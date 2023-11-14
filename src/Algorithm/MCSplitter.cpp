#include "C4Hex/Algorithm/MCSplitter.hpp"

#include <MC3D/Mesh/MCMeshManipulator.hpp>

#include <MC3D/Algorithm/MotorcycleTracer.hpp>
#include "C4Hex/Algorithm/PathRouter.hpp"
#include "C4Hex/Algorithm/SurfaceRouter.hpp"


#include <fstream>

namespace c4hex
{

MCSplitter::MCSplitter(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

bool MCSplitter::bisectPatchAcrossDir(const FH& p, UVWDir dir, bool allowLoop, vector<FH>& psSub)
{
    psSub.clear();

    // Quantization must be present
    assert(mcMeshProps().isAllocated<ARC_INT_LENGTH>());

    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    HFH hpAny = mcMesh.halfface_handle(p, 0);
    if (!mcMesh.incident_cell(hpAny).is_valid())
        hpAny = mcMesh.halfface_handle(p, 1);

    auto side2has = halfpatchHalfarcsByDir(hpAny);

    UVWDir dir1 = findMatching(side2has,
                               [&](const pair<const UVWDir, vector<HEH>>& kv)
                               { return (dim(dir) != 1 || (kv.first & dir) == UVWDir::NONE) && kv.second.size() > 1; })
                      .first;

    if (dir1 == UVWDir::NONE)
        return false;

    for (UVWDir dir2 : {dir1, -dir1})
        if (containsMatching(side2has.at(dir2),
                             [&](const HEH& ha)
                             {
                                 return mcMeshProps().isAllocated<ARC_INT_LENGTH>()
                                        && isZeroArc(mcMesh.edge_handle(ha));
                             }))
        {
            DLOG(INFO) << "Can not refine patch " << p << " safely yet, has 0-arc";
            return false;
        }

    DLOG(INFO) << "Bisecting patch " << p << " side along dir " << toVec(dir1);

    UVWDir dir1opp = -dir1;

    auto& hasDir1 = side2has.at(dir1);
    auto& hasDir1opp = side2has.at(dir1opp);

    HEH haFirst1 = hasDir1.front();
    HEH haLast1opp = mcMesh.opposite_halfedge_handle(hasDir1opp.back());

    cutToMatch(haFirst1, haLast1opp);

    VH n1 = mcMesh.to_vertex_handle(haFirst1);
    VH n1opp = mcMesh.to_vertex_handle(haLast1opp);
    list<HEH> pathHesStart = mcMeshProps().haHalfedges(mcMesh.opposite_halfedge_handle(haLast1opp));
    list<HEH> pathHesEnd = mcMeshProps().haHalfedges(haFirst1);

    // Find the side connecting haLast1opp -> ... -> haFirst1 and store the full path
    list<HEH> pathHes = pathHesStart;
    int orthIntLength = 0;
    double orthDblLength = 0.0;
    for (HEH ha : findMatching(side2has,
                               [&](const pair<const UVWDir, vector<HEH>>& kv)
                               {
                                   return mcMesh.to_vertex_handle(kv.second.back())
                                              == mcMesh.from_vertex_handle(haFirst1)
                                          && kv.first != dir1 && kv.first != -dir1;
                               })
                      .second)
    {
        EH a = mcMesh.edge_handle(ha);
        auto hesHa = mcMeshProps().haHalfedges(ha);
        pathHes.insert(pathHes.end(), hesHa.begin(), hesHa.end());
        orthIntLength += mcMeshProps().get<ARC_INT_LENGTH>(a);
        orthDblLength += mcMeshProps().isAllocated<ARC_DBL_LENGTH>() ? mcMeshProps().get<ARC_DBL_LENGTH>(a) : 0.0;
    }

    if (!allowLoop && mcMeshProps().isAllocated<ARC_INT_LENGTH>() && orthIntLength == 0 && n1 == n1opp)
    {
        DLOG(WARNING) << "Can not refine patch " << p << " safely, as it is (will be) a cigar patch!";
        return false;
    }

    // Find patch bisector path in embedding
    pathHes.insert(pathHes.end(), pathHesEnd.begin(), pathHesEnd.end());
    {
        set<HFH> hfsTransferred;
        auto ret = PathRouter(meshProps()).reroutePathThroughPatch(p, pathHes, hfsTransferred);
        if (ret != PathRouter::SUCCESS)
            throw std::logic_error("Can not reroute path through patch, programming error");
        assert(!hfsTransferred.empty());
    }

    // Add new face bisector arc and give it all needed properties
    EH aSplitting = mcMesh.add_edge(n1opp, n1, true);

    if (mcMeshProps().isAllocated<IS_SINGULAR>())
        mcMeshProps().set<IS_SINGULAR>(aSplitting, false);
    mcMeshProps().set<ARC_INT_LENGTH>(aSplitting, orthIntLength);
    if (mcMeshProps().isAllocated<ARC_DBL_LENGTH>())
    {
        if (orthDblLength > 0)
            mcMeshProps().set<ARC_DBL_LENGTH>(aSplitting, orthDblLength);
        else
        {
            double length = 0.0;
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(aSplitting))
                length += edgeLengthUVW<CHART>(tetMesh.edge_handle(he));
            mcMeshProps().set<ARC_DBL_LENGTH>(aSplitting, length);
        }
    }
    if (mcMeshProps().isAllocated<ARC_MESH_HALFEDGES>())
        mcMeshProps().set<ARC_MESH_HALFEDGES>(aSplitting, pathHes);

    if (meshProps().isAllocated<MC_ARC>())
        for (HEH he : pathHes)
            meshProps().set<MC_ARC>(tetMesh.edge_handle(he), aSplitting);
    if (meshProps().isAllocated<IS_ARC>())
        for (HEH he : pathHes)
            meshProps().set<IS_ARC>(tetMesh.edge_handle(he), true);

    set<CH> affectedBs;

    // Use bisector arc to bisect patch
    psSub = splitPatch(p, aSplitting, affectedBs);

    if (mcMeshProps().isAllocated<ARC_INT_LENGTH>() && orthIntLength == 0)
        if (n1 == n1opp)
            cutToMatchLength(mcMesh.halfedge_handle(aSplitting, 0), 0);

    return true;
}

bool MCSplitter::bisectBlockOrPatch(const CH& b, vector<CH>& subBlocks)
{
    assert(!mcMeshProps().mesh().is_deleted(b));
    assert(mcMeshProps().ref<BLOCK_FACE_PATCHES>(b).size() == DIM_1_DIRS.size());

    bool changing = true;
    bool diff = false;

    while (changing)
    {
        assert(!mcMeshProps().mesh().is_deleted(b));
        assert(mcMeshProps().ref<BLOCK_FACE_PATCHES>(b).size() == DIM_1_DIRS.size());

        map<int, vector<FH>> nPatchesPerSide2ps;
        for (UVWDir dir : DIM_1_DIRS)
        {
            auto& facePatches = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b).at(dir);
            int nPatchesPerSide = -((int)facePatches.size());
            for (FH p : facePatches)
                nPatchesPerSide2ps[nPatchesPerSide].push_back(p);
        }
        changing = containsMatching(nPatchesPerSide2ps,
                                    [&](const pair<const int, vector<FH>>& kv)
                                    {
                                        return containsMatching(
                                            kv.second,
                                            [&](const FH& p)
                                            {
                                                vector<FH> trash;
                                                return bisectPatchAcrossDir(p, UVWDir::ANY, false, trash); });
                                    });
        if (changing)
            diff = true;
    }

    bool boundaryFinished = !containsMatching(mcMeshProps().mesh().cell_halffaces(b),
                                              [&](const HFH& hp)
                                              {
                                                  return containsMatching(halfpatchHalfarcsByDir(hp),
                                                                          [&](const pair<const UVWDir, vector<HEH>>& kv)
                                                                          { return kv.second.size() > 1; });
                                              });

    if (boundaryFinished)
        diff |= bisectBlockInterior(b, subBlocks);

    return diff;
}

FH MCSplitter::cutBlock(const CH& b, const vector<HEH>& hasRing, vector<CH>& bsCut)
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    set<EH> ringAs;
    for (HEH ha : hasRing)
        ringAs.insert(mcMesh.edge_handle(ha));

    HFH seedHp = findMatching(mcMesh.halfedge_halffaces(hasRing.front()),
                              [&](const HFH& hp) { return mcMesh.incident_cell(hp) == b; });

    // floodfill halfpatches on one side of the block:
    list<HFH> hpQ({seedHp});
    set<HFH> sideHps({seedHp});
    while (!hpQ.empty())
    {
        HFH hp = hpQ.front();
        hpQ.pop_front();
        for (HEH ha3 : mcMesh.halfface_halfedges(hp))
        {
            if (ringAs.count(mcMesh.edge_handle(ha3)) != 0)
                continue;
            HFH hpNext = mcMesh.adjacent_halfface_in_cell(hp, ha3);
            if (!hpNext.is_valid() || sideHps.find(hpNext) != sideHps.end())
                continue;
            hpQ.emplace_back(hpNext);
            sideHps.insert(hpNext);
        }
    }

    // Gather halffaces in a set
    set<HFH> pNewHfs;
    if (mcMeshProps().isAllocated<PATCH_MESH_HALFFACES>())
        for (HFH hp : sideHps)
        {
            auto hfsHp = mcMeshProps().hpHalffaces(hp);
            pNewHfs.insert(hfsHp.begin(), hfsHp.end());
        }

    // Reroute halffaces
    set<CH> transferredTets;
    if (SurfaceRouter(meshProps()).rerouteSurfaceThroughBlock(b, pNewHfs, transferredTets) != SurfaceRouter::SUCCESS)
        throw std::logic_error("Block not properly connected");

    FH pNew = mcMesh.add_face(hasRing);
    // Default transition is identity, no need to adjust
    if (mcMeshProps().isAllocated<PATCH_MIN_DIST>())
        mcMeshProps().set<PATCH_MIN_DIST>(pNew, DBL_MAX);
    if (mcMeshProps().isAllocated<PATCH_MESH_HALFFACES>())
        mcMeshProps().set<PATCH_MESH_HALFFACES>(pNew, pNewHfs);

    if (meshProps().isAllocated<MC_PATCH>())
        for (HFH hf : pNewHfs)
            meshProps().set<MC_PATCH>(tetMesh.face_handle(hf), pNew);
    if (meshProps().isAllocated<IS_WALL>())
        for (HFH hf : pNewHfs)
            meshProps().set<IS_WALL>(tetMesh.face_handle(hf), true);
    if (meshProps().isAllocated<WALL_DIST>())
        for (HFH hf : pNewHfs)
            meshProps().set<WALL_DIST>(tetMesh.face_handle(hf), DBL_MAX);

    bsCut = splitBlock(b, pNew);

    return pNew;
}

bool MCSplitter::bisectBlockInterior(const CH& b, vector<CH>& subBlocks)
{
    auto& mcMesh = mcMeshProps().mesh();

    for (auto& dir2faceArcs : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
        if (!dir2faceArcs.second.empty())
        {
            HEH haSeed = mcMesh.halfedge_handle(*dir2faceArcs.second.begin(), 0);

            auto hasRing = findArcLoopOnBoundary(b, haSeed);
            if (hasRing.empty())
            {
                DLOG(INFO) << "Could not close surface loop on block " << b << " delaying refinement";
                continue;
            }

            set<EH> asRing;
            for (HEH ha : hasRing)
                asRing.insert(mcMesh.edge_handle(ha));

            // See if loop really partitions block boundary into 2 parts
            int i = 0;
            set<HFH> hpsVisited;
            for (HFH hpStart : mcMesh.cell_halffaces(b))
            {
                if (hpsVisited.count(hpStart) != 0)
                    continue;

                hpsVisited.insert(hpStart);
                list<HFH> hpQ({hpStart});
                while (!hpQ.empty())
                {
                    HFH hp = hpQ.front();
                    hpQ.pop_front();
                    for (HEH ha3 : mcMesh.halfface_halfedges(hp))
                    {
                        if (asRing.count(mcMesh.edge_handle(ha3)) != 0)
                            continue;
                        HFH hpNext = mcMesh.adjacent_halfface_in_cell(hp, ha3);
                        if (!hpNext.is_valid() || hpsVisited.find(hpNext) != hpsVisited.end())
                            continue;
                        hpQ.emplace_back(hpNext);
                        hpsVisited.insert(hpNext);
                    }
                }
                i++;
            }
            if (i != 2)
            {
                DLOG(INFO) << "Block surface loop is nonmanifold and splits surface into != 2 patches";
                continue;
            }

            vector<CH> bsCut;
            cutBlock(b, hasRing, bsCut);

            subBlocks = {bsCut[0], bsCut[1]};
            return true;
        }

    return false;
}

vector<HEH> MCSplitter::findArcLoopOnBoundary(const CH& b, const HEH& haSeed)
{
    auto& mcMesh = mcMeshProps().mesh();

    // Collect ring
    vector<HEH> hasRing;
    UVWDir allDirs = UVWDir::NONE;
    HEH haCurr = haSeed;
    set<EH> asVisited;
    do
    {
        hasRing.emplace_back(haCurr);
        asVisited.insert(mcMesh.edge_handle(haCurr));
        VH nCurr = mcMesh.to_vertex_handle(haCurr);
        bool onCorner = containsMatching(mcMeshProps().ref<BLOCK_CORNER_NODES>(b),
                                         [&](const pair<const UVWDir, VH>& corner) { return corner.second == nCurr; });
        if (onCorner)
            break;
        bool onEdge
            = containsMatching(mcMeshProps().ref<BLOCK_EDGE_NODES>(b),
                               [&](const pair<const UVWDir, set<VH>>& kv) { return kv.second.count(nCurr) != 0; });

        haCurr = findMatching(
            mcMesh.outgoing_halfedges(nCurr),
            [&](const HEH& ha)
            {
                if (ha == haSeed || asVisited.count(mcMesh.edge_handle(ha)) == 0)
                {
                    if (onEdge)
                    {
                        UVWDir sideCurr = findMatching(mcMeshProps().ref<BLOCK_FACE_ARCS>(b),
                                                       [&](const pair<const UVWDir, set<EH>>& kv)
                                                       { return kv.second.count(mcMesh.edge_handle(haCurr)) != 0; })
                                              .first;
                        UVWDir sideNew = findMatching(mcMeshProps().ref<BLOCK_FACE_ARCS>(b),
                                                      [&](const pair<const UVWDir, set<EH>>& kv)
                                                      { return kv.second.count(mcMesh.edge_handle(ha)) != 0; })
                                             .first;
                        if (sideNew != sideCurr && sideNew != UVWDir::NONE)
                            return !containsMatching(mcMeshProps().ref<BLOCK_EDGE_ARCS>(b),
                                                     [&](const pair<const UVWDir, set<EH>>& kv)
                                                     { return kv.second.count(mcMesh.edge_handle(ha)) != 0; });
                    }
                    else if (halfarcDirInBlock(ha, b) == halfarcDirInBlock(haCurr, b))
                        return true;
                }
                return false;
            });
    } while (haCurr.is_valid() && haCurr != haSeed);

    if (dim(allDirs) != 3 && mcMesh.to_vertex_handle(hasRing.back()) == mcMesh.from_vertex_handle(hasRing.front()))
        return hasRing;

    return {};
}

void MCSplitter::cutToMatch(HEH& ha1, HEH& ha2)
{
    auto& mcMesh = mcMeshProps().mesh();

    int length1 = mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha1));
    int length2 = mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha2));

    if (length1 < length2)
        ha2 = cutToMatchLength(ha2, length1)[0];
    else if (length2 < length1)
        ha1 = cutToMatchLength(ha1, length2)[0];
}

VH MCSplitter::findVertexAtLengthFraction(const EH& a, double frac)
{
    auto& tetMesh = meshProps().mesh();

    auto& aHes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);
    if (aHes.size() == 1)
        splitHalfEdge(aHes.front(), *tetMesh.hec_iter(aHes.front()), Q(0.5));
    assert(aHes.size() >= 2);
    double fullLen = 0;
    for (HEH he : aHes)
        fullLen += edgeLengthUVW<CHART>(tetMesh.edge_handle(he));
    double len = 0;
    VH vSplit;
    for (HEH he : aHes)
    {
        if (len / fullLen > frac || he == aHes.back())
        {
            vSplit = tetMesh.from_vertex_handle(he);
            break;
        }
        len += edgeLengthUVW<CHART>(tetMesh.edge_handle(he));
    }
    return vSplit;
}

vector<HEH> MCSplitter::cutToMatchLength(const HEH& ha, int lengthToMatch)
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    EH a = mcMesh.edge_handle(ha);
    int lengthHa = mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));

    bool isHalfarc0 = (ha.idx() % 2) == 0;
    double frac = lengthHa == 0 ? 0.5 : (double)lengthToMatch / lengthHa;
    if (!isHalfarc0)
        frac = 1. - frac;

    frac = std::min(0.9, std::max(0.1, frac));

    VH vSplit = findVertexAtLengthFraction(mcMesh.edge_handle(ha), frac);

    VH nSplit = mcMesh.add_vertex(tetMesh.vertex(vSplit));
    if (meshProps().isAllocated<MC_NODE>())
        meshProps().set<MC_NODE>(vSplit, nSplit);
    if (mcMeshProps().isAllocated<NODE_MESH_VERTEX>())
        mcMeshProps().set<NODE_MESH_VERTEX>(nSplit, vSplit);

    set<FH> affectedPs;
    set<CH> affectedBs;
    auto subArcs = splitArc(a, nSplit, affectedPs, affectedBs);
    if (mcMeshProps().isAllocated<ARC_INT_LENGTH>())
    {
        DLOG(INFO) << "Split arc " << a << " into child arcs " << subArcs[0] << " and " << subArcs[1];
        mcMeshProps().set<ARC_INT_LENGTH>(subArcs[0], isHalfarc0 ? lengthToMatch : lengthHa - lengthToMatch);
        mcMeshProps().set<ARC_INT_LENGTH>(subArcs[1], isHalfarc0 ? lengthHa - lengthToMatch : lengthToMatch);
    }

    if (isHalfarc0)
        return {mcMesh.halfedge_handle(subArcs[0], 0), mcMesh.halfedge_handle(subArcs[1], 0)};
    return {mcMesh.halfedge_handle(subArcs[1], 1), mcMesh.halfedge_handle(subArcs[0], 1)};
}

} // namespace c4hex
