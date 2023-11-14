#include "C4Hex/Algorithm/IGMInitializer.hpp"

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#define UNIFORM_SPACING

namespace c4hex
{

IGMInitializer::IGMInitializer(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

void IGMInitializer::allocateIGM()
{
    mcMeshProps().allocate<PATCH_IGM_TRANSITION>();
    meshProps().allocate<CHART_IGM>();
    meshProps().allocate<TRANSITION_IGM>();

    for (CH tet : meshProps().mesh().cells())
    {
        auto& chart = meshProps().ref<CHART_IGM>(tet);
        for (VH v : meshProps().mesh().tet_vertices(tet))
            chart[v] = Vec3Q(0, 0, 0);
    }
}

IGMInitializer::RetCode IGMInitializer::initializeFromQuantization()
{
    allocateIGM();

    initializeOnNodes();

    determineIGMTransitions();

    initializeOnArcs();

    splitAllNecessaryForInjectiveBoundary();

    initializeOnPatches();

    // We dont split the sufficient amount to avoid as much unnecessary refinement as possible.
    // Splitting the rest of elements is left to the untangling routine.
    splitSomeNecessaryForInjectiveInterior();

    initializeInBlocks();

    for (CH tet : meshProps().mesh().cells())
    {
        if (doubleVolumeIGM(tet) <= 1e-6 && rationalVolumeIGM(tet) <= 0)
        {
            LOG(WARNING) << "Initial IGM obtained via naive 3D Tutte has invalid elements, IGM untangling is needed";
            return INVALID_ELEMENTS;
        }
    }

    return SUCCESS;
}

bool IGMInitializer::splitAllNecessaryForInjectiveBoundary()
{
    const MCMesh& mcMesh = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();
    set<EH> cutEdges;
    set<FH> cutFaces;

    for (EH a : mcMesh.edges())
        for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
            for (VH v : mesh.halfedge_vertices(he))
                for (HEH he2 : mesh.outgoing_halfedges(v))
                {
                    if (cutEdges.count(mesh.edge_handle(he2)) != 0)
                        continue;
                    if (!meshProps().isInArc(he2))
                    {
                        VH v2 = mesh.to_vertex_handle(he2);
                        for (EH e3 : mesh.vertex_edges(v2))
                        {
                            EH a3 = meshProps().get<MC_ARC>(e3);
                            if (a3.is_valid())
                            {
                                if (a3 == a)
                                {
                                    cutEdges.insert(mesh.edge_handle(he2));
                                    break;
                                }
                                for (CH tet : mesh.halfedge_cells(he2))
                                {
                                    CH b = meshProps().get<MC_BLOCK>(tet);
                                    UVWDir dirHa = halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b);
                                    if (dirHa == UVWDir::NONE)
                                        continue;
                                    UVWDir dirHa3 = halfarcDirInBlock(mcMesh.halfedge_handle(a3, 0), b);
                                    if (dirHa3 == UVWDir::NONE)
                                        continue;
                                    dirHa = dirHa | -dirHa;
                                    dirHa3 = dirHa3 | -dirHa3;
                                    if (dirHa == dirHa3)
                                    {
                                        Vec3Q igm1 = meshProps().ref<CHART_IGM>(tet).at(v);
                                        Vec3Q igm2 = meshProps().ref<CHART_IGM>(tet).at(v2);
                                        if (dim(toDir(igm1 - igm2)) <= 1
                                            && toCoord(toDir(igm1 - igm2)) == toCoord(dirHa))
                                        {
                                            cutEdges.insert(mesh.edge_handle(he2));
                                            break;
                                        }
                                    }
                                }
                                if (cutEdges.count(mesh.edge_handle(he2)) != 0)
                                    break;
                            }
                        }
                    }
                }
    for (FH p : mcMesh.faces())
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            for (VH v : mesh.halfface_vertices(hf))
            {
                for (HEH he2 : mesh.outgoing_halfedges(v))
                {
                    if (cutEdges.count(mesh.edge_handle(he2)) != 0 || meshProps().isInPatch(he2))
                        continue;

                    VH v2 = mesh.to_vertex_handle(he2);
                    for (CH tet : mesh.halfedge_cells(he2))
                    {
                        CH b = meshProps().get<MC_BLOCK>(tet);
                        HFH hp = mcMesh.halfface_handle(p, 0);
                        if (mcMesh.incident_cell(hp) != b)
                            hp = mcMesh.opposite_halfface_handle(hp);
                        if (mcMesh.incident_cell(hp) != b)
                            continue;
                        bool bothHpsOnB = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp)) == b;
                        UVWDir dirP = halfpatchNormalDir(hp);
                        for (FH f2 : mesh.vertex_faces(v2))
                        {
                            FH p2 = meshProps().get<MC_PATCH>(f2);
                            if (p2.is_valid())
                            {
                                HFH hp2 = mcMesh.halfface_handle(p2, 0);
                                if (mcMesh.incident_cell(hp2) != b)
                                    hp2 = mcMesh.opposite_halfface_handle(hp2);
                                if (mcMesh.incident_cell(hp2) != b)
                                    continue;
                                UVWDir dirP2 = halfpatchNormalDir(hp2);
                                if (bothHpsOnB)
                                    dirP2 = dirP2 | -dirP2;
                                if ((dirP2 & dirP) != UVWDir::NONE)
                                    cutEdges.insert(mesh.edge_handle(he2));
                                break;
                            }
                        }
                        if (cutEdges.count(mesh.edge_handle(he2)) != 0)
                            break;
                    }
                }
                for (FH f2 : mesh.vertex_faces(v))
                {
                    if (cutFaces.count(f2) != 0 || meshProps().isInPatch(f2))
                        continue;
                    for (CH tet : mesh.face_cells(f2))
                    {
                        if (!tet.is_valid())
                            continue;

                        CH b = meshProps().get<MC_BLOCK>(tet);
                        HFH hp = mcMesh.halfface_handle(p, 0);
                        if (mcMesh.incident_cell(hp) != b)
                            hp = mcMesh.opposite_halfface_handle(hp);
                        if (mcMesh.incident_cell(hp) != b)
                            continue;

                        UVWDir dirP = halfpatchNormalDir(hp);

                        int nAligned = 0;
                        for (EH e : mesh.face_edges(f2))
                        {
                            bool aligned = false;
                            for (FH f3 : mesh.edge_faces(e))
                            {
                                FH p3 = meshProps().get<MC_PATCH>(f3);
                                if (p3.is_valid())
                                {
                                    HFH hp3 = mcMesh.halfface_handle(p3, 0);
                                    if (mcMesh.incident_cell(hp3) != b)
                                        hp3 = mcMesh.opposite_halfface_handle(hp3);
                                    if (mcMesh.incident_cell(hp3) != b)
                                        continue;
                                    UVWDir dirP3 = halfpatchNormalDir(hp3);
                                    if (dirP3 == dirP)
                                    {
                                        aligned = true;
                                        break;
                                    }
                                }
                            }
                            if (aligned)
                                nAligned++;
                        }
                        if (nAligned == 3)
                            cutFaces.insert(f2);
                    }
                }
            }

    for (EH e : cutEdges)
        splitHalfEdge(mesh.halfedge_handle(e, 0), *mesh.ec_iter(e), Q(0.5));
    for (FH f : cutFaces)
        if (!mesh.is_deleted(f))
            splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});
    LOG(INFO) << "Split " << cutEdges.size() << " edges and " << cutFaces.size()
              << " faces, to prevent degenerate boundary IGM";
    return !cutEdges.empty() || !cutFaces.empty();
}

bool IGMInitializer::splitSomeNecessaryForInjectiveInterior()
{
    // const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    set<EH> cutEdges;
    set<FH> cutFaces;
    // Care only about tets which are inverted and all 4 vertices constrained on boundary
    for (CH tet : mesh.cells())
        if (rationalVolumeIGM(tet) <= 0
            && !containsMatching(mesh.tet_vertices(tet), [this](const VH& v) { return !meshProps().isInPatch(v); }))
            for (EH e : mesh.cell_edges(tet))
                if (!meshProps().isInPatch(e))
                    cutEdges.insert(e);
    for (EH e : cutEdges)
        splitHalfEdge(mesh.halfedge_handle(e, 0), *mesh.ec_iter(e), Q(0.5));
    for (FH f : cutFaces)
        splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});
    LOG(INFO) << "Split " << cutEdges.size() << " edges and " << cutFaces.size()
              << " faces, to prevent overconstrained inverted interior IGM";
    return !cutEdges.empty() || !cutFaces.empty();
}

bool IGMInitializer::splitAllSufficientForInjectiveInterior()
{
    // const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    set<EH> cutEdges;
    set<FH> cutFaces;
    for (EH e : mesh.edges())
    {
        if (!meshProps().isInPatch(e)
            && !containsMatching(mesh.edge_vertices(e), [this](const VH& v) { return !meshProps().isInPatch(v); }))
            cutEdges.insert(e);
    }
    for (FH f : mesh.faces())
    {
        if (meshProps().isInPatch(f))
            continue;
        int nEonPatch = 0;
        for (EH e : mesh.face_edges(f))
            if (meshProps().isInPatch(e))
                nEonPatch++;
        if (nEonPatch == 3)
            cutFaces.insert(f);
    }
    for (EH e : cutEdges)
        splitHalfEdge(mesh.halfedge_handle(e, 0), *mesh.ec_iter(e), Q(0.5));
    for (FH f : cutFaces)
        splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});
    LOG(INFO) << "Split " << cutEdges.size() << " edges and " << cutFaces.size()
              << " faces, to guarantee existence of bijective map";
    return !cutEdges.empty() || !cutFaces.empty();
}

void IGMInitializer::minimizeCutSurface()
{
    const MCMesh& mc = mcMeshProps().mesh();
    set<CH> bsVisited({*mc.c_iter()});
    list<CH> bQ({*mc.c_iter()});

    while (!bQ.empty())
    {
        CH b = bQ.front();
        bQ.pop_front();

        for (HFH hpOut : mc.cell_halffaces(b))
        {
            HFH hpIn = mc.opposite_halfface_handle(hpOut);
            CH bNext = mc.incident_cell(hpIn);
            if (bNext.is_valid() && bsVisited.count(bNext) == 0)
            {
                bsVisited.insert(bNext);
                bQ.push_back(bNext);
                Transition trans = mcMeshProps().hpTransition<PATCH_IGM_TRANSITION>(hpIn);
                Transition transNoTranslation = trans;
                transNoTranslation.translation = Vec3Q(0, 0, 0);
                applyTransitionToBlock(transNoTranslation, bNext);
                applyTransitionIGMToBlock(trans, bNext, false);
            }
        }
    }
}

void IGMInitializer::initializeOnNodes()
{
    const MCMesh& mc = mcMeshProps().mesh();

    for (CH b : mc.cells())
    {
        Vec3Q startUVW(0, 0, 0);
        VH n0 = *mc.cv_iter(b);
        VH v0 = mcMeshProps().get<NODE_MESH_VERTEX>(n0);
        CH tet0 = anyIncidentTetOfBlock(v0, b);

        set<pair<VH, CH>> tetsVisited;
        forVertexNeighbourTetsInBlock(v0,
                                      tet0,
                                      [this, &tetsVisited, n0](const CH& tet)
                                      {
                                          tetsVisited.insert({n0, tet});
                                          return false;
                                      });
        list<pair<pair<VH, CH>, Vec3Q>> nQ({{{n0, tet0}, startUVW}});

        set<HEH> bHas;
        for (HEH ha : mc.cell_halfedges(b))
            bHas.insert(ha);

        while (!nQ.empty())
        {
            auto nodeTetUVW = nQ.front();
            nQ.pop_front();
            VH n = nodeTetUVW.first.first;
            Vec3Q UVW = nodeTetUVW.second;
            VH v = mcMeshProps().get<NODE_MESH_VERTEX>(n);
            CH tet = nodeTetUVW.first.second;

            // Flood UVW to all tets around v
            forVertexNeighbourTetsInBlock(v,
                                          tet,
                                          [this, &v, &UVW](const CH& tet2)
                                          {
                                              meshProps().ref<CHART_IGM>(tet2)[v] = UVW;
                                              return false; // dont break
                                          });

            // Find next nodes in MC
            for (HEH ha : mc.outgoing_halfedges(n))
            {
                if (bHas.count(ha) == 0)
                    continue;

                VH nNext = mc.to_vertex_handle(ha);
                VH vNext = mcMeshProps().get<NODE_MESH_VERTEX>(nNext);

                // Find volumetric path from current vertex to target vertex along ha
                set<CH> tetsPath;
                CH tetNext = traceToEndWithinBlock(tet, ha, tetsPath);
                if (tetNext.is_valid() && tetsVisited.count({nNext, tetNext}) == 0)
                {
                    forVertexNeighbourTetsInBlock(vNext,
                                                  tetNext,
                                                  [this, &tetsVisited, nNext](const CH& tet2)
                                                  {
                                                      tetsVisited.insert({nNext, tet2});
                                                      return false;
                                                  });
                    Vec3Q UVWnext = UVW
                                    + Vec3Q(toVec(halfarcDirInBlock(ha, b))
                                            * mcMeshProps().get<ARC_INT_LENGTH>(mc.edge_handle(ha)));
                    nQ.push_back({{nNext, tetNext}, UVWnext});
                }
            }
        }
    }
}

void IGMInitializer::determineIGMTransitions()
{
    const MCMesh& mc = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    for (FH p : mc.faces())
    {
        if (mc.is_boundary(p))
            continue;

        auto& hfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);

        auto itPairPNs = mc.face_vertices(p);
        auto ns = vector<VH>(itPairPNs.first, itPairPNs.second);
        auto UVWs0 = map<VH, Vec3Q>();
        auto UVWs1 = map<VH, Vec3Q>();
        for (VH n : ns)
        {
            VH v = mcMeshProps().get<NODE_MESH_VERTEX>(n);
            CH tet0, tet1;
            for (bool front : {true, false})
                for (HFH hf : tetMesh.vertex_halffaces(v))
                    if (hfs.count(front ? hf : tetMesh.opposite_halfface_handle(hf)) != 0)
                    {
                        (front ? tet0 : tet1) = tetMesh.incident_cell(hf);
                        break;
                    }
            assert(tet0.is_valid() && tet1.is_valid());
            UVWs0[n] = meshProps().ref<CHART_IGM>(tet0).at(v);
            UVWs1[n] = meshProps().ref<CHART_IGM>(tet1).at(v);
        }
        Transition transIGM = mcMeshProps().get<PATCH_TRANSITION>(p);
        Vec3Q deltaIGM = transIGM.rotate(UVWs0.begin()->second) - UVWs1.begin()->second;
        transIGM.translation = -deltaIGM;

        mcMeshProps().set<PATCH_IGM_TRANSITION>(p, transIGM);
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            meshProps().setTransition<TRANSITION_IGM>(hf, transIGM);
    }
}

void IGMInitializer::initializeOnArcs()
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    for (EH a : mc.edges())
    {
        HEH ha = mc.halfedge_handle(a, 0);
        const auto& hes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);

        VH nStart = mc.from_vertex_handle(ha);
        VH nEnd = mc.to_vertex_handle(ha);
        VH vStart = mcMeshProps().get<NODE_MESH_VERTEX>(nStart);
        VH vEnd = mcMeshProps().get<NODE_MESH_VERTEX>(nEnd);
        assert(mesh.from_vertex_handle(*hes.begin()) == vStart && mesh.to_vertex_handle(*hes.rbegin()) == vEnd);

        set<CH> pathTets;
        CH tetStart = *mesh.hec_iter(hes.front());
        CH tetEnd = traceToEndWithinBlock(tetStart, ha, pathTets);
        Vec3Q startIGM = meshProps().ref<CHART_IGM>(tetStart).at(vStart);
        Vec3Q endIGM = meshProps().ref<CHART_IGM>(tetEnd).at(vEnd);

        double totalDist = 0.0;
        map<VH, Q> v2dist;
        for (HEH he : hes)
        {
            VH v2 = mesh.to_vertex_handle(he);
#ifndef UNIFORM_SPACING
            VH v1 = mesh.from_vertex_handle(he);
            double heLength = std::max(1e-4, (mesh.vertex(v2) - mesh.vertex(v1)).length());
            totalDist += heLength;
#else
            totalDist += 1;
#endif
            v2dist[v2] = totalDist;
        }

        for (const auto& kv : v2dist)
        {
            VH v = kv.first;
            Q dist = kv.second;
            if (v == vEnd)
                continue;
            for (CH tet0 : mesh.vertex_cells(v))
                if (pathTets.count(tet0) != 0)
                {
                    Vec3Q igm0 = startIGM + dist / totalDist * (endIGM - startIGM);
                    auto tet2trans = determineTransitionsAroundVertex<TRANSITION_IGM>(v, tet0);
                    for (auto& kv2 : tet2trans)
                        meshProps().ref<CHART_IGM>(kv2.first)[v] = kv2.second.apply(igm0);
                    break;
                }
        }
    }
}

void IGMInitializer::initializeOnPatches()
{
#ifdef UNIFORM_SPACING
    initializeOnPatches2DTutte(false, false);
#else
    initializeOnPatches2DTutte(true, true);
#endif
}

void IGMInitializer::initializeOnPatches2DTutte(bool meanValWeights, bool xyzTarget)
{

    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    for (FH p : mc.faces())
    {
        CH b = mc.incident_cell(mc.halfface_handle(p, 0));
        if (!b.is_valid())
            b = mc.incident_cell(mc.halfface_handle(p, 1));
        Transition transIGM = mcMeshProps().get<PATCH_IGM_TRANSITION>(p);

        auto& hfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);

        int isoCoord = toCoord(halfpatchNormalDir(mc.halfface_handle(p, 0)));
        Q isoValue(0);
        bool front = true;
        {
            VH vCorner = mcMeshProps().get<NODE_MESH_VERTEX>(*mc.fv_iter(p));
            CH tet0, tet1;
            for (bool checkFront : {true, false})
                for (HFH hf : mesh.vertex_halffaces(vCorner))
                    if (hfs.count(checkFront ? hf : mesh.opposite_halfface_handle(hf)) != 0)
                    {
                        (checkFront ? tet0 : tet1) = mesh.incident_cell(hf);
                        break;
                    }
            assert(tet0.is_valid() || tet1.is_valid());
            if (tet0.is_valid())
            {
                auto chart = meshProps().ref<CHART_IGM>(tet0);
                isoValue = chart.at(vCorner)[isoCoord];
            }
            else
            {
                front = false;
                auto chart = meshProps().ref<CHART_IGM>(tet1);
                isoValue = chart.at(vCorner)[isoCoord];
            }
        }

        map<VH, int> vtx2index;
        set<EH> innerEdges;
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
        {
            for (VH v : mesh.halfface_vertices(hf))
            {
                if (meshProps().isInArc(v))
                    continue;
                int nextIdx = vtx2index.size();
                if (vtx2index.find(v) == vtx2index.end())
                    vtx2index[v] = nextIdx;
            }
            for (EH e : mesh.halfface_edges(hf))
                if (!meshProps().isInArc(e))
                    innerEdges.insert(e);
        }
        if (vtx2index.size() == 0)
            continue;

        map<EH, double> e2weight;
        if (meanValWeights)
            e2weight = calc2DmeanValueWeights(b, innerEdges, xyzTarget);

        solve2DTutte<double>(vtx2index, e2weight, isoCoord, isoValue, front, hfs);

        bool hasFlipsInDoublePrecision = false;
        {
            HFH hp = mcMeshProps().mesh().halfface_handle(p, 0);
            if (!mcMeshProps().mesh().incident_cell(hp).is_valid())
                hp = mcMeshProps().mesh().opposite_halfface_handle(hp);
            VH nMin = mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::NEG_U_NEG_V_NEG_W);
            VH nMax = mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::POS_U_POS_V_POS_W);
            Vec3Q minIGM = nodeIGMinBlock(nMin, b);
            Vec3Q maxIGM = nodeIGMinBlock(nMax, b);
            Vec3Q midPt = (minIGM + maxIGM) * 0.5;
            for (HFH hf : mcMeshProps().hpHalffaces(hp))
            {
                CH tet = meshProps().mesh().incident_cell(hf);
                vector<Vec3Q> uvws;
                for (VH v : meshProps().mesh().halfface_vertices(hf))
                    uvws.push_back(meshProps().ref<CHART_IGM>(tet).at(v));
                if (dot(cross(uvws[1] - uvws[0], uvws[2] - uvws[0]), midPt - uvws[0]) <= 0)
                {
                    LOG(WARNING) << "Boundary hf " + std::to_string(hf.idx()) + " of block " + std::to_string(b.idx())
                                        + " inverted! Have to recompute tutte in exact coords";
                    hasFlipsInDoublePrecision = true;
                }
            }
        }

        if (hasFlipsInDoublePrecision)
            solve2DTutte<Q>(vtx2index, e2weight, isoCoord, isoValue, front, hfs);
    }
}

map<EH, double> IGMInitializer::calc2DmeanValueWeights(const CH& b, const set<EH>& edges, bool xyzTarget) const
{
    const TetMesh& mesh = meshProps().mesh();

    map<EH, double> e2meanValWeight;

    for (EH e : edges)
    {
        HEH heOut = mesh.halfedge_handle(e, 0);

        if (e2meanValWeight.find(e) != e2meanValWeight.end())
            continue;
        vector<HFH> adjPatchHfs;
        for (HFH hf : mesh.halfedge_halffaces(heOut))
            if (meshProps().isInPatch(mesh.face_handle(hf)))
                adjPatchHfs.emplace_back(hf);
        if (adjPatchHfs.size() != 2)
            continue;

        vector<CH> adjBlockCells;
        for (HFH hf : adjPatchHfs)
            for (CH t : mesh.face_cells(mesh.face_handle(hf)))
                if (meshProps().get<MC_BLOCK>(t) == b)
                {
                    adjBlockCells.emplace_back(t);
                    break;
                }

        VH vI = mesh.from_vertex_handle(heOut);
        VH vJ = mesh.to_vertex_handle(heOut);
        double weightIJ = 0;
        if (xyzTarget)
        {
            Vec3d ptI = mesh.vertex(vI);
            Vec3d ptJ = mesh.vertex(vJ);
            Vec3d vecIJ = ptJ - ptI;
            double length = vecIJ.length();
            vecIJ.normalize();

            for (int plusminus = 0; plusminus < 2; plusminus++)
            {
                HEH heNext = mesh.next_halfedge_in_halfface(heOut, adjPatchHfs[plusminus]);
                Vec3d ptJpm1 = mesh.vertex(mesh.to_vertex_handle(heNext));
                Vec3d vecIJpm1 = ptJpm1 - ptI;
                assert(vecIJpm1.sqrnorm() > 0);
                vecIJpm1.normalize();
                double angle = std::acos(std::max(-1.0, std::min(1.0, vecIJ | vecIJpm1))) / 2.0;
                weightIJ += std::tan(angle);
            }
            weightIJ /= 2.0 * std::max((double)FLT_MIN, length);
            weightIJ = std::max((double)FLT_MIN, weightIJ);
        }
        else
        {
            Vec3Q ptI = meshProps().ref<CHART>(adjBlockCells[0]).at(vI);
            Vec3Q ptJ = meshProps().ref<CHART>(adjBlockCells[0]).at(vJ);
            Vec3d vecIJ = Vec3Q2d(ptJ - ptI);
            double length = vecIJ.length();
            vecIJ.normalize();

            for (int plusminus = 0; plusminus < 2; plusminus++)
            {
                HEH heNext = mesh.next_halfedge_in_halfface(heOut, adjPatchHfs[plusminus]);
                Vec3Q ptJpm1 = meshProps().ref<CHART>(adjBlockCells[plusminus]).at(mesh.to_vertex_handle(heNext));
                Vec3d vecIJpm1 = Vec3Q2d(ptJpm1 - ptI);
                vecIJpm1.normalize();
                double angle = std::acos(std::max(-1.0, std::min(1.0, vecIJ | vecIJpm1))) / 2.0;
                weightIJ += std::tan(angle);
            }
            weightIJ /= 2.0 * std::max((double)FLT_MIN, length);
            weightIJ = std::max((double)FLT_MIN, weightIJ);
        }
        e2meanValWeight[e] = weightIJ;
    }
    return e2meanValWeight;
}

void IGMInitializer::initializeInBlocks()
{
#ifdef UNIFORM_SPACING
    initializeInBlocks3DTutte(false, false);
#else
    initializeInBlocks3DTutte(true, true);
#endif
}

void IGMInitializer::initializeInBlocks3DTutte(bool cotWeights, bool xyzTarget)
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    vector<double> edgeWeights;
    if (cotWeights)
        edgeWeights = calc3DcotangentWeights(xyzTarget);

    for (CH b : mc.cells())
    {
        map<VH, int> vtx2index;
        for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
        {
            for (VH v : mesh.cell_vertices(tet))
            {
                if (meshProps().isInPatch(v))
                    continue;
                int nextIdx = vtx2index.size();
                if (vtx2index.find(v) == vtx2index.end())
                    vtx2index[v] = nextIdx;
            }
        }

        Eigen::SparseMatrix<double, Eigen::RowMajor> A(vtx2index.size(), vtx2index.size());
        Eigen::VectorXd rhsU(Eigen::VectorXd::Zero(vtx2index.size()));
        Eigen::VectorXd rhsV(Eigen::VectorXd::Zero(vtx2index.size()));
        Eigen::VectorXd rhsW(Eigen::VectorXd::Zero(vtx2index.size()));

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(50 * vtx2index.size());

        for (const auto& kv : vtx2index)
        {
            auto vI = kv.first;
            auto i = kv.second;

            double aII = 0;
            for (HEH heOut : mesh.outgoing_halfedges(vI))
            {
                CH tet = *mesh.hec_iter(heOut);
                VH vJ = mesh.to_vertex_handle(heOut);
                double weightIJ = cotWeights ? edgeWeights[mesh.edge_handle(heOut).idx()] : 1.0;
                aII -= weightIJ;

                auto jIt = vtx2index.find(vJ);
                if (jIt != vtx2index.end())
                {
                    // vJ Inner
                    auto j = jIt->second;
                    triplets.push_back(Eigen::Triplet<double>(i, j, weightIJ));
                }
                else
                {
                    // vJ Boundary
                    rhsU(i) -= weightIJ * meshProps().ref<CHART_IGM>(tet).at(vJ)[0].get_d();
                    rhsV(i) -= weightIJ * meshProps().ref<CHART_IGM>(tet).at(vJ)[1].get_d();
                    rhsW(i) -= weightIJ * meshProps().ref<CHART_IGM>(tet).at(vJ)[2].get_d();
                }
            }

            triplets.push_back(Eigen::Triplet<double>(i, i, aII));
        }
        A.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);

        assert(solver.info() == Eigen::Success);

        Eigen::VectorXd solutionU = solver.solve(rhsU).eval();
        Eigen::VectorXd solutionV = solver.solve(rhsV).eval();
        Eigen::VectorXd solutionW = solver.solve(rhsW).eval();

        for (const auto& kv : vtx2index)
        {
            Vec3Q scaledUVW(solutionU[kv.second], solutionV[kv.second], solutionW[kv.second]);
            for (CH tet : mesh.vertex_cells(kv.first))
                meshProps().ref<CHART_IGM>(tet)[kv.first] = scaledUVW;
        }
    }
}

vector<double> IGMInitializer::calc3DcotangentWeights(bool xyzTarget) const
{
    const TetMesh& mesh = meshProps().mesh();
    vector<double> edgeWeights(mesh.n_edges(), DBL_MAX);

    for (EH e : mesh.edges())
    {
        if (edgeWeights[e.idx()] != DBL_MAX || meshProps().isInPatch(e))
            continue;

        double w = 0;

        for (CH tet : mesh.edge_cells(e))
        {
            TetElements elems = getTetElements(tet, e);

            HEH heCB = mesh.opposite_halfedge_handle(elems.heBC);

            HEH heCD = mesh.next_halfedge_in_halfface(elems.heBC, elems.hfBCD);
            HFH hfDCA = mesh.adjacent_halfface_in_cell(elems.hfBCD, heCD);
            HEH heDC = mesh.opposite_halfedge_handle(heCD);
            HEH heAD = mesh.prev_halfedge_in_halfface(heDC, hfDCA);
            HEH heBA = mesh.next_halfedge_in_halfface(heCB, elems.hfCBA);
            HFH hfBDA = mesh.adjacent_halfface_in_cell(elems.hfCBA, heBA);

            double L_kl = 0;
            double theta = 0;

            if (xyzTarget)
            {
                L_kl = (mesh.vertex(mesh.to_vertex_handle(heAD)) - mesh.vertex(mesh.from_vertex_handle(heAD))).length();
                theta = dihedralAngleXYZ(hfDCA, hfBDA);
            }
            else
            {
                L_kl = Vec3Q2d(meshProps().ref<CHART>(tet).at(mesh.to_vertex_handle(heAD))
                               - meshProps().ref<CHART>(tet).at(mesh.from_vertex_handle(heAD)))
                           .length();
                theta = dihedralAngleUVW(hfDCA, hfBDA);
            }

            if (L_kl != 0)
            {
                w += L_kl * std::cos(theta) / std::sin(theta);
            }
            assert(std::isfinite(w));
        }

        w /= 6;
        if (w > 1e10)
            w = 1e10;
        else if (w < 0)
            w = 0;
        edgeWeights[e.idx()] = w;
    }

    return edgeWeights;
}

CH IGMInitializer::traceToEndWithinBlock(const CH& startTet, const HEH& haTrace, set<CH>& tetsVisited) const
{
    auto& tetMesh = meshProps().mesh();
    tetsVisited.clear();

    auto pathHes = mcMeshProps().haHalfedges(haTrace);
    set<VH> haVs;
    for (HEH he : pathHes)
        haVs.insert(tetMesh.from_vertex_handle(he));
    haVs.insert(tetMesh.to_vertex_handle(pathHes.back()));
    bool correctSector = false;

    forVertexNeighbourTetsInBlock(tetMesh.from_vertex_handle(pathHes.front()),
                                  startTet,
                                  [this, &correctSector, eStart = tetMesh.edge_handle(pathHes.front())](const CH& tet2)
                                  {
                                      correctSector = contains(meshProps().mesh().cell_edges(tet2), eStart);
                                      return correctSector;
                                  });
    if (!correctSector)
        return CH();

    tetsVisited.insert(startTet);
    list<CH> tetQ({{startTet}});
    while (!tetQ.empty())
    {
        CH tet = tetQ.front();
        tetQ.pop_front();

        for (HFH hf : tetMesh.cell_halffaces(tet))
        {
            if (meshProps().isBlockBoundary(hf))
                continue;
            bool hasV = false;
            for (VH v : meshProps().get_halfface_vertices(hf))
                if (haVs.count(v) != 0)
                {
                    hasV = true;
                    break;
                }
            if (!hasV)
                continue;
            CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
            if (tetsVisited.count(tetNext) != 0)
                continue;
            tetsVisited.insert(tetNext);
            tetQ.emplace_back(tetNext);
        }
    }
    for (CH tet : tetMesh.halfedge_cells(pathHes.back()))
        if (tetsVisited.count(tet) != 0)
            return tet;
    return CH();
}

template <typename Scalar>
void IGMInitializer::solve2DTutte(const map<VH, int>& vtx2index,
                                  const map<EH, double>& e2weight,
                                  int isoCoord,
                                  Q isoValue,
                                  bool front,
                                  const set<HFH>& hfs)
{
    auto& tetMesh = meshProps().mesh();
    vector<int> nonIsoCoords;
    for (int i = 0; i < 3; i++)
        if (i != isoCoord)
            nonIsoCoords.emplace_back(i);

    using VectorX = Eigen::Matrix<Scalar, -1, 1>;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> A(vtx2index.size(), vtx2index.size());
    VectorX rhsCoord1(VectorX::Zero(vtx2index.size()));
    VectorX rhsCoord2(VectorX::Zero(vtx2index.size()));

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(6 * vtx2index.size());

    for (const auto& kv : vtx2index)
    {
        VH vI = kv.first;
        int i = kv.second;

        Scalar aII = 0;
        for (HEH heOut : tetMesh.outgoing_halfedges(vI))
        {
            vector<HFH> adjPatchHfs;
            for (HFH hf : tetMesh.edge_halffaces(tetMesh.edge_handle(heOut)))
                if (hfs.count(front ? hf : tetMesh.opposite_halfface_handle(hf)))
                    adjPatchHfs.emplace_back(hf);
            assert(adjPatchHfs.size() == 2 || adjPatchHfs.size() == 0);
            if (adjPatchHfs.size() != 2)
                continue;

            CH tet = tetMesh.incident_cell(adjPatchHfs.front());

            VH vJ = tetMesh.to_vertex_handle(heOut);
            Scalar weightIJ = e2weight.empty() ? 1.0 : e2weight.at(tetMesh.edge_handle(heOut));

            aII -= weightIJ;

            auto jIt = vtx2index.find(vJ);
            if (jIt != vtx2index.end())
            {
                // vJ Inner
                int j = jIt->second;
                triplets.push_back(Eigen::Triplet<Scalar>(i, j, weightIJ));
            }
            else
            {
                // vJ Boundary
                if constexpr (std::is_same_v<Scalar, Q>)
                {
                    rhsCoord1(i) -= weightIJ * meshProps().ref<CHART_IGM>(tet).at(vJ)[nonIsoCoords[0]];
                    rhsCoord2(i) -= weightIJ * meshProps().ref<CHART_IGM>(tet).at(vJ)[nonIsoCoords[1]];
                }
                else
                {
                    rhsCoord1(i) -= weightIJ * meshProps().ref<CHART_IGM>(tet).at(vJ)[nonIsoCoords[0]].get_d();
                    rhsCoord2(i) -= weightIJ * meshProps().ref<CHART_IGM>(tet).at(vJ)[nonIsoCoords[1]].get_d();
                }
            }
        }

        triplets.push_back(Eigen::Triplet<Scalar>(i, i, aII));
    }

    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<Scalar>> solver(A);

    if (solver.info() != Eigen::Success)
        throw std::logic_error("Could not solve tutte");

    VectorX solutionCoord1 = solver.solve(rhsCoord1).eval();
    VectorX solutionCoord2 = solver.solve(rhsCoord2).eval();

    for (const auto& kv : vtx2index)
    {
        Vec3Q scaledUVW(0, 0, 0);
        scaledUVW[isoCoord] = isoValue;
        scaledUVW[nonIsoCoords[0]] = solutionCoord1[kv.second];
        scaledUVW[nonIsoCoords[1]] = solutionCoord2[kv.second];
        CH tet;
        for (HFH hf : tetMesh.vertex_halffaces(kv.first))
            if (hfs.count(front ? hf : tetMesh.opposite_halfface_handle(hf)) != 0)
            {
                tet = tetMesh.incident_cell(hf);
                break;
            }
        assert(tet.is_valid());

        auto tet2trans = determineTransitionsAroundVertex<TRANSITION_IGM>(kv.first, tet);
        for (auto& kv2 : tet2trans)
            meshProps().ref<CHART_IGM>(kv2.first)[kv.first] = kv2.second.apply(scaledUVW);
    }
}

template void IGMInitializer::solve2DTutte<Q>(const map<VH, int>& vtx2index,
                                              const map<EH, double>& e2weight,
                                              int isoCoord,
                                              Q isoValue,
                                              bool front,
                                              const set<HFH>& hfs);
template void IGMInitializer::solve2DTutte<double>(const map<VH, int>& vtx2index,
                                                   const map<EH, double>& e2weight,
                                                   int isoCoord,
                                                   Q isoValue,
                                                   bool front,
                                                   const set<HFH>& hfs);


} // namespace c4hex
