#include "C4Hex/Algorithm/PathRouter.hpp"
#include "C4Hex/Algorithm/SurfaceRouter.hpp"

#include <fstream>
#include <queue>


namespace c4hex
{

namespace
{

struct VtxHeuristic
{
    VtxHeuristic(const VH& _v, double _dist, double _h) : v(_v), dist(_dist), heuristic(_h)
    {
    }

    VH v;
    double dist;
    double heuristic;
};

template <typename HEURISTIC>
struct LeastHeuristicComp
{
    bool operator()(const HEURISTIC& p1, const HEURISTIC& p2) const
    {
        return p1.heuristic > p2.heuristic;
    }
};
} // namespace

PathRouter::PathRouter(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

PathRouter::RetCode PathRouter::reroutePathThroughSurface(list<HEH>& path,
                                                          set<FH>& surface,
                                                          set<FH>& forbiddenFs,
                                                          set<EH>& forbiddenEs,
                                                          set<VH>& forbiddenVs,
                                                          set<HFH>& hfsTransferred)
{
    _forbiddenFs = forbiddenFs;
    _forbiddenEs = forbiddenEs;
    _forbiddenVs = forbiddenVs;

    auto& tetMesh = meshProps().mesh();
    VH vFrom = tetMesh.from_vertex_handle(path.front());
    VH vTo = tetMesh.to_vertex_handle(path.back());

    assert(!path.empty());
    assert(!surface.empty());
    assert(_forbiddenVs.count(vTo) == 0);
    assert(_forbiddenVs.count(vFrom) == 0);
    if (vFrom == vTo)
        LOG(WARNING) << "Trying to reroute path between identical from/to vertex";

    auto pathRerouted = path;
    {
        TemporaryPropAllocator<TetMeshProps, CHILD_HALFEDGES> propGuard(meshProps());
        auto ret = shortestPathThroughSurface(vFrom, vTo, pathRerouted, surface, forbiddenFs, forbiddenEs, forbiddenVs);

        // Replace deleted halfedges in path
        for (auto it = path.begin(); it != path.end();)
        {
            HEH he = *it;
            if (meshProps().mesh().is_deleted(he))
            {
                path.erase(it++);
                auto children = meshProps().get<CHILD_HALFEDGES>(he);
                it = path.insert(it, children.begin(), children.end());
            }
            else
                it++;
        }

        if (ret != SUCCESS)
        {
            LOG(ERROR) << "shortest path failed from vertex " << vFrom << " to " << vTo;
            forbiddenFs = _forbiddenFs;
            forbiddenEs = _forbiddenEs;
            forbiddenVs = _forbiddenVs;
            return ret;
        }
    }

    // Determine where pathRerouted branches off path to then determine the spanned surface portions between branches
    list<list<HEH>> branchesPath;
    list<list<HEH>> branchesPathRerouted;
    auto ret = determineBranches(path, pathRerouted, branchesPath, branchesPathRerouted);
    if (ret != SUCCESS)
        return ret;

    ret = getEnclosedHalffaces(branchesPath, branchesPathRerouted, surface, hfsTransferred);
    if (ret != SUCCESS)
        return ret;

    path = pathRerouted;

    forbiddenFs = _forbiddenFs;
    forbiddenEs = _forbiddenEs;
    forbiddenVs = _forbiddenVs;

    return SUCCESS;
}

PathRouter::RetCode PathRouter::shortestPathThroughSurface(const VH& vFrom,
                                                           const VH& vTo,
                                                           list<HEH>& path,
                                                           set<FH>& surface,
                                                           set<FH>& forbiddenFs,
                                                           set<EH>& forbiddenEs,
                                                           set<VH>& forbiddenVs)
{
    _forbiddenFs = forbiddenFs;
    _forbiddenEs = forbiddenEs;
    _forbiddenVs = forbiddenVs;
    if (vFrom == vTo)
        LOG(WARNING) << "WARNING: Trying to reroute path between identical from/to vertex, this probably does not work";
    path.clear();

    auto& tetMesh = meshProps().mesh();

    set<EH> allowedEs;
    set<VH> allowedVs;
    for (FH f : surface)
    {
        for (EH e : tetMesh.face_edges(f))
            if (forbiddenEs.count(e) == 0)
                allowedEs.insert(e);
        for (VH v : tetMesh.face_vertices(f))
            if (forbiddenVs.count(v) == 0)
                allowedVs.insert(v);
    }
    auto ret = aStarShortestPath(vFrom, vTo, path, allowedEs, allowedVs);
    if (ret != SUCCESS)
    {
        ret = refineSurfaceToAllowReroute(surface, vFrom, vTo);
        if (ret != SUCCESS)
            return ret;
        allowedEs.clear();
        allowedVs.clear();
        for (FH f : surface)
        {
            for (EH e : tetMesh.face_edges(f))
                if (forbiddenEs.count(e) == 0)
                    allowedEs.insert(e);
            for (VH v : tetMesh.face_vertices(f))
                if (forbiddenVs.count(v) == 0)
                    allowedVs.insert(v);
        }

        ret = aStarShortestPath(vFrom, vTo, path, allowedEs, allowedVs);
        if (ret != SUCCESS)
            LOG(WARNING) << "shortest path failed from vertex " << vFrom << " to " << vTo;
    }

    forbiddenFs = _forbiddenFs;
    forbiddenEs = _forbiddenEs;
    forbiddenVs = _forbiddenVs;
    return ret;
}

PathRouter::RetCode PathRouter::reroutePathThroughVolume(list<HEH>& path,
                                                         set<CH>& volume,
                                                         set<FH>& forbiddenFs,
                                                         set<EH>& forbiddenEs,
                                                         set<VH>& forbiddenVs,
                                                         set<HFH>* hfsTransferred)
{
    _forbiddenFs = forbiddenFs;
    _forbiddenEs = forbiddenEs;
    _forbiddenVs = forbiddenVs;

    auto& tetMesh = meshProps().mesh();
    VH vFrom = tetMesh.from_vertex_handle(path.front());
    VH vTo = tetMesh.to_vertex_handle(path.back());

    if (vFrom == vTo)
        LOG(WARNING) << "WARNING: Trying to reroute path between identical from/to vertex, this probably does not work";
    assert(forbiddenVs.count(vTo) == 0);
    assert(forbiddenVs.count(vFrom) == 0);

    auto pathRerouted = path;
    {
        TemporaryPropAllocator<TetMeshProps, CHILD_HALFEDGES> propGuard(meshProps());
        auto ret = shortestPathThroughVolume(vFrom, vTo, pathRerouted, volume, forbiddenFs, forbiddenEs, forbiddenVs);

        // Replace deleted halfedges in path
        for (auto it = path.begin(); it != path.end();)
        {
            HEH he = *it;
            if (meshProps().mesh().is_deleted(he))
            {
                path.erase(it++);
                auto children = meshProps().get<CHILD_HALFEDGES>(he);
                it = path.insert(it, children.begin(), children.end());
            }
            else
                it++;
        }

        if (ret != SUCCESS)
        {
            forbiddenFs = _forbiddenFs;
            forbiddenEs = _forbiddenEs;
            forbiddenVs = _forbiddenVs;
            return ret;
        }
    }

    if (hfsTransferred != nullptr)
    {
        // Determine where pathRerouted branches off path
        list<list<HEH>> branchesPath;
        list<list<HEH>> branchesPathRerouted;
        auto ret = determineBranches(path, pathRerouted, branchesPath, branchesPathRerouted);
        if (ret != SUCCESS)
        {
            forbiddenFs = _forbiddenFs;
            forbiddenEs = _forbiddenEs;
            forbiddenVs = _forbiddenVs;
            return ret;
        }
        ret = getEnclosedHalffaces(branchesPath, branchesPathRerouted, volume, *hfsTransferred);
        if (ret != SUCCESS)
        {
            forbiddenFs = _forbiddenFs;
            forbiddenEs = _forbiddenEs;
            forbiddenVs = _forbiddenVs;
            return ret;
        }
    }

    path = pathRerouted;

    forbiddenFs = _forbiddenFs;
    forbiddenEs = _forbiddenEs;
    forbiddenVs = _forbiddenVs;

    return SUCCESS;
}

PathRouter::RetCode PathRouter::shortestPathThroughVolume(const VH& vFrom,
                                                          const VH& vTo,
                                                          list<HEH>& path,
                                                          set<CH>& volume,
                                                          set<FH>& forbiddenFs,
                                                          set<EH>& forbiddenEs,
                                                          set<VH>& forbiddenVs)
{
    _forbiddenFs = forbiddenFs;
    _forbiddenEs = forbiddenEs;
    _forbiddenVs = forbiddenVs;
    if (vFrom == vTo)
        LOG(WARNING) << "WARNING: Trying to reroute path between identical from/to vertex, this probably does not work";
    path.clear();

    auto& tetMesh = meshProps().mesh();

    set<CH> confinedVolume;

    Vec3Q bboxMin;
    Vec3Q bboxMax;

    CH seedTet = findSomeOf(tetMesh.vertex_cells(vFrom), volume);
    assert(seedTet.is_valid());

    {
        set<CH> tetVisited({seedTet});
        list<pair<CH, Transition>> tetQ({{seedTet, Transition()}});
        Vec3Q uvwFrom = meshProps().ref<CHART>(seedTet).at(vFrom);
        Vec3Q uvwTo;
        while (!tetQ.empty())
        {
            auto tetTrans = tetQ.front();
            tetQ.pop_front();

            if (contains(tetMesh.tet_vertices(tetTrans.first), vTo))
            {
                uvwTo = tetTrans.second.invert().apply(meshProps().ref<CHART>(tetTrans.first).at(vTo));
                break;
            }

            for (HFH hf : tetMesh.cell_halffaces(tetTrans.first))
            {
                CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                if (tetNext.is_valid() && tetVisited.count(tetNext) == 0 && volume.count(tetNext) != 0)
                {
                    tetVisited.insert(tetNext);
                    Transition transNext = tetTrans.second.chain(meshProps().hfTransition<TRANSITION>(hf));
                    tetQ.push_back({tetNext, transNext});
                }
            }
        }
        bboxMin = uvwFrom;
        bboxMax = bboxMin;
        for (int i = 0; i < 3; i++)
            if (bboxMin[i] > uvwTo[i])
                bboxMin[i] = uvwTo[i];
            else if (bboxMax[i] < uvwTo[i])
                bboxMax[i] = uvwTo[i];
        Vec3Q diff = (bboxMax - bboxMin) * Q(1, 20) + Vec3Q(Q(1, 1000), Q(1, 1000), Q(1, 1000));
        bboxMin -= diff;
        bboxMax += diff;
    }

    {
        confinedVolume.insert(seedTet);
        set<CH> tetVisited({seedTet});
        list<pair<CH, Transition>> tetQ({{seedTet, Transition()}});
        while (!tetQ.empty())
        {
            auto tetTrans = tetQ.front();
            tetQ.pop_front();

            for (HFH hf : tetMesh.cell_halffaces(tetTrans.first))
            {
                CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                if (tetNext.is_valid() && tetVisited.count(tetNext) == 0 && volume.count(tetNext) != 0)
                {
                    tetVisited.insert(tetNext);
                    Transition transNext = tetTrans.second.chain(meshProps().hfTransition<TRANSITION>(hf));
                    bool withinBbox = containsMatching(
                        tetMesh.tet_vertices(tetNext),
                        [&, this](const VH& v)
                        {
                            Vec3Q uvw = transNext.invert().apply(meshProps().ref<CHART>(tetNext).at(v));
                            return uvw[0] <= bboxMax[0] && uvw[0] >= bboxMin[0] && uvw[1] <= bboxMax[1]
                                   && uvw[1] >= bboxMin[1] && uvw[2] <= bboxMax[2] && uvw[2] >= bboxMin[2];
                        });
                    if (withinBbox)
                    {
                        confinedVolume.insert(tetNext);
                        tetQ.push_back({tetNext, transNext});
                    }
                }
            }
        }
    }

    for (bool confined : {true, false})
    {
        if (confined && confinedVolume.empty())
            continue;
        if (!confined)
        {
            auto ret = refineVolumeToAllowReroute(volume, vFrom, vTo);
            if (ret != SUCCESS)
            {
                forbiddenFs = _forbiddenFs;
                forbiddenEs = _forbiddenEs;
                forbiddenVs = _forbiddenVs;
                return ret;
            }
        }
        if (confined)
        {
            // Update main volume
            list<CH> children(volume.begin(), volume.end());
            for (auto it = children.begin(); it != children.end();)
            {
                CH tet = *it;
                if (tetMesh.is_deleted(tet))
                {
                    for (CH child : meshProps().get<CHILD_CELLS>(tet))
                        children.push_back(child);
                    children.erase(it++);
                }
                else
                    it++;
            }
            volume = {children.begin(), children.end()};
        }

        assert(forbiddenVs.count(vTo) == 0);
        assert(forbiddenVs.count(vFrom) == 0);

        set<EH> allowedEs;
        set<VH> allowedVs;
        for (CH tet : confined ? confinedVolume : volume)
        {
            for (EH e : tetMesh.cell_edges(tet))
                if (forbiddenEs.count(e) == 0)
                    allowedEs.insert(e);
            for (VH v : tetMesh.cell_vertices(tet))
                if (forbiddenVs.count(v) == 0)
                    allowedVs.insert(v);
        }
        auto ret = aStarShortestPath(vFrom, vTo, path, allowedEs, allowedVs);
        if (ret != SUCCESS)
        {
            if (!confined)
            {
                forbiddenFs = _forbiddenFs;
                forbiddenEs = _forbiddenEs;
                forbiddenVs = _forbiddenVs;
                return ret;
            }
        }
        else
        {
            forbiddenFs = _forbiddenFs;
            forbiddenEs = _forbiddenEs;
            forbiddenVs = _forbiddenVs;
            return ret;
        }
    }

    forbiddenFs = _forbiddenFs;
    forbiddenEs = _forbiddenEs;
    forbiddenVs = _forbiddenVs;
    return SUCCESS;
}

PathRouter::RetCode PathRouter::determineBranches(const list<HEH>& path,
                                                  const list<HEH>& pathRerouted,
                                                  list<list<HEH>>& branchesPath,
                                                  list<list<HEH>>& branchesPathRerouted)
{
    auto& tetMesh = meshProps().mesh();

    set<VH> pathVs;
    for (HEH he : path)
        pathVs.insert(tetMesh.from_vertex_handle(he));
    pathVs.insert(tetMesh.to_vertex_handle(path.back()));

    auto itPath = path.begin();
    auto itPathRerouted = pathRerouted.begin();

    set<VH> vsPathPassed;

    bool branched = false;
    while (itPathRerouted != pathRerouted.end())
    {
        if (!branched)
        {
            if (*itPathRerouted == *itPath)
                itPath++;
            else
            {
                branched = true;
                branchesPath.emplace_back();
                branchesPathRerouted.emplace_back();
            }
        }

        if (branched)
        {
            branchesPathRerouted.back().emplace_back(*itPathRerouted);
            VH vCurr = tetMesh.to_vertex_handle(*itPathRerouted);
            if (vsPathPassed.count(vCurr) != 0 && vCurr != tetMesh.from_vertex_handle(pathRerouted.front()))
                return NOT_CONNECTED;
            if (pathVs.count(vCurr) != 0)
            {
                branched = false;
                while (itPath != path.end() && tetMesh.from_vertex_handle(*itPath) != vCurr)
                {
                    vsPathPassed.insert(tetMesh.from_vertex_handle(*itPath));
                    branchesPath.back().emplace_back(*(itPath++));
                }
            }
        }
        itPathRerouted++;
    }
    assert(branchesPath.size() == branchesPathRerouted.size());
    return SUCCESS;
}
PathRouter::RetCode PathRouter::getEnclosedHalffaces(const list<list<HEH>>& branchesPath,
                                                     const list<list<HEH>>& branchesPathRerouted,
                                                     const set<FH>& surface,
                                                     set<HFH>& hfsEnclosed) const
{
    auto& tetMesh = meshProps().mesh();

    for (auto itPath = branchesPath.begin(), itRerouted = branchesPathRerouted.begin();
         itPath != branchesPath.end() && itRerouted != branchesPathRerouted.end();
         itPath++, itRerouted++)
    {
        auto& branchPath = *itPath;
        auto& branchRerouted = *itRerouted;
        assert(tetMesh.to_vertex_handle(branchPath.back()) == tetMesh.to_vertex_handle(branchRerouted.back()));
        assert(tetMesh.from_vertex_handle(branchPath.front()) == tetMesh.from_vertex_handle(branchRerouted.front()));
        set<HEH> boundaryHes(branchPath.begin(), branchPath.end());
        set<EH> boundaryEs;
        for (HEH he : branchRerouted)
            boundaryHes.insert(tetMesh.opposite_halfedge_handle(he));
        for (HEH he : boundaryHes)
            boundaryEs.insert(tetMesh.edge_handle(he));

        set<HFH> hfsEnclosedInBranch;
        set<HFH> hfsEnclosedInBranchAlt;

        list<HFH> hfQ;
        list<HFH> hfQAlt;
        for (HFH hf : tetMesh.halfedge_halffaces(*boundaryHes.begin()))
        {
            if (surface.count(tetMesh.face_handle(hf)) != 0)
            {
                if (hfQ.empty())
                {
                    hfQ.push_back(hf);
                    hfsEnclosedInBranch.insert(hf);
                }
                else
                {
                    hfQAlt.push_back(hf);
                    hfsEnclosedInBranchAlt.insert(hf);
                }
            }
        }
        assert(hfQ.size() == 1 && hfQAlt.size() < 2);
        bool valid = !hfQ.empty();
        bool validAlt = !hfQAlt.empty();
        bool alt = false;
        while (!(valid && hfQ.empty()) && !(validAlt && hfQAlt.empty()))
        {
            if (!valid || hfQ.empty())
                alt = true;
            else if (!validAlt || hfQAlt.empty())
                alt = false;
            else
                alt = !alt;
            auto& QQ = alt ? hfQAlt : hfQ;
            auto& hfs = alt ? hfsEnclosedInBranchAlt : hfsEnclosedInBranch;
            bool& val = alt ? validAlt : valid;

            HFH hf = QQ.front();
            QQ.pop_front();
            for (HEH he : tetMesh.halfface_halfedges(hf))
            {
                if (boundaryHes.count(he) != 0)
                    continue;
                he = tetMesh.opposite_halfedge_handle(he);
                val = false;
                // if we run into the boundary of "surface" while floodfilling, we are floodfilling the wrong side
                // mark this by -> val = false
                for (HFH hfNext : tetMesh.halfedge_halffaces(he))
                {
                    if (tetMesh.opposite_halfface_handle(hfNext) == hf
                        || surface.count(tetMesh.face_handle(hfNext)) == 0)
                        continue;
                    val = true;
                    if (hfs.count(hfNext) != 0)
                        continue;
                    QQ.push_back(hfNext);
                    hfs.insert(hfNext);
                }
                if (!val)
                    break;
            }
        }
        if (!valid && !validAlt)
        {
            LOG(ERROR) << "Branch does not enclose disk-topology surface";
            return NOT_CONNECTED;
        }
        assert(valid || validAlt);
        assert((valid && hfQ.empty()) != (validAlt && hfQAlt.empty()));

        if (valid && hfQ.empty())
            hfsEnclosed.insert(hfsEnclosedInBranch.begin(), hfsEnclosedInBranch.end());
        else
            hfsEnclosed.insert(hfsEnclosedInBranchAlt.begin(), hfsEnclosedInBranchAlt.end());
    }
    return SUCCESS;
}

PathRouter::RetCode PathRouter::getEnclosedHalffaces(const list<list<HEH>>& branchesPath,
                                                     const list<list<HEH>>& branchesPathRerouted,
                                                     set<CH>& volume,
                                                     set<HFH>& hfsEnclosed)
{
    auto& tetMesh = meshProps().mesh();
    SurfaceRouter sr(meshProps());
    hfsEnclosed.clear();
    for (auto itPath = branchesPath.begin(), itRerouted = branchesPathRerouted.begin();
            itPath != branchesPath.end() && itRerouted != branchesPathRerouted.end();
            itPath++, itRerouted++)
    {
        auto& branchPath = *itPath;
        auto& branchRerouted = *itRerouted;
        assert(tetMesh.to_vertex_handle(branchPath.back()) == tetMesh.to_vertex_handle(branchRerouted.back()));
        assert(tetMesh.from_vertex_handle(branchPath.front())
                == tetMesh.from_vertex_handle(branchRerouted.front()));
        set<HEH> boundaryHes(branchPath.begin(), branchPath.end());
        for (HEH he : branchRerouted)
            boundaryHes.insert(tetMesh.opposite_halfedge_handle(he));
        set<EH> boundaryEs;
        set<VH> boundaryVs;
        for (HEH he : boundaryHes)
        {
            boundaryEs.insert(tetMesh.edge_handle(he));
            boundaryVs.insert(tetMesh.from_vertex_handle(he));
        }

        set<VH> forbiddenVsMinusBoundary;
        set<EH> forbiddenEsMinusBoundary;
        for (VH v : _forbiddenVs)
            if (boundaryVs.count(v) == 0)
                forbiddenVsMinusBoundary.insert(v);
        for (EH e : _forbiddenEs)
            if (boundaryEs.count(e) == 0)
                forbiddenEsMinusBoundary.insert(e);

        set<HFH> hfsEnclosedInBranch;
        if (sr.calcMinimalSurfaceByLP(volume,
                                      _forbiddenFs,
                                      forbiddenEsMinusBoundary,
                                      forbiddenVsMinusBoundary,
                                      boundaryEs,
                                      boundaryHes,
                                      boundaryVs,
                                      hfsEnclosedInBranch)
            != SurfaceRouter::SUCCESS)
            return NOT_CONNECTED;
        hfsEnclosed.insert(hfsEnclosedInBranch.begin(), hfsEnclosedInBranch.end());
    }
    return SUCCESS;
}

#define SPLIT_IF_ALL_VS_FORBIDDEN(FORBIDDENSET, SPLITSET, ELEMENT, VERTEXRANGE)                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        if (FORBIDDENSET.count(ELEMENT) != 0 || SPLITSET.count(ELEMENT) != 0)                                          \
            continue;                                                                                                  \
        bool allVsForbidden = true;                                                                                    \
        for (VH v : tetMesh.VERTEXRANGE(ELEMENT))                                                                      \
            if (vFrom != v && vTo != v && _forbiddenVs.count(v) == 0)                                                  \
                allVsForbidden = false;                                                                                \
        if (allVsForbidden)                                                                                            \
            SPLITSET.insert(ELEMENT);                                                                                  \
    } while (0)

PathRouter::RetCode PathRouter::refineSurfaceToAllowReroute(set<FH>& surface, const VH& vFrom, const VH& vTo)
{
    auto& tetMesh = meshProps().mesh();

    TemporaryPropAllocator<TetMeshProps, CHILD_FACES> propGuard(meshProps());

    set<FH> splitFs;
    set<EH> splitEs;
    for (FH f : surface)
        for (EH e : tetMesh.face_edges(f))
            SPLIT_IF_ALL_VS_FORBIDDEN(_forbiddenEs, splitEs, e, edge_vertices);

    for (FH f : surface)
        if (!containsSomeOf(tetMesh.face_edges(f), splitEs))
            SPLIT_IF_ALL_VS_FORBIDDEN(_forbiddenFs, splitFs, f, face_vertices);

    for (FH f : splitFs)
        splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});

    for (EH e : splitEs)
        if (!containsSomeOf(tetMesh.edge_faces(e), _forbiddenFs))
            splitHalfEdge(tetMesh.halfedge_handle(e, 0), *tetMesh.ec_iter(e), 0.5);

    list<FH> children(surface.begin(), surface.end());
    for (auto it = children.begin(); it != children.end();)
    {
        FH f = *it;
        if (tetMesh.is_deleted(f))
        {
            for (FH child : meshProps().get<CHILD_FACES>(f))
                children.push_back(child);
            children.erase(it++);
        }
        else
            it++;
    }
    surface = {children.begin(), children.end()};

    return SUCCESS;
}

PathRouter::RetCode PathRouter::refineVolumeToAllowReroute(set<CH>& volume, const VH& vFrom, const VH& vTo)
{
    auto& tetMesh = meshProps().mesh();

    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());

    set<CH> dummy;
    set<CH> splitTets;
    set<FH> splitFs;
    set<EH> splitEs;
    for (CH tet : volume)
        for (EH e : tetMesh.cell_edges(tet))
            SPLIT_IF_ALL_VS_FORBIDDEN(_forbiddenEs, splitEs, e, edge_vertices);
    for (CH tet : volume)
        for (FH f : tetMesh.cell_faces(tet))
            if (!containsSomeOf(tetMesh.face_edges(f), splitEs))
                SPLIT_IF_ALL_VS_FORBIDDEN(_forbiddenFs, splitFs, f, face_vertices);
    for (CH tet : volume)
        if (!containsSomeOf(tetMesh.cell_edges(tet), splitEs) && !containsSomeOf(tetMesh.cell_faces(tet), splitFs))
            SPLIT_IF_ALL_VS_FORBIDDEN(dummy, splitTets, tet, tet_vertices);

    for (CH tet : splitTets)
        splitTet(tet, {Q(1, 4), Q(1, 4), Q(1, 4), Q(1, 4)});

    for (FH f : splitFs)
        splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});

    for (EH e : splitEs)
        if (!containsMatching(tetMesh.edge_faces(e), [this](const FH& f) { return _forbiddenFs.count(f) != 0; }))
            splitHalfEdge(tetMesh.halfedge_handle(e, 0), *tetMesh.ec_iter(e), 0.5);

    list<CH> children(volume.begin(), volume.end());
    for (auto it = children.begin(); it != children.end();)
    {
        CH tet = *it;
        if (tetMesh.is_deleted(tet))
        {
            for (CH child : meshProps().get<CHILD_CELLS>(tet))
                children.push_back(child);
            children.erase(it++);
        }
        else
            it++;
    }
    volume = {children.begin(), children.end()};

    return SUCCESS;
}

PathRouter::RetCode
PathRouter::aStarShortestPath(const VH& vFrom, const VH& vTo, list<HEH>& path, set<EH>& allowedEs, set<VH>& allowedVs)
{
    using VtxQueue = std::priority_queue<VtxHeuristic, std::deque<VtxHeuristic>, LeastHeuristicComp<VtxHeuristic>>;
    auto& tetMesh = meshProps().mesh();

#ifdef MINIMIZE_XYZ
    Vec3d XYZto = tetMesh.vertex(vTo);
#endif
    VtxQueue vQ;
    vQ.push({VtxHeuristic(vFrom,
                          0,
#ifdef MINIMIZE_XYZ
                          (XYZto - tetMesh.vertex(vFrom)).length()
#else
                          0
#endif
                              )});

    map<VH, std::pair<double, HEH>> v2minDistAndHe;
    for (VH v : allowedVs)
        v2minDistAndHe[v] = {DBL_MAX, HEH()};
    set<EH> esVisited;
    set<VH> vsExpanded;
    while (!vQ.empty())
    {
        auto minHeuristic = vQ.top();
        vQ.pop();
        VH v = minHeuristic.v;
        if (v2minDistAndHe.at(v).second.is_valid())
        {
            vsExpanded.insert(v);
            if (v == vTo)
                break;
        }

        double dist = minHeuristic.dist;

        for (HEH he : tetMesh.outgoing_halfedges(v))
        {
            VH vNext = tetMesh.to_vertex_handle(he);
            auto itMinDist = v2minDistAndHe.find(vNext);
            if (itMinDist == v2minDistAndHe.end() || vsExpanded.find(vNext) != vsExpanded.end()
                || allowedEs.count(tetMesh.edge_handle(he)) == 0 || esVisited.count(tetMesh.edge_handle(he)) != 0)
                continue;

            double heLength =
#ifdef MINIMIZE_XYZ
                tetMesh.length(he);
#else
                edgeLengthUVW<CHART>(tetMesh.edge_handle(he));
#endif
            double distNext = dist + heLength;

            if (itMinDist->second.first <= distNext)
                continue;
            esVisited.insert(tetMesh.edge_handle(he));

            itMinDist->second = {distNext, he};
            vQ.push({vNext,
                     distNext,
                     distNext +
#ifdef MINIMIZE_XYZ
                         (XYZto - tetMesh.vertex(vNext)).length()
#else
                         0
#endif
            });
        }
    }
    // assert(vsExpanded.find(vTo) != vsExpanded.end());

    if (vsExpanded.find(vTo) == vsExpanded.end())
        return NOT_CONNECTED;

    VH vCurr = vTo;
    HEH he = v2minDistAndHe.at(vTo).second;
    do
    {
        assert(tetMesh.to_vertex_handle(he) == vCurr);
        path.emplace_front(he);
        vCurr = tetMesh.from_vertex_handle(he);
        he = v2minDistAndHe.at(vCurr).second;
    } while (vCurr != vFrom);
    assert(tetMesh.from_vertex_handle(path.front()) == vFrom);
    assert(!he.is_valid());

    return SUCCESS;
}

PathRouter::RetCode PathRouter::reroutePathThroughPatch(const FH& p, list<HEH>& path, set<HFH>& hfsTransferred)
{
    _forbiddenFs.clear();
    _forbiddenVs.clear();
    _forbiddenEs.clear();

    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    VH vFrom = tetMesh.from_vertex_handle(path.front());
    VH vTo = tetMesh.to_vertex_handle(path.back());
    if (vFrom == vTo)
        LOG(WARNING) << "WARNING: Trying to reroute path between identical from/to vertex, this probably does not work";
    _p = p;
    auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);

    set<VH> forbiddenVs;
    set<EH> forbiddenEs;
    set<FH> forbiddenFs;
    (void)mcMesh;
    // First try minimal path
    {
        assert(!path.empty());
        assert(!pHfs.empty());
        auto pathRerouted = path;
        // Determine shift/flip direction
        bool forbiddenFront = meshProps().isInArc(path.front());
        _inFirstHp = pHfs.count(findMatching(
                         tetMesh.halfedge_halffaces(forbiddenFront ? path.front()
                                                                   : tetMesh.opposite_halfedge_handle(path.back())),
                         [&, this](const HFH& hf) { return meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)) == p; }))
                     != 0;

        for (EH a : mcMesh.face_edges(p))
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
            {
                forbiddenEs.insert(tetMesh.edge_handle(he));
                for (VH v : tetMesh.halfedge_vertices(he))
                    if (v != vFrom && v != vTo)
                        forbiddenVs.insert(v);
            }
        set<FH> surface;
        for (HFH hf : pHfs)
            surface.insert(tetMesh.face_handle(hf));
        auto ret = shortestPathThroughSurface(vFrom, vTo, pathRerouted, surface, forbiddenFs, forbiddenEs, forbiddenVs);
        if (ret == SUCCESS)
        {
            // LOG(INFO) << "succeeded!";
            set<EH> boundary;
            for (HEH he : path)
                boundary.insert(tetMesh.edge_handle(he));
            for (HEH he : pathRerouted)
                boundary.insert(tetMesh.edge_handle(he));

            HFH hfSeed = findMatching(tetMesh.halfedge_halffaces(forbiddenFront ? path.front()
                                                                    : tetMesh.opposite_halfedge_handle(path.back())),

                         [&, this](const HFH& hf) { return meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)) == p; }
                                                                    );
            _inFirstHp = pHfs.find(hfSeed) != pHfs.end();
            if (!_inFirstHp)
                hfSeed = tetMesh.opposite_halfface_handle(hfSeed);

            list<HFH> hfQ({hfSeed});
            vector<bool> hfVisited(tetMesh.n_halffaces(), false);
            hfVisited[hfSeed.idx()] = true;
            while (!hfQ.empty())
            {
                HFH hf = hfQ.back();
                hfQ.pop_back();

                hfsTransferred.insert(hf);

                for (HEH he : tetMesh.halfface_halfedges(hf))
                {
                    EH e = tetMesh.edge_handle(he);
                    // Do not spread beyond boundary
                    if (boundary.count(e) > 0)
                        continue;

                    for (HFH hfNext : tetMesh.halfedge_halffaces(he))
                    {
                        if (hfNext == hf)
                            continue;
                        HFH hfAdj = tetMesh.opposite_halfface_handle(hfNext);
                        if (!hfVisited[hfAdj.idx()] && pHfs.count(hfAdj) > 0)
                        {
                            hfVisited[hfAdj.idx()] = true;
                            hfQ.push_back(hfAdj);
                        }
                    }
                }
            }
            path = pathRerouted;
            return SUCCESS;
        }
    }

    // Backup: incremental rerouting
    {
        LOG(INFO) << "Backup: Incrementally Rerouting path through patch " << p << " incident on blocks "
                  << mcMeshProps().mesh().incident_cell(mcMeshProps().mesh().halfface_handle(p, 0)) << " "
                  << mcMeshProps().mesh().incident_cell(mcMeshProps().mesh().halfface_handle(p, 1));

        _forbiddenFs = forbiddenFs;
        _forbiddenEs = forbiddenEs;
        _forbiddenVs = forbiddenVs;
        {
            auto ret = refineToAllowReroute(p, vFrom, vTo);
            if (ret != SUCCESS)
                return ret;
        }

        bool invert = _forbiddenEs.find(tetMesh.edge_handle(path.front())) != _forbiddenEs.end();
        assert(invert || _forbiddenEs.find(tetMesh.edge_handle(path.back())) != _forbiddenEs.end());
        list<HEH> pathRerouted = path;
        if (invert)
        {
            pathRerouted.reverse();
            for (auto& he : pathRerouted)
                he = tetMesh.opposite_halfedge_handle(he);
        }

        set<HEH> currentHes(pathRerouted.begin(), pathRerouted.end());
        set<VH> currentVs;
        int nForbiddenEs = 0;
        int nForbiddenVs = 0;
        {
            auto ret = gatherCurrentElements(currentHes, currentVs, nForbiddenEs, nForbiddenVs);
            if (ret != SUCCESS)
                return ret;
        }

        // Determine shift/flip direction
        _inFirstHp =  pHfs.count(findMatching(
                         tetMesh.halfedge_halffaces(pathRerouted.front()),
                         [&, this](const HFH& hf) { return meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)) == p; }))
                     != 0;

        bool change = true;

        while (change)
        {
            change = false;
            set<VH> vsPassed({tetMesh.from_vertex_handle(pathRerouted.front())});
            for (auto it = pathRerouted.begin(); it != pathRerouted.end();)
            {

                HEH heCurrent = *it;
                HFH transferredHf = getTransferredHf(heCurrent);
                // To be filled later
                list<HEH> removedHes;
                list<HEH> addedHes;
                set<VH> removedVs;
                set<VH> addedVs;

                bool doReroute = transferredHf.is_valid();
                if (doReroute)
                {
                    bool force = incidentOnForbiddenElement(heCurrent);

                    auto ret = exchangedElements(transferredHf, currentHes, removedHes, addedHes, removedVs, addedVs);
                    if (ret != SUCCESS)
                        return ret;

                    doReroute = force || rerouteRecommended(heCurrent, removedHes, addedHes, removedVs, addedVs);
                }

                if (doReroute)
                {
                    reroute(pathRerouted,
                            it,
                            hfsTransferred,
                            currentVs,
                            currentHes,
                            vsPassed,
                            removedHes,
                            addedHes,
                            removedVs,
                            addedVs);
                }
                else
                {
                    vsPassed.insert(tetMesh.to_vertex_handle(*it));
                    it++;
                }
            }
        }

        if (invert)
        {
            pathRerouted.reverse();
            for (auto& he : pathRerouted)
                he = tetMesh.opposite_halfedge_handle(he);
        }

        path = pathRerouted;
    }

    return SUCCESS;
}

PathRouter::RetCode PathRouter::gatherCurrentElements(const set<HEH>& currentHes,
                                                      set<VH>& currentVs,
                                                      int& nForbiddenEs,
                                                      int& nForbiddenVs) const
{
    auto& tetMesh = meshProps().mesh();

    nForbiddenEs = 0;
    nForbiddenVs = 0;

    for (HEH he : currentHes)
    {
        EH e = tetMesh.edge_handle(he);
        if (_forbiddenEs.find(e) != _forbiddenEs.end())
            nForbiddenEs++;
        for (VH v : tetMesh.edge_vertices(e))
            currentVs.insert(v);
    }

    for (VH v : currentVs)
        if (_forbiddenVs.find(v) != _forbiddenVs.end())
            nForbiddenVs++;

    return SUCCESS;
}

PathRouter::RetCode PathRouter::refineToAllowReroute(const FH& p, const VH& vFrom, const VH& vTo)
{
    auto& tetMesh = meshProps().mesh();

    // Split all edges that connect 2 arc vertices but are not arcs themselves
    auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);
    for (HFH hf : pHfs)
        for (EH e : tetMesh.halfface_edges(hf))
            if (meshProps().isInArc(e))
                _forbiddenEs.insert(e);

    for (EH e : _forbiddenEs)
        for (VH v : tetMesh.edge_vertices(e))
            if (v != vFrom && v != vTo)
                _forbiddenVs.insert(v);

    set<EH> esCut;
    for (HFH hf : pHfs)
        for (EH e : tetMesh.halfface_edges(hf))
            if (_forbiddenEs.count(e) == 0)
            {
                auto vs = tetMesh.edge_vertices(e);
                bool fromForbidden = _forbiddenVs.find(vs[0]) != _forbiddenVs.end();
                bool toForbidden = _forbiddenVs.find(vs[1]) != _forbiddenVs.end();
                if (fromForbidden && toForbidden)
                    esCut.insert(e);
            }

    for (EH e : esCut)
    {
        HEH he = tetMesh.halfedge_handle(e, 0);
        splitHalfEdge(he, *tetMesh.hec_iter(he), 0.5);
    }

    return SUCCESS;
}

bool PathRouter::incidentOnForbiddenElement(const HEH& he) const
{
    auto& tetMesh = meshProps().mesh();

    if (_forbiddenEs.find(tetMesh.edge_handle(he)) != _forbiddenEs.end())
        return true;
    for (VH v : tetMesh.halfedge_vertices(he))
        if (_forbiddenVs.find(v) != _forbiddenVs.end())
            return true;

    return false;
}

bool PathRouter::rerouteRecommended(const HEH& he,
                                    const list<HEH>& removedHes,
                                    const list<HEH>& addedHes,
                                    const set<VH>& removedVs,
                                    const set<VH>& addedVs) const
{
    auto& tetMesh = meshProps().mesh();
    HFH hf = getTransferredHf(he);

    (void)removedVs;

    if (!hf.is_valid())
        return false;

    for (HEH heAdded : addedHes)
        if (_forbiddenEs.count(tetMesh.edge_handle(heAdded)))
            return false;
    for (VH v : addedVs)
        if (_forbiddenVs.count(v))
            return false;

    return removedHes.size() > addedHes.size();
}

PathRouter::RetCode PathRouter::exchangedElements(const HFH& hf,
                                                  const set<HEH>& currentHes,
                                                  list<HEH>& removedHes,
                                                  list<HEH>& addedHes,
                                                  set<VH>& removedVs,
                                                  set<VH>& addedVs) const
{
    auto& tetMesh = meshProps().mesh();

    // Get the first removed halfedge along the path
    HEH heStart = *tetMesh.hfhe_iter(hf);
    while (currentHes.find(heStart) == currentHes.end()
           || currentHes.find(tetMesh.prev_halfedge_in_halfface(heStart, hf)) != currentHes.end())
    {
        heStart = tetMesh.next_halfedge_in_halfface(heStart, hf);
    }

    removedHes.emplace_back(heStart);
    bool previousRemoved = true;
    HEH he = heStart;
    for (he = tetMesh.next_halfedge_in_halfface(heStart, hf); he != heStart;
         he = tetMesh.next_halfedge_in_halfface(he, hf))
    {
        if (currentHes.find(he) == currentHes.end())
        {
            if (!previousRemoved)
                addedVs.insert(tetMesh.from_vertex_handle(he));
            addedHes.emplace_front(tetMesh.opposite_halfedge_handle(he));
            previousRemoved = false;
        }
        else
        {
            if (previousRemoved)
                removedVs.insert(tetMesh.from_vertex_handle(he));
            removedHes.emplace_back(he);
            previousRemoved = true;
        }
    }
    assert(previousRemoved == false);

    return SUCCESS;
}

HFH PathRouter::getTransferredHf(const HEH& pathHe) const
{
    auto& tetMesh = meshProps().mesh();
    const auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(_p);
    for (HFH hf : tetMesh.halfedge_halffaces(pathHe))
    {
        FH f = tetMesh.face_handle(hf);
        if (meshProps().get<MC_PATCH>(f) == _p && _inFirstHp == (pHfs.find(hf) != pHfs.end()))
            return hf;
    }
    return HFH(-1);
}

void PathRouter::reroute(list<HEH>& pathRerouted,
                         list<HEH>::iterator& pathIt,
                         set<HFH>& hfsTransferred,
                         set<VH>& currentVs,
                         set<HEH>& currentHes,
                         set<VH>& vsPassed,
                         list<HEH>& removedHes,
                         list<HEH>& addedHes,
                         set<VH>& removedVs,
                         set<VH>& addedVs)
{
    auto& tetMesh = meshProps().mesh();
    HEH heCurrent = *pathIt;
    auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(_p);
    hfsTransferred.insert(getTransferredHf(heCurrent));

    if (addedHes.size() == 2 && currentVs.find(tetMesh.to_vertex_handle(addedHes.front())) != currentVs.end())
    {
        assert(removedHes.front() == heCurrent);
        set<HEH> cycleHes;

        VH vPinch = tetMesh.to_vertex_handle(addedHes.front());
        addedVs.erase(vPinch);
        if (vsPassed.find(vPinch) != vsPassed.end())
        {
            // Pinch -> remove backward cycle
            cycleHes.insert(addedHes.front());
            addedHes.pop_front();
            while (tetMesh.from_vertex_handle(heCurrent) != vPinch)
            {
                pathIt--;
                heCurrent = *pathIt;
                removedHes.emplace_front(heCurrent);
                cycleHes.insert(heCurrent);
                removedVs.insert(tetMesh.to_vertex_handle(heCurrent));
            }
            assert(removedHes.front() == heCurrent);
        }
        else
        {
            // Pinch -> remove forward cycle
            cycleHes.insert(addedHes.back());
            addedHes.pop_back();
            auto itLocal = pathIt;
            while (tetMesh.to_vertex_handle(heCurrent) != vPinch)
            {
                itLocal++;
                heCurrent = *itLocal;
                removedHes.emplace_back(heCurrent);
                cycleHes.insert(heCurrent);
                removedVs.insert(tetMesh.from_vertex_handle(heCurrent));
            }
            heCurrent = *pathIt;
        }

        // Floodfill pocket hfs -> enclosedHfs
        list<HFH> hfQ;
        set<HFH> hfsToTransfer;

        HFH hfSeed = getTransferredHf(*cycleHes.begin());
        hfQ.emplace_back(hfSeed);
        hfsToTransfer.insert(hfSeed);

        while (!hfQ.empty())
        {
            HFH hf = hfQ.front();
            hfQ.pop_front();

            for (HEH he : tetMesh.halfface_halfedges(hf))
            {
                if (cycleHes.find(he) != cycleHes.end())
                    continue;
                HEH heOpp = tetMesh.opposite_halfedge_handle(he);

                for (HFH hfNext : tetMesh.halfedge_halffaces(heOpp))
                {
                    FH f = tetMesh.face_handle(hfNext);
                    if (meshProps().get<MC_PATCH>(f) == _p && _inFirstHp == (pHfs.find(hfNext) != pHfs.end())
                        && hfsToTransfer.find(hfNext) == hfsToTransfer.end())
                    {
                        hfQ.emplace_back(hfNext);
                        hfsToTransfer.insert(hfNext);
                        break;
                    }
                }
            }
        }

        hfsTransferred.insert(hfsToTransfer.begin(), hfsToTransfer.end());
    }
    else if (removedHes.front() != heCurrent)
    {
        pathIt--;
        heCurrent = *pathIt;
        assert(removedHes.front() == heCurrent);
    }
    for (HEH he : removedHes)
    {
        currentHes.erase(he);
        assert(heCurrent == he);
        pathIt = pathRerouted.erase(pathIt);
        heCurrent = *pathIt;
    }
    pathIt = pathRerouted.insert(pathIt, addedHes.begin(), addedHes.end());
    for (HEH he : addedHes)
        currentHes.insert(he);
    for (VH v : removedVs)
        currentVs.erase(v);
    for (VH v : addedVs)
        currentVs.insert(v);
}

} // namespace c4hex
