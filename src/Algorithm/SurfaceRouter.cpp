#include "C4Hex/Algorithm/SurfaceRouter.hpp"

#include <fstream>

#include <gurobi_c++.h>

namespace c4hex
{

SurfaceRouter::SurfaceRouter(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

// This function is so unwieldy because I want to keep GUROBI out of the header so I cant separate
// too much of the function into other functions.
SurfaceRouter::RetCode SurfaceRouter::calcMinimalSurfaceByLP(set<CH>& space,
                                                             const set<FH>& forbiddenFs,
                                                             const set<EH>& forbiddenEs,
                                                             const set<VH>& forbiddenVs,
                                                             const set<EH>& boundaryEs,
                                                             const set<HEH>& boundaryHes,
                                                             const set<VH>& boundaryVs,
                                                             set<HFH>& surfaceNew,
                                                             int depth)
{
    _forbiddenVs = forbiddenVs;
    _forbiddenEs = forbiddenEs;
    _forbiddenFs = forbiddenFs;

    auto& tetMesh = meshProps().mesh();
    if (depth >= 3)
        return NO_CONVERGENCE;

    assert(!space.empty());

    map<VH, int> v2valence;
    for (EH e : boundaryEs)
        for (VH v : tetMesh.edge_vertices(e))
            v2valence[v]++;
    set<VH> nonMfVs;
    for (auto kv : v2valence)
        if (kv.second != 2)
        {
            DLOG(WARNING) << "NON-2 valence at vertex " << kv.first << ", can not guarantee surface manifoldness";
            nonMfVs.insert(kv.first);
        }

    set<CH> confinedVolume = confineVolume(boundaryVs, space);

    for (bool confined : {true, false})
    {
        if (confined && (confinedVolume.empty() || confinedVolume == space))
            continue;

        surfaceNew.clear();
        {
            TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());
            {
                auto ret = refineVolumeToAllowReroute(confined ? confinedVolume : space, boundaryVs, boundaryEs);
                if (ret != SUCCESS)
                    return ret;
            }
            if (confined)
            {
                // Also update main space
                list<CH> children(space.begin(), space.end());
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
                space = {children.begin(), children.end()};
            }
        }

        try
        {
            GRBEnv env = GRBEnv(true);
            env.set(GRB_IntParam_Threads, 2);
            env.set(GRB_IntParam_LogToConsole, false);
            env.set(GRB_DoubleParam_TimeLimit, 5 * 60);
            env.start();

            GRBModel model = GRBModel(env);

            // Gather set of faces and edges
            map<HFH, GRBVar> hf2var;
            set<EH> edges;
            set<FH> faces;
            for (CH tet : confined ? confinedVolume : space)
                for (FH f : tetMesh.cell_faces(tet))
                {
                    if (_forbiddenFs.count(f) != 0 || containsSomeOf(tetMesh.face_edges(f), _forbiddenEs)
                        || containsSomeOf(tetMesh.face_vertices(f), _forbiddenVs))
                        continue;
                    faces.insert(f);
                    for (HFH hf : tetMesh.face_halffaces(f))
                        if (hf2var.count(hf) == 0)
                        {
                            for (EH e : tetMesh.halfface_edges(hf))
                                edges.insert(e);
                            hf2var[hf] = model.addVar(
                                0.0, 1.0, 0., GRB_CONTINUOUS, std::string("Halfface") + std::to_string(hf.idx()));
                        }
                }

            // Formulate objective: min z * Area(z)^T
            GRBQuadExpr objective = 0.;
            for (auto& kv : hf2var)
            {
                HFH hf = kv.first;
                CH tet = tetMesh.incident_cell(hf);
                if (!tet.is_valid())
                {
                    hf = tetMesh.opposite_halfface_handle(hf);
                    tet = tetMesh.incident_cell(hf);
                }
                auto vertices = meshProps().get_halfface_vertices(hf);

                vector<Vec3d> pos;
                for (VH v : vertices)
#ifdef MINIMIZE_XYZ
                    pos.emplace_back(tetMesh.vertex(v));
#else
                    pos.emplace_back(Vec3Q2d(meshProps().ref<CHART>(tet).at(v)));
#endif
                double area = ((pos[2] - pos[0]) % (pos[1] - pos[0])).length();
                area = std::max(area, 1e-6 * (1 + (double)rand() / RAND_MAX));

                objective += area * kv.second;
            }

            model.setObjective(objective, GRB_MINIMIZE);

            // Formulate constraints B * z = r
            // Formulate manifoldness constraints
            for (EH e : edges)
            {
                HEH he = tetMesh.halfedge_handle(e, 0);
                GRBLinExpr Bzi = 0;
                for (FH f : tetMesh.edge_faces(e))
                {
                    for (HFH hf : tetMesh.face_halffaces(f))
                    {
                        auto it = hf2var.find(hf);
                        if (it == hf2var.end())
                            continue;
                        bool coherent = false;
                        for (HEH heOther : tetMesh.halfface_halfedges(hf))
                            if (he == heOther)
                                coherent = true;
                        Bzi += coherent ? it->second : -it->second;
                    }
                }

                GRBLinExpr ri = 0;
                if (boundaryEs.count(e) != 0)
                {
                    ri = boundaryHes.count(he) == 0 ? -1 : 1;
                }

                std::string constraintName = "Edge" + std::to_string(e.idx());
                model.addConstr(Bzi, GRB_EQUAL, ri, constraintName);
            }

            bool nonManifoldSolution = true;

            set<EH> complexEs;
            set<VH> complexVs;

            int nIter = 0;

            while (nonManifoldSolution)
            {
                if (nIter++ >= 10)
                {
                    if (confined)
                        break; // next iteration -> non-confined
                    else
                        return NO_CONVERGENCE;
                }
                DLOG(INFO) << "Gurobi solving surface minimization...";
                model.optimize();

                int status = model.get(GRB_IntAttr_Status);
                if (status != 2 && status != 9)
                {
                    if (!surfaceNew.empty())
                    {
                        TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());

                        auto doNotSplitTheseEdges = _forbiddenEs;
                        doNotSplitTheseEdges.insert(boundaryEs.begin(), boundaryEs.end());
                        std::swap(doNotSplitTheseEdges, _forbiddenEs);
                        // Split
                        refineVolumeToAvoidComplex(confined ? confinedVolume : space, complexEs, complexVs, surfaceNew);
                        std::swap(doNotSplitTheseEdges, _forbiddenEs);

                        DLOG(INFO) << "SurfaceRouter: After refining complex elements, trying to resolve";
                        // resolve
                        auto ret = calcMinimalSurfaceByLP(confined ? confinedVolume : space,
                                                          _forbiddenFs,
                                                          _forbiddenEs,
                                                          _forbiddenVs,
                                                          boundaryEs,
                                                          boundaryHes,
                                                          boundaryVs,
                                                          surfaceNew,
                                                          depth + 1);
                        if (confined)
                        {
                            // Update main space
                            list<CH> children(space.begin(), space.end());
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
                            space = {children.begin(), children.end()};
                        }
                        if (ret != SUCCESS && confined)
                            break; // next iteration -> non-confined
                        else
                            return ret;
                    }
                    else
                    {
                        if (confined)
                            break; // next iteration -> non-confined
                        else
                        {
                            DLOG(ERROR) << "Bad status return by GUROBI solver";
                            return NO_CONVERGENCE;
                        }
                    }
                }
                auto obj = model.getObjective();
                DLOG(INFO) << "Solved with final surface area " << obj.getValue();

                bool nonIntegerSolution
                    = containsMatching(hf2var,
                                       [&](const pair<const HFH, GRBVar>& kv)
                                       {
                                           double val = kv.second.get(GRB_DoubleAttr_X);
                                           if (std::abs(val - std::round(val)) > 1e-3)
                                           {
                                               DLOG(WARNING)
                                                   << "Non integer solution " << val << " on hf " << kv.first
                                                   << ", can not solve by continuous LP, adding constraint to force 0";
                                               std::string constraintName = "NonIntF" + std::to_string(kv.first.idx());
                                               model.addConstr(kv.second, GRB_EQUAL, 0, constraintName);
                                               return true;
                                           }
                                           return false;
                                       });
                if (nonIntegerSolution)
                    continue;

                surfaceNew.clear();
                for (auto kv : hf2var)
                    if ((int)std::round(kv.second.get(GRB_DoubleAttr_X)) == 1)
                        surfaceNew.insert(kv.first);

                complexEs.clear();
                complexVs.clear();

                nonManifoldSolution = false;

                // Check and try to constrain complex edges
                {
                    map<EH, set<HFH>> esToFaces;
                    for (HFH hf : surfaceNew)
                        for (EH e : tetMesh.halfface_edges(hf))
                            esToFaces[e].insert(hf);

                    for (auto& kv : esToFaces)
                        if ((int)kv.second.size() != (2 - (int)boundaryEs.count(kv.first)))
                        {
                            DLOG(WARNING) << "Minimal surface has non-manifold edge " << kv.first
                                          << ", constraining to not include all current incident hfs";
                            complexEs.insert(kv.first);
                            nonManifoldSolution = true;

                            GRBLinExpr sum = 0;
                            for (HFH hf : kv.second)
                            {
                                auto var = hf2var.at(hf);
                                sum += var;
                            }

                            if (boundaryEs.count(kv.first) == 0)
                            {
                                DLOG(INFO) << "Adding constraint over " << kv.second.size()
                                           << " halffaces, that sum must be less/equal than 2";
                                std::string constraintName = "ComplexE" + std::to_string(kv.first.idx());
                                model.addConstr(sum, GRB_LESS_EQUAL, 2, constraintName);
                            }
                            else
                            {
                                DLOG(INFO) << "Adding constraint over " << kv.second.size()
                                           << " halffaces, that sum must be less/equal than 1";
                                std::string constraintName = "ComplexBoundaryE" + std::to_string(kv.first.idx());
                                model.addConstr(sum, GRB_LESS_EQUAL, 1, constraintName);
                            }
                        }
                    if (nonManifoldSolution)
                        continue;
                }

                // Check and try to constrain complex vertices
                // This is cumbersome as we have to account for (allowed) selfadjacency on the surface boundary
                {
                    map<VH, set<HFH>> vsToFaces;
                    for (HFH hf : surfaceNew)
                        for (VH v : meshProps().get_halfface_vertices(hf))
                            vsToFaces[v].insert(hf);

                    for (const auto& kv : vsToFaces)
                    {
                        if (nonMfVs.count(kv.first) != 0)
                            continue;
                        int nHfs = 0;
                        VH v = kv.first;
                        HFH hfSeed = findSomeOf(tetMesh.vertex_halffaces(v), surfaceNew);
                        set<HFH> hfVisited({{hfSeed}});
                        list<HFH> hfQ({{hfSeed}});
                        nHfs++;
                        while (!hfQ.empty())
                        {
                            HFH hf = hfQ.front();
                            hfQ.pop_front();
                            for (HEH he : tetMesh.halfface_halfedges(hf))
                            {
                                auto vs = tetMesh.halfedge_vertices(he);
                                if (vs[0] != v && vs[1] != v)
                                    continue;
                                HEH heOpp = tetMesh.opposite_halfedge_handle(he);
                                for (HFH hfNext : tetMesh.halfedge_halffaces(heOpp))
                                {
                                    if (surfaceNew.find(hfNext) != surfaceNew.end()
                                        && hfVisited.find(hfNext) == hfVisited.end())
                                    {
                                        hfVisited.insert(hfNext);
                                        hfQ.emplace_back(hfNext);
                                        nHfs++;
                                        break;
                                    }
                                }
                            }
                        }
                        if (nHfs != (int)kv.second.size())
                        {
                            DLOG(WARNING) << "Minimal surface has non-manifold vertex " << v
                                          << ", constraining to not include all current incident hfs";

                            complexVs.insert(v);
                            nonManifoldSolution = true;

                            GRBLinExpr sum = 0;
                            for (HFH hf : kv.second)
                            {
                                auto var = hf2var.at(hf);
                                sum += var;
                            }

                            DLOG(INFO) << "Adding constraint over " << kv.second.size()
                                       << " halffaces, that sum must be less/equal than "
                                       << (double)((int)kv.second.size() - 1);
                            std::string constraintName = "ComplexV" + std::to_string(v.idx());
                            model.addConstr(sum, GRB_LESS_EQUAL, (double)((int)kv.second.size() - 1), constraintName);
                        }
                    }
                    if (nonManifoldSolution)
                        continue;
                }

                // Check for genus > 0
                // This is cumbersome as we have to account for (allowed) selfadjacency on the surface boundary
                bool nonGenus0solution = false;
                {
                    map<VH, list<HFH>> vIncidence;
                    for (HFH hf : surfaceNew)
                        for (VH v : meshProps().get_halfface_vertices(hf))
                            vIncidence[v].push_back(hf);
                    map<VH, set<HFH>> v2firstSector;
                    for (auto kv : vIncidence)
                    {
                        auto& v = kv.first;
                        set<HFH> incidentHfs(kv.second.begin(), kv.second.end());
                        assert(incidentHfs.size() == kv.second.size());
                        set<HFH> localNeighbors;
                        for (int i = 0; i < 2; i++)
                        {
                            set<HFH> localDisk;
                            list<HFH> hfQ;
                            auto hfSeed = findNoneOf(incidentHfs, localNeighbors);
                            if (!hfSeed.is_valid())
                                continue;
                            localDisk.insert(hfSeed);
                            hfQ.push_back(hfSeed);
                            while (!hfQ.empty())
                            {
                                HFH hf = hfQ.front();
                                hfQ.pop_front();

                                for (HEH he : tetMesh.halfface_halfedges(hf))
                                {
                                    auto vs = tetMesh.halfedge_vertices(he);
                                    if ((vs[0] != v && vs[1] != v) || boundaryHes.count(he) != 0)
                                        continue;
                                    HEH heOpp = tetMesh.opposite_halfedge_handle(he);
                                    for (HFH hfNext : tetMesh.halfedge_halffaces(heOpp))
                                        if (localDisk.count(hfNext) == 0 && incidentHfs.count(hfNext) != 0)
                                        {
                                            localDisk.insert(hfNext);
                                            hfQ.push_back(hfNext);
                                        }
                                }
                            }
                            v2firstSector[kv.first] = localDisk;
                            localNeighbors.insert(localDisk.begin(), localDisk.end());
                        }
                        assert(localNeighbors == incidentHfs);
                    }

                    {
                        set<HEH> pHes;
                        set<VH> pVs;
                        set<VH> pVsAlt;
                        for (HFH hf2 : surfaceNew)
                        {
                            for (HEH he : tetMesh.halfface_halfedges(hf2))
                            {
                                if (boundaryHes.count(he) != 0)
                                    pHes.insert(he);
                                else if ((he.idx() % 2) == 0)
                                    pHes.insert(he);
                            }
                            for (VH v : meshProps().get_halfface_vertices(hf2))
                            {
                                if (v2firstSector.at(v).count(hf2) != 0)
                                    pVs.insert(v);
                                else
                                    pVsAlt.insert(v);
                            }
                        }
                        if ((int)surfaceNew.size() - (int)pHes.size() + (int)pVs.size() + (int)pVsAlt.size() != 1)
                        {
                            DLOG(WARNING) << "Genus of minimal surface > 0"
                                          << ", can not solve by continuous LP";
                            nonGenus0solution = true;
                            nonManifoldSolution = true;
                        }
                    }
                }
                if (nonGenus0solution)
                {
                    if (confined)
                        break; // next iteration -> non-confined
                    else
                        return NO_CONVERGENCE;
                }
            }
            if (!nonManifoldSolution)
                return SUCCESS;
        }
        catch (GRBException e)
        {
            DLOG(ERROR) << "Gurobi exception, errcode: " << e.getErrorCode();
            DLOG(ERROR) << "Gurobi error message: " << e.getMessage();
            if (confined)
                continue; // next iteration -> non-confined
            return NO_CONVERGENCE;
        }
    }
    return SUCCESS;
}

SurfaceRouter::RetCode
SurfaceRouter::rerouteSurfaceThroughBlock(const CH& b, set<HFH>& surface, set<CH>& transferredTets)
{
    DLOG(INFO) << "Rerouting surface through block " << b;

    // Try solving with minimal surface
    {
        auto surfaceNew = surface;
        if (minimalSurfaceThroughBlock(b, surfaceNew, transferredTets) == SUCCESS)
        {
            surface = surfaceNew;
            return SUCCESS;
        }
        else
        {
            LOG(INFO) << "minimal surface through block didnt work, now incrementally rerouting";
            transferredTets.clear();
        }
    }

    // Else fall back to incremental shift
    auto ret = shiftSurfaceThroughBlock(b, surface, transferredTets);
    return ret;
}

SurfaceRouter::RetCode
SurfaceRouter::minimalSurfaceThroughBlock(const CH& b, set<HFH>& surface, set<CH>& transferredTets)
{
    DLOG(INFO) << "Determining minimal surface through block " << b;
    auto& tetMesh = meshProps().mesh();

    CH tetAny = tetMesh.incident_cell(*surface.begin());
    bool frontHalfface = tetAny.is_valid() && meshProps().get<MC_BLOCK>(tetAny) == b;

    // Gather contour:
    set<HEH> boundaryHes;
    set<EH> boundaryEs;
    set<VH> boundaryVs;
    getSurfaceBoundary(surface, boundaryVs, boundaryEs, boundaryHes);

    setForbiddenElements(b);
    for (EH e : boundaryEs)
        _forbiddenEs.erase(e);
    for (VH v : boundaryVs)
        _forbiddenVs.erase(v);

    auto space = mcMeshProps().get<BLOCK_MESH_TETS>(b);
    set<HFH> surfaceNew;
    auto ret = calcMinimalSurfaceByLP(
        space, _forbiddenFs, _forbiddenEs, _forbiddenVs, boundaryEs, boundaryHes, boundaryVs, surfaceNew);
    if (ret != SUCCESS)
        return ret;
    assert(surfaceNew.size() != 0);

    // Floodfill space between old and new surface to find transferred tets
    vector<bool> tetVisited(tetMesh.n_cells(), false);
    for (CH tetSeed : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        if (tetVisited[tetSeed.idx()])
            continue;
        list<CH> tetQ({tetSeed});
        tetVisited[tetSeed.idx()] = true;

        bool wrongSide = false;
        bool onlyNewSurfaceTouched = true;
        vector<CH> tetsFlooded;
        while (!tetQ.empty())
        {
            CH tet = tetQ.back();
            tetQ.pop_back();
            assert(!tetMesh.is_deleted(tet));

            assert(meshProps().get<MC_BLOCK>(tet) == b);
            tetsFlooded.push_back(tet);

            for (HFH hf : tetMesh.cell_halffaces(tet))
            {
                HFH hfOpp = tetMesh.opposite_halfface_handle(hf);

                bool touchedNewSurfaceFront = surfaceNew.count(hf) > 0;
                bool touchedNewSurfaceBack = surfaceNew.count(hfOpp) > 0;
                // Do not spread beyond boundary
                if (touchedNewSurfaceFront || touchedNewSurfaceBack)
                {
                    if (touchedNewSurfaceFront == frontHalfface)
                        wrongSide = true;
                    continue;
                }

                // Register if correct flooded half
                if (surface.count(hf) > 0 || surface.count(hfOpp) > 0)
                {
                    onlyNewSurfaceTouched = false;
                    continue;
                }

                CH tetNext = tetMesh.incident_cell(hfOpp);
                if (!tetNext.is_valid() || meshProps().get<MC_BLOCK>(tetNext) != b)
                    onlyNewSurfaceTouched = false;
                if (tetNext.is_valid() && !tetVisited[tetNext.idx()] && meshProps().get<MC_BLOCK>(tetNext) == b)
                {
                    tetVisited[tetNext.idx()] = true;
                    tetQ.push_back(tetNext);
                }
            }
        }
        if (!wrongSide)
        {
            if (onlyNewSurfaceTouched)
            {
                LOG(ERROR) << "Touched only new surface with " << surfaceNew.size() << " hfs, flooded "
                           << tetsFlooded.size() << " tets";
                throw std::logic_error("Pocket created, should be impossible");
            }
            transferredTets.insert(tetsFlooded.begin(), tetsFlooded.end());
        }
    }

    surface = surfaceNew;

    return SUCCESS;
}

SurfaceRouter::RetCode SurfaceRouter::shiftSurfaceThroughBlock(const CH& b, set<HFH>& surface, set<CH>& transferredTets)
{
    auto& tetMesh = meshProps().mesh();

    setForbiddenElements(b);
    {
        auto ret = refineToAllowReroute(b, surface);
        if (ret != SUCCESS)
            return ret;
    }

    set<HEH> boundaryHes;
    set<EH> boundaryEs;
    set<VH> boundaryVs;
    getSurfaceBoundary(surface, boundaryVs, boundaryEs, boundaryHes);

    for (EH e : boundaryEs)
        _forbiddenEs.erase(e);
    for (VH v : boundaryVs)
        _forbiddenVs.erase(v);

    set<HFH> currentHfs = surface;
    set<EH> currentEs;
    set<VH> currentVs;
    int nForbiddenFs = 0;
    int nForbiddenEs = 0;
    int nForbiddenVs = 0;

    {
        auto ret = gatherCurrentElements(currentHfs, currentEs, currentVs, nForbiddenFs, nForbiddenEs, nForbiddenVs);
        if (ret != SUCCESS)
            return ret;
    }

    // Determine shift/flip direction
    CH tetAny = tetMesh.incident_cell(*currentHfs.begin());
    bool frontHalfface = tetAny.is_valid() && meshProps().get<MC_BLOCK>(tetAny) == b;

    set<HFH> surfaceHfsToFlip;
    list<HFH> hfsList;

    for (HFH hf : surface)
    {
        if (!frontHalfface)
            hf = tetMesh.opposite_halfface_handle(hf);
        if (incidentOnForbiddenElement(hf))
        {
            surfaceHfsToFlip.insert(hf);
            hfsList.push_back(hf);
        }
    }

    int nRefined = 0;
    int iter = 0;
    int nForbiddenLastV = nForbiddenVs;
    int nForbiddenLastE = nForbiddenEs;
    int nForbiddenLastF = nForbiddenFs;
    int nInArow = 0;
    while (!surfaceHfsToFlip.empty())
    {
        iter++;
        if ((iter % 10000) == 0)
        {
            DLOG(INFO) << "Still have to lift surface off of: " << nForbiddenFs << " fs , " << nForbiddenEs << " es, "
                       << nForbiddenVs << " vs";
            if (nForbiddenVs == nForbiddenLastV && nForbiddenEs == nForbiddenLastE && nForbiddenFs == nForbiddenLastF)
            {
                nInArow++;
                if (nInArow > 10)
                    throw std::logic_error("Infinite patch pulling");
            }
            else
                nInArow = 0;
            nForbiddenLastV = nForbiddenVs;
            nForbiddenLastE = nForbiddenEs;
            nForbiddenLastF = nForbiddenFs;
        }

        int listLength = hfsList.size();
        // This should not happen, but if it does:
        // Try to salvage by refining a lot
        if (nInArow == 3)
        {
            set<CH> tets;
            for (HFH hf : surfaceHfsToFlip)
                tets.insert(tetMesh.incident_cell(hf));
            for (CH tet : tets)
                splitTet(tet, {Q(1, 4), Q(1, 4), Q(1, 4), Q(1, 4)});
            for (HFH hf : currentHfs)
            {
                if (!frontHalfface)
                    hf = tetMesh.opposite_halfface_handle(hf);
                if (incidentOnForbiddenElement(hf))
                {
                    if (surfaceHfsToFlip.count(hf) == 0)
                    {
                        surfaceHfsToFlip.insert(hf);
                        hfsList.push_back(hf);
                    }
                }
            }
        }
        else if (nInArow == 6)
        {
            set<CH> tets;
            for (HFH hf : hfsList)
                tets.insert(tetMesh.incident_cell(hf));
            for (CH tet : tets)
                for (HFH hf : tetMesh.cell_halffaces(tet))
                    if (currentHfs.count(hf) == 0 && currentHfs.count(tetMesh.opposite_halfface_handle(hf)) == 0
                        && _forbiddenFs.count(tetMesh.face_handle(hf)) == 0)
                        splitFace(tetMesh.face_handle(hf), {Q(1, 4), Q(1, 4), Q(1, 4)});
            for (HFH hf : currentHfs)
            {
                if (!frontHalfface)
                    hf = tetMesh.opposite_halfface_handle(hf);
                if (incidentOnForbiddenElement(hf))
                {
                    if (surfaceHfsToFlip.count(hf) == 0)
                    {
                        surfaceHfsToFlip.insert(hf);
                        hfsList.push_back(hf);
                    }
                }
            }
        }

        // First: only transfer tets that do not create nonmanifold surface
        CH transferredTet;
        for (int i = 0; i < listLength; i++)
        {
            HFH hf = hfsList.front();
            hfsList.pop_front();
            if (surfaceHfsToFlip.count(hf) != 0)
            {
                assert(transferredTets.count(tetMesh.incident_cell(hf)) == 0);
                if (transferCreatesNonmanifold(hf, currentHfs, currentEs, currentVs) == NONMF_NONE)
                {
                    CH tet = tetMesh.incident_cell(hf);
                    assert(tet.is_valid() && !tetMesh.is_deleted(tet) && meshProps().get<MC_BLOCK>(tet) == b);
                    for (HFH hfTet : tetMesh.cell_halffaces(tet))
                        surfaceHfsToFlip.erase(hfTet);
                    transferredTet = tet;
                    break;
                }
                else
                    hfsList.push_back(hf);
            }
        }

        // ...else transfer tet that would create complex vertex but refine before
        if (!transferredTet.is_valid())
        {
            for (int i = 0; i < listLength; i++)
            {
                HFH hf = hfsList.front();
                assert(transferredTets.count(tetMesh.incident_cell(hf)) == 0);
                hfsList.pop_front();
                auto nonMfType = transferCreatesNonmanifold(hf, currentHfs, currentEs, currentVs);
                assert(nonMfType != NONMF_NONE);
                if (nonMfType == NONMF_VERTEX)
                {
                    DLOG(INFO) << "Refining to avoid non-manifold vertex";
                    auto ret = refineToAvoidNonmanifoldVertex(hf);
                    if (ret != SUCCESS)
                        return ret;
                    assert(transferCreatesNonmanifold(hf, currentHfs, currentEs, currentVs) == NONMF_NONE);
                    assert(!tetMesh.is_deleted(hf));
                    CH tet = tetMesh.incident_cell(hf);
                    assert(tet.is_valid() && !tetMesh.is_deleted(tet) && meshProps().get<MC_BLOCK>(tet) == b);
                    for (HFH hfTet : tetMesh.cell_halffaces(tet))
                        surfaceHfsToFlip.erase(hfTet);
                    transferredTet = tet;
                    nRefined++;
                    break;
                }
                else
                    hfsList.push_back(hf);
            }
        }

        // ...else try to find an entire pocket with triangular opening, that can not otherwise be resolved,
        // or transfer a tet that would create a complex edge but refine before
        if (!transferredTet.is_valid())
        {
            // There should either be a pocket which we can find and eliminate
            if (findAndTransferPocket(boundaryEs,
                                      frontHalfface,
                                      transferredTets,
                                      surfaceHfsToFlip,
                                      hfsList,
                                      currentHfs,
                                      currentEs,
                                      currentVs,
                                      nForbiddenFs,
                                      nForbiddenEs,
                                      nForbiddenVs))
            {
                continue;
            }
            // ... or the nonManifold transfer is avoidable through refinement
            for (int i = 0; i < listLength; i++)
            {
                HFH hf = hfsList.front();
                hfsList.pop_front();
                assert(transferredTets.count(tetMesh.incident_cell(hf)) == 0);
                auto nonMfType = transferCreatesNonmanifold(hf, currentHfs, currentEs, currentVs);
                assert(nonMfType != NONMF_NONE);
                if (nonMfType == NONMF_EDGE)
                {
                    DLOG(INFO) << "Refining to avoid non-manifold edge";
                    auto ret = refineToAvoidNonmanifoldEdge(tetMesh.incident_cell(hf), currentHfs);
                    if (ret != SUCCESS)
                        return ret;
                    assert(transferCreatesNonmanifold(hf, currentHfs, currentEs, currentVs) == NONMF_NONE);
                    assert(!tetMesh.is_deleted(hf));
                    CH tet = tetMesh.incident_cell(hf);
                    assert(tet.is_valid() && !tetMesh.is_deleted(tet) && meshProps().get<MC_BLOCK>(tet) == b);
                    for (HFH hfTet : tetMesh.cell_halffaces(tet))
                        surfaceHfsToFlip.erase(hfTet);
                    transferredTet = tet;
                    break;
                }
                else
                    hfsList.push_back(hf);
            }
        }
        transferredTets.insert(transferredTet);

        // Now that the tet to transfer is chosen: execute the traversal/transfer

        set<HFH> removedHfs;
        set<HFH> addedHfs;
        set<EH> removedEs;
        set<EH> addedEs;
        set<VH> removedVs;
        set<VH> addedVs;
        {
            auto ret = exchangedElements(transferredTet,
                                         currentHfs,
                                         frontHalfface,
                                         removedHfs,
                                         addedHfs,
                                         removedEs,
                                         addedEs,
                                         removedVs,
                                         addedVs);
            if (ret != SUCCESS)
                return ret;
        }
        assert(removedHfs.size() <= 3);
        assert(addedHfs.size() <= 3);
        assert(removedEs.size() <= 3);
        assert(addedEs.size() <= 3);
        assert(removedVs.size() <= 1);
        assert(addedVs.size() <= 1);

        for (VH v : removedVs)
        {
            auto it = currentVs.find(v);
            assert(it != currentVs.end());
            currentVs.erase(it);
            if (_forbiddenVs.find(v) != _forbiddenVs.end())
                nForbiddenVs--;
        }
        for (VH v : addedVs)
        {
            assert(currentVs.find(v) == currentVs.end());
            currentVs.insert(v);
            assert(currentVs.find(v) != currentVs.end() || _forbiddenVs.find(v) == _forbiddenVs.end());
        }
        for (EH e : removedEs)
        {
            auto it = currentEs.find(e);
            assert(it != currentEs.end());
            currentEs.erase(it);
            if (_forbiddenEs.find(e) != _forbiddenEs.end())
                nForbiddenEs--;
        }
        for (EH e : addedEs)
        {
            assert(currentEs.find(e) == currentEs.end());
            currentEs.insert(e);
            assert(currentEs.find(e) != currentEs.end() || _forbiddenEs.find(e) == _forbiddenEs.end());
        }
        for (HFH hf : removedHfs)
        {
            auto it = currentHfs.find(hf);
            assert(it != currentHfs.end());
            currentHfs.erase(it);
            if (_forbiddenFs.find(tetMesh.face_handle(hf)) != _forbiddenFs.end())
                nForbiddenFs--;
        }
        for (HFH hf : addedHfs)
        {
            assert(currentHfs.find(hf) == currentHfs.end());
            currentHfs.insert(hf);
            assert(_forbiddenFs.find(tetMesh.face_handle(hf)) == _forbiddenFs.end());
        }

        for (HFH hf : addedHfs)
        {
            if (!frontHalfface)
                hf = tetMesh.opposite_halfface_handle(hf);
            if (incidentOnForbiddenElement(hf))
            {
                surfaceHfsToFlip.insert(hf);
                hfsList.push_back(hf);
            }
        }
    }

    surface = currentHfs;

    return SUCCESS;
}

SurfaceRouter::RetCode SurfaceRouter::gatherCurrentElements(const set<HFH>& currentHfs,
                                                            set<EH>& currentEs,
                                                            set<VH>& currentVs,
                                                            int& nForbiddenFs,
                                                            int& nForbiddenEs,
                                                            int& nForbiddenVs) const
{
    auto& tetMesh = meshProps().mesh();

    nForbiddenFs = 0;
    nForbiddenEs = 0;
    nForbiddenVs = 0;

    for (HFH hf : currentHfs)
    {
        FH f = tetMesh.face_handle(hf);
        if (_forbiddenFs.find(f) != _forbiddenFs.end())
            nForbiddenFs++;
        for (EH e : tetMesh.face_edges(f))
            currentEs.insert(e);
        for (VH v : tetMesh.face_vertices(f))
            currentVs.insert(v);
    }

    for (EH e : currentEs)
        if (_forbiddenEs.find(e) != _forbiddenEs.end())
            nForbiddenEs++;
    for (VH v : currentVs)
        if (_forbiddenVs.find(v) != _forbiddenVs.end())
            nForbiddenVs++;

    return SUCCESS;
}

bool SurfaceRouter::incidentOnForbiddenElement(const HFH& hf) const
{
    auto& tetMesh = meshProps().mesh();

    if (_forbiddenFs.find(tetMesh.face_handle(hf)) != _forbiddenFs.end())
        return true;
    for (EH e : tetMesh.halfface_edges(hf))
        if (_forbiddenEs.find(e) != _forbiddenEs.end())
            return true;
    for (VH v : meshProps().get_halfface_vertices(hf))
        if (_forbiddenVs.find(v) != _forbiddenVs.end())
            return true;

    return false;
}

SurfaceRouter::NonMF SurfaceRouter::transferCreatesNonmanifold(const HFH& hfFlip,
                                                               const set<HFH>& currentHfs,
                                                               const set<EH>& currentEs,
                                                               const set<VH>& currentVs) const
{
    auto& tetMesh = meshProps().mesh();

    {
        HFH hf = hfFlip;
        HEH he = *tetMesh.hfhe_iter(hf);
        hf = tetMesh.adjacent_halfface_in_cell(hfFlip, he);
        he = tetMesh.next_halfedge_in_halfface(tetMesh.opposite_halfedge_handle(he), hf);

        for (int i = 0; i < 3; i++)
        {
            HFH hfOpp = tetMesh.opposite_halfface_handle(hf);

            bool eIsCurrent = currentEs.find(tetMesh.edge_handle(he)) != currentEs.end();
            if (eIsCurrent)
            {
                HFH hfNext = tetMesh.adjacent_halfface_in_cell(hf, he);
                HFH hfNextOpp = tetMesh.opposite_halfface_handle(hfNext);
                bool fIsRemovedNotAdded
                    = currentHfs.find(hf) != currentHfs.end() || currentHfs.find(hfOpp) != currentHfs.end();
                bool fNextIsRemovedNotAdded
                    = currentHfs.find(hfNext) != currentHfs.end() || currentHfs.find(hfNextOpp) != currentHfs.end();
                if (!fIsRemovedNotAdded && !fNextIsRemovedNotAdded)
                    return NONMF_EDGE;
            }

            hf = tetMesh.adjacent_halfface_in_cell(hf, he);
            he = tetMesh.prev_halfedge_in_halfface(tetMesh.opposite_halfedge_handle(he), hf);
        }
    }

    VH v = tetMesh.halfface_opposite_vertex(hfFlip);
    bool vIsCurrent = currentVs.find(v) != currentVs.end();
    if (!vIsCurrent)
        return NONMF_NONE;

    vector<bool> removedNotAdded;
    for (HEH he : tetMesh.halfface_halfedges(hfFlip))
    {
        HFH hfAdj = tetMesh.adjacent_halfface_in_cell(hfFlip, he);
        HFH hfAdjOpp = tetMesh.opposite_halfface_handle(hfAdj);
        removedNotAdded.emplace_back(currentHfs.find(hfAdj) != currentHfs.end()
                                     || currentHfs.find(hfAdjOpp) != currentHfs.end());
    }
    assert(removedNotAdded.size() == 3);

    if (!removedNotAdded[0] && !removedNotAdded[1] && !removedNotAdded[2])
        return NONMF_VERTEX;

    return NONMF_NONE;
}

bool SurfaceRouter::findAndTransferPocket(const set<EH>& boundaryEs,
                                          bool frontHalfface,
                                          set<CH>& transferredTets,
                                          set<HFH>& surfaceHfsToFlip,
                                          list<HFH>& hfsList,
                                          set<HFH>& currentHfs,
                                          set<EH>& currentEs,
                                          set<VH>& currentVs,
                                          int& nForbiddenHfs,
                                          int& nForbiddenEs,
                                          int& nForbiddenVs)
{
    auto& tetMesh = meshProps().mesh();
    // Find pocket entrance adjacent to a surfaceHf that must be flipped
    HFH pocketHf;
    for (HFH hf : surfaceHfsToFlip)
    {
        for (HEH he : tetMesh.halfface_halfedges(hf))
        {
            HFH hfAdj = tetMesh.adjacent_halfface_in_cell(hf, he);
            assert(currentHfs.find(frontHalfface ? tetMesh.opposite_halfface_handle(hfAdj) : hfAdj)
                   == currentHfs.end());
            if (currentHfs.count(frontHalfface ? hfAdj : tetMesh.opposite_halfface_handle(hfAdj)) == 0
                &&!containsMatching(tetMesh.halfface_edges(hfAdj),
                                    [&](const EH& e) { return currentEs.count(e) == 0; }))
            {
                pocketHf = hfAdj;
                break;
            }
        }
        if (pocketHf.is_valid())
            break;
    }
    if (!pocketHf.is_valid())
    {
        DLOG(INFO) << "No resolvable pocket found";
        return false;
    }

    bool oppHf = false;

    auto itPairEs = tetMesh.halfface_edges(pocketHf);
    vector<EH> lockedEs(itPairEs.first, itPairEs.second);
    list<HFH> hfQ;
    set<HFH> hfVisited;
    for (HEH he : tetMesh.halfface_halfedges(pocketHf))
    {
        for (HFH hfNext : tetMesh.halfedge_halffaces(tetMesh.opposite_halfedge_handle(he)))
        {
            assert(hfNext != pocketHf);
            if (currentHfs.find(frontHalfface ? hfNext : tetMesh.opposite_halfface_handle(hfNext)) != currentHfs.end())
            {
                hfVisited.insert(hfNext);
                hfQ.emplace_back(hfNext);
            }
        }
    }

    while (!hfQ.empty())
    {
        HFH hf = hfQ.front();
        hfQ.pop_front();

        for (HEH he : tetMesh.halfface_halfedges(hf))
        {
            EH e = tetMesh.edge_handle(he);
            if (e == lockedEs[0] || e == lockedEs[1] || e == lockedEs[2])
                continue;
            if (boundaryEs.find(e) != boundaryEs.end())
            {
                DLOG(INFO) << "Encountered boundary edge in floodfill search: " << e;
                oppHf = true;
                break;
            }
            for (HFH hfNext : tetMesh.halfedge_halffaces(tetMesh.opposite_halfedge_handle(he)))
            {
                if (currentHfs.find(frontHalfface ? hfNext : tetMesh.opposite_halfface_handle(hfNext))
                        != currentHfs.end()
                    && hfVisited.find(hfNext) == hfVisited.end())
                {
                    hfQ.emplace_front(hfNext);
                    hfVisited.insert(hfNext);
                }
            }
        }
        if (oppHf)
            break;
    }
    assert(oppHf || hfQ.empty());

    if (oppHf)
        pocketHf = tetMesh.opposite_halfface_handle(pocketHf);

    assert(tetMesh.incident_cell(pocketHf).is_valid());
    list<CH> tetQ({tetMesh.incident_cell(pocketHf)});
    set<CH> pocketTets({tetQ.front()});
    vector<HFH> removedHfs;
    while (!tetQ.empty())
    {
        CH tetPocket = tetQ.front();
        tetQ.pop_front();

        for (HFH hf : tetMesh.cell_halffaces(tetPocket))
        {
            if (hf != pocketHf)
            {
                auto it = currentHfs.find(frontHalfface ? hf : tetMesh.opposite_halfface_handle(hf));
                assert(currentHfs.find(frontHalfface ? tetMesh.opposite_halfface_handle(hf) : hf) == currentHfs.end());
                if (it == currentHfs.end())
                {
                    CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                    assert(tetNext.is_valid());
                    if (pocketTets.find(tetNext) == pocketTets.end())
                    {
                        pocketTets.insert(tetNext);
                        tetQ.emplace_back(tetNext);
                    }
                }
                else
                {
                    removedHfs.emplace_back(frontHalfface ? hf : tetMesh.opposite_halfface_handle(hf));
                }
            }
        }
    }
    assert(pocketTets.size() >= 1);
    assert(removedHfs.size() >= 3);
    for (CH tet : pocketTets)
        transferredTets.insert(tet);

    array<VH, 3> lockedVs(meshProps().get_halfface_vertices(pocketHf));
    set<VH> removedVs;
    set<EH> removedEs;

    for (HFH hf : removedHfs)
    {
        if (_forbiddenFs.find(tetMesh.face_handle(hf)) != _forbiddenFs.end())
            nForbiddenHfs--;
        auto it = currentHfs.find(hf);
        assert(it != currentHfs.end());
        currentHfs.erase(it);
        surfaceHfsToFlip.erase(frontHalfface ? hf : tetMesh.opposite_halfface_handle(hf));

        for (EH e : tetMesh.halfface_edges(hf))
            if (e != lockedEs[0] && e != lockedEs[1] && e != lockedEs[2])
                removedEs.insert(e);
    }

    for (EH e : removedEs)
    {
        if (_forbiddenEs.find(e) != _forbiddenEs.end())
            nForbiddenEs--;
        auto it = currentEs.find(e);
        assert(it != currentEs.end());
        currentEs.erase(it);
        for (VH v : tetMesh.edge_vertices(e))
            if (v != lockedVs[0] && v != lockedVs[1] && v != lockedVs[2])
                removedVs.insert(v);
    }
    for (VH v : removedVs)
    {
        if (_forbiddenVs.find(v) != _forbiddenVs.end())
            nForbiddenVs--;
        auto it = currentVs.find(v);
        assert(it != currentVs.end());
        currentVs.erase(it);
    }

    assert(nForbiddenHfs >= 0);
    assert(nForbiddenEs >= 0);
    assert(nForbiddenVs >= 0);

    currentHfs.insert(frontHalfface ? tetMesh.opposite_halfface_handle(pocketHf) : pocketHf);

    if (incidentOnForbiddenElement(tetMesh.opposite_halfface_handle(pocketHf)))
    {
        surfaceHfsToFlip.insert(tetMesh.opposite_halfface_handle(pocketHf));
        hfsList.push_back(tetMesh.opposite_halfface_handle(pocketHf));
    }

    return true;
}

SurfaceRouter::RetCode SurfaceRouter::exchangedElements(const CH& tet,
                                                        const set<HFH>& currentHfs,
                                                        bool frontHalfface,
                                                        set<HFH>& removedHfs,
                                                        set<HFH>& addedHfs,
                                                        set<EH>& removedEs,
                                                        set<EH>& addedEs,
                                                        set<VH>& removedVs,
                                                        set<VH>& addedVs) const
{
    DLOG(INFO) << "Checking the transferred elements of tet " << tet;
    auto& tetMesh = meshProps().mesh();

    for (HFH hf : tetMesh.cell_halffaces(tet))
    {
        HFH hfOpp = tetMesh.opposite_halfface_handle(hf);
        bool current = currentHfs.find(frontHalfface ? hf : hfOpp) != currentHfs.end();
        if (current)
            removedHfs.insert(frontHalfface ? hf : hfOpp);
        else
            addedHfs.insert(frontHalfface ? hfOpp : hf);
    }

    {
        HFH hf = *tetMesh.chf_iter(tet);
        HEH he = *tetMesh.hfhe_iter(hf);

        for (int i = 0; i < 6; i++)
        {
            EH e = tetMesh.edge_handle(he);
            HFH hfOpp = tetMesh.opposite_halfface_handle(hf);

            HFH hfNext = tetMesh.adjacent_halfface_in_cell(hf, he);
            HFH hfNextOpp = tetMesh.opposite_halfface_handle(hfNext);
            bool fIsRemovedNotAdded = removedHfs.find(frontHalfface ? hf : hfOpp) != removedHfs.end();
            bool fNextIsRemovedNotAdded = removedHfs.find(frontHalfface ? hfNext : hfNextOpp) != removedHfs.end();
            if (fIsRemovedNotAdded && fNextIsRemovedNotAdded)
                removedEs.insert(e);
            else if (!fIsRemovedNotAdded && !fNextIsRemovedNotAdded)
                addedEs.insert(e);

            if (i < 2)
                he = tetMesh.next_halfedge_in_halfface(he, hf);
            else if (i == 2)
            {
                hf = tetMesh.adjacent_halfface_in_cell(hf, he);
                he = tetMesh.next_halfedge_in_halfface(tetMesh.opposite_halfedge_handle(he), hf);
            }
            else
            {
                hf = tetMesh.adjacent_halfface_in_cell(hf, he);
                he = tetMesh.prev_halfedge_in_halfface(tetMesh.opposite_halfedge_handle(he), hf);
            }
        }
    }

    for (VH v : tetMesh.tet_vertices(tet))
    {
        vector<bool> removedNotAdded;
        for (HFH hf : tetMesh.vertex_halffaces(v))
            if (tetMesh.incident_cell(hf) == tet)
            {
                HFH hfOpp = tetMesh.opposite_halfface_handle(hf);
                removedNotAdded.emplace_back(removedHfs.find(frontHalfface ? hf : hfOpp) != removedHfs.end());
            }

        if (removedNotAdded[0] && removedNotAdded[1] && removedNotAdded[2])
            removedVs.insert(v);
        else if (!removedNotAdded[0] && !removedNotAdded[1] && !removedNotAdded[2])
            addedVs.insert(v);
    }

    return SUCCESS;
}

SurfaceRouter::RetCode SurfaceRouter::refineToAllowReroute(const CH& b, set<HFH>& surface)
{
    auto& tetMesh = meshProps().mesh();

    set<EH> surfaceEs;
    set<VH> surfaceVs;
    for (HFH hf : surface)
    {
        for (EH e : tetMesh.halfface_edges(hf))
            surfaceEs.insert(e);
        for (VH v : meshProps().get_halfface_vertices(hf))
            surfaceVs.insert(v);
    }

    set<EH> esCut;
    set<FH> fsSplit;
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        for (EH e : tetMesh.cell_edges(tet))
            if (esCut.find(e) == esCut.end() && _forbiddenEs.find(e) == _forbiddenEs.end()
                && surfaceEs.find(e) == surfaceEs.end())
            {
                auto vs = tetMesh.edge_vertices(e);
                bool fromForbidden = _forbiddenVs.find(vs[0]) != _forbiddenVs.end();
                bool toForbidden = _forbiddenVs.find(vs[1]) != _forbiddenVs.end();
                bool oneCurrent = surfaceVs.find(vs[0]) != surfaceVs.end() || surfaceVs.find(vs[1]) != surfaceVs.end();
                if (fromForbidden && toForbidden && oneCurrent)
                    esCut.insert(e);
            }
        for (HFH hf : tetMesh.cell_halffaces(tet))
        {
            FH f = tetMesh.face_handle(hf);
            if (fsSplit.find(f) == fsSplit.end() && _forbiddenFs.find(f) == _forbiddenFs.end()
                && surface.find(hf) == surface.end()
                && surface.find(tetMesh.opposite_halfface_handle(hf)) == surface.end())
            {
                vector<bool> forbidden;
                vector<bool> current;
                for (EH e : tetMesh.face_edges(f))
                {
                    forbidden.emplace_back(_forbiddenEs.find(e) != _forbiddenEs.end());
                    current.emplace_back(surfaceEs.find(e) != surfaceEs.end());
                }
                bool someForbidden = (forbidden[0] || forbidden[1] || forbidden[2]);
                bool someCurrent = (current[0] || current[1] || current[2]);
                bool allCurrentOrForbidden
                    = (forbidden[0] || current[0]) && (forbidden[1] || current[1]) && (forbidden[2] || current[2]);
                if (someForbidden && someCurrent && allCurrentOrForbidden)
                    fsSplit.insert(f);
            }
        }
    }

    for (EH e : esCut)
    {
        HEH he = tetMesh.halfedge_handle(e, 0);
        splitHalfEdge(he, *tetMesh.hec_iter(he), 0.5);
    }

    Vec3Q barCoords(1, 1, 1);
    barCoords /= 3;
    for (FH f : fsSplit)
        splitFace(f, barCoords);

    return SUCCESS;
}

SurfaceRouter::RetCode
SurfaceRouter::refineVolumeToAllowReroute(set<CH>& space, const set<VH>& boundaryVs, const set<EH>& boundaryEs)
{
    auto& tetMesh = meshProps().mesh();

    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());

#define SPLIT_IF_ALL_VS_FORBIDDEN(FORBIDDENSET, SPLITSET, ELEMENT, VERTEXRANGE)                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        if (FORBIDDENSET.count(ELEMENT) == 0                                                                           \
            && !containsMatching(tetMesh.VERTEXRANGE(ELEMENT),                                                         \
                                 [&](const VH& v) { return boundaryVs.count(v) == 0 && _forbiddenVs.count(v) == 0; })) \
            SPLITSET.insert(ELEMENT);                                                                                  \
    } while (0)

    set<CH> dummy;
    set<CH> splitTets;
    set<FH> splitFs;
    set<EH> splitEs;
    for (CH tet : space)
        for (EH e : tetMesh.cell_edges(tet))
            if (boundaryEs.count(e) == 0)
                SPLIT_IF_ALL_VS_FORBIDDEN(_forbiddenEs, splitEs, e, edge_vertices);
    for (CH tet : space)
        for (FH f : tetMesh.cell_faces(tet))
            if (!containsSomeOf(tetMesh.face_edges(f), splitEs))
                SPLIT_IF_ALL_VS_FORBIDDEN(_forbiddenFs, splitFs, f, face_vertices);
    for (CH tet : space)
        if (!containsSomeOf(tetMesh.cell_edges(tet), splitEs) && !containsSomeOf(tetMesh.cell_faces(tet), splitFs))
            SPLIT_IF_ALL_VS_FORBIDDEN(dummy, splitTets, tet, tet_vertices);

    for (CH tet : splitTets)
        splitTet(tet, {Q(1, 4), Q(1, 4), Q(1, 4), Q(1, 4)});
    for (FH f : splitFs)
        splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});
    for (EH e : splitEs)
        if (!containsSomeOf(tetMesh.edge_faces(e), _forbiddenFs))
            splitHalfEdge(tetMesh.halfedge_handle(e, 0), *tetMesh.ec_iter(e), 0.5);

    list<CH> children(space.begin(), space.end());
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
    space = {children.begin(), children.end()};

#undef SPLIT_IF_ALL_VS_FORBIDDEN
    return SUCCESS;
}

SurfaceRouter::RetCode SurfaceRouter::refineVolumeToAvoidComplex(set<CH>& space,
                                                                 const set<EH>& complexEs,
                                                                 const set<VH>& complexVs,
                                                                 const set<HFH>& surface)
{
    auto& tetMesh = meshProps().mesh();

    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());

    (void)surface;
    set<EH> splitEs;

    for (EH e1 : complexEs)
        for (VH v : tetMesh.edge_vertices(e1))
            for (EH e : tetMesh.vertex_edges(v))
                if (e != e1)
                    splitEs.insert(e);

    for (VH v : complexVs)
        for (EH e : tetMesh.vertex_edges(v))
            splitEs.insert(e);

    for (EH e : splitEs)
        if (!containsSomeOf(tetMesh.edge_faces(e), _forbiddenFs) && _forbiddenEs.count(e) == 0)
            splitHalfEdge(tetMesh.halfedge_handle(e, 0), *tetMesh.ec_iter(e), 0.5);

    list<CH> children(space.begin(), space.end());
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
    space = {children.begin(), children.end()};

    return SUCCESS;
}

SurfaceRouter::RetCode SurfaceRouter::refineToAvoidNonmanifoldVertex(const HFH& hf, set<CH>* space)
{
    auto& tetMesh = meshProps().mesh();

    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());

    set<EH> esToCut;
    set<FH> fsCut;

    for (HEH he : tetMesh.halfface_halfedges(hf))
    {
        HFH hfAdj = tetMesh.adjacent_halfface_in_cell(hf, he);
        HEH heNext = tetMesh.next_halfedge_in_halfface(tetMesh.opposite_halfedge_handle(he), hfAdj);
        if (_forbiddenVs.find(tetMesh.from_vertex_handle(heNext)) != _forbiddenVs.end())
        {
            esToCut.insert(tetMesh.edge_handle(heNext));
            fsCut.insert(tetMesh.face_handle(hfAdj));
            fsCut.insert(tetMesh.face_handle(tetMesh.adjacent_halfface_in_cell(hfAdj, heNext)));
        }
    }

    set<FH> fsToCut;
    for (HEH he : tetMesh.halfface_halfedges(hf))
    {
        FH fAdj = tetMesh.face_handle(tetMesh.adjacent_halfface_in_cell(hf, he));
        if (_forbiddenEs.count(tetMesh.edge_handle(he)) != 0 && fsCut.count(fAdj) == 0)
            fsToCut.insert(fAdj);
    }

    for (EH e : esToCut)
    {
        HEH he = tetMesh.halfedge_handle(e, 0);
        splitHalfEdge(he, *tetMesh.hec_iter(he), 0.5);
    }

    for (FH f : fsToCut)
        splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});

    if (esToCut.empty() && fsToCut.empty())
        splitTet(tetMesh.incident_cell(hf), {Q(1, 4), Q(1, 4), Q(1, 4), Q(1, 4)});

    if (space != nullptr)
    {
        list<CH> children(space->begin(), space->end());
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
        *space = {children.begin(), children.end()};
    }

    return SUCCESS;
}

SurfaceRouter::RetCode
SurfaceRouter::refineToAvoidNonmanifoldEdge(const CH& tet, const set<HFH>& currentHfs, set<CH>* space)
{
    auto& tetMesh = meshProps().mesh();

    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS> propGuard(meshProps());

    set<FH> fsToCut;
    for (HFH hf : tetMesh.cell_halffaces(tet))
        for (HEH he : tetMesh.halfface_halfedges(hf))
        {
            HFH hfAdj = tetMesh.adjacent_halfface_in_cell(hf, he);
            if (currentHfs.count(hfAdj) != 0)
                continue;
            FH fAdj = tetMesh.face_handle(hfAdj);
            if (_forbiddenEs.count(tetMesh.edge_handle(he)) != 0
                || _forbiddenVs.count(tetMesh.from_vertex_handle(he)) != 0
                || _forbiddenVs.count(tetMesh.to_vertex_handle(he)))
                fsToCut.insert(fAdj);
        }
    for (FH f : fsToCut)
        splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});
    if (fsToCut.empty())
        splitTet(tet, {Q(1, 4), Q(1, 4), Q(1, 4), Q(1, 4)});

    if (space != nullptr)
    {
        list<CH> children(space->begin(), space->end());
        for (auto it = children.begin(); it != children.end();)
        {
            CH tet2 = *it;
            if (tetMesh.is_deleted(tet2))
            {
                for (CH child : meshProps().get<CHILD_CELLS>(tet2))
                    children.push_back(child);
                children.erase(it++);
            }
            else
                it++;
        }
        *space = {children.begin(), children.end()};
    }

    return SUCCESS;
}

set<CH> SurfaceRouter::confineVolume(const set<VH>& vsBoundary, const set<CH>& volume) const
{
    auto& tetMesh = meshProps().mesh();

    set<CH> confinedVolume;

    Vec3Q bboxMin;
    Vec3Q bboxMax;

    CH seedTet = findSomeOf(tetMesh.vertex_cells(*vsBoundary.begin()), volume);
    assert(seedTet.is_valid());
    {
        map<CH, Transition> tet2trans({{seedTet, Transition()}});
        list<pair<CH, Transition>> tetQ({{seedTet, Transition()}});
        while (!tetQ.empty())
        {
            auto tetTrans = tetQ.front();
            tetQ.pop_front();

            for (HFH hf : tetMesh.cell_halffaces(tetTrans.first))
            {
                if (!containsSomeOf(tetMesh.halfface_vertices(hf), vsBoundary))
                    continue;
                CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                if (tetNext.is_valid() && tet2trans.count(tetNext) == 0 && volume.count(tetNext) != 0)
                {
                    Transition transNext = tetTrans.second.chain(meshProps().hfTransition<TRANSITION>(hf));
                    tet2trans.insert({tetNext, transNext});
                    tetQ.push_back({tetNext, transNext});
                }
            }
        }
        bboxMin = Vec3Q(DBL_MAX, DBL_MAX, DBL_MAX);
        bboxMax = Vec3Q(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        for (VH v : vsBoundary)
        {
            for (CH tet : tetMesh.vertex_cells(v))
            {
                auto it = tet2trans.find(tet);
                if (it != tet2trans.end())
                {
                    Vec3Q uvw = it->second.invert().apply(meshProps().ref<CHART>(it->first).at(v));

                    for (int i = 0; i < 3; i++)
                    {
                        if (bboxMin[i] > uvw[i])
                            bboxMin[i] = uvw[i];
                        if (bboxMax[i] < uvw[i])
                            bboxMax[i] = uvw[i];
                    }
                }
            }
        }
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
                        [&](const VH& v)
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
    return confinedVolume;
}

void SurfaceRouter::setForbiddenElements(const CH& b)
{
    auto& tetMesh = meshProps().mesh();

    _forbiddenFs.clear();
    _forbiddenEs.clear();
    _forbiddenVs.clear();

    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
        for (HFH hf : tetMesh.cell_halffaces(tet))
            if (meshProps().isBlockBoundary(hf))
                _forbiddenFs.insert(tetMesh.face_handle(hf));

    for (FH f : _forbiddenFs)
    {
        for (EH e : tetMesh.face_edges(f))
            _forbiddenEs.insert(e);
        for (VH v : tetMesh.face_vertices(f))
            _forbiddenVs.insert(v);
    }
}

void SurfaceRouter::getSurfaceBoundary(const set<HFH>& surface,
                                       set<VH>& vsBoundary,
                                       set<EH>& esBoundary,
                                       set<HEH>& hesBoundary) const
{
    auto& tetMesh = meshProps().mesh();
    vsBoundary.clear();
    esBoundary.clear();
    hesBoundary.clear();
    for (HFH hf : surface)
        for (HEH he : tetMesh.halfface_halfedges(hf))
        {
            EH e = tetMesh.edge_handle(he);
            auto it = esBoundary.find(e);
            if (it != esBoundary.end())
            {
                auto it2 = hesBoundary.find(tetMesh.opposite_halfedge_handle(he));
                if (it2 == hesBoundary.end())
                {
                    LOG(ERROR) << "Non-manifold surface or inverse halffaces inside the surface";
                    throw std::logic_error("");
                }
                assert(it2 != hesBoundary.end());
                hesBoundary.erase(it2);
                esBoundary.erase(it);
            }
            else
            {
                hesBoundary.insert(he);
                esBoundary.insert(e);
            }
        }
    for (EH e : esBoundary)
    {
        auto vs = tetMesh.edge_vertices(e);
        vsBoundary.insert(vs.begin(), vs.end());
    }
}

} // namespace c4hex
