#include "C4Hex/Algorithm/MCSmoother.hpp"

#include "C4Hex/Algorithm/PathRouter.hpp"
#include "C4Hex/Algorithm/SurfaceRouter.hpp"
#include <MC3D/Algorithm/TetRemesher.hpp>

namespace c4hex
{

MCSmoother::MCSmoother(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

bool MCSmoother::isUVWAligned(const EH& a) const
{
    auto& tetMesh = meshProps().mesh();

    CH b = *mcMeshProps().mesh().ec_iter(a);
    UVWDir dir = UVWDir::NONE;
    return !containsMatching(mcMeshProps().ref<ARC_MESH_HALFEDGES>(a),
                             [&, this](const HEH& he)
                             {
                                 CH tetB = findMatching(tetMesh.halfedge_cells(he),
                                                        [&, this](const CH& tet)
                                                        { return meshProps().get<MC_BLOCK>(tet) == b; });
                                 if (tetB.is_valid())
                                     dir = dir | edgeDirection(tetMesh.edge_handle(he), tetB);
                                 return dim(dir) > 1;
                             });
}

bool MCSmoother::isUVWAligned(const FH& p) const
{
    UVWDir dir = UVWDir::NONE;
    return !containsMatching(mcMeshProps().get<PATCH_MESH_HALFFACES>(p),
                             [&, this](const HFH& hf)
                             {
                                 dir = dir | normalDirUVW(hf);
                                 return dim(dir) != 1;
                             });
}

void MCSmoother::smoothMC()
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();
    TetRemesher remesher(meshProps());
    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS, CHILD_FACES, CHILD_EDGES, CHILD_HALFFACES, CHILD_HALFEDGES>
        propGuard(meshProps());

    set<FH> psOptimized;
    list<EH> as;
    for (EH a : mcMesh.edges())
        as.push_back(a);
    set<EH> asOnQ(as.begin(), as.end());
    map<EH, int> a2timesOptimized;
    while (!as.empty())
    {
        EH a = as.front();
        as.pop_front();
        asOnQ.erase(a);

        if (a2timesOptimized[a] >= 3 || mcMeshProps().get<IS_SINGULAR>(a)
            || (mcMeshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps().get<IS_FEATURE_E>(a))
            || (!mcMesh.is_boundary(a) && mcMeshProps().isAllocated<IS_FEATURE_F>()
                && containsMatching(mcMesh.edge_faces(a),
                                    [this](const FH& p) { return mcMeshProps().get<IS_FEATURE_F>(p); }))
            || isUVWAligned(a))
            continue;

        if (tetMesh.n_logical_cells() > 2000000)
        {
            remesher.collapseAllPossibleEdges(false, true, true, true, 10.0);
            remesher.remeshToImproveAngles(true, false, TetRemesher::QualityMeasure::ANGLES);
        }

        LOG(INFO) << "Smoothing arc " << a;

        _patchesRerouted.clear();
        _blocksChanged.clear();
        _traversedHfs.clear();
        _aReroute = a;
        collectAllowedRegionAroundArc(a);

        auto hesA = mcMeshProps().get<ARC_MESH_HALFEDGES>(a);

        bool change = false;
        auto ret = reroute(a, &change);
        if (ret != SUCCESS)
        {
            updateArcEmbedding(a, hesA);
            assertValidMC(true, true);
            continue;
        }
        if (!change)
        {
            LOG(INFO) << "Arc " << a << " already as smooth as possible";
            continue;
        }
        map<FH, set<HFH>> p2hfs;
        map<FH, Transition> p2trans;
        for (FH p : mcMesh.edge_faces(a))
        {
            p2trans[p] = mcMeshProps().ref<PATCH_TRANSITION>(p);
            p2hfs[p] = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);
        }

        // Order patches incident on a starting from hp
        auto psToReroute = cyclicOrderPatches(*mcMesh.ehf_iter(a), a);

        // for all ordered patches p incident on a (circulating around a)
        for (FH p : psToReroute)
        {
            if (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p)
                && !mcMesh.is_boundary(p))
                continue;
            DLOG(INFO) << "Smoothing patch " << p;
            ret = reroute(p);
            if (ret != SUCCESS)
                break;
            DLOG(INFO) << "Smoothed patch " << p;
        }

        if (ret != SUCCESS)
        {
            updateArcEmbedding(a, hesA);
            for (auto& kv : p2hfs)
                updatePatchEmbedding(kv.first, kv.second, &p2trans.at(kv.first));
            assertValidMC(true, true);
            continue;
        }

        map<CH, set<CH>> b2tets;
        set<CH> bs;
        for (FH p : mcMesh.edge_faces(a))
            for (CH b : mcMesh.face_cells(p))
                if (b.is_valid())
                    bs.insert(b);
        for (CH b : bs)
            b2tets[b] = mcMeshProps().ref<BLOCK_MESH_TETS>(b);
        _blocksChanged = bs;
        if (refloodFillBlocks() != SUCCESS)
        {
            updateArcEmbedding(a, hesA);
            for (auto& kv : p2hfs)
                updatePatchEmbedding(kv.first, kv.second, &p2trans.at(kv.first));
            refloodFillBlocks();
            assertValidMC(false, false);
            continue;
        }

        LOG(INFO) << "Successfully smoothed arc " << a;
        assertValidMC(true, true);
        a2timesOptimized[a]++;
        for (FH p : psToReroute)
        {
            psOptimized.insert(p);
            for (EH aNext : mcMesh.face_edges(p))
            {
                if (aNext != a && asOnQ.count(aNext) == 0 && a2timesOptimized[a] < 3)
                {
                    as.push_back(aNext);
                    asOnQ.insert(aNext);
                }
            }
        }
    }

    for (FH p : mcMesh.faces())
    {
        if (mcMesh.is_boundary(p) || psOptimized.count(p) != 0
            || (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p)) || isUVWAligned(p))
            continue;

        _patchesRerouted.clear();
        _blocksChanged.clear();
        _traversedHfs.clear();
        _aReroute = EH(-1);
        for (auto a : mcMesh.face_edges(p))
        {
            if (a2timesOptimized[a] >= 3 || mcMeshProps().get<IS_SINGULAR>(a)
                || (mcMeshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps().get<IS_FEATURE_E>(a))
                || (!mcMesh.is_boundary(a) && mcMeshProps().isAllocated<IS_FEATURE_F>()
                    && containsMatching(mcMesh.edge_faces(a),
                                        [this](const FH& p2) { return mcMeshProps().get<IS_FEATURE_F>(p2); }))
                || isUVWAligned(a))
                continue;
            _aReroute = a;
            break;
        }
        if (!_aReroute.is_valid())
            continue;
        EH a = _aReroute;
        collectAllowedRegionAroundArc(_aReroute);

        auto hesA = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);

        bool change = false;
        auto ret = reroute(a, &change);
        if (ret != SUCCESS)
        {
            updateArcEmbedding(a, hesA);
            assertValidMC(true, true);
            continue;
        }
        map<FH, set<HFH>> p2hfs;
        map<FH, Transition> p2trans;
        for (FH p2 : mcMesh.edge_faces(a))
        {
            p2trans[p2] = mcMeshProps().ref<PATCH_TRANSITION>(p2);
            p2hfs[p2] = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p2);
        }

        // Order patches incident on a starting from hp
        auto psToReroute = cyclicOrderPatches(*mcMesh.ehf_iter(a), a);

        // for all ordered patches p incident on a (circulating around a)
        for (FH p2 : psToReroute)
        {
            if (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p)
                && !mcMesh.is_boundary(p2))
                continue;
            LOG(INFO) << "Smoothing patch " << p2;
            ret = reroute(p2);
            if (ret != SUCCESS)
                break;
            psOptimized.insert(p2);
            LOG(INFO) << "Smoothed patch " << p2;
        }

        if (ret != SUCCESS)
        {
            updateArcEmbedding(a, hesA);
            for (auto& kv : p2hfs)
                updatePatchEmbedding(kv.first, kv.second, &p2trans.at(kv.first));
            assertValidMC(true, true);
            continue;
        }

        map<CH, set<CH>> b2tets;
        set<CH> bs;
        for (FH p2 : mcMesh.edge_faces(a))
            for (CH b : mcMesh.face_cells(p2))
                if (b.is_valid())
                    bs.insert(b);
        for (CH b : bs)
            b2tets[b] = mcMeshProps().ref<BLOCK_MESH_TETS>(b);
        _blocksChanged = bs;
        if (refloodFillBlocks() != SUCCESS)
        {
            updateArcEmbedding(a, hesA);
            for (auto& kv : p2hfs)
                updatePatchEmbedding(kv.first, kv.second, &p2trans.at(kv.first));
            refloodFillBlocks();
            assertValidMC(true, true);
            continue;
        }
        assertValidMC(true, true);
    }
}

MCSmoother::RetCode MCSmoother::makeVolumeTransitionFree(const set<CH>& space, const CH& tetSeed)
{
    set<CH> tetVisited({tetSeed});
    list<CH> tetStack({tetSeed});

#ifndef NDEBUG
    set<FH> innerFaces;
#endif

    auto pushTransition = [this,
                           &space,
#ifndef NDEBUG
                           &innerFaces,
#endif
                           &tetVisited](const CH& tet1)
    {
        for (HFH hf1to2 : meshProps().mesh().cell_halffaces(tet1))
        {
            HFH hf2to1 = meshProps().mesh().opposite_halfface_handle(hf1to2);
            CH tet2 = meshProps().mesh().incident_cell(hf2to1);
            if (tet2.is_valid() && space.count(tet2) != 0)
            {
#ifndef NDEBUG
                FH f = meshProps().mesh().face_handle(hf1to2);
                innerFaces.insert(f);
#endif
                if (tetVisited.count(tet2) != 0)
                    continue;
                Transition tr2to1 = meshProps().hfTransition<TRANSITION>(hf2to1);
                for (auto& kv : meshProps().ref<CHART>(tet2))
                {
                    auto& v = kv.first;
                    auto& uvw = kv.second;
                    (void)v;
                    uvw = tr2to1.apply(uvw);
                }
                for (HFH hf2to3 : meshProps().mesh().cell_halffaces(tet2))
                {
                    HFH hf3to2 = meshProps().mesh().opposite_halfface_handle(hf2to3);
                    if (!meshProps().mesh().is_boundary(hf3to2))
                    {
                        Transition tr3to2 = meshProps().hfTransition<TRANSITION>(hf3to2);
                        meshProps().setTransition<TRANSITION>(hf3to2, tr3to2.chain(tr2to1));
                    }
                }
            }
        }
    };

    while (!tetStack.empty())
    {
        CH tet = tetStack.back();
        tetStack.pop_back();

        pushTransition(tet);

        for (HFH hf : meshProps().mesh().cell_halffaces(tet))
        {
            FH f = meshProps().mesh().face_handle(hf);
            CH tetOpp = meshProps().mesh().incident_cell(meshProps().mesh().opposite_halfface_handle(hf));
            if (!meshProps().isBlockBoundary(f) && tetVisited.count(tetOpp) == 0)
            {
                tetVisited.insert(tetOpp);
                tetStack.push_back(tetOpp);
            }
        }
    }
    assert(space == tetVisited);

#ifndef NDEBUG
    for (FH f : innerFaces)
        assert(meshProps().get<TRANSITION>(f).isIdentity());
#endif
    return SUCCESS;
}

void MCSmoother::collectAllowedRegionAroundArc(const EH& aIn)
{
    _allowedVolume.clear();
    _forbiddenFs.clear();
    _forbiddenEs.clear();
    _forbiddenVs.clear();
    _b2unaffectedPs.clear();
    _b2unaffectedAs.clear();
    _b2unaffectedNs.clear();
    _torusSplitter.clear();
    _torusSplitterE.clear();
    _torusSplitterV.clear();
    _a2sectorBack.clear();
    _a2sectorFront.clear();
    _aBoundary2sectorBack.clear();
    _aBoundary2sectorFront.clear();
    _p2boundary.clear();
    _p2sector.clear();

    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    set<CH> bs;
    set<FH> ps;
    set<FH> boundaryPs;
    set<EH> as;
    set<VH> ns;

    for (FH p : mcMesh.edge_faces(aIn))
        for (CH b : mcMesh.face_cells(p))
            if (b.is_valid())
                bs.insert(b);
    for (CH b : bs)
    {
        for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
        {
            if (tetMesh.is_deleted(tet))
                throw std::logic_error("Tet deleted before inserting into allowedspace");
            _allowedVolume.insert(tet);
            if (meshProps().isAllocated<TOUCHED>())
            {
                for (VH v : tetMesh.tet_vertices(tet))
                {
                    meshProps().set<TOUCHED>(v, true);
                    for (VH v2 : tetMesh.vertex_vertices(v))
                        meshProps().set<TOUCHED>(v2, true);
                }
            }
        }
        for (FH p : mcMesh.cell_faces(b))
        {
            ps.insert(p);
            if (!contains(mcMesh.face_edges(p), aIn))
            {
                boundaryPs.insert(p);
                _b2unaffectedPs[b].push_back(p);
                for (HFH element : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                    _forbiddenFs.insert(tetMesh.face_handle(element));
            }
        }
        // assert(!_b2unaffectedPs[b].empty());
        for (EH a : mcMesh.cell_edges(b))
        {
            if (a != aIn)
                _b2unaffectedAs[b].push_back(a);
            as.insert(a);
        }
        for (VH n : mcMesh.cell_vertices(b))
        {
            _b2unaffectedNs[b].push_back(n);
            ns.insert(n);
        }
    }

    for (FH p : ps)
        if (!contains(mcMesh.face_edges(p), aIn))
            for (HFH element : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                _forbiddenFs.insert(tetMesh.face_handle(element));
    for (EH a : as)
        if (a != aIn)
            for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                _forbiddenEs.insert(tetMesh.edge_handle(he));
    for (VH n : ns)
        _forbiddenVs.insert(mcMeshProps().get<NODE_MESH_VERTEX>(n));

    // Check if non-ball topology and register splitting faces in _torusSplitter
    for (auto& kv : _b2unaffectedPs)
    {
        CH b = kv.first;
        auto& ps2 = kv.second;
        for (FH p : ps2)
        {
            auto bsP = mcMesh.face_cells(p);
            CH bNext;
            if (bsP[0] == b)
                bNext = bsP[1];
            else
                bNext = bsP[0];
            if (bs.count(bNext) != 0)
            {
                DLOG(INFO) << "Splitting toroidal allowed space by marking some faces splitters";
                for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                {
                    for (VH v : meshProps().get_halfface_vertices(hf))
                        _torusSplitterV.insert(v);
                    for (EH e : tetMesh.halfface_edges(hf))
                        _torusSplitterE.insert(e);
                    _torusSplitter.insert(tetMesh.face_handle(hf));
                }
            }
        }
    }
    // Check if non-ball topology and register splitting edges in _torusSplitterE
    for (auto& kv : _b2unaffectedAs)
    {
        CH b = kv.first;
        auto& as2 = kv.second;
        for (EH a : as2)
        {
            set<CH> allBs;
            for (FH p : mcMesh.edge_faces(a))
                for (CH b2 : mcMesh.face_cells(p))
                    if (b2.is_valid() && bs.count(b2) != 0)
                        allBs.insert(b2);
            list<CH> bQ({b});
            set<CH> bVisited({b});
            while (!bQ.empty())
            {
                CH b2 = bQ.front();
                bQ.pop_front();

                for (HFH hp : mcMesh.cell_halffaces(b2))
                {
                    if (boundaryPs.count(mcMesh.face_handle(hp)) != 0 || !contains(mcMesh.halfface_edges(hp), a))
                        continue;
                    CH bNext = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp));
                    if (bNext.is_valid() && bVisited.count(bNext) == 0)
                    {
                        assert(bs.count(bNext) != 0);
                        bVisited.insert(bNext);
                        bQ.push_back(bNext);
                    }
                }
            }
            if (bVisited.size() != allBs.size())
            {
                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                {
                    for (VH v : tetMesh.halfedge_vertices(he))
                        _torusSplitterV.insert(v);
                    _torusSplitterE.insert(tetMesh.edge_handle(he));
                }
            }
        }
    }
    // Check if non-ball topology and register splitting vertices in _torusSplitterV
    for (auto& kv : _b2unaffectedNs)
    {
        CH b = kv.first;
        auto& ns2 = kv.second;
        for (VH n : ns2)
        {
            set<CH> allBs;
            for (CH b2 : mcMesh.vertex_cells(n))
                if (bs.count(b2) != 0)
                    allBs.insert(b2);
            list<CH> bQ({b});
            set<CH> bVisited({b});
            while (!bQ.empty())
            {
                CH b2 = bQ.front();
                bQ.pop_front();

                for (HFH hp : mcMesh.cell_halffaces(b2))
                {
                    if (boundaryPs.count(mcMesh.face_handle(hp)) != 0 || !contains(mcMesh.halfface_vertices(hp), n))
                        continue;
                    CH bNext = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp));
                    if (bNext.is_valid() && bVisited.count(bNext) == 0)
                    {
                        assert(bs.count(bNext) != 0);
                        bVisited.insert(bNext);
                        bQ.push_back(bNext);
                    }
                }
            }
            if (bVisited.size() != allBs.size())
                _torusSplitterV.insert(mcMeshProps().ref<NODE_MESH_VERTEX>(n));
        }
    }

    for (FH p : ps)
    {
        set<VH> boundaryVs;
        set<FH> surfaceFs;
        set<HEH> pBoundary;
        for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
            for (HEH he : mcMeshProps().haHalfedges(ha))
                pBoundary.insert(he);
        for (HEH he : pBoundary)
            boundaryVs.insert(tetMesh.from_vertex_handle(he));
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            surfaceFs.insert(tetMesh.face_handle(hf));
        _p2boundary[p] = pBoundary;
        if (!mcMesh.is_boundary(p))
        {
            auto& tetsVisited = _p2sector[p];
            list<CH> tetQ;
            for (VH v : boundaryVs)
                for (FH f : tetMesh.vertex_faces(v))
                    if (surfaceFs.count(f) != 0)
                        for (CH tet : tetMesh.face_cells(f))
                            if (_allowedVolume.count(tet) != 0)
                            {
                                tetsVisited.insert(tet);
                                tetQ.push_back(tet);
                            }
            assert(!tetQ.empty());
            while (!tetQ.empty())
            {
                CH tet = tetQ.front();
                tetQ.pop_front();

                for (HFH hf : tetMesh.cell_halffaces(tet))
                {
                    if (!containsMatching(tetMesh.halfface_vertices(hf),
                                          [&](const VH& v) { return boundaryVs.count(v) != 0; })
                        || _forbiddenFs.count(tetMesh.face_handle(hf)) != 0)
                        continue;
                    CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                    if (tetNext.is_valid() && _allowedVolume.count(tetNext) != 0 && tetsVisited.count(tetNext) == 0)
                    {
                        tetsVisited.insert(tetNext);
                        tetQ.push_back(tetNext);
                    }
                }
            }
            assert(!tetsVisited.empty());
        }
    }
}

list<FH> MCSmoother::cyclicOrderPatches(const HFH& hp, const EH& aReroute) const
{
    auto& mcMesh = mcMeshProps().mesh();

    FH p0 = mcMesh.face_handle(hp);
    list<FH> ps({p0});

    HFH hpFront = hp;
    HFH hpBack = mcMesh.opposite_halfface_handle(hpFront);
    HEH haFront = findMatching(mcMesh.halfface_halfedges(hpFront),
                               [&](const HEH& ha) { return mcMesh.edge_handle(ha) == aReroute; });
    HEH haBack = mcMesh.opposite_halfedge_handle(haFront);

    bool front = false;
    while (mcMesh.incident_cell(hpBack) != mcMesh.incident_cell(hpFront))
    {
        if (mcMesh.is_boundary(hpBack))
            front = true;
        else if (mcMesh.is_boundary(hpFront))
            front = false;
        else
            front = !front;
        (front ? hpFront : hpBack) = mcMesh.opposite_halfface_handle(
            mcMesh.adjacent_halfface_in_cell((front ? hpFront : hpBack), (front ? haFront : haBack)));
        ps.push_back(mcMesh.face_handle(front ? hpFront : hpBack));
    }
    return ps;
}

MCSmoother::RetCode MCSmoother::refloodFillBlocks()
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    for (CH b : _blocksChanged)
    {
        auto& unaffectedPs = _b2unaffectedPs[b];
        FH pSeed;
        if (!unaffectedPs.empty())
        {
            pSeed = unaffectedPs.front();
        }
        else
            pSeed = *mcMesh.cf_iter(b);
        set<FH> ps;
        for (FH p : mcMesh.cell_faces(b))
            ps.insert(p);
        assert(ps.count(pSeed) != 0);
        for (FH p : ps)
            for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                if (tetMesh.is_deleted(hf))
                {
                    LOG(ERROR) << "Floodfilling block failed, because a boundary patch halfface is deleted";
                    return REROUTE_ERROR;
                }

        HFH hpSeed
            = findMatching(mcMesh.cell_halffaces(b), [&](const HFH& hp) { return mcMesh.face_handle(hp) == pSeed; });
        assert(hpSeed.is_valid());

        CH tetSeed = tetMesh.incident_cell(*mcMeshProps().hpHalffaces(hpSeed).begin());
        assert(unaffectedPs.empty() || meshProps().get<MC_BLOCK>(tetSeed) == b);
        set<CH> newBlockTets({tetSeed});
        list<CH> tetQ({tetSeed});

        while (!tetQ.empty())
        {
            CH tet = tetQ.front();
            tetQ.pop_front();

            for (HFH hf : tetMesh.cell_halffaces(tet))
            {
                FH p = meshProps().get<MC_PATCH>(tetMesh.face_handle(hf));
                if (p.is_valid() && ps.count(p) == 0)
                {
                    LOG(ERROR) << "Floodfilling block failed, because block not closed";
                    return REROUTE_ERROR;
                }
                assert(!p.is_valid() || ps.count(p) != 0);
                CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
                if (p.is_valid() || newBlockTets.count(tetNext) != 0)
                    continue;
                if (!tetNext.is_valid())
                {
                    LOG(ERROR) << "Floodfilling block failed, because boundary uncovered";
                    return REROUTE_ERROR;
                }
                newBlockTets.insert(tetNext);
                tetQ.push_back(tetNext);
            }
        }

        makeVolumeTransitionFree(newBlockTets, tetSeed);
        mcMeshProps().set<BLOCK_MESH_TETS>(b, newBlockTets);
        for (CH tet : newBlockTets)
            if (tetMesh.is_deleted(tet))
                throw std::logic_error("Deleted tet floodfilled at collapse");
        for (CH tet : newBlockTets)
        {
            assert(_blocksChanged.count(meshProps().get<MC_BLOCK>(tet)) != 0);
            meshProps().set<MC_BLOCK>(tet, b);
        }
    }
    for (CH b : _blocksChanged)
        for (HFH hp : mcMesh.cell_halffaces(b))
            mcMeshProps().setHpTransition<PATCH_TRANSITION>(
                hp, meshProps().hfTransition<TRANSITION>(*mcMeshProps().hpHalffaces(hp).begin()));

    for (CH b : _blocksChanged)
        for (HFH hp : mcMesh.cell_halffaces(b))
        {
            Transition trans = mcMeshProps().hpTransition<PATCH_TRANSITION>(hp);
            for (HFH hf : mcMeshProps().hpHalffaces(hp))
            {
                if (!(meshProps().hfTransition<TRANSITION>(hf) == trans))
                {
                    LOG(ERROR) << "Floodfilling block failed, because transition doesnt match";
                    return REROUTE_ERROR;
                }
            }
        }
    return SUCCESS;
}

MCSmoother::RetCode MCSmoother::reroute(const EH& a, bool* change)
{
    _aReroute = a;
    if (change != nullptr)
        *change = false;

    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    bool aIsBoundary = mcMesh.is_boundary(a);
    bool aIsSingular = mcMeshProps().get<IS_SINGULAR>(a);

    auto oldHes = mcMeshProps().get<ARC_MESH_HALFEDGES>(a);

    list<HEH> pathRerouted = oldHes;

    VH vFrom = tetMesh.from_vertex_handle(pathRerouted.front());
    VH vTo = tetMesh.to_vertex_handle(pathRerouted.back());
    // Force to keep first halfedge at non-shifted end

    // Reroute only if not singular (if singular, just keep the path resulting after appending)
    if (!aIsSingular)
    {
        // allowedVolumeArc = aIsBoundary ? allowedVolume : boundaryFaces of allowedVolume
        if (aIsBoundary)
        {
            set<FH> forbiddenFsArc = _forbiddenFs;
            set<EH> forbiddenEsArc = _forbiddenEs;
            set<VH> forbiddenVsArc = _forbiddenVs;
            set<FH> allowedFacesArc;
            for (CH tet : _allowedVolume)
                if (tetMesh.is_deleted(tet))
                    throw std::logic_error("Deleted tet");
            for (CH tet : _allowedVolume)
                for (FH f : tetMesh.cell_faces(tet))
                    if (tetMesh.is_deleted(f))
                        throw std::logic_error("Deleted face");
            for (CH tet : _allowedVolume)
                for (FH f : tetMesh.cell_faces(tet))
                    if (tetMesh.is_boundary(f))
                        allowedFacesArc.insert(f);
            for (FH f : forbiddenFsArc)
                for (EH e : tetMesh.face_edges(f))
                    forbiddenEsArc.insert(e);
            for (EH e : forbiddenEsArc)
                for (VH v : tetMesh.edge_vertices(e))
                    forbiddenVsArc.insert(v);

            forbiddenVsArc.erase(vFrom);
            forbiddenVsArc.erase(vTo);
            if (!_aBoundary2sectorFront.empty() && !_aBoundary2sectorBack.empty())
            {
                for (FH f : tetMesh.vertex_faces(vFrom))
                    if (_aBoundary2sectorFront[a].count(f) == 0 && _aBoundary2sectorBack[a].count(f) == 0)
                        allowedFacesArc.erase(f);
                for (FH f : tetMesh.vertex_faces(vTo))
                    if (_aBoundary2sectorFront[a].count(f) == 0 && _aBoundary2sectorBack[a].count(f) == 0)
                        allowedFacesArc.erase(f);
            }

            // meshedges[a] = shortest path (from[a] -> nTo) through allowedVolumeArc
            // traversedfaces[a] = minimal surface (boundary = oldHes + meshedges[a] + collapseedges
            if (PathRouter(meshProps())
                    .reroutePathThroughSurface(
                        pathRerouted, allowedFacesArc, forbiddenFsArc, forbiddenEsArc, forbiddenVsArc, _traversedHfs[a])
                != PathRouter::SUCCESS)
            {
                LOG(ERROR) << "Could not determine shortest arc path through boundary";
                return REROUTE_ERROR;
            }
        }
        else
        {
            set<CH> allowedVolumeArc = _allowedVolume;
            set<FH> forbiddenFsArc = _forbiddenFs;
            set<EH> forbiddenEsArc = _forbiddenEs;
            set<VH> forbiddenVsArc = _forbiddenVs;
            for (FH f : forbiddenFsArc)
                for (EH e : tetMesh.face_edges(f))
                    forbiddenEsArc.insert(e);
            for (EH e : forbiddenEsArc)
                for (VH v : tetMesh.edge_vertices(e))
                    forbiddenVsArc.insert(v);

            for (CH tet : allowedVolumeArc)
            {
                for (FH f : tetMesh.cell_faces(tet))
                    if (tetMesh.is_boundary(f))
                        forbiddenFsArc.insert(f);
                for (EH e : tetMesh.cell_edges(tet))
                    if (tetMesh.is_boundary(e))
                        forbiddenEsArc.insert(e);
                for (VH v : tetMesh.cell_vertices(tet))
                    if (tetMesh.is_boundary(v))
                        forbiddenVsArc.insert(v);
            }

            forbiddenVsArc.erase(vFrom);
            forbiddenVsArc.erase(vTo);

            if (!_a2sectorFront.empty() && !_a2sectorBack.empty())
            {
                for (CH tet : tetMesh.vertex_cells(vFrom))
                    if (_a2sectorFront[a].count(tet) == 0 && _a2sectorBack[a].count(tet) == 0)
                        allowedVolumeArc.erase(tet);
                for (CH tet : tetMesh.vertex_cells(vTo))
                    if (_a2sectorFront[a].count(tet) == 0 && _a2sectorBack[a].count(tet) == 0)
                        allowedVolumeArc.erase(tet);
            }

            // meshedges[a] = shortest path (from[a] -> nTo) through allowedVolumeArc
            // traversedfaces[a] = minimal surface (boundary = oldHes + meshedges[a] + collapseedges
            if (PathRouter(meshProps())
                    .reroutePathThroughVolume(
                        pathRerouted, allowedVolumeArc, forbiddenFsArc, forbiddenEsArc, forbiddenVsArc)
                != PathRouter::SUCCESS)
            {
                LOG(ERROR) << "Could not determine shortest arc path through volume";
                return REROUTE_ERROR;
            }
        }

        // Replace deleted tets in allowedVolume, faces in forbiddenFs, halfedges in oldHes
        meshProps().replaceAllByChildren(_allowedVolume, _forbiddenFs, _forbiddenEs, oldHes);
        for (auto& kv : _traversedHfs)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _a2sectorFront)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _a2sectorBack)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _aBoundary2sectorFront)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _aBoundary2sectorBack)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _p2sector)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _p2boundary)
            meshProps().replaceByChildren(kv.second);

        for (auto it = oldHes.begin(); it != oldHes.end();)
        {
            HEH he = *it;
            if (tetMesh.is_deleted(he))
            {
                oldHes.erase(it++);
                auto children = mcMeshProps().get<CHILD_HALFEDGES>(he);
                it = oldHes.insert(it, children.begin(), children.end());
            }
            else
                it++;
        }
    }

    pathRerouted = pathRerouted;

    if (!aIsSingular)
        for (HEH he : pathRerouted)
        {
            if (_forbiddenEs.count(tetMesh.edge_handle(he)) != 0)
            {
                LOG(ERROR) << "Arc rerouted through forbidden edge";
                return REROUTE_ERROR;
            }
            if (he != pathRerouted.back())
            {
                VH v = tetMesh.to_vertex_handle(he);
                if (_forbiddenVs.count(v) != 0)
                {
                    LOG(ERROR) << "Arc rerouted through forbidden vertex";
                    return REROUTE_ERROR;
                }
            }
        }

    updateArcEmbedding(a, pathRerouted);

    // Forbid rerouted arc:
    for (HEH he : pathRerouted)
        _forbiddenEs.insert(tetMesh.edge_handle(he));

    if (change != nullptr && oldHes != pathRerouted)
        *change = true;

    return SUCCESS;
}

MCSmoother::RetCode MCSmoother::reroute(const FH& p)
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    if (!mcMesh.is_boundary(p))
    {
        set<VH> boundaryVs;
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            for (EH e : tetMesh.halfface_edges(hf))
            {
                int n = 0;
                for (HFH hf2 : tetMesh.edge_halffaces(e))
                    if (mcMeshProps().ref<PATCH_MESH_HALFFACES>(p).count(hf2) != 0)
                        n++;
                if (n != 2)
                    for (VH v : tetMesh.edge_vertices(e))
                        boundaryVs.insert(v);
            }
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            for (VH v : meshProps().get_halfface_vertices(hf))
                if (boundaryVs.count(v) == 0 && tetMesh.is_boundary(v))
                    throw std::logic_error("Non-boundary patch touches boundary before rerouting");
    }

    bool pIsBoundary = mcMesh.is_boundary(p);

    set<HEH> pBoundary;
    for (HEH ha : mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0)))
    {
        for (HEH he : mcMeshProps().haHalfedges(ha))
        {
            auto it = pBoundary.find(tetMesh.opposite_halfedge_handle(he));
            if (it != pBoundary.end())
                pBoundary.erase(it);
            else
                pBoundary.insert(he);
        }
    }

    map<VH, int> incidenceIn, incidenceOut;
    for (HEH he : pBoundary)
    {
        incidenceIn[tetMesh.to_vertex_handle(he)]++;
        incidenceOut[tetMesh.from_vertex_handle(he)]++;
    }
    for (auto* collPtr : {&incidenceIn, &incidenceOut})
        for (auto kv : *collPtr)
            if (kv.second != 1)
            {
                LOG(WARNING) << "Nonmanifold patch boundary, can not guarantee manifoldness of surface at v "
                             << kv.first;
            }

    auto oldPatchHalffaces = mcMeshProps().get<PATCH_MESH_HALFFACES>(p);
    auto newPatchHalffaces = oldPatchHalffaces;

    if (pIsBoundary)
    {
        assert(mcMesh.is_boundary(_aReroute));
        set<HFH> finalTraversedHfs;
        for (EH a : mcMesh.face_edges(p))
        {
            if (!mcMesh.is_boundary(a))
                continue;
            auto& traversedHfs = _traversedHfs[a];
            if (!traversedHfs.empty())
                DLOG(INFO) << "Inserting/erasing traversedhfs of arc " << a << "into/from patch " << p;
            for (HFH hf : traversedHfs)
            {
                bool flip = tetMesh.is_boundary(*oldPatchHalffaces.begin()) != tetMesh.is_boundary(hf);
                bool isShrinking = finalTraversedHfs.count(flip ? tetMesh.opposite_halfface_handle(hf) : hf) != 0;
                if (isShrinking)
                    finalTraversedHfs.erase(flip ? tetMesh.opposite_halfface_handle(hf) : hf);
                else
                    finalTraversedHfs.insert(flip ? tetMesh.opposite_halfface_handle(hf) : hf);
            }
        }
        for (HFH hf : finalTraversedHfs)
        {
            if (oldPatchHalffaces.count(hf) != 0)
                newPatchHalffaces.erase(hf);
            else
            {
                if (_forbiddenFs.count(tetMesh.face_handle(hf)) != 0)
                {
                    LOG(ERROR) << "Inserting forbidden traversed hfs into patch " + std::to_string(p.idx());
                    return REROUTE_ERROR;
                }
                newPatchHalffaces.insert(hf);
            }
        }
    }
    else
    {
        // meshfaces[p] = minimal surface (connecting boundary[p])
        //  through allowedVolumePatch
        auto allowedVolumePatch = _allowedVolume;
        auto forbiddenFsPatch = _forbiddenFs;

        set<EH> boundaryEs;
        for (HEH he : pBoundary)
        {
            assert(boundaryEs.count(tetMesh.edge_handle(he)) == 0);
            boundaryEs.insert(tetMesh.edge_handle(he));
        }
        set<VH> boundaryVs;
        for (HEH he : pBoundary)
            boundaryVs.insert(tetMesh.from_vertex_handle(he));

        set<VH> oldBoundaryVs;
        for (HEH he : _p2boundary.at(p))
            oldBoundaryVs.insert(tetMesh.from_vertex_handle(he));

        for (VH v : oldBoundaryVs)
            for (CH tet : tetMesh.vertex_cells(v))
                if (_p2sector.at(p).count(tet) == 0)
                    allowedVolumePatch.erase(tet);

        for (CH tet : allowedVolumePatch)
            for (FH f : tetMesh.cell_faces(tet))
                if (tetMesh.is_boundary(f))
                    forbiddenFsPatch.insert(f);

        set<EH> forbiddenEsPatch;
        set<VH> forbiddenVsPatch;
        for (EH e : _forbiddenEs)
        {
            if (boundaryEs.count(e) == 0)
            {
                forbiddenEsPatch.insert(e);
                for (VH v : tetMesh.edge_vertices(e))
                    if (boundaryVs.count(v) == 0)
                        forbiddenVsPatch.insert(v);
            }
        }
        for (VH v : _forbiddenVs)
            if (boundaryVs.count(v) == 0)
                forbiddenVsPatch.insert(v);
        for (FH f : forbiddenFsPatch)
        {
            for (EH e : tetMesh.face_edges(f))
                if (boundaryEs.count(e) == 0)
                    forbiddenEsPatch.insert(e);
            for (VH v : tetMesh.face_vertices(f))
                if (boundaryVs.count(v) == 0)
                    forbiddenVsPatch.insert(v);
        }

        if (SurfaceRouter(meshProps())
                .calcMinimalSurfaceByLP(allowedVolumePatch,
                                        forbiddenFsPatch,
                                        forbiddenEsPatch,
                                        forbiddenVsPatch,
                                        boundaryEs,
                                        pBoundary,
                                        boundaryVs,
                                        newPatchHalffaces)
            != SurfaceRouter::SUCCESS)
        {
            LOG(ERROR) << "Could not determine minimal surface by LP for patch " << p;
            return REROUTE_ERROR;
        }

        // Replace deleted tets in allowedVolume, faces in forbiddenFs, halfedges in oldHes
        meshProps().replaceAllByChildren(_allowedVolume, _forbiddenFs, _forbiddenEs);
        for (auto& kv : _traversedHfs)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _a2sectorFront)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _a2sectorBack)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _aBoundary2sectorFront)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _aBoundary2sectorBack)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _p2sector)
            meshProps().replaceByChildren(kv.second);
        for (auto& kv : _p2boundary)
            meshProps().replaceByChildren(kv.second);

        for (HFH hf : newPatchHalffaces)
            for (VH v : meshProps().get_halfface_vertices(hf))
                if (boundaryVs.count(v) == 0 && tetMesh.is_boundary(v))
                {
                    LOG(ERROR) << "Non-boundary patch touches boundary after rerouting";
                    return REROUTE_ERROR;
                }
    }

    for (HFH hf : newPatchHalffaces)
    {
        FH p2 = meshProps().get<MC_PATCH>(tetMesh.face_handle(hf));
        if (p2.is_valid() && _patchesRerouted.count(p2) != 0)
        {
            LOG(ERROR) << "Patch " + std::to_string(p.idx()) + " rerouted onto an already rerouted patch "
                              + std::to_string(p2.idx());
            return REROUTE_ERROR;
        }
        if (p2.is_valid() && !contains(mcMesh.face_edges(p2), _aReroute))
        {
            LOG(ERROR) << "Pillow patch rerouted onto an unaffected patch";
            return REROUTE_ERROR;
        }
    }

    updatePatchEmbedding(p, newPatchHalffaces);

    for (CH b : mcMesh.face_cells(p))
        if (b.is_valid())
            _blocksChanged.insert(b);

    for (HFH hf : newPatchHalffaces)
        if (_forbiddenFs.count(tetMesh.face_handle(hf)) != 0)
        {
            LOG(ERROR) << "Patch rerouted through forbidden face";
            return REROUTE_ERROR;
        }

    // allowedVolumePatch = allowedVolumePatch - meshfaces[p]
    for (HFH hf : newPatchHalffaces)
        _forbiddenFs.insert(tetMesh.face_handle(hf));

    _patchesRerouted.insert(p);

    return SUCCESS;
}

void MCSmoother::updateArcEmbedding(const EH& a, list<HEH>& hesA)
{
    auto& tetMesh = meshProps().mesh();

    meshProps().replaceByChildren(hesA);
    for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
    {
        if (meshProps().get<MC_ARC>(tetMesh.edge_handle(he)) == a)
        {
            meshProps().reset<MC_ARC>(tetMesh.edge_handle(he));
            meshProps().reset<IS_ARC>(tetMesh.edge_handle(he));
        }
    }
    mcMeshProps().set<ARC_MESH_HALFEDGES>(a, hesA);
    for (HEH he : hesA)
    {
        meshProps().set<MC_ARC>(tetMesh.edge_handle(he), a);
        meshProps().set<IS_ARC>(tetMesh.edge_handle(he), true);
    }
}

void MCSmoother::updatePatchEmbedding(const FH& p, set<HFH>& hfsP, Transition* transP)
{
    auto& tetMesh = meshProps().mesh();

    meshProps().replaceByChildren(hfsP);
    for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
    {
        if (meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)) == p)
        {
            meshProps().reset<MC_PATCH>(tetMesh.face_handle(hf));
            meshProps().reset<IS_WALL>(tetMesh.face_handle(hf));
            if (transP != nullptr)
                meshProps().reset<TRANSITION>(tetMesh.face_handle(hf));
        }
    }
    mcMeshProps().set<PATCH_MESH_HALFFACES>(p, hfsP);
    for (HFH hf : hfsP)
    {
        meshProps().set<MC_PATCH>(tetMesh.face_handle(hf), p);
        meshProps().set<IS_WALL>(tetMesh.face_handle(hf), true);
        if (transP != nullptr)
            meshProps().setTransition<TRANSITION>(hf, *transP);
    }
}

} // namespace c4hex
