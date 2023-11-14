#include "C4Hex/Algorithm/EmbeddingCollapser.hpp"

#include <list>

namespace c4hex
{

EmbeddingCollapser::EmbeddingCollapser(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

void EmbeddingCollapser::resetVars()
{
    _layerFaces.clear();
    _layerTets.clear();
    _collapseHes.clear();

    _haCollapse = {};
    _pCollapse = {};
    _haReroute = {};
    _haRerouteInto = {};
    _heCurrent = {};
    _vPullUp = {};

    _aToAppend = {};

    _a2pToAppend.clear();
}

void EmbeddingCollapser::setVars(const HEH& ha)
{
    resetVars();
    _a2pToAppend.clear();
    _haCollapse = ha;
    _aToAppend = determineLastSuccessorArc();
    _collapseHes = mcMeshProps().haHalfedges(_haCollapse);
}

void EmbeddingCollapser::setVars(const FH& p, const HEH& haMoving, const HEH& haStationary)
{
    resetVars();
    _pCollapse = p;
    _haReroute = haMoving;
    _haRerouteInto = haStationary;
    _a2pToAppend.clear();
    _a2pToAppend[mcMeshProps().mesh().edge_handle(_haReroute)]
        = determineLastSuccessorPatch(_pCollapse, mcMeshProps().mesh().edge_handle(_haReroute));
}

bool EmbeddingCollapser::currentlyCollapsingPatch() const
{
    return _pCollapse.is_valid();
}

EmbeddingCollapser::RetCode EmbeddingCollapser::collapseArcEmbedding(const HEH& ha)
{
    DLOG(INFO) << "PRECOLLAPSE of halfarc " << ha << " which has "
               << mcMeshProps().ref<ARC_MESH_HALFEDGES>(mcMeshProps().mesh().edge_handle(ha)).size() << " halfedges";

    setVars(ha);

    // Mark vertices of collapse arc as touched to guide local collapsing/remeshing
    if (meshProps().isAllocated<TOUCHED>())
        for (HEH he : _collapseHes)
        {
            VH vFrom = meshProps().mesh().from_vertex_handle(he);
            meshProps().set<TOUCHED>(vFrom, true);
            for (VH vRing : meshProps().mesh().vertex_vertices(vFrom))
                meshProps().set<TOUCHED>(vRing, true);
        }

    return collapseEdgeByEdge();
}

EmbeddingCollapser::RetCode
EmbeddingCollapser::collapsePillowPatchEmbedding(const FH& p, const HEH& haMoving, const HEH& haStationary)
{
    DLOG(INFO) << "PRECOLLAPSE of patch " << p << " which has " << mcMeshProps().ref<PATCH_MESH_HALFFACES>(p).size()
               << " halffaces";

    setVars(p, haMoving, haStationary);

    // Mark vertices of collapse patch as touched to guide local collapsing/remeshing
    if (meshProps().isAllocated<TOUCHED>())
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(_pCollapse))
            for (VH v : meshProps().get_halfface_vertices(hf))
            {
                meshProps().set<TOUCHED>(v, true);
                for (VH v2 : meshProps().mesh().vertex_vertices(v))
                    meshProps().set<TOUCHED>(v2, true);
            }

    return collapseFaceByFace();
}

EmbeddingCollapser::RetCode
EmbeddingCollapser::collapsePillowBlockEmbedding(const CH& b, const HFH& hpMoving, const HFH& hpStationary)
{
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();

    // Mark vertices of collapse block as touched to guide local collapsing/remeshing
    if (meshProps().isAllocated<TOUCHED>())
        for (CH tet : mcMeshProps().get<BLOCK_MESH_TETS>(b))
            for (VH v : tetMesh.tet_vertices(tet))
            {
                meshProps().set<TOUCHED>(v, true);
                for (VH vRing : tetMesh.vertex_vertices(v))
                    meshProps().set<TOUCHED>(vRing, true);
            }

    FH p1 = mcMesh.face_handle(hpStationary);
    FH p2 = mcMesh.face_handle(hpMoving);

    CH bOpp = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hpMoving));
    applyTransitionToBlock(mcMeshProps().hpTransition<PATCH_TRANSITION>(hpMoving), b);
    assert(mcMeshProps().hpTransition<PATCH_TRANSITION>(hpMoving).isIdentity());

    auto newPatchHalffaces = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p1);
    bool flip1 = (hpStationary.idx() % 2) == (hpMoving.idx() % 2);
    if (flip1)
    {
        set<HFH> swapHalffaces;
        for (HFH hf : newPatchHalffaces)
            swapHalffaces.insert(tetMesh.opposite_halfface_handle(hf));
        newPatchHalffaces = std::move(swapHalffaces);
    }
    Transition transFromTo = mcMeshProps().hpTransition<PATCH_TRANSITION>(hpMoving);
    // Transfer tets, remove p2 mapping
    auto& tets = mcMeshProps().ref<BLOCK_MESH_TETS>(b);
    auto& tetsOpp = mcMeshProps().ref<BLOCK_MESH_TETS>(bOpp);
    for (CH tet : tets)
    {
        for (auto& kv : meshProps().ref<CHART>(tet))
            kv.second = transFromTo.apply(kv.second);
        meshProps().set<MC_BLOCK>(tet, bOpp);
        assert(!tetMesh.is_deleted(tet));
    }
    tetsOpp.insert(tets.begin(), tets.end());
    tets.clear();

    for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p2))
    {
        assert(meshProps().get<MC_PATCH>(tetMesh.face_handle(hf)) == p2);
        meshProps().reset<MC_PATCH>(tetMesh.face_handle(hf));
        meshProps().reset<IS_WALL>(tetMesh.face_handle(hf));
        meshProps().reset<TRANSITION>(tetMesh.face_handle(hf));
        if (meshProps().isAllocated<TOUCHED>())
        {
            for (VH v : meshProps().get_halfface_vertices(hf))
            {
                meshProps().set<TOUCHED>(v, true);
                for (VH v2 : tetMesh.vertex_vertices(v))
                    meshProps().set<TOUCHED>(v2, true);
            }
        }
    }
    mcMeshProps().set<PATCH_MESH_HALFFACES>(p2, newPatchHalffaces);
    mcMeshProps().setHpTransition<PATCH_TRANSITION>(
        hpStationary, transFromTo.invert().chain(mcMeshProps().hpTransition<PATCH_TRANSITION>(hpStationary)));

    return SUCCESS;
}

EmbeddingCollapser::RetCode EmbeddingCollapser::collapseCigarBlockEmbedding(const CH& b)
{
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();

    // Mark vertices of collapse block as touched to guide local collapsing/remeshing
    if (meshProps().isAllocated<TOUCHED>())
        for (CH tet : mcMeshProps().get<BLOCK_MESH_TETS>(b))
            for (VH v : tetMesh.tet_vertices(tet))
            {
                meshProps().set<TOUCHED>(v, true);
                for (VH vRing : tetMesh.vertex_vertices(v))
                    meshProps().set<TOUCHED>(vRing, true);
            }

    HFH hpCigar = *mcMesh.chf_iter(b);
    CH bOther = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hpCigar));

    auto tets = mcMeshProps().get<BLOCK_MESH_TETS>(b);
    mcMeshProps().ref<BLOCK_MESH_TETS>(b).clear();
    assert(bOther.is_valid());
    auto& tetsOther = mcMeshProps().ref<BLOCK_MESH_TETS>(bOther);
    tetsOther.insert(tets.begin(), tets.end());
    for (CH tet : tets)
    {
        meshProps().set<MC_BLOCK>(tet, bOther);
        assert(!tetMesh.is_deleted(tet));
    }

    Transition transFromTo = mcMeshProps().hpTransition<PATCH_TRANSITION>(hpCigar);
    for (CH tet : tets)
        for (auto& kv : meshProps().ref<CHART>(tet))
            kv.second = transFromTo.apply(kv.second);
    auto hfsCigar = mcMeshProps().get<PATCH_MESH_HALFFACES>(mcMesh.face_handle(hpCigar));
    mcMeshProps().ref<PATCH_MESH_HALFFACES>(mcMesh.face_handle(hpCigar)).clear();
    for (HFH hf : hfsCigar)
    {
        FH f = tetMesh.face_handle(hf);
#ifndef NDEBUG
        for (auto tet : tetMesh.face_cells(f))
        {
            auto block = meshProps().get<MC_BLOCK>(tet);
            if (block != bOther)
                LOG(ERROR) << "Non-patch face " << f << " incident to blocks " << block << " and " << bOther
                           << std::endl;
            assert(block == bOther);
        }
#endif
        meshProps().reset<TRANSITION>(f);
        meshProps().reset<MC_PATCH>(f);
        meshProps().reset<IS_WALL>(f);
        if (meshProps().isAllocated<TOUCHED>())
        {
            for (VH v : meshProps().get_halfface_vertices(hf))
            {
                meshProps().set<TOUCHED>(v, true);
                for (VH v2 : tetMesh.vertex_vertices(v))
                    meshProps().set<TOUCHED>(v2, true);
            }
        }
    }

    return SUCCESS;
}

void EmbeddingCollapser::nodeShift(const HEH& he)
{
    auto& tetMesh = meshProps().mesh();
    EH e0 = tetMesh.edge_handle(he);
    VH vFrom = tetMesh.from_vertex_handle(he);
    VH vTo = tetMesh.to_vertex_handle(he);

    VH n = meshProps().get<MC_NODE>(vFrom);

    meshProps().reset<MC_NODE>(vFrom);
    mcMeshProps().set<NODE_MESH_VERTEX>(n, vTo);

    if (!meshProps().get<MC_NODE>(vTo).is_valid())
        meshProps().set<MC_NODE>(vTo, n);

    EH a = meshProps().get<MC_ARC>(e0);
    meshProps().reset<MC_ARC>(e0);
    meshProps().reset<IS_ARC>(e0);
    auto& aHes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);
    if (aHes.front() == he)
        aHes.pop_front();
    else
        aHes.pop_back();
}

set<CH> EmbeddingCollapser::getCurrentLayerTets() const
{
    auto& tetMesh = meshProps().mesh();

    set<CH> layerTets;
    for (CH tetSeed : tetMesh.halfedge_cells(_heCurrent))
    {
        if (layerTets.count(tetSeed) != 0)
            continue;
        auto tets = getDomeTets(tetSeed);
        layerTets.insert(tets.begin(), tets.end());
    }
    return layerTets;
}

set<FH> EmbeddingCollapser::getCurrentLayerPatchFaces() const
{
    auto& tetMesh = meshProps().mesh();

    set<FH> layerFaces;
    for (FH fPatchSeed : tetMesh.halfedge_faces(_heCurrent))
    {
        if (!meshProps().isInPatch(fPatchSeed))
            continue;
        layerFaces.insert(fPatchSeed);
    }
    return layerFaces;
}

pair<HEH, vector<HFH>> EmbeddingCollapser::patchCorner(const HEH& heStart, const HFH& patchHf) const
{
    auto& tetMesh = meshProps().mesh();

    // Traverse triangle dome to find next arc in same patch
    pair<HEH, vector<HFH>> dome;
    HEH heCurr = heStart;
    HFH hfCurr = patchHf;
    do
    {
        dome.second.push_back(hfCurr);
        heCurr = tetMesh.prev_halfedge_in_halfface(heCurr, hfCurr);
        hfCurr = adjacentHfOnWall(hfCurr, heCurr);
        heCurr = tetMesh.opposite_halfedge_handle(heCurr);
    } while (!meshProps().isInArc(heCurr) && heCurr != heStart);
    dome.first = heCurr;
    return dome;
}

bool EmbeddingCollapser::reconnectNextArcToNodeThroughPatch()
{
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    HEH he0 = _heCurrent;
    EH e0 = tetMesh.edge_handle(he0);

    // find arc-a1 halfedge e1 incident on from_vertex(he) layer1-patch-connected to he0
    for (HFH hfPatch : tetMesh.halfedge_halffaces(he0))
    {
        FH fPatch = tetMesh.face_handle(hfPatch);

        // only search within patch
        if (!meshProps().isInPatch(fPatch))
            continue;
        // only search within layer
        if (_layerFaces.count(fPatch) == 0)
            continue;

        // Only incident to volume layer
        if (!containsMatching(tetMesh.face_cells(fPatch),
                              [&_layerTets = _layerTets](const CH& tet) { return _layerTets.count(tet) != 0; }))
            continue;

        auto sector = patchCorner(he0, hfPatch);
        // Skip if sector spans whole disk
        if (sector.first == he0)
        {
            for (HFH hf : sector.second)
                _layerFaces.erase(tetMesh.face_handle(hf));
            continue;
        }
        HEH heArc = tetMesh.opposite_halfedge_handle(sector.first); // Flip so heArc points to fan center vertex
        EH a = meshProps().get<MC_ARC>(tetMesh.edge_handle(heArc));
        vector<HFH>& triangleFan = sector.second;

        // Only reroute boundary arcs within boundary
        if ((tetMesh.is_boundary(heArc) && !tetMesh.is_boundary(fPatch)))
            continue;

        // Only reroute feature-patch-arcs within feature-patch
        bool inFeature = mcMeshProps().isAllocated<IS_FEATURE_F>()
                         && containsMatching(mcMesh.edge_faces(a),
                                             [this](const FH& p) { return mcMeshProps().get<IS_FEATURE_F>(p); });
        if (inFeature && !mcMeshProps().get<IS_FEATURE_F>(meshProps().get<MC_PATCH>(fPatch)))
            continue;

        // Special handling of last successor
        if (a == _aToAppend)
        {
            if (meshProps().isInArc(e0)) // was already appended
                continue;
            appendHalfedge(heArc);
            return true;
        }

        arcShift(heArc, triangleFan);

        return true;
    }
    return false;
}

bool EmbeddingCollapser::reconnectNextArcToNodeThroughBlock()
{
    auto& tetMesh = meshProps().mesh();
    auto& mcMesh = mcMeshProps().mesh();
    HEH he0 = _heCurrent;
    EH e0 = tetMesh.edge_handle(he0);
    VH v0 = tetMesh.from_vertex_handle(he0);

    set<CH> tetsVisited;
    for (CH tetSeed : tetMesh.halfedge_cells(he0))
    {
        if (tetsVisited.count(tetSeed) != 0 || _layerTets.count(tetSeed) == 0)
            continue;
        auto domeTets = getDomeTets(tetSeed);
        tetsVisited.insert(domeTets.begin(), domeTets.end());

        // Find incoming arc on boundary of this cell dome
        for (HEH heArc : tetMesh.incoming_halfedges(v0))
            if (heArc != he0 && !tetMesh.is_boundary(heArc) && meshProps().isInArc(heArc)
                && containsSomeOf(tetMesh.halfedge_cells(heArc), domeTets))
            {
                // Dont reroute feature-arc
                EH a = meshProps().get<MC_ARC>(tetMesh.edge_handle(heArc));
                if (mcMeshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps().get<IS_FEATURE_E>(a))
                    continue;

                // Dont reroute feature-patch-arc through block
                if (mcMeshProps().isAllocated<IS_FEATURE_F>()
                    && containsMatching(mcMesh.edge_faces(a),
                                        [this](const FH& p) { return mcMeshProps().get<IS_FEATURE_F>(p); }))
                    continue;

                // Special handling of last successor
                if (meshProps().get<MC_ARC>(tetMesh.edge_handle(heArc)) == _aToAppend)
                {
                    if (meshProps().isInArc(e0)) // Already appended
                        continue;
                    appendHalfedge(heArc);
                    return true;
                }

                VH from = tetMesh.to_vertex_handle(he0);
                VH to = tetMesh.from_vertex_handle(heArc);

                // Find any path between from and to
                list<HEH> path = findPathOnDomeHull(from, to, domeTets);
                assert(!path.empty());

                // Transform path into a triangle fan
                // Refinement of triangle fan is done in arcShift()
                vector<HFH> triangleFan;
                for (HEH he : path)
                    triangleFan.push_back(findMatching(
                        tetMesh.halfedge_halffaces(he),
                        [&, this](const HFH& hf) { return contains(meshProps().get_halfface_vertices(hf), v0); }));

                arcShift(heArc, triangleFan);

                return true;
            }
    }
    return false;
}

bool EmbeddingCollapser::isShiftableDome(const DomeElements& domeElems) const
{
    auto& tetMesh = meshProps().mesh();
    set<FH> ps;
    for (HFH hf : domeElems.floor)
    {
        FH p = meshProps().get<MC_PATCH>(tetMesh.face_handle(hf));
        ps.insert(p);
    }
    if (ps.size() > 1)
        return false;
    if (ps.size() == 1 && (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(*ps.begin())))
        return false;

    // CeilingBorderOrdered may be empty if dome is not manifold (i.e. invalid)
    if (domeElems.ceilingBorderOrdered.empty())
        return false;

    if (containsMatching(domeElems.floorInnerEdges,
                         [&, this](const EH& e) {
                             return (currentlyCollapsingPatch() || tetMesh.edge_handle(_heCurrent) != e)
                                    && meshProps().isInArc(e);
                         }))
        return false;

    return true;
}

bool EmbeddingCollapser::pullAwayNextPatch()
{
    auto& tetMesh = meshProps().mesh();
    VH vCenter = currentlyCollapsingPatch() ? _vPullUp : tetMesh.from_vertex_handle(_heCurrent);

    // No need to unravel single patch
    {
        set<FH> ps;
        for (FH f : tetMesh.vertex_faces(vCenter))
        {
            FH p = meshProps().get<MC_PATCH>(f);
            if (p.is_valid())
                ps.insert(p);
        }
        if (ps.size() <= 1)
            return false;
    }

    for (CH tetSeed : tetMesh.vertex_cells(vCenter))
    {
        if (_layerTets.count(tetSeed) == 0)
            continue;

        auto domeTets = getDomeTets(tetSeed);
        auto domeElems = getDomeConnectivity(domeTets);

        if (!isShiftableDome(domeElems))
        {
            for (CH tet : domeTets)
                _layerTets.erase(tet);
            continue;
        }

        refineAndUpdateDome(domeTets, domeElems);

        patchDomeShift(domeTets, domeElems.floor, domeElems.ceiling);

        return true;
    }
    return false;
}

set<CH> EmbeddingCollapser::getDomeTets(const CH& tetSeed) const
{
    auto& tetMesh = meshProps().mesh();
    VH vCenter = currentlyCollapsingPatch() ? _vPullUp : tetMesh.from_vertex_handle(_heCurrent);

    // all cells incident to vcenter volume-connected to tetSeed, stopping at walls
    set<CH> domeTets({tetSeed});
    list<CH> tetQ({tetSeed});
    while (!tetQ.empty())
    {
        CH tet = tetQ.front();
        tetQ.pop_front();
        for (HFH hf : tetMesh.cell_halffaces(tet))
        {
            if (meshProps().isInPatch(hf) || !contains(tetMesh.halfface_vertices(hf), vCenter))
                continue;
            CH tetNext = tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf));
            assert(tetNext.is_valid());
            if (domeTets.count(tetNext) == 0)
            {
                domeTets.insert(tetNext);
                tetQ.push_back(tetNext);
            }
        }
    }
    return domeTets;
}

void EmbeddingCollapser::arcShift(const HEH& heArc, vector<HFH>& triangleFan)
{
    auto& tetMesh = meshProps().mesh();
    HEH he0 = _heCurrent;
    VH vCollapseFrom = tetMesh.from_vertex_handle(_heCurrent);

    EH aReroute = meshProps().get<MC_ARC>(tetMesh.edge_handle(heArc));
    auto& aHes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(aReroute);
    bool revert = heArc != aHes.back();

    // Remove all fan hfs from layerFaces
    for (HFH hf : triangleFan)
        _layerFaces.erase(tetMesh.face_handle(hf));

    refineAndUpdateTriangleFan(heArc, triangleFan);

    auto stackOfHfs = triangleFan;
    while (!stackOfHfs.empty())
    {
        HFH hf = stackOfHfs.back();
        stackOfHfs.pop_back();

        // Order so that edge not touching vCollapseFrom is always first
        list<HEH> hes;
        for (HEH he : tetMesh.halfface_halfedges(hf))
            hes.push_back(he);
        while (tetMesh.from_vertex_handle(hes.front()) == vCollapseFrom
               || tetMesh.to_vertex_handle(hes.front()) == vCollapseFrom)
        {
            hes.push_back(hes.front());
            hes.pop_front();
        }

        // Toggle arc status for each halfedge of current hf
        for (HEH he : hes)
        {
            EH e = tetMesh.edge_handle(he);

            if (tetMesh.edge_handle(he0) == e)
                continue;

            if (meshProps().isInArc(e))
            {
                meshProps().reset<IS_ARC>(e);
                meshProps().reset<MC_ARC>(e);
                if (revert)
                    aHes.erase(
                        std::next(std::find(aHes.rbegin(), aHes.rend(), tetMesh.opposite_halfedge_handle(he))).base());
                else
                    aHes.erase(std::find(aHes.begin(), aHes.end(), he));
            }
            else
            {
                meshProps().set<IS_ARC>(e, true);
                meshProps().set<MC_ARC>(e, aReroute);
                if (revert)
                    aHes.push_front(he);
                else
                    aHes.push_back(tetMesh.opposite_halfedge_handle(he));
            }
        }
        // ... and reconnect patches to arc
        for (HEH he : hes)
        {
            EH e = tetMesh.edge_handle(he);
            if (tetMesh.edge_handle(he0) == e)
                continue;
            if (!meshProps().isInArc(e))
            {
                // Reconnect all incident patch faces in a cycle around he starting from hf
                if (tetMesh.is_boundary(hf))
                    reconnectPatchesToArc(
                        tetMesh.opposite_halfedge_handle(he), aReroute, tetMesh.opposite_halfface_handle(hf));
                else
                    reconnectPatchesToArc(he, aReroute, hf);
            }
        }
    }
}

void EmbeddingCollapser::reconnectPatchesToArc(const HEH& he, const EH& a, const HFH& hfFlip)
{
    auto& tetMesh = meshProps().mesh();
    VH vCenter = currentlyCollapsingPatch() ? _vPullUp : tetMesh.from_vertex_handle(_heCurrent);
    FH fFlip = tetMesh.face_handle(hfFlip);

    FH fJoin;
    // choose a patch to take up previous patch space (if one exists)
    FH pFirst = meshProps().get<MC_PATCH>(fFlip);
    FH pToAppend;
    if (_a2pToAppend.count(a) != 0)
        pToAppend = _a2pToAppend.at(a);
    else if (pFirst.is_valid())
        pToAppend = determineLastSuccessorPatch(pFirst, a);

    if (pToAppend.is_valid())
    {
        fJoin = findMatching(tetMesh.halfedge_faces(he),
                             [&](const FH& f) { return f != fFlip && meshProps().get<MC_PATCH>(f) == pToAppend; });
        if (pToAppend != pFirst && !fJoin.is_valid())
            throw std::logic_error("Cant append pToAppend");
    }

    // Get the ring of block domes around he in an ordered fashin (via propeller of halffaces around he)
    auto hfsPair = cyclicOrderHalfFaces(he, hfFlip, fJoin);
    auto hfsFront = hfsPair.first;
    auto hfsBack = hfsPair.second;

    {
        // Last patch halfface can occupy the space of hfFlip
        HFH hfLastPatch = hfsFront.back().front();
        FH fLastPatch = tetMesh.face_handle(hfLastPatch);
        assert(hfsBack.empty() || fLastPatch == tetMesh.face_handle(hfsBack.back().front()));
        if (fJoin.is_valid())
            assert(fLastPatch == fJoin);
        hfsFront.pop_back();
        if (!hfsBack.empty())
            hfsBack.pop_back();

        if (_pCollapse.is_valid())
            assert(pFirst == _pCollapse);

        Transition transLast = meshProps().hfTransition<TRANSITION>(hfLastPatch);

        // Cache the appended patch, if none has been determined before
        FH pLast = meshProps().get<MC_PATCH>(fLastPatch);
        assert(pLast.is_valid());
        if (_a2pToAppend.count(a) == 0)
            _a2pToAppend[a] = pLast;

        if (pLast == _pCollapse)
            throw std::logic_error("Trying to append pCollapse into pCollapse");

        // Transfer properties of hfFlip between the two patches
        if (pFirst != pLast)
        {
            if (pFirst.is_valid())
            {
                auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(pFirst);
                bool revert = pHfs.count(hfFlip) == 0;
                pHfs.erase(revert ? tetMesh.opposite_halfface_handle(hfFlip) : hfFlip);
            }
            auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(pLast);
            bool revert = pHfs.count(hfLastPatch) == 0;
            pHfs.insert(revert ? tetMesh.opposite_halfface_handle(hfFlip) : hfFlip);
            meshProps().set<IS_WALL>(fFlip, true);
            meshProps().set<MC_PATCH>(fFlip, pLast);
            meshProps().setTransition<TRANSITION>(hfFlip, transLast);
        }
    }

    // Now shift the block dome walls back onto the shifted arc in peel by peel (first cw then ccw)
    for (bool front : {true, false})
    {
        auto& hfs = front ? hfsFront : hfsBack;
        HEH hePivot = front ? tetMesh.opposite_halfedge_handle(he) : he;
        for (auto it = hfs.begin(); it != hfs.end(); it++)
        {
            auto& hfsPropeller = *it;
            HFH hfPatch = hfsPropeller.front();
            FH p = meshProps().get<MC_PATCH>(tetMesh.face_handle(hfPatch));
            if (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p))
            {
                assert(false);
                throw std::logic_error("Error, have to move feature patch");
            }

            auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);
            bool revert = pHfs.count(hfPatch) == 0;

            CH bFrom = meshProps().get<MC_BLOCK>(tetMesh.incident_cell(hfPatch));
            CH bTo = meshProps().get<MC_BLOCK>(tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hfPatch)));

            Transition transFromTo = meshProps().hfTransition<TRANSITION>(hfPatch);

            refineAndUpdateTetFan(hePivot, hfFlip, hfPatch, hfsPropeller);

            meshProps().replaceByChildren(_layerTets);
            list<CH> tetsToErase;
            for (CH tet : _layerTets)
                if (!contains(tetMesh.tet_vertices(tet), vCenter))
                    tetsToErase.push_back(tet);
            for (CH tet : tetsToErase)
                _layerTets.erase(tet);

            meshProps().replaceByChildren(_layerFaces);
            list<FH> fsToErase;
            for (FH f : _layerFaces)
                if (!contains(tetMesh.face_vertices(f), vCenter))
                    fsToErase.push_back(f);
            for (FH f : fsToErase)
                _layerFaces.erase(f);

            // Now shift patch across tets and transfer tets between blocks
            for (HFH hfPropeller : hfsPropeller)
            {
                assert(!tetMesh.is_deleted(hfPropeller));
                CH tet = tetMesh.incident_cell(hfPropeller);
                for (HFH hf : tetMesh.cell_halffaces(tet))
                {
                    FH f = tetMesh.face_handle(hf);
                    if (f == fFlip)
                        continue;
                    if (meshProps().isInPatch(f))
                    {
                        assert(p == meshProps().get<MC_PATCH>(f));
                        pHfs.erase(revert ? tetMesh.opposite_halfface_handle(hf) : hf);
                        meshProps().reset<IS_WALL>(f);
                        meshProps().reset<MC_PATCH>(f);
                        meshProps().reset<TRANSITION>(f);
                    }
                    else
                    {
                        pHfs.insert(revert ? hf : tetMesh.opposite_halfface_handle(hf));
                        meshProps().set<IS_WALL>(f, true);
                        meshProps().set<MC_PATCH>(f, p);
                        meshProps().setTransition<TRANSITION>(hf, transFromTo.invert());
                    }
                }
                assert(meshProps().get<MC_BLOCK>(tet) == bFrom);
                assert(!tetMesh.is_deleted(tet));
                meshProps().set<MC_BLOCK>(tet, bTo);
                mcMeshProps().ref<BLOCK_MESH_TETS>(bFrom).erase(tet);
                mcMeshProps().ref<BLOCK_MESH_TETS>(bTo).insert(tet);
                _layerTets.erase(tet);

                for (auto& kv : meshProps().ref<CHART>(tet))
                    kv.second = transFromTo.apply(kv.second);
            }

            auto itNext = it;
            itNext++;
            if (itNext != hfs.end())
                for (HFH hf : hfsPropeller)
                    itNext->push_back(hf);
        }
    }
}

pairTT<list<list<HFH>>> EmbeddingCollapser::cyclicOrderHalfFaces(const HEH& he, const HFH& hf, const FH& fJoin) const
{
    auto& tetMesh = meshProps().mesh();

    assert(!tetMesh.is_boundary(hf));

    FH f0 = tetMesh.face_handle(hf);
    list<FH> fs({f0});

    HFH hfFront = hf;
    HFH hfBack = tetMesh.opposite_halfface_handle(hfFront);
    HEH heFront = he;
    HEH heBack = tetMesh.opposite_halfedge_handle(heFront);

    assert(contains(tetMesh.halfface_halfedges(hfFront), heFront));
    assert(contains(tetMesh.halfface_halfedges(hfBack), heBack));

    pairTT<list<list<HFH>>> hfLists;
    bool front = false;
    do
    {
        if (tetMesh.is_boundary(hfBack) || (fJoin.is_valid() && tetMesh.face_handle(hfBack) == fJoin))
            front = true;
        else if (tetMesh.is_boundary(hfFront) || (fJoin.is_valid() && tetMesh.face_handle(hfFront) == fJoin))
            front = false;
        else
            front = !front;

        assert(!tetMesh.is_boundary(front ? hfFront : hfBack));
        (front ? hfLists.first : hfLists.second).emplace_back();
        do
        {
            (front ? hfFront : hfBack) = tetMesh.opposite_halfface_handle(
                tetMesh.adjacent_halfface_in_cell((front ? hfFront : hfBack), (front ? heFront : heBack)));
            (front ? hfLists.first : hfLists.second)
                .back()
                .push_front(front ? tetMesh.opposite_halfface_handle(hfFront)
                                  : tetMesh.opposite_halfface_handle(hfBack));
        } while (!meshProps().isInPatch(tetMesh.face_handle(front ? hfFront : hfBack)));
    } while ((!tetMesh.is_boundary(hfBack) || !tetMesh.is_boundary(hfFront))
             && hfFront != tetMesh.opposite_halfface_handle(hfBack));

    return hfLists;
}

void EmbeddingCollapser::appendHalfedge(const HEH& he)
{
    auto& tetMesh = meshProps().mesh();
    HEH he0 = _heCurrent;
    EH e0 = tetMesh.edge_handle(he0);
    assert(!meshProps().isInArc(e0));

    EH a = meshProps().get<MC_ARC>(tetMesh.edge_handle(he));
    auto& aHes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);
    bool reverse = std::find(aHes.begin(), aHes.end(), he) == aHes.end();
    if (reverse)
        aHes.push_front(tetMesh.opposite_halfedge_handle(he0));
    else
        aHes.push_back(he0);
    meshProps().set<IS_ARC>(e0, true);
    meshProps().set<MC_ARC>(e0, a);
}

EmbeddingCollapser::DomeElements EmbeddingCollapser::getDomeConnectivity(const set<CH>& domeTets) const
{
    auto& tetMesh = meshProps().mesh();
    VH vCenter = currentlyCollapsingPatch() ? _vPullUp : tetMesh.from_vertex_handle(_heCurrent);

    DomeElements elems;

    for (CH tet : domeTets)
        for (HFH hf : tetMesh.cell_halffaces(tet))
            if (domeTets.count(tetMesh.incident_cell(tetMesh.opposite_halfface_handle(hf))) == 0)
            {
                if (contains(meshProps().get_halfface_vertices(hf), vCenter))
                    elems.floor.insert(hf);
                else
                    elems.ceiling.insert(hf);
            }
    for (HFH hf : elems.ceiling)
        for (HEH he : tetMesh.halfface_halfedges(hf))
        {
            bool hasHf = false;
            for (HFH hf2 : tetMesh.halfedge_halffaces(tetMesh.opposite_halfedge_handle(he)))
                if (elems.ceiling.count(hf2) != 0)
                {
                    hasHf = true;
                    break;
                }
            if (!hasHf)
                elems.ceilingBorder.insert(he);
            else
                elems.ceilingInnerEdges.insert(tetMesh.edge_handle(he));
        }
    for (HFH hf : elems.floor)
        for (HEH he : tetMesh.halfface_halfedges(hf))
            if (elems.ceilingBorder.count(he) == 0
                && elems.ceilingBorder.count(tetMesh.opposite_halfedge_handle(he)) == 0)
                elems.floorInnerEdges.insert(tetMesh.edge_handle(he));

    bool nonManifold = false;
    {
        map<VH, int> incidenceIn;
        map<VH, int> incidenceOut;
        for (HEH he : elems.ceilingBorder)
        {
            incidenceIn[tetMesh.to_vertex_handle(he)]++;
            incidenceOut[tetMesh.from_vertex_handle(he)]++;
        }
        for (auto kv : incidenceIn)
            if (kv.second != 1)
            {
                nonManifold = true;
                break;
            }
        if (!nonManifold)
            for (auto kv : incidenceOut)
                if (kv.second != 1)
                {
                    nonManifold = true;
                    break;
                }
    }
    for (HEH he : elems.ceilingBorder)
        elems.ceilingBorderVertices.insert(tetMesh.to_vertex_handle(he));

    if (!elems.ceilingBorder.empty() && !nonManifold)
    {
        HEH startHe = *elems.ceilingBorder.begin();
        HEH currHe = startHe;
        do
        {
            elems.ceilingBorderOrdered.push_back(currHe);
            VH vTo = tetMesh.to_vertex_handle(currHe);
            for (HEH he : elems.ceilingBorder)
                if (tetMesh.from_vertex_handle(he) == vTo)
                {
                    currHe = he;
                    break;
                }
        } while (currHe != startHe);
    }
    for (EH e : elems.ceilingInnerEdges)
        for (VH v : tetMesh.edge_vertices(e))
            if (elems.ceilingBorderVertices.count(v) == 0)
                elems.ceilingInnerVertices.insert(v);

    return elems;
}

list<HEH> EmbeddingCollapser::findPathOnDomeHull(const VH& from, const VH& to, set<CH>& domeTets)
{
    auto& tetMesh = meshProps().mesh();

    // Get dome connectivity
    auto domeElems = getDomeConnectivity(domeTets);

    set<VH> border1, border2;
    bool first = true;
    for (HEH he : domeElems.ceilingBorderOrdered)
    {
        for (VH v : tetMesh.halfedge_vertices(he))
            if (v != from && v != to)
                (first ? border1 : border2).insert(v);
        VH v = tetMesh.to_vertex_handle(he);
        if (v == from || v == to)
            first = !first;
    }

    // Refine dome boundary
    if (domeElems.ceilingInnerEdges.empty())
    {
        assert(domeElems.ceiling.size() == 1);
        VH vNew = splitFace(tetMesh.face_handle(*domeElems.ceiling.begin()), {Q(1, 3), Q(1, 3), Q(1, 3)});
        domeElems.ceilingInnerVertices.insert(vNew);
        for (VH v : domeElems.ceilingBorderVertices)
            domeElems.ceilingInnerEdges.insert(tetMesh.edge_handle(tetMesh.find_halfedge(v, vNew)));
        meshProps().replaceByChildren(domeElems.ceiling);
        meshProps().replaceByChildren(_layerTets);
        meshProps().replaceByChildren(domeTets);
    }
    else
    {
        bool found = false;
        bool opp = false;
        HEH he = tetMesh.find_halfedge(from, to);
        if (he.is_valid())
        {
            if (domeElems.ceilingBorder.count(he) != 0)
            {
                found = true;
                opp = false;
            }
            else if (domeElems.ceilingBorder.count(tetMesh.opposite_halfedge_handle(he)) != 0)
            {
                found = true;
                opp = true;
            }
        }
        if (found)
        {
            if (opp)
                he = tetMesh.opposite_halfedge_handle(he);
            HFH hfIncident = findSomeOf(tetMesh.halfedge_halffaces(he), domeElems.ceiling);
            assert(hfIncident.is_valid());
            int nBoundaryVs = 0;
            list<VH> hfVs;
            for (VH v : meshProps().get_halfface_vertices(hfIncident))
            {
                hfVs.push_back(v);
                if (domeElems.ceilingBorderVertices.count(v) != 0)
                    nBoundaryVs++;
            }
            if (nBoundaryVs == 3)
            {
                VH vNew = splitFace(tetMesh.face_handle(hfIncident), {Q(1, 3), Q(1, 3), Q(1, 3)});
                domeElems.ceilingInnerVertices.insert(vNew);
                for (VH v : hfVs)
                    domeElems.ceilingInnerEdges.insert(tetMesh.edge_handle(tetMesh.find_halfedge(v, vNew)));
                meshProps().replaceByChildren(domeElems.ceiling);
                meshProps().replaceByChildren(_layerTets);
                meshProps().replaceByChildren(domeTets);
            }
        }
        else
        {
            set<EH> splitEs;
            for (EH e : domeElems.ceilingInnerEdges)
            {
                auto vs = tetMesh.edge_vertices(e);
                if (domeElems.ceilingBorderOrdered.empty())
                {
                    if (domeElems.ceilingBorderVertices.count(vs[0]) != 0
                        && domeElems.ceilingBorderVertices.count(vs[1]) != 0)
                        splitEs.insert(e);
                }
                else if ((border1.count(vs[0]) != 0 && border2.count(vs[1]) != 0)
                         || (border1.count(vs[1]) != 0 && border2.count(vs[0]) != 0))
                    splitEs.insert(e);
            }
            for (EH e : splitEs)
            {
                auto eVs = tetMesh.edge_vertices(e);
                list<VH> vsOpp;
                for (HFH hf : tetMesh.edge_halffaces(e))
                    if (domeElems.ceiling.count(hf) != 0)
                        vsOpp.push_back(findNoneOf(tetMesh.halfface_vertices(hf), tetMesh.edge_vertices(e)));
                VH vNew = splitHalfEdge(tetMesh.halfedge_handle(e, 0), *tetMesh.ec_iter(e), Q(0.5));

                if (domeElems.ceilingInnerEdges.count(e) != 0)
                {
                    domeElems.ceilingInnerEdges.erase(e);
                    domeElems.ceilingInnerEdges.insert(tetMesh.edge_handle(tetMesh.find_halfedge(eVs[0], vNew)));
                    domeElems.ceilingInnerEdges.insert(tetMesh.edge_handle(tetMesh.find_halfedge(eVs[1], vNew)));
                    domeElems.ceilingInnerVertices.insert(vNew);
                }

                for (VH v : vsOpp)
                    domeElems.ceilingInnerEdges.insert(tetMesh.edge_handle(tetMesh.find_halfedge(v, vNew)));

                meshProps().replaceByChildren(domeElems.ceiling);
                meshProps().replaceByChildren(_layerTets);
                meshProps().replaceByChildren(domeTets);
            }
        }
    }

    return findPath(from, to, domeElems.ceilingInnerEdges, domeElems.ceilingInnerVertices);
}

list<HEH>
EmbeddingCollapser::findPath(const VH& from, const VH& to, const set<EH>& esAllowed, const set<VH>& vsAllowed) const
{
    auto& tetMesh = meshProps().mesh();

    assert(!esAllowed.empty());

    set<VH> vVisited({from});
    list<list<HEH>> heQ;
    for (HEH he : tetMesh.outgoing_halfedges(from))
        if (esAllowed.count(tetMesh.edge_handle(he)) != 0)
        {
            VH vNext = tetMesh.to_vertex_handle(he);
            if (tetMesh.to_vertex_handle(he) == to)
                return {he};
            if (vsAllowed.count(vNext) != 0)
            {
                vVisited.insert(vNext);
                heQ.push_back({he});
            }
        }

    while (!heQ.empty())
    {
        list<HEH> hes = heQ.front();
        heQ.pop_front();

        VH vTo = tetMesh.to_vertex_handle(hes.back());
        for (HEH heNext : tetMesh.outgoing_halfedges(vTo))
        {
            VH vNext = tetMesh.to_vertex_handle(heNext);
            if (esAllowed.count(tetMesh.edge_handle(heNext)) != 0 && vVisited.count(vNext) == 0)
            {
                auto hesNext = hes;
                hesNext.push_back(heNext);
                if (vNext == to)
                    return hesNext;
                if (vsAllowed.count(vNext) != 0)
                {
                    heQ.push_back(hesNext);
                    vVisited.insert(vNext);
                }
            }
        }
    }

    assert(false);
    return {};
}

FH EmbeddingCollapser::determineLastSuccessorPatch(const FH& pFirst, const EH& aShift) const
{
    auto& mcMesh = mcMeshProps().mesh();
    if (!pFirst.is_valid())
        return FH();
    if (mcMesh.is_boundary(pFirst))
        return findMatching(mcMesh.edge_faces(aShift),
                            [&](const FH& p) { return p != pFirst && mcMesh.is_boundary(p); });
    else if (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(pFirst))
        return findMatching(mcMesh.edge_faces(aShift),
                            [&, this](const FH& p) { return p != pFirst && mcMeshProps().get<IS_FEATURE_F>(p); });
    else
    {
        auto b2trans = determineTransitionsAroundArc(aShift, *mcMesh.ec_iter(aShift), Transition());
        HFH hpNonBoundary = mcMesh.halfface_handle(pFirst, 0);
        if (mcMesh.is_boundary(hpNonBoundary))
            hpNonBoundary = mcMesh.opposite_halfface_handle(hpNonBoundary);
        UVWDir hpFirstNormal
            = b2trans.at(mcMesh.incident_cell(hpNonBoundary)).invert().rotate(halfpatchNormalDir(hpNonBoundary));
        for (FH p : mcMesh.edge_faces(aShift))
        {
            if (p == pFirst || mcMesh.is_boundary(p) != mcMesh.is_boundary(pFirst)
                || (mcMeshProps().isAllocated<IS_FEATURE_F>()
                    && !mcMeshProps().get<IS_FEATURE_F>(pFirst) != !mcMeshProps().get<IS_FEATURE_F>(p)))
                continue;
            HFH hp = mcMesh.halfface_handle(p, 0);
            if (mcMesh.is_boundary(hp))
                hp = mcMesh.opposite_halfface_handle(hp);
            UVWDir normalHp = b2trans.at(mcMesh.incident_cell(hp)).invert().rotate(halfpatchNormalDir(hp));
            if (((normalHp | -normalHp) & hpFirstNormal) != UVWDir::NONE)
                return p;
        }
        return findMatching(mcMesh.edge_faces(aShift),
                            [&, this](const FH& p)
                            {
                                return p != pFirst && mcMesh.is_boundary(p) == mcMesh.is_boundary(pFirst)
                                       && (!mcMeshProps().isAllocated<IS_FEATURE_F>()
                                           || !mcMeshProps().get<IS_FEATURE_F>(pFirst)
                                                  == !mcMeshProps().get<IS_FEATURE_F>(p));
                            });
    }
    return FH();
}

EH EmbeddingCollapser::determineLastSuccessorArc() const
{
    auto& mcMesh = mcMeshProps().mesh();
    if (mcMeshProps().get<IS_SINGULAR>(mcMesh.edge_handle(_haCollapse)))
        return findMatching(mcMesh.vertex_edges(mcMesh.from_vertex_handle(_haCollapse)),
                            [&, this](const EH& a)
                            { return a != mcMesh.edge_handle(_haCollapse) && mcMeshProps().get<IS_SINGULAR>(a); });
    else if (mcMeshProps().isAllocated<IS_FEATURE_E>()
             && mcMeshProps().get<IS_FEATURE_E>(mcMesh.edge_handle(_haCollapse)))
        return findMatching(mcMesh.vertex_edges(mcMesh.from_vertex_handle(_haCollapse)),
                            [&, this](const EH& a)
                            { return a != mcMesh.edge_handle(_haCollapse) && mcMeshProps().get<IS_FEATURE_E>(a); });
    else
    {
        bool collapseInFeature = mcMeshProps().isAllocated<IS_FEATURE_F>()
                                 && containsMatching(mcMesh.halfedge_faces(_haCollapse),
                                                     [&](const FH& p) { return mcMeshProps().get<IS_FEATURE_F>(p); });

        CH bRef = *mcMesh.hec_iter(_haCollapse);
        // Choose a well aligned boundary arc
        auto b2trans = determineTransitionsAroundNode(mcMesh.from_vertex_handle(_haCollapse), bRef, Transition());
        UVWDir dirHa = halfarcDirInBlock(_haCollapse, bRef);

        for (bool aligned : {true, false})
            for (HEH ha : mcMesh.incoming_halfedges(mcMesh.from_vertex_handle(_haCollapse)))
                if (ha != mcMesh.opposite_halfedge_handle(_haCollapse)
                    && (mcMesh.is_boundary(_haCollapse) == mcMesh.is_boundary(ha))
                    && (!mcMeshProps().isAllocated<IS_FEATURE_E>()
                        || (!mcMeshProps().get<IS_FEATURE_E>(mcMesh.edge_handle(_haCollapse))
                            == !mcMeshProps().get<IS_FEATURE_E>(mcMesh.edge_handle(ha))))
                    && (!aligned
                        || dirHa
                               == b2trans.at(*mcMesh.hec_iter(ha))
                                      .invert()
                                      .rotate(halfarcDirInBlock(ha, *mcMesh.hec_iter(ha))))
                    && collapseInFeature
                           == (mcMeshProps().isAllocated<IS_FEATURE_F>()
                               && containsMatching(mcMesh.halfedge_faces(ha),
                                                   [&](const FH& p) { return mcMeshProps().get<IS_FEATURE_F>(p); })))
                    return mcMesh.edge_handle(ha);
    }
    return EH();
}

EmbeddingCollapser::RetCode EmbeddingCollapser::collapseEdgeByEdge()
{
    while (!_collapseHes.empty())
    {
        bool done = false;
        _heCurrent = _collapseHes.front();
        _collapseHes.pop_front();
        nodeShift(_heCurrent);

        while (!done)
        {
            done = true;

            _layerTets = getCurrentLayerTets();
            _layerFaces = getCurrentLayerPatchFaces();

            while (!_layerFaces.empty())
                if (reconnectNextArcToNodeThroughPatch())
                    done = false;
                else
                    break;
            while (!_layerTets.empty())
                if (reconnectNextArcToNodeThroughBlock())
                    done = false;
                else
                    break;
            while (!_layerTets.empty())
                if (pullAwayNextPatch())
                    done = false;
                else
                    break;
        }
        _heCurrent = {};
    }
    return SUCCESS;
}

EmbeddingCollapser::RetCode EmbeddingCollapser::collapseFaceByFace()
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    EH aReroute = mcMesh.edge_handle(_haReroute);
    HFH hpCollapse = mcMesh.halfface_handle(_pCollapse, 0);
    if ((_haReroute.idx() % 2) != 0)
    {
        _haReroute = mcMesh.opposite_halfedge_handle(_haReroute);
        _haRerouteInto = mcMesh.opposite_halfedge_handle(_haRerouteInto);
        hpCollapse = mcMesh.opposite_halfface_handle(hpCollapse);
    }

    set<HFH> hfsCollapse = mcMeshProps().hpHalffaces(hpCollapse);

    auto& aHes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(aReroute);

    set<HEH> hesCurve;
    set<HEH> hesLimit;
    set<HFH> hfsBoundary;
    for (HEH he : aHes)
    {
        hesCurve.insert(he);
        for (HFH hf : tetMesh.halfedge_halffaces(he))
            if (hfsCollapse.count(hf) != 0)
                hfsBoundary.insert(hf);
    }

    // TODO Bad fix/workaround for MCMesh connectivity breaking for some extreme selfadjacency
    if (hfsBoundary.empty())
    {
        hpCollapse = mcMesh.opposite_halfface_handle(hpCollapse);
        hfsCollapse = mcMeshProps().hpHalffaces(hpCollapse);
        for (HEH he : aHes)
            for (HFH hf : tetMesh.halfedge_halffaces(he))
                if (hfsCollapse.count(hf) != 0)
                    hfsBoundary.insert(hf);
    }

    for (HEH he : mcMeshProps().haHalfedges(_haRerouteInto))
        hesLimit.insert(tetMesh.opposite_halfedge_handle(he));

    bool skippedFull = true;
    list<HFH> hfs(hfsBoundary.begin(), hfsBoundary.end());
    if (hfs.empty())
    {
        LOG(ERROR) << "Patch with no boundary";
        return REROUTE_ERROR;
    }

    // Zipper algorithm
    for (auto it = hfs.begin(); it != hfs.end();)
    {

        HFH hf = *it;
        vector<VH> vs;
        list<HEH> hes;
        for (HEH he : tetMesh.halfface_halfedges(hf))
        {
            hes.push_back(he);
            vs.push_back(tetMesh.to_vertex_handle(he));
        }

        if (shiftMakesArcNonManifold(hesCurve, hes))
        {
            ++it;
            if (it == hfs.end())
            {
                if (skippedFull)
                {
                    LOG(ERROR) << "Patch to reroute is not genus 0 with boundary";
                    return REROUTE_ERROR;
                }
                skippedFull = true;
                it = hfs.begin();
            }
            continue;
        }
        skippedFull = false;

        it = hfs.erase(it);
        hfsBoundary.erase(hf);
        hfsCollapse.erase(hf);

        // Handle similar as arcShift
        // Order halfedges so that first halfedge is on path and last is not
        while (hesCurve.count(hes.front()) == 0 || hesCurve.count(hes.back()) != 0)
        {
            hes.push_back(hes.front());
            hes.pop_front();
        }

        if (hesCurve.count(*(++hes.begin())) != 0)
        {
            _heCurrent = *(++hes.begin());
            _vPullUp = tetMesh.to_vertex_handle(hes.front());
            _layerTets = getCurrentLayerTets();
        }
        else
        {
            _heCurrent = HEH();
            _vPullUp = VH();
        }

        auto itHes = std::find(aHes.begin(), aHes.end(), hes.front());

        int skip = 0;
        // Toggle arc status for each halfedge of current hf
        for (HEH he : hes)
        {
            EH e = tetMesh.edge_handle(he);
            if (hesCurve.count(he) != 0)
            {
                meshProps().reset<IS_ARC>(e);
                meshProps().reset<MC_ARC>(e);
                itHes = aHes.erase(itHes);
                hesCurve.erase(he);
            }
            else
            {
                HEH heOpp = tetMesh.opposite_halfedge_handle(he);
                itHes = aHes.insert(itHes, heOpp);
                hesCurve.insert(heOpp);
                HFH hfNewBoundary = findMatching(
                    tetMesh.halfedge_halffaces(heOpp),
                    [&](const HFH& hf2) { return hfsCollapse.count(hf2) != 0 && hfsBoundary.count(hf2) == 0; });
                if (hfNewBoundary.is_valid())
                {
                    skip++;
                    hfsBoundary.insert(hfNewBoundary);
                    it = hfs.insert(it, hfNewBoundary);
                }
                if (hesLimit.count(heOpp) == 0)
                {
                    meshProps().set<IS_ARC>(e, true);
                    meshProps().set<MC_ARC>(e, aReroute);
                }
            }
        }
        for (int i = 0; i < skip; i++)
            it++;

        // ... and reconnect patches to arc
        // if (!_vPullUp.is_valid())
        {
            // Reconnect all incident patch faces in a cycle around he starting from hf
            if (tetMesh.is_boundary(hf))
                reconnectPatchesToArc(
                    tetMesh.opposite_halfedge_handle(hes.front()), aReroute, tetMesh.opposite_halfface_handle(hf));
            else
                reconnectPatchesToArc(hes.front(), aReroute, hf);
        }

        // ... and remove patches from possible overlaps at vertex
        if (_vPullUp.is_valid())
        {
            bool done = false;
            while (!done)
            {
                done = true;
                while (!_layerTets.empty())
                    if (pullAwayNextPatch())
                        done = false;
                    else
                        break;
                if (!done)
                    _layerTets = getCurrentLayerTets();
            }
        }

        if (it == hfs.end())
            it = hfs.begin();
    }
    return SUCCESS;
}

bool EmbeddingCollapser::shiftMakesArcNonManifold(const set<HEH>& hesCurve, const list<HEH>& hfHes) const
{
    auto& tetMesh = meshProps().mesh();
    // If 2 edges are not in curve yet, check if enclosed vertex is in curve
    list<HEH> nonCurveHes;
    for (HEH he : hfHes)
        if (hesCurve.count(he) == 0)
            nonCurveHes.push_back(he);
    assert(!nonCurveHes.empty());
    assert(nonCurveHes.size() != 3);
    if (nonCurveHes.size() == 1)
        return false;

    auto vsA = tetMesh.halfedge_vertices(nonCurveHes.front());
    auto vsB = tetMesh.halfedge_vertices(nonCurveHes.back());
    VH vShared = vsA[1] == vsB[0] ? vsA[1] : vsA[0];
    return containsSomeOf(tetMesh.outgoing_halfedges(vShared), hesCurve)
           || containsSomeOf(tetMesh.incoming_halfedges(vShared), hesCurve);
}

void EmbeddingCollapser::refineAndUpdateDome(set<CH>& domeTets, DomeElements& elems)
{
    auto& tetMesh = meshProps().mesh();
    VH vCenter = _pCollapse.is_valid() ? _vPullUp : tetMesh.from_vertex_handle(_heCurrent);

    list<HEH> splitHes;
    for (VH v : elems.ceilingInnerVertices)
        if (meshProps().isInPatch(v))
            splitHes.push_back(tetMesh.find_halfedge(v, vCenter));

    for (HEH heSplit : splitHes)
    {
        VH vCeiling = tetMesh.from_vertex_handle(heSplit);
        list<HEH> hesOpp;
        for (HFH hf : tetMesh.vertex_halffaces(vCeiling))
        {
            if (elems.ceiling.count(hf) != 0)
            {
                elems.ceiling.erase(hf);
                _layerTets.erase(tetMesh.incident_cell(hf));
                domeTets.erase(tetMesh.incident_cell(hf));
                hesOpp.push_back(findMatching(tetMesh.halfface_halfedges(hf),
                                              [&](const HEH& he) {
                                                  return tetMesh.from_vertex_handle(he) != vCeiling
                                                         && tetMesh.to_vertex_handle(he) != vCeiling;
                                              }));
            }
        }
        for (EH e : tetMesh.vertex_edges(vCeiling))
            if (elems.ceilingInnerEdges.count(e) != 0)
                elems.ceilingInnerEdges.erase(e);
        elems.ceilingInnerVertices.erase(vCeiling);

        VH vNew = splitHalfEdge(heSplit, *tetMesh.hec_iter(heSplit), Q(0.5));
        elems.ceilingInnerVertices.insert(vNew);
        for (HEH he : hesOpp)
        {
            elems.ceilingInnerEdges.insert(
                tetMesh.edge_handle(tetMesh.find_halfedge(tetMesh.from_vertex_handle(he), vNew)));
            HFH hf = tetMesh.find_halfface({vNew, tetMesh.from_vertex_handle(he), tetMesh.to_vertex_handle(he)});
            assert(hf.is_valid());
            elems.ceiling.insert(hf);
            _layerTets.insert(tetMesh.incident_cell(hf));
            domeTets.insert(tetMesh.incident_cell(hf));
        }
    }

    list<EH> splitHfs;
    for (EH e : elems.ceilingInnerEdges)
        if (meshProps().isInPatch(e))
            splitHfs.push_back(e);

    for (EH eCeiling : splitHfs)
    {
        auto eVs = tetMesh.edge_vertices(eCeiling);
        HFH hfSplit = tetMesh.find_halfface({eVs[0], eVs[1], vCenter});
        elems.ceilingInnerEdges.erase(eCeiling);
        list<array<VH, 3>> hfVs;
        for (HFH hf : tetMesh.edge_halffaces(eCeiling))
        {
            if (elems.ceiling.count(hf) != 0)
            {
                elems.ceiling.erase(hf);
                _layerTets.erase(tetMesh.incident_cell(hf));
                domeTets.erase(tetMesh.incident_cell(hf));
                hfVs.push_back(meshProps().get_halfface_vertices(hf));
            }
        }

        VH vNew = splitFace(tetMesh.face_handle(hfSplit), {Q(1, 3), Q(1, 3), Q(1, 3)});

        elems.ceilingInnerVertices.insert(vNew);
        for (auto vs : hfVs)
        {
            VH vOpp;
            vector<vector<VH>> hfsNew(2);
            for (VH v : vs)
            {
                if (v == eVs[0])
                {
                    hfsNew[0].push_back(v);
                    hfsNew[1].push_back(vNew);
                }
                else if (v == eVs[1])
                {
                    hfsNew[0].push_back(vNew);
                    hfsNew[1].push_back(v);
                }
                else
                {
                    hfsNew[0].push_back(v);
                    hfsNew[1].push_back(v);
                    vOpp = v;
                }
            }
            assert(tetMesh.find_halfedge(vOpp, vNew).is_valid());
            elems.ceilingInnerEdges.insert(tetMesh.edge_handle(tetMesh.find_halfedge(vOpp, vNew)));
            for (auto vsNew : hfsNew)
            {
                HFH hfNew = tetMesh.find_halfface(vsNew);
                assert(hfNew.is_valid());
                assert(!tetMesh.is_boundary(hfNew));
                elems.ceiling.insert(hfNew);
                _layerTets.insert(tetMesh.incident_cell(hfNew));
                domeTets.insert(tetMesh.incident_cell(hfNew));
            }
        }
    }

    list<HFH> splitTets;
    for (HFH hf : elems.ceiling)
        if (meshProps().isInPatch(tetMesh.face_handle(hf)))
            splitTets.push_back(hf);

    for (HFH hfCeiling : splitTets)
    {
        elems.ceiling.erase(hfCeiling);
        _layerTets.erase(tetMesh.incident_cell(hfCeiling));
        domeTets.erase(tetMesh.incident_cell(hfCeiling));
        CH tetSplit = tetMesh.incident_cell(hfCeiling);

        VH vNew = splitTet(tetSplit, {Q(1, 4), Q(1, 4), Q(1, 4), Q(1, 4)});

        for (VH v : meshProps().get_halfface_vertices(hfCeiling))
            elems.ceilingInnerEdges.insert(tetMesh.edge_handle(tetMesh.find_halfedge(v, vNew)));
        for (HEH he : tetMesh.halfface_halfedges(hfCeiling))
        {
            HFH hfNew = tetMesh.find_halfface({tetMesh.from_vertex_handle(he), tetMesh.to_vertex_handle(he), vNew});
            assert(hfNew.is_valid());
            _layerTets.insert(tetMesh.incident_cell(hfNew));
            domeTets.insert(tetMesh.incident_cell(hfNew));
            elems.ceiling.insert(hfNew);
        }
        elems.ceilingInnerVertices.insert(vNew);
    }
}

void EmbeddingCollapser::refineAndUpdateTriangleFan(const HEH& heArc, vector<HFH>& triangleFan)
{
    auto& tetMesh = meshProps().mesh();
    HEH he0 = _heCurrent;
    VH vCollapseFrom = tetMesh.from_vertex_handle(he0);

    bool throughPatch = meshProps().isInPatch(tetMesh.face_handle(triangleFan.front()));

    // if v (except first and last) on outer side of triangle fan touches an arc/patch
    // split halfedge between vCollapseFrom and v and update triangle fan
    list<VH> vs;
    for (HFH hf : triangleFan)
    {
        for (HEH he : tetMesh.halfface_halfedges(hf))
            if (tetMesh.from_vertex_handle(he) != vCollapseFrom && tetMesh.to_vertex_handle(he) != vCollapseFrom
                && tetMesh.to_vertex_handle(he) != tetMesh.from_vertex_handle(heArc))
                vs.push_back(tetMesh.to_vertex_handle(he));
    }

    for (VH& v : vs)
    {
        bool needsSplit = (throughPatch && meshProps().isInArc(v)) || (!throughPatch && meshProps().isInPatch(v));
        if (needsSplit)
        {
            VH vOld = v;
            HEH he = tetMesh.find_halfedge(v, vCollapseFrom);
            v = splitHalfEdge(he, *tetMesh.hec_iter(he), Q(0.5));
            meshProps().replaceByChildren(_layerFaces);
            meshProps().replaceByChildren(_layerTets);
            list<FH> fsToErase;
            for (FH f : _layerFaces)
                if (contains(tetMesh.face_vertices(f), vOld))
                    fsToErase.push_back(f);
            for (FH f : fsToErase)
                _layerFaces.erase(f);
            list<CH> tetsToErase;
            for (CH tet : _layerTets)
                if (contains(tetMesh.tet_vertices(tet), vOld))
                    tetsToErase.push_back(tet);
            for (CH tet : tetsToErase)
                _layerTets.erase(tet);
        }
    }

    vs.push_front(tetMesh.to_vertex_handle(he0));
    vs.push_back(tetMesh.from_vertex_handle(heArc));

    if (vs.size() == 2)
    {
        assert(triangleFan.size() == 1);
        HEH he = tetMesh.find_halfedge(vs.front(), vs.back());
        bool doSplit = false;
        if (throughPatch)
            doSplit = meshProps().isInArc(he);
        else
            doSplit = containsMatching(tetMesh.halfedge_faces(he),
                                       [this](const FH& f) { return meshProps().isInPatch(f); });
        if (doSplit)
        {
            VH vFront = vs.front();
            vs.pop_front();
            vs.push_front(splitFace(tetMesh.face_handle(triangleFan.front()), {Q(1, 3), Q(1, 3), Q(1, 3)}));
            vs.push_front(vFront);

            meshProps().replaceByChildren(_layerFaces);
            meshProps().replaceByChildren(_layerTets);
            list<FH> fsToErase;
            for (FH f : _layerFaces)
                if (!contains(tetMesh.face_vertices(f), vCollapseFrom))
                    fsToErase.push_back(f);
            for (FH f : fsToErase)
                _layerFaces.erase(f);
            list<CH> tetsToErase;
            for (CH tet : _layerTets)
                if (!contains(tetMesh.tet_vertices(tet), vCollapseFrom))
                    tetsToErase.push_back(tet);
            for (CH tet : tetsToErase)
                _layerTets.erase(tet);
        }
    }

    triangleFan.clear();
    for (auto it = vs.begin(); it != (--vs.end());)
    {
        VH v = *it;
        VH vNext = *(++it);
        assert(tetMesh.find_halfface({v, vNext, vCollapseFrom}).is_valid());
        triangleFan.push_back(tetMesh.find_halfface({v, vNext, vCollapseFrom}));
    }
}

void EmbeddingCollapser::refineAndUpdateTetFan(const HEH& hePivot,
                                               const HFH& hfFlip,
                                               HFH& hfPatch,
                                               list<HFH>& hfsBetween)
{
    auto& tetMesh = meshProps().mesh();
    FH fPatch = tetMesh.face_handle(hfPatch);
    FH fFlip = tetMesh.face_handle(hfFlip);
    // Split all faces, insert smaller halffaces that are still incident on he
    for (HFH& hf : hfsBetween)
    {
        VH vOpp = tetMesh.to_vertex_handle(tetMesh.next_halfedge_in_halfface(hePivot, hf));
        bool onWall = false;
        CH tet = tetMesh.incident_cell(hf);
        if (hf == hfPatch && hfsBetween.size() == 1)
            onWall = containsMatching(tetMesh.cell_faces(tet),
                                      [&, this](const FH& f)
                                      { return f != fPatch && f != fFlip && meshProps().isInPatch(f); })
                     || containsMatching(tetMesh.cell_edges(tet),
                                         [&, this](const EH& e)
                                         {
                                             if (!contains(tetMesh.halfface_edges(hf), e)
                                                 && !contains(tetMesh.halfface_edges(hfFlip), e))
                                             {
                                                 int nPatches = 0;
                                                 bool onTet = false;
                                                 for (FH f : tetMesh.edge_faces(e))
                                                     if (meshProps().isInPatch(f))
                                                     {
                                                         nPatches++;
                                                         auto tets = tetMesh.face_cells(f);
                                                         if (tets[0] == tet || tets[1] == tet)
                                                             onTet = true;
                                                     }
                                                 if (!onTet && nPatches > 0)
                                                     return true;
                                             }
                                             return false;
                                         });
        else
            onWall = meshProps().isInPatch(vOpp);
        if (onWall)
        {
            VH v = splitFace(tetMesh.face_handle(hf), {Q(1, 3), Q(1, 3), Q(1, 3)});
            HFH hfNew
                = tetMesh.find_halfface({tetMesh.from_vertex_handle(hePivot), tetMesh.to_vertex_handle(hePivot), v});
            if (hf == hfPatch)
                hfPatch = hfNew;
            hf = hfNew;
        }
    }
}

void EmbeddingCollapser::patchDomeShift(const set<CH>& domeTets, const set<HFH>& floor, const set<HFH>& ceiling)
{
    auto& tetMesh = meshProps().mesh();

    FH p = meshProps().get<MC_PATCH>(tetMesh.face_handle(*floor.begin()));
    auto& pHfs = mcMeshProps().ref<PATCH_MESH_HALFFACES>(p);
    bool revert = pHfs.count(*floor.begin()) == 0;

    CH bFrom = meshProps().get<MC_BLOCK>(tetMesh.incident_cell(*floor.begin()));
    CH bTo = meshProps().get<MC_BLOCK>(tetMesh.incident_cell(tetMesh.opposite_halfface_handle(*floor.begin())));
    Transition transFromTo = meshProps().hfTransition<TRANSITION>(*floor.begin());

    for (HFH hf : floor)
    {
        FH f = tetMesh.face_handle(hf);
        pHfs.erase(revert ? tetMesh.opposite_halfface_handle(hf) : hf);
        meshProps().reset<IS_WALL>(f);
        meshProps().reset<MC_PATCH>(f);
        meshProps().reset<TRANSITION>(f);
    }
    for (HFH hf : ceiling)
    {
        FH f = tetMesh.face_handle(hf);
        pHfs.insert(revert ? hf : tetMesh.opposite_halfface_handle(hf));
        meshProps().set<IS_WALL>(f, true);
        meshProps().set<MC_PATCH>(f, p);
        meshProps().setTransition<TRANSITION>(hf, transFromTo.invert());
    }

    for (CH tet : domeTets)
    {
        assert(!tetMesh.is_deleted(tet));
        meshProps().set<MC_BLOCK>(tet, bTo);
        mcMeshProps().ref<BLOCK_MESH_TETS>(bFrom).erase(tet);
        mcMeshProps().ref<BLOCK_MESH_TETS>(bTo).insert(tet);
        _layerTets.erase(tet);

        for (auto& kv : meshProps().ref<CHART>(tet))
            kv.second = transFromTo.apply(kv.second);
    }
}

} // namespace c4hex
