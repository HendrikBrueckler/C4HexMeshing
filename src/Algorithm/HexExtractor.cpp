#include "C4Hex/Algorithm/HexExtractor.hpp"

namespace c4hex
{

template <typename MESHPROPS>
HexExtractor<MESHPROPS>::HexExtractor(const TetMeshProps& meshProps, MESHPROPS& hexMeshProps)
    : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps), _hexMeshProps(hexMeshProps)
{
}

template <typename MESHPROPS>
typename HexExtractor<MESHPROPS>::RetCode HexExtractor<MESHPROPS>::HexExtractor::extractHexMesh(int subdiv)
{
    if (!std::is_same<MESHPROPS, PolyMeshProps>::value)
        subdiv = 0;

    int stepsPerInt = std::max(1, subdiv + 1);
    Q stepSize = Q(1, stepsPerInt);

    _hexMeshProps.mesh().clear();

    _hexMeshProps.template allocate<MC_NODE_ID>(-1);
    _hexMeshProps.template allocate<MC_ARC_ID>(-1);
    _hexMeshProps.template allocate<MC_PATCH_ID>(-1);
    _hexMeshProps.template allocate<MC_BLOCK_ID>(-1);
    _hexMeshProps.template allocate<IS_SINGULAR>(false);
    _hexMeshProps.template allocate<IS_FEATURE_F>(0);
    _hexMeshProps.template allocate<IS_FEATURE_E>(0);
    _hexMeshProps.template allocate<IS_FEATURE_V>(0);
    _hexMeshProps.template allocate<MARK_N>(0);
    _hexMeshProps.template allocate<MARK_A>(0);
    _hexMeshProps.template allocate<MARK_P>(0);
    _hexMeshProps.template allocate<MARK_B>(0);

    auto n2hexV = createNodeHexV();
    auto a2delta2hexV = createArcHexVE(n2hexV, subdiv);
    auto p2igm2hexV = createPatchHexVEF(a2delta2hexV, subdiv);
    createBlockHexVEFC(p2igm2hexV, subdiv);

    return SUCCESS;
}

template <typename MESHPROPS>
map<VH, VH> HexExtractor<MESHPROPS>::createNodeHexV()
{
    map<VH, VH> n2hexV;

    for (VH n : mcMeshProps().mesh().vertices())
    {
        VH v = _hexMeshProps.mesh().add_vertex(meshProps().mesh().vertex(mcMeshProps().get<NODE_MESH_VERTEX>(n)));
        n2hexV[n] = v;
        if (mcMeshProps().isAllocated<IS_FEATURE_V>())
            _hexMeshProps.template set<IS_FEATURE_V>(v, mcMeshProps().get<IS_FEATURE_V>(n));
        if (mcMeshProps().isAllocated<MARK_N>())
            _hexMeshProps.template set<MARK_N>(v, mcMeshProps().get<MARK_N>(n));
        _hexMeshProps.template set<MC_NODE_ID>(v, n.idx());
    }

    assert(n2hexV.size() == mcMeshProps().mesh().n_logical_vertices());

    return n2hexV;
}

template <typename MESHPROPS>
map<EH, vector<VH>> HexExtractor<MESHPROPS>::createArcHexVE(const map<VH, VH>& n2hexV, int subdiv)
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& tetMesh = meshProps().mesh();

    int stepsPerInt = std::max(1, subdiv + 1);
    Q stepSize = Q(1, stepsPerInt);

    map<EH, vector<VH>> a2delta2hexV;

    for (EH a : mc.edges())
    {
        HEH ha = mc.halfedge_handle(a, 0);
        VH nFrom = mc.from_vertex_handle(ha);
        VH nTo = mc.to_vertex_handle(ha);

        auto& hes = mcMeshProps().ref<ARC_MESH_HALFEDGES>(a);
        auto itHeCurrent = hes.begin();
        CH b = meshProps().get<MC_BLOCK>(*tetMesh.hec_iter(*itHeCurrent));
        UVWDir dir = halfarcDirInBlock(ha, b);
        bool negDir = isNeg(dir);
        int deltaCoord = toCoord(dir);
        VH vStart = tetMesh.from_vertex_handle(*itHeCurrent);
        CH tetStart = anyIncidentTetOfBlock(vStart, b);
        Vec3Q igmStart = meshProps().ref<CHART_IGM>(tetStart).at(vStart);

        int arcLen = mcMeshProps().get<ARC_INT_LENGTH>(a);
        auto& delta2hexV = a2delta2hexV[a];
        delta2hexV = vector<VH>(stepsPerInt * arcLen + 1);
        delta2hexV[0] = n2hexV.at(nFrom);
        delta2hexV[arcLen * stepsPerInt] = n2hexV.at(nTo);
        for (int delta = 0; delta < arcLen; delta++)
        {
            for (int step = 0; step < stepsPerInt; step++)
            {
                if (delta == 0 && step == 0)
                    continue;

                Vec3Q igmCurrent = igmStart;
                if (negDir)
                    igmCurrent[deltaCoord] -= (delta + step * stepSize);
                else
                    igmCurrent[deltaCoord] += (delta + step * stepSize);
                VH vFrom(-1);
                VH vTo(-1);
                Q t;
                do
                {
                    HEH heCurrent = *itHeCurrent;
                    vFrom = tetMesh.from_vertex_handle(heCurrent);
                    vTo = tetMesh.to_vertex_handle(heCurrent);
                    CH tetFrom = anyIncidentTetOfBlock(vFrom, b);
                    CH tetTo = anyIncidentTetOfBlock(vTo, b);
                    Q igmFrom = meshProps().ref<CHART_IGM>(tetFrom).at(vFrom)[deltaCoord];
                    Q igmTo = meshProps().ref<CHART_IGM>(tetTo).at(vTo)[deltaCoord];
                    t = (Q(igmCurrent[deltaCoord]) - igmFrom) / (igmTo - igmFrom);
                    if (t > 1)
                        itHeCurrent++;
                    assert(itHeCurrent != hes.end());
                } while (t > 1);
                assert(t >= 0 && t <= 1);
                auto xyz = t * Vec3Q(tetMesh.vertex(vTo)) + (Q(1) - t) * Vec3Q(tetMesh.vertex(vFrom));
                delta2hexV[delta * stepsPerInt + step] = _hexMeshProps.mesh().add_vertex(Vec3Q2d(xyz));
            }
        }
        assert((int)delta2hexV.size() == stepsPerInt * arcLen + 1);
    }

#ifndef NDEBUG
    for (const auto& kv : a2delta2hexV)
        for (VH v : kv.second)
            assert(_hexMeshProps.mesh().vertex(v)[0] != DBL_MAX);
#endif

    for (auto a2delta2v : a2delta2hexV)
    {
        EH a = a2delta2v.first;
        auto delta2v = a2delta2v.second;
        for (int i = 0; i < (int)delta2v.size() - 1; i++)
        {
            EH e = _hexMeshProps.mesh().add_edge(delta2v[i], delta2v[i + 1]);
            if (mcMeshProps().isAllocated<IS_FEATURE_E>())
                _hexMeshProps.template set<IS_FEATURE_E>(e, mcMeshProps().get<IS_FEATURE_E>(a));
            if (mcMeshProps().isAllocated<IS_SINGULAR>() && mcMeshProps().get<IS_SINGULAR>(a))
                _hexMeshProps.template set<IS_SINGULAR>(e, true);
            if (mcMeshProps().isAllocated<MARK_A>())
                _hexMeshProps.template set<MARK_A>(e, mcMeshProps().get<MARK_A>(a));
            _hexMeshProps.template set<MC_ARC_ID>(e, a.idx());
        }
    }

    return a2delta2hexV;
}

template <typename MESHPROPS>
map<FH, map<Vec3Q, VH>> HexExtractor<MESHPROPS>::createPatchHexVEF(const map<EH, vector<VH>>& a2delta2hexV, int subdiv)
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& tetMesh = meshProps().mesh();
    map<FH, map<Vec3Q, VH>> p2igm2hexV;

    int stepsPerInt = std::max(1, subdiv + 1);
    Q stepSize(1, stepsPerInt);

    for (FH p : mc.faces())
    {
        HFH hp = mc.halfface_handle(p, 0);
        assert(!mc.is_boundary(hp));
        bool flipped = false;
        if (mc.is_boundary(hp))
        {
            flipped = true;
            hp = mc.opposite_halfface_handle(hp);
        }
        CH b = mc.incident_cell(hp);

        auto cornersHp = orderedHalfpatchCorners(hp);
        auto dir2orderedHas = halfpatchHalfarcsByDir(hp);

        auto& igm2hexV = p2igm2hexV[p];
        for (const auto& kv : dir2orderedHas)
        {
            auto& has = kv.second;
            UVWDir dir = kv.first;
            assert(dim(dir) == 1);
            int deltaCoord = toCoord(dir);
            for (HEH ha : has)
            {
                if (ha.idx() % 2 != 0)
                    ha = mc.opposite_halfedge_handle(ha);
                EH a = mc.edge_handle(ha);
                auto arcLen = mcMeshProps().get<ARC_INT_LENGTH>(a);

                VH nFrom = mc.from_vertex_handle(ha);
                VH vFrom = mcMeshProps().get<NODE_MESH_VERTEX>(nFrom);
                CH tetFrom = anyIncidentTetOfBlock(vFrom, b);
                Vec3Q igmFrom = (meshProps().ref<CHART_IGM>(tetFrom).at(vFrom));

                VH nTo = mc.to_vertex_handle(ha);
                VH vTo = mcMeshProps().get<NODE_MESH_VERTEX>(nTo);
                CH tetTo = anyIncidentTetOfBlock(vTo, b);
                Vec3Q igmTo = (meshProps().ref<CHART_IGM>(tetTo).at(vTo));

                const auto& delta2hexv = a2delta2hexV.at(a);
                Vec3Q igmBetween = igmFrom;

                igm2hexV[igmTo] = delta2hexv.at(arcLen * stepsPerInt);
                for (int delta = 0; delta < arcLen; delta++)
                {
                    for (int step = 0; step < stepsPerInt; step++)
                    {
                        igm2hexV[igmBetween] = delta2hexv.at(delta * stepsPerInt + step);
                        igmBetween[deltaCoord] += (igmTo[deltaCoord] > igmFrom[deltaCoord] ? stepSize : -stepSize);
                    }
                }
                assert(igmBetween == igmTo);
            }
        }

        {
            pairTT<int> sideLengths;

            VH vCorner0 = mcMeshProps().get<NODE_MESH_VERTEX>(cornersHp[0]);
            CH tetCorner0 = anyIncidentTetOfBlock(vCorner0, b);
            Vec3i igmCorner0 = Vec3Q2i(meshProps().ref<CHART_IGM>(tetCorner0).at(vCorner0));
            VH vCorner1 = mcMeshProps().get<NODE_MESH_VERTEX>(cornersHp[1]);
            CH tetCorner1 = anyIncidentTetOfBlock(vCorner1, b);
            Vec3i igmCorner1 = Vec3Q2i(meshProps().ref<CHART_IGM>(tetCorner1).at(vCorner1));
            VH vCorner3 = mcMeshProps().get<NODE_MESH_VERTEX>(cornersHp[3]);
            CH tetCorner3 = anyIncidentTetOfBlock(vCorner3, b);
            Vec3i igmCorner3 = Vec3Q2i(meshProps().ref<CHART_IGM>(tetCorner3).at(vCorner3));

            int deltaCoord1 = -1;
            int deltaCoord3 = -1;
            for (int coord = 0; coord < 3; coord++)
            {
                if (igmCorner1[coord] - igmCorner0[coord] != 0)
                    deltaCoord1 = coord;
                else if (igmCorner3[coord] - igmCorner0[coord] != 0)
                    deltaCoord3 = coord;
            }
            assert(deltaCoord1 != -1);
            assert(deltaCoord3 != -1);

            sideLengths.first = std::abs((igmCorner1 - igmCorner0)[deltaCoord1]);
            sideLengths.second = std::abs((igmCorner3 - igmCorner0)[deltaCoord3]);
            assert((int)igm2hexV.size() == sideLengths.first * stepsPerInt * 2 + sideLengths.second * stepsPerInt * 2);

            Vec3Q igmBetween(igmCorner0);
            for (int deltaDir1 = 0; deltaDir1 < sideLengths.first; deltaDir1++)
            {
                for (int stepDir1 = (deltaDir1 == 0 ? 1 : 0); stepDir1 < stepsPerInt; stepDir1++)
                {
                    igmBetween[deltaCoord1] = igmCorner0[deltaCoord1];
                    if (igmCorner1[deltaCoord1] > igmCorner0[deltaCoord1])
                        igmBetween[deltaCoord1] += (deltaDir1 + stepDir1 * stepSize);
                    else
                        igmBetween[deltaCoord1] -= (deltaDir1 + stepDir1 * stepSize);
                    for (int deltaDir3 = 0; deltaDir3 < sideLengths.second; deltaDir3++)
                    {
                        for (int stepDir3 = (deltaDir3 == 0 ? 1 : 0); stepDir3 < stepsPerInt; stepDir3++)
                        {
                            igmBetween[deltaCoord3] = igmCorner0[deltaCoord3];
                            if (igmCorner3[deltaCoord3] > igmCorner0[deltaCoord3])
                                igmBetween[deltaCoord3] += (deltaDir3 + stepDir3 * stepSize);
                            else
                                igmBetween[deltaCoord3] -= (deltaDir3 + stepDir3 * stepSize);
                            igm2hexV[igmBetween] = _hexMeshProps.mesh().add_vertex(Vec3d(DBL_MAX, DBL_MAX, DBL_MAX));
                        }
                    }
                }
            }
            assert((int)igm2hexV.size()
                   == (sideLengths.first * stepsPerInt + 1) * (sideLengths.second * stepsPerInt + 1));

            {
                for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
                {
                    if (flipped)
                        hf = tetMesh.opposite_halfface_handle(hf);
                    CH tet = tetMesh.incident_cell(hf);
                    assert(tet.is_valid());
                    Vec3d bboxMin(DBL_MAX, DBL_MAX, DBL_MAX);
                    Vec3d bboxMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
                    for (VH v : meshProps().get_halfface_vertices(hf))
                    {
                        Vec3d igm = Vec3Q2d(meshProps().ref<CHART_IGM>(tet).at(v));
                        for (int coord = 0; coord < 3; coord++)
                        {
                            bboxMin[coord] = std::min(igm[coord], bboxMin[coord]);
                            bboxMax[coord] = std::max(igm[coord], bboxMax[coord]);
                        }
                    }
                    Vec3Q igm(igmCorner0);
                    for (int valDir1 = std::ceil(bboxMin[deltaCoord1] * stepsPerInt);
                         valDir1 <= std::floor(bboxMax[deltaCoord1] * stepsPerInt);
                         valDir1++)
                    {
                        igm[deltaCoord1] = Q(valDir1, stepsPerInt);
                        // Assign xyz only to patch-interior vertices, skip edge vertices (already assigned)
                        if (igm[deltaCoord1] == igmCorner0[deltaCoord1] || igm[deltaCoord1] == igmCorner1[deltaCoord1])
                            continue;
                        for (int valDir3 = std::ceil(bboxMin[deltaCoord3] * stepsPerInt);
                             valDir3 <= std::floor(bboxMax[deltaCoord3] * stepsPerInt);
                             valDir3++)
                        {
                            igm[deltaCoord3] = Q(valDir3, stepsPerInt);
                            // Assign xyz only to patch-interior vertices, skip edge vertices (already assigned)
                            if (igm[deltaCoord3] == igmCorner0[deltaCoord3]
                                || igm[deltaCoord3] == igmCorner3[deltaCoord3])
                                continue;

                            Vec3Q barCoords(0, 0, 0);
                            if (barycentricCoordsIGM(hf, igm, deltaCoord1, deltaCoord3, barCoords))
                            {
                                auto vs = meshProps().get_halfface_vertices(hf);
                                Vec3Q xyz(0, 0, 0);
                                for (int i = 0; i < 3; i++)
                                    xyz += barCoords[i] * Vec3Q(tetMesh.vertex(vs[i]));
                                _hexMeshProps.mesh().template set_vertex(igm2hexV.at(igm), Vec3Q2d(xyz));
                            }
                        }
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    for (const auto& kv : p2igm2hexV)
        for (const auto& kv2 : kv.second)
            assert(_hexMeshProps.mesh().vertex(kv2.second)[0] != DBL_MAX);
#endif

    for (auto p2igm2v : p2igm2hexV)
    {
        FH p = p2igm2v.first;
        auto igm2v = p2igm2v.second;
        Vec3Q minIGM = igm2v.begin()->first;
        Vec3Q maxIGM = igm2v.rbegin()->first;

        vector<int> varCoords;
        for (int i = 0; i < 3; i++)
            if (minIGM[i] != maxIGM[i])
                varCoords.emplace_back(i);

        for (Q varVal1 = minIGM[varCoords[0]]; varVal1 < maxIGM[varCoords[0]]; varVal1 += stepSize)
            for (Q varVal2 = minIGM[varCoords[1]]; varVal2 < maxIGM[varCoords[1]]; varVal2 += stepSize)
            {
                Vec3Q UVW = minIGM;
                UVW[varCoords[0]] = varVal1;
                UVW[varCoords[1]] = varVal2;
                vector<Vec3Q> uvws({UVW, UVW, UVW, UVW});
                uvws[1][varCoords[0]] += stepSize;
                uvws[2][varCoords[0]] += stepSize;
                uvws[2][varCoords[1]] += stepSize;
                uvws[3][varCoords[1]] += stepSize;
                vector<VH> vs;
                for (Vec3Q uvw : uvws)
                    vs.emplace_back(igm2v.at(uvw));
                FH f = _hexMeshProps.mesh().add_face(vs);
                if (mcMeshProps().isAllocated<IS_FEATURE_F>())
                    _hexMeshProps.template set<IS_FEATURE_F>(f, mcMeshProps().get<IS_FEATURE_F>(p));
                if (mcMeshProps().isAllocated<MARK_P>())
                    _hexMeshProps.template set<MARK_P>(f, mcMeshProps().get<MARK_P>(p));
                _hexMeshProps.template set<MC_PATCH_ID>(f, p.idx());
            }
    }

    return p2igm2hexV;
}

template <typename MESHPROPS>
map<CH, map<Vec3Q, VH>> HexExtractor<MESHPROPS>::createBlockHexVEFC(const map<FH, map<Vec3Q, VH>>& p2igm2hexV,
                                                                    int subdiv)
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& tetMesh = meshProps().mesh();
    map<CH, map<Vec3Q, VH>> b2igm2hexV;

    int stepsPerInt = std::max(1, subdiv + 1);
    Q stepSize(1, stepsPerInt);

    for (CH b : mc.cells())
    {
        auto& igm2hexV = b2igm2hexV[b];

        Vec3i minIGM
            = Vec3Q2i(nodeIGMinBlock(mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::NEG_U_NEG_V_NEG_W), b));
        Vec3i maxIGM
            = Vec3Q2i(nodeIGMinBlock(mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::POS_U_POS_V_POS_W), b));

        for (int U = minIGM[0]; U <= maxIGM[0]; U++)
            for (int stepU = 0; stepU < (U == maxIGM[0] ? 1 : stepsPerInt); stepU++)
                for (int V = minIGM[1]; V <= maxIGM[1]; V++)
                    for (int stepV = 0; stepV < (V == maxIGM[1] ? 1 : stepsPerInt); stepV++)
                        for (int W = minIGM[2]; W <= maxIGM[2]; W++)
                            for (int stepW = 0;
                                 stepW < (W == maxIGM[2] || (stepU != 0 && stepV != 0) ? 1 : stepsPerInt);
                                 stepW++)
                            {
                                Vec3Q igm(U + stepU * stepSize, V + stepV * stepSize, W + stepW * stepSize);
                                bool patchLookup = true;
                                UVWDir lookupDir = UVWDir::NONE;
                                if (U == minIGM[0] && stepU == 0)
                                    lookupDir = UVWDir::NEG_U;
                                else if (U == maxIGM[0] && stepU == 0)
                                    lookupDir = UVWDir::POS_U;
                                else if (V == minIGM[1] && stepV == 0)
                                    lookupDir = UVWDir::NEG_V;
                                else if (V == maxIGM[1] && stepV == 0)
                                    lookupDir = UVWDir::POS_V;
                                else if (W == minIGM[2] && stepW == 0)
                                    lookupDir = UVWDir::NEG_W;
                                else if (W == maxIGM[2] && stepW == 0)
                                    lookupDir = UVWDir::POS_W;
                                else
                                    patchLookup = false;

                                if (!patchLookup)
                                    igm2hexV[igm] = _hexMeshProps.mesh().add_vertex(Vec3d(DBL_MAX, DBL_MAX, DBL_MAX));
                                else
                                {
                                    const auto& patches = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b).at(lookupDir);
#ifndef NDEBUG
                                    bool found = false;
#endif
                                    for (FH p : patches)
                                    {
                                        assert(p2igm2hexV.find(p) != p2igm2hexV.end());
                                        auto& patchigm2hexV = p2igm2hexV.find(p)->second;
                                        bool mainHP = mc.incident_cell(mc.halfface_handle(p, 0)) == b;
                                        Vec3Q transformedIGM = igm;
                                        if (!mainHP)
                                        {
                                            Transition trans = mcMeshProps().get<PATCH_IGM_TRANSITION>(p);
                                            transformedIGM = trans.invert().apply(Vec3Q(igm));
                                        }
                                        auto it = patchigm2hexV.find(transformedIGM);
                                        if (it != patchigm2hexV.end())
                                        {
#ifndef NDEBUG
                                            found = true;
#endif
                                            igm2hexV[igm] = it->second;
                                        }
                                    }
                                    assert(found);
                                }
                            }

        for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
        {
            if (doubleVolumeIGM(tet) < 1e-6 && rationalVolumeIGM(tet) <= 0)
                continue;

            Vec3Q bboxMin(DBL_MAX, DBL_MAX, DBL_MAX);
            Vec3Q bboxMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
            for (VH v : tetMesh.tet_vertices(tet))
            {
                Vec3Q igm = meshProps().ref<CHART_IGM>(tet).at(v);
                for (int coord = 0; coord < 3; coord++)
                {
                    bboxMin[coord] = std::min(igm[coord], bboxMin[coord]);
                    bboxMax[coord] = std::max(igm[coord], bboxMax[coord]);
                }
            }
            if (std::ceil(bboxMin[0].get_d()) > std::floor(bboxMax[0].get_d())
                && std::ceil(bboxMin[1].get_d()) > std::floor(bboxMax[1].get_d())
                && std::ceil(bboxMin[2].get_d()) > std::floor(bboxMax[2].get_d()))
                continue;
            Q minU = Q(std::ceil(Q(bboxMin[0] * stepsPerInt).get_d())) / stepsPerInt;
            Q maxU = Q(std::floor(Q(bboxMax[0] * stepsPerInt).get_d())) / stepsPerInt;
            Q minV = Q(std::ceil(Q(bboxMin[1] * stepsPerInt).get_d())) / stepsPerInt;
            Q maxV = Q(std::floor(Q(bboxMax[1] * stepsPerInt).get_d())) / stepsPerInt;
            Q minW = Q(std::ceil(Q(bboxMin[2] * stepsPerInt).get_d())) / stepsPerInt;
            Q maxW = Q(std::floor(Q(bboxMax[2] * stepsPerInt).get_d())) / stepsPerInt;
            for (Q U = minU; U <= maxU; U += stepSize)
            {
                if (U == minIGM[0] || U == maxIGM[0])
                    continue;
                for (Q V = minV; V <= maxV; V += stepSize)
                {
                    if (V == minIGM[1] || V == maxIGM[1])
                        continue;
                    for (Q W = minW; W <= maxW; W += stepSize)
                    {
                        if (W == minIGM[2] || W == maxIGM[2])
                            continue;
                        double Ud = U.get_d();
                        double Vd = V.get_d();
                        double Wd = W.get_d();
                        if (std::trunc(Ud) != Ud && std::trunc(Vd) != Vd && std::trunc(Wd) != Wd)
                            continue;
                        Vec3Q igm = Vec3Q(U, V, W);

                        Vec4Q barCoords(0, 0, 0, 0);
                        if (barycentricCoordsIGM(tet, igm, barCoords))
                        {
                            auto vsItPair = tetMesh.tet_vertices(tet);
                            vector<VH> vs(vsItPair.first, vsItPair.second);
                            Vec3Q xyz(0, 0, 0);
                            for (int i = 0; i < 4; i++)
                                xyz += barCoords[i] * Vec3Q(tetMesh.vertex(vs[i]));
                            _hexMeshProps.mesh().template set_vertex(igm2hexV.at(igm), Vec3Q2d(xyz));
                        }
                    }
                }
            }
        }

#ifndef NDEBUG
        for (const auto& kv : igm2hexV)
            assert(_hexMeshProps.mesh().vertex(kv.second)[0] != DBL_MAX);
#endif
    }

    for (const auto& kv : b2igm2hexV)
    {
        CH b = kv.first;
        auto& igm2hexV = kv.second;

        Vec3i minIGM
            = Vec3Q2i(nodeIGMinBlock(mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::NEG_U_NEG_V_NEG_W), b));
        Vec3i maxIGM
            = Vec3Q2i(nodeIGMinBlock(mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::POS_U_POS_V_POS_W), b));

        for (int U = minIGM[0]; U < maxIGM[0]; U++)
            for (int V = minIGM[1]; V < maxIGM[1]; V++)
                for (int W = minIGM[2]; W < maxIGM[2]; W++)
                {
                    Vec3Q UVW(U, V, W);
                    vector<HFH> hfs;
                    for (int quantCoord : {0, 1, 2})
                    {
                        for (int quantVal : {0, 1})
                        {
                            int nonQuantCoord1 = (quantCoord + 2) % 3;
                            int nonQuantCoord2 = (quantCoord + 1) % 3;
                            if (quantVal == 0)
                                std::swap(nonQuantCoord1, nonQuantCoord2);
                            vector<Vec3Q> corners({{UVW}, {UVW}, {UVW}, {UVW}});
                            corners[1][nonQuantCoord1] += stepSize;
                            corners[2][nonQuantCoord1] += stepSize;
                            corners[2][nonQuantCoord2] += stepSize;
                            corners[3][nonQuantCoord2] += stepSize;
                            for (int sub1 = 0; sub1 < stepsPerInt; sub1++)
                            {
                                for (int sub2 = 0; sub2 < stepsPerInt; sub2++)
                                {
                                    Vec3Q offset(0, 0, 0);
                                    offset[quantCoord] = quantVal;
                                    offset[nonQuantCoord1] = sub1 * stepSize;
                                    offset[nonQuantCoord2] = sub2 * stepSize;
                                    vector<VH> hexVs;
                                    // int i = 0;
                                    for (Vec3Q corner : corners)
                                    {
                                        hexVs.emplace_back(igm2hexV.at(corner + offset));
                                    }

                                    HFH hf = _hexMeshProps.mesh().find_halfface(hexVs);
                                    if (!hf.is_valid())
                                        hf = _hexMeshProps.mesh().halfface_handle(_hexMeshProps.mesh().add_face(hexVs),
                                                                                  0);
                                    assert(hf.is_valid());
                                    hfs.emplace_back(hf);
                                }
                            }
                        }
                    }
                    assert((int)hfs.size() == stepsPerInt * stepsPerInt * 6);
                    CH c = _hexMeshProps.mesh().add_cell(hfs);
                    assert(c.is_valid());
                    _hexMeshProps.template set<MC_BLOCK_ID>(c, b.idx());
                    if (mcMeshProps().isAllocated<MARK_B>())
                        _hexMeshProps.template set<MARK_B>(c, mcMeshProps().get<MARK_B>(b));
                }
    }

    return b2igm2hexV;
}

template <typename MESHPROPS>
bool HexExtractor<MESHPROPS>::barycentricCoordsIGM(
    const HFH& hf, const Vec3Q& igmUVW, int coord1, int coord3, Vec3Q& barCoords) const
{
    const TetMesh& tetMesh = meshProps().mesh();

    CH tet = tetMesh.incident_cell(hf);
    auto vs = meshProps().get_halfface_vertices(hf);

    vector<Vec3Q> cornerIGM;
    vector<Vec3Q> edgeVecs;
    for (int corner = 0; corner < 3; corner++)
    {
        const auto& igmUVWfrom = meshProps().ref<CHART_IGM>(tet).at(vs[corner]);
        const auto& igmUVWto = meshProps().ref<CHART_IGM>(tet).at(vs[(corner + 1) % 3]);
        cornerIGM.emplace_back(igmUVWfrom);
        edgeVecs.emplace_back(igmUVWto - igmUVWfrom);
    }

    for (int corner = 0; corner < 2; corner++)
    {
        int edge = (corner + 1) % 3;
        assert((cornerIGM[corner][coord1] - cornerIGM[edge][coord1]) * edgeVecs[edge][coord3]
                   - (cornerIGM[corner][coord3] - cornerIGM[edge][coord3]) * edgeVecs[edge][coord1]
               != 0);
        if ((cornerIGM[corner][coord1] - cornerIGM[edge][coord1]) * edgeVecs[edge][coord3]
                - (cornerIGM[corner][coord3] - cornerIGM[edge][coord3]) * edgeVecs[edge][coord1]
            == 0)
            return false;
        barCoords[corner] = ((igmUVW[coord1] - cornerIGM[edge][coord1]) * edgeVecs[edge][coord3]
                             - (igmUVW[coord3] - cornerIGM[edge][coord3]) * edgeVecs[edge][coord1])
                            / ((cornerIGM[corner][coord1] - cornerIGM[edge][coord1]) * edgeVecs[edge][coord3]
                               - (cornerIGM[corner][coord3] - cornerIGM[edge][coord3]) * edgeVecs[edge][coord1]);
        if (barCoords[corner] < 0 || barCoords[corner] > 1)
            return false;
    }

    barCoords[2] = Q(1) - barCoords[0] - barCoords[1];
    if (barCoords[2] < 0)
        return false;

    assert(barCoords[0] + barCoords[1] + barCoords[2] == 1);
    assert(barCoords[0] >= 0 && barCoords[1] >= 0 && barCoords[2] >= 0);
    return true;
}

template <typename MESHPROPS>
bool HexExtractor<MESHPROPS>::barycentricCoordsIGM(const CH& tet, const Vec3Q& igmUVW, Vec4Q& barCoords) const
{
    const TetMesh& tetMesh = meshProps().mesh();

    vector<Vec3Q> cornerIGMs;
    vector<Vec3Q> oppHfNormals;
    for (VH v : tetMesh.tet_vertices(tet))
        cornerIGMs.emplace_back(meshProps().ref<CHART_IGM>(tet).at(v));

    for (int corner = 0; corner < 3; corner++)
    {
        Vec3Q normal = cross(cornerIGMs[(corner + 2) % 4] - cornerIGMs[(corner + 1) % 4],
                                cornerIGMs[(corner + 3) % 4] - cornerIGMs[(corner + 1) % 4]);
        oppHfNormals.emplace_back(normal);
    }

    for (int corner = 0; corner < 3; corner++)
    {
        const auto& cornerIGM = cornerIGMs[corner];
        const auto& normal = oppHfNormals[corner];
        assert(dot(cornerIGM - cornerIGMs[(corner + 1) % 4], normal) != 0);
        barCoords[corner] = dot(igmUVW - cornerIGMs[(corner + 1) % 4], normal)
                            / dot(cornerIGM - cornerIGMs[(corner + 1) % 4], normal);
        if (barCoords[corner] < 0 || barCoords[corner] > 1)
            return false;
    }

    barCoords[3] = Q(1) - barCoords[0] - barCoords[1] - barCoords[2];
    if (barCoords[3] < 0)
        return false;

    assert(barCoords[0] + barCoords[1] + barCoords[2] + barCoords[3] == 1);
    assert(barCoords[0] >= 0 && barCoords[1] >= 0 && barCoords[2] >= 0 && barCoords[3] >= 0);

    return true;
}

template class HexExtractor<HexMeshProps>;
template class HexExtractor<PolyMeshProps>;

}; // namespace c4hex
