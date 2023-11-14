#include "C4Hex/Algorithm/IGMUntangler.hpp"

#include "C4Hex/Algorithm/IGMInitializer.hpp"

#include "Eigen/Geometry"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"
#include "Eigen/StdVector"

#include <LBFGS.h>
#include <LBFGSB.h>
#include <LBFGSpp/LineSearchMoreThuente.h>

#include <chrono>
#include <iostream>

#ifdef C4HEX_WITH_TLC
#include <TLC.h>
#include <nlopt.hpp>
#endif

namespace c4hex
{

inline double calcChi(double eps, double det)
{
    eps *= eps;
    double det2 = det * det;
    if (det > 0)
        return (det + std::sqrt(eps + det2)) * .5;
    return .5 * eps / (std::sqrt(eps + det2) - det);
}

inline double calcChiDeriv(double eps, double det)
{
    return .5 + det / (2. * std::sqrt(eps * eps + det * det));
}

IGMUntangler::FoldoverEnergy::FoldoverEnergy(double areaVsAngles_,
                                             double e_,
                                             const Eigen::Matrix4Xi& tetVtxIndices_,
                                             const Eigen::Matrix4Xd& tetZ_,
                                             const Eigen::VectorXd& tetVol_,
                                             const int nInteriorVs_,
                                             const Eigen::VectorXd& tetWeights_)
    : areaVsAngles(areaVsAngles_), e(e_), tetVtxIndices(tetVtxIndices_), tetZ(tetZ_), tetVol(tetVol_),
      nInteriorVs(nInteriorVs_), tetJ(Eigen::Matrix3Xd::Zero(3, 3 * tetVtxIndices.cols())),
      tetDetJ(Eigen::VectorXd::Zero(tetVtxIndices.cols())), tetWeights(tetWeights_)
{
}

double IGMUntangler::FoldoverEnergy::operator()(const Eigen::VectorXd& UVWflat, Eigen::VectorXd& gradF)
{
    if (!updateJ(UVWflat) || !updateGradF(gradF))
        return DBL_MAX;

    return calcF();
}

bool IGMUntangler::FoldoverEnergy::updateJ(const Eigen::VectorXd& UVWflat)
{

    auto ret = IGMUntangler::updateJ(tetVtxIndices, UVWflat, tetZ, tetJ, tetDetJ, minDetJ, nFlipped);
    if (ret != SUCCESS)
        return false;
    return true;
}

bool IGMUntangler::FoldoverEnergy::updateGradF(Eigen::VectorXd& gradF)
{
    gradF.setZero();
    auto ret = IGMUntangler::updateGradF(
        tetVtxIndices, tetVol, tetWeights, tetZ, tetJ, tetDetJ, nInteriorVs, e, areaVsAngles, gradF);
    if (ret != SUCCESS)
        return false;
    return true;
}

double IGMUntangler::FoldoverEnergy::calcF()
{
    double F = 0.;
    IGMUntangler::calcF(tetVtxIndices, tetVol, tetWeights, tetJ, tetDetJ, e, areaVsAngles, F);
    if (!std::isfinite(F))
        return DBL_MAX;
    return F;
}

IGMUntangler::IGMUntangler(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps),
      MCMeshManipulator(meshProps)
{
}

IGMUntangler::RetCode IGMUntangler::untangleIGM(double areaVsAngles, int maxIter, int secondsTimeLimit)
{
    const MCMesh& mc = mcMeshProps().mesh();
    const TetMesh& mesh = meshProps().mesh();

    // Sanity check on input before trying to untangle it
    for (CH b : mcMeshProps().mesh().cells())
    {
        bool selfadjacent = false;
        for (VH n : mc.cell_vertices(b))
        {
            VH v = mcMeshProps().get<NODE_MESH_VERTEX>(n);
            CH tet = anyIncidentTetOfBlock(v, b);
            set<CH> tets;
            forVertexNeighbourTetsInBlock(v,
                                          tet,
                                          [&tets](const CH& tet2)
                                          {
                                              tets.insert(tet2);
                                              return false;
                                          });
            for (CH tet2 : mesh.vertex_cells(v))
                if (meshProps().get<MC_BLOCK>(tet2) == b && tets.count(tet2) == 0)
                    selfadjacent = true;

            if (selfadjacent)
                break;
        }
        // Sanity check not working for selfadjacent blocks
        if (selfadjacent)
            continue;

        VH nMin = mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::NEG_U_NEG_V_NEG_W);
        VH nMax = mcMeshProps().ref<BLOCK_CORNER_NODES>(b).at(UVWDir::POS_U_POS_V_POS_W);
        Vec3Q minIGM = nodeIGMinBlock(nMin, b);
        Vec3Q maxIGM = nodeIGMinBlock(nMax, b);
        Vec3Q midPt = (minIGM + maxIGM) * 0.5;
        for (HFH hp : mcMeshProps().mesh().cell_halffaces(b))
            for (HFH hf : mcMeshProps().hpHalffaces(hp))
            {
                CH tet = meshProps().mesh().incident_cell(hf);
                vector<Vec3Q> uvws;
                for (VH v : meshProps().mesh().halfface_vertices(hf))
                    uvws.push_back(meshProps().ref<CHART_IGM>(tet).at(v));
                if (dot(cross(uvws[1] - uvws[0], uvws[2] - uvws[0]), midPt - uvws[0]) <= 0)
                {
                    LOG(ERROR) << "Boundary hf " + std::to_string(hf.idx()) + " of block " + std::to_string(b.idx())
                                      + " inverted!";
                    LOG(ERROR) << "UVW coords: ";
                    for (Vec3Q uvw : uvws)
                        LOG(ERROR) << Vec3Q2d(uvw);
                    LOG(ERROR) << "Area: " << Vec3Q2d(cross(uvws[1] - uvws[0], uvws[2] - uvws[0])).norm();
                    return NUMERICAL_ISSUE;
                }
            }
    }

    // Variables for statistics
    int nBadBlocks = 0;
    int nFlippedTotal = 0;
    double totalVol = 0.0;
    double negVolTotal = 0.0;
    set<EH> cutEdges;
    set<FH> cutFaces;

    for (CH b : mc.cells())
    {
        double blockVol = 0.0;
        double negVolPre = 0.0;
        double minVolPre = 0.0;
        set<CH> tetsFlippedPre;
        determineIGMStats(b, blockVol, negVolPre, minVolPre, tetsFlippedPre);
        totalVol += blockVol;

        if (tetsFlippedPre.empty())
            continue;

        map<CH, map<VH, Vec3Q>> cell2igm;
        for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
            cell2igm[tet] = meshProps().ref<CHART_IGM>(tet);

        set<CH> flippedTetsTLC;
        double negVolTLC = DBL_MAX;
        double minVolTLC = DBL_MAX;
#ifdef C4HEX_WITH_TLC
        if (_blocksOptimizedByTLC.count(b) == 0)
        {
            _blocksOptimizedByTLC.insert(b);
            untangleByTLC(b, maxIter, secondsTimeLimit);
            determineIGMStats(b, blockVol, negVolTLC, minVolTLC, flippedTetsTLC);

            LOG(INFO) << "BLOCK " << b << ": after TLC, bad/all tets: " << flippedTetsTLC.size() << "/"
                      << mcMeshProps().ref<BLOCK_MESH_TETS>(b).size() << " ("
                      << (double)flippedTetsTLC.size() / mcMeshProps().ref<BLOCK_MESH_TETS>(b).size() * 100.0
                      << "%); neg/total volume: " << negVolTLC << "/" << blockVol << " ("
                      << negVolTLC / blockVol * 100.0 << "%)";

            if (flippedTetsTLC.size() == 0)
                continue;

            if (flippedTetsTLC.size() > tetsFlippedPre.size()
                || (flippedTetsTLC.size() == tetsFlippedPre.size() && minVolTLC < minVolPre))
            {
                for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
                    meshProps().ref<CHART_IGM>(tet) = cell2igm[tet];
                negVolTLC = negVolPre;
                minVolTLC = minVolPre;
                flippedTetsTLC = tetsFlippedPre;
            }
            else
            {
                for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
                    cell2igm[tet] = meshProps().ref<CHART_IGM>(tet);
            }
        }
#endif

        untangleByFFM(b, areaVsAngles, maxIter, secondsTimeLimit);

        double negVolFFM = 0.0;
        double minVolFFM = DBL_MAX;
        set<CH> flippedTetsFFM;
        determineIGMStats(b, blockVol, negVolFFM, minVolFFM, flippedTetsFFM);
        LOG(INFO) << "BLOCK " << b << ": after FFM, bad/all tets: " << flippedTetsFFM.size() << "/"
                  << mcMeshProps().ref<BLOCK_MESH_TETS>(b).size() << " ("
                  << (double)flippedTetsFFM.size() / mcMeshProps().ref<BLOCK_MESH_TETS>(b).size() * 100.0
                  << "%); neg/total volume: " << negVolFFM << "/" << blockVol << " (" << negVolFFM / blockVol * 100.0
                  << "%)";

        if (flippedTetsFFM.size() == 0)
            continue;

        nBadBlocks++;

        set<CH> flippedTets;
        if (minVolTLC == DBL_MAX || flippedTetsFFM.size() < flippedTetsTLC.size()
            || (flippedTetsFFM.size() == flippedTetsTLC.size() && minVolFFM > minVolTLC))
        {
            negVolTotal += negVolFFM;
            nFlippedTotal += flippedTetsFFM.size();
            flippedTets = flippedTetsFFM;
        }
        else
        {
            negVolTotal += negVolTLC;
            nFlippedTotal += flippedTetsTLC.size();
            flippedTets = flippedTetsTLC;
            for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
                meshProps().ref<CHART_IGM>(tet) = cell2igm[tet];
        }

        for (CH tet : flippedTets)
        {
            for (EH e : mesh.cell_edges(tet))
            {
                if (!meshProps().isInPatch(e)
                    && !containsMatching(mesh.edge_vertices(e),
                                         [this](const VH& v) { return !meshProps().isInPatch(v); }))
                    cutEdges.insert(e);
            }
            for (FH f : mesh.cell_faces(tet))
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
        }
    }
    if (!cutFaces.empty() || !cutEdges.empty())
    {
        for (EH e : cutEdges)
            splitHalfEdge(mesh.halfedge_handle(e, 0), *mesh.ec_iter(e), Q(0.5));
        for (FH f : cutFaces)
            if (!mesh.is_deleted(f))
                splitFace(f, {Q(1, 3), Q(1, 3), Q(1, 3)});
        reset();
        LOG(INFO) << "Found splittable edges, splitting and trying to untangle again...";
        return untangleIGM(areaVsAngles, maxIter, secondsTimeLimit);
    }

    LOG(INFO) << "After igm optimization, bad/all blocks: " << nBadBlocks << "/" << mc.n_logical_cells()
              << "; bad/all tets: " << nFlippedTotal << "/" << mesh.n_logical_cells() << " ("
              << (double)nFlippedTotal / mesh.n_logical_cells() * 100.0 << "%); neg/total volume: " << negVolTotal
              << "/" << totalVol << " (" << negVolTotal / totalVol * 100.0 << "%)";

    if (nFlippedTotal > 0)
        return NO_CONVERGENCE;
#ifndef NDEBUG
    else
    {
        for (CH tet : meshProps().mesh().cells())
            if (rationalVolumeIGM(tet) <= 0)
                throw std::logic_error("Flipped tet after successful per block optimization");
    }
#endif
    return SUCCESS;
}

void IGMUntangler::determineIGMStats(
    const CH& b, double& blockVol, double& volInverted, double& minVol, set<CH>& tetsInverted) const
{
    volInverted = 0.0;
    blockVol = 0.0;
    minVol = DBL_MAX;
    tetsInverted.clear();
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        Q vol = rationalVolumeIGM(tet);
        minVol = std::min(minVol, vol.get_d());
        blockVol += vol.get_d();
        if (vol <= 0)
        {
            tetsInverted.insert(tet);
            volInverted += -vol.get_d();
        }
    }
}

IGMUntangler::RetCode IGMUntangler::untangleByTLC(const CH& b, int maxIter, int secondsTimeLimit)
{
#ifndef C4HEX_WITH_TLC
    (void)b;
    (void)maxIter;
    (void)secondsTimeLimit;
    LOG(ERROR) << "Compiled without nlopt/TLC, unable to untangle by tlc";
    return NO_CONVERGENCE;
#else
    const TetMesh& mesh = meshProps().mesh();

    LOG(INFO) << "Fixing/optimizing by TLC block " << b;

    // import data

    int maxIdx = 0;
    map<VH, map<CH, int>> v2corner2idx;
    vector<pair<VH, vector<CH>>> idx2v;
    vector<CH> idx2tet;
    vector<vector<unsigned>> F(mcMeshProps().ref<BLOCK_MESH_TETS>(b).size());
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        for (VH v : mesh.tet_vertices(tet))
        {
            if (v2corner2idx[v].count(tet) == 0)
            {
                set<CH> tetsLocal;
                forVertexNeighbourTetsInBlock(v,
                                              tet,
                                              [&tetsLocal](const CH& tet2)
                                              {
                                                  tetsLocal.insert(tet2);
                                                  return false;
                                              });
                int nextIdx = maxIdx++;
                idx2v.push_back({v, {}});
                for (CH tet2 : tetsLocal)
                {
                    v2corner2idx[v][tet2] = nextIdx;
                    idx2v.back().second.push_back(tet2);
                }
            }
            F[idx2tet.size()].push_back(v2corner2idx[v].at(tet));
        }
        idx2tet.push_back(tet);
    }

    vector<vector<double>> XYZs(idx2v.size());
    vector<vector<double>> initUVWs(idx2v.size());
    vector<unsigned> fixedVsIndices;
    for (int idx = 0; idx < (int)idx2v.size(); idx++)
    {
        VH v = idx2v[idx].first;
        auto& tets = idx2v[idx].second;
        if (meshProps().isInPatch(v))
            fixedVsIndices.push_back(idx);
        Vec3d pos = mesh.vertex(v);
        XYZs[idx] = {pos[0], pos[1], pos[2]};
        Vec3d uvw = Vec3Q2d(meshProps().ref<CHART_IGM>(tets.front()).at(v));
        initUVWs[idx] = {uvw[0], uvw[1], uvw[2]};
    }

    // import options
    tlc::LiftedData data(XYZs, initUVWs, F, fixedVsIndices, "Tutte", 1e-6, -1);
    data.stopCode = "all_good";

    unsigned nv = XYZs.size();
    unsigned nfree = nv - fixedVsIndices.size();
    unsigned embed_dim = initUVWs[0].size();
    unsigned problem_dim = nfree * embed_dim;

    // set algorithm
    nlopt::algorithm algo = nlopt::LD_LBFGS;
    nlopt::opt opt(algo, problem_dim);

    // set stop criteria
    opt.set_ftol_abs(1e-15);
    opt.set_ftol_rel(1e-15);
    opt.set_xtol_abs(1e-15);
    opt.set_xtol_rel(1e-15);
    opt.set_maxeval(maxIter * 100);
    opt.set_maxtime(secondsTimeLimit);
    opt.set_min_objective(tlc::lifted_func, &data);

    vector<double> x = data.x0;
    double minf;

    nlopt::result result;
    // untangle
    try
    {
        result = opt.optimize(x, minf);
    }
    catch (std::exception& e)
    {
        std::cout << "nlopt failed: " << e.what() << std::endl;
        result = nlopt::FAILURE;
    }

    std::cout << "result: ";
    switch (result)
    {
    case nlopt::SUCCESS:
        std::cout << "SUCCESS" << std::endl;
        break;
    case nlopt::STOPVAL_REACHED:
        std::cout << "STOPVAL_REACHED" << std::endl;
        break;
    case nlopt::FTOL_REACHED:
        std::cout << "FTOL_REACHED" << std::endl;
        break;
    case nlopt::XTOL_REACHED:
        std::cout << "XTOL_REACHED" << std::endl;
        break;
    case nlopt::MAXEVAL_REACHED:
        std::cout << "MAXEVAL_REACHED" << std::endl;
        break;
    case nlopt::MAXTIME_REACHED:
        std::cout << "MAXTIME_REACHED" << std::endl;
        break;
    default:
        std::cout << "unexpected return code!" << std::endl;
        break;
    }

    for (int i = 0; i < (int)idx2tet.size(); i++)
    {
        CH tet = idx2tet[i];
        auto js = F[i];
        for (auto j : js)
        {
            VH v = idx2v[j].first;
            if (meshProps().isInPatch(v))
                continue;
            auto uvw = data.V[j];
            meshProps().ref<CHART_IGM>(tet).at(v) = Vec3Q(uvw[0], uvw[1], uvw[2]);
        }
    }

    return SUCCESS;
#endif
}

void IGMUntangler::reset()
{
    _blocksOptimizedByTLC.clear();
}

IGMUntangler::RetCode IGMUntangler::untangleByFFM(const CH& b, double areaVsAngles, int maxIter, int secondsTimeLimit)
{
    const TetMesh& mesh = meshProps().mesh();

    // Gather block vertices
    map<VH, map<CH, int>> v2corner2idx;
    Eigen::VectorXd lowerBounds, upperBounds;
    int nInteriorVs = 0;
    Eigen::Matrix4Xi tetVtxIndices = fixBlockVertices(b, v2corner2idx, lowerBounds, upperBounds, nInteriorVs);

    if (v2corner2idx.size() == 0)
        return SUCCESS;

    LOG(INFO) << "Fixing/optimizing by FFM block " << b;

    int nVertices = lowerBounds.size() / 3;
    int nTets = mcMeshProps().ref<BLOCK_MESH_TETS>(b).size();

    Eigen::VectorXd UVWflat(Eigen::VectorXd::Zero(3 * nVertices));
    Eigen::VectorXd XYZflat(Eigen::VectorXd::Zero(3 * nVertices));
    Eigen::VectorXd tetWeights(Eigen::VectorXd::Ones(nTets));

    for (const auto& kv : v2corner2idx)
    {
        VH v = kv.first;
        map<int, vector<CH>> indices;
        for (auto& kv2 : kv.second)
            indices[kv2.second * 3].push_back(kv2.first);

        for (auto& kv2 : indices)
        {
            int idx = kv2.first;
            auto& tets = kv2.second;
            Vec3d igm = Vec3Q2d(meshProps().ref<CHART_IGM>(tets.front()).at(v));
            Vec3d xyz = mesh.vertex(v);
            for (int coord = 0; coord < 3; coord++)
            {
                UVWflat(idx + coord) = igm[coord];
                XYZflat(idx + coord) = xyz[coord];
            }
        }
    }

    Eigen::Matrix4Xd tetZ(Eigen::Matrix4Xd::Zero(4, 3 * nTets));
    Eigen::VectorXd tetVol(Eigen::VectorXd::Zero(nTets));
    auto ret = precompute(tetVtxIndices, XYZflat, UVWflat, tetZ, tetVol);
    if (ret != SUCCESS)
        return ret;

    Eigen::Matrix3Xd tetJ(Eigen::Matrix3Xd::Zero(3, 3 * nTets));
    Eigen::VectorXd tetDetJ(Eigen::VectorXd::Zero(nTets));
    double minDetJ0 = DBL_MAX;
    int nFlipped = INT_MAX;
    ret = updateJ(tetVtxIndices, UVWflat, tetZ, tetJ, tetDetJ, minDetJ0, nFlipped);
    if (ret != SUCCESS)
        return ret;

    double e0 = 10.0;
    double F0 = 0.0;

    ret = calcF(tetVtxIndices, tetVol, tetWeights, tetJ, tetDetJ, e0, areaVsAngles, F0);
    if (ret != SUCCESS)
        return ret;

    LOG(INFO) << "Initial F " << F0 << ", MinDetJ " << minDetJ0 << ", nFlipped " << nFlipped;

    double minDetJ = minDetJ0;
    double e = e0;
    double F = F0;
    double lastF = DBL_MAX;
    double lastE = DBL_MAX;

    int bestNFlipped = nFlipped;
    double bestMinDetJ = minDetJ0;
    Eigen::VectorXd bestUVWflat = UVWflat;
    Eigen::VectorXd bestTetDetJ = tetDetJ;

    Eigen::VectorXd gradFflat(3 * nInteriorVs);
    Eigen::SparseMatrix<double> hessF(3 * nInteriorVs, 3 * nInteriorVs);

    int iter = 0;
    int iterNoImprovement = 0;
    int maxIterNoImprovement = 100;
    int mode = 0;
    // modes:
    // 0 : lbfgs
    // 1 : gradient descent
    // 2 : newton
    // 3 : lbfgs with many inner iterations

    auto start_time = std::chrono::high_resolution_clock::now();
    while ((minDetJ <= 1e-6 || std::abs(lastF - F) / F > 1e-5) && iter++ < maxIter
           && iterNoImprovement < maxIterNoImprovement
           && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time)
                      .count()
                  < secondsTimeLimit)
    {
        bool numericalIssue = false;

        {
            if (iterNoImprovement >= 0 && iterNoImprovement <= 0.25 * maxIterNoImprovement)
                mode = 2;
            else if (iterNoImprovement > 0.25 * maxIterNoImprovement
                     && iterNoImprovement <= 0.45 * maxIterNoImprovement)
                mode = 2;
            else if (iterNoImprovement > 0.45 * maxIterNoImprovement && iterNoImprovement <= 0.5 * maxIterNoImprovement)
                mode = 2;
            else if (iterNoImprovement > 0.5 * maxIterNoImprovement && iterNoImprovement <= 0.75 * maxIterNoImprovement)
                mode = 2;
            else if (iterNoImprovement > 0.75 * maxIterNoImprovement
                     && iterNoImprovement <= 0.95 * maxIterNoImprovement)
                mode = 3;
            else
                mode = 3;
        }

        if (mode != 0 && !numericalIssue)
        {
            ret = updateGradF(
                tetVtxIndices, tetVol, tetWeights, tetZ, tetJ, tetDetJ, nInteriorVs, e, areaVsAngles, gradFflat);
            if (ret != SUCCESS)
                numericalIssue = true;
        }

        if (mode == 2 && !numericalIssue)
        {
            ret = updateHessF(
                tetVtxIndices, tetVol, tetWeights, tetZ, tetJ, tetDetJ, nInteriorVs, e, areaVsAngles, hessF);
            if (ret != SUCCESS)
                mode = 3;
        }

        if (!numericalIssue)
        {

            ret = updateUVW(tetVtxIndices,
                            tetVol,
                            tetWeights,
                            gradFflat,
                            hessF,
                            tetZ,
                            nInteriorVs,
                            lowerBounds,
                            upperBounds,
                            e,
                            areaVsAngles,
                            tetJ,
                            tetDetJ,
                            UVWflat,
                            iterNoImprovement,
                            mode);

            if (ret != SUCCESS || iterNoImprovement >= (int)(0.5 * maxIterNoImprovement))
                numericalIssue = true;
        }

        ret = updateJ(tetVtxIndices, UVWflat, tetZ, tetJ, tetDetJ, minDetJ, nFlipped);
        if (ret != SUCCESS)
            break;

        if (nFlipped < bestNFlipped || (nFlipped == bestNFlipped && minDetJ > bestMinDetJ))
        {
            bestNFlipped = nFlipped;
            bestMinDetJ = minDetJ;
            bestUVWflat = UVWflat;
            bestTetDetJ = tetDetJ;
        }

        lastF = F;
        ret = calcF(tetVtxIndices, tetVol, tetWeights, tetJ, tetDetJ, e, areaVsAngles, F);
        if (ret != SUCCESS)
            numericalIssue = true;

        lastE = e;
        if (numericalIssue)
        {
            // RESET
            DLOG(WARNING) << "Increasing untangling eps due to numerical issue";
            e = std::min(1.0, e * 1e3);
            lastF = F;
            UVWflat = bestUVWflat;
            updateJ(tetVtxIndices, UVWflat, tetZ, tetJ, tetDetJ, minDetJ, nFlipped);
            calcF(tetVtxIndices, tetVol, tetWeights, tetJ, tetDetJ, e, areaVsAngles, F);
            if (ret != SUCCESS)
                break;
        }
        else
        {
            ret = updateE(minDetJ, lastF, F, lastE, e);
            if (ret != SUCCESS)
                break;
        }
        if (iter % std::max(1, (int)(0.2 * maxIter)) == 0)
            LOG(INFO) << "Iter " << iter << " e " << e << " F " << F << ", MinDetJ " << minDetJ << ", nFlipped "
                      << nFlipped;
        for (int tetIdx = 0; tetIdx < tetVtxIndices.cols(); tetIdx++)
        {
            if (iterNoImprovement > 0.25 * maxIterNoImprovement && tetDetJ(tetIdx) < 0)
                tetWeights(tetIdx) *= 1.1;
            else if (tetWeights(tetIdx) > 1.0)
                tetWeights(tetIdx) -= (tetWeights(tetIdx) - 1) / 2.0;
        }
    }

    LOG(INFO) << "Final MinDetJ " << bestMinDetJ << ", nFlipped " << bestNFlipped;
    applyUntangling(v2corner2idx, nInteriorVs, bestUVWflat);

    int tetIdx = 0;
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        bool valid = rationalVolumeIGM(tet) > 0;
        if (valid != (bestTetDetJ(tetIdx) > 0))
        {
            LOG(INFO) << "Tet " << tet << " volume is " << rationalVolumeIGM(tet).get_d() << ", but detJ is "
                      << bestTetDetJ(tetIdx);
        }
        tetIdx++;
    }

    if (bestMinDetJ > 0)
        return SUCCESS;
    else
        return NO_CONVERGENCE;
}

Eigen::Matrix4Xi IGMUntangler::fixBlockVertices(const CH& b,
                                                map<VH, map<CH, int>>& v2corner2idx,
                                                Eigen::VectorXd& lowerBounds,
                                                Eigen::VectorXd& upperBounds,
                                                int& nInteriorVertices)
{
    DLOG(INFO) << "Fixing block vertices";
    v2corner2idx.clear();
    Eigen::Matrix4Xi tetVtxIndices(Eigen::Matrix4Xi::Zero(4, mcMeshProps().ref<BLOCK_MESH_TETS>(b).size()));
    const TetMesh& mesh = meshProps().mesh();
    Vec3d minIGMblock(DBL_MAX, DBL_MAX, DBL_MAX);
    Vec3d maxIGMblock(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    int maxIdx = 0;
    nInteriorVertices = 0;

    // First non-boundary
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        for (VH v : mesh.tet_vertices(tet))
        {
            if (v2corner2idx[v].count(tet) == 0)
            {
                if (meshProps().isInPatch(v))
                    continue;

                set<CH> tetsLocal;
                forVertexNeighbourTetsInBlock(v,
                                              tet,
                                              [&tetsLocal](const CH& tet2)
                                              {
                                                  tetsLocal.insert(tet2);
                                                  return false;
                                              });
                int nextIdx = maxIdx++;
                for (CH tet2 : tetsLocal)
                    v2corner2idx[v][tet2] = nextIdx;
                Vec3d igm = Vec3Q2d(meshProps().ref<CHART_IGM>(tet).at(v));
                for (int coord = 0; coord < 3; coord++)
                {
                    minIGMblock[coord] = std::min(minIGMblock[coord], igm[coord]);
                    maxIGMblock[coord] = std::max(maxIGMblock[coord], igm[coord]);
                }
            }
        }
    }
    nInteriorVertices = maxIdx;

    // Then boundary
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        for (VH v : mesh.tet_vertices(tet))
        {
            if (v2corner2idx[v].count(tet) == 0)
            {
                if (!meshProps().isInPatch(v))
                    continue;
                set<CH> tetsLocal;
                forVertexNeighbourTetsInBlock(v,
                                              tet,
                                              [&tetsLocal](const CH& tet2)
                                              {
                                                  tetsLocal.insert(tet2);
                                                  return false;
                                              });
                int nextIdx = maxIdx++;
                for (CH tet2 : tetsLocal)
                    v2corner2idx[v][tet2] = nextIdx;
                Vec3d igm = Vec3Q2d(meshProps().ref<CHART_IGM>(tet).at(v));
                for (int coord = 0; coord < 3; coord++)
                {
                    minIGMblock[coord] = std::min(minIGMblock[coord], igm[coord]);
                    maxIGMblock[coord] = std::max(maxIGMblock[coord], igm[coord]);
                }
            }
        }
    }

    lowerBounds = Eigen::VectorXd(3 * maxIdx);
    lowerBounds.setConstant(-DBL_MAX);
    upperBounds = Eigen::VectorXd(3 * maxIdx);
    upperBounds.setConstant(DBL_MAX);

    for (int i = 0; i < lowerBounds.size(); i++)
    {
        lowerBounds(i) = minIGMblock[i % 3] + 1e-6;
        upperBounds(i) = maxIGMblock[i % 3] - 1e-6;
    }

    int tetIdx = 0;
    for (CH tet : mcMeshProps().ref<BLOCK_MESH_TETS>(b))
    {
        int vtxNum = 0;
        for (VH v : mesh.tet_vertices(tet))
        {
            tetVtxIndices(vtxNum, tetIdx) = v2corner2idx.at(v).at(tet);
            vtxNum++;
        }
        tetIdx++;
    }

    return tetVtxIndices;
}

IGMUntangler::RetCode IGMUntangler::precompute(const Eigen::Matrix4Xi& tetVtxIndices,
                                               const Eigen::VectorXd XYZflat,
                                               const Eigen::VectorXd UVWflat,
                                               Eigen::Matrix4Xd& tetZ,
                                               Eigen::VectorXd& tetVol)
{
    Eigen::Vector3d minXYZ(DBL_MAX, DBL_MAX, DBL_MAX);
    Eigen::Vector3d maxXYZ(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    Eigen::Vector3d minUVW(DBL_MAX, DBL_MAX, DBL_MAX);
    Eigen::Vector3d maxUVW(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    for (int i = 0; i < XYZflat.size(); i += 3)
    {
        Eigen::Vector3d point = XYZflat.segment<3>(i);
        Eigen::Vector3d param = UVWflat.segment<3>(i);
        for (int coord = 0; coord < 3; coord++)
        {
            minXYZ(coord) = std::min(point(coord), minXYZ(coord));
            maxXYZ(coord) = std::max(point(coord), maxXYZ(coord));
            minUVW(coord) = std::min(param(coord), minUVW(coord));
            maxUVW(coord) = std::max(param(coord), maxUVW(coord));
        }
    }

    double xyzDiag = (maxXYZ - minXYZ).norm();
    double uvwDiag = (maxUVW - minUVW).norm();
    double scaling = uvwDiag / xyzDiag * 3.0;

    for (int tetIdx = 0; tetIdx < tetVtxIndices.cols(); tetIdx++)
    {
        Eigen::Matrix3d S;
        Eigen::Vector3d xyz0((XYZflat.segment<3>(3 * tetVtxIndices(0, tetIdx)) - minXYZ) * scaling);
        for (int vtxNum = 1; vtxNum < 4; vtxNum++)
        {
            Eigen::Vector3d xyz((XYZflat.segment<3>(3 * tetVtxIndices(vtxNum, tetIdx)) - minXYZ) * scaling);
            S.col(vtxNum - 1) = xyz - xyz0;
        }
        tetVol(tetIdx) = S.determinant() / 6.;
        double clampVol = 1e-6;
        if (tetVol(tetIdx) <= clampVol)
        {
            DLOG(WARNING) << "Cant use XYZ for target shape of tet " << tetIdx << " because volume " << tetVol(tetIdx)
                          << " is too small, choosing uniform target shape as fallback";
            double factor = std::pow(clampVol, 1.0 / 3.0);
            S = Eigen::Matrix3d::Zero();
            S << 1, 0, 1, 1, 1, 0, 0, 1, 1;
            S *= 1.5 * factor;
            tetVol(tetIdx) = S.determinant() / 6.;
        }
        Eigen::Matrix3d Sinv = S.inverse();
        tetZ.block<3, 3>(1, 3 * tetIdx) = Sinv;
        tetZ.block<1, 3>(0, 3 * tetIdx) = Eigen::RowVector3d(-1, -1, -1) * Sinv;
    }

    return SUCCESS;
}

IGMUntangler::RetCode IGMUntangler::updateJ(const Eigen::Matrix4Xi& tetVtxIndices,
                                            const Eigen::VectorXd UVWflat,
                                            const Eigen::Matrix4Xd& tetZ,
                                            Eigen::Matrix3Xd& tetJ,
                                            Eigen::VectorXd& tetDetJ,
                                            double& minDetJ,
                                            int& nFlipped)
{
    nFlipped = 0;
    // DLOG(INFO) << "Updating J";
    minDetJ = DBL_MAX;
    for (int tetIdx = 0; tetIdx < tetVtxIndices.cols(); tetIdx++)
    {
        Eigen::Matrix<double, 3, 4> UVWmat;
        for (int vtxNum = 0; vtxNum < 4; vtxNum++)
            UVWmat.col(vtxNum) = UVWflat.segment<3>(3 * tetVtxIndices(vtxNum, tetIdx));
        tetJ.block<3, 3>(0, 3 * tetIdx) = UVWmat * tetZ.block<4, 3>(0, 3 * tetIdx);
        Eigen::Matrix3d J = tetJ.block<3, 3>(0, 3 * tetIdx);
        double det = J.determinant();
        if (det < 1e-4)
        {
            Vec3Q v1(J(0, 0), J(1, 0), J(2, 0)), v2(J(0, 1), J(1, 1), J(2, 1)), v3(J(0, 2), J(1, 2), J(2, 2));
            det = dot(v1, cross(v2, v3)).get_d();
            if (det <= 0)
                nFlipped++;
        }
        tetDetJ(tetIdx) = det;
        if (!std::isfinite(tetDetJ(tetIdx)))
        {
            DLOG(WARNING) << "Non-finite det J encountered, probably due to numerical issues";
            return NUMERICAL_ISSUE;
        }
        minDetJ = std::min(minDetJ, tetDetJ(tetIdx));
    }
    // DLOG(INFO) << "New min det J: " << minDetJ;
    return SUCCESS;
}

IGMUntangler::RetCode IGMUntangler::calcF(const Eigen::Matrix4Xi& tetVtxIndices,
                                          const Eigen::VectorXd& tetVol,
                                          const Eigen::VectorXd& tetWeights,
                                          const Eigen::Matrix3Xd& tetJ,
                                          const Eigen::VectorXd& tetDetJ,
                                          double e,
                                          double areaVsAngles,
                                          double& F)
{
    // DLOG(INFO) << "Updating F: ";
    F = 0.;
    for (int tetIdx = 0; tetIdx < tetVtxIndices.cols(); tetIdx++)
    {
        double detJ = tetDetJ(tetIdx);
        double chi = calcChi(e, detJ);
        double f = 0.0;
        for (int i = 0; i < 3; i++)
            f += tetJ.col(3 * tetIdx + i).squaredNorm();
        f /= std::pow(chi, 2.0 / 3.0);
        double g = (detJ * detJ + 1) / chi;
        F += ((1 - areaVsAngles) * f + areaVsAngles * g) * tetVol(tetIdx) * tetWeights(tetIdx);
    }
    if (!std::isfinite(F))
    {
        DLOG(WARNING) << "Non-finite F encountered, probably due to numerical issues";
        return NUMERICAL_ISSUE;
    }
    // DLOG(INFO) << "New F: " << F;
    return SUCCESS;
}

IGMUntangler::RetCode IGMUntangler::updateGradF(const Eigen::Matrix4Xi& tetVtxIndices,
                                                const Eigen::VectorXd& tetVol,
                                                const Eigen::VectorXd& tetWeights,
                                                const Eigen::Matrix4Xd& tetZ,
                                                const Eigen::Matrix3Xd& tetJ,
                                                const Eigen::VectorXd& tetDetJ,
                                                const int nInteriorVs,
                                                double e,
                                                double areaVsAngles,
                                                Eigen::VectorXd& gradFflat)
{
    // DLOG(INFO) << "Updating gradient";
    gradFflat.setZero();
    for (int tetIdx = 0; tetIdx < tetVtxIndices.cols(); tetIdx++)
    {
        Eigen::Matrix3d J = tetJ.block<3, 3>(0, 3 * tetIdx);
        Eigen::Matrix3d B;
        for (int i = 0; i < 3; i++)
            B.col(i) = J.col((i + 1) % 3).cross(J.col((i + 2) % 3));
        double detJ = tetDetJ(tetIdx);
        double chi = calcChi(e, detJ);
        double chi23 = std::pow(chi, 2.0 / 3.0);
        double dchi = calcChiDeriv(e, detJ);
        double f = 0.0;
        for (int i = 0; i < 3; i++)
            f += tetJ.col(3 * tetIdx + i).squaredNorm();
        f /= chi23;
        double g = (detJ * detJ + 1) / chi;
        Eigen::Matrix3d dfda = 2 / chi23 * J - (2 * f * dchi / (3 * chi)) * B;
        Eigen::Matrix3d dgda = (2 * detJ - g * dchi) / chi * B;
        Eigen::Matrix<double, 4, 3> Z = tetZ.block<4, 3>(0, 3 * tetIdx);

        for (int vtxNum = 0; vtxNum < 4; vtxNum++)
        {
            int vtxIdx = tetVtxIndices(vtxNum, tetIdx);
            if (vtxIdx >= nInteriorVs)
                continue;

            Eigen::Vector3d sum(0, 0, 0);
            for (int i = 0; i < 3; i++)
                sum += Z(vtxNum, i) * ((1 - areaVsAngles) * dfda.col(i) + areaVsAngles * dgda.col(i));

            gradFflat.segment<3>(3 * vtxIdx) += tetVol(tetIdx) * tetWeights(tetIdx) * sum;

            if (!std::isfinite(sum(0)) || !std::isfinite(sum(1)) || !std::isfinite(sum(2)))
            {
                DLOG(WARNING) << "Non-finite gradF encountered, probably due to numerical issues";
                return NUMERICAL_ISSUE;
            }
        }
    }
    // DLOG(INFO) << "Updated gradient, gradient squared length " << gradFflat.squaredNorm();
    return SUCCESS;
}

IGMUntangler::RetCode IGMUntangler::updateHessF(const Eigen::Matrix4Xi& tetVtxIndices,
                                                const Eigen::VectorXd& tetVol,
                                                const Eigen::VectorXd& tetWeights,
                                                const Eigen::Matrix4Xd& tetZ,
                                                const Eigen::Matrix3Xd& tetJ,
                                                const Eigen::VectorXd& tetDetJ,
                                                const int nInteriorVs,
                                                double e,
                                                double areaVsAngles,
                                                Eigen::SparseMatrix<double>& hessF)
{
    hessF = Eigen::SparseMatrix<double>(hessF.rows(), hessF.cols());
    vector<Eigen::Triplet<double>> triplets;
    for (int tetIdx = 0; tetIdx < tetVtxIndices.cols(); tetIdx++)
    {
        Eigen::Matrix3d J = tetJ.block<3, 3>(0, 3 * tetIdx);
        Eigen::Matrix3d B;
        for (int i = 0; i < 3; i++)
            B.col(i) = J.col((i + 1) % 3).cross(J.col((i + 2) % 3));
        double detJ = tetDetJ(tetIdx);
        double chi = calcChi(e, detJ);
        double chi23 = std::pow(chi, 2.0 / 3.0);
        double dchi = calcChiDeriv(e, detJ);
        double dchi2 = dchi * dchi;
        double chi2 = chi * chi;

        Eigen::Matrix<double, 4, 3> Z = tetZ.block<4, 3>(0, 3 * tetIdx);

        Eigen::Matrix<double, 9, 10> bStretch(Eigen::Matrix<double, 9, 10>::Zero(9, 10));
        for (int ii = 0; ii < 9; ii++)
            bStretch(ii, ii) = 1.0;

        for (int i = 0; i < 3; i++)
            bStretch.block<3, 1>(3 * i, 9) = B.col(i);

        Eigen::Matrix<double, 10, 10> ddInner(Eigen::Matrix<double, 10, 10>::Zero(10, 10));
        for (int ii = 0; ii < 9; ii++)
            ddInner(ii, ii) = (1. - areaVsAngles) * 2. / chi23;

        ddInner(9, 9)
            = (1. - areaVsAngles) * 2. / 3. * (1 + 2. / 3.) * J.squaredNorm() * (dchi2 / (chi2 * chi23))
              + areaVsAngles * (2. / chi - 4 * detJ * dchi / (chi2) + 2 * (1 + detJ * detJ) * (dchi2 / (chi2 * chi)));

        double factor = (1. - areaVsAngles) * -4. / 3. * dchi / (chi * chi23);
        for (int i = 0; i < 9; i++)
            ddInner(i, 9) = ddInner(9, i) = factor * J(i % 3, i / 3);

        Eigen::Matrix<double, 9, 9> Mplus = bStretch * ddInner * bStretch.transpose();
        for (int vtxNum = 0; vtxNum < 4; vtxNum++)
        {
            int vtxIdx = tetVtxIndices(vtxNum, tetIdx);
            if (vtxIdx >= nInteriorVs)
                continue;

            for (int vtxNum2 = 0; vtxNum2 < 4; vtxNum2++)
            {
                int vtxIdx2 = tetVtxIndices(vtxNum2, tetIdx);
                if (vtxIdx2 >= nInteriorVs)
                    continue;
                Eigen::Matrix3d hessSum(Eigen::Matrix3d::Zero(3, 3));
                for (int m = 0; m < 3; m++)
                    for (int l = 0; l < 3; l++)
                        hessSum += tetVol(tetIdx) * tetWeights(tetIdx) * Z(vtxNum2, m) * Mplus.block<3, 3>(3 * m, 3 * l)
                                   * Z(vtxNum, l);

                if (!hessSum.array().isFinite().all())
                {
                    DLOG(WARNING) << "Non-finite hessF encountered, probably due to numerical issues";
                    return NUMERICAL_ISSUE;
                }
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        triplets.emplace_back(Eigen::Triplet<double>(3 * vtxIdx2 + i, 3 * vtxIdx + j, hessSum(i, j)));
            }
        }
    }
    hessF.setFromTriplets(triplets.begin(), triplets.end());

    return SUCCESS;
}

IGMUntangler::RetCode IGMUntangler::updateUVW(const Eigen::Matrix4Xi& tetVtxIndices,
                                              const Eigen::VectorXd& tetVol,
                                              const Eigen::VectorXd& tetWeights,
                                              const Eigen::VectorXd& gradF,
                                              const Eigen::SparseMatrix<double>& hessF,
                                              const Eigen::Matrix4Xd& tetZ,
                                              const int nInteriorVs,
                                              const Eigen::VectorXd& lowerBounds,
                                              const Eigen::VectorXd& upperBounds,
                                              double e,
                                              double areaVsAngles,
                                              Eigen::Matrix3Xd& tetJ,
                                              Eigen::VectorXd& tetDetJ,
                                              Eigen::VectorXd& UVWflat,
                                              int& iterNoImprovement,
                                              int mode)
{
    double startMinDetJ = DBL_MAX;
    double startF = DBL_MAX;
    int nFlipped = INT_MAX;
    if (updateJ(tetVtxIndices, UVWflat, tetZ, tetJ, tetDetJ, startMinDetJ, nFlipped) != SUCCESS
        || calcF(tetVtxIndices, tetVol, tetWeights, tetJ, tetDetJ, e, areaVsAngles, startF) != SUCCESS)
        return NUMERICAL_ISSUE;

    double bestF = startF;
    Eigen::VectorXd bestUVWflat;

    if (mode == 0 || mode == 3)
    {
        (void)gradF;
        (void)hessF;
        LBFGSpp::LBFGSBParam<double> param;
        param.epsilon = 1e-6;
        param.max_iterations = mode == 0 ? 10 : 100;
        param.max_linesearch = 30;
        LBFGSpp::LBFGSBSolver<double> solver(param);

        auto tetVolCopy = tetVol;
        tetVolCopy.array() *= tetWeights.array();
        FoldoverEnergy foldEnergy(areaVsAngles, e, tetVtxIndices, tetZ, tetVolCopy, nInteriorVs, tetWeights);

        bestUVWflat = UVWflat;

        try
        {
            solver.minimize(foldEnergy, bestUVWflat, bestF, lowerBounds, upperBounds);
        }
        catch (std::exception& exc)
        {
            bestF = DBL_MAX;
        }
    }
    else
    {
        Eigen::VectorXd descentDirShort(Eigen::VectorXd::Zero(gradF.size()));
        if (mode == 1)
        {
            (void)hessF;
            descentDirShort = -gradF;
        }
        else
        {
            // DLOG(INFO) << "Newton to find optimal UVWflat";
            if (nInteriorVs > 0)
            {
                double maxAbsVal = -DBL_MAX;
                for (int k = 0; k < hessF.outerSize(); ++k)
                    for (Eigen::SparseMatrix<double>::InnerIterator it(hessF, k); it; ++it)
                        maxAbsVal = std::max(std::abs(it.value()), maxAbsVal);
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(hessF / maxAbsVal);
                if (chol.info() != Eigen::Success)
                {
                    DLOG(INFO) << "Cholesky solving failed";
                    return updateUVW(tetVtxIndices,
                                     tetVol,
                                     tetWeights,
                                     gradF,
                                     hessF,
                                     tetZ,
                                     nInteriorVs,
                                     lowerBounds,
                                     upperBounds,
                                     e,
                                     areaVsAngles,
                                     tetJ,
                                     tetDetJ,
                                     UVWflat,
                                     iterNoImprovement,
                                     3);
                }
                else
                {
                    descentDirShort = -chol.solve(gradF / gradF.cwiseAbs().maxCoeff());
                }
            }
        }

        descentDirShort /= descentDirShort.cwiseAbs().maxCoeff();

        Eigen::VectorXd descentDir(Eigen::VectorXd::Zero(UVWflat.size()));
        Eigen::VectorXd nextGrad(Eigen::VectorXd::Zero(UVWflat.size()));

        double maxStep = DBL_MAX;
        for (int i = 0; i < descentDirShort.size(); i++)
        {
            descentDir(i) = descentDirShort(i);
            nextGrad(i) = gradF(i);
            if (descentDir(i) > 0.0)
            {
                double stepLimit = (upperBounds(i) - UVWflat(i)) / descentDir(i);
                if (stepLimit < 1e-5)
                    descentDir(i) = (upperBounds(i) - UVWflat(i)) / maxStep;
                else
                    maxStep = std::min(maxStep, stepLimit);
            }
            else if (descentDir(i) < 0.0)
            {
                double stepLimit = (lowerBounds(i) - UVWflat(i)) / descentDir(i);
                if (stepLimit < 1e-5)
                    descentDir(i) = (lowerBounds(i) - UVWflat(i)) / maxStep;
                else
                    maxStep = std::min(maxStep, stepLimit);
            }
        }
        maxStep *= 0.99;

        auto tetVolCopy = tetVol;
        tetVolCopy.array() *= tetWeights.array();
        FoldoverEnergy foldEnergy(areaVsAngles, e, tetVtxIndices, tetZ, tetVolCopy, nInteriorVs, tetWeights);
        LBFGSpp::LBFGSBParam<double> param;
        param.epsilon = 1e-6;
        param.max_linesearch = 30;
        param.max_step = maxStep;
        bestUVWflat = UVWflat;
        double dg = descentDir.dot(nextGrad);

        double step = std::min(1.0, maxStep);
        try
        {
            LBFGSpp::LineSearchMoreThuente<double>::LineSearch(
                foldEnergy, param, UVWflat, descentDir, maxStep, step, bestF, nextGrad, dg, bestUVWflat);
        }
        catch (std::exception& exc)
        {
            bestF = DBL_MAX;
            DLOG(INFO) << "LineSearch failed: " << exc.what();
        }
    }

    if (bestF >= startF)
    {
        iterNoImprovement++;
    }
    else
    {
        iterNoImprovement = 0;
        for (int i = 0; i < 3 * nInteriorVs; i++)
            if (i < 3 * nInteriorVs)
                UVWflat(i) = std::min(upperBounds(i), std::max(lowerBounds(i), bestUVWflat(i)));
    }

    return SUCCESS;
}

IGMUntangler::RetCode IGMUntangler::updateE(double minDetJ, double lastF, double F, double lastE, double& e)
{
    double o = std::max(0.1, 1 - F / lastF);
    if (minDetJ < 1e-6)
    {
        double sqrteeDD = std::sqrt(lastE * lastE + minDetJ * minDetJ);
        e = (1 - o * sqrteeDD / (std::abs(minDetJ) + sqrteeDD)) * lastE;
    }
    else
        e = (1 - o) * lastE;

    if (!std::isfinite(e))
    {
        LOG(ERROR) << "Non-finite e encountered, probably due to numerical issues";
        return NUMERICAL_ISSUE;
    }
    return SUCCESS;
}

void IGMUntangler::applyUntangling(const map<VH, map<CH, int>>& v2corner2idx,
                                   const int nInteriorVs,
                                   const Eigen::VectorXd& UVWflat)
{
    const TetMesh& mesh = meshProps().mesh();
    for (const auto& kv : v2corner2idx)
    {
        VH v = kv.first;

        map<int, vector<CH>> indices;
        for (auto& kv2 : kv.second)
            indices[kv2.second].push_back(kv2.first);

        for (auto& kv2 : indices)
        {
            int idx = kv2.first;
            auto tets = kv2.second;
            Vec3Q igm = meshProps().ref<CHART_IGM>(tets.front()).at(v);
            if (idx >= nInteriorVs)
                continue;

            for (int coord = 0; coord < 3; coord++)
                igm[coord] = UVWflat(3 * idx + coord);

            CH tetRef = tets.front();
            auto tet2trans = determineTransitionsAroundVertex<TRANSITION_IGM>(v, tetRef);

            for (CH tet : mesh.vertex_cells(v))
                meshProps().ref<CHART_IGM>(tet).at(v) = tet2trans.at(tet).apply(igm);
        }
    }
}

} // namespace c4hex
