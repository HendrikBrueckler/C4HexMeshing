#include "C4Hex/Interface/IGMGenerator.hpp"

#include "C4Hex/Algorithm/IGMInitializer.hpp"
#include "C4Hex/Algorithm/IGMUntangler.hpp"

#include <MC3D/Algorithm/TetRemesher.hpp>

namespace c4hex
{
using namespace mc3d;

IGMGenerator::IGMGenerator(TetMeshProps& meshProps)
    : TetMeshNavigator(meshProps), TetMeshManipulator(meshProps), MCMeshNavigator(meshProps)
{
}

IGMGenerator::RetCode IGMGenerator::generateBlockwiseIGM(bool simplifyBaseMesh, int maxInnerIter, int maxUntanglingIter)
{
    TetRemesher remesher(meshProps());

    if (simplifyBaseMesh)
    {
        LOG(INFO) << "Simplifying base mesh to simplify IGM computation";
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 30);
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 20);
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 10);
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 5);
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 2);
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, false, true, 5);
    }
    else
        remesher.collapseAllPossibleEdges(true, true, true, true, 5);

    IGMInitializer init(meshProps());
    auto retInit = init.initializeFromQuantization();
    if (retInit != IGMInitializer::SUCCESS && retInit != IGMInitializer::INVALID_ELEMENTS)
        return INITIALIZATION_ERROR;
    if (simplifyBaseMesh)
    {
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 20);
        // Initialize again, hoping for less refinement
        retInit = init.initializeFromQuantization();
        if (retInit != IGMInitializer::SUCCESS && retInit != IGMInitializer::INVALID_ELEMENTS)
            return INITIALIZATION_ERROR;

        if (retInit != IGMInitializer::INVALID_ELEMENTS)
            return SUCCESS;
    }

    LOG(INFO) << "After initializing the IGM via naive 3D tutte, the following blocks have inversions:";
    for (CH b : mcMeshProps().mesh().cells())
    {
        int nInvalidUVW = 0;
        for (CH tet : mcMeshProps().get<BLOCK_MESH_TETS>(b))
            if (rationalVolumeIGM(tet) <= 0)
                nInvalidUVW++;
        if (nInvalidUVW > 0)
            LOG(INFO) << "Block " << b << " nTetsInvalidIGM " << nInvalidUVW;
    }
    LOG(INFO) << "...proceeding to IGM untangling";

    IGMUntangler optimizer(meshProps());
    auto retUntangling = IGMUntangler::NO_CONVERGENCE;
    for (int i = 0; i < maxUntanglingIter && retUntangling != IGMUntangler::SUCCESS; i++)
    {
        double areaVsAngles = 0.5 + (((i + 1) % 3) - 1) * (0.49);
        retUntangling = optimizer.untangleIGM(
            areaVsAngles,
            i < 0.75 * maxUntanglingIter ? maxInnerIter
                                         : (i < 0.9 * maxUntanglingIter ? 3 * maxInnerIter : 10 * maxInnerIter));
        if (retUntangling != IGMUntangler::SUCCESS && retUntangling != IGMUntangler::NO_CONVERGENCE)
            return UNTANGLING_ERROR;
        if (retUntangling != IGMUntangler::SUCCESS)
        {
            // Mark valid blocks as excluded
            set<CH> excludedBlocks;
            for (CH b : mcMeshProps().mesh().cells())
                if (!containsMatching(meshProps().get<MC_MESH_PROPS>()->get<BLOCK_MESH_TETS>(b),
                                      [&](const CH& tet) { return rationalVolumeIGM(tet) <= 0; }))
                    excludedBlocks.insert(b);
            meshProps().allocate<TOUCHED>(true);
            // Mark vertices not incident on non-excluded block excluded
            for (auto v : meshProps().mesh().vertices())
                if (findNoneOf(meshProps().mesh().vertex_cells(v), excludedBlocks).is_valid())
                    meshProps().set<TOUCHED>(v, false);

            // Try different remeshing/collapsing variations to improve condition of untangling problem
            if (simplifyBaseMesh)
            {
                if (i % 3 == 0)
                {
                    if (i < 0.75 * maxUntanglingIter)
                        remesher.remeshToImproveAngles(
                            true, true, TetRemesher::QualityMeasure::ANGLES, 0, excludedBlocks);
                    else
                        remesher.remeshToImproveAngles(
                            true, true, TetRemesher::QualityMeasure::VL_RATIO, 0, excludedBlocks);
                }
                else if (i % 3 == 1)
                    remesher.collapseAllPossibleEdges(simplifyBaseMesh, true, true, true, 5);
                else
                    remesher.collapseAllPossibleEdges(simplifyBaseMesh, true, true, true, 0);
            }
            else
            {
                if (i % 3 == 2)
                    remesher.collapseAllPossibleEdges(true, true, true, true, 0);
                else
                    remesher.collapseAllPossibleEdges(true, true, true, true, 5);
            }

            if (i > 5)
            {
                if ((i % 2) == 0)
                    init.splitAllSufficientForInjectiveInterior();
                else if ((i % 2) == 1)
                    init.splitSomeNecessaryForInjectiveInterior();
            }
            else
            {
                if ((i % 3) == 2)
                    init.splitAllSufficientForInjectiveInterior();
                else if ((i % 3) == 1)
                    init.splitSomeNecessaryForInjectiveInterior();
                else
                {
                    if (i < 3)
                        remesher.collapseAllPossibleEdges(!simplifyBaseMesh, true, true, true, 10);
                    else
                        remesher.collapseAllPossibleEdges(!simplifyBaseMesh, true, true, false, 0);
                }
            }
            optimizer.reset();
        }
    }
    meshProps().allocate<TOUCHED>(true);
    if (simplifyBaseMesh)
        remesher.collapseAllPossibleEdges(false, true, true);
    else
        remesher.collapseAllPossibleEdges(true, true, true, true, 5);

    if (retUntangling == IGMUntangler::NO_CONVERGENCE)
    {
        LOG(INFO) << "Some inverted tets remain in IGM in the following blocks:";
        init.minimizeCutSurface();
        for (CH b : mcMeshProps().mesh().cells())
        {
            int nInvalidUVW = 0;
            for (CH tet : mcMeshProps().get<BLOCK_MESH_TETS>(b))
                if (rationalVolumeIGM(tet) <= 0)
                    nInvalidUVW++;
            if (nInvalidUVW > 0)
                LOG(INFO) << "Block " << b << " nTetsInvalidIGM " << nInvalidUVW;
        }
        return NO_CONVERGENCE;
    }

    return SUCCESS;
}

} // namespace c4hex
