#include <MC3D/Interface/MCGenerator.hpp>
#include <MC3D/Interface/Reader.hpp>
#include <MC3D/Interface/Writer.hpp>

#include <MC3D/Algorithm/MCBuilder.hpp>
#include <MC3D/Algorithm/MotorcycleSpawner.hpp>
#include <MC3D/Algorithm/MotorcycleTracer.hpp>
#include <MC3D/Algorithm/SingularityInitializer.hpp>
#include <MC3D/Algorithm/TetRemesher.hpp>

#include <QGP3D/ConstraintExtractor.hpp>
#include <QGP3D/ConstraintWriter.hpp>
#include <QGP3D/IQP/IQPQuantizer.hpp>
#include <QGP3D/ISP/ISPQuantizer.hpp>
#include <QGP3D/ObjectiveBuilder.hpp>
#include <QGP3D/StructurePreserver.hpp>

#include <C4Hex/Algorithm/MCCollapser.hpp>
#include <C4Hex/Interface/HexRemesher.hpp>
#include <C4Hex/Interface/IGMGenerator.hpp>

#include <C4Hex/Algorithm/IGMInitializer.hpp>
#include <C4Hex/Algorithm/IGMUntangler.hpp>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include <CLI/CLI.hpp>

#include <chrono>
#include <iomanip>
#include <string>

using namespace mc3d;
using namespace qgp3d;
using namespace c4hex;

#define ASSERT_SUCCESS(stage, call)                                                                                    \
    do                                                                                                                 \
    {                                                                                                                  \
        LOG(INFO) << stage << "...";                                                                                   \
        if (auto _ERROR_CODE_ = call; _ERROR_CODE_ != 0)                                                               \
        {                                                                                                              \
            LOG(ERROR) << stage << " failed with error code " << _ERROR_CODE_ << ", aborting...";                      \
            return (_ERROR_CODE_);                                                                                     \
        }                                                                                                              \
        LOG(INFO) << stage << " was successful!";                                                                      \
    } while (0)

int main(int argc, char** argv)
{
    // Manage cli options
    CLI::App app{"MC3D"};
    std::string inputFile = "";
    std::string wallsFile = "";
    bool simulateBC = false;
    bool inputHasMCwalls = false;
    bool splitSelfadjacent = false;
    bool reduceSingularWalls = false;
    bool exactOutput = false;
    bool forceSanitization = false;
    bool doCollapse = false;

    double individualArcFactor = 0.01;
    std::string constraintFile = "";
    std::string outputIGMFile = "";
    int untanglingIter = 500;

    std::string outputHexFile = "";

    bool blockStructured = false;
    bool debugCollapse = false;
    int timesMinimalHexes = 100;
    bool randomOrder = false;
    int direction = 0;

    bool optimizeBaseMesh = false;

    // 1st bit: 0 flexible / 1 fixed singularities
    // 2nd bit: 0 distortion / 1 arcdeviation objective
    // 3rd bit: 1 for singularity padding (only has effect when 1st bit is 0)
    uint8_t variant = 0;

    double nTargetHexes = 10000;

    auto optInput
        = app.add_option("--input", inputFile, "Specify the input mesh & seamless parametrization file.")->required();
    app.add_flag("--input-has-walls",
                 inputHasMCwalls,
                 "Use, if the input already contains precomputed MC walls. Only use this, if you are sure the input is "
                 "numerically sane!");
    app.add_flag("--force-sanitization",
                 forceSanitization,
                 "Whether input sanitization should be forced even for exact rational input");
    app.add_option("--output-walls",
                   wallsFile,
                   "Specify an output file to write the refined"
                   " mesh & parametrization & MC walls to.");
    app.add_flag("--output-exact, !--param-output-double",
                 exactOutput,
                 "Whether the parametrization should be output in rational numbers"
                 " (not conforming to standard .hexex format, but numerically safe!) or doubles (numerically unsafe)");
    app.add_flag("--bc, !--mc", simulateBC, "Whether BC or MC is to be computed");
    app.add_flag("--split-selfadjacent, !--allow-selfadjacent",
                 splitSelfadjacent,
                 "Whether selfadjacent blocks should be split");
    app.add_flag("--reduce-singularity-walls, !--keep-singularity-walls",
                 reduceSingularWalls,
                 "Whether walls at singularities may be removed");
    app.add_flag("--collapse", doCollapse, "Collapse 0 arcs (only works with non-negative arc lengths!");
    app.add_option("--output-constraint-file",
                   constraintFile,
                   "Set this string to generate a quantization constraint file (optional)");
    app.add_option(
        "--output-igm-file", outputIGMFile, "Specify a file to generate an IGM from the quantization into (optional)");
    auto optOutput = app.add_option(
        "--output-hex-file", outputHexFile, "Specify a file to generate a hex mesh from the IGM into (optional)");
    // app.add_option("--scaling", scaling, "Set this double to affect relative quantization scaling");
    app.add_option(
        "--arc-factor", individualArcFactor, "Percent weight of individual arc lengths (vs critical path lengths)");
    app.add_option(
        "--untangling-iter", untanglingIter, "Number of IGM foldover-removal untangling iterations (default 500)");
    app.add_flag("--block-structured",
                 blockStructured,
                 "Whether to collapse first and then requantize to an estimated reasonable amount of hexes");
    app.add_option(
        "--times-minimal", timesMinimalHexes, "By which factor to scale the minimum number of hexes for block structured output");
    app.add_flag("--random-order", randomOrder, "Ordering of collapses");
    app.add_option("--collapse-direction", direction, "Direction of collapses");
    app.add_flag(
        "--optimize-base-mesh",
        optimizeBaseMesh,
        "Optimize the base mesh for IGM generation. More time consuming but better IGM quality and less inversions.");
    app.add_option("--variant",
                   variant,
                   "Pick algorithm variant: \n"
                   "1st bit: 0 flexible / 1 fixed singularities \n 2nd bit: 0 distortion / 1 arcdeviation objective \n "
                   "3rd bit: 1 for singularity padding (only has effect when 1st bit is 0)");
    app.add_option("--num-hexes", nTargetHexes, "Target number of hexes");
#ifdef MC3D_WITH_VIEWER
    bool useViewer = true;
    app.add_flag("--use-viewer", useViewer, "Use viewer to visualize some steps of the pipeline. \n Close the viewer to let algorithm continue.");
#endif
#ifndef QGP3D_WITHOUT_IQP
    int iqpTimeLimit = 180;
    app.add_option("--iqp-time-limit",
                   iqpTimeLimit,
                   "Time limit for quantization IQP solvers in seconds. 0 disables use of IQP solvers.");
#endif

    // Parse cli options
    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    // Create base meshes and add property wrapper
    TetMesh meshRaw;
    MCMesh mcMeshRaw;
    TetMeshProps meshProps(meshRaw, mcMeshRaw);

    meshProps.allocate<ALGO_VARIANT>(variant);
    LOG(INFO) << "Running algorithm variant " << meshProps.get<ALGO_VARIANT>() << " with "
              << (variant % 2 ? "fixed " : ("flexible " + variant / 4 % 2 ? "but padded " : "")) << "singularities"
              << " and " << (variant / 2 % 2 ? "arc-length deviation " : "feature distortion ") << "objective";
    meshProps.allocate<TOUCHED>(true);
    Reader reader(meshProps, inputFile, forceSanitization);
    if (inputHasMCwalls)
        ASSERT_SUCCESS("Reading precomputed MC walls", reader.readSeamlessParamWithWalls());
    else
        ASSERT_SUCCESS("Reading seamless map", reader.readSeamlessParam());

    double vol = 0;
    for (CH tet : meshRaw.cells())
        vol += reader.doubleVolumeUVW(tet);
    double newScaling = std::max(0.001, std::cbrt(nTargetHexes / vol));

    LOG(INFO) << "Scaling parametrization by factor " << newScaling << " to match target hex count";
#ifdef MC3D_WITH_VIEWER
    if (useViewer)
    {
        SingularityInitializer init(meshProps);
        ASSERT_SUCCESS("Determining transitions", init.initTransitions());
        ASSERT_SUCCESS("Determining singularities", init.initSingularities());
        LOG(INFO) << "Visualizing scaled seamless parametrization (close viewer to continue)";
        reader.visualizeParametrization<CHART>(newScaling);
    }
#endif

    MCGenerator mcgen(meshProps);
    if (!inputHasMCwalls)
    {
        // For default usage, the interface is simple to use and requires no property management
        ASSERT_SUCCESS("Tracing and connecting the raw MC",
                       mcgen.traceMC(true, splitSelfadjacent || doCollapse, simulateBC, !constraintFile.empty()));
        MCMeshNavigator(meshProps).assertValidMC(true, true);
        if (!simulateBC)
            ASSERT_SUCCESS("Reducing the raw MC",
                           mcgen.reduceMC(!reduceSingularWalls || doCollapse, splitSelfadjacent || doCollapse, true));
    }
    else
    {
        LOG(INFO) << "Connecting a precomputed MC";

        // For advanced usage, some property management is required.
        // Required/Generated properties are documented for each callable function
        meshProps.allocate<CHILD_CELLS>({});
        meshProps.allocate<CHILD_EDGES>({});
        meshProps.allocate<CHILD_FACES>({});

        // For advanced usage, the library provides specialized classes for each algorithmic step
        SingularityInitializer init(meshProps);
        ASSERT_SUCCESS("Determining transitions", init.initTransitions());
        ASSERT_SUCCESS("Determining singularities", init.initSingularities());
        ASSERT_SUCCESS("Making features consistent", init.makeFeaturesConsistent());

        if (!constraintFile.empty())
        {
            meshProps.allocate<CHART_ORIG>();
            for (CH tet : meshRaw.cells())
                meshProps.set<CHART_ORIG>(tet, meshProps.ref<CHART>(tet));
            meshProps.allocate<TRANSITION_ORIG>();
            for (FH f : meshRaw.faces())
                meshProps.set<TRANSITION_ORIG>(f, meshProps.ref<TRANSITION>(f));
        }

        MCBuilder builder(meshProps);
        ASSERT_SUCCESS("Discovering block structure", builder.discoverBlocks());

        MotorcycleQueue mQ;
        MotorcycleSpawner spawner(meshProps, mQ);
        MotorcycleTracer tracer(meshProps, mQ, simulateBC);

        for (int n = builder.nToroidalBlocks(); n > 0; n = builder.nToroidalBlocks())
        {
            LOG(INFO) << "Splitting toroidal blocks. " << n << " remaining";

            ASSERT_SUCCESS("Spawning torus splitting motorcycles", spawner.spawnTorusSplitMotorcycle());
            ASSERT_SUCCESS("Tracing torus splitting motorcycles", tracer.traceAllMotorcycles());

            auto newWalls = tracer.getNewWalls();
            HFH anyHf = meshRaw.halfface_handle(newWalls.front(), 0);
            ASSERT_SUCCESS("Updating split block", builder.updateSingleBlock(meshRaw.incident_cell(anyHf)));
            tracer.clearNewWalls();
        }

        for (int n = builder.nSelfadjacentBlocks(); splitSelfadjacent && n > 0; n = builder.nSelfadjacentBlocks())
        {
            LOG(INFO) << "Splitting selfadjacent blocks. " << n << " remaining";

            ASSERT_SUCCESS("Spawning selfadjacency splitting motorcycles", spawner.spawnSelfadjacencySplitMotorcycle());
            ASSERT_SUCCESS("Tracing selfadjacency splitting motorcycles", tracer.traceAllMotorcycles());

            auto newWalls = tracer.getNewWalls();
            HFH anyHf = meshRaw.halfface_handle(newWalls.front(), 0);
            HFH anyHfOpp = meshRaw.halfface_handle(newWalls.front(), 1);
            ASSERT_SUCCESS("Updating split block", builder.updateSingleBlock(meshRaw.incident_cell(anyHf)));
            ASSERT_SUCCESS("Updating split block", builder.updateSingleBlock(meshRaw.incident_cell(anyHfOpp)));
            tracer.clearNewWalls();
        }

        ASSERT_SUCCESS("Connecting the MC", builder.connectMCMesh(true, splitSelfadjacent || doCollapse));
        LOG(INFO) << "The MC has " << mcMeshRaw.n_logical_vertices() << " nodes, " << mcMeshRaw.n_logical_edges()
                  << " arcs, " << mcMeshRaw.n_logical_faces() << " patches and " << mcMeshRaw.n_logical_cells()
                  << " blocks";

        MCMeshNavigator(meshProps).assertValidMC(true, true);
        if (!simulateBC)
            ASSERT_SUCCESS("Reducing the raw MC",
                           mcgen.reduceMC(!reduceSingularWalls || doCollapse, splitSelfadjacent || doCollapse));
    }

    MCMeshNavigator(meshProps).assertValidMC(true, true);
    TetRemesher remesher(meshProps);

    if (!inputHasMCwalls)
    {
        LOG(INFO) << "Derefining the base mesh after MC computation...";
        meshProps.allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(true, true, true, false);
        if (optimizeBaseMesh)
        {
            LOG(INFO) << "Optimizing the base mesh before MC collapsing...";
            meshProps.allocate<TOUCHED>(true);
            remesher.collapseAllPossibleEdges(false, true, true, true, 10);
        }
    }

    if (!wallsFile.empty())
    {
        ASSERT_SUCCESS("Writing walls", Writer(meshProps, wallsFile, exactOutput).writeSeamlessParamAndWalls());
    }

    double lowerBound = -DBL_MAX;
    if (doCollapse)
        lowerBound = 0;
    else if (constraintFile.empty())
        lowerBound = 1;

    int minimalHexes = 1;
    if (!constraintFile.empty() || doCollapse || !outputIGMFile.empty() || !outputHexFile.empty())
    {
        ObjectiveBuilder builder(meshProps);
        {
            if (blockStructured)
            {
                StructurePreserver sep(meshProps);
                QuadraticObjective coarsestObjective = builder.simplifiedDistortionObjective(0.001, 0.01);
                ASSERT_SUCCESS("IQP Quantization",
                               IQPQuantizer(meshProps, sep, coarsestObjective).quantize(lowerBound, 60));
                minimalHexes = sep.numHexesInQuantization();
                meshProps.get<MC_MESH_PROPS>()->release<ARC_INT_LENGTH>();
            }
            {
                QuadraticObjective obj
                    = ((variant / 2) % 2)
                          ? ObjectiveBuilder(meshProps).arcLengthDeviationObjective(newScaling)
                          : ObjectiveBuilder(meshProps).simplifiedDistortionObjective(newScaling, individualArcFactor);
                StructurePreserver sep(meshProps);
                ASSERT_SUCCESS("Quantization (ISP initialize)", ISPQuantizer(meshProps, sep, obj).quantize(0.0));
#ifndef QGP3D_WITHOUT_IQP
                if (iqpTimeLimit > 0)
                    ASSERT_SUCCESS("Quantization (IQP improve)",
                                   IQPQuantizer(meshProps, sep, obj).quantize(0.0, iqpTimeLimit));
#endif
            }
        }
    }

    if (doCollapse)
    {
        if (optimizeBaseMesh)
        {
            LOG(INFO) << "Optimizing the base mesh before MC collapsing...";
            meshProps.allocate<TOUCHED>(true);
            remesher.collapseAllPossibleEdges(false, true, true, true, 10);
        }
#ifdef MC3D_WITH_VIEWER
        if (useViewer)
            ASSERT_SUCCESS("Collapsing 0-arcs interactively (close viewer to continue automatically)",
                           MCCollapser(meshProps).interactiveViewCollapse(optimizeBaseMesh, randomOrder, direction));
        else
#endif
            ASSERT_SUCCESS("Collapsing 0-arcs",
                           MCCollapser(meshProps).collapseAllZeroElements(optimizeBaseMesh, randomOrder, direction));
        if (blockStructured && timesMinimalHexes >= 1)
        {
            Q paramVol = 0;
            for (CH tet : meshRaw.cells())
                paramVol += mcgen.rationalVolumeUVW(tet);
            double optimalScaling = std::pow(timesMinimalHexes * minimalHexes / paramVol.get_d(), 1.0 / 3);

            for (EH a : mcMeshRaw.edges())
            {
                double length = 0.0;
                for (HEH he : meshProps.get<MC_MESH_PROPS>()->ref<ARC_MESH_HALFEDGES>(a))
                    length += reader.edgeLengthUVW<CHART>(meshRaw.edge_handle(he));
                meshProps.get<MC_MESH_PROPS>()->set<ARC_DBL_LENGTH>(a, length);
            }
            meshProps.get<MC_MESH_PROPS>()->release<ARC_INT_LENGTH>();
            meshProps.set<ALGO_VARIANT>(3);
            StructurePreserver sep(meshProps);
            ObjectiveBuilder builder(meshProps);
            QuadraticObjective obj = builder.arcLengthDeviationObjective(optimalScaling);
            ASSERT_SUCCESS("Quantization block structured", IQPQuantizer(meshProps, sep, obj).quantize(1.0, 120));
        }
    }
    else if (!constraintFile.empty())
        ASSERT_SUCCESS("Writing constraints", ConstraintWriter(meshProps, constraintFile).writeTetPathConstraints());

    // if (!outputIGMFile.empty() || !outputHexFile.empty())
    {
        IGMGenerator igmgen(meshProps);
        LOG(INFO) << "Generating IGM...";
        auto ret = igmgen.generateBlockwiseIGM(optimizeBaseMesh, untanglingIter, 8);
        if (ret == IGMGenerator::SUCCESS)
            LOG(INFO) << "Generating IGM was successful";
        else if (ret == IGMGenerator::NO_CONVERGENCE)
            LOG(INFO) << "Generated IGM has some inversions";
        else
            LOG(ERROR) << "Generating IGM failed with error code " << ret;
        std::string success = (ret == IGMGenerator::SUCCESS) ? "_clean" : "_dirty";

#ifdef MC3D_WITH_VIEWER
        if (useViewer)
        {
            LOG(INFO) << "Visualizing output IGM (without further distortion optimization!)";
            LOG(INFO) << "Close viewer to continue";
            igmgen.visualizeParametrization<CHART_IGM>();
        }
#endif

        if (!outputIGMFile.empty())
        {
            ASSERT_SUCCESS("Writing IGM in exact rational numbers",
                           Writer(meshProps,
                                  outputIGMFile + "_" + std::to_string(timesMinimalHexes) + success + "_exact.igm",
                                  true)
                               .writeIGMAndWalls());
            ASSERT_SUCCESS(
                "Writing IGM in double precision",
                Writer(meshProps, outputIGMFile + "_" + std::to_string(timesMinimalHexes) + success + ".igm", false)
                    .writeIGM());
        }

        {
            size_t nHexes = StructurePreserver(meshProps).numHexesInQuantization();

            int nSmooth = 1;
            if (nHexes > 10000)
                nSmooth = 10;
            else if (nHexes > 1000)
                nSmooth = 5;
            else if (nHexes > 50)
                nSmooth = 2;

            HexRemesher hexer(meshProps);
            HexMesh hexMeshRaw;
            HexMeshProps hexMeshProps(hexMeshRaw);
            ASSERT_SUCCESS("Extracting hex mesh", hexer.extractHexMesh(hexMeshProps));
            ASSERT_SUCCESS("Smoothing hex mesh (preserves features, projects back to boundary)",
                           hexer.smooth(hexMeshProps, nSmooth));
            if (!outputHexFile.empty())
                ASSERT_SUCCESS("Writing hex mesh",
                               hexer.writeHexMesh(
                                   hexMeshProps, outputHexFile + "_" + std::to_string(timesMinimalHexes) + "_hex.ovm"));
#ifdef MC3D_WITH_VIEWER
            hexer.viewHexMesh(hexMeshProps);
#endif
        }
    }

    return 0;
}
