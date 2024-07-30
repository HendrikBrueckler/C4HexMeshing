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
#include <QGP3D/ISP/ISPQuantizer.hpp>
#include <QGP3D/SeparationChecker.hpp>

#include <C4Hex/Algorithm/MCCollapser.hpp>
#include <C4Hex/Interface/HexRemesher.hpp>
#include <C4Hex/Interface/IGMGenerator.hpp>

#include <C4Hex/Algorithm/IGMInitializer.hpp>
#include <C4Hex/Algorithm/IGMUntangler.hpp>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include <CLI/CLI.hpp>

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

    double scaling = 0.0;
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
    app.add_option("--scaling", scaling, "Set this double to affect relative quantization scaling");
    app.add_option(
        "--untangling-iter", untanglingIter, "Number of IGM foldover-removal untangling iterations (default 500)");
    app.add_flag("--block-structured",
                 blockStructured,
                 "Whether to collapse first and then requantize to an estimated reasonable amount of hexes");
    app.add_option(
        "--times-minimal", timesMinimalHexes, "By which factor to scale number of hexes for block structured");
    app.add_flag("--random-order", randomOrder, "Ordering of collapses");
    app.add_option("--collapse-direction", direction, "Direction of collapses");
    app.add_flag(
        "--optimize-base-mesh",
        optimizeBaseMesh,
        "Optimize the base mesh for IGM generation. More time consuming but better IGM quality and less inversions.");

    // Parse cli options
    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError& e)
    {
        return app.exit(e);
    }

    scaling = std::max(scaling, 0.001);
    // Create base meshes and add property wrapper
    TetMesh meshRaw;
    MCMesh mcMeshRaw;
    TetMeshProps meshProps(meshRaw, mcMeshRaw);

    meshProps.allocate<TOUCHED>(true);
    Reader reader(meshProps, inputFile, forceSanitization);
    if (inputHasMCwalls)
        ASSERT_SUCCESS("Reading precomputed MC walls", reader.readSeamlessParamWithWalls());
    else
        ASSERT_SUCCESS("Reading seamless map", reader.readSeamlessParam());

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
        remesher.collapseAllPossibleEdges(true, true, false, false);
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

    double newScaling = scaling;
    int minimalHexes = 1;
    if (!constraintFile.empty() || doCollapse || !outputIGMFile.empty() || !outputHexFile.empty())
    {
        {
            SeparationChecker sep(meshProps);
            ASSERT_SUCCESS("Quantization", ISPQuantizer(meshProps, sep).quantize(0.0001, lowerBound));
            minimalHexes = sep.numHexesInQuantization();
            vector<double> aLengths;
            for (auto a : mcMeshRaw.edges())
                aLengths.emplace_back(meshProps.get<MC_MESH_PROPS>()->get<ARC_DBL_LENGTH>(a));
            std::sort(aLengths.begin(), aLengths.end());
            double percentileArcLength = aLengths.at(
                std::max(0lu, std::min(aLengths.size() - 1, (size_t)(aLengths.size() * (1.0 - scaling)))));
            newScaling = 0.4 / percentileArcLength;
            LOG(INFO) << "To keep " << scaling * 100 << "% of arcs above length 0.5, a scaling factor of " << newScaling
                      << " for collapsing was chosen";
        }
        SeparationChecker sep(meshProps);
        ASSERT_SUCCESS("Quantization", ISPQuantizer(meshProps, sep).quantize(newScaling, lowerBound));
        if (doCollapse && !MCCollapser(meshProps).hasZeroLengthArcs())
        {
            LOG(INFO) << "No 0-arcs, nothing to collapse, exiting...";
            return 0;
        }
        MCCollapser(meshProps).markZeros();
    }

    if (doCollapse)
    {
        if (optimizeBaseMesh)
        {
            LOG(INFO) << "Optimizing the base mesh before MC collapsing...";
            meshProps.allocate<TOUCHED>(true);
            remesher.collapseAllPossibleEdges(false, true, true, true, 10);
        }
        ASSERT_SUCCESS("Collapsing 0-arcs",
                       MCCollapser(meshProps).collapseAllZeroElements(optimizeBaseMesh, randomOrder, direction));
        if (blockStructured)
        {
            Q paramVol = 0;
            for (CH tet : meshRaw.cells())
                paramVol += mcgen.rationalVolumeUVW(tet);
            double optimalScaling = std::pow(timesMinimalHexes * minimalHexes / paramVol.get_d(), 1.0 / 3);

            SeparationChecker sep(meshProps);
            ASSERT_SUCCESS("Quantization block structured", ISPQuantizer(meshProps, sep).quantize(optimalScaling, 1.0));
        }
    }
    else if (!constraintFile.empty())
        ASSERT_SUCCESS("Writing constraints", ConstraintWriter(meshProps, constraintFile).writeTetPathConstraints());

    if (!outputIGMFile.empty() || !outputHexFile.empty())
    {
        IGMGenerator igmgen(meshProps);
        LOG(INFO) << "Generating IGM...";
        auto ret = igmgen.generateBlockwiseIGM(optimizeBaseMesh, untanglingIter, 40);
        if (ret == IGMGenerator::SUCCESS)
            LOG(INFO) << "Generating IGM was successful";
        else if (ret == IGMGenerator::NO_CONVERGENCE)
            LOG(INFO) << "Generated IGM has some inversions";
        else
            LOG(ERROR) << "Generating IGM failed with error code " << ret << ", aborting...";

        if (!outputIGMFile.empty())
        {
            // remesher.remeshToImproveAngles(true, true, TetMeshManipulator::QualityMeasure::VL_RATIO);
            ASSERT_SUCCESS("Writing IGM", Writer(meshProps, outputIGMFile + "_exact.igm", true).writeIGMAndWalls());
            ASSERT_SUCCESS("Writing IGM", Writer(meshProps, outputIGMFile + ".igm", false).writeIGM());
        }

        if (!outputHexFile.empty())
        {
            int nHexes = SeparationChecker(meshProps).numHexesInQuantization();
            int nsub = 0;
            if (nHexes > 100000)
                nsub = 2;
            else if (nHexes > 10000)
                nsub = 3;
            else if (nHexes > 5000)
                nsub = 4;
            else if (nHexes > 1000)
                nsub = 5;
            else if (nHexes > 500)
                nsub = 6;
            else
                nsub = 10;

            int nSmooth = 0;
            if (nHexes > 10000)
                nSmooth = 3;
            else if (nHexes > 1000)
                nSmooth = 2;
            else if (nHexes > 50)
                nSmooth = 1;

            HexRemesher hexer(meshProps);
            {
                PolyMesh polyMCMesh;
                PolyMeshProps polyMCMeshProps(polyMCMesh);
                ASSERT_SUCCESS("Extracting MC mesh", hexer.extractMCMesh(polyMCMeshProps));
                ASSERT_SUCCESS("Writing MC mesh", hexer.writePolyHexMesh(polyMCMeshProps, outputHexFile + "_MC.ovm"));
            }

            {
                HexMesh hexMeshRaw;
                HexMeshProps hexMeshProps(hexMeshRaw);
                ASSERT_SUCCESS("Extracting hex mesh", hexer.extractHexMesh(hexMeshProps));
                ASSERT_SUCCESS("Smoothing hex mesh", hexer.smoothSurface(hexMeshProps, 0));
                ASSERT_SUCCESS("Writing hex mesh", hexer.writeHexMesh(hexMeshProps, outputHexFile + "_hex.ovm"));
            }

            PolyMesh polyHexMesh;
            PolyMeshProps polyMeshProps(polyHexMesh);

            ASSERT_SUCCESS("Extracting poly hex mesh", hexer.extractPolyHexMesh(polyMeshProps, nsub));
            ASSERT_SUCCESS("Smoothing poly hex mesh", hexer.smoothSurface(polyMeshProps, nSmooth));
            ASSERT_SUCCESS("Writing poly hex mesh", hexer.writePolyHexMesh(polyMeshProps, outputHexFile + "_poly.ovm"));
        }
    }

    return 0;
}
