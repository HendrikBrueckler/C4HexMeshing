#include "C4Hex/Algorithm/MCCollapser.hpp"

#include "C4Hex/Algorithm/EmbeddingCollapser.hpp"
#include "C4Hex/Algorithm/MCSmoother.hpp"
#include "C4Hex/Algorithm/PathRouter.hpp"
#include "C4Hex/Algorithm/SurfaceRouter.hpp"

#include <MC3D/Algorithm/MCReducer.hpp>
#include <MC3D/Algorithm/TetRemesher.hpp>
#include <MC3D/Interface/Writer.hpp>

#include <chrono>
#include <iostream>

#ifdef MC3D_WITH_VIEWER
#include <QGP3D/IQP/IQPQuantizer.hpp>
#include <QGP3D/ISP/ISPQuantizer.hpp>
#include <QGP3D/ObjectiveBuilder.hpp>
#include <QGP3D/StructurePreserver.hpp>
#include <misc/cpp/imgui_stdlib.h>
#include <settings/AppState.h>
#include <util/ImGuiUtil.h>
#include <volumeshOS.h>
#endif

namespace c4hex
{
#ifdef MC3D_WITH_VIEWER
using namespace qgp3d;
#endif

MCCollapser::MCCollapser(TetMeshProps& meshProps_)
    : TetMeshNavigator(meshProps_), TetMeshManipulator(meshProps_), MCMeshNavigator(meshProps_),
      MCMeshManipulator(meshProps_), _refiner(meshProps_)
{
}

bool MCCollapser::hasZeroLengthArcs() const
{
    for (EH a : mcMeshProps().mesh().edges())
        if (isZeroArc(a))
            return true;
    return false;
}

bool MCCollapser::collapseIsLocked(const HEH& ha) const
{
    auto& mcMesh = mcMeshProps().mesh();
    VH nFrom = mcMesh.from_vertex_handle(ha);

    if (mcMesh.is_boundary(nFrom) && !mcMesh.is_boundary(ha))
        return true;

    if (mcMeshProps().get<IS_CRITICAL_N>(nFrom))
        return true;

    bool nodeInCriticalLinks = containsMatching(
        mcMesh.vertex_edges(nFrom), [&, this](const EH& a2) { return mcMeshProps().get<IS_CRITICAL_A>(a2); });
    EH a = mcMesh.edge_handle(ha);
    if (nodeInCriticalLinks && !mcMeshProps().get<IS_CRITICAL_A>(a))
        return true;

    bool nodeInCriticalRegion = containsMatching(
        mcMesh.vertex_faces(nFrom), [&, this](const FH& p) { return mcMeshProps().get<IS_CRITICAL_P>(p); });
    bool arcInCriticalRegion = containsMatching(mcMesh.halfedge_faces(ha),
                                                [&, this](const FH& p) { return mcMeshProps().get<IS_CRITICAL_P>(p); });
    if (nodeInCriticalRegion && !arcInCriticalRegion)
        return true;

    if (_direction == 0 && !_randomOrder && mcMeshProps().isAllocated<BLOCK_COLLAPSE_DIR>())
    {
        // Necessary to iterate like this, because selfadjacency messes with edgecelliter
        set<CH> bs;
        for (FH p : mcMesh.edge_faces(mcMesh.edge_handle(ha)))
            for (CH b : mcMesh.face_cells(p))
                if (b.is_valid())
                    bs.insert(b);

        if (!containsMatching(bs,
                              [&, this](const CH& b)
                              {
                                  UVWDir dir = halfarcDirInBlock(ha, b);
                                  return ((dir | -dir) & mcMeshProps().get<BLOCK_COLLAPSE_DIR>(b)) != UVWDir::NONE;
                              }))
            return true;
    }

    // If there is a second arc that has the same vertices => locked
    auto ns = mcMesh.halfedge_vertices(ha);
    return containsMatching(mcMesh.outgoing_halfedges(ns[0]),
                            [&](const HEH& ha2)
                            {
                                return mcMesh.edge_handle(ha2) != mcMesh.edge_handle(ha)
                                       && mcMesh.to_vertex_handle(ha2) == ns[1] && isZeroArc(mcMesh.edge_handle(ha2));
                            });
}

bool MCCollapser::collapseMovesSingularity(const HEH& ha) const
{
    auto& mcMesh = mcMeshProps().mesh();
    VH nFrom = mcMesh.from_vertex_handle(ha);

    auto type = mcMeshProps().nodeType(nFrom);

    if (type.first == SingularNodeType::SINGULAR)
        return true;

    if (type.first == SingularNodeType::SEMI_SINGULAR && !mcMeshProps().get<IS_SINGULAR>(mcMesh.edge_handle(ha)))
        return true;

    return false;
}

HEH MCCollapser::preferredCollapseHalfarc(const EH& a) const
{
    auto& mcMesh = mcMeshProps().mesh();
    HEH ha0 = mcMesh.halfedge_handle(a, 0);
    HEH ha1 = mcMesh.halfedge_handle(a, 1);
    bool locked0 = collapseIsLocked(ha0);
    bool locked1 = collapseIsLocked(ha1);
    if (locked0 && locked1)
        return HEH();
    if (locked0)
        return ha1;
    if (locked1)
        return ha0;

    bool movesSing0 = collapseMovesSingularity(ha0);
    bool movesSing1 = collapseMovesSingularity(ha0);
    if (movesSing0 && !movesSing1)
        return ha1;
    if (movesSing1 && !movesSing0)
        return ha0;

    // If none is locked estimate which direction causes less work/MC distortion
    if (_direction == 0 && mcMeshProps().isAllocated<BLOCK_COLLAPSE_DIR>())
    {
        CH b = *mcMesh.ec_iter(a);
        UVWDir dir0 = halfarcDirInBlock(ha0, b);
        UVWDir dir1 = -dir0;
        UVWDir dirB = mcMeshProps().get<BLOCK_COLLAPSE_DIR>(b);
        if (dir0 == dirB)
            return ha0;
        else if (dir1 == dirB)
            return ha1;
    }

    if (_direction == 2)
    {
        // Variant 1: lower arc valence on source node
        if (mcMesh.valence(mcMesh.from_vertex_handle(ha0)) >= mcMesh.valence(mcMesh.to_vertex_handle(ha0)))
            return ha1;
    }
    else if (_direction == 1)
    {
        return ((double)rand() / INT_MAX > 0.5 ? ha0 : ha1);
    }

    return ha0;
}

MCCollapser::RetCode MCCollapser::interactiveViewCollapse(bool optimize, bool randomOrder, int direction)
{
#ifndef MC3D_WITH_VIEWER
    return collapseAllZeroElements(optimize, randomOrder, direction);
#else
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    LOG(INFO) << "Collapsing randomOrder: " << randomOrder << " and direction: " << direction;
    _randomOrder = randomOrder;
    _direction = direction;

    MCReducer reducer(meshProps());
    TetRemesher remesher(meshProps());

    {
        for (EH a : mcMesh.edges())
        {
            if (mcMeshProps().get<ARC_INT_LENGTH>(a) != 0)
                continue;
            _numZeroAs++;
        }
        for (FH p : mcMesh.faces())
        {
            auto hasByDir = halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0));
            for (auto& kv : hasByDir)
            {
                bool zeroSide = true;
                for (HEH ha : kv.second)
                {
                    if (mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha)) != 0)
                    {
                        zeroSide = false;
                        break;
                    }
                }
                if (zeroSide)
                {
                    _numZeroPs++;
                    break;
                }
            }
        }
        for (CH b : mcMesh.cells())
        {
            for (UVWDir dir : {UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W})
            {
                bool zeroSide = true;
                auto dir2arcs = mcMeshProps().get<BLOCK_ALL_ARCS>(b);
                for (auto dir2 : {dir, -dir})
                {
                    for (EH a : dir2arcs[dir2])
                        if (mcMeshProps().get<ARC_INT_LENGTH>(a) != 0)
                        {
                            zeroSide = false;
                            break;
                        }
                    if (!zeroSide)
                        break;
                }
                if (zeroSide)
                {
                    _numZeroBs++;
                    break;
                }
            }
        }
        LOG(INFO) << _numZeroAs << " arcs, " << _numZeroPs << " patches, " << _numZeroBs << " blocks to collapse";
    }

    assignCollapseDirs();

    OVM::Vec4f color0A(0.24706, 0.58039, 0.67059, 1.0);
    OVM::Vec4f color1A(0.58039, 0.00000, 1.00000, 1.0);
    OVM::Vec4f color2A(0.99216, 0.14118, 0.72941, 1.0);
    OVM::Vec4f color3A(0.99608, 0.63529, 0.00000, 1.0);
    OVM::Vec4f color0P(0.870, 1, 0.996, 0.31);
    OVM::Vec4f color4P(1, 0.011, 0.749, 0.74);
    OVM::Vec4f color5P(1, 0.917, 0.376, 0.8);
    OVM::Vec4f color6P(1, 0.917, 0.376, 0.4);
    OVM::Vec4f colorHighlight(1.0, 0.3, 0.3, 0.5);
    auto setZeroColors = [&, this](bool highlight = false)
    {
        if (meshProps().isAllocated<COLOR_F>())
            meshProps().release<COLOR_F>();
        if (meshProps().isAllocated<COLOR_E>())
            meshProps().release<COLOR_E>();
        if (meshProps().isAllocated<COLOR_V>())
            meshProps().release<COLOR_V>();
        meshProps().allocate<COLOR_F>(Vec4f(1.0f, 1.0f, 1.0f, 0.0f));
        meshProps().allocate<COLOR_E>(Vec4f(1.0f, 1.0f, 1.0f, 0.0f));
        meshProps().allocate<COLOR_V>(Vec4f(1.0f, 1.0f, 1.0f, 0.0f));

        for (FH f : tetMesh.faces())
        {
            bool hasTouched = false;
            for (VH v : tetMesh.face_vertices(f))
                if (meshProps().get<TOUCHED>(v))
                {
                    hasTouched = true;
                    break;
                }
            FH p = meshProps().get<MC_PATCH>(f);
            if (!p.is_valid())
                meshProps().set<COLOR_F>(f, Vec4f(1.0f, 1.0f, 1.0f, 0.0f));
            else
            {
                if (highlight && hasTouched)
                {
                    meshProps().set<COLOR_F>(f, colorHighlight);
                }
                else
                    switch (highlight ? mcMeshProps().get<MARK_P>(p) : (mcMeshProps().get<MARK_P>(p) % 1000))
                    {
                    case 0:
                        meshProps().set<COLOR_F>(f, color0P);
                        break;
                    case 4:
                        meshProps().set<COLOR_F>(f, color4P);
                        break;
                    case 5:
                        meshProps().set<COLOR_F>(f, color5P);
                        break;
                    case 6:
                        meshProps().set<COLOR_F>(f, color6P);
                        break;
                    default:
                        meshProps().set<COLOR_F>(f, Vec4f(1.0f, 1.0f, 1.0f, 0.0f));
                    }
            }
        }
        for (EH e : tetMesh.edges())
        {
            bool hasTouched = false;
            for (VH v : tetMesh.edge_vertices(e))
                if (meshProps().get<TOUCHED>(v))
                {
                    hasTouched = true;
                    break;
                }
            EH a = meshProps().get<MC_ARC>(e);
            if (!a.is_valid())
                meshProps().reset<COLOR_E>(e);
            else if (highlight && hasTouched)
            {
                meshProps().set<COLOR_E>(e, colorHighlight);
            }
            else
                switch (highlight ? mcMeshProps().get<MARK_A>(a) : (mcMeshProps().get<MARK_A>(a) % 1000))
                {
                case 0:
                    meshProps().set<COLOR_E>(e, color0A);
                    break;
                case 1:
                    meshProps().set<COLOR_E>(e, color1A);
                    break;
                case 2:
                    meshProps().set<COLOR_E>(e, color2A);
                    break;
                case 3:
                    meshProps().set<COLOR_E>(e, color3A);
                    break;
                default:
                    meshProps().set<COLOR_E>(e, colorHighlight);
                }
        }
        for (VH v : tetMesh.vertices())
        {
            VH n = meshProps().get<MC_NODE>(v);
            if (!n.is_valid())
                meshProps().reset<COLOR_V>(v);
            else
            {
                if (highlight && meshProps().get<TOUCHED>(v))
                {
                    meshProps().set<COLOR_V>(v, colorHighlight);
                }
                else
                    switch (highlight ? mcMeshProps().get<MARK_N>(n) : (mcMeshProps().get<MARK_N>(n) % 1000))
                    {
                    case 0:
                        meshProps().set<COLOR_V>(v, color0A);
                        break;
                    case 3:
                        meshProps().set<COLOR_V>(v, color3A);
                        break;
                    default:
                        meshProps().set<COLOR_V>(v, colorHighlight);
                    }
            }
        }
    };

    bool newCollapseDir = true;
    bool change = true;
    auto executeUntilNextChange = [&, this]()
    {
        do
        {
            if (collapseNextPillowBlock() || collapseNextPillowPatch() || bisectNextBlockByPillowPatch()
                || collapseNextZeroArc() || bisectNextAlmostPillowPatch(false) || bisectNextAlmostPillowBlock())
            {
                change = true;
                return;
            }

            if (!newCollapseDir && hasZeroLengthArcs())
            {
                assignCollapseDirs();
                newCollapseDir = true;
            }
            else if (bisectNextZeroPatch(true))
            {
                newCollapseDir = false;
                change = true;
                return;
            }
            else
                newCollapseDir = false;
        } while (newCollapseDir);
    };

    double cylinderRadius = 0.0;
    double sphereRadius = 0.0;
    {
        double avgEdgeLength = 0.0;
        double minEdgeLength = DBL_MAX;
        vector<double> edgeLengths;
        edgeLengths.reserve(tetMesh.n_logical_edges());
        for (EH e : tetMesh.edges())
        {
            double length = tetMesh.length(e);
            minEdgeLength = std::min(minEdgeLength, length);
            avgEdgeLength += length;
            edgeLengths.emplace_back(length);
        }
        avgEdgeLength /= tetMesh.n_logical_edges();
        std::sort(edgeLengths.begin(), edgeLengths.end());
        double percentileEdgeLength = edgeLengths.at((size_t)(edgeLengths.size() * 0.01));

        Vec3d bboxMin(DBL_MAX, DBL_MAX, DBL_MAX);
        Vec3d bboxMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        for (VH v : tetMesh.vertices())
        {
            Vec3d pos = tetMesh.vertex(v);
            for (int i = 0; i < 3; i++)
            {
                if (pos[i] > bboxMax[i])
                    bboxMax[i] = pos[i];
                if (pos[i] < bboxMin[i])
                    bboxMin[i] = pos[i];
            }
        }
        Vec3d diagonal = (bboxMax - bboxMin);
        double bboxDiameter = (diagonal).length();
        // double scalefactor = 7.0 / std::max(diagonal[0], std::max(diagonal[1], diagonal[2]));
        percentileEdgeLength = std::max(percentileEdgeLength, bboxDiameter / 150.0);
        cylinderRadius = percentileEdgeLength * 0.4;
        sphereRadius = percentileEdgeLength * 0.8;
    }

    markZeros();
    setZeroColors();
    volumeshOS::VMesh mesh;
    auto updateVMesh = [&, this](bool onlyColors = false)
    {
        if (!onlyColors)
        {
            if (!mesh.is_valid())
            {
                LOG(INFO) << "loading vismesh";
                mesh = volumeshOS::load(&tetMesh);
                volumeshOS::Internal::AppState::settings.camera.position = {3.8, 7.5, 20.0};
                volumeshOS::Internal::AppState::settings.camera.mode = volumeshOS::CameraMode::FLY;
                volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                volumeshOS::set_camera_position(3.8, 7.5, 20.0);

                mesh.set_cell_rounding(0.0);
                mesh.set_scale(1.00);
                mesh.use_base_color(false);
                mesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
                mesh.set_shading_mode(volumeshOS::ShadingMode::FLAT);
                mesh.set_ambient(0.7f);
                mesh.set_diffuse(0.3f);
                mesh.set_specular(0.1f);
                mesh.set_specular_coefficient(8.0f);
                mesh.use_two_sided_lighting(true);
            }
            else
            {
                volumeshOS::update(mesh, &tetMesh);
                // mesh.set_cell_rounding(0.0);
                // mesh.set_scale(1.00);
                // mesh.use_base_color(false);
                // mesh.set_lighting_mode(volumeshOS::LightingMode::PHONG);
                // mesh.set_shading_mode(volumeshOS::ShadingMode::FLAT);
                // mesh.set_ambient(0.7f);
                // mesh.set_diffuse(0.3f);
                // mesh.set_specular(0.1f);
                // mesh.set_specular_coefficient(8.0f);
                // mesh.use_two_sided_lighting(true);
            }
        }
        volumeshOS::remove_shapes();
        // Actual mesh: ignore block 3
        for (HFH hf : tetMesh.halffaces())
        {
            FH f = tetMesh.face_handle(hf);
            mesh.set_color(hf, meshProps().get<COLOR_F>(f));
        }
        for (EH e : tetMesh.edges())
        {
            if (meshProps().get<COLOR_E>(e)[3] > 0.0)
            {
                auto highlight = mesh.add_shape<volumeshOS::VCylinder>();
                assert(!tetMesh.is_deleted(e));
                auto vs = tetMesh.edge_vertices(e);
                Vec3d dir = tetMesh.vertex(vs[1]) - tetMesh.vertex(vs[0]);
                highlight.set_position(tetMesh.vertex(vs[0]) + 0.5 * dir);
                highlight.set_color(meshProps().get<COLOR_E>(e));
                highlight.set_direction(dir.normalized());
                highlight.set_scale(OVM::Vec3d({cylinderRadius, dir.norm(), cylinderRadius}));
            }
        }
        for (VH v : tetMesh.vertices())
        {
            if (meshProps().get<COLOR_V>(v)[3] > 0.0)
            {
                Vec3d pos = tetMesh.vertex(v);
                auto highlight = mesh.add_shape<volumeshOS::VSphere>();
                highlight.set_position(pos);
                highlight.set_color(meshProps().get<COLOR_V>(v));
                highlight.set_scale((float)sphereRadius);
            }
        }
    };

    std::string screenshotFilename = "screenshot";
    struct ExportOptions
    {
        int width = -1;                  // viewport width by default
        int height = -1;                 // viewport height by default
        bool include_background = true;  // include background in image
        bool include_shapes = true;      // include shapes in image
        bool include_ground = true;      // include ground in image
        bool ground_shadow_only = false; // if ground and shadows are active, export only the area in shadow
    };
    static int scrnum = 0;
    auto makeScreenShot = [&]()
    {
        volumeshOS::ExportOptions exopt;
        exopt.width = volumeshOS::get_viewport_width();   // viewport width by default
        exopt.height = volumeshOS::get_viewport_height(); // viewport height by default
        exopt.include_background = true;                  // include background in image
        exopt.include_shapes = true;                      // include shapes in image
        exopt.include_ground = true;                      // include ground in image
        exopt.ground_shadow_only = true; // if ground and shadows are active, export only the area in shadow

        auto str = std::to_string(scrnum++);
        auto num = std::string(5 - std::min(5ul, str.length()), '0') + str;
        volumeshOS::export_image("./" + screenshotFilename + num + ".png", exopt);
    };

    volumeshOS::use_transparency(true);
    volumeshOS::use_shadows(true);
    volumeshOS::set_shape_lighting_mode(volumeshOS::LightingMode::PHONG);
    volumeshOS::set_shape_ambient(0.75f);
    volumeshOS::set_shape_diffuse(0.25f);
    volumeshOS::set_shape_specular(0.3f);
    volumeshOS::set_shape_specular_coefficient(8.0f);
    volumeshOS::set_sky_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::set_ground_color(OVM::Vec3f{1.0f, 1.0f, 1.0f});
    volumeshOS::Internal::AppState::settings.light.color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.sky_color = {1.0f, 1.0f, 1.0f};
    volumeshOS::Internal::AppState::settings.sky.fog_density = 0.0f;
    volumeshOS::Internal::AppState::settings.post_processing.active = false;
    volumeshOS::Internal::AppState::settings.ground.use_pbr = false;
    volumeshOS::Internal::AppState::settings.shadow.shadow_strength = 0.24f;
    volumeshOS::Internal::AppState::settings.shadow.penumbra_scale = 22.0f;
    updateVMesh();

    bool animating = false;
    int stage = 0;
    set<FH> psOptimized;
    set<EH> asOnQ;
    list<EH> as;
    map<EH, int> a2timesOptimized;

    int quantEnabled = 1;
    bool justFinished = false;
    auto animate = [&]()
    {
        if (stage == 0)
        {
            markZeros();

            float tempAlpha = 1.0;
            std::swap(tempAlpha, color0P[3]);
            setZeroColors();
            updateVMesh(true);
            makeScreenShot();
            std::swap(tempAlpha, color0P[3]);
            stage++;
        }
        else if (stage == 1)
        {
            setZeroColors();
            updateVMesh(true);
            makeScreenShot();
            change = true;
            stage++;
        }
        else if (stage == 2)
        {
            if (change)
                newCollapseDir = false;
            change = false;

            remesher.collapseAllPossibleEdges(true, true, false, false);
            executeUntilNextChange();
            setZeroColors();
            markZeros();
            updateVMesh();
            makeScreenShot();
            if (!change)
            {
                stage++;
                reducer.init(true, true, true);
            }
        }
        else if (stage == 3)
        {
            if (reducer.isReducible())
            {
                reducer.removeNextPatch();
                setZeroColors();
                markZeros();
                updateVMesh();
                makeScreenShot();
            }
            else
            {
                psOptimized.clear();
                as.clear();
                asOnQ.clear();
                a2timesOptimized.clear();
                for (EH a : mcMesh.edges())
                    as.push_back(a);
                asOnQ = set<EH>(as.begin(), as.end());
                stage++;
            }
        }
        else if (stage == 4)
        {
            MCSmoother(meshProps()).smoothMC();
            setZeroColors();
            markZeros();
            updateVMesh();
            makeScreenShot();
            animating = false;
            stage = 0;
            quantEnabled = 2 * !zeroElementsRemain();
            justFinished = true;
        }
    };

    meshProps().allocate<TOUCHED>(true);
    for (VH v : tetMesh.vertices())
        meshProps().set<TOUCHED>(v, false);
    double vol = 0.0;
    for (CH tet : tetMesh.cells())
        vol += doubleVolumeUVW(tet);
    int nHexes = std::max(1,
                          mcMeshProps().isAllocated<ARC_INT_LENGTH>()
                              ? (int)StructurePreserver(meshProps()).numHexesInQuantization()
                              : (int)std::round(vol));
    int focusMode = 0;
    int selectionID = -1;
    bool init = false;
    volumeshOS::on_gui_render(
        [&, this]()
        {
            if (!init)
            {
                markZeros();
                setZeroColors();
                updateVMesh(true);
                init = true;
            }
            if (animating)
                animate();
            ImGui::Begin("MyPanel");
            if (ImGui::Button("UPDATE COLORING") || justFinished)
            {
                markZeros();
                setZeroColors();
                updateVMesh(true);
                justFinished = false;
            }
            // if (ImGui::Button("SCREENSHOT"))
            // {
            //     makeScreenShot();
            // }
            if (quantEnabled)
            {
                if (ImGui::Button("QUANTIZE (using below parameters"))
                {
                    if (quantEnabled == 2 && meshProps().isAllocated<ALGO_VARIANT>())
                        meshProps().set<ALGO_VARIANT>(3);
                    StructurePreserver sep(meshProps());
                    ObjectiveBuilder builder(meshProps());
                    double scaling = std::max(0.001, std::cbrt(nHexes / vol));
                    QuadraticObjective obj
                        = (quantEnabled == 2
                           || (meshProps().isAllocated<ALGO_VARIANT>() && ((meshProps().get<ALGO_VARIANT>() / 2) % 2)))
                              ? ObjectiveBuilder(meshProps()).arcLengthDeviationObjective(scaling)
                              : ObjectiveBuilder(meshProps()).simplifiedDistortionObjective(scaling, 0.01);
                    ISPQuantizer(meshProps(), sep, obj).quantize(0.00);
#ifndef QGP3D_WITHOUT_IQP
                    IQPQuantizer(meshProps(), sep, obj).quantize(0.0, 30);
#endif
                    markZeros();
                    setZeroColors();
                    updateVMesh();
                }
            }
            else
            {
                if (ImGui::Button("[Quantization disabled after collapsing]"))
                {
                }
            }
            if (ImGui::Button("SINGLE COLLAPSING STEP"))
            {
                if (change)
                    newCollapseDir = false;
                change = false;

                remesher.collapseAllPossibleEdges(true, true, false, false);
                executeUntilNextChange();
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                quantEnabled = 2 * !zeroElementsRemain();
                if (quantEnabled == 2)
                    justFinished = true;
            }
            if (ImGui::Button("10 STEPS"))
            {
                change = true;
                for (int i = 0; i < 10 && change; i++)
                {
                    if (change)
                        newCollapseDir = false;
                    change = false;

                    remesher.collapseAllPossibleEdges(true, true, false, false);
                    executeUntilNextChange();
                }
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                quantEnabled = 2 * !zeroElementsRemain();
                if (quantEnabled == 2)
                    justFinished = true;
            }
            if (ImGui::Button("50 STEPS"))
            {
                change = true;
                for (int i = 0; i < 50 && change; i++)
                {
                    if (change)
                        newCollapseDir = false;
                    change = false;

                    remesher.collapseAllPossibleEdges(true, true, false, false);
                    executeUntilNextChange();
                }
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                quantEnabled = 2 * !zeroElementsRemain();
                if (quantEnabled == 2)
                    justFinished = true;
            }
            if (ImGui::Button("FINISH COLLAPSING"))
            {
                change = true;
                while (change)
                {
                    if (change)
                        newCollapseDir = false;
                    change = false;

                    remesher.collapseAllPossibleEdges(true, true, false, false);
                    executeUntilNextChange();
                }
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                quantEnabled = 2 * !zeroElementsRemain();
                if (quantEnabled == 2)
                    justFinished = true;
            }
            if (ImGui::Button("REDUCE (Postprocess after collapsing)"))
            {
                change = true;
                reducer.init(true, true, true);
                while (reducer.isReducible())
                {
                    reducer.removeNextPatch();
                }
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                justFinished = true;
            }
            if (ImGui::Button("DECIMATE_TETMESH (Reduce cell count)"))
            {
                assertValidMC(false, false);
                meshProps().allocate<TOUCHED>(true);
                remesher.collapseAllPossibleEdges(true, true, false, false);
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                justFinished = true;
            }
            if (ImGui::Button("REMESH_TETMESH (Improve cell geometry)"))
            {
                remesher.remeshToImproveAngles(true, false, TetRemesher::QualityMeasure::ANGLES, -1);
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                justFinished = true;
            }
            if (ImGui::Button("SMOOTH MC (Postprocess for less jagged MC geometry)"))
            {
                if (!zeroElementsRemain())
                {
                    reducer.init(true, true, true);
                    while (reducer.isReducible())
                        reducer.removeNextPatch();
                }
                MCSmoother(meshProps()).smoothMC();
                setZeroColors();
                markZeros();
                updateVMesh();
                assertValidMC(false, false);
                justFinished = true;
            }
            if (ImGui::Button("ANIMATE COLLAPSING (and make screenshots)"))
            {
                animating = true;
                quantEnabled = 0;
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Visualization", 12))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "Target hexes for quantization",
                    [&] { ImGui::InputInt("##Target hexes for quantization", &nHexes); });

                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "ScreenshotFilename", [&] { ImGui::InputText("##ScreenshotFilename", &screenshotFilename); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "CylinderRadius", [&] { ImGui::InputDouble("##CylinderRadius", &cylinderRadius); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "SphereRadius", [&] { ImGui::InputDouble("##SphereRadius", &sphereRadius); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color0A", color0A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color1A", color1A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color2A", color2A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color3A", color3A.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color0P", color0P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color4P", color4P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color5P", color5P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "", [&] { ImGui::ColorEdit4("Color6P", color6P.data(), ImGuiColorEditFlags_NoInputs); });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }
            if (volumeshOS::Internal::ImGuiUtil::begin_menu_with_background("Focus", 3))
            {
                volumeshOS::Internal::ImGuiUtil::menu_item_filled("Select by ID",
                                                                  [&]
                                                                  {
                                                                      constexpr const char* selection_modes[] = {
                                                                          "Off", "Vertices", "Edges", "Faces", "Cells"};
                                                                      ImGui::Combo("##SelectionMode",
                                                                                   &focusMode,
                                                                                   selection_modes,
                                                                                   IM_ARRAYSIZE(selection_modes),
                                                                                   IM_ARRAYSIZE(selection_modes));
                                                                  });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "ID", [&] { ImGui::InputInt("##ManualSelectionID", &selectionID); });
                volumeshOS::Internal::ImGuiUtil::menu_item_filled(
                    "",
                    [&]
                    {
                        if (ImGui::Button("Focus"))
                        {
                            volumeshOS::set_camera_mode(volumeshOS::CameraMode::ORBIT);
                            Vec3d target;
                            switch (focusMode)
                            {
                            case 0:
                                volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 1:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_vertices())
                                    target = tetMesh.vertex(VH(selectionID));
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 2:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_edges())
                                {
                                    auto vs = tetMesh.edge_vertices(EH(selectionID));
                                    target = 0.5 * tetMesh.vertex(vs[0]) + 0.5 * tetMesh.vertex(vs[1]);
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 3:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_faces())
                                {
                                    vector<Vec3d> vs;
                                    for (VH v : tetMesh.face_vertices(FH(selectionID)))
                                        vs.push_back(tetMesh.vertex(v));
                                    for (Vec3d pos : vs)
                                        target += pos;
                                    target /= vs.size();
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            case 4:
                                if (selectionID >= 0 && selectionID <= (int)tetMesh.n_cells())
                                {
                                    vector<Vec3d> vs;
                                    for (VH v : tetMesh.cell_vertices(CH(selectionID)))
                                        vs.push_back(tetMesh.vertex(v));
                                    for (Vec3d pos : vs)
                                        target += pos;
                                    target /= vs.size();
                                }
                                else
                                    volumeshOS::set_camera_mode(volumeshOS::CameraMode::FLY);
                                break;
                            default:
                                break;
                            }
                            volumeshOS::set_camera_target(mesh.get_transformed_point(target));
                        }
                    });
                volumeshOS::Internal::ImGuiUtil::end_menu();
            }

            ImGui::End();
        });

    volumeshOS::open();

    if (zeroElementsRemain())
    {
        change = true;
        while (change)
        {
            if (change)
                newCollapseDir = false;
            change = false;

            remesher.collapseAllPossibleEdges(true, true, false, false);
            executeUntilNextChange();
        }
    }

    if (zeroElementsRemain())
        return COLLAPSE_ERROR;

    reducer.init(true, true, true);
    while (reducer.isReducible())
    {
        reducer.removeNextPatch();
        assertValidMC(false, true);
    }

    assertValidMC(true, true);
    meshProps().allocate<TOUCHED>(true);

    if (optimize)
    {
        remesher.collapseAllPossibleEdges(false, true, true, true, 10.0);
        assertValidMC(true, true);
        remesher.remeshToImproveAngles(true, false, TetRemesher::QualityMeasure::ANGLES);
        assertValidMC(true, true);
        MCSmoother(meshProps()).smoothMC();
    }

    assertValidMC(true, true);
    return SUCCESS;
#endif
}

MCCollapser::RetCode MCCollapser::collapseAllZeroElements(bool optimize, bool randomOrder, int direction)
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    TetRemesher remesher(meshProps());

    assertValidMC(true, true);

    LOG(INFO) << "Collapsing randomOrder: " << randomOrder << " and direction: " << direction;
    _randomOrder = randomOrder;
    _direction = direction;

    MCReducer reducer(meshProps());

    countZeroElements(_numZeroAs, _numZeroPs, _numZeroBs);
    LOG(INFO) << _numZeroBs << " blocks, " << _numZeroPs << " patches, " << _numZeroAs << " arcs to collapse";
    if (_numZeroAs == 0)
        return SUCCESS;

    _nTetsPre = (int)tetMesh.n_logical_cells();
    _nFsPre = (int)tetMesh.n_logical_faces();
    _nEsPre = (int)tetMesh.n_logical_edges();
    _nVsPre = (int)tetMesh.n_logical_vertices();

    _nBsPre = (int)mcMesh.n_logical_cells();
    _nPsPre = (int)mcMesh.n_logical_faces();
    _nAsPre = (int)mcMesh.n_logical_edges();
    _nNsPre = (int)mcMesh.n_logical_vertices();

    size_t initial = tetMesh.n_logical_cells();

    std::chrono::nanoseconds totalTime{};
    std::chrono::nanoseconds decimationTime{};
    auto start_time = std::chrono::high_resolution_clock::now();

    if (_numZeroAs > 0)
        assignCollapseDirs();

    bool newCollapseDir = true;
    bool change = true;

    while (change || newCollapseDir)
    {
        if (change)
            newCollapseDir = false;

        // if (tetMesh.n_logical_cells() > 1.0 * initial)
        if (change)
        {
            auto delta = std::chrono::high_resolution_clock::now() - start_time;
            totalTime += delta;
            start_time = std::chrono::high_resolution_clock::now();
            if (optimize)
                remesher.collapseAllPossibleEdges(false, true, true, true, 10.0);
            else
                remesher.collapseAllPossibleEdges(true, true, false, false);
            delta = std::chrono::high_resolution_clock::now() - start_time;
            totalTime += delta;
            decimationTime += delta;
            start_time = std::chrono::high_resolution_clock::now();
        }
        change = false;

        if (optimize && tetMesh.n_logical_cells() > 1.3 * initial)
        {
            auto delta = std::chrono::high_resolution_clock::now() - start_time;
            totalTime += delta;
            start_time = std::chrono::high_resolution_clock::now();
            remesher.remeshToImproveAngles(true, false, TetRemesher::QualityMeasure::ANGLES);
            initial = tetMesh.n_logical_cells();
            delta = std::chrono::high_resolution_clock::now() - start_time;
            totalTime += delta;
            decimationTime += delta;
            start_time = std::chrono::high_resolution_clock::now();
        }

        if (collapseNextPillowBlock() || collapseNextPillowPatch() || bisectNextBlockByPillowPatch()
            || collapseNextZeroArc() || bisectNextAlmostPillowPatch(false) || bisectNextAlmostPillowBlock())
        {
            change = true;
            continue;
        }

        if (!newCollapseDir && hasZeroLengthArcs())
        {
            assignCollapseDirs();
            newCollapseDir = true;
        }
        else if (bisectNextZeroPatch(true))
            change = true;
        else
            newCollapseDir = false;
    }

    if (zeroElementsRemain())
        return COLLAPSE_ERROR;

    reducer.init(true, true, true);
    while (reducer.isReducible())
        reducer.removeNextPatch();

    auto delta = std::chrono::high_resolution_clock::now() - start_time;
    totalTime += delta;
    start_time = std::chrono::high_resolution_clock::now();
    assertValidMC(true, true);
    delta = std::chrono::high_resolution_clock::now() - start_time;
    totalTime += delta;
    decimationTime += delta;

    auto msTotal = std::chrono::duration_cast<std::chrono::milliseconds>(totalTime);
    auto msDecimation = std::chrono::duration_cast<std::chrono::milliseconds>(decimationTime);

    LOG(INFO) << "POST-COLLAPSE stats: " << _nTetsPre << "; " << (int)tetMesh.n_logical_cells() << "; " << _nFsPre
              << "; " << (int)tetMesh.n_logical_faces() << "; " << _nEsPre << "; " << (int)tetMesh.n_logical_edges()
              << "; " << _nVsPre << "; " << (int)tetMesh.n_logical_vertices() << "; " << _nBsPre << "; "
              << (int)mcMesh.n_logical_cells() << "; " << _nPsPre << "; " << (int)mcMesh.n_logical_faces() << "; "
              << _nAsPre << "; " << (int)mcMesh.n_logical_edges() << "; " << _nNsPre << "; "
              << (int)mcMesh.n_logical_vertices() << "; " << _numZeroBs << "; " << _numZeroPs << "; " << _numZeroAs
              << "; " << _nCollapsedBs << "; " << _nCollapsedPs << "; " << _nCollapsedAs << "; "
              << _refiner.nBisectedBlocks() << "; " << _refiner.nBisectedPatches() << "; " << _refiner.nBisectedArcs()
              << "; " << msTotal.count() << "; " << msDecimation.count() << "; " << 0;

    start_time = std::chrono::high_resolution_clock::now();
    if (optimize)
    {
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(true, true, true, true, 10.0);
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 10.0);
        assertValidMC(true, true);
        remesher.remeshToImproveAngles(true, false, TetRemesher::QualityMeasure::ANGLES);
        assertValidMC(true, true);
        start_time = std::chrono::high_resolution_clock::now();
        MCSmoother(meshProps()).smoothMC();
        auto optimizationTime = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_time);
        meshProps().allocate<TOUCHED>(true);
        remesher.collapseAllPossibleEdges(false, true, true, true, 20);
        remesher.remeshToImproveAngles(true, false, TetRemesher::QualityMeasure::ANGLES);
        LOG(INFO) << "POST-COLLAPSE stats2: " << _nTetsPre << "; " << (int)tetMesh.n_logical_cells() << "; " << _nFsPre
                  << "; " << (int)tetMesh.n_logical_faces() << "; " << _nEsPre << "; " << (int)tetMesh.n_logical_edges()
                  << "; " << _nVsPre << "; " << (int)tetMesh.n_logical_vertices() << "; " << _nBsPre << "; "
                  << (int)mcMesh.n_logical_cells() << "; " << _nPsPre << "; " << (int)mcMesh.n_logical_faces() << "; "
                  << _nAsPre << "; " << (int)mcMesh.n_logical_edges() << "; " << _nNsPre << "; "
                  << (int)mcMesh.n_logical_vertices() << "; " << _numZeroBs << "; " << _numZeroPs << "; " << _numZeroAs
                  << "; " << _nCollapsedBs << "; " << _nCollapsedPs << "; " << _nCollapsedAs << "; "
                  << _refiner.nBisectedBlocks() << "; " << _refiner.nBisectedPatches() << "; "
                  << _refiner.nBisectedArcs() << "; " << msTotal.count() << "; " << msDecimation.count() << "; "
                  << optimizationTime.count();
    }

    assertValidMC(true, true);
    return SUCCESS;
}

MCCollapser::RetCode MCCollapser::collapseHalfarc(const HEH& haCollapse)
{
    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS, CHILD_FACES, CHILD_EDGES, CHILD_HALFFACES, CHILD_HALFEDGES>
        propGuard(meshProps());

    LOG(INFO) << "Collapsing halfarc " << haCollapse;
    EmbeddingCollapser(meshProps()).collapseArcEmbedding(haCollapse);
    collapseArcConnectivity(haCollapse);

    assertValidMC(false, false);
    return SUCCESS;
}

MCCollapser::RetCode MCCollapser::collapsePillowPatch(const FH& pCollapse, const set<HEH>& has)
{
    LOG(INFO) << "Collapsing pillow patch " << pCollapse;
    HEH haMoving = *has.begin();
    HEH haStationary = *(++has.begin());

    determineStationaryEnd(pCollapse, haMoving, haStationary);
    auto ret = EmbeddingCollapser(meshProps()).collapsePillowPatchEmbedding(pCollapse, haMoving, haStationary);
    if (ret != EmbeddingCollapser::SUCCESS)
        throw std::logic_error("Rerouting failed");
    collapsePillowPatchConnectivity(pCollapse, haMoving, haStationary);

    assertValidMC(false, false);
    return SUCCESS;
}

MCCollapser::RetCode MCCollapser::collapsePillowBlock(const CH& bCollapse)
{
    LOG(INFO) << "Collapsing minimal pillow block " << bCollapse;

    set<HFH> hps;
    for (HFH hp : mcMeshProps().mesh().cell_halffaces(bCollapse))
        hps.insert(hp);

    auto itHps = hps.begin();
    HFH hpStationary = *itHps;
    HFH hpMoving = *(++itHps);

    determineStationaryEnd(hpMoving, hpStationary);
    EmbeddingCollapser(meshProps()).collapsePillowBlockEmbedding(bCollapse, hpMoving, hpStationary);
    collapsePillowBlockConnectivity(bCollapse, hpMoving, hpStationary);

    assertValidMC(false, false);
    return SUCCESS;
}

MCCollapser::RetCode MCCollapser::collapseCigarBlock(const CH& bCollapse)
{
    LOG(INFO) << "Collapsing cigar block " << bCollapse;

    HFH hpCigar = *mcMeshProps().mesh().chf_iter(bCollapse);
    EmbeddingCollapser(meshProps()).collapseCigarBlockEmbedding(bCollapse);
    collapsePillowBlockConnectivity(bCollapse, hpCigar, HFH());

    assertValidMC(false, false);
    return SUCCESS;
}

void MCCollapser::collapseArcConnectivity(const HEH& haCollapse)
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    EH aCollapse = mcMesh.edge_handle(haCollapse);

    VH nFrom = mcMesh.from_vertex_handle(haCollapse);
    VH nTo = mcMesh.to_vertex_handle(haCollapse);

    vector<EH> asOnFrom;
    for (EH a : mcMesh.vertex_edges(nFrom))
        asOnFrom.emplace_back(a);

    set<CH> bsOnFrom;
    for (CH b : mcMesh.vertex_cells(nFrom))
        bsOnFrom.insert(b);
    set<CH> bsOnTo;
    for (CH b : mcMesh.vertex_cells(nTo))
        bsOnTo.insert(b);

    set<CH> bsOnCollapse;
    for (FH p : mcMesh.edge_faces(aCollapse))
        for (CH b : mcMesh.face_cells(p))
            if (b.is_valid())
                bsOnCollapse.insert(b);

    DLOG(INFO) << "DELETING NODE " << nFrom << " AND REPLACING BY " << nTo;

    // UPDATE REFERENCES TO NODES IN BLOCKS
    for (CH b : bsOnFrom)
    {
        bool containsNTo = contains(mcMesh.cell_vertices(b), nTo);
        bool oneIsCorner = containsMatching(mcMeshProps().ref<BLOCK_CORNER_NODES>(b),
                                            [&](const pair<const UVWDir, VH>& kv)
                                            { return kv.second == nTo || kv.second == nFrom; });

        for (auto& dir2n : mcMeshProps().ref<BLOCK_CORNER_NODES>(b))
            if (dir2n.second == nFrom)
            {
                dir2n.second = nTo;
                if (containsNTo)
                {
                    for (UVWDir dirEdge : DIM_2_DIRS)
                        mcMeshProps().ref<BLOCK_EDGE_NODES>(b).at(dirEdge).erase(nTo);
                    for (UVWDir dirFace : DIM_1_DIRS)
                        mcMeshProps().ref<BLOCK_FACE_NODES>(b).at(dirFace).erase(nTo);
                }
            }

        for (auto& dir2ns : mcMeshProps().ref<BLOCK_EDGE_NODES>(b))
        {
            auto it = dir2ns.second.find(nFrom);
            if (it != dir2ns.second.end())
            {
                dir2ns.second.erase(it);
                if (!containsNTo)
                    dir2ns.second.insert(nTo);
                else if (mcMeshProps().ref<BLOCK_EDGE_ARCS>(b).at(dir2ns.first).count(aCollapse) != 0 || oneIsCorner)
                {
                    // Do nothing
                }
                else
                {
                    dir2ns.second.insert(nTo);
                    auto dirsFace = decompose(dir2ns.first, DIM_1_DIRS);
                    for (int i = 0; i != 2; i++)
                    {
                        UVWDir dirFace = dirsFace[i];
                        if (mcMeshProps().ref<BLOCK_FACE_ARCS>(b).at(dirFace).count(aCollapse) != 0)
                        {
                            auto& fns = mcMeshProps().ref<BLOCK_FACE_NODES>(b).at(dirFace);
                            auto itf = fns.find(nTo);

                            if (itf != fns.end())
                                fns.erase(itf);
                        }
                    }
                    assert(dir2ns.second.count(nTo) != 0);
                }
            }
        }
        // Test if that properly works for all models
        bool edgesHaveTo
            = containsMatching(mcMeshProps().ref<BLOCK_EDGE_NODES>(b),
                               [&](const pair<const UVWDir, set<VH>>& kv) { return kv.second.count(nTo) != 0; });
        for (auto& dir2ns : mcMeshProps().ref<BLOCK_FACE_NODES>(b))
        {
            auto it = dir2ns.second.find(nFrom);
            if (it != dir2ns.second.end())
            {
                dir2ns.second.erase(it);
                if (!containsNTo && !edgesHaveTo)
                    dir2ns.second.insert(nTo);
            }
            if (edgesHaveTo)
            {
                it = dir2ns.second.find(nTo);
                if (it != dir2ns.second.end())
                    dir2ns.second.erase(it);
            }
        }
    }

    // Reconnect arcs incident on nFrom to nTo
    for (EH a : asOnFrom)
    {
        auto vs = mcMesh.edge_vertices(a);
        if (vs[0] == nFrom)
            vs[0] = nTo;
        if (vs[1] == nFrom)
            vs[1] = nTo;
        mcMesh.set_edge(a, vs[0], vs[1]);
    }

    deferredDeleteNode(nFrom);

    // DELETE REFERENCES TO ARC IN BLOCKS
    set<FH> psOnA;
    set<CH> bsOnA;
    for (FH p : mcMesh.edge_faces(aCollapse))
        psOnA.insert(p);
    for (FH p : psOnA)
        for (CH b : mcMesh.face_cells(p))
            if (b.is_valid())
                bsOnA.insert(b);

    DLOG(INFO) << "DELETING ARC " << aCollapse;
    replaceArcIncidentPatches({{haCollapse, {}}, {tetMesh.opposite_halfedge_handle(haCollapse), {}}}, psOnA);
    map<EH, vector<EH>> aReplacements({{aCollapse, {}}});
    updateBlockArcReferences(aReplacements, bsOnA);
    deferredDeleteArc(aCollapse);

    _nCollapsedAs++;
}

void MCCollapser::collapsePillowPatchConnectivity(const FH& p, const HEH& haRemoved, const HEH& haRemaining)
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    EH aRemaining = mcMesh.edge_handle(haRemaining);
    EH aRemoved = mcMesh.edge_handle(haRemoved);

    DLOG(INFO) << "DELETING PATCH " << p;
    DLOG(INFO) << "REPLACING ARC " << aRemoved << " BY " << aRemaining;
    DLOG(INFO) << "AND HALFARC " << mcMesh.opposite_halfedge_handle(haRemoved) << " BY " << haRemaining;
    DLOG(INFO) << "AND HALFARC " << haRemoved << " BY " << mcMesh.opposite_halfedge_handle(haRemaining);

    set<CH> bsOnP;
    for (CH b : mcMesh.face_cells(p))
        if (b.is_valid())
            bsOnP.insert(b);

    set<FH> psOnA2;
    set<CH> bsOnA2;
    for (FH p2 : mcMesh.edge_faces(aRemoved))
        psOnA2.insert(p2);
    for (FH p2 : psOnA2)
        for (CH b : mcMesh.face_cells(p2))
            if (b.is_valid())
                bsOnA2.insert(b);

    bool flipA1AfterReplacement = (haRemaining.idx() % 2) == (haRemoved.idx() % 2);
    map<CH, UVWDir> b2dirA2;
    if (flipA1AfterReplacement)
        for (CH b : bsOnA2)
            b2dirA2[b] = halfarcDirInBlock(mcMesh.halfedge_handle(aRemoved, 0), b);

    for (CH b : bsOnA2)
    {
        // bool containsP = bsOnP.count(b) != 0;
        bool containsRemaining = contains(mcMesh.cell_edges(b), aRemaining);

        for (auto& dir2as : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b))
        {
            auto it = dir2as.second.find(aRemoved);
            if (it != dir2as.second.end())
            {
                dir2as.second.erase(it);
                dir2as.second.insert(aRemaining);
                if (containsRemaining)
                {
                    for (UVWDir dirFace : decompose(dir2as.first, DIM_1_DIRS))
                    {
                        auto& fas = mcMeshProps().ref<BLOCK_FACE_ARCS>(b).at(dirFace);
                        auto itf = fas.find(aRemaining);
                        if (itf != fas.end())
                            fas.erase(itf);
                    }
                }
            }
        }

        for (auto& dir2as : mcMeshProps().ref<BLOCK_FACE_ARCS>(b))
        {
            auto it = dir2as.second.find(aRemoved);
            if (it != dir2as.second.end())
            {
                dir2as.second.erase(it);
                if (!containsRemaining)
                    dir2as.second.insert(aRemaining);
            }
        }
        UVWDir dirA2 = halfarcDirInBlock(mcMesh.halfedge_handle(aRemoved, 0), b);
        assert(dirA2 != UVWDir::NONE);
        bool hasA1 = false;
        for (auto& kv : mcMeshProps().ref<BLOCK_ALL_ARCS>(b))
        {
            if (kv.second.count(aRemoved) != 0)
                kv.second.erase(aRemoved);
            if (kv.second.count(aRemaining) != 0)
                hasA1 = true;
        }
        if (!hasA1)
            mcMeshProps().ref<BLOCK_ALL_ARCS>(b).at(flipA1AfterReplacement ? -dirA2 : dirA2).insert(aRemaining);
    }

    replaceArcIncidentPatches({{haRemoved, {tetMesh.opposite_halfedge_handle(haRemaining)}},
                               {tetMesh.opposite_halfedge_handle(haRemoved), {haRemaining}}},
                              psOnA2);
    map<EH, EH> aRemovedReplacements({{aRemoved, aRemaining}});

    int singularityIdx = mcMeshProps().get<IS_SINGULAR>(aRemaining) + mcMeshProps().get<IS_SINGULAR>(aRemoved);
    mcMeshProps().set<IS_SINGULAR>(aRemaining, singularityIdx);

    replacePatchIncidentBlocks({{tetMesh.halfface_handle(p, 0), {}}, {tetMesh.halfface_handle(p, 1), {}}}, bsOnP);
    map<FH, vector<FH>> pReplacements({{p, {}}});
    updateBlockPatchReferences(pReplacements, bsOnP);

    deferredDeletePatch(p);
    deferredDeleteArc(aRemoved);

    _nCollapsedPs++;
}

void MCCollapser::collapsePillowBlockConnectivity(const CH& b, const HFH& hpRemoved, const HFH& hpRemaining)
{
    auto& mcMesh = mcMeshProps().mesh();

    FH pRemaining = hpRemaining.is_valid() ? mcMesh.face_handle(hpRemaining) : FH(-1);
    FH pRemoved = mcMesh.face_handle(hpRemoved);

    vector<HEH> has;
    if (hpRemaining.is_valid())
        for (HEH ha : mcMesh.face_halfedges(pRemaining))
            has.push_back(ha);

    set<CH> bsOnP2;
    for (CH b2 : mcMesh.face_cells(pRemoved))
        if (b2.is_valid() && b2 != b)
            bsOnP2.insert(b2);

    if (!hpRemaining.is_valid())
    {
        CH b2 = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hpRemoved));
        vector<EH> as;
        for (EH a : mcMesh.halfface_edges(hpRemoved))
            as.push_back(a);

        for (EH a : as)
        {
            bool inEdgeArcs = false;
            for (auto& dir2as : mcMeshProps().ref<BLOCK_EDGE_ARCS>(b2))
                inEdgeArcs = inEdgeArcs || dir2as.second.count(a) != 0;

            if (inEdgeArcs)
                for (auto& dir2as : mcMeshProps().ref<BLOCK_FACE_ARCS>(b2))
                    dir2as.second.erase(a);
        }
    }

    DLOG(INFO) << "DELETING BLOCK " << b;
    deferredDeleteBlock(b);

    map<HFH, vector<HFH>> hpReplacements;
    map<FH, vector<FH>> pReplacements;
    if (hpRemaining.is_valid())
    {
        hpReplacements = {{hpRemoved, {mcMesh.opposite_halfface_handle(hpRemaining)}},
                          {mcMesh.opposite_halfface_handle(hpRemoved), {hpRemaining}}};
        pReplacements = {{pRemoved, {pRemaining}}};
    }
    else
    {
        hpReplacements = {{hpRemoved, {}}, {mcMesh.opposite_halfface_handle(hpRemoved), {}}};
        pReplacements = {{pRemoved, {}}};
    }
    replacePatchIncidentBlocks(hpReplacements, bsOnP2);
    updateBlockPatchReferences(pReplacements, bsOnP2);

    if (hpRemaining.is_valid() && mcMeshProps().isAllocated<PATCH_MIN_DIST>())
        mcMeshProps().set<PATCH_MIN_DIST>(
            pRemoved,
            std::min(mcMeshProps().get<PATCH_MIN_DIST>(pRemaining), mcMeshProps().get<PATCH_MIN_DIST>(pRemoved)));

    if (hpRemaining.is_valid())
        mcMesh.set_face(pRemaining, has);

    DLOG(INFO) << "REPLACING PATCH " << pRemoved << " BY " << pRemaining;
    deferredDeletePatch(pRemoved);

    _nCollapsedBs++;
}

bool MCCollapser::collapseNextZeroArc()
{
    auto& mcMesh = mcMeshProps().mesh();

    DLOG(INFO) << "Looking for 0-arc to collapse";

    vector<EH> asZero;
    for (EH a : mcMesh.edges())
        if (isZeroArc(a))
            asZero.push_back(a);
    if (_randomOrder)
        std::random_shuffle(asZero.begin(), asZero.end());
    else
        std::sort(asZero.begin(),
                  asZero.end(),
                  [&](const EH& a, const EH& b) -> bool
                  { return mcMeshProps().get<ARC_DBL_LENGTH>(a) < mcMeshProps().get<ARC_DBL_LENGTH>(b); });
    if (!asZero.empty())
        assert(mcMeshProps().get<ARC_DBL_LENGTH>(asZero.front()) <= mcMeshProps().get<ARC_DBL_LENGTH>(asZero.back()));
    for (EH a : asZero)
    {
        HEH haPreferred = preferredCollapseHalfarc(a);
        if (!haPreferred.is_valid())
            continue;

        collapseHalfarc(haPreferred);
        return true;
    }
    return false;
}

bool MCCollapser::collapseNextPillowPatch()
{
    auto& mcMesh = mcMeshProps().mesh();
    TemporaryPropAllocator<TetMeshProps, CHILD_CELLS, CHILD_FACES, CHILD_EDGES, CHILD_HALFFACES, CHILD_HALFEDGES>
        propGuard(meshProps());

    DLOG(INFO) << "Looking for pillow patch to collapse";

    for (FH p : mcMesh.faces())
    {
        HFH hp = mcMesh.halfface_handle(p, 0);
        set<HEH> has;
        for (HEH ha : mcMesh.halfface_halfedges(hp))
            has.insert(ha);
        if (has.size() == 2 && *(has.begin()) != mcMesh.opposite_halfedge_handle(*(++has.begin())))
        {
            HEH ha1 = *has.begin();
            HEH ha2 = *(++has.begin());
            determineStationaryEnd(p, ha1, ha2);
            if (ha1.is_valid())
            {
                collapsePillowPatch(p, has);
                return true;
            }
            else
            {
                LOG(WARNING) << "Patch " << p << " skipped (this might be an algorithm error)";
            }
        }
    }
    return false;
}

bool MCCollapser::collapseNextPillowBlock()
{
    auto& mcMesh = mcMeshProps().mesh();

    DLOG(INFO) << "Looking for pillow block to collapse";

    for (CH b : mcMesh.cells())
    {
        set<HFH> hps;
        for (HFH hp : mcMesh.cell_halffaces(b))
            hps.insert(hp);
        if (hps.size() == 1)
        {
            collapseCigarBlock(b);
            return true;
        }
        else if (hps.size() == 2)
        {
            auto itHps = hps.begin();
            HFH hp1 = *itHps;
            HFH hp2 = *(++itHps);
            FH p1 = mcMesh.face_handle(hp1);
            FH p2 = mcMesh.face_handle(hp2);

            set<EH> as1, as2;
            for (EH a : mcMesh.face_edges(p1))
                as1.insert(a);
            for (EH a : mcMesh.face_edges(p2))
                as2.insert(a);
            // IN AN ALMOST-CIGAR AS1 != AS2, CAN NOT MERGE SUCH THING
            if (as1 != as2)
                continue;

            determineStationaryEnd(hp1, hp2);
            if (hp1.is_valid())
            {
                collapsePillowBlock(b);
                return true;
            }
            else
                LOG(WARNING) << "Block " << b << " skipped (this might be an algorithm error)";
        }
    }
    return false;
}

bool MCCollapser::bisectNextAlmostPillowPatch(bool allowZeroLoop)
{
    auto& mcMesh = mcMeshProps().mesh();

    DLOG(INFO) << "Looking for quasi-pillow-patches to bisect";

    auto quasiPillowPatches = findAlmostPillowPatches();

    for (FH p : quasiPillowPatches)
    {
        auto itPair = mcMesh.halfface_halfedges(mcMesh.halfface_handle(p, 0));
        size_t numHas = std::distance(itPair.first, itPair.second);
        vector<FH> psSub;
        if (numHas > 3 && _refiner.bisectPatchAcrossDir(p, UVWDir::ANY, allowZeroLoop, psSub))
        {
            LOG(INFO) << "Bisected patch " << p;
            _nBisectionsP++;
            assertValidMC(false, false);
            return true;
        }
    }
    return false;
}

bool MCCollapser::bisectNextZeroPatch(bool allowZeroLoop)
{
    auto& mcMesh = mcMeshProps().mesh();

    DLOG(INFO) << "Looking for 0-patches to bisect";

    for (FH p : mcMesh.faces())
    {
        HFH hp = mcMesh.halfface_handle(p, 0);
        if (mcMesh.is_boundary(hp))
            hp = mcMesh.opposite_halfface_handle(hp);
        auto dirNormal = halfpatchNormalDir(hp);
        UVWDir dirsHp = ~(dirNormal | -dirNormal);
        auto dir2has = halfpatchHalfarcsByDir(hp);
        UVWDir zeroDir = UVWDir::NONE;
        for (UVWDir dir : decompose(dirsHp, DIM_1_DIRS))
            if (dir2has.count(dir) == 0)
            {
                zeroDir = dir | -dir;
                break;
            }
        if (zeroDir == UVWDir::NONE)
        {
            for (auto& kv : dir2has)
            {
                int length = 0;
                for (HEH ha : kv.second)
                    length += mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha));
                if (length == 0)
                {
                    zeroDir = kv.first | -kv.first;
                    break;
                }
            }
        }
        vector<FH> psSub;
        if (zeroDir != UVWDir::NONE && _refiner.bisectPatchAcrossDir(p, zeroDir, allowZeroLoop, psSub))
        {
            LOG(INFO) << "Bisected patch " << p;
            _nBisectionsP++;
            assertValidMC(false, false);
            return true;
        }
    }
    return false;
}

bool MCCollapser::bisectNextAlmostPillowBlock()
{
    auto& mcMesh = mcMeshProps().mesh();

    DLOG(INFO) << "Looking for 0-blocks to bisect";
    auto pillowBlocks = findAlmostPillowBlocks();
    if (!pillowBlocks.empty())
    {
        for (CH b : pillowBlocks)
        {
            set<HFH> hps;
            for (HFH hp : mcMesh.cell_halffaces(b))
                hps.insert(hp);
            if (hps.size() <= 2)
                continue;
            vector<CH> subBlocks;
            bool refined = _refiner.bisectBlockOrPatch(b, subBlocks);
            if (refined)
            {
                LOG(INFO) << "Bisected block " << b;
                _nBisectionsB++;
                assertValidMC(false, false);
                return true;
            }
        }
    }
    return false;
}

set<FH> MCCollapser::findAlmostPillowPatches() const
{
    auto& mcMesh = mcMeshProps().mesh();
    set<FH> quasiPillowPatches;

    for (FH p : mcMesh.faces())
    {
        bool hasZero = false;
        auto hasByDir = halfpatchHalfarcsByDir(mcMesh.halfface_handle(p, 0));
        if (hasByDir.size() == 2)
        {
            hasZero = containsMatching(hasByDir,
                                       [&](const pair<const UVWDir, vector<HEH>>& kv)
                                       {
                                           return containsMatching(kv.second,
                                                                   [&, this](const HEH& ha)
                                                                   { return isZeroArc(mcMesh.edge_handle(ha)); });
                                       });
            if (!hasZero)
                quasiPillowPatches.insert(p);
        }
    }
    return quasiPillowPatches;
}

bool MCCollapser::bisectNextBlockByPillowPatch()
{
    auto& mcMesh = mcMeshProps().mesh();
    auto& tetMesh = meshProps().mesh();

    DLOG(INFO) << "Looking for 2-zero-arc-cycle to bisect";
    for (CH b : mcMesh.cells())
        for (HEH ha : mcMesh.cell_halfedges(b))
        {
            if (!isZeroArc(mcMesh.edge_handle(ha)))
                continue;
            for (HEH ha2 : mcMesh.cell_halfedges(b))
            {
                if (mcMesh.edge_handle(ha) != mcMesh.edge_handle(ha2)
                    && mcMesh.from_vertex_handle(ha) == mcMesh.to_vertex_handle(ha2)
                    && mcMesh.to_vertex_handle(ha) == mcMesh.from_vertex_handle(ha2))
                {
                    if (mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha))
                        != mcMeshProps().get<ARC_INT_LENGTH>(mcMesh.edge_handle(ha2)))
                        continue;
                    bool hasP = containsMatching(mcMesh.cell_faces(b),
                                                 [&, this](const FH& p)
                                                 {
                                                     auto itPair = mcMesh.face_edges(p);
                                                     return std::distance(itPair.first, itPair.second) == 2
                                                            && contains(itPair, mcMesh.edge_handle(ha))
                                                            && contains(itPair, mcMesh.edge_handle(ha2));
                                                 });
                    if (!hasP)
                    {
                        vector<HEH> ringHas = {ha, ha2};
                        set<EH> ringAs;
                        for (HEH ha3 : ringHas)
                            ringAs.insert(mcMesh.edge_handle(ha3));
                        HFH seedHp = findMatching(mcMesh.halfedge_halffaces(ringHas.front()),
                                                  [&](const HFH& hp) { return mcMesh.incident_cell(hp) == b; });
                        list<HFH> hpQ({seedHp});
                        set<HFH> sideHps({seedHp});
                        while (!hpQ.empty())
                        {
                            HFH hp = hpQ.front();
                            hpQ.pop_front();
                            for (HEH ha3 : mcMesh.halfface_halfedges(hp))
                            {
                                if (ringAs.count(mcMesh.edge_handle(ha3)) != 0)
                                    continue;
                                HFH hpNext = safeAdjacentHalffaceInBlock(hp, ha3);
                                if (!hpNext.is_valid() || sideHps.find(hpNext) != sideHps.end())
                                    continue;
                                hpQ.emplace_back(hpNext);
                                sideHps.insert(hpNext);
                            }
                        }
                        set<HFH> hpsCheck;
                        for (HFH hp : mcMesh.cell_halffaces(b))
                            hpsCheck.insert(hp);

                        if (hpsCheck.size() > sideHps.size())
                        {
                            // For safety against really weird MC connectivity: see if ha, ha2 edges also partition the
                            // blocks boundary faces into two
                            vector<HFH> hfs;
                            for (HFH hp : hpsCheck)
                            {
                                auto hfsHp = mcMeshProps().hpHalffaces(hp);
                                hfs.insert(hfs.end(), hfsHp.begin(), hfsHp.end());
                            }
                            set<EH> ringEs;
                            for (EH a : ringAs)
                                for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(a))
                                    ringEs.insert(tetMesh.edge_handle(he));
                            HFH hfSeed = hfs.front();
                            set<HFH> hfsVisited({hfSeed});
                            list<HFH> hfQ({hfSeed});
                            while (!hfQ.empty())
                            {
                                HFH hf = hfQ.front();
                                hfQ.pop_front();
                                for (HEH he : tetMesh.halfface_halfedges(hf))
                                {
                                    if (ringEs.count(tetMesh.edge_handle(he)) != 0)
                                        continue;
                                    HFH hfNext = adjacentHfOnWall(hf, he);
                                    if (!hfNext.is_valid() || hfsVisited.find(hfNext) != hfsVisited.end())
                                        continue;
                                    hfQ.emplace_back(hfNext);
                                    hfsVisited.insert(hfNext);
                                }
                            }
                            if (hfsVisited.size() < hfs.size())
                            {
                                std::stringstream sstr1;
                                for (HFH hp : mcMesh.cell_halffaces(b))
                                {
                                    sstr1 << hp << " (";
                                    for (HEH ha3 : mcMesh.halfface_halfedges(hp))
                                        sstr1 << ha3 << ", ";
                                    sstr1 << "), ";
                                }
                                LOG(INFO) << "Bisecting block " << b << " which consists of halfpatches " << sstr1.str()
                                          << " by a pillow patch";
                                std::stringstream sstr2;
                                for (HEH ringHa : ringHas)
                                    sstr2 << ringHa << ", ";
                                LOG(INFO) << "...which has the following halfarcs on its boundary: " << sstr2.str();
                                vector<CH> subBlocks;
                                _refiner.cutBlock(b, ringHas, subBlocks);
                                _nBisectionsP++;
                                assertValidMC(false, false);
                                return true;
                            }
                        }
                    }
                }
            }
        }
    return false;
}

set<CH> MCCollapser::findAlmostPillowBlocks() const
{
    auto& mcMesh = mcMeshProps().mesh();
    set<CH> pillowBlocks;

    for (CH b : mcMesh.cells())
    {
        auto& dir2ps = mcMeshProps().ref<BLOCK_FACE_PATCHES>(b);
        int nFaces = 0;
        int nPatches = 0;
        for (auto& kv : dir2ps)
            if (!kv.second.empty())
            {
                nPatches += (int)kv.second.size();
                nFaces++;
            }
        if (nFaces == 2)
            pillowBlocks.insert(b);
    }
    return pillowBlocks;
}

void MCCollapser::assignCollapseDirs()
{
    if (_direction != 0)
        return;
    auto& mcMesh = mcMeshProps().mesh();

    if (mcMeshProps().isAllocated<BLOCK_COLLAPSE_DIR>())
        mcMeshProps().release<BLOCK_COLLAPSE_DIR>();

    // find the shortest zero arc and floodfill collapsedir into all blocks and onwards
    vector<EH> asZero;
    for (EH a : mcMesh.edges())
        if (isZeroArc(a))
            asZero.push_back(a);
    std::sort(asZero.begin(),
              asZero.end(),
              [&](const EH& a, const EH& b) -> bool
              { return mcMeshProps().get<ARC_DBL_LENGTH>(a) < mcMeshProps().get<ARC_DBL_LENGTH>(b); });

    EH aZero = findMatching(
        asZero,
        [&](const EH& a)
        { return !collapseIsLocked(mcMesh.halfedge_handle(a, 0)) || !collapseIsLocked(mcMesh.halfedge_handle(a, 1)); });

    if (aZero.is_valid())
    {
        HEH haZero = preferredCollapseHalfarc(aZero);
        list<pair<CH, UVWDir>> bQ;
        mcMeshProps().allocate<BLOCK_COLLAPSE_DIR>(UVWDir::NONE);
        set<CH> bs;
        for (FH p : mcMesh.edge_faces(aZero))
            for (CH b : mcMesh.face_cells(p))
                if (b.is_valid())
                    bs.insert(b);
        for (CH b : bs)
        {
            UVWDir dir = halfarcDirInBlock(haZero, b);
            mcMeshProps().set<BLOCK_COLLAPSE_DIR>(b, dir);
            bQ.push_back({b, dir});
        }
        while (!bQ.empty())
        {
            auto bDir = bQ.front();
            bQ.pop_front();
            CH b = bDir.first;
            UVWDir dirCollapse = bDir.second;

            for (HFH hp : mcMesh.cell_halffaces(b))
            {
                CH bNext = mcMesh.incident_cell(mcMesh.opposite_halfface_handle(hp));
                if (bNext.is_valid() && mcMeshProps().get<BLOCK_COLLAPSE_DIR>(b) == UVWDir::NONE)
                {
                    UVWDir dirNext = mcMeshProps().hpTransition<PATCH_TRANSITION>(hp).rotate(dirCollapse);
                    auto dirs2as = mcMeshProps().get<BLOCK_ALL_ARCS>(b);
                    bool hasZeroInDirNext
                        = containsMatching(mcMesh.cell_edges(bNext),
                                           [&, this](const EH& a)
                                           {
                                               return isZeroArc(a)
                                                      && (halfarcDirInBlock(mcMesh.halfedge_handle(a, 0), b)
                                                          & (dirNext | -dirNext))
                                                             != UVWDir::NONE;
                                           })
                          || containsMatching(mcMesh.cell_halffaces(b),
                                              [&, this](const HFH& hp2)
                                              {
                                                  UVWDir dirNormal = halfpatchNormalDir(hp2);
                                                  auto dir2has = halfpatchHalfarcsByDir(hp2);
                                                  return containsMatching(
                                                      decompose(~(dirNormal | -dirNormal), DIM_1_DIRS),
                                                      [&](const UVWDir& dir)
                                                      {
                                                          return dir2has.count(dir) == 0
                                                                 && (dir & (dirNext | -dirNext)) != UVWDir::NONE;
                                                      });
                                              })
                          || containsMatching((UVWDir[]){UVWDir::POS_U, UVWDir::POS_V, UVWDir::POS_W},
                                              [&](const UVWDir& dir)
                                              {
                                                  return dirs2as.count(dir) == 0 && dirs2as.count(-dir) == 0
                                                         && (dir & (dirNext | -dirNext)) != UVWDir::NONE;
                                              });
                    if (hasZeroInDirNext)
                    {
                        mcMeshProps().set<BLOCK_COLLAPSE_DIR>(bNext, dirNext);
                        bQ.push_back({bNext, dirNext});
                    }
                }
            }
        }
        int n = 0, nOpp = 0;

        for (HEH ha : mcMesh.halfedges())
        {
            if (isZeroArc(mcMesh.edge_handle(ha)))
            {
                CH b = *mcMesh.hec_iter(ha);
                UVWDir dir0 = halfarcDirInBlock(ha, b);
                UVWDir dirB = mcMeshProps().get<BLOCK_COLLAPSE_DIR>(b);
                if (dir0 == dirB && collapseIsLocked(ha))
                    n++;
                else if (-dir0 == dirB && collapseIsLocked(mcMesh.opposite_halfedge_handle(ha)))
                    nOpp++;
            }
        }
        if (n > nOpp)
            for (CH b : mcMesh.cells())
                mcMeshProps().set<BLOCK_COLLAPSE_DIR>(b, -mcMeshProps().get<BLOCK_COLLAPSE_DIR>(b));
    }
}

void MCCollapser::countZeroElements(int& numZeroAs, int& numZeroPs, int& numZeroBs)
{
    auto& mcMesh = mcMeshProps().mesh();

    numZeroAs = numZeroPs = numZeroBs = 0;

    for (EH a : mcMesh.edges())
        if (isZeroArc(a))
            numZeroAs++;
    for (FH p : mcMesh.faces())
        if (isZeroPatch(p))
            numZeroPs++;
    for (CH b : mcMesh.cells())
        if (isZeroBlock(b))
            numZeroBs++;
}

bool MCCollapser::zeroElementsRemain()
{
    int nAs, nPs, nBs;
    countZeroElements(nAs, nPs, nBs);

    return nAs > 0 || nPs > 0 || nBs > 0;
}

void MCCollapser::determineStationaryEnd(const FH& pCollapse, HEH& haMoving, HEH& haStationary) const
{
    auto& mcMesh = mcMeshProps().mesh();
    EH aStationary = mcMesh.edge_handle(haStationary);
    EH aMoving = mcMesh.edge_handle(haMoving);

    bool movingCritical = mcMeshProps().get<IS_CRITICAL_A>(aMoving)
                          || (containsMatching(mcMesh.edge_faces(aMoving),
                                               [&, this](const FH& p) { return mcMeshProps().get<IS_CRITICAL_P>(p); })
                              && !mcMeshProps().get<IS_CRITICAL_P>(pCollapse));
    bool stationaryCritical
        = mcMeshProps().get<IS_CRITICAL_A>(aStationary)
          || (containsMatching(mcMesh.edge_faces(aStationary),
                               [&, this](const FH& p) { return mcMeshProps().get<IS_CRITICAL_P>(p); })
              && !mcMeshProps().get<IS_CRITICAL_P>(pCollapse));
    if (movingCritical && stationaryCritical)
    {
        haMoving = HEH();
        haStationary = HEH();
        return;
    }
    if (stationaryCritical)
        return;
    if (movingCritical)
    {
        std::swap(haStationary, haMoving);
        return;
    }

    // Soft criterion: valences
    if (mcMesh.valence(aMoving) > mcMesh.valence(aStationary))
    {
        std::swap(haStationary, haMoving);
        return;
    }

    if (mcMesh.valence(aMoving) == mcMesh.valence(aStationary))
    {
        double lengthMoving = 0.0;
        double lengthStationary = 0.0;
        for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(aMoving))
            lengthMoving += edgeLengthUVW<CHART>(meshProps().mesh().edge_handle(he));
        for (HEH he : mcMeshProps().ref<ARC_MESH_HALFEDGES>(aStationary))
            lengthStationary += edgeLengthUVW<CHART>(meshProps().mesh().edge_handle(he));

        if (lengthStationary > lengthMoving)
            std::swap(haStationary, haMoving);
    }
}

void MCCollapser::determineStationaryEnd(HFH& hpMoving, HFH& hpStationary) const
{
    auto& mcMesh = mcMeshProps().mesh();

    FH pMoving = mcMesh.face_handle(hpMoving);
    FH pStationary = mcMesh.face_handle(hpStationary);

    bool movingCritical = mcMesh.is_boundary(pMoving) || mcMeshProps().get<IS_CRITICAL_P>(pMoving);
    bool stationaryCritical = mcMesh.is_boundary(pStationary) || mcMeshProps().get<IS_CRITICAL_P>(pStationary);

    if (movingCritical && stationaryCritical)
    {
        hpMoving = HFH();
        hpStationary = HFH();
        return;
    }

    if (movingCritical)
    {
        std::swap(hpStationary, hpMoving);
        return;
    }

    if (stationaryCritical)
        return;

    double areaMoving = 0.0;
    double areaStationary = 0.0;
    for (HFH hf : mcMeshProps().hpHalffaces(hpMoving))
        areaMoving += areaUVW<CHART>(meshProps().mesh().face_handle(hf));
    for (HFH hf : mcMeshProps().hpHalffaces(hpStationary))
        areaStationary += areaUVW<CHART>(meshProps().mesh().face_handle(hf));

    if (areaStationary > areaMoving)
        std::swap(hpStationary, hpMoving);
}

void MCCollapser::markZeros()
{
    mcMeshProps().allocate<MARK_N>(0);
    mcMeshProps().allocate<MARK_A>(0);
    mcMeshProps().allocate<MARK_P>(0);
    mcMeshProps().allocate<MARK_B>(0);

    for (EH a : mcMeshProps().mesh().edges())
        if (isZeroArc(a))
        {
            bool hasZeroPatch = false;
            for (FH p : mcMeshProps().mesh().edge_faces(a))
                if (isZeroPatch(p))
                {
                    hasZeroPatch = true;
                    break;
                }
            bool hasZeroBlock = false;

            set<CH> bs;
            for (FH p : mcMeshProps().mesh().edge_faces(a))
                for (CH b : mcMeshProps().mesh().face_cells(p))
                    if (b.is_valid())
                        bs.insert(b);
            for (CH b : bs)
                if (isZeroBlock(b))
                {
                    hasZeroBlock = true;
                    break;
                }
            if (!hasZeroPatch)
                mcMeshProps().set<MARK_A>(a, 1);
            else if (!hasZeroBlock)
                mcMeshProps().set<MARK_A>(a, 2);
            else
                mcMeshProps().set<MARK_A>(a, 3);
        }

    for (FH p : mcMeshProps().mesh().faces())
    {
        bool zeroPatch = isZeroPatch(p);
        bool hasZeroBlock = false;
        for (CH b : mcMeshProps().mesh().face_cells(p))
            if (b.is_valid() && isZeroBlock(b))
            {
                hasZeroBlock = true;
                break;
            }
        if (zeroPatch)
        {
            if (!hasZeroBlock)
                mcMeshProps().set<MARK_P>(p, 4);
            else
                mcMeshProps().set<MARK_P>(p, 5);
        }
        else if (hasZeroBlock)
            mcMeshProps().set<MARK_P>(p, 6);
    }
    for (CH b : mcMeshProps().mesh().cells())
        if (isZeroBlock(b))
            mcMeshProps().set<MARK_B>(b, 7);
}

} // namespace c4hex
