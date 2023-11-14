#include "C4Hex/Interface/HexRemesher.hpp"

#include "C4Hex/Algorithm/HexExtractor.hpp"

#include "OpenVolumeMesh/FileManager/FileManager.hh"

#include <fstream>
#include <iomanip>

namespace c4hex
{

HexRemesher::HexRemesher(const TetMeshProps& meshProps) : TetMeshNavigator(meshProps), MCMeshNavigator(meshProps)
{
}

HexRemesher::RetCode HexRemesher::extractHexMesh(HexMeshProps& hexMeshProps)
{
    HexExtractor<HexMeshProps> hexex(meshProps(), hexMeshProps);
    hexex.extractHexMesh(0);
    return SUCCESS;
}

HexRemesher::RetCode HexRemesher::extractPolyHexMesh(PolyMeshProps& hexMeshProps, int subdiv)
{
    HexExtractor<PolyMeshProps> hexex(meshProps(), hexMeshProps);
    hexex.extractHexMesh(subdiv);
    return SUCCESS;
}

HexRemesher::RetCode HexRemesher::extractMCMesh(PolyMeshProps& polyMeshProps)
{
    map<VH, VH> v2v;
    map<EH, EH> e2e;
    map<FH, FH> f2f;

    PolyMesh& polyMesh = polyMeshProps.mesh();

    polyMeshProps.allocate<MC_NODE_ID>(-1);
    polyMeshProps.allocate<MC_ARC_ID>(-1);
    polyMeshProps.allocate<MC_PATCH_ID>(-1);
    polyMeshProps.allocate<MC_BLOCK_ID>(-1);
    polyMeshProps.allocate<IS_SINGULAR>(false);
    polyMeshProps.allocate<IS_FEATURE_F>(0);
    polyMeshProps.allocate<IS_FEATURE_E>(0);
    polyMeshProps.allocate<IS_FEATURE_V>(0);
    polyMeshProps.template allocate<MARK_N>(0);
    polyMeshProps.template allocate<MARK_A>(0);
    polyMeshProps.template allocate<MARK_P>(0);
    polyMeshProps.template allocate<MARK_B>(0);

    map<FH, set<HFH>> p2hfsNew;
    for (FH p : mcMeshProps().mesh().faces())
    {
        set<VH> vs;
        set<EH> es;
        set<FH> fs;
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
            fs.insert(meshProps().mesh().face_handle(hf));
        for (FH f : fs)
            for (EH e : meshProps().mesh().face_edges(f))
                es.insert(e);
        for (EH e : es)
            for (VH v : meshProps().mesh().edge_vertices(e))
                vs.insert(v);
        for (VH v : vs)
            if (v2v.count(v) == 0)
            {
                VH vNew = polyMesh.add_vertex(meshProps().mesh().vertex(v));
                v2v[v] = vNew;
                VH n = meshProps().get<MC_NODE>(v);
                if (n.is_valid())
                {
                    polyMeshProps.set<MC_NODE_ID>(vNew, n.idx());
                    if (mcMeshProps().isAllocated<IS_FEATURE_V>() && mcMeshProps().get<IS_FEATURE_V>(n))
                        polyMeshProps.set<IS_FEATURE_V>(vNew, mcMeshProps().get<IS_FEATURE_V>(n));
                    if (mcMeshProps().isAllocated<MARK_N>())
                        polyMeshProps.set<MARK_N>(vNew, mcMeshProps().get<MARK_N>(n));
                }
            }
        for (EH e : es)
            if (e2e.count(e) == 0)
            {
                auto evs = meshProps().mesh().edge_vertices(e);
                EH eNew = polyMesh.add_edge(v2v.at(evs[0]), v2v.at(evs[1]));
                e2e[e] = eNew;
                EH a = meshProps().get<MC_ARC>(e);
                if (a.is_valid())
                {
                    polyMeshProps.set<MC_ARC_ID>(eNew, a.idx());
                    if (mcMeshProps().isAllocated<IS_FEATURE_E>() && mcMeshProps().get<IS_FEATURE_E>(a))
                        polyMeshProps.set<IS_FEATURE_E>(eNew, mcMeshProps().get<IS_FEATURE_E>(a));
                    if (mcMeshProps().isAllocated<IS_SINGULAR>() && mcMeshProps().get<IS_SINGULAR>(a))
                        polyMeshProps.set<IS_SINGULAR>(eNew, true);
                    if (mcMeshProps().isAllocated<MARK_A>())
                        polyMeshProps.set<MARK_A>(eNew, mcMeshProps().get<MARK_A>(a));
                }
            }
        for (FH f : fs)
            if (f2f.count(f) == 0)
            {
                auto fvs = meshProps().mesh().get_halfface_vertices(meshProps().mesh().halfface_handle(f, 0));
                FH fNew = polyMesh.add_face({v2v.at(fvs[0]), v2v.at(fvs[1]), v2v.at(fvs[2])});
                f2f[f] = fNew;
                {
                    polyMeshProps.set<MC_PATCH_ID>(fNew, p.idx());
                    if (mcMeshProps().isAllocated<IS_FEATURE_F>() && mcMeshProps().get<IS_FEATURE_F>(p))
                        polyMeshProps.set<IS_FEATURE_F>(fNew, mcMeshProps().get<IS_FEATURE_F>(p));
                    if (mcMeshProps().isAllocated<MARK_P>())
                        polyMeshProps.set<MARK_P>(fNew, mcMeshProps().get<MARK_P>(p));
                }
            }
        auto& phfsNew = p2hfsNew[p];
        for (HFH hf : mcMeshProps().ref<PATCH_MESH_HALFFACES>(p))
        {
            vector<VH> hfvs;
            for (VH v : meshProps().mesh().halfface_vertices(hf))
                hfvs.push_back(v2v.at(v));
            HFH hfNew = polyMesh.find_halfface(hfvs);
            phfsNew.insert(hfNew);
        }
    }
    for (CH b : mcMeshProps().mesh().cells())
    {
        vector<HFH> hfsNew;
        for (HFH hp : mcMeshProps().mesh().cell_halffaces(b))
        {
            for (HFH hfNew : p2hfsNew.at(mcMeshProps().mesh().face_handle(hp)))
            {
                if ((hp.idx() % 2) == 0)
                    hfsNew.push_back(hfNew);
                else
                    hfsNew.push_back(polyMesh.opposite_halfface_handle(hfNew));
            }
        }
        CH bNew = polyMesh.add_cell(hfsNew);
        polyMeshProps.set<MC_BLOCK_ID>(bNew, b.idx());
        if (mcMeshProps().isAllocated<MARK_B>())
            polyMeshProps.set<MARK_B>(bNew, mcMeshProps().get<MARK_B>(b));
    }

    return SUCCESS;
}

HexRemesher::RetCode HexRemesher::smoothSurface(HexMeshProps& hexMeshProps, int iter)
{
    HexMesh& hexMesh = hexMeshProps.mesh();

    vector<VH> vs;
    vs.reserve(hexMesh.n_logical_vertices());
    for (VH v : hexMesh.vertices())
        if (hexMesh.is_boundary(v))
            vs.push_back(v);
    for (VH v : hexMesh.vertices())
        if (!hexMesh.is_boundary(v))
            vs.push_back(v);
    for (int i = 0; i < iter; i++)
    {
        for (VH v : vs)
        {
            bool isBoundary = hexMesh.is_boundary(v);

            bool onFeature = hexMeshProps.get<IS_FEATURE_V>(v);
            if (!onFeature)
                for (FH f : hexMesh.vertex_faces(v))
                    if (hexMeshProps.get<IS_FEATURE_F>(f))
                        onFeature = true;
            if (!onFeature)
                for (EH e : hexMesh.vertex_edges(v))
                    if (hexMeshProps.get<IS_FEATURE_E>(e))
                        onFeature = true;
            if (!onFeature)
                for (EH e : hexMesh.vertex_edges(v))
                    if (hexMesh.is_boundary(e) && hexMeshProps.get<IS_SINGULAR>(e))
                        onFeature = true;
            if (onFeature)
                continue;

            {
                Vec3d avgPos(0, 0, 0);
                int numNeighbors = 0;
                for (HEH he : hexMesh.outgoing_halfedges(v))
                {
                    VH otherVert = hexMesh.to_vertex_handle(he);
                    // For each neighbour
                    if (!isBoundary || hexMesh.is_boundary(hexMesh.edge_handle(he)))
                    {
                        Vec3d pos = hexMesh.vertex(otherVert);
                        avgPos += pos;
                        numNeighbors++;
                    }
                }
                avgPos /= numNeighbors;
                Vec3d newPos = (avgPos + hexMesh.vertex(v)) / 2.0;

                if (isBoundary)
                {
                    Vec3d avgNormal(0, 0, 0);
                    for (HFH hf : hexMesh.vertex_halffaces(v))
                    {
                        if (hexMesh.is_boundary(hf))
                        {
                            std::vector<Vec3d> vertices;
                            for (VH vHF : hexMesh.halfface_vertices(hf))
                                vertices.emplace_back(hexMesh.vertex(vHF));
                            Vec3d normal = (vertices[1] - vertices[0]) % (vertices[3] - vertices[0]);
                            normal.normalize();
                            avgNormal += normal;
                        }
                    }
                    avgNormal.normalize();

                    bool tooSpread = false;
                    for (HFH hf : hexMesh.vertex_halffaces(v))
                    {
                        if (hexMesh.is_boundary(hf))
                        {
                            std::vector<Vec3d> vertices;
                            for (VH vHF : hexMesh.halfface_vertices(hf))
                                vertices.emplace_back(hexMesh.vertex(vHF));
                            Vec3d normal = (vertices[1] - vertices[0]) % (vertices[3] - vertices[0]);
                            normal.normalize();
                            if ((avgNormal | normal) < 0.8)
                            {
                                tooSpread = true;
                                break;
                            }
                        }
                    }
                    if (tooSpread)
                        continue;

                    Vec3d diff = newPos - hexMesh.vertex(v);
                    double dist = diff | avgNormal;
                    newPos = newPos - dist * avgNormal;
                }

                hexMesh.set_vertex(v, newPos);
            }
        }
    }
    return SUCCESS;
}

HexRemesher::RetCode HexRemesher::smoothSurface(PolyMeshProps& hexMeshProps, int iter)
{
    PolyMesh& hexMesh = hexMeshProps.mesh();

    vector<VH> vs;
    vs.reserve(hexMesh.n_logical_vertices());
    for (VH v : hexMesh.vertices())
        if (hexMesh.is_boundary(v))
            vs.push_back(v);
    for (VH v : hexMesh.vertices())
        if (!hexMesh.is_boundary(v))
            vs.push_back(v);
    for (int i = 0; i < iter; i++)
    {
        for (VH v : vs)
        {
            bool isBoundary = hexMesh.is_boundary(v);

            bool onFeature = hexMeshProps.get<IS_FEATURE_V>(v);
            if (!onFeature)
                for (FH f : hexMesh.vertex_faces(v))
                    if (hexMeshProps.get<IS_FEATURE_F>(f))
                        onFeature = true;
            if (!onFeature)
                for (EH e : hexMesh.vertex_edges(v))
                    if (hexMeshProps.get<IS_FEATURE_E>(e))
                        onFeature = true;
            if (!onFeature)
                for (EH e : hexMesh.vertex_edges(v))
                    if (hexMesh.is_boundary(e) && hexMeshProps.get<IS_SINGULAR>(e))
                        onFeature = true;
            if (onFeature)
                continue;

            {
                Vec3d avgPos(0, 0, 0);
                int numNeighbors = 0;
                for (HEH he : hexMesh.outgoing_halfedges(v))
                {
                    VH otherVert = hexMesh.to_vertex_handle(he);
                    // For each neighbour
                    if (!isBoundary || hexMesh.is_boundary(hexMesh.edge_handle(he)))
                    {
                        Vec3d pos = hexMesh.vertex(otherVert);
                        avgPos += pos;
                        numNeighbors++;
                    }
                }
                avgPos /= numNeighbors;
                Vec3d newPos = (avgPos + hexMesh.vertex(v)) / 2.0;

                if (isBoundary)
                {
                    Vec3d avgNormal(0, 0, 0);
                    for (HFH hf : hexMesh.vertex_halffaces(v))
                    {
                        if (hexMesh.is_boundary(hf))
                        {
                            std::vector<Vec3d> vertices;
                            for (VH vHF : hexMesh.halfface_vertices(hf))
                                vertices.emplace_back(hexMesh.vertex(vHF));
                            Vec3d normal = (vertices[1] - vertices[0]) % (vertices[3] - vertices[0]);
                            normal.normalize();
                            avgNormal += normal;
                        }
                    }
                    avgNormal.normalize();
                    bool tooSpread = false;
                    for (HFH hf : hexMesh.vertex_halffaces(v))
                    {
                        if (hexMesh.is_boundary(hf))
                        {
                            std::vector<Vec3d> vertices;
                            for (VH vHF : hexMesh.halfface_vertices(hf))
                                vertices.emplace_back(hexMesh.vertex(vHF));
                            Vec3d normal = (vertices[1] - vertices[0]) % (vertices[3] - vertices[0]);
                            normal.normalize();
                            if ((avgNormal | normal) < 0.8)
                            {
                                tooSpread = true;
                                break;
                            }
                        }
                    }
                    if (tooSpread)
                        continue;

                    Vec3d diff = newPos - hexMesh.vertex(v);
                    double dist = diff | avgNormal;
                    newPos = newPos - dist * avgNormal;
                }

                hexMesh.set_vertex(v, newPos);
            }
        }
    }
    return SUCCESS;
}

HexRemesher::RetCode HexRemesher::writeHexMesh(const HexMeshProps& hexMeshProps, const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out.good())
    {
        LOG(ERROR) << "Could not write to file " << filename;
        return FILE_INACCESSIBLE;
    }
    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    LOG(INFO) << "Writing extracted hex mesh to " << filename;

    OVM::IO::FileManager().writeStream(out, hexMeshProps.mesh());

    return SUCCESS;
}

HexRemesher::RetCode HexRemesher::writePolyHexMesh(const PolyMeshProps& hexMeshProps, const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out.good())
    {
        LOG(ERROR) << "Could not write to file " << filename;
        return FILE_INACCESSIBLE;
    }
    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    LOG(INFO) << "Writing extracted hex mesh to " << filename;

    OVM::IO::FileManager().writeStream(out, hexMeshProps.mesh());

    return SUCCESS;
}

} // namespace c4hex
