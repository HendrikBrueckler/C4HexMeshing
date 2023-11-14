#ifndef C4HEX_HEXEXTRACTORBASE_HPP
#define C4HEX_HEXEXTRACTORBASE_HPP

#include "C4Hex/Mesh/HexMeshProps.hpp"
#include "C4Hex/Mesh/PolyMeshProps.hpp"
#include <MC3D/Mesh/MCMeshNavigator.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages the extraction of a hex mesh from a tet mesh equipped with an MC and an IGM
 *
 * @tparam MESHPROPS the type of mesh to extract to. May be either PolyMeshProps or HexMeshProps
 */
template <typename MESHPROPS>
class HexExtractor : public MCMeshNavigator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
    };

    /**
     * @brief Create an instance that manages extraction of a hex mesh from a tet mesh equipped with MC and IGM
     *
     * @param tetMeshProps IN: tet mesh equipped with MC and IGM
     * @param hexMeshProps OUT: hex mesh will be extracted here
     */
    HexExtractor(const TetMeshProps& tetMeshProps, MESHPROPS& hexMeshProps);

    /**
     * @brief Extract the hex mesh.
     *
     * @param how IN: many times to divide the integer grid facets (1x1) along the two facet axes.
     *                only the integer grid facets is refined, not the inside of the cells!
     *                This only has an effect if MESHPROPS == PolyMeshProps!
     * @return RetCode SUCCESS
     */
    RetCode extractHexMesh(int subdiv);

  protected:
    /**
     * @brief Add a hex vertex for every MC node
     *
     * @return map<VH, VH> mapping of nodes to hex vertices
     */
    map<VH, VH> createNodeHexV();

    /**
     * @brief Add hex vertices for integer grid points on arcs (incorporates node mappings) and connect by edges
     *
     * @param n2hexV IN: mapping of nodes to hex vertices
     * @param subdiv IN: how many times to divide integer grid segments (length 1) along each arc
     * @return map<EH, vector<VH>> mapping of arcs to ordered hex vertices on that arc
     */
    map<EH, vector<VH>> createArcHexVE(const map<VH, VH>& n2hexV, int subdiv);

    /**
     * @brief Add hex vertices for integer grid points on patches (incorporates arc mappings) and connect by edges and
     *        faces
     *
     * @param a2delta2hexV IN: mapping of arcs to ordered hex vertices on that arc
     * @param subdiv IN: how many times to divide the integer grid facets (1x1) along the two facet axes
     * @return map<FH, map<Vec3Q, VH>> mapping of patches to hex vertices by IGM (in coord
     *                                                             system of first halfpatches block)
     */
    map<FH, map<Vec3Q, VH>> createPatchHexVEF(const map<EH, vector<VH>>& a2delta2hexV, int subdiv);

    /**
     * @brief Add hex vertices for integer grid points in blocks (incorporates patch mappings) and connect by edges,
     *        faces and cells
     *
     * @param p2igm2hexV mapping of patches to hex vertices by IGM (in coord system of first halfpatches block)
     * @param subdiv IN: how many times to divide the integer grid facets (1x1) along the two facet axes.
     *                   only the integer grid facets is refined, not the inside of the cells!
     * @return map<CH, map<Vec3Q, VH>> mapping of blocks to hex vertices by IGM
     */
    map<CH, map<Vec3Q, VH>> createBlockHexVEFC(const map<FH, map<Vec3Q, VH>>& p2igm2hexV, int subdiv);

    /**
     * @brief Determine the barycentric coordinates of \p igmUVW wrt \p hf given that \p coord1 and \p coord3
     *        are the variable coordinates
     *
     * @param hf IN: triangle wrt which the barycentric coordinates are to be determined
     * @param igmUVW IN: igm of the point for which to determine barycentric coordinates
     * @param coord1 IN: first variable coordinate
     * @param coord3 IN: second variable coordinate
     * @param barCoords OUT: barycentric coordinates of \p igmUVW wrt \p hf (length of 3 if inside, else empty)
     * @return true if inside
     * @return false else
     */
    bool barycentricCoordsIGM(const HFH& hf, const Vec3Q& igmUVW, int coord1, int coord3, Vec3Q& barCoords) const;

    /**
     * @brief Determine the barycentric coordinates of \p igmUVW wrt \p tet
     *
     * @param tet IN: tet wrt which the barycentric coordinates are to be determined
     * @param igmUVW IN: igm of the point for which to determine barycentric coordinates
     * @param barCoords OUT: barycentric coordinates of \p igmUVW wrt \p hf (length of 4 if inside, else empty)
     * @return true if inside
     * @return false else
     */
    bool barycentricCoordsIGM(const CH& tet, const Vec3Q& igmUVW, Vec4Q& barCoords) const;

    MESHPROPS& _hexMeshProps;
};

} // namespace c4hex

#endif
