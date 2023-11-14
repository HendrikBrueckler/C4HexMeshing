#ifndef C4HEX_HEXREMESHER_HPP
#define C4HEX_HEXREMESHER_HPP

#include "C4Hex/Mesh/HexMeshProps.hpp"
#include "C4Hex/Mesh/PolyMeshProps.hpp"
#include <MC3D/Mesh/MCMeshNavigator.hpp>
#include <MC3D/Mesh/TetMeshProps.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages the extraction of a hex mesh from a tet mesh equipped with an MC and an IGM
 */
class HexRemesher : public virtual MCMeshNavigator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        FILE_INACCESSIBLE = 1,
    };

    /**
     * @brief Create an instance that manages extraction of a hex mesh from a tet mesh equipped with MC and IGM
     *
     * @param hexMeshProps OUT: hex mesh will be extracted here
     * @param meshProps IN: tet mesh equipped with MC and IGM
     */
    HexRemesher(const TetMeshProps& meshProps);

    /**
     * @brief Extract the hex mesh
     *
     * @param hexMeshProps OUT: hex mesh will be extracted here
     * @return RetCode SUCCESS
     */
    RetCode extractHexMesh(HexMeshProps& hexMeshProps);

    /**
     * @brief Extract the (poly)-hex mesh
     *
     * @param hexMeshProps OUT: hex mesh will be extracted here
     * @return RetCode SUCCESS
     */
    RetCode extractPolyHexMesh(PolyMeshProps& hexMeshProps, int subDiv);

    /**
     * @brief Extract the MC mesh
     *
     * @param hexMeshProps OUT: MC mesh will be extracted here
     * @return RetCode SUCCESS
     */
    RetCode extractMCMesh(PolyMeshProps& hexMeshProps);

    /**
     * @brief Write the extracted hex mesh to \p filename
     *
     * @param filename IN: name of the file to write to
     * @return RetCode SUCCESS or FILE_INACCESSIBLE
     */
    RetCode writeHexMesh(const HexMeshProps& hexMeshProps, const std::string& filename) const;

    /**
     * @brief Write the extracted hex mesh to \p filename
     *
     * @param filename IN: name of the file to write to
     * @return RetCode SUCCESS or FILE_INACCESSIBLE
     */
    RetCode writePolyHexMesh(const PolyMeshProps& hexMeshProps, const std::string& filename) const;

    /**
     * @brief Smooth the surface of the hex mesh by shifting vertices towards the center of their neighbors.
     *        Respects curvature by projecting into average plane before shift.
     *
     * @param iter how often each vertex should be shifted
     * @return RetCode SUCCESS
     */
    RetCode smoothSurface(HexMeshProps& hexMeshProps, int iter);

    /**
     * @brief Smooth the surface of the hex mesh by shifting vertices towards the center of their neighbors.
     *        Respects curvature by projecting into average plane before shift.
     *
     * @param iter how often each vertex should be shifted
     * @return RetCode SUCCESS
     */
    RetCode smoothSurface(PolyMeshProps& hexMeshProps, int iter);
};

} // namespace c4hex

#endif
