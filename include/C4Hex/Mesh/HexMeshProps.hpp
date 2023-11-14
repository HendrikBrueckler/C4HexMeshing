#ifndef MC3D_HEXMESHPROPS_HPP
#define MC3D_HEXMESHPROPS_HPP

#include <MC3D/Data/BlockData.hpp>
#include <MC3D/Mesh/MCMeshProps.hpp>
#include <MC3D/Mesh/TetMeshProps.hpp>
#include <MC3D/Mesh/MeshPropsInterface.hpp>
#include <MC3D/Types.hpp>

namespace mc3d
{

using HexMeshPropsBase = MeshPropsInterface<HexMesh,
                                            MC_BLOCK_ID,
                                            MC_PATCH_ID,
                                            MC_ARC_ID,
                                            MC_NODE_ID,
                                            IS_SINGULAR,
                                            IS_FEATURE_F,
                                            IS_FEATURE_E,
                                            IS_FEATURE_V,
                                            MARK_N,
                                            MARK_P,
                                            MARK_A,
                                            MARK_B>;

/**
 * @brief Class/struct to manage predefined properties of a raw hex mesh
 */
class HexMeshProps : public HexMeshPropsBase
{
  public:
    /**
     * @brief Create a property wrapper around \p mesh_
     *
     * @param mesh_ IN/OUT: raw mesh to augment by a set of predefined properties
     */
    HexMeshProps(HexMesh& mesh_) : HexMeshPropsBase(mesh_)
    {
    }

};

} // namespace mc3d

#endif
