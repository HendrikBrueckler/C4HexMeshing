#ifndef C4HEX_IGMGENERATOR_HPP
#define C4HEX_IGMGENERATOR_HPP

#include <MC3D/Mesh/MCMeshNavigator.hpp>
#include <MC3D/Mesh/TetMeshManipulator.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages the generation of an integer grid map (IGM)
 *        from a seamless parametrization on a tet mesh equipped with an MC.
 *
 */
class IGMGenerator : public virtual TetMeshManipulator, public virtual MCMeshNavigator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        QUANTIZATION_ERROR = 20,
        INITIALIZATION_ERROR = 21,
        UNTANGLING_ERROR = 22,
        NO_CONVERGENCE = 23,
    };

    /**
     * @brief Create an instance that manages the generation of an IGM from a seamless param and an MC
     *
     * @param meshProps IN/OUT: tet mesh equipped with a seamless param and an MC
     */
    IGMGenerator(TetMeshProps& meshProps);

    /**
     * @brief From a quantized MC with nonzero arc lengths, generate IGM for each block.
     *
     * @param simplifyBaseMesh IN: whether to decimate and remesh base mesh to improve condition of param problem
     * @param maxUntanglingIter IN: maximum iterations of outer untangling iterations to eliminate parametric inversions
     * @param maxInnerIter IN: maximum iterations of inner untangling iterations to eliminate parametric inversions
     * @return RetCode SUCCESS, QUANTIZATION_ERROR or RESCALING_ERROR
     */
    RetCode generateBlockwiseIGM(bool simplifyBaseMesh, int maxInnerIter = 500, int maxUntanglingIter = 40);
};

} // namespace c4hex

#endif
