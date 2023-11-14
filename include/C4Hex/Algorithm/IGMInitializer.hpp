#ifndef C4HEX_PARAMRESCALER_HPP
#define C4HEX_PARAMRESCALER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages initialization of the integer grid map from a given quantization.
 *        Not guaranteed inversion-free.
 *
 */
class IGMInitializer : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        INVALID_ELEMENTS = 26,
    };

    /**
     * @brief Create an instance initializing the IGM on given mesh. Requires a quantization.
     *
     * @param meshProps IN/OUT: mesh for which an IGM parametrization should be generated
     */
    IGMInitializer(TetMeshProps& meshProps);

    /**
     * @brief Actually rescale the seamless UVW param into an IGM parametrization.
     *
     * @return RetCode SUCCESS
     */
    RetCode initializeFromQuantization();

    /**
     * @brief Initialize IGM for all node-vertices
     */
    void initializeOnNodes();

    /**
     * @brief Initialize IGM for all vertices on arcs (without nodes)
     */
    void initializeOnArcs();

    /**
     * @brief Initialize IGM for all vertices on patches (without arcs)
     */
    void initializeOnPatches();

    /**
     * @brief Initialize IGM for all vertices within blocks patches (without arcs)
     */
    void initializeInBlocks();

    /**
     * @brief Determine the transitions between blocks as implied by the MC quantization
     *        and the initial rescaled node IGMs
     */
    void determineIGMTransitions();

    /**
     * @brief Split block-inner edges of tets that are fully constrained and inverted
     *
     * @return true if any edges were split
     * @return false else
     */
    bool splitSomeNecessaryForInjectiveInterior();

    /**
     * @brief Split block-inner edges of tets that are fully constrained
     *
     * @return true if any edges were split
     * @return false else
     */
    bool splitAllSufficientForInjectiveInterior();

    /**
     * @brief Split all edges/faces on block boundary that prevent a bijective boundary map
     *
     * @return true if any edges or faces were split
     * @return false else
     */
    bool splitAllNecessaryForInjectiveBoundary();

    /**
     * @brief Walking along halfedges of \p haTrace trace the volumetric path to the to-node.
     *        Necessary to handle selfadjacency.
     *
     * @param startTet IN: tet in volumetric sector on from-node of \p haTrace
     * @param haTrace IN: halfarc to trace along
     * @param tetsVisited OUT: tets visited along the volumetric path trace
     * @return CH tet on to-node of \p haTrace in same volumetric sector as \p startTet
     */
    CH traceToEndWithinBlock(const CH& startTet, const HEH& haTrace, set<CH>& tetsVisited) const;

    /**
     * @brief Minimize the cut surface of the IGM. Keeps coordinate system orientations consistent with seamless param.
     */
    void minimizeCutSurface();

    /**
     * @brief Allocate needed mesh properties to store IGM and transitions
     */
    void allocateIGM();

  private:
    /**
     * @brief Use 2D-Tutte to rescale the vertex IGMs for patches
     *
     * @param meanValWeights IN: whether to use mean value weights instead of uniform weights
     * @param xyzTarget IN: whether to calculate weights based on XYZ instead of seamless param UVW
     */
    void initializeOnPatches2DTutte(bool meanValWeights, bool xyzTarget);

    /**
     * @brief Calculate the mean value weights for block \p b for all edges in \p edges
     *
     * @param b IN: block
     * @param edges IN: patch edges of \p b for which mean 2D value weights should be computed
     * @param xyzTarget IN: whether to calculate weights based on XYZ instead of seamless param UVW
     * @return map<EH, double>
     */
    map<EH, double> calc2DmeanValueWeights(const CH& b, const set<EH>& edges, bool xyzTarget) const;

    /**
     * @brief Use 3D-Tutte to rescale the vertex IGMs in block interiors
     *
     * @param cotWeights IN: whether to use cotangent weights instead of uniform weights
     * @param xyzTarget IN: whether to calculate weights based on XYZ instead of seamless param UVW
     */
    void initializeInBlocks3DTutte(bool cotWeights, bool xyzTarget);

    /**
     * @brief Calculate the cotangent weights for all mesh edges
     *
     * @param xyzTarget IN: whether to calculate weights based on XYZ instead of seamless param UVW
     * @return vector<double> cotangent weight per mesh edge (in order of idx())
     */
    vector<double> calc3DcotangentWeights(bool xyzTarget) const;

    /**
     * @brief Computes 2D tutte embedding for patch-interior vertices given fixed patch boundary.
     *
     * @tparam Scalar double or Q
     * @param vtx2index IN: entry for each patch-interior vertex and its index (in range 0 to vtx2index.size() - 1)
     * @param e2weight IN: edge weights
     * @param isoCoord IN: which of the 3D coordinates is const
     * @param front IN: wether to compute for the front-facing halfface of patch
     * @param isoValue IN: value of the const coordinate
     * @param hfs IN: halffaces of patch
     */
    template <typename Scalar>
    void solve2DTutte(const map<VH, int>& vtx2index,
                      const map<EH, double>& e2weight,
                      int isoCoord,
                      Q isoValue,
                      bool front,
                      const set<HFH>& hfs);
};

} // namespace c4hex

#endif
