#ifndef C4HEX_MCBISECTOR_HPP
#define C4HEX_MCBISECTOR_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>
#include <MC3D/Mesh/MCMeshProps.hpp>
#include <MC3D/Mesh/TetMeshProps.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages the simplification of a given MC. All methods require a quantization to be present
 *
 */
class MCSplitter : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        PATH_ROUTE_FAILURE = 1,
    };

    /**
     * @brief Create an instance that handles cell splitting of the MC meta-mesh in \p meshProps .
     *
     * @param meshProps IN/OUT: mesh equipped with an MC
     */
    MCSplitter(TetMeshProps& meshProps);

    /**
     * @brief Bisect a patch \p p across direction \p dir relative to coord system of block incident on halfpatch 0.
     *
     * @param p IN: patch
     * @param dir IN: direction. must be one of the two arc directions of \p p . Bisector is parallel to dir.
     * @param allowLoop IN: whether to allow cutting an annular patch along its annular direction
     * @param psSub OUT: pairTT<FH> child patches if successful, else invalid handles
     * @return true if bisection successful, else false
     */
    bool bisectPatchAcrossDir(const FH& p, UVWDir dir, bool allowLoop, vector<FH>& psSub);

    /**
     * @brief Bisect \p b . Will first bisect any non-minimal patch on its boundary if present.
     *
     * @param b IN: block to refine
     * @param subBlocks IN: should be empty, OUT: all subblocks remaining after simplification
     * @return true if anything was actually simplified, false if block was already simple
     */
    bool bisectBlockOrPatch(const CH& b, vector<CH>& subBlocks);

    /**
     * @brief Cut block into two parts given an arc loop on its boundary. A bisecting patch will be inserted
     *        having the arc loop as its boundary and a minimal surface as its embedding.
     *
     * @param b IN: block
     * @param ringHas IN: arc loop on boundary of \p b
     * @param bsCut OUT: the two child blocks
     * @return FH the inserted bisecting patch
     */
    FH cutBlock(const CH& b, const vector<HEH>& ringHas, vector<CH>& bsCut);

  private:
    /**
     * @brief Bisect the interior of block \p b . Assumes boundary only contains simple patches and no t-junctions
     *        between arcs.
     *
     * @param b IN: block
     * @param subBlocks OUT: the two child blocks
     * @return true if bisection successful
     * @return false else
     */
    bool bisectBlockInterior(const CH& b, vector<CH>& subBlocks);

    /**
     * @brief Find a loop on the boundary of \p b that has no kink around its perimeter (i.e. lies in a plane under the
     *        quantization)
     *
     * @param b IN: block
     * @param haSeed IN: halfarc to start loop from
     * @return vector<HEH> loop of arcs, if successful, else empty
     */
    vector<HEH> findArcLoopOnBoundary(const CH& b, const HEH& haSeed);

    /**
     * @brief If one of the two (parallel) halfarcs is longer than the other under the quantization, cut it and replace
     *        the reference by the from-half
     *
     * @param ha1 IN/OUT: halfarc, updated if cut
     * @param ha2 IN/OUT: halfarc, updated if cut
     */
    void cutToMatch(HEH& ha1, HEH& ha2);

    /**
     * @brief Find a vertex in the embedding of \p ha approximately at a fraction \p frac of its length.
     *
     * @param ha IN: halfarc
     * @param frac IN: value between 0 and 1
     * @return VH background mesh vertex
     */
    VH findVertexAtLengthFraction(const EH& ha, double frac);

    /**
     * @brief Split halfarc \p ha so that its from-half has a quantized length of \p
     *
     * @param ha IN: halfarc
     * @param lengthToMatch IN: length to cut ha to
     * @return vector<HEH> the sub-halfarcs of ha (form-half first, to-half second)
     */
    vector<HEH> cutToMatchLength(const HEH& ha, int lengthToMatch);
};

} // namespace c4hex

#endif
