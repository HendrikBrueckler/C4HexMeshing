#ifndef C4HEX_MCCOLLAPSER_HPP
#define C4HEX_MCCOLLAPSER_HPP

#include "C4Hex/Algorithm/MCSplitter.hpp"

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages the simplification of a given MC
 */
class MCCollapser : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        INVALID_QUANTIZATION = 27,
        REROUTE_ERROR = 28,
        COLLAPSE_ERROR = 29
    };

    /**
     * @brief Create an instance that collapses the MC meta-mesh and its embedding for a \p meshProps .
     *
     * @param meshProps IN/OUT: mesh equipped with an MC
     */
    MCCollapser(TetMeshProps& meshProps);

    /**
     * @brief Whether the MC has 0-arcs
     */
    bool hasZeroLengthArcs() const;

    /**
     * @brief Collapse all zero elements in the MC and update its embedding accordingly
     *
     * @param optimize IN: whether to smooth the resulting MC afterwards
     * @param randomOrder IN: which order of collapses (random or smalles elements first)
     * @param direction IN: direction of collapses (0: globally coordinated, 1: random, 2: from min valence)
     * @return RetCode SUCCESS or error code
     */
    RetCode collapseAllZeroElements(bool optimize = false, bool randomOrder = false, int direction = 0);

    /**
     * @brief Count the number of zero elements according to the underlying quantization
     *
     * @param numZeroAs OUT: number of zero-arcs
     * @param numZeroPs OUT: number of zero-patches
     * @param numZeroBs OUT: number of zero-blocks
     */
    void countZeroElements(int& numZeroAs, int& numZeroPs, int& numZeroBs);

    /**
     * @brief Whether any uncollapsed zero-elements remain (fully or partially)
     *
     * @return true if some remain
     * @return false else
     */
    bool zeroElementsRemain();

    /**
     * @brief Mark all zero elements (and subelements of those) using properties MARK_N, MARK_A, MARK_P, MARK_B
     *
     */
    void markZeros();

  private:
    /**
     * @brief Assign globally preferred collapse directions to each block
     */
    void assignCollapseDirs();

    /**
     * @brief Whether the collapse of from[ha] onto to[ha] is locked.
     *        A collapse is locked if from[ha] is a singular node,
     *        or if from[ha] is a node on a singular arc but ha is not a singular arc
     *        or if from[ha] is a boundary node but ha is not a boundary arc.
     *
     *
     * @param ha IN: halfarc to collapse
     * @return true if locked
     * @return false else
     */
    bool collapseIsLocked(const HEH& ha) const;

    /**
     * @brief Get the preferred collapse direction of \p a
     *
     * @param a IN: zero-arc
     * @return HEH directed halfarc
     */
    HEH preferredCollapseHalfarc(const EH& a) const;

    /**
     * @brief Determine the stationary side of a pillow patch collapse
     *
     * @param haMoving IN: any halfedge of collapse patch, OUT: moving side of patch
     * @param haStationary IN: other halfedge of collapse patch, OUT: stationary side of patch
     */
    void determineStationaryEnd(HEH& haMoving, HEH& haStationary) const;

    /**
     * @brief Determine the stationary side of a pillow block collapse
     *
     * @param hpMoving IN: any halfpatch of collapse block, OUT: moving side of block
     * @param hpStationary IN: other halfpatch of collapse block, OUT: stationary side of block
     */
    void determineStationaryEnd(HFH& hpMoving, HFH& hpStationary) const;

    /**
     * @brief Directed arc collapse
     *
     * @param haCollapse IN: zero-halfarc
     * @return RetCode SUCCESS or error code
     */
    RetCode collapseHalfarc(const HEH& haCollapse);

    /**
     * @brief Directed pillow patch collapse. Ideal direction is chosen automatically.
     *
     * @param pCollapse IN: patch to collapse
     * @param has IN: halfarcs of frontside of \p pCollapse
     * @return RetCode SUCCESS or error code
     */
    RetCode collapsePillowPatch(const FH& pCollapse, const set<HEH>& has);

    /**
     * @brief Directed pillow block collapse. Ideal direction is chosen automatically.
     *
     * @param bCollapse IN: block to collapse
     * @return RetCode SUCCESS or error code
     */
    RetCode collapsePillowBlock(const CH& bCollapse);

    /**
     * @brief Cigar block collapse.
     *
     * @param bCollapse IN: block to collapse
     * @return RetCode SUCCESS or error code
     */
    RetCode collapseCigarBlock(const CH& bCollapse);

    /**
     * @brief Handle the MC connectivity update of the arc collapse operation
     *
     * @param haCollapse IN: collapse halfarc
     */
    void collapseArcConnectivity(const HEH& haCollapse);

    /**
     * @brief Handle the MC connectivity update of the patch collapse operation
     *
     * @param p IN: collapse patch
     * @param haRemoved IN: halfarc of \p p to be removed
     * @param haRemaining IN: halfarc of \p p to remain
     */
    void collapsePillowPatchConnectivity(const FH& p, const HEH& haRemoved, const HEH& haRemaining);

    /**
     * @brief Handle the MC connectivity update of the block collapse operation
     *
     * @param b IN: collapse block
     * @param hpRemoved IN: halfpatch of \p b to be removed
     * @param hpRemaining IN: halfpatch of \p b to remain
     */
    void collapsePillowBlockConnectivity(const CH& b, const HFH& hpRemoved, const HFH& hpRemaining);

    /**
     * @brief Collapse the next encountered zero arc
     *
     * @return true if one was found and collapsed
     * @return false else
     */
    bool collapseNextZeroArc();

    /**
     * @brief Collapse the next encountered pillow patch
     *
     * @return true if one was found and collapsed
     * @return false else
     */
    bool collapseNextPillowPatch();

    /**
     * @brief Collapse the next encountered pillow block
     *
     * @return true if one was found and collapsed
     * @return false else
     */
    bool collapseNextPillowBlock();

    /**
     * @brief Bisect the next encountered pillow block, that has a pillow patch bisector.
     *        This operation is needed to avoid some extreme selfadjacency which OVM data
     *        structure can not handle.
     *
     * @return true if one was found and bisected
     * @return false else
     */
    bool bisectNextBlockByPillowPatch();

    /**
     * @brief Bisect the next encountered zero patch.
     *
     * @param allowZeroLoop IN: whether to allow to bisect an arc-selfincident pillow patch (creating a temporary 0-loop
     * @return true if one was found and bisected
     * @return false else
     */
    bool bisectNextZeroPatch(bool allowZeroLoop);

    /**
     * @brief Bisect the next encountered almost-pillow patch.
     *
     * @param allowZeroLoop IN: whether to allow to bisect an arc-selfincident pillow patch (creating a temporary
     * 0-loop)
     * @return true if one was found and bisected
     * @return false else
     */
    bool bisectNextAlmostPillowPatch(bool allowZeroLoop);

    /**
     * @brief  Bisect the next encountered almost-pillow block.
     *
     * @return true if one was found and bisected
     * @return false else
     */
    bool bisectNextAlmostPillowBlock();

    /**
     * @brief Get all almost-pillow patches
     *
     * @return set<FH> almost-pillow patches
     */
    set<FH> findAlmostPillowPatches() const;

    /**
     * @brief Get all almost-pillow blocks
     *
     * @return set<FH> almost-pillow blocks
     */
    set<CH> findAlmostPillowBlocks() const;

    int _nCollapsedAs = 0; // Accumulated number of arc collapses
    int _nCollapsedPs = 0; // Accumulated number of patch collapses
    int _nCollapsedBs = 0; // Accumulated number of block collapses

    int _nTetsPre = 0; // Number of tets before collapsing
    int _nFsPre = 0;   // Number of faces before collapsing
    int _nEsPre = 0;   // Number of edges before collapsing
    int _nVsPre = 0;   // Number of vertices before collapsing

    int _nBsPre = 0; // Number of blocks before collapsing
    int _nPsPre = 0; // Number of patches before collapsing
    int _nAsPre = 0; // Number of arcs before collapsing
    int _nNsPre = 0; // Number of nodes before collapsing

    int _numZeroBs = 0; // Number of zero-arcs before collapsing
    int _numZeroPs = 0; // Number of zero-patches before collapsing
    int _numZeroAs = 0; // Number of zero-blocks before collapsing

    MCSplitter _refiner; // Internal refiner for bisection operations

    bool _randomOrder = false; // Whether to execute collapses in random order
    int _direction = 0;        // Directedness of collapse process
};

} // namespace c4hex

#endif
