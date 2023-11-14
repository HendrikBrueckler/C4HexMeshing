#ifndef C4HEX_MCSMOOTHER_HPP
#define C4HEX_MCSMOOTHER_HPP

#include "C4Hex/Algorithm/MCCollapser.hpp"

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages smoothing of jaggy MC embedding
 */
class MCSmoother : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        REROUTE_ERROR = 28,
        COLLAPSE_ERROR = 29
    };

    /**
     * @brief Create an instance that smoothes the MC meta-meshs embedding for a \p meshProps .
     *
     * @param meshProps IN/OUT: mesh equipped with an MC
     */
    MCSmoother(TetMeshProps& meshProps);

    /**
     * @brief Execute the smoothing (via shortest path, minimal surface in UVW).
     *        Success is not guaranteed but failsafes are in place and a best-effort
     *        smoothing is always executed.
     */
    void smoothMC();

  private:
    /**
     * @brief Tries to push transitions to the boundary of \p space .
     *
     * @param space IN: volume-connected portion of object
     * @param tetSeed IN: tet whose coordinate system should be propagated into entire \p space
     * @return RetCode SUCCESS or error code
     */
    RetCode makeVolumeTransitionFree(const set<CH>& space, const CH& tetSeed);

    /**
     * @brief Initialize internal variables, describing the available space for rerouting around arc \p a
     *
     * @param a IN: arc around which rerouting will take place
     */
    void collectAllowedRegionAroundArc(const EH& a);

    /**
     * @brief Cyclically order the patches around \p aReroute starting from \p hp
     *
     * @param hp IN: halfpatch to start cycle from (in both cw and ccw)
     * @param aReroute IN: arc
     * @return list<FH> List of halfpatches ordered by circular cell distance from \p hp
     */
    list<FH> cyclicOrderPatches(const HFH& hp, const EH& aReroute) const;

    /**
     * @brief Determine block embedding by floodfilling blocks.
     *        Requires consistent patch embeddings for all involved blocks, else fails.
     *
     * @return RetCode SUCCESS or error code
     */
    RetCode refloodFillBlocks();

    /**
     * @brief Reroute arc \p a by computing a shortest path between its endpoints.
     *        Success is not guaranteed for tricky MC connectivity.
     *
     * @param a IN: arc
     * @param change OUT: whether the arcs embedding changed at all
     * @return RetCode SUCCESS or error code
     */
    RetCode reroute(const EH& a, bool* change = nullptr);

    /**
     * @brief Reroute patch \p p by computing a minimal surface within its boundary.
     *        Success is not guaranteed for tricky MC connectivity.
     *
     * @param p IN: patch
     * @return RetCode SUCCESS or error code
     */
    RetCode reroute(const FH& p);

    /**
     * @brief Whether arc \p a is already perfectly U-, V- or W-aligned
     *
     * @param a IN: arc
     * @return true if completely aligned
     * @return false else
     */
    bool isUVWAligned(const EH& a) const;

    /**
     * @brief Whether patch \p a is already perfectly U-, V- or W-aligned
     *
     * @param p IN: patch
     * @return true if completely aligned
     * @return false else
     */
    bool isUVWAligned(const FH& p) const;

    /**
     * @brief Update the embedding of arc \p a to lie on halfedges \p hesA
     *
     * @param a IN: arc
     * @param hesA IN: sequence of halfedges, OUT: after replacing by children in case of refinement
     */
    void updateArcEmbedding(const EH& a, list<HEH>& hesA);

    /**
     * @brief Update the embedding of patch \p p to lie on halffaces \p hfsP
     *
     * @param p IN: patch
     * @param hfsP IN: surface of halffaces, OUT: after replacing by children in case of refinement
     * @param transP IN: optional transition to set on the newly embedded patch
     */
    void updatePatchEmbedding(const FH& p, set<HFH>& hfsP, Transition* transP = nullptr);

    set<CH> _allowedVolume;              // The tetrahedra available for rerouting the currently processed MC elements
    set<FH> _forbiddenFs;                // Faces forbidden for rerouting (no other element may be embedded here)
    set<EH> _forbiddenEs;                // Edges forbidden for rerouting (no other element may be embedded here)
    set<VH> _forbiddenVs;                // Vertices forbidden for rerouting (no other element may be embedded here)
    map<CH, vector<FH>> _b2unaffectedPs; // Patches unaffected by reroute for each involved block
    map<CH, vector<EH>> _b2unaffectedAs; // Arcs unaffected by reroute for each involved block
    map<CH, vector<VH>> _b2unaffectedNs; // Nodes unaffected by reroute for each involved block

    set<FH> _patchesRerouted;        // Patches already rerouted (in the process of an arc reroute)
    set<CH> _blocksChanged;          // Blocks affected in the process of an ongoing reroute
    EH _aReroute;                    // Arc currently being rerouted
    map<EH, set<HFH>> _traversedHfs; // Halffaces traversed by boundary arcs during their reroute

    set<FH> _torusSplitter;  // Marked to prevent reroute through these faces (for toroidal allowedspace)
    set<EH> _torusSplitterE; // Marked to prevent reroute through these edges (for toroidal allowedspace)
    set<VH> _torusSplitterV; // Marked to prevent reroute through these vertices (for toroidal allowedspace)

    map<EH, set<FH>> _aBoundary2sectorFront; // To force arc reroutes to pass through a surface sector on frontside
    map<EH, set<CH>> _a2sectorFront;         // To force arc reroutes to pass through a volume sector on frontside
    map<EH, set<FH>> _aBoundary2sectorBack;  // To force arc reroutes to pass through a surface sector on backside
    map<EH, set<CH>> _a2sectorBack;          // To force arc reroutes to pass through a volume sector on backside

    map<FH, set<CH>> _p2sector;    // To force patch reroutes to pass through a volume sector
    map<FH, set<HEH>> _p2boundary; // Store the original boundary cycle of a patch
};

} // namespace c4hex

#endif
