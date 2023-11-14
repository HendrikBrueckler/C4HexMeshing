#ifndef C4HEX_EMBEDDINGCOLLAPSER_HPP
#define C4HEX_EMBEDDINGCOLLAPSER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages collapsing of an MCs embedding
 */
class EmbeddingCollapser : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        REROUTE_ERROR = 28
    };

    /**
     * @brief Create an instance managing collapsing of an MCs embedding in \p meshProps
     *
     * @param meshProps IN: mesh for which collapse operations should be performed
     */
    EmbeddingCollapser(TetMeshProps& meshProps);

    /**
     * @brief One-directional incremental collapse of the embedding of halfarc \p ha
     *
     * @param ha IN: halfarc whose embedding should be collapsed (i.e. vanish)
     * @return RetCode SUCCESS or error code
     */
    RetCode collapseArcEmbedding(const HEH& ha);

    /**
     * @brief One-directional incremental collapse of the embedding of patch \p p
     *
     * @param p IN: patch to collapse
     * @param haMoving IN: halfarc bounding \p p which is shifted towards \p haStationary
     * @param haStationary IN: halfarc bounding \p p which is stationary
     * @return RetCode SUCCESS or error code
     */
    RetCode collapsePillowPatchEmbedding(const FH& p, const HEH& haMoving, const HEH& haStationary);

    /**
     * @brief One-directional all-at-once collapse of the embedding of block \p p
     *
     * @param b IN: block to collapse
     * @param hpMoving IN: halfpatch bounding \p p which is shifted towards \p hpStationary
     * @param hpStationary IN: halfpatch bounding \p p which is stationary
     * @return RetCode SUCCESS or error code
     */
    RetCode collapsePillowBlockEmbedding(const CH& b, const HFH& hpMoving, const HFH& hpStationary);

    /**
     * @brief All-at-once collapse of the embedding of block \p b
     *
     * @param b IN: cigar block to collapse
     * @return RetCode SUCCESS or error code
     */
    RetCode collapseCigarBlockEmbedding(const CH& b);

  private:
    /**
     * @brief Represents a dome over a vertex. This is the maximum volume-connected set of tetrahedra
     *        bounded by block walls which are all incident to one vertex on the block's boundary
     *
     */
    struct DomeElements
    {
        set<HFH> floor;                 // all halffaces on the boundary of the dome that touch the center vertex
        set<HFH> ceiling;               // all halffaces on the boundary of the dome that do NOT touch the center vertex
        set<EH> floorInnerEdges;        // all inner edges of floor (without those shared with ceiling)
        set<EH> ceilingInnerEdges;      // all inner edges of ceiling (without those shared with floor)
        set<VH> ceilingInnerVertices;   // all inner vertices of ceiling (without those shared with floor)
        set<HEH> ceilingBorder;         // all border edges of ceiling (those shared with ceiling)
        set<VH> ceilingBorderVertices;  // all border vertices of ceiling (those shared with ceiling)
        list<HEH> ceilingBorderOrdered; // all border edges of ceiling in cyclic order (empty if no single cycle exists)
    };

    /**
     * @brief Shift node sitting on from[he] across \p he
     *
     * @param he IN: halfedge with a node associated to from[he]
     */
    void nodeShift(const HEH& he);

    /**
     * @brief Get the fan of triangles starting from \p heStart within the same patch as \p patchHf and get the other
     *        limit edge. All triangles in the fan will be incident on from[heStart].
     *
     * @param heStart IN: halfedge from which fan grows
     * @param patchHf IN: patch within which to grow fan
     * @return pair<HEH, vector<HFH>> first: next arc-edge on other side of the fan, second: halffaces in fan
     */
    pair<HEH, vector<HFH>> patchCorner(const HEH& heStart, const HFH& patchHf) const;

    /**
     * @brief Whether this instance is currently working on a patch collapse
     *
     * @return true if this instance collapses a patch currently
     * @return false if it instead collapses an arc
     */
    bool currentlyCollapsingPatch() const;

    /**
     * @brief Reset state variables
     */
    void resetVars();

    /**
     * @brief Set the internal state variables for a collapse of \p ha
     *
     * @param ha IN: halfarc to be collapsed
     */
    void setVars(const HEH& ha);

    /**
     * @brief Set the internal state variables for a collapse of \p p
     *
     * @param p IN: patch to be collapsed
     * @param haMoving IN: non-stationary patch side
     * @param haStationary IN: stationary patch side
     */
    void setVars(const FH& p, const HEH& haMoving, const HEH& haStationary);

    /**
     * @brief Reconnect next left-behind arc to shifted node by shifting it through a patch surface
     *
     * @return true if any arc was reconnected
     * @return false else
     */
    bool reconnectNextArcToNodeThroughPatch();

    /**
     * @brief Reconnect next left-behind arc to shifted node by shifting it through a block volume
     *
     * @return true if any arc was reconnected
     * @return false else
     */
    bool reconnectNextArcToNodeThroughBlock();

    /**
     * @brief Pull patch away from a vertex where multiple patches are entangled
     *
     * @return true if any arc was reconnected
     * @return false else
     */
    bool pullAwayNextPatch();

    /**
     * @brief Shift an arc across a fan of triangles
     *
     * @param he IN: Halfedge on one side of the fan. The arc sitting on this halfedge will be shifted
     * @param triangleFan IN: halffaces in the triangle fan, OUT: after refinement and update
     */
    void arcShift(const HEH& he, vector<HFH>& triangleFan);

    /**
     * @brief Order halffaces around he.
     *
     * Fat drawn faces are block walls
     *
     *      1     hf     7
     *       \    hf    /
     *        \   hf   /
     *         \  hf  /
     *          \ hf /
     *           \hf/
     * 2 =========he========= 6
     *           /||\
     *          / || \
     *         /  ||  \
     *        /   ||   \
     *       /    ||    \
     *      3      4     5
     *
     * Order is:
     * Front: [2,1], [4,3], Back: [6, 7], [4, 5]
     *
     * Inserted halffaces point back towards hf
     * Last halfface in return.first is opposite of last halfface in return.second
     *
     * @param he IN: halfedge
     * @param hf IN: halfface to start peeling from (in both directions)
     * @param fJoin IN: optional prescribed face separating cw and ccw expansion domains (e.g. 4 in above example)
     * @return pairTT<list<list<HFH>>> ring of halffaces around \p he ordered as shown above
     */
    pairTT<list<list<HFH>>> cyclicOrderHalfFaces(const HEH& he, const HFH& hf, const FH& fJoin = FH()) const;

    /**
     * @brief Reconnect all left-behind patches to shifted arc by shifting them through block volume
     *
     * @param he IN: halfedge previously belonging to shifted arc
     * @param a IN: arc that was shifted
     * @param hf IN: halfface \p a was shifted across
     */
    void reconnectPatchesToArc(const HEH& he, const EH& a, const HFH& hf);

    /**
     * @brief Reconnect arc of \p he to shifted node by appending the halfedge that node was shifted across.
     *
     * @param he IN: halfedge belonging to arc previously connected to shifted node
     */
    void appendHalfedge(const HEH& arcHe);

    /**
     * @brief Get the tets in a dome around a central vertex (internal variable).
     *        Dome is bounded by block walls.
     *
     * @param tetSeed IN: tet to grow dome from
     * @return set<CH> tets in dome
     */
    set<CH> getDomeTets(const CH& tetSeed) const;

    /**
     * @brief Get the tets in the currently processing layer
     * @return set<CH> tets in current layer
     */
    set<CH> getCurrentLayerTets() const;

    /**
     * @brief Get the patch faces in the currently processing layer
     * @return set<FH> patch faces in current layer
     */
    set<FH> getCurrentLayerPatchFaces() const;

    /**
     * @brief Whether a patch on the boundary of this tet dome can be shifted across it
     *
     * @param domeElems IN: connectivity of the tet dome
     * @return true if shifting is feasible
     * @return false else
     */
    bool isShiftableDome(const DomeElements& domeElems) const;

    /**
     * @brief Get the connectivity of the dome by elements
     *
     * @param domeTets IN: tets forming the dome
     * @return DomeElements partitioned elements of the dome
     */
    DomeElements getDomeConnectivity(const set<CH>& domeTets) const;

    /**
     * @brief Find the path between two vertices within the ceiling of a tet dome
     *
     * @param from IN: from vertex
     * @param to IN: to vertex
     * @param domeTets IN: tets forming the dome
     * @return list<HEH> path connecting \p from to \p to
     */
    list<HEH> findPathOnDomeHull(const VH& from, const VH& to, set<CH>& domeTets);

    /**
     * @brief Find path between two vertices given a set of allowed vertices and edges
     *
     * @param from IN: from vertex
     * @param to IN: to vertex
     * @param esAllowed IN: allowed edges
     * @param vsAllowed IN: allowed vertices
     * @return list<HEH> path connecting \p from to \p to using only alllowed elements
     */
    list<HEH> findPath(const VH& from, const VH& to, const set<EH>& esAllowed, const set<VH>& vsAllowed) const;

    /**
     * @brief Determine the patch to take up the space of the opposite patch across a shifted arc
     *
     * @param pFirst IN: patch that shrinks
     * @param aShift IN: arc that shifts into shrinking patch
     * @return FH optimal patch to take up the space of \p pFirst
     */
    FH determineLastSuccessorPatch(const FH& pFirst, const EH& aShift) const;

    /**
     * @brief Determine arc to take up the space of the currently collapsing arc
     *
     * @return EH optimal patch to take up the space of collapsing arc
     */
    EH determineLastSuccessorArc() const;

    /**
     * @brief Collapse current arcs embedding edge by edge
     *
     * @return RetCode SUCCESS or error code
     */
    RetCode collapseEdgeByEdge();

    /**
     * @brief Collapse current patches embedding face by face
     *
     * @return RetCode SUCCESS or error code
     */
    RetCode collapseFaceByFace();

    /**
     * @brief Whether shift of an arc across a triangle described by the halfedges \p hfHes makes arc nonmanifold
     *
     * @param hesCurve IN: halfedges currently part of shifting arc
     * @param hfHes IN: halfedges of a triangle to shift arc across
     * @return true if resulting arc would be nonmanifold
     * @return false else
     */
    bool shiftMakesArcNonManifold(const set<HEH>& hesCurve, const list<HEH>& hfHes) const;

    /**
     * @brief Refine dome so a patch shift across it is possible without overlappings
     *
     * @param domeTets IN: tets in dome, OUT: after refinement and update
     * @param elems IN: dome elements, OUT: after refinement and update
     */
    void refineAndUpdateDome(set<CH>& domeTets, DomeElements& elems);

    /**
     * @brief Shift patch across tet dome
     *
     * @param domeTets IN: tets of dome
     * @param floor IN: floor halffaces of dome, patch starts here
     * @param ceiling IN: ceiling halffaces of dome, patch stops here
     */
    void patchDomeShift(const set<CH>& domeTets, const set<HFH>& floor, const set<HFH>& ceiling);

    /**
     * @brief Refine triangle fan so an arc shift across it is possible without overlaps
     *
     * @param he IN: Halfedge on one side of the fan. The arc sitting on this halfedge will be shifted
     * @param triangleFan IN: halffaces in the triangle fan, OUT: after refinement and update
     */
    void refineAndUpdateTriangleFan(const HEH& he, vector<HFH>& triangleFan);

    /**
     * @brief Refine a sector of tets incident on \p hePivot so a patch shift across it is possible without overlaps
     *
     * @param hePivot IN: edge around which sector unfolds
     * @param hfFlip IN: halfface bounding the sector from one side
     * @param hfPatch IN: halfface bounding the sector from the other side, has patch to be shifted
     * @param hfsBetween IN: fan of halffaces between the two boundaries
     */
    void refineAndUpdateTetFan(const HEH& hePivot, const HFH& hfFlip, HFH& hfPatch, list<HFH>& hfsBetween);

    // Patch faces in layer that is currently being processed (i.e. peeled off)
    // this will be only refreshed after a whole layer has been peeled off
    set<FH> _layerFaces;

    // Tets in layer that is currently being processed (i.e. peeled off)
    // this will be only refreshed after a whole layer has been peeled off
    set<CH> _layerTets;

    list<HEH> _collapseHes;   // Remaining halfedges of collapsing arc

    HEH _haCollapse;          // arc that is currently being collapsed (invalid if patch collapse)
    FH _pCollapse;            // patch that is currently being collapsed (invalid if arc collapse)
    HEH _haReroute;           // moving side of collapse patch
    HEH _haRerouteInto;       // stationary side of collapse patch
    HEH _heCurrent;           // the last halfedge across which the collapsing node was shifted
    VH _vPullUp;              // vertex to pull patches away from (only valid sometimes during patch collapse)

    EH _aToAppend;            // Arc to append to currently collapsing halfarc

    map<EH, FH> _a2pToAppend; // for each arc around the shifting node, which patch should take up the traversed space
};

} // namespace c4hex

#endif
