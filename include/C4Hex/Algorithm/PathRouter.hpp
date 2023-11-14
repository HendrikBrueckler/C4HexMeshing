#ifndef C4HEX_PATHROUTER_HPP
#define C4HEX_PATHROUTER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages finding (tet-mesh embedded) paths between two vertices constrained to some surface or
 *        volume. Can refine the tet mesh if needed.
 */
class PathRouter : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        NOT_CONNECTED = 26, // No path between given elements exists in the mesh
    };

    /**
     * @brief Create an instance that manages finding paths on the tet mesh associated with \p meshProps
     *
     * @param meshProps IN: tet mesh
     */
    PathRouter(TetMeshProps& meshProps);

    /**
     * @brief Reroute a path within a given surface, defined by a set of triangles, so that it no
     *        longer coincides with any elements marked forbidden.
     *
     * @param path IN: initial path, OUT: rerouted path
     * @param surface IN: surface to which path finding is confined, OUT: possibly refined input surface
     * @param forbiddenFs IN: faces marked as forbidden OUT: possibly refined input
     * @param forbiddenEs IN: edges marked as forbidden OUT: possibly refined input
     * @param forbiddenVs IN: vertices marked as forbidden OUT: possibly refined input
     * @param hfsTransferred OUT: halffaces of \p surface traversed by path rerouting
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode reroutePathThroughSurface(list<HEH>& path,
                                      set<FH>& surface,
                                      set<FH>& forbiddenFs,
                                      set<EH>& forbiddenEs,
                                      set<VH>& forbiddenVs,
                                      set<HFH>& hfsTransferred);

    /**
     * @brief Find the shortest path between two vertices within a set of allowed faces.
     *
     * @param vFrom IN: from vertex
     * @param vTo IN: to vertex
     * @param path OUT: path
     * @param surface IN: surface within which to route OUT: possibly refined input surface
     * @param forbiddenFs IN: explicitly forbidden faces
     * @param forbiddenEs IN: explicitly forbidden edges
     * @param forbiddenVs IN: explicitly forbidden vertices
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode shortestPathThroughSurface(const VH& vFrom,
                                       const VH& vTo,
                                       list<HEH>& path,
                                       set<FH>& surface,
                                       set<FH>& forbiddenFs,
                                       set<EH>& forbiddenEs,
                                       set<VH>& forbiddenVs);

    /**
     * @brief Reroute a path within a given volume, defined by a set of tets, so that it no
     *        longer coincides with any elements marked forbidden.
     * @param path IN: initial path, OUT: rerouted path
     * @param volume IN: volume to which path finding is confined, OUT: possibly refined input volume
     * @param forbiddenFs IN: faces marked as forbidden OUT: possibly refined input
     * @param forbiddenEs IN: edges marked as forbidden OUT: possibly refined input
     * @param forbiddenVs IN: vertices marked as forbidden OUT: possibly refined input
     * @param hfsTransferred OUT: halffaces of \p volume traversed by path rerouting
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode reroutePathThroughVolume(list<HEH>& path,
                                     set<CH>& volume,
                                     set<FH>& forbiddenFs,
                                     set<EH>& forbiddenEs,
                                     set<VH>& forbiddenVs,
                                     set<HFH>* hfsTransferred = nullptr);

    /**
     * @brief Find the shortest path between two vertices within a set of allowed tets.
     *
     * @param vFrom IN: from vertex
     * @param vTo IN: to vertex
     * @param path OUT: path
     * @param volume IN: volume within which to route OUT: possibly refined input volume
     * @param forbiddenFs IN: explicitly forbidden faces
     * @param forbiddenEs IN: explicitly forbidden edges
     * @param forbiddenVs IN: explicitly forbidden vertices
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode shortestPathThroughVolume(const VH& vFrom,
                                      const VH& vTo,
                                      list<HEH>& path,
                                      set<CH>& volume,
                                      set<FH>& forbiddenFs,
                                      set<EH>& forbiddenEs,
                                      set<VH>& forbiddenVs);

    /**
     * @brief Reroute a path incident on forbidden (already occupied) edges
     *        so that is no longer incident on any forbidden elements.
     *        Mesh is refined in the process.
     *
     * @param p IN: patch through which the path is to be routed
     * @param path IN: initial invalid path, OUT: rerouted path
     * @param hfsTransferred OUT: halffaces traversed by path reroute
     * @param minimize IN: whether to minimize path length
     * @return RetCode SUCCESS or error code
     */
    RetCode reroutePathThroughPatch(const FH& p, list<HEH>& path, set<HFH>& hfsTransferred);

  protected:
    /**
     * @brief Determine the segments of \p pathRerouted that branch off of the original \p path (i.e. do not coincide
     * with)
     *
     * @param path IN: original path
     * @param pathRerouted IN: rerouted path
     * @param branchesPath OUT: list of segments where of \p path that do not coincide with \p pathRerouted
     * @param branchesPathRerouted OUT: list of segments where of \p pathRerouted that do not coincide with \p path
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode determineBranches(const list<HEH>& path,
                              const list<HEH>& pathRerouted,
                              list<list<HEH>>& branchesPath,
                              list<list<HEH>>& branchesPathRerouted);

    /**
     * @brief Get all halffaces enclosed between the path branches within a given surface
     *
     * @param branchesPath IN: branch curves of original path
     * @param branchesPathRerouted IN: branch curves of rerouted path
     * @param surface IN: surface within which rerouting was done
     * @param hfsEnclosed OUT: halffaces crossed during rerouting (= halffaces between branches)
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode getEnclosedHalffaces(const list<list<HEH>>& branchesPath,
                                 const list<list<HEH>>& branchesPathRerouted,
                                 const set<FH>& surface,
                                 set<HFH>& hfsEnclosed) const;

    /**
     * @brief Get all halffaces enclosed between the path branches within a given volume
     *
     * @param branchesPath IN: branch curves of original path
     * @param branchesPathRerouted IN: branch curves of rerouted path
     * @param volume IN: volume within which rerouting was done OUT: refined
     * @param hfsEnclosed OUT: halffaces crossed during rerouting (= halffaces between branches)
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode getEnclosedHalffaces(const list<list<HEH>>& branchesPath,
                                 const list<list<HEH>>& branchesPathRerouted,
                                 set<CH>& volume,
                                 set<HFH>& hfsEnclosed);

    /**
     * @brief Refine surface to guarantee the existence of an edge-embedded path within allowed surface.
     *
     * @param surface IN: surface within which reroute is confined OUT: refined surface
     * @param vFrom IN: from vertex
     * @param vTo IN: to vertex
     * @return RetCode SUCCESS
     */
    RetCode refineSurfaceToAllowReroute(set<FH>& surface, const VH& vFrom, const VH& vTo);

    /**
     * @brief Refine volume to guarantee the existence of an edge-embedded path within allowed volume.
     *
     * @param volume IN: volume within which reroute is confined OUT: refined volume
     * @param vFrom IN: from vertex
     * @param vTo IN: to vertex
     * @return RetCode SUCCESS
     */
    RetCode refineVolumeToAllowReroute(set<CH>& volume, const VH& vFrom, const VH& vTo);

    /**
     * @brief A star algorithm to find shortest path within sets of allowed edges and vertices
     *
     * @param vFrom IN: from vertex
     * @param vTo IN: to vertex
     * @param path OUT: shortest path
     * @param allowedEs IN: all edges allowed to use
     * @param allowedVs IN: all vertices allowed to use
     * @return RetCode SUCCESS or NOT_CONNECTED
     */
    RetCode aStarShortestPath(const VH& vFrom, const VH& vTo, list<HEH>& path, set<EH>& allowedEs, set<VH>& allowedVs);

    /**
     * @brief Refine patch \p p to allow a reroute of a path between \p vFrom and \p vTo without
     *        using any edges or vertices of _forbiddenEs or _forbiddenVs respectively.
     *
     * @param p IN: patch to refine
     * @param vFrom IN: path start
     * @param vTo IN: path end
     * @return RetCode SUCCESS
     */
    RetCode refineToAllowReroute(const FH& p, const VH& vFrom, const VH& vTo);

    /**
     * @brief Given a set of halfedges \p currentHes which form the current path, get
     *        the corresponding set of vertices and count the forbidden elements along the path.
     *
     * @param currentHes IN: path halfedges
     * @param currentVs OUT: path vertices
     * @param nForbiddenEs OUT: number of forbidden path edges
     * @param nForbiddenVs OUT: number of forbidden path vertices
     * @return RetCode SUCCESS
     */
    RetCode
    gatherCurrentElements(const set<HEH>& currentHes, set<VH>& currentVs, int& nForbiddenEs, int& nForbiddenVs) const;

    /**
     * @brief Check whether halfedge \p he is incident on a forbidden element on the path
     *
     * @param he IN: halfedge
     * @return true if incident on forbidden element
     * @return false else
     */
    bool incidentOnForbiddenElement(const HEH& he) const;

    /**
     * @brief Whether the transfer of the halfface incident on the given halfedge
     *        from one side of the path to the other is recommended under the path
     *        length minimization heuristic.
     *
     * @param he IN: halfedge
     * @param removedHes IN: halfedges to be removed from path by transfer
     * @param addedHes IN: halfedges to be added to path by transfer
     * @param removedVs IN: vertices to be removed from path by transfer
     * @param addedVs IN: vertices to be added to path by transfer
     * @return true if recommended
     * @return false else
     */
    bool rerouteRecommended(const HEH& he,
                            const list<HEH>& removedHes,
                            const list<HEH>& addedHes,
                            const set<VH>& removedVs,
                            const set<VH>& addedVs) const;

    /**
     * @brief Gather the elements added to or removed from the path by transferring \p hf to the other side
     *        of the path
     *
     * @param hf IN: halfface
     * @param currentHes IN: path halfedges
     * @param removedHes OUT: halfedges removed from path
     * @param addedHes OUT: halfedges added to path
     * @param removedVs OUT: vertices removed from path
     * @param addedVs OUT: vertices added to path
     * @return RetCode SUCCESS
     */
    RetCode exchangedElements(const HFH& hf,
                              const set<HEH>& currentHes,
                              list<HEH>& removedHes,
                              list<HEH>& addedHes,
                              set<VH>& removedVs,
                              set<VH>& addedVs) const;

    /**
     * @brief Get the halfface transferred when removing the given halfedge
     *
     * @param pathHe IN: halfedge on path to be removed (and replaced)
     * @return HFH halfface transferred to other path side during replacement
     */
    HFH getTransferredHf(const HEH& pathHe) const;

    /**
     * @brief Execute the path reroute given by precomputed element transfers/toggles.
     *
     * @param pathRerouted IN: path, OUT: path after single reroute operation
     * @param pathIt IN: iterator to current path edge, OUT: iterator to current path edge after reroute
     * @param hfsTransferred IN: halffaces previously transferred, OUT: updated after reroute
     * @param currentVs IN: vertices previously part of path, OUT: updated after reroute
     * @param currentHes IN: halfedges previously part of path, OUT: updated after reroute
     * @param vsPassed IN: vertices previously visited, OUT: updated after reroute
     * @param removedHes IN: halfedges to be removed from path
     * @param addedHes IN: halfedges to be added to path
     * @param removedVs IN: vertices to be removed from path
     * @param addedVs IN: vertices to be added to path
     */
    void reroute(list<HEH>& pathRerouted,
                 list<HEH>::iterator& pathIt,
                 set<HFH>& hfsTransferred,
                 set<VH>& currentVs,
                 set<HEH>& currentHes,
                 set<VH>& vsPassed,
                 list<HEH>& removedHes,
                 list<HEH>& addedHes,
                 set<VH>& removedVs,
                 set<VH>& addedVs);

  private:
    bool _inFirstHp;
    FH _p;

    set<FH> _forbiddenFs;
    set<EH> _forbiddenEs;
    set<VH> _forbiddenVs;
};

} // namespace c4hex

#endif
