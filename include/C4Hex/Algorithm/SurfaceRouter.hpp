#ifndef C4HEX_SURFACEROUTER_HPP
#define C4HEX_SURFACEROUTER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages finding (geometric, i.e. tet mesh) surfaces spanning a given boundary
 */
class SurfaceRouter : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        NO_CONVERGENCE = 25, // No convergence to a valid surface
        NOT_CONNECTED = 26,  // the boundaries are not connected by any surface in the tet mesh
    };

    /**
     * @brief Create an instance that manages finding surfaces through the tet mesh associated with \p meshProps
     *
     * @param meshProps IN: tet mesh on which to find surfaces spanning boundaries
     */
    SurfaceRouter(TetMeshProps& meshProps);

    /**
     * @brief Given an initial surface \p surface in the tet mesh, try to iteratively shift it away from
     *        occupied elements so that it shares no elements except the boundary with other patches.
     *
     * @param b IN: block through which \p surface passes / should pass
     * @param surface IN: surface through \p b with boundary \p boundaryHes with some already-occupied elements,
     *                OUT: surface through \p b with boundary \p boundaryHes with no already-occupied elements except
     *                     boundary
     * @param transferredTets OUT: tets shifted across the moving boundary
     * @return RetCode SUCCESS or error code
     */
    RetCode rerouteSurfaceThroughBlock(const CH& b, set<HFH>& surface, set<CH>& transferredTets);

    /**
     * @brief Find the minimal discrete surface through a given block, connecting the given boundary
     *
     * @param b IN: block
     * @param surface IN: initial surface (must be a proper disk-like surface with ring-like boundary!), OUT: rerouted
     *                    surface
     * @param transferredTets OUT: tets transferred from one side of the surface to the other during reroute
     * @return RetCode SUCCESS or error code
     */
    RetCode minimalSurfaceThroughBlock(const CH& b, set<HFH>& surface, set<CH>& transferredTets);

    /**
     * @brief Incrementally shift the surface from the blocks boundary into its interior by transferring tets across.
     *
     * @param b IN: block
     * @param surface IN: initial surface touching the blocks boundary, OUT: final surface no longer touching the block
     * boundary
     * @param transferredTets OUT: tets transferred across the shifting surface
     * @return RetCode SUCCESS or error code
     */
    RetCode shiftSurfaceThroughBlock(const CH& b, set<HFH>& surface, set<CH>& transferredTets);

    /**
     * @brief Calculate (by linear program) a minimal discrete surface through a given connected volume, avoiding
     *        forbidden elements.
     *
     * @param space IN: connected volume, OUT: volume after refinement
     * @param forbiddenFs IN: forbidden faces
     * @param forbiddenEs IN: forbidden edges
     * @param forbiddenVs IN: forbidden vertices
     * @param boundaryEs IN: boundary edges of the surface
     * @param boundaryHes IN: boundary halfedges of the surface (inwards)
     * @param boundaryVs IN: boundary vertices of the surface
     * @param surfaceNew OUT: minimal discrete surface connecting the given boundary
     * @param depth IN: current attempt index. Internally after recursively refining, reconstraining and retrying 3
     *                  times without success, NO_CONVERGENCE is returned
     * @return RetCode SUCCESS or error code
     */
    RetCode calcMinimalSurfaceByLP(set<CH>& space,
                                   const set<FH>& forbiddenFs,
                                   const set<EH>& forbiddenEs,
                                   const set<VH>& forbiddenVs,
                                   const set<EH>& boundaryEs,
                                   const set<HEH>& boundaryHes,
                                   const set<VH>& boundaryVs,
                                   set<HFH>& surfaceNew,
                                   int depth = 0);

  protected:
    /**
     * @brief Non-manifoldness classificator
     */
    enum NonMF
    {
        NONMF_NONE = 0,
        NONMF_EDGE = 1,
        NONMF_VERTEX = 2
    };

    /**
     * @brief Restrict the given \p volume to a subvolume that does not exceed the axis-aligned bounding box of \p
     * vsBoundary
     *
     * @param vsBoundary IN: vertices on the boundary of a surface
     * @param volume IN: full volume
     * @return set<CH> subset of \p volume overlapping with the bounding box of \p vsBoundary
     */
    set<CH> confineVolume(const set<VH>& vsBoundary, const set<CH>& volume) const;

    /**
     * @brief Get the boundary elements of the given surface
     *
     * @param surface IN: discrete surface
     * @param vsBoundary OUT: boundary vertices
     * @param esBoundary OUT: boundary edges
     * @param hesBoundary OUT: boundary halfedges (inwards)
     */
    void
    getSurfaceBoundary(const set<HFH>& surface, set<VH>& vsBoundary, set<EH>& esBoundary, set<HEH>& hesBoundary) const;

    /**
     * @brief Gather the sets of elements (faces, edges, vertices) that are part of the
     *        current surface given as \p currentHfs . Also, count the number of current
     *        forbidden elements covered by \p currentHfs
     *
     * @param currentHfs IN: halffaces currently part of the (oriented) surface
     * @param currentEs OUT: edges currently part of the (oriented) surface
     * @param currentVs OUT: vertices currently part of the (oriented) surface
     * @param nForbiddenHfs OUT: number of forbidden halffaces currently part of the (oriented) surface
     * @param nForbiddenEs OUT: number of forbidden edges currently part of the (oriented) surface
     * @param nForbiddenVs OUT: number of forbidden vertices currently part of the (oriented) surface
     * @return RetCode SUCCESS or error code
     */
    RetCode gatherCurrentElements(const set<HFH>& currentHfs,
                                  set<EH>& currentEs,
                                  set<VH>& currentVs,
                                  int& nForbiddenHfs,
                                  int& nForbiddenEs,
                                  int& nForbiddenVs) const;

    /**
     * @brief Check whether \p tet is incident to any forbidden elements
     *
     * @param tet IN: tet for which the priority should be updated
     * @param vPrioTets IN: last tier with no subtiers
     * @return true if incident on forbidden element
     * @return false else
     */
    bool incidentOnForbiddenElement(const HFH& hf) const;

    /**
     * @brief Whether transferring a tet incident on the surface from one side of the surface to the other
     *        would create a non-manifold surface.
     *
     * @param tet IN: tet that should be transferred
     * @param currentHfs IN: halffaces currently part of the (oriented) surface
     * @param currentEs IN: edges currently part of the (oriented) surface
     * @param currentVs IN: vertices currently part of the (oriented) surface
     * @return true if a non-manifold surface would be created
     * @return false else
     */
    NonMF transferCreatesNonmanifold(const HFH& hf,
                                     const set<HFH>& currentHfs,
                                     const set<EH>& currentEs,
                                     const set<VH>& currentVs) const;

    /**
     * @brief If surface shifting can not proceed due to non-manifold intermediate configurations,
     *        this can be used to resolve such configurations, where the surface forms a pocket
     *        with triangular "entrance". The whole pocket (i.e. the contained tets) is first detected
     *        and then transferred to the other side of the surface, so that the triangular open entrance
     *        is "sealed" and the pocket is "removed".
     *
     * @param forbiddenFs IN: forbidden faces
     * @param forbiddenEs IN: forbidden edges
     * @param forbiddenVs IN: forbidden vertices
     * @param boundaryEs IN: prescribed boundary
     * @param frontHalfface IN: whether the front halffaces of the surface are inside the transferred tets
     * @param transferredTets IN: tets already shifted across the moving boundary OUT: input + shifted pocket tets
     * @param surfaceHfsToFlip IN/OUT: halffaces of the surface, which are still incident on forbidden elements
     * @param currentHfs IN/OUT: halffaces currently part of the (oriented) surface
     * @param currentEs IN/OUT: edges currently part of the (oriented) surface
     * @param currentVs IN/OUT: vertices currently part of the (oriented) surface
     * @param nForbiddenHfs OUT: number of forbidden halffaces currently part of the (oriented) surface
     * @param nForbiddenEs OUT: number of forbidden edges currently part of the (oriented) surface
     * @param nForbiddenVs OUT: number of forbidden vertices currently part of the (oriented) surface
     * @return true if a pocket could be found and resolved
     * @return false else
     */
    bool findAndTransferPocket(const set<EH>& boundaryEs,
                               bool frontHalfface,
                               set<CH>& transferredTets,
                               set<HFH>& surfaceHfsToFlip,
                               list<HFH>& hfsList,
                               set<HFH>& currentHfs,
                               set<EH>& currentEs,
                               set<VH>& currentVs,
                               int& nForbiddenHfs,
                               int& nForbiddenEs,
                               int& nForbiddenVs);

    /**
     * @brief Get the elements transferred between the surface and the rest of the mesh
     *        when transferring one tet incident on the surface from one side of the surface to the other
     *
     * @param tet IN: tet transferred from one side of the surface to the other
     * @param currentHfs IN: halffaces currently part of the (oriented) surface
     * @param frontHalfface IN: whether the front halffaces of the surface are inside the transferred tets
     * @param removedHfs OUT: halffaces removed from the surface
     * @param addedHfs OUT: halffaces added to the surface
     * @param removedEs OUT: edges removed from the surface
     * @param addedEs OUT: edges added to the surface
     * @param removedVs OUT: vertices removed from the surface
     * @param addedVs OUT: vertices added to the surface
     * @return RetCode SUCCESS or error code
     */
    RetCode exchangedElements(const CH& tet,
                              const set<HFH>& currentHfs,
                              bool frontHalfface,
                              set<HFH>& removedHfs,
                              set<HFH>& addedHfs,
                              set<EH>& removedEs,
                              set<EH>& addedEs,
                              set<VH>& removedVs,
                              set<VH>& addedVs) const;

    /**
     * @brief Refine block \p b to allow a reroute of a surface \p surface
     *        using any edges or vertices of _forbiddenEs or _forbiddenVs respectively.
     *
     * @param b IN: block
     * @param surface IN: initial surface
     * @return RetCode SUCCESS or error code
     */
    RetCode refineToAllowReroute(const CH& b, set<HFH>& surface);

    /**
     * @brief Refine the given volume so that a discrete surface connecting the given boundary not touching
     *        any forbidden elements is guaranteed to exist.
     *
     * @param space IN: allowed volume OUT: refined volume
     * @param boundaryVs IN: boundary vertices of surface
     * @param boundaryEs IN: boundary edges of surface
     * @return RetCode SUCCESS or error code
     */
    RetCode refineVolumeToAllowReroute(set<CH>& space, const set<VH>& boundaryVs, const set<EH>& boundaryEs);

    /**
     * @brief Refine the given volume so that there are more degrees of freedom around a set of complex elements
     *        from a previous solution of discrete minimal surface computation.
     *
     * @param space IN: allowed volume OUT: refined volume
     * @param complexEs IN: previously complex edges
     * @param complexVs IN: previously complex vertices
     * @param surface IN: previously computed surface with complex elements
     * @return RetCode SUCCESS or error code
     */
    RetCode refineVolumeToAvoidComplex(set<CH>& space,
                                       const set<EH>& complexEs,
                                       const set<VH>& complexVs,
                                       const set<HFH>& surface);

    /**
     * @brief Refine the tet incident to \p hf to avoid creating a nonmanifold vertex when transferring the tet
     *
     * @param hf IN: halfface for which the incident tet should be transferred across the surface
     * @param forbiddenVs IN: forbidden vertices
     * @return RetCode SUCCESS or error code
     */
    RetCode refineToAvoidNonmanifoldVertex(const HFH& hf, set<CH>* space = nullptr);

    /**
     * @brief Refine the tet incident to \p hf to avoid creating a nonmanifold edge when transferring the tet
     *
     * @param tet IN: tet
     * @param currentHfs IN: current surface
     * @param forbiddenVs IN: forbidden vertices
     * @param forbiddenEs IN: forbidden edges
     * @return RetCode SUCCESS or error code
     */
    RetCode refineToAvoidNonmanifoldEdge(const CH& tet, const set<HFH>& currentHfs, set<CH>* space = nullptr);

    /**
     * @brief Internally set all forbidden elements from the block's boundary
     *
     * @param b IN: block whose boundary should be set to forbidden
     */
    void setForbiddenElements(const CH& b);

    set<VH> _forbiddenVs;
    set<EH> _forbiddenEs;
    set<FH> _forbiddenFs;
};

} // namespace c4hex

#endif
