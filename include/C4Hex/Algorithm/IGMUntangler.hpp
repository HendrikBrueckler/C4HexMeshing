#ifndef C4HEX_PARAMOPTIMIZER_HPP
#define C4HEX_PARAMOPTIMIZER_HPP

#include <MC3D/Mesh/MCMeshManipulator.hpp>

namespace c4hex
{
using namespace mc3d;

/**
 * @brief Class that manages IGM untangling
 */
class IGMUntangler : public virtual MCMeshManipulator
{
  public:
    enum RetCode
    {
        SUCCESS = 0,
        NUMERICAL_ISSUE = 23,
        NO_CONVERGENCE = 24,
    };

    /**
     * @brief Create an instance that manages the untangling of the IGM associated with \p meshProps
     *
     * @param meshProps IN/OUT: mesh whose IGM should be optimized
     */
    IGMUntangler(TetMeshProps& meshProps);

    /**
     * @brief Optimize the IGM of the passed mesh using relative weighting of \p areaVsAngles to enforce area/angle
     *        preservation. This will incrementally split fully-constrained edges adjacent to inverted tets, that
     *        are not untangleable.
     *
     * @param areaVsAngles IN: factor in [0, 1] determining the weight of area/volume preservation relative to the
     *                         weight of angle preservation
     * @param maxIter IN: maximum inner iterations
     * @param secondsTimeLimit IN: maximum time per block and cycle (multiple cycles happen only when edges are split)
     * @return RetCode SUCCESS, NUMERICAL_ISSUE or NO_CONVERGENCE
     */
    RetCode untangleIGM(double areaVsAngles, int maxIter, int secondsTimeLimit = 300);

    /**
     * @brief Reset optimization history. Specifically this allows to reuse methods that have been discarded because of
     * previous failure
     *
     */
    void reset();

  private:
    /**
     * @brief Determine some stats for IGM in given block
     *
     * @param b IN: block
     * @param vol OUT: IGM volume of \p b
     * @param volInverted OUT: inverted IGM volume of \p b
     * @param minVol OUT: minimal IGM volume of any tet of \p b
     * @param tetsInverted OUT: inverted tets of \p b
     */
    void determineIGMStats(const CH& b, double& vol, double& volInverted, double& minVol, set<CH>& tetsInverted) const;

    /**
     * @brief Index vertices, set the coordinates of vertices on block boundaries fixed and return
     *        tet-vertex indices
     *
     * @param b IN: block
     * @param v2corner2idx OUT: mapping of block vertices/tetcorners to their array/vector index
     * @param lb OUT: lower bounds (block/patch/arc bounds) for each vertex coordinate
     * @param ub OUT: upper bounds (block/patch/arc bounds) for each vertex coordinate
     * @param nInteriorVs OUT: Number of interior vertices
     * @return Eigen::Matrix4Xi tet-vertex indices
     */
    Eigen::Matrix4Xi fixBlockVertices(const CH& b,
                                      map<VH, map<CH, int>>& v2corner2idx,
                                      Eigen::VectorXd& lb,
                                      Eigen::VectorXd& ub,
                                      int& nInteriorVs);

    /**
     * @brief Compute the Z (shape) matrix and volume for each tet
     *
     * @param tetVtxIndices IN: indices of each tets vertices
     * @param XYZflat IN: xyz coordinates
     * @param UVWflat IN: initial UVW coordinates
     * @param tetZ OUT: Z matrix per tet
     * @param tetVol OUT: tet volume per tet
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    static RetCode precompute(const Eigen::Matrix4Xi& tetVtxIndices,
                              const Eigen::VectorXd XYZflat,
                              const Eigen::VectorXd UVWflat,
                              Eigen::Matrix4Xd& tetZ,
                              Eigen::VectorXd& tetVol);

    /**
     * @brief Update jacobian stats for each tet
     *
     * @param tetVtxIndices IN: indices of each tets vertices
     * @param UVWflat IN: current UVW coordinates
     * @param tetZ IN: Z matrix per tet
     * @param tetJ OUT: J matrix per tet
     * @param tetDetJ OUT: determinant per tet
     * @param minDetJ OUT: minimum determinant
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    static RetCode updateJ(const Eigen::Matrix4Xi& tetVtxIndices,
                           const Eigen::VectorXd UVWflat,
                           const Eigen::Matrix4Xd& tetZ,
                           Eigen::Matrix3Xd& tetJ,
                           Eigen::VectorXd& tetDetJ,
                           double& minDetJ,
                           int& nFlipped);

    /**
     * @brief Calculate the foldover energy F
     *
     * @param tetVtxIndices IN: indices of each tets vertices
     * @param tetVol IN: tet volume per tet
     * @param tetJ IN: J matrix per tet
     * @param tetDetJ IN: determinant per tet
     * @param e IN: epsilon value (check the paper)
     * @param areaOverAngles IN: weighting of area vs angle preservation
     * @param F OUT: total energy
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    static RetCode calcF(const Eigen::Matrix4Xi& tetVtxIndices,
                         const Eigen::VectorXd& tetVol,
                         const Eigen::VectorXd& weight,
                         const Eigen::Matrix3Xd& tetJ,
                         const Eigen::VectorXd& tetDetJ,
                         double e,
                         double areaOverAngles,
                         double& F);

    /**
     * @brief Update the gradient of foldover energy function F
     *
     * @param tetVtxIndices IN: indices of each tets vertices
     * @param tetVol IN: tet volume per tet
     * @param tetZ IN: Z matrix per tet
     * @param tetJ IN: J matrix per tet
     * @param tetDetJ IN: determinant per tet
     * @param nInteriorVs IN: Number of interior vertices
     * @param e IN: epsilon value (check the paper)
     * @param areaOverAngles IN: weighting of area vs angle preservation
     * @param gradF OUT: gradient of foldover energy function F
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    static RetCode updateGradF(const Eigen::Matrix4Xi& tetVtxIndices,
                               const Eigen::VectorXd& tetVol,
                               const Eigen::VectorXd& weight,
                               const Eigen::Matrix4Xd& tetZ,
                               const Eigen::Matrix3Xd& tetJ,
                               const Eigen::VectorXd& tetDetJ,
                               const int nInteriorVs,
                               double e,
                               double areaOverAngles,
                               Eigen::VectorXd& gradF);

    /**
     * @brief Update the hessian of foldover energy function F
     *
     * @param tetVtxIndices IN: indices of each tets vertices
     * @param tetVol IN: tet volume per tet
     * @param tetZ IN: Z matrix per tet
     * @param tetJ IN: J matrix per tet
     * @param tetDetJ IN: determinant per tet
     * @param nInteriorVs IN: Number of interior vertices
     * @param e IN: epsilon value (check the paper)
     * @param areaOverAngles IN: weighting of area vs angle preservation
     * @param hessF OUT: hessian matrix of foldover energy function F
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    static RetCode updateHessF(const Eigen::Matrix4Xi& tetVtxIndices,
                               const Eigen::VectorXd& tetVol,
                               const Eigen::VectorXd& weight,
                               const Eigen::Matrix4Xd& tetZ,
                               const Eigen::Matrix3Xd& tetJ,
                               const Eigen::VectorXd& tetDetJ,
                               const int nInteriorVs,
                               double e,
                               double areaOverAngles,
                               Eigen::SparseMatrix<double>& hessF);

    /**
     * @brief Update the current UVW values by one of: gradient descent, L-BFGS or Newton
     *
     *
     * @param tetVtxIndices IN: indices of each tets vertices
     * @param tetVol IN: tet volume per tet
     * @param gradF IN: gradient of foldover energy function F
     * @param hessF IN: hessian matrix of foldover energy function F
     * @param tetZ IN: Z matrix per tet
     * @param nInteriorVs IN: Number of interior vertices
     * @param lowerBounds IN: lower block/patch/arc bounds for clamping values
     * @param upperBounds IN: upper block/patch/arc bounds for clamping values
     * @param e IN: epsilon value (check the paper)
     * @param areaOverAngles IN: weighting of area vs angle preservation
     * @param tetJ IN: J matrix per tet
     * @param tetDetJ IN: determinant per tet
     * @param UVWflat IN: current UVW coordinates, OUT: updated UVW coordinates
     * @param iterNoImprovement IN: number of outer iterations without improvement OUT: new number of outer iteration
     * without improvement
     * @param mode IN: 0: L-BFGS, 1: gradient descent, 2: Newton, 3: L-BFGS with more inner iterations
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    static RetCode updateUVW(const Eigen::Matrix4Xi& tetVtxIndices,
                             const Eigen::VectorXd& tetVol,
                             const Eigen::VectorXd& weight,
                             const Eigen::VectorXd& gradF,
                             const Eigen::SparseMatrix<double>& hessF,
                             const Eigen::Matrix4Xd& tetZ,
                             const int nInteriorVs,
                             const Eigen::VectorXd& lowerBounds,
                             const Eigen::VectorXd& upperBounds,
                             double e,
                             double areaOverAngles,
                             Eigen::Matrix3Xd& tetJ,
                             Eigen::VectorXd& tetDetJ,
                             Eigen::VectorXd& UVWflat,
                             int& iterNoImprovement,
                             int mode);

    /**
     * @brief Update the epsilon value
     *
     * @param minDetJ IN: minimum determinant
     * @param lastF IN: energy before UVW update
     * @param F IN: energy after UVW update
     * @param lastE IN: previous epsilon
     * @param e OUT: new epsilon
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    static RetCode updateE(double minDetJ, double lastF, double F, double lastE, double& e);

    /**
     * @brief Apply the untangled local UVW values to the parametrization mesh property
     *
     * @param vtx2index IN: mapping of block vertex to array/matrix inde
     * @param nInteriorVs IN: Number of interior vertices
     * @param UVWflat IN: untangled UVW values per vertex
     */
    void applyUntangling(const map<VH, map<CH, int>>& v2corner2idx,
                         const int nInteriorVs,
                         const Eigen::VectorXd& UVWflat);

    /**
     * @brief Optimize the IGM of the passed block (inner vertices only) using total lifted content method
     *
     * @param b IN: block to untangle
     * @param maxIter IN: max number of iterations
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    RetCode untangleByTLC(const CH& b, int maxIter, int secondsTimeLimit);

    /**
     * @brief Optimize the IGM of the passed block (inner vertices only) using foldover-free maps method
     *
     * @param b IN: block to untangle
     * @param areaVsAngles IN: factor in [0, 1] determining the weight of area/volume preservation relative to the
     *                         weight of angle preservation
     * @param maxIter IN: max number of iterations
     * @return RetCode SUCCESS or NUMERICAL_ISSUE
     */
    RetCode untangleByFFM(const CH& b, double areaVsAngles, int maxIter, int secondsTimeLimit);

    set<CH> _blocksOptimizedByTLC; // Blocks already untangled by TLC method (to avoid retrying in case of failure)

    // Interface for L-BFGS
    class FoldoverEnergy
    {
      public:
        FoldoverEnergy(double areaVsAngles_,
                       double e_,
                       const Eigen::Matrix4Xi& tetVtxIndices_,
                       const Eigen::Matrix4Xd& tetZ_,
                       const Eigen::VectorXd& tetVol_,
                       const int nInteriorVs,
                       const Eigen::VectorXd& tetWeights_);

        double operator()(const Eigen::VectorXd& UVWflat, Eigen::VectorXd& gradF);

      private:
        bool updateJ(const Eigen::VectorXd& UVWflat);

        bool updateGradF(Eigen::VectorXd& gradF);

        double calcF();

        double areaVsAngles;
        double minDetJ;
        int nFlipped;
        double e;
        const Eigen::Matrix4Xi& tetVtxIndices;
        const Eigen::Matrix4Xd& tetZ;
        const Eigen::VectorXd& tetVol;
        const int nInteriorVs;
        Eigen::Matrix3Xd tetJ;
        Eigen::VectorXd tetDetJ;
        Eigen::VectorXd tetWeights;
    };
};

} // namespace c4hex

#endif
