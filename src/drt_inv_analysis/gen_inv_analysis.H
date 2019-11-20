/*----------------------------------------------------------------------*/
/*! \file
\brief gen inv analysis

\level 2

\maintainer Harald Willmann
*/
/*----------------------------------------------------------------------*/
#ifndef GEN_INV_ANALYSIS_H
#define GEN_INV_ANALYSIS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_inpar/inpar_problemtype.H"

// forward declarations

namespace COMM_UTILS
{
  class NestedParGroup;
}  // namespace COMM_UTILS

namespace STR
{
  class TimInt;
}  // namespace STR

namespace IO
{
  class DiscretizationWriter;
}  // namespace IO

namespace LINALG
{
  class Solver;
  class SparseMatrix;
}  // namespace LINALG

namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace STR
{
  /*!
  \brief Set material parameters

  @param prob ??
  @param p_cur ??
  @param mymatset ??

  \note This routine is to be called also from supporting procs in case of nested parallel
  multi-scale inverse analysis.
  */
  void SetMaterialParameters(int prob, Epetra_SerialDenseVector& p_cur, std::set<int>& mymatset);

  /*!
  \brief Inverse Analysis of Structures

  To optimize the material parameters \f$ \mathbf p\f$ , given by

  \f[ \mathbf{p} = [p_1, ..., p_k]^T \f]

  of each constitutive model, an inverse analysis is performed.

  Beginning with an initial guess of the material parameters, displacements
  are calculated and compared to those obtained from the experiment.
  For the comparison, selected dofs of selected nodes, at every
  time step of the simulation are measured.

  As long as the difference between calculated and measured dofs
  are higher than a tolerance, a new set of parameters is calculated by
  utilizing the Levenberg- Marquardt Algorithm, developed by
  [Levenberg1944] and [Marquardt1963]. This algorithm is based
  on a least square fit, and has already shown good results for soft
  biological tissue [Kauer2002, Mahnken1996, Moulton1995, Seshaiyer2003].

  In this approach, the material parameters are updated according to
  \f[
  \Delta \mathbf{p} = - [\mathbf{J}^T(\mathbf {p})\mathbf{J}(\mathbf {p}) +
  \mu ~ diag(\mathbf{J}^T(\mathbf{p})\mathbf{J}(\mathbf{p}))]^{-1}\mathbf{J}^T(\mathbf{p})\mathbf{R}
  \f]

  with \f$\mu\f$ being a dimensionless training parameter and \f$\mathbf{J}\f$ being the
  rectangular Jacobian-Matrix of the dofs change with respect to the parameter change, i. e.

  \f[
  \begin{array}{c}
  \ \mathbf J= \left( \begin{array}{ccc}
  \frac {\partial \mathbf{u}_1} {\partial p_1} & ... & \frac {\partial \mathbf{u}_1}{\partial p_k}
  \\
  \vdots  &   &   \vdots \\
  \frac {\partial \mathbf{u}_{n}} {\partial p_1} &  ... & \frac {\partial
  \mathbf{u}_{n}}{\partial p_k} \end{array} \right). \end{array}
  \f]

  At each of the \f$ n\f$ time steps the dofs \f$\mathbf{u}\f$ were measured.
  \f$\mathbf{R}\f$ is the residual vector given by the difference between
  calculated and experimentally determined displacements
  \f[
         \begin{array}{c}
         \mathbf R =  \left( \begin{array}{c}
         \mathbf {r}_{1}\\
         \vdots \\
         \mathbf {r}_{n}
         \end{array} \right)
          =  \left( \begin{array}{c}
         \mathbf{u}_{1, cal} - \mathbf{u}_{1, exp}\\
         \vdots \\
         \mathbf{u}_{n,cal} - \mathbf{u}_{n, exp}
         \end{array} \right).
         \end{array}
  \f]

  The standardized error \f$ err\f$ used is the root mean square of
  \f$\mathbf R \f$
  \f[
  err =  \sqrt {\frac { \sum _{i=1}^{n} (\mathbf{r}_i\cdotp\mathbf{r}_i)} {n}}
  \f]

  \author rausch
  */
  class GenInvAnalysis
  {
   private:
    /// my structure discretization to work at
    Teuchos::RCP<DRT::Discretization> discret_;

    /// linear algebraic solver
    Teuchos::RCP<LINALG::Solver> solver_;

    /// context for output and restart
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    ///  number of measurement points the
    double tstep_;  // timestep

    /// number of fitted parameters (default 3, 4 with kappa)
    int np_;

    /// total number of measured dofs or points
    int nmp_;

    /// forward problem type given in dat-file FORWARD_PROBLEMTYP
    ProblemType forward_prb_type_;

    /// number of dimensions of forward problem
    int ndim_;

    //! @name Variables of multiple inverse analysis
    //! @{

    //! define if multiple inverse analysis is used
    bool spec_inv_ana_mult_;

    //! number of all measured dofs of every experiment
    int nmp_mult_;

    //! number of all experiments
    int nexp_;

    //! measured displacement of all the experiments (target value of multiple inverse analysis)
    Epetra_SerialDenseVector mcurve_mult_;

    //! array of nmp of each experiment
    Epetra_SerialDenseVector nmp_array_mult_;

    //! array of facnorms of each experiment
    Epetra_SerialDenseVector facnorm_mult_;

    //! @}

    //! @name Variables of coupled inverse analysis
    //! @{

    //! define if coupled inverse analysis is used
    bool spec_inv_ana_coup_;

    //! data points of pressure-volume-change experiment
    int nmp_volexp_;

    //! @}

    /// measured displacement of the experiments (target value)
    Epetra_SerialDenseVector mcurve_;

    /// time values where displacements are measured
    std::vector<double> timesteps_;

    //! Damping paramete for Levenberg-Marquardt in current iteration
    double mu_;

    //! Damping paramete for Levenberg-Marquardt in previous iteration
    double mu_o_;

    //! @name Parameter for convergence check of fitting iteration
    //! @{

    //! Error in current iteration
    double error_;

    //! Error in previous iteration
    double error_o_;

    /// Gradient based error measure in current iteration
    double error_grad_;

    /// Gradient based error measure in previous iteration
    double error_grad_o_;

    //! keeps the current error for the tolerance check depending on the update strategy
    double error_i_;

    //! User-given tolerance for convergence check of fitting iteration
    double tol_;

    //! @}

    /*! stores if inverse analysis does not longer converge (due to no descent direction) to stop
     * the inverse analysis
     */
    int nodescentdirection_;

    /// Update Strategy for mu
    enum INV_REG_UPDATE_STRATEGY
    {
      grad_based,  ///< gradient based update of regularization parameter
      res_based    ///< residual based update of regularization parameter
    };
    INV_REG_UPDATE_STRATEGY reg_update_;

    //! don't accept optimization step for negative material parameters
    int check_neg_params_;


    std::string filename_;

    double kappa_multi_;

    //! Counter for fitting iterations
    int numb_run_;

    //! parameter for finding starting point in #p_ of material parameters of coupexppol-material
    int coupexppol_parameter_position_;

    //! stores position of parameter alpha in #p_ of material viscofract
    int viscofract_alpha_position_;

    //! parameter vector
    Epetra_SerialDenseVector p_;

    //! parameters of the last step
    Epetra_SerialDenseVector p_o_;

    //! output parameters
    Epetra_SerialDenseVector p_print_;

    //! @name Storage
    //!@{

    Epetra_SerialDenseMatrix p_s_;
    Epetra_SerialDenseMatrix delta_p_s_;
    Epetra_SerialDenseMatrix par_s_;
    Epetra_SerialDenseVector mu_s_;
    Epetra_SerialDenseVector error_s_;
    Epetra_SerialDenseVector error_grad_s_;
    // Epetra_SerialDenseMatrix ccurve_s_;

    //! @}

    std::vector<std::set<int>> matset_;

    /// number of nodes or points observed
    int nnodes_;

    /// total number of dofs monitored
    int ndofs_;

    /// number of time steps observed
    int nsteps_;

    /// list of node gids observed
    std::vector<int> nodes_;

    /// list of dofs on each node that are monitored
    std::vector<std::vector<int>> dofs_;

    //! @name Patch stuff
    //! @{

    bool patches_;
    int numpatches_;
    int startindexpatches_;
    std::set<int> matset_patches_;
    Teuchos::RCP<LINALG::SparseMatrix> smoother_;
    int smoothingsteps_;

    enum DEFINE_PATCHES
    {
      material,
      uniform
    };
    DEFINE_PATCHES which_patches_;

    //! @}

    /// check whether we have to reset history parameters with this bool
    bool materialhashistory_;

    /// Strategy for evaluating distance from experiment to simulation
    enum MEASUREMENT_TYPE
    {
      dof_based,  /// distance measured based on individual node displacements
      point_based /*!< distance measured at specific points in prescribed directions, representing
                     the interface shape in experiment*/
    };
    MEASUREMENT_TYPE meas_type_;

   public:
    /// standard constructor
    GenInvAnalysis(Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
        Teuchos::RCP<LINALG::Solver> solver,               ///< linear solver
        Teuchos::RCP<IO::DiscretizationWriter> output      ///< output writer
    );

    /// Read in routine for the material parameters
    void ReadInParameters();

    /*!
    \brief Read the dof-based monitor file (default format)

    A dof-based monitor file provides measurement data that is linked to nodes of the finite element
    mesh and their respective degrees of freedom.

    <h3> Layout of the monitor file</h3>

    The monitor file has to be formatted as follows:

    \verbatim

    steps <numberOfTimeSteps> nnodes <numberOfNodesWithMeasurements>
    <nodeIdOne> <numberOfDofsWithMeasurements> <directionIndicators>
    <nodeIdTwo> <numberOfDofsWithMeasurements> <directionIndicators>
    ...
    <nodeIdN> <numberOfDofsWithMeasurements> <directionIndicators>
    # You can provide any number of comment lines here, but only at this position
    # After the comments, list the actual measurement data.
    # time | disp_1_x | disp_1_y | disp_2_x | disp_2_y | ... | disp_n_x | disp_n_y
    0.1   0.123   0.32165  1.254  1.3481
    0.2   0.224   0.564    2.8158 2.18

    \endverbatim

    with
    - <numberOfTimeSteps> denoting the number of time/load steps
    - <numberOfNodesWithMeasurements> denoting the total number of nodes where measurement data
    is available
    - <nodeIdOne>/<nodeIdTwo>/<nodeIdN> denoting the node ID from the Baci dat-file of the
      firsts/second/.../n-th node to be monitored
    - <numberOfDofsWithMeasurements> denotes the number of dimensions, where measurements are
    given at this node
    - <directionIndicators> is a list of directions (0 for x, 1 for y, 2 for z) where measurements
      are given

    After the comments, the data is given in column format. The first column lists the time/load
    step. Then, the first node lists its data in all directions in <directionIndicators>, followed
    by the second/third/... node.

    If measurements are only available for selected directions at a node, then only these directions
    are included into the <directionIndicators> list and the directions without measurement data
    are omitted from the <directionIndicators> list.

    <h3> Examples</h3>

    - Example for the file header of a 2D problem with x- and y-data at node 7:
        - <nodeID> = 7
        - <numberOfDofsWithMeasurements> = 2
        - <directionIndicators>  = 0 1

        \verbatim
        7 2 0 1
        \endverbatim

    - Example for the file header of a 3D problem with x-, y-, and z-data at node 11:
        - <nodeID> = 11
        - <numberOfDofsWithMeasurements> = 3
        - <directionIndicators>  = 0 1 2

        \verbatim
        11 3 0 1 2
        \endverbatim

    - Example for the file header of a 3D problem with only y-data at node 13:
        - <nodeID> = 13
        - <numberOfDofsWithMeasurements> = 1
        - <directionIndicators>  = 1

        \verbatim
        13 1 1
        \endverbatim

    @param[in] myrank Proc ID

    \note Constructor sequence has been moved into its own method.
    */
    void ReadMonitorDofBased(const int myrank);


    /// Read in routine for the point based monitor file
    void ReadMonitorPointBased(int myrank);

    /*! Set material parameters
     *
     * ToDo Pass in 'const Epetra_SerialDenseVector&' to avoid copy?
     *
     * @param p_cur
     */
    void SetParameters(Epetra_SerialDenseVector p_cur);

    /// get current state of calculated curve
    Epetra_SerialDenseVector GetCalculatedCurve(const Epetra_Vector& disp);

    /// get Distance of target points to respective interface at matching time
    Epetra_SerialDenseVector GetDistancePointsToInterfaceContour(
        const Epetra_Vector& disp, int writestep);

    /// coarse fitting of parameters
    void PrecondParameters();

    /*!
    \brief Fitting the material parameters using the Levenberg-Marquardt algorithm

    This implements the update of the parameters based on the Levenberg-Marquardt algorithm. At its
    core, we solve
    \[
    \left(J^T J + \mu \diag\left(j^T J\right)\right)\delta p
    = J^T r
    \]
    with the residual \f$r\f$ denoting the deviation of measurements and simulation results and its
    derivative \f$J\f$. The derivative \f$J\f$ is approximated via finite differencing.

    \note Since the number of parameters to be fitted, i.e. the size of \f$\delta p\f$, is low, this
    can be implemented with serial data structures and can be solved with a direct solver.

    @param cmatrix
    @param perturb
    */
    void CalcNewParameters(Epetra_SerialDenseMatrix& cmatrix, std::vector<double>& perturb);


    /// Calculate the new dof vector for comparison with the experimental data
    Epetra_SerialDenseVector CalcCvector(bool outputtofile);

    /// Calculate the new dof vector for comparison with the experimental data from a FSI problem
    Epetra_SerialDenseVector CalcCvectorFSI(bool outputtofile, int iLMstep);

    /// NP Calculate the new dof vector for comparison with the experimental data
    Epetra_SerialDenseVector CalcCvector(
        bool outputtofile, Teuchos::RCP<COMM_UTILS::NestedParGroup> group);

    /// actual inverse analysis loop
    void Integrate();

    /// additional calculations for integrate of multiple inverse analysis
    void NPIntegrateMult();

    /// acual inverse analysis loop (nested parallel version)
    void NPIntegrate();

    /*!
    \brief Output of Inverse Analysis

    ToDo Pass delta_p as 'const Epetra_SerialDenseVector&' to avoid copy?
    */
    void PrintStorage(Epetra_SerialDenseVector delta_p);
    void PrintFile();

    /// SetMaterialParameters for patch materials
    void SetMaterialParametersPatches(Epetra_SerialDenseVector& p_cur);

    /// initialize patches and prolongation operator
    void InitPatches();

    /// compute elementparameters from patch parameters
    void ComputeParametersFromPatches(
        Epetra_Vector p_patches, Teuchos::RCP<Epetra_Vector>& eleparameters);

    /// Check if mat params are out of bounds
    void CheckOptStep();

    /// Re-initialized micro-scale
    void MultiInvAnaInit();

    /// pressure-volume-change data of lung parenchyma
    void pVData(std::vector<double>& volexp_deltap, std::vector<double>& volexp_deltaV);

    /// compute pressure value of volume-change and material
    void ComputePressureNHIso1Coup3(std::vector<double> perturb, std::vector<double> volexp_deltaV,
        Epetra_SerialDenseMatrix& volexp_p_comp);

    /// check if pressure-volume-change-relation is physiological
    void CheckPhysiologicalPVRelationNHIso1Coup3();

    /// check if difference in dat-files of multiple inverse analysis occur
    void CheckDiffDat(int& expid);

    /// Constrain alpha of viscofract-material
    void ConstrainAlpha();


  };  // class GenInvAnalysis

}  // namespace STR

#endif /*GEN_INV_ANALYSIS_H_*/
