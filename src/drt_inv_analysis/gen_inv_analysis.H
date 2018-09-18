/*----------------------------------------------------------------------*/
/*!
 * \file gen_inv_analysis.H

\brief gen inv analysis

\level 2

\maintainer Harald Willmann
            willmann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
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
}

namespace STR
{
  class TimInt;
}

namespace IO
{
  class DiscretizationWriter;
}

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
  /*====================================================================*/
  /*!
   * Inverse Analysis of Structures
   *
   * To optimize the material parameters \f$ \mathbf p\f$ , given by
   *
   * \f[ \mathbf{p} = [p_1, ..., p_k]^T \f]
   *
   * of each constitutive model, an inverse analysis was performed.
   *
   *
   *
   * Beginning with an initial guess of the material parameters, displacements
   * were calculated and compared to those obtained from the experiment.
   * For the comparison, selected dofs of selected nodes, at every
   * time step of the simulation were measured.
   *
   * As long as the difference between calculated and measured dofs
   * was higher than a tolerance, a new set of parameters was calculated by
   * utilizing the Levenberg- Marquardt Algorithm, developed by
   * [Levenberg1944] and [Marquardt1963]. This algorithm is based
   * on a least square fit, and has already shown good results for soft
   * biological tissue [Kauer2002, Mahnken1996, Moulton1995, Seshaiyer2003].
   *
   * In this approach the material parameters were updated according to
   * \f[
   *       \Delta \mathbf{p} = [\mathbf{J}^T(\mathbf {p})\mathbf{J}(\mathbf {p}) +
   *       \mu ~ diag(\mathbf{J}^T(\mathbf {p})\mathbf{J}(\mathbf {p}))]^{-1} \mathbf{J}^T(\mathbf
   *{p}) \mathbf {R} \f]
   *
   * with \f$\mu\f$ being a dimensionless training parameter and \f$\mathbf{J}\f$ being the
   *rectangular Jacobian-Matrix of the dofs change with respect to the parameter change, i. e.
   *
   * \f[
   *     \begin{array}{c}
   *     \ \mathbf J= \left( \begin{array}{ccc}
   *     \frac {\partial \mathbf{u}_1} {\partial p_1} & ... & \frac {\partial \mathbf{u}_1}{\partial
   *p_k} \\
   *     \vdots  &   &   \vdots \\
   *     \frac {\partial \mathbf{u}_{n}} {\partial p_1} &  ... & \frac {\partial
   *\mathbf{u}_{n}}{\partial p_k} \end{array} \right). \end{array} \f]
   *
   * At each of the \f$ n\f$ time steps the dofs \f$\mathbf{u}\f$ were measured.
   * \f$\mathbf{R}\f$ is the residual vector given by the difference between
   * calculated and experimentally determined displacements
   *\f[
   *        \begin{array}{c}
   *        \mathbf R =  \left( \begin{array}{c}
   *        \mathbf {r}_{1}\\
   *        \vdots \\
   *        \mathbf {r}_{n}
   *        \end{array} \right)
   *         =  \left( \begin{array}{c}
   *        \mathbf{u}_{1, cal} - \mathbf{u}_{1, exp}\\
   *        \vdots \\
   *        \mathbf{u}_{n,cal} - \mathbf{u}_{n, exp}
   *        \end{array} \right).
   *        \end{array}
   *\f]
   *
   * The standardized error \f$ err\f$ used is the root mean square of
   * \f$\mathbf R \f$
   *\f[
   *     err =  \sqrt {\frac { \sum _{i=1}^{n} (\mathbf{r}_i\cdotp\mathbf{r}_i)} {n}}
   *\f]
   *
   *
   * \author rausch
   */

  /// Set material parameters (to be called also from supporting procs in case
  /// of nested parallel multi-scale inverse analysis)
  void SetMaterialParameters(int prob, Epetra_SerialDenseVector& p_cur, std::set<int>& mymatset);


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

    int problem_type_;  // Neumann or Dirichlet Problem
    int np_;            // number of fitted parameters (default 3, 4 with kappa)
    int nmp_;           // total number of measured dofs

    /// forward problem type given in dat-file FORWARD_PROBLEMTYP
    PROBLEM_TYP forward_prb_type_;
    /// number of dimensions of forward problem
    int ndim_;

    /// variables of multiple inverse analysis
    bool spec_inv_ana_mult_;                // define if multiple inverse analysis is used
    int nmp_mult_;                          // number of all measured dofs of every experiment
    int nexp_;                              // number of all experiments
    Epetra_SerialDenseVector mcurve_mult_;  // measured displacement of all the experiments (target
                                            // value of multiple inverse analysis)
    Epetra_SerialDenseVector nmp_array_mult_;  // array of nmp of each experiment
    Epetra_SerialDenseVector facnorm_mult_;    // array of facnorms of each experiment

    /// variables of coupled inverse analysis
    bool spec_inv_ana_coup_;  // define if coupled inverse analysis is used
    int nmp_volexp_;          // data points of pressure-volume-change experiment

    Epetra_SerialDenseVector mcurve_;  // measured displacement of the experiments (target value)
    std::vector<double> timesteps_;    // time values where displacements are measured

    double mu_;
    double mu_o_;

    double error_;
    double error_o_;
    // Gradient based error measure
    double error_grad_;
    double error_grad_o_;

    double error_i_;  // keeps the current error for the tolerance check depending on the update
                      // strategy
    double tol_;
    int nodescentdirection_;  // stores if inverse analysis does not longer
    // converge (due to no descent direction) to stop the inverse analysis

    // Update Strategy for mu
    enum INV_REG_UPDATE_STRATEGY
    {
      grad_based,
      res_based
    };
    INV_REG_UPDATE_STRATEGY reg_update_;

    // don't accept optimization step for negative material parameters
    int check_neg_params_;


    std::string filename_;

    double kappa_multi_;

    int numb_run_;

    // parameter for finding starting point in p_ of material parameters of coupexppol-material
    int coupexppol_parameter_position_;
    // stores position of parameter alpha in p_ of material viscofract
    int viscofract_alpha_position_;

    Epetra_SerialDenseVector p_;        // parameter vector
    Epetra_SerialDenseVector p_o_;      // parameters of the last step
    Epetra_SerialDenseVector p_print_;  // output parameters

    /// @name Storage
    //@{
    Epetra_SerialDenseMatrix p_s_;
    Epetra_SerialDenseMatrix delta_p_s_;
    Epetra_SerialDenseMatrix par_s_;
    Epetra_SerialDenseVector mu_s_;
    Epetra_SerialDenseVector error_s_;
    Epetra_SerialDenseVector error_grad_s_;
    // Epetra_SerialDenseMatrix ccurve_s_;
    //@}

    std::vector<std::set<int>> matset_;

    // observed nodes, dofs, timesteps
    int nnodes_;                          // number of nodes observed
    int ndofs_;                           // total number of dofs monitored
    int nsteps_;                          // number of time steps observed
    std::vector<int> nodes_;              // list of node gids observed
    std::vector<std::vector<int>> dofs_;  // list of dofs on each node that are monitored

    // patch stuff
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

    // check whether we have to reset history parameters with this bool
    bool materialhashistory_;

    enum MEASUREMENT_TYPE
    {
      dof_based,
      point_based
    };
    MEASUREMENT_TYPE meas_type_;

   public:
    /// standard constructor
    GenInvAnalysis(Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
        Teuchos::RCP<LINALG::Solver> solver,               ///< algebraic solver
        Teuchos::RCP<IO::DiscretizationWriter> output      ///< outputer
    );

    /// Read in routine for the material parameters
    void ReadInParameters();

    /// Read in routine for the dof based monitor file (default format)
    // constructor sequence moved into own method
    void ReadMonitorDofBased(int myrank);


    /// Read in routine for the point based monitor file
    void ReadMonitorPointBased(int myrank);

    /// Set material parameters
    void SetParameters(Epetra_SerialDenseVector p_cur);

    /// get current state of calculated curve
    Epetra_SerialDenseVector GetCalculatedCurve(const Epetra_Vector& disp);

    /// get Distance of target points to respective interface at matching time
    Epetra_SerialDenseVector GetDistancePointsToInterfaceContour(
        const Epetra_Vector& disp, int writestep);

    /// coarse fitting of parameters
    void PrecondParameters();

    /// fitting the material parameters using the Levenberg-Marquart algorithm
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

    /// Output of Inverse Analysis
    void PrintStorage(Epetra_SerialDenseVector delta_p);
    void PrintFile();

    // SetMaterialParameters for patch materials
    void SetMaterialParametersPatches(Epetra_SerialDenseVector& p_cur);
    // initialize patches and prolongation operator
    void InitPatches();
    // compute elementparameters from patch parameters
    void ComputeParametersFromPatches(
        Epetra_Vector p_patches, Teuchos::RCP<Epetra_Vector>& eleparameters);

    // Check if mat params are out of bounds
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
