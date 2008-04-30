/*!----------------------------------------------------------------------
\file strugenalpha.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUGENALPHA_H
#define STRUGENALPHA_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_io/io.H"
#include "../drt_contact/drt_contact_manager.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_surfstress/drt_surfstress_manager.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include <NOX.H>
#include <NOX_Epetra.H>
#include "ml_common.h"
#include "ml_include.h"
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"

// definitions in algorithm
#define STRUGENALPHA_INCRUPDT  /* incremental update of vel and acc */
//#define STRUGENALPHA_STRONGDBC  /* strong Dirichlet BCs, i.e. DBC reigns  */
                                /* its displacements, velocities and acceleration */
#ifdef STRUGENALPHA_STRONGDBC
#ifdef STRUGENALPHA_INCRUPDT
#undef STRUGENALPHA_INCRUPDT
#endif
#endif
//#define STRUGENALPHA_BE  /* enable backward Euler with alpham=0,alphaf=0,beta=1/2,delta=1,gamma=1 */
//#define STRUGENALPHA_FINTLIKETR  /* enable TR-like mid-internal force vector F_{int}*/
                                 /* TR-like mid-point internal force:
                                  *    F_{int,m} = (1-alpha_f)*F_{int,n+1} + alpha_f*F_{int,n}
                                  *              = (1-alpha_f)*F_{int}(D_{n+1}) + alpha_f*F_{int}(D_n)
                                  * in contrast to the IMR-like mid-point internal force (default)
                                  *    F_{int,m} = F_{int}(D_m)
                                  *              = F_{int}( (1-alpha_f)*D_{n+1} + alpha_f*D_n ) */

using namespace std;
using namespace Teuchos;
class NoxInterface;

/*!
\brief Generalized Alpha time integration for structural problems

\author gee (gee@lnm.mw.tum.de)
*/
class StruGenAlpha
{
  friend class StruResultTest;
public:

  /*!
  \brief The nox interface is a friend of this class

  */
  friend class NoxInterface;


  /*!
  \brief Standard Constructor

  */
  explicit StruGenAlpha(ParameterList& params,
                        DRT::Discretization& dis,
                        LINALG::Solver& solver,
                        IO::DiscretizationWriter& output);

  /*!
  \brief Destructor

  */
  virtual ~StruGenAlpha() {};

  /*!
  \brief Create default parameter list

  */
  static void SetDefaults(ParameterList& params);



  /*!
  \brief Do time integration

  */
  virtual void Integrate();

  /*!
  \brief Integrate single step in time

  */
  virtual void IntegrateStep();

  /*!
  \brief do constant predictor step

  */
  virtual void ConstantPredictor();

  /*!
  \brief do consistent constant predictor step

  */
  virtual void ConsistentPredictor();

  /*!
    \brief build linear system matrix and rhs

    Monolithic FSI needs to access the linear fluid problem.

    \param disp new guess at displacement
  */
  void Evaluate(Teuchos::RCP<const Epetra_Vector> disp);

  /*!
  \brief do full Newton iteration

  */
  virtual void FullNewton();



  /*!
    \brief do classical augmented lagrange for volume constraint

    */
  virtual void NonLinearUzawaFullNewton(int predictor);


  /*!
    \brief do full Newton iteration respecting volume constraint

    */
  virtual void FullNewtonLinearUzawa();

  /*!
  \brief do modified Newton iteration (makes sense with direct solvers only)

  */
  virtual void ModifiedNewton();

  /*!
  \brief do nonlinear cg iteration using nox

  */
  virtual void NonlinearCG();

  /*!
  \brief do pseudo transient continuation nonlinear iteration

  Pseudo transient continuation is a variant of full newton that has a
  larger convergence radius than newton and is therefore more stable and/or
  can do larger time steps

  */
  virtual void PTC();

  /*!
  \brief Compute the residual to a given state x

  This method is used by the nox user interface class and by the linalg matrix free operator

  \param x (in/out) : Current solution state (dirichlet dofs will be blanked)
  \param F (out): residual matching current solution state (dirichlet dofs will be blanked)
  */
  virtual void computeF(const Epetra_Vector& x, Epetra_Vector& F);

  /*!
  \brief Compute the Jacobian to a given state x

  This method is used by the nox user interface

  \param x (in/out) : Current solution state (dirichlet dofs will be blanked)
  */
  virtual void computeJacobian(const Epetra_Vector& x);

  /*!
  \brief determine new state at t_{n+1}

  */
  virtual void ExtrapolateEndState();

  /*!
  \brief update configuration and output to file/screen

  */
  virtual void UpdateandOutput();

  /// update configuration
  virtual void Update();

  /// output to file/screen
  virtual void Output();

  /*!
  \brief restart calculation from a given step

  */
  virtual void ReadRestart(int step);

  /*!
  \brief return displacement increments at midpoint

  */
  virtual const Epetra_Vector& Getdu() { return *disi_; }

  /*!
  \brief return displacements at midpoint

  */
  virtual const Epetra_Vector& Getdm() { return *dism_; }

  /*!
  \brief return residual vector at midpoint

  */
  virtual const Epetra_Vector& Getfresm() { return *fresm_; }

  /*!
  \brief return rowmap

  */
  virtual const Epetra_Map* GetMap() { return discret_.DofRowMap(); }

  /*!
  \brief Get time integration scheme (TIS) parameters

  */
  virtual void GetTISPara(double& beta,double& gamma,double& alpham,double& alphaf);

  /*!
  \brief set time step size

  */
  virtual void SetTimeStepSize(const double& timstpsiz);

  /*!
  \brief Set time

  */
  virtual void SetTime(const double& tim);

  /*!
  \brief Set time step

  */
  virtual void SetTimeStep(const int& stp);

  /*!
  \brief Check convergence of Newton iteration

  */
  virtual bool Converged(const string type, const double disinorm,
                         const double resnorm, const double toldisp,
                         const double tolres);

  /*!
  \brief Check convergence of Newton iteration
  */
  virtual bool Converged(const string type, const double disinorm,
                         const double resnorm, const double volnorm,
                         const double toldisp, const double tolres,
                         const double tolvol);

  /*!
  \brief Determine reference norm of displacements and maximum of
  norms of internal, external and inertial forces

  */
  virtual void CalcRefNorms();

  /*!
  \brief Print Newton output to screen and/or error file

  */
  virtual void PrintNewton(bool printscreen, bool printerr, bool print_unconv,
                           FILE* errfile, Epetra_Time timer, int numiter,
                           int maxiter, double fresmnorm, double disinorm,
                           string convcheck);

  /*!
    \brief Print Newton output to screen and/or error file considering volume constraints

    */
    virtual void PrintNewton(bool printscreen, bool printerr, bool print_unconv,
                             FILE* errfile, Epetra_Time timer, int numiter,
                             int maxiter, double fresmnorm, double disinorm,
                             string convcheck, double volnorm, double UPara);

  /*!
  \brief Print PTC output to screen and/or error file

  \param dti (in): Psi tc stabilization parameter
  */
  virtual void PrintPTC(bool printscreen, bool printerr, bool print_unconv,
                        FILE* errfile, Epetra_Time timer, int numiter,
                        int maxiter, double fresmnorm, double disinorm,
                        string convcheck,double dti);

   /*!
  \brief Print to screen and/or error file

  */
  virtual void PrintPredictor(string convcheck, double fresmnorm);

  /*!
  \brief Get last converged accelerations

  */
  virtual const RCP<Epetra_Vector>& GetAcc()
  { return acc_; }

  /*!
  \brief Get last converged accelerations

  */
  virtual const RCP<Epetra_Vector>& GetAccn()
  { return accn_; }

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear structure problem

  Teuchos::RCP<Epetra_Vector> Residual() { return fresm_; }
  Teuchos::RCP<Epetra_Vector> Dispn() { return disn_; }
  Teuchos::RCP<Epetra_Vector> Dispm() { return dism_; }
  Teuchos::RCP<Epetra_Vector> Disp() { return dis_; }
  Teuchos::RCP<Epetra_Vector> Vel() { return vel_; }
  Teuchos::RCP<Epetra_Vector> Acc() { return acc_; }
  Teuchos::RCP<Epetra_Vector> FRobin() { return frobin_; }
  Teuchos::RCP<Epetra_Vector> FExtn() { return fextn_; }
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return stiff_; }
  Teuchos::RCP<Epetra_Vector> Dirichlet() { return dirichtoggle_; }
  Teuchos::RCP<DRT::Discretization> Discretization() { return Teuchos::rcp(&discret_,false); }
  double DispIncrFactor() { return 1./(1.-params_.get<double>("alpha f",0.459)); }
  double AlphaF() { return params_.get<double>("alpha f",0.459); }
  const Epetra_Map& DomainMap() { return mass_->DomainMap(); }
  bool HaveConstraint() {return constrMan_->HaveConstraint(); }

  /// linear solve with prescribed loads and without history
  /*!
    This is the linear solve as needed for steepest descent FSI.
   */
  Teuchos::RCP<Epetra_Vector> LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  /// setup equilibrium with additional external forces
  void ApplyExternalForce(const LINALG::MapExtractor& extractor,
                          Teuchos::RCP<Epetra_Vector> iforce);

  /// introduce (robin) fsi surface extractor object
  /*!
    Introduce the FSI surface to the fluid algorithm. This way special case
    algorithms (robin coupling) can work on the surface.
   */
  void SetFSISurface(const LINALG::MapExtractor* fsisurface) { fsisurface_ = fsisurface; }

  //@}

protected:

  // don't want = operator and cctor
  StruGenAlpha operator = (const StruGenAlpha& old);
  StruGenAlpha(const StruGenAlpha& old);

  //! @name (very) general purpose algorithm parameters

  ParameterList&                params_;
  DRT::Discretization&          discret_;
  LINALG::Solver&               solver_;
  IO::DiscretizationWriter&     output_;
  int                           myrank_;

  //@}

  //! @name The usual structural dynamics matrix triple

  RCP<LINALG::SparseMatrix> stiff_;
  RCP<LINALG::SparseMatrix> mass_;
  RCP<LINALG::SparseMatrix> damp_;

  //@}

  //! @name Vectors for the people

  RCP<Epetra_Vector>    dirichtoggle_;
  RCP<Epetra_Vector>    invtoggle_;
  RCP<Epetra_Vector>    zeros_;
  RCP<Epetra_Vector>    dis_;
  RCP<Epetra_Vector>    vel_;
  RCP<Epetra_Vector>    acc_;
  RCP<Epetra_Vector>    disn_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    accn_;
  RCP<Epetra_Vector>    dism_;
  RCP<Epetra_Vector>    velm_;
  RCP<Epetra_Vector>    accm_;
  RCP<Epetra_Vector>    disi_;
  RCP<Epetra_Vector>    fint_;
#ifdef STRUGENALPHA_FINTLIKETR
  RCP<Epetra_Vector>    fintn_;
#endif
  RCP<Epetra_Vector>    finert_;
  RCP<Epetra_Vector>    fvisc_;
  RCP<Epetra_Vector>    fext_;
  RCP<Epetra_Vector>    fextm_;
  RCP<Epetra_Vector>    fextn_;
  RCP<Epetra_Vector>    fresm_;
  RCP<Epetra_Vector>    frobin_;

  //@}

  double                        ref_disnorm_;
  double                        ref_fnorm_;

  double                        V0_;

  /// whatever constraints
  RCP<ConstrManager>   constrMan_;

  /// surface stresses
  RCP<DRT::SurfStressManager>   surf_stress_man_;

  // For a relative convergence check, reference norms have to be
  // calculated. We need to know whether we are in the first step
  // of a simulation (either step=0 or step=genprob.restart) in
  // order to determine these norms at the right time.
  bool firststep_;

  //! @name stuff needed to use NOX as nonlinear solver
  RCP<NLNML::NLNML_FineLevelNoxInterface> fineinterface_;
  RCP<NLNML::NLNML_Preconditioner>        prec_;
  RCP<NOX::Epetra::MatrixFree>            matfreeoperator_;
  RCP<NOX::Epetra::LinearSystemAztecOO>   rcpazlinsys_;
  RCP<NOX::StatusTest::Combo>             combo_;
  //@}

  /// (robin) fsi coupling row dof map extractor
  const LINALG::MapExtractor* fsisurface_;

}; // class StruGenAlpha



/*!
\brief An interface to Nox/ML for structural gen alpha

\author gee (gee@lnm.mw.tum.de)
*/
class NoxInterface : public NLNML::NLNML_FineLevelNoxInterface
{
public:


  /*!
  \brief Standard Constructor

  */
  NoxInterface(StruGenAlpha& integrator, int printlevel) :
    NLNML_FineLevelNoxInterface(),
    integrator_(integrator)
  {
    isnewJacobian_   = false;
    numJacobian_     = 0;
    t_               = 0.;
    ncalls_computeF_ = 0;
    printlevel_      = printlevel;
    return;
  }


  /*!
  \brief Destructor

  */
  virtual ~NoxInterface()
  {
    return;
  }



  //! @name calculation methods

  /*!
  \brief For some given solution vector x, compute residual forces F

  */
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F,
			const FillType fillFlag)
  {
    isnewJacobian_ = false;
    ncalls_computeF_++;
    Epetra_Time timer(integrator_.discret_.Comm());
    double t0 = timer.ElapsedTime();
    //========================================================================
    integrator_.computeF(x,F);
    //========================================================================
    double t1 = timer.ElapsedTime();
    t_ += (t1-t0);
    return true;
  }



  /*!
  \brief For some given solution state x, compute tangent matrix J

  */
  virtual bool computeJacobian(const Epetra_Vector& x,
                               Epetra_Operator& Jac)
  {
    isnewJacobian_ = true;
    numJacobian_++;
    ncalls_computeF_++;
    Epetra_Time timer(integrator_.discret_.Comm());
    double t0 = timer.ElapsedTime();
    //========================================================================
    integrator_.computeJacobian(x);
    //========================================================================
    double t1 = timer.ElapsedTime();
    t_ += (t1-t0);
    return true;
  }

  /*!
  \brief Return pointer to last tangent matrix J

  */
  virtual Epetra_CrsMatrix* getJacobian()
  {
    if (integrator_.stiff_==null) dserror("integrator_.stiff_ = null!");
    return &*integrator_.stiff_->EpetraMatrix();
  }

  /*!
  \brief Return pointer to graph of problem (graph of J)

  */
  virtual const Epetra_CrsGraph* getGraph()
  {
    if (integrator_.stiff_==null) dserror("integrator_.stiff_ = null!");
    return &(integrator_.stiff_->EpetraMatrix()->Graph());
  }

  /*!
  \brief Return pointer to graph of problem (graph of J)

  */
  virtual const Epetra_CrsGraph* getModifiedGraph()
  {
    if (integrator_.stiff_==null) dserror("integrator_.stiff_ = null!");
    return &(integrator_.stiff_->EpetraMatrix()->Graph());
  }

  /*!
  \brief Return pointer to current solution vector

  */
  virtual const Epetra_Vector* getSolution()
  {
    return integrator_.disi_.get();
  }

  /*!
  \brief Return dofrowmap of problem

  */
  virtual const Epetra_Map& getMap()
  {
    return integrator_.mass_->RowMap();
  }

  /*!
  \brief Return nullspace of tangent matrix J.

         The nullspace was
         previously computed for ML in the linear solver section
         and is extracted from there. A copy is made as the
         nonlinear ml precondition will take ownership of it.

  */
  virtual double* Get_Nullspace(const int nummyrows, const int numpde,
                                const int dim_nullsp)
  {
    ParameterList& p = integrator_.solver_.Params();
    bool hasml = p.isSublist("ML Parameters");
    if (!hasml) dserror("Solver does not have ML configured from input file");
    ParameterList&  mllist = p.sublist("ML Parameters");
    bool hasns = mllist.isParameter("null space: vectors");
    if (!hasns) dserror("Ml ParameterList does not have a nullspace inside");
    double* ns = mllist.get<double*>("null space: vectors",NULL);
    if (!ns) dserror("Nullspace vector is NULL");

    double* nullspace = new double[nummyrows*dim_nullsp];
    for (int i=0; i<nummyrows*dim_nullsp; ++i) nullspace[i] = ns[i];
    return nullspace;
  }


  /*!
  \brief For problems with non-constant nodal block size one can use the
         VBMETIS coarsening scheme (currently not supported).

  */
  virtual bool getBlockInfo(int *nblocks, vector<int>& blocks, vector<int>& block_pde)
  {
    cout << "getBlockInfo not impl.\n"; fflush(stdout);
    exit(0);
    return false;
  }

  //@}

private:

  // don't want = operator and cctor
  NoxInterface operator = (const NoxInterface& old);
  NoxInterface(const NoxInterface& old);

  StruGenAlpha& integrator_;


}; // class NoxInterface


#endif  // #ifndef STRUGENALPHA_H
#endif  // #ifdef CCADISCRET
