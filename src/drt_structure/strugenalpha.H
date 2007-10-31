/*!----------------------------------------------------------------------
\file strugenalpha.H

\class StruGenAlpha

\brief Generalized Alpha time integration for structural problems

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef STRUGENALPHA_H
#define STRUGENALPHA_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../io/io_drt.H"
#include "../drt_contact/drt_contact_manager.H"
#include "../drt_lib/linalg_mfoperator.H"

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

using namespace std;
using namespace Teuchos;
class NoxInterface;

/*!
\class StruGenAlpha

\brief Generalized Alpha time integration for structural problems

\author gee (gee@lnm.mw.tum.de)
*/
class StruGenAlpha
{
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
  void Integrate();

  /*!
  \brief Integrate single step in time

  */
  void IntegrateStep();

  /*!
  \brief do constant predictor step

  */
  void ConstantPredictor();

  /*!
  \brief do matrix free constant predictor

  */
  void MatrixFreeConstantPredictor();

  /*!
  \brief do consistent constant predictor step

  */
  void ConsistentPredictor();

  /*!
  \brief do full Newton iteration

  */
  void FullNewton();

  /*!
  \brief do modified Newton iteration (makes sense with direct solvers only)

  */
  void ModifiedNewton();

  /*!
  \brief do matrix free Newton iteration (only with Aztec)

  */
  void MatrixFreeNewton();

  /*!
  \brief do nonlinear cg iteration using nox

  */
  void NonlinearCG();

  /*!
  \brief do pseudo transient continuation nonlinear iteration

  Pseudo transient continuation is a variant of full newton that has a
  larger convergence radius than newton and is therefore more stable and/or
  can do larger time steps

  */
  void PTC();

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
  \brief update configuration and output to file/screen

  */
  void UpdateandOutput();

  /*!
  \brief restart calculation from a given step

  */
  void ReadRestart(int step);

  /*!
  \brief return displacements at midpoint

  */
  const Epetra_Vector& Getdu() { return *disi_; }

  /*!
  \brief return rowmap

  */
  const Epetra_Map* GetMap() { return discret_.DofRowMap(); }

  /*!
  \brief Get time integration scheme (TIS) parameters
  
  */
  void GetTISPara(double& beta, 
                  double& gamma, 
                  double& alpham, 
                  double& alphaf);

  /*
  \brief set time step size

  */
  void SetTimeStepSize(const double& timstpsiz);

  /*
  \brief Get last converged accelerations

  */
  const RefCountPtr<Epetra_Vector>& GetAcc() 
  { return acc_; }

  /*
  \brief Get last converged accelerations

  */
  const RefCountPtr<Epetra_Vector>& GetAccn() 
  { return accn_; }

protected:

  // don't want = operator and cctor
  StruGenAlpha operator = (const StruGenAlpha& old);
  StruGenAlpha(const StruGenAlpha& old);

  ParameterList&                params_;
  DRT::Discretization&          discret_;
  LINALG::Solver&               solver_;
  IO::DiscretizationWriter&     output_;
  int                           myrank_;
  int                           maxentriesperrow_;
  double                        norm_;

  bool                          havecontact_;
  RCP<CONTACT::Manager>         contactmanager_;

  RefCountPtr<Epetra_CrsMatrix> stiff_;
  RefCountPtr<Epetra_CrsMatrix> mass_;
  RefCountPtr<Epetra_CrsMatrix> damp_;

  RefCountPtr<Epetra_Vector>    dirichtoggle_;
  RefCountPtr<Epetra_Vector>    invtoggle_;
  RefCountPtr<Epetra_Vector>    zeros_;
  RefCountPtr<Epetra_Vector>    dis_;
  RefCountPtr<Epetra_Vector>    vel_;
  RefCountPtr<Epetra_Vector>    acc_;
  RefCountPtr<Epetra_Vector>    disn_;
  RefCountPtr<Epetra_Vector>    veln_;
  RefCountPtr<Epetra_Vector>    accn_;
  RefCountPtr<Epetra_Vector>    dism_;
  RefCountPtr<Epetra_Vector>    velm_;
  RefCountPtr<Epetra_Vector>    accm_;
  RefCountPtr<Epetra_Vector>    disi_;
  RefCountPtr<Epetra_Vector>    fint_;
  RefCountPtr<Epetra_Vector>    fext_;
  RefCountPtr<Epetra_Vector>    fextm_;
  RefCountPtr<Epetra_Vector>    fextn_;
  RefCountPtr<Epetra_Vector>    fresm_;

  // stuff needed to use NOX as nonlinear solver
  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> fineinterface_;
  RefCountPtr<NLNML::NLNML_Preconditioner>        prec_;
  RefCountPtr<NOX::Epetra::MatrixFree>            matfreeoperator_;
  RefCountPtr<NOX::Epetra::LinearSystemAztecOO>   rcpazlinsys_;
  RefCountPtr<NOX::StatusTest::Combo>             combo_;

}; // class StruGenAlpha



/*!
\class NoxInterface

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
    return integrator_.stiff_.get();
  }

  /*!
  \brief Return pointer to graph of problem (graph of J)

  */
  virtual const Epetra_CrsGraph* getGraph()
  {
    if (integrator_.stiff_==null) dserror("integrator_.stiff_ = null!");
    return &(integrator_.stiff_->Graph());
  }

  /*!
  \brief Return pointer to graph of problem (graph of J)

  */
  virtual const Epetra_CrsGraph* getModifiedGraph()
  {
    if (integrator_.stiff_==null) dserror("integrator_.stiff_ = null!");
    return &(integrator_.stiff_->Graph());
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
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
