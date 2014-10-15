/*----------------------------------------------------------------------------*/
/*!
\file ale.H

\brief

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef ALENEW_H_
#define ALENEW_H_

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_adapter/ad_ale.H"

#include "../drt_inpar/inpar_ale.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace LINALG
{
  class Solver;
  class Preconditioner;
  class SparseOperator;
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
}

namespace DRT
{
  class ResultTest;
  class Discretization;

  namespace UTILS
  {
    class LocsysManager;
  }
}

namespace IO
{
  class DiscretizationWriter;
}




namespace ALENEW {

namespace UTILS
{
  class MapExtractor;
}

/// forward declaration

/// General ALE interface
/*!
  An abstract base class meant to be derived for each ALE
  implementation we want to support.

  The purpose of the ALE field is to take an interface displacement and to
  return the deformed interface. All in all pretty simple.

  \warning Further cleanup is still needed.

  \sa Structure, Fluid
  \author u.kue
  \date 11/07
 */
  class Ale : public ::ADAPTER::Ale
  {
    //friend class AleResultTest;

  public:

  Ale(Teuchos::RCP<DRT::Discretization> actdis,
      Teuchos::RCP<LINALG::Solver> solver,
      Teuchos::RCP<Teuchos::ParameterList> params,
      Teuchos::RCP<IO::DiscretizationWriter> output,
      bool dirichletcond);

  virtual ~Ale() {}

  /*! \brief Create Systemmatrix
   *
   * We allocate the LINALG object just once, the result is an empty LINALG object.
   * Evaluate has to be called separately.
   *
   */
  void CreateSystemMatrix(bool full = true ///< Sparsematrix (true) or Blocksparsematrix (false)
      );

  /*! \brief evaluate and assemble residual #residual_ and jacobian matrix #sysmat_
   *
   *  use this as evaluate routine for pure ALE problems as well as for coupled problems.
   *  Update in case of monolithic coupling is done by passing stepinc, Teuchos::null is assumed for
   *  non monolithic case.
   */
   virtual void Evaluate(
       Teuchos::RCP<const Epetra_Vector> stepinc = Teuchos::null ///< step increment such that \f$ x_{n+1}^{k+1} = x_{n}^{converged}+ stepinc \f$
       );

   /// solve this time step with currently applied interface displacements
   virtual int Solve();

   /// get the linear solver object used for this field
   Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

   //! update displacement with iterative increment
   virtual void UpdateIter(); // TODO (Biehler) make private later

   /// take the current solution to be the final one for this time step
   virtual void Update();

   /// convergence test for newton
   virtual const bool Converged(const int iter);

   virtual void EvaluateElements();

   const std::string ElementActionString(const enum INPAR::ALE::AleDynamic name ///< enum to convert
       );

   Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner();


   virtual void CreatePreconditioner(bool full = true ///< Sparsematrix (true) or Blocksparsematrix (false)
       );
   //__________________________________________________________________________________________________CLEANUP

  //! @name Apply interface displacements

  /// displacements at the interface are prescribed
  virtual void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp);

  /// displacements at the free surface are prescribed
  virtual void ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp);

  /// displacements at the ale update condition nodes are prescribed
  virtual void ApplyAleUpdateDisplacements(Teuchos::RCP<Epetra_Vector> audisp);

  //@}

  //! @name Time step helpers

  /// a very simple time loop to be used for standalone ALE problems
  virtual int Integrate();

  /// start a new time step
  virtual void PrepareTimeStep();

  virtual void TimeStep();



  /// write output
  virtual void Output();

  /*!
  \brief Reset time step

  In case of time step size adaptivity, time steps might have to be repeated.
  Therefore, we need to reset the solution back to the initial solution of the
  time step.

  \author mayr.mt
  \date 08/2013
  */
  virtual void ResetStep();

  /*!
  \brief Reset time and step in case that a time step has to be repeated

  ALE field increments time and step at the beginning of a time step. If a time
  step has to be repeated, we need to take this into account and decrease time and
  step beforehand. They will be incremented right at the beginning of the repetition
  and, thus, everything will be fine. Currently, this is needed for time step size
  adaptivity in FSI.

  \author mayr.mt
  \date 08/2013
   */
  virtual void ResetTime(const double dtold);

  /// Get current simulation time
  virtual double Time() const { return time_; }

  /// Get the time step size
  virtual double Dt() const { return dt_; }

  /// set time step step size
  virtual void SetDt(const double dtnew);

  /// read restart for given step
  virtual void ReadRestart(const int step);

  //@}

  //! @name Reading access to displacement
  //@{

  /// get the whole displacement field at time step \f$t^{n+1}\f$
  virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const {return dispnp_;}

  /// get the whole displacement field at time step \f$t^{n}\f$
  virtual Teuchos::RCP<const Epetra_Vector> Dispn() const {return dispn_;}

  //@}

  //! @name Writing access to displacement

  /// write access to whole displacement field at time step \f$t^{n+1}\f$
  virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const {return dispnp_;}

  //@}

  //! @name Vector access

  /// initial guess of Newton's method
  Teuchos::RCP<const Epetra_Vector> InitialGuess() const { return zeros_; }

  /// rhs of Newton's method
  Teuchos::RCP<const Epetra_Vector> RHS() const { return rhs_; }

  //@}

  //! @name Misc

  /// dof map of vector of unknowns
  Teuchos::RCP<const Epetra_Map> DofRowMap() const;

  /// direct access to system matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

  /// direct access to system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

  /// direct access to discretization
  Teuchos::RCP<const DRT::Discretization> Discretization() const { return discret_; }

  /// writing access to discretization
  Teuchos::RCP<DRT::Discretization> WriteAccessDiscretization() { return discret_; }

  /// communication object at the interface
  Teuchos::RCP<const ALENEW::UTILS::MapExtractor> Interface() const { return interface_; }

  /// setup Dirichlet boundary condition map extractor
  void SetupDBCMapEx(bool dirichletcond);

  /// create result test for encapsulated algorithm
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return dbcmaps_; }

  //! Return (rotatory) transformation matrix of local co-ordinate systems
  Teuchos::RCP<const LINALG::SparseMatrix> GetLocSysTrafo() const;

  //! Return locsys manager
  Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() {return locsysman_;}

  //! Apply Dirichlet boundary conditions on provided state vectors
  void ApplyDirichletBC
  (
    Teuchos::ParameterList& params,
    Teuchos::RCP<Epetra_Vector> systemvector,   //!< (may be Teuchos::null)
    Teuchos::RCP<Epetra_Vector> systemvectord,  //!< (may be Teuchos::null)
    Teuchos::RCP<Epetra_Vector> systemvectordd, //!< (may be Teuchos::null)
    bool recreatemap  //!< recreate mapextractor/toggle-vector
                      //!< which stores the DOF IDs subjected
                      //!< to Dirichlet BCs
                      //!< This needs to be true if the bounded DOFs
                      //!< have been changed.
  );

  /// Reset state vectors to zero
  void Reset();


  //@}

  private:

  //! @name Misc

  Teuchos::RCP<DRT::Discretization> discret_;
  Teuchos::RCP<LINALG::Solver> solver_;
  Teuchos::RCP<Teuchos::ParameterList> params_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! Dirichlet BCs with local co-ordinate system
  Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

  //@}

  //! @name Algorithm core variables
  int step_;              ///< step counter
  int numstep_;           ///< max number of steps
  double time_;           ///< simulation time
  double maxtime_;        ///< max simulation time
  double dt_;             ///< time step size
  int writerestartevery_; ///< write restart every n steps
  int writeresultsevery_; ///< write results every n steps
  bool incremental_;      ///< incremental ale yes/no?
  //@}

  //! @name matrices, vectors
  //@{

  Teuchos::RCP<LINALG::SparseOperator> sysmat_; ///< stiffness matrix

  /*! \brief residual vector
   *
   *  This is the "mechanical" residual \f$res = - f_{int}\f$ as it comes
   *  from the discret_->Evaluate() call.
   *
   *  \author mayr.mt \date 10/2014
   */
  Teuchos::RCP<Epetra_Vector> residual_;

  /*! \brief right hand side of Newton-type algorithm
   *
   *  Use this as the right hand side for a Newton algorithm. It should equal
   *  the negative residual: #rhs_ = - #residual_
   *
   *  We update this variabel only right after the discret_->Evaluate() call.
   *
   *  \warning DO NOT TOUCH THIS VARIBALE AT OHTER PLACES!!!
   *
   *  \author mayr.mt \date 10/2014
   */
  Teuchos::RCP<Epetra_Vector> rhs_;

  Teuchos::RCP<Epetra_Vector>    dispnp_; ///< unknown solution at \f$t_{n+1}\f$
  Teuchos::RCP<Epetra_Vector>    dispn_;  ///< known solution at \f$t_{n}\f$
  Teuchos::RCP<Epetra_Vector>    disi_;  ///< iterative displacement increment
  Teuchos::RCP<const Epetra_Vector>    zeros_;  ///< zero vector for dbc handling
  //@}

  //! @name map extractors
  //@{

  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_; ///< map extractor for extracting Dirichlet and free DOF sets

  Teuchos::RCP<ALENEW::UTILS::MapExtractor> interface_; ///< interface map setup for interface <-> full translation

  //@}

  /// print info about current time step to screen
  virtual void PrintTimeStepHeader() const;

  /// write restart data
  virtual void OutputRestart(bool& datawritten);

  /// write output data
  virtual void OutputState(bool& datawritten);

  /// cached preconditioner on constant system matrix
  Teuchos::RCP<LINALG::Preconditioner> precond_;

  /// ale formulation read from inputfile
  const INPAR::ALE::AleDynamic aletype_;

  //! @name solver parameters
  //@{
  //! maximum number of newton iterations
  const int maxiter_;

  //! tolerance of length scaled L2 residual norm
  const double tolres_;

  //! tolerance of length scaled L2 increment norm
  const double toldisp_;

  //! error handling in case of unconverged nonlinear solver
  const INPAR::ALE::DivContAct divercont_;

  //@}

};

}

#endif /* ALENEW_H_ */
