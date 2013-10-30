/*----------------------------------------------------------------------*/
/*!
\file ale.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef ALE_H
#define ALE_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

/// forward declarations
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

namespace ALE {

/// forward declaration
namespace UTILS
{
  class MapExtractor;
}


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
class Ale
{
  friend class AleResultTest;

public:

  Ale(Teuchos::RCP<DRT::Discretization> actdis,
      Teuchos::RCP<LINALG::Solver> solver,
      Teuchos::RCP<Teuchos::ParameterList> params,
      Teuchos::RCP<IO::DiscretizationWriter> output,
      bool dirichletcond);

  virtual ~Ale() {}

  //! @name Apply interface displacements

  /// displacements at the interface are prescribed
  void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp);

  /// displacements at the free surface are prescribed
  void ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp);

  //@}

  //! @name Time step helpers

  /// a very simple time loop to be used for standalone ALE problems
  void Integrate();

  /// start a new time step
  void PrepareTimeStep();

  /// calculate linear system
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp, std::string incrementtype = std::string("step")) = 0;

  /// take the current solution to be the final one for this time step
  virtual void Update() = 0;

  /// write output
  virtual void Output() = 0;

  /*!
  \brief Reset time step

  In case of time step size adaptivity, time steps might have to be repeated.
  Therefore, we need to reset the solution back to the initial solution of the
  time step.

  \author: mayr.mt
  \date: 08/2013
  */
  void ResetStep();

  /*!
  \brief Reset time and step in case that a time step has to be repeated

  ALE field increments time and step at the beginning of a time step. If a time
  step has to be repeated, we need to take this into account and decrease time and
  step beforehand. They will be incremented right at the beginning of the repetition
  and, thus, everything will be fine. Currently, this is needed for time step size
  adaptivity in FSI.

  \author: mayr.mt
  \date: 08/2013
   */
  void ResetTime(const double dtold);

  /// set time step step size
  void SetDt(const double dtnew);

  /// read restart for given step
  void ReadRestart(const int step);

  //@}

  //! @name Solver calls

  /// solve this time step with currently applied interface displacements
  virtual void Solve() = 0;

  // Solve call for biofilm growth
  virtual void SolveBioGr() = 0;

  /// get the linear solver object used for this field
  Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

  /// get the linear solver object used for this field
  virtual Teuchos::RCP<LINALG::Preconditioner> ConstPreconditioner() = 0;

  //@}

  //! @name Extract displacement

  /// write access to whole displacement field at time step \f$t^{n+1}\f$
  virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const = 0;

  /// write access to whole displacement field at time step \f$t^{n}\f$
  virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispn() const = 0;

  //@}

  //! @name Vector access

  /// initial guess of Newton's method
  Teuchos::RCP<Epetra_Vector> InitialGuess() { return residual_; }

  /// rhs of Newton's method
  Teuchos::RCP<Epetra_Vector> RHS() const { return residual_; }

  //@}

  virtual void SolveAleXFluidFluidFSI() { };

  //! @name Misc

  /// dof map of vector of unknowns
  Teuchos::RCP<const Epetra_Map> DofRowMap();

  /// direct access to system matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

  /// direct access to system matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

  /// direct access to discretization
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  /// communication object at the interface
  Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const { return interface_; }

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true) = 0;

  /// setup Dirichlet boundary condition map extractor
  void SetupDBCMapEx(bool dirichletcond);

  /// create result test for encapsulated algorithm
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return dbcmaps_; }

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

  /// Reset state vectors
  void Reset();


  //@}

protected:

  //! @name Misc

  Teuchos::RCP<DRT::Discretization> discret_;
  Teuchos::RCP<LINALG::Solver> solver_;
  Teuchos::RCP<Teuchos::ParameterList> params_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! Dirichlet BCs with local co-ordinate system
  Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;

  //@}

  //! @name Algorithm core variables
  int step_;
  int numstep_;
  double time_;
  double maxtime_;
  double dt_;
  int uprestart_;
  bool incremental_;
  //@}

  //! @name matrices, vectors, map extractors

  Teuchos::RCP<LINALG::SparseOperator> sysmat_;
  Teuchos::RCP<Epetra_Vector>    residual_;

  Teuchos::RCP<Epetra_Vector>    dispnp_;
  Teuchos::RCP<Epetra_Vector>    dispn_;

  //! maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  //! the interface map setup for interface <-> full translation
  Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;

  //@}

};


/// Base class of algorithms that use an ale field
class AleBaseAlgorithm
{
public:

  /// constructor
  explicit AleBaseAlgorithm(const Teuchos::ParameterList& prbdyn, int disnum = 0);

  /// virtual destructor to support polymorph destruction
  virtual ~AleBaseAlgorithm();

  /// ale field solver
  Ale& AleField() { return *ale_; }
  const Ale& AleField() const { return *ale_; }
  Teuchos::RCP<Ale> AleFieldrcp(){ return ale_; }

private:

  /// setup ale algorithm (overriding some parameters with
  /// values specified in given problem-dependent ParameterList)
  void SetupAle(const Teuchos::ParameterList& prbdyn, int disnum = 0);

  /// ale field solver
  Teuchos::RCP<Ale> ale_;
};

}

#endif
