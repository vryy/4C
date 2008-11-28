/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_strugenalpha.H

\brief Structure field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_STRUGENALPHA_H
#define ADAPTER_STRUCTURE_STRUGENALPHA_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/drt_utils.H"
#include "../drt_structure/strugenalpha.H"
#include "../drt_structure/stru_resulttest.H"



namespace ADAPTER {

  /// adapter to generalized alpha algorithm
  class StructureGenAlpha : public Structure
  {
  public:

    StructureGenAlpha(Teuchos::RCP<Teuchos::ParameterList> params,
                     Teuchos::RCP<DRT::Discretization> dis,
                     Teuchos::RCP<LINALG::Solver> solver,
                     Teuchos::RCP<IO::DiscretizationWriter> output);

    Teuchos::RCP<const Epetra_Vector> InitialGuess();
    Teuchos::RCP<const Epetra_Vector> RHS();
    Teuchos::RCP<const Epetra_Vector> Dispnp();
    Teuchos::RCP<const Epetra_Vector> Dispn();
    /// original comment: 'unknown displacements at t(n-1)'
    /// actually Dispm()-called, mid-displacement? 
    Teuchos::RCP<const Epetra_Vector> Dispnm();
    Teuchos::RCP<const Epetra_Map> DofRowMap();
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();
    Teuchos::RCP<DRT::Discretization> Discretization();
    Teuchos::RCP<const Epetra_Vector> FExtn();
    Teuchos::RCP<const Epetra_Vector> FRobin();
    //virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS() const;
    //virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS() const;

    //Teuchos::RCP<Epetra_Vector> InterfaceDispn();

    /// @name Time step helpers
    //@{

    /// Get the current time
    virtual const double GetTime() const { dserror("Not impl."); return -1.0; }
    
    /// Get upper limit of time range of interest
    virtual const double GetTimeEnd() const { dserror("Not impl."); return -1.0; }

    /// Get time step size \f$\Delta t_n\f$
    virtual const double GetTimeStepSize() const { dserror("Not impl."); return -1.0; }

    /// Get current time step \f$n\f$ (sometimes called \f$k\f$ [but not in Baci])
    virtual const int GetTimeStep() const { dserror("Not impl."); return -1; }

    /// Get number of time steps
    virtual const int GetTimeNumStep() const { dserror("Not impl."); return -1; }

    void PrepareTimeStep();
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disp);
    void Update();
    void UpdateIterIncrConstr(Teuchos::RCP<Epetra_Vector> lagrincr)
    {
      structure_.UpdateIterIncrConstr(lagrincr);
    }
    void Output();
    
    //@}

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const { return interface_; }

    virtual const Epetra_Map& DomainMap();

    virtual void ReadRestart(int step);

    virtual void Solve();
    
    /// are there any algebraic constraints?
    virtual bool HaveConstraint()
    {
      return structure_.HaveConstraint();
    };

    /// get constraint manager defined in the structure
    virtual const RCP<const UTILS::ConstrManager> GetConstraintManager()
    {
      return structure_.GetConstraintManager();
    };
    
    //! Return MapExtractor for Dirichlet boundary conditions 
    virtual const RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return LINALG::ConvertDirichletToggleVectorToMaps(structure_.Dirichlet());
    };
    
    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// linear structure solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid structural solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce);

    //! @name Extract interface values

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn();

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp();

    /// extract interface forces (necessary for robin-BC)
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDispnp();

    //@}

    //! @name Apply interface forces

    /// apply interface forces to structural solver
    /*!
      This prepares a new solve of the structural field within one time
      step. The middle values are newly created.

      \note This is not yet the most efficient implementation.
    */
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce);

    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> iforce,
                                          Teuchos::RCP<Epetra_Vector> ifluidvel);

    //@}

    /// create result test for encapulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  private:

    /// the actual structure algorithm
    StruGenAlpha structure_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor interface_;

    /// sum of displacement vectors already applied
    Teuchos::RCP<Epetra_Vector> sumdisi_;

    //! @name local copies of input parameters

    Teuchos::RCP<DRT::Discretization> dis_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<LINALG::Solver> solver_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //@}
  };

}

#endif
#endif
