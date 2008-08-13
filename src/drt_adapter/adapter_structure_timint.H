/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_timint.H

\brief Structural time integration field adapter

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_TIMINT_H
#define ADAPTER_STRUCTURE_TIMINT_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_lib/drt_utils.H"
#include "../drt_structure/strugenalpha.H"
#include "../drt_structure/stru_resulttest.H"

#include "adapter_structure.H"



namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to general structural time integration
   *
   * \date 08/08
   */
  class StructureTimInt : public Structure
  {
  public:

    StructureTimInt(Teuchos::RCP<Teuchos::ParameterList> params,
                    Teuchos::RCP<DRT::Discretization> dis,
                    Teuchos::RCP<LINALG::Solver> solver,
                    Teuchos::RCP<IO::DiscretizationWriter> output);

    Teuchos::RCP<const Epetra_Vector> InitialGuess();
    Teuchos::RCP<const Epetra_Vector> RHS();
    Teuchos::RCP<const Epetra_Vector> Dispnp();
    Teuchos::RCP<const Epetra_Vector> Dispn();
    Teuchos::RCP<const Epetra_Vector> Dispnm();
    Teuchos::RCP<const Epetra_Map> DofRowMap();
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();
    Teuchos::RCP<DRT::Discretization> Discretization();
    double DispIncrFactor();
    Teuchos::RCP<const Epetra_Vector> FExtn();
    Teuchos::RCP<const Epetra_Vector> FRobin();
    //virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS() const;
    //virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS() const;

    //Teuchos::RCP<Epetra_Vector> InterfaceDispn();

    void PrepareTimeStep();
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disp);
    void Update();
    void Output();

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const { return interface_; }

    virtual const Epetra_Map& DomainMap();

    virtual void ReadRestart(int step);

    virtual void Solve();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// linear structure solve with just a interface load
    ///
    ///  The very special solve done in steepest descent relaxation
    ///  calculation (and matrix free Newton Krylov).
    ///
    ///  \note Can only be called after a valid structural solve.
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
    ///
    ///  This prepares a new solve of the structural field within one time
    ///  step. The middle values are newly created.
    ///
    ///  \note This is not yet the most efficient implementation.
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
  };  // class StructureTimInt

  /*====================================================================*/
  /*!
   * \brief Structure field solver
   * 
   * \date 08/08
   */
  class StructureTimIntBaseAlgorithm
  {
  public:

    /// constructor
    explicit StructureTimIntBaseAlgorithm(const Teuchos::ParameterList& prbdyn);

    /// virtual destructor to support polymorph destruction
    virtual ~StructureTimIntBaseAlgorithm();

    /// structural field solver
    Structure& StructureField() { return *structure_; }
    const Structure& StructureField() const { return *structure_; }

  private:

    /// setup structure algorithm
    void SetupStructure(const Teuchos::ParameterList& prbdyn);

    /// structural field solver
    Teuchos::RCP<Structure> structure_;
  };  // class StructureTimIntBaseAlgorithm

}  // namespace ADAPTER

#endif
#endif
