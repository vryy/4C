/*----------------------------------------------------------------------*/
/*!
\file fsi_structure.H

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

#ifndef MFSI_STRUCTURE_H
#define MFSI_STRUCTURE_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_structure/strugenalpha.H"

#include "fsi_utils.H"


namespace FSI {

  /// general structural field interface
  /*!

  The point is to keep FSI as far apart from our field solvers as
  possible. Each structure field solver we want to use should get its own
  subclass of this. The FSI algorithm should be able to extract all the
  information from the structure field it needs using this interface.

  \warning A cleanup is needed badly.

  \author u.kue
  \date 11/07
  */
  class Structure
  {
  public:

    virtual ~Structure();

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() = 0;
    virtual Teuchos::RCP<Epetra_Vector> RHS() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Dispnp() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Disp() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Dispm() = 0;
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;
    virtual Teuchos::RCP<Epetra_Operator> SysMat() = 0;
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;
    virtual double DispIncrFactor() = 0;

    virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS() = 0;
    virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS() = 0;

    virtual void PrepareTimeStep() = 0;
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disp) = 0;
    virtual void Update() = 0;
    virtual void Output() = 0;

    virtual Teuchos::RCP<Epetra_Map> InterfaceMap() = 0;
    virtual Teuchos::RCP<Epetra_Map> InnerDisplacementRowMap() = 0;

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const = 0;

    virtual const Epetra_Map& DomainMap() = 0;

    virtual void ReadRestart(int step) = 0;

    virtual void UpdateandOutput() = 0;

    virtual void Solve() = 0;

    /// linear structure solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid structural solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    //! @name Extract interface values

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn() = 0;

    /// extract displacements
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDisplacement() = 0;

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDisplacement() = 0;

    //@}

    //! @name Apply interface forces

    /// apply interface forces to structural solver
    /*!
      This prepares a new solve of the structural field within one time
      step. The middle values are newly created.

      \note This is not yet the most efficient implementation.
    */
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce) = 0;

    //@}
  };


  /// adapter to generalized alpha algorithm
  class StructureAdapter : public Structure
  {
  public:

    StructureAdapter(Teuchos::RCP<Teuchos::ParameterList> params,
                     Teuchos::RCP<DRT::Discretization> dis,
                     Teuchos::RCP<LINALG::Solver> solver,
                     Teuchos::RCP<IO::DiscretizationWriter> output);

    Teuchos::RCP<const Epetra_Vector> InitialGuess();
    Teuchos::RCP<Epetra_Vector> RHS();
    Teuchos::RCP<Epetra_Vector> Dispnp();
    Teuchos::RCP<Epetra_Vector> Disp();
    Teuchos::RCP<Epetra_Vector> Dispm();
    Teuchos::RCP<const Epetra_Map> DofRowMap();
    Teuchos::RCP<Epetra_Operator> SysMat();
    Teuchos::RCP<DRT::Discretization> Discretization();
    double DispIncrFactor();

    virtual Teuchos::RCP<Epetra_Vector> FluidCondRHS();
    virtual Teuchos::RCP<Epetra_Vector> MeshCondRHS();

    //Teuchos::RCP<Epetra_Vector> InterfaceDispn();

    void PrepareTimeStep();
    void Evaluate(Teuchos::RCP<const Epetra_Vector> disp);
    void Update();
    void Output();

    Teuchos::RCP<Epetra_Map> InterfaceMap() { return interface_.CondMap(); }
    Teuchos::RCP<Epetra_Map> InnerDisplacementRowMap() { return interface_.OtherMap(); }

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const { return interface_; }

    virtual const Epetra_Map& DomainMap();

    virtual void ReadRestart(int step);

    virtual void UpdateandOutput();

    virtual void Solve();

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
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDisplacement();

    /// Predictor for interface displacements
    virtual Teuchos::RCP<Epetra_Vector> PredictInterfaceDisplacement();

    //@}

    //! @name Apply interface forces

    /// apply interface forces to structural solver
    /*!
      This prepares a new solve of the structural field within one time
      step. The middle values are newly created.

      \note This is not yet the most efficient implementation.
    */
    virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce);

    //@}

  private:

    /// the actual structure algorithm
    StruGenAlpha structure_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor interface_;

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
