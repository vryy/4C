/*----------------------------------------------------------------------*/
/*!
\brief Partitioned FSI base

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FSI_PARTITIONED_H
#define FSI_PARTITIONED_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <AztecOO.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

#include <vector>

#include "fsi_algorithm.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_mortar.H"


/// FSI: Fluid Structure Interaction
/*!
  This is the place where we bring fluid and structure together.
 */
namespace FSI {

  class DebugWriter;

/// Partitioned FSI base
/*!

  This is the framework for partitioned FSI. The actual work is done by
  subclasses.

  \author u.kue
  \date 04/08
 */
class Partitioned : public Algorithm,
                    public NOX::Epetra::Interface::Required
{
public:

  /// setup FSI algorithm
  explicit Partitioned(Epetra_Comm& comm);

  /// outer level FSI time loop
  void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

  /// compute FSI interface residual S^{-1}(F(d)) - d
  bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

protected:

  //! @name Transfer helpers

  Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv);
  Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv);

  //@}

  //! @name Operators implemented by subclasses

  /// composed FSI operator
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                                              const FillType fillFlag);

  /// interface structural operator
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                               const FillType fillFlag);

  //@}

  //! @name Encapsulation of interface unknown
  /// default is displacement, but subclasses might change that

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  //@}

  //! @name Access methods for subclasses

  /// interface velocity calculation given the interface displacements
  Teuchos::RCP<Epetra_Vector> InterfaceVelocity(const Teuchos::RCP<Epetra_Vector> idispnp) const;

  /// current interface displacements
  /*!
    Extract structural displacement at t(n+1)
   */
  Teuchos::RCP<Epetra_Vector> InterfaceDisp();

  /// current interface forces
  /*!
    Extract fluid force at t(n+1)
   */
  Teuchos::RCP<Epetra_Vector> InterfaceForce();

  //@}

  /// create convergance tests
  virtual void CreateStatusTest(ParameterList& nlParams,
                                Teuchos::RCP<NOX::Epetra::Group> grp,
                                Teuchos::RCP<NOX::StatusTest::Combo> converged);

private:

  /// create linear solver framework
  Teuchos::RCP<NOX::Epetra::LinearSystem>
  CreateLinearSystem(ParameterList& nlParams,
                     const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface,
                     NOX::Epetra::Vector& noxSoln,
                     Teuchos::RCP<NOX::Utils> utils);

  /// create convergance tests including testing framework
  Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(ParameterList& nlParams,
                                                        Teuchos::RCP<NOX::Epetra::Group> grp);

  /// setup list with default parameters
  void SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

  /// coupling of structure and fluid at the interface, with mortar. Experimental!
  ADAPTER::CouplingMortar         coupsfm_;

  /// interface displacement from time step begin
  Teuchos::RCP<Epetra_Vector> idispn_;

  //! connection of interface dofs for finite differences
  Teuchos::RCP<Epetra_CrsGraph> rawGraph_;

  //! counters on how many times the residuum was called in a time step
  /*!
    NOX knows different types of residuum calls depending on
    circumstances (normal, finite difference, matrix free
    jacobi). It is possible to do approximations depending on the
    type.
   */
  std::vector<int> counter_;

  //! number of residuum calculations per nonlinear solve in one time step
  std::vector<int> linsolvcount_;

  /// print parameters and stuff
  /*!
    \warning This variable is only valid during when the time loop runs.
   */
  Teuchos::RCP<NOX::Utils> utils_;

  /// parameters handed in to NOX
  Teuchos::ParameterList noxparameterlist_;

  int mfresitemax_;

  /// nodes at the fluid-structure interface match
  bool matchingnodes_;

  /// special debugging output
  Teuchos::RCP<DebugWriter> debugwriter_;
};

}

#endif
#endif
