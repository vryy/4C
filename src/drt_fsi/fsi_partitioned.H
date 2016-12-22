/*----------------------------------------------------------------------*/
/*!
 * \file fsi_partitioned.H
\brief Partitioned FSI base

\level 1

\maintainer Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362

*/
/*----------------------------------------------------------------------*/



#ifndef FSI_PARTITIONED_H
#define FSI_PARTITIONED_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <AztecOO.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_LinearSystem_AztecOO.H>

#include "fsi_algorithm.H"

// forward declarations
namespace ADAPTER
{
  class CouplingMortar;
}


/// FSI: Fluid Structure Interaction
/*!
  This is the place where we bring fluid and structure together.
 */
namespace FSI {

namespace UTILS {
  class DebugWriter;
}

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
  explicit Partitioned(const Epetra_Comm& comm);

  /// setup this object
  virtual void Setup();

  /// outer level FSI time loop
  virtual void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

  /// compute FSI interface residual S^{-1}(F(d)) - d
  bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// return true if nodes at interface are matching
  bool matchingnodes() {return matchingnodes_; }

  /// Open door in the time loop for sliding ale algo to do remeshing
  virtual void Remeshing();

  /// Setup of coupling at fsi interface
  virtual void SetupCoupling(const Teuchos::ParameterList& fsidyn ,const Epetra_Comm& comm);

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

  /*! \brief Calculate interface velocity based on given interface displacements
   *
   *  Two options to transform the structural interface displacement into a fluid interface velocity:
   *  - second order (cf. eq. (6.1.2) in [2])               --> set input parameter SECONDORDER = Yes
   *  - fisrt order Backward Euler (cf. eq. (6.1.2) in [2]) --> set input parameter SECONDORDER = No
   *
   *  A derivation of these kinematic coupling conditions is given in chapter 6.2.3 in [1].
   *
   *  References:
   *  - [1] C Foerster, Robust methods for fluid-structure interaction with stabilised finite elemtes,
   *    PhD-Thesis, 2007
   *  - [2] U Kuettler, Effiziente Loesungsverfahren fuer Fluid-Struktur-Interaktions-Probleme,
   *    PhD-Thesis, 2009
   */
  Teuchos::RCP<Epetra_Vector> InterfaceVelocity(Teuchos::RCP<const Epetra_Vector> idispnp) const;

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

  /// create convergence tests
  virtual void CreateStatusTest(Teuchos::ParameterList& nlParams,
                                Teuchos::RCP<NOX::Epetra::Group> grp,
                                Teuchos::RCP<NOX::StatusTest::Combo> converged);

  Teuchos::RCP<UTILS::DebugWriter> MyDebugWriter() const { return debugwriter_; }

  /// return coupsfm_
  ADAPTER::CouplingMortar& StructureFluidCouplingMortar();

  /// return coupsfm_
  const ADAPTER::CouplingMortar& StructureFluidCouplingMortar() const;

  /// access to iteration counter
  virtual std::vector<int> IterationCounter(){return counter_;};

  /// extract idispn_ iveln_
  virtual void ExtractPreviousInterfaceSolution();

  /// setup list with default parameters
  virtual void SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

private:

  /// create linear solver framework
  Teuchos::RCP<NOX::Epetra::LinearSystem>
  CreateLinearSystem(Teuchos::ParameterList& nlParams,
                     const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface,
                     NOX::Epetra::Vector& noxSoln,
                     Teuchos::RCP<NOX::Utils> utils);

  /// create convergence tests including testing framework
  Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(Teuchos::ParameterList& nlParams,
                                                        Teuchos::RCP<NOX::Epetra::Group> grp);

  /// interface displacement from time step begin
  Teuchos::RCP<Epetra_Vector> idispn_;

  /// interface velocity from time step begin
  Teuchos::RCP<Epetra_Vector> iveln_;

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

protected:
  int mfresitemax_;


  /// coupling of structure and fluid at the interface, with mortar.
  Teuchos::RCP<ADAPTER::CouplingMortar>         coupsfm_;

  /// nodes at the fluid-structure interface match
  bool matchingnodes_;

  /// parameters handed in to NOX
  Teuchos::ParameterList noxparameterlist_;

  /// special debugging output
  Teuchos::RCP<UTILS::DebugWriter> debugwriter_;

};

}

#endif
