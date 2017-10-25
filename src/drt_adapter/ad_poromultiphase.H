/*----------------------------------------------------------------------*/
/*!
 \file ad_poromultiphase.H

 \brief

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_POROMULTIPHASE_H_
#define SRC_DRT_ADAPTER_AD_POROMULTIPHASE_H_

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>
#include "../linalg/linalg_utils.H"

namespace ADAPTER
{
  class PoroFluidMultiphaseWrapper;
  class Structure;
}

namespace LINALG
{
  class Solver;
}

namespace ADAPTER
{
class PoroMultiPhase
{
public:

  /// constructor
  PoroMultiPhase(){};

  /// virtual destructor
  virtual ~PoroMultiPhase(){};

  /// initialization
  virtual void Init(
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& algoparams,
      const Teuchos::ParameterList& structparams,
      const Teuchos::ParameterList& fluidparams,
      const std::string& struct_disname,
      const std::string& fluid_disname,
      bool isale,
      int nds_disp,
      int nds_vel,
      int nds_solidpressure,
      int ndsporofluid_scatra) =0;

  /// read restart
  virtual void ReadRestart( int restart ) =0;

  /// test results (if necessary)
  virtual void CreateFieldTest() =0;

  /// setup
  virtual void SetupSystem() = 0;

  /// setup the solver (only for monolithic system)
  virtual bool SetupSolver() = 0;

  /// prepare timeloop of coupled problem
  virtual void PrepareTimeLoop() = 0;

  /// timeloop of coupled problem
  virtual void Timeloop() = 0;

  /// time step of coupled problem
  virtual void TimeStep() = 0;

  /// time step of coupled problem
  virtual void PrepareTimeStep() = 0;

  //! update time step and print to screen
  virtual void UpdateAndOutput() = 0;

  /// set structure solution on scatra field
  virtual void SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                          Teuchos::RCP<const Epetra_Vector> vel ) =0;

  /// set scatra solution on fluid field
  virtual void SetScatraSolution(unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars ) =0;

  /// dof map of vector of unknowns
  virtual Teuchos::RCP<const Epetra_Map> StructDofRowMap() const = 0;

  /// unknown displacements at \f$t_{n+1}\f$
  virtual Teuchos::RCP<const Epetra_Vector> StructDispnp() const =0;

  /// unknown velocity at \f$t_{n+1}\f$
  virtual Teuchos::RCP<const Epetra_Vector> StructVelnp() const =0;

  /// dof map of vector of unknowns
  virtual Teuchos::RCP<const Epetra_Map> FluidDofRowMap() const = 0;

  /// return fluid flux
  virtual Teuchos::RCP<const Epetra_MultiVector> FluidFlux() const =0;

  /// return fluid solution variable
  virtual Teuchos::RCP<const Epetra_Vector> FluidPhinp() const =0;

  /// return fluid solution variable
  virtual Teuchos::RCP<const Epetra_Vector> FluidSaturation() const =0;

  /// return fluid solution variable
  virtual Teuchos::RCP<const Epetra_Vector> FluidPressure() const =0;

  /// return fluid solution variable
  virtual Teuchos::RCP<const Epetra_Vector> SolidPressure() const =0;

  //! unique map of all dofs that should be constrained with DBC
  virtual Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const = 0;

  //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
  virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> sx,
      Teuchos::RCP<const Epetra_Vector> fx,
      const bool firstcall) = 0;

  //! access to monolithic right-hand side vector
  virtual Teuchos::RCP<const Epetra_Vector> RHS() const = 0;

  //! update all fields after convergence (add increment on displacements and fluid primary variables)
  virtual void UpdateFieldsAfterConvergence(
      Teuchos::RCP<const Epetra_Vector>& sx,
      Teuchos::RCP<const Epetra_Vector>& fx
      ) = 0;

  //! get the extractor
  virtual Teuchos::RCP<const LINALG::MultiMapExtractor> Extractor() const = 0;

  //! get the monolithic system matrix
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const = 0;

  //! get structure field
  virtual const Teuchos::RCP<ADAPTER::Structure>& StructureField() = 0;

  //! get fluid field
  virtual const Teuchos::RCP<ADAPTER::PoroFluidMultiphaseWrapper>& FluidField() = 0;

  //! build the block null spaces
  virtual void BuildBlockNullSpaces(Teuchos::RCP<LINALG::Solver>& solver) = 0;

};
}


#endif /* SRC_DRT_ADAPTER_AD_POROMULTIPHASE_H_ */
