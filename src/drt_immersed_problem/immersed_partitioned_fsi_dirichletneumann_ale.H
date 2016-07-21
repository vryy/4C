  /*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann_ale.H

\brief partitioned immersed fsi algorithm for dirichlet-neumann coupling combined with ALE

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_ALE_H
#define IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_ALE_H

#include "immersed_partitioned_fsi_dirichletneumann.H"


namespace IMMERSED {

class ImmersedPartitionedFSIDirichletNeumannALE : public ImmersedPartitionedFSIDirichletNeumann
{

public:

  /// setup partitioned ImmersedFSI algorithm
  explicit ImmersedPartitionedFSIDirichletNeumannALE(const Epetra_Comm& comm);

  /// here we need to call the standard FSI setup
  virtual void SetupCoupling(const Teuchos::ParameterList& fsidyn ,const Epetra_Comm& comm);

  /// set necessary pointers to state vectors
  virtual void SetStatesFluidOP();

  /// set state necessary state vectors
  virtual void SetStatesVelocityCorrection();

  /// set necessary pointers to state vectors
  virtual void SetStatesStructOP();

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> bforce,
                                              const FillType fillFlag);

  /// interface structural operator
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> struct_bdry_traction,
                                               const FillType fillFlag);

  /// initial guess
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  /// update nodal positions of fluid nodes (stored in this class)
  virtual void UpdateCurrentPositionsFluidNodes();

  /// don't call the nonlinear fluid solver here
  virtual void SolveFluid();

  /// don't call the nonlinear structural solver here
  virtual void SolveStruct(){ return;};

  /// call to special extraction method
  virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp();

  /// call to special application of interface forces
  virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> full_traction_vec);

  /// call to FSI ExtractPreviousInterfaceSolution()
  virtual void ExtractPreviousInterfaceSolution();

  /// call to special routine that adds dirichlet values to fluid field
  virtual void AddDirichCond();

  /// call to special routine that removes dirichlet values from fluid field
  virtual void RemoveDirichCond();

  /// calc the fsi residual
  virtual int CalcResidual(Epetra_Vector &F, const Teuchos::RCP<Epetra_Vector> newstate , const Teuchos::RCP<Epetra_Vector> oldstate);

  /// calc the current fluid tractions interpolated to the structural surface
  virtual void CalcFluidTractionsOnStructure();

  /// setup list with default parameters
  virtual void SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

  /// Immersed + FSI interface vector
  Teuchos::RCP<Epetra_Vector> combined_newstate_;

  //! @name Various interface states
  //@{
  Teuchos::RCP<Epetra_Vector> idispnp_;
  Teuchos::RCP<Epetra_Vector> ivelnp_;
  //@}

};// class ImmersedFSI
} // namspace IMMERSED

#endif
