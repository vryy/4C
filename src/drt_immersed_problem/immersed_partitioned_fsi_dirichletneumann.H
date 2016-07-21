  /*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann.H

\brief partitioned immersed fsi algorithm for neumann-(dirichlet-neumann) like coupling

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_H
#define IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_H

#include "immersed_base.H"
#include "fsi_partitioned_immersed.H"


namespace FLD
{
  class FluidImplicitTimeInt;
}

namespace GEO
{
  class SearchTree;
}

namespace Teuchos
{
  class Time;
}

namespace ADAPTER
{
  class FSIStructureWrapperImmersed;
}

namespace IMMERSED {

class ImmersedPartitionedFSIDirichletNeumann : public ImmersedBase, public FSI::PartitionedImmersed
{

public:

  /// setup partitioned ImmersedFSI algorithm
  explicit ImmersedPartitionedFSIDirichletNeumann(const Epetra_Comm& comm);

  /// initialize search tree for structure discretization and ghost structure redundantly
  void SetupStructuralDiscretization();

protected:

  /// composed FSI operator
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> bforce,
                                              const FillType fillFlag);

  /// interface structural operator
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> struct_bdry_traction,
                                               const FillType fillFlag);
  /// initial guess
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  /// get immersed nodes and determine their dofs
  void BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& dirichmap, const Teuchos::RCP<const Epetra_Map>& dirichmap_original);

  /// add immersed dirichlet values from immersed dis to systemvector of background dis
  void DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap);

  /// set state necessary state vectors
  virtual void SetStatesFluidOP();

  /// set state necessary state vectors
  virtual void SetStatesVelocityCorrection();

  /// set state necessary state vectors
  virtual void SetStatesStructOP();

  /// update nodal positions of fluid nodes (stored in this class)
  virtual void UpdateCurrentPositionsFluidNodes(){/* fixed grid */ return;};

  /// call the nonlinear fluid solver
  virtual void SolveFluid();

  /// call the nonlinear structural solver
  virtual void SolveStruct();

  /// determine elements cut by the boundary
  void PrepareFluidOp();

  /// call to special extraction method
  virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp();

  /// call to special application of interface forces
  virtual void ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> full_traction_vec);

  /// call to special routine that adds dirichlet values to fluid field
  virtual void AddDirichCond();

  /// call to special routine that removes dirichlet values from fluid field
  virtual void RemoveDirichCond();

  /// calc the fsi residual
  virtual int CalcResidual(Epetra_Vector &F, const Teuchos::RCP<Epetra_Vector> newstate , const Teuchos::RCP<Epetra_Vector> oldstate);

  /// calc the current fluid tractions interpolated to the structural surface
  virtual void CalcFluidTractionsOnStructure();


  //! @name Various global forces
  //@{
  Teuchos::RCP<Epetra_Vector> struct_bdry_traction_;     //!< bdry traction rhs on struct FSI surface
  Teuchos::RCP<Epetra_Vector> fluid_artificial_velocity_;//!< background velocity interpolated from immersed dis
  Teuchos::RCP<Epetra_Vector> fluid_artificial_velocity_old_;//!< background velocity interpolated from immersed dis
  Teuchos::RCP<Epetra_Map> dbcmap_immersed_;             //!< dirichlet bc map of immersed values
  Teuchos::RCP<GEO::SearchTree> fluid_SearchTree_;       //!< search tree for fluid domain
  Teuchos::RCP<GEO::SearchTree> structure_SearchTree_;   //!< search tree for structure domain
  std::map<int,std::set<int> > curr_subset_of_fluiddis_;    //!< fluid elements to evaluate the dirichlet interpolation
  std::map<int,LINALG::Matrix<3,1> > currpositions_struct_; //!< map of vectors for search tree containing current structural positions
  std::map<int,LINALG::Matrix<3,1> > currpositions_fluid_;  //!< map of vectors for search tree containing current fluid positions
  //@}

  //! @name  basic information for parallelism
  int myrank_;
  int numproc_;
  //@}

  /// pointer to global problem
  DRT::Problem* globalproblem_;

  //! @name various bools and switches
  //@{
  bool displacementcoupling_;
  bool multibodysimulation_;
  bool output_evry_nlniter_;
  int coupalgo_;
  int correct_boundary_velocities_;
  int degree_gp_fluid_bound_;
  //@}

  //! @name pointer to discretizations
  //@{
  Teuchos::RCP<DRT::Discretization> fluiddis_;
  Teuchos::RCP<DRT::Discretization> structdis_;
  //@}

  /// pointer to immersed structure adapter
  Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> immersedstructure_;

};// class ImmersedFSI
} // namspace IMMERSED

#endif
