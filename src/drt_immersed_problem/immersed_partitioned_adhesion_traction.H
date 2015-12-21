/*!----------------------------------------------------------------------
\file immersed_partitioned_adhesion_traction.H

\brief partitioned immersed cell-ecm interaction via adhesion traction

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_ADHESION_TRACTION_H
#define IMMERSED_PARTITIONED_ADHESION_TRACTION_H

#include "immersed_partitioned.H"
#include "immersed_field_exchange_manager.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  //class Monolithic;
  class PoroScatraBase;
}

namespace ADAPTER
{
  class FPSIStructureWrapper;
  class FSIStructureWrapperImmersed;
}

/*----------------------------------------------------------------------*/

namespace IMMERSED {

class ImmersedPartitionedAdhesionTraction : public ImmersedPartitioned
{

public:

  /// setup partitioned immersed cell migration algorithm
  explicit ImmersedPartitionedAdhesionTraction(const Teuchos::ParameterList& params, const Epetra_Comm& comm);

  /// read restart
  virtual void ReadRestart(int step);

protected:

  /// composed Field-Field-Interaction operator
  void CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// background operator
  void BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
      const FillType fillFlag);

  /// immersed operator
  Teuchos::RCP<Epetra_Vector> ImmersedOp(Teuchos::RCP<Epetra_Vector> bdry_traction,
      const FillType fillFlag);

  /// preparation of background field solve
  void PrepareBackgroundOp();

  /// prepare output of data to HD
  virtual void PrepareOutput();

  /// output of data to HD
  virtual void Output();

  /// update fields
  virtual void Update();

  /// prepare time step of single fields and print header
  virtual void PrepareTimeStep();

  /// initial guess
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  /// calculate force from cell adhesion node inside ecm element to ecm element nodes
  void CalcAdhesionForce();

  /// distribute force from cell adhesion node inside ecm element to ecm element nodes
  void DistributeAdhesionForce(Teuchos::RCP<Epetra_Vector> adhesionforce);

  /// interpolate ecm element displacements to cell adhesion node inside element
  void CalcAdhesionDisplacements();

  /// apply ecm element displacements to cell adhesion node inside element
  void ApplyAdhesionDisplacements();


private:

  /// basic information for parallelism
  int myrank_;
  int numproc_;

  /// pointer to global problem
  DRT::Problem* globalproblem_;

  /// immersed coupling variable
  bool displacementcoupling_;

  //! @name Various global vectors
  //@{
  Teuchos::RCP<Epetra_Vector> ecm_adhesion_forces_; //!< current nodal adhesion forces on ecm
  Teuchos::RCP<Epetra_Vector> cell_adhesion_disp_;  //!< current nodal cell adhesion displacements
  Teuchos::RCP<Epetra_Vector> Freact_cell_;         //!< current nodal reaction force because of dirichlet at adhesions
  //@}

  //! @name adhesion dynamics related map to save parameter space coords. of adhesion nodes
  //@{
  std::map<int, LINALG::Matrix<3,1> > adh_nod_param_coords_in_backgrd_ele_;
  std::map<int, int> adh_nod_backgrd_ele_mapping_;
  //@}

  //! @name Pointers to Discretizations
  //@{
  Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;
  Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;
  Teuchos::RCP<DRT::Discretization> immerseddis_;
  //@}

  //! @name Manager for data exchange between fields
  //@{
  DRT::ImmersedFieldExchangeManager* exchange_manager_ ;
  //@}

  //! @name Search Trees
  //@{
  Teuchos::RCP<GEO::SearchTree> fluid_SearchTree_;  //!< search tree for fluid domain
  Teuchos::RCP<GEO::SearchTree> cell_SearchTree_;   //!< search tree for cell domain
  //@}

  //! @name current nodal positions of the discretizations
  //@{
  std::map<int,LINALG::Matrix<3,1> >* currpositions_cell_; //!< pointer to map of vectors for search tree containing current structural positions
  std::map<int,LINALG::Matrix<3,1> >* currpositions_ECM_;  //!< pointer to map of vectors for search tree containing current fluid positions
  //@}

  //! @name current subsets of discretizations
  //@{
  std::map<int,std::set<int> > curr_subset_of_backgrounddis_;    //!< background elements covered by immersed discretisation
  //@}

  /// true if problem is modeled as a quasi 2D problem
  bool isPseudo2D_;

  /// cell structure of the immersed problem
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapperImmersed> cellstructure_;


  /// pointer to poroelast subproblem (poro-scatra interaction)
  Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

  /// poro structure of the immersed problem
  Teuchos::RCP< ::ADAPTER::FPSIStructureWrapper> porostructure_;

  /// reinitialize the transfer vectors
  void ReinitTransferVectors(){  ecm_adhesion_forces_->PutScalar(0.0);
                                 cell_adhesion_disp_ ->PutScalar(0.0);
                                 return;
                              };

  /// get immersed nodes and determine their dofs
  void BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis,
                              Teuchos::RCP<Epetra_Map>& dirichmap,
                              const Teuchos::RCP<const Epetra_Map>& dirichmap_original,
                              int dofsetnum);

  /// add immersed dirichlet values from immersed dis to systemvector of background dis
  void DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector,
                               Teuchos::RCP<Epetra_Vector> dirichvals,
                               Teuchos::RCP<const Epetra_Map> dbcmap);

};// class ImmersedPartitionedAdhesionTraction
} // namespace IMMERSED

#endif
