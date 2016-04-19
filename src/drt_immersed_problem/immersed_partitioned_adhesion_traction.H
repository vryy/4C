/*!----------------------------------------------------------------------
\file immersed_partitioned_adhesion_traction.H
\level 2

\brief partitioned immersed cell-ecm interaction via adhesion traction

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

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

  /*!
  \brief Distribute force from cell adhesion node inside ecm element to ecm element nodes.

   Requires immersed cell discretization to have 'CellFocalAdhesion' defined in the input
   file (via '--DESIGN FOCAL ADHESION SURF CONDITIONS'). Requires the same surface/s to
   have also defined Dirichlet conditions to be able to retrieve the reaction forces from
   the immersed cell.

   1) Loop over all nodes with 'CellFocalAdhesion' condition.
     -> 1.1 get displacement of current node.

     -> 1.2 Loop over curr_subset_of_backgrounddis_ and search only in elements for which IsBoundaryImmersed() == true.
        -> 1.2.1 find backgrd. ECM element in which current node is immersed.
        -> 1.2.2 get parameter space coord. of current node in backgrd. element.
        -> 1.2.3 write pair 'cell adhesion node id' <-> 'xi' in adh_nod_param_coords_in_backgrd_ele_.
        -> 1.2.4 write pair 'adhesion node id' <-> 'backgrdele id' in map adh_nod_backgrd_ele_mapping_.
        -> 1.2.5 spread reaction force at current node to nodes of backgrd. element by the ratio of
                 shapefcts. value at current 'xi' and add those contributions to ecm_adhesion_forces_.

  \param adhesionforce  (in) : Force at adhesion nodes of immersed cell discretization

  \note This method uses member variables and cannot be used stand-alone

  \warning none
  \return void
  */
  void DistributeAdhesionForce(Teuchos::RCP<Epetra_Vector> adhesionforce);

  /*!
  \brief Interpolate backgrd. ecm element displacements to immersed cell adhesion node inside element.

   Requires immersed cell discretization to have 'CellFocalAdhesion' defined in the input
   file (via '--DESIGN FOCAL ADHESION SURF CONDITIONS'). Requires the same surface/s to
   have also defined Dirichlet conditions to be able to retrieve the reaction forces from
   the immersed cell.

   1) Loop over all nodes in adh_nod_param_coords_in_backgrd_ele_.
     -> 1.1 extract background ECM element nodal displacements.
     -> 1.2 interpolate ECM displacement to parameter space coordinate occupied by the cell adhesion node
            prior to the ECM deformation. Thereby the adhesion node is supposed to be fixed to the same
            material point.
     -> 1.3 save the interpolated displacements in cell_adhesion_disp_.

  \note this method uses member variables and cannot be used stand-alone

  \warning none
  \return void
  */
  void CalcAdhesionDisplacements();


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

  /// do we need an update of our immersed information (e.g. isimmersed_ in nodes)
  bool immersed_information_invalid_;

  /// write output after every fixed-point iteration
  bool output_evry_nlniter_;

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


  /// reset immersed information like e.g.: IsImmersed()
  void ResetImmersedInformation();

};// class ImmersedPartitionedAdhesionTraction
} // namespace IMMERSED

#endif
