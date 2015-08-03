  /*!----------------------------------------------------------------------
\file immersed_partitioned_cellmigration.H

\brief partitioned immersed fsi algorithm for dirichlet-neumann coupling (volume dirichlet coupling)

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_CELLMIGRATION_H
#define IMMERSED_PARTITIONED_CELLMIGRATION_H

#include "immersed_partitioned.H"
#include "immersed_field_exchange_manager.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  //class Monolithic;
  class PORO_SCATRA_Base;
}

namespace ADAPTER
{
  class FSIStructureWrapper;
}

/*----------------------------------------------------------------------*/

namespace IMMERSED {

class ImmersedPartitionedCellMigration : public ImmersedPartitioned
{

public:

  /// setup partitioned ImmersedFSI algorithm
  explicit ImmersedPartitionedCellMigration(const Epetra_Comm& comm);

  /// initialize search tree for structure discretization and ghost structure redundantly
  void SetupImmersedDiscretization();

  /// initialize search tree for poro discretization
  void SetupBackgroundDiscretization();

  /// access to field cell
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapper> CellField(){return cellstructure_;};

  /// access to field Poro
  //Teuchos::RCP<POROELAST::Monolithic> PoroField(){return poroelast_subproblem_;};
  Teuchos::RCP<POROELAST::PORO_SCATRA_Base> PoroField(){return poroelast_subproblem_;};

  /// read restart
  virtual void ReadRestart(int step);

protected:

  /// composed FSI operator
  void CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// background operator
  void BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
                    const FillType fillFlag);

  /// immersed operator
  Teuchos::RCP<Epetra_Vector> ImmersedOp(Teuchos::RCP<Epetra_Vector> bdry_traction,
                                         const FillType fillFlag);
  /// initial guess
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  /// get immersed nodes and determine their dofs
  void BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis,
                              Teuchos::RCP<Epetra_Map>& dirichmap,
                              const Teuchos::RCP<const Epetra_Map>& dirichmap_original,
                              int dofsetnum);

  /// get immersed nodes and determine their dofs
  void BuildImmersedScaTraDirichMap(Teuchos::RCP<DRT::Discretization> immersedinfodis, //< dis carrying information on artificial domain
                                    Teuchos::RCP<DRT::Discretization> dis,             //< dis for which the dirichlet map is to be modified
                                    Teuchos::RCP<Epetra_Map>& dirichmap,               //< new dirichmap to be filled
                                    const Teuchos::RCP<const Epetra_Map>& dirichmap_original, //< original dirich
                                    int dofsetnum);                                           //< dofsetnumber of dis node dofs to be retrieved

  /// prepare time step of single fields and print header
  virtual void PrepareTimeStep();

  /// add immersed dirichlet values from immersed dis to systemvector of background dis
  void DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap);

  /// zero out the cell-ecm inteaction forces at the leading edge
  //void DeleteForceAtLeadingEdge(Teuchos::RCP<Epetra_Vector> force);

  /// prepare the immersed solve
  void PrepareImmersedOp();

  /// prepare the background solve
  void PrepareBackgroundOp();

  /// prepare output of data to HD
  virtual void PrepareOutput();

  /// output of data to HD
  virtual void Output();

  /// update fields
  virtual void Update();

  /// initialize the relaxation parameters and print output
  void SetupRelaxation();

  /// setup poro algorithm (coupling, combined dofmaps, etc.)
  void SetupPoroAlgorithm();

  /// extract the pressure part
  void ExtractPressurePart();

  /// extract the pressure part
  virtual void SetFieldDt();

  /// check input vector for negatives and set them zero
  //void PreventInterfacePullingTraction(Teuchos::RCP<Epetra_Vector> increment_vector,
  //                                     Teuchos::RCP<Epetra_Vector> total_vector);

  //! @name Various global vectors
  //@{
  Teuchos::RCP<Epetra_Vector> porofluid_artificial_velocity_;//!< background velocity interpolated from immersed dis
  Teuchos::RCP<Epetra_Vector> porofluid_artificial_source_;  //!< background mass change interpolated from immersed dis
  Teuchos::RCP<Epetra_Vector> poroscatra_segregated_phi_;    //!< background phi interpolated from immersed dis
  Teuchos::RCP<Epetra_Vector> cell_bdry_traction_;           //!< boundary traction rhs on immersed cell
  Teuchos::RCP<Epetra_Vector> cell_penalty_traction_;        //!< boundary traction rhs on immersed cell from cell-ecm interaction
  Teuchos::RCP<Epetra_Vector> cell_delta_penalty_traction_;  //!< boundary traction rhs on immersed cell from cell-ecm interaction
  Teuchos::RCP<Epetra_Vector> cell_delta_penalty_traction_old_;  //!< boundary traction rhs on immersed cell from cell-ecm interaction of previous iteration
  Teuchos::RCP<Epetra_Vector> cell_bdry_traction_old_;       //!< old boundary traction rhs on immersed cell
  Teuchos::RCP<Epetra_Vector> penalty_gap_;
  Teuchos::RCP<Epetra_Vector> curr_nodal_normals_;
  //@}

  //! @name Pointers to Discretizations
    //@{
  Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;
  Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;
  Teuchos::RCP<DRT::Discretization> immerseddis_;
  Teuchos::RCP<DRT::Discretization> scatradis_;
  //@}

  //! @name parameters and variables for relaxation
    //@{
  Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory_;
  bool relaxforceglobally_;
  bool relaxvelglobally_;
  double forcerelax_;
  double velrelax_;
  int coupalgo_;
  //@}

  /// cell structure of the immersed problem
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    cellstructure_;

  /// poroelast subproblem
  //Teuchos::RCP<POROELAST::Monolithic>       poroelast_subproblem_;
  Teuchos::RCP<POROELAST::PORO_SCATRA_Base> poroelast_subproblem_;

  /// basic information for parallelism
  int myrank_;
  int numproc_;

  /// pointer to global problem
  DRT::Problem* globalproblem_;

  /// various other parameters
  bool displacementcoupling_;
  bool multicellmigration_;
  int migrationtype_;
  double segregationconstant_;
  int segregationtype_;
  int segregationby_;
  double penalty_start_; //!< start penalty parameter at the beginning of each timestep
  double penalty_init_;  //!< init penalty parameter
  bool initialize_cell_;
  double timestep_;
  double initial_timestep_;
  int initialization_steps_;
  bool ecm_interaction_;
  int continued_steps_;
  double nu_; //<- AITKEN relaxation parameter


  DRT::ImmersedFieldExchangeManager* exchange_manager_ ;

  //! @name Search Trees
    //@{
    Teuchos::RCP<GEO::SearchTree> fluid_SearchTree_;  //!< search tree for fluid domain
    Teuchos::RCP<GEO::SearchTree> cell_SearchTree_;   //!< search tree for cell domain
    //@}

  //! @name current nodal positions of the discretizations
    //@{
    std::map<int,LINALG::Matrix<3,1> > currpositions_cell_; //!< map of vectors for search tree containing current structural positions
    std::map<int,LINALG::Matrix<3,1> > currpositions_ECM_;  //!< map of vectors for search tree containing current fluid positions
   //@}

    //! @name current subsets of discretizations
      //@{
    std::map<int,std::set<int> > curr_subset_of_backgrounddis_;    //!< backgrd elements to evaluate the dirichlet interpolation
    //@}

    //! @name current subsets of discretizations
      //@{
    Teuchos::RCP<const Epetra_Map> dbcmap_ordinary_;       //!< dirichlet bc map before consideration of immersed values
    Teuchos::RCP<Epetra_Map> dbcmap_immersed_;             //!< dirichlet bc map of immersed values
    Teuchos::RCP<Epetra_Map> dbcmap_immersed_scatra_;      //!< dirichlet bc map of immersed scatra values
    //@}

};// class ImmersedFSI
} // namspace IMMERSED

#endif
