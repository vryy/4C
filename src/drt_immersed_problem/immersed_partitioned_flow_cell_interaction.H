  /*!----------------------------------------------------------------------
\file immersed_partitioned_flow_cell_interaction.H

\brief partitioned immersed cell - interstitial flow interaction algorithm

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_FLOW_CELL_INTERACTION_H
#define IMMERSED_PARTITIONED_FLOW_CELL_INTERACTION_H

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
  class AdapterScatraWrapperCellMigration;
  class FPSIStructureWrapper;
  class FSIStructureWrapperImmersed;
}

namespace SSI
{
  class SSI_Base;
}

/*----------------------------------------------------------------------*/

namespace IMMERSED {

class ImmersedPartitionedFlowCellInteraction : public ImmersedPartitioned
{

public:

  /// setup partitioned immersed cell migration algorithm
  explicit ImmersedPartitionedFlowCellInteraction(const Teuchos::ParameterList& params, const Epetra_Comm& comm);

  /*! \brief Initialize this object

  \warning none
  \return int
  \date 08/16
  \author rauch  */
  virtual int Init(const Teuchos::ParameterList& params);

  /*! \brief Setup all class internal objects and members

  \warning none
  \return void
  \date 08/16
  \author rauch  */
  virtual void Setup();

  /// access to field cell
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapperImmersed> CellField(){return cellstructure_;};

  /// access to field Poro
  Teuchos::RCP<POROELAST::PoroScatraBase> PoroField(){return poroscatra_subproblem_;};

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

  /// extract the pressure part
  void ExtractPressurePart();

  /// calc the current fluid tractions interpolated to the structural surface
  void CalcFluidTractionOnStructure();

  //! @name Various global vectors
  //@{
  Teuchos::RCP<Epetra_Vector> porofluid_artificial_velocity_;//!< background velocity interpolated from immersed dis
  Teuchos::RCP<Epetra_Vector> poroscatra_segregated_phi_;    //!< background phi interpolated from immersed dis
  Teuchos::RCP<Epetra_Vector> cell_bdry_traction_;           //!< boundary traction rhs on immersed cell
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
  Teuchos::RCP< ::ADAPTER::FSIStructureWrapperImmersed> cellstructure_;

  /// poro structure of the immersed problem
  Teuchos::RCP< ::ADAPTER::FPSIStructureWrapper>        porostructure_;

  /// pointer to cell subproblem (structure-scatra interaction)
  Teuchos::RCP<SSI::SSI_Base> cell_subproblem_;

  /// pointer to poroelast subproblem (poro-scatra interaction)
  Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

  /// scatra adapter for poro-scatra time integration
  Teuchos::RCP<ADAPTER::AdapterScatraWrapperCellMigration> poroscatra_wrapper_;

  /// basic information for parallelism
  int myrank_;
  int numproc_;

  /// pointer to global problem
  DRT::Problem* globalproblem_;

  /// coupling scheme
  bool displacementcoupling_;

  /// switch that determines if there are more than one distinct immersed cells
  bool multicellmigration_;  //! distinct cells must be marked with ImmersedSearchbox condition

  /// true if problem is modeled as a quasi 2D problem
  bool isPseudo2D_;

  /// amoeboid or proteolytic migration
  int migrationtype_;

  /// global timstep
  double timestep_;

  //! @name proteolytic migration control parameters
  //@{
  double segregationconstant_;
  int segregationtype_;
  int segregationby_;
  //@}

  /// counter for number of unconverged timesteps
  int continued_steps_;

  /// number of integration points related to mass conserving adaption of equations
  int degree_gp_fluid_bound_; //!< number of integration points in fluid elements cut by structural boundary

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
  std::map<int,std::set<int> > curr_subset_of_backgrounddis_;    //!< backgrd elements to evaluate the dirichlet interpolation
  //@}

  //! @name Maps for Dirichlet BCs
  //@{
  Teuchos::RCP<Epetra_Map> dbcmap_immersed_;             //!< dirichlet bc map of immersed values
  Teuchos::RCP<Epetra_Map> dbcmap_immersed_scatra_;      //!< dirichlet bc map of immersed scatra values
  //@}

private:

  /// reinitialize the transfer vectors
  void ReinitTransferVectors(){  porofluid_artificial_velocity_->PutScalar(0.0);
                                 cell_bdry_traction_->PutScalar(0.0);
                                 poroscatra_segregated_phi_->PutScalar(0.0);
                                 return;
                              };

  /// reset element and node information about immersed method
  void ResetImmersedInformation();

};// class ImmersedPartitionedFlowCellInteraction
} // namespace IMMERSED

#endif
