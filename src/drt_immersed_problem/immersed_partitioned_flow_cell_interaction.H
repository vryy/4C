/*!----------------------------------------------------------------------

\brief partitioned immersed cell - interstitial flow interaction algorithm

\level 2

\maintainer Jonas Eichinger
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
  // class Monolithic;
  class PoroScatraBase;
}  // namespace POROELAST

namespace ADAPTER
{
  class AdapterScatraWrapperCellMigration;
  class FPSIStructureWrapper;
  class FSIStructureWrapperImmersed;
}  // namespace ADAPTER

namespace SSI
{
  class SSI_Base;
}

/*----------------------------------------------------------------------*/

namespace IMMERSED
{
  class ImmersedPartitionedFlowCellInteraction : public ImmersedPartitioned
  {
   public:
    /// setup partitioned immersed cell migration algorithm
    explicit ImmersedPartitionedFlowCellInteraction(
        const Teuchos::ParameterList& params, const Epetra_Comm& comm);

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
    Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed> CellField() { return cellstructure_; };

    /// access to field Poro
    Teuchos::RCP<POROELAST::PoroScatraBase> PoroField() { return poroscatra_subproblem_; };

    /// read restart
    virtual void ReadRestart(int step);

    /// write drag on immersed solid to file
    void WriteDrag();

   protected:
    /// composed Field-Field-Interaction operator
    void CouplingOp(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

    /// background operator
    void BackgroundOp(
        Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values, const FillType fillFlag);

    /// immersed operator
    Teuchos::RCP<Epetra_Vector> ImmersedOp(
        Teuchos::RCP<Epetra_Vector> bdry_traction, const FillType fillFlag);
    /// initial guess
    virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

    /// get immersed nodes and determine their dofs
    void BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map>& dirichmap,
        const Teuchos::RCP<const Epetra_Map>& dirichmap_original, int dofsetnum);

    /// get immersed nodes and determine their dofs
    void BuildImmersedScaTraDirichMap(
        Teuchos::RCP<DRT::Discretization>
            immersedinfodis,  //< dis carrying information on artificial domain
        Teuchos::RCP<DRT::Discretization>
            dis,                              //< dis for which the dirichlet map is to be modified
        Teuchos::RCP<Epetra_Map>& dirichmap,  //< new dirichmap to be filled
        const Teuchos::RCP<const Epetra_Map>& dirichmap_original,  //< original dirich
        int dofsetnum);  //< dofsetnumber of dis node dofs to be retrieved

    /// prepare time step of single fields and print header
    virtual void PrepareTimeStep();

    /// add immersed dirichlet values from immersed dis to systemvector of background dis
    void DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector,
        Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap);

    /// set state necessary state vectors
    virtual void SetStatesBackgroundOP();

    /// set state necessary state vectors
    virtual void SetStatesVelocityCorrection();

    /// set state necessary state vectors
    virtual void SetStatesImmersedOP();

    /// update nodal positions of ECM nodes
    void UpdateCurrentPositionsBackgroundNodes();

    /// prepare the background solve
    void PrepareBackgroundOp();

    /// prepare output of data to HD
    virtual void PrepareOutput();

    /// output of data to HD
    virtual void Output();

    /// update fields
    virtual void Update();

    /// apply given vector as Dirichlet to artificial fluid field
    virtual void ApplyImmersedDirichlet(Teuchos::RCP<Epetra_Vector> artificial_velocity);

    /// apply given vector as Dirichlet to scatra field
    virtual void ApplyImmersedDirichletScatra(Teuchos::RCP<Epetra_Vector> scatra_values);

    /// call to special routine that removes dirichlet values from fluid field
    virtual void RemoveDirichCond();

    /// call to special routine that removes scatra values from scatra field
    virtual void RemoveDirichCondScatra();

    /*!
    \brief calc the current fluid tractions integrated over structural surface

    Calc new tractions on immersed boundary. After leaving this method the vector
    struct_bdry_traction_ contains the current tractions on the immersed boundary.

    */
    void CalcFluidTractionOnStructure();

    /*!
    \brief calc the current artificial velocity by projection from structure velocity onto fluid

    Calc new velocity in artificial fluid domain. After leaving this method the vector
    fluid_artificial_velocity_ contains the current velocity projected from the structure onto the
    fluid.

     \note Artificial velocity is only projected if the member artificial_velocity_isvalid_ is=false
    on entry. if artificial_velocity_isvalid_=true the current artificial velocity stored in
    fluid_artificial_velocity_ is returned. This is important when restart is requested since then
    the we need to perform the projection in the first restart step. This is done in \ref
    InitialGuess(). In case of no restart, we can just return the velocity since it had recently
    been projected after the last structural solve. This work would be done twice, else.
    */
    virtual Teuchos::RCP<Epetra_Vector> CalcArtificialVelocity();

    /// return vector containing information for inter-module coupling
    Teuchos::RCP<Epetra_Vector> ReturnCouplingInfo();

    /// print information on the current solution step
    virtual void PrintStepInfo();

    //! @name Various global vectors
    //@{
    Teuchos::RCP<Epetra_Vector>
        porofluid_artificial_velocity_;  //!< background velocity interpolated from immersed dis
    Teuchos::RCP<Epetra_Vector>
        poroscatra_segregated_phi_;  //!< background phi interpolated from immersed dis
    Teuchos::RCP<Epetra_Vector> cell_bdry_traction_;  //!< boundary traction rhs on immersed cell
    Teuchos::RCP<Epetra_Vector>
        cell_bdry_traction_pressure_part_;  //!< boundary traction rhs on immersed cell (only
                                            //!< pressure part)
    //@}

    //! @name Pointers to Discretizations
    //@{
    Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;
    Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;
    Teuchos::RCP<DRT::Discretization> immerseddis_;
    Teuchos::RCP<DRT::Discretization> scatradis_;
    //@}

    /// pointer to cell structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed> cellstructure_;

    /// pointer to poro structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FPSIStructureWrapper> porostructure_;

    /// pointer to cell subproblem (structure-scatra interaction)
    Teuchos::RCP<SSI::SSI_Base> cell_subproblem_;

    /// pointer to poroelast subproblem (poro-scatra interaction)
    Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

    /// scatra adapter for poro-scatra time integration
    Teuchos::RCP<ADAPTER::AdapterScatraWrapperCellMigration> poroscatra_wrapper_;

    //! @name  basic information for parallelism
    int myrank_;
    int numproc_;
    //@}

    /// pointer to global problem
    DRT::Problem* globalproblem_;

    //! @name various bools and switches
    //@{
    bool displacementcoupling_;  //!< true if increment of displacement is chosen for convergence
                                 //!< check of partitioned scheme.
    bool multicellmigration_;    //!< true if more than one cell is present. Each cell needs to be
                                 //!< labeled with "ImmersedSearchbox" condition.
    bool isPseudo2D_;            //!< true if problem is modeled as a quasi 2D problem
    bool isALE_;                 //!< true if background mesh is suppose to be deformable
    //@}

    //! @name proteolytic migration control parameters
    //@{
    int migrationtype_;  //!< amoeboid or proteolytic migration
    double segregationconstant_;
    int segregationtype_;
    int segregationby_;
    //@}

    /// counter for number of unconverged timesteps
    int continued_steps_;

    /// number of integration points related to mass conserving adaption of equations
    int degree_gp_fluid_bound_;  //!< number of integration points in fluid elements cut by
                                 //!< structural boundary

    //! @name Manager for data exchange between fields
    //@{
    DRT::ImmersedFieldExchangeManager* exchange_manager_;
    //@}

    //! @name Search Trees
    //@{
    Teuchos::RCP<GEO::SearchTree> fluid_SearchTree_;  //!< search tree for fluid domain
    Teuchos::RCP<GEO::SearchTree> cell_SearchTree_;   //!< search tree for cell domain
    //@}

    //! @name current nodal positions of the discretizations
    //@{
    std::map<int, LINALG::Matrix<3, 1>>*
        currpositions_cell_;  //!< pointer to map of vectors for search tree containing current
                              //!< structural positions
    std::map<int, LINALG::Matrix<3, 1>>*
        currpositions_ECM_;  //!< pointer to map of vectors for search tree containing current fluid
                             //!< positions
    //@}

    //! @name current subsets of discretizations
    //@{
    std::map<int, std::set<int>> curr_subset_of_backgrounddis_;  //!< backgrd elements to evaluate
                                                                 //!< the dirichlet interpolation
    //@}

    //! @name Maps for Dirichlet BCs
    //@{
    Teuchos::RCP<Epetra_Map> dbcmap_immersed_;  //!< dirichlet bc map of immersed values
    Teuchos::RCP<Epetra_Map>
        dbcmap_immersed_scatra_;  //!< dirichlet bc map of immersed scatra values
                                  //@}

   private:
    /// improve quality of solution near the interface
    virtual void CorrectInterfaceVelocity();

    /// reset immersed information in fluid dis
    /// e.g. isimmersed_, isboundaryimmersed_
    void ResetImmersedInformation();

    //! @name various validation switches
    //@{
    bool artificial_velocity_isvalid_;  //!< this flag is true if the current vel. in artificial
                                        //!< fluid domain is consistent with current state of
                                        //!< structure
    bool boundary_traction_isvalid_;    //!< this flag is true if the current bdry. traction on the
                                        //!< structure is consistent with current state of fluid
    bool immersed_info_isvalid_;  //!< true if struct. configuration and fluid configuration are
                                  //!< consistently reflected in immersed infos, like e.g.
                                  //!< isimmersed_
    bool is_relaxation_;          //!< true if relaxation is used
    int correct_boundary_velocities_;  //!< true if interface accuracy correction is requested
    //@}

  };  // class ImmersedPartitionedFlowCellInteraction
}  // namespace IMMERSED

#endif
