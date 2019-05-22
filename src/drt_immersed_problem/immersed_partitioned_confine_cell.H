/*!----------------------------------------------------------------------

\brief partitioned immersed cell-ecm interaction via cell confinement in ecm pore

\level 3

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_CONFINE_CELL_H
#define IMMERSED_PARTITIONED_CONFINE_CELL_H

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
  class FPSIStructureWrapper;
  class FSIStructureWrapperImmersed;
}  // namespace ADAPTER

/*----------------------------------------------------------------------*/

namespace IMMERSED
{
  class ImmersedPartitionedConfineCell : public ImmersedPartitioned
  {
   public:
    /// setup partitioned immersed cell migration algorithm
    explicit ImmersedPartitionedConfineCell(
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

    /// read restart
    virtual void ReadRestart(int step);

   protected:
    /// composed Field-Field-Interaction operator
    void CouplingOp(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

    /// background operator
    void BackgroundOp(
        Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values, const FillType fillFlag);

    /// immersed operator
    Teuchos::RCP<Epetra_Vector> ImmersedOp(
        Teuchos::RCP<Epetra_Vector> bdry_traction, const FillType fillFlag);

    /// preparation of background field solve
    void PrepareBackgroundOp();

    /// preparation of immersed field solve
    void PrepareImmersedOp();

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


   private:
    /// basic information for parallelism
    int myrank_;
    int numproc_;

    /// global timstep
    double timestep_;

    /// pointer to global problem
    DRT::Problem* globalproblem_;

    /// immersed coupling variable
    bool displacementcoupling_;

    /// start penalty parameter at the beginning of each timestep
    double penalty_start_;

    //! @name Various global vectors
    //@{
    Teuchos::RCP<Epetra_Vector> penalty_gap_;  //!< gap integrated over surface
    Teuchos::RCP<Epetra_Vector>
        curr_nodal_normals_;  //!< current nodal normals included in cell output for visualization
    Teuchos::RCP<Epetra_Vector> cell_penalty_traction_;  //!< boundary traction rhs on immersed cell
                                                         //!< from cell-ecm interaction
    Teuchos::RCP<Epetra_Vector> cell_delta_penalty_traction_;  //!< force increment
    //@}

    //! @name adhesion dynamics related map to save parameter space coords. of adhesion nodes
    //@{
    std::map<int, LINALG::Matrix<3, 1>> adh_nod_param_coords_in_backgrd_ele_;
    std::map<int, int> adh_nod_backgrd_ele_mapping_;
    //@}

    //! @name Pointers to Discretizations
    //@{
    Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;
    Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;
    Teuchos::RCP<DRT::Discretization> immerseddis_;
    //@}

    //! @name cell initialization
    //@{
    bool initialize_cell_;
    double initial_timestep_;
    double penalty_init_;  //!< init penalty parameter
    int initialization_steps_;
    //@}

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
    std::map<int, std::set<int>>
        curr_subset_of_backgrounddis_;  //!< background elements covered by immersed discretisation
    //@}

    /// true if problem is modeled as a quasi 2D problem
    bool isPseudo2D_;

    /// cell structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed> cellstructure_;


    /// pointer to poroelast subproblem (poro-scatra interaction)
    Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

    /// poro structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FPSIStructureWrapper> porostructure_;

    /// reinitialize the transfer vectors
    void ReinitTransferVectors()
    {
      cell_delta_penalty_traction_->PutScalar(0.0);
      penalty_gap_->PutScalar(0.0);
      return;
    };

    /// calc the current tractions from the ecm to the cell due to confinement
    void CalcECMTractionOnStructure();

    /// reset element and node information about immersed method
    void ResetImmersedInformation();

    /// set time step in single fields
    virtual void SetFieldDt();

  };  // class ImmersedPartitionedConfineCell
}  // namespace IMMERSED

#endif /* IMMERSED_PARTITIONED_CONFINE_CELL_H */
