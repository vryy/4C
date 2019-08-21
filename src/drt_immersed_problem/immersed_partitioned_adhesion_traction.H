/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed cell-ecm interaction via adhesion traction

\level 2

\maintainer Jonas Eichinger

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
  // class Monolithic;
  class PoroScatraBase;
}  // namespace POROELAST

namespace ADAPTER
{
  class FPSIStructureWrapper;
  class FSIStructureWrapperImmersed;
}  // namespace ADAPTER

namespace SSI
{
  class SSI_Part2WC;
}

/*----------------------------------------------------------------------*/

namespace IMMERSED
{
  class ImmersedPartitionedAdhesionTraction : public ImmersedPartitioned
  {
   public:
    /// setup partitioned immersed cell migration algorithm
    explicit ImmersedPartitionedAdhesionTraction(
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

    /*!
    \brief Update the set of coupled cell-ECM nodes.

     Updates the coupled ECM fixpoints depending on whether there are bound surface proteins or not.

    \param none
    \note none
    \warning none
    \return void
    */
    void UpdateAdhesionInformation();


    /// print information on the current solution step
    virtual void PrintStepInfo();


    /*!
    \brief Distribute force from cell adhesion node inside ecm element to ecm element nodes.

     Description: Fills members \ref adh_nod_param_coords_in_backgrd_ele_ and \ref
    adh_nod_backgrd_ele_mapping_ . Performs the force spreading from the cell adhesion nodes onto
    the adjacent ecm nodes.

     Parallelism: The adhesion nodes may be on different processors than the corresponding
    background elements. We therefore have to perform a round-robin loop which tries to match the
    adhesion nodes belonging to each processor to a background element on any proc.

     1) Loop over all nodes with 'CellFocalAdhesion' condition.
       -> 1.1 get displacement of current node.

       -> 1.2 Loop over curr_subset_of_backgrounddis_ and search only in elements for which
    IsBoundaryImmersed() == true.
          -> 1.2.1 find backgrd. ECM element in which current node is immersed.
          -> 1.2.2 get parameter space coord. of current node in backgrd. element.
          -> 1.2.3 write pair 'cell adhesion node id' <-> 'xi' in
    adh_nod_param_coords_in_backgrd_ele_.
          -> 1.2.4 write pair 'adhesion node id' <-> 'backgrdele id' in map
    adh_nod_backgrd_ele_mapping_.
          -> 1.2.5 spread reaction force at current node to nodes of backgrd. element by the ratio
    of shapefcts. value at current 'xi' and add those contributions to ecm_adhesion_forces_.

    \param adhesionforce  (out) : Force at adhesion nodes of immersed cell discretization

    \note Requires immersed cell discretization to have 'CellFocalAdhesion' defined in the input
          file (via '--DESIGN FOCAL ADHESION SURF CONDITIONS').

    \warning This method uses member variables and cannot be used stand-alone.

    \return void
    */
    void DistributeAdhesionForce(const Teuchos::RCP<const Epetra_Vector>& force_to_distribute);

    /*!
    \brief Interpolate backgrd. ecm element displacements to immersed cell adhesion node inside
    element.

     Requires immersed cell discretization to have 'CellFocalAdhesion' defined in the input
     file (via '--DESIGN FOCAL ADHESION SURF CONDITIONS'). Requires the same surface/s to
     have also defined Dirichlet conditions to be able to retrieve the reaction forces from
     the immersed cell.

     1) Loop over all nodes in adh_nod_param_coords_in_backgrd_ele_.
       -> 1.1 extract background ECM element nodal displacements.
       -> 1.2 interpolate ECM displacement to parameter space coordinate occupied by the cell
    adhesion node prior to the ECM deformation. Thereby the adhesion node is supposed to be fixed to
    the same material point.
       -> 1.3 save the interpolated displacements in cell_adhesion_disp_.

    \note this method uses member variables and cannot be used stand-alone.

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

    /// adhesion coupling by Dirichlet or penalty
    bool dirichletcoupling_;

    /// concentration defining reset concentration for adhesion sites
    double reset_conc_;

    //! @name Various global vectors
    //@{
    Teuchos::RCP<Epetra_Vector> ecm_adhesion_forces_;   //!< current nodal adhesion forces on ecm
    Teuchos::RCP<Epetra_Vector> cell_adhesion_forces_;  //!< current nodal adhesion forces on cell
    Teuchos::RCP<Epetra_Vector> cell_adhesion_disp_;  //!< current nodal cell adhesion displacements
    Teuchos::RCP<Epetra_Vector>
        cell_adhesion_nod_coords_;  //!< current nodal cell adhesion coordinates
    Teuchos::RCP<Epetra_Vector>
        Freact_cell_;  //!< current nodal reaction force because of dirichlet at adhesions
    //@}

    //! @name adhesion dynamics related map to save parameter space coords. of adhesion nodes
    //@{
    std::map<int, LINALG::Matrix<3, 1>>
        adh_nod_param_coords_in_backgrd_ele_;  //!< map adhesion node id <-> xi in backgrd. ele
    std::map<int, int> adh_nod_backgrd_ele_mapping_;  //!< map adhesion node id <-> backgrdele id
    //@}

    //! @name Pointers to Discretizations
    //@{
    Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;
    Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;
    Teuchos::RCP<DRT::Discretization> immerseddis_;
    Teuchos::RCP<DRT::Discretization> immersedscatradis_;
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

    /// method for treatment of artificial ecm solid phase
    int artificial_ecm_treatment_;

    /// do we need an update of our immersed information (e.g. isimmersed_ in nodes)
    bool immersed_information_invalid_;

    /// write output after every fixed-point iteration
    bool output_evry_nlniter_;

    /// switch deciding whether adhesion state has to be (re-)determined
    bool determine_adhesion_information_;

    /// cell structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed> cellstructure_;


    /// pointer to poroelast subproblem (poro-scatra interaction)
    Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

    /// poro structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FPSIStructureWrapper> porostructure_;

    /// pointer to cellscatra subproblem (structure-scatra interaction)
    Teuchos::RCP<SSI::SSI_Part2WC> cellscatra_subproblem_;

    /// reinitialize the transfer vectors
    void ReinitTransferVectors()
    {
      ecm_adhesion_forces_->Scale(0.0);
      if (dirichletcoupling_) cell_adhesion_disp_->Scale(0.0);
      return;
    };

    /// get immersed nodes and determine their dofs
    void BuildImmersedDirichMap(const Teuchos::RCP<const DRT::Discretization>& dis,
        const Teuchos::RCP<const Epetra_Map>& dirichmap_orig,
        Teuchos::RCP<Epetra_Map>& dirichmap_new, std::map<int, int>& adh_nod_backgrd_ele_mapping);

    /// add immersed dirichlet values from immersed dis to systemvector of background dis
    void DoImmersedDirichletCond(const Teuchos::RCP<Epetra_Vector>& statevector,
        const Teuchos::RCP<Epetra_Vector>& dirichvals,
        const Teuchos::RCP<const Epetra_Map>& dbcmap);

    /// reset immersed information like e.g.: IsImmersed()
    void ResetImmersedInformation();

  };  // class ImmersedPartitionedAdhesionTraction
}  // namespace IMMERSED

#endif
