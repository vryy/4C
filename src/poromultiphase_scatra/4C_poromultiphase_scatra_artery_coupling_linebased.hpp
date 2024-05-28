/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for line-based (non-conforming) coupling between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_LINEBASED_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_ARTERY_COUPLING_LINEBASED_HPP

#include "4C_config.hpp"

#include "4C_inpar_bio.hpp"
#include "4C_poromultiphase_scatra_artery_coupling_nonconforming.hpp"

// forward declaration
class Epetra_IntVector;

FOUR_C_NAMESPACE_OPEN

namespace POROMULTIPHASESCATRA
{
  // forward declaration
  class PoroMultiPhaseScatraArteryCouplingPairBase;

  //! Line based coupling between artery network and poromultiphasescatra algorithm
  class PoroMultiPhaseScaTraArtCouplLineBased : public PoroMultiPhaseScaTraArtCouplNonConforming
  {
   public:
    //! create using a Epetra_Comm
    PoroMultiPhaseScaTraArtCouplLineBased(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname);

    //! set-up linear system of equations of coupled problem
    void SetupSystem(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs, Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art,
        Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_cont,
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_art) override;

    //! setup the strategy
    void Setup() override;

    //! apply mesh movement (on artery elements)
    void ApplyMeshMovement() override;

    //! access to blood vessel volume fraction
    Teuchos::RCP<const Epetra_Vector> blood_vessel_volume_fraction() override;

   private:
    //! pre-evaluate the pairs and sort out duplicates
    void pre_evaluate_coupling_pairs();

    //! fill the length not changed through deformation and initialize curr length
    void fill_unaffected_artery_length();

    //! fill the integrated diameter not changed through varying blood vessel diameter
    void fill_unaffected_integrated_diam();

    //! calculate the volume fraction occupied by blood vessels
    void calculate_blood_vessel_volume_fraction();

    //! create the GID to segment vector
    void create_gid_to_segment_vector();

    //! fill the GID to segment vector
    void fill_gid_to_segment_vector(
        const std::vector<Teuchos::RCP<
            POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
        std::map<int, std::vector<double>>& gid_to_seglength);

    //! set the artery diameter in column based vector
    void fill_artery_ele_diam_col();

    //! (re-)set the artery diameter in material to be able to use it on 1D discretization
    void set_artery_diam_in_material() override;

    //! reset the integrated diameter vector to zero
    void reset_integrated_diam_to_zero() override;

    /*!
     * @brief Utility function for depth-first search for the connected components of the 1D artery
     * discretization
     *
     * @param actnode : currently checked node
     * @param visited : vector that marks visited nodes
     * @param artsearchdis : artery-discretization in fully-overlapping format
     * @param ele_diams_artery_full_overlap : vector of element diameters in fully-overlapping
     * format
     * @param this_connected_comp : current connected component
     */
    void depth_first_search_util(DRT::Node* actnode, Teuchos::RCP<Epetra_IntVector> visited,
        Teuchos::RCP<DRT::Discretization> artconncompdis,
        Teuchos::RCP<const Epetra_Vector> ele_diams_artery_full_overlap,
        std::vector<int>& this_connected_comp);

    /*!
     * @brief find free-hanging 1D elements
     *
     * Find the free hanging 1D elements which have to be taken out of simulation during blood
     * vessel collapse. This is realized by getting the connected components of the 1D graph. If no
     * node of a connected component has a DBC or if it is smaller than a user-specified
     * threshold, all its elements are taken out
     * @param : eles_to_be_deleted vector of free-hanging elements
     */
    void find_free_hanging1_d_elements(std::vector<int>& eles_to_be_deleted);

    //! evaluate additional linearization of (integrated) element diameter dependent terms
    //! (Hagen-Poiseuille)
    void evaluate_additional_linearizationof_integrated_diam() override;

    /*!
     * @brief apply additional Dirichlet boundary condition for collapsed 1D elements to avoid
     * singular stiffness matrix
     *
     * apply additional dirichlet boundary conditions of zero pressure or mass fraction on nodes
     * which only border collapsed 1D elements, i.e., free-hanging nodes with zero row in global
     * stiffness matrix to avoid singularity of this matrix
     * \note this procedure is equivalent to taking collapsed elements out of the simulation
     * entirely
     *
     * @param[in]        dbcmap_art map of nodes with DBC of 1D discretization
     * @param[in,out]   rhs_art_with_collapsed right hand side of artery subpart
     * @returns dbcmap, also containing additional boundary condition for collapsed eles
     */
    Teuchos::RCP<Epetra_Map> get_additional_dbc_for_collapsed_eles(
        Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_art,
        Teuchos::RCP<Epetra_Vector> rhs_art_with_collapsed);

    //! FE-assemble into global force and stiffness
    void fe_assemble_ele_force_stiff_into_system_vector_matrix(const int& ele1gid,
        const int& ele2gid, const double& integrated_diam,
        std::vector<CORE::LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector> rhs) override;

    //! get the segment lengths of element 'artelegid'
    std::vector<double> get_ele_segment_lengths(const int artelegid) override;

    //! check for duplicate segment
    bool is_duplicate_segment(
        const std::vector<Teuchos::RCP<
            POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
        const Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
            possible_duplicate);

    //! check for identical segment
    bool is_identical_segment(
        const std::vector<Teuchos::RCP<
            POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
        const int& ele1gid, const double& etaA, const double& etaB, int& elepairID);

    //! set flag if varying diameter has to be calculated
    void set_varying_diam_flag() override;

    //! print output of mesh tying pairs
    void output_summary() const;

    //! print out the coupling method
    void print_out_coupling_method() const override;

    //! maximum number of segments per artery element
    int maxnumsegperartele_;

    //! length of artery elements unaffected by deformation
    Teuchos::RCP<Epetra_FEVector> unaffected_seg_lengths_artery_;

    //! length of artery elements in current configuration
    Teuchos::RCP<Epetra_FEVector> current_seg_lengths_artery_;

    //! diameter of the artery element integrated over the length of the artery element (row format
    //! and FE vector due to non-local assembly)
    Teuchos::RCP<Epetra_FEVector> integrated_diams_artery_row_;

    //! diameter of artery element integrated over the length of the artery element (col format)
    Teuchos::RCP<Epetra_Vector> integrated_diams_artery_col_;

    //! diameter of artery element (col format)
    Teuchos::RCP<Epetra_Vector> ele_diams_artery_col_;

    //! unaffected diameter integrated over the length of the artery element
    //! (protruding elements for which diameter does not change)
    Teuchos::RCP<Epetra_Vector> unaffected_integrated_diams_artery_col_;

    //! volume fraction of blood vessels (for output)
    Teuchos::RCP<Epetra_Vector> bloodvesselvolfrac_;

    //! gid to segment: stores [GID; [eta_a eta_b]_1, [eta_a eta_b]_2, ..., [eta_a eta_b]_n]
    //!  of artery elements in column format, i.e. fully overlapping
    std::map<int, std::vector<double>> gid_to_segment_;

    //! gid to segment length: stores [GID; seglength_1, seglength_2, ..., seglength_n]
    //!  of artery elements in column format, i.e. fully overlapping (only used for
    //!  porofluid-problems)
    std::map<int, std::vector<double>> gid_to_seglength_;
  };
}  // namespace POROMULTIPHASESCATRA

FOUR_C_NAMESPACE_CLOSE

#endif
