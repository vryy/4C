/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
       with volume coupling

\level 3

*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_DIRICHLETNEUMANN_VOLCOUPL_HPP
#define FOUR_C_FSI_DIRICHLETNEUMANN_VOLCOUPL_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fsi_dirichletneumann_disp.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEO
{
  class SearchTree;
}

namespace FSI
{
  class InterfaceCorrector;

  /// Dirichlet-Neumann VolCoupled system
  class DirichletNeumannVolCoupl : public DirichletNeumannDisp
  {
    friend class DirichletNeumannFactory;

   protected:
    /**
     *  \brief constructor
     *
     * You will have to use the FSI::DirichletNeumannFactory to create an instance of this class
     */
    explicit DirichletNeumannVolCoupl(const Epetra_Comm& comm);

   public:
    /// setup this object
    void Setup() override;

   protected:
    /// setup
    void setup_coupling_struct_ale(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm);

    /// setup
    void setup_interface_corrector(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm);

    /** \brief interface fluid operator
     *
     * Solve the fluid field problem.  Since the fluid field is the Dirichlet partition, the
     * interface displacement is prescribed as a Dirichlet boundary condition.
     *
     * \param[in] idisp interface displacement
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns interface force
     */
    Teuchos::RCP<Epetra_Vector> fluid_op(
        Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag) final;


    void extract_previous_interface_solution() override;

    /// structure to ale mapping
    Teuchos::RCP<Epetra_Vector> stucture_to_ale(Teuchos::RCP<Epetra_Vector> iv) const;

    /// structure to ale mapping
    Teuchos::RCP<Epetra_Vector> structure_to_ale(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// ale to structure mapping
    Teuchos::RCP<Epetra_Vector> ale_to_structure(Teuchos::RCP<Epetra_Vector> iv) const;

    /// ale to structure
    Teuchos::RCP<Epetra_Vector> ale_to_structure(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> coupsa_;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<InterfaceCorrector> icorrector_;
  };

  class VolCorrector;

  class InterfaceCorrector
  {
   public:
    /// constructor
    InterfaceCorrector() : idisp_(Teuchos::null){};

    /// destructor
    virtual ~InterfaceCorrector() = default;

    virtual void Setup(Teuchos::RCP<ADAPTER::FluidAle> fluidale);

    void set_interface_displacements(
        Teuchos::RCP<Epetra_Vector>& idisp_struct, CORE::ADAPTER::Coupling& icoupfs);

    virtual void correct_interface_displacements(Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

   private:
    Teuchos::RCP<const Epetra_Vector> idisp_;
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoupfs_;

    Teuchos::RCP<Epetra_Vector> deltadisp_;
    Teuchos::RCP<ADAPTER::FluidAle> fluidale_;

    Teuchos::RCP<VolCorrector> volcorrector_;
  };


  class VolCorrector
  {
   public:
    /// constructor
    VolCorrector() : dim_(-1){};

    /// destructor
    virtual ~VolCorrector() = default;

    virtual void Setup(const int dim, Teuchos::RCP<ADAPTER::FluidAle> fluidale);

    virtual void correct_vol_displacements(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
        Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

   private:
    virtual void correct_vol_displacements_para_space(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
        Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

    virtual void correct_vol_displacements_phys_space(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
        Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

    void init_dop_normals();

    std::map<int, CORE::LINALG::Matrix<9, 2>> calc_background_dops(
        Teuchos::RCP<DRT::Discretization> searchdis);

    CORE::LINALG::Matrix<9, 2> calc_dop(DRT::Element& ele);

    std::vector<int> search(
        DRT::Element& ele, std::map<int, CORE::LINALG::Matrix<9, 2>>& currentKDOPs);

    //! Spatial dimension of the problem
    int dim_;

    //! Searchtree for mortar evaluations
    Teuchos::RCP<CORE::GEO::SearchTree> search_tree_;

    //! Dop normals for search algorithm
    CORE::LINALG::Matrix<9, 3> dopnormals_;

    std::map<int, std::vector<int>> fluidaleelemap_;

    std::map<int, std::vector<int>> fluidalenodemap_;
    std::map<int, std::vector<int>> fluidalenode_fs_imap_;
  };


}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
