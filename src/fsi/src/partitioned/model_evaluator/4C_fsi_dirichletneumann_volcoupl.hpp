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

namespace Core::Geo
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
    void setup() override;

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
    Teuchos::RCP<Core::LinAlg::Vector<double>> fluid_op(
        Teuchos::RCP<Core::LinAlg::Vector<double>> idisp, const FillType fillFlag) final;


    void extract_previous_interface_solution() override;

    /// structure to ale mapping
    Teuchos::RCP<Core::LinAlg::Vector<double>> stucture_to_ale(
        Teuchos::RCP<Core::LinAlg::Vector<double>> iv) const;

    /// structure to ale mapping
    Teuchos::RCP<Core::LinAlg::Vector<double>> structure_to_ale(
        const Core::LinAlg::Vector<double>& iv) const;

    /// ale to structure mapping
    Teuchos::RCP<Core::LinAlg::Vector<double>> ale_to_structure(
        Core::LinAlg::Vector<double>& iv) const;

    /// ale to structure
    Teuchos::RCP<Core::LinAlg::Vector<double>> ale_to_structure(
        Teuchos::RCP<const Core::LinAlg::Vector<double>> iv) const;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<Coupling::Adapter::MortarVolCoupl> coupsa_;

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

    virtual void setup(Teuchos::RCP<Adapter::FluidAle> fluidale);

    void set_interface_displacements(Teuchos::RCP<Core::LinAlg::Vector<double>>& idisp_struct,
        Coupling::Adapter::Coupling& icoupfs);

    virtual void correct_interface_displacements(
        Teuchos::RCP<Core::LinAlg::Vector<double>> idisp_fluid,
        Teuchos::RCP<FLD::Utils::MapExtractor> const& finterface);

   private:
    Teuchos::RCP<const Core::LinAlg::Vector<double>> idisp_;
    Teuchos::RCP<Coupling::Adapter::Coupling> icoupfs_;

    Teuchos::RCP<Core::LinAlg::Vector<double>> deltadisp_;
    Teuchos::RCP<Adapter::FluidAle> fluidale_;

    Teuchos::RCP<VolCorrector> volcorrector_;
  };


  class VolCorrector
  {
   public:
    /// constructor
    VolCorrector() : dim_(-1){};

    /// destructor
    virtual ~VolCorrector() = default;

    virtual void setup(const int dim, Teuchos::RCP<Adapter::FluidAle> fluidale);

    virtual void correct_vol_displacements(Teuchos::RCP<Adapter::FluidAle> fluidale,
        Teuchos::RCP<Core::LinAlg::Vector<double>> deltadisp,
        Teuchos::RCP<Core::LinAlg::Vector<double>> idisp_fluid,
        Teuchos::RCP<FLD::Utils::MapExtractor> const& finterface);

   private:
    virtual void correct_vol_displacements_para_space(Teuchos::RCP<Adapter::FluidAle> fluidale,
        Teuchos::RCP<Core::LinAlg::Vector<double>> deltadisp,
        Teuchos::RCP<Core::LinAlg::Vector<double>> idisp_fluid,
        Teuchos::RCP<FLD::Utils::MapExtractor> const& finterface);

    virtual void correct_vol_displacements_phys_space(Teuchos::RCP<Adapter::FluidAle> fluidale,
        Teuchos::RCP<Core::LinAlg::Vector<double>> deltadisp,
        Teuchos::RCP<Core::LinAlg::Vector<double>> idisp_fluid,
        Teuchos::RCP<FLD::Utils::MapExtractor> const& finterface);

    void init_dop_normals();

    std::map<int, Core::LinAlg::Matrix<9, 2>> calc_background_dops(
        Core::FE::Discretization& searchdis);

    Core::LinAlg::Matrix<9, 2> calc_dop(Core::Elements::Element& ele);

    std::vector<int> search(
        Core::Elements::Element& ele, std::map<int, Core::LinAlg::Matrix<9, 2>>& currentKDOPs);

    //! Spatial dimension of the problem
    int dim_;

    //! Searchtree for mortar evaluations
    Teuchos::RCP<Core::Geo::SearchTree> search_tree_;

    //! Dop normals for search algorithm
    Core::LinAlg::Matrix<9, 3> dopnormals_;

    std::map<int, std::vector<int>> fluidaleelemap_;

    std::map<int, std::vector<int>> fluidalenodemap_;
    std::map<int, std::vector<int>> fluidalenode_fs_imap_;
  };


}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
