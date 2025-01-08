// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
    explicit DirichletNeumannVolCoupl(MPI_Comm comm);

   public:
    /// setup this object
    void setup() override;

   protected:
    /// setup
    void setup_coupling_struct_ale(const Teuchos::ParameterList& fsidyn, MPI_Comm comm);

    /// setup
    void setup_interface_corrector(const Teuchos::ParameterList& fsidyn, MPI_Comm comm);

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
    std::shared_ptr<Core::LinAlg::Vector<double>> fluid_op(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp, const FillType fillFlag) final;


    void extract_previous_interface_solution() override;

    /// structure to ale mapping
    std::shared_ptr<Core::LinAlg::Vector<double>> stucture_to_ale(
        std::shared_ptr<Core::LinAlg::Vector<double>> iv) const;

    /// structure to ale mapping
    std::shared_ptr<Core::LinAlg::Vector<double>> structure_to_ale(
        const Core::LinAlg::Vector<double>& iv) const;

    /// ale to structure mapping
    std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_structure(
        Core::LinAlg::Vector<double>& iv) const;

    /// ale to structure
    std::shared_ptr<Core::LinAlg::Vector<double>> ale_to_structure(
        std::shared_ptr<const Core::LinAlg::Vector<double>> iv) const;

    /// coupling of structure and ale at the interface
    std::shared_ptr<Coupling::Adapter::MortarVolCoupl> coupsa_;

    /// coupling of structure and ale at the interface
    std::shared_ptr<InterfaceCorrector> icorrector_;
  };

  class VolCorrector;

  class InterfaceCorrector
  {
   public:
    /// constructor
    InterfaceCorrector() : idisp_(nullptr) {};

    /// destructor
    virtual ~InterfaceCorrector() = default;

    virtual void setup(std::shared_ptr<Adapter::FluidAle> fluidale);

    void set_interface_displacements(std::shared_ptr<Core::LinAlg::Vector<double>>& idisp_struct,
        Coupling::Adapter::Coupling& icoupfs);

    virtual void correct_interface_displacements(
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp_fluid,
        std::shared_ptr<FLD::Utils::MapExtractor> const& finterface);

   private:
    std::shared_ptr<const Core::LinAlg::Vector<double>> idisp_;
    std::shared_ptr<Coupling::Adapter::Coupling> icoupfs_;

    std::shared_ptr<Core::LinAlg::Vector<double>> deltadisp_;
    std::shared_ptr<Adapter::FluidAle> fluidale_;

    std::shared_ptr<VolCorrector> volcorrector_;
  };


  class VolCorrector
  {
   public:
    /// constructor
    VolCorrector() : dim_(-1) {};

    /// destructor
    virtual ~VolCorrector() = default;

    virtual void setup(const int dim, std::shared_ptr<Adapter::FluidAle> fluidale);

    virtual void correct_vol_displacements(std::shared_ptr<Adapter::FluidAle> fluidale,
        std::shared_ptr<Core::LinAlg::Vector<double>> deltadisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp_fluid,
        std::shared_ptr<FLD::Utils::MapExtractor> const& finterface);

   private:
    virtual void correct_vol_displacements_para_space(std::shared_ptr<Adapter::FluidAle> fluidale,
        std::shared_ptr<Core::LinAlg::Vector<double>> deltadisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp_fluid,
        std::shared_ptr<FLD::Utils::MapExtractor> const& finterface);

    virtual void correct_vol_displacements_phys_space(std::shared_ptr<Adapter::FluidAle> fluidale,
        std::shared_ptr<Core::LinAlg::Vector<double>> deltadisp,
        std::shared_ptr<Core::LinAlg::Vector<double>> idisp_fluid,
        std::shared_ptr<FLD::Utils::MapExtractor> const& finterface);

    void init_dop_normals();

    std::map<int, Core::LinAlg::Matrix<9, 2>> calc_background_dops(
        Core::FE::Discretization& searchdis);

    Core::LinAlg::Matrix<9, 2> calc_dop(Core::Elements::Element& ele);

    std::vector<int> search(
        Core::Elements::Element& ele, std::map<int, Core::LinAlg::Matrix<9, 2>>& currentKDOPs);

    //! Spatial dimension of the problem
    int dim_;

    //! Searchtree for mortar evaluations
    std::shared_ptr<Core::Geo::SearchTree> search_tree_;

    //! Dop normals for search algorithm
    Core::LinAlg::Matrix<9, 3> dopnormals_;

    std::map<int, std::vector<int>> fluidaleelemap_;

    std::map<int, std::vector<int>> fluidalenodemap_;
    std::map<int, std::vector<int>> fluidalenode_fs_imap_;
  };


}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
