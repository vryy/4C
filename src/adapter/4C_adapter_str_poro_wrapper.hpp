// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_PORO_WRAPPER_HPP
#define FOUR_C_ADAPTER_STR_PORO_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_field_wrapper.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Global
{
  class Problem;
}

namespace PoroElast
{
  class Monolithic;
}

namespace Adapter
{
  class FluidPoro;

  /// Just wrap, do nothing new, provides methods which are not available for Base Class
  /// Adapter::Field!
  class StructurePoroWrapper : public FieldWrapper
  {
   public:
    /// constructor
    explicit StructurePoroWrapper(
        std::shared_ptr<Field> field, FieldWrapper::Fieldtype type, bool NOXCorrection = false);

    /// setup
    void setup();

    /// communication object at the interface
    virtual std::shared_ptr<const Solid::MapExtractor> interface() const
    {
      return structure_->interface();
    }

    /// direct access to discretization
    virtual std::shared_ptr<Core::FE::Discretization> discretization()
    {
      return structure_->discretization();
    }

    /// return time integration factor
    virtual double tim_int_param() const { return structure_->tim_int_param(); }

    /// Access to output object
    virtual std::shared_ptr<Core::IO::DiscretizationWriter> disc_writer()
    {
      return structure_->disc_writer();
    }

    /// unknown displacements at \f$t_{n+1}\f$
    virtual std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp() const
    {
      return structure_->dispnp();
    }

    /// unknown displacements at \f$t_{n+1}\f$
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> write_access_dispnp() const
    {
      return structure_->write_access_dispnp();
    }

    /// get constraint manager defined in the structure
    virtual std::shared_ptr<CONSTRAINTS::ConstrManager> get_constraint_manager()
    {
      return structure_->get_constraint_manager();
    }

    // access to contact/meshtying bridge
    virtual std::shared_ptr<CONTACT::MeshtyingContactBridge> meshtying_contact_bridge()
    {
      return structure_->meshtying_contact_bridge();
    }

    /// extract interface displacements at \f$t_{n}\f$
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_dispn()
    {
      return structure_->extract_interface_dispn();
    }

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_interface_dispnp()
    {
      return structure_->extract_interface_dispnp();
    }

    /// extract interface displacements at \f$t_{n}\f$
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_fpsi_interface_dispn()
    {
      return structure_->interface()->extract_fpsi_cond_vector(*structure_->dispn());
    }

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_fpsi_interface_dispnp()
    {
      return structure_->interface()->extract_fpsi_cond_vector(*structure_->dispnp());
    }

    //! unique map of all dofs that should be constrained with DBC
    virtual std::shared_ptr<const Epetra_Map> combined_dbc_map();

    //! perform result test
    void test_results(Global::Problem* problem);

    //! return poro poro_field
    const std::shared_ptr<PoroElast::Monolithic>& poro_field();

    //! return poro structure_field
    const std::shared_ptr<FSIStructureWrapper>& structure_field();

    //! return poro fluid_field
    const std::shared_ptr<Adapter::FluidPoro>& fluid_field();

    //! Insert FSI Condition Vector
    std::shared_ptr<Core::LinAlg::Vector<double>> insert_fsi_cond_vector(
        const Core::LinAlg::Vector<double>& cond);

    //! Recover Lagrange Multiplier during iteration (does nothing for structure)
    void recover_lagrange_multiplier_after_newton_step(
        std::shared_ptr<Core::LinAlg::Vector<double>> iterinc);

    bool is_poro() { return (type_ == FieldWrapper::type_PoroField); }

   protected:
    std::shared_ptr<PoroElast::Monolithic> poro_;     ///< underlying poro time integration
    std::shared_ptr<FSIStructureWrapper> structure_;  ///< underlying structural time integration
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
