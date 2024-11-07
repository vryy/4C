// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_FSIWRAPPER_IMMERSED_HPP
#define FOUR_C_ADAPTER_STR_FSIWRAPPER_IMMERSED_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_structure_new_dbc.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Solid
{
  class Dbc;

  namespace Aux
  {
    class MapExtractor;
  }
}  // namespace Solid


namespace Adapter
{
  class FSIStructureWrapperImmersed : public FPSIStructureWrapper
  {
   public:
    /// constructor
    explicit FSIStructureWrapperImmersed(std::shared_ptr<Structure> structure);

    /// extract interface displacements at \f$t_{n+1}\f$ of immersed interface
    virtual std::shared_ptr<Core::LinAlg::Vector<double>> extract_immersed_interface_dispnp();

    /// Get mutable reference to DBC object
    Solid::Dbc& get_dbc();

    /// expand dirichlet bc map
    /// old struct. time integration version
    void add_dirich_dofs(const std::shared_ptr<const Epetra_Map> maptoadd) override;

    /// contract dirichlet bc map
    /// old struct. time integration version
    void remove_dirich_dofs(const std::shared_ptr<const Epetra_Map> maptoremove) override;

    /// set the state of the nox group and the global state data container
    void set_state(const std::shared_ptr<Core::LinAlg::Vector<double>>& x) override;

    /// @name Apply interface forces

    /// apply interface forces to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    void apply_immersed_interface_forces(std::shared_ptr<Core::LinAlg::Vector<double>> iforce_fsi,
        std::shared_ptr<Core::LinAlg::Vector<double>> iforce_immersed);

    /*!
      \brief Write extra output for specified step and time.
             Useful if you want to write output every iteration in partitioned schemes.
             If no step and time is provided, standard Output of structure field is invoked.

      \param forced_writerestart (in) : Force to write restart
      \param step (in) : Pseudo-step for which extra output is written
      \param time (in) : Pseudo-time for which extra output is written

      \note This is a pure DEBUG functionality. Originally used in immersed method development.

      \warning This method partly re-implements redundantly few lines of the common structure
      output() routine. \return void
    */
    virtual void output(
        bool forced_writerestart = false, const int step = -1, const double time = -1.0);

   protected:
    /// the interface map setup for immersed interface <-> fsi interface distinction
    std::shared_ptr<Core::LinAlg::MapExtractor> combinedinterface_;

    /// combined matching FSI - IMMERSED interface map
    std::shared_ptr<Epetra_Map> combinedmap_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
