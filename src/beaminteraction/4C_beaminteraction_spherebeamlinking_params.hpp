// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_SPHEREBEAMLINKING_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_SPHEREBEAMLINKING_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_inpar_beaminteraction.hpp"
#include "4C_utils_exceptions.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN



// forward declaration

namespace Solid
{
  namespace TimeInt
  {
    class BaseDataGlobalState;
  }
}  // namespace Solid
namespace Mat
{
  class CrosslinkerMat;
}
namespace BeamInteraction
{
  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction */
  class SphereBeamLinkingParams
  {
   public:
    //! constructor
    SphereBeamLinkingParams();

    //! destructor
    virtual ~SphereBeamLinkingParams() = default;

    //! initialize with the stuff coming from input file
    void init(Solid::TimeInt::BaseDataGlobalState const& gstate);

    //! setup member variables
    void setup();

    //! reset time step in case structure time is adapted during simulation time
    void reset_time_step(double structure_delta_time);

    //! returns the isinit_ flag
    inline const bool& is_init() const { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& is_setup() const { return issetup_; };

    //! Checks the init and setup status
    inline void check_init_setup() const
    {
      if (!is_init() or !is_setup()) FOUR_C_THROW("Call init() and setup() first!");
    }

    //! Checks the init status
    inline void check_init() const
    {
      if (!is_init()) FOUR_C_THROW("init() has not been called, yet!");
    }

    /// linker material id
    std::shared_ptr<Mat::CrosslinkerMat> get_linker_material() const
    {
      /// HACK: FIX IF MORE THAN ONE CROSSLINKER TYPE
      check_init_setup();
      return mat_.back();
    };

    /// time step for stochastic events concerning crosslinking
    double const& delta_time() const
    {
      check_init_setup();
      return deltatime_;
    };

    /// contraction rate of cell (integrin linker) in [microm/s]
    double contraction_rate(Inpar::BeamInteraction::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return contractionrate_.at(linkertype);
    };

    /// number of linker per type
    std::vector<int> const& max_num_linker_per_type() const
    {
      check_init_setup();
      return maxnumlinkerpertype_;
    };

    /// material number for linker types
    std::vector<int> const& mat_linker_per_type() const
    {
      check_init_setup();
      return matlinkerpertype_;
    };

    /// get all active linker types
    std::vector<Inpar::BeamInteraction::CrosslinkerType> const& linker_types() const
    {
      check_init_setup();
      return linkertypes_;
    };

    // distance between two binding spots on a filament
    double filament_bspot_interval_global(Inpar::BeamInteraction::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotintervalglobal_.at(linkertype);
    };

    // distance between two binding spots on a filament
    double filament_bspot_interval_local(Inpar::BeamInteraction::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotintervallocal_.at(linkertype);
    };

    // start and end arc parameter for binding spots on a filament
    std::pair<double, double> const& filament_bspot_range_local(
        Inpar::BeamInteraction::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotrangelocal_.at(linkertype);
    };

    // start and end arc parameter for binding spots on a filament
    std::pair<double, double> const& filament_bspot_range_global(
        Inpar::BeamInteraction::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotrangeglobal_.at(linkertype);
    };

   private:
    bool isinit_;

    bool issetup_;

    /// time step for stochastic events concerning integrins, e.g. catch-slip-bond behavior
    double deltatime_;
    bool own_deltatime_;
    /// contraction rate of cell (integrin linker) in [microm/s]
    std::map<Inpar::BeamInteraction::CrosslinkerType, double> contractionrate_;
    /// crosslinker material
    std::vector<std::shared_ptr<Mat::CrosslinkerMat>> mat_;
    /// number of crosslinkers in the simulated volume
    std::vector<int> maxnumlinkerpertype_;
    /// material numbers for crosslinker types
    std::vector<int> matlinkerpertype_;
    /// linker and therefore binding spot types
    std::vector<Inpar::BeamInteraction::CrosslinkerType> linkertypes_;
    /// distance between two binding spots on each filament
    std::map<Inpar::BeamInteraction::CrosslinkerType, double> filamentbspotintervalglobal_;
    /// distance between two binding spots on a filament as percentage of filament reference length
    std::map<Inpar::BeamInteraction::CrosslinkerType, double> filamentbspotintervallocal_;
    /// start and end arc parameter for binding spots on a filament
    std::map<Inpar::BeamInteraction::CrosslinkerType, std::pair<double, double>>
        filamentbspotrangeglobal_;
    /// start and end arc parameter for binding spots on a filament
    /// in percent of filament reference length
    std::map<Inpar::BeamInteraction::CrosslinkerType, std::pair<double, double>>
        filamentbspotrangelocal_;
  };

}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
