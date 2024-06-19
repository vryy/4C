/*---------------------------------------------------------------------*/
/*! \file
\brief data container holding all crosslinking input parameters

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_CROSSLINKING_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_CROSSLINKING_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_inpar_beaminteraction.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace STR
{
  namespace TimeInt
  {
    class BaseDataGlobalState;
  }
}  // namespace STR
namespace BEAMINTERACTION
{
  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction
   * author eichinger*/
  class CrosslinkingParams
  {
   public:
    //! constructor
    CrosslinkingParams();

    //! destructor
    virtual ~CrosslinkingParams() = default;

    //! initialize with the stuff coming from input file
    void init(STR::TimeInt::BaseDataGlobalState const& gstate);

    //! setup member variables
    void setup();

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

    /// number of crosslinkers per type
    std::vector<int> const& num_crosslinker_per_type() const
    {
      check_init_setup();
      return numcrosslinkerpertype_;
    };

    /// number of crosslinkers per type
    int num_init_crosslinker_per_crosslinker_mat_id(int matid) const
    {
      check_init_setup();
      return maxnum_init_crosslinker_pertype_.at(matid);
    };

    /// number of crosslinkers per type
    int total_num_init_crosslinker() const
    {
      check_init_setup();
      int sum = 0;
      for (auto const& iter : maxnum_init_crosslinker_pertype_) sum += iter.second;
      return sum;
    };

    /// material number for crosslinker types
    std::vector<int> const& mat_crosslinker_per_type() const
    {
      check_init_setup();
      return matcrosslinkerpertype_;
    };

    /// get all active crosslinker types
    std::vector<Inpar::BEAMINTERACTION::CrosslinkerType> const& LinkerTypes() const
    {
      check_init_setup();
      return linkertypes_;
    };

    /// number of different crosslinker types in simulation volume
    int number_of_crosslinker_types() const
    {
      check_init_setup();
      return static_cast<int>(numcrosslinkerpertype_.size());
    };

    /// ~ 1e-3 / 2.27 according to cyron2011 eq 52 ff, viscosity of surrounding fluid
    double const& Viscosity() const
    {
      check_init_setup();
      return viscosity_;
    };

    /// thermal energy
    double const& KT() const
    {
      check_init_setup();
      return kt_;
    };

    /// time step for stochastic events concerning crosslinking
    double const& DeltaTime() const
    {
      check_init_setup();
      return deltatime_;
    };

    /// time step for stochastic events concerning crosslinking
    Core::LinAlg::Matrix<3, 2> const& linker_initialization_box() const
    {
      check_init_setup();
      return init_box_;
    };

    // distance between two binding spots on a filament
    int max_number_of_bonds_per_filament_bspot(
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return max_num_bonds_per_filament_bspot_.at(linkertype);
    };

    // distance between two binding spots on a filament
    double filament_bspot_interval_global(Inpar::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotintervalglobal_.at(linkertype);
    };

    // distance between two binding spots on a filament
    double filament_bspot_interval_local(Inpar::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotintervallocal_.at(linkertype);
    };

    // start and end arc parameter for binding spots on a filament
    std::pair<double, double> const& filament_bspot_range_local(
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotrangelocal_.at(linkertype);
    };

    // start and end arc parameter for binding spots on a filament
    std::pair<double, double> const& filament_bspot_range_global(
        Inpar::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return filamentbspotrangeglobal_.at(linkertype);
    };

   private:
    bool isinit_;

    bool issetup_;

    /// viscosity
    double viscosity_;
    /// thermal energy
    double kt_;
    /// time step for stochastic events concerning crosslinking
    double deltatime_;
    /// box corners
    Core::LinAlg::Matrix<3, 2> init_box_;
    /// number of crosslinker that are initially set
    std::map<int, int> maxnum_init_crosslinker_pertype_;
    /// number of crosslinkers in the simulated volume
    std::vector<int> numcrosslinkerpertype_;
    /// material numbers for crosslinker types
    std::vector<int> matcrosslinkerpertype_;
    /// linker and therefore binding spot types
    std::vector<Inpar::BEAMINTERACTION::CrosslinkerType> linkertypes_;
    /// maximal number of bonds per filament binding spot
    std::map<Inpar::BEAMINTERACTION::CrosslinkerType, int> max_num_bonds_per_filament_bspot_;
    /// distance between two binding spots on each filament
    std::map<Inpar::BEAMINTERACTION::CrosslinkerType, double> filamentbspotintervalglobal_;
    /// distance between two binding spots on a filament as percentage of filament reference length
    std::map<Inpar::BEAMINTERACTION::CrosslinkerType, double> filamentbspotintervallocal_;
    /// start and end arc parameter for binding spots on a filament
    std::map<Inpar::BEAMINTERACTION::CrosslinkerType, std::pair<double, double>>
        filamentbspotrangeglobal_;
    /// start and end arc parameter for binding spots on a filament
    /// in percent of filament reference length
    std::map<Inpar::BEAMINTERACTION::CrosslinkerType, std::pair<double, double>>
        filamentbspotrangelocal_;
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
