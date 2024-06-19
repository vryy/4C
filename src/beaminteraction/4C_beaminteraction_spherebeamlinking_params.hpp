/*---------------------------------------------------------------------*/
/*! \file


\brief data container holding all contractile cells input parameters

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_SPHEREBEAMLINKING_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_SPHEREBEAMLINKING_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_inpar_beaminteraction.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN



// forward declaration

namespace STR
{
  namespace TimeInt
  {
    class BaseDataGlobalState;
  }
}  // namespace STR
namespace Mat
{
  class CrosslinkerMat;
}
namespace BEAMINTERACTION
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
    void init(STR::TimeInt::BaseDataGlobalState const& gstate);

    //! setup member variables
    void setup();

    //! reset time step in case structure time is adapted during simulation time
    void ResetTimeStep(double structure_delta_time);

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
    Teuchos::RCP<Mat::CrosslinkerMat> GetLinkerMaterial() const
    {
      /// HACK: FIX IF MORE THAN ONE CROSSLINKER TYPE
      check_init_setup();
      return mat_.back();
    };

    /// time step for stochastic events concerning crosslinking
    double const& DeltaTime() const
    {
      check_init_setup();
      return deltatime_;
    };

    /// contraction rate of cell (integrin linker) in [microm/s]
    double ContractionRate(Inpar::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      check_init_setup();
      return contractionrate_.at(linkertype);
    };

    /// number of linker per type
    std::vector<int> const& MaxNumLinkerPerType() const
    {
      check_init_setup();
      return maxnumlinkerpertype_;
    };

    /// material number for linker types
    std::vector<int> const& MatLinkerPerType() const
    {
      check_init_setup();
      return matlinkerpertype_;
    };

    /// get all active linker types
    std::vector<Inpar::BEAMINTERACTION::CrosslinkerType> const& LinkerTypes() const
    {
      check_init_setup();
      return linkertypes_;
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

    /// time step for stochastic events concerning integrins, e.g. catch-slip-bond behavior
    double deltatime_;
    bool own_deltatime_;
    /// contraction rate of cell (integrin linker) in [microm/s]
    std::map<Inpar::BEAMINTERACTION::CrosslinkerType, double> contractionrate_;
    /// crosslinker material
    std::vector<Teuchos::RCP<Mat::CrosslinkerMat>> mat_;
    /// number of crosslinkers in the simulated volume
    std::vector<int> maxnumlinkerpertype_;
    /// material numbers for crosslinker types
    std::vector<int> matlinkerpertype_;
    /// linker and therefore binding spot types
    std::vector<Inpar::BEAMINTERACTION::CrosslinkerType> linkertypes_;
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
