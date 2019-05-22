/*---------------------------------------------------------------------*/
/*!
\brief data container holding all crosslinking input parameters

\level 3

\maintainer Jonas Eichinger
*/
/*---------------------------------------------------------------------*/

#ifndef CROSSLINKINGPARAMS_H
#define CROSSLINKINGPARAMS_H

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include "../linalg/linalg_fixedsizematrix.H"

// forward declaration
namespace STR
{
  namespace TIMINT
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
    virtual ~CrosslinkingParams(){};

    //! initialize with the stuff coming from input file
    void Init(STR::TIMINT::BaseDataGlobalState const& gstate);

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline const bool& IsInit() const { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const { return issetup_; };

    //! Checks the init and setup status
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
    }

    //! Checks the init status
    inline void CheckInit() const
    {
      if (!IsInit()) dserror("Init() has not been called, yet!");
    }

    /// number of crosslinkers per type
    std::vector<int> const& NumCrosslinkerPerType() const
    {
      CheckInitSetup();
      return numcrosslinkerpertype_;
    };

    /// number of crosslinkers per type
    int NumInitCrosslinkerPerCrosslinkerMatId(int matid) const
    {
      CheckInitSetup();
      return maxnum_init_crosslinker_pertype_.at(matid);
    };

    /// number of crosslinkers per type
    int TotalNumInitCrosslinker() const
    {
      CheckInitSetup();
      int sum = 0;
      for (auto const& iter : maxnum_init_crosslinker_pertype_) sum += iter.second;
      return sum;
    };

    /// material number for crosslinker types
    std::vector<int> const& MatCrosslinkerPerType() const
    {
      CheckInitSetup();
      return matcrosslinkerpertype_;
    };

    /// get all active crosslinker types
    std::vector<INPAR::BEAMINTERACTION::CrosslinkerType> const& LinkerTypes() const
    {
      CheckInitSetup();
      return linkertypes_;
    };

    /// number of different crosslinker types in simulation volume
    int NumberOfCrosslinkerTypes() const
    {
      CheckInitSetup();
      return static_cast<int>(numcrosslinkerpertype_.size());
    };

    /// ~ 1e-3 / 2.27 according to cyron2011 eq 52 ff, viscosity of surrounding fluid
    double const& Viscosity() const
    {
      CheckInitSetup();
      return viscosity_;
    };

    /// thermal energy
    double const& KT() const
    {
      CheckInitSetup();
      return kt_;
    };

    /// time step for stochastic events concerning crosslinking
    double const& DeltaTime() const
    {
      CheckInitSetup();
      return deltatime_;
    };

    /// time step for stochastic events concerning crosslinking
    LINALG::Matrix<3, 2> const& LinkerInitializationBox() const
    {
      CheckInitSetup();
      return init_box_;
      ;
    };

    // distance between two binding spots on a filament
    int MaxNumberOfBondsPerFilamentBspot(INPAR::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      CheckInitSetup();
      return max_num_bonds_per_filament_bspot_.at(linkertype);
    };

    // distance between two binding spots on a filament
    double FilamentBspotIntervalGlobal(INPAR::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      CheckInitSetup();
      return filamentbspotintervalglobal_.at(linkertype);
    };

    // distance between two binding spots on a filament
    double FilamentBspotIntervalLocal(INPAR::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      CheckInitSetup();
      return filamentbspotintervallocal_.at(linkertype);
    };

    // start and end arc parameter for binding spots on a filament
    std::pair<double, double> const& FilamentBspotRangeLocal(
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      CheckInitSetup();
      return filamentbspotrangelocal_.at(linkertype);
    };

    // start and end arc parameter for binding spots on a filament
    std::pair<double, double> const& FilamentBspotRangeGlobal(
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype) const
    {
      CheckInitSetup();
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
    LINALG::Matrix<3, 2> init_box_;
    /// number of crosslinker that are initially set
    std::map<int, int> maxnum_init_crosslinker_pertype_;
    /// number of crosslinkers in the simulated volume
    std::vector<int> numcrosslinkerpertype_;
    /// material numbers for crosslinker types
    std::vector<int> matcrosslinkerpertype_;
    /// linker and therefore binding spot types
    std::vector<INPAR::BEAMINTERACTION::CrosslinkerType> linkertypes_;
    /// maximal number of bonds per filament binding spot
    std::map<INPAR::BEAMINTERACTION::CrosslinkerType, int> max_num_bonds_per_filament_bspot_;
    /// distance between two binding spots on each filament
    std::map<INPAR::BEAMINTERACTION::CrosslinkerType, double> filamentbspotintervalglobal_;
    /// distance between two binding spots on a filament as percentage of filament reference length
    std::map<INPAR::BEAMINTERACTION::CrosslinkerType, double> filamentbspotintervallocal_;
    /// start and end arc parameter for binding spots on a filament
    std::map<INPAR::BEAMINTERACTION::CrosslinkerType, std::pair<double, double>>
        filamentbspotrangeglobal_;
    /// start and end arc parameter for binding spots on a filament
    /// in percent of filament reference length
    std::map<INPAR::BEAMINTERACTION::CrosslinkerType, std::pair<double, double>>
        filamentbspotrangelocal_;
  };
}  // namespace BEAMINTERACTION

#endif
