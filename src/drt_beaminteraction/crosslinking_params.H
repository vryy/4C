/*---------------------------------------------------------------------*/
/*!
\file crosslinking_params.H

\brief data container holding all crosslinking input parameters

\level 3

\maintainer Jonas Eichinger, Maximilian Grill
*/
/*---------------------------------------------------------------------*/

#ifndef CROSSLINKINGPARAMS_H
#define CROSSLINKINGPARAMS_H

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_beaminteraction.H"

// forward declaration
namespace STR {
namespace TIMINT {
  class BaseDataGlobalState;
}
}
namespace BEAMINTERACTION {
  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction
   * author eichinger*/
  class CrosslinkingParams
  {
  public:
    //! constructor
    CrosslinkingParams();

    //! destructor
    virtual ~CrosslinkingParams() {};

    //! initialize with the stuff coming from input file
    void Init( STR::TIMINT::BaseDataGlobalState const& gstate );

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline const bool& IsInit() const
    { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const
    { return issetup_; };

    //! Checks the init and setup status
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    //! Checks the init status
    inline void CheckInit() const
    {
      if (!IsInit())
        dserror("Init() has not been called, yet!");
    }

    /// number of crosslinkers per type
    std::vector<int> const& NumCrosslinkerPerType() const
    {
      CheckInitSetup();
      return numcrosslinkerpertype_;
    };

    /// material number for crosslinker types
    std::vector<int> const& MatCrosslinkerPerType() const
    {
      CheckInitSetup();
      return matcrosslinkerpertype_;
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

    // distance between two binding spots on a filament
    double const& FilamentBspotInterval() const
    {
      CheckInitSetup();
      return filamentbspotinterval_;
    };

    // start and end arc parameter for binding spots on a filament
    std::vector<double> const& FilamentBspotRange() const
    {
      CheckInitSetup();
      return filamentbspotrange_;
    };

  private:

    bool isinit_;

    bool issetup_;

    /// number of crosslinkers in the simulated volume
    std::vector<int> numcrosslinkerpertype_;
    /// material numbers for crosslinker types
    std::vector<int> matcrosslinkerpertype_;
    /// viscosity
    double viscosity_;
    /// thermal energy
    double kt_;
    /// time step for stochastic events concering crosslinking
    double deltatime_;
    /// distance between two binding spots on a filament
    double filamentbspotinterval_;
    /// start and end arc parameter for binding spots on a filament
    std::vector<double> filamentbspotrange_;

  };
}

#endif
