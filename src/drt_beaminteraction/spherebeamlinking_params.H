/*---------------------------------------------------------------------*/
/*!
\file spherebeamlinking_params.H


\brief data container holding all contractile cells input parameters

\level 3

\maintainer Jonas Eichinger
*/
/*---------------------------------------------------------------------*/

#ifndef SPHEREBEAMLINKING_PARAMS_H
#define SPHEREBEAMLINKING_PARAMS_H

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include<Teuchos_RCP.hpp>



// forward declaration

namespace STR {
namespace TIMINT {
  class BaseDataGlobalState;
}
}
namespace MAT
{
  class CrosslinkerMat;
}
namespace BEAMINTERACTION {

  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction */
  class SphereBeamLinkingParams
  {
  public:
    //! constructor
    SphereBeamLinkingParams();

    //! destructor
    virtual ~SphereBeamLinkingParams() {};

    //! initialize with the stuff coming from input file
    void Init( STR::TIMINT::BaseDataGlobalState const & gstate );

    //! setup member variables
    void Setup();

    //! reset time step in case structure time is adapted during simulation time
    void ResetTimeStep( double structure_delta_time );

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

    /// linker material id
    Teuchos::RCP< MAT::CrosslinkerMat > GetLinkerMaterial() const
    {
      /// HACK: FIX IF MORE THAN ONE CROSSLINKER TYPE
      CheckInitSetup();
      return mat_.back();
    };

    /// time step for stochastic events concerning crosslinking
    double const & DeltaTime() const
    {
      CheckInitSetup();
      return deltatime_;
    };

    /// contraction rate of cell (integrin linker) in [μm/s]
    double ContractionRate( INPAR::BEAMINTERACTION::CrosslinkerType linkertype ) const
    {
      CheckInitSetup();
      return contractionrate_.at(linkertype);
    };

    /// number of linker per type
    std::vector<int> const & MaxNumLinkerPerType() const
    {
      CheckInitSetup();
      return maxnumlinkerpertype_;
    };

    /// material number for linker types
    std::vector<int> const & MatLinkerPerType() const
    {
      CheckInitSetup();
      return matlinkerpertype_;
    };

    /// get all active linker types
    std::vector< INPAR::BEAMINTERACTION::CrosslinkerType > const & LinkerTypes() const
    {
      CheckInitSetup();
      return linkertypes_;
    };

    // distance between two binding spots on a filament
    double FilamentBspotIntervalGlobal( INPAR::BEAMINTERACTION::CrosslinkerType linkertype ) const
    {
      CheckInitSetup();
      return filamentbspotintervalglobal_.at( linkertype );
    };

    // distance between two binding spots on a filament
    double FilamentBspotIntervalLocal( INPAR::BEAMINTERACTION::CrosslinkerType linkertype ) const
    {
      CheckInitSetup();
      return filamentbspotintervallocal_.at( linkertype );
    };

    // start and end arc parameter for binding spots on a filament
    std::pair< double, double > const &
    FilamentBspotRangeLocal( INPAR::BEAMINTERACTION::CrosslinkerType linkertype ) const
    {
      CheckInitSetup();
      return filamentbspotrangelocal_.at( linkertype );
    };

    // start and end arc parameter for binding spots on a filament
    std::pair< double, double > const &
    FilamentBspotRangeGlobal( INPAR::BEAMINTERACTION::CrosslinkerType linkertype ) const
    {
      CheckInitSetup();
      return filamentbspotrangeglobal_.at( linkertype );
    };

  private:

    bool isinit_;

    bool issetup_;

    /// time step for stochastic events concerning integrins, e.g. catch-slip-bond behavior
    double deltatime_;
    bool own_deltatime_;
    /// contraction rate of cell (integrin linker) in [μm/s]
    std::map< INPAR::BEAMINTERACTION::CrosslinkerType, double > contractionrate_;
    /// crosslinker material
    std::vector< Teuchos::RCP< MAT::CrosslinkerMat > > mat_;
    /// number of crosslinkers in the simulated volume
    std::vector<int> maxnumlinkerpertype_;
    /// material numbers for crosslinker types
    std::vector<int> matlinkerpertype_;
    /// linker and therefore binding spot types
    std::vector< INPAR::BEAMINTERACTION::CrosslinkerType > linkertypes_;
    /// distance between two binding spots on each filament
    std::map< INPAR::BEAMINTERACTION::CrosslinkerType, double > filamentbspotintervalglobal_;
    /// distance between two binding spots on a filament as percentage of filament reference length
    std::map< INPAR::BEAMINTERACTION::CrosslinkerType, double > filamentbspotintervallocal_;
    /// start and end arc parameter for binding spots on a filament
    std::map< INPAR::BEAMINTERACTION::CrosslinkerType, std::pair< double, double > > filamentbspotrangeglobal_;
    /// start and end arc parameter for binding spots on a filament
    /// in percent of filament reference length
    std::map< INPAR::BEAMINTERACTION::CrosslinkerType, std::pair< double, double > > filamentbspotrangelocal_;

  };

}

#endif
