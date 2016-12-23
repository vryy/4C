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
#include "../drt_inpar/inpar_crosslinking.H"

// forward declaration
namespace BEAMINTERACTION {

  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction */
  class CrosslinkingParams
  {
  public:
    //! constructor
    CrosslinkingParams();

    //! destructor
    virtual ~CrosslinkingParams() {};

    //! initialize with the stuff coming from input file
    void Init();

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

    /// number of crosslinkers in the simulated volume
    int const& NumCrosslink() const
    {
      CheckInitSetup();
      return numcrosslink_;
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

    /// distance between the two binding domains of a linker
    double const& LinkingLength() const
    {
      CheckInitSetup();
      return linkinglength_;
    };

    /// tolerance for linker length in the sense: length +- tolerance
    double const& LinkingLengthTolerance() const
    {
      CheckInitSetup();
      return linkinglengthtol_;
    };

    /// assuming circular cross section of crosslinker with diameter according to phd tharmann
    double const& ALink() const
    {
      CheckInitSetup();
      return alink_;
    };

    /// must set ilink=0.0 ->truss!
    double const& ILink() const
    {
      CheckInitSetup();
      return ilink_;
    };

    /// 2xilink
    double const& IpLink() const
    {
      CheckInitSetup();
      return iplink_;
    };

    /// preferred binding angle enclosed by two filaments' axes in radians
    double const& LinkingAngle() const
    {
     CheckInitSetup();
     return linkingangle_;
    };

    /// tolerance for preferred binding angle in radians in the sense of: angle +- tolerance
    double const& LinkingAngleTolerance() const
    {
     CheckInitSetup();
     return linkingangletol_;
    };

    /// chemical dissociation-rate after ktswitchtime // 100.0 1000.0
    double const& KOn() const
    {
      CheckInitSetup();
      return k_on_;
    };

    /// in analogy to k_on
    double const& KOff() const
    {
      CheckInitSetup();
      return k_off_;
    };

  private:

    bool isinit_;

    bool issetup_;

    // number of crosslinkers in the simulated volume
    int numcrosslink_;
    // viscosity
    double viscosity_;
    // thermal energy
    double kt_;
    // distance between the two binding domains of a linker
    double linkinglength_;
    // deviation from r_link
    double linkinglengthtol_;
    // assuming circular cross section of crosslinker with diameter according to phd tharmann
    double alink_;
    // must set ilink=0.0 ->truss!
    double ilink_;
    // 2xilink
    double iplink_;
    // preferred binding angle betw. two filaments
    double linkingangle_;
    // deviation from phi0
    double linkingangletol_;
    // chemical association-rate before ktswitchtime
    double k_on_;
    // chemical dissociation-rate after ktswitchtime // 100.0 1000.0
    // in analogy to k_on
    double k_off_;

  };

}

#endif
