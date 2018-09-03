/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_sphere_contact_params.H

\brief data container holding all beam to sphere contact input parameters

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BEAMTOSPHERECONTACTPARAMS_H
#define BEAMTOSPHERECONTACTPARAMS_H

#include "../drt_lib/drt_dserror.H"

namespace BEAMINTERACTION
{
  class BeamToSphereContactParams
  {
   public:
    //! constructor
    BeamToSphereContactParams();

    //! destructor
    virtual ~BeamToSphereContactParams(){};

    //! initialize with the stuff coming from input file
    void Init();

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

    inline double BeamToSpherePenaltyParam() const { return BTSPH_penalty_param_; }

   private:
    bool isinit_;

    bool issetup_;

    //! beam-to-sphere penalty parameter
    double BTSPH_penalty_param_;
  };

}  // namespace BEAMINTERACTION

#endif
