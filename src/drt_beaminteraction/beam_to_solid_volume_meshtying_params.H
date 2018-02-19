/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_solid_volume_meshtying_params.H

\brief data container holding all beam to solid volume meshtying contact input parameters

\level 3

\maintainer Alexander Popp
*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BEAMTOSOLIDVOLUMEMESHTYINGPARAMS_H
#define BEAMTOSOLIDVOLUMEMESHTYINGPARAMS_H

#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_beamcontact.H"

namespace BEAMINTERACTION {

  class BeamToSolidVolumeMeshtyingParams
  {
  public:

    //! constructor
    BeamToSolidVolumeMeshtyingParams();

    //! destructor
    virtual ~BeamToSolidVolumeMeshtyingParams() {};

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

    inline double BeamToSolidVolumeMeshtyingPenaltyParam() const
    {
      return BTSVOLMT_penalty_param_;
    }


  private:

    bool isinit_;

    bool issetup_;

    double BTSVOLMT_penalty_param_;


  };

}

#endif
