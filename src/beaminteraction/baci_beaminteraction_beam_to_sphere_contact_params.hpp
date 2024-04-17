/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all beam to sphere contact input parameters

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PARAMS_HPP

#include "baci_config.hpp"

#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  class BeamToSphereContactParams
  {
   public:
    //! constructor
    BeamToSphereContactParams();

    //! destructor
    virtual ~BeamToSphereContactParams() = default;

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

    inline double BeamToSpherePenaltyParam() const { return penalty_parameter_; }

   private:
    bool isinit_;

    bool issetup_;

    //! beam-to-sphere penalty parameter
    double penalty_parameter_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
