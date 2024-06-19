/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all beam to sphere contact input parameters

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

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
    void init();

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

    inline double beam_to_sphere_penalty_param() const { return penalty_parameter_; }

   private:
    bool isinit_;

    bool issetup_;

    //! beam-to-sphere penalty parameter
    double penalty_parameter_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
