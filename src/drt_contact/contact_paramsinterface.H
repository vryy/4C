/*---------------------------------------------------------------------*/
/*!
\file contact_paramsinterface.H

\brief Contact parameter interface. Necessary for the communication
       between the structural time integration framework and the
       contact strategies.

\level 3

\maintainer Michael Hiermeier

\date May 2, 2016

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_PARAMSINTERFACE_H_
#define SRC_DRT_CONTACT_CONTACT_PARAMSINTERFACE_H_

#include "../drt_mortar/mortar_paramsinterface.H"

namespace CONTACT {
  class ParamsInterface : public MORTAR::ParamsInterface
  {
  public:
    //! destructor
    virtual ~ParamsInterface() {};

    //! get the predictor status
    virtual const bool& IsPredictor() const = 0;

    //! get the current step length
    virtual double GetStepLength() const = 0;

    //! get the is_default_step indicator
    virtual bool IsDefaultStep() const = 0;

    //! get the current time step [derived]
    virtual double GetDeltaTime() const = 0;
  };
}


#endif /* SRC_DRT_CONTACT_CONTACT_PARAMSINTERFACE_H_ */
