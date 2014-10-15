/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_crack.H

 \brief Crack Wrapper for the ALE time integration

 <pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            089 - 289-15257
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_CRACK_H_
#define AD_ALE_CRACK_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace ALENEW
{
namespace UTILS
{
  class MapExtractor;
}
}

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{

/*! \brief Just a wrapper that does nothing, meant to be derived from
 *
 *  This wrapper just encapsulated the ADAPTER::Ale and implements all
 *  routines that are pure virtual in ADAPTER::Ale. For a specific ALE adapter
 *  just derive from this one and overload those routines you need with your
 *  problem specific routine.
 *
 *  \author mayr.mt \date 10/2014
 */
class AleCrackWrapper : public AleWrapper
{
public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleCrackWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleCrackWrapper(){}

  //@}

private:

}; // class AleCrackWrapper
} // namespace ADAPTER

#endif /* AD_ALE_CRACK_H_ */
