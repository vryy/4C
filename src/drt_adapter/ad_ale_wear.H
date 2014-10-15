/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_wear.H

 \brief Wrapper for the ALE time integration

 <pre>
       Maintainer: Philip Farah
       farah@lnm.mw.tum.de
       089 - 289-15257
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_WEAR_H_
#define AD_ALE_WEAR_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{


class AleWearWrapper : public AleWrapper
{
public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleWearWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleWearWrapper(){}

  //@}



}; // class AleWearWrapper
} // namespace ADAPTER

#endif /* AD_ALE_WEAR_H_ */
