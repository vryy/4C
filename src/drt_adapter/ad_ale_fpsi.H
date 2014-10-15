/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fpsi.H

 \brief FPSI wrapper for the ALE time integration

 <pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            089 - 289-15240
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_FPSI_H_
#define AD_ALE_FPSI_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fsi.H"

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{

class AleFpsiWrapper : public AleFsiWrapper
{
public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleFpsiWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleFpsiWrapper(){}

  //@}


}; // class AleFpsiWrapper
} // namespace ADAPTER

#endif /* AD_ALE_FPSI_H_ */
