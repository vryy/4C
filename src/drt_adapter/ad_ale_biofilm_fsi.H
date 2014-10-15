/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_biofilm_fsi.H

 \brief Biofilm FSI Wrapper for the ALE time integration

 <pre>
 Maintainer: Christoph Ager
            ager@lnm.mw.tum.de
            089 - 289-15249
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_BIOFILM_FSI_H_
#define AD_ALE_BIOFILM_FSI_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fsi.H"


/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{


class AleBiofilmFsiWrapper : public AleFsiWrapper
{
public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleBiofilmFsiWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleBiofilmFsiWrapper(){}

  //@}

  virtual void SolveBioGr();


}; // class AleBiofilmFsiWrapper
} // namespace ADAPTER

#endif /* AD_ALE_BIOFILM_FSI_H_ */
