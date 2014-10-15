/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fsi.H

 \brief Wrapper for the ALE time integration

 <pre>
       Maintainer: Matthias Mayr
       mayr@mhpc.mw.tum.de
       089 - 289-10362
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_XFFSI_H_
#define AD_ALE_XFFSI_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fsi.H"

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
class AleXFFsiWrapper : public AleFsiWrapper
{
public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleXFFsiWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleXFFsiWrapper(){}

  //@}

  virtual void SolveAleXFluidFluidFSI();

private:

}; // class AleXFFsiWrapper
} // namespace ADAPTER

#endif /* AD_ALE_XFFSI_H_ */
