/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fsi.H

 \brief FSI Wrapper for the ALE time integration

 <pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_FSI_H_
#define AD_ALE_FSI_H_


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
class AleFsiWrapper : public AleWrapper
{
public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleFsiWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleFsiWrapper(){}

  //@}

  //! communicate object at the interface
  Teuchos::RCP<const ALENEW::UTILS::MapExtractor> Interface() const;

  //! ToDo (mayr) move this to XFluidFluid adapter
  virtual void SolveAleXFluidFluidFSI() { };

  //! ToDo (mayr) move this to XFluidFluid adapter
  virtual void SolveBioGr() { };

  //! apply interface displacements
  void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
  {
    interface_->InsertFSICondVector(idisp,WriteAccessDispnp());
  }

  //! get Dirichlet map extractor
  virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
  {
    return AleWrapper::GetDBCMapExtractor();
  }

private:

  Teuchos::RCP<ALENEW::UTILS::MapExtractor> interface_;

}; // class AleFsiWrapper
} // namespace ADAPTER

#endif /* AD_ALE_FSI_H_ */
