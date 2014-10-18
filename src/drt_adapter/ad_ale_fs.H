/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fs.H

 \brief Wrapper for the ALE time integration

 <pre>
 Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_FS_H_
#define AD_ALE_FS_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"

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

class AleFSWrapper : public AleWrapper
{

public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleFSWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleFSWrapper(){}

  //@}

  //! communicate object at the interface
  Teuchos::RCP<const ALENEW::UTILS::MapExtractor> Interface() const;

  /// displacements at the free surface are prescribed
  void ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp);

  /// displacements at the ale update condition nodes are prescribed
  void ApplyAleUpdateDisplacements(Teuchos::RCP<Epetra_Vector> audisp);

private:

  Teuchos::RCP<ALENEW::UTILS::MapExtractor> interface_;

}; // class AleFSWrapper
} // namespace ADAPTER

#endif /* AD_ALE_FS_H_ */
