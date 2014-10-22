/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fs.H

 \brief Wrapper for the ALE time integration in free-surface problem

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

/*! \brief ALE Wrapper for free-surface problems
 *
 *  \sa ADAPTER::Ale, ADAPTER::AleWrapper
 *
 *  \author kruse \date 10/2014
 */
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

  //! apply displacements at the free surface nodes
  void ApplyFreeSurfaceDisplacements(Teuchos::RCP<const Epetra_Vector> fsdisp);

  //! apply displacements at the ale update condition nodes
  void ApplyAleUpdateDisplacements(Teuchos::RCP<const Epetra_Vector> audisp);

private:

  //! map extractor for conditioned nodesets of the free-surface problem
  Teuchos::RCP<ALENEW::UTILS::MapExtractor> interface_;

}; // class AleFSWrapper
} // namespace ADAPTER

#endif /* AD_ALE_FS_H_ */
