/*----------------------------------------------------------------------------*/
/*!
\file ad_ale_fluid.H

\brief Wrapper for the ALE time integration for fluid problems with moving boundaries

\level 2

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_FLUID_H_
#define AD_ALE_FLUID_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"

/* forward declarations */
namespace ALE
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

/*! \brief ALE Wrapper for fluid problems with moving boundaries
 *
 *  Provide ALE functionalities for problems involving moving boundaries
 *  here by overloading the respective routines from ADAPTER::AleWrapper
 *
 *  \sa ADAPTER::Ale, ADAPTER::AleWrapper
 *
 *  \author kruse \date 10/2014
 */
class AleFluidWrapper : public AleWrapper
{

public:

  //! @name Construction / Destruction
  //@{

  //! constructor
  explicit AleFluidWrapper(Teuchos::RCP<Ale> ale);

  //! Destructor
  virtual ~AleFluidWrapper(){}

  //@}

  //! communicate object at the interface
  Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

  //! solve
  int Solve();

  //! apply displacements at the free surface nodes
  void ApplyFreeSurfaceDisplacements(Teuchos::RCP<const Epetra_Vector> fsdisp);

  //! apply displacements at the ale update condition nodes
  void ApplyAleUpdateDisplacements(Teuchos::RCP<const Epetra_Vector> audisp);

  //! apply FSI interface displacements (required for staggered FSI)
  void ApplyInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp);

  //! add additional FSI interface displacements (required for cell migration)
  void AddInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp);

private:

  Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;

}; // class AleFluidWrapper
} // namespace ADAPTER

#endif /* AD_ALE_FLUID_H_ */
