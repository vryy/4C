/*----------------------------------------------------------------------------*/
/*!

 \brief Wrapper for the ALE time integration

 \maintainer Matthias Mayr

 \level 2
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef AD_ALE_WEAR_H_
#define AD_ALE_WEAR_H_


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace ALE
{
  namespace UTILS
  {
    class MapExtractor;
  }
}  // namespace ALE

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
    virtual ~AleWearWrapper() {}

    //@}

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

    //! add ALE wear condition vector
    void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
    {
      interface_->AddAleWearCondVector(idisp, WriteAccessDispnp());
    }

   private:
    Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;

  };  // class AleWearWrapper
}  // namespace ADAPTER

#endif /* AD_ALE_WEAR_H_ */
