/*----------------------------------------------------------------------------*/
/*! \file

 \brief Wrapper for the ALE time integration


 \level 2
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_ALE_WEAR_HPP
#define FOUR_C_ADAPTER_ALE_WEAR_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_config.hpp"

#include "baci_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

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

FOUR_C_NAMESPACE_CLOSE

#endif
