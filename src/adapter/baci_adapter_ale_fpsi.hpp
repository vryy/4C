/*----------------------------------------------------------------------------*/
/*! \file

 \brief FPSI wrapper for the ALE time integration

 \level 2

 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_ALE_FPSI_HPP
#define FOUR_C_ADAPTER_ALE_FPSI_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_config.hpp"

#include "baci_adapter_ale_wrapper.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{
  class AleFpsiWrapper : public AleWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleFpsiWrapper(Teuchos::RCP<Ale> ale);

    //! specialized method to apply displacements to fpsi interface
    void ApplyInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp);

    //! specialized method to apply displacements to fsi interface
    void ApplyFSIInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp);

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

    //@}

   private:
    //! interface map extractor
    Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;


  };  // class AleFpsiWrapper
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_ALE_FPSI_H
