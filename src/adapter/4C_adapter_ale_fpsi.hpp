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
#include "4C_config.hpp"

#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

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
    void apply_interface_displacements(Teuchos::RCP<const Epetra_Vector> idisp);

    //! specialized method to apply displacements to fsi interface
    void apply_fsi_interface_displacements(Teuchos::RCP<const Epetra_Vector> idisp);

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

    //@}

   private:
    //! interface map extractor
    Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;


  };  // class AleFpsiWrapper
}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif
