/*----------------------------------------------------------------------------*/
/*! \file

 \brief FSI Wrapper for the ALE time integration

 \level 1

 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_ALE_FSI_HPP
#define BACI_ADAPTER_ALE_FSI_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_config.hpp"

#include "baci_adapter_ale_wrapper.hpp"

BACI_NAMESPACE_OPEN

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
  /*! \brief ALE Wrapper for FSI Problems
   *
   *  Provide FSI specific ALE functionalities here by overloading the respective
   *  routines from ADAPTER::AleWrapper
   *
   *  \sa ADAPTER::Ale, ADAPTER::AleWrapper
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

    //@}

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

    //! apply interface displacements
    void ApplyInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp)
    {
      interface_->InsertFSICondVector(idisp, WriteAccessDispnp());
    }

    //! get Dirichlet map extractor
    Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() override
    {
      return AleWrapper::GetDBCMapExtractor();
    }

   private:
    Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;

  };  // class AleFsiWrapper
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_ALE_FSI_H
