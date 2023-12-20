/*--------------------------------------------------------------------------*/
/*! \file

\brief FSI Wrapper for the ALE time integration with internal mesh tying or mesh sliding interface


\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_ADAPTER_ALE_FSI_MSHT_H
#define BACI_ADAPTER_ALE_FSI_MSHT_H



/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_config.H"

#include "baci_adapter_ale_fsi.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace ALE
{
  namespace UTILS
  {
    class FsiMapExtractor;
  }
}  // namespace ALE

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace ADAPTER
{
  /*! \brief ALE Wrapper for FSI Problems ith internal mesh tying or mesh sliding interface
   *
   *  Provide an additional map extractor
   *
   *  \sa ADAPTER::Ale, ADAPTER::AleWrapper, ADAPTER::AleFsiWrapper
   *
   *  \author wirtz \date 02/2016
   */
  class AleFsiMshtWrapper : public AleFsiWrapper
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor
    explicit AleFsiMshtWrapper(Teuchos::RCP<Ale> ale);

    //@}

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::FsiMapExtractor> FsiInterface() const;

   private:
    Teuchos::RCP<ALE::UTILS::FsiMapExtractor> fsiinterface_;

  };  // class AleFsiWrapper
}  // namespace ADAPTER


BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_ALE_FSI_MSHT_H
