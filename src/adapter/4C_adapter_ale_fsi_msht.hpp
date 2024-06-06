/*--------------------------------------------------------------------------*/
/*! \file

\brief FSI Wrapper for the ALE time integration with internal mesh tying or mesh sliding interface


\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_ALE_FSI_MSHT_HPP
#define FOUR_C_ADAPTER_ALE_FSI_MSHT_HPP



/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_adapter_ale_fsi.hpp"

FOUR_C_NAMESPACE_OPEN

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
namespace Adapter
{
  /*! \brief ALE Wrapper for FSI Problems ith internal mesh tying or mesh sliding interface
   *
   *  Provide an additional map extractor
   *
   *  \sa Adapter::Ale, Adapter::AleWrapper, Adapter::AleFsiWrapper
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
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
