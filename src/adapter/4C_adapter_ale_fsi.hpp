#ifndef FOUR_C_ADAPTER_ALE_FSI_HPP
#define FOUR_C_ADAPTER_ALE_FSI_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace ALE
{
  namespace Utils
  {
    class MapExtractor;
  }
}  // namespace ALE

/*----------------------------------------------------------------------------*/
/* class definitions */
namespace Adapter
{
  /*! \brief ALE Wrapper for FSI Problems
   *
   *  Provide FSI specific ALE functionalities here by overloading the respective
   *  routines from Adapter::AleWrapper
   *
   *  \sa Adapter::Ale, Adapter::AleWrapper
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
    Teuchos::RCP<const ALE::Utils::MapExtractor> interface() const;

    //! apply interface displacements
    void apply_interface_displacements(Teuchos::RCP<const Core::LinAlg::Vector<double>> idisp)
    {
      interface_->insert_fsi_cond_vector(*idisp, *write_access_dispnp());
    }

    //! get Dirichlet map extractor
    Teuchos::RCP<const Core::LinAlg::MapExtractor> get_dbc_map_extractor() override
    {
      return AleWrapper::get_dbc_map_extractor();
    }

   private:
    Teuchos::RCP<ALE::Utils::MapExtractor> interface_;

  };  // class AleFsiWrapper
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
