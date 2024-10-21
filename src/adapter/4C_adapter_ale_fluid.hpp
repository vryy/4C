#ifndef FOUR_C_ADAPTER_ALE_FLUID_HPP
#define FOUR_C_ADAPTER_ALE_FLUID_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

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
  /*! \brief ALE Wrapper for fluid problems with moving boundaries
   *
   *  Provide ALE functionalities for problems involving moving boundaries
   *  here by overloading the respective routines from Adapter::AleWrapper
   *
   *  \sa Adapter::Ale, Adapter::AleWrapper
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

    //@}

    //! communicate object at the interface
    Teuchos::RCP<const ALE::Utils::MapExtractor> interface() const;

    //! solve
    int solve() override;

    //! apply displacements at the free surface nodes
    void apply_free_surface_displacements(Teuchos::RCP<const Core::LinAlg::Vector<double>> fsdisp);

    //! apply displacements at the ale update condition nodes
    void apply_ale_update_displacements(Teuchos::RCP<const Core::LinAlg::Vector<double>> audisp);

    //! apply FSI interface displacements (required for staggered FSI)
    void apply_interface_displacements(Teuchos::RCP<const Core::LinAlg::Vector<double>> idisp);

   private:
    Teuchos::RCP<ALE::Utils::MapExtractor> interface_;

  };  // class AleFluidWrapper
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
