/*----------------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the ALE time integration for fluid problems with moving boundaries

\level 2


 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_ALE_FLUID_HPP
#define BACI_ADAPTER_ALE_FLUID_HPP


/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_config.hpp"

#include "baci_adapter_ale_wrapper.hpp"

BACI_NAMESPACE_OPEN

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

    //@}

    //! communicate object at the interface
    Teuchos::RCP<const ALE::UTILS::MapExtractor> Interface() const;

    //! solve
    int Solve() override;

    //! apply displacements at the free surface nodes
    void ApplyFreeSurfaceDisplacements(Teuchos::RCP<const Epetra_Vector> fsdisp);

    //! apply displacements at the ale update condition nodes
    void ApplyAleUpdateDisplacements(Teuchos::RCP<const Epetra_Vector> audisp);

    //! apply FSI interface displacements (required for staggered FSI)
    void ApplyInterfaceDisplacements(Teuchos::RCP<const Epetra_Vector> idisp);

   private:
    Teuchos::RCP<ALE::UTILS::MapExtractor> interface_;

  };  // class AleFluidWrapper
}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_ALE_FLUID_H
