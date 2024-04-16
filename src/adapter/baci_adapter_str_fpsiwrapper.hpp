/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for FPSI problems containing the interface
       and methods dependent on the interface


\level 3

*/

#ifndef FOUR_C_ADAPTER_STR_FPSIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_FPSIWRAPPER_HPP

#include "baci_config.hpp"

#include "baci_adapter_str_fsiwrapper.hpp"

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class FPSIStructureWrapper : public FSIStructureWrapper
  {
   public:
    /// constructor
    explicit FPSIStructureWrapper(Teuchos::RCP<Structure> structure);

    /*!
    \brief extract interface displacements at \f$t_{n}\f$

    \param FPSI (in) : if true perform some pre-stress and fpsi specific stuff

    \note if param FPSI = false the base class version of this method is called

    */
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn(bool FPSI = false);

    /*!
    \brief  extract interface displacements at \f$t_{n+1}\f$

    \param FPSI (in) : if true perform some pre-stress and fpsi specific stuff

    \note if param FPSI = false the base class version of this method is called

    */
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp(bool FPSI = false);
  };
}  // namespace ADAPTER
BACI_NAMESPACE_CLOSE

#endif
