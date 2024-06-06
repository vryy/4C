/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for FPSI problems containing the interface
       and methods dependent on the interface


\level 3

*/

#ifndef FOUR_C_ADAPTER_STR_FPSIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_FPSIWRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
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
    virtual Teuchos::RCP<Epetra_Vector> extract_interface_dispn(bool FPSI = false);

    /*!
    \brief  extract interface displacements at \f$t_{n+1}\f$

    \param FPSI (in) : if true perform some pre-stress and fpsi specific stuff

    \note if param FPSI = false the base class version of this method is called

    */
    virtual Teuchos::RCP<Epetra_Vector> extract_interface_dispnp(bool FPSI = false);
  };
}  // namespace Adapter
FOUR_C_NAMESPACE_CLOSE

#endif
