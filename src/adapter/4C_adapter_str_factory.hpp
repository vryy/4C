/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_FACTORY_HPP
#define FOUR_C_ADAPTER_STR_FACTORY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class StructureBaseAlgorithmNew;
  class StructureFactory
  {
   public:
    //! constructor
    StructureFactory();

    //! destructor
    virtual ~StructureFactory() = default;

    //! Build the structural adapter object
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> build_structure_algorithm(
        const Teuchos::ParameterList& sdyn) const;
  };  // class Factory

  // non-member function
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> build_structure_algorithm(
      const Teuchos::ParameterList& sdyn);
}  // namespace ADAPTER


FOUR_C_NAMESPACE_CLOSE

#endif
