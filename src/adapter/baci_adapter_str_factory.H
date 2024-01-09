/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_ADAPTER_STR_FACTORY_H
#define BACI_ADAPTER_STR_FACTORY_H

#include "baci_config.H"

#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

BACI_NAMESPACE_OPEN

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
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> BuildStructureAlgorithm(
        const Teuchos::ParameterList& sdyn) const;
  };  // class Factory

  // non-member function
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> BuildStructureAlgorithm(
      const Teuchos::ParameterList& sdyn);
}  // namespace ADAPTER


BACI_NAMESPACE_CLOSE

#endif  // ADAPTER_STR_FACTORY_H
