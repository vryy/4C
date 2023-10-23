/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_ADAPTER_STR_FACTORY_H
#define BACI_ADAPTER_STR_FACTORY_H

#include <Teuchos_RCP.hpp>

// forward declaration
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace ADAPTER
{
  class StructureBaseAlgorithmNew;
  namespace STR
  {
    class Factory
    {
     public:
      //! constructor
      Factory();

      //! destructor
      virtual ~Factory() = default;

      //! Build the structural adapter object
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> BuildStructureAlgorithm(
          const Teuchos::ParameterList& sdyn) const;
    };  // class Factory

    // non-member function
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> BuildStructureAlgorithm(
        const Teuchos::ParameterList& sdyn);
  }  // namespace STR
}  // namespace ADAPTER


#endif  // ADAPTER_STR_FACTORY_H
