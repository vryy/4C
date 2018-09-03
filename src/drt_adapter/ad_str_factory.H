/*-----------------------------------------------------------*/
/*!
\file ad_str_factory.H

\maintainer Michael Hiermeier

\date Jan 5, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_STR_FACTORY_H_
#define SRC_DRT_ADAPTER_AD_STR_FACTORY_H_

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
      virtual ~Factory(){};

      //! Build the structural adapter object
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> BuildStructureAlgorithm(
          const Teuchos::ParameterList& sdyn) const;
    };  // class Factory

    // non-member function
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> BuildStructureAlgorithm(
        const Teuchos::ParameterList& sdyn);
  }  // namespace STR
}  // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_AD_STR_FACTORY_H_ */
