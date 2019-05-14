/*-----------------------------------------------------------*/
/*!

\maintainer Anh-Tu Vuong

\brief  structure adapter using the LOCA library


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_STR_LOCA_H_
#define SRC_DRT_ADAPTER_AD_STR_LOCA_H_

#include "ad_str_structure_new.H"

namespace ADAPTER
{
  class StructureLocaAlgorithm : public StructureBaseAlgorithmNew
  {
   public:
    //! constructor
    StructureLocaAlgorithm();

    //! destructor
    virtual ~StructureLocaAlgorithm(){};

    //! setup (derived)
    virtual void Setup();

   protected:
    //! // Set different time integrator specific parameters in the different parameter lists
    //! (derived)
    virtual void SetParams(Teuchos::ParameterList& ioflags, Teuchos::ParameterList& xparams,
        Teuchos::ParameterList& taflags);

    /// set the final structure time integrator object (derived)
    virtual void SetStructureWrapper(const Teuchos::ParameterList& ioflags,
        const Teuchos::ParameterList& sdyn, const Teuchos::ParameterList& xparams,
        const Teuchos::ParameterList& taflags, Teuchos::RCP<::STR::TIMINT::Base> ti_strategy);
  };  // class StructureLocaAlgorithm
}  // namespace ADAPTER



#endif /* SRC_DRT_ADAPTER_AD_STR_LOCA_H_ */
