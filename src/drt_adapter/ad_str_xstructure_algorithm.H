/*----------------------------------------------------------------------------*/
/**
\file ad_str_xstructure_algorithm.H

\brief Adaption for the xstructure case


\maintainer Michael Hiermeier

\date Jun 29, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_STR_XSTRUCTURE_ALGORITHM_H_
#define SRC_DRT_ADAPTER_AD_STR_XSTRUCTURE_ALGORITHM_H_


#include "ad_str_structure_new.H"

namespace ADAPTER
{
  class XStructureAlgorithm : public StructureBaseAlgorithmNew
  {
   public:
    /// constructor
    XStructureAlgorithm(){};

    /// destructor
    virtual ~XStructureAlgorithm(){};

   protected:
    /** \brief Create and initiate the global state data container
     *
     *  Please note, that we do NOT setup the global state data container.
     *  This is due to the fact, that we need the full dof row map to create
     *  the different state variables. The full dof row map is unknown at
     *  this point. First we need the merged FEM/xFEM discretization to do
     *  it. These things are done in the Setup() routine of the
     *  ADAPTER::StructureXContact class. */
    virtual void SetGlobalState(Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& dataglobalstate,
        const Teuchos::RCP<const STR::TIMINT::BaseDataSDyn>& datasdyn);

    /** \brief Create and initiate the time integration strategy object
     *
     *  Please note, that we do not Setup() the time integration object,
     *  because of the uncompleted global state object. */
    virtual void SetTimeIntegrationStrategy(Teuchos::RCP<STR::TIMINT::Base>& ti_strategy,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& dataio,
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& datasdyn,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& dataglobalstate, const int& restart);
  };

}  // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_AD_STR_XSTRUCTURE_ALGORITHM_H_ */
