/*-----------------------------------------------------------*/
/*! \file


\brief wrapper for structure adapter using the LOCA library

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_STR_LOCA_WRAPPER_H_
#define SRC_DRT_ADAPTER_AD_STR_LOCA_WRAPPER_H_


#include "ad_str_wrapper.H"

namespace ADAPTER
{
  class StructureLocaWrapper : public StructureWrapper
  {
   public:
    //! constructor
    explicit StructureLocaWrapper(Teuchos::RCP<Structure> structure)
        : StructureWrapper(structure){};

    //! destructor
    virtual ~StructureLocaWrapper(){};

    //! actual time loop
    virtual int Integrate();

    //! wrapper for things that should be done before PrepareTimeStep is called
    virtual void PrePredict(){};

    //! wrapper for things that should be done after PrepareTimeStep is called
    virtual void PostPredict(){};

    //! wrapper for things that should be done before solving the nonlinear iterations
    virtual void PreSolve(){};

    //! wrapper for things that should be done after solving the nonlinear iterations
    virtual void PostSolve(){};

    //! wrapper for things that should be done before updating
    virtual void PreUpdate(){};

    //! wrapper for things that should be done after solving the update
    virtual void PostUpdate(){};

    //! wrapper for things that should be done before the output
    virtual void PreOutput(){};

    //! wrapper for things that should be done after the output
    virtual void PostOutput(){};
  };
}  // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_AD_STR_LOCA_WRAPPER_H_ */
