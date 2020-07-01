/*-----------------------------------------------------------*/
/*! \file

\brief Factory class to build predictor objects


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_PREDICT_FACTORY_H_
#define STR_PREDICT_FACTORY_H_

#include "../drt_inpar/inpar_structure.H"

// forward declaration
namespace STR
{
  namespace PREDICT
  {
    class Generic;

    /*! \brief Factory to build the desired predictor
     *
     *  \author Michael Hiermeier */
    class Factory
    {
     public:
      //! constructor
      Factory();

      //! destructor
      virtual ~Factory(){};

      //! build the desired predictor
      Teuchos::RCP<STR::PREDICT::Generic> BuildPredictor(
          const enum INPAR::STR::PredEnum& predType) const;
    };

    /*! \brief Non-member function, which relates to the STR::PREDICT::Factory class
     *
     * \note Call this method from outside!
     */
    Teuchos::RCP<STR::PREDICT::Generic> BuildPredictor(const enum INPAR::STR::PredEnum& predType);

  }  // namespace PREDICT
}  // namespace STR


#endif /* STR_PREDICT_FACTORY_H_ */
