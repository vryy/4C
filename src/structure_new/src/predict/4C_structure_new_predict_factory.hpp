/*-----------------------------------------------------------*/
/*! \file

\brief Factory class to build predictor objects


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_PREDICT_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_PREDICT_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

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
      virtual ~Factory() = default;

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


FOUR_C_NAMESPACE_CLOSE

#endif
