// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_NEW_PREDICT_FACTORY_HPP
#define FOUR_C_STRUCTURE_NEW_PREDICT_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Solid
{
  namespace Predict
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
      Teuchos::RCP<Solid::Predict::Generic> build_predictor(
          const enum Inpar::Solid::PredEnum& predType) const;
    };

    /*! \brief Non-member function, which relates to the Solid::Predict::Factory class
     *
     * \note Call this method from outside!
     */
    Teuchos::RCP<Solid::Predict::Generic> build_predictor(
        const enum Inpar::Solid::PredEnum& predType);

  }  // namespace Predict
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
