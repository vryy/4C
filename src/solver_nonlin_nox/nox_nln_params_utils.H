/*----------------------------------------------------------------------------*/
/*!
 * \file nox_nln_params_utils.H
 *
 * \brief Utility methods to handle and parse (input) parameter lists
 *
 * \level 3
 *
 * \maintainer Michael Hiermeier
 *
 * \date Feb 13, 2018
 */
/*----------------------------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_PARAMS_UTILS_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_PARAMS_UTILS_H_

#include <Teuchos_StandardParameterEntryValidators.hpp>

namespace NOX
{
  namespace NLN
  {
    namespace PUTILS
    {
      /// \brief set a new set of validate integral parameters in the parameter list
      /** This function is a wrapper for the Teuchos method. The main difference
       *  is that the parameter-list contains already the input parameter and,
       *  therefore, it is necessary to store a copy before the validator can
       *  be set-up. Subsequently, the input parameter is tested and if valid,
       *  the corresponding integral value will be returned.
       *
       *  \author hiermeier \date 02/18 */
      template <typename IntegralType>
      IntegralType SetAndValidate(Teuchos::ParameterList& p, const std::string& param_name,
          const std::string& default_value, const std::string& documentation,
          const Teuchos::ArrayView<std::string>& value_list,
          const Teuchos::ArrayView<IntegralType>& ivalue_list)
      {
        // store copy of the input value
        const std::string input_value(p.get<std::string>(param_name, default_value));

        // setup validator
        Teuchos::setStringToIntegralParameter<IntegralType>(
            param_name, default_value, documentation, value_list, ivalue_list, &p);

        // set and validate input value
        p.set(param_name, input_value);

        // convert the input value and return
        return Teuchos::getIntegralValue<IntegralType>(p, param_name);
      }
    }  // namespace PUTILS
  }    // namespace NLN
}  // namespace NOX

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_PARAMS_UTILS_H_ */
