/*---------------------------------------------------------------------*/
/*! \file
\brief Global utility file for the augmented contact formulation

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_UTILS_HPP
#define FOUR_C_CONTACT_AUG_UTILS_HPP

#include "4C_config.hpp"

#include "4C_inpar_contact.hpp"
#include "4C_utils_pairedmatrix.hpp"

#include <ext/mt_allocator.h>

// #define INSERT_AND_SORT
#define QUICK_INSERT

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Solver;
}  // namespace Core::LinAlg
namespace CONTACT
{
  class Interface;
  class AbstractStrategy;
  namespace Aug
  {
    typedef std::map<int, double> plain_double_map;
    typedef std::map<int, plain_double_map> plain_map_map;

    typedef std::vector<Teuchos::RCP<CONTACT::AbstractStrategy>> plain_strategy_set;
    typedef std::vector<Teuchos::RCP<CONTACT::Interface>> plain_interface_set;
    typedef std::vector<Teuchos::RCP<Core::LinAlg::Solver>> plain_lin_solver_set;
    typedef std::vector<plain_interface_set> plain_interface_sets;
    typedef std::vector<FourC::Inpar::CONTACT::SolvingStrategy> plain_strattype_set;

#if defined(INSERT_AND_SORT)
    typedef Core::Gen::Pairedvector<int, double, Gen::insert_and_sort_policy<int, double>>
        Deriv1stMap;
    typedef Core::Gen::pairedmatrix<int, double, Gen::insert_and_sort_policy<int, double>>
        Deriv2ndMap;
#elif defined(QUICK_INSERT)
    typedef Core::Gen::quick_pairedvector<int, double> Deriv1stMap;
    typedef Core::Gen::quick_pairedmatrix<int, double> Deriv2ndMap;
#else
    typedef Gen::default_pairedvector<int, double> Deriv1stMap;
    typedef Gen::default_pairedmatrix<int, double> Deriv2ndMap;
#endif

    typedef std::vector<Deriv1stMap> Deriv1stVecMap;
    typedef std::vector<Deriv2ndMap> Deriv2ndVecMap;

    inline std::string contact_func_name(const std::string& pretty_func_name)
    {
      // skip return parameters
      const unsigned b = pretty_func_name.find("CONTACT", 0);
      // find brackets for input parameters
      const unsigned e_1 = pretty_func_name.find("(", b);
      // find brackets for template parameters
      const unsigned e_2 = pretty_func_name.substr(0, e_1).find("<", b);

      return pretty_func_name.substr(b, std::min(e_1, e_2) - b);
    }

/** macro which returns the full function name (incl. namespaces) of any
 *  function starting with the namespace CONTACT
 *
 *  \author hiermeier \date 04/17 */
#define CONTACT_FUNC_NAME contact_func_name(__PRETTY_FUNCTION__)

  }  // namespace Aug
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
