/*---------------------------------------------------------------------*/
/*!
\file contact_aug_utils.H

\brief Global utility file for the augmented contact formulation

\level 3

\maintainer Matthias Mayr

\date Mar 23, 2017

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_UTILS_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_UTILS_H_

#include "../drt_inpar/inpar_contact.H"
#include <ext/mt_allocator.h>

#include "../headers/pairedmatrix.H"

namespace GEN_DATA = GEN;

//#define INSERT_AND_SORT
#define QUICK_INSERT


namespace LINALG
{
  class Solver;
}  // namespace LINALG
namespace CONTACT
{
  class CoInterface;
  class CoAbstractStrategy;
  namespace AUG
  {
    typedef std::map<int, double> plain_double_map;
    typedef std::map<int, plain_double_map> plain_map_map;

    typedef std::vector<Teuchos::RCP<CONTACT::CoAbstractStrategy>> plain_strategy_set;
    typedef std::vector<Teuchos::RCP<CONTACT::CoInterface>> plain_interface_set;
    typedef std::vector<Teuchos::RCP<LINALG::Solver>> plain_lin_solver_set;
    typedef std::vector<plain_interface_set> plain_interface_sets;
    typedef std::vector<INPAR::CONTACT::SolvingStrategy> plain_strattype_set;

#if defined(INSERT_AND_SORT)
    typedef GEN_DATA::pairedvector<int, double, GEN_DATA::insert_and_sort_policy<int, double>>
        Deriv1stMap;
    typedef GEN_DATA::pairedmatrix<int, double, GEN_DATA::insert_and_sort_policy<int, double>>
        Deriv2ndMap;
#elif defined(QUICK_INSERT)
    typedef GEN_DATA::quick_pairedvector<int, double> Deriv1stMap;
    typedef GEN_DATA::quick_pairedmatrix<int, double> Deriv2ndMap;
#else
    typedef GEN_DATA::default_pairedvector<int, double> Deriv1stMap;
    typedef GEN_DATA::default_pairedmatrix<int, double> Deriv2ndMap;
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

  }  // namespace AUG
}  // namespace CONTACT

#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_UTILS_H_ */
