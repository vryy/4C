/*-----------------------------------------------------------*/
/*! \file
\brief Input parameters for the eXtended contact approach

\maintainer Matthias Mayr

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_INPAR_INPAR_CONTACT_XCONTACT_H_
#define SRC_DRT_INPAR_INPAR_CONTACT_XCONTACT_H_

#include <string.h>
#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace INPAR
{
  namespace XCONTACT
  {
    /// How shall we couple the level-set and the contact simulation?
    enum FieldCouplingAlgorithm
    {
      field_coupl_algo_monolithic,   ///< monolithic coupling algorithm
      field_coupl_algo_partitioned,  ///< partitioned coupling algorithm
      field_coupl_algo_undefined     ///< undefined coupling algorithm
    };

    /// map field coupling scheme enumerator to std::string
    static inline const std::string FieldCouplingAlgorithm2String(
        const enum FieldCouplingAlgorithm& scheme)
    {
      switch (scheme)
      {
        case field_coupl_algo_monolithic:
          return "field_coupl_algo_monolithic";
          break;
        case field_coupl_algo_partitioned:
          return "field_coupl_algo_partitioned";
          break;
        case field_coupl_algo_undefined:
          return "field_coupl_algo_undefined";
          break;
        default:
          return "unknown field coupling algorithm";
          break;
      }
      return "";
    }

    /// set the eXtended contact parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace XCONTACT
}  // namespace INPAR

#endif /* SRC_DRT_INPAR_INPAR_CONTACT_XCONTACT_H_ */
