/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_PASI_HPP
#define FOUR_C_INPAR_PASI_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | Input parameters for particle structure interaction                       |
 *---------------------------------------------------------------------------*/
namespace Inpar
{
  namespace PaSI
  {
    //! type of partitioned coupling
    enum PartitionedCouplingType
    {
      partitioned_onewaycoup,                 //!< one-way coupling
      partitioned_twowaycoup,                 //!< two-way coupling
      partitioned_twowaycoup_disprelax,       //!< two-way coupling with constant relaxation
      partitioned_twowaycoup_disprelaxaitken  //!< two-way coupling with dynamic aitken relaxation
    };

    //! set valid parameters for particle structure interaction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PaSI

}  // namespace Inpar

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
