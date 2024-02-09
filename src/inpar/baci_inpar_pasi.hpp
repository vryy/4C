/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_INPAR_PASI_HPP
#define BACI_INPAR_PASI_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_inpar_parameterlist_utils.hpp"

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | Input parameters for particle structure interaction                       |
 *---------------------------------------------------------------------------*/
namespace INPAR
{
  namespace PASI
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

  }  // namespace PASI

}  // namespace INPAR

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
