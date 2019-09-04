/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
#ifndef INPAR_PASI_H
#define INPAR_PASI_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
#include "inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | Input parameters for particle structure interaction        sfuchs 01/2017 |
 *---------------------------------------------------------------------------*/
namespace INPAR
{
  namespace PASI
  {
    //! type of partitioned coupling
    enum PartitionedCouplingType
    {
      partitioned_onewaycoup,
      partitioned_twowaycoup,
      partitioned_twowaycoup_forcerelax,
      partitioned_twowaycoup_forcerelaxaitken
    };

    //! set valid parameters for particle structure interaction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PASI

}  // namespace INPAR

/*---------------------------------------------------------------------------*/
#endif
