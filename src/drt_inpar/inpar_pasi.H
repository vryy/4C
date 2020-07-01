/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef INPAR_PASI_H
#define INPAR_PASI_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "inpar_parameterlist_utils.H"

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
#endif
