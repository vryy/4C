/*---------------------------------------------------------------------------*/
/*! \file
\brief input parameters for particle structure interaction problems

\level 3

\maintainer Sebastian Fuchs
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
      partitioned_onewaycoup,                 //! one-way structure to particle coupling
      partitioned_twowaycoup,                 //! two-way Dirichlet-Neumann coupling
      partitioned_twowaycoup_disprelax,       //! two-way Dirichlet-Neumann coupling with constant
                                              //! displacement relaxation
      partitioned_twowaycoup_disprelaxaitken  //! two-way Dirichlet-Neumann coupling with dynamic
                                              //! displacement relaxation
    };

    //! set valid parameters for particle structure interaction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PASI

}  // namespace INPAR

/*---------------------------------------------------------------------------*/
#endif
