/*!----------------------------------------------------------------------
\file pasi_partitioned_onewaycoup.H

\brief one way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                           sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_ONEWAYCOUP_H
#define PASI_PARTITIONED_ONEWAYCOUP_H

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
#include "pasi_partitioned.H"

/*----------------------------------------------------------------------*
 | one way coupled pasi algorithm                        sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
namespace PASI
{

//! one way coupled partitioned algorithm for particle structure interaction
class PASI_PartOneWayCoup: public PartitionedAlgo
{

public:

  //! constructor
  explicit PASI_PartOneWayCoup(
      const Epetra_Comm&              comm,         //! communicator
      const Teuchos::ParameterList&   pasi_params   //! input parameters for particle structure interaction
      );

  //! virtual destructor
  virtual ~PASI_PartOneWayCoup(){};

  //! timeloop of coupled problem
  virtual void Timeloop();

}; // class PASI_PartOneWayCoup: public PartitionedAlgo

} // namespace PASI

/*----------------------------------------------------------------------*/
#endif // PASI_PARTITIONED_ONEWAYCOUP_H
