/*---------------------------------------------------------------------------*/
/*!
\brief one way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_ONEWAYCOUP_H
#define PASI_PARTITIONED_ONEWAYCOUP_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
namespace PASI
{
  //! one way coupled partitioned algorithm
  class PASI_PartOneWayCoup : public PartitionedAlgo
  {
   public:
    //! constructor
    explicit PASI_PartOneWayCoup(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~PASI_PartOneWayCoup(){};

    //! setup pasi algorithm
    virtual void Setup() override;

    //! partitioned one way coupled timeloop
    virtual void Timeloop() override;
  };

}  // namespace PASI

/*---------------------------------------------------------------------------*/
#endif
