/*---------------------------------------------------------------------------*/
/*! \file
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

    //! destructor
    ~PASI_PartOneWayCoup() override{};

    //! setup pasi algorithm
    void Setup() override;

    //! partitioned one way coupled timeloop
    void Timeloop() override;

   protected:
    //! output of fields
    void Output() override;
  };

}  // namespace PASI

/*---------------------------------------------------------------------------*/
#endif
