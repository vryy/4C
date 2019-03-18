/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part.H

 \brief  base class for partitioned poroelasticity scalar transport interaction algorithms

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/

#ifndef PORO_SCATRA_PART_H
#define PORO_SCATRA_PART_H

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include "poro_scatra_base.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
}

namespace POROELAST
{
  class PoroBase;
}

/*----------------------------------------------------------------------*
 |                                                                       |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  /// partitioned algorithm for scalar transport in porous media
  class PoroScatraPart : public PoroScatraBase
  {
   public:
    /// Constructor
    explicit PoroScatraPart(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// destructor
    virtual ~PoroScatraPart(){};

    // Methods

    //! solve one time/incremental step of porous media problem (depending on coupling algorithm)
    virtual void DoPoroStep() = 0;
    //! solve one time/incremental step of scalar transport problem (depending on coupling
    //! algorithm)
    virtual void DoScatraStep() = 0;
  };
}  // namespace POROELAST

#endif
