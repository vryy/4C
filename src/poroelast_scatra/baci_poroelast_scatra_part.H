/*----------------------------------------------------------------------*/
/*! \file

 \brief  base class for partitioned poroelasticity scalar transport interaction algorithms

\level 2

 *----------------------------------------------------------------------*/

#ifndef BACI_POROELAST_SCATRA_PART_H
#define BACI_POROELAST_SCATRA_PART_H

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include "baci_config.H"

#include "baci_poroelast_scatra_base.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

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
namespace POROELASTSCATRA
{
  /// partitioned algorithm for scalar transport in porous media
  class PoroScatraPart : public PoroScatraBase
  {
   public:
    /// Constructor
    explicit PoroScatraPart(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    // Methods

    //! solve one time/incremental step of porous media problem (depending on coupling algorithm)
    virtual void DoPoroStep() = 0;
    //! solve one time/incremental step of scalar transport problem (depending on coupling
    //! algorithm)
    virtual void DoScatraStep() = 0;
  };
}  // namespace POROELASTSCATRA

BACI_NAMESPACE_CLOSE

#endif
