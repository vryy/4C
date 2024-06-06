/*----------------------------------------------------------------------*/
/*! \file

 \brief  base class for partitioned poroelasticity scalar transport interaction algorithms

\level 2

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROELAST_SCATRA_PART_HPP
#define FOUR_C_POROELAST_SCATRA_PART_HPP

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_poroelast_scatra_base.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace Adapter
{
  class ScaTraBaseAlgorithm;
}

namespace PoroElast
{
  class PoroBase;
}

/*----------------------------------------------------------------------*
 |                                                                       |
 *----------------------------------------------------------------------*/
namespace PoroElastScaTra
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
    virtual void do_scatra_step() = 0;
  };
}  // namespace PoroElastScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
