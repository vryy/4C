/*!----------------------------------------------------------------------
 \file poro_scatra.H
 \brief

 <pre>
 Maintainer: Anh-Tu Vuong
 vuong@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15264
 </pre>

 *----------------------------------------------------------------------*/

#ifndef PORO_SCATRA_PART_H
#define PORO_SCATRA_PART_H

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include "poro_scatra_base.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
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
  class PORO_SCATRA_Part: public PORO_SCATRA_Base
  {
  public:

    ///Constructor
    explicit PORO_SCATRA_Part(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    ///destructor
    virtual ~PORO_SCATRA_Part()
    {
    }
    ;

    //Methods

    //! solve one time step of porous media problem
    void DoPoroStep();
    //! solve one time step of scalar transport problem
    void DoScatraStep();

    //! apply solution of poro-problem to scatra
    void SetPoroSolution();
    //! apply displacement fields to scatra
    void SetMeshDisp();
    //! apply velocity fields to scatra
    void SetVelocityFields();

    //! Main time loop.
    virtual void Timeloop()=0;
  };
}

#endif
