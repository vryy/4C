/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_base.H

 \brief  base class for all poroelasticity scalar transport interaction algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef PORO_SCATRA_BASE_H_
#define PORO_SCATRA_BASE_H_

#include "../drt_adapter/adapter_algorithmbase.H"

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
  /// base class of algorithms for scalar transport in porous media
  class PORO_SCATRA_Base : public ADAPTER::AlgorithmBase
  {
  public:
    /// create using a Epetra_Comm
    explicit PORO_SCATRA_Base(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);// Problem builder

    /// virtual destructor
    virtual ~PORO_SCATRA_Base(){};

    //! setup discretization, includes cloning the structure discretization
    void SetupDiscretizations(const Epetra_Comm& comm);

    //! Main time loop.
    virtual void Timeloop()=0;

    //! read and set fields needed for restart
    virtual void ReadRestart(int restart);

    //! setup for single fields
    void SetupSystem();

    //! perform result test
    void TestResults(const Epetra_Comm& comm);

  protected:
    //Members

    //! Pointer to the porous media problem. (drt_poroelastic)
    Teuchos::RCP<POROELAST::PoroBase> poro_;
    //! Pointer to the ScaTra problem.     (drt_scatra)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;
  };
}


#endif /* PORO_SCATRA_BASE_H_ */
