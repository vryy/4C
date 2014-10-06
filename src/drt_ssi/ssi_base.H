/*!------------------------------------------------------------------------------------------------*
 \file ssi_base.H

 \brief base class for all scalar structure algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_ALGO_H_
#define SSI_ALGO_H_

#include "../drt_adapter/adapter_algorithmbase.H"

#include <Epetra_Vector.h>

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class Structure;
  class ScaTraTimIntImpl;
}

namespace SSI
{

  class SSI_Base : public ADAPTER::AlgorithmBase
  {
  public:

    /// create using a Epetra_Comm
    explicit SSI_Base(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams);// Problem builder

    /// virtual destructor
    virtual ~SSI_Base(){};

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    void TestResults(const Epetra_Comm& comm);

    /// read restart
    void ReadRestart( int restart );

    /// read restart from given time
    void ReadRestartfromTime( double restarttime );

    //! access to structural field
    //const Teuchos::RCP<ADAPTER::Structure>& StructureField(){return structure_;}

    //! access to scatra field
    //const Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>&  ScaTraField(){return scatra_;}

  protected:

    /// underlying structure of the SSI problem
    Teuchos::RCP< ::ADAPTER::Structure>    structure_;

    /// underlying scatra problem of the SSI problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    /// a zero vector of full length
    Teuchos::RCP<Epetra_Vector> zeros_;

  private:
    /// clone scatra discretization form structure discretization
    void SetupDiscretizations(const Epetra_Comm& comm);

    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    //void CheckTimeStepping(double dt1, double dt2);
  };
}

#endif /* SSI_ALGO_H_ */
