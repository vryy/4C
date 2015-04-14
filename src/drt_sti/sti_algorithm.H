/*----------------------------------------------------------------------*/
/*!
\file sti_algorithm.H

\brief monolithic algorithm for scatra-thermo interaction

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STI_ALGORITHM_H
#define STI_ALGORITHM_H

#include "../drt_adapter/adapter_algorithmbase.H"

// forward declaration
namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace STI
{
  //! monolithic algorithm for scatra-thermo interaction
  class Algorithm : public ADAPTER::AlgorithmBase
  {
    public:

      //! constructor
      explicit Algorithm(
          const Epetra_Comm&              comm,          //! communicator
          const Teuchos::ParameterList&   scatradyn,     //! scalar transport parameter list
          const Teuchos::ParameterList&   solverparams   //! solver parameter list
          );

      //! read restart data
      void ReadRestart(
          int step
          )
      {return;};

      //! time loop
      void TimeLoop(){return;};

      //! access scatra time integrator
      const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Scatra() const {return scatra_;}

      //! access thermo time integrator
      const Teuchos::RCP<SCATRA::ScaTraTimIntImpl>& Thermo() const {return thermo_;}

    private:

      //! scatra time integrator
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

      //! thermo time integrator
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> thermo_;
  }; // class Algorithm : public ADAPTER::AlgorithmBase
} // namespace STI
#endif
