/*----------------------------------------------------------------------*/
/*!
\file sti_algorithm.cpp

\brief monolithic algorithm for scatra-thermo interaction

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "sti_algorithm.H"

/*--------------------------------------------------------------------------------*
 | constructor                                                         fang 04/15 |
 *--------------------------------------------------------------------------------*/
STI::Algorithm::Algorithm(
    const Epetra_Comm&              comm,          //! communicator
    const Teuchos::ParameterList&   scatradyn,     //! scalar transport parameter list
    const Teuchos::ParameterList&   solverparams   //! solver parameter list
    ) :
    // instantiate base class
    AlgorithmBase(comm,scatradyn),

    // initialize scatra time integrator
    scatra_(Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,false,"scatra",solverparams))->ScaTraField()),

    // initialize thermo time integrator
    thermo_(Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,false,"thermo",solverparams))->ScaTraField())
{
  return;
}
