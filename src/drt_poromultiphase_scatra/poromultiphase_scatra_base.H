/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_base.H

 \brief base algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_BASE_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_BASE_H_


#include "../drt_adapter/adapter_algorithmbase.H"

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
}

namespace DRT
{
  class Discretization;
}

namespace POROMULTIPHASE
{
 class PoroMultiPhaseBase;
}

namespace POROMULTIPHASESCATRA
{

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseScaTraBase : public ADAPTER::AlgorithmBase
  {
  public:

    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraBase(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);// Problem builder

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraBase(){};

    /// initialization
    virtual void Init(
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams,
        const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams,
        const std::string& struct_disname,
        const std::string& fluid_disname,
        const std::string& scatra_disname,
        bool isale,
        int nds_disp,
        int nds_vel,
        int nds_solidpressure) =0;

    /// read restart
    void ReadRestart( int restart );

    /// create result test for subproblems
    void CreateFieldTest();

    /// setup
    virtual void SetupSystem() = 0;

    /// prepare timeloop of coupled problem
    virtual void PrepareTimeLoop() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// time step of coupled problem
    virtual void TimeStep() = 0;

    /// time step of coupled problem
    virtual void PrepareTimeStep() = 0;

    //! update time step and print to screen
    virtual void UpdateAndOutput() = 0;

    //! apply solution of poro-problem to scatra
    void SetPoroSolution();

    //! apply solution of scatra to poro
    void SetScatraSolution();

  protected:
    /// underlying poroelast multi phase
    Teuchos::RCP<POROMULTIPHASE::PoroMultiPhaseBase> poromulti_;

    /// underlying scatra problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

  private:

    //! apply displacement fields to scatra
    void SetMeshDisp();
    //! apply solution fields (except displacements) to scatra
    void SetSolutionFields();

  }; //PoroMultiPhaseScaTraBase


}



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_BASE_H_ */
