/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_base.H

 \brief  base class for all poroelasticity scalar transport interaction algorithms

 \level 2

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15251
 *----------------------------------------------------------------------*/


#ifndef PORO_SCATRA_BASE_H_
#define PORO_SCATRA_BASE_H_

#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_poroelast/poro_base.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class MortarVolCoupl;
  class FluidPoro;
  class FPSIStructureWrapper;
}

//namespace POROELAST
//{
//  class PoroBase;
//}

/*----------------------------------------------------------------------*
 |                                                                       |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  /// base class of algorithms for scalar transport in porous media
  class PoroScatraBase : public ADAPTER::AlgorithmBase
  {
  public:
    /// create using a Epetra_Comm
    explicit PoroScatraBase(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);// Problem builder

    /// virtual destructor
    virtual ~PoroScatraBase(){};

    //! Main time loop.
    virtual void Timeloop()=0;

    //! prepare time step for single fields
    virtual void PrepareTimeStep(bool printheader = true){dserror("not implemented in base class. override in subclass.");};

    //! perform iteration loop between fields
    virtual void Solve(){dserror("not implemented in base class. override in subclass.");};;

    //! prepare output
    virtual void PrepareOutput(){dserror("not implemented in base class. override in subclass.");};;

    //! update time step
    virtual void Update(){dserror("not implemented in base class. override in subclass.");};;

    //! write output print to screen
    virtual void Output(){dserror("not implemented in base class. override in subclass.");};;

    //! read and set fields needed for restart
    virtual void ReadRestart(int restart)=0;

    //! setup for single fields
    virtual void SetupSystem();

    //! Build the combined dirichlet map of the monolithic poro problem
    virtual void BuildCombinedDBCMap(){poro_->BuildCombinedDBCMap();};

    //! perform result test
    void TestResults(const Epetra_Comm& comm);

    //! apply solution of poro-problem to scatra
    void SetPoroSolution();

    //! apply solution of scatra to poro
    void SetScatraSolution();

    //! return pointer to porous medium problem
    const Teuchos::RCP<POROELAST::PoroBase>& PoroField(){return poro_;};

    //! return pointer to interstitial fluid
    const Teuchos::RCP<ADAPTER::FluidPoro>& FluidField(){return poro_->FluidField();};

    //! return pointer to porous structure
    const Teuchos::RCP<ADAPTER::FPSIStructureWrapper>& StructureField(){return poro_->StructureField();};

    //! return pointer to scalar transport problem
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() { return scatra_->ScaTraField(); };

    //! return pointer to scalar problem adapter base class
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> ScaTraFieldBase() { return scatra_; };

    //! setup solver (for monolithic only)
    virtual bool SetupSolver(){return true;};

  protected:
    //!setup up of dofsets for two way coupling
    void ReplaceDofSets(Teuchos::RCP<DRT::Discretization> structdis,
                    Teuchos::RCP<DRT::Discretization> fluiddis,
                    Teuchos::RCP<DRT::Discretization> scatradis);

    //!setup up coupling objects if necessary
    void SetupCoupling(Teuchos::RCP<DRT::Discretization> structdis,
                       Teuchos::RCP<DRT::Discretization> fluiddis,
                       Teuchos::RCP<DRT::Discretization> scatradis);

    //! Pointer to the porous media problem. (drt_poroelastic)
    Teuchos::RCP<POROELAST::PoroBase> poro_;
    //! Pointer to the ScaTra problem.     (drt_scatra)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    //! @name Volume Mortar stuff

    //! flag for matchinggrid
    const bool matchinggrid_;

    //! volume coupling (using mortar) adapter
    Teuchos::RCP<ADAPTER::MortarVolCoupl> volcoupl_structurescatra_;
    Teuchos::RCP<ADAPTER::MortarVolCoupl> volcoupl_fluidscatra_;
    //@}

  private:
    //! apply displacement fields to scatra
    void SetMeshDisp();
    //! apply velocity fields to scatra
    void SetVelocityFields();

  };
}


#endif /* PORO_SCATRA_BASE_H_ */
