/*!----------------------------------------------------------------------
\file fpsi_paritioned.H


<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de

</pre>

*----------------------------------------------------------------------*/


#ifndef FPSI_PARTITIONED_H
#define FPSI_PARTITIONED_H

#include "fpsi.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_adapter/ad_fld_fluid_ale.H"
#include "../drt_poroelast/poro_base.H"
#include "../drt_adapter/adapter_coupling.H"

namespace FPSI
{
  class Partitioned : public FPSI_Base
  {
  public:

    //ctor
    Partitioned(const Epetra_Comm& comm,
                const Teuchos::ParameterList& timeparams);

    //dtor
    virtual ~Partitioned(){};

//////////////// FPSI_DYN
//
    void ReadRestart(int restartstep);
//
    // redistribute FPSI interface if running on parallel
    void RedistributeInterface();

    void SetupSystem();
//

//
      ////////// TIMELOOP
      void Timeloop();
//    //
      void InitialCalculations();
//    //
      void IncrementTimeAndStep();
//    //
      void PrepareTimeStep();
//    //
      void OuterLoop();
//    //
      void TimeUpdateAndOutput();
//    //

    void TestResults(const Epetra_Comm& comm);

    bool ConvergenceCheck(int itnum);
//
///////////////// END ALGORITHM


    bool matchingnodes() {return matchingnodes_; }

    void SetDefaultParameters(const Teuchos::ParameterList& fsidyn);

    Teuchos::RCP<Epetra_Vector> InterfaceDisp();

    Teuchos::RCP<Epetra_Vector> InterfaceForce();

    void InterfaceVelocity(Teuchos::RCP<const Epetra_Vector> idispnp);


    const Teuchos::RCP<POROELAST::PoroBase>& PoroField(){return poroelast_subproblem_;}
    const Teuchos::RCP<ADAPTER::FluidMovingBoundary>& FluidField(){return fluid_subproblem_;}

  protected:

    int step_;
    int numstep_;

    double dt_;
    double time_;
    double timemax_;

    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;

  private:

    Teuchos::RCP<POROELAST::PoroBase>             poroelast_subproblem_;
    Teuchos::RCP< ::ADAPTER::FluidMovingBoundary> fluid_subproblem_;

    double  ittol_;

    bool matchingnodes_;
    Teuchos::RCP<ADAPTER::Coupling>  coupsf_;

    /// interface displacement and velocity from time step begin (d^n and v^n)
    Teuchos::RCP<Epetra_Vector> idispn_;
    Teuchos::RCP<Epetra_Vector> iveln_;

  };
}

#endif // FPSI_PARTITIONED
