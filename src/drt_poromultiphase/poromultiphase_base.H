/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_base.H

 \brief base class for all porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_BASE_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_BASE_H_


#include "../drt_adapter/adapter_algorithmbase.H"

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>

namespace ADAPTER
{
  class PoroFluidMultiphase;
  class Structure;
}

namespace DRT
{
  class Discretization;
}

namespace POROMULTIPHASE
{

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseBase : public ADAPTER::AlgorithmBase
  {
  public:

    /// create using a Epetra_Comm
    explicit PoroMultiPhaseBase(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);// Problem builder

    /// virtual destructor
    virtual ~PoroMultiPhaseBase(){};

    /// initialization
    virtual void Init(
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams,
        const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams,
        const std::string& struct_disname,
        const std::string& fluid_disname,
        bool isale,
        int nds_disp,
        int nds_vel,
        int nds_solidpressure,
        int ndsporofluid_scatra) =0;

    /// read restart
    void ReadRestart( int restart );

    /// test results (if necessary)
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

    /// set structure solution on scatra field
    void SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                            Teuchos::RCP<const Epetra_Vector> vel ) ;

    /// set solid pressure on structure field
    void SetSolidPressure( Teuchos::RCP<const Epetra_Vector> phi ) ;

    /// set scatra solution on fluid field
    void SetScatraSolution(unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars ) ;

    /// unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> StructDispnp() const;

    /// unknown velocity at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> StructVelnp() const;

    /// return fluid flux
    Teuchos::RCP<const Epetra_MultiVector> FluidFlux() const;

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> FluidPhinp() const;

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> FluidSaturation() const;

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> FluidPressure() const;

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> SolidPressure() const;

  private:

    /// set structure mesh displacement on fluid field
    void SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp ) ;

    /// set structure velocity field on fluid field
    void SetVelocityFields( Teuchos::RCP<const Epetra_Vector> vel ) ;

  protected:
    /// underlying structure of the PoroMultiPhase problem
    Teuchos::RCP< ::ADAPTER::Structure>    structure_;

    /// underlying fluid problem of the PoroMultiPhase problem
    Teuchos::RCP<ADAPTER::PoroFluidMultiphase> fluid_;

    /// a zero vector of full length of structure dofs
    Teuchos::RCP<Epetra_Vector> zeros_;

  }; //PoroMultiPhaseBase


}



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_BASE_H_ */
