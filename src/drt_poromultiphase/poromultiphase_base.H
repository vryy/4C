/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_base.H

 \brief base class for all porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_BASE_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_BASE_H_


#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_adapter/ad_poromultiphase.H"
#include "../drt_lib/drt_dserror.H"

namespace ADAPTER
{
  class PoroFluidMultiphaseWrapper;
  class Structure;
}

namespace DRT
{
  class Discretization;
}

namespace POROMULTIPHASE
{

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseBase : public ADAPTER::AlgorithmBase, public ADAPTER::PoroMultiPhase
  {
  public:

    /// create using a Epetra_Comm
    PoroMultiPhaseBase(
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
    void PrepareTimeLoop();

    /// timeloop of coupled problem
    void Timeloop();

    /// time step of coupled problem
    virtual void TimeStep() = 0;

    /// prepare time step of coupled problem
    void PrepareTimeStep();

    //! update fields after convergence
    void UpdateAndOutput();

    /// dof map of vector of unknowns of structure field
    virtual Teuchos::RCP<const Epetra_Map> StructDofRowMap() const;

    /// dof map of vector of unknowns of fluid field
    virtual Teuchos::RCP<const Epetra_Map> FluidDofRowMap() const;

    //! access to structural field
    const Teuchos::RCP<ADAPTER::Structure>& StructureField(){return structure_;}

    //! access to fluid field
    const Teuchos::RCP<ADAPTER::PoroFluidMultiphaseWrapper>& FluidField(){return fluid_;}

    /// set structure solution on scatra field
    void SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                            Teuchos::RCP<const Epetra_Vector> vel ) ;

    /// set solid pressure on structure field
    void SetSolidPressure( Teuchos::RCP<const Epetra_Vector> phi ) ;

    /// set scatra solution on fluid field
    void SetScatraSolution(unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars ) ;

    //! setup solver (for monolithic only)
    virtual bool SetupSolver(){return false;};

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

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const
        {dserror("CombinedDBCMap() only available for monolithic schemes!"); return Teuchos::null;};

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> sx,
        Teuchos::RCP<const Epetra_Vector> fx,
        const bool firstcall
        )
        {dserror("Evaluate() only available for monolithic schemes!"); return;};

    //! update all fields after convergence (add increment on displacements and fluid primary variables)
    virtual void UpdateFieldsAfterConvergence(
        Teuchos::RCP<const Epetra_Vector>& sx,
        Teuchos::RCP<const Epetra_Vector>& fx
        )
        {dserror("UpdateFieldsAfterConvergence() only available for monolithic schemes!"); return;};

    //! get monolithic rhs vector
    virtual Teuchos::RCP<const Epetra_Vector> RHS() const
        {dserror("RHS() only available for monolithic schemes!"); return Teuchos::null;};

    //! get extractor
    virtual Teuchos::RCP<const LINALG::MultiMapExtractor> Extractor() const
       {dserror("Extractor() only available for monolithic schemes!"); return Teuchos::null;};

    //! get monolithic block system matrix
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const
       {dserror("BlockSystemMatrix() only available for monolithic schemes!"); return Teuchos::null;};

  private:

    /// set structure mesh displacement on fluid field
    void SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp ) ;

    /// set structure velocity field on fluid field
    void SetVelocityFields( Teuchos::RCP<const Epetra_Vector> vel ) ;

    /// underlying structure of the PoroMultiPhase problem
    Teuchos::RCP< ::ADAPTER::Structure>    structure_;

    /// underlying fluid problem of the PoroMultiPhase problem
    Teuchos::RCP<ADAPTER::PoroFluidMultiphaseWrapper> fluid_;

  protected:

    /// a zero vector of full length of structure dofs
    Teuchos::RCP<Epetra_Vector> zeros_;
    //! here the computation of the structure can be skipped, this is helpful if only fluid-scatra coupling should be calculated
    bool solve_structure_;

    /// Print user output that structure field is disabled
    void PrintStructureDisabledInfo() ;

  }; //PoroMultiPhaseBase


}



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_BASE_H_ */
