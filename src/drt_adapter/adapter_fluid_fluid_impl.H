/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_fluid_impl.H


<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_FLUID_IMPL_H_
#define ADAPTER_FLUID_FLUID_IMPL_H_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/xfluidfluid.H"
#include "../drt_fluid/xfluidfluidresulttest.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "adapter_fluid.H"

namespace ADAPTER {
  class FluidFluidImpl : public Fluid
  {
  public:
   FluidFluidImpl(Teuchos::RCP<DRT::Discretization> embfluiddis,
             Teuchos::RCP<DRT::Discretization> bgfluiddis,
             Teuchos::RCP<LINALG::Solver> solver,
             Teuchos::RCP<ParameterList> params,
             bool isale,
             bool dirichletcond,
             bool monolithicfluidfluidfsi);


    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess();
    virtual Teuchos::RCP<const Epetra_Vector> RHS();
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual();
    virtual Teuchos::RCP<const Epetra_Vector> Velnp();
    virtual Teuchos::RCP<const Epetra_Vector> Velaf();
    virtual Teuchos::RCP<const Epetra_Vector> Veln();
    virtual Teuchos::RCP<const Epetra_Vector> Velnm() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accn() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accam() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Scaam() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> GridVel() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<const Epetra_Vector> Hist(){ dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp();
    virtual Teuchos::RCP<const Epetra_Vector> Dispn(){dserror("not implemented"); return null;}
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel();
    virtual Teuchos::RCP<const Epetra_Vector> FsVel(){ dserror("not implemented"); return Teuchos::null; };
    virtual const Teuchos::RCP<Epetra_Vector> StdVeln()  { dserror("not implemented"); return Teuchos::null; };
    virtual const Teuchos::RCP<Epetra_Vector> StdVelnp() { dserror("not implemented"); return Teuchos::null; };
    virtual const Teuchos::RCP<Epetra_Vector> StdVelaf() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap(unsigned nds) { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<LINALG::SparseMatrix>SystemMatrix() { return fluid_.SystemMatrix(); }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() { return fluid_.ShapeDerivatives(); }
    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual Teuchos::RCP<const DRT::DofSet> DofSet() {dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();
    virtual Teuchos::RCP<const Epetra_Map>  FluidDirichMaps() {return fluid_.FluidFluidDirichMaps();}

    virtual void TimeLoop();
    virtual void PrepareTimeStep();
    virtual void PrepareTimeandStep() { dserror("not implemented!"); return; }
    virtual void PrepareSolve() { dserror("not implemented!"); return; }
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);
    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol) { dserror("not implemented!"); return false; }
    virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment) { dserror("not implemented!"); return; }
    virtual void Update();
    virtual void StatisticsAndOutput();
    virtual void Output();
    virtual void StatisticsOutput() { dserror("not implemented!"); return; }
    virtual IO::DiscretizationWriter& DiscWriter() {return fluid_.DiscWriter();}

    virtual void NonlinearSolve();
    virtual void Predictor() {dserror("not implemented!");}
    virtual void MultiCorrector() {dserror("not implemented!");}


    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    virtual Teuchos::RCP<const Epetra_Map>    InnerVelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    VelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    PressureRowMap();

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    /// scaling factor needed to convert the residual to real forces
    virtual double ResidualScaling() const;

    virtual double TimeScaling() const;

    /// communication object at the interface
    virtual const FLD::UTILS::MapExtractor& Interface() const { return interface_; }

    virtual void ReadRestart(int step);
    virtual void SetRestart(const int step,
                    const double time,
                    Teuchos::RCP<const Epetra_Vector> readvelnp,
                    Teuchos::RCP<const Epetra_Vector> readveln,
                    Teuchos::RCP<const Epetra_Vector> readvelnm,
                    Teuchos::RCP<const Epetra_Vector> readaccnp,
                    Teuchos::RCP<const Epetra_Vector> readaccn) { dserror("not implemented!"); return; }
    virtual double Time() const;
    virtual int Step() const;
    virtual double Dt() const;
    virtual void LiftDrag();
    virtual const double EvalTime() const { dserror("not implemented!"); return 0.0; }
    virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph) { dserror("not implemented!"); return; }
    virtual void TransferVectorsToNewDistribution(const Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle) { dserror("not implemented!"); return; }

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin() {dserror("not implemented!"); return null;}
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln()  {dserror("not implemented!"); return Teuchos::null;}
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce) {dserror("not implemented");}
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);
    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");};
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx) {dserror("not implemented");}
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx) {dserror("not implemented");}

    virtual FLD::UTILS::FluidXFluidMapExtractor XFluidFluidMapExtractor() {return fluid_.XFluidFluidMapExtractor();}

    virtual Teuchos::RCP<const Epetra_Vector> Stepinc() {return fluid_.Stepinc();}

    virtual int  Itemax() const;
    virtual void SetItemax(int itemax);

    /*!
    \brief return type of time integration scheme

    */
    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return fluid_.TimIntScheme();}

    virtual Teuchos::RCP<Epetra_Vector>   IntegrateInterfaceShape();

    virtual void UseBlockMatrix(bool splitmatrix);

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    virtual Teuchos::RCP<Epetra_Vector>    RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) {dserror("not implemented"); return null;}

    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();

    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

    virtual void SetTopOptData(RCP<Epetra_Vector> porosity,RCP<TOPOPT::Optimizer> optimizer)
    { dserror("not implemented!"); };

    /// set initial flow field
    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno);

    virtual void SetInitialPorosityField(const INPAR::POROELAST::InitialField initfield,const int startfuncno)
    { dserror("not implemented!");return; };

    /// update velocity increment after Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel)  { dserror("not implemented!");return; };

    ///  set fields for low-Mach-number flow within outer iteration loop
    virtual void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                                   RCP<const Epetra_Vector> scalaram,
                                   RCP<const Epetra_Vector> scalardtam,
                                   const double             thermpressaf,
                                   const double             thermpressam,
                                   const double             thermpressdtaf,
                                   const double             thermpressdtam,
                                   Teuchos::RCP<DRT::Discretization> scatradis)
    {dserror("not implemented");}

    /// set fields for low-Mach-number flow at end of time step
    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                                   const double             thermpressnp,
                                   RCP<const Epetra_Vector> scatraresidual,
                                   Teuchos::RCP<DRT::Discretization> scatradis,
                                   const int whichscalar = -1)
    {dserror("not implemented");};

    /// set velocity field (separate computation)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp)
    { dserror("not implemented!");return; };

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
        {dserror("not implemented"); return Teuchos::null;}


  private:

    /// the actual fluid algorithm
    FLD::XFluidFluid        fluid_;

    /// the interface map setup for fsi interface, free surface, interior translation
    FLD::UTILS::MapExtractor               interface_;

    /// ALE dof map
    LINALG::MapExtractor                   meshmap_;

    // FSI-Cond Map of ALE (used for monolithic-fluidfluid-fsi)
    LINALG::MapExtractor                   fsiale_;

    //! @name local copies of input parameters
    // fluid discretization
    Teuchos::RCP<DRT::Discretization>      embfluiddis_;
    // xfluid discretization
    Teuchos::RCP<DRT::Discretization>      bgfluiddis_;
//     // fluid-xfluid-boundary discretization
//     Teuchos::RCP<DRT::Discretization>      fluidxfluidboundarydis_;

//     Teuchos::RCP<IO::DiscretizationWriter> fluidxfluidboundaryoutput_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RefCountPtr<Epetra_Map>       innervelmap_;

    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_; // xfluidfluid monolithic approach

    bool monolithicfluidfluidfsi_;
  };
}

#endif /* ADAPTER_FLUID_FLUID_IMPL_H_ */
#endif
