/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_projection.H

\brief Fluid field adapter (pressure correction algorithm)

Created on: Jun 16, 2009

<pre>
Maintainer: tw
wiesner@lnm.mw.tum.de
</pre>
*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_PROJECTION_H_
#define ADAPTER_FLUID_PROJECTION_H_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/fluid_projectionmethod.H"
#include "../drt_fluid/fluidresulttest.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "adapter_fluid.H"

namespace ADAPTER
{
/// adapter to pressure correction fluid algorithm
class FluidProjection : public Fluid
{
public:
    FluidProjection(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<LINALG::Solver> solverp, Teuchos::RCP<ParameterList> params, Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond);

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess();
    virtual Teuchos::RCP<const Epetra_Vector> RHS();
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual();
    virtual Teuchos::RCP<const Epetra_Vector> Velnp();
    virtual Teuchos::RCP<const Epetra_Vector> Velaf();
    virtual Teuchos::RCP<const Epetra_Vector> Veln();
    virtual Teuchos::RCP<const Epetra_Vector> Velnm() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accn() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accam();
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Scaam() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Hist();
    virtual Teuchos::RCP<const Epetra_Vector> GridVel() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp();
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() {dserror("not implemented"); return null;};
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel();
    virtual Teuchos::RCP<const Epetra_Vector> FsVel() { dserror("not implemented"); return Teuchos::null; };
    virtual const Teuchos::RCP<Epetra_Vector> StdVeln()  { dserror("not implemented"); return Teuchos::null; };
    virtual const Teuchos::RCP<Epetra_Vector> StdVelnp() { dserror("not implemented"); return Teuchos::null; };
    virtual const Teuchos::RCP<Epetra_Vector> StdVelaf() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap(unsigned nds) { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<LINALG::SparseMatrix>SystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives();
    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual Teuchos::RCP<const DRT::DofSet> DofSet() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();
    /// Return dirichlet boundray condition (used for fluidfluid)
    virtual Teuchos::RCP<const Epetra_Map> FluidDirichMaps(){dserror("not implemented"); return Teuchos::null;};

    virtual FLD::UTILS::FluidXFluidMapExtractor XFluidFluidMapExtractor() {dserror("not implemented");
      // this is a hack to compile the code! danger of dangling reference!
      FLD::UTILS::FluidXFluidMapExtractor  object;
      FLD::UTILS::FluidXFluidMapExtractor& reftoobject=object;
      return reftoobject;}

    virtual Teuchos::RCP<const Epetra_Vector> Stepinc() {dserror("not implemented");return Teuchos::null;}

    virtual void TimeLoop();
    virtual void PrepareTimeStep();
    virtual void PrepareTimeandStep() { dserror("not implemented!"); return; }
    virtual void PrepareSolve() { dserror("not implemented!"); return; }

    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol) { dserror("not implemented!"); return false; }
    virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment) { dserror("not implemented!"); return; }
    virtual void Update();
    virtual void Output();
    virtual void StatisticsOutput() { dserror("not implemented!"); return; }
    virtual IO::DiscretizationWriter& DiscWriter() { return fluid_.DiscWriter();}

    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc);
    virtual void StatisticsAndOutput();

    virtual void NonlinearSolve();
    virtual void Predictor();
    virtual void MultiCorrector();

    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }
    Teuchos::RCP<LINALG::Solver> LinearSolverP() { return solverp_; }

    virtual Teuchos::RCP<const Epetra_Map>    InnerVelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    VelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    PressureRowMap();

    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    virtual double ResidualScaling() const;

    virtual double TimeScaling() const;

    virtual const FLD::UTILS::MapExtractor&  Interface() const { return interface_; }

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

    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno);

    virtual void SetInitialPorosityField(const INPAR::POROELAST::InitialField initfield,const int startfuncno)
    { dserror("not implemented!");return; };

    /// update velocity increment after Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel)  { dserror("not implemented!");return; };

    virtual void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                                   RCP<const Epetra_Vector> scalaram,
                                   RCP<const Epetra_Vector> scalardtam,
                                   RCP<const Epetra_Vector> fsscalaraf,
                                   const double             thermpressaf,
                                   const double             thermpressam,
                                   const double             thermpressdtaf,
                                   const double             thermpressdtam,
                                   Teuchos::RCP<DRT::Discretization> scatradis);

    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                                   const double             thermpressnp,
                                   RCP<const Epetra_Vector> scatraresidual,
                                   Teuchos::RCP<DRT::Discretization> scatradis,
                                   const int whichscalar = 1);

    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp)
    { dserror("not implemented!");return; };

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
       {dserror("not implemented"); return Teuchos::null;}

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln();
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce);
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);
    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");};
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    virtual int  Itemax() const;
    virtual void SetItemax(int itemax);

    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {/*return fluid_.TimIntScheme();*/ dserror("not implemented");return INPAR::FLUID::timeint_gen_alpha;}

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);
    virtual void SetTopOptData(RCP<Epetra_Vector> porosity,RCP<TOPOPT::Optimizer> optimizer){ dserror("not implemented!"); };
    virtual Teuchos::RCP<Epetra_Vector>   IntegrateInterfaceShape();
    virtual void UseBlockMatrix(bool splitmatrix);
    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();
private:
    /// the actual fluid algorithm
    FLD::FluidProjectionMethod			fluid_;

    /// the interface map setup for interface <-> full translation
    FLD::UTILS::MapExtractor                   interface_;

    LINALG::MapExtractor                   	meshmap_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      	dis_;
    Teuchos::RCP<LINALG::Solver>           	solver_;
    Teuchos::RCP<LINALG::Solver>       		solverp_;
    Teuchos::RCP<Teuchos::ParameterList>   	params_;
    Teuchos::RCP<IO::DiscretizationWriter> 	output_;
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RefCountPtr<Epetra_Map>       innervelmap_;
};
}


#endif /* ADAPTER_FLUID_PROJECTION_H_ */
#endif /* CCADISCRET */
