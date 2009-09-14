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
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp();
    virtual Teuchos::RCP<const Epetra_Vector> SgVelVisc();
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel();
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap();
    virtual Teuchos::RCP<LINALG::SparseMatrix>SystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives();
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    virtual void PrepareTimeStep();

    virtual void Update();
    virtual void Output();

    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);
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
    virtual double Time() const;
    virtual int Step() const;
    virtual void LiftDrag();

    virtual void SetInitialFlowField(int whichinitialfield,int startfuncno);

    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                         RCP<const Epetra_Vector> scalarn,
                         RCP<const Epetra_Vector> scalardtn,
                         const double             thermpressnp,
                         const double             thermpressn,
                         const double             thermpressdtn,
                         RCP<const Epetra_Vector> scatraresidual,
                         const int                numscal);

    virtual void SetIterLomaFields(RCP<const Epetra_Vector> scalarnp,
                         RCP<const Epetra_Vector> scalardtnp,
                         const double             thermpressnp,
                         const double             thermpressdtnp,
                         const int                numscal);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln();
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce);
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp);
    virtual void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel);

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    virtual int  Itemax() const;
    virtual void SetItemax(int itemax);

    const FLUID_TIMEINTTYPE TimIntScheme() const {/*return fluid_.TimIntScheme();*/ dserror("not implemented");return timeint_gen_alpha;}

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);
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
