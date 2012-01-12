/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_impl.H

\brief Fluid field adapter

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_IMPL_H
#define ADAPTER_FLUID_IMPL_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/fluidimplicitintegration.H"

#include "adapter_fluid.H"



// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// adapter to normal fluid algorithm
  /*!
    <h3>Scaling factors in monolithic FSI</h3>

    For monolithic FSI we need to build one block system with overlapping
    equations at the interface. In order to do so we have to scale all blocks
    such that we the terms we add are of the same units. We go for real forces
    and displacement increments at the interface.

    It is the duty of the adapter to provide the necessary scaling
    factors.

    We assume a very simple velocity -- displacement relation at the interface
    \f$ \mathbf{u}^{n+1}_{i+1} = \frac{1}{\Delta t} \Delta
    \mathbf{d}^{n+1}_{i+1} \f$

    The fluid interface equations have to be scaled with \f$
    \frac{\rho}{\theta\Delta t} \f$ to be physical forces. The matrix columns
    that go with \f$ \mathbf{u}^{n+1}_{\Gamma,i+1} \f$ have to be scaled with
    \f$ \frac{1}{\Delta t} \f$. This results in the scaled blocks:

    - \f$ \frac{\rho}{\theta\Delta t^2} \mathbf{F}_{\Gamma\Gamma} \f$

    - \f$ \frac{1}{\Delta t} \mathbf{F}_{I\Gamma} \f$

    - \f$ \frac{\rho}{\theta\Delta t} \mathbf{F}_{\Gamma I} \f$

    \note As we couple a generalized-\f$\alpha\f$ structure to a
    one-step-\f$\theta\f$ fluid here all we can do is backward Euler time
    discretization on both fields.

    \author u.kue
    \date 02/08
   */
  class FluidImpl : public Fluid
  {
  public:
    FluidImpl(Teuchos::RCP<DRT::Discretization> dis,
              Teuchos::RCP<LINALG::Solver> solver,
              Teuchos::RCP<ParameterList> params,
              Teuchos::RCP<IO::DiscretizationWriter> output,
              bool isale,
              bool dirichletcond);

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess();
    virtual Teuchos::RCP<const Epetra_Vector> RHS();
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual();
    virtual Teuchos::RCP<const Epetra_Vector> Velnp();
    virtual Teuchos::RCP<const Epetra_Vector> Velaf();
    virtual Teuchos::RCP<const Epetra_Vector> Veln();
    virtual Teuchos::RCP<const Epetra_Vector> Velnm();
    virtual Teuchos::RCP<const Epetra_Vector> Accnp();
    virtual Teuchos::RCP<const Epetra_Vector> Accn();
    virtual Teuchos::RCP<const Epetra_Vector> Accam();
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf();
    virtual Teuchos::RCP<const Epetra_Vector> Scaam();
    virtual Teuchos::RCP<const Epetra_Vector> Hist();
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp();
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel();
    virtual Teuchos::RCP<const Epetra_Vector> FsVel();
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap(unsigned nds);
    virtual Teuchos::RCP<LINALG::SparseMatrix>SystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives();
    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual Teuchos::RCP<const DRT::DofSet> DofSet() { dserror("not implemented"); return Teuchos::null; };
    virtual FLD::UTILS::FluidXFluidMapExtractor XFluidFluidMapExtractor() {dserror("not implemented");
      // this is a hack to compile the code! danger of dangling reference!
      FLD::UTILS::FluidXFluidMapExtractor  object;
      FLD::UTILS::FluidXFluidMapExtractor& reftoobject=object;
      return reftoobject;}
    virtual Teuchos::RCP<const Epetra_Vector> Stepinc() {dserror("not implemented");return Teuchos::null;}
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();
    /// Return dirichlet boundray condition (used for fluidfluid)
    virtual Teuchos::RCP<const Epetra_Map> FluidDirichMaps(){dserror("not implemented"); return Teuchos::null;};

    //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const;

    virtual void TimeLoop();
    virtual void PrepareTimeStep();
    virtual void PrepareTimeandStep();
    virtual void PrepareSolve();
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);
    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol);
    virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment);
    virtual void Update();
    virtual void StatisticsAndOutput();
    virtual void Output();
    virtual void StatisticsOutput();
    virtual IO::DiscretizationWriter& DiscWriter() { return fluid_.DiscWriter();}

    virtual void NonlinearSolve();
    virtual void Predictor();
    virtual void MultiCorrector();


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
                    Teuchos::RCP<const Epetra_Vector> readaccn);
    virtual double Time() const;
    virtual int Step() const;
    virtual double Dt() const;
    virtual void LiftDrag();
    virtual const double EvalTime() const { dserror("not implemented!"); return 0.0; }
    virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph) { dserror("not implemented!"); return; }
    virtual void TransferVectorsToNewDistribution(const Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle) { dserror("not implemented!"); return; }

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
    virtual Teuchos::RCP<Epetra_Vector>    RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();

    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

    virtual void SetTopOptPorosityField(Teuchos::RCP<Epetra_Vector> porosity);

    /// set initial flow field
    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno);

    ///  set fields for low-Mach-number flow within outer iteration loop
    virtual void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                                   RCP<const Epetra_Vector> scalaram,
                                   RCP<const Epetra_Vector> scalardtam,
                                   const double             thermpressaf,
                                   const double             thermpressam,
                                   const double             thermpressdtaf,
                                   const double             thermpressdtam,
                                   Teuchos::RCP<DRT::Discretization> scatradis);

    /// set fields for low-Mach-number flow at end of time step
    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                                   const double             thermpressnp,
                                   RCP<const Epetra_Vector> scatraresidual,
                                   Teuchos::RCP<DRT::Discretization> scatradis,
                                   const int                whichscalar = -1);

    /// set velocity field (separate computation)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp);

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
       {return fluid_.TurbulenceStatisticManager();}

  private:

    /// the actual fluid algorithm
    FLD::FluidImplicitTimeInt              fluid_;

    /// the interface map setup for fsi interface, free surface, interior translation
    FLD::UTILS::MapExtractor               interface_;
    Teuchos::RCP<Epetra_Vector>            interfaceforcen_;  /// interface force at old time level t_n

    /// ALE dof map
    LINALG::MapExtractor                   meshmap_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      dis_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RefCountPtr<Epetra_Map>       innervelmap_;
  };

}

#endif
#endif
