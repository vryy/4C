/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_xfluid_impl.H


<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_XFLUID_IMPL_H_
#define ADAPTER_FLUID_XFLUID_IMPL_H_

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/fluidxfluidimplicitintegration.H"
#include "../drt_fluid/fluidxfluidresulttest.H"
#include "../drt_lib/drt_utils.H"

#include "adapter_fluid.H"

namespace ADAPTER {
  class FluidXFluidImpl : public Fluid
  {
  public:
   FluidXFluidImpl(Teuchos::RCP<DRT::Discretization> fluiddis,
             Teuchos::RCP<DRT::Discretization> xfluiddis,
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
    virtual Teuchos::RCP<const Epetra_Vector> Velnm() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accn() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accam() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<const Epetra_Vector> Hist() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp();
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel();
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<LINALG::SparseMatrix>SystemMatrix() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() {dserror("not implemented"); return null;}
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() {dserror("not implemented!"); return null;}
    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual Teuchos::RCP<const DRT::DofSet> DofSet() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor();

    virtual Teuchos::RCP<Epetra_Vector> ExtractEmbVector(Teuchos::RCP<const Epetra_Vector> f) {dserror("not implemented"); return Teuchos::null;};

    virtual void TimeLoop();
    virtual void PrepareTimeStep();
    virtual void PrepareTimeandStep() { dserror("not implemented!"); return; }
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel) {dserror("not implemented");}
    virtual void Update();
    virtual void StatisticsAndOutput();
    virtual void Output();
    virtual void StatisticsOutput() { dserror("not implemented!"); return; }
    virtual IO::DiscretizationWriter& DiscWriter() { return fluid_.DiscWriter();}

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

    virtual int  Itemax() const;
    virtual void SetItemax(int itemax);

    /*!
    \brief return type of time integration scheme

    */
    const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return fluid_.TimIntScheme();}

    virtual Teuchos::RCP<Epetra_Vector>   IntegrateInterfaceShape();

    virtual void UseBlockMatrix(bool splitmatrix) {dserror("not implemented");}

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    virtual Teuchos::RCP<Epetra_Vector>    RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) {dserror("not implemented"); return null;}

    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();

    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

    /// set initial flow field
    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno);

    ///  set fields for low-Mach-number flow within outer iteration loop
    virtual void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                                   RCP<const Epetra_Vector> scalaram,
                                   RCP<const Epetra_Vector> scalardtam,
                                   const double             thermpressaf,
                                   const double             thermpressam,
                                   const double             thermpressdtam,
                                   Teuchos::RCP<DRT::Discretization> scatradis) {dserror("not implemented");}

    /// set fields for low-Mach-number flow at end of time step
    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                                   const double             thermpressnp,
                                   RCP<const Epetra_Vector> scatraresidual,
                                   Teuchos::RCP<DRT::Discretization> scatradis,
                                   const int                whichscalar = -1)
                                   {dserror("not implemented");};

    /// set velocity field (separate computation)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp)
    { dserror("not implemented!");return; };

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
       {dserror("not implemented"); return Teuchos::null;}



  private:

    void PrepareFluidXFluidBoundaryDis();

    /// the actual fluid algorithm
    FLD::FluidXFluidImplicitTimeInt        fluid_;

    /// the interface map setup for fsi interface, free surface, interior translation
    FLD::UTILS::MapExtractor               interface_;

    /// ALE dof map
    LINALG::MapExtractor                   meshmap_;

    //! @name local copies of input parameters
    // fluid discretization
    Teuchos::RCP<DRT::Discretization>      fluiddis_;
    // xfluid discretization
    Teuchos::RCP<DRT::Discretization>      xfluiddis_;
    // fluid-xfluid-boundary discretization
    Teuchos::RCP<DRT::Discretization>      fluidxfluidboundarydis_;

    Teuchos::RCP<IO::DiscretizationWriter> fluidxfluidboundaryoutput_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    //! @name state vectors using the fluid-xfluid boundary parallel distribution (row vectors)
    Teuchos::RCP<Epetra_Vector>            fxfidispnp_;   ///< interface displacement  n+1
    Teuchos::RCP<Epetra_Vector>            fxfivelnp_;    ///< interface velocity      n+1
    Teuchos::RCP<Epetra_Vector>            fxfitrueresnp_;///< interface residual      n+1
    Teuchos::RCP<Epetra_Vector>            fxfidispn_;    ///< interface displacement  n
    Teuchos::RCP<Epetra_Vector>            fxfiveln_;     ///< interface velocity      n
    Teuchos::RCP<Epetra_Vector>            fxfivelnm_;    ///< interface velocity      n-1
    Teuchos::RCP<Epetra_Vector>            fxfiaccnp_;    ///< interface acceleration  n+1
    Teuchos::RCP<Epetra_Vector>            fxfiaccn_;     ///< interface acceleration  n
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RefCountPtr<Epetra_Map>       innervelmap_;
  };
}

#endif /* ADAPTER_FLUID_XFLUID_IMPL_H_ */
#endif
