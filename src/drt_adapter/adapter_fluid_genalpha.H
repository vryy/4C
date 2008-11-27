/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_genalpha.H

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

#ifndef ADAPTER_FLUID_GENALPHA_H
#define ADAPTER_FLUID_GENALPHA_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluidresulttest.H"
#include "../drt_lib/drt_utils.H"

#include "adapter_fluid.H"


/// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// adapter to generalized alpha fluid algorithm
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
    \rho \f$ to be physical forces. The matrix columns
    that go with \f$ \mathbf{u}^{n+1}_{\Gamma,i+1} \f$ have to be scaled with
    \f$ \frac{1}{\Delta t} \f$. This results in the scaled blocks:

    - \f$ \frac{\rho}{\Delta t} \mathbf{F}_{\Gamma\Gamma} \f$

    - \f$ \frac{1}{\Delta t} \mathbf{F}_{I\Gamma} \f$

    - \f$ \rho \mathbf{F}_{\Gamma I} \f$

    \author pg
    \date 02/08
   */
  class FluidGenAlpha : public Fluid
  {
  public:
    FluidGenAlpha(
      Teuchos::RCP<DRT::Discretization>      dis,
      Teuchos::RCP<LINALG::Solver>           solver,
      Teuchos::RCP<ParameterList>            params,
      Teuchos::RCP<IO::DiscretizationWriter> output,
      bool isale
      );

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess();
    virtual Teuchos::RCP<const Epetra_Vector> RHS();
    virtual Teuchos::RCP<const Epetra_Vector> Velnp();
    virtual Teuchos::RCP<const Epetra_Vector> Velaf();
    virtual Teuchos::RCP<const Epetra_Vector> Veln();
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp();
    virtual Teuchos::RCP<const Epetra_Map>    DofRowMap();
    virtual Teuchos::RCP<LINALG::SparseMatrix>SystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> MeshMoveMatrix();
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const;

    virtual void PrepareTimeStep();
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> dacc);
    virtual void Update();
    virtual void StatisticsAndOutput();
    virtual void Output();

    virtual void NonlinearSolve();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// the interface map contains all velocity dofs at the interface
    virtual Teuchos::RCP<const Epetra_Map>    InnerVelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    VelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map>    PressureRowMap();

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    /// scaling factor needed to convert the residual to real forces
    virtual double ResidualScaling() const;

    virtual double TimeScaling() const;

    /// communication object at the interface
    virtual const LINALG::MapExtractor&  Interface() const { return interface_; }

    virtual const LINALG::MapExtractor& FreeSurface() const { return freesurface_; }

    virtual void ReadRestart(int step);
    virtual double Time();
    virtual int Step();
    virtual void LiftDrag();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce);
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp);
    virtual void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel);

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    virtual int  Itemax() const;
    virtual void SetItemax(int itemax);

    /*!
    \brief return type of time integration scheme

    */
    virtual const FLUID_TIMEINTTYPE TimIntScheme() const {return fluid_.TimIntScheme();}

    virtual Teuchos::RCP<Epetra_Vector>   IntegrateInterfaceShape();

    virtual void UseBlockMatrix(const LINALG::MultiMapExtractor& domainmaps,
                                const LINALG::MultiMapExtractor& rangemaps,
                                bool splitmatrix);

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    virtual Teuchos::RCP<Epetra_Vector>    RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();

    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

    /// set initial flow field
    virtual void SetInitialFlowField(int whichinitialfield,int startfuncno);

    ///  set time-step-related fields for low-Mach-number flow
    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> densnp,
                                   RCP<const Epetra_Vector> densn,
                                   RCP<const Epetra_Vector> densdtn){ dserror("not implemented!"); return; };

    /// set outer-iteration-related fields for low-Mach-number flow
    virtual void SetIterLomaFields(RCP<const Epetra_Vector> densnp,
                                   RCP<const Epetra_Vector> densdtnp){ dserror("not implemented!"); return; };

  private:

    /// the actual fluid algorithm
    FLD::FluidGenAlphaIntegration          fluid_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor                   interface_;

    LINALG::MapExtractor                   freesurface_;

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
