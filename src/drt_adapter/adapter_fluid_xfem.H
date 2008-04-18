/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_xfem.H

\brief 

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_XFEM_H
#define ADAPTER_FLUID_XFEM_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "adapter_fluid_moving_boundary.H"
#include "adapter_fluid_base_algorithm.H"
#include "../drt_fsi/fsi_coupling.H"

namespace ADAPTER
{

/// linear ALE algorithm
/*!
  The most basic, most simple implementation.
 */
class XFluidFSIBoundary : public FluidMovingBoundary
{
  friend class BoundaryResultTest;

public:

    XFluidFSIBoundary(Teuchos::RCP<DRT::Discretization> actdis,
            Teuchos::RCP<LINALG::Solver> solver,
            Teuchos::RCP<ParameterList> params,
            Teuchos::RCP<IO::DiscretizationWriter> output,
            int aletype,
            bool dirichletcond=true);

  /// build system matrix
  virtual void BuildSystemMatrix(bool full=true);

  /// a very simple time loop to be used for standalone ALE problems
  void Integrate();

  //! @name Algorithm core methods
  virtual void PrepareTimeStep();
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp) const;
  virtual void Solve();
  virtual void Update();
  virtual void Output();
  //@}

  /// get the linear solver object used for this field
  Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

  void SetInterfaceMap(Teuchos::RCP<Epetra_Map> im);

  void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp);
  void ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp);

  virtual Teuchos::RCP<Epetra_Vector> ExtractDisplacement() const;

  //! @name access functions for monolithic FSI

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess()
    { return residual_; }
  virtual Teuchos::RCP<Epetra_Vector> RHS() const
    { return residual_; }
  virtual Teuchos::RCP<Epetra_Vector> Disp()
    { return dispnp_; }
  virtual Teuchos::RCP<const Epetra_Map> DofRowMap() const
    { return Teuchos::rcp(discret_->DofRowMap(),false); }
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const
  { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }

  /// Aii matrix block
  virtual const LINALG::SparseMatrix* InteriorMatrixBlock() const;

  /// Aig matrix block
  virtual const LINALG::SparseMatrix* InterfaceMatrixBlock() const;

  virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return discret_; }
  //virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() const;

  //@}

  /// read restart for given step
  virtual double ReadRestart(int step);

  /// communication object at the interface
  const LINALG::MapExtractor& Interface() const { return interface_; }

  const LINALG::MapExtractor& FreeSurface() const { return freesurface_; }

private:

  /// the interface map setup for interface <-> full translation
  LINALG::MapExtractor interface_;

  LINALG::MapExtractor freesurface_;

  //! @name local copies of input parameters
  Teuchos::RCP<DRT::Discretization> discret_;
  Teuchos::RCP<LINALG::Solver> solver_;
  Teuchos::RCP<ParameterList> params_;
  Teuchos::RCP<IO::DiscretizationWriter> output_;
  //@}

  //! @name Algorithm core variables
  int step_;
  int numstep_;
  double time_;
  double maxtime_;
  double dt_;

  int aletype_;

  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  Teuchos::RCP<Epetra_Vector>    dirichtoggle_;
  //Teuchos::RCP<Epetra_Vector>    zeros_;

  Teuchos::RCP<Epetra_Vector>    dispnp_;
  Teuchos::RCP<Epetra_Vector>    dispn_;
  //Teuchos::RCP<Epetra_Vector>    dispnm_;

  Teuchos::RCP<Epetra_Vector>    residual_;
  //@}

private:

  /// element call
  void EvaluateElements();

  //! @name restart
  int restartstep_;
  int uprestart_;
  //@}

};

  /// fluid with moving interfaces implemented by the XFEM
  class FluidXFEM : public FluidMovingBoundary
  {
  public:

    FluidXFEM(const Teuchos::ParameterList& prbdyn,
                    std::string condname);

    /// fluid field
    Fluid& FluidField() { return boundary_.FluidField(); }
    const Fluid& FluidField() const { return boundary_.FluidField(); }

    /// ale field
//    Boundary& BoundaryField() { return boundary_.BoundaryField(); }
//    const Boundary& BoundaryField() const { return boundary_.BoundaryField(); }

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual const LINALG::MapExtractor& Interface() const;

    virtual void PrepareTimeStep();
    virtual void Update();
    virtual void Output();
    virtual double ReadRestart(int step);

    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,Teuchos::RCP<Epetra_Vector> ivel);
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    virtual int Itemax() const { return FluidField().Itemax(); }
    virtual void SetItemax(int itemax) { FluidField().SetItemax(itemax); }

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  protected:

    //! @name Transfer helpers

//    /// field transform
//    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;
//
//    /// field transform
//    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const;
//
//    /// interface transform
//    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;
//
//    /// interface transform
//    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    //Coupling& FluidAleFieldCoupling()             { return coupfa_; }
    //const Coupling& FluidAleFieldCoupling() const { return coupfa_; }

  private:

    FluidBaseAlgorithm fluid_;
    FluidBaseAlgorithm boundary_;

    /// coupling of fluid and boundary (interface only)
    FSI::Coupling               icoupfb_;

  };

}

#endif
#endif
