/*----------------------------------------------------------------------*/
/*!
\file mfsi_fluid.H

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

#ifndef MFSI_FLUID_H
#define MFSI_FLUID_H

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_fluid/fluidimplicitintegration.H"

#include "../drt_fsi/fsi_conditiondofmap.H"


namespace MFSI {

  /// general fluid field interface for monolithic FSI
  /*!

  The point is to keep monolithic FSI as far apart from our field
  solvers as possible. Each fluid field solver we want to use should
  get its own subclass of this. The monolithic FSI algorithm should be
  able to extract all the information from the fluid field it needs
  using this interface.

  \warning A cleanup is needed badly.

  \author u.kue
  \date 11/07
  */
  class Fluid
  {
  public:

    virtual ~Fluid();

    virtual Teuchos::RCP<Epetra_Vector> InitialGuess() = 0;
    virtual Teuchos::RCP<Epetra_Vector> RHS() = 0;
    virtual Teuchos::RCP<Epetra_Vector> Vel() = 0;
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() = 0;
    virtual Teuchos::RCP<Epetra_Operator> SysMat() = 0;
    virtual Teuchos::RCP<DRT::Discretization> Discretization() = 0;

    virtual Teuchos::RCP<Epetra_Vector> StructCondRHS() = 0;

    virtual void PrepareTimeStep() = 0;
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel) = 0;
    virtual void Update() = 0;
    virtual void Output() = 0;

    /// the interface map contains all velocity dofs at the interface
    virtual void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im) = 0;
    virtual Teuchos::RCP<Epetra_Map> InterfaceMap() = 0;
    virtual Teuchos::RCP<Epetra_Map> InnerVelocityRowMap() = 0;
    virtual Teuchos::RCP<Epetra_Map> PressureRowMap() = 0;

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<Epetra_Map> mm) = 0;

    /// tell the mesh displacement to the fluid solver
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) = 0;

    virtual double ResidualScaling() = 0;

    /// communication object at the interface
    virtual const FSI::ConditionDofMap& Interface() const = 0;
  };


  /// adapter to normal fluid algorithm
  class FluidAdapter : public Fluid
  {
  public:
    FluidAdapter(Teuchos::RCP<DRT::Discretization> dis,
                 Teuchos::RCP<LINALG::Solver> solver,
                 Teuchos::RCP<ParameterList> params,
                 Teuchos::RCP<IO::DiscretizationWriter> output);

    virtual Teuchos::RCP<Epetra_Vector> InitialGuess();
    virtual Teuchos::RCP<Epetra_Vector> RHS();
    virtual Teuchos::RCP<Epetra_Vector> Vel();
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();
    virtual Teuchos::RCP<Epetra_Operator> SysMat();
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    virtual Teuchos::RCP<Epetra_Vector> StructCondRHS();

    virtual void PrepareTimeStep();
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);
    virtual void Update();
    virtual void Output();

    /// the interface map contains all velocity dofs at the interface
    virtual void SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im);
    virtual Teuchos::RCP<Epetra_Map> InterfaceMap();
    virtual Teuchos::RCP<Epetra_Map> InnerVelocityRowMap();
    virtual Teuchos::RCP<Epetra_Map> PressureRowMap();

    /// the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<Epetra_Map> mm);

    /// tell the mesh displacement to the fluid solver
    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp);

    virtual double ResidualScaling();

    /// communication object at the interface
    virtual const FSI::ConditionDofMap& Interface() const { return interface_; }

  private:

    /// the actual fluid algorithm
    FluidImplicitTimeInt fluid_;

    /// the interface map setup for interface <-> full translation
    FSI::ConditionDofMap interface_;

    /// ALE dof map
    FSI::ConditionDofMap meshmap_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization> dis_;
    Teuchos::RCP<LINALG::Solver> solver_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RefCountPtr<Epetra_Map> innervelmap_;

    /// what a hack
    Teuchos::RCP<Epetra_Vector> sumincvel_;
  };

}

#endif
#endif
