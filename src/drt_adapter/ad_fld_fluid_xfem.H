/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_xfem.H

\brief Fluid field adapter for xfem-fluids with moving boundaries

\level 1

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
 </pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_XFEM_H
#define ADAPTER_FLUID_XFEM_H

#include "ad_fld_moving_boundary.H"
#include "ad_fld_base_algorithm.H"

namespace ADAPTER
{
  // forward declarations
  class Coupling;

  /// fluid with moving interfaces implemented by the XFEM
  class FluidXFEM : public FluidMovingBoundary
  {
   public:
    /// constructor
    explicit FluidXFEM(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_; }

    /// return the boundary discretization that matches the structure discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    /// return the boundary discretization that matches the structure discretization
    Teuchos::RCP<DRT::Discretization> BoundaryDiscretization();

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const
    {
      return fluid_->Interface();
    }

    /// communication object at the struct interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& StructInterface();

    //@}

    /*========================================================================*/
    //! @name Time step helpers
    /*========================================================================*/

    /// start new time step
    virtual void PrepareTimeStep();

    /// update at time step end
    virtual void Update();

    /// output results
    virtual void Output();

    /// read restart information for given time step
    virtual double ReadRestart(int step);

    /*========================================================================*/
    //! @name Solver calls
    /*========================================================================*/

    /// nonlinear solve
    virtual void NonlinearSolve(
        Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel);

    /// relaxation solve
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(
        Teuchos::RCP<Epetra_Vector> idisp, double dt);
    //@}

    /*========================================================================*/
    //! @name Extract interface forces
    /*========================================================================*/

    /// After the fluid solve we need the forces at the FSI interface.
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    //@}

    /*========================================================================*/
    //! @name extract helpers
    /*========================================================================*/

    /// extract the interface velocity at time t^(n+1)
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    /// extract the interface velocity at time t^n
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    //@}

    /*========================================================================*/
    //! @name Number of Newton iterations
    /*========================================================================*/

    //! For simplified FD MFNK solve we want to temporally limit the
    /// number of Newton steps inside the fluid solver

    /// get the maximum number of iterations from the fluid field
    virtual int Itemax() const { return fluid_->Itemax(); }

    /// set the maximum number of iterations for the fluid field
    virtual void SetItemax(int itemax) { fluid_->SetItemax(itemax); }

    //@}

    /*========================================================================*/
    //! @name others
    /*========================================================================*/

    /// integrate the interface shape functions
    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    /// create the testing of fields
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();



   private:
    /// fluid base algorithm object
    Teuchos::RCP<Fluid> fluid_;
  };

}  // namespace ADAPTER

#endif
