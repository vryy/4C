/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for immersed-fluids

\level 3

*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_IMMERSED_H
#define ADAPTER_FLUID_IMMERSED_H

#include "ad_fld_moving_boundary.H"
#include "ad_fld_base_algorithm.H"

namespace ADAPTER
{
  // forward declarations
  class Coupling;

  /// fluid with moving interfaces
  class FluidImmersed : public FluidMovingBoundary
  {
   public:
    /// constructor
    explicit FluidImmersed(const Teuchos::ParameterList& prbdyn, std::string condname);

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    virtual const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_->FluidField(); }

    /// return the boundary discretization that matches the structure discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const;

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
    /// Not used for IMMERSED yet !

    /// get the maximum number of iterations from the fluid field
    virtual int Itemax() const { return fluidadapter_->Itemax(); }

    /// set the maximum number of iterations for the fluid field
    virtual void SetItemax(int itemax) { FluidField()->SetItemax(itemax); }

    /// add dirichlet conditions to dirichlet condmap before next fluid solve
    virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    /// remove dirichlet conditions from dirichlet condmap after last fluid solve
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

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
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid_;
    Teuchos::RCP<ADAPTER::Fluid> fluidadapter_;
  };

}  // namespace ADAPTER

#endif
