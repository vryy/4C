/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_xfem.H

\brief Fluid field adapter for xfem-fluids with moving boundaries

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_FLUID_XFEM_H
#define ADAPTER_FLUID_XFEM_H

#include "ad_fld_moving_boundary.H"
#include "ad_fld_base_algorithm.H"

namespace ADAPTER
{
  /// forward declarations

  class Coupling;

  /// fluid with moving interfaces implemented by the XFEM
  class FluidXFEM : public FluidMovingBoundary
  {
  public:

    /// constructor
    explicit FluidXFEM(const Teuchos::ParameterList& prbdyn,
                       std::string condname
                       );

    /*========================================================================*/
    //! @name Misc
    /*========================================================================*/

    /// fluid field
    const Fluid& FluidField() const { return fluid_.FluidField(); }
    Fluid& FluidField() { return fluid_.FluidField(); }

    /// return the boundary discretization that matches the structure discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization();

    /// communication object at the interface
    virtual const FLD::UTILS::MapExtractor& Interface() const;

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
    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,Teuchos::RCP<Epetra_Vector> ivel);

    /// relaxation solve
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt);
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
    virtual int Itemax() const { return FluidField().Itemax(); }

    /// set the maximum number of iterations for the fluid field
    virtual void SetItemax(int itemax) { FluidField().SetItemax(itemax); }

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
    FluidBaseAlgorithm fluid_;

    /// coupling of interface in to parallel distributions
    Teuchos::RCP<Coupling> icoupsf_;

  };

}

#endif
