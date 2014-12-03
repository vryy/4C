/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.H

\brief Fluid field adapter for fsi

Can only be used in conjunction with #FluidImplicitTimeInt

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_FSI_H
#define ADAPTER_FLUID_FSI_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "ad_fld_wrapper.H"

namespace DRT
{
  class Condition;
}

namespace LINALG
{
  class Solver;
  class MapExtractor;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace FLD
{
  class FluidImplicitTimeInt;
  namespace UTILS
  {
    class MapExtractor;
  }
}

namespace ADAPTER {

  /*! \brief Fluid field adapter for fsi
   *
   *
   *  Can only be used in conjunction with #FLD::FluidImplicitTimeInt
   */
  class FluidFSI : public FluidWrapper
  {
  public:

    /// Constructor
    FluidFSI(Teuchos::RCP<Fluid> fluid,
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output,
        bool isale,
        bool dirichletcond);

    /// initialize algorithm
    virtual void Init();

    virtual Teuchos::RCP<const Epetra_Map> DofRowMap();

    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds);

    /// Velocity-displacement conversion at the fsi interface
    virtual double TimeScaling() const;

    /// take current results for converged and save for next time step
    virtual void Update();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver();

    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel);

    /// communication object at the interface
    virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }

    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    /// Return interface velocity at new time level n+1
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp();

    /// Return interface velocity at old time level n
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();

    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln();

    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");};

    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel);

    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm);

    //! @name Conversion between displacement and velocity at interface

    //! Conversion of displacement to velocity at the interface without predictors or inhomogeneous DBCs
    //!
    //! All input vectors have to live on the fluid field map.
    virtual void DisplacementToVelocity(
        Teuchos::RCP<Epetra_Vector> fcx     ///< interface displacement step increment
    );

    //! Conversion of velocity to displacement at the interface without predictors or inhomogeneous DBCs
    //!
    //! All input vectors have to live on the fluid field map.
    virtual void VelocityToDisplacement(
        Teuchos::RCP<Epetra_Vector> fcx     ///< interface velocity step increment at interface
    );

    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx);

    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx);

    //@}

    virtual Teuchos::RCP<Epetra_Vector>   IntegrateInterfaceShape();

    virtual void UseBlockMatrix(bool splitmatrix);

    //! Project the velocity field into a divergence free subspace
    //!
    //! Project the velocity field into a divergence free subspace
    //! while interface and Dirichlet DOFS are not affected.
    //! The projection is done by the following operation:
    //!
    //! \$f v_{divfree} = \left(I - B(B^TB)^{-1}B^T)\right) v + B(B^TB)^{-1} R\$f
    //!
    //! The vector \$f R \$f ensures that interface and Dirichlet DOFs are not modified.
    //!
    //! \autor: mayr.mt
    //! \date  06/2012
    void ProjVelToDivZero();

    /// reset state vectors
    virtual void Reset(
        bool completeReset = false,
        bool newFiles = false,
        int iter = -1);

    /// calculate error in comparison to analytical solution
    virtual void CalculateError();

    //! @name Time step size adaptivity in monolithic FSI
    //@{

    /*! Do one step with auxiliary time integration scheme
     *
     *  Do a single time step with the user given auxiliary time integration
     *  scheme. Result is stored in #locerrvelnp_ and is used later to estimate
     *  the local discretization error of the marching time integration scheme.
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual void TimeStepAuxiliar();

    /*! Indicate norms of temporal discretization error
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual void IndicateErrorNorms(
        double& err,        ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,    ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,   ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,     ///< L-inf-norm of temporal discretization error based on all DOFs
        double& errinfcond, ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother ///< L-inf-norm of temporal discretization error based on interior DOFs
        );

    /// provide access to fluid time integrator
    Teuchos::RCP<FLD::FluidImplicitTimeInt> FluidImplTimeInt()
    {return fluidimpl_;};

    ///Calculate WSS vector
    virtual Teuchos::RCP<Epetra_Vector> CalculateWallShearStresses();

    //@}

  protected:

    /// create conditioned dof-map extractor for the fluid
    void SetupInterface();

    /// build inner velocity dof-map
    void BuildInnerVelMap();

    /// A casted pointer to the fluid itself
    Teuchos::RCP<FLD::FluidImplicitTimeInt> fluidimpl_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      dis_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    bool dirichletcond_;
    //@}

    //! \brief interface map setup for fsi interface, free surface, interior translation
    //!
    //! Note: full map contains velocity AND pressure DOFs
    Teuchos::RCP<FLD::UTILS::MapExtractor> interface_;

    /// interface force at old time level t_n
    Teuchos::RCP<Epetra_Vector>            interfaceforcen_;

    /// ALE dof map
    Teuchos::RCP<LINALG::MapExtractor>     meshmap_;

    /// all velocity dofs not at the interface
    Teuchos::RCP<Epetra_Map>               innervelmap_;

  private:

    //! Time step size adaptivity in monolithic FSI
    //@{

    /*! \brief Do an explicit Euler step as auxiliary time integrator
     *
     *  Based on state vector \f$x_n\f$ and its time derivative \f$\dot{x}_{n}\f$
     *  at time \f$t_{n}\f$, we calculate \f$x_{n+1}\f$ using an explicit Euler
     *  step:
     *
     *  \f[
     *    x_{n+1} = x_{n} + \Delta t_{n} \dot{x}_{n}
     *  \f]
     *
     *  \author mayr.mt
     *  \date 10/2013
     */
    void ExplicitEuler(const Epetra_Vector& veln, ///< velocity at \f$t_n\f$
                       const Epetra_Vector& accn, ///< acceleration at \f$t_n\f$
                       Epetra_Vector& velnp       ///< velocity at \f$t_{n+1}\f$
                       ) const;

    /*! \brief Do an Adams-Bashforth 2 step as auxiliary time integrator
     *
     *  Based on state vector \f$x_n\f$ and its time derivatives \f$\dot{x}_{n}\f$
     *  and \f$\dot{x}_{n-1}\f$ at time steps \f$t_{n}\f$ and \f$t_{n-1}\f$,
     *  respectively, we calculate \f$x_{n+1}\f$ using an Adams-Bashforth 2 step
     *  with varying time step sizes:
     *
     *  \f[
     *    x_{n+1} = x_{n}
     *            + \frac{2\Delta t_{n} \Delta t_{n-1} + \Delta t_{n}^2}{2\Delta t_{n-1}} \dot{x}_{n}
     *            - \frac{\Delta t_{n}^2}{2\Delta t_{n-1}} \dot{x}_{n-1}
     *  \f]
     *
     *  <h3> References: </h3>
     *
     *  - Wall, W. A.: Fluid-Struktur-Interaktion mit stabilisierten Finiten
     *  Elementen, PhD Thesis, Universitaet Stuttgart, 1999
     *
     *  - Bornemann, B.: Time Integration Algorithms for the Steady States of
     *  Dissipative Non-Linear Dynamic Systems, PhD Thesis, Imperial College
     *  London, 2003
     *
     *  - Gresho, P. M., Griffiths, D. F., Silvester, D. J.: Adaptive Time-Stepping
     *  for Incompressible Flow Part I: Scalar Advection-Diffusion,
     *  SIAM J. Sci. Comput. (30), pp. 2018-2054, 2008
     *
     *  - Kay, D. A., Gresho, P. M., Griffiths, D. F., Silvester, D. J.: Adaptive
     *  Time Stepping for Incompressible Flow Part II: Navier-Stokes Equations,
     *  SIAM J. Sci. Comput. (32), pp. 111-128, 2010
     *
     *  \author mayr.mt
     *  \date 11/2013
     */
    void AdamsBashforth2(const Epetra_Vector& veln,   ///< velocity at \f$t_n\f$
                         const Epetra_Vector& accn,   ///< acceleration at \f$t_n\f$
                         const Epetra_Vector& accnm,  ///< acceleration at \f$t_{n-1}\f$
                         Epetra_Vector& velnp         ///< velocity at \f$t_{n+1}\f$
                         ) const;

    //! Compute length-scaled L2-norm of a vector
    virtual double CalculateErrorNorm(const Epetra_Vector& vec, ///< vector to compute norm of
                                      const int numneglect = 0  ///< number of DOFs that have to be neglected for length scaling
                                       ) const;

    Teuchos::RCP<Epetra_Vector> locerrvelnp_; ///< vector of temporal local discretization error

    int auxintegrator_; ///< auxiliary time integrator in fluid field

    int numfsidbcdofs_; ///< number of interface DOFs with Dirichlet boundary condition

    //@}

  };
}

#endif
