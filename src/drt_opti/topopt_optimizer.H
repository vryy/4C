/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer.H

\brief optimizer of the topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_OPTIMIZER_H_
#define TOPOPT_OPTIMIZER_H_


#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_topopt.H"


/// forward declarations
class Epetra_Map;
class Epetra_Vector;

namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;
}


namespace TOPOPT
{
  /*
  \brief The optimizer for topology optimization

  This class shall perform the optimization step(s) due to an appropriate
  optimization approach. It shall therefore use the primary solution (e.g.
  Navier-Stokes equation) as well as the adjoint solution. With these data
  the optimization gradient has to be built, followed by the optimization step.

  This class shall hold all the data an optimizer requires:
  - primal solutions (stationary or instationary)
  - dual solutions (the same)
  - the ingredients of the objective function
  - the optimization variables with its values
  It shall not have data which it does not really need. Such data shall be hold
  by the main function of the optimization, e.g. topopt_algorithm.H

  If the examples become larger, it might be sensible to write the solutions
  into files and read these files here. Especially in the case of instationary
  simulations.

  Possibly optimization algorithms are the MMA (method of moving asymptotes),
  the GCMMA (globally convergent MMA), see Svanberg or the Trust-Region method.

  \author winklmaier
  \date 12/11
  */
  class Optimizer
  {
  public:
    /// constructor
    explicit Optimizer(Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& params);

    /// desctructor
    virtual ~Optimizer(){};

    /// compute the value of the objective funtion
    double ComputeObjectiveValue(Teuchos::RCP<Epetra_Vector> porosity);

    /// import the fluid solution of one time step
    void ImportFluidData(Teuchos::RCP<Epetra_Vector> vel,int step);

    /// import the adjoint fluid solution of one time step
    void ImportAdjointFluidData(Teuchos::RCP<Epetra_Vector> vel, int step);

    /// export the complete fluid solution of all time steps
    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > ExportFluidData(){return vel_;} /// fluid velocity

    /// density field of last optimization step
    Teuchos::RCP<Epetra_Vector> DensityI(){return dens_i_;}

    /// density field of new optimization step
    Teuchos::RCP<Epetra_Vector> DensityIp(){return dens_ip_;}

    /// discretization of fluid field
    Teuchos::RCP<DRT::Discretization> Discret(){return discret_;}

    /// row map of optimization variable
    const Epetra_Map* RowMap();

    /// column map of optimization variable
    const Epetra_Map* ColMap();

  protected:

    /// set the initial density field of the optimization variable
    void SetInitialDensityField(const INPAR::TOPOPT::InitialDensityField init,const int startfuncno);

    /// evaluate the objective function
    double EvaluateObjective(
        Teuchos::RCP<Epetra_Vector> porosity,
        const std::string& action);

  private:

    /// check if all required data for an optimization step is present
    bool DataComplete() const;

    Teuchos::RCP<Epetra_Vector> dens_i_;/// old density in row map
    Teuchos::RCP<Epetra_Vector> dens_ip_; /// new density in row map
    Teuchos::RCP<DRT::Discretization> discret_; /// discretization of optimization field
    const Teuchos::ParameterList& params_; /// parameters

    //! @name fields for optimization

    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > vel_; /// velocities of all time steps
    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > adjointvel_; /// adjoint velocities of all time steps

    size_t num_timesteps_; /// number of time steps
    //@}

    //! @name parameter of the objective function

    bool dissipation_; //! optimize dissipation
    bool pressure_inlet_; //! optimize pressure at inlet
    bool pressure_drop_; //! optimize pressure drop from inlet to outlet

    double dissipation_fac_; //! dissipation weighting factor
    double pressure_inlet_fac_; //! inlet pressure weighting factor
    double pressure_drop_fac_; //! inlet pressure weighting factor
    //@}

  }; // end class Optimizer
} // end namespace TOPOPT

#endif // #ifndef TOPOPT_OPTIMIZER_H_
