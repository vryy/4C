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
class Epetra_MultiVector;
class Epetra_Vector;

namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidAdjoint3ImplParameter;
  }
}

namespace OPTI
{
  class GCMMA;
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

  The optimizer contains global as well as element routines. All functions
  containing element routines end up with a "3" and there code is split to a
  separate file.

  \author winklmaier
  \date 12/11
  */
  class Optimizer
  {
  public:
    /// constructor
    explicit Optimizer(
        Teuchos::RCP<DRT::Discretization> optidis,
        Teuchos::RCP<DRT::Discretization> fluiddis,
        const Teuchos::ParameterList& params
    );

    /// desctructor
    virtual ~Optimizer(){};

    /// compute the value of the objective funtion
    void ComputeValues();

    /// compute the gradient of the objective function
    void ComputeGradients();

    /// import time integration parameters
    void ImportFlowParams(Teuchos::RCP<Teuchos::ParameterList>& fluidParams);

    /// import the fluid solution of one time step
    void ImportFluidData(Teuchos::RCP<Epetra_Vector> vel,int step);

    /// import the adjoint fluid solution of one time step
    void ImportAdjointFluidData(Teuchos::RCP<Epetra_Vector> vel,int step);

    /// clear the field data for a new optimization step
    void ClearFieldData()
    {
      fluidvel_->clear();
      adjointvel_->clear();
    }

    void Iterate(
        bool& doGradient
    );

    void FinishIteration(
        bool& doGradient
    );

    bool Converged(
        bool& doGradient
    );

    Teuchos::RCP<const DRT::Discretization> OptiDis(){return optidis_;}
    Teuchos::RCP<const std::map<int,Teuchos::RCP<Epetra_Vector> > > ExportFluidData() const {return fluidvel_;} /// fluid velocity
    Teuchos::RCP<const std::map<int,Teuchos::RCP<Epetra_Vector> > > ExportAdjointData() const {return adjointvel_;} /// fluid velocity
    Teuchos::RCP<const Epetra_Vector> Density() const {return dens_;}
    const double ObjectiveValue(bool compute = false)
    {
      if (compute) ComputeValues();
      return obj_value_;
    }
    const Teuchos::RCP<const Epetra_Vector> ObjectiveGradient(bool compute)
    {
      if (compute) ComputeGradients();
      return obj_grad_;
    }
    const Epetra_Map* RowMap();
    const Epetra_Map* ColMap();

  protected:

    /// set the initial density field of the optimization variable
    void SetInitialDensityField(const INPAR::TOPOPT::InitialDensityField init,const int startfuncno);

    /// check if all required data for an optimization step is present
    bool DataComplete() const;

    void TransformFlowFields();

    //! @name general stuff
    Teuchos::RCP<DRT::Discretization> optidis_; /// discretization of optimization field
    Teuchos::RCP<DRT::Discretization> fluiddis_; /// discretization of optimization field
    const Teuchos::ParameterList& params_; /// parameters
    Teuchos::RCP<Teuchos::ParameterList> fluidParams_; /// fluid parameters

    //@}

    //! @name fluid fields

    /// velocities of all time steps
    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > fluidvel_;
    /// adjoint velocities of all time steps
    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > adjointvel_;

    //@}

    //! @name fields for optimization
    Teuchos::RCP<Epetra_Vector> dens_; /// new density in row map

    double obj_value_; /// value of the objective (adjoint solution not required)
    Teuchos::RCP<Epetra_Vector> obj_grad_; /// gradient of the objective (with help of the adjoint solution)

    double* constr_; /// values of the constraints
    Teuchos::RCP<Epetra_MultiVector> constr_deriv_; /// derivatives of the constraints w.r.t. the optimization variables
    int num_constr_; /// number of constraints

    Teuchos::RCP<OPTI::GCMMA> optimizer_; /// optimizer
    //@}

  }; // end class Optimizer
} // end namespace TOPOPT

#endif // #ifndef TOPOPT_OPTIMIZER_H_
