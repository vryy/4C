/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer.H

\brief optimizer of the topology optimization

\level 3

\maintainer Martin Winklmaier
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_OPTIMIZER_H_
#define TOPOPT_OPTIMIZER_H_


#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_topopt.H"
#include <Epetra_MultiVector.h>

// forward declarations
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseVector;
class Epetra_Vector;

namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;
  class ResultTest;

  namespace ELEMENTS
  {
    class FluidAdjoint3ImplParameter;
  }
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
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
    explicit Optimizer(Teuchos::RCP<DRT::Discretization> optidis,
        Teuchos::RCP<DRT::Discretization> fluiddis, const Teuchos::ParameterList& params,
        Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// desctructor
    virtual ~Optimizer(){};

    /// setting of objective function and constraints
    void Values() { ComputeValues(obj_, *constr_); };

    /// setting of gradients of objective function and constraints
    void Gradients() { ComputeGradients(obj_der_, constr_der_); };

    /// compute the gradient of the objective function in one direction
    void ComputeGradientDirectionForFD(
        const double fac, const int GID, const int index, const int numFDPoints);

    /// write output
    void Output();

    /// import time integration parameters
    void ImportFlowParams(Teuchos::RCP<Teuchos::ParameterList>& fluidParams);

    /// import the fluid solution of one time step
    void ImportFluidData(Teuchos::RCP<Epetra_Vector> vel, int step);

    /// import the adjoint fluid solution of one time step
    void ImportAdjointFluidData(Teuchos::RCP<Epetra_Vector> vel, int step);

    /// clear the field data for a new optimization step
    void ClearFieldData()
    {
      fluidvel_->clear();
      adjointvel_->clear();
    };

    /// perform one optimization iteration
    void Iterate(bool& doGradient);

    /// finish one optimization iteration
    void FinishIteration(bool& doGradient);

    /// check convergence of the optimization
    bool Converged(bool& doGradient);

    /// create result test object for optimization variables
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    Teuchos::RCP<const DRT::Discretization> OptiDis() { return optidis_; }
    const Teuchos::ParameterList& OptiParams() { return params_; }
    Teuchos::RCP<const std::map<int, Teuchos::RCP<Epetra_Vector>>> ExportFluidData() const
    {
      return fluidvel_;
    }  /// fluid velocity
    Teuchos::RCP<const std::map<int, Teuchos::RCP<Epetra_Vector>>> ExportAdjointData() const
    {
      return adjointvel_;
    }  /// fluid velocity
    Teuchos::RCP<Epetra_Vector> Density() const { return dens_; }

    double ObjectiveValue() { return obj_; };
    const Teuchos::RCP<const Epetra_Vector> ObjectiveGradient() { return obj_der_; };
    const Teuchos::RCP<const Epetra_SerialDenseVector> ConstraintValues() { return constr_; };
    const Teuchos::RCP<const Epetra_MultiVector> ConstraintGradients() { return constr_der_; };

    const Epetra_Map* RowMap();
    const Epetra_Map* ColMap();

    void UpdateOptimizationParameter();

    void AdoptDensityForFD(const double value, const int GID);

    int Iter() const;

    void ReadRestart(const int step);

   protected:
    /// computation of objective function and constraints
    void ComputeValues(double& objective, Epetra_SerialDenseVector& constraints);

    /// computation of gradients of objective function and constraints
    void ComputeGradients(
        Teuchos::RCP<Epetra_Vector> obj_der, Teuchos::RCP<Epetra_MultiVector> constr_der);

    /// set the initial density field of the optimization variable
    void SetInitialDensityField(
        const INPAR::TOPOPT::InitialDensityField init, const int startfuncno);

    /// check if all required data for an optimization step is present
    bool CheckData(const bool doAdjoint);

    void TransformFlowFields(const bool doAdjoint, const bool rowToCol);

    //! @name general stuff
    Teuchos::RCP<DRT::Discretization> optidis_;         /// discretization of optimization field
    Teuchos::RCP<DRT::Discretization> fluiddis_;        /// discretization of optimization field
    const Teuchos::ParameterList& params_;              /// parameters
    INPAR::TOPOPT::GradientType gradienttype_;          /// type of gradient computation
    Teuchos::RCP<Teuchos::ParameterList> fluidParams_;  /// fluid parameters

    //@}

    //! @name fluid fields

    /// velocities of all time steps
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>> fluidvel_;
    /// adjoint velocities of all time steps
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_Vector>>> adjointvel_;

    //@}

    //! @name fields for optimization
    Teuchos::RCP<Epetra_Vector> dens_;  /// new density in col map

    double obj_;  /// value of the objective (adjoint solution not required)
    Teuchos::RCP<Epetra_Vector>
        obj_der_;  /// gradient of the objective (with help of the adjoint solution)

    Teuchos::RCP<Epetra_SerialDenseVector> constr_;  /// values of the constraints
    Teuchos::RCP<Epetra_MultiVector>
        constr_der_;  /// derivatives of the constraints w.r.t. the optimization variables
    int num_constr_;  /// number of constraints

    Teuchos::RCP<Epetra_SerialDenseVector> objective_FD_;
    Teuchos::RCP<Epetra_SerialDenseMatrix> constraints_FD_;

    Teuchos::RCP<OPTI::GCMMA> optimizer_;  /// optimizer
    //@}

  };  // end class Optimizer
}  // end namespace TOPOPT

#endif  // #ifndef TOPOPT_OPTIMIZER_H_
