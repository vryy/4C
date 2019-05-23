/*!----------------------------------------------------------------------
\brief Control routine for acoustic explicit time integration.

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*----------------------------------------------------------------------*/

#ifndef ACOU_EXPL_H
#define ACOU_EXPL_H

#include "acou_timeint.H"
#include "../drt_lib/drt_discret_hdg.H"

#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/base/index_set.h>

#include "acou_expl_worker.H"
#include "deal_tria_wrapper.H"

namespace ACOU
{
  class AcouExplicitTimeInt;

  using namespace dealii;

  template <typename Operator>
  class ExplicitIntegrator;

  template <int dim, typename Number>
  class WaveEquationProblem
  {
   public:
    typedef typename WaveEquationOperationBase<dim, Number>::value_type value_type;
    WaveEquationProblem(Teuchos::RCP<DRT::DiscretizationHDG> discretin,
        const Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<ACOU::AcouExplicitTimeInt> timeintin);

    ~WaveEquationProblem();

    // evaluator is set
    void set_evaluator(const Teuchos::RCP<Teuchos::ParameterList> params, const int fe_degree);

    // Integration, i.e. time loop
    void run();

    // preparation for output in baci format
    void write_deal_cell_values();

    // evaluator
    Teuchos::RCP<WaveEquationOperationBase<dim, Number>> evaluator;

    // actual and previous pressure and velocity solution
    std::vector<parallel::distributed::Vector<value_type>> solutions, previous_solutions;

    // vector of up_res solution vectors of pressure and velocity (for adjoint with acouopt)
    std::vector<std::vector<parallel::distributed::Vector<value_type>>> stored_forward_solutions;

    // time step (public for adaption in base class)
    double time_step;

    // set time and step, required for restart
    void set_time_and_step(double timein, int stepin);

    // get sos gradient value for row element with id rowid
    double get_SoS_gradient(int rowid);

    // get density gradient value for row element with id rowid
    double get_density_gradient(int rowid);

   private:
    void make_grid_and_dofs(Teuchos::RCP<DRT::DiscretizationHDG> discret);
    void output_results(const unsigned int timestep_number);
    void write_restart(const unsigned int timestep_number);
    void read_restart(const unsigned int timestep_number, double& restarttime, bool fw);
    void compute_post_pressure();
    void compute_post_velocity();
    Number determine_pressure_value(Point<dim> point, unsigned int cell);

    std::vector<parallel::distributed::Vector<value_type>>
        post_quantity;  // pressure for fluid, velocity for solid

    ConditionalOStream pcout;

    DistributedTriangulation<dim> triangulation;
    FE_DGQArbitraryNodes<dim> fe, fe_spectral, fe_post_disp;
    DoFHandler<dim> dof_handler, dof_handler_spectral, dof_handler_post_disp;

    IndexSet locally_relevant_dofs, loc_disp;

    // actual time, time step size, final time
    double time, final_time;

    // maximal step, write restart, write output
    int step_max, up_res, step;

    // cfl number (determined during construction)
    const double cfl_number;

    // time integration scheme
    INPAR::ACOU::DynamicType dyna;

    // exact solution given by number
    int exactsolutionfuncno;

    // pointer to the baci discretization
    Teuchos::RCP<DRT::DiscretizationHDG> discret;

    // pointer to the baci time integration (for baci output)
    Teuchos::RCP<ACOU::AcouExplicitTimeInt> bacitimeint;

    // flags needed for image reconstruction
    bool invana, adjoint, acouopt, timereversal, solid, reduction, writebacioutput;

    // forward time integration in the adjoint run
    void intermediate_integrate(
        Teuchos::RCP<ExplicitIntegrator<WaveEquationOperationBase<dim, Number>>> integrator);
  };

  class AcouExplicitTimeInt : public AcouTimeInt
  {
   public:
    /// constructor
    AcouExplicitTimeInt(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// virtual destructor
    virtual ~AcouExplicitTimeInt();

    /// read restart information for given time step
    virtual void ReadRestart(int step);

    /// set initial field to zero (needed for adjoint problem)
    void SetInitialZeroField();

    /// set initial field by given function
    void SetInitialField(int startfuncno);

    /// set initial field by given light distribution
    // void SetInitialPhotoAcousticField(Teuchos::RCP<SCATRA::TimIntStationary> scatraalgo);
    void SetInitialPhotoAcousticField(Teuchos::RCP<Epetra_Vector> light,
        Teuchos::RCP<DRT::Discretization> scatradis, bool meshconform);
    /// run a complete simulation
    void Integrate();

    /// Output nodal internal pressure field at new time step
    virtual void NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec);

    /// return the name of the time integration scheme
    std::string Name();

    /// return the protected parameter list
    const Teuchos::RCP<Teuchos::ParameterList> Params() { return params_; }

    /// return the protected processor id
    int ProcId() { return myrank_; }

    // get values for speed of sound gradient
    double GetSoSGradient(int rowid);

    // get values for density gradient
    double GetDensityGradient(int rowid);

    virtual void NodalPsiField(Teuchos::RCP<Epetra_Vector> outvec);

   private:
    /// in case the time step size has to be smaller, we have to tell the base class
    void UpdateTimeStepSize();

    // the two possible wave equation problems
    Teuchos::RCP<WaveEquationProblem<2, double>> wave2dd_;
    Teuchos::RCP<WaveEquationProblem<3, double>> wave3dd_;
    Teuchos::RCP<WaveEquationProblem<2, float>> wave2df_;
    Teuchos::RCP<WaveEquationProblem<3, float>> wave3df_;

    bool doubleorfloat_;

    bool writebacioutput_;
  };

}  // namespace ACOU

#else  // #ifdef HAVE_DEAL_II

namespace ACOU
{
  // Dummy class that cannot work without deal.II
  class AcouExplicitTimeInt : public AcouTimeInt
  {
   public:
    /// constructor
    AcouExplicitTimeInt(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output)
        : AcouTimeInt(actdis, solver, params, output)
    {
      dserror("Explicit time integration requires interface to deal.II package");
    }

    /// set initial field by given function
    virtual void SetInitialField(int startfuncno)
    {
      dserror("Explicit time integration requires interface to deal.II package");
    }

    /// return the name of the time integration scheme
    std::string Name()
    {
      dserror("Explicit time integration requires interface to deal.II package");
      return std::string();
    }

    /// Output nodal internal pressure field at new time step
    virtual void NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec)
    {
      dserror("Explicit time integration requires interface to deal.II package");
    }

    /// run a complete simulation
    virtual void Integrate()
    {
      dserror("Explicit time integration requires interface to deal.II package");
    }

    /// output routine: not inherited since it does something completely different than the one for
    /// the implicit integration
    virtual void Output()
    {
      dserror("Explicit time integration requires interface to deal.II package");
    }

    // get values for speed of sound gradient
    double GetSoSGradient(int rowid)
    {
      dserror("Explicit time integration requires interface to deal.II package");
      return 0.0;
    }

    // get values for density gradient
    double GetDensityGradient(int rowid)
    {
      dserror("Explicit time integration requires interface to deal.II package");
      return 0.0;
    }
  };

}  // namespace ACOU

#endif  // HAVE_DEAL_II

#endif
