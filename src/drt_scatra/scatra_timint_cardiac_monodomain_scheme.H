/*!----------------------------------------------------------------------

\brief  connecting time-integration schemes (OST, BDF2, GenAlpha, Stationary) with
        Cardiac-monodomain-specific implementation (class TimIntCardiacMonodomain)

\level 2

\maintainer Anh-Tu Vuong

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_H
#define SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_H

#include "scatra_timint_cardiac_monodomain.H"
#include "scatra_timint_ost.H"
#include "scatra_timint_bdf2.H"
#include "scatra_timint_genalpha.H"
//#include "scatra_timint_stat.H"


namespace SCATRA
{
  class TimIntCardiacMonodomainOST : public virtual TimIntCardiacMonodomain,
                                     public virtual TimIntOneStepTheta
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainOST(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~TimIntCardiacMonodomainOST();

    //! setup time integration scheme
    virtual void Setup();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

   protected:
    //! write additional data required for restart
    void OutputRestart() const;

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false);

    //! do not calculate initial scalar time derivatives for ep
    virtual void CalcInitialTimeDerivative() { return; };

  };  // class TimIntCardiacMonodomainOST


  class TimIntCardiacMonodomainBDF2 : public virtual TimIntCardiacMonodomain,
                                      public virtual TimIntBDF2
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainBDF2(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~TimIntCardiacMonodomainBDF2();

    //! setup time integration scheme
    virtual void Setup();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

   protected:
    //! write additional data required for restart
    void OutputRestart() const;

    //! do not calculate initial scalar time derivatives for ep
    virtual void CalcInitialTimeDerivative() { return; };

  };  // class TimIntCardiacMonodomainBDF2


  class TimIntCardiacMonodomainGenAlpha : public virtual TimIntCardiacMonodomain,
                                          public virtual TimIntGenAlpha
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainGenAlpha(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~TimIntCardiacMonodomainGenAlpha();

    //! Setup time integration scheme
    virtual void Setup();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num = 0);

    //! read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

    /// add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false);

   protected:
    //! write additional data required for restart
    void OutputRestart() const;

    //! do not calculate initial scalar time derivatives for ep
    virtual void CalcInitialTimeDerivative() { return; };

  };  // class TimIntCardiacMonodomainGenAlpha

  /*
  class TimIntCardiacMonodomainStationary : public virtual TimIntCardiacMonodomain, public virtual
  TimIntStationary
  {

  public:

    //! Standard Constructor
    TimIntCardiacMonodomainStationary(
      Teuchos::RCP<DRT::Discretization>       dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

    //! Destructor
    virtual ~TimIntCardiacMonodomainStationary();

    //! initialize time integration scheme
    virtual void Init();

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update(const int num=0);

    //! read restart data
    void ReadRestart(int step,Teuchos::RCP<IO::InputControl> input = Teuchos::null);

  protected:

    //! write additional data required for restart
    void OutputRestart() const;

  }; // class TimIntCardiacMonodomainStatationary
  */

}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_H
