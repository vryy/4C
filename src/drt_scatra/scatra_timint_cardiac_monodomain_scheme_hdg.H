/*!----------------------------------------------------------------------
\file scatra_timint_cardiac_monodomain_scheme_hdg.H

\brief  Connecting time-integration schemes for HDG with
        cardiac-monodomain-specific implementation

<pre>
\level 3
\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_H
#define SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_H

#include "scatra_timint_cardiac_monodomain.H"
#include "scatra_timint_hdg.H"


namespace SCATRA
{

class TimIntCardiacMonodomainHDG : public virtual TimIntCardiacMonodomain, public virtual TimIntHDG
{
  public:

  //! Standard Constructor
  TimIntCardiacMonodomainHDG(
    Teuchos::RCP<DRT::Discretization>       dis,
    Teuchos::RCP<LINALG::Solver>            solver,
    Teuchos::RCP<Teuchos::ParameterList>    params,
    Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList>    extraparams,
    Teuchos::RCP<IO::DiscretizationWriter>  output);

  //! Destructor
  virtual ~TimIntCardiacMonodomainHDG();

  //! setup time integration scheme
  virtual void Setup();

  //! Update the solution after convergence of the nonlinear iteration.
  //! Current solution becomes old solution of next timestep.
  void Update(const int num=0);

  //! write current state to BINIO
  void OutputState();

protected:

  //! update gating variables
  virtual void ElementMaterialTimeUpdate();

  //! write problem specific output
  virtual void WriteProblemSpecificOutput(
      Teuchos::RCP<Epetra_Vector>    interpolatedPhi
       );

  //! adapt material
  virtual void PackMaterial();

  //! adapt material
  virtual void UnpackMaterial();

  //! project material field
  virtual void ProjectMaterial();

  //! read restart
  virtual void ReadRestart(const int step,Teuchos::RCP<IO::InputControl> input = Teuchos::null);

private:

  //! activation time
  Teuchos::RCP<Epetra_Vector> activation_time_interpol_;

  //! element data
  Teuchos::RCP<std::vector<char> > data_;


}; // class TimIntCardiacMonodomainHDG

} // namespace SCATRA

#endif /* SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_H */
