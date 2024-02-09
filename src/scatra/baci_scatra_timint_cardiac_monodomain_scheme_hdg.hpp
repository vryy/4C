/*----------------------------------------------------------------------*/
/*! \file

\brief  Connecting time-integration schemes for HDG with
        cardiac-monodomain-specific implementation

\level 3

*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HPPDG_HPP
#define BACI_SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HPPDG_HPP

#include "baci_config.hpp"

#include "baci_scatra_timint_cardiac_monodomain.hpp"
#include "baci_scatra_timint_hdg.hpp"

BACI_NAMESPACE_OPEN


namespace SCATRA
{
  class TimIntCardiacMonodomainHDG : public virtual TimIntCardiacMonodomain,
                                     public virtual TimIntHDG
  {
   public:
    //! Standard Constructor
    TimIntCardiacMonodomainHDG(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    //! setup time integration scheme
    void Setup() override;

    //! Update the solution after convergence of the nonlinear iteration.
    //! Current solution becomes old solution of next timestep.
    void Update() override;

    //! write current state to BINIO
    void OutputState() override;

   protected:
    //! update gating variables
    virtual void ElementMaterialTimeUpdate();

    //! write problem specific output
    void WriteProblemSpecificOutput(Teuchos::RCP<Epetra_Vector> interpolatedPhi) override;

    //! adapt material
    void PackMaterial() override;

    //! adapt material
    void UnpackMaterial() override;

    //! project material field
    void ProjectMaterial() override;

    //! read restart
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override;

   private:
    //! activation time
    Teuchos::RCP<Epetra_Vector> activation_time_interpol_;

    //! element data
    Teuchos::RCP<std::vector<char>> data_;


  };  // class TimIntCardiacMonodomainHDG

}  // namespace SCATRA

BACI_NAMESPACE_CLOSE

#endif  // SCATRA_TIMINT_CARDIAC_MONODOMAIN_SCHEME_HDG_H
