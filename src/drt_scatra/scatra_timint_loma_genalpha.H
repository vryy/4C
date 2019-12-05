/*----------------------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time-integration scheme with extensions for
       loma problems

\level 2

\maintainer Anh-Tu Vuong

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_LOMA_GENALPHA_H
#define SCATRA_TIMINT_LOMA_GENALPHA_H

#include "scatra_timint_loma.H"
#include "scatra_timint_genalpha.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace SCATRA
{
  class TimIntLomaGenAlpha : public virtual ScaTraTimIntLoma, public virtual TimIntGenAlpha
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// Destructor
    virtual ~TimIntLomaGenAlpha();

    /// initialize time integration scheme
    virtual void Init();

    /// setup time integration scheme
    virtual void Setup();

    /// predict thermodynamic pressure and time derivative for low-Mach-number flow
    void PredictThermPressure();

    /// compute values of thermodynamic pressure at intermediate time steps
    void ComputeThermPressureIntermediateValues();

    /// compute thermodynamic pressure and time derivative for low-Mach-number flow
    void ComputeThermPressure();

    ///  compute time derivative of thermodynamic pressure
    void ComputeThermPressureTimeDerivative();

    /// update thermodynamic pressure and time derivative for low-Mach-number flow
    void UpdateThermPressure();

    /// read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

    /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
    double ThermPressAf() { return thermpressaf_; }

    /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
    double ThermPressAm() { return thermpressam_; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_F for
    /// low-Mach-number flow
    double ThermPressDtAf() { return thermpressdtaf_; }

    /// routine to return time derivative of thermo. press. at time step n+alpha_M for
    /// low-Mach-number flow
    double ThermPressDtAm() { return thermpressdtam_; }

   protected:
    /// write additional data required for restart
    void OutputRestart() const;

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs();

    /// dynamic Vreman model
    void DynamicComputationOfCv();

    /// add thermodynamic pressure to parameter list for element evaluation
    void AddThermPressToParameterList(Teuchos::ParameterList& params  //!< parameter list
    );

   private:
    /// LOMA-specific parameter: thermodynamic pressure at n+alpha_F and n+alpha_M
    /// and its time derivative at n+alpha_F and n+alpha_M
    double thermpressaf_;
    double thermpressam_;
    double thermpressdtaf_;
    double thermpressdtam_;


  };  // class TimIntLomaGenAlpha

}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_LOMA_GENALPHA_H
