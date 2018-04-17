/*!----------------------------------------------------------------------
\file fluid_timint_loma.H

\brief TimIntLoma

\level 2
<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_LOMA_H
#define FLD_TIMINT_LOMA_H


#include "fluidimplicitintegration.H"


namespace FLD
{

class TimIntLoma : public virtual FluidImplicitTimeInt
{
  friend class TimIntGenAlpha;
public:
  /// Standard Constructor
  TimIntLoma(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntLoma();

  /*!
  \brief initialization

  */
  virtual void Init();

  /*!
  \brief set scalar fields within outer iteration loop

  */
  void SetLomaIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
                         Teuchos::RCP<const Epetra_Vector> scalaram,
                         Teuchos::RCP<const Epetra_Vector> scalardtam,
                         Teuchos::RCP<const Epetra_Vector> fsscalaraf,
                         const double             thermpressaf,
                         const double             thermpressam,
                         const double             thermpressdtaf,
                         const double             thermpressdtam,
                         Teuchos::RCP<DRT::Discretization> scatradis);

  /*!
  \brief set scalar fields

  */
  void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
                         const double             thermpressnp,
                         Teuchos::RCP<const Epetra_Vector> scatraresidual,
                         Teuchos::RCP<DRT::Discretization> scatradis,
                         const int                whichscalar = -1);

  /*!
  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element*/
  void SetElementCustomParameter();

  /*!
  \brief print turbulence model

  */
  void PrintTurbulenceModel();

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  void SetCustomEleParamsApplyNonlinearBoundaryConditions(Teuchos::ParameterList& eleparams);

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams);

  /*!
  \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

  */
  void SetCustomEleParamsLinearRelaxationSolve(Teuchos::ParameterList& eleparams);

  /*!
  \brief Call statistics manager (special case in TimIntLoma)

  */
  void CallStatisticsManager();

  /// prepare AVM3-based scale separation
  void AVM3Preparation();

  /*!
  \brief return thermpressaf_ in TimIntLoma

  */
  double ReturnThermpressaf(){return thermpressaf_;}


protected:

private:

  /// for low-Mach-number flow solver: thermodynamic pressure at n+alpha_F/n+1
  /// and at n+alpha_M/n as well as its time derivative at n+alpha_F/n+1 and n+alpha_M/n
  double thermpressaf_;
  double thermpressam_;
  double thermpressdtaf_;
  double thermpressdtam_;


}; // class TimIntLoma

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_LOMA_H
