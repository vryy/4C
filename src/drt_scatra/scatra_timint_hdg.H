/*!----------------------------------------------------------------------
\file scatra_timint_hdg.H
\brief Time integration class for HDG discretization for scalar transport

<pre>
\level 3
\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_HDG_H_
#define SCATRA_TIMINT_HDG_H_

#include "scatra_timint_genalpha.H"

namespace SCATRA
{
/*!
\brief time integration for HDG scatra
*/
class TimIntHDG : public TimIntGenAlpha
{
public:

  //! standard constructor
  TimIntHDG(
      const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<Teuchos::ParameterList>&   extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>        output
      );

  //! destructor
  virtual ~TimIntHDG() {}

  //! setup
  virtual void Setup();

  //! set theta_ to its value, dependent on integration method for GenAlpha and BDF2
  virtual void SetTheta();

  //! set states in the time integration schemes: additional vectors for HDG
  virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver=false);

  //! set the part of the right hand side belonging to the last time step
  virtual void SetOldPartOfRighthandside();

  //! update the solution after convergence of the nonlinear iteration,
  //! current solution becomes old solution of next time step
  virtual void Update(const int num=0);

//  //! Initialization procedure before the first time step is done
//  void PrepareFirstTimeStep ();

  //! update configuration and output to file/screen
  virtual void OutputState();

  //! write output for restart
  virtual void OutputRestart();

  //! read restart
  virtual void ReadRestart(const int step);

  //! set the initial scalar field phi
  void SetInitialField(
    const INPAR::SCATRA::InitialField   init,         //!< type of initial field
    const int                           startfuncno   //!< number of spatial function
  );

  //! accessor to interior concentrations
  virtual Teuchos::RCP<Epetra_Vector>  ReturnIntPhinp(){return intphinp_;}
  virtual Teuchos::RCP<Epetra_Vector>  ReturnIntPhin(){return intphin_;}
//  virtual Teuchos::RCP<Epetra_Vector>  ReturnIntPhinm(){return intphinm_;}

  virtual Teuchos::RCP<Epetra_Vector>  InterpolatedPhinp() const {return interpolatedPhinp_;}

  //! prepare time loop
  void PrepareTimeLoop();

protected:

  //! copy constructor
  TimIntHDG (const TimIntHDG& old);

  //! update time derivative for generalized-alpha time integration
  virtual void GenAlphaComputeTimeDerivative();

  //! compute values at intermediate time steps for gen.-alpha
  virtual void GenAlphaIntermediateValues() ;

  //! number of dofset for interior variables
  int   nds_intvar_;

  //! @name concentration and concentration gradient at different times for element interior for HDG
  //@{
  Teuchos::RCP<Epetra_Vector> intphinp_;   //!< concentration at time \f$t^{n+1}\f$
  Teuchos::RCP<Epetra_Vector> intphin_;    //!< concentration at time \f$t^{n}\f$
//  Teuchos::RCP<Epetra_Vector> intphinm_;   //!< concentration at time \f$t^{n-1}\f$
//  Teuchos::RCP<Epetra_Vector> intphiaf_;   //!< concentration at time \f$t^{n+\alpha_F}\f$
//  Teuchos::RCP<Epetra_Vector> intphiam_;   //!< concentration at time \f$t^{n+\alpha_M}\f$
  //@}

  //! @name scalar time derivative of concentration and concentration gradient
  //! at time n+1, n and n+alpha_M/(n+alpha_M/n) and n-1 for element interior in HDG
  //@{
//  Teuchos::RCP<Epetra_Vector> intphidtnp_;   //!< time derivative at time \f$t^{n+1}\f$
//  Teuchos::RCP<Epetra_Vector> intphidtn_;    //!< time derivative at time \f$t^{n}\f$
//  Teuchos::RCP<Epetra_Vector> intphidtnm_;   //!< time derivative at time \f$t^{n-1}\f$
//  Teuchos::RCP<Epetra_Vector> intphidtam_;   //!< time derivative at time \f$t^{n+\alpha_M}\f$
  //@}

  //! @name other HDG-specific auxiliary vectors
  //@{
  Teuchos::RCP<Epetra_Vector> interpolatedPhinp_;  //!< concentrations for output at time \f$t^{n+1}\f$
  //@}

  //! calculate intermediate solution
  virtual void ComputeIntermediateValues();

  //! compute values at the interior of the elements
  virtual void ComputeInteriorValues();

  //! update interior variables
  virtual void UpdateInteriorVariables(
      Teuchos::RCP<Epetra_Vector>          updatevector
      );

  //! write problem specific output
  virtual void WriteProblemSpecificOutput(
      Teuchos::RCP<Epetra_Vector>    interpolatedPhi
       ) {return;}

  //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
  void CalcInitialTimeDerivative() {return;};

  //! fd check
  virtual void FDCheck();

  //! adapt degree of test functions and change dofsets accordingly
  virtual void AdaptDegree();

  //! adapt variable vectors required due to the change of the degrees of the test functions
  virtual void AdaptVariableVector(
      Teuchos::RCP<Epetra_Vector>                phi_new,
      Teuchos::RCP<Epetra_Vector>                phi_old,
      Teuchos::RCP<Epetra_Vector>                intphi_new,
      Teuchos::RCP<Epetra_Vector>                intphi_old,
      int                                        nds_var_old,
      int                                        nds_intvar_old,
      std::vector<DRT::Element::LocationArray>   la_old,
      const Teuchos::RCP<Epetra_Map>&            dofrowmap_old,
      const Teuchos::RCP<Epetra_Map>&            intdofrowmap_old
      );

  //! calculate matrices on element
  virtual void CalcMatInitial();

  //! chooses the assemble process (assemble matrix and rhs or only rhs)
  void AssembleMatAndRHS();

  //! contains the assembly process only for rhs
  void AssembleRHS();

private:

  //! time algorithm flag actually set (we internally reset it)
  INPAR::SCATRA::TimeIntegrationScheme timealgoset_;

  //! @name time stepping variable
  bool   startalgo_;   //!< flag for starting algorithm

  //! time-integration-scheme factor theta
  double theta_;

  //! activation_time at times n+1
  Teuchos::RCP<Epetra_Vector> activation_time_interpol_np_;

  //! HDG discretization
  DRT::DiscretizationHDG* hdgdis_;

  //! p-adativitity
  bool padaptivity_;

  //! error tolerance for p-adativitity
  double padapterrortol_;

  //! error tolerance base for p-adativitity (delta p=log_padapterrorbase_(E/padapterrortol_))
  double padapterrorbase_;


  //! max. degree of shape functions for p-adativitity
  double padaptdegreemax_;

  //! element degree
  Teuchos::RCP<Epetra_Vector> elementdegree_;

}; // class TimIntHDG
} // namespace SCATRA

#endif // #ifndef SCATRA_TIMINT_HDG_H
