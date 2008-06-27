/*----------------------------------------------------------------------*/
/*!
\file strutimint_impl.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_IMPL_H
#define STRUTIMINT_IMPL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint.H"

/*======================================================================*/
/*!
\brief Front-end for structural dynamics 
       by implicit structure time integration. 
\author bborn
\date 06/08
*/
class StruTimIntImpl : public StruTimInt
{

public:

  //! NOX interface is friend of this class
  //friend class NoxInterface;

  //! type of solution techniques
  enum SolTechEnum
  {
    soltech_vague,  //!< undefined
    soltech_newtonfull,  //!< full Newton-Raphson iteration
    soltech_newtonmod,  //!< modified Newton-Raphson iteration
    soltech_newtonuzawalin,  //!< linear Uzawa iteration for 
                             //!< constraint system
    soltech_newtonuzawanonlin  //!< non-linear Uzawa iteration 
                               //!< for constraint system
  };

  //! map solution technique identification string to enum term
  static enum SolTechEnum MapSolTechStringToEnum
  (
    const std::string name  //!< name identification string
  );

  //! kind of iteration convergence checks
  //! \author lw \date 12/07
  enum ConvCheckEnum
  {
    convcheck_absres_or_absdis,  //!< absolute norms of residual forces
                                 //!< OR iterative displacement increments
    convcheck_absres_and_absdis,  //!< absolute norms of residual forces
                                  //!< AND iterative displacement increments
    convcheck_relres_or_absdis,  //!< relative norm of residual forces
                                 //!< OR absolute norm if iterative
                                 //!< displacement increments
    convcheck_relres_and_absdis,  //!< relative norm of residual forces
                                  //!< AND absolute norm if iterative
                                  //!< displacement increments
    convcheck_relres_or_reldis,  //!< relative norms of residual forces
                                 //!< OR iterative displacement increments
    convcheck_relres_and_reldis,  //!< relative norms of residual forces
                                  //!< AND iterative displacement increments
    convcheck_vague  //!< undetermined
  };

  //! map convergence check to enum term
  static enum ConvCheckEnum MapConvCheckStringToEnum
  (
    const std::string name  //!< name identification string
  );

  //! kind of vector norms available for convergence check
  enum VectorNormEnum
  {
    vectornorm_vague,  //!< undetermined
    vectornorm_l1,  //!< L1/linear/average norm
    vectornorm_l2,  //!< L2/Euclidean norm
    vectornorm_inf  //!< Maximum/infinity norm
  };

  //! return norm of vector according to requested norm
  static double EvaluateNorm
  (
    const enum VectorNormEnum norm,  //!< requested norm
    const Teuchos::RCP<Epetra_Vector> vect  //!< input vector
  );

  //! Type of predictor
  enum PredEnum
  {
    pred_vague,  //!< undetermined
    pred_constdis,  //!< constant displacements
    pred_constdisvelacc  //!< constant displacements, 
                         //!< velocities and accelerations
  };

  //! Map enum term to input string
  static enum PredEnum MapPredictorStringToEnum
  (
    const std::string name  //!< identifier
  );

  //! Type of stress
  enum StressEnum
  {
    stress_none,  //!< undetermined
    stress_cauchy,  //!< Cauchy stress
    stress_pk1,  //!< 1st Piola-Kirchhoff stress
    stress_pk2  //! 2nd Piola-Kirchhoff stress
  };

  //! Map stress input string to enum
  static enum StressEnum MapStressStringToEnum
  (
    const std::string name  //!< identifier
  );

  //! Type of strain
  enum StrainEnum
  {
    strain_none,  //!< undetermined
    strain_ea,  //!< Euler-Almansi strain
    strain_gl  //!< Green-Lagrange strain
  };

  //! Map strain input string to enum
  static enum StrainEnum MapStrainStringToEnum
  (
    const std::string name  //!< identifier
  );

  //! Constructor
  StruTimIntImpl
  (
    const Teuchos::ParameterList& ioparams,  //!< ioflags
    const Teuchos::ParameterList& sdynparams,  //!< input parameters
    const Teuchos::ParameterList& xparams,  //!< extra flags
    DRT::Discretization& actis,  //!< current discretisation
    LINALG::Solver& solver,  //!< the solver
    IO::DiscretizationWriter& output  //!< the output
  );

  //! Destructor
  ~StruTimIntImpl() { ; }

  //! Equilibrate the initial state by identifying the consistent 
  //! initial accelerations and (if applicable) internal variables
  //! Make damping and mass matrix
  void DetermineMassDampConsistAccel();

  //! Do time integration of multiple steps
  void Integrate();

  //! Do time integration of single step
  void IntegrateStep();

  //! @name Prediction
  //@{

  //! Predict target solution and identify residual
  void Predict();

  //! Predict constant displacements, velocities and accelerations,
  //! i.e. the initial guess is equal to the last converged step
  //! except Dirichlet BCs
  void PredictConstDisVelAcc();

  //! Predict constant displacements, however the velocities
  //! and accelerations are consistent to the time integration
  //! if the constant displacements are taken as correct displacement
  //! solution.
  //! This method has to be implemented by the individual time
  //! integrator.
  virtual void PredictConstDisConsistVelAcc() = 0;

  //@}

  //! Apply Dirichlet boundary conditions on provided state vectors
  void ApplyDirichletBC
  (  
    const double time,  //!< at time
    Teuchos::RCP<Epetra_Vector>& dis,  //!< displacements (may be Teuchos::null) 
    Teuchos::RCP<Epetra_Vector>& vel,  //!< velocities (may be Teuchos::null) 
    Teuchos::RCP<Epetra_Vector>& acc  //!< accelerations (may be Teuchos::null)
  );

  //! @name Forces and stiffnesses
  //@{

  //! Do residual force due to global balance of momentum
  //!
  //! The residual #fres_ is determined readily negative
  virtual void EvaluateForceStiffResidual() = 0;

  //! Apply external force
  void ApplyForceExternal
  (
    const double time,  //!< evaluation time
    const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
    Teuchos::RCP<Epetra_Vector>& fext  //!< external force
  );

  //! Evaluate ordinary internal force, its stiffness at state
  void ApplyForceStiffInternal
  (
    const double time,  //!< evaluation time
    const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
    const Teuchos::RCP<Epetra_Vector> disi,  //!< residual displacements
    Teuchos::RCP<Epetra_Vector>& fint,  //!< internal force
    Teuchos::RCP<LINALG::SparseMatrix>& stiff  //!< stiffness matrix
  );

  //! Evaluate forces due to constraints
  void ApplyForceStiffConstraint
  (
    const double time,
    const Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<Epetra_Vector>& fint,
    Teuchos::RCP<LINALG::SparseMatrix>& stiff
  );

  //! Evaluate _certain_ surface stresses and stiffness
  //!
  //! The Evaluation happens internal-force like.
  void ApplyForceStiffSurfstress
  (
    const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
    Teuchos::RCP<Epetra_Vector>& fint,  //!< internal force
    Teuchos::RCP<LINALG::SparseMatrix>& stiff  //! stiffness 
  );

  //! evaluate certain potential forces and its stiffness
  //!
  //! Evaluation happens internal-force like.
  //! Which potentials are actually meant can be found 
  //! at DRT::PotentialManager
  void ApplyForceStiffPotential
  (
    const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state
    Teuchos::RCP<Epetra_Vector>& fint,  //!< internal force vector
    Teuchos::RCP<LINALG::SparseMatrix>& stiff  //!< stiffness
  );

  //@}

  //! @name Solution
  //@{

  //! determine characteristic norms for relative
  //! error checks of residual displacements
  //! \author lw  \date 12/07
  virtual double CalcRefNormDisplacement() = 0;

  //! determine characteristic norms for relative
  //! error checks of residual forces
  //! \author lw  \date 12/07
  virtual double CalcRefNormForce() = 0;

  //! Determine norm of force residual
  static double CalculateNorm
  (
    const enum VectorNormEnum norm,  //!< norm to use, see #VectorNormEnum
    const RCP<Epetra_Vector> vect  //!< the vector
  );

  //! Test convergence of iterative solution technique
  //! \author lw  \date 12/07
  bool Converged();

  //! Solve dynamic equilibrium
  //!
  //! This is a general wrapper around the specific techniques.
  void Solve();

  //! Do full Newton-Raphson iteration
  //! 
  //! This routines expects a prepared negative reisdual force #fres_
  //! and associated effective stiffness matrix #stiff_
  void NewtonFull();

  //! Do (so-called) modified Newton-Raphson iteration in which
  //! the initial tangent is kept and not adapted to the current
  //! state of the displacement solution
  void NewtonModified() { }

  //! Do classical augmented lagrange for volume constraint
  //! DEAR SIRS AND MADAMS, PLEASE PROVIDE REFERENCE
  void NewtonFullUzawaNonLinear() { }

  //! do full Newton iteration respecting volume constraint
  void NewtonFullUzawaLinear() { }

  //! Do non-linear CG iteration using NOX
  void NonlinearCG() { }

  //! Do pseudo transient continuation non-linear iteration
  //! 
  //! Pseudo transient continuation is a variant of full newton that has a
  //! larger convergence radius than newton and is therefore more stable 
  //! and/or can do larger time steps
  void PTC() { }

  //! Update iteration
  //!
  //! This handles the iterative update of the current 
  //! displacements \f$D_{n+1}\f$ with the residual dispplacements
  //! The velocities and accelerations follow on par.
  virtual void UpdateIteration() = 0;

  //@}


  //! Update configuration after time step
  //! This means, the state set
  //! \f$D_{n} := D_{n+1}\f$, \f$V_{n} := V_{n+1}\f$ and \f$A_{n} := A_{n+1}\f$
  //! Thus the 'last' converged is lost and a reset of the time step becomes
  //! impossible. We are ready and keen awating the next time step.
  virtual void UpdateStep() = 0;

  //! Update step for surface stress
  void UpdateStepSurfstress()
  {
    if (surfstressman_ != Teuchos::null) surfstressman_->Update();
  }

  //! Update step for potential forces
  void UpdateStepPotential()
  {
    if (potman_ != Teuchos::null) potman_->Update();
  }

  //! @name Output
  //@{

  //! Print to screen predictor informations about residual norm etc.
  //! \author lw (originally) \date 12/07
  void PrintPredictor();

  //! Print to screen information about residual forces and displacements
  //! \author lw (originally) \date 12/07
  void PrintNewtonIter();

  //! Contains text to PrintNewtonIter
  //! \author lw (originally) \date 12/07
  void PrintNewtonIterText
  (
    FILE* ofile  //!< output file handle
  );

  //! print statistics of converged Newton-Raphson iteration
  void PrintNewtonConv();

  //! print summary after step
  void PrintStep();

  //! The text for summary print, see #PrintStep
  void PrintStepText
  (
    FILE* ofile  //!< output file handle
  );

  //! Output to file
  //! This routine always prints the last converged state, i.e.
  //! \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
  //! upon object prior to writing stuff here.
  //! \author mwgee (originally) \date 03/07
  void Output();

  //! Write restart
  //! \author mwgee (originally) \date 03/07
  void OutputRestart
  (
    bool& datawritten  //!< (in/out) read and append if
                       //!< it was written at this time step
  );

  //! Output displacements, velocities and accelerations
  //! \author mwgee (originally) \date 03/07
  void OutputState
  (
    bool& datawritten  //!< (in/out) read and append if
                       //!< it was written at this time step
  );

  //! Stress & strain output
  //! \author lw (originally)
  void OutputStressStrain
  (
    bool& datawritten  //!< (in/out) read and append if
                       //!< it was written at this time step
  );

  //@}


  //! @name Attribute access functions
  //@{
  Teuchos::RCP<DRT::Discretization> Discretization() 
  { 
    return Teuchos::rcp(&discret_, false); 
  }
  //@}

protected:

  //! equal operator is NOT wanted
  //StruTimIntImpl operator = (const StruTimIntImpl& old);

  //! copy constructor is NOT wanted
  //StruTimIntImpl(const StruTimIntImpl& old);

  //! @name General purpose algorithm parameters
  //@{
  DRT::Discretization& discret_;  //!< attached discretisation
  int myrank_;  //!< ID of actual processor in parallel
  const Epetra_Map* dofrowmap_;  //!< vector layout from the discretization 
                           //!< to construct matching
                           //!< vectors and matrices
  LINALG::Solver& solver_;  //!< linear algebraic solver
  bool solveradapttol_;  //!< adapt solver tolerance
  double solveradaptolbetter_;  //!< tolerance to which is adapated ????
  enum PredEnum pred_;  //!< predictor
  //@}

  //! @name Printing and output
  //@{
  IO::DiscretizationWriter& output_;  //!< binary output
  bool printscreen_;  //!< print infos to standard out
  FILE* errfile_;  //!< error file handle
  bool printerrfile_;  //!< print infos to error file #errfile_
  bool printiter_;  //!< print intermediate iterations during solution
  int writerestartevery_;  //!< write restart every given step;
                      //!< if 0, restart is not written
  bool writestate_;  //!< write state on/off
  int writestateevery_;  //!< write state every given step
  int writestrevery_;  //!< write stress/strain every given step
  enum StressEnum writestress_;  //!< stress output type
  enum StrainEnum writestrain_;  //!< strain output type
  //@}

  //! whatever constraints
  Teuchos::RCP<ConstrManager> constrman_;  //!< constraint manager
  Teuchos::RCP<UzawaSolver> uzawasolv_;  //!< Uzawa constraint solver

  //! surface stresses
  Teuchos::RCP<DRT::SurfStressManager> surfstressman_;

  //! potential conditions
  Teuchos::RCP<DRT::PotentialManager> potman_;

  //! @name Rayleigh damping
  //!       \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
  //@{
  bool damping_;  //!< damping active/inactive
  double dampk_;  //!< damping factor for stiffness \f$c_\text{K}\f$
  double dampm_;  //!< damping factor for mass \f$c_\text{M}\f$
  //@}

  //! @name General control parameters
  //@{
  double time_;  //!< time \f$t_{n}\f$ of last converged step
  double timen_;  //!< target time \f$t_{n+1}\f$
  double dt_;  //!< time step size \f$\Delta t\f$
  double timemax_;  //!< final time \f$t_\text{fin}\f$
  int stepmax_;  //!< final step \f$N\f$
  int step_;  //!< time step index \f$n\f$
  int stepn_;  //!< time step index \f$n+1\f$
  //@}

  //! @name Iterative solution technique
  //@{
  enum SolTechEnum itertype_;  //!< kind of iteration technique
  enum ConvCheckEnum itercnvchk_;  //!< convergence check
  enum VectorNormEnum iternorm_;  //!< vector norm to check with
  int itermax_;  //!< maximal permitted iterations
  double toldisi_;  //!< tolerance residual displacements
  double tolfres_;  //!< tolerance force residual
  double tolcon_;  //!< tolerance constraint
  int iter_;  //!< iteration step
  double normcharforce_;  //!< characteristic norm for residual force
  double normchardis_;  //!< characteristic norm for residual displacements
  double normfres_;  //!< norm of residual forces
  double normdisi_;  //!< norm of residual displacements
  double normcon_;  //!< norm of constraint
  Teuchos::RCP<Epetra_Vector> disi_;  //!< residual displacements
                                      //!< \f$\Delta{D}^{<k>}_{n+1}\f$
  Epetra_Time timer_;  //!< timer for solution technique
  //@}
 
  //! @name Global vectors
  //@{
  Teuchos::RCP<Epetra_Vector> dirichtoggle_;  //!< At a free DOF this vector
                                     //!< holds a 0; at a supported DOF
                                     //!< it is a 1.
  Teuchos::RCP<Epetra_Vector> invtoggle_;  //!< Opposite of 
                                           //!< #dirichtoggle_ vector:
                                           //!< i.e.\ for each component 
                                           //!< a 0 indicates a supported 
                                           //!< corresponding DOF; a 1 indicates
                                           //!< a free/unsupported DOF.
  Teuchos::RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
  //@}

  //! @name Global state vectors
  //! (Allright, accelerations are not truely state vectors.
  //! Moreover, the state is also constituted by the elementwise stored 
  //! hybrid/mixed quantities.)
  //@{
  Teuchos::RCP<Epetra_Vector> dis_;  //!< global displacements
                                     //!< \f${D}_{n}\f$  at \f$t_n\f$
  Teuchos::RCP<Epetra_Vector> vel_;  //!< global velocities
                                     //!< \f${V}_{n}\f$ at \f$t_n\f$
  Teuchos::RCP<Epetra_Vector> acc_;  //!< global accelerations
                                     //!< \f${A}_{n}\f$ at \f$t_n\f$
  Teuchos::RCP<Epetra_Vector> disn_;  //!< global displacements
                                      //!< \f${D}_{n+1}\f$  
                                      //!< at \f$t_{n+1}\f$
  Teuchos::RCP<Epetra_Vector> veln_;  //!< global velocities
                                      //!< \f${V}_{n+1}\f$  
                                      //!< at \f$t_{n+1}\f$
  Teuchos::RCP<Epetra_Vector> accn_;  //!< global accelerations
                                      //!< \f${A}_{n+1}\f$  
                                      //!< at \f$t_{n+1}\f$
  //@}
  
  //! @name Various global forces
  //@{
  //Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force vector
  //                           //!< \f${f}_\text{int;n}\f$
  //                            //!< at \f$t_{n}\f$
  //Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal forces at \f$t_{n+1}\f$
  //Teuchos::RCP<Epetra_Vector> finert_;  //!< inertial forces
  //Teuchos::RCP<Epetra_Vector> fvisc_;  //!< viscous Rayleigh forces
  //Teuchos::RCP<Epetra_Vector> fext_;  //!< external force at \f$t_{n+1}\f$
  //Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_n\f$
  //Teuchos::RCP<Epetra_Vector> frobin_;  //!< Robin forces due to I-DON'T-KNOW-WHAT
  Teuchos::RCP<Epetra_Vector> fres_;  //!< force residual used for solution
  //@}

  //! @name System matrices
  //@{
  Teuchos::RCP<LINALG::SparseMatrix> stiff_;  //!< holds eventually effective stiffness
  Teuchos::RCP<LINALG::SparseMatrix> mass_;  //!< mass matrix (constant)
  Teuchos::RCP<LINALG::SparseMatrix> damp_;  //!< Rayleigh damping matrix (constant)
  //@}

};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_IMPL_H
#endif  // #ifdef CCADISCRET
