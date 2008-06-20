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
*/
class StruTimIntImpl : public StruTimInt
{

public:

  //! NOX interface is friend of this class
  friend class NoxInterface;

  //! Constructor
  StruTimIntImpl
  (
    const Teuchos::ParameterList& sdynparams,  //!< input parameters
    DRT::Discretization& actis,  //!< current discretisation
    LINALG::Solver& solver,  //!< the solver
    IO::DiscretizationWriter& output  //!< the output
  );

  //! Destructor
  ~StruTimIntImpl() { }

  //! Initialise system:
  //! equilibrate the initial state by identifying the consistent 
  //! initial accelerations and (if applicable) internal variables
  void Initialise() { }

  //! Predict target solution
  void Predict() { }

  //! Apply Dirichlet boundary conditions
  void ApplyDirichletBC() { }

  /*
  //! Compute the residual to a given state x
  //!
  //! This method is used by the NOX user interface class 
  //! and by the linalg matrix free operator
  void computeF
  (
    const Epetra_Vector& x,  //!< Current solution state 
                             //!< (Dirichlet DOFs will be blanked)
    Epetra_Vector& F         //!< residual matching current solution state 
                             //!< (Dirichlet DOFs will be blanked)
  );
  */

  //! Do residual force due to global balance of momentum
  virtual void ResidualForce() { }

  //! Do tangent of residual force with respect to displacements
  virtual void Tangent() { }

  //! Test convergence of iterative solution technique
  bool Converged() { return true; }
  /*
  (
    const string type,  //!< type
    const double disinorm,  //!< current displacement norm
    const double resnorm,  //!< current force residual norm
    const double toldisp,  //!< displacement tolerance
    const double tolres  //!< force residual tolerance
  )
  */

  /*
  //! Test convergence of iterative solution technique
  bool Converged
  (
    const string type, 
    const double disinorm,
    const double resnorm,
    const double volnorm,
    const double toldisp, 
    const double tolres,
    const double tolvol
  );
  */

  //! Solve dynamic equilibrium
  void Solve() { }

  //! Do full Newton-Raphson iteration
  void NewtonFull() { }

  //! Do (so-called) modified Newton-Raphson iteration in which
  //! the initial tangent is kept and not adapted to the current
  //! state of the displacement solution
  void NewtonModified() { }

  //! Do classical augmented lagrange for volume constraint
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

  //! Do time integration of multiple steps
  void Integrate()
  { 
    cout << "You are here." << endl
         << "You have a time step size of " << dt_ << endl
         << endl;
  }

  //! Do time integration of single step
  void IntegrateStep() { }

  //! update configuration
  virtual void Update() { }

  //! output to file/screen
  void Output() { }

protected:

  //! type of solution techniques
  typedef enum
  {
    vague,  //!< undefined
    newtonfull,  //!< full Newton-Raphson iteration
    newtonmod  //!< modified Newton-Raphson iteration
  } SolTech;

  //! equal operator is NOT wanted
  StruTimIntImpl operator = (const StruTimIntImpl& old);

  //! copy constructor is NOT wanted
  StruTimIntImpl(const StruTimIntImpl& old);

  //! @name General purpose algorithm parameters
  //@{
  DRT::Discretization& discret_;  //!< attached discretisation
  LINALG::Solver& solver_;  //!< linear algebraic solver
  IO::DiscretizationWriter& output_;  //!< binary output
  int myrank_;  //!< ID of actual processor in parallel
  //@}

  //! @name General control parameters
  //@{
  double time_;  //!< time \f$t_{n}\f$ of last converged step
  double timen_;  //!< target time \f$t_{n+1}\f$
  double dt_;  //!< time step size \f$\Delta t\f$
  double timemax_;  //!< final time \f$t_\text{fin}\f$
  int step_;  //!< time step index \f$n\f$
  int stepmax_;  //!< final step \f$N\f$
  //@}

  //! whatever constraints
  RCP<ConstrManager> constrman_;

  //! surface stresses
  RCP<DRT::SurfStressManager> surfstressman_;

  //! @name Rayleigh damping
  //!       \f$\boldsymbol{C} = c_\text{K} \boldsymbol{K} + c_\text{M} \boldsymbol{M}\f$
  //@{
  bool damping_;  //!< damping active/inactive
  double dampk_;  //!< damping factor for stiffness \f$c_\text{K}\f$
  double dampm_;  //!< damping factor for mass \f$c_\text{M}\f$
  //@}

  //! @name Iterative solution technique
  //@{
  //int itertype_;  //!< kind of iteration technique
  int iter_;  //!< iteration step
  int itermax_;  //!< maximal permitted iterations
  double toldis_;  //!< tolerance residual displacements
  double tolres_;  //!< tolerance force residual
  double tolcon_;  //!< tolerance constraint
  //@}
 
  //! @name Global vectors
  //@{
  RCP<Epetra_Vector> dirichtoggle_;  //!< At a free DOF this vector
                                     //!< holds a 0; at a supported DOF
                                     //!< it is a 1.
  RCP<Epetra_Vector> invtoggle_;  //!< Opposite of \c dirichtoggle_ vector:
                                  //!< i.e.\ for each component 
                                  //!< a 0 indicates a supported 
                                  //!< corresponding DOF; a 1 indicates
                                  //!< a free/unsupported DOF.
  RCP<Epetra_Vector> zeros_;  //!< a zero vector of full length
  RCP<Epetra_Vector> dis_;  //!< global displacements
                            //!< \f$\boldsymbol{D}_{n}\f$  at \f$t_n\f$
  RCP<Epetra_Vector> vel_;  //!< global velocities
                            //!< \f$\boldsymbol{V}_{n}\f$ at \f$t_n\f$
  RCP<Epetra_Vector> acc_;  //!< global accelerations
                            //!< \f$\boldsymbol{A}_{n}\f$ at \f$t_n\f$
  RCP<Epetra_Vector> disn_;  //!< global displacements
                             //!< \f$\boldsymbol{D}_{n+1}\f$  
                             //!< at \f$t_{n+1}\f$
  RCP<Epetra_Vector> veln_;  //!< global velocities
                             //!< \f$\boldsymbol{V}_{n+1}\f$  
                             //!< at \f$t_{n+1}\f$
  RCP<Epetra_Vector> accn_;  //!< global accelerations
                             //!< \f$\boldsymbol{A}_{n+1}\f$  
                             //!< at \f$t_{n+1}\f$
  RCP<Epetra_Vector> disi_;  //!< residual displacements
                             //!< \f$\Delta\boldsymbol{D}^{<k>}_{n+1}\f$
  RCP<Epetra_Vector> fint_;  //!< internal force vector
                             //!< \f$\boldsymbol{f}_\text{int;n}\f$
                             //!< at \f$t_{n}\f$
  RCP<Epetra_Vector> fintn_;  //!< internal forces at \f$t_{n+1}\f$
  RCP<Epetra_Vector> finert_;  //!< inertial forces
  RCP<Epetra_Vector> fvisc_;  //!< viscous Rayleigh forces
  RCP<Epetra_Vector> fext_;  //!< external force at \f$t_{n+1}\f$
  RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_n\f$
  RCP<Epetra_Vector> fres_;  //!< force residual
  RCP<Epetra_Vector> frobin_;
  //@}

  //! @name Global matrices
  //@{
  RCP<LINALG::SparseMatrix> stiff_;
  RCP<LINALG::SparseMatrix> mass_;
  RCP<LINALG::SparseMatrix> damp_;
  //@}
};

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_IMPL_H
#endif  // #ifdef CCADISCRET
