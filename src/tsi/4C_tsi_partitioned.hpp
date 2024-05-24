/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all TSI algorithms that perform a coupling between the
       structural field equation and temperature field equations


\level 2
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_TSI_PARTITIONED_HPP
#define FOUR_C_TSI_PARTITIONED_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_inpar_tsi.hpp"
#include "4C_tsi_algorithm.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                           dano 12/09 |
 *----------------------------------------------------------------------*/
namespace TSI
{
  //! TSI algorithm base
  //!
  //!
  //!  Base class of TSI algorithms. Derives from structure_base_algorithm and
  //!  ThermoBaseAlgorithm with temperature field.
  //!  There can (and will) be different subclasses that implement different
  //!  coupling schemes.
  //!
  //!  \warning The order of calling the two BaseAlgorithm-constructors (that
  //!  is the order in which we list the base classes) is important here! In the
  //!  constructors control file entries are written. And these entries define
  //!  the order in which the filters handle the Discretizations, which in turn
  //!  defines the dof number ordering of the Discretizations... Don't get
  //!  confused. Just always list structure, thermo. In that order.
  //!
  //!  \author u.kue
  //!  \date 02/08
  class Partitioned : public Algorithm
  {
   public:
    //! create using a Epetra_Comm
    explicit Partitioned(const Epetra_Comm& comm);


    //! outer level time loop (to be implemented by deriving classes)
    void TimeLoop() override;

    //! non-linear solve, i.e. (multiple) corrector
    void Solve() override;

    //! initialise internal variables needed as guess for the partitioned TSI algorithm
    void SetupSystem() override;

    //! time loop for TSI algorithm with one-way coupling
    void TimeLoopSequStagg();

    //! time loop for TSI algorithm with one-way coupling
    void TimeLoopOneWay();

    //! time loop for TSI algorithm with iteration between fields (full coupling)
    void TimeLoopFull();

    //! read restart data
    void read_restart(int step  //!< step number where the calculation is continued
        ) override;


   protected:
    //! @name Time loop building blocks

    //! start a new time step
    void prepare_time_step() override;

    //! calculate stresses, strains, energies
    void prepare_output() override;

    //! take current results for converged and save for next time step
    void Update() override;
    //@}

    void prepare_contact_strategy() override;

    //! @name Solve

    //! solve temperature equations for current time step
    void DoThermoStep();

    //! solve displacement equations for current time step
    void DoStructureStep();

    // two-way coupling (iterative staggered)

    //! outer iteration loop
    void OuterIterationLoop();

    //@}

    //! convergence check for iterative staggered TSI solver
    bool convergence_check(int itnum,  //!< index of current iteration
        int itmax,                     //!< maximal number of iterations
        double ittol                   //!< iteration tolerance
    );

    enum INPAR::TSI::ConvNorm normtypeinc_;  //!< convergence check for residual temperatures

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;

    //! temperature increment of the outer loop
    Teuchos::RCP<Epetra_Vector> tempincnp_;
    //! displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

   private:
    //! displacements at time step (t_n) or (t_n+1)
    Teuchos::RCP<const Epetra_Vector> disp_;
    //! velocities at time step (t_n) or (t_n+1)
    Teuchos::RCP<const Epetra_Vector> vel_;

    //! temperature at time step (t_n) or (t_n+1)
    Teuchos::RCP<Epetra_Vector> temp_;

    //! @name Aitken relaxation

    //! difference of last two solutions
    // del = r^{i+1}_{n+1} = d^{i+1}_{n+1} - d^i_{n+1}
    Teuchos::RCP<Epetra_Vector> del_;
    //! difference of difference of last two pair of solutions
    // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
    Teuchos::RCP<Epetra_Vector> delhist_;
    //! Aitken factor
    double mu_;

    //@}

    //! coupling algorithm
    INPAR::TSI::SolutionSchemeOverFields coupling_;
    //! we couple based on displacements
    bool displacementcoupling_;
    //! quasi-static solution of the mechanical equation
    bool quasistatic_;

  };  // Partitioned
}  // namespace TSI


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
