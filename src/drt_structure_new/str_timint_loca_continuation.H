/*-----------------------------------------------------------*/
/*!
\file str_timint_loca_continuation.H

\brief Library of Continuation Algorithms (LOCA) implementation.

\maintainer Michael Hiermeier

\date Nov 18, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_TIMINT_LOCA_CONTINUATION_H_
#define SRC_DRT_STRUCTURE_NEW_STR_TIMINT_LOCA_CONTINUATION_H_


#include "str_timint_implicitbase.H"  // base class

#include "../loca_continuation/loca_nln_stepper.H"  // enum
#include "../solver_nonlin_nox/nox_nln_enum_lists.H" // enum

//forward declaration
namespace LOCA {
  class ParameterVector;
namespace MultiContinuation {
  class ConstraintInterface;
} // namespace MultiContinuation
namespace Epetra {
  class Factory;
} // namespace Epetra
} // namespace LOCA
namespace STR {
namespace TIMINT {
  /** \brief LOCA continuation strategy
   *
   * \author Michael Hiermeier */
  class LOCAContinuation : public ImplicitBase
  {
  public:
    //! constructor
    LOCAContinuation();

    //! destructor
    virtual ~LOCAContinuation() {};

    //! derived setup routine
    virtual void Setup();

    //! call of the stepper run() routine (derived)
    virtual int Integrate();

    /// integrate the current step (derived)
    virtual int IntegrateStep();

    ///  Evaluate routine for coupled problems with monolithic approach
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> disiterinc);

    /// don't update displacement but evaluate elements (implicit only)
    virtual void Evaluate() ;

    /// set the state of the nox group and the global state data container
    virtual void SetState(const Teuchos::RCP<Epetra_Vector> & x);

    /*! \brief nonlinear solve
     *
     *  Do the nonlinear solve, i.e. (multiple) corrector,
     *  for the time step. All boundary conditions have
     *  been set. */
    virtual INPAR::STR::ConvergenceStatus Solve();

    //! start new time step (derived)
    virtual void PrepareTimeStep();

    //! start new partitioned step
    virtual void PreparePartitionStep();


  protected:
    //! Returns the current solution group
    virtual const NOX::Abstract::Group& GetSolutionGroup() const;

    //! Check for necessary parameters which have to be supplied by the user.
    void CheckForNecessaryParameters() const;

    //! combine the LOCA and NOX parameter lists
    void CombineParameterLists();

    //! create the LOCA parameter vector (continuation informations)
    Teuchos::RCP<LOCA::ParameterVector> CreateParameterVector(const Teuchos::ParameterList& p_loca);

  protected:
    //! parameter list which holds the loca and the nox parameters
    Teuchos::RCP<Teuchos::ParameterList> p_loca_nox_ptr_;

    //! LOCA constraints pointer (e.g. arc-length constraints object)
    Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> loca_constraints_ptr_;

    Teuchos::RCP< std::vector<std::string> > loca_constraint_param_names_ptr_;

    //! LOCA stepper (doing the main work)
    Teuchos::RCP<LOCA::NLN::Stepper> loca_stepper_ptr_;

    //! LOCA stepper status
    LOCA::Abstract::Iterator::IteratorStatus loca_stepper_status_;
  };  // class LOCAContinuation
} // namespace TIMINT
} // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_TIMINT_LOCA_CONTINUATION_H_ */
