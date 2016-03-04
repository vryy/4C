/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_springdashpot.H

\maintainer Martin Pfaller

\date Feb 29, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_SPRINGDASHPOT_H
#define STR_MODEL_EVALUATOR_SPRINGDASHPOT_H

#include "str_model_evaluator_generic.H"
#include "../drt_constraint/springdashpot_new.H"
#include <Teuchos_RCP.hpp>

//forward declaration ...
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace STR{
namespace MODELEVALUATOR{

  class SpringDashpot : public Generic
  {
  public:

    //! constructor
    SpringDashpot();

    //! destructor
    virtual ~SpringDashpot(){};

    virtual void Setup();

    //! derived
    virtual bool ApplyForce(
        const Epetra_Vector& x,
        Epetra_Vector& f);

    //! derived
    virtual bool ApplyStiff(
        const Epetra_Vector& x,
        LINALG::SparseOperator& jac);

    //! derived
    virtual bool ApplyForceStiff(
        const Epetra_Vector& x,
        Epetra_Vector& f,
        LINALG::SparseOperator& jac);

    //! derived
    virtual void UpdateStepState();

    //! derived
    virtual void UpdateStepElement();

    //! derived
    virtual void OutputStepState();

    void OutputRestart();

//    void ReadRestart(
//        IO::DiscretizationReader& reader,
//        const double& time);

  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);
    //! derived
    virtual void Reset(const Epetra_Vector& x, LINALG::SparseOperator& jac);

  private:
    //! all spring dashpot instances
    std::vector<Teuchos::RCP<UTILS::SpringDashpotNew> > springs_;

    //! number of spring dashpot conditions
    int n_conds_;

    //! structural displacement at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

    //! structural stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr_;
  };

} // namespace ModelEvaluator
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_SPRINGDASHPOT_H */
