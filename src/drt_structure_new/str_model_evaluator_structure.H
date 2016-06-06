/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_structure.H

\brief Evaluation and assembly of all structure terms

\maintainer Michael Hiermeier

\date Aug 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_STRUCTURE_H_
#define STR_MODEL_EVALUATOR_STRUCTURE_H_

#include "str_model_evaluator_generic.H"       // base class
#include "../drt_inpar/inpar_structure.H"      // enumerators
#include <Teuchos_RCP.hpp>

#include "str_elements_paramsinterface.H"            // interface to the element evaluation

//forward declaration ...
namespace LINALG {
  class SparseMatrix;
} // namespace LINALG
namespace STR{
namespace MODELEVALUATOR{

  class Structure : public Generic
  {
  public:
    //! constructor
    Structure();

    //! destructor
    virtual ~Structure(){};

    virtual void Setup();

    //! @name Derived public STR::MODELEVALUATOR::Generic methods
    //! @{
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
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! derived
    virtual void ReadRestart(IO::DiscretizationReader& ioreader);

    //! derived
    virtual void RecoverState(
            const Epetra_Vector& xold,
            const Epetra_Vector& dir,
            const Epetra_Vector& xnew);

    //! derived
    virtual void UpdateStepState();

    //! derived
    virtual void UpdateStepElement();

    //! derived
    void DetermineStressStrain();

    //! derived
    void DetermineEnergy();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;
    //! @}

  protected:
    //! derived
    virtual void Reset(const Epetra_Vector& x);
    //! derived
    virtual void Reset(const Epetra_Vector& x, LINALG::SparseOperator& jac);

  private:

    bool ApplyForceInternal();

    bool ApplyForceExternal();

    bool ApplyForceStiffInternal();

    bool ApplyForceStiffExternal();

    void EvaluateNeumann(
        Teuchos::RCP<Epetra_Vector> eval_vec,
        Teuchos::RCP<LINALG::SparseOperator> eval_mat);

    void EvaluateNeumann(
        Teuchos::ParameterList& p,
        Teuchos::RCP<Epetra_Vector> eval_vec,
        Teuchos::RCP<LINALG::SparseOperator> eval_mat);

    void EvaluateInternal(
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    void EvaluateInternal(
        Teuchos::ParameterList& p,
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

  private:
    //! @name pointers to the data container content
    //! @{
    //! global internal force at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fintnp_ptr_;

    //! global external force at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fextnp_ptr_;

    //! inertial force at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> finertialnp_ptr_;

    //! viscous force at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> fvisconp_ptr_;

    //! structural displacement at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> disnp_ptr_;

    //! structural stiffness matrix
    Teuchos::RCP<LINALG::SparseMatrix> stiff_ptr_;

    //! structural mass matrix
    Teuchos::RCP<LINALG::SparseOperator> mass_ptr_;

    //! structural damping matrix
    Teuchos::RCP<LINALG::SparseOperator> damp_ptr_;

    //! structural element evaluation time
    double* dt_ele_ptr_;

    //! mass linearization type
    enum INPAR::STR::MassLin masslin_type_;
    //! @}

    //! @name class only variables @{
    /*! \brief displacement increment
     *  Necessary for the EAS reconstruction, incremental strain evaluation,
     *  etc.. */
    Teuchos::RCP<Epetra_Vector> dis_incr_ptr_;
    //! @}
  };

} // namespace ModelEvaluator
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_STRUCTURE_H_ */
