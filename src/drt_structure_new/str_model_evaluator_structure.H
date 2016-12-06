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

#include "str_elements_paramsinterface.H"     // interface to the element evaluation

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
    INPAR::STR::ModelType Type() const
    { return INPAR::STR::model_structure; }

    //! derived
    virtual void Reset(const Epetra_Vector& x);

    //! derived
    virtual bool EvaluateForce();

    //! derived
    virtual bool EvaluateStiff();

    //! derived
    virtual bool EvaluateForceStiff();

    /*! \brief Initialize viscous and inertial matrices
     *
     *  This is the place where we calculate the default mass matrix and the
     *  Rayleigh damping matrix only once during the EquilibriateInitialState routine.
     *
     *  \date 09/16
     *  \author hiermeier */
    bool InitializeInertiaAndDamping();

    //! derived
    virtual bool AssembleForce(Epetra_Vector& f,
        const double & timefac_np) const;

    //! derived
    virtual bool AssembleJacobian(LINALG::SparseOperator& jac,
        const double & timefac_np) const;

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
    virtual void UpdateStepState(const double& timefac_n);

    //! derived
    virtual void UpdateStepElement();

    //! derived
    void DetermineStressStrain();

    //! derived
    void DetermineEnergy();

    //! derived
    virtual void ResetStepState();

    //! derived
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;

    //! [derived]
    virtual void PostOutput();

    //! @}

  protected:

    //! pre-operator for \ref EvaluateInternal
    virtual void PreEvaluateInternal(){/* empty */};

  private:
    //! apply the internal force contributions
    bool ApplyForceInternal();

    //! apply the external force contributions
    bool ApplyForceExternal();

    //! apply the internal force contributions and the evaluate the structural stiffness terms
    bool ApplyForceStiffInternal();

    //! apply the external force contributions and evaluate possible linearization contributions
    bool ApplyForceStiffExternal();

    //! Set the ParamsInterface in the parameter list and call the other EvaluateNeumann routine
    void EvaluateNeumann(
        Teuchos::RCP<Epetra_Vector> eval_vec,
        Teuchos::RCP<LINALG::SparseOperator> eval_mat);

    /*! \brief Check if the given parameter list is valid and call the
     *  EvaluateNeumann routine of the discretization
     *
     *  \param eval_vec (out) : external force vector
     *  \param eval_mat (out) : linearization of the external force (optional)
     *
     *  \date 08/15
     *  \author hiermeier */
    void EvaluateNeumann(
        Teuchos::ParameterList& p,
        Teuchos::RCP<Epetra_Vector> eval_vec,
        Teuchos::RCP<LINALG::SparseOperator> eval_mat);

    //! Set the ParamsInterface in the parameter list and call the other EvaluateInternal routine
    void EvaluateInternal(
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    /*! \brief Check if the given parameter list is valid and call the
     *  Evaluate routine of the discretization
     *
     *  \param eval_vec (out) : array of different internal forces (f_int, f_inertial)
     *  \param eval_mat (out) : array of different matrices (stiffness, mass, damping)
     *
     *  \date 08/15
     *  \author hiermeier */
    void EvaluateInternal(
        Teuchos::ParameterList& p,
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    /*! \brief Add static structural internal force and stiffness matrix to the
     *         evaluate call (default)
     *
     *  Set matrix, vector and default action type.
     *
     *  \param eval_mat (out): pointer to the evaluation matrix array, which is
     *                         changed.
     *  \param eval_vec (out): pointer to the evaluation vector array, which is
     *                         changed.
     *
     *  \date 09/16
     *  \author hiermeier */
    void StaticContributions(
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    /*! \brief Add static structural internal force to the evaluate call (default)
     *
     *  Set vector and default action type.
     *
     *  \param eval_vec (out): pointer to the evaluation vector array, which is
     *                         changed.
     *
     *  \date 09/16
     *  \author hiermeier */
    void StaticContributions(
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    /*! \brief Add material damping matrix  to the evaluate call (optional)
     *
     *  \param eval_mat (out): pointer to the evaluation matrix array,
     *                         which will be augmented with the damping matrix
     *                         if desired.
     *
     *  \warning Material damping and non-linear mass effects cannot be
     *  considered at the same time at the moment!
     *
     *  \date 09/16
     *  \author hiermeier */
    void MaterialDampingContributions(
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat);

    /*! \brief Add mass matrix and inertial force to the evaluate call (optional)
     *
     *  Set vector, matrix and default mass or lumped mass action type.
     *
     *  \param eval_mat (out): pointer to the evaluation matrix array,
     *                         which will be augmented with the mass matrix
     *                         if desired.
     *
     *  \param eval_vec (out): pointer to the evaluation vector array,
     *                         which will be augmented with the inertial force
     *                         vector if desired.
     *
     *  \warning Material damping and non-linear mass effects cannot be
     *  considered at the same time at the moment!
     *
     *  \date 09/16
     *  \author hiermeier */
    void InertialContributions(
        Teuchos::RCP<LINALG::SparseOperator>* eval_mat,
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    /*! \brief Add inertial force to the evaluate call (optional)
     *
     *  Set vectors and internal inertial force action type.
     *
     *  \param eval_vec (out): pointer to the evaluation vector array,
     *                         which will be augmented with the inertial force
     *                         vector if desired.
     *
     *  \warning Material damping and non-linear mass effects cannot be
     *  considered at the same time at the moment!
     *
     *  \date 09/16
     *  \author hiermeier */
    void InertialContributions(
        Teuchos::RCP<Epetra_Vector>* eval_vec);

    /*! \brief Evaluate the inertial forces (for the standard case) and
     *         any viscous damping forces
     *
     *  \date 09/16
     *  \author hiermeier */
    void InertialAndViscousForces();

    /*! Check the FillComplete status of the stiffness and mass matrix
     *  and complete them, if necessary */
    void FillComplete();

    /*! \biref Assemble the Rayleigh damping matrix
     *
     *  Please note, that this has to been done only once during the
     *  STR::Integrator::EquilibriateInitialState routine!
     *
     *  \date 09/16
     *  \author hiermeier */
    void RayleighDampingMatrix();

    /*! \brief Returns the interial force vector for non-linear mass problems
     *
     *  This function zeros the inertial force vector and returns it,
     *  if a non-linear mass problem is solved. Otherwise, a Teuchos::null
     *  pointer is returned.
     *
     *  \date 09/16
     *  \author hiermeier */
    Teuchos::RCP<Epetra_Vector> GetInertialForce();

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
