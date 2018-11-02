/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_partitioned_ssi_adhesiondynamics.H

\brief Model evaluator for partitioned ssi adhesion dynamics.

\maintainer Andreas Rauch

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SSI_STR_MODEL_EVALUATOR_PARTITIONED_SSI_ADHESIONDYNAMICS_H_
#define SSI_STR_MODEL_EVALUATOR_PARTITIONED_SSI_ADHESIONDYNAMICS_H_

#include "../drt_ssi/ssi_str_model_evaluator_partitioned.H"

// forward declarations
namespace ADAPTER
{
  class Structure;
}  // namespace ADAPTER

namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG

namespace STR
{
  namespace MODELEVALUATOR
  {
    class PartitionedSSIAdhesionDynamics : public PartitionedSSI
    {
     public:
      //! constructor
      PartitionedSSIAdhesionDynamics();

      //! destructor
      virtual ~PartitionedSSIAdhesionDynamics(){};

      //! setup class variables [derived]
      virtual void Setup();

      //! @name Functions which are derived from the base generic class
      //! @{
      //! [derived]
      INPAR::STR::ModelType Type() const { return INPAR::STR::model_partitioned_coupling; }

      //! derived
      virtual void Reset(const Epetra_Vector& x);

      //! [derived]
      virtual bool EvaluateForce();

      //! [derived]
      virtual bool EvaluateStiff();

      //! [derived]
      virtual bool EvaluateForceStiff();

      //! derived
      virtual bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const;

      //! Assemble the jacobian at \f$t_{n+1}\f$
      virtual bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

      //! [derived]
      virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;
      //! @}

      //! Returns a pointer to the adhesion forces
      virtual Teuchos::RCP<Epetra_Vector> GetAdhesionForcesPtr() { return adhesion_force_ptr_; };

      //! Sets a pointer to the adhesion fixpoints
      void SetAdhesionFixpointPtr(Teuchos::RCP<Epetra_Vector> me)
      {
        fixpoint_coord_ptr_ = me;
        return;
      }

     private:
      //! true if adhesion dynamics are turned on
      bool adhesion_dynamics_;

      //! adhesion penalty parameter
      double penalty_;

      //! number of different bound species
      int num_bound_species_;

      //! adhesion forces on cell
      Teuchos::RCP<Epetra_Vector> adhesion_force_ptr_;

      //! structural stiffness matrix
      Teuchos::RCP<LINALG::SparseMatrix> stiff_adhesion_force_ptr_;

      //! structural displacement at \f$t_{n+1}\f$
      Teuchos::RCP<const Epetra_Vector> disnp_ptr_;

      //! structural coordinates of adhesion fixpoints
      Teuchos::RCP<const Epetra_Vector> fixpoint_coord_ptr_;

    };  // class PartitionedSSIAdhesionDynamics

  }  // namespace MODELEVALUATOR
}  // namespace STR


#endif /* SSI_STR_MODEL_EVALUATOR_PARTITIONED_SSI_ADHESIONDYNAMICS_H_ */
