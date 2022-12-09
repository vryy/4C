/*-----------------------------------------------------------*/
/*! \file
\brief Model evaluator for structure part of partitioned ssi.

\level 3


*/
/*-----------------------------------------------------------*/


#ifndef SSI_STR_MODEL_EVALUATOR_PARTITIONED_H_
#define SSI_STR_MODEL_EVALUATOR_PARTITIONED_H_

#include "ssi_str_model_evaluator_base.H"

namespace ADAPTER
{
  class Structure;
}  // namespace ADAPTER

namespace SSI
{
  class SSIPart;
}

namespace STR
{
  namespace MODELEVALUATOR
  {
    class PartitionedSSI : public BaseSSI
    {
     public:
      //! constructor
      PartitionedSSI(const Teuchos::RCP<const SSI::SSIPart>
              ssi_part  //!< partitioned algorithm for scalar-structure interaction
      );

      void Setup() override;

      //! @name Functions which are derived from the base generic class
      //! @{
      [[nodiscard]] INPAR::STR::ModelType Type() const override
      {
        return INPAR::STR::model_partitioned_coupling;
      }

      bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const override;

      bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const override;

      void DetermineStressStrain() override{};

      void RunPreComputeX(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
          const NOX::NLN::Group& curr_grp) override;
      //! @}

     private:
      //! partitioned algorithm for scalar-structure interaction
      const Teuchos::RCP<const SSI::SSIPart> ssi_part_;
    };  // class PartitionedSSI

  }  // namespace MODELEVALUATOR
}  // namespace STR


#endif /* SSI_STR_MODEL_EVALUATOR_PARTITIONED_H_ */
