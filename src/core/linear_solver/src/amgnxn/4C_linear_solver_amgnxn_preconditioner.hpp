/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_AMGNXN_PRECONDITIONER_HPP
#define FOUR_C_LINEAR_SOLVER_AMGNXN_PRECONDITIONER_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linear_solver_amgnxn_hierarchies.hpp"
#include "4C_linear_solver_amgnxn_smoothers.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_Utilities.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace CORE::LINEAR_SOLVER
{
  class AmGnxnPreconditioner : public PreconditionerTypeBase
  {
   public:
    AmGnxnPreconditioner(Teuchos::ParameterList &params);

    void Setup(bool create, Epetra_Operator *matrix, Epetra_MultiVector *x,
        Epetra_MultiVector *b) override;

    virtual void Setup(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A);

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override;

   private:
    // Private variables
    Teuchos::RCP<Epetra_Operator> p_;                      // The underlying preconditioner object
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a_;  // A own copy of the system matrix
    Teuchos::ParameterList &params_;

  };  // AMGnxn_Preconditioner

  class AmGnxnInterface
  {
   public:
    AmGnxnInterface(Teuchos::ParameterList &params, int NumBlocks);

    std::vector<std::string> GetMueLuXmlFiles() { return xml_files_; }
    std::vector<int> GetNumPdes() { return num_pdes_; }
    std::vector<int> GetNullSpacesDim() { return null_spaces_dim_; }
    std::vector<Teuchos::RCP<std::vector<double>>> GetNullSpacesData() { return null_spaces_data_; }
    // int GetNumLevelAMG(){return NumLevelAMG_;}
    Teuchos::ParameterList get_preconditioner_params() { return prec_params_; }
    Teuchos::ParameterList GetSmoothersParams() { return smoo_params_; }
    std::string get_preconditioner_type() { return prec_type_; }

   private:
    std::vector<std::string> xml_files_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    // int NumLevelAMG_;
    Teuchos::ParameterList prec_params_;
    Teuchos::ParameterList smoo_params_;
    std::string prec_type_;

    void Params_TSI_AMG_BGS(Teuchos::ParameterList &params);

    // Helper function to convert int to std::string
    std::string convert_int(int number)
    {
      std::stringstream ss;
      ss << number;
      return ss.str();
    }
  };

  class AmGnxnOperator : virtual public Epetra_Operator
  {
   public:
    AmGnxnOperator(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A, std::vector<int> num_pdes,
        std::vector<int> null_spaces_dim,
        std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data,
        const Teuchos::ParameterList &amgnxn_params, const Teuchos::ParameterList &smoothers_params,
        const Teuchos::ParameterList &muelu_params);

    // virtual functions given by Epetra_Operator. The only one to be used is ApplyInverse()
    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override;

    int SetUseTranspose(bool UseTranspose) override
    {
      // default to false
      return 0;
    }

    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override
    {
      FOUR_C_THROW("Function not implemented");
      return -1;
    }

    double NormInf() const override
    {
      FOUR_C_THROW("Function not implemented");
      return -1.0;
    }

    const char *Label() const override { return "AMG(BlockSmoother)"; }

    bool UseTranspose() const override
    {
      // default to false
      return false;
    }

    bool HasNormInf() const override
    {
      FOUR_C_THROW("Function not implemented");
      return false;
    }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Comm &Comm() const override { return a_->Comm(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorDomainMap() const override { return a_->OperatorDomainMap(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorRangeMap() const override { return a_->OperatorRangeMap(); }

    void Setup();

   private:
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a_;
    std::vector<Teuchos::ParameterList> muelu_lists_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    Teuchos::ParameterList amgnxn_params_;
    Teuchos::ParameterList smoothers_params_;
    Teuchos::ParameterList muelu_params_;

    bool is_setup_flag_;

    Teuchos::RCP<AMGNXN::CoupledAmg> v_;

  };  // class AMGnxn_Operator

  class BlockSmootherOperator : virtual public Epetra_Operator
  {
   public:
    BlockSmootherOperator(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A,
        std::vector<int> num_pdes, std::vector<int> null_spaces_dim,
        std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data,
        const Teuchos::ParameterList &amgnxn_params,
        const Teuchos::ParameterList &smoothers_params);

    // virtual functions given by Epetra_Operator. The only one to be used is ApplyInverse()
    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override;

    int SetUseTranspose(bool UseTranspose) override
    {
      // default to false
      return 0;
    }

    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override
    {
      FOUR_C_THROW("Function not implemented");
      return -1;
    }

    double NormInf() const override
    {
      FOUR_C_THROW("Function not implemented");
      return -1.0;
    }

    const char *Label() const override { return "BlockSmoother(X)"; }

    bool UseTranspose() const override
    {
      // default to false
      return false;
    }

    bool HasNormInf() const override
    {
      FOUR_C_THROW("Function not implemented");
      return false;
    }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Comm &Comm() const override { return a_->Comm(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorDomainMap() const override { return a_->OperatorDomainMap(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorRangeMap() const override { return a_->OperatorRangeMap(); }

    void Setup();  // TODO

   private:
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    Teuchos::ParameterList amgnxn_params_;
    Teuchos::ParameterList smoothers_params_;

    bool is_setup_flag_;
    Teuchos::RCP<AMGNXN::BlockedSmoother> s_;
    Teuchos::RCP<AMGNXN::GenericSmoother> sbase_;

  };  // class BlockSmoother_Operator

  class MergedOperator : virtual public Epetra_Operator
  {
   public:
    MergedOperator(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A,
        const Teuchos::ParameterList &amgnxn_params,
        const Teuchos::ParameterList &smoothers_params);

    // virtual functions given by Epetra_Operator. The only one to be used is ApplyInverse()
    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override;

    int SetUseTranspose(bool UseTranspose) override
    {
      // default to false
      return 0;
    }

    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const override
    {
      FOUR_C_THROW("Function not implemented");
      return -1;
    }

    double NormInf() const override
    {
      FOUR_C_THROW("Function not implemented");
      return -1.0;
    }

    const char *Label() const override { return "Merged matrix plus smoother"; }

    bool UseTranspose() const override
    {
      // default to false
      return false;
    }

    bool HasNormInf() const override
    {
      FOUR_C_THROW("Function not implemented");
      return false;
    }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Comm &Comm() const override { return a_->Comm(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorDomainMap() const override { return a_->OperatorDomainMap(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorRangeMap() const override { return a_->OperatorRangeMap(); }

    void Setup();  // TODO

   private:
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> asp_;
    Teuchos::ParameterList amgnxn_params_;
    Teuchos::ParameterList smoothers_params_;

    bool is_setup_flag_;
    Teuchos::RCP<AMGNXN::IfpackWrapper> s_;

  };  // class Merged_Operator

  void PrintMap(const Epetra_Map &Map, std::string prefix);

}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
