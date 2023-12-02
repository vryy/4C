/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_AMGNXN_PRECONDITIONER_H
#define BACI_LINEAR_SOLVER_AMGNXN_PRECONDITIONER_H

#include "baci_linalg_blocksparsematrix.H"
#include "baci_linear_solver_amgnxn_hierarchies.H"
#include "baci_linear_solver_amgnxn_smoothers.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_linear_solver_preconditioner_type.H"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_Utilities.hpp>
#include <Teuchos_RCP.hpp>


namespace CORE::LINEAR_SOLVER
{
  class AMGnxn_Preconditioner : public PreconditionerType
  {
   public:
    AMGnxn_Preconditioner(Teuchos::ParameterList &params);

    void Setup(bool create, Epetra_Operator *matrix, Epetra_MultiVector *x,
        Epetra_MultiVector *b) override;

    virtual void Setup(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A);

    Epetra_Operator *PrecOperator() const override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const override;

    void Print(std::ostream &stream) const override { stream << "AMGnxn"; }
    std::string getParameterListName() const override { return "AMGnxn Parameters"; }

   private:
    // Private variables
    Teuchos::RCP<Epetra_Operator> P_;                      // The underlying preconditioner object
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A_;  // A own copy of the system matrix
    Teuchos::ParameterList &params_;

  };  // AMGnxn_Preconditioner

  class AMGnxn_Interface
  {
   public:
    AMGnxn_Interface(Teuchos::ParameterList &params, int NumBlocks);

    std::vector<std::string> GetMueLuXmlFiles() { return xml_files_; }
    std::vector<int> GetNumPdes() { return num_pdes_; }
    std::vector<int> GetNullSpacesDim() { return null_spaces_dim_; }
    std::vector<Teuchos::RCP<std::vector<double>>> GetNullSpacesData() { return null_spaces_data_; }
    // int GetNumLevelAMG(){return NumLevelAMG_;}
    Teuchos::ParameterList GetPreconditionerParams() { return prec_params_; }
    Teuchos::ParameterList GetSmoothersParams() { return smoo_params_; }
    std::string GetPreconditionerType() { return prec_type_; }

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
    std::string ConvertInt(int number)
    {
      std::stringstream ss;
      ss << number;
      return ss.str();
    }
  };

  class AMGnxn_Operator : virtual public Epetra_Operator
  {
   public:
    AMGnxn_Operator(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A, std::vector<int> num_pdes,
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
      dserror("Function not implemented");
      return -1;
    }

    double NormInf() const override
    {
      dserror("Function not implemented");
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
      dserror("Function not implemented");
      return false;
    }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Comm &Comm() const override { return A_->Comm(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorDomainMap() const override { return A_->OperatorDomainMap(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorRangeMap() const override { return A_->OperatorRangeMap(); }

    void Setup();

   private:
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A_;
    std::vector<Teuchos::ParameterList> muelu_lists_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    Teuchos::ParameterList amgnxn_params_;
    Teuchos::ParameterList smoothers_params_;
    Teuchos::ParameterList muelu_params_;

    bool is_setup_flag_;

    Teuchos::RCP<AMGNXN::CoupledAmg> V_;

  };  // class AMGnxn_Operator

  class BlockSmoother_Operator : virtual public Epetra_Operator
  {
   public:
    BlockSmoother_Operator(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A,
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
      dserror("Function not implemented");
      return -1;
    }

    double NormInf() const override
    {
      dserror("Function not implemented");
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
      dserror("Function not implemented");
      return false;
    }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Comm &Comm() const override { return A_->Comm(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorDomainMap() const override { return A_->OperatorDomainMap(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorRangeMap() const override { return A_->OperatorRangeMap(); }

    void Setup();  // TODO

   private:
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    Teuchos::ParameterList amgnxn_params_;
    Teuchos::ParameterList smoothers_params_;

    bool is_setup_flag_;
    Teuchos::RCP<AMGNXN::BlockedSmoother> S_;
    Teuchos::RCP<AMGNXN::GenericSmoother> Sbase_;

  };  // class BlockSmoother_Operator

  class Merged_Operator : virtual public Epetra_Operator
  {
   public:
    Merged_Operator(Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A,
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
      dserror("Function not implemented");
      return -1;
    }

    double NormInf() const override
    {
      dserror("Function not implemented");
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
      dserror("Function not implemented");
      return false;
    }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Comm &Comm() const override { return A_->Comm(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorDomainMap() const override { return A_->OperatorDomainMap(); }

    // Only required to properly define an Epetra_Operator, not should be used!
    const Epetra_Map &OperatorRangeMap() const override { return A_->OperatorRangeMap(); }

    void Setup();  // TODO

   private:
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> Asp_;
    Teuchos::ParameterList amgnxn_params_;
    Teuchos::ParameterList smoothers_params_;

    bool is_setup_flag_;
    Teuchos::RCP<AMGNXN::IfpackWrapper> S_;

  };  // class Merged_Operator

  void PrintMap(const Epetra_Map &Map, std::string prefix);

}  // namespace CORE::LINEAR_SOLVER

#endif  // SOLVER_AMGNXN_PRECONDITIONER_H
