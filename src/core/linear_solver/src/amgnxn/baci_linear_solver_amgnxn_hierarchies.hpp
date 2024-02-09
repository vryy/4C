/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_AMGNXN_HPPIERARCHIES_HPP
#define BACI_LINEAR_SOLVER_AMGNXN_HPPIERARCHIES_HPP

#include "baci_config.hpp"

#include "baci_linalg_blocksparsematrix.hpp"
#include "baci_linear_solver_amgnxn_smoothers.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_linear_solver_preconditioner_type.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_Utilities.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER::AMGNXN
{
  class Hierarchies
  {
   public:
    Hierarchies(Teuchos::RCP<AMGNXN::BlockedMatrix> A,
        std::vector<Teuchos::ParameterList> muelu_params, std::vector<int> num_pdes,
        std::vector<int> null_spaces_dim,
        std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data, int NumLevelAMG,
        std::string verbosity = "off");

    int GetNumLevelMin();
    int GetNumBlocks();
    int GetNumLevels(int block);

    Teuchos::RCP<AMGNXN::BlockedMatrix> GetBlockMatrix();
    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> GetH(int block);

    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetA(int block, int level);
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetP(int block, int level);
    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetR(int block, int level);
    Teuchos::RCP<AMGNXN::MueluSmootherWrapper> GetSPre(int block, int level);
    Teuchos::RCP<AMGNXN::MueluSmootherWrapper> GetSPos(int block, int level);

    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> GetA(int block);
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> GetP(int block);
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> GetR(int block);
    std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>> GetSPre(int block);
    std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>> GetSPos(int block);

    int GetNumPDEs(int block);
    int GetNullSpaceDim(int block);
    Teuchos::RCP<std::vector<double>> GetNullSpaceData(int block);

   private:
    Teuchos::RCP<AMGNXN::BlockedMatrix> A_;
    std::vector<Teuchos::ParameterList> muelu_params_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    int NumBlocks_;
    int NumLevelMax_;
    int NumLevelMin_;
    int NumLevelAMG_;
    std::vector<Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> H_block_;
    std::vector<std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>>> A_block_level_;
    std::vector<std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>>> P_block_level_;
    std::vector<std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>>> R_block_level_;
    std::vector<std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>>> SPre_block_level_;
    std::vector<std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>>> SPos_block_level_;
    std::string verbosity_;

    void Setup();

    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BuildMueLuHierarchy(
        Teuchos::ParameterList paramListFromXml, int numdf, int dimns,
        Teuchos::RCP<std::vector<double>> nsdata, Teuchos::RCP<Epetra_Operator> A_eop, int block,
        int NumBlocks, std::vector<int>& offsets, int offsetFineLevel);

    std::string ConvertInt(int number)
    {
      std::stringstream ss;
      ss << number;
      return ss.str();
    }
  };

  class Vcycle;

  class MonolithicHierarchy
  {
   public:
    MonolithicHierarchy(Teuchos::RCP<AMGNXN::Hierarchies> H, const Teuchos::ParameterList& params,
        const Teuchos::ParameterList& params_smoothers);

    int GetNumLevels();

    Teuchos::RCP<AMGNXN::Hierarchies> GetHierarchies();

    Teuchos::RCP<AMGNXN::BlockedMatrix> GetA(int level);

    Teuchos::RCP<Vcycle> BuildVCycle();

   private:
    Teuchos::RCP<AMGNXN::Hierarchies> H_;
    int NumLevels_;
    int NumBlocks_;
    std::vector<Teuchos::RCP<AMGNXN::BlockedMatrix>> A_;
    std::vector<Teuchos::RCP<AMGNXN::BlockedMatrix>> P_;
    std::vector<Teuchos::RCP<AMGNXN::BlockedMatrix>> R_;
    std::vector<Teuchos::RCP<AMGNXN::GenericSmoother>> Spre_;
    std::vector<Teuchos::RCP<AMGNXN::GenericSmoother>> Spos_;
    Teuchos::ParameterList params_;
    Teuchos::ParameterList params_smoothers_;

    void Setup();
    Teuchos::RCP<AMGNXN::GenericSmoother> BuildSmoother(int level);
  };
}  // namespace CORE::LINEAR_SOLVER::AMGNXN

BACI_NAMESPACE_CLOSE

#endif  // SOLVER_AMGNXN_HIERARCHIES_H
