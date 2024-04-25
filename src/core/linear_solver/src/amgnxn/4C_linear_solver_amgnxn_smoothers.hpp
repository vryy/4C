/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 1
Created on: Feb 27, 2014
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_AMGNXN_SMOOTHERS_HPP
#define FOUR_C_LINEAR_SOLVER_AMGNXN_SMOOTHERS_HPP

// Trilinos includes
#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linear_solver_amgnxn_objects.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Ifpack.h>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_Utilities.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace CORE::LINEAR_SOLVER::AMGNXN
{
  class GenericSmoother
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~GenericSmoother() = default;

    virtual void Solve(
        const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const = 0;

    void Richardson(Teuchos::RCP<GenericSmoother> Ainv, const BlockedMatrix& A,
        const BlockedVector& X, BlockedVector& Y, int iters, double omega,
        bool InitialGuessIsZero) const;
    // if InitialGuessIsZero == true we can input any random initial guess and the smoother will
    // take care of making the final result be as if the initial guess would be zero. This
    // avoids to scale to zero the initial guess, and make a little more efficient the smoother
  };

  class SingleFieldSmoother : public GenericSmoother
  {
   public:
    void Solve(
        const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const override
    {
      CheckSingleFieldVector(X);
      CheckSingleFieldVector(Y);
      Apply(*(X.GetVector(0)), *(Y.GetVector(0)), InitialGuessIsZero);
      return;
    }

    virtual void Apply(
        const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const = 0;

   protected:
    void CheckSingleFieldVector(const BlockedVector& V) const
    {
      if (not V.HasOnlyOneBlock()) FOUR_C_THROW("We need here a single field vector");
      return;
    }
  };

  class BlockedSmoother : public GenericSmoother
  {
  };

  class BgsSmoother : public BlockedSmoother
  {
   public:
    BgsSmoother(Teuchos::RCP<BlockedMatrix> A, std::vector<Teuchos::RCP<GenericSmoother>> smoothers,
        std::vector<std::vector<int>> superblocks, unsigned iter, double omega,
        std::vector<unsigned> iters, std::vector<double> omegas)
        : a_(A),
          smoothers_(smoothers),
          superblocks_(superblocks),
          iter_(iter),
          omega_(omega),
          iters_(iters),
          omegas_(omegas)
    {
    }

    void Solve(
        const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const override;

   private:
    Teuchos::RCP<BlockedMatrix> a_;
    std::vector<Teuchos::RCP<GenericSmoother>> smoothers_;
    std::vector<std::vector<int>> superblocks_;
    unsigned iter_;
    double omega_;
    std::vector<unsigned> iters_;
    std::vector<double> omegas_;
  };

  class SimpleSmoother : public BlockedSmoother
  {
   public:
    SimpleSmoother(Teuchos::RCP<BlockedMatrix> A, Teuchos::RCP<BlockedMatrix> invApp,
        Teuchos::RCP<BlockedMatrix> Schur, Teuchos::RCP<GenericSmoother> SmooApp,
        Teuchos::RCP<GenericSmoother> SmooSchur, std::vector<int> BlocksPred,
        std::vector<int> BlocksSchur, unsigned iter, double alpha)
        : a_(A),
          inv_app_(invApp),
          schur_(Schur),
          smoo_app_(SmooApp),
          smoo_schur_(SmooSchur),
          blocks_pred_(BlocksPred),
          blocks_schur_(BlocksSchur),
          iter_(iter),
          alpha_(alpha)
    {
    }

    void Solve(
        const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const override;

   private:
    Teuchos::RCP<BlockedMatrix> a_;
    Teuchos::RCP<BlockedMatrix> inv_app_;
    Teuchos::RCP<BlockedMatrix> schur_;
    Teuchos::RCP<GenericSmoother> smoo_app_;
    Teuchos::RCP<GenericSmoother> smoo_schur_;
    std::vector<int> blocks_pred_;
    std::vector<int> blocks_schur_;
    unsigned iter_;
    double alpha_;
    mutable Teuchos::RCP<BlockedVector> xp_tmp_;
    mutable Teuchos::RCP<BlockedVector> xs_tmp_;
    mutable Teuchos::RCP<BlockedVector> yp_tmp_;
    mutable Teuchos::RCP<BlockedVector> d_ys_;
    mutable Teuchos::RCP<BlockedVector> d_xp_;
    mutable Teuchos::RCP<BlockedVector> d_xs_;
  };

  class MergeAndSolve : public BlockedSmoother
  {
   public:
    MergeAndSolve()
        : solver_(Teuchos::null),
          sparse_matrix_(Teuchos::null),
          block_sparse_matrix_(Teuchos::null),
          a_(Teuchos::null),
          x_(Teuchos::null),
          b_(Teuchos::null),
          is_set_up_(false)
    {
    }

    void Setup(BlockedMatrix matrix);

    void Solve(
        const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const override;

   private:
    Teuchos::RCP<CORE::LINALG::Solver> solver_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sparse_matrix_;
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> block_sparse_matrix_;
    Teuchos::RCP<Epetra_Operator> a_;
    mutable Teuchos::RCP<Epetra_MultiVector> x_;
    mutable Teuchos::RCP<Epetra_MultiVector> b_;
    bool is_set_up_;
  };

  // Forward declarations
  class Hierarchies;
  class MonolithicHierarchy;
  class Vcycle;
  class VcycleSingle;

  class CoupledAmg : public BlockedSmoother
  {
   public:
    CoupledAmg(Teuchos::RCP<AMGNXN::BlockedMatrix> A, std::vector<int> num_pdes,
        std::vector<int> null_spaces_dim,
        std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data,
        const Teuchos::ParameterList& amgnxn_params, const Teuchos::ParameterList& smoothers_params,
        const Teuchos::ParameterList& muelu_params);

    void Solve(
        const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const override;

   private:
    void Setup();

    Teuchos::RCP<AMGNXN::BlockedMatrix> a_;
    std::vector<Teuchos::ParameterList> muelu_lists_;
    std::vector<int> num_pdes_;
    std::vector<int> null_spaces_dim_;
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
    Teuchos::ParameterList amgnxn_params_;
    Teuchos::ParameterList smoothers_params_;
    Teuchos::ParameterList muelu_params_;

    bool is_setup_flag_;
    Teuchos::RCP<AMGNXN::Hierarchies> h_;
    Teuchos::RCP<AMGNXN::MonolithicHierarchy> m_;
    Teuchos::RCP<AMGNXN::Vcycle> v_;
  };

  class MueluSmootherWrapper : public SingleFieldSmoother
  {
   public:
    MueluSmootherWrapper(
        Teuchos::RCP<MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>> S)
        : s_(S)
    {
    }

    void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
        bool InitialGuessIsZero = false) const override;

   private:
    Teuchos::RCP<MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>> s_;
  };

  class MueluHierarchyWrapper : public SingleFieldSmoother  // Not used
  {
   public:
    MueluHierarchyWrapper(
        Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H);

    void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
        bool InitialGuessIsZero = false) const override;

   private:
    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> h_;
    Teuchos::RCP<Epetra_Operator> p_;
  };

  class MueluAMGWrapper : public SingleFieldSmoother
  {
   public:
    MueluAMGWrapper(Teuchos::RCP<CORE::LINALG::SparseMatrix> A, int num_pde, int null_space_dim,
        Teuchos::RCP<std::vector<double>> null_space_data,
        const Teuchos::ParameterList& muelu_list);

    void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
        bool InitialGuessIsZero = false) const override;

    void Setup();

   protected:
    Teuchos::RCP<CORE::LINALG::SparseMatrix> A_;
    int num_pde_;
    int null_space_dim_;
    Teuchos::RCP<std::vector<double>> null_space_data_;
    Teuchos::ParameterList muelu_list_;
    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;
    void BuildHierarchy();

   private:
    Teuchos::RCP<Epetra_Operator> p_;
  };

  class SingleFieldAMG : public MueluAMGWrapper
  {
   public:
    SingleFieldAMG(Teuchos::RCP<CORE::LINALG::SparseMatrix> A, int num_pde, int null_space_dim,
        Teuchos::RCP<std::vector<double>> null_space_data, const Teuchos::ParameterList& muelu_list,
        const Teuchos::ParameterList& fine_smoother_list);

    void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
        bool InitialGuessIsZero = false) const override;

   private:
    Teuchos::ParameterList fine_smoother_list_;
    Teuchos::RCP<VcycleSingle> v_;
    void Setup();
  };

  class IfpackWrapper : public SingleFieldSmoother
  {
   public:
    IfpackWrapper(Teuchos::RCP<CORE::LINALG::SparseMatrixBase> A, Teuchos::ParameterList& list);
    ~IfpackWrapper() override { delete prec_; }
    void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
        bool InitialGuessIsZero = false) const override;

   private:
    Ifpack_Preconditioner* prec_;
    Teuchos::RCP<CORE::LINALG::SparseMatrixBase> a_;
    Teuchos::RCP<Epetra_RowMatrix> arow_;
    Teuchos::ParameterList list_;
    std::string type_;
  };

  class DirectSolverWrapper : public SingleFieldSmoother
  {
   public:
    DirectSolverWrapper();
    void Setup(Teuchos::RCP<CORE::LINALG::SparseMatrix> matrix,
        Teuchos::RCP<Teuchos::ParameterList> params);

    void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
        bool InitialGuessIsZero = false) const override;

   private:
    Teuchos::RCP<CORE::LINALG::Solver> solver_;
    Teuchos::RCP<Epetra_Operator> a_;
    mutable Teuchos::RCP<Epetra_MultiVector> x_;
    mutable Teuchos::RCP<Epetra_MultiVector> b_;
    bool is_set_up_;
  };

  // Auxiliary class to wrap the null space data to be used within the smoothers
  class NullSpaceInfo
  {
   public:
    NullSpaceInfo() {}
    NullSpaceInfo(
        int num_pdes, int null_space_dim, Teuchos::RCP<std::vector<double>> null_space_data)
        : num_pdes_(num_pdes), null_space_dim_(null_space_dim), null_space_data_(null_space_data)
    {
    }

    int GetNumPDEs() { return num_pdes_; }
    int GetNullSpaceDim() { return null_space_dim_; }
    Teuchos::RCP<std::vector<double>> GetNullSpaceData() { return null_space_data_; }

   private:
    int num_pdes_;
    int null_space_dim_;
    Teuchos::RCP<std::vector<double>> null_space_data_;
  };

  class Hierarchies;  // forward declaration

  class SmootherManager  // TODO: this is quite lengthy. This can be done with a ParameterList
  {
   public:
    SmootherManager();
    Teuchos::RCP<BlockedMatrix> GetOperator();
    Teuchos::ParameterList GetParams();
    Teuchos::ParameterList GetParamsSmoother();
    Teuchos::RCP<Hierarchies> GetHierarchies();
    int GetLevel();
    int GetBlock();
    std::vector<int> GetBlocks();
    std::string GetSmootherName();
    std::string GetType();
    std::string GetVerbosity();
    NullSpaceInfo GetNullSpace();
    std::vector<NullSpaceInfo> GetNullSpaceAllBlocks();

    void SetOperator(Teuchos::RCP<BlockedMatrix> in);
    void SetParams(const Teuchos::ParameterList& in);
    void SetParamsSmoother(const Teuchos::ParameterList& in);
    void SetHierarchies(Teuchos::RCP<Hierarchies> in);
    void SetLevel(int in);
    void SetBlock(int in);
    void SetBlocks(std::vector<int> in);
    void SetSmootherName(std::string in);
    void SetType(std::string in);
    void SetVerbosity(std::string in);
    void SetNullSpace(const NullSpaceInfo& in);
    void SetNullSpaceAllBlocks(const std::vector<NullSpaceInfo>& in);

    bool IsSetOperator();
    bool IsSetParams();
    bool IsSetParamsSmoother();
    bool IsSetHierarchies();
    bool IsSetLevel();
    bool IsSetBlock();
    bool IsSetBlocks();
    bool IsSetSmootherName();
    bool IsSetType();
    bool IsSetVerbosity();
    bool IsSetNullSpace();
    bool IsSetNullSpaceAllBlocks();

   private:
    Teuchos::RCP<BlockedMatrix> operator_;
    Teuchos::ParameterList params_;
    Teuchos::ParameterList params_subsolver_;
    Teuchos::RCP<Hierarchies> hierarchies_;
    int level_;
    int block_;
    std::vector<int> blocks_;
    std::string subsolver_name_;
    std::string type_;
    std::string verbosity_;
    NullSpaceInfo null_space_;
    std::vector<NullSpaceInfo> null_space_all_blocks_;

    bool set_operator_;
    bool set_params_;
    bool set_params_subsolver_;
    bool set_hierarchies_;
    bool set_level_;
    bool set_block_;
    bool set_blocks_;
    bool set_subsolver_name_;
    bool set_type_;
    bool set_verbosity_;
    bool set_null_space_;
    bool set_null_space_all_blocks_;
  };

  class SmootherFactoryBase : public SmootherManager
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~SmootherFactoryBase() = default;

    virtual Teuchos::RCP<GenericSmoother> Create() = 0;
  };

  // This class is able to create any smoother. The smoother to be created is given in a
  // parameter list
  class SmootherFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;

   private:
    void SetTypeAndParams();
  };

  class BgsSmootherFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;

   private:
    void ParseSmootherNames(const std::string& smoothers_string,
        std::vector<std::string>& smoothers_vector, std::vector<std::vector<int>> superblocks);
  };

  class CoupledAmgFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };

  class SimpleSmootherFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;

   private:
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ApproximateInverse(
        const CORE::LINALG::SparseMatrixBase& A, const std::string& method);
    Teuchos::RCP<BlockedMatrix> ComputeSchurComplement(const BlockedMatrix& invApp,
        const BlockedMatrix& Aps, const BlockedMatrix& Asp, const BlockedMatrix& Ass);
  };

  class MergeAndSolveFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };

  class IfpackWrapperFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };

  class MueluSmootherWrapperFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };

  class HierarchyRemainderWrapperFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };

  class MueluAMGWrapperFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };

  class SingleFieldAMGFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };

  class DirectSolverWrapperFactory : public SmootherFactoryBase
  {
   public:
    Teuchos::RCP<GenericSmoother> Create() override;
  };
}  // namespace CORE::LINEAR_SOLVER::AMGNXN

FOUR_C_NAMESPACE_CLOSE

#endif
