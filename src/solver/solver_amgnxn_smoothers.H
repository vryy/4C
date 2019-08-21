/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
*/
/*----------------------------------------------------------------------*/
#ifndef AMGNXN_SMOOTHERS_H
#define AMGNXN_SMOOTHERS_H
#ifdef HAVE_MueLu

// Trilinos includes
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_UseDefaultTypes.hpp>

// Baci includes
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "solver_preconditionertype.H"
#include "solver_amgnxn_objects.H"
//#include "solver_amgnxn_vcycle.H"


namespace LINALG
{
  namespace SOLVER
  {
    namespace AMGNXN
    {
      class GenericSmoother
      {
       public:
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
        virtual void Solve(
            const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const
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
          if (not V.HasOnlyOneBlock()) dserror("We need here a single field vector");
          return;
        }
      };

      class BlockedSmoother : public GenericSmoother
      {
      };


      class BgsSmoother : public BlockedSmoother
      {
       public:
        BgsSmoother(Teuchos::RCP<BlockedMatrix> A,
            std::vector<Teuchos::RCP<GenericSmoother>> smoothers,
            std::vector<std::vector<int>> superblocks, int iter, double omega,
            std::vector<int> iters, std::vector<double> omegas)
            : A_(A),
              smoothers_(smoothers),
              superblocks_(superblocks),
              iter_(iter),
              omega_(omega),
              iters_(iters),
              omegas_(omegas)
        {
        }

        virtual void Solve(
            const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const;

       private:
        Teuchos::RCP<BlockedMatrix> A_;
        std::vector<Teuchos::RCP<GenericSmoother>> smoothers_;
        std::vector<std::vector<int>> superblocks_;
        int iter_;
        double omega_;
        std::vector<int> iters_;
        std::vector<double> omegas_;
      };

      class SimpleSmoother : public BlockedSmoother
      {
       public:
        SimpleSmoother(Teuchos::RCP<BlockedMatrix> A, Teuchos::RCP<BlockedMatrix> invApp,
            Teuchos::RCP<BlockedMatrix> Schur, Teuchos::RCP<GenericSmoother> SmooApp,
            Teuchos::RCP<GenericSmoother> SmooSchur, std::vector<int> BlocksPred,
            std::vector<int> BlocksSchur, int iter, double alpha)
            : A_(A),
              invApp_(invApp),
              Schur_(Schur),
              SmooApp_(SmooApp),
              SmooSchur_(SmooSchur),
              BlocksPred_(BlocksPred),
              BlocksSchur_(BlocksSchur),
              iter_(iter),
              alpha_(alpha)
        {
        }

        virtual void Solve(
            const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const;

       private:
        Teuchos::RCP<BlockedMatrix> A_;
        Teuchos::RCP<BlockedMatrix> invApp_;
        Teuchos::RCP<BlockedMatrix> Schur_;
        Teuchos::RCP<GenericSmoother> SmooApp_;
        Teuchos::RCP<GenericSmoother> SmooSchur_;
        std::vector<int> BlocksPred_;
        std::vector<int> BlocksSchur_;
        int iter_;
        double alpha_;
        mutable Teuchos::RCP<BlockedVector> Xp_tmp_;
        mutable Teuchos::RCP<BlockedVector> Xs_tmp_;
        mutable Teuchos::RCP<BlockedVector> Yp_tmp_;
        mutable Teuchos::RCP<BlockedVector> DYs_;
        mutable Teuchos::RCP<BlockedVector> DXp_;
        mutable Teuchos::RCP<BlockedVector> DXs_;
      };

      class MergeAndSolve : public BlockedSmoother
      {
       public:
        MergeAndSolve()
            : solver_(Teuchos::null),
              sparse_matrix_(Teuchos::null),
              block_sparse_matrix_(Teuchos::null),
              A_(Teuchos::null),
              x_(Teuchos::null),
              b_(Teuchos::null),
              isSetUp_(false)
        {
        }

        void Setup(BlockedMatrix matrix);

        virtual void Solve(
            const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const;

       private:
        Teuchos::RCP<LINALG::Solver> solver_;
        Teuchos::RCP<LINALG::SparseMatrix> sparse_matrix_;
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> block_sparse_matrix_;
        Teuchos::RCP<Epetra_Operator> A_;
        mutable Teuchos::RCP<Epetra_MultiVector> x_;
        mutable Teuchos::RCP<Epetra_MultiVector> b_;
        bool isSetUp_;
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
            const Teuchos::ParameterList& amgnxn_params,
            const Teuchos::ParameterList& smoothers_params,
            const Teuchos::ParameterList& muelu_params);

        virtual void Solve(
            const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero = false) const;


       private:
        void Setup();

        Teuchos::RCP<AMGNXN::BlockedMatrix> A_;
        std::vector<Teuchos::ParameterList> muelu_lists_;
        std::vector<int> num_pdes_;
        std::vector<int> null_spaces_dim_;
        std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data_;
        Teuchos::ParameterList amgnxn_params_;
        Teuchos::ParameterList smoothers_params_;
        Teuchos::ParameterList muelu_params_;

        bool is_setup_flag_;
        Teuchos::RCP<AMGNXN::Hierarchies> H_;
        Teuchos::RCP<AMGNXN::MonolithicHierarchy> M_;
        Teuchos::RCP<AMGNXN::Vcycle> V_;
      };


      class MueluSmootherWrapper : public SingleFieldSmoother
      {
       public:
        MueluSmootherWrapper(
            Teuchos::RCP<MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>> S)
            : S_(S)
        {
        }
        virtual void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
            bool InitialGuessIsZero = false) const;

       private:
        Teuchos::RCP<MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>> S_;
      };


      class MueluHierarchyWrapper : public SingleFieldSmoother  // Not used
      {
       public:
        MueluHierarchyWrapper(
            Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H);
        virtual void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
            bool InitialGuessIsZero = false) const;

       private:
        Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;
        Teuchos::RCP<Epetra_Operator> P_;
      };

      class MueluAMGWrapper : public SingleFieldSmoother
      {
       public:
        MueluAMGWrapper(Teuchos::RCP<SparseMatrix> A, int num_pde, int null_space_dim,
            Teuchos::RCP<std::vector<double>> null_space_data, Teuchos::ParameterList muelu_list);
        virtual void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
            bool InitialGuessIsZero = false) const;
        void Setup();

       protected:
        Teuchos::RCP<SparseMatrix> A_;
        int num_pde_;
        int null_space_dim_;
        Teuchos::RCP<std::vector<double>> null_space_data_;
        Teuchos::ParameterList muelu_list_;
        Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;
        void BuildHierarchy();

       private:
        Teuchos::RCP<Epetra_Operator> P_;
      };

      class SingleFieldAMG : public MueluAMGWrapper
      {
       public:
        SingleFieldAMG(Teuchos::RCP<SparseMatrix> A, int num_pde, int null_space_dim,
            Teuchos::RCP<std::vector<double>> null_space_data, Teuchos::ParameterList muelu_list,
            Teuchos::ParameterList fine_smoother_list);
        virtual void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
            bool InitialGuessIsZero = false) const;

       private:
        Teuchos::ParameterList fine_smoother_list_;
        Teuchos::RCP<VcycleSingle> V_;
        void Setup();
      };

      class IfpackWrapper : public SingleFieldSmoother
      {
       public:
        IfpackWrapper(Teuchos::RCP<SparseMatrixBase> A, Teuchos::ParameterList& list);
        ~IfpackWrapper() { delete prec_; }
        virtual void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
            bool InitialGuessIsZero = false) const;

       private:
        Ifpack_Preconditioner* prec_;
        Teuchos::RCP<SparseMatrixBase> A_;
        Teuchos::RCP<Epetra_RowMatrix> Arow_;
        Teuchos::ParameterList list_;
        std::string type_;
      };



      class DirectSolverWrapper : public SingleFieldSmoother
      {
       public:
        DirectSolverWrapper();
        void Setup(Teuchos::RCP<LINALG::SparseMatrix> matrix);
        virtual void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,
            bool InitialGuessIsZero = false) const;

       private:
        Teuchos::RCP<LINALG::Solver> solver_;
        Teuchos::RCP<Epetra_Operator> A_;
        mutable Teuchos::RCP<Epetra_MultiVector> x_;
        mutable Teuchos::RCP<Epetra_MultiVector> b_;
        bool isSetUp_;
      };

      // Auxiliary class to wrap the null space data to be used within the smoothers
      class NullSpaceInfo
      {
       public:
        NullSpaceInfo() {}
        NullSpaceInfo(
            int num_pdes, int null_space_dim, Teuchos::RCP<std::vector<double>> null_space_data)
            : num_pdes_(num_pdes),
              null_space_dim_(null_space_dim),
              null_space_data_(null_space_data)
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
        virtual Teuchos::RCP<GenericSmoother> Create() = 0;
      };

      // This class is able to create any smoother. The smoother to be created is given in a
      // parameter list
      class SmootherFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();

       private:
        void SetTypeAndParams();
      };


      class BgsSmootherFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();

       private:
        void ParseSmootherNames(const std::string& smoothers_string,
            std::vector<std::string>& smoothers_vector, std::vector<std::vector<int>> superblocks);
      };

      class CoupledAmgFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };

      class SimpleSmootherFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();

       private:
        Teuchos::RCP<SparseMatrix> ApproximateInverse(
            const SparseMatrixBase& A, std::string method);
        Teuchos::RCP<BlockedMatrix> ComputeSchurComplement(const BlockedMatrix& invApp,
            const BlockedMatrix& Aps, const BlockedMatrix& Asp, const BlockedMatrix& Ass);
      };

      class MergeAndSolveFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };

      class IfpackWrapperFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };

      class MueluSmootherWrapperFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };

      class HierarchyRemainderWrapperFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };

      class MueluAMGWrapperFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };

      class SingleFieldAMGFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };

      class DirectSolverWrapperFactory : public SmootherFactoryBase
      {
       public:
        virtual Teuchos::RCP<GenericSmoother> Create();
      };


    }  // namespace AMGNXN
  }    // namespace SOLVER
}  // namespace LINALG


#endif  // HAVE_MueLu
#endif  // AMGNXN_OPERATOR_H
