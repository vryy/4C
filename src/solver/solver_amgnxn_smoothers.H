/*!----------------------------------------------------------------------
\file solver_amgnxn_smoothers.H

<pre>
Maintainer: Francesc Verdugo
            verdugo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/
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
#include <MueLu_UseShortNames.hpp>

// Baci includes
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "solver_preconditionertype.H"
#include "solver_amgnxn_vcycle.H"


namespace LINALG
{


  Teuchos::RCP<SparseOperator> Multiply(
      Teuchos::RCP<SparseOperator> A, bool transA,
      Teuchos::RCP<SparseOperator> B, bool transB,
      bool complite);

  namespace SOLVER
  {
    class BlockSparseMatrix_Creator
    {
      public:
        Teuchos::RCP<BlockSparseMatrixBase> CreateBlockSparseMatrix(
            std::vector< Teuchos::RCP<SparseMatrix> > blocks,
            int rows,
            int cols,
            Epetra_DataAccess access=Copy,
            bool explicitdirichlet=true,
            bool savegraph=false);
    };

    class BlockAggrupator
    {
      public:
        BlockAggrupator(
            Teuchos::RCP<BlockSparseMatrixBase> A,
            std::vector<std::vector<int> > super_blocks_to_blocks);


        std::vector< Teuchos::RCP<SparseOperator> > GetSuperBlocksRowOrder();
        Teuchos::RCP<MultiMapExtractor> GetRangeMapExtractor();
        Teuchos::RCP<MultiMapExtractor> GetDomainMapExtractor();
        int GetNumSuperBlocks();
        std::vector<std::vector<int> > GetSuperBlocksToBlocks();
        Teuchos::RCP<SparseOperator> GetSuperBlock(int srow, int scol);

      private:
        std::vector<std::vector<int> > superblocks_to_blocks_;
        Teuchos::RCP<BlockSparseMatrixBase> A_;
        std::vector< Teuchos::RCP<LINALG::SparseOperator> > superblocks_; // row major order
        Teuchos::RCP<MultiMapExtractor> range_ex_;
        Teuchos::RCP<MultiMapExtractor> domain_ex_;
        int num_superblocks_;

        void Setup();
    };


    class AMGnxn_SmootherBase
    {
      public:
        virtual ~AMGnxn_SmootherBase(){};
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const=0;
        void Richardson(
            Teuchos::RCP<const AMGnxn_SmootherBase> Ainv,
            Teuchos::RCP<const SparseOperator> A,
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            int iters,
            double omega) const;
    };


    class AMGnxn_BlockSmootherBase : public AMGnxn_SmootherBase
    {
      public:
        virtual ~AMGnxn_BlockSmootherBase(){};
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const=0;
    };


    class AMGnxn_NonBlockSmootherBase : public AMGnxn_SmootherBase
    {
      public:
        virtual ~AMGnxn_NonBlockSmootherBase(){};
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const=0;
    };


    class MueluSmootherWrapper : public AMGnxn_NonBlockSmootherBase
    {
      public:
        MueluSmootherWrapper(Teuchos::RCP<SmootherBase> S): S_(S){}
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        Teuchos::RCP<SmootherBase> S_;
    };


    class MueluHierarchyWrapper : public AMGnxn_NonBlockSmootherBase // Not used
    {
      public:
        MueluHierarchyWrapper(Teuchos::RCP<Hierarchy> H);
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        Teuchos::RCP<Hierarchy> H_;
        Teuchos::RCP<Epetra_Operator> P_;
    };

    class MueluAMGWrapper : public AMGnxn_NonBlockSmootherBase
    {
      public:
        MueluAMGWrapper(
            Teuchos::RCP<SparseMatrix> A,
            int num_pde,
            int null_space_dim,
            Teuchos::RCP<std::vector<double> > null_space_data,
            Teuchos::ParameterList muelu_list);
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        Teuchos::RCP<SparseMatrix> A_;
        int num_pde_;
        int null_space_dim_;
        Teuchos::RCP<std::vector<double> > null_space_data_;
        Teuchos::ParameterList muelu_list_;
        Teuchos::RCP<Hierarchy> H_;
        Teuchos::RCP<Epetra_Operator> P_;
        void Setup();
    };


    class HierarchyRemainderWrapper : public AMGnxn_NonBlockSmootherBase
    {
      public:
        HierarchyRemainderWrapper(Teuchos::RCP<Richardson_Vcycle_Operator> S,int start_level=0);
        virtual void Apply
          (const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero=false) const;
      private:
        int start_level_;
        Teuchos::RCP<Richardson_Vcycle_Operator> S_;
    };


    class IfpackWrapper : public AMGnxn_NonBlockSmootherBase
    {
      public:
        IfpackWrapper(
            Teuchos::RCP<SparseMatrixBase> A,
            Teuchos::ParameterList& list);
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        Ifpack_Preconditioner* prec_;
        Teuchos::RCP<SparseMatrixBase> A_;
        Teuchos::RCP<Epetra_RowMatrix> Arow_;
        Teuchos::ParameterList list_;
        std::string type_;
    };



    class DirectSolverWrapper : public AMGnxn_BlockSmootherBase
    {
      public:
        DirectSolverWrapper();
        void Setup(Teuchos::RCP<LINALG::SparseMatrix>     matrix);
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        Teuchos::RCP<LINALG::Solver>                    solver_;
        Teuchos::RCP<Epetra_Operator>                   A_;
        mutable Teuchos::RCP<Epetra_MultiVector>        x_;
        mutable Teuchos::RCP<Epetra_MultiVector>        b_;
        bool isSetUp_;
    };

    class MergeAndSolve : public AMGnxn_BlockSmootherBase
    {
      public:
        MergeAndSolve();
        void Setup(Teuchos::RCP<LINALG::BlockSparseMatrixBase>     matrix);
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        Teuchos::RCP<LINALG::Solver>                    solver_;
        Teuchos::RCP<LINALG::SparseMatrix>              sparse_matrix_;
        Teuchos::RCP<Epetra_Operator>                   A_;
        mutable Teuchos::RCP<Epetra_MultiVector>        x_;
        mutable Teuchos::RCP<Epetra_MultiVector>        b_;
        bool isSetUp_;
    };


    class SIMPLE_BlockSmoother: public AMGnxn_BlockSmootherBase
    {
      public:
        SIMPLE_BlockSmoother(
            Teuchos::RCP<LINALG::SparseOperator>  App,
            Teuchos::RCP<LINALG::SparseOperator>  Ass,
            Teuchos::RCP<LINALG::SparseOperator>  Aps,
            Teuchos::RCP<LINALG::SparseOperator>  Asp,
            Teuchos::RCP<LINALG::SparseOperator>  invApp,
            Teuchos::RCP<LINALG::SparseOperator>  S,
            Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_App,
            Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_S,
            Teuchos::RCP<MultiMapExtractor> range_ex,
            Teuchos::RCP<MultiMapExtractor> domain_ex,
            int p,
            int s,
            int iter,
            double alpha=1.0,
            std::string correction = "approximated inverse");
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
        void Setup();
      private:
        Teuchos::RCP<LINALG::SparseOperator>  App_;
        Teuchos::RCP<LINALG::SparseOperator>  Ass_;
        Teuchos::RCP<LINALG::SparseOperator>  Aps_;
        Teuchos::RCP<LINALG::SparseOperator>  Asp_;
        Teuchos::RCP<LINALG::SparseOperator>  invApp_;
        Teuchos::RCP<LINALG::SparseOperator>  S_;
        Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_App_;
        Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_S_;
        Teuchos::RCP<MultiMapExtractor> range_ex_;
        Teuchos::RCP<MultiMapExtractor> domain_ex_;
        int p_;
        int s_;
        int iter_;
        double alpha_;
        std::string correction_;
    };


    class BGS_BlockSmoother : public AMGnxn_BlockSmootherBase
    {
      public:
        BGS_BlockSmoother(
            std::vector< Teuchos::RCP<LINALG::SparseOperator> > blocks,
            std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > smoothers,
            std::vector<int> indices,
            Teuchos::RCP<MultiMapExtractor> range_ex,
            Teuchos::RCP<MultiMapExtractor> domain_ex,
            int size,
            int iter,
            double omega,
            std::vector<int> iters,
            std::vector<double> omegas
            );
        virtual void Apply(
            const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        std::vector< Teuchos::RCP<LINALG::SparseOperator> > blocks_; // row major order
        std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > smoothers_; // only for diagonal blocks
        std::vector<int> indices_; // permutation for recover info in range_ex_ and domain_ex_
        Teuchos::RCP<MultiMapExtractor> range_ex_;
        Teuchos::RCP<MultiMapExtractor> domain_ex_;
        int size_;
        int iter_;
        double omega_;
        std::vector<int> iters_;
        std::vector<double> omegas_;
    };

    class AMGnxn_Operator; // Forward declaration

    class AMG_BlockSmoother : public AMGnxn_BlockSmootherBase
    {

      public:

        AMG_BlockSmoother(
            Teuchos::RCP<BlockSparseMatrixBase> A,
            std::vector<int> num_pdes,
            std::vector<int> null_spaces_dim,
            std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
            const Teuchos::ParameterList& amgnxn_params,
            const Teuchos::ParameterList& smoothers_params);

        virtual void Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;
      private:
        Teuchos::RCP<LINALG::SOLVER::AMGnxn_Operator>  P_;



    };


    // Auxiliary class to wrap the null space data to be used within the smoothers
    class NullSpaceInfo
    {
      public:
        NullSpaceInfo(){}
        NullSpaceInfo(
            int num_pdes,
            int null_space_dim,
            Teuchos::RCP<std::vector<double> >  null_space_data);
        int GetNumPDEs(){return num_pdes_;};
        int GetNullSpaceDim(){return null_space_dim_;};
        Teuchos::RCP<std::vector<double> > GetNullSpaceData(){return null_space_data_;};
      private:
        int num_pdes_;
        int null_space_dim_;
        Teuchos::RCP<std::vector<double> >  null_space_data_;
    };

    class AMGnxn_Hierarchies;

    class AMGnxn_SmootherFactoryBase
    {
      public:
        AMGnxn_SmootherFactoryBase();
        Teuchos::RCP<SparseOperator>     GetOperator();
        Teuchos::ParameterList           GetParams();
        Teuchos::ParameterList           GetParamsSmoother();
        Teuchos::RCP<AMGnxn_Hierarchies> GetHierarchies();
        int                              GetLevel();
        int                              GetBlock();
        std::vector<int>                 GetBlocks();
        std::string                      GetSmootherName();
        std::string                      GetType();
        std::string                      GetVerbosity();
        NullSpaceInfo                    GetNullSpace();
        std::vector<NullSpaceInfo>       GetNullSpaceAllBlocks();

        void SetOperator           (Teuchos::RCP<SparseOperator>      in);
        void SetParams             (const Teuchos::ParameterList&     in);
        void SetParamsSmoother     (const Teuchos::ParameterList&     in);
        void SetHierarchies        (Teuchos::RCP<AMGnxn_Hierarchies>  in);
        void SetLevel              (int                               in);
        void SetBlock              (int                               in);
        void SetBlocks             (std::vector<int>                  in);
        void SetSmootherName       (std::string                       in);
        void SetType               (std::string                       in);
        void SetVerbosity          (std::string                       in);
        void SetNullSpace          (const NullSpaceInfo&              in);
        void SetNullSpaceAllBlocks (const std::vector<NullSpaceInfo>& in);

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

        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create()=0;

      private:
        Teuchos::RCP<SparseOperator>     operator_;
        Teuchos::ParameterList           params_;
        Teuchos::ParameterList           params_subsolver_;
        Teuchos::RCP<AMGnxn_Hierarchies> hierarchies_;
        int                              level_;
        int                              block_;
        std::vector<int>                 blocks_;
        std::string                      subsolver_name_;
        std::string                      type_;
        std::string                      verbosity_;
        NullSpaceInfo                    null_space_;
        std::vector<NullSpaceInfo>       null_space_all_blocks_;

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

    // This class creates any type of smoother.
    // The smoother to be constructed is specified given a list containing parameters
    // for one or several smoothers using SetParamsSmoother()
    // The specific smoother to be constructed is defined by the function SetSmootherName()
    class AMGnxn_SmootherFactory : public AMGnxn_SmootherFactoryBase {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
      private:
        void SetTypeAndParams();
    };

    // The classes derived from this one, create each one only a particular type of smoother
    class AMGnxn_ParticularSmootherFactory : public AMGnxn_SmootherFactoryBase {
      protected:
        void ParseBlocks(
            const std::string& block_string,
            const std::vector<int>& blocks,
            std::vector<std::vector<int> >& superblocks_to_blocks,
            std::vector<std::vector<int> >& superblocks_to_blocks_local);

    };

    class BGS_BlockSmootherFactory : public  AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
      private:
        void ParseSmootherNames(
            const std::string& smoothers_string,
            std::vector<std::string>& smoothers_vector,std::vector<std::vector<int> > superblocks);
    };

    class AMG_BlockSmootherFactory : public  AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
    };

    class SIMPLE_BlockSmootherFactory : public AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
      private:
        Teuchos::RCP<SparseMatrix> ApproximateInverse(const SparseMatrix& A,std::string method);
        Teuchos::RCP<LINALG::SparseOperator> ComputeSchurComplement(
            Teuchos::RCP<SparseOperator> invApp,
            Teuchos::RCP<SparseOperator> Aps,
            Teuchos::RCP<SparseOperator> Asp,
            Teuchos::RCP<SparseOperator> Ass);
    };

    class MergeAndSolveFactory : public AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
    };

    class IfpackWrapperFactory : public AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
    };


    class MueluSmootherWrapperFactory : public AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
    };

    class HierarchyRemainderWrapperFactory : public AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
    };

    class MueluAMGWrapperFactory : public AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
    };

    class DirectSolverWrapperFactory : public AMGnxn_ParticularSmootherFactory
    {
      public:
        virtual Teuchos::RCP<AMGnxn_SmootherBase> Create();
    };

  } //namespace LINALG
} // namespace SOLVER

#endif // HAVE_MueLu
#endif // AMGNXN_OPERATOR_H
