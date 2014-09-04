/*!----------------------------------------------------------------------
\file solver_amgnxn_hierarchies.H

<pre>
Maintainer: Francesc Verdugo
            verdugo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Aug 11, 2014
</pre>
*----------------------------------------------------------------------*/
#ifndef AMGNXN_HIERARCHIES_H
#define AMGNXN_HIERARCHIES_H
#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q1_2014

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Teuchos_RCP.hpp"
#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "solver_preconditionertype.H"
#include "solver_amgnxn_smoothers.H"

namespace LINALG
{
  namespace SOLVER
  {
    class AMGnxn_Hierarchies
    {
      public:
        AMGnxn_Hierarchies(
            Teuchos::RCP<BlockSparseMatrixBase> A,
            std::vector<Teuchos::ParameterList> muelu_params,
            std::vector<int> num_pdes,
            std::vector<int> null_spaces_dim,
            std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
            int NumLevelAMG,
            std::string verbosity="off");

        int GetNumLevelMin();
        int GetNumBlocks();
        int GetNumLevels(int block);

        Teuchos::RCP<BlockSparseMatrixBase> GetBlockMatrix();
        Teuchos::RCP<Hierarchy> GetH(int block);

        Teuchos::RCP<SparseMatrix> GetA(int block, int level);
        Teuchos::RCP<SparseMatrix> GetP(int block, int level);
        Teuchos::RCP<SparseMatrix> GetR(int block, int level);
        Teuchos::RCP<MueluSmootherWrapper> GetSPre(int block, int level);
        Teuchos::RCP<MueluSmootherWrapper> GetSPos(int block, int level);

        std::vector<Teuchos::RCP<SparseMatrix> > GetA(int block);
        std::vector<Teuchos::RCP<SparseMatrix> > GetP(int block);
        std::vector<Teuchos::RCP<SparseMatrix> > GetR(int block);
        std::vector<Teuchos::RCP<MueluSmootherWrapper> > GetSPre(int block);
        std::vector<Teuchos::RCP<MueluSmootherWrapper> > GetSPos(int block);

        int GetNumPDEs(int block);
        int GetNullSpaceDim(int block);
        Teuchos::RCP<std::vector<double> > GetNullSpaceData(int block);

      private:
        Teuchos::RCP<BlockSparseMatrixBase> A_;
        std::vector<Teuchos::ParameterList> muelu_params_;
        std::vector<int> num_pdes_;
        std::vector<int> null_spaces_dim_;
        std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data_;
        int NumBlocks_;
        int NumLevelMax_;
        int NumLevelMin_;
        int NumLevelAMG_;
        std::vector<Teuchos::RCP<Hierarchy> > H_block_;
        std::vector<std::vector<Teuchos::RCP<SparseMatrix> > > A_block_level_;
        std::vector<std::vector<Teuchos::RCP<SparseMatrix> > > P_block_level_;
        std::vector<std::vector<Teuchos::RCP<SparseMatrix> > > R_block_level_;
        std::vector<std::vector<Teuchos::RCP<MueluSmootherWrapper> > > SPre_block_level_;
        std::vector<std::vector<Teuchos::RCP<MueluSmootherWrapper> > > SPos_block_level_;
        std::string verbosity_;

        void Setup();

        Teuchos::RCP<Hierarchy> BuildMueLuHierarchy(
            Teuchos::ParameterList paramListFromXml,
            int num_pde,
            int null_space_dim,
            Teuchos::RCP<std::vector<double> > null_space_data,
            Teuchos::RCP<Epetra_Operator> A_eop,
            int block,
            int NumBlocks,
            std::vector<int>& offsets,
            int offsetFineLevel);

        std::string ConvertInt(int number)
        {
          std::stringstream ss;
          ss << number;
          return ss.str();
        }

    }; //AMGnxn_Hierarchies


    class AMGnxn_MonolithicHierarchy
    {
      public:
        AMGnxn_MonolithicHierarchy(
            Teuchos::RCP<AMGnxn_Hierarchies> H,
            const Teuchos::ParameterList& params,
            const Teuchos::ParameterList& params_smoothers);

        int GetNumLevels();

        Teuchos::RCP<AMGnxn_Hierarchies> GetHierarchies();

        Teuchos::RCP<BlockSparseMatrixBase> GetA(int level);

        Teuchos::RCP<Richardson_Vcycle_Operator> BuildVCycle();

      private:
        Teuchos::RCP<AMGnxn_Hierarchies> H_;
        int NumLevels_;
        int NumBlocks_;
        std::vector< Teuchos::RCP<BlockSparseMatrixBase> > A_;
        std::vector< Teuchos::RCP<BlockSparseMatrixBase> > P_;
        std::vector< Teuchos::RCP<BlockSparseMatrixBase> > R_;
        std::vector< Teuchos::RCP<AMGnxn_BlockSmootherBase> > Spre_;
        std::vector< Teuchos::RCP<AMGnxn_BlockSmootherBase> > Spos_;
        Teuchos::ParameterList params_;
        Teuchos::ParameterList params_smoothers_;

        void Setup();
        Teuchos::RCP<AMGnxn_BlockSmootherBase> BuildSmoother(int level);

    };





  } //namespace LINALG
} // namespace SOLVER

#endif // HAVE_MueLu
#endif // HAVE_Trilinos_Q1_2014
#endif // AMGNXN_HIERARCHIES_H
