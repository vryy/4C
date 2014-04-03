/*!----------------------------------------------------------------------
\file solver_amgnxn_preconditioner.H

<pre>
Maintainer: Francesc Verdugo 
            verdugo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/
#ifndef AMGNXN_OPERATOR_H
#define AMGNXN_OPERATOR_H

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

#include "../linalg/linalg_blocksparsematrix.H"
#include "solver_preconditionertype.H"

//=========================================================================
// This might be moved to another place (begin)
//=========================================================================

namespace LINALG
{
namespace SOLVER 
{

/// Very simple matrix container implemented using std::vector(s).
//! You are welcome to add further member functions if you need them!
template<class T>
class matrix
{
  public:
  matrix(){};
  matrix(int n,int m);
  matrix(int n,int m,const T &val);
  void resize(int n,int m);
  void assign(int n,int m,const T &val);
  int nrow() const;
  int ncol() const;
  std::vector<T>& operator[](int i);
  const std::vector<T>& operator[](int i) const;
  matrix& operator=(const matrix& A);
  private:
  std::vector<std::vector<T> >  matrix_;
}; // class matrix

} //SOLVER 
}//LINALG


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
LINALG::SOLVER::matrix<T>::matrix(int n,int m)
{
  resize(n,m);
}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
LINALG::SOLVER::matrix<T>::matrix(int n,int m,const T &val)
{
  assign(n,m,val);
}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
void LINALG::SOLVER::matrix<T>::resize(int n,int m)
{
  matrix_.resize(n);
  for(int i=0;i<n;i++)
    matrix_[i].resize(m);
  return;
}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
void LINALG::SOLVER::matrix<T>::assign(int n,int m,const T &val)
{
  matrix_.resize(n);
  for(int i=0;i<n;i++)
    matrix_[i].assign(m,val);
}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
int LINALG::SOLVER::matrix<T>::nrow() const
{
  return matrix_.size();
}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
int LINALG::SOLVER::matrix<T>::ncol() const
{
  int ncol;
  if(matrix_.size()>0)
    ncol=matrix_[0].size();
  else
    ncol=0;
  return ncol;
}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
std::vector<T>& LINALG::SOLVER::matrix<T>::operator[](int i)
{
  return matrix_[i];
}
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
const std::vector<T>& LINALG::SOLVER::matrix<T>::operator[](int i) const
{
  return matrix_[i];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
template<class T>
LINALG::SOLVER::matrix<T>& LINALG::SOLVER::matrix<T>::operator=(const LINALG::SOLVER::matrix<T>& A)
{
  if(this == &A)
       return *this;

  int n=A.nrow();
  int m=A.ncol();
  this->resize(n,m);
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      matrix_[i][j]=A[i][j];

  return *this;

}

//=========================================================================
// This might be moved to another place (end)
//=========================================================================


namespace LINALG
{
  namespace SOLVER
  {

    //==================================================================
    /// Auxiliary class implementing a minimal Blocked sparse matrix.
    /*!
     * This class allows to build up a block matrix from simple matrices whose maps may begin with 0 
     * (or any other arbitrary value). This class is used only to bypass the fact that muelu 
     * furnishes matrices  whose maps begin with 0 in some cases. 
     * Therefore we cannot use these matrices to build a 
     * standard "BlockSparseMatrixBase" object.
     *
     * Don't use this class. We want to remove it in the future.
     */
    //==================================================================

    class BlockSparseMatrixAUX : virtual public Epetra_Operator
    {

      public:

        BlockSparseMatrixAUX( 
            std::vector< Teuchos::RCP<SparseMatrix> > blocks,
            int rows,
            int cols,
            Epetra_DataAccess access=Copy);

        int Cols() const {return cols_;}

        int Rows() const {return rows_;}

        const SparseMatrix& Matrix(int i, int j) const {return *blocks_[i*Cols()+j];}

        int ApplyBlock(int r,int c, const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        const MultiMapExtractor& DomainExtractor() const { return *domainmaps_; }

        const MultiMapExtractor& RangeExtractor() const { return *rangemaps_; }

        Teuchos::RCP<Epetra_Import> DomainImporters(int i) const {return domainimporters_[i];}

        Teuchos::RCP<Epetra_Import> RangeImporters(int i) const {return rangeimporters_[i];}

        virtual int   SetUseTranspose (bool UseTranspose)
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual int   Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        virtual int   ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual double   NormInf () const
        {
          dserror("Function not implemented");
          return -1.0;
        }

        virtual const char *   Label () const
        {
          return "BlockSparseMatrixAUX";
        }

        virtual bool   UseTranspose () const
        {
          dserror("Function not implemented");
          return false;
        }
        virtual bool   HasNormInf () const
        {
          dserror("Function not implemented");
          return false;
        }
        virtual const Epetra_Comm &   Comm () const {return blocks_[0]->Comm();} 

        virtual const Epetra_Map &   OperatorDomainMap () const { return *domainmaps_->FullMap(); } 

        virtual const Epetra_Map &   OperatorRangeMap () const { return *rangemaps_->FullMap(); }


      private:

        void Setup(std::vector< Teuchos::RCP<SparseMatrix> > blocks,Epetra_DataAccess access);
        void FillWithZeroOfDiagonalBlocks(
            const std::vector<Teuchos::RCP<SparseMatrix> >& blocks,
            std::vector<Teuchos::RCP<SparseMatrix> >& blocks_all);
        Teuchos::RCP<Epetra_Map> ComputeShiftedMap(const Epetra_Map& Map,int shift);
        const Epetra_Map& DomainMapUnShifted(int col) const {return Matrix(0,col).DomainMap();}
        const Epetra_Map& RangeMapUnShifted (int row) const {return Matrix(row,0).RangeMap() ;}
        Epetra_Map CopyMapByHand(const Epetra_Map& MapIn);

        int rows_;
        int cols_;
        std::vector< Teuchos::RCP<SparseMatrix>  > blocks_;
        Teuchos::RCP<MultiMapExtractor> domainmaps_;
        Teuchos::RCP<MultiMapExtractor> rangemaps_;
        std::vector<Teuchos::RCP<Epetra_Import> > domainimporters_;
        std::vector<Teuchos::RCP<Epetra_Import> > rangeimporters_;

    };


    //==================================================================
    /*! \brief
     * Helper class used to create an instance of BlockSparseMatrixBase
     * from a vector containing pointers to the individual blocks.
     */
    //! Not used (yet)
    //==================================================================

    class FactoryBlockSparseMatrix
    {
      public:
        /*!
          \param blocks Vector containing the individual blocks using a ROW MAJOR 
          matrix storage (i.e sweeping by rows). The number of given
          blocks should be rows*cols. The range and domain maps
          of the blocks should be consistent! If rows==cols and only $(rows) blocks
          are given, the blocks are interpreted as the diagonal blocks. The off
          diagonal blocks are assumed to be zero.
          \param rows Number of block rows
          \param cols Number of block columns
          \param access Allows to create the block matrix as deep or shalow copies of the blocks
          \param explicitdirichlet whether to remove Dirichlet zeros from the
          matrix graphs in each block
          \param savegraph whether to save the matrix graphs of each block and
          recreate filled matrices the next time
          */
        Teuchos::RCP<BlockSparseMatrixBase> CreateBlockSparseMatrix
          (
           std::vector< Teuchos::RCP<SparseMatrix> > blocks,
           int rows,
           int cols,
           Epetra_DataAccess access=Copy,
           bool explicitdirichlet=true,
           bool savegraph=false
          );

    };

    //==================================================================
    /*! \brief
     * Interface allowing to use different objects as smoothers 
     * in the Richardson_Vcycle_Operator class
     */
    //==================================================================

    class SmootherWrapperBase
    {
      public:
        virtual ~SmootherWrapperBase(){};
        virtual void Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero) const=0;
    };

    //==================================================================
    /*! \brief
     * Interface allowing to use different block smoothers in the same
     * Richardson_Vcycle_Operator
     */
    //==================================================================

    class BlockSmootherWrapperBase : public SmootherWrapperBase
    {
      public:
        virtual ~BlockSmootherWrapperBase(){};
        virtual void Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero) const=0;
    };

    //==================================================================
    /*! \brief
     * Interface allowing to use different single field smoothers 
     * inside a block smoother 
     */
    //==================================================================

    class NonBlockSmootherWrapperBase : public SmootherWrapperBase
    {
      public:
        virtual ~NonBlockSmootherWrapperBase(){};
        virtual void Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero) const=0;
    };

    //==================================================================
    /*! \brief
     * Auxiliary class to bypass the problem with the maps in the
     * smoothers
     */
    //==================================================================

    class NonBlockSmootherAUX : public NonBlockSmootherWrapperBase
    {
      public:
        NonBlockSmootherAUX(
            Teuchos::RCP<NonBlockSmootherWrapperBase> S,
            Teuchos::RCP<Epetra_Import> DomainImporter,
            Teuchos::RCP<Epetra_Import> RangeImporter)
          : S_(S), DomainImporter_(DomainImporter), RangeImporter_(RangeImporter) {}

        virtual void Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero) const
        {
          Epetra_MultiVector XU(DomainImporter_->TargetMap(),X.NumVectors());
          Epetra_MultiVector YU( RangeImporter_->SourceMap(),Y.NumVectors());
          XU.Import(X,*DomainImporter_,Insert);
          S_->Apply(XU,YU,InitialGuessIsZero);
          Y.Import(YU,*RangeImporter_,Insert);
        }

      private:
        Teuchos::RCP<NonBlockSmootherWrapperBase> S_;
        Teuchos::RCP<Epetra_Import> DomainImporter_;
        Teuchos::RCP<Epetra_Import> RangeImporter_;
    };


    //==================================================================
    /// Auxiliary class to hide Xpetra in the MueLu smoother
    //==================================================================

    class NonBlockSmootherWrapperMueLu : public NonBlockSmootherWrapperBase
    {
      public:
        NonBlockSmootherWrapperMueLu(Teuchos::RCP<SmootherBase> S): S_(S){}

        virtual void Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero) const
        {

          // Convert to Xpetra
          Teuchos::RCP<Epetra_MultiVector> X_rcp =
            Teuchos::rcp(new Epetra_MultiVector(X));
          // TODO Deep copy!!!!!!! May be there is a better way to do this
          Teuchos::RCP<Xpetra::EpetraMultiVector> Xex =
            Teuchos::rcp(new Xpetra::EpetraMultiVector(X_rcp));
          Teuchos::RCP<MultiVector> Xx =
            Teuchos::rcp_dynamic_cast<MultiVector>(Xex);
          Teuchos::RCP<Epetra_MultiVector> Y_rcp =
            Teuchos::rcp(new Epetra_MultiVector(Y));
          //TODO Deep copy!!!!!!! May be there is a better way to do this
          Teuchos::RCP<Xpetra::EpetraMultiVector> Yex =
            Teuchos::rcp(new Xpetra::EpetraMultiVector(Y_rcp));
          Teuchos::RCP<MultiVector> Yx = Teuchos::rcp_dynamic_cast<MultiVector>(Yex);
          // Apply underlying smoother
          S_->Apply(*Yx,*Xx,InitialGuessIsZero);
          // Convert to Epetra
          const Teuchos::RCP<Epetra_MultiVector>& Ye =
            MueLu::Utils<double,int,int,Node,LocalMatOps>::MV2NonConstEpetraMV(Yx);
          Y = *Ye;

          return;
        }

      private:
        Teuchos::RCP<SmootherBase> S_;

    }; // class NonBlockSmootherWrapperMueLu 


    /// Class implementing a block Gauss-Seidel smoother
    class BlockSmootherWrapperBGS : public BlockSmootherWrapperBase
    {

      public:

        BlockSmootherWrapperBGS
          (
           //Teuchos::RCP<BlockSparseMatrixBase> A,
           Teuchos::RCP<BlockSparseMatrixAUX> A,
           std::vector< Teuchos::RCP<NonBlockSmootherWrapperBase> > S,
           int global_iter,
           double global_omega,
           bool flip_order=false
          );

        virtual void Apply(const Epetra_MultiVector& X, 
            Epetra_MultiVector& Y, bool InitialGuessIsZero) const;

      private:

        /// A Richardson iteration using BGS
        void RichardsonBGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

        /// Bock matrix
        //Teuchos::RCP<BlockSparseMatrixBase>                       A_;  
        Teuchos::RCP<BlockSparseMatrixAUX>                        A_;  

        /// Smoothers for the diagonal blocks 
        std::vector< Teuchos::RCP<NonBlockSmootherWrapperBase> >  S_;  

        /// Number of GLOBAL Richardson sweeps
        int                            global_iter_;                   

        /// Damping in the GLOBAL Richardson update
        double                         global_omega_;                  

        /// 0,1,...,n-1 or n-1,n-2,...,0 depending on the flip order
        std::vector<int>               index_order_;                   

        /// Number of matrix blocks
        int                            NumBlocks_;                     

    }; //class BlockSmootherWrapperBGS


    //==================================================================
    ///  A general Vcycle
    //==================================================================

    class Richardson_Vcycle_Operator
    {
      public:

        Richardson_Vcycle_Operator(int NumLevels,int NumSweeps,double omega);

        void SetOperators  (std::vector< Teuchos::RCP<Epetra_Operator> > Avec); 
        void SetProjectors (std::vector< Teuchos::RCP<Epetra_Operator> > Pvec);
        void SetRestrictors(std::vector< Teuchos::RCP<Epetra_Operator> > Rvec);
        void SetPreSmoothers (std::vector< Teuchos::RCP<SmootherWrapperBase> > SvecPre);
        void SetPosSmoothers (std::vector< Teuchos::RCP<SmootherWrapperBase> > SvecPos);

        void Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y,int start_level=0) const;
        void Richardson_Vcycle
          (const Epetra_MultiVector& X, Epetra_MultiVector& Y,int start_level) const;
        void Vcycle(const Epetra_MultiVector& X, 
            Epetra_MultiVector& Y, int level, bool InitialGuessIsZero) const;


      private:

        int NumLevels_;
        int NumSweeps_;
        double omega_;

        std::vector< Teuchos::RCP<Epetra_Operator> > Avec_;
        std::vector< Teuchos::RCP<Epetra_Operator> > Pvec_;
        std::vector< Teuchos::RCP<Epetra_Operator> > Rvec_;
        std::vector< Teuchos::RCP<SmootherWrapperBase> > SvecPre_;
        std::vector< Teuchos::RCP<SmootherWrapperBase> > SvecPos_;

        bool flag_set_up_A_;
        bool flag_set_up_P_;
        bool flag_set_up_R_;
        bool flag_set_up_Pre_;
        bool flag_set_up_Pos_;

    }; // class Richardson_Vcycle_Operator


    //==================================================================
    ///  A Vcycle seen as an smoother
    //==================================================================

    class SmootherWrapperVcycle : public NonBlockSmootherWrapperBase
    {
      public:
        SmootherWrapperVcycle(Teuchos::RCP<Richardson_Vcycle_Operator> S,int start_level=0):
          start_level_(start_level),S_(S) {}
        virtual void Apply
          (const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
          {

            S_->Apply(X,Y,start_level_);
            return;
          }
      private:
        int start_level_;
        Teuchos::RCP<Richardson_Vcycle_Operator> S_;
    };

    //==================================================================
    /// The AMGnxn preconditioner operator
    //==================================================================

    class AMGnxn_Operator : virtual public Epetra_Operator
    {
      public:
        AMGnxn_Operator
          (
           std::vector< Teuchos::RCP<Hierarchy>  > H,
           Teuchos::RCP<BlockSparseMatrixBase>     A,
           int                                     NumLevelAMG,
           int                                     NumSweepsAMG,
           double                                  omegaAMG,
           std::vector<int>                        NumSweepsPreSmoo,
           std::vector<double>                     omegaPreSmoo,     
           std::vector<bool>                       flipPreSmoo,      
           std::vector<int>                        NumSweepsPosSmoo,  
           std::vector<double>                     omegaPosSmoo,     
           std::vector<bool>                       flipPosSmoo      
          );


        void SetUp();

        // virtual functions given by Epetra_Operator. The only one to be used is ApplyInverse()
        virtual int   ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        virtual int   SetUseTranspose (bool UseTranspose)
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual int   Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
        {
          dserror("Function not implemented");
          return -1;
        }

        virtual double   NormInf () const
        {
          dserror("Function not implemented");
          return -1.0;
        }

        virtual const char *   Label () const
        {
          return "AMGnxn_Operator";
        }

        virtual bool   UseTranspose () const
        {
          dserror("Function not implemented");
          return false;
        }
        virtual bool   HasNormInf () const
        {
          dserror("Function not implemented");
          return false;
        }
        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Comm &   Comm () const {return A_->Comm();} 

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorDomainMap () const {return A_->OperatorDomainMap();} 

        // Only required to properly define an Epetra_Operator, not should be used!
        virtual const Epetra_Map &   OperatorRangeMap () const {return A_->OperatorRangeMap();}


      private:
        /// MueLu Hierarchies for each one of the fields
        std::vector< Teuchos::RCP<Hierarchy>  >   H_;                
        /// System matrix
        Teuchos::RCP<BlockSparseMatrixBase>       A_;                
        /// Number of field or blocks in the block matrix
        int                                       NumBlocks_;        
        /// Number of sweeps in the global AMG preconditioner
        int                                       NumSweepsAMG_;     
        /// Damping parameter for the Richardson iteration using the global AMG preconditioner
        double                                    omegaAMG_;         
        /// Number of sweeps for block pre-smoothers for each level
        std::vector<int>                          NumSweepsPreSmoo_;  
        /// Damping for Block pre-smoother for each level
        std::vector<double>                       omegaPreSmoo_;     
        /// Flip the BGS pre-smoother in each level
        std::vector<bool>                         flipPreSmoo_;      
        /// Number of sweeps for block post-smoothers for each level
        std::vector<int>                          NumSweepsPosSmoo_;  
        /// Damping for Block post-smoother for each level
        std::vector<double>                       omegaPosSmoo_;     
        /// Flip the BGS post-smoother in each level
        std::vector<bool>                         flipPosSmoo_;      
        /// Flag indicating if the preconditioner can be used
        bool                                      is_setup_flag_;    
        /// Where the AMG things happen
        Teuchos::RCP<Richardson_Vcycle_Operator>  P_;                
        int                                       NumLevelMax_;
        int                                       NumLevelMin_;
        int                                       NumLevelAMG_;
        matrix<Teuchos::RCP<SparseMatrix> >    ALocal_;
        matrix<Teuchos::RCP<SparseMatrix> >    PLocal_;
        matrix<Teuchos::RCP<SparseMatrix> >    RLocal_;    
        matrix<Teuchos::RCP<LINALG::SOLVER::NonBlockSmootherWrapperMueLu> > SPreLocal_;
        matrix<Teuchos::RCP<LINALG::SOLVER::NonBlockSmootherWrapperMueLu> > SPosLocal_;
        Teuchos::RCP<Richardson_Vcycle_Operator> 
          CreateRemainingHierarchy(int start_level,int num_levels,int block);

    }; // class AMGnxn_Operator



    //==================================================================
    /// The AMGnxn preconditioner. 
    //==================================================================

    class AMGnxn_Preconditioner : public PreconditionerType
    {
      public:
        AMGnxn_Preconditioner( FILE * outfile, Teuchos::ParameterList & params);
        // TODO remove this destructor.

        virtual void Setup
          (
           bool create,
           Epetra_Operator * matrix,
           Epetra_MultiVector * x,
           Epetra_MultiVector * b 
          );

        virtual Epetra_Operator * PrecOperator() const ;
        virtual void Print( std::ostream & stream ){stream << "AMGnxn";}
        virtual const std::string getParameterListName() const { return "AMGnxn Parameters"; }

      private:

        // Helper function to extract muelu hierarchy 
        Teuchos::RCP<Hierarchy> BuildMueLuHierarchy
          (
           Teuchos::ParameterList& mllist,
           Teuchos::RCP<Epetra_Operator> A_eop,
           int block,
           int NumBlocks,
           std::vector<int>& offsets
          );

        // Helper function to convert int to std::string 
        std::string ConvertInt(int number)
        {
          std::stringstream ss;
          ss << number;
          return ss.str();  
        }

        // Private variables
        Teuchos::RCP<Epetra_Operator> P_; // The underlying preconditioner object
        Teuchos::RCP<BlockSparseMatrixBase> A_; // A own copy of the system matrix
        Teuchos::ParameterList & params_;

    };//AMGnxn_Preconditioner




  } //namespace LINALG
} // namespace SOLVER

#endif // HAVE_MueLu
#endif // HAVE_Trilinos_Q1_2014
#endif // AMGNXN_OPERATOR_H
