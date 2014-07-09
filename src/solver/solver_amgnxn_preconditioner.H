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

#include "../linalg/linalg_solver.H"
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
    /*! \brief
     * Helper class used to create an instance of BlockSparseMatrixBase
     * from a vector containing pointers to the individual blocks.
     */
    //! 
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
            bool InitialGuessIsZero=false) const=0;
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
            bool InitialGuessIsZero=false) const;

      private:
        Teuchos::RCP<SmootherBase> S_;

    }; // class NonBlockSmootherWrapperMueLu 


    //==================================================================
    /// Auxiliary class to hide an Ifpack smoother 
    //==================================================================

    class NonBlockSmootherWrapperIfpack : public NonBlockSmootherWrapperBase
    {
      public:
        NonBlockSmootherWrapperIfpack(
            Teuchos::RCP<SparseMatrixBase> A, 
            Teuchos::ParameterList& list);

        virtual void Apply(const Epetra_MultiVector& X,
            Epetra_MultiVector& Y,
            bool InitialGuessIsZero=false) const;

      private:

        Ifpack_Preconditioner* prec_;
        Teuchos::RCP<SparseMatrixBase> A_; 
        Teuchos::RCP<Epetra_RowMatrix> Arow_;
        Teuchos::ParameterList list_;
        std::string type_;



    }; // class NonBlockSmootherWrapperMueLu 


    //==================================================================
    /// Class implementing a block Gauss-Seidel smoother
    //==================================================================

    class BlockSmootherWrapperBGS : public BlockSmootherWrapperBase
    {

      public:

        BlockSmootherWrapperBGS
          (
           Teuchos::RCP<BlockSparseMatrixBase> A,
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
        Teuchos::RCP<BlockSparseMatrixBase>                       A_;  

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
    /// Class implementing a SIMPLE(C) smoother (only for 2x2) block matrices
    //==================================================================

    class BlockSmootherWrapperSIMPLE : public BlockSmootherWrapperBase
    {
      public:

        BlockSmootherWrapperSIMPLE
          (
           Teuchos::RCP<BlockSparseMatrixBase> A,
           Teuchos::RCP<NonBlockSmootherWrapperBase> Smoother_A00,
           const Teuchos::ParameterList & ParamsSmoother_S,
           int iter,
           double omega,
           bool flip
          ) :
          A_(A),Smoother_A00_(Smoother_A00),
          ParamsSmoother_S_(ParamsSmoother_S),
          iter_(iter),omega_(omega),flip_(flip),flip_vector_(2,0){ Setup();}


        virtual void Apply(const Epetra_MultiVector& X, 
            Epetra_MultiVector& Y, bool InitialGuessIsZero) const;

        void Setup();

      private:

        /// Bock matrix
        Teuchos::RCP<BlockSparseMatrixBase> A_;  

        /// Smoother for the block 0,0 (fliped order)
        Teuchos::RCP<NonBlockSmootherWrapperBase>  Smoother_A00_;  

        /// Cheap approximation of the inverse of block 0,0 (fliped order)
        Teuchos::RCP<Epetra_Vector> invA00_;

        /// Smoother for the Schur complement 
        Teuchos::RCP<NonBlockSmootherWrapperBase>  Smoother_S_;  

        /// Schur complement
        Teuchos::RCP<SparseMatrix> S_;

        /// Parameters for creating the smother for the Schur complement
        Teuchos::ParameterList  ParamsSmoother_S_;

        /// Number of sweeps
        int iter_;                   

        /// Damping 
        double omega_;                  

        /// flip order ( if false, the Schur complement is located at the last diagonal block. The other way arround if flip_==true    )
        bool flip_; 

        std::vector<int> flip_vector_;


    }; //class BlockSmootherWrapperBGS


    //==================================================================
    /// Class wrapping a linear solver for a BlockSparseMatrix 
    // It is used to solve in the coarsest level of the AMGnxn preconditioner
    //==================================================================
    
    class BlockSmootherWrapperKLU : public BlockSmootherWrapperBase
    {

      public:

        BlockSmootherWrapperKLU();

        void Setup(Teuchos::RCP<LINALG::BlockSparseMatrixBase>     matrix);

        virtual void Apply(const Epetra_MultiVector& X, 
            Epetra_MultiVector& Y, bool InitialGuessIsZero) const;

      private:

        Teuchos::RCP<LINALG::Solver>                    solver_;
        Teuchos::RCP<LINALG::SparseMatrix>              sparse_matrix_;
        Teuchos::RCP<Epetra_Operator>                   A_;
        mutable Teuchos::RCP<Epetra_MultiVector>        x_;
        mutable Teuchos::RCP<Epetra_MultiVector>        b_;
        bool isSetUp_;

    };

    //==================================================================
    ///  A general Vcycle
    //==================================================================

    class Richardson_Vcycle_Operator
    {
      public:

        Richardson_Vcycle_Operator
          (int NumLevels,int NumSweeps,double omega,
           bool activate_analysis=false, std::string file_for_analysis="AMGnxn_analysis.txt");

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

        // For analysis
        bool activate_analysis_;
        mutable Teuchos::RCP<std::ofstream> out_analysis_;
        void compute_residual_and_write
          (const Epetra_MultiVector& X,const Epetra_MultiVector& Y,
           int lev,const  std::string& str) const;

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
           std::vector<bool>                       flipPosSmoo,
           std::string                             BlockSmoother="BGS",
           bool                                    klu_on_coarsest_level=false,
           bool                                    analysis_AMG=false,
           std::string                             file_analysis_AMG="amgnxn_analysis.txt",
           Teuchos::RCP<Teuchos::ParameterList>    ParamsSmoother_Schur=Teuchos::null
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
        // Use a direct solver on the coarsest level
        bool                                      klu_on_coarsest_level_;
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
        Teuchos::ParameterList   ParamsSmoother_Schur_;
        std::string BlockSmoother_;

        /// For analysis
        bool analysis_AMG_;
        static int instance_id_;
        std::string file_analysis_AMG_;

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

        virtual void Setup
          (
           bool create,
           Epetra_Operator * matrix,
           Epetra_MultiVector * x,
           Epetra_MultiVector * b 
          );

        virtual void Setup(Teuchos::RCP<BlockSparseMatrixBase> A);


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
           std::vector<int>& offsets,
           int offsetFineLevel
          );

        // Helper to extract the Schur complement parameters
        Teuchos::RCP<Teuchos::ParameterList> Read_Schur_Smoother_params
          (std::string xmlFileName,std::string  BlockSmoother);

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
