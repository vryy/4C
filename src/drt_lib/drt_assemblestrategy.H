#ifndef DRT_ASSEMBLESTRATEGY_H
#define DRT_ASSEMBLESTRATEGY_H

#include <vector>

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_Vector.h>

namespace LINALG
{
  class SparseOperator;
}

namespace DRT
{
  class Discretization;

  /// General assembly strategy
  class AssembleStrategy
  {
  public:

    /// Construct
    AssembleStrategy( Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                      Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                      Teuchos::RCP<Epetra_Vector>          systemvector1,
                      Teuchos::RCP<Epetra_Vector>          systemvector2,
                      Teuchos::RCP<Epetra_Vector>          systemvector3 );

    virtual ~AssembleStrategy();

    /// Create matrix graph to create a filled matrix
    virtual Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( DRT::Discretization & dis, Teuchos::RCP<const Epetra_Map> dbcmap );

    //! @name Assembly Flags
    bool Assemblemat1() { return systemmatrix1_!=Teuchos::null; }
    bool Assemblemat2() { return systemmatrix2_!=Teuchos::null; }
    bool Assemblevec1() { return systemvector1_!=Teuchos::null; }
    bool Assemblevec2() { return systemvector2_!=Teuchos::null; }
    bool Assemblevec3() { return systemvector3_!=Teuchos::null; }
    //@}

    //! @name Access Methods to Global Object
    Teuchos::RCP<LINALG::SparseOperator> Systemmatrix1() { return systemmatrix1_; }
    Teuchos::RCP<LINALG::SparseOperator> Systemmatrix2() { return systemmatrix2_; }
    Teuchos::RCP<Epetra_Vector>          Systemvector1() { return systemvector1_; }
    Teuchos::RCP<Epetra_Vector>          Systemvector2() { return systemvector2_; }
    Teuchos::RCP<Epetra_Vector>          Systemvector3() { return systemvector3_; }
    //@}

    //! @name Access Methods to Element Local Object
    Epetra_SerialDenseMatrix & Elematrix1() { return elematrix1_; }
    Epetra_SerialDenseMatrix & Elematrix2() { return elematrix2_; }
    Epetra_SerialDenseVector & Elevector1() { return elevector1_; }
    Epetra_SerialDenseVector & Elevector2() { return elevector2_; }
    Epetra_SerialDenseVector & Elevector3() { return elevector3_; }
    //@}

    /// zero element memory
    void ClearElementStorage( int eledim );

    //! @name Assembly Methods

    void AssembleMatrix1(int eid,
                         const std::vector<int>& lm,
                         const std::vector<int>& lmowner)
    {
      if ( Assemblemat1() )
      {
        Assemble( *systemmatrix1_, eid, elematrix1_, lm, lmowner );
      }
    }

    void AssembleMatrix2(int eid,
                         const std::vector<int>& lm,
                         const std::vector<int>& lmowner)
    {
      if ( Assemblemat2() )
      {
        Assemble( *systemmatrix2_, eid, elematrix2_, lm, lmowner );
      }
    }

    void AssembleVector1( const std::vector<int>& lm, const std::vector<int>& lmowner )
    {
      if ( Assemblevec1() )
      {
        Assemble( *systemvector1_, elevector1_, lm, lmowner );
      }
    }

    void AssembleVector2( const std::vector<int>& lm, const std::vector<int>& lmowner )
    {
      if ( Assemblevec2() )
      {
        Assemble( *systemvector2_, elevector2_, lm, lmowner );
      }
    }

    void AssembleVector3( const std::vector<int>& lm, const std::vector<int>& lmowner )
    {
      if ( Assemblevec3() )
      {
        Assemble( *systemvector3_, elevector3_, lm, lmowner );
      }
    }

    //@}

    //! @name Assembly Methods

    virtual void Assemble(LINALG::SparseOperator& sysmat,
                          int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner);

    virtual void Assemble(LINALG::SparseOperator& sysmat,
                          int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol);

    virtual void Assemble(LINALG::SparseOperator& sysmat, double val, int rgid, int cgid);

    virtual void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
                          const std::vector<int>& lm, const std::vector<int>& lmowner);

    virtual void Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
                          const std::vector<int>& lm, const std::vector<int>& lmowner);

    //@}

  private:

    //! @name Global Objects

    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1_;
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2_;
    Teuchos::RCP<Epetra_Vector>          systemvector1_;
    Teuchos::RCP<Epetra_Vector>          systemvector2_;
    Teuchos::RCP<Epetra_Vector>          systemvector3_;

    //@}

    //! @name Element Local Objects
    /// define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1_;
    Epetra_SerialDenseMatrix elematrix2_;
    Epetra_SerialDenseVector elevector1_;
    Epetra_SerialDenseVector elevector2_;
    Epetra_SerialDenseVector elevector3_;

    //@}
  };

}

#endif
