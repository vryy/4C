/*----------------------------------------------------------------------*/
/*!
\file drt_assemblestrategy.H

\brief Setup of the list of valid conditions for input

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/


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
  /*!
    The global matrices and vectors that are filled during the assembly
    process. This class manages the element assembling. This simplifies the
    element loop in DRT::Discretization::Evaluate(). Furthermore, the strategy
    can be exchanged in the assembly process needs to be modified.

    \author u.kue
   */
  class AssembleStrategy
  {
   public:
    /// Construct with allocated global objects or Teuchos::null.
    AssembleStrategy(int firstdofset, int seconddofset,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3);

    /// Destruct
    virtual ~AssembleStrategy();

    int FirstDofSet() const { return firstdofset_; }

    int SecondDofSet() const { return seconddofset_; }

    /// Create matrix graph to create a filled matrix
    /*!
      The graph will already include the Dirichlet lines. That is there will
      be just the diagonal entry on each line that is a Dirichlet line.
      LINALG::SparseMatrix can assemble to filled matrices with ready made
      Dirichlet lines. Thus the idea is to explicitly create a graph and
      creating filled matrices from it. It will never be required to assemble
      to unfilled matrices or reallocate any matrices as long as the mesh
      remains the same.

      \note It is much more expensive to assemble to an unfilled matrix than a
      filled one.
     */
    virtual Teuchos::RCP<Epetra_CrsGraph> MatrixGraph(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Map> dbcmap);

    //! @name Assembly Flags
    /// Tell which matrix and vector is available to be assembled.

    bool Assemblemat1() { return systemmatrix1_ != Teuchos::null; }
    bool Assemblemat2() { return systemmatrix2_ != Teuchos::null; }
    bool Assemblevec1() { return systemvector1_ != Teuchos::null; }
    bool Assemblevec2() { return systemvector2_ != Teuchos::null; }
    bool Assemblevec3() { return systemvector3_ != Teuchos::null; }
    //@}

    //! @name Access Methods to Global Object
    Teuchos::RCP<LINALG::SparseOperator> Systemmatrix1() { return systemmatrix1_; }
    Teuchos::RCP<LINALG::SparseOperator> Systemmatrix2() { return systemmatrix2_; }
    Teuchos::RCP<Epetra_Vector> Systemvector1() { return systemvector1_; }
    Teuchos::RCP<Epetra_Vector> Systemvector2() { return systemvector2_; }
    Teuchos::RCP<Epetra_Vector> Systemvector3() { return systemvector3_; }
    //@}

    //! @name Access Methods to Element Local Object
    Epetra_SerialDenseMatrix& Elematrix1() { return elematrix1_; }
    Epetra_SerialDenseMatrix& Elematrix2() { return elematrix2_; }
    Epetra_SerialDenseVector& Elevector1() { return elevector1_; }
    Epetra_SerialDenseVector& Elevector2() { return elevector2_; }
    Epetra_SerialDenseVector& Elevector3() { return elevector3_; }
    //@}

    /// zero global storage
    void Zero();

    /// complete any global objects
    virtual void Complete();

    /// zero element memory
    void ClearElementStorage(int rdim, int cdim);

    //! @name Specific Assembly Methods

    /// Asseble to matrix 1
    void AssembleMatrix1(int eid, const std::vector<int>& rlm, const std::vector<int>& clm,
        const std::vector<int>& lmowner, const std::vector<int>& lmstride)
    {
      if (Assemblemat1())
      {
        Assemble(*systemmatrix1_, eid, lmstride, elematrix1_, rlm, lmowner, clm);
      }
    }

    /// Asseble to matrix 2
    void AssembleMatrix2(int eid, const std::vector<int>& rlm, const std::vector<int>& clm,
        const std::vector<int>& lmowner, const std::vector<int>& lmstride)
    {
      if (Assemblemat2())
      {
        Assemble(*systemmatrix2_, eid, lmstride, elematrix2_, rlm, lmowner, clm);
      }
    }

    /// Asseble to vector 1
    void AssembleVector1(const std::vector<int>& lm, const std::vector<int>& lmowner)
    {
      if (Assemblevec1())
      {
        Assemble(*systemvector1_, elevector1_, lm, lmowner);
      }
    }

    /// Asseble to vector 2
    void AssembleVector2(const std::vector<int>& lm, const std::vector<int>& lmowner)
    {
      if (Assemblevec2())
      {
        Assemble(*systemvector2_, elevector2_, lm, lmowner);
      }
    }

    /// Asseble to vector 3
    void AssembleVector3(const std::vector<int>& lm, const std::vector<int>& lmowner)
    {
      if (Assemblevec3())
      {
        Assemble(*systemvector3_, elevector3_, lm, lmowner);
      }
    }

    //@}

    //! @name General Purpose Assembly Methods

    /// Assemble to given matrix using nodal stride
    virtual void Assemble(LINALG::SparseOperator& sysmat, int eid, const std::vector<int>& lmstride,
        const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lm,
        const std::vector<int>& lmowner);

    /// Assemble to given matrix
    virtual void Assemble(LINALG::SparseOperator& sysmat, int eid, const std::vector<int>& lmstride,
        const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
        const std::vector<int>& lmrowowner, const std::vector<int>& lmcol);

    /// Assemble to given matrix
    virtual void Assemble(LINALG::SparseOperator& sysmat, double val, int rgid, int cgid);

    /// Assemble to given vector
    virtual void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
        const std::vector<int>& lm, const std::vector<int>& lmowner);

    /// Assemble to given vector
    virtual void Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
        const std::vector<int>& lm, const std::vector<int>& lmowner);

    //@}

   private:
    int firstdofset_;
    int seconddofset_;

    //! @name Global Objects

    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1_;
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2_;
    Teuchos::RCP<Epetra_Vector> systemvector1_;
    Teuchos::RCP<Epetra_Vector> systemvector2_;
    Teuchos::RCP<Epetra_Vector> systemvector3_;

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

}  // namespace DRT

#endif
