/*----------------------------------------------------------------------*/
/*! \file

\brief Routines for handing a collection of element matrices and vectors to the actual assembly
calls into one global sparse matrix and global load vector

\level 0


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_FEM_GENERAL_ASSEMBLESTRATEGY_HPP
#define FOUR_C_FEM_GENERAL_ASSEMBLESTRATEGY_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SparseOperator;
}

namespace Core::FE
{
  /// General assembly strategy
  /*!
    The global matrices and vectors that are filled during the assembly
    process. This class manages the element assembling. This simplifies the
    element loop in Core::FE::Discretization::evaluate(). Furthermore, the strategy
    can be exchanged in the assembly process needs to be modified.

    \author u.kue
   */
  class AssembleStrategy
  {
   public:
    /// Construct with allocated global objects or Teuchos::null.
    AssembleStrategy(int firstdofset, int seconddofset,
        Teuchos::RCP<LinAlg::SparseOperator> systemmatrix1,
        Teuchos::RCP<LinAlg::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3);

    /// Destruct
    virtual ~AssembleStrategy() = default;

    int first_dof_set() const { return firstdofset_; }

    int second_dof_set() const { return seconddofset_; }

    //! @name Assembly Flags
    /// Tell which matrix and vector is available to be assembled.

    bool assemblemat1() { return systemmatrix1_ != Teuchos::null; }
    bool assemblemat2() { return systemmatrix2_ != Teuchos::null; }
    bool assemblevec1() { return systemvector1_ != Teuchos::null; }
    bool assemblevec2() { return systemvector2_ != Teuchos::null; }
    bool assemblevec3() { return systemvector3_ != Teuchos::null; }
    //@}

    //! @name Access Methods to Global Object
    Teuchos::RCP<LinAlg::SparseOperator> systemmatrix1() { return systemmatrix1_; }
    Teuchos::RCP<LinAlg::SparseOperator> systemmatrix2() { return systemmatrix2_; }
    Teuchos::RCP<Epetra_Vector> systemvector1() { return systemvector1_; }
    Teuchos::RCP<Epetra_Vector> systemvector2() { return systemvector2_; }
    Teuchos::RCP<Epetra_Vector> systemvector3() { return systemvector3_; }
    //@}

    //! @name Access Methods to Element Local Object
    LinAlg::SerialDenseMatrix& elematrix1() { return elematrix1_; }
    LinAlg::SerialDenseMatrix& elematrix2() { return elematrix2_; }
    LinAlg::SerialDenseVector& elevector1() { return elevector1_; }
    LinAlg::SerialDenseVector& elevector2() { return elevector2_; }
    LinAlg::SerialDenseVector& elevector3() { return elevector3_; }
    //@}

    /// zero global storage
    void zero();

    /// complete any global objects
    virtual void complete();

    /// zero element memory
    void clear_element_storage(int rdim, int cdim);

    //! @name Specific Assembly Methods

    /// Asseble to matrix 1
    void assemble_matrix1(int eid, const std::vector<int>& rlm, const std::vector<int>& clm,
        const std::vector<int>& lmowner, const std::vector<int>& lmstride)
    {
      if (assemblemat1())
      {
        assemble(*systemmatrix1_, eid, lmstride, elematrix1_, rlm, lmowner, clm);
      }
    }

    /// Asseble to matrix 2
    void assemble_matrix2(int eid, const std::vector<int>& rlm, const std::vector<int>& clm,
        const std::vector<int>& lmowner, const std::vector<int>& lmstride)
    {
      if (assemblemat2())
      {
        assemble(*systemmatrix2_, eid, lmstride, elematrix2_, rlm, lmowner, clm);
      }
    }

    /// Asseble to vector 1
    void assemble_vector1(const std::vector<int>& lm, const std::vector<int>& lmowner)
    {
      if (assemblevec1())
      {
        assemble(*systemvector1_, elevector1_, lm, lmowner);
      }
    }

    /// Asseble to vector 2
    void assemble_vector2(const std::vector<int>& lm, const std::vector<int>& lmowner)
    {
      if (assemblevec2())
      {
        assemble(*systemvector2_, elevector2_, lm, lmowner);
      }
    }

    /// Asseble to vector 3
    void assemble_vector3(const std::vector<int>& lm, const std::vector<int>& lmowner)
    {
      if (assemblevec3())
      {
        assemble(*systemvector3_, elevector3_, lm, lmowner);
      }
    }

    //@}

    //! @name General Purpose Assembly Methods

    /// Assemble to given matrix using nodal stride
    virtual void assemble(LinAlg::SparseOperator& sysmat, int eid, const std::vector<int>& lmstride,
        const LinAlg::SerialDenseMatrix& Aele, const std::vector<int>& lm,
        const std::vector<int>& lmowner);

    /// Assemble to given matrix
    virtual void assemble(LinAlg::SparseOperator& sysmat, int eid, const std::vector<int>& lmstride,
        const LinAlg::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
        const std::vector<int>& lmrowowner, const std::vector<int>& lmcol);

    /// Assemble to given matrix
    virtual void assemble(LinAlg::SparseOperator& sysmat, double val, int rgid, int cgid);

    /// Assemble to given vector
    virtual void assemble(Epetra_Vector& V, const LinAlg::SerialDenseVector& Vele,
        const std::vector<int>& lm, const std::vector<int>& lmowner);

    /// Assemble to given vector
    virtual void assemble(Epetra_MultiVector& V, const int n, const LinAlg::SerialDenseVector& Vele,
        const std::vector<int>& lm, const std::vector<int>& lmowner);

    //@}

   private:
    int firstdofset_;
    int seconddofset_;

    //! @name Global Objects

    Teuchos::RCP<LinAlg::SparseOperator> systemmatrix1_;
    Teuchos::RCP<LinAlg::SparseOperator> systemmatrix2_;
    Teuchos::RCP<Epetra_Vector> systemvector1_;
    Teuchos::RCP<Epetra_Vector> systemvector2_;
    Teuchos::RCP<Epetra_Vector> systemvector3_;

    //@}

    //! @name Element Local Objects
    /// define element matrices and vectors
    LinAlg::SerialDenseMatrix elematrix1_;
    LinAlg::SerialDenseMatrix elematrix2_;
    LinAlg::SerialDenseVector elevector1_;
    LinAlg::SerialDenseVector elevector2_;
    LinAlg::SerialDenseVector elevector3_;

    //@}
  };

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
