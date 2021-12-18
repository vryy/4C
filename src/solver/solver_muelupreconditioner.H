/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef SOLVER_MUELUPRECONDITIONER_H_
#define SOLVER_MUELUPRECONDITIONER_H_

#include "solver_preconditionertype.H"

// Xpetra
#include <Xpetra_MultiVector.hpp>

// MueLu
#include <MueLu_Hierarchy.hpp>

// Header files for default types must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>

// Baci
#include "../linalg/linalg_blocksparsematrix.H"

namespace LINALG
{
  namespace SOLVER
  {
    /*! \brief General layout of the MueLu preconditioner working with sparse matrices
     * in CRS format.
     *
     * For more specialized versions (e.g. block systems) derive from this
     * class as it already owns all core components of a MueLu preconditioner.
     */
    class MueLuPreconditioner : public PreconditionerType
    {
     public:
      MueLuPreconditioner(FILE* outfile, Teuchos::ParameterList& muelulist);

      /*! \brief Create and compute the preconditioner
       *
       * The MueLu preconditioner only works for matrices of the type
       * Epetra_CrsMatrix. We check whether the input matrix is of proper type
       * and throw an error if not!
       *
       * This routine either re-creates the entire preconditioner from scratch or
       * it re-uses the existing preconditioner and only updates the fine level matrix
       * for the Krylov solver.
       *
       * It maintains backward compability to the ML interface!
       *
       * @param create Boolean flag to enforce (re-)creation of the preconditioner
       * @param matrix Epetra_CrsMatrix to be used as input for the preconditioner
       * @param x Solution of the linear system
       * @param b Right-hand side of the linear system
       */
      void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      //! Get the actual operator used for preconditioning
      Epetra_Operator* PrecOperator() const final { return &*P_; }

      //! linear operator used for preconditioning
      Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const final { return P_; }

      void Print(std::ostream& stream) { stream << "MueLuPreconditioner"; }

      //! return name of sublist in paramterlist which contains parameters for preconditioner
      const std::string getParameterListName() const { return "MueLu Parameters"; }

     private:
      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix_;

     protected:
      //! MueLu parameter list
      Teuchos::ParameterList& muelulist_;

      //! preconditioner
      Teuchos::RCP<Epetra_Operator> P_;

      //! MueLu hierarchy
      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;

    };  // class MueLuPreconditioner


    /*! \brief Specialized MueLu preconditioner for blocked fluid systems.
     */
    class MueLuFluidBlockPreconditioner : public MueLuPreconditioner
    {
     public:
      MueLuFluidBlockPreconditioner(FILE* outfile, Teuchos::ParameterList& muelulist);

      /*! \brief Create and compute the preconditioner
       *
       * The Muelu block-preconditioner only works for matrices of the type
       * Epetra_BlockCrsMatrix. We check whether the input matrix is of proper type
       * and throw an error if not!
       *
       * @param create Boolean flag to enforce (re-)creation of the preconditioner
       * @param matrix BlockSparseMatrixBase to be used as input for the preconditioner
       * @param x Solution of the linear system
       * @param b Right-hand side of the linear system
       */
      void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
          Epetra_MultiVector* b) override;

      void Print(std::ostream& stream) override { stream << "MueLuFlBlockPreconditioner"; }

      //! return name of sublist in paramterlist which contains parameters for preconditioner
      const std::string getParameterListName() const override { return "MueLu (Fluid) Parameters"; }

     private:
      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<BlockSparseMatrixBase> Pmatrix_;

    };  // class MueLuFlBlockPreconditioner


    /*! \brief Specialized MueLu preconditioner for blocked tsi systems.
     */
    class MueLuTsiBlockPreconditioner : public MueLuPreconditioner
    {
     public:
      MueLuTsiBlockPreconditioner(FILE* outfile, Teuchos::ParameterList& muelulist);

      /*! \brief Create and compute the preconditioner
       *
       * The Muelu block-preconditioner only works for matrices of the type
       * Epetra_BlockCrsMatrix. We check whether the input matrix is of proper type
       * and throw an error if not!
       *
       * @param create Boolean flag to enforce (re-)creation of the preconditioner
       * @param matrix BlockSparseMatrixBase to be used as input for the preconditioner
       * @param x Solution of the linear system
       * @param b Right-hand side of the linear system
       */
      void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
          Epetra_MultiVector* b) override;

      void Print(std::ostream& stream) override { stream << "MueLuTsiBlockPreconditioner"; }

      // return name of sublist in paramterlist which contains parameters for preconditioner
      const std::string getParameterListName() const override { return "MueLu (TSI) Parameters"; }

     private:
      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<BlockSparseMatrixBase> Pmatrix_;

    };  // class MueLuTsiBlockPreconditioner


    /*! \brief MueLu preconditioner for blocked linear systems of equations for contact problems
     * in saddlepoint formulation.
     */
    class MueLuContactSpPreconditioner : public MueLuPreconditioner
    {
     public:
      MueLuContactSpPreconditioner(FILE* outfile, Teuchos::ParameterList& muelulist);

      /*! \brief Create and compute the preconditioner
       *
       * The saddle-point preconditioner only works for matrices of the type
       * BlockSparseMatrix. We assume the structure and block indices to be
       * \f[
       * A = \left[\begin{array}{cc}
       *       A_{11} & A_{12}\\
       *       A_{21} & A_{22}
       *     \end{array}\right]
       * \f]
       * We check whether the input \c matrix is of proper type and throw an
       * error if not.
       *
       * This routine either re-create the entire preconditioner from scratch or
       * it re-uses the existing preconditioner and only updates the fine level matrix
       * for the Krylov solver.
       *
       * @param create Boolean flag to enforce (re-)creation of the preconditioner
       * @param matrix BlockSparseMatrix to be used as input for the preconditioner
       * @param x Solution of the linear system
       * @param b Right-hand side of the linear system
       */
      void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
          Epetra_MultiVector* b) override;

      void Print(std::ostream& stream) override { stream << "MueLuContactSpPreconditioner"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      const std::string getParameterListName() const override
      {
        return "MueLu (Contact) Parameters";
      }

     private:
      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<BlockSparseMatrixBase> Pmatrix_;

    };  // class MueLuContactSpPreconditioner


    namespace MUELU
    {
      namespace UTILS
      {
        /*!
         * \brief Extract nullspace from Baci-provided ML list and convert to Xpetra::MultiVector
         *
         * \pre The input parameter list \c mllist needs to contain these entries:
         *   - "null space: dimension" (type: \c int )
         *   - "nullspace" (type: \c RCP<std::vector<double>> )
         *
         * @param[in] rowMap Map to be used to create the null space vectors
         * @param[in] mllist Parameter list, where Baci has stored the null space data in ML style
         * @return MultiVector with nullspace data
         */
        Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
        ExtractNullspaceFromMLList(
            const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
            Teuchos::ParameterList& muelulist);

        /*!
         * \brief Sets the striding maps of a matrix
         *
         * \note this function is copied from Xpetra and can be
         * removed after the next trilinos update
         *
         * @param matrix Matrix with striding maps to be set
         * @param rangeStridingInfo Range map striding information
         * @param domainStridingInfo Domain map striding information
         */
        void convertMatrixToStridedMaps(
            Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> matrix,
            std::vector<size_t>& rangeStridingInfo, std::vector<size_t>& domainStridingInfo);

      }  // namespace UTILS
    }    // namespace MUELU
  }      // namespace SOLVER
}  // namespace LINALG

#endif /* SOLVER_MUELUPRECONDITIONER_H_ */
