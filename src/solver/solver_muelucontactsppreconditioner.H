/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 2

\maintainer Matthias Mayr
*/
/*-----------------------------------------------------------------------------*/

#ifndef SOLVER_MUELUCONTACTSPPRECONDITIONER_H_
#define SOLVER_MUELUCONTACTSPPRECONDITIONER_H_

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

#include "solver_preconditionertype.H"

#include <Xpetra_BlockedCrsMatrix_fwd.hpp>

#include <MueLu_Hierarchy.hpp>  // for MueLu
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>
#include <MueLu_SchurComplementFactory_fwd.hpp>


// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include "muelu/MueLu_ContactSPAggregationFactory_fwd.hpp"
#include "muelu/MueLu_MeshtyingSPAmalgamationFactory_fwd.hpp"
#include "muelu/MueLu_ContactSPRepartitionInterface_fwd.hpp"

// BACI includes
#include "../linalg/linalg_blocksparsematrix.H"

namespace LINALG
{
  namespace SOLVER
  {
    /*! \brief MueLu preconditioner for contact problems in saddle-point formulation

     Set of single-matrix algebraic multi-grid preconditioners for contact/meshtying problems
     */
    class MueLuContactSpPreconditioner : public PreconditionerType
    {
     public:
      // some typedefs
      typedef MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> FactoryManager;
      typedef MueLu::FactoryBase FactoryBase;
      typedef Scalar SC;
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;
      typedef Node NO;

      MueLuContactSpPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist);

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
      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      //! Get the actual operator used for preconditioning
      virtual Epetra_Operator* PrecOperator() const { return &*P_; }

      //! Get the actual operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print(std::ostream& stream) { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const
      {
        return "MueLu (Contact) Parameters";
      }  // TODO check me
     private:
      void Write(
          const Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& H);

      //! @name Create block smoothers
      //!@{

      /*! \brief Get block smoother factory on fine and intermediate levels
       *
       * The parameter list stores the smoother factory configuration. We extract this information
       * and then call subroutines to actually create a concrete instance of the smoother factory.
       *
       * @param paramList Parameter list with smoother configuration
       * @param level Current level ID
       * @return Instance of a concrete MueLu::SmootherFactory
       *
       * \sa GetCoarsestBlockSmootherFactory(), CreateBlockSmootherFactory()
       */
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetBlockSmootherFactory(const Teuchos::ParameterList& paramList, int level) const;

      /*! \brief Get block smoother factory on coarsest level
       *
       * The parameter list stores the smoother factory configuration. We extract this information
       * and then call subroutines to actually create a concrete instance of the smoother factory.
       *
       * @param paramList paramList Parameter list with smoother configuration
       * @return Instance of a concrete MueLu::SmootherFactory
       *
       * \sa GetBlockSmootherFactory(), CreateBlockSmootherFactory()
       */
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetCoarsestBlockSmootherFactory(const Teuchos::ParameterList& paramList) const;

      /*! \brief Create block smoother factory from parameter list
       *
       * The parameter list stores the smoother factory configuration. We extract this information
       * and then call subroutines to actually create a concrete instance of the smoother factory.
       *
       * @param type Type of block smoother to be created
       * @param paramList Parameter list with smoother configuration
       * @param bCoarse Boolean flag to indicate if this is the coarsest level
       * @return
       *
       * \sa GetBlockSmootherFactory(), GetCoarsestBlockSmootherFactory(),
       * GetUzawaSmootherFactory(), GetIndefBlockedDiagonalSmootherFactory(),
       * GetSIMPLESmootherFactory(), GetBraessSarazinSmootherFactory()
       */
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      CreateBlockSmootherFactory(const std::string& type, const Teuchos::ParameterList& paramList,
          const bool bCoarse) const;

      /*! \brief Interpret ML parameter lists for creating an Uzawa smoother
       *
       * @param paramList Parameter list with smoother configuration
       * @param bCoarse Boolean flag to indicate if this is the coarsest level
       * @return Instance of a Uzawa smoother
       *
       * \sa GetBlockSmootherFactory(), GetCoarsestBlockSmootherFactory()
       */
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetUzawaSmootherFactory(
          const Teuchos::ParameterList& paramList, const bool bCoarse = false) const;

      /*! \brief Interpret ML parameter lists for creating a indefinite block diagonal level
       * smoother
       *
       * @param paramList Parameter list with smoother configuration
       * @param bCoarse Boolean flag to indicate if this is the coarsest level
       * @return Instance of a indefinite block diagonal level smoother
       *
       * \sa GetBlockSmootherFactory(), GetCoarsestBlockSmootherFactory()
       */
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetIndefBlockedDiagonalSmootherFactory(
          const Teuchos::ParameterList& paramList, const bool bCoarse = false) const;

      /*! \brief Interpret ML parameter lists for creating a SIMPLE-type level smoother
       *
       * @param paramList Parameter list with smoother configuration
       * @param bCoarse Boolean flag to indicate if this is the coarsest level
       * @return Instance of a SIMPLE-type smoother
       *
       * \sa GetBlockSmootherFactory(), GetCoarsestBlockSmootherFactory()
       */
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetSIMPLESmootherFactory(
          const Teuchos::ParameterList& paramList, const bool bCoarse = false) const;

      /*! \brief Interprets ML parameter lists for creating a Braess Sarazin level smoother
       *
       * @param paramList Parameter list with smoother configuration
       * @param bCoarse Boolean flag to indicate if this is the coarsest level
       * @return Instance of a Braess-Sarazin smoother
       *
       * \sa GetBlockSmootherFactory(), GetCoarsestBlockSmootherFactory()
       */
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetBraessSarazinSmootherFactory(
          const Teuchos::ParameterList& paramList, const bool bCoarse = false) const;

      //!@}

      //! function interprets ML parameter lists for creating a Braess Sarazin level smoother
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      InterpretBACIList2MueLuSmoother(
          const Teuchos::ParameterList& paramList, Teuchos::RCP<FactoryBase> AFact) const;

      //! ML parameter list
      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<BlockSparseMatrixBase> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;

      /// MueLu hierarchy
      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;
    };  // class MueLuContactSpPreconditioner
  }     // namespace SOLVER
}  // namespace LINALG

#endif  // HAVE_MueLuContact

#endif /* SOLVER_MUELUCONTACTSPPRECONDITIONER_H_ */
