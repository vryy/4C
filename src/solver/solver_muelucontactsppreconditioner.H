/*!-----------------------------------------------------------------------------
\file solver_muelucontactsppreconditioner.H

\brief Declaration

\level 2

\maintainer Matthias Mayr
*/
/*-----------------------------------------------------------------------------*/

#ifndef SOLVER_MUELUCONTACTSPPRECONDITIONER_H_
#define SOLVER_MUELUCONTACTSPPRECONDITIONER_H_

#ifdef HAVE_MueLu

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
       *  \warning The saddle-point preconditioner only works for matrices of the type
       * BlockSparseMatrix. We check whether the input \c matrix is of proper type and throw an
       * error if not.
       *
       * @param create Boolean flag to enforce creation of the preconditioner
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

      //! function interprets ML parameter lists
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetBlockSmootherFactory(const Teuchos::ParameterList& paramList, int level);

      //! function interprets coarse ML parameter list
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetCoarsestBlockSmootherFactory(const Teuchos::ParameterList& paramList);

      //! function interprets ML parameter lists for creating a indefinite block diagonal level
      //! smoother
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetUzawaSmootherFactory(const Teuchos::ParameterList& paramList, bool bCoarse = false);

      //! function interprets ML parameter lists for creating a indefinite block diagonal level
      //! smoother
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetIndefBlockedDiagonalSmootherFactory(
          const Teuchos::ParameterList& paramList, bool bCoarse = false);

      //! function interprets ML parameter lists for creating a SIMPLE level smoother
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetSIMPLESmootherFactory(const Teuchos::ParameterList& paramList, bool bCoarse = false);

      //! function interprets ML parameter lists for creating a Braess Sarazin level smoother
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetBraessSarazinSmootherFactory(
          const Teuchos::ParameterList& paramList, bool bCoarse = false);

      //! function interprets ML parameter lists for creating a Braess Sarazin level smoother
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      InterpretBACIList2MueLuSmoother(
          const Teuchos::ParameterList& paramList, Teuchos::RCP<FactoryBase> AFact);

      //! ML parameter list
      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<BlockSparseMatrixBase> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;

      /// MueLu hierarchy
      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H_;
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif  // HAVE_MueLu

#endif /* SOLVER_MUELUCONTACTSPPRECONDITIONER_H_ */
