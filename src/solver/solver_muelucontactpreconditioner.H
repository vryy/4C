/*!----------------------------------------------------------------------

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jan 11, 2012
*/
/*----------------------------------------------------------------------*/

#ifndef SOLVER_MUELUCONTACTPRECONDITIONER_H_
#define SOLVER_MUELUCONTACTPRECONDITIONER_H_

#ifdef HAVE_MueLu

#include "solver_preconditionertype.H"


#include <MueLu_Hierarchy.hpp>  // for MueLu
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_DirectSolver_fwd.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

namespace LINALG
{
  namespace SOLVER
  {
    /// ml preconditioners
    /*!
      Set of single-matrix algebraic multi-grid preconditioners
      for contact/meshtying problems
     */
    class MueLuContactPreconditioner : public PreconditionerType
    {
     public:
      MueLuContactPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist);

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual Epetra_Operator* PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print(std::ostream& stream) { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const
      {
        return "MueLu (Contact) Parameters";
      }

     private:
      //! function interprets ML parameter lists for creating contact level smoothers
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetContactSmootherFactory(const Teuchos::ParameterList& paramList, int level,
          const Teuchos::RCP<MueLu::FactoryBase>& AFact);

      //! function interprets ML parameter lists for creating contact coarse solver
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetContactCoarsestSolverFactory(const Teuchos::ParameterList& paramList,
          const Teuchos::RCP<MueLu::FactoryBase>& AFact = Teuchos::null);

      // private Setup routine for MueLu Multigrid hierarchy
      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> SetupHierarchy(
          const Teuchos::ParameterList& params,
          const Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
          const Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nsp =
              Teuchos::null);

      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif  // HAVE_MueLu


#endif /* SOLVER_MUELUCONTACTPRECONDITIONER_H_ */
