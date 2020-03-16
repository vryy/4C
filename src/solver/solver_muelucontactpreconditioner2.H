/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 2

\maintainer Martin Kronbichler
*/
/*-----------------------------------------------------------------------------*/

#ifndef SOLVER_MUELUCONTACTPRECONDITIONER2_H_
#define SOLVER_MUELUCONTACTPRECONDITIONER2_H_

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

#include "solver_preconditionertype.H"


#include <MueLu_Hierarchy.hpp>  // for MueLu
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_DirectSolver_fwd.hpp>
#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_PermutingSmoother_fwd.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

#include "muelu/MueLu_ContactAFilterFactory_fwd.hpp"

namespace LINALG
{
  namespace SOLVER
  {
    /// ml preconditioners
    /*!
      Set of single-matrix algebraic multi-grid preconditioners
      for contact/meshtying problems
     */
    class MueLuContactPreconditioner2 : public PreconditionerType
    {
     public:
      // some typedefs
      typedef MueLu::FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> FactoryManager;
      typedef MueLu::FactoryBase FactoryBase;
      typedef Scalar SC;
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;
      typedef Node NO;

      MueLuContactPreconditioner2(FILE* outfile, Teuchos::ParameterList& mllist);

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual Epetra_Operator* PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print(std::ostream& stream) { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const
      {
        return "MueLu (Contact2) Parameters";
      }

     private:
      //! function interprets ML parameter lists for creating contact level smoothers
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetContactSmootherFactory(const Teuchos::ParameterList& paramList, int level,
          const Teuchos::RCP<FactoryBase>& AFact);

      //! function interprets ML parameter lists for creating contact coarse solver
      Teuchos::RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
      GetContactCoarsestSolverFactory(const Teuchos::ParameterList& paramList,
          const Teuchos::RCP<FactoryBase>& AFact = Teuchos::null);

      // void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* AFact,
      // const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm, int
      // numDofsPerNode);

      // private Initialization routine for MueLu Multigrid hierarchy
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> InitializeHierarchy(
          const Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>& h,
          const Teuchos::ParameterList& params,
          const Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>>& A,
          const Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nsp = Teuchos::null);

      // private Setup routine for MueLu Multigrid factories
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> SetupFactories(
          const Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>& h,
          const Teuchos::ParameterList& params);

      // private Setup routine for smoothers
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> SetupSmoothers(
          const Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>& h,
          const Teuchos::ParameterList& params);

      // private Setup routine for MueLu Multigrid hierarchy
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> SetupHierarchy(
          const Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>>& h,
          const Teuchos::ParameterList& params,
          const Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>>& A,
          const Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nsp = Teuchos::null);

      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;

      /// hierarchy object
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H_;

      /// multigrid levels
      Teuchos::Array<Teuchos::RCP<FactoryManager>> vecManager_;

      /// A Factory which shall be reused by MyTrilinosSmoother
      Teuchos::RCP<MueLu::Factory> singleNodeAFact_;
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif  // HAVE_MueLuContact

#endif /* SOLVER_MUELUCONTACTPRECONDITIONER2_H_ */
