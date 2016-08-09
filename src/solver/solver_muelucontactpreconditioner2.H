/*!-----------------------------------------------------------------------------
\file solver_muelucontactpreconditioner2.H

\brief Declaration

\level 2

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------------*/

#ifndef SOLVER_MUELUCONTACTPRECONDITIONER2_H_
#define SOLVER_MUELUCONTACTPRECONDITIONER2_H_

#ifdef HAVE_MueLu

#include "solver_preconditionertype.H"


#include <MueLu_Hierarchy.hpp> // for MueLu
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_DirectSolver_fwd.hpp>
#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_PermutingSmoother_fwd.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>

#include "muelu/muelu_ContactAFilterFactory_fwd.hpp"

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

      MueLuContactPreconditioner2( FILE * outfile, Teuchos::ParameterList & mllist );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print( std::ostream & stream ) { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "MueLu (Contact2) Parameters"; }
    private:

      //! function interprets ML parameter lists for creating contact level smoothers
      Teuchos::RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetContactSmootherFactory(const Teuchos::ParameterList & paramList, int level, const Teuchos::RCP<FactoryBase> & AFact);

      //! function interprets ML parameter lists for creating contact coarse solver
      Teuchos::RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node> > GetContactCoarsestSolverFactory(const Teuchos::ParameterList & paramList, const Teuchos::RCP<FactoryBase> & AFact = Teuchos::null);

      //void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* AFact, const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm, int numDofsPerNode);

      // private Initialization routine for MueLu Multigrid hierarchy
      Teuchos::RCP<Hierarchy> InitializeHierarchy(const Teuchos::RCP<Hierarchy> & h, const Teuchos::ParameterList & params, const Teuchos::RCP<Matrix> & A, const Teuchos::RCP<MultiVector> nsp = Teuchos::null);

      // private Setup routine for MueLu Multigrid factories
      Teuchos::RCP<Hierarchy> SetupFactories(const Teuchos::RCP<Hierarchy> & h, const Teuchos::ParameterList & params);

      // private Setup routine for smoothers
      Teuchos::RCP<Hierarchy> SetupSmoothers(const Teuchos::RCP<Hierarchy> & h, const Teuchos::ParameterList & params);

      // private Setup routine for MueLu Multigrid hierarchy
      Teuchos::RCP<Hierarchy> SetupHierarchy(const Teuchos::RCP<Hierarchy> & h, const Teuchos::ParameterList & params, const Teuchos::RCP<Matrix> & A, const Teuchos::RCP<MultiVector> nsp = Teuchos::null);

      Teuchos::ParameterList & mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix>                 Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator>                  P_;

      /// hierarchy object
      Teuchos::RCP<Hierarchy> H_;

      /// multigrid levels
      Teuchos::Array<Teuchos::RCP<FactoryManager> > vecManager_;

      /// A Factory which shall be reused by MyTrilinosSmoother
      Teuchos::RCP<Factory> singleNodeAFact_;
    };
  }
}

#endif // HAVE_MueLu

#endif /* SOLVER_MUELUCONTACTPRECONDITIONER2_H_ */
