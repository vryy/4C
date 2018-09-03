/*!----------------------------------------------------------------------
\file solver_muelucontactpreconditioner_old.H

\brief Interface class for MueLu contact preconditioner

<pre>
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jan 11, 2012
</pre>
*----------------------------------------------------------------------*/

#ifndef SOLVER_MUELUCONTACTPRECONDITIONER_H_
#define SOLVER_MUELUCONTACTPRECONDITIONER_H_

#ifdef HAVE_MueLu

#include "solver_preconditionertype.H"


#include <MueLu_Hierarchy.hpp>  // for MueLu

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

/*class HHierarchy;
class OOperator;
class MMultiVector;*/

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

     private:
      // void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* AFact,
      // const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm, int
      // numDofsPerNode);

      // private Setup routine for MueLu Multigrid hierarchy
      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> SetupHierarchy(
          const Teuchos::ParameterList& params,
          const Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
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
