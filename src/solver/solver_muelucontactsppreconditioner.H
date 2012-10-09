/*
 * solver_muelucontactsppreconditioner.H
 *
 *  Created on: Sep 23, 2012
 *      Author: tobias
 */

#ifndef SOLVER_MUELUCONTACTSPPRECONDITIONER_H_
#define SOLVER_MUELUCONTACTSPPRECONDITIONER_H_

#ifdef HAVE_MueLu

#include "solver_preconditionertype.H"


#include <Xpetra_BlockedCrsMatrix_fwd.hpp>

#include <MueLu_Hierarchy.hpp> // for MueLu
#include <MueLu_SmootherFactory.hpp>



// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>



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
    class MueLuContactSpPreconditioner : public PreconditionerType
    {
    public:

      MueLuContactSpPreconditioner( FILE * outfile, Teuchos::ParameterList & mllist );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      virtual void Print( std::ostream & stream ) { stream << "MueLu"; }

    private:

      //! function interprets ML parameter lists for creating contact level smoothers
      //Teuchos::RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > GetContactSmootherFactory(const Teuchos::ParameterList & paramList, int level, const Teuchos::RCP<FactoryBase> & AFact);

      //void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* AFact, const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm, int numDofsPerNode);

      // private Setup routine for MueLu Multigrid hierarchy
      Teuchos::RCP<Hierarchy> SetupHierarchy(const Teuchos::ParameterList & params, const Teuchos::RCP<Matrix> & A, const Teuchos::RCP<MultiVector> nsp = Teuchos::null);

      Teuchos::ParameterList & mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix>                 Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator>                  P_;

    };
  }
}

#endif // HAVE_MueLu


#endif /* SOLVER_MUELUCONTACTSPPRECONDITIONER_H_ */
