/*
 * solver_muelucontactpreconditioner3.H
 *
 *  Created on: Dec 6, 2012
 *      Author: tobias
 */

#ifndef SOLVER_MUELUCONTACTPRECONDITIONER3_H_
#define SOLVER_MUELUCONTACTPRECONDITIONER3_H_

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q3_2013

#include "solver_preconditionertype.H"


#include <MueLu_Hierarchy.hpp> // for MueLu
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_DirectSolver_fwd.hpp>
#include <MueLu_PermutationFactory_fwd.hpp>
#include <MueLu_AdaptiveSaMLParameterListInterpreter_fwd.hpp>
//#include <MueLu_PermutingSmoother_fwd.hpp>
#include <MueLu_MapTransferFactory_fwd.hpp>

// header files for default types, must be included after all other MueLu/Xpetra headers
#include <MueLu_UseDefaultTypes.hpp> // => Scalar=double, LocalOrdinal=GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>

namespace LINALG
{
  namespace SOLVER
  {
    /// ml preconditioners
    /*!
      Set of single-matrix algebraic multi-grid preconditioners
      for contact/meshtying problems
     */
    class MueLuContactPreconditioner3 : public PreconditionerType
    {
    public:

      MueLuContactPreconditioner3( FILE * outfile, Teuchos::ParameterList & mllist );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      virtual void Print( std::ostream & stream ) { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "MueLu (Contact3) Parameters"; }

    private:

      Teuchos::ParameterList & mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix>                 Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator>                  P_;

    };
  }
}

#endif // HAVE_Trilinos_Q1_2013
#endif // HAVE_MueLu


#endif /* SOLVER_MUELUCONTACTPRECONDITIONER3_H_ */
