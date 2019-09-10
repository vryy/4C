/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Dec 1, 2011
*/
/*----------------------------------------------------------------------*/

#ifndef SOLVER_MUELUPRECONDITIONER_H_
#define SOLVER_MUELUPRECONDITIONER_H_

#ifdef HAVE_MueLu

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// ml preconditioners
    /*!
      Set of single-matrix algebraic multi-grid preconditioners.
     */
    class MueLuPreconditioner : public PreconditionerType
    {
     public:
      MueLuPreconditioner(FILE* outfile, Teuchos::ParameterList& mllist);

      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      virtual Epetra_Operator* PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print(std::ostream& stream) { stream << "MueLu"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "MueLu Parameters"; }

     private:
      // some typedefs
      typedef Scalar SC;
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;
      typedef Node NO;

      void registerEpetraMapinMueLuLevel(
          Teuchos::ParameterList& paramList, std::string name, MueLu::Hierarchy<SC, LO, GO, NO>& H)
      {
        Teuchos::RCP<Epetra_Map> myMap = Teuchos::null;
        if (paramList.isSublist("Linear System properties"))
        {
          const Teuchos::ParameterList& linSystemProps =
              paramList.sublist("Linear System properties");
          // extract information provided by solver (e.g. PermutedAztecSolver)
          if (linSystemProps.isParameter(name))
            myMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>(name);
        }

        // check if myMap is an Epetra_Map
        Teuchos::RCP<Xpetra::EpetraMap> xmyMap = Teuchos::null;
        Teuchos::RCP<Epetra_Map> epTest = Teuchos::rcp_dynamic_cast<Epetra_Map>(myMap);
        if (epTest != Teuchos::null)
        {
          xmyMap = Teuchos::rcp(new Xpetra::EpetraMap(myMap));
          H.GetLevel(0)->Set(
              name, Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(xmyMap));
        }
        return;
      }

      void registerMapinMueLuLevel(
          Teuchos::ParameterList& paramList, std::string name, MueLu::Hierarchy<SC, LO, GO, NO>& H)
      {
        Teuchos::RCP<Xpetra::Map<LO, GO, NO>> myMap = Teuchos::null;
        if (paramList.isSublist("Linear System properties"))
        {
          const Teuchos::ParameterList& linSystemProps =
              paramList.sublist("Linear System properties");
          // extract information provided by solver (e.g. PermutedAztecSolver)
          if (linSystemProps.isParameter(name))
            myMap = linSystemProps.get<Teuchos::RCP<Xpetra::Map<LO, GO, NO>>>(name);
        }
        H.GetLevel(0)->Set(name, Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, Node>>(myMap));
        return;
      }

      void replaceAll(std::string& str, const std::string& from, const std::string& to);

      Teuchos::ParameterList& mllist_;

      //! system of equations used for preconditioning used by P_ only
      Teuchos::RCP<Epetra_CrsMatrix> Pmatrix_;

      /// preconditioner
      Teuchos::RCP<Epetra_Operator> P_;
    };
  }  // namespace SOLVER
}  // namespace LINALG

#endif  // HAVE_MueLu

#endif /* SOLVER_MUELUPRECONDITIONER_H_ */
