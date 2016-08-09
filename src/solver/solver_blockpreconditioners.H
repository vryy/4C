/*!----------------------------------------------------------------------
\file solver_blockpreconditioners.H

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*----------------------------------------------------------------------*/
#ifndef SOLVER_BLOCKPRECONDITIONERS_H_
#define SOLVER_BLOCKPRECONDITIONERS_H_

#include "solver_preconditionertype.H"

namespace LINALG
{
  namespace SOLVER
  {
    /// MueLu block preconditioner
    /*!
     * Multigrid block preconditioner based on MueLu
     * Block preconditioners assume the Epetra_Operator to be a
     * LINALG::BlockSparseMatrix.
     */
    class MueLuBlockPreconditioner : public PreconditionerType
    {
    public:

      MueLuBlockPreconditioner( FILE * outfile,
          Teuchos::ParameterList & params );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print( std::ostream & stream ) { stream << "MueLuBlockPreconditioner"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      /*! This function is not needed since we use the XML interface */
      virtual const std::string getParameterListName() const { return ""; }

    private:
      Teuchos::ParameterList & params_;

      Teuchos::RCP<Epetra_Operator>    P_;
    };

    /// SIMPLE(R) block preconditioner
    /*!
      Block preconditioners assume the Epetra_Operator to be a
      LINALG::BlockSparseMatrix.
     */
    class SimplePreconditioner : public PreconditionerType
    {
    public:

      SimplePreconditioner( FILE * outfile,
          Teuchos::ParameterList & params );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print( std::ostream & stream ) { stream << "CheapSimple"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "CheapSIMPLE Parameters"; }

    private:
      Teuchos::ParameterList & params_;

      Teuchos::RCP<Epetra_Operator>    P_;
    };

    /// General purpose block gauss-seidel preconditioner
    /*!
      2x2 block preconditioner
     */
    class BGSPreconditioner : public PreconditionerType
    {
    public:

      BGSPreconditioner( FILE * outfile,
          Teuchos::ParameterList & params, Teuchos::ParameterList & bgslist );

      virtual void Setup( bool create,
                          Epetra_Operator * matrix,
                          Epetra_MultiVector * x,
                          Epetra_MultiVector * b );

      virtual Epetra_Operator * PrecOperator() const { return &*P_; }

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const { return P_; }

      virtual void Print( std::ostream & stream ) { stream << "BGS"; }

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      virtual const std::string getParameterListName() const { return "BGS Parameters"; }

    private:

      Teuchos::ParameterList & params_;
      Teuchos::ParameterList & bgslist_;
      Teuchos::RCP<Epetra_Operator>    P_;
    };
  }
}


#endif /* SOLVER_BLOCKPRECONDITIONERS_H_ */
