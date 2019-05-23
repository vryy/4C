/*!----------------------------------------------------------------------

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_CHEAPSIMPLE_H
#define SOLVER_CHEAPSIMPLE_H

#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_mapextractor.H"  // This is for mmex_ (reference, not Teuchos::RCP by intention)

using Teuchos::ParameterList;  // note, we do not include a whole namespace!

namespace LINALG
{
  ////////////////////////////////////
  // forward declarations
  class Solver;
  class SparseMatrixBase;
  class MultiMapExtractor;
  class DownwindMatrix;
  namespace ANA
  {
    class Vector;
  }
  ////////////////////////////////////

  namespace SOLVER
  {
    /*!
    \brief A Semi-implicit Method for Pressure Linked Equations (SIMPLE)
           preconditioner

    This Operator implements the family of SIMPLE methods such as
    SIMPLE, SIMPLER, SIMPLEC and CheapSIMPLE (a m.gee's variation of SIMPLE)

    Literature:<br>

    Elman, H., Howle, V.E., Shadid, J., Shuttleworth, R., Tuminaro, R.:
    A taxonomy and comparison of parallel block multi-level
    preconditioners for the incomp. Navier-Stokes equations.
    Sandia technical report SAND2007-2761, 2007,
    Also appeared in JCP

    Pernice, M., Tocci, M.D.:
    A Multigrid Preconditioned Newton-Krylov method for the incomp.
    Navier-Stokes equations, Siam, J. Sci. Comp. 23, pp. 398-418 (2001)

    \author gee (gee@lnm.mw.tum.de)
     */
    class CheapSIMPLE_BlockPreconditioner : public virtual Epetra_Operator
    {
     public:
      /*!
      \brief Standard Constructor
       */
      explicit CheapSIMPLE_BlockPreconditioner(Teuchos::RCP<Epetra_Operator> A,
          const Teuchos::ParameterList& velocitylist, const Teuchos::ParameterList& pressurelist,
          FILE* outfile = NULL);

      /*!
      \brief Destructor
       */
      virtual ~CheapSIMPLE_BlockPreconditioner()
      {
        vsolver_ = Teuchos::null;  // solver objects only generated for SIMPLE and SIMPLER (not for
                                   // CheapSIMPLE)
        psolver_ = Teuchos::null;
        Ppredict_ = Teuchos::null;
        Pschur_ = Teuchos::null;
      }

      /*!
      \brief Setup the label of this class.
       */
      const std::string SetupLabel()
      {
        std::stringstream strLabel;

        strLabel << "CheapSIMPLE ";

        const bool visml = predictSolver_list_.isSublist("ML Parameters");
        const bool vismuelu = predictSolver_list_.isSublist("MueLu Parameters");
        const bool visifpack = predictSolver_list_.isSublist("IFPACK Parameters");
        strLabel << "[Inv1=";
        if (visml)
          strLabel << "ML-type";
        else if (vismuelu)
          strLabel << "MueLu-type";
        else if (visifpack)
        {
          const Teuchos::ParameterList& ifpackParams =
              predictSolver_list_.sublist("IFPACK Parameters");
          if (ifpackParams.isParameter("relaxation: type"))
          {
            strLabel << ifpackParams.get<std::string>("relaxation: type");
            strLabel << " (" << ifpackParams.get<int>("relaxation: sweeps") << ", "
                     << ifpackParams.get<double>("relaxation: damping factor") << ")";
          }
          else
            strLabel << "ILU";
        }
        else
          strLabel << "unknown";

        strLabel << ",Inv2=";
        const bool pisml = schurSolver_list_.isSublist("ML Parameters");
        const bool pismuelu = predictSolver_list_.isSublist("MueLu Parameters");
        const bool pisifpack = schurSolver_list_.isSublist("IFPACK Parameters");
        if (pisml)
          strLabel << "ML-type";
        else if (pismuelu)
          strLabel << "MueLu-type";
        else if (pisifpack)
        {
          const Teuchos::ParameterList& ifpackParams =
              schurSolver_list_.sublist("IFPACK Parameters");
          if (ifpackParams.isParameter("relaxation: type"))
          {
            strLabel << ifpackParams.get<std::string>("relaxation: type");
            strLabel << " (" << ifpackParams.get<int>("relaxation: sweeps") << ", "
                     << ifpackParams.get<double>("relaxation: damping factor") << ")";
          }
          else
            strLabel << "ILU";
        }
        else
          strLabel << "unknown";

        strLabel << "]";

        return strLabel.str();
      }
      /*
       *  \brief Returns the label of this class.
       */
      const char* Label() const { return &label_[0]; }

      /*!
      \brief get Comm of this class

      Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class

       */
      const Epetra_Comm& Comm() const { return (A_->Comm()); }


      /*!
      \brief Get fine level OperatorDomainMap

      Derived from Epetra_Operator, get fine level OperatorDomainMap

       */
      const Epetra_Map& OperatorDomainMap() const { return A_->FullDomainMap(); }

      /*!
      \brief Get fine level OperatorRangeMap

      Derived from Epetra_Operator, get fine level OperatorRangeMap

       */
      const Epetra_Map& OperatorRangeMap() const { return A_->FullRangeMap(); }

      /*!
      \brief ApplyInverse the preconditioner

       ApplyInverse the preconditioner. Method is derived from Epetra_Operator.


       \param X   (In) : Epetra_MultiVector matching the fine level map of this
                         preconditioner
       \param Y (Out)  : Epetra_MultiVector containing the result on output
       */
      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      /*!
      \brief not implemented
       */
      int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        dserror("Apply does not make sense for LINALG::SIMPLER_Operator");
        return (-1);
      }

      /*!
      \brief not implemented
       */
      int SetUseTranspose(bool UseTranspose)
      {
        dserror("SetUseTranspose not impl.");
        return -1;
      }

      /*!
      \brief not implemented
       */
      double NormInf() const
      {
        dserror("NormInf not impl.");
        return (-1.0);
      }


      /*!
      \brief not implemented
       */
      bool UseTranspose() const
      {
        // TAW: we can safely remove this. The (new) Belos now always checks the UseTransposed flag.
        // dserror("UseTranspose not impl.");
        return false;
      }

      /*!
      \brief not implemented
       */
      bool HasNormInf() const
      {
        dserror("HasNormInf not impl.");
        return false;
      }

     private:
      // don't want copy-ctor and = operator
      CheapSIMPLE_BlockPreconditioner(CheapSIMPLE_BlockPreconditioner& old);
      CheapSIMPLE_BlockPreconditioner operator=(const CheapSIMPLE_BlockPreconditioner& old);

      /*!
      \brief setup phase of preconditioner
       */
      void Setup(Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& origvlist,
          const Teuchos::ParameterList& origplist);

      /*!
      \brief do one sweep of simple or simplec preconditioning
       */
      void Simple(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px, LINALG::ANA::Vector& vb,
          LINALG::ANA::Vector& pb) const;

      /*!
      \brief do one sweep of simpler preconditioning
       */
      void Simpler(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px, LINALG::ANA::Vector& vb,
          LINALG::ANA::Vector& pb) const;


      /*!
      \brief do one sweep of simple or simplec preconditioning without subsolves
       */
      void CheapSimple(LINALG::ANA::Vector& vx, LINALG::ANA::Vector& px, LINALG::ANA::Vector& vb,
          LINALG::ANA::Vector& pb) const;

      FILE* outfile_;                     // error file to write output to or NULL
      ParameterList predictSolver_list_;  // list for primary solver
      ParameterList schurSolver_list_;    // list for secondary solver
      double alpha_;                      // pressure damping \in (0,1]

      MultiMapExtractor mmex_;  // a multimapetxractor to handle extracts (reference not an
                                // Teuchos::RCP by intention!)
      Teuchos::RCP<BlockSparseMatrixBase> A_;  // 2x2 block matrix
      Teuchos::RCP<SparseMatrix> diagAinv_;    // inverse of main diagonal of A(0,0)
      Teuchos::RCP<SparseMatrix> S_;           // Approximate Schur complement on the pressure space

      Teuchos::RCP<Epetra_Operator> Ppredict_;  // preconditioner for primary prediction subproblem
      Teuchos::RCP<Epetra_Operator> Pschur_;    // preconditioner for Schur-complement subproblem

      Teuchos::RCP<LINALG::ANA::Vector> vx_;      // velocity solution
      Teuchos::RCP<LINALG::ANA::Vector> px_;      // pressure solution
      Teuchos::RCP<LINALG::ANA::Vector> vb_;      // velocity rhs
      Teuchos::RCP<LINALG::ANA::Vector> pb_;      // pressure rhs
      Teuchos::RCP<LINALG::ANA::Vector> vwork1_;  // working vector velocity dimension
      Teuchos::RCP<LINALG::ANA::Vector> vwork2_;  // working vector velocity dimension
      Teuchos::RCP<LINALG::ANA::Vector> pwork1_;  // working vector pressure dimension
      Teuchos::RCP<LINALG::ANA::Vector> pwork2_;  // working vector pressure dimension

      Teuchos::RCP<LINALG::Solver> vsolver_;  // velocity solver // not used for CheapSIMPLE
      Teuchos::RCP<LINALG::Solver> psolver_;  // pressure solver // not used for CheapSIMPLE

      bool vdw_;  // indicate downwinding for velocity subproblem
      Teuchos::RCP<LINALG::DownwindMatrix>
          vdwind_;  // downwinding reindexer for velocity subproblem
      Teuchos::RCP<Epetra_CrsMatrix> dwA00_;
      Teuchos::RCP<LINALG::ANA::Vector> vdwin_;   // working vector
      Teuchos::RCP<LINALG::ANA::Vector> vdwout_;  // working vector

      bool pdw_;  // indicate downwinding for velocity subproblem
      Teuchos::RCP<LINALG::DownwindMatrix>
          pdwind_;  // downwinding reindexer for velocity subproblem
      Teuchos::RCP<Epetra_CrsMatrix> dwS_;
      Teuchos::RCP<LINALG::ANA::Vector> pdwin_;   // working vector
      Teuchos::RCP<LINALG::ANA::Vector> pdwout_;  // working vector

      std::string label_;  // label
    };                     // class  CheapSIMPLE_Operator : public virtual Epetra_Operator
  }                        // namespace SOLVER
}  // namespace LINALG


#endif  // CHEAPSIMPEL_BLOCKPRECONDITIONER_H
