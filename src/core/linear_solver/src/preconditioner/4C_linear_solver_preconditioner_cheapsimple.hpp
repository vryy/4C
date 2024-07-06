/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_CHEAPSIMPLE_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_CHEAPSIMPLE_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_linalg_ana.hpp"
#include "4C_linear_solver_preconditioner_type.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinearSolver
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

   */
  class CheapSimpleBlockPreconditioner : public virtual Epetra_Operator
  {
   public:
    /*!
        \brief Standard Constructor
    */
    explicit CheapSimpleBlockPreconditioner(Teuchos::RCP<Epetra_Operator> A,
        const Teuchos::ParameterList& predict_list, const Teuchos::ParameterList& correct_list);

    /*!
        \brief Destructor
    */
    ~CheapSimpleBlockPreconditioner() override
    {
      vsolver_ = Teuchos::null;  // solver objects only generated for SIMPLE and SIMPLER (not for
                                 // CheapSIMPLE)
      psolver_ = Teuchos::null;
      ppredict_ = Teuchos::null;
      pschur_ = Teuchos::null;
    }

    /*!
     \brief Setup the label of this class.
    */
    std::string setup_label()
    {
      std::stringstream strLabel;

      strLabel << "CheapSIMPLE ";

      const bool visml = predict_solver_list_.isSublist("ML Parameters");
      const bool vismuelu = predict_solver_list_.isSublist("MueLu Parameters");
      const bool visifpack = predict_solver_list_.isSublist("IFPACK Parameters");
      strLabel << "[Inv1=";
      if (visml)
        strLabel << "ML-type";
      else if (vismuelu)
        strLabel << "MueLu-type";
      else if (visifpack)
      {
        const Teuchos::ParameterList& ifpackParams =
            predict_solver_list_.sublist("IFPACK Parameters");
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
      const bool pisml = schur_solver_list_.isSublist("ML Parameters");
      const bool pismuelu = predict_solver_list_.isSublist("MueLu Parameters");
      const bool pisifpack = schur_solver_list_.isSublist("IFPACK Parameters");
      if (pisml)
        strLabel << "ML-type";
      else if (pismuelu)
        strLabel << "MueLu-type";
      else if (pisifpack)
      {
        const Teuchos::ParameterList& ifpackParams =
            schur_solver_list_.sublist("IFPACK Parameters");
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
    const char* Label() const override { return label_.data(); }

    /*!
     \brief get Comm of this class

     Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class

     */
    const Epetra_Comm& Comm() const override { return (a_->Comm()); }


    /*!
     \brief Get fine level OperatorDomainMap

     Derived from Epetra_Operator, get fine level OperatorDomainMap

     */
    const Epetra_Map& OperatorDomainMap() const override { return a_->full_domain_map(); }

    /*!
      \brief Get fine level OperatorRangeMap

      Derived from Epetra_Operator, get fine level OperatorRangeMap

    */
    const Epetra_Map& OperatorRangeMap() const override { return a_->full_range_map(); }

    /*!
      \brief ApplyInverse the preconditioner

       ApplyInverse the preconditioner. Method is derived from Epetra_Operator.


       \param X   (In) : Epetra_MultiVector matching the fine level map of this
                         preconditioner
       \param Y (Out)  : Epetra_MultiVector containing the result on output
    */
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /*!
      \brief not implemented
     */
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
    {
      FOUR_C_THROW("Apply does not make sense for Core::LinAlg::SIMPLER_Operator");
      return (-1);
    }

    /*!
      \brief not implemented
    */
    int SetUseTranspose(bool UseTranspose) override
    {
      FOUR_C_THROW("SetUseTranspose not impl.");
      return -1;
    }

    /*!
      \brief not implemented
    */
    double NormInf() const override
    {
      FOUR_C_THROW("NormInf not impl.");
      return (-1.0);
    }

    /*!
     \brief not implemented
    */
    bool UseTranspose() const override
    {
      // TAW: we can safely remove this. The (new) Belos now always checks the UseTransposed flag.
      // FOUR_C_THROW("UseTranspose not impl.");
      return false;
    }

    /*!
      \brief not implemented
    */
    bool HasNormInf() const override
    {
      FOUR_C_THROW("HasNormInf not impl.");
      return false;
    }

   private:
    // don't want copy-ctor and = operator
    CheapSimpleBlockPreconditioner(CheapSimpleBlockPreconditioner& old);
    CheapSimpleBlockPreconditioner operator=(const CheapSimpleBlockPreconditioner& old);

    /*!
      \brief setup phase of preconditioner
    */
    void setup(Teuchos::RCP<Epetra_Operator> A, const Teuchos::ParameterList& origvlist,
        const Teuchos::ParameterList& origplist);

    /*!
      \brief do one sweep of simple or simplec preconditioning
    */
    void simple(Core::LinAlg::Ana::Vector& vx, Core::LinAlg::Ana::Vector& px,
        Core::LinAlg::Ana::Vector& vb, Core::LinAlg::Ana::Vector& pb) const;

    /*!
      \brief do one sweep of simpler preconditioning
    */
    void simpler(Core::LinAlg::Ana::Vector& vx, Core::LinAlg::Ana::Vector& px,
        Core::LinAlg::Ana::Vector& vb, Core::LinAlg::Ana::Vector& pb) const;


    /*!
      \brief do one sweep of simple or simplec preconditioning without subsolves
    */
    void cheap_simple(Core::LinAlg::Ana::Vector& vx, Core::LinAlg::Ana::Vector& px,
        Core::LinAlg::Ana::Vector& vb, Core::LinAlg::Ana::Vector& pb) const;

    Teuchos::ParameterList predict_solver_list_;  // list for primary solver
    Teuchos::ParameterList schur_solver_list_;    // list for secondary solver
    double alpha_;                                // pressure damping \in (0,1]

    Core::LinAlg::MultiMapExtractor mmex_;  // a multimapetxractor to handle extracts (reference not
                                            // an Teuchos::RCP by intention!)
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> a_;  // 2x2 block matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> diag_ainv_;   // inverse of main diagonal of A(0,0)
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        s_;  // Approximate Schur complement on the pressure space

    Teuchos::RCP<Epetra_Operator> ppredict_;  // preconditioner for primary prediction subproblem
    Teuchos::RCP<Epetra_Operator> pschur_;    // preconditioner for Schur-complement subproblem

    Teuchos::RCP<Core::LinAlg::Ana::Vector> vx_;      // velocity solution
    Teuchos::RCP<Core::LinAlg::Ana::Vector> px_;      // pressure solution
    Teuchos::RCP<Core::LinAlg::Ana::Vector> vb_;      // velocity rhs
    Teuchos::RCP<Core::LinAlg::Ana::Vector> pb_;      // pressure rhs
    Teuchos::RCP<Core::LinAlg::Ana::Vector> vwork1_;  // working vector velocity dimension
    Teuchos::RCP<Core::LinAlg::Ana::Vector> vwork2_;  // working vector velocity dimension
    Teuchos::RCP<Core::LinAlg::Ana::Vector> pwork1_;  // working vector pressure dimension
    Teuchos::RCP<Core::LinAlg::Ana::Vector> pwork2_;  // working vector pressure dimension

    Teuchos::RCP<Core::LinAlg::Solver> vsolver_;  // velocity solver // not used for CheapSIMPLE
    Teuchos::RCP<Core::LinAlg::Solver> psolver_;  // pressure solver // not used for CheapSIMPLE

    bool vdw_;  // indicate downwinding for velocity subproblem
    Teuchos::RCP<Core::LinAlg::DownwindMatrix>
        vdwind_;  // downwinding reindexer for velocity subproblem
    Teuchos::RCP<Epetra_CrsMatrix> dw_a00_;
    Teuchos::RCP<Core::LinAlg::Ana::Vector> vdwin_;   // working vector
    Teuchos::RCP<Core::LinAlg::Ana::Vector> vdwout_;  // working vector

    bool pdw_;  // indicate downwinding for velocity subproblem
    Teuchos::RCP<Core::LinAlg::DownwindMatrix>
        pdwind_;  // downwinding reindexer for velocity subproblem
    Teuchos::RCP<Epetra_CrsMatrix> dw_s_;
    Teuchos::RCP<Core::LinAlg::Ana::Vector> pdwin_;   // working vector
    Teuchos::RCP<Core::LinAlg::Ana::Vector> pdwout_;  // working vector

    std::string label_;  // label
  };
}  // namespace Core::LinearSolver

FOUR_C_NAMESPACE_CLOSE

#endif
