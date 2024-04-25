/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all multipoint constraint equation terms

\level 3
    */
/*-----------------------------------------------------------*/
#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_EQUATION_MPC_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_EQUATION_MPC_HPP

#include "4C_config.hpp"

#include "4C_inpar_mpc_rve.hpp"
#include "4C_linalg_sparsematrix.hpp"

#include <Epetra_CrsMatrix.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN
namespace CONSTRAINTS::SUBMODELEVALUATOR
{
  /*! \brief The MultiPointConstraintEquationBase class serves as a base interface
   * for managing multi-point constraint equations within the constraint framework
   */
  class MultiPointConstraintEquationBase
  {
   public:
    virtual ~MultiPointConstraintEquationBase() = default;
    //! Constructor
    MultiPointConstraintEquationBase() = default;

    /*! \brief Add the penalty stiffness contribution to the constraint_vector and the
     * coupling-stiffness
     *
     * @param [in]  \f$Q_{dd}\f$ coupling-stiffnes matrix
     * @param [in]  \f$Q_{dL}\f$ coupling-stiffnes matrix
     * @param [in]  \f$Q_{Ld}\f$ coupling-stiffnes matrix
     * @param [in] constraint_vector constraint vector
     * @param [in] displacements \f$D_{n+1}\f$
     */
    virtual void EvaluateEquation(CORE::LINALG::SparseMatrix& Q_dd,
        CORE::LINALG::SparseMatrix& Q_dL, CORE::LINALG::SparseMatrix& Q_Ld,
        Epetra_Vector& constraint_vector, const Epetra_Vector& D_np1) = 0;

    /*! \brief Return the number of multi point constraints the object contains
     *
     * @return [out] number of multi point constraint equation the object contains
     */
    int GetNumberOfMPCs() const;

    /*! \brief Return the global id of the affected row of this equation
     *
     * @return [out] global id of the affected row of this equation
     */
    int GetFirstRowId() const;

    /*! \brief Return the global id of the affected row of this equation
     *
     * @param [in] global_row_id global id of the affected row of this equation
     */
    void SetFirstRowId(int global_row_id);

   private:
    //! Number of dof coupled per Object (= Number of MPCs per Obj.)
    int n_dof_coupled_ = 1;

    //! ID of the first constraint in the set
    int first_row_id_;
  };
  /*! \brief The class provides the method for evaluating linear coupled
   *  equations and manages associated coefficients and data.
   */
  class LinearCoupledEquation : public MultiPointConstraintEquationBase
  {
   public:
    ~LinearCoupledEquation() override = default;

    //! Default Constructor
    LinearCoupledEquation() = default;

    /*!
        \brief Standard Constructor
    */
    LinearCoupledEquation(int id, const std::vector<int>& dofs, std::vector<double> coefficients);

    //! derived
    void EvaluateEquation(CORE::LINALG::SparseMatrix& Q_dd, CORE::LINALG::SparseMatrix& Q_dL,
        CORE::LINALG::SparseMatrix& Q_Ld, Epetra_Vector& constraint_vector,
        const Epetra_Vector& D_np1) override;

   private:
    //! Struct with Term data: Coef, RowID, DofID
    struct TermData
    {
      double Coef;
      int RowId;  // is equal to the id of the MPC
      int DofId;  // dof the coefficient is multiplied with.
    };

    //! Vector with the data of the terms of a single equation
    std::vector<TermData> equation_data_;
  };
}  // namespace CONSTRAINTS::SUBMODELEVALUATOR

FOUR_C_NAMESPACE_CLOSE
#endif
