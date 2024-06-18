/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
coupling the displacements


\level 1
*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_DIRICHLETNEUMANN_DISP_HPP
#define FOUR_C_FSI_DIRICHLETNEUMANN_DISP_HPP

#include "4C_config.hpp"

#include "4C_fsi_dirichletneumann.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  /**
   *  \brief Dirichlet-Neumann displacement-based algorithm
   *
   *  This class implements the abstract interface FSI::DirichletNeumann
   *  for the algorithm class of Dirichlet-Neumann partitioned FSI problems.
   *  Specifically, this algorithm performs coupling of the solid and ALE displacements.
   *
   *  FluidOp() takes an interface displacement, applies it to the ale
   *  field, solves the ale field, calculates the interface velocities,
   *  applies them to the fluid field, solves the fluid field on the
   *  newly deformed fluid mesh and returns the interface forces.
   *
   *  StructOp() takes interface forces, applies them to the structural
   *  field, solves the field and returns the interface displacements.
   */
  class DirichletNeumannDisp : public DirichletNeumann
  {
    friend class DirichletNeumannFactory;

   protected:
    /**
     *  \brief constructor
     *
     * You will have to use the FSI::DirichletNeumannFactory to create an instance of this class
     *
     * \param[in] comm Communicator
     */
    explicit DirichletNeumannDisp(const Epetra_Comm &comm);

   public:
    /// setup this object
    void setup() override;

   protected:
    /** \brief interface fluid operator
     *
     * Solve the fluid field problem.  Since the fluid field is the Dirichlet partition, the
     * interface displacement is prescribed as a Dirichlet boundary condition.
     *
     * \param[in] idisp interface displacement
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns interface force
     */
    Teuchos::RCP<Epetra_Vector> fluid_op(
        Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag) override;

    /** \brief interface structural operator
     *
     * Solve the structure field problem.  Since the structure field is the Neumann partition, the
     * interface forces are prescribed as a Neumann boundary condition.
     *
     * \param[in] iforce interface force
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns interface displacement
     */
    Teuchos::RCP<Epetra_Vector> struct_op(
        Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag) final;

    //! Predictor
    Teuchos::RCP<Epetra_Vector> initial_guess() override;
  };

}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
