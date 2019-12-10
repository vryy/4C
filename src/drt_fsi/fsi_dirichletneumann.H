/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_DIRICHLETNEUMANN_H
#define FSI_DIRICHLETNEUMANN_H

#include "fsi_partitioned.H"

namespace FSI
{
  /**
   *  \brief Dirichlet-Neumann coupled interface system
   *
   * This is an abstract interface for the algorithm class of Dirichlet-Neumann
   * partitioned FSI problems. The fields themselves are solved using appropriate
   * field algorithms (that are used for standalone solvers as well.)
   * The FSI interface problem is solved using NOX.
   *
   * The FSIOp() method does one FSI cycle, that is one
   * solve of all participating fields. But for sake of clarity
   * this cycle is expressed via the two operator methods FluidOp() and
   * StructOp().
   *
   * FluidOp() takes a kinematic interface variable, applies it in some
   * shape or form to the fluid field, solves the fluid field on the
   * newly deformed fluid mesh and returns the interface forces. This needs
   * to be implemented in a derived class.
   *
   * StructOp() takes interface forces, applies them to the structural
   * field, solves the field and returns the kinematic interface variable.
   * This needs to be implemented in a derived class
   *
   * This coupling process builds on the available field solvers. However,
   * the independent parallel distribution of the fields complicates the
   * exchange of coupling information. Therefore three instances of the
   * Coupling class are used that couple those fields. On top of these
   * there are helper methods StructToAle(), StructToFluid(),
   * FluidToStruct() and AleToFluid() to easily exchange distributed
   * interface vectors between fields.
   */
  class DirichletNeumann : public Partitioned
  {
   public:
    /** \brief Constructor
     *
     * @param[in] comm Communicator
     */
    explicit DirichletNeumann(const Epetra_Comm &comm);

    /// Creates the appropriate DirichletNeumann algortihm
    //    Teuchos::RCP<DirichletNeumann> FSI::DirichletNeumann::Factory(_PROBLEM_TYP type); \FIXME

    /// setup this object
    virtual void Setup();

    /// Access function for the #kinematiccoupling_ flag
    bool GetKinematicCoupling() const { return kinematiccoupling_; }

    /** \brief Set function for the kinematiccoupling variable
     *  \param[in] variable value to write into #kinematiccoupling_
     */
    void SetKinematicCoupling(const bool variable) { kinematiccoupling_ = variable; }

   protected:
    /** \brief composed FSI operator
     *
     * The FSI operator performs the coupling iteration between the solid and the fluid field.
     *
     * \param[in] x Interface coupling variable (kinematic quantity or force), depending on the
     * value of #kinematiccoupling_.
     *
     * \param[in, out] F residual vector
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     */
    virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

    /** \brief interface fluid operator
     * \param[in] icoup kinematic interface variable
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns interface force
     */
    virtual Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> icoup, const FillType fillFlag) = 0;

    /** \brief interface structural operator
     * \param[in] iforce interface force
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns kinematic interface variable
     */
    virtual Teuchos::RCP<Epetra_Vector> StructOp(
        Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag) = 0;

    virtual Teuchos::RCP<Epetra_Vector> InitialGuess() = 0;

   private:
    /**
     * \brief  Flag to switch between kinematic and force coupling.
     *
     * - When performing kinematic coupling, we perform a fix-point scheme for the kinematic
     * variable obtained from the structure. Hence, we solve the fluid problem first.
     * - When performing force coupling, we perform a fix-point scheme for the interface force
     * obtained from the fluid field. Hence, we solve the structure problem first.
     */
    bool kinematiccoupling_;
  };

}  // namespace FSI

#endif
