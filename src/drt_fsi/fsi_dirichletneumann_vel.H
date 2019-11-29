/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann_vel.H

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

\maintainer Nora Hagmeyer

\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_DIRICHLETNEUMANN_DIRICH_H
#define FSI_DIRICHLETNEUMANN_DIRICH_H

#include "fsi_dirichletneumann.H"

namespace FSI
{
  /**
   * \brief Dirichlet-Neumann interface velocity based algorithm
   *
   * This is the algorithm class of partitioned FSI problems. Here we
   * do the time loop and the coupling between fields. The fields
   * themselves are solved using appropriate field algorithms (that are
   * used for stand-alone solvers as well.) The FSI interface problem is
   * solved using NOX.
   *
   * Many different things come together in this class. First and
   * foremost there is the Timeloop() method, that contains the overall
   * FSI time stepping scheme. The time loop is build around the FSI
   * interface problem, that is the problem of finding the correct interface
   * coupling quantities that satisfy the coupled problem. The time loop knows
   * nothing about the fields themselves.
   *
   * Inside the time loop the interface problem is solved using NOX. To
   * do so Timeloop() needs to know an object which defines the nonlinear
   * residual of the FSI problem. This is always a Teuchos::RCP to the
   * DirichletNeumannCoupling object itself!
   *
   * So the second part of this class consists of the interface residual
   * evaluation computeF(). This method does one FSI cycle, that is one
   * solve of all participating fields. But for sake of clarity
   * this cycle is expressed via the two operator methods FluidOp() and
   * StructOp().
   *
   * FluidOp() takes an kinematic interface variable, applies it in some
   * shape or form to the fluid field, solves the fluid field on the
   * newly deformed fluid mesh and returns the interface forces.
   *
   * StructOp() takes interface forces, applies them to the structural
   * field, solves the field and returns the interface displacements.
   *
   * This coupling process builds on the available field solvers. However,
   * the independent parallel distribution of the fields complicates the
   * exchange of coupling information. Therefore three instances of the
   * Coupling class are used that couple those fields. On top of these
   * there are helper methods StructToAle(), StructToFluid(),
   * FluidToStruct() and AleToFluid() to easily exchange distributed
   * interface vectors between fields.
   *
   * The FSI algorithm requires repeated evaluations of the interface
   * residual via computeF(). So the field solvers themselves must be
   * clean, subsequent calls with the same interface input must yield the
   * same results. The time stepping therefore needs the further methods
   * PrepareTimeStep() to start a new time step as well as Update() and
   * Output() to finish the current step, save the result and write the
   * files.
   */
  class DirichletNeumannVel : public DirichletNeumann
  {
    friend class DirichletNeumannFactory;

   protected:
    /**
     *  \brief constructor
     *
     * You will have to use the FSI::DirichletNeumannFactory to create an instance of this class
     */
    explicit DirichletNeumannVel(const Epetra_Comm &comm);

   public:
    /// setup this object
    virtual void Setup() override;

   protected:
    /** \brief interface fluid operator
     *
     * In here, the nonlinear solve for the fluid field is prepared and called and the resulting
     * interface force is returned
     *
     * \param[in] ivel The interface velocity
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     */
    virtual Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> ivel, const FillType fillFlag) override;

    /** \brief interface structural operator
     *
     * In here, the nonlinear solve for the structure field is prepared and called and the resulting
     * interface velocity is returned
     *
     * \param[in] iforce The interface force
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     */
    virtual Teuchos::RCP<Epetra_Vector> StructOp(
        Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag) override;

    /// Computes initial guess for the next iteration
    virtual Teuchos::RCP<Epetra_Vector> InitialGuess() override;

    // todo use or throw error
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) override
    {
      return Teuchos::null;
    };

    // todo use or throw error
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) override
    {
      return Teuchos::null;
    };

   private:
    /**
     * \brief Object that allows to capsule the different constraint enforcement strategies and
     * effectively separating it from the actual algorithm
     */
    Teuchos::RCP<ADAPTER::FBIConstraintenforcer> constraint_manager_;
  };
}  // namespace FSI

#endif
