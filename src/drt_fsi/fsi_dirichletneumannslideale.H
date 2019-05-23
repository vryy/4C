/*----------------------------------------------------------------------*/
/*!

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
       with sliding ALE-structure interfaces


\maintainer Matthias Mayr

\level 1

*/
/*----------------------------------------------------------------------*/



#ifndef FSI_DIRICHLETNEUMANNSLIPALE_H
#define FSI_DIRICHLETNEUMANNSLIPALE_H

#include "fsi_dirichletneumann.H"

namespace FSI
{
  namespace UTILS
  {
    class SlideAleUtils;
  }
  /**
   *  \brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
   *   with sliding ALE-structure interfaces
   *
   * This class implements the abstract interface FSI::DirichletNeumann
   * for the algorithm class of Dirichlet-Neumann partitioned FSI problems.
   * Specifically, this class enables tangential sliding between the solid and the ALE mesh, while
   * the tangential relative motion between solid and fluid particles adheres to a perfect stick
   * condition.
   *
   * FluidOp() takes a interface displacement, applies it to the ale
   * field, solves the ale field, calculates the interface velocities,
   * applies them to the fluid field, solves the fluid field on the
   * newly deformed fluid mesh and returns the interface forces.
   *
   * StructOp() takes interface forces, applies them to the structural
   * field, solves the field and returns the interface displacements.
   *
   * Furthermore this class contains a Remeshing() method containing the
   * the computation of rotation free ALE displacement values, remeshing
   * of the fluid field and reevaluation of the Mortar interface.
   */
  class DirichletNeumannSlideale : public DirichletNeumann
  {
    friend class DirichletNeumannFactory;

   protected:
    /**
     *  \brief almost empty ctor within saving of important things
     *
     * You will have to use the FSI::DirichletNeumannFactory to create an instance of this class
     *
     * @param[in] comm Communicator
     */
    explicit DirichletNeumannSlideale(const Epetra_Comm& comm);

   public:
    /// setup this object
    virtual void Setup();

    /**
     *  \brief Perform remeshing to account for mesh sliding at the interface
     *
     * This encompasses computation of rotation free ALE displacement values, remeshing
     * of the fluid field and reevaluation of the Mortar interface terms.
     */
    void Remeshing();

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
    virtual Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag) final;

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
    virtual Teuchos::RCP<Epetra_Vector> StructOp(
        Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag) final;

    /// predictor
    virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

   private:
    //! Sliding Ale helper class
    Teuchos::RCP<FSI::UTILS::SlideAleUtils> slideale_;

    //! Displacement of slave side of the interface
    Teuchos::RCP<Epetra_Vector> islave_;

    //! Slave displacement on master side at every time step begin
    Teuchos::RCP<Epetra_Vector> FTStemp_;
  };

}  // namespace FSI

#endif
