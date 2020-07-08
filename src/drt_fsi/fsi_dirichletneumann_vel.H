/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach


\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_DIRICHLETNEUMANN_DIRICH_H
#define FSI_DIRICHLETNEUMANN_DIRICH_H

#include "fsi_dirichletneumann.H"

// Forward declarations

namespace ADAPTER
{
  class FBIConstraintenforcer;
}

namespace BEAMINTERACTION
{
  class BeamToFluidMeshtyingVtkOutputWriter;
}

namespace FSI
{
  /**
   * \brief Dirichlet-Neumann interface velocity based algorithm
   *
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
    explicit DirichletNeumannVel(const Epetra_Comm& comm);

   public:
    /*! \brief Outer level FSI time loop
     *
     * We overload this interface here in order to carry out operations that only have to be done
     * once at the start of the simulation, but need information which are not available during
     * setup in the case of a restart.
     *
     *  \param[in] interface Our interface to NOX
     */
    void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface) override;

    /** \brief Here we decide which type of coupling we are going to use
     *
     * Here we check the input for the coupling variable
     *
     */
    void Setup() override;

    /** \brief Here the base class writes output for each field and in addition we write coupling
     * related output
     *
     */
    void Output() override;

   protected:
    /** \brief interface fluid operator
     *
     * In here, the nonlinear solve for the fluid field is prepared and called and the resulting
     * interface force is returned
     *
     * \param[in] ivel The interface velocity
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     */
    Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> ivel, const FillType fillFlag) override;

    /** \brief interface structural operator
     *
     * In here, the nonlinear solve for the structure field is prepared and called and the resulting
     * interface velocity is returned
     *
     * \param[in] iforce The interface force
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     */
    Teuchos::RCP<Epetra_Vector> StructOp(
        Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag) override;

    /// Computes initial guess for the next iteration
    Teuchos::RCP<Epetra_Vector> InitialGuess() override;

    /**
     * \brief In here all coupling related quantities are given to the fluid solver
     *
     * \param[in] iv In our case, the input parameter is not used!
     * \returns ivel The fluid velocity on the (whole) fluid domain
     */
    Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) override;

    /**
     * \brief In here all coupling related quantities are assembled for the structure solver
     *
     * \param[in] iv In our case, the input parameter is not used!
     * \returns iforce The fsi force acting on the structure
     */

    Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) override;


   private:
    /**
     * \brief Object that allows to capsule the different constraint enforcement strategies and
     * effectively separating it from the actual algorithm
     */
    Teuchos::RCP<ADAPTER::FBIConstraintenforcer> constraint_manager_;

    Teuchos::RCP<BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter> vtk_output_writer_;
  };
}  // namespace FSI

#endif
