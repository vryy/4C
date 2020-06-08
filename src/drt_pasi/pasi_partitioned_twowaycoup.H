/*---------------------------------------------------------------------------*/
/*! \file
\brief two way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_TWOWAYCOUP_H
#define PASI_PARTITIONED_TWOWAYCOUP_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned.H"

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PASI
{
  /*!
   * \brief two way coupled partitioned algorithm
   *
   * Two way coupled partitioned particle structure interaction algorithm following a
   * Dirichlet-Neumann coupling scheme with particle field as Dirichlet partition and Structure
   * field as Neumann partition.
   *
   * \author Sebastian Fuchs \date 02/2017
   */
  class PASI_PartTwoWayCoup : public PartitionedAlgo
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 02/2017
     *
     * \param[in] comm   communicator
     * \param[in] params particle structure interaction parameter list
     */
    explicit PASI_PartTwoWayCoup(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PASI_PartTwoWayCoup() override{};

    /*!
     * \brief init pasi algorithm
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void Init() override;

    /*!
     * \brief setup pasi algorithm
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void Setup() override;

    /*!
     * \brief read restart information for given time step
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] restartstep restart step
     */
    void ReadRestart(int restartstep) override;

    /*!
     * \brief partitioned two way coupled timeloop
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void Timeloop() final;

   protected:
    /*!
     * \brief iteration loop between coupled fields
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    virtual void Outerloop();

    /*!
     * \brief output of fields
     *
     * \author Sebastian Fuchs \date 03/2017
     */
    void Output() override;

    /*!
     * \brief reset increment states
     *
     * Reset the interface displacement increment and the interface force increment states to the
     * interface displacement and the interface force. The increments are build after the structure
     * and particle field are solved.
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] intfdispnp  interface displacement
     * \param[in] intfforcenp interface force
     */
    void ResetIncrementStates(Teuchos::RCP<const Epetra_Vector> intfdispnp,
        Teuchos::RCP<const Epetra_Vector> intfforcenp);

    /*!
     * \brief build increment states
     *
     * Finalize the interface displacement increment and the interface force increment states.
     *
     * \author Sebastian Fuchs \date 03/2017
     */
    void BuildIncrementStates();

    /*!
     * \brief set interface forces
     *
     * Apply the interface forces as handed in to the structural field.
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] intfforcenp interface force
     */
    void SetInterfaceForces(Teuchos::RCP<const Epetra_Vector> intfforcenp);

    /*!
     * \brief reset particle states
     *
     * Reset the particle states to the converged states of the last time step.
     *
     * \author Sebastian Fuchs \date 03/2017
     */
    void ResetParticleStates();

    /*!
     * \brief clear interface forces
     *
     * Clear the interface forces in the particle wall handler.
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void ClearInterfaceForces();

    /*!
     * \brief get interface forces
     *
     * Get the interface forces via assemblation of the forces from the particle wall handler. This
     * includes communication, since the structural discretization and the particle wall
     * discretization are in general distributed independently of each other to all processors.
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void GetInterfaceForces();

    /*!
     * \brief convergence check of the outer loop
     *
     * Convergence check of the partitioned coupling outer loop based on relative and scaled
     * interface displacement and force increment norms.
     *
     * \author Sebastian Fuchs \date 02/2017
     *
     * \param[in] itnum iteration counter
     *
     * \return flag indicating status of convergence check
     */
    bool ConvergenceCheck(int itnum);

    /*!
     * \brief save particle states
     *
     * Save the converged particle states of the last time step.
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    void SaveParticleStates();

    //! interface force acting
    Teuchos::RCP<Epetra_Vector> intfforcenp_;

    //! interface displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> intfdispincnp_;

    //! interface force increment of the outer loop
    Teuchos::RCP<Epetra_Vector> intfforceincnp_;

    //! maximum iteration steps
    const int itmax_;

    //! tolerance of relative interface displacement increments in partitioned iterations
    const double convtolrelativedisp_;

    //! tolerance of dof and dt scaled interface displacement increments in partitioned iterations
    const double convtolscaleddisp_;

    //! tolerance of relative interface force increments in partitioned iterations
    const double convtolrelativeforce_;

    //! tolerance of dof and dt scaled interface force increments in partitioned iterations
    const double convtolscaledforce_;

    //! ignore convergence check and proceed simulation
    const bool ignoreconvcheck_;

    //! write restart every n steps
    const int writerestartevery_;
  };

  /*!
   * \brief two way coupled partitioned algorithm with constant interface displacement relaxation
   *
   * Two way coupled partitioned particle structure interaction algorithm following a
   * Dirichlet-Neumann coupling scheme with particle field as Dirichlet partition and Structure
   * field as Neumann partition and constant interface displacement relaxation.
   *
   * \author Sebastian Fuchs \date 03/2017
   */
  class PASI_PartTwoWayCoup_DispRelax : public PASI_PartTwoWayCoup
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] comm   communicator
     * \param[in] params particle structure interaction parameter list
     */
    explicit PASI_PartTwoWayCoup_DispRelax(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PASI_PartTwoWayCoup_DispRelax() override{};

   protected:
    /*!
     * \brief iteration loop between coupled fields with relaxed displacements
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void Outerloop() override;

    /*!
     * \brief calculate relaxation parameter
     *
     * No computation of the relaxation parameter necessary in the constant case.
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] omega relaxation parameter
     * \param[in] itnum iteration counter
     */
    virtual void CalcOmega(double& omega, const int itnum);

    //! relaxation parameter
    double omega_;

   private:
    /*!
     * \brief perform relaxation of interface states
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in,out] intfdispnp relaxed interface displacement
     * \param[out]    intfvelnp  relaxed interface velocity
     * \param[out]    intfaccnp  relaxed interface acceleration
     */
    void PerformRelaxationInterfaceStates(Teuchos::RCP<Epetra_Vector> intfdispnp,
        Teuchos::RCP<Epetra_Vector> intfvelnp, Teuchos::RCP<Epetra_Vector> intfaccnp);
  };

  /*!
   * \brief two way coupled partitioned algorithm with dynamic interface displacement relaxation
   *
   * Two way coupled partitioned particle structure interaction algorithm following a
   * Dirichlet-Neumann coupling scheme with particle field as Dirichlet partition and Structure
   * field as Neumann partition and dynamic interface displacement relaxation following Aitken's
   * delta^2 method.
   *
   * \author Sebastian Fuchs \date 03/2017
   */
  class PASI_PartTwoWayCoup_DispRelaxAitken : public PASI_PartTwoWayCoup_DispRelax
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] comm   communicator
     * \param[in] params particle structure interaction parameter list
     */
    PASI_PartTwoWayCoup_DispRelaxAitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PASI_PartTwoWayCoup_DispRelaxAitken() override{};

    /*!
     * \brief init pasi algorithm
     *
     * \author Sebastian Fuchs \date 03/2017
     */
    void Init() override;

    /*!
     * \brief read restart information for given time step
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] restartstep restart step
     */
    void ReadRestart(int restartstep) override;

   protected:
    /*!
     * \brief output of fields
     *
     * \author Sebastian Fuchs \date 03/2017
     */
    void Output() override;

    /*!
     * \brief calculate relaxation parameter
     *
     * Computation of the relaxation parameter following Aitken's delta^2 method.
     *
     * Refer to PhD thesis U. Kuettler, equation (3.5.29).
     *
     * \author Sebastian Fuchs \date 03/2017
     *
     * \param[in] omega relaxation parameter
     * \param[in] itnum iteration counter
     */
    void CalcOmega(double& omega, const int itnum) override;

    //! old interface displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> intfdispincnpold_;

    //! maximal relaxation parameter
    double maxomega_;

    //! minimal relaxation parameter
    double minomega_;
  };

}  // namespace PASI

/*---------------------------------------------------------------------------*/
#endif
