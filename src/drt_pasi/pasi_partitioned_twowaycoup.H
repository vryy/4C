/*---------------------------------------------------------------------------*/
/*! \file
\brief two way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_TWOWAYCOUP_H
#define PASI_PARTITIONED_TWOWAYCOUP_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
namespace PASI
{
  //! two way coupled partitioned algorithm
  class PASI_PartTwoWayCoup : public PartitionedAlgo
  {
   public:
    //! constructor
    explicit PASI_PartTwoWayCoup(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PASI_PartTwoWayCoup() override{};

    //! init pasi algorithm
    void Init() override;

    //! setup pasi algorithm
    void Setup() override;

    //! read restart information for given time step
    void ReadRestart(int restartstep) override;

    //! partitioned two way coupled timeloop
    void Timeloop() final;

   protected:
    //! iteration loop between coupled fields
    virtual void Outerloop();

    //! output of fields
    void Output() override;

    //! reset increment states
    void ResetIncrementStates(Teuchos::RCP<const Epetra_Vector> intfdispnp,
        Teuchos::RCP<const Epetra_Vector> intfforceincnp);

    //! build increment states
    void BuildIncrementStates();

    //! set interface forces
    void SetInterfaceForces(Teuchos::RCP<const Epetra_Vector> intfforcenp);

    //! reset particle states
    void ResetParticleStates();

    //! clear interface forces
    void ClearInterfaceForces();

    //! get interface forces
    void GetInterfaceForces();

    //! convergence check for structure and particles fields
    bool ConvergenceCheck(int itnum);

    //! save particle states
    void SaveParticleStates();

    //! interface force acting on structural boundary
    Teuchos::RCP<Epetra_Vector> intfforcenp_;

    //! interface displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> intfdispincnp_;

    //! interface force increment of the outer loop
    Teuchos::RCP<Epetra_Vector> intfforceincnp_;

    //! maximum iteration steps
    int itmax_;

    //! convergence tolerance
    double ittol_;

    //! ignore convergence check and proceed simulation
    bool ignoreconvcheck_;

    //! write restart every n steps
    int writerestartevery_;
  };

  //! two way coupled partitioned algorithm with displacement relaxation
  class PASI_PartTwoWayCoup_DispRelax : public PASI_PartTwoWayCoup
  {
   public:
    //! constructor
    explicit PASI_PartTwoWayCoup_DispRelax(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PASI_PartTwoWayCoup_DispRelax() override{};

   protected:
    //! iteration loop between coupled fields with relaxed forces
    void Outerloop() override;

    //! calculate relaxation parameter
    virtual void CalcOmega(double& omega, const int itnum);

    //! relaxation parameter
    double omega_;

   private:
    //! perform relaxation of interface states
    void PerformRelaxationInterfaceStates(Teuchos::RCP<Epetra_Vector> intfdispnp,
        Teuchos::RCP<Epetra_Vector> intfvelnp, Teuchos::RCP<Epetra_Vector> intfaccnp);
  };

  //! two way coupled partitioned algorithm with displacement relaxation via aitken
  class PASI_PartTwoWayCoup_DispRelaxAitken : public PASI_PartTwoWayCoup_DispRelax
  {
   public:
    //! constructor
    PASI_PartTwoWayCoup_DispRelaxAitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~PASI_PartTwoWayCoup_DispRelaxAitken() override{};

    //! init pasi algorithm
    void Init() override;

    //! read restart information for given time step
    void ReadRestart(int restartstep) override;

   protected:
    //! output of fields
    void Output() override;

    //! calculate relaxation parameter
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
