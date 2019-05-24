/*---------------------------------------------------------------------------*/
/*!
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

    //! virtual destructor
    virtual ~PASI_PartTwoWayCoup(){};

    //! init pasi algorithm
    virtual void Init() override;

    //! setup pasi algorithm
    virtual void Setup() override;

    //! read restart information for given time step
    virtual void ReadRestart(int restartstep) override;

    //! partitioned two way coupled timeloop
    virtual void Timeloop() override;

   protected:
    //! iteration loop between coupled fields
    virtual void Outerloop();

    //! output of fields
    virtual void Output();

    //! update the current states in every iteration
    void IterUpdateStates(
        Teuchos::RCP<const Epetra_Vector> dispnp, Teuchos::RCP<const Epetra_Vector> forcenp);

    //! set wall forces
    void SetWallForces(Teuchos::RCP<const Epetra_Vector> forcenp);

    //! reset particle states
    void ResetParticleStates();

    //! clear wall forces
    void ClearWallForces();

    //! get wall forces
    void GetWallForces();

    //! convergence check for structure and particles fields
    bool ConvergenceCheck(int itnum);

    //! save particle states
    void SaveParticleStates();

    //! particle force acting on structural boundary
    Teuchos::RCP<Epetra_Vector> forcenp_;

    //! structural displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

    //! particle force increment of the outer loop
    Teuchos::RCP<Epetra_Vector> forceincnp_;

    //! maximum iteration steps
    int itmax_;

    //! convergence tolerance
    double ittol_;

    //! ignore convergence check and proceed simulation
    bool ignoreconvcheck_;

    //! write restart every n steps
    int writerestartevery_;
  };

  //! two way coupled partitioned algorithm with force relaxation
  class PASI_PartTwoWayCoup_ForceRelax : public PASI_PartTwoWayCoup
  {
   public:
    //! constructor
    explicit PASI_PartTwoWayCoup_ForceRelax(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~PASI_PartTwoWayCoup_ForceRelax(){};

   protected:
    //! iteration loop between coupled fields with relaxed forces
    virtual void Outerloop() override;

    //! calculate relaxation parameter
    virtual void CalcOmega(double& omega, const int itnum);

    //! relaxation parameter
    double omega_;
  };

  //! two way coupled partitioned algorithm with force relaxation via aitken
  class PASI_PartTwoWayCoup_ForceRelaxAitken : public PASI_PartTwoWayCoup_ForceRelax
  {
   public:
    //! constructor
    PASI_PartTwoWayCoup_ForceRelaxAitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~PASI_PartTwoWayCoup_ForceRelaxAitken(){};

    //! init pasi algorithm
    virtual void Init() override;

    //! read restart information for given time step
    virtual void ReadRestart(int restartstep) override;

   protected:
    //! output of fields
    virtual void Output() override;

    //! calculate relaxation parameter
    virtual void CalcOmega(double& omega, const int itnum) override;

    //! old particle force increment of the outer loop
    Teuchos::RCP<Epetra_Vector> forceincnpold_;

    //! maximal relaxation parameter
    double maxomega_;

    //! minimal relaxation parameter
    double minomega_;
  };

}  // namespace PASI

/*---------------------------------------------------------------------------*/
#endif
