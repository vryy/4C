/*!----------------------------------------------------------------------
\file pasi_partitioned_twowaycoup.H

\brief two way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                           sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
#ifndef PASI_PARTITIONED_TWOWAYCOUP_H
#define PASI_PARTITIONED_TWOWAYCOUP_H

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
#include "pasi_partitioned.H"

/*----------------------------------------------------------------------*
 | two way coupled pasi algorithm                        sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
namespace PASI
{
  //! two way coupled partitioned algorithm for particle structure interaction
  class PASI_PartTwoWayCoup : public PartitionedAlgo
  {
   public:
    //! constructor
    explicit PASI_PartTwoWayCoup(const Epetra_Comm& comm,  //! communicator
        const Teuchos::ParameterList&
            pasi_params  //! input parameters for particle structure interaction
    );

    //! virtual destructor
    virtual ~PASI_PartTwoWayCoup(){};

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return int
    \date 01/17
    \author rauch  */
    virtual void Init(const Epetra_Comm& comm  //! communicator
    );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 01/17
    \author rauch  */
    virtual void Setup();

    //! read restart data
    virtual void ReadRestart(int step  //! time step for restart
    );

    //! timeloop of coupled problem
    virtual void Timeloop();

   protected:
    //! iteration loop between coupled fields
    virtual void Outerloop();

    //! update and output
    virtual void UpdateOutput();

    //! reset particle states in outerloop iteration
    virtual void ResetParticleStates();

    //! get particle forces
    virtual void GetParticleForces();

    //! set particle forces
    virtual void SetParticleForces(
        Teuchos::RCP<const Epetra_Vector> forcenp  //! particle wall forces at \f$t_{n+1}\f$
    );

    //! update the current states in every iteration
    virtual void IterUpdateStates(
        Teuchos::RCP<const Epetra_Vector> dispnp,  //! structural displacements at \f$t_{n+1}\f$
        Teuchos::RCP<const Epetra_Vector> forcenp  //! particle wall forces at \f$t_{n+1}\f$
    );

    //! convergence check for structure and particles fields
    virtual bool ConvergenceCheck(int itnum  //! number of iterations
    );

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

  };  // class PASI_PartTwoWayCoup: public PartitionedAlgo

  //! two way coupled partitioned algorithm for particle structure interaction with force relaxation
  class PASI_PartTwoWayCoup_ForceRelax : public PASI_PartTwoWayCoup
  {
   public:
    //! constructor
    explicit PASI_PartTwoWayCoup_ForceRelax(const Epetra_Comm& comm,  //! communicator
        const Teuchos::ParameterList&
            pasi_params  //! input parameters for particle structure interaction
    );

    //! virtual destructor
    virtual ~PASI_PartTwoWayCoup_ForceRelax(){};

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return int
    \date 01/17
    \author rauch  */
    virtual void Init(const Epetra_Comm& comm  //! communicator
    );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 01/17
    \author rauch  */
    virtual void Setup();

    //! read restart data
    virtual void ReadRestart(int step  //! time step for restart
    );

   protected:
    //! iteration loop between coupled fields with relaxed forces
    virtual void Outerloop();

    //! update and output
    virtual void UpdateOutput();

    //! calculate relaxation parameter
    virtual void CalcOmega(double& omega,  //! relaxation parameter
        const int itnum                    //! number of partitioned coupling iterations
    );

    //! relaxation parameter
    double omega_;

  };  // class PASI_PartTwoWayCoup_ForceRelax: public PASI_PartTwoWayCoup

  //! two way coupled partitioned algorithm for particle structure interaction with force relaxation
  //! via aitken
  class PASI_PartTwoWayCoup_ForceRelaxAitken : public PASI_PartTwoWayCoup_ForceRelax
  {
   public:
    //! constructor
    PASI_PartTwoWayCoup_ForceRelaxAitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    virtual ~PASI_PartTwoWayCoup_ForceRelaxAitken(){};

    /*!
    \brief Setup this object

     Initializes members and performs problem specific setup.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch
    */
    virtual void Init(const Epetra_Comm& comm  //! communicator
    );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 01/17
    \author rauch  */
    virtual void Setup();

   protected:
    //! calculate relaxation parameter
    virtual void CalcOmega(double& omega,  //! relaxation parameter
        const int itnum                    //! number of partitioned coupling iterations
    );

    //! old particle force increment of the outer loop
    Teuchos::RCP<Epetra_Vector> forceincnpold_;

    //! maximal relaxation parameter
    double maxomega_;

    //! minimal relaxation parameter
    double minomega_;

  };  // class PASI_PartTwoWayCoup_ForceRelaxAitken : public PASI_PartTwoWayCoup_ForceRelax

}  // namespace PASI

/*----------------------------------------------------------------------*/
#endif  // PASI_PARTITIONED_TWOWAYCOUP_H
