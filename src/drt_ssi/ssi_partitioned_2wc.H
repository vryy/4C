/*!------------------------------------------------------------------------------------------------*
 \brief two way coupled partitioned scalar structure interaction

 \level 2

 \maintainer Christoph Schmidt

 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED2_WC_H_
#define SSI_PARTITIONED2_WC_H_

#include "ssi_partitioned.H"
#include "../drt_inpar/inpar_ssi.H"

namespace SSI
{
  //! base class to deal with partioned 2WC SSI. Mainly it is the same as for thermo-structure
  //! interaction (TSI)
  class SSI_Part2WC : public SSI_Part
  {
   public:
    //! constructor
    SSI_Part2WC(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    virtual ~SSI_Part2WC(){};

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
    \date 08/16
    \author rauch  */
    virtual int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle);

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
    \date 08/16
    \author rauch  */
    virtual void Setup();

    /*! \brief Initialize discretizations

    \date 08/16
    \author rauch  */
    virtual void InitDiscretizations(const Epetra_Comm& comm, const std::string& struct_disname,
        const std::string& scatra_disname);

    //! full time loop
    virtual void Timeloop();

    //! perform iteration loop between fields
    virtual void OuterLoop();

    //! prepare time loop
    virtual void PrepareTimeLoop();

    //! prepare time step for single fields
    virtual void PrepareTimeStep(bool printheader = true);

    //! update time step and print to screen
    virtual void UpdateAndOutput();

   protected:
    //! solve field 1
    virtual void Operator1() { DoStructStep(); };

    //! solve field 2
    virtual void Operator2() { DoScatraStep(); };

    //! pre operator called before first field operator
    virtual void PreOperator1();

    //! pre operator called before second field operator
    virtual void PreOperator2(){};

    //! post operator called after first field operator
    virtual void PostOperator1(){};

    //! post operator called after second field operator
    virtual void PostOperator2(){};

    //! perform iteration step of structure field and set the new disp and vel states in the scatra
    //! field
    virtual void DoStructStep();

    //! perform iteration step of scatra field and set the new phi state in the structure field
    virtual void DoScatraStep();

    //! update the current states in every iteration
    //! states are set to the last solutions obtained
    virtual void IterUpdateStates();

    //! convergence check of outer loop
    virtual bool ConvergenceCheck(int itnum);

    /// velocity calculation given the displacements
    Teuchos::RCP<Epetra_Vector> CalcVelocity(Teuchos::RCP<const Epetra_Vector> dispnp);

    //! scalar increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnp_;
    //! displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

    //! convergence tolerance
    double ittol_;
    //! maximum iteration steps
    int itmax_;
  };

  //! class to deal with displacement relaxated partioned 2WC SSI. Relaxation parameter is constant
  class SSI_Part2WC_SolidToScatra_Relax : public SSI_Part2WC
  {
   public:
    //! constructor
    SSI_Part2WC_SolidToScatra_Relax(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    virtual ~SSI_Part2WC_SolidToScatra_Relax(){};

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
    virtual int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle);

   protected:
    //! perform iteration loop between fields with relaxed displacements
    void OuterLoop();

    //! //! calculate relaxation parameter
    virtual void CalcOmega(double& omega, const int itnum);

    //! relaxation parameter
    double omega_;
  };

  //! class to deal with displacement relaxated 2WC partioned SSI. Relaxation parameter is
  //! calculated via Aitken
  class SSI_Part2WC_SolidToScatra_Relax_Aitken : public SSI_Part2WC_SolidToScatra_Relax
  {
   public:
    //! constructor
    SSI_Part2WC_SolidToScatra_Relax_Aitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    virtual ~SSI_Part2WC_SolidToScatra_Relax_Aitken(){};

    /*!
    \brief Setup all class internal objects and members

    Setup() is not supposed to have any input arguments !

    Must only be called after Init().

    Construct all objects depending on the parallel distribution and
    relying on valid maps like, e.g. the state vectors, system matrices, etc.

    Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 01/18
    \author fang
    */
    void Setup();

   protected:
    //! Calculate relaxation parameter via Aitken
    virtual void CalcOmega(double& omega, const int itnum);

    //! old displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnpold_;
  };

  //! class to deal with scalar relaxated 2WC partioned SSI. Relaxation parameter is constant
  class SSI_Part2WC_ScatraToSolid_Relax : public SSI_Part2WC
  {
   public:
    //! constructor
    SSI_Part2WC_ScatraToSolid_Relax(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    virtual ~SSI_Part2WC_ScatraToSolid_Relax(){};

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
    virtual int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle);

   protected:
    //! perform iteration loop between fields with relaxed scalar
    void OuterLoop();

    //! calculate relaxation parameter
    virtual void CalcOmega(double& omega, const int itnum);

    //! relaxation parameter
    double omega_;
  };

  //! class to deal with scalar relaxated 2WC partioned SSI. Relaxation parameter is calculated via
  //! Aitken
  class SSI_Part2WC_ScatraToSolid_Relax_Aitken : public SSI_Part2WC_ScatraToSolid_Relax
  {
   public:
    //! constructor
    SSI_Part2WC_ScatraToSolid_Relax_Aitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    virtual ~SSI_Part2WC_ScatraToSolid_Relax_Aitken(){};

    /*!
    \brief Setup all class internal objects and members

    Setup() is not supposed to have any input arguments !

    Must only be called after Init().

    Construct all objects depending on the parallel distribution and
    relying on valid maps like, e.g. the state vectors, system matrices, etc.

    Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 01/18
    \author fang
    */
    void Setup();

   protected:
    //! Calculate relaxation parameter via Aitken
    virtual void CalcOmega(double& omega, const int itnum);

    //! old scatra increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnpold_;
  };


}  // namespace SSI

#endif /* SSI_PARTITIONED2_WC_H_ */
