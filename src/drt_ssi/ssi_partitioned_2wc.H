/*----------------------------------------------------------------------*/
/*! \file
 \brief two way coupled partitioned scalar structure interaction

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED2_WC_H_
#define SSI_PARTITIONED2_WC_H_

#include "ssi_partitioned.H"
#include "../drt_inpar/inpar_ssi.H"

namespace SSI
{
  //! base class to deal with partioned 2WC SSI. Mainly it is the same as for thermo-structure
  //! interaction (TSI)
  class SSIPart2WC : public SSIPart
  {
   public:
    //! constructor
    SSIPart2WC(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart2WC() override = default;

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
    int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle) override;

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
    void Setup() override;

    /*! \brief Initialize discretizations

    \date 08/16
    \author rauch  */
    void InitDiscretizations(const Epetra_Comm& comm, const std::string& struct_disname,
        const std::string& scatra_disname) override;

    //! full time loop
    void Timeloop() override;

    //! perform iteration loop between fields
    virtual void OuterLoop();

    //! prepare time loop
    virtual void PrepareTimeLoop();

    //! prepare time step for single fields
    void PrepareTimeStep(bool printheader = true) override;

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
    void DoStructStep() override;

    //! perform iteration step of scatra field and set the new phi state in the structure field
    void DoScatraStep() override;

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
  class SSIPart2WCSolidToScatraRelax : public SSIPart2WC
  {
   public:
    //! constructor
    SSIPart2WCSolidToScatraRelax(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart2WCSolidToScatraRelax() override = default;

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
    int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle) override;

   protected:
    //! perform iteration loop between fields with relaxed displacements
    void OuterLoop() override;

    //! //! calculate relaxation parameter
    virtual void CalcOmega(double& omega, const int itnum);

    //! relaxation parameter
    double omega_;
  };

  //! class to deal with displacement relaxated 2WC partioned SSI. Relaxation parameter is
  //! calculated via Aitken
  class SSIPart2WCSolidToScatraRelaxAitken : public SSIPart2WCSolidToScatraRelax
  {
   public:
    //! constructor
    SSIPart2WCSolidToScatraRelaxAitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart2WCSolidToScatraRelaxAitken() override = default;

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
    void Setup() override;

   protected:
    //! Calculate relaxation parameter via Aitken
    void CalcOmega(double& omega, const int itnum) override;

    //! old displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnpold_;
  };

  //! class to deal with scalar relaxated 2WC partioned SSI. Relaxation parameter is constant
  class SSIPart2WCScatraToSolidRelax : public SSIPart2WC
  {
   public:
    //! constructor
    SSIPart2WCScatraToSolidRelax(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart2WCScatraToSolidRelax() override = default;

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
    int Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string struct_disname, const std::string scatra_disname, bool isAle) override;

   protected:
    //! perform iteration loop between fields with relaxed scalar
    void OuterLoop() override;

    //! calculate relaxation parameter
    virtual void CalcOmega(double& omega, const int itnum);

    //! relaxation parameter
    double omega_;
  };

  //! class to deal with scalar relaxated 2WC partioned SSI. Relaxation parameter is calculated via
  //! Aitken
  class SSIPart2WCScatraToSolidRelaxAitken : public SSIPart2WCScatraToSolidRelax
  {
   public:
    //! constructor
    SSIPart2WCScatraToSolidRelaxAitken(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    //! virtual destructor
    ~SSIPart2WCScatraToSolidRelaxAitken() override = default;

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
    void Setup() override;

   protected:
    //! Calculate relaxation parameter via Aitken
    void CalcOmega(double& omega, const int itnum) override;

    //! old scatra increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnpold_;
  };
}  // namespace SSI

#endif /* SSI_PARTITIONED2_WC_H_ */
