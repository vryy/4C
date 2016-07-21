/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned_2wc.H

 \brief two way coupled partitioned scalar structure interaction

\maintainer    Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_PARTITIONED2_WC_H_
#define SSI_PARTITIONED2_WC_H_

#include "ssi_partitioned.H"
#include "../drt_inpar/inpar_ssi.H"

namespace SSI
{
  //! base class to deal with partioned 2WC SSI. Mainly it is the same as for thermo-structure interaction (TSI)
  class SSI_Part2WC : public SSI_Part
  {
  public:

    //!explicit constructor
    explicit SSI_Part2WC(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //!virtual destructor
    virtual ~SSI_Part2WC(){};

    //!full time loop
    virtual void Timeloop();

    //! perform iteration loop between fields
    virtual void OuterLoop();

    //! prepare time step for single fields
    virtual void PrepareTimeStep(bool printheader = true);

    // return current iteration step
    int Itnum(){return itnum_;};

    //! update time step and print to screen
    virtual void UpdateAndOutput();

  protected:

    //! perform iteration step of structure field and set the new disp and vel states in the scatra field
    virtual void DoStructStep();

    //! perform iteration step of scatra field and set the new phi state in the structure field
    virtual void DoScatraStep();

    //! update the current states in every iteration
    //! states are set to the last solutions obtained
    virtual void IterUpdateStates();

    //! convergence check of outer loop
    virtual bool ConvergenceCheck(int itnum);

    /// velocity calculation given the displacements
    Teuchos::RCP<Epetra_Vector> CalcVelocity( Teuchos::RCP<const Epetra_Vector> dispnp );

    //! scalar increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnp_;
    //! displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

    //! convergence tolerance
    double ittol_;
    //! maximum iteration steps
    int itmax_;

  private:

    //! current iteration step
    int itnum_;

  };

  //! class to deal with displacement relaxated partioned 2WC SSI. Relaxation parameter is constant
  class SSI_Part2WC_SolidToScatra_Relax : public SSI_Part2WC
  {
  public:

    //!explicit constructor
    explicit SSI_Part2WC_SolidToScatra_Relax(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //!virtual destructor
    virtual ~SSI_Part2WC_SolidToScatra_Relax(){};

  protected:
    //! perform iteration loop between fields with relaxed displacements
    void OuterLoop();

    //! //! calculate relaxation parameter
    virtual void CalcOmega( double& omega , const int itnum );

    //! relaxation parameter
    double omega_;
  };

  //! class to deal with displacement relaxated 2WC partioned SSI. Relaxation parameter is calculated via Aitken
  class SSI_Part2WC_SolidToScatra_Relax_Aitken : public SSI_Part2WC_SolidToScatra_Relax
  {
  public:

    //!explicit constructor
    explicit SSI_Part2WC_SolidToScatra_Relax_Aitken(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //!virtual destructor
    virtual ~SSI_Part2WC_SolidToScatra_Relax_Aitken(){};

  protected:

    //! Calculate relaxation parameter via Aitken
    virtual void CalcOmega( double& omega , const int itnum );

    //! old displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnpold_;

  };

  //! class to deal with scalar relaxated 2WC partioned SSI. Relaxation parameter is constant
  class SSI_Part2WC_ScatraToSolid_Relax : public SSI_Part2WC
  {
  public:

    //!explicit constructor
    explicit SSI_Part2WC_ScatraToSolid_Relax(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //!virtual destructor
    virtual ~SSI_Part2WC_ScatraToSolid_Relax(){};

  protected:
    //! perform iteration loop between fields with relaxed scalar
    void OuterLoop();

    //! calculate relaxation parameter
    virtual void CalcOmega( double& omega , const int itnum );

    //! relaxation parameter
    double omega_;
  };

  //! class to deal with scalar relaxated 2WC partioned SSI. Relaxation parameter is calculated via Aitken
  class SSI_Part2WC_ScatraToSolid_Relax_Aitken : public SSI_Part2WC_ScatraToSolid_Relax
  {
  public:

    //!explicit constructor
    explicit SSI_Part2WC_ScatraToSolid_Relax_Aitken(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //!virtual destructor
    virtual ~SSI_Part2WC_ScatraToSolid_Relax_Aitken(){};

  protected:

    //! Calculate relaxation parameter via Aitken
    virtual void CalcOmega( double& omega , const int itnum );

    //! old scatra increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnpold_;

  };


}

#endif /* SSI_PARTITIONED2_WC_H_ */
