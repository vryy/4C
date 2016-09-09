/*!------------------------------------------------------------------------------------------------*
 \file ssi_base.H

 \brief base class for all scalar structure algorithms

 \level 1

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15264

 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_BASE_H_
#define SSI_BASE_H_

#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_lib/drt_discret.H"

#include <Epetra_Vector.h>

#include "../drt_inpar/inpar_ssi.H"

namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class Structure;
  class ScaTraTimIntImpl;
}

namespace SSI
{
  // forward declaration
  class SSICouplingBase;

  //! Base class of all solid-scatra algorithms
  class SSI_Base : public ADAPTER::AlgorithmBase
  {
  public:

    /// create using a Epetra_Comm
    explicit SSI_Base(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~SSI_Base(){};

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
    \return bool
    \date 08/16
    \author rauch  */
    virtual bool Init(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname) = 0;

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

    /*! \brief Setup discretizations and dofsets

     Init coupling object \ref ssicoupling_ and
     other possible coupling objects in derived
     classes

    \return void
    \date 08/16
    \author vuong, rauch  */
    virtual bool InitFieldCoupling(
        const Epetra_Comm& comm,
        const std::string& struct_disname,
        const std::string& scatra_disname);

    /*! \brief Setup discretizations

    \date 08/16
    \author rauch  */
    virtual void InitDiscretizations(
        const Epetra_Comm& comm,
        const std::string& struct_disname,
        const std::string& scatra_disname);

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    void TestResults(const Epetra_Comm& comm);

    /// read restart
    void ReadRestart( int restart );

    /// read restart from given time
    void ReadRestartfromTime( double restarttime );

    //! access to structural field
    const Teuchos::RCP<ADAPTER::Structure>& StructureField() const {return structure_;};

    //! access to scatra field
    const Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>&  ScaTraField() const {return scatra_;};

    /// set structure solution on scatra field
    void SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                            Teuchos::RCP<const Epetra_Vector> vel ) ;

    /// set scatra solution on structure field
    void SetScatraSolution( Teuchos::RCP<const Epetra_Vector> phi ) ;

  protected:

    /// interpolate structure quantity to a scatra one (e.g. via volmortar or meshtying)
    const Teuchos::RCP<const Epetra_Vector> StructureToScatra(const Teuchos::RCP<const Epetra_Vector> structurevector) const;

    /// interpolate scatra quantity to a structure one (e.g. via volmortar or meshtying)
    const Teuchos::RCP<const Epetra_Vector> ScatraToStructure(const Teuchos::RCP<const Epetra_Vector> scatravector) const;

    /// underlying structure of the SSI problem
    Teuchos::RCP< ::ADAPTER::Structure>    structure_;

    /// underlying scatra problem of the SSI problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

  private:

    ///Set material pointers for matching structure and scatra grids
    void SetMaterialPointersMatchingGrid( Teuchos::RCP<const DRT::Discretization> sourcedis, Teuchos::RCP<const DRT::Discretization> targetdis);

    /// a zero vector of full length
    Teuchos::RCP<Epetra_Vector> zeros_;

    //! Type of coupling strategy between the two fields of the SSI problems
    const INPAR::SSI::FieldCoupling fieldcoupling_;

  private:

    /// set structure mesh displacement on scatra field
    void SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp ) ;

    /// set structure velocity field on scatra field
    void SetVelocityFields( Teuchos::RCP<const Epetra_Vector> vel ) ;

    /// helper class for applying SSI couplings
    Teuchos::RCP< SSICouplingBase> ssicoupling_;

  private:

    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

  protected:

    //! returns true if Setup() was called and is still valid
    bool IsSetup()
    {return issetup_;};

    //! returns true if Init(..) was called and is still valid
    bool IsInit()
    {return isinit_;};

    //! check if \ref Setup() was called
    void CheckIsSetup()
    {if(not IsSetup()) dserror("Setup() was not called.");};

    //! check if \ref Init() was called
    void CheckIsInit()
    {if(not IsInit()) dserror("Init(...) was not called.");};

  public:

    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse){issetup_ = trueorfalse;};

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse){isinit_ = trueorfalse;};


  }; //SSI_Base


}

#endif /* SSI_BASE_H_ */
