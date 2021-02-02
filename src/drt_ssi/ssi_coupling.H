/*----------------------------------------------------------------------*/
/*! \file
 \brief helper classes for  scalar-structure coupling

 \level 3


 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SSI_SSI_COUPLING_H_
#define SRC_DRT_SSI_SSI_COUPLING_H_

#include "../drt_lib/drt_discret.H"


namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class Structure;
  class ScaTraTimIntImpl;
  class CouplingMortar;
  class MortarVolCoupl;
}  // namespace ADAPTER

namespace SSI
{
  //! Base class of solid-scatra coupling helper classes
  class SSICouplingBase
  {
   public:
    SSICouplingBase() = default;

    virtual ~SSICouplingBase() = default;

    //! Init
    virtual void Init(const int ndim,                 //!< dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis,  //!< underlying scatra discretization
        Teuchos::RCP<DRT::Discretization> scatra_manifold_dis) = 0;

    //! \brief setup this class
    virtual void Setup() = 0;


    //! exchange material pointers of both discratizations
    virtual void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis   //!< underlying scatra discretization
        ) = 0;

    //! set structure mesh displacement on scatra field
    virtual void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                                 scatra,        //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp  //!< displacement field to set
        ) = 0;

    //! set structure velocity fields on scatra field
    virtual void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                                       scatra,  //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel       //!< velocity field to set
        ) = 0;

    //! set scatra solution on structure field
    virtual void SetScalarField(DRT::Discretization& structdis,  //!< structural discretization
        Teuchos::RCP<const Epetra_Vector> phi                    //!< scalar field
        ) = 0;

    //! set temperature field on structure field
    virtual void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) = 0;
  };

  //! solid-scatra coupling for matching volume meshes
  class SSICouplingMatchingVolume : public SSICouplingBase
  {
   public:
    //! constructor
    SSICouplingMatchingVolume() : issetup_(false), isinit_(false){};

    //! destructor
    ~SSICouplingMatchingVolume() override = default;

    //! Init
    void Init(const int ndim,                         //!< dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis,  //!< underlying scatra discretization
        Teuchos::RCP<DRT::Discretization> scatra_manifold_dis) override;

    //! \brief setup this class
    void Setup() override;


    //! exchange material pointers of both discretizations
    void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis   //!< underlying scatra discretization
        ) override;

    //! set structure mesh displacement on scatra field
    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                         scatra,                //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp  //!< displacement field to set
        ) override;

    //! set structure velocity fields on scatra field
    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                               scatra,  //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel       //!< velocity field to set
        ) override;

    //! set scatra solution on structure field
    void SetScalarField(DRT::Discretization& structdis,  //!< structural discretization
        Teuchos::RCP<const Epetra_Vector> phi            //!< scalar field
        ) override;

    void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) override;

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() { return isinit_; };

    //! returns true if class was setup and setup is still valid
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! returns true if class was init and init is still valid
    void CheckIsInit()
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };
  };

  //! solid-scatra coupling for matching boundary meshes
  class SSICouplingNonMatchingBoundary : public SSICouplingBase
  {
   public:
    SSICouplingNonMatchingBoundary()
        : adaptermeshtying_(Teuchos::null),
          extractor_(Teuchos::null),
          issetup_(false),
          isinit_(false){};

    ~SSICouplingNonMatchingBoundary() override = default;

    //! Setup
    void Init(const int ndim,                         //!< dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis,  //!< underlying scatra discretization
        Teuchos::RCP<DRT::Discretization> scatra_manifold_dis) override;

    //! \brief setup this class
    void Setup() override;

    //! exchange material pointers of both discretizations
    void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis   //!< underlying scatra discretization
        ) override;

    //! set structure mesh displacement on scatra field
    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                         scatra,                //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp  //!< displacement field to set
        ) override;

    //! set structure velocity field on scatra field
    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                               scatra,  //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel       //!< velocity field to set
        ) override;

    //! set scatra solution on structure field
    void SetScalarField(DRT::Discretization& structdis,  //!< structural discretization
        Teuchos::RCP<const Epetra_Vector> phi            //!< scalar field
        ) override;

    void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) override
    {
      dserror("only for matching nodes");
    };

   private:
    //! adapter to mortar framework
    Teuchos::RCP<ADAPTER::CouplingMortar> adaptermeshtying_;

    //! extractor for coupled surface of structure discretization with surface scatra
    Teuchos::RCP<LINALG::MapExtractor> extractor_;

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

    //! spatial dimension of the global problem
    int problem_dimension_;

    //! pointer to structdis_
    Teuchos::RCP<DRT::Discretization> structdis_;

    //! pointer to scatradis_
    Teuchos::RCP<DRT::Discretization> scatradis_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() { return isinit_; };

    //! returns true if class was setup and setup is still valid
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! returns true if class was init and init is still valid
    void CheckIsInit()
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };
  };

  //! solid-scatra coupling for non-matching boundary meshes
  class SSICouplingNonMatchingVolume : public SSICouplingBase
  {
   public:
    SSICouplingNonMatchingVolume()
        : volcoupl_structurescatra_(Teuchos::null), issetup_(false), isinit_(false){};

    ~SSICouplingNonMatchingVolume() override = default;

    //! Init
    void Init(const int ndim,                         //!< dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis,  //!< underlying scatra discretization
        Teuchos::RCP<DRT::Discretization> scatra_manifold_dis) override;

    //! \brief setup this class
    void Setup() override;

    //! exchange material pointers of both discretizations
    void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis   //!< underlying scatra discretization
        ) override;

    //! set structure mesh displacement on scatra field
    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                         scatra,                //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp  //!< displacement field to set
        ) override;

    //! set structure velocity field on scatra field
    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                               scatra,  //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel       //!< velocity field to set
        ) override;

    //! set scatra solution on structure field
    void SetScalarField(DRT::Discretization& structdis,  //!< structural discretization
        Teuchos::RCP<const Epetra_Vector> phi            //!< scalar field
        ) override;

    void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) override
    {
      dserror("only for matching nodes");
    };

   private:
    //! volume coupling (using mortar) adapter
    Teuchos::RCP<::ADAPTER::MortarVolCoupl> volcoupl_structurescatra_;

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() { return isinit_; };

    //! returns true if class was setup and setup is still valid
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! returns true if class was init and init is still valid
    void CheckIsInit()
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };
  };

  //! solid-scatra coupling for matching volume and boundary meshes
  class SSICouplingMatchingVolumeAndBoundary : public SSICouplingBase
  {
   public:
    SSICouplingMatchingVolumeAndBoundary() : issetup_(false), isinit_(false){};

    ~SSICouplingMatchingVolumeAndBoundary() override = default;

    //! \brief Setup
    //! \note  We need to make sure that the parallel distribution of Volume and Boundary
    //!        is the same externally! The best thing is if you do this in your *_dyn.cpp,
    //!        i.e., your global control algorithm.
    void Init(const int ndim,                         //!< dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis,  //!< underlying scatra discretization
        Teuchos::RCP<DRT::Discretization> scatra_manifold_dis) override;

    //! \brief setup this class
    void Setup() override;

    //! exchange material pointers of both discretizations
    void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,  //!< underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis   //!< underlying scatra discretization
        ) override;

    //! set structure mesh displacement on scatra field
    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                         scatra,                //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp  //!< displacement field to set
        ) override;

    //! set structure velocity field on scatra field
    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>
                               scatra,  //!< underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel       //!< velocity field to set
        ) override;

    //! set scatra solution on structure field
    void SetScalarField(DRT::Discretization& structdis,  //!< structural discretization
        Teuchos::RCP<const Epetra_Vector> phi            //!< scalar field
        ) override;

    void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) override;

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() { return isinit_; };

    //! returns true if class was setup and setup is still valid
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! returns true if class was init and init is still valid
    void CheckIsInit()
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

   public:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };
  };
}  // namespace SSI


#endif /* SRC_DRT_SSI_SSI_COUPLING_H_ */
