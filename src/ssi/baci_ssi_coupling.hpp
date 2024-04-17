/*----------------------------------------------------------------------*/
/*! \file
 \brief helper classes for  scalar-structure coupling

 \level 3


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SSI_COUPLING_HPP
#define FOUR_C_SSI_COUPLING_HPP

#include "baci_config.hpp"

#include "baci_coupling_adapter_mortar.hpp"
#include "baci_coupling_adapter_volmortar.hpp"
#include "baci_lib_discret.hpp"
#include "baci_scatra_timint_implicit.hpp"
#include "baci_ssi_base.hpp"

FOUR_C_NAMESPACE_OPEN


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

    //! \brief init this class
    //!
    //! \param ndim                        dimension of the problem
    //! \param structdis                   underlying structure discretization
    //! \param ssi_base                    underlying scatra-structure time integrator
    virtual void Init(const int ndim, Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<SSI::SSIBase> ssi_base) = 0;

    //! \brief setup this class
    virtual void Setup() = 0;


    //! \brief exchange material pointers of both discratizations
    //!
    //! \param structdis   underlying structure discretization
    //! \param scatradis   underlying scatra discretization
    virtual void AssignMaterialPointers(Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> scatradis) = 0;

    //!
    //! \param scatradis      underlying scatra discretization
    //! \param stress_state   mechanical stress state vector to set
    //! \param nds            number of dofset to write state on
    virtual void SetMechanicalStressState(DRT::Discretization& scatradis,
        Teuchos::RCP<const Epetra_Vector> stress_state, unsigned nds) = 0;

    //! \brief set structure mesh displacement on other field
    //!
    //! \param scatra    underlying scatra problem of the SSI problem
    //! \param disp      displacement field to set
    virtual void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> disp) = 0;

    //! \brief set structure velocity fields on other field
    //!
    //! \param scatra    underlying scatra problem of the SSI problem
    //! \param convvel   convective velocity field to set
    //! \param vel       velocity field to set
    virtual void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> convvel, Teuchos::RCP<const Epetra_Vector> vel) = 0;

    //! \brief set scatra solution on other field
    //!
    //! \param dis    discretization to write scatra solution on
    //! \param phi    scalar field solution
    //! \param nds    number of dofset to write state on
    virtual void SetScalarField(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) = 0;

    //! \brief set micro soultion of scatra field other field
    //!
    //! \param dis     discretization to write micro scatra solution on
    //! \param phi     micro scatra solution
    //! \param nds     number of dofset to write micro scatra solution on
    virtual void SetScalarFieldMicro(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) = 0;

    //! set temperature field on structure field
    virtual void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) = 0;
  };

  //! solid-scatra coupling for matching volume meshes
  class SSICouplingMatchingVolume : public SSICouplingBase
  {
   public:
    SSICouplingMatchingVolume() : issetup_(false), isinit_(false){};
    void Init(const int ndim, Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<SSI::SSIBase> ssi_base) override;

    void Setup() override;

    void AssignMaterialPointers(Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> scatradis) override;

    void SetMechanicalStressState(DRT::Discretization& scatradis,
        Teuchos::RCP<const Epetra_Vector> stress_statetemp, unsigned nds) override;

    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> disp) override;

    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> convvel, Teuchos::RCP<const Epetra_Vector> vel) override;

    void SetScalarField(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

    void SetScalarFieldMicro(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

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
    void Init(const int ndim, Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<SSI::SSIBase> ssi_base) override;

    void Setup() override;

    void AssignMaterialPointers(Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> scatradis) override;

    void SetMechanicalStressState(DRT::Discretization& scatradis,
        Teuchos::RCP<const Epetra_Vector> stress_state, unsigned nds) override
    {
      dserror("only implemented for 'SSICouplingMatchingVolume'");
    }

    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> disp) override;

    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> convvel, Teuchos::RCP<const Epetra_Vector> vel) override;

    void SetScalarField(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

    void SetScalarFieldMicro(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

    void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) override
    {
      dserror("only for matching nodes");
    };

   private:
    //! adapter to mortar framework
    Teuchos::RCP<CORE::ADAPTER::CouplingMortar> adaptermeshtying_;

    //! extractor for coupled surface of structure discretization with surface scatra
    Teuchos::RCP<CORE::LINALG::MapExtractor> extractor_;

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
    void Init(const int ndim, Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<SSI::SSIBase> ssi_base) override;

    void Setup() override;

    void AssignMaterialPointers(Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> scatradis) override;

    void SetMechanicalStressState(DRT::Discretization& scatradis,
        Teuchos::RCP<const Epetra_Vector> stress_state, unsigned nds) override
    {
      dserror("only implemented for 'SSICouplingMatchingVolume'");
    }

    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> disp) override;

    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> convvel, Teuchos::RCP<const Epetra_Vector> vel) override;

    void SetScalarField(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

    void SetScalarFieldMicro(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

    void SetTemperatureField(
        DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp) override
    {
      dserror("only for matching nodes");
    };

   private:
    //! volume coupling (using mortar) adapter
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> volcoupl_structurescatra_;

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
    void Init(const int ndim, Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<SSI::SSIBase> ssi_base) override;

    void Setup() override;


    void AssignMaterialPointers(Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> scatradis) override;

    void SetMechanicalStressState(DRT::Discretization& scatradis,
        Teuchos::RCP<const Epetra_Vector> stress_state, unsigned nds) override
    {
      dserror("only implemented for 'SSICouplingMatchingVolume'");
    }

    void SetMeshDisp(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> disp) override;

    void SetVelocityFields(Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
        Teuchos::RCP<const Epetra_Vector> convvel, Teuchos::RCP<const Epetra_Vector> vel) override;

    void SetScalarField(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

    void SetScalarFieldMicro(
        DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds) override;

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


FOUR_C_NAMESPACE_CLOSE

#endif
