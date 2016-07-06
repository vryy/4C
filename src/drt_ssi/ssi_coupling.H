/*----------------------------------------------------------------------*/
/*!
 \file ssi_coupling.H

 \brief helper classes for  scalar-structure coupling

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
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
}

namespace SSI
{
  //! Base class of solid-scatra coupling helper classes
  class SSICouplingBase
  {
  public:

    //!constructor
    explicit SSICouplingBase(){};

    //!destructor
    virtual ~SSICouplingBase(){};

    //! Setup
    virtual void Setup(
        const int                               ndim,    /// dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        )=0;

    //! exchange material pointers of both discratizations
    virtual void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        )=0;

    /// set structure mesh displacement on scatra field
    virtual void SetMeshDisp(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
        )=0;

    /// set structure velocity fields on scatra field
    virtual void SetVelocityFields(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
        )=0;

    /// set scatra solution on structure field
    virtual void SetScalarField(
        Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
        Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
        ) =0;
  };

  //! solid-scatra coupling for matching volume meshes
  class SSICouplingMatchingVolume : public SSICouplingBase
  {
  public:

    //!constructor
    explicit SSICouplingMatchingVolume(){};

    //!destructor
    virtual ~SSICouplingMatchingVolume(){};

    //! Setup
    virtual void Setup(
        const int                              ndim,    // dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    //! exchange material pointers of both discretizations
    virtual void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    /// set structure mesh displacement on scatra field
    virtual void SetMeshDisp(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
        );

    /// set structure velocity fields on scatra field
    virtual void SetVelocityFields(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
        );

    /// set scatra solution on structure field
    virtual void SetScalarField(
        Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
        Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
        );
  };

  //! solid-scatra coupling for matching boundary meshes
  class SSICouplingNonMatchingBoundary : public SSICouplingBase
  {
  public:

    //!constructor
    explicit SSICouplingNonMatchingBoundary():
    adaptermeshtying_(Teuchos::null),
    extractor_(Teuchos::null)
    {};

    //!destructor
    virtual ~SSICouplingNonMatchingBoundary(){};

    //! Setup
    virtual void Setup(
        const int                         ndim,          // dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    //! exchange material pointers of both discretizations
    virtual void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    /// set structure mesh displacement on scatra field
    virtual void SetMeshDisp(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
        );

    /// set structure velocity field on scatra field
    virtual void SetVelocityFields(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
        );

    /// set scatra solution on structure field
    virtual void SetScalarField(
        Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
        Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
        );

  private:

    //! adapter to mortar framework
    Teuchos::RCP<ADAPTER::CouplingMortar> adaptermeshtying_;

    //! extractor for coupled surface of structure discretization with surface scatra
    Teuchos::RCP<LINALG::MapExtractor> extractor_;

  };

  //! solid-scatra coupling for non-matching boundary meshes
  class SSICouplingNonMatchingVolume : public SSICouplingBase
  {
  public:

    //!constructor
    explicit SSICouplingNonMatchingVolume():volcoupl_structurescatra_(Teuchos::null)
    {};

    //!destructor
    virtual ~SSICouplingNonMatchingVolume(){};

    //! Setup
    virtual void Setup(
        const int                         ndim,          /// dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    //! exchange material pointers of both discretizations
    virtual void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    /// set structure mesh displacement on scatra field
    virtual void SetMeshDisp(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
        );

    /// set structure velocity field on scatra field
    virtual void SetVelocityFields(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
        );

    /// set scatra solution on structure field
    virtual void SetScalarField(
        Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
        Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
        );

  private:
    //! volume coupling (using mortar) adapter
    Teuchos::RCP< ::ADAPTER::MortarVolCoupl> volcoupl_structurescatra_;
  };

  //! solid-scatra coupling for matching volume and boundary meshes
  class SSICouplingMatchingVolumeAndBoundary : public SSICouplingBase
  {
  public:

    //!constructor
    explicit SSICouplingMatchingVolumeAndBoundary()
    {};

    //!destructor
    virtual ~SSICouplingMatchingVolumeAndBoundary(){};

    //! Setup
    virtual void Setup(
        const int                         ndim,          /// dimension of the problem
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    //! exchange material pointers of both discretizations
    virtual void AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        );

    /// set structure mesh displacement on scatra field
    virtual void SetMeshDisp(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
        );

    /// set structure velocity field on scatra field
    virtual void SetVelocityFields(
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
        Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
        Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
        );

    /// set scatra solution on structure field
    virtual void SetScalarField(
        Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
        Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
        );

  };
}


#endif /* SRC_DRT_SSI_SSI_COUPLING_H_ */
