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
    virtual void Setup(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname) = 0;

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

    /// setup adapters for transport on boundary if necessary
    void SetupFieldCoupling(const std::string struct_disname, const std::string scatra_disname);

    /// setup discretizations and dofsets
    void SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname, const std::string scatra_disname);

    ///Set material pointers for matching structure and scatra grids
    void SetMaterialPointersMatchingGrid( Teuchos::RCP<const DRT::Discretization> sourcedis, Teuchos::RCP<const DRT::Discretization> targetdis);

    /// a zero vector of full length
    Teuchos::RCP<Epetra_Vector> zeros_;

    //! Type of coupling strategy between the two fields of the SSI problems
    const INPAR::SSI::FieldCoupling fieldcoupling_;

    /// setup discretizations and dofsets
    void SetupDiscretizationsAndFieldCoupling(
        const Epetra_Comm& comm,
        const std::string& struct_disname,
        const std::string& scatra_disname);

  private:

    /// set structure mesh displacement on scatra field
    void SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp ) ;

    /// set structure velocity field on scatra field
    void SetVelocityFields( Teuchos::RCP<const Epetra_Vector> vel ) ;

    /// helper class for applying SSI couplings
    Teuchos::RCP< SSICouplingBase> ssicoupling_;


  }; //SSI_Base


}

#endif /* SSI_BASE_H_ */
