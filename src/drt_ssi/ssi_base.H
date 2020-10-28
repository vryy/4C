/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure algorithms

 \level 1


 *------------------------------------------------------------------------------------------------*/

#ifndef SSI_BASE_H_
#define SSI_BASE_H_

#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_lib/drt_discret.H"

#include <Epetra_Vector.h>

#include "../drt_inpar/inpar_ssi.H"

// forward declarations
namespace ADAPTER
{
  class Coupling;
  class CouplingSlaveConverter;
  class Structure;
  class ScaTraBaseAlgorithm;
  class SSIStructureWrapper;
  class StructureBaseAlgorithmNew;
}  // namespace ADAPTER

namespace LINALG
{
  class MultiMapExtractor;
}

namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace SSI
{
  // forward declaration
  class SSICouplingBase;

  enum RedistributionType
  {
    none,     //!< unknown redistribution type
    binning,  //!< redistribute by binning
    match     //!< redistribute by node matching
  };

  //! Base class of all solid-scatra algorithms
  class SSIBase : public ADAPTER::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit SSIBase(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    ~SSIBase() override = default;

    //! return scatra-scatra interface coupling adapter for structure field
    const Teuchos::RCP<ADAPTER::Coupling>& InterfaceCouplingAdapterStructure() const
    {
      return icoup_structure_;
    };

    //! return scatra-scatra interface coupling adapter slave converter for structure field
    ADAPTER::CouplingSlaveConverter& InterfaceCouplingAdapterStructureSlaveConverter() const
    {
      return *icoup_structure_slave_converter_;
    }

    //! return counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm)
    int IterationCount() const { return iter_; }

    //! reset the counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm)
    void ResetIterationCount() { iter_ = 0; }

    //! increment the counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm) by 1
    void IncrementIterationCount() { iter_ += 1; }

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
        const std::string struct_disname, const std::string scatra_disname, bool isAle) = 0;

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

    //! returns true if Setup() was called and is still valid
    bool IsSetup() const { return issetup_; };

    /*! \brief Setup discretizations and dofsets

     Init coupling object \ref ssicoupling_ and
     other possible coupling objects in derived
     classes

    \return int
    \date 08/16
    \author vuong, rauch  */
    virtual int InitFieldCoupling(const Epetra_Comm& comm, const std::string& struct_disname,
        const std::string& scatra_disname);

    /*! \brief Setup discretizations

    \date 08/16
    \author rauch  */
    virtual void InitDiscretizations(const Epetra_Comm& comm, const std::string& struct_disname,
        const std::string& scatra_disname);

    //! return structural map extractor (0: other degrees of freedom, 1: slave-side degrees of
    //! freedom, 2: master-side degrees of freedom)
    const Teuchos::RCP<const LINALG::MultiMapExtractor>& MapsStructure() const
    {
      return maps_structure_;
    };

    //! return map for interior and master-side structural degrees of freedom
    const Teuchos::RCP<const Epetra_Map>& MapStructureCondensed() const
    {
      return map_structure_condensed_;
    };

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm) const;

    /// read restart
    void ReadRestart(int restart) override;

    /// read restart from given time
    void ReadRestartfromTime(double restarttime) override;

    //! access to structural field
    const Teuchos::RCP<::ADAPTER::SSIStructureWrapper>& StructureField() const
    {
      return structure_;
    };

    /// pointer to the underlying structure problem base algorithm
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> StructureBaseAlgorithm() const
    {
      return struct_adapterbase_ptr_;
    }

    //! access the scalar transport base algorithm
    const Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>& ScaTraBaseAlgorithm() const
    {
      return scatra_base_algorithm_;
    }

    //! access the scalar transport field
    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() const;

    /// set structure solution on scatra field
    void SetStructSolution(
        Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel);

    /// set scatra solution on structure field
    void SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi);

    /// set temperature field  by evaluating time dependent function
    void EvaluateAndSetTemperatureField();

    //! get bool indicating if we have at least one ssi interface meshtying condition
    bool SSIInterfaceMeshtying() const { return ssiinterfacemeshtying_; };

   protected:
    //! get bool indicating if old structural time integration is used
    bool UseOldStructureTimeInt() const { return use_old_structure_; }

    //! check if \ref Setup() was called
    void CheckIsSetup() const
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! check if \ref Init() was called
    void CheckIsInit() const
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

    //! copies modified time step from scatra to structure and to this SSI algorithm
    void SetDtFromScaTraToStructure();

    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

    //! set up structural model evaluator for scalar-structure interaction
    virtual void SetupModelEvaluator() const = 0;

   private:
    //! returns true if Init(..) was called and is still valid
    bool IsInit() const { return isinit_; }

    /// set structure mesh displacement on scatra field
    void SetMeshDisp(Teuchos::RCP<const Epetra_Vector> disp);

    /// set structure velocity field on scatra field
    void SetVelocityFields(Teuchos::RCP<const Epetra_Vector> vel);

    //! Type of coupling strategy between the two fields of the SSI problems
    const INPAR::SSI::FieldCoupling fieldcoupling_;

    //! scatra-scatra interface coupling adapter for structure field
    Teuchos::RCP<ADAPTER::Coupling> icoup_structure_;

    //! scatra-scatra interface coupling adapter for structure field
    Teuchos::RCP<ADAPTER::CouplingSlaveConverter> icoup_structure_slave_converter_;

    //! flag indicating if class is initialized
    bool isinit_;

    //! flag indicating if class is setup
    bool issetup_;

    //! counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling iterations
    //! (partitioned algorithm)
    int iter_;

    //! structural map extractor (0: other degrees of freedom, 1: slave-side degrees of freedom, 2:
    //! master-side degrees of freedom)
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_structure_;

    //! map for interior and master-side structural degrees of freedom
    Teuchos::RCP<const Epetra_Map> map_structure_condensed_;

    /// ptr underlying scatra problem base algorithm
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_base_algorithm_;

    /// helper class for applying SSI couplings
    Teuchos::RCP<SSICouplingBase> ssicoupling_;

    //! bool indicating if we have at least one ssi interface meshtying condition
    const bool ssiinterfacemeshtying_;

    /// ptr to underlying structure
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure_;

    /// ptr to the underlying structure problem base algorithm
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> struct_adapterbase_ptr_;

    //! number of function for prescribed temperature
    const int temperature_funct_num_;

    //! vector of temperatures
    Teuchos::RCP<Epetra_Vector> temperature_vector_;

    //! Temporary flag to indicate whether old strutural time integration is used.
    //! todo Remove as soon as structural elements are adapted to new structure.
    bool use_old_structure_;

    /// a zero vector of full length
    Teuchos::RCP<Epetra_Vector> zeros_;

  };  // SSI_Base
}  // namespace SSI
#endif /* SSI_BASE_H_ */
