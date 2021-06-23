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

namespace INPAR
{
  namespace SSI
  {
    enum class FieldCoupling;
  }

}  // namespace INPAR

namespace LINALG
{
  class MultiMapExtractor;
}

namespace SCATRA
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace SCATRA

namespace SSI
{
  // forward declaration
  class SSICouplingBase;

  namespace UTILS
  {
    class SSIMeshTyingMaps;
    class SSISlaveSideConverter;
    class SSIStructureMeshTying;
  }  // namespace UTILS

  enum class RedistributionType
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
    \return void
    \date 08/16
    \author rauch  */
    virtual void Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string& struct_disname, const std::string& scatra_disname, bool isAle) = 0;

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

    /*!
     * @brief checks whether simulation is restarted or not
     *
     * @return  flag indicating if simulation is restarted
     */
    bool IsRestart() const;

    /*! \brief Setup discretizations and dofsets

     Init coupling object \ref ssicoupling_ and
     other possible coupling objects in derived
     classes

    \return int
    \date 08/16
    \author vuong, rauch  */
    virtual RedistributionType InitFieldCoupling(const Epetra_Comm& comm,
        const std::string& struct_disname, const std::string& scatra_disname);

    /*! \brief Setup discretizations

    \date 08/16
    \author rauch  */
    virtual void InitDiscretizations(const Epetra_Comm& comm, const std::string& struct_disname,
        const std::string& scatra_disname);

    //! structural map extractor for coupling of two entities (0: other degrees of freedom, 1:
    //! slave-side degrees of freedom, 2: master-side degrees of freedom)
    Teuchos::RCP<const LINALG::MultiMapExtractor> MapsCoupStruct() const;

    //! Interior dofs of structure mesh tying
    Teuchos::RCP<const Epetra_Map> MapStructureInterior() const;

    //! Master side dofs of structure mesh tying
    Teuchos::RCP<const Epetra_Map> MapStructureMaster() const;

    //! Slave side dofs of structure mesh tying
    Teuchos::RCP<const Epetra_Map> MapStructureSlave() const;

    //! structural map extractor for 3 domains intersecting (0: other degrees of freedom, 1:
    //! slave-side degrees  of freedom, 2: master-side degrees of freedom)
    Teuchos::RCP<const LINALG::MultiMapExtractor> MapsCoupStruct3DomainIntersection() const;

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
    }

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

    //! access the scalar transport base algorithm on manifolds
    const Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>& ScaTraManifoldBaseAlgorithm() const
    {
      return scatra_manifold_base_algorithm_;
    }

    //! access the scalar transport field
    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() const;

    //! access the scalar transport field on manifolds
    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraManifold() const;

    /// set structure solution on other fields
    void SetStructSolution(
        Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel);

    /// set scatra solution on other fields
    virtual void SetScatraSolution(Teuchos::RCP<const Epetra_Vector> phi) const;

    /// set scatra on manifold solution on other fields
    void SetScatraManifoldSolution(Teuchos::RCP<const Epetra_Vector> phi);

    /// set temperature field  by evaluating time dependent function
    void EvaluateAndSetTemperatureField();

    //! get bool indicating if we have at least one ssi interface meshtying condition
    bool SSIInterfaceMeshtying() const { return ssiinterfacemeshtying_; }

    //! indicates, if meshtying includes triple points
    bool Meshtying3DomainIntersection() const { return meshtying_3_domain_intersection_; }

    //! return the scatra-scatra interface meshtying strategy
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingStrategyS2I() const
    {
      return meshtying_strategy_s2i_;
    }

    //! returns whether calculation of the initial potential field is performed
    bool DoCalculateInitialPotentialField() const;

    //! returns if the scalar transport time integration is of type electrochemistry
    bool IsElchScaTraTimIntType() const;

    //! solve additional scatra field on manifolds
    bool IsScaTraManifold() const { return is_scatra_manifold_; }

    //! Redistribute nodes and elements on processors
    void Redistribute(RedistributionType redistribution_type);

    //! get bool indicating if we have at least one ssi interface contact condition
    bool SSIInterfaceContact() const { return ssiinterfacecontact_; }

    //! SSI structure meshtying object containing coupling adapters, converters and maps
    Teuchos::RCP<SSI::UTILS::SSIStructureMeshTying> SSIStructureMeshTying() const
    {
      return ssi_structure_meshtying_;
    }

   protected:
    //! get bool indicating if old structural time integration is used
    bool UseOldStructureTimeInt() const { return use_old_structure_; }

    //! check if \ref Setup() was called
    void CheckIsSetup() const
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    }

    //! check if \ref Init() was called
    void CheckIsInit() const
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    }

    //! copy modified time step from scatra to scatra manifold field
    void SetDtFromScaTraToManifold();

    //! copy modified time step from scatra to this SSI algorithm
    void SetDtFromScaTraToSSI();

    //! copy modified time step from scatra to structure field
    void SetDtFromScaTraToStructure();

    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; }

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; }

    //! set up structural model evaluator for scalar-structure interaction
    virtual void SetupModelEvaluator() const = 0;

   private:
    /*!
     * @brief Checks whether flags for adaptive time stepping in ssi have been set consistently
     *
     * @param[in] scatraparams  parameter list containing the SCALAR TRANSPORT DYNAMIC parameters
     * @param[in] structparams  parameter list containing the STRUCTURAL DYNAMIC parameters
     */
    static void CheckAdaptiveTimeStepping(
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams);

    /*!
     * @brief Time integrators for the scalar and structure fields are instantiated and initialized
     *
     * @param[in] globaltimeparams  parameter list containing the SSI CONTROL parameters
     * @param[in] scatraparams      parameter list containing the SCALAR TRANSPORT DYNAMIC
     *                              parameters
     * @param[in] structparams      parameter list containing the STRUCTURAL DYNAMIC parameters
     * @param[in] struct_disname    name of structure discretization
     * @param[in] scatra_disname    name of scalar transport discretization
     * @param[in] isAle             flag indicating if ALE is activated
     */
    void InitTimeIntegrators(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string& struct_disname, const std::string& scatra_disname, const bool isAle);

    //! check whether scatra-structure interaction flags are set correctly
    void CheckSSIFlags() const;

    /*!
     * @brief SSI interface condition definition is checked
     *
     * @param[in] struct_disname  name of structure discretization
     */
    void CheckSSIInterfaceConditions(const std::string& struct_disname) const;

    //! returns true if Init(..) was called and is still valid
    bool IsInit() const { return isinit_; }

    /// set structure mesh displacement on scatra field
    void SetMeshDisp(Teuchos::RCP<const Epetra_Vector> disp);

    /// set structure velocity field on scatra field
    void SetVelocityFields(Teuchos::RCP<const Epetra_Vector> vel);

    //! Type of coupling strategy between the two fields of the SSI problems
    const INPAR::SSI::FieldCoupling fieldcoupling_;

    //! flag indicating if class is initialized
    bool isinit_;

    //! flag indicating if class is setup
    bool issetup_;

    //! solve additional scatra field on manifolds
    const bool is_scatra_manifold_;

    //! counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling
    //! iterations (partitioned algorithm)
    int iter_;

    //! indicates, if meshtying includes three entities
    const bool meshtying_3_domain_intersection_;

    //! meshtying strategy for scatra-scatra interface coupling on scatra discretization
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i_;

    //! underlying scatra problem base algorithm
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_base_algorithm_;

    //! underlying scatra problem base algorithm on manifolds
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_manifold_base_algorithm_;

    //! SSI structure mesh tying object containing coupling adapters, converters and maps
    Teuchos::RCP<SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying_;

    /// helper class for applying SSI couplings
    Teuchos::RCP<SSICouplingBase> ssicoupling_;

    //! bool indicating if we have at least one ssi interface contact condition
    const bool ssiinterfacecontact_;

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

    //! Flag to indicate whether old structural time integration is used.
    bool use_old_structure_;

    //! a zero vector of full length with structure dofs
    Teuchos::RCP<Epetra_Vector> zeros_structure_;
  };  // SSI_Base
}  // namespace SSI
#endif /* SSI_BASE_H_ */
