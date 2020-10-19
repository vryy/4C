/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure thermo algorithms

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef SSTI_ALGORITHM_H_
#define SSTI_ALGORITHM_H_

#include "../drt_adapter/adapter_algorithmbase.H"

#include <Epetra_Vector.h>

// forward declarations
namespace ADAPTER
{
  class Coupling;
  class ScaTraBaseAlgorithm;
  class SSIStructureWrapper;
  class Structure;
  class StructureBaseAlgorithmNew;
}  // namespace ADAPTER

namespace INPAR
{
  namespace SSTI
  {
    enum class SolutionScheme;
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

namespace SSTI
{
  //! Base class of all solid-scatra algorithms
  class SSTIAlgorithm : public ADAPTER::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit SSTIAlgorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& sstiparams);

    /// virtual destructor
    ~SSTIAlgorithm() override = default;

    //! Setup of algorithm
    //! Clone Discretizations, init and setup subproblems, setup coupling adapters at interfaces,
    //! setup submatrices for coupling between fields
    //@{
    virtual void Init(const Epetra_Comm& comm, const Teuchos::ParameterList& sstitimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams,
        const Teuchos::ParameterList& structparams) = 0;
    virtual void Setup();
    virtual void SetupSystem() = 0;
    //@}

    //! return coupling
    //@{
    const Teuchos::RCP<ADAPTER::Coupling> CouplingAdapterStructure() const
    {
      return icoup_structure_;
    };
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingScatra() const
    {
      return meshtying_strategy_scatra_;
    };
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingThermo() const
    {
      return meshtying_strategy_thermo_;
    };
    //@}

    //! return subproblems
    //@{
    const Teuchos::RCP<::ADAPTER::SSIStructureWrapper> StructureField() const
    {
      return structure_;
    };
    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() const;
    const Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ThermoField() const;
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> ScaTraFieldBase() { return scatra_; };
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> ThermoFieldBase() { return thermo_; };
    //@}

    //! get bool indicating if we have at least one ssi interface meshtying condition
    bool InterfaceMeshtying() const { return interfacemeshtying_; };

    //! read restart
    void ReadRestart(int restart) override;

    //! timeloop of coupled problem
    virtual void Timeloop() = 0;

    //! test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm) const;

   protected:
    //! clone scatra from structure and then thermo from scatra
    virtual void CloneDiscretizations(const Epetra_Comm& comm);

    //! copies modified time step from scatra to structure and to this SSI algorithm
    void DistributeDtFromScaTra();

    //! copy states between subproblems
    //@{
    void SetSolutionAllFields();
    void SetScatraSolution();
    void SetStructureSolution();
    void SetThermoSolution();
    //@}

    //! counter for Newton-Raphson iterations (monolithic algorithm) or outer coupling iterations
    //! (partitioned algorithm)
    unsigned int iter_;

   private:
    //! excahnge materials between discretizations
    void AssignMaterialPointers();

    void CheckIsInit();

    //! clone thermo parameters from scatra parameters and adjust where needed
    Teuchos::ParameterList CloneThermoParams(
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams);

    //! Setup meshtying for structure
    void SetupCouplingAdapterStructure();

    //! Pointers to subproblems
    //@{
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure_;
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> struct_adapterbase_ptr_;
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo_;
    //@}

    //! Pointers to coupling strategies
    //@{
    Teuchos::RCP<ADAPTER::Coupling> icoup_structure_;
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra_;
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo_;
    //@}

    //! bool indicating if we have at least one ssi interface meshtying condition
    const bool interfacemeshtying_;

    //! flag indicating if class is initialized
    bool isinit_;

    //! flag indicating if class is setup
    bool issetup_;
  };  // SSTI_Algorithm

  //! Construct specific SSTI algorithm
  Teuchos::RCP<SSTI::SSTIAlgorithm> BuildSSTI(INPAR::SSTI::SolutionScheme coupling,
      const Epetra_Comm& comm, const Teuchos::ParameterList& sstiparams);
}  // namespace SSTI
#endif /* SSTI_ALGORITHM_H_ */
