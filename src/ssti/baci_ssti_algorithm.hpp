/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all scalar structure thermo algorithms

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SSTI_ALGORITHM_HPP
#define FOUR_C_SSTI_ALGORITHM_HPP

#include "baci_config.hpp"

#include "baci_adapter_algorithmbase.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

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

namespace CORE::LINALG
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
  namespace UTILS
  {
    class SSIMeshTying;
  }
}  // namespace SSI

namespace SSTI
{
  //! Base class of all solid-scatra algorithms
  class SSTIAlgorithm : public ADAPTER::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    explicit SSTIAlgorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

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

    //! increment the counter for Newton-Raphson iterations (monolithic algorithm)
    void IncrementIter() { ++iter_; }

    //! return the counter for Newton-Raphson iterations (monolithic algorithm)
    unsigned int Iter() const { return iter_; }

    //! reset the counter for Newton-Raphson iterations (monolithic algorithm)
    void ResetIter() { iter_ = 0; }

    //! return coupling
    //@{
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingScatra() const
    {
      return meshtying_strategy_scatra_;
    }
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> MeshtyingThermo() const
    {
      return meshtying_strategy_thermo_;
    }
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> SSTIStructureMeshTying() const
    {
      return ssti_structure_meshtying_;
    }
    //@}

    //! return subproblems
    //@{
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> StructureField() const { return structure_; };
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() const;
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ThermoField() const;
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

    //! distribute states between subproblems
    //@{
    void DistributeSolutionAllFields();
    void DistributeScatraSolution();
    void DistributeStructureSolution();
    void DistributeThermoSolution();
    //@}

   private:
    //! counter for Newton-Raphson iterations (monolithic algorithm)
    unsigned int iter_;

    //! exchange materials between discretizations
    void AssignMaterialPointers();

    void CheckIsInit();

    //! clone thermo parameters from scatra parameters and adjust where needed
    Teuchos::ParameterList CloneThermoParams(
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& thermoparams);

    //! Pointers to subproblems
    //@{
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure_;
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> struct_adapterbase_ptr_;
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo_;
    //@}

    //! Pointers to coupling strategies
    //@{
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra_;
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo_;
    Teuchos::RCP<SSI::UTILS::SSIMeshTying> ssti_structure_meshtying_;
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
FOUR_C_NAMESPACE_CLOSE

#endif
