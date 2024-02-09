/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef BACI_SSTI_UTILS_HPP
#define BACI_SSTI_UTILS_HPP

#include "baci_config.hpp"

#include "baci_coupling_adapter.hpp"
#include "baci_ssi_clonestrategy.hpp"
#include "baci_sti_clonestrategy.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class Coupling;
  class SSIStructureWrapper;
}  // namespace ADAPTER

namespace INPAR
{
  namespace SCATRA
  {
    enum class MatrixType;
  }
}  // namespace INPAR

namespace CORE::LINALG
{
  class BlockSparseMatrixBase;
  enum class MatrixType;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace SCATRA

namespace SSTI
{
  class SSTIAlgorithm;
  class SSTIMono;

  //! holds all maps in context of SSTI simulations
  class SSTIMaps
  {
   public:
    SSTIMaps(const SSTI::SSTIMono& ssti_mono_algorithm);

    //! get maps of subproblems
    //@{
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapScatra() const
    {
      return block_map_scatra_;
    }
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapStructure() const
    {
      return block_map_structure_;
    }
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapThermo() const
    {
      return block_map_thermo_;
    }
    //@}

    /*!
     * @brief global map extractor
     * @note only access with GetProblemPosition method
     */
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> MapsSubProblems() const
    {
      return maps_subproblems_;
    }

    //! return map with dofs on both sides of interface
    Teuchos::RCP<Epetra_Map> MapInterface(
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy) const;

    //! return block map with dofs on both sides of interface
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> MapsInterfaceBlocks(
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
        CORE::LINALG::MatrixType scatramatrixtype, unsigned nummaps) const;

    //! return block map with dofs on slave side of interface
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> MapsInterfaceBlocksSlave(
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
        CORE::LINALG::MatrixType scatramatrixtype, unsigned nummaps) const;

   private:
    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure_;

    //! map extractor associated with all degrees of freedom inside thermo field
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo_;

    //! global map extractor (0: scalar transport, 1: structure, 2: thermo)
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> maps_subproblems_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //! holds all maps in context of SSTI monolithic simulations
  class SSTIMapsMono : public SSTIMaps
  {
   public:
    SSTIMapsMono(const SSTI::SSTIMono& ssti_mono_algorithm);

    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMapSystemMatrix() const
    {
      return block_map_system_matrix_;
    };

   private:
    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_system_matrix_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //! sets up and holds all sub blocks of system matrices and system matrix for SSTI simulations
  class SSTIMatrices
  {
   public:
    SSTIMatrices(Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono,
        const CORE::LINALG::MatrixType matrixtype_global,
        const CORE::LINALG::MatrixType matrixtype_scatra, bool interfacemeshtying);

    //! method that clears all ssi matrices
    void ClearMatrices();

    //! call complete on all coupling matrices
    void CompleteCouplingMatrices();

    //! call uncomplete on all coupling matrices
    void UnCompleteCouplingMatrices();

    Teuchos::RCP<CORE::LINALG::SparseOperator> SystemMatrix() { return systemmatrix_; };

    //! return sub blocks of system matrix
    //@{
    Teuchos::RCP<CORE::LINALG::SparseOperator> ScaTraStructureDomain()
    {
      return scatrastructuredomain_;
    };
    Teuchos::RCP<CORE::LINALG::SparseOperator> ScaTraStructureInterface()
    {
      return scatrastructureinterface_;
    };
    Teuchos::RCP<CORE::LINALG::SparseOperator> ScaTraThermoDomain() { return scatrathermodomain_; };
    Teuchos::RCP<CORE::LINALG::SparseOperator> ScaTraThermoInterface()
    {
      return scatrathermointerface_;
    };
    Teuchos::RCP<CORE::LINALG::SparseOperator> StructureScaTraDomain()
    {
      return structurescatradomain_;
    };
    Teuchos::RCP<CORE::LINALG::SparseOperator> StructureThermoDomain()
    {
      return structurethermodomain_;
    };
    Teuchos::RCP<CORE::LINALG::SparseOperator> ThermoScaTraDomain() { return thermoscatradomain_; };
    Teuchos::RCP<CORE::LINALG::SparseOperator> ThermoScaTraInterface()
    {
      return thermoscatrainterface_;
    };
    Teuchos::RCP<CORE::LINALG::SparseOperator> ThermoStructureDomain()
    {
      return thermostructuredomain_;
    };
    Teuchos::RCP<CORE::LINALG::SparseOperator> ThermoStructureInterface()
    {
      return thermostructureinterface_;
    };
    //@}

   private:
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> SetupBlockMatrix(
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> row_map,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> col_map);

    Teuchos::RCP<CORE::LINALG::SparseMatrix> SetupSparseMatrix(
        const Teuchos::RCP<const Epetra_Map> row_map);

    //! scalar transport matrix type
    const CORE::LINALG::MatrixType matrixtype_scatra_;

    //! maps for monolithic treatment of scalar transport-structure-thermo-interaction
    Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono_;

    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix_;
    //! subblocks of system matrix
    //@{
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrastructuredomain_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrastructureinterface_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrathermodomain_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrathermointerface_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> structurescatradomain_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> structurethermodomain_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> thermoscatradomain_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> thermoscatrainterface_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> thermostructuredomain_;
    Teuchos::RCP<CORE::LINALG::SparseOperator> thermostructureinterface_;
    //@}

    //! bool indicating if we have at least one ssi interface meshtying condition
    const bool interfacemeshtying_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  class ConvCheckMono
  {
   public:
    ConvCheckMono(const Teuchos::ParameterList params);

    //! Is this Newton step converged
    bool Converged(const SSTI::SSTIMono& ssti_mono);

   private:
    //! maximum number of Newton-Raphson iteration steps
    const unsigned itermax_;

    //! relative tolerance for Newton-Raphson iteration
    const double itertol_;

    //! absolute tolerance for residual vectors
    const double restol_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  class SSTIScatraStructureCloneStrategy : public SSI::ScatraStructureCloneStrategy
  {
   public:
    /// returns condition names to be copied (source and target name)
    std::map<std::string, std::string> ConditionsToCopy() const override;

   protected:
    //! provide cloned element with element specific data (material etc.)
    void SetElementData(
        Teuchos::RCP<DRT::Element> newele,  //! current cloned element on target discretization
        DRT::Element* oldele,               //! current element on source discretization
        const int matid,                    //! material of cloned element
        const bool isnurbs                  //! nurbs flag
        ) override;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  class SSTIScatraThermoCloneStrategy : public STI::ScatraThermoCloneStrategy
  {
   protected:
    /// returns condition names to be copied (source and target name)
    std::map<std::string, std::string> ConditionsToCopy() const override;
  };
}  // namespace SSTI

BACI_NAMESPACE_CLOSE

#endif  // SSTI_UTILS_H
