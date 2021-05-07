/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef SSTI_UTILS_H_
#define SSTI_UTILS_H_

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

#include "../drt_ssi/ssi_clonestrategy.H"
#include "../drt_sti/sti_clonestrategy.H"

#include "../drt_adapter/adapter_coupling.H"

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

namespace LINALG
{
  class BlockSparseMatrixBase;
  enum class MatrixType;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

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

    virtual ~SSTIMaps() = default;

    //! get maps of subproblems
    //@{
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsScatra() const { return maps_scatra_; };
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsStructure() const
    {
      return maps_structure_;
    };
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsThermo() const { return maps_thermo_; };
    //@}

    //! global map extractor (0: scalar transport, 1: structure, 2: thermo)
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsSubproblems() const
    {
      return maps_subproblems_;
    };

    //! return map with dofs on both sides of interface
    const Teuchos::RCP<Epetra_Map> MapInterface(
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy) const;

    //! return block map with dofs on both sides of interface
    const Teuchos::RCP<LINALG::MultiMapExtractor> MapsInterfaceBlocks(
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
        LINALG::MatrixType scatramatrixtype, unsigned nummaps) const;

    //! return block map with dofs on slave side of interface
    const Teuchos::RCP<LINALG::MultiMapExtractor> MapsInterfaceBlocksSlave(
        Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategy,
        LINALG::MatrixType scatramatrixtype, unsigned nummaps) const;

   private:
    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_structure_;

    //! global map extractor (0: scalar transport, 1: structure, 2: thermo)
    Teuchos::RCP<LINALG::MultiMapExtractor> maps_subproblems_;

    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_thermo_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //! holds all maps in context of SSTI monolithic simulations
  class SSTIMapsMono : public SSTIMaps
  {
   public:
    SSTIMapsMono(const SSTI::SSTIMono& ssti_mono_algorithm);

    ~SSTIMapsMono() override = default;

    //! map extractor associated with blocks of global system matrix
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsSystemMatrixSubblocks() const
    {
      return maps_systemmatrix_subblocks_;
    };

   private:
    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_systemmatrix_subblocks_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //! sets up and holds all sub blocks of system matrices and system matrix for SSTI simulations
  class SSTIMatrices
  {
   public:
    SSTIMatrices(Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono,
        const LINALG::MatrixType matrixtype_global, const LINALG::MatrixType matrixtype_scatra,
        Teuchos::RCP<Epetra_Map> interface_map_scatra,
        Teuchos::RCP<Epetra_Map> interface_map_thermo,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapscatrainterface,
        Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermointerface, bool isinterfacemeshtying);

    //! method that clears all ssi matrices
    void ClearMatrices();

    //! call complete on the scalar transport - structure off-diagonal matrices
    void CompleteScaTraStructureMatrices();

    Teuchos::RCP<LINALG::SparseOperator> SystemMatrix() { return systemmatrix_; };

    //! return sub blocks of system matrix
    //@{
    Teuchos::RCP<LINALG::SparseOperator> ScaTraStructureDomain() { return scatrastructuredomain_; };
    Teuchos::RCP<LINALG::SparseOperator> ScaTraStructureInterface()
    {
      return scatrastructureinterface_;
    };
    Teuchos::RCP<LINALG::SparseOperator> ScaTraThermoDomain() { return scatrathermodomain_; };
    Teuchos::RCP<LINALG::SparseOperator> ScaTraThermoInterface() { return scatrathermointerface_; };
    Teuchos::RCP<LINALG::SparseOperator> StructureScaTraDomain() { return structurescatradomain_; };
    Teuchos::RCP<LINALG::SparseOperator> StructureThermoDomain() { return structurethermodomain_; };
    Teuchos::RCP<LINALG::SparseOperator> ThermoScaTraDomain() { return thermoscatradomain_; };
    Teuchos::RCP<LINALG::SparseOperator> ThermoScaTraInterface() { return thermoscatrainterface_; };
    Teuchos::RCP<LINALG::SparseOperator> ThermoStructureDomain() { return thermostructuredomain_; };
    Teuchos::RCP<LINALG::SparseOperator> ThermoStructureInterface()
    {
      return thermostructureinterface_;
    };
    //@}

   private:
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SetupBlockMatrix(
        Teuchos::RCP<const LINALG::MultiMapExtractor> row_map,
        Teuchos::RCP<const LINALG::MultiMapExtractor> col_map);

    Teuchos::RCP<LINALG::SparseMatrix> SetupSparseMatrix(
        const Teuchos::RCP<const Epetra_Map> row_map);

    //! scalar transport interface map
    Teuchos::RCP<Epetra_Map> interface_map_scatra_;

    //! scalar transport matrix type
    const LINALG::MatrixType matrixtype_scatra_;

    //! maps for monolithic treatment of scalar transport-structure-thermo-interaction
    Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono_;

    Teuchos::RCP<LINALG::SparseOperator> systemmatrix_;
    //! subblocks of system matrix
    //@{
    Teuchos::RCP<LINALG::SparseOperator> scatrastructuredomain_;
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface_;
    Teuchos::RCP<LINALG::SparseOperator> scatrathermodomain_;
    Teuchos::RCP<LINALG::SparseOperator> scatrathermointerface_;
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain_;
    Teuchos::RCP<LINALG::SparseOperator> structurethermodomain_;
    Teuchos::RCP<LINALG::SparseOperator> thermoscatradomain_;
    Teuchos::RCP<LINALG::SparseOperator> thermoscatrainterface_;
    Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain_;
    Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface_;
    //@}
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
    std::map<std::string, std::string> ConditionsToCopy() override;

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
    std::map<std::string, std::string> ConditionsToCopy() override;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  class SSTIStructuralMeshtying
  {
   public:
    SSTIStructuralMeshtying();

    void Setup(const SSTI::SSTIAlgorithm& ssti_algorithm);

    //! coupling adapters for structral meshtying
    //@{
    const Teuchos::RCP<ADAPTER::Coupling>& InterfaceCouplingAdapterStructure() const
    {
      return icoup_structure_;
    }
    const Teuchos::RCP<ADAPTER::Coupling>& InterfaceCouplingAdapterStructure3DomainIntersection()
        const
    {
      return icoup_structure_3_domain_intersection_;
    }
    //@}

    //! Slave structure converter for structure coupling
    //@{
    const Teuchos::RCP<ADAPTER::CouplingSlaveConverter>
    InterfaceCouplingAdapterStructureSlaveConverter() const
    {
      return icoup_structure_slave_converter_;
    }
    const Teuchos::RCP<ADAPTER::CouplingSlaveConverter>
    InterfaceCouplingAdapterStructureSlaveConverter3DomainIntersection() const
    {
      return icoup_structure_slave_converter_3_domain_intersection_;
    }
    //@}

    //! return map for interior and master-side structural degrees of freedom
    const Teuchos::RCP<const Epetra_Map> MapStructureCondensed() const
    {
      return map_structure_condensed_;
    }

    //! structural map extractor (0: slave-side degrees of freedom, 1: master-side degrees of
    //! freedom, 3: others[unused])
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsInterfaceStructure() const
    {
      return maps_interface_structure_;
    }

    //! structural map extractor (0: slave-side degrees of freedom, 1: master-side degrees of
    //! freedom
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsInterfaceStructure3DomainIntersection()
        const
    {
      return maps_interface_structure_3_domain_intersection_;
    }

    //! indicates, if meshtying includes triple points
    bool Meshtying3DomainIntersection() const { return meshtying_3_domain_intersection_; };

   private:
    //! coupling adapters for structral meshtying
    //@{
    Teuchos::RCP<ADAPTER::Coupling> icoup_structure_;
    Teuchos::RCP<ADAPTER::Coupling> icoup_structure_3_domain_intersection_;
    //@}

    //! Slave structure converter for structure coupling
    //@{
    Teuchos::RCP<ADAPTER::CouplingSlaveConverter> icoup_structure_slave_converter_;
    Teuchos::RCP<ADAPTER::CouplingSlaveConverter>
        icoup_structure_slave_converter_3_domain_intersection_;
    //@}

    //! map for interior and master-side structural degrees of freedom
    Teuchos::RCP<const Epetra_Map> map_structure_condensed_;

    //! structural map extractor (0: slave-side degrees of freedom, 1: master-side degrees of
    //! freedom, 3: others[unused])
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_interface_structure_;

    //! structural map extractor (0: slave-side degrees of freedom, 1: master-side degrees of
    //! freedom
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_interface_structure_3_domain_intersection_;

    //! indicates, if meshtying includes three entities
    const bool meshtying_3_domain_intersection_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  void PrintSSTILogo(int pid);
}  // namespace SSTI

#endif /* SSTI_UTILS_H_ */
