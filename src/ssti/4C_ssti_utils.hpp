/*----------------------------------------------------------------------*/
/*! \file
 \brief Utility methods for SSI

 \level 2


 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SSTI_UTILS_HPP
#define FOUR_C_SSTI_UTILS_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_ssi_clonestrategy.hpp"
#include "4C_sti_clonestrategy.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class Coupling;
  class SSIStructureWrapper;
}  // namespace Adapter

namespace Inpar
{
  namespace ScaTra
  {
    enum class MatrixType;
  }
}  // namespace Inpar

namespace Core::LinAlg
{
  class BlockSparseMatrixBase;
  enum class MatrixType;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;
}  // namespace Core::LinAlg

namespace ScaTra
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace ScaTra

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
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_scatra() const
    {
      return block_map_scatra_;
    }
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_structure() const
    {
      return block_map_structure_;
    }
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo() const
    {
      return block_map_thermo_;
    }
    //@}

    /*!
     * @brief global map extractor
     * @note only access with GetProblemPosition method
     */
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> maps_sub_problems() const
    {
      return maps_subproblems_;
    }

    //! return map with dofs on both sides of interface
    Teuchos::RCP<Epetra_Map> map_interface(
        Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtyingstrategy) const;

    //! return block map with dofs on both sides of interface
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> maps_interface_blocks(
        Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtyingstrategy,
        Core::LinAlg::MatrixType scatramatrixtype, unsigned nummaps) const;

    //! return block map with dofs on slave side of interface
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> maps_interface_blocks_slave(
        Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtyingstrategy,
        Core::LinAlg::MatrixType scatramatrixtype, unsigned nummaps) const;

   private:
    //! map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_scatra_;

    //! map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_structure_;

    //! map extractor associated with all degrees of freedom inside thermo field
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo_;

    //! global map extractor (0: scalar transport, 1: structure, 2: thermo)
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> maps_subproblems_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //! holds all maps in context of SSTI monolithic simulations
  class SSTIMapsMono : public SSTIMaps
  {
   public:
    SSTIMapsMono(const SSTI::SSTIMono& ssti_mono_algorithm);

    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_system_matrix() const
    {
      return block_map_system_matrix_;
    };

   private:
    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_system_matrix_;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  //! sets up and holds all sub blocks of system matrices and system matrix for SSTI simulations
  class SSTIMatrices
  {
   public:
    SSTIMatrices(Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono,
        const Core::LinAlg::MatrixType matrixtype_global,
        const Core::LinAlg::MatrixType matrixtype_scatra, bool interfacemeshtying);

    //! method that clears all ssi matrices
    void clear_matrices();

    //! call complete on all coupling matrices
    void complete_coupling_matrices();

    //! call uncomplete on all coupling matrices
    void un_complete_coupling_matrices();

    Teuchos::RCP<Core::LinAlg::SparseOperator> system_matrix() { return systemmatrix_; };

    //! return sub blocks of system matrix
    //@{
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_structure_domain()
    {
      return scatrastructuredomain_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_structure_interface()
    {
      return scatrastructureinterface_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_thermo_domain()
    {
      return scatrathermodomain_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatra_thermo_interface()
    {
      return scatrathermointerface_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> structure_scatra_domain()
    {
      return structurescatradomain_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> structure_thermo_domain()
    {
      return structurethermodomain_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermo_scatra_domain()
    {
      return thermoscatradomain_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermo_scatra_interface()
    {
      return thermoscatrainterface_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermo_structure_domain()
    {
      return thermostructuredomain_;
    };
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermo_structure_interface()
    {
      return thermostructureinterface_;
    };
    //@}

   private:
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> setup_block_matrix(
        Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> row_map,
        Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> col_map);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> setup_sparse_matrix(
        const Teuchos::RCP<const Epetra_Map> row_map);

    //! scalar transport matrix type
    const Core::LinAlg::MatrixType matrixtype_scatra_;

    //! maps for monolithic treatment of scalar transport-structure-thermo-interaction
    Teuchos::RCP<SSTI::SSTIMapsMono> ssti_maps_mono_;

    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix_;
    //! subblocks of system matrix
    //@{
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrastructuredomain_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrastructureinterface_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermodomain_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermointerface_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> structurescatradomain_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> structurethermodomain_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermoscatradomain_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermoscatrainterface_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermostructuredomain_;
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermostructureinterface_;
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
    bool converged(const SSTI::SSTIMono& ssti_mono);

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
    std::map<std::string, std::string> conditions_to_copy() const override;

   protected:
    //! provide cloned element with element specific data (material etc.)
    void set_element_data(Teuchos::RCP<Core::Elements::Element>
                              newele,     //! current cloned element on target discretization
        Core::Elements::Element* oldele,  //! current element on source discretization
        const int matid,                  //! material of cloned element
        const bool isnurbs                //! nurbs flag
        ) override;
  };

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  class SSTIScatraThermoCloneStrategy : public STI::ScatraThermoCloneStrategy
  {
   protected:
    /// returns condition names to be copied (source and target name)
    std::map<std::string, std::string> conditions_to_copy() const override;
  };
}  // namespace SSTI

FOUR_C_NAMESPACE_CLOSE

#endif
