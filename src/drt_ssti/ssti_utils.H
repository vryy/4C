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
  enum class MatrixType;
  class MultiMapExtractor;
}  // namespace LINALG

namespace SCATRA
{
  class MeshtyingStrategyS2I;
  class ScaTraTimIntImpl;
}  // namespace SCATRA

namespace SSTI
{
  class SSTIMono;

  class SSTIMaps
  {
   public:
    SSTIMaps(Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> thermo,
        Teuchos::RCP<ADAPTER::Coupling> icoupstructure, bool interfacemeshyting);

    virtual ~SSTIMaps() = default;

    //! getmaps of subproblems
    //@{
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsScatra() const { return maps_scatra_; };
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsStructure() const
    {
      return maps_structure_;
    };
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsThermo() const { return maps_thermo_; };
    //@}

    //! return map for interior and master-side structural degrees of freedom
    const Teuchos::RCP<const Epetra_Map> MapStructureCondensed() const
    {
      return map_structure_condensed_;
    };

    //! global map extractor (0: scalar transport, 1: structure, 2: thermo)
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsSubproblems() const
    {
      return maps_subproblems_;
    };

    //! structural map extractor (0: slave-side degrees of freedom, 1: master-side degrees of
    //! freedom, 3: others[unused])
    const Teuchos::RCP<const LINALG::MultiMapExtractor> MapsInterfaceStructure() const
    {
      return maps_interface_structure_;
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
    //! map for interior and master-side structural degrees of freedom
    Teuchos::RCP<const Epetra_Map> map_structure_condensed_;

    //! structural map extractor (0: slave-side degrees of freedom, 1: master-side degrees of
    //! freedom, 3: others[unused])
    Teuchos::RCP<const LINALG::MultiMapExtractor> maps_interface_structure_;

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
  class SSTIMapsMono : public SSTIMaps
  {
   public:
    SSTIMapsMono(Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
        Teuchos::RCP<SCATRA::ScaTraTimIntImpl> thermo,
        Teuchos::RCP<ADAPTER::Coupling> icoupstructure, bool interfacemeshyting);

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
   protected:
    /// returns condition names to be copied (source and target name)
    std::map<std::string, std::string> ConditionsToCopy() override;
  };  // class SSTIScatraCloneStrategy

  /*---------------------------------------------------------------------------------*
   *---------------------------------------------------------------------------------*/
  void PrintSSTILogo(int pid);
}  // namespace SSTI

#endif /* SSTI_UTILS_H_ */
