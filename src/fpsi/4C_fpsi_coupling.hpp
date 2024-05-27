/*----------------------------------------------------------------------*/
/*! \file

\brief \brief FPSI Coupling Object: Holds all objects on the Fluid-Poro-Interface and evaluates the
Fluid-Poro-Coupling Matrixes!


\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FPSI_COUPLING_HPP
#define FOUR_C_FPSI_COUPLING_HPP

// ALE includes
#include "4C_config.hpp"

#include "4C_ale_utils_mapextractor.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class SparseMatrix;
}

namespace POROELAST
{
  class Monolithic;
}

namespace ADAPTER
{
  class Fluid;
  class AleFpsiWrapper;
}  // namespace ADAPTER

namespace CORE::LINALG
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace CORE::LINALG

/*----------------------------------------------------------------------*/

namespace FPSI
{
  namespace UTILS
  {
    class MapExtractor;
  }

  class FpsiCoupling
  {
   public:
    // ctor
    explicit FpsiCoupling(Teuchos::RCP<POROELAST::Monolithic> poro,
        Teuchos::RCP<ADAPTER::Fluid> fluid, Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale,
        Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap,
        Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap);

    // Setup the Coupling Objects
    void setup_interface_coupling();

    // Method reinitializes the matrix transformation objects
    void re_init_coupling_matrix_transform();

    // Evaluate Coupling Matrixes and Coupling RHS
    void evaluate_coupling_matrixes_rhs();

    //! @name access coupling matrixes

    // Poro-Poro Coupling Matrix
    CORE::LINALG::BlockSparseMatrixBase& C_pp() { return *c_pp_; }
    // Fluid-Fluid Coupling Matrix
    CORE::LINALG::BlockSparseMatrixBase& C_ff()
    {
      return *c_ff_;
    }  // blockmatrix for condensation!!!
    // Poro-Fluid Coupling Matrix
    CORE::LINALG::BlockSparseMatrixBase& C_pf() { return *c_pf_; }
    // Fluid-Poro Coupling Matrix
    CORE::LINALG::BlockSparseMatrixBase& C_fp() { return *c_fp_; }
    // Poro-Ale Coupling Matrix
    CORE::LINALG::BlockSparseMatrixBase& C_pa() { return *c_pa_; }
    // Fluid-Ale Coupling Matrix
    CORE::LINALG::SparseMatrix& C_fa() { return *c_fa_; }

    //@}

    // Poro Coupling RHS (structure)
    Teuchos::RCP<Epetra_Vector>& RHS_s() { return c_rhs_s_; }
    // Poro Coupling RHS (fluid)
    Teuchos::RCP<Epetra_Vector>& RHS_pf() { return c_rhs_pf_; }
    // Fluid Coupling RHS
    Teuchos::RCP<Epetra_Vector>& RHS_f() { return c_rhs_f_; }

    //! @name transform helpers

    // Vector Transform
    Teuchos::RCP<Epetra_Vector> iFluidToPorofluid(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_pf_f_->SlaveToMaster(iv);
    }

    Teuchos::RCP<Epetra_Vector> iPorofluidToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_pf_f_->MasterToSlave(iv);
    }

    Teuchos::RCP<Epetra_Vector> iFluidToPorostruct(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_f_->SlaveToMaster(iv);
    }

    Teuchos::RCP<Epetra_Vector> iPorostructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_f_->MasterToSlave(iv);
    }

    Teuchos::RCP<Epetra_Vector> iAleToPorostruct(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_a_->SlaveToMaster(iv);
    }

    Teuchos::RCP<Epetra_Vector> iPorostructToAle(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_a_->MasterToSlave(iv);
    }

    //@}

    //! @name access coupling objects

    CORE::ADAPTER::Coupling& poro_fluid_fluid_coupling() { return *icoup_pf_f_; }

    CORE::ADAPTER::Coupling& poro_structure_fluid_coupling() { return *icoup_ps_f_; }

    CORE::ADAPTER::Coupling& poro_structure_ale_coupling() { return *icoup_ps_a_; }

    //@}

    //! @name access extractors

    const Teuchos::RCP<CORE::LINALG::MapExtractor>& fluid_fpsi_vel_pres_extractor() const
    {
      return fluidvelpres_extractor_;
    }
    const Teuchos::RCP<CORE::LINALG::MapExtractor>& fluid_fpsi_vel_extractor() const
    {
      return fluidvel_extractor_;
    }
    const Teuchos::RCP<CORE::LINALG::MapExtractor>& poro_fluid_fpsi_vel_pres_extractor() const
    {
      return porofluid_extractor_;
    }
    const Teuchos::RCP<CORE::LINALG::MultiMapExtractor>& PoroExtractor() const
    {
      return poro_extractor_;
    }
    const Teuchos::RCP<FPSI::UTILS::MapExtractor>& fluid_fsi_fpsi_extractor() const
    {
      return fluid_fsifpsi_extractor_;
    }

    //@}

    // set hydraulic conductivity
    void SetConductivity(double conduct);

   private:
    // access to the fields
    const Teuchos::RCP<POROELAST::Monolithic>& poro_field() { return poro_; }
    const Teuchos::RCP<ADAPTER::Fluid>& fluid_field() { return fluid_; }
    const Teuchos::RCP<ADAPTER::AleFpsiWrapper>& ale_field() { return ale_; }

    // Initialize Coupling Matrixes and Coupling RHS
    void init_coupling_matrixes_rhs();

    // underlying poroelast problem
    Teuchos::RCP<POROELAST::Monolithic> poro_;
    // underlying fluid of the FPSI problem
    Teuchos::RCP<ADAPTER::Fluid> fluid_;
    // underlying ale of the FPSI problem
    Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale_;

    // Poro-Poro Coupling Matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> c_pp_;
    // Fluid-Fluid Coupling Matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> c_ff_;
    // Poro-Fluid Coupling Matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> c_pf_;
    // Fluid-Poro Coupling Matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> c_fp_;
    // Poro-Ale Coupling Matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> c_pa_;
    // Fluid-Ale Coupling Matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        c_fa_;  // block matrix to cut out just ale other block (-->interface (fpsi & fsi) is
                // condensed to structural dofs!)

    // Poro Coupling RHS
    Teuchos::RCP<Epetra_Vector> c_rhs_s_;
    Teuchos::RCP<Epetra_Vector> c_rhs_pf_;
    // Fluid Coupling RHS
    Teuchos::RCP<Epetra_Vector> c_rhs_f_;

    // Interface Coupling PoroFluid - Fluid velocities and pressure are/is coupled
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoup_pf_f_;
    // Interface Coupling PoroStructure - Fluid
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoup_ps_f_;
    // Interface Coupling PoroStructure - Ale
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoup_ps_a_;

    // extractor for fpsi condition from fluid
    Teuchos::RCP<CORE::LINALG::MapExtractor> fluidvelpres_extractor_;
    Teuchos::RCP<CORE::LINALG::MapExtractor> fluidvel_extractor_;
    // extractor for fpsi condition from poro fluid
    Teuchos::RCP<CORE::LINALG::MapExtractor> porofluid_extractor_;
    // extractor for fpsi condition from (poro) structure
    Teuchos::RCP<CORE::LINALG::MapExtractor> porostruct_extractor_;
    //! dof row map splitted in inner structure (0), struct interface (1)
    //! inner porofluid (2) and porofluidinterface (3)
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> poro_extractor_;
    Teuchos::RCP<FPSI::UTILS::MapExtractor> fluid_fsifpsi_extractor_;

    // Evaluate is called first time!
    bool isfirstcall_;

    Teuchos::RCP<std::map<int, int>> fluid_poro_fluid_interface_map_;
    Teuchos::RCP<std::map<int, int>> poro_fluid_fluid_interface_map_;

    Teuchos::RCP<CORE::LINALG::MatrixRowTransform>
        couplingrowtransform_;  /// g_fpsi || F->PF transform (FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform>
        couplingrowtransform2_;  /// g_fpsi || PF->F transform (FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform>
        couplingrowtransform3_;  /// g_fpsi || PF->F transform (FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform>
        couplingrowtransform4_;  /// g_fpsi || F->PS transform (FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> couplingrowtransform5_;
    ;  /// g_fpsi || F->PS transform (FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixColTransform>
        couplingcoltransform_;  /// for Row/Col-Map for Full - fluid_field & F->PS transform (FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixColTransform>
        couplingcoltransform2_;  /// for Row/Col-Map for Full - ale_field & A->PS transform (FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform>
        couplingrowcoltransform_;  /// g_fpsi/g_fpsi || F->PS/F->PS transform (FPSI/FPSI)
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform>
        couplingrowcoltransform2_;  /// g_fpsi/g_fpsi || F->PF/F->PS transform (FPSI/FPSI)

    // hydraulic conductivity (needed for coupling in case of probtype fps3i)
    double conductivity_;
  };  // fpsi_coupling
}  // namespace FPSI

FOUR_C_NAMESPACE_CLOSE

#endif
