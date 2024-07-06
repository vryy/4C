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
namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace PoroElast
{
  class Monolithic;
}

namespace Adapter
{
  class Fluid;
  class AleFpsiWrapper;
}  // namespace Adapter

namespace Core::LinAlg
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace Core::LinAlg

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
    explicit FpsiCoupling(Teuchos::RCP<PoroElast::Monolithic> poro,
        Teuchos::RCP<Adapter::Fluid> fluid, Teuchos::RCP<Adapter::AleFpsiWrapper> ale,
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
    Core::LinAlg::BlockSparseMatrixBase& c_pp() { return *c_pp_; }
    // Fluid-Fluid Coupling Matrix
    Core::LinAlg::BlockSparseMatrixBase& c_ff()
    {
      return *c_ff_;
    }  // blockmatrix for condensation!!!
    // Poro-Fluid Coupling Matrix
    Core::LinAlg::BlockSparseMatrixBase& c_pf() { return *c_pf_; }
    // Fluid-Poro Coupling Matrix
    Core::LinAlg::BlockSparseMatrixBase& c_fp() { return *c_fp_; }
    // Poro-Ale Coupling Matrix
    Core::LinAlg::BlockSparseMatrixBase& c_pa() { return *c_pa_; }
    // Fluid-Ale Coupling Matrix
    Core::LinAlg::SparseMatrix& c_fa() { return *c_fa_; }

    //@}

    // Poro Coupling RHS (structure)
    Teuchos::RCP<Epetra_Vector>& rhs_s() { return c_rhs_s_; }
    // Poro Coupling RHS (fluid)
    Teuchos::RCP<Epetra_Vector>& rhs_pf() { return c_rhs_pf_; }
    // Fluid Coupling RHS
    Teuchos::RCP<Epetra_Vector>& rhs_f() { return c_rhs_f_; }

    //! @name transform helpers

    // Vector Transform
    Teuchos::RCP<Epetra_Vector> i_fluid_to_porofluid(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_pf_f_->slave_to_master(iv);
    }

    Teuchos::RCP<Epetra_Vector> i_porofluid_to_fluid(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_pf_f_->master_to_slave(iv);
    }

    Teuchos::RCP<Epetra_Vector> i_fluid_to_porostruct(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_f_->slave_to_master(iv);
    }

    Teuchos::RCP<Epetra_Vector> i_porostruct_to_fluid(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_f_->master_to_slave(iv);
    }

    Teuchos::RCP<Epetra_Vector> i_ale_to_porostruct(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_a_->slave_to_master(iv);
    }

    Teuchos::RCP<Epetra_Vector> i_porostruct_to_ale(Teuchos::RCP<const Epetra_Vector> iv) const
    {
      return icoup_ps_a_->master_to_slave(iv);
    }

    //@}

    //! @name access coupling objects

    Core::Adapter::Coupling& poro_fluid_fluid_coupling() { return *icoup_pf_f_; }

    Core::Adapter::Coupling& poro_structure_fluid_coupling() { return *icoup_ps_f_; }

    Core::Adapter::Coupling& poro_structure_ale_coupling() { return *icoup_ps_a_; }

    //@}

    //! @name access extractors

    const Teuchos::RCP<Core::LinAlg::MapExtractor>& fluid_fpsi_vel_pres_extractor() const
    {
      return fluidvelpres_extractor_;
    }
    const Teuchos::RCP<Core::LinAlg::MapExtractor>& fluid_fpsi_vel_extractor() const
    {
      return fluidvel_extractor_;
    }
    const Teuchos::RCP<Core::LinAlg::MapExtractor>& poro_fluid_fpsi_vel_pres_extractor() const
    {
      return porofluid_extractor_;
    }
    const Teuchos::RCP<Core::LinAlg::MultiMapExtractor>& poro_extractor() const
    {
      return poro_extractor_;
    }
    const Teuchos::RCP<FPSI::UTILS::MapExtractor>& fluid_fsi_fpsi_extractor() const
    {
      return fluid_fsifpsi_extractor_;
    }

    //@}

    // set hydraulic conductivity
    void set_conductivity(double conduct);

   private:
    // access to the fields
    const Teuchos::RCP<PoroElast::Monolithic>& poro_field() { return poro_; }
    const Teuchos::RCP<Adapter::Fluid>& fluid_field() { return fluid_; }
    const Teuchos::RCP<Adapter::AleFpsiWrapper>& ale_field() { return ale_; }

    // Initialize Coupling Matrixes and Coupling RHS
    void init_coupling_matrixes_rhs();

    // underlying poroelast problem
    Teuchos::RCP<PoroElast::Monolithic> poro_;
    // underlying fluid of the FPSI problem
    Teuchos::RCP<Adapter::Fluid> fluid_;
    // underlying ale of the FPSI problem
    Teuchos::RCP<Adapter::AleFpsiWrapper> ale_;

    // Poro-Poro Coupling Matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> c_pp_;
    // Fluid-Fluid Coupling Matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> c_ff_;
    // Poro-Fluid Coupling Matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> c_pf_;
    // Fluid-Poro Coupling Matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> c_fp_;
    // Poro-Ale Coupling Matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> c_pa_;
    // Fluid-Ale Coupling Matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        c_fa_;  // block matrix to cut out just ale other block (-->interface (fpsi & fsi) is
                // condensed to structural dofs!)

    // Poro Coupling RHS
    Teuchos::RCP<Epetra_Vector> c_rhs_s_;
    Teuchos::RCP<Epetra_Vector> c_rhs_pf_;
    // Fluid Coupling RHS
    Teuchos::RCP<Epetra_Vector> c_rhs_f_;

    // Interface Coupling PoroFluid - Fluid velocities and pressure are/is coupled
    Teuchos::RCP<Core::Adapter::Coupling> icoup_pf_f_;
    // Interface Coupling PoroStructure - Fluid
    Teuchos::RCP<Core::Adapter::Coupling> icoup_ps_f_;
    // Interface Coupling PoroStructure - Ale
    Teuchos::RCP<Core::Adapter::Coupling> icoup_ps_a_;

    // extractor for fpsi condition from fluid
    Teuchos::RCP<Core::LinAlg::MapExtractor> fluidvelpres_extractor_;
    Teuchos::RCP<Core::LinAlg::MapExtractor> fluidvel_extractor_;
    // extractor for fpsi condition from poro fluid
    Teuchos::RCP<Core::LinAlg::MapExtractor> porofluid_extractor_;
    // extractor for fpsi condition from (poro) structure
    Teuchos::RCP<Core::LinAlg::MapExtractor> porostruct_extractor_;
    //! dof row map splitted in inner structure (0), struct interface (1)
    //! inner porofluid (2) and porofluidinterface (3)
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> poro_extractor_;
    Teuchos::RCP<FPSI::UTILS::MapExtractor> fluid_fsifpsi_extractor_;

    // Evaluate is called first time!
    bool isfirstcall_;

    Teuchos::RCP<std::map<int, int>> fluid_poro_fluid_interface_map_;
    Teuchos::RCP<std::map<int, int>> poro_fluid_fluid_interface_map_;

    Teuchos::RCP<Core::LinAlg::MatrixRowTransform>
        couplingrowtransform_;  /// g_fpsi || F->PF transform (FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixRowTransform>
        couplingrowtransform2_;  /// g_fpsi || PF->F transform (FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixRowTransform>
        couplingrowtransform3_;  /// g_fpsi || PF->F transform (FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixRowTransform>
        couplingrowtransform4_;  /// g_fpsi || F->PS transform (FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixRowTransform> couplingrowtransform5_;
    ;  /// g_fpsi || F->PS transform (FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixColTransform>
        couplingcoltransform_;  /// for Row/Col-Map for Full - fluid_field & F->PS transform (FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixColTransform>
        couplingcoltransform2_;  /// for Row/Col-Map for Full - ale_field & A->PS transform (FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixRowColTransform>
        couplingrowcoltransform_;  /// g_fpsi/g_fpsi || F->PS/F->PS transform (FPSI/FPSI)
    Teuchos::RCP<Core::LinAlg::MatrixRowColTransform>
        couplingrowcoltransform2_;  /// g_fpsi/g_fpsi || F->PF/F->PS transform (FPSI/FPSI)

    // hydraulic conductivity (needed for coupling in case of probtype fps3i)
    double conductivity_;
  };  // fpsi_coupling
}  // namespace FPSI

FOUR_C_NAMESPACE_CLOSE

#endif
