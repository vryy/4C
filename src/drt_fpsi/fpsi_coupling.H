/*----------------------------------------------------------------------*/
/*! \file

\brief \brief FPSI Coupling Object: Holds all objects on the Fluid-Poro-Interface and evaluates the
Fluid-Poro-Coupling Matrixes!

\maintainer Johannes Kremheller

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FPSI_COUPLING_H
#define FPSI_COUPLING_H

// ALE includes
#include "../drt_ale/ale_utils_mapextractor.H"

// ADAPTER includes
#include "../drt_adapter/adapter_coupling.H"

// LINALG includes
#include "../linalg/linalg_utils_densematrix_manipulation.H"

// FLUID includes
#include "../drt_fluid/fluid_utils_mapextractor.H"

// forward declarations
namespace LINALG
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

namespace FSI
{
  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

/*----------------------------------------------------------------------*/

namespace FPSI
{
  namespace UTILS
  {
    class MapExtractor;
  }

  class FPSICoupling
  {
   public:
    // ctor
    explicit FPSICoupling(Teuchos::RCP<POROELAST::Monolithic> poro,
        Teuchos::RCP<::ADAPTER::Fluid> fluid, Teuchos::RCP<::ADAPTER::AleFpsiWrapper> ale,
        Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap,
        Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap);

    // Setup the Coupling Objects
    void SetupInterfaceCoupling();

    // Method reinitializes the matrix transformation objects
    void ReInitCouplingMatrixTransform();

    // Evaluate Coupling Matrixes and Coupling RHS
    void EvaluateCouplingMatrixesRHS();

    //! @name access coupling matrixes

    // Poro-Poro Coupling Matrix
    LINALG::BlockSparseMatrixBase& C_pp() { return *c_pp_; }
    // Fluid-Fluid Coupling Matrix
    LINALG::BlockSparseMatrixBase& C_ff() { return *c_ff_; }  // blockmatrix for condensation!!!
    // Poro-Fluid Coupling Matrix
    LINALG::BlockSparseMatrixBase& C_pf() { return *c_pf_; }
    // Fluid-Poro Coupling Matrix
    LINALG::BlockSparseMatrixBase& C_fp() { return *c_fp_; }
    // Poro-Ale Coupling Matrix
    LINALG::BlockSparseMatrixBase& C_pa() { return *c_pa_; }
    // Fluid-Ale Coupling Matrix
    LINALG::SparseMatrix& C_fa() { return *c_fa_; }

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

    ADAPTER::Coupling& PoroFluidFluidCoupling() { return *icoup_pf_f_; }

    ADAPTER::Coupling& PoroStructureFluidCoupling() { return *icoup_ps_f_; }

    ADAPTER::Coupling& PoroStructureAleCoupling() { return *icoup_ps_a_; }

    //@}

    //! @name access extractors

    const Teuchos::RCP<LINALG::MapExtractor>& FluidFpsiVelPresExtractor() const
    {
      return fluidvelpres_extractor_;
    }
    const Teuchos::RCP<LINALG::MapExtractor>& FluidFpsiVelExtractor() const
    {
      return fluidvel_extractor_;
    }
    const Teuchos::RCP<LINALG::MapExtractor>& PoroFluidFpsiVelPresExtractor() const
    {
      return porofluid_extractor_;
    }
    const Teuchos::RCP<LINALG::MultiMapExtractor>& PoroExtractor() const { return poro_extractor_; }
    const Teuchos::RCP<FPSI::UTILS::MapExtractor>& FluidFsiFpsiExtractor() const
    {
      return fluid_fsifpsi_extractor_;
    }

    //@}

    // set hydraulic conductivity
    void SetConductivity(double conduct);

   private:
    // access to the fields
    const Teuchos::RCP<POROELAST::Monolithic>& PoroField() { return poro_; }
    const Teuchos::RCP<ADAPTER::Fluid>& FluidField() { return fluid_; }
    const Teuchos::RCP<ADAPTER::AleFpsiWrapper>& AleField() { return ale_; }

    // Initialize Coupling Matrixes and Coupling RHS
    void InitCouplingMatrixesRHS();

    // underlying poroelast problem
    Teuchos::RCP<POROELAST::Monolithic> poro_;
    // underlying fluid of the FPSI problem
    Teuchos::RCP<::ADAPTER::Fluid> fluid_;
    // underlying ale of the FPSI problem
    Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale_;

    // Poro-Poro Coupling Matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> c_pp_;
    // Fluid-Fluid Coupling Matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> c_ff_;
    // Poro-Fluid Coupling Matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> c_pf_;
    // Fluid-Poro Coupling Matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> c_fp_;
    // Poro-Ale Coupling Matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> c_pa_;
    // Fluid-Ale Coupling Matrix
    Teuchos::RCP<LINALG::SparseMatrix>
        c_fa_;  // block matrix to cut out just ale other block (-->interface (fpsi & fsi) is
                // condensed to structural dofs!)

    // Poro Coupling RHS
    Teuchos::RCP<Epetra_Vector> c_rhs_s_;
    Teuchos::RCP<Epetra_Vector> c_rhs_pf_;
    // Fluid Coupling RHS
    Teuchos::RCP<Epetra_Vector> c_rhs_f_;

    // Interface Coupling PoroFluid - Fluid velocities and pressure are/is coupled
    Teuchos::RCP<ADAPTER::Coupling> icoup_pf_f_;
    // Interface Coupling PoroStructure - Fluid
    Teuchos::RCP<ADAPTER::Coupling> icoup_ps_f_;
    // Interface Coupling PoroStructure - Ale
    Teuchos::RCP<ADAPTER::Coupling> icoup_ps_a_;

    // extractor for fpsi condition from fluid
    Teuchos::RCP<LINALG::MapExtractor> fluidvelpres_extractor_;
    Teuchos::RCP<LINALG::MapExtractor> fluidvel_extractor_;
    // extractor for fpsi condition from poro fluid
    Teuchos::RCP<LINALG::MapExtractor> porofluid_extractor_;
    // extractor for fpsi condition from (poro) structure
    Teuchos::RCP<LINALG::MapExtractor> porostruct_extractor_;
    //! dof row map splitted in inner structure (0), struct interface (1)
    //! inner porofluid (2) and porofluidinterface (3)
    Teuchos::RCP<LINALG::MultiMapExtractor> poro_extractor_;
    Teuchos::RCP<FPSI::UTILS::MapExtractor> fluid_fsifpsi_extractor_;

    // Evaluate is called first time!
    bool isfirstcall_;

    Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap_;
    Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap_;

    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>
        couplingrowtransform_;  /// g_fpsi || F->PF transform (FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>
        couplingrowtransform2_;  /// g_fpsi || PF->F transform (FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>
        couplingrowtransform3_;  /// g_fpsi || PF->F transform (FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>
        couplingrowtransform4_;  /// g_fpsi || F->PS transform (FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> couplingrowtransform5_;
    ;  /// g_fpsi || F->PS transform (FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixColTransform>
        couplingcoltransform_;  /// for Row/Col-Map for Full - FluidField & F->PS transform (FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixColTransform>
        couplingcoltransform2_;  /// for Row/Col-Map for Full - AleField & A->PS transform (FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform>
        couplingrowcoltransform_;  /// g_fpsi/g_fpsi || F->PS/F->PS transform (FPSI/FPSI)
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform>
        couplingrowcoltransform2_;  /// g_fpsi/g_fpsi || F->PF/F->PS transform (FPSI/FPSI)

    // hydraulic conductivity (needed for coupling in case of probtype fps3i)
    double conductivity_;
  };  // FPSICoupling
}  // namespace FPSI

#endif  // FPSI_COUPLING_H
