/*----------------------------------------------------------------------*/
/*! \file

 \brief  Basis of all porous media algorithms

 \level 2

 *-----------------------------------------------------------------------*/
#ifndef FOUR_C_POROELAST_BASE_HPP
#define FOUR_C_POROELAST_BASE_HPP


#include "baci_config.hpp"

#include "baci_adapter_algorithmbase.hpp"
#include "baci_adapter_field.hpp"
#include "baci_coupling_adapter.hpp"
#include "baci_coupling_adapter_volmortar.hpp"
#include "baci_linalg_mapextractor.hpp"
#include "baci_poroelast_utils.hpp"
#include "baci_utils_exceptions.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

namespace CORE::LINALG
{
  class MapExtractor;
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace CORE::LINALG

namespace ADAPTER
{
  class FluidPoro;
  class Coupling;
  class FPSIStructureWrapper;
  class MortarVolCoupl;
}  // namespace ADAPTER

namespace POROELAST
{
  class NoPenetrationConditionHandle;

  //! base class for porous media algorithm
  class PoroBase : public ADAPTER::AlgorithmBase, public ADAPTER::Field
  {
   public:
    //! create using a Epetra_Comm
    explicit PoroBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter);

    //! read restart data
    void ReadRestart(const int step) override;

    //! outer level time loop
    virtual void TimeLoop();

    //! initialize system
    virtual void SetupSystem() = 0;

    //! perform result tests
    virtual void TestResults(const Epetra_Comm& comm);

    //! build combined dirichlet map for the monolithic problem
    virtual void BuildCombinedDBCMap()
    {
      dserror(
          "BuildCombinedDBCMap() not implemented in base class. must be implemented in sub "
          "classes.");
    }

    //! @name access methods

    //! access to structural field
    const Teuchos::RCP<ADAPTER::FPSIStructureWrapper>& StructureField() { return structure_; }

    //! access to fluid field
    const Teuchos::RCP<ADAPTER::FluidPoro>& FluidField() { return fluid_; }

    //! composed system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override
    {
      dserror("SystemMatrix() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! block system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      dserror("BlockSystemMatrix() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap() override
    {
      dserror("DofRowMap() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! dof row map of Structure field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapStructure() = 0;

    //! dof row map of Fluid field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapFluid() = 0;

    //! extractor to communicate between full monolithic map and block maps
    virtual Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> Extractor() const
    {
      dserror("ExtractorPointer only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> CombinedDBCMap() const
    {
      dserror("CombinedDBCMap() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! return rhs of poro problem
    Teuchos::RCP<const Epetra_Vector> RHS() override
    {
      dserror("RHS() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //!@}

    //! update all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void UpdateStateIncrementally(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
        ) override
    {
      dserror("UpdateStateIncrementally() only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
        ) override
    {
      dserror("Evaluate() only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx)
    {
      dserror("Evaluate(sx,fx) only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc,  //!< increment between iteration i and i+1
        bool firstiter) override
    {
      dserror("Evaluate() only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx, bool firstiter)
    {
      dserror("Evaluate(sx,fx) only available for monolithic schemes!");
    }

    //! solve time step (depending on algorithm)
    virtual void Solve() = 0;

    //! perform one time step (setup + solve + output)
    virtual void DoTimeStep() = 0;

    //! setup solver (for monolithic only)
    virtual bool SetupSolver() { return false; }

    virtual void SetupRHS(bool firstcall = false)
    {
      dserror("SetupRHS() only available for monolithic schemes!");
    }

    //! @name Time loop building blocks

    //! start a new time step
    void PrepareTimeStep() override;

    //! take current results for converged and save for next time step
    void Update() override;

    //! calculate stresses, strains, energies
    void PrepareOutput(bool force_prepare_timestep) override;

    //! output
    void Output(bool forced_writerestart = false) override;

    //!@}

    //! return whether the poro discretization contains submeshes (i.e. it is coupled with a pure
    //! solid)
    bool HasSubmeshes() { return submeshes_; }

    //! return coupling object
    CORE::ADAPTER::Coupling& FluidStructureCoupling() { return *coupling_fluid_structure_; }

   protected:
    //! @name Transfer helpers

    //! field transform
    Teuchos::RCP<Epetra_Vector> StructureToFluidField(Teuchos::RCP<const Epetra_Vector> iv);

    //!@}

    //! @name Transfer methods

    //! set fluid solution on structure
    void SetFluidSolution();

    //! set structure solution on fluid
    void SetStructSolution();

    //!@}

    //! Extractor used for constraint structure
    Teuchos::RCP<CORE::LINALG::MapExtractor> cond_splitter_;

    //! true if the poroelast problem is only part of a larger problem (like e.g. in FPSI)
    bool is_part_of_multifield_problem_;

    //! Extractor used for poro structure interaction
    Teuchos::RCP<CORE::LINALG::MapExtractor> psi_extractor_;

    //! helper class for no penetration condition
    Teuchos::RCP<NoPenetrationConditionHandle> nopen_handle_;

    //! flag for partial integration condition of porous fluid continuity equation
    bool part_int_cond_;

    //! flag for pressure coupling condition
    bool pres_int_cond_;

    //! flag for additional porosity degree of freedom
    bool porosity_dof_;

    Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter_;



    //! @name Volume Mortar stuff

    //! flag for matchinggrid
    const bool matchinggrid_;

    //! volume coupling (using mortar) adapter
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> volcoupl_;
    //!@}

    //! flag for old time integration
    const bool oldstructimint_;  // delete this once unused!

    //! flag indicating if nitsche contact is active
    bool nit_contact_;

   private:
    //! setup of everything, that is needed for the volumetric coupling
    void SetupCoupling();

    //! add dof set of structure/fluid discretization to fluid/structure discretization
    void ReplaceDofSets();

    //! check for special poro conditions and set bools
    void CheckForPoroConditions();

    //! flag indicating not fully overlapping fluid and structure discretization
    bool submeshes_;

    //! coupling of fluid and structure (whole field)
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupling_fluid_structure_;

    //! @name Underlying fields

    //! underlying structure of the poro problem
    Teuchos::RCP<ADAPTER::FPSIStructureWrapper> structure_;

    //! underlying fluid of the poro problem
    Teuchos::RCP<ADAPTER::FluidPoro> fluid_;

    //!@}
  };

  //! helper class for no penetration condition
  class NoPenetrationConditionHandle
  {
   public:
    //! create using a Epetra_Comm
    explicit NoPenetrationConditionHandle(std::vector<DRT::Condition*> nopencond)
        : cond_ids_(Teuchos::null),
          cond_dofs_(Teuchos::null),
          cond_rhs_(Teuchos::null),
          nopencond_(std::move(nopencond)),
          nopenetration_(Teuchos::null),
          has_cond_(false),
          fluid_fluid_constraint_matrix_(Teuchos::null),
          fluid_structure_constraint_matrix_(Teuchos::null),
          structure_vel_constraint_matrix_(Teuchos::null)
    {
      if (nopencond_.size())
      {
        has_cond_ = true;
        cond_ids_ = Teuchos::rcp(new std::set<int>());
      }
    }

    //! build map containing dofs with no penetration condition (fluid)
    void BuidNoPenetrationMap(const Epetra_Comm& comm, Teuchos::RCP<const Epetra_Map> dofRowMap);

    //! apply rhs terms of no penetration condition to global rhs vector
    void ApplyCondRHS(Teuchos::RCP<Epetra_Vector> iterinc, Teuchos::RCP<Epetra_Vector> rhs);

    //! return no penetration map extractor
    Teuchos::RCP<const CORE::LINALG::MapExtractor> Extractor() { return nopenetration_; }

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int>> CondIDs() { return cond_ids_; }

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> CondVector() { return cond_dofs_; }

    //! check if a no penetration condition exists
    bool HasCond() { return has_cond_; }

    //! return condrhs
    Teuchos::RCP<Epetra_Vector> RHS() { return cond_rhs_; }

    //! clear everything that is needed for coupling
    void Clear(POROELAST::coupltype coupltype = POROELAST::undefined);

    //! setup coupling matrixes and vecors
    void Setup(Teuchos::RCP<const Epetra_Map> dofRowMap, const Epetra_Map* dofRowMapFluid);

    //! return constraint matrix, that fits to coupling type
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ConstraintMatrix(POROELAST::coupltype coupltype);

    //! return constraint matrix for structure velocity coupling
    Teuchos::RCP<CORE::LINALG::SparseMatrix> StructVelConstraintMatrix(
        POROELAST::coupltype coupltype);

   private:
    //! set containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int>> cond_ids_;

    //! vector marking dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> cond_dofs_;

    //! vector containing rhs terms from no penetration condition
    Teuchos::RCP<Epetra_Vector> cond_rhs_;

    //! vector containing no penetration - conditions
    std::vector<DRT::Condition*> nopencond_;

    //! Extractor used for no penetration condition
    Teuchos::RCP<CORE::LINALG::MapExtractor> nopenetration_;

    //! flag indicating if a no penetration condition exists
    bool has_cond_;

    //! @name coupling matrices
    Teuchos::RCP<CORE::LINALG::SparseMatrix> fluid_fluid_constraint_matrix_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> fluid_structure_constraint_matrix_;
    Teuchos::RCP<CORE::LINALG::SparseMatrix> structure_vel_constraint_matrix_;
    //!@}
  };
}  // namespace POROELAST

FOUR_C_NAMESPACE_CLOSE

#endif
