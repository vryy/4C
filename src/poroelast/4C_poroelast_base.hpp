/*----------------------------------------------------------------------*/
/*! \file

 \brief  Basis of all porous media algorithms

 \level 2

 *-----------------------------------------------------------------------*/
#ifndef FOUR_C_POROELAST_BASE_HPP
#define FOUR_C_POROELAST_BASE_HPP


#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_adapter_field.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fem_condition.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::LinAlg
{
  class MapExtractor;
  class SparseMatrix;
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace Adapter
{
  class FluidPoro;
  class Coupling;
  class FPSIStructureWrapper;
  class MortarVolCoupl;
}  // namespace Adapter

namespace PoroElast
{
  class NoPenetrationConditionHandle;

  //! base class for porous media algorithm
  class PoroBase : public Adapter::AlgorithmBase, public Adapter::Field
  {
   public:
    //! create using a Epetra_Comm
    explicit PoroBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter);

    //! read restart data
    void read_restart(const int step) override;

    //! outer level time loop
    virtual void TimeLoop();

    //! initialize system
    virtual void SetupSystem() = 0;

    //! perform result tests
    virtual void TestResults(const Epetra_Comm& comm);

    //! build combined dirichlet map for the monolithic problem
    virtual void build_combined_dbc_map()
    {
      FOUR_C_THROW(
          "build_combined_dbc_map() not implemented in base class. must be implemented in sub "
          "classes.");
    }

    //! @name access methods

    //! access to structural field
    const Teuchos::RCP<Adapter::FPSIStructureWrapper>& structure_field() { return structure_; }

    //! access to fluid field
    const Teuchos::RCP<Adapter::FluidPoro>& fluid_field() { return fluid_; }

    //! composed system matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> system_matrix() override
    {
      FOUR_C_THROW("system_matrix() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! block system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override
    {
      FOUR_C_THROW("block_system_matrix() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! full monolithic dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() override
    {
      FOUR_C_THROW("dof_row_map() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! dof row map of Structure field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapStructure() = 0;

    //! dof row map of Fluid field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapFluid() = 0;

    //! extractor to communicate between full monolithic map and block maps
    virtual Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> Extractor() const
    {
      FOUR_C_THROW("ExtractorPointer only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> combined_dbc_map() const
    {
      FOUR_C_THROW("combined_dbc_map() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //! return rhs of poro problem
    Teuchos::RCP<const Epetra_Vector> RHS() override
    {
      FOUR_C_THROW("RHS() only available for monolithic schemes!");
      return Teuchos::null;
    }

    //!@}

    //! update all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void update_state_incrementally(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
        ) override
    {
      FOUR_C_THROW("update_state_incrementally() only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc  //!< increment between iteration i and i+1
        ) override
    {
      FOUR_C_THROW("Evaluate() only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx)
    {
      FOUR_C_THROW("Evaluate(sx,fx) only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc,  //!< increment between iteration i and i+1
        bool firstiter) override
    {
      FOUR_C_THROW("Evaluate() only available for monolithic schemes!");
    }

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx, bool firstiter)
    {
      FOUR_C_THROW("Evaluate(sx,fx) only available for monolithic schemes!");
    }

    //! solve time step (depending on algorithm)
    virtual void Solve() = 0;

    //! perform one time step (setup + solve + output)
    virtual void do_time_step() = 0;

    //! setup solver (for monolithic only)
    virtual bool SetupSolver() { return false; }

    virtual void setup_rhs(bool firstcall = false)
    {
      FOUR_C_THROW("setup_rhs() only available for monolithic schemes!");
    }

    //! @name Time loop building blocks

    //! start a new time step
    void prepare_time_step() override;

    //! take current results for converged and save for next time step
    void Update() override;

    //! calculate stresses, strains, energies
    void prepare_output(bool force_prepare_timestep) override;

    //! output
    void Output(bool forced_writerestart = false) override;

    //!@}

    //! return whether the poro discretization contains submeshes (i.e. it is coupled with a pure
    //! solid)
    bool HasSubmeshes() { return submeshes_; }

    //! return coupling object
    Core::Adapter::Coupling& fluid_structure_coupling() { return *coupling_fluid_structure_; }

   protected:
    //! @name Transfer helpers

    //! field transform
    Teuchos::RCP<Epetra_Vector> structure_to_fluid_field(Teuchos::RCP<const Epetra_Vector> iv);

    //!@}

    //! @name Transfer methods

    //! set fluid solution on structure
    void set_fluid_solution();

    //! set structure solution on fluid
    void set_struct_solution();

    //!@}

    //! Extractor used for constraint structure
    Teuchos::RCP<Core::LinAlg::MapExtractor> cond_splitter_;

    //! true if the poroelast problem is only part of a larger problem (like e.g. in FPSI)
    bool is_part_of_multifield_problem_;

    //! Extractor used for poro structure interaction
    Teuchos::RCP<Core::LinAlg::MapExtractor> psi_extractor_;

    //! helper class for no penetration condition
    Teuchos::RCP<NoPenetrationConditionHandle> nopen_handle_;

    //! flag for partial integration condition of porous fluid continuity equation
    bool part_int_cond_;

    //! flag for pressure coupling condition
    bool pres_int_cond_;

    //! flag for additional porosity degree of freedom
    bool porosity_dof_;

    Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter_;



    //! @name Volume Mortar stuff

    //! flag for matchinggrid
    const bool matchinggrid_;

    //! volume coupling (using mortar) adapter
    Teuchos::RCP<Core::Adapter::MortarVolCoupl> volcoupl_;
    //!@}

    //! flag for old time integration
    const bool oldstructimint_;  // delete this once unused!

    //! flag indicating if nitsche contact is active
    bool nit_contact_;

   private:
    //! setup of everything, that is needed for the volumetric coupling
    void setup_coupling();

    //! add dof set of structure/fluid discretization to fluid/structure discretization
    void replace_dof_sets();

    //! check for special poro conditions and set bools
    void check_for_poro_conditions();

    //! flag indicating not fully overlapping fluid and structure discretization
    bool submeshes_;

    //! coupling of fluid and structure (whole field)
    Teuchos::RCP<Core::Adapter::Coupling> coupling_fluid_structure_;

    //! @name Underlying fields

    //! underlying structure of the poro problem
    Teuchos::RCP<Adapter::FPSIStructureWrapper> structure_;

    //! underlying fluid of the poro problem
    Teuchos::RCP<Adapter::FluidPoro> fluid_;

    //!@}
  };

  //! helper class for no penetration condition
  class NoPenetrationConditionHandle
  {
   public:
    //! create using a Epetra_Comm
    explicit NoPenetrationConditionHandle(std::vector<Core::Conditions::Condition*> nopencond)
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
    void buid_no_penetration_map(const Epetra_Comm& comm, Teuchos::RCP<const Epetra_Map> dofRowMap);

    //! apply rhs terms of no penetration condition to global rhs vector
    void ApplyCondRHS(Teuchos::RCP<Epetra_Vector> iterinc, Teuchos::RCP<Epetra_Vector> rhs);

    //! return no penetration map extractor
    Teuchos::RCP<const Core::LinAlg::MapExtractor> Extractor() { return nopenetration_; }

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int>> CondIDs() { return cond_ids_; }

    //! return vector containing global IDs of dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> CondVector() { return cond_dofs_; }

    //! check if a no penetration condition exists
    bool HasCond() { return has_cond_; }

    //! return condrhs
    Teuchos::RCP<Epetra_Vector> RHS() { return cond_rhs_; }

    //! clear everything that is needed for coupling
    void Clear(PoroElast::Coupltype coupltype = PoroElast::undefined);

    //! setup coupling matrixes and vecors
    void Setup(Teuchos::RCP<const Epetra_Map> dofRowMap, const Epetra_Map* dofRowMapFluid);

    //! return constraint matrix, that fits to coupling type
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ConstraintMatrix(PoroElast::Coupltype coupltype);

    //! return constraint matrix for structure velocity coupling
    Teuchos::RCP<Core::LinAlg::SparseMatrix> struct_vel_constraint_matrix(
        PoroElast::Coupltype coupltype);

   private:
    //! set containing global IDs of dofs with no penetration condition
    Teuchos::RCP<std::set<int>> cond_ids_;

    //! vector marking dofs with no penetration condition
    Teuchos::RCP<Epetra_Vector> cond_dofs_;

    //! vector containing rhs terms from no penetration condition
    Teuchos::RCP<Epetra_Vector> cond_rhs_;

    //! vector containing no penetration - conditions
    std::vector<Core::Conditions::Condition*> nopencond_;

    //! Extractor used for no penetration condition
    Teuchos::RCP<Core::LinAlg::MapExtractor> nopenetration_;

    //! flag indicating if a no penetration condition exists
    bool has_cond_;

    //! @name coupling matrices
    Teuchos::RCP<Core::LinAlg::SparseMatrix> fluid_fluid_constraint_matrix_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> fluid_structure_constraint_matrix_;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> structure_vel_constraint_matrix_;
    //!@}
  };
}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
