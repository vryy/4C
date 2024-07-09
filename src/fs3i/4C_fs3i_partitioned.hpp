/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with general algorithmic routines for
       partitioned solution approaches to fluid-structure-scalar-scalar
       interaction (FS3I), that is, algorithmic routines not specifically
       related to partitioned solution approaches to one -or
       two-way-coupled problem configurations, respectively

\level 2



*----------------------------------------------------------------------*/


#ifndef FOUR_C_FS3I_PARTITIONED_HPP
#define FOUR_C_FS3I_PARTITIONED_HPP


#include "4C_config.hpp"

#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fs3i.hpp"
#include "4C_inpar_fs3i.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Adapter
{
  class Coupling;
  class MortarVolCoupl;
  class ScaTraBaseAlgorithm;
}  // namespace Adapter

namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

namespace Core::LinAlg
{
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class SparseMatrix;
  class Solver;
}  // namespace Core::LinAlg


namespace FS3I
{
  class PartFS3I : public FS3IBase
  {
   public:
    //! constructor of base class for partitioned FS3I
    PartFS3I(const Epetra_Comm& comm);

    //! initialize this class
    void init() override;

    //! setup this class
    void setup() override;

    //! @name overall FS3I system
    //@{

    //! time loop to be defined in inherited classes (structure depends on
    //! considered coupling, i.e. one-way or two-way)
    void timeloop() override = 0;

    //! flag whether time loop should be finished
    bool not_finished() const { return ((step_ < numstep_) and ((time_ + 1e-14) < timemax_)); };

    //! read and set fields needed for restart
    void read_restart() override;

    /// redistribute the  FPSI interface, if running on parallel. Just needed in the case of FPS3I
    void redistribute_interface() override { return; };

    /// create a volmortar object
    Teuchos::RCP<Core::Adapter::MortarVolCoupl> create_vol_mortar_object(
        Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<Core::FE::Discretization> slavedis);

    //! set-up of FSI and ScaTra systems
    void setup_system() override;

    //! test results for individual fields
    void test_results(const Epetra_Comm& comm) override;

    //! information transfer FSI -> ScaTra
    void set_fsi_solution();

    /// set scatra solution on structure field
    void set_struct_scatra_solution() const;

    //! check convergence of monolithic ScaTra problem (depends on which
    // coupling is considered)
    virtual bool scatra_convergence_check(int itnum) = 0;

    //! return communicator
    const Epetra_Comm& get_comm() const { return comm_; }

    /// extract fluid convective and structure convective velocities
    void extract_vel(std::vector<Teuchos::RCP<const Epetra_Vector>>& vel,
        std::vector<Teuchos::RCP<const Epetra_Vector>>& convel) const;

    void set_velocity_fields() const;

    //! routine for preparing time step to be defined in inherited classes
    //! (structure depends on coupling, that is, either one- or two-way)
    virtual void prepare_time_step() = 0;

    void set_mesh_disp() const;

    /// provide wall shear stresses from FS3I subproblem for scatra subproblem
    virtual void set_wall_shear_stresses() const;

    /// extract Wall Shear Stresses at the interface
    void extract_wss(std::vector<Teuchos::RCP<const Epetra_Vector>>& wss) const;

    Teuchos::ParameterList& manipulate_fsi_time_params(const Teuchos::ParameterList& fs3idyn) const;

    //@}

    /// transport quantity from fluid to fluid-scalar
    Teuchos::RCP<const Epetra_Vector> fluid_to_fluid_scalar(
        const Teuchos::RCP<const Epetra_Vector> fluidvector) const;

    /// transport quantity from fluid-scalar to fluid
    Teuchos::RCP<const Epetra_Vector> fluid_scalar_to_fluid(
        const Teuchos::RCP<const Epetra_Vector> fluidscalarvector) const;

    /// transport quantity from structure to structure-scalar
    Teuchos::RCP<const Epetra_Vector> structure_to_structure_scalar(
        const Teuchos::RCP<const Epetra_Vector> structurevector) const;

    /// transport quantity from structure-scalar to structure
    Teuchos::RCP<const Epetra_Vector> structure_scalar_to_structure(
        const Teuchos::RCP<const Epetra_Vector> structurescalavector) const;

   private:
    /// transport quantity from i-th volmortar master to i-th volmortar slave
    Teuchos::RCP<const Epetra_Vector> vol_mortar_master_to_slavei(
        const int i, const Teuchos::RCP<const Epetra_Vector> mastervector) const;

    /// transport quantity from i-th volmortar slave to i-th volmortar master
    Teuchos::RCP<const Epetra_Vector> vol_mortar_slave_to_masteri(
        const int i, const Teuchos::RCP<const Epetra_Vector> slavevector) const;

   protected:
    /// fsi algorithm
    Teuchos::RCP<FSI::Monolithic> fsi_;

    /// vector of scatra volume couplings (i.e. fluid to fluid-scalar and structure to
    /// structure-scalar)
    std::vector<Inpar::FS3I::VolumeCoupling> volume_fieldcouplings_;

   private:
    //! volume coupling (using mortar) adapter
    std::vector<Teuchos::RCP<Core::Adapter::MortarVolCoupl>> volume_coupling_objects_;

    /// communication (mainly for screen output)
    const Epetra_Comm& comm_;

    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> fluidscatra_;

    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> structscatra_;
  };
}  // namespace FS3I

FOUR_C_NAMESPACE_CLOSE

#endif
