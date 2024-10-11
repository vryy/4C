/*----------------------------------------------------------------------*/
/*! \file
\brief mortar coupling terms of ehl

\level 3

*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_COUPLING_EHL_MORTAR_HPP
#define FOUR_C_ADAPTER_COUPLING_EHL_MORTAR_HPP

/*---------------------------------------------------------------------*
 | headers                                                             |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_adapter_coupling_nonlin_mortar.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                                |
 *---------------------------------------------------------------------*/
namespace Core::LinAlg
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
}  // namespace Core::LinAlg

namespace Core::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace Core::IO

namespace Adapter
{
  class CouplingEhlMortar : public CouplingNonLinMortar
  {
   public:
    //! @name Construction / Destruction
    //@{
    /*!
    \brief Constructor.
    */
    CouplingEhlMortar(int spatial_dimension, Teuchos::ParameterList mortar_coupling_params,
        Teuchos::ParameterList contact_dynamic_params,
        Core::FE::ShapeFunctionType shape_function_type);

    //@}


    //! @name Evaluation
    //@{
    /*!
    \brief Read Mortar Condition

    */
    void read_mortar_condition(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond, Teuchos::ParameterList& input,
        std::map<int, Core::Nodes::Node*>& mastergnodes,
        std::map<int, Core::Nodes::Node*>& slavegnodes,
        std::map<int, Teuchos::RCP<Core::Elements::Element>>& masterelements,
        std::map<int, Teuchos::RCP<Core::Elements::Element>>& slaveelements) override;

    /*!
    \brief initialize routine

    */
    void setup(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond) override;

    /// perform interface integration
    virtual void integrate(Teuchos::RCP<const Core::LinAlg::Vector<double>> disp, const double dt);

    /// perform condensation of contact Lagrange multipliers
    virtual void condense_contact(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Core::LinAlg::Vector<double>>& combined_RHS,
        Teuchos::RCP<const Core::LinAlg::Vector<double>> disp, const double dt);
    virtual void recover_coupled(Teuchos::RCP<Core::LinAlg::Vector<double>> sinc,
        Teuchos::RCP<Core::LinAlg::Vector<double>> tinc);

    void evaluate_rel_mov();
    void store_dirichlet_status(const Core::LinAlg::MapExtractor& dbcmaps);

    /// check whether this displacement state has already been evaluated
    virtual bool already_evaluated(Teuchos::RCP<const Core::LinAlg::Vector<double>> disp);

    /// Assemble linearization D_{ij,k}*x_i
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> assemble_ehl_lin_d(
        const Teuchos::RCP<Core::LinAlg::Vector<double>> x);

    /// Assemble linearization M_{il,k}*x_i
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> assemble_ehl_lin_m(
        const Teuchos::RCP<Core::LinAlg::Vector<double>> x);

    /// Assemble linearization G_{ij,k}*x_i
    Teuchos::RCP<Core::LinAlg::SparseMatrix> assemble_surf_grad_deriv(
        const Core::LinAlg::Vector<double>& x);

    virtual void write_restart(Core::IO::DiscretizationWriter& output);
    virtual void read_restart(Core::IO::DiscretizationReader& reader);
    void create_active_slip_toggle(Teuchos::RCP<Core::LinAlg::Vector<double>>* active,
        Teuchos::RCP<Core::LinAlg::Vector<double>>* slip,
        Teuchos::RCP<Core::LinAlg::Vector<double>>* active_old = nullptr);
    void create_force_vec(Teuchos::RCP<Core::LinAlg::Vector<double>>& n,
        Teuchos::RCP<Core::LinAlg::Vector<double>>& t);
    bool has_contact() { return contact_regularization_; }
    double contact_res() { return contact_rhs_norm_; }
    double contact_incr() { return contact_LM_incr_norm_; }
    int active_contact();
    int slip_contact();
    //@}

    //! @name Access
    //@{

    /// relative tangential velocity
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> rel_tang_vel() { return relTangVel_; }
    /// relative tangential velocity derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> rel_tang_vel_deriv()
    {
      return relTangVel_deriv_;
    }
    /// average tangential velocity
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> av_tang_vel() { return avTangVel_; }
    /// average tangential velocity derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> av_tang_vel_deriv()
    {
      return avTangVel_deriv_;
    }
    /// nodal gap (not weighted)
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> nodal_gap() { return nodal_gap_; }
    /// nodal gap (not weighted) derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> nodal_gap_deriv() { return deriv_nodal_gap_; }
    /// nodal normals
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> normals() { return normals_; }
    /// nodal normals derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> nderiv_matrix() { return Nderiv_; }
    /// surfrace gradient operator
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> surf_grad_matrix() { return SurfGrad_; }
    /// slave+master dof map
    virtual Teuchos::RCP<const Epetra_Map> s_mdof_map() { return smdofrowmap_; }
    //@}

   private:
    //! @name Unwanted parent functions
    //@{
    Teuchos::RCP<Core::LinAlg::Vector<double>> gap() override
    {
      FOUR_C_THROW("stop");
      return Teuchos::null;
    }
    Teuchos::RCP<Core::LinAlg::SparseMatrix> n_matrix() override
    {
      FOUR_C_THROW("stop");
      return Teuchos::null;
    }
    void evaluate() override { FOUR_C_THROW("stop"); }
    void evaluate(Teuchos::RCP<Core::LinAlg::Vector<double>> idisp) override
    {
      FOUR_C_THROW("stop");
    }
    void evaluate_with_mesh_relocation(
        Teuchos::RCP<Core::FE::Discretization> slavedis,    ///< slave discretization
        Teuchos::RCP<Core::FE::Discretization> aledis,      ///< ALE discretization
        Teuchos::RCP<Core::LinAlg::Vector<double>>& idisp,  ///< ALE displacements
        const Epetra_Comm& comm,                            ///< communicator
        bool slavewithale                                   ///< flag defining if slave is ALE
        ) override
    {
      FOUR_C_THROW("stop");
    }

    void evaluate_geometry(std::vector<Teuchos::RCP<Mortar::IntCell>>&
            intcells  //!< vector of mortar integration cells
        ) override
    {
      FOUR_C_THROW("stop");
    }
    //@}

   protected:
    //! @name Assemlby of global vectors/matrices
    //@{

    /// real (not weighted) gap
    void assemble_real_gap();
    /// real (not weighted) gap derivative wrt displacements
    void assemble_real_gap_deriv();
    /// slave sided nodal normals
    void assemble_normals();
    /// slave sided nodal normals, derivative wrt displacements
    void assemble_normals_deriv();
    /// slave sided tangential gradient operator
    void assemble_surf_grad();
    /// relative and average tangential velocities and their derivatives
    void assemble_interface_velocities(const double dt);
    //@}

    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        SurfGrad_;  ///< matrix to compute surface gradient at nodes: Grad(x) \approx TangGrad * x
    Teuchos::RCP<Core::LinAlg::Vector<double>> relTangVel_;  ///< relative tangential velocity
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        relTangVel_deriv_;  ///< derivative of relative tangential velocity
    Teuchos::RCP<Core::LinAlg::Vector<double>> avTangVel_;  ///< average tangential velocity
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        avTangVel_deriv_;  ///< derivative of average tangential velocity
    Teuchos::RCP<Core::LinAlg::Vector<double>>
        nodal_gap_;  ///< NOT the weighted gap but the real distance
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        deriv_nodal_gap_;  ///< derivative of nodal_gap_ w.r.t. displacements
    Teuchos::RCP<Core::LinAlg::Vector<double>> normals_;  ///< nodal normals (slave side)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> Nderiv_;  ///< derivtative of nodal slave-sided normals

    Teuchos::RCP<Core::LinAlg::Vector<double>>
        evaluated_state_;  ///< displacement state that has last been evaluated

    bool as_converged_;
    double contact_rhs_norm_;
    double contact_LM_incr_norm_;
    Teuchos::RCP<Core::LinAlg::Vector<double>>
        fscn_;  // structural contact forces of last time step (needed for time integration)
    Teuchos::RCP<Core::LinAlg::Vector<double>> z_;  //  contact lagrange mulitplier

    // recovery of contact LM
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        dinvA_;  // dinv on active displacement dofs (for recovery)
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        kss_a_;  // Part of structure-stiffness (kss) that corresponds to active slave rows
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        kst_a_;  // Part of coupling-stiffness  (kst) that corresponds to active slave rows
    Teuchos::RCP<Core::LinAlg::Vector<double>>
        rs_a_;  // Part of structural residual that corresponds to active slave rows
    Teuchos::RCP<Core::LinAlg::Vector<double>>
        sdirichtoggle_;  // global dirichlet toggle of all slave dofs

    bool contact_regularization_;
    double regularization_thickness_;
    double regularization_compliance_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
