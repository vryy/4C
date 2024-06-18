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
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
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
    void Setup(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<Core::FE::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond) override;

    /// perform interface integration
    virtual void Integrate(Teuchos::RCP<const Epetra_Vector> disp, const double dt);

    /// perform condensation of contact Lagrange multipliers
    virtual void CondenseContact(Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector>& combined_RHS, Teuchos::RCP<const Epetra_Vector> disp,
        const double dt);
    virtual void RecoverCoupled(Teuchos::RCP<Epetra_Vector> sinc, Teuchos::RCP<Epetra_Vector> tinc);

    void EvaluateRelMov();
    void store_dirichlet_status(Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmaps);

    /// check whether this displacement state has already been evaluated
    virtual bool AlreadyEvaluated(Teuchos::RCP<const Epetra_Vector> disp);

    /// Assemble linearization D_{ij,k}*x_i
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> AssembleEHLLinD(
        const Teuchos::RCP<Epetra_Vector> x);

    /// Assemble linearization M_{il,k}*x_i
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> AssembleEHLLinM(
        const Teuchos::RCP<Epetra_Vector> x);

    /// Assemble linearization G_{ij,k}*x_i
    Teuchos::RCP<Core::LinAlg::SparseMatrix> assemble_surf_grad_deriv(
        const Teuchos::RCP<const Epetra_Vector> x);

    virtual void write_restart(Core::IO::DiscretizationWriter& output);
    virtual void read_restart(Core::IO::DiscretizationReader& reader);
    void create_active_slip_toggle(Teuchos::RCP<Epetra_Vector>* active,
        Teuchos::RCP<Epetra_Vector>* slip, Teuchos::RCP<Epetra_Vector>* active_old = nullptr);
    void CreateForceVec(Teuchos::RCP<Epetra_Vector>& n, Teuchos::RCP<Epetra_Vector>& t);
    bool HasContact() { return contact_regularization_; }
    double ContactRes() { return contact_rhs_norm_; }
    double ContactIncr() { return contact_LM_incr_norm_; }
    int ActiveContact();
    int SlipContact();
    //@}

    //! @name Access
    //@{

    /// relative tangential velocity
    virtual Teuchos::RCP<Epetra_Vector> RelTangVel() { return relTangVel_; }
    /// relative tangential velocity derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> RelTangVelDeriv() { return relTangVel_deriv_; }
    /// average tangential velocity
    virtual Teuchos::RCP<Epetra_Vector> AvTangVel() { return avTangVel_; }
    /// average tangential velocity derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> AvTangVelDeriv() { return avTangVel_deriv_; }
    /// nodal gap (not weighted)
    virtual Teuchos::RCP<Epetra_Vector> Nodal_Gap() { return nodal_gap_; }
    /// nodal gap (not weighted) derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> Nodal_GapDeriv() { return deriv_nodal_gap_; }
    /// nodal normals
    virtual Teuchos::RCP<Epetra_Vector> Normals() { return normals_; }
    /// nodal normals derivative
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> NderivMatrix() { return Nderiv_; }
    /// surfrace gradient operator
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> SurfGradMatrix() { return SurfGrad_; }
    /// slave+master dof map
    virtual Teuchos::RCP<const Epetra_Map> SMdofMap() { return smdofrowmap_; }
    //@}

   private:
    //! @name Unwanted parent functions
    //@{
    Teuchos::RCP<Epetra_Vector> Gap() override
    {
      FOUR_C_THROW("stop");
      return Teuchos::null;
    }
    Teuchos::RCP<Core::LinAlg::SparseMatrix> NMatrix() override
    {
      FOUR_C_THROW("stop");
      return Teuchos::null;
    }
    void evaluate() override { FOUR_C_THROW("stop"); }
    void evaluate(Teuchos::RCP<Epetra_Vector> idisp) override { FOUR_C_THROW("stop"); }
    void evaluate_with_mesh_relocation(
        Teuchos::RCP<Core::FE::Discretization> slavedis,  ///< slave discretization
        Teuchos::RCP<Core::FE::Discretization> aledis,    ///< ALE discretization
        Teuchos::RCP<Epetra_Vector>& idisp,               ///< ALE displacements
        const Epetra_Comm& comm,                          ///< communicator
        bool slavewithale                                 ///< flag defining if slave is ALE
        ) override
    {
      FOUR_C_THROW("stop");
    }

    void EvaluateGeometry(std::vector<Teuchos::RCP<Mortar::IntCell>>&
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
    Teuchos::RCP<Epetra_Vector> relTangVel_;  ///< relative tangential velocity
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        relTangVel_deriv_;                   ///< derivative of relative tangential velocity
    Teuchos::RCP<Epetra_Vector> avTangVel_;  ///< average tangential velocity
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        avTangVel_deriv_;                    ///< derivative of average tangential velocity
    Teuchos::RCP<Epetra_Vector> nodal_gap_;  ///< NOT the weighted gap but the real distance
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        deriv_nodal_gap_;                  ///< derivative of nodal_gap_ w.r.t. displacements
    Teuchos::RCP<Epetra_Vector> normals_;  ///< nodal normals (slave side)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> Nderiv_;  ///< derivtative of nodal slave-sided normals

    Teuchos::RCP<Epetra_Vector>
        evaluated_state_;  ///< displacement state that has last been evaluated

    bool as_converged_;
    double contact_rhs_norm_;
    double contact_LM_incr_norm_;
    Teuchos::RCP<Epetra_Vector>
        fscn_;  // structural contact forces of last time step (needed for time integration)
    Teuchos::RCP<Epetra_Vector> z_;  //  contact lagrange mulitplier

    // recovery of contact LM
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        dinvA_;  // dinv on active displacement dofs (for recovery)
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        kss_a_;  // Part of structure-stiffness (kss) that corresponds to active slave rows
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        kst_a_;  // Part of coupling-stiffness  (kst) that corresponds to active slave rows
    Teuchos::RCP<Epetra_Vector>
        rs_a_;  // Part of structural residual that corresponds to active slave rows
    Teuchos::RCP<Epetra_Vector> sdirichtoggle_;  // global dirichlet toggle of all slave dofs

    bool contact_regularization_;
    double regularization_thickness_;
    double regularization_compliance_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
