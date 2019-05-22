/*!----------------------------------------------------------------------
\brief mortar coupling terms of ehl

\level 3
\maintainer Mostafa Faraji

*----------------------------------------------------------------------*/

#ifndef ADAPTER_COUPLING_EHL_MORTAR_H_
#define ADAPTER_COUPLING_EHL_MORTAR_H_

/*---------------------------------------------------------------------*
 | headers                                                             |
 *---------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include "../drt_lib/drt_dserror.H"
#include "adapter_coupling_nonlin_mortar.H"
#include "../drt_contact/contact_interface.H"

/*---------------------------------------------------------------------*
 | forward declarations                                                |
 *---------------------------------------------------------------------*/
namespace LINALG
{
  class SparseMatrix;
  class BlockSparseMatrixBase;
  class MapExtractor;
}  // namespace LINALG

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

namespace ADAPTER
{
  class CouplingEhlMortar : public CouplingNonLinMortar
  {
   public:
    //! @name Construction / Destruction
    //@{
    /*!
    \brief Empty constructor

    */
    CouplingEhlMortar();

    /*!
    \brief Virtual destructor

    */
    virtual ~CouplingEhlMortar(){};
    //@}


    //! @name Evaluation
    //@{
    /*!
    \brief Read Mortar Condition

    */
    virtual void ReadMortarCondition(Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<DRT::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond, Teuchos::ParameterList& input,
        std::map<int, DRT::Node*>& mastergnodes, std::map<int, DRT::Node*>& slavegnodes,
        std::map<int, Teuchos::RCP<DRT::Element>>& masterelements,
        std::map<int, Teuchos::RCP<DRT::Element>>& slaveelements);

    /*!
    \brief initialize routine

    */
    virtual void Setup(Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<DRT::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond);

    /// perform interface integration
    virtual void Integrate(Teuchos::RCP<const Epetra_Vector> disp, const double dt);

    /// perform condensation of contact Lagrange multipliers
    virtual void CondenseContact(Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
        Teuchos::RCP<Epetra_Vector>& combined_RHS, Teuchos::RCP<const Epetra_Vector> disp,
        const double dt);
    virtual void RecoverCoupled(Teuchos::RCP<Epetra_Vector> sinc, Teuchos::RCP<Epetra_Vector> tinc);

    void EvaluateRelMov();
    void StoreDirichletStatus(Teuchos::RCP<const LINALG::MapExtractor> dbcmaps);

    /// check whether this displacement state has already been evaluated
    virtual bool AlreadyEvaluated(Teuchos::RCP<const Epetra_Vector> disp);

    /// Assemble linearization D_{ij,k}*x_i
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleEHLLinD(const Teuchos::RCP<Epetra_Vector> x);

    /// Assemble linearization M_{il,k}*x_i
    virtual Teuchos::RCP<LINALG::SparseMatrix> AssembleEHLLinM(const Teuchos::RCP<Epetra_Vector> x);

    /// Assemble linearization G_{ij,k}*x_i
    Teuchos::RCP<LINALG::SparseMatrix> AssembleSurfGradDeriv(
        const Teuchos::RCP<const Epetra_Vector> x);

    virtual void WriteRestart(IO::DiscretizationWriter& output);
    virtual void ReadRestart(IO::DiscretizationReader& reader);
    void CreateActiveSlipToggle(Teuchos::RCP<Epetra_Vector>* active,
        Teuchos::RCP<Epetra_Vector>* slip, Teuchos::RCP<Epetra_Vector>* active_old = NULL);
    void CreateForceVec(Teuchos::RCP<Epetra_Vector>& n, Teuchos::RCP<Epetra_Vector>& t);
    bool HasContact() { return contact_regularization_; }
    double ContactRes() { return contact_rhs_norm_; }
    double ContactIncr() { return contact_LM_incr_norm_; }
    int ActiveContact() { return interface_->ActiveNodes()->NumGlobalElements(); }
    int SlipContact() { return interface_->SlipNodes()->NumGlobalElements(); }
    //@}

    //! @name Access
    //@{

    /// relative tangential velocity
    virtual Teuchos::RCP<Epetra_Vector> RelTangVel() { return relTangVel_; }
    /// relative tangential velocity derivative
    virtual Teuchos::RCP<LINALG::SparseMatrix> RelTangVelDeriv() { return relTangVel_deriv_; }
    /// average tangential velocity
    virtual Teuchos::RCP<Epetra_Vector> AvTangVel() { return avTangVel_; }
    /// average tangential velocity derivative
    virtual Teuchos::RCP<LINALG::SparseMatrix> AvTangVelDeriv() { return avTangVel_deriv_; }
    /// nodal gap (not weighted)
    virtual Teuchos::RCP<Epetra_Vector> Nodal_Gap() { return nodal_gap_; }
    /// nodal gap (not weighted) derivative
    virtual Teuchos::RCP<LINALG::SparseMatrix> Nodal_GapDeriv() { return deriv_nodal_gap_; }
    /// nodal normals
    virtual Teuchos::RCP<Epetra_Vector> Normals() { return normals_; }
    /// nodal normals derivative
    virtual Teuchos::RCP<LINALG::SparseMatrix> NderivMatrix() { return Nderiv_; }
    /// surfrace gradient operator
    virtual Teuchos::RCP<LINALG::SparseMatrix> SurfGradMatrix() { return SurfGrad_; }
    /// slave+master dof map
    virtual Teuchos::RCP<const Epetra_Map> SMdofMap() { return smdofrowmap_; }
    //@}

   private:
    //! @name Unwanted parent functions
    //@{
    virtual Teuchos::RCP<Epetra_Vector> Gap()
    {
      dserror("stop");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<LINALG::SparseMatrix> NMatrix()
    {
      dserror("stop");
      return Teuchos::null;
    }
    virtual void Evaluate() { dserror("stop"); }
    virtual void Evaluate(Teuchos::RCP<Epetra_Vector> idisp) { dserror("stop"); }
    virtual void EvaluateWithMeshRelocation(
        Teuchos::RCP<DRT::Discretization> slavedis,  ///< slave discretization
        Teuchos::RCP<DRT::Discretization> aledis,    ///< ALE discretization
        Teuchos::RCP<Epetra_Vector>& idisp,          ///< ALE displacements
        const Epetra_Comm& comm,                     ///< communicator
        bool slavewithale                            ///< flag defining if slave is ALE
    )
    {
      dserror("stop");
    }
    virtual void SetupForUQAbuseNormalCalculation(
        Teuchos::RCP<DRT::Discretization> slavedis,  ///< slave discretization
        const Epetra_Comm& comm                      ///< communicator
    )
    {
      dserror("stop");
    }
    virtual void EvaluateGeometry(std::vector<Teuchos::RCP<MORTAR::IntCell>>&
            intcells  //!< vector of mortar integration cells
    )
    {
      dserror("stop");
    }
    //@}

   protected:
    //! @name Assemlby of global vectors/matrices
    //@{

    /// real (not weighted) gap
    void AssembleRealGap();
    /// real (not weighted) gap derivative wrt displacements
    void AssembleRealGapDeriv();
    /// slave sided nodal normals
    void AssembleNormals();
    /// slave sided nodal normals, derivative wrt displacements
    void AssembleNormalsDeriv();
    /// slave sided tangential gradient operator
    void AssembleSurfGrad();
    /// relative and average tangential velocities and their derivatives
    void AssembleInterfaceVelocities(const double dt);
    //@}

    Teuchos::RCP<LINALG::SparseMatrix>
        SurfGrad_;  ///< matrix to compute surface gradient at nodes: Grad(x) \approx TangGrad * x
    Teuchos::RCP<Epetra_Vector> relTangVel_;  ///< relative tangential velocity
    Teuchos::RCP<LINALG::SparseMatrix>
        relTangVel_deriv_;                   ///< derivative of relative tangential velocity
    Teuchos::RCP<Epetra_Vector> avTangVel_;  ///< average tangential velocity
    Teuchos::RCP<LINALG::SparseMatrix>
        avTangVel_deriv_;                    ///< derivative of average tangential velocity
    Teuchos::RCP<Epetra_Vector> nodal_gap_;  ///< NOT the weighted gap but the real distance
    Teuchos::RCP<LINALG::SparseMatrix>
        deriv_nodal_gap_;                        ///< derivative of nodal_gap_ w.r.t. displacements
    Teuchos::RCP<Epetra_Vector> normals_;        ///< nodal normals (slave side)
    Teuchos::RCP<LINALG::SparseMatrix> Nderiv_;  ///< derivtative of nodal slave-sided normals

    Teuchos::RCP<Epetra_Vector>
        evaluated_state_;  ///< displacement state that has last been evaluated

    bool as_converged_;
    double contact_rhs_norm_;
    double contact_LM_incr_norm_;
    Teuchos::RCP<Epetra_Vector>
        fscn_;  // structural contact forces of last time step (needed for time integration)
    Teuchos::RCP<Epetra_Vector> z_;  //  contact lagrange mulitplier

    // recovery of contact LM
    Teuchos::RCP<LINALG::SparseMatrix> dinvA_;  // dinv on active displacement dofs (for recovery)
    Teuchos::RCP<LINALG::SparseMatrix>
        kss_a_;  // Part of structure-stiffness (kss) that corresponds to active slave rows
    Teuchos::RCP<LINALG::SparseMatrix>
        kst_a_;  // Part of coupling-stiffness  (kst) that corresponds to active slave rows
    Teuchos::RCP<Epetra_Vector>
        rs_a_;  // Part of structural residual that corresponds to active slave rows
    Teuchos::RCP<Epetra_Vector> sdirichtoggle_;  // global dirichlet toggle of all slave dofs

    bool contact_regularization_;
    double regularization_thickness_;
    double regularization_compliance_;
  };
}  // namespace ADAPTER

#endif /* ADAPTER_COUPLING_EHL_MORTAR_H_ */
