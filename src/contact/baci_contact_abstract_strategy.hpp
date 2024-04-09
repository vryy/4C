/*---------------------------------------------------------------------*/
/*! \file
\brief Main abstract class for contact solution strategies

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_ABSTRACT_STRATEGY_HPP
#define FOUR_C_CONTACT_ABSTRACT_STRATEGY_HPP

#include "baci_config.hpp"

#include "baci_contact_paramsinterface.hpp"
#include "baci_contact_utils.hpp"
#include "baci_inpar_contact.hpp"
#include "baci_inpar_mortar.hpp"
#include "baci_mortar_strategy_base.hpp"

#include <Epetra_Operator.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN

// forward declarations
namespace NOX::NLN
{
  class Group;
}  // namespace NOX::NLN

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace CORE::LINALG
{
  class MultiMapExtractor;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace CONTACT
{
  // forward declarations
  class Interface;
  class NoxInterface;

  /*! \brief Data container object for the abstract strategy
   *
   *  This object makes it possible to interchange and share the current state of the
   *  contact simulation between different strategy objects. By using this the
   *  actual strategy stays stateless!
   *
   *  \author  hiermeier
   *  \date 05/16 */
  class AbstractStratDataContainer : public MORTAR::StratDataContainer
  {
   public:
    //! constructor
    AbstractStratDataContainer();

    //! @name Accessors
    //!@{

    //! Return parallel unbalance factors (evaluation time) for current time step \f$t_{n+1}\f$
    std::vector<double>& UnbalanceTimeFactors() { return unbalanceEvaluationTime_; };
    const std::vector<double>& UnbalanceTimeFactors() const { return unbalanceEvaluationTime_; };

    //! Return parallel unbalance factors (number of slave elements) for current time step
    //! \f$t_{n+1}\f$
    std::vector<int>& UnbalanceElementFactors() { return unbalanceNumSlaveElements_; };
    const std::vector<int>& UnbalanceElementFactors() const { return unbalanceNumSlaveElements_; };

    //! return global Lagrange mult. dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GLmDofRowMapPtr() { return glmdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> GLmDofRowMapPtr() const { return glmdofrowmap_; };

    //! return global reference dof row map for self contact Lagr. multipliers (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GSelfContactRefDofRowMapPtr() { return gscrefdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> GSelfContactRefDofRowMapPtr() const { return gscrefdofrowmap_; };

    //! return global self-contact Lagrange mult. dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GSelfContactLmDofRowMapPtr() { return gsclmdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> GSelfContactLmDofRowMapPtr() const { return gsclmdofrowmap_; };

    //! return global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GSlNodeRowMapPtr() { return gsnoderowmap_; };
    Teuchos::RCP<const Epetra_Map> GSlNodeRowMapPtr() const { return gsnoderowmap_; };

    //! return global master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GMaNodeRowMapPtr() { return gmnoderowmap_; };
    Teuchos::RCP<const Epetra_Map> GMaNodeRowMapPtr() const { return gmnoderowmap_; };

    //! return global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GSlDofRowMapPtr() { return gsdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> GSlDofRowMapPtr() const { return gsdofrowmap_; };

    //! return global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GMaDofRowMapPtr() { return gmdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> GMaDofRowMapPtr() const { return gmdofrowmap_; };

    //! return global internal dof row map
    Teuchos::RCP<Epetra_Map>& GInternalDofRowMapPtr() { return gndofrowmap_; };
    Teuchos::RCP<const Epetra_Map> GInternalDofRowMapPtr() const { return gndofrowmap_; };

    //! return global slave and master dof row map (s+m map)
    Teuchos::RCP<Epetra_Map>& GSlMaDofRowMapPtr() { return gsmdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> GSlMaDofRowMapPtr() const { return gsmdofrowmap_; };

    //! return global displacement dof row map (s+m+n map)
    Teuchos::RCP<Epetra_Map>& GDispDofRowMapPtr() { return gdisprowmap_; };
    Teuchos::RCP<const Epetra_Map> GDispDofRowMapPtr() const { return gdisprowmap_; };

    //! return global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GActiveNodeRowMapPtr() { return gactivenodes_; };
    Teuchos::RCP<const Epetra_Map> GActiveNodeRowMapPtr() const { return gactivenodes_; };
    Epetra_Map& GActiveNodeRowMap()
    {
      if (gactivenodes_.is_null()) dserror("The gactivenodes_ is not initialized!");
      return *gactivenodes_;
    }

    //! return global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GActiveDofRowMapPtr() { return gactivedofs_; };
    Teuchos::RCP<const Epetra_Map> GActiveDofRowMapPtr() const { return gactivedofs_; };
    Epetra_Map& GActiveDofRowMap()
    {
      if (gactivedofs_.is_null()) dserror("The gAugActiveSlaveDofsPtr_ is not initialized!");
      return *gactivedofs_;
    }


    //! return global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GInActiveNodeRowMapPtr() { return ginactivenodes_; };
    Teuchos::RCP<const Epetra_Map> GInActiveNodeRowMapPtr() const { return ginactivenodes_; };
    Epetra_Map& GInActiveNodeRowMap()
    {
      if (ginactivenodes_.is_null()) dserror("The ginactivenodes_ is not initialized!");
      return *ginactivenodes_;
    }

    //! return global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GInActiveDofRowMapPtr() { return ginactivedofs_; };
    Teuchos::RCP<const Epetra_Map> GInActiveDofRowMapPtr() const { return ginactivedofs_; };
    Epetra_Map& GInActiveDofRowMap()
    {
      if (ginactivedofs_.is_null()) dserror("The gAugActiveSlaveDofsPtr_ is not initialized!");
      return *ginactivedofs_;
    }


    //! return global active slave dof row map in normal direction (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GActiveNDofRowMapPtr() { return gactiven_; };
    Teuchos::RCP<const Epetra_Map> GActiveNDofRowMapPtr() const { return gactiven_; };
    Epetra_Map& GActiveNDofRowMap()
    {
      if (gactiven_.is_null()) dserror("The gactiven_ is not initialized!");
      return *gactiven_;
    }

    //! return global active slave dof row map in tangential direction (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GActiveTDofRowMapPtr() { return gactivet_; };
    Teuchos::RCP<const Epetra_Map> GActiveTDofRowMapPtr() const { return gactivet_; };
    Epetra_Map& GActiveTDofRowMap()
    {
      if (gactivet_.is_null()) dserror("The gactivet_ is not initialized!");
      return *gactivet_;
    }

    //! return global slip slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GSlipNodeRowMapPtr() { return gslipnodes_; };
    Teuchos::RCP<const Epetra_Map> GSlipNodeRowMapPtr() const { return gslipnodes_; };

    //! return global slip slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GSlipDofRowMapPtr() { return gslipdofs_; };
    Teuchos::RCP<const Epetra_Map> GSlipDofRowMapPtr() const { return gslipdofs_; };

    //! return global slip slave dof row map in tangential direction (of all interfaces)
    Teuchos::RCP<Epetra_Map>& GSlipTDofRowMapPtr() { return gslipt_; };
    Teuchos::RCP<const Epetra_Map> GSlipTDofRowMapPtr() const { return gslipt_; };

    //! return global slave dof row map associated with vertex nodes
    Teuchos::RCP<Epetra_Map>& GSDofVertexRowMapPtr() { return gsdofVertex_; };
    Teuchos::RCP<const Epetra_Map> GSDofVertexRowMapPtr() const { return gsdofVertex_; };

    //! return global slave dof row map associated with edge nodes
    Teuchos::RCP<Epetra_Map>& GSDofEdgeRowMapPtr() { return gsdofEdge_; };
    Teuchos::RCP<const Epetra_Map> GSDofEdgeRowMapPtr() const { return gsdofEdge_; };

    //! return global slave dof row map associated with surface nodes
    Teuchos::RCP<Epetra_Map>& GSDofSurfRowMapPtr() { return gsdofSurf_; };
    Teuchos::RCP<const Epetra_Map> GSDofSurfRowMapPtr() const { return gsdofSurf_; };

    //! return global LM dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map>& PGLmDofRowMapPtr() { return pglmdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> PGLmDofRowMapPtr() const { return pglmdofrowmap_; };

    //! return global slave dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map>& PGSlDofRowMapPtr() { return pgsdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> PGSlDofRowMapPtr() const { return pgsdofrowmap_; };

    //! return global master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map>& PGMaDofRowMapPtr() { return pgmdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> PGMaDofRowMapPtr() const { return pgmdofrowmap_; };

    //! return global slave and master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map>& PGSlMaDofRowMapPtr() { return pgsmdofrowmap_; };
    Teuchos::RCP<const Epetra_Map> PGSlMaDofRowMapPtr() const { return pgsmdofrowmap_; };

    //! return global dirichlet toggle of all slave dofs (before parallel redistribution)
    Teuchos::RCP<Epetra_Vector>& PGSlDirichToggleDofRowMapPtr() { return pgsdirichtoggle_; };
    Teuchos::RCP<const Epetra_Vector> PGSlDirichToggleDofRowMapPtr() const
    {
      return pgsdirichtoggle_;
    };

    //! return initial col ele map for binning strategy (s m)
    std::vector<Teuchos::RCP<Epetra_Map>>& InitialSlMaEleColMap() { return initial_elecolmap_; };
    const std::vector<Teuchos::RCP<Epetra_Map>>& InitialSlMaEleColMap() const
    {
      return initial_elecolmap_;
    };

    //! return global Mortar matrix D
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& DMatrixPtr() { return dmatrix_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> DMatrixPtr() const { return dmatrix_; };
    CORE::LINALG::SparseMatrix& DMatrix()
    {
      if (dmatrix_.is_null()) dserror("The dmatrix_ is not initialized!");
      return *dmatrix_;
    }

    //! return global Mortar matrix M
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& MMatrixPtr() { return mmatrix_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> MMatrixPtr() const { return mmatrix_; };
    CORE::LINALG::SparseMatrix& MMatrix()
    {
      if (mmatrix_.is_null()) dserror("The mmatrix_ is not initialized!");
      return *mmatrix_;
    }

    //! return global weighted gap vector g
    Teuchos::RCP<Epetra_Vector>& WGapPtr() { return g_; };
    Teuchos::RCP<const Epetra_Vector> WGapPtr() const { return g_; };
    Epetra_Vector& WGap()
    {
      if (g_.is_null()) dserror("The wGapRhsPtr_ is not initialized!");
      return *g_;
    }

    //! return global tangential rhs vector
    Teuchos::RCP<Epetra_Vector>& TangRhsPtr() { return tangrhs_; };
    Teuchos::RCP<const Epetra_Vector> TangRhsPtr() const { return tangrhs_; };

    //! return gloabl inactive rhs vector
    Teuchos::RCP<Epetra_Vector>& InactiveRhsPtr() { return inactiverhs_; };
    Teuchos::RCP<const Epetra_Vector> InactiveRhsPtr() const { return inactiverhs_; };
    Epetra_Vector& InactiveRhs()
    {
      if (inactiverhs_.is_null()) dserror("The inactiverhs_ is not initialized!");
      return *inactiverhs_;
    }

    //! Return the structural contact right-hand-side contributions of the current time step
    //! \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& StrContactRhsPtr() { return strContactRhsPtr_; }
    Teuchos::RCP<const Epetra_Vector> StrContactRhsPtr() const { return strContactRhsPtr_; }
    Epetra_Vector& StrContactRhs()
    {
      if (strContactRhsPtr_.is_null()) dserror("The strContactRhsPtr_ is not initialized!");
      return *strContactRhsPtr_;
    }

    //! return global constraint rhs vector (only for saddlepoint problems)
    Teuchos::RCP<Epetra_Vector>& ConstrRhsPtr() { return constrrhs_; };
    Teuchos::RCP<const Epetra_Vector> ConstrRhsPtr() const { return constrrhs_; };
    Epetra_Vector& ConstrRhs()
    {
      if (constrrhs_.is_null()) dserror("The constrrhs_ is not initialized!");
      return *constrrhs_;
    }

    //! return global Matrix LinD containing slave fc derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& DLinMatrixPtr() { return lindmatrix_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> DLinMatrixPtr() const { return lindmatrix_; };
    CORE::LINALG::SparseMatrix& DLinMatrix()
    {
      if (lindmatrix_.is_null()) dserror("The augDnLinMatrixPtr_ is not initialized!");
      return *lindmatrix_;
    }

    //! return global Matrix LinM containing master fc derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& MLinMatrixPtr() { return linmmatrix_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> MLinMatrixPtr() const { return linmmatrix_; };
    CORE::LINALG::SparseMatrix& MLinMatrix()
    {
      if (linmmatrix_.is_null()) dserror("The augMnLinMatrixPtr_ is not initialized!");
      return *linmmatrix_;
    }

    //! return global Matrix kteffnew containing modified jacobian
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& kteffnewMatrixPtr() { return kteffnew_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> kteffnewMatrixPtr() const { return kteffnew_; };
    CORE::LINALG::SparseMatrix& kteffnewMatrix()
    {
      if (kteffnew_.is_null()) dserror("The kteffnewMatrixPtr is not initialized!");
      return *kteffnew_;
    }

    //! return global Mortar matrix D (last end-point \f$t_{n}\f$)
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& OldDMatrixPtr() { return dold_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> OldDMatrixPtr() const { return dold_; };

    //! return global Mortar matrix M (last end-point \f$t_{n}\f$)
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& OldMMatrixPtr() { return mold_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> OldMMatrixPtr() const { return mold_; };

    //! return current vector of Lagrange multipliers at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& LmPtr() { return z_; };
    Teuchos::RCP<const Epetra_Vector> LmPtr() const { return z_; };

    //! return old vector of Lagrange multipliers at \f$t_{n}\f$
    Teuchos::RCP<Epetra_Vector>& OldLmPtr() { return zold_; };
    Teuchos::RCP<const Epetra_Vector> OldLmPtr() const { return zold_; };

    /*! \brief Return Lagrange multiplier vector increment
     *
     *  \remark This is NOT the increment of z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!) */
    Teuchos::RCP<Epetra_Vector>& LmIncrPtr() { return zincr_; };
    Teuchos::RCP<const Epetra_Vector> LmIncrPtr() const { return zincr_; };

    //! return vector of Lagrange multipliers from last Uzawa step
    Teuchos::RCP<Epetra_Vector>& LmUzawaPtr() { return zuzawa_; };
    Teuchos::RCP<const Epetra_Vector> LmUzawaPtr() const { return zuzawa_; };

    //! return vector of normal contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& StressNormalPtr() { return stressnormal_; };
    Teuchos::RCP<const Epetra_Vector> StressNormalPtr() const { return stressnormal_; };

    //! return vector of tangential contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& StressTangentialPtr() { return stresstangential_; };
    Teuchos::RCP<const Epetra_Vector> StressTangentialPtr() const { return stresstangential_; };

    //! return vector of normal contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& ForceNormalPtr() { return forcenormal_; };
    Teuchos::RCP<const Epetra_Vector> ForceNormalPtr() const { return forcenormal_; };

    //! return vector of tangential contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& ForceTangentialPtr() { return forcetangential_; };
    Teuchos::RCP<const Epetra_Vector> ForceTangentialPtr() const { return forcetangential_; };

    //! return time step index at \f$t_{n+1}\f$
    int& StepNp() { return stepnp_; };
    int StepNp() const { return stepnp_; };

    //! return non-linear (Newton) iteration index
    int& NlnIter() { return iter_; };
    int NlnIter() const { return iter_; };

    //! return flag indicating global contact status
    bool& IsInContact() { return isincontact_; };
    bool IsInContact() const { return isincontact_; };

    //! return flag indicating global contact status of this time step (history)
    bool& WasInContact() { return wasincontact_; };
    bool WasInContact() const { return wasincontact_; };

    //! return flag indicating global contact status of last time step
    bool& WasInContactLastTimeStep() { return wasincontactlts_; };
    bool WasInContactLastTimeStep() const { return wasincontactlts_; };

    //! return flag indicating potential self contact
    bool& IsSelfContact() { return isselfcontact_; };
    bool IsSelfContact() const { return isselfcontact_; };

    //! return flag for frictional contact
    bool& IsFriction() { return friction_; };
    bool IsFriction() const { return friction_; };

    //! return flag for nonsmooth contact
    bool& IsNonSmoothContact() { return nonSmoothContact_; };
    const bool& IsNonSmoothContact() const { return nonSmoothContact_; };

    //! return flag for regularized contact
    bool& IsRegularized() { return regularized_; };
    bool IsRegularized() const { return regularized_; };

    //! return flag indicating whether trafo should be applied
    bool& IsDualQuadSlaveTrafo() { return dualquadslavetrafo_; };
    bool IsDualQuadSlaveTrafo() const { return dualquadslavetrafo_; };

    //! return transformation matrix T for dual quad 3D case
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& TrafoPtr() { return trafo_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> TrafoPtr() const { return trafo_; };

    //! return inverse trafo matrix T^(-1) for dual quad 3D case
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& InvTrafoPtr() { return invtrafo_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> InvTrafoPtr() const { return invtrafo_; };

    //! return modified global Mortar matrix D
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& ModifiedDMatrixPtr() { return dmatrixmod_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> ModifiedDMatrixPtr() const
    {
      return dmatrixmod_;
    };

    //! return modified global Mortar matrix Dold
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& OldModifiedDMatrixPtr() { return doldmod_; };
    Teuchos::RCP<const CORE::LINALG::SparseMatrix> OldModifiedDMatrixPtr() const
    {
      return doldmod_;
    };

    //! return integration time
    double& IntTime() { return inttime_; };
    double IntTime() const { return inttime_; };

    //! return mean interface velocity
    std::vector<double>& MeanInterfaceVels() { return ivel_; };
    const std::vector<double>& MeanInterfaceVels() const { return ivel_; };

    //! return current used solving strategy
    INPAR::CONTACT::SolvingStrategy& SolType() { return stype_; };
    INPAR::CONTACT::SolvingStrategy SolType() const { return stype_; };

    //! return direction in which the contact constraints are formulated
    INPAR::CONTACT::ConstraintDirection& ConstrDirection() { return constr_direction_; };
    INPAR::CONTACT::ConstraintDirection ConstrDirection() const { return constr_direction_; };

    INPAR::MORTAR::ParallelRedist& ParType() { return partype_; };
    INPAR::MORTAR::ParallelRedist ParType() const { return partype_; };

    //!@}

   private:
    //! global Lagrange multiplier dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> glmdofrowmap_;

    //! global reference dof row map for self contact Lagrange multipliers (of all interfaces)
    Teuchos::RCP<Epetra_Map> gscrefdofrowmap_;

    //! global Lagrange mult. dof row map for self contact (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsclmdofrowmap_;

    //! global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsnoderowmap_;

    //! global master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmnoderowmap_;

    //! global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsdofrowmap_;

    //! global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmdofrowmap_;

    //! global internal dof row map
    Teuchos::RCP<Epetra_Map> gndofrowmap_;

    //! global slave and master dof row map (s+m map)
    Teuchos::RCP<Epetra_Map> gsmdofrowmap_;

    //! global displacement dof row map (s+m+n map)
    Teuchos::RCP<Epetra_Map> gdisprowmap_;

    //! @name Active set
    //!@{

    //! global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gactivenodes_;

    //! global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gactivedofs_;

    //! global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> ginactivenodes_;

    //! global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> ginactivedofs_;

    //! global active slave dof row map in normal direction (of all interfaces)
    Teuchos::RCP<Epetra_Map> gactiven_;

    //! global dof row map of matrix T (of all interfaces)
    Teuchos::RCP<Epetra_Map> gactivet_;

    //! global slip slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gslipnodes_;

    //! global slip slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gslipdofs_;

    //! global slip slave dof row map in tangential direction (of all interfaces)
    Teuchos::RCP<Epetra_Map> gslipt_;

    //!@}

    //! global slave dof row map of vertex nodes
    Teuchos::RCP<Epetra_Map> gsdofVertex_;

    //! global slave dof row map of edge nodes
    Teuchos::RCP<Epetra_Map> gsdofEdge_;

    //! global slave dof row map of surface nodes
    Teuchos::RCP<Epetra_Map> gsdofSurf_;

    //! @name Parallel redistribution
    //!@{

    /*! Max-to-min ratio of evaluation time across all processes for currnet time step \f$t_{n+1}\f$
     */
    std::vector<double> unbalanceEvaluationTime_;

    /*! Max-to-min ratio of number of row slave elements across all processes for current time step
     * \f$t_{n+1}\f$
     */
    std::vector<int> unbalanceNumSlaveElements_;

    //! global LM dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pglmdofrowmap_;

    //! global slave dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsdofrowmap_;

    //! global master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgmdofrowmap_;

    //! global slave and master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsmdofrowmap_;

    //! global dirichlet toggle of all slave dofs (before parallel redistribution)
    Teuchos::RCP<Epetra_Vector> pgsdirichtoggle_;

    //! parallel redistribution type
    INPAR::MORTAR::ParallelRedist partype_;

    //!@}

    //! @name Binning strategy
    //!@{

    //! initial col ele map for binning strategy (s m)
    std::vector<Teuchos::RCP<Epetra_Map>> initial_elecolmap_;

    //!@}

    //! global Mortar matrix \f$D\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dmatrix_;

    //! global Mortar matrix \f$M\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix> mmatrix_;

    //! global weighted gap vector \f$g\f$
    Teuchos::RCP<Epetra_Vector> g_;

    //! global tangential right-hand side vector (formulation with incremental #z_)
    Teuchos::RCP<Epetra_Vector> tangrhs_;

    /*! \brief Gloabl inactive right-hand side vector
     *
     * This is used for the formulation with incremental #z_ and saddle point system.
     */
    Teuchos::RCP<Epetra_Vector> inactiverhs_;

    //! structural contact right-hand-side vector at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> strContactRhsPtr_;

    //! global constraint right-hand side vector (only for saddlepoint problems)
    Teuchos::RCP<Epetra_Vector> constrrhs_;

    //! global Matrix LinD containing slave fc derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix> lindmatrix_;

    //! global Matrix LinM containing master fc derivatives
    Teuchos::RCP<CORE::LINALG::SparseMatrix> linmmatrix_;

    //! global K matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> kteffnew_;

    //! global Mortar matrix D (last end-point \f$t_{n}\f$)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dold_;

    //! global Mortar matrix M (last end-point \f$t_{n}\f$)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> mold_;

    //! current vector of Lagrange multipliers at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> z_;

    //! old vector of Lagrange multipliers at \f$t_{n}\f$
    Teuchos::RCP<Epetra_Vector> zold_;

    /*! \brief Lagrange multiplier vector increment within SaddlePointSolve
     *
     *  \remark This is \em not the increment of #z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!)
     */
    Teuchos::RCP<Epetra_Vector> zincr_;

    //! vector of Lagrange multipliers from last Uzawa step
    Teuchos::RCP<Epetra_Vector> zuzawa_;

    //! vector of normal contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> stressnormal_;

    //! vector of tangential contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> stresstangential_;

    //! vector of normal contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> forcenormal_;

    //! vector of tangential contact forces at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> forcetangential_;

    //! @name Counters and indices
    //!@{

    //! time step index at \f$t_{n+1}\f$
    int stepnp_;

    //! Nonlinear iteration index, e.g. Newton iteration
    int iter_;

    //!@}

    //! @name Status flags
    //!@{

    //! flag indicating global contact status
    bool isincontact_;

    //! flag indicating global contact status of this time step (history)
    bool wasincontact_;

    //! flag indicating global contact status of last time step
    bool wasincontactlts_;

    //! flag indicating potential self contact
    bool isselfcontact_;

    //! flag for frictional contact
    bool friction_;

    //! flag for non-smooth contact
    bool nonSmoothContact_;

    //! flag for regularized contact
    bool regularized_;

    //! flag indicating whether trafo should be applied
    bool dualquadslavetrafo_;

    //!@}

    //! transformation matrix T for dual quad 3D case
    Teuchos::RCP<CORE::LINALG::SparseMatrix> trafo_;

    //! inverse trafo matrix T^(-1) for dual quad 3D case
    Teuchos::RCP<CORE::LINALG::SparseMatrix> invtrafo_;

    //! modified global Mortar matrix D
    Teuchos::RCP<CORE::LINALG::SparseMatrix> dmatrixmod_;

    //! modified global Mortar matrix Dold
    Teuchos::RCP<CORE::LINALG::SparseMatrix> doldmod_;

    /*! \brief Integration time
     *
     * \todo Is this the wall clock time required to perform the mortar integration?
     */
    double inttime_;

    //! mean interface velocity
    std::vector<double> ivel_;

    //! current used solving strategy
    INPAR::CONTACT::SolvingStrategy stype_;

    //! direction in which the contact constraints are formulated
    INPAR::CONTACT::ConstraintDirection constr_direction_;

  };  // class AbstractStratDataContainer


  /*! \brief Main abstract class for contact solution strategies
   *
   * This is the templating abstract class for all contact solution algorithms.
   * Every solution algorithm has to fit into the set of functions and calls defined herein
   * and has to be specified in a corresponding subclass defining the concrete algorithmic steps.
   *
   * This class it itself derived from the MORTAR::StrategyBase class, which is an even
   * more abstract framework for any solution strategies involving mortar coupling.
   *
   * \remark Please add no new member variables to the abstract strategy! Use
   * the corresponding data container instead (--> CONTACT::AbstractStratDataContainer).
   *
   * Refer also to the Semesterarbeit of Bernd Budich, 2009
   *
   */
  class AbstractStrategy : public MORTAR::StrategyBase
  {
   public:
    /*!
    \brief Standard constructor

    This constructor uses the given DataContainer to store and share all its
    member variables. The declared member variables are just references to
    the container content!

    Creates the strategy base object and initializes all global variables.


    \param[in] stratData Data container object
    \param[in] DofRowMap Dof row map of underlying problem
    \param[in] NodeRowMap Node row map of underlying problem
    \param[in] params List of contact/parameters
    \param[in] spatialDim Spatial dimension of the problem
    \param[in] comm Communicator
    \param[in] alphaf Mid-point for Generalized-alpha time integration
    \param[in] maxdof Highest DOF number in global problem
    */
    AbstractStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params, const int spatialDim,
        const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, int const maxdof);

    /*! \brief Setup this strategy object (maps, vectors, etc.)

     All global maps and vectors are initialized by collecting
     the necessary information from all interfaces. In the case
     of a parallel redistribution, this method is called again
     to re-setup the above mentioned quantities. In this case
     we set the input parameter redistributed=TRUE. Moreover,
     when called for the first time (in the constructor) this
     method is given the input parameter init=TRUE to account
     for initialization of the active set. */
    virtual void Setup(bool redistributed, bool init);


    //! return the current solution type
    virtual INPAR::CONTACT::SolvingStrategy Type() const { return stype_; }

    //! @name Access methods
    //!@{

    //! Return the NOX::NLN::CONSTRAINT::Interface::Required member object
    const Teuchos::RCP<CONTACT::NoxInterface>& NoxInterfacePtr() { return noxinterface_ptr_; };

    /*! \brief Return the Lagrange multiplier dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> LMDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist()) return Data().PGLmDofRowMapPtr();

      return Data().GLmDofRowMapPtr();
    };
    virtual const Epetra_Map& LMDoFRowMap(const bool& redist) const
    {
      return *LMDoFRowMapPtr(redist);
    }

    /*! \brief Return the Lagrange multiplier dof row map for the global linear
     *  system
     *
     *  \note This map is NOT used internally. Its only purpose is to provide a
     *  map as meaningful upper bound for potentially acquired LM dofs.
     *
     *  \date 04/2018
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> LinSystemLMDoFRowMapPtr() const
    {
      if (SystemType() != INPAR::CONTACT::system_saddlepoint) return Teuchos::null;

      if (IsSelfContact())
      {
        if (ParRedist()) dserror("Parallel redistribution is not supported for self contact!");
        return Data().GSelfContactLmDofRowMapPtr();
      }
      else
        return LMDoFRowMapPtr(false);
    };
    virtual const Epetra_Map& LinSystemLMDoFRowMap() const { return *LinSystemLMDoFRowMapPtr(); }

    /*! \brief Return the slave dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> SlDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist()) return Data().PGSlDofRowMapPtr();

      return Data().GSlDofRowMapPtr();
    };
    virtual const Epetra_Map& SlDoFRowMap(const bool& redist) const
    {
      return *SlDoFRowMapPtr(redist);
    }

    /*! \brief Return the slave dof row map in normal direction
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> SlNormalDoFRowMapPtr(const bool& redist) const
    {
      dserror("Map not available in abstract strategy!");
      if ((not redist) and ParRedist())
        dserror("The original / not redistributed slave normal row map is not available!");

      return Teuchos::null;
    };
    virtual const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const
    {
      // currently not supported for the abstract strategy
      dserror("SlNormalDoFRowMap() seems currently unsupported!");
      exit(EXIT_FAILURE);
    }

    /*! \brief Return the slave dof row map in the tangential directions
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> SlTangentialDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist())
        dserror("The original / not redistributed slave tangential row map is not available!");

      return Teuchos::null;
    };
    virtual const Epetra_Map& SlTangentialDoFRowMap(const bool& redist) const
    {
      return *gslipdofs_;
    }

    /*! \brief Return the master dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> MaDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist()) return Data().PGMaDofRowMapPtr();

      return Data().GMaDofRowMapPtr();
    };
    virtual const Epetra_Map& MaDoFRowMap(const bool& redist) const
    {
      return *MaDoFRowMapPtr(redist);
    }

    /*! \brief Return the combined slave/master dof row map
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> SlMaDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist()) return Data().PGSlMaDofRowMapPtr();

      return Data().GSlMaDofRowMapPtr();
    };
    virtual const Epetra_Map& SlMaDoFRowMap(const bool& redist) const
    {
      return *SlMaDoFRowMapPtr(redist);
    }


    /*! \brief Return the desired right-hand-side block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum CONTACT::VecBlockType& bt) const
    {
      dserror("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    /*! \brief Return the desired right-hand side block pointer for norm check
     *  (read-only)
     *
     *  In the default case this method returns the standard right-hand side block,
     *  i.e. the same as for the assembly procedure. Anyway, in some cases it is
     *  meaningful to use a modified right-hand side, e.g. without penalty
     *  contributions in an augmented framework.
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired vector block type, e.g. displ, constraint, ...
     *
     *  \author hiermeier \date 08/17  */
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtrForNormCheck(
        const enum CONTACT::VecBlockType& bt) const
    {
      return GetRhsBlockPtr(bt);
    }

    /*! Return the condensed right-hand-side (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Vector> GetCondensedRhsPtr(
        Epetra_Vector& f, const double& timefac_np) const
    {
      dserror("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    /*! \brief Return the desired matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired matrix block type, e.g. displ_displ, displ_lm, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum CONTACT::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = nullptr) const
    {
      dserror("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    //! Apply modifications (e.g. condensation) directly before linear solve
    virtual void RunPreApplyJacobianInverse(
        Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff, Epetra_Vector& rhs)
    { /* do nothing */
    }

    /*! Return the condensed matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     */
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> GetCondensedMatrixBlockPtr(
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& kteff, const double& timefac_np) const
    {
      dserror("Not yet implemented!");
      exit(EXIT_FAILURE);

      return Teuchos::null;
    };

    //! Return global slave node row map
    Teuchos::RCP<Epetra_Map> SlaveRowNodes() override { return Data().GSlNodeRowMapPtr(); }
    Teuchos::RCP<const Epetra_Map> SlRowNodesPtr() const { return Data().GSlNodeRowMapPtr(); }
    const Epetra_Map& SlRowNodes() const { return *Data().GSlNodeRowMapPtr(); }

    //! Return global slave node row map
    Teuchos::RCP<const Epetra_Map> MaRowNodesPtr() const { return Data().GMaNodeRowMapPtr(); }
    const Epetra_Map& MaRowNodes() const { return *Data().GMaNodeRowMapPtr(); }

    //! Return global active node row map
    Teuchos::RCP<Epetra_Map> ActiveRowNodes() override { return Data().GActiveNodeRowMapPtr(); };
    virtual Teuchos::RCP<const Epetra_Map> ActiveRowNodes() const
    {
      return Data().GActiveNodeRowMapPtr();
    };

    //! Return global slip node row map
    Teuchos::RCP<Epetra_Map> SlipRowNodes() override { return Data().GSlipNodeRowMapPtr(); };
    Teuchos::RCP<const Epetra_Map> SlipRowNodes() const { return Data().GSlipNodeRowMapPtr(); };

    //! Return global slave dof row map
    Teuchos::RCP<Epetra_Map> SlaveRowDofs() { return Data().GSlDofRowMapPtr(); }

    //! Return global active dof row map
    Teuchos::RCP<Epetra_Map> ActiveRowDofs() override { return Data().GActiveDofRowMapPtr(); }

    //! Return global master dof row map
    Teuchos::RCP<Epetra_Map> MasterRowDofs() { return Data().GMaDofRowMapPtr(); }

    //! Return global slave dof row map
    Teuchos::RCP<Epetra_Map> SlaveMasterRowDofs() { return Data().GSlMaDofRowMapPtr(); }

    //! Return not redistributed global slave dof row map
    Teuchos::RCP<Epetra_Map> NotReDistSlaveRowDofs() override { return Data().PGSlDofRowMapPtr(); }

    //! Return not redistributed global master dof row map
    Teuchos::RCP<Epetra_Map> NotReDistMasterRowDofs() override { return Data().PGMaDofRowMapPtr(); }

    /*!
    \brief Gather maps needed for contact/meshtying specific multigrid preconditioners

    @param MasterDofMap Dof row map of master interface
    @param SlaveDofMap Dof row map of slave interface
    @param InnerDofMap Dof row map of interior volume
    @param ActiveDofMap Dof row map of active slave contact interface
    */
    void CollectMapsForPreconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap,
        Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap,
        Teuchos::RCP<Epetra_Map>& ActiveDofMap) const override;

    //! Return Lagrange multiplier vector (\f$t_{n+1}\f$)
    Teuchos::RCP<Epetra_Vector> LagrMult() override { return z_; }

    /*! \brief Return Lagrange multiplier vector \f$(t_{n+1})\f$
     *
     *  \param redist (in): If TRUE, the redistributed vector is returned,
     *                      otherwise the vector with the original map before
     *                      any redistribution took place.
     *
     *  \warning The vector is returned with the slave dof row map, i.e. actually the wrong map!
     *
     *  \author hiermeier
     *  \date 05/16 */
    virtual Teuchos::RCP<const Epetra_Vector> GetLagrMultNp(const bool& redist) const;

    //! Return old Lagrange multiplier vector (\f$t_{n}\f$)
    Teuchos::RCP<Epetra_Vector> LagrMultOld() override { return Data().OldLmPtr(); }

    /*! \brief Return old Lagrange multiplier vector \f$(t_n)\f$
     *
     *  \param redist (in): If TRUE, the redistributed vector is returned,
     *                      otherwise the vector with the original map before
     *                      any redistribution took place.
     *
     *  \warning The vector is returned with the slave dof row map, i.e. actually the wrong map!
     *
     *  \author hiermeier
     *  \date 05/16 */
    virtual Teuchos::RCP<const Epetra_Vector> GetLagrMultN(const bool& redist) const;

    //! Return Lagrange multiplier vector from last Uzawa step
    Teuchos::RCP<Epetra_Vector> LagrMultUzawa() { return Data().LmUzawaPtr(); }

    //! Return constraint rhs vector (only in saddle-point formulation
    Teuchos::RCP<Epetra_Vector> ConstrRhs() override { return Data().ConstrRhsPtr(); }

    //! Returns increment of LagrangeMultiplier solution vector in SaddlePointSolve routine
    Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() override { return Data().LmIncrPtr(); }
    Teuchos::RCP<const Epetra_Vector> GetLagrMultSolveIncr() const { return Data().LmIncrPtr(); };

    //! Return mortar matrix D
    Teuchos::RCP<CORE::LINALG::SparseMatrix> DMatrix() override { return Data().DMatrixPtr(); }

    //! Return mortar matrix M
    Teuchos::RCP<CORE::LINALG::SparseMatrix> MMatrix() override { return Data().MMatrixPtr(); }

    //! Return vector of normal contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactNorStress() override { return Data().StressNormalPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactNorStress() const { return Data().StressNormalPtr(); }
    //! Return weighted gap
    Teuchos::RCP<Epetra_Vector> ContactWGap() { return Data().WGapPtr(); }

    //! Return vector of tangential contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactTanStress() override { return Data().StressTangentialPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactTanStress() const
    {
      return Data().StressTangentialPtr();
    }

    //! Return vector of normal contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactNorForce() override { return Data().ForceNormalPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactNorForce() const { return Data().ForceNormalPtr(); }

    //! Return vector of tangential contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactTanForce() override { return Data().ForceTangentialPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactTanForce() const
    {
      return Data().ForceTangentialPtr();
    }


    //! Return required Integration time
    double Inttime() override { return Data().IntTime(); };

    //! Set integration time to zero
    void Inttime_init() override { Data().IntTime() = 0.0; };

    //! Return current global contact status
    bool IsInContact() const override { return Data().IsInContact(); }

    /*! \brief Return old global contact status (this time step)

     True if there has been contact in any nonlinear iteration
     step of the current time step. */
    bool WasInContact() const override { return Data().WasInContact(); }

    /*!
    \brief Return old global contact status (last time step)

    True if there has been contact at the end of the last
    time step (last converged state)
    */
    bool WasInContactLastTimeStep() const override { return Data().WasInContactLastTimeStep(); }

    /*!
    \brief Return global self contact status

    Note that at the moment this only gives information about the
    POTENTIAL self contact of the global problem and not about
    an actual self contact occurring.

    TODO: automatically recognize ACTUAL self contact
    */
    bool& IsSelfContact() { return Data().IsSelfContact(); }
    bool IsSelfContact() const { return Data().IsSelfContact(); };

    //! Return global frictional status
    bool Friction() const override { return Data().IsFriction(); }

    //! Return contact interfaces
    const std::vector<Teuchos::RCP<CONTACT::Interface>>& ContactInterfaces() const
    {
      return Interfaces();
    }

    /*!
    \brief Get dual quadratic 3d slave element flag

    Returns TRUE if at least one higher-order 3d slave element with
    dual Lagrange mutliplier shape functions in any interface.
    */
    virtual bool Dualquadslavetrafo() const { return Data().IsDualQuadSlaveTrafo(); };

    //! Return parallel redistribution status (yes or no)
    inline bool ParRedist() const
    {
      return (Data().ParType() != INPAR::MORTAR::ParallelRedist::redist_none);
    }


    //! Return specific parallel redistribution status
    inline INPAR::MORTAR::ParallelRedist WhichParRedist() const { return Data().ParType(); }

    //! Return matrix T
    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> TMatrix() { return Teuchos::null; }

    //! Return number of active nodes
    int NumberOfActiveNodes() const override
    {
      if (not Data().GActiveNodeRowMapPtr().is_null())
        return Data().GActiveNodeRowMapPtr()->NumGlobalElements();
      return 0;
    }

    //! Return number of frictional slip nodes
    int NumberOfSlipNodes() const override
    {
      if (not Data().GSlipNodeRowMapPtr().is_null())
        return Data().GSlipNodeRowMapPtr()->NumGlobalElements();
      return 0;
    }

    //!@}

    //! @name Parallel redistribution
    //!@{

    /*!
    \brief Redistribute all contact interfaces in parallel

    We have two code paths to perform contact load balancing:
    - Using RedistributeWithSafeGhosting() will guarantee, that the master-sided ghosting is
    sufficiently far and no master elements will be missed in the subsequent contact search.
    Applicability of this code path is limited to some contact scenarios.
    - RedistributeContactOld() provides the legacy implementation to be used with all specialized
    contact features. However, master-sided interface ghosting might be insufficient.

    \post Each contact interface is FillComplete().

    \param[in] dis Current displacement state
    \param[in] vel Current velocity state

    \return TRUE if the interface has been redistributed. Return FALSE otherwise.
    */
    bool RedistributeContact(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel) override;

    /** \brief Redistribute all contact interfaces in parallel
     *
     *  In contrast to RedistributeContact this routine takes place at a different
     *  point during the simulation. For example, the redistribution can be initiated
     *  each time a certain amount of Newton steps per load step has been reached.
     *  In this way an adaption can be made quicker directly after a large predictor
     *  step or another unforeseen scenario which might have changed the contact
     *  situation severely. */
    virtual bool DynRedistributeContact(const Teuchos::RCP<const Epetra_Vector>& dis,
        Teuchos::RCP<const Epetra_Vector> vel, const int nlniter)
    {
      return false;
    };

    //!@}

    //! @name Evaluation methods
    //! @{

    /*!
    \brief Global evaluation method called from time integrator

    This routine handles the evaluation of all contact terms. This is time consuming. To assess
    timing, detailed time measurements can be performed. This requires synchronization among all MPI
    ranks. By default, detailed time measurements are turned off.

    @param dis Current displacement state
    @param[in/out] kt Global Jacobian matrix
    @param[in/out] f Global residual vector
    @param[in] timeStep Current time step
    @param[in] nonlinearIteration Current nonlinear iteration step
    @param[in] predictor Is this called during the predictor?
    */
    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
        const int timeStep, const int nonlinearIteration, bool predictor = false) override;

    /*! \brief Reset the internal state variables
     *
     *  \date 02/2016
     *  \author hiermeier */
    virtual void Reset(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp,
        const Epetra_Vector& xnew);

    /*! \brief Global evaluation method called from STR::MODELEVALUATOR::Contact class
     *
     *  \date 03/2016
     *  \author hiermeier */
    void Evaluate(CONTACT::ParamsInterface& cparams) { Evaluate(cparams, nullptr); }

    /*! \brief Global evaluation method called from STR::MODELEVALUATOR::Contact class
     *
     *  \date 03/2016
     *  \author hiermeier */

    void Evaluate(CONTACT::ParamsInterface& cparams,
        const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec)
    {
      Evaluate(cparams, eval_vec, nullptr);
    }

    /*! \brief Global evaluation method called from STR::MODELEVALUATOR::Contact class
     *
     * This is the central place to enter contact evaluation.
     * The actual evaluation operation is governed by the MORTAR::ActionType in the
     * CONTACT::ParamsInterface. We use a switch on the MORTAR::ActionType to call the actual
     * evaluation routine.
     *
     * \note This routine is \em not virtual as it is not supposed to be overloaded.
     *
     * \date 03/2016
     * \author hiermeier */
    void Evaluate(CONTACT::ParamsInterface& cparams,
        const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec,
        const std::vector<Teuchos::RCP<Epetra_Vector>>* eval_vec_mutable);

    /*! \brief Set current deformation state

    All interfaces are called to set the current deformation state
    (u, xspatial) in their nodes. Additionally, the new contact
    element areas are computed.

    \param statename (in): std::string defining which quantity to set (either "displacement" or
                           "olddisplacement")
    \param vec (in): current global state of the quantity defined by statename
    */
    void SetState(const enum MORTAR::StateType& statetype, const Epetra_Vector& vec) override;

    /*! \brief Evaluate reference state

     for frictional contact we need history values (relative velocity) and
     therefore we store the nodal entries of mortar matrices (reference
     configuration) before the first time step

     \pre SetState() has been called.
     */
    void EvaluateReferenceState() override;

    /*! \brief Evaluate matrix of nodal normals

     This is needed for energy-conserving time integration (Velocity-Update) */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> EvaluateNormals(
        Teuchos::RCP<Epetra_Vector> dis) override;

    //!@}

    //! @name Merit function methods
    //!@{

    /// return the potential contributions of the active contact strategy
    virtual double GetPotentialValue(
        const enum NOX::NLN::MeritFunction::MeritFctName mrt_type) const;

    /// return contributions of the active contact strategy to the linear model
    virtual double GetLinearizedPotentialValueTerms(const Epetra_Vector& dir,
        const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
        const enum NOX::NLN::MeritFunction::LinOrder linorder,
        const enum NOX::NLN::MeritFunction::LinType lintype) const;

    //!@}

    //! @name Preconditioner methods
    //!@{

    //! Is this a saddle-point system?
    bool IsSaddlePointSystem() const override;

    //! Is this a condensed system?
    bool IsCondensedSystem() const override;

    /*! \brief Fill the maps vector for the linear solver preconditioner

    The following order is pre-defined:
    (0) masterDofMap
    (1) slaveDofMap
    (2) innerDofMap
    (3) activeDofMap

    \author hiermeier
    */
    void FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map>>& maps) const override;

    //! compute the preconditioner operator
    bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
        Teuchos::ParameterList* precParams = nullptr) override;

    //!@}

    //! @name Quantity control methods
    //!@{

    /*! \brief Get some nodal quantity globally and store into Nodes

     The enum input parameter defines, which quantity is be updated.
     Currently the possibilities "lmold", "lmcurrent", "lmupdate" and
     "lmuzawa" exist. Note that "lmold" means the converged value LM_n
     of the last time / load step, whereas "lmcurrent" adresses the current
     (not necessarily converged) value of the LM_n+1. "lmupdate" is a special
     option called only in Recover() after the update of the Lagr. multipliers.
     It basically does the same as "lmcurrent", but also checks for D.B.C.
     problems. Finally, "lmuzawa" addresses the LM update within an
     Uzawa augmented Lagrangian scheme.

     \param type (in): enum defining which quantity to store into Nodes

     */
    void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type) override;

    /*! \brief Evaluate contact stresses in normal direction and tangential plane

     This is called at the end of each time or load step. It calculates
     the stress vector in normal direction and the stress vector in the
     tangential plane. */
    void ComputeContactStresses() override;

    /*! \brief Get dirichlet B.C. status and store into Nodes

     This is called once at the beginning of the simulation
     to set the D.B.C. status in each CNode.

     \param dbcmaps (in): MapExtractor carrying global dbc map */
    void StoreDirichletStatus(Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmaps) override;

    virtual void SetParentState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec,
        const Teuchos::RCP<DRT::Discretization> dis){
        /* standard contact methods don't need the corresponding bulk element */};

    /*! \brief Update contact at end of time step

     \param dis (in):  current displacements (-> old displacements)

     */
    void Update(Teuchos::RCP<const Epetra_Vector> dis) override;

    /*! \brief Perform a write restart

     A write restart is initiated by the contact manager. However, the manager has no
     direct access to the nodal quantities. Hence, a portion of the restart has to be
     performed on the level of the contact algorithm, for short: here's the right place.
     */
    void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const override;

    /*!
    \brief Read restart data from disk

    @param reader Discretization reader to be used for reading the restart data
    @param dis Displacement vector of the solid field
    */
    void DoReadRestart(
        IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis) override
    {
      DoReadRestart(reader, dis, Teuchos::null);
    };

    /*!
    \brief Read restart data from disk

    @param reader Discretization reader to be used for reading the restart data
    @param dis Displacement vector of the solid field
    @param cparams_ptr ??
    */
    virtual void DoReadRestart(IO::DiscretizationReader& reader,
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    //!@}

    //! @name Output
    //!@{

    /*! \brief Write strategy specific output
     *
     *  \param(in) writer: output writer */
    virtual void WriteOutput(IO::DiscretizationWriter& writer) const { return; }

    /*! \brief Compute interface forces and moments
     *
     * Compute current interface forces and moments at n+1-alphaf using current
     * Lagrange multiplier values and current Mortar matrices D and M at n+1. When
     * doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
     * matrices of the last converged time / load step n (TR-like interpolation).
     *
     *\param output (in): flag indicating whether force output shall be written
     */
    void InterfaceForces(bool output = false) override;

    //! Print interfaces
    void Print(std::ostream& os) const override;

    //! Print summary of active set status to screen
    void PrintActiveSet() const override;

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    void PostprocessQuantitiesPerInterface(Teuchos::RCP<Teuchos::ParameterList> outputParams) final;

    //!@}

    //! @name Debugging methods
    //!@{

    /*! \brief Visualize contact stuff with gmsh

     \param step (in): current time step index
     \param iter (in): current iteration index
     */
    void VisualizeGmsh(const int step, const int iter) override;

    //!@}

    /*! @name Purely virtual functions
     *
     * All these functions are defined in one or more specific derived classes,
     * i.e CONTACT::LagrangeStrategy or CONTACT::PenaltyStrategy.
     * As the base class MORTAR::StrategyBase is always called from the control routine
     * (time integrator), these functions need to be defined purely virtual here.
     */
    //!@{

    bool ActiveSetSemiSmoothConverged() const override = 0;
    bool ActiveSetConverged() override = 0;
    virtual int ActiveSetSteps() = 0;
    virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const = 0;
    virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const = 0;
    double ConstraintNorm() const override = 0;
    virtual void EvaluateContact(
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual void EvaluateFriction(
        Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    void EvaluateRelMovPredict() override = 0;
    double InitialPenalty() override = 0;
    void Initialize() override = 0;
    void InitializeUzawa(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff) override = 0;
    void Recover(Teuchos::RCP<Epetra_Vector> disi) override = 0;
    void ResetActiveSet() override = 0;
    void ResetPenalty() override = 0;
    void ModifyPenalty() override = 0;
    void BuildSaddlePointSystem(Teuchos::RCP<CORE::LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override = 0;
    void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override = 0;
    virtual void EvalConstrRHS() = 0;
    void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override = 0;
    void UpdateActiveSet() override = 0;
    void UpdateActiveSetSemiSmooth(const bool firstStepPredictor = false) override = 0;
    void UpdateUzawaAugmentedLagrange() override = 0;
    void UpdateConstraintNorm(int uzawaiter = 0) override = 0;

    //!@}

    /*! @name Empty functions (meshtying)
     *
     * All these functions only have functionality in meshtying simulations, thus they
     * are defined as empty here in the case of contact. They can be called from the
     * control routine (time integrator), whenever you like.
     */
    //!@{

    void RedistributeMeshtying() final {}
    void RestrictMeshtyingZone() override {}
    void EvaluateMeshtying(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) override
    {
    }
    Teuchos::RCP<const Epetra_Vector> MeshInitialization() override { return Teuchos::null; };

    void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis) override {}

    //!@}

   protected:
    //! @name Pre/Postoperators
    //!@{

    //! Run after the StoreDirichletStatus() routine has been called
    virtual void PostStoreDirichletStatus(Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmaps){};

    /*! \brief Run at the beginning of the Evaluate() routine
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void PreEvaluate(CONTACT::ParamsInterface& cparams){};

    /*! \brief Run in the end of the Evaluate() routine
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void PostEvaluate(CONTACT::ParamsInterface& cparams){};

    /*! \brief Run in the end of the Setup() routine
     *
     *  Can be used to redistribute member variables of derived classes, if necessary.
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void PostSetup(bool redistributed, bool init){};

    //!@}

    //! @Internal evaluate routines
    //!@{

    /*! \brief Compute force and stiffness terms
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time
     * integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void EvalForceStiff(CONTACT::ParamsInterface& cparams);

    /*! \brief Compute force terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration
     *
     *  \author hiermeier \date 03/2016 */
    virtual void EvalForce(CONTACT::ParamsInterface& cparams);

    /*! \brief Compute the constraint rhs
     *
     *  \param(in) cparams: parameter interface between the contact objects and
     *                      the structural time integrator
     *
     *  \author hiermeier \date 12/17 */
    virtual void EvalStaticConstraintRHS(CONTACT::ParamsInterface& cparams);

    /** \brief Run at the very beginning of a call to STR::ModelEvaluator::Evalute*
     *
     *  \param cparams (in): parameter interface between the contact objects and
     *                       the structural time integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void RunPreEvaluate(CONTACT::ParamsInterface& cparams);

    /** \brief Run in the end of a call to STR::ModelEvaluator::EvaluteForce/Stiff/ForceStiff
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     *                       integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void RunPostEvaluate(CONTACT::ParamsInterface& cparams);

    /*! \brief Recover the current state
     *
     *  The main task of this method is to recover the Lagrange multiplier solution.
     *  The Lagrange multiplier solution will be stored inside the corresponding strategy
     *  and is necessary for different internal evaluation methods. If the Lagrange multiplier
     *  is condensed, this method is the right place to recover it from the displacement solution.
     *  If it is not condensed (saddle-point system) use the ResetLagrangeMultiplier routine
     *  instead.
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     *                       integration
     *  \param xold    (in): old solution vector of the NOX solver
     *  \param dir     (in): current search direction (in general NOT the actual step, keep in mind
     *                       that the step length can differ from 1.0)
     *  \param xnew    (in): new solution vector of the NOX solver
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual void RunPostComputeX(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
        const Epetra_Vector& dir, const Epetra_Vector& xnew);

    /*! \brief run pre-compute x routine for contact
     *
     *  This method is called at the very beginning of the NOX::NLN::Group::ComputeX()
     *  routine and gives you the opportunity to modify/augment the current Newton
     *  direction.
     *
     *  \param cparams (in)    : parameter interface between the contact objects
     *                           and the structural time integration
     *  \param xold    (in)    : old solution vector of the NOX solver
     *  \param dir     (in/out): current search direction (in general NOT the actual
     *                           step, keep in mind that the step length can differ from 1.0)
     *
     *  \date 03/2017
     *  \author hiermeier */
    virtual void RunPreComputeX(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable);

    /*! \brief Executed at the end of the NOX::NLN::Group::applyJacobianInverse()
     *  method
     *
     *  \param cparams: parameter interface between the contact objects and the
     *                  structural time integration
     *  \param rhs    : read-only access to the rhs vector
     *  \param result : full access to the result vector
     *  \param xold   : read-only access to the jacobian
     *  \param grp    : read only access to the group object
     *
     *  \author hiermeier \date 12/2017 */
    virtual void RunPostApplyJacobianInverse(const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& rhs, Epetra_Vector& result, const Epetra_Vector& xold,
        const NOX::NLN::Group& grp);

    /*! \brief run pre-compute x routine for contact
     *
     *  This routine is called in the end of a ::NOX::Solver::step() call.
     *
     *  \param cparams (in)    : parameter interface between the contact objects
     *                           and the structural time integration
     *
     *  \author hiermeier \date 03/2017  */
    virtual void RunPostIterate(const CONTACT::ParamsInterface& cparams);

    /// run before before the nonlinear solver starts
    virtual void RunPreSolve(const Teuchos::RCP<const Epetra_Vector>& curr_disp,
        const CONTACT::ParamsInterface& cparams);

    /*! \brief Reset the internal stored Lagrange multipliers
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     *                       integration
     *  \param xnew    (in): new solution vector of the NOX solver
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual void ResetLagrangeMultipliers(
        const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew);

    //! Evaluate the weighted gap gradient error
    virtual void EvalWeightedGapGradientError(CONTACT::ParamsInterface& cparams);

    virtual void CorrectParameters(
        CONTACT::ParamsInterface& cparams, const NOX::NLN::CorrectionType type);

    /*! \brief Remove condensed contact contributions from the structural right-hand side
     *
     *  \param(in) str_rhs: reference to the structural right-hand side
     *  \author hiermeier \date 03/18 */
    virtual void RemoveCondensedContributionsFromRhs(Epetra_Vector& str_rhs) const;

    //!@}

   protected:
    //! access the contact interfaces of the concrete strategies (read and write)
    virtual std::vector<Teuchos::RCP<CONTACT::Interface>>& Interfaces() = 0;

    //! access the contact interfaces of the concrete strategies (read-only)
    virtual const std::vector<Teuchos::RCP<CONTACT::Interface>>& Interfaces() const = 0;

    /*! \brief Evaluate contact

     This is just a tiny control routine, deciding which Evaluate-routine
     of those listed below is to be called (based on input-file information).
     Note that into ALL derived Evaluate() routines, a REFERENCE to the pointer
     on the effective stiffness matrix is handed in. This way, after building the
     new effective stiffness matrix with contact, we can simply let the pointer
     kteff point onto the new object. The same is true for the effective force
     vector feff. Be careful: kteff is of type Teuchos::RCP<CORE::LINALG::SparseOperator>&.

     \param kteff (in/out): effective stiffness matrix (without -> with contact)
     \param feff (in/out): effective residual / force vector (without -> with contact)

     */
    void Evaluate(Teuchos::RCP<CORE::LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) override;

    /*! \brief Evaluate relative movement of contact bodies

     This is for evaluating the relative movement of contact bodies. This
     can either be done with regarding the different movement of material points
     or regarding the change of mortar projection. The second possibility
     is definitely objective whereas the first possibility is objective
     only when the gap is zero. */
    void EvaluateRelMov() override;

    /*! \brief Initialize and evaluate Mortar stuff for the next Newton step

     This method first checks if we are dealing with self contact and updates
     the interface slave and master sets if so. Then it resets the global
     Mortar matrices D and M and the global gap vector g accordingly.

     The nodal quantities computed in InitEvalInterface() are then assembled
     to global matrices and vectors respectively. No setup of the global system
     is to be done here yet, so there is no need to pass in the effective
     stiffness K or the effective load vector f. */
    void InitMortar() override;
    void AssembleMortar() override;

    /*! \brief Initialize and evaluate interface for the next Newton step

     This method calls Initialize() on all contact interfaces, which
     resets all kind of nodal quantities like normal vector, weighted
     gap or Mortar and linearization maps. It then calls Evaluate() on
     all contact interfaces, which does all the geometric contact stuff.
     Concretely, this is an evaluation of all involved quantities at nodal
     level plus the setup of all corresponding linearizations.
     It includes the nodal normal calculations, contact search, projection
     and overlap detection, integration of the  Mortar terms D, M and of the
     weighted gap. Additionally, the linearizations of geometric quantities
     (delta_n, delta_t, delta_D, delta_M) are calculated. */
    void InitEvalInterface() override { InitEvalInterface(Teuchos::null); };
    virtual void InitEvalInterface(Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    /*! check the parallel distribution and initialize a possible
     *  redistribution */
    void CheckParallelDistribution(const double& t_start);

    /// update the parallel distribution status
    void UpdateParallelDistributionStatus(const double& my_total_time);

    /*! \brief Update Mortar matrices D and M

     The std::string input parameter defines in which direction the conversion
     is to be performed. Currently only the possibilities "old" and "current"
     exist, with "old" meaning the Mortar matrices of the last time / load step
     will be set to the current values D_n+1 and M_n+1 (this happens after
     completion of a time / load step!). The std::string "current" addresses the
     current Mortar matrices, which when called will be reset to the last
     converged values D_n and M_n (this happens in the predictor step when
     the active set has not yet converged!).

     \param state (in): std::string defining in which direction to convert D and M
     */
    void StoreDM(const std::string& state);

    /*! \brief Store current (contact) nodal entries to old ones

     Contact nodes own their current entries and old ones (last converged
     state) from. p.e. the mortar matrices D and M. This function writes the
     current ones to the old ones. */
    void StoreToOld(MORTAR::StrategyBase::QuantityType type);

    /*! \brief Update global self contact state

     This becomes necessary for self contact simulations, because in a
     self contact interface master and slave status are assigned dynamically
     and thus the global maps change constantly.

     */
    void UpdateGlobalSelfContactState();

    /// access global self contact lagrange multiplier map (read & write)
    inline const Epetra_Map& GSelfContactLmMap() const
    {
      return *Data().GSelfContactLmDofRowMapPtr();
    }

    inline const Epetra_Map& GSelfContactRefMap() const
    {
      return *Data().GSelfContactRefDofRowMapPtr();
    }

   private:
    /*!
    \brief Check if this is the first time step of the simulation

    As we don't have the time step counter available here, let's check for the size of some member
    variables: a size of zero indicates the first time step.

    \warning This checks relies on the proper (re-)initialization of some member variables. Behavior
    could change, if these member variables are (re-)initialized differently.

    \return Boolean flag to indicate the first time step (true) or not (false)
    */
    bool IsFirstTimeStep() const;

    //! @name Parallel redistribution and ghosting
    //! @{

    /*!
    \brief Decide whether interface discretizations need to be rebalanced

    The decision to perform rebalancing is based on user input as well as history of
    - the max-to-min ratio of contact evaluation time across all processes
    - the max-to-min ratio of the number of row slave elements across all processes

    averaged over all contact evaluations of the previous time step.

    Naturally, serial runs do never require rebalancing.

    \sa CheckParallelDistribution(), UpdateParallelDistributionStatus()

    @param[in] Flag to indicate first time step after start/restart of simulation
    @return True if rebalancing is necessary, false otherwise.
    */
    bool IsRebalancingNecessary(const bool first_time_step);

    /*!
    \brief Compute and reset indicators for necessity of parallel rebalancing

    Unbalance is measured as the max-to-min ratio of eveluation time / number of row slave elements
    over all processes.

    We average the unbalance of interface evaluation time and interface element count over all
    contact evaluations of a time step. These will be used to decide, whether rebalancing of the
    interface discretization is necessary. At the end, reset the indicators in preparation for the
    next time step.

    \sa CheckParallelDistribution(), IsRebalancingNecessary(), UpdateParallelDistributionStatus()

    @param[in/out] time_average Average max-to-min ratio of evlation time accross procs over all
                                evaluations of previous time step
    @param[in/out] elements_average Average max-to-min ratio of row slave elements accross procs
                                    over all evaluations of previous time step
    @param[in] Flag to indicate first time step after start/restart of simulation
    */
    void ComputeAndResetParallelBalanceIndicators(double& time_average, double& elements_average);

    /*!
    \brief Print indicators for current status of parallel load balancing

    Indicators will be printed to screen on proc 0.

    @param[in/out] time_average Average max-to-min ratio of evlation time accross procs over all
                                evaluations of previous time step
    @param[in/out] elements_average Average max-to-min ratio of row slave elements accross procs
                                    over all evaluations of previous time step
    @param[in] max_time_unbalance Upper bound for imbalance in evaluation time given in input file
    */
    void PrintParallelBalanceIndicators(
        double& time_average, double& elements_average, const double& max_time_unbalance) const;

    /*!
    \brief Is an update of the interface ghosting necessary?

    It depends on the actual interface ghosting strategy, if the interface ghosting needs to be
    update in this time step:
    - Any redundant storage does not require an update. It just has to be done in the very first
    time step.
    - In case of any non-redundant storage, the interface ghosting needs to be updated in every time
    step to guarantee sufficient ghosting for a correct and safe contact search.

    @param ghosting_strategy User-chosen strategy to extend the interface ghosting
    @param[in] Flag to indicate first time step after start/restart of simulation

    @return Flag to indicate, whether ghosting needs to be updated (true) or not (false)
    */
    bool IsUpdateOfGhostingNecessary(
        const INPAR::MORTAR::ExtendGhosting& ghosting_strategy, const bool first_time_step) const;

    /*!
    \brief Calculate absolute value of mean velocity of interface for binning

    For each interface, extract its velocity DOFs and compute their mean value.
    Then, take the absolute value and store it to #ivel_.

    \param[in] velocity Vector with velocities for all solid DOFs
    */
    void CalcMeanVelocityForBinning(const Epetra_Vector& velocity);

    /*!
    \brief Update parallel load balancing of each contact interface and guarantee correct ghosting

    We hand in the current global displacement state, so that a contact search can be performed and
    we can SetState the displacement field.

    The current velocity state is required in case of extending the ghosting via binning to account
    for relative motion between interfaces.

    \post Each contact interface is FillComplete().

    @param[in] displacement
    @param[in] velocity
    @return TRUE if the interface has been redistributed. Return FALSE otherwise.
    */
    bool RedistributeWithSafeGhosting(
        const Epetra_Vector& displacement, const Epetra_Vector& velocity);

    /*!
    \brief Redistribute all contact interfaces in parallel (legacy implementation)

    We hand in the current global displacement state, so that a contact search can be performed and
    we can SetState the displacement field.

    The current velocity state is required in case of extending the ghosting via binning to account
    for relative motion between interfaces.

    \post Each contact interface is FillComplete().

    \warning The interplay of parallel redistribution and extension of the interface ghosting is
    somehow fragile. The interface ghosting is only updated after an actual redistribution. However,
    it can happen that redistribution is not necessary (or disabled by the user), but the ghosting
    still needs to be updated due to changes in the contact area topology (e.g. large sliding). Such
    cases are not captured properly. As a result, the ghosting does not include all necessary master
    nodes and, thus, the contact search fails to detect all close slave/master element pairs. Use
    RedistributeWithSafeGhosting() instead!

    \param[in] dis Current displacement state
    \param[in] vel Current velocity state

    \return TRUE if the interface has been redistributed. Return FALSE otherwise.
    */
    bool RedistributeContactOld(
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel);

    //! @}

    /*! \brief Create the global Lagrange multiplier DoF row map
     *
     *  The global Lagrange multiplier DoF row map is created in a deterministic
     *  manner based on the previously created global slave DoF row map. This is
     *  necessary for the later ReplaceMap calls. Especially, the std::sort during
     *  a CORE::LINALG::MergeMap call would otherwise destroy the correlation. This becomes
     *  obvious if more than one interface is considered.
     *
     *  \pre The method UpdateLagMultSets() has to be called on each involved
     *  interface before this method is executed.
     *
     *  \param[in] gsdofrowmap: Already new global slave DoF row map.
     *
     *  \return New Lagrange multiplier DoF row map in correlation to the given
     *          global slave DoF row map.
     *
     *  \author hiermeier \date 10/17 */
    Teuchos::RCP<Epetra_Map> CreateDeterministicLMDofRowMap(const Epetra_Map& gsdofrowmap) const;

    /*! return the mutable contact abstract data container
     *
     * \remark This has to stay PRIVATE, otherwise the function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    CONTACT::AbstractStratDataContainer& Data()
    {
      if (data_ptr_.is_null()) dserror("The AbstractStratDataContainer is not initialized!");
      return *data_ptr_;
    };

    /*! return the read-only abstract contact data container
     *
     * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
     *
     * \author hiermeier
     * \date 05/16 */
    const CONTACT::AbstractStratDataContainer& Data() const
    {
      if (data_ptr_.is_null()) dserror("The AbstractStratDataContainer is not initialized!");
      return *data_ptr_;
    };

   protected:
    //! Global Lagrange multiplier dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& glmdofrowmap_;

    //! Global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gsnoderowmap_;

    //! Global master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gmnoderowmap_;

    //! Global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gsdofrowmap_;

    //! Global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gmdofrowmap_;

    //! Global internal dof row map
    Teuchos::RCP<Epetra_Map>& gndofrowmap_;

    //! Global slave and master dof row map (salve+master map)
    Teuchos::RCP<Epetra_Map>& gsmdofrowmap_;

    //! Global displacement dof row map (s+m+n map)
    Teuchos::RCP<Epetra_Map>& gdisprowmap_;

    //! @name Active set and slip set
    //!@{

    //! Global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gactivenodes_;

    //! Global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gactivedofs_;

    //! Global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& ginactivenodes_;

    //! Global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& ginactivedofs_;

    /*! \brief Global dof row map of matrix \f$N\f$ (of all interfaces)
     *
     * \todo What is the matrix N?
     */
    Teuchos::RCP<Epetra_Map>& gactiven_;

    /*! \brief Global dof row map of matrix \f$T\f$ (of all interfaces)
     *
     * \todo What is the matrix T?
     */
    Teuchos::RCP<Epetra_Map>& gactivet_;

    //! Global slip slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gslipnodes_;

    //! Global slip slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gslipdofs_;

    /*! \brief Global row map of matrix \f$T\f$ for slip dofs (of all interfaces)
     *
     * \todo What is the matrix T?
     */
    Teuchos::RCP<Epetra_Map>& gslipt_;

    //!@}

    //! Global slave row map of vertex nodes
    Teuchos::RCP<Epetra_Map>& gsdofVertex_;

    //! Global slave row map of edge nodes
    Teuchos::RCP<Epetra_Map>& gsdofEdge_;

    //! Global slave row map of surface nodes
    Teuchos::RCP<Epetra_Map>& gsdofSurf_;

    //! @name Parallel redistribution and ghosting
    //!@{

    //! Parallel unbalance factors (evaluation time) for current time step \f$t_{n+1}\f$
    std::vector<double>& unbalanceEvaluationTime_;

    //! Parallel unbalance factors (num. of slave elements) for current time step \f$t_{n+1}\f$
    std::vector<int>& unbalanceNumSlaveElements_;

    //! Global Lagrange multiplier dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pglmdofrowmap_;

    //! Global slave dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pgsdofrowmap_;

    //! Global master dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pgmdofrowmap_;

    //! Global slave and master dof row map before parallel redistribution
    Teuchos::RCP<Epetra_Map>& pgsmdofrowmap_;

    //!< Global dirichlet toggle of all slave dofs before parallel redistribution
    Teuchos::RCP<Epetra_Vector>& pgsdirichtoggle_;

    //!@}

    //! @name Binning strategy
    //!@{

    //!< Initial element columns map for binning strategy (slave and master)
    std::vector<Teuchos::RCP<Epetra_Map>>& initial_elecolmap_;

    //!@}

    //! Global Mortar matrix \f$D\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& dmatrix_;

    //! Global Mortar matrix \f$M\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& mmatrix_;

    //! Global weighted gap vector \f$g\f$
    Teuchos::RCP<Epetra_Vector>& g_;

    //! Global tangential right-hand side vector (formulation with incremental #z_)
    Teuchos::RCP<Epetra_Vector>& tangrhs_;

    /*! \brief Global inactive right-hand side vector
     *
     * This is used for the formulation with incremental #z_ and saddle point system.
     */
    Teuchos::RCP<Epetra_Vector>& inactiverhs_;

    //! Global structural contact contributions to right-hand side vector at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& strcontactrhs_;

    //! Global constraint right-hand side vector (only for saddlepoint problems)
    Teuchos::RCP<Epetra_Vector>& constrrhs_;

    /*! \brief Global Matrix LinD containing slave fc derivatives
     *
     * \todo What is fc?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& lindmatrix_;

    /*! \brief Global Matrix LinM containing master fc derivatives
     *
     * \todo What is fc?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& linmmatrix_;

    Teuchos::RCP<CORE::LINALG::SparseMatrix>& kteffnew_;

    //! Global Mortar matrix \f$D\f$ at end of last time step \f$t_{n}\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& dold_;

    //! Global Mortar matrix \f$M\f$ at end of last time step \f$t_{n}\f$
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& mold_;

    //!< Current vector of Lagrange multipliers at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& z_;

    //! Old vector of Lagrange multipliers at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector>& zold_;

    /*! \brief Lagrange multiplier vector increment within SaddlePointSolve
     *
     * \note This is \em not the increment of #z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!
     */
    Teuchos::RCP<Epetra_Vector>& zincr_;

    //! Vector of Lagrange multipliers from last Uzawa step
    Teuchos::RCP<Epetra_Vector>& zuzawa_;

    /*! \brief Vector of normal contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #forcenormal_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& stressnormal_;

    /*! \brief Vector of tangential contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #forcetangential_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& stresstangential_;

    /*! \brief Vector of normal contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #stressnormal_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& forcenormal_;

    /*! \brief Vector of tangential contact forces at \f$t_{n+1}\f$
     *
     * \todo What's the difference to #stresstangential_? Update documentation!
     */
    Teuchos::RCP<Epetra_Vector>& forcetangential_;

    //! @name Counters and indices
    //!@{

    //! Time step index
    int& step_;

    //! Nonlinear iteration index, e.g. Newton iteration
    int& iter_;

    //!@}

    //! @name Status flags
    //!@{

    //! Flag indicating global contact status
    bool& isincontact_;

    //! Flag indicating global contact status of this time step (history)
    bool& wasincontact_;

    //! Flag indicating global contact status of last time step
    bool& wasincontactlts_;

    //! Flag indicating potential self contact
    bool& isselfcontact_;

    //! Flag for frictional contact
    bool& friction_;

    //! Flag for non-smooth contact algorithm
    bool& nonSmoothContact_;

    //! Flag for regularized contact
    bool& regularized_;

    /*! \brief Flag indicating whether transformation should be applied
     *
     * \todo Which transformation?
     */
    bool& dualquadslavetrafo_;

    //!@}

    /*! \brief Transformation matrix \f$T\f$ for dual quad 3D case
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& trafo_;

    /*! \brief Transformation matrix \f$T\f$ for dual quad 3D case (all problem dofs)
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> systrafo_;

    /*! \brief Inverse transformation matrix \f$T\f$ for dual quad 3D case (all problem dofs)
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix> invsystrafo_;

    /*! \brief Inverse transformation matrix \f$T^{-1}\f$ for dual quad 3D case
     *
     * \todo What does quad refer to? Quadratic or quadrilateral elements?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& invtrafo_;

    /*! \brief Modified global Mortar matrix \f$D\d$
     *
     * \todo What modifications?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& dmatrixmod_;

    /*! \brief Modified global Mortar matrix Dold
     *
     * \todo What modifications?
     */
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& doldmod_;

    /*! \brief Integration time
     *
     * \todo Is this the wall clock time required to perform the mortar integration?
     */
    double& inttime_;

    //! Mean velocity of each interface
    std::vector<double>& ivel_;

    //! Current used solving strategy
    INPAR::CONTACT::SolvingStrategy& stype_;

    //! Direction in which the contact constraints are formulated
    INPAR::CONTACT::ConstraintDirection& constr_direction_;

   private:
    /*!
    \brief Copy constructor

    @param old Instance of this class to be copied
    */
    AbstractStrategy(const AbstractStrategy& old) = delete;

    //! pointer to the data container object
    Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr_;

    //! pointer to the NOX::NLN::CONSTRAINT::Interface::Required object
    Teuchos::RCP<CONTACT::NoxInterface> noxinterface_ptr_;

  };  // namespace CONTACT
}  // namespace CONTACT

//! << operator
std::ostream& operator<<(std::ostream& os, const CONTACT::AbstractStrategy& strategy);

BACI_NAMESPACE_CLOSE

#endif
