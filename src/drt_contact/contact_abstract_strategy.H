/*---------------------------------------------------------------------*/
/*!
\file contact_abstract_strategy.H

\brief Main abstract class for contact solution strategies

\level 2

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_STRATEGY_H
#define CONTACT_STRATEGY_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>

#include <Epetra_Operator.h>

#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_utils_matrix_vector_enums.H"

// forward declarations
namespace NOX
{
  namespace NLN
  {
    class Group;
  }  // namespace NLN
}  // namespace NOX

namespace DRT
{
  class Discretization;
}  // namespace DRT

namespace LINALG
{
  class MultiMapExtractor;
  class SparseMatrix;
}  // namespace LINALG

namespace CONTACT
{
  // forward declarations
  class CoInterface;
  class ParamsInterface;
  class NoxInterface;
  // class WearInterface;

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

    //! destructor
    virtual ~AbstractStratDataContainer(){};

    //! @name Accessors
    //! @{

    //! return parallel unbalance factors (time) for current time step \f$t_{n+1}\f$
    std::vector<double>& UnbalanceTimeFactors() { return tunbalance_; };
    const std::vector<double>& UnbalanceTimeFactors() const { return tunbalance_; };

    //! return parallel unbalance factors (elements) for current time step \f$t_{n+1}\f$
    std::vector<int>& UnbalanceElementFactors() { return eunbalance_; };
    const std::vector<int>& UnbalanceElementFactors() const { return eunbalance_; };

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
    Teuchos::RCP<LINALG::SparseMatrix>& DMatrixPtr() { return dmatrix_; };
    Teuchos::RCP<const LINALG::SparseMatrix> DMatrixPtr() const { return dmatrix_; };
    LINALG::SparseMatrix& DMatrix()
    {
      if (dmatrix_.is_null()) dserror("The dmatrix_ is not initialized!");
      return *dmatrix_;
    }

    //! return global Mortar matrix M
    Teuchos::RCP<LINALG::SparseMatrix>& MMatrixPtr() { return mmatrix_; };
    Teuchos::RCP<const LINALG::SparseMatrix> MMatrixPtr() const { return mmatrix_; };
    LINALG::SparseMatrix& MMatrix()
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
    Teuchos::RCP<LINALG::SparseMatrix>& DLinMatrixPtr() { return lindmatrix_; };
    Teuchos::RCP<const LINALG::SparseMatrix> DLinMatrixPtr() const { return lindmatrix_; };
    LINALG::SparseMatrix& DLinMatrix()
    {
      if (lindmatrix_.is_null()) dserror("The augDnLinMatrixPtr_ is not initialized!");
      return *lindmatrix_;
    }

    //! return global Matrix LinM containing master fc derivatives
    Teuchos::RCP<LINALG::SparseMatrix>& MLinMatrixPtr() { return linmmatrix_; };
    Teuchos::RCP<const LINALG::SparseMatrix> MLinMatrixPtr() const { return linmmatrix_; };
    LINALG::SparseMatrix& MLinMatrix()
    {
      if (linmmatrix_.is_null()) dserror("The augMnLinMatrixPtr_ is not initialized!");
      return *linmmatrix_;
    }

    //! return global Matrix kteffnew containing modified jacobian
    Teuchos::RCP<LINALG::SparseMatrix>& kteffnewMatrixPtr() { return kteffnew_; };
    Teuchos::RCP<const LINALG::SparseMatrix> kteffnewMatrixPtr() const { return kteffnew_; };
    LINALG::SparseMatrix& kteffnewMatrix()
    {
      if (kteffnew_.is_null()) dserror("The kteffnewMatrixPtr is not initialized!");
      return *kteffnew_;
    }

    //! return global Mortar matrix D (last end-point \f$t_{n}\f$)
    Teuchos::RCP<LINALG::SparseMatrix>& OldDMatrixPtr() { return dold_; };
    Teuchos::RCP<const LINALG::SparseMatrix> OldDMatrixPtr() const { return dold_; };

    //! return global Mortar matrix M (last end-point \f$t_{n}\f$)
    Teuchos::RCP<LINALG::SparseMatrix>& OldMMatrixPtr() { return mold_; };
    Teuchos::RCP<const LINALG::SparseMatrix> OldMMatrixPtr() const { return mold_; };

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
    Teuchos::RCP<LINALG::SparseMatrix>& TrafoPtr() { return trafo_; };
    Teuchos::RCP<const LINALG::SparseMatrix> TrafoPtr() const { return trafo_; };

    //! return inverse trafo matrix T^(-1) for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix>& InvTrafoPtr() { return invtrafo_; };
    Teuchos::RCP<const LINALG::SparseMatrix> InvTrafoPtr() const { return invtrafo_; };

    //! return modified global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix>& ModifiedDMatrixPtr() { return dmatrixmod_; };
    Teuchos::RCP<const LINALG::SparseMatrix> ModifiedDMatrixPtr() const { return dmatrixmod_; };

    //! return modified global Mortar matrix Dold
    Teuchos::RCP<LINALG::SparseMatrix>& OldModifiedDMatrixPtr() { return doldmod_; };
    Teuchos::RCP<const LINALG::SparseMatrix> OldModifiedDMatrixPtr() const { return doldmod_; };

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

    INPAR::MORTAR::ParRedist& ParType() { return partype_; };
    INPAR::MORTAR::ParRedist ParType() const { return partype_; };
    //! @}

   private:
    //! parallel unbalance factors (time) for current time step \f$t_{n+1}\f$
    std::vector<double> tunbalance_;

    //! parallel unbalance factors (eles) for current time step \f$t_{n+1}\f$
    std::vector<int> eunbalance_;

    //! global Lagrange mult. dof row map (of all interfaces)
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

    //! global slave dof row map of vertex nodes
    Teuchos::RCP<Epetra_Map> gsdofVertex_;

    //! global slave dof row map of edge nodes
    Teuchos::RCP<Epetra_Map> gsdofEdge_;

    //! global slave dof row map of surface nodes
    Teuchos::RCP<Epetra_Map> gsdofSurf_;

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

    //! initial col ele map for binning strategy (s m)
    std::vector<Teuchos::RCP<Epetra_Map>> initial_elecolmap_;

    //! global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix> dmatrix_;

    //! global Mortar matrix M
    Teuchos::RCP<LINALG::SparseMatrix> mmatrix_;

    //! global weighted gap vector g
    Teuchos::RCP<Epetra_Vector> g_;

    //! global tangential rhs vector (formulation with incremental z_)
    Teuchos::RCP<Epetra_Vector> tangrhs_;

    //! gloabl inactive rhs vector (formulation with incremental z_ and saddle point system)
    Teuchos::RCP<Epetra_Vector> inactiverhs_;

    //! structural contact right-hand-side vector at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> strContactRhsPtr_;

    //! global constraint rhs vector (only for saddlepoint problems)
    Teuchos::RCP<Epetra_Vector> constrrhs_;

    //! global Matrix LinD containing slave fc derivatives
    Teuchos::RCP<LINALG::SparseMatrix> lindmatrix_;

    //! global Matrix LinM containing master fc derivatives
    Teuchos::RCP<LINALG::SparseMatrix> linmmatrix_;

    //! global K matrix
    Teuchos::RCP<LINALG::SparseMatrix> kteffnew_;

    //! global Mortar matrix D (last end-point \f$t_{n}\f$)
    Teuchos::RCP<LINALG::SparseMatrix> dold_;

    //! global Mortar matrix M (last end-point \f$t_{n}\f$)
    Teuchos::RCP<LINALG::SparseMatrix> mold_;

    //! current vector of Lagrange multipliers at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> z_;

    //! old vector of Lagrange multipliers at \f$t_{n}\f$
    Teuchos::RCP<Epetra_Vector> zold_;

    /*! \brief Lagrange multiplier vector increment within SaddlePointSolve
     *
     *  \remark This is NOT the increment of z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!) */
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

    //! time step index at \f$t_{n+1}\f$
    int stepnp_;

    //! nonlinear (Newton) iteration index
    int iter_;

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

    //! transformation matrix T for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix> trafo_;

    //! inverse trafo matrix T^(-1) for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix> invtrafo_;

    //! modified global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix> dmatrixmod_;

    //! modified global Mortar matrix Dold
    Teuchos::RCP<LINALG::SparseMatrix> doldmod_;

    //! integration time
    double inttime_;

    //! mean interface velocity
    std::vector<double> ivel_;

    //! current used solving strategy
    INPAR::CONTACT::SolvingStrategy stype_;

    //! direction in which the contact constraints are formulated
    INPAR::CONTACT::ConstraintDirection constr_direction_;

    //! parallel redistribution type
    INPAR::MORTAR::ParRedist partype_;

  };  // class AbstractStratDataContainer


  /*! \brief Main abstract class for contact solution strategies
   *
   *  This is the templating abstract class for all contact solution algorithms.
   *  Every solution algorithm has to fit into the set of functions and calls defined herein
   *  and has to be specified in a corresponding subclass defining the concrete algorithmic steps.
   *
   *  This class it itself derived from the MORTAR::StrategyBase class, which is an even
   *  more abstract framework for any solution strategies involving mortar coupling.
   *
   *  Refer also to the Semesterarbeit of Bernd Budich, 2009
   *
   *  \author popp (popp@lnm.mw.tum.de) */
  class CoAbstractStrategy : public MORTAR::StrategyBase
  {
   public:
    /*! Standard constructor
     *
     *  This constructor uses the given DataContainer to store and share all its
     *  member variables. The declared member variables are just references to
     *  the container content! */
    CoAbstractStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& stratData,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params, int dim, const Teuchos::RCP<const Epetra_Comm>& comm,
        double alphaf, int maxdof);

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

    //! Destructor
    virtual ~CoAbstractStrategy(){};

    //! return the current solution type
    virtual INPAR::CONTACT::SolvingStrategy Type() const { return stype_; }

    //! @name Access methods
    //! @{
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
      dserror("ERROR: Map not available in abstract strategy!");
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
     *  \param bt (in): Desired vector block type, e.g. block_displ, block_constraint, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum DRT::UTILS::VecBlockType& bt) const
    {
      dserror("Not yet implemented!");
      exit(EXIT_FAILURE);
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
     *  \param bt (in): Desired vector block type, e.g. block_displ, block_constraint, ...
     *
     *  \author hiermeier \date 08/17  */
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtrForNormCheck(
        const enum DRT::UTILS::VecBlockType& bt) const
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
    };

    /*! \brief Return the desired matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     *
     *  \param bt (in): Desired matrix block type, e.g. block_displ_displ, block_displ_lm, ...
     *
     *  \date 05/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = NULL) const
    {
      dserror("Not yet implemented!");
      exit(EXIT_FAILURE);
    };

    //! Apply modifications (e.g. condensation) directly before linear solve
    virtual void RunPreApplyJacobianInverse(
        Teuchos::RCP<LINALG::SparseMatrix> kteff, Epetra_Vector& rhs)
    { /* do nothing */
      ;
    }

    /*! Return the condensed matrix block pointer (read-only)
     *
     *  \remark Please note, that a Teuchos::null pointer is returned, if no active contact
     *  contributions are present.
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetCondensedMatrixBlockPtr(
        Teuchos::RCP<LINALG::SparseMatrix>& kteff, const double& timefac_np) const
    {
      dserror("Not yet implemented!");
      exit(EXIT_FAILURE);
    };

    //! Return global slave node row map
    const Teuchos::RCP<Epetra_Map> SlaveRowNodes() { return Data().GSlNodeRowMapPtr(); }
    Teuchos::RCP<const Epetra_Map> SlRowNodesPtr() const { return Data().GSlNodeRowMapPtr(); }
    const Epetra_Map& SlRowNodes() const { return *Data().GSlNodeRowMapPtr(); }

    //! Return global slave node row map
    Teuchos::RCP<const Epetra_Map> MaRowNodesPtr() const { return Data().GMaNodeRowMapPtr(); }
    const Epetra_Map& MaRowNodes() const { return *Data().GMaNodeRowMapPtr(); }

    //! Return global active node row map
    virtual const Teuchos::RCP<Epetra_Map> ActiveRowNodes()
    {
      return Data().GActiveNodeRowMapPtr();
    };
    virtual Teuchos::RCP<const Epetra_Map> ActiveRowNodes() const
    {
      return Data().GActiveNodeRowMapPtr();
    };

    //! Return global slip node row map
    const Teuchos::RCP<Epetra_Map> SlipRowNodes() { return Data().GSlipNodeRowMapPtr(); };
    Teuchos::RCP<const Epetra_Map> SlipRowNodes() const { return Data().GSlipNodeRowMapPtr(); };

    //! Return global slave dof row map
    const Teuchos::RCP<Epetra_Map> SlaveRowDofs() { return Data().GSlDofRowMapPtr(); }

    //! Return global active dof row map
    virtual const Teuchos::RCP<Epetra_Map> ActiveRowDofs() { return Data().GActiveDofRowMapPtr(); }

    //! Return global master dof row map
    const Teuchos::RCP<Epetra_Map> MasterRowDofs() { return Data().GMaDofRowMapPtr(); }

    //! Return global slave dof row map
    const Teuchos::RCP<Epetra_Map> SlaveMasterRowDofs() { return Data().GSlMaDofRowMapPtr(); }

    //! Return not redistributed global slave dof row map
    const Teuchos::RCP<Epetra_Map> NotReDistSlaveRowDofs() { return Data().PGSlDofRowMapPtr(); }

    //! Return not redistributed global master dof row map
    const Teuchos::RCP<Epetra_Map> NotReDistMasterRowDofs() { return Data().PGMaDofRowMapPtr(); }

    //! provides maps needed for contact/meshtying specific multigrid preconditioners
    void CollectMapsForPreconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap,
        Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap,
        Teuchos::RCP<Epetra_Map>& ActiveDofMap);

    //! Return Lagrange multiplier vector (\f$t_{n+1}\f$)
    virtual Teuchos::RCP<Epetra_Vector> LagrMult() { return z_; }

    /*! \brief Return Lagrange multiplier vector \f$(t_{n+1})\f$
     *
     *  \param redist (in): If TRUE, the redistributed vector is returned,
     *                      otherwise the vector with the original map before
     *                      any redistribution took place.
     *
     *  Please note, that the vector is returned with the slave dof row map,
     *  i.e. actually the wrong map!
     *
     *  \author hiermeier
     *  \date 05/16 */
    virtual Teuchos::RCP<const Epetra_Vector> GetLagrMultNp(const bool& redist) const;

    //! Return old Lagrange multiplier vector (\f$t_{n}\f$)
    virtual Teuchos::RCP<Epetra_Vector> LagrMultOld() { return Data().OldLmPtr(); }

    /*! \brief Return old Lagrange multiplier vector \f$(t_n)\f$
     *
     *  \param redist (in): If TRUE, the redistributed vector is returned,
     *                      otherwise the vector with the original map before
     *                      any redistribution took place.
     *
     *  Please note, that the vector is returned with the slave dof row map,
     *  i.e. actually the wrong map!
     *
     *  \author hiermeier
     *  \date 05/16 */
    virtual Teuchos::RCP<const Epetra_Vector> GetLagrMultN(const bool& redist) const;

    //! Return Lagrange multiplier vector from last Uzawa step
    Teuchos::RCP<Epetra_Vector> LagrMultUzawa() { return Data().LmUzawaPtr(); }

    //! Return constraint rhs vector (only in saddle-point formulation
    Teuchos::RCP<Epetra_Vector> ConstrRhs() { return Data().ConstrRhsPtr(); }

    //! Returns increment of LagrangeMultiplier solution vector in SaddlePointSolve routine
    Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() { return Data().LmIncrPtr(); }
    Teuchos::RCP<const Epetra_Vector> GetLagrMultSolveIncr() const { return Data().LmIncrPtr(); };

    //! Return mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix> DMatrix() { return Data().DMatrixPtr(); }

    //! Return mortar matrix M
    Teuchos::RCP<LINALG::SparseMatrix> MMatrix() { return Data().MMatrixPtr(); }

    //! Return vector of normal contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactNorStress() { return Data().StressNormalPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactNorStress() const { return Data().StressNormalPtr(); }

    //! Return vector of tangential contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactTanStress() { return Data().StressTangentialPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactTanStress() const
    {
      return Data().StressTangentialPtr();
    }

    //! Return vector of normal contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactNorForce() { return Data().ForceNormalPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactNorForce() const { return Data().ForceNormalPtr(); }

    //! Return vector of tangential contact stresses \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> ContactTanForce() { return Data().ForceTangentialPtr(); }
    Teuchos::RCP<const Epetra_Vector> ContactTanForce() const
    {
      return Data().ForceTangentialPtr();
    }


    //! Return required Integration time
    double Inttime() { return Data().IntTime(); };

    //! Set integration time to zero
    void Inttime_init() { Data().IntTime() = 0.0; };

    //! Return current global contact status
    bool IsInContact() const { return Data().IsInContact(); }

    /*! \brief Return old global contact status (this time step)

     True if there has been contact in any nonlinear iteration
     step of the current time step. */
    bool WasInContact() const { return Data().WasInContact(); }

    /*!
     \brief Return old global contact status (last time step)

     True if there has been contact at the end of the last
     time step (last converged state)
     */
    bool WasInContactLastTimeStep() const { return Data().WasInContactLastTimeStep(); }

    /*! \brief Return global self contact status

     Note that at the moment this only gives information about the
     POTENTIAL self contact of the global problem and not about
     an actual self contact occurring.

     TODO: automatically recognize ACTUAL self contact */
    bool& IsSelfContact() { return Data().IsSelfContact(); }
    bool IsSelfContact() const { return Data().IsSelfContact(); };

    //! Return global frictional status
    bool Friction() const { return Data().IsFriction(); }

    //! Return contact interfaces
    const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& ContactInterfaces() const
    {
      return Interfaces();
    }

    /*! \brief Get dual quadratic 3d slave element flag

     Returns TRUE if at least one higher-order 3d slave element with
     dual Lagrange mutliplier shape functions in any interface. */
    virtual bool Dualquadslavetrafo() const { return Data().IsDualQuadSlaveTrafo(); };

    //! Return parallel redistribution status (yes or no)
    inline bool ParRedist() const { return (Data().ParType() != INPAR::MORTAR::parredist_none); }


    //! Return specific parallel redistribution status
    inline INPAR::MORTAR::ParRedist WhichParRedist() const { return Data().ParType(); }

    //! Return matrix T
    virtual Teuchos::RCP<LINALG::SparseMatrix> TMatrix() { return Teuchos::null; }

    //! Return number of active nodes
    int NumberOfActiveNodes() const
    {
      if (not Data().GActiveNodeRowMapPtr().is_null())
        return Data().GActiveNodeRowMapPtr()->NumGlobalElements();
      return 0;
    }

    //! Return number of frictional slip nodes
    int NumberOfSlipNodes() const
    {
      if (not Data().GSlipNodeRowMapPtr().is_null())
        return Data().GSlipNodeRowMapPtr()->NumGlobalElements();
      return 0;
    }

    //! @}

    //! @name Evaluation methods
    //! @{

    /*! \brief Redistribute all contact interfaces in parallel

     We hand in the current global displacement state so that
     a contact search can be performed and set state called.

     \return TRUE if the interface has been redistributed. Return FALSE otherwise.*/
    virtual bool RedistributeContact(Teuchos::RCP<const Epetra_Vector> dis);

    /** \brief Redistribute all contact interfaces in parallel
     *
     *  In contrast to RedistributeContact this routine takes place at a different
     *  point during the simulation. For example, the redistribution can be initiated
     *  each time a certain amount of Newton steps per load step has been reached.
     *  In this way an adaption can be made quicker directly after a large predictor
     *  step or another unforeseen scenario which might have changed the contact
     *  situation severely. */
    virtual bool DynRedistributeContact(
        const Teuchos::RCP<const Epetra_Vector>& dis, const int nlniter)
    {
      return false;
    };

    //! Global evaluation method called from time integrator
    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, const int step,
        const int iter, bool predictor = false);

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
    void Evaluate(CONTACT::ParamsInterface& cparams) { Evaluate(cparams, NULL); }
    void Evaluate(CONTACT::ParamsInterface& cparams,
        const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec)
    {
      Evaluate(cparams, eval_vec, NULL);
    }
    void Evaluate(CONTACT::ParamsInterface& cparams,
        const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec,
        const std::vector<Teuchos::RCP<Epetra_Vector>>* eval_vec_mutable);

    /*! \brief Set current deformation state

     All interfaces are called to set the current deformation state
     (u, xspatial) in their nodes. Additionally, the new contact
     element areas are computed.

     \param statename (in): std::string defining which quantity to set (either "displacement" or
     "olddisplacement") \param vec (in): current global state of the quantity defined by statename
   */
    virtual void SetState(const enum MORTAR::StateType& statename, const Epetra_Vector& vec);

    /*! \brief Evaluate reference state

     for frictional contact we need history values (relative velocity) and
     therefore we store the nodal entries of mortar matrices (reference
     configuration) before the first time step */
    virtual void EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> vec);

    void InitBinStrategyforTimestep(Teuchos::RCP<const Epetra_Vector> vel);

    /*! \brief Evaluate matrix of nodal normals

     This is needed for energy-conserving time integration (Velocity-Update) */
    Teuchos::RCP<LINALG::SparseMatrix> EvaluateNormals(Teuchos::RCP<Epetra_Vector> dis);

    /** \brief Write strategy specific output
     *
     *  \param(in) writer: output writer */
    virtual void WriteOutput(IO::DiscretizationWriter& writer) const { return; }

    //! @}

    //! @name merit function methods
    //! @{

    /// return the potential contributions of the active contact strategy
    virtual double GetPotentialValue(
        const enum NOX::NLN::MeritFunction::MeritFctName mrt_type) const;

    /// return contributions of the active contact strategy to the linear model
    virtual double GetLinearizedPotentialValueTerms(const Epetra_Vector& dir,
        const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
        const enum NOX::NLN::MeritFunction::LinOrder linorder,
        const enum NOX::NLN::MeritFunction::LinType lintype) const;

    //! @}

    //! @name Preconditioner methods
    //! @{

    virtual bool IsSaddlePointSystem() const;

    virtual bool IsCondensedSystem() const;

    /* Fill the maps vector for the linear solver preconditioner
     *
     * The following order is pre-defined:
     * (0) masterDofMap
     * (1) slaveDofMap
     * (2) innerDofMap
     * (3) activeDofMap
     *
     * \author hiermeier */
    virtual void FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map>>& maps) const;

    //! compute the preconditioner operator
    virtual bool computePreconditioner(
        const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams = 0);

    //! @}

    //! @name Quantity control methods
    //! @{
    /*! \brief Get some nodal quantity globally and store into CoNodes

     The enum input parameter defines, which quantity is be updated.
     Currently the possibilities "lmold", "lmcurrent", "lmupdate" and
     "lmuzawa" exist. Note that "lmold" means the converged value LM_n
     of the last time / load step, whereas "lmcurrent" adresses the current
     (not necessarily converged) value of the LM_n+1. "lmupdate" is a special
     option called only in Recover() after the update of the Lagr. multipliers.
     It basically does the same as "lmcurrent", but also checks for D.B.C.
     problems. Finally, "lmuzawa" addresses the LM update within an
     Uzawa augmented Lagrangian scheme.

     \param type (in): enum defining which quantity to store into CoNodes

     */
    virtual void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type);

    /*! \brief Evaluate contact stresses in normal direction and tangential plane

     This is called at the end of each time or load step. It calculates
     the stress vector in normal direction and the stress vector in the
     tangential plane. */
    virtual void OutputStresses();

    /*! \brief Get dirichlet B.C. status and store into CoNodes

     This is called once at the beginning of the simulation
     to set the D.B.C. status in each CNode.

     \param dbcmaps (in): MapExtractor carrying global dbc map */
    virtual void StoreDirichletStatus(Teuchos::RCP<const LINALG::MapExtractor> dbcmaps);

    virtual void SetParentState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec,
        const Teuchos::RCP<DRT::Discretization> dis){
        /* standard contact methods don't need the corresponding bulk element */};

    /*! \brief Update contact at end of time step

     \param dis (in):  current displacements (-> old displacements)

     */
    virtual void Update(Teuchos::RCP<const Epetra_Vector> dis);

    /*! \brief Perform a write restart

     A write restart is initiated by the contact manager. However, the manager has no
     direct access to the nodal quantities. Hence, a portion of the restart has to be
     performed on the level of the contact algorithm, for short: here's the right place.
     */
    virtual void DoWriteRestart(std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors,
        bool forcedrestart = false) const;

    /*! \brief Perform a write restart

     A write restart is initiated by the contact manager. However, the manager has no
     direct access to the nodal quantities. Hence, all the restart action has to be
     performed on the level of the contact algorithm, for short: here's the right place.
     */
    virtual void DoReadRestart(
        IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis)
    {
      DoReadRestart(reader, dis, Teuchos::null);
    };
    virtual void DoReadRestart(IO::DiscretizationReader& reader,
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    //! @}

    //! @name Debugging methods
    //! @{
    /*! \brief Compute interface forces and moments

     Compute current interface forces and moments at n+1-alphaf using current
     Lagrange multiplier values and current Mortar matrices D and M at n+1. When
     doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
     matrices of the last converged time / load step n (TR-like interpolation).

     \param output (in): flag indicating whether force output shall be written
     */
    virtual void InterfaceForces(bool output = false);

    //! Print interfaces
    virtual void Print(std::ostream& os) const;

    //! Print current active set to screen for debugging purposes
    virtual void PrintActiveSet() const;

    /*! \brief Visualize contact stuff with gmsh

     \param step (in): current time step index
     \param iter (in): current iteration index
     */
    void VisualizeGmsh(const int step, const int iter);

    //! @}

    //! @name Purely virtual functions
    //! @{
    // All these functions are defined in one or more specific derived classes,
    // i.e CONTACT::CoLagrangeStrategy or CONTACT::CoPenaltyStrategy.
    // As the base class MORTAR::StrategyBase is always called from the control routine
    // (time integrator), these functions need to be defined purely virtual here.
    virtual bool ActiveSetSemiSmoothConverged() const = 0;
    virtual bool ActiveSetConverged() = 0;
    virtual int ActiveSetSteps() = 0;
    virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const = 0;
    virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const = 0;
    virtual double ConstraintNorm() const = 0;
    virtual void EvaluateContact(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual void EvaluateFriction(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual void EvaluateRelMovPredict() = 0;
    virtual double InitialPenalty() = 0;
    virtual void Initialize() = 0;
    virtual void InitializeUzawa(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi) = 0;
    virtual void ResetActiveSet() = 0;
    virtual void ResetPenalty() = 0;
    virtual void ModifyPenalty() = 0;
    virtual void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,
        Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
        Teuchos::RCP<Epetra_Vector>& blockrhs) = 0;
    virtual void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol) = 0;
    virtual void EvalConstrRHS() = 0;
    virtual void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) = 0;
    virtual void UpdateActiveSet() = 0;
    virtual void UpdateActiveSetSemiSmooth(const bool firstStepPredictor = false) = 0;
    virtual void UpdateUzawaAugmentedLagrange() = 0;
    virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;

    //! @}

    //! @name Empty functions (meshtying)
    //! @{
    // All these functions only have functionality in meshtying simulations, thus they
    // are defined as empty here in the case of contact. They can be called from the
    // control routine (time integrator), whenever you like.

    void RedistributeMeshtying() {}
    void RestrictMeshtyingZone() {}
    void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis)
    {
    }
    Teuchos::RCP<Epetra_Vector> MeshInitialization() { return Teuchos::null; };

    void MortarCoupling(const Teuchos::RCP<const Epetra_Vector>& dis) {}
    //! @}

   protected:
    //! @name Pre/Postoperators
    //! @{
    //! Run after the StoreDirichletStatus() routine has been called
    virtual void PostStoreDirichletStatus(Teuchos::RCP<const LINALG::MapExtractor> dbcmaps){};

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
    //! @}

    //! @Internal evaluate routines
    //! @{
    /*! \brief Evaluate only the weighted gap values
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time
     * integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void EvalWeightedGap(CONTACT::ParamsInterface& cparams);

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
     * integration
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
     * instead.
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time
     * integration \param xold    (in): old solution vector of the NOX solver \param dir     (in):
     * current search direction (in general NOT the actual step, keep in mind that the step length
     * can differ from 1.0) \param xnew    (in): new solution vector of the NOX solver
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
     *  This routine is called in the end of a NOX::Solver::step() call.
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
     * integration \param xnew    (in): new solution vector of the NOX solver
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual void ResetLagrangeMultipliers(
        const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew);

    /// \brief Evaluate the weighted gap gradient error
    virtual void EvalWeightedGapGradientError(CONTACT::ParamsInterface& cparams);

    virtual void CorrectParameters(
        CONTACT::ParamsInterface& cparams, const NOX::NLN::CorrectionType type);

    /// remove condensed contact contributions from the structural rhs
    /**
     *  \param(in) str_rhs: reference to the structural rhs
     *  \author hiermeier \date 03/18 */
    virtual void RemoveCondensedContributionsFromRhs(Epetra_Vector& str_rhs) const;

    //! @}
   protected:
    //! access the contact interfaces of the concrete strategies (read and write)
    virtual std::vector<Teuchos::RCP<CONTACT::CoInterface>>& Interfaces() = 0;

    //! access the contact interfaces of the concrete strategies (read-only)
    virtual const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& Interfaces() const = 0;

    /*! \brief Evaluate contact

     This is just a tiny control routine, deciding which Evaluate-routine
     of those listed below is to be called (based on input-file information).
     Note that into ALL derived Evaluate() routines, a REFERENCE to the pointer
     on the effective stiffness matrix is handed in. This way, after building the
     new effective stiffness matrix with contact, we can simply let the pointer
     kteff point onto the new object. The same is true for the effective force
     vector feff. Be careful: kteff is of type Teuchos::RCP<LINALG::SparseOperator>&.

     \param kteff (in/out): effective stiffness matrix (without -> with contact)
     \param feff (in/out): effective residual / force vector (without -> with contact)

     */
    void Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff,
        Teuchos::RCP<Epetra_Vector> dis);

    /*! \brief Evaluate relative movement of contact bodies

     This is for evaluating the relative movement of contact bodies. This
     can either be done with regarding the different movement of material points
     or regarding the change of mortar projection. The second possibility
     is definitely objective wheras the first possibility is objective
     only when the gap is zero. */
    void EvaluateRelMov();

    /*! \brief Initialize and evaluate Mortar stuff for the next Newton step

     This method first checks if we are dealing with self contact and updates
     the interface slave and master sets if so. Then it resets the global
     Mortar matrices D and M and the global gap vector g accordingly.

     The nodal quantites computed in InitEvalInterface() are then assembled
     to global matrices and vectors respectively. No setup of the global system
     is to be done here yet, so there is no need to pass in the effective
     stiffness K or the effective load vector f. */
    virtual void InitMortar();
    virtual void AssembleMortar();

    /*! \brief Initialize and evaluate interface for the next Newton step

     This method calls Initialize() on all contact interfaces, which
     resets all kind of nodal quantities like normal vector, weighted
     gap or Mortar and linearization maps. It then calls Evaluate() on
     all contact interfaces, which does all the geometric contact stuff.
     Concretely, this is an evaluation of all involved quantites at nodal
     level plus the setup of all corresponding linearizations.
     It includes the nodal normal calculations, contact search, projection
     and overlap detection, integration of the  Mortar terms D, M and of the
     weighted gap. Additionally, the linearizations of geometrical quantities
     (delta_n, delta_t, delta_D, delta_M) are calculated. */
    void InitEvalInterface() { InitEvalInterface(Teuchos::null); };
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
     completion of a time / load step!). The std::string "current" adresses the
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
    /*! \brief Create the global Lagrange multiplier DoF row map
     *
     *  The global Lagrange multiplier DoF row map is created in a deterministic
     *  manner based on the previously created global slave DoF row map. This is
     *  necessary for the later ReplaceMap calls. Especially, the std::sort during
     *  a LINALG::MergeMap call would otherwise destroy the correlation. This becomes
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

    //! \brief Calculate mean velocity of interface for binning range
    void CalcMeanVelforBinning(Teuchos::RCP<const Epetra_Vector> vel);

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
    // don't want cctor (= operator impossible anyway for abstract class)
    CoAbstractStrategy(const CoAbstractStrategy& old);

    /*! @name References to the data container content
     *
     * \remark Please add no new member variables to the abstract strategy! Use
     *  the corresponding data container instead (--> CONTACT::AbstractStratDataContainer).
     *  If you have any questions concerning this, do not hesitate and ask me.
     *                                                                   hiermeier 05/16 */
    //! @{
    std::vector<double>& tunbalance_;  //!< ref. to parallel unbalance factors (time) for current
                                       //!< time step \f$t_{n+1}\f$
    std::vector<int>& eunbalance_;  //!< ref. to parallel unbalance factors (eles) for current time
                                    //!< step \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Map>&
        glmdofrowmap_;  //!< ref. to global Lagrange mult. dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gsnoderowmap_;  //!< ref. to global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gmnoderowmap_;  //!< ref. to global master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gsdofrowmap_;  //!< ref. to global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gmdofrowmap_;  //!< ref. to global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>& gndofrowmap_;  //!< ref. to global internal dof row map
    Teuchos::RCP<Epetra_Map>&
        gsmdofrowmap_;  //!< ref. to global slave and master dof row map (s+m map)
    Teuchos::RCP<Epetra_Map>&
        gdisprowmap_;  //!< ref. to global displacement dof row map (s+m+n map)

    Teuchos::RCP<Epetra_Map>&
        gactivenodes_;  //!< ref. to global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gactivedofs_;  //!< ref. to global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        ginactivenodes_;  //!< ref. to global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        ginactivedofs_;  //!< ref. to global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gactiven_;  //!< ref. to global dof row map of matrix N (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gactivet_;  //!< ref. to global dof row map of matrix T (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gslipnodes_;  //!< ref. to global slip slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gslipdofs_;  //!< ref. to global slip slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map>&
        gslipt_;  //!< ref. to global row map of matrix T for slip dofs (of all interfaces)

    Teuchos::RCP<Epetra_Map>& gsdofVertex_;  //!< ref. to global slave row map of vertex nodes
    Teuchos::RCP<Epetra_Map>& gsdofEdge_;    //!< ref. to global slave row map of edge nodes
    Teuchos::RCP<Epetra_Map>& gsdofSurf_;    //!< ref. to global slave row map of surface nodes

    Teuchos::RCP<Epetra_Map>&
        pglmdofrowmap_;  //!< ref. to global LM dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map>&
        pgsdofrowmap_;  //!< ref. to global slave dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map>&
        pgmdofrowmap_;  //!< ref. to global master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map>& pgsmdofrowmap_;       //!< ref. to global slave and master dof row map
                                                    //!< (before parallel redistribution)
    Teuchos::RCP<Epetra_Vector>& pgsdirichtoggle_;  //!< ref. to global dirichlet toggle of all
                                                    //!< slave dofs (before parallel redistribution)

    std::vector<Teuchos::RCP<Epetra_Map>>&
        initial_elecolmap_;  //!< initial col ele map for binning strategy (s m)

    Teuchos::RCP<LINALG::SparseMatrix>& dmatrix_;  //!< ref. to global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix>& mmatrix_;  //!< ref. to global Mortar matrix M
    Teuchos::RCP<Epetra_Vector>& g_;               //!< ref. to global weighted gap vector g
    Teuchos::RCP<Epetra_Vector>&
        tangrhs_;  //!< ref. to global tangential rhs vector (formulation with incremental z_)
    Teuchos::RCP<Epetra_Vector>& inactiverhs_;  //!< ref. to gloabl inactive rhs vector (formulation
                                                //!< with incremental z_ and saddle point system)
    Teuchos::RCP<Epetra_Vector>&
        strcontactrhs_;  //!< ref. to global structural contact rhs contributions at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>&
        constrrhs_;  //!< ref. to global constraint rhs vector (only for saddlepoint problems)

    Teuchos::RCP<LINALG::SparseMatrix>&
        lindmatrix_;  //!< ref. to global Matrix LinD containing slave fc derivatives
    Teuchos::RCP<LINALG::SparseMatrix>&
        linmmatrix_;  //!< ref. to global Matrix LinM containing master fc derivatives
    Teuchos::RCP<LINALG::SparseMatrix>& kteffnew_;

    Teuchos::RCP<LINALG::SparseMatrix>&
        dold_;  //!< ref. to global Mortar matrix D (last end-point \f$t_{n}\f$)
    Teuchos::RCP<LINALG::SparseMatrix>&
        mold_;  //!< ref. to global Mortar matrix D (last end-point \f$t_{n}\f$)

    Teuchos::RCP<Epetra_Vector>&
        z_;  //!< ref. to current vector of Lagrange multipliers at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& zold_;  //!< ref. to old vector of Lagrange multipliers at t_n
    Teuchos::RCP<Epetra_Vector>&
        zincr_;  //!< ref. to Lagrange multiplier vector increment within SaddlePointSolve (this is
                 //!< NOT the increment of z_ between \f$t_{n+1}\f$ and \f$t_{n}\f$!)
    Teuchos::RCP<Epetra_Vector>&
        zuzawa_;  //!< ref. to vector of Lagrange multipliers from last Uzawa step

    Teuchos::RCP<Epetra_Vector>&
        stressnormal_;  //!< ref. to vector of normal contact forces at t_n+1
    Teuchos::RCP<Epetra_Vector>&
        stresstangential_;  //!< ref. to vector of tangential contact forces at t_n+1
    Teuchos::RCP<Epetra_Vector>&
        forcenormal_;  //!< ref. to vector of normal contact forces at t_n+1
    Teuchos::RCP<Epetra_Vector>&
        forcetangential_;  //!< ref. to vector of tangential contact forces at t_n+1

    int& step_;  //!< ref. to time step index
    int& iter_;  //!< ref. to nonlinear (Newton) iteration index

    bool& isincontact_;       //!< ref. to flag indicating global contact status
    bool& wasincontact_;      //!< ref. to flag indicating global contact status of this time step
                              //!< (history)
    bool& wasincontactlts_;   //!< ref. to flag indicating global contact status of last time step
    bool& isselfcontact_;     //!< ref. to flag indicating potential self contact
    bool& friction_;          //!< ref. to flag for frictional contact
    bool& nonSmoothContact_;  //!< ref. to flag for nonsmooth contact algorithm
    bool& regularized_;       //!< ref. to flag for regularized contact

    bool& dualquadslavetrafo_;  //!< ref. to flag indicating whether trafo should be applied
    Teuchos::RCP<LINALG::SparseMatrix>&
        trafo_;  //!< ref. to transformation matrix T for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix>
        systrafo_;  //! transformation matrix T for dual quad 3D case (all problem dofs)
    Teuchos::RCP<LINALG::SparseMatrix>
        invsystrafo_;  //! inverse transformation matrix T for dual quad 3D case (all problem dofs)
    Teuchos::RCP<LINALG::SparseMatrix>&
        invtrafo_;  //!< ref. to inverse trafo matrix T^(-1) for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix>& dmatrixmod_;  //!< ref. to modified global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix>& doldmod_;  //!< ref. to modified global Mortar matrix Dold

    double& inttime_;            //!< ref. to integration time
    std::vector<double>& ivel_;  //!< ref. to mean interface velocity

    INPAR::CONTACT::SolvingStrategy& stype_;  //!< ref. to current used solving strategy
    INPAR::CONTACT::ConstraintDirection&
        constr_direction_;  //!< ref. to direction in which the contact constraints are formulated
   private:
    //! pointer to the data container object
    Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr_;

    //! pointer to the NOX::NLN::CONSTRAINT::Interface::Required object
    Teuchos::RCP<CONTACT::NoxInterface> noxinterface_ptr_;

  };  // class CoAbstractStrategy
}  // namespace CONTACT

// << operator
std::ostream& operator<<(std::ostream& os, const CONTACT::CoAbstractStrategy& strategy);

#endif  // #ifndef CONTACT_STRATEGY_H
