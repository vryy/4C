/*----------------------------------------------------------------------*/
/*!
\file contact_augmented_strategy.H

\brief Augmented contact solving strategy with standard Lagrangian
       multipliers.

\level 2

\maintainer Michael Hiermeier

\date Apr 7, 2014

*/
/*----------------------------------------------------------------------*/

#ifndef CONTACT_AUGMENTED_STRATEGY_H_
#define CONTACT_AUGMENTED_STRATEGY_H_

#include "../drt_contact/contact_lagrange_strategy.H"
#include "contact_aug_utils.H"
#include "contact_aug_enum_lists.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_singletondestruction.H"
#include "contact_aug_timemonitor.H"

// forward declaration
namespace DRT
{
  class Node;
}  // namespace DRT
namespace MORTAR
{
  class MatrixRowColTransformer;
}  // namespace MORTAR
namespace CONTACT
{
  namespace AUG
  {
    namespace STEEPESTASCENT
    {
      class DataContainer;
    }  // namespace STEEPESTASCENT
    class Interface;
    class Potential;
    namespace POTENTIAL
    {
      enum class Type : int;
      enum class SetType : int;
    }  // namespace POTENTIAL
    class ComboStrategy;
    class ActiveSet;
    class ParallelDistributionController;

    /*--------------------------------------------------------------------------*/
    class DataContainer : public ::CONTACT::AbstractStratDataContainer
    {
     public:
      //! @name Constructors and destructors and related methods
      //! @{
      //! Standard constructor
      DataContainer();

      //! Destructor
      virtual ~DataContainer() { return; }

      //! initialize optional sub-data container
      void InitSubDataContainer(const INPAR::CONTACT::SolvingStrategy strat_type);

      //! @}

      //! @name Access methods
      //! @{

      //! @name Sub-data container accessors
      //! @{

      CONTACT::AUG::STEEPESTASCENT::DataContainer& SaData()
      {
        if (sa_data_ptr_.is_null())
          dserror("The steepest ascent sub data container was not initialized!");
        return *sa_data_ptr_;
      }

      const CONTACT::AUG::STEEPESTASCENT::DataContainer& SaData() const
      {
        if (sa_data_ptr_.is_null())
          dserror("The steepest ascent sub data container was not initialized!");
        return *sa_data_ptr_;
      }

      //! @}

      //! @name Booleans accessors
      //! @{

      //! Return WasInContactLastIter indicator
      bool& WasInContactLastIter() { return wasincontactlastiter_; }
      bool WasInContactLastIter() const { return wasincontactlastiter_; }

      //! Return TRUE if the active set is converged, otherwise false
      bool& IsActiveSetConverged() { return isactivesetconverged_; };
      bool IsActiveSetConverged() const { return isactivesetconverged_; }

      //! print the linear conservation check results if TRUE
      bool& PrintLinearMomConservation() { return printlinearconservation_; };
      bool PrintLinearMomConservation() const { return printlinearconservation_; };

      //! print the angular conservation check results if TRUE
      bool& PrintAngularMomConservation() { return printangularconservation_; };
      bool PrintAngularMomConservation() const { return printangularconservation_; };

      bool& AddInactivForceContributions() { return add_inactive_force_; };
      bool AddInactivForceContributions() const { return add_inactive_force_; };

      //! set the the semi smooth newton flag
      void SetIsSemiSmoothNewton(bool is_semi_smooth) { is_semi_smooth_newton_ = is_semi_smooth; }

      //! get the semi smooth newton flag
      bool IsSemiSmoothNewton() const { return is_semi_smooth_newton_; }

      //! Set the matrix maps status
      void SetMatrixMapsValid(bool isvalid) { matrix_maps_valid_ = isvalid; }

      //! Are the matrix maps valid or did they change?
      bool MatrixMapsValid() const { return matrix_maps_valid_; }

      //! Set the vector maps status
      void SetVectorMapsValid(bool isvalid) { vector_maps_valid_ = isvalid; }

      //! Are the vector maps valid or did they change?
      bool VectorMapsValid() const { return vector_maps_valid_; }

      //! @}

      //! set constant semi smooth cn (coming from the Input file)
      void SetConstantCn(double cn) { cn_ = cn; }

      //! access constant semi smooth cn (coming from the Input file)
      double ConstantCn() const
      {
        if (cn_ < 0.0) dserror("cn was not set correctly! (cn_ = %d)", cn_);
        return cn_;
      };

      void SetCurrentEvalState(MORTAR::ActionType eval_state) { eval_state_ = eval_state; }

      MORTAR::ActionType GetCurrentEvalState() const { return eval_state_; }

      //! get the parallel strategy enumerator
      void SetParallelStrategy(INPAR::MORTAR::ParallelStrategy parallel_strat)
      {
        parallel_strategy_ = parallel_strat;
      }

      //! get the parallel strategy enumerator
      INPAR::MORTAR::ParallelStrategy ParallelStrategy() const { return parallel_strategy_; }

      void SetVariationalApproachType(const enum INPAR::CONTACT::VariationalApproach var_type)
      {
        var_type_ = var_type;
      }

      enum INPAR::CONTACT::VariationalApproach VariationalApproachType() const { return var_type_; }

      void SetFDCheckType(const enum INPAR::CONTACT::FDCheck fd_check_type)
      {
        fd_check_type_ = fd_check_type;
      }

      enum INPAR::CONTACT::FDCheck FDCheckType() const { return fd_check_type_; }

      void SetPotential(const Teuchos::RCP<Potential>& potential) { potentialPtr_ = potential; }

      CONTACT::AUG::Potential& Potential() { return *potentialPtr_; }

      const CONTACT::AUG::Potential& Potential() const { return *potentialPtr_; }

      CONTACT::AUG::Potential& MutablePotential() const { return *potentialPtr_; }

      void SetPDController(
          const Teuchos::RCP<CONTACT::AUG::ParallelDistributionController>& pd_control)
      {
        pd_control_ = pd_control;
      }

      CONTACT::AUG::ParallelDistributionController& PDController()
      {
        if (pd_control_.is_null()) dserror("ParallelDistributionController is not initialized!");

        return *pd_control_;
      }

      void SetMatrixRowColTransformer(
          const Teuchos::RCP<MORTAR::MatrixRowColTransformer>& transformer)
      {
        mat_row_col_transformer_ = transformer;
      }

      void InitMatrixRowColTransformer();

      MORTAR::MatrixRowColTransformer& MatrixRowColTransformer() const
      {
        if (mat_row_col_transformer_.is_null())
          dserror("The mat_row_col_transformer was not initialized correctly!");

        return *mat_row_col_transformer_;
      }

      double TotalGradientError() const { return grad_error_.total_; }

      double& TotalGradientError() { return grad_error_.total_; }

      const std::vector<std::pair<int, double>>& NodalGradientError_MaProj() const
      {
        return grad_error_.master_;
      }

      std::vector<std::pair<int, double>>& NodalGradientError_MaProj()
      {
        return grad_error_.master_;
      }

      const std::vector<std::pair<int, double>>& NodalGradientError_Jacobian() const
      {
        return grad_error_.jacobian_;
      }

      std::vector<std::pair<int, double>>& NodalGradientError_Jacobian()
      {
        return grad_error_.jacobian_;
      }

      Teuchos::RCP<Epetra_Vector>& GSeleEvalTimesPtr() { return gSeleEvalTimesPtr_; }

      //! @name Matrix accessors
      //! @{

      /**! \brief Return BMatrix (combination of D and M)
       *
       *   For the complete variational approach the B-Matrix represents the
       *   transpose of the gradient of the weighted gap vector and is
       *   equal to dLmNWGapLin. In the incomplete case the B-Matrix represents
       *   an incomplete estimate of the gradient. */
      Teuchos::RCP<LINALG::SparseMatrix>& BMatrixPtr() { return BMatrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> BMatrixPtr() const { return BMatrixPtr_; }
      LINALG::SparseMatrix& BMatrix()
      {
        if (BMatrixPtr_.is_null()) dserror("The BMatrixPtr_ is not initialized!");
        return *BMatrixPtr_;
      }

      //! Return dGLmLinMatrix
      Teuchos::RCP<LINALG::SparseMatrix>& DGLmLinMatrixPtr() { return dGLmLinMatrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> DGLmLinMatrixPtr() const
      {
        return dGLmLinMatrixPtr_;
      }
      LINALG::SparseMatrix& DGLmLinMatrix()
      {
        if (dGLmLinMatrixPtr_.is_null()) dserror("The dGLmLinMatrixPtr_ is not initialized!");
        return *dGLmLinMatrixPtr_;
      }

      //! Return dGGLinMatrix
      Teuchos::RCP<LINALG::SparseMatrix>& DGGLinMatrixPtr() { return dGGLinMatrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> DGGLinMatrixPtr() const { return dGGLinMatrixPtr_; }
      LINALG::SparseMatrix& DGGLinMatrix()
      {
        if (dGGLinMatrixPtr_.is_null()) dserror("The dGGLinMatrixPtr_ is not initialized!");
        return *dGGLinMatrixPtr_;
      }

      //! Return dLmNWGapLinMatrix
      Teuchos::RCP<LINALG::SparseMatrix>& DLmNWGapLinMatrixPtr() { return dLmNWGapLinMatrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> DLmNWGapLinMatrixPtr() const
      {
        return dLmNWGapLinMatrixPtr_;
      }
      LINALG::SparseMatrix& DLmNWGapLinMatrix()
      {
        if (dLmNWGapLinMatrixPtr_.is_null())
          dserror("The dLmNWGapLinMatrixPtr_ is not initialized!");
        return *dLmNWGapLinMatrixPtr_;
      }

      //! Return dLmTLmTMatrix
      Teuchos::RCP<LINALG::SparseMatrix>& DLmTLmTMatrixPtr() { return dLmTLmTMatrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> DLmTLmTMatrixPtr() const
      {
        return dLmTLmTMatrixPtr_;
      }
      LINALG::SparseMatrix& DLmTLmTMatrix()
      {
        if (dLmTLmTMatrixPtr_.is_null()) dserror("The dLmTLmTMatrixPtr_ is not initialized!");
        return *dLmTLmTMatrixPtr_;
      }

      //! Return dLmTLmTLinMatrix
      Teuchos::RCP<LINALG::SparseMatrix>& DLmTLmTLinMatrixPtr() { return dLmTLmTLinMatrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> DLmTLmTLinMatrixPtr() const
      {
        return dLmTLmTLinMatrixPtr_;
      }
      LINALG::SparseMatrix& DLmTLmTLinMatrix()
      {
        if (dLmTLmTLinMatrixPtr_.is_null()) dserror("The dLmTLmTLinMatrixPtr_ is not initialized!");
        return *dLmTLmTLinMatrixPtr_;
      }

      Teuchos::RCP<LINALG::SparseMatrix>& InactiveLinMatrixPtr() { return inactiveLinMatrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> InactiveLinMatrixPtr() const
      {
        return inactiveLinMatrixPtr_;
      }
      LINALG::SparseMatrix& InactiveLinMatrix()
      {
        if (inactiveLinMatrixPtr_.is_null())
          dserror("The augInactiveLinMatrixPtr_ is not initialized!");
        return *inactiveLinMatrixPtr_;
      }

      Teuchos::RCP<LINALG::SparseMatrix>& InactiveDDMatrixPtr() { return inactive_dd_matrixPtr_; }
      Teuchos::RCP<const LINALG::SparseMatrix> InactiveDDMatrixPtr() const
      {
        return inactive_dd_matrixPtr_;
      }
      LINALG::SparseMatrix& InactiveDDMatrix()
      {
        if (inactive_dd_matrixPtr_.is_null())
          dserror("The inactive_dd_matrixPtr_ is not initialized!");
        return *inactive_dd_matrixPtr_;
      }
      //! @}

      //! @name Vector accessors
      //! @{

      /*! \brief Return a vector of the tributary (weighted) nodal areas
       *
       *  These quantities are integrated over the whole slave side! */
      Teuchos::RCP<Epetra_Vector>& AVecPtr() { return aPtr_; }
      Teuchos::RCP<const Epetra_Vector> AVecPtr() const { return aPtr_; }
      Epetra_Vector& AVec()
      {
        if (aPtr_.is_null()) dserror("The augAPtr_ is not initialized!");
        return *aPtr_;
      }

      //! This vector contributes as diagonal inactive matrix to the kzz block
      Teuchos::RCP<Epetra_Vector>& InactiveDiagMatrixPtr() { return inactiveDiagMatrixPtr_; }
      Teuchos::RCP<const Epetra_Vector> InactiveDiagMatrixPtr() const
      {
        return inactiveDiagMatrixPtr_;
      }
      Epetra_Vector& InactiveDiagMatrix()
      {
        if (inactiveDiagMatrixPtr_.is_null())
          dserror("The inactiveDiagMatrixPtr_ is not initialized!");
        return *inactiveDiagMatrixPtr_;
      }

      /*! \brief Return a vector of the tributary (weighted) nodal areas
       *
       *  These quantities are integrated over the projectable slave side! */
      Teuchos::RCP<Epetra_Vector>& KappaVecPtr() { return kappaPtr_; }
      Teuchos::RCP<const Epetra_Vector> KappaVecPtr() const { return kappaPtr_; }
      Epetra_Vector& KappaVec()
      {
        if (kappaPtr_.is_null()) dserror("The kappaPtr_ is not initialized!");
        return *kappaPtr_;
      }

      //! Return a vector of the lagrange multipliers for the frictionless/normal constraint
      Teuchos::RCP<Epetra_Vector>& LmNPtr() { return lmNPtr_; }
      Teuchos::RCP<const Epetra_Vector> LmNPtr() const { return lmNPtr_; }
      Epetra_Vector& LmN()
      {
        if (lmNPtr_.is_null()) dserror("The lmNPtr_ is not initialized!");
        return *lmNPtr_;
      }

      //! Return the averaged weighted gap
      Teuchos::RCP<Epetra_Vector>& AWGapPtr() { return aWGapPtr_; }
      Teuchos::RCP<const Epetra_Vector> AWGapPtr() const { return aWGapPtr_; }
      Epetra_Vector& AWGap()
      {
        if (aWGapPtr_.is_null()) dserror("The aWGapRhsPtr_ is not initialized!");
        return *aWGapPtr_;
      }

      //! return global weighted gap vector of all slave nodes
      Teuchos::RCP<Epetra_Vector>& WGapAllSlNodesPtr() { return wGapAllPtr_; };
      Teuchos::RCP<const Epetra_Vector> WGapAllSlNodesPtr() const { return wGapAllPtr_; };
      Epetra_Vector& WGapAllSlNodes()
      {
        if (wGapAllPtr_.is_null()) dserror("The wgap_all_ vector is not initialized!");
        return *wGapAllPtr_;
      }

      /*! \brief Return the tangential Lagrange multiplier right hand side
       *
       *  In the frictionless case this is supposed to do nothing and just keeps
       *  the tangential Lagrange multipliers at a constant value of zero. */
      Teuchos::RCP<Epetra_Vector>& DLmTLmTRhsPtr() { return dLmTLmTRhsPtr_; }
      Teuchos::RCP<const Epetra_Vector> DLmTLmTRhsPtr() const { return dLmTLmTRhsPtr_; }
      Epetra_Vector& DLmTLmTRhs()
      {
        if (dLmTLmTRhsPtr_.is_null()) dserror("The dLmTLmTRhsPtr_ is not initialized!");
        return *dLmTLmTRhsPtr_;
      }

      /*! \brief Return the slave force due to active Lagrange multiplier values */
      Teuchos::RCP<Epetra_Vector>& SlForceLmPtr() { return slForceLmPtr_; }
      Teuchos::RCP<const Epetra_Vector> SlForceLmPtr() const { return slForceLmPtr_.getConst(); }
      Epetra_Vector& SlForceLm()
      {
        if (slForceLmPtr_.is_null()) dserror("The slForceLmPtr_ is not initialized!");
        return *slForceLmPtr_;
      }

      /*! \brief Return the slave force due to inactive Lagrange multiplier values */
      Teuchos::RCP<Epetra_Vector>& SlForceLmInactivePtr() { return slForceLmInactivePtr_; }
      Teuchos::RCP<const Epetra_Vector> SlForceLmInactivePtr() const
      {
        return slForceLmInactivePtr_.getConst();
      }
      Epetra_Vector& SlForceLmInactive()
      {
        if (slForceLmInactivePtr_.is_null())
          dserror("The slForceLmInactivePtr_ is not initialized!");
        return *slForceLmInactivePtr_;
      }

      /*! \brief Return the slave force due to gap values */
      Teuchos::RCP<Epetra_Vector>& SlForceGPtr() { return slForceGPtr_; }
      Teuchos::RCP<const Epetra_Vector> SlForceGPtr() const { return slForceGPtr_.getConst(); }
      Epetra_Vector& SlForceG()
      {
        if (slForceGPtr_.is_null()) dserror("The slForceGPtr_ is not initialized!");
        return *slForceGPtr_;
      }

      /*! \brief Return the master force due to Lagrange multiplier values */
      Teuchos::RCP<Epetra_Vector>& MaForceLmPtr() { return maForceLmPtr_; }
      Teuchos::RCP<const Epetra_Vector> MaForceLmPtr() const { return maForceLmPtr_.getConst(); }
      Epetra_Vector& MaForceLm()
      {
        if (maForceLmPtr_.is_null()) dserror("The maForceLmPtr_ is not initialized!");
        return *maForceLmPtr_;
      }

      /*! \brief Return the master force due to gap values */
      Teuchos::RCP<Epetra_Vector>& MaForceGPtr() { return maForceGPtr_; }
      Teuchos::RCP<const Epetra_Vector> MaForceGPtr() const { return maForceGPtr_.getConst(); }
      Epetra_Vector& MaForceG()
      {
        if (maForceGPtr_.is_null()) dserror("The maForceGPtr_ is not initialized!");
        return *maForceGPtr_;
      }

      //! Returns the nodal c_n vector
      Teuchos::RCP<Epetra_Vector>& CnPtr() { return cnPtr_; }
      Teuchos::RCP<const Epetra_Vector> CnPtr() const { return cnPtr_; }
      Epetra_Vector& Cn()
      {
        if (cnPtr_.is_null()) dserror("The cnPtr_ is not initialized!");
        return *cnPtr_;
      }
      const Epetra_Vector& Cn() const
      {
        if (cnPtr_.is_null()) dserror("The cnPtr_ is not initialized!");
        return *cnPtr_;
      }

      //! @}

      //! @name Global map accessors
      //! @{
      Teuchos::RCP<Epetra_Map>& GSlNormalDofRowMapPtr() { return gsndofrowmapPtr_; }
      Teuchos::RCP<const Epetra_Map> GSlNormalDofRowMapPtr() const { return gsndofrowmapPtr_; }
      Epetra_Map& GSlNormalDofRowMap()
      {
        if (gsndofrowmapPtr_.is_null()) dserror("The gsndofrowmapPtr_ is not initialized!");
        return *gsndofrowmapPtr_;
      }

      Teuchos::RCP<Epetra_Map>& GSlTangentialDofRowMapPtr() { return gstdofrowmapPtr_; }
      Teuchos::RCP<const Epetra_Map> GSlTangentialDofRowMapPtr() const { return gstdofrowmapPtr_; }
      Epetra_Map& GSlTangentialDofRowMap()
      {
        if (gstdofrowmapPtr_.is_null()) dserror("The gstdofrowmapPtr_ is not initialized!");
        return *gstdofrowmapPtr_;
      }

      Teuchos::RCP<Epetra_Map>& GOldActiveSlaveNodesPtr() { return gOldActiveSlaveNodesPtr_; }
      Teuchos::RCP<const Epetra_Map> GOldActiveSlaveNodesPtr() const
      {
        return gOldActiveSlaveNodesPtr_;
      }
      Epetra_Map& GOldActiveSlaveNodes()
      {
        if (gOldActiveSlaveNodesPtr_.is_null())
          dserror("The gAugOldActiveSlaveNodesPtr_ is not initialized!");
        return *gOldActiveSlaveNodesPtr_;
      }

      Teuchos::RCP<Epetra_Map>& GSeleRowMapPtr() { return gSeleRowMapPtr_; }
      Teuchos::RCP<const Epetra_Map> GSeleRowMapPtr() const { return gSeleRowMapPtr_; }
      Epetra_Map& GSeleRowMap()
      {
        if (gSeleRowMapPtr_.is_null()) dserror("The gSeleRowMapPtr_ is not initialized!");
        return *gSeleRowMapPtr_;
      }

      Teuchos::RCP<Epetra_Map>& GSeleColMapPtr() { return gSeleColMapPtr_; }
      Teuchos::RCP<const Epetra_Map> GSeleColMapPtr() const { return gSeleColMapPtr_; }
      Epetra_Map& GSeleColMap()
      {
        if (gSeleColMapPtr_.is_null()) dserror("The gSeleColMapPtr_ is not initialized!");
        return *gSeleColMapPtr_;
      }

      Teuchos::RCP<Epetra_Map>& GMeleRowMapPtr() { return gMeleRowMapPtr_; }
      Teuchos::RCP<const Epetra_Map> GMeleRowMapPtr() const { return gMeleRowMapPtr_; }
      Epetra_Map& GMeleRowMap()
      {
        if (gMeleRowMapPtr_.is_null()) dserror("The gMeleRowMapPtr_ is not initialized!");
        return *gMeleRowMapPtr_;
      }

      Teuchos::RCP<Epetra_Map>& GMeleColMapPtr() { return gMeleColMapPtr_; }
      Teuchos::RCP<const Epetra_Map> GMeleColMapPtr() const { return gMeleColMapPtr_; }
      Epetra_Map& GMeleColMap()
      {
        if (gMeleColMapPtr_.is_null()) dserror("The gMeleColMapPtr_ is not initialized!");
        return *gMeleColMapPtr_;
      }
      //! @}
      //! @}

     protected:
      // don't want = operator and cctor
      DataContainer operator=(const DataContainer& old);
      DataContainer(const DataContainer& old);

      //! flag indicating global contact status of the last nonlinear iteration step
      bool wasincontactlastiter_;

      //! flag indicating whether the active set is converged or not (semi-smooth Newton case only)
      bool isactivesetconverged_;

      //! print the linear conservation check results if TRUE
      bool printlinearconservation_;

      //! print the angular conservation check results if TRUE
      bool printangularconservation_;

      //! flag indicating if a semi smooth newton is used
      bool is_semi_smooth_newton_;

      //! consider forces from Lagrange multiplier values of inactive nodes
      bool add_inactive_force_ = false;

      //! flag indicating if the slave and master matrix maps changed ( e.g. due to redistribution )
      bool matrix_maps_valid_;

      //! flag indicating if the slave and master vector maps changed ( e.g. due to redistribution )
      bool vector_maps_valid_;

      bool ispredict_ = false;

      //! constant cn value from the input file
      double cn_;

      /// access the current evaluation state of the contact strategy
      enum MORTAR::ActionType eval_state_;

      /// defines the type of the parallel evaluation strategy
      enum INPAR::MORTAR::ParallelStrategy parallel_strategy_;

      /// which variational approach shall be followed? (complete, incomplete, ...)
      enum INPAR::CONTACT::VariationalApproach var_type_;

      /// global or gp/nodal finite difference check?
      enum INPAR::CONTACT::FDCheck fd_check_type_;

      /// pointer to the potential object
      Teuchos::RCP<CONTACT::AUG::Potential> potentialPtr_;

      /// row column parallel transformer for matrices
      Teuchos::RCP<MORTAR::MatrixRowColTransformer> mat_row_col_transformer_;

      /// parallel redistribution controller
      Teuchos::RCP<CONTACT::AUG::ParallelDistributionController> pd_control_;

      //! global matrix B-matrix (slave + master)
      Teuchos::RCP<LINALG::SparseMatrix> BMatrixPtr_;

      //! global matrix dGLmLin (slave + master)
      Teuchos::RCP<LINALG::SparseMatrix> dGLmLinMatrixPtr_;

      //! global matrix dGGLin (slave + master)
      Teuchos::RCP<LINALG::SparseMatrix> dGGLinMatrixPtr_;

      //! global matrix dLmNWGapLin
      Teuchos::RCP<LINALG::SparseMatrix> dLmNWGapLinMatrixPtr_;

      //! global matrix dLmTLmT
      Teuchos::RCP<LINALG::SparseMatrix> dLmTLmTMatrixPtr_;

      //! global matrix dLmTLmTLin
      Teuchos::RCP<LINALG::SparseMatrix> dLmTLmTLinMatrixPtr_;

      //! global inactive linearization matrix
      Teuchos::RCP<LINALG::SparseMatrix> inactiveLinMatrixPtr_;

      //! global inactive 2-nd order derivative matrix w.r.t. the displacements
      Teuchos::RCP<LINALG::SparseMatrix> inactive_dd_matrixPtr_;

      //! global inactive diagonal matrix ( represented as vector )
      Teuchos::RCP<Epetra_Vector> inactiveDiagMatrixPtr_;

      //! tributary slave area vector
      Teuchos::RCP<Epetra_Vector> aPtr_;

      //! kappa vector
      Teuchos::RCP<Epetra_Vector> kappaPtr_;

      //! lagrange multiplier vector in normal direction
      Teuchos::RCP<Epetra_Vector> lmNPtr_;

      //! averaged weighted gap vector
      Teuchos::RCP<Epetra_Vector> aWGapPtr_;

      //! weighted gap vector containing values for all slave nodes
      Teuchos::RCP<Epetra_Vector> wGapAllPtr_;

      //! inactive tangential lagrange multiplier right-hand-side
      Teuchos::RCP<Epetra_Vector> dLmTLmTRhsPtr_;

      //! interface force due to Lagrange multiplier values on the slave-side
      Teuchos::RCP<Epetra_Vector> slForceLmPtr_;

      //! interface force due to inactive Lagrange multiplier values on the slave-side
      Teuchos::RCP<Epetra_Vector> slForceLmInactivePtr_;

      //! interface force due to gap values on the slave-side
      Teuchos::RCP<Epetra_Vector> slForceGPtr_;

      //! interface force due to Lagrange multiplier values on the master-side
      Teuchos::RCP<Epetra_Vector> maForceLmPtr_;

      //! interface force due to gap values on the master-side
      Teuchos::RCP<Epetra_Vector> maForceGPtr_;

      //! cn-values of each node
      Teuchos::RCP<Epetra_Vector> cnPtr_;

      //! global slave dof row map in normal direction (for all interfaces)
      Teuchos::RCP<Epetra_Map> gsndofrowmapPtr_;

      //! global slave dof row map in tangential direction (for all interfaces)
      Teuchos::RCP<Epetra_Map> gstdofrowmapPtr_;

      //! global row map of all active slave nodes of the previous Newton step
      Teuchos::RCP<Epetra_Map> gOldActiveSlaveNodesPtr_;

      /// global row map of all slave elements
      Teuchos::RCP<Epetra_Map> gSeleRowMapPtr_;

      /// global column map of all slave elements
      Teuchos::RCP<Epetra_Map> gSeleColMapPtr_;

      /// global row map of all master elements
      Teuchos::RCP<Epetra_Map> gMeleRowMapPtr_;

      /// global column map of all master elements
      Teuchos::RCP<Epetra_Map> gMeleColMapPtr_;

      /// global vector containing the evaluation times of the slave elements
      Teuchos::RCP<Epetra_Vector> gSeleEvalTimesPtr_;

      /// gradient error container
      struct GradError
      {
        GradError() = default;
        ~GradError() = default;

        /// total error
        double total_ = 0.0;

        /// error wrt to the master parametric projection variation
        std::vector<std::pair<int, double>> master_;

        /// error due to the missing jacobian determinant variation
        std::vector<std::pair<int, double>> jacobian_;
      };

      /// see container description
      GradError grad_error_;

      /// @name sub-data container
      /// @{

      /// steepest ascent data container
      Teuchos::RCP<CONTACT::AUG::STEEPESTASCENT::DataContainer> sa_data_ptr_;

      /// @}
    };  // class DataContainer

    /*--------------------------------------------------------------------------*/
    /*! \brief Augmented Lagrangian strategy
     *
     *  \author hiermeier
     *  \date 04/15*/
    class Strategy : public ::CONTACT::CoAbstractStrategy
    {
      /** The combo_strategy is a wrapper class for a set of augmented Lagrangian
       *  strategies and needs access to all methods. */
      friend class CONTACT::AUG::ComboStrategy;

      friend class CONTACT::AUG::ParallelDistributionController;
      friend class CONTACT::AUG::ActiveSet;

     public:
      //! Standard constructor
      Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
          const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
          const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
          const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof);

      //! Destructor
      virtual ~Strategy(){};

      /// derived
      INPAR::CONTACT::SolvingStrategy Type() const override
      {
        return INPAR::CONTACT::solution_augmented;
      }

      /// derived
      bool IsSaddlePointSystem() const override;

      //! reset active set convergence flags
      void ResetActiveSet() override { Data().IsActiveSetConverged() = false; };

      //! Return the L2-norm of the constraint right-hand-side [deprecated]
      double ConstraintNorm() const override;

      //! Save the reference state
      void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override{};

      //! Read the restart information and adjust the members accordingly
      void DoReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<const Epetra_Vector> dis,
          Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr) override;

      /*! \brief Reset the internal state variables [derived]
       *
       *  \date 02/2016
       *  \author hiermeier */
      void Reset(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp,
          const Epetra_Vector& xnew) override;

      //! @name Accessors
      //! @{
      /*! \brief Return convergence status of semi-smooth active set search
       *
       *  The augmented Lagrangian strategy supports only the semi-smooth Newton
       *  check at the moment, so no difference is made! */
      bool ActiveSetSemiSmoothConverged() const override { return Data().IsActiveSetConverged(); }

      //! Return the active node row map of the previous Newton step
      Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const override
      {
        return Data().GOldActiveSlaveNodesPtr();
      };

      //! Return the slip node row map of the previous Newton step
      Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const override
      {
        dserror("No frictional contact for the augmented Lagrangian contact formulation!");
        exit(EXIT_FAILURE);
      };

      /*! \brief Return the slave dof row map in normal direction
       *
       *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
       *                      original map before any redistribution took place.
       *
       *  \date 04/2016
       *  \author hiermeier */
      Teuchos::RCP<const Epetra_Map> SlNormalDoFRowMapPtr(const bool& redist) const override
      {
        if ((not redist) and ParRedist())
          dserror("The original / not redistributed slave normal row map is not available!");

        return Data().GSlNormalDofRowMapPtr();
      };
      const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const override
      {
        return *SlNormalDoFRowMapPtr(redist);
      }

      /*! \brief Return the slave dof row map in the tangential directions
       *
       *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
       *                      original map before any redistribution took place.
       *
       *  \date 04/2016
       *  \author hiermeier */
      Teuchos::RCP<const Epetra_Map> SlTangentialDoFRowMapPtr(const bool& redist) const override
      {
        if ((not redist) and ParRedist())
          dserror("The original / not redistributed slave tangential row map is not available!");

        return Data().GSlTangentialDofRowMapPtr();
      };
      const Epetra_Map& SlTangentialDoFRowMap(const bool& redist) const override
      {
        return *SlTangentialDoFRowMapPtr(redist);
      }

      //! derived
      Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
          const enum DRT::UTILS::VecBlockType& bt) const override;

      //! Return the condensed right hand side vector (currently unsupported!)
      Teuchos::RCP<const Epetra_Vector> GetCondensedRhsPtr(
          Epetra_Vector& f, const double& timefac_np) const override
      {
        dserror("There is no condensed rhs pointer in the augmented Lagrangian case!");
        exit(EXIT_FAILURE);
      };

      //! Return the desired matrix block pointer
      Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(const enum DRT::UTILS::MatBlockType& bt,
          const CONTACT::ParamsInterface* cparams = NULL) const override;

      //! Return the condensed matrix block pointer (currently unsupported!)
      Teuchos::RCP<LINALG::SparseMatrix> GetCondensedMatrixBlockPtr(
          Teuchos::RCP<LINALG::SparseMatrix>& kteff, const double& timefac_np) const override
      {
        dserror("There is no condensed matrix block in the augmented Lagrangian case!");
        exit(EXIT_FAILURE);
      };

      //! Return the contact status of the last nonlinear solver iteration
      bool WasInContactLastIter() const { return Data().WasInContactLastIter(); }

      //! Return augmented constraint rhs vector
      Teuchos::RCP<Epetra_Vector> ConstrRhs() override { return Data().ConstrRhsPtr(); }

      /*! \brief return the weighted gap vector (slave normal dof row map layout)
       *
       *  \param[in] type: Choose the map type. Default is the weighted gap of
       *                   the active nodes only. Alternatively, the active
       *                   weighted gap of all slave nodes can be accessed. In
       *                   the second case non-projectable contributions will
       *                   be filled with a value of 1.0e+12.
       *
       * \author hiermeier \date 10/17 */
      const Epetra_Vector& GetWeightedGap(
          const enum MapType type = MapType::active_slave_nodes) const;

      /*! \brief return the weighted gap gradient
       *
       *  \note If map_type == all_slave_nodes: In general, not all contributions
       *  can be considered, if non-projectable slave elements are involved. Be
       *  aware of possible consequences in your calculations.
       *
       *  \param[in] grad_type  Choose between the gradient of the weighted gap
       *                        used for the force balance and the always
       *                        consistently computed gradient used for the
       *                        constraint linearization.
       *  \param[in] map_type   Specify the range map of the weighted gap gradient.
       *                        If the weighted gap gradient is assembled for
       *                        all slave nodes, non-projectable contributions
       *                        will contain zero values.
       *
       *  \author hiermeier \date 09/17 */
      Teuchos::RCP<const LINALG::SparseMatrix> GetWeightedGapGradient(
          const enum WGapGradientType grad_type, const enum MapType map_type) const;

      void EvalWeightedGapGradientError(CONTACT::ParamsInterface& cparams) override;

      //! @}

      double CharacteristicInterfaceElementLength(const enum CONTACT::AUG::SideType stype) const;

      inline double GetTotalGradientError() const { return Data().TotalGradientError(); }

      const std::vector<std::pair<int, double>>& GetNodalGradientError_Jacobian() const
      {
        return Data().NodalGradientError_Jacobian();
      }

      const std::vector<std::pair<int, double>>& GetNodalGradientError_MaProj() const
      {
        return Data().NodalGradientError_MaProj();
      }

      //! @name Evaluate routines
      //! @{

      //! Evaluate the augmented contact forces on slave and master side
      void AugForces(Epetra_Vector& augfs_lm, Epetra_Vector& augfs_g, Epetra_Vector& augfm_lm,
          Epetra_Vector& augfm_g) const;

      //! @}

      //! @name merit function methods
      //! @{

      /// return the potential function value
      double GetPotentialValue(
          const enum NOX::NLN::MeritFunction::MeritFctName mrt_type) const override;

      /// return the desired contributions to the linear model of the potential function
      double GetLinearizedPotentialValueTerms(const Epetra_Vector& dir,
          const enum NOX::NLN::MeritFunction::MeritFctName mrt_type,
          const enum NOX::NLN::MeritFunction::LinOrder linorder,
          const enum NOX::NLN::MeritFunction::LinType lintype) const override;

      /** \brief Split state vector (e.g. direction or x) into the displacement
       *  part as well as into its active and inactive normal Lagrange multiplier
       *  parts
       *
       *  \note If NULL pointers are used as input arguments, the pointers are
       *  initialized automatically with the correct maps.
       *
       *  \author hiermeier \date 08/17 */
      void SplitStateVector(const Epetra_Vector& full_state,
          Teuchos::RCP<Epetra_Vector>& displ_state_slma_ptr,
          Teuchos::RCP<Epetra_Vector>& z_state_active_ptr,
          Teuchos::RCP<Epetra_Vector>& z_state_inactive_ptr) const;

      /** \brief Split state vector (e.g. direction or x) into the displacement
       *  part as well as into its active and inactive normal Lagrange multiplier
       *  parts
       *
       *  \author hiermeier \date 08/17 */
      void SplitStateVector(const Epetra_Vector& full_state, Epetra_Vector& displ_state_slma,
          Epetra_Vector& z_state_active, Epetra_Vector& z_state_inactive) const;

     protected:
      /// helper function to access the right model terms
      double GetLinearizedPotentialModelTerms(const Epetra_Vector& dir,
          const enum POTENTIAL::Type pottype, const enum POTENTIAL::SetType potset,
          const enum NOX::NLN::MeritFunction::LinOrder linorder,
          const enum NOX::NLN::MeritFunction::LinType lintype) const;

      /// helper function to access the correct 1st order derivative model terms
      void GetLinearizedPotentialModelTerms_1stOrder(const enum POTENTIAL::Type pottype,
          const enum POTENTIAL::SetType potset, const enum NOX::NLN::MeritFunction::LinType lintype,
          double& linval) const;

      // helper function to access the correct 2nd order derivative model terms
      void GetLinearizedPotentialModelTerms_2ndOrder(const enum POTENTIAL::Type pottype,
          const enum POTENTIAL::SetType potset, const enum NOX::NLN::MeritFunction::LinType lintype,
          double& linval) const;
      //! @}

      /// return the scaling factor for the inactive part
      virtual double InactiveScaleFactor() const { return 0.5; }

      /// derived
      void EvaluateReferenceState(Teuchos::RCP<const Epetra_Vector> vec) override;

      /// set current evaluation state in the \c cparams object
      void SetCurrentEvalState(const CONTACT::ParamsInterface& cparams);

      /// perform parallel distribution check
      void CheckParallelDistribution(const GlobalTimeMonitor& global_timer);

      /// handle dynamic redistribution of the contact interface elements
      bool DynRedistributeContact(
          const Teuchos::RCP<const Epetra_Vector>& dis, const int nlniter) override;

      /** spread the global evaluation times of the slave elements to the different
       *  contact interfaces */
      void SpreadGlobalSeleEvalTimesToInterfaces();

      /*! \brief Update routine called at the end of a load/time step [derived]
       *
       *  \param[in] dis  Converged displacemnt state.
       *
       *  \author hiermeier \date 06/17 */
      void Update(Teuchos::RCP<const Epetra_Vector> dis) override;

      //! derived
      plain_interface_set& Interfaces() override;

      //! derived
      const plain_interface_set& Interfaces() const override;

      /*! \brief  Assemble the global sets gsndofrowmap_ and gstdofrowmap_
       *
       *   These maps contain all slave dofs in normal and tangential direction,
       *   assembled over all interfaces. */
      void AssembleGlobalSlNTDofRowMaps();

      /// assemble global element maps
      void AssembleGlobalEleMaps();

      /*! \brief Setup the cn-vector
       *
       *  We allow theoretically different cn values for each node in the
       *  augmented Lagrangian  case. */
      void InitializeCn(const double cn_init);

      /// redistribute the cn vector
      void RedistributeCn();

      /*! \brief Initialize and evaluate augmented Mortar stuff for the next Newton step
       *
       *  This method first checks if we are dealing with self contact and updates
       *  the interface slave and master sets if necessary. Then it resets the global
       *  Mortar matrices Dn and Mn.
       *
       *  The nodal quantities computed in InitEvalInterface() are then assembled
       *  to global matrices. No setup of the global system is to be done here yet,
       *  so there is no need to pass in the effective stiffness K or the effective
       *  load vector f. (-->Evaluate routine) */
      void InitMortar() override;

      //! Assemble the mortar matrices
      void AssembleMortar() override;

      /*! \brief Split of the Dn and Mn matrices into an active and inactive part
       *
       *  This split becomes necessary for the least squares update routine and can
       *  be prevented in all other cases, if the evaluation order is changed (1st:
       *  update of the active set, 2nd: evaluation of Dn/Mn). Unfortunately, it's not
       *  possible to keep the complete Dn/Mn matrices (slave map), because the inactive
       *  Lagrange multipliers are not reduced to zero in one Newton step (because of
       *  the consistent linearization of the area, see also the hint in the EvaluateContact
       *  method). A direct consequence is, that we have to distinguish in a hard way
       *  between inactive and active quantities in contrast to the standard Lagrangian case,
       *  where the inactive quantities vanish in the force balance, due to the one-step
       *  reduction of the inactive Lagrange multipliers.
       *
       *  \remark The split is quite expensive and should be used carefully (Complete calls).
       *
       *  \author hiermeier */
      void SplitMortar();

      /*! \brief Update augmented active set and check for convergence
       *
       *  In this function we loop over all interfaces and then over all
       *  slave nodes to check, whether the assumption of them being active
       *  or inactive respectively has been correct. If a single node changes
       *  state, the active set is adapted accordingly and the convegence
       *  flag is kept on false. Here we have the semi-smooth Newton case
       *  with one combined iteration loop for active set search and large
       *  deformations. As a consequence this method is called AFTER each
       *  (not yet converged) Newton step. If there is a change in the active
       *  set or the residual and disp norm are still above their limits,
       *  another Newton step has to be performed.
       *
       *  \author hiermeier */
      void UpdateActiveSetSemiSmooth(const CONTACT::ParamsInterface& cparams);

      //! Initialize all matrices
      void Initialize() override { Initialize(MORTAR::eval_force_stiff); }
      //! Initialize only the necessary member variables
      void Initialize(enum MORTAR::ActionType actiontype);

      /*! \brief Projection of the nodal LM values in the nodal normal and tangential direction.
       *
       *  Pay attention, that the Lagrange multipliers
       *  don't have any directions, they are scalar values. This function is just
       *  used to visualize the LM values.
       *  A better and more realistic visualization is the corresponding contact force.
       *
       *  Function overwrites the basic function in the abstract lagrange strategy. */
      void OutputStresses() override;

      /*! \brief Write contact force output
       *
       *  \author hiermeier \date 12/17 */
      void WriteOutput(IO::DiscretizationWriter& writer) const override;

      /*! @name Auxiliary routines, debugging and visualization methods
       *        (all these methods are defined in the contact_augmented_strategy_tools.cpp) */
      //! @{

      //! Check linear and angular momentum conservation
      void CheckConservationLaws(CONTACT::ParamsInterface& cparams);

      //! \brief adapt the cn-value (optional)
      /** The adaptive cn routine is in addition to the default regularization
       *  parameter correction. Think of an meaningful initial choice of the
       *  regularization parameter for example.
       *
       * \param(in) cparams  Contact parameter interface
       *
       *  \author hiermeier \date 03/18 */
      virtual void AdaptiveCn(CONTACT::ParamsInterface& cparams);

      //! Finite Difference check of the augmented Lagrange terms at GP level
      void AugFDCheckGP(CONTACT::ParamsInterface& cparams);

      //! Finite Difference check of the augmented Lagrange terms at global level
      void AugFDCheckGlobal(CONTACT::ParamsInterface& cparams);

      /*! \brief Evaluate the global FD check w.r.t. the displacements
       *
       *  This function is only suitable for the augmented Lagrange formulation due to
       *  the special structure of the EvaluateContact implementation! */
      class FD_Debug : public DRT::SingletonDestruction
      {
       public:
        /// singleton creation
        static FD_Debug* Instance(
            Strategy* strat, const double delta = 1.0e-8, const bool delete_me = false);

        /// delete singleton pointer
        virtual void Done();

        /// perform the FD evaluate calls
        void Evaluate(const Teuchos::RCP<LINALG::SparseMatrix>& derivMatrixPtr,
            Teuchos::RCP<Epetra_Vector>& rhsVector, CONTACT::ParamsInterface& cparams);

       private:
        /// constructor
        FD_Debug() : strat_(NULL), delta_(0.0), is_fd_check_(false)
        { /*empty*/
        }

        /// destructor
        ~FD_Debug(){};

        /// initialize this object
        void Init(Strategy* strat, const double delta = 1.0e-8);

        /// apply the perturbation
        void DoPerturbation(const int gid, const int dof);

        /// undo the perturbation
        void UndoPerturbation(const int gid, const int dof) const;

        DRT::Node* FindINode(const int gid) const;

       private:
        /// call back to the wrapping strategy
        Strategy* strat_;

        /// finite difference delta
        double delta_;

        /// avoid recursive calls
        bool is_fd_check_;

        /// reference values to undo a previous perturbation
        std::map<int, LINALG::Matrix<3, 1>> ref_x_;
      };
      //@}

      //! @name (Derived) internal evaluate routines
      //! @{
      /*! \brief Redistribute and setup augmented Lagrangian members
       *
       *  \param redistributed (in) : currently a redistribution is taking place
       *  \param init          (in) : initialization is running
       *
       *  \author hiermeier \date 03/2016 */
      void PostSetup(bool redistributed, bool init) override;

      /*! \brief Run at the beginning of a call to EvalForce
       *
       *  Prepare the evaluation and assembly and integrate all necessary quantities.
       *
       *  \author hiermeier \date 03/17 */
      void PreEvalForce(CONTACT::ParamsInterface& cparams);

      /*! \brief Compute force terms
       *
       *  \param cparams (in): parameter interface between the contact objects and
       *                       the structural time integration
       *
       *  \date 03/2016
       *  \author hiermeier */
      void EvalForce(CONTACT::ParamsInterface& cparams) override;

      /*! \brief Run in the end of a call to EvalForce
       *
       *  Post-process force terms.
       *
       *  \author hiermeier \date 04/17 */
      void PostEvalForce(CONTACT::ParamsInterface& cparams);

      /// evaluate augmented forces w.r.t. Lagrange multiplier and weighted gap
      void EvalAugmentedForces();

      /// evaluate only the forces originating from the weight gap values
      void EvalConstraintForces();

      /// set all lm forces to zero
      void ZeroizeLMForces();

      /*! \brief Compute force and stiffness terms
       *
       * \param cparams (in): parameter interface between the contact objects and the structural
       * time integration
       *
       *  \date 03/2016
       *  \author hiermeier */
      void EvalForceStiff(CONTACT::ParamsInterface& cparams) override;

      /// evaluate only all contributions corresponding to the weighted gap
      void EvalStaticConstraintRHS(CONTACT::ParamsInterface& cparams) override;

      /*! \brief Run in the end of a call to EvalForceStiff
       *
       *  \author hiermeier \date 03/17 */
      void PostEvalForceStiff(CONTACT::ParamsInterface& cparams);

      /*! \brief Recover contact specific solution variables
       *
       * \param cparams (in): parameter interface between the contact objects and the structural
       * time integration \param xold    (in): old solution vector of the NOX solver \param dir
       * (in): current search direction (in general NOT the actual step, keep in mind that the step
       * length can differ from 1.0) \param xnew    (in): new solution vector of the NOX solver
       *
       * (see the CONTACT::CoAbstractStrategy for more information)
       *
       * \date 03/2016
       * \author hiermeier */
      void RunPostComputeX(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      /*! \brief Reset the internal stored Lagrange multipliers
       *
       *  \param cparams (in): parameter interface between the contact objects and the structural
       * time integration \param xnew    (in): new solution vector of the NOX solver
       *
       *  \date 07/2016
       *  \author hiermeier */
      void ResetLagrangeMultipliers(
          const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew) override;

      //! \brief Derived function
      inline void InitEvalInterface(CONTACT::ParamsInterface& cparams)
      {
        Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::rcpFromRef(cparams);
        InitEvalInterface(cparams_ptr);
      }
      void InitEvalInterface(Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr) override;

      //! evaluate interface
      void EvalInterface(CONTACT::AUG::Interface& interface, const int rriter,
          const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

      //! The contributions to the structural right-hand-side block are calculated.
      virtual void EvalStrContactRHS();

      //! The entries of the constraint right-hand side are calculated.
      void EvalConstrRHS() override;

      //! All necessary contributions are added to the constraint right-hand side.
      virtual void AddContributionsToConstrRHS(Epetra_Vector& augConstrRhs) const;

      //! @}

      //! Assemble contact contributions to the rhs (frictionless)
      void AssembleGap();

      //! Assemble contact contributions to the rhs (frictionless)
      bool AssembleContactRHS();

      //! Evaluate contact contributions to stiff matrix (frictionless)
      void AssembleContactStiff();

      /*! \brief Create new stiffness state
       *
       *  \author hiermeier \date 03/17 */
      void CreateStiffnessState(const Epetra_Map& gAugInactiveSlaveDofs);

      /*! \brief Create new rhs state
       *
       *  \author hiermeier \date 03/17 */
      void CreateRhsState(const Epetra_Map& gAugInactiveSlaveDofs);

      /*! \brief Zeroize existing stiffness state
       *
       *  Zeroize all quantities which are directly related to the evaluation
       *  of the stiffness contributions. This covers matrices as well as
       *  vectors representing diagonal matrices.
       *
       *  \author hiermeier \date 03/17 */
      void ZeroizeStiffnessState();

      /*! \brief Zeroize existing rhs state
       *
       *  \author hiermeier \date 03/17 */
      void ZeroizeRhsState();

      /** \brief All contributions to the final system matrix block:
       *  ROW => Displacement, COLUMN => Displacement
       *
       *  \author hiermeier \date 03/17 */
      virtual void AddContributionsToMatrixBlockDisplDispl(
          LINALG::SparseMatrix& kdd, const CONTACT::ParamsInterface* cparams = NULL) const;

      /** \brief All contributions to the final system matrix block:
       *  ROW => Displacement, COLUMN => Lagrange Multiplier
       *
       *  \author hiermeier \date 03/17 */
      void AddContributionsToMatrixBlockDisplLm(LINALG::SparseMatrix& kdz) const;

      /** \brief All contributions to the final system matrix block:
       *  ROW => Lagrange Multiplier, COLUMN => Displacement
       *
       *  \author hiermeier \date 03/17 */
      void AddContributionsToMatrixBlockLmDispl(LINALG::SparseMatrix& kzd) const;

      /** \brief All contributions to the final system matrix block:
       *  ROW => Lagrange Multiplier, COLUMN => Lagrange Multiplier
       *
       *  \author hiermeier \date 03/17 */
      virtual void AddContributionsToMatrixBlockLmLm(LINALG::SparseMatrix& kzz) const;

     protected:
      /*! \brief Get access to the internal data container of the strategy (mutable)
       *
       * \author hiermeier
       * \date 05/16 */
      inline CONTACT::AUG::DataContainer& Data() { return augData_; }

      /*! \brief Get access to the internal data container of the strategy (read-only)
       *
       * \author hiermeier
       * \date 05/16 */
      inline const CONTACT::AUG::DataContainer& Data() const { return augData_; }

     public:
      //! @name Unsupported derived routines (dead-end)
      //! @{
      //! @name No support for nested active set loops
      //! @{
      bool ActiveSetConverged() override { return true; };
      int ActiveSetSteps() override { return -1; };
      void UpdateActiveSet() override
      {
        dserror("No support for fixed nested active set strategy!");
      };
      //! @}

      //! @name Deprecated methods
      //! @{
      void EvaluateContact(
          Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) override
      {
        dserror("Deprecated function call!");
      };
      void EvaluateFriction(
          Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) override
      {
        dserror("Deprecated function call!");
      };
      void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
          Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
          Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,
          Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
          Teuchos::RCP<Epetra_Vector>& blockrhs) override
      {
        dserror("Deprecated function call!");
      };
      void UpdateDisplacementsAndLMincrements(
          Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol) override
      {
        dserror("Deprecated function call!");
      };
      void Recover(Teuchos::RCP<Epetra_Vector> disi) override
      {
        dserror("Deprecated function call! Replaced by RunPostComputeX().");
      };
      void UpdateActiveSetSemiSmooth() override
      {
        dserror(
            "Deprecated function call! Replaced by "
            "UpdateActiveSetSemiSmooth( const CONTACT::ParamsInterface& ).");
      };
      //! @}

      //! @name No frictional support at the moment
      //! @{
      void EvaluateRelMovPredict() override
      {
        if (Data().IsFriction()) dserror("No frictional contact support at the moment!");
        return;
      };
      //! @}

      /*! @name Dead-end for penalty and Uzawa methods (wrong strategy)
       *
       * Please note, that the definition of these functions seems completely unnecessary here.
       * Actually it would be a much better idea to cast the object to the right strategy at the
       * place where it is needed.                                            hiermeier 05/16 */
      //! @{
      double InitialPenalty() override
      {
        dserror("Wrong strategy!");
        exit(EXIT_FAILURE);
      };
      void InitializeUzawa(
          Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) override
      {
        dserror("Wrong strategy!");
      };
      void ResetPenalty() override { dserror("Wrong strategy!"); };
      void ModifyPenalty() override { dserror("Wrong strategy!"); };
      void UpdateUzawaAugmentedLagrange() override { dserror("Wrong strategy!"); };
      void UpdateConstraintNorm(int uzawaiter = 0) override { dserror("Wrong strategy!"); };
      //! @}
      //! @}


     protected:
      // Don't want = operator and cctor
      Strategy operator=(const Strategy& old) = delete;
      Strategy(const Strategy& old) = delete;

     private:
      /// augmented contact strategy data container pointer
      Teuchos::RCP<CONTACT::AUG::DataContainer> augDataPtr_;

      /// reference to the augmented contact strategy data container
      CONTACT::AUG::DataContainer& augData_;

      /// contact interface set
      plain_interface_set interface_;

    };  // class Strategy

    /// compute the l2-error of the values in the \c error_map
    void L2ErrorNormPerNode(const std::unordered_map<int, Deriv1stMap>& error_map,
        std::vector<std::pair<int, double>>& error_norm_per_node);

    /// return the sum of all squared errors in \c error_map_vec
    double MyTotalSquareError(
        const std::unordered_map<int, Deriv1stMap>* const* error_map_vec, const unsigned num_vecs);

    /*--------------------------------------------------------------------------*/
    /** \brief Extract a sub-sparse matrix from a source sparse matrix
     *
     *  \param[in] source         Given sparse matrix
     *  \param[in] target_row_map row/range map of the target matrix
     *  \param[in] target_col_map column/domain map of the target matrix
     *
     *  \return A RCP to the completed sub-sparse matrix will be returned.
     *
     *  \author hiermeier \date 06/17 */
    Teuchos::RCP<LINALG::SparseMatrix> ExtractMatrix(const LINALG::SparseMatrix& source,
        const Epetra_Map& target_row_map, const Epetra_Map& target_col_map);

    /*--------------------------------------------------------------------------*/
    /*! \brief Scale the vector entries of the target vector with the entries of
     *  the source vector
     *
     *  The length of the source2targetMap has to be as long as the length
     *  of the target vector map (on each processor). The boolean <inverse> gives
     *  the possibility to scale the target vector with the reciprocal values of the
     *  source vector entries.
     *
     *  The routine does exactly the same as
     *  \f[
     *                    target = diag(source) \cdot target.
     *  \f]
     *
     *  \author hiermeier */
    void MultiplyElementwise(const Epetra_Vector& source, const Epetra_Map& source2targetMap,
        Epetra_Vector& target, const bool inverse);
    inline void MultiplyElementwise(
        const Epetra_Vector& source, const Epetra_Map& source2targetMap, Epetra_Vector& target)
    {
      MultiplyElementwise(source, source2targetMap, target, false);
    };

    /*--------------------------------------------------------------------------*/
    /*! \brief Redistribute a row map by using a reference map
     *
     *  The reference map is used as blue-print for the redistribution of the
     *  red_map entries.
     *  The reference map must have the right distribution and the %red_map
     *  has to be a subset of the reference map. In this way we can look up
     *  the right processor in the reference map and the redistribution is
     *  achieved easily.
     *
     *  \author hiermeier \date 03/17 */
    void RedistributeRowMap(const Epetra_Map& ref_map, Epetra_Map& red_map);

  }  // namespace AUG
}  // namespace CONTACT


#endif /* #ifndef CONTACT_AUGMENTED_STRATEGY_H_ */
