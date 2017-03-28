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

// forward declaration
namespace CONTACT {
namespace AUG {
namespace STEEPESTASCENT {
  class DataContainer;
} // namespace STEEPESTASCENT
  class Interface;
  class Potential;
  class ComboStrategy;

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
    void InitSubDataContainer( const INPAR::CONTACT::SolvingStrategy strat_type);

    //! @}

    //! @name Access methods
    //! @{

    //! @name Sub-data container accessors
    //! @{

    CONTACT::AUG::STEEPESTASCENT::DataContainer& SaData()
    {
      if ( sa_data_ptr_.is_null() )
        dserror( "The steepest ascent sub data container was not initialized!" );
      return *sa_data_ptr_;
    }

    const CONTACT::AUG::STEEPESTASCENT::DataContainer& SaData() const
    {
      if ( sa_data_ptr_.is_null() )
        dserror( "The steepest ascent sub data container was not initialized!" );
      return *sa_data_ptr_;
    }

    //! @}

    //! @name Booleans accessors
    //! @{

    //! Return gFdCheck_ indicator
    bool& FiniteDifferenceIndicator() { return gFdCheck_; }
    bool FiniteDifferenceIndicator() const { return gFdCheck_; }

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

    //! set the the semi smooth newton flag
    void SetIsSemiSmoothNewton( bool is_semi_smooth )
    {
      is_semi_smooth_newton_ = is_semi_smooth;
    }

    //! get the semi smooth newton flag
    bool IsSemiSmoothNewton() const
    {
      return is_semi_smooth_newton_;
    }

    //! Set the matrix maps status
    void SetMatrixMapsValid( bool isvalid )
    {
      matrix_maps_valid_ = isvalid;
    }

    //! Are the matrix maps valid or did they change?
    bool MatrixMapsValid() const
    {
      return matrix_maps_valid_;
    }

    //! Set the vector maps status
    void SetVectorMapsValid( bool isvalid )
    {
      vector_maps_valid_ = isvalid;
    }

    //! Are the vector maps valid or did they change?
    bool VectorMapsValid() const
    {
      return vector_maps_valid_;
    }

    //! @}

    //! set constant semi smooth cn (coming from the Input file)
    void SetConstantCn( double cn )
    {
      cn_ = cn;
    }

    //! access constant semi smooth cn (coming from the Input file)
    double ConstantdCn() const
    {
      if ( cn_ < 0.0 )
        dserror( "cn was not set correctly! (cn_ = %d)", cn_ );
      return cn_;
    };

    //! get the parallel strategy enumerator
    void SetParallelStrategy( INPAR::MORTAR::ParallelStrategy parallel_strat )
    {
      parallel_strategy_ = parallel_strat;
    }

    //! get the parallel strategy enumerator
    INPAR::MORTAR::ParallelStrategy ParallelStrategy() const
    {
      return parallel_strategy_;
    }

    void SetPotential( const Teuchos::RCP<Potential>& potential )
    {
      potentialPtr_ = potential;
    }

    CONTACT::AUG::Potential& Potential()
    {
      return *potentialPtr_;
    }

    const CONTACT::AUG::Potential& Potential() const
    {
      return *potentialPtr_;
    }

    //! @name Matrix accessors
    //! @{

    //! Return dGLmSlLinMatrix
    Teuchos::RCP<LINALG::SparseMatrix>& DGLmSlLinMatrixPtr() { return dGLmSlLinMatrixPtr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> DGLmSlLinMatrixPtr() const { return dGLmSlLinMatrixPtr_; }
    LINALG::SparseMatrix& DGLmSlLinMatrix()
    {
      if (dGLmSlLinMatrixPtr_.is_null())
        dserror("The dGLmSlLinMatrixPtr_ is not initialized!");
      return *dGLmSlLinMatrixPtr_;
    }

    //! Return dGLmMaLinMatrix
    Teuchos::RCP<LINALG::SparseMatrix>& DGLmMaLinMatrixPtr() { return dGLmMaLinMatrixPtr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> DGLmMaLinMatrixPtr() const { return dGLmMaLinMatrixPtr_; }
    LINALG::SparseMatrix& DGLmMaLinMatrix()
    {
      if (dGLmMaLinMatrixPtr_.is_null())
        dserror("The dGLmMaLinMatrixPtr_ is not initialized!");
      return *dGLmMaLinMatrixPtr_;
    }

    //! Return dGGSlLinMatrix
    Teuchos::RCP<LINALG::SparseMatrix>& DGGSlLinMatrixPtr() { return dGGSlLinMatrixPtr_;}
    Teuchos::RCP<const LINALG::SparseMatrix> DGGSlLinMatrixPtr() const { return dGGSlLinMatrixPtr_;}
    LINALG::SparseMatrix& DGGSlLinMatrix()
    {
      if (dGGSlLinMatrixPtr_.is_null())
        dserror("The dGGSlLinMatrixPtr_ is not initialized!");
      return*dGGSlLinMatrixPtr_;
    }

    //! Return dGGMaLinMatrix
    Teuchos::RCP<LINALG::SparseMatrix>& DGGMaLinMatrixPtr() { return dGGMaLinMatrixPtr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> DGGMaLinMatrixPtr() const { return dGGMaLinMatrixPtr_; }
    LINALG::SparseMatrix& DGGMaLinMatrix()
    {
      if (dGGMaLinMatrixPtr_.is_null())
        dserror("The dGGMaLinMatrixPtr_ is not initialized!");
      return *dGGMaLinMatrixPtr_;
    }

    //! Return dLmNWGapLinMatrix
    Teuchos::RCP<LINALG::SparseMatrix>& DLmNWGapLinMatrixPtr() { return dLmNWGapLinMatrixPtr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> DLmNWGapLinMatrixPtr() const { return dLmNWGapLinMatrixPtr_; }
    LINALG::SparseMatrix& DLmNWGapLinMatrix()
    {
      if (dLmNWGapLinMatrixPtr_.is_null())
        dserror("The dLmNWGapLinMatrixPtr_ is not initialized!");
      return *dLmNWGapLinMatrixPtr_;
    }

    //! Return dLmTLmTMatrix
    Teuchos::RCP<LINALG::SparseMatrix>& DLmTLmTMatrixPtr() { return dLmTLmTMatrixPtr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> DLmTLmTMatrixPtr() const { return dLmTLmTMatrixPtr_; }
    LINALG::SparseMatrix& DLmTLmTMatrix()
    {
      if (dLmTLmTMatrixPtr_.is_null())
        dserror("The dLmTLmTMatrixPtr_ is not initialized!");
      return *dLmTLmTMatrixPtr_;
    }

    //! Return dLmTLmTLinMatrix
    Teuchos::RCP<LINALG::SparseMatrix>& DLmTLmTLinMatrixPtr() { return dLmTLmTLinMatrixPtr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> DLmTLmTLinMatrixPtr() const { return dLmTLmTLinMatrixPtr_; }
    LINALG::SparseMatrix& DLmTLmTLinMatrix()
    {
      if (dLmTLmTLinMatrixPtr_.is_null())
        dserror("The dLmTLmTLinMatrixPtr_ is not initialized!");
      return*dLmTLmTLinMatrixPtr_;
    }

    Teuchos::RCP<LINALG::SparseMatrix>& InactiveLinMatrixPtr() { return inactiveLinMatrixPtr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> InactiveLinMatrixPtr() const { return inactiveLinMatrixPtr_; }
    LINALG::SparseMatrix& InactiveLinMatrix()
    {
      if (inactiveLinMatrixPtr_.is_null())
        dserror("The augInactiveLinMatrixPtr_ is not initialized!");
      return *inactiveLinMatrixPtr_;
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
      if (aPtr_.is_null())
        dserror("The augAPtr_ is not initialized!");
      return *aPtr_;
    }

    //! This vector contributes as diagonal inactive matrix to the kzz block
    Teuchos::RCP<Epetra_Vector>& InactiveDiagMatrixPtr() { return inactiveDiagMatrixPtr_; }
    Teuchos::RCP<const Epetra_Vector> InactiveDiagMatrixPtr() const { return inactiveDiagMatrixPtr_; }
    Epetra_Vector& InactiveDiagMatrix()
    {
      if (inactiveDiagMatrixPtr_.is_null())
        dserror("The inactiveDiagMatrixPtr_ is not initialized!");
      return *inactiveDiagMatrixPtr_;
    }

    /*! \brief Return a vector of the tributary (weighted) nodal areas
     *
     *  These quantities are integrated over the projectable slave side! */
    Teuchos::RCP<Epetra_Vector>& KappaVecPtr() {return kappaPtr_; }
    Teuchos::RCP<const Epetra_Vector> KappaVecPtr() const {return kappaPtr_; }
    Epetra_Vector& KappaVec()
    {
      if (kappaPtr_.is_null())
        dserror("The kappaPtr_ is not initialized!");
      return *kappaPtr_;
    }

    //! Return a vector of the lagrange multipliers for the frictionless/normal constraint
    Teuchos::RCP<Epetra_Vector>& LmNPtr() {return lmNPtr_; }
    Teuchos::RCP<const Epetra_Vector> LmNPtr() const {return lmNPtr_; }
    Epetra_Vector& LmN()
    {
      if (lmNPtr_.is_null())
        dserror("The lmNPtr_ is not initialized!");
      return *lmNPtr_;
    }

    //! Return the averaged weighted gap
    Teuchos::RCP<Epetra_Vector>& AWGapPtr() { return aWGapPtr_; }
    Teuchos::RCP<const Epetra_Vector> AWGapPtr() const { return aWGapPtr_; }
    Epetra_Vector& AWGap()
    {
      if (aWGapPtr_.is_null())
        dserror("The aWGapRhsPtr_ is not initialized!");
      return *aWGapPtr_;
    }

    /*! \brief Return the tangential Lagrange multiplier right hand side
     *
     *  In the frictionless case this is supposed to do nothing and just keeps
     *  the tangential Lagrange multipliers at a constant value of zero. */
    Teuchos::RCP<Epetra_Vector>& DLmTLmTRhsPtr() {return dLmTLmTRhsPtr_; }
    Teuchos::RCP<const Epetra_Vector> DLmTLmTRhsPtr() const {return dLmTLmTRhsPtr_; }
    Epetra_Vector& DLmTLmTRhs()
    {
      if (dLmTLmTRhsPtr_.is_null())
        dserror("The dLmTLmTRhsPtr_ is not initialized!");
      return *dLmTLmTRhsPtr_;
    }

    /*! \brief Return the slave force due to Lagrange multiplier values */
    Teuchos::RCP<Epetra_Vector>& SlForceLmPtr() { return slForceLmPtr_; }
    Teuchos::RCP<const Epetra_Vector> SlForceLmPtr() const { return slForceLmPtr_.getConst(); }
    Epetra_Vector& SlForceLm()
    {
      if ( slForceLmPtr_.is_null() )
        dserror( "The slForceLmPtr_ is not initialized!" );
      return *slForceLmPtr_;
    }

    /*! \brief Return the slave force due to gap values */
    Teuchos::RCP<Epetra_Vector>& SlForceGPtr() { return slForceGPtr_; }
    Teuchos::RCP<const Epetra_Vector> SlForceGPtr() const { return slForceGPtr_.getConst(); }
    Epetra_Vector& SlForceG()
    {
      if ( slForceGPtr_.is_null() )
        dserror( "The slForceGPtr_ is not initialized!" );
      return *slForceGPtr_;
    }

    /*! \brief Return the master force due to Lagrange multiplier values */
    Teuchos::RCP<Epetra_Vector>& MaForceLmPtr() { return maForceLmPtr_; }
    Teuchos::RCP<const Epetra_Vector> MaForceLmPtr() const { return maForceLmPtr_.getConst(); }
    Epetra_Vector& MaForceLm()
    {
      if ( maForceLmPtr_.is_null() )
        dserror( "The maForceLmPtr_ is not initialized!" );
      return *maForceLmPtr_;
    }

    /*! \brief Return the master force due to gap values */
    Teuchos::RCP<Epetra_Vector>& MaForceGPtr() { return maForceGPtr_; }
    Teuchos::RCP<const Epetra_Vector> MaForceGPtr() const { return maForceGPtr_.getConst(); }
    Epetra_Vector& MaForceG()
    {
      if ( maForceGPtr_.is_null() )
        dserror( "The maForceGPtr_ is not initialized!" );
      return *maForceGPtr_;
    }

    //! Returns the nodal c_n vector
    Teuchos::RCP<Epetra_Vector>& CnPtr() { return cnPtr_; }
    Teuchos::RCP<const Epetra_Vector> CnPtr() const { return cnPtr_; }
    Epetra_Vector& Cn()
    {
      if (cnPtr_.is_null())
        dserror("The cnPtr_ is not initialized!");
      return *cnPtr_;
    }
    const Epetra_Vector& Cn() const
    {
      if (cnPtr_.is_null())
        dserror("The cnPtr_ is not initialized!");
      return *cnPtr_;
    }

    //! @}

    //! @name Global map accessors
    //! @{
    Teuchos::RCP<Epetra_Map>& GSlNormalDofRowMapPtr() { return gsndofrowmapPtr_; }
    Teuchos::RCP<const Epetra_Map> GSlNormalDofRowMapPtr() const { return gsndofrowmapPtr_; }
    Epetra_Map& GSlNormalDofRowMap()
    {
      if (gsndofrowmapPtr_.is_null())
        dserror("The gsndofrowmapPtr_ is not initialized!");
      return *gsndofrowmapPtr_;
    }

    Teuchos::RCP<Epetra_Map>& GSlTangentialDofRowMapPtr() { return gstdofrowmapPtr_; }
    Teuchos::RCP<const Epetra_Map> GSlTangentialDofRowMapPtr() const { return gstdofrowmapPtr_; }
    Epetra_Map& GSlTangentialDofRowMap()
    {
      if (gstdofrowmapPtr_.is_null())
        dserror("The gstdofrowmapPtr_ is not initialized!");
      return *gstdofrowmapPtr_;
    }

    Teuchos::RCP<Epetra_Map>& GOldActiveSlaveNodesPtr() { return gOldActiveSlaveNodesPtr_; }
    Teuchos::RCP<const Epetra_Map> GOldActiveSlaveNodesPtr() const { return gOldActiveSlaveNodesPtr_; }
    Epetra_Map& GOldActiveSlaveNodes()
    {
      if (gOldActiveSlaveNodesPtr_.is_null())
        dserror("The gAugOldActiveSlaveNodesPtr_ is not initialized!");
      return *gOldActiveSlaveNodesPtr_;
    }
    //! @}
    //! @}

  protected:
    // don't want = operator and cctor
    DataContainer operator =(const DataContainer& old);
    DataContainer(const DataContainer& old);

    //! indicates if the global finite difference check is active
    bool gFdCheck_;

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

    //! flag indicating if the slave and master matrix maps changed ( e.g. due to redistribution )
    bool matrix_maps_valid_;

    //! flag indicating if the slave and master vector maps changed ( e.g. due to redistribution )
    bool vector_maps_valid_;

    //! constant cn value from the input file
    double cn_;

    enum INPAR::MORTAR::ParallelStrategy parallel_strategy_;

    Teuchos::RCP<CONTACT::AUG::Potential> potentialPtr_;

    //! global matrix dGLmSlLin
    Teuchos::RCP<LINALG::SparseMatrix> dGLmSlLinMatrixPtr_;

    //! global matrix dGLmMaLin
    Teuchos::RCP<LINALG::SparseMatrix> dGLmMaLinMatrixPtr_;

    //! global matrix dGGSlLin
    Teuchos::RCP<LINALG::SparseMatrix> dGGSlLinMatrixPtr_;

    //! global matrix dGGMaLin
    Teuchos::RCP<LINALG::SparseMatrix> dGGMaLinMatrixPtr_;

    //! global matrix dLmNWGapLin
    Teuchos::RCP<LINALG::SparseMatrix> dLmNWGapLinMatrixPtr_;

    //! global matrix dLmTLmT
    Teuchos::RCP<LINALG::SparseMatrix> dLmTLmTMatrixPtr_;

    //! global matrix dLmTLmTLin
    Teuchos::RCP<LINALG::SparseMatrix> dLmTLmTLinMatrixPtr_;

    //! global inactive linearization matrix
    Teuchos::RCP<LINALG::SparseMatrix> inactiveLinMatrixPtr_;

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

    //! inactive tangential lagrange multiplier right-hand-side
    Teuchos::RCP<Epetra_Vector> dLmTLmTRhsPtr_;

    //! interface force due to Lagrange multiplier values on the slave-side
    Teuchos::RCP<Epetra_Vector> slForceLmPtr_;

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

    //! global map of all active slave nodes of the previous Newton step
    Teuchos::RCP<Epetra_Map> gOldActiveSlaveNodesPtr_;

    /// @name sub-data container
    /// @{

    Teuchos::RCP<CONTACT::AUG::STEEPESTASCENT::DataContainer> sa_data_ptr_;

    /// @}
  };  // class DataContainer

  /*--------------------------------------------------------------------------*/
  /*! \brief AUGMENTED LAGRANGIAN STRATEGY for a full step non-linear solver method
   *
   *  \remark If you are interested in using more sophisticated non-linear solver
   *  strategies please see the CONTACT::NLN::Augmented, CONTACT::NLN::SteepestAscent
   *  and CONTACT::NLN::Lagrange classes.
   *
   *  \author hiermeier
   *  \date 04/15*/
  class Strategy : public ::CONTACT::CoAbstractStrategy
  {
    /** The combo_strategy is a wrapper class for a set of augmented Lagrangian
     *  strategies and needs access to all methods. */
    friend class CONTACT::AUG::ComboStrategy;

  public:
    //! Standard constructor
    Strategy(
        const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        const std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
        int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm,
        int maxdof);

    //! Destructor
    virtual ~Strategy() {};

    virtual INPAR::CONTACT::SolvingStrategy Type() const
    {
      return INPAR::CONTACT::solution_augmented;
    }

    virtual bool IsSaddlePointSystem() const;

    //! reset active set convergence flags
    virtual void ResetActiveSet()
    { Data().IsActiveSetConverged() = false; };

    //! Return the L2-norm of the constraint right-hand-side [deprecated]
    virtual double ConstraintNorm() const;

    //! Save the reference state
    virtual void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis)
    { };

    //! Read the restart information and adjust the members accordingly
    virtual void DoReadRestart(
        IO::DiscretizationReader& reader,
        Teuchos::RCP<const Epetra_Vector> dis,
        Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr );

    //! @name Accessors
    //! @{
    /*! \brief Return convergence status of semi-smooth active set search
     *
     *  The augmented Lagrangian strategy supports only the semi-smooth Newton
     *  check at the moment, so no difference is made! */
    virtual bool ActiveSetSemiSmoothConverged() const
    {
      return Data().IsActiveSetConverged();
    }

    //! Return the active node row map of the previous Newton step
    virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const
    { return Data().GOldActiveSlaveNodesPtr(); };

    //! Return the slip node row map of the previous Newton step
    virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const
    { dserror("No frictional contact for the augmented Lagrangian contact formulation!");
      exit(EXIT_FAILURE); };

    /*! \brief Return the slave dof row map in normal direction
     *
     *  \param redist (in): If TRUE, the redistributed map is returned, otherwise the
     *                      original map before any redistribution took place.
     *
     *  \date 04/2016
     *  \author hiermeier */
    virtual Teuchos::RCP<const Epetra_Map> SlNormalDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist())
        dserror("The original / not redistributed slave normal row map is not available!");

      return Data().GSlNormalDofRowMapPtr();
    };
    virtual const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const
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
    virtual Teuchos::RCP<const Epetra_Map> SlTangentialDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist())
        dserror("The original / not redistributed slave tangential row map is not available!");

      return Data().GSlTangentialDofRowMapPtr();
    };
    virtual const Epetra_Map& SlTangentialDoFRowMap(const bool& redist) const
    {
      return *SlTangentialDoFRowMapPtr(redist);
    }

    //! derived
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum DRT::UTILS::VecBlockType& bt) const;

    //! Return the condensed right hand side vector (currently unsupported!)
    virtual Teuchos::RCP<const Epetra_Vector> GetCondensedRhsPtr() const
    { dserror("There is no condensed rhs pointer in the augmented Lagrangian case!");
      exit(EXIT_FAILURE); };

    //! Return the desired matrix block pointer
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt) const;

    //! Return the condensed matrix block pointer (currently unsupported!)
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetCondensedMatrixBlockPtr() const
    { dserror("There is no condensed matrix block in the augmented Lagrangian case!");
      exit(EXIT_FAILURE); };

    //! Return the contact status of the last nonlinear solver iteration
    bool WasInContactLastIter() const { return Data().WasInContactLastIter(); }

    //! Return augmented constraint rhs vector
    virtual Teuchos::RCP<Epetra_Vector> ConstrRhs() { return Data().ConstrRhsPtr(); }

    //! @}

    //! @name Evaluate routines
    //! @{

    //! Evaluate the augmented contact forces on slave and master side
    void AugForces(Epetra_Vector& augfs_lm, Epetra_Vector& augfs_g,
                   Epetra_Vector& augfm_lm, Epetra_Vector& augfm_g) const;

    //! @}


  protected:
    //! derived
    virtual plain_interface_set& Interfaces();

    //! derived
    virtual const plain_interface_set& Interfaces() const;

    /*! \brief  Assemble the global sets gsndofrowmap_ and gstdofrowmap_
     *
     *   These maps contain all slave dofs in normal and tangential direction, over all
     *   interfaces, respectively. */
    void AssembleGlobalSlNTDofRowMaps();

    /*! \brief Setup the cn-vector
     *
     *  We allow different cn values for each node in the augmented Lagrangian case. */
    void InitializeCn(bool redistributed, bool init);

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
    virtual void InitMortar();

    //! Assemble the mortar matrices
    virtual void AssembleMortar();

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
    virtual void UpdateActiveSetSemiSmooth()
    { UpdateActiveSetSemiSmooth(false); };
    void UpdateActiveSetSemiSmooth(const bool& correction);

    //! Initialize all matrices
    virtual void Initialize()
    {
      Initialize( MORTAR::eval_force_stiff );
    }
    //! Initialize only the necessary member variables
    void Initialize( enum MORTAR::ActionType actiontype );

    /*! \brief Projection of the nodal LM values in the nodal normal and tangential direction.
     *
     *  Pay attention, that the Lagrange multipliers
     *  don't have any directions, they are scalar values. This function is just
     *  used to visualize the LM values.
     *  A better and more realistic visualization is the corresponding contact force.
     *
     *  Function overwrites the basic function in the abstract lagrange strategy. */
    virtual void OutputStresses();

    /*! @name Auxiliary routines, debugging and visualization methods
     *        (all these methods are defined in the contact_augmented_strategy_tools.cpp) */
    //! @{

     //! Check linear and angular momentum conservation
    void CheckConservationLaws(const Epetra_Vector& augfs, const Epetra_Vector& augfm);

    //! Finite Difference check of the augmented Lagrange terms at GP level
    void AugFDCheckGP(CONTACT::ParamsInterface& cparams);

    //! Finite Difference check of the augmented Lagrange terms at global level
    void AugFDCheckGlobal(CONTACT::ParamsInterface& cparams);

    /*! \brief Evaluate the global FD check w.r.t. the displacements
     *
     *  This function is only suitable for the augmented Lagrange formulation due to
     *  the special structure of the EvaluateContact implementation! */
    void EvalFDCheckGlobalDispl(
        LINALG::SparseMatrix derivMatrix,
        Teuchos::RCP<Epetra_Vector>& rhsVector,
        CONTACT::ParamsInterface& cparams);
    //@}

    //! @name (Derived) internal evaluate routines
    //! @{
    /*! \brief Redistribute and setup augmented Lagrangian members
     *
     *  \param redistributed (in) : currently a redistribution is taking place
     *  \param init          (in) : initialization is running
     *
     *  \author hiermeier \date 03/2016 */
    virtual void PostSetup( bool redistributed, bool init );

    /*! \brief Compute force terms
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void EvalForce(CONTACT::ParamsInterface& cparams);

    void EvalAugmentedForces();

    /*! \brief Run at the beginning of a call to EvalForce
     *
     *  Prepare the evaluation and assembly and integrate all necessary quantities.
     *
     *  \author hiermeier \date 03/17 */
    void PreEvalForce( CONTACT::ParamsInterface& cparams );

    /*! \brief Compute force and stiffness terms
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time integration
     *
     *  \date 03/2016
     *  \author hiermeier */
    virtual void EvalForceStiff(CONTACT::ParamsInterface& cparams);

    /*! \brief Run in the end of a call to EvalForceStiff
     *
     *  \author hiermeier \date 03/17 */
    void PostEvalForceStiff( CONTACT::ParamsInterface& cparams );

    /*! \brief Recover contact specific solution variables
     *
     * \param cparams (in): parameter interface between the contact objects and the structural time integration
     * \param xold    (in): old solution vector of the NOX solver
     * \param dir     (in): current search direction (in general NOT the actual step, keep in mind that the step length can differ from 1.0)
     * \param xnew    (in): new solution vector of the NOX solver
     *
     * (see the CONTACT::CoAbstractStrategy for more information)
     *
     * \date 03/2016
     * \author hiermeier */
    virtual void RecoverState(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew);

    /*! \brief Reset the internal stored Lagrange multipliers
     *
     *  \param cparams (in): parameter interface between the contact objects and the structural time integration
     *  \param xnew    (in): new solution vector of the NOX solver
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual void ResetLagrangeMultipliers(
        const CONTACT::ParamsInterface& cparams,
        const Epetra_Vector& xnew);

    //! \brief Derived function
    virtual void InitEvalInterface(Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    //! The contributions to the structural right-hand-side block are calculated.
    virtual void EvalStrContactRHS();

    //! The entries of the constraint right-hand side are calculated.
    virtual void EvalConstrRHS();

    //! All necessary contributions are added to the constraint right-hand side.
    virtual void AddContributionsToConstrRHS( Epetra_Vector& augConstrRhs ) const;

    //! @}

    //! Assemble contact contributions to the rhs (frictionless)
    void AssembleGap();

    //! Assemble contact contributions to the rhs (frictionless)
    virtual bool AssembleContactRHS();

    //! Evaluate contact contributions to stiff matrix (frictionless)
    virtual void AssembleContactStiff();

    /*! \brief Create new state matrices
     *
     *  \author hiermeier \date 03/17 */
    void CreateMatrices( const Epetra_Map & gAugInactiveSlaveDofs );

    /*! \brief Create new vectors
     *
     *  \author hiermeier \date 03/17 */
    void CreateVectors( const Epetra_Map & gAugInactiveSlaveDofs );

    /*! \brief Zeroize existing state matrices
     *
     *  \author hiermeier \date 03/17 */
    void ZeroizeMatrices();

    /*! \brief Zeroize existing state vectors
     *
     *  \author hiermeier \date 03/17 */
    void ZeroizeVectors();

    /** \brief All contributions to the final system matrix block:
     *  ROW => Displacement, COLUMN => Displacement
     *
     *  \author hiermeier \date 03/17 */
    virtual void AddContributionsToMatrixBlockDisplDispl( LINALG::SparseMatrix& kdd ) const;

    /** \brief All contributions to the final system matrix block:
     *  ROW => Displacement, COLUMN => Lagrange Multiplier
     *
     *  \author hiermeier \date 03/17 */
    virtual void AddContributionsToMatrixBlockDisplLm( LINALG::SparseMatrix& kdz ) const;

    /** \brief All contributions to the final system matrix block:
     *  ROW => Lagrange Multiplier, COLUMN => Displacement
     *
     *  \author hiermeier \date 03/17 */
    virtual void AddContributionsToMatrixBlockLmDispl( LINALG::SparseMatrix& kzd ) const;

    /** \brief All contributions to the final system matrix block:
     *  ROW => Lagrange Multiplier, COLUMN => Lagrange Multiplier
     *
     *  \author hiermeier \date 03/17 */
    virtual void AddContributionsToMatrixBlockLmLm( LINALG::SparseMatrix& kzz ) const;

  protected:
    /*! \brief Get access to the internal data container of the strategy (mutable)
     *
     * \author hiermeier
     * \date 05/16 */
    inline CONTACT::AUG::DataContainer& Data()
    {
      return augData_;
    }

    /*! \brief Get access to the internal data container of the strategy (read-only)
     *
     * \author hiermeier
     * \date 05/16 */
    inline const CONTACT::AUG::DataContainer& Data() const
    {
      return augData_;
    }

  public:
    //! @name Unsupported derived routines (dead-end)
    //! @{
    //! @name No support for nested active set loops
    //! @{
    virtual bool ActiveSetConverged() {return true; };
    virtual int ActiveSetSteps() {return -1; };
    virtual void UpdateActiveSet()
    { dserror("No support for fixed nested active set strategy!"); };
    //! @}

    //! @name Deprecated methods
    //! @{
    virtual void EvaluateContact(
        Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff)
    { dserror("Deprecated function call!"); };
    virtual void EvaluateFriction(
        Teuchos::RCP<LINALG::SparseOperator>& kteff,
        Teuchos::RCP<Epetra_Vector>& feff)
    { dserror("Deprecated function call!"); };
    virtual void BuildSaddlePointSystem(
        Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd,
        Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps,
        int numiter,Teuchos::RCP<Epetra_Operator>& blockMat,
        Teuchos::RCP<Epetra_Vector>& blocksol,
        Teuchos::RCP<Epetra_Vector>& blockrhs)
    { dserror("Deprecated function call!"); };
    virtual void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<Epetra_Vector> blocksol)
    { dserror("Deprecated function call!"); };
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi)
    { dserror("Deprecated function call! Replaced by RecoverState()."); };
    //! @}

    //! @name No frictional support at the moment
    //! @{
    virtual void EvaluateRelMovPredict()
    {
      if (Data().IsFriction())
        dserror("No frictional contact support at the moment!");
      return;
    };
    //! @}

    /*! @name Dead-end for penalty and Uzawa methods (wrong strategy)
     *
     * Please note, that the definition of these functions seems completely unnecessary here.
     * Actually it would be a much better idea to cast the object to the right strategy at the
     * place where it is needed.                                            hiermeier 05/16 */
    //! @{
    virtual double InitialPenalty()
    { dserror("Wrong strategy!"); exit(EXIT_FAILURE); };
    virtual void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
    { dserror("Wrong strategy!"); };
    virtual void ResetPenalty()
    { dserror("Wrong strategy!"); };
    virtual void ModifyPenalty()
    { dserror("Wrong strategy!"); };
    virtual void UpdateUzawaAugmentedLagrange()
    { dserror("Wrong strategy!"); };
    virtual void UpdateConstraintNorm(int uzawaiter = 0)
    { dserror("Wrong strategy!"); };
    //! @}
    //! @}


  protected:
    // Don't want = operator and cctor
    Strategy operator = (const Strategy& old);
    Strategy(const Strategy& old);

  private:
    /// augmented contact strategy data container pointer
    Teuchos::RCP<CONTACT::AUG::DataContainer> augDataPtr_;

    /// reference to the augmented contact strategy data container
    CONTACT::AUG::DataContainer& augData_;

    /// contact interface set
    plain_interface_set interface_;

  };  // class Strategy

  /*--------------------------------------------------------------------------*/
  /*! \brief Scale the vector entries of the target vector with the entries of
   *  the source vector
   *
   *  The length of the source2tragetMap has to be as long as the length
   *  of the target vector map (on each processor). The boolean <inverse> gives
   *  the possibility to scale the target vector with the reciprocal values of the
   *  source vector entries.
   *
   *  The routine does exactly the same as
   *  \f[
   *                    target = diag(source) * target.
   *  \f]
   *
   *  \author hiermeier */
  void MultiplyElementwise(const Epetra_Vector& source,
                           const Epetra_Map& source2targetMap,
                           Epetra_Vector& target,
                           const bool inverse);
  inline void MultiplyElementwise(const Epetra_Vector& source,
                             const Epetra_Map& source2targetMap,
                             Epetra_Vector& target)
  {
    MultiplyElementwise(source,source2targetMap,target,false);
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
   void RedistributeRowMap(const Epetra_Map& ref_map,
       Epetra_Map& red_map);

} // namespace AUG
} // namespace CONTACT


#endif /* #ifndef CONTACT_AUGMENTED_STRATEGY_H_ */
