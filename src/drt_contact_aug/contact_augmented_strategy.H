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

namespace CONTACT
{
// forward declarartion
class AugmentedInterface;

class AugStratDataContainer : public AbstractStratDataContainer
{
public:
  //! @name Constructors and destructors and related methods
  //! @{
  //! Standard constructor
  AugStratDataContainer();

  //! Destructor
  virtual ~AugStratDataContainer() { return; }
  //! @}

  //! @name Access methods
  //! @{

  //! @name Booleans accessors
  //! @{

  //! Return gFdCheck_ indicator
  bool& FiniteDifferenceIndicator() { return gFdCheck_; }
  const bool& FiniteDifferenceIndicator() const { return gFdCheck_; }

  //! Return WasInContactLastIter indicator
  bool& WasInContactLastIter() { return wasincontactlastiter_; }
  const bool& WasInContactLastIter() const { return wasincontactlastiter_; }

  //! Return TRUE if the active set is converged, otherwise false
  bool& IsActiveSetConverged() { return isactivesetconverged_; };
  const bool& IsActiveSetConverged() const { return isactivesetconverged_; }

  //! print the linear conservation check results if TRUE
  bool& PrintLinearMomConservation() { return printlinearconservation_; };
  const bool& PrintLinearMomConservation() const { return printlinearconservation_; };

  //! print the angular conservation check results if TRUE
  bool& PrintAngularMomConservation() { return printangularconservation_; };
  const bool& PrintAngularMomConservation() const { return printangularconservation_; };

  //! @}

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

  //!
  Teuchos::RCP<LINALG::SparseMatrix>& InactiveMatrixPtr() { return inactiveMatrixPtr_; }
  Teuchos::RCP<const LINALG::SparseMatrix> InactiveMatrixPtr() const { return inactiveMatrixPtr_; }
  LINALG::SparseMatrix& InactiveMatrix()
  {
    if (inactiveMatrixPtr_.is_null())
      dserror("The augInactiveMatrixPtr_ is not initialized!");
    return *inactiveMatrixPtr_;
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

  //! Returns the nodal c_n vector
  Teuchos::RCP<Epetra_Vector>& CnPtr() { return cnPtr_; }
  Teuchos::RCP<const Epetra_Vector> CnPtr() const { return cnPtr_; }
  Epetra_Vector& Cn()
  {
    if (cnPtr_.is_null())
      dserror("The cnPtr_ is not initialized!");
    return *cnPtr_;
  }

  //! Returns the nodal c_n update vector
  Teuchos::RCP<Epetra_Vector>& UCnPtr() { return uCnPtr_; }
  Teuchos::RCP<const Epetra_Vector> UCnPtr() const { return uCnPtr_; }
  Epetra_Vector& UCn()
  {
    if (uCnPtr_.is_null())
      dserror("The uCnPtr_ is not initialized!");
    return *uCnPtr_;
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
  AugStratDataContainer operator =(const AugStratDataContainer& old);
  AugStratDataContainer(const AugStratDataContainer& old);

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

  //! global inactive matrix
  Teuchos::RCP<LINALG::SparseMatrix> inactiveMatrixPtr_;

  //! global inactive linearization matrix
  Teuchos::RCP<LINALG::SparseMatrix> inactiveLinMatrixPtr_;

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

  //! cn-values of each node
  Teuchos::RCP<Epetra_Vector> cnPtr_;

  //! Vector which is used to decide if the cn-value of the node should be updated or not.
  Teuchos::RCP<Epetra_Vector> uCnPtr_;

  //! global slave dof row map in normal direction (for all interfaces)
  Teuchos::RCP<Epetra_Map> gsndofrowmapPtr_;

  //! global slave dof row map in tangential direction (for all interfaces)
  Teuchos::RCP<Epetra_Map> gstdofrowmapPtr_;

  //! global map of all active slave nodes of the previous Newton step
  Teuchos::RCP<Epetra_Map> gOldActiveSlaveNodesPtr_;
};

/*! \brief AUGMENTED LAGRANGIAN STRATEGY for a full step non-linear solver method
 *
 *  \remark If you are interested in using more sophisticated non-linear solver
 *  strategies please see the CONTACT::NLN::Augmented, CONTACT::NLN::SteepestAscent
 *  and CONTACT::NLN::Lagrange classes.
 *
 *  \author hiermeier
 *  \date 04/15*/
class AugmentedLagrangeStrategy : public CoAbstractStrategy
{
public:
  //! Standard constructor
  AugmentedLagrangeStrategy(
      const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
      const Epetra_Map* DofRowMap,
      const Epetra_Map* NodeRowMap,
      const Teuchos::ParameterList& params,
      std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
      const int& dim,
      const Teuchos::RCP<const Epetra_Comm>& comm,
      const int& maxdof);

  //! Destructor
  virtual ~AugmentedLagrangeStrategy() {};

  //! reset active set convergence flags
  virtual void ResetActiveSet()
  { Data().IsActiveSetConverged() = false; };

  //! Return the L2-norm of the constraint right-hand-side
  virtual double ConstraintNorm() const;

  //! Save the reference state
  virtual void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis)
  { };

  //! @name Accessors
  //! @{
  /*! \brief Return convergence status of semi-smooth active set search
   *
   *  The augmented Lagrangian strategy supports only the semi-smooth Newton
   *  check at the moment, so no difference is made! */
  bool ActiveSetSemiSmoothConverged() const
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

  //! Dummy function which is necessary in a derived nln-strategy class
  virtual bool IsPreAsymptotic() const {return false; }
  //! @}

  //! @name Evaluate routines
  //! @{

  //! Evaluate the augmented contact forces on slave and master side
  void AugForces(Epetra_Vector& augfs_lm, Epetra_Vector& augfs_g,
                 Epetra_Vector& augfm_lm, Epetra_Vector& augfm_g) const;

  //! @}

protected:
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
  void Initialize();

  /*! \brief Projection of the nodal LM values in the nodal normal and tangential direction.
   *
   *  Pay attention, that the Lagrange multipliers
   *  don't have any directions, they are scalar values. This function is just
   *  used to visualize the LM values.
   *  A better and more realistic visualization is the corresponding contact force.
   *
   *  Function overwrites the basic function in the abstract lagrange strategy. */
  void OutputStresses();

  /*! @name Auxiliary routines, debugging and visualization methods
   *        (all these methods are defined in the contact_augmented_strategy_tools.cpp) */
  //! @{

  /*! \brief Scale the vector entries of the target vector with the entries of the source vector
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
                             Epetra_Vector& target) const
  { MultiplyElementwise(source,source2targetMap,target,false); };
  void MultiplyElementwise(const Epetra_Vector& source,
                           const Epetra_Map& source2targetMap,
                           Epetra_Vector& target,
                           const bool inverse) const;

  /*! \brief Redistribute a row map by using a reference map
   *
   *  The reference map must have the right distribution and the modMap
   *  has to be a subset of the reference map. In this way we can look up
   *  the right processor in the reference map and the redistribution is easily
   *  performed. */
  void RedistributeRowMap(const Epetra_Map& refMap, Epetra_Map& modMap);

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
  /*! \brief Redistribute augmented Lagrangian stuff ... */
  virtual void PostSetup(const bool& redistributed, const bool& init);

  /*! \brief Compute force terms
   *
   *  \param cparams (in): parameter interface between the contact objects and the structural time integration
   *
   *  \date 03/2016
   *  \author hiermeier */
  virtual void EvalForce(CONTACT::ParamsInterface& cparams);

  /*! \brief Compute force and stiffness terms
   *
   * \param cparams (in): parameter interface between the contact objects and the structural time integration
   *
   *  \date 03/2016
   *  \author hiermeier */
  virtual void EvalForceStiff(CONTACT::ParamsInterface& cparams);

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

  //! @}

  //! Assemble contact contributions to the rhs (frictionless)
  virtual void AssembleContactRHS();

  //! Evaluate contact contributions to stiff matrix (frictionless)
  virtual void AssembleContactStiff();

  //! @name Functions which are implemented in some of the derived classes
  //! @{

  //! Assemble the global slave/master dof row map without DBC DoFs
  virtual void AssembleSlMaNoDbcDofRowMap(Teuchos::RCP<LINALG::MapExtractor> dbcmaps) {};
  //! @}
private:
  /*! \brief Get access to the internal data container of the strategy (mutable)
   *
   * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
   *
   * \author hiermeier
   * \date 05/16 */
  CONTACT::AugStratDataContainer& Data()
  {
    if (augDataPtr_.is_null())
      dserror("The augmented strategy data container is not initialized!");
    return *augDataPtr_;
  }

  /*! \brief Get access to the internal data container of the strategy (read-only)
   *
   * \remark This has to stay PRIVATE, otherwise this function becomes ambiguous.
   *
   * \author hiermeier
   * \date 05/16 */
  const CONTACT::AugStratDataContainer& Data() const
  {
    if (augDataPtr_.is_null())
      dserror("The augmented strategy data container is not initialized!");
    return *augDataPtr_;
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
  AugmentedLagrangeStrategy operator = (const AugmentedLagrangeStrategy& old);
  AugmentedLagrangeStrategy(const AugmentedLagrangeStrategy& old);

private:
  Teuchos::RCP<CONTACT::AugStratDataContainer> augDataPtr_;
  std::vector<Teuchos::RCP<CONTACT::AugmentedInterface> > interface_;
};
}


#endif /* #ifndef CONTACT_AUGMENTED_STRATEGY_H_ */
