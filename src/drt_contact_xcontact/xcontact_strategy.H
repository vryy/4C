/*---------------------------------------------------------------------------*/
/*!
\file xcontact_strategy.H

\brief Extended contact solving strategy with Level-Set and X-FEM

\level 3

\maintainer Michael Hiermeier

\date Aug 18, 2016
*/
/*---------------------------------------------------------------------------*/

#ifndef XCONTACT_STRATEGY_H
#define XCONTACT_STRATEGY_H

#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_xfem/xfem_enum_lists.H"

namespace XCONTACT
{
  // Forward declaration
  class Interface;

  class DataContainer : public CONTACT::AbstractStratDataContainer
  {
   public:
    //! @name Constructors and destructors and related methods
    //! @{

    //! Standard constructor
    DataContainer();

    //! Destructor
    virtual ~DataContainer() {}

    //! @}


    //! @name Accessors
    //! @{

    //! @name Boolean accessors
    //! @{
    //! TODO
    //! @}

    //! @name Vector accessors
    //! @{

    //! Return strContactRhsPtr_
    Teuchos::RCP<Epetra_Vector>& StrContactRhsPtr() { return strContactRhsPtr_; }
    Teuchos::RCP<const Epetra_Vector> StrContactRhsPtr() const { return strContactRhsPtr_; }
    Epetra_Vector& StrContactRhs()
    {
      if (strContactRhsPtr_.is_null())
      {
        dserror("XContact data container: Pointer to vector strContactRhsPtr_ not initialized!");
      }
      return *strContactRhsPtr_;
    }

    //! Return lmN_ptr_
    Teuchos::RCP<Epetra_Vector>& LmNPtr() { return lmN_ptr_; }
    Teuchos::RCP<const Epetra_Vector> LmNPtr() const { return lmN_ptr_; }
    Epetra_Vector& LmN()
    {
      if (lmN_ptr_.is_null())
      {
        dserror("XContact data container: Pointer to vector lmN_ptr_ not initialized!");
      }
      return *lmN_ptr_;
    }

    //! Return Wc_lm_
    Teuchos::RCP<Epetra_Vector>& WcLmPtr() { return Wc_lm_ptr_; }
    Teuchos::RCP<const Epetra_Vector> WcLmPtr() const { return Wc_lm_ptr_; }
    Epetra_Vector& WcLm()
    {
      if (Wc_lm_ptr_.is_null())
      {
        dserror("XContact data container: Pointer to vector Wc_lm_ not initialized.");
      }
      return *Wc_lm_ptr_;
    }

    //! @}

    //! @name Matrix accessors
    //! @{

    //! Return Wc_su_u_
    Teuchos::RCP<LINALG::SparseMatrix>& WcSuUPtr() { return Wc_su_u_ptr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> WcSuUPtr() const { return Wc_su_u_ptr_; }
    LINALG::SparseMatrix& WcSuU()
    {
      if (Wc_su_u_ptr_.is_null())
      {
        dserror("XContact data container: Pointer to matrix Wc_su_u_ not initialized.");
      }
      return *Wc_su_u_ptr_;
    }

    //! Return Wc_mu_u_
    Teuchos::RCP<LINALG::SparseMatrix>& WcMuUPtr() { return Wc_mu_u_ptr_; }
    Teuchos::RCP<const LINALG::SparseMatrix> WcMuUPtr() const { return Wc_mu_u_ptr_; }
    LINALG::SparseMatrix& WcMuU()
    {
      if (Wc_mu_u_ptr_.is_null())
      {
        dserror("XContact data container: Pointer to matrix Wc_mu_u_ not initialized.");
      }
      return *Wc_mu_u_ptr_;
    }

    //! @}

    //! @name Global map accessors
    //! @{

    //! Return gsndofrowmapPtr_
    Teuchos::RCP<Epetra_Map>& GSlNormalDofRowMapPtr() { return gsndofrowmapPtr_; }
    Teuchos::RCP<const Epetra_Map> GSlNormalDofRowMapPtr() const { return gsndofrowmapPtr_; }
    Epetra_Map& GSlNormalDofRowMap()
    {
      if (gsndofrowmapPtr_.is_null())
      {
        dserror("XContact data container: Pointer to map gsndofrowmapPtr_ not initialized.");
      }
      return *gsndofrowmapPtr_;
    }

    //! Return gstdofrowmapPtr_
    Teuchos::RCP<Epetra_Map>& GSlTangentialDofRowMapPtr() { return gstdofrowmapPtr_; }
    Teuchos::RCP<const Epetra_Map> GSlTangentialDofRowMapPtr() const { return gstdofrowmapPtr_; }
    Epetra_Map& GSlTangentialDofRowMap()
    {
      if (gstdofrowmapPtr_.is_null())
      {
        dserror("XContact data container: Pointer to map gstdofrowmapPtr_ not initialized.");
      }
      return *gstdofrowmapPtr_;
    }

    //! @}

    //! @}


   protected:
    //! @name Member variables
    //! @{

    //! @name Global vectors
    //! @{

    //! Pointer to structural contact residuum vector
    Teuchos::RCP<Epetra_Vector> strContactRhsPtr_;

    //! Pointer to Lagrange multiplier vector in normal direction
    Teuchos::RCP<Epetra_Vector> lmN_ptr_;

    //! Pointer to constraint residuum vector
    Teuchos::RCP<Epetra_Vector> Wc_lm_ptr_;

    //! @}


    //! @name Global matrices
    //! @{

    //! Pointer to linearization of slave varied contact potential
    Teuchos::RCP<LINALG::SparseMatrix> Wc_su_u_ptr_;

    //! Pointer to linearization of master varied contact potential
    Teuchos::RCP<LINALG::SparseMatrix> Wc_mu_u_ptr_;

    //! @}


    //! @name Global maps
    //! @{

    //! Pointer to global slave DOF row map in normal direction (for all interfaces)
    Teuchos::RCP<Epetra_Map> gsndofrowmapPtr_;

    //! Pointer to global slave DOF row map in tangential direction (for all interfaces)
    Teuchos::RCP<Epetra_Map> gstdofrowmapPtr_;

    //! @}

    //! @}

  };  // class DataContainer

  /*-------------------------------------------------------------------------*/

  class Strategy : public CONTACT::CoAbstractStrategy
  {
   public:
    //! @name Constructors and destructors and related methods
    //! @{

    //! Standard constructor
    Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap,
        const Teuchos::ParameterList& params,
        const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& interfaces, const int& dim,
        const Teuchos::RCP<const Epetra_Comm>& comm, const int& maxdof);

    //! Destructor
    virtual ~Strategy() {}

    //! @}


    //! @name Definition of purely virtual functions
    //! @{

    //! Compute L2-norm of constraint right-hand-side
    virtual double ConstraintNorm() const;

    //! Save reference state
    virtual void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis){};

    //! @}

    //! @name Set methods
    //! @{

    //! set contact status (from partitioned algorithm)
    void SetContactStatus(const bool& is_in_contact);

    //! @}

    //! @name Accessors
    //! @{

    const GEN::pairedvector<XFEM::FieldName, DRT::Discretization*>& GetIDiscrets() const;

    //! return the weighted gap vector (slave normal dof row map layout)
    const Epetra_Vector& GetWeightedGap() const;

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
        const enum DRT::UTILS::VecBlockType& bt) const;

    //! [derived]
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMatrixBlockPtr(
        const enum DRT::UTILS::MatBlockType& bt,
        const CONTACT::ParamsInterface* cparams = NULL) const;

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Map> SlNormalDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist())  // TODO: ! and &&
      {
        dserror("The original / not redistributed slave normal row map is not available!");
      }
      return Data().GSlNormalDofRowMapPtr();
    }
    virtual const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const
    {
      return *SlNormalDoFRowMapPtr(redist);
    }

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Map> SlTangentialDoFRowMapPtr(const bool& redist) const
    {
      if ((not redist) and ParRedist())
      {
        dserror("The original / not redistributed slave tangential row map is not available!");
      }
      return Data().GSlTangentialDofRowMapPtr();
    }
    virtual const Epetra_Map& SlTangentialDoFRowMap(const bool& redist) const
    {
      return *SlTangentialDoFRowMapPtr(redist);
    }

    //! @}

    //! @name Evaluate routines
    //! @{
    //! TODO
    //! @}

    //! @name Preconditioner methods
    //! @{

    //! [derived]
    virtual bool IsSaddlePointSystem() const;

    //! @}


   protected:
    // TODO: Combine purely virtual functions and derived internal evaluate routines?

    //! @name Definition of purely virtual functions
    //! @{

    //! derived
    virtual std::vector<Teuchos::RCP<CONTACT::CoInterface>>& Interfaces();

    //! derived
    virtual const std::vector<Teuchos::RCP<CONTACT::CoInterface>>& Interfaces() const;

    //! Initialize all matrices
    virtual void Initialize() { Initialize(MORTAR::eval_force_stiff); }

    //! Initialize only the necessary member variables
    void Initialize(MORTAR::ActionType actiontype);

    //! Evaluate constraint right-hand side (constraint residuum)
    virtual void EvalConstrRHS();

    //! Update active set (TODO: Just for testing.)
    // virtual void UpdateActiveSetSemiSmooth();

    //! @}


    //! @name Derived internal evaluate routines
    //! @{

    //! [derived]
    virtual void EvalWeightedGap(CONTACT::ParamsInterface& cparams);

    //! [derived]
    virtual void EvalForceStiff(CONTACT::ParamsInterface& cparams);

    //! [derived]
    virtual void EvalForce(CONTACT::ParamsInterface& cparams);

    //! [derived]
    virtual void InitEvalInterface(Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

    //! The contributions to the structural right-hand-side block are calculated.
    virtual void EvalStrContactRHS();

    //! [derived]
    virtual void RunPostComputeX(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
        const Epetra_Vector& dir, const Epetra_Vector& xnew);

    //! [derived]
    virtual void ResetLagrangeMultipliers(
        const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew);

    //! @}


    //! @name Internal evaluate routines
    //! @{

    //! Assemble global maps gsndofrowmap_ and gstdofrowmap_
    void AssembleGlobalSlNTDofRowMaps();

    //! Assemble mortar matrices
    virtual void AssembleMortar();

    //! Assemble the nodal weighted gap vector
    void AssembleWeightedGap();

    //! Assemble contact contributions to the rhs (frictionless)
    virtual void AssembleContactRHS();

    // Assemble contact stiffness matrix (contact tangent matrix)
    virtual void AssembleContactStiff();

    //! @}


    //! @name Auxiliary routines, debugging and visualization methods
    //! @{
    //! TODO
    //! @}


   private:
    //! store a paired vector with pointers to the interface discretizations
    void StoreIDiscrets();

    //! Get access to internal strategy data container (mutable)
    DataContainer& Data()
    {
      if (xDataPtr_.is_null())
      {
        dserror("The XContact strategy data container is not initialized!");
      }
      return *xDataPtr_;
    }

    //! Get access to internal strategy data container (read-only)
    const DataContainer& Data() const
    {
      if (xDataPtr_.is_null())
      {
        dserror("The XContact strategy data container is not initialized!");
      }
      return *xDataPtr_;
    }

   private:
    //! @name Member variables
    //! @{

    //! Data container
    Teuchos::RCP<DataContainer> xDataPtr_;

    //! Interfaces
    std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface_;

    //! pointer to the interface discretizations
    GEN::pairedvector<XFEM::FieldName, DRT::Discretization*> idiscrets_;

    //! @}

   public:
    //! @name Unsupported purely virtual functions
    //! @{

    //! @name Friction not supported at the moment
    //! @{
    virtual void EvaluateRelMovPredict()
    {
      if (Data().IsFriction())
      {
        dserror("No frictional contact support at the moment!");
      }
      return;
    };
    //! @}
    //! @name Active set strategy not used
    //! @{
    virtual void UpdateActiveSet() {}
    virtual void ResetActiveSet() {}
    virtual void UpdateActiveSetSemiSmooth() {}
    virtual int ActiveSetSteps() { return -1; }
    virtual bool ActiveSetConverged() { return true; }
    virtual bool ActiveSetSemiSmoothConverged() const { return true; }
    virtual Teuchos::RCP<const Epetra_Map> GetOldActiveRowNodes() const { return Teuchos::null; }
    virtual Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const { return Teuchos::null; }
    //! @}

    //! @name Penalty and Uzawa methods (wrong strategy) (TODO: See note in
    //! AugmentedLagrangeStrategy)
    //! @{
    virtual double InitialPenalty()
    {
      dserror("Wrong strategy!");
      exit(EXIT_FAILURE);
    }
    virtual void InitializeUzawa(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
    {
      dserror("Wrong strategy!");
    }
    virtual void ResetPenalty() { dserror("Wrong strategy!"); }
    virtual void ModifyPenalty() { dserror("Wrong strategy!"); }
    virtual void UpdateUzawaAugmentedLagrange() { dserror("Wrong strategy!"); }
    virtual void UpdateConstraintNorm(int uzawaiter = 0) { dserror("Wrong strategy!"); }
    //! @}

    //! @name Deprecated methods
    //! @{
    virtual void EvaluateContact(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
    {
      dserror("Deprecated function call!");
    }
    virtual void EvaluateFriction(
        Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff)
    {
      dserror("Deprecated function call!");
    }
    virtual void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd,
        Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
        Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,
        Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol,
        Teuchos::RCP<Epetra_Vector>& blockrhs)
    {
      dserror("Deprecated function call!");
    }
    virtual void UpdateDisplacementsAndLMincrements(
        Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol)
    {
      dserror("Deprecated function call!");
    }
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi)
    {
      dserror("Deprecated function call! Replaced by RunPostComputeX().");
    }
    //! @}

    //! @}

  };  // class Strategy
}  // namespace XCONTACT

#endif  // #ifndef XCONTACT_STRATEGY_H
