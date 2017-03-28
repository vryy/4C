/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_interface.H

\brief Augmented contact interface.

\level 2

\maintainer Michael Hiermeier

\date Apr 16, 2014

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_AUGMENTED_INTERFACE_H_
#define CONTACT_AUGMENTED_INTERFACE_H_

#include "../drt_contact/contact_interface.H"

namespace CONTACT {
  // forward declaration
  class ParamsInterface;

namespace AUG {

  /*----------------------------------------------------------------------------*/
  /** \brief Augmented contact interface data container
   *
   *  This class is supposed to contain all relevant members for the augmented
   *  contact interfaces. The external storage in this object, instead of the
   *  actual interface class itself, makes it possible to share the interface
   *  data between different interface objects w/o the need of copying them.
   *
   *  \author hiermeier \date 03/17 */
  class IDataContainer : public ::CONTACT::IDataContainer
  {
  public:
    /// constructor
    IDataContainer();

    /// destructor
    virtual ~IDataContainer() {};

    /// @name Accessors
    /// @{

    inline double& PenBound()
    {
      return penBound_;
    }

    inline double PenBound() const
    {
      return penBound_;
    }

    inline double& Ct()
    {
      return ct_;
    }

    inline double Ct() const
    {
      return ct_;
    }

    inline int& MaxNumMasterElements()
    {
      return maxNumMasterElements_;
    }

    inline int MaxNumMasterElements() const
    {
      return maxNumMasterElements_;
    }

    inline bool& IsSetup()
    {
      return issetup_;
    }

    inline bool IsSetup() const
    {
      return issetup_;
    }

    inline Teuchos::RCP<Epetra_Map>& SNDofRowMap()
    {
      return sndofrowmap_;
    }

    inline Teuchos::RCP<const Epetra_Map> SNDofRowMap() const
    {
      return sndofrowmap_;
    }

    inline Teuchos::RCP<Epetra_Map>& STDofRowMap()
    {
      return stdofrowmap_;
    }

    inline Teuchos::RCP<const Epetra_Map> STDofRowMap() const
    {
      return stdofrowmap_;
    }

    /// @}

  private:
    //! interface penetration bound
    double penBound_;

    double ct_;

    /** Approximated maximal number of master elements, in which one slave element
     *  will project. */
    int maxNumMasterElements_;

    bool issetup_;

    //! slave dofs in normal direction
    Teuchos::RCP<Epetra_Map>    sndofrowmap_;

    //! slave dofs in tangential direction
    Teuchos::RCP<Epetra_Map>    stdofrowmap_;

  }; // class CONTACT::AUG::IDataContainer

  /*--------------------------------------------------------------------------*/
  /** \brief Augmented contact interface class
   *
   *  \author hiermeier \date 03/17 */
  class Interface : public ::CONTACT::CoInterface
  {
  public:
    /** \brief Alternative constructor
     *
     *  A prerequisite for this constructor is, that the passed
     *  shared interface data object has been filled/initialized already.
     *
     *  \param idata_ptr (in) : filled shared augmented contact interface
     *                          data container object
     *
     *  \author hiermeier \date 03/17 */
    Interface( const Teuchos::RCP<CONTACT::AUG::IDataContainer>& idata_ptr );

    //! Constructor
    Interface(
        const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr,
        int id,
        const Epetra_Comm& comm,
        int dim,
        const Teuchos::ParameterList& icontact,
        bool selfcontact,
        INPAR::MORTAR::RedundantStorage redundant);

    //! Destructor
    virtual ~Interface() {}

    //! @name Accessors
    //! @{

    /// share the data with other derived interfaces via copy constructor
    Teuchos::RCP<AUG::IDataContainer> SharedInterfaceDataPtr() const
    {
      return idata_ptr_;
    }

    //! Get row map of slave normal dofs
    const Teuchos::RCP<Epetra_Map> SlaveRowNDofs() const
    {
      if (Filled()) return idata_.SNDofRowMap();
      else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
      exit(EXIT_FAILURE);
    }

    //! Get row map of slave tangential dofs
    const Teuchos::RCP<Epetra_Map> SlaveRowTDofs() const
    {
      if (Filled()) return idata_.STDofRowMap();
      else dserror("CONTACT::AugmentedInterface::FillComplete was not called");
      exit(EXIT_FAILURE);
    }

    //! Returns the penetration bound of the current interface.
    inline double PenetrationBound() const
    {
      return idata_.PenBound();
    };

    //! @}

    //! @name Initialize and evaluate interface element contributions
    //! @{
    /*! \brief Finalize the construction of an augmented contact interface
     *
     *  We look up the minimal edge length only once at the construction point
     *  and call the standard mortar_interface FillComplete routine. */
    virtual void FillComplete(int maxdof=0, bool newghosting=true);

    //! Initialization of all augmented contact related quantities
    void Initialize();

    /*! \brief Reduced evaluate of the contact interface. We don't need any
     *  search algorithm and we don't have to build the nodal normals again. */
    void RedEvaluate(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    //@}

    //! @name Assemble partial global matrices
    //! @{
    //! @name Assemble auxiliary terms
    //! @{
    /*! \brief Assemble the area vector
     *
     *  We need it for scaling purposes. */
    void AssembleAugAVector(Epetra_Vector& augAVec, Epetra_Vector& kappaVec) const;


    //! Assemble the directly calculated Uzawa update increment
    void AssembleUzawaUpdateIncr(Epetra_Vector& znincrVec,
                                 const Epetra_Vector& cnVec);
    //! @}

    //! @name Assemble scalar quantities
    //! @{

    void AssembleContactPotentialTerms(
        const Epetra_Vector& cnVec,
        double& zn_gn,
        double& gn_gn,
        double& zn_zn,
        double& zt_zt ) const;
    //! @}

    //! @name Assemble right hand side vectors, Dn and Mn
    //! @{

    //! Assemble Dn and Mn matrices (equivalent to D and M in the Lagrange formulation)
    void AssembleAugDnMnMatrix(
        LINALG::SparseMatrix& augDnMatrix,
        LINALG::SparseMatrix& augMnMatrix) const;

    /*! Assemble  the averaged weighted gap vector and the weighted gap
     *  vector  */
    void AssembleGapVectors(
        Epetra_Vector& aWGapVec,
        Epetra_Vector& wGapVec) const;

    /*! Assemble the normal Lagrange multiplier vector */
    void AssembleLmNVector(
        Epetra_Vector& lmNVec) const;

    //! Assemble tangential constraint equation for active nodes (dLmTLmTrhs, ACTIVE)
    void AssembleDLmTLmTRhs(Epetra_Vector& dLmTLmTRhs) const;

    /*! \brief Assemble inactive constraint rhs
     *
     *  Combination of tangential and normal parts. */
    void AssembleAugInactiveRhs(
        Epetra_Vector& augInactiveRhs,
        Epetra_Vector& cnVec) const;
    //! @}

    //! @name Assemble linearization matrices
    //! @{
    /*! \brief Assemble DGLmLinMatrix
     *
     *  Linearization w.r.t. the displacements */
    void AssembleDGLmLinMatrix(
        LINALG::SparseMatrix& dGLmSlLinMatrix,
        LINALG::SparseMatrix& dGLmMaLinMatrix) const;

    /*! \brief Assemble DGGLinMatrix
     *
     *  Linearization w.r.t. the displacements */
    virtual void AssembleDGGLinMatrix(
        LINALG::SparseMatrix& dGGSlLinMatrix,
        LINALG::SparseMatrix& dGGMaLinMatrix,
        const Epetra_Vector& cnVec) const;

    /*! \brief Assemble DLmNWGapLinMatrix
     *
     *  Linearization w.r.t. the displ. */
    void AssembleDLmNWGapLinMatrix(
        LINALG::SparseMatrix& dLmNWGapLinMatrix) const;

    /*! \brief Assemble DLmTLmTMatrix
     *
     *  Linearization w.r.t. the LM */
    void AssembleDLmTLmTMatrix(
        LINALG::SparseMatrix& dLmTLmTMatrix) const;

    /*! \brief Assemble DLmTLmTLinMatrix
     *
     *  Linearization w.r.t. the displ. */
    void AssembleDLmTLmTLinMatrix(
        LINALG::SparseMatrix& dLmTLmTLinMatrix) const;

    /*! \brief Assemble AugInactiveMatrix
     *
     *  Linearization w.r.t. the LM */
    void AssembleAugInactiveDiagMatrix(
        Epetra_Vector& augInactiveDiagMatrix,
        const Epetra_Vector& cnVec) const;

    /*! \brief Assemble AugInactiveLinMatrix
     *
     *  Linearization w.r.t. the displ. */
    void AssembleAugInactiveLinMatrix(
        LINALG::SparseMatrix& augInactiveLinMatrix,
        const Epetra_Vector& cnVec) const;
    //! @}

    //! @name Calculate nodal quantities
    //! @{
    //! Calculate the nodal weighted gap values
    void WGap() const;

    //! Calculate the linearization of the averaged weighted gap
    void AWGapLin() const;

    //! @}

    //! @name Augmented active set functions
    //! @{

    //! Update augmented active set for each interface
    virtual bool BuildActiveSet(bool init=false);

    //! Split augmented active dof set for each interface
    void SplitAugActiveDofs();

    Teuchos::RCP<Epetra_Map> BuildActiveForceMap( const Epetra_Vector& force,
        const double threshold = 0.0 ) const;
    //! @}

    //! @name Basic setup (only necessary in the augmented lagrange case)
    //! @{

    //! Split the sdofrowmap_ into a normal and tangential part
    void SplitSlaveDofs();
    //! @}

    //! @name Debugging methods of the augmented Lagrange formulation
    //! @{
    //! Check linearization of kappa with finite differences
    void FDCheckKappaLin(CONTACT::ParamsInterface& cparams);

    //! Check linearization of augA with finite differences
    void FDCheckAugALin(CONTACT::ParamsInterface& cparams);

    //! Check linearization of the averaged weighted gap with finite differences
    void FDCheckAWGapLin(CONTACT::ParamsInterface& cparams);

    //! Check linearization of varWGapSl with finite differences
    void FDCheckVarWGapLinSl(CONTACT::ParamsInterface& cparams);

    //! Check linearization of varWGapMa with finite differences
    void FDCheckVarWGapLinMa(CONTACT::ParamsInterface& cparams);

    //! Update Interfaces for global finite difference-check
    bool UpdateInterfaces(
        int gid,int dof,
        double delta,bool forward,
        CONTACT::ParamsInterface& cparams);

    //! @}
  protected:
    //! derived
    virtual void UpdateMasterSlaveSets();

    /** \brief Assemble a map object of type T into a SparseMatrix
     *
     *  Standard as well as FE sparse matrices are supported. The global column
     *  ID is supposed to be the key of the given map.
     *
     *  \param row       (in) : global row ID
     *  \param scal      (in) : scalar for scaling the map data values
     *  \param values    (in) : map object holding the data
     *  \param mat      (out) : Sparse Matrix object (standard or FE)
     *  \param threshold (in) : values with an absolute value below this threshold
     *                          will not be assembled.
     *
     *  \author hiermeier \date 03/17 */
    template< class T >
    void AssembleMapIntoMatrix(
        int row,
        double scal,
        const T& values,
        LINALG::SparseMatrix& mat,
        double threshold = 0.0 ) const;

    virtual void AssembleDGGLinMatrixOnSlaveSide(
        const CoNode& cnode,
        const GEN::pairedvector<int,std::pair<int,double> >& varWGapSlMap,
        const std::map<int,double>& aWGapLinMap,
        double cn,
        double aWGap,
        LINALG::SparseMatrix& dGGSlLinMatrix ) const;

    virtual void AssembleDGGLinMatrixOnMasterSide(
        const CoNode& cnode,
        const std::map<int,std::pair<int,double> >& varWGapMaMap,
        const std::map<int,double>& aWGapLinMap,
        double cn,
        double aWGap,
        LINALG::SparseMatrix& dGGMaLinMatrix ) const;

    virtual void AddKappaLinToGapLinearization(
        const std::map<int,double>& kappaLinMap,
        double x,
        double kappainv,
        double varWGap,
        double scale,
        std::map<int,double>& aWGapLinMap ) const;

    void Setup();

  protected:
    // don't want = operator and cctor
    Interface( const CONTACT::AUG::Interface& source );
    Interface operator = (const Interface& old);

  private:
    /** \remark Please add no new member variables to this class and use the
     *  corresponding data container, instead! If you have any questions
     *  concerning this, you can ask me.
     *                                                        hiermeier 03/17 */

    /// pointer to the interface data object
    Teuchos::RCP<AUG::IDataContainer> idata_ptr_;

    /// reference to the interface data object
    AUG::IDataContainer&              idata_;

  };  // class Interface

} // namespace AUG
} // namespace CONTACT


#endif /* CONTACT_AUGMENTED_INTERFACE_H_ */
