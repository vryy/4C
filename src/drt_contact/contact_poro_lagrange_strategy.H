/*---------------------------------------------------------------------*/
/*! \file
\brief Contact Strategy handling the porous no penetraction condition on the active contact
interface

\level 3

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_PORO_LAGRANGE_STRATEGY_H
#define CONTACT_PORO_LAGRANGE_STRATEGY_H

#include "contact_monocoupled_lagrange_strategy.H"
#include "../drt_lib/drt_dserror.H"

namespace ADAPTER
{
  class Coupling;
}


namespace LINALG
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace LINALG

namespace CONTACT
{
  class PoroLagrangeStrategy : public MonoCoupledLagrangeStrategy
  {
   public:
    /*!
    \brief Standard Constructor

    */
    PoroLagrangeStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* DofRowMap, const Epetra_Map* NodeRowMap, Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface>> interface, int dim,
        Teuchos::RCP<Epetra_Comm> comm, double alphaf, int maxdof, bool poroslave, bool poromaster);

    /*!
    \brief Destructor

    */
    virtual ~PoroLagrangeStrategy(){};

    //! @name Access methods

    //@}

    //! @name Evaluation methods

    virtual void DoReadRestart(IO::DiscretizationReader& reader,
        Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr);

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

    /*!
    /brief Activate the No Penetration Condition for the active contact surface!
     */
    void SetupNoPenetrationCondition();

    /*!
    \brief Initialize poro contact variables for next Newton step

    For a poro lagrangian strategy this includes the global normal / tangent matrices N and T,
    //Todo to be updated the global derivative matrices S and P and Tresca friction matrix L +
    vector r.
    */
    void PoroInitialize(
        ADAPTER::Coupling& coupfs, Teuchos::RCP<const Epetra_Map> fluiddofs, bool fullinit = true);

    /*!
    \brief as D and M Matrices are initialized here
    */
    void PoroMtInitialize();

    /*!
    \brief Prepare Matrices D and M, that are not Computed coming from the mortar adapter
    */
    void PoroMtPrepareFluidCoupling();

    /*!
    \brief
    set some coupling matrices for the proro meshtying case mhataam-, dhat_ and invda_
    */
    void PoroMtSetCouplingMatrices();

    /*!
    \brief
    set old matrices dold_ mold_ and the old lagrange multiplier in case of poro meshtying
    */
    void PoroMtUpdate();

    /*!
    \brief Evaluate poro no penetration contact

     Evaluate poro coupling contact matrices for no penetration
     condition on contact surface
    */
    void EvaluatePoroNoPenContact(Teuchos::RCP<LINALG::SparseMatrix>&
                                      k_fseff,  // global poro Coupling Matrix Fluid Structure K_FS
        Teuchos::RCP<LINALG::SparseMatrix>& Feff,  // global fluid Matrix in poro
        Teuchos::RCP<Epetra_Vector>& feff);        // global RHS of fluid in poro

    /*!
    \brief Evaluate poro no penetration contact

     Evaluate poro coupling contact matrices for no penetration
     condition on contact surface
    */
    void EvaluatePoroNoPenContact(Teuchos::RCP<LINALG::SparseMatrix>& k_fseff,
        std::map<int, Teuchos::RCP<LINALG::SparseMatrix>*>& Feff,
        Teuchos::RCP<Epetra_Vector>& feff);

    /*!
    \brief Evaluate poro no penetration contact

     Evaluate all other matrixes here!
    */
    void EvaluateMatPoroNoPen(Teuchos::RCP<LINALG::SparseMatrix>&
                                  k_fseff,   // global poro Coupling Matrix Fluid Structure K_FS
        Teuchos::RCP<Epetra_Vector>& feff);  // global RHS of fluid in poro

    /*!
    \brief Evaluate poro no penetration contact

     Evaluate all other matrixes here!
    */
    void EvaluateOtherMatPoroNoPen(Teuchos::RCP<LINALG::SparseMatrix>& Feff, int Column_Block_Id);

    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers for poro no penetration condition here, which had been
    statically condensated during the setup of the global problem!*/

    void RecoverPoroNoPen(Teuchos::RCP<Epetra_Vector> disi, Teuchos::RCP<Epetra_Vector> inc);

    void RecoverPoroNoPen(
        Teuchos::RCP<Epetra_Vector> disi, std::map<int, Teuchos::RCP<Epetra_Vector>> inc);


    void UpdatePoroContact();

    /*!
    \brief Set current state
    ...Standard Implementation in Abstract Strategy:
    All interfaces are called to set the current deformation state
    (u, xspatial) in their nodes. Additionally, the new contact
    element areas are computed.

    ... + Overloaded Implementation in Poro Lagrange Strategy
    Set structure & fluid velocity and lagrangean multiplier to Contact nodes data container!!!

    \param statetype (in): enumerator defining which quantity to set (see mortar_interface.H for an
    overview) \param vec (in): current global state of the quantity defined by statetype
    */
    virtual void SetState(const enum MORTAR::StateType& statetype, const Epetra_Vector& vec);

    virtual void SetParentState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec,
        const Teuchos::RCP<DRT::Discretization> dis);

    // Flag for Poro No Penetration Condition
    bool HasPoroNoPenetration() const { return no_penetration_; }

    // Return Lagrange Multiplier for No Penetration Condition!
    Teuchos::RCP<Epetra_Vector> LambdaNoPen() { return lambda_; }
    Teuchos::RCP<const Epetra_Vector> LambdaNoPen() const { return lambda_; }

    // Return all active fluid slave dofs
    Teuchos::RCP<Epetra_Map> FluidActiveNDofMap() { return fgactiven_; };

   protected:
    // don't want = operator and cctor
    PoroLagrangeStrategy operator=(const PoroLagrangeStrategy& old) = delete;
    PoroLagrangeStrategy(const PoroLagrangeStrategy& old) = delete;

    // flag activation poro contact no penetration condition
    // h.Willmann the name is misleading as the bool is also used for other cases to access some
    // methods
    bool no_penetration_;

    // time integration
    double nopenalpha_;  // 1-theta!!!

    Teuchos::RCP<LINALG::SparseMatrix>
        csfsn_;  // poro coupling stiffness block Csf_sn (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix>
        csfsm_;  // poro coupling stiffness block Csf_sm (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix>
        csfss_;  // poro coupling stiffness block Csf_ss (needed for LM)

    // For Recovery of no penetration lagrange multiplier!!!
    Teuchos::RCP<Epetra_Vector> ffs_;  // poro fluid RHS (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix>
        cfssn_;  // poro coupling stiffness block Cfs_sn (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix>
        cfssm_;  // poro coupling stiffness block Cfs_sm (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix>
        cfsss_;  // poro coupling stiffness block Cfs_ss (needed for no pen LM)

    Teuchos::RCP<LINALG::SparseMatrix>
        fsn_;  // poro fluid stiffness block F_sn (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix>
        fsm_;  // poro fluid stiffness block F_sm (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix>
        fss_;  // poro fluid stiffness block F_ss (needed for no pen LM)

    std::map<int, Teuchos::RCP<LINALG::SparseOperator>>
        cfx_s_;  // offdiagonal coupling stiffness blocks on slave side!

    // Matrices transformed to the fluid dofs!!!
    Teuchos::RCP<LINALG::SparseMatrix> fdhat_;
    Teuchos::RCP<LINALG::SparseMatrix> fmhataam_;
    Teuchos::RCP<LINALG::SparseMatrix> finvda_;
    Teuchos::RCP<LINALG::SparseMatrix> ftanginvD_;

    Teuchos::RCP<LINALG::SparseMatrix>
        fdoldtransp_;  // global transposed Mortar matrix D (last end-point t_n)
    Teuchos::RCP<LINALG::SparseMatrix>
        fmoldtransp_;  // global transposed Mortar matrix M (last end-point t_n)


    Teuchos::RCP<Epetra_Map> fgsdofrowmap_;   // fluid slave dofs
    Teuchos::RCP<Epetra_Map> fgmdofrowmap_;   // fluid master dofs
    Teuchos::RCP<Epetra_Map> fgsmdofrowmap_;  // fluid slave + master dofs
    Teuchos::RCP<Epetra_Map> fgndofrowmap_;   // fluid other dofs
    Teuchos::RCP<Epetra_Map> fgactivedofs_;   // fluid active slave dofs
    Teuchos::RCP<Epetra_Map> falldofrowmap_;  // all fluid dofs
    Teuchos::RCP<Epetra_Map> fgactiven_;      // active normal fluid dofs
    Teuchos::RCP<Epetra_Map> fgactivet_;      // active tangential fluid dofs

    /// @name matrix transformation
    //! transform object for linearized ncoup matrix \f$linncoup\f$
    Teuchos::RCP<LINALG::MatrixRowTransform> linncoupveltransform_;
    //! transform object for linearized ncoup matrix \f$linncoup\f$
    Teuchos::RCP<LINALG::MatrixRowTransform> linncoupdisptransform_;
    //! transform object for tangential times Dinv matrix \f$T*D^-1\f$
    Teuchos::RCP<LINALG::MatrixRowColTransform> tanginvtransform_;
    //! transform object for linearized tangentlambda matrix \f$lintanglambda\f$
    Teuchos::RCP<LINALG::MatrixRowTransform> lintangentlambdatransform_;
    //! transform object for linearized Dlambda matrix \f$linDlambda\f$
    Teuchos::RCP<LINALG::MatrixRowTransform> porolindmatrixtransform_;
    //! transform object for linearized Mlambda matrix \f$linMlambda\f$
    Teuchos::RCP<LINALG::MatrixRowTransform> porolinmmatrixtransform_;  // h.Willmann
    //! transform object for mhataam = invda * mmatrixa
    Teuchos::RCP<LINALG::MatrixRowColTransform> mhataamtransform_;
    //! transform object for dhat
    Teuchos::RCP<LINALG::MatrixRowTransform> dhattransform_;
    //! transform object for mold
    Teuchos::RCP<LINALG::MatrixRowTransform> doldtransform_;
    //! transform object for dold
    Teuchos::RCP<LINALG::MatrixRowTransform> moldtransform_;
    //! transform object for active part of inverse D matrix \f$invDa\f$
    Teuchos::RCP<LINALG::MatrixRowTransform> invDatransform_;


    Teuchos::RCP<Epetra_Vector>
        lambda_;  // current vector of Lagrange multipliers(for poro no pen.) at t_n+1
    Teuchos::RCP<Epetra_Vector>
        lambdaold_;  // old vector of Lagrange multipliers(for poro no pen.) at t_n

    //... add the relevant matrices !!!
    Teuchos::RCP<Epetra_Vector> NCoup_;  ///< normal coupling vector (for RHS)
    Teuchos::RCP<LINALG::SparseMatrix>
        NCoup_lindisp_;  ///< linearisation of normal coupling w.r.t. displacements
    Teuchos::RCP<LINALG::SparseMatrix>
        NCoup_linvel_;  ///< linearisation of normal coupling w.r.t. fluid velocity

    Teuchos::RCP<Epetra_Vector>
        fNCoup_;  ///< normal coupling vector (for RHS) -- transformed to fluid dofs
    Teuchos::RCP<LINALG::SparseMatrix>
        fNCoup_lindisp_;  ///< linearisation of normal coupling w.r.t. displacements -- transformed
                          ///< to fluid dofs
    Teuchos::RCP<LINALG::SparseMatrix>
        fNCoup_linvel_;  ///< linearisation of normal coupling w.r.t. fluid velocity -- transformed
                         ///< to fluid dofs

    Teuchos::RCP<LINALG::SparseMatrix> Tangential_;  ///< matrix with tangential vectors inside

    Teuchos::RCP<LINALG::SparseMatrix> linTangentiallambda_;  ///< linearized tangential times
                                                              ///< lambda

    Teuchos::RCP<LINALG::SparseMatrix>
        flinTangentiallambda_;  ///< linearized tangential times lambda -- transformed to fluid dofs

    Teuchos::RCP<LINALG::SparseMatrix>
        porolindmatrix_;  // global Matrix LinD containing slave fc derivatives (with lm from poro
                          // no penetration)
    Teuchos::RCP<LINALG::SparseMatrix>
        porolinmmatrix_;  // global Matrix LinM containing master fc derivatives (with lm from poro
                          // no penetration)

    Teuchos::RCP<LINALG::SparseMatrix>
        fporolindmatrix_;  // global Matrix LinD containing slave fc derivatives (with lm from poro
                           // no penetration) -- transformed to fluid dofs
    Teuchos::RCP<LINALG::SparseMatrix>
        fporolinmmatrix_;  // global Matrix LinM containing master fc derivatives (with lm from poro
                           // no penetration) -- transformed to fluid dofs

    bool poroslave_;   // true if interface slave side is purely poroelastic
    bool poromaster_;  // true if interface master sided is purely poroelastic
    // it must be assured that these two are previously set correctly and that there is no mixed
    // master or slave interface with both structural and poroelastic elements
  };  // class POROLagrangeStrategy
}  // namespace CONTACT


#endif  // #ifndef CONTACT_PORO_LAGRANGE_STRATEGY_H
