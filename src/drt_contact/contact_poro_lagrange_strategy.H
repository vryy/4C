/*!----------------------------------------------------------------------
\file contact_poro_lagrange_strategy.H

<pre>
Maintainer: Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/
#ifndef CONTACT_PORO_LAGRANGE_STRATEGY_H
#define CONTACT_PORO_LAGRANGE_STRATEGY_H

#include "contact_monocoupled_lagrange_strategy.H"
#include "../drt_lib/drt_dserror.H"

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }
}

namespace CONTACT
{

// forward declaration
//class WearInterface;
/*!
 \brief Contact solving strategy with (standard/dual) Lagrangian multipliers.

 This is a specialization of the abstract contact algorithm as defined in CoAbstractStrategy.
 For a more general documentation of the involved functions refer to CoAbstract Strategy.

 \author popp (popp@lnm.mw.tum.de)
 */
class PoroLagrangeStrategy : public MonoCoupledLagrangeStrategy
{
  public:

    /*!
    \brief Standard Constructor

    */
    PoroLagrangeStrategy(
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface,
        int dim,
        Teuchos::RCP<Epetra_Comm> comm,
        double alphaf,
        int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~PoroLagrangeStrategy() {};

    //! @name Access methods

    //@}

    //! @name Evaluation methods

    /*!
    /brief Activate the No Penetration Condition for the active contact surface!
     */
    void SetupNoPenetrationCondition();

    /*!
    \brief Initialize poro contact variables for next Newton step

    For a poro lagrangian strategy this includes the global normal / tangent matrices N and T, //Todo to be updated
    the global derivative matrices S and P and Tresca friction matrix L + vector r.

    */
   void PoroInitialize(ADAPTER::Coupling& coupfs , Teuchos::RCP<const Epetra_Map> fluiddofs);

    /*!
    \brief Evaluate poro no penetration contact

     Evaluate poro coupling contact matrices for no penetration
     condition on contact surface
    */
    void EvaluatePoroNoPenContact(Teuchos::RCP<LINALG::SparseMatrix>& k_fseff, //global poro Coupling Matrix Fluid Structure K_FS
                                  Teuchos::RCP<LINALG::SparseMatrix>& Feff, //global fluid Matrix in poro
                                  Teuchos::RCP<Epetra_Vector>& feff); //global RHS of fluid in poro

    /*!
    \brief Recovery method

    We only recover the Lagrange multipliers for poro no penetration condition here, which had been
    statically condensated during the setup of the global problem!*/

    void RecoverPoroNoPen(Teuchos::RCP<Epetra_Vector> disi,Teuchos::RCP<Epetra_Vector> veli);

    void UpdatePoroContact();

    /*!
    \brief Set current state
    ...Standard Implementation in Abstract Strategy:
    All interfaces are called to set the current deformation state
    (u, xspatial) in their nodes. Additionally, the new contact
    element areas are computed.

    ... + Overloaded Implementation in Poro Lagrange Strategy
    Set structure & fluid velocity and lagrangean multiplier to Contact nodes data container!!!

    \param statename (in): std::string defining which quantity to set (either "displacement", "olddisplacement", "svelocity", "fvelocity", "lm")
    \param vec (in): current global state of the quantity defined by statename
    */
    void SetState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec);

    void SetParentState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec, const Teuchos::RCP<DRT::Discretization> dis);

    //Flag for Poro No Penetration Condition
    bool HasPoroNoPenetration(){return no_penetration_;}

    //Return Lagrange Multiplier for No Penetration Condition!
    Teuchos::RCP<Epetra_Vector> LambdaNoPen(){return lambda_;}

    protected:

    // don't want = operator and cctor
    PoroLagrangeStrategy operator = (const PoroLagrangeStrategy& old);
    PoroLagrangeStrategy(const PoroLagrangeStrategy& old);

    // flag activation poro contact no penetration condition
    bool no_penetration_;

    //time integration
    double nopenalpha_; //1-theta!!!

    Teuchos::RCP<LINALG::SparseMatrix> csfsn_;  //poro coupling stiffness block Csf_sn (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> csfsm_;  //poro coupling stiffness block Csf_sm (needed for LM)
    Teuchos::RCP<LINALG::SparseMatrix> csfss_;  //poro coupling stiffness block Csf_ss (needed for LM)

    //For Recovery of no penetration lagrange multiplier!!!
    Teuchos::RCP<Epetra_Vector> ffs_; //poro fluid RHS (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix> cfssn_; //poro coupling stiffness block Cfs_sn (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix> cfssm_; //poro coupling stiffness block Cfs_sm (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix> cfsss_; //poro coupling stiffness block Cfs_ss (needed for no pen LM)

    Teuchos::RCP<LINALG::SparseMatrix> fsn_; //poro fluid stiffness block F_sn (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix> fsm_; //poro fluid stiffness block F_sm (needed for no pen LM)
    Teuchos::RCP<LINALG::SparseMatrix> fss_; //poro fluid stiffness block F_ss (needed for no pen LM)

    //Matrices transformed to the fluid dofs!!!
    Teuchos::RCP<LINALG::SparseMatrix> fdhat_;
    Teuchos::RCP<LINALG::SparseMatrix> fmhataam_;
    Teuchos::RCP<LINALG::SparseMatrix> finvda_;
    Teuchos::RCP<LINALG::SparseMatrix> ftanginvD_;

    Teuchos::RCP<LINALG::SparseMatrix> fdoldtransp_;       // global transposed Mortar matrix D (last end-point t_n)
    Teuchos::RCP<LINALG::SparseMatrix> fmoldtransp_;       // global transposed Mortar matrix M (last end-point t_n)


    Teuchos::RCP<Epetra_Map> fgsdofrowmap_; //fluid slave dofs
    Teuchos::RCP<Epetra_Map> fgmdofrowmap_; //fluid master dofs
    Teuchos::RCP<Epetra_Map> fgsmdofrowmap_; //fluid slave + master dofs
    Teuchos::RCP<Epetra_Map> fgndofrowmap_; //fluid other dofs
    Teuchos::RCP<Epetra_Map> fgactivedofs_; //fluid active slave dofs
    Teuchos::RCP<Epetra_Map> falldofrowmap_; //all fluid dofs
    Teuchos::RCP<Epetra_Map> fgactiven_; //active normal fluid dofs
    Teuchos::RCP<Epetra_Map> fgactivet_; //active tangential fluid dofs

    /// @name matrix transformation
    //! transform object for linearized ncoup matrix \f$linncoup\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform> linncoupveltransform_;
    //! transform object for linearized ncoup matrix \f$linncoup\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> linncoupdisptransform_;
    //! transform object for linearized ncoup matrix \f$linncoup\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform> tanginvtransform_;
    //! transform object for linearized tangentlambda matrix \f$lintanglambda\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> lintangentlambdatransform_;
    //! transform object for linearized Dlambda matrix \f$linDlambda\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> porolindmatrixtransform_;
    //! transform object for mhataam = invda * mmatrixa
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> mhataamtransform_;
    //! transform object for dhat
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> dhattransform_;
    //! transform object for mold
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> doldtransform_;
    //! transform object for dold
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> moldtransform_;
     //! transform object for active part of inverse D matrix \f$invDa\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> invDatransform_;


    Teuchos::RCP<Epetra_Vector> lambda_;                 // current vector of Lagrange multipliers(for poro no pen.) at t_n+1
    Teuchos::RCP<Epetra_Vector> lambdaold_;              // old vector of Lagrange multipliers(for poro no pen.) at t_n

    //... add the relevant matrices !!!
    Teuchos::RCP<Epetra_Vector>            NCoup_;           ///<normal coupling vector (for RHS)
    Teuchos::RCP<LINALG::SparseMatrix>     NCoup_lindisp_;   ///<linearisation of normal coupling w.r.t. displacements
    Teuchos::RCP<LINALG::SparseMatrix>     NCoup_linvel_;    ///<linearisation of normal coupling w.r.t. fluid velocity

    Teuchos::RCP<Epetra_Vector>            fNCoup_;           ///<normal coupling vector (for RHS) -- transformed to fluid dofs
    Teuchos::RCP<LINALG::SparseMatrix>     fNCoup_lindisp_;   ///<linearisation of normal coupling w.r.t. displacements -- transformed to fluid dofs
    Teuchos::RCP<LINALG::SparseMatrix>     fNCoup_linvel_;    ///<linearisation of normal coupling w.r.t. fluid velocity -- transformed to fluid dofs

    Teuchos::RCP<LINALG::SparseMatrix>     Tangential_;      ///<matrix with tangential vectors inside

    Teuchos::RCP<LINALG::SparseMatrix>     linTangentiallambda_;  ///<linearized tangential times lambda

    Teuchos::RCP<LINALG::SparseMatrix>     flinTangentiallambda_;  ///<linearized tangential times lambda -- transformed to fluid dofs

    Teuchos::RCP<LINALG::SparseMatrix> porolindmatrix_; // global Matrix LinD containing slave fc derivatives (with lm from poro no penetration)
    Teuchos::RCP<LINALG::SparseMatrix> porolinmmatrix_; // global Matrix LinM containing master fc derivatives (with lm from poro no penetration)

    Teuchos::RCP<LINALG::SparseMatrix> fporolindmatrix_; // global Matrix LinD containing slave fc derivatives (with lm from poro no penetration) -- transformed to fluid dofs
    Teuchos::RCP<LINALG::SparseMatrix> fporolinmmatrix_; // global Matrix LinM containing master fc derivatives (with lm from poro no penetration) -- transformed to fluid dofs

}; // class POROLagrangeStrategy
} // namespace CONTACT


#endif  // #ifndef CONTACT_PORO_LAGRANGE_STRATEGY_H
