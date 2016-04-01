/*!----------------------------------------------------------------------
\file meshtying_abstract_strategy.H

<pre>
\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*-----------------------------------------------------------------------*/
#ifndef MESHTYING_ABSTRACT_STRATEGY_H
#define MESHTYING_ABSTRACT_STRATEGY_H

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_inpar/inpar_mortar.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace MORTAR
{
  class MortarInterface;
}

namespace LINALG
{
  class SparseMatrix;
}

namespace CONTACT
{

/*!
\brief Main abstract class for meshtying solution strategies

This is the templating abstract class for all meshtying solution algorithms.
Every solution algorithm has to fit into the set of functions and calls defined herein
and has to be specified in a corresponding subclass defining the concrete algorithmic steps.

This class it itself derived from the MORTAR::StrategyBase class, which is an even
more abstract framework for any solution strategies involving mortar coupling.

\author popp (popp@lnm.mw.tum.de)
*/
class MtAbstractStrategy : public MORTAR::StrategyBase
{
  public:

    /*!
    \brief Standard Constructor

    Creates the strategy object and initializes all global variables, including
    all necessary Epetra_Maps and global vector and matrix quantities.

    */
    MtAbstractStrategy(
        const Epetra_Map* DofRowMap,
        const Epetra_Map* NodeRowMap,
        Teuchos::ParameterList params,
        std::vector<Teuchos::RCP<MORTAR::MortarInterface> > interface,
        int dim,
        Teuchos::RCP<Epetra_Comm> comm,
        double alphaf,
        int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~MtAbstractStrategy() {};

    //! @name Access methods

    /*!
    \brief Return Lagrange multiplier vector (t_n+1)

    */
    Teuchos::RCP<Epetra_Vector> LagrMult() { return z_; }

    /*!
    \brief Return old Lagrange multiplier vector (t_n)

    */
    Teuchos::RCP<Epetra_Vector> LagrMultOld() { return zold_; }

    /*!
    \brief Return old rescaled Lagrange multiplier vector (t_n)

    */
    Teuchos::RCP<Epetra_Vector> LagrMultOldRescaled() ;

    /*!
    \brief Return Lagrange multiplier vector from last Uzawa step

    */
    Teuchos::RCP<Epetra_Vector> LagrMultUzawa() { return zuzawa_; }

    /*!
    \brief Return constraint rhs vector (only in saddle-point formulation

    */
    Teuchos::RCP<Epetra_Vector> ConstrRhs() { return constrrhs_; }

    /*!
    \brief Returns increment of LagrangeMultiplier solution vector in SaddlePointSolve routine

    */
    Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() { return zincr_; }

    /*!
    \brief provides maps needed for contact/meshtying specific multigrid preconditioners

    */
    void CollectMapsForPreconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap, Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap, Teuchos::RCP<Epetra_Map>& ActiveDofMap );

    /*!
    \brief Return mortar matrix D

    */
    Teuchos::RCP<LINALG::SparseMatrix> DMatrix() { return dmatrix_; }

    /*!
    \brief Return mortar matrix M

    */
    Teuchos::RCP<LINALG::SparseMatrix> MMatrix() { return mmatrix_; }

    /*!
    \brief Get dual quadratic 3d slave element flag

    Returns TRUE if at least one higher-order 3d slave element with
    dual Lagrange mutliplier shape functions in any interface.

    */
    virtual const bool& Dualquadslave3d() { return dualquadslave3d_; };

    /*!
    \brief Return parallel redistribution status (yes or no)

    */
    bool ParRedist()
    {
      INPAR::MORTAR::ParRedist partype = DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(Params(),"PARALLEL_REDIST");
      if (partype!=INPAR::MORTAR::parredist_none) return true;
      else return false;
    }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Redistribute all meshtying interfaces in parallel

    */
    void RedistributeMeshtying();

    /*!
    \brief Global evaluation method called from time integrator

    */
    void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt,
                            Teuchos::RCP<Epetra_Vector>& f,
                            const int step, const int iter,
                            bool predictor = false);

    /*!
    \brief Set current deformation state

    All interfaces are called to set the current deformation state.

    \param statename (in): std::string defining which quantity to set (only "displacement" applicable)
    \param vec (in): current global state of the quantity defined by statename

    */
    void SetState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec);

    /*!
    \brief Do mortar coupling in reference configuration

    Only do this ONCE for meshtying upon initialization!
    This method calls Initialize() on all contact interfaces, which
    resets all kind of nodal quantities. It then calls Evaluate() on
    all meshtying interfaces, which does all the geometric coupling stuff.
    Concretely, this is an evaluation of all involved quantites at nodal
    level. It includes the nodal normal calculations, search, projection
    and overlap detection and integration of the Mortar terms D and M.

    Then - on global level - it resets the Mortar matrices D and M accordingly.
    The nodal quantites computed before are assembled to global matrices. No
    setup of the global system is to be done here yet, so there is no need to
    pass in the effective stiffness K or the effective load vector f.

    Note: Only quantities common to all subsequent solving strategies (Lagrange,
    Penalty) are computed here. In case they need additional mortar variables,
    use the overloaded function call in the derived class and refer back to this function.

    */
    virtual void MortarCoupling(const Teuchos::RCP<Epetra_Vector> dis);

    //@}

    //! @name Quantity control methods

    /*!
    \brief Get some nodal quantity globally and store into MortarNodes

    The enum input parameter defines, which quantity is be updated.
    Currently the possibilities "lmold", "lmcurrent", "lmupdate" and
    "lmuzawa" exist. Note that "lmold" means the converged value LM_n
    of the last time / load step, whereas "lmcurrent" adresses the current
    (not necessarily converged) value of the LM_n+1. "lmupdate" is a special
    option called only in Recover() after the update of the Lagr. multipliers.
    It basically does the same as "lmcurrent", but also checks for D.B.C.
    problems. Finally, "lmuzawa" addresses the LM update within an Uzawa
    augmented Lagrangian scheme.

    \param type (in): enum defining which quantity to store into MortarNodes

    */
    void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type);

    /*!
    \brief Get dirichlet B.C. status and store into MortarNodes

    This is called once at the beginning of the simulation
    to set the D.B.C. status in each MortarNodes.

    \param dbcmaps (in): MapExtractor carrying global dbc map

    */
    void StoreDirichletStatus(Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

    /*!
    \brief Update meshtying at end of time step

    \param dis (in):  current displacements (-> old displacements)

    */
    void Update(Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief Perform a write restart

    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Different from writing a restart step, now
    all the restart action has to be performed on the level of the meshtying algorithm,
    for short: here's the right place.

    */
    void DoReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis);

    //@}

    //! @name Debugging methods

    /*!
    \brief Compute interface forces and moments

    Compute current interface forces and moments at n+1-alphaf using current
    Lagrange multiplier values and current Mortar matrices D and M at n+1. When
    doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
    matrices of the last converged time / load step n (TR-like interpolation).

    \param fresm (in): residual / force vector at state n+1 of current Newton step
    \param output (in): flag indicating whether force output shall be written

    */
    void InterfaceForces(bool output = false);

    /*!
    \brief Print interfaces

    */
    void Print(std::ostream& os) const;

    /*!
    \brief Print current active set to screen for debugging purposes

    */
    void PrintActiveSet();

    /*!
    \brief Visualize contact stuff with gmsh

    \param step (in): current time step index
    \param iter (in): current iteration index

    */
    void VisualizeGmsh(const int step, const int iter);

    //@}

    //! @name Purely virtual functions

    // All these functions are defined in one or more specific derived classes,
    // i.e CONTACT::MeshtyingLagrangeStrategy or CONTACT::MeshtyingPenaltyStrategy.
    // As the base class MORTAR::StrategyBase is always called from the control routine
    // (time integrator), these functions need to be defined purely virtual here.

    virtual double ConstraintNorm() = 0;
    virtual void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual double InitialPenalty() = 0;
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi) = 0;
    virtual void ResetPenalty() = 0;
    virtual void ModifyPenalty() = 0;
    virtual void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) = 0;
    virtual void UpdateDisplacementsAndLMincrements(Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol) = 0;
    virtual void UpdateUzawaAugmentedLagrange() = 0;
    virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;

    //@}

    //! @name Empty functions (contact)

    // All these functions only have functionality in contact simulations, thus they
    // are defined as empty here in the case of meshtying. They can be called from the
    // control routine (time integrator), whenever you like.

    bool ActiveSetConverged() {return true;}
    bool ActiveSetSemiSmoothConverged() {return true;}
    bool Friction() {return false;}
    bool WearBothDiscrete() {return false;}
    bool IsInContact() { return true; }
    bool WasInContact() { return true; }
    bool WasInContactLastTimeStep() { return true; }
    Teuchos::RCP<Epetra_Vector> ContactNorStress() {return Teuchos::null;}
    Teuchos::RCP<Epetra_Vector> ContactTanStress() {return Teuchos::null;}
    void AssembleMortar() {}
    void InitBinStrategyforTimestep(Teuchos::RCP<Epetra_Vector> vel) {}
    void DoWriteRestart(Teuchos::RCP<Epetra_Vector>& activetoggle, Teuchos::RCP<Epetra_Vector>& sliptoggle, Teuchos::RCP<Epetra_Vector>& weightedwear,Teuchos::RCP<Epetra_Vector>& realwear, bool forcedrestart = false) {}
    void InitEvalInterface() {}
    void InitMortar() {}
    void Initialize() {}
    double Inttime() {return inttime_;};
    void Inttime_init(){inttime_=0.0;};
    int NumberOfActiveNodes() {return 0;}
    int NumberOfSlipNodes() {return 0;}
    void OutputStresses() {}
    void AugForces(Epetra_Vector& augfs_lm,Epetra_Vector& augfs_g,Epetra_Vector& augfm_lm,Epetra_Vector& augfm_g){};
    void RedistributeContact(Teuchos::RCP<Epetra_Vector> dis) {}
    void ResetActiveSet() {}
    void SaveReferenceState(Teuchos::RCP<Epetra_Vector> dis) {}
    void UpdateActiveSet() {}
    void UpdateActiveSetSemiSmooth() {}
    Teuchos::RCP<LINALG::SparseMatrix> EvaluateNormals(Teuchos::RCP<Epetra_Vector> dis) {return Teuchos::null;}
    void EvaluateReferenceState(const Teuchos::RCP<Epetra_Vector> vec){}
    void EvaluateRelMov() {}
    void EvaluateRelMovPredict() {}
    const Teuchos::RCP<Epetra_Map> SlaveRowNodes() { return gsnoderowmap_; }
    const Teuchos::RCP<Epetra_Map> ActiveRowNodes() {return Teuchos::null;}
    const Teuchos::RCP<Epetra_Map> ActiveRowDofs() {return Teuchos::null;}
    const Teuchos::RCP<Epetra_Map> NotReDistSlaveRowDofs() { return pgsdofrowmap_; }
    const Teuchos::RCP<Epetra_Map> NotReDistMasterRowDofs() { return pgmdofrowmap_; }
    const Teuchos::RCP<Epetra_Map> SlipRowNodes() {return Teuchos::null;}
    //@}

  protected:
    /*!
    \brief Assemble global coordinate vector

    \param sidename (in): std::string indicating slave or master side
    \param ref (in): boolean indicating evaluation in reference configuration
    \param vec (in/out)):  empty global vetcor to be assembled to

    */
    void AssembleCoords(
        const std::string& sidename,
        bool ref,
        Teuchos::RCP<Epetra_Vector> vec);

    /*!
    \brief Do mesh intialization for rotational invariance

    Only do this ONCE for meshtying upon initialization!
    This method relocates the slave nodes such that the meshtying constraint
    is satisifed in the reference condiguration, which is a prerequisite for
    ensuring both rotational invariance and absence of initial stresses at the
    same time. Basically the constraint equation needs to be solved for this,
    which is specific to the applied solving strategy (dual Lagrange or Penalty).
    In the dual LM, matrix D is diagonal, thus its inversion is trivial and no
    linear system needs to be solved. In the penalty case, matrix D is not diagonal
    and we apply a default LINALG::Solver to solve for the modified slave poitions.
    Thus, this linear system solve is done in the derived method FIRST and then
    we refer back to this base class function.

    \param Xslavemod (in): modified slave reference configuration

    */
    virtual void MeshInitialization(
        Teuchos::RCP<Epetra_Vector> Xslavemod);

  private:
    /*!
    \brief Evaluate contact

    This is just a tiny control routine, deciding which Evaluate-routine
    of those listed below is to be called (based on input-file information)
    Note that into ALL derived Evaluate() routines, a REFERENCE to the pointer
    on the effective stiffness matrix is handed in. This way, after building the
    new effective stiffness matrix with contact, we can simply let the pointer
    kteff point onto the new object. The same is true for the effective force
    vector feff. Be careful: kteff is of type Teuchos::RCP<LINALG::SparseOperator>&.

    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)

    */
    void Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff,
                  Teuchos::RCP<Epetra_Vector>& feff,
                  Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief Restrict slave boundary to actual meshtying zone

    Only do this ONCE for meshtying upon initialization!
    This method first detects for each interface the actually tied part
    of the slave surface (i.e. the nodes that carry a D/M contribution).
    Then all slave maps on interface level and on global level are
    re-initialized and re-setup according to the the above defined
    actual slave meshtying zone. This is necessary for problems in which
    the slave surface does not fully project onto the master surface
    and thus the actual meshtying zone cannot be defined within the
    input file. Thus, it is computed here.

    */
    virtual void RestrictMeshtyingZone();

    /*!
    \brief Setup this strategy object (maps, vectors, etc.)

    All global maps and vectors are initialized by collecting
    the necessary information from all interfaces. In the case
    of a parallel redistribution, this method is called again
    to re-setup the above mentioned quantities. In this case
    the input parameter is set to TRUE.

    */
    void Setup(bool redistributed);

  protected:

    // don't want cctor (= operator impossible anyway for abstract class)
    MtAbstractStrategy(const MtAbstractStrategy& old);

    std::vector<Teuchos::RCP<MORTAR::MortarInterface> > interface_; // meshtying interfaces

    Teuchos::RCP<Epetra_Map> glmdofrowmap_;         // global Lagrange mult. dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsdofrowmap_;          // global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmdofrowmap_;          // global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gndofrowmap_;          // global internal dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsmdofrowmap_;         // global slave and master dof row map (s+m map)
    Teuchos::RCP<Epetra_Map> gdisprowmap_;          // global displacement dof row map (s+m+n map)

    Teuchos::RCP<Epetra_Map> gsnoderowmap_;         // global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmnoderowmap_;         // global master node row map (of all interfaces)

    Teuchos::RCP<Epetra_Map> pglmdofrowmap_;        // global LM dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsdofrowmap_;         // global slave dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgmdofrowmap_;         // global master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsmdofrowmap_;        // global slave and master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Vector> pgsdirichtoggle_;   // global dirichlet toggle of all slave dofs (before parallel redistribution)

    std::vector<Teuchos::RCP<Epetra_Map> > initial_elecolmap_; // init a col ele map for binning strategy (s m)

    Teuchos::RCP<LINALG::SparseMatrix> dmatrix_;    // global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix> mmatrix_;    // global Mortar matrix M
    Teuchos::RCP<Epetra_Vector> g_;                 // global weighted gap vector g
    Teuchos::RCP<Epetra_Vector> constrrhs_;         // global constraint rhs vector (only for saddlepoint problems)

    Teuchos::RCP<Epetra_Vector> z_;                 // current vector of Lagrange multipliers at t_n+1
    Teuchos::RCP<Epetra_Vector> zold_;              // old vector of Lagrange multipliers at t_n
    Teuchos::RCP<Epetra_Vector> zincr_;             // Lagrange multiplier vector increment within SaddlePointSolve (this is NOT the increment of z_ between t_{n+1} and t_{n}!)
    Teuchos::RCP<Epetra_Vector> zuzawa_;            // vector of Lagrange multipliers from last Uzawa step

    bool                      dualquadslave3d_;     // flag indicating whether trafo should be applied
    Teuchos::RCP<LINALG::SparseMatrix> trafo_;      // transformation matrix T for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix> invtrafo_;   // inverse trafo matrix T^(-1) for dual quad 3D case
    double inttime_;                                // Integration time

}; // class MtAbstractStrategy
} // namespace CONTACT

// << operator
std::ostream& operator << (std::ostream& os, const CONTACT::MtAbstractStrategy& strategy);

#endif  // #ifndef MESHTYING_ABSTRACT_STRATEGY_H
