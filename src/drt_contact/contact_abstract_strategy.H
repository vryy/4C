/*!----------------------------------------------------------------------
\file contact_abstract_strategy.H

<pre>
\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*-----------------------------------------------------------------------*/
#ifndef CONTACT_STRATEGY_H
#define CONTACT_STRATEGY_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>

#include <Epetra_Operator.h>

#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"

// forward declarations
namespace DRT { class Discretization; }

namespace LINALG { class MultiMapExtractor; class SparseMatrix; }

namespace CONTACT
{

// forward declarations
class CoInterface;
//class WearInterface;

/*!
 \brief Main abstract class for contact solution strategies

 This is the templating abstract class for all contact solution algorithms.
 Every solution algorithm has to fit into the set of functions and calls defined herein
 and has to be specified in a corresponding subclass defining the concrete algorithmic steps.

 This class it itself derived from the MORTAR::StrategyBase class, which is an even
 more abstract framework for any solution strategies involving mortar coupling.

 Refer also to the Semesterarbeit of Bernd Budich, 2009

 \author popp (popp@lnm.mw.tum.de)
 */
class CoAbstractStrategy: public MORTAR::StrategyBase
{
public:

  /*!
   \brief Standard Constructor

   Creates the strategy object and initializes all global variables, including
   all necessary Epetra_Maps and global vector and matrix quantities.

   */
  CoAbstractStrategy(
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
  virtual ~CoAbstractStrategy() {};

  //! @name Access methods

  /*!
   \brief Return global slave node row map

   */
  const Teuchos::RCP<Epetra_Map> SlaveRowNodes() { return gsnoderowmap_; }

  /*!
   \brief Return global active node row map

   */
  virtual const Teuchos::RCP<Epetra_Map> ActiveRowNodes() { return gactivenodes_; }

  /*!
   \brief Return global slip node row map

   */
  const Teuchos::RCP<Epetra_Map> SlipRowNodes() { return gslipnodes_; }

  /*!
   \brief Return global slave dof row map

   */
  const Teuchos::RCP<Epetra_Map> SlaveRowDofs() { return gsdofrowmap_; }

  /*!
   \brief Return global active dof row map

   */
  virtual const Teuchos::RCP<Epetra_Map> ActiveRowDofs() { return gactivedofs_; }

  /*!
   \brief Return global master dof row map

   */
  const Teuchos::RCP<Epetra_Map> MasterRowDofs() { return gmdofrowmap_; }

  /*!
   \brief Return global slave dof row map

   */
  const Teuchos::RCP<Epetra_Map> SlaveMasterRowDofs() { return gsmdofrowmap_; }

  /*!
   \brief Return not redistributed global slave dof row map

   */
  const Teuchos::RCP<Epetra_Map> NotReDistSlaveRowDofs() { return pgsdofrowmap_; }

  /*!
   \brief Return not redistributed global master dof row map

   */
  const Teuchos::RCP<Epetra_Map> NotReDistMasterRowDofs() { return pgmdofrowmap_; }

  /*!
   \brief provides maps needed for contact/meshtying specific multigrid preconditioners

   */
  void CollectMapsForPreconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap,
      Teuchos::RCP<Epetra_Map>& SlaveDofMap,
      Teuchos::RCP<Epetra_Map>& InnerDofMap,
      Teuchos::RCP<Epetra_Map>& ActiveDofMap);

  /*!
   \brief Return Lagrange multiplier vector (t_n+1)

   */
  Teuchos::RCP<Epetra_Vector> LagrMult() { return z_; }

  /*!
   \brief Return old Lagrange multiplier vector (t_n)

   */
  Teuchos::RCP<Epetra_Vector> LagrMultOld() { return zold_; }

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
   \brief Return mortar matrix D

   */
  Teuchos::RCP<LINALG::SparseMatrix> DMatrix() { return dmatrix_; }

  /*!
   \brief Return mortar matrix M

   */
  Teuchos::RCP<LINALG::SparseMatrix> MMatrix() { return mmatrix_; }

  /*!
   \brief Return vector of normal contact stresses (t_n+1)

   */
  Teuchos::RCP<Epetra_Vector> ContactNorStress() { return stressnormal_; }

  /*!
   \brief Return vector of tangential contact stresses (t_n+1)

   */
  Teuchos::RCP<Epetra_Vector> ContactTanStress() { return stresstangential_; }

  /*!
   \brief Return required Integration time

   */
  double Inttime() { return inttime_; }
  ;

  /*!
   \brief Set integration time to zero

   */
  void Inttime_init() { inttime_ = 0.0; }
  ;

  /*!
   \brief Return current global contact status

   */
  bool IsInContact() { return isincontact_; }

  /*!
   \brief Return old global contact status (this time step)

   True if there has been contact in any nonlinear iteration
   step of the current time step.
   */
  bool WasInContact() { return wasincontact_; }

  /*!
   \brief Return old global contact status (last time step)

   True if there has been contact at the end of the last
   time step (last converged state)
   */
  bool WasInContactLastTimeStep() { return wasincontactlts_; }

  /*!
   \brief Return global self contact status

   Note that at the moment this only gives information about the
   POTENTIAL self contact of the global problem and not about
   an actual self contact occurring.

   TODO: automatically recognize ACTUAL self contact

   */
  bool& IsSelfContact() { return isselfcontact_; }

  /*!
   \brief Return global frictional status

   */
  bool Friction() { return friction_; }

  /*!
   \brief Return contact interfaces

   */
  const std::vector<Teuchos::RCP<CONTACT::CoInterface> > ContactInterfaces() { return interface_; }

  /*!
   \brief Get dual quadratic 3d slave element flag

   Returns TRUE if at least one higher-order 3d slave element with
   dual Lagrange mutliplier shape functions in any interface.

   */
  virtual const bool& Dualquadslave3d() { return dualquadslave3d_; }
  ;

  /*!
   \brief Return parallel redistribution status (yes or no)

  */
  bool ParRedist()
  {
    INPAR::MORTAR::ParRedist partype = DRT::INPUT::IntegralValue<
        INPAR::MORTAR::ParRedist>(Params(), "PARALLEL_REDIST");
    if (partype != INPAR::MORTAR::parredist_none)
      return true;
    else
      return false;
  }


  /*!
   \brief Return specific parallel redistribution status

   */
  INPAR::MORTAR::ParRedist WhichParRedist()
  {
    return DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(Params(),
        "PARALLEL_REDIST");
  }

  /*!
   \brief Return matrix T

   */
  virtual Teuchos::RCP<LINALG::SparseMatrix> TMatrix() { return Teuchos::null; }

  /*!
   \brief Return number of active nodes

   */
  int NumberOfActiveNodes()
  {
    if (gactivenodes_ != Teuchos::null)
      return gactivenodes_->NumGlobalElements();
    return 0;
  }

  /*!
   \brief Return number of frictional slip nodes

   */
  int NumberOfSlipNodes()
  {
    if (gslipnodes_ != Teuchos::null)
      return gslipnodes_->NumGlobalElements();
    return 0;
  }

  //@}

  //! @name Evaluation methods

  /*!
   \brief Redistribute all contact interfaces in parallel

   We hand in the current global displacement state so that
   a contact search can be performed and set state called.

   */
  void RedistributeContact(Teuchos::RCP<Epetra_Vector> dis);

  /*!
   \brief Global evaluation method called from time integrator

   */
  void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis,
      Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
      const int step, const int iter, bool predictor = false);

  /*!
   \brief Set current deformation state

   All interfaces are called to set the current deformation state
   (u, xspatial) in their nodes. Additionally, the new contact
   element areas are computed.

   \param statename (in): std::string defining which quantity to set (either "displacement" or "olddisplacement")
   \param vec (in): current global state of the quantity defined by statename

   */
  virtual void SetState(const std::string& statename,
      const Teuchos::RCP<Epetra_Vector> vec);

  /*!
   \brief Evaluate reference state

   for frictional contact we need history values (relative velocity) and
   therefore we store the nodal entries of mortar matrices (reference
   configuration) before the first time step

   */
  void EvaluateReferenceState(const Teuchos::RCP<Epetra_Vector> vec);

  void InitBinStrategyforTimestep(Teuchos::RCP<Epetra_Vector> vel);

  /*!
   \brief Evaluate matrix of nodal normals

   This is needed for energy-conserving time integration (Velocity-Update)
   */

  Teuchos::RCP<LINALG::SparseMatrix> EvaluateNormals(
      Teuchos::RCP<Epetra_Vector> dis);

  //@}

  //! @name Quantity control methods

  /*!
   \brief Get some nodal quantity globally and store into CoNodes

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
  void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type);

  /*!
   \brief Evaluate contact stresses in normal direction and tangential plane

   This is called at the end of each time or load step. It calculates
   the stress vector in normal direction and the stress vector in the
   tangential plane.

   */
  void OutputStresses();

  /*!
   \brief Get dirichlet B.C. status and store into CoNodes

   This is called once at the beginning of the simulation
   to set the D.B.C. status in each CNode.

   \param dbcmaps (in): MapExtractor carrying global dbc map

   */
  void StoreDirichletStatus(Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

  virtual void SetParentState(
      const std::string& statename,
      const Teuchos::RCP<Epetra_Vector> vec,
      const Teuchos::RCP<DRT::Discretization> dis)
  {/* standard contact methods don't need the corresponding bulk element */};

  /*!
   \brief Update contact at end of time step

   \param dis (in):  current displacements (-> old displacements)

   */
  void Update(Teuchos::RCP<Epetra_Vector> dis);

  /*!
   \brief Perform a write restart

   A write restart is initiated by the contact manager. However, the manager has no
   direct access to the nodal quantities. Hence, a portion of the restart has to be
   performed on the level of the contact algorithm, for short: here's the right place.

   */
  void DoWriteRestart(Teuchos::RCP<Epetra_Vector>& activetoggle,
      Teuchos::RCP<Epetra_Vector>& sliptoggle,
      Teuchos::RCP<Epetra_Vector>& weightedwear,
      Teuchos::RCP<Epetra_Vector>& realwear, bool forcedrestart = false);

  /*!
   \brief Perform a write restart

   A write restart is initiated by the contact manager. However, the manager has no
   direct access to the nodal quantities. Hence, all the restart action has to be
   performed on the level of the contact algorithm, for short: here's the right place.

   */
  void DoReadRestart(IO::DiscretizationReader& reader,
      Teuchos::RCP<Epetra_Vector> dis);

    //@}

    //! @name Empty functions if no augmented Lagrange formulation is used
    virtual void BuildGlobalAugActiveSet(const int it) {return;}
    virtual void UpdateStructuralRHS(Teuchos::RCP<Epetra_Vector>& feff) {return;}
    virtual void UpdateStructuralStiff(Teuchos::RCP<LINALG::SparseOperator>& kteff) {return;}
    virtual void AugForces(Epetra_Vector& augfs_lm, Epetra_Vector& augfs_g,
                           Epetra_Vector& augfm_lm, Epetra_Vector& augfm_g) {return;}
    //@}



  //! @name Debugging methods
  /*!
   \brief Compute interface forces and moments

   Compute current interface forces and moments at n+1-alphaf using current
   Lagrange multiplier values and current Mortar matrices D and M at n+1. When
   doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
   matrices of the last converged time / load step n (TR-like interpolation).

   \param output (in): flag indicating whether force output shall be written

   */
  virtual void InterfaceForces(bool output = false);

  /*!
   \brief Print interfaces

   */
  virtual void Print(std::ostream& os) const;

  /*!
   \brief Print current active set to screen for debugging purposes

   */
  virtual void PrintActiveSet();

  /*!
   \brief Visualize contact stuff with gmsh

   \param step (in): current time step index
   \param iter (in): current iteration index

   */
  void VisualizeGmsh(const int step, const int iter);

  //@}

  //! @name Purely virtual functions

  // All these functions are defined in one or more specific derived classes,
  // i.e CONTACT::CoLagrangeStrategy or CONTACT::CoPenaltyStrategy.
  // As the base class MORTAR::StrategyBase is always called from the control routine
  // (time integrator), these functions need to be defined purely virtual here.
  virtual bool ActiveSetSemiSmoothConverged() = 0;
  virtual bool ActiveSetConverged() = 0;
  virtual int ActiveSetSteps() = 0;
  virtual double ConstraintNorm() = 0;
  virtual void EvaluateContact(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
  virtual void EvaluateFriction(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
  virtual void EvaluateRelMovPredict() = 0;
  virtual double InitialPenalty() = 0;
  virtual void Initialize() = 0;
  virtual void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
  virtual void Recover(Teuchos::RCP<Epetra_Vector> disi) = 0;
  virtual void ResetActiveSet() = 0;
  virtual void ResetPenalty() = 0;
  virtual void ModifyPenalty() = 0;
  virtual void BuildSaddlePointSystem(Teuchos::RCP<LINALG::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<LINALG::MapExtractor> dbcmaps, int numiter,Teuchos::RCP<Epetra_Operator>& blockMat, Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) = 0;
  virtual void UpdateDisplacementsAndLMincrements(Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> blocksol) = 0;
  virtual void EvalConstrRHS() = 0;
  virtual void SaveReferenceState(const Teuchos::RCP<Epetra_Vector> dis) = 0;
  virtual void UpdateActiveSet() = 0;
  virtual void UpdateActiveSetSemiSmooth() = 0;
  virtual void UpdateUzawaAugmentedLagrange() = 0;
  virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;

  //@}

  //! @name Empty functions (meshtying)

  // All these functions only have functionality in meshtying simulations, thus they
  // are defined as empty here in the case of contact. They can be called from the
  // control routine (time integrator), whenever you like.

  void RedistributeMeshtying() {}
  void RestrictMeshtyingZone() {}
  void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) {}
  Teuchos::RCP<Epetra_Vector> MeshInitialization() {return Teuchos::null;};

  void MortarCoupling(const Teuchos::RCP<Epetra_Vector> dis) {}
  Teuchos::RCP<Epetra_Vector> LagrMultOldRescaled() {return Teuchos::null; }

protected:

  /*!
   \brief Evaluate contact

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
  void Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff,
      Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis);

  /*!
   \brief Evaluate relative movement of contact bodies

   This is for evaluating the relative movement of contact bodies. This
   can either be done with regarding the different movement of material points
   or regarding the change of mortar projection. The second possibility
   is definitely objective wheras the first possibility is objective
   only when the gap is zero.

   */

  void EvaluateRelMov();

  /*!
   \brief Initialize and evaluate Mortar stuff for the next Newton step

   This method first checks if we are dealing with self contact and updates
   the interface slave and master sets if so. Then it resets the global
   Mortar matrices D and M and the global gap vector g accordingly.

   The nodal quantites computed in InitEvalInterface() are then assembled
   to global matrices and vectors respectively. No setup of the global system
   is to be done here yet, so there is no need to pass in the effective
   stiffness K or the effective load vector f.

   */
  void InitMortar();
  void AssembleMortar();

  /*!
   \brief Initialize and evaluate interface for the next Newton step

   This method calls Initialize() on all contact interfaces, which
   resets all kind of nodal quantities like normal vector, weighted
   gap or Mortar and linearization maps. It then calls Evaluate() on
   all contact interfaces, which does all the geometric contact stuff.
   Concretely, this is an evaluation of all involved quantites at nodal
   level plus the setup of all corresponding linearizations.
   It includes the nodal normal calculations, contact search, projection
   and overlap detection, integration of the  Mortar terms D, M and of the
   weighted gap. Additionally, the linearizations of geometrical quantities
   (delta_n, delta_t, delta_D, delta_M) are calculated.

   */
  void InitEvalInterface();

  /*!
   \brief Setup this strategy object (maps, vectors, etc.)

   All global maps and vectors are initialized by collecting
   the necessary information from all interfaces. In the case
   of a parallel redistribution, this method is called again
   to re-setup the above mentioned quantities. In this case
   we set the input parameter redistributed=TRUE. Moreover,
   when called for the first time (in the constructor) this
   method is given the input parameter init=TRUE to account
   for initialization of the active set.

   */
  void Setup(bool redistributed, bool init);

  /*!
   \brief Update Mortar matrices D and M

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

  /*!
   \brief Store current (contact) nodal entries to old ones

   Contact nodes own their current entries and old ones (last converged
   state) from. p.e. the mortar matrices D and M. This function writes the
   current ones to the old ones.

   */
  void StoreToOld(MORTAR::StrategyBase::QuantityType type);

  /*!
   \brief Update global master and slave sets

   This becomes necessary for self contact simulations, because in a
   self contact interface master and slave status are assigned dynamically
   and thus the global maps change constantly.

   */
  void UpdateMasterSlaveSetsGlobal();

private:
  /*!
   \brief Calculate mean velocity of interface for binning range

   */
  void CalcMeanVelforBinning(Teuchos::RCP<Epetra_Vector> vel);

protected:

    // don't want cctor (= operator impossible anyway for abstract class)
    CoAbstractStrategy(const CoAbstractStrategy& old);

    std::vector<Teuchos::RCP<CONTACT::CoInterface> > interface_;  // contact interfaces

    std::vector<double>      tunbalance_;           // parallel unbalance factors (time) for current timestep
    std::vector<int>         eunbalance_;           // parallel unbalance factors (eles) for current timestep

    Teuchos::RCP<Epetra_Map> glmdofrowmap_;         // global Lagrange mult. dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsnoderowmap_;         // global slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmnoderowmap_;         // global master node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gsdofrowmap_;          // global slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gmdofrowmap_;          // global master dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gndofrowmap_;          // global internal dof row map
    Teuchos::RCP<Epetra_Map> gsmdofrowmap_;         // global slave and master dof row map (s+m map)
    Teuchos::RCP<Epetra_Map> gdisprowmap_;          // global displacement dof row map (s+m+n map)

    Teuchos::RCP<Epetra_Map> gactivenodes_;         // global active slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gactivedofs_;          // global active slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gactiven_;             // global dof row map of matrix N (of all interfaces)
    Teuchos::RCP<Epetra_Map> gactivet_;             // global dof row map of matrix T (of all interfaces)
    Teuchos::RCP<Epetra_Map> gslipnodes_;           // global slip slave node row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gslipdofs_;            // global slip slave dof row map (of all interfaces)
    Teuchos::RCP<Epetra_Map> gslipt_;               // global row map of matrix T for slip dofs (of all interfaces)

    Teuchos::RCP<Epetra_Map> pglmdofrowmap_;        // global LM dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsdofrowmap_;         // global slave dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgmdofrowmap_;         // global master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Map> pgsmdofrowmap_;        // global slave and master dof row map (before parallel redistribution)
    Teuchos::RCP<Epetra_Vector> pgsdirichtoggle_;   // global dirichlet toggle of all slave dofs (before parallel redistribution)

    std::vector<Teuchos::RCP<Epetra_Map> > initial_elecolmap_; // initla col ele map for binning strategy (s m)

    Teuchos::RCP<LINALG::SparseMatrix> dmatrix_;    // global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix> mmatrix_;    // global Mortar matrix M
    Teuchos::RCP<Epetra_Vector> g_;                 // global weighted gap vector g
    Teuchos::RCP<Epetra_Vector> tangrhs_;           // global tangential rhs vector (formulation with incremental z_)
    Teuchos::RCP<Epetra_Vector> inactiverhs_;       // gloabl inactive rhs vector (formulation with incremental z_ and saddle point system)
    Teuchos::RCP<Epetra_Vector> constrrhs_;         // global constraint rhs vector (only for saddlepoint problems)

    Teuchos::RCP<LINALG::SparseMatrix> lindmatrix_; // global Matrix LinD containing slave fc derivatives
    Teuchos::RCP<LINALG::SparseMatrix> linmmatrix_; // global Matrix LinM containing master fc derivatives

    Teuchos::RCP<LINALG::SparseMatrix> dold_;       // global Mortar matrix D (last end-point t_n)
    Teuchos::RCP<LINALG::SparseMatrix> mold_;       // global Mortar matrix D (last end-point t_n)

    Teuchos::RCP<Epetra_Vector> z_;                 // current vector of Lagrange multipliers at t_n+1
    Teuchos::RCP<Epetra_Vector> zold_;              // old vector of Lagrange multipliers at t_n
    Teuchos::RCP<Epetra_Vector> zincr_;             // Lagrange multiplier vector increment within SaddlePointSolve (this is NOT the increment of z_ between t_{n+1} and t_{n}!)
    Teuchos::RCP<Epetra_Vector> zuzawa_;            // vector of Lagrange multipliers from last Uzawa step

    Teuchos::RCP<Epetra_Vector> stressnormal_;      // vector of normal contact forces at t_n+1
    Teuchos::RCP<Epetra_Vector> stresstangential_;  // vector of tangential contact forces at t_n+1

    // global storage fields
    int step_;                // time step index
    int iter_;                // nonlinear (Newton) iteration index
    int iterls_;              // nonlinear (Line Search) iteration index

    bool isincontact_;        // flag indicating global contact status
    bool wasincontact_;       // flag indicating global contact status of this time step (history)
    bool wasincontactlts_;    // flag indicating global contact status of last time step
    bool isselfcontact_;      // flag indicating potential self contact
    bool friction_;           // flag for frictional contact
    bool regularized_;        // flag for regularized contact

    bool                                dualquadslave3d_; // flag indicating whether trafo should be applied
    Teuchos::RCP<LINALG::SparseMatrix>  trafo_;           // transformation matrix T for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix>  invtrafo_;        // inverse trafo matrix T^(-1) for dual quad 3D case
    Teuchos::RCP<LINALG::SparseMatrix>  dmatrixmod_;      // modified global Mortar matrix D
    Teuchos::RCP<LINALG::SparseMatrix>  doldmod_;         // modified global Mortar matrix Dold

    double inttime_;                                      // integration time
    std::vector<double > ivel_;                           // mean interface velocity

    INPAR::CONTACT::SolvingStrategy stype_;               // current used solving strategy
    INPAR::CONTACT::ConstraintDirection constr_direction_;// direction in which the contact constraints are formulated
}; // class CoAbstractStrategy
} // namespace CONTACT

// << operator
std::ostream& operator << (std::ostream& os, const CONTACT::CoAbstractStrategy& strategy);

#endif  // #ifndef CONTACT_STRATEGY_H
