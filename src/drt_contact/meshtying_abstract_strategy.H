/*!----------------------------------------------------------------------
\file meshtying_abstract_strategy.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MESHTYING_ABSTRACT_STRATEGY_H
#define MESHTYING_ABSTRACT_STRATEGY_H

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_mortar.H"

using namespace std;
using namespace Teuchos;

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

// forward declarations

/*!
\brief Main abstract class for meshtying solution strategies
 
This is the templating abstract class for all meshyting solution algorithms.
Every solution algorithm has to fit into the set of functions and calls definied herein 
and has to be specified in a corresponding subclass defining the concrete algorithmic steps.

This class it itself derived from the MORTAR::StrategyBase class, which is an even
more abstract framework for any solution strategies involing mortar coupling.

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
    MtAbstractStrategy(DRT::Discretization& discret, RCP<Epetra_Map> problemrowmap,
                       Teuchos::ParameterList params, vector<RCP<MORTAR::MortarInterface> > interface,
                       int dim, RCP<Epetra_Comm> comm, double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~MtAbstractStrategy() {};
    
    //! @name Access methods
    
    /*!
    \brief Get the underlying problem discretization

    */
    DRT::Discretization& ProblemDiscret() const { return probdiscret_; }

    
    /*!
    \brief Return Lagrange multiplier vector (t_n+1)
     
    */
    RCP<Epetra_Vector> LagrMult() { return z_; }

    /*!
    \brief Return old Lagrange multiplier vector (t_n)
     
    */
    RCP<Epetra_Vector> LagrMultOld() { return zold_; }
    
    /*!
    \brief Return Lagrange multiplier vector from last Uzawa step
     
    */
    RCP<Epetra_Vector> LagrMultUzawa()
    {
      INPAR::CONTACT::SolvingStrategy soltype = Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");
      if (soltype != INPAR::CONTACT::solution_auglag) dserror("LagrMultUzawa() only defined for augmented strategy");
      return zuzawa_;
    }  
    
    /*!
	  \brief Return mortar matrix D

	  */
	  RCP<LINALG::SparseMatrix> DMatrix() { return dmatrix_; }
	  
	  /*!
	  \brief Return mortar matrix M

	  */
	  RCP<LINALG::SparseMatrix> MMatrix() { return mmatrix_; }

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
	  	INPAR::MORTAR::ParRedist partype = Teuchos::getIntegralValue<INPAR::MORTAR::ParRedist>(Params(),"PARALLEL_REDIST");
      if (partype!=INPAR::MORTAR::parredist_none) return true;
      else return false;
	  }

    /*!
    \brief Assemble global coordinate vector

    \param sidename (in): string indicating slave or master side
    \param ref (in): boolean indicating evaluation in reference configuration
    \param vec (in/out)):  empty global vetcor to be assembled to

    */
    void AssembleCoords(const string& sidename, bool ref, RCP<Epetra_Vector> vec);

    //@}

    //! @name Evaluation methods

    /*!
    \brief Redistribute all meshtying interfaces in parallel

    */
    void RedistributeMeshtying();

    /*!
    \brief Setup this strategy object (maps, vectors, etc.)

    All global maps and vectors are initialized by collecting
    the necessary information from all interfaces. In the case
    of a parallel redistribution, this method is called again
    to re-setup the above mentioned quantities. In this case
    the input parameter is set to TRUE.
    
    */
    void Setup(bool redistributed);

    /*!
    \brief Global evaluation method called from time integrator
    
    */
    void ApplyForceStiffCmt(RCP<Epetra_Vector> dis, RCP<LINALG::SparseOperator>& kt,
                            RCP<Epetra_Vector>& f, bool predictor = false)
    {
      // mortar initialization and evaluation
      SetState("displacement",dis);
      InitEvalInterface();
      InitEvalMortar();
      
      // evaluate relative movement for friction
      if (predictor) EvaluateRelMovPredict();
      else           EvaluateRelMov();
      
      // update active set
      UpdateActiveSetSemiSmooth();
      
      // apply contact forces and stiffness
      Initialize();
      Evaluate(kt,f,dis);
      InterfaceForces();
        
      return;
    }
    
    /*!
    \brief Set current deformation state
     
    All interfaces are called to set the current deformation state.
     
    \param statename (in): string defining which quantity to set (only "displacement" applicable)
    \param vec (in): current global state of the quantity defined by statename
     
    */
    void SetState(const string& statename, const RCP<Epetra_Vector> vec);
   
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
    virtual void MortarCoupling(const RCP<Epetra_Vector> dis);

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
    virtual void MeshInitialization(RCP<Epetra_Vector> Xslavemod);
    
    /*!
    \brief Evaluate contact
     
    This is just a tiny control routine, deciding which Evaluate-routine
    of those listed below is to be called (based on input-file information)
    Note that into ALL derived Evaluate() routines, a REFERENCE to the pointer
    on the effective stiffness matrix is handed in. This way, after building the
    new effective stiffness matrix with contact, we can simply let the pointer
    kteff point onto the new object. The same is true for the effective force
    vector feff. Be careful: kteff is of type RCP<LINALG::SparseOperator>&.
    
    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)
     
    */
    void Evaluate(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis);

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
    problems. Finally, "lmuzawa" addresses the LM update within an
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
    void StoreDirichletStatus(RCP<LINALG::MapExtractor> dbcmaps);
    
    /*!
    \brief Update meshtying at end of time step

    \param iter (in): step index for gmsh file creation
    \param dis (in):  current displacements (-> old displacements)

    */
    void Update(int istep, RCP<Epetra_Vector> dis);
      
    /*!
    \brief Perform a write restart
    
    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Different from writing a restart step, now
    all the restart action has to be performed on the level of the meshtying algorithm,
    for short: here's the right place. 
    
    */
    void DoReadRestart(IO::DiscretizationReader& reader, RCP<Epetra_Vector> dis);

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
    void Print(ostream& os) const;
    
    /*!
    \brief Print current active set to screen for debugging purposes
     
    */
    void PrintActiveSet();
    
    /*!
    \brief Visualize contact stuff with gmsh

    \param step (in): current time step index
    \param iter (in): current iteration index

    */
    void VisualizeGmsh(const int step, const int iter = 0);
    
    //@}

    //! @name Purely virtual functions
    
    // All these functions are defined in one or more specific derived classes,
    // i.e CONTACT::MeshtyingLagrangeStrategy or CONTACT::MeshtyingPenaltyStrategy.
    // As the base class MORTAR::StrategyBase is always called from the control routine
    // (time integrator), these functions need to be defined purely virtual here.
    
    virtual double ConstraintNorm() = 0;
    virtual void EvaluateMeshtying(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis) = 0;
    virtual void InitializeUzawa(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff) = 0;
    virtual double InitialPenalty() = 0;
    virtual void Recover(RCP<Epetra_Vector> disi) = 0;
    virtual void ResetPenalty() = 0;
    virtual void SaddlePointSolve(LINALG::Solver& solver, RCP<LINALG::SparseOperator> kdd, RCP<Epetra_Vector> fd, RCP<Epetra_Vector> sold, RCP<Epetra_Vector> dirichtoggle, int numiter) = 0;
    virtual void UpdateAugmentedLagrange() = 0;
    virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;
    
    //@}
    
    //! @name Empty functions (contact)
        
    // All these functions only have functionality in contact simulations, thus they
    // are defined as empty here in the case of meshtying. They can be called from the
    // control routine (time integrator), whenever you like.
    
    bool ActiveSetConverged() {return true;}
    bool ActiveSetSemiSmoothConverged() {return true;}
    bool Friction() {return false;}
    bool Wear() {return false;}
    RCP<Epetra_Vector> ContactNorStress() {return Teuchos::null;}
    RCP<Epetra_Vector> ContactTanStress() {return Teuchos::null;}
    RCP<Epetra_Vector> ContactWear() {return Teuchos::null;}
    void DoWriteRestart(RCP<Epetra_Vector>& activetoggle, RCP<Epetra_Vector>& sliptoggle, RCP<Epetra_Vector>& weightedwear) {}
    void InitEvalInterface() {}
    void InitEvalMortar() {}
    void Initialize() {}
    void OutputStresses() {}
    void OutputWear() {}
    void RedistributeContact(RCP<Epetra_Vector> dis) {}
    void ResetActiveSet() {}
    void ForceRefConfig() {}
    void SaveReferenceState(RCP<Epetra_Vector> dis) {}
    void UpdateActiveSet() {}
    void UpdateActiveSetSemiSmooth() {}
    void EvaluateReferenceState (int step,const RCP<Epetra_Vector> vec){}
    void EvaluateRelMov() {}
    void EvaluateRelMovPredict() {}
    
    //@}
        
  protected:

    // don't want cctor (= operator impossible anyway for abstract class)
    MtAbstractStrategy(const MtAbstractStrategy& old);

    DRT::Discretization&  probdiscret_;               // the underlying problem discretization
    vector<RCP<MORTAR::MortarInterface> > interface_; // contact interfaces
    
    RCP<Epetra_Map> glmdofrowmap_;         // global Lagrange mult. dof row map (of all interfaces)
    RCP<Epetra_Map> gsdofrowmap_;          // global slave dof row map (of all interfaces)
    RCP<Epetra_Map> gmdofrowmap_;          // global master dof row map (of all interfaces)
    RCP<Epetra_Map> gndofrowmap_;          // global internal dof row map (of all interfaces)
    RCP<Epetra_Map> gsmdofrowmap_;         // global slave and master dof row map (s+m map)
    RCP<Epetra_Map> gdisprowmap_;          // global displacement dof row map (s+m+n map)

    RCP<Epetra_Map> gsnoderowmap_;         // global slave node row map (of all interfaces)
    RCP<Epetra_Map> gmnoderowmap_;         // global master node row map (of all interfaces)

    RCP<Epetra_Map> pglmdofrowmap_;        // global LM dof row map (before parallel redistribution)
    RCP<Epetra_Map> pgsdofrowmap_;         // global slave dof row map (before parallel redistribution)
    RCP<Epetra_Map> pgmdofrowmap_;         // global master dof row map (before parallel redistribution)
    RCP<Epetra_Map> pgsmdofrowmap_;        // global slave and master dof row map (before parallel redistribution)

    RCP<LINALG::SparseMatrix> dmatrix_;    // global Mortar matrix D
    RCP<LINALG::SparseMatrix> mmatrix_;    // global Mortar matrix M
    RCP<Epetra_Vector> g_;                 // global weighted gap vector g

    RCP<Epetra_Vector> z_;                 // current vector of Lagrange multipliers at t_n+1
    RCP<Epetra_Vector> zold_;              // old vector of Lagrange multipliers at t_n
    RCP<Epetra_Vector> zuzawa_;            // vector of Lagrange multipliers from last Uzawa step

    bool                      dualquadslave3d_; // flag indicating whether trafo should be applied
    RCP<LINALG::SparseMatrix> trafo_;           // transformation matrix T for dual quad 3D case
    RCP<LINALG::SparseMatrix> invtrafo_;        // inverse trafo matrix T^(-1) for dual quad 3D case

}; // class MtAbstractStrategy
} // namespace CONTACT

// << operator
ostream& operator << (ostream& os, const CONTACT::MtAbstractStrategy& strategy);

#endif  // #ifndef MESHTYING_ABSTRACT_STRATEGY_H
#endif  // #ifdef CCADISCRET
