/*!----------------------------------------------------------------------
\file drt_contact_abstract_strategy.H

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
#ifndef DRT_CONTACT_STRATEGY_H
#define DRT_CONTACT_STRATEGY_H

#include "Teuchos_RefCountPtr.hpp"
#include "drt_contact_interface.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_sparsematrix.H"

using namespace std;
using namespace Teuchos;

namespace CONTACT
{

// forward declarations

/*!
\brief Main abstract class for contact solution strategies
 
This is the templating abstract class for all contact solution algorithms.
Every solution algorithm has to fit into the set of functions and calls definied herein 
and has to be specified in a corresponding subclass defining the concrete algorithmic steps.

Refer also to the Semesterarbeit of Bernd Budich, 2009

\author popp (popp@lnm.mw.tum.de)
*/
class AbstractStrategy
{
  public:

    //! @name Enums and Friends
    enum QuantityType // types that can be called by StoreNodalQuantities()
    {
      lmcurrent,
      lmold,
      activeold,
      lmupdate,
      lmuzawa,
      jump
    };
    //@}

    /*!
    \brief Standard Constructor
     
    Creates the strategy object and initializes all global variables, including 
    all necessary Epetra_Maps and global vector and matrix quantities. 

    */
    explicit AbstractStrategy(RCP<Epetra_Map> problemrowmap, Teuchos::ParameterList params,
                              vector<RCP<CONTACT::Interface> > interface, int dim,
                              RCP<Epetra_Comm> comm, double alphaf);

    /*!
    \brief Destructor

    */
    virtual ~AbstractStrategy() {};
    
    //! @name Access methods
    
    /*!
    \brief Get contact parameter list

    */
    inline Teuchos::ParameterList& Params() { return scontact_; }
    
    /*!
    \brief Get problem dimension
     
    */
    inline const int& Dim() { return dim_; }

    /*!
    \brief Get Epetra communicator

    */
    inline const Epetra_Comm& Comm() const { return *comm_; }

    /*!
    \brief Return global slave dof row map
     
    */
    inline const RCP<Epetra_Map> SlaveRowDofs() { return gsdofrowmap_; }
    
    /*!
    \brief Return global slave node row map
     
    */
    inline const RCP<Epetra_Map> SlaveRowNodes() { return gsnoderowmap_; }

    /*!
    \brief Return convergence status of active set loop
     
    */
    inline bool& ActiveSetConverged() { return activesetconv_; }
     
    /*!
    \brief Return no. of active sets in this time/load step
     
    */
    inline int& ActiveSetSteps() { return activesetsteps_; }

    /*!
    \brief Return contact force vector
     
    */
    inline RCP<Epetra_Vector> GetContactForces() { return fc_; }

    /*!
    \brief Return Lagrange multiplier vector (t_n+1)
     
    */
    inline RCP<Epetra_Vector> LagrMult() { return z_; }

    /*!
    \brief Return old Lagrange multiplier vector (t_n)
     
    */
    inline RCP<Epetra_Vector> LagrMultOld() { return zold_; }
    
    /*!
    \brief Return Lagrange multiplier vector from last Uzawa step
     
    */
    inline RCP<Epetra_Vector> LagrMultUzawa()
    {
      INPAR::CONTACT::SolvingStrategy soltype = Teuchos::getIntegralValue<INPAR::CONTACT::SolvingStrategy>(Params(),"STRATEGY");
      if (soltype != INPAR::CONTACT::solution_auglag) dserror("LagrMultUzawa() only defined for augmented strategy");
      return zuzawa_;
    }

    /*!
    \brief Return jump vector (t_n+1)
     
    */
    inline RCP<Epetra_Vector> Jump() { return jump_; }

    /*!
    \brief Return current global contact status
     
    */
    inline bool& IsInContact() { return isincontact_; }

    //@}

    //! @name Evaluation methods

    /*!
    \brief Set current deformation state
     
    All interfaces are called to set the current deformation state
    (u, xspatial) in their nodes. Additionally, the new contact
    element areas are computed.
     
    \param statename (in): string defining which quantity to set (either "displacement" or "olddisplacement")
    \param vec (in): current global state of the quantity defined by statename
     
    */
    void SetState(const string& statename, const RCP<Epetra_Vector> vec);

    /*!
    \brief Initialize general Mortar related quantities for the next Newton step
     
    This method calls Initialize() on all contact interfaces, which
    resets all kind of nodal quantities like normal vector, weighted
    gap or Mortar and linearization maps. It also resets the global
    Mortar matrices D and M and the global gap vector g.
     
    Note: Only quantities common to all subsequent solving strategies are
    initialized here. In case they need additional Mortar variables,
    use the overloaded function call in the derived class and refer
    back to this function.
         
    */
    virtual void InitializeMortar();

    /*!
    \brief Initialize general contact variables for next Newton step
     
    This method initializes all global quantities that have to do with
    contact directly, or to be more precise with the active set.
    
    Note: Since this is specific for different solving strategies,
    this has to be implemented in derived classes. 
     
    */
    virtual void Initialize() = 0;

    /*!
    \brief Evaluate contact
     
    This is just a tiny control routine, deciding which Evaluate-routine
    of those listed below is to be called (based on input-file information)
     
    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)
     
    */
    void Evaluate(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff);
    
    /*!
    \brief Evaluate Mortar related quantites
     
    This method does all the geometric contact stuff that has not
    yet to do with the active / inactive node sets. Concretely it
    calls Evaluate() on all interfaces, which does a evaluation of
    all involved quantites at a nodal level plus
    the setup of all corresponding linearizations. 
    These nodal quantites are then assembled to global  matrices and vectors respectively. 
    No setup of the global system is to be done here yet, so there is no need to pass
    in the effective stiffness K or the effective load vector f.
     
    */
    virtual void EvaluateMortar() = 0;
    
    /*!
    \brief Evaluate contact
     
    This is the main routine of our contact algorithms on a global level.
    It contains the setup of the global linear system including contact.
    This makes modifications to the effective stiffness K and the effective
    load vector f necessary.
     
    Note: No basis transformation is applied to the displacement variables:
    delta_d_slave and delta_d_master are used directly as unknowns!
     
    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)
     
    */
    virtual void EvaluateContact(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff) = 0;

    /*!
    \brief Evaluate relative movement of contact bodies
     
    This is for evaluating the relative movement of contact bodies. This 
    can either be done with regarding the different movement of material points
    or regarding the change of mortar projection. The second possibility 
    is definitely objective wheras the first possibility is objective 
    only when the gap is zero.         
    
    */
    virtual void EvaluateRelMov(RCP<Epetra_Vector> disi) = 0;

    /*!
    \brief Evaluate frictional contact
     
    Note that this includes Perfect Stick, Tresca friction and Coulomb
    friction, of course! But it also includes the mesh tying case, as
    this is very similar to Perfect Stick from a conceptual point of view.
     
    \param kteff (in/out): effective stiffness matrix (without -> with contact)
    \param feff (in/out): effective residual / force vector (without -> with contact)
     
    */
    virtual void EvaluateFriction(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff) = 0;

    /*!
    \brief Recovery method
     
    Recovery is needed for the specialization of lagrangian multiplier approach.
    See LagrangeStrategy::Recover for details.
     
    \param disi (in): displacement increment of current Newton step solution
     
    */
    virtual void Recover(RCP<Epetra_Vector> disi) = 0;

    /*!
    \brief Update active set and check for convergence
     
    ActiveSet technology is needed for the specialization of lagrangian multiplier approach.
    See LAgrangeStrategy::UpdateActiveSet for details.
     
    */
    virtual void UpdateActiveSet() = 0;

    /*!
    \brief Update active set and check for convergence
     
    ActiveSet technology is needed for the specialization of lagrangian multiplier approach.
    See LAgrangeStrategy::UpdateActiveSet for details.
     
    */
    virtual void UpdateActiveSetSemiSmooth() = 0;

    //@}
    
    //! @name Quantity control methods
    
    /*!
    \brief Get some nodal quantity globally and store into CNodes

    The enum input parameter defines, which quantity is be updated.
    Currently the possibilities "lmold", "lmcurrent", "lmupdate", "jump"
    and "activeold" exist. Note that "lmold" means the converged value LM_n
    of the last time / load step, whereas "lmcurrent" adresses the current
    (not necessarily converged) value of the LM_n+1 and "jump" the current
    (not necessarily converged) value of d_n+1-d_n. "lmupdate" is a special
    option called only in Recover() after the update of the Lagr. multipliers.
    It basically does the same as "lmcurrent", but also checks for D.B.C.
    problems.

    \param type (in): enum defining which quantity to store into CNodes

    */
    void StoreNodalQuantities(AbstractStrategy::QuantityType type);

    /*!
    \brief Get dirichlet B.C. status and store into CNodes

    This is called once at the beginning of the simulation
    to set the D.B.C. status in each CNode.

    \param dbcmaps (in): MapExtractor carrying global dbc map

    */
    void StoreDirichletStatus(RCP<LINALG::MapExtractor> dbcmaps);
      
    /*!
    \brief Store current (contact) nodal entries of D and M to old ones
     
    Contact nodes own their current entries and old ones (last converged
    state) from the mortar matrices D and M. This function writes the
    current ones to the old ones.  
     
    */
    void StoreDMToNodes();

    /*!
    \brief Update Mortar matrices D and M
     
    The string input parameter defines in which direction the conversion
    is to be performed. Currently only the possibilities "old" and "current"
    exist, with "old" meaning the Mortar matrices of the last time / load step
    will be set to the current values D_n+1 and M_n+1 (this happens after
    completion of a time / load step!). The string "current" adresses the
    current Mortar matrices, which when called will be reset to the last
    converged values D_n and M_n (this happens in the predictor step when
    the active set has not yet converged!).
     
    \param state (in): string defining in which direction to convert D and M
     
    */
    void StoreDM(const string& state);
    
    /*!
    \brief Update contact at end of time step

    \param iter (in): step index for gmsh file creation

    */
    void Update(int istep);
      
    /*!
    \brief Perform a write restart
    
    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Hence, a portion of the restart has to be
    performed on the level of the contact algorithm, for short: here's the right place. 
    
    */
    void DoWriteRestart(RCP<Epetra_Vector>& activetoggle, RCP<Epetra_Vector>& sliptoggle);

    /*!
    \brief Perform a write restart
    
    A write restart is initiated by the contact manager. However, the manager has no
    direct access to the nodal quantities. Hence, a portion of the restart has to be
    performed on the level of the contact algorithm, for short: here's the right place. 
    
    */
    void DoReadRestart(RCP<Epetra_Vector> activetoggle, RCP<Epetra_Vector> sliptoggle,
                       RCP<Epetra_Vector> dis);

    //@}
    
    //! @name Debugging methods
    
    /*!
    \brief Compute contact force vector
     
    Compute current contact forces at n+1-alphaf using current Lagrange
    multiplier values and current Mortar matrices D and M at n+1. When
    doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
    matrices of the last converged time / load step n (TR-like interpolation).
     
    \param fresm (in): residual / force vector at state n+1 of current Newton step
    Note that we only pass in fresm for debugging purposes!
     
    */
    virtual void ContactForces(RCP<Epetra_Vector> fresm);

    /*!
    \brief Print interfaces

    */
    virtual void Print(ostream& os) const;
    
    /*!
    \brief Print current active set to screen for debugging purposes
     
    */
    virtual void PrintActiveSet();

    /*!
    \brief Visualize contact stuff with gmsh
     
    \param step (in): current time step index
    \param iter (in): current iteration index
    
    */
    void VisualizeGmsh(const int step, const int iter = 0);
    
    //@}

  protected:

    vector<RCP<CONTACT::Interface> > interface_; // contact interfaces
    Teuchos::ParameterList scontact_;            // containing contact input parameters
    int dim_;                                    // dimension of problem (2D or 3D)
    RCP<Epetra_Comm> comm_;                      // communicator
    double alphaf_;                              // Generalized-alpha parameter (0.0 for statics)

    RCP<Epetra_Map> problemrowmap_; // global problem dof row map (not only interfaces)
    RCP<Epetra_Map> gsnoderowmap_;  // global slave node row map (of all interfaces)
    RCP<Epetra_Map> gsdofrowmap_;   // global slave dof row map (of all interfaces)
    RCP<Epetra_Map> gmdofrowmap_;   // global master dof row map (of all interfaces)
    RCP<Epetra_Map> gndofrowmap_;   // global internal dof row map

    RCP<Epetra_Map> gactivenodes_;  // global active slave node row map (of all interfaces)
    RCP<Epetra_Map> gactivedofs_;   // global active slave dof row map (of all interfaces)
    RCP<Epetra_Map> gactiven_;      // global row map of matrix N (of all interfaces)
    RCP<Epetra_Map> gactivet_;      // global row map of matrix T (of all interfaces)
    RCP<Epetra_Map> gslipnodes_;    // global slip slave node row map (of all interfaces)
    RCP<Epetra_Map> gslipdofs_;     // global slip slave dof row map (of all interfaces)
    RCP<Epetra_Map> gslipt_;        // global row map of matrix T for slip dofs (of all interfaces)

    RCP<LINALG::SparseMatrix> dmatrix_; // global Mortar matrix D
    RCP<LINALG::SparseMatrix> mmatrix_; // global Mortar matrix M
    RCP<Epetra_Vector> g_;              // global weighted gap vector g

    RCP<LINALG::SparseMatrix> lindmatrix_; // global Matrix LinD containing slave fc derivatives
    RCP<LINALG::SparseMatrix> linmmatrix_; // global Matrix LinM containing master fc derivatives

    RCP<LINALG::SparseMatrix> dold_; // global Mortar matrix D (last end-point t_n)
    RCP<LINALG::SparseMatrix> mold_; // global Mortar matrix D (last end-point t_n)

    RCP<Epetra_Vector> z_;        // current vector of Lagrange multipliers at t_n+1
    RCP<Epetra_Vector> zold_;     // old vector of Lagrange multipliers at t_n
    RCP<Epetra_Vector> zuzawa_;   // vector of Lagrange multipliers from last Uzawa step
    RCP<Epetra_Vector> incrjump_; // incremental displacement jump (slave dofs)
    RCP<Epetra_Vector> jump_;     // displacement jump (slave dofs)
    RCP<Epetra_Vector> fc_;       // vector of current contact forces

    // global storage fields
    bool activesetconv_; // convergence flag for active set
    int activesetsteps_; // no. of active sets in this time/load step
    bool isincontact_;   // flag indicating global contact status

}; // class AbstractStrategy
} // namespace CONTACT

// << operator
ostream& operator << (ostream& os, const CONTACT::AbstractStrategy& strategy);

#endif  // #ifndef DRT_CONTACT_STRATEGY_H
#endif  // #ifdef CCADISCRET
