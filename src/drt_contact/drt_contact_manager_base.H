/*!----------------------------------------------------------------------
\file drt_contact_manager_base.H

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
<\pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_MANAGER_BASE_H
#define DRT_CONTACT_MANAGER_BASE_H

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
\brief Main class to control all contact

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class ManagerBase
{
public:

  //! @name Enums and Friends
  enum QuantityType // types that can be called by StoreNodalQuantities()
  {
    lmcurrent,
    lmold,
    lmupdate,
    jump,
    dirichlet
  };
  //@}

  /*!
  \brief Standard Constructor

  The base class constructor is empty.
  
  One needs a derived class for a concrete implementation of the contact
  manager into a given FE code environment (see e.g. drt_contact_manager.H
  and drt_contact_manager.cpp for the BACI implementation). This constructor
  then has to be fed with a discretization that is expected to carry at least
  two contact boundary conditions. It extracts all contact boundary conditions
  and constructs one or multiple contact interfaces from them and stores them.
  
  It calls CONTACT::Interface::FillComplete() on all contact interfaces which
  makes the nodes and elements of a contact interfaces redundant on all
  processors that either own a node or an element on the interfaces in the
  input discretization.
  
  \param alphaf (in): Generalized-alpha parameter (set to 0.0 by default)
  
  */
  explicit ManagerBase();

  /*!
  \brief Destructor

  */
  virtual ~ManagerBase() {};

  //! @name Access methods
  
  /*!
  \brief Get contact parameter list

  */
  inline Teuchos::ParameterList Params() const { return scontact_;};
    
  /*!
  \brief Get problem dimension
  
  */
  inline const int& Dim() {return dim_;};
  
  /*!
  \brief Get Epetra communicator

  */
  inline const Epetra_Comm& Comm() const { return *comm_;};
    
  /*!
  \brief Return global slave dof row map
  
  */
  inline const RCP<Epetra_Map> SlaveRowDofs() {return gsdofrowmap_; }
  
  /*!
  \brief Return global slave dof row map
    
  */
  inline const RCP<Epetra_Map> SlaveRowNodes() {return gsnoderowmap_; }
  
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
  
  \param statename (in): string defining which quantity to set (e.g. displacement)
  \param vec (in): current global state of the quantity defined by statename
  
  */
  void SetState(const string& statename, const RCP<Epetra_Vector> vec);
  
  /*!
  \brief Initialize Mortar stuff for next Newton step
  
  This method calls Initialize() on all contact interfaces, which
  resets all kind of nodal quantities like normal vector, weighted
  gap or Mortar and linearization maps. It also resets the global
  Mortar matrices D and M and the global gap vector g.
    
  */
  void InitializeMortar();
    
  /*!
  \brief Initialize contact for next Newton step
  
  This method initializes all global quantities that have to do with
  contact directly, or to be more precise with the active set. This
  includes the global normal / tangent matrices N and T, the global
  derivative matrices S and P and Tresca friction matrix L + vector r.
  
  */
  void Initialize();
  
  /*!
  \brief Evaluate Mortar matrices D,M and weighted gap
  
  This method does all the geometric contact stuff that has not
  yet to do with the active / inactive node sets. Concretely it
  calls Evaluate() on all interfaces, which does the nodal normal
  calculations, contact search, projection and overlap detection,
  integration of the Mortar terms D+M and of the weighted gap plus
  the setup of all linearizations of geometrical quantities (delta_n,
  delta_t, delta_D, delta_M). The nodal D, M and g values are then
  assembled into global matrices and vectors respectively. No setup of
  the global system is to be done here yet, so there is no need to pass
  in the effective stiffness K or the effective load vector f.
      
  */
  void EvaluateMortar();
    
  /*!
  \brief Evaluate contact
  
  This is just a tiny control routine, deciding which Evaluate-routine
  of those listed below is to be called (based on input-file information)
  
  \param kteff (in/out): effective stiffness matrix (without -> with contact)
  \param feff (in/out): effective residual / force vector (without -> with contact)
        
  */
  void Evaluate(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff);
  
  /*!
  \brief Evaluate Tresca friction (with basis transformation)
  
  \param kteff (in/out): effective stiffness matrix (without -> with contact)
  \param feff (in/out): effective residual / force vector (without -> with contact)
    
  */
  void EvaluateTrescaBasisTrafo(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff);
  
  /*!
  \brief Evaluate Tresca friction (without basis transformation)
  
  \param kteff (in/out): effective stiffness matrix (without -> with contact)
  \param feff (in/out): effective residual / force vector (without -> with contact)
      
  */
  void EvaluateTrescaNoBasisTrafo(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff);
   
  /*!
  \brief Evaluate contact (with basis transformation)
  
  This is the main routine of our contact algorithms on a global level.
  It contains the setup of the global linear system including contact.
  This makes modifications to the effective stiffness K and the effective
  load vector f necessary. As these modifications are usually applied
  to groups of dofs (slave, master, active etc.) some matrix and vector
  splitting and a lot of matrix-vector calculation is done! 
  A basis transformation is applied to the displacement variables: 
  [delta_d] := delta_d_slave - inv(D) * M * delta_d_master !
  
  \param kteff (in/out): effective stiffness matrix (without -> with contact)
  \param feff (in/out): effective residual / force vector (without -> with contact)
  
  */
  void EvaluateBasisTrafo(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff);
  
  /*!
  \brief Evaluate contact (without basis transformation)
  
  This is the main routine of our contact algorithms on a global level.
  It contains the setup of the global linear system including contact.
  This makes modifications to the effective stiffness K and the effective
  load vector f necessary. As these modifications are usually applied
  to groups of dofs (slave, master, active etc.) some matrix and vector
  splitting and a lot of matrix-vector calculation is done!
  No basis transformation is applied to the displacement variables:
  delta_d_slave and delta_d_master are used directly as unknowns !
  
  \param kteff (in/out): effective stiffness matrix (without -> with contact)
  \param feff (in/out): effective residual / force vector (without -> with contact)
       
  */
  void EvaluateNoBasisTrafo(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff);
    
  /*!
  \brief Recovery method (for displacements and Lagrange multipliers)
  
  Again, this is just a tiny control routine, deciding which Recover-routine
  of those listed below is to be called (based on input-file information)
  
  \param disi (in/out): displacement increment of current Newton step solution
  
  */
  void Recover(RCP<Epetra_Vector> disi);
    
  /*!
  \brief Recovery method (with basis transformation)
  
  During the setup of the global tangent stiffness the so-called basis
  transformation introduced the incremental jump of displacements
  [d] = d_slave - inv(D) * M * d_master into the vector of unknowns
  instead of the incremental slave displacments. This step is reversed
  here by recovering the incremental slave displacements.
  Moreover we recover the Lagrange multipliers here, which had been
  statically condensated during the setup of the global problem!
  
  \param disi (in/out): displacement increment of current Newton step solution
  
  */
  void RecoverBasisTrafo(RCP<Epetra_Vector> disi);
  
  /*!
  \brief Recovery method (without basis transformation)
  
  Since no basis transformation has been applied, we only recover
  the Lagrange multipliers here, which had been statically condensated
  during the setup of the global problem!
  
  \param disi (in): displacement increment of current Newton step solution
  
  */
  void RecoverNoBasisTrafo(RCP<Epetra_Vector> disi);
  
  /*!
  \brief Update active set and check for convergence
  
  In this function we loop over all interfaces and then over all
  slave nodes to check, whether the assumption of them being active
  or inactive respectively has been correct. If a single node changes
  state, the active set is adapted accordingly and the convegence
  flag is kept on false. Here we have the "standard" case of two
  nested iteration loops, and as a consequence this method is
  called AFTER convergence of the inner Newton iteration. If there
  is a change in the active set, another full Newton iteration has
  to be performed for the current time / load step.
  
  */
  void UpdateActiveSet();
  
  /*!
  \brief Update active set and check for convergence
  
  In this function we loop over all interfaces and then over all
  slave nodes to check, whether the assumption of them being active
  or inactive respectively has been correct. If a single node changes
  state, the active set is adapted accordingly and the convegence
  flag is kept on false. Here we have the semi-smooth Newton case
  with one combined iteration loop for active set search and large
  deformations. As a consequence this method is called AFTER each
  (not yet converged) Newton step. If there is a change in the active
  set or the residual and disp norm are still above their limits,
  another Newton step has to be performed.

  */
  void UpdateActiveSetSemiSmooth();
  
  /*!
  \brief Compute contact force vector
  
  Compute current contact forces at n+1-alphaf using current Lagrange
  multiplier values and current Mortar matrices D and M at n+1. When
  doing dynamics with alpha_f > 0, this also uses the old LM and Mortar
  matrices of the last converged time / load step n (TR-like interpolation).
  
  \param fresm (in): residual / force vector at state n+1 of current Newton step
                     Note that we only pass in fresm for debugging purposes!
  
  */
  void ContactForces(RCP<Epetra_Vector> fresm);
  
  /*!
  \brief Get some nodal quantity globally and store into CNodes
  
  The enum input parameter defines, which quantity is be updated.
  Currently the possibilities "lmold", "lmcurrent", "lmupdate", "jump"
  and "dirichlet" exist. Note that "lmold" means the converged value LM_n
  of the last time / load step, whereas "lmcurrent" adresses the current
  (not necessarily converged) value of the LM_n+1 and "jump" the current
  (not necessarily converged) value of d_n+1-d_n. "lmupdate" is a special
  option called only in Recover() after the update of the Lagr. multipliers.
  It basically does the same as "lmcurrent", but also checks for D.B.C.
  problems. Finally, "dirichlet" is called once at the beginning of the
  simulation to set the D.B.C. status in each CNode.
  
  \param type (in): enum defining which quantity to store into CNodes
  \param vec (in): optional, vector carrying global (only for "dirichlet")
  
  */
  void StoreNodalQuantities(ManagerBase::QuantityType type, RCP<Epetra_Vector> vec = null);
  
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
  \brief Print current active set to screen
  
  */
  void PrintActiveSet();
  
  /*!
  \brief Visualize contact stuff with gmsh
  
  \param step (in): current time step index
  \param iter (in): current iteration index
    
  */
  void VisualizeGmsh(const int step, const int iter = 0);
  
  //@}
  
  
protected:

  // don't want = operator and cctor
  ManagerBase operator = (const ManagerBase& old);
  ManagerBase(const ManagerBase& old);

  vector<RCP<CONTACT::Interface> >   interface_;       // contact interfaces
  Teuchos::ParameterList             scontact_;        // containing contact input parameters
  int                                dim_;             // dimension of problem (2D or 3D)
  double                             alphaf_;          // Generalized-alpha parameter (0.0 for statics)
  RCP<Epetra_Comm>                   comm_;            // communicator
  
  RCP<Epetra_Map>                    problemrowmap_;   // global problem dof row map (not only interfaces)
  RCP<Epetra_Map>                    gsnoderowmap_;    // global slave node row map (of all interfaces)
  RCP<Epetra_Map>                    gsdofrowmap_;     // global slave dof row map (of all interfaces)
  RCP<Epetra_Map>                    gmdofrowmap_;     // global master dof row map (of all interfaces)
  RCP<Epetra_Map>                    gndofrowmap_;     // global internal dof row map
  
  RCP<Epetra_Map>                    gactivenodes_;    // global active slave node row map (of all interfaces)
  RCP<Epetra_Map>                    gactivedofs_;     // global active slave dof row map (of all interfaces)
  RCP<Epetra_Map>                    gactiven_;        // global row map of matrix N (of all interfaces)
  RCP<Epetra_Map>                    gactivet_;        // global row map of matrix T (of all interfaces)
  RCP<Epetra_Map>                    gslipnodes_;      // global slip slave node row map (of all interfaces)
  RCP<Epetra_Map>                    gslipdofs_;       // global slip slave dof row map (of all interfaces)
  RCP<Epetra_Map>                    gslipt_;          // global row map of matrix T for slip dofs (of all interfaces)
  bool                               activesetconv_;   // convergence flag for active set
  int                                activesetsteps_;  // no. of active sets in this time/load step
    
  RCP<LINALG::SparseMatrix>          dmatrix_;         // global Mortar matrix D
  RCP<LINALG::SparseMatrix>          mmatrix_;         // global Mortar matrix M
  RCP<LINALG::SparseMatrix>          mhatmatrix_;      // product of global Mortar matrices inv(D)*M
  RCP<Epetra_Vector>                 g_;               // global weighted gap vector g
  
  RCP<LINALG::SparseMatrix>          mhata_;           // active part of mhatmatrix_ (needed for HüeberCheck)
  
  RCP<LINALG::SparseMatrix>          nmatrix_;         // global Matrix N containing active node normals
  RCP<LINALG::SparseMatrix>          tmatrix_;         // global Matrix T containing active node tangents
  RCP<LINALG::SparseMatrix>          smatrix_;         // global Matrix S containing normal+D+M derivatives
  RCP<LINALG::SparseMatrix>          pmatrix_;         // global Matrix P containing tangent derivatives
  RCP<LINALG::SparseMatrix>          lindmatrix_;      // global Matrix LinD containing slave fc derivatives
  RCP<LINALG::SparseMatrix>          linmmatrix_;      // global Matrix LinM containing master fc derivatives
  
  RCP<LINALG::SparseMatrix>          dold_;            // global Mortar matrix D (last end-point t_n)
  RCP<LINALG::SparseMatrix>          mold_;            // global Mortar matrix D (last end-point t_n)
  
  RCP<Epetra_Vector>                 z_;               // current vector of Lagrange multipliers at t_n+1
  RCP<Epetra_Vector>                 zold_;            // old vector of Lagrange multipliers at t_n
  RCP<Epetra_Vector>                 incrjump_;        // incremental displacement jump (slave dofs)
  RCP<Epetra_Vector>                 jump_;            // displacement jump (slave dofs)
  RCP<Epetra_Vector>                 fc_;              // vector of current contact forces
  bool                               isincontact_;     // flag indicating global contact status
  
  RCP<Epetra_Vector>                 fs_;              // slave side effective forces (needed for LM)
  RCP<LINALG::SparseMatrix>          invd_;            // inverse of Mortar matrix D (needed for LM)
  RCP<LINALG::SparseMatrix>          ksn_;             // stiffness block K_sn (needed for LM)
  RCP<LINALG::SparseMatrix>          ksm_;             // stiffness block K_sm (needed for LM)
  RCP<LINALG::SparseMatrix>          kss_;             // stiffness block K_ss (needed for LM)
  
  RCP<Epetra_Map>                    zigzagone_;       // active node set of last active set try
  RCP<Epetra_Map>                    zigzagtwo_;       // active node set of second-last active set try
  RCP<Epetra_Map>                    zigzagthree_;     // active node set of third-last active set try
  
  RCP<Epetra_Map>                    zigzagslipone_;   // slip node set of last slip set try
  RCP<Epetra_Map>                    zigzagsliptwo_;   // slip node set of second-last slip set try
  RCP<Epetra_Map>                    zigzagslipthree_; // slip node set of third-last slip set try
    
  RCP<LINALG::SparseMatrix>          lmatrix_;        // global matrix L for Tresca friction slip nodes
  RCP<Epetra_Vector>                 r_;              // r.h.s vector for Tresca friction slip nodes
  RCP<Epetra_Vector>                 rstick_;         // r.h.s vector for fricton stick nodes   
  
}; // class ManagerBase
} // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_MANAGER_BASE_H
#endif  // #ifdef CCADISCRET
