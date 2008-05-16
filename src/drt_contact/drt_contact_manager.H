/*!----------------------------------------------------------------------
\file drt_contact_manager.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_MANAGER_H
#define DRT_CONTACT_MANAGER_H

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "drt_contact_interface.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_lib/linalg_ana.H"

using namespace std;
using namespace Teuchos;

namespace CONTACT
{

// forward declarations

/*!
\brief Main class to control all contact


\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class Manager
{
public:

  //! @name Enums and Friends
  //@}

  /*!
  \brief Standard Constructor

  The constructor takes a discretization 
  that is expected to have an even number of contact boundary conditions.
  It extracts all contact boundary conditions and constructs one or multiple
  contact interfaces from them and stores them.
  
  It calls CONTACT::Interface::FillComplete() on all contact interfaces which
  makes the nodes and elements of a contact interfaces redundant on all
  processors that either own a node or an element on the interfaces in the
  input discretization.
  
  \param discret (in): A discretization containing contact boundary conditions
  */
  explicit Manager(DRT::Discretization& discret);


  /*!
  \brief Destructor

  */
  virtual ~Manager() {};


  /*!
  \brief Print this Manager

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Get Epetra communicator

  */
  inline const Epetra_Comm& Comm() const { return discret_.Comm();};

  //! @name Access methods
  
  /*!
  \brief Get discretization

  */
  inline DRT::Discretization& Discret() const { return discret_;};
  
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
  \brief Return contact force vector
  
  */
  inline RCP<Epetra_Vector> GetContactForces() { return fc_; }
  
  /*!
  \brief Return Lagrange multiplier vector (mid-point)
    
  */
  inline RCP<Epetra_Vector> LagrMult() { return z_; }
  
  /*!
  \brief Return old Lagrange multiplier vector (mid-point)
      
  */
  inline RCP<Epetra_Vector> LagrMultOld() { return zold_; }
  
  /*!
  \brief Return Lagrange multiplier vector (end-point)
        
  */
  inline RCP<Epetra_Vector> LagrMultEnd() { return zn_; }
  
  /*!
  \brief Return current global contact status
        
  */
  inline bool& IsInContact() { return isincontact_; }
  
  //@}
  
  //! @name Evaluation methods
  
  /*!
  \brief Read and check contact input parameters
  
  */
  bool ReadAndCheckInput();
  
  /*!
  \brief Write restart information for contact
    
  */
  RCP<Epetra_Vector> WriteRestart();
    
  /*!
  \brief Read restart information for contact
  
  */
  void ReadRestart(const RCP<Epetra_Vector> activetoggle);
  
  /*!
  \brief Initialize contact for next Newton step
  
  */
  void Initialize(int numiter);
  
  /*!
  \brief Set current deformation state
  
  */
  void SetState(const string& statename, const RCP<Epetra_Vector> vec);
  
  /*!
  \brief Evaluate contact Mortar matrices D,M only
    
  This is a light version of the Evaluate() method, which include the
  evluation of the current Mortar operators D and M only. No setup of
  the global system is to be done here, so there is no need to pass
  in the effective stiffness K or the effective load vector f.
      
  */
  void EvaluateMortar();
    
  /*!
  \brief Evaluate contact
        
  */
  void Evaluate(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff, int numiter);
  
  /*!
  \brief Evaluate Tresca friction (with basis transformation)
    
  */
  void EvaluateTrescaBasisTrafo(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff, int numiter);
   
  /*!
  \brief Evaluate contact (with basis transformation)
  
  This is the main routine of our contact algorithms on a global level.
  It includes the evaluation of the current Mortar operators D and M
  as well as the setup of the global linear system including contact.
  This makes modifications to the effective stiffness K and the effective
  load vector f necessary. As these modifications are usually applied
  to groups of dofs (slave, master, etc.) some matrix and vector
  splitting and a lot of matrix-vector calculation is done! 
   
  */
  void EvaluateBasisTrafo(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff, int numiter);
  
  /*!
  \brief Evaluate contact (without basis transformation)
      
  */
  void EvaluateNoBasisTrafo(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff, int numiter);
    
  /*!
  \brief Transform displacement increment vector
  
  */
  void RecoverDisp(RCP<Epetra_Vector> disi);
    
  /*!
  \brief Transform displacement increment vector (with basis transformation)
  
  During the setup of the global tangent stiffness the so-called basis
  transformation introduced the incremental jump of displacements
  [d] = d_slave - inv(D) * M * d_master into the vector of unknowns
  instead of the incremental slave displacments. This step is reversed
  here by recovering the incremental slave displacements. But as we
  still need the incremental jump for the active set strategy, it is
  extracted at the beginning and stored in incrjump_.
  Moreover we recover the Lagrange multipliers here, which had been
  statically condensated during the setup of the global problem!
  
  */
  void RecoverDispBasisTrafo(RCP<Epetra_Vector> disi);
  
  /*!
  \brief RecoverDisp (without basis transformation)
  
  */
  void RecoverDispNoBasisTrafo(RCP<Epetra_Vector> disi);
  
  /*!
  \brief Update active set and check for convergence
  
  In this function we loop over all interfaces and then over all
  slave nodes to check, whether the assumption of them being active
  or inactive respectively has been correct. If a single node changes
  state, the active set is adapted accordingly and the convegence
  flag is kept on false -> another Newton iteration!
  
  */
  void UpdateActiveSet(int numiteractive, RCP<Epetra_Vector> dism);
  
  /*!
  \brief Compute contact force vector
    
  */
  void ContactForces(RCP<Epetra_Vector> fresm);
  
  //@}
  
  
protected:

  // don't want = operator and cctor
  Manager operator = (const Manager& old);
  Manager(const Manager& old);


  DRT::Discretization&               discret_;        // the underlying problem discretization
  vector<RCP<CONTACT::Interface> >   interface_;      // contact interfaces
  Teuchos::ParameterList             scontact_;       // containing contact input parameters
  
  RCP<Epetra_Map>                    gsnoderowmap_;   // global slave node row map (of all interfaces)
  RCP<Epetra_Map>                    gsdofrowmap_;    // global slave dof row map (of all interfaces)
  RCP<Epetra_Map>                    gmdofrowmap_;    // global master dof row map (of all interfaces)
  RCP<Epetra_Map>                    gndofrowmap_;    // global internal dof row map
  
  RCP<Epetra_Map>                    gactivenodes_;   // global active slave node row map (of all interfaces)
  RCP<Epetra_Map>                    gactivedofs_;    // global active slave dof row map (of all interfaces)
  RCP<Epetra_Map>                    gactiven_;       // global row map of matrix N (of all interfaces)
  RCP<Epetra_Map>                    gactivet_;       // global row map of matrix T (of all interfaces)
  bool                               activesetconv_;  // convergence flag for active set
  
  RCP<LINALG::SparseMatrix>          dmatrix_;        // global Mortar matrix D
  RCP<LINALG::SparseMatrix>          mmatrix_;        // global Mortar matrix M
  RCP<LINALG::SparseMatrix>          mhatmatrix_;     // product of global Mortar matrices inv(D)*M
  RCP<Epetra_Vector>                 g_;              // global weighted gap vector g
  
  RCP<LINALG::SparseMatrix>          nmatrix_;        // global Matrix N containing active node normals
  RCP<LINALG::SparseMatrix>          tmatrix_;        // global Matrix T containing active node tangents
  RCP<LINALG::SparseMatrix>          smatrix_;        // global Matrix S containing active node normal derivatives
  
  RCP<Epetra_Vector>                 z_;              // vector of Lagrange multipliers (mid-point)
  RCP<Epetra_Vector>                 zold_;           // vector of last time step Lagrange multipliers (mid-point)
  RCP<Epetra_Vector>                 zn_;             // vector of Lagrange multipliers (end-point)
  RCP<Epetra_Vector>                 incrjump_;       // incremental displacement jump (slave dofs)
  RCP<Epetra_Vector>                 fc_;             // vector of current contact forces
  bool                               isincontact_;    // flag indicatoing global contact status
  
  RCP<LINALG::SparseMatrix>          invd_;           // inverse of Mortar matrix D (needed for LM)
  RCP<LINALG::SparseMatrix>          ksn_;            // stiffness block K_sn (needed for LM)
  RCP<LINALG::SparseMatrix>          ksm_;            // stiffness block K_sm (needed for LM)
  RCP<LINALG::SparseMatrix>          kss_;            // stiffness block K_ss (needed for LM)
  
  RCP<Epetra_Map>                    zigzagone_;      // active node set of last active set try
  RCP<Epetra_Map>                    zigzagtwo_;      // active node set of second-last active set try
  RCP<Epetra_Map>                    zigzagthree_;    // active node set of third-last active set try

  RCP<LINALG::SparseMatrix>          lmatrix_;        // global matrix L for Tresca friction slip nodes
  RCP<Epetra_Vector>                 r_;              // r.h.s vector for Tresca friction slip nodes

}; // class Manager
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::Manager& manager);









#endif  // #ifndef DRT_CONTACT_MANAGER_H
#endif  // #ifdef CCADISCRET
