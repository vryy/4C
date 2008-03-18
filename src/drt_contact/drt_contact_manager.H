/*!----------------------------------------------------------------------
\file drt_contact_manager.H

\class CONTACT::Manager

\brief Main class to control all contact

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

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace CONTACT
{

// forward declarations

/*!
\class Manager

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
  inline RCP<Epetra_Vector>& GetContactForces() { return fc_; }
  
  /*!
  \brief Return Lagrange multiplier vector (mid-point)
    
  */
  inline RCP<Epetra_Vector>& LagrMult() { return z_; }
  
  /*!
  \brief Return old Lagrange multiplier vector (mid-point)
      
  */
  inline RCP<Epetra_Vector>& LagrMultOld() { return zold_; }
  
  /*!
  \brief Return Lagrange multiplier vector (end-point)
        
  */
  inline RCP<Epetra_Vector>& LagrMultEnd() { return zn_; }
  
  /*!
  \brief Return current global contact status
        
  */
  inline bool& IsInContact() { return isincontact_; }
  
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
      
  */
  void EvaluateMortar();
    
  /*!
  \brief Evaluate contact
    
  */
  void Evaluate(RCP<LINALG::SparseMatrix> kteff, RCP<Epetra_Vector> feff, int numiter);
  
  /*!
  \brief Transform displacement increment vector
  
  During the setup of the global tangent stiffness the so-called basis
  transformation introduced the incremental jump of displacements
  [d] = d_slave - inv(D) * M * d_master into the vector of unknowns
  instead of the incremental slave displacments. This step is reversed
  here by recovering the incremental slave displacements. But as we
  still need the incremental jump for the active set strategy, it is
  extracted at the beginning and stored in incrjump_
  
  */
  void RecoverDisp(RCP<Epetra_Vector>& disi);
  
  /*!
  \brief Update active set and check for convergence
  
  In this function we loop over all interfaces and then over all
  slave nodes to check, whether the assumption of them being active
  or inactive respectively has been correct. If a single node changes
  state, the active set is adapted accordingly and the convegence
  flag is kept on false -> another Newton iteration!
  
  */
  void UpdateActiveSet();
  
  /*!
  \brief Compute contact force vector
    
  */
  void ContactForces(RCP<Epetra_Vector>& fresm);
  
  

protected:

  // don't want = operator and cctor
  Manager operator = (const Manager& old);
  Manager(const Manager& old);


  DRT::Discretization&               discret_;        // the underlying problem discretization
  vector<RCP<CONTACT::Interface> >   interface_;      // contact interfaces
  
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
  
  RCP<LINALG::SparseMatrix>          nmatrix_;        // global Matrix N of active node normals
  RCP<LINALG::SparseMatrix>          tmatrix_;        // global Matrix T of active node tangents
  
  RCP<Epetra_Vector>                 z_;              // vector of Lagrange multipliers (mid-point)
  RCP<Epetra_Vector>                 zold_;           // vector of last time step Lagrange multipliers (mid-point)
  RCP<Epetra_Vector>                 zn_;             // vector of agrange multipliers (end-point)
  RCP<Epetra_Vector>                 incrjump_;       // incremental displacement jump (slave dofs)
  RCP<Epetra_Vector>                 fc_;             // vector of current contact forces
  bool                               isincontact_;    // flag indicatoing global contact status
  
  RCP<LINALG::SparseMatrix>          invd_;           // inverse of Mortar matrix D (needed for LM)
  RCP<LINALG::SparseMatrix>          ksn_;            // stiffness block K_sn (needed for LM)
  RCP<LINALG::SparseMatrix>          ksm_;            // stiffness block K_sm (needed for LM)
  RCP<LINALG::SparseMatrix>          kss_;            // stiffness block K_ss (needed for LM)


}; // class Manager
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::Manager& manager);









#endif  // #ifndef DRT_CONTACT_MANAGER_H
#endif  // #ifdef CCADISCRET
