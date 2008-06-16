/*!----------------------------------------------------------------------
\file drt_contact_interface.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_INTERFACE_H
#define DRT_CONTACT_INTERFACE_H

#include "Epetra_Comm.h"
#include "../drt_lib/drt_discret.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "../drt_lib/linalg_systemmatrix.H"

using namespace std;
using namespace Teuchos;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief One contact interface


\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class Interface
{
public:

  //! @name Enums and Friends
  //@}

  /*!
  \brief Standard constructor creating empty contact interface

  \param id (in): Unique interface id
  \param comm (in): A communicator object
  */
  explicit Interface(const int id, const Epetra_Comm& comm);

  /*!
  \brief Destructor

  */
  virtual ~Interface() {}

  /*!
  \brief Print this Interface

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return comm_;}

  //! @name Access methods
    
  /*!
  \brief Check whether interface was called FillComplete

  */
  virtual inline bool Filled() const { return idiscret_->Filled();}

  /*!
  \brief Get discretization of this interface

  */
  virtual inline DRT::Discretization& Discret() const { return *idiscret_; }

  /*!
  \brief Get row map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveRowNodes() const 
  { 
    if (Filled()) return snoderowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterRowNodes() const 
  { 
    if (Filled()) return mnoderowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColNodes() const 
  { 
    if (Filled()) return snodecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterColNodes() const 
  { 
    if (Filled()) return mnodecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveFullNodes() const 
  { 
    if (Filled()) return snodefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of master nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterFullNodes() const 
  { 
    if (Filled()) return mnodefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of slave elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveRowElements() const 
  { 
    if (Filled()) return selerowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterRowElements() const 
  { 
    if (Filled()) return melerowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColElements() const 
  { 
    if (Filled()) return selecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterColElements() const 
  { 
    if (Filled()) return melecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of slave elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveFullElements() const 
  { 
    if (Filled()) return selefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of master elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterFullElements() const 
  { 
    if (Filled()) return melefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveRowDofs() const 
  { 
    if (Filled()) return sdofrowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master dofs (Filled()==true is prerequisite)
  
  */
  virtual const RCP<Epetra_Map> MasterRowDofs() const 
  { 
    if (Filled()) return mdofrowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of active nodes (Filled()==true is prerequisite)
    
  */
  virtual const RCP<Epetra_Map> ActiveNodes() const 
  { 
    if (Filled()) return activenodes_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of active dofs (Filled()==true is prerequisite)
      
  */
  virtual const RCP<Epetra_Map> ActiveDofs() const 
  { 
    if (Filled()) return activedofs_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of matrix N (Filled()==true is prerequisite)
      
  */
  virtual const RCP<Epetra_Map> ActiveNDofs() const 
  { 
    if (Filled()) return activen_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
    
  /*!
  \brief Get row map of matrix T (Filled()==true is prerequisite)
      
  */
  virtual const RCP<Epetra_Map> ActiveTDofs() const 
  { 
    if (Filled()) return activet_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of slip nodes (Filled()==true is prerequisite)
    
  */
  virtual const RCP<Epetra_Map> SlipNodes() const 
  { 
    if (Filled()) return slipnodes_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of slip node dofs (Filled()==true is prerequisite)
      
  */
  virtual const RCP<Epetra_Map> SlipDofs() const 
  { 
    if (Filled()) return slipdofs_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of matrix T for slip nodes (Filled()==true is prerequisite)
      
  */
  virtual const RCP<Epetra_Map> SlipTDofs() const 
  { 
    if (Filled()) return slipt_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  //@}
  
  //! @name Evlauation methods
  
  /*!
  \brief Add a CONTACT::CNode to the interface (Filled()==true NOT prerequisite)
  
  \param cnode (in): rcp to a contact node
  
  \return Filled()==false
  */
  inline virtual void AddCNode(RCP<CONTACT::CNode> cnode) 
  { idiscret_->AddNode(cnode); return; }

  /*!
  \brief Add a CONTACT::CElement to the interface (Filled()==true is prerequisite)

  \param cele (in): rcp to a contact element

  \return Filled()==false
  */
  inline virtual void AddCElement(RCP<CONTACT::CElement> cele) 
  { idiscret_->AddElement(cele); return; }

  /*!
  \brief Finalize construction of contact interface

  The methods completes construction phase of a contact interface.
  It creates complete overlap ghosting of all elements and nodes
  on all processors that own at least one node or element on the interface.
  On return, row maps of nodes and elements are the same as on input but column maps
  are redundant on all participating processors.
  also creates row and column maps of slave and master nodes separately such that
  these can be easily accessed.
  
  \return Filled()==true
  */
  virtual void FillComplete();
  
  /*!
  \brief Initialize / reset interface for contact
  
  */
  virtual void Initialize();
  
  /*!
  \brief Set current deformation state
    
  */
  virtual void SetState(const string& statename,
                        const RCP<Epetra_Vector> vec);

  /*!
  \brief Evaluate contact interface
      
  */
  virtual void Evaluate();
  
  /*!
  \brief Search algorithm for limitation of potentially contacting slave / master pairs
  
  */
  virtual bool EvaluateContactSearch();
  
  /*!
  \brief Integrate Mortar matrix D on slave element (2D)
  
  This method is only called, if the flag CONTACTONEMORTARLOOP is
  not chosen. Then it integrates the Mortar matrix D on the current
  slave CElement sele (full integration from -1 to 1).
          
  */
  virtual bool IntegrateSlave2D(CONTACT::CElement& sele);
  
  /*!
  \brief Projection of slave / master pair (2D)
  
  This method projects the nodes of the slave CElement sele onto
  the master CElement mele and vice versa. The parameter hasproj
  stores a boolean variable for each of the 4 end nodes, indicating
  whether a feasible projection was found or not. The local element
  coordinates of the 4 projection points are stored in xiproj.
  
  */
  virtual bool Project2D(CONTACT::CElement& sele,
                         CONTACT::CElement& mele,
                         vector<bool>& hasproj,
                         vector<double>& xiproj);
  
  /*!
  \brief Detect overlap of slave / master pair (2D)
  
  This method evaluates the overlap of the current CElement pair
  sele / mele based on the projection status of the 4 end nodes
  (hasproj) and the coordinates of the projection points (xiproj).
  According to the detected overlap case, the integration limits
  are determined and written into xiproj and the overlap status
  is returned by a boolean variable.
  
  */
  virtual bool DetectOverlap2D(CONTACT::CElement& sele,
                               CONTACT::CElement& mele,
                               vector<bool>& hasproj,
                               vector<double>& xiproj);
  
  /*!
  \brief Integrate overlap of slave / master pair (2D)
  
  This method integrates the overlap of the current CElement
  pair sele / mele based on the integration limits (xiproj). The
  integration always includes the Mortar matrix M and the gap g.
  In some cases it also includes the Mortar matrix D (if the
  flag CONTACTONEMORTARLOOP is chosen) and the Mortar matrix
  Mmod (if the respective conditions are fulfilled).
        
  */
  virtual bool IntegrateOverlap2D(CONTACT::CElement& sele,
                                  CONTACT::CElement& mele,
                                  vector<double>& xiproj);
  
  /*!
  \brief Assemble Mortar matrices and weighted gap
        
  */
  virtual void AssembleDMG(LINALG::SparseMatrix& dglobal,
                           LINALG::SparseMatrix& mglobal,
                           Epetra_Vector& gglobal);
  
  /*!
  \brief Assemble matrices N / T containing nodal normals / tangents
  
  */
  virtual void AssembleNT(LINALG::SparseMatrix& nglobal,
                          LINALG::SparseMatrix& tglobal);
  
  /*!
  \brief Assemble matrix L and vector R for Tresca friction (slip)
    
  */
  virtual void AssembleTresca(LINALG::SparseMatrix& lglobal,
                              Epetra_Vector& rglobal,
                              double frbound);
    
  /*!
  \brief Assemble matrix S containing linearizations
  
  This method builds an algebraic form of the FULL linearization
  of the normal contact condition g~ = 0. Concretely, this
  includes assembling the linearizations of the slave side
  nodal normals and of the Mortar matrices D  and M.
  */
  virtual void AssembleS(LINALG::SparseMatrix& sglobal);
  
  /*!
  \brief Assemble matrix P containing linearizations
    
  This method builds an algebraic form of the FULL linearization
  of the tangential contact condition (frictionless). Concretely,
  this means assembling the linearization of the slave side
  nodal tangents and the current Lagrange multipliers.
  */
  virtual void AssembleP(LINALG::SparseMatrix& pglobal);
  
  /*!
  \brief Assemble matrices LinD, LinM containing linearizations
      
  This method builds an algebraic form of the FULL linearization
  of the contact force vector. Concretely, this means assembling
  the linearization of the Mortar matrices D and M and the
  current Lagrange multipliers.
  */
  virtual void AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                             LINALG::SparseMatrix& linmglobal);
    
  /*!
  \brief Initialize active set of this interface (t=0)
  
  \param initialcontact (in): true, if contact at t=0 assumed
        
  */
  virtual bool InitializeActiveSet(bool initialcontact);
    
  /*!
  \brief Build active set (nodes / dofs) of this interface
      
  */
  virtual bool BuildActiveSet();
    
  /*!
  \brief Split active dofs into N- and T-part
    
  */
  virtual bool SplitActiveDofs();
  
  //@}
  
  //! @name Visualization and Debugging methods
  
  /*!
  \brief Visualize node normals and projections with gmsh
          
  */
  virtual void VisualizeGmsh(const Epetra_SerialDenseMatrix& csegs);
  
  /*!
  \brief Return matrix of contact segment coordinates
  
  */
  virtual Epetra_SerialDenseMatrix& CSegs() {return contactsegs_;};
  
  /*!
  \brief Check normal/tangent derivatives with finite differences
    
  */
  virtual void FDCheckNormalDeriv();
  
  /*!
  \brief Check Mortar matrix D derivatives with finite differences
     
  */
  virtual void FDCheckMortarDDeriv();
  
  /*!
  \brief Check Mortar matrix M derivatives with finite differences
       
  */
  virtual void FDCheckMortarMDeriv();
  
  /*!
  \brief Check weighted gap g derivatives with finite differences
         
  */
  virtual void FDCheckGapDeriv();
  
  /*!
  \brief Check tangential LM derivatives with finite differences
           
  */
  virtual void FDCheckTangLMDeriv();
    
  //@}

protected:

  // don't want = operator and cctor
  Interface operator = (const Interface& old);
  Interface(const Interface& old);


  int                         id_;            // unique interface id
  const Epetra_Comm&          comm_;          // communicator
  RCP<DRT::Discretization>    idiscret_;      // the discretization of the contact interface
  
  RCP<Epetra_Map>             snoderowmap_;   // row map of all slave nodes
  RCP<Epetra_Map>             snodecolmap_;   // column map of all slave nodes
  RCP<Epetra_Map>             snodefullmap_;  // fully overlapping column map of all slave nodes
  RCP<Epetra_Map>             mnoderowmap_;   // row map of all master nodes
  RCP<Epetra_Map>             mnodecolmap_;   // column map of all master nodes
  RCP<Epetra_Map>             mnodefullmap_;  // fully overlapping column map of all master nodes
  
  RCP<Epetra_Map>             selerowmap_;    // row map of all slave elements
  RCP<Epetra_Map>             selecolmap_;    // column map of all slave elements
  RCP<Epetra_Map>             selefullmap_;   // fully overlapping column map of all slave elements
  RCP<Epetra_Map>             melerowmap_;    // row map of all master elements
  RCP<Epetra_Map>             melecolmap_;    // column map of all master elements
  RCP<Epetra_Map>             melefullmap_;   // fully overlapping column map of all master elements
  
  RCP<Epetra_Map>             sdofrowmap_;    // row map of all slave dofs
  RCP<Epetra_Map>             mdofrowmap_;    // row map of all master dofs
  
  RCP<Epetra_Map>             snodecolmapbound_;    // col map of slave nodes (+ boundary nodes)
  RCP<Epetra_Map>             snodefullmapbound_;   // full map of slave nodes (+ boundary nodes)
  RCP<Epetra_Map>             mnodefullmapnobound_; // full map of master nodes (- boundary nodes)
  
  RCP<Epetra_Map>             activenodes_;   // row map of all active slave nodes
  RCP<Epetra_Map>             activedofs_;    // row map of all active slave dofs
  RCP<Epetra_Map>             slipnodes_;     // row map of all slip slave nodes
  RCP<Epetra_Map>             slipdofs_;      // row map of all slip slave dofs  
  RCP<Epetra_Map>             activen_;       // row map of global N-matrix
  RCP<Epetra_Map>             activet_;       // row map of global T-matrix
  RCP<Epetra_Map>             slipt_;         // row map of part of T-matrix (slip nodes)
  
  Epetra_SerialDenseMatrix    contactsegs_;   // contains coordinates of all contact segments (gmsh)
  int                         counter_;       // counts the number of calls to this interface (gmsh)
  
}; // class Interface
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::Interface& interface);









#endif  // #ifndef DRT_CONTACT_INTERFACE_H
#endif  // #ifdef CCADISCRET
