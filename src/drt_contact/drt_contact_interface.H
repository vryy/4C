/*!----------------------------------------------------------------------
\file drt_contact_interface.H

\class CONTACT::Interface

\brief One contact interface

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

using namespace std;
using namespace Teuchos;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\class Interface

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
  virtual const Epetra_Map& SlaveRowNodes() const 
  { 
    if (Filled()) return *snoderowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master nodes (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& MasterRowNodes() const 
  { 
    if (Filled()) return *mnoderowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& SlaveColNodes() const 
  { 
    if (Filled()) return *snodecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master nodes (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& MasterColNodes() const 
  { 
    if (Filled()) return *mnodecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& SlaveFullNodes() const 
  { 
    if (Filled()) return *snodefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of master nodes (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& MasterFullNodes() const 
  { 
    if (Filled()) return *mnodefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of slave elements (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& SlaveRowElements() const 
  { 
    if (Filled()) return *selerowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master elements (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& MasterRowElements() const 
  { 
    if (Filled()) return *melerowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave elements (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& SlaveColElements() const 
  { 
    if (Filled()) return *selecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master elements (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& MasterColElements() const 
  { 
    if (Filled()) return *melecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of slave elements (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& SlaveFullElements() const 
  { 
    if (Filled()) return *selefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get fully overlapping column map of master elements (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& MasterFullElements() const 
  { 
    if (Filled()) return *melefullmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get row map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const Epetra_Map& SlaveRowDofs() const 
  { 
    if (Filled()) return *sdofrowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master dofs (Filled()==true is prerequisite)
  
  */
  virtual const Epetra_Map& MasterRowDofs() const 
  { 
    if (Filled()) return *mdofrowmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }  
  
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
  virtual void SetState(const string& statename, const RCP<Epetra_Vector> vec);

  /*!
  \brief Evaluate contact interface
      
  */
  virtual void Evaluate(Epetra_CrsMatrix& Dglobal,
  											Epetra_CrsMatrix& Mglobal,
  											Epetra_Vector& gglobal);
  
  /*!
  \brief Search algorithm for limitation of potentially contacting slave / master pairs
  
  */
  virtual bool EvaluateContactSearch();
  
  /*!
  \brief Integrate Mortar matrix D on slave element (2D)
          
  */
  virtual bool IntegrateSlave_2D(CONTACT::CElement& sele);
  
  /*!
  \brief Determine overlap and integrate slave / master pair (2D)
        
  */
  virtual bool IntegrateOverlap_2D(CONTACT::CElement& sele,
  													       CONTACT::CElement& mele);
  
  /*!
  \brief Assemble Mortar matrices and weighted gap
        
  */
  virtual void Assemble(Epetra_CrsMatrix& Dglobal,
  											Epetra_CrsMatrix& Mglobal,
  											Epetra_Vector& gglobal);
  
  /*!
  \brief Visualize node normals and projections with gmsh
          
  */
  virtual void VisualizeGmsh(const Epetra_SerialDenseMatrix& csegs);
  
  /*!
  \brief Return matrix of contact segment coordinates
  
  */
  virtual Epetra_SerialDenseMatrix& CSegs() {return contactsegs_;};

protected:

  // don't want = operator and cctor
  Interface operator = (const Interface& old);
  Interface(const Interface& old);


  int                           id_;            // unique interface id
  const Epetra_Comm&            comm_;          // communicator
  RCP<DRT::Discretization>      idiscret_;      // the discretization of the contact interface
  
  RCP<Epetra_Map>               snoderowmap_;   // row map of all slave nodes
  RCP<Epetra_Map>               snodecolmap_;   // column map of all slave nodes
  RCP<Epetra_Map>               snodefullmap_;  // fully overlapping column map of all slave nodes
  RCP<Epetra_Map>               mnoderowmap_;   // row map of all master nodes
  RCP<Epetra_Map>               mnodecolmap_;   // column map of all master nodes
  RCP<Epetra_Map>               mnodefullmap_;  // fully overlapping column map of all master nodes
  
  RCP<Epetra_Map>               selerowmap_;   // row map of all slave elements
  RCP<Epetra_Map>               selecolmap_;   // column map of all slave elements
  RCP<Epetra_Map>               selefullmap_;  // fully overlapping column map of all slave elements
  RCP<Epetra_Map>               melerowmap_;   // row map of all master elements
  RCP<Epetra_Map>               melecolmap_;   // column map of all master elements
  RCP<Epetra_Map>               melefullmap_;  // fully overlapping column map of all master elements
  
  RCP<Epetra_Map>               sdofrowmap_;   // row map of all slave dofs
  RCP<Epetra_Map>               mdofrowmap_;   // row map of all master dofs
  
  Epetra_SerialDenseMatrix		  contactsegs_;	 // contains coordinates of all contact segments (gmsh)
  int														counter_;			 // counts the number of calls to this interface (gmsh)
  
}; // class Interface
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::Interface& interface);









#endif  // #ifndef DRT_CONTACT_INTERFACE_H
#endif  // #ifdef CCADISCRET
