/*!----------------------------------------------------------------------
\file mortar_interface.H

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
#ifndef MORTAR_INTERFACE_H
#define MORTAR_INTERFACE_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_Comm.h"
#include "mortar_node.H"
#include "mortar_element.H"
#include "mortar_coupling3d_classes.H"
#include "mortar_binarytree.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_inpar/inpar_mortar.H"

using namespace std;
using namespace Teuchos;

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations

/*!
\brief One mortar coupling interface

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class MortarInterface
{
public:

  //! @name Enums and Friends
  
  /*!
  \brief Enum for shape function used for lagrange multiplier discretization
  
  */
  enum ShapeFcnType
  {
    StandardFunctions,
    DualFunctions,
    Undefined
  };
  
  //@}

  /*!
  \brief Standard constructor creating empty mortar interface
 
  This initializes the employed shape function set for lagrangian mutlipliers
  to a specific setting. Throughout the evaluation process, this set will be employed
  for the field of lagrangian multipliers.

  \param id (in): Unique interface id
  \param comm (in): A communicator object
  \param dim (in): Global problem dimension
  \param icontact (in): Global contact parameter list
 
  */
  explicit MortarInterface(const int id, const Epetra_Comm& comm, const int dim,
                           const Teuchos::ParameterList& icontact);
    
  /*!
  \brief Destructor

  */
  virtual ~MortarInterface() {}

  /*!
  \brief Print this Interface

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return comm_;}

  /*!
  \brief Get local interface communicator (Filled()==true is prerequisite)
         This Epetra_Comm only contains the participating processors of the interface

  */
  virtual const Epetra_Comm* lComm() const { return lcomm_.get();}


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
  \brief Get problem dimension

  Note that only 2D and 3D are possible here as this refers to the global
  problem dimension. On interface level this corresponds to 1D interfaces
  (dim_==2) and 2D interfaces (dim_==3)!

  */
  virtual inline const int& Dim() {return dim_;};

  /*!
  \brief Get interface contact parameter list

  */
  virtual inline Teuchos::ParameterList& IParams() {return icontact_;};

  /*!
  \brief Get type of search algorithm

  */
  virtual inline const INPAR::MORTAR::SearchAlgorithm SearchAlg()
  {return Teuchos::getIntegralValue<INPAR::MORTAR::SearchAlgorithm>(icontact_,"SEARCH_ALGORITHM");};

  /*!
  \brief Get search algorithm parameter

  */
  virtual inline const double SearchParam()
  {return icontact_.get<double>("SEARCH_PARAM");};
  
  /*!
  \brief Get column map of all interface nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> OldColNodes() const
  {
    if (Filled()) return oldnodecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get column map of all interface elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> OldColElements() const
  {
    if (Filled()) return oldelecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
      
  /*!
  \brief Get row map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveRowNodes() const
  {
    if (Filled()) return snoderowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterRowNodes() const
  {
    if (Filled()) return mnoderowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColNodes() const
  {
    if (Filled()) return snodecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave nodes + boundary nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColNodesBound() const
  {
    if (Filled()) return snodecolmapbound_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get column map of master nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterColNodes() const
  {
    if (Filled()) return mnodecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveFullNodes() const
  {
    if (Filled()) return snodefullmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of master nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterFullNodes() const
  {
    if (Filled()) return mnodefullmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of slave elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveRowElements() const
  {
    if (Filled()) return selerowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterRowElements() const
  {
    if (Filled()) return melerowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColElements() const
  {
    if (Filled()) return selecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterColElements() const
  {
    if (Filled()) return melecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of slave elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveFullElements() const
  {
    if (Filled()) return selefullmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of master elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterFullElements() const
  {
    if (Filled()) return melefullmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveRowDofs() const
  {
    if (Filled()) return sdofrowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColDofs() const
  {
    if (Filled()) return sdofcolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveFullDofs() const
  {
    if (Filled()) return sdoffullmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterRowDofs() const
  {
    if (Filled()) return mdofrowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterColDofs() const
  {
    if (Filled()) return mdofcolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of master dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterFullDofs() const
  {
    if (Filled()) return mdoffullmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  //@}

  //! @name Evlauation methods

  /*!
  \brief Add a MORTAR::MortarNode to the interface (Filled()==true NOT prerequisite)

  \param mrtrnode (in): rcp to a mortar node

  \return Filled()==false

  */
  inline virtual void AddMortarNode(RCP<MORTAR::MortarNode> mrtrnode)
  { idiscret_->AddNode(mrtrnode); return; }

  /*!
  \brief Add a MORTAR::MortarElement to the interface (Filled()==true is prerequisite)

  \param mrtrele (in): rcp to a mortar element

  \return Filled()==false

  */
  inline virtual void AddMortarElement(RCP<MORTAR::MortarElement> mrtrele)
  { idiscret_->AddElement(mrtrele); return; }

  /*!
  \brief Finalize construction of mortar interface

  The methods completes construction phase of a mortar interface.
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
  \brief Update interface master and slave sets
  
  This update is usually only done ONCE during FillComplete() in
  the intialization phase and sets up the slave and master sets
  (elements, nodes, dofs) for the whole simulation. Yet, in the
  case of self contact the sets need to be updated again and again
  during simulation time, as the slave/master status is assigned dynamically.

  */
  virtual void UpdateMasterSlaveSets();
  
  /*!
  \brief Initialize / reset mortar interface

  */
  virtual void Initialize();

  /*!
  \brief Set current deformation state

  */
  virtual void SetState(const string& statename,
                        const RCP<Epetra_Vector> vec);

  /*!
  \brief Evaluate mortar interface

  This is the main routine of the MortarInterface class, where nodal normals
  are computed, search is performed, mortar segments are set up and the entries
  of the mortar matrices D and M are integrated.

  */
  virtual void Evaluate();

  /*!
  \brief Brute force search algorithm for potentially coupling
         slave / master pairs (node-based algorithm)

  */
  virtual bool EvaluateSearch();

  /*!
  \brief Brute force search algorithm for potentially coupling
         slave / master pairs (element-based algorithm)

  */
  virtual void EvaluateSearchBruteForce(const double& eps);

  /*!
  \brief Binary tree search algorithm for potentially coupling
         slave / master pairs (element-based algorithm)

  */
  virtual bool EvaluateSearchBinarytree();

  /*!
  \brief Integrate Mortar matrix D on slave element

  This method is only called, if the flag MORTARONELOOP is
  not chosen. Then it integrates the Mortar matrix D on the current
  slave MortarElement sele (full integration from -1 to 1).

  */
  virtual bool IntegrateSlave(MORTAR::MortarElement& sele);

  /*!
  \brief Integrate Mortar matrix M and gap g on slave/master overlap

  */
  virtual bool IntegrateCoupling(MORTAR::MortarElement& sele,
                                 MORTAR::MortarElement& mele);

  /*!
  \brief Split MortarElements into IntElements for 3D quadratic coupling

  */
  virtual bool SplitIntElements(MORTAR::MortarElement& ele,
                                vector<RCP<MORTAR::IntElement> >& auxele);
  
  /*!
  \brief Assemble lagrange multipliers into global z vector (penalty strategy)
  
  */
  virtual void AssembleLM(Epetra_Vector& zglobal);
  
  /*!
  \brief Assemble Mortar matrices D and M

  */
  virtual void AssembleDM(LINALG::SparseMatrix& dglobal,
                          LINALG::SparseMatrix& mglobal);
  
  /*!
  \brief Assemble weighted gap g

  */
  virtual void AssembleG(Epetra_Vector& gglobal);

  //@}

  //! @name Visualization and Debugging methods

  /*!
  \brief Visualize contact stuff with gmsh

  */
  //virtual void VisualizeGmsh(const int step, const int iter, const bool fric);

  /*!
  \brief Print shape function type (enum)

  */
  inline void PrintShapeFcn() { cout << shapefcn_ << endl; };
  
  //@}

protected:

  // don't want = operator and cctor
  MortarInterface operator = (const MortarInterface& old);
  MortarInterface(const MortarInterface& old);


  int                         id_;            // unique interface id
  const Epetra_Comm&          comm_;          // communicator
  RCP<Epetra_Comm>            lcomm_;         // interface local communicator
  map<int,int>                procmap_;       // mapping global -> local communicator PIDs

  RCP<DRT::Discretization>    idiscret_;      // the discretization of the mortar interface
  int                         dim_;           // dimension of problem (2D or 3D)
  Teuchos::ParameterList      icontact_;      // containing contact input parameters of interface
  ShapeFcnType                shapefcn_;      // employed type of shape function set
  
  RCP<Epetra_Map>             oldnodecolmap_; // column map of all interface nodes (overlap=1)
  RCP<Epetra_Map>             oldelecolmap_;  // column map of all interface elements (overlap=1)
    
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
  RCP<Epetra_Map>             sdofcolmap_;    // column map of all slave dofs
  RCP<Epetra_Map>             sdoffullmap_;   // fully overlapping column map of all slave dofs
  RCP<Epetra_Map>             mdofrowmap_;    // row map of all master dofs
  RCP<Epetra_Map>             mdofcolmap_;    // column map of all master dofs
  RCP<Epetra_Map>             mdoffullmap_;   // fully overlapping column map of all master dofs

  RCP<Epetra_Map>             snodecolmapbound_;    // col map of slave nodes (+ boundary nodes)
  RCP<Epetra_Map>             snodefullmapbound_;   // full map of slave nodes (+ boundary nodes)
  RCP<Epetra_Map>             mnodefullmapnobound_; // full map of master nodes (- boundary nodes)

  RCP<MORTAR::BinaryTree>     binarytree_;    // binary searchtree

}; // class MortarInterface
} // namespace MORTAR


// << operator
ostream& operator << (ostream& os, const MORTAR::MortarInterface& interface);


#endif  // #ifndef MORTAR_INTERFACE_H
#endif  // #ifdef CCADISCRET
