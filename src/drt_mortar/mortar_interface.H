/*!----------------------------------------------------------------------
\file mortar_interface.H

\brief One mortar coupling interface

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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_INTERFACE_H
#define MORTAR_INTERFACE_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "../drt_inpar/inpar_mortar.H"
#include "../drt_lib/drt_dserror.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
}

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations
class MortarNode;
class MortarElement;
class IntElement;
class BinaryTree;

/*!
\brief One mortar coupling interface

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class MortarInterface
{
public:

  /*!
  \brief Standard constructor creating empty mortar interface
 
  This initializes the employed shape function set for lagrangian mutlipliers
  to a specific setting. Throughout the evaluation process, this set will be employed
  for the field of lagrangian multipliers.

  \param id (in): Unique interface id
  \param comm (in): A communicator object
  \param dim (in): Global problem dimension
  \param imortar (in): Global contact parameter list
  \param redundant (in): type of redundancy in interface storage
 
  */
  MortarInterface(const int id, const Epetra_Comm& comm, const int dim,
                  const Teuchos::ParameterList& imortar,
                  INPAR::MORTAR::RedundantStorage redundant);
    
  /*!
  \brief Destructor

  */
  virtual ~MortarInterface() {}

  /*!
  \brief Print this Interface

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Check whether interface was called FillComplete

  */
  virtual bool Filled() const;

  /*!
  \brief Print parallel distribution of this Interface

  */
  virtual void PrintParallelDistribution(int index);

  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return comm_;}

  /*!
  \brief Get local interface communicator (Filled()==true is prerequisite)
         This Epetra_Comm only contains the participating processors of the interface

  */
  virtual const Epetra_Comm* lComm() const { return lcomm_.get();}

  /*!
  \brief Get global -> local processor map

  */
  virtual std::map<int,int>& Procmap(){ return procmap_;}

  //! @name Access methods

  /*!
  \brief Get discretization of this interface

  */
  virtual DRT::Discretization& Discret() const { return *idiscret_; }

  /*!
  \brief Get problem dimension

  Note that only 2D and 3D are possible here as this refers to the global
  problem dimension. On interface level this corresponds to 1D interfaces
  (dim_==2) and 2D interfaces (dim_==3)!

  */
  virtual const int& Dim() {return dim_;};

  /*!
  \brief Get interface contact parameter list

  */
  virtual Teuchos::ParameterList& IParams() {return imortar_;};

  /*!
  \brief Get quadratic 2d/3d slave element flag

  Returns TRUE if at least one higher-order 2d/3d slave element in interface.

  */
  virtual const bool& Quadslave() { return quadslave_; };

  /*!
  \brief Get redundant storage flag

  Returns status of redundancy in interface storage (nodes, elements).

  */
  virtual const INPAR::MORTAR::RedundantStorage& Redundant() { return redundant_; };

  /*!
  \brief Get type of search algorithm

  */
  virtual const INPAR::MORTAR::SearchAlgorithm SearchAlg()
  {return DRT::INPUT::IntegralValue<INPAR::MORTAR::SearchAlgorithm>(imortar_,"SEARCH_ALGORITHM");};

  /*!
  \brief Get search algorithm parameter

  */
  virtual const double SearchParam()
  {return imortar_.get<double>("SEARCH_PARAM");};
  
  /*!
  \brief Get column map of all interface nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> OldColNodes() const
  {
    if (Filled()) return oldnodecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get column map of all interface elements (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> OldColElements() const
  {
    if (Filled()) return oldelecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
      
  /*!
  \brief Get row map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveRowNodes() const
  {
    if (Filled()) return snoderowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterRowNodes() const
  {
    if (Filled()) return mnoderowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveColNodes() const
  {
    if (Filled()) return snodecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterColNodes() const
  {
    if (Filled()) return mnodecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of slave nodes + boundary nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveRowNodesBound() const
  {
    if (Filled()) return snoderowmapbound_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave nodes + boundary nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveColNodesBound() const
  {
    if (Filled()) return snodecolmapbound_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master nodes - boundary nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterRowNodesNoBound() const
  {
    if (Filled()) return mnoderowmapnobound_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master nodes - boundary nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterColNodesNoBound() const
  {
    if (Filled()) return mnodecolmapnobound_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of slave elements (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveRowElements() const
  {
    if (Filled()) return selerowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master elements (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterRowElements() const
  {
    if (Filled()) return melerowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave elements (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveColElements() const
  {
    if (Filled()) return selecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master elements (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterColElements() const
  {
    if (Filled()) return melecolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveRowDofs() const
  {
    if (Filled()) return sdofrowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlaveColDofs() const
  {
    if (Filled()) return sdofcolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of master dofs (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterRowDofs() const
  {
    if (Filled()) return mdofrowmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get column map of master dofs (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> MasterColDofs() const
  {
    if (Filled()) return mdofcolmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get map of Lagrange multiplier dofs (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> LagMultDofs() const
  {
    if (Filled()) return lmdofmap_;
    else dserror("MORTAR::MortarInterface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get maxmimum dof ID in global discretization

  */
  virtual const int& MaxDofGlobal()
  {
    if (maxdofglobal_ < 0) dserror("MaxDofGlobal not yet initialized");
    return maxdofglobal_;
  }
  
  //@}

  //! @name Evlauation methods

  /*!
  \brief Add a MORTAR::MortarNode to the interface (Filled()==true NOT prerequisite)

  \param mrtrnode (in): Teuchos::rcp to a mortar node

  \return Filled()==false

  */
  virtual void AddMortarNode(Teuchos::RCP<MORTAR::MortarNode> mrtrnode);

  /*!
  \brief Add a MORTAR::MortarElement to the interface (Filled()==true is prerequisite)

  \param mrtrele (in): Teuchos::rcp to a mortar element

  \return Filled()==false

  */
  virtual void AddMortarElement(Teuchos::RCP<MORTAR::MortarElement> mrtrele);

  /*!
  \brief Finalize construction of mortar interface

  The methods completes construction phase of a mortar interface.
  It creates complete overlap ghosting of all elements and nodes
  on all processors that own at least one node or element on the interface.
  On return, row maps of nodes and elements are the same as on input but column maps
  are redundant on all participating processors. Also creates row and column maps
  of slave and master nodes separately such that these can be easily accessed.

  \param maxdof (in): maximum dof ID in global discretization
                      (if default=0, then no lmdofrowmap is created)

  \return Filled()==true

  */
  virtual void FillComplete(int maxdof = 0);
  
  /*!
  \brief Redistribute interface among all procs

  When first creating a mortar interface, its parallel distribution
  is simply copied from the underlying problem discretization. This,
  of course, is not the optimal parallel distribution for evaluating
  the mortar coupling terms, as the interface ownership might be
  restricted to only very fee processors. Moreover, no parallel
  scalability can be achieved with this procedure, because adding
  processors to the problem discretization does not automatically
  mean adding processors to the interface discretization.

  Thus, an independent parallel distribution of the interface is
  desirable, which divides the interface among all available
  processors. Redistribute() is the method to achieve this.
  Internally, we call ParMetis to re-partition both slave and
  master side of the interface independently. This results in new
  "optimal" node/element maps of the interface discretization.
  Note that after Redistribute(), we must call FillComplete() again.

  */
  virtual void Redistribute();

  /*!
  \brief Create interface ghosting

  To ease our search algorithms we'll afford the luxury to ghost all nodes
  on all processors in the general mortar coupling framework. Currently, this
  still holds for both slave and master side of the interface.

  */
  virtual void CreateInterfaceGhosting();

  /*!
  \brief Create binary search tree

  The methods creates a binary tree object for efficient search.

  */
  virtual void CreateSearchTree();

  /*!
  \brief Update interface master and slave sets
  
  This update is usually only done ONCE in the initialization phase
  and sets up the slave and master sets (elements, nodes, dofs) for
  the whole simulation. Yet, in the case of self contact the sets
  need to be updated again and again during simulation time, as the
  slave/master status is assigned dynamically.

  */
  virtual void UpdateMasterSlaveSets();
  
  /*!
  \brief Restrict slave sets to actual meshtying zone

  This update is done ONCE in the initialization phase, when the
  slave surface only partially overlaps with the master surface.
  It restricts the slave sets (nodes, dofs, but NOT the elements)
  to the actual meshtying zone for the whole simulation.

  */
  virtual void RestrictSlaveSets();

  /*!
  \brief Update interface Lagrange multiplier sets

  This update is usually only done ONCE in the initialization phase
  and sets up the Lagrange mutliplier set (only dofs) for the whole
  simulation. Yet, in the case of self contact the sets need to be
  updated again and again during simulation time, as the slave/master
  status and thus the LM set is assigned dynamically.

  */
  virtual void UpdateLagMultSets(int offset_if);

  /*!
  \brief Initialize / reset mortar interface

  */
  virtual void Initialize();

  /*!
  \brief Set current deformation state

  */
  virtual void SetState(const string& statename,
                        const Teuchos::RCP<Epetra_Vector> vec);

  /*!
  \brief Set element areas

  */
  virtual void SetElementAreas();

  /*!
  \brief Evaluate mortar interface

  This is the main routine of the MortarInterface class, where nodal normals
  are computed, search is performed, mortar segments are set up and the entries
  of the mortar matrices D and M are integrated. If the boolean flag "nonlinear"
  is set to true (only for contact), then nonlinear mortar coupling is performed
  including evaluation of the weighted gap g and linearizations of all mortar quantities.
  
  \param nonlinear (in): nonlinear mortar coupling flag

  */
  virtual void Evaluate();

  /*!
  \brief Evaluate nodal normals

  */
  virtual void EvaluateNodalNormals();

  /*!
  \brief Export nodal normals

  This method exports / communicates the nodal normal vector from row to
  column map layout (needed for coupling evaluation). The reason behind this
  is that only the row nodes can correctly evaluate the nodal normal based on
  averaging of the adjacent element normals. A column node might not know of
  all adjacent elements and thus would compute wrong nodal normals by itself.

  */
  virtual void ExportNodalNormals();

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

  This method is outdated!!! We now obtain the D-matrix by
  integrating the mortar segments (just like the M-matrix entries).

  */
  virtual bool IntegrateSlave(MORTAR::MortarElement& sele);

  /*!
  \brief Integrate Mortar matrix M and gap g on slave/master overlaps
  
  */
  virtual bool IntegrateCoupling(MORTAR::MortarElement* sele,
                                 std::vector<MORTAR::MortarElement*> mele);

  /*!
  \brief Split MortarElements into IntElements for 3D quadratic coupling

  */
  virtual bool SplitIntElements(MORTAR::MortarElement& ele,
                                std::vector<Teuchos::RCP<MORTAR::IntElement> >& auxele);
  
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
  \brief Assemble matrix of normals N

  */
  virtual void AssembleNormals(LINALG::SparseMatrix& nglobal);

  /*!
  \brief Assemble transformation matrices T and T^(-1)

  These matrices need to be applied to the slave displacements
  in the cases of dual LM interpolation for tet10/hex20 meshes
  in 3D. Here, the displacement basis functions have been modified
  in order to assure positivity of the D matrix entries and at
  the same time biorthogonality. Thus, to scale back the modified
  discrete displacements \hat{d} to the nodal discrete displacements
  {d}, we have to apply the transformation matrix T and vice
  versa with the transformation matrix T^(-1).
  */
  virtual void AssembleTrafo(LINALG::SparseMatrix& trafo,
                             LINALG::SparseMatrix& invtrafo,
                             std::set<int>& donebefore);

  /*!
  \brief Detect actual meshtying zone (node by node)

  */
  virtual void DetectTiedSlaveNodes(int& founduntied);

  /*!
  \brief return integration time for current interface

  */
  virtual double Inttime(){return inttime_interface_;};

  //@}

  //! @name Visualization and Debugging methods

  /*!
  \brief Visualize mortar stuff with gmsh (may be overloaded)

  */
  virtual void VisualizeGmsh(const int step, const int iter);

  /*!
  \brief Print shape function type (enum)

  */
  void PrintShapeFcn() { std::cout << shapefcn_ << endl; };
  
  //@}

protected:

  // don't want = operator and cctor
  MortarInterface operator = (const MortarInterface& old);
  MortarInterface(const MortarInterface& old);

  int                                id_;            // unique interface id
  const Epetra_Comm&                 comm_;          // communicator
  Teuchos::RCP<Epetra_Comm>          lcomm_;         // interface local communicator
  std::map<int,int>                  procmap_;       // mapping global -> local communicator PIDs

  Teuchos::RCP<DRT::Discretization>  idiscret_;      // the discretization of the mortar interface
  int                                dim_;           // dimension of problem (2D or 3D)
  Teuchos::ParameterList             imortar_;      // containing contact input parameters of interface
  INPAR::MORTAR::ShapeFcn            shapefcn_;      // employed type of shape function set
  bool                               quadslave_;     // flag indicating quadratic 2d/3d slave elements
  INPAR::MORTAR::RedundantStorage    redundant_;     // employed type of redundancy in storage of interface
  
  Teuchos::RCP<Epetra_Map>           oldnodecolmap_; // column map of all interface nodes (overlap=1)
  Teuchos::RCP<Epetra_Map>           oldelecolmap_;  // column map of all interface elements (overlap=1)
    
  Teuchos::RCP<Epetra_Map>           snoderowmap_;   // row map of all slave nodes
  Teuchos::RCP<Epetra_Map>           snodecolmap_;   // column map of all slave nodes
  Teuchos::RCP<Epetra_Map>           mnoderowmap_;   // row map of all master nodes
  Teuchos::RCP<Epetra_Map>           mnodecolmap_;   // column map of all master nodes

  Teuchos::RCP<Epetra_Map>           snoderowmapbound_;    // row map of slave nodes (+ boundary nodes)
  Teuchos::RCP<Epetra_Map>           snodecolmapbound_;    // col map of slave nodes (+ boundary nodes)
  Teuchos::RCP<Epetra_Map>           mnoderowmapnobound_;  // row map of master nodes (- boundary nodes)
  Teuchos::RCP<Epetra_Map>           mnodecolmapnobound_;  // col map of master nodes (- boundary nodes)

  Teuchos::RCP<Epetra_Map>           selerowmap_;    // row map of all slave elements
  Teuchos::RCP<Epetra_Map>           selecolmap_;    // column map of all slave elements
  Teuchos::RCP<Epetra_Map>           melerowmap_;    // row map of all master elements
  Teuchos::RCP<Epetra_Map>           melecolmap_;    // column map of all master elements

  Teuchos::RCP<Epetra_Map>           sdofrowmap_;    // row map of all slave dofs
  Teuchos::RCP<Epetra_Map>           sdofcolmap_;    // column map of all slave dofs
  Teuchos::RCP<Epetra_Map>           mdofrowmap_;    // row map of all master dofs
  Teuchos::RCP<Epetra_Map>           mdofcolmap_;    // column map of all master dofs

  Teuchos::RCP<Epetra_Map>           lmdofmap_;      // row map of all Lagrange multiplier dofs
  int                                maxdofglobal_;  // maximum dof ID in global discretization
  
  Teuchos::RCP<MORTAR::BinaryTree>   binarytree_;    // binary searchtree

  double                             inttime_interface_;// integration time

}; // class MortarInterface
} // namespace MORTAR


// << operator
ostream& operator << (ostream& os, const MORTAR::MortarInterface& interface);


#endif  // #ifndef MORTAR_INTERFACE_H
