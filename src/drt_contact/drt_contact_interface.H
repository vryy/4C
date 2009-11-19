/*!----------------------------------------------------------------------
\file drt_contact_interface.H

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
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_INTERFACE_H
#define DRT_CONTACT_INTERFACE_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_Comm.h"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "drt_contact_coupling3d_classes.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_contact/drt_contact_binarytree.H"
#include "../drt_contact/drt_contact_binarytree_self.H"
#include "../drt_inpar/inpar_contact.H"

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
  \brief Standard constructor creating empty contact interface
 
  This initializes the employed shape function set for lagrangian mutlipliers
  to a specific setting. Throughout the evaluation process, this set will be employed
  for the field of lagrangian multipliers.

  \param id (in): Unique interface id
  \param comm (in): A communicator object
  \param dim (in): Global problem dimension
  \param icontact (in): Global contact parameter list
  \param selfcontact (in): Flag for self contact status (default=false)
 
  */
  explicit Interface(const int id, const Epetra_Comm& comm, const int dim,
                     const Teuchos::ParameterList& icontact, bool selfcontact = false);
    
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
  \brief Get type of contact search algorithm

  */
  virtual inline const INPAR::CONTACT::ContactSearchAlgorithm SearchAlg()
  {return Teuchos::getIntegralValue<INPAR::CONTACT::ContactSearchAlgorithm>(icontact_,"SEARCH_ALGORITHM");};

  /*!
  \brief Get contact search algorithm parameter

  */
  virtual inline const double SearchParam()
  {return icontact_.get<double>("SEARCH_PARAM");};

  /*!
  \brief Get self contact status of this interface

  */
  virtual inline const bool& SelfContact() { return selfcontact_;}
  
  /*!
  \brief Get column map of all interface nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> OldColNodes() const
  {
    if (Filled()) return oldnodecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
  
  /*!
  \brief Get column map of all interface elements (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> OldColElements() const
  {
    if (Filled()) return oldelecolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }
      
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
  \brief Get column map of slave nodes + boundary nodes (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColNodesBound() const
  {
    if (Filled()) return snodecolmapbound_;
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
  \brief Get column map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveColDofs() const
  {
    if (Filled()) return sdofcolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of slave dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> SlaveFullDofs() const
  {
    if (Filled()) return sdoffullmap_;
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
  \brief Get column map of master dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterColDofs() const
  {
    if (Filled()) return mdofcolmap_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get fully overlapping column map of master dofs (Filled()==true is prerequisite)

  */
  virtual const RCP<Epetra_Map> MasterFullDofs() const
  {
    if (Filled()) return mdoffullmap_;
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
  \brief Update interface master and slave sets
  
  This update is usually only done ONCE during FillComplete() in
  the intialization phase and sets up the slave and master sets
  (elements, nodes, dofs) for the whole simulation. Yet, in the
  case of self contact the sets need to be updated again and again
  during simulation time, as the slave/master status is assigned dynamically.

  */
  virtual void UpdateMasterSlaveSets();
  
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

  This is the main routine of the Interface class, where nodal normals
  are computed, the contact search is performed, contact segments are
  set up and the entries of the mortar matrices D and M are integrated.

  */
  virtual void Evaluate();

  /*!
  \brief Brute force contact search algorithm for potentially contacting
         slave / master pairs (node-based algorithm)

  */
  virtual bool EvaluateContactSearch();

  /*!
  \brief Brute force contact search algorithm for potentially contacting
         slave / master pairs (element-based algorithm)

  */
  virtual void EvaluateContactSearchBruteForce(const double& eps);

  /*!
  \brief Binary tree contact search algorithm for potentially contacting
         slave / master pairs (element-based algorithm, SA T. Eberl)

  */
  virtual bool EvaluateContactSearchBinarytree();

  /*!
  \brief Integrate Mortar matrix D on slave element

  This method is only called, if the flag CONTACTONEMORTARLOOP is
  not chosen. Then it integrates the Mortar matrix D on the current
  slave CElement sele (full integration from -1 to 1).

  */
  virtual bool IntegrateSlave(CONTACT::CElement& sele);

  /*!
  \brief Integrate Mortar matrix M and gap g on slave/master overlap

  */
  virtual bool IntegrateCoupling(CONTACT::CElement& sele,
                                 CONTACT::CElement& mele);

  // (THIS IS A PURE FINITE DIFFERENCE VERSION!!!)
  virtual bool IntegrateCoupling(CONTACT::CElement& sele,
                                 CONTACT::CElement& mele,
                                 vector<vector<double> >& testv,
                                 bool printderiv=false);

  // (THIS IS A PURE FINITE DIFFERENCE VERSION!!!)
  virtual bool IntegrateCoupling(CONTACT::CElement& sele,
                                 CONTACT::CElement& mele,
                                 vector<vector<double> >& testgps,
                                 vector<vector<double> >& testgpm,
                                 vector<vector<double> >& testjs,
                                 vector<vector<double> >& testji,
                                 bool printderiv=false);

  /*!
  \brief Split CElements into IntElements for 3D quadratic coupling

  */
  virtual bool SplitIntElements(CONTACT::CElement& ele,
                                vector<RCP<CONTACT::IntElement> >& auxele);
  
  /*!
  \brief Integrate penalty sclaing factor kappa on slave element

  This method is only called, if a penalty strategy is applied. It is
  called ONCE at the beginning of the simulation and evaluates the
  penalty scaling factor kappa_j = int_{slave} (N_j) dslave. The
  correct interpolation N_j is chosen for any case (2D, 3D, linear
  quadratic, CONTACTPETROVGALERKIN...)
  
  TODO: maybe update kappa each time step?

  */
  virtual bool IntegrateKappaPenalty(CONTACT::CElement& sele);
    
  /*!
  \brief Evaluate relative movement (jump) of slave nodes
  
  In the case of frictional contact, an important geometric measure is
  the relative movement (jump) of the contacting bodies. Here, this is evaluated
  over change of mortar projection. Also, the directional derivatives are
  evaluated here.
  
  */ 
  virtual void EvaluateRelMov();
  
  /*!
  \brief Assemble jump vector into global jump vector 
  
  */
  virtual void AssembleRelMov(Epetra_Vector& jumpglobal);

  /*!
  \brief Assemble gap-computed lagrange multipliers and nodal linlambda derivatives into nodal quantities
  using the Macauley bracket
  
  When dealing with penalty methods, the lagrange multipliers are not independent variables anymore.
  Instead, they can be computed in terms of the weighted gap and the penalty parameter.
  This is done here so every node stores the correct lm and thus we integrate smoothly into the overlaying algorithm.
  
  Additionally, we use the performed loop over all nodes to store the nodal derivlambda_j matrix right there.
  
  As a result, the function notifies the calling routine if any negative gap was detected 
  and thus wether the interface is in contact or not. In consequence, after calling this routine
  from within the penalty strategy object, the contact status is known at a global level.
  
  Note: To be able to perform this computation, weighted gaps and normals have to be available within every node!
  Since this computation is done via Interface::Evaluate() in the Integrator class,
  these corresponding methods have to be called before AssembleMacauley()!
  
  */ 
  virtual void AssembleRegNormalForces(bool& localisincontact, bool& localactivesetchange);
  
  /*!
  \brief Assemble geometry dependent, tangential lagrange multipliers 
  and their derivatives in the penalty case
  */ 
  virtual void AssembleRegTangentForces();
    
  /*!
  \brief Assemble lagrange multipliers into global z vector (penalty strategy)
  
  */
  virtual void AssembleLM(Epetra_Vector& zglobal);
  
  /*
  \brief Assemble LM derivatives into global matrix (penalty strategy)
 
  */
  virtual void AssembleLinZ(LINALG::SparseMatrix& lambdaglobal);
  
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
                              double& frbound, double& ct);

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
  \brief Assemble matrix LinStick containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the tangential stick condition delta tg = 0. Concretely, this
  includes assembling the linearizations of the slave side
  nodal tangents and of the Mortar matrices D  and M.

  */
  virtual void AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
  		                          LINALG::SparseMatrix& linstickDISglobal,
                                Epetra_Vector& linstickRHSglobal);
  /*!
  \brief Assemble matrix LinSlip containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the tangential slip condition. Concretely, this
  includes assembling the linearizations of the slave side
  nodal tangents and of the Mortar matrices D  and M.

  */
  virtual void AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
  		                         LINALG::SparseMatrix& linslipDISglobal,
  		                         Epetra_Vector& linslipRHSglobal);

  /*!
  \brief Initialize active set of this interface (t=0)

  */
  virtual bool InitializeActiveSet();

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
  \brief Visualize contact stuff with gmsh

  */
  virtual void VisualizeGmsh(const Epetra_SerialDenseMatrix& csegs,
                             const int step, const int iter, const bool fric);

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

  /*!
  \brief Check stick condition derivatives with finite differences

  */
  virtual void FDCheckStickDeriv();

  /*!
  \brief Check slip condition derivatives with finite differences

  */
  virtual void FDCheckSlipDeriv();

  /*!
  \brief Check 3D coupling vertex derivatives with finite differences

  */
  virtual void FDCheckVertex3DDeriv(vector<vector<double> >& testv);

  /*!
  \brief Check 3D Gauss point derivatives and slave and intcell
  Jacobian derivatives with finite differences

  */
  virtual void FDCheckGP3DDeriv(vector<vector<double> >& testgps,
                                vector<vector<double> >& testgpm,
                                vector<vector<double> >& testjs,
                                vector<vector<double> >& testjm);

  /*!
  \brief Check penalty approach with finite differences

  */
  virtual void FDCheckPenaltyTracNor();  
  virtual void FDCheckPenaltyKTeff(RCP<LINALG::SparseMatrix> kc, RCP<Epetra_Vector> fc, RCP<LINALG::SparseMatrix> dold,
                                   RCP<LINALG::SparseMatrix> mold, RCP<Epetra_Vector> zold,
                                   RCP<Epetra_Map> problemrowmap, double alphaf);
  virtual void FDCheckPenaltyLinD(RCP<LINALG::SparseMatrix> lindmatrix, RCP<Epetra_Vector> reffcd);
  virtual void FDCheckPenaltyLinM(RCP<LINALG::SparseMatrix> linmmatrix, RCP<Epetra_Vector> reffcm);

  /*!
  \brief Check frictional penalty traction with finite differences

  */
  virtual void FDCheckPenaltyTracFric();
  
  /*!
  \brief Print shape function type (enum)

  */
  inline void PrintShapeFcn() { cout << shapefcn_ << endl; };
  
  //@}

protected:

  // don't want = operator and cctor
  Interface operator = (const Interface& old);
  Interface(const Interface& old);


  int                         id_;            // unique interface id
  const Epetra_Comm&          comm_;          // communicator
  RCP<Epetra_Comm>            lcomm_;         // interface local communicator
  map<int,int>                procmap_;       // mapping global -> local communicator PIDs

  RCP<DRT::Discretization>    idiscret_;      // the discretization of the contact interface
  int                         dim_;           // dimension of problem (2D or 3D)
  Teuchos::ParameterList      icontact_;      // containing contact input parameters of interface
  ShapeFcnType                shapefcn_;      // employed type of shape function set
  bool                        selfcontact_;   // flag indicating if this is a self contact interface
  
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

  RCP<Epetra_Map>             activenodes_;   // row map of all active slave nodes
  RCP<Epetra_Map>             activedofs_;    // row map of all active slave dofs
  RCP<Epetra_Map>             slipnodes_;     // row map of all slip slave nodes
  RCP<Epetra_Map>             slipdofs_;      // row map of all slip slave dofs
  RCP<Epetra_Map>             activen_;       // row map of global N-matrix
  RCP<Epetra_Map>             activet_;       // row map of global T-matrix
  RCP<Epetra_Map>             slipt_;         // row map of part of T-matrix (slip nodes)

  RCP<CONTACT::BinaryTree>       binarytree_;     // binary tree for contact search
  RCP<CONTACT::BinaryTreeSelf>   binarytreeself_; // binary tree for self contact search

  Epetra_SerialDenseMatrix       contactsegs_;    // contains coordinates of all contact segments (gmsh)

}; // class Interface
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::Interface& interface);


#endif  // #ifndef DRT_CONTACT_INTERFACE_H
#endif  // #ifdef CCADISCRET
