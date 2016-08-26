/*---------------------------------------------------------------------*/
/*!
\file contact_interface.H

\brief One contact interface

\level 2

\maintainer Philipp Farah, Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_INTERFACE_H
#define CONTACT_INTERFACE_H

#include "../drt_mortar/mortar_interface.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_sparseoperator.H"
#include <Epetra_FEVector.h>


/*!
\brief CONTACT: namespace of the BACI contact and meshtying algorithms

*/
namespace ADAPTER
{
// forward declarations
  class Coupling;
}

namespace CONTACT
{

// forward declarations
class CoNode;
class CoElement;
class SelfBinaryTree;

/*!
\brief One contact interface

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class CoInterface : public MORTAR::MortarInterface
{
public:

  /*!
  \brief Standard constructor creating empty contact interface

  This initializes the employed shape function set for lagrangian mutlipliers
  to a specific setting. Throughout the evaluation process, this set will be employed
  for the field of lagrangian multipliers.

  \param id (in): Unique interface id
  \param comm (in): A communicator object
  \param dim (in): Global problem dimension
  \param icontact (in): Global contact parameter list
  \param selfcontact (in): Flag for self contact status
  \param redundant (in): type of redundancy in interface storage

  */
  CoInterface(const int id, const Epetra_Comm& comm, const int dim,
              const Teuchos::ParameterList& icontact,
              bool selfcontact, INPAR::MORTAR::RedundantStorage redundant);

  /*!
  \brief Destructor

  */
  virtual ~CoInterface() {}


  /*!
  \brief Print this Interface

  */
  virtual void Print(std::ostream& os) const;

  //! @name Access methods

  /*!
  \brief Get self contact status of this interface

  */
  virtual const bool& SelfContact() {return selfcontact_;}

  /*!
  \brief Get row map of active nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> ActiveNodes() const
  {
    if (Filled()) return activenodes_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }


  /*!
  \brief Get row map of active dofs (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> ActiveDofs() const
  {
    if (Filled()) return activedofs_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of matrix N (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> ActiveNDofs() const
  {
    if (Filled()) return activen_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of matrix T (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> ActiveTDofs() const
  {
    if (Filled()) return activet_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of slip nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlipNodes() const
  {
    if (Filled()) return slipnodes_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of slip node dofs (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlipDofs() const
  {
    if (Filled()) return slipdofs_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of matrix T for slip nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SlipTDofs() const
  {
    if (Filled()) return slipt_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of nonsmooth node (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> NonSmoothNodes() const
  {
    if (Filled()) return nonsmoothnodes_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }

  /*!
  \brief Get row map of smooth nodes (Filled()==true is prerequisite)

  */
  virtual const Teuchos::RCP<Epetra_Map> SmoothNodes() const
  {
    if (Filled()) return smoothnodes_;
    else dserror("CONTACT::Interface::FillComplete was not called");
    exit(EXIT_FAILURE); // calm down the compiler
  }



  /*!
  \brief Get number of slave / master pairs of this interface (proc local)

  */
  virtual const int& SlaveMasterPairs() { return smpairs_;}

  /*!
  \brief Get number of slave / master integration pairs of this interface (proc local)

  */
  virtual const int& SlaveMasterIntPairs() { return smintpairs_;}

  /*!
  \brief Get number of integration cells of this interface (proc local)

  */
  virtual const int& IntegrationCells() { return intcells_;}

  //@}

  //! @name Evlauation methods

  /*!
  \brief Add a CONTACT::CoNode to the interface (Filled()==true NOT prerequisite)

  \param cnode (in): Teuchos::rcp to a contact node

  \return Filled()==false

  */
  virtual void AddCoNode(Teuchos::RCP<CONTACT::CoNode> cnode);

  /*!
  \brief Add a CONTACT::CoElement to the interface (Filled()==true is prerequisite)

  \param cele (in): Teuchos::rcp to a contact element

  \return Filled()==false

  */
  virtual void AddCoElement(Teuchos::RCP<CONTACT::CoElement> cele);

  /*!
  \brief Redistribute contact interface among all procs

  Derived version!

  When first creating a contact interface, its parallel distribution
  is simply copied from the underlying problem discretization. This,
  of course, is not the optimal parallel distribution for evaluating
  the contact coupling terms, as the interface ownership might be
  restricted to only very few processors. Moreover, no parallel
  scalability can be achieved with this procedure, because adding
  processors to the problem discretization does not automatically
  mean adding processors to the interface discretization.

  Thus, an independent parallel distribution of the interface is
  desirable, which divides the interface among all available
  processors. Redistribute() is the method to achieve this.
  Moreover, for contact problems we have to account for the fact
  that only parts of the slave surface actually need to evaluate
  contact terms (those parts that are "close" to the master side).

  Internally, we call ZOLTAN to re-partition the contact interfacs
  in three independent parts: (1) close slave part, (2) non-close
  slave part, (3) master part. This results in new "optimal" node/element
  maps of the interface discretization. Note that after Redistribute(),
  we must call FillComplete() again. Note also that for contact
  simulations Redistribute() might be called dynamically again and
  again to account for changes of the contact zone.

  Two special cases are treated seperately: First, if ALL slave
  elements of the interface have some "close" neighbors, we do not
  need to distinguish the two different slave parts. Thus, we
  simply call the base class method Redistribute() also used for
  meshtying. Second, if NO slave element of the interface has any
  "close" neighbors, we do not need to redistribute at all. This
  is indicated by returning with a boolean return value FALSE.

  \param index (in): global index of current interface
  \return flag indicating whether redistribution has been performed

  */
  bool Redistribute(int index);

  void RoundRobinEvaluate();

  void RoundRobinChangeOwnership();

  void RoundRobinDetectGhosting();

  void RoundRobinExtendGhosting(bool firstevaluation);

  /*!
  \brief Collect data concering load balance and parallel distribution

  */
  void CollectDistributionData(int& loadele, int& crowele);

  /*!
  \brief Create binary search tree

  The method creates a binary tree object for efficient search. This is
  an overloaded method specific for contact, as in this case we have to
  consider the possibility of SELF-contact.

  Derived version!

  */
  void CreateSearchTree();

  /*!
  \brief Initialize / reset interface for contact

  Derived version!

  */
  void Initialize();

  /*!
  \brief Set element areas

  Derived version!

  */
  void SetElementAreas();

  /*!
  \brief Export nodal normals

  This method exports / communicates the nodal normal vector and all
  associated information (nodal tangent vectors, normal and tangent
  linearizations) from row to column map layout.

  Derived version!

  */
  void ExportNodalNormals();

  /*!
  \brief Binary tree search algorithm for potentially coupling slave /
  master pairs (element-based algorithm) including self-contact

  Derived version!

  */
  bool EvaluateSearchBinarytree();

  /*!
  \brief Integrate Mortar matrix M and gap g on slave/master overlaps

  Derived verision!

  */
  virtual bool MortarCoupling(MORTAR::MortarElement* sele,
      std::vector<MORTAR::MortarElement*> mele,
      const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

  /*!
  \brief evaluate coupling terms for nts coupling + lin

  */
  void EvaluateNTS();

  /*!
  \brief evaluate coupling terms for lts coupling + lin

  */
  void EvaluateLTS();

  /*!
  \brief evaluate coupling terms for ltl coupling + lin

  */
  void EvaluateLTL();

  /*!
  \brief evaluate coupling terms for stl coupling + lin

  */
  void EvaluateSTL();

  /*!
  \brief Integrate penalty sclaing factor kappa on slave element

  This method is only called, if a penalty strategy is applied. It is
  called ONCE at the beginning of the simulation and evaluates the
  penalty scaling factor kappa_j = int_{slave} (N_j) dslave. The
  correct interpolation N_j is chosen for any case (2D, 3D, linear
  quadratic, piecewise linear...)

  TODO: maybe update kappa each time step?

  */
  virtual bool IntegrateKappaPenalty(CONTACT::CoElement& sele);

  /*!
  \brief Evaluate relative movement (jump) of slave nodes

  In the case of frictional contact, an important geometric measure is
  the relative movement (jump) of the contacting bodies. Here, this is evaluated
  over change of mortar projection. Also, the directional derivatives are
  evaluated here.

  */
  virtual void EvaluateRelMov(const Teuchos::RCP<Epetra_Vector> xsmod,
                              const Teuchos::RCP<LINALG::SparseMatrix> dmatrixmod,
                              const Teuchos::RCP<LINALG::SparseMatrix> doldmod);

  /*!
    \brief Evaluate nodal distances and linearization

  */
  virtual void EvaluateDistances(const Teuchos::RCP<const Epetra_Vector>& vec,
      std::map<int, std::vector<double> >& mynormals,
      std::map<int,std::vector<GEN::pairedvector<int,double> > >& dmynormals,
      std::map<int,double>& mygap,
      std::map<int,std::map<int,double> >& dmygap);

  /*!
  \brief Assemble slave coordinates (xs)

  */
  virtual void AssembleSlaveCoord(Teuchos::RCP<Epetra_Vector>& xsmod);

  /*!
  \brief Evaluate L2 Norm of tangential contact conditions

  */
  virtual void EvaluateTangentNorm(double& cnormtan);

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
  virtual void AssembleRegTangentForcesPenalty();

  /*!
  \brief Assemble geometry dependent, tangential lagrange multipliers
  and their derivatives in the Uzawa augmented lagrange case
  */
  virtual void AssembleRegTangentForcesUzawa();

  /*
  \brief Assemble LM derivatives into global matrix (penalty strategy)

  */
  virtual void AssembleLinZ(LINALG::SparseMatrix& lambdaglobal);

  /*!
  \brief Assemble matrix T containing nodal tangents and/or matrix N containing nodal normals!

  */
  virtual void AssembleTN(Teuchos::RCP<LINALG::SparseMatrix> tglobal = Teuchos::null,
                          Teuchos::RCP<LINALG::SparseMatrix> nglobal = Teuchos::null);

  /*!
  \brief Assemble matrix S containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the normal contact condition g~ = 0. Concretely, this
  includes assembling the linearizations of the slave side
  nodal normals and of the Mortar matrices D  and M.

  */
  virtual void AssembleS(LINALG::SparseMatrix& sglobal);

  /*!
  \brief Assemble matrix Tderiv and Nderiv containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the tangential contact condition (frictionless) and/or normal
  condition (tractionlss).
  Concretely, this means assembling the linearization of the slave side
  nodal tangents / nodal normals and the current Lagrange multipliers.

    usePoroLM: linearisation will be multiplied with ...
   - true ->  poro no penetration lagrange multiplier!
   - false -> standard contact lagrange multiplier!

  */
  virtual void AssembleTNderiv(Teuchos::RCP<LINALG::SparseMatrix> tderivglobal = Teuchos::null,
      Teuchos::RCP<LINALG::SparseMatrix> nderivglobal = Teuchos::null, bool usePoroLM = false);

  /*!
  \brief Assemble matrices LinD, LinM containing linearizations

  This method builds an algebraic form of the FULL linearization
  of the contact force vector. Concretely, this means assembling
  the linearization of the Mortar matrices D and M and the
  current Lagrange multipliers.

  usePoroLM: linearisation will be multiplied with ...
   - true ->  poro no penetration lagrange multiplier!
   - false -> standard contact lagrange multiplier!

  */
  virtual void AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                             LINALG::SparseMatrix& linmglobal,
                             bool usePoroLM = false);

  /*!
  \brief subroutine assemble lin d
  */
  virtual void AssembleLinD(LINALG::SparseMatrix& lindglobal,
                             bool usePoroLM = false);

  /*!
  \brief subroutine assemble lin m
  */
  virtual void AssembleLinM(LINALG::SparseMatrix& linmglobal,
                             bool usePoroLM = false);


  /*!
  \brief Assemble weighted gap g

  Derived version! It is very important to note that g has a different
  meaning here in contact than in standard mortar meshtying applications,
  thus we need a derived method. Referring to MORTAR::Interface::AssembleG(),
  we notice that g is a vector-quantity at each node there. Yet, in
  (frictionless) we are only interested in the normal part, which makes
  g a scalar quantity here. Compare also the different definitions of g_
  in CONTACT::MtAbstractStrategy::MortarCoupling() -> gsdofrowmap_ and
  in CONTACT::CoAbstractStrategy::InitMortar()/AssembleMortar() -> gsnoderowmap_!!!

  */
  virtual void AssembleG(Epetra_Vector& gglobal);

  /*!
  \brief Assemble inactive rhs (incremental delta_z_)
  */
  virtual void AssembleInactiverhs(Epetra_Vector& inactiverhs);

  /*!
  \brief Assemble tangential rhs (incremental delta_z_)
  */
  virtual void AssembleTangrhs(Epetra_Vector& tangrhs);

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
  \brief Build active set (nodes / dofs) of this interface

  If the flag init==true, the active set is initialized (for t=0)
  according to the contact initialization defined in the input file.

  */
  virtual bool BuildActiveSet(bool init=false);

  /*!
  \brief Split active dofs into N- and T-part

  */
  virtual bool SplitActiveDofs();

  /*!
  \brief Assemble normal coupling weighted condition for poro contact

  */
  virtual void AssembleNCoup(Epetra_Vector& gglobal);

  /*!
  \brief Assemble linearisation of normal coupling weighted condition for poro contact

  */
  virtual void AssembleNCoupLin(LINALG::SparseMatrix& sglobal,
      ADAPTER::Coupling& coupfs,
      bool AssembleVelocityLin = false //if true velocity linearisation will be assembled into sglobal, otherwise lin. w.r.t. displacements!
      );

  //@}

  //! @name Visualization and Debugging methods

  /*!
  \brief Calculate the resultant angular momentum
  Necessary to check the angular momentum conservation

  */
  void EvalResultantMoment(const Epetra_Vector& fs,
                           const Epetra_Vector& fm);

  /*!
  \brief Visualize contact stuff with gmsh

  */
  void VisualizeGmsh(const int step, const int iter);

  /*!
  \brief Check normal/tangent derivatives with finite differences

  */
  void FDCheckNormalDeriv();

  /*!
  \brief Check normal/tangent derivatives with finite differences

  */
  void FDCheckNormalCPPDeriv();

  /*!
  \brief Check Mortar matrix D derivatives with finite differences

  */
  void FDCheckMortarDDeriv();

  /*!
  \brief Check Mortar matrix M derivatives with finite differences

  */
  void FDCheckMortarMDeriv();

  /*!
  \brief Check weighted gap g derivatives with finite differences

  */
  void FDCheckGapDeriv();


  /*!
  \brief Check alpha derivatives with finite differences (for hybrid formulation)

  */
  void FDCheckAlphaDeriv();


  /*!
  \brief Check weighted slip increment derivatives with finite differences (gp-wise calculated)

  */
  void FDCheckSlipIncrDerivTXI(); //- TXI
  void FDCheckSlipIncrDerivTETA(); //- TETA

  /*!
  \brief Check tangential LM derivatives with finite differences

  */
  void FDCheckTangLMDeriv();

  /*!
  \brief Check stick condition derivatives with finite differences

  */
  virtual void FDCheckStickDeriv(LINALG::SparseMatrix& linstickLMglobal,
                                 LINALG::SparseMatrix& linstickDISglobal);

  /*!
  \brief Check slip condition derivatives with finite differences

  */
  virtual void FDCheckSlipDeriv(LINALG::SparseMatrix& linslipLMglobal,
                                LINALG::SparseMatrix& linslipDISglobal);

  /*!
  \brief Check penalty approach with finite differences

  */
  void FDCheckPenaltyTracNor();

  /*!
  \brief Check frictional penalty traction with finite differences

  */
  virtual void FDCheckPenaltyTracFric();

  void AddGPTSforces(Teuchos::RCP<Epetra_FEVector> feff);
  void AddGPTSstiffness(Teuchos::RCP<LINALG::SparseMatrix> kteff);

  Teuchos::RCP<Epetra_Vector>& GetCn()
  {
    return cnValues_;
  };

  Epetra_Vector& GetCnRef()
  {
    if (cnValues_.is_null())
      dserror("The cnValues_ is not initialized!");
    return *cnValues_;
  }

  Teuchos::RCP<Epetra_Vector>& GetCt()
  {
    return ctValues_;
  };

  Epetra_Vector& GetCtRef()
  {
    if (ctValues_.is_null())
      dserror("The ctValues_ is not initialized!");
    return *ctValues_;
  }
  //@}

protected:

  /*!
  \brief initialize node and element data container

  Derived version!

  */
  virtual void InitializeDataContainer();

  /*!
  \brief initialize slave/master node status for corner/edge modification

  Derived version!

  */
  virtual void InitializeCornerEdge();

  /*!
  \brief set cpp normal and tangent to data container

  Derived version!

  */
  virtual void SetCPPNormal(
      MORTAR::MortarNode& snode,
      double* normal,
      std::vector<GEN::pairedvector<int,double> >& normallin);

  /*!
  \brief do calculations which are required for contact term evaluation:
         for example: nodal normal calculation

  Derived version!

  */
  virtual void PreEvaluate();

  /*!
  \brief Routine to control contact term evaluation. Here, we decide if mortar, nts
         etc. is evaluated

  Derived version!

  */
  virtual void EvaluateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

  /*!
  \brief export master nodal normals for cpp calculation

  */
  virtual void ExportMasterNodalNormals();

  /*!
  \brief evaluate cpp normals on slave side based on averaged normal field on master side

  */
  virtual void EvaluateCPPNormals();

  /*!
  \brief do calculations which are required after contact term evaluation:
         for example: scale nodal entries

  Derived version!

  */
  virtual void PostEvaluate(const int step, const int iter);

  /*!
  \brief Compute cpp normal based on averaged nodal normal field on master side.

  */
  virtual double ComputeCPPNormal(
      MORTAR::MortarNode& mrtrnode,
      std::vector<MORTAR::MortarElement*> meles,
      double* normal,
      std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

  /*!
  \brief 2D routine for cpp normal

  */
  virtual double ComputeCPPNormal2D(
      MORTAR::MortarNode& mrtrnode,
      std::vector<MORTAR::MortarElement*> meles,
      double* normal,
      std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

  /*!
  \brief 3D routine for cpp normal

  */
  virtual double ComputeCPPNormal3D(
      MORTAR::MortarNode& mrtrnode,
      std::vector<MORTAR::MortarElement*> meles,
      double* normal,
      std::vector<GEN::pairedvector<int,double> >& normaltolineLin);

  /*!
  \brief Compute normal between slave and master node

  */
  virtual double ComputeNormalNodeToNode(
      MORTAR::MortarNode& snode,
      MORTAR::MortarNode& mnode,
      double* normal,
      std::vector<GEN::pairedvector<int,double> >& normaltonodelin);
  /*!
  \brief Compute scaling factors for transition between nts, mortar etc.

  */
  virtual void ComputeScaling();

  /*!
  \brief 2D version of scaling computation

  */
  virtual void ComputeScaling2D();

  /*!
  \brief 3D version of scaling computation

  */
  virtual void ComputeScaling3D();

  /*!
  \brief Set new cn and ct values (global interface vector)

  */
  virtual void SetCnCtValues();

  /*!
  \brief Routine scale between nts, mortar, lts etc. This is required for non-smooth
         contact geometries

  */
  virtual void ScaleTerms();

  /*!
  \brief 2D version of scaling routine

  */
  virtual void ScaleTerms2D();

  /*!
  \brief 3D version of scaling routine

  */
  virtual void ScaleTerms3D();

  /*!
  \brief Routine to scale nodal normals between nodal averaged normal and cpp normal

  */
  virtual void ScaleNormals();

  /*!
  \brief 2D version of nodal normal scaling

  */
  virtual void ScaleNormals2D();

  /*!
  \brief 3D version of nodal normal scaling

  */
  virtual void ScaleNormals3D();

  /*!
  \brief Routine which stores entries from nts algorithm into mortar nodes to reuse
         standard assemble functions

  */
  virtual void StoreNTSvalues();

  /*!
  \brief Routine which stores entries from lts algorithm into mortar nodes to reuse
         standard assemble functions

  */
  virtual void StoreLTSvalues();

  /*!
  \brief Routine which stores entries from ltl algorithm into mortar nodes to reuse
         standard assemble functions

  */
  virtual void StoreLTLvalues();

  /*!
  \brief These functions are not properly implemented/used!!!!

  */
  virtual void DetectNonSmoothGeometries();
  virtual void EvaluateAveragedNodalNormals();

protected:
  // don't want = operator and cctor
  CoInterface operator = (const CoInterface& old);
  CoInterface(const CoInterface& old);

  bool                          selfcontact_;           // flag indicating if this is a self contact interface
  bool                          friction_;              // flag for frictional contact
  bool                          nonSmoothContact_;      // flag for non-smooth contact algorithm
  INPAR::CONTACT::ConstraintDirection constr_direction_;// direction in which the contact constraints are formulated

  Teuchos::RCP<Epetra_Map>      activenodes_;    // row map of all active slave nodes
  Teuchos::RCP<Epetra_Map>      activedofs_;     // row map of all active slave dofs
  Teuchos::RCP<Epetra_Map>      activen_;        // row map of global N-matrix
  Teuchos::RCP<Epetra_Map>      activet_;        // row map of global T-matrix
  Teuchos::RCP<Epetra_Map>      slipnodes_;      // row map of all slip slave nodes
  Teuchos::RCP<Epetra_Map>      slipdofs_;       // row map of all slip slave dofs
  Teuchos::RCP<Epetra_Map>      slipt_;          // row map of part of T-matrix (slip nodes)

  Teuchos::RCP<Epetra_Map>      nonsmoothnodes_; // row map of all nonsmooth slave nodes
  Teuchos::RCP<Epetra_Map>      smoothnodes_;    // row map of all smooth slave nodes

  Teuchos::RCP<Epetra_Map>      nextendedghosting_;
  Teuchos::RCP<Epetra_Map>      eextendedghosting_;

  Teuchos::RCP<SelfBinaryTree>  binarytreeself_; // binary tree for self contact search

  //! cn-values of each node
  Teuchos::RCP<Epetra_Vector>   cnValues_;         // cn
  Teuchos::RCP<Epetra_Vector>   ctValues_;         // ct

  int                           smpairs_;        // proc local number of slave/master pairs
  int                           smintpairs_;     // proc local number of slave/master integration pairs
  int                           intcells_;       // proc local number of integration cells

private:
  static bool abs_compare(int a, int b)
  {
      return (std::abs(a) < std::abs(b));
  }
}; // class CoInterface
} // namespace CONTACT


// << operator
std::ostream& operator << (std::ostream& os, const CONTACT::CoInterface& interface);


#endif  // #ifndef CONTACT_INTERFACE_H
