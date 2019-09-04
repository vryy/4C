/*---------------------------------------------------------------------------*/
/*! \file
\brief Interface class for the eXtended contact evaluation

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------------*/

#include "../drt_contact_xcontact/xcontact_interface.H"
#include "../drt_structure_xstructure/xstr_multi_discretization_wrapper.H"

#include "../drt_contact/contact_node.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
XCONTACT::Interface::Interface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr, const int id,
    const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
    bool selfcontact, INPAR::MORTAR::RedundantStorage redundant)
    : CONTACT::CoInterface(interfaceData_ptr, id, comm, dim, icontact, selfcontact, redundant),
      sndofrowmap_(Teuchos::null),
      stdofrowmap_(Teuchos::null),
      parent_discret_(*icontact.get<Teuchos::RCP<XSTR::MultiDiscretizationWrapper::cXDisPair>>(
          "ParentDiscretPair"))
{
  // Empty constructor body
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::Initialize()
{
  // Return if not participating in interface
  if (!lComm())
  {
    return;
  }

  // ========================================================================
  // Initialize all general contact related quantities
  // ========================================================================

  // Call initialization routine of contact interface
  CoInterface::Initialize();


  // ========================================================================
  // Initialize all XContact related quantities
  // ========================================================================

  // Loop over all slave nodes to reset quantities (standard column map)
  // (include slave side boundary nodes/crosspoints)
  for (int i = 0; i < SlaveColNodesBound()->NumMyElements(); ++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node)
    {
      dserror("Cannot find node with gid %d.", gid);
    }

    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    cnode->CoData().GetWcLm() = 0.0;
    cnode->CoData().GetWcSuLm().clear();
    cnode->CoData().GetWcMuLm().clear();
    cnode->CoData().GetWcSuU().clear();
    cnode->CoData().GetWcMuU().clear();
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::PreMortarCoupling(const MORTAR::MortarElement* sele,
    const std::vector<MORTAR::MortarElement*> mele,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr) const
{
  // ToDo check if it is possible to delete the nodal data containers
  // create the element vectors and element matrices at this point
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::PostMortarCoupling(const MORTAR::MortarElement* sele,
    const std::vector<MORTAR::MortarElement*> mele,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr) const
{
  // ToDo check if it is possible to delete the nodal data containers
  // assemble the element vectors and element matrices at this point
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::AssembleWeightedGap(Epetra_Vector& wgap) const
{
  // Return if not participating in interface
  if (!lComm()) return;

  // Loop over procs active slave nodes of the interface for assembly
  // (use standard row map to assemble each node only once)
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
    {
      dserror("Cannot find slave node with gid %d", gid);
    }
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    std::vector<int> rgid(1, sndofrowmap_->GID(i));
    std::vector<int> rowner(1, cnode->Owner());

    // Get Lagrange multiplier in normal direction of
    Epetra_SerialDenseVector wgap_node(1);
    wgap_node[0] = cnode->CoData().Getg();

    // Assemble nodal quantities to global matrices
    LINALG::Assemble(wgap, wgap_node, rgid, rowner);
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::AssembleContactRHS(Epetra_Vector& Wc_lm, Epetra_Vector& lmN) const
{
  // Return if not participating in interface
  if (!lComm())
  {
    return;
  }

  // Loop over procs active slave nodes of the interface for assembly
  // (use standard row map to assemble each node only once)
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
    {
      dserror("Cannot find slave node with gid %d", gid);
    }
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    std::vector<int> rgid(1, cnode->Dofs()[0]);
    std::vector<int> rowner(1, cnode->Owner());

    // Get Lagrange multiplier in normal direction of
    Epetra_SerialDenseVector lmN_node(1);
    lmN_node[0] = cnode->MoData().lm()[0];

    // Get contact residuum of node
    Epetra_SerialDenseVector Wc_lm_node(1);
    Wc_lm_node[0] = cnode->CoData().GetWcLm();

    // Assemble nodal quantities to global matrices
    LINALG::Assemble(Wc_lm, Wc_lm_node, rgid, rowner);
    LINALG::Assemble(lmN, lmN_node, rgid, rowner);
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::AssembleMortar(LINALG::SparseMatrix& D, LINALG::SparseMatrix& M) const
{
  // Return if not participating in interface
  if (!lComm())
  {
    return;
  }

  // Loop over procs active slave nodes of the interface for assembly
  // (use standard row map to assemble each node only once)
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    const int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
    {
      dserror("Cannot find node with gid %d.", gid);
    }

    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (cnode->Owner() != Comm().MyPID())
    {
      dserror("Node ownership inconsistency.");
    }

    // get the corresponding slave normal dof gid
    const int rowId = sndofrowmap_->GID(i);

    // Define type of map iterators
    typedef GEN::pairedvector<int, double>::const_iterator CI;
    typedef std::map<int, double>::const_iterator CIP;


    // ========================================================================
    // Assemble mortar matrix D (slave)
    // ========================================================================

    GEN::pairedvector<int, double>& Wc_su_lm = cnode->CoData().GetWcSuLm();

    for (CI p = Wc_su_lm.begin(); p != Wc_su_lm.end(); ++p)
    {
      const int cDofId = p->first;
      const double cval = p->second;
      D.Assemble(cval, rowId, cDofId);
    }


    // ========================================================================
    // Assemble mortar matrix M (master)
    // ========================================================================

    std::map<int, double>& Wc_mu_lm = cnode->CoData().GetWcMuLm();

    for (CIP p = Wc_mu_lm.begin(); p != Wc_mu_lm.end(); ++p)
    {
      const int cDofId = p->first;
      const double cval = p->second;
      M.Assemble(cval, rowId, cDofId);
    }
  }
}

/*---------------------------------------------------------------------------*
 | Assemble structural contact tangent matrix                    Hofer 08/16 |
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::AssembleWcUU(
    LINALG::SparseMatrix& Wc_su_u, LINALG::SparseMatrix& Wc_mu_u) const
{
  // Return if not participating in interface
  if (!lComm())
  {
    return;
  }

  // Loop over procs active slave nodes of the interface for assembly
  // (use standard row map to assemble each node only once)
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(idiscret_->gNode(gid));
    if (!cnode)
    {
      dserror("Cannot find node with gid %d.", gid);
    }

    // Get Lagrange multiplier in normal direction
    double lmN = cnode->MoData().lm()[0];

    // Define type of map iterators
    typedef std::map<int, double>::const_iterator CI;
    typedef std::map<int, std::map<int, double>>::const_iterator CII;


    // ========================================================================
    // Slave varied part
    // ========================================================================

    std::map<int, std::map<int, double>>& Wc_su_u_node = cnode->CoData().GetWcSuU();

    // iteration over ALL slave Dof Ids
    for (CII p = Wc_su_u_node.begin(); p != Wc_su_u_node.end(); ++p)
    {
      int sRow = p->first;

      for (CI pp = Wc_su_u_node[sRow].begin(); pp != Wc_su_u_node[sRow].end(); ++pp)
      {
        int col = pp->first;
        double val = pp->second * lmN;
        if (abs(val) > 1.0e-12)
        {
          Wc_su_u.FEAssemble(val, sRow, col);
        }
      }
    }

    // ========================================================================
    // Master varied part
    // ========================================================================

    std::map<int, std::map<int, double>>& Wc_mu_u_node = cnode->CoData().GetWcMuU();

    // iteration over ALL master Dof Ids
    for (CII p = Wc_mu_u_node.begin(); p != Wc_mu_u_node.end(); ++p)
    {
      const int mRow = p->first;

      for (CI pp = Wc_mu_u_node[mRow].begin(); pp != Wc_mu_u_node[mRow].end(); ++pp)
      {
        const int col = pp->first;
        const double val = pp->second * lmN;
        if (abs(val) > 1.0e-12)
        {
          Wc_mu_u.FEAssemble(val, mRow, col);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::UpdateMasterSlaveSets()
{
  // TODO: Copied from CONTACT::AugmentedInterface::UpdateMasterSlaveSets

  //
  MORTAR::MortarInterface::UpdateMasterSlaveSets();

  //
  SplitSlaveDofs();
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XCONTACT::Interface::SplitSlaveDofs()
{
  // Return if there is no contact
  if (snoderowmap_ == Teuchos::null || snoderowmap_->NumGlobalElements() == 0)
  {
    sndofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    stdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    return;
  }

  // Define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNGids(snoderowmap_->NumMyElements());
  std::vector<int> myTGids((Dim() - 1) * snoderowmap_->NumMyElements());

  // Check dimension
  double dimcheck = sdofrowmap_->NumGlobalElements() / snoderowmap_->NumGlobalElements();
  if (dimcheck != Dim())
  {
    dserror("Split slave DOFs: Nodes <-> DOFs dimension mismatch.");
  }

  // Loop over all contact nodes
  for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
  {
    int gid = snoderowmap_->GID(k);
    CONTACT::CoNode* cnode = static_cast<CONTACT::CoNode*>(idiscret_->gNode(gid));
    if (!cnode)
    {
      dserror("Cannot find slave node with gid %d", gid);
    }

    // Add first DOF to normal map
    myNGids[countN] = cnode->Dofs()[0];
    ++countN;

    // Add remaining DOFs to tangential map
    for (int i = 1; i < cnode->NumDof(); ++i)
    {
      myTGids[countT] = cnode->Dofs()[i];
      ++countT;
    }
  }

  // Resize the temporary vectors
  myNGids.resize(countN);
  myTGids.resize(countT);

  // Communicate countN and countT among processors
  int gCountN;
  int gCountT;
  Comm().SumAll(&countN, &gCountN, 1);
  Comm().SumAll(&countT, &gCountT, 1);

  // Check global dimensions
  if (gCountN + gCountT != sdofrowmap_->NumGlobalElements())
    dserror("Split slave DOFs: Splitting went wrong.");

  // Create normal map and tangential map objects
  sndofrowmap_ = Teuchos::rcp(new Epetra_Map(gCountN, countN, &myNGids[0], 0, Comm()));
  stdofrowmap_ = Teuchos::rcp(new Epetra_Map(gCountT, countT, &myTGids[0], 0, Comm()));
}
