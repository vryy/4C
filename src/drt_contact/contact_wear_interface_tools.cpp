/*-----------------------------------------------------------------------*/
/*! \file
\level 2

\maintainer Matthias Mayr

\brief Some tools for wear problems
*/
/*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    farah 09/13|
 *----------------------------------------------------------------------*/
#include "contact_wear_interface.H"
#include "contact_integrator.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "selfcontact_binarytree.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_dofset.H"
#include "../drt_mortar/mortar_integrator.H"
#include "../drt_mortar/mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_io/io_control.H"


/*----------------------------------------------------------------------*
 | Finite difference check for normal gap derivatives        farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckGapDeriv()
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for gap values
  int nrow = snoderowmap_->NumMyElements();
  std::vector<double> refG(nrow);
  std::vector<double> newG(nrow);

  // problem dimension (2D or 3D)
  int dim = Dim();

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    // store gap-values into refG
    refG[i] = cnode->CoData().Getg();
  }

  // global loop to apply FD scheme to all slave dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];
    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // apply finite difference scheme
    /*if (Comm().MyPID()==snode->Owner())
    {
      std::cout << "\nBuilding FD for Slave Node: " << snode->Id() << " Dof(l): " << fd%dim
           << " Dof(g): " << snode->Dofs()[fd%dim] << std::endl;
    }*/

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);

      // store gap-values into newG
      newG[k] = kcnode->CoData().Getg();

      if (abs(newG[k] - refG[k]) > 1e-12 && newG[k] != 1.0e12 && refG[k] != 1.0e12)
      {
        double finit = (newG[k] - refG[k]) / delta;
        double analy = kcnode->CoData().GetDerivG()[snode->Dofs()[fd % dim]];
        double dev = finit - analy;

        // kgid: id of currently tested slave node
        // snode->Dofs()[fd%dim]: currently modified slave dof
        std::cout << "(" << kgid << "," << snode->Dofs()[fd % dim] << ") : fd=" << finit
                  << " derivg=" << analy << " DEVIATION " << dev;

        if (abs(dev) > 1e-4)
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if (abs(dev) > 1e-5)
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // global loop to apply FD scheme to all master dofs (=dim*nodes)
  for (int fd = 0; fd < dim * mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find master node with gid %", gid);
    CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(node);

    int mdof = mnode->Dofs()[fd % dim];
    std::cout << "\nDERIVATIVE FOR M-NODE # " << gid << " DOF: " << mdof << std::endl;

    // apply finite difference scheme
    /*if (Comm().MyPID()==mnode->Owner())
    {
      std::cout << "\nBuilding FD for Master Node: " << mnode->Id() << " Dof(l): " << fd%dim
           << " Dof(g): " << mnode->Dofs()[fd%dim] << std::endl;
    }*/

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      mnode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      mnode->xspatial()[1] += delta;
    }
    else
    {
      mnode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);

      if (kcnode->Active())
      {
        // check two versions of weighted gap
        double defgap = 0.0;
        double wii = (kcnode->MoData().GetD())[kcnode->Id()];

        for (int j = 0; j < dim; ++j)
          defgap -= (kcnode->MoData().n()[j]) * wii * (kcnode->xspatial()[j]);

        std::map<int, double>& mmap = kcnode->MoData().GetM();
        std::map<int, double>::const_iterator mcurr;

        for (int m = 0; m < mnodefullmap->NumMyElements(); ++m)
        {
          int gid = mnodefullmap->GID(m);
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %", gid);
          CONTACT::CoNode* cmnode = dynamic_cast<CONTACT::CoNode*>(mnode);
          bool hasentry = false;

          // look for this master node in M-map of the active slave node
          for (mcurr = mmap.begin(); mcurr != mmap.end(); ++mcurr)
            if ((mcurr->first) == cmnode->Id())
            {
              hasentry = true;
              break;
            }

          double mik = mmap[cmnode->Id()];
          double* mxi = cmnode->xspatial();

          // get out of here, if master node not adjacent or coupling very weak
          if (!hasentry || abs(mik) < 1.0e-12) continue;

          for (int j = 0; j < dim; ++j) defgap += (kcnode->MoData().n()[j]) * mik * mxi[j];
        }

        // std::cout << "SNode: " << kcnode->Id() << " IntGap: " << kcnode->CoData().Getg << "
        // DefGap: " << defgap << std::endl; kcnode->CoData().Getg = defgap;
      }

      // store gap-values into newG
      newG[k] = kcnode->CoData().Getg();

      if (abs(newG[k] - refG[k]) > 1e-12 && newG[k] != 1.0e12 && refG[k] != 1.0e12)
      {
        double finit = (newG[k] - refG[k]) / delta;
        double analy = kcnode->CoData().GetDerivG()[mnode->Dofs()[fd % dim]];
        double dev = finit - analy;

        // kgid: id of currently tested slave node
        // mnode->Dofs()[fd%dim]: currently modified slave dof
        std::cout << "(" << kgid << "," << mnode->Dofs()[fd % dim] << ") : fd=" << finit
                  << " derivg=" << analy << " DEVIATION " << dev;

        if (abs(dev) > 1e-4)
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if (abs(dev) > 1e-5)
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      mnode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      mnode->xspatial()[1] -= delta;
    }
    else
    {
      mnode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // back to normal...
  Initialize();
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for normal gap derivatives        farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckGapDeriv_W()
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for gap values
  int nrow = snoderowmap_->NumMyElements();
  std::vector<double> refG(nrow);
  std::vector<double> newG(nrow);

  // problem dimension (2D or 3D)
  int dim = Dim();

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    // store gap-values into refG
    refG[i] = cnode->CoData().Getg();
  }

  // global loop to apply FD scheme to all slave dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    int sdof = snode->Dofs()[fd % dim];
    std::cout << "\nW --- DERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    double delta = 1e-8;
    if (fd % dim == 0) snode->WearData().wcurr()[0] += delta;


    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::CoNode* kcnode = dynamic_cast<CONTACT::CoNode*>(knode);

      // store gap-values into newG
      newG[k] = kcnode->CoData().Getg();

      if (abs(newG[k] - refG[k]) > 1e-12 && newG[k] != 1.0e12 && refG[k] != 1.0e12)
      {
        double finit = (newG[k] - refG[k]) / delta;
        double analy = kcnode->CoData().GetDerivGW()[snode->Dofs()[0]];
        double dev = finit - analy;

        // kgid: id of currently tested slave node
        // snode->Dofs()[fd%dim]: currently modified slave dof
        std::cout << "(" << kgid << "," << snode->Dofs()[fd % dim] << ") : fd=" << finit
                  << " derivg=" << analy << " DEVIATION " << dev;

        if (abs(dev) > 1e-4)
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if (abs(dev) > 1e-5)
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0) snode->WearData().wcurr()[0] -= delta;


    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // back to normal...
  Initialize();
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check of wear condition derivatives     farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckDerivE_D(LINALG::SparseMatrix& linedis)
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for values of complementary function C
  int nrow = snoderowmap_->NumMyElements();
  std::vector<double> reft(nrow);
  std::vector<double> newt(nrow);

  int dim = Dim();
  typedef std::map<int, double>::const_iterator CI;

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    reft[i] = 0.0;

    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->WearData().GetE().size() > 0)
    {
      std::map<int, double> scurr = cnode->WearData().GetE()[0];

      for (CI p = scurr.begin(); p != scurr.end(); ++p)
      {
        int gid2 = (int)((p->first) / (dim));
        DRT::Node* node2 = idiscret_->gNode(gid2);
        if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
        CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

        double w = cnode2->WearData().wcurr()[0];

        reft[i] += (p->second) * w;
      }
    }
    else
    {
      reft[i] = 0.0;
    }
  }  // loop over procs slave nodes

  // ********************************************************************************************
  // global loop to apply FD scheme to all slave dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->xspatial()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      newt[k] = 0.0;

      int gid3 = snoderowmap_->GID(k);
      DRT::Node* node3 = idiscret_->gNode(gid3);
      if (!node3) dserror("ERROR: Cannot find node with gid %", gid3);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(node3);

      if (kcnode->WearData().GetE().size() > 0)
      {
        std::map<int, double> sscurr = kcnode->WearData().GetE()[0];

        for (CI p = sscurr.begin(); p != sscurr.end(); ++p)
        {
          int gid2 = (int)((p->first) / (dim));
          DRT::Node* node2 = idiscret_->gNode(gid2);
          if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
          CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

          double w = cnode2->WearData().wcurr()[0];

          newt[k] += (p->second) * w;
        }
      }
      else
        newt[k] = 0.0;

      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linedis.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[0]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[0], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;


      // std::cout << "ref= " << reft[k] << "\t new= " << newt[k] << std::endl;

      if (abs(newt[k] - reft[k]) > 1e-12)
      {
        std::cout << "E WEAR DIS-Deriv: " << kcnode->Dofs()[0]
                  << "\t w.r.t Slave: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newt[k] - reft[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newt[k] - reft[k]) / delta);
        if (abs(analyt_txi - (newt[k] - reft[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }
  }  // loop over procs slave nodes

  // ********************************************************************************************
  // global loop to apply FD scheme to all master dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * mnodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->xspatial()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      newt[k] = 0.0;

      int gid3 = snoderowmap_->GID(k);
      DRT::Node* node3 = idiscret_->gNode(gid3);
      if (!node3) dserror("ERROR: Cannot find node with gid %", gid3);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(node3);

      if (kcnode->WearData().GetE().size() > 0)
      {
        std::map<int, double> sscurr = kcnode->WearData().GetE()[0];

        for (CI p = sscurr.begin(); p != sscurr.end(); ++p)
        {
          int gid2 = (int)((p->first) / (dim));
          DRT::Node* node2 = idiscret_->gNode(gid2);
          if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
          CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

          double w = cnode2->WearData().wcurr()[0];

          newt[k] += (p->second) * w;
        }
      }
      else
        newt[k] = 0.0;

      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linedis.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[0]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[0], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;


      // std::cout << "ref= " << reft[k] << "\t new= " << newt[k] << std::endl;

      if (abs(newt[k] - reft[k]) > 1e-12)
      {
        std::cout << "E WEAR DIS-Deriv: " << kcnode->Dofs()[0]
                  << "\t w.r.t Master: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newt[k] - reft[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newt[k] - reft[k]) / delta);
        if (abs(analyt_txi - (newt[k] - reft[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }
  }  // loop over procs slave nodes
  // back to normal...
  Initialize();
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check of wear condition derivatives     farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckDerivE_D_Master(LINALG::SparseMatrix& linedis)
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for values of complementary function C
  int nrow = snoderowmap_->NumMyElements();
  std::vector<double> reft(nrow);
  std::vector<double> newt(nrow);

  int dim = Dim();
  typedef std::map<int, double>::const_iterator CI;

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
  {
    reft[i] = 0.0;

    int gid = mnoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->WearData().GetE().size() > 0)
    {
      std::map<int, double> scurr = cnode->WearData().GetE()[0];

      for (CI p = scurr.begin(); p != scurr.end(); ++p)
      {
        int gid2 = (int)((p->first) / (dim));
        DRT::Node* node2 = idiscret_->gNode(gid2);
        if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
        CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

        double w = cnode2->WearData().wcurr()[0];

        reft[i] += (p->second) * w;
      }
    }
    else
    {
      reft[i] = 0.0;
    }
  }  // loop over procs slave nodes

  // ********************************************************************************************
  // global loop to apply FD scheme to all slave dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->xspatial()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < mnoderowmap_->NumMyElements(); ++k)
    {
      newt[k] = 0.0;

      int gid3 = mnoderowmap_->GID(k);
      DRT::Node* node3 = idiscret_->gNode(gid3);
      if (!node3) dserror("ERROR: Cannot find node with gid %", gid3);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(node3);

      if (kcnode->WearData().GetE().size() > 0)
      {
        std::map<int, double> sscurr = kcnode->WearData().GetE()[0];

        for (CI p = sscurr.begin(); p != sscurr.end(); ++p)
        {
          int gid2 = (int)((p->first) / (dim));
          DRT::Node* node2 = idiscret_->gNode(gid2);
          if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
          CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

          double w = cnode2->WearData().wcurr()[0];

          newt[k] += (p->second) * w;
        }
      }
      else
        newt[k] = 0.0;

      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linedis.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[0]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[0], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;


      // std::cout << "ref= " << reft[k] << "\t new= " << newt[k] << std::endl;

      if (abs(newt[k] - reft[k]) > 1e-12)
      {
        std::cout << "E WEAR DIS-Deriv: " << kcnode->Dofs()[0]
                  << "\t w.r.t Slave: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newt[k] - reft[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newt[k] - reft[k]) / delta);
        if (abs(analyt_txi - (newt[k] - reft[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }
  }  // loop over procs slave nodes

  // ********************************************************************************************
  // global loop to apply FD scheme to all master dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * mnodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->xspatial()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < mnoderowmap_->NumMyElements(); ++k)
    {
      newt[k] = 0.0;

      int gid3 = mnoderowmap_->GID(k);
      DRT::Node* node3 = idiscret_->gNode(gid3);
      if (!node3) dserror("ERROR: Cannot find node with gid %", gid3);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(node3);

      if (kcnode->WearData().GetE().size() > 0)
      {
        std::map<int, double> sscurr = kcnode->WearData().GetE()[0];

        for (CI p = sscurr.begin(); p != sscurr.end(); ++p)
        {
          int gid2 = (int)((p->first) / (dim));
          DRT::Node* node2 = idiscret_->gNode(gid2);
          if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
          CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

          double w = cnode2->WearData().wcurr()[0];

          newt[k] += (p->second) * w;
        }
      }
      else
        newt[k] = 0.0;

      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linedis.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[0]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[0], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;


      // std::cout << "ref= " << reft[k] << "\t new= " << newt[k] << std::endl;

      if (abs(newt[k] - reft[k]) > 1e-12)
      {
        std::cout << "E WEAR DIS-Deriv: " << kcnode->Dofs()[0]
                  << "\t w.r.t Master: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newt[k] - reft[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newt[k] - reft[k]) / delta);
        if (abs(analyt_txi - (newt[k] - reft[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }
  }  // loop over procs slave nodes
  // back to normal...
  Initialize();
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check of wear condition derivatives     farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckDerivT_D(LINALG::SparseMatrix& lintdis)
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements() == 0) return;

  // create storage for values of complementary function C
  int nrow = slipnodes_->NumMyElements();
  std::vector<double> reft(nrow);
  std::vector<double> newt(nrow);

  int dim = Dim();
  typedef std::map<int, double>::const_iterator CI;

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < slipnodes_->NumMyElements(); ++i)
  {
    reft[i] = 0.0;

    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->WearData().GetT().size() > 0)
    {
      std::map<int, double> scurr = cnode->WearData().GetT()[0];

      for (CI p = scurr.begin(); p != scurr.end(); ++p)
      {
        int gid2 = (int)((p->first) / (dim));
        DRT::Node* node2 = idiscret_->gNode(gid2);
        if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
        CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

        double lmn = 0.0;
        for (int g = 0; g < dim; ++g) lmn += (cnode2->MoData().n()[g]) * (cnode2->MoData().lm()[g]);

        reft[i] += (p->second) * lmn;
      }
    }
    else
    {
      reft[i] = 0.0;
    }
  }  // loop over procs slave nodes

  // ********************************************************************************************
  // global loop to apply FD scheme to all slave dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->xspatial()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < slipnodes_->NumMyElements(); ++k)
    {
      newt[k] = 0.0;

      int gid3 = slipnodes_->GID(k);
      DRT::Node* node3 = idiscret_->gNode(gid3);
      if (!node3) dserror("ERROR: Cannot find node with gid %", gid3);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(node3);

      if (kcnode->WearData().GetT().size() > 0)
      {
        std::map<int, double> sscurr = kcnode->WearData().GetT()[0];

        for (CI p = sscurr.begin(); p != sscurr.end(); ++p)
        {
          int gid2 = (int)((p->first) / (dim));
          DRT::Node* node2 = idiscret_->gNode(gid2);
          if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
          CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

          double lmn = 0.0;
          for (int g = 0; g < dim; ++g)
            lmn += (cnode2->MoData().n()[g]) * (cnode2->MoData().lm()[g]);

          newt[k] += (p->second) * lmn;
        }
      }
      else
        newt[k] = 0.0;

      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = lintdis.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[0]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[0], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;


      if (abs(newt[k] - reft[k]) > 1e-12)
      {
        std::cout << "WEAR DIS-Deriv: " << kcnode->Dofs()[0]
                  << "\t w.r.t Slave: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newt[k] - reft[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newt[k] - reft[k]) / delta);
        if (abs(analyt_txi - (newt[k] - reft[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }
  }  // loop over procs slave nodes

  // back to normal...
  Initialize();
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check of wear condition derivatives     farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckDerivT_D_Master(LINALG::SparseMatrix& lintdis)
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // nothing to do if no slip nodes
  if (slipmasternodes_->NumMyElements() == 0) return;

  // create storage for values of complementary function C
  int nrow = slipmasternodes_->NumMyElements();
  std::vector<double> reft(nrow);
  std::vector<double> newt(nrow);

  int dim = Dim();
  typedef std::map<int, double>::const_iterator CI;

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < slipmasternodes_->NumMyElements(); ++i)
  {
    reft[i] = 0.0;

    int gid = slipmasternodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->WearData().GetT().size() > 0)
    {
      std::map<int, double> scurr = cnode->WearData().GetT()[0];

      for (CI p = scurr.begin(); p != scurr.end(); ++p)
      {
        int gid2 = (int)((p->first) / (dim));
        DRT::Node* node2 = idiscret_->gNode(gid2);
        if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
        CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

        double lmn = 0.0;
        for (int g = 0; g < dim; ++g) lmn += (cnode2->MoData().n()[g]) * (cnode2->MoData().lm()[g]);

        reft[i] += (p->second) * lmn;
      }
    }
    else
    {
      reft[i] = 0.0;
    }
  }  // loop over procs slave nodes

  // ********************************************************************************************
  // global loop to apply FD scheme to all slave dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->xspatial()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < slipmasternodes_->NumMyElements(); ++k)
    {
      newt[k] = 0.0;

      int gid3 = slipmasternodes_->GID(k);
      DRT::Node* node3 = idiscret_->gNode(gid3);
      if (!node3) dserror("ERROR: Cannot find node with gid %", gid3);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(node3);

      if (kcnode->WearData().GetT().size() > 0)
      {
        std::map<int, double> sscurr = kcnode->WearData().GetT()[0];

        for (CI p = sscurr.begin(); p != sscurr.end(); ++p)
        {
          int gid2 = (int)((p->first) / (dim));
          DRT::Node* node2 = idiscret_->gNode(gid2);
          if (!node2) dserror("ERROR: Cannot find node with gid %", gid2);
          CONTACT::FriNode* cnode2 = dynamic_cast<CONTACT::FriNode*>(node2);

          double lmn = 0.0;
          for (int g = 0; g < dim; ++g)
            lmn += (cnode2->MoData().n()[g]) * (cnode2->MoData().lm()[g]);

          newt[k] += (p->second) * lmn;
        }
      }
      else
        newt[k] = 0.0;

      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = lintdis.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[0]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[0], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;


      if (abs(newt[k] - reft[k]) > 1e-12)
      {
        std::cout << "WEAR DIS-Deriv: " << kcnode->Dofs()[0]
                  << "\t w.r.t Slave: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newt[k] - reft[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newt[k] - reft[k]) / delta);
        if (abs(analyt_txi - (newt[k] - reft[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }
  }  // loop over procs slave nodes

  // back to normal...
  Initialize();
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check of slip condition derivatives     farah 08/13|
 | Not for Wear Lin. or modifications concerning the compl.             |
 | fnc. !!! See flags CONSISTENTSTICK / CONSISTENTSLIP                  |
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckSlipDeriv(LINALG::SparseMatrix& linslipLMglobal,
    LINALG::SparseMatrix& linslipDISglobal, LINALG::SparseMatrix& linslipWglobal)
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");
  double frbound = InterfaceParams().get<double>("FRBOUND");
  double frcoeff = InterfaceParams().get<double>("FRCOEFF");
  double ct = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  double cn = InterfaceParams().get<double>("SEMI_SMOOTH_CN");

  // create storage for values of complementary function C
  int nrow = snoderowmap_->NumMyElements();
  std::vector<double> refCtxi(nrow);
  std::vector<double> refCteta(nrow);
  std::vector<double> newCtxi(nrow);
  std::vector<double> newCteta(nrow);

  int dim = Dim();

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    double jumptxi = 0;
    double jumpteta = 0;
    double ztxi = 0;
    double zteta = 0;
    double znor = 0;
    double euclidean = 0;

    if (cnode->FriData().Slip())
    {
      // calculate value of C-function
      double D = cnode->MoData().GetD()[cnode->Id()];
      double Dold = cnode->FriData().GetDOld()[cnode->Id()];

      for (int dim = 0; dim < cnode->NumDof(); ++dim)
      {
        jumptxi -= (cnode->CoData().txi()[dim]) * (D - Dold) * (cnode->xspatial()[dim]);
        jumpteta -= (cnode->CoData().teta()[dim]) * (D - Dold) * (cnode->xspatial()[dim]);
        ztxi += (cnode->CoData().txi()[dim]) * (cnode->MoData().lm()[dim]);
        zteta += (cnode->CoData().teta()[dim]) * (cnode->MoData().lm()[dim]);
        znor += (cnode->MoData().n()[dim]) * (cnode->MoData().lm()[dim]);
      }

      std::map<int, double>& mmap = cnode->MoData().GetM();
      std::map<int, double>& mmapold = cnode->FriData().GetMOld();

      std::map<int, double>::const_iterator colcurr;
      std::set<int> mnodes;

      for (colcurr = mmap.begin(); colcurr != mmap.end(); colcurr++) mnodes.insert(colcurr->first);

      for (colcurr = mmapold.begin(); colcurr != mmapold.end(); colcurr++)
        mnodes.insert(colcurr->first);

      std::set<int>::iterator mcurr;

      // loop over all master nodes (find adjacent ones to this slip node)
      for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %", gid);
        CONTACT::FriNode* cmnode = dynamic_cast<CONTACT::FriNode*>(mnode);

        double mik = mmap[cmnode->Id()];
        double mikold = mmapold[cmnode->Id()];

        std::map<int, double>::iterator mcurr;

        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          jumptxi += (cnode->CoData().txi()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
          jumpteta += (cnode->CoData().teta()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
        }
      }  //  loop over master nodes

      // gp-wise slip !!!!!!!
      if (DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == true)
      {
        jumptxi = cnode->FriData().jump_var()[0];
        jumpteta = 0.0;

        if (Dim() == 3) jumpteta = cnode->FriData().jump_var()[1];
      }

      // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
      std::vector<double> sum1(Dim() - 1, 0);
      sum1[0] = ztxi + ct * jumptxi;
      if (Dim() == 3) sum1[1] = zteta + ct * jumpteta;
      if (Dim() == 2) euclidean = abs(sum1[0]);
      if (Dim() == 3) euclidean = sqrt(sum1[0] * sum1[0] + sum1[1] * sum1[1]);
    }  // if cnode == Slip

    // store C in vector
    if (ftype == INPAR::CONTACT::friction_tresca)
    {
      refCtxi[i] = euclidean * ztxi - frbound * (ztxi + ct * jumptxi);
      refCteta[i] = euclidean * zteta - frbound * (zteta + ct * jumpteta);
    }
    else if (ftype == INPAR::CONTACT::friction_coulomb)
    {
      refCtxi[i] = euclidean * ztxi - (frcoeff * znor) * (ztxi + ct * jumptxi);
      refCteta[i] = euclidean * zteta - (frcoeff * znor) * (zteta + ct * jumpteta);
    }
    else
      dserror("ERROR: Friction law is neiter Tresca nor Coulomb");

    refCtxi[i] =
        euclidean * ztxi - (frcoeff * (znor - cn * cnode->CoData().Getg())) * (ztxi + ct * jumptxi);
    refCteta[i] = euclidean * zteta -
                  (frcoeff * (znor - cn * cnode->CoData().Getg())) * (zteta + ct * jumpteta);

  }  // loop over procs slave nodes

  // **********************************************************************************
  // global loop to apply FD scheme for LM to all slave dofs (=3*nodes)
  // **********************************************************************************
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->MoData().lm()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->MoData().lm()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->MoData().lm()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!node) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      double jumptxi = 0;
      double jumpteta = 0;
      double ztxi = 0;
      double zteta = 0;
      double znor = 0;
      double euclidean = 0;

      if (kcnode->FriData().Slip())
      {
        // check two versions of weighted gap
        double D = kcnode->MoData().GetD()[kcnode->Id()];
        double Dold = kcnode->FriData().GetDOld()[kcnode->Id()];
        for (int dim = 0; dim < kcnode->NumDof(); ++dim)
        {
          jumptxi -= (kcnode->CoData().txi()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          jumpteta -= (kcnode->CoData().teta()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          ztxi += (kcnode->CoData().txi()[dim]) * (kcnode->MoData().lm()[dim]);
          zteta += (kcnode->CoData().teta()[dim]) * (kcnode->MoData().lm()[dim]);
          znor += (kcnode->MoData().n()[dim]) * (kcnode->MoData().lm()[dim]);
        }

        std::map<int, double> mmap = kcnode->MoData().GetM();
        std::map<int, double> mmapold = kcnode->FriData().GetMOld();

        std::map<int, double>::iterator colcurr;
        std::set<int> mnodes;

        for (colcurr = mmap.begin(); colcurr != mmap.end(); colcurr++)
          mnodes.insert(colcurr->first);

        for (colcurr = mmapold.begin(); colcurr != mmapold.end(); colcurr++)
          mnodes.insert(colcurr->first);

        std::set<int>::iterator mcurr;

        // loop over all master nodes (find adjacent ones to this stick node)
        for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %", gid);
          CONTACT::FriNode* cmnode = dynamic_cast<CONTACT::FriNode*>(mnode);
          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          std::map<int, double>::iterator mcurr;

          for (int dim = 0; dim < kcnode->NumDof(); ++dim)
          {
            jumptxi += (kcnode->CoData().txi()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
            jumpteta += (kcnode->CoData().teta()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
          }
        }  //  loop over master nodes

        // gp-wise slip !!!!!!!
        if (DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == true)
        {
          jumptxi = kcnode->FriData().jump_var()[0];
          jumpteta = 0.0;

          if (Dim() == 3) jumpteta = kcnode->FriData().jump_var()[1];
        }

        // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
        std::vector<double> sum1(Dim() - 1, 0);
        sum1[0] = ztxi + ct * jumptxi;
        if (Dim() == 3) sum1[1] = zteta + ct * jumpteta;
        if (Dim() == 2) euclidean = abs(sum1[0]);
        if (Dim() == 3) euclidean = sqrt(sum1[0] * sum1[0] + sum1[1] * sum1[1]);
      }  // if cnode == Slip

      // store C in vector
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        newCtxi[k] = euclidean * ztxi - frbound * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - frbound * (zteta + ct * jumpteta);
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        newCtxi[k] = euclidean * ztxi - (frcoeff * znor) * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - (frcoeff * znor) * (zteta + ct * jumpteta);
      }
      else
        dserror("ERROR: Friction law is neiter Tresca nor Coulomb");

      newCtxi[k] = euclidean * ztxi -
                   (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (ztxi + ct * jumptxi);
      newCteta[k] = euclidean * zteta -
                    (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (zteta + ct * jumpteta);

      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linslipLMglobal.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[1]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[1], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;

      // ********************************* TETA
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs2 = linslipLMglobal.EpetraMatrix();
      sparse_crs2->FillComplete();
      double sparse_2 = 0.0;
      int sparsenumentries2 = 0;
      int sparselength2 = sparse_crs2->NumGlobalEntries(kcnode->Dofs()[2]);
      std::vector<double> sparsevalues2(sparselength2);
      std::vector<int> sparseindices2(sparselength2);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(kcnode->Dofs()[2], sparselength2, sparsenumentries2,
          &sparsevalues2[0], &sparseindices2[0]);

      for (int h = 0; h < sparselength2; ++h)
      {
        if (sparseindices2[h] == coldof)
        {
          sparse_2 = sparsevalues2[h];
          break;
        }
        else
          sparse_2 = 0.0;
      }
      double analyt_teta = sparse_2;

      // print results (derivatives) to screen
      if (abs(newCtxi[k] - refCtxi[k]) > 1e-12)
      {
        std::cout << "SLIP LM-Deriv_xi: " << kcnode->Id() << "\t w.r.t: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newCtxi[k] - refCtxi[k]) / delta
                  << "\t analyt= " << std::setprecision(4) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newCtxi[k] - refCtxi[k]) / delta);
        if (abs(analyt_txi - (newCtxi[k] - refCtxi[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }

      // print results (derivatives) to screen
      if (abs(newCteta[k] - refCteta[k]) > 1e-12)
      {
        std::cout << "SLIP LM-Deriv_eta: " << kcnode->Id()
                  << "\t w.r.t: " << snode->Dofs()[fd % dim] << "\t FD= " << std::setprecision(4)
                  << (newCteta[k] - refCteta[k]) / delta << "\t analyt= " << std::setprecision(4)
                  << analyt_teta
                  << "\t Error= " << analyt_teta - ((newCteta[k] - refCteta[k]) / delta);
        if (abs(analyt_teta - (newCteta[k] - refCteta[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->MoData().lm()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->MoData().lm()[1] -= delta;
    }
    else
    {
      snode->MoData().lm()[2] -= delta;
    }
  }  // loop over procs slave nodes


  // ********************************************************************************************
  // global loop to apply FD scheme to all slave dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::FriNode* snode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
      coldof = snode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
      coldof = snode->Dofs()[1];
    }
    else
    {
      snode->xspatial()[2] += delta;
      coldof = snode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!node) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      double jumptxi = 0;
      double jumpteta = 0;
      double ztxi = 0;
      double zteta = 0;
      double znor = 0;
      double euclidean = 0;

      if (kcnode->FriData().Slip())
      {
        // check two versions of weighted gap
        double D = kcnode->MoData().GetD()[kcnode->Id()];
        double Dold = kcnode->FriData().GetDOld()[kcnode->Id()];

        for (int dim = 0; dim < kcnode->NumDof(); ++dim)
        {
          jumptxi -= (kcnode->CoData().txi()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          jumpteta -= (kcnode->CoData().teta()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          ztxi += (kcnode->CoData().txi()[dim]) * (kcnode->MoData().lm()[dim]);
          zteta += (kcnode->CoData().teta()[dim]) * (kcnode->MoData().lm()[dim]);
          znor += (kcnode->MoData().n()[dim]) * (kcnode->MoData().lm()[dim]);
        }

        std::map<int, double> mmap = kcnode->MoData().GetM();
        std::map<int, double> mmapold = kcnode->FriData().GetMOld();

        std::map<int, double>::iterator colcurr;
        std::set<int> mnodes;

        for (colcurr = mmap.begin(); colcurr != mmap.end(); colcurr++)
          mnodes.insert(colcurr->first);

        for (colcurr = mmapold.begin(); colcurr != mmapold.end(); colcurr++)
          mnodes.insert(colcurr->first);

        std::set<int>::iterator mcurr;

        // loop over all master nodes (find adjacent ones to this stick node)
        for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %", gid);
          CONTACT::FriNode* cmnode = dynamic_cast<CONTACT::FriNode*>(mnode);

          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          std::map<int, double>::iterator mcurr;

          for (int dim = 0; dim < kcnode->NumDof(); ++dim)
          {
            jumptxi += (kcnode->CoData().txi()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
            jumpteta += (kcnode->CoData().teta()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
          }
        }  //  loop over master nodes

        // gp-wise slip !!!!!!!
        if (DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == true)
        {
          jumptxi = kcnode->FriData().jump_var()[0];
          jumpteta = 0.0;

          if (Dim() == 3) jumpteta = kcnode->FriData().jump_var()[1];
        }

        // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
        std::vector<double> sum1(Dim() - 1, 0);
        sum1[0] = ztxi + ct * jumptxi;
        if (Dim() == 3) sum1[1] = zteta + ct * jumpteta;
        if (Dim() == 2) euclidean = abs(sum1[0]);
        if (Dim() == 3) euclidean = sqrt(sum1[0] * sum1[0] + sum1[1] * sum1[1]);

      }  // if cnode == Slip

      // store C in vector
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        newCtxi[k] = euclidean * ztxi - frbound * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - frbound * (zteta + ct * jumpteta);
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        newCtxi[k] = euclidean * ztxi - (frcoeff * znor) * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - (frcoeff * znor) * (zteta + ct * jumpteta);
      }
      else
        dserror("ERROR: Friction law is neiter Tresca nor Coulomb");

      newCtxi[k] = euclidean * ztxi -
                   (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (ztxi + ct * jumptxi);
      newCteta[k] = euclidean * zteta -
                    (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (zteta + ct * jumpteta);



      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linslipDISglobal.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[1]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[1], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;

      // ********************************* TETA
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs2 = linslipDISglobal.EpetraMatrix();
      sparse_crs2->FillComplete();
      double sparse_2 = 0.0;
      int sparsenumentries2 = 0;
      int sparselength2 = sparse_crs2->NumGlobalEntries(kcnode->Dofs()[2]);
      std::vector<double> sparsevalues2(sparselength2);
      std::vector<int> sparseindices2(sparselength2);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(kcnode->Dofs()[2], sparselength2, sparsenumentries2,
          &sparsevalues2[0], &sparseindices2[0]);

      for (int h = 0; h < sparselength2; ++h)
      {
        if (sparseindices2[h] == coldof)
        {
          sparse_2 = sparsevalues2[h];
          break;
        }
        else
          sparse_2 = 0.0;
      }
      double analyt_teta = sparse_2;


      if (abs(newCtxi[k] - refCtxi[k]) > 1e-12)
      {
        std::cout << "SLIP DIS-Deriv_xi: " << kcnode->Id()
                  << "\t w.r.t Slave: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newCtxi[k] - refCtxi[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newCtxi[k] - refCtxi[k]) / delta);
        if (abs(analyt_txi - (newCtxi[k] - refCtxi[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }

      // print results (derivatives) to screen
      if (abs(newCteta[k] - refCteta[k]) > 1e-12)
      {
        std::cout << "SLIP DIS-Deriv_eta: " << kcnode->Id()
                  << "\t w.r.t Slave: " << snode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newCteta[k] - refCteta[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_teta
                  << "\t Error= " << analyt_teta - ((newCteta[k] - refCteta[k]) / delta);
        if (abs(analyt_teta - (newCteta[k] - refCteta[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }
  }  // loop over procs slave nodes

  // ********************************************************************************************
  // global loop to apply FD scheme to all master dofs (=3*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * mnodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find master node with gid %", gid);
    CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      mnode->xspatial()[0] += delta;
      coldof = mnode->Dofs()[0];
    }
    else if (fd % dim == 1)
    {
      mnode->xspatial()[1] += delta;
      coldof = mnode->Dofs()[1];
    }
    else
    {
      mnode->xspatial()[2] += delta;
      coldof = mnode->Dofs()[2];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      double jumptxi = 0;
      double jumpteta = 0;
      double ztxi = 0;
      double zteta = 0;
      double znor = 0;
      double euclidean = 0;

      if (kcnode->FriData().Slip())
      {
        // check two versions of weighted gap
        double D = kcnode->MoData().GetD()[kcnode->Id()];
        double Dold = kcnode->FriData().GetDOld()[kcnode->Id()];

        for (int dim = 0; dim < kcnode->NumDof(); ++dim)
        {
          jumptxi -= (kcnode->CoData().txi()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          jumpteta -= (kcnode->CoData().teta()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          ztxi += (kcnode->CoData().txi()[dim]) * (kcnode->MoData().lm()[dim]);
          zteta += (kcnode->CoData().teta()[dim]) * (kcnode->MoData().lm()[dim]);
          znor += (kcnode->MoData().n()[dim]) * (kcnode->MoData().lm()[dim]);
        }

        std::map<int, double> mmap = kcnode->MoData().GetM();
        std::map<int, double> mmapold = kcnode->FriData().GetMOld();

        std::map<int, double>::iterator colcurr;
        std::set<int> mnodes;

        for (colcurr = mmap.begin(); colcurr != mmap.end(); colcurr++)
          mnodes.insert(colcurr->first);

        for (colcurr = mmapold.begin(); colcurr != mmapold.end(); colcurr++)
          mnodes.insert(colcurr->first);

        std::set<int>::iterator mcurr;

        // loop over all master nodes (find adjacent ones to this stick node)
        for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %", gid);
          CONTACT::FriNode* cmnode = dynamic_cast<CONTACT::FriNode*>(mnode);

          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          std::map<int, double>::iterator mcurr;

          for (int dim = 0; dim < kcnode->NumDof(); ++dim)
          {
            jumptxi += (kcnode->CoData().txi()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
            jumpteta += (kcnode->CoData().teta()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
          }
        }  //  loop over master nodes

        // gp-wise slip !!!!!!!
        if (DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == true)
        {
          jumptxi = kcnode->FriData().jump_var()[0];
          jumpteta = 0.0;

          if (Dim() == 3) jumpteta = kcnode->FriData().jump_var()[1];
        }

        // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
        std::vector<double> sum1(Dim() - 1, 0);
        sum1[0] = ztxi + ct * jumptxi;
        if (Dim() == 3) sum1[1] = zteta + ct * jumpteta;
        if (Dim() == 2) euclidean = abs(sum1[0]);
        if (Dim() == 3) euclidean = sqrt(sum1[0] * sum1[0] + sum1[1] * sum1[1]);

      }  // if cnode == Slip

      // store C in vector
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        newCtxi[k] = euclidean * ztxi - frbound * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - frbound * (zteta + ct * jumpteta);
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        newCtxi[k] = euclidean * ztxi - (frcoeff * znor) * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - (frcoeff * znor) * (zteta + ct * jumpteta);
      }
      else
        dserror("ERROR: Friction law is neiter Tresca nor Coulomb");

      newCtxi[k] = euclidean * ztxi -
                   (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (ztxi + ct * jumptxi);
      newCteta[k] = euclidean * zteta -
                    (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (zteta + ct * jumpteta);



      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linslipDISglobal.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[1]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[1], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;

      // ********************************* TETA
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs2 = linslipDISglobal.EpetraMatrix();
      sparse_crs2->FillComplete();
      double sparse_2 = 0.0;
      int sparsenumentries2 = 0;
      int sparselength2 = sparse_crs2->NumGlobalEntries(kcnode->Dofs()[2]);
      std::vector<double> sparsevalues2(sparselength2);
      std::vector<int> sparseindices2(sparselength2);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(kcnode->Dofs()[2], sparselength2, sparsenumentries2,
          &sparsevalues2[0], &sparseindices2[0]);

      for (int h = 0; h < sparselength2; ++h)
      {
        if (sparseindices2[h] == coldof)
        {
          sparse_2 = sparsevalues2[h];
          break;
        }
        else
          sparse_2 = 0.0;
      }
      double analyt_teta = sparse_2;

      // print results (derivatives) to screen
      if (abs(newCtxi[k] - refCtxi[k]) > 1e-12)
      {
        std::cout << "SLIP DIS-Deriv_xi: " << kcnode->Id()
                  << "\t w.r.t Master: " << mnode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newCtxi[k] - refCtxi[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newCtxi[k] - refCtxi[k]) / delta);
        if (abs(analyt_txi - (newCtxi[k] - refCtxi[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }

      if (abs(newCteta[k] - refCteta[k]) > 1e-12)
      {
        std::cout << "SLIP DIS-Deriv_eta: " << kcnode->Id()
                  << "\t w.r.t Master: " << mnode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newCteta[k] - refCteta[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_teta
                  << "\t Error= " << analyt_teta - ((newCteta[k] - refCteta[k]) / delta);
        if (abs(analyt_teta - (newCteta[k] - refCteta[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      mnode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      mnode->xspatial()[1] -= delta;
    }
    else
    {
      mnode->xspatial()[2] -= delta;
    }
  }
  // ********************************************************************************************
  // global loop to apply FD scheme to all Wear dofs (=1*nodes)
  // ********************************************************************************************
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    int coldof = 0;
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find master node with gid %", gid);
    CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      mnode->WearData().wcurr()[0] += delta;
      coldof = mnode->Dofs()[0];
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      double jumptxi = 0;
      double jumpteta = 0;
      double ztxi = 0;
      double zteta = 0;
      double znor = 0;
      double euclidean = 0;

      if (kcnode->FriData().Slip())
      {
        // check two versions of weighted gap
        double D = kcnode->MoData().GetD()[kcnode->Id()];
        double Dold = kcnode->FriData().GetDOld()[kcnode->Id()];

        for (int dim = 0; dim < kcnode->NumDof(); ++dim)
        {
          jumptxi -= (kcnode->CoData().txi()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          jumpteta -= (kcnode->CoData().teta()[dim]) * (D - Dold) * (kcnode->xspatial()[dim]);
          ztxi += (kcnode->CoData().txi()[dim]) * (kcnode->MoData().lm()[dim]);
          zteta += (kcnode->CoData().teta()[dim]) * (kcnode->MoData().lm()[dim]);
          znor += (kcnode->MoData().n()[dim]) * (kcnode->MoData().lm()[dim]);
        }

        std::map<int, double> mmap = kcnode->MoData().GetM();
        std::map<int, double> mmapold = kcnode->FriData().GetMOld();

        std::map<int, double>::iterator colcurr;
        std::set<int> mnodes;

        for (colcurr = mmap.begin(); colcurr != mmap.end(); colcurr++)
          mnodes.insert(colcurr->first);

        for (colcurr = mmapold.begin(); colcurr != mmapold.end(); colcurr++)
          mnodes.insert(colcurr->first);

        std::set<int>::iterator mcurr;

        // loop over all master nodes (find adjacent ones to this stick node)
        for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %", gid);
          CONTACT::FriNode* cmnode = dynamic_cast<CONTACT::FriNode*>(mnode);

          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          std::map<int, double>::iterator mcurr;

          for (int dim = 0; dim < kcnode->NumDof(); ++dim)
          {
            jumptxi += (kcnode->CoData().txi()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
            jumpteta += (kcnode->CoData().teta()[dim]) * (mik - mikold) * (cmnode->xspatial()[dim]);
          }
        }  //  loop over master nodes

        // gp-wise slip !!!!!!!
        if (DRT::INPUT::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == true)
        {
          jumptxi = kcnode->FriData().jump_var()[0];
          jumpteta = 0.0;

          if (Dim() == 3) jumpteta = kcnode->FriData().jump_var()[1];
        }

        // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
        std::vector<double> sum1(Dim() - 1, 0);
        sum1[0] = ztxi + ct * jumptxi;
        if (Dim() == 3) sum1[1] = zteta + ct * jumpteta;
        if (Dim() == 2) euclidean = abs(sum1[0]);
        if (Dim() == 3) euclidean = sqrt(sum1[0] * sum1[0] + sum1[1] * sum1[1]);

      }  // if cnode == Slip

      // store C in vector
      if (ftype == INPAR::CONTACT::friction_tresca)
      {
        newCtxi[k] = euclidean * ztxi - frbound * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - frbound * (zteta + ct * jumpteta);
      }
      else if (ftype == INPAR::CONTACT::friction_coulomb)
      {
        newCtxi[k] = euclidean * ztxi - (frcoeff * znor) * (ztxi + ct * jumptxi);
        newCteta[k] = euclidean * zteta - (frcoeff * znor) * (zteta + ct * jumpteta);
      }
      else
        dserror("ERROR: Friction law is neiter Tresca nor Coulomb");

      newCtxi[k] = euclidean * ztxi -
                   (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (ztxi + ct * jumptxi);
      newCteta[k] = euclidean * zteta -
                    (frcoeff * (znor - cn * kcnode->CoData().Getg())) * (zteta + ct * jumpteta);


      // ************************************************************************
      // Extract linearizations from sparse matrix !!!
      // ************************************************************************

      // ********************************* TXI
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = linslipWglobal.EpetraMatrix();
      sparse_crs->FillComplete();
      double sparse_ij = 0.0;
      int sparsenumentries = 0;
      int sparselength = sparse_crs->NumGlobalEntries(kcnode->Dofs()[1]);
      std::vector<double> sparsevalues(sparselength);
      std::vector<int> sparseindices(sparselength);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(
          kcnode->Dofs()[1], sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);

      for (int h = 0; h < sparselength; ++h)
      {
        if (sparseindices[h] == coldof)
        {
          sparse_ij = sparsevalues[h];
          break;
        }
        else
          sparse_ij = 0.0;
      }
      double analyt_txi = sparse_ij;

      // ********************************* TETA
      Teuchos::RCP<Epetra_CrsMatrix> sparse_crs2 = linslipWglobal.EpetraMatrix();
      sparse_crs2->FillComplete();
      double sparse_2 = 0.0;
      int sparsenumentries2 = 0;
      int sparselength2 = sparse_crs2->NumGlobalEntries(kcnode->Dofs()[2]);
      std::vector<double> sparsevalues2(sparselength2);
      std::vector<int> sparseindices2(sparselength2);
      // int sparseextractionstatus =
      sparse_crs->ExtractGlobalRowCopy(kcnode->Dofs()[2], sparselength2, sparsenumentries2,
          &sparsevalues2[0], &sparseindices2[0]);

      for (int h = 0; h < sparselength2; ++h)
      {
        if (sparseindices2[h] == coldof)
        {
          sparse_2 = sparsevalues2[h];
          break;
        }
        else
          sparse_2 = 0.0;
      }
      double analyt_teta = sparse_2;

      // print results (derivatives) to screen
      if (abs(newCtxi[k] - refCtxi[k]) > 1e-12)
      {
        std::cout << "SLIP W-Deriv_xi: " << kcnode->Id()
                  << "\t w.r.t Master: " << mnode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newCtxi[k] - refCtxi[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_txi
                  << "\t Error= " << analyt_txi - ((newCtxi[k] - refCtxi[k]) / delta);
        if (abs(analyt_txi - (newCtxi[k] - refCtxi[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }

      if (abs(newCteta[k] - refCteta[k]) > 1e-12)
      {
        std::cout << "SLIP W-Deriv_eta: " << kcnode->Id()
                  << "\t w.r.t Master: " << mnode->Dofs()[fd % dim]
                  << "\t FD= " << std::setprecision(4) << (newCteta[k] - refCteta[k]) / delta
                  << "\t analyt= " << std::setprecision(5) << analyt_teta
                  << "\t Error= " << analyt_teta - ((newCteta[k] - refCteta[k]) / delta);
        if (abs(analyt_teta - (newCteta[k] - refCteta[k]) / delta) > 1.0e-4)
          std::cout << "*** WARNING ***" << std::endl;
        else
          std::cout << " " << std::endl;
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      mnode->WearData().wcurr()[0] -= delta;
    }
  }

  // back to normal...
  Initialize();
  Evaluate();

  return;

}  // FDCheckSlipTrescaDeriv

/*----------------------------------------------------------------------*
 | Finite difference check for T-Mortar derivatives          farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckMortarTDeriv()
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for D-Matrix entries
  std::map<int, std::map<int, double>> refT;  // stores dof-wise the entries of D
  std::map<int, std::map<int, double>> newT;

  std::map<int, std::map<int, std::map<int, double>>>
      refDerivT;  // stores old derivm for every node

  // problem dimension (2D or 3D)
  int dim = Dim();

  // print reference to screen (D-derivative-maps) and store them for later comparison
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    // typedef std::map<int,std::map<int,double> >::const_iterator CID;
    // typedef std::map<int,double>::const_iterator CI;

    if ((int)(cnode->WearData().GetT().size()) == 0) continue;

    //    for( int d=0; d<dim; d++ )
    //    {
    int dof = cnode->Dofs()[0];
    refT[dof] = cnode->WearData().GetT()[0];
    //    }

    refDerivT[gid] = cnode->WearData().GetDerivTw();
  }

  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      if ((int)(kcnode->WearData().GetT().size()) == 0) continue;

      typedef std::map<int, double>::const_iterator CI;

      for (int d = 0; d < 1; d++)
      {
        int dof = kcnode->Dofs()[d];

        // store D-values into refD
        newT[dof] = kcnode->WearData().GetT()[d];

        // print results (derivatives) to screen
        for (CI p = newT[dof].begin(); p != newT[dof].end(); ++p)
        {
          if (abs(newT[dof][p->first] - refT[dof][p->first]) > 1e-12)
          {
            double finit = (newT[dof][p->first] - refT[dof][p->first]) / delta;
            double analy = ((refDerivT[kgid])[(p->first) / Dim()])[sdof];
            double dev = finit - analy;

            // kgid: currently tested dof of slave node kgid
            // (p->first)/Dim(): paired master
            // sdof: currently modified slave dof
            std::cout << "(" << dof << "," << (p->first) << "," << sdof << ") : fd=" << finit
                      << " derivT=" << analy << " DEVIATION " << dev;

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // contents of Evaluate()
  // *******************************************************************
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for T-Mortar derivatives (Master) farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckMortarT_Master_Deriv()
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for D-Matrix entries
  std::map<int, std::map<int, double>> refT;  // stores dof-wise the entries of D
  std::map<int, std::map<int, double>> newT;

  std::map<int, std::map<int, std::map<int, double>>>
      refDerivT;  // stores old derivm for every node

  // problem dimension (2D or 3D)
  int dim = Dim();

  // print reference to screen (D-derivative-maps) and store them for later comparison
  // loop over proc's slave nodes
  for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
  {
    int gid = mnoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    // typedef std::map<int,std::map<int,double> >::const_iterator CID;
    // typedef std::map<int,double>::const_iterator CI;

    if ((int)(cnode->WearData().GetT().size()) == 0) continue;

    int dof = cnode->Dofs()[0];
    refT[dof] = cnode->WearData().GetT()[0];

    refDerivT[gid] = cnode->WearData().GetDerivTw();
  }

  //**************************************************************
  //                            SLAVE
  //**************************************************************
  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < mnoderowmap_->NumMyElements(); ++k)
    {
      int kgid = mnoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      if ((int)(kcnode->WearData().GetT().size()) == 0) continue;

      typedef std::map<int, double>::const_iterator CI;

      for (int d = 0; d < 1; d++)
      {
        int dof = kcnode->Dofs()[d];

        // store D-values into refD
        newT[dof] = kcnode->WearData().GetT()[d];

        // print results (derivatives) to screen
        for (CI p = newT[dof].begin(); p != newT[dof].end(); ++p)
        {
          if (abs(newT[dof][p->first] - refT[dof][p->first]) > 1e-12)
          {
            double finit = (newT[dof][p->first] - refT[dof][p->first]) / delta;
            double analy = ((refDerivT[kgid])[(p->first) / Dim()])[sdof];
            double dev = finit - analy;

            // kgid: currently tested dof of slave node kgid
            // (p->first)/Dim(): paired master
            // sdof: currently modified slave dof
            std::cout << "(" << dof << "," << (p->first) << "," << sdof << ") : fd=" << finit
                      << " derivT=" << analy << " DEVIATION " << dev;

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  //**************************************************************
  //                            MASTER
  //**************************************************************
  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd = 0; fd < dim * mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];

    std::cout << "\nDERIVATIVE FOR M-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < mnoderowmap_->NumMyElements(); ++k)
    {
      int kgid = mnoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      if ((int)(kcnode->WearData().GetT().size()) == 0) continue;

      typedef std::map<int, double>::const_iterator CI;

      for (int d = 0; d < 1; d++)
      {
        int dof = kcnode->Dofs()[d];

        // store D-values into refD
        newT[dof] = kcnode->WearData().GetT()[d];

        // print results (derivatives) to screen
        for (CI p = newT[dof].begin(); p != newT[dof].end(); ++p)
        {
          if (abs(newT[dof][p->first] - refT[dof][p->first]) > 1e-12)
          {
            double finit = (newT[dof][p->first] - refT[dof][p->first]) / delta;
            double analy = ((refDerivT[kgid])[(p->first) / Dim()])[sdof];
            double dev = finit - analy;

            // kgid: currently tested dof of slave node kgid
            // (p->first)/Dim(): paired master
            // sdof: currently modified slave dof
            std::cout << "(" << dof << "," << (p->first) << "," << sdof << ") : fd=" << finit
                      << " derivT=" << analy << " DEVIATION " << dev;

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // contents of Evaluate()
  // *******************************************************************
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for E-Mortar derivatives          farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckMortarEDeriv()
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for D-Matrix entries
  std::map<int, std::map<int, double>> refE;  // stores dof-wise the entries of D
  std::map<int, std::map<int, double>> newE;

  std::map<int, std::map<int, std::map<int, double>>>
      refDerivE;  // stores old derivm for every node

  // problem dimension (2D or 3D)
  int dim = Dim();

  // print reference to screen (D-derivative-maps) and store them for later comparison
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    // typedef std::map<int,std::map<int,double> >::const_iterator CID;
    // typedef std::map<int,double>::const_iterator CI;

    if ((int)(cnode->WearData().GetE().size()) == 0) continue;

    //    for( int d=0; d<dim; d++ )
    //    {
    int dof = cnode->Dofs()[0];
    refE[dof] = cnode->WearData().GetE()[0];
    //    }

    refDerivE[gid] = cnode->WearData().GetDerivE();
  }

  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      if ((int)(kcnode->WearData().GetE().size()) == 0) continue;

      typedef std::map<int, double>::const_iterator CI;

      for (int d = 0; d < 1; d++)
      {
        int dof = kcnode->Dofs()[d];

        // store D-values into refD
        newE[dof] = kcnode->WearData().GetE()[d];

        // print results (derivatives) to screen
        for (CI p = newE[dof].begin(); p != newE[dof].end(); ++p)
        {
          if (abs(newE[dof][p->first] - refE[dof][p->first]) > 1e-12)
          {
            double finit = (newE[dof][p->first] - refE[dof][p->first]) / delta;
            double analy = ((refDerivE[kgid])[(p->first) / Dim()])[sdof];
            double dev = finit - analy;

            // kgid: currently tested dof of slave node kgid
            // (p->first)/Dim(): paired master
            // sdof: currently modified slave dof
            std::cout << "(" << dof << "," << (p->first) / Dim() << "," << sdof
                      << ") : fd=" << finit << " derivE=" << analy << " DEVIATION " << dev;

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // contents of Evaluate()
  // *******************************************************************
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for E-Mortar derivatives (Master) farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckMortarE_Master_Deriv()
{
  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for D-Matrix entries
  std::map<int, std::map<int, double>> refE;  // stores dof-wise the entries of D
  std::map<int, std::map<int, double>> newE;

  std::map<int, std::map<int, std::map<int, double>>>
      refDerivE;  // stores old derivm for every node

  // problem dimension (2D or 3D)
  int dim = Dim();

  // print reference to screen (D-derivative-maps) and store them for later comparison
  // loop over proc's slave nodes
  for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
  {
    int gid = mnoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    // typedef std::map<int,std::map<int,double> >::const_iterator CID;
    // typedef std::map<int,double>::const_iterator CI;

    if ((int)(cnode->WearData().GetE().size()) == 0) continue;

    int dof = cnode->Dofs()[0];
    refE[dof] = cnode->WearData().GetE()[0];

    refDerivE[gid] = cnode->WearData().GetDerivE();
  }

  //**************************************************************
  //                            SLAVE
  //**************************************************************
  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < mnoderowmap_->NumMyElements(); ++k)
    {
      int kgid = mnoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      if ((int)(kcnode->WearData().GetE().size()) == 0) continue;

      typedef std::map<int, double>::const_iterator CI;

      for (int d = 0; d < 1; d++)
      {
        int dof = kcnode->Dofs()[d];

        // store D-values into refD
        newE[dof] = kcnode->WearData().GetE()[d];

        // print results (derivatives) to screen
        for (CI p = newE[dof].begin(); p != newE[dof].end(); ++p)
        {
          if (abs(newE[dof][p->first] - refE[dof][p->first]) > 1e-12)
          {
            double finit = (newE[dof][p->first] - refE[dof][p->first]) / delta;
            double analy = ((refDerivE[kgid])[(p->first) / Dim()])[sdof];
            double dev = finit - analy;

            // kgid: currently tested dof of slave node kgid
            // (p->first)/Dim(): paired master
            // sdof: currently modified slave dof
            std::cout << "(" << dof << "," << (p->first) / Dim() << "," << sdof
                      << ") : fd=" << finit << " derivE=" << analy << " DEVIATION " << dev;

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }


  //**************************************************************
  //                            MASTER
  //**************************************************************
  // global loop to apply FD scheme to all SLAVE dofs (=dim*nodes)
  for (int fd = 0; fd < dim * mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];

    std::cout << "\nDERIVATIVE FOR M-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < mnoderowmap_->NumMyElements(); ++k)
    {
      int kgid = mnoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      if ((int)(kcnode->WearData().GetE().size()) == 0) continue;

      typedef std::map<int, double>::const_iterator CI;

      for (int d = 0; d < 1; d++)
      {
        int dof = kcnode->Dofs()[d];

        // store D-values into refD
        newE[dof] = kcnode->WearData().GetE()[d];

        // print results (derivatives) to screen
        for (CI p = newE[dof].begin(); p != newE[dof].end(); ++p)
        {
          if (abs(newE[dof][p->first] - refE[dof][p->first]) > 1e-12)
          {
            double finit = (newE[dof][p->first] - refE[dof][p->first]) / delta;
            double analy = ((refDerivE[kgid])[(p->first) / Dim()])[sdof];
            double dev = finit - analy;

            // kgid: currently tested dof of slave node kgid
            // (p->first)/Dim(): paired master
            // sdof: currently modified slave dof
            std::cout << "(" << dof << "," << (p->first) / Dim() << "," << sdof
                      << ") : fd=" << finit << " derivE=" << analy << " DEVIATION " << dev;

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  // back to normal...

  // Initialize
  Initialize();

  // compute element areas
  SetElementAreas();

  // *******************************************************************
  // contents of Evaluate()
  // *******************************************************************
  Evaluate();

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for normal Wear-lm derivatives    farah 07/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckWearDerivLm()
{
  double wcoeff = InterfaceParams().get<double>("WEARCOEFF");

  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for gap values
  int nrow = snoderowmap_->NumMyElements();
  std::vector<double> refW(nrow);
  std::vector<double> newW(nrow);

  // problem dimension (2D or 3D)
  int dim = Dim();

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    // store wear-values into refW
    refW[i] = cnode->WearData().WeightedWear();
  }

  // global loop to apply FD scheme to all slave dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize -- set WearData.Wear() = 0 !!!
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->MoData().lm()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->MoData().lm()[1] += delta;
    }
    else
    {
      snode->MoData().lm()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();


    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      // store gap-values into newG
      newW[k] = wcoeff * kcnode->WearData().DeltaWeightedWear();

      //      if (abs(newW[k]-refW[k]) > 1e-12 && newW[k]!=1.0e12 && refW[k] != 1.0e12)
      //      {
      double finit = (newW[k] - refW[k]) / delta;
      double analy = kcnode->CoData().GetDerivWlm()[snode->Dofs()[fd % dim]];
      double dev = finit - analy;

      // kgid: id of currently tested slave node
      // snode->Dofs()[fd%dim]: currently modified slave dof
      std::cout << "(" << kgid << "," << snode->Dofs()[fd % dim] << ") : fd=" << finit
                << " derivw=" << analy << " DEVIATION " << dev;

      if (abs(dev) > 1e-4)
      {
        std::cout << " ***** WARNING ***** ";
        w++;
      }
      else if (abs(dev) > 1e-5)
      {
        std::cout << " ***** warning ***** ";
        w++;
      }

      std::cout << std::endl;
      //      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->MoData().lm()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->MoData().lm()[1] -= delta;
    }
    else
    {
      snode->MoData().lm()[2] -= delta;
    }


    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
    //    if (w>0)
    //      dserror("WARNING!!!");
  }

  // LM only on slave nodes!!!

  // back to normal...
  Initialize();
  Evaluate();

  // write the init. wear back into the node
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    // store gap-values into refG
    cnode->WearData().WeightedWear() = refW[i];
  }


  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for normal Wear-d derivatives     farah 07/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::FDCheckWearDeriv()
{
  double wcoeff = InterfaceParams().get<double>("WEARCOEFF");

  // FD checks only for serial case
  Teuchos::RCP<Epetra_Map> snodefullmap = LINALG::AllreduceEMap(*snoderowmap_);
  Teuchos::RCP<Epetra_Map> mnodefullmap = LINALG::AllreduceEMap(*mnoderowmap_);
  if (Comm().NumProc() > 1) dserror("ERROR: FD checks only for serial case");

  // get out of here if not participating in interface
  if (!lComm()) return;

  // create storage for gap values
  int nrow = snoderowmap_->NumMyElements();
  std::vector<double> refW(nrow);
  std::vector<double> newW(nrow);

  // problem dimension (2D or 3D)
  int dim = Dim();

  // store reference
  // loop over proc's slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    // store wear-values into refW
    refW[i] = cnode->WearData().WeightedWear();
  }

  // global loop to apply FD scheme to all slave dofs (=dim*nodes)
  for (int fd = 0; fd < dim * snodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize -- set WearData.Wear() = 0 !!!
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(node);

    int sdof = snode->Dofs()[fd % dim];
    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd % dim == 0)
    {
      snode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();


    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      // store gap-values into newG
      newW[k] = wcoeff * kcnode->WearData().DeltaWeightedWear();

      if (abs(newW[k] - refW[k]) > 1e-12 && newW[k] != 1.0e12 && refW[k] != 1.0e12)
      {
        double finit = (newW[k] - refW[k]) / delta;
        double analy = kcnode->CoData().GetDerivW()[snode->Dofs()[fd % dim]];
        double dev = finit - analy;

        // kgid: id of currently tested slave node
        // snode->Dofs()[fd%dim]: currently modified slave dof
        std::cout << "(" << kgid << "," << snode->Dofs()[fd % dim] << ") : fd=" << finit
                  << " derivw=" << analy << " DEVIATION " << dev;

        if (abs(dev) > 1e-4)
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if (abs(dev) > 1e-5)
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }
    }
    // undo finite difference modification
    if (fd % dim == 0)
    {
      snode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      snode->xspatial()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
    //    if (w>0)
    //      dserror("WARNING!!!");
  }

  // global loop to apply FD scheme to all master dofs (=dim*nodes)
  for (int fd = 0; fd < dim * mnodefullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap->GID(fd / dim);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find master node with gid %", gid);
    CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(node);

    int mdof = mnode->Dofs()[fd % dim];
    std::cout << "\nDERIVATIVE FOR M-NODE # " << gid << " DOF: " << mdof << std::endl;

    // do step forward (modify nodal displacement)
    double delta = 1e-6;
    if (fd % dim == 0)
    {
      mnode->xspatial()[0] += delta;
    }
    else if (fd % dim == 1)
    {
      mnode->xspatial()[1] += delta;
    }
    else
    {
      mnode->xspatial()[2] += delta;
    }

    // compute element areas
    SetElementAreas();

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    Evaluate();

    // compute finite difference derivative
    for (int k = 0; k < snoderowmap_->NumMyElements(); ++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CONTACT::FriNode* kcnode = dynamic_cast<CONTACT::FriNode*>(knode);

      // store gap-values into newG
      newW[k] = wcoeff * kcnode->WearData().DeltaWeightedWear();

      if (abs(newW[k] - refW[k]) > 1e-12 && newW[k] != 1.0e12 && refW[k] != 1.0e12)
      {
        double finit = (newW[k] - refW[k]) / delta;
        double analy = kcnode->CoData().GetDerivW()[mnode->Dofs()[fd % dim]];
        double dev = finit - analy;

        // kgid: id of currently tested slave node
        // mnode->Dofs()[fd%dim]: currently modified slave dof
        std::cout << "(" << kgid << "," << mnode->Dofs()[fd % dim] << ") : fd=" << finit
                  << " derivw=" << analy << " DEVIATION " << dev;

        if (abs(dev) > 1e-4)
        {
          std::cout << " ***** WARNING ***** ";
          w++;
        }
        else if (abs(dev) > 1e-5)
        {
          std::cout << " ***** warning ***** ";
          w++;
        }

        std::cout << std::endl;
      }
    }

    // undo finite difference modification
    if (fd % dim == 0)
    {
      mnode->xspatial()[0] -= delta;
    }
    else if (fd % dim == 1)
    {
      mnode->xspatial()[1] -= delta;
    }
    else
    {
      mnode->xspatial()[2] -= delta;
    }

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
    //    if (w>0)
    //      dserror("WARNING!!!");
  }

  // back to normal...
  Initialize();
  Evaluate();

  // write the init. wear back into the node
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    cnode->WearData().WeightedWear() = refW[i];
  }

  return;
}
