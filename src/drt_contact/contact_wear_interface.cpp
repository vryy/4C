/*!----------------------------------------------------------------------
\file contact_wear_interface.H

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

Questions? Contact Prof. Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    farah 09/13|
 *----------------------------------------------------------------------*/
#include "contact_defines.H"
#include "contact_wear_interface.H"
#include "contact_interface.H"
#include "contact_node.H"
#include "contact_element.H"
#include "friction_node.H"

#include "../drt_mortar/mortar_dofset.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_mortar/mortar_element.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 09/13|
 *----------------------------------------------------------------------*/
CONTACT::WearInterface::WearInterface(const int id, const Epetra_Comm& comm,
                                  const int dim,
                                  const Teuchos::ParameterList& icontact,
                                  bool selfcontact,
                                  INPAR::MORTAR::RedundantStorage redundant) :
CONTACT::CoInterface(id,comm,dim,icontact,selfcontact,redundant)
{
  // set wear contact status
  INPAR::CONTACT::WearType wtype = DRT::INPUT::IntegralValue<INPAR::CONTACT::WearType>(icontact,"WEARTYPE");
  if (wtype == INPAR::CONTACT::wear_impl)
    wearimpl_ = true;
  else
    wearimpl_=false;

  // set wear contact discretization
  if (wtype == INPAR::CONTACT::wear_discr)
    weardiscr_ = true;
  else
    weardiscr_=false;


  return;
}

/*----------------------------------------------------------------------*
 |  Assemble Mortar wear matrices                            farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleTE(LINALG::SparseMatrix& tglobal,
                                      LINALG::SparseMatrix& eglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (slipnodes_==Teuchos::null)
    return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* fnode = static_cast<FriNode*>(node);

    if (fnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleTE: Node ownership inconsistency!");

    /**************************************************** T-matrix ******/
    if ((fnode->FriDataPlus().GetT()).size()>0)
    {
      std::vector<std::map<int,double> >& tmap = fnode->FriDataPlus().GetT();
      int colsize = (int)tmap[0].size();

      std::map<int,double>::iterator colcurr;

      int row = fnode->Dofs()[0];
      int k = 0;

      for (colcurr=tmap[0].begin();colcurr!=tmap[0].end();++colcurr)
      {
        int col = colcurr->first;
        double val = colcurr->second;

        // don't check for diagonality
        // since for standard shape functions, as in general when using
        // arbitrary shape function types, this is not the case

        // create the d matrix, do not assemble zeros
        if (abs(val)>1.0e-12) tglobal.Assemble(val, row, col);

        ++k;
      }

      if (k!=colsize)
        dserror("ERROR: AssembleTE: k = %i but colsize = %i",k,colsize);

    }

    /**************************************************** E-matrix ******/
    if ((fnode->FriDataPlus().GetE()).size()>0)
    {
      std::vector<std::map<int,double> >& emap = fnode->FriDataPlus().GetE();
      int rowsize = 1;//fnode->NumDof();
      int colsize = (int)emap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)emap[j].size() != (int)emap[j+1].size())
          dserror("ERROR: AssembleTE: Column dim. of nodal E-map is inconsistent!");

      std::map<int,double>::iterator colcurr;

      int row = fnode->Dofs()[0];
      int k = 0;

      for (colcurr=emap[0].begin();colcurr!=emap[0].end();++colcurr)
      {
        int col = colcurr->first;
        double val = colcurr->second;

        // do not assemble zeros into m matrix
        if (abs(val)>1.0e-12) eglobal.Assemble(val,row,col);
        ++k;
      }

      if (k!=colsize)
        dserror("ERROR: AssembleTE: k = %i but colsize = %i",k,colsize);

    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinT containing disp derivatives         farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinT_D(LINALG::SparseMatrix& lintglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (slipnodes_==Teuchos::null)
    return;

  /**********************************************************************/
  // we have: T_wj,c with j = Lagrange multiplier slave dof
  //                 with w = wear slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinT)_kc = T_wj,c * z_j
  /**********************************************************************/

  for (int j=0;j<slipnodes_->NumMyElements();++j)
  {
    int gid = slipnodes_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* fnode = static_cast<FriNode*>(node);

    // Mortar matrix Tw derivatives
    std::map<int,std::map<int,double> >& tderiv = fnode->FriDataPlus().GetDerivTw();

    // get sizes and iterator start
    int slavesize = (int)tderiv.size(); // column size
    std::map<int,std::map<int,double> >::iterator scurr = tderiv.begin();

    /********************************************** LinTMatrix **********/
    // loop over all DISP slave nodes in the DerivT-map of the current LM slave node
    for (int k=0;k<slavesize;++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = idiscret_->gNode(sgid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
      FriNode* csnode = static_cast<FriNode*>(snode);

      // current Lagrange multipliers
      double lmn = 0.0;
      if (Dim()==2)
        lmn=(csnode->MoData().lm()[0]) * (csnode->MoData().n()[0]) + (csnode->MoData().lm()[1]) * (csnode->MoData().n()[1]);
      else if (Dim()==3)
        lmn=(csnode->MoData().lm()[0]) * (csnode->MoData().n()[0]) + (csnode->MoData().lm()[1]) * (csnode->MoData().n()[1]) + (csnode->MoData().lm()[2]) * (csnode->MoData().n()[2]);
      else
        dserror("False Dimension!");

      // Mortar matrix T derivatives
      std::map<int,double>& thisdderive = fnode->FriDataPlus().GetDerivTw()[sgid];
      int mapsize = (int)(thisdderive.size());

      // we choose the first node dof as wear dof
      int row = fnode->Dofs()[0];
      std::map<int,double>::iterator scolcurr = thisdderive.begin();

      // loop over all directional derivative entries
      for (int c=0;c<mapsize;++c)
      {
        int col = scolcurr->first;
        double val = lmn * (scolcurr->second);
        ++scolcurr;

        // owner of LM slave node can do the assembly, although it actually
        // might not own the corresponding rows in lindglobal (DISP slave node)
        // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
        //std::cout << "Assemble LinE: " << row << " " << col << " " << val << std::endl;
        if (abs(val)>1.0e-12) lintglobal.FEAssemble(val,row,col);
      }

      // check for completeness of DerivD-Derivatives-iteration
      if (scolcurr!=thisdderive.end())
        dserror("ERROR: AssembleLinE_D: Not all derivative entries of DerivE considered!");
    }
    // check for completeness of DerivD-Slave-iteration
    if (scurr!=tderiv.end())
      dserror("ERROR: AssembleLinE_D: Not all DISP slave entries of DerivE considered!");
    /******************************** Finished with LinTmatrix for delta T **********/
  }

  // *******************************************************************************
  //            Considering linearization of nodal normal vectors                 //
  // *******************************************************************************
  // loop over all LM slave nodes (row map)
  for (int j=0;j<slipnodes_->NumMyElements();++j)
  {
    int gid = slipnodes_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* fnode = static_cast<FriNode*>(node);

    if (fnode->FriDataPlus().GetT().size()>0)
    {
      // map iterator
      typedef std::map<int,double>::const_iterator CI;
      std::map<int,double>& nmap = fnode->FriDataPlus().GetT()[0];

      // loop over col entries
      for (CI z=nmap.begin();z!=nmap.end();++z)
      {
        //std::cout << "t-irst= " << z->first << std::endl;
        int gid3 = (int)((z->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        FriNode* csnode = static_cast<FriNode*>(snode);

        for (int u=0;u<Dim();++u)
        {
          std::map<int,double>& numap = csnode->CoData().GetDerivN()[u];
          double lmu = csnode->MoData().lm()[u];

          // multiply T-column entry with lin n*lambda
          for (CI b=numap.begin(); b!=numap.end(); ++b)
          {
            int row     =   fnode->Dofs()[0];
            int col     =   (b->first);
            double val  =   (z->second) * (b->second) * lmu;
            // owner of LM slave node can do the assembly, although it actually
            // might not own the corresponding rows in lindglobal (DISP slave node)
            // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
            //std::cout << "Assemble LinT N: " << row << " " << col << " " << val << std::endl;
            if (abs(val)>1.0e-12) lintglobal.FEAssemble(val,row,col);
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinE containing disp derivatives         farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinE_D(LINALG::SparseMatrix& lineglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements()==0)
    return;

  /**********************************************************************/
  // we have: E_wj,c with j = Lagrange multiplier slave dof
  //                 with w = wear slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinE)_kc = T_wj,c * z_j
  /**********************************************************************/

  // loop over all LM slave nodes (row map)
  for (int j=0;j<slipnodes_->NumMyElements();++j)
  {
    int gid = slipnodes_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* fnode = static_cast<FriNode*>(node);

    // Mortar matrix Tw derivatives
    std::map<int,std::map<int,double> >& ederiv = fnode->FriDataPlus().GetDerivE();

    // get sizes and iterator start
    int slavesize = (int)ederiv.size(); // column size
    std::map<int,std::map<int,double> >::iterator scurr = ederiv.begin();

    /********************************************** LinTMatrix **********/
    // loop over all DISP slave nodes in the DerivT-map of the current LM slave node
    for (int k=0;k<slavesize;++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = idiscret_->gNode(sgid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
      FriNode* csnode = static_cast<FriNode*>(snode);

      // current Lagrange multipliers
      double w = 0.0;
      w=csnode->FriDataPlus().wcurr()[0];

      // Mortar matrix T derivatives
      std::map<int,double>& thisdderive = fnode->FriDataPlus().GetDerivE()[sgid];
      int mapsize = (int)(thisdderive.size());

      // we choose the first node dof as wear dof
      int row = fnode->Dofs()[0];//csnode->Dofs()[0];
      std::map<int,double>::iterator scolcurr = thisdderive.begin();

      // loop over all directional derivative entries
      for (int c=0;c<mapsize;++c)
      {
        int col = scolcurr->first;
        double val = w * (scolcurr->second);
        ++scolcurr;

        // owner of LM slave node can do the assembly, although it actually
        // might not own the corresponding rows in lindglobal (DISP slave node)
        // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
        //std::cout << "Assemble LinE: " << row << " " << col << " " << val << std::endl;
        if (abs(val)>1.0e-12) lineglobal.FEAssemble(val,row,col);
      }

      // check for completeness of DerivD-Derivatives-iteration
      if (scolcurr!=thisdderive.end())
        dserror("ERROR: AssembleLinE_D: Not all derivative entries of DerivE considered!");
    }
    // check for completeness of DerivD-Slave-iteration
    if (scurr!=ederiv.end())
      dserror("ERROR: AssembleLinE_D: Not all DISP slave entries of DerivE considered!");
    /******************************** Finished with LinTmatrix for delta T **********/
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinT containing lm derivatives           farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinT_LM(LINALG::SparseMatrix& lintglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (slipnodes_==Teuchos::null)
    return;

  typedef std::map<int,double>::const_iterator CI;

  // loop over all LM slave nodes (row map)
  for (int j=0;j<slipnodes_->NumMyElements();++j)
  {
    int gid = slipnodes_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* fnode = static_cast<FriNode*>(node);

    typedef std::map<int,double>::const_iterator CI;

    if (fnode->FriDataPlus().GetT().size()>0)
    {
      // column entries for row f
      std::map<int,double>& fmap = fnode->FriDataPlus().GetT()[0];

      for (CI p=fmap.begin();p!=fmap.end();++p)
      {
        int gid2 = (int)((p->first)/Dim());
        DRT::Node* node2 = idiscret_->gNode(gid2);
        if (!node2) dserror("ERROR: Cannot find node with gid %",gid2);
        FriNode* jnode = static_cast<FriNode*>(node2);

        for (int iter=0;iter<Dim();++iter)
        {
          double n= 0.0;
          n=jnode->MoData().n()[iter];

          int row= fnode->Dofs()[0];
          int col = jnode->Dofs()[iter];
          double val = n * (p->second);
          //std::cout << "Assemble LinT: " << row << " " << col << " " << val << std::endl;
          if (abs(val)>1.0e-12) lintglobal.FEAssemble(val,row,col);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ derivatives          farah 09/13|
 |  PS: "AssembleS" is an outdated name which could make                |
 |  you confused.                                                       |
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleS(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    // prepare assembly
    std::map<int,double>& dgmap = cnode->CoData().GetDerivG();
    std::map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;
      if (DRT::INPUT::IntegralValue<int>(imortar_,"LM_NODAL_SCALE")==true &&
          cnode->MoData().GetScale() != 0.0)
        val /= cnode->MoData().GetScale();
      //std::cout << "Assemble S: " << row << " " << col << " " << val << std::endl;
      // do not assemble zeros into s matrix
      if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }
    if (DRT::INPUT::IntegralValue<int>(imortar_,"LM_NODAL_SCALE")==true &&
        cnode->MoData().GetScale() != 0.0)
    {
      std::map<int,double>& dscalemap = cnode->CoData().GetDerivScale();
      double scalefac = cnode->MoData().GetScale();
      for (colcurr=dscalemap.begin(); colcurr!=dscalemap.end(); ++colcurr)
      {
        int col = colcurr->first;
        double val = -cnode->CoData().Getg()*(colcurr->second) / (scalefac*scalefac);
        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
      }
    }

    /*************************************************************************
     * Wear implicit linearization   --> obviously, we need a new linear.    *
     *************************************************************************/
    if(wearimpl_)
    {
      // prepare assembly
      std::map<int,double>& dwmap = cnode->CoData().GetDerivW();

      for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
      {
        int col = colcurr->first;
        double val = colcurr->second;
        //std::cout << "Assemble S: " << row << " " << col << " " << val << endl;
        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
      }
    }
  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ lm derivatives       farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinG_W(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    // prepare assembly
    std::map<int,double>& dgmap = cnode->CoData().GetDerivGW();
    std::map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;
      //std::cout << "Assemble S: " << row << " " << col << " " << val << std::endl;
      // do not assemble zeros into s matrix
      if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }
  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}
/*----------------------------------------------------------------------*
 |  Assemble matrix LinStick with tangential+D+M derivatives  mgit 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
                                          LINALG::SparseMatrix& linstickDISglobal,
                                          Epetra_Vector& linstickRHSglobal)
{

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
  Teuchos::RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements()==0)
    return;

  INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  bool consistent = false;

#if defined(CONSISTENTSTICK) && defined(CONSISTENTSLIP)
  dserror("It's not reasonable to activate both, the consistent stick and slip branch, "
      "because both together will lead again to an inconsistent formulation!");
#endif
#ifdef CONSISTENTSTICK
  consistent = true;
#endif

  if (consistent && ftype == INPAR::CONTACT::friction_coulomb)
  {
    // loop over all stick nodes of the interface
    for (int i=0;i<sticknodes->NumMyElements();++i)
    {
      int gid = sticknodes->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

      // prepare assembly, get information from node
      std::vector<std::map<int,double> > dnmap = cnode->CoData().GetDerivN();
      std::vector<std::map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<std::map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
      std::map<int,double> dgmap = cnode->CoData().GetDerivG();

      // check for Dimension of derivative maps
      for (int j=0;j<Dim()-1;++j)
        if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivN-map is inconsistent!");

      for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim()==3)
      {
        for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      // information from interface contact parameter list
      double frcoeff = 1.0;//IParams().get<double>("FRCOEFF");
      double ct = 1.0;//IParams().get<double>("SEMI_SMOOTH_CT");
      double cn = IParams().get<double>("SEMI_SMOOTH_CN");

      // more information from node
      double* n = cnode->MoData().n();
      double* z = cnode->MoData().lm();
      double& wgap = cnode->CoData().Getg();

      // iterator for maps
      std::map<int,double>::iterator colcurr;

      // row number of entries
      std::vector<int> row (Dim()-1);
      if (Dim()==2)
      {
        row[0] = stickt->GID(i);
      }
      else if (Dim()==3)
      {
        row[0] = stickt->GID(2*i);
        row[1] = stickt->GID(2*i)+1;
      }
      else
        dserror("ERROR: AssemblelinStick: Dimension not correct");

      // evaluation of specific components of entries to assemble
      // evaluation of specific components of entries to assemble
      double znor = 0;
      double jumptxi = 0;
      double jumpteta = 0;

      for (int i=0;i<Dim();i++)
        znor += n[i]*z[i];

#ifdef OBJECTVARSLIPINCREMENT

        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];

#else
        double* jump = cnode->FriData().jump();
        double* txi = cnode->CoData().txi();
        double* teta = cnode->CoData().teta();

        // more information from node
        for (int i=0;i<Dim();i++)
        {
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
#endif


      // check for dimensions
      if(Dim()==2 and (jumpteta != 0.0))
        dserror ("ERROR: AssembleLinStick: jumpteta must be zero in 2D");

      //**************************************************************
       // calculation of matrix entries of linearized stick condition
       //**************************************************************

       // 1) Entries from differentiation with respect to LM
       /******************************************************************/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        double valtxi = 0.0;
        double valteta = 0.0;
        int col = cnode->Dofs()[dim];
        valtxi = -frcoeff * ct * jumptxi * n[dim];

        if (Dim()==3)
        {
          valteta = -frcoeff * ct * jumpteta * n[dim];
        }
        // do not assemble zeros into matrix
        if (abs(valtxi)>1.0e-12) linstickLMglobal.Assemble(valtxi,row[0],col);
        if (Dim()==3)
        if (abs(valteta)>1.0e-12) linstickLMglobal.Assemble(valteta,row[1],col);
      }

      // Entries on right hand side ****************************
      Epetra_SerialDenseVector rhsnode(Dim()-1);
      std::vector<int> lm(Dim()-1);
      std::vector<int> lmowner(Dim()-1);
      rhsnode(0) = frcoeff * (znor - cn * wgap) * ct * jumptxi;

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();
      if (Dim()==3)
      {
        rhsnode(1) = frcoeff * (znor - cn * wgap) * ct * jumpteta;
        lm[1] = cnode->Dofs()[2];
        lmowner[1] = cnode->Owner();
      }

      LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);

      // 3) Entries from differentiation with respect to displacements
      /******************************************************************/

#ifdef OBJECTVARSLIPINCREMENT
      std::vector<std::map<int,double> > derivjump_ = cnode->FriData().GetDerivVarJump();

      //txi
      for (colcurr=derivjump_[0].begin();colcurr!=derivjump_[0].end();++colcurr)
      {
        int col = colcurr->first;
        double valtxi = - frcoeff * (znor - cn * wgap) * ct * (colcurr->second);
        if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
      }
      //teta
      for (colcurr=derivjump_[1].begin();colcurr!=derivjump_[1].end();++colcurr)
      {
        int col = colcurr->first;
        double valteta = - frcoeff * (znor - cn * wgap) * ct * (colcurr->second);
        if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
      }

      // ... old slip
      // get linearization of jump vector
     std::vector<std::map<int,double> > derivjump = cnode->FriData().GetDerivJump();

     // loop over dimensions
     for (int dim=0;dim<cnode->NumDof();++dim)
     {
#else
     // ... old slip
     // get linearization of jump vector
    std::vector<std::map<int,double> > derivjump = cnode->FriData().GetDerivJump();

    // loop over dimensions
    for (int dim=0;dim<cnode->NumDof();++dim)
    {
      // loop over all entries of the current derivative map (jump)
      for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
      {
        int col = colcurr->first;

        double valtxi=0.0;
        valtxi = - frcoeff * (znor - cn * wgap) * ct * txi[dim] * (colcurr->second);
        // do not assemble zeros into matrix
        if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

        if (Dim()==3)
        {
          double valteta=0.0;
          valteta = - frcoeff * (znor - cn * wgap) * ct * teta[dim] * (colcurr->second);
          // do not assemble zeros into matrix
          if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }
      }

      // linearization first tangential direction *********************************
      // loop over all entries of the current derivative map (txi)
      for (colcurr=dtximap[dim].begin();colcurr!=dtximap[dim].end();++colcurr)
      {
        int col = colcurr->first;
        double valtxi=0.0;
        valtxi = - frcoeff*(znor-cn*wgap) * ct * jump[dim] * colcurr->second;


        // do not assemble zeros into matrix
        if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
      }
      // linearization second tangential direction *********************************
      if (Dim()==3)
      {
        // loop over all entries of the current derivative map (teta)
        for (colcurr=dtetamap[dim].begin();colcurr!=dtetamap[dim].end();++colcurr)
        {
          int col = colcurr->first;
          double valteta=0.0;
          valteta = - frcoeff * (znor-cn*wgap) * ct * jump[dim] * colcurr->second;

          // do not assemble zeros into matrix
          if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }
      }
#endif
        // linearization of normal direction *****************************************
        // loop over all entries of the current derivative map
        for (colcurr=dnmap[dim].begin();colcurr!=dnmap[dim].end();++colcurr)
        {
          int col = colcurr->first;
          double valtxi=0.0;
          valtxi = - frcoeff * z[dim] * colcurr->second * ct * jumptxi;
          // do not assemble zeros into matrix
          if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

          if (Dim()==3)
          {
            double valteta=0.0;
            valteta = - frcoeff * z[dim] * colcurr->second * ct * jumpteta;
            // do not assemble zeros into matrix
            if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
          }
        }
    } // loop over all dimensions

    // linearization of weighted gap**********************************************
    // loop over all entries of the current derivative map fixme
    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
    {
      int col = colcurr->first;
      double valtxi=0.0;
      valtxi = frcoeff * colcurr->second * ct * cn * jumptxi;
      // do not assemble zeros into matrix
      if (abs(valtxi)>1e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
      if (Dim()==3)
      {
        double valteta=0.0;
        valteta = frcoeff * colcurr->second * ct * cn * jumpteta;
        // do not assemble zeros into matrix
        if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
      }
    }
    if(wearimpl_)
    {
      // linearization of weighted wear w.r.t. displacements **********************
      std::map<int,double>& dwmap = cnode->CoData().GetDerivW();
      for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
      {
        int col = colcurr->first;
        double valtxi=0.0;
        valtxi = frcoeff * colcurr->second * ct * cn * jumptxi;
        // do not assemble zeros into matrix
        if (abs(valtxi)>1e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
        if (Dim()==3)
        {
          double valteta=0.0;
          valteta = frcoeff * colcurr->second * ct * cn * jumpteta;
          // do not assemble zeros into matrix
          if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }
      }
    } // if wearimpl_
    } // loop over stick nodes
  }
  else
  {
    // loop over all stick nodes of the interface
    for (int i=0;i<sticknodes->NumMyElements();++i)
    {
      int gid = sticknodes->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

      // prepare assembly, get information from node
      std::vector<std::map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<std::map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();

      for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
        dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim()==3)
      {
        for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      // iterator for maps
      std::map<int,double>::iterator colcurr;

      // row number of entries
      std::vector<int> row (Dim()-1);
      if (Dim()==2)
      {
        row[0] = stickt->GID(i);
      }
      else if (Dim()==3)
      {
        row[0] = stickt->GID(2*i);
        row[1] = stickt->GID(2*i)+1;
      }
      else
        dserror("ERROR: AssemblelinStick: Dimension not correct");

      // evaluation of specific components of entries to assemble
      double jumptxi=0;
      double jumpteta=0;

#ifdef OBJECTVARSLIPINCREMENT

      jumptxi=cnode->FriData().jump_var()[0];

      if (Dim()==3)
        jumpteta=cnode->FriData().jump_var()[1];

#else
      // more information from node
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* jump = cnode->FriData().jump();

      for (int i=0;i<Dim();i++)
      {
        jumptxi += txi[i]*jump[i];
        jumpteta += teta[i]*jump[i];
      }
#endif
      // check for dimensions
      if(Dim()==2 and (jumpteta != 0.0))
        dserror ("ERROR: AssembleLinStick: jumpteta must be zero in 2D");

      // Entries on right hand side
      /************************************************ (-utxi, -uteta) ***/
      Epetra_SerialDenseVector rhsnode(Dim()-1);
      std::vector<int> lm(Dim()-1);
      std::vector<int> lmowner(Dim()-1);

      // modification to stabilize the convergence of the lagrange multiplier incr (hiermeier 08/13)
      if (abs(jumptxi)<1e-15)
        rhsnode(0) = 0.0;
      else
        rhsnode(0) = -jumptxi;

      lm[0] = cnode->Dofs()[1];
      lmowner[0] = cnode->Owner();

      if (Dim()==3)
      {
        // modification to stabilize the convergence of the lagrange multiplier incr (hiermeier 08/13)
        if (abs(jumpteta)<1e-15)
          rhsnode(1) = 0.0;
        else
          rhsnode(1) = -jumpteta;

        lm[1] = cnode->Dofs()[2];
        lmowner[1] = cnode->Owner();
      }

      LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);

      // Entries from differentiation with respect to displacements
      /*** 1 ************************************** tangent.deriv(jump) ***/
#ifdef OBJECTVARSLIPINCREMENT
      std::map<int,double> derivjump1 = cnode->FriData().GetDerivVarJump()[0];
      std::map<int,double> derivjump2 = cnode->FriData().GetDerivVarJump()[1];

      for (colcurr=derivjump1.begin();colcurr!=derivjump1.end();++colcurr)
      {
        int col = colcurr->first;
        double valtxi = colcurr->second;

        if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
      }

      if(Dim()==3)
      {
        for (colcurr=derivjump2.begin();colcurr!=derivjump2.end();++colcurr)
        {
          int col = colcurr->first;
          double valteta = colcurr->second;

          if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }
      }

#else
      // get linearization of jump vector
      std::vector<std::map<int,double> > derivjump = cnode->FriData().GetDerivJump();

      if (derivjump.size()<1)
        dserror ("AssembleLinStick: Derivative of jump is not exiting!");

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map (jump)
        for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
        {
          int col = colcurr->first;
          double valtxi = txi[dim]*colcurr->second;

          // do not assemble zeros into matrix
          if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

          if(Dim()==3)
          {
            double valteta = teta[dim]*colcurr->second;
            if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
          }
        }
      }

      /*** 2 ************************************** deriv(tangent).jump ***/
      // loop over dimensions
      for (int j=0;j<Dim();++j)
      {
        // loop over all entries of the current derivative map (txi)
        for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = jump[j]*colcurr->second;

          // do not assemble zeros into s matrix
          if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
        }

        if(Dim()==3)
        {
          // loop over all entries of the current derivative map (teta)
          for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = jump[j]*colcurr->second;

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
          }
        }
      }
#endif
    }
  }
  return;
}

  /*----------------------------------------------------------------------*
  |  Assemble matrix LinSlip with W derivatives               farah 09/13|
  *----------------------------------------------------------------------*/
  void CONTACT::WearInterface::AssembleLinSlip_W(LINALG::SparseMatrix& linslipWglobal)
  {
    // get out of here if not participating in interface
    if (!lComm())
      return;

    // nothing to do if no slip nodes
    if (slipnodes_->NumMyElements()==0)
      return;

    // information from interface contact parameter list
    INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
    double frcoeff = IParams().get<double>("FRCOEFF");
    double ct = IParams().get<double>("SEMI_SMOOTH_CT");
    double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  #ifdef CONTACTFRICTIONLESSFIRST
    // get systemtype if CONTACTFRICTIONLESSFIRST is active
    INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(IParams(),"SYSTEM");
  #endif

    //**********************************************************************
    //**********************************************************************
    //**********************************************************************
    // Coulomb Friction
    //**********************************************************************
    //**********************************************************************
    //**********************************************************************
    if (ftype == INPAR::CONTACT::friction_coulomb)
    {
      // loop over all slip nodes of the interface
      for (int i=0;i<slipnodes_->NumMyElements();++i)
      {
        int gid = slipnodes_->GID(i);
        DRT::Node* node = idiscret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        FriNode* cnode = static_cast<FriNode*>(node);

        if (cnode->Owner() != Comm().MyPID())
          dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

        // prepare assembly, get information from node
        std::vector<std::map<int,double> > dnmap = cnode->CoData().GetDerivN();
        std::vector<std::map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
        std::vector<std::map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
        double scalefac=1.;
        std::map<int,double> dscmap = cnode->CoData().GetDerivScale();

        // check for Dimension of derivative maps
        for (int j=0;j<Dim()-1;++j)
          if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

         for (int j=0;j<Dim()-1;++j)
            if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
              dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

         if (Dim()==3)
         {
           for (int j=0;j<Dim()-1;++j)
            if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
              dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
         }

        // more information from node
        double* n = cnode->MoData().n();
        double* txi = cnode->CoData().txi();
        double* teta = cnode->CoData().teta();
        double* z = cnode->MoData().lm();
        double& wgap = cnode->CoData().Getg();
        wgap /= scalefac;

        // iterator for maps
        std::map<int,double>::iterator colcurr;

        // row number of entries
        std::vector<int> row (Dim()-1);
        if (Dim()==2)
        {
          row[0] = slipt_->GID(i);
        }
        else if (Dim()==3)
        {
          row[0] = slipt_->GID(2*i);
          row[1] = slipt_->GID(2*i)+1;
        }
        else
          dserror("ERROR: AssemblelinSlip: Dimension not correct");

        // boolean variable if flag "CONTACTFRICTIONLESSFIRST" AND
        // ActiveOld = true
        bool friclessandfirst = false;

        // evaluation of specific components of entries to assemble
        double znor = 0;
        double ztxi = 0;
        double zteta = 0;
        double jumptxi = 0;
        double jumpteta = 0;
        double euclidean = 0;

  #ifdef OBJECTVARSLIPINCREMENT

        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];

        for (int i=0;i<Dim();i++)
        {
          znor += n[i]*z[i];
          ztxi += txi[i]*z[i];
          zteta += teta[i]*z[i];
        }
  #else
        double* jump = cnode->FriData().jump();
        for (int i=0;i<Dim();i++)
        {
          znor += n[i]*z[i];
          ztxi += txi[i]*z[i];
          zteta += teta[i]*z[i];
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
  #endif

        // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
        std::vector<double> sum1 (Dim()-1,0);
        sum1[0] = ztxi+ct*jumptxi;
        if (Dim()==3) sum1[1] = zteta+ct*jumpteta;
        if (Dim()==2) euclidean = abs(sum1[0]);
        if (Dim()==3) euclidean = sqrt(sum1[0]*sum1[0]+sum1[1]*sum1[1]);

        // check of dimensions
        if(Dim()==2 and (zteta != 0.0 or jumpteta != 0.0))
          dserror ("ERROR: AssemblelinSlip: zteta and jumpteta must be zero in 2D");

        // check of euclidean norm
        if (euclidean==0.0)
          dserror ("ERROR: AssemblelinSlip: Euclidean norm is zero");

  #ifdef CONTACTFRICTIONLESSFIRST

        // in the case of frictionless contact for nodes just coming into
        // contact, the frictionless contact condition is applied.
        if (cnode->FriData().ActiveOld()==false)
        {
          friclessandfirst=true;
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cnode->Dofs()[dim];
            double valtxi = txi[dim];
            double valteta = 0;
            if (Dim()==3) valteta = teta[dim];

            if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
            if (Dim()==3)
              if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);

          }

          Epetra_SerialDenseVector rhsnode(Dim()-1);
          std::vector<int> lm(Dim() - 1);
          std::vector<int> lmowner(Dim() - 1);

          rhsnode(0)  = 0.0;
          lm[0]     = cnode->Dofs()[1];
          lmowner[0]  = cnode->Owner();


          if (systype == INPAR::CONTACT::system_condensed)
            rhsnode[0] = 0.0;
          else
            rhsnode[0] = -ztxi;   // already negative rhs!!!

          if(Dim()==3)
          {
            if (systype == INPAR::CONTACT::system_condensed)
              rhsnode[1] = 0.0;
            else
                rhsnode[1] = -zteta;    // already negative rhs!!!

              lm[1] = cnode->Dofs()[2];
              lmowner[1] = cnode->Owner();
          }
          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);


          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            for (colcurr=dtximap[dim].begin();colcurr!=dtximap[dim].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (colcurr->second)*z[dim];
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            }

            if(Dim()==3)
            {
              for (colcurr=dtetamap[dim].begin();colcurr!=dtetamap[dim].end();++colcurr)
              {
                int col = colcurr->first;
                double valteta = (colcurr->second)*z[dim];
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }
        }
  #endif

        // this is not evaluated if "FRICTIONLESSFIRST" is flaged on AND the node
        // is just coming into contact
        if(friclessandfirst==false)
        {
          //****************************************************************
          // CONSISTENT TREATMENT OF CASE FRCOEFF=0 (FRICTIONLESS)
          //****************************************************************
          // popp 08/2012
          //
          // There is a problem with our frictional nonlinear complementarity
          // function when applied to the limit case frcoeff=0 (frictionless).
          // In this case, the simple frictionless sliding condition should
          // be consistently recovered, which unfortunately is not the case.
          // This fact is well-known (see PhD thesis S. Hüeber) and now
          // taken care of by a special treatment as can be seen below
          //
          //****************************************************************
          if (frcoeff==0.0)
          {
            //coming soon....
          }

          //****************************************************************
          // STANDARD TREATMENT OF CASE FRCOEFF!=0 (FRICTIONAL)
          //****************************************************************
          else
          {
            //**************************************************************
            // calculation of matrix entries of linearized slip condition
            //**************************************************************

            /*** 1 ****************** frcoeff*cn*deriv (g).(ztan+ct*utan) ***/
            // prepare assembly
            std::map<int,double>& dgwmap = cnode->CoData().GetDerivGW();

            // loop over all entries of the current derivative map
            for (colcurr=dgwmap.begin();colcurr!=dgwmap.end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = frcoeff*cn*(colcurr->second)*(ztxi+ct*jumptxi);
              double valteta = frcoeff*cn*(colcurr->second)*(zteta+ct*jumpteta);

              // do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipWglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipWglobal.Assemble(valteta,row[1],col);
            }
          } // if (frcoeff==0.0)
        } // if (frictionlessandfirst == false)
      } // loop over all slip nodes of the interface
    } // Coulomb friction
    else
      dserror("linslip wear only for coulomb friction!");

    return;
  }

/*----------------------------------------------------------------------*
|  Assemble matrix LinSlip with tangential+D+M derivatives    mgit 02/09|
*----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
                                         LINALG::SparseMatrix& linslipDISglobal,
                                         Epetra_Vector& linslipRHSglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements()==0)
    return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
  double frbound = IParams().get<double>("FRBOUND");
  double frcoeff = IParams().get<double>("FRCOEFF");
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
#ifdef CONTACTFRICTIONLESSFIRST
  // get systemtype if CONTACTFRICTIONLESSFIRST is active
  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(IParams(),"SYSTEM");
#endif

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  // Coulomb Friction
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  if (ftype == INPAR::CONTACT::friction_coulomb)
  {
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // prepare assembly, get information from node
      std::vector<std::map<int,double> > dnmap = cnode->CoData().GetDerivN();
      std::vector<std::map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<std::map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
      double scalefac=1.;
      std::map<int,double> dscmap = cnode->CoData().GetDerivScale();
      bool scderiv=false;
      if (DRT::INPUT::IntegralValue<int>(imortar_,"LM_NODAL_SCALE")==true &&
                cnode->MoData().GetScale() != 0.)
      {
        scderiv=true;
        scalefac=cnode->MoData().GetScale();

      }

      // check for Dimension of derivative maps
      for (int j=0;j<Dim()-1;++j)
        if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

       for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

       if (Dim()==3)
       {
         for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
       }

      // more information from node
      double* n = cnode->MoData().n();
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* z = cnode->MoData().lm();
      double& wgap = cnode->CoData().Getg();
      wgap /= scalefac;

      // iterator for maps
      std::map<int,double>::iterator colcurr;

      // row number of entries
      std::vector<int> row (Dim()-1);
      if (Dim()==2)
      {
        row[0] = slipt_->GID(i);
      }
      else if (Dim()==3)
      {
        row[0] = slipt_->GID(2*i);
        row[1] = slipt_->GID(2*i)+1;
      }
      else
        dserror("ERROR: AssemblelinSlip: Dimension not correct");

      // boolean variable if flag "CONTACTFRICTIONLESSFIRST" AND
      // ActiveOld = true
      bool friclessandfirst = false;

      // evaluation of specific components of entries to assemble
      double znor = 0;
      double ztxi = 0;
      double zteta = 0;
      double jumptxi = 0;
      double jumpteta = 0;
      double euclidean = 0;

#ifdef OBJECTVARSLIPINCREMENT

      jumptxi=cnode->FriData().jump_var()[0];

      if (Dim()==3)
        jumpteta=cnode->FriData().jump_var()[1];

      for (int i=0;i<Dim();i++)
      {
        znor += n[i]*z[i];
        ztxi += txi[i]*z[i];
        zteta += teta[i]*z[i];
      }
#else
      double* jump = cnode->FriData().jump();
      for (int i=0;i<Dim();i++)
      {
        znor += n[i]*z[i];
        ztxi += txi[i]*z[i];
        zteta += teta[i]*z[i];
        jumptxi += txi[i]*jump[i];
        jumpteta += teta[i]*jump[i];
      }
#endif

      // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
      std::vector<double> sum1 (Dim()-1,0);
      sum1[0] = ztxi+ct*jumptxi;
      if (Dim()==3) sum1[1] = zteta+ct*jumpteta;
      if (Dim()==2) euclidean = abs(sum1[0]);
      if (Dim()==3) euclidean = sqrt(sum1[0]*sum1[0]+sum1[1]*sum1[1]);

      // check of dimensions
      if(Dim()==2 and (zteta != 0.0 or jumpteta != 0.0))
        dserror ("ERROR: AssemblelinSlip: zteta and jumpteta must be zero in 2D");

      // check of euclidean norm
      if (euclidean==0.0)
      {
        std::cout << "owner= " << cnode->Owner() <<  "  " << lComm()->MyPID() <<std::endl;
        dserror ("ERROR: AssemblelinSlip: Euclidean norm is zero");
      }

#ifdef CONTACTFRICTIONLESSFIRST

      // in the case of frictionless contact for nodes just coming into
      // contact, the frictionless contact condition is applied.
      if (cnode->FriData().ActiveOld()==false)
      {
        friclessandfirst=true;
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];
          double valteta = 0;
          if (Dim()==3) valteta = teta[dim];

          if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);

        }

        Epetra_SerialDenseVector rhsnode(Dim()-1);
        std::vector<int> lm(Dim() - 1);
        std::vector<int> lmowner(Dim() - 1);

        rhsnode(0)  = 0.0;
        lm[0]     = cnode->Dofs()[1];
        lmowner[0]  = cnode->Owner();


        if (systype == INPAR::CONTACT::system_condensed)
          rhsnode[0] = 0.0;
        else
          rhsnode[0] = -ztxi;   // already negative rhs!!!

        if(Dim()==3)
        {
          if (systype == INPAR::CONTACT::system_condensed)
            rhsnode[1] = 0.0;
          else
              rhsnode[1] = -zteta;    // already negative rhs!!!

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
        }
        LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);


        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          for (colcurr=dtximap[dim].begin();colcurr!=dtximap[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (colcurr->second)*z[dim];
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
          }

          if(Dim()==3)
          {
            for (colcurr=dtetamap[dim].begin();colcurr!=dtetamap[dim].end();++colcurr)
            {
              int col = colcurr->first;
              double valteta = (colcurr->second)*z[dim];
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }
        }
      }
#endif

      // this is not evaluated if "FRICTIONLESSFIRST" is flaged on AND the node
      // is just coming into contact
      if(friclessandfirst==false)
      {
        //****************************************************************
        // CONSISTENT TREATMENT OF CASE FRCOEFF=0 (FRICTIONLESS)
        //****************************************************************
        // popp 08/2012
        //
        // There is a problem with our frictional nonlinear complementarity
        // function when applied to the limit case frcoeff=0 (frictionless).
        // In this case, the simple frictionless sliding condition should
        // be consistently recovered, which unfortunately is not the case.
        // This fact is well-known (see PhD thesis S. Hüeber) and now
        // taken care of by a special treatment as can be seen below
        //
        //****************************************************************
        if (frcoeff==0.0)
        {
          //**************************************************************
          // calculation of matrix entries of linearized slip condition
          //**************************************************************

          // 1) Entries from differentiation with respect to LM
          /******************************************************************/

          // loop over the dimension
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cnode->Dofs()[dim];
            double valtxi = txi[dim];

            double valteta = 0.0;
            if (Dim()==3) valteta = teta[dim];

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
            if (Dim()==3)
              if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);
          }

          // 2) Entries on right hand side
          /******************************************************************/

          Epetra_SerialDenseVector rhsnode(Dim()-1);
          std::vector<int> lm(Dim()-1);
          std::vector<int> lmowner(Dim()-1);

          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          rhsnode[0] = -ztxi;   // already negative rhs!!!

          if(Dim()==3)
          {
            rhsnode[1] = -zteta;    // already negative rhs!!!

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }

          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

          // 3) Entries from differentiation with respect to displacements
          /******************************************************************/

          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = (colcurr->second)*z[j];

              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if (Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = (colcurr->second)*z[j];

                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }
        }

        //****************************************************************
        // STANDARD TREATMENT OF CASE FRCOEFF!=0 (FRICTIONAL)
        //****************************************************************
        else
        {
          //**************************************************************
          // calculation of matrix entries of linearized slip condition
          //**************************************************************

          // 1) Entries from differentiation with respect to LM
          /******************************************************************/

          // loop over the dimension
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            double valtxi = 0.0;
            int col = cnode->Dofs()[dim];

            double valtxi0 = euclidean*txi[dim];
            double valtxi1 = ((ztxi+ct*jumptxi)/euclidean*ztxi)*txi[dim];
            double valtxi3 = (zteta+ct*jumpteta)/euclidean*ztxi*teta[dim];

#ifdef CONSISTENTSLIP
            valtxi0 = valtxi0 / (znor - cn * wgap);
            valtxi1 = valtxi1 / (znor - cn * wgap);
            valtxi3 = valtxi3 / (znor - cn * wgap);

            // Additional term
            valtxi0 -= euclidean * ztxi / pow(znor -cn * wgap,2.0) * n[dim];

            double valtxi2 = -frcoeff * txi[dim];
#else
            double valtxi2 = -frcoeff*(znor-cn*wgap)*txi[dim]-frcoeff*(ztxi+ct*jumptxi)*n[dim];
#endif
            valtxi = valtxi0 + valtxi1 + valtxi2 + valtxi3;

            double valteta = 0.0;
            if (Dim()==3)
            {
              double valteta0 = euclidean*teta[dim];
              double valteta1 = ((ztxi+ct*jumptxi)/euclidean*zteta)*txi[dim];
              double valteta3 = (zteta+ct*jumpteta)/euclidean*zteta*teta[dim];
#ifdef CONSISTENTSLIP
              valteta0 = valteta0 / (znor - cn * wgap);
              valteta1 = valteta1 / (znor - cn * wgap);
              valteta3 = valteta3 / (znor - cn * wgap);

              // Additional term
              valteta0 -= euclidean * zteta / pow(znor -cn * wgap,2.0) * n[dim];

              double valteta2 = -frcoeff * teta[dim];
#else
              double valteta2 = -frcoeff*(znor-cn*wgap)*teta[dim]-frcoeff*(zteta+ct*jumpteta)*n[dim];
#endif
              valteta = valteta0 + valteta1 + valteta2 + valteta3;
            }

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
            if (Dim()==3)
              if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);
          }

          // 2) Entries on right hand side
          /******************************************************************/
          Epetra_SerialDenseVector rhsnode(Dim()-1);
          std::vector<int> lm(Dim()-1);
          std::vector<int> lmowner(Dim()-1);

#ifdef CONSISTENTSLIP
          double valuetxi1 = -(euclidean)*ztxi / (znor - cn * wgap) + frcoeff*(ztxi+ct*jumptxi);
#else
          double valuetxi1 = -(euclidean)*ztxi+(frcoeff*(znor-cn*wgap))*(ztxi+ct*jumptxi);
#endif
          rhsnode(0) = valuetxi1;
          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          if(Dim()==3)
          {
#ifdef CONSISTENTSLIP
          double valueteta1 = -(euclidean)*zteta / (znor - cn * wgap) + frcoeff * (zteta + ct * jumpteta);
#else
          double valueteta1 = -(euclidean)*zteta+(frcoeff*(znor-cn*wgap))*(zteta+ct*jumpteta);
#endif

          rhsnode(1) = valueteta1;

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }

          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

          // 3) Entries from differentiation with respect to displacements
          /******************************************************************/

          /*** 01  ********* -Deriv(euclidean).ct.tangent.deriv(u)*ztan ***/
#ifdef OBJECTVARSLIPINCREMENT
          std::map<int,double> derivjump1 = cnode->FriData().GetDerivVarJump()[0];
          std::map<int,double> derivjump2 = cnode->FriData().GetDerivVarJump()[1];

          for (colcurr=derivjump1.begin();colcurr!=derivjump1.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi1 = (ztxi+ct*jumptxi)/euclidean*ct*colcurr->second*ztxi;
            double valteta1 = (ztxi+ct*jumptxi)/euclidean*ct*colcurr->second*zteta;

            if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
            if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);

          }

          if (Dim()==3)
          {
            for (colcurr=derivjump2.begin();colcurr!=derivjump2.end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi2 = (zteta+ct*jumpteta)/euclidean*ct*colcurr->second*ztxi;
              double valteta2 = (zteta+ct*jumpteta)/euclidean*ct*colcurr->second*zteta;

              if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
              if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
            }
          }

#else
          // get linearization of jump vector
          std::vector<std::map<int,double> > derivjump = cnode->FriData().GetDerivJump();

          // loop over dimensions
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            // loop over all entries of the current derivative map (jump)
            for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
            {
              int col = colcurr->first;

              double valtxi1 = (ztxi+ct*jumptxi)/euclidean*ct*txi[dim]*colcurr->second*ztxi;
              double valteta1 = (ztxi+ct*jumptxi)/euclidean*ct*txi[dim]*colcurr->second*zteta;
              double valtxi2 = (zteta+ct*jumpteta)/euclidean*ct*teta[dim]*colcurr->second*ztxi;
              double valteta2 = (zteta+ct*jumpteta)/euclidean*ct*teta[dim]*colcurr->second*zteta;

#ifdef CONSISTENTSLIP
              valtxi1   = valtxi1  / (znor - cn * wgap);
              valteta1  = valteta1 / (znor - cn * wgap);
              valtxi2   = valtxi2  / (znor - cn * wgap);
              valteta2  = valteta2 / (znor - cn * wgap);
#endif

              // do not assemble zeros into matrix
              if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
              if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
              if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
              if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
            }

#ifdef CONSISTENTSLIP
            /*** Additional Terms ***/
            // normal derivative
            for (colcurr=dnmap[dim].begin();colcurr!=dnmap[dim].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi   = - euclidean * ztxi * z[dim] / pow(znor - cn * wgap, 2.0) * (colcurr->second);
              double valteta  = - euclidean * zteta * z[dim] / pow(znor - cn * wgap, 2.0) * (colcurr->second);

              // do not assemble zeros into s matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
#endif
          }
#endif


#ifdef CONSISTENTSLIP
          /*** Additional Terms ***/
          // wgap derivative
          std::map<int,double>& dgmap = cnode->CoData().GetDerivG();

          for (colcurr=dgmap.begin(); colcurr!=dgmap.end(); ++colcurr)
          {
            int col = colcurr->first;
            double valtxi  = + euclidean * ztxi  / pow(znor - cn * wgap, 2.0) * cn * (colcurr->second)/scalefac;
            double valteta = + euclidean * zteta / pow(znor - cn * wgap, 2.0) * cn * (colcurr->second)/scalefac;

            //do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
<<<<<<< .mine
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
=======
          if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
>>>>>>> .r18285
          }
          if (scderiv)
            for (colcurr=dscmap.begin(); colcurr!=dscmap.end(); ++colcurr)
            {
              int col =colcurr->first;
              double valtxi  = + euclidean * ztxi  / pow(znor - cn * wgap, 2.0) * cn * wgap/scalefac*colcurr->second;
              double valteta = + euclidean * zteta / pow(znor - cn * wgap, 2.0) * cn * wgap/scalefac*colcurr->second;
              //do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
#endif


          /*** 02 ***************** frcoeff*znor*ct*tangent.deriv(jump) ***/
#ifdef OBJECTVARSLIPINCREMENT

          for (colcurr=derivjump1.begin();colcurr!=derivjump1.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = -frcoeff*(znor-cn*wgap)*ct*colcurr->second;

            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
          }

          if (Dim()==3)
          {
            for (colcurr=derivjump2.begin();colcurr!=derivjump2.end();++colcurr)
            {
              int col = colcurr->first;
              double valteta = -frcoeff*(znor-cn*wgap)*ct*colcurr->second;

              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }

#else
          // loop over dimensions
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            // loop over all entries of the current derivative map (jump)
            for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
           {
              int col = colcurr->first;

              //std::cout << "val " << colcurr->second << std::endl;
#ifdef CONSISTENTSLIP
              double valtxi = - frcoeff * ct * txi[dim] * colcurr->second;
              double valteta = - frcoeff * ct * teta[dim] * colcurr->second;
#else
              double valtxi = (-1)*(frcoeff*(znor-cn*wgap))*ct*txi[dim]*colcurr->second;
              double valteta = (-1)*(frcoeff*(znor-cn*wgap))*ct*teta[dim]*colcurr->second;
#endif
              // do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);

              if (Dim()==3)
              {
               if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }
#endif
          /*** 1 ********************************* euclidean.deriv(T).z ***/
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = euclidean*(colcurr->second)*z[j];

#ifdef CONSISTENTSLIP
              val = val / (znor - cn * wgap);
#endif

              // do not assemble zeros into s matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if (Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double val = euclidean*(colcurr->second)*z[j];

#ifdef CONSISTENTSLIP
              val = val / (znor - cn * wgap);
#endif

                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 2 ********************* deriv(euclidean).deriv(T).z.ztan ***/
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*zteta;

#ifdef CONSISTENTSLIP
              valtxi  = valtxi / (znor - cn*wgap);
              valteta   = valteta / (znor - cn*wgap);
#endif

             // do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (Dim()==3)
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }

            if(Dim()==3)
            {
              // 3D loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*zteta;

#ifdef CONSISTENTSLIP
                valtxi  = valtxi / (znor - cn*wgap);
                valteta = valteta / (znor - cn*wgap);
#endif

                // do not assemble zeros into matrix
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }

          /*** 3 ****************** deriv(euclidean).deriv(T).jump.ztan ***/
#ifdef OBJECTVARSLIPINCREMENT
          //!!!!!!!!!!!!!!! DO NOTHING !!!!!!!
#else
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

#ifdef CONSISTENTSLIP
              valtxi  = valtxi / (znor - cn*wgap);
              valteta   = valteta / (znor - cn*wgap);
#endif

              // do not assemble zeros into s matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

#ifdef CONSISTENTSLIP
                valtxi  = valtxi / (znor - cn*wgap);
                valteta   = valteta / (znor - cn*wgap);
#endif

                // do not assemble zeros into matrix
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }
#endif
          /*** 4 ************************** (frcoeff*znor).deriv(T).z ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
#ifdef CONSISTENTSLIP
              double val = - frcoeff * (colcurr->second)*z[j];
#else
              double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];
#endif
              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
#ifdef CONSISTENTSLIP
                double val = - frcoeff * (colcurr->second) * z[j];
#else
                double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];
#endif
                // do not assemble zeros into matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }

          /*** 5 *********************** (frcoeff*znor).deriv(T).jump ***/
#ifdef OBJECTVARSLIPINCREMENT
          //!!!!!!!!!!!!!!! DO NOTHING
#else
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
#ifdef CONSISTENTSLIP
              double val = - frcoeff * ct * (colcurr->second) * jump[j];
#else
              double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];
#endif
              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
#ifdef CONSISTENTSLIP
                double val = - frcoeff * ct * (colcurr->second) * jump[j];
#else
                double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];
#endif

                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }
#endif

#ifndef CONSISTENTSLIP
          /*** 6 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map
            for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (-1)*(ztxi+ct*jumptxi)*frcoeff*(colcurr->second)*z[j];
              double valteta = (-1)*(zteta+ct*jumpteta)*frcoeff*(colcurr->second)*z[j];

              // do not assemble zeros into s matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }

          /*** 7 ****************** frcoeff*cn*deriv (g).(ztan+ct*utan) ***/
          // prepare assembly
          std::map<int,double>& dgmap = cnode->CoData().GetDerivG();

          // loop over all entries of the current derivative map
          for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = frcoeff*cn*(colcurr->second)*(ztxi+ct*jumptxi);
            double valteta = frcoeff*cn*(colcurr->second)*(zteta+ct*jumpteta);

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
          }

          /*** 8 ****************** scale factor ***/
          // loop over all entries of the current derivative map
          if (scderiv)
            for (colcurr=dscmap.begin();colcurr!=dscmap.end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = frcoeff*cn*wgap/scalefac*(ztxi+ct*jumptxi)*colcurr->second;
              double valteta = frcoeff*cn*wgap/scalefac*(zteta+ct*jumpteta)*colcurr->second;

              // do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
#endif

          /*************************************************************************
           * Wear implicit linearization w.r.t. displ.                                          *
           *************************************************************************/
          if(wearimpl_)
          {
            std::map<int,double>& dwmap = cnode->CoData().GetDerivW();
#ifdef CONSISTENTSLIP
          // loop over all entries of the current derivative map
          for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = euclidean*ztxi*(pow(znor-cn*wgap,-2.0) * cn*(colcurr->second));
            double valteta = euclidean*zteta*(pow(znor-cn*wgap,-2.0) * cn*(colcurr->second));

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
          }
#else
            // loop over all entries of the current derivative map
            for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = frcoeff*cn*(colcurr->second)*(ztxi+ct*jumptxi);
              double valteta = frcoeff*cn*(colcurr->second)*(zteta+ct*jumpteta);

              // do not assemble zeros into matrix
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
#endif // end if no consistentslip
          }// end wearimplicit
        } // if (frcoeff==0.0)
      } // if (frictionlessandfirst == false)
    } // loop over all slip nodes of the interface
  } // Coulomb friction

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  // Tresca Friction
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  if (ftype == INPAR::CONTACT::friction_tresca)
  {
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = static_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      std::vector<std::map<int,double> > dnmap = cnode->CoData().GetDerivN();

      // iterator
      std::map<int,double>::iterator colcurr;

      std::vector <std::map<int,double> > dtmap(Dim());

      for (colcurr=dnmap[0].begin(); colcurr!=dnmap[0].end(); colcurr++)
        dtmap[1].insert(std::pair<int,double>(colcurr->first,colcurr->second));

      for (colcurr=dnmap[1].begin(); colcurr!=dnmap[1].end(); colcurr++)
        dtmap[0].insert(std::pair<int,double>(colcurr->first,(-1)*colcurr->second));

      // get more information from node
      double* jump = cnode->FriData().jump();
      double* txi = cnode->CoData().txi();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];
      //double temp = ztan + ct*jumptan;
      //double epk = frbound/abs(temp);
      //double Fpk = ztan*temp/(frbound*abs(temp));
      //double Mpk = epk*(1-Fpk);
      //double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      //****************************************************************
      // CONSISTENT TREATMENT OF CASE FRBOUND=0 (FRICTIONLESS)
      //****************************************************************
      // popp 08/2012
      //
      // There is a problem with our frictional nonlinear complementarity
      // function when applied to the limit case frbound=0 (frictionless).
      // In this case, the simple frictionless sliding condition should
      // be consistently recovered, which unfortunately is not the case.
      // This fact is well-known (see PhD thesis S. Hüeber) and now
      // taken care of by a special treatment as can be seen below
      //
      //****************************************************************
      if (frbound==0.0)
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];

          // do not assemble zeros into matrix
          if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row,col);
        }

        // 2) Entries on right hand side
        /******************************************************************/

        Epetra_SerialDenseVector rhsnode(1);
        std::vector<int> lm(1);
        std::vector<int> lmowner(1);

        rhsnode(0) = -ztan;
        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/

        // loop over dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (colcurr->second)*z[j];

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }
      }

      //****************************************************************
      // STANDARD TREATMENT OF CASE FRBOUND!=0 (FRICTIONAL)
      //****************************************************************
      else
      {
        //****************************************************************
        // calculation of matrix entries of the linearized slip condition
        //****************************************************************

        // 1) Entries from differentiation with respect to LM
        /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frbound).tan ***/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = (prefactor*ztan+sum-frbound)*txi[dim];

  #ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = txi[dim];
  #endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
        }

        // 2) Entries on right hand side
        /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/
        Epetra_SerialDenseVector rhsnode(1);

        // -C and remaining terms
        double value1= -(abs(ztan+ct*jumptan))*ztan+frbound*(ztan+ct*jumptan);

        rhsnode(0) = value1;

        std::vector<int> lm(1);
        std::vector<int> lmowner(1);

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false)
        {
          lm[0]     = cnode->Dofs()[1];
          lmowner[0]  = cnode->Owner();


          if (systype == INPAR::CONTACT::system_condensed)
            rhsnode[0] = 0.0;
          else
            rhsnode[0] = -ztan;   // already negative rhs!!!
        }
#endif

        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

        // 3) Entries from differentiation with respect to displacements
        /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

        // we need the nodal entries of the D-matrix and the old one
        double D= (cnode->MoData().GetD()[0])[cnode->Dofs()[0]];
        double Dold= (cnode->FriData().GetDOld()[0])[cnode->Dofs()[0]];

        if (abs(Dold)<0.0001)
          dserror ("Error:No entry for Dold");

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
         //std::cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;
  #ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = 0;
  #endif

         // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

        // we need the nodal entries of the M-matrix and the old one
        std::vector<std::map<int,double> > mmap = cnode->MoData().GetM();
        std::vector<std::map<int,double> > mmapold = cnode->FriData().GetMOld();

        // create a set of nodes including nodes according to M entries
        // from current and previous time step
        std::set <int> mnodes;

        // iterator
        std::set<int>::iterator mcurr;

        std::set <int> mnodescurrent = cnode->FriData().GetMNodes();
        std::set <int> mnodesold = cnode->FriData().GetMNodesOld();

        for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
          mnodes.insert(*mcurr);

        for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
          mnodes.insert(*mcurr);

        // loop over all master nodes (find adjacent ones to this stick node)
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          const int* mdofs = cmnode->Dofs();

          double mik = (mmap[0])[mdofs[0]];
          double mikold = (mmapold[0])[mdofs[0]];

          // compute linstick-matrix entry of the current active node / master node pair
          // loop over all derivative maps (=dimensions)
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cmnode->Dofs()[dim];
            double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
            //std::cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

  #ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = 0;
  #endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /************************************** frbound*ct*tan.(D-Dn-1) ***/

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = frbound*ct*txi[dim]*(D-Dold);
          //std::cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

  #ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = 0;
  #endif

          // do not assemble zeros into matrix
         if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /********************************** -frbound*ct*tan.(M-Mn-1).xm ***/

        // loop over all master nodes
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          const int* mdofs = cmnode->Dofs();

          double mik = (mmap[0])[mdofs[0]];
          double mikold = (mmapold[0])[mdofs[0]];

          // loop over all derivative maps (=dimensions)
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cmnode->Dofs()[dim];
            double val = frbound*(-1)*ct*txi[dim]*(mik-mikold);
            //std::cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

  #ifdef CONTACTFRICTIONLESSFIRST
            if (cnode->FriData().ActiveOld()==false) val = 0;
  #endif
            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //std::cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
            if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //std::cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
          if (cnode->FriData().ActiveOld()==false) val = (colcurr->second)*z[j];
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //std::cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

        if (k!=colsize)
          dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        std::map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //std::cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

          // do not assemble zeros into matrix
          if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        std::map<int,std::map<int,double> >& dmmap = cnode->CoData().GetDerivM();
        std::map<int,std::map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = static_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          std::map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //std::cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /****************************************** -frbound.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*(colcurr->second)*z[j];
            //std::cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************ -frbound.ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*ct*(colcurr->second)*jump[j];
            //std::cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************* +frbound.ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
         tdotx = 0.0;
         for (int dim=0;dim<cnode->NumDof();++dim)
           tdotx += txi[dim]*xi[dim];

         // loop over all entries of the current derivative map
         for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
         {
           int col = colcurr->first;
           double val = (-1)*(-1)*frbound*ct*tdotx*colcurr->second;
           //std::cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

           // do not assemble zeros into matrix
           if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
         }

         /********************************  -frbound.ct.T.DerivM.x ******/

         // loop over all master nodes in the DerivM-map of the active slave node
         for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
         {
           int gid = dmcurr->first;
           DRT::Node* mnode = idiscret_->gNode(gid);
           if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
           FriNode* cmnode = static_cast<FriNode*>(mnode);
           double* mxi = cmnode->xspatial();

           // we need the dot product ns*xm of this node pair
           double tdotx = 0.0;
           for (int dim=0;dim<cnode->NumDof();++dim)
             tdotx += txi[dim]*mxi[dim];

           // compute entry of the current active node / master node pair
           std::map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

           // loop over all entries of the current derivative map
           for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
           {
             int col = colcurr->first;
             double val = (-1)*frbound*ct*tdotx*colcurr->second;
            //std::cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

#ifdef CONTACTFRICTIONLESSFIRST
        if (cnode->FriData().ActiveOld()==false) val = 0;
#endif

             // do not assemble zeros into matrix
             if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
           }
         }
       }
     }
   }// Tresca friction

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix W_lm containing wear w~ derivatives      farah 07/13|
 |  w.r.t. lm       -- impl wear                                                     |
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinWLm(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  if (!wearimpl_)
    dserror("This matrix deriv. is only required for implicit wear algorithm!");

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i) //(int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleWLm: Node ownership inconsistency!");

    // prepare assembly
    std::map<int,double>& dwmap = cnode->CoData().GetDerivWlm();
    std::map<int,double>::iterator colcurr;
    int row = activen_->GID(i);
    // row number of entries

    for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;

      if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix W_lmsl containing wear w~ derivatives    farah 07/13|
 |  w.r.t. lm  --> for consistent stick                                 |
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinWLmSt(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  if (!wearimpl_)
    dserror("This matrix deriv. is only required for implicit wear algorithm!");

  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
  Teuchos::RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements()==0)
    return;

  // get input params
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  double frcoeff = IParams().get<double>("FRCOEFF");


  // loop over all stick slave nodes of the interface
  for (int i=0;i<sticknodes->NumMyElements();++i)
  {
    int gid = sticknodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    // prepare assembly, get information from node
    std::vector<std::map<int,double> > dnmap = cnode->CoData().GetDerivN();
    std::map<int,double>& dwmap = cnode->CoData().GetDerivWlm();

    double* jump = cnode->FriData().jump();
    double* n = cnode->MoData().n();
    double* txi = cnode->CoData().txi();
    double* teta = cnode->CoData().teta();
    double* z = cnode->MoData().lm();

    // iterator for maps
    std::map<int,double>::iterator colcurr;

    // row number of entries
    std::vector<int> row (Dim()-1);
    if (Dim()==2)
    {
      row[0] = stickt->GID(i);
    }
    else if (Dim()==3)
    {
      row[0] = stickt->GID(2*i);
      row[1] = stickt->GID(2*i)+1;
    }
    else
      dserror("ERROR: AssemblelinSlip: Dimension not correct");

    // evaluation of specific components of entries to assemble
    double znor = 0;
    double ztxi = 0;
    double zteta = 0;
    double jumptxi = 0;
    double jumpteta = 0;
    for (int i=0;i<Dim();i++)
    {
      znor += n[i]*z[i];
      ztxi += txi[i]*z[i];
      zteta += teta[i]*z[i];
      jumptxi += txi[i]*jump[i];
      jumpteta += teta[i]*jump[i];
    }

    // loop over all entries of the current derivative map
    for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
    {
      int col = colcurr->first;
      double valtxi = frcoeff*cn*ct*jumptxi*(colcurr->second);

      // do not assemble zeros into matrix
      if (abs(valtxi)>1.0e-12) sglobal.Assemble(valtxi,row[0],col);

      if (Dim()==3)
      {
        double valteta = frcoeff*cn*ct*jumpteta*(colcurr->second);
        if (abs(valteta)>1.0e-12) sglobal.Assemble(valteta,row[1],col);
      }
    }

  }
  return;
}
/*----------------------------------------------------------------------*
 |  Assemble matrix W_lmsl containing wear w~ derivatives    farah 07/13|
 |  w.r.t. lm  --> for slip                                             |
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleLinWLmSl(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements()==0)
    return;

  if (!wearimpl_)
    dserror("This matrix deriv. is only required for implicit wear algorithm!");

  // get input params
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");

  // loop over all active slave nodes of the interface
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

    // prepare assembly, get information from node
    std::vector<std::map<int,double> > dnmap = cnode->CoData().GetDerivN();
    std::vector<std::map<int,double> > dtximap = cnode->CoData().GetDerivTxi();
    //std::vector<std::map<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
    std::map<int,double>& dwmap = cnode->CoData().GetDerivWlm();

    double* jump = cnode->FriData().jump();
    double* n = cnode->MoData().n();
    double* txi = cnode->CoData().txi();
    double* teta = cnode->CoData().teta();
    double* z = cnode->MoData().lm();

    // iterator for maps
    std::map<int,double>::iterator colcurr;

    // row number of entries
    std::vector<int> row (Dim()-1);
    if (Dim()==2)
    {
      row[0] = slipt_->GID(i);
    }
    else if (Dim()==3)
    {
      row[0] = slipt_->GID(2*i);
      row[1] = slipt_->GID(2*i)+1;
    }
    else
      dserror("ERROR: AssemblelinSlip: Dimension not correct");

    // evaluation of specific components of entries to assemble
    double znor = 0;
    double ztxi = 0;
    double zteta = 0;
    double jumptxi = 0;
    double jumpteta = 0;
    //double euclidean = 0;
    for (int i=0;i<Dim();i++)
    {
      znor += n[i]*z[i];
      ztxi += txi[i]*z[i];
      zteta += teta[i]*z[i];
      jumptxi += txi[i]*jump[i];
      jumpteta += teta[i]*jump[i];
    }

#ifdef CONSISTENTSLIP

    double euclidean = 0;
    // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
    std::vector<double> sum1 (Dim()-1,0);
    sum1[0] = ztxi+ct*jumptxi;
    if (Dim()==3) sum1[1] = zteta+ct*jumpteta;
    if (Dim()==2) euclidean = abs(sum1[0]);
    if (Dim()==3) euclidean = sqrt(sum1[0]*sum1[0]+sum1[1]*sum1[1]);

    // loop over all entries of the current derivative map
    for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
    {
      int col = colcurr->first;
      double valtxi = euclidean*ztxi*(pow(znor-cn*wgap,-2.0) * cn*(colcurr->second));

      // do not assemble zeros into matrix
      if (abs(valtxi)>1.0e-12) sglobal.Assemble(valtxi,row[0],col);

      if (Dim()==3)
      {
        double valteta = euclidean*zteta*(pow(znor-cn*wgap,-2.0) * cn*(colcurr->second));
        if (abs(valteta)>1.0e-12) sglobal.Assemble(valteta,row[1],col);
      }
    }
#else
    double frcoeff = IParams().get<double>("FRCOEFF");
    // loop over all entries of the current derivative map
    for (colcurr=dwmap.begin();colcurr!=dwmap.end();++colcurr)
    {
      int col = colcurr->first;
      double valtxi = frcoeff*cn*(colcurr->second)*(ztxi+ct*jumptxi);

      // do not assemble zeros into matrix
      if (abs(valtxi)>1.0e-12) sglobal.Assemble(valtxi,row[0],col);

      if (Dim()==3)
      {
        double valteta = frcoeff*cn*(colcurr->second)*(zteta+ct*jumpteta);
        if (abs(valteta)>1.0e-12) sglobal.Assemble(valteta,row[1],col);
      }
    }
#endif

  }
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble wear                                         gitterle 12/10|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleWear(Epetra_Vector& gglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* frinode = static_cast<FriNode*>(node);

    if (frinode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleWear: Node ownership inconsistency!");

    /**************************************************** w-vector ******/
    double wear = frinode->FriDataPlus().Wear();

    Epetra_SerialDenseVector wnode(1);
    std::vector<int> lm(1);
    std::vector<int> lmowner(1);

    wnode(0) = wear;
    lm[0] = frinode->Id();
    lmowner[0] = frinode->Owner();

    LINALG::Assemble(gglobal,wnode,lm,lmowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble Mortar matrice for both sided wear              farah 06/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleD2(LINALG::SparseMatrix& dglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(imortar_,"BOTH_SIDED_WEAR") == INPAR::CONTACT::wear_slave)
    dserror("ERROR: AssembleD2 only for mapped both-sided wear!");

  //*******************************************************
  // assemble second D matrix for both-sided wear
  //*******************************************************
  for (int i=0;i<mnoderowmap_->NumMyElements();++i)
  {
    int gid = mnoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDM: Node ownership inconsistency!");

    /**************************************************** D-matrix ******/
    if ((cnode->CoData().GetD2()).size()>0)
    {
      std::vector<std::map<int,double> > dmap = cnode->CoData().GetD2();
      int rowsize = cnode->NumDof();
      int colsize = (int)dmap[0].size();

      for (int j=0;j<rowsize-1;++j)
        if ((int)dmap[j].size() != (int)dmap[j+1].size())
          dserror("ERROR: AssembleDM: Column dim. of nodal D-map is inconsistent!");

      std::map<int,double>::iterator colcurr;

      for (int j=0;j<rowsize;++j)
      {
        int row = cnode->Dofs()[j];
        int k = 0;

        for (colcurr=dmap[j].begin();colcurr!=dmap[j].end();++colcurr)
        {
          int col = colcurr->first;
          double val = colcurr->second;

          // do the assembly into global D matrix
          if (shapefcn_ == INPAR::MORTAR::shape_dual || shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin)
          {
            // check for diagonality
            if (row!=col && abs(val)>1.0e-12)
              dserror("ERROR: AssembleDM: D-Matrix is not diagonal!");

            // create an explicitly diagonal d matrix
            if (row==col)
              dglobal.Assemble(val, row, col);
          }
          else if (shapefcn_ == INPAR::MORTAR::shape_standard)
          {
            dserror("both-sided wear only for dual lagr. mult.");
          }

          ++k;
        }

        if (k!=colsize)
          dserror("ERROR: AssembleDM: k = %i but colsize = %i",k,colsize);
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                           popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::WearInterface::BuildActiveSet(bool init)
{
  // define local variables
  std::vector<int> mynodegids(0);
  std::vector<int> mydofgids(0);
  std::vector<int> myslipnodegids(0);
  std::vector<int> myslipdofgids(0);
  std::vector<int> mymnodegids(0);
  std::vector<int> mymdofgids(0);

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // *******************************************************************
    // This is given by the CoNode member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding CoNodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (init)
    {
      // flag for initialization of init active nodes with nodal gaps
      bool initcontactbygap = DRT::INPUT::IntegralValue<int>(IParams(),"INITCONTACTBYGAP");
      // value
      double initcontactval = IParams().get<double>("INITCONTACTGAPVALUE");

      // Either init contact by definition or by gap
      if(cnode->IsInitActive() and initcontactbygap)
        dserror("Init contact either by definition in condition or by gap!");

      // check if node is initially active or, if initialization with nodal, gap,
      // the gap is smaller than the prescribed value
      if (cnode->IsInitActive() or (initcontactbygap and cnode->CoData().Getg() < initcontactval))
      {
/*
      // **********************************************************************************
      // ***                     CONSISTENT INITIALIZATION STATE                        ***
      // **********************************************************************************
      // hiermeier 08/2013
      //
      // The consistent nonlinear complementarity function C_tau is for the special initialization case
      //
      //                zn == 0   and   gap == 0
      //
      // in conjunction with the Coulomb friction model no longer unique and a special treatment is
      // necessary. For this purpose we identify the critical nodes here and set the slip-state to true.
      // In the next step the tangential part of the Lagrange multiplier vector will be set to zero. Hence,
      // we treat these nodes as frictionless nodes! (see AssembleLinSlip)
      INPAR::CONTACT::FrictionType ftype =
          DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
      if (ftype == INPAR::CONTACT::friction_coulomb)
      {
      static_cast<FriNode*>(cnode)->FriData().Slip() = true;
      static_cast<FriNode*>(cnode)->FriData().InconInit() = true;
      myslipnodegids.push_back(cnode->Id());

      for (int j=0;j<numdof;++j)
        myslipdofgids.push_back(cnode->Dofs()[j]);
      }

*/
        cnode->Active()=true;
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is initially in slip state
      if (friction_)
      {
        // do nothing: we always assume STICK at t=0
      }
    }

    // *******************************************************************
    // RE-BUILDING OF THE ACTIVE SET
    // *******************************************************************
    else
    {
      // check if node is active
      if (cnode->Active())
      {
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is in slip state
      if (friction_)
      {
        if (static_cast<FriNode*>(cnode)->FriData().Slip())
        {
          myslipnodegids.push_back(cnode->Id());

          for (int j=0;j<numdof;++j)
            myslipdofgids.push_back(cnode->Dofs()[j]);
        }
      }
    }
  }

  // loop over all master nodes - both-sided wear specific
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(imortar_,"BOTH_SIDED_WEAR") != INPAR::CONTACT::wear_slave)
  {
    for (int k=0;k<mnoderowmap_->NumMyElements();++k)
    {
      int gid = mnoderowmap_->GID(k);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);
      const int numdof = cnode->NumDof();

      if (cnode->InvolvedM()==true)
      {
        mymnodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
        {
          mymdofgids.push_back(cnode->Dofs()[j]);
        }
        //reset it
        cnode->InvolvedM()=false;
      }
    }
  }

  // create active node map and active dof map
  activenodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)mynodegids.size(),&mynodegids[0],0,Comm()));
  activedofs_  = Teuchos::rcp(new Epetra_Map(-1,(int)mydofgids.size(),&mydofgids[0],0,Comm()));

  // create map for all involved master nodes -- both-sided wear specific
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(imortar_,"BOTH_SIDED_WEAR") != INPAR::CONTACT::wear_slave)
  {
    involvednodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)mymnodegids.size(),&mymnodegids[0],0,Comm()));
    involveddofs_  = Teuchos::rcp(new Epetra_Map(-1,(int)mymdofgids.size(),&mymdofgids[0],0,Comm()));
  }

  if (friction_)
  {
    // create slip node map and slip dof map
    slipnodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)myslipnodegids.size(),&myslipnodegids[0],0,Comm()));
    slipdofs_  = Teuchos::rcp(new Epetra_Map(-1,(int)myslipdofgids.size(),&myslipdofgids[0],0,Comm()));
  }

  // split active dofs and slip dofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  finalize construction of interface (public)              mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::FillComplete(int maxdof)
{
  // store maximum global dof ID handed in
  // this ID is later needed when setting up the Lagrange multiplier
  // dof map, which of course must not overlap with existing dof ranges
  maxdofglobal_ = maxdof;

  // we'd like to call idiscret_.FillComplete(true,false,false) but this
  // will assign all nodes new degrees of freedom which we don't want.
  // We would like to use the degrees of freedom that were stored in the
  // mortar nodes. To do so, we have to create and set our own
  // version of a DofSet class before we call FillComplete on the
  // interface discretization.
  // Our special dofset class will not assign new dofs but will assign the
  // dofs stored in the nodes.
  {
    Teuchos::RCP<MORTAR::MortarDofSet> mrtrdofset = Teuchos::rcp(new MORTAR::MortarDofSet());
    Discret().ReplaceDofSet(mrtrdofset);
    // do not assign dofs yet, we'll do this below after
    // shuffling around of nodes and elements (saves time)
    Discret().FillComplete(false,false,false);
  }

  //**********************************************************************
  // check whether crosspoints / edge nodes shall be considered or not
  bool crosspoints = DRT::INPUT::IntegralValue<int>(IParams(),"CROSSPOINTS");

  // modify crosspoints / edge nodes
  if (crosspoints)
  {
    // only applicable for 2D problems up to now
    if (Dim()==3) dserror("ERROR: Crosspoint / edge node modification not yet impl. for 3D");

    // ---------------------------------------------------------------------
    // Detect relevant nodes on slave side
    // ---------------------------------------------------------------------
    // A typical application are so-called crosspoints within mortar mesh
    // tying, where this approach is necessary to avoid over-constraint.
    // Otherwise these crosspoints would be active with respect to more
    // than one interface and thus the LM cannot sufficiently represent
    // all geometrical constraints. Another typical application is mortar
    // contact, when we want to make use of symmetry boundary conditions.
    // In this case, we deliberately modify so-called edge nodes of the
    // contact boundary and thus free them from any contact constraint.
    // ---------------------------------------------------------------------
    // Basically, the status of the crosspoints / edge nodes is simply
    // changed to MASTER and consequently they will NOT carry Lagrange
    // multipliers later on. In order to sustain the partition of unity
    // property of the LM shape functions on the adjacent slave elements,
    // the LM shape functions of the adjacent nodes will be modified! This
    // way, the mortar operator entries of the crosspoints / edge nodes are
    // transfered to the neighboring slave nodes!
    // ---------------------------------------------------------------------

    for (int i=0; i<(Discret().NodeRowMap())->NumMyElements();++i)
    {
      MORTAR::MortarNode* node = static_cast<MORTAR::MortarNode*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with only 1 adjacent MortarElement
      if (node->IsSlave() && node->NumElement()==1)
      {
        //case1: linear shape functions, boundary nodes already found
        if ((node->Elements()[0])->NumNode() == 2)
        {
          node->SetBound()=true;
          node->SetSlave()=false;
        }
        //case2: quad. shape functions, middle nodes must be sorted out
        else if (node->Id() != (node->Elements()[0])->NodeIds()[2])
        {
          node->SetBound()=true;
          node->SetSlave()=false;
        }
      }
    }
  }
  //**********************************************************************

  //**********************************************************************
  // check for linear interpolation of 2D/3D quadratic Lagrange multipliers
  bool lagmultlin = (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(IParams(),"LAGMULT_QUAD")
                     == INPAR::MORTAR::lagmult_lin_lin);

  // modify crosspoints / edge nodes
  if (lagmultlin)
  {
    // modified treatment of vertex nodes and edge nodes
    // detect middle nodes (quadratic nodes) on slave side
    // set status of middle nodes -> MASTER
    // set status of vertex nodes -> SLAVE

    // loop over all elements
    for (int i=0; i<Discret().NodeRowMap()->NumMyElements(); ++i)
    {
      // get node and cast to cnode
      MORTAR::MortarNode* node = static_cast<MORTAR::MortarNode*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with shape line3 (2D), tri6 and quad8/9 (3D)
      if (node->IsSlave())
      {
        //search the first adjacent element
        MORTAR::MortarElement::DiscretizationType shape = (node->Elements()[0])->Shape();

        // which discretization type
        switch(shape)
        {
          // line3 contact elements (= quad8/9 or tri6 discretizations)
          case MORTAR::MortarElement::line3:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0]
             || node->Id() == (node->Elements()[0])->NodeIds()[1])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

          // tri6 contact elements (= tet10 discretizations)
          case MORTAR::MortarElement::tri6:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0]
             || node->Id() == (node->Elements()[0])->NodeIds()[1]
             || node->Id() == (node->Elements()[0])->NodeIds()[2])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

          // quad8 contact elements (= hex20 discretizations)
          case MORTAR::MortarElement::quad8:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0]
             || node->Id() == (node->Elements()[0])->NodeIds()[1]
             || node->Id() == (node->Elements()[0])->NodeIds()[2]
             || node->Id() == (node->Elements()[0])->NodeIds()[3])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

          // quad9 contact elements (= hex27 discretizations)
          case MORTAR::MortarElement::quad9:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0]
             || node->Id() == (node->Elements()[0])->NodeIds()[1]
             || node->Id() == (node->Elements()[0])->NodeIds()[2]
             || node->Id() == (node->Elements()[0])->NodeIds()[3])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

          // other cases
          default:
          {
            dserror("ERROR: Lin/Lin interpolation of LM only for line3/tri6/quad8/quad9 mortar elements");
            break;
          }
        } // switch(Shape)
      } // if (IsSlave())
    } // for-loop
  }
  //**********************************************************************

  // later we will export node and element column map to FULL overlap,
  // thus store the standard column maps first
  // get standard nodal column map (overlap=1)
  oldnodecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().NodeColMap())));
  // get standard element column map (overlap=1)
  oldelecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().ElementColMap())));

  // create interface local communicator
  // find all procs that have business on this interface (own or ghost nodes/elements)
  // build a Epetra_Comm that contains only those procs
  // this intra-communicator will be used to handle most stuff on this
  // interface so the interface will not block all other procs
  {
#ifdef PARALLEL
    std::vector<int> lin(Comm().NumProc());
    std::vector<int> gin(Comm().NumProc());
    for (int i=0; i<Comm().NumProc(); ++i)
      lin[i] = 0;

    // check ownership or ghosting of any elements / nodes
    //const Epetra_Map* nodemap = Discret().NodeColMap();
    //const Epetra_Map* elemap  = Discret().ElementColMap();

    //********************************************************************
    // NOTE: currently we choose local=global communicator, but we have
    // all structures present in the code to change this assignment any time.
    //********************************************************************
    //if (nodemap->NumMyElements() || elemap->NumMyElements())
      lin[Comm().MyPID()] = 1;

    Comm().MaxAll(&lin[0],&gin[0],Comm().NumProc());
    lin.clear();

    // build global -> local communicator PID map
    // we need this when calling Broadcast() on lComm later
    int counter = 0;
    for (int i=0; i<Comm().NumProc(); ++i)
    {
      if (gin[i])
        procmap_[i]=counter++;
      else
        procmap_[i]=-1;
    }

    // typecast the Epetra_Comm to Epetra_MpiComm
    Teuchos::RCP<Epetra_Comm> copycomm = Teuchos::rcp(Comm().Clone());
    Epetra_MpiComm* epetrampicomm = dynamic_cast<Epetra_MpiComm*>(copycomm.get());
    if (!epetrampicomm)
      dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");

    // split the communicator into participating and none-participating procs
    int color;
    int key = Comm().MyPID();
    // I am taking part in the new comm if I have any ownership
    if (gin[Comm().MyPID()])
      color = 0;
    // I am not taking part in the new comm
    else
      color = MPI_UNDEFINED;

    // tidy up
    gin.clear();

    // create the local communicator
    MPI_Comm  mpi_global_comm = epetrampicomm->GetMpiComm();
    MPI_Comm  mpi_local_comm;
    MPI_Comm_split(mpi_global_comm,color,key,&mpi_local_comm);

    // create the new Epetra_MpiComm
    if (mpi_local_comm == MPI_COMM_NULL)
      lcomm_ = Teuchos::null;
    else
      lcomm_ = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));

#else  // the easy serial case
    Teuchos::RCP<Epetra_Comm> copycomm = Teuchos::rcp(Comm().Clone());
    Epetra_SerialComm* serialcomm = dynamic_cast<Epetra_SerialComm*>(copycomm.get());
    if (!serialcomm)
      dserror("ERROR: casting Epetra_Comm -> Epetra_SerialComm failed");
    lcomm_ = Teuchos::rcp(new Epetra_SerialComm(*serialcomm));
#endif // #ifdef PARALLEL
  }

  // create interface ghosting
  // (currently, the slave is kept with the standard overlap of one,
  // but the master is made fully redundant, i.e. it is exported to
  // fully overlapping column layout, for the ease of interface search)
  // (the only exceptions are self contact and coupled problems, where
  // also the slave is still made fully redundant)
  CreateInterfaceGhosting();

  // make sure discretization is complete
  Discret().FillComplete(true,false,false);

  // need row and column maps of slave and master nodes / elements / dofs
  // separately so we can easily address them
  UpdateMasterSlaveSets();

  // initialize node data container
  // (include slave side boundary nodes / crosspoints)
  for (int i=0; i<SlaveColNodesBound()->NumMyElements(); ++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %i",gid);
    MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(node);

    //********************************************************
    // NOTE: depending on which kind of node this really is,
    // i.e. mortar, contact or friction node, several derived
    // versions of the InitializeDataContainer() methods will
    // be called here, apart from the base class version.
    //********************************************************

    // initialize container if not yet initialized before
    mnode->InitializeDataContainer();
  }

  //***********************************************************
  // both-sided wear
  // here we need a datacontainer for the masternodes too
  // they have to know their involved nodes/dofs
  //***********************************************************
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(imortar_,"BOTH_SIDED_WEAR") != INPAR::CONTACT::wear_slave)
  {
    for (int i=0; i<MasterRowNodes()->NumMyElements(); ++i) //col //for (int i=0; i<MasterRowNodes()->NumMyElements(); ++i)
    {
      int gid = MasterRowNodes()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %i",gid);
      MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(node);

      //********************************************************
      // NOTE: depending on which kind of node this really is,
      // i.e. mortar, contact or friction node, several derived
      // versions of the InitializeDataContainer() methods will
      // be called here, apart from the base class version.
      //********************************************************

      // initialize container if not yet initialized before
      mnode->InitializeDataContainer();
    }
  }

  // initialize element data container
  for (int i=0; i<SlaveColElements()->NumMyElements(); ++i)
  {
    int gid = SlaveColElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
    MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);

    // initialize container if not yet initialized before
    mele->InitializeDataContainer();
  }

  // communicate quadslave status among ALL processors
  // (not only those participating in interface)
  int localstatus = (int)(quadslave_);
  int globalstatus = 0;
  Comm().SumAll(&localstatus,&globalstatus,1);
  quadslave_ = (bool)(globalstatus);

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble inactive wear right hand side                   farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleInactiveWearRhs(Epetra_Vector& inactiverhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // FIXME It's possible to improve the performance, if only recently active nodes of the inactive node set,
  // i.e. nodes, which were active in the last iteration, are considered. Since you know, that the lagrange
  // multipliers of former inactive nodes are still equal zero.

  Teuchos::RCP<Epetra_Map> inactivenodes  = LINALG::SplitMap(*snoderowmap_, *activenodes_);
  Teuchos::RCP<Epetra_Map> inactivedofs   = LINALG::SplitMap(*sdofrowmap_, *activedofs_);

  for (int i=0;i<inactivenodes->NumMyElements();++i)
  {
    int gid = inactivenodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = static_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleInactiverhs: Node ownership inconsistency!");

    if (Dim() == 2)
    {
      std::vector<int> w_gid(1);
      std::vector<int> w_owner(1);

      // calculate the tangential rhs
      Epetra_SerialDenseVector w_i(1);

      w_owner[0] = cnode->Owner();
      w_i[0]     = - cnode->FriDataPlus().wcurr()[0];    // already negative rhs!!!
      w_gid[0]   = cnode->Dofs()[0];//inactivedofs->GID(2*i);


      if (abs(w_i[0])>1e-12) LINALG::Assemble(inactiverhs, w_i, w_gid, w_owner);
    }
    else if (Dim() == 3)
    {
      std::vector<int> w_gid(1);
      std::vector<int> w_owner(1);

      // calculate the tangential rhs
      Epetra_SerialDenseVector w_i(1);

      w_owner[0] = cnode->Owner();
      w_i[0]     = - cnode->FriDataPlus().wcurr()[0];    // already negative rhs!!!
      w_gid[0]   = inactivedofs->GID(3*i);

      if (abs(w_i[0])>1e-12) LINALG::Assemble(inactiverhs, w_i, w_gid, w_owner);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Assemble wear-cond. right hand side                      farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::AssembleWearCondRhs(Epetra_Vector& rhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // nothing to do if no active nodes
  if (slipnodes_==Teuchos::null)
    return;

  double wcoeff = IParams().get<double>("WEARCOEFF");

  typedef std::map<int,double>::const_iterator CI;

  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* fnode = static_cast<FriNode*>(node);

    if (fnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleWearCondRhs: Node ownership inconsistency!");

    /**************************************************** E-matrix ******/
    if ((fnode->FriDataPlus().GetE()).size()>0)
    {
      std::map<int,double>  emap = fnode->FriDataPlus().GetE()[0];

      for (CI p=emap.begin();p!=emap.end();++p)
      {
        int gid3 = (int)((p->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        FriNode* csnode = static_cast<FriNode*>(snode);

        std::vector<int> w_gid(1);
        std::vector<int> w_owner(1);

        Epetra_SerialDenseVector w_i(1);

        w_owner[0] = fnode->Owner();
        w_i[0]     = - (csnode->FriDataPlus().wcurr()[0]) * (p->second);
        w_gid[0]   = fnode->Dofs()[0];

        if (abs(w_i[0])>1e-12) LINALG::Assemble(rhs, w_i, w_gid, w_owner);
      }
    }

    /**************************************************** T-matrix ******/
    if ((fnode->FriDataPlus().GetT()).size()>0)
    {
      std::map<int,double> tmap = fnode->FriDataPlus().GetT()[0];

      for (CI p=tmap.begin();p!=tmap.end();++p)
      {
        int gid3 = (int)((p->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        FriNode* csnode = static_cast<FriNode*>(snode);

        double lmn = 0.0;
        for (int u=0;u<Dim();++u)
          lmn += (csnode->MoData().n()[u]) * (csnode->MoData().lm()[u]);

        std::vector<int> w_gid(1);
        std::vector<int> w_owner(1);

        Epetra_SerialDenseVector w_i(1);

        w_owner[0] = fnode->Owner();
        w_i[0]     =  wcoeff * lmn * (p->second);
        w_gid[0]   = fnode->Dofs()[0];

        if (abs(w_i[0])>1e-12) LINALG::Assemble(rhs, w_i, w_gid, w_owner);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 |  initialize / reset interface for wear                    farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)

  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CoNode* node = static_cast<CONTACT::CoNode*>(idiscret_->lColNode(i));

    // reset feasible projection and segmentation status
    node->HasProj() = false;
    node->HasSegment() = false;

    if (friction_)
    {
      FriNode* frinode = static_cast<FriNode*>(node);

      // reset nodal mechanical dissipation
      frinode->MechDiss() = 0.0;

      // reset matrix B quantities
      frinode->GetBNodes().clear();

      // reset nodal B maps
      for (int j=0;j<(int)((frinode->GetB()).size());++j)
        (frinode->GetB())[j].clear();

      (frinode->GetB()).resize(0);

    }
  }

  //******************************************************
  // for both-sided wear
  //******************************************************
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(imortar_,"BOTH_SIDED_WEAR") != INPAR::CONTACT::wear_slave)
  {
    for (int i=0;i<MasterRowNodes()->NumMyElements();++i) //col //for (int i=0;i<MasterRowNodes()->NumMyElements();++i)
    {
      int gid = MasterRowNodes()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = static_cast<CoNode*>(node);

      if (cnode->IsSlave() ==false)
      {
        // reset nodal Mortar maps
        for (int j=0;j<(int)((cnode->CoData().GetD2()).size());++j)
          (cnode->CoData().GetD2())[j].clear();

        (cnode->CoData().GetD2()).resize(0);
      }
    }
  }

  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i=0;i<SlaveColNodesBound()->NumMyElements();++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);


    // reset nodal Mortar maps
    for (int j=0;j<(int)((cnode->MoData().GetD()).size());++j)
      (cnode->MoData().GetD())[j].clear();
    for (int j=0;j<(int)((cnode->MoData().GetM()).size());++j)
      (cnode->MoData().GetM())[j].clear();
    for (int j=0;j<(int)((cnode->MoData().GetMmod()).size());++j)
      (cnode->MoData().GetMmod())[j].clear();

    (cnode->MoData().GetD()).resize(0);
    (cnode->MoData().GetM()).resize(0);
    (cnode->MoData().GetMmod()).resize(0);

    // reset nodal scaling factor
    (cnode->MoData().GetScale())=0.0;
    (cnode->CoData().GetDerivScale()).clear();

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((cnode->CoData().GetDerivN()).size());++j)
      (cnode->CoData().GetDerivN())[j].clear();
    (cnode->CoData().GetDerivN()).resize(0);

    // reset derivative maps of tangent vectors
    for (int j=0;j<(int)((cnode->CoData().GetDerivTxi()).size());++j)
      (cnode->CoData().GetDerivTxi())[j].clear();
    (cnode->CoData().GetDerivTxi()).resize(0);
    for (int j=0;j<(int)((cnode->CoData().GetDerivTeta()).size());++j)
      (cnode->CoData().GetDerivTeta())[j].clear();
    (cnode->CoData().GetDerivTeta()).resize(0);

    // reset derivative map of Mortar matrices
    (cnode->CoData().GetDerivD()).clear();
    (cnode->CoData().GetDerivM()).clear();

    // reset nodal weighted gap and derivative
    cnode->CoData().Getg() = 1.0e12;
    (cnode->CoData().GetDerivG()).clear();

    // reset weighted wear increment and derivative
    // only for implicit wear algorithm
    if(wearimpl_)
    {
      (cnode->CoData().GetDerivW()).clear();
      (cnode->CoData().GetDerivWlm()).clear();
    }

    // reset derivative map of lagrange multipliers
    for (int j=0; j<(int)((cnode->CoData().GetDerivZ()).size()); ++j)
      (cnode->CoData().GetDerivZ())[j].clear();
    (cnode->CoData().GetDerivZ()).resize(0);

    if (friction_)
    {
      FriNode* frinode = static_cast<FriNode*>(cnode);

      // reset SNodes and Mnodes
      frinode->FriData().GetSNodes().clear();
      frinode->FriData().GetMNodes().clear();

#ifdef OBJECTVARSLIPINCREMENT
      // reset jump deriv.
      for (int j=0 ;j<(int)((frinode->FriData().GetDerivVarJump()).size()); ++j )
        (frinode->FriData().GetDerivVarJump())[j].clear();

      (frinode->FriData().GetDerivVarJump()).resize(2);

      // reset jumps
      frinode->FriData().jump_var()[0]=0.0;
      frinode->FriData().jump_var()[1]=0.0;
#endif
      if(weardiscr_)
      {
        // reset nodal Mortar wear maps
        for (int j=0;j<(int)((frinode->FriDataPlus().GetT()).size());++j)
          (frinode->FriDataPlus().GetT())[j].clear();
        for (int j=0;j<(int)((frinode->FriDataPlus().GetE()).size());++j)
          (frinode->FriDataPlus().GetE())[j].clear();

        (frinode->FriDataPlus().GetT()).resize(0);
        (frinode->FriDataPlus().GetE()).resize(0);

        (frinode->FriDataPlus().GetDerivTw()).clear();
        (frinode->FriDataPlus().GetDerivE()).clear();

        (frinode->CoData().GetDerivGW()).clear();
      }

      if (tsi_)
      {
        // reset matrix A quantities
        frinode->FriDataPlus().GetANodes().clear();

        // reset nodal A maps
        for (int j=0;j<(int)((frinode->FriDataPlus().GetA()).size());++j)
          (frinode->FriDataPlus().GetA())[j].clear();

        (frinode->FriDataPlus().GetA()).resize(0);
      }

      if (wear_)
      {
        // reset wear increment
        frinode->FriDataPlus().DeltaWear() = 0.0;

        // reset abs. wear.
        // for impl. wear algor. the abs. wear equals the
        // delta-wear
        if(wearimpl_)
          frinode->FriDataPlus().Wear() = 0.0;
      }
    }
  }

  //**********************************************************************
  // In general, it is sufficient to reset search candidates only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need to reset the search candidates of
  // all slave elements in the fully overlapping column map there. This
  // is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact())
  {
    // loop over all elements to reset candidates / search lists
    // (use fully overlapping column map of S+M elements)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      DRT::Element* ele = idiscret_->lColElement(i);
      MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);

      mele->MoData().SearchElements().resize(0);

      // dual shape function coefficient matrix
      mele->MoData().ResetDualShape();
      mele->MoData().ResetDerivDualShape();
    }
  }
  else
  {
    // loop over all elements to reset candidates / search lists
    // (use standard slave column map)
    for (int i=0;i<SlaveColElements()->NumMyElements();++i)
    {
      int gid = SlaveColElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
      MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);

      mele->MoData().SearchElements().resize(0);

      // dual shape function coefficient matrix
      mele->MoData().ResetDualShape();
      mele->MoData().ResetDerivDualShape();
    }
  }

  // reset s/m pairs and intcell counters
  smpairs_    = 0;
  smintpairs_ = 0;
  intcells_   = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  create snode n_                                          farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::SplitSlaveDofs()
{
  // get out of here if active set is empty
  if (snoderowmap_==Teuchos::null)
  {
    sndofmap_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return;
  }

  else if (snoderowmap_->NumGlobalElements()==0)
  {
    sndofmap_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return;
  }

  // define local variables
  int countN=0;
  std::vector<int> myNgids(snoderowmap_->NumMyElements());

  // dimension check
  double dimcheck =(sdofrowmap_->NumGlobalElements())/(snoderowmap_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitSlaveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    // add first dof to Nmap
    myNgids[countN] = cnode->Dofs()[0];
    ++countN;
  }

  // resize the temporary vectors
  myNgids.resize(countN);

  // communicate countN and countT among procs
  int gcountN;
  Comm().SumAll(&countN,&gcountN,1);

  // check global dimensions
  if ((gcountN)!=snoderowmap_->NumGlobalElements())
    dserror("ERROR: SplitSlaveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  sndofmap_ = Teuchos::rcp(new Epetra_Map(gcountN,countN,&myNgids[0],0,Comm()));

  return;
}

/*----------------------------------------------------------------------*
 |  compute element areas (public)                            popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::SetElementAreas()
{
  //**********************************************************************
  // In general, it is sufficient to compute element areas only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need the element areas of all elements
  // (slave and master) in the fully overlapping column map there. At the
  // same time we initialize the element data containers for self contact.
  // This is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact())
  {
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      MORTAR::MortarElement* element = static_cast<MORTAR::MortarElement*>(idiscret_->lColElement(i));
      element->InitializeDataContainer();
      element->MoData().Area()=element->ComputeArea();
    }
  }
  else
  {
    // for both-sided wear we need element areas of master elements --> colelements
    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(imortar_,"BOTH_SIDED_WEAR") != INPAR::CONTACT::wear_slave)
    {
      for (int i=0;i<idiscret_->NumMyColElements();++i)
      {
        MORTAR::MortarElement* element = static_cast<MORTAR::MortarElement*>(idiscret_->lColElement(i));
        element->InitializeDataContainer();
        element->MoData().Area()=element->ComputeArea();
      }
    }
    else
    {
      //refer call back to base class version
      MORTAR::MortarInterface::SetElementAreas();
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  update wear set (dofs)                                   farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::WearInterface::UpdateWSets(int offset_if, int maxdofwear)
{
  //********************************************************************
  // WEAR DOFS --  one per node
  //********************************************************************
  // NOTE: we want no gap between the displacement dofs and the newly
  // defined Lagrange multiplier dofs!! Thus, if the maximum displacement
  // dof is 12.345, we want the LM dofs to start with 12.346. This can
  // be readily achieved, because we know that the lmdofmap will have
  // the same parallel distribution as the slavedofrowmap. The only
  // thing we need to take care of is to avoid overlapping of the LM
  // dofs among different processors. Therefore, the total number of
  // slave nodes (and thus LM nodes) of each processor is communicated
  // to ALL other processors and an offset is then determined for each
  // processor based on this information.
  //********************************************************************
  // temporary vector of W dofs
  std::vector<int> wdof;

  // gather information over all procs
  std::vector<int> localnumwdof(Comm().NumProc());
  std::vector<int> globalnumlmdof(Comm().NumProc());
  localnumwdof[Comm().MyPID()] = (int)( (sdofrowmap_->NumMyElements())/Dim() );
  Comm().SumAll(&localnumwdof[0],&globalnumlmdof[0],Comm().NumProc());

  // compute offet for LM dof initialization for all procs
  int offset = 0;
  for (int k=0;k<Comm().MyPID();++k)
    offset += globalnumlmdof[k];

  // loop over all slave dofs and initialize LM dofs
  for (int i=0; i<(int)( (sdofrowmap_->NumMyElements())/Dim() ); ++i)
    wdof.push_back(maxdofwear + 1 + offset_if + offset + i);

  // create interface w map
  // (if maxdofglobal_ == 0, we do not want / need this)
  if (maxdofwear>0)
    wdofmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)wdof.size(),&wdof[0],0,Comm()));

  return;
}
