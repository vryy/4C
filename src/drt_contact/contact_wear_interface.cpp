/*----------------------------------------------------------------------*/
/*!
\file contact_wear_interface.cpp

\brief Wear interface implementation.

\level 2

\maintainer Philipp Farah
*/
/*----------------------------------------------------------------------*/

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
#include <Epetra_FEVector.h>

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 09/13|
 *----------------------------------------------------------------------*/
WEAR::WearInterface::WearInterface(
    const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr,
    const int id, const Epetra_Comm& comm,
    const int dim,
    const Teuchos::ParameterList& icontact,
    bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant) :
CONTACT::CoInterface(idata_ptr,id,comm,dim,icontact,selfcontact,redundant),
wear_(false),
wearimpl_(false),
wearpv_(false),
wearboth_(false),
sswear_(DRT::INPUT::IntegralValue<int>(icontact,"SSWEAR"))
{
  // set wear contact status
  INPAR::WEAR::WearType wtype =
      DRT::INPUT::IntegralValue<INPAR::WEAR::WearType>(icontact,"WEARTYPE");

  INPAR::WEAR::WearTimInt wtimint =
      DRT::INPUT::IntegralValue<INPAR::WEAR::WearTimInt>(icontact,"WEARTIMINT");

  INPAR::WEAR::WearSide wside =
      DRT::INPUT::IntegralValue<INPAR::WEAR::WearSide>(icontact,"WEAR_SIDE");

  INPAR::WEAR::WearLaw wlaw =
      DRT::INPUT::IntegralValue<INPAR::WEAR::WearLaw>(icontact,"WEARLAW");

  if (wlaw != INPAR::WEAR::wear_none)
    wear_ = true;

  if (wtimint == INPAR::WEAR::wear_impl)
    wearimpl_ = true;

  // set wear contact discretization
  if (wtype == INPAR::WEAR::wear_primvar)
    wearpv_ = true;

  // set wear contact discretization
  if (wside == INPAR::WEAR::wear_both)
    wearboth_ = true;

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble Mortar wear matrices                            farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleTE(LINALG::SparseMatrix& tglobal,
                                      LINALG::SparseMatrix& eglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  //nodes for loop
  Teuchos::RCP<Epetra_Map> considerednodes;

  // nothing to do if no active nodes
  if(sswear_)
  {
    if (activenodes_==Teuchos::null)
      return;
      considerednodes = activenodes_;
  }
  else
  {
    if (slipnodes_==Teuchos::null)
      return;
      considerednodes = slipnodes_;
  }

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<considerednodes->NumMyElements();++i)
  {
    int gid = considerednodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (fnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleTE: Node ownership inconsistency!");

    /**************************************************** T-matrix ******/
    if ((fnode->WearData().GetT()).size()>0)
    {
      std::vector<std::map<int,double> >& tmap = fnode->WearData().GetT();
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
    if ((fnode->WearData().GetE()).size()>0)
    {
      std::vector<std::map<int,double> >& emap = fnode->WearData().GetE();
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
        if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
        {
          if (abs(val)>1.0e-12) eglobal.Assemble(val,row,col);
          ++k;
        }
        else if (WearShapeFcn() == INPAR::WEAR::wear_shape_dual)
        {
          if (col==row)
            if (abs(val)>1.0e-12) eglobal.Assemble(val,row,col);
          ++k;
        }
        else
          dserror("Chosen wear shape function not supported!");
      }

      if (k!=colsize)
        dserror("ERROR: AssembleTE: k = %i but colsize = %i",k,colsize);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble Mortar wear matrices (for master side)          farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleTE_Master(LINALG::SparseMatrix& tglobal,
                                               LINALG::SparseMatrix& eglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  if (!(wearboth_ and wearpv_))
    dserror("ERROR: AssembleTE_Master only for discr both-sided wear!");

  //*******************************************************
  // assemble second D matrix for both-sided wear :
  // Loop over allreduced map, so that any proc join
  // the loop. Do non-local assembly by FEAssemble to
  // allow owned and ghosted nodes assembly. In integrator
  // only the nodes associated to the owned slave element
  // are filled with entries. Therefore, the integrated
  // entries are unique on nodes!
  //*******************************************************

  // nothing to do if no active nodes
  if (slipmasternodes_==Teuchos::null)
    return;

  const Teuchos::RCP<Epetra_Map> slmasternodes = LINALG::AllreduceEMap(*(slipmasternodes_));

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<slmasternodes->NumMyElements();++i)
  {
    int gid = slmasternodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    /**************************************************** T-matrix ******/
    if ((fnode->WearData().GetT()).size()>0)
    {
      std::vector<std::map<int,double> >& tmap = fnode->WearData().GetT();
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
        if (abs(val)>1.0e-12) tglobal.FEAssemble(val, row, col);

        ++k;
      }

      if (k!=colsize)
        dserror("ERROR: AssembleTE: k = %i but colsize = %i",k,colsize);

    }

    /**************************************************** E-matrix ******/
    if ((fnode->WearData().GetE()).size()>0)
    {
      std::vector<std::map<int,double> >& emap = fnode->WearData().GetE();
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
        if (WearShapeFcn() == INPAR::WEAR::wear_shape_standard)
        {
          if (abs(val)>1.0e-12) eglobal.FEAssemble(val,row,col);
          ++k;
        }
        else if (WearShapeFcn() == INPAR::WEAR::wear_shape_dual)
        {
          if (col==row)
            if (abs(val)>1.0e-12) eglobal.FEAssemble(val,row,col);
          ++k;
        }
        else
          dserror("Choosen wear shape function not supported!");
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
void WEAR::WearInterface::AssembleLinT_D(LINALG::SparseMatrix& lintglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  //nodes for loop
  Teuchos::RCP<Epetra_Map> considerednodes;

  // nothing to do if no active nodes
  if(sswear_)
  {
    if (activenodes_==Teuchos::null)
      return;
      considerednodes = activenodes_;
  }
  else
  {
    if (slipnodes_==Teuchos::null)
      return;
      considerednodes = slipnodes_;
  }

  /**********************************************************************/
  // we have: T_wj,c with j = Lagrange multiplier slave dof
  //                 with w = wear slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinT)_kc = T_wj,c * z_j
  /**********************************************************************/

  for (int j=0;j<considerednodes->NumMyElements();++j)
  {
    int gid = considerednodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    // Mortar matrix Tw derivatives
    std::map<int,std::map<int,double> >& tderiv = fnode->WearData().GetDerivTw();

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
      CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

      // current Lagrange multipliers
      double lmn = 0.0;
      if (Dim()==2)
        lmn=(csnode->MoData().lm()[0]) * (csnode->MoData().n()[0]) + (csnode->MoData().lm()[1]) * (csnode->MoData().n()[1]);
      else if (Dim()==3)
        lmn=(csnode->MoData().lm()[0]) * (csnode->MoData().n()[0]) + (csnode->MoData().lm()[1]) * (csnode->MoData().n()[1]) + (csnode->MoData().lm()[2]) * (csnode->MoData().n()[2]);
      else
        dserror("False Dimension!");

      // Mortar matrix T derivatives
      std::map<int,double>& thisdderive = fnode->WearData().GetDerivTw()[sgid];
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
  for (int j=0;j<considerednodes->NumMyElements();++j)
  {
    int gid = considerednodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (fnode->WearData().GetT().size()>0)
    {
      // map iterator
      typedef std::map<int,double>::const_iterator           CI;
      typedef GEN::pairedvector<int,double>::const_iterator _CI;

      std::map<int,double>& nmap = fnode->WearData().GetT()[0];

      // loop over col entries
      for (CI z=nmap.begin();z!=nmap.end();++z)
      {
        //std::cout << "t-irst= " << z->first << std::endl;
        int gid3 = (int)((z->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

        for (int u=0;u<Dim();++u)
        {
          GEN::pairedvector<int,double>& numap = csnode->CoData().GetDerivN()[u];
          double lmu = csnode->MoData().lm()[u];

          // multiply T-column entry with lin n*lambda
          for (_CI b=numap.begin(); b!=numap.end(); ++b)
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
 |  Assemble matrix LinT containing disp derivatives         farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleLinT_D_Master(LINALG::SparseMatrix& lintglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (slipmasternodes_==Teuchos::null)
    return;

  /**********************************************************************/
  // we have: T_wj,c with j = Lagrange multiplier slave dof
  //                 with w = wear slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinT)_kc = T_wj,c * z_j
  /**********************************************************************/
  const Teuchos::RCP<Epetra_Map> slmasternodes = LINALG::AllreduceEMap(*(slipmasternodes_));

  for (int j=0;j<slmasternodes->NumMyElements();++j)
  {
    int gid = slmasternodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    // Mortar matrix Tw derivatives
    std::map<int,std::map<int,double> >& tderiv = fnode->WearData().GetDerivTw();

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
      CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

      // current Lagrange multipliers
      double lmn = 0.0;
      if (Dim()==2)
        lmn=(csnode->MoData().lm()[0]) * (csnode->MoData().n()[0]) + (csnode->MoData().lm()[1]) * (csnode->MoData().n()[1]);
      else if (Dim()==3)
        lmn=(csnode->MoData().lm()[0]) * (csnode->MoData().n()[0]) + (csnode->MoData().lm()[1]) * (csnode->MoData().n()[1]) + (csnode->MoData().lm()[2]) * (csnode->MoData().n()[2]);
      else
        dserror("False Dimension!");

      // Mortar matrix T derivatives
      std::map<int,double>& thisdderive = fnode->WearData().GetDerivTw()[sgid];
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
  for (int j=0;j<slmasternodes->NumMyElements();++j)
  {
    int gid = slmasternodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (fnode->WearData().GetT().size()>0)
    {
      // map iterator
      typedef std::map<int,double>::const_iterator           CI;
      typedef GEN::pairedvector<int,double>::const_iterator _CI;

      std::map<int,double>& nmap = fnode->WearData().GetT()[0];

      // loop over col entries
      for (CI z=nmap.begin();z!=nmap.end();++z)
      {
        //std::cout << "t-irst= " << z->first << std::endl;
        int gid3 = (int)((z->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

        for (int u=0;u<Dim();++u)
        {
          GEN::pairedvector<int,double>& numap = csnode->CoData().GetDerivN()[u];
          double lmu = csnode->MoData().lm()[u];

          // multiply T-column entry with lin n*lambda
          for (_CI b=numap.begin(); b!=numap.end(); ++b)
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
void WEAR::WearInterface::AssembleLinE_D(LINALG::SparseMatrix& lineglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  //nodes for loop
  Teuchos::RCP<Epetra_Map> considerednodes;

  // nothing to do if no active nodes
  if(sswear_)
  {
    if (activenodes_==Teuchos::null)
      return;
    considerednodes = activenodes_;
  }
  else
  {
    if (slipnodes_==Teuchos::null)
      return;
    considerednodes = slipnodes_;
  }


  /**********************************************************************/
  // we have: E_wj,c with j = Lagrange multiplier slave dof
  //                 with w = wear slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinE)_kc = T_wj,c * z_j
  /**********************************************************************/

  // loop over all LM slave nodes (row map)
  for (int j=0;j<considerednodes->NumMyElements();++j)
  {
    int gid = considerednodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    // Mortar matrix Tw derivatives
    std::map<int,std::map<int,double> >& ederiv = fnode->WearData().GetDerivE();

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
      CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

      // current wear - wear from last converged iteration step (partitioned scheme)
      double w = 0.0;
      w=(csnode->WearData().wcurr()[0] + csnode->WearData().wold()[0]);

      // Mortar matrix T derivatives
      std::map<int,double>& thisdderive = fnode->WearData().GetDerivE()[sgid];
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
        if (abs(val)>1.0e-15) lineglobal.FEAssemble(val,row,col);
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
 |  Assemble matrix LinE containing disp derivatives         farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleLinE_D_Master(LINALG::SparseMatrix& lineglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipmasternodes_->NumGlobalElements()==0)
    return;

  /**********************************************************************/
  // we have: E_wj,c with j = Lagrange multiplier slave dof
  //                 with w = wear slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinE)_kc = T_wj,c * z_j
  /**********************************************************************/
  const Teuchos::RCP<Epetra_Map> slmasternodes = LINALG::AllreduceEMap(*(slipmasternodes_));

  // loop over all LM slave nodes (row map)
  for (int j=0;j<slmasternodes->NumMyElements();++j)
  {
    int gid = slmasternodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    // Mortar matrix Tw derivatives
    std::map<int,std::map<int,double> >& ederiv = fnode->WearData().GetDerivE();

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
      CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

      // current wear - wear from last converged iteration step (partitioned scheme)
      double w = 0.0;
      w=(csnode->WearData().wcurr()[0] + csnode->WearData().wold()[0]);

      // Mortar matrix T derivatives
      std::map<int,double>& thisdderive = fnode->WearData().GetDerivE()[sgid];
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
void WEAR::WearInterface::AssembleLinT_LM(LINALG::SparseMatrix& lintglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  //nodes for loop
  Teuchos::RCP<Epetra_Map> considerednodes;

  // nothing to do if no active nodes
  if(sswear_)
  {
    if (activenodes_==Teuchos::null)
      return;
    considerednodes = activenodes_;
  }
  else
  {
    if (slipnodes_==Teuchos::null)
      return;
    considerednodes = slipnodes_;
  }

  //typedef std::map<int,double>::const_iterator CI;

  // loop over all LM slave nodes (row map)
  for (int j=0;j<considerednodes->NumMyElements();++j)
  {
    int gid = considerednodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    typedef std::map<int,double>::const_iterator CI;

    if (fnode->WearData().GetT().size()>0)
    {
      // column entries for row f
      std::map<int,double>& fmap = fnode->WearData().GetT()[0];

      for (CI p=fmap.begin();p!=fmap.end();++p)
      {
        int gid2 = (int)((p->first)/Dim());
        DRT::Node* node2 = idiscret_->gNode(gid2);
        if (!node2) dserror("ERROR: Cannot find node with gid %",gid2);
        CONTACT::FriNode* jnode = dynamic_cast<CONTACT::FriNode*>(node2);

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
 |  Assemble matrix LinT containing lm derivatives           farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleLinT_LM_Master(LINALG::SparseMatrix& lintglobal)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (slipmasternodes_==Teuchos::null)
    return;

  const Teuchos::RCP<Epetra_Map> slmasternodes = LINALG::AllreduceEMap(*(slipmasternodes_));


  //typedef std::map<int,double>::const_iterator CI;

  // loop over all LM slave nodes (row map)
  for (int j=0;j<slmasternodes->NumMyElements();++j)
  {
    int gid = slmasternodes->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    typedef std::map<int,double>::const_iterator CI;

    if (fnode->WearData().GetT().size()>0)
    {
      // column entries for row f
      std::map<int,double>& fmap = fnode->WearData().GetT()[0];

      for (CI p=fmap.begin();p!=fmap.end();++p)
      {
        int gid2 = (int)((p->first)/Dim());
        DRT::Node* node2 = idiscret_->gNode(gid2);
        if (!node2) dserror("ERROR: Cannot find node with gid %",gid2);
        CONTACT::FriNode* jnode = dynamic_cast<CONTACT::FriNode*>(node2);

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
 |  evaluate nodal normals (public)                          farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::EvaluateNodalNormals() const
{
  // call mortar function
  MORTAR::MortarInterface::EvaluateNodalNormals();

  // for both-sided discrete wear
  if(wearboth_ == true and
     wearpv_   == true)
  {
    for(int i=0; i<mnoderowmap_->NumMyElements();++i)
    {
      int gid = mnoderowmap_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(node);

      // build averaged normal at each master node
      mrtrnode->BuildAveragedNormal();
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  export nodal normals (public)                            farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::ExportNodalNormals()
{
  // call contact function
  CONTACT::CoInterface::ExportNodalNormals();

  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > triad;

  std::map<int,std::vector<int> > n_x_key;
  std::map<int,std::vector<int> > n_y_key;
  std::map<int,std::vector<int> > n_z_key;
  std::map<int,std::vector<int> > txi_x_key;
  std::map<int,std::vector<int> > txi_y_key;
  std::map<int,std::vector<int> > txi_z_key;
  std::map<int,std::vector<int> > teta_x_key;
  std::map<int,std::vector<int> > teta_y_key;
  std::map<int,std::vector<int> > teta_z_key;

  std::map<int,std::vector<double> > n_x_val;
  std::map<int,std::vector<double> > n_y_val;
  std::map<int,std::vector<double> > n_z_val;
  std::map<int,std::vector<double> > txi_x_val;
  std::map<int,std::vector<double> > txi_y_val;
  std::map<int,std::vector<double> > txi_z_val;
  std::map<int,std::vector<double> > teta_x_val;
  std::map<int,std::vector<double> > teta_y_val;
  std::map<int,std::vector<double> > teta_z_val;

  GEN::pairedvector<int,double>::iterator iter;

  // --------------------------------------------------------------------------------------
  // for both-sided discrete wear we need the same normal information on the master side:
  // --------------------------------------------------------------------------------------
  if(wearboth_ and wearpv_)
  {
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(mnoderowmap_));

    // build info on row map
     for(int i=0; i<mnoderowmap_->NumMyElements();++i)
     {
       int gid = mnoderowmap_->GID(i);
       DRT::Node* node = idiscret_->gNode(gid);
       if (!node) dserror("ERROR: Cannot find node with gid %",gid);
       CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

       // fill nodal matrix
       Teuchos::RCP<Epetra_SerialDenseMatrix> loc = Teuchos::rcp(new Epetra_SerialDenseMatrix(3,3));
       (*loc)(0,0) = cnode->MoData().n()[0];
       (*loc)(1,0) = cnode->MoData().n()[1];
       (*loc)(2,0) = cnode->MoData().n()[2];
       (*loc)(0,1) = cnode->CoData().txi()[0];
       (*loc)(1,1) = cnode->CoData().txi()[1];
       (*loc)(2,1) = cnode->CoData().txi()[2];
       (*loc)(0,2) = cnode->CoData().teta()[0];
       (*loc)(1,2) = cnode->CoData().teta()[1];
       (*loc)(2,2) = cnode->CoData().teta()[2];

       triad[gid] = loc;

       // fill nodal derivative vectors
       std::vector<GEN::pairedvector<int,double> >& derivn    = cnode->CoData().GetDerivN();
       std::vector<GEN::pairedvector<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
       std::vector<GEN::pairedvector<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

       for(iter=derivn[0].begin();iter!=derivn[0].end();++iter)
       {
         n_x_key[gid].push_back(iter->first);
         n_x_val[gid].push_back(iter->second);
       }
       for(iter=derivn[1].begin();iter!=derivn[1].end();++iter)
       {
         n_y_key[gid].push_back(iter->first);
         n_y_val[gid].push_back(iter->second);
       }
       for(iter=derivn[2].begin();iter!=derivn[2].end();++iter)
       {
         n_z_key[gid].push_back(iter->first);
         n_z_val[gid].push_back(iter->second);
       }

       for(iter=derivtxi[0].begin();iter!=derivtxi[0].end();++iter)
       {
         txi_x_key[gid].push_back(iter->first);
         txi_x_val[gid].push_back(iter->second);
       }
       for(iter=derivtxi[1].begin();iter!=derivtxi[1].end();++iter)
       {
         txi_y_key[gid].push_back(iter->first);
         txi_y_val[gid].push_back(iter->second);
       }
       for(iter=derivtxi[2].begin();iter!=derivtxi[2].end();++iter)
       {
         txi_z_key[gid].push_back(iter->first);
         txi_z_val[gid].push_back(iter->second);
       }

       for(iter=derivteta[0].begin();iter!=derivteta[0].end();++iter)
       {
         teta_x_key[gid].push_back(iter->first);
         teta_x_val[gid].push_back(iter->second);
       }
       for(iter=derivteta[1].begin();iter!=derivteta[1].end();++iter)
       {
         teta_y_key[gid].push_back(iter->first);
         teta_y_val[gid].push_back(iter->second);
       }
       for(iter=derivteta[2].begin();iter!=derivteta[2].end();++iter)
       {
         teta_z_key[gid].push_back(iter->first);
         teta_z_val[gid].push_back(iter->second);
       }
     }

     // communicate from master node row to column map
     DRT::Exporter ex(*mnoderowmap_,*masternodes,Comm());
     ex.Export(triad);

     ex.Export(n_x_key);
     ex.Export(n_x_val);
     ex.Export(n_y_key);
     ex.Export(n_y_val);
     ex.Export(n_z_key);
     ex.Export(n_z_val);

     ex.Export(txi_x_key);
     ex.Export(txi_x_val);
     ex.Export(txi_y_key);
     ex.Export(txi_y_val);
     ex.Export(txi_z_key);
     ex.Export(txi_z_val);

     ex.Export(teta_x_key);
     ex.Export(teta_x_val);
     ex.Export(teta_y_key);
     ex.Export(teta_y_val);
     ex.Export(teta_z_key);
     ex.Export(teta_z_val);

     // extract info on column map
     for(int i=0; i<masternodes->NumMyElements();++i)
     {
       // only do something for ghosted nodes
       int gid = masternodes->GID(i);
       DRT::Node* node = idiscret_->gNode(gid);
       if (!node) dserror("ERROR: Cannot find node with gid %",gid);
       CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
       int linsize = cnode->GetLinsize()+(int)(n_x_key[gid].size());

       if (cnode->Owner()==Comm().MyPID())
         continue;

       // extract info
       Teuchos::RCP<Epetra_SerialDenseMatrix> loc = triad[gid];
       cnode->MoData().n()[0]    = (*loc)(0,0);
       cnode->MoData().n()[1]    = (*loc)(1,0);
       cnode->MoData().n()[2]    = (*loc)(2,0);
       cnode->CoData().txi()[0]  = (*loc)(0,1);
       cnode->CoData().txi()[1]  = (*loc)(1,1);
       cnode->CoData().txi()[2]  = (*loc)(2,1);
       cnode->CoData().teta()[0] = (*loc)(0,2);
       cnode->CoData().teta()[1] = (*loc)(1,2);
       cnode->CoData().teta()[2] = (*loc)(2,2);

       // extract derivative info
       std::vector<GEN::pairedvector<int,double> >& derivn    = cnode->CoData().GetDerivN();
       std::vector<GEN::pairedvector<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
       std::vector<GEN::pairedvector<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

       for (int k=0;k<(int)(derivn.size());++k)
         derivn[k].clear();
       derivn.resize(3,linsize);
       for (int k=0;k<(int)(derivtxi.size());++k)
         derivtxi[k].clear();
       derivtxi.resize(3,linsize);
       for (int k=0;k<(int)(derivteta.size());++k)
         derivteta[k].clear();
       derivteta.resize(3,linsize);

       cnode->CoData().GetDerivN()[0].resize(linsize);
       cnode->CoData().GetDerivN()[1].resize(linsize);
       cnode->CoData().GetDerivN()[2].resize(linsize);

       cnode->CoData().GetDerivTxi()[0].resize(linsize);
       cnode->CoData().GetDerivTxi()[1].resize(linsize);
       cnode->CoData().GetDerivTxi()[2].resize(linsize);

       cnode->CoData().GetDerivTeta()[0].resize(linsize);
       cnode->CoData().GetDerivTeta()[1].resize(linsize);
       cnode->CoData().GetDerivTeta()[2].resize(linsize);

       for (int k=0;k<(int)(n_x_key[gid].size());++k)
         (cnode->CoData().GetDerivN()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
       for (int k=0;k<(int)(n_y_key[gid].size());++k)
         (cnode->CoData().GetDerivN()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
       for (int k=0;k<(int)(n_z_key[gid].size());++k)
         (cnode->CoData().GetDerivN()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

       for (int k=0;k<(int)(txi_x_key[gid].size());++k)
         (cnode->CoData().GetDerivTxi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
       for (int k=0;k<(int)(txi_y_key[gid].size());++k)
         (cnode->CoData().GetDerivTxi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
       for (int k=0;k<(int)(txi_z_key[gid].size());++k)
         (cnode->CoData().GetDerivTxi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

       for (int k=0;k<(int)(teta_x_key[gid].size());++k)
         (cnode->CoData().GetDerivTeta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
       for (int k=0;k<(int)(teta_y_key[gid].size());++k)
         (cnode->CoData().GetDerivTeta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
       for (int k=0;k<(int)(teta_z_key[gid].size());++k)
         (cnode->CoData().GetDerivTeta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
     }

     // free memory
     triad.clear();

     n_x_key.clear();
     n_y_key.clear();
     n_z_key.clear();
     txi_x_key.clear();
     txi_y_key.clear();
     txi_z_key.clear();
     teta_x_key.clear();
     teta_y_key.clear();
     teta_z_key.clear();

     n_x_val.clear();
     n_y_val.clear();
     n_z_val.clear();
     txi_x_val.clear();
     txi_y_val.clear();
     txi_z_val.clear();
     teta_x_val.clear();
     teta_y_val.clear();
     teta_z_val.clear();

  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ derivatives          farah 09/13|
 |  PS: "AssembleS" is an outdated name which could make                |
 |  you confused.                                                       |
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleS(LINALG::SparseMatrix& sglobal)
{
  // call contact function
  CONTACT::CoInterface::AssembleS(sglobal);

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
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    // prepare assembly
    std::map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    /*************************************************************************
     * Wear implicit linearization   --> obviously, we need a new linear.    *
     *************************************************************************/
    if(wearimpl_ and !wearpv_)
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
 |  Assemble matrix S containing gap g~ w derivatives       farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleLinG_W(LINALG::SparseMatrix& sglobal)
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
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

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
void WEAR::WearInterface::AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
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

  // information from interface contact parameter list
  const double frcoeff = IParams().get<double>("FRCOEFF");
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");

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
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

      // prepare assembly, get information from node
      std::vector<GEN::pairedvector<int,double> > dnmap = cnode->CoData().GetDerivN();
      std::vector<GEN::pairedvector<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<GEN::pairedvector<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
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

      // more information from node
      double* n = cnode->MoData().n();
      double* z = cnode->MoData().lm();
      double& wgap = cnode->CoData().Getg();

      // iterator for maps
      std::map<int,double>::iterator           colcurr;
      GEN::pairedvector<int,double>::iterator _colcurr;

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
      double znor = 0;
      double jumptxi = 0;
      double jumpteta = 0;
      double* jump = cnode->FriData().jump();
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();

      for (int i=0;i<Dim();i++)
        znor += n[i]*z[i];

      // for slip
      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];
      }
      else
      {
        // more information from node
        for (int i=0;i<Dim();i++)
        {
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
      }

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

      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
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
          // linearization of normal direction *****************************************
          // loop over all entries of the current derivative map
          for (_colcurr=dnmap[dim].begin();_colcurr!=dnmap[dim].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi=0.0;
            valtxi = - frcoeff * z[dim] * _colcurr->second * ct * jumptxi;
            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
              double valteta=0.0;
              valteta = - frcoeff * z[dim] * _colcurr->second * ct * jumpteta;
              // do not assemble zeros into matrix
              if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }
        } // loop over all dimensions
      }
      else // std slip
      {
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
          for (_colcurr=dtximap[dim].begin();_colcurr!=dtximap[dim].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi=0.0;
            valtxi = - frcoeff*(znor-cn*wgap) * ct * jump[dim] * _colcurr->second;

            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
          }
          // linearization second tangential direction *********************************
          if (Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr=dtetamap[dim].begin();_colcurr!=dtetamap[dim].end();++_colcurr)
            {
              int col = _colcurr->first;
              double valteta=0.0;
              valteta = - frcoeff * (znor-cn*wgap) * ct * jump[dim] * _colcurr->second;

              // do not assemble zeros into matrix
              if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }
          // linearization of normal direction *****************************************
          // loop over all entries of the current derivative map
          for (_colcurr=dnmap[dim].begin();_colcurr!=dnmap[dim].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi=0.0;
            valtxi = - frcoeff * z[dim] * _colcurr->second * ct * jumptxi;
            // do not assemble zeros into matrix
            if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
              double valteta=0.0;
              valteta = - frcoeff * z[dim] * _colcurr->second * ct * jumpteta;
              // do not assemble zeros into matrix
              if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }
        } // loop over all dimensions
      }

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
      if(wearimpl_ and !wearpv_)
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
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

      // prepare assembly, get information from node
      std::vector<GEN::pairedvector<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<GEN::pairedvector<int,double> > dtetamap = cnode->CoData().GetDerivTeta();

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
      std::map<int,double>::iterator           colcurr;
      GEN::pairedvector<int,double>::iterator _colcurr;

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

      // more information from node
      double* txi = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* jump = cnode->FriData().jump();

      // slip
      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];
      }
      else
      {
        for (int i=0;i<Dim();i++)
        {
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
      }

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
      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
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
      }
      else // std slip
      {
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
          for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
          {
            int col = _colcurr->first;
            double val = jump[j]*_colcurr->second;

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
          }

          if(Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
              double val = jump[j]*_colcurr->second;

              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
            }
          }
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
|  Assemble matrix LinSlip with W derivatives               farah 09/13|
*----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleLinSlip_W(LINALG::SparseMatrix& linslipWglobal)
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
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // prepare assembly, get information from node
      std::vector<GEN::pairedvector<int,double> > dnmap = cnode->CoData().GetDerivN();
      std::vector<GEN::pairedvector<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<GEN::pairedvector<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
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
      double* jump = cnode->FriData().jump();

      // for slip
      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];

        for (int i=0;i<Dim();i++)
        {
          znor += n[i]*z[i];
          ztxi += txi[i]*z[i];
          zteta += teta[i]*z[i];
        }
      }
      else
      {
        for (int i=0;i<Dim();i++)
        {
          znor += n[i]*z[i];
          ztxi += txi[i]*z[i];
          zteta += teta[i]*z[i];
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
      }

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
void WEAR::WearInterface::AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
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
  double frcoeff = IParams().get<double>("FRCOEFF");
  double ct = IParams().get<double>("SEMI_SMOOTH_CT");
  double cn = IParams().get<double>("SEMI_SMOOTH_CN");

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
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // prepare assembly, get information from node
      std::vector<GEN::pairedvector<int,double> > dnmap = cnode->CoData().GetDerivN();
      std::vector<GEN::pairedvector<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<GEN::pairedvector<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
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
      std::map<int,double>::iterator           colcurr;
      GEN::pairedvector<int,double>::iterator _colcurr;

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
      double euclidean = 0;
      double* jump = cnode->FriData().jump();

      // for gp slip
      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];

        for (int i=0;i<Dim();i++)
        {
          znor += n[i]*z[i];
          ztxi += txi[i]*z[i];
          zteta += teta[i]*z[i];
        }
      }
      else
      {
        for (int i=0;i<Dim();i++)
        {
          znor += n[i]*z[i];
          ztxi += txi[i]*z[i];
          zteta += teta[i]*z[i];
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
      }

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
          for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
          {
            int col = _colcurr->first;
            double val = (_colcurr->second)*z[j];

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
          }

          if (Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
              double val = (_colcurr->second)*z[j];

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
        std::map<int,double> derivjump1, derivjump2;  // for gp slip
        std::vector<std::map<int,double> > derivjump; // for dm slip

        /*** 01  ********* -Deriv(euclidean).ct.tangent.deriv(u)*ztan ***/
        if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
        {
          derivjump1 = cnode->FriData().GetDerivVarJump()[0];
          derivjump2 = cnode->FriData().GetDerivVarJump()[1];

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
        }
        else //std slip
        {
          // get linearization of jump vector
          derivjump = cnode->FriData().GetDerivJump();

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
        }


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
          if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
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
        if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
        {
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
        }
        else // std. slip
        {
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
        }

        /*** 1 ********************************* euclidean.deriv(T).z ***/
        // loop over dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
          {
            int col = _colcurr->first;
            double val = euclidean*(_colcurr->second)*z[j];

#ifdef CONSISTENTSLIP
            val = val / (znor - cn * wgap);
#endif

            // do not assemble zeros into s matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
          }

          if (Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
              double val = euclidean*(_colcurr->second)*z[j];

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
          for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi = (ztxi+ct*jumptxi)/euclidean*(_colcurr->second)*z[j]*ztxi;
            double valteta = (ztxi+ct*jumptxi)/euclidean*(_colcurr->second)*z[j]*zteta;

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
            for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
              double valtxi = (zteta+ct*jumpteta)/euclidean*(_colcurr->second)*z[j]*ztxi;
              double valteta = (zteta+ct*jumpteta)/euclidean*(_colcurr->second)*z[j]*zteta;

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
        if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
        {
          //!!!!!!!!!!!!!!! DO NOTHING !!!!!!!
        }
        else
        {
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*ct*(_colcurr->second)*jump[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*ct*(_colcurr->second)*jump[j]*zteta;

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
              for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
              {
                int col = _colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*ct*(_colcurr->second)*jump[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*ct*(_colcurr->second)*jump[j]*zteta;

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
        }

        /*** 4 ************************** (frcoeff*znor).deriv(T).z ***/
        // loop over all dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
          {
            int col = _colcurr->first;
#ifdef CONSISTENTSLIP
            double val = - frcoeff * (_colcurr->second)*z[j];
#else
            double val = (-1)*(frcoeff*(znor-cn*wgap))*(_colcurr->second)*z[j];
#endif
            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
          }

          if(Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
#ifdef CONSISTENTSLIP
              double val = - frcoeff * (_colcurr->second) * z[j];
#else
              double val = (-1.0)*(frcoeff*(znor-cn*wgap))*(_colcurr->second)*z[j];
#endif
              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
            }
          }
        }

        /*** 5 *********************** (frcoeff*znor).deriv(T).jump ***/
        if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
        {
          //!!!!!!!!!!!!!!! DO NOTHING !!!!!!!!!!!!!!!!!!!!!!
        }
        else
        {
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
  #ifdef CONSISTENTSLIP
              double val = - frcoeff * ct * (_colcurr->second) * jump[j];
  #else
              double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(_colcurr->second)*jump[j];
  #endif
              // do not assemble zeros into matrix
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
              {
                int col = _colcurr->first;
  #ifdef CONSISTENTSLIP
                double val = - frcoeff * ct * (_colcurr->second) * jump[j];
  #else
                double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(_colcurr->second)*jump[j];
  #endif
                // do not assemble zeros into s matrix
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }
        }

#ifndef CONSISTENTSLIP
        /*** 6 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
        // loop over all dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map
          for (_colcurr=dnmap[j].begin();_colcurr!=dnmap[j].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi = (-1.0)*(ztxi+ct*jumptxi)*frcoeff*(_colcurr->second)*z[j];
            double valteta = (-1.0)*(zteta+ct*jumpteta)*frcoeff*(_colcurr->second)*z[j];

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
        if(wearimpl_ and !wearpv_)
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
    dserror("Tresca friction not implemented for wear !!!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix W_lm containing wear w~ derivatives      farah 07/13|
 |  w.r.t. lm       -- impl wear                                        |
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleLinWLm(LINALG::SparseMatrix& sglobal)
{
  /************************************************
   *  This function is only for Implicit Wear !!! *
   ************************************************/
  if (!wearimpl_)
    dserror("This matrix deriv. is only required for implicit wear algorithm!");

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i) //(int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

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
 |  w.r.t. lm  --> (ONLY !!!) for consistent stick                      |
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleLinWLmSt(LINALG::SparseMatrix& sglobal)
{
  /************************************************
   *  This function is only for Implicit Wear !!! *
   ************************************************/
  if (!wearimpl_)
    dserror("This matrix deriv. is only required for implicit wear algorithm!");

  // get out of here if not participating in interface
  if (!lComm())
    return;

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
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    // prepare assembly, get information from node
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
void WEAR::WearInterface::AssembleLinWLmSl(LINALG::SparseMatrix& sglobal)
{
  /************************************************
   *  This function is only for Implicit Wear !!! *
   ************************************************/

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
  double frcoeff = IParams().get<double>("FRCOEFF");

  // loop over all active slave nodes of the interface
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

    // prepare assembly, get information from node
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
void WEAR::WearInterface::AssembleWear(Epetra_Vector& wglobal)
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
    CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

    if (frinode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleWear: Node ownership inconsistency!");

    /**************************************************** w-vector ******/
    double wear = frinode->WearData().WeightedWear();

    Epetra_SerialDenseVector wnode(1);
    std::vector<int> lm(1);
    std::vector<int> lmowner(1);

    wnode(0) = wear;
    lm[0] = frinode->Id();
    lmowner[0] = frinode->Owner();

    LINALG::Assemble(wglobal,wnode,lm,lmowner);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble Mortar matrice for both sided wear              farah 06/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleD2(LINALG::SparseMatrix& dglobal)
{
  /**************************************************************
   *  This function is only for wear mapping by enforcing weak  *
   *  dirichlet bound. cond on master side !!! *                *
   **************************************************************/

  //*******************************************************
  // assemble second D matrix for both-sided wear :
  // Loop over allreduced map, so that any proc join
  // the loop. Do non-local assembly by FEAssemble to
  // allow owned and ghosted nodes assembly. In integrator
  // only the nodes associated to the owned slave element
  // are filled with entries. Therefore, the integrated
  // entries are unique on nodes!
  //*******************************************************
  const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(mnoderowmap_));

  for (int i=0;i<masternodes->NumMyElements();++i) //mnoderowmap_
  {
    int gid = masternodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    /**************************************************** D2-matrix ******/
    if ((cnode->WearData().GetD2()).size()>0 )
    {
      std::vector<std::map<int,double> > dmap = cnode->WearData().GetD2();
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
          // check for diagonality
          if (row!=col && abs(val)>1.0e-12)
            dserror("ERROR: AssembleDM: D-Matrix is not diagonal!");

          // create an explicitly diagonal d matrix
          if (row==col)
            dglobal.FEAssemble(val, row, col);

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
 |  build active set (nodes / dofs) for Master               farah 11/13|
 *----------------------------------------------------------------------*/
bool WEAR::WearInterface::BuildActiveSetMaster()
{
  //****************************************
  //for both-sided discr wear
  //****************************************

  // spread active and slip information to all procs
  std::vector<int> a;
  std::vector<int> sl;
  // loop over all slave nodes on the current interface
  for (int j=0;j<SlaveRowNodes()->NumMyElements();++j)
  {
    int gid = SlaveRowNodes()->GID(j);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

    if (frinode->Active())
      a.push_back(frinode->Id());

    if (frinode->FriData().Slip())
      sl.push_back(frinode->Id());
  }
  Teuchos::RCP<Epetra_Map> auxa  = Teuchos::rcp(new Epetra_Map(-1,(int)a.size(), &a[0], 0,Comm()));
  Teuchos::RCP<Epetra_Map> auxsl = Teuchos::rcp(new Epetra_Map(-1,(int)sl.size(),&sl[0],0,Comm()));

  const Teuchos::RCP<Epetra_Map> ara  = LINALG::AllreduceEMap(*(auxa));
  const Teuchos::RCP<Epetra_Map> arsl = LINALG::AllreduceEMap(*(auxsl));

  for (int j=0;j<SlaveColNodes()->NumMyElements();++j)
  {
    int gid = SlaveColNodes()->GID(j);

    if(ara->LID(gid)==-1)
      continue;

    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);;
    CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

    if(frinode->Owner()!=Comm().MyPID())
    {
      frinode->Active()=true;
    }
  }
  for (int j=0;j<SlaveColNodes()->NumMyElements();++j)
  {
    int gid = SlaveColNodes()->GID(j);

    if(arsl->LID(gid)==-1)
      continue;

    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);;
    CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

    if(frinode->Owner()!=Comm().MyPID())
      frinode->FriData().Slip()=true;
  }

  // spread info for attached status...
  const Teuchos::RCP<Epetra_Map> meleall = LINALG::AllreduceEMap(*(MasterRowElements()));
  std::vector<int> eleatt;

  for (int j=0;j<meleall->NumMyElements();++j)
  {
    int gid = meleall->GID(j);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarElement* moele = dynamic_cast<MORTAR::MortarElement*>(ele);

    if(moele->IsAttached()==true and moele->Owner()==Comm().MyPID())
      eleatt.push_back(moele->Id());
  }

  Teuchos::RCP<Epetra_Map> auxe = Teuchos::rcp(new Epetra_Map(-1,(int)eleatt.size(),&eleatt[0],0,Comm()));
  const Teuchos::RCP<Epetra_Map> att = LINALG::AllreduceEMap(*(auxe));

  for (int j=0;j<att->NumMyElements();++j)
  {
    int gid = att->GID(j);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find node with gid %",gid);
    MORTAR::MortarElement* moele = dynamic_cast<MORTAR::MortarElement*>(ele);

    if(moele->IsAttached()==false)
      moele->SetAttached()=true;
  }

  //Detect maps
  std::vector<int> wa;
  std::vector<int> wsl;
  std::vector<int> wad;
  std::vector<int> wsln;

  // loop over all slave nodes on the current interface
  for (int j=0;j<SlaveRowNodes()->NumMyElements();++j)
  {
    int gid = SlaveRowNodes()->GID(j);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

    // get elements from node (SLAVE)
    for(int u=0;u<(int)frinode->NumElement();++u)
    {
      // all found MASTER elements:
      for(int k=0;k<(int)dynamic_cast<MORTAR::MortarElement*>(frinode->Elements()[u])->MoData().NumSearchElements();++k)
      {
        int gid2 = dynamic_cast<MORTAR::MortarElement*>(frinode->Elements()[u])->MoData().SearchElements()[k];
        DRT::Element* ele2 = Discret().gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        MORTAR::MortarElement* celement = dynamic_cast<MORTAR::MortarElement*>(ele2);

        // nodes cor. to this master element
        if(celement->IsAttached()==true)
        {
          for(int p=0;p<celement->NumNode();++p)
          {
            CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(celement->Nodes()[p]);

            if(mnode->IsDetected()==false)
            {
              // active master nodes!
              if(frinode->Active())
              {
                wa.push_back(celement->Nodes()[p]->Id());
                wad.push_back(dynamic_cast<MORTAR::MortarNode*>(celement->Nodes()[p])->Dofs()[0]);
              }

              //slip master nodes!
              if(frinode->FriData().Slip())
              {
                wsl.push_back(celement->Nodes()[p]->Id());
                wsln.push_back(dynamic_cast<MORTAR::MortarNode*>(celement->Nodes()[p])->Dofs()[0]);
              }

              //set detection status
//              if(frinode->Active() || frinode->FriData().Slip())
//                mnode->SetDetected()=true;
            }
          }
        }
      }
    }
  }//node loop

  //reset nodes
  // loop over all master nodes on the current interface
  const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(mnoderowmap_));
  const Teuchos::RCP<Epetra_Map> mastereles = LINALG::AllreduceEMap(*(melerowmap_));

  for (int j=0;j<masternodes->NumMyElements();++j)
  {
    int gid = masternodes->GID(j);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(node);

    mnode->SetDetected()=false;
  }
  for (int j=0;j<mastereles->NumMyElements();++j)
  {
    int gid = mastereles->GID(j);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find element with gid %",gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    mele->SetAttached()=false;
  }

  Teuchos::RCP<Epetra_Map> actmn = Teuchos::rcp(new Epetra_Map(-1,(int)wa.size(),&wa[0],0,Comm()));
  Teuchos::RCP<Epetra_Map> slimn = Teuchos::rcp(new Epetra_Map(-1,(int)wsl.size(),&wsl[0],0,Comm()));
  Teuchos::RCP<Epetra_Map> slimd = Teuchos::rcp(new Epetra_Map(-1,(int)wsln.size(),&wsln[0],0,Comm()));

  const Teuchos::RCP<Epetra_Map> ARactmn = LINALG::AllreduceOverlappingEMap(*(actmn));
  const Teuchos::RCP<Epetra_Map> ARslimn = LINALG::AllreduceOverlappingEMap(*(slimn));

  std::vector<int> ga;
  std::vector<int> gs;
  std::vector<int> gsd;

  for (int j=0;j<ARactmn->NumMyElements();++j)
  {
    int gid = ARactmn->GID(j);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(node);

    if(mnode->Owner()==Comm().MyPID())
    {
      bool isin=false;
      for (int k=0;k<(int)ga.size();++k)
      {
        if(ga[k] == mnode->Id())
          isin=true;
      }
      if(!isin)
        ga.push_back(mnode->Id());
    }
  }


  for (int j=0;j<ARslimn->NumMyElements();++j)
  {
    int gid = ARslimn->GID(j);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* mnode = dynamic_cast<CONTACT::FriNode*>(node);

    if(mnode->Owner()==Comm().MyPID())
    {
      bool isin=false;
      for (int k=0;k<(int)gs.size();++k)
      {
        if(gs[k] == mnode->Id())
          isin=true;
      }
      if(!isin)
      {
        gs.push_back(mnode->Id());
        gsd.push_back(mnode->Dofs()[0]);
      }
    }
  }

  activmasternodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)ga.size(),&ga[0],0,Comm()));
  slipmasternodes_  = Teuchos::rcp(new Epetra_Map(-1,(int)gs.size(),&gs[0],0,Comm()));
  slipmn_           = Teuchos::rcp(new Epetra_Map(-1,(int)gsd.size(),&gsd[0],0,Comm()));

  for (int j=0;j<SlaveColNodes()->NumMyElements();++j)
  {
    int gid = SlaveColNodes()->GID(j);

    if(ara->LID(gid)==-1)
      continue;

    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);;
    CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

    if(frinode->Owner()!=Comm().MyPID())
      frinode->Active()=false;
  }
  for (int j=0;j<SlaveColNodes()->NumMyElements();++j)
  {
    int gid = SlaveColNodes()->GID(j);

    if(arsl->LID(gid)==-1)
      continue;

    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);;
    CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

    if(frinode->Owner()!=Comm().MyPID())
      frinode->FriData().Slip()=false;
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                          farah 02/16|
 *----------------------------------------------------------------------*/
bool WEAR::WearInterface::BuildActiveSet(bool init)
{
  // call contact function
  CONTACT::CoInterface::BuildActiveSet(init);

  // define local variables
  std::vector<int> mynodegids(0);
  std::vector<int> mydofgids(0);
  std::vector<int> myslipnodegids(0);
  std::vector<int> myslipdofgids(0);
  std::vector<int> mymnodegids(0);
  std::vector<int> mymdofgids(0);

  // *******************************************************
  // loop over all master nodes - both-sided wear specific
  if (wearboth_ and !wearpv_)
  {
    // The node and dof map of involved master nodes for both-sided wear should have an
    // analog distribution as the master node row map. Therefore, we loop over
    // allreduced map (all procs involved). If a node has the corresponding bool
    // "InvolvedM()" then this information will be communicated to all other procs.
    // Note, this node could be ghosted! At the end, the owning proc will take the
    // information to build the map.

    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(mnoderowmap_));

    for (int k=0;k<masternodes->NumMyElements();++k) //mnoderowmap_
    {
      int gid = masternodes->GID(k);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
      const int numdof = cnode->NumDof();

      int inv = 0;
      if (cnode->InvolvedM()==true)
        inv=1;

      //TODO: not nice... alternative to sumall?
      int invglobal =0;
      Comm().SumAll(&inv,&invglobal,1);
      Comm().Barrier();

      if (cnode->Owner()==Comm().MyPID() && invglobal>0)
      {
        mymnodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
        {
          mymdofgids.push_back(cnode->Dofs()[j]);
        }
      }
      //reset it
      cnode->InvolvedM()=false;
    }

    // create map for all involved master nodes -- both-sided wear specific
    involvednodes_ = Teuchos::rcp(new Epetra_Map(-1,(int)mymnodegids.size(),&mymnodegids[0],0,Comm()));
    involveddofs_  = Teuchos::rcp(new Epetra_Map(-1,(int)mymdofgids.size(),&mymdofgids[0],0,Comm()));
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Initialize Data Container for nodes and elements         farah 02/16|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::InitializeDataContainer()
{
  // call contact class function which calls mortar function
  CONTACT::CoInterface::InitializeDataContainer();

  //***********************************************************
  // both-sided wear
  // here we need a datacontainer for the masternodes too
  // they have to know their involved nodes/dofs
  //***********************************************************
  if (wearboth_)
  {
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(MasterRowNodes()));

    for (int i=0; i<masternodes->NumMyElements(); ++i) //MasterRowNodes()
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %i",gid);
      MORTAR::MortarNode* mnode = dynamic_cast<MORTAR::MortarNode*>(node);

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

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble inactive wear right hand side                   farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleInactiveWearRhs(Epetra_Vector& inactiverhs)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm()) return;

  // FIXME It's possible to improve the performance, if only recently active nodes of the inactive node set,
  // i.e. nodes, which were active in the last iteration, are considered. Since you know, that the lagrange
  // multipliers of former inactive nodes are still equal zero.

  Teuchos::RCP<Epetra_Map> inactivenodes;

  if(sswear_)
    inactivenodes = LINALG::SplitMap(*snoderowmap_, *activenodes_);
  else
    inactivenodes = LINALG::SplitMap(*snoderowmap_, *slipnodes_);

  for (int i=0;i<inactivenodes->NumMyElements();++i)
  {
    int gid = inactivenodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleInactiverhs: Node ownership inconsistency!");

    if (Dim() == 2)
    {
      std::vector<int> w_gid(1);
      std::vector<int> w_owner(1);

      // calculate the tangential rhs
      Epetra_SerialDenseVector w_i(1);

      w_owner[0] = cnode->Owner();
      w_i[0]     = - cnode->WearData().wold()[0] - cnode->WearData().wcurr()[0];    // already negative rhs!!!
      w_gid[0]   = cnode->Dofs()[0];//inactivedofs->GID(2*i);

      if (abs(w_i[0])>1e-15) LINALG::Assemble(inactiverhs, w_i, w_gid, w_owner);
    }
    else if (Dim() == 3)
    {
      std::vector<int> w_gid(1);
      std::vector<int> w_owner(1);

      // calculate the tangential rhs
      Epetra_SerialDenseVector w_i(1);

      w_owner[0] = cnode->Owner();
      w_i[0]     = - cnode->WearData().wold()[0] - cnode->WearData().wcurr()[0];    // already negative rhs!!!
      w_gid[0]   = cnode->Dofs()[0];//inactivedofs->GID(3*i);

      if (abs(w_i[0])>1e-15) LINALG::Assemble(inactiverhs, w_i, w_gid, w_owner);
    }
  }
}


/*----------------------------------------------------------------------*
 |  Assemble inactive wear right hand side                   farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleInactiveWearRhs_Master(Epetra_FEVector& inactiverhs)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm()) return;

  Teuchos::RCP<Epetra_Map> inactivenodes  = LINALG::SplitMap(*mnoderowmap_, *slipmasternodes_);
  Teuchos::RCP<Epetra_Map> inactivedofs  = LINALG::SplitMap(*(MNDofs()), *slipmn_);



  const Teuchos::RCP<Epetra_Map> allredi = LINALG::AllreduceEMap(*(inactivedofs));

  Teuchos::RCP<Epetra_Vector> rhs= LINALG::CreateVector(*allredi, true);

  for (int i=0;i<inactivenodes->NumMyElements();++i)
  {
    int gid = inactivenodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleInactiverhs: Node ownership inconsistency!");

    if (Dim() == 2)
    {
      std::vector<int> w_gid(1);
      std::vector<int> w_owner(1);

      // calculate the tangential rhs
      Epetra_SerialDenseVector w_i(1);

      w_owner[0] = Comm().MyPID();//cnode->Owner();
      w_i[0]     = - cnode->WearData().wold()[0] - cnode->WearData().wcurr()[0];    // already negative rhs!!!
      w_gid[0]   = cnode->Dofs()[0];//inactivedofs->GID(2*i);

      if (abs(w_i[0])>1e-12) LINALG::Assemble(*rhs, w_i, w_gid, w_owner);
    }
    else if (Dim() == 3)
    {
      std::vector<int> w_gid(1);
      std::vector<int> w_owner(1);

      // calculate the tangential rhs
      Epetra_SerialDenseVector w_i(1);

      w_owner[0] = Comm().MyPID();//cnode->Owner();
      w_i[0]     = - cnode->WearData().wold()[0] - cnode->WearData().wcurr()[0];    // already negative rhs!!!
      w_gid[0]   = cnode->Dofs()[0];//inactivedofs->GID(3*i);

      if (abs(w_i[0])>1e-12) LINALG::Assemble(*rhs, w_i, w_gid, w_owner);
    }
  }

  Teuchos::RCP<Epetra_Export> exp = Teuchos::rcp(new Epetra_Export(*allredi,*inactivedofs));
  inactiverhs.Export(*rhs,*exp,Add);


  return;
}


/*----------------------------------------------------------------------*
 |  Assemble wear-cond. right hand side (discr)              farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleWearCondRhs(Epetra_Vector& rhs)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm()) return;

  //nodes for loop
  Teuchos::RCP<Epetra_Map> considerednodes;

  // nothing to do if no active nodes
  if(sswear_)
  {
    if (activenodes_==Teuchos::null)
      return;
    considerednodes = activenodes_;
  }
  else
  {
    if (slipnodes_==Teuchos::null)
      return;
    considerednodes = slipnodes_;
  }

  INPAR::CONTACT::SystemType systype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(IParams(),"SYSTEM");

  double wcoeff = IParams().get<double>("WEARCOEFF");

  typedef std::map<int,double>::const_iterator CI;

  for (int i=0;i<considerednodes->NumMyElements();++i)
  {
    int gid = considerednodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    if (fnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleWearCondRhs: Node ownership inconsistency!");

    /**************************************************** E-matrix ******/
    if ((fnode->WearData().GetE()).size()>0)
    {
      std::map<int,double>  emap = fnode->WearData().GetE()[0];

      for (CI p=emap.begin();p!=emap.end();++p)
      {
        int gid3 = (int)((p->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

        std::vector<int> w_gid(1);
        std::vector<int> w_owner(1);

        Epetra_SerialDenseVector w_i(1);

        w_owner[0] = fnode->Owner();
        w_i[0]     = (- (csnode->WearData().wold()[0]) - (csnode->WearData().wcurr()[0])) * (p->second);
        w_gid[0]   = fnode->Dofs()[0];

        if (abs(w_i[0])>1e-15) LINALG::Assemble(rhs, w_i, w_gid, w_owner);
      }
    }

    /**************************************************** T-matrix ******/
    // for condensation of lm and wear we condense the system with absol. lm
    // --> therefore we do not need the lm^i term...
    if (((fnode->WearData().GetT()).size()>0) && systype != INPAR::CONTACT::system_condensed)
    {
      std::map<int,double> tmap = fnode->WearData().GetT()[0];

      for (CI p=tmap.begin();p!=tmap.end();++p)
      {
        int gid3 = (int)((p->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

        double lmn = 0.0;
        for (int u=0;u<Dim();++u)
          lmn += (csnode->MoData().n()[u]) * (csnode->MoData().lm()[u]);

        std::vector<int> w_gid(1);
        std::vector<int> w_owner(1);

        Epetra_SerialDenseVector w_i(1);

        w_owner[0] = fnode->Owner();
        w_i[0]     =  wcoeff * lmn * (p->second);
        w_gid[0]   = fnode->Dofs()[0];

        if (abs(w_i[0])>1e-15) LINALG::Assemble(rhs, w_i, w_gid, w_owner);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble wear-cond. right hand side (discr)              farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::AssembleWearCondRhs_Master(Epetra_FEVector& RHS)
{
  /************************************************
   *  This function is only for discrete Wear !!! *
   ************************************************/

  // get out of here if not participating in interface
  if (!lComm()) return;

  // nothing to do if no active nodes
  if (slipmasternodes_==Teuchos::null)
    return;

  INPAR::CONTACT::SystemType systype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(IParams(),"SYSTEM");

  double wcoeff = IParams().get<double>("WEARCOEFF_MASTER");

  typedef std::map<int,double>::const_iterator CI;

  const Teuchos::RCP<Epetra_Map> slmasternodes = LINALG::AllreduceEMap(*(slipmasternodes_));
  const Teuchos::RCP<Epetra_Map> slmastern = LINALG::AllreduceEMap(*(slipmn_));

  Teuchos::RCP<Epetra_Vector> rhs= LINALG::CreateVector(*slmastern, true);

  for (int i=0;i<slmasternodes->NumMyElements();++i)
  {
    int gid = slmasternodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::FriNode* fnode = dynamic_cast<CONTACT::FriNode*>(node);

    /**************************************************** E-matrix ******/
    if ((fnode->WearData().GetE()).size()>0)
    {
      std::map<int,double>  emap = fnode->WearData().GetE()[0];

      for (CI p=emap.begin();p!=emap.end();++p)
      {
        int gid3 = (int)((p->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

        std::vector<int> w_gid(1);
        std::vector<int> w_owner(1);

        Epetra_SerialDenseVector w_i(1);

        w_owner[0] = Comm().MyPID();//fnode->Owner();
        w_i[0]     = (- (csnode->WearData().wold()[0]) - (csnode->WearData().wcurr()[0])) * (p->second);
        w_gid[0]   = fnode->Dofs()[0];

        if (abs(w_i[0])>1e-15) LINALG::Assemble(*rhs, w_i, w_gid, w_owner);
      }
    }

    /**************************************************** T-matrix ******/
    // for condensation of lm and wear we condense the system with absol. lm
    // --> therefore we do not need the lm^i term...
    if (((fnode->WearData().GetT()).size()>0) && systype != INPAR::CONTACT::system_condensed)
    {
      std::map<int,double> tmap = fnode->WearData().GetT()[0];

      for (CI p=tmap.begin();p!=tmap.end();++p)
      {
        int gid3 = (int)((p->first)/Dim());
        DRT::Node* snode = idiscret_->gNode(gid3);
        if (!snode) dserror("ERROR: Cannot find node with gid");
        CONTACT::FriNode* csnode = dynamic_cast<CONTACT::FriNode*>(snode);

        double lmn = 0.0;
        for (int u=0;u<Dim();++u)
          lmn += (csnode->MoData().n()[u]) * (csnode->MoData().lm()[u]);

        std::vector<int> w_gid(1);
        std::vector<int> w_owner(1);

        Epetra_SerialDenseVector w_i(1);

        w_owner[0] = Comm().MyPID();//fnode->Owner();
        w_i[0]     = wcoeff * lmn * (p->second);
        w_gid[0]   = fnode->Dofs()[0];

        if (abs(w_i[0])>1e-15) LINALG::Assemble(*rhs, w_i, w_gid, w_owner);
      }
    }
  }

  Teuchos::RCP<Epetra_Export> exp = Teuchos::rcp(new Epetra_Export(*slmastern,*slipmn_));
  RHS.Export(*rhs,*exp,Add);

  return;
}


/*----------------------------------------------------------------------*
 |  initialize / reset interface for wear                    farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)

  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscret_->lColNode(i));

    // reset feasible projection and segmentation status
    node->HasProj()    = false;
    node->HasSegment() = false;
  }

  //**************************************************
  // for both-sided wear
  //**************************************************
  if (wearboth_ and !wearpv_)
  {
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(MasterRowNodes()));

    for (int i=0;i<masternodes->NumMyElements();++i) //for (int i=0;i<MasterRowNodes()->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);

      if (cnode->IsSlave() == false)
      {
        // reset nodal Mortar maps
        for (int j=0;j<(int)((cnode->WearData().GetD2()).size());++j)
          (cnode->WearData().GetD2())[j].clear();

        (cnode->WearData().GetD2()).resize(0);
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
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    // reset nodal Mortar maps
    cnode->MoData().GetD().clear();
    cnode->MoData().GetM().clear();
    cnode->MoData().GetMmod().clear();

    // reset nodal scaling factor
    (cnode->MoData().GetScale())=0.0;
    (cnode->CoData().GetDerivScale()).clear();

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((cnode->CoData().GetDerivN()).size());++j)
      (cnode->CoData().GetDerivN())[j].clear();
    (cnode->CoData().GetDerivN()).resize(0,0);

    // reset derivative maps of tangent vectors
    for (int j=0;j<(int)((cnode->CoData().GetDerivTxi()).size());++j)
      (cnode->CoData().GetDerivTxi())[j].clear();
    (cnode->CoData().GetDerivTxi()).resize(0,0);
    for (int j=0;j<(int)((cnode->CoData().GetDerivTeta()).size());++j)
      (cnode->CoData().GetDerivTeta())[j].clear();
    (cnode->CoData().GetDerivTeta()).resize(0,0);

    // reset derivative map of Mortar matrices
    (cnode->CoData().GetDerivD()).clear();
    (cnode->CoData().GetDerivM()).clear();

    // reset nodal weighted gap and derivative
    cnode->CoData().Getg() = 1.0e12;
    (cnode->CoData().GetDerivG()).clear();

    // reset derivative map of lagrange multipliers
    for (int j=0; j<(int)((cnode->CoData().GetDerivZ()).size()); ++j)
      (cnode->CoData().GetDerivZ())[j].clear();
    (cnode->CoData().GetDerivZ()).resize(0);

    //************************************
    //              friction
    //*************************************
    if (friction_)
    {
      CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(cnode);

      // reset SNodes and Mnodes
      frinode->FriData().GetSNodes().clear();
      frinode->FriData().GetMNodes().clear();

      // for gp slip
      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
        // reset jump deriv.
        for (int j=0 ;j<(int)((frinode->FriData().GetDerivVarJump()).size()); ++j )
          (frinode->FriData().GetDerivVarJump())[j].clear();

        (frinode->FriData().GetDerivVarJump()).resize(2);

        // reset jumps
        frinode->FriData().jump_var()[0]=0.0;
        frinode->FriData().jump_var()[1]=0.0;
      }

      //************************************
      //              wear
      //*************************************
      // reset weighted wear increment and derivative
      // only for implicit wear algorithm
      if(wearimpl_ and !wearpv_)
      {
        (cnode->CoData().GetDerivW()).clear();
        (cnode->CoData().GetDerivWlm()).clear();
      }

      if(wearpv_)
      {
        // reset nodal Mortar wear maps
        for (int j=0;j<(int)((frinode->WearData().GetT()).size());++j)
          (frinode->WearData().GetT())[j].clear();
        for (int j=0;j<(int)((frinode->WearData().GetE()).size());++j)
          (frinode->WearData().GetE())[j].clear();

        (frinode->WearData().GetT()).resize(0);
        (frinode->WearData().GetE()).resize(0);

        (frinode->WearData().GetDerivTw()).clear();
        (frinode->WearData().GetDerivE()).clear();

        (frinode->CoData().GetDerivGW()).clear();
      }
      if (wear_)
      {
        // reset wear increment
        if(!wearpv_)
          frinode->WearData().DeltaWeightedWear() = 0.0;

        // reset abs. wear.
        // for impl. wear algor. the abs. wear equals the
        // delta-wear
        if(wearimpl_ and !wearpv_)
          frinode->WearData().WeightedWear() = 0.0;
      }
    }
  }

  // for both-sided wear with discrete wear
  if (wearboth_ and wearpv_)
  {
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(MasterRowNodes()));

    for (int i=0;i<masternodes->NumMyElements();++i) //for (int i=0;i<MasterRowNodes()->NumMyElements();++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(node);

      // reset nodal Mortar wear maps
      for (int j=0;j<(int)((frinode->WearData().GetT()).size());++j)
        (frinode->WearData().GetT())[j].clear();
      for (int j=0;j<(int)((frinode->WearData().GetE()).size());++j)
        (frinode->WearData().GetE())[j].clear();

      (frinode->WearData().GetT()).resize(0);
      (frinode->WearData().GetE()).resize(0);

      (frinode->WearData().GetDerivTw()).clear();
      (frinode->WearData().GetDerivE()).clear();
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
      MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

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
      MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

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
void WEAR::WearInterface::SplitSlaveDofs()
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
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

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
 |  create mnode n_                                          farah 11/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::SplitMasterDofs()
{
  // get out of here if active set is empty
  if (mnoderowmap_==Teuchos::null)
  {
    mndofmap_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return;
  }

  else if (mnoderowmap_->NumGlobalElements()==0)
  {
    mndofmap_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return;
  }

  // define local variables
  int countN=0;
  std::vector<int> myNgids(mnoderowmap_->NumMyElements());

  // dimension check
  double dimcheck =(mdofrowmap_->NumGlobalElements())/(mnoderowmap_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitMasterDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slave nodes
  for (int i=0;i<mnoderowmap_->NumMyElements();++i)
  {
    int gid = mnoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

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
  if ((gcountN)!=mnoderowmap_->NumGlobalElements())
    dserror("ERROR: SplitSlaveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  mndofmap_ = Teuchos::rcp(new Epetra_Map(gcountN,countN,&myNgids[0],0,Comm()));

  return;
}


/*----------------------------------------------------------------------*
 |  compute element areas (public)                           farah 02/16|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::SetElementAreas()
{
  //**********************************************************************
  // In general, it is sufficient to compute element areas only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need the element areas of all elements
  // (slave and master) in the fully overlapping column map there. At the
  // same time we initialize the element data containers for self contact.
  // This is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact() or wearboth_)
  {
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      MORTAR::MortarElement* element = dynamic_cast<MORTAR::MortarElement*>(idiscret_->lColElement(i));
      element->InitializeDataContainer();
      element->MoData().Area()=element->ComputeArea();
    }
  }
  else
  {
    //refer call back to base class version
    MORTAR::MortarInterface::SetElementAreas();
  }

  return;
}


/*----------------------------------------------------------------------*
 |  update wear set (dofs)                                   farah 09/13|
 *----------------------------------------------------------------------*/
void WEAR::WearInterface::UpdateWSets(int offset_if, int maxdofwear, bool bothdiscr)
{
  //********************************************************************
  // WEAR DOFS --  one per node
  //********************************************************************

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

  //********************************************************************
  // For discrete both-sided wear
  //********************************************************************
  if(bothdiscr)
  {
    // temporary vector of W dofs
    std::vector<int> wmdof;

    maxdofwear+=wdofmap_->NumGlobalElements();

    // gather information over all procs
    std::vector<int> localnumwdof(Comm().NumProc());
    std::vector<int> globalnumlmdof(Comm().NumProc());
    localnumwdof[Comm().MyPID()] = (int)( (mdofrowmap_->NumMyElements())/Dim() );
    Comm().SumAll(&localnumwdof[0],&globalnumlmdof[0],Comm().NumProc());

    // compute offet for LM dof initialization for all procs
    int offset = 0;
    for (int k=0;k<Comm().MyPID();++k)
      offset += globalnumlmdof[k];

    // loop over all slave dofs and initialize LM dofs
    for (int i=0; i<(int)( (mdofrowmap_->NumMyElements())/Dim() ); ++i)
      wmdof.push_back(maxdofwear + 1 + offset_if + offset + i);

    // create interface w map
    // (if maxdofglobal_ == 0, we do not want / need this)
    if (maxdofwear>0)
      wmdofmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)wmdof.size(),&wmdof[0],0,Comm()));
  }

  return;
}
