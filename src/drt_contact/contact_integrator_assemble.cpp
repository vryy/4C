/*!----------------------------------------------------------------------
\file contact_integrator_tools.cpp
\brief ...

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
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_integrator.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_defines.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

#include "../drt_inpar/inpar_contact.H"
#include "contact_manager.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"

/*----------------------------------------------------------------------*
 | Assemble third mortar matrix D2                           farah 06/13|
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleD2(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& mele,
                                    Epetra_SerialDenseMatrix& d2seg)
{
  // get adjacent nodes to assemble to
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleD2: Null pointer for mnodes!");

  if (mele.IsSlave()==true )
    dserror("This function is only for master elements");

  // loop over all master nodes
  for (int master=0;master<mele.NumNode();++master)
  {
    CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);

    // only process slave node rows that belong to this proc
    //std::cout << "owner= " << mnode->Owner() << "  current proc= " << comm.MyPID()<< std::endl;
//    if (mnode->Owner() != comm.MyPID())
//      continue;

    int mndof = mnode->NumDof();

    // loop over all dofs of the slave node
    for (int mdof=0;mdof<mndof;++mdof)
    {
      // loop over all master nodes again ("master2 nodes")
      for (int master2=0;master2<mele.NumNode();++master2)
      {
        CONTACT::CoNode* m2node = static_cast<CONTACT::CoNode*>(mnodes[master2]);
        const int* m2dofs = m2node->Dofs();
        int m2ndof = m2node->NumDof();

        // loop over all dofs of the slave node again ("master dofs")
        for (int m2dof=0;m2dof<m2ndof;++m2dof)
        {
          int col = m2dofs[m2dof];
          double val = d2seg(master*mndof+mdof,master2*m2ndof+m2dof);

          if(abs(val)>1e-12)
          {
            mnode->AddD2Value(mdof,col,val);
            //set this master node as "involved" for both sided wear
            mnode->InvolvedM()=true;
            //std::cout << "node= " << mnode->Id()<< std::endl;
            //std::cout << "DOF FROM NODE " << *mnode->Dofs() << std::endl;
          }
        }
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave             |
 |  element to the D map of the adjacent slave nodes.                   |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleD(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    Epetra_SerialDenseMatrix& dseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleD: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
    int sndof = snode->NumDof();

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all slave nodes again ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(snodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the slave node again ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = dseg(slave*sndof+sdof,master*mndof+mdof);

          // BOUNDARY NODE MODIFICATION **********************************
          // We have modified their neighbors' dual shape functions, so we
          // now have a problem with off-diagonal entries occuring in D.
          // Of course we want to keep the diagonality property of the D
          // matrix, but still we may not modify the whole Mortar coupling
          // setting! We achieve both by appling a quite simple but very
          // effective trick: The boundary nodes have already been defined
          // as being master nodes, so all we have to do here, is to shift
          // the off-diagonal terms from D to the resepective place in M,
          // which is not diagonal anyway! (Mind the MINUS sign!!!)
          // *************************************************************
          if (mnode->IsOnBound())
          {
            double minusval = -val;
            if(abs(val)>1e-12) snode->AddMValue(sdof,col,minusval);
            if(abs(val)>1e-12) snode->AddMNode(mnode->Id()); // only for friction!
          }
          else
          {
            if(abs(val)>1e-12) snode->AddDValue(sdof,col,val);
            if(abs(val)>1e-12) snode->AddSNode(mnode->Id()); // only for friction!
          }
        }
      }
    }
    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << std::endl;
    std::map<int, double> nodemap0 = (snode->GetD())[0];
    std::map<int, double> nodemap1 = (snode->GetD())[1];
    typedef std::map<int,double>::const_iterator CI;

    std::cout << "Row dof id: " << sdofs[0] << std::endl;;
    for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
      std::cout << p->first << '\t' << p->second << std::endl;

    std::cout << "Row dof id: " << sdofs[1] << std::endl;
    for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
      std::cout << p->first << '\t' << p->second << std::endl;
#endif // #ifdef DEBUG
    */
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble T contribution (2D / 3D)                        farah 09/13|
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleT(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    Epetra_SerialDenseMatrix& tseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleT: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::FriNode* snode = static_cast<CONTACT::FriNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // loop over all slave nodes again ("master nodes")
    for (int master=0;master<sele.NumNode();++master)
    {
      CONTACT::FriNode* mnode = static_cast<CONTACT::FriNode*>(snodes[master]);
      const int* mdofs = mnode->Dofs();

      int col = mdofs[0];
      int row = 0;
      double val = tseg(slave,master);

      if(abs(val)>1e-12) snode->AddTValue(row,col,val);
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble E contribution (2D / 3D)                        farah 09/13|
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleE(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    Epetra_SerialDenseMatrix& eseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleE: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::FriNode* snode = static_cast<CONTACT::FriNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all slave nodes again ("master nodes")
    for (int master=0;master<sele.NumNode();++master)
    {
      CONTACT::FriNode* mnode = static_cast<CONTACT::FriNode*>(snodes[master]);
      const int* mdofs = mnode->Dofs();

      int row = 0;
      int col = mdofs[0];
      double val = eseg(slave,master);

      if(abs(val)>1e-12) snode->AddEValue(row,col,val);
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution (2D / 3D)                         popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleD(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    MORTAR::IntElement& sintele,
                                    Epetra_SerialDenseMatrix& dseg)
{
  // get adjacent int nodes to assemble to
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("ERROR: AssembleD: Null pointer for sintnodes!");
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleD: Null pointer for snodes!");

  // loop over all slave int nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    CONTACT::CoNode* sintnode = static_cast<CONTACT::CoNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all slave nodes ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(snodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the slave node ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = dseg(slave*sintndof+sdof,master*mndof+mdof);

          // assembly
          if(abs(val)>1e-12) sintnode->AddDValue(sdof,col,val);
          if(abs(val)>1e-12) sintnode->AddSNode(mnode->Id()); // only for friction!
        }
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave and master  |
 |  overlap pair to the M map of the adjacent slave nodes.              |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleM(const Epetra_Comm& comm,
                                    MORTAR::MortarElement& sele,
                                    MORTAR::MortarElement& mele,
                                    Epetra_SerialDenseMatrix& mseg)
{
  // get adjacent slave nodes and master nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleM: Null pointer for snodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleM: Null pointer for mnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
    int sndof = snode->NumDof();

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mseg(slave*sndof+sdof,master*mndof+mdof);
          if(abs(val)>1e-12) snode->AddMValue(sdof,col,val);
          if(abs(val)>1e-12) snode->AddMNode(mnode->Id());  // only for friction!
        }
      }
    }
    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << std::endl;
    std::map<int, double> nodemap0 = (snode->GetM())[0];
    std::map<int, double> nodemap1 = (snode->GetM())[1];
    typedef std::map<int,double>::const_iterator CI;

    std::cout << "Row dof id: " << sdofs[0] << std::endl;;
    for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
      std::cout << p->first << '\t' << p->second << std::endl;

    std::cout << "Row dof id: " << sdofs[1] << std::endl;
    for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
      std::cout << p->first << '\t' << p->second << std::endl;
#endif // #ifdef DEBUG
     */
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Assemble M contribution for eleb. int.(2D / 3D)          farah 01/13|
 |  Assemblle also the contribution for friction                        |
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleM_EleBased(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                         std::vector<MORTAR::MortarElement*> meles,
                                         Epetra_SerialDenseMatrix& mseg)
{
  //problem dimension via numdof of an arbitrary slave node
  DRT::Node** testnodes = sele.Nodes();
  CONTACT::CoNode* tnode = static_cast<CONTACT::CoNode*>(testnodes[0]);
  int dim = tnode->NumDof();

  for ( int nummasterele=0;nummasterele<(int)meles.size();++nummasterele )
  {
    //if all entries zero-valued --> no assembly
    bool nozero=false;
    for ( int i=0;i<dim*meles[nummasterele]->NumNode();++i  )
    {
      for (int j=0;j<dim*sele.NumNode();++j)
      {
        if(mseg(j,i+nummasterele*meles[nummasterele]->NumNode()*dim) != 0.0 )
        {
          nozero=true;
          break;
        }
      }
    }

    if (nozero==true)
    {
      // get adjacent slave nodes and master nodes
      DRT::Node** snodes = sele.Nodes();
      if (!snodes)
      dserror("ERROR: AssembleM_EleBased: Null pointer for snodes!");
      DRT::Node** mnodes = meles[nummasterele]->Nodes();
      if (!mnodes)
      dserror("ERROR: AssembleM_EleBased: Null pointer for mnodes!");

      // loop over all slave nodes
      for (int slave=0;slave<sele.NumNode();++slave)
      {
        CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
        int sndof = snode->NumDof();

        // only process slave node rows that belong to this proc
        if (snode->Owner() != comm.MyPID())
          continue;

        // do not process slave side boundary nodes
        // (their row entries would be zero anyway!)
        if (snode->IsOnBound())
          continue;

        // loop over all dofs of the slave node
        for (int sdof=0;sdof<sndof;++sdof)
        {
          // loop over all master nodes
          for (int master=0;master<meles[nummasterele]->NumNode();++master)
          {
            CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
            const int* mdofs = mnode->Dofs(); //global dofs
            int mndof = mnode->NumDof();

            // loop over all dofs of the master node
            for (int mdof=0;mdof<mndof;++mdof)
            {
              int col = mdofs[mdof];

              double val = mseg(slave*sndof+sdof,master*mndof+mdof+nummasterele*mndof*meles[nummasterele]->NumNode());
              //snode->AddMValue(sdof,col,val);

              if(abs(val)>1e-12) snode->AddMValue(sdof,col,val);
              if(abs(val)>1e-12) snode->AddMNode(mnode->Id());  // only for friction!
            }
          }
        }
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution (2D / 3D)                         popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleM(const Epetra_Comm& comm,
                                    MORTAR::IntElement& sintele,
                                    MORTAR::MortarElement& mele,
                                    Epetra_SerialDenseMatrix& mseg)
{
  // get adjacent slave int nodes and master nodes
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("ERROR: AssembleM: Null pointer for sintnodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: AssembleM: Null pointer for mnodes!");

  // loop over all slave int nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    CONTACT::CoNode* sintnode = static_cast<CONTACT::CoNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mseg(slave*sintndof+sdof,master*mndof+mdof);
          if(abs(val)>1e-12) sintnode->AddMValue(sdof,col,val);
          if(abs(val)>1e-12) sintnode->AddMNode(mnode->Id()); // only for friction!
        }
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 01/08|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the weighted gap of the adjacent slave nodes.       |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleG(const Epetra_Comm& comm,
                                      MORTAR::MortarElement& sele,
                                      Epetra_SerialDenseVector& gseg)
{
  // get adjacent slave nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleG: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = gseg(slave);
    snode->AddgValue(val);

    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << std::endl;
    std::cout << "Weighted gap: " << snode->Getg() << std::endl;
#endif // #ifdef DEBUG
    */
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble u~ contribution (2D / 3D)                       farah 08/13|
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleU(const Epetra_Comm& comm,
                                      MORTAR::MortarElement& sele,
                                      Epetra_SerialDenseMatrix& useg)
{
  // get adjacent slave nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleG: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::FriNode* snode = static_cast<CONTACT::FriNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = 0.0;

    val=useg(slave,0);
    snode->AddJumpValue(val,0);

    if (Dim()==3)
    {
      val=useg(slave,1);
      snode->AddJumpValue(val,1);
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleG(const Epetra_Comm& comm,
                                      MORTAR::IntElement& sintele,
                                      Epetra_SerialDenseVector& gseg)
{
  // get adjacent slave int nodes to assemble to
  DRT::Node** snodes = sintele.Nodes();
  if (!snodes) dserror("ERROR: AssembleG: Null pointer for sintnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = gseg(slave);
    snode->AddgValue(val);
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble  mechanical dissipation to slave nodes       gitterle 08/10|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the mechanical dissipation of adjacent nodes        |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleMechDissSlave(const Epetra_Comm& comm,
                                                  MORTAR::MortarElement& sele,
                                                  Epetra_SerialDenseVector& mdissseg)
{
  // get adjacent slave node to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleMechDissSlave: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::FriNode* snode = static_cast<CONTACT::FriNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = mdissseg(slave);
    snode->AddMechDissValue(val);
   }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble  mechanical dissipation to master nodes      gitterle 10/10|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the mechanical dissipation of adjacent nodes        |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleMechDissMaster(const Epetra_Comm& comm,
                                                   MORTAR::MortarElement& mele,
                                                   Epetra_SerialDenseVector& mdissseg)
{
  // get adjacent master node to assemble to
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: AssembleMechDissMaster: Null pointer for snodes!");

  // loop over all master nodes
  for (int master=0;master<mele.NumNode();++master)
  {
    CONTACT::FriNode* mnode = static_cast<CONTACT::FriNode*>(mnodes[master]);

    double val = mdissseg(master);
    mnode->AddMechDissValue(val);
   }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble A contribution (2D / 3D)                     gitterle 10/10|
 |  This method assembles the contrubution of a 1D/2D slave             |
 |  element to the A map of the adjacent slave nodes.                   |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleA(const Epetra_Comm& comm,
                                     MORTAR::MortarElement& sele,
                                     Epetra_SerialDenseMatrix& aseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleA: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(snodes[slave]);
    int sndof = snode->NumDof();

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all slave nodes again ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(snodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the slave node again ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = aseg(slave*sndof+sdof,master*mndof+mdof);

          if(abs(val)>1e-12) snode->AddAValue(sdof,col,val);
          if(abs(val)>1e-12) snode->AddANode(mnode->Id());
        }
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble B contribution (2D / 3D)                     gitterle 10/10|
 |  This method assembles the contribution of a 1D/2D master            |
 |  element to the B map of the adjacent master node.                   |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleB(const Epetra_Comm& comm,
                                     MORTAR::MortarElement& mele,
                                     Epetra_SerialDenseMatrix& bseg)
{
  // get adjacent nodes to assemble to
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleB: Null pointer for mnodes!");

  // loop over all master nodes
  for (int slave=0;slave<mele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = static_cast<CONTACT::CoNode*>(mnodes[slave]);
    int sndof = snode->NumDof();

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all master nodes again ("master nodes")
      for (int master=0;master<mele.NumNode();++master)
      {
        CONTACT::CoNode* mnode = static_cast<CONTACT::CoNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node again ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = bseg(slave*sndof+sdof,master*mndof+mdof);

          if(abs(val)>1e-12) snode->AddBValue(sdof,col,val);
          if(abs(val)>1e-12) snode->AddBNode(mnode->Id());
        }
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
  |  Assemble wear contribution (2D / 3D)                  gitterle 12/10|
  |  This method assembles the contribution of a 1D/2D slave and master  |
  |  overlap pair to the wear of the adjacent slave nodes.               |
  *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleWear(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                         Epetra_SerialDenseVector& wseg)
{
  // get adjacent slave nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleWear: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::FriNode* snode = static_cast<CONTACT::FriNode*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = wseg(slave);
    snode->AddDeltaWearValue(val);
  }

  return true;
}

