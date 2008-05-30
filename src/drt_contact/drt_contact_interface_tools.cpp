/*!----------------------------------------------------------------------
\file drt_contact_interface_tools.cpp
\brief

<pre>
Maintainer: Alexander Popp
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_interface.H"
#include "drt_cdofset.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io_gmsh.H"
#include "drt_contact_projector.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"
#include "../drt_lib/drt_globalproblem.H"

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |  Visualize node projections with gmsh                      popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::VisualizeGmsh(const Epetra_SerialDenseMatrix& csegs)
{
  // increase counter variable by one
  counter_+=1;

  // construct unique filename for gmsh output
  std::ostringstream filename;
  filename << allfiles.outputfile_kenner << "_";
  if (counter_<10)
    filename << 0 << 0 << 0 << 0;
  else if (counter_<100)
    filename << 0 << 0 << 0;
  else if (counter_<1000)
    filename << 0 << 0;
  else if (counter_<10000)
    filename << 0;

  filename << counter_ << ".pos";

  // do output to file in c-style
  FILE* fp = NULL;

  for (int proc=0;proc<comm_.NumProc();++proc)
  {
    if (proc==comm_.MyPID())
    {
      // open file (overwrite if proc==0, else append)
      if (proc==0) fp = fopen(filename.str().c_str(), "w");
      else fp = fopen(filename.str().c_str(), "a");

      // write output to temporary stringstream
      std::stringstream gmshfilecontent;
      if (proc==0) gmshfilecontent << "View \" Slave and master side CElements \" {" << endl;

      // plot elements
      for (int i=0; i<idiscret_->NumMyRowElements(); ++i)
      {
        CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lRowElement(i));
        int nnodes = element->NumNode();
        LINALG::SerialDenseMatrix coord(3,nnodes);
        double area = element->Area();

        // 2D linear case (2noded line elements)
        if (element->Shape()==DRT::Element::line2)
        {
          coord = element->GetNodalCoords();

          gmshfilecontent << "SL(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << ")";
          gmshfilecontent << "{" << scientific << area << "," << area << "};" << endl;
        }

        // 2D quadratic case (3noded line elements)
        if (element->Shape()==DRT::Element::line3)
        {
          coord = element->GetNodalCoords();

          gmshfilecontent << "SL2(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                              << coord(2,2) << ")";
          gmshfilecontent << "{" << scientific << area << "," << area << "," << area << "};" << endl;
        }
      }

      // plot normal vectors
      for (int i=0; i<snodecolmap_->NumMyElements(); ++i)
      {
        int gid = snodecolmap_->GID(i);
        DRT::Node* node = idiscret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cnode = static_cast<CNode*>(node);
        if (!cnode) dserror("ERROR: Static Cast to CNode* failed");

         double nc[3];
         double nn[3];

         for (int j=0;j<3;++j)
         {
           nc[j]=cnode->xspatial()[j];
           nn[j]=3*cnode->n()[j]; // 2.9 because of gmsh color bar (3 procs(
         }

         gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
         gmshfilecontent << "{" << scientific << nn[0] << "," << nn[1] << "," << nn[2] << "};" << endl;
      }

      // plot contact segments (slave and master projections)
      if (csegs.M()!=0)
      {
        for (int i=0; i<csegs.M(); ++i)
        {
          gmshfilecontent << "SQ(" << scientific << csegs(i,0) << "," << csegs(i,1) << ","
                                   << csegs(i,2) << "," << csegs(i,3) << "," << csegs(i,4) << ","
                                   << csegs(i,5) << "," << csegs(i,6) << "," << csegs(i,7) << ","
                                   << csegs(i,8) << "," << csegs(i,9) << "," << csegs(i,10) << ","
                                   << csegs(i,11) << ")";
          gmshfilecontent << "{" << scientific << proc << "," << proc << "," << proc << "," << proc << "};" << endl;

          gmshfilecontent << "SL(" << scientific << csegs(i,0) << "," << csegs(i,1) << ","
                          << csegs(i,2) << "," << csegs(i,3) << "," << csegs(i,4) << ","
                          << csegs(i,5) << ")";
          gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;

          gmshfilecontent << "SL(" << scientific << csegs(i,6) << "," << csegs(i,7) << ","
                          << csegs(i,8) << "," << csegs(i,9) << "," << csegs(i,10) << ","
                          << csegs(i,11) << ")";
          gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
        }
      }

      if (proc==comm_.NumProc()-1) gmshfilecontent << "};" << endl;

      // move everything to gmsh post-processing file and close it
      fprintf(fp,gmshfilecontent.str().c_str());
      fclose(fp);
    }
    comm_.Barrier();
  }

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for normal/tangent deriv.          popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckNormalDeriv()
{
  // global loop to apply FD scheme to all slave dofs (=2*nodes)
  for (int i=0; i<2*snodefullmap_->NumMyElements();++i)
  {
    // reset normal etc.
    Initialize();
    
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int j=0;j<idiscret_->NumMyColElements();++j)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(j));
      element->Area()=element->ComputeArea();
      //cout << "Element: " << element->Id() << " Nodes: " << (element->Nodes()[0])->Id() << " "
      //     << (element->Nodes()[1])->Id() << " Area: " << element->Area() << endl;
    }  
    
    // create storage for normals / tangents
    vector<double> refnx(int(snodecolmapbound_->NumMyElements()));
    vector<double> refny(int(snodecolmapbound_->NumMyElements()));
    vector<double> newnx(int(snodecolmapbound_->NumMyElements()));
    vector<double> newny(int(snodecolmapbound_->NumMyElements()));
    
    vector<double> reftx(int(snodecolmapbound_->NumMyElements()));
    vector<double> refty(int(snodecolmapbound_->NumMyElements()));
    vector<double> newtx(int(snodecolmapbound_->NumMyElements()));
    vector<double> newty(int(snodecolmapbound_->NumMyElements()));
        
    // compute and print all nodal normals / derivatives (reference)
    for(int j=0; j<snodecolmapbound_->NumMyElements();++j)
    {
      int jgid = snodecolmapbound_->GID(j);
      DRT::Node* jnode = idiscret_->gNode(jgid);
      if (!jnode) dserror("ERROR: Cannot find node with gid %",jgid);
      CNode* jcnode = static_cast<CNode*>(jnode);

      // build averaged normal at each slave node
      jcnode->BuildAveragedNormal();
      
      typedef map<int,double>::const_iterator CI;
      
      // print reference data only once
      if (i==0)
      {
        cout << endl << "Node: " << jcnode->Id() << "  Owner: " << jcnode->Owner() << endl;
        cout << "Normal-derivative-maps: " << endl;
        cout << "Row dof id: " << jcnode->Dofs()[0] << endl;
        for (CI p=(jcnode->GetDerivN()[0]).begin();p!=(jcnode->GetDerivN()[0]).end();++p)
          cout << p->first << '\t' << p->second << endl;    
        cout << "Row dof id: " << jcnode->Dofs()[1] << endl;
        for (CI p=(jcnode->GetDerivN()[1]).begin();p!=(jcnode->GetDerivN()[1]).end();++p)
          cout << p->first << '\t' << p->second << endl;
        cout << "Tangent-derivative-maps: " << endl;
        cout << "Row dof id: " << jcnode->Dofs()[0] << endl;
        for (CI p=(jcnode->GetDerivT()[0]).begin();p!=(jcnode->GetDerivT()[0]).end();++p)
          cout << p->first << '\t' << p->second << endl;    
        cout << "Row dof id: " << jcnode->Dofs()[1] << endl;
        for (CI p=(jcnode->GetDerivT()[1]).begin();p!=(jcnode->GetDerivT()[1]).end();++p)
          cout << p->first << '\t' << p->second << endl;
      }
      
      // store reference normals / tangents
      refnx[j] = jcnode->n()[0];
      refny[j] = jcnode->n()[1];
      reftx[j] =-jcnode->n()[1];
      refty[j] = jcnode->n()[0];
    }
    
    // now fincally get the node we want to apply the FD scheme to
    int gid = snodefullmap_->GID(i/2);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);
    
    // apply finite difference scheme
    cout << "\nBuilding FD for Slave Node: " << snode->Id() << " Dof(l): " << i%2
         << " Dof(g): " << snode->Dofs()[i%2] << endl;
    
    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (i%2==0)
    {
      snode->xspatial()[0] += delta;
      snode->u()[0] += delta;
    }
    else
    {
      snode->xspatial()[1] += delta;
      snode->u()[1] += delta;
    }
  
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int j=0;j<idiscret_->NumMyColElements();++j)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(j));
      element->Area()=element->ComputeArea();
      //cout << "Element: " << element->Id() << " Nodes: " << (element->Nodes()[0])->Id() << " "
      //     << (element->Nodes()[1])->Id() << " Area: " << element->Area() << endl;
    }           
        
    // compute finite difference derivative
    for(int k=0; k<snodecolmapbound_->NumMyElements();++k)
    {
      int kgid = snodecolmapbound_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %",kgid);
      CNode* kcnode = static_cast<CNode*>(knode);
      
      // build NEW averaged normal at each slave node
      kcnode->BuildAveragedNormal();
            
      newnx[k] = kcnode->n()[0];
      newny[k] = kcnode->n()[1];
      newtx[k] =-kcnode->n()[1];
      newty[k] = kcnode->n()[0];
      
      // get reference normal / tangent
      double refn[2] = {0.0, 0.0};
      double reft[2] = {0.0, 0.0};
      refn[0] = refnx[k];
      refn[1] = refny[k];
      reft[0] = reftx[k];
      reft[1] = refty[k];
      
      // get modified normal / tangent
      double newn[2] = {0.0, 0.0};
      double newt[2] = {0.0, 0.0};
      newn[0] = newnx[k];
      newn[1] = newny[k];
      newt[0] = newtx[k];
      newt[1] = newty[k];
      
      // print results (derivatives) to screen
      if (abs(newn[0]-refn[0])>1e-12 || abs(newn[1]-refn[1])>1e-12)
      {
        cout << "Node: " << kcnode->Id() << "  Owner: " << kcnode->Owner() << endl;
        cout << "Normal derivative (FD):" << endl;
        if (abs(newn[0]-refn[0])>1e-12)
        { 
          double val = (newn[0]-refn[0])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[0] << endl;
          cout << snode->Dofs()[i%2] << '\t' << val << endl;
        }
      
        if (abs(newn[1]-refn[1])>1e-12)
        {
          double val = (newn[1]-refn[1])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[1] << endl;
          cout << snode->Dofs()[i%2] << '\t' << val << endl;
        }
      }
      
      if (abs(newt[0]-reft[0])>1e-12 || abs(newt[1]-reft[1])>1e-12)
      {
        cout << "Node: " << kcnode->Id() << "  Owner: " << kcnode->Owner() << endl;
        cout << "Tangent derivative (FD):" << endl;
        if (abs(newt[0]-reft[0])>1e-12)
        { 
          double val = (newt[0]-reft[0])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[0] << endl;
          cout << snode->Dofs()[i%2] << '\t' << val << endl;
        }
      
        if (abs(newt[1]-reft[1])>1e-12)
        {
          double val = (newt[1]-reft[1])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[1] << endl;
          cout << snode->Dofs()[i%2] << '\t' << val << endl;
        }
      }
      
    }
    
    // undo finite difference modification
    if (i%2==0)
    {
      snode->xspatial()[0] -= delta;
      snode->u()[0] -= delta;
    }
    else
    {
      snode->xspatial()[1] -= delta;
      snode->u()[1] -= delta;
    }
  }
 
  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for D-Mortar derivatives           popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckMortarDDeriv()
{
  // create storage for D-Matrix entries
  vector<double> refD(int(snoderowmap_->NumMyElements()));
  vector<double> newD(int(snoderowmap_->NumMyElements()));
      
  // print reference to screen (D-derivative-maps)
  // loop over proc's slave nodes
  for (int i=0; i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    
    typedef map<int,double>::const_iterator CI;
    map<int,double> derivdmap = cnode->GetDerivD();
    map<int,double > dmap = cnode->GetD()[0];
    
    cout << endl << "Node: " << cnode->Id() << "  Owner: " << cnode->Owner() << endl;
    
    // store D-value into refD
    refD[i]=dmap[cnode->Dofs()[0]];
   
    cout << "D-derivative-map: " << endl;
    for (CI p=derivdmap.begin();p!=derivdmap.end();++p)
      cout << p->first << '\t' << p->second << endl;
  }
      
  // global loop to apply FD scheme to all slave dofs (=2*nodes)
  for (int i=0; i<2*snodefullmap_->NumMyElements();++i)
  {
    // reset Mortar map D
    // loop over all nodes (use fully overlapping column map)
    for (int k=0;k<idiscret_->NumMyColNodes();++k)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(k));
      
      // reset nodal Mortar map D
      for (int j=0;j<(int)((node->GetD()).size());++j)
        (node->GetD())[j].clear();

      (node->GetD()).resize(0);
    }  
    
    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap_->GID(i/2);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);
    
    // apply finite difference scheme
    if (Comm().MyPID()==snode->Owner())
    {
      cout << "\nBuilding FD for Slave Node: " << snode->Id() << " Dof(l): " << i%2
           << " Dof(g): " << snode->Dofs()[i%2] << endl;
    }
    
    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (i%2==0)
    {
      snode->xspatial()[0] += delta;
      snode->u()[0] += delta;
    }
    else
    {
      snode->xspatial()[1] += delta;
      snode->u()[1] += delta;
    }
    
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int j=0;j<idiscret_->NumMyColElements();++j)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(j));
      element->Area()=element->ComputeArea();
    } 
        
    // compute new D-Matrix entries
    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int j=0; j<selecolmap_->NumMyElements();++j)
    {
      int gid1 = selecolmap_->GID(j);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);
      
  #ifndef CONTACTONEMORTARLOOP
      // integrate Mortar matrix D (lives on slave side only!)
      IntegrateSlave2D(*selement);
  #endif // #ifndef CONTACTONEMORTARLOOP
    }
    
    // compute finite difference derivative
    for (int k=0;k<snoderowmap_->NumMyElements();++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %",kgid);
      CNode* kcnode = static_cast<CNode*>(knode);
      map<int,double > dmap = kcnode->GetD()[0];
      
      newD[k]=dmap[kcnode->Dofs()[0]];
              
      // print results (derivatives) to screen
      if (abs(newD[k]-refD[k])>1e-12)
      {
        cout << "Node: " << kcnode->Id() << "  Owner: " << kcnode->Owner() << endl;
        //cout << "Ref-D: " << refD[k] << endl;
        //cout << "New-D: " << newD[k] << endl;
        cout << "Deriv: " << snode->Dofs()[i%2] << " " << (newD[k]-refD[k])/delta << endl;
      }
    }
    
    // undo finite difference modification
    if (i%2==0)
    {
      snode->xspatial()[0] -= delta;
      snode->u()[0] -= delta;
    }
    else
    {
      snode->xspatial()[1] -= delta;
      snode->u()[1] -= delta;
    }
  }
  
  // back to normal(1)...
  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));
    
    // reset nodal Mortar maps
    for (int j=0;j<(int)((node->GetD()).size());++j)
      (node->GetD())[j].clear();

    (node->GetD()).resize(0);
  }
      
  // back to normal(2)...
  // loop over all elements to set current element length / area
  // (use fully overlapping column map)
  for (int j=0;j<idiscret_->NumMyColElements();++j)
  {
    CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(j));
    element->Area()=element->ComputeArea();
  }
  
  // back to normal(3)...
  // compute new D-Matrix entries
  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int j=0; j<selecolmap_->NumMyElements();++j)
  {
    int gid1 = selecolmap_->GID(j);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CElement* selement = static_cast<CElement*>(ele1);
    
#ifndef CONTACTONEMORTARLOOP
    // integrate Mortar matrix D (lives on slave side only!)
    IntegrateSlave2D(*selement);
#endif // #ifndef CONTACTONEMORTARLOOP
  }
      
  //exit(0);
  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for M-Mortar derivatives           popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckMortarMDeriv()
{
  // create storage for M-Matrix entries
  int nrow = snoderowmap_->NumMyElements();
  vector<map<int,double> > refM(nrow);
  vector<map<int,double> > newM(nrow);
  
  // print reference to screen (M-derivative-maps)
  // loop over proc's slave nodes
  for (int i=0; i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    
    typedef map<int,map<int,double> >::const_iterator CIM;
    typedef map<int,double>::const_iterator CI;
    map<int,map<int,double> >& derivmmap = cnode->GetDerivM();
    
    map<int,double > mmap = cnode->GetM()[0];
    
    cout << endl << "Node: " << cnode->Id() << "  Owner: " << cnode->Owner() << endl;
    
    // store M-values into refM
    refM[i]=mmap;
   
    // print to screen
    for (CIM p=derivmmap.begin();p!=derivmmap.end();++p)
    {
      cout << "M-derivative-map for pair S" << cnode->Id() << " and M" << p->first << endl;
      map<int,double>& currmap = derivmmap[p->first];
      for (CI q=currmap.begin();q!=currmap.end();++q)
        cout << q->first << "\t" << q->second << endl;
    }
  }
  
  
  // global loop to apply FD scheme to all slave dofs (=2*nodes)
  for (int fd=0; fd<2*snodefullmap_->NumMyElements();++fd)
  {
    // Initialize
    Initialize();
    
    // now get the node we want to apply the FD scheme to
    int gid = snodefullmap_->GID(fd/2);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);
    
    // apply finite difference scheme
    if (Comm().MyPID()==snode->Owner())
    {
      cout << "\nBuilding FD for Slave Node: " << snode->Id() << " Dof(l): " << fd%2
           << " Dof(g): " << snode->Dofs()[fd%2] << endl;
    }
    
    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%2==0)
    {
      snode->xspatial()[0] += delta;
      snode->u()[0] += delta;
    }
    else
    {
      snode->xspatial()[1] += delta;
      snode->u()[1] += delta;
    }
    
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int j=0;j<idiscret_->NumMyColElements();++j)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(j));
      element->Area()=element->ComputeArea();
    }
    
    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************    
    // loop over proc's slave nodes of the interface
    // use standard column map to include processor's ghosted nodes
    // use boundary map to include slave side boundary nodes
    for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
    {
      int gid1 = snodecolmapbound_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid1);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid1);
      CNode* cnode = static_cast<CNode*>(node);
      
      // build averaged normal at each slave node
      cnode->BuildAveragedNormal();
    }

    // contact search algorithm
    EvaluateContactSearch();

    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);

  #ifndef CONTACTONEMORTARLOOP
      // integrate Mortar matrix D (lives on slave side only!)
      IntegrateSlave2D(*selement);
  #else
  #ifdef CONTACTFULLLIN
      dserror("ERROR: Full linearization not yet implemented for 1 mortar loop case!");
  #endif // #ifdef CONTACTFULLLIN
  #endif // #ifndef CONTACTONEMORTARLOOP

      // loop over the contact candidate master elements
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);

        // prepare overlap integration
        vector<bool> hasproj(4);
        vector<double> xiproj(4);
        bool overlap = false;

        // project the element pair
        Project2D(*selement,*melement,hasproj,xiproj);

        // check for element overlap
        overlap = DetectOverlap2D(*selement,*melement,hasproj,xiproj);

        // integrate the element overlap
        if (overlap) IntegrateOverlap2D(*selement,*melement,xiproj);
      }
    }
    // *******************************************************************
    
    // compute finite difference derivative
    for (int k=0;k<snoderowmap_->NumMyElements();++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %",kgid);
      CNode* kcnode = static_cast<CNode*>(knode);
      map<int,double> mmap = kcnode->GetM()[0];
      typedef map<int,double>::const_iterator CI;
      
      // store M-values into refM
      newM[k]=mmap;
              
      // print results (derivatives) to screen
      for (CI p=newM[k].begin();p!=newM[k].end();++p)
      {
        if (abs(newM[k][p->first]-refM[k][p->first]) > 1e-12)
        {
          cout << "M-FD-derivative for pair S" << kcnode->Id() << " and M" << (p->first)/2 << endl;
          //cout << "Ref-M: " << refM[k][p->first] << endl;
          //cout << "New-M: " << newM[k][p->first] << endl;
          cout << "Deriv: " << snode->Dofs()[fd%2] << " " << (newM[k][p->first]-refM[k][p->first])/delta << endl;
        }
      }
    }
    
    // undo finite difference modification
    if (fd%2==0)
    {
      snode->xspatial()[0] -= delta;
      snode->u()[0] -= delta;
    }
    else
    {
      snode->xspatial()[1] -= delta;
      snode->u()[1] -= delta;
    }       
  }
  
  // global loop to apply FD scheme to all master dofs (=2*nodes)
  for (int fd=0; fd<2*mnodefullmap_->NumMyElements();++fd)
  {
    // Initialize
    Initialize();
    
    // now get the node we want to apply the FD scheme to
    int gid = mnodefullmap_->GID(fd/2);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find master node with gid %",gid);
    CNode* mnode = static_cast<CNode*>(node);
    
    // apply finite difference scheme
    if (Comm().MyPID()==mnode->Owner())
    {
      cout << "\nBuilding FD for Master Node: " << mnode->Id() << " Dof(l): " << fd%2
           << " Dof(g): " << mnode->Dofs()[fd%2] << endl;
    }
    
    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (fd%2==0)
    {
      mnode->xspatial()[0] += delta;
      mnode->u()[0] += delta;
    }
    else
    {
      mnode->xspatial()[1] += delta;
      mnode->u()[1] += delta;
    }
    
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int j=0;j<idiscret_->NumMyColElements();++j)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(j));
      element->Area()=element->ComputeArea();
    }
    
    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************    
    // loop over proc's slave nodes of the interface
    // use standard column map to include processor's ghosted nodes
    // use boundary map to include slave side boundary nodes
    for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
    {
      int gid1 = snodecolmapbound_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid1);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid1);
      CNode* cnode = static_cast<CNode*>(node);
      
      // build averaged normal at each slave node
      cnode->BuildAveragedNormal();
    }

    // contact search algorithm
    EvaluateContactSearch();

    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);

  #ifndef CONTACTONEMORTARLOOP
      // integrate Mortar matrix D (lives on slave side only!)
      IntegrateSlave2D(*selement);
  #else
  #ifdef CONTACTFULLLIN
      dserror("ERROR: Full linearization not yet implemented for 1 mortar loop case!");
  #endif // #ifdef CONTACTFULLLIN
  #endif // #ifndef CONTACTONEMORTARLOOP

      // loop over the contact candidate master elements
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);

        // prepare overlap integration
        vector<bool> hasproj(4);
        vector<double> xiproj(4);
        bool overlap = false;

        // project the element pair
        Project2D(*selement,*melement,hasproj,xiproj);

        // check for element overlap
        overlap = DetectOverlap2D(*selement,*melement,hasproj,xiproj);

        // integrate the element overlap
        if (overlap) IntegrateOverlap2D(*selement,*melement,xiproj);
      }
    }
    // *******************************************************************
    
    // compute finite difference derivative
    for (int k=0;k<snoderowmap_->NumMyElements();++k)
    {
      int kgid = snoderowmap_->GID(k);
      DRT::Node* knode = idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %",kgid);
      CNode* kcnode = static_cast<CNode*>(knode);
      map<int,double> mmap = kcnode->GetM()[0];
      typedef map<int,double>::const_iterator CI;
      
      // store M-values into refM
      newM[k]=mmap;
              
      // print results (derivatives) to screen
      for (CI p=newM[k].begin();p!=newM[k].end();++p)
      {
        if (abs(newM[k][p->first]-refM[k][p->first]) > 1e-12)
        {
          cout << "M-FD-derivative for pair S" << kcnode->Id() << " and M" << (p->first)/2 << endl;
          //cout << "Ref-M: " << refM[k][p->first] << endl;
          //cout << "New-M: " << newM[k][p->first] << endl;
          cout << "Deriv: " << mnode->Dofs()[fd%2] << " " << (newM[k][p->first]-refM[k][p->first])/delta << endl;
        }
      }
    }
    
    // undo finite difference modification
    if (fd%2==0)
    {
      mnode->xspatial()[0] -= delta;
      mnode->u()[0] -= delta;
    }
    else
    {
      mnode->xspatial()[1] -= delta;
      mnode->u()[1] -= delta;
    }       
  }
  
  // back to normal...
  
  // Initialize
  Initialize();
  
  // loop over all elements to set current element length / area
  // (use fully overlapping column map)
  for (int j=0;j<idiscret_->NumMyColElements();++j)
  {
    CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(j));
    element->Area()=element->ComputeArea();
  }
    
  // *******************************************************************
  // contents of Evaluate()
  // *******************************************************************    
  // loop over proc's slave nodes of the interface
  // use standard column map to include processor's ghosted nodes
  // use boundary map to include slave side boundary nodes
  for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
  {
    int gid1 = snodecolmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid1);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid1);
    CNode* cnode = static_cast<CNode*>(node);
    
    // build averaged normal at each slave node
    cnode->BuildAveragedNormal();
  }

  // contact search algorithm
  EvaluateContactSearch();

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CElement* selement = static_cast<CElement*>(ele1);

#ifndef CONTACTONEMORTARLOOP
    // integrate Mortar matrix D (lives on slave side only!)
    IntegrateSlave2D(*selement);
#else
#ifdef CONTACTFULLLIN
    dserror("ERROR: Full linearization not yet implemented for 1 mortar loop case!");
#endif // #ifdef CONTACTFULLLIN
#endif // #ifndef CONTACTONEMORTARLOOP

    // loop over the contact candidate master elements
    // use slave element's candidate list SearchElements !!!
    for (int j=0;j<selement->NumSearchElements();++j)
    {
      int gid2 = selement->SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      CElement* melement = static_cast<CElement*>(ele2);

      // prepare overlap integration
      vector<bool> hasproj(4);
      vector<double> xiproj(4);
      bool overlap = false;

      // project the element pair
      Project2D(*selement,*melement,hasproj,xiproj);

      // check for element overlap
      overlap = DetectOverlap2D(*selement,*melement,hasproj,xiproj);

      // integrate the element overlap
      if (overlap) IntegrateOverlap2D(*selement,*melement,xiproj);
    }
  }
  // *******************************************************************
  
  //exit(0);
  return;
}

#endif  // #ifdef CCADISCRET
