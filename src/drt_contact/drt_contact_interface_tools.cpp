/*!----------------------------------------------------------------------
\file drt_contact_interface_tools.cpp
\brief

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
#include "drt_contact_coupling.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"


/*----------------------------------------------------------------------*
 |  Visualize contact stuff with gmsh                         popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::VisualizeGmsh(const Epetra_SerialDenseMatrix& csegs,
                                       const int step, const int iter,
                                       const bool fric)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
  // construct unique filename for gmsh output
  // first index = time step index
  std::ostringstream filename;
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename << "o/gmsh_output/" << filebase << "_";
  if (step<10)
    filename << 0 << 0 << 0 << 0;
  else if (step<100)
    filename << 0 << 0 << 0;
  else if (step<1000)
    filename << 0 << 0;
  else if (step<10000)
    filename << 0;
  else if (step>99999)
    dserror("Gmsh output implemented for a maximum of 99.999 time steps");
  filename << step;
  
  // construct unique filename for gmsh output
  // second index = Newton iteration index
#ifdef CONTACTGMSH2
  filename << "_";
  if (iter<10)
    filename << 0;
  else if (iter>99)
    dserror("Gmsh output implemented for a maximum of 99 iterations");
  filename << iter;
#endif // #ifdef CONTACTGMSH2
  
  filename << ".pos";

  // do output to file in c-style
  FILE* fp = NULL;

  //**********************************************************************
  // Start GMSH output
  //**********************************************************************
  for (int proc=0;proc<lComm()->NumProc();++proc)
  {
    if (proc==lComm()->MyPID())
    {
      // open file (overwrite if proc==0, else append)
      if (proc==0) fp = fopen(filename.str().c_str(), "w");
      else fp = fopen(filename.str().c_str(), "a");

      // write output to temporary stringstream
      std::stringstream gmshfilecontent;
      if (proc==0) gmshfilecontent << "View \" Step " << step << " Iter " << iter << " \" {" << endl;

      //******************************************************************
      // plot elements
      //******************************************************************
      for (int i=0; i<idiscret_->NumMyRowElements(); ++i)
      {
        CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lRowElement(i));
        int nnodes = element->NumNode();
        LINALG::SerialDenseMatrix coord(3,nnodes);
        double color = (double)element->IsSlave();

        //local center
        double xi[2] = {0.0, 0.0};
        
        // 2D linear case (2noded line elements)
        if (element->Shape()==DRT::Element::line2)
        {
          element->GetNodalCoords(coord);

          gmshfilecontent << "SL(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "};" << endl;
        }

        // 2D quadratic case (3noded line elements)
        if (element->Shape()==DRT::Element::line3)
        {
          element->GetNodalCoords(coord);

          gmshfilecontent << "SL2(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                              << coord(2,2) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "," << color << "};" << endl;
        }
        
        // 3D linear case (3noded triangular elements)
        if (element->Shape()==DRT::Element::tri3)
        {
          element->GetNodalCoords(coord);
  
          gmshfilecontent << "ST(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                              << coord(2,2) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "," << color << "};" << endl;
          xi[0] = 1.0/3; xi[1] = 1.0/3;
        }
        
        // 3D bilinear case (4noded quadrilateral elements)
        if (element->Shape()==DRT::Element::quad4)
        {
          element->GetNodalCoords(coord);
  
          gmshfilecontent << "SQ(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                              << coord(2,2) << "," << coord(0,3) << "," << coord(1,3) << ","
                              << coord(2,3) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "," << color << "," << color << "};" << endl;
        }
        
        // plot element number in element center
        double elec[3];
        element->LocalToGlobal(xi,elec,0);
        gmshfilecontent << "T3(" << scientific << elec[0] << "," << elec[1] << "," << elec[2] << "," << 17 << ")";
        if (element->IsSlave())
          gmshfilecontent << "{" << "S" << element->Id() << "};" << endl;
        else
          gmshfilecontent << "{" << "M" << element->Id() << "};" << endl;
      }

      //******************************************************************
      // plot normal, tangent, contact status and Lagr. mutlipliers
      //******************************************************************
      for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
      {
        int gid = snoderowmap_->GID(i);
        DRT::Node* node = idiscret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cnode = static_cast<CNode*>(node);
        if (!cnode) dserror("ERROR: Static Cast to CNode* failed");

        double nc[3];
        double nn[3];
        double nt[3];
        double lmn = 0.0;
        double lmt = 0.0;
        
        for (int j=0;j<3;++j)
        {
          nc[j]=cnode->xspatial()[j];
          nn[j]=cnode->n()[j];
          nt[j]=cnode->txi()[j];
          lmn +=  (cnode->Active())*nn[j]* cnode->lm()[j];
          lmt +=  (cnode->Active())*nt[j]* cnode->lm()[j];
        }

        //******************************************************************
        // plot normal and tangent vectors (only 2D)
        //******************************************************************
        gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
        gmshfilecontent << "{" << scientific << nn[0] << "," << nn[1] << "," << nn[2] << "};" << endl;
         
        if (fric)
        {
          gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
          gmshfilecontent << "{" << scientific << nt[0] << "," << nt[1] << "," << nt[2] << "};" << endl;
        }
         
        //******************************************************************
        // plot contact status of slave nodes (inactive, active, stick, slip)
        //******************************************************************
        // frictionless contact, active node = {A}
        if (!fric && cnode->Active())
        {
          gmshfilecontent << "T3(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << "," << 17 << ")";
          gmshfilecontent << "{" << "A" << "};" << endl;
        }
        
        // frictionless contact, inactive node = { }
        else if (!fric && !cnode->Active())
        {
          //do nothing
        }
        
        // frictional contact, slip node = {G}
        else if (fric && cnode->Active() && cnode->Slip())
        {
          gmshfilecontent << "T3(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << "," << 17 << ")";
          gmshfilecontent << "{" << "G" << "};" << endl;
        }
        
        //frictional contact, stick node = {H}
        else if (fric && cnode->Active() && !cnode->Slip())
        {
          gmshfilecontent << "T3(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << "," << 17 << ")";
          gmshfilecontent << "{" << "H" << "};" << endl;
        }
        
        //******************************************************************
        // plot Lagrange multipliers (normal+tang. contact stresses)
        //******************************************************************
        gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
        gmshfilecontent << "{" << scientific << lmn*nn[0] << "," << lmn*nn[1] << "," << lmn*nn[2] << "};" << endl;
        
        if (fric)
        {
          gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
          gmshfilecontent << "{" << scientific << lmt*nn[0] << "," << lmt*nn[1] << "," << lmt*nn[2] << "};" << endl;
        }
      }
      
      //******************************************************************
      // plot contact segments (search + projections)
      //******************************************************************
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
      
      if (proc==lComm()->NumProc()-1) gmshfilecontent << "};" << endl;

      // move everything to gmsh post-processing file and close it
      fprintf(fp,gmshfilecontent.str().c_str());
      fclose(fp);
    }
    lComm()->Barrier();
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Visualize contact stuff with gmsh                         popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::VisualizeGmshLight()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
  // construct unique filename for gmsh output
  // first index = time step index
  std::ostringstream filename;
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename << "o/gmsh_output/" << "normals_" << filebase;
  filename << ".pos";

  // do output to file in c-style
  FILE* fp = NULL;

  //**********************************************************************
  // Start GMSH output
  //**********************************************************************
  for (int proc=0;proc<lComm()->NumProc();++proc)
  {
    if (proc==lComm()->MyPID())
    {
      // open file (overwrite if proc==0, else append)
      if (proc==0) fp = fopen(filename.str().c_str(), "w");
      else fp = fopen(filename.str().c_str(), "a");

      // write output to temporary stringstream
      std::stringstream gmshfilecontent;
      if (proc==0) gmshfilecontent << "View \" Normals \" {" << endl;

      //******************************************************************
      // plot elements
      //******************************************************************
      for (int i=0; i<idiscret_->NumMyRowElements(); ++i)
      {
        CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lRowElement(i));
        int nnodes = element->NumNode();
        LINALG::SerialDenseMatrix coord(3,nnodes);
        double color = (double)element->IsSlave();

        // 2D linear case (2noded line elements)
        if (element->Shape()==DRT::Element::line2)
        {
          element->GetNodalCoords(coord);

          gmshfilecontent << "SL(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "};" << endl;
        }

        // 2D quadratic case (3noded line elements)
        if (element->Shape()==DRT::Element::line3)
        {
          element->GetNodalCoords(coord);

          gmshfilecontent << "SL2(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                              << coord(2,2) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "," << color << "};" << endl;
        }
        
        // 3D linear case (3noded triangular elements)
        if (element->Shape()==DRT::Element::tri3)
        {
          element->GetNodalCoords(coord);
  
          gmshfilecontent << "ST(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                              << coord(2,2) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "," << color << "};" << endl;
        }
        
        // 3D bilinear case (4noded quadrilateral elements)
        if (element->Shape()==DRT::Element::quad4)
        {
          element->GetNodalCoords(coord);
  
          gmshfilecontent << "SQ(" << scientific << coord(0,0) << "," << coord(1,0) << ","
                              << coord(2,0) << "," << coord(0,1) << "," << coord(1,1) << ","
                              << coord(2,1) << "," << coord(0,2) << "," << coord(1,2) << ","
                              << coord(2,2) << "," << coord(0,3) << "," << coord(1,3) << ","
                              << coord(2,3) << ")";
          gmshfilecontent << "{" << scientific << color << "," << color << "," << color << "," << color << "};" << endl;
        }
      }

      //******************************************************************
      // plot normal + tangent vectors (2D or 3D)
      //******************************************************************
      for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
      {
        int gid = snoderowmap_->GID(i);
        DRT::Node* node = idiscret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cnode = static_cast<CNode*>(node);
        if (!cnode) dserror("ERROR: Static Cast to CNode* failed");

        double nc[3];
        double nn[3];
        double ntxi[3];
        double nteta[3];
        
        for (int j=0;j<3;++j)
        {
          nc[j]=cnode->xspatial()[j];
          nn[j]=cnode->n()[j];
          ntxi[j]=cnode->txi()[j];
          nteta[j]=cnode->teta()[j];
        }

        gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
        gmshfilecontent << "{" << scientific << 4*nn[0] << "," << 4*nn[1] << "," << 4*nn[2] << "};" << endl;
        
        gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
        gmshfilecontent << "{" << scientific << 2*ntxi[0] << "," << 2*ntxi[1] << "," << 2*ntxi[2] << "};" << endl;
        
        gmshfilecontent << "VP(" << scientific << nc[0] << "," << nc[1] << "," << nc[2] << ")";
        gmshfilecontent << "{" << scientific << 3*nteta[0] << "," << 3*nteta[1] << "," << 3*nteta[2] << "};" << endl;
      }
      
      if (proc==lComm()->NumProc()-1) gmshfilecontent << "};" << endl;

      // move everything to gmsh post-processing file and close it
      fprintf(fp,gmshfilecontent.str().c_str());
      fclose(fp);
    }
    lComm()->Barrier();
  }

  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for normal/tangent deriv.          popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckNormalDeriv()
{
  /****************************************************/
  /* NOTE: This is a combined 2D / 3D method already! */
  /****************************************************/
  
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
  // global loop to apply FD scheme to all slave dofs (=3*nodes)
  for (int i=0; i<3*snodefullmap_->NumMyElements();++i)
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
    vector<double> refnz(int(snodecolmapbound_->NumMyElements()));
    vector<double> newnx(int(snodecolmapbound_->NumMyElements()));
    vector<double> newny(int(snodecolmapbound_->NumMyElements()));
    vector<double> newnz(int(snodecolmapbound_->NumMyElements()));
    
    vector<double> reftxix(int(snodecolmapbound_->NumMyElements()));
    vector<double> reftxiy(int(snodecolmapbound_->NumMyElements()));
    vector<double> reftxiz(int(snodecolmapbound_->NumMyElements()));
    vector<double> newtxix(int(snodecolmapbound_->NumMyElements()));
    vector<double> newtxiy(int(snodecolmapbound_->NumMyElements()));
    vector<double> newtxiz(int(snodecolmapbound_->NumMyElements()));
    
    vector<double> reftetax(int(snodecolmapbound_->NumMyElements()));
    vector<double> reftetay(int(snodecolmapbound_->NumMyElements()));
    vector<double> reftetaz(int(snodecolmapbound_->NumMyElements()));
    vector<double> newtetax(int(snodecolmapbound_->NumMyElements()));
    vector<double> newtetay(int(snodecolmapbound_->NumMyElements()));
    vector<double> newtetaz(int(snodecolmapbound_->NumMyElements()));
        
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
        cout << "Row dof id: " << jcnode->Dofs()[2] << endl;
        for (CI p=(jcnode->GetDerivN()[2]).begin();p!=(jcnode->GetDerivN()[2]).end();++p)
          cout << p->first << '\t' << p->second << endl;
        
        cout << "Tangent txi-derivative-maps: " << endl;
        cout << "Row dof id: " << jcnode->Dofs()[0] << endl;
        for (CI p=(jcnode->GetDerivTxi()[0]).begin();p!=(jcnode->GetDerivTxi()[0]).end();++p)
          cout << p->first << '\t' << p->second << endl;    
        cout << "Row dof id: " << jcnode->Dofs()[1] << endl;
        for (CI p=(jcnode->GetDerivTxi()[1]).begin();p!=(jcnode->GetDerivTxi()[1]).end();++p)
          cout << p->first << '\t' << p->second << endl;
        cout << "Row dof id: " << jcnode->Dofs()[2] << endl;
        for (CI p=(jcnode->GetDerivTxi()[2]).begin();p!=(jcnode->GetDerivTxi()[2]).end();++p)
          cout << p->first << '\t' << p->second << endl;
        
        cout << "Tangent teta-derivative-maps: " << endl;
        cout << "Row dof id: " << jcnode->Dofs()[0] << endl;
        for (CI p=(jcnode->GetDerivTeta()[0]).begin();p!=(jcnode->GetDerivTeta()[0]).end();++p)
          cout << p->first << '\t' << p->second << endl;    
        cout << "Row dof id: " << jcnode->Dofs()[1] << endl;
        for (CI p=(jcnode->GetDerivTeta()[1]).begin();p!=(jcnode->GetDerivTeta()[1]).end();++p)
          cout << p->first << '\t' << p->second << endl;
        cout << "Row dof id: " << jcnode->Dofs()[2] << endl;
        for (CI p=(jcnode->GetDerivTeta()[2]).begin();p!=(jcnode->GetDerivTeta()[2]).end();++p)
          cout << p->first << '\t' << p->second << endl;
      }
      
      // store reference normals / tangents
      refnx[j] = jcnode->n()[0];
      refny[j] = jcnode->n()[1];
      refnz[j] = jcnode->n()[2];
      reftxix[j] = jcnode->txi()[0];
      reftxiy[j] = jcnode->txi()[1];
      reftxiz[j] = jcnode->txi()[2];
      reftetax[j] = jcnode->teta()[0];
      reftetay[j] = jcnode->teta()[1];
      reftetaz[j] = jcnode->teta()[2];
      
    }
    
    // now fincally get the node we want to apply the FD scheme to
    int gid = snodefullmap_->GID(i/3);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);
    
    // apply finite difference scheme
    cout << "\nBuilding FD for Slave Node: " << snode->Id() << " Dof(l): " << i%3
         << " Dof(g): " << snode->Dofs()[i%3] << endl;
    
    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (i%3==0)
    {
      snode->xspatial()[0] += delta;
      snode->u()[0] += delta;
    }
    else if (i%3==1)
    {
      snode->xspatial()[1] += delta;
      snode->u()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
      snode->u()[2] += delta;
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
      newnz[k] = kcnode->n()[2];
      newtxix[k] = kcnode->txi()[0];
      newtxiy[k] = kcnode->txi()[1];
      newtxiz[k] = kcnode->txi()[2];
      newtetax[k] = kcnode->teta()[0];
      newtetay[k] = kcnode->teta()[1];
      newtetaz[k] = kcnode->teta()[2];
      
      // get reference normal / tangent
      double refn[3] = {0.0, 0.0, 0.0};
      double reftxi[3] = {0.0, 0.0, 0.0};
      double refteta[3] = {0.0, 0.0, 0.0};
      refn[0] = refnx[k];
      refn[1] = refny[k];
      refn[2] = refnz[k];
      reftxi[0] = reftxix[k];
      reftxi[1] = reftxiy[k];
      reftxi[2] = reftxiz[k];
      refteta[0] = reftetax[k];
      refteta[1] = reftetay[k];
      refteta[2] = reftetaz[k];
      
      // get modified normal / tangent
      double newn[3] = {0.0, 0.0, 0.0};
      double newtxi[3] = {0.0, 0.0, 0.0};
      double newteta[3] = {0.0, 0.0, 0.0};
      newn[0] = newnx[k];
      newn[1] = newny[k];
      newn[2] = newnz[k];
      newtxi[0] = newtxix[k];
      newtxi[1] = newtxiy[k];
      newtxi[2] = newtxiz[k];
      newteta[0] = newtetax[k];
      newteta[1] = newtetay[k];
      newteta[2] = newtetaz[k];
      
      // print results (derivatives) to screen
      if (abs(newn[0]-refn[0])>1e-12 || abs(newn[1]-refn[1])>1e-12 || abs(newn[2]-refn[2])>1e-12)
      {
        cout << "Node: " << kcnode->Id() << "  Owner: " << kcnode->Owner() << endl;
        cout << "Normal derivative (FD):" << endl;
        if (abs(newn[0]-refn[0])>1e-12)
        { 
          double val = (newn[0]-refn[0])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[0] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
      
        if (abs(newn[1]-refn[1])>1e-12)
        {
          double val = (newn[1]-refn[1])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[1] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
        
        if (abs(newn[2]-refn[2])>1e-12)
        {
          double val = (newn[2]-refn[2])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[2] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
      }
      
      if (abs(newtxi[0]-reftxi[0])>1e-12 || abs(newtxi[1]-reftxi[1])>1e-12 | abs(newtxi[2]-reftxi[2])>1e-12)
      {
        cout << "Node: " << kcnode->Id() << "  Owner: " << kcnode->Owner() << endl;
        cout << "Tangent txi derivative (FD):" << endl;
        if (abs(newtxi[0]-reftxi[0])>1e-12)
        { 
          double val = (newtxi[0]-reftxi[0])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[0] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
      
        if (abs(newtxi[1]-reftxi[1])>1e-12)
        {
          double val = (newtxi[1]-reftxi[1])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[1] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
        
        if (abs(newtxi[2]-reftxi[2])>1e-12)
        {
          double val = (newtxi[2]-reftxi[2])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[2] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
      }
      
      if (abs(newteta[0]-refteta[0])>1e-12 || abs(newteta[1]-refteta[1])>1e-12 | abs(newteta[2]-refteta[2])>1e-12)
      {
        cout << "Node: " << kcnode->Id() << "  Owner: " << kcnode->Owner() << endl;
        cout << "Tangent teta derivative (FD):" << endl;
        if (abs(newteta[0]-refteta[0])>1e-12)
        { 
          double val = (newteta[0]-refteta[0])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[0] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
      
        if (abs(newteta[1]-refteta[1])>1e-12)
        {
          double val = (newteta[1]-refteta[1])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[1] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
        
        if (abs(newteta[2]-refteta[2])>1e-12)
        {
          double val = (newteta[2]-refteta[2])/delta;
          cout << "Row dof id: " << kcnode->Dofs()[2] << endl;
          cout << snode->Dofs()[i%3] << '\t' << val << endl;
        }
      }
      
    }
    
    // undo finite difference modification
    if (i%3==0)
    {
      snode->xspatial()[0] -= delta;
      snode->u()[0] -= delta;
    }
    else if (i%3==1)
    {
      snode->xspatial()[1] -= delta;
      snode->u()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
      snode->u()[2] -= delta;
    }
  }
 
  // back to normal...
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
      
  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for D-Mortar derivatives           popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckMortarDDeriv()
{
  /****************************************************/
  /* NOTE: This is a combined 2D / 3D method already! */
  /****************************************************/
  
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
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
      
  // global loop to apply FD scheme to all slave dofs (=3*nodes)
  for (int i=0; i<3*snodefullmap_->NumMyElements();++i)
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
    int gid = snodefullmap_->GID(i/3);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %",gid);
    CNode* snode = static_cast<CNode*>(node);
    
    // apply finite difference scheme
    if (Comm().MyPID()==snode->Owner())
    {
      cout << "\nBuilding FD for Slave Node: " << snode->Id() << " Dof(l): " << i%3
           << " Dof(g): " << snode->Dofs()[i%3] << endl;
    }
    
    // do step forward (modify nodal displacement)
    double delta = 1e-8;
    if (i%3==0)
    {
      snode->xspatial()[0] += delta;
      snode->u()[0] += delta;
    }
    else if (i%3==1)
    {
      snode->xspatial()[1] += delta;
      snode->u()[1] += delta;
    }
    else
    {
      snode->xspatial()[2] += delta;
      snode->u()[2] += delta;
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
      IntegrateSlave(*selement);
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
        cout << "Deriv: " << snode->Dofs()[i%3] << " " << (newD[k]-refD[k])/delta << endl;
      }
    }
    
    // undo finite difference modification
    if (i%3==0)
    {
      snode->xspatial()[0] -= delta;
      snode->u()[0] -= delta;
    }
    else if (i%3==1)
    {
      snode->xspatial()[1] -= delta;
      snode->u()[1] -= delta;
    }
    else
    {
      snode->xspatial()[2] -= delta;
      snode->u()[2] -= delta;
    }
  }
  
  // back to normal...
  
  // Initialize
  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

    //reset nodal normal vector
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
      node->txi()[j]=0.0;
      node->teta()[j]=0.0;
    }

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((node->GetDerivN()).size());++j)
      (node->GetDerivN())[j].clear();
    (node->GetDerivN()).resize(0);
    
    // reset derivative maps of tangent vector
    for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
      (node->GetDerivTxi())[j].clear();
    (node->GetDerivTxi()).resize(0);
    for (int j=0;j<(int)((node->GetDerivTeta()).size());++j)
      (node->GetDerivTeta())[j].clear();
    (node->GetDerivTeta()).resize(0);
        
    // reset closest node
    // (FIXME: at the moment we do not need this info. in the next
    // iteration, but it might be helpful for accelerated search!!!)
    node->ClosestNode() = -1;

    // reset nodal Mortar maps
    for (int j=0;j<(int)((node->GetD()).size());++j)
      (node->GetD())[j].clear();
    for (int j=0;j<(int)((node->GetM()).size());++j)
      (node->GetM())[j].clear();
    for (int j=0;j<(int)((node->GetMmod()).size());++j)
      (node->GetMmod())[j].clear();

    (node->GetD()).resize(0);
    (node->GetM()).resize(0);
    (node->GetMmod()).resize(0);

    // reset derivative map of Mortar matrices
    (node->GetDerivD()).clear();
    (node->GetDerivM()).clear();
    
    // reset nodal weighted gap
    node->Getg() = 1.0e12;

    // reset feasible projection status
    node->HasProj() = false;
  }

  // loop over all elements to reset contact candidates / search lists
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
    CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
    element->SearchElements().resize(0);
  }

  // reset matrix containing interface contact segments (gmsh)
  CSegs().Shape(0,0);
  
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
    IntegrateSlave(*selement);
#endif // #ifndef CONTACTONEMORTARLOOP
  }
  
  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CElement* selement = static_cast<CElement*>(ele1);
    
    // loop over the contact candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j=0;j<selement->NumSearchElements();++j)
    {
      int gid2 = selement->SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      CElement* melement = static_cast<CElement*>(ele2);
      
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pair)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      IntegrateCoupling(*selement,*melement);
    }
  }
  // *******************************************************************
  
  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for M-Mortar derivatives           popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckMortarMDeriv()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
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
    
    if ((int)(cnode->GetM().size())==0)
      break;
    
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
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

      //reset nodal normal vector
      for (int j=0;j<3;++j)
      {
        node->n()[j]=0.0;
        node->txi()[j]=0.0;
        node->teta()[j]=0.0;
      }

      // reset derivative maps of normal vector
      for (int j=0;j<(int)((node->GetDerivN()).size());++j)
        (node->GetDerivN())[j].clear();
      (node->GetDerivN()).resize(0);
      
      // reset derivative maps of tangent vector
      for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
        (node->GetDerivTxi())[j].clear();
      (node->GetDerivTxi()).resize(0);
          
      // reset closest node
      // (FIXME: at the moment we do not need this info. in the next
      // iteration, but it might be helpful for accelerated search!!!)
      node->ClosestNode() = -1;

      // reset nodal Mortar maps
      for (int j=0;j<(int)((node->GetD()).size());++j)
        (node->GetD())[j].clear();
      for (int j=0;j<(int)((node->GetM()).size());++j)
        (node->GetM())[j].clear();
      for (int j=0;j<(int)((node->GetMmod()).size());++j)
        (node->GetMmod())[j].clear();

      (node->GetD()).resize(0);
      (node->GetM()).resize(0);
      (node->GetMmod()).resize(0);

      // reset derivative map of Mortar matrices
      (node->GetDerivD()).clear();
      (node->GetDerivM()).clear();
      
      // reset nodal weighted gap
      node->Getg() = 1.0e12;

      // reset feasible projection status
      node->HasProj() = false;
    }

    // loop over all elements to reset contact candidates / search lists
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->SearchElements().resize(0);
    }

    // reset matrix containing interface contact segments (gmsh)
    CSegs().Shape(0,0);
      
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
      IntegrateSlave(*selement);
  #endif // #ifndef CONTACTONEMORTARLOOP
    }
    
    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);
      
      // loop over the contact candidate master elements of sele_
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);
        
        //********************************************************************
        // 1) perform coupling (projection + overlap detection for sl/m pair)
        // 2) integrate Mortar matrix M and weighted gap g
        // 3) compute directional derivative of M and g and store into nodes
        //********************************************************************
        IntegrateCoupling(*selement,*melement);
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
      
      if ((int)(kcnode->GetM().size())==0)
        break;
      
      map<int,double> mmap = kcnode->GetM()[0];
      
      typedef map<int,double>::const_iterator CI;
      
      // store M-values into refM
      newM[k]=mmap;
              
      // print results (derivatives) to screen
      for (CI p=newM[k].begin();p!=newM[k].end();++p)
      {
        if (abs(newM[k][p->first]-refM[k][p->first]) > 1e-12)
        {
          // if we want to use the FD for DerivM
          //map<int,map<int,double> >& derivm = kcnode->GetDerivM();
          //double val = (newM[k][p->first]-refM[k][p->first])/delta;
          //int outercol = p->first;
          //int innercol = snode->Dofs()[fd%2];
          //(derivm[outercol/2])[innercol] += val;
          
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
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

      //reset nodal normal vector
      for (int j=0;j<3;++j)
      {
        node->n()[j]=0.0;
        node->txi()[j]=0.0;
        node->teta()[j]=0.0;
      }

      // reset derivative maps of normal vector
      for (int j=0;j<(int)((node->GetDerivN()).size());++j)
        (node->GetDerivN())[j].clear();
      (node->GetDerivN()).resize(0);
      
      // reset derivative maps of tangent vector
      for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
        (node->GetDerivTxi())[j].clear();
      (node->GetDerivTxi()).resize(0);
          
      // reset closest node
      // (FIXME: at the moment we do not need this info. in the next
      // iteration, but it might be helpful for accelerated search!!!)
      node->ClosestNode() = -1;

      // reset nodal Mortar maps
      for (int j=0;j<(int)((node->GetD()).size());++j)
        (node->GetD())[j].clear();
      for (int j=0;j<(int)((node->GetM()).size());++j)
        (node->GetM())[j].clear();
      for (int j=0;j<(int)((node->GetMmod()).size());++j)
        (node->GetMmod())[j].clear();

      (node->GetD()).resize(0);
      (node->GetM()).resize(0);
      (node->GetMmod()).resize(0);

      // reset derivative map of Mortar matrices
      (node->GetDerivD()).clear();
      (node->GetDerivM()).clear();
      
      // reset nodal weighted gap
      node->Getg() = 1.0e12;

      // reset feasible projection status
      node->HasProj() = false;
    }

    // loop over all elements to reset contact candidates / search lists
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->SearchElements().resize(0);
    }

    // reset matrix containing interface contact segments (gmsh)
    CSegs().Shape(0,0);
      
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
      IntegrateSlave(*selement);
  #endif // #ifndef CONTACTONEMORTARLOOP
    }
    
    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);
      
      // loop over the contact candidate master elements of sele_
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);
        
        //********************************************************************
        // 1) perform coupling (projection + overlap detection for sl/m pair)
        // 2) integrate Mortar matrix M and weighted gap g
        // 3) compute directional derivative of M and g and store into nodes
        //********************************************************************
        IntegrateCoupling(*selement,*melement);
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
      
      if ((int)(kcnode->GetM().size())==0)
        break;
      
      map<int,double> mmap = kcnode->GetM()[0];
      typedef map<int,double>::const_iterator CI;
      
      // store M-values into refM
      newM[k]=mmap;
              
      // print results (derivatives) to screen
      for (CI p=newM[k].begin();p!=newM[k].end();++p)
      {
        if (abs(newM[k][p->first]-refM[k][p->first]) > 1e-12)
        {
          // if we want to use the FD for DerivM
          //map<int,map<int,double> >& derivm = kcnode->GetDerivM();
          //double val = (newM[k][p->first]-refM[k][p->first])/delta;
          //int outercol = p->first;
          //int innercol = mnode->Dofs()[fd%2];
          //(derivm[outercol/2])[innercol] += val;
          
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
  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

    //reset nodal normal vector
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
      node->txi()[j]=0.0;
      node->teta()[j]=0.0;
    }

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((node->GetDerivN()).size());++j)
      (node->GetDerivN())[j].clear();
    (node->GetDerivN()).resize(0);
    
    // reset derivative maps of tangent vector
    for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
      (node->GetDerivTxi())[j].clear();
    (node->GetDerivTxi()).resize(0);
        
    // reset closest node
    // (FIXME: at the moment we do not need this info. in the next
    // iteration, but it might be helpful for accelerated search!!!)
    node->ClosestNode() = -1;

    // reset nodal Mortar maps
    for (int j=0;j<(int)((node->GetD()).size());++j)
      (node->GetD())[j].clear();
    for (int j=0;j<(int)((node->GetM()).size());++j)
      (node->GetM())[j].clear();
    for (int j=0;j<(int)((node->GetMmod()).size());++j)
      (node->GetMmod())[j].clear();

    (node->GetD()).resize(0);
    (node->GetM()).resize(0);
    (node->GetMmod()).resize(0);

    // reset derivative map of Mortar matrices
    (node->GetDerivD()).clear();
    (node->GetDerivM()).clear();
    
    // reset nodal weighted gap
    node->Getg() = 1.0e12;

    // reset feasible projection status
    node->HasProj() = false;
  }

  // loop over all elements to reset contact candidates / search lists
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
    CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
    element->SearchElements().resize(0);
  }

  // reset matrix containing interface contact segments (gmsh)
  CSegs().Shape(0,0);
  
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
    IntegrateSlave(*selement);
#endif // #ifndef CONTACTONEMORTARLOOP
  }
  
  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CElement* selement = static_cast<CElement*>(ele1);
    
    // loop over the contact candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j=0;j<selement->NumSearchElements();++j)
    {
      int gid2 = selement->SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      CElement* melement = static_cast<CElement*>(ele2);
      
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pair)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      IntegrateCoupling(*selement,*melement);
    }
  }
  // *******************************************************************
  
  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for normal gap derivatives         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckGapDeriv()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
  // create storage for gap values
  int nrow = snoderowmap_->NumMyElements();
  vector<double> refG(nrow);
  vector<double> newG(nrow);
  
  // store reference
  // loop over proc's slave nodes
  for (int i=0; i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    
    if (cnode->Active())
    {
      // check two versions of weighted gap
      double defgap = 0.0;
      double wii = (cnode->GetD()[0])[cnode->Dofs()[0]];
      
      for (int j=0;j<3;++j)
        defgap-= (cnode->n()[j])*wii*(cnode->xspatial()[j]);
      
      vector<map<int,double> > mmap = cnode->GetM();
      map<int,double>::iterator mcurr;
          
      for (int m=0;m<mnodefullmap_->NumMyElements();++m)
      {
        int gid = mnodefullmap_->GID(m);
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CNode* cmnode = static_cast<CNode*>(mnode);
        const int* mdofs = cmnode->Dofs();
        bool hasentry = false;
        
        // look for this master node in M-map of the active slave node
        for (mcurr=mmap[0].begin();mcurr!=mmap[0].end();++mcurr)
          if ((mcurr->first)==mdofs[0])
          {
            hasentry=true;
            break;
          }
        
        double mik = (mmap[0])[mdofs[0]];
        double* mxi = cmnode->xspatial();
        
        // get out of here, if master node not adjacent or coupling very weak
        if (!hasentry || abs(mik)<1.0e-12) continue;
              
        for (int j=0;j<3;++j)
          defgap+= (cnode->n()[j]) * mik * mxi[j];
      }
      
      //cout << "SNode: " << cnode->Id() << " IntGap: " << gap << " DefGap: " << defgap << endl;
      //cnode->Getg() = defgap;
    }
    
    // store gap-values into refG
    refG[i]=cnode->Getg();
  }
  
  // global loop to apply FD scheme to all slave dofs (=2*nodes)
  for (int fd=0; fd<2*snodefullmap_->NumMyElements();++fd)
  {
    // Initialize
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

      //reset nodal normal vector
      for (int j=0;j<3;++j)
      {
        node->n()[j]=0.0;
        node->txi()[j]=0.0;
        node->teta()[j]=0.0;
      }

      // reset derivative maps of normal vector
      for (int j=0;j<(int)((node->GetDerivN()).size());++j)
        (node->GetDerivN())[j].clear();
      (node->GetDerivN()).resize(0);
      
      // reset derivative maps of tangent vector
      for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
        (node->GetDerivTxi())[j].clear();
      (node->GetDerivTxi()).resize(0);
          
      // reset closest node
      // (FIXME: at the moment we do not need this info. in the next
      // iteration, but it might be helpful for accelerated search!!!)
      node->ClosestNode() = -1;

      // reset nodal Mortar maps
      for (int j=0;j<(int)((node->GetD()).size());++j)
        (node->GetD())[j].clear();
      for (int j=0;j<(int)((node->GetM()).size());++j)
        (node->GetM())[j].clear();
      for (int j=0;j<(int)((node->GetMmod()).size());++j)
        (node->GetMmod())[j].clear();

      (node->GetD()).resize(0);
      (node->GetM()).resize(0);
      (node->GetMmod()).resize(0);

      // reset derivative map of Mortar matrices
      (node->GetDerivD()).clear();
      (node->GetDerivM()).clear();
      
      // reset nodal weighted gap
      node->Getg() = 1.0e12;

      // reset feasible projection status
      node->HasProj() = false;
    }

    // loop over all elements to reset contact candidates / search lists
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->SearchElements().resize(0);
    }

    // reset matrix containing interface contact segments (gmsh)
    CSegs().Shape(0,0);
      
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
      IntegrateSlave(*selement);
  #endif // #ifndef CONTACTONEMORTARLOOP
    }
    
    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);
      
      // loop over the contact candidate master elements of sele_
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);
        
        //********************************************************************
        // 1) perform coupling (projection + overlap detection for sl/m pair)
        // 2) integrate Mortar matrix M and weighted gap g
        // 3) compute directional derivative of M and g and store into nodes
        //********************************************************************
        IntegrateCoupling(*selement,*melement);
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
      
      if (kcnode->Active())
      {
        // check two versions of weighted gap
        double defgap = 0.0;
        double wii = (kcnode->GetD()[0])[kcnode->Dofs()[0]];
        
        for (int j=0;j<3;++j)
          defgap-= (kcnode->n()[j])*wii*(kcnode->xspatial()[j]);
        
        vector<map<int,double> > mmap = kcnode->GetM();
        map<int,double>::iterator mcurr;
            
        for (int m=0;m<mnodefullmap_->NumMyElements();++m)
        {
          int gid = mnodefullmap_->GID(m);
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          CNode* cmnode = static_cast<CNode*>(mnode);
          const int* mdofs = cmnode->Dofs();
          bool hasentry = false;
          
          // look for this master node in M-map of the active slave node
          for (mcurr=mmap[0].begin();mcurr!=mmap[0].end();++mcurr)
            if ((mcurr->first)==mdofs[0])
            {
              hasentry=true;
              break;
            }
          
          double mik = (mmap[0])[mdofs[0]];
          double* mxi = cmnode->xspatial();
          
          // get out of here, if master node not adjacent or coupling very weak
          if (!hasentry || abs(mik)<1.0e-12) continue;
                
          for (int j=0;j<3;++j)
            defgap+= (kcnode->n()[j]) * mik * mxi[j];
        }
        
        //cout << "SNode: " << kcnode->Id() << " IntGap: " << gap << " DefGap: " << defgap << endl;
        //kcnode->Getg() = defgap;
      }
      
      // store gap-values into newG
      newG[k]=kcnode->Getg();
              
      // print results (derivatives) to screen
      if (abs(newG[k]-refG[k]) > 1e-12)
      {
        cout << "G-FD-derivative for node S" << kcnode->Id() << endl;
        //cout << "Ref-G: " << refG[k] << endl;
        //cout << "New-G: " << newG[k] << endl;
        cout << "Deriv: " << snode->Dofs()[fd%2] << " " << (newG[k]-refG[k])/delta << endl;
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
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

      //reset nodal normal vector
      for (int j=0;j<3;++j)
      {
        node->n()[j]=0.0;
        node->txi()[j]=0.0;
        node->teta()[j]=0.0;
      }

      // reset derivative maps of normal vector
      for (int j=0;j<(int)((node->GetDerivN()).size());++j)
        (node->GetDerivN())[j].clear();
      (node->GetDerivN()).resize(0);
      
      // reset derivative maps of tangent vector
      for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
        (node->GetDerivTxi())[j].clear();
      (node->GetDerivTxi()).resize(0);
          
      // reset closest node
      // (FIXME: at the moment we do not need this info. in the next
      // iteration, but it might be helpful for accelerated search!!!)
      node->ClosestNode() = -1;

      // reset nodal Mortar maps
      for (int j=0;j<(int)((node->GetD()).size());++j)
        (node->GetD())[j].clear();
      for (int j=0;j<(int)((node->GetM()).size());++j)
        (node->GetM())[j].clear();
      for (int j=0;j<(int)((node->GetMmod()).size());++j)
        (node->GetMmod())[j].clear();

      (node->GetD()).resize(0);
      (node->GetM()).resize(0);
      (node->GetMmod()).resize(0);

      // reset derivative map of Mortar matrices
      (node->GetDerivD()).clear();
      (node->GetDerivM()).clear();
      
      // reset nodal weighted gap
      node->Getg() = 1.0e12;

      // reset feasible projection status
      node->HasProj() = false;
    }

    // loop over all elements to reset contact candidates / search lists
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->SearchElements().resize(0);
    }

    // reset matrix containing interface contact segments (gmsh)
    CSegs().Shape(0,0);
      
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
      IntegrateSlave(*selement);
  #endif // #ifndef CONTACTONEMORTARLOOP
    }
    
    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);
      
      // loop over the contact candidate master elements of sele_
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);
        
        //********************************************************************
        // 1) perform coupling (projection + overlap detection for sl/m pair)
        // 2) integrate Mortar matrix M and weighted gap g
        // 3) compute directional derivative of M and g and store into nodes
        //********************************************************************
        IntegrateCoupling(*selement,*melement);
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
      
      if (kcnode->Active())
      {
        // check two versions of weighted gap
        double defgap = 0.0;
        double wii = (kcnode->GetD()[0])[kcnode->Dofs()[0]];
        
        for (int j=0;j<3;++j)
          defgap-= (kcnode->n()[j])*wii*(kcnode->xspatial()[j]);
        
        vector<map<int,double> > mmap = kcnode->GetM();
        map<int,double>::iterator mcurr;
            
        for (int m=0;m<mnodefullmap_->NumMyElements();++m)
        {
          int gid = mnodefullmap_->GID(m);
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          CNode* cmnode = static_cast<CNode*>(mnode);
          const int* mdofs = cmnode->Dofs();
          bool hasentry = false;
          
          // look for this master node in M-map of the active slave node
          for (mcurr=mmap[0].begin();mcurr!=mmap[0].end();++mcurr)
            if ((mcurr->first)==mdofs[0])
            {
              hasentry=true;
              break;
            }
          
          double mik = (mmap[0])[mdofs[0]];
          double* mxi = cmnode->xspatial();
          
          // get out of here, if master node not adjacent or coupling very weak
          if (!hasentry || abs(mik)<1.0e-12) continue;
                
          for (int j=0;j<3;++j)
            defgap+= (kcnode->n()[j]) * mik * mxi[j];
        }
        
        //cout << "SNode: " << kcnode->Id() << " IntGap: " << gap << " DefGap: " << defgap << endl;
        //kcnode->Getg() = defgap;
      }
      
      // store gap-values into newG
      newG[k]=kcnode->Getg();
              
      // print results (derivatives) to screen
      if (abs(newG[k]-refG[k]) > 1e-12)
      {
        cout << "G-FD-derivative for node S" << kcnode->Id() << endl;
        //cout << "Ref-G: " << refG[k] << endl;
        //cout << "New-G: " << newG[k] << endl;
        cout << "Deriv: " << mnode->Dofs()[fd%2] << " " << (newG[k]-refG[k])/delta << endl;
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
  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

    //reset nodal normal vector
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
      node->txi()[j]=0.0;
      node->teta()[j]=0.0;
    }

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((node->GetDerivN()).size());++j)
      (node->GetDerivN())[j].clear();
    (node->GetDerivN()).resize(0);
    
    // reset derivative maps of tangent vector
    for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
      (node->GetDerivTxi())[j].clear();
    (node->GetDerivTxi()).resize(0);
        
    // reset closest node
    // (FIXME: at the moment we do not need this info. in the next
    // iteration, but it might be helpful for accelerated search!!!)
    node->ClosestNode() = -1;

    // reset nodal Mortar maps
    for (int j=0;j<(int)((node->GetD()).size());++j)
      (node->GetD())[j].clear();
    for (int j=0;j<(int)((node->GetM()).size());++j)
      (node->GetM())[j].clear();
    for (int j=0;j<(int)((node->GetMmod()).size());++j)
      (node->GetMmod())[j].clear();

    (node->GetD()).resize(0);
    (node->GetM()).resize(0);
    (node->GetMmod()).resize(0);

    // reset derivative map of Mortar matrices
    (node->GetDerivD()).clear();
    (node->GetDerivM()).clear();
    
    // reset nodal weighted gap
    node->Getg() = 1.0e12;

    // reset feasible projection status
    node->HasProj() = false;
  }

  // loop over all elements to reset contact candidates / search lists
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
    CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
    element->SearchElements().resize(0);
  }

  // reset matrix containing interface contact segments (gmsh)
  CSegs().Shape(0,0);
  
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
    IntegrateSlave(*selement);
#endif // #ifndef CONTACTONEMORTARLOOP
  }
  
  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CElement* selement = static_cast<CElement*>(ele1);
    
    // loop over the contact candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j=0;j<selement->NumSearchElements();++j)
    {
      int gid2 = selement->SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      CElement* melement = static_cast<CElement*>(ele2);
      
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pair)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      IntegrateCoupling(*selement,*melement);
    }
  }
  // *******************************************************************
  
  return;
}

/*----------------------------------------------------------------------*
 | Finite difference check for tang. LM derivatives           popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FDCheckTangLMDeriv()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;
    
  // create storage for tangential LM values
  int nrow = snoderowmap_->NumMyElements();
  vector<double> refTLM(nrow);
  vector<double> newTLM(nrow);
  
  // store reference
  // loop over proc's slave nodes
  for (int i=0; i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CNode* cnode = static_cast<CNode*>(node);
    
    double val = 0.0;
    for (int dim=0;dim<3;++dim)
      val += (cnode->txi()[dim])*(cnode->lm()[dim]);
    
    // store gap-values into refTLM
    refTLM[i]=val;
  }
  
  // global loop to apply FD scheme to all slave dofs (=2*nodes)
  for (int fd=0; fd<2*snodefullmap_->NumMyElements();++fd)
  {
    // Initialize
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

      //reset nodal normal vector
      for (int j=0;j<3;++j)
      {
        node->n()[j]=0.0;
        node->txi()[j]=0.0;
        node->teta()[j]=0.0;
      }

      // reset derivative maps of normal vector
      for (int j=0;j<(int)((node->GetDerivN()).size());++j)
        (node->GetDerivN())[j].clear();
      (node->GetDerivN()).resize(0);
      
      // reset derivative maps of tangent vector
      for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
        (node->GetDerivTxi())[j].clear();
      (node->GetDerivTxi()).resize(0);
          
      // reset closest node
      // (FIXME: at the moment we do not need this info. in the next
      // iteration, but it might be helpful for accelerated search!!!)
      node->ClosestNode() = -1;

      // reset nodal Mortar maps
      for (int j=0;j<(int)((node->GetD()).size());++j)
        (node->GetD())[j].clear();
      for (int j=0;j<(int)((node->GetM()).size());++j)
        (node->GetM())[j].clear();
      for (int j=0;j<(int)((node->GetMmod()).size());++j)
        (node->GetMmod())[j].clear();

      (node->GetD()).resize(0);
      (node->GetM()).resize(0);
      (node->GetMmod()).resize(0);

      // reset derivative map of Mortar matrices
      (node->GetDerivD()).clear();
      (node->GetDerivM()).clear();
      
      // reset nodal weighted gap
      node->Getg() = 1.0e12;

      // reset feasible projection status
      node->HasProj() = false;
    }

    // loop over all elements to reset contact candidates / search lists
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->SearchElements().resize(0);
    }

    // reset matrix containing interface contact segments (gmsh)
    CSegs().Shape(0,0);
      
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
      IntegrateSlave(*selement);
  #endif // #ifndef CONTACTONEMORTARLOOP
    }
    
    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);
      
      // loop over the contact candidate master elements of sele_
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);
        
        //********************************************************************
        // 1) perform coupling (projection + overlap detection for sl/m pair)
        // 2) integrate Mortar matrix M and weighted gap g
        // 3) compute directional derivative of M and g and store into nodes
        //********************************************************************
        IntegrateCoupling(*selement,*melement);
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
      
      double val = 0.0;
      for (int dim=0;dim<3;++dim)
        val += (kcnode->txi()[dim])*(kcnode->lm()[dim]);
      
      // store gap-values into newTLM
      newTLM[k]=val;
              
      // print results (derivatives) to screen
      if (abs(newTLM[k]-refTLM[k]) > 1e-12)
      {
        cout << "TLM-FD-derivative for node S" << kcnode->Id() << endl;
        //cout << "Ref-TLM: " << refTLM[k] << endl;
        //cout << "New-TLM: " << newTLM[k] << endl;
        cout << "Deriv: " << snode->Dofs()[fd%2] << " " << (newTLM[k]-refTLM[k])/delta << endl;
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
    // loop over all nodes to reset normals, closestnode and Mortar maps
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColNodes();++i)
    {
      CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

      //reset nodal normal vector
      for (int j=0;j<3;++j)
      {
        node->n()[j]=0.0;
        node->txi()[j]=0.0;
        node->teta()[j]=0.0;
      }

      // reset derivative maps of normal vector
      for (int j=0;j<(int)((node->GetDerivN()).size());++j)
        (node->GetDerivN())[j].clear();
      (node->GetDerivN()).resize(0);
      
      // reset derivative maps of tangent vector
      for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
        (node->GetDerivTxi())[j].clear();
      (node->GetDerivTxi()).resize(0);
          
      // reset closest node
      // (FIXME: at the moment we do not need this info. in the next
      // iteration, but it might be helpful for accelerated search!!!)
      node->ClosestNode() = -1;

      // reset nodal Mortar maps
      for (int j=0;j<(int)((node->GetD()).size());++j)
        (node->GetD())[j].clear();
      for (int j=0;j<(int)((node->GetM()).size());++j)
        (node->GetM())[j].clear();
      for (int j=0;j<(int)((node->GetMmod()).size());++j)
        (node->GetMmod())[j].clear();

      (node->GetD()).resize(0);
      (node->GetM()).resize(0);
      (node->GetMmod()).resize(0);

      // reset derivative map of Mortar matrices
      (node->GetDerivD()).clear();
      (node->GetDerivM()).clear();
      
      // reset nodal weighted gap
      node->Getg() = 1.0e12;

      // reset feasible projection status
      node->HasProj() = false;
    }

    // loop over all elements to reset contact candidates / search lists
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
      element->SearchElements().resize(0);
    }

    // reset matrix containing interface contact segments (gmsh)
    CSegs().Shape(0,0);
      
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
      IntegrateSlave(*selement);
  #endif // #ifndef CONTACTONEMORTARLOOP
    }
    
    // loop over proc's slave elements of the interface for integration
    // use standard column map to include processor's ghosted elements
    for (int i=0; i<selecolmap_->NumMyElements();++i)
    {
      int gid1 = selecolmap_->GID(i);
      DRT::Element* ele1 = idiscret_->gElement(gid1);
      if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
      CElement* selement = static_cast<CElement*>(ele1);
      
      // loop over the contact candidate master elements of sele_
      // use slave element's candidate list SearchElements !!!
      for (int j=0;j<selement->NumSearchElements();++j)
      {
        int gid2 = selement->SearchElements()[j];
        DRT::Element* ele2 = idiscret_->gElement(gid2);
        if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
        CElement* melement = static_cast<CElement*>(ele2);
        
        //********************************************************************
        // 1) perform coupling (projection + overlap detection for sl/m pair)
        // 2) integrate Mortar matrix M and weighted gap g
        // 3) compute directional derivative of M and g and store into nodes
        //********************************************************************
        IntegrateCoupling(*selement,*melement);
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
      
      double val = 0.0;
      for (int dim=0;dim<3;++dim)
        val += (kcnode->txi()[dim])*(kcnode->lm()[dim]);
      
      // store gap-values into newTLM
      newTLM[k]=val;
              
      // print results (derivatives) to screen
      if (abs(newTLM[k]-refTLM[k]) > 1e-12)
      {
        cout << "TLM-FD-derivative for node S" << kcnode->Id() << endl;
        //cout << "Ref-TLM: " << refTLM[k] << endl;
        //cout << "New-TLM: " << newTLM[k] << endl;
        cout << "Deriv: " << mnode->Dofs()[fd%2] << " " << (newTLM[k]-refTLM[k])/delta << endl;
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
  // loop over all nodes to reset normals, closestnode and Mortar maps
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CNode* node = static_cast<CONTACT::CNode*>(idiscret_->lColNode(i));

    //reset nodal normal vector
    for (int j=0;j<3;++j)
    {
      node->n()[j]=0.0;
      node->txi()[j]=0.0;
      node->teta()[j]=0.0;
    }

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((node->GetDerivN()).size());++j)
      (node->GetDerivN())[j].clear();
    (node->GetDerivN()).resize(0);
    
    // reset derivative maps of tangent vector
    for (int j=0;j<(int)((node->GetDerivTxi()).size());++j)
      (node->GetDerivTxi())[j].clear();
    (node->GetDerivTxi()).resize(0);
        
    // reset closest node
    // (FIXME: at the moment we do not need this info. in the next
    // iteration, but it might be helpful for accelerated search!!!)
    node->ClosestNode() = -1;

    // reset nodal Mortar maps
    for (int j=0;j<(int)((node->GetD()).size());++j)
      (node->GetD())[j].clear();
    for (int j=0;j<(int)((node->GetM()).size());++j)
      (node->GetM())[j].clear();
    for (int j=0;j<(int)((node->GetMmod()).size());++j)
      (node->GetMmod())[j].clear();

    (node->GetD()).resize(0);
    (node->GetM()).resize(0);
    (node->GetMmod()).resize(0);

    // reset derivative map of Mortar matrices
    (node->GetDerivD()).clear();
    (node->GetDerivM()).clear();
    
    // reset nodal weighted gap
    node->Getg() = 1.0e12;

    // reset feasible projection status
    node->HasProj() = false;
  }

  // loop over all elements to reset contact candidates / search lists
  // (use fully overlapping column map)
  for (int i=0;i<idiscret_->NumMyColElements();++i)
  {
    CONTACT::CElement* element = static_cast<CONTACT::CElement*>(idiscret_->lColElement(i));
    element->SearchElements().resize(0);
  }

  // reset matrix containing interface contact segments (gmsh)
  CSegs().Shape(0,0);
  
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
    IntegrateSlave(*selement);
#endif // #ifndef CONTACTONEMORTARLOOP
  }
  
  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CElement* selement = static_cast<CElement*>(ele1);
    
    // loop over the contact candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j=0;j<selement->NumSearchElements();++j)
    {
      int gid2 = selement->SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      CElement* melement = static_cast<CElement*>(ele2);
      
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pair)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      IntegrateCoupling(*selement,*melement);
    }
  }
  // *******************************************************************
  
  return;
}

#endif  // #ifdef CCADISCRET
