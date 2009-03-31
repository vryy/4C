/*!----------------------------------------------------------------------
\file drt_contact_coupling_tools.cpp
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
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_coupling3d.H"
#include "drt_contact_projector.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
CONTACT::Coupling3d::Coupling3d(DRT::Discretization& idiscret, int dim, bool quad,
                                bool auxplane,
                                CONTACT::CElement& sele, CONTACT::CElement& mele,
                                vector<vector<double> >& testv,
                                bool printderiv) :
idiscret_(idiscret),
dim_(dim),
quad_(quad),
auxplane_(auxplane),
sele_(sele),
mele_(mele)
{
  // *********************************************************************
  // the three-dimensional case
  // *********************************************************************
  // check for quadratic elements
  if (sele.Shape()!=DRT::Element::tri3 && sele.Shape()!=DRT::Element::quad4)
    dserror("ERROR: FD check for 3D mortar coupling not yet impl. for quadratic elements");
  if (mele.Shape()!=DRT::Element::tri3 && mele.Shape()!=DRT::Element::quad4)
    dserror("ERROR: FD check for 3D mortar coupling not yet impl. for quadratic elements");

  // rough check whether elements are "near"
  bool near = RoughCheck();
  if (!near) return;
  
  // map to store projection parameter alpha for each master node
  map<int,double> projpar;
  
  // *******************************************************************
  // ************ Coupling with or without auxiliary plane *************
  // *******************************************************************
  if (CouplingInAuxPlane())
  {
    // compute auxiliary plane for 3D coupling
    AuxiliaryPlane();
    
    // project slave element nodes onto auxiliary plane
    ProjectSlave();
    
    // project master element nodes onto auxiliary plane
    ProjectMaster();
    
    // tolerance for polygon clipping
    double sminedge = sele.MinEdgeSize();
    double mminedge = mele.MinEdgeSize(); 
    double tol = CONTACTCLIPTOL * min(sminedge,mminedge);
    
    
    // do clipping in auxiliary plane
    PolygonClipping(SlaveVertices(),MasterVertices(),Clip(),tol);
  }
  // *******************************************************************
  else //(!CouplingInAuxPlane())
  {
    // get some data
    int nsnodes = SlaveElement().NumNode();
    int nmnodes = MasterElement().NumNode();
    
    // get slave vertices in slave element parameter space (direct)
    // additionally get slave vertex Ids for later linearization
    vector<vector<double> > svertices(nsnodes,vector<double>(3));
    vector<int> snodeids(1);
    
    for (int i=0;i<nsnodes;++i)
    {     
      double xi[2] = {0.0, 0.0};
      SlaveElement().LocalCoordinatesOfNode(i,xi);
      svertices[i][0] = xi[0];
      svertices[i][1] = xi[1];
      svertices[i][2] = 0.0;
   
      // relevant ids (here only slave node id itself)
      snodeids[0] = SlaveElement().NodeIds()[i];
      
      // store into vertex data structure
      SlaveVertices().push_back(Vertex(svertices[i],Vertex::slave,snodeids,NULL,NULL,false,false,NULL,-1.0));
    }
    
    // get master vertices in slave element parameter space (project)
    // additionally get master vertex Ids for later linearization
    vector<vector<double> > mvertices(nmnodes,vector<double>(3));
    vector<int> mnodeids(1);
    for (int i=0;i<nmnodes;++i)
    {
      int gid = MasterElement().NodeIds()[i];
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* mnode = static_cast<CNode*>(node);
      
      // do the projection
      double sxi[2] = {0.0, 0.0};
      double alpha = 0.0;
      CONTACT::Projector projector(3);
      //cout << "Projecting master node ID: " << mnode->Id() << endl;
      projector.ProjectElementNormal3D(*mnode,SlaveElement(),sxi,alpha);
      
      mvertices[i][0] = sxi[0];
      mvertices[i][1] = sxi[1];
      mvertices[i][2] = 0.0;
      
      // relevant ids (here only master node id itself)
      mnodeids[0] = gid;
      
      // store proj. parameter for later linearization
      projpar[gid] = alpha;
      
      // store into vertex data structure
      MasterVertices().push_back(Vertex(mvertices[i],Vertex::projmaster,mnodeids,NULL,NULL,false,false,NULL,-1.0));
    }
    
    // normal is (0,0,1) in slave element parameter space
    Auxn()[0] = 0.0; Auxn()[1] = 0.0; Auxn()[2] = 1.0;
    Lauxn() = 1.0;
    
    // tolerance for polygon clipping
    // minimum edge size in parameter space is 1
    double tol = CONTACTCLIPTOL;
      
    // do clipping in slave element parameter space
    PolygonClipping(SlaveVertices(),MasterVertices(),Clip(),tol);
  }
  // *******************************************************************
  
  // proceed only if clipping polygon is at least a triangle
  int clipsize = (int)(Clip().size());
  bool overlap = false;
  if (clipsize>=3) overlap = true;
  if (overlap)
  {
    // fill testv vector for FDVertex check
    for (int k=0;k<clipsize;++k)
    {
      vector<double> testventry(6);
      testventry[0] = (Clip()[k]).Coord()[0];
      testventry[1] = (Clip()[k]).Coord()[1];
      testventry[2] = (Clip()[k]).Coord()[2];
      testventry[3] = SlaveElement().Id();
      testventry[4] = MasterElement().Id();
      testventry[5] = (Clip()[k]).VType();
      testv.push_back(testventry);
    }
    
    // check / set  projection status of slave nodes
    HasProjStatus();
    
    // do linearization + triangulation of clip polygon
    Triangulation(projpar,testv,printderiv);
    
    // do integration of integration cells
    IntegrateCells();
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
CONTACT::Coupling3d::Coupling3d(DRT::Discretization& idiscret, int dim, bool quad,
                                bool auxplane,
                                CONTACT::CElement& sele, CONTACT::CElement& mele,
                                vector<vector<double> >& testgps,
                                vector<vector<double> >& testgpm,
                                vector<vector<double> >& testjs,
                                vector<vector<double> >& testji,
                                bool printderiv) :
idiscret_(idiscret),
dim_(dim),
quad_(quad),
auxplane_(auxplane),
sele_(sele),
mele_(mele)
{
  // *********************************************************************
  // the three-dimensional case
  // *********************************************************************
  // check for quadratic elements
  if (sele.Shape()!=DRT::Element::tri3 && sele.Shape()!=DRT::Element::quad4)
    dserror("ERROR: FD check for 3D mortar coupling not yet impl. for quadratic elements");
  if (mele.Shape()!=DRT::Element::tri3 && mele.Shape()!=DRT::Element::quad4)
    dserror("ERROR: FD check for 3D mortar coupling not yet impl. for quadratic elements");

  // rough check whether elements are "near"
  bool near = RoughCheck();
  if (!near) return;
  
  // map to store projection parameter alpha for each master node
  map<int,double> projpar;
  
  // *******************************************************************
  // ************ Coupling with or without auxiliary plane *************
  // *******************************************************************
  if (CouplingInAuxPlane())
  {
    // compute auxiliary plane for 3D coupling
    AuxiliaryPlane();
    
    // project slave element nodes onto auxiliary plane
    ProjectSlave();
    
    // project master element nodes onto auxiliary plane
    ProjectMaster();
    
    // tolerance for polygon clipping
    double sminedge = sele.MinEdgeSize();
    double mminedge = mele.MinEdgeSize(); 
    double tol = CONTACTCLIPTOL * min(sminedge,mminedge);
    
    // do clipping in auxiliary plane
    PolygonClipping(SlaveVertices(),MasterVertices(),Clip(),tol);
  }
  // *******************************************************************
  else //(!CouplingInAuxPlane())
  {
    // get some data
    int nsnodes = SlaveElement().NumNode();
    int nmnodes = MasterElement().NumNode();
    
    // get slave vertices in slave element parameter space (direct)
    // additionally get slave vertex Ids for later linearization
    vector<vector<double> > svertices(nsnodes,vector<double>(3));
    vector<int> snodeids(1);
    
    for (int i=0;i<nsnodes;++i)
    {     
      double xi[2] = {0.0, 0.0};
      SlaveElement().LocalCoordinatesOfNode(i,xi);
      svertices[i][0] = xi[0];
      svertices[i][1] = xi[1];
      svertices[i][2] = 0.0;
   
      // relevant ids (here only slave node id itself)
      snodeids[0] = SlaveElement().NodeIds()[i];
      
      // store into vertex data structure
      SlaveVertices().push_back(Vertex(svertices[i],Vertex::slave,snodeids,NULL,NULL,false,false,NULL,-1.0));
    }
    
    // get master vertices in slave element parameter space (project)
    // additionally get master vertex Ids for later linearization
    vector<vector<double> > mvertices(nmnodes,vector<double>(3));
    vector<int> mnodeids(1);
    for (int i=0;i<nmnodes;++i)
    {
      int gid = MasterElement().NodeIds()[i];
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CNode* mnode = static_cast<CNode*>(node);
      
      // do the projection
      // the third component of sxi will be the proj. parameter alpha!
      double sxi[2] = {0.0, 0.0};
      double alpha = 0.0;
      CONTACT::Projector projector(3);
      //cout << "Projecting master node ID: " << mnode->Id() << endl;
      projector.ProjectElementNormal3D(*mnode,SlaveElement(),sxi,alpha);
      
      mvertices[i][0] = sxi[0];
      mvertices[i][1] = sxi[1];
      mvertices[i][2] = 0.0;
      
      // relevant ids (here only master node id itself)
      mnodeids[0] = gid;
      
      // store proj. parameter for later linearization
      projpar[gid] = alpha;
      
      // store into vertex data structure
      MasterVertices().push_back(Vertex(mvertices[i],Vertex::projmaster,mnodeids,NULL,NULL,false,false,NULL,-1.0));
    }
    
    // normal is (0,0,1) in slave element parameter space
    Auxn()[0] = 0.0; Auxn()[1] = 0.0; Auxn()[2] = 1.0;
    
    // tolerance for polygon clipping
    // minimum edge size in parameter space is 1
    double tol = CONTACTCLIPTOL;
      
    // do clipping in slave element parameter space
    PolygonClipping(SlaveVertices(),MasterVertices(),Clip(),tol);
  }
  // *******************************************************************
  
  // proceed only if clipping polygon is at least a triangle
  int clipsize = (int)(Clip().size());
  bool overlap = false;
  if (clipsize>=3) overlap = true;
  if (overlap)
  {
    // check / set  projection status of slave nodes
    HasProjStatus();
    
    // do linearization + triangulation of clip polygon
    Triangulation(projpar);
    
    // do integration of integration cells
    IntegrateCells(testgps,testgpm,testjs,testji,printderiv);
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon (3D)                        popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::Triangulation(map<int,double>& projpar,
                                        vector<vector<double> >& testv, bool printderiv)
{
  // preparations
  int clipsize = (int)(Clip().size());
  vector<vector<map<int,double> > > linvertex(clipsize,vector<map<int,double> >(3));
  vector<map<int,double> > lincenter(3);
  
  //**********************************************************************
  // (1) Linearization of clip vertex coordinates
  //**********************************************************************
  VertexLinearization(linvertex,projpar,printderiv);
  
  //**********************************************************************
  // (2) Find center of clipping polygon (centroid formulas)
  //**********************************************************************
  vector<double> clipcenter(3);
  for (int k=0;k<3;++k) clipcenter[k] = 0.0;
  double fac = 0.0;
  
  // first we need node averaged center
  double nac[3] = {0.0, 0.0, 0.0};
  for (int i=0;i<clipsize;++i)
    for (int k=0;k<3;++k)
      nac[k] += (Clip()[i].Coord()[k] / clipsize);
  
  // loop over all triangles of polygon
  for (int i=0; i<clipsize; ++i)
  {
    double xi_i[3] = {0.0, 0.0, 0.0};
    double xi_ip1[3] = {0.0, 0.0, 0.0};
    
    // standard case    
    if (i<clipsize-1)
    {
      for (int k=0;k<3;++k) xi_i[k] = Clip()[i].Coord()[k];
      for (int k=0;k<3;++k) xi_ip1[k] = Clip()[i+1].Coord()[k];
    }
    // last vertex of clip polygon
    else
    {
      for (int k=0;k<3;++k) xi_i[k] = Clip()[clipsize-1].Coord()[k];
      for (int k=0;k<3;++k) xi_ip1[k] = Clip()[0].Coord()[k];
    }

    // triangle area
    double diff1[3] = {0.0, 0.0, 0.0};
    double diff2[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k) diff1[k] = xi_ip1[k] - xi_i[k];
    for (int k=0;k<3;++k) diff2[k] = xi_i[k] - nac[k];
    
    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = diff1[1]*diff2[2] - diff1[2]*diff2[1];
    cross[1] = diff1[2]*diff2[0] - diff1[0]*diff2[2];
    cross[2] = diff1[0]*diff2[1] - diff1[1]*diff2[0];
    
    double Atri = 0.5 * sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
    
    // add contributions to clipcenter and fac
    fac += Atri;
    for (int k=0;k<3;++k) clipcenter[k] += 1.0/3.0 * (xi_i[k] + xi_ip1[k] + nac[k]) * Atri;
  }
  
  //final clipcenter
  for (int k=0;k<3;++k) clipcenter[k] /= fac; 
  //cout << "Clipcenter: " << clipcenter[0] << " " << clipcenter[1] << " " << clipcenter[2] << endl;
    
  // fill testv vector for FDVertex check with center
  // we use 3 for VType (because 0,1,2 are already in use)
  vector<double> testventry(6);
  testventry[0] = clipcenter[0];
  testventry[1] = clipcenter[1];
  testventry[2] = clipcenter[2];
  testventry[3] = SlaveElement().Id();
  testventry[4] = MasterElement().Id();
  testventry[5] = 3;
  testv.push_back(testventry);
  
  //**********************************************************************
  // (3) Linearization of clip center coordinates
  //**********************************************************************
  CenterLinearization(linvertex,lincenter);
  
  if (printderiv)
  {
    static int counter = 0;
    typedef map<int,double>::const_iterator CI;
    for (int i=0;i<clipsize;++i)
    {
      cout << "\nAnalytical derivative for Vertex " << counter << " (x-component). VType=" << Clip()[i].VType() << endl;
      for (CI p=linvertex[i][0].begin();p!=linvertex[i][0].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      cout << "Analytical derivative for Vertex " << counter << " (y-component). VType=" << Clip()[i].VType() << endl;
      for (CI p=linvertex[i][1].begin();p!=linvertex[i][1].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      cout << "Analytical derivative for Vertex " << counter << " (z-component). VType=" << Clip()[i].VType() << endl;
      for (CI p=linvertex[i][2].begin();p!=linvertex[i][2].end();++p)
        cout << "Dof: " << p->first << "\t" << p->second << endl;
      counter = counter+1;
    }
    cout << "\nAnalytical derivative for Vertex " << counter << " (x-component). VType=" << 3 << endl;
    for (CI p=lincenter[0].begin();p!=lincenter[0].end();++p)
      cout << "Dof: " << p->first << "\t" << p->second << endl;
    cout << "Analytical derivative for Vertex " << counter << " (y-component). VType=" << 3 << endl;
    for (CI p=lincenter[1].begin();p!=lincenter[1].end();++p)
      cout << "Dof: " << p->first << "\t" << p->second << endl;
    cout << "Analytical derivative for Vertex " << counter << " (z-component). VType=" << 3 << endl;
    for (CI p=lincenter[2].begin();p!=lincenter[2].end();++p)
      cout << "Dof: " << p->first << "\t" << p->second << endl;
    counter = counter+1;
  }
  
  //**********************************************************************
  // (4)Triangulation -> Intcells
  //**********************************************************************
  vector<Intcell> cells;
  
  // easy if clip polygon = triangle: 1 Intcell
  if (clipsize==3)
  {
    // Intcell vertices = clip polygon vertices
    Epetra_SerialDenseMatrix coords(3,clipsize);
    for (int i=0;i<clipsize;++i)
      for (int k=0;k<3;++k)
        coords(k,i) = Clip()[i].Coord()[k];
    
    // create Intcell object and push back
    Cells().push_back(rcp(new Intcell(0,3,coords,Auxn(),DRT::Element::tri3,
      CouplingInAuxPlane(),linvertex[0],linvertex[1],linvertex[2],GetDerivAuxn())));
    
  }
  
  // triangulation if clip polygon > triangle
  else
  {
    // No. of Intcells is equal to no. of clip polygon vertices
    for (int num=0;num<clipsize;++num)
    {
      // the first vertex is always the clip center
      // the second vertex is always the current clip vertex
      Epetra_SerialDenseMatrix coords(3,3);
      for (int k=0;k<3;++k)
      {
        coords(k,0) = clipcenter[k];
        coords(k,1) = Clip()[num].Coord()[k];
      }
      
      // the third vertex is the next vertex on clip polygon
      int numplus1 = num+1;
      if (num==clipsize-1)
      {
        for (int k=0;k<3;++k) coords(k,2) = Clip()[0].Coord()[k];
        numplus1 = 0;
      }
      else
        for (int k=0;k<3;++k) coords(k,2) = Clip()[num+1].Coord()[k];
      
      // create Intcell object and push back
      Cells().push_back(rcp(new Intcell(num,3,coords,Auxn(),DRT::Element::tri3,
        CouplingInAuxPlane(),lincenter,linvertex[num],linvertex[numplus1],GetDerivAuxn())));
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling3d::IntegrateCells(vector<vector<double> >& testgps,
                                         vector<vector<double> >& testgpm,
                                         vector<vector<double> >& testjs,
                                         vector<vector<double> >& testji,
                                         bool printderiv)
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Integrate the Mortar matrix M and the weighted gap function g~ on  */
  /* the current integration cell of the slave / master CElement pair   */
  /**********************************************************************/
  
  // create an integrator instance with correct NumGP and Dim
  // it is sufficient to do this once as all Intcells are triangles
  CONTACT::Integrator integrator(Cells()[0]->Shape());
    
  // loop over all integration cells
  for (int i=0;i<(int)(Cells().size());++i)
  {
    // compare intcell area with slave element area
    double intcellarea = Cells()[i]->Area();
    double selearea = 0.0;
    DRT::Element::DiscretizationType dt = SlaveElement().Shape();
    if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
      selearea = 4.0;
    else if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
      selearea = 0.5;
    else dserror("ERROR: IntegrateCells: Invalid 3D slave element type");
    
    // integrate cell only if not neglectable
    if (intcellarea<CONTACTINTLIM*selearea) continue;
    
    // do the two integrations
    // *******************************************************************
    // ************ Coupling with or without auxiliary plane *************
    // *******************************************************************
    int nrow = sele_.NumNode();
    int ncol = mele_.NumNode();
    RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),ncol*Dim()));
    RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));
    
    if (CouplingInAuxPlane())
      integrator.IntegrateDerivCell3DAuxPlane(sele_,mele_,Cells()[i],Auxn(),mseg,gseg,testgps,testgpm,testjs,testji,printderiv);
    else //(!CouplingInAuxPlane())
      integrator.IntegrateDerivCell3D(sele_,mele_,Cells()[i],mseg,gseg,testgps,testgpm,testjs,testji,printderiv);
    // *******************************************************************
    
    // do the two assemblies into the slave nodes
    // if CONTACTONEMORTARLOOP defined, then AssembleM does M AND D matrices !!!
    integrator.AssembleM(Comm(),sele_,mele_,*mseg);
    integrator.AssembleG(Comm(),sele_,*gseg);
  }
  return true;
}

#endif //#ifdef CCADISCRET
