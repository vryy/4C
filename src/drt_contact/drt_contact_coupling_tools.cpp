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

#include "drt_contact_coupling.H"
#include "drt_contact_projector.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
CONTACT::Coupling::Coupling(DRT::Discretization& idiscret,
                            CONTACT::CElement& sele, CONTACT::CElement& mele,
                            int dim, Epetra_SerialDenseMatrix& csegs,
                            vector<vector<double> >& testv,
                            bool printderiv) :
idiscret_(idiscret),
sele_(sele),
mele_(mele),
dim_(dim),
contactsegs_(csegs)
{
  // check sanity of dimension
  if (Dim()!=2 && Dim()!=3) dserror("Dim. must be 2D or 3D!");
  
  // *********************************************************************
  // the two-dimensional case
  // *********************************************************************
  if (Dim()==2)
  {
    // prepare overlap integration
    vector<bool> hasproj(4);
    vector<double> xiproj(4);
    bool overlap = false;
  
    // project the element pair
    Project2D(hasproj,xiproj);
  
    // check for element overlap
    overlap = DetectOverlap2D(hasproj,xiproj);
  
    // integrate the element overlap
    if (overlap) IntegrateOverlap2D(xiproj);
  }
  
  // *********************************************************************
  // the three-dimensional case
  // *********************************************************************
  else if (Dim()==3)
  {
    // check for quadratic elements
    if (sele.Shape()!=DRT::Element::tri3 && sele.Shape()!=DRT::Element::quad4)
      dserror("ERROR: 3D mortar coupling not yet impl. for quadratic elements");
    if (mele.Shape()!=DRT::Element::tri3 && mele.Shape()!=DRT::Element::quad4)
      dserror("ERROR: 3D mortar coupling not yet impl. for quadratic elements");

    // rough check whether elements are "near"
    bool near = RoughCheck3D();
    if (!near) return;
    
    // map to store projection parameter alpha for each master node
    map<int,double> projpar;
    
    // *******************************************************************
    // ************ Coupling with or without auxiliary plane *************
    // *******************************************************************
#ifdef CONTACTAUXPLANE
    // compute auxiliary plane for 3D coupling
    AuxiliaryPlane3D();
    
    // project slave element nodes onto auxiliary plane
    ProjectSlave3D();
    
    // project master element nodes onto auxiliary plane
    ProjectMaster3D();
    
    // tolerance for polygon clipping
    double sminedge = sele.MinEdgeSize();
    double mminedge = mele.MinEdgeSize(); 
    double tol = CONTACTCLIPTOL * min(sminedge,mminedge);
    
    
    // do clipping in auxiliary plane
    PolygonClipping3D(SlaveVertices(),MasterVertices(),Clip(),tol);
    int clipsize = (int)(Clip().size());
#else
    
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
    
    // tolerance for polygon clipping
    // minimum edge size in parameter space is 1
    double tol = CONTACTCLIPTOL;
      
    // do clipping in slave element parameter space
    PolygonClipping3D(SlaveVertices(),MasterVertices(),Clip(),tol);
    int clipsize = (int)(Clip().size());
#endif // #ifdef CONTACTAUXPLANE
    // *******************************************************************
    
    // proceed only if clipping polygon is at least a triangle
    bool overlap = false;
    if (clipsize>=3) overlap = true;
    if (overlap)
    {
      // fill testv vector for FDVertex check
      for (int k=0;k<clipsize;++k)
      {
        vector<double> testventry(5);
        testventry[0] = (Clip()[k]).Coord()[0];
        testventry[1] = (Clip()[k]).Coord()[1];
        testventry[2] = SlaveElement().Id();
        testventry[3] = MasterElement().Id();
        testventry[4] = (Clip()[k]).VType();
        testv.push_back(testventry);
      }
      // check / set  projection status of slave nodes
      HasProjStatus3D();
      
      // do linearization + triangulation of clip polygon
      Triangulation3D(projpar,testv,printderiv);
      
      // do integration of integration cells
      IntegrateCells3D();
    }
  }
  
  // *********************************************************************
  // invalid cases for dim (other than 2D or 3D)
  // *********************************************************************
  else dserror("ERROR: Coupling can only be called for 2D or 3D!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
CONTACT::Coupling::Coupling(DRT::Discretization& idiscret,
                            CONTACT::CElement& sele, CONTACT::CElement& mele,
                            int dim, Epetra_SerialDenseMatrix& csegs,
                            vector<vector<double> >& testgps,
                            vector<vector<double> >& testgpm,
                            vector<vector<double> >& testjs,
                            vector<vector<double> >& testji,
                            bool printderiv) :
idiscret_(idiscret),
sele_(sele),
mele_(mele),
dim_(dim),
contactsegs_(csegs)
{
  // check sanity of dimension
  if (Dim()!=2 && Dim()!=3) dserror("Dim. must be 2D or 3D!");
  
  // *********************************************************************
  // the two-dimensional case
  // *********************************************************************
  if (Dim()==2)
  {
    // prepare overlap integration
    vector<bool> hasproj(4);
    vector<double> xiproj(4);
    bool overlap = false;
  
    // project the element pair
    Project2D(hasproj,xiproj);
  
    // check for element overlap
    overlap = DetectOverlap2D(hasproj,xiproj);
  
    // integrate the element overlap
    if (overlap) IntegrateOverlap2D(xiproj);
  }
  
  // *********************************************************************
  // the three-dimensional case
  // *********************************************************************
  else if (Dim()==3)
  {
    // check for quadratic elements
    if (sele.Shape()!=DRT::Element::tri3 && sele.Shape()!=DRT::Element::quad4)
      dserror("ERROR: 3D mortar coupling not yet impl. for quadratic elements");
    if (mele.Shape()!=DRT::Element::tri3 && mele.Shape()!=DRT::Element::quad4)
      dserror("ERROR: 3D mortar coupling not yet impl. for quadratic elements");

    // rough check whether elements are "near"
    bool near = RoughCheck3D();
    if (!near) return;
    
    // map to store projection parameter alpha for each master node
    map<int,double> projpar;
    
    // *******************************************************************
    // ************ Coupling with or without auxiliary plane *************
    // *******************************************************************
#ifdef CONTACTAUXPLANE
    // compute auxiliary plane for 3D coupling
    AuxiliaryPlane3D();
    
    // project slave element nodes onto auxiliary plane
    ProjectSlave3D();
    
    // project master element nodes onto auxiliary plane
    ProjectMaster3D();
    
    // tolerance for polygon clipping
    double sminedge = sele.MinEdgeSize();
    double mminedge = mele.MinEdgeSize(); 
    double tol = CONTACTCLIPTOL * min(sminedge,mminedge);
    
    // do clipping in auxiliary plane
    PolygonClipping3D(SlaveVertices(),MasterVertices(),Clip(),tol);
    int clipsize = (int)(Clip().size());
#else
    
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
    PolygonClipping3D(SlaveVertices(),MasterVertices(),Clip(),tol);
    int clipsize = (int)(Clip().size());
#endif // #ifdef CONTACTAUXPLANE
    // *******************************************************************
    
    // proceed only if clipping polygon is at least a triangle
    bool overlap = false;
    if (clipsize>=3) overlap = true;
    if (overlap)
    {
      // check / set  projection status of slave nodes
      HasProjStatus3D();
      
      // do linearization + triangulation of clip polygon
      Triangulation3D(projpar);
      
      // do integration of integration cells
      IntegrateCells3D(testgps,testgpm,testjs,testji,printderiv);
    }
  }
  
  // *********************************************************************
  // invalid cases for dim (other than 2D or 3D)
  // *********************************************************************
  else dserror("ERROR: Coupling can only be called for 2D or 3D!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon (3D)                        popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::Triangulation3D(map<int,double>& projpar,
                                        vector<vector<double> >& testv, bool printderiv)
{
  // preparations
  vector<double> clipcenter(3);
  for (int k=0;k<3;++k) clipcenter[k] = 0.0;
  int clipsize = (int)(Clip().size());
  vector<vector<map<int,double> > > linvertex(clipsize,vector<map<int,double> >(3));
  vector<map<int,double> > lincenter(3);
  
#ifndef CONTACTAUXPLANE
  //**********************************************************************
  // (1) Linearization of clip vertex coordinates
  //**********************************************************************
  VertexLinearization3D(linvertex,projpar,printderiv);
#endif // #ifndef CONTACTAUXPLANE
  
  //**********************************************************************
  // (2) Find center of clipping polygon
  //**********************************************************************
  
  // ************* Coupling with or without auxiliary plane **************
#ifdef CONTACTAUXPLANE
  // as a first shot we use simple node averaging here
  // arithmetic center, NOT geometric center
  for (int i=0;i<clipsize;++i)
    for (int k=0;k<3;++k)
      clipcenter[k] += (Clip()[i].Coord()[k] / clipsize);
  //cout << "Clipcenter (simple): " << clipcenter[0] << " " << clipcenter[1] << " " << clipcenter[2] << endl;
  
  // fill testv vector for FDVertex check with center
  // we use 3 for VType (because 0,1,2 are already in use)
  vector<double> testventry(5);
  testventry[0] = clipcenter[0];
  testventry[1] = clipcenter[1];
  testventry[2] = SlaveElement().Id();
  testventry[3] = MasterElement().Id();
  testventry[4] = 3;
  testv.push_back(testventry);
#else
  // as an improved version we use centroid formulas here
  // this really yields the geometric center
  // (taken from: Moertel package, Trilinos, Sandia NL)
  double A = 0.0;
  
  for (int i=0; i<clipsize; ++i)
  {
    // check for 2D clip polygon in parameter space
    if (abs(Clip()[i].Coord()[2])>1.0e-12) dserror("ERROR: Clip polygon point with z!=0");
    double xi_i[2] = {0.0, 0.0};
    double xi_ip1[2] = {0.0, 0.0};
    
    // standard case    
    if (i<clipsize-1)
    {
      xi_i[0] = Clip()[i].Coord()[0]; xi_i[1] = Clip()[i].Coord()[1];
      xi_ip1[0] = Clip()[i+1].Coord()[0]; xi_ip1[1] = Clip()[i+1].Coord()[1];
    }
    // last vertex of clip polygon
    else
    {
      xi_i[0] = Clip()[clipsize-1].Coord()[0]; xi_i[1] = Clip()[clipsize-1].Coord()[1];
      xi_ip1[0] = Clip()[0].Coord()[0]; xi_ip1[1] = Clip()[0].Coord()[1];
    }
    
    // add contribution to area and centroid coords
    A     += xi_ip1[0]*xi_i[1] - xi_i[0]*xi_ip1[1];
    clipcenter[0] += (xi_i[0]+xi_ip1[0])*(xi_ip1[0]*xi_i[1]-xi_i[0]*xi_ip1[1]);
    clipcenter[1] += (xi_i[1]+xi_ip1[1])*(xi_ip1[0]*xi_i[1]-xi_i[0]*xi_ip1[1]);
  }
  
  // final centroid coords
  clipcenter[0] /= (3.0*A);
  clipcenter[1] /= (3.0*A);
  //cout << "Clipcenter (centroid): " << clipcenter[0] << " " << clipcenter[1] << " " << clipcenter[2] << endl;
  
  // fill testv vector for FDVertex check with center
  // we use 3 for VType (because 0,1,2 are already in use)
  vector<double> testventry(5);
  testventry[0] = clipcenter[0];
  testventry[1] = clipcenter[1];
  testventry[2] = SlaveElement().Id();
  testventry[3] = MasterElement().Id();
  testventry[4] = 3;
  testv.push_back(testventry);
        
#endif // #ifdef CONTACTAUXPLANE
  // *********************************************************************
  
#ifndef CONTACTAUXPLANE
  //**********************************************************************
  // (3) Linearization of clip center coordinates
  //**********************************************************************
  CenterLinearization3D(linvertex,lincenter);
  
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
      counter = counter+1;
    }
    cout << "\nAnalytical derivative for Vertex " << counter << " (x-component). VType=" << 3 << endl;
    for (CI p=lincenter[0].begin();p!=lincenter[0].end();++p)
      cout << "Dof: " << p->first << "\t" << p->second << endl;
    cout << "Analytical derivative for Vertex " << counter << " (y-component). VType=" << 3 << endl;
    for (CI p=lincenter[1].begin();p!=lincenter[1].end();++p)
      cout << "Dof: " << p->first << "\t" << p->second << endl;
    counter = counter+1;
  }
#endif // #ifndef CONTACTAUXPLANE
  
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
    Cells().push_back(rcp(new Intcell(0,3,coords,DRT::Element::tri3,
                      linvertex[0],linvertex[1],linvertex[2])));
    
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
      Cells().push_back(rcp(new Intcell(num,3,coords,DRT::Element::tri3,
                        lincenter,linvertex[num],linvertex[numplus1])));
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells (3D)                                 popp 02/09|
 |  THIS IS A PURE FINITE DIFFERENCE VERSION!!!                         |
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::IntegrateCells3D(vector<vector<double> >& testgps,
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
    else dserror("ERROR: IntegrateCells3D: Invalid 3D slave element type");
    
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
#ifdef CONTACTAUXPLANE
    integrator.IntegrateDerivCell3DAuxPlane(sele_,mele_,Cells()[i],Auxn(),mseg,gseg,testgps,testgpm,testjs,testji,printderiv);
#else
    integrator.IntegrateDerivCell3D(sele_,mele_,Cells()[i],mseg,gseg,testgps,testgpm,testjs,testji,printderiv);
#endif // #ifdef CONTACTAUXPLANE
    // *******************************************************************
    
    // do the two assemblies into the slave nodes
    // if CONTACTONEMORTARLOOP defined, then AssembleM does M AND D matrices !!!
    integrator.AssembleM(Comm(),sele_,mele_,*mseg);
    integrator.AssembleG(Comm(),sele_,*gseg);
  }
  return true;
}

#endif //#ifdef CCADISCRET
