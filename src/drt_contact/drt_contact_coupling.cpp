/*!----------------------------------------------------------------------
\file drt_contact_coupling.cpp
\brief A class for mortar coupling of ONE slave element and ONE master
       element of a contact interface in 2D and 3D.

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
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Intcell::Intcell(int id, int nvertices,
    Epetra_SerialDenseMatrix& coords, const DRT::Element::DiscretizationType& shape) :
id_(id),
nvertices_(nvertices),
coords_(coords),
shape_(shape)
{
   // check nvertices_ and shape_
  if (nvertices_!=3) dserror("ERROR: Integration cell must have 3 vertices");
  if (shape_!=DRT::Element::tri3) dserror("ERROR: Integration cell must be tri3");
  
  // check dimensions of coords_
  if (coords_.M() != 3) dserror("ERROR: Inconsistent coord matrix");
  if (coords_.N() != nvertices_) dserror("ERROR: Inconsistent coord matrix");
  
  // compute area of Intcell
  double t1[3] = {0.0, 0.0, 0.0};
  double t2[3] = {0.0, 0.0, 0.0};
  for (int k=0;k<3;++k)
  {
    t1[k]=Coords()(k,1)-Coords()(k,0);
    t2[k]=Coords()(k,2)-Coords()(k,0);
  }
  
  double t1xt2[3] = {0.0, 0.0, 0.0};
  t1xt2[0] = t1[1]*t2[2]-t1[2]*t2[1];
  t1xt2[1] = t1[2]*t2[0]-t1[0]*t2[2];
  t1xt2[2] = t1[0]*t2[1]-t1[1]*t2[0]; 
  area_ = 0.5*sqrt(t1xt2[0]*t1xt2[0]+t1xt2[1]*t1xt2[1]+t1xt2[2]*t1xt2[2]);
  
  return;
}

/*----------------------------------------------------------------------*
 |  cctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Intcell::Intcell(const Intcell& old) :
id_(old.id_),
nvertices_(old.nvertices_),
area_(old.area_),
coords_(old.coords_),
shape_(old.shape_)
{
  // empty copy constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  Get global coords for given local coords (Intcell)        popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Intcell::LocalToGlobal(const double* xi,
                                               double* globcoord,
                                               int inttype)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");
  
  // collect fundamental data
  int nnodes = NumVertices();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);
  
  // Evaluate shape, get nodal coords and interpolate global coords
  EvaluateShape(xi, val, deriv);
  for (int i=0;i<3;++i) globcoord[i]=0.0;
  
  for (int i=0;i<nnodes;++i)
  {
    if (inttype==0)
    {
      // use shape function values for interpolation
      globcoord[0]+=val[i]*Coords()(0,i);
      globcoord[1]+=val[i]*Coords()(1,i);
      globcoord[2]+=val[i]*Coords()(2,i);
    }
    else if (inttype==1)
    {
      // use shape function derivatives xi for interpolation
      globcoord[0]+=deriv(i,0)*Coords()(0,i);
      globcoord[1]+=deriv(i,0)*Coords()(1,i);
      globcoord[2]+=deriv(i,0)*Coords()(2,i);
    }
    else if (inttype==2)
    {
      // use shape function derivatives eta for interpolation
      globcoord[0]+=deriv(i,1)*Coords()(0,i);
      globcoord[1]+=deriv(i,1)*Coords()(1,i);
      globcoord[2]+=deriv(i,1)*Coords()(2,i);
    }
    else
      dserror("ERROR: Invalid interpolation type requested, only 0,1,2!");
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate shape functions (Intcell)                        popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Intcell::EvaluateShape(const double* xi,
    LINALG::SerialDenseVector& val, LINALG::SerialDenseMatrix& deriv)
{
  if (!xi)
    dserror("ERROR: EvaluateShape called with xi=NULL");
  
  // 3noded triangular element
  if(Shape()==DRT::Element::tri3)
  {
    val[0] = 1-xi[0]-xi[1]; 
    val[1] = xi[0];
    val[2] = xi[1];
    deriv(0,0) = -1.0; deriv(0,1) = -1.0;
    deriv(1,0) =  1.0; deriv(1,1) =  0.0;
    deriv(2,0) =  0.0; deriv(2,1) =  1.0;
  }
  
  // unknown case
  else dserror("ERROR: EvaluateShape (Intcell) called for type != tri3");
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant (Intcell)                   popp 11/08|
 *----------------------------------------------------------------------*/
double CONTACT::Intcell::Jacobian(double* xi)
{
  double jac = 0.0;
  vector<double> gxi(3);
  vector<double> geta(3);
    
  // 2D linear case (2noded line element)
  if (Shape()==DRT::Element::tri3)
    jac = Area()*2;
  
  // unknown case
  else dserror("ERROR: Jacobian called for unknown element type!");
  
  return jac;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Vertex::Vertex(vector<double> coord, Vertex* next, Vertex* prev,
                        bool intersect, bool entryexit, Vertex* neighbor,
                        double alpha) :
coord_(coord),
next_(next),
prev_(prev),
intersect_(intersect),
entryexit_(entryexit),
neighbor_(neighbor),
alpha_(alpha)
{
   // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  cctor (public)                                            popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Vertex::Vertex(const Vertex& old) :
coord_(old.coord_),
next_(old.next_),
prev_(old.prev_),
intersect_(old.intersect_),
entryexit_(old.entryexit_),
neighbor_(old.neighbor_),
alpha_(old.alpha_)
{
   // empty copy constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling::Coupling(DRT::Discretization& idiscret,
                            CONTACT::CElement& sele, CONTACT::CElement& mele,
                            int dim, Epetra_SerialDenseMatrix& csegs) :
idiscret_(idiscret),
sele_(sele),
mele_(mele),
dim_(dim),
contactsegs_(csegs)
{
  // check sanity of dimension
  if (Dim()!=2 && Dim()!=3) dserror("Dim. must be 2D or 3D!");
  
  // the two-dimensional case
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
  
  // the three-dimensional case
  else if (Dim()==3)
  {
    // check for quadratic elements
    if (sele.Shape()!=DRT::Element::tri3 && sele.Shape()!=DRT::Element::quad4)
      dserror("ERROR: 3D mortar coupling not yet impl. for quadratic elements");
    if (mele.Shape()!=DRT::Element::tri3 && mele.Shape()!=DRT::Element::quad4)
      dserror("ERROR: 3D mortar coupling not yet impl. for quadratic elements");
        
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
    Clip() = PolygonClipping3D(SlaveVertices(),MasterVertices(),tol);
    int clipsize = (int)(Clip().size());
    
    // proceed only if clipping polygon is at least triangle
    bool overlap = false;
    if (clipsize>=3) overlap = true;
    if (overlap)
    {
      // do triangulation of clip polygon
      Triangulation3D();
      
      //cout << "\nNo. of integration cells: " << (int)(Cells().size()) << endl;
      //for (int i=0;i<(int)(Cells().size());++i)
      //  cout << "->Cell " << i << ": " << endl << Cells()[i]->Coords() << endl;
      
      // do integration of integration cells
      IntegrateCells3D();
    }
  }
  
  // invalid cases for dim
  else dserror("ERROR: Coupling can only be called for 2D or 3D!");
  
  return;
}

/*----------------------------------------------------------------------*
 |  Project slave / master element pair (public)              popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::Project2D(vector<bool>& hasproj,
                                  vector<double>& xiproj)
{
  // initialize projection status
  hasproj[0] = false;   // slave 0 end node
  hasproj[1] = false;   // slave 1 end node
  hasproj[2] = false;   // master 0 end node
  hasproj[3] = false;   // master 1 end node
  
  // get slave and master element nodes
  DRT::Node** mysnodes = sele_.Nodes();
  if (!mysnodes)
    dserror("ERROR: IntegrateOverlap2D: Null pointer for mysnodes!");
  DRT::Node** mymnodes = mele_.Nodes();
  if (!mymnodes)
      dserror("ERROR: IntegrateOverlap2D: Null pointer for mymnodes!");

  // create a projector instance of problem dimension Dim()
  CONTACT::Projector projector(Dim());

  // project slave nodes onto master element
  for (int i=0;i<sele_.NumNode();++i)
  {
    CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(mysnodes[i]);
    double xi[2] = {0.0, 0.0};
    projector.ProjectNodalNormal(*snode,mele_,xi);

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!
    if ((-1.0-CONTACTPROJTOL<=xi[0]) && (xi[0]<=1.0+CONTACTPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i<2)
      {
        hasproj[i]=true;
        xiproj[i]=xi[0];
      }
      // nevertheless we need the inner node projection status later (weighted gap)
      snode->HasProj()=true;
    }
  }

  // project master nodes onto slave element
  for (int i=0;i<2;++i)
  {
    CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(mymnodes[i]);
    double xi[2] = {0.0, 0.0};
    projector.ProjectElementNormal(*mnode,sele_,xi);

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!!!
    if ((-1.0-CONTACTPROJTOL<=xi[0]) && (xi[0]<=1.0+CONTACTPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i<2)
      {
        hasproj[i+2]=true;
        xiproj[i+2]=xi[0];
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Detect overlap of slave / master pair (public)            popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::DetectOverlap2D(vector<bool>& hasproj,
                                        vector<double>& xiproj)
{
  /**********************************************************************/
  /* OVERLAP CASES                                                      */
  /* Depending on mxi and sxi overlap will be decided!                  */
  /* Even for 3noded CElements only the two end nodes matter in 2D!     */
  /* There are several cases how the 2 elements can overlap. Handle all */
  /* of them, including the ones that they don't overlap at all!        */
  /**********************************************************************/

  // For the non-overlapping cases, the possibility of an identical local
  // node numbering direction for both sides is taken into account!!
  // (this can happen, when elements far from each other are projected,
  // which actually should be impossible due to the CONTACTCRITDIST
  // condition in the potential contact pair search above!
  // But you never know...)

  // For the overlapping cases, it is a prerequisite that the two local
  // node numbering directions are opposite!!
  // (this is the case, when the elements are sufficiently near each other,
  // which is ensured by only processing nodes that fulfill the
  // CONTACTCRITDIST condition above!)

  // CAUTION: The bool output variable in this method is a REAL output
  // variable, determining whether there is an overlap or not!

  // initialize local working variables
  bool overlap = false;
  double sxia = 0.0;
  double sxib = 0.0;
  double mxia = 0.0;
  double mxib = 0.0;

  // local working copies of input variables
  bool s0hasproj = hasproj[0];
  bool s1hasproj = hasproj[1];
  bool m0hasproj = hasproj[2];
  bool m1hasproj = hasproj[3];

  vector<double> sprojxi(2);
  sprojxi[0] = xiproj[0];
  sprojxi[1] = xiproj[1];

  vector<double> mprojxi(2);
  mprojxi[0] = xiproj[2];
  mprojxi[1] = xiproj[3];


  /* CASE 1 (NO OVERLAP):
     no feasible projection found for any of the 4 outer element nodes  */

  if (!s0hasproj && !s1hasproj && !m0hasproj && !m1hasproj)
  {
    //do nothing
  }

  /* CASES 2-5 (NO OVERLAP):
     feasible projection found only for 1 of the 4 outer element nodes
     (this can happen due to the necessary projection tolerance!!!)     */

  else if  (s0hasproj && !s1hasproj && !m0hasproj && !m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=sprojxi[0]) && (sprojxi[0]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap2D: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  else if  (!s0hasproj && s1hasproj && !m0hasproj && !m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=sprojxi[1]) && (sprojxi[1]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap2D: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  else if  (!s0hasproj && !s1hasproj && m0hasproj && !m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=mprojxi[0]) && (mprojxi[0]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap2D: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  else if  (!s0hasproj && !s1hasproj && !m0hasproj && m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=mprojxi[1]) && (mprojxi[1]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap2D: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  /* CASE 6 (OVERLAP):
     feasible projection found for all 4 outer element nodes
     (this can happen due to the necessary projection tolerance!!!)     */

  else if (s0hasproj && s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;

    // internal case 1 for global CASE 6
    // (equivalent to global CASE 7, slave fully projects onto master)
    if ((sprojxi[0]<1.0) && (sprojxi[1]>-1.0))
    {
      sxia = -1.0;
      sxib = 1.0;
      mxia = sprojxi[1];      // local node numbering always anti-clockwise!!!
      mxib = sprojxi[0];
      //cout << "Problem solved with internal case 1!" << endl;
    }

    // internal case 2 for global CASE 6
    // (equivalent to global CASE 8, master fully projects onto slave)
    else if ((mprojxi[0]<1.0) && (mprojxi[1]>-1.0))
    {
      mxia = -1.0;
      mxib = 1.0;
      sxia = mprojxi[1];      // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
      //cout << "Problem solved with internal case 2!" << endl;
    }

    // internal case 3 for global CASE 6
    // (equivalent to global CASE 9, both nodes no. 0 project successfully)
    else if ((sprojxi[0]<1.0+CONTACTPROJLIM) && (mprojxi[0]<1.0+CONTACTPROJLIM))
    {
      sxia = -1.0;
      sxib = mprojxi[0];      // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
      //cout << "Problem solved with internal case 3!" << endl;
    }

    // internal case 4 for global CASE 6
    // (equivalent to global CASE 10, both nodes no. 1 project successfully)
    else if ((sprojxi[1]>-1.0-CONTACTPROJLIM) && (mprojxi[1]>-1.0-CONTACTPROJLIM))
    {
      sxia = mprojxi[1];
      sxib = 1.0;            // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
      //cout << "Problem solved with internal case 4!" << endl;
    }

    // unknown internal case for global CASE 6
    else
    {
      cout << "CONTACT::Interface::IntegrateOverlap2D "<< endl << "has detected '4 projections'-case for Sl./Ma. pair "
           << sele_.Id() << "/" << mele_.Id() << endl;
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap2D: Unknown overlap case found in global case 6!");
    }
  }

  /* CASES 7-8 (OVERLAP):
     feasible projections found for both nodes of one element, this
      means one of the two elements is projecting fully onto the other!  */

  else if (s0hasproj && s1hasproj && !m0hasproj && !m1hasproj)
  {
    overlap = true;
    sxia = -1.0;
    sxib = 1.0;
    mxia = sprojxi[1];      // local node numbering always anti-clockwise!!!
    mxib = sprojxi[0];
  }

  else if (!s0hasproj && !s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;
    mxia = -1.0;
    mxib = 1.0;
    sxia = mprojxi[1];      // local node numbering always anti-clockwise!!!
    sxib = mprojxi[0];
  }

  /* CASES 9-10 (OVERLAP):
     feasible projections found for one node of each element, due to
     node numbering only identical local node ID pairs possible!        */

  else if (s0hasproj && !s1hasproj && m0hasproj && !m1hasproj)
  {
    // do the two elements really have an overlap?
    if ((sprojxi[0]>-1.0+CONTACTPROJLIM) && (mprojxi[0]>-1.0+CONTACTPROJLIM))
    {
      overlap = true;
      sxia = -1.0;
      sxib = mprojxi[0];      // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
    }
  }

  else if (!s0hasproj && s1hasproj && !m0hasproj && m1hasproj)
  {
    // do the two elements really have an overlap?
    if ((sprojxi[1]<1.0-CONTACTPROJLIM) && (mprojxi[1]<1.0-CONTACTPROJLIM))
    {
      overlap = true;
      sxia = mprojxi[1];
      sxib = 1.0;            // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
    }
  }

  /* CASES 11-14 (OVERLAP):
     feasible projections found for 3 out of the total 4 nodes,
     this can either lead to cases 7/8 or 9/10!                         */
  else if (s0hasproj && s1hasproj && m0hasproj && !m1hasproj)
  {
    overlap = true;
    // equivalent to global case 7
    if (mprojxi[0]>1.0)
    {
      sxia = -1.0;
      sxib = 1.0;
      mxia = sprojxi[1];    // local node numbering always anti-clockwise!!!
      mxib = sprojxi[0];
    }
    // equivalent to global case 9
    else
    {
      sxia = -1.0;
      sxib = mprojxi[0];    // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
    }
  }

  else if (s0hasproj && s1hasproj && !m0hasproj && m1hasproj)
  {
    overlap = true;
    // equivalent to global case 7
    if (mprojxi[1]<-1.0)
    {
      sxia = -1.0;
      sxib = 1.0;
      mxia = sprojxi[1];  // local node numbering always anti-clockwise!!!
      mxib = sprojxi[0];
    }
    // equivalent to global case 10
    else
    {
      sxia = mprojxi[1];
      sxib = 1.0;          // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
    }
  }

  else if (s0hasproj && !s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;
    // equivalent to global case 8
    if (sprojxi[0]>1.0)
    {
      mxia = -1.0;
      mxib = 1.0;
      sxia = mprojxi[1];      // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
    }
    // equivalent to global case 9
    else
    {
      sxia = -1.0;
      sxib = mprojxi[0];    // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
    }
  }

  else if (!s0hasproj && s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;
    // equivalent to global case 8
    if (sprojxi[1]<-1.0)
    {
      mxia = -1.0;
      mxib = 1.0;
      sxia = mprojxi[1];  // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
    }
    // equivalent to global case 10
    else
    {
      sxia = mprojxi[1];
      sxib = 1.0;          // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
    }
  }

  /* CASE DEFAULT: unknown overlap case                                  */
  else
  {
    cout << "SElement: " << sele_.NodeIds()[0] << " " << sele_.NodeIds()[1] << endl;
    cout << "MElement: " << mele_.NodeIds()[0] << " " << mele_.NodeIds()[1] << endl;
    cout << "s0: " << s0hasproj << " s1: " << s1hasproj << endl;
    cout << "m0: " << m0hasproj << " m1: " << m1hasproj << endl;
    dserror("ERROR: IntegrateOverlap2D: Unknown overlap case found!");
  }

  // check for 1:1 node projections and for infeasible limits
  if ((sxia<-1.0) || (sxib>1.0) || (mxia<-1.0) || (mxib>1.0))
  {
    if (abs(sxia+1.0)<CONTACTPROJLIM) sxia=-1.0;
    if (abs(sxib-1.0)<CONTACTPROJLIM) sxib= 1.0;
    if (abs(mxia+1.0)<CONTACTPROJLIM) mxia=-1.0;
    if (abs(mxib-1.0)<CONTACTPROJLIM) mxib= 1.0;

    if ((sxia<-1.0) || (sxib>1.0) || (mxia<-1.0) || (mxib>1.0))
    {
      cout << "Slave: " << sxia << " " << sxib << endl;
      cout << "Master: " << mxia << " " << mxib << endl;
      dserror("ERROR: IntegrateOverlap2D: Determined infeasible limits!");
    }
  }

  // update integration limits in xiproj
  xiproj[0]=sxia;
  xiproj[1]=sxib;
  xiproj[2]=mxia;
  xiproj[3]=mxib;

  // prepare gmsh visualization
#ifdef DEBUG
  if (overlap)
  {
    double sxialoc[2] = {sxia, 0.0};
    double sxibloc[2] = {sxib, 0.0};
    double mxialoc[2] = {mxia, 0.0};
    double mxibloc[2] = {mxib, 0.0};

    double sxiaglob[3] = {0.0, 0.0, 0.0};
    double sxibglob[3] = {0.0, 0.0, 0.0};
    double mxiaglob[3] = {0.0, 0.0, 0.0};
    double mxibglob[3] = {0.0, 0.0, 0.0};

    sele_.LocalToGlobal(sxialoc,sxiaglob,0);
    sele_.LocalToGlobal(sxibloc,sxibglob,0);
    mele_.LocalToGlobal(mxialoc,mxiaglob,0);
    mele_.LocalToGlobal(mxibloc,mxibglob,0);

    Epetra_SerialDenseMatrix& segs = CSegs();
    segs.Reshape(segs.M()+1,12);
    segs(segs.M()-1,0)  = sxiaglob[0];
    segs(segs.M()-1,1)  = sxiaglob[1];
    segs(segs.M()-1,2)  = sxiaglob[2];
    segs(segs.M()-1,3)  = mxibglob[0];
    segs(segs.M()-1,4)  = mxibglob[1];
    segs(segs.M()-1,5)  = mxibglob[2];
    segs(segs.M()-1,6)  = mxiaglob[0];
    segs(segs.M()-1,7)  = mxiaglob[1];
    segs(segs.M()-1,8)  = mxiaglob[2];
    segs(segs.M()-1,9)  = sxibglob[0];
    segs(segs.M()-1,10) = sxibglob[1];
    segs(segs.M()-1,11) = sxibglob[2];
  }
#endif // #ifdef DEBUG

  return overlap;
}

/*----------------------------------------------------------------------*
 |  Integrate slave / master overlap (public)                 popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::IntegrateOverlap2D(vector<double>& xiproj)
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj entries integrate the Mortar   */
  /* matrix M and the weighted gap function g~ on the overlap of the    */
  /* current slave / master CElement pair                               */
  /**********************************************************************/

  //local working copies of input variables
  double sxia = xiproj[0];
  double sxib = xiproj[1];
  double mxia = xiproj[2];
  double mxib = xiproj[3];

  // create an integrator instance with correct NumGP and Dim
  CONTACT::Integrator integrator(sele_.Shape());

  // do the two integrations
  RCP<Epetra_SerialDenseMatrix> mseg = integrator.IntegrateM(sele_,sxia,sxib,mele_,mxia,mxib);
  RCP<Epetra_SerialDenseVector> gseg = integrator.IntegrateG(sele_,sxia,sxib,mele_,mxia,mxib);

  // compute directional derivative of M and store into nodes
  // if CONTACTONEMORTARLOOP defined, then DerivM does linearization of M AND D matrices !!!
  integrator.DerivM(sele_,sxia,sxib,mele_,mxia,mxib);
    
  // do the two assemblies into the slave nodes
  // if CONTACTONEMORTARLOOP defined, then AssembleM does M AND D matrices !!!
  integrator.AssembleM(Comm(),sele_,mele_,*mseg);
  integrator.AssembleG(Comm(),sele_,*gseg);

  
  /*----------------------------------------------------------------------
  // check for the modification of the M matrix for curved interfaces
  // (based on the paper by M. Puso / B. Wohlmuth, IJNME, 2005)
  // (note: we assume that the modification is not useful for mortar CONTACT,
  // but only for mortar MESH TYING as published!)
  //----------------------------------------------------------------------
  bool modification = false;

  // conditions for modification
  // (1) linear shape functions for the slave elements
  // (2) 2D problems (3D unknown / unpublished ??)
  // (3) curved interface, n1!=n2
  // (4) use of dual shape functions for LM (always true in our case !!)
  if (sele_.Shape()==DRT::Element::line2)
  {
    if (integrator.Dim()==2)
    {
      CNode* snode0 = static_cast<CNode*>(sele_.Nodes()[0]);
      CNode* snode1 = static_cast<CNode*>(sele_.Nodes()[1]);

      const double* n0 = snode0->n();
      const double* n1 = snode1->n();
      double delta = (n0[0]-n1[0])*(n0[0]-n1[0]) + (n0[1]-n1[1])*(n0[1]-n1[1]);

      if (delta>1.0e-8)
        modification = true;
    }
  }

  // integrate and assemble the modification, if necessary

  if (modification)
  {
    RCP<Epetra_SerialDenseMatrix> mmodseg = integrator.IntegrateMmod(sele_,sxia,sxib,mele_,mxia,mxib);
    integrator.AssembleMmod(Comm(),sele_,mele_,*mmodseg);
  }
  //--------------------------------------------------------------------*/
  
  return true;
}


/*----------------------------------------------------------------------*
 |  Build auxiliary plane from slave element (public)         popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::AuxiliaryPlane3D()
{
  // we first need the element center:
  // for quad4, quad8, quad9 elements: xi = eta = 0.0
  // for tri3, tri6 elements: xi = eta = 1/3
  double loccenter[2];
    
  DRT::Element::DiscretizationType dt = SlaveElement().Shape();
  if (dt==CElement::tri3 || dt==CElement::tri6)
  {
    loccenter[0] = 1/3;
    loccenter[1] = 1/3;
  }
  else if (dt==CElement::quad4 || dt==CElement::quad8 || dt==CElement::quad9)
  {
    loccenter[0] = 0.0;
    loccenter[1] = 0.0;
  }
  else dserror("ERROR: AuxiliaryPlane3D called for unknown element type");
  
  // compute element center via shape fct. interpolation
  SlaveElement().LocalToGlobal(loccenter,Auxc(),0);
  
  // we then compute the unit normal vector at the element center
  SlaveElement().ComputeUnitNormalAtXi(loccenter,Auxn());
  
  //cout << "Slave Element: " << SlaveElement().Id() << endl;
  //cout << "->Center: " << Auxc()[0] << " " << Auxc()[1] << " " << Auxc()[2] << endl;
  //cout << "->Normal: " << Auxn()[0] << " " << Auxn()[1] << " " << Auxn()[2] << endl;
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Project slave element onto auxiliary plane (public)       popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::ProjectSlave3D()
{
  // project slave nodes onto auxiliary plane
  int nnodes = SlaveElement().NumNode();
  DRT::Node** mynodes = SlaveElement().Nodes();
  if (!mynodes) dserror("ERROR: ProjectSlave3D: Null pointer!");
  
  // initialize storage for slave vertices
  vector<vector<double> > vertices(nnodes,vector<double>(3));
  
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: ProjectSlave3D: Null pointer!");
    
    // first build difference of point and element center
    // and then dot product with unit normal at center
    double dist = (mycnode->xspatial()[0]-Auxc()[0])*Auxn()[0]
                + (mycnode->xspatial()[1]-Auxc()[1])*Auxn()[1]
                + (mycnode->xspatial()[2]-Auxc()[2])*Auxn()[2];
    
    // compute projection
    for (int k=0;k<3;++k) vertices[i][k] = mycnode->xspatial()[k] - dist * Auxn()[k];
 
    //cout << "->RealNode(S) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " << mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << endl; 
    //cout << dist << endl;
    //cout << "->ProjNode(S) " << mycnode->Id() << ": " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << endl; 
  }
  
  // store slave vertices into svertices_
  svertices_ = vertices;
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Project master element onto auxiliary plane (public)      popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::ProjectMaster3D()
{
  // project master nodes onto auxiliary plane
  int nnodes = MasterElement().NumNode();
  DRT::Node** mynodes = MasterElement().Nodes();
  if (!mynodes) dserror("ERROR: ProjectMaster3D: Null pointer!");
  
  // initialize storage for master vertices
  vector<vector<double> > vertices(nnodes,vector<double>(3));
    
  for (int i=0;i<nnodes;++i)
  {
    CNode* mycnode = static_cast<CNode*> (mynodes[i]);
    if (!mycnode) dserror("ERROR: ProjectMaster3D: Null pointer!");
    
    // first build difference of point and element center
    // and then dot product with unit normal at center
    double dist = (mycnode->xspatial()[0]-Auxc()[0])*Auxn()[0]
                + (mycnode->xspatial()[1]-Auxc()[1])*Auxn()[1]
                + (mycnode->xspatial()[2]-Auxc()[2])*Auxn()[2];
    
    // compute projection
    for (int k=0;k<3;++k) vertices[i][k] = mycnode->xspatial()[k] - dist * Auxn()[k];
    
    //cout << "->RealNode(M) " << mycnode->Id() << ": " << mycnode->xspatial()[0] << " " << mycnode->xspatial()[1] << " " << mycnode->xspatial()[2] << endl; 
    //cout << dist << endl;
    //cout << "->ProjNode(M) " << mycnode->Id() << ": " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << endl; 
  }
  
  // store master vertices into mvertices_
  mvertices_ = vertices;
    
  return true;
}

/*----------------------------------------------------------------------*
 |  Clipping of two polygons                                  popp 11/08|
 *----------------------------------------------------------------------*/
vector<vector<double> > CONTACT::Coupling::PolygonClipping3D(
    vector<vector<double> > poly1, vector<vector<double> > poly2, double& tol)
{
  // print to screen
  //cout << "\n\n*****************************************************";
  //cout << "\n*          P O L Y G O N   C L I P P I N G          *";
  //cout << "\n*****************************************************\n\n";
  
  //**********************************************************************
  // STEP1: Input check
  // - input polygons must consist of min. 3 vertices each
  // - rotation of polygon 1 must be c-clockwise w.r.t. Auxn()
  // - rotation of polygon 2 is changed to c-clockwise w.r.t. Auxn()
  // - both input polygons must be convex
  //**********************************************************************
  
  // check input variables
  if ((int)poly1.size()<3 || (int)poly2.size()<3)
    dserror("ERROR: Input Polygons must consist of min. 3 vertices each");
  
  // check for rotation of polygon1 (slave) and polgon 2 (master)
  // note that we implicitly already rely on convexity here!
  // first get geometric centers of polygon1 and polygon2
  double center1[3] = {0.0, 0.0, 0.0};
  double center2[3] = {0.0, 0.0, 0.0};
  
  for (int i=0;i<(int)poly1.size();++i)
    for (int k=0;k<3;++k)
      center1[k] += poly1[i][k]/((int)poly1.size());
  
  for (int i=0;i<(int)poly2.size();++i)
    for (int k=0;k<3;++k)
      center2[k] += poly2[i][k]/((int)poly2.size());
  
  //cout << "Center 1: " << center1[0] << " " << center1[1] << " " << center1[2] << endl;
  //cout << "Center 2: " << center2[0] << " " << center2[1] << " " << center2[2] << endl;
  
  // then we compute the counter-clockwise plane normal
  double diff1[3] = {0.0, 0.0, 0.0};
  double edge1[3] = {0.0, 0.0, 0.0};
  double diff2[3] = {0.0, 0.0, 0.0};
  double edge2[3] = {0.0, 0.0, 0.0};
  
  for (int k=0;k<3;++k)
  {
    diff1[k] = poly1[0][k]-center1[k];
    edge1[k] = poly1[1][k]-poly1[0][k];
    diff2[k] = poly2[0][k]-center2[k];
    edge2[k] = poly2[1][k]-poly2[0][k];
  }
  
  double cross1[3] = {0.0, 0.0, 0.0};
  double cross2[3] = {0.0, 0.0, 0.0};
  
  cross1[0] = diff1[1]*edge1[2]-diff1[2]*edge1[1];
  cross1[1] = diff1[2]*edge1[0]-diff1[0]*edge1[2];
  cross1[2] = diff1[0]*edge1[1]-diff1[1]*edge1[0];
  
  cross2[0] = diff2[1]*edge2[2]-diff2[2]*edge2[1];
  cross2[1] = diff2[2]*edge2[0]-diff2[0]*edge2[2];
  cross2[2] = diff2[0]*edge2[1]-diff2[1]*edge2[0];
  
  // check against auxiliary plane normal
  double check1 = cross1[0]*Auxn()[0]+cross1[1]*Auxn()[1]+cross1[2]*Auxn()[2];
  double check2 = cross2[0]*Auxn()[0]+cross2[1]*Auxn()[1]+cross2[2]*Auxn()[2];
  
  // check polygon 1 and throw dserror if not c-clockwise
  if (check1<=0) dserror("ERROR: Polygon 1 (slave) not ordered counter-clockwise!");
  
  // check polygon 2 and reorder in c-clockwise direction
  if (check2<0)
  {
    //cout << "Polygon 2 (master) not ordered counter-clockwise -> reordered!" << endl;
    vector<vector<double> > newpoly2((int)poly2.size(),vector<double>(3));
    for (int i=0;i<(int)poly2.size();++i)
      newpoly2[(int)poly2.size()-1-i] = poly2[i];
    poly2 = newpoly2;
  }
   
  // check if the two input polygons are convex
  // a polygon is convex if the scalar product of an edge normal and the
  // next edge direction is negative for all edges
  for (int i=0;i<(int)poly1.size();++i)
  {
    // we need the edge vector first
    double edge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i!=(int)poly1.size()-1) edge[k] = poly1[i+1][k] - poly1[i][k];
      else edge[k] = poly1[0][k] - poly1[i][k];
    }
    
    // edge normal is result of cross product
    double n[3] = {0.0, 0.0, 0.0};
    n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
    n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
    n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
    
    // we need the next edge vector now
    double nextedge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i<(int)poly1.size()-2) nextedge[k] = poly1[i+2][k] - poly1[i+1][k];
      else if (i==(int)poly1.size()-2) nextedge[k] = poly1[0][k] - poly1[i+1][k];
      else nextedge[k] = poly1[1][k] - poly1[0][k];
    }
    
    // check scalar product
    double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
    if (check>0) dserror("ERROR: Input polygon 1 not convex");
  }
  
  for (int i=0;i<(int)poly2.size();++i)
  {
    // we need the edge vector first
    double edge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i!=(int)poly2.size()-1) edge[k] = poly2[i+1][k] - poly2[i][k];
      else edge[k] = poly2[0][k] - poly2[i][k];
    }
    
    // edge normal is result of cross product
    double n[3] = {0.0, 0.0, 0.0};
    n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
    n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
    n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
    
    // we need the next edge vector now
    double nextedge[3] = {0.0, 0.0, 0.0};
    for (int k=0;k<3;++k)
    {
      if (i<(int)poly2.size()-2) nextedge[k] = poly2[i+2][k] - poly2[i+1][k];
      else if (i==(int)poly2.size()-2) nextedge[k] = poly2[0][k] - poly2[i+1][k];
      else nextedge[k] = poly2[1][k] - poly2[0][k];
    }
    
    // check scalar product
    double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
    if (check>0) dserror("ERROR: Input polygon 2 not convex");
  }
  
  // print final input polygons to screen
  //cout << "\nInput Poylgon 1:";
  //for (int i=0;i<(int)poly1.size();++i)
  //  cout << "\nVertex " << i << ":\t" << scientific << poly1[i][0] << "\t" << poly1[i][1] << "\t" << poly1[i][2];
 
  //cout << "\nInput Poylgon 2:";
  //for (int i=0;i<(int)poly2.size();++i)
  //  cout << "\nVertex " << i << ":\t" << scientific << poly2[i][0] << "\t" << poly2[i][1] << "\t" << poly2[i][2];
  
  //cout << endl << endl;
  
  //**********************************************************************
  // STEP2: Create Vertex data structures
  // - convert the input vectors into a vector of Vertex objects each
  // - assign Next() and Prev() pointers to initialize linked structure
  //**********************************************************************
  vector<Vertex> poly1list;
  vector<Vertex> poly2list;
  
  // create Vertex types from input data
  for (int i=0;i<(int)poly1.size();++i)
    poly1list.push_back(Vertex(poly1[i],NULL,NULL,false,false,NULL,-1.0));

  for (int i=0;i<(int)poly2.size();++i)
    poly2list.push_back(Vertex(poly2[i],NULL,NULL,false,false,NULL,-1.0));

  // set previous and next Vertex pointer for all elements in lists
  for (int i=0;i<(int)poly1list.size();++i)
  {
    // standard case
    if (i!=0 && i!=(int)poly1list.size()-1)
    {
      poly1list[i].AssignNext(&poly1list[i+1]);
      poly1list[i].AssignPrev(&poly1list[i-1]);
    }
    // first element in list
    else if (i==0)
    {
      poly1list[i].AssignNext(&poly1list[i+1]);
      poly1list[i].AssignPrev(&poly1list[(int)poly1list.size()-1]);
    }
    // last element in list
    else
    {
      poly1list[i].AssignNext(&poly1list[0]);
      poly1list[i].AssignPrev(&poly1list[i-1]);
    }
  }
  for (int i=0;i<(int)poly2list.size();++i)
  {
    // standard case
    if (i!=0 && i!=(int)poly2list.size()-1)
    {
      poly2list[i].AssignNext(&poly2list[i+1]);
      poly2list[i].AssignPrev(&poly2list[i-1]);
    }
    // first element in list
    else if (i==0)
    {
      poly2list[i].AssignNext(&poly2list[i+1]);
      poly2list[i].AssignPrev(&poly2list[(int)poly2list.size()-1]);
    }
    // last element in list
    else
    {
      poly2list[i].AssignNext(&poly2list[0]);
      poly2list[i].AssignPrev(&poly2list[i-1]);
    }
  }
  
  //**********************************************************************
  // STEP3: Avoid degenerate cases
  // - if a point of poly1 is close (<tol) to a edge of poly2 or vice
  //   versa we move this point away from the edge by tol
  //**********************************************************************
  for (int i=0;i<(int)poly1list.size();++i)
  {
    for (int j=0;j<(int)poly2list.size();++j)
    {
      // we need diff vector and edge2 first
      double diff1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diff1[k] = poly1list[i].Coord()[k] - poly2list[j].Coord()[k];
        edge2[k] = (poly2list[j].Next())->Coord()[k] - poly2list[j].Coord()[k];
      }
      
      // check if point of poly1 lies within [0,1] for edge2
      double checkalpha = diff1[0]*edge2[0]+diff1[1]*edge2[1]+diff1[2]*edge2[2];
      checkalpha /= (edge2[0]*edge2[0]+edge2[1]*edge2[1]+edge2[2]*edge2[2]);
      
      // proceed only if inside [0,1] with tolerance tol
      if (checkalpha<-tol || checkalpha>1+tol) continue;
      
      // compute distance from point on poly1 to edge2
      double n2[3] = {0.0, 0.0, 0.0};
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];
      double ln = sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
      for (int k=0;k<3;++k) n2[k] /= ln;
      
      double dist = diff1[0]*n2[0]+diff1[1]*n2[1]+diff1[2]*n2[2];
      
      // move point away if very close to edge 2
      if (dist > -tol && dist < 0)
      {
        //cout << "Vertex " << i << " on poly1 is very close to edge " << j << " of poly2 -> moved inside!" << endl;
        poly1list[i].Coord()[0] -= tol*n2[0];
        poly1list[i].Coord()[1] -= tol*n2[1];
        poly1list[i].Coord()[2] -= tol*n2[2];                                   
      }
      else if (dist < tol && dist >= 0)
      {
        //cout << "Vertex " << i << " on poly1 is very close to edge " << j << " of poly2 -> moved outside!" << endl;
        poly1list[i].Coord()[0] += tol*n2[0];
        poly1list[i].Coord()[1] += tol*n2[1];
        poly1list[i].Coord()[2] += tol*n2[2];     
      }
      else
      {
        // do nothing, point is not very close
      }    
    }
  }
  
  for (int i=0;i<(int)poly2list.size();++i)
  {
    for (int j=0;j<(int)poly1list.size();++j)
    {
      // we need diff vector and edge1 first
      double diff2[3] = {0.0, 0.0, 0.0};
      double edge1[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diff2[k] = poly2list[i].Coord()[k] - poly1list[j].Coord()[k];
        edge1[k] = (poly1list[j].Next())->Coord()[k] - poly1list[j].Coord()[k];
      }
      
      // check if point of poly2 lies within [0,1] for edge1
      double checkalpha = diff2[0]*edge1[0]+diff2[1]*edge1[1]+diff2[2]*edge1[2];
      checkalpha /= (edge1[0]*edge1[0]+edge1[1]*edge1[1]+edge1[2]*edge1[2]);
      
      // proceed only if inside [0,1] with tolerance tol
      if (checkalpha<-tol || checkalpha>1+tol) continue;
      
      // compute distance from point on poly2 to edge1
      double n1[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      double ln = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
      for (int k=0;k<3;++k) n1[k] /= ln;
      
      double dist = diff2[0]*n1[0]+diff2[1]*n1[1]+diff2[2]*n1[2];
      
      // move point away if very close to edge 2
      if (dist > -tol && dist < 0)
      {
        //cout << "Vertex " << i << " on poly2 is very close to edge " << j << " of poly1 -> moved inside!" << endl;
        poly2list[i].Coord()[0] -= tol*n1[0];
        poly2list[i].Coord()[1] -= tol*n1[1];
        poly2list[i].Coord()[2] -= tol*n1[2];                                   
      }
      else if (dist < tol && dist >= 0)
      {
        //cout << "Vertex " << i << " on poly2 is very close to edge " << j << " of poly1 -> moved outside!" << endl;
        poly2list[i].Coord()[0] += tol*n1[0];
        poly2list[i].Coord()[1] += tol*n1[1];
        poly2list[i].Coord()[2] += tol*n1[2];     
      }
      else
      {
        // do nothing, point is not very close
      }    
    }
  }
  
  
  //**********************************************************************
  // STEP4: Perform line intersection of all edge pairs
  // - this yields two new vectors of intersection vertices
  // - by default the respective edge end vertices are assumed to be
  //   the next/prev vertices and connectivity is set up accordingly
  //**********************************************************************
  vector<Vertex> intersec1;
  vector<Vertex> intersec2;
  
  for (int i=0;i<(int)poly1list.size();++i)
  {
    for (int j=0;j<(int)poly2list.size();++j)
    {
      // we need two diff vectors and edges first
      double diffp1[3] = {0.0, 0.0, 0.0};
      double diffp2[3] = {0.0, 0.0, 0.0};
      double diffq1[3] = {0.0, 0.0, 0.0};
      double diffq2[3] = {0.0, 0.0, 0.0};
      double edge1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        diffp1[k] = poly1list[i].Coord()[k] - poly2list[j].Coord()[k];
        diffp2[k] = (poly1list[i].Next())->Coord()[k] - poly2list[j].Coord()[k];
        diffq1[k] = poly2list[j].Coord()[k] - poly1list[i].Coord()[k];
        diffq2[k] = (poly2list[j].Next())->Coord()[k] - poly1list[i].Coord()[k];
        edge1[k] = (poly1list[i].Next())->Coord()[k] - poly1list[i].Coord()[k];
        edge2[k] = (poly2list[j].Next())->Coord()[k] - poly2list[j].Coord()[k];
      }
      
      // outward edge normals of polygon 1 and 2 edges
      double n1[3] = {0.0, 0.0, 0.0};
      double n2[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];
      
      // check for parallelity of edges
      double parallel = edge1[0]*n2[0]+edge1[1]*n2[1]+edge1[2]*n2[2];
      if(abs(parallel)<1.0e-12)
      {
        //cout << "WARNING: Detected two parallel edges! (" << i << "," << j << ")" << endl;
        continue;
      }
      
      //cout << "Searching intersection (" << i << "," << j << ")" << endl;
      
      // check for intersection of non-parallel edges
      double wec_p1 = 0.0;
      double wec_p2 = 0.0;
      for (int k=0;k<3;++k)
      {
        wec_p1 += (poly1list[i].Coord()[k] - poly2list[j].Coord()[k]) * n2[k];
        wec_p2 += ((poly1list[i].Next())->Coord()[k] - poly2list[j].Coord()[k]) * n2[k];
      }
     
      if (wec_p1*wec_p2<=0)
      {
        double wec_q1 = 0.0;
        double wec_q2 = 0.0;
        for (int k=0;k<3;++k)
        {
          wec_q1 += (poly2list[j].Coord()[k] - poly1list[i].Coord()[k]) * n1[k];
          wec_q2 += ((poly2list[j].Next())->Coord()[k] - poly1list[i].Coord()[k]) * n1[k];
        }
        
        if (wec_q1*wec_q2<=0)
        {
          double alphap = wec_p1/(wec_p1-wec_p2);
          double alphaq = wec_q1/(wec_q1-wec_q2);
          vector<double> ip(3);
          vector<double> iq(3);
          for (int k=0;k<3;++k)
          {
            ip[k] = (1-alphap) * poly1list[i].Coord()[k] + alphap * (poly1list[i].Next())->Coord()[k];
            iq[k] = (1-alphaq) * poly2list[j].Coord()[k] + alphaq * (poly2list[j].Next())->Coord()[k];
            if (abs(ip[k])<1.0e-12) ip[k] = 0.0;
            if (abs(iq[k])<1.0e-12) iq[k] = 0.0;
          }
          
          //cout << "Found intersection! (" << i << "," << j << ") " << alphap << " " << alphaq << endl;
          //cout << "On Polygon 1: " << ip[0] << " " << ip[1] << " " << ip[2] << endl;
          //cout << "On Polygon 2: " << iq[0] << " " << iq[1] << " " << iq[2] << endl;
          
          intersec1.push_back(Vertex(ip,poly1list[i].Next(),&poly1list[i],true,false,NULL,alphap));
          intersec2.push_back(Vertex(iq,poly2list[j].Next(),&poly2list[j],true,false,NULL,alphaq));
        }
      }
    }
  }
  
  // do clipping
  vector<vector<double> > respoly;
    
  //**********************************************************************
  // STEP5: Find result polygon for no intersection case
  // - if there are no intersections polygon 1 could lie within polygon 2,
  //   polygon 2 could lie within polygon 1 or they are fully adjacent
  // - by default the respective edge end vertices are assumed to be
  //   the next/prev vertices and connectivity is set up accordingly
  //**********************************************************************
  if ((int)intersec1.size()==0 && (int)intersec2.size()==0)
  {
    // if they are nested the clip polygon is the inner input polygon
    bool poly1inner = true;
    bool poly2inner = true;
    
    // (A) check if poly1list[0] inside poly2
    for (int i=0;i<(int)poly2list.size();++i)
    {
      double edge[3] = {0.0, 0.0, 0.0};
      double diff[3] = {0.0, 0.0, 0.0};
      
      for (int k=0;k<3;++k)
      {
        edge[k] = (poly2list[i].Next())->Coord()[k] - poly2list[i].Coord()[k];
        diff[k] = poly1list[0].Coord()[k] - poly2list[i].Coord()[k];
      }
      
      double n[3] = {0.0, 0.0, 0.0};
      n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
      n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
      n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
      
      double check = diff[0]*n[0] + diff[1]*n[1] + diff[2]*n[2];
      
      // if check>0 then poly1list[0] NOT inside poly2
      if (check>0)
      {
        poly1inner = false;
        break;
      }
    }
    
    if (poly1inner==true)
    {
      //cout << "Polygon S lies fully inside polygon M!" << endl; 
      respoly = poly1;
    }
    
    else
    {
      // (A) check if poly2list[0] inside poly1
      for (int i=0;i<(int)poly1list.size();++i)
      {
        double edge[3] = {0.0, 0.0, 0.0};
        double diff[3] = {0.0, 0.0, 0.0};
        
        for (int k=0;k<3;++k)
        {
          edge[k] = (poly1list[i].Next())->Coord()[k] - poly1list[i].Coord()[k];
          diff[k] = poly2list[0].Coord()[k] - poly1list[i].Coord()[k];
        }
        
        double n[3] = {0.0, 0.0, 0.0};
        n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
        n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
        n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
        
        double check = diff[0]*n[0] + diff[1]*n[1] + diff[2]*n[2];
        
        // if check>0 then poly2list[0] NOT inside poly1
        if (check>0)
        {
          poly2inner = false;
          break;
        }
      }
      
      if (poly2inner==true)
      {
        //cout << "Polygon M lies fully inside polygon S!" << endl; 
        respoly = poly2;
      }
      
      // fully adjacent case
      else
      {
        //cout << "Polygons S and M are fully adjacent!" << endl; 
        vector<vector<double> > empty(0,vector<double>(3));
        respoly = empty;
      }
    }
  }
  
  // invalid case
  else if ((int)intersec1.size()*(int)intersec2.size()==0)
  {
    dserror("ERROR: Found intersection points only on one input polygon...?");
  }
  
  //**********************************************************************
  // STEP6: Find result polygon for intersection case
  // - assign neighbor connectivity for intersection points
  // - check for edges where 2 intersection points have been found
  // - establish new connectivity (next/prev) accordingly: first for
  //   the intersection points, then for the adjacent edge nodes
  // - perform entry exit classification of intersection points
  // - build result polygon (path finding through linked data structures)
  // - check for sanity of result polygon (last==first)
  // - reorder result polygon in c-clockwise direction if necessary
  // - check if result polygon is convex
  //**********************************************************************
  else
  {
    // assign neighbor intersection nodes on the respective other polygon
    for (int i=0;i<(int)intersec1.size();++i)
      intersec1[i].AssignNeighbor(&intersec2[i]);
    for (int i=0;i<(int)intersec2.size();++i)
      intersec2[i].AssignNeighbor(&intersec1[i]);
        
    // check all edges for double intersections
    for (int i=0;i<(int)poly1list.size();++i)
    {
      vector<Vertex*> dis;
      for (int z=0;z<(int)intersec1.size();++z)
      {
       if (intersec1[z].Next()==poly1list[i].Next() && intersec1[z].Prev()==&poly1list[i])
       {
         dis.push_back(&intersec1[z]);
       }
      }
      
      if ((int)dis.size()<2) continue;
      if ((int)dis.size()>2) dserror("ERROR: More than 2 intersections on 1 edge impossible!");
      
      double alpha1 = dis[0]->Alpha();
      double alpha2 = dis[1]->Alpha();
      
      if (alpha1<alpha2)
      {
        // ordering is polylist1[i] -> dis[0] -> dis[1] -> polylist1[i].Next()
        dis[0]->AssignNext(dis[1]);
        dis[1]->AssignPrev(dis[0]);
      }
      else if (alpha1==alpha2)
      {
        dserror("ERROR: Two identical intersection points on 1 edge!");
      }
      else
      {
        // ordering is polylist1[i] -> dis[1] -> dis[0] -> polylist1[i].Next()
        dis[1]->AssignNext(dis[0]);
        dis[0]->AssignPrev(dis[1]);
      }
    }
    
    for (int i=0;i<(int)poly2list.size();++i)
    {
      vector<Vertex*> dis;
      for (int z=0;z<(int)intersec2.size();++z)
      {
       if (intersec2[z].Next()==poly2list[i].Next() && intersec2[z].Prev()==&poly2list[i])
       {
         dis.push_back(&intersec2[z]);
       }
      }
      
      if ((int)dis.size()<2) continue;
      if ((int)dis.size()>2) dserror("ERROR: More than 2 intersections on 1 edge impossible!");
      
      double alpha1 = dis[0]->Alpha();
      double alpha2 = dis[1]->Alpha();
      
      if (alpha1<alpha2)
      {
        // ordering is polylist2[i] -> dis[0] -> dis[1] -> polylist2[i].Next()
        dis[0]->AssignNext(dis[1]);
        dis[1]->AssignPrev(dis[0]);
      }
      else if (alpha1==alpha2)
      {
        dserror("ERROR: Two identical intersection points on 1 edge!");
      }
      else
      {
        // ordering is polylist2[i] -> dis[1] -> dis[0] -> polylist2[i].Next()
        dis[1]->AssignNext(dis[0]);
        dis[0]->AssignPrev(dis[1]);
      }
    }
    
    // assign new next / previous nodes for vertices near intersections
    for (int i=0;i<(int)intersec1.size();++i)
    {
      (intersec1[i].Prev())->AssignNext(&intersec1[i]);
      (intersec1[i].Next())->AssignPrev(&intersec1[i]);
    }
    for (int i=0;i<(int)intersec2.size();++i)
    {
      (intersec2[i].Prev())->AssignNext(&intersec2[i]);
      (intersec2[i].Next())->AssignPrev(&intersec2[i]);
    }  

    // perform entry / exit classification of intersections
    // we move along both polygons and determine whether each intersection
    // point is an entry or exit point with respect to the other polygon.
    // this status is then stored into the vertex data structure.
    
    for (int i=0;i<(int)intersec1.size();++i)
    {
      // check if previous vertex is inside for first intersection
      double edge1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        edge1[k] = (intersec1[i].Next())->Coord()[k] - (intersec1[i].Prev())->Coord()[k];
        edge2[k] = ((intersec1[i].Neighbor())->Next())->Coord()[k] - ((intersec1[i].Neighbor())->Prev())->Coord()[k];
      }
      double n1[3] = {0.0, 0.0, 0.0};
      double n2[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];
      
      double check = edge1[0]*n2[0] + edge1[1]*n2[1] + edge1[2]*n2[2];
      if (check<0) intersec1[i].EntryExit()=true;
    }
    
    for (int i=0;i<(int)intersec1.size();++i)
    {
      // check if previous vertex is inside for first intersection
      double edge1[3] = {0.0, 0.0, 0.0};
      double edge2[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        edge1[k] = (intersec2[i].Next())->Coord()[k] - (intersec2[i].Prev())->Coord()[k];
        edge2[k] = ((intersec2[i].Neighbor())->Next())->Coord()[k] - ((intersec2[i].Neighbor())->Prev())->Coord()[k];
      }
      double n1[3] = {0.0, 0.0, 0.0};
      double n2[3] = {0.0, 0.0, 0.0};
      n1[0] = edge1[1]*Auxn()[2]-edge1[2]*Auxn()[1];
      n1[1] = edge1[2]*Auxn()[0]-edge1[0]*Auxn()[2];
      n1[2] = edge1[0]*Auxn()[1]-edge1[1]*Auxn()[0];
      n2[0] = edge2[1]*Auxn()[2]-edge2[2]*Auxn()[1];
      n2[1] = edge2[2]*Auxn()[0]-edge2[0]*Auxn()[2];
      n2[2] = edge2[0]*Auxn()[1]-edge2[1]*Auxn()[0];
      
      double check = edge1[0]*n2[0] + edge1[1]*n2[1] + edge1[2]*n2[2];
      if (check<0) intersec2[i].EntryExit()=true;
    }
    
    // print intersection points and their status
    //cout << endl;
    for (int i=0;i<(int)intersec1.size();++i)
    {
      //cout << "Intersec1: " << i << " " << intersec1[i].Coord()[0] << " " << intersec1[i].Coord()[1] << " " << intersec1[i].Coord()[2];
      //cout << " EntryExit: " << intersec1[i].EntryExit() << endl;
    }
    
    //cout << endl;
    for (int i=0;i<(int)intersec2.size();++i)
    {
     // cout << "Intersec2: " << i << " " <<  intersec2[i].Coord()[0] << " " << intersec2[i].Coord()[1] << " " << intersec2[i].Coord()[2];
      //cout << " EntryExit: " << intersec2[i].EntryExit() << endl;
    }
    
    // create clipped polygon by filtering
    // We simply have to find our way through the linked data structures of
    // poly1, poly2 and intersection vertices. For this we start at an
    // intersection point and move on according to the entry / exit status.
    // When we reach the next intersection point we jump to the other polygon
    // according to the neighboring vertex pointer.
    // The result will be an ordered list of vertices of the clipped polygon!
    Vertex* current = &intersec1[0];
    
    // push_back start Vertex coords into result polygon
    //cout << "\nStart loop on Slave at " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
    vector<double> coords = current->Coord();
    respoly.push_back(coords);
    
    do {
      // find next Vertex / Vertices (path)
      if (current->EntryExit()==true)
      {
        //cout << "Intersection was Entry, so move to Next() on same polygon!" << endl;
        do {
          current = current->Next();
          //cout << "Current vertex is " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
          vector<double> coords = current->Coord();
          respoly.push_back(coords);
        } while (current->Intersect()==false);
        //cout << "Found intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      }
      else
      {
        //cout << "Intersection was Exit, so move to Prev() on same polygon!" << endl;
        do {
          current = current->Prev();
          //cout << "Current vertex is " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
          vector<double> coords = current->Coord();
          respoly.push_back(coords);
        } while (current->Intersect()==false);
        //cout << "Found intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      }
      
      // jump to the other input polygon
      current = current->Neighbor();
      //cout << "Jumping to other polygon at intersection: " << current->Coord()[0] << " " << current->Coord()[1] << " " << current->Coord()[2] << endl;
      //cout << "Length of result list so far: " << (int)respoly.size() << endl;
      
    } while (current!=&intersec1[0] && current!=&intersec2[0]);
    
    // check if last entry is identical to first entry
    double fldiff[3] = {0.0, 0.0, 0.0};
    bool identical = true;
    for (int k=0;k<3;++k)
    {
      fldiff[k] = respoly[(int)respoly.size()-1][k] - respoly[0][k];
      if (abs(fldiff[k]>1.0e-12)) identical = false;
    }
    // remove last entry if so, throw dserror if not so
    if (identical) respoly.pop_back();
    else dserror("ERROR: We did not arrive at the staring point again...?");
    
    
    // collapse respoly points that are very close
    vector<vector<double> > collapsedrespoly;
    for (int i=0;i<(int)respoly.size();++i)
    {
      // find distance between two consecutive points
      // first point of respoly
      if (i==0)
        collapsedrespoly.push_back(respoly[i]);
      
      // last point of respoly
      else if (i==(int)respoly.size()-1)
      {
        double diff[3] = {0.0, 0.0, 0.0};
        double diff2[3] = {0.0, 0.0, 0.0};
        
        for (int k=0;k<3;++k) diff[k] = respoly[i][k] - respoly[i-1][k];
        for (int k=0;k<3;++k) diff2[k] = respoly[0][k] - respoly[i][k];
        
        double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
        double dist2 = sqrt(diff2[0]*diff2[0]+diff2[1]*diff2[1]+diff2[2]*diff2[2]);
        double tolcollapse = 10*tol;
        
        if (abs(dist) >= tolcollapse && abs(dist2) >= tolcollapse)
          collapsedrespoly.push_back(respoly[i]);
        else {}
         // cout << "Collapsed two points in result polygon!" << endl;
      }
      
      // standard case
      else
      {
        double diff[3] = {0.0, 0.0, 0.0};
        for (int k=0;k<3;++k) diff[k] = respoly[i][k] - respoly[i-1][k];

        double dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
        double tolcollapse = 10*tol;
        
        if (abs(dist) >= tolcollapse)
          collapsedrespoly.push_back(respoly[i]);
        else {}
          //cout << "Collapsed two points in result polygon!" << endl;
      }
    }
    
    // replace respoly by collapsed respoly
    respoly = collapsedrespoly;
    //cout << "Final length of result list: " << (int)respoly.size() << endl;
        
    // check if respoly collapsed to nothing
    if ((int)collapsedrespoly.size()<3)
    {
     // cout << "Collapsing of result polygon led to < 3 vertices -> no respoly!" << endl;
      vector<vector<double> > empty(0,vector<double>(3));
      respoly = empty;
    }
 
    // check for rotation of result polygon (must be clockwise!!!)
    // first get geometric center
    double center[3] = {0.0, 0.0, 0.0};
    
    for (int i=0;i<(int)respoly.size();++i)
      for (int k=0;k<3;++k)
        center[k] += respoly[i][k]/((int)respoly.size());
    
    //cout << "\nCenter ResPoly: " << center[0] << " " << center[1] << " " << center[2] << endl;
    
    // then we compute the clockwise plane normal
    double diff[3] = {0.0, 0.0, 0.0};
    double edge[3] = {0.0, 0.0, 0.0};
    
    for (int k=0;k<3;++k)
    {
      diff[k] = respoly[0][k]-center[k];
      edge[k] = respoly[1][k]-respoly[0][k];
    }
    
    double cross[3] = {0.0, 0.0, 0.0};
    
    cross[0] = diff[1]*edge[2]-diff[2]*edge[1];
    cross[1] = diff[2]*edge[0]-diff[0]*edge[2];
    cross[2] = diff[0]*edge[1]-diff[1]*edge[0];

    // check against auxiliary plane normal
    double check = cross[0]*Auxn()[0]+cross[1]*Auxn()[1]+cross[2]*Auxn()[2];
    vector<vector<double> > newrespoly((int)respoly.size(),vector<double>(3));
    if (check<0)
    {
      // reorder result polygon in clockwise direction
     // cout << "Result polygon not ordered counter-clockwise -> reordered!" << endl;
      for (int i=0;i<(int)respoly.size();++i)
        newrespoly[(int)respoly.size()-1-i] = respoly[i];
      respoly = newrespoly;
    }
    
    // print final input polygons to screen
      //cout << "\nResult Poylgon:";
      //for (int i=0;i<(int)respoly.size();++i)
      //  cout << "\nVertex " << i << ":\t" << respoly[i][0] << "\t" << respoly[i][1] << "\t" << respoly[i][2];
      
    // check if result polygon is convex
    // a polygon is convex if the scalar product of an edge normal and the
    // next edge direction is negative for all edges
    for (int i=0;i<(int)respoly.size();++i)
    {
      // we need the edge vector first
      double edge[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        if (i!=(int)respoly.size()-1) edge[k] = respoly[i+1][k] - respoly[i][k];
        else edge[k] = respoly[0][k] - respoly[i][k];
      }
      // edge normal is result of cross product
      double n[3] = {0.0, 0.0, 0.0};
      n[0] = edge[1]*Auxn()[2]-edge[2]*Auxn()[1];
      n[1] = edge[2]*Auxn()[0]-edge[0]*Auxn()[2];
      n[2] = edge[0]*Auxn()[1]-edge[1]*Auxn()[0];
      
      // we need the next edge vector now
      double nextedge[3] = {0.0, 0.0, 0.0};
      for (int k=0;k<3;++k)
      {
        if (i<(int)respoly.size()-2) nextedge[k] = respoly[i+2][k] - respoly[i+1][k];
        else if (i==(int)respoly.size()-2) nextedge[k] = respoly[0][k] - respoly[i+1][k];
        else nextedge[k] = respoly[1][k] - respoly[0][k];
      }
      // check scalar product
      double check = n[0]*nextedge[0]+n[1]*nextedge[1]+n[2]*nextedge[2];
      if (check>0) dserror("ERROR: Result polygon not convex!");
    }
  }
  
  //cout << "\nRESULT POLYGON FINISHED!!!!!!!!!!!!!!!\n" << endl;

  /*
  // **********************************************************************
  // STEP6: Result visualization with GMSH
  // - plot the two input polygons and their vertex numbering
  // - plot the result polygon and its vertex numbering
  // **********************************************************************
  std::ostringstream filename;
  static int gmshcount=0;
  filename << "o/gmsh_output/" << "clipping_";
  if (gmshcount<10)
    filename << 0 << 0 << 0 << 0;
  else if (gmshcount<100)
    filename << 0 << 0 << 0;
  else if (gmshcount<1000)
    filename << 0 << 0;
  else if (gmshcount<10000)
    dserror("Gmsh output implemented for a maximum of 9.999 clip polygons");
  filename << gmshcount << ".pos";
  gmshcount++;

  // do output to file in c-style
  FILE* fp = NULL;
  fp = fopen(filename.str().c_str(), "w");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Clipping \" {" << endl;
        
  for (int i=0;i<(int)poly1.size();++i)
  {
    if (i!=(int)poly1.size()-1)
    {
      gmshfilecontent << "SL(" << scientific << poly1[i][0] << "," << poly1[i][1] << ","
                               << poly1[i][2] << "," << poly1[i+1][0] << "," << poly1[i+1][1] << ","
                               << poly1[i+1][2] << ")";
      gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
    }
    else
    {
      gmshfilecontent << "SL(" << scientific << poly1[i][0] << "," << poly1[i][1] << ","
                               << poly1[i][2] << "," << poly1[0][0] << "," << poly1[0][1] << ","
                               << poly1[0][2] << ")";
      gmshfilecontent << "{" << scientific << 1.0 << "," << 1.0 << "};" << endl;
      
    }
    gmshfilecontent << "T3(" << scientific << poly1[i][0] << "," << poly1[i][1] << "," << poly1[i][2] << "," << 17 << ")";
    gmshfilecontent << "{" << "S" << i << "};" << endl;
  }
  
  for (int i=0;i<(int)poly2.size();++i)
  {
    if (i!=(int)poly2.size()-1)
    {
      gmshfilecontent << "SL(" << scientific << poly2[i][0] << "," << poly2[i][1] << ","
                               << poly2[i][2] << "," << poly2[i+1][0] << "," << poly2[i+1][1] << ","
                               << poly2[i+1][2] << ")";
      gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
    }
    else
    {
      gmshfilecontent << "SL(" << scientific << poly2[i][0] << "," << poly2[i][1] << ","
                               << poly2[i][2] << "," << poly2[0][0] << "," << poly2[0][1] << ","
                               << poly2[0][2] << ")";
      gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
      
    }
    gmshfilecontent << "T3(" << scientific << poly2[i][0] << "," << poly2[i][1] << "," << poly2[i][2] << "," << 17 << ")";
    gmshfilecontent << "{" << "M" << i << "};" << endl;
  }
  
  for (int i=0;i<(int)respoly.size();++i)
  {
    if (i!=(int)respoly.size()-1)
    {
      gmshfilecontent << "SL(" << scientific << respoly[i][0] << "," << respoly[i][1] << ","
                               << respoly[i][2] << "," << respoly[i+1][0] << "," << respoly[i+1][1] << ","
                               << respoly[i+1][2] << ")";
      gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
    }
    else
    {
      gmshfilecontent << "SL(" << scientific << respoly[i][0] << "," << respoly[i][1] << ","
                               << respoly[i][2] << "," << respoly[0][0] << "," << respoly[0][1] << ","
                               << respoly[0][2] << ")";
      gmshfilecontent << "{" << scientific << 0.0 << "," << 0.0 << "};" << endl;
      
    }
    gmshfilecontent << "T3(" << scientific << respoly[i][0] << "," << respoly[i][1] << "," << respoly[i][2] << "," << 27 << ")";
    gmshfilecontent << "{" << "R" << i << "};" << endl;
  }
  
//  for (int i=0;i<(int)intersec1.size();++i)
//  {
//    gmshfilecontent << "T3(" << scientific << intersec1[i].Coord()[0] << "," << intersec1[i].Coord()[1] << "," << intersec1[i].Coord()[2] << "," << 17 << ")";
//    if (intersec1[i].EntryExit()==true && intersec2[i].EntryExit()==true) gmshfilecontent << "{" << "SEME" << "};" << endl;
//    else if (intersec1[i].EntryExit()==false && intersec2[i].EntryExit()==true) gmshfilecontent << "{" << "SXME" << "};" << endl;
//    else if (intersec1[i].EntryExit()==true && intersec2[i].EntryExit()==false) gmshfilecontent << "{" << "SEMX" << "};" << endl;
//    else gmshfilecontent << "{" << "SXMX" << "};" << endl;
//  }
  
  gmshfilecontent << "};" << endl;

  // move everything to gmsh post-processing file and close it
  fprintf(fp,gmshfilecontent.str().c_str());
  fclose(fp);
  */
  // return result
  return respoly;
}

/*----------------------------------------------------------------------*
 |  Triangulation of clip polygon in auxiliary plane          popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::Triangulation3D()
{
  // find geometric center of clipping polygon
  // as a first shot we use simple node averaging here
  vector<double> clipcenter(3);
  for (int k=0;k<3;++k) clipcenter[k] = 0.0;
  int clipsize = (int)(Clip().size());

  for (int i=0;i<clipsize;++i)
    for (int k=0;k<3;++k)
      clipcenter[k] += (Clip()[i][k] / clipsize);
 
  //cout << "-> This is a hack: WE NEED TO FIND GEOMETRIC CENTER OF CLIP POLYGON..." << endl;
  //cout << "Clipcenter (simple): " << clipcenter[0] << " " << clipcenter[1] << " " << clipcenter[2] << endl;
  
  // do triangulization (create Intcells)
  vector<Intcell> cells;
  
  // easy if clip polygon = triangle: 1 Intcell
  if (clipsize==3)
  {
    // Intcell vertices = clip polygon vertices
    Epetra_SerialDenseMatrix coords(3,clipsize);
    for (int i=0;i<clipsize;++i)
      for (int k=0;k<3;++k)
        coords(k,i) = Clip()[i][k];
    
    // create Intcell object and push back
    Cells().push_back(rcp(new Intcell(0,3,coords,DRT::Element::tri3)));    
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
        coords(k,1) = Clip()[num][k];
      }
      
      // the third vertex is the next vertex on clip polygon
      if (num==clipsize-1)
        for (int k=0;k<3;++k) coords(k,2) = Clip()[0][k];
      else
        for (int k=0;k<3;++k) coords(k,2) = Clip()[num+1][k];
      
      // create Intcell object and push back
      Cells().push_back(rcp(new Intcell(num,3,coords,DRT::Element::tri3)));
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Integration of cells in aux. plane (3D)                   popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling::IntegrateCells3D()
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
    // do the two integrations
    RCP<Epetra_SerialDenseMatrix> mseg = integrator.IntegrateMAuxPlane3D(sele_,mele_,Cells()[i],Auxn());
    RCP<Epetra_SerialDenseVector> gseg = integrator.IntegrateGAuxPlane3D(sele_,mele_,Cells()[i],Auxn());
  
    // compute directional derivative of M and store into nodes
    // if CONTACTONEMORTARLOOP defined, then DerivM does linearization of M AND D matrices !!!
    //integrator.DerivM(sele_,sxia,sxib,mele_,mxia,mxib);
      
    // do the two assemblies into the slave nodes
    // if CONTACTONEMORTARLOOP defined, then AssembleM does M AND D matrices !!!
    integrator.AssembleM(Comm(),sele_,mele_,*mseg);
    integrator.AssembleG(Comm(),sele_,*gseg);
  }
  return true;
}

#endif //#ifdef CCADISCRET
