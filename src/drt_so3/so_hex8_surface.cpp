/*!----------------------------------------------------------------------
\file so_hex8_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "so_hex8.H"
#include "../discret/linalg_utils.H"
#include "../discret/drt_utils.H"
#include "../discret/drt_discret.H"
#include "../discret/drt_dserror.H"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../discret/dstrc.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Surface::Soh8Surface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::So_hex8* parent,
                              const int lsurface) :
DRT::Element(id,element_soh8surface,owner),
parent_(parent),
lsurface_(lsurface)
{
  DSTraceHelper dst("Soh8Surface::Soh8Surface");
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Surface::Soh8Surface(const DRT::Elements::Soh8Surface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  DSTraceHelper dst("Soh8Surface::Soh8Surface");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Soh8Surface::Clone() const
{
  DSTraceHelper dst("Soh8Surface::Clone");
  DRT::Elements::Soh8Surface* newelement = new DRT::Elements::Soh8Surface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::Soh8Surface::Shape() const
{
  return quad4;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Surface::Pack(vector<char>& data) const
{
  DSTraceHelper dst("Soh8Surface::Pack");
  data.resize(0);
  dserror("this Soh8Surface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Surface::Unpack(const vector<char>& data)
{
  DSTraceHelper dst("Soh8Surface::Unpack");
  dserror("this Soh8Surface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::Soh8Surface::~Soh8Surface()
{
  DSTraceHelper dst("Soh8Surface::~Soh8Surface");
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Soh8Surface::Print(ostream& os) const
{
  DSTraceHelper dst("Soh8Surface::Print");
  os << "Soh8Surface ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     maf 04/07*
 * ---------------------------------------------------------------------*/
int DRT::Elements::Soh8Surface::EvaluateNeumann(ParameterList&           params,
                                                DRT::Discretization&     discretization,
                                                DRT::Condition&          condition,
                                                vector<int>&             lm,
                                                Epetra_SerialDenseVector& elevec1)
{
  DSTraceHelper dst("Soh8Surface::EvaluateNeumann");

  // get values and switches from the condition
  vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0]; 
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    dyn_facfromcurve(curvenum,time,&curvefac); 
  // **

  // element geometry
  const int numnod = 4;
  Epetra_SerialDenseMatrix xsrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  for (int i=0; i<numnod; ++i){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];
  }
  
  /*
  ** Here, we integrate a 4-node surface with 2x2 Gauss Points
  */
  const int ngp = 4;
  
  // gauss parameters
  const double gpweight = 1.0;
  const double gploc    = 1.0/sqrt(3.0);
  Epetra_SerialDenseMatrix gpcoord (ngp,2);
  gpcoord(0,0) = - gploc;
  gpcoord(0,1) = - gploc;
  gpcoord(1,0) =   gploc;
  gpcoord(1,1) = - gploc;
  gpcoord(2,0) = - gploc;
  gpcoord(2,1) =   gploc;
  gpcoord(3,0) =   gploc;
  gpcoord(3,1) =   gploc;

  for (int gpid = 0; gpid < ngp; ++gpid) {    // loop over intergration points
    // get shape functions and derivatives of element surface
    vector<double> funct(ngp);                // 4 shape function values
    double drs;                               // surface area factor
    soh8_surface_integ(&funct,&drs,&xsrefe,gpcoord(gpid,0),gpcoord(gpid,1));
    double fac = gpweight * drs * curvefac;   // integration factor

    // distribute over element load vector
    for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
      for(int dim=0; dim<NUMDIM_SOH8; dim++) {
        elevec1[nodid+dim] += funct[nodid] * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }
  }
  return 0;
}
                                                   
/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      maf 05/07*
 * ---------------------------------------------------------------------*/
void DRT::Elements::Soh8Surface::soh8_surface_integ(
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      const Epetra_SerialDenseMatrix* xsrefe,// (i) material element coords
      const double r,                        // (i) coord in r-direction 
      const double s)                        // (i) coord in s-direction
{
  DSTraceHelper dst("Soh8Surface::soh8_surface_metric");
  
  // shape functions for 4 nodes
  (*funct)[0] = 0.25 * (1.0-r) * (1.0-s);
  (*funct)[1] = 0.25 * (1.0+r) * (1.0-s);
  (*funct)[2] = 0.25 * (1.0+r) * (1.0+s);
  (*funct)[3] = 0.25 * (1.0-r) * (1.0+s);
  // derivatives of 4 shape functions wrt 2 directions
  Epetra_SerialDenseMatrix deriv(4,2);
  deriv(0,0) = -0.25 * (1.0-s);
  deriv(0,1) = -0.25 * (1.0-r);
  deriv(1,0) =  0.25 * (1.0-s);
  deriv(1,1) = -0.25 * (1.0+r);
  deriv(2,0) =  0.25 * (1.0+s);
  deriv(2,1) =  0.25 * (1.0+r);
  deriv(3,0) = -0.25 * (1.0+s);
  deriv(3,1) =  0.25 * (1.0-r);
                        
  // compute dXYZ / drs
  Epetra_SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('T','N',1.0,deriv,xsrefe,1.0);
  
  /* compute covariant metric tensor G for surface element
  **                        | g11   g12 |
  **                    G = |           |
  **                        | g12   g22 |
  ** where (o denotes the inner product, xyz a vector)
  **
  **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
  ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
  **        dr     dr            dr     ds            ds     ds
  */
  Epetra_SerialDenseMatrix metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,1.0);
  (*sqrtdetg) = sqrt( metrictensor(0,0)*metrictensor(1,1)
                     -metrictensor(0,1)*metrictensor(1,0));
  return;
}                     
  

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
