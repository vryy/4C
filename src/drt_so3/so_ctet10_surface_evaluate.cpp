/*!----------------------------------------------------------------------*###
\file so_ctet10_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by : Alexander Volf
			alexander.volf@mytum.de  
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET
#ifdef CCADISCRET

#include "so_ctet10.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     vlf 04/07*
 * ---------------------------------------------------------------------*/
int DRT::Elements::Soctet10Surface::EvaluateNeumann(ParameterList&           params,
                                                DRT::Discretization&     discretization,
                                                DRT::Condition&          condition,
                                                vector<int>&             lm,
                                                Epetra_SerialDenseVector& elevec1)
{
  //cout << "DRT::Elements::Soctet10Surface::EvaluateNeumann" << endl;
  //getchar();
  // OBSOLETE!!! needs change
  Epetra_SerialDenseMatrix* shapefct;
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOTET10_FACE]
  sotet4_surface_shapefunc(&shapefct,&weights);
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // element geometry
  const int numnod = 6;
  Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOCTET10);  // material coord. of element
  for (int i=0; i<numnod; i++){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];
  }

  SUB_STRUCTURE_SURF sub(xsrefe);
  
  for (int num_ele = 0; num_ele < 4; num_ele++)
  {  
    /*
    ** Here, we integrate a 6-node surface with 3 Gauss Points
    */
    double fac = (*weights)(0) * (sub.my_surface[num_ele]).my_area() * curvefac;   // integration factor
 
    // gauss parameters
    for (int gpid = 0; gpid < 1; gpid++) {    // loop over intergration points
      // get shape functions and derivatives of element surface
      // distribute over element load vector
      for (int nodid=0; nodid < 3; nodid++) {
         for(int dim=0; dim < 3; dim++) {
           elevec1[((sub.my_surface[num_ele]).my_nodes[nodid]).global_id * NUMDIM_SOCTET10 + dim] +=\
        	  (*shapefct)(nodid,gpid) * (*onoff)[dim] * (*val)[dim] * fac;
         }
      }
  
    } 
  }  
  return 0;
} //Sotet10Surface::EvaluateNeumann(..)


DRT::Elements::Soctet10Surface::SUB_NODE::SUB_NODE()
{
	//nothing to do
}

void DRT::Elements::Soctet10Surface::SUB_NODE::init(
	const int in_local_id,
	const int in_global_id,
	const Epetra_SerialDenseMatrix& xrefe)
{	
	local_id =in_local_id;
  	global_id=in_global_id;
  	
  	my_x[0]=xrefe(global_id,0);
	my_x[1]=xrefe(global_id,1);
	my_x[2]=xrefe(global_id,2);
}

DRT::Elements::Soctet10Surface::SUB_NODE::~SUB_NODE()
{
	//nothing to do
}

DRT::Elements::Soctet10Surface::TRI3_SUB::TRI3_SUB()
{
	//nothing to do
}

void DRT::Elements::Soctet10Surface::TRI3_SUB::init(
	const int& node1,
	const int& node2,
	const int& node3,
	const Epetra_SerialDenseMatrix& xrefe)
{
	my_nodes[0].init(0,node1,xrefe);
	my_nodes[1].init(1,node2,xrefe);
	my_nodes[2].init(2,node3,xrefe);
}

double DRT::Elements::Soctet10Surface::TRI3_SUB::my_area()
{
    Epetra_SerialDenseMatrix xsrefe(3,NUMDIM_SOCTET10);  // material coord. of element
    for (int i=0; i<3; i++){
      xsrefe(i,0) = my_nodes[i].my_x[0];
      xsrefe(i,1) = my_nodes[i].my_x[1];
      xsrefe(i,2) = my_nodes[i].my_x[2];
    }
   
    Epetra_SerialDenseVector A(NUMDIM_SOCTET10);
    Epetra_SerialDenseVector B(NUMDIM_SOCTET10);
    Epetra_SerialDenseVector C(NUMDIM_SOCTET10);
  
    /*
     * to compute the Jacobian i compute the area of the triangle 
     * for that i get the boundary vectors A,B of the triangle 
     * the area is then |A x B| 
     */
  
    A(0)=xsrefe(1,0)-xsrefe(0,0);
    A(1)=xsrefe(1,1)-xsrefe(0,1);
    A(2)=xsrefe(1,2)-xsrefe(0,2);
  
    B(0)=xsrefe(2,0)-xsrefe(0,0);
    B(1)=xsrefe(2,1)-xsrefe(0,1);
    B(2)=xsrefe(2,2)-xsrefe(0,2);
  
    /* C = A x B  */
    C(0)=A(0)*B(1) - A(1)*B(0);
    C(1)=A(1)*B(2) - A(2)*B(1);
    C(2)=A(2)*B(0) - A(0)*B(2);
  
    /* return = |A x B|*/
    return C.Norm2()/2;
}

DRT::Elements::Soctet10Surface::TRI3_SUB::~TRI3_SUB()
{
	//nothing to do
}

DRT::Elements::Soctet10Surface::SUB_STRUCTURE_SURF::SUB_STRUCTURE_SURF(const Epetra_SerialDenseMatrix& xrefe)
{
	my_surface[ 0].init( 0, 3, 5, xrefe);
	my_surface[ 1].init( 3, 1, 4, xrefe);
	my_surface[ 2].init( 5, 4, 2, xrefe);
	my_surface[ 3].init( 3, 4, 5, xrefe);
	
}

DRT::Elements::Soctet10Surface::SUB_STRUCTURE_SURF::~SUB_STRUCTURE_SURF()
{
	//nothing to do
}
/*----------------------------------------------------------------------*
 * Get shape functions for a tet4 face					       vlf 08/07*
 * ---------------------------------------------------------------------*/
void DRT::Elements::Soctet10Surface::sotet4_surface_shapefunc(
      Epetra_SerialDenseMatrix** shapefct,  // pointer to pointer of shapefct
      Epetra_SerialDenseVector** weights)   // pointer to pointer of weights
{

  static Epetra_SerialDenseMatrix  f(3,1);  // shape functions
  static Epetra_SerialDenseVector weightfactors(1);   // weights for each gp

 //Quadrature rule from Carlos A. Felippa: Adv. FEM  §17 
  const double gploc_alpha    = (double)1/3;    // gp sampling point value for liner. fct
  const double w			  = (double)1;

  const double ksi1[1] = {gploc_alpha };
  const double ksi2[1] = {gploc_alpha };
  const double ksi3[1] = {gploc_alpha };
  
  for (int i=0; i<1; i++) {
      f(0,i) = ksi1[i];
      f(1,i) = ksi2[i];
      f(2,i) = ksi3[i];
      weightfactors[i] = w; // just for clarity how to get weight factors
   } 
   
   *weights  = &weightfactors;
   *shapefct = &f;

   return;

}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOTET
