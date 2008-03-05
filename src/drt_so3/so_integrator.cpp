/*!----------------------------------------------------------------------
\file so_tet10.H

\class DRT::ELEMENTS::So_tet10

\brief A C++ version of the 10-node tet solid element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by: Alexander Volf
			alexander.volf@mytum.de     
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3

#include "so_integrator.H"
#include "so_tet4.H"
#include "so_tet10.H"
#include "so_ctet10.H"
#include "so_hex8.H"
#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_utils_integration.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"

/*----------------------------------------------------------------------*
 | destructor for the Discrete_integrator class               volf 09/07|
 | makes sure the arrays shapefct_gp and deriv_gp are correctly			|
 | destroyed															|				
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_integrator::~So_integrator()
{
	return;
}

/*----------------------------------------------------------------------*
 |  shape functions and derivatives for So_hex8                vlf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::Integrator_So_hex8::Integrator_So_hex8()
{

    num_gp = NUMGPT_SOH8;
    num_nodes = NUMNOD_SOH8 ;
    num_coords = 3;
    shapefct_gp.resize(NUMGPT_SOH8);
    deriv_gp.resize(NUMGPT_SOH8);
    weights.Size(NUMGPT_SOH8);
    
#if 1 //OUR SHAPE FUNCTIONS

  #if 1 // long double GP
    const long double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
    const long double gpw      = 1.0;              // weight at every gp for linear fct

    // (r,s,t) gp-locations of fully integrated linear 8-node Hex
    const long double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
    const long double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
    const long double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};
    const long double w[NUMGPT_SOH8] = {   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw};

    // fill up nodal f at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
    (shapefct_gp[i]).Size(NUMNOD_SOH8);    
      (shapefct_gp[i])(0) = ((long double)1.0-r[i])*((long double)1.0-s[i])*((long double)1.0-t[i])*(long double)0.125;
      (shapefct_gp[i])(1) = ((long double)1.0+r[i])*((long double)1.0-s[i])*((long double)1.0-t[i])*(long double)0.125;
      (shapefct_gp[i])(2) = ((long double)1.0+r[i])*((long double)1.0+s[i])*((long double)1.0-t[i])*(long double)0.125;
      (shapefct_gp[i])(3) = ((long double)1.0-r[i])*((long double)1.0+s[i])*((long double)1.0-t[i])*(long double)0.125;
      (shapefct_gp[i])(4) = ((long double)1.0-r[i])*((long double)1.0-s[i])*((long double)1.0+t[i])*(long double)0.125;
      (shapefct_gp[i])(5) = ((long double)1.0+r[i])*((long double)1.0-s[i])*((long double)1.0+t[i])*(long double)0.125;
      (shapefct_gp[i])(6) = ((long double)1.0+r[i])*((long double)1.0+s[i])*((long double)1.0+t[i])*(long double)0.125;
      (shapefct_gp[i])(7) = ((long double)1.0-r[i])*((long double)1.0+s[i])*((long double)1.0+t[i])*(long double)0.125;
      weights[i] = w[i]*w[i]*w[i]; // just for clarity how to get weight factors
    }

    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
      (deriv_gp[i]).Shape(3,NUMNOD_SOH8);
        // df wrt to r "+0" for each node(0..7) at each gp [i]
        (deriv_gp[i])(0,0) = -((long double)1.0-s[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(0,1) =  ((long double)1.0-s[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(0,2) =  ((long double)1.0+s[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(0,3) = -((long double)1.0+s[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(0,4) = -((long double)1.0-s[i])*((long double)1.0+t[i])*(long double)0.125;
        (deriv_gp[i])(0,5) =  ((long double)1.0-s[i])*((long double)1.0+t[i])*(long double)0.125;
        (deriv_gp[i])(0,6) =  ((long double)1.0+s[i])*((long double)1.0+t[i])*(long double)0.125;
        (deriv_gp[i])(0,7) = -((long double)1.0+s[i])*((long double)1.0+t[i])*(long double)0.125;

        // df wrt to s "+1" for each node(0..7) at each gp [i]
        (deriv_gp[i])(1,0) = -((long double)1.0-r[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(1,1) = -((long double)1.0+r[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(1,2) =  ((long double)1.0+r[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(1,3) =  ((long double)1.0-r[i])*((long double)1.0-t[i])*(long double)0.125;
        (deriv_gp[i])(1,4) = -((long double)1.0-r[i])*((long double)1.0+t[i])*(long double)0.125;
        (deriv_gp[i])(1,5) = -((long double)1.0+r[i])*((long double)1.0+t[i])*(long double)0.125;
        (deriv_gp[i])(1,6) =  ((long double)1.0+r[i])*((long double)1.0+t[i])*(long double)0.125;
        (deriv_gp[i])(1,7) =  ((long double)1.0-r[i])*((long double)1.0+t[i])*(long double)0.125;

        // df wrt to t "+2" for each node(0..7) at each gp [i]
        (deriv_gp[i])(2,0) = -((long double)1.0-r[i])*((long double)1.0-s[i])*(long double)0.125;
        (deriv_gp[i])(2,1) = -((long double)1.0+r[i])*((long double)1.0-s[i])*(long double)0.125;
        (deriv_gp[i])(2,2) = -((long double)1.0+r[i])*((long double)1.0+s[i])*(long double)0.125;
        (deriv_gp[i])(2,3) = -((long double)1.0-r[i])*((long double)1.0+s[i])*(long double)0.125;
        (deriv_gp[i])(2,4) =  ((long double)1.0-r[i])*((long double)1.0-s[i])*(long double)0.125;
        (deriv_gp[i])(2,5) =  ((long double)1.0+r[i])*((long double)1.0-s[i])*(long double)0.125;
        (deriv_gp[i])(2,6) =  ((long double)1.0+r[i])*((long double)1.0+s[i])*(long double)0.125;
        (deriv_gp[i])(2,7) =  ((long double)1.0-r[i])*((long double)1.0+s[i])*(long double)0.125;
    }
  
  #else // double GP 
  const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
    const double gpw      = 1.0;              // weight at every gp for linear fct

    // (r,s,t) gp-locations of fully integrated linear 8-node Hex
    const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
    const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
    const double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};
    const double w[NUMGPT_SOH8] = {   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw};
    
    // fill up nodal f at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
    (shapefct_gp[i]).Size(NUMNOD_SOH8);    
      (shapefct_gp[i])(0) = (1.0-r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
      (shapefct_gp[i])(1) = (1.0+r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
      (shapefct_gp[i])(2) = (1.0+r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
      (shapefct_gp[i])(3) = (1.0-r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
      (shapefct_gp[i])(4) = (1.0-r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
      (shapefct_gp[i])(5) = (1.0+r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
      (shapefct_gp[i])(6) = (1.0+r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
      (shapefct_gp[i])(7) = (1.0-r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
      weights[i] = w[i]*w[i]*w[i]; // just for clarity how to get weight factors
    }

    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    for (int i=0; i<NUMGPT_SOH8; ++i) {
      (deriv_gp[i]).Shape(3,NUMNOD_SOH8);
        // df wrt to r "+0" for each node(0..7) at each gp [i]
        (deriv_gp[i])(0,0) = -(1.0-s[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(0,1) =  (1.0-s[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(0,2) =  (1.0+s[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(0,3) = -(1.0+s[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(0,4) = -(1.0-s[i])*(1.0+t[i])*0.125;
        (deriv_gp[i])(0,5) =  (1.0-s[i])*(1.0+t[i])*0.125;
        (deriv_gp[i])(0,6) =  (1.0+s[i])*(1.0+t[i])*0.125;
        (deriv_gp[i])(0,7) = -(1.0+s[i])*(1.0+t[i])*0.125;

        // df wrt to s "+1" for each node(0..7) at each gp [i]
        (deriv_gp[i])(1,0) = -(1.0-r[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(1,1) = -(1.0+r[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(1,2) =  (1.0+r[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(1,3) =  (1.0-r[i])*(1.0-t[i])*0.125;
        (deriv_gp[i])(1,4) = -(1.0-r[i])*(1.0+t[i])*0.125;
        (deriv_gp[i])(1,5) = -(1.0+r[i])*(1.0+t[i])*0.125;
        (deriv_gp[i])(1,6) =  (1.0+r[i])*(1.0+t[i])*0.125;
        (deriv_gp[i])(1,7) =  (1.0-r[i])*(1.0+t[i])*0.125;

        // df wrt to t "+2" for each node(0..7) at each gp [i]
        (deriv_gp[i])(2,0) = -(1.0-r[i])*(1.0-s[i])*0.125;
        (deriv_gp[i])(2,1) = -(1.0+r[i])*(1.0-s[i])*0.125;
        (deriv_gp[i])(2,2) = -(1.0+r[i])*(1.0+s[i])*0.125;
        (deriv_gp[i])(2,3) = -(1.0-r[i])*(1.0+s[i])*0.125;
        (deriv_gp[i])(2,4) =  (1.0-r[i])*(1.0-s[i])*0.125;
        (deriv_gp[i])(2,5) =  (1.0+r[i])*(1.0-s[i])*0.125;
        (deriv_gp[i])(2,6) =  (1.0+r[i])*(1.0+s[i])*0.125;
        (deriv_gp[i])(2,7) =  (1.0-r[i])*(1.0+s[i])*0.125;
    }
    #endif //#if // long double GP
    
#else  // DRT SHAPE FUNCTIONS
  // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
  // fill up nodal f at each gp
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule_hex_8point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule_);
 
  for (int igp = 0; igp < intpoints.nquad; ++igp) {
    const double r = intpoints.qxg[igp][0];
    const double s = intpoints.qxg[igp][1];
    const double t = intpoints.qxg[igp][2];

    shapefct_gp[igp].Size(NUMNOD_SOH8);
    deriv_gp[igp].Shape(NUMDIM_SOH8, NUMNOD_SOH8);
    DRT::UTILS::shape_function_3D(shapefct_gp[igp], r, s, t, hex8);
    DRT::UTILS::shape_function_3D_deriv1(deriv_gp[igp], r, s, t, hex8);
    weights[igp] = intpoints.qwgt[igp];
 }
#endif  // #if 0 //OUR SHAPE FUNCTIONS

#if 0 // COMPARE DRT SHAPE FUNCTIONS
  // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
  // fill up nodal f at each gp
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule_hex_8point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule_);
 
  for (int igp = 0; igp < intpoints.nquad; ++igp) {
    const double r = intpoints.qxg[igp][0];
    const double s = intpoints.qxg[igp][1];
    const double t = intpoints.qxg[igp][2];

    Epetra_SerialDenseVector shapefct_gpmy(NUMNOD_SOH8);
    Epetra_SerialDenseVector temp_vec(NUMNOD_SOH8);
    Epetra_SerialDenseMatrix deriv_gpmy(NUMDIM_SOH8, NUMNOD_SOH8);
    Epetra_SerialDenseMatrix temp_mat(NUMDIM_SOH8, NUMNOD_SOH8);
    DRT::UTILS::shape_function_3D(shapefct_gpmy, r, s, t, hex8);
    DRT::UTILS::shape_function_3D_deriv1(deriv_gpmy, r, s, t, hex8);
    
    temp_vec+=shapefct_gp[igp];
    temp_vec.Scale(-1);
    temp_vec+=shapefct_gpmy;
    
    temp_mat+=deriv_gp[igp];
    temp_mat.Scale(-1);
    temp_mat+=deriv_gpmy;
    
    cout << temp_vec;
    cout << temp_mat;
    getchar();
  }
#endif //#if 0 // COMPARE DRT SHAPE FUNCTIONS

  return;
}  // of DRT::ELEMENTS::So_hex8::Integrator_So_hex8::Integrator_So_hex8


/*----------------------------------------------------------------------*
 |  shape functions and derivatives for SoDisp                 vlf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDisp::Integrator_SoDisp::Integrator_SoDisp(DRT::ELEMENTS::SoDisp& this_element)
{

  const DRT::Element::DiscretizationType distype = this_element.Shape();
    // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
    // fill up nodal f at each gp
    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    // forward initialization of necessary attributes

    const DRT::UTILS::IntegrationPoints3D intpoints = \
        DRT::UTILS::getIntegrationPoints3D(this_element.gaussrule_);
    
    num_gp = this_element.numgpt_disp_;
    num_nodes = this_element.numnod_disp_ ;
    num_coords = 3;
    shapefct_gp.resize(num_gp);
    deriv_gp.resize(num_gp);
    weights.Size(this_element.numgpt_disp_);
    
    for (int igp = 0; igp < intpoints.nquad; ++igp) {
      const double r = intpoints.qxg[igp][0];
      const double s = intpoints.qxg[igp][1];
      const double t = intpoints.qxg[igp][2];

      Epetra_SerialDenseVector funct(num_nodes);
      Epetra_SerialDenseMatrix deriv(num_coords, num_nodes);
      DRT::UTILS::shape_function_3D(funct, r, s, t, distype);
      DRT::UTILS::shape_function_3D_deriv1(deriv, r, s, t, distype);
        
      // return adresses of just evaluated matrices
      shapefct_gp[igp] = funct;
      deriv_gp[igp]    = deriv;
      weights[igp] = intpoints.qwgt[igp]; // return adress of static object to target of pointer
    }

    return;
}


DRT::ELEMENTS::Integrator_tet4_1point::Integrator_tet4_1point(void)
{
	 // forward initialization of necessary attributes
  num_gp = NUMGPT_SOTET4;
  num_nodes = NUMNOD_SOTET4;
  num_coords = NUMCOORD_SOTET4;
  shapefct_gp.resize(NUMGPT_SOTET4);
  deriv_gp.resize(NUMGPT_SOTET4);
  weights.Size(num_gp);
  
  //Quadrature rule from Carlos A. Felippa: Adv. FEM  §16.4 
  const double gploc_alpha    =  0.25;    // gp sampling point value for linear. fct
  const double gpw = 1.0;
  // (xsi1, xsi2, xsi3 ,xsi4) gp-locations of fully integrated linear 4-node Tet
  const double xsi1[NUMGPT_SOTET4] = {gploc_alpha};
  const double xsi2[NUMGPT_SOTET4] = {gploc_alpha};
  const double xsi3[NUMGPT_SOTET4] = {gploc_alpha};
  const double xsi4[NUMGPT_SOTET4] = {gploc_alpha};
  const double w[NUMGPT_SOTET4]    = {gpw};
	
   // fill up nodal f at each gp 
  for (int gp=0; gp<num_gp; gp++) {
      (shapefct_gp[gp]).Size(num_nodes);
      (shapefct_gp[gp])[0] = xsi1[gp];
      (shapefct_gp[gp])[1] = xsi2[gp];
      (shapefct_gp[gp])[2] = xsi3[gp];
      (shapefct_gp[gp])[3] = xsi4[gp];
      weights[gp] = w[gp]; 	// just for clarity how to get weight factors    
  }
 
  // fill up df xsi1, xsi2, xsi3, xsi4 directions (NUMDIM) at each gp
  for (int gp=0; gp<num_gp; gp++) {
  	(deriv_gp[gp]).Shape(num_nodes,num_coords);
  	
  	// deriv_gp wrt to xsi1 "(0,..)" for each node(0..9) at each gp [i]
	(deriv_gp[gp])(0,0) = 1;
	(deriv_gp[gp])(1,0) = 0;
    (deriv_gp[gp])(2,0) = 0;
    (deriv_gp[gp])(3,0) = 0;
    
    // deriv_gp wrt to xsi2 "(1,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,1) = 0;
    (deriv_gp[gp])(1,1) = 1;
    (deriv_gp[gp])(2,1) = 0;
    (deriv_gp[gp])(3,1) = 0;
     	
    // deriv_gp wrt to xsi3 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,2) = 0;
    (deriv_gp[gp])(1,2) = 0;
    (deriv_gp[gp])(2,2) = 1;
    (deriv_gp[gp])(3,2) = 0;
      
    // deriv_gp wrt to xsi4 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,3) = 0;
    (deriv_gp[gp])(1,3) = 0;
    (deriv_gp[gp])(2,3) = 0;
    (deriv_gp[gp])(3,3) = 1;
  }
}

DRT::ELEMENTS::Integrator_tet4_4point::Integrator_tet4_4point(void)
{
  // forward initialization of necessary attributes
  num_gp = 4;
  num_nodes = 4;
  num_coords = 4;
  shapefct_gp.resize(num_gp);
  deriv_gp.resize(num_gp);
  weights.Size(num_gp);
  
 //Quadrature rule from Carlos A. Felippa: Adv. FEM  §16.4 
  double gploc_alpha    = (5.0 + 3.0*sqrt(5.0))/20.0;    // gp sampling point value for quadr. fct
  double gploc_beta     = (5.0 - sqrt(5.0))/20.0;
  double gpw      = 0.25;              // weight at every gp for linear fct
 
  // (xsi1, xsi2, xsi3 ,xsi4) gp-locations of fully integrated linear 10-node Tet
  const double xsi1[4] = {gploc_alpha, gploc_beta , gploc_beta , gploc_beta };
  const double xsi2[4] = {gploc_beta , gploc_alpha, gploc_beta , gploc_beta };
  const double xsi3[4] = {gploc_beta , gploc_beta , gploc_alpha, gploc_beta };
  const double xsi4[4] = {gploc_beta , gploc_beta , gploc_beta , gploc_alpha};
  const double w[4]    = {   gpw,   gpw,   gpw,   gpw};
  
   // fill up nodal f at each gp 
  for (int gp=0; gp<num_gp; gp++) {
      (shapefct_gp[gp]).Size(num_nodes);
      (shapefct_gp[gp])[0] = xsi1[gp];
      (shapefct_gp[gp])[1] = xsi2[gp];
      (shapefct_gp[gp])[2] = xsi3[gp];
      (shapefct_gp[gp])[3] = xsi4[gp];
      weights[gp] = w[gp]; 	// just for clarity how to get weight factors    
  }
 
  // fill up df xsi1, xsi2, xsi3, xsi4 directions (NUMDIM) at each gp
  for (int gp=0; gp<num_gp; gp++) {
  	(deriv_gp[gp]).Shape(num_nodes,num_coords);
  	
  	// deriv_gp wrt to xsi1 "(0,..)" for each node(0..9) at each gp [i]
	(deriv_gp[gp])(0,0) = 1;
	(deriv_gp[gp])(1,0) = 0;
    (deriv_gp[gp])(2,0) = 0;
    (deriv_gp[gp])(3,0) = 0;
    
    // deriv_gp wrt to xsi2 "(1,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,1) = 0;
    (deriv_gp[gp])(1,1) = 1;
    (deriv_gp[gp])(2,1) = 0;
    (deriv_gp[gp])(3,1) = 0;
     	
    // deriv_gp wrt to xsi3 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,2) = 0;
    (deriv_gp[gp])(1,2) = 0;
    (deriv_gp[gp])(2,2) = 1;
    (deriv_gp[gp])(3,2) = 0;
      
    // deriv_gp wrt to xsi4 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,3) = 0;
    (deriv_gp[gp])(1,3) = 0;
    (deriv_gp[gp])(2,3) = 0;
    (deriv_gp[gp])(3,3) = 1;
  }
}

/*----------------------------------------------------------------------*
 | constructor for a integrator class              			  volf 09/07|
 | uses shape functions of a quadratic tetrahedra using so-called       |
 | "natural coordinates" as described by Carlos A. Felippa in Adv. FEM  |
 | Aerospace Engineering Sciences - University of Colorado at Boulder   |
 *----------------------------------------------------------------------*/
 
DRT::ELEMENTS::Integrator_tet10_4point::Integrator_tet10_4point(void)
{
  // forward initialization of necessary attributes
  num_gp = NUMGPT_SOTET10;
  num_nodes = NUMNOD_SOTET10 ;
  num_coords = NUMCOORD_SOTET10;
  shapefct_gp.resize(4);
  deriv_gp.resize(4);
  weights.Size(num_gp);

  //Quadrature rule from Carlos A. Felippa: Adv. FEM  §16.4 
  const double gploc_alpha    = (5.0 + 3.0*sqrt(5.0))/20.0;    // gp sampling point value for quadr. fct
  const double gploc_beta     = (5.0 - sqrt(5.0))/20.0;
  const double gpw      = 0.25;              // weight at every gp for linear fct
 
  // (xsi1, xsi2, xsi3 ,xsi4) gp-locations of fully integrated linear 10-node Tet
  const double xsi1[NUMGPT_SOTET10] = {gploc_alpha, gploc_beta , gploc_beta , gploc_beta };
  const double xsi2[NUMGPT_SOTET10] = {gploc_beta , gploc_alpha, gploc_beta , gploc_beta };
  const double xsi3[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_alpha, gploc_beta };
  const double xsi4[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_beta , gploc_alpha};
  const double w[NUMGPT_SOTET10]    = {   gpw,   gpw,   gpw,   gpw};
	
   // fill up nodal f at each gp 
  for (int gp=0; gp<num_gp; gp++) {
      (shapefct_gp[gp]).Size(num_nodes);
      (shapefct_gp[gp])[0] = xsi1[gp] * (2*xsi1[gp] -1);
      (shapefct_gp[gp])[1] = xsi2[gp] * (2*xsi2[gp] -1);
      (shapefct_gp[gp])[2] = xsi3[gp] * (2*xsi3[gp] -1);
      (shapefct_gp[gp])[3] = xsi4[gp] * (2*xsi4[gp] -1);
      (shapefct_gp[gp])[4] = 4 * xsi1[gp] * xsi2[gp];
      (shapefct_gp[gp])[5] = 4 * xsi2[gp] * xsi3[gp];
      (shapefct_gp[gp])[6] = 4 * xsi3[gp] * xsi1[gp];
      (shapefct_gp[gp])[7] = 4 * xsi1[gp] * xsi4[gp];  
      (shapefct_gp[gp])[8] = 4 * xsi2[gp] * xsi4[gp];
      (shapefct_gp[gp])[9] = 4 * xsi3[gp] * xsi4[gp];
      weights[gp] = w[gp]; 	// just for clarity how to get weight factors    
  }
 
  // fill up df xsi1, xsi2, xsi3, xsi4 directions (NUMDIM) at each gp
  for (int gp=0; gp<num_gp; gp++) {
  	(deriv_gp[gp]).Shape(num_nodes,num_coords);
  	// deriv_gp wrt to xsi1 "(0,..)" for each node(0..9) at each gp [i]
	(deriv_gp[gp])(0,0) = 4 * xsi1[gp]-1;
	(deriv_gp[gp])(1,0) = 0;
    (deriv_gp[gp])(2,0) = 0;
    (deriv_gp[gp])(3,0) = 0;
    
    (deriv_gp[gp])(4,0) = 4 * xsi2[gp];
    (deriv_gp[gp])(5,0) = 0;
    (deriv_gp[gp])(6,0) = 4 * xsi3[gp];
    (deriv_gp[gp])(7,0) = 4 * xsi4[gp];  
    (deriv_gp[gp])(8,0) = 0;
    (deriv_gp[gp])(9,0) = 0;

    // deriv_gp wrt to xsi2 "(1,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,1) = 0;
    (deriv_gp[gp])(1,1) = 4 * xsi2[gp] - 1;
    (deriv_gp[gp])(2,1) = 0;
    (deriv_gp[gp])(3,1) = 0;
     	
    (deriv_gp[gp])(4,1) = 4 * xsi1[gp];
    (deriv_gp[gp])(5,1) = 4 * xsi3[gp];
    (deriv_gp[gp])(6,1) = 0;
    (deriv_gp[gp])(7,1) = 0;  
    (deriv_gp[gp])(8,1) = 4 * xsi4[gp];
    (deriv_gp[gp])(9,1) = 0;

    // deriv_gp wrt to xsi3 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,2) = 0;
    (deriv_gp[gp])(1,2) = 0;
    (deriv_gp[gp])(2,2) = 4 * xsi3[gp] - 1;
    (deriv_gp[gp])(3,2) = 0;
      
    (deriv_gp[gp])(4,2) = 0;
    (deriv_gp[gp])(5,2) = 4 * xsi2[gp];
    (deriv_gp[gp])(6,2) = 4 * xsi1[gp];
    (deriv_gp[gp])(7,2) = 0;  
    (deriv_gp[gp])(8,2) = 0;
    (deriv_gp[gp])(9,2) = 4 * xsi4[gp];
        
    // deriv_gp wrt to xsi4 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,3) = 0;
    (deriv_gp[gp])(1,3) = 0;
    (deriv_gp[gp])(2,3) = 0;
    (deriv_gp[gp])(3,3) = 4 * xsi4[gp] - 1;
      
    (deriv_gp[gp])(4,3) = 0;
    (deriv_gp[gp])(5,3) = 0;
    (deriv_gp[gp])(6,3) = 0;
    (deriv_gp[gp])(7,3) = 4 * xsi1[gp];  
    (deriv_gp[gp])(8,3) = 4 * xsi2[gp];
    (deriv_gp[gp])(9,3) = 4 * xsi3[gp];
  }
 }
 
/*----------------------------------------------------------------------*
 | constructor for a integrator class              			  volf 09/07|
 | uses shape functions of a quadratic tetrahedra using so-called       |
 | "natural coordinates" as described by Carlos A. Felippa in Adv. FEM  |
 | Aerospace Engineering Sciences - University of Colorado at Boulder   |
 *----------------------------------------------------------------------*/  
DRT::ELEMENTS::Integrator_tet10_14point::Integrator_tet10_14point(void)
{
  const int number_gp =14;
  // forward initialization of necessary attributes
  num_gp = number_gp;
  num_nodes = NUMNOD_SOTET10 ;
  num_coords = NUMCOORD_SOTET10;
  shapefct_gp.resize(number_gp);
  deriv_gp.resize(number_gp);
  weights.Size(number_gp);

/*	14-point  Quadrature rule from Carlos A. Felippa: Adv. FEM  §17.4 
    If [rule==14,   (* g1,g2 +roots of P(g)=0, P=9+96*g-
      1712*g^2-30464*g^3-127232*g^4+86016*g^5+1060864*g^6 *)
    g1=0.09273525031089122640232391373703060;
	g2=0.31088591926330060979734573376345783;
	g3=0.45449629587435035050811947372066056;
	If [!numer,{g1,g2,g3}=Rationalize[{g1,g2,g3},eps]];
	w1=(-1+6*g2*(2+g2*(-7+8*g2))+14*g3-60*g2*(3+4*g2*
    	(-3+4*g2))*g3+4*(-7+30*g2*(3+4*g2*(-3+4*g2)))*g3^2)/
    	(120*(g1-g2)*(g2*(-3+8*g2)+6*g3+8*g2*(-3+4*g2)*g3-4*
	    (3+4*g2*(-3+4*g2))*g3^2+8*g1^2*(1+12*g2*
    	(-1+2*g2)+4*g3-8*g3^2)+g1*(-3-96*g2^2+24*g3*(-1+2*g3)+
    	g2*(44+32*(1-2*g3)*g3))));
	w2=(-1-20*(1+12*g1*(2*g1-1))*w1+20*g3*(2*g3-1)*(4*w1-1))/
   		(20*(1+12*g2*(2*g2-1)+4*g3-8*g3^2));
	If [i<5,      info={{g1,g1,g1,g1},w1};info[[1,i]]=1-3*g1];
	If [i>4&&i<9, info={{g2,g2,g2,g2},w2};info[[1,i-4]]=1-3*g2];
	If [i>8,      info={{g3,g3,g3,g3},1/6-2*(w1+w2)/3};
    	{j,k}=jk6[[i-8]]; info[[1,j]]=info[[1,k]]=1/2-g3] ];
    	      jk6= {{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}} 
*/
  
  const long double g1=0.09273525031089122640232391373703060;
  const long double g2=0.31088591926330060979734573376345783;
  const long double g3=0.45449629587435035050811947372066056;
  const long double gi1=1-3*g1;
  const long double gi2=1-3*g2;
  const long double gi3=1/2-g3; 
  
  const long double 	w1=(-1+6*g2*(2+g2*(-7+8*g2))+14*g3-60*g2*(3+4*g2*\
    	(-3+4*g2))*g3+4*(-7+30*g2*(3+4*g2*(-3+4*g2)))*g3*g3)/\
    	(120*(g1-g2)*(g2*(-3+8*g2)+6*g3+8*g2*(-3+4*g2)*g3-4*\
	    (3+4*g2*(-3+4*g2))*g3*g3+8*g1*g1*(1+12*g2*\
    	(-1+2*g2)+4*g3-8*g3*g3)+g1*(-3-96*g2*g3+24*g3*(-1+2*g3)+\
    	g2*(44+32*(1-2*g3)*g3))));
  const long double w2=(-1-20*(1+12*g1*(2*g1-1))*w1+20*g3*(2*g3-1)*(4*w1-1))/\
   		(20*(1+12*g2*(2*g2-1)+4*g3-8*g3*g3));
  const long double w3= 1/6-2*(w1+w2)/3;
   			
  // (xsi1, xsi2, xsi3 ,xsi4) gp-locations of fully integrated linear 10-node Tet
  const double xsi1[number_gp] = {gi1,  g1,  g1,  g1, gi2,  g2,  g2,  g2, gi3, gi3, gi3,  g3,  g3,  g3};
  const double xsi2[number_gp] = { g1, gi1,  g1,  g1,  g2, gi2,  g2,  g2, gi3,  g3,  g3, gi3, gi3,  g3};
  const double xsi3[number_gp] = { g1,  g1, gi1,  g1,  g2,  g2, gi2,  g2,  g3, gi3,  g3, gi3,  g3, gi3};
  const double xsi4[number_gp] = { g1,  g1,  g1, gi1,  g2,  g2,  g2, gi2,  g3,  g3, gi3,  g3, gi3, gi3};
  const double w[number_gp]    = { w1,  w1,  w1,  w1,  w2,  w2,  w2,  w2,  w3,  w3,  w3,  w3,  w3,  w3};
	
   // fill up nodal f at each gp 
  for (int gp=0; gp<num_gp; gp++) {
      (shapefct_gp[gp]).Size(num_nodes);
      (shapefct_gp[gp])[0] = xsi1[gp] * (2*xsi1[gp] -1);
      (shapefct_gp[gp])[1] = xsi2[gp] * (2*xsi2[gp] -1);
      (shapefct_gp[gp])[2] = xsi3[gp] * (2*xsi3[gp] -1);
      (shapefct_gp[gp])[3] = xsi4[gp] * (2*xsi4[gp] -1);
      (shapefct_gp[gp])[4] = 4 * xsi1[gp] * xsi2[gp];
      (shapefct_gp[gp])[5] = 4 * xsi2[gp] * xsi3[gp];
      (shapefct_gp[gp])[6] = 4 * xsi3[gp] * xsi1[gp];
      (shapefct_gp[gp])[7] = 4 * xsi1[gp] * xsi4[gp];  
      (shapefct_gp[gp])[8] = 4 * xsi2[gp] * xsi4[gp];
      (shapefct_gp[gp])[9] = 4 * xsi3[gp] * xsi4[gp];
      weights[gp] = w[gp]; 	// just for clarity how to get weight factors    
  }
 
  // fill up df xsi1, xsi2, xsi3, xsi4 directions (NUMDIM) at each gp
  for (int gp=0; gp<num_gp; gp++) {
  	(deriv_gp[gp]).Shape(num_nodes,num_coords);
  	// deriv_gp wrt to xsi1 "(0,..)" for each node(0..9) at each gp [i]
	(deriv_gp[gp])(0,0) = 4 * xsi1[gp]-1;
	(deriv_gp[gp])(1,0) = 0;
    (deriv_gp[gp])(2,0) = 0;
    (deriv_gp[gp])(3,0) = 0;
    
    (deriv_gp[gp])(4,0) = 4 * xsi2[gp];
    (deriv_gp[gp])(5,0) = 0;
    (deriv_gp[gp])(6,0) = 4 * xsi3[gp];
    (deriv_gp[gp])(7,0) = 4 * xsi4[gp];  
    (deriv_gp[gp])(8,0) = 0;
    (deriv_gp[gp])(9,0) = 0;

    // deriv_gp wrt to xsi2 "(1,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,1) = 0;
    (deriv_gp[gp])(1,1) = 4 * xsi2[gp] - 1;
    (deriv_gp[gp])(2,1) = 0;
    (deriv_gp[gp])(3,1) = 0;
     	
    (deriv_gp[gp])(4,1) = 4 * xsi1[gp];
    (deriv_gp[gp])(5,1) = 4 * xsi3[gp];
    (deriv_gp[gp])(6,1) = 0;
    (deriv_gp[gp])(7,1) = 0;  
    (deriv_gp[gp])(8,1) = 4 * xsi4[gp];
    (deriv_gp[gp])(9,1) = 0;

    // deriv_gp wrt to xsi3 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,2) = 0;
    (deriv_gp[gp])(1,2) = 0;
    (deriv_gp[gp])(2,2) = 4 * xsi3[gp] - 1;
    (deriv_gp[gp])(3,2) = 0;
      
    (deriv_gp[gp])(4,2) = 0;
    (deriv_gp[gp])(5,2) = 4 * xsi2[gp];
    (deriv_gp[gp])(6,2) = 4 * xsi1[gp];
    (deriv_gp[gp])(7,2) = 0;  
    (deriv_gp[gp])(8,2) = 0;
    (deriv_gp[gp])(9,2) = 4 * xsi4[gp];
        
    // deriv_gp wrt to xsi4 "(2,..)" for each node(0..9) at each gp [gp]
    (deriv_gp[gp])(0,3) = 0;
    (deriv_gp[gp])(1,3) = 0;
    (deriv_gp[gp])(2,3) = 0;
    (deriv_gp[gp])(3,3) = 4 * xsi4[gp] - 1;
      
    (deriv_gp[gp])(4,3) = 0;
    (deriv_gp[gp])(5,3) = 0;
    (deriv_gp[gp])(6,3) = 0;
    (deriv_gp[gp])(7,3) = 4 * xsi1[gp];  
    (deriv_gp[gp])(8,3) = 4 * xsi2[gp];
    (deriv_gp[gp])(9,3) = 4 * xsi3[gp];
    
  }
  
}

/*----------------------------------------------------------------------*
 * Get shape functions for a tet4 face					       vlf 08/07*
 * ---------------------------------------------------------------------*/
DRT::ELEMENTS::Integrator_tri3_1point::Integrator_tri3_1point(void)
{
  const int number_gp = 1;
  // forward initialization of necessary attributes
  num_gp = number_gp;
  num_nodes = NUMNOD_SOTET4 ;
  num_coords = NUMCOORD_SOTET4;
  shapefct_gp.resize(number_gp);
  deriv_gp.resize(number_gp);
  weights.Size(number_gp);

 //Quadrature rule from Carlos A. Felippa: Adv. FEM  §17 
  const double gploc_alpha    = (double)1/3;    // gp sampling point value for liner. fct
  const double w			  = (double)1;

  const double ksi1[NUMGPT_SOTET4_FACE] = {gploc_alpha };
  const double ksi2[NUMGPT_SOTET4_FACE] = {gploc_alpha };
  const double ksi3[NUMGPT_SOTET4_FACE] = {gploc_alpha };
  
  for (int gp=0; gp<NUMGPT_SOTET4_FACE; gp++) {
  	  (shapefct_gp[gp]).Size(num_nodes);
      shapefct_gp[gp](0) = ksi1[gp];
      shapefct_gp[gp](1) = ksi2[gp];
      shapefct_gp[gp](2) = ksi3[gp];
      weights[gp] = w; // just for clarity how to get weight factors
   }
}

/*----------------------------------------------------------------------*
 |                                                         vlf 06/07    |
 | recursive computation of determinant of a  matrix using Sarrus rule  |
 |                              										|
 *----------------------------------------------------------------------*/
long double det_volf(Epetra_SerialDenseMatrix& in_matrix)
{
	//Epetra_SerialDenseMatrix temp_matrix(in_matrix.N()-1,in_matrix.N()-1);
	
	if (in_matrix.N()==1)
	{
		return (long double) in_matrix(0,0);	
	}	
	else if (in_matrix.N()==2)
	{
		return 	((long double)(in_matrix(0,0)*(long double)in_matrix(1,1))-\
					((long double)in_matrix(0,1)*(long double)in_matrix(1,0)));
	}
	else if (in_matrix.N()>2)
	{
		long double out_det=0;
		int sign=1;
		for (int i_col=0;i_col < in_matrix.N();i_col++)
		{
			Epetra_SerialDenseMatrix temp_matrix(in_matrix.N()-1,in_matrix.N()-1);
			
			for (int c_col=0;c_col < i_col;c_col++)
			{				
				for(int row=1;row<in_matrix.N();row++)
				temp_matrix(row-1,c_col)=in_matrix(row,c_col);							
			}
			for (int c_col=i_col+1;c_col <  in_matrix.N();c_col++)
			{
			     for(int row=1;row<in_matrix.N();row++)
		              temp_matrix(row-1,c_col-1)=in_matrix(row,c_col);	
			}
 
			out_det=out_det+((long double)sign *\
					(long double)in_matrix(0,i_col) * (long double)det_volf(temp_matrix));
			sign*=-1;
		}
		return out_det;
	}	
	else return 0;
}

#endif

#if 0
/*----------------------------------------------------------------------*
 |                                                         vlf 06/07    |
 | recursive computation of determinant of a  matrix using Sarrus rule  |
 |                              										|
 *----------------------------------------------------------------------*/
double det_volf(Epetra_SerialDenseMatrix& in_matrix)
{
	//Epetra_SerialDenseMatrix temp_matrix(in_matrix.N()-1,in_matrix.N()-1);
	
	if (in_matrix.N()==1)
	{
		return in_matrix(0,0);	
	}	
	else if (in_matrix.N()==2)
	{
		return 	((in_matrix(0,0)*in_matrix(1,1))-(in_matrix(0,1)*in_matrix(1,0)));
	}
	else if (in_matrix.N()>2)
	{
		double out_det=0;
		int sign=1;
		for (int i_col=0;i_col < in_matrix.N();i_col++)
		{
			
			Epetra_SerialDenseMatrix temp_matrix(in_matrix.N()-1,in_matrix.N()-1);
			for (int c_col=0;c_col < i_col;c_col++)
			{				
				for(int row=1;row<in_matrix.N();row++)
				temp_matrix(row-1,c_col)=in_matrix(row,c_col);							
			}
			for (int c_col=i_col+1;c_col <  in_matrix.N();c_col++)
			{
			for(int row=1;row<in_matrix.N();row++)
		temp_matrix(row-1,c_col-1)=in_matrix(row,c_col);	
			}
 
			out_det=out_det+(sign* in_matrix(0,i_col)*det_volf(temp_matrix));
			sign*=-1;
		}
		return out_det;
	}	
	else return 0;
}
#endif


#if 0
/*----------------------------------------------------------------------*
 |                                                         vlf 06/07    |
 | recursive computation of determinant of a  matrix using Sarrus rule  |
 |                              										|
 *----------------------------------------------------------------------*/
long double LINALG::SerialDenseMatrix::Det_long()
{
 if (N()==1)
	 {
		return (*this)(0,0);	
	}	
	else if (N()==2)
	{
		return 	((this->(0,0)*this->(1,1))-(this->(0,1)*this->(1,0)));
	}
	else if (N()>2)
	{
		long double (out_det)=0;
		int sign=1;
		for (int i_col=0;i_col < N();i_col++)
		{
			
			SerialDenseMatrix temp_matrix(N()-1,N()-1);
			for (int c_col=0;c_col < i_col;c_col++)
			{				
				for(int row=1;row<N();row++)
				   temp_matrix(row-1,c_col)=in_matrix(row,c_col);							
			}
			for (int c_col=i_col+1;c_col < N();c_col++)
			{
			   for(int row=1;row<N();row++)
			     temp_matrix(row-1,c_col-1)=in_matrix(row,c_col);	
			} 
			out_det = out_det + (long double ( sign)* long double (in_matrix(0,i_col))\
					* long double (temp_matrix.det));
			sign*=-1;
		}
		return out_det;
	}	
	else return 0;
}

#endif
