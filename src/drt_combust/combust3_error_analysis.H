/*----------------------------------------------------------------------*/
/*!
\file combust3_error_analysis.H

\brief comparison to analytical solution for error estimation
       for premixed combustion and two-phase flow problems / 3D combustion element

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef COMBUST3_ERROR_ESTIMATION_H
#define COMBUST3_ERROR_ESTIMATION_H

#include <Teuchos_TimeMonitor.hpp>

#include "combust3_sysmat.H"
#include "combust_defines.H"
#include "combust3_local_assembler.H"
#include "combust3_utils.H"
#include "combust3_interpolation.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_f3/xfluid3_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../drt_xfem/enrichment.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"


using namespace XFEM::PHYSICS;

namespace COMBUST
{


/*!
  Transforms Gaussian point of cell from 3D-Eta-coordinates in 3D-XYZ-coordinates
 */
template <DRT::Element::DiscretizationType CELL_DISTYPE,
          class M1, class M2>
void GetGPinXYZ(
    GEO::DomainIntCell          cell,
    M1&                         pos_eta_domain,
    M2&                         posGPinXYZ)
{
  const int nsd = 3;

  // number of nodes of this cell
  const size_t numnode3Dcell = DRT::UTILS::DisTypeToNumNodePerEle<CELL_DISTYPE>::numNodePerElement;

  // cel 3D-shape functions
  LINALG::Matrix<numnode3Dcell,1> funct3Dcell;
  DRT::UTILS::shape_function_3D(funct3Dcell,pos_eta_domain(0),pos_eta_domain(1),pos_eta_domain(2),CELL_DISTYPE);

  LINALG::SerialDenseMatrix xyze_cell_tmp = cell.CellNodalPosXYZ();
  const LINALG::Matrix<nsd,numnode3Dcell> xyze_cell(xyze_cell_tmp);

  posGPinXYZ.Multiply(xyze_cell,funct3Dcell);
  return;
}



/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: couette_20x20x1_stat.dat
 */
void AnalyticalSol_Couette_20x20x1_plain(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: x+0.025
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 20.0*y+10.0 |
   * u_x(x,y,z)= | 0             |(linear in y, konst in x,z)
   *             | 0             |
   * du_x/dy = 20.0 other derivatives are zero
   *
   * p(x,y,z)= 	7.0 for x<-0.025
   * 				4.0 for x>-0.025
   * gradp = 0.0
   */

  double exact_levelset_gp = posXYZgp(0)+0.025;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


  // set analytical velocity
  analyt_vel(0) = 20.0*posXYZgp(1)+10.0;
  analyt_vel(1) = 0.0;
  analyt_vel(2) = 0.0;

  // set analytical vel_derivatives
  analyt_vel_deriv.Clear();
  analyt_vel_deriv(0,1) = 20.0;

  // set analytical pressure
  if(exact_levelset_gp<0) analyt_pres(0,0) = 7.0;
  else if(exact_levelset_gp>0) analyt_pres(0,0) = 4.0;
  else dserror("domain-Gaussian point lies on Interface");

  // set analytical pressure derivatives
  analyt_pres_deriv.Clear();

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: couette_20x20x1_stat.dat
 */
void AnalyticalSol_Couette_20x20x1_plain_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: x+0.025
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 20.0*y+10.0 |
   * u_x(x,y,z)= | 0             |(linear in y, konst in x,z)
   *             | 0             |
   * du_x/dy = 20.0 other derivatives are zero
   *
   * p(x,y,z)= 	7.0 for x<-0.025
   * 				4.0 for x>-0.025
   * gradp = 0.0
   */


  // set analytical velocity
  analyt_vel_plus(0) = 20.0*posXYZgp(1)+10.0;
  analyt_vel_plus(1) = 0.0;
  analyt_vel_plus(2) = 0.0;

  analyt_vel_minus(0) = 20.0*posXYZgp(1)+10.0;
  analyt_vel_minus(1) = 0.0;
  analyt_vel_minus(2) = 0.0;

  // set analytical vel_derivatives
  analyt_vel_deriv_plus.Clear();
  analyt_vel_deriv_plus(0,1) = 20.0;

  analyt_vel_deriv_minus.Clear();
  analyt_vel_deriv_minus(0,1) = 20.0;


  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: shear.dat
 */
void AnalyticalSol_shear(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{

  // evaluate the analytical solution of the flow problem at posXYZgp
  /*
   *    * !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * FLAGS:
   * #define COMBUST_TESTCOUETTEFLOW
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   *
   * 	Levelset: x+y+0.05
   *
   * [|u|]:=5, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 20.0/0.6*y+10.0 |
   * u_x(x,y,z)= | 0               |(linear in y, konst in x,z) für Omega^-
   *             | 0               |
   *
   *             | 20.0/0.6*y+15.0 |
   * u_x(x,y,z)= | 0               |(linear in y, konst in x,z) für Omega^+
   *             | 0               |
   * du_x/dy = 20.0/0.6 other derivatives are zero
   *
   * p(x,y,z)= 	5.0 for x+y<-0.05
   * 				2.0 for x+y>-0.05
   * gradp = 0.0
   *
   *
   */

  double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


  // set analytical velocity
  if(exact_levelset_gp < 0){
    analyt_vel(0) = 20.0/0.6*posXYZgp(1)+10.0;
    analyt_vel(1) = 0.0;
    analyt_vel(2) = 0.0;
  }
  else if(exact_levelset_gp > 0){
    analyt_vel(0) = 20.0/0.6*posXYZgp(1)+15.0;
    analyt_vel(1) = 0.0;
    analyt_vel(2) = 0.0;
  }

  // set analytical vel_derivatives
  analyt_vel_deriv.Clear();
  analyt_vel_deriv(0,1) = 20.0/0.6;

  // set analytical pressure
  if(exact_levelset_gp < 0) analyt_pres(0,0) = 5.0;
  else if(exact_levelset_gp > 0) analyt_pres(0,0) = 2.0;
  else dserror("domain-Gaussian point lies on Interface");

  // set analytical pressure derivatives
  analyt_pres_deriv.Clear();

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: shear.dat
 */
void AnalyticalSol_shear_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: x+y+0.05
   *
   * [|u|]:=5, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 20.0/0.6*y+10.0 |
   * u_x(x,y,z)= | 0               |(linear in y, konst in x,z) für Omega^-
   *             | 0               |
   *
   *             | 20.0/0.6*y+15.0 |
   * u_x(x,y,z)= | 0               |(linear in y, konst in x,z) für Omega^+
   *             | 0               |
   * du_x/dy = 20.0/0.6 other derivatives are zero
   *
   * p(x,y,z)= 	5.0 for x+y<-0.05
   * 				2.0 for x+y>-0.05
   * gradp = 0.0
   *
   *
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * FLAGS:
   * #define COMBUST_TESTCOUETTEFLOW
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */

  //	double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  //	// check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  //	// matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  //	if((matid == 3 && exact_levelset_gp > 0) ||
  //	   (matid == 4 && exact_levelset_gp < 0)); // do nothing
  //	else cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


  // set analytical velocity

  analyt_vel_minus(0) = 20.0/0.6*posXYZgp(1)+10.0;
  analyt_vel_minus(1) = 0.0;
  analyt_vel_minus(2) = 0.0;

  analyt_vel_plus(0) = 20.0/0.6*posXYZgp(1)+15.0;
  analyt_vel_plus(1) = 0.0;
  analyt_vel_plus(2) = 0.0;


  // set analytical vel_derivatives
  analyt_vel_deriv_plus.Clear();
  analyt_vel_deriv_plus(0,1) = 20.0/0.6;

  analyt_vel_deriv_minus.Clear();
  analyt_vel_deriv_minus(0,1) = 20.0/0.6;


  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: straight_bodyforce.dat
 */
void AnalyticalSol_straight_bodyforce_inclined(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: x+y+0.05
   *
   * bodyforce(0) = -10.0 =>   dp/dx= -10.0
   *
   * [|u|]:=5, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 2.0 |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^-
   *             | 0               |
   *
   *             | 7.0 |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^+
   *             | 0               |
   * derivatives are zero
   *
   * p(x,y,z)= 	5.0 for x+y<-0.05
   * 				2.0 for x+y>-0.05
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


  // set analytical velocity
  if(exact_levelset_gp < 0){
    analyt_vel(0) = 2.0;
    analyt_vel(1) = 0.0;
    analyt_vel(2) = 0.0;
  }
  else if(exact_levelset_gp > 0){
    analyt_vel(0) = 7.0;
    analyt_vel(1) = 0.0;
    analyt_vel(2) = 0.0;
  }

  // set analytical vel_derivatives
  analyt_vel_deriv.Clear();


  // set analytical pressure
  if(exact_levelset_gp < 0) analyt_pres(0,0) = 17.0-10.0*(posXYZgp(0)+0.5);
  else if(exact_levelset_gp > 0) analyt_pres(0,0) = 17.0-10.0*(posXYZgp(0)+0.5) - 3.0;
  else dserror("domain-Gaussian point lies on Interface");

  // set analytical pressure derivatives
  analyt_pres_deriv.Clear();
  analyt_pres_deriv(0,0)= -10.0;
  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: shear.dat
 */
void AnalyticalSol_straight_bodyforce_inclined_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: x+y+0.05
   *
   * bodyforce(0) = 10.0 =>   dp/dx= -10.0
   *
   * [|u|]:=5, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 2.0 |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^-
   *             | 0               |
   *
   *             | 7.0 |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^+
   *             | 0               |
   * derivatives are zero
   *
   * p(x,y,z)= 	5.0 for x+y<-0.05
   * 				2.0 for x+y>-0.05
   * gradp = 0.0
   *
   */

  //	double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  //	// check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  //	// matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  //	if((matid == 3 && exact_levelset_gp > 0) ||
  //	   (matid == 4 && exact_levelset_gp < 0)); // do nothing
  //	else cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


  // set analytical velocity

  analyt_vel_minus(0) = 2.0;
  analyt_vel_minus(1) = 0.0;
  analyt_vel_minus(2) = 0.0;

  analyt_vel_plus(0) = 7.0;
  analyt_vel_plus(1) = 0.0;
  analyt_vel_plus(2) = 0.0;


  // set analytical vel_derivatives
  analyt_vel_deriv_plus.Clear();

  analyt_vel_deriv_minus.Clear();


  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   * 				2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = 0.5;
  double curvature = 1.0/rad_1;
  const double surftens = 70.0*1e-003;



  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  if(exact_levelset_gp>0) analyt_pres=0.0;
  else if(exact_levelset_gp<0) analyt_pres= 0.0+curvature*surftens;

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   * 				2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnx1(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   * 				2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.16;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = 0.4;
  double curvature = 1.0/rad_1;
  const double surftens = 70.0*1e-003;



  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //	if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04>0) analyt_pres=0.0;
  //	else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if(matid==3) analyt_pres=0.0;
  else if(matid==4) analyt_pres= 0.0+curvature*surftens;



  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnx1_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   * 				2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnxn(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09<0
   * 				2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09<0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.09;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else
  {

  } //cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = 0.3;
  double rad_2 = 0.3;
  double curvature = 1.0/rad_1 + 1.0/rad_2;
  const double surftens = 70.0*1e-003;


  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //	if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04>0) analyt_pres=0.0;
  //	else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if(matid==3) analyt_pres=0.0;
  else if(matid==4) analyt_pres= 0.0+curvature*surftens;
  //		if(matid==3) analyt_pres=4.0;
  //		else if(matid==4) analyt_pres= 4.0+0.0;


  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnxn_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.09<0
   * 				2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.09>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!!!!!!!!!!

  return;
}


/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnxn_instat(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   * 				0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else
  {

  } //cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = 0.232401774;
  double rad_2 = 0.232401774;
  double curvature = 1.0/rad_1 + 1.0/rad_2;
  const double surftens = 70.0*1e-003;


  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //	if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04>0) analyt_pres=0.0;
  //	else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if(matid==3) analyt_pres=0.0;
  else if(matid==4) analyt_pres= 0.0+curvature*surftens;
  //		if(matid==3) analyt_pres=4.0;
  //		else if(matid==4) analyt_pres= 4.0+0.0;


  return;
}

/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnxn_instat_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)(z-0.01)*(z-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   * 				0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!!!!!!!!!!

  return;
}



/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnx1_instat(LINALG::Matrix<3,1>& posXYZgp, int matid,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.01
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   * 				0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.01;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((matid == 3 && exact_levelset_gp > 0) ||
      (matid == 4 && exact_levelset_gp < 0)); // do nothing
  else
  {

  } //cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = sqrt(0.1);
  //double rad_2 = 0.1;
  double curvature = 1.0/rad_1;
  const double surftens = 70.0*1e-003;


  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //	if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04>0) analyt_pres=0.0;
  //	else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if(matid==3) analyt_pres=0.0;
  else if(matid==4) analyt_pres= 0.0+curvature*surftens;
  //		if(matid==3) analyt_pres=4.0;
  //		else if(matid==4) analyt_pres= 4.0+0.0;


  return;
}

/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnx1_instat_Interface(LINALG::Matrix<3,1>& posXYZgp,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*	Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)(z-0.01)*(z-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)= 	2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   * 				0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!!!!!!!!!!

  return;
}


/*!
  Calculate Nitsche errors with domain integrals
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class V2>
void Nitsche_BuildDomainIntegratedErrors(
    ParameterList&                 eleparams,       /// list of parameters with actions to be calculated
    const INPAR::COMBUST::NitscheError&  NitscheErrorType, /// which solution is the reference
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih, ///< information about the interface
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const M1&                      evelnp,       ///< nodal velocity n+1 values for this element
    const V1&                      eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                      ephi,         ///< nodal G-function values for this element
    Teuchos::RCP<const MAT::Material> material,  ///< fluid material
    double&                        ele_meas_plus,
    double&                        ele_meas_minus
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion Nitsche errors - domain");

  //==================================================================================================  
  // get actual values of norms^2 of domain errors,stored in eleparams, then update the values with error part of this element
  // at the end of this function we set the eleparams new
  double eleVelDomErr 			= eleparams.get<double>("L2 integrated velocity domain error");
  double eleGradVelDomErr 		= eleparams.get<double>("L2 integrated grad_velocity domain error");
  double elePresDomErr			= eleparams.get<double>("L2 integrated pressure domain error");
  double eleGradPresDomErr		= eleparams.get<double>("L2 integrated grad_pressure domain error");
  double eleWeightPresDomErr	= eleparams.get<double>("L2 integrated weighted pressure domain error");
  double eleNitscheErr			= eleparams.get<double>("Nitsche integrated error");
  // here are stored the !!!squared norms!!!
  //==================================================================================================

  // function to evaluate analytical solution at Gaussian point -> return velx,vely,velz,pres at Gaussian point
  // pointer to function
  void (*analytical)(	LINALG::Matrix<3,1>& posXYZgp,
      int matid,
      LINALG::Matrix<3,1>& analyt_vel,
      LINALG::Matrix<3,3>& analyt_vel_deriv,
      LINALG::Matrix<1,1>& analyt_pres,
      LINALG::Matrix<3,1>& analyt_pres_deriv);
  analytical=NULL;

  // set pointer to function which evaluates the analytical solution

  if     (NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnx1)
    analytical = &AnalyticalSol_static_bubble_nxnx1;	// static bubble with n x n x 1 Elements
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnxn)
    analytical = &AnalyticalSol_static_bubble_nxnxn;  // static bubble with n x n x n Elements
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_couette_20x20x1)
    analytical = &AnalyticalSol_Couette_20x20x1_plain;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_shear)
    analytical = &AnalyticalSol_shear;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_straight_bodyforce)
    analytical = &AnalyticalSol_straight_bodyforce_inclined;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_ellipsoid_bubble_3D)
    analytical = &AnalyticalSol_ellipsoid_bubble_nxnxn_instat;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_ellipsoid_bubble_2D)
    analytical = &AnalyticalSol_ellipsoid_bubble_nxnx1_instat;
  else dserror("there is no analytical solution for 'flowproblem' stored in eleparams");

  //==================================================================================================
  //update the errors with contribution of this element ele

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // get node coordinates of this element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif
  // get material list for this element
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());


  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele));

  //============================= get parameters for norm calculation ==========

  // get both viscosities for weighted pressure norm -> mu_max is needed
  Teuchos::RCP<const MAT::Material> matptr_plus = matlist->MaterialById(3);
  Teuchos::RCP<const MAT::Material> matptr_minus = matlist->MaterialById(4);

  const MAT::NewtonianFluid* mat_plus = static_cast<const MAT::NewtonianFluid*>(matptr_plus.get());
  const MAT::NewtonianFluid* mat_minus = static_cast<const MAT::NewtonianFluid*>(matptr_minus.get());

  const double dynvisc_plus = mat_plus->Viscosity();
  const double dynvisc_minus = mat_minus->Viscosity();

  const double densnp_plus = mat_plus->Density();
  const double densnp_minus = mat_minus->Density();


  const double dynvisc_max = max(dynvisc_plus,dynvisc_minus);

  //=================================================================================
  //----------------------------------------------------------------------------------------------
  // loop over domain integration cells
  //----------------------------------------------------------------------------------------------
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {

    double densnp = 0.0;
    double dynvisc = 0.0;

    //------------------------------------------------------------------------------------------
    // get material parameters for this integration cell
    //------------------------------------------------------------------------------------------
    int matid = 777;
    // check on which side of the interface the cell is located
    if(cell->getDomainPlus()) // cell belongs to burnt domain
    {
      // burnt material
      matid = 3;
      densnp = densnp_plus;
      dynvisc = dynvisc_plus;
    }
    else // cell belongs to unburnt domain
    {
      // unburnt material
      matid = 4;
      densnp = densnp_minus;
      dynvisc = dynvisc_minus;
    }


    //------------------------------------------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //------------------------------------------------------------------------------------------
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);

    //------------------------------------------------------------------------------------------
    // get Gaussian points for this integration cell
    //------------------------------------------------------------------------------------------
    // TODO: Check Gaußrule
    // we only need to integrate u_h and p_h exact,
    // if the exact solution is smoother we wont get a better
    // order of convergence because u_h and p_h approximates u and p only with this order of convergence
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Intersected(), cell->Shape());
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);


    //------------------------------------------------------------------------------------------
    // integration loop over Gaussian points
    //------------------------------------------------------------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //--------------------------------------------------------------------------------------
      // transform coordinates of this Gaussian point
      //--------------------------------------------------------------------------------------
      // coordinates of the current integration point in cell coordinates \eta^domain
      const LINALG::Matrix<nsd,1> pos_eta_domain(intpoints.qxg[iquad]);

      // coordinates of the current integration point in element coordinates \xi
      static LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
      const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
#ifdef DEBUG
      if (detcell < 0.0)
      {
        cout << "detcell :  " << detcell << endl;
        dserror("negative detcell!");
      }
#endif
      //--------------------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //--------------------------------------------------------------------------------------
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      //--------------------------------------------------------------------------------------
      // procedures involving Jacobian matrix
      //--------------------------------------------------------------------------------------
      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
      xjm.MultiplyNT(deriv,xyze);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det*detcell;
      //cout << "domain spatial factor " << fac << endl;

      if (det < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
      }

      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Invert(xjm);

      //--------------------------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //--------------------------------------------------------------------------------------
      static LINALG::Matrix<3,numnode> derxy;
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);

      // flag for higher order elements
      const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();
      // compute second global derivative
      static LINALG::Matrix<6,numnode> derxy2;
      if (higher_order_ele)
      {
        static LINALG::Matrix<6,numnode> deriv2;
        DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
        DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }


      //--------------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        //cout << "this is a fully or partially enriched element" << endl;
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSize,1> enr_funct;
        enr_funct.Clear();
        static LINALG::Matrix<3,shpVecSize> enr_derxy;
        enr_derxy.Clear();
        static LINALG::Matrix<6,shpVecSize> enr_derxy2;
        enr_derxy2.Clear();

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            Velx,
            funct,
            derxy,
            derxy2,
            enr_funct,
            enr_derxy,
            enr_derxy2);

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          shp.d0(iparam) = enr_funct(iparam);
          shp.dx(iparam) = enr_derxy(0,iparam);
          shp.dy(iparam) = enr_derxy(1,iparam);
          shp.dz(iparam) = enr_derxy(2,iparam);
          shp.dxdx(iparam) = enr_derxy2(0,iparam);
          shp.dxdy(iparam) = enr_derxy2(3,iparam);
          shp.dxdz(iparam) = enr_derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = enr_derxy2(1,iparam);
          shp.dydz(iparam) = enr_derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = enr_derxy2(2,iparam);
        }
      }
      else // not xfem_assembly i.e. standard assembly
      {
        // fill approximation functions for standard FEM
        // remark: numparamvelx == numnode, for standard FEM
        for (size_t iparam = 0; iparam < numnode; ++iparam)
        {
          shp.d0(iparam) = funct(iparam);
          shp.dx(iparam) = derxy(0,iparam);
          shp.dy(iparam) = derxy(1,iparam);
          shp.dz(iparam) = derxy(2,iparam);
          shp.dxdx(iparam) = derxy2(0,iparam);
          shp.dxdy(iparam) = derxy2(3,iparam);
          shp.dxdz(iparam) = derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = derxy2(1,iparam);
          shp.dydz(iparam) = derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = derxy2(2,iparam);
        }
      }

      if(cell->getDomainPlus()){ele_meas_plus += fac;}
      else {ele_meas_minus += fac;}

      //--------------------------------------------------------------------------------------
      // get velocities of solution at integration point
      //--------------------------------------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);

      //--------------------------------------------------------------------------------------
      // get velocity derivatives of solution at integration point
      //--------------------------------------------------------------------------------------
      static LINALG::Matrix<3,nsd> vderxy;
      vderxy.Clear();

      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          vderxy(isd,0) += evelnp(isd,iparam) * shp.dx(iparam);
          vderxy(isd,1) += evelnp(isd,iparam) * shp.dy(iparam);
          vderxy(isd,2) += evelnp(isd,iparam) * shp.dz(iparam);
        }
      }


      //--------------------------------------------------------------------------------------
      // get pressure at integration point
      //--------------------------------------------------------------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shp.d0(iparam);

      //--------------------------------------------------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shp.dx(iparam);
        gradp(1) += eprenp(iparam) * shp.dy(iparam);
        gradp(2) += eprenp(iparam) * shp.dz(iparam);
      }


      //=========================== ANALYTICAL SOLUTION ===============
      // get the analytical solution at Gaussian point
      // variables for analytical solution
      static LINALG::Matrix<nsd,1> analyt_vel;
      static LINALG::Matrix<nsd,nsd> analyt_vel_derxy;
      static LINALG::Matrix<1,1> analyt_pres;
      static LINALG::Matrix<nsd,1> analyt_pres_derxy;
      analyt_vel.Clear();
      analyt_vel_derxy.Clear();
      analyt_pres.Clear();
      analyt_pres_derxy.Clear();

      // evaluate the analytical solution
      static LINALG::Matrix<nsd,1> posGPinXYZ;
      posGPinXYZ.Clear();


      // get posGPinXYZ and the domain (according to phi_h) in which the Gaussian point lies

      switch (cell->Shape())
      {
      case DRT::Element::tet4:
        COMBUST::GetGPinXYZ<DRT::Element::tet4>(
            *cell, pos_eta_domain, posGPinXYZ);
        break;
      case DRT::Element::hex8:
        COMBUST::GetGPinXYZ<DRT::Element::hex8>(
            *cell, pos_eta_domain, posGPinXYZ);
        break;
      default:
        dserror("cell distype not implemented yet");
      }

      analytical(posGPinXYZ, matid, analyt_vel, analyt_vel_derxy, analyt_pres, analyt_pres_derxy);

      //============================== NITSCHE ERRORS =================================
      // calculate difference between analytical solution and converged solution in function and derivatives

      static LINALG::Matrix<nsd,1> err_vel;			// v-v_h
      static LINALG::Matrix<nsd,nsd> err_vel_derxy;	// grad(v-v_h)
      static LINALG::Matrix<1,1> err_pres;			// p-p_h
      static LINALG::Matrix<nsd,1> err_pres_derxy;	// grad(p-p_h)
      err_vel.Clear();
      err_vel_derxy.Clear();
      err_pres.Clear();
      err_pres_derxy.Clear();

      // loop over Velx,Vely,Velz
      // get difference between exact and approximative solution
      for (size_t i = 0; i< nsd; ++i)
      {
        err_vel(i) = analyt_vel(i)-gpvelnp(i);
        // 0=d/dx,1=d/dy,2=d/dz
        err_vel_derxy(i,0) = analyt_vel_derxy(i,0) - vderxy(i,0); // d/dx
        err_vel_derxy(i,1) = analyt_vel_derxy(i,1) - vderxy(i,1); // d/dy
        err_vel_derxy(i,2) = analyt_vel_derxy(i,2) - vderxy(i,2); // d/dz
      }

      err_pres = analyt_pres(0,0) - pres;
      err_pres_derxy(0,0) = analyt_pres_derxy(0,0) - gradp(0); // d/dx
      err_pres_derxy(1,0) = analyt_pres_derxy(1,0) - gradp(1); // d/dy
      err_pres_derxy(2,0) = analyt_pres_derxy(2,0) - gradp(2); // d/dz


      // update element errors
      for (size_t i = 0; i< nsd; ++i)
      {
        eleVelDomErr 	+= err_vel(i)*err_vel(i)*fac;
        eleNitscheErr 	+= err_vel(i)*err_vel(i)*fac;

        eleGradVelDomErr 	+= dynvisc * err_vel_derxy(0,i) * err_vel_derxy(0,i) * fac;
        eleGradVelDomErr 	+= dynvisc * err_vel_derxy(1,i) * err_vel_derxy(1,i) * fac;
        eleGradVelDomErr 	+= dynvisc * err_vel_derxy(2,i) * err_vel_derxy(2,i) * fac;
        eleNitscheErr 		+= dynvisc * err_vel_derxy(0,i) * err_vel_derxy(0,i) * fac;
        eleNitscheErr 		+= dynvisc * err_vel_derxy(1,i) * err_vel_derxy(1,i) * fac;
        eleNitscheErr 		+= dynvisc * err_vel_derxy(2,i) * err_vel_derxy(2,i) * fac;

        eleGradPresDomErr   += err_pres_derxy(i,0) * err_pres_derxy(i,0) * fac;

      }

      elePresDomErr 			+= err_pres(0,0) * err_pres(0,0) * fac;

      eleWeightPresDomErr 	+= 1.0/dynvisc_max * err_pres(0,0) * err_pres(0,0) * fac;
      eleNitscheErr			+= 1.0/dynvisc_max * err_pres(0,0) * err_pres(0,0) * fac;
    } // end loop over Gauß points
  } // end loop over integration cells

  //==================================================================================================  
  //  cout << eleGradPresDomErr << endl;

  //==================================================================================================
  // set the domain errors with the updated values
  eleparams.set<double>("L2 integrated velocity domain error", eleVelDomErr);
  eleparams.set<double>("L2 integrated grad_velocity domain error", eleGradVelDomErr);
  eleparams.set<double>("L2 integrated pressure domain error", elePresDomErr);
  eleparams.set<double>("L2 integrated grad_pressure domain error", eleGradPresDomErr);
  eleparams.set<double>("L2 integrated weighted pressure domain error", eleWeightPresDomErr);
  eleparams.set<double>("Nitsche integrated error", eleNitscheErr);
  //==================================================================================================

  return;
}


/*!
  Calculate boundary Nitsche errors for current cell
 */
template <DRT::Element::DiscretizationType DISTYPE,
          DRT::Element::DiscretizationType CELLDISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          class M1, class M2, class M3, class V1, class V2>
void Nitsche_BuildBoundaryIntegratedErrors_Cell(
    void (*analytical_Interface)( LINALG::Matrix<3,1>& posXYZgp,
                                  LINALG::Matrix<3,1>& analyt_vel_plus,
                                  LINALG::Matrix<3,1>& analyt_vel_minus,
                                  LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
                                  LINALG::Matrix<3,3>& analyt_vel_deriv_minus),
    double&                           eleViscInterfErr,      /// "H-1/2 integrated viscosity interface error"
    double&                           eleVelJumpInterfErr,   /// "H1/2 integrated velocity jump interface error"
    double&                           eleNitscheErr,         /// "Nitsche integrated error"
    const double&                     hk_eleDiam,
    GEO::BoundaryIntCell              cell,
    const DRT::UTILS::GaussRule2D&    intrule2D,
    M2&                               xyze,
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const bool&                       higher_order_ele,
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M3&                         grad_ephi_smoothed,
    const XFEM::ElementEnrichmentValues&       enrvals_plus,
    const XFEM::ElementEnrichmentValues&       enrvals_minus,
    const double&                              kappa_plus,    /// weights of average in plus domain
    const double&                              kappa_minus,    /// weights of average in minus domain
    const double&                              visc_plus,    /// viscosity in plus domain
    const double&                              visc_minus,   /// viscosity in minus domain
    const size_t&                              numparamvelx,
    const bool                                 smoothed_boundary_integration
)
{
  const size_t nsd =3;

  // number of nodes of element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  // TODO: are 6 Gaussian points enough?
  //const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);
  const DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

  // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
  LINALG::SerialDenseMatrix cellXiDomaintmp = cell.CellNodalPosXiDomain();
  // get coordinates of vertices of boundary integration cell in physical coordinates \x^domain
  LINALG::SerialDenseMatrix cellXYZDomaintmp = cell.CellNodalPosXYZ();
  // transform to fixed size format
  const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);
  const LINALG::Matrix<nsd,numvertices> cellXYZDomain(cellXYZDomaintmp);

  //----------------------------------------------------------------------------------------------
  // integration loop over Gaussian points
  //----------------------------------------------------------------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {

    // new transformation for boundary integrals
    // 1. define a coupled transformation x_3D(xi_3D(eta_2D)): transformation from 2D->3D
    // 2. compute the corresponding Jacobian J_eta2D->x_3D and
    // 3. the corresponding surface integral factor sqrt(det(J_eta2D->x_3D^T * J_eta2D->x_3D))
    // 4. approximate integral with Gauss rule in eta coordinate system
    // 5. evaluate the transformed integrand f(x(xi(eta)))

    const LINALG::Matrix<2,1> gpinEta2D(intpoints.qxg[iquad]);

    // jacobian for coupled transformation
    // get derivatives dxi_3D/deta_2D
    static LINALG::Matrix<2,numvertices> deriv_eta2D;
    DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),CELLDISTYPE);

    // calculate dxi3Ddeta2D
    static LINALG::Matrix<3,2> dXi3Ddeta2D;
    dXi3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < (int)numvertices; k++)
          dXi3Ddeta2D(i,j) += cellXiDomain(i,k)*deriv_eta2D(j,k);

    // transform Gauss point to xi3D space (element parameter space)
    static LINALG::Matrix<nsd,1> gpinXi3D;
    gpinXi3D.Clear();
    // coordinates of this integration point in element coordinates \xi^domain
    GEO::mapEtaBToXiD(cell, gpinEta2D, gpinXi3D);


    static LINALG::Matrix<3,numnode> deriv_xi3D;
    DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DISTYPE);


    // calculate dx3Ddxi3D
    static LINALG::Matrix<3,3> dX3DdXi3D;
    dX3DdXi3D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 3; j++) // derivatives
        for (int k = 0; k < (int)numnode; k++)
          dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);

    // get the coupled Jacobian dx3Ddeta2D
    static LINALG::Matrix<3,2> dx3Ddeta2D;
    dx3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < 3; k++)
          dx3Ddeta2D(i,j) += dX3DdXi3D(i,k) * dXi3Ddeta2D(k,j);

    // get deformation factor
    static LINALG::Matrix<2,2> Jac_tmp; // J^T*J
    Jac_tmp.Clear();
    Jac_tmp.MultiplyTN(dx3Ddeta2D,dx3Ddeta2D);

    if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");

    const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

    const double fac = intpoints.qwgt[iquad]*deform_factor;



    LINALG::Matrix<2,1> posEtaBoundary;
    posEtaBoundary.Clear();
    for (int i= 0; i< 2; i++){
      posEtaBoundary(i,0) = gpinEta2D(i,0);
    }

    LINALG::Matrix<nsd,1> posXiDomain;
    posXiDomain.Clear();
    for (int i= 0; i< 3; i++){
      posXiDomain(i,0) = gpinXi3D(i,0);
    }

    //--------------------------------------------------------------------------------------------
    // compute normal vector (in physical coordinates)
    // remark: for linear boundary integrastion cells this could be done before the loop over all
    //         Gaussian points
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> normal(true);
    normal.Clear();

    // schott Sep 2, 2010
    // call smoothed normals for an accurate evaluation of boundary integration terms
    if(smoothed_boundary_integration == false)
    {
      // here normal points from Omega^+ to Omega^-
      GEO::computeNormalToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal);
      // but we need a normal vector pointing from Omega^- to Omega^+ -> Scale(-1.0) => n = n^- = -n^+
      normal.Scale(-1.0);
    }
    else // smoothed_boundary_integration == true
    {
      // get shape functions at gaussian point
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
      DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);


      normal.Clear();
      LINALG::Matrix<nsd,1> grad_xi;
      grad_xi.Clear();

      // loop over nodes
      for(size_t i = 0; i< numnode; i++)
      {
        LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
        LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
        nodal_grad_xyz_tmp.Clear();
        nodal_grad_xi_tmp.Clear();


        normal(0)+=grad_ephi_smoothed(0,i)*funct_gp(i);
        normal(1)+=grad_ephi_smoothed(1,i)*funct_gp(i);
        normal(2)+=grad_ephi_smoothed(2,i)*funct_gp(i);
      }

      // 4. normalize final gradient
      if (normal.Norm2() == 0.0 ) dserror("normal norm is zero !!! ");
      normal.Scale(1.0/normal.Norm2());
    }



    //--------------------------------------------------------------------------------------------
    // evaluate shape functions and their first derivatives at this Gaussian point
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<numnode,1> funct;
    funct.Clear();
    static LINALG::Matrix<nsd,numnode> deriv;
    deriv.Clear();
    DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
    DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    //--------------------------------------------------------------------------------------------
    // procedures involving Jacobian matrix for domain mapping
    //--------------------------------------------------------------------------------------------
    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
    static LINALG::Matrix<nsd,nsd> xjm;
    xjm.Clear();
    xjm.MultiplyNT(deriv,xyze);
    // determinant for mapping from physical space (X^3D/domain) to element (Xi^3D/domain)
    const double detXtoXi = xjm.Determinant();

    if (detXtoXi < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), detXtoXi);
    }

    // inverse of jacobian
    static LINALG::Matrix<nsd,nsd> xji;
    xji.Clear();
    xji.Invert(xjm);

    //--------------------------------------------------------------------------------------------
    // compute global derivates of shape functions at this Gaussian point
    //--------------------------------------------------------------------------------------------
    // compute first global derivative
    static LINALG::Matrix<3,numnode> derxy;
    derxy.Clear();
    // derxy(i,j) = xji(i,k) * deriv(k,j)
    derxy.Multiply(xji,deriv);

    // compute second global derivative
    static LINALG::Matrix<6,numnode> derxy2;
    if (higher_order_ele) // if (higher_order_ele)
    {
      static LINALG::Matrix<6,numnode> deriv2;
      DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
    }
    else
    {
      derxy2.Clear();
    }



    //--------------------------------------------------------------------------------------------
    // rearrange (enriched) shape functions and derivatives as approximation functions
    //--------------------------------------------------------------------------------------------
    const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
    static LINALG::Matrix<shpVecSize,1>       enrfunct_plus;
    enrfunct_plus.Clear();
    static LINALG::Matrix<shpVecSize,1>       enrfunct_minus;
    enrfunct_minus.Clear();

    static LINALG::Matrix<3,shpVecSize> enrderxy_plus;
    enrderxy_plus.Clear();
    static LINALG::Matrix<3,shpVecSize> enrderxy_minus;
    enrderxy_minus.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_plus;
    enrderxy2_plus.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_minus;
    enrderxy2_minus.Clear();

    // shape functions for nodal parameters (dofs) on plus and minus side
    enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
        enrfunct_plus, enrderxy_plus, enrderxy2_plus);
    enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
        enrfunct_minus, enrderxy_minus, enrderxy2_minus);
    // schott Apr 29, 2010:   Velx okay, all shapefunctions for Velx, Vely, Velz, Pres are the same
    // okay schott Apr 29 2010

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<0,shpVecSize> shp_jump;
    static XFEM::ApproxFunc<1,shpVecSize> shp_mean_visc;


#ifdef TPF_NITSCHE_LPS
    static XFEM::ApproxFunc<0,shpVecSize> shp_presstab;

#endif

    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_jump.d0(iparam) = enrfunct_plus(iparam) - enrfunct_minus(iparam);

      shp_mean_visc.dx(iparam) = kappa_plus*visc_plus*enrderxy_plus(0,iparam) + kappa_minus*visc_minus*enrderxy_minus(0,iparam);
      shp_mean_visc.dy(iparam) = kappa_plus*visc_plus*enrderxy_plus(1,iparam) + kappa_minus*visc_minus*enrderxy_minus(1,iparam);
      shp_mean_visc.dz(iparam) = kappa_plus*visc_plus*enrderxy_plus(2,iparam) + kappa_minus*visc_minus*enrderxy_minus(2,iparam);
    }

    //--------------------------------------------------------------------------------------------
    // compute data at Gaussian point for rhs
    //--------------------------------------------------------------------------------------------

    // velocity jump
    static LINALG::Matrix<nsd,1> vjump(true);
    vjump.Clear();
    vjump = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_jump.d0, numparamvelx);


    // get velocity (np,i) derivatives at integration point
    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    static LINALG::Matrix<nsd,nsd> vderxy_mean_visc(true);
    vderxy_mean_visc.Clear();
    for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      for (size_t isd = 0; isd < nsd; ++isd)
      {
        vderxy_mean_visc(isd,0) += evelnp(isd,iparam) * shp_mean_visc.dx(iparam);
        vderxy_mean_visc(isd,1) += evelnp(isd,iparam) * shp_mean_visc.dy(iparam);
        vderxy_mean_visc(isd,2) += evelnp(isd,iparam) * shp_mean_visc.dz(iparam);
      }



    //=========================== ANALYTICAL SOLUTION ===============
    // get the analytical solution at Gaussian point
    // variables for analytical solution at each side (values and derivatives)
    static LINALG::Matrix<nsd,1> analyt_vel_plus;
    static LINALG::Matrix<nsd,1> analyt_vel_minus;
    static LINALG::Matrix<nsd,nsd> analyt_deriv_plus;
    static LINALG::Matrix<nsd,nsd> analyt_deriv_minus;

    analyt_vel_plus.Clear();
    analyt_vel_minus.Clear();
    analyt_deriv_plus.Clear();
    analyt_deriv_minus.Clear();


    // evaluate the analytical solution
    static LINALG::Matrix<nsd,1> posEtaBoundaryinXYZ;
    posEtaBoundaryinXYZ.Multiply(xyze,funct);

    analytical_Interface(posEtaBoundaryinXYZ, analyt_vel_plus, analyt_vel_minus, analyt_deriv_plus, analyt_deriv_minus);

    //============================== NITSCHE ERRORS =================================

    // get [v]
    static LINALG::Matrix<nsd,1> analyt_vel_jump;
    analyt_vel_jump.Clear();
    for (size_t i = 0; i< nsd; i++)
    {
      analyt_vel_jump(i) = analyt_vel_plus(i)-analyt_vel_minus(i);
    }

    // get [v_h]
    // => [v_h] is stored in vjump

    // get deriv_visc_mean = {2*mu*E(v_h)}
    // vderxy_mean_visc stores {mu*v_x,y} (3x3)
    static LINALG::Matrix<nsd,nsd> deriv_visc_mean;
    deriv_visc_mean.Clear();
    deriv_visc_mean(0,0) = 2.0* (vderxy_mean_visc(0,0));
    deriv_visc_mean(0,1) =      (vderxy_mean_visc(1,0) + vderxy_mean_visc(0,1));
    deriv_visc_mean(0,2) =      (vderxy_mean_visc(2,0) + vderxy_mean_visc(0,2));
    deriv_visc_mean(1,1) = 2.0* (vderxy_mean_visc(1,1));
    deriv_visc_mean(1,2) =      (vderxy_mean_visc(1,2) + vderxy_mean_visc(2,1));
    deriv_visc_mean(2,2) = 2.0* (vderxy_mean_visc(2,2));
    deriv_visc_mean(1,0) = deriv_visc_mean(0,1);
    deriv_visc_mean(2,0) = deriv_visc_mean(0,2);
    deriv_visc_mean(2,1) = deriv_visc_mean(1,2);


    // get {mu * v_x,y} => 3x3
    static LINALG::Matrix<nsd,nsd> analyt_vderxy_mean_visc;
    analyt_vderxy_mean_visc.Clear();
    for(size_t isd = 0; isd < nsd; isd++)
    {
      analyt_vderxy_mean_visc(isd,0) = 0.5* (visc_plus * analyt_deriv_plus(isd,0) + visc_minus * analyt_deriv_minus(isd,0));
      analyt_vderxy_mean_visc(isd,1) = 0.5* (visc_plus * analyt_deriv_plus(isd,1) + visc_minus * analyt_deriv_minus(isd,1));
      analyt_vderxy_mean_visc(isd,2) = 0.5* (visc_plus * analyt_deriv_plus(isd,2) + visc_minus * analyt_deriv_minus(isd,2));
    }


    // get analyt_deriv_visc_mean = {2*mu*E(v)}
    static LINALG::Matrix<nsd,nsd> analyt_deriv_visc_mean;
    analyt_deriv_visc_mean.Clear();
    analyt_deriv_visc_mean(0,0) = 2.0* (analyt_vderxy_mean_visc(0,0));
    analyt_deriv_visc_mean(0,1) =      (analyt_vderxy_mean_visc(1,0) + analyt_vderxy_mean_visc(0,1));
    analyt_deriv_visc_mean(0,2) =      (analyt_vderxy_mean_visc(2,0) + analyt_vderxy_mean_visc(0,2));
    analyt_deriv_visc_mean(1,1) = 2.0* (analyt_vderxy_mean_visc(1,1));
    analyt_deriv_visc_mean(1,2) =      (analyt_vderxy_mean_visc(1,2) + analyt_vderxy_mean_visc(2,1));
    analyt_deriv_visc_mean(2,2) = 2.0* (analyt_vderxy_mean_visc(2,2));
    analyt_deriv_visc_mean(1,0) = analyt_deriv_visc_mean(0,1);
    analyt_deriv_visc_mean(2,0) = analyt_deriv_visc_mean(0,2);
    analyt_deriv_visc_mean(2,1) = analyt_deriv_visc_mean(1,2);

    // calculate difference between analytical solution and converged solution in function and derivatives

    static LINALG::Matrix<nsd,1> err_vel_jump;			// [v]-[v_h]
    static LINALG::Matrix<nsd,nsd> err_deriv_vel_visc_mean;	// {2*mu*E(v-v_h)} = {2*mu*E(v)}-{2*mu*E(v_h)}
    static LINALG::Matrix<nsd,1> err_visc2_vel_mean_normal; // {2*mu*E(v-v_h)}*n

    err_vel_jump.Clear();
    err_deriv_vel_visc_mean.Clear();
    err_visc2_vel_mean_normal.Clear();


    // get difference between exact and approximative solution

    //loop over Velx, Vely, Velz (size_t i)
    for (size_t i = 0; i < nsd; i ++ )
    {
      // get [v]-[v_h]
      err_vel_jump(i) = analyt_vel_jump(i)-vjump(i);

      // get {2*mu*E(v)-2*mu*E(v_h)}
      // 0=d/dx,1=d/dy,2=d/dz
      err_deriv_vel_visc_mean(i,0) = analyt_deriv_visc_mean(i,0) - deriv_visc_mean(i,0); // error of d/dx
      err_deriv_vel_visc_mean(i,1) = analyt_deriv_visc_mean(i,1) - deriv_visc_mean(i,1); // error of d/dy
      err_deriv_vel_visc_mean(i,2) = analyt_deriv_visc_mean(i,2) - deriv_visc_mean(i,2); // error of d/dz
    }



    // get {2*mu*E(v)-2*mu*E(v_h)}*n
    for (size_t i = 0; i< nsd; i++){
      err_visc2_vel_mean_normal(i,0) = err_deriv_vel_visc_mean(i,0)*normal(0)
                                         + err_deriv_vel_visc_mean(i,1)*normal(1)
                                         + err_deriv_vel_visc_mean(i,2)*normal(2);
    }

    // update element errors
    for (size_t i = 0; i< nsd; ++i)
    {
      eleViscInterfErr		+= hk_eleDiam * err_visc2_vel_mean_normal(i,0) * err_visc2_vel_mean_normal(i,0) * fac;
      eleNitscheErr			+= hk_eleDiam * err_visc2_vel_mean_normal(i,0) * err_visc2_vel_mean_normal(i,0) * fac;
      eleVelJumpInterfErr	+= 1.0/hk_eleDiam * err_vel_jump(i) * err_vel_jump(i) * fac;
      eleNitscheErr		    += 1.0/hk_eleDiam * err_vel_jump(i) * err_vel_jump(i) * fac;
    }

  } // end loop over Gauß points

  return;
}


/*!
  Calculate Nitsche errors with domain integrals
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2>
void Nitsche_BuildBoundaryIntegratedErrors(
    ParameterList&                 eleparams,       /// list of parameters with actions to be calculated
    const INPAR::COMBUST::NitscheError&  NitscheErrorType, /// which solution is the reference
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih, ///< information about the interface
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const M1&                      evelnp,       ///< nodal velocity n+1 values for this element
    const V1&                      eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                      ephi,         ///< nodal G-function values for this element
    const M2&                      grad_ephi_smoothed,
    Teuchos::RCP<const MAT::Material> material,  ///< fluid material
    double&                        ele_meas_plus,
    double&                        ele_meas_minus,
    const bool                     smoothed_boundary_integration
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion Nitsche errors - boundary");

  //==================================================================================================  
  // get actual values of norms^2 of domain errors,stored in eleparams, then update the values with error part of this element
  // at the end of this function we set the eleparams new
  double eleViscInterfErr 		= eleparams.get<double>("H-1/2 integrated viscosity interface error");
  double eleVelJumpInterfErr 	= eleparams.get<double>("H1/2 integrated velocity jump interface error");
  double eleNitscheErr			= eleparams.get<double>("Nitsche integrated error");


  // here are stored the !!!squared norms!!!
  //==================================================================================================

  // function to evaluate analytical solution at Gaussian point -> return velx,vely,velz,pres at Gaussian point
  // pointer to function
  void (*analytical_Interface)(	LINALG::Matrix<3,1>& posXYZgp,
      LINALG::Matrix<3,1>& analyt_vel_plus,
      LINALG::Matrix<3,1>& analyt_vel_minus,
      LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
      LINALG::Matrix<3,3>& analyt_vel_deriv_minus);
  analytical_Interface=NULL;

  // set pointer to function which evaluates the analytical solution
  if     (NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnx1)
    analytical_Interface = &AnalyticalSol_static_bubble_nxnx1_Interface;	// static bubble with n x n x 1 Elements
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnxn)
    analytical_Interface = &AnalyticalSol_static_bubble_nxnxn_Interface;  // static bubble with n x n x n Elements
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_couette_20x20x1)
    analytical_Interface = &AnalyticalSol_Couette_20x20x1_plain_Interface;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_shear)
    analytical_Interface = &AnalyticalSol_shear_Interface;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_straight_bodyforce)
    analytical_Interface = &AnalyticalSol_straight_bodyforce_inclined_Interface;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_ellipsoid_bubble_3D)
    analytical_Interface = &AnalyticalSol_ellipsoid_bubble_nxnxn_instat_Interface;
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_ellipsoid_bubble_2D)
    analytical_Interface = &AnalyticalSol_ellipsoid_bubble_nxnx1_instat_Interface;
  else dserror("there is no analytical solution for 'flowproblem' stored in eleparams");


  //==================================================================================================
  //update the errors with contribution of this element ele

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  //------------------------------------------------------------------------------------------------
  // get material parameters for all boundary integration cell of this element
  //------------------------------------------------------------------------------------------------
#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif

  // get material list for this element
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
  // index plus == burnt material
  const int matid_plus = 3;
  // get material from list of materials
  Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid_plus);
  // check if we really have a fluid material
  dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
  // get the kinematic viscosity \mu
  const double visc_plus = mat->Viscosity();


  // index minus = unburnt material
  const int matid_minus = 4;
  // get material from list of materials
  matptr = matlist->MaterialById(matid_minus);
  // check if we really have a fluid material
  dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
  // get the dynamic viscosity \mu
  const double visc_minus = mat->Viscosity();


  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // calculate element diameter for Nitsche norms
  const double hk_eleDiam = COMBUST::getEleDiameter<DISTYPE>(xyze);

  //=====================================================================================================
  // get weights for average operators { () } = kappa_plus * ()_plus + kappa_minus * ()_minus
  //                               and < () > = kappa_minus* ()_plus + kappa_plus  * ()_minus

  double ele_meas = ele_meas_plus + ele_meas_minus;
  if(ele_meas == 0.0) dserror("element measure is zero!!! Weighted average operators can not get defined!");
  double kappa_plus  = ele_meas_plus/ele_meas;
  double kappa_minus = ele_meas_minus/ele_meas;

  //=====================================================================================================



  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  //  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());

  //------------------------------------------------------------------------------------------
  // evaluate the enrichment function at the interface (boundary integration cells)
  //------------------------------------------------------------------------------------------
  const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi);
  const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi);


  //------------------------------------------------------------------------------------------------
  // loop over boundary integration cells
  //------------------------------------------------------------------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    switch (cell->Shape())
    {

    case DRT::Element::tri3:
      COMBUST::Nitsche_BuildBoundaryIntegratedErrors_Cell<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
          analytical_Interface,
          eleViscInterfErr, eleVelJumpInterfErr, eleNitscheErr, hk_eleDiam,
          *cell, DRT::UTILS::intrule_tri_37point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi_smoothed,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, visc_plus, visc_minus,
          numparamvelx,smoothed_boundary_integration
      );
      break;
    case DRT::Element::quad4:
      // TODO check this intrule !!!
      COMBUST::Nitsche_BuildBoundaryIntegratedErrors_Cell<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          analytical_Interface,
          eleViscInterfErr, eleVelJumpInterfErr, eleNitscheErr, hk_eleDiam,
          *cell, DRT::UTILS::intrule_quad_9point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi_smoothed,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, visc_plus, visc_minus,
          numparamvelx,smoothed_boundary_integration
      );
      break;
    default:
      dserror("cell distype not implemented yet");
    }


  } // end loop over integration cells

  //==================================================================================================


  //==================================================================================================
  // set the boundary errors and Nitsche integrated error with the updated values

  eleparams.set<double>("H-1/2 integrated viscosity interface error", eleViscInterfErr);
  eleparams.set<double>("H1/2 integrated velocity jump interface error", eleVelJumpInterfErr);
  eleparams.set<double>("Nitsche integrated error", eleNitscheErr);
  //==================================================================================================
  return;
}


} // namespace COMBUST


#endif
#endif
#endif
