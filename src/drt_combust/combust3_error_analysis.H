/*----------------------------------------------------------------------*/
/*!
\file combust3_error_analysis.H

\brief comparison to analytical solution for error estimation
       for premixed combustion and two-phase flow problems / 3D combustion element

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_ERROR_ESTIMATION_H
#define COMBUST3_ERROR_ESTIMATION_H


#include "combust_defines.H"

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
void AnalyticalSol_Couette_20x20x1_plain(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: x+0.025
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 20.0*y+10.0 |
   * u_x(x,y,z)= | 0             |(linear in y, konst in x,z)
   *             | 0             |
   * du_x/dy = 20.0 other derivatives are zero
   *
   * p(x,y,z)=   7.0 for x<-0.025
   *         4.0 for x>-0.025
   * gradp = 0.0
   */

  double exact_levelset_gp = posXYZgp(0)+0.025;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0)); // do nothing
  else std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


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
void AnalyticalSol_Couette_20x20x1_plain_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: x+0.025
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 20.0*y+10.0 |
   * u_x(x,y,z)= | 0             |(linear in y, konst in x,z)
   *             | 0             |
   * du_x/dy = 20.0 other derivatives are zero
   *
   * p(x,y,z)=   7.0 for x<-0.025
   *         4.0 for x>-0.025
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

  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: shear.dat
 */
void AnalyticalSol_shear(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
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
   *   Levelset: x+y+0.05
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
   * p(x,y,z)=   5.0 for x+y<-0.05
   *         2.0 for x+y>-0.05
   * gradp = 0.0
   *
   *
   */

  double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0)); // do nothing
  else std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


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
void AnalyticalSol_shear_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: x+y+0.05
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
   * p(x,y,z)=   5.0 for x+y<-0.05
   *         2.0 for x+y>-0.05
   * gradp = 0.0
   *
   *
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   * FLAGS:
   * #define COMBUST_TESTCOUETTEFLOW
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */

  //  double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  //  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  //  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  //  if((domainplus && exact_levelset_gp > 0) ||
  //     (!domainplus && exact_levelset_gp < 0)); // do nothing
  //  else std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


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

  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: straight_bodyforce.dat
 */
void AnalyticalSol_straight_bodyforce_inclined(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: x+y+0.05
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
   * p(x,y,z)=   5.0 for x+y<-0.05
   *         2.0 for x+y>-0.05
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0)); // do nothing
  else std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


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
void AnalyticalSol_straight_bodyforce_inclined_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: x+y+0.05
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
   * p(x,y,z)=   5.0 for x+y<-0.05
   *         2.0 for x+y>-0.05
   * gradp = 0.0
   *
   */

  //  double exact_levelset_gp = posXYZgp(0)+posXYZgp(1)+0.05;
  //  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  //  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  //  if((domainplus && exact_levelset_gp > 0) ||
  //     (!domainplus && exact_levelset_gp < 0)); // do nothing
  //  else std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";


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

  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   *         2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0)); // do nothing
  else std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

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
void AnalyticalSol_static_bubble_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   *         2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!
  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnx1(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   *         2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.16;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0))
    {} // do nothing
//  else std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = 0.4;
  double curvature = 1.0/rad_1;
  const double surftens = 70.0*1e-003;



  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //  if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04>0) analyt_pres=0.0;
  //  else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if(domainplus)
    analyt_pres=0.0;
  else
    analyt_pres= 0.0+curvature*surftens;



  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnx1_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04<0
   *         2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!
  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnxn(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09<0
   *         2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09<0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.09;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0)); // do nothing
  else
  {

  } //std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = 0.3;
  double rad_2 = 0.3;
  double curvature = 1.0/rad_1 + 1.0/rad_2;
  const double surftens = 70.0*1e-003;


  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //  if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04>0) analyt_pres=0.0;
  //  else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if(domainplus)
    analyt_pres=0.0;
  else
    analyt_pres= 0.0+curvature*surftens;

  return;
}


/*
 * author schott
 * May 27, 2010
 *
 * implements the analytical solution for test case: static_bubble
 */
void AnalyticalSol_static_bubble_nxnxn_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.09
   *
   * [|u|]:=0, [|sigma*n|]:=3, g=0 -> [|-p|]=3
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   5.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.09<0
   *         2.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.09>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!!!!!!!!!!
  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}


/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnxn_instat(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   *         0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0)); // do nothing
  else
  {

  } //std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = 0.232401774;
  double rad_2 = 0.232401774;
  double curvature = 1.0/rad_1 + 1.0/rad_2;
  const double surftens = 70.0*1e-003;


  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //  if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04>0) analyt_pres=0.0;
  //  else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if (domainplus)
    analyt_pres=0.0;
  else
    analyt_pres= 0.0+curvature*surftens;
  // if(domainplus) analyt_pres=4.0;
  // else analyt_pres= 4.0+0.0;


  return;
}

/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnxn_instat_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)(z-0.01)*(z-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   *         0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!!!!!!!!!!
  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}



/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnx1_instat(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)-0.01
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   *         0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  double exact_levelset_gp = (posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)-0.01;
  // check if Gaussian point (calculate in Phi_h-approx Interface) lies at the same side of the real!!! Interface
  // matid = 3 -> burnt domain (+ domain), matid = 4 -> unburnt domain (- domain)
  if((domainplus && exact_levelset_gp > 0) ||
      (!domainplus && exact_levelset_gp < 0)); // do nothing
  else
  {

  } //std::cout << "\n !!! real Interface side (via PHI) and approximativ Interface side (via PHI_h) are not the same for Gaussian point !!! \n";

  double rad_1 = sqrt(0.1);
  //double rad_2 = 0.1;
  double curvature = 1.0/rad_1;
  const double surftens = 70.0*1e-003;


  analyt_vel.Clear();
  analyt_vel_deriv.Clear();

  analyt_pres_deriv.Clear();

  //  if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04>0) analyt_pres=0.0;
  //  else if((posXYZgp(0)-0.01)*(posXYZgp(0)-0.01)+(posXYZgp(1)-0.01)*(posXYZgp(1)-0.01)+(posXYZgp(2)-0.01)*(posXYZgp(2)-0.01)-0.04<0) analyt_pres= 0.0+curvature*surftens;
  //
  // get the solution which is the best on the approximativ interface
  // take the constant curvature of the not approximated interface (1/rad)
  if (domainplus)
    analyt_pres=0.0;
  else
    analyt_pres= 0.0+curvature*surftens;

  return;
}

/*
 * author schott
 * June 22, 2010
 *
 * implements the analytical solution for test case: ellipsoid_bubble
 */
void AnalyticalSol_ellipsoid_bubble_nxnx1_instat_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // evaluate the analytical solution of the flow problem at posXYZgp
  /*  Levelset: (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)(z-0.01)*(z-0.01)-0.04
   *
   * [|u|]:=0, [|sigma*n|]:=2/rad * surftens, g=0 -> [|-p|]=2/rad * surftens
   *             | 0               |
   * u_x(x,y,z)= | 0               |(konst in x,y,z) für Omega^- und Omega^+
   *             | 0               |
   *
   * vel derivatives are zero
   *
   * p(x,y,z)=   2/rad * surftens for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04<0
   *         0.0 for (x-0.01)*(x-0.01)+(y-0.01)*(y-0.01)+(z-0.01)*(z-0.01)-0.04>0
   * gradp = 0.0
   *
   */

  // do nothing, all values are zero!!!!!!!!!!!!!!!

  analyt_pres_plus.Clear();
  analyt_pres_minus.Clear();

  return;
}


/*
 * implements the analytical solution for the Betrami flow
 */
void AnalyticalSol_Beltrami(
    const LINALG::Matrix<3,1>& posXYZgp,
    const bool domainplus,
    const double time,
    const double dynvisc,
    const double dens,
    LINALG::Matrix<3,1>& analyt_vel,
    LINALG::Matrix<3,3>& analyt_vel_deriv,
    LINALG::Matrix<1,1>& analyt_pres,
    LINALG::Matrix<3,1>& analyt_pres_deriv)
{
  const double visc = dynvisc/dens;
  // set constants for analytical solution
  const double a = M_PI/4.0;
  const double d = M_PI/2.0;

  const double t = time;

  // compute analytical pressure
  analyt_pres(0,0) = -a*a/2.0 *
      ( exp(2.0*a*posXYZgp(0))
      + exp(2.0*a*posXYZgp(1))
      + exp(2.0*a*posXYZgp(2))
      + 2.0 * sin(a*posXYZgp(0) + d*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) * exp(a*(posXYZgp(1)+posXYZgp(2)))
      + 2.0 * sin(a*posXYZgp(1) + d*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) * exp(a*(posXYZgp(2)+posXYZgp(0)))
      + 2.0 * sin(a*posXYZgp(2) + d*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) * exp(a*(posXYZgp(0)+posXYZgp(1)))
      )* exp(-2.0*visc*d*d*t);

  // compute analytical velocities
  analyt_vel(0) = -a * ( exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) +
                         exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);
  analyt_vel(1) = -a * ( exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) +
                         exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);
  analyt_vel(2) = -a * ( exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) +
                         exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);

  // du_x/dx
  analyt_vel_deriv(0,0) = -a * ( a*exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) -
                                 a*exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);
  // du_x/dy
  analyt_vel_deriv(0,1) = -a * ( a*exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) -
                                 d*exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);
  // du_x/dz
  analyt_vel_deriv(0,2) = -a * ( d*exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) +
                                 a*exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);

  // du_y/dx
  analyt_vel_deriv(1,0) = -a * ( d*exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) +
                                 a*exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);
  // du_y/dy
  analyt_vel_deriv(1,1) = -a * ( a*exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) -
                                 a*exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);
  // du_y/dz
  analyt_vel_deriv(1,2) = -a * ( a*exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) -
                                 d*exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);

  // du_z/dx
  analyt_vel_deriv(2,0) = -a * ( a*exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) -
                                 d*exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);
  // du_z/dy
  analyt_vel_deriv(2,1) = -a * ( d*exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) +
                                 a*exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);
  // du_z/dz
  analyt_vel_deriv(2,2) = -a * ( a*exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) -
                                 a*exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);

  // not implemented, the pressure gradient is not 0!
  analyt_pres_deriv.Clear();

  return;
}


/*
 * implements the analytical solution for the Betrami flow
 */
void AnalyticalSol_Beltrami_Interface(
    const LINALG::Matrix<3,1>& posXYZgp,
    const double time,
    const double dynviscplus,
    const double dynviscminus,
    const double densplus,
    const double densminus,
    LINALG::Matrix<3,1>& analyt_vel_plus,
    LINALG::Matrix<3,1>& analyt_vel_minus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
    LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
    LINALG::Matrix<1,1>& analyt_pres_plus,
    LINALG::Matrix<1,1>& analyt_pres_minus)
{
  // we assume that material parameters do not jump across the interface (zero jumps)
  if (dynviscplus != dynviscminus)
    dserror("identical dynamic viscosities expected");
  if (densplus != densminus)
    dserror("identical densities expected");
  // chose the plus side to compute kinematic viscosity
  const double visc = dynviscplus/densplus;

  // set constants for analytical solution
  const double a = M_PI/4.0;
  const double d = M_PI/2.0;

  const double t = time;

  // compute analytical velocities
  analyt_vel_plus(0) = -a * ( exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) +
                              exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);
  analyt_vel_plus(1) = -a * ( exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) +
                              exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);
  analyt_vel_plus(2) = -a * ( exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) +
                              exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);

  // du_x/dx
  analyt_vel_deriv_plus(0,0) = -a * ( a*exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) -
                                      a*exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);
  // du_x/dy
  analyt_vel_deriv_plus(0,1) = -a * ( a*exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) -
                                      d*exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);
  // du_x/dz
  analyt_vel_deriv_plus(0,2) = -a * ( d*exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) +
                                      a*exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) ) * exp(-visc*d*d*t);

  // du_y/dx
  analyt_vel_deriv_plus(1,0) = -a * ( d*exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) +
                                      a*exp(a*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);
  // du_y/dy
  analyt_vel_deriv_plus(1,1) = -a * ( a*exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) -
                                      a*exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);
  // du_y/dz
  analyt_vel_deriv_plus(1,2) = -a * ( a*exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) -
                                      d*exp(a*posXYZgp(0)) * sin(a*posXYZgp(1) + d*posXYZgp(2)) ) * exp(-visc*d*d*t);

  // du_z/dx
  analyt_vel_deriv_plus(2,0) = -a * ( a*exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) -
                                      d*exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);
  // du_z/dy
  analyt_vel_deriv_plus(2,1) = -a * ( d*exp(a*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) +
                                      a*exp(a*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);
  // du_z/dz
  analyt_vel_deriv_plus(2,2) = -a * ( a*exp(a*posXYZgp(2)) * sin(a*posXYZgp(0) + d*posXYZgp(1)) -
                                      a*exp(a*posXYZgp(1)) * sin(a*posXYZgp(2) + d*posXYZgp(0)) ) * exp(-visc*d*d*t);

  // assume zero jumps at interface (identical material parameters)
  analyt_vel_minus = analyt_vel_plus;
  analyt_vel_deriv_minus = analyt_vel_deriv_plus;

  analyt_pres_plus(0,0) = -a*a/2.0 *
      ( exp(2.0*a*posXYZgp(0))
      + exp(2.0*a*posXYZgp(1))
      + exp(2.0*a*posXYZgp(2))
      + 2.0 * sin(a*posXYZgp(0) + d*posXYZgp(1)) * cos(a*posXYZgp(2) + d*posXYZgp(0)) * exp(a*(posXYZgp(1)+posXYZgp(2)))
      + 2.0 * sin(a*posXYZgp(1) + d*posXYZgp(2)) * cos(a*posXYZgp(0) + d*posXYZgp(1)) * exp(a*(posXYZgp(2)+posXYZgp(0)))
      + 2.0 * sin(a*posXYZgp(2) + d*posXYZgp(0)) * cos(a*posXYZgp(1) + d*posXYZgp(2)) * exp(a*(posXYZgp(0)+posXYZgp(1)))
      )* exp(-2.0*visc*d*d*t);

  analyt_pres_minus = analyt_pres_plus;

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
    Teuchos::ParameterList&                         eleparams,        ///< list of parameters with actions to be calculated
    const INPAR::COMBUST::NitscheError&    NitscheErrorType, ///< which solution is the reference
    const DRT::ELEMENTS::Combust3*         ele,              ///< element whose matrix is calculated
    const COMBUST::InterfaceHandleCombust* ih,               ///< information about the interface
    const XFEM::ElementDofManager&         dofman,           ///< dofmanager of this element
    const M1&                              evelnp,           ///< nodal velocity n+1 values for this element
    const V1&                              eprenp,           ///< nodal pressure n+1 values for this element
    const V2&                              ephi,             ///< nodal G-function values for this element
    Teuchos::RCP<const MAT::Material>      material,         ///< fluid material
    const double                           time,             ///< current time
    double&                                ele_meas_plus,
    double&                                ele_meas_minus,
    const bool                             withkinks         ///< distinguish kinks and jumps
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion Nitsche errors - domain");

  std::cout << "#################################################################################\n";
  std::cout << "## WARNING: You should use the calculation of errorrs with caution! Check all  ##\n";
  std::cout << "##          functionality before using it!                                     ##\n";
  std::cout << "#################################################################################";
  // I am not sure whether this function has been carefully considered when further developing the combustion
  // module. It has been introduced in the DA of Schott and also been applied by Flo Henke for his dissertation. Since
  // I fixed a bug here (February 19th 2015), I would not rely on the results obtained with this function and strongly
  // recommend to check it before using!! Ursula

  //==================================================================================================
  // get actual values of norms^2 of domain errors,stored in eleparams, then update the values with error part of this element
  // at the end of this function we set the eleparams new
  double eleVelDomErr           = eleparams.get<double>("L2 integrated velocity domain error");
  double eleGradVelDomErr       = eleparams.get<double>("L2 integrated grad_velocity domain error");
  double elePresDomErr          = eleparams.get<double>("L2 integrated pressure domain error");
  double eleGradPresDomErr      = eleparams.get<double>("L2 integrated grad_pressure domain error");
  double eleWeightPresDomErr    = eleparams.get<double>("L2 integrated weighted pressure domain error");
  double eleNitscheErr          = eleparams.get<double>("Nitsche integrated error");
  // here are stored the !!!squared norms!!!
  //==================================================================================================

  // function to evaluate analytical solution at Gaussian point -> return velx,vely,velz,pres at Gaussian point
  // pointer to function
  void (*analytical)(
      const LINALG::Matrix<3,1>& posXYZgp,
      const bool domainplus,
      const double time,
      const double dynvisc,
      const double dens,
      LINALG::Matrix<3,1>& analyt_vel,
      LINALG::Matrix<3,3>& analyt_vel_deriv,
      LINALG::Matrix<1,1>& analyt_pres,
      LINALG::Matrix<3,1>& analyt_pres_deriv);
  analytical = NULL;

  // set pointer to function which evaluates the analytical solution

  // static bubble with n x n x 1 Elements
  if (NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnx1)
    analytical = &AnalyticalSol_static_bubble_nxnx1;
  // static bubble with n x n x n Elements
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnxn)
    analytical = &AnalyticalSol_static_bubble_nxnxn;
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
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_beltrami)
    analytical = &AnalyticalSol_Beltrami;
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

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  //=================================================================================
  //-----------------------------------
  // loop over domain integration cells
  //-----------------------------------
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    //--------------------------------------------------
    // get material parameters for this integration cell
    //--------------------------------------------------
    // density \rho^{n+1}
    double densaf  = -1.0;
    // dynamic viscosity \mu^{n+1}
    double dynvisc = -1.0;
    double dynvisc_max = -1.0;
    // set density and dynamic viscosity
    {
      double densplus = -1.0;
      double densminus = -1.0;
      double dynviscplus = -1.0;
      double dynviscminus = -1.0;
      COMBUST::GetMaterialParams(material, true, densplus, dynviscplus);
      COMBUST::GetMaterialParams(material, false, densminus, dynviscminus);

      // blend smoothly to minus (gas bubble) domain
      COMBUST::BlendMaterial(ele, cell->getDomainPlus(), time,
                             densaf, dynvisc, densplus, dynviscplus,
                             densminus, dynviscminus);

      dynvisc_max = std::max(dynviscplus,dynviscminus);
    }

    //-----------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //-----------------------------------------------------------
    Teuchos::RCP<XFEM::ElementEnrichmentValues> enrvals;
    // caution: kink enrichments are evaluated in gauss point loop
    if (!withkinks)
      enrvals = Teuchos::rcp(new XFEM::ElementEnrichmentValues(*ele,dofman,*cell,ephi));

    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
    // TODO: Check Gaussrule
    // we only need to integrate u_h and p_h exact, if the exact solution is smoother we wont get a better
    // order of convergence because u_h and p_h approximates u and p only with this order of convergence
    DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;
    if (!withkinks)
      gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Splited(), cell->Shape());
    // special getXFEMGaussruleKinkEnr for kink enrichment is called as parabolic shape functions are obtained
    // after multipying N and Psi
    else
      gaussrule = XFEM::getXFEMGaussruleKinkEnr<DISTYPE>(ele, xyze, ele->Splited(),cell->Shape());
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //---------------------------------------------
      // transform coordinates of this Gaussian point
      //---------------------------------------------
      // coordinates of the current integration point in cell coordinates \eta^domain
      const LINALG::Matrix<nsd,1> pos_eta_domain(intpoints.qxg[iquad]);

      // coordinates of the current integration point in element coordinates \xi
      static LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
      const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
#ifdef DEBUG
      if (detcell < 0.0)
      {
        std::cout << "detcell :  " << detcell << std::endl;
        dserror("negative detcell!");
      }
#endif
      //----------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //----------------------------------------------------------------------------
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      //-------------------------------------
      // procedures involving Jacobian matrix
      //-------------------------------------
      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
      xjm.MultiplyNT(deriv,xyze);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det*detcell;

      if (det < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
      }

      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Invert(xjm);

      //-------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //-------------------------------------------------------------------
      static LINALG::Matrix<3,numnode> derxy;
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);

      // compute second global derivative
      static LINALG::Matrix<6,numnode> derxy2;
      if (higher_order_ele)
      {
        static LINALG::Matrix<6,numnode> deriv2;
        DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
        DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }

      //-----------------------------------------------------------
      // evaluate the enrichment function for this integration cell
      //-----------------------------------------------------------
      // kink enrichments are called with level-set values
      // jump enrichments are called with domain cells!
      if (withkinks)
        enrvals = Teuchos::rcp(new XFEM::ElementEnrichmentValues(*ele,dofman,ephi,*cell,funct,derxy,derxy2));

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
//      static XFEM::ApproxFunc<2,shpVecSize> shp;
      static XFEM::ApproxFunc<2,shpVecSize> shpvel;
      static XFEM::ApproxFunc<2,shpVecSize> shppres;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSize,1> enr_funct_vel;
        enr_funct_vel.Clear();
        static LINALG::Matrix<3,shpVecSize> enr_derxy_vel;
        enr_derxy_vel.Clear();
        static LINALG::Matrix<6,shpVecSize> enr_derxy2_vel;
        enr_derxy2_vel.Clear();
        static LINALG::Matrix<shpVecSize,1> enr_funct_pres;
        enr_funct_pres.Clear();
        static LINALG::Matrix<3,shpVecSize> enr_derxy_pres;
        enr_derxy_pres.Clear();
        static LINALG::Matrix<6,shpVecSize> enr_derxy2_pres;
        enr_derxy2_pres.Clear();

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals->ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Velx,
            funct,
            derxy,
            derxy2,
            enr_funct_vel,
            enr_derxy_vel,
            enr_derxy2_vel);

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          shpvel.d0(iparam) = enr_funct_vel(iparam);
          shpvel.dx(iparam) = enr_derxy_vel(0,iparam);
          shpvel.dy(iparam) = enr_derxy_vel(1,iparam);
          shpvel.dz(iparam) = enr_derxy_vel(2,iparam);
          shpvel.dxdx(iparam) = enr_derxy2_vel(0,iparam);
          shpvel.dxdy(iparam) = enr_derxy2_vel(3,iparam);
          shpvel.dxdz(iparam) = enr_derxy2_vel(4,iparam);
          shpvel.dydx(iparam) = shpvel.dxdy(iparam);
          shpvel.dydy(iparam) = enr_derxy2_vel(1,iparam);
          shpvel.dydz(iparam) = enr_derxy2_vel(5,iparam);
          shpvel.dzdx(iparam) = shpvel.dxdz(iparam);
          shpvel.dzdy(iparam) = shpvel.dydz(iparam);
          shpvel.dzdz(iparam) = enr_derxy2_vel(2,iparam);
        }

        if(withkinks)
        {
          enrvals->ComputeModifiedEnrichedNodalShapefunction(
              XFEM::PHYSICS::Pres,
              funct,
              derxy,
              derxy2,
              enr_funct_pres,
              enr_derxy_pres,
              enr_derxy2_pres);
          for (size_t iparam = 0; iparam < numparampres; ++iparam)
          {
            shppres.d0(iparam) = enr_funct_pres(iparam);
            shppres.dx(iparam) = enr_derxy_pres(0,iparam);
            shppres.dy(iparam) = enr_derxy_pres(1,iparam);
            shppres.dz(iparam) = enr_derxy_pres(2,iparam);
            shppres.dxdx(iparam) = enr_derxy2_pres(0,iparam);
            shppres.dxdy(iparam) = enr_derxy2_pres(3,iparam);
            shppres.dxdz(iparam) = enr_derxy2_pres(4,iparam);
            shppres.dydx(iparam) = shppres.dxdy(iparam);
            shppres.dydy(iparam) = enr_derxy2_pres(1,iparam);
            shppres.dydz(iparam) = enr_derxy2_pres(5,iparam);
            shppres.dzdx(iparam) = shppres.dxdz(iparam);
            shppres.dzdy(iparam) = shppres.dydz(iparam);
            shppres.dzdz(iparam) = enr_derxy2_pres(2,iparam);
          }
        }
        else
        {
          dsassert(numparamvelx == numparampres, "assumption violation");
          for (size_t iparam = 0; iparam < numparampres; ++iparam)
          {
            shppres.d0(iparam) = enr_funct_vel(iparam);
            shppres.dx(iparam) = enr_derxy_vel(0,iparam);
            shppres.dy(iparam) = enr_derxy_vel(1,iparam);
            shppres.dz(iparam) = enr_derxy_vel(2,iparam);
            shppres.dxdx(iparam) = enr_derxy2_vel(0,iparam);
            shppres.dxdy(iparam) = enr_derxy2_vel(3,iparam);
            shppres.dxdz(iparam) = enr_derxy2_vel(4,iparam);
            shppres.dydx(iparam) = shpvel.dxdy(iparam);
            shppres.dydy(iparam) = enr_derxy2_vel(1,iparam);
            shppres.dydz(iparam) = enr_derxy2_vel(5,iparam);
            shppres.dzdx(iparam) = shpvel.dxdz(iparam);
            shppres.dzdy(iparam) = shpvel.dydz(iparam);
            shppres.dzdz(iparam) = enr_derxy2_vel(2,iparam);
          }
        }
      }
      else // not xfem_assembly i.e. standard assembly
      {
        // fill approximation functions for standard FEM
        // remark: numparamvelx == numnode, for standard FEM
        for (size_t iparam = 0; iparam < numnode; ++iparam)
        {
          shpvel.d0(iparam) = funct(iparam);
          shpvel.dx(iparam) = derxy(0,iparam);
          shpvel.dy(iparam) = derxy(1,iparam);
          shpvel.dz(iparam) = derxy(2,iparam);
          shpvel.dxdx(iparam) = derxy2(0,iparam);
          shpvel.dxdy(iparam) = derxy2(3,iparam);
          shpvel.dxdz(iparam) = derxy2(4,iparam);
          shpvel.dydx(iparam) = shpvel.dxdy(iparam);
          shpvel.dydy(iparam) = derxy2(1,iparam);
          shpvel.dydz(iparam) = derxy2(5,iparam);
          shpvel.dzdx(iparam) = shpvel.dxdz(iparam);
          shpvel.dzdy(iparam) = shpvel.dydz(iparam);
          shpvel.dzdz(iparam) = derxy2(2,iparam);

          shppres.d0(iparam) = funct(iparam);
          shppres.dx(iparam) = derxy(0,iparam);
          shppres.dy(iparam) = derxy(1,iparam);
          shppres.dz(iparam) = derxy(2,iparam);
          shppres.dxdx(iparam) = derxy2(0,iparam);
          shppres.dxdy(iparam) = derxy2(3,iparam);
          shppres.dxdz(iparam) = derxy2(4,iparam);
          shppres.dydx(iparam) = shppres.dxdy(iparam);
          shppres.dydy(iparam) = derxy2(1,iparam);
          shppres.dydz(iparam) = derxy2(5,iparam);
          shppres.dzdx(iparam) = shppres.dxdz(iparam);
          shppres.dzdy(iparam) = shppres.dydz(iparam);
          shppres.dzdz(iparam) = derxy2(2,iparam);
        }
      }
      // set measure of element in plus and minus domain
      if(cell->getDomainPlus() == true) {ele_meas_plus  += fac;}
      else  {ele_meas_minus += fac;}
      //------------------------------------------------
      // get velocities of solution at integration point
      //------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shpvel.d0, numparamvelx);

      //----------------------------------------------------------
      // get velocity derivatives of solution at integration point
      //----------------------------------------------------------
      static LINALG::Matrix<3,nsd> vderxy;
      vderxy.Clear();
      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          vderxy(isd,0) += evelnp(isd,iparam) * shpvel.dx(iparam);
          vderxy(isd,1) += evelnp(isd,iparam) * shpvel.dy(iparam);
          vderxy(isd,2) += evelnp(isd,iparam) * shpvel.dz(iparam);
        }
      }

      //----------------------------------
      // get pressure at integration point
      //----------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shppres.d0(iparam);

      //--------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shppres.dx(iparam);
        gradp(1) += eprenp(iparam) * shppres.dy(iparam);
        gradp(2) += eprenp(iparam) * shppres.dz(iparam);
      }

      //=========================== ANALYTICAL SOLUTION ===============
      // get the analytical solution at Gaussian point
      // variables for analytical solution
      static LINALG::Matrix<nsd,1>   analyt_vel;
      static LINALG::Matrix<nsd,nsd> analyt_vel_derxy;
      static LINALG::Matrix<1,1>     analyt_pres;
      static LINALG::Matrix<nsd,1>   analyt_pres_derxy;
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
        COMBUST::GetGPinXYZ<DRT::Element::tet4>(*cell, pos_eta_domain, posGPinXYZ);
        break;
      case DRT::Element::hex8:
        COMBUST::GetGPinXYZ<DRT::Element::hex8>(*cell, pos_eta_domain, posGPinXYZ);
        break;
      default:
        dserror("cell distype not implemented yet");
      }

      analytical(posGPinXYZ, cell->getDomainPlus(), time, dynvisc, densaf, analyt_vel, analyt_vel_derxy, analyt_pres, analyt_pres_derxy);

      //============================== NITSCHE ERRORS =================================
      // calculate difference between analytical solution and converged solution in function and derivatives

      static LINALG::Matrix<nsd,1>   err_vel;       // v-v_h
      static LINALG::Matrix<nsd,nsd> err_vel_derxy; // grad(v-v_h)
      static LINALG::Matrix<1,1>     err_pres;      // p-p_h
      static LINALG::Matrix<nsd,1> err_pres_derxy;  // grad(p-p_h)
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
        eleVelDomErr  += err_vel(i) * err_vel(i)*fac;
        eleNitscheErr += err_vel(i) * err_vel(i)*fac;

        eleGradVelDomErr += dynvisc * err_vel_derxy(0,i) * err_vel_derxy(0,i) * fac;
        eleGradVelDomErr += dynvisc * err_vel_derxy(1,i) * err_vel_derxy(1,i) * fac;
        eleGradVelDomErr += dynvisc * err_vel_derxy(2,i) * err_vel_derxy(2,i) * fac;
        eleNitscheErr    += dynvisc * err_vel_derxy(0,i) * err_vel_derxy(0,i) * fac;
        eleNitscheErr    += dynvisc * err_vel_derxy(1,i) * err_vel_derxy(1,i) * fac;
        eleNitscheErr    += dynvisc * err_vel_derxy(2,i) * err_vel_derxy(2,i) * fac;

        eleGradPresDomErr += err_pres_derxy(i,0) * err_pres_derxy(i,0) * fac;
      }

      elePresDomErr += err_pres(0,0) * err_pres(0,0) * fac;

      eleWeightPresDomErr += 1.0/dynvisc_max * err_pres(0,0) * err_pres(0,0) * fac;
      eleNitscheErr       += 1.0/dynvisc_max * err_pres(0,0) * err_pres(0,0) * fac;
    } // end loop over Gauß points
  } // end loop over integration cells

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
    void (*analytical_Interface)( const LINALG::Matrix<3,1>& posXYZgp,
                                  const double time,
                                  const double dynviscplus,
                                  const double dynviscminus,
                                  const double densplus,
                                  const double densminus,
                                  LINALG::Matrix<3,1>& analyt_vel_plus,
                                  LINALG::Matrix<3,1>& analyt_vel_minus,
                                  LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
                                  LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
                                  LINALG::Matrix<1,1>& analyt_pres_plus,
                                  LINALG::Matrix<1,1>& analyt_pres_minus),
    double&                              eleViscInterfErr,    /// "H-1/2 integrated viscosity interface error"
    double&                              eleVelJumpInterfErr, /// "H1/2 integrated velocity jump interface error"
    double&                              eleFluxJumpInterfErr,/// "H-1/2 integrated flux jump interface error"
    double&                              eleNitscheErr,       /// "Nitsche integrated error"
    const double&                        hk_eleDiam,
    GEO::BoundaryIntCell                 cell,
    const DRT::UTILS::GaussRule2D&       intrule2D,
    M2&                                  xyze,
    const DRT::ELEMENTS::Combust3*       ele,                ///< the element those matrix is calculated
    const bool&                          higher_order_ele,
    const M1&                            evelnp,
    const V1&                            eprenp,
    const V2&                            ephi,
    const M3&                            grad_ephi_smoothed,
    const XFEM::ElementEnrichmentValues& enrvals_plus,
    const XFEM::ElementEnrichmentValues& enrvals_minus,
    const double&                        kappa_plus,     /// weights of average in plus domain
    const double&                        kappa_minus,    /// weights of average in minus domain
    const double&                        dynvisc_plus,   /// viscosity in plus domain
    const double&                        dynvisc_minus,  /// viscosity in minus domain
    const double&                        dens_plus,      /// density in plus domain
    const double&                        dens_minus,     /// density in minus domain
    const double                         time,
    const size_t&                        numparamvelx,
    const bool                           smoothed_boundary_integration,
    const bool                           withkinks
)
{
  const size_t nsd =3;

  // number of nodes of element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  //----------------------------------------------
  // get Gaussian points for this integration cell
  //----------------------------------------------
  const DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

  // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
  LINALG::SerialDenseMatrix cellXiDomaintmp = cell.CellNodalPosXiDomain();
  // get coordinates of vertices of boundary integration cell in physical coordinates \x^domain
  LINALG::SerialDenseMatrix cellXYZDomaintmp = cell.CellNodalPosXYZ();
  // transform to fixed size format
  const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);
  const LINALG::Matrix<nsd,numvertices> cellXYZDomain(cellXYZDomaintmp);

  //--------------------------------------
  // integration loop over Gaussian points
  //--------------------------------------
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
    // remark: for linear boundary integration cells this could be done before the loop over all
    //         Gaussian points
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> normal(true);
    normal.Clear();
    // normal vector points from Omega^+ to Omega^- (n = n^+ = -n^-)
    GEO::computeNormalToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal);

    static LINALG::Matrix<nsd,1> normal_smooth(true);
    normal_smooth.Clear();

    // call smoothed normals for an accurate evaluation of boundary integration terms
    if (smoothed_boundary_integration == true)
    {
      // get shape functions at gaussian point
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
      DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      // loop over nodes
      for(size_t i = 0; i< numnode; i++)
      {
        LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
        nodal_grad_xyz_tmp.Clear();

        // compute smoothed normal vector
        normal_smooth(0)+=grad_ephi_smoothed(0,i)*funct_gp(i);
        normal_smooth(1)+=grad_ephi_smoothed(1,i)*funct_gp(i);
        normal_smooth(2)+=grad_ephi_smoothed(2,i)*funct_gp(i);
      }
#ifdef COMBUST_2D
      normal_smooth(2) = 0.0;
#endif
      // normalize final gradient
      // remark: n = -grad(\phi)
      const double norm = normal_smooth.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0) dserror("norm of normal vector is zero!");
      normal_smooth.Scale(-1.0/norm);
    }
    else
    {
      normal_smooth = normal;
    }

    //----------------------------------------------------------------------------
    // evaluate shape functions and their first derivatives at this Gaussian point
    //----------------------------------------------------------------------------
    static LINALG::Matrix<numnode,1> funct;
    funct.Clear();
    static LINALG::Matrix<nsd,numnode> deriv;
    deriv.Clear();
    DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
    DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    //--------------------------------------------------------
    // procedures involving Jacobian matrix for domain mapping
    //--------------------------------------------------------
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

    //-------------------------------------------------------------------
    // compute global derivates of shape functions at this Gaussian point
    //-------------------------------------------------------------------
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
      DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
    }
    else
    {
      derxy2.Clear();
    }

    //--------------------------------------------------------------------------------
    // rearrange (enriched) shape functions and derivatives as approximation functions
    //--------------------------------------------------------------------------------
    const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
    static LINALG::Matrix<shpVecSize,1> enrfunct_plus;
    enrfunct_plus.Clear();
    static LINALG::Matrix<shpVecSize,1> enrfunct_minus;
    enrfunct_minus.Clear();

    static LINALG::Matrix<3,shpVecSize> enrderxy_plus;
    enrderxy_plus.Clear();
    static LINALG::Matrix<3,shpVecSize> enrderxy_minus;
    enrderxy_minus.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_plus;
    enrderxy2_plus.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_minus;
    enrderxy2_minus.Clear();

    if (!withkinks)
    {
      // shape functions for nodal parameters (dofs) on plus and minus side
      enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_plus, enrderxy_plus, enrderxy2_plus);
      enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_minus, enrderxy_minus, enrderxy2_minus);
    }
    else
    {
      // shape functions for nodal parameters (dofs) on plus and minus side
      enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Pres, funct, derxy, derxy2,
          enrfunct_plus, enrderxy_plus, enrderxy2_plus);
      enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Pres, funct, derxy, derxy2,
          enrfunct_minus, enrderxy_minus, enrderxy2_minus);
    }

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<0,shpVecSize> shp_jump;
    static XFEM::ApproxFunc<1,shpVecSize> shp_jump_visc;
    static XFEM::ApproxFunc<1,shpVecSize> shp_mean_visc;

    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_jump.d0(iparam) = enrfunct_minus(iparam) - enrfunct_plus(iparam);

      // not used, but filled for security
      shp_jump_visc.d0(iparam) = dynvisc_minus*enrfunct_minus(iparam) - dynvisc_plus*enrfunct_plus(iparam);

      shp_jump_visc.dx(iparam) = dynvisc_minus*enrderxy_minus(0,iparam) - dynvisc_plus*enrderxy_plus(0,iparam);
      shp_jump_visc.dy(iparam) = dynvisc_minus*enrderxy_minus(1,iparam) - dynvisc_plus*enrderxy_plus(1,iparam);
      shp_jump_visc.dz(iparam) = dynvisc_minus*enrderxy_minus(2,iparam) - dynvisc_plus*enrderxy_plus(2,iparam);

      shp_mean_visc.dx(iparam) = kappa_plus * dynvisc_plus*enrderxy_plus(0,iparam) + kappa_minus * dynvisc_minus*enrderxy_minus(0,iparam);
      shp_mean_visc.dy(iparam) = kappa_plus * dynvisc_plus*enrderxy_plus(1,iparam) + kappa_minus * dynvisc_minus*enrderxy_minus(1,iparam);
      shp_mean_visc.dz(iparam) = kappa_plus * dynvisc_plus*enrderxy_plus(2,iparam) + kappa_minus * dynvisc_minus*enrderxy_minus(2,iparam);
    }

    //---------------------------------------
    // compute data at Gaussian point for rhs
    //---------------------------------------

    // evaluate velocity jumps for jump enrichments only
    // if pressure jumps should be evaluated, they have to be added here
    if (!withkinks)
    {
      // velocity jump
      static LINALG::Matrix<nsd,1> vjump(true);
      vjump.Clear();
      vjump = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_jump.d0, numparamvelx);

      // pressure jump
      static double presjump;
      presjump = 0.0;
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
         presjump += eprenp(iparam)*shp_jump.d0(iparam);

      // get velocity (np,i) derivatives at integration point
      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      static LINALG::Matrix<nsd,nsd> vderxy_jump_visc(true);
      vderxy_jump_visc.Clear();
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          vderxy_jump_visc(isd,0) += evelnp(isd,iparam) * shp_jump_visc.dx(iparam);
          vderxy_jump_visc(isd,1) += evelnp(isd,iparam) * shp_jump_visc.dy(iparam);
          vderxy_jump_visc(isd,2) += evelnp(isd,iparam) * shp_jump_visc.dz(iparam);
        }

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
      static LINALG::Matrix<nsd,1>   analyt_vel_plus;
      static LINALG::Matrix<nsd,1>   analyt_vel_minus;
      static LINALG::Matrix<nsd,nsd> analyt_deriv_plus;
      static LINALG::Matrix<nsd,nsd> analyt_deriv_minus;
      static LINALG::Matrix<1,1>     analyt_pres_plus;
      static LINALG::Matrix<1,1>     analyt_pres_minus;

      analyt_vel_plus.Clear();
      analyt_vel_minus.Clear();
      analyt_deriv_plus.Clear();
      analyt_deriv_minus.Clear();
      analyt_pres_plus.Clear();
      analyt_pres_minus.Clear();

      // evaluate the analytical solution
      static LINALG::Matrix<nsd,1> posEtaBoundaryinXYZ;
      posEtaBoundaryinXYZ.Multiply(xyze,funct);

      analytical_Interface(posEtaBoundaryinXYZ, time, dynvisc_plus, dynvisc_minus, dens_plus, dens_minus,
                           analyt_vel_plus, analyt_vel_minus, analyt_deriv_plus, analyt_deriv_minus,
                           analyt_pres_plus, analyt_pres_minus);

      //============================== NITSCHE ERRORS =================================

      // get velocity jump [[v]]
      static LINALG::Matrix<nsd,1> analyt_vel_jump;
      analyt_vel_jump.Clear();
      for (size_t i = 0; i< nsd; i++)
      {
        analyt_vel_jump(i) = analyt_vel_minus(i)-analyt_vel_plus(i);
      }

      // get flux jump [[sigma]]*normal
      static LINALG::Matrix<nsd,1> analyt_flux_jump;
      analyt_flux_jump.Clear();

      analyt_flux_jump(0) = - (analyt_pres_minus(0,0) - analyt_pres_plus(0,0))*normal(0)
                            + (analyt_deriv_minus(0,0)+analyt_deriv_minus(0,0) - (analyt_deriv_plus(0,0)+analyt_deriv_plus(0,0)))*normal(0)
                            + (analyt_deriv_minus(0,1)+analyt_deriv_minus(1,0) - (analyt_deriv_plus(0,1)+analyt_deriv_plus(1,0)))*normal(1)
                            + (analyt_deriv_minus(0,2)+analyt_deriv_minus(2,0) - (analyt_deriv_plus(0,2)+analyt_deriv_plus(2,0)))*normal(2);
      analyt_flux_jump(1) = - (analyt_pres_minus(0,0) - analyt_pres_plus(0,0))*normal(1)
                            + (analyt_deriv_minus(1,0)+analyt_deriv_minus(0,1) - (analyt_deriv_plus(1,0)+analyt_deriv_plus(0,1)))*normal(0)
                            + (analyt_deriv_minus(1,1)+analyt_deriv_minus(1,1) - (analyt_deriv_plus(1,1)+analyt_deriv_plus(1,1)))*normal(1)
                            + (analyt_deriv_minus(1,2)+analyt_deriv_minus(2,1) - (analyt_deriv_plus(1,2)+analyt_deriv_plus(2,1)))*normal(2);
      analyt_flux_jump(2) = - (analyt_pres_minus(0,0) - analyt_pres_plus(0,0))*normal(2)
                            + (analyt_deriv_minus(2,0)+analyt_deriv_minus(0,2) - (analyt_deriv_plus(2,0)+analyt_deriv_plus(0,2)))*normal(0)
                            + (analyt_deriv_minus(2,1)+analyt_deriv_minus(1,2) - (analyt_deriv_plus(2,1)+analyt_deriv_plus(1,2)))*normal(1)
                            + (analyt_deriv_minus(2,2)+analyt_deriv_minus(2,2) - (analyt_deriv_plus(2,2)+analyt_deriv_plus(2,2)))*normal(2);

      // flux jump
      static LINALG::Matrix<nsd,1> fluxjump(true);
      fluxjump.Clear();

      fluxjump(0) = - presjump*normal(0)
                    + (vderxy_jump_visc(0,0)+vderxy_jump_visc(0,0))*normal(0)
                    + (vderxy_jump_visc(0,1)+vderxy_jump_visc(1,0))*normal(1)
                    + (vderxy_jump_visc(0,2)+vderxy_jump_visc(2,0))*normal(2);
      fluxjump(1) = - presjump*normal(1)
                    + (vderxy_jump_visc(1,0)+vderxy_jump_visc(0,1))*normal(0)
                    + (vderxy_jump_visc(1,1)+vderxy_jump_visc(1,1))*normal(1)
                    + (vderxy_jump_visc(1,2)+vderxy_jump_visc(2,1))*normal(2);
      fluxjump(2) = - presjump*normal(2)
                    + (vderxy_jump_visc(2,0)+vderxy_jump_visc(0,2))*normal(0)
                    + (vderxy_jump_visc(2,1)+vderxy_jump_visc(1,2))*normal(1)
                    + (vderxy_jump_visc(2,2)+vderxy_jump_visc(2,2))*normal(2);

      // get {mu * v_x,y} => 3x3
      static LINALG::Matrix<nsd,nsd> analyt_vderxy_mean_visc;
      analyt_vderxy_mean_visc.Clear();
      for(size_t isd = 0; isd < nsd; isd++)
      {
        analyt_vderxy_mean_visc(isd,0) = 0.5* (dynvisc_plus * analyt_deriv_plus(isd,0) + dynvisc_minus * analyt_deriv_minus(isd,0));
        analyt_vderxy_mean_visc(isd,1) = 0.5* (dynvisc_plus * analyt_deriv_plus(isd,1) + dynvisc_minus * analyt_deriv_minus(isd,1));
        analyt_vderxy_mean_visc(isd,2) = 0.5* (dynvisc_plus * analyt_deriv_plus(isd,2) + dynvisc_minus * analyt_deriv_minus(isd,2));
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

     static LINALG::Matrix<nsd,1> err_vel_jump;              // [v]-[v_h]
     static LINALG::Matrix<nsd,1> err_flux_jump;             // [sigma]-[sigma_h]
     static LINALG::Matrix<nsd,nsd> err_deriv_vel_visc_mean; // {2*mu*E(v-v_h)} = {2*mu*E(v)}-{2*mu*E(v_h)}
     static LINALG::Matrix<nsd,1> err_visc2_vel_mean_normal; // {2*mu*E(v-v_h)}*n

     err_vel_jump.Clear();
     err_flux_jump.Clear();
     err_deriv_vel_visc_mean.Clear();
     err_visc2_vel_mean_normal.Clear();

      // get difference between exact and approximative solution

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

     //loop over Velx, Vely, Velz (size_t i)
     for (size_t i = 0; i < nsd; i ++ )
     {
       // get [[v]]-[[v_h]]
       err_vel_jump(i) = analyt_vel_jump(i)-vjump(i);
       err_flux_jump(i) = analyt_flux_jump(i)-fluxjump(i);
       // get {2*mu*E(v)-2*mu*E(v_h)}
       // 0=d/dx,1=d/dy,2=d/dz
       err_deriv_vel_visc_mean(i,0) = analyt_deriv_visc_mean(i,0) - deriv_visc_mean(i,0); // error of d/dx
       err_deriv_vel_visc_mean(i,1) = analyt_deriv_visc_mean(i,1) - deriv_visc_mean(i,1); // error of d/dy
       err_deriv_vel_visc_mean(i,2) = analyt_deriv_visc_mean(i,2) - deriv_visc_mean(i,2); // error of d/dz
     }

//std::cout << "analyt " << analyt_flux_jump << std::endl;

      // get {2*mu*E(v)-2*mu*E(v_h)}*n
      for (size_t i = 0; i< nsd; i++)
      {
        err_visc2_vel_mean_normal(i,0) = err_deriv_vel_visc_mean(i,0)*normal(0)
                                       + err_deriv_vel_visc_mean(i,1)*normal(1)
                                       + err_deriv_vel_visc_mean(i,2)*normal(2);
      }

      // update element errors
      for (size_t i = 0; i< nsd; ++i)
      {
        eleViscInterfErr     += hk_eleDiam * err_visc2_vel_mean_normal(i,0) * err_visc2_vel_mean_normal(i,0) * fac;
        eleNitscheErr        += hk_eleDiam * err_visc2_vel_mean_normal(i,0) * err_visc2_vel_mean_normal(i,0) * fac;
        eleVelJumpInterfErr  += 1.0/hk_eleDiam * err_vel_jump(i) * err_vel_jump(i) * fac;
        eleFluxJumpInterfErr += hk_eleDiam * err_flux_jump(i) * err_flux_jump(i) * fac;
        eleNitscheErr        += 1.0/hk_eleDiam * err_vel_jump(i) * err_vel_jump(i) * fac;
      }
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
    Teuchos::ParameterList&                eleparams,         ///< list of parameters with actions to be calculated
    const INPAR::COMBUST::NitscheError&    NitscheErrorType,  ///< which solution is the reference
    const DRT::ELEMENTS::Combust3*         ele,               ///< element whose matrix is calculated
    const COMBUST::InterfaceHandleCombust* ih,                ///< information about the interface
    const XFEM::ElementDofManager&         dofman,            ///< dofmanager of this element
    const M1&                              evelnp,            ///< nodal velocity n+1 values for this element
    const V1&                              eprenp,            ///< nodal pressure n+1 values for this element
    const V2&                              ephi,              ///< nodal G-function values for this element
    const M2&                              grad_ephi,         ///<
    Teuchos::RCP<const MAT::Material>      material,          ///< fluid material
    const double                           time,
    const double                           ele_meas_plus,
    const double                           ele_meas_minus,
    const bool                             smoothed_boundary_integration,
    const bool                             withkinks
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion Nitsche errors - boundary");

  //==================================================================================================
  // get actual values of norms^2 of domain errors,stored in eleparams, then update the values with error part of this element
  // at the end of this function we set the eleparams new
  double eleViscInterfErr     = eleparams.get<double>("H-1/2 integrated viscosity interface error");
  double eleVelJumpInterfErr  = eleparams.get<double>("H1/2 integrated velocity jump interface error");
  double eleFluxJumpInterfErr = eleparams.get<double>("H-1/2 integrated flux jump interface error");
  double eleNitscheErr        = eleparams.get<double>("Nitsche integrated error");

  // here are stored the !!!squared norms!!!
  //==================================================================================================

  // function to evaluate analytical solution at Gaussian point -> return velx,vely,velz,pres at Gaussian point
  // pointer to function
  void (*analytical_Interface)(
      const LINALG::Matrix<3,1>& posXYZgp,
      const double time,
      const double dynviscplus,
      const double dynviscminus,
      const double densplus,
      const double densminus,
      LINALG::Matrix<3,1>& analyt_vel_plus,
      LINALG::Matrix<3,1>& analyt_vel_minus,
      LINALG::Matrix<3,3>& analyt_vel_deriv_plus,
      LINALG::Matrix<3,3>& analyt_vel_deriv_minus,
      LINALG::Matrix<1,1>& analyt_pres_plus,
      LINALG::Matrix<1,1>& analyt_pres_minus);
  analytical_Interface=NULL;

  // set pointer to function which evaluates the analytical solution
  // static bubble with n x n x 1 Elements
  if (NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnx1)
    analytical_Interface = &AnalyticalSol_static_bubble_nxnx1_Interface;
  // static bubble with n x n x n Elements
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_static_bubble_nxnxn)
    analytical_Interface = &AnalyticalSol_static_bubble_nxnxn_Interface;
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
  else if(NitscheErrorType == INPAR::COMBUST::nitsche_error_beltrami)
    analytical_Interface = &AnalyticalSol_Beltrami_Interface;
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

  //---------------------------------------------------------------------------
  // get material parameters for all boundary integration cells of this element
  //---------------------------------------------------------------------------
  // density \rho
  double dens_plus = -1.0;
  // dynamic viscosity \mu
  double dynvisc_plus = -1.0;
  // density \rho
  double dens_minus = -1.0;
  // dynamic viscosity \mu
  double dynvisc_minus = -1.0;
  // set density and dynamic viscosity
  COMBUST::GetMaterialParams(material, dens_plus, dynvisc_plus, dens_minus, dynvisc_minus);

  // calculate element diameter for Nitsche norms
  const double hk_eleDiam = COMBUST::getEleDiameter<DISTYPE>(xyze);

  //-----------------------------------------------------------------------------------------
  // get weights for average operators { () } = kappa_plus * ()_plus + kappa_minus * ()_minus
  //                               and < () > = kappa_minus* ()_plus + kappa_plus  * ()_minus
  //-----------------------------------------------------------------------------------------
  double ele_meas = ele_meas_plus + ele_meas_minus;
  if(ele_meas == 0.0) dserror("element measure is zero!!! Weighted average operators can not get defined!");
  double kappa_plus  = ele_meas_plus/ele_meas;
  double kappa_minus = ele_meas_minus/ele_meas;

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->ElementBoundaryIntCells(ele->Id());

  //-------------------------------------------------------------------------------
  // evaluate the enrichment function at the interface (boundary integration cells)
  //-------------------------------------------------------------------------------
  const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi);
  const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi);

  //-------------------------------------
  // loop over boundary integration cells
  //-------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    switch (cell->Shape())
    {
    case DRT::Element::tri3:
    {
      if (!withkinks)
        COMBUST::Nitsche_BuildBoundaryIntegratedErrors_Cell<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
          analytical_Interface,
          eleViscInterfErr, eleVelJumpInterfErr, eleFluxJumpInterfErr, eleNitscheErr, hk_eleDiam,
          *cell, DRT::UTILS::intrule_tri_37point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus, dens_plus, dens_minus,
          time, numparamvelx, smoothed_boundary_integration, withkinks);
      else
         COMBUST::Nitsche_BuildBoundaryIntegratedErrors_Cell<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
           analytical_Interface,
           eleViscInterfErr, eleVelJumpInterfErr, eleFluxJumpInterfErr, eleNitscheErr, hk_eleDiam,
           *cell, DRT::UTILS::intrule_tri_37point, xyze, ele, higher_order_ele,
           evelnp, eprenp, ephi, grad_ephi,
           enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus, dens_plus, dens_minus,
           time, numparampres, smoothed_boundary_integration, withkinks);
      break;
    }
    case DRT::Element::quad4:
    {
      if (!withkinks)
        COMBUST::Nitsche_BuildBoundaryIntegratedErrors_Cell<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          analytical_Interface,
          eleViscInterfErr, eleVelJumpInterfErr, eleFluxJumpInterfErr, eleNitscheErr, hk_eleDiam,
          *cell, DRT::UTILS::intrule_quad_25point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus, dens_plus, dens_minus,
          time, numparamvelx, smoothed_boundary_integration, withkinks);
      else
        COMBUST::Nitsche_BuildBoundaryIntegratedErrors_Cell<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          analytical_Interface,
          eleViscInterfErr, eleVelJumpInterfErr, eleFluxJumpInterfErr, eleNitscheErr, hk_eleDiam,
          *cell, DRT::UTILS::intrule_quad_25point, xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, grad_ephi,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus, dens_plus, dens_minus,
          time, numparampres, smoothed_boundary_integration, withkinks);
      break;
    }
    default:
      dserror("cell distype not implemented yet");
    }
  } // end loop over integration cells

  //==================================================================================================
  // set the boundary errors and Nitsche integrated error with the updated values
  eleparams.set<double>("H-1/2 integrated viscosity interface error", eleViscInterfErr);
  eleparams.set<double>("H1/2 integrated velocity jump interface error", eleVelJumpInterfErr);
  eleparams.set<double>("H-1/2 integrated flux jump interface error", eleFluxJumpInterfErr);
  eleparams.set<double>("Nitsche integrated error", eleNitscheErr);
  //==================================================================================================
  return;
}


} // namespace COMBUST


#endif
