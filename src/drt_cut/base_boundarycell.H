/*---------------------------------------------------------------------*/
/*!
\file base_boundarycell.H

\brief used in direct divergence line integration

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef BASE_BOUNDARYCELL_H
#define BASE_BOUNDARYCELL_H

#include<vector>
#include<cmath>
#include<iostream>

#include "../linalg/linalg_fixedsizematrix.H"

/*!
\brief Returns the base function when the boundaycell is projected in y-z plane
*/
//pt(0,0) --> y
//pt(1,0) --> z
double base_func_surfX(const LINALG::Matrix<2,1>& pt, int inte_num,std::vector<double> alfa)
{
  double basef_surf=0.0;
  if(inte_num==1) //f(x,y,z) = 1
  {
    basef_surf = pt(0,0);
    return basef_surf;
  }
  double y=pt(0,0),y2=pow(pt(0,0),2),z=pt(1,0);
  double a0=alfa[0],a1=alfa[1],a2=alfa[2];
  if(inte_num==2) //f(x,y,z) = x
  {
    basef_surf = a2*y*z+0.5*a1*y2+a0*y;
    return basef_surf;
  }
  if(inte_num==3) //f(x,y,z) = y
  {
    basef_surf = 0.5*y*y;
    return basef_surf;
  }
  if(inte_num==4) //f(x,y,z) = z
  {
    basef_surf = z*y;
    return basef_surf;
  }

  double a22=pow(a2,2),a02=pow(a0,2),a12=pow(a1,2);
  double z2=pow(z,2),y3=pow(y,3);
  if(inte_num==5) //f(x,y,z) = x^2
  {
    basef_surf = a22*y*z2+2*a0*(a2*y*z+0.5*a1*y2)+a1*a2*y2*z+a12*y3/3.0+a02*y;
    return basef_surf;

  }
  if(inte_num==6) //f(x,y,z) = xy
  {
    basef_surf = y2*(3*a2*z+3*a0)+2*a1*y3;
    return basef_surf/6.0;
  }
  if(inte_num==7) //f(x,y,z) = xz
  {
    basef_surf = (a2*y*z+0.5*a1*y2+a0*y)*z;
    return basef_surf;
  }
  if(inte_num==8) //f(x,y,z) = y^2
  {
    basef_surf = y3/3.0;
    return basef_surf;
  }
  if(inte_num==9) //f(x,y,z) = yz
  {
    basef_surf = 0.5*y2*z;
    return basef_surf;
  }
  if(inte_num==10) //f(x,y,z) = z^2
  {
    basef_surf = z2*y;
    return basef_surf;
  }

  double a03=pow(a0,3),a13=pow(a1,3),a23=pow(a2,3);
  double y4=pow(y,4),z3=pow(z,3);
  if(inte_num==11) //f(x,y,z) = x^3
  {
    basef_surf = a23*y*z3+3.0*a0*(a22*y*z2+a1*a2*y2*z+a12*y3/3.0)+1.5*a1*a22*y2*z2+3.0*a02*(a2*y*z+0.5*a1*y2)+
      a12*a2*y3*z+0.25*a13*y4+a03*y;
    return basef_surf;
  }
  if(inte_num==12) //f(x,y,z) = x^2*y
  {
    basef_surf = y2*(6.0*a22*z2+12.0*a0*a2*z+6.0*a02)+y3*(8.0*a1*a2*z+8.0*a0*a1)+3.0*a12*y4;
    return basef_surf/12.0;
  }
  if(inte_num==13) //f(x,y,z) = x^2*z
  {
    basef_surf = z*(a22*y*z2+2.0*a0*(a2*y*z+0.5*a1*y2)+a1*a2*y2*z+a12*y3/3.0+a02*y);
    return basef_surf;
  }
  if(inte_num==14) //f(x,y,z) = xy^2
  {
    basef_surf = y3*(4.0*a2*z+4.0*a0)+3.0*a1*y4;
    return basef_surf/12.0;
  }
  if(inte_num==15) //f(x,y,z) = xyz
  {
    basef_surf = z*(y2*(3.0*a2*z+3.0*a0)+2.0*a1*y3);
    return basef_surf/6.0;
  }
  if(inte_num==16) //f(x,y,z) = xz^2
  {
    basef_surf = z2*(a2*y*z+0.5*a1*y2+a0*y);
    return basef_surf;
  }
  if(inte_num==17) //f(x,y,z) = y^3
  {
    basef_surf = 0.25*y4;
    return basef_surf;
  }
  if(inte_num==18) //f(x,y,z) = y^2*z
  {
    basef_surf = y3*z/3.0;
    return basef_surf;
  }
  if(inte_num==19) //f(x,y,z) = yz^2
  {
    basef_surf = y2*z2*0.5;
    return basef_surf;
  }
  if(inte_num==20) //f(x,y,z) = z^3
  {
    basef_surf = y*z3;
    return basef_surf;
  }

  double a04=pow(a0,4),a14=pow(a1,4),a24=pow(a2,4);
  double y5=pow(y,5),z4=pow(z,4);
  if(inte_num==21) //f(x,y,z) = x^4
  {
    basef_surf = a24*y*z4+4.0*a0*(a23*y*z3+1.5*a1*a22*y2*z2+a12*a2*y3*z+0.25*a13*y4)+2.0*a1*a23*y2*z3+
      6.0*a02*(a22*y*z2+a1*a2*y2*z+a12*y3/3.0)+2.0*a12*a22*y3*z2+4.0*a03*(a2*y*z+0.5*a1*y2)+
      a13*a2*y4*z+0.2*a14*y5+a04*y;
    return basef_surf;
  }
  if(inte_num==22) //f(x,y,z) = x^3*y
  {
    basef_surf = y2*(10.0*a23*z3+30.0*a0*a22*z2+30.0*a02*a2*z+10.0*a03)+y3*(20.0*a1*a22*z2+40.0*a0*a1*a2*z+
      20.0*a02*a1)+y4*(15.0*a12*a2*z+15.0*a0*a12)+4.0*a13*y5;
    return basef_surf*0.05;
  }
  if(inte_num==23) //f(x,y,z) = x^3*z
  {
    basef_surf = z*(a23*y*z3+3.0*a0*(a22*y*z2+a1*a2*y2*z+a12*y3/3.0)+1.5*a1*a22*y2*z2+3.0*a02*(a2*y*z+0.5*a1*y2)+
      a12*a2*y3*z+0.25*a13*y4+a03*y);
    return basef_surf;
  }
  if(inte_num==24) //f(x,y,z) = x^2*y^2
  {
    basef_surf = y3*(10.0*a22*z2+20.0*a0*a2*z+10.0*a02)+y4*(15.0*a1*a2*z+15.0*a0*a1)+6.0*a12*y5;
    return basef_surf/30.0;
  }
  if(inte_num==25) //f(x,y,z) = x^2*yz
  {
    basef_surf = z*(y2*(6.0*a22*z2+12.0*a0*a2*z+6.0*a02)+y3*(8.0*a1*a2*z+8.0*a0*a1)+3.0*a12*y4);
    return basef_surf/12.0;
  }
  if(inte_num==26) //f(x,y,z) = x^2*z^2
  {
    basef_surf = z2*(a22*y*z2+2.0*a0*(a2*y*z+0.5*a1*y2)+a1*a2*y2*z+a12*y3/3.0+a02*y);
    return basef_surf;
  }
  if(inte_num==27) //f(x,y,z) = x*y^3
  {
    basef_surf = y4*(5.0*a2*z+5.0*a0)+4.0*a1*y5;
    return basef_surf*0.05;
  }
  if(inte_num==28) //f(x,y,z) = x*y^2*z
  {
    basef_surf = z*(y3*(4.0*a2*z+4.0*a0)+3.0*a1*y4);
    return basef_surf/12.0;
  }
  if(inte_num==29) //f(x,y,z) = xyz^2
  {
    basef_surf = z2*(y2*(3.0*a2*z+3.0*a0)+2.0*a1*y3);
    return basef_surf/6.0;
  }
  if(inte_num==30) //f(x,y,z) = xz^3
  {
    basef_surf = z3*(a2*y*z+0.5*a1*y2+a0*y);
    return basef_surf;
  }
  if(inte_num==31) //f(x,y,z) = y^4
  {
    basef_surf = 0.2*y5;
    return basef_surf;
  }
  if(inte_num==32) //f(x,y,z) = y^3*z
  {
    basef_surf = 0.25*y4*z;
    return basef_surf;
  }
  if(inte_num==33) //f(x,y,z) = y^2*z^2
  {
    basef_surf = y3*z2/3.0;
    return basef_surf;
  }
  if(inte_num==34) //f(x,y,z) = yz^3
  {
    basef_surf = 0.5*y2*z3;
    return basef_surf;
  }
  if(inte_num==35) //f(x,y,z) = z^4
  {
    basef_surf = y*z4;
    return basef_surf;
  }

  double a05=pow(a0,5),a15=pow(a1,5),a25=pow(a2,5);
  double y6=pow(y,6),z5=pow(z,5);
  if(inte_num==36) //f(x,y,z) = x^5
  {
    basef_surf = a25*y*z5+5.0*a0*(a24*y*z4+2.0*a1*a23*y2*z3+2.0*a12*a22*y3*z2+a13*a2*y4*z+0.2*a14*y5)+2.5*a1*a24*y2*z4+
      10.0*a02*(a23*y*z3+1.5*a1*a22*y2*z2+a12*a2*y3*z+0.25*a13*y4)+10.0*a12*a23*y3*z3/3.0+10.0*a03*(a22*y*z2+
      a1*a2*y2*z+a12*y3/3.0)+2.5*a13*a22*y4*z2+5.0*a04*(a2*y*z+0.5*a1*y2)+a14*a2*y5*z+a15*y6/6.0+a05*y;
    return basef_surf;
  }
  if(inte_num==37) //f(x,y,z) = x^4*y
  {
    basef_surf = y2*(15.0*a24*z4+60.0*a0*a23*z3+90.0*a02*a22*z2+60.0*a03*a2*z+15.0*a04)+y3*(40.0*a1*a23*z3+120.0*a0*a1*a22*z2+
      120.0*a02*a1*a2*z+40.0*a03*a1)+y4*(45.0*a12*a22*z2+90.0*a0*a12*a2*z+45.0*a02*a12)+y5*(24.0*a13*a2*z+24.0*a0*a13)+
      5.0*a14*y6;
    return basef_surf/30.0;
  }
  if(inte_num==38) //f(x,y,z) = x^4*z
  {
    basef_surf = z*(a24*y*z4+4.0*a0*(a23*y*z3+1.5*a1*a22*y2*z2+a12*a2*y3*z+0.25*a13*y4)+2.0*a1*a23*y2*z3+6.0*a02*(a22*y*z2+
      a1*a2*y2*z+a12*y3/3.0)+2.0*a12*a22*y3*z2+4.0*a03*(a2*y*z+0.5*a1*y2)+a13*a2*y4*z+0.2*a14*y5+a04*y);
    return basef_surf;
  }
  if(inte_num==39) //f(x,y,z) = x^3*y^2
  {
    basef_surf = y3*(20.0*a23*z3+60.0*a0*a22*z2+60.0*a02*a2*z+20.0*a03)+y4*(45.0*a1*a22*z2+90.0*a0*a1*a2*z+45.0*a02*a1)+
      y5*(36.0*a12*a2*z+36.0*a0*a12)+10.0*a13*y6;
    return basef_surf/60.0;
  }
  if(inte_num==40) //f(x,y,z) = x^3*yz
  {
    basef_surf = z*(y2*(10.0*a23*z3+30.0*a0*a22*z2+30.0*a02*a2*z+10.0*a03)+y3*(20.0*a1*a22*z2+40.0*a0*a1*a2*z+20.0*a02*a1)+
      y4*(15.0*a12*a2*z+15.0*a0*a12)+4.0*a13*y5);
    return basef_surf*0.05;
  }
  if(inte_num==41) //f(x,y,z) = x^3*z^2
  {
    basef_surf = z2*(a23*y*z3+3.0*a0*(a22*y*z2+a1*a2*y2*z+a12*y3/3.0)+1.5*a1*a22*y2*z2+3.0*a02*(a2*y*z+0.5*a1*y2)+a12*a2*y3*z+
      0.25*a13*y4+a03*y);
    return basef_surf;
  }
  if(inte_num==42) //f(x,y,z) = x^2*y^3
  {
    basef_surf = y4*(15.0*a22*z2+30.0*a0*a2*z+15.0*a02)+y5*(24.0*a1*a2*z+24.0*a0*a1)+10.0*a12*y6;
    return basef_surf/60.0;
  }
  if(inte_num==43) //f(x,y,z) = x^2*y^2*z
  {
    basef_surf = z*(y3*(10.0*a22*z2+20.0*a0*a2*z+10.0*a02)+y4*(15.0*a1*a2*z+15.0*a0*a1)+6.0*a12*y5);
    return basef_surf/30.0;
  }
  if(inte_num==44) //f(x,y,z) = x^2*yz^2
  {
    basef_surf = z2*(y2*(6.0*a22*z2+12.0*a0*a2*z+6.0*a02)+y3*(8.0*a1*a2*z+8.0*a0*a1)+3.0*a12*y4);
    return basef_surf/12.0;
  }
  if(inte_num==45) //f(x,y,z) = x^2*z^3
  {
    basef_surf = z3*(a22*y*z2+2.0*a0*(a2*y*z+0.5*a1*y2)+a1*a2*y2*z+a12*y3/3.0+a02*y);
    return basef_surf;
  }
  if(inte_num==46) //f(x,y,z) = xy^4
  {
    basef_surf = y5*(6.0*a2*z+6.0*a0)+5.0*a1*y6;
    return basef_surf/30.0;
  }
  if(inte_num==47) //f(x,y,z) = xy^3*z
  {
    basef_surf = z*(y4*(5.0*a2*z+5.0*a0)+4.0*a1*y5);
    return basef_surf*0.05;
  }
  if(inte_num==48) //f(x,y,z) = xy^2*z^2
  {
    basef_surf = z2*(y3*(4.0*a2*z+4.0*a0)+3.0*a1*y4);
    return basef_surf/12.0;
  }
  if(inte_num==49) //f(x,y,z) = xyz^3
  {
    basef_surf = z3*(y2*(3.0*a2*z+3.0*a0)+2.0*a1*y3);
    return basef_surf/6.0;
  }
  if(inte_num==50) //f(x,y,z) = xz^4
  {
    basef_surf = z4*(a2*y*z+0.5*a1*y2+a0*y);
    return basef_surf;
  }
  if(inte_num==51) //f(x,y,z) = y^5
  {
    basef_surf = y6/6.0;
    return basef_surf;
  }
  if(inte_num==52) //f(x,y,z) = y^4*z
  {
    basef_surf = 0.2*y5*z;
    return basef_surf;
  }
  if(inte_num==53) //f(x,y,z) = y^3*z^2
  {
    basef_surf = 0.25*y4*z2;
    return basef_surf;
  }
  if(inte_num==54) //f(x,y,z) = y^2*z^3
  {
    basef_surf = y3*z3/3.0;
    return basef_surf;
  }
  if(inte_num==55) //f(x,y,z) = yz^4
  {
    basef_surf = y2*z4*0.5;
    return basef_surf;
  }
  if(inte_num==56) //f(x,y,z) = z^5
  {
    basef_surf = y*z5;
    return basef_surf;
  }

  double y7=pow(y,7),z6=pow(z,6);
  if(inte_num==57) //f(x,y,z) = x^6
  {
    if(fabs(a1)>0.00000001)
      basef_surf = std::pow((a2*z+a1*y+a0),7)/7.0/a1;
    else
      basef_surf = std::pow((a0+a2*z),6)*y;
    return basef_surf;
  }
  if(inte_num==58) //f(x,y,z) = x^5*y
  {
    basef_surf = y2*(21.0*a25*z5+105.0*a0*a24*z4+210.0*a02*a23*z3+210.0*a03*a22*z2+105.0*a04*a2*z+21.0*a05)+
      y3*(70.0*a1*a24*z4+280.0*a0*a1*a23*z3+420.0*a02*a1*a22*z2+280.0*a03*a1*a2*z+70.0*a04*a1)+
      y4*(105.0*a12*a23*z3+315.0*a0*a12*a22*z2+315.0*a02*a12*a2*z+105.0*a03*a12)+y5*(84.0*a13*a22*z2+
      168.0*a0*a13*a2*z+84.0*a02*a13)+y6*(35.0*a14*a2*z+35.0*a0*a14)+6.0*a15*y7;
    return basef_surf/42.0;
  }
  if(inte_num==59) //f(x,y,z) = x^5*z
  {
    basef_surf = z*(a25*y*z5+5.0*a0*(a24*y*z4+2.0*a1*a23*y2*z3+2.0*a12*a22*y3*z2+a13*a2*y4*z+0.2*a14*y5)+
      2.5*a1*a24*y2*z4+10.0*a02*(a23*y*z3+1.5*a1*a22*y2*z2+a12*a2*y3*z+0.25*a13*y4)+10.0/3.0*a12*a23*y3*z3+
      10.0*a03*(a22*y*z2+a1*a2*y2*z+a12*y3/3.0)+2.5*a13*a22*y4*z2+5.0*a04*(a2*y*z+0.5*a1*y2)+
      a14*a2*y5*z+a15*y6/6.0+a05*y);
    return basef_surf;
  }
  if(inte_num==60) //f(x,y,z) = x^4*y^2
  {
    basef_surf = y3*(35.0*a24*z4+140.0*a0*a23*z3+210.0*a02*a22*z2+140.0*a03*a2*z+35.0*a04)+y4*(105.0*a1*a23*z3+315.0*a0*a1*a22*z2+
      315.0*a02*a1*a2*z+105.0*a03*a1)+y5*(126.0*a12*a22*z2+252.0*a0*a12*a2*z+126.0*a02*a12)+y6*(70.0*a13*a2*z+70.0*a0*a13)+
      15.0*a14*y7;
    return basef_surf/105.0;
  }
  if(inte_num==61) //f(x,y,z) = x^4*yz
  {
    basef_surf = z*(y2*(15.0*a24*z4+60.0*a0*a23*z3+90.0*a02*a22*z2+60.0*a03*a2*z+15.0*a04)+y3*(40.0*a1*a23*z3+120.0*a0*a1*a22*z2+
      120.0*a02*a1*a2*z+40.0*a03*a1)+y4*(45.0*a12*a22*z2+90.0*a0*a12*a2*z+45.0*a02*a12)+y5*(24.0*a13*a2*z+24.0*a0*a13)+5.0*a14*y6);
    return basef_surf/30.0;
  }
  if(inte_num==62) //f(x,y,z) = x^4*z^2
  {
    basef_surf = z2*(a24*y*z4+4.0*a0*(a23*y*z3+1.5*a1*a22*y2*z2+a12*a2*y3*z+0.25*a13*y4)+2.0*a1*a23*y2*z3+6.0*a02*(a22*y*z2+
      a1*a2*y2*z+a12*y3/3.0)+2.0*a12*a22*y3*z2+4.0*a03*(a2*y*z+0.5*a1*y2)+a13*a2*y4*z+0.2*a14*y5+a04*y);
    return basef_surf;
  }
  if(inte_num==63) //f(x,y,z) = x^3*y^3
  {
    basef_surf = y4*(35.0*a23*z3+105.0*a0*a22*z2+105.0*a02*a2*z+35.0*a03)+y5*(84.0*a1*a22*z2+168.0*a0*a1*a2*z+84.0*a02*a1)+
      y6*(70.0*a12*a2*z+70.0*a0*a12)+20.0*a13*y7;
    return basef_surf/140.0;
  }
  if(inte_num==64) //f(x,y,z) = x^3*y^2*z
  {
    basef_surf = z*(y3*(20.0*a23*z3+60.0*a0*a22*z2+60.0*a02*a2*z+20.0*a03)+y4*(45.0*a1*a22*z2+90.0*a0*a1*a2*z+45.0*a02*a1)+
      y5*(36.0*a12*a2*z+36.0*a0*a12)+10.0*a13*y6);
    return basef_surf/60.0;
  }
  if(inte_num==65) //f(x,y,z) = x^3*y*z^2
  {
    basef_surf = z2*(y2*(10.0*a23*z3+30.0*a0*a22*z2+30.0*a02*a2*z+10.0*a03)+y3*(20.0*a1*a22*z2+40.0*a0*a1*a2*z+20.0*a02*a1)+
      y4*(15.0*a12*a2*z+15.0*a0*a12)+4.0*a13*y5);
    return basef_surf*0.05;
  }
  if(inte_num==66) //f(x,y,z) = x^3*z^3
  {
    basef_surf = z3*(a23*y*z3+3.0*a0*(a22*y*z2+a1*a2*y2*z+a12*y3/3.0)+1.5*a1*a22*y2*z2+3.0*a02*(a2*y*z+0.5*a1*y2)+
      a12*a2*y3*z+0.25*a13*y4+a03*y);
    return basef_surf;
  }
  if(inte_num==67) //f(x,y,z) = x^2*y^4
  {
    basef_surf = y5*(21.0*a22*z2+42.0*a0*a2*z+21.0*a02)+y6*(35.0*a1*a2*z+35.0*a0*a1)+15.0*a12*y7;
    return basef_surf/105.0;
  }
  if(inte_num==68) //f(x,y,z) = x^2*y^3*z
  {
    basef_surf = z*(y4*(15.0*a22*z2+30.0*a0*a2*z+15.0*a02)+y5*(24.0*a1*a2*z+24.0*a0*a1)+10.0*a12*y6);
    return basef_surf/60.0;
  }
  if(inte_num==69) //f(x,y,z) = x^2*y^2*z^2
  {
    basef_surf = z2*(y3*(10.0*a22*z2+20.0*a0*a2*z+10.0*a02)+y4*(15.0*a1*a2*z+15.0*a0*a1)+6.0*a12*y5);
    return basef_surf/30.0;
  }
  if(inte_num==70) //f(x,y,z) = x^2*yz^3
  {
    basef_surf = z3*(y2*(6.0*a22*z2+12.0*a0*a2*z+6.0*a02)+y3*(8.0*a1*a2*z+8.0*a0*a1)+3.0*a12*y4);
    return basef_surf/12.0;
  }
  if(inte_num==71) //f(x,y,z) = x^2*z^4
  {
    basef_surf = z4*(a22*y*z2+2.0*a0*(a2*y*z+0.5*a1*y2)+a1*a2*y2*z+a12*y3/3.0+a02*y);
    return basef_surf;
  }
  if(inte_num==72) //f(x,y,z) = xy^5
  {
    basef_surf = y6*(7.0*a2*z+7.0*a0)+6.0*a1*y7;
    return basef_surf/42.0;
  }
  if(inte_num==73) //f(x,y,z) = xy^4*z
  {
    basef_surf = z*(y5*(6.0*a2*z+6.0*a0)+5.0*a1*y6);
    return basef_surf/30.0;
  }
  if(inte_num==74) //f(x,y,z) = xy^3*z^2
  {
    basef_surf = z2*(y4*(5.0*a2*z+5.0*a0)+4.0*a1*y5);
    return basef_surf*0.05;
  }
  if(inte_num==75) //f(x,y,z) = xy^2*z^3
  {
    basef_surf = z3*(y3*(4.0*a2*z+4.0*a0)+3.0*a1*y4);
    return basef_surf/12.0;
  }
  if(inte_num==76) //f(x,y,z) = xyz^4
  {
    basef_surf = z4*(y2*(3.0*a2*z+3.0*a0)+2.0*a1*y3);
    return basef_surf/6.0;
  }
  if(inte_num==77) //f(x,y,z) = xz^5
  {
    basef_surf = z5*(a2*y*z+0.5*a1*y2+a0*y);
    return basef_surf;
  }
  if(inte_num==78) //f(x,y,z) = y^6
  {
    basef_surf = y7/7.0;
    return basef_surf;
  }
  if(inte_num==79) //f(x,y,z) = y^5*z
  {
    basef_surf = y6*z/6.0;
    return basef_surf;
  }
  if(inte_num==80) //f(x,y,z) = y^4*z^2
  {
    basef_surf = 0.2*y5*z2;
    return basef_surf;
  }
  if(inte_num==81) //f(x,y,z) = y^3*z^3
  {
    basef_surf = 0.25*y4*z3;
    return basef_surf;
  }
  if(inte_num==82) //f(x,y,z) = y^2*z^4
  {
    basef_surf = y3*z4/3.0;
    return basef_surf;
  }
  if(inte_num==83) //f(x,y,z) = yz^5
  {
    basef_surf = 0.5*y2*z5;
    return basef_surf;
  }
  if(inte_num==84) //f(x,y,z) = z^6
  {
    basef_surf = y*z6;
    return basef_surf;
  }

  dserror("The base function for boundarycell integration undefined");
  exit(1);
  return 0.0;

}

/*!
\brief Returns the base function when the boundaycell is projected in z-x plane
*/
//pt(0,0) --> z
//pt(1,0) --> x
double base_func_surfY(const LINALG::Matrix<2,1>& pt, int inte_num,std::vector<double> alfa)
{
  double basef_surf=0.0;

  if(inte_num==1) //f(x,y,z) = 1
  {
    basef_surf = pt(0,0);
    return basef_surf;
  }

  double a0=alfa[0],a1=alfa[1],a2=alfa[2];
  double z=pt(0,0),x=pt(1,0),z2=pow(pt(0,0),2);
  if(inte_num==2) //f(x,y,z) = x
  {
    basef_surf = x*z;
    return basef_surf;
  }
  if(inte_num==3) //f(x,y,z) = y
  {
    basef_surf = 0.5*a1*z2+a2*x*z+a0*z;
    return basef_surf;
  }
  if(inte_num==4) //f(x,y,z) = z
  {
    basef_surf = 0.5*z*z;
    return basef_surf;
  }

  double a22=pow(a2,2),a02=pow(a0,2),a12=pow(a1,2);
  double x2=pow(x,2),z3=pow(z,3);
  if(inte_num==5) //f(x,y,z) = x^2
  {
    basef_surf = x2*z;
    return basef_surf;
  }
  if(inte_num==6) //f(x,y,z) = xy
  {
    basef_surf = x*(a1*z2*0.5+a2*x*z+a0*z);
    return basef_surf;
  }
  if(inte_num==7) //f(x,y,z) = xz
  {
    basef_surf = 0.5*x*z2;
    return basef_surf;
  }
  if(inte_num==8) //f(x,y,z) = y^2
  {
    basef_surf = a12*z3/3.0+2.0*a0*(0.5*a1*z2+a2*x*z)+a1*a2*x*z2+a22*x2*z+a02*z;
    return basef_surf;
  }
  if(inte_num==9) //f(x,y,z) = yz
  {
    basef_surf = 2.0*a1*z3+(3.0*a2*x+3.0*a0)*z2;
    return basef_surf/6.0;
  }
  if(inte_num==10) //f(x,y,z) = z^2
  {
    basef_surf = z3/3.0;
    return basef_surf;
  }

  double a03=pow(a0,3),a13=pow(a1,3),a23=pow(a2,3);
  double x3=pow(x,3),z4=pow(z,4);
  if(inte_num==11) //f(x,y,z) = x^3
  {
    basef_surf = x3*z;
    return basef_surf;
  }
  if(inte_num==12) //f(x,y,z) = x^2*y
  {
    basef_surf = x2*(0.5*a1*z2+a2*x*z+a0*z);
    return basef_surf;
  }
  if(inte_num==13) //f(x,y,z) = x^2*z
  {
    basef_surf = x2*z2*0.5;
    return basef_surf;
  }
  if(inte_num==14) //f(x,y,z) = xy^2
  {
    basef_surf = x*(a12*z3/3.0+2.0*a0*(a1*z2*0.5+a2*x*z)+a1*a2*x*z2+a22*x2*z+a02*z);
    return basef_surf;
  }
  if(inte_num==15) //f(x,y,z) = xyz
  {
    basef_surf = x*(2.0*a1*z3+(3.0*a2*x+3.0*a0)*z2);
    return basef_surf/6.0;
  }
  if(inte_num==16) //f(x,y,z) = xz^2
  {
    basef_surf = x*z3/3.0;
    return basef_surf;
  }
  if(inte_num==17) //f(x,y,z) = y^3
  {
    basef_surf = 0.25*a13*z4+3.0*a0*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+a12*a2*x*z3+3.0*a02*(0.5*a1*z2+a2*x*z)+
      1.5*a1*a22*x2*z2+a23*x3*z+a03*z;
    return basef_surf;
  }
  if(inte_num==18) //f(x,y,z) = y^2*z
  {
    basef_surf = 3.0*a12*z4+(8.0*a1*a2*x+8.0*a0*a1)*z3+(6.0*a22*x2+12.0*a0*a2*x+6.0*a02)*z2;
    return basef_surf/12.0;
  }
  if(inte_num==19) //f(x,y,z) = yz^2
  {
    basef_surf = 3.0*a1*z4+(4.0*a2*x+4.0*a0)*z3;
    return basef_surf/12.0;
  }
  if(inte_num==20) //f(x,y,z) = z^3
  {
    basef_surf = 0.25*z4;
    return basef_surf;
  }

  double a04=pow(a0,4),a14=pow(a1,4),a24=pow(a2,4);
  double x4=pow(x,4),z5=pow(z,5);
  if(inte_num==21) //f(x,y,z) = x^4
  {
    basef_surf = x4*z;
    return basef_surf;
  }
  if(inte_num==22) //f(x,y,z) = x^3*y
  {
    basef_surf = x3*(0.5*a1*z2+a2*x*z+a0*z);
    return basef_surf;
  }
  if(inte_num==23) //f(x,y,z) = x^3*z
  {
    basef_surf = x3*z2*0.5;
    return basef_surf;
  }
  if(inte_num==24) //f(x,y,z) x^2*y^2
  {
    basef_surf = x2*(a12*z3/3.0+2.0*a0*(0.5*a1*z2+a2*x*z)+a1*a2*x*z2+a22*x2*z+a02*z);
    return basef_surf;
  }
  if(inte_num==25) //f(x,y,z) = x^2*yz
  {
    basef_surf = x2*(2.0*a1*z3+(3.0*a2*x+3.0*a0)*z2);
    return basef_surf/6.0;
  }
  if(inte_num==26) //f(x,y,z) = x^2*z^2
  {
    basef_surf = x2*z3/3.0;
    return basef_surf;
  }
  if(inte_num==27) //f(x,y,z) = xy^3
  {
    basef_surf = x*(a13*z4*0.25+3.0*a0*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+a12*a2*x*z3+3.0*a02*(0.5*a1*z2+a2*x*z)+
      1.5*a1*a22*x2*z2+a23*x3*z+a03*z);
    return basef_surf;
  }
  if(inte_num==28) //f(x,y,z) = xy^2*z
  {
    basef_surf = x*(3.0*a12*z4+(8.0*a1*a2*x+8.0*a0*a1)*z3+(6.0*a22*x2+12.0*a0*a2*x+6.0*a02)*z2);
    return basef_surf/12.0;
  }
  if(inte_num==29) //f(x,y,z) = xyz^2
  {
    basef_surf = x*(3.0*a1*z4+(4.0*a2*x+4.0*a0)*z3);
    return basef_surf/12.0;
  }
  if(inte_num==30) //f(x,y,z) = xz^3
  {
    basef_surf = 0.25*x*z4;
    return basef_surf;
  }
  if(inte_num==31) //f(x,y,z) = y^4
  {
    basef_surf = 0.2*a14*z5+4.0*a0*(0.25*a13*z4+a12*a2*x*z3+1.5*a1*a22*x2*z2+a23*x3*z)+a13*a2*x*z4+6.0*a02*(a12*z3/3.0+
      a1*a2*x*z2+a22*x2*z)+2.0*a12*a22*x2*z3+4.0*a03*(0.5*a1*z2+a2*x*z)+2.0*a1*a23*x3*z2+a24*x4*z+a04*z;
    return basef_surf;
  }
  if(inte_num==32) //f(x,y,z) = y^3*z
  {
    basef_surf = 4.0*a13*z5+(15.0*a12*a2*x+15.0*a0*a12)*z4+(20.0*a1*a22*x2+40.0*a0*a1*a2*x+20.0*a02*a1)*z3+(10.0*a23*x3+
      30.0*a0*a22*x2+30.0*a02*a2*x+10.0*a03)*z2;
    return basef_surf*0.05;
  }
  if(inte_num==33) //f(x,y,z) = y^2*z^2
  {
    basef_surf = 6.0*a12*z5+(15.0*a1*a2*x+15.0*a0*a1)*z4+(10.0*a22*x2+20.0*a0*a2*x+10.0*a02)*z3;
    return basef_surf/30.0;
  }
  if(inte_num==34) //f(x,y,z) = yz^3
  {
    basef_surf = 4.0*a1*z5+(5.0*a2*x+5.0*a0)*z4;
    return basef_surf*0.05;
  }
  if(inte_num==35) //f(x,y,z) = z^4
  {
    basef_surf = 0.2*z5;
    return basef_surf;
  }

  double a05=pow(a0,5),a15=pow(a1,5),a25=pow(a2,5);
  double x5=pow(x,5),z6=pow(z,6);
  if(inte_num==36)  //f(x,y,z) = x^5
  {
    basef_surf = x5*z;
    return basef_surf;
  }
  if(inte_num==37) //f(x,y,z) = x^4*y
  {
    basef_surf = x4*(0.5*a1*z2+a2*x*z+a0*z);
    return basef_surf;
  }
  if(inte_num==38) //f(x,y,z) = x^4*z
  {
    basef_surf = x4*z2*0.5;
    return basef_surf;
  }
  if(inte_num==39) //f(x,y,z) = x^3*y^2
  {
    basef_surf = x3*(a12*z3/3.0+2.0*a0*(0.5*a1*z2+a2*x*z)+a1*a2*x*z2+a22*x2*z+a02*z);
    return basef_surf;
  }
  if(inte_num==40) //f(x,y,z) = x^3*yz
  {
    basef_surf = x3*(2.0*a1*z3+(3.0*a2*x+3.0*a0)*z2);
    return basef_surf/6.0;
  }
  if(inte_num==41) //f(x,y,z) = x^3*z^2
  {
    basef_surf = x3*z3/3.0;
    return basef_surf;
  }
  if(inte_num==42) //f(x,y,z) = x^2*y^3
  {
    basef_surf = x2*(0.25*a13*z4+3.0*a0*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+a12*a2*x*z3+3.0*a02*(0.5*a1*z2+a2*x*z)+
      1.5*a1*a22*x2*z2+a23*x3*z+a03*z);
    return basef_surf;
  }
  if(inte_num==43) //f(x,y,z) = x^2*y^2*z
  {
    basef_surf = x2*(3.0*a12*z4+(8.0*a1*a2*x+8.0*a0*a1)*z3+(6.0*a22*x2+12.0*a0*a2*x+6.0*a02)*z2);
    return basef_surf/12.0;
  }
  if(inte_num==44) //f(x,y,z) = x^2*yz^2
  {
    basef_surf = x2*(3.0*a1*z4+(4.0*a2*x+4.0*a0)*z3);
    return basef_surf/12.0;
  }
  if(inte_num==45) //f(x,y,z) = x^2*z^3
  {
    basef_surf = x2*z4*0.25;
    return basef_surf;
  }
  if(inte_num==46) //f(x,y,z) = xy^4
  {
    basef_surf = x*(0.2*a14*z5+4.0*a0*(0.25*a13*z4+a12*a2*x*z3+1.5*a1*a22*x2*z2+a23*x3*z)+a13*a2*x*z4+
      6.0*a02*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+2.0*a12*a22*x2*z3+4.0*a03*(0.5*a1*z2+a2*x*z)+
      2.0*a1*a23*x3*z2+a24*x4*z+a04*z);
    return basef_surf;
  }
  if(inte_num==47) //f(x,y,z) = xy^3*z
  {
    basef_surf = x*(4.0*a13*z5+(15.0*a12*a2*x+15.0*a0*a12)*z4+(20.0*a1*a22*x2+40.0*a0*a1*a2*x+20.0*a02*a1)*z3+
      (10.0*a23*x3+30.0*a0*a22*x2+30.0*a02*a2*x+10.0*a03)*z2);
    return basef_surf*0.05;
  }
  if(inte_num==48) //f(x,y,z) = xy^2*z^2
  {
    basef_surf = x*(6.0*a12*z5+(15.0*a1*a2*x+15.0*a0*a1)*z4+(10.0*a22*x2+20.0*a0*a2*x+10.0*a02)*z3);
    return basef_surf/30.0;
  }
  if(inte_num==49) //f(x,y,z) = xyz^3
  {
    basef_surf = x*(4.0*a1*z5+(5.0*a2*x+5.0*a0)*z4);
    return basef_surf*0.05;
  }
  if(inte_num==50) //f(x,y,z) = xz^4
  {
    basef_surf = x*z5*0.2;
    return basef_surf;
  }
  if(inte_num==51) //f(x,y,z) = y^5
  {
    basef_surf = a15*z6/6.0+5.0*a0*(0.2*a14*z5+a13*a2*x*z4+2.0*a12*a22*x2*z3+2.0*a1*a23*x3*z2+a24*x4*z)+
      a14*a2*x*z5+10.0*a02*(0.25*a13*z4+a12*a2*x*z3+1.5*a1*a22*x2*z2+a23*x3*z)+2.5*a13*a22*x2*z4+
      10.0*a03*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+10.0/3.0*a12*a23*x3*z3+5.0*a04*(0.5*a1*z2+a2*x*z)+
      2.5*a1*a24*x4*z2+a25*x5*z+a05*z;
    return basef_surf;
  }
  if(inte_num==52) //f(x,y,z) = y^4*z
  {
    basef_surf = 5.0*a14*z6+(24.0*a13*a2*x+24.0*a0*a13)*z5+(45.0*a12*a22*x2+90.0*a0*a12*a2*x+45.0*a02*a12)*z4+
      (40.0*a1*a23*x3+120.0*a0*a1*a22*x2+120.0*a02*a1*a2*x+40.0*a03*a1)*z3+(15.0*a24*x4+60.0*a0*a23*x3+
      90.0*a02*a22*x2+60.0*a03*a2*x+15.0*a04)*z2;
    return basef_surf/30.0;
  }
  if(inte_num==53) //f(x,y,z) = y^3*z^2
  {
    basef_surf = 10.0*a13*z6+(36.0*a12*a2*x+36.0*a0*a12)*z5+(45.0*a1*a22*x2+90.0*a0*a1*a2*x+45.0*a02*a1)*z4+
      (20.0*a23*x3+60.0*a0*a22*x2+60.0*a02*a2*x+20.0*a03)*z3;
    return basef_surf/60.0;
  }
  if(inte_num==54) //f(x,y,z) = y^2*z^3
  {
    basef_surf = 10.0*a12*z6+(24.0*a1*a2*x+24.0*a0*a1)*z5+(15.0*a22*x2+30.0*a0*a2*x+15.0*a02)*z4;
    return basef_surf/60.0;
  }
  if(inte_num==55) //f(x,y,z) = yz^4
  {
    basef_surf = 5.0*a1*z6+(6.0*a2*x+6.0*a0)*z5;
    return basef_surf/30.0;
  }
  if(inte_num==56) //f(x,y,z) = z^5
  {
    basef_surf = z6/6.0;
    return basef_surf;
  }

  double x6=pow(x,6),z7=pow(z,7);
  if(inte_num==57) //f(x,y,z) = x^6
  {
    basef_surf = x6*z;
    return basef_surf;
  }
  if(inte_num==58) //f(x,y,z) = x^5*y
  {
    basef_surf = x5*(0.5*a1*z2+a2*x*z+a0*z);
    return basef_surf;
  }
  if(inte_num==59) //f(x,y,z) = x^5*z
  {
    basef_surf = x5*z2*0.5;
    return basef_surf;
  }
  if(inte_num==60) //f(x,y,z) = x^4*y^2
  {
    basef_surf = x4*(a12*z3/3.0+2.0*a0*(0.5*a1*z2+a2*x*z)+a1*a2*x*z2+a22*x2*z+a02*z);
    return basef_surf;
  }
  if(inte_num==61) //f(x,y,z) = x^4*yz
  {
    basef_surf = x4*(2.0*a1*z3+(3.0*a2*x+3.0*a0)*z2);
    return basef_surf/6.0;
  }
  if(inte_num==62) //f(x,y,z) = x^4*z^2
  {
    basef_surf = x4*z3/3.0;
    return basef_surf;
  }
  if(inte_num==63) //f(x,y,z) = x^3*y^3
  {
    basef_surf = x3*(0.25*a13*z4+3.0*a0*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+a12*a2*x*z3+3.0*a02*(0.5*a1*z2+a2*x*z)+
      1.5*a1*a22*x2*z2+a23*x3*z+a03*z);
    return basef_surf;
  }
  if(inte_num==64) //f(x,y,z) = x^3*y^2*z
  {
    basef_surf = x3*(3.0*a12*z4+(8.0*a1*a2*x+8.0*a0*a1)*z3+(6.0*a22*x2+12.0*a0*a2*x+6.0*a02)*z2);
    return basef_surf/12.0;
  }
  if(inte_num==65) //f(x,y,z) = x^3*y*z^2
  {
    basef_surf = x3*(3.0*a1*z4+(4.0*a2*x+4.0*a0)*z3);
    return basef_surf/12.0;
  }
  if(inte_num==66) //f(x,y,z) = x^3*z^3
  {
    basef_surf = x3*z4*0.25;
    return basef_surf;
  }
  if(inte_num==67) //f(x,y,z) = x^2*y^4
  {
    basef_surf = x2*(0.2*a14*z5+4.0*a0*(0.25*a13*z4+a12*a2*x*z3+1.5*a1*a22*x2*z2+a23*x3*z)+a13*a2*x*z4+
      6.0*a02*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+2.0*a12*a22*x2*z3+4.0*a03*(0.5*a1*z2+a2*x*z)+
      2.0*a1*a23*x3*z2+a24*x4*z+a04*z);
    return basef_surf;
  }
  if(inte_num==68) //f(x,y,z) = x^2*y^3*z
  {
    basef_surf = x2*(4.0*a13*z5+(15.0*a12*a2*x+15.0*a0*a12)*z4+(20.0*a1*a22*x2+40.0*a0*a1*a2*x+20.0*a02*a1)*z3+
      (10.0*a23*x3+30.0*a0*a22*x2+30.0*a02*a2*x+10.0*a03)*z2);
    return basef_surf*0.05;
  }
  if(inte_num==69) //f(x,y,z) = x^2*y^2*z^2
  {
    basef_surf = x2*(6.0*a12*z5+(15.0*a1*a2*x+15.0*a0*a1)*z4+(10.0*a22*x2+20.0*a0*a2*x+10.0*a02)*z3);
    return basef_surf/30.0;
  }
  if(inte_num==70) //f(x,y,z) = x^2*yz^3
  {
    basef_surf = x2*(4.0*a1*z5+(5.0*a2*x+5.0*a0)*z4);
    return basef_surf*0.05;
  }
  if(inte_num==71) //f(x,y,z) = x^2*z^4
  {
    basef_surf = x2*z5*0.2;
    return basef_surf;
  }
  if(inte_num==72) //f(x,y,z) = xy^5
  {
    basef_surf = x*(a15*z6/6.0+5.0*a0*(a14*z5*0.2+a13*a2*x*z4+2.0*a12*a22*x2*z3+2.0*a1*a23*x3*z2+
      a24*x4*z)+a14*a2*x*z5+10.0*a02*(0.25*a13*z4+a12*a2*x*z3+1.5*a1*a22*x2*z2+a23*x3*z)+
      2.5*a13*a22*x2*z4+10.0*a03*(a12*z3/3.0+a1*a2*x*z2+a22*x2*z)+10.0/3.0*a12*a23*x3*z3+
      5.0*a04*(0.5*a1*z2+a2*x*z)+2.5*a1*a24*x4*z2+a25*x5*z+a05*z);
    return basef_surf;
  }
  if(inte_num==73) //f(x,y,z) = xy^4*z
  {
    basef_surf = x*(5.0*a14*z6+(24.0*a13*a2*x+24.0*a0*a13)*z5+(45.0*a12*a22*x2+90.0*a0*a12*a2*x+45.0*a02*a12)*z4+
      (40.0*a1*a23*x3+120.0*a0*a1*a22*x2+120.0*a02*a1*a2*x+40.0*a03*a1)*z3+(15.0*a24*x4+60.0*a0*a23*x3+
      90.0*a02*a22*x2+60.0*a03*a2*x+15.0*a04)*z2);
    return basef_surf/30.0;
  }
  if(inte_num==74) //f(x,y,z) = xy^3*z^2
  {
    basef_surf = x*(10.0*a13*z6+(36.0*a12*a2*x+36.0*a0*a12)*z5+(45.0*a1*a22*x2+90.0*a0*a1*a2*x+45.0*a02*a1)*z4+
      (20.0*a23*x3+60.0*a0*a22*x2+60.0*a02*a2*x+20.0*a03)*z3);
    return basef_surf/60.0;
  }
  if(inte_num==75) //f(x,y,z) = xy^2*z^3
  {
    basef_surf = x*(10.0*a12*z6+(24.0*a1*a2*x+24.0*a0*a1)*z5+(15.0*a22*x2+30.0*a0*a2*x+15.0*a02)*z4);
    return basef_surf/60.0;
  }
  if(inte_num==76) //f(x,y,z) = xyz^4
  {
    basef_surf = x*(5.0*a1*z6+(6.0*a2*x+6.0*a0)*z5);
    return basef_surf/30.0;
  }
  if(inte_num==77) //f(x,y,z) = xz^5
  {
    basef_surf = x*z6/6.0;
    return basef_surf;
  }
  if(inte_num==78) //f(x,y,z) = y^6
  {
    if(fabs(a1)>0.00000001)
      basef_surf = std::pow((a0+a2*x+a1*z),7)/7.0/a1;
    else
      basef_surf = std::pow((a0+a1*x),6)*z;
    return basef_surf;
  }
  if(inte_num==79) //f(x,y,z) = y^5*z
  {
    basef_surf = 6.0*a15*z7+(35.0*a14*a2*x+35.0*a0*a14)*z6+(84.0*a13*a22*x2+168.0*a0*a13*a2*x+84.0*a02*a13)*z5+
      (105.0*a12*a23*x3+315.0*a0*a12*a22*x2+315.0*a02*a12*a2*x+105.0*a03*a12)*z4+(70.0*a1*a24*x4+280.0*a0*a1*a23*x3+
      420.0*a02*a1*a22*x2+280.0*a03*a1*a2*x+70.0*a04*a1)*z3+(21.0*a25*x5+105.0*a0*a24*x4+210.0*a02*a23*x3+
      210.0*a03*a22*x2+105.0*a04*a2*x+21.0*a05)*z2;
    return basef_surf/42.0;
  }
  if(inte_num==80) //f(x,y,z) = y^4*z^2
  {
    basef_surf = 15.0*a14*z7+(70.0*a13*a2*x+70.0*a0*a13)*z6+(126.0*a12*a22*x2+252.0*a0*a12*a2*x+126.0*a02*a12)*z5+
      (105.0*a1*a23*x3+315.0*a0*a1*a22*x2+315.0*a02*a1*a2*x+105.0*a03*a1)*z4+(35.0*a24*x4+140.0*a0*a23*x3+
      210.0*a02*a22*x2+140.0*a03*a2*x+35.0*a04)*z3;
    return basef_surf/105.0;
  }
  if(inte_num==81) //f(x,y,z) = y^3*z^3
  {
    basef_surf = 20.0*a13*z7+(70.0*a12*a2*x+70.0*a0*a12)*z6+(84.0*a1*a22*x2+168.0*a0*a1*a2*x+84.0*a02*a1)*z5+
      (35.0*a23*x3+105.0*a0*a22*x2+105.0*a02*a2*x+35.0*a03)*z4;
    return basef_surf/140.0;
  }
  if(inte_num==82) //f(x,y,z) = y^2*z^4
  {
    basef_surf = 15.0*a12*z7+(35.0*a1*a2*x+35.0*a0*a1)*z6+(21.0*a22*x2+42.0*a0*a2*x+21.0*a02)*z5;
    return basef_surf/105.0;
  }
  if(inte_num==83) //f(x,y,z) = yz^5
  {
    basef_surf = 6.0*a1*z7+(7.0*a2*x+7.0*a0)*z6;
    return basef_surf/42.0;
  }
  if(inte_num==84) //f(x,y,z) = z^6
  {
    basef_surf = z7/7.0;
    return basef_surf;
  }
  dserror("The base function for boundarycell integration undefined");
  exit(1);
  return 0.0;

}

/*!
\brief Returns the base function when the boundaycell is projected in x-y plane
*/
//pt(0,0) --> x
//pt(1,0) --> y
double base_func_surfZ(const LINALG::Matrix<2,1>& pt, int inte_num,std::vector<double> alfa)
{
  double basef_surf=0.0;
  if(inte_num==1) //f(x,y,z) = 1
  {
    basef_surf = pt(0,0);
    return basef_surf;
  }

  double a0=alfa[0],a1=alfa[1],a2=alfa[2];
  double x=pt(0,0),y=pt(1,0),x2=pow(pt(0,0),2);
  if(inte_num==2) //f(x,y,z) = x
  {
    basef_surf = 0.5*x*x;
    return basef_surf;
  }
  if(inte_num==3) //f(x,y,z) = y
  {
    basef_surf = y*x;
    return basef_surf;
  }
  if(inte_num==4) //f(x,y,z) = z
  {
    basef_surf = a2*x*y+a1*x2*0.5+a0*x;
    return basef_surf;
  }

  double a22=pow(a2,2),a02=pow(a0,2),a12=pow(a1,2);
  double x3=pow(x,3),y2=pow(y,2);
  if(inte_num==5) //f(x,y,z) = x^2
  {
    basef_surf = x3/3.0;
    return basef_surf;
  }
  if(inte_num==6) //f(x,y,z) = xy
  {
    basef_surf = 0.5*x2*y;
    return basef_surf;
  }
  if(inte_num==7) //f(x,y,z) = xz
  {
    basef_surf = x2*(3*a2*y+3*a0)+2*a1*x3;
    return basef_surf/6.0;
  }
  if(inte_num==8) //f(x,y,z) = y^2
  {
    basef_surf = x*y2;
    return basef_surf;
  }
  if(inte_num==9) //f(x,y,z) = yz
  {
    basef_surf = y*(a2*x*y+0.5*a1*x2+a0*x);
    return basef_surf;
  }
  if(inte_num==10) //f(x,y,z) = z^2
  {
    basef_surf = a22*x*y2+2.0*a0*(a2*x*y+0.5*a1*x2)+a1*a2*x2*y+a12*x3/3.0+a02*x;
    return basef_surf;
  }

  double x4=pow(x,4),y3=pow(y,3);
  double a03=pow(a0,3),a13=pow(a1,3),a23=pow(a2,3);
  if(inte_num==11) //f(x,y,z) = x^3
  {
    basef_surf = 0.25*x4;
    return basef_surf;
  }
  if(inte_num==12) //f(x,y,z) = x^2*y
  {
    basef_surf = x3*y/3.0;
    return basef_surf;
  }
  if(inte_num==13) //f(x,y,z) = x^2*z
  {
    basef_surf = x3*(4.0*a2*y+4.0*a0)+3.0*a1*x4;
    return basef_surf/12.0;
  }
  if(inte_num==14) //f(x,y,z) = xy^2
  {
    basef_surf = 0.5*x2*y2;
    return basef_surf;
  }
  if(inte_num==15) //f(x,y,z) = xyz
  {
    basef_surf = y*(x2*(3.0*a2*y+3.0*a0)+2.0*a1*x3);
    return basef_surf/6.0;
  }
  if(inte_num==16) //f(x,y,z) = xz^2
  {
    basef_surf = x2*(6.0*a22*y2+12.0*a0*a2*y+6.0*a02)+x3*(8.0*a1*a2*y+8.0*a0*a1)+3.0*a12*x4;
    return basef_surf/12.0;
  }
  if(inte_num==17) //f(x,y,z) = y^3
  {
    basef_surf = y3*x;
    return basef_surf;
  }
  if(inte_num==18) //f(x,y,z) = y^2*z
  {
    basef_surf = y2*(a2*x*y+0.5*a1*x2+a0*x);
    return basef_surf;
  }
  if(inte_num==19) //f(x,y,z) = yz^2
  {
    basef_surf = y*(a22*x*y2+2.0*a0*(a2*x*y+0.5*a1*x2)+a1*a2*x2*y+a12*x3/3.0+a02*x);
    return basef_surf;
  }
  if(inte_num==20) //f(x,y,z) = z^3
  {
    basef_surf = a23*x*y3+3.0*a0*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+1.5*a1*a22*x2*y2+3.0*a02*(a2*x*y+0.5*a1*x2)+
      a12*a2*x3*y+0.25*a13*x4+a03*x;
    return basef_surf;
  }

  double x5=pow(x,5),y4=pow(y,4);
  double a04=pow(a0,4),a14=pow(a1,4),a24=pow(a2,4);
  if(inte_num==21) //f(x,y,z) = x^4
  {
    basef_surf = 0.2*x5;
    return basef_surf;
  }
  if(inte_num==22) //f(x,y,z) = x^3*y
  {
    basef_surf = 0.25*x4*y;
    return basef_surf;
  }
  if(inte_num==23) //f(x,y,z) = x^3*z
  {
    basef_surf = x4*(5.0*a2*y+5.0*a0)+4.0*a1*x5;
    return basef_surf*0.05;
  }
  if(inte_num==24) //f(x,y,z) = x^2*y^2
  {
    basef_surf = x3*y2/3.0;
    return basef_surf;
  }
  if(inte_num==25) //f(x,y,z) = x^2*yz
  {
    basef_surf = y*(x3*(4.0*a2*y+4.0*a0)+3.0*a1*x4);
    return basef_surf/12.0;
  }
  if(inte_num==26) //f(x,y,z) = x^2*z^2
  {
    basef_surf = x3*(10.0*a22*y2+20.0*a0*a2*y+10.0*a02)+x4*(15.0*a1*a2*y+15.0*a0*a1)+6.0*a12*x5;
    return basef_surf/30.0;
  }
  if(inte_num==27) //f(x,y,z) = x*y^3
  {
    basef_surf = x2*y3*0.5;
    return basef_surf;
  }
  if(inte_num==28) //f(x,y,z) = xy^2*z
  {
    basef_surf = y2*(x2*(3.0*a2*y+3.0*a0)+2.0*a1*x3);
    return basef_surf/6.0;
  }
  if(inte_num==29) //f(x,y,z) = xyz^2
  {
    basef_surf = y*(x2*(6.0*a22*y2+12.0*a0*a2*y+6.0*a02)+x3*(8.0*a1*a2*y+8.0*a0*a1)+3.0*a12*x4);
    return basef_surf/12.0;
  }
  if(inte_num==30) //f(x,y,z) = xz^3
  {
    basef_surf = x2*(10.0*a23*y3+30.0*a0*a22*y2+30.0*a02*a2*y+10.0*a03)+x3*(20.0*a1*a22*y2+40.0*a0*a1*a2*y+
      20.0*a02*a1)+x4*(15.0*a12*a2*y+15.0*a0*a12)+4.0*a13*x5;
    return basef_surf*0.05;
  }
  if(inte_num==31) //f(x,y,z) = y^4
  {
    basef_surf = x*y4;
    return basef_surf;
  }
  if(inte_num==32) //f(x,y,z) = y^3*z
  {
    basef_surf = y3*(a2*x*y+0.5*a1*x2+a0*x);
    return basef_surf;
  }
  if(inte_num==33) //f(x,y,z) = y^2*z^2
  {
    basef_surf = y2*(a22*x*y2+2.0*a0*(a2*x*y+0.5*a1*x2)+a1*a2*x2*y+a12*x3/3.0+a02*x);
    return basef_surf;
  }
  if(inte_num==34) //f(x,y,z) = yz^3
  {
    basef_surf = y*(a23*x*y3+3.0*a0*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+1.5*a1*a22*x2*y2+3.0*a02*(a2*x*y+0.5*a1*x2)+
      a12*a2*x3*y+0.25*a13*x4+a03*x);
    return basef_surf;
  }
  if(inte_num==35) //f(x,y,z) = z^4
  {
    basef_surf = a24*x*y4+4.0*a0*(a23*x*y3+1.5*a1*a22*x2*y2+a12*a2*x3*y+0.25*a13*x4)+2.0*a1*a23*x2*y3+6.0*a02*(a22*x*y2+
      a1*a2*x2*y+a12*x3/3.0)+2.0*a12*a22*x3*y2+4.0*a03*(a2*x*y+0.5*a1*x2)+a13*a2*x4*y+0.2*a14*x5+a04*x;
    return basef_surf;
  }

  double x6=pow(x,6),y5=pow(y,5);
  double a05=pow(a0,5),a15=pow(a1,5),a25=pow(a2,5);
  if(inte_num==36) //f(x,y,z) = x^5
  {
    basef_surf = x6/6.0;
    return basef_surf;
  }
  if(inte_num==37) //f(x,y,z) = x^4*y
  {
    basef_surf = 0.2*x5*y;
    return basef_surf;
  }
  if(inte_num==38) //f(x,y,z) = x^4*z
  {
    basef_surf = x5*(6.0*a2*y+6.0*a0)+5.0*a1*x6;
    return basef_surf/30.0;
  }
  if(inte_num==39) //f(x,y,z) = x^3*y^2
  {
    basef_surf = 0.25*x4*y2;
    return basef_surf;
  }
  if(inte_num==40) //f(x,y,z) = x^3*yz
  {
    basef_surf = y*(x4*(5.0*a2*y+5.0*a0)+4.0*a1*x5);
    return basef_surf*0.05;
  }
  if(inte_num==41) //f(x,y,z) = x^3*z^2
  {
    basef_surf = x4*(15.0*a22*y2+30.0*a0*a2*y+15.0*a02)+x5*(24.0*a1*a2*y+24.0*a0*a1)+10.0*a12*x6;
    return basef_surf/60.0;
  }
  if(inte_num==42) //f(x,y,z) = x^2*y^3
  {
    basef_surf = x3*y3/3.0;
    return basef_surf;
  }
  if(inte_num==43) //f(x,y,z) = x^2*y^2*z
  {
    basef_surf = y2*(x3*(4.0*a2*y+4.0*a0)+3.0*a1*x4);
    return basef_surf/12.0;
  }
  if(inte_num==44) //f(x,y,z) = x^2*yz^2
  {
    basef_surf = y*(x3*(10.0*a22*y2+20.0*a0*a2*y+10.0*a02)+x4*(15.0*a1*a2*y+15.0*a0*a1)+6.0*a12*x5);
    return basef_surf/30.0;
  }
  if(inte_num==45) //f(x,y,z) = x^2*z^3
  {
    basef_surf = x3*(20.0*a23*y3+60.0*a0*a22*y2+60.0*a02*a2*y+20.0*a03)+x4*(45.0*a1*a22*y2+90.0*a0*a1*a2*y+45.0*a02*a1)+
      x5*(36.0*a12*a2*y+36.0*a0*a12)+10.0*a13*x6;
    return basef_surf/60.0;
  }
  if(inte_num==46) //f(x,y,z) = xy^4
  {
    basef_surf = x2*y4*0.5;
    return basef_surf;
  }
  if(inte_num==47) //f(x,y,z) = xy^3*z
  {
    basef_surf = y3*(x2*(3.0*a2*y+3.0*a0)+2.0*a1*x3);
    return basef_surf/6.0;
  }
  if(inte_num==48) //f(x,y,z) = xy^2*z^2
  {
    basef_surf = y2*(x2*(6.0*a22*y2+12.0*a0*a2*y+6.0*a02)+x3*(8.0*a1*a2*y+8.0*a0*a1)+3.0*a12*x4);
    return basef_surf/12.0;
  }
  if(inte_num==49) //f(x,y,z) = xyz^3
  {
    basef_surf = y*(x2*(10.0*a23*y3+30.0*a0*a22*y2+30.0*a02*a2*y+10.0*a03)+x3*(20.0*a1*a22*y2+40.0*a0*a1*a2*y+
      20.0*a02*a1)+x4*(15.0*a12*a2*y+15.0*a0*a12)+4.0*a13*x5);
    return basef_surf*0.05;
  }
  if(inte_num==50) //f(x,y,z) = xz^4
  {
    basef_surf = x2*(15.0*a24*y4+60.0*a0*a23*y3+90.0*a02*a22*y2+60.0*a03*a2*y+15.0*a04)+x3*(40.0*a1*a23*y3+120.0*a0*a1*a22*y2+
      120.0*a02*a1*a2*y+40.0*a03*a1)+x4*(45.0*a12*a22*y2+90.0*a0*a12*a2*y+45.0*a02*a12)+x5*(24.0*a13*a2*y+24.0*a0*a13)+
      5.0*a14*x6;
    return basef_surf/30.0;
  }
  if(inte_num==51) //f(x,y,z) = y^5
  {
    basef_surf = y5*x;
    return basef_surf;
  }
  if(inte_num==52) //f(x,y,z) = y^4*z
  {
    basef_surf = y4*(a2*x*y+0.5*a1*x2+a0*x);
    return basef_surf;
  }
  if(inte_num==53) //f(x,y,z) = y^3*z^2
  {
    basef_surf = y3*(a22*x*y2+2.0*a0*(a2*x*y+0.5*a1*x2)+a1*a2*x2*y+a12*x3/3.0+a02*x);
    return basef_surf;
  }
  if(inte_num==54) //f(x,y,z) = y^2*z^3
  {
    basef_surf = y2*(a23*x*y3+3.0*a0*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+1.5*a1*a22*x2*y2+3.0*a02*(a2*x*y+
      0.5*a1*x2)+a12*a2*x3*y+0.25*a13*x4+a03*x);
    return basef_surf;
  }
  if(inte_num==55) //f(x,y,z) = yz^4
  {
    basef_surf = y*(a24*x*y4+4.0*a0*(a23*x*y3+1.5*a1*a22*x2*y2+a12*a2*x3*y+0.25*a13*x4)+2.0*a1*a23*x2*y3+
      6.0*a02*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+2.0*a12*a22*x3*y2+4.0*a03*(a2*x*y+0.5*a1*x2)+
      a13*a2*x4*y+0.2*a14*x5+a04*x);
    return basef_surf;
  }
  if(inte_num==56) //f(x,y,z) = z^5
  {
    basef_surf = a25*x*y5+5.0*a0*(a24*x*y4+2.0*a1*a23*x2*y3+2.0*a12*a22*x3*y2+a13*a2*x4*y+a14*x5*0.2)+
      2.5*a1*a24*x2*y4+10.0*a02*(a23*x*y3+1.5*a1*a22*x2*y2+a12*a2*x3*y+0.25*a13*x4)+10.0/3.0*a12*a23*x3*y3+
      10.0*a03*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+2.5*a13*a22*x4*y2+5.0*a04*(a2*x*y+0.5*a1*x2)+a14*a2*x5*y+
      a15*x6/6.0+a05*x;
    return basef_surf;
  }

  double x7=pow(x,7),y6=pow(y,6);
  if(inte_num==57) //f(x,y,z) = x^6
  {
    basef_surf = x7/7.0;
    return basef_surf;
  }
  if(inte_num==58) //f(x,y,z) = x^5*y
  {
    basef_surf = x6*y/6.0;
    return basef_surf;
  }
  if(inte_num==59) //f(x,y,z) = x^5*z
  {
    basef_surf = x6*(7.0*a2*y+7.0*a0)+6.0*a1*x7;
    return basef_surf/42.0;
  }
  if(inte_num==60) //f(x,y,z) = x^4*y^2
  {
    basef_surf = 0.2*x5*y2;
    return basef_surf;
  }
  if(inte_num==61) //f(x,y,z) = x^4*yz
  {
    basef_surf = y*(x5*(6.0*a2*y+6.0*a0)+5.0*a1*x6);
    return basef_surf/30.0;
  }
  if(inte_num==62) //f(x,y,z) = x^4*z^2
  {
    basef_surf = x5*(21.0*a22*y2+42.0*a0*a2*y+21.0*a02)+x6*(35.0*a1*a2*y+35.0*a0*a1)+15.0*a12*x7;
    return basef_surf/105.0;
  }
  if(inte_num==63) //f(x,y,z) = x^3*y^3
  {
    basef_surf = 0.25*x4*y3;
    return basef_surf;
  }
  if(inte_num==64) //f(x,y,z) = x^3*y^2*z
  {
    basef_surf = y2*(x4*(5.0*a2*y+5.0*a0)+4.0*a1*x5);
    return basef_surf*0.05;
  }
  if(inte_num==65) //f(x,y,z) = x^3*y*z^2
  {
    basef_surf = y*(x4*(15.0*a22*y2+30.0*a0*a2*y+15.0*a02)+x5*(24.0*a1*a2*y+24.0*a0*a1)+10.0*a12*x6);
    return basef_surf/60.0;
  }
  if(inte_num==66) //f(x,y,z) = x^3*z^3
  {
    basef_surf = x4*(35.0*a23*y3+105.0*a0*a22*y2+105.0*a02*a2*y+35.0*a03)+x5*(84.0*a1*a22*y2+168.0*a0*a1*a2*y+
      84.0*a02*a1)+x6*(70.0*a12*a2*y+70.0*a0*a12)+20.0*a13*x7;
    return basef_surf/140.0;
  }
  if(inte_num==67) //f(x,y,z) = x^2*y^4
  {
    basef_surf = x3*y4/3.0;
    return basef_surf;
  }
  if(inte_num==68) //f(x,y,z) = x^2*y^3*z
  {
    basef_surf = 4.0*a2*x3*y4+(3.0*a1*x4+4.0*a0*x3)*y3;
    return basef_surf/12.0;
  }
  if(inte_num==69) //f(x,y,z) = x^2*y^2*z^2
  {
    basef_surf = y2*(x3*(10.0*a22*y2+20.0*a0*a2*y+10.0*a02)+x4*(15.0*a1*a2*y+15.0*a0*a1)+6.0*a12*x5);
    return basef_surf/30.0;
  }
  if(inte_num==70) //f(x,y,z) = x^2*yz^3
  {
    basef_surf = y*(x3*(20.0*a23*y3+60.0*a0*a22*y2+60.0*a02*a2*y+20.0*a03)+x4*(45.0*a1*a22*y2+90.0*a0*a1*a2*y+
      45.0*a02*a1)+x5*(36.0*a12*a2*y+36.0*a0*a12)+10.0*a13*x6);
    return basef_surf/60.0;
  }
  if(inte_num==71) //f(x,y,z) = x^2*z^4
  {
    basef_surf = x3*(35.0*a24*y4+140.0*a0*a23*y3+210.0*a02*a22*y2+140.0*a03*a2*y+35.0*a04)+x4*(105.0*a1*a23*y3+
      315.0*a0*a1*a22*y2+315.0*a02*a1*a2*y+105.0*a03*a1)+x5*(126.0*a12*a22*y2+252.0*a0*a12*a2*y+126.0*a02*a12)+
      x6*(70.0*a13*a2*y+70.0*a0*a13)+15.0*a14*x7;
    return basef_surf/105.0;
  }
  if(inte_num==72) //f(x,y,z) = xy^5
  {
    basef_surf = 0.5*x2*y5;
    return basef_surf;
  }
  if(inte_num==73) //f(x,y,z) = xy^4*z
  {
    basef_surf = y4*(x2*(3.0*a2*y+3.0*a0)+2.0*a1*x3);
    return basef_surf/6.0;
  }
  if(inte_num==74) //f(x,y,z) = xy^3*z^2
  {
    basef_surf = y3*(x2*(6.0*a22*y2+12.0*a0*a2*y+6.0*a02)+x3*(8.0*a1*a2*y+8.0*a0*a1)+3.0*a12*x4);
    return basef_surf/12.0;
  }
  if(inte_num==75) //f(x,y,z) = xy^2*z^3
  {
    basef_surf = y2*(x2*(10.0*a23*y3+30.0*a0*a22*y2+30.0*a02*a2*y+10.0*a03)+x3*(20.0*a1*a22*y2+
      40.0*a0*a1*a2*y+20.0*a02*a1)+x4*(15.0*a12*a2*y+15.0*a0*a12)+4.0*a13*x5);
    return basef_surf*0.05;
  }
  if(inte_num==76) //f(x,y,z) = xyz^4
  {
    basef_surf = y*(x2*(15.0*a24*y4+60.0*a0*a23*y3+90.0*a02*a22*y2+60.0*a03*a2*y+15.0*a04)+x3*(40.0*a1*a23*y3+
      120.0*a0*a1*a22*y2+120.0*a02*a1*a2*y+40.0*a03*a1)+x4*(45.0*a12*a22*y2+90.0*a0*a12*a2*y+45.0*a02*a12)+
      x5*(24.0*a13*a2*y+24.0*a0*a13)+5.0*a14*x6);
    return basef_surf/30.0;
  }
  if(inte_num==77) //f(x,y,z) = xz^5
  {
    basef_surf = x2*(21.0*a25*y5+105.0*a0*a24*y4+210.0*a02*a23*y3+210.0*a03*a22*y2+105.0*a04*a2*y+21.0*a05)+
      x3*(70.0*a1*a24*y4+280.0*a0*a1*a23*y3+420.0*a02*a1*a22*y2+280.0*a03*a1*a2*y+70.0*a04*a1)+x4*(105.0*a12*a23*y3+
      315.0*a0*a12*a22*y2+315.0*a02*a12*a2*y+105.0*a03*a12)+x5*(84.0*a13*a22*y2+168.0*a0*a13*a2*y+
      84.0*a02*a13)+x6*(35.0*a14*a2*y+35.0*a0*a14)+6.0*a15*x7;
    return basef_surf/42.0;
  }
  if(inte_num==78) //f(x,y,z) = y^6
  {
    basef_surf = y6*x;
    return basef_surf;
  }
  if(inte_num==79) //f(x,y,z) = y^5*z
  {
    basef_surf = y5*(a2*x*y+0.5*a1*x2+a0*x);
    return basef_surf;
  }
  if(inte_num==80) //f(x,y,z) = y^4*z^2
  {
    basef_surf = y4*(a22*x*y2+2.0*a0*(a2*x*y+0.5*a1*x2)+a1*a2*x2*y+a12*x3/3.0+a02*x);
    return basef_surf;
  }
  if(inte_num==81) //f(x,y,z) = y^3*z^3
  {
    basef_surf = y3*(a23*x*y3+3.0*a0*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+1.5*a1*a22*x2*y2+
      3.0*a02*(a2*x*y+0.5*a1*x2)+a12*a2*x3*y+0.25*a13*x4+a03*x);
    return basef_surf;
  }
  if(inte_num==82) //f(x,y,z) = y^2*z^4
  {
    basef_surf = y2*(a24*x*y4+4.0*a0*(a23*x*y3+1.5*a1*a22*x2*y2+a12*a2*x3*y+0.25*a13*x4)+2.0*a1*a23*x2*y3+
      6.0*a02*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+2.0*a12*a22*x3*y2+4.0*a03*(a2*x*y+0.5*a1*x2)+
      a13*a2*x4*y+0.2*a14*x5+a04*x);
    return basef_surf;
  }
  if(inte_num==83) //f(x,y,z) = yz^5
  {
    basef_surf = y*(a25*x*y5+5.0*a0*(a24*x*y4+2.0*a1*a23*x2*y3+2.0*a12*a22*x3*y2+a13*a2*x4*y+0.2*a14*x5)+
      2.5*a1*a24*x2*y4+10.0*a02*(a23*x*y3+1.5*a1*a22*x2*y2+a12*a2*x3*y+0.25*a13*x4)+10.0/3.0*a12*a23*x3*y3+
      10.0*a03*(a22*x*y2+a1*a2*x2*y+a12*x3/3.0)+2.5*a13*a22*x4*y2+5.0*a04*(a2*x*y+0.5*a1*x2)+
      a14*a2*x5*y+a15*x6/6.0+a05*x);
    return basef_surf;
  }
  if(inte_num==84) //f(x,y,z) = z^6
  {
    if(fabs(a1)>0.00000001)
      basef_surf = std::pow((a0+a1*x+a2*y),7)/7.0/a1;
    else
      basef_surf = std::pow((a0+a2*y),6)*x;
    return basef_surf;
  }
  dserror("The base function for boundarycell integration undefined");
  exit(1);
  return 0.0;

}
#endif
