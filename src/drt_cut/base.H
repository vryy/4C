#ifndef BASE_H
#define BASE_H

#include<vector>
#include<cmath>
#include<iostream>

//compute the base functions that are integrated along the boundary lines of the irregular volume 
double base_func_line_int(std::vector<double> pt, int inte_num,std::vector<double> alfa)
{
        double basef_line=0.0;
        if(inte_num==1) //f(x,y,z) = 1
        {
                basef_line = alfa[0]*pt[0]+alfa[1]*pt[0]*pt[0]*0.5+alfa[2]*pt[0]*pt[1]; 
                return basef_line;
        }
        if(inte_num==2) //f(x,y,z) = x
        {
                if(fabs(alfa[1])<0.0000001)
                        basef_line = 0.5*pow((alfa[0]+alfa[2]*pt[1]),2)*pt[0];
                else
                        basef_line = pow((alfa[0]+alfa[1]*pt[0]+alfa[2]*pt[1]),3)/6.0/alfa[1];
                return basef_line;
        }
        if(inte_num==3) //f(x,y,z) = y
        {
                basef_line = (0.5*alfa[0]+alfa[1]*pt[0]/3.0+0.5*alfa[2]*pt[1])*pt[0]*pt[0];
                return basef_line;
        }
        if(inte_num==4) //f(x,y,z) = z
        {
                basef_line = (alfa[0]+0.5*alfa[1]*pt[0]+alfa[2]*pt[1])*pt[0]*pt[1];
                return basef_line;
        }

        if(inte_num==5) //f(x,y,z) = x^2
        {
                if(fabs(alfa[1])<0.0000001)
                        basef_line = pow((alfa[0]+alfa[2]*pt[1]),3)*pt[0]/3.0;
                else
                        basef_line = pow((alfa[0]+alfa[1]*pt[0]+alfa[2]*pt[1]),4)/12.0/alfa[1];
                return basef_line;
        }
        if(inte_num==6) //f(x,y,z) = xy 
        {
                basef_line = 6*pow((alfa[2]*pt[0]*pt[1]),2)+(8*alfa[1]*alfa[2]*pow(pt[0],3)+12*alfa[0]*alfa[2]*pow(pt[0],2))*pt[1]+
                        3*pow(alfa[1],2)*pow(pt[0],4)+8*alfa[0]*alfa[1]*pow(pt[0],3)+6*pow((alfa[0]*pt[0]),2);
                return basef_line/24.0;
        }
        if(inte_num==7) //f(x,y,z) = xz
        {
                basef_line = pt[0]*pow((alfa[2]*pt[1]),2)+2*alfa[0]*(alfa[2]*pt[0]*pt[1]+0.5*alfa[1]*pt[0]*pt[0])+alfa[1]*alfa[2]*pt[0]*pt[0]*pt[1]+
                        pow((alfa[1]*pt[0]),2)*pt[0]/3.0+alfa[0]*alfa[0]*pt[0];
               return basef_line*pt[1]*0.5;     
        }
        if(inte_num==8) //f(x,y,z) = y^2
        {
                basef_line = pow(pt[0],3)*(4*alfa[2]*pt[1]+4*alfa[0])+3*alfa[1]*pow(pt[0],4);
                return basef_line/12.0;
        }
        if(inte_num==9) //f(x,y,z) = yz
        {
                basef_line = (pt[0]*pt[0]*(3*alfa[2]*pt[1]+3*alfa[0])+2*alfa[1]*pow(pt[0],3))*pt[1];
                return basef_line/6.0;
        }
        if(inte_num==10) //f(x,y,z) = z^2
        {
                basef_line = pow(pt[1],2)*(alfa[2]*pt[0]*pt[1]+0.5*alfa[1]*pt[0]*pt[0]+alfa[0]*pt[0]);
                return basef_line;
        }

        if(inte_num==11) //f(x,y,z) = x^3
        {
                basef_line = pt[0]*pow((alfa[2]*pt[1]),4)+4*alfa[0]*(pt[0]*pow((alfa[2]*pt[1]),3)+1.5*alfa[1]*pow((alfa[2]*pt[0]*pt[1]),2)+alfa[1]*alfa[1]*alfa[2]*pow((pt[0]),3)*pt[1]
                                +pow((alfa[1]*pt[0]),3)*0.25*pt[0])+2*alfa[1]*pt[0]*pt[0]*pow((alfa[2]*pt[1]),3)+6*alfa[0]*alfa[0]*(pow((alfa[2]*pt[1]),2)*pt[0]+alfa[1]*alfa[2]*pt[0]*pt[0]*pt[1]+
                                pow((alfa[1]*pt[0]),2)*pt[0]/3.0)+2*pow((alfa[1]*alfa[2]*pt[0]*pt[1]),2)*pt[0]+4*pow((alfa[0]),3)*(alfa[2]*pt[0]*pt[1]+0.5*alfa[1]*pt[0]*pt[0])+pow((alfa[1]),3)*
                                alfa[2]*pow((pt[0]),4)*pt[1]+pow((alfa[1]*pt[0]),4)*pt[0]*0.2+pow(alfa[0],4)*pt[0];
                return basef_line*0.25;
        }
        if(inte_num==12) //f(x,y,z) = x^2 y
        {
                basef_line = pt[0]*pt[0]*(10*pow((alfa[2]*pt[1]),3)+30*alfa[0]*pow((alfa[2]*pt[1]),2)+30*alfa[0]*alfa[0]*alfa[2]*pt[1]+10*pow((alfa[0]),3))+
                        pow((pt[0]),3)*(20*alfa[1]*pow((alfa[2]*pt[1]),2)+40*alfa[0]*alfa[1]*alfa[2]*pt[1]+20*alfa[0]*alfa[0]*alfa[1])+
                        pow((pt[0]),4)*(15*alfa[1]*alfa[1]*alfa[2]*pt[1]+15*alfa[0]*alfa[1]*alfa[1])+4*pow(alfa[1],3)*pow(pt[0],5);
                return basef_line/60.0;
        }
        if(inte_num==13) //f(x,y,z) = x^2 z
        {
                basef_line = pt[1]*(pt[0]*pow((alfa[2]*pt[1]),3)+3*alfa[0]*(pt[0]*pow((alfa[2]*pt[1]),2)+alfa[1]*alfa[2]*pt[0]*pt[0]*pt[1]+pow((alfa[1]*pt[0]),2)*pt[0]/3.0)+
                               1.5*alfa[1]*pow((alfa[2]*pt[0]*pt[1]),2)+3*alfa[0]*alfa[0]*(alfa[2]*pt[0]*pt[1]+0.5*alfa[1]*pt[0]*pt[0])+alfa[2]*pow((alfa[1]*pt[0]),2)*pt[0]*pt[1]+
                               0.25*pow(alfa[1],3)*pow(pt[0],4)+pow(alfa[0],3)*pt[0]);
                return basef_line/3.0;
        }
        if(inte_num==14) //f(x,y,z) = xy^2
        {
                basef_line = pow(pt[0],3)*(10*pow((alfa[2]*pt[1]),2)+20*alfa[0]*alfa[2]*pt[1]+10*alfa[0]*alfa[0])+pow(pt[0],4)*(15*alfa[1]*alfa[2]*pt[1]+15*alfa[0]*alfa[1])+
                        6*alfa[1]*alfa[1]*pow(pt[0],5);
                return basef_line/60.0;
        }
        if(inte_num==15) //f(x,y,z) = xyz
        {
                basef_line = pt[1]*(pt[0]*pt[0]*(6*pow((alfa[2]*pt[1]),2)+12*alfa[0]*alfa[2]*pt[1]+6*alfa[0]*alfa[0])+pow(pt[0],3)*8*(alfa[1]*alfa[2]*pt[1]+alfa[0]*alfa[1])+
                                3*alfa[1]*alfa[1]*pow(pt[0],4));
                return basef_line/24.0;
        }
        if(inte_num==16) //f(x,y,z) = xz^2
        {
                basef_line = 3*alfa[2]*alfa[2]*pt[0]*pow(pt[1],4)+pow(pt[1],3)*(3*alfa[1]*alfa[2]*pt[0]*pt[0]+6*alfa[0]*alfa[2]*pt[0])+pt[1]*pt[1]*(alfa[1]*alfa[1]*pow(pt[0],3)+
                                3*alfa[0]*alfa[1]*pt[0]*pt[0]+3*alfa[0]*alfa[0]*pt[0]);
                return basef_line/6.0;
        }
        if(inte_num==17) //f(x,y,z) = y^3
        {
                basef_line = pow(pt[0],4)*5*(alfa[2]*pt[1]+alfa[0])+4*alfa[1]*pow(pt[0],5);
                return basef_line*0.05;
        }
        if(inte_num==18) //f(x,y,z) = y^2 z
        {
                basef_line = 4*alfa[2]*pow((pt[0]*pt[1]),2)*pt[0]+pt[1]*(3*alfa[1]*pow(pt[0],4)+4*alfa[0]*pow(pt[0],3));
                return basef_line/12.0;
        }
        if(inte_num==19) //f(x,y,z) = yz^2
        {
                basef_line = pt[1]*pt[1]*(pt[0]*pt[0]*(3*alfa[2]*pt[1]+3*alfa[0])+2*alfa[1]*pow(pt[0],3));
                return basef_line/6.0;
        }
        if(inte_num==20) //f(x,y,z) = z^3
        {
                basef_line = pow(pt[1],3)*(alfa[2]*pt[0]*pt[1]+0.5*alfa[1]*pt[0]*pt[0]+alfa[0]*pt[0]);
                return basef_line;
        }

        double a02=pow(alfa[0],2),a12=pow(alfa[1],2),a22=pow(alfa[2],2),y2=pow(pt[0],2),z2=pow(pt[1],2);
        double a03=pow(alfa[0],3),a13=pow(alfa[1],3),a23=pow(alfa[2],3),y3=pow(pt[0],3),z3=pow(pt[1],3);
        double a04=pow(alfa[0],4),a14=pow(alfa[1],4),a24=pow(alfa[2],4),y4=pow(pt[0],4),z4=pow(pt[1],4);
        double a05=pow(alfa[0],5),a15=pow(alfa[1],5),a25=pow(alfa[2],5),y5=pow(pt[0],5),z5=pow(pt[1],5);
        double y6=pow(pt[0],6);
        if(inte_num==21) //f(x,y,z) = x^4
        {
                basef_line = 6*a25*pt[0]*z5+(15*alfa[1]*y2+30*alfa[0]*pt[0])*a24*z4+(20*a12*y3+60*alfa[0]*alfa[1]*y2+60*a02*pt[0])*a23*z3+
                        (15*a13*y4+60*alfa[0]*a12*y3+90*a02*alfa[1]*y2+60*a03*pt[0])*a22*z2+(6*a14*y5+30*alfa[0]*a13*y4+60*a02*a12*y3+
                        60*a03*alfa[1]*y2+30*a04*pt[0])*alfa[2]*pt[1]+a15*y6+6*alfa[0]*a14*y5+15*a02*a13*y4+20*a03*a12*y3+15*a04*alfa[1]*y2+6*a05*pt[0];
                return basef_line/30.0;
        }
        if(inte_num==22) //f(x,y,z) = x^3*y
        {
                basef_line = y2*(15*a24*z4+60*alfa[0]*a23*z3+90*a02*a22*z2+60*a03*alfa[2]*pt[1]+15*a04)+y3*(40*alfa[1]*a23*z3+120*alfa[0]*alfa[1]*a22*z2
                        +120*a02*alfa[1]*alfa[2]*pt[1]+40*a03*alfa[1])+y4*(45*a12*a22*z2+90*alfa[0]*a12*alfa[2]*pt[1]+45*a02*a12)+y5*(24*a13*alfa[2]*pt[1]
                        +24*alfa[0]*a13)+5*a14*y6;
                return basef_line/120.0;
        }
        if(inte_num==23) //f(x,y,z) = x^3*z
        {
                basef_line = pt[1]*(a24*pt[0]*z4+4*alfa[0]*(a23*pt[0]*z3+1.5*alfa[1]*a22*y2*z2+a12*alfa[2]*y3*pt[1]+0.25*a13*y4)+2*alfa[1]*a23*y2*z3
                                +6*a02*(a22*pt[0]*z2+alfa[1]*alfa[2]*y2*pt[1]+a12*y3/3.0)+2*a12*a22*y3*z2+4*a03*(alfa[2]*pt[0]*pt[1]+alfa[1]*y2*0.5)
                                +a13*alfa[2]*y4*pt[1]+0.2*a14*y5+a04*pt[0]);
                return 0.25*basef_line;
        }
        if(inte_num==24) //f(x,y,z) = x^2*y^2
        {
                basef_line = y3*(20*a23*z3+60*alfa[0]*a22*z2+60*a02*alfa[2]*pt[1]+20*a03)+y4*(45*alfa[1]*a22*z2+90*alfa[0]*alfa[1]*alfa[2]*pt[1]+45*a02*alfa[1])
                                +y5*(36*a12*alfa[2]*pt[1]+36*alfa[0]*a12)+10*a13*y6;
                return basef_line/180.0;
        }
        if(inte_num==25) //f(x,y,z) = x^2*yz
        {
                basef_line = pt[1]*(y2*(10*a23*z3+30*alfa[0]*a22*z2+30*a02*alfa[2]*pt[1]+10*a03)+y3*(20*alfa[1]*a22*z2+40*alfa[0]*alfa[1]*alfa[2]*pt[1]
                                +20*a02*alfa[1])+y4*(15*a12*alfa[2]*pt[1]+15*alfa[0]*a12)+4*a13*y5);
                return basef_line/60.0;
        }
        if(inte_num==26) //f(x,y,z) = x^2*z^2
        {
                basef_line = z2*(a23*pt[0]*z3+3*alfa[0]*(a22*pt[0]*z2+alfa[1]*alfa[2]*y2*pt[1]+a12*y3/3.0)+1.5*alfa[1]*a22*y2*z2+3*a02*(alfa[2]*pt[0]*pt[1]+
                                0.5*alfa[1]*y2)+a12*alfa[2]*y3*pt[1]+0.25*a13*y4+a03*pt[0]);
                return basef_line/3.0;
        }
        if(inte_num==27) //f(x,y,z) = x*y^3
        {
                basef_line = y4*(15*a22*z2+30*alfa[0]*alfa[2]*pt[1]+15*a02)+y5*(24*alfa[1]*alfa[2]*pt[1]+24*alfa[0]*alfa[1])+10*a12*y6;
                return basef_line/120.0;
        }
        if(inte_num==28) //f(x,y,z) = x*y^2*z
        {
                basef_line = pt[1]*(y3*(10*a22*z2+20*alfa[0]*alfa[2]*pt[1]+10*a02)+y4*(15*alfa[1]*alfa[2]*pt[1]+15*alfa[0]*alfa[1])+6*a12*y5);
                return basef_line/60.0;
        }
        if(inte_num==29) //f(x,y,z) = x*y*z^2
        {
                basef_line = z2*(y2*(6*a22*z2+12*alfa[0]*alfa[2]*pt[1]+6*a02)+y3*(8*alfa[1]*alfa[2]*pt[1]+8*alfa[0]*alfa[1])+3*a12*y4);
                return basef_line/24.0;
        }
        if(inte_num==30) //f(x,y,z) = x*z^3
        {
                basef_line = z3*(a22*pt[0]*z2+2*alfa[0]*(alfa[2]*pt[0]*pt[1]+0.5*alfa[1]*y2)+alfa[1]*alfa[2]*y2*pt[1]+a12*y3/3.0+a02*pt[0]);
                return basef_line*0.5;
        }
        if(inte_num==31) //f(x,y,z) = y^4
        {
                basef_line = y5*(6*alfa[2]*pt[1]+6*alfa[0])+5*alfa[1]*y6;
                return basef_line/30.0;
        }
        if(inte_num==32) //f(x,y,z) = y^3*z
        {
                basef_line = y4*(5*alfa[2]*pt[1]+5*alfa[0])+4*alfa[1]*y5;
                return basef_line*0.05*pt[1];
        }
        if(inte_num==33) //f(x,y,z) = y^2*z^2
        {
                basef_line = 4*y3*(alfa[2]*pt[1]+alfa[0])+3*alfa[1]*y4;
                return basef_line/12.0*z2;
        }
        if(inte_num==34) //f(x,y,z) = y*z^3
        {
                basef_line = y2*(3*alfa[2]*pt[1]+3*alfa[0])+2*alfa[1]*y3;
                return basef_line*z3/6.0;
        }
        if(inte_num==35) //f(x,y,z) = y*z^3
        {
                basef_line = alfa[2]*pt[0]*pt[1]+0.5*alfa[1]*y2+alfa[0]*pt[0];
                return basef_line*z4;
        }
}

/*//returns function which is used to form the matrix of moment fitting equations
double base_function(std::vector<double> coordi,int base_num)
{
        if(base_num==1) //f(x,y,z) = 1.0
                return 1.0;
        if(base_num==2) //f(x,y,z) = x
                return coordi[0];
        if(base_num==3) //f(x,y,z) = y
                return coordi[1];
        if(base_num==4) //f(x,y,z) = z
                return coordi[2];

        if(base_num==5) //f(x,y,z) = x^2
                return coordi[0]*coordi[0];
        if(base_num==6) //f(x,y,z) = xy
                return coordi[0]*coordi[1];
        if(base_num==7) //f(x,y,z) = xz
                return coordi[0]*coordi[2];
        if(base_num==8) //f(x,y,z) = y^2
                return coordi[1]*coordi[1];
        if(base_num==9) //f(x,y,z) = yz
                return coordi[1]*coordi[2];
        if(base_num==10) //f(x,y,z) = z^2
                return coordi[2]*coordi[2];

        if(base_num==11) //f(x,y,x) = x^3
                return coordi[0]*coordi[0]*coordi[0];
        if(base_num==12) //f(x,y,x) = x^2*y
                return coordi[0]*coordi[0]*coordi[1];
        if(base_num==13) //f(x,y,x) = x^2*z
                return coordi[0]*coordi[0]*coordi[2];
        if(base_num==14) //f(x,y,x) = xy^2
                return coordi[0]*coordi[1]*coordi[1];
        if(base_num==15) //f(x,y,x) = xyz
                return coordi[0]*coordi[1]*coordi[2];
        if(base_num==16) //f(x,y,x) = xz^2
                return coordi[0]*coordi[2]*coordi[2];
        if(base_num==17) //f(x,y,x) = y^3
                return coordi[1]*coordi[1]*coordi[1];
        if(base_num==18) //f(x,y,x) = y^2*z
                return coordi[1]*coordi[1]*coordi[2];
        if(base_num==19) //f(x,y,x) = yz^2
                return coordi[1]*coordi[2]*coordi[2];
        if(base_num==20) //f(x,y,x) = z^3
                return coordi[2]*coordi[2]*coordi[2];

        if(base_num==21) //f(x,y,z) = x^4
                return coordi[0]*coordi[0]*coordi[0]*coordi[0];
        if(base_num==22) //f(x,y,z) = x^3*y
                return coordi[0]*coordi[0]*coordi[0]*coordi[1];
        if(base_num==23) //f(x,y,z) = x^3*z
                return coordi[0]*coordi[0]*coordi[0]*coordi[2];
        if(base_num==24) //f(x,y,z) = x^2*y^2
                return coordi[0]*coordi[0]*coordi[1]*coordi[1];
        if(base_num==25) //f(x,y,z) = x^2*yz
                return coordi[0]*coordi[0]*coordi[1]*coordi[2];
        if(base_num==26) //f(x,y,z) = x^2*y^2
                return coordi[0]*coordi[0]*coordi[1]*coordi[1];
        if(base_num==27) //f(x,y,z) = xy^3
                return coordi[0]*coordi[1]*coordi[1]*coordi[1];
        if(base_num==28) //f(x,y,z) = xy^2*z
                return coordi[0]*coordi[1]*coordi[1]*coordi[2];
        if(base_num==29) //f(x,y,z) = xyz^2
                return coordi[0]*coordi[1]*coordi[2]*coordi[2];
        if(base_num==30) //f(x,y,z) = xz^3
                return coordi[0]*coordi[2]*coordi[2]*coordi[2];
        if(base_num==31) //f(x,y,z) = y^4
                return coordi[1]*coordi[1]*coordi[1]*coordi[1];
        if(base_num==32) //f(x,y,z) = y^3*z
                return coordi[1]*coordi[1]*coordi[1]*coordi[2];
        if(base_num==33) //f(x,y,z) = y^2*z^2
                return coordi[1]*coordi[1]*coordi[2]*coordi[2];
        if(base_num==34) //f(x,y,z) = yz^3
                return coordi[1]*coordi[2]*coordi[2]*coordi[2];
        if(base_num==35) //f(x,y,z) = z^4
                return coordi[2]*coordi[2]*coordi[2]*coordi[2];
}*/
#endif

