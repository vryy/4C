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

