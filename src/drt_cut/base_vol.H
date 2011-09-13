#ifndef BASE_VOL_H
#define BASE_VOL_H

#include<vector>
//returns function which is used to form the matrix of moment fitting equations
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
        if(base_num==26) //f(x,y,z) = x^2*z^2
                return coordi[0]*coordi[0]*coordi[2]*coordi[2];
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
}
#endif

