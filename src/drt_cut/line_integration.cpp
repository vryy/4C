#include "line_integration.H"
#include<cmath>
#include<iostream>
#include "base.H"

std::vector<double> LineIntegration::compute_normal()
{
        std::vector<double> normal;
        normal.resize(2);
        double dy = point_end_[1]-point_begin_[1];
        double dx = -point_end_[0]+point_begin_[0];
        double modd = sqrt(dx*dx+dy*dy);

        normal[0] = dy/modd;
        normal[1] = dx/modd;
//      std::cout<<normal[0]<<"\t"<<normal[1]<<std::endl;
        return normal;
}

std::vector<double> LineIntegration::get_Gauss_weights()
{
        std::vector<double> line_wei;
        line_wei.push_back(0.1012285362903762591525314);
        line_wei.push_back(0.2223810344533744705443560);
        line_wei.push_back(0.3137066458778872873379622);
        line_wei.push_back(0.3626837833783619829651504);
        line_wei.push_back(0.3626837833783619829651504);
        line_wei.push_back(0.3137066458778872873379622);
        line_wei.push_back(0.2223810344533744705443560);
        line_wei.push_back(0.1012285362903762591525314);

        return line_wei;
}

std::vector<double> LineIntegration::get_Gauss_line_pts()
{
        std::vector<double> line_tau;
        line_tau.push_back(-0.9602898564975362316835609);
        line_tau.push_back(-0.7966664774136267395915539);
        line_tau.push_back(-0.5255324099163289858177390);
        line_tau.push_back(-0.1834346424956498049394761);
        line_tau.push_back(0.1834346424956498049394761);
        line_tau.push_back(0.5255324099163289858177390);
        line_tau.push_back(0.7966664774136267395915539);
        line_tau.push_back(0.9602898564975362316835609);

        return line_tau;
}

//Obtain the actual integration points from the points available in (-1,1) interval
std::vector<std::vector<double> > LineIntegration::find_line_integration_pts()
{
        std::vector<double> line_tau = get_Gauss_line_pts();

/*      for(int i=0;i<3;i++)
                std::cout<<line_tau[i]<<std::endl;*/

        std::vector<std::vector<double> > line_int_pts;
        line_int_pts.resize(2);
        for(int i=0;i<2;i++)
                line_int_pts[i].resize(8);
        double xmid[2];
        //middle point in all 2 coordinates
        for(int i=0;i<2;i++)
        {
                xmid[i] = 0.5*(point_begin_[i]+point_end_[i]);
        }
        //Finding the 3 integration points used in Gaussian quadrature
        for(int i=0;i<8;i++)
        {
                double tau_fac = line_tau[i];
                for(int j=0;j<2;j++)
                {
                        line_int_pts[j][i] = (xmid[j]-point_begin_[j])*tau_fac+xmid[j];
                }
        }
        
/*      for(int i=0;i<3;i++)
        {
                for(int j=0;j<2;j++)
                        std::cout<<line_int_pts[j][i]<<std::endl;
        }*/
        return line_int_pts;
}

//The half length of the line appears in the Gaussian quadrature when the integration
//is scaled from (a,b) to (-1,1) for which the weights are available 
double LineIntegration::half_length()
{
        double x1 = point_begin_[0];
        double y1 = point_begin_[1];
        
        double x2 = point_end_[0];
        double y2 = point_end_[1];

        double hal_len = 0.5*sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        return hal_len;
}

double LineIntegration::integrate_line()
{
        std::vector<double> normal;
        normal = compute_normal();
//      std::cout<<normal[0]<<"\t"<<normal[1]<<std::endl;

        if (fabs(normal[0])<0.000000001)
                return 0.0;

        std::vector<std::vector<double> > line_int_pts;
        line_int_pts = find_line_integration_pts();

        std::vector<double> line_wei = get_Gauss_weights();

        double half_len = half_length();

//      std::cout<<half_len<<std::endl;
//
        double inte = 0.0;
        for (int i=0;i<8;i++)
        {
                std::vector<double> pt;
                pt.push_back(line_int_pts[0][i]);
                pt.push_back(line_int_pts[1][i]);

                double linein = base_func_line_int(pt, inte_num_,alpha_);
        //      std::cout<<"base = "<<linein<<std::endl;
                inte = inte+line_wei[i]*linein;
        }
        inte = inte*normal[0]*half_len;
//      std::cout<<"line_inte = "<<inte<<std::endl;

        return inte;
}
