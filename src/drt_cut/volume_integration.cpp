#include<iostream>
#include<fstream>
#include <cmath>
#include<algorithm>
#include "volume_integration.H"
#include "base_vol.H"
#include "least_squares.H"

//compute the rhs of the moment fitting equations
//Integration of base functions take place inside this
std::vector<double> GEO::CUT::VolumeIntegration::compute_rhs_moment()
{
    std::vector<double> rhs_mom;

    const plain_facet_set & facete = volcell_->Facets();
    for(int fnc=1;fnc<=num_func_;fnc++)
    {
        double mome = 0.0;

        for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
        {
            Facet *fe = *i;
            FacetIntegration faee1(fe,elem1_,position_);
            faee1.set_integ_number(fnc);
            mome += faee1.integrate_facet();
        }
        rhs_mom.push_back(mome);
//      std::cout<<mome<<std::endl;
    }
    return rhs_mom;
}

//compute the gaussian points of the volumecell
void GEO::CUT::VolumeIntegration::compute_Gaussian_points()
{
    
/*  for(std::vector<FacetIntegration>::iterator i=facet_list_.begin();
        i!=facet_list_.end();i++)
    {
        FacetIntegration facee = *i;
    }*/


/*************** probably need to remove this*****************************************/
//create a bounding box over the volume and generate gaussian points over that volume   
    const plain_facet_set & facete = volcell_->Facets();
    std::vector<double> x1,y1,z1;
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
        Facet *fac = *i;
        const std::vector<Point*> & corners = fac->CornerPoints();
        for(std::vector<Point*>::const_iterator k=corners.begin();k!=corners.end();k++)
        {
            const Point* po = *k;
            const double * coords = po->X();
            x1.push_back(coords[0]);
            y1.push_back(coords[1]);
            z1.push_back(coords[2]);
        }

    }
    std::sort(x1.begin(),x1.end());
    std::sort(y1.begin(),y1.end());
    std::sort(z1.begin(),z1.end());

    unsigned numb = x1.size();
    double minnx=x1[0],maxxx=x1[numb-1],minny=y1[0],maxxy=y1[numb-1],minnz=z1[0],maxxz=z1[numb-1];  
//    double tauu[] = {0.1,0.2,0.5,0.6};
//    double tauuz[] = {0.1,0.2,0.5,0.6};
    double tauux[] = {0.2,0.4, 0.55,0.7,0.9};
    double tauuy[] = {0.2,0.4, 0.55,0.7,0.9};
    double tauuz[] = {0.2,0.4, 0.55,0.7,0.9};
    for(int i=0;i<5;i++)
    {
        
        std::vector<double> ptt;
        ptt.resize(3);
        ptt[2] = minnz+(maxxz-minnz)*tauuz[i];
        for(int j=0;j<5;j++)
        {

            ptt[1] = minny+(maxxy-minny)*tauuy[j];
            for(int k=0;k<5;k++)
            {
                ptt[0] = minnx+(maxxx-minnx)*tauux[k];
                gaus_pts_.push_back(ptt);
            }
        }
    }
//    std::cout<<"size"<<gaus_pts_.size()<<std::endl;
/*  for(int i=0;i<27;i++)
    {
        std::cout<<gau_pts[i][0]<<"\t"<<gau_pts[i][1]<<"\t"<<gau_pts[i][2]<<std::endl;
    }*/
/*******************until this*******************************************/
}

//form the moment fitting matrix
void GEO::CUT::VolumeIntegration::moment_fitting_matrix(std::vector<std::vector<double> >&mom)
{
    for(int i=0;i<num_func_;i++)
    {
        int k=0;
        for(std::vector<std::vector<double> >::iterator j=gaus_pts_.begin();j!=gaus_pts_.end();j++)
        {
            std::vector<double> cordi = *j;
            mom[i][k] = base_function(cordi,i+1);
            k++;
        }
    }

/*  std::ofstream myfile;
      myfile.open ("example.txt");
      for(int i=0;i<mom.size();i++)
      {
          for(int j=0;j<mom[0].size();j++)
          {
          myfile << mom[i][j]<<"\t";
          }
          myfile<<"\n";
      }
      myfile.close();*/

}

//compute the weights and returns the coordinates of weights and the weighing points
std::vector<double> GEO::CUT::VolumeIntegration::compute_weights()
{
    std::vector<double> rhs_moment;
    rhs_moment = compute_rhs_moment();
    
/*  for(std::vector<double>::iterator i=rhs_moment.begin();i!=rhs_moment.end();i++)
        std::cout<<*i<<std::endl;*/

    compute_Gaussian_points();
/*  for(int i=0;i<27;i++)
    {
        std::cout<<gaus_pts_[i][0]<<"\t"<<gaus_pts_[i][1]<<"\t"<<gaus_pts_[i][2]<<std::endl;
    }*/

    std::vector<std::vector<double> > moment_matrix(num_func_,std::vector<double>(gaus_pts_.size()));
    moment_fitting_matrix(moment_matrix);

    std::vector<double> weights;
    LeastSquares least(moment_matrix,rhs_moment);
    weights = least.linear_least_square();

/*   for(int j=0;j<num_func_;j++)
    {
    double chek = 0.0;
    for(int i=0;i<weights.size();i++)
    {
        std::vector<double> coorrr;
        coorrr.push_back(gaus_pts_[i][0]);
        coorrr.push_back(gaus_pts_[i][1]);
        coorrr.push_back(gaus_pts_[i][2]);
        chek += weights[i]*base_function(coorrr,j+1);
    }
    std::cout<<"check"<<rhs_moment[j]-chek<<std::endl;
    }*/

/*  double chek=0.0;
    for(int i=0;i<weights.size();i++)
    {
    //  chek += weights[i]*(pow(gaus_pts_[i][0],2)+pow(gaus_pts_[i][1],2)+5.0);     
    //  chek += weights[i]*(gaus_pts_[i][0]*gaus_pts_[i][0]*gaus_pts_[i][0]+gaus_pts_[i][1]*gaus_pts_[i][1]*gaus_pts_[i][1]+5.0);       
        chek += pow(gaus_pts_[i][0],4)*weights[i];
    }
    std::cout<<"check"<<chek<<std::endl;
    std::cout<<"check"<<chek-(rhs_moment[10]+rhs_moment[16]+5*rhs_moment[0])<<std::endl;*/

/*  for(int i=0;i<weights.size();i++)
        std::cout<<"wei"<<weights[i]<<std::endl;*/

    return weights;
}
