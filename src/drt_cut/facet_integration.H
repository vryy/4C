#ifndef FACET_INTEGRATION_H
#define FACET_INTEGRATION_H

#include<vector>
#include "cut_element.H"
#include "cut_facet.H"

namespace GEO
{
    namespace CUT
    {
        class Element;
        class Facet;

//While performing the facet integration, the points of the facet should be arranged in anti-clockwise manner when looking the surface away from the volume 
//this ensures outward normal vector
class FacetIntegration
{
    public:
        FacetIntegration(Facet *face1,Element *element1,const GEO::CUT::Point::PointPosition posi)
            :face1_(face1),
             elem1_(element1),
             position_(posi)
        {
        }

        ~FacetIntegration(){};
        void set_integ_number(int inte_num){inte_num_ = inte_num;}
        double integrate_facet();
        std::vector<double> get_equation() {return eqn_plane_;}
//      std::vector<std::pair<int,int> > get_connectivity(){return connecti_;};

    private:
        Facet *face1_;
        Element *elem1_;
        const GEO::CUT::Point::PointPosition position_;
        int inte_num_;
        bool clockwise_;
        std::vector<double> eqn_plane_;
        std::vector<double> equation_plane(const std::vector<std::vector<double> > cornersLocal);
        void IsClockwise(const std::vector<double> eqn_plane);
        std::vector<double> compute_alpha(std::vector<double> eqn_plane);
        
    };
}
}

#endif
