#ifndef VOLUME_INTEGRATION_H
#define VOLUME_INTEGRATION_H

#include<vector>
#include "facet_integration.H"
#include "cut_volumecell.H"
#include "cut_element.H"
#include "cut_facet.H"
namespace GEO
{
    namespace CUT
    {
        class Element;
        class Facet;
        class VolumeCell;
//While performing the volume integration, the points of the facet should be arranged in anti-clockwise manner when looking the surface away from the body
//this ensures outward normal vector
class VolumeIntegration
{
    public:
        VolumeIntegration(VolumeCell* volcell, Element* elem, const GEO::CUT::Point::PointPosition posi, int num_func)
            :volcell_(volcell),
            elem1_(elem),
            position_(posi),
            num_func_(num_func)
        {
        }

        ~VolumeIntegration(){};
        std::vector<double> compute_weights();

    private:
        VolumeCell* volcell_;
        Element* elem1_;
        const GEO::CUT::Point::PointPosition position_;
        int num_func_;
        std::vector<double> compute_rhs_moment();
        std::vector<std::vector<double> > gaus_pts_;
        void compute_Gaussian_points();
        void moment_fitting_matrix(std::vector<std::vector<double> >&mom);
};
}
}

#endif
