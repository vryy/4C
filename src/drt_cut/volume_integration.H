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
        std::vector<std::vector<double> > eqn_facets_;
        bool compute_Gaussian_points(int numeach);
        void moment_fitting_matrix(std::vector<std::vector<double> >&mom);
        bool IsIntersect(double *pt, double *mini, double *maxi,vector<vector<double> >& linePts,
                        vector<vector<double> >zcoord,vector<vector<double> >ycoord,double toler,int numeach);
        bool IsContainArea(double minn[3],double maxx[3], double&zmin,vector<vector<double> > &pts,
                vector<vector<double> >zcoord,vector<vector<double> >ycoord,double toler,int numeach);
        void GaussPointGmsh();
        void OnLine(vector<double>inter1,vector<double>inter2,vector<vector<double> >&linePts,int num);
        void get_zcoordinates(vector<vector<double> >& zcoord,vector<vector<double> >& ycoord);
        int pnpoly(int npol, vector<double>xp, vector<double>yp, double x, double y);
	bool traceVolume(vector<vector<double> >& zcoord,vector<vector<double> >& ycoord,double*minn,double*maxx,double toler);
};
}
}

#endif
