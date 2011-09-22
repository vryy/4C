#include "facet_integration.H"
#include "line_integration.H"
#include<iostream>
#include<cmath>

//compute the equation of the plane Ax+By+Cz=D
std::vector<double> GEO::CUT::FacetIntegration::equation_plane(const std::vector<Point*>& corners)
{
        std::vector<double> eqn_plane;

        double x1[3],y1[3],z1[3];
        int mm=0;
        for(std::vector<Point*>::const_iterator k=corners.begin();k!=corners.end();k++)
        {
            const Point* po = *k;
            const double * coords = po->X();
//          std::cout<<coords[0]<<"\t"<<coords[1]<<"\t"<<coords[2]<<std::endl;
            x1[mm] = coords[0];
            y1[mm] = coords[1];
            z1[mm] = coords[2];
            mm++;
           //3 points are sufficient to find the equation of plane
            if(mm==3)
                    break;
//
        }
        double a=0.0,b=0.0,c=0.0,d=0.0;
        a = y1[0]*(z1[1]-z1[2])+y1[1]*(z1[2]-z1[0])+y1[2]*(z1[0]-z1[1]);
        b = z1[0]*(x1[1]-x1[2])+z1[1]*(x1[2]-x1[0])+z1[2]*(x1[0]-x1[1]);
        c = x1[0]*(y1[1]-y1[2])+x1[1]*(y1[2]-y1[0])+x1[2]*(y1[0]-y1[1]);
        d = x1[0]*(y1[1]*z1[2]-y1[2]*z1[1])+x1[1]*(y1[2]*z1[0]-y1[0]*z1[2])+x1[2]*(y1[0]*z1[1]-y1[1]*z1[0]);

        eqn_plane.push_back(a);
        eqn_plane.push_back(b);
        eqn_plane.push_back(c);
        eqn_plane.push_back(d);
//      std::cout<<eqn_plane[0]<<"\t"<<eqn_plane[1]<<"\t"<<eqn_plane[2]<<"\t"<<eqn_plane[3]<<std::endl;

        return eqn_plane;
}

//compute only the x-component of unit-normal vector which is used in further computations
//also determine whether the plane is numbered in clockwise or anticlockwise sense when seen away from the face
void GEO::CUT::FacetIntegration::IsClockwise(const std::vector<double> eqn_plane)
{
        clockwise_ = 0;
        Side* parent = face1_->ParentSide();
        std::vector<Point*>  par_pts;
        const std::vector<Node*> &par_nodes = parent->Nodes();
        for(std::vector<Node*>::const_iterator i=par_nodes.begin();i!=par_nodes.end();i++)
        {
                Node *nod = *i;
                Point *pt = nod->point();
                par_pts.push_back(pt);
        }
        std::vector<double> eqn_par = equation_plane(par_pts);

        bool iscut = face1_->OnCutSide();
        if(iscut)
        {
                if(position_==-2)
                {
                        if(eqn_plane[0]*eqn_par[0]<0.0)
                                clockwise_ = 1;
                }
                else
                {
                        if(eqn_plane[0]*eqn_par[0]>0.0)
                                clockwise_ = 1;
                }
        }       
        else
        {
        //        const std::vector<Side*> &ele_sides = elem1_->Sides();
        /*      for(std::vector<Side*>::const_iterator i=ele_sides.begin();i!=ele_sides.end();i++)              
                {
                        Side*sss = *i;
                        const std::vector<Node*> &nnn = sss->Nodes();
                        for(std::vector<Node*>::const_iterator j=nnn.begin();j!=nnn.end();j++)
                        {
                                Node *nn = *j;
                                double*chh;
                                nn->Coordinates( chh );
                        std::cout<<chh[0]<<"\t"<<chh[1]<<"\t"<<chh[2]<<std::endl;       
                        }
                }*/

//should check whether this is sufficient or do we need to find the number of side in the element and
//use their orientation to get clockwise ordering
                if(eqn_plane[0]*eqn_par[0]<0.0 ||eqn_plane[1]*eqn_par[1]<0.0||eqn_plane[2]*eqn_par[2]<0.0)
                        clockwise_ = 1;
        }
//      std::cout<<clockwise_<<std::endl;
}

//computes x=a1+a2y+a3z from the plane equation
//equation of this form is used to replace x in the line integral
std::vector<double> GEO::CUT::FacetIntegration::compute_alpha(std::vector<double> eqn_plane)
{
        std::vector<double> alfa;
        alfa.resize(3);
        double a = eqn_plane[0];
        double b = eqn_plane[1];
        double c = eqn_plane[2];
        double d = eqn_plane[3];

        alfa[0] = d/a;
        alfa[1] = -1.0*b/a;
        alfa[2] = -1.0*c/a;

//      std::cout<<alfa[0]<<"\t"<<alfa[1]<<"\t"<<alfa[2]<<std::endl;
        return alfa;
}

//perform integration over the facet
double GEO::CUT::FacetIntegration::integrate_facet()
{
        const std::vector<Point*> & corners = face1_->CornerPoints();
        eqn_plane_ = equation_plane(corners);
//      std::cout<<eqn_plane[0]<<"\t"<<eqn_plane[1]<<"\t"<<eqn_plane[2]<<"\t"<<eqn_plane[3]<<std::endl;

// the face is in the x-y or in y-z plane which gives zero facet integral
        if(fabs(eqn_plane_[0])<0.0000001)
                return 0.0;
//x=0 plane which again do not contribute to facet integral
        if(fabs(eqn_plane_[1])<0.0000001 && fabs(eqn_plane_[2])<0.0000001 && fabs(eqn_plane_[3])<0.0000001)
                return 0.0;

        IsClockwise(eqn_plane_);
//      std::cout<<normals<<std::endl;

        std::vector<double> alpha;
        alpha = compute_alpha(eqn_plane_);
//      std::cout<<alpha[0]<<"\t"<<alpha[1]<<"\t"<<alpha[2]<<std::endl;

        //integrating over each line of the facet
        double facet_integ = 0.0;
        for(std::vector<Point*>::const_iterator k=corners.begin();k!=corners.end();k++)
        {
                const Point* po1 = *k;
                const Point* po2;
        //for the last line the end point is the first point of the facet
                if(k!=(corners.end()-1))
                        po2 = *(k+1);
                else
                        po2 = *(corners.begin());
                const double * coords1 = po1->X();
                const double * coords2 = po2->X();

                std::vector<double> coord1,coord2;
                coord1.resize(2);
                coord2.resize(2);
//The facet is projected over y-z plane and then the integration is performed
//so only y- and z-coordinates are passed to make the lines
//[0]-- indicates y and [1] indicate z because we are now in y-z plane
                coord1[0] = coords1[1];
                coord1[1] = coords1[2];
                coord2[0] = coords2[1];
                coord2[1] = coords2[2];

//              std::cout<<coord1[0]<<"\t"<<coord1[1]<<"\t"<<coord2[0]<<"\t"<<coord2[1]<<std::endl;
//
                LineIntegration line1(coord1,coord2,inte_num_,alpha);
                facet_integ += line1.integrate_line();
        }

//this condition results in negative normal for all the lines in the line integral
        if(clockwise_)
        {
                facet_integ = -1.0*facet_integ;
                for(int i=0;i<4;i++)
                       eqn_plane_[i] = -1.0*eqn_plane_[i];
        }
        return facet_integ;
}
