#include "facet_integration.H"
#include "line_integration.H"
#include<iostream>
#include<cmath>

//compute the equation of the plane Ax+By+Cz=D with the local coordinates of corner points
std::vector<double> GEO::CUT::FacetIntegration::equation_plane(const std::vector<std::vector<double> > cornersLocal)
{
  std::vector<double> eqn_plane(4);

  double x1[3]={0.0,0.0,0.0},y1[3]={0.0,0.0,0.0},z1[3]={0.0,0.0,0.0};
  unsigned mm=0;
  for(std::vector<std::vector<double> >::const_iterator k=cornersLocal.begin();k!=cornersLocal.end();k++)
  {
    const std::vector<double> coords = *k;
//          std::cout<<coords[0]<<"\t"<<coords[1]<<"\t"<<coords[2]<<std::endl;
    x1[mm] = coords[0];
    y1[mm] = coords[1];
    z1[mm] = coords[2];

//make sure all three points do not lie on a line
//construct a vector pt2-pt1 and pt3-pt1. If the cross product is zero then all are in same line
//go to next point in such cases
    if(mm==2)
    {
      double pt1pt2[3],pt1pt3[3],cross[3];
      pt1pt2[0] = x1[1]-x1[0];pt1pt3[0]=x1[2]-x1[0];
      pt1pt2[1] = y1[1]-y1[0];pt1pt3[1]=y1[2]-y1[0];
      pt1pt2[2] = z1[1]-z1[0];pt1pt3[2]=z1[2]-z1[0];
      cross[0] = fabs(pt1pt2[1]*pt1pt3[2]-pt1pt2[2]*pt1pt3[1]);
      cross[1] = fabs(pt1pt2[0]*pt1pt3[2]-pt1pt2[2]*pt1pt3[0]);
      cross[2] = fabs(pt1pt2[1]*pt1pt3[0]-pt1pt2[0]*pt1pt3[1]);
        if(cross[0]<0.000001 && cross[1]<0.000001 && cross[2]<0.000001)
          continue;
        else
          mm++;
      }
      else
        mm++;
      //3 points are sufficient to find the equation of plane
      if(mm==3)
        break;
      if(mm==cornersLocal.size())
        dserror("All points are on a line");
  }
  eqn_plane[0] = y1[0]*(z1[1]-z1[2])+y1[1]*(z1[2]-z1[0])+y1[2]*(z1[0]-z1[1]);
  eqn_plane[1] = z1[0]*(x1[1]-x1[2])+z1[1]*(x1[2]-x1[0])+z1[2]*(x1[0]-x1[1]);
  eqn_plane[2] = x1[0]*(y1[1]-y1[2])+x1[1]*(y1[2]-y1[0])+x1[2]*(y1[0]-y1[1]);
  eqn_plane[3] = x1[0]*(y1[1]*z1[2]-y1[2]*z1[1])+x1[1]*(y1[2]*z1[0]-y1[0]*z1[2])+x1[2]*(y1[0]*z1[1]-y1[1]*z1[0]);
//      std::cout<<eqn_plane[0]<<"\t"<<eqn_plane[1]<<"\t"<<eqn_plane[2]<<"\t"<<eqn_plane[3]<<std::endl;

  return eqn_plane;
}

//compute only the x-component of unit-normal vector which is used in further computations
//also determine whether the plane is numbered in clockwise or anticlockwise sense when seen away from the face
void GEO::CUT::FacetIntegration::IsClockwise(const std::vector<double> eqn_plane,const std::vector<std::vector<double> > cornersLocal)
{
  std::string ordering;

  if(cornersLocal.size()==3 || cornersLocal.size()==4)
  {
    if(eqn_plane[0]>0.0)
      ordering = "acw";
    else
      ordering = "cw";
  }
#if 1
  else
  {
    double crossProd = 0.0;
    for(unsigned i=0;i<cornersLocal.size();i++)
    {
      unsigned j = (i+1)%cornersLocal.size();
      crossProd += (cornersLocal[j][1]-cornersLocal[i][1])*(cornersLocal[j][2]+cornersLocal[i][2]);
    }
    if(crossProd>0)
      ordering = "cw";
    else if(crossProd<0)
      ordering = "acw";
    else
      dserror("the points in the facet are neither ordered in clockwise not anti-clockwise or all collinear"
          " point in this facet");
  }
#else
  else
  {
    double crossProd = 0.0;
    int count=0;
    for(unsigned i=0;i<cornersLocal.size();i++)
    {
      unsigned j = (i+1)%cornersLocal.size();
      unsigned k = (i+2)%cornersLocal.size();
      crossProd = (cornersLocal[j][1]-cornersLocal[i][1])*(cornersLocal[k][2]-cornersLocal[j][2]);
      crossProd -= (cornersLocal[j][2]-cornersLocal[i][2])*(cornersLocal[k][1]-cornersLocal[j][1]);
      std::cout<<"cross = "<<crossProd<<"\n";
      if(crossProd<0)
        count--;
      else if(crossProd>0)
        count++;

    }
    if(count>0)
      ordering = "acw";
    else if(count<0)
      ordering = "cw";
    else
      dserror("the points in the facet are neither ordered in clockwise not anti-clockwise or all collinear"
          " point in this facet");
  }
#endif

  clockwise_ = 0;
  Side* parent = face1_->ParentSide();
  const std::vector<Node*> &par_nodes = parent->Nodes();
  std::vector<vector<double> > corners(par_nodes.size());
  int mm=0;
  for(std::vector<Node*>::const_iterator i=par_nodes.begin();i!=par_nodes.end();i++)
  {
    Node *nod = *i;
    double x1[3];
    nod->Coordinates(x1);
    LINALG::Matrix<3,1> glo,loc;
    std::vector<double> pt_local(3);
    for(int nodno=0;nodno<3;nodno++)
      glo(nodno,0) = x1[nodno];

     elem1_->LocalCoordinates(glo,loc);

     pt_local[0] = loc(0,0);
     pt_local[1] = loc(1,0);
     pt_local[2] = loc(2,0);

  //			 std::cout<<loc(0,0)<<"\t"<<loc(1,0)<<"\t"<<loc(2,0)<<"\t";//blockkk

     corners[mm] = pt_local;
     mm++;
    }
//        std::cout<<"\n";//blockkk
    std::vector<double> eqn_par = equation_plane(corners);

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
			const std::vector<Side*> &ele_sides = elem1_->Sides();
//				std::cout<<"element nodes"<<"\n";//blockkk or remove until the next comment
			int parentSideno = 0;
			for(std::vector<Side*>::const_iterator i=ele_sides.begin();i!=ele_sides.end();i++)
			{
					Side*sss = *i;
					if(sss==parent)
					{
						break;	//unblockkk
					}
					parentSideno++;
/*						std::cout<<"side nodes"<<"\n";//blockkk
					const std::vector<Node*> &nnn = sss->Nodes();
					for(std::vector<Node*>::const_iterator j=nnn.begin();j!=nnn.end();j++)
					{
						Node *nn = *j;
						double*chh;
						nn->Coordinates( chh );
						LINALG::Matrix<3,1> glo,loc;
						glo(0,0) = chh[0];
						glo(1,0) = chh[1];
						glo(2,0) = chh[2];
						elem1_->LocalCoordinates(glo,loc);
//						std::cout<<glo(0,0)<<"\t"<<glo(1,0)<<"\t"<<glo(2,0)<<"\t";//blockkk or remove
						std::cout<<loc(0,0)<<"\t"<<loc(1,0)<<"\t"<<loc(2,0)<<"\n";
//                       std::cout<<chh[0]<<"\t"<<chh[1]<<"\t"<<chh[2]<<std::endl;
					}*/
			}

//should check whether this is sufficient or do we need to find the number of side in the element and
//use their orientation to get clockwise ordering
//                std::cout<<"parent side no = "<<parentSideno<<"\n";
//			std::cout<<"parentSideno = "<<parentSideno<<"corner = "<<corners[0][0]<<"\n";
//ParentSideno=1 is x=1 face and 3 is x=-1 face
#if 0
			if(parentSideno==1 && eqn_plane_[0]<0.0)
				clockwise_ = 1;
			if(parentSideno==3 && eqn_plane_[0]>0.0)
				clockwise_ = 1;
#endif
    if(parentSideno==1 && ordering=="cw")
      clockwise_ = 1;
    if(parentSideno==3 && ordering=="acw")
      clockwise_ = 1;

			/*if(corners[0][0]==1 && eqn_plane_[0]<0.0)
				clockwise_ = 1;
			if(corners[0][0]==-1 && eqn_plane_[0]>0.0)
				clockwise_ = 1;*/

  }
}

//computes x=a1+a2y+a3z from the plane equation
//equation of this form is used to replace x in the line integral
std::vector<double> GEO::CUT::FacetIntegration::compute_alpha(std::vector<double> eqn_plane,std::string intType)
{
  std::vector<double> alfa(3);
  double a = eqn_plane[0];
  double b = eqn_plane[1];
  double c = eqn_plane[2];
  double d = eqn_plane[3];


  if(intType=="x")
  {
    alfa[0] = d/a;
    alfa[1] = -1.0*b/a;
    alfa[2] = -1.0*c/a;
  }
  else if(intType=="y")
  {
        alfa[0] = d/b;
    alfa[1] = -1.0*c/b;
    alfa[2] = -1.0*a/b;
  }
  else if(intType=="z")
  {
    alfa[0] = d/c;
    alfa[1] = -1.0*a/c;
    alfa[2] = -1.0*b/c;
  }
  else
  {
    dserror("The facet integration type undefined");
    exit(1);
  }
//      std::cout<<alfa[0]<<"\t"<<alfa[1]<<"\t"<<alfa[2]<<std::endl;
  return alfa;
}

//perform integration over the facet
double GEO::CUT::FacetIntegration::integrate_facet()
{
//		std::cout<<"face"<<"\n";
		const std::vector<std::vector<double> > cornersLocal = face1_->CornerPointsLocal(elem1_,1);//blockkk !!!!remve the const
		/*if(bcellInt_==true)//blockkk
		{
			cornersLocal.clear();
			const std::vector<Point*> c1 = face1_->CornerPoints();
			for(std::vector<Point*>::const_iterator mm=c1.begin();mm!=c1.end();mm++)
			{
				Point* p1 = *mm;
				double x3[3];
				p1->Coordinates(x3);
				std::vector<double> dd;
				dd.push_back(x3[0]);
				dd.push_back(x3[1]);
				dd.push_back(x3[2]);
				cornersLocal.push_back(dd);
			}
		}//blockkk until this*/
        eqn_plane_ = equation_plane(cornersLocal);
//        std::cout<<"equation "<<eqn_plane_[0]<<"\t"<<eqn_plane_[1]<<"\t"<<eqn_plane_[2]<<"\t"<<eqn_plane_[3]<<std::endl;

// the face is in the x-y or in y-z plane which gives zero facet integral
        if(fabs(eqn_plane_[0])<0.0000001 && bcellInt_==false)
                return 0.0;
//x=0 plane which again do not contribute to facet integral
        if(fabs(eqn_plane_[1])<0.0000001 && fabs(eqn_plane_[2])<0.0000001 && fabs(eqn_plane_[3])<0.0000001 && bcellInt_==false)
                return 0.0;

/*        std::cout<<"facet points"<<"\n";
		for(unsigned i=0;i<cornersLocal.size();i++)
		{
			std::cout<<cornersLocal[i][0]<<"\t"<<cornersLocal[i][1]<<"\t"<<cornersLocal[i][2]<<"\n";
		}*/

        IsClockwise(eqn_plane_,cornersLocal);
//        std::cout<<"equation "<<eqn_plane_[0]<<"\t"<<eqn_plane_[1]<<"\t"<<eqn_plane_[2]<<"\t"<<eqn_plane_[3]<<std::endl;
  //      std::cout<<"clockwise = "<<clockwise_<<std::endl;


//        std::cout<<"alpha "<<alpha[0]<<"\t"<<alpha[1]<<"\t"<<alpha[2]<<std::endl;//blockkkk

        //integrating over each line of the facet
        double facet_integ = 0.0;
        if(bcellInt_==false)
        {
        	std::vector<double> alpha;
        	alpha = compute_alpha(eqn_plane_,"x");
			for(std::vector<std::vector<double> >::const_iterator k=cornersLocal.begin();k!=cornersLocal.end();k++)
			{
				const std::vector<double> coords1 = *k;
				std::vector<double> coords2;
				//for the last line the end point is the first point of the facet
				if(k!=(cornersLocal.end()-1))
						coords2 = *(k+1);
				else
						coords2= *(cornersLocal.begin());

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
				LineIntegration line1(coord1,coord2,inte_num_,alpha,bcellInt_);
				facet_integ += line1.integrate_line();
			}
        }
        else
        {
//to reduce the truncation error introduced during the projection of plane,
//the plane, along which the normal component is maximum, is chosen
        	std::string projType;
			if(fabs(eqn_plane_[0])<1e-8)
			{
				if(fabs(eqn_plane_[1])<1e-8)
					projType = "z";
				else if(fabs(eqn_plane_[2])<1e-8)
					projType = "y";
				else
				{
					if(fabs(eqn_plane_[1])>fabs(eqn_plane_[2]))
						projType = "y";
					else
						projType = "z";
				}
			}
			else if(fabs(eqn_plane_[1])<1e-8)
			{
				if(fabs(eqn_plane_[2])<1e-8)
					projType = "x";
				else
				{
					if(fabs(eqn_plane_[0])>fabs(eqn_plane_[2]))
						projType = "x";
					else
						projType = "z";
				}
			}
			else if(fabs(eqn_plane_[2])<1e-8)
			{
				if(fabs(eqn_plane_[1])>fabs(eqn_plane_[0]))
					projType = "y";
				else
					projType = "x";
			}
			else
			{
				if(fabs(eqn_plane_[0])>=fabs(eqn_plane_[1]) && fabs(eqn_plane_[0])>=fabs(eqn_plane_[2]))
					projType = "x";
				else if(fabs(eqn_plane_[1])>=fabs(eqn_plane_[2]) && fabs(eqn_plane_[1])>=fabs(eqn_plane_[0]))
					projType = "y";
				else
					projType = "z";
			}

			if(projType!="x" && projType!="y" && projType!="z")
			{
				dserror("projection plane is not defined");
				exit(1);
			}

			BoundaryFacetIntegration(cornersLocal,facet_integ,projType);

  //      	dserror("not yet implemented");
        }

//this condition results in negative normal for all the lines in the line integral
        if(clockwise_)
        {
                facet_integ = -1.0*facet_integ;
                for(int i=0;i<4;i++)
                       eqn_plane_[i] = -1.0*eqn_plane_[i];
        }
//	std::cout<<"facet_in"<<facet_integ<<"\n";
        return facet_integ;
}

void GEO::CUT::FacetIntegration::BoundaryFacetIntegration(const std::vector<std::vector<double> > cornersLocal,double &facet_integ,
		std::string intType)
{
	std::vector<double> alpha;
	double abs_normal=0.0,normalScale=0.0;//define this correspondingly
	for(unsigned i=0;i<3;i++)
		normalScale += eqn_plane_[i]*eqn_plane_[i];
	normalScale = sqrt(normalScale);
	for(std::vector<std::vector<double> >::const_iterator k=cornersLocal.begin();k!=cornersLocal.end();k++)
	{
		const std::vector<double> coords1 = *k;
		std::vector<double> coords2;
		//for the last line the end point is the first point of the facet
		if(k!=(cornersLocal.end()-1))
				coords2 = *(k+1);
		else
				coords2= *(cornersLocal.begin());

//		std::cout<<coords1[0]<<"\t"<<coords1[1]<<"\t"<<coords1[2]<<"\t"<<coords2[0]<<"\t"<<coords2[1]<<"\t"<<coords2[2]<<"\n";//blockkk

		std::vector<double> coord1,coord2;
		coord1.resize(2);
		coord2.resize(2);
		if(intType=="x")
		{
			if(k==cornersLocal.begin())
			{
				alpha = compute_alpha(eqn_plane_,"x");
				abs_normal = eqn_plane_[0]/normalScale;
			}
	//The facet is projected over y-z plane and then the integration is performed
	//so only y- and z-coordinates are passed to make the lines
	//[0]-- indicates y and [1] indicate z because we are now in y-z plane
			coord1[0] = coords1[1];
			coord1[1] = coords1[2];
			coord2[0] = coords2[1];
			coord2[1] = coords2[2];

			LineIntegration line1(coord1,coord2,inte_num_,alpha,bcellInt_);
			line1.set_integ_type("x");
			facet_integ += line1.integrate_line();
		}
		else if(intType=="y")
		{
			if(k==cornersLocal.begin())
			{
				alpha = compute_alpha(eqn_plane_,"y");
				abs_normal = eqn_plane_[1]/normalScale;
			}
	//The facet is projected over y-z plane and then the integration is performed
	//so only y- and z-coordinates are passed to make the lines
	//[0]-- indicates y and [1] indicate z because we are now in y-z plane
			coord1[0] = coords1[2];
			coord1[1] = coords1[0];
			coord2[0] = coords2[2];
			coord2[1] = coords2[0];
			LineIntegration line1(coord1,coord2,inte_num_,alpha,bcellInt_);
			line1.set_integ_type("y");
			facet_integ += line1.integrate_line();
		}

		else if(intType=="z")
		{
			if(k==cornersLocal.begin())
			{
				alpha = compute_alpha(eqn_plane_,"z");
				abs_normal = eqn_plane_[2]/normalScale;
			}
	//The facet is projected over y-z plane and then the integration is performed
	//so only y- and z-coordinates are passed to make the lines
	//[0]-- indicates y and [1] indicate z because we are now in y-z plane
			coord1[0] = coords1[0];
			coord1[1] = coords1[1];
			coord2[0] = coords2[0];
			coord2[1] = coords2[1];
			LineIntegration line1(coord1,coord2,inte_num_,alpha,bcellInt_);
			line1.set_integ_type("z");
			facet_integ += line1.integrate_line();
		}
		else
		{
			dserror("The facet integration type not supported");
			exit(1);
		}

//              std::cout<<coord1[0]<<"\t"<<coord1[1]<<"\t"<<coord2[0]<<"\t"<<coord2[1]<<std::endl;
//

	}
	//std::cout<<"abs normal = "<<abs_normal<<"\n";//blockkk
	facet_integ = facet_integ/abs_normal;
}
