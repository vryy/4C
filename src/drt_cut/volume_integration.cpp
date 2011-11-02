#include<iostream>
#include<fstream>
#include <cmath>
#include<algorithm>
#include "volume_integration.H"
#include "base_vol.H"
#include "least_squares.H"
#include "cut_boundingbox.H"
#include "cut_options.H"

using namespace std;

//compute the rhs of the moment fitting equations
//Integration of base functions take place inside this
Epetra_SerialDenseVector GEO::CUT::VolumeIntegration::compute_rhs_moment()
{
     Epetra_SerialDenseVector rhs_mom(num_func_);

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
            if(fnc==1)
            {
                std::vector<double> eqn = faee1.get_equation();
                eqn_facets_.push_back(eqn);
//		std::cout<<"eqn_plane"<<eqn[0]<<"\t"<<eqn[1]<<"\t"<<eqn[2]<<"\t"<<eqn[3]<<"\n"; //blockkk or remove
            }
        }
        rhs_mom(fnc-1) = mome;
//      std::cout<<mome<<std::endl;
//
    }
    return rhs_mom;
}

//compute the gaussian points of the volumecell with "numeach" points in each 3-directions
//numeach should be more than 1
bool GEO::CUT::VolumeIntegration::compute_Gaussian_points(int numeach)
{
    
/*  for(std::vector<FacetIntegration>::iterator i=facet_list_.begin();
        i!=facet_list_.end();i++)
    {
        FacetIntegration facee = *i;
    }*/

/*        for(unsigned i=0;i<eqn_facets_.size();i++)      
               std::cout<<eqn_facets_[i][0]<<"\t"<<eqn_facets_[i][1]<<"\t"<<eqn_facets_[i][2]<<"\t"<<eqn_facets_[i][3]<<std::endl;*/

    /*********geneate Gaussian points by checking the intersection of an arbitrary line
     ********* with the sides of the volumecell**************************************/
//construct a bounding box within which the Gaussian points are distributed and checked whether
//they are inside the volume cell or not
    BoundingBox box1(*volcell_,elem1_);
    double minn[3],maxx[3];
    minn[0] = box1.minx();
    maxx[0] = box1.maxx();
    minn[1] = box1.miny();
    maxx[1] = box1.maxy();
    minn[2] = box1.minz();
    maxx[2] = box1.maxz();
    /*std::cout<<minn[0]<<"\t"<<maxx[0]<<"\t";
    std::cout<<minn[1]<<"\t"<<maxx[1]<<"\t";
    std::cout<<minn[2]<<"\t"<<maxx[2]<<std::endl;*/
 
    vector<vector<double> > zcoord;
    vector<vector<double> > ycoord;
    get_zcoordinates(zcoord,ycoord);

//min z plane which contains significant area
    double zmin=minn[2]+0.01*(maxx[2]-minn[2]),zmax=maxx[2]-0.01*(maxx[2]-minn[2]);
    while(1)
    {
        bool area=false;
        vector<vector<double> > InPlane;
        area = IsContainArea(minn,maxx,zmin,InPlane,zcoord,ycoord,0.01,numeach);
        if(area)
        {
              gaus_pts_.insert(gaus_pts_.end(),InPlane.begin(),InPlane.end());
              break;
        }
        zmin += 0.01*(maxx[2]-minn[2]);
        if((zmax-zmin)<0.01*(maxx[2]-minn[2]))
                break;
    }
//    std::cout<<"intermediate gauss points = "<<gaus_pts_.size()<<"\n";
//max z plane which contains significant area
    while(1)
    {
        bool area=false;
        vector<vector<double> > InPlane;
        area = IsContainArea(minn,maxx,zmax,InPlane,zcoord,ycoord,0.01,numeach);
        if(area)
        {
              gaus_pts_.insert(gaus_pts_.end(),InPlane.begin(),InPlane.end());
              break;
        }
        zmax -= 0.01*(maxx[2]-minn[2]);
        if((zmax-zmin)<0.01*(maxx[2]-minn[2]))
                break;
    }
    bool wei=true;
    if((zmax-zmin)<0.01*(maxx[2]-minn[2]))
    {
//check for very thin volumecells
//works fine but time taken is slightly higher
        double zmin=minn[2]+0.005*(maxx[2]-minn[2]),zmax=maxx[2]-0.005*(maxx[2]-minn[2]);
        while(1)
        {
			bool area=false;
			vector<vector<double> > InPlane;
			area = IsContainArea(minn,maxx,zmin,InPlane,zcoord,ycoord,0.001,numeach);
			if(area)
			{
				  gaus_pts_.insert(gaus_pts_.end(),InPlane.begin(),InPlane.end());
				  break;
			}
			zmin += 0.001*(maxx[2]-minn[2]);
			if((zmax-zmin)<0.001*(maxx[2]-minn[2]))
					break;
        }
        while(1)
        {
			bool area=false;
			vector<vector<double> > InPlane;
			area = IsContainArea(minn,maxx,zmax,InPlane,zcoord,ycoord,0.001,numeach);
			if(area)
			{
				  gaus_pts_.insert(gaus_pts_.end(),InPlane.begin(),InPlane.end());
				  break;
			}
			zmax -= 0.001*(maxx[2]-minn[2]);
			if((zmax-zmin)<0.001*(maxx[2]-minn[2]))
					break;
        }

    //    cout<<"very small area"<<endl;
        if(gaus_pts_.size()==0)
                wei = false;
        else
                wei = true;
/*        if(gaus_pts_.size()==0)
                wei = false;*/
        cout<<"number of Gauss points"<<gaus_pts_.size()<<endl;//blockkk
        return wei;
    }
//find z-planes in between zmin and zmax to generate Gauss points
    else
    {
        int num = numeach;   
        double dzz = (zmax-zmin)/(num-1);
        std::vector<double> zplane;
        for(int i=0;i<num;i++)
                zplane.push_back(zmin+i*dzz);
         bool previous = true;
         for(int i=1;i<num-1;i++)
         {
                bool area = false;
                vector<vector<double> > InPlane;
                area = IsContainArea(minn,maxx,zplane[i],InPlane,zcoord,ycoord,0.01,numeach);
                if(area)
                {
                        gaus_pts_.insert(gaus_pts_.end(),InPlane.begin(),InPlane.end());
                        InPlane.clear();
                        previous = true;
                        continue;
                }
//if the considered z-plane does not contain significant area, the interval is subdivided to check whether
//any plane lies in between
                else
                {
                        if(previous==false)
                                continue;
                        double dzzm = dzz;
                        for(int k=0;k<4;k++)
                        {
                                dzzm += 0.5*dzzm;
                                bool innerarea = false;
                                double zz = zplane[i]-dzzm;
                                innerarea = IsContainArea(minn,maxx,zz,InPlane,zcoord,ycoord,0.01,numeach);
                                if(innerarea)
                                {
                                        gaus_pts_.insert(gaus_pts_.end(),InPlane.begin(),InPlane.end());
                                        InPlane.clear();
                                        break; 
                                }
                        }
                        previous = false;
                }
        }
    }
    cout<<"number of Gauss points"<<gaus_pts_.size()<<endl;//blockkk
    return wei;
     
/*************** probably need to remove this*****************************************/
//create a bounding box over the volume and generate gaussian points over that volume   
/*    const plain_facet_set & facete = volcell_->Facets();
    std::vector<double> x1,y1,z1;
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
        Facet *fac = *i;
        const std::vector<Point*> & corners = fac->CornerPoints();
//        std::cout<<"facet"<<std::endl;
        for(std::vector<Point*>::const_iterator k=corners.begin();k!=corners.end();k++)
        {
            const Point* po = *k;
            const double * coords = po->X();
            x1.push_back(coords[0]);
            y1.push_back(coords[1]);
            z1.push_back(coords[2]);
//            std::cout<<coords[0]<<"\t"<<coords[1]<<"\t"<<coords[2]<<std::endl;
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
    }*/
//    std::cout<<"size"<<gaus_pts_.size()<<std::endl;
/*  for(int i=0;i<27;i++)
    {
        std::cout<<gau_pts[i][0]<<"\t"<<gau_pts[i][1]<<"\t"<<gau_pts[i][2]<<std::endl;
    }*/
/*******************until this*******************************************/
}

//store the z- and y-coordinates of the all corner points which will be used to find whether the intersection
//point lies inside the volume or not
void GEO::CUT::VolumeIntegration::get_zcoordinates(vector<vector<double> >& zcoord,
                vector<vector<double> >& ycoord)
{
    const plain_facet_set & facete = volcell_->Facets();
    int faceno=0;
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
            vector<double> thisplane1,thisplane2;
//since the imaginary line never intersect this plane
            if(fabs(eqn_facets_[faceno][0])<0.00000001)
            {
                    thisplane1.push_back(0.0);
                    zcoord.push_back(thisplane1);
                    ycoord.push_back(thisplane1);
                    faceno++;
                    continue;
            }
            Facet *face1 = *i;
            const std::vector<vector<double> > corLocal = face1->CornerPointsLocal(elem1_,0);
            for(std::vector<std::vector<double> >::const_iterator k=corLocal.begin();k!=corLocal.end();k++)
            {
				std::vector<double> coords1 = *k;

                thisplane1.push_back(coords1[1]);
                thisplane2.push_back(coords1[2]);
            }
            ycoord.push_back(thisplane1);
            zcoord.push_back(thisplane2);
            faceno++;
    }

    /*std::cout<<"zcoord\n";
    for(vector<vector<double> >::iterator i=zcoord.begin();i!=zcoord.end();i++)
    {
    	vector<double> zz=*i;
    	for(vector<double>::iterator j=zz.begin();j!=zz.end();j++)
    	{
    		std::cout<<*j<<"\t";
    	}
    	std::cout<<"\n";
    }
    std::cout<<"ycoord\n";
	for(vector<vector<double> >::iterator i=ycoord.begin();i!=ycoord.end();i++)
	{
		vector<double> zz=*i;
		for(vector<double>::iterator j=zz.begin();j!=zz.end();j++)
		{
			std::cout<<*j<<"\t";
		}
		std::cout<<"\n";
	}*/
}

//check whether the considered imaginary line intersect any of the facets
//if so generate gauss points on that line
bool GEO::CUT::VolumeIntegration::IsIntersect(double *pt, double *mini, double *maxi,vector<vector<double> >& linePts,
                vector<vector<double> >zcoord,vector<vector<double> >ycoord,double toler,int numeach)
{
        bool intersect = false;
        vector<int> planeinter,InterFaces;

//stores all the facets which is not in x-y or x-z
//only in the remaining planes a horizontal line possibly intersect
        for(unsigned i=0;i<eqn_facets_.size();i++)
        {
                if(fabs(eqn_facets_[i][0])>0.0000001)
                        planeinter.push_back(i);
        }

//stores the facets which are actually cut by the line
        for(unsigned i=0;i<planeinter.size();i++)
        {
                int faceno = planeinter[i];
                vector<double>planez = zcoord[faceno];
                vector<double>planey = ycoord[faceno];
                int cutno=0;
//check whether the intersection point lies inside the facet area
                cutno = pnpoly(planez.size(), planey, planez, pt[1], pt[2]);
                if(cutno==1)
                        InterFaces.push_back(faceno);
        }


//        cout<<InterFaces.size()<<endl;
        if(InterFaces.size()<2)
        {
                intersect = false;
                return intersect;
        }
//simple geometries result in two intersections
        if(InterFaces.size()==2)
        {
                int fa1 = InterFaces[0];
                int fa2 = InterFaces[1];

                double int1 = (eqn_facets_[fa1][3]-eqn_facets_[fa1][1]*pt[1]-eqn_facets_[fa1][2]*pt[2])/eqn_facets_[fa1][0];
                double int2 = (eqn_facets_[fa2][3]-eqn_facets_[fa2][1]*pt[1]-eqn_facets_[fa2][2]*pt[2])/eqn_facets_[fa2][0];

                if(fabs(int2-int1)<toler*(maxi[0]-mini[0]))
                {
                        intersect = false;
                        return intersect;
                }
                else
                        intersect = true;

 //               cout<<int1<<"\t"<<int2<<endl;
                vector<double> inter1,inter2;
                inter1.push_back(int1);
                inter1.push_back(pt[1]);
                inter1.push_back(pt[2]);
                inter2.push_back(int2);
                inter2.push_back(pt[1]);
                inter2.push_back(pt[2]);
                if(fabs(inter2[0]-inter1[0])<0.025*(maxi[0]-mini[0]))
                {
                        vector<double> middle(3);
                        middle[0] = 0.5*(inter2[0]+inter1[0]);
                        middle[1] = inter2[1];
                        middle[2] = inter2[2];
                        linePts.push_back(middle);
                        intersect = true;
                }
                else if(fabs(inter2[0]-inter1[0])<0.05*(maxi[0]-mini[0]))
                {
                        OnLine(inter1,inter2,linePts,2/*numeach/2+1*/);//blockkk or remove
                        intersect = true;
                }
                else
                        OnLine(inter1,inter2,linePts,numeach);
                return intersect;
        } 
        else
        {
//map is useful since we need to arrange the elements from minimum x-cut value
                cout<<"in more than two cuts"<<endl;
                map<vector<double>,int> interPoints;
                for(vector<int>::iterator i=InterFaces.begin();i!=InterFaces.end();i++)
                {
                        int fa1 = *i;
                        double int1 = (eqn_facets_[fa1][3]-eqn_facets_[fa1][1]*pt[1]-eqn_facets_[fa1][2]*pt[2])/eqn_facets_[fa1][0];
                        vector<double> inter1;
                        inter1.push_back(int1);
                        inter1.push_back(pt[1]);
                        inter1.push_back(pt[2]);
                        interPoints[inter1] = fa1;
                }

/*                for(map<vector<double>,int>::iterator i=interPoints.begin();i!=interPoints.end();i++)
                {
                        vector<double> mmm = i->first;
                        int kkk = i->second;
                        cout<<mmm[0]<<"\t"<<mmm[1]<<"\t"<<mmm[2]<<"\t"<<kkk<<"\n";
                }*/

                const plain_facet_set & facete = volcell_->Facets();
                plain_facet_set::const_iterator f=facete.begin();
                unsigned cut_count=0;
                for(map<vector<double>,int>::iterator i=interPoints.begin();
                                i!=interPoints.end();i++)
                {
                        bool ptsInside = false;
                        unsigned face1 = i->second;
                        Facet *facet1 = *(f+face1);
                        vector<double> inter1 = i->first;
                        i++;
                        unsigned face2 = i->second;
                        Facet *facet2 = *(f+face2);
                        vector<double> inter2 = i->first;
                        i--;
//Among two consequetive facets which are cut by the line, one must be the cut surface
                        if(facet1->OnCutSide())
                        {
                                if(eqn_facets_[face1][0]<0.0)
                                        ptsInside=true;
           //                     cout<<"face1"<<endl;
             //                   cout<<inter1[0]<<"\t"<<inter2[0]<<"\t"<<eqn_facets_[face1][0]<<endl;
                        }
                        else if(facet2->OnCutSide())
                        {
                                if(eqn_facets_[face2][0]>0.0)
                                        ptsInside=true;
               //                 cout<<"face2"<<endl;
                 //               cout<<inter1[0]<<"\t"<<inter2[0]<<"\t"<<eqn_facets_[face2][0]<<endl;
                        }
                        else
                                cout<<"The assumption that one must be a cut surface is false"<<endl;
                        if(ptsInside)
                        {
                                if(interPoints.size()==2)//if intersection is on the triangulated line, results only in two face intersections
                                {
                                        OnLine(inter1,inter2,linePts,numeach);
                                        intersect = true;
                                        break;
                                }
                                if((inter2[0]-inter1[0])<toler*(maxi[0]-mini[0]))
                                {
                                        vector<double> middle(3);
                                        middle[0] = 0.5*(inter2[0]+inter1[0]);
                                        middle[1] = inter2[1];
                                        middle[2] = inter2[2];
                                        linePts.push_back(middle);
                                        intersect = true;
                                }
                                else if((inter2[0]-inter1[0])<0.05*(maxi[0]-mini[0]))
                                {
                                        OnLine(inter1,inter2,linePts,2/*numeach/3+1*/);//blockkk or remove
                                        intersect = true;
                                }
                                else
                                {
                                        OnLine(inter1,inter2,linePts,numeach/2+1);
                                        intersect = true;
                                }
                        }
                        cut_count++;
                        if(cut_count==(interPoints.size()-1))
                                break;
                }
        }
        return intersect;
}

//find whether the point is inside the polygon
int GEO::CUT::VolumeIntegration::pnpoly(int npol, vector<double>xp, vector<double>yp, double x, double y)
{
      int i, j, c = 0;
      for (i = 0, j = npol-1; i < npol; j = i++) {
        if ((((yp[i]<=y) && (y<yp[j])) ||
             ((yp[j]<=y) && (y<yp[i]))) &&
            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
          c = !c;
      }
      return c;
}

//check whether at a particular z-plane significant area exists, and if so returns the Gauss points
//in the particular plane
bool GEO::CUT::VolumeIntegration::IsContainArea(double minn[3],double maxx[3], double&zmin,vector<vector<double> > &pts,
                vector<vector<double> >zcoord,vector<vector<double> >ycoord,double toler,int numeach)
{
     bool isArea=true;
     double dx1[3];
     dx1[0]=maxx[0]-minn[0],dx1[1]=maxx[1]-minn[1],dx1[2]=maxx[2]-minn[2];
     double xmin=minn[0]+toler*dx1[0],ymin=minn[1],ymax=maxx[1];

//checl for the lowest line
     while(1)
     {
		 vector<vector<double> > linePts;
		 ymin += toler*dx1[1];
		 bool intersec = false;
		 double a[] = {xmin,ymin,zmin};
		 intersec = IsIntersect(a, minn, maxx, linePts,zcoord,ycoord,toler,numeach);

// the area of the solid volume in this particular z-plane is almost negligible
// if the area of volumecell is less than 1% of bounding box
		 if((ymax-ymin)<toler*dx1[1])
		 {
			  isArea = false;
			  return isArea;
		 }
		 if(intersec)
		 {
//                 cout<<ymin<<"\t"<<ymax<<endl;
			  pts.insert(pts.end(),linePts.begin(),linePts.end());
			  break;
		 }
     }
//   cout<<"line1"<<endl; 
//check for the topmost line
     while(1)
     {
		 vector<vector<double> > linePts;
		 ymax -= toler*dx1[1];
//               cout<<ymax<<endl;
		 bool intersec = false;
		 double a[] = {xmin,ymax,zmin};
		 intersec = IsIntersect(a, minn, maxx, linePts,zcoord,ycoord,toler,numeach);
		 if((ymax-ymin)<toler*dx1[1])
		 {
			 isArea = false;
			 return isArea;
		 }
		 if(intersec)
		 {
			 pts.insert(pts.end(),linePts.begin(),linePts.end());
			 break;
		 }
     }
//     cout<<"in this plane"<<ymin<<"\t"<<ymax<<endl;//blockkk or remove

//generate points in between the topmost and lowest line
    int num = numeach;   
    double dyy = (ymax-ymin)/(num-1);
    std::vector<double> yplane;
    for(int i=0;i<num;i++)
                yplane.push_back(ymin+i*dyy);
    bool previous = true;
    for(int i=1;i<num-1;i++)
    {
		double a[] = {xmin,yplane[i],zmin};
		bool intersec = false;
		vector<vector<double> > linePts;
		intersec = IsIntersect(a, minn, maxx, linePts,zcoord,ycoord,toler,numeach);
		if(intersec)
		{
				pts.insert(pts.end(),linePts.begin(),linePts.end());
				linePts.clear();
				previous = true;
				continue;
		}
		else
		{
				if(previous==false)
						continue;
				double dyym = dyy;
				for(int k=0;k<4;k++)
				{
						dyym += 0.5*dyym;
						double a[] = {xmin,yplane[i]-dyym,zmin};
						bool innerintersec = false;
						vector<vector<double> > linePts;
						innerintersec = IsIntersect(a, minn, maxx, linePts,zcoord,ycoord,toler,numeach);

						 if(innerintersec)
						 {
								pts.insert(pts.end(),linePts.begin(),linePts.end());
								linePts.clear();
								break;
						 }
				}
				previous = false;
		}
    }
    return isArea;  
}

//generates specified number of points (num) on the imaginary line
void GEO::CUT::VolumeIntegration::OnLine(vector<double>inter1,vector<double>inter2,
                vector<vector<double> >&linePts,int num)
{
		vector<double> left,right;
		if(inter1[0]<inter2[0])
		{
				left = inter1;
				right = inter2;
		}
		else
		{
				left = inter2;
				right = inter1;
		}
		double xlen = right[0]-left[0];
		inter1[0] = left[0]+0.05*xlen;
		inter2[0] = right[0]-0.05*xlen;
		double xdiff = (inter2[0]-inter1[0])/(num-1);
		for(int i=0;i<num;i++)
		{
				vector<double> temp(3);
				temp[0] = inter1[0]+i*xdiff;
				temp[1] = inter1[1];
				temp[2] = inter1[2];
				linePts.push_back(temp);
		}
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
Epetra_SerialDenseVector GEO::CUT::VolumeIntegration::compute_weights()
{

    Epetra_SerialDenseVector rhs_moment(num_func_);
    rhs_moment = compute_rhs_moment();
    
/*  for(std::vector<double>::iterator i=rhs_moment.begin();i!=rhs_moment.end();i++)
        std::cout<<*i<<std::endl;*/
    bool wei;
//we should ask for more than 1 point in each direction
    wei = compute_Gaussian_points(6);
/*  for(int i=0;i<27;i++)
    {
        std::cout<<gaus_pts_[i][0]<<"\t"<<gaus_pts_[i][1]<<"\t"<<gaus_pts_[i][2]<<std::endl;
    }*/

    Epetra_SerialDenseVector weights;
    if(wei)
    {
        std::vector<std::vector<double> > moment_matrix(num_func_,std::vector<double>(gaus_pts_.size()));
        moment_fitting_matrix(moment_matrix);

        LeastSquares least(moment_matrix,rhs_moment);
        weights.Size(moment_matrix[0].size());
        weights = least.linear_least_square();
    }
    else
    {
        gaus_pts_.clear();
        vector<double> zer(3);
        zer[0]=0.0;zer[1]=0.0;zer[2]=0.0;
        gaus_pts_.push_back(zer);
        weights.Size(1);
        weights(0) = 0.0;
    }

/*    for(int j=0;j<num_func_;j++)
    {
    double chek = 0.0;
    for(int i=0;i<gaus_pts_.size();i++)
    {
        std::vector<double> coorrr;
        coorrr.push_back(gaus_pts_[i][0]);
        coorrr.push_back(gaus_pts_[i][1]);
        coorrr.push_back(gaus_pts_[i][2]);
        chek += weights(i)*base_function(coorrr,j+1);
    }
    std::cout<<"check"<<(rhs_moment(j)-chek)<<std::endl;
//    std::cout<<"value"<<chek<<std::endl;
    }*/

/*    double chek=0.0;
    for(unsigned i=0;i<weights.size();i++)
    {
         chek += weights(i)*(pow(gaus_pts_[i][0],4)+pow(gaus_pts_[i][1],4)+5.0);     
    //  chek += weights(i)*(gaus_pts_[i][0]*gaus_pts_[i][0]*gaus_pts_[i][0]+gaus_pts_[i][1]*gaus_pts_[i][1]*gaus_pts_[i][1]+5.0);       
//        chek += pow(gaus_pts_[i][0],4)*weights(i];
    }
    std::cout<<scientific<<"check"<<chek<<std::endl;
    std::cout<<scientific<<"error"<<chek-(rhs_moment(20)+rhs_moment(30)+5*rhs_moment(0))<<std::endl;*/

/*  for(int i=0;i<weights.size();i++)
        std::cout<<"wei"<<weights(i)<<std::endl;*/

    std::cout<<"volume = "<<rhs_moment(0)<<"\t"<<rhs_moment(1)<<"\t"<<rhs_moment(4)<<std::endl;

#ifdef DEBUGCUTLIBRARY
    GaussPointGmsh();
#endif

    return weights;
}

void GEO::CUT::VolumeIntegration::GaussPointGmsh()
{
        volcell_->DumpGmshGaussPoints(gaus_pts_);
}
