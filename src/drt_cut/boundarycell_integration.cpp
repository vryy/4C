/*---------------------------------------------------------------------*/
/*!
\file boundarycell_integration.cpp

\brief used in boundary cell integration

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "boundarycell_integration.H"
#include "least_squares.H"
#include "base_vol.H"

#include "../drt_io/io_pstream.H"

/*--------------------------------------------------------------------------------------------*
  Generate quadrature rule for boundarycells
  Unlike volume integration the facets whose normal-x if zero, cannot be eliminated
  facet is projected into appropriate coordinate plane to obtain quadrature
*---------------------------------------------------------------------------------------------*/
Epetra_SerialDenseVector GEO::CUT::BoundarycellIntegration::GenerateBoundaryCellIntegrationRule()
{
  std::vector<std::vector<double> > corners1;
  bcell_->CornerPointsLocal(elem1_,corners1);
  Epetra_SerialDenseVector rhs_bcell_temp(num_func_);
  FacetIntegration faee1(bcell_,elem1_,position_,true,false);
  for(int fnc=1;fnc<=num_func_;fnc++)
  {
    faee1.set_integ_number(fnc);
    rhs_bcell_temp(fnc-1) = faee1.integrate_facet();
  }

  /*if(rhs_bcell_temp(0)<0.0)
    dserror("Negative area found in base function integration. Is ordering of vertices a problem?");*/

  Epetra_SerialDenseVector Bcellweights;

  int ptsEachLine = 14;//14 points gave min error for the test cases

  while(1)
  {
#if 0
    std::cout<<"pts on each line = "<<ptsEachLine<<"\n";
#endif

    std::vector<double> eqn = faee1.get_equation();
    DistributeBoundaryCellGaussPoints(eqn,corners1,BcellgausPts_,ptsEachLine);


    std::vector<std::vector<double> > moment_matbc(num_func_,std::vector<double>(BcellgausPts_.size()));
    momentFittingMatrix(moment_matbc,BcellgausPts_);


//if all the elements in a row of the moment fitting matrix are zero, then the row has to be deleted
//this ensures non-zero diagonal elements in the matrix
//it always occurs when the boundarycell is oriented along any of the coordinate planes
    std::vector<int> deleteRowNos;
    for(unsigned row=0;row<moment_matbc.size();row++)
    {
      bool deleteRow=true;
      for(unsigned col=0;col!=moment_matbc[0].size();col++)
      {
        if(fabs(moment_matbc[row][col])>1e-8)
        {
          deleteRow = false;
          break;
        }
      }
      if(deleteRow==true)
        deleteRowNos.push_back(row);
    }

    if(deleteRowNos.size()!=0)
    {
      for(unsigned row=0;row<deleteRowNos.size();row++)
      {
        int delno = deleteRowNos[row]-row;

        moment_matbc.erase(moment_matbc.begin()+delno);
      }
    }

    Epetra_SerialDenseVector rhs_bcell(num_func_-deleteRowNos.size());
    if(deleteRowNos.size()==0)
    {
      for(int m=0;m<rhs_bcell_temp.Length();m++)
        rhs_bcell(m) = rhs_bcell_temp(m);
    }
    else
    {
      int rowno=0,rhsno=0;
      for(int m=0;m<rhs_bcell_temp.Length();m++)
      {
        int deleteNo = deleteRowNos[rowno];
        if(m==deleteNo)
        {
          rowno++;
          continue;
        }
        else
        {
          rhs_bcell(rhsno) = rhs_bcell_temp(m);
          rhsno++;
        }
      }
    }

    LeastSquares least(moment_matbc,rhs_bcell);
    Bcellweights.Size(moment_matbc[0].size());
    Bcellweights = least.linear_least_square();

    Epetra_SerialDenseVector err(num_func_);
    for(int i=0;i<num_func_;i++)
    {
      err(i) = 0.0;
      for(unsigned j=0;j<BcellgausPts_.size();j++)
      {
        err(i) += Bcellweights(j)*base_function(BcellgausPts_[j],i+1);
      }

      if(fabs(rhs_bcell_temp(i))>1e-8)
        err(i) = (err(i)-rhs_bcell_temp(i))/rhs_bcell_temp(i);
      else
        err(i) = err(i)-rhs_bcell_temp(i);
    }

    double maxerr = err.InfNorm();
#if 0
    std::cout<<"numpts = "<<ptsEachLine<<"\tmax error = "<<maxerr<<"\n";
#endif

    if(maxerr<1e-10 || ptsEachLine>25)
      break;
    else
    {
      ptsEachLine++;
      BcellgausPts_.clear();
    }
  }

#ifdef DEBUGCUTLIBRARY
    BcellGaussPointGmsh(BcellgausPts_,corners1);
#endif

  return Bcellweights;
}

/*--------------------------------------------------------------------------------------------*
   The arbitrarily oriented plane is first projected into one of the coordinate planes
   the gauss points are generated in the coordinate plane - now it has only two coordinates
   value of third coordinate can be calculated from the equation of arbitrary plane
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::BoundarycellIntegration::DistributeBoundaryCellGaussPoints(std::vector<double> eqn, std::vector<std::vector<double> > corners,
                    std::vector<std::vector<double> >& bcGausspts,int ptNos)
{
  std::vector<double> co1(corners.size()), co2(corners.size());

  std::string projType;

  int dim1=0,dim2=1,dim3=2;

//to reduce the truncation error introduced during the projection of plane,
//the plane along which the normal component is maximum, is chosen

//initial check whether the cell is already in the coordinate plane
  if(fabs(eqn[0])<1e-8)
  {
    if(fabs(eqn[1])<1e-8)
      projType = "xy";
    else if(fabs(eqn[2])<1e-8)
      projType = "zx";
    else
    {
      if(fabs(eqn[1])>fabs(eqn[2]))
        projType = "zx";
      else
        projType = "xy";
    }
  }
  else if(fabs(eqn[1])<1e-8)
  {
    if(fabs(eqn[2])<1e-8)
      projType = "yz";
    else
    {
      if(fabs(eqn[0])>fabs(eqn[2]))
        projType = "yz";
      else
        projType = "xy";
    }
  }
  else if(fabs(eqn[2])<1e-8)
  {
    if(fabs(eqn[1])>fabs(eqn[0]))
      projType = "zx";
    else
      projType = "yz";
  }
// if cell is not in the coordinate plane, project it along the plane of maximum normal
  else
  {
    if(fabs(eqn[0])>=fabs(eqn[1]) && fabs(eqn[0])>=fabs(eqn[2]))
      projType = "yz";
    else if(fabs(eqn[1])>=fabs(eqn[2]) && fabs(eqn[1])>=fabs(eqn[0]))
      projType = "zx";
    else
      projType = "xy";
  }

  if(projType=="xy"){
    dim1=0,dim2=1,dim3=2;}
  else if(projType=="yz"){
    dim1=1,dim2=2,dim3=0;}
  else if(projType=="zx"){
    dim1=2,dim2=0,dim3=1;}
  else
  {
    dserror("the projection plane is not set");
    exit(1);
  }

  for(unsigned i=0;i<corners.size();i++)
  {
    co1[i] = corners[i][dim1];
    co2[i] = corners[i][dim2];
  }

  //store the equation of each line in the form of ax+by=c
  std::vector<std::vector<double> > eqnLines(co1.size());
  for(unsigned i=0;i<co1.size();i++)
  {
    double x1 = co1[i];
    double x2 = co1[(i+1)%co1.size()];
    double y1 = co2[i];
    double y2 = co2[(i+1)%co2.size()];
    std::vector<double> line(3);
    if(fabs(x2-x1)<0.00000001)
    {
      line[0] = 1.0;
      line[1] = 0.0;
      line[2] = x1;
    }
    else if(fabs(y2-y1)<0.00000001)
    {
      line[0] = 0.0;
      line[1] = 1.0;
      line[2] = y1;
    }
    else
    {
      line[0] = 1.0;
      line[1] = (x1-x2)/(y2-y1);
      line[2] = x1+line[1]*y1;
    }
    eqnLines[i] = line;
  }

//  double max1 = *(std::max_element(co1.begin(), co1.end()));
//  double min1 = *(std::min_element(co1.begin(), co1.end()));
  double max2 = *(std::max_element(co2.begin(), co2.end()));
  double min2 = *(std::min_element(co2.begin(), co2.end()));

  int gausno=0;//number of gauss points in each line
  bcGausspts.resize(ptNos*ptNos,std::vector<double> (3));
  double interval = (max2-min2)/ptNos;
  std::map<double,int> intersections;
  for(int i=0;i<ptNos;i++)
  {
    double val2 = min2+(i+0.5)*interval;

    intersections.clear();
    for(unsigned j=0;j<co2.size();j++)
    {
      double y1 = co2[j];
      double y2 = co2[(j+1)%co2.size()];
      if(fabs(y1-y2)<1e-8)
        continue;
      if((y1<=val2 && y2>=val2) || (y1>=val2 && y2<=val2))
      {
        double val1 = (eqnLines[j][2]-eqnLines[j][1]*val2)/eqnLines[j][0];
        intersections[val1] = j;
      }
    }
    if(intersections.size()==2)
    {
      int count=0;
      double x1[2];
      for(std::map<double,int>::iterator j=intersections.begin();j!=intersections.end();j++)
      {
        x1[count] = j->first;
        count++;
      }

      int actNo = ptNos;
      /*if(fabs(x1[1]-x1[0])<0.05*(max1-min1))
        actNo = 5;
      else if(fabs(x1[1]-x1[0])<0.1*(max1-min1))
        actNo = 10;*/
      double intervalx = (x1[1]-x1[0])/actNo;
      for(int j=0;j<actNo;j++)
      {
        bcGausspts[gausno][dim1] = x1[0]+(j+0.5)*intervalx;
        bcGausspts[gausno][dim2] = val2;
        gausno++;
      }
    }
    else //more than 2 intersections for bcell is unlikely (???). Not tested enough
    {
      int count=0, numcut = intersections.size();;
      std::vector<double> x1(numcut);
      for(std::map<double,int>::iterator j=intersections.begin();j!=intersections.end();j++)
      {
        x1[count] = j->first;
        count++;
      }

      int actNo=ptNos/numcut+1;
      for(int i=0;i<numcut/2;i++)
      {
        double x11 = x1[i*2];
        double x21 = x1[i*2+1];
        double intervalx = (x21-x11)/actNo;
        for(int j=0;j<actNo;j++)
        {
          bcGausspts[gausno][dim1] = x11+(j+0.5)*intervalx;
          bcGausspts[gausno][dim2] = val2;
          gausno++;
        }
      }
    }
  }

//  bcGausspts.erase (bcGausspts.begin()+gausno,bcGausspts.end());

//calculation of third coordinate from equation of plane
//every point in the bcell should satisfy the plane equation
  for(unsigned i=0;i<bcGausspts.size();i++)
  {
    if(projType=="xy")
      bcGausspts[i][dim3] = (eqn[3]-eqn[1]*bcGausspts[i][dim2]-eqn[0]*bcGausspts[i][dim1])/eqn[2];
    else if(projType=="yz")
      bcGausspts[i][dim3] = (eqn[3]-eqn[1]*bcGausspts[i][dim1]-eqn[2]*bcGausspts[i][dim2])/eqn[0];
    else if(projType=="zx")
      bcGausspts[i][dim3] = (eqn[3]-eqn[0]*bcGausspts[i][dim2]-eqn[2]*bcGausspts[i][dim1])/eqn[1];
    else
    {
      dserror("projection type not assigned");
      exit(1);
    }
  }
}

/*--------------------------------------------------------------------------------------------*
                                   form the moment fitting matrix
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::BoundarycellIntegration::momentFittingMatrix(std::vector<std::vector<double> >&mom,std::vector<std::vector<double> >gauspts)
{
    for(int i=0;i<num_func_;i++)
    {
        int k=0;
        for(std::vector<std::vector<double> >::iterator j=gauspts.begin();j!=gauspts.end();j++)
        {
            std::vector<double> cordi = *j;
            mom[i][k] = base_function(cordi,i+1);
            k++;
        }
    }
}

/*--------------------------------------------------------------------------------------------*
  The geometry of boundarycell and the location of Gaussian points are written in GMSH output
  file for visualization and for debugging
 *--------------------------------------------------------------------------------------------*/
void GEO::CUT::BoundarycellIntegration::BcellGaussPointGmsh(const std::vector<std::vector<double> > bcGausspts,
    const std::vector<std::vector<double> > corners)
{
  std::string filename="bcell";
  std::ofstream file;

  static int bcellno = 0;
  bcellno++;
  std::stringstream out;
  out <<"bcell"<<bcellno<<".pos";
  filename = out.str();
  file.open(filename.c_str());

  int pointno=1,point_end,lineno=1;
  for(unsigned i=0;i<corners.size();i++)
  {
     file<<"Point("<<pointno<<")={"<<corners[i][0]<<","<<corners[i][1]<<","<<corners[i][2]<<","<<"1"<<"};"<<std::endl;
     pointno++;
  }
  point_end = pointno;
  for(int i=1;i!=point_end;i++)
  {
     if(i!=point_end-1)
         file<<"Line("<<lineno<<")={"<<i<<","<<i+1<<"};"<<std::endl;
     else
         file<<"Line("<<lineno<<")={"<<i<<","<<1<<"};"<<std::endl;
     lineno++;
  }

  for(unsigned i=0;i<bcGausspts.size();i++)
  {
       file<<"Point("<<pointno<<")={"<<bcGausspts[i][0]<<","<<bcGausspts[i][1]<<","<<bcGausspts[i][2]<<","<<"1"<<"};"<<std::endl;
       pointno++;
  }
  file.close();
}
