/*---------------------------------------------------------------------*/
/*! \file

\brief Integrates base functions over volume, distribute Gauss points and solve moment fitting
equations

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/

#include <fstream>
#include <cmath>
#include <algorithm>
#include "volume_integration.H"
#include "base_vol.H"
#include "least_squares.H"
#include "cut_boundingbox.H"
#include "cut_options.H"
#include "cut_enum.H"
#include <Epetra_MultiVector.h>

#include "../drt_lib/drt_colors.H"

/*--------------------------------------------------------------------------*
         compute the rhs of the moment fitting equations
         Integration of base functions take place inside this
*---------------------------------------------------------------------------*/
Epetra_SerialDenseVector GEO::CUT::VolumeIntegration::compute_rhs_moment()
{
  Epetra_SerialDenseVector rhs_mom(num_func_);

  const plain_facet_set &facete = volcell_->Facets();

  // integrate each base function over this volumecell
  for (int fnc = 1; fnc <= num_func_; fnc++)
  {
    double mome = 0.0;

    for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
    {
      Facet *fe = *i;

      // an equivalent function corresponding to every base function will be
      // integrated over facet of vcell
      FacetIntegration faee1(fe, elem1_, position_, false, false);
      faee1.set_integ_number(fnc);
      mome += faee1.integrate_facet();

      if (fnc == 1)
      {
        std::vector<double> eqn = faee1.get_equation();
        eqn_facets_.push_back(eqn);
      }
    }
    rhs_mom(fnc - 1) = mome;

    if (fnc == 1)
    {
      if (rhs_mom.InfNorm() > 1e-5 && rhs_mom(0) < 0.0)
        dserror("negaive volume in base function integration. is ordering of vertices right?");
    }
  }

  /*double intee = rhs_mom(56)+rhs_mom(68)+rhs_mom(72)+rhs_mom(83);
  std::cout<<std::setprecision(20)<<"the INTEGRAL = "<<intee<<"\n";*/

  // set the volume of this volumecell
  // the volume from local coordinates is converted in terms of global coordinates
  double volGlobal = 0.0;
  switch (elem1_->Shape())
  {
    case DRT::Element::hex8:
    {
      volGlobal =
          elem1_->ScalarFromLocalToGlobal<3, DRT::Element::hex8>(rhs_mom(0), "LocalToGlobal");
      break;
    }
    case DRT::Element::tet4:
    {
      volGlobal =
          elem1_->ScalarFromLocalToGlobal<3, DRT::Element::tet4>(rhs_mom(0), "LocalToGlobal");
      break;
    }
    case DRT::Element::wedge6:
    {
      volGlobal =
          elem1_->ScalarFromLocalToGlobal<3, DRT::Element::wedge6>(rhs_mom(0), "LocalToGlobal");
      break;
    }
    case DRT::Element::pyramid5:
    {
      volGlobal =
          elem1_->ScalarFromLocalToGlobal<3, DRT::Element::pyramid5>(rhs_mom(0), "LocalToGlobal");
      break;
    }
    default:
      dserror("unsupported integration cell type ( cell type = %s )",
          DRT::DistypeToString(elem1_->Shape()).c_str());
      exit(EXIT_FAILURE);
  }
  volcell_->SetVolume(volGlobal);

  return rhs_mom;
}

/*----------------------------------------------------------------------------*
    compute the gaussian points of the volumecell with "numeach" points in each
    3-directions
    numeach should be more than 1
    uses ray tracing method
*-----------------------------------------------------------------------------*/
bool GEO::CUT::VolumeIntegration::compute_Gaussian_points(int numeach)
{
  Teuchos::RCP<BoundingBox> box1 = Teuchos::rcp(BoundingBox::Create(*volcell_, elem1_));
  double minn[3], maxx[3];
  minn[0] = box1->minx();
  maxx[0] = box1->maxx();
  minn[1] = box1->miny();
  maxx[1] = box1->maxy();
  minn[2] = box1->minz();
  maxx[2] = box1->maxz();

  std::vector<std::vector<double>> zcoord;
  std::vector<std::vector<double>> ycoord;
  get_zcoordinates(zcoord, ycoord);

  // min z plane which contains significant area
  double zmin = minn[2] + 0.01 * (maxx[2] - minn[2]), zmax = maxx[2] - 0.01 * (maxx[2] - minn[2]);
  while (1)
  {
    bool area = false;
    std::vector<std::vector<double>> InPlane;
    area = IsContainArea(minn, maxx, zmin, InPlane, zcoord, ycoord, 0.01, numeach);
    if (area)
    {
      gaus_pts_.insert(gaus_pts_.end(), InPlane.begin(), InPlane.end());
      break;
    }
    zmin += 0.01 * (maxx[2] - minn[2]);
    if ((zmax - zmin) < 0.01 * (maxx[2] - minn[2])) break;
  }

  // max z plane which contains significant area
  while (1)
  {
    bool area = false;
    std::vector<std::vector<double>> InPlane;
    area = IsContainArea(minn, maxx, zmax, InPlane, zcoord, ycoord, 0.01, numeach);
    if (area)
    {
      gaus_pts_.insert(gaus_pts_.end(), InPlane.begin(), InPlane.end());
      break;
    }
    zmax -= 0.01 * (maxx[2] - minn[2]);
    if ((zmax - zmin) < 0.01 * (maxx[2] - minn[2])) break;
  }
  bool wei = true;
  if ((zmax - zmin) < 0.01 * (maxx[2] - minn[2]))
  {
    // check for very thin volumecells
    // works fine but time taken is slightly higher
    double zmin = minn[2] + 0.005 * (maxx[2] - minn[2]),
           zmax = maxx[2] - 0.005 * (maxx[2] - minn[2]);
    while (1)
    {
      bool area = false;
      std::vector<std::vector<double>> InPlane;
      area = IsContainArea(minn, maxx, zmin, InPlane, zcoord, ycoord, 0.001, numeach);
      if (area)
      {
        gaus_pts_.insert(gaus_pts_.end(), InPlane.begin(), InPlane.end());
        break;
      }
      zmin += 0.001 * (maxx[2] - minn[2]);
      if ((zmax - zmin) < 0.001 * (maxx[2] - minn[2])) break;
    }
    while (1)
    {
      bool area = false;
      std::vector<std::vector<double>> InPlane;
      area = IsContainArea(minn, maxx, zmax, InPlane, zcoord, ycoord, 0.001, numeach);
      if (area)
      {
        gaus_pts_.insert(gaus_pts_.end(), InPlane.begin(), InPlane.end());
        break;
      }
      zmax -= 0.001 * (maxx[2] - minn[2]);
      if ((zmax - zmin) < 0.001 * (maxx[2] - minn[2])) break;
    }

    if (gaus_pts_.size() == 0)
      wei = false;
    else
      wei = true;
#ifdef DEBUGCUTLIBRARY
    std::cout << "number of Gauss points " << gaus_pts_.size() << std::endl;
#endif
    return wei;
  }

  // find z-planes in between zmin and zmax to generate Gauss points
  else
  {
    int num = numeach;
    double dzz = (zmax - zmin) / (num - 1);
    std::vector<double> zplane;
    for (int i = 0; i < num; i++) zplane.push_back(zmin + i * dzz);
    bool previous = true;
    for (int i = 1; i < num - 1; i++)
    {
      bool area = false;
      std::vector<std::vector<double>> InPlane;
      area = IsContainArea(minn, maxx, zplane[i], InPlane, zcoord, ycoord, 0.01, numeach);
      if (area)
      {
        gaus_pts_.insert(gaus_pts_.end(), InPlane.begin(), InPlane.end());
        InPlane.clear();
        previous = true;
        continue;
      }

      // if the considered z-plane does not contain significant area, the interval is subdivided to
      // check whether any plane lies in between
      else
      {
        if (previous == false) continue;
        double dzzm = dzz;
        for (int k = 0; k < 4; k++)
        {
          dzzm += 0.5 * dzzm;
          bool innerarea = false;
          double zz = zplane[i] - dzzm;
          innerarea = IsContainArea(minn, maxx, zz, InPlane, zcoord, ycoord, 0.01, numeach);
          if (innerarea)
          {
            gaus_pts_.insert(gaus_pts_.end(), InPlane.begin(), InPlane.end());
            InPlane.clear();
            break;
          }
        }
        previous = false;
      }
    }
  }

#ifdef DEBUGCUTLIBRARY
  std::cout << "number of Gauss points" << gaus_pts_.size() << std::endl;
#endif
  return wei;
}

/*-------------------------------------------------------------------------------------------------------------------*
    Store the z- and y-coordinates of the all corner points which will be used to find whether the
intersection point lies inside the volume or not
*--------------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeIntegration::get_zcoordinates(
    std::vector<std::vector<double>> &zcoord, std::vector<std::vector<double>> &ycoord)
{
  const plain_facet_set &facete = volcell_->Facets();
  int faceno = 0;
  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    std::vector<double> thisplane1, thisplane2;

    // since the imaginary line never intersect this plane
    if (fabs(eqn_facets_[faceno][0]) < 0.00000001)
    {
      thisplane1.push_back(0.0);
      zcoord.push_back(thisplane1);
      ycoord.push_back(thisplane1);
      faceno++;
      continue;
    }
    Facet *face1 = *i;
    std::vector<std::vector<double>> corLocal;
    face1->CornerPointsLocal(elem1_, corLocal);
    for (std::vector<std::vector<double>>::iterator k = corLocal.begin(); k != corLocal.end(); k++)
    {
      std::vector<double> &coords1 = *k;

      thisplane1.push_back(coords1[1]);
      thisplane2.push_back(coords1[2]);
    }
    ycoord.push_back(thisplane1);
    zcoord.push_back(thisplane2);
    faceno++;
  }
}

/*-----------------------------------------------------------------------------------------*
              check whether the generated ray intersect any of the facets
              if so generate gauss points along the ray
*------------------------------------------------------------------------------------------*/
bool GEO::CUT::VolumeIntegration::IsIntersect(double *pt, double *mini, double *maxi,
    std::vector<std::vector<double>> &linePts, std::vector<std::vector<double>> zcoord,
    std::vector<std::vector<double>> ycoord, double toler, int numeach)
{
  bool intersect = false;
  std::vector<int> planeinter, InterFaces;

  // stores all the facets which is not in x-y or x-z
  // only in the remaining planes a horizontal line possibly intersect
  for (unsigned i = 0; i < eqn_facets_.size(); i++)
  {
    if (fabs(eqn_facets_[i][0]) > 0.0000001) planeinter.push_back(i);
  }

  // stores the facets which are actually cut by the line
  for (unsigned i = 0; i < planeinter.size(); i++)
  {
    int faceno = planeinter[i];
    std::vector<double> planez = zcoord[faceno];
    std::vector<double> planey = ycoord[faceno];
    int cutno = 0;
    // check whether the intersection point lies inside the facet area
    cutno = pnpoly(planez.size(), planey, planez, pt[1], pt[2]);
    if (cutno == 1) InterFaces.push_back(faceno);
  }

  if (InterFaces.size() < 2)
  {
    intersect = false;
    return intersect;
  }

  // simple geometries result in two intersections
  if (InterFaces.size() == 2)
  {
    int fa1 = InterFaces[0];
    int fa2 = InterFaces[1];

    double int1 =
        (eqn_facets_[fa1][3] - eqn_facets_[fa1][1] * pt[1] - eqn_facets_[fa1][2] * pt[2]) /
        eqn_facets_[fa1][0];
    double int2 =
        (eqn_facets_[fa2][3] - eqn_facets_[fa2][1] * pt[1] - eqn_facets_[fa2][2] * pt[2]) /
        eqn_facets_[fa2][0];

    if (fabs(int2 - int1) < toler * (maxi[0] - mini[0]))
    {
      intersect = false;
      return intersect;
    }
    else
      intersect = true;

    std::vector<double> inter1, inter2;
    inter1.push_back(int1);
    inter1.push_back(pt[1]);
    inter1.push_back(pt[2]);
    inter2.push_back(int2);
    inter2.push_back(pt[1]);
    inter2.push_back(pt[2]);

#if 1  // old way of point distribution method in intersection lines
    if (fabs(inter2[0] - inter1[0]) < 0.025 * (maxi[0] - mini[0]))
    {
      std::vector<double> middle(3);
      middle[0] = 0.5 * (inter2[0] + inter1[0]);
      middle[1] = inter2[1];
      middle[2] = inter2[2];
      linePts.push_back(middle);
      intersect = true;
    }
    else if (fabs(inter2[0] - inter1[0]) < 0.1 * (maxi[0] - mini[0]))
    {
      OnLine(inter1, inter2, linePts, 2);
      intersect = true;
    }
    else if (fabs(inter2[0] - inter1[0]) < 0.2 * (maxi[0] - mini[0]))
    {
      OnLine(inter1, inter2, linePts, 3);
      intersect = true;
    }
    else if (fabs(inter2[0] - inter1[0]) < 0.3 * (maxi[0] - mini[0]))
    {
      OnLine(inter1, inter2, linePts, 4);
      intersect = true;
    }
    else
      OnLine(inter1, inter2, linePts, numeach);
#endif
#if 0  // scaled point distribution method in intersection lines (produced more error than old
       // method)
    int numX = (int)((fabs(inter2[0]-inter1[0])/(maxi[0]-mini[0]))*numeach+1);
    if(numX==1)
    {
      std::vector<double> middle(3);
      middle[0] = 0.5*(inter2[0]+inter1[0]);
      middle[1] = inter2[1];
      middle[2] = inter2[2];
      linePts.push_back(middle);
      intersect = true;
    }
    else
    {
      OnLine(inter1,inter2,linePts,numX);
      intersect = true;
    }
#endif
    return intersect;
  }
  else
  {
// map is useful since we need to arrange the elements from minimum x-cut value
#ifdef DEBUGCUTLIBRARY
    std::cout << "in more than two cuts" << std::endl;
#endif
    std::map<std::vector<double>, int> interPoints;
    for (std::vector<int>::iterator i = InterFaces.begin(); i != InterFaces.end(); i++)
    {
      int fa1 = *i;
      double int1 =
          (eqn_facets_[fa1][3] - eqn_facets_[fa1][1] * pt[1] - eqn_facets_[fa1][2] * pt[2]) /
          eqn_facets_[fa1][0];
      std::vector<double> inter1;
      inter1.push_back(int1);
      inter1.push_back(pt[1]);
      inter1.push_back(pt[2]);
      interPoints[inter1] = fa1;
    }

    const plain_facet_set &facete = volcell_->Facets();
    plain_facet_set::const_iterator f = facete.begin();
    unsigned cut_count = 0;
    for (std::map<std::vector<double>, int>::iterator i = interPoints.begin();
         i != interPoints.end(); i++)
    {
      bool ptsInside = false;
      unsigned face1 = i->second;
      plain_facet_set::const_iterator f1 = f;
      for (int k = 0; k < i->second; k++) f1++;
      Facet *facet1 = *f1;
      std::vector<double> inter1 = i->first;
      i++;
      unsigned face2 = i->second;
      plain_facet_set::const_iterator f2 = f;
      for (int k = 0; k < i->second; k++) f2++;
      Facet *facet2 = *f2;
      std::vector<double> inter2 = i->first;
      i--;

      // Among two consequetive facets which are cut by the line, one must be the cut surface
      if (facet1->OnCutSide())
      {
        if (eqn_facets_[face1][0] < 0.0) ptsInside = true;
      }
      else if (facet2->OnCutSide())
      {
        if (eqn_facets_[face2][0] > 0.0) ptsInside = true;
      }
      else
      {
#ifdef DEBUGCUTLIBRARY
        std::cout << "The assumption that one must be a cut surface is false" << std::endl;
#endif
      }

      if (ptsInside)
      {
        if (interPoints.size() == 2)  // if intersection is on the triangulated line, results only
                                      // in two face intersections
        {
          OnLine(inter1, inter2, linePts, numeach);
          intersect = true;
          break;
        }
        if ((inter2[0] - inter1[0]) < toler * (maxi[0] - mini[0]))
        {
          std::vector<double> middle(3);
          middle[0] = 0.5 * (inter2[0] + inter1[0]);
          middle[1] = inter2[1];
          middle[2] = inter2[2];
          linePts.push_back(middle);
          intersect = true;
        }
        else if ((inter2[0] - inter1[0]) < 0.05 * (maxi[0] - mini[0]))
        {
          OnLine(inter1, inter2, linePts, 2 /*numeach/3+1*/);
          intersect = true;
        }
        else
        {
          OnLine(inter1, inter2, linePts, numeach / 2 + 1);
          intersect = true;
        }
      }
      cut_count++;
      if (cut_count == (interPoints.size() - 1)) break;
    }
  }
  return intersect;
}

/*-------------------------------------------------------------------------------------------------------------------*
         Check whether the intersection point, which is in the plane containing the facet, actually
         lies with in the facet area
*--------------------------------------------------------------------------------------------------------------------*/
int GEO::CUT::VolumeIntegration::pnpoly(
    int npol, std::vector<double> xp, std::vector<double> yp, double x, double y)
{
  // check whether given point is one of the corner points
  for (int i = 0; i < npol; i++)
  {
    if (fabs(xp[i] - x) < 1e-8 && fabs(yp[i] - y) < 1e-8) return 1;
  }

  // check whether given point is on the edge
  for (int i = 0; i < npol; i++)
  {
    double x1 = xp[i], y1 = yp[i], x2 = xp[(i + 1) % npol], y2 = yp[(i + 1) % npol];
    // line with constant ind1
    if (fabs(x2 - x1) < 1e-8 && fabs(x2 - x) < 1e-8)
    {
      if (y > std::min(y2, y1) && y < std::max(y2, y1)) return 1;
    }
    // line with constant ind2
    if (fabs(y2 - y1) < 1e-8 && fabs(y2 - y) < 1e-8)
    {
      if (x > std::min(x2, x1) && x < std::max(x2, x1)) return 1;
    }

    // pt fall on edge if it satisfies eqn of edge
    if (fabs(y2 - y1) < 1e-8 && fabs(y2 - y) < 1e-8)
    {
      double xpre = x1 + (x2 - x1) / (y2 - y1) * (y - y1);
      if (fabs(xpre - x) < 1e-8) return 1;
    }
  }

  // check inside
  int i, j, c = 0;
  for (i = 0, j = npol - 1; i < npol; j = i++)
  {
    if ((((yp[i] <= y) && (y < yp[j])) || ((yp[j] <= y) && (y < yp[i]))) &&
        (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
      c = !c;
  }
  return c;
}

/*-------------------------------------------------------------------------------------------------------------------*
???
*--------------------------------------------------------------------------------------------------------------------*/
int GEO::CUT::VolumeIntegration::pnpoly(const std::vector<std::vector<double>> &xp,
    const LINALG::Matrix<3, 1> &pt, GEO::CUT::ProjectionDirection projType)
{
  int npol = xp.size();
  int ind1 = 1, ind2 = 2;
  if (projType == GEO::CUT::proj_y)
  {
    ind1 = 2;
    ind2 = 0;
  }
  else if (projType == GEO::CUT::proj_z)
  {
    ind1 = 0;
    ind2 = 1;
  }

  // check whether given point is one of the corner points
  for (int i = 0; i < npol; i++)
  {
    if (fabs(xp[i][ind1] - pt(ind1, 0)) < 1e-8 && fabs(xp[i][ind2] - pt(ind2, 0)) < 1e-8) return 1;
  }

  // check whether given point is on the edge
  for (int i = 0; i < npol; i++)
  {
    std::vector<double> end1 = xp[i];
    std::vector<double> end2 = xp[(i + 1) % npol];

    // line with constant ind1
    if (fabs(end1[ind1] - end2[ind1]) < 1e-8 &&
        fabs(end1[ind1] - pt(ind1, 0)) < 1e-8)  // pt is on infinite line check
    {
      // pt is within limits check
      if (pt(ind2, 0) > std::min(end1[ind2], end2[ind2]) &&
          pt(ind2, 0) < std::max(end1[ind2], end2[ind2]))
        return 1;
      /*else
        return 0;*/
    }
    // line with constant ind2
    if (fabs(end1[ind2] - end2[ind2]) < 1e-8 && fabs(end1[ind2] - pt(ind2, 0)) < 1e-8)
    {
      if (pt(ind1, 0) > std::min(end1[ind1], end2[ind1]) &&
          pt(ind1, 0) < std::max(end1[ind1], end2[ind1]))
        return 1;
      /*else
        return 0;*/
    }

    // pt fall on edge if it satisfies eqn of edge
    if (pt(ind2, 0) > std::min(end1[ind2], end2[ind2]) &&  // check within limits
        pt(ind2, 0) < std::max(end1[ind2], end2[ind2]))
    {
      double xpre = end1[ind1] + (end2[ind1] - end1[ind1]) / (end2[ind2] - end1[ind2]) *
                                     (pt(ind2, 0) - end1[ind2]);
      if (fabs(xpre - pt(ind1, 0)) < 1e-8)  // check eqn of line is satisfied with given pt
        return 1;
    }
  }

  // check for inside
  int i, j, c = 0;
  for (i = 0, j = npol - 1; i < npol; j = i++)
  {
    if ((((xp[i][ind2] <= pt(ind2, 0)) && (pt(ind2, 0) < xp[j][ind2])) ||
            ((xp[j][ind2] <= pt(ind2, 0)) && (pt(ind2, 0) < xp[i][ind2]))) &&
        (pt(ind1, 0) < (xp[j][ind1] - xp[i][ind1]) * (pt(ind2, 0) - xp[i][ind2]) /
                               (xp[j][ind2] - xp[i][ind2]) +
                           xp[i][ind1]))
      c = !c;
  }
  return c;
}
/*-------------------------------------------------------------------------------------------------------------------*
        Check whether the particular z-plane of the volumecell contains significant area so as to
 distribute the Gauss points in that plane
 *--------------------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::VolumeIntegration::IsContainArea(double minn[3], double maxx[3], double &zmin,
    std::vector<std::vector<double>> &pts, std::vector<std::vector<double>> zcoord,
    std::vector<std::vector<double>> ycoord, double toler, int numeach)
{
  bool isArea = true;
  double dx1[3];
  dx1[0] = maxx[0] - minn[0], dx1[1] = maxx[1] - minn[1], dx1[2] = maxx[2] - minn[2];
  double xmin = minn[0] + toler * dx1[0], ymin = minn[1], ymax = maxx[1];

  // check for the lowest line
  while (1)
  {
    std::vector<std::vector<double>> linePts;
    ymin += toler * dx1[1];
    bool intersec = false;
    double a[] = {xmin, ymin, zmin};
    intersec = IsIntersect(a, minn, maxx, linePts, zcoord, ycoord, toler, numeach);

    // the area of the solid volume in this particular z-plane is almost negligible
    // if the area of volumecell is less than 1% of bounding box
    if ((ymax - ymin) < toler * dx1[1])
    {
      isArea = false;
      return isArea;
    }
    if (intersec)
    {
      pts.insert(pts.end(), linePts.begin(), linePts.end());
      break;
    }
  }

  // check for the topmost line
  while (1)
  {
    std::vector<std::vector<double>> linePts;
    ymax -= toler * dx1[1];

    bool intersec = false;
    double a[] = {xmin, ymax, zmin};
    intersec = IsIntersect(a, minn, maxx, linePts, zcoord, ycoord, toler, numeach);
    if ((ymax - ymin) < toler * dx1[1])
    {
      isArea = false;
      return isArea;
    }
    if (intersec)
    {
      pts.insert(pts.end(), linePts.begin(), linePts.end());
      break;
    }
  }

  // generate points in between the topmost and lowest line
  int num = numeach;
  double dyy = (ymax - ymin) / (num - 1);
  std::vector<double> yplane;
  for (int i = 0; i < num; i++) yplane.push_back(ymin + i * dyy);
  bool previous = true;
  for (int i = 1; i < num - 1; i++)
  {
    double a[] = {xmin, yplane[i], zmin};
    bool intersec = false;
    std::vector<std::vector<double>> linePts;
    intersec = IsIntersect(a, minn, maxx, linePts, zcoord, ycoord, toler, numeach);
    if (intersec)
    {
      pts.insert(pts.end(), linePts.begin(), linePts.end());
      linePts.clear();
      previous = true;
      continue;
    }
    else
    {
      if (previous == false) continue;
      double dyym = dyy;
      for (int k = 0; k < 4; k++)
      {
        dyym += 0.5 * dyym;
        double a[] = {xmin, yplane[i] - dyym, zmin};
        bool innerintersec = false;
        std::vector<std::vector<double>> linePts;
        innerintersec = IsIntersect(a, minn, maxx, linePts, zcoord, ycoord, toler, numeach);

        if (innerintersec)
        {
          pts.insert(pts.end(), linePts.begin(), linePts.end());
          linePts.clear();
          break;
        }
      }
      previous = false;
    }
  }
  return isArea;
}

/*-------------------------------------------------------------------------------------------------------------------*
        Generates equally spaced "num" number of points on the line whose end points are specified
by inter1 and inter2
*-------------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeIntegration::OnLine(std::vector<double> inter1, std::vector<double> inter2,
    std::vector<std::vector<double>> &linePts, int num)
{
  std::vector<double> left, right;
  if (inter1[0] < inter2[0])
  {
    left = inter1;
    right = inter2;
  }
  else
  {
    left = inter2;
    right = inter1;
  }
  double xlen = right[0] - left[0];
  inter1[0] = left[0] + 0.05 * xlen;
  inter2[0] = right[0] - 0.05 * xlen;
  double xdiff = (inter2[0] - inter1[0]) / (num - 1);
  for (int i = 0; i < num; i++)
  {
    std::vector<double> temp(3);
    temp[0] = inter1[0] + i * xdiff;
    temp[1] = inter1[1];
    temp[2] = inter1[2];
    linePts.push_back(temp);
  }
}

/*-------------------------------------------------------------------------------------------------------------------*
                                                form the moment fitting matrix
*--------------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeIntegration::moment_fitting_matrix(
    std::vector<std::vector<double>> &mom, std::vector<std::vector<double>> gauspts)
{
  for (int i = 0; i < num_func_; i++)
  {
    int k = 0;
    for (std::vector<std::vector<double>>::iterator j = gauspts.begin(); j != gauspts.end(); j++)
    {
      std::vector<double> cordi = *j;
      mom[i][k] = GEO::CUT::base_function(cordi, i + 1);
      k++;
    }
  }
}

/*-------------------------------------------------------------------------------------------------------------------*
    Compute Gauss point weights by solving the moment fitting equations and returns the coordinates
of Gauss points and their corresponding weights
*--------------------------------------------------------------------------------------------------------------------*/
Epetra_SerialDenseVector GEO::CUT::VolumeIntegration::compute_weights()
{
  Epetra_SerialDenseVector rhs_moment(num_func_);
  rhs_moment = compute_rhs_moment();

  bool wei;
  // we should ask for more than 1 point in each direction
  int numeach = 7;
  Epetra_SerialDenseVector weights;
  while (1)
  {
    gaus_pts_.clear();
    wei = compute_Gaussian_points(numeach);

    if (wei)
    {
      std::vector<std::vector<double>> moment_matrix(
          num_func_, std::vector<double>(gaus_pts_.size()));
      moment_fitting_matrix(moment_matrix, gaus_pts_);

#if 0  // addition of linear combination of monomials. found no improvement even if it is added
        if(num_func_>1)
          FirstOrderAdditionalTerms(moment_matrix,rhs_moment);
        if(num_func_>4)
          SecondOrderAdditionalTerms(moment_matrix,rhs_moment);
        if(num_func_>10)
          ThirdOrderAdditionalTerms(moment_matrix,rhs_moment);
        if(num_func_>20)
          FourthOrderAdditionalTerms(moment_matrix,rhs_moment);
        if(num_func_>35)
          FifthOrderAdditionalTerms(moment_matrix,rhs_moment);
        if(num_func_>56)
          SixthOrderAdditionalTerms(moment_matrix,rhs_moment);
#endif

      // if all the elements in a row of the moment fitting matrix are zero, then the row has to be
      // deleted this ensures non-zero diagonal elements in the matrix
      // REquires further checking
      /*vector<int> deleteRowNos;
      for(unsigned row=0;row<moment_matrix.size();row++)
      {
        bool deleteRow=true;
        for(unsigned col=0;col!=moment_matrix[0].size();col++)
        {
          if(fabs(moment_matrix[row][col])>1e-8)
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
  //        std::cout<<delno<<"\n";
          moment_matrix.erase(moment_matrix.begin()+delno);
        }
      }*/

      LeastSquares least(moment_matrix, rhs_moment);
      weights.Size(moment_matrix[0].size());
      weights = least.linear_least_square();
    }
    else  // the considered volumecell has negligible volume and can be eliminated
    {
      gaus_pts_.clear();
      std::vector<double> zer(3);
      zer[0] = 0.0;
      zer[1] = 0.0;
      zer[2] = 0.0;
      gaus_pts_.push_back(zer);
      weights.Size(1);
      weights(0) = 0.0;
      break;
    }

    Epetra_SerialDenseVector err(num_func_);
    for (int i = 0; i < num_func_; i++)
    {
      err(i) = 0.0;
      for (unsigned j = 0; j < gaus_pts_.size(); j++)
      {
        err(i) += weights(j) * base_function(gaus_pts_[j], i + 1);
      }
      if (fabs(rhs_moment(i)) > 1e-8)
        err(i) = (err(i) - rhs_moment(i)) / rhs_moment(i);
      else
        err(i) = err(i) - rhs_moment(i);
    }

#if 0  // call the computation of error when integrating specific functions
    ErrorForSpecificFunction(rhs_moment,weights,numeach);
#endif

    const double maxError = err.InfNorm();
#ifdef DEBUGCUTLIBRARY
    std::cout << "max error = " << maxError << "\n";
#endif
    if (maxError < 1e-10 || numeach == 13)
      break;
    else
      numeach++;
    break;
  }

#ifdef DEBUGCUTLIBRARY
  std::cout << "volume = " << rhs_moment(0) << "\t" << rhs_moment(1) << "\t" << rhs_moment(4)
            << std::endl;
#endif

#ifdef DEBUGCUTLIBRARY
  GaussPointGmsh();
#endif

  return weights;
}

// Writes the Geometry of volumecell and location of Gauss points in GMSH format output file
void GEO::CUT::VolumeIntegration::GaussPointGmsh()
{
  volcell_->DumpGmshGaussPointsMomFit(gaus_pts_);
}

// compute integration of x+y,y+z and x+z values from the integration of x, y and z values
void GEO::CUT::VolumeIntegration::FirstOrderAdditionalTerms(
    std::vector<std::vector<double>> &mat, Epetra_SerialDenseVector &rhs)
{
  unsigned int i = mat.size(), kk = mat[0].size();
  // no of additional elements is n(n+1)/2 where n=(no_of_monomials-1). Here no_of
  // monomials=3=>(x,y,z)
  rhs.Resize(i + 3);
  mat.resize(i + 3, std::vector<double>(kk));
  unsigned int ibegin = 1, iend = 3;  // the row numbers that contain the first order terms in rhs
  for (unsigned j = ibegin; j <= iend; j++)
  {
    for (unsigned k = j + 1; k <= iend; k++)
    {
      rhs(i) = rhs(j) + rhs(k);
      for (unsigned m = 0; m < mat[0].size(); m++) mat[i][m] = mat[j][m] + mat[k][m];
      i++;
    }
  }
}

// integration of linear combination of second order terms like x^2+xy+y^2+yz
void GEO::CUT::VolumeIntegration::SecondOrderAdditionalTerms(
    std::vector<std::vector<double>> &mat, Epetra_SerialDenseVector &rhs)
{
  unsigned int i = mat.size(), kk = mat[0].size();
  rhs.Resize(i + 15);
  mat.resize(i + 15, std::vector<double>(kk));
  unsigned int ibegin = 4, iend = 9;  // the row numbers that contain the second order terms in rhs
  for (unsigned j = ibegin; j <= iend; j++)
  {
    for (unsigned k = j + 1; k <= iend; k++)
    {
      rhs(i) = rhs(j) + rhs(k);
      for (unsigned m = 0; m < mat[0].size(); m++) mat[i][m] = mat[j][m] + mat[k][m];
      i++;
    }
  }
}

// integration of linear combination of third order terms x^3+xyz
void GEO::CUT::VolumeIntegration::ThirdOrderAdditionalTerms(
    std::vector<std::vector<double>> &mat, Epetra_SerialDenseVector &rhs)
{
  unsigned int i = mat.size(), kk = mat[0].size();
  rhs.Resize(i + 45);
  mat.resize(i + 45, std::vector<double>(kk));
  unsigned int ibegin = 10, iend = 19;  // the row numbers that contain the third order terms in rhs
  for (unsigned j = ibegin; j <= iend; j++)
  {
    for (unsigned k = j + 1; k <= iend; k++)
    {
      rhs(i) = rhs(j) + rhs(k);
      for (unsigned m = 0; m < mat[0].size(); m++) mat[i][m] = mat[j][m] + mat[k][m];
      i++;
    }
  }
}

void GEO::CUT::VolumeIntegration::FourthOrderAdditionalTerms(
    std::vector<std::vector<double>> &mat, Epetra_SerialDenseVector &rhs)
{
  unsigned int i = mat.size(), kk = mat[0].size();
  rhs.Resize(i + 105);
  mat.resize(i + 105, std::vector<double>(kk));
  unsigned int ibegin = 20, iend = 34;  // the row numbers that contain the fourth order terms in
                                        // rhs
  for (unsigned j = ibegin; j <= iend; j++)
  {
    for (unsigned k = j + 1; k <= iend; k++)
    {
      rhs(i) = rhs(j) + rhs(k);
      for (unsigned m = 0; m < mat[0].size(); m++) mat[i][m] = mat[j][m] + mat[k][m];
      i++;
    }
  }
}

void GEO::CUT::VolumeIntegration::FifthOrderAdditionalTerms(
    std::vector<std::vector<double>> &mat, Epetra_SerialDenseVector &rhs)
{
  unsigned int i = mat.size(), kk = mat[0].size();
  rhs.Resize(i + 210);
  mat.resize(i + 210, std::vector<double>(kk));
  unsigned int ibegin = 35, iend = 55;  // the row numbers that contain the fifth order terms in rhs
  for (unsigned j = ibegin; j <= iend; j++)
  {
    for (unsigned k = j + 1; k <= iend; k++)
    {
      rhs(i) = rhs(j) + rhs(k);
      for (unsigned m = 0; m < mat[0].size(); m++) mat[i][m] = mat[j][m] + mat[k][m];
      i++;
    }
  }
}

void GEO::CUT::VolumeIntegration::SixthOrderAdditionalTerms(
    std::vector<std::vector<double>> &mat, Epetra_SerialDenseVector &rhs)
{
  unsigned int i = mat.size(), kk = mat[0].size();
  rhs.Resize(i + 561);
  mat.resize(i + 561, std::vector<double>(kk));
  unsigned int ibegin = 56, iend = 83;  // the row numbers that contain the fifth order terms in rhs
  for (unsigned j = ibegin; j <= iend; j++)
  {
    for (unsigned k = j + 1; k <= iend; k++)
    {
      rhs(i) = rhs(j) + rhs(k);
      for (unsigned m = 0; m < mat[0].size(); m++) mat[i][m] = mat[j][m] + mat[k][m];
      i++;
    }
  }
}

/*  Computes the error introduced by the generated integration rule for integrating some specific
   functions Used only in post-processing    */
void GEO::CUT::VolumeIntegration::ErrorForSpecificFunction(
    Epetra_SerialDenseVector rhs_moment, Epetra_SerialDenseVector weights, int numeach)
{
  static std::vector<int> gausSize;
  gausSize.push_back(gaus_pts_.size());
  static std::vector<std::vector<double>> errAccu;
  std::vector<double> error(7);

  double chek = 0.0, val = 0.0;
  /*for(unsigned i=0;i<gaus_pts_.size();i++)
  {
    double xx = gaus_pts_[i][0];
    double yy = gaus_pts_[i][1];
    double zz = gaus_pts_[i][2];
    chek += (pow(xx,6)+xx*pow(yy,4)*zz+xx*xx*yy*yy*zz*zz+pow(zz,6))*weights(i);
  }
  std::cout<<"MOMENT FITTING :: Integral value = "<<chek<<"\n";*/

  chek = 0.0;
  for (unsigned i = 0; i < gaus_pts_.size(); i++)
  {
    chek += 5 * weights(i);
  }
  val = 5 * rhs_moment(0);
  error[0] = (val - chek) / val;

  chek = 0.0;
  for (unsigned i = 0; i < gaus_pts_.size(); i++)
  {
    chek += 2 * weights(i) * gaus_pts_[i][0] + 3 * weights(i) * gaus_pts_[i][1] +
            5 * weights(i) * gaus_pts_[i][2];
  }
  val = 2 * rhs_moment(1) + 3 * rhs_moment(2) + 5 * rhs_moment(3);
  error[1] = (val - chek) / val;

  chek = 0.0;
  for (unsigned i = 0; i < gaus_pts_.size(); i++)
  {
    chek += weights(i) * (gaus_pts_[i][0] * gaus_pts_[i][0] - gaus_pts_[i][0] * gaus_pts_[i][1] +
                             gaus_pts_[i][1] * gaus_pts_[i][2] + gaus_pts_[i][2] * gaus_pts_[i][2]);
  }
  val = rhs_moment(4) - rhs_moment(5) + rhs_moment(8) + rhs_moment(9);
  error[2] = (val - chek) / val;

  chek = 0.0;
  for (unsigned i = 0; i < gaus_pts_.size(); i++)
  {
    chek += weights(i) *
            (pow(gaus_pts_[i][0], 3) + gaus_pts_[i][0] * gaus_pts_[i][1] * gaus_pts_[i][2] +
                pow(gaus_pts_[i][1], 3));
  }
  val = rhs_moment(10) + rhs_moment(14) + rhs_moment(16);
  error[3] = (val - chek) / val;

  chek = 0.0;
  for (unsigned i = 0; i < gaus_pts_.size(); i++)
  {
    chek += weights(i) *
            (pow(gaus_pts_[i][0], 4) - 7 * gaus_pts_[i][0] * gaus_pts_[i][1] * gaus_pts_[i][1] +
                pow(gaus_pts_[i][2], 3));
  }
  val = rhs_moment(20) - 7 * rhs_moment(13) + rhs_moment(19);
  error[4] = (val - chek) / val;

  chek = 0.0;
  for (unsigned i = 0; i < gaus_pts_.size(); i++)
  {
    chek += weights(i) * (pow(gaus_pts_[i][1], 5) + gaus_pts_[i][0] * pow(gaus_pts_[i][2], 4) +
                             pow(gaus_pts_[i][2], 5));
  }
  val = rhs_moment(50) + rhs_moment(49) + rhs_moment(55);
  error[5] = (val - chek) / val;

  chek = 0.0;
  for (unsigned i = 0; i < gaus_pts_.size(); i++)
  {
    chek += weights(i) *
            (pow(gaus_pts_[i][0], 6) + pow(gaus_pts_[i][1], 6) + pow(gaus_pts_[i][2], 6) +
                pow(gaus_pts_[i][0], 2) * pow(gaus_pts_[i][1], 2) * pow(gaus_pts_[i][2], 2) +
                gaus_pts_[i][0] * pow(gaus_pts_[i][1], 5));
  }
  val = rhs_moment(56) + rhs_moment(77) + rhs_moment(83) + rhs_moment(68) + rhs_moment(71);
  error[6] = (val - chek) / val;

  errAccu.push_back(error);

  if (numeach == 18)
  {
    std::string filename = "wrong";
    std::ofstream file;

    std::stringstream out;
    out << "funcOutput.dat";
    filename = out.str();
    file.open(filename.c_str());

    for (unsigned i = 0; i < errAccu.size(); i++)
    {
      file << gausSize[i] << "\t";
      for (unsigned j = 0; j < errAccu[0].size(); j++) file << errAccu[i][j] << "\t";
      file << "\n";
    }
    file.close();
  }
}

/*-------------------------------------------------------------------------------------*
 * Check whether the point with this element Local coordinates is inside,              *
 * outside or on boundary of this volumecell                            sudhakar 07/12 *
 *-------------------------------------------------------------------------------------*/
std::string GEO::CUT::VolumeIntegration::IsPointInside(LINALG::Matrix<3, 1> &rst)
{
  //  std::cout << RED << "** GEO::CUT::VolumeIntegration::IsPointInside -- This function does not
  //  work for general volumecells (concave, facets with holes, ...)! **" << END_COLOR << std::endl;
  const plain_facet_set &facete = volcell_->Facets();

  //--------------------------------------------------------------------------------
  // Step 1: Classify facets into facets having zero normal in x-direction and other
  //--------------------------------------------------------------------------------
  std::vector<plain_facet_set::const_iterator> XFacets;     // facets with non-zero nx
  std::vector<plain_facet_set::const_iterator> NotXFacets;  // facets with zero nx
  std::vector<std::vector<double>> Eqnplane;                // eqn of plane for all facets

  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet *fe = *i;
    std::vector<std::vector<double>> cornersLocal;
    fe->CornerPointsLocal(elem1_, cornersLocal);

    FacetIntegration faee1(fe, elem1_, position_, false, false);
    std::vector<double> eqnFacet = faee1.equation_plane(cornersLocal);
    Eqnplane.push_back(eqnFacet);

    if (fabs(eqnFacet[0]) > 1e-10)
      XFacets.push_back(i);
    else
      NotXFacets.push_back(i);
  }

  //-------------------------------------------------------------------------
  // Step 2: Shoot a ray along x-direction and find all intersection points
  //-------------------------------------------------------------------------
  std::map<double, int> Xintersect;
  // ray intersects only XFacets
  for (unsigned i = 0; i < XFacets.size(); i++)
  {
    // check whether the infinite ray thru given (y,z) intersect the facet
    Facet *fe = *XFacets[i];
    std::vector<std::vector<double>> cornersLocal;
    fe->CornerPointsLocal(elem1_, cornersLocal);
    int cutno = pnpoly(cornersLocal, rst, GEO::CUT::proj_x);
    if (cutno == 1)
    {
      // find x-value of intersection point, (yInt,zInt) = (y,z) of given pt
      const int idx = std::distance(facete.begin(), XFacets[i]);
      std::vector<double> eqn = Eqnplane[idx];
      double x = (eqn[3] - eqn[1] * rst(1, 0) - eqn[2] * rst(2, 0)) / eqn[0];
      if (fabs((x - rst(0, 0)) / x) < 1e-8)  // pt is on one of the XFacets
        return "onBoundary";
      Xintersect[x] = i;
    }
  }

  //-------------------------------------------------------------------------------------
  // Step 3: based on relative location of intersection points w.r to given point, decide
  //-------------------------------------------------------------------------------------
  int numInter = Xintersect.size();

  // Check to see if point is in x-z or x-y oriented facets
  for (unsigned i = 0; i < NotXFacets.size(); i++)
  {
    const int idx = std::distance(facete.begin(), NotXFacets[i]);
    std::vector<double> eqn = Eqnplane[idx];
    // check pt on x-y plane facet
    if (fabs(eqn[2]) < 1e-10)  // make sure it is x-y facet
    {
      if (fabs((eqn[3] / eqn[1] - rst(1, 0)) / rst(1, 0)) < 1e-10)  // make sure pt is on same plane
      {
        Facet *fe = *NotXFacets[i];
        std::vector<std::vector<double>> cornersLocal;
        fe->CornerPointsLocal(elem1_, cornersLocal);
        int cutno = pnpoly(cornersLocal, rst, GEO::CUT::proj_y);
        if (cutno == 1)  // make sure pt is within facet area
        {
          return "onBoundary";
        }
      }
    }
    // check pt on x-z plane facet
    if (fabs(eqn[1]) < 1e-10)  // make sure it is x-z facet
    {
      if (fabs((eqn[3] / eqn[2] - rst(2, 0)) / rst(2, 0)) < 1e-10)
      {
        Facet *fe = *NotXFacets[i];
        std::vector<std::vector<double>> cornersLocal;
        fe->CornerPointsLocal(elem1_, cornersLocal);
        int cutno = pnpoly(cornersLocal, rst, GEO::CUT::proj_z);
        if (cutno == 1)
        {
          return "onBoundary";
        }
      }
    }
  }

  // if pt is not on x-y or x-z facet and no intersection --> outside
  if (numInter == 0) return "outside";

  // add given point --> to sort intersection pts along with given pt
  Xintersect[rst(0, 0)] = -1;
  numInter++;


  std::map<double, int>::iterator it = Xintersect.find(rst(0, 0));

  // all intersecting facets are right side to given pt
  if (it == Xintersect.begin()) return "outside";

  std::map<double, int>::iterator itEnd = Xintersect.end();
  --itEnd;
  // all intersecting facets are left side to given pt
  if (it == itEnd) return "outside";

  // if the no of facets right side of given pt is even number --> outside
  int rightFacets = std::distance(it, itEnd);
  if (rightFacets % 2 == 0) return "outside";
  return "inside";
}
