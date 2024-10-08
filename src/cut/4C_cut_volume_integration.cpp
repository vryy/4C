/*---------------------------------------------------------------------*/
/*! \file

\brief Integrates base functions over volume, distribute Gauss points and solve moment fitting
equations

\level 3


*----------------------------------------------------------------------*/

#include "4C_cut_volume_integration.hpp"

#include "4C_cut_boundingbox.hpp"
#include "4C_cut_enum.hpp"
#include "4C_cut_options.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
         compute the rhs of the moment fitting equations
         Integration of base functions take place inside this
*---------------------------------------------------------------------------*/
Core::LinAlg::SerialDenseVector Cut::VolumeIntegration::compute_rhs_moment()
{
  Core::LinAlg::SerialDenseVector rhs_mom(num_func_);

  const plain_facet_set &facete = volcell_->facets();

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
      if (rhs_mom.normInf() > 1e-5 && rhs_mom(0) < 0.0)
        FOUR_C_THROW("negaive volume in base function integration. is ordering of vertices right?");
    }
  }

  /*double intee = rhs_mom(56)+rhs_mom(68)+rhs_mom(72)+rhs_mom(83);
  std::cout<<std::setprecision(20)<<"the INTEGRAL = "<<intee<<"\n";*/

  // set the volume of this volumecell
  // the volume from local coordinates is converted in terms of global coordinates
  double volGlobal = 0.0;
  switch (elem1_->shape())
  {
    case Core::FE::CellType::hex8:
    {
      volGlobal = elem1_->scalar_from_local_to_global<3, Core::FE::CellType::hex8>(
          rhs_mom(0), "local_to_global");
      break;
    }
    case Core::FE::CellType::tet4:
    {
      volGlobal = elem1_->scalar_from_local_to_global<3, Core::FE::CellType::tet4>(
          rhs_mom(0), "local_to_global");
      break;
    }
    case Core::FE::CellType::wedge6:
    {
      volGlobal = elem1_->scalar_from_local_to_global<3, Core::FE::CellType::wedge6>(
          rhs_mom(0), "local_to_global");
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      volGlobal = elem1_->scalar_from_local_to_global<3, Core::FE::CellType::pyramid5>(
          rhs_mom(0), "local_to_global");
      break;
    }
    default:
      FOUR_C_THROW("unsupported integration cell type ( cell type = %s )",
          Core::FE::cell_type_to_string(elem1_->shape()).c_str());
      exit(EXIT_FAILURE);
  }
  volcell_->set_volume(volGlobal);

  return rhs_mom;
}

/*----------------------------------------------------------------------------*
    compute the gaussian points of the volumecell with "numeach" points in each
    3-directions
    numeach should be more than 1
    uses ray tracing method
*-----------------------------------------------------------------------------*/
bool Cut::VolumeIntegration::compute_gaussian_points(int numeach)
{
  Teuchos::RCP<BoundingBox> box1 = Teuchos::RCP(BoundingBox::create(*volcell_, elem1_));
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
    area = is_contain_area(minn, maxx, zmin, InPlane, zcoord, ycoord, 0.01, numeach);
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
    area = is_contain_area(minn, maxx, zmax, InPlane, zcoord, ycoord, 0.01, numeach);
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
      area = is_contain_area(minn, maxx, zmin, InPlane, zcoord, ycoord, 0.001, numeach);
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
      area = is_contain_area(minn, maxx, zmax, InPlane, zcoord, ycoord, 0.001, numeach);
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
      area = is_contain_area(minn, maxx, zplane[i], InPlane, zcoord, ycoord, 0.01, numeach);
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
          innerarea = is_contain_area(minn, maxx, zz, InPlane, zcoord, ycoord, 0.01, numeach);
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

  return wei;
}

/*-------------------------------------------------------------------------------------------------------------------*
    Store the z- and y-coordinates of the all corner points which will be used to find whether the
intersection point lies inside the volume or not
*--------------------------------------------------------------------------------------------------------------------*/
void Cut::VolumeIntegration::get_zcoordinates(
    std::vector<std::vector<double>> &zcoord, std::vector<std::vector<double>> &ycoord)
{
  const plain_facet_set &facete = volcell_->facets();
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
    face1->corner_points_local(elem1_, corLocal);
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
bool Cut::VolumeIntegration::is_intersect(double *pt, double *mini, double *maxi,
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

    // old way of point distribution method in intersection lines
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
      on_line(inter1, inter2, linePts, 2);
      intersect = true;
    }
    else if (fabs(inter2[0] - inter1[0]) < 0.2 * (maxi[0] - mini[0]))
    {
      on_line(inter1, inter2, linePts, 3);
      intersect = true;
    }
    else if (fabs(inter2[0] - inter1[0]) < 0.3 * (maxi[0] - mini[0]))
    {
      on_line(inter1, inter2, linePts, 4);
      intersect = true;
    }
    else
      on_line(inter1, inter2, linePts, numeach);
    return intersect;
  }
  else
  {
    // map is useful since we need to arrange the elements from minimum x-cut value
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

    const plain_facet_set &facete = volcell_->facets();
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
      if (facet1->on_cut_side())
      {
        if (eqn_facets_[face1][0] < 0.0) ptsInside = true;
      }
      else if (facet2->on_cut_side())
      {
        if (eqn_facets_[face2][0] > 0.0) ptsInside = true;
      }
      else
      {
        std::cout << "The assumption that one must be a cut surface is false" << std::endl;
      }

      if (ptsInside)
      {
        if (interPoints.size() == 2)  // if intersection is on the triangulated line, results only
                                      // in two face intersections
        {
          on_line(inter1, inter2, linePts, numeach);
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
          on_line(inter1, inter2, linePts, 2 /*numeach/3+1*/);
          intersect = true;
        }
        else
        {
          on_line(inter1, inter2, linePts, numeach / 2 + 1);
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
int Cut::VolumeIntegration::pnpoly(
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
int Cut::VolumeIntegration::pnpoly(const std::vector<std::vector<double>> &xp,
    const Core::LinAlg::Matrix<3, 1> &pt, Cut::ProjectionDirection projType)
{
  int npol = xp.size();
  int ind1 = 1, ind2 = 2;
  if (projType == Cut::proj_y)
  {
    ind1 = 2;
    ind2 = 0;
  }
  else if (projType == Cut::proj_z)
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
bool Cut::VolumeIntegration::is_contain_area(double minn[3], double maxx[3], double &zmin,
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
    std::array<double, 3> a = {xmin, ymin, zmin};
    intersec = is_intersect(a.data(), minn, maxx, linePts, zcoord, ycoord, toler, numeach);

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
    std::array<double, 3> a = {xmin, ymax, zmin};
    intersec = is_intersect(a.data(), minn, maxx, linePts, zcoord, ycoord, toler, numeach);
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
    std::array<double, 3> a = {xmin, yplane[i], zmin};
    bool intersec = false;
    std::vector<std::vector<double>> linePts;
    intersec = is_intersect(a.data(), minn, maxx, linePts, zcoord, ycoord, toler, numeach);
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
        std::array<double, 3> a = {xmin, yplane[i] - dyym, zmin};
        bool innerintersec = false;
        std::vector<std::vector<double>> linePts;
        innerintersec = is_intersect(a.data(), minn, maxx, linePts, zcoord, ycoord, toler, numeach);

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
void Cut::VolumeIntegration::on_line(std::vector<double> inter1, std::vector<double> inter2,
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

/*-------------------------------------------------------------------------------------*
 * Check whether the point with this element Local coordinates is inside,              *
 * outside or on boundary of this volumecell                            sudhakar 07/12 *
 *-------------------------------------------------------------------------------------*/
std::string Cut::VolumeIntegration::is_point_inside(Core::LinAlg::Matrix<3, 1> &rst)
{
  const plain_facet_set &facete = volcell_->facets();

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
    fe->corner_points_local(elem1_, cornersLocal);

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
    fe->corner_points_local(elem1_, cornersLocal);
    int cutno = pnpoly(cornersLocal, rst, Cut::proj_x);
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
        fe->corner_points_local(elem1_, cornersLocal);
        int cutno = pnpoly(cornersLocal, rst, Cut::proj_y);
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
        fe->corner_points_local(elem1_, cornersLocal);
        int cutno = pnpoly(cornersLocal, rst, Cut::proj_z);
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

FOUR_C_NAMESPACE_CLOSE
