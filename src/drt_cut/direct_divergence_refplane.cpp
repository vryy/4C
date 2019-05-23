/*!-----------------------------------------------------------------------------------------------*
\file direct_divergence_refplane.cpp

\brief Construct reference plane for direct divergence method when used in global
coordinate system

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
 *------------------------------------------------------------------------------------------------*/
#include "direct_divergence_refplane.H"
#include "cut_kernel.H"
#include "cut_side.H"
#include "cut_volumecell.H"
#include "cut_options.H"
#include "cut_output.H"

/*-----------------------------------------------------------------------------------*
 * Perform all the operations related to computing reference plane           sudhakar 06/15
 *-----------------------------------------------------------------------------------*/
std::vector<double> GEO::CUT::DirectDivergenceGlobalRefplane::GetReferencePlane()
{
  if (elem1_->Shape() != DRT::Element::hex8)
  {
    std::cout << "Element Shape: " << elem1_->Shape() << std::endl;
    throw std::runtime_error("Currently can handle only hexagonal family\n");
    // dserror("Currently can handle only hexagonal family\n");
  }

  std::vector<double> RefPlaneEqn(4, 0.0);

  bool comp_ref_plane = false;

  // In the first round check if the projection of element corner points is inside
  std::vector<Point*> points = elem1_->Points();

  static const int max_attempts = 9;
  double tol = 1e-8;

  for (int refplaneattempt = 0; refplaneattempt < max_attempts; ++refplaneattempt)
  {
    //---
    // First estimate -- Compute reference plane based on the facets of the volumecell information
    //---
    comp_ref_plane = FacetBasedRef(RefPlaneEqn, points, tol);
    if (comp_ref_plane)
    {
      return RefPlaneEqn;
    }
    refPtsGmsh_.clear();

    //---
    // Second estimate -- Compute reference plane based on the diagonal information
    //---
    comp_ref_plane = DiagonalBasedRef(RefPlaneEqn, points, tol);
    if (comp_ref_plane)
    {
      return RefPlaneEqn;
    }
    refPtsGmsh_.clear();

    //---
    // Third estimate -- Compute reference plane based on Side information
    //---
    comp_ref_plane = SideBasedRef(RefPlaneEqn, points, tol);
    if (comp_ref_plane)
    {
      return RefPlaneEqn;
    }

    if (not comp_ref_plane)
    {
      std::stringstream str;
      str << ".no_refplane_" << refplaneattempt << "_CUTFAIL.pos";
      std::string filename(GEO::CUT::OUTPUT::GenerateGmshOutputFilename(str.str()));
      std::ofstream file(filename.c_str());

      GEO::CUT::OUTPUT::GmshCompleteCutElement(file, elem1_);
      GEO::CUT::OUTPUT::GmshNewSection(file, "VolumeCell");
      GEO::CUT::OUTPUT::GmshVolumecellDump(file, volcell_);
      GEO::CUT::OUTPUT::GmshEndSection(file, true);
    }

    if (!refplaneattempt)
    {
      // We couldn't find a reference plane where gauss points from the whole element are inside.
      // But anyway this is not required, since we are just integrating a vc here --> therefore
      // project all points from the vc to check if it is a good reference plane --> this is more
      // expensive depending on the number of points describing the vc but will just be neccesary in
      // some cases ...
      points.clear();
      for (plain_facet_set::const_iterator fit = volcell_->Facets().begin();
           fit != volcell_->Facets().end(); ++fit)
      {
        for (uint p = 0; p < (*fit)->Points().size(); ++p)
        {
          bool insert = true;
          for (uint ap = 0; ap < points.size(); ++ap)
            if (points[ap]->Id() == (*fit)->Points()[p]->Id())
            {
              insert = false;
              break;
            }
          if (insert) points.push_back((*fit)->Points()[p]);
        }
      }
    }
    else
    {
      tol *= 10;
      std::cout << "Warning: Increasing my tolerance to find a proper reference plane to tol = "
                << tol << std::endl;
    }
  }
  dserror("Proper reference plane not found");

  return RefPlaneEqn;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the diagonal of background element sudhakar 06/15 This
 *considers all 6 diagonals of background Hex element, and choose the one that has maximum normal
 *component in x-direction
 *-------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::DirectDivergenceGlobalRefplane::DiagonalBasedRef(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  if (options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_all &&
      options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_diagonal &&
      options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_diagonal_side)
    return false;

  // Any method should involve the following two steps

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh
  // output. Here we take all possible 6 diagonals of the Hex element, and choose the diagonal which
  // has maximum normal component in x-direction as the reference plane
  //---
  std::vector<Point*> ptslist = elem1_->Points();

  std::vector<Point*> diag;
  std::vector<std::vector<Point*>> diagonals;
  for (unsigned i = 0; i < 24; ++i)  // 24 is the number of considered diagonals in
                                     // tri_diags_[24][3]
  {
    diag.clear();
    for (uint plid = 0; plid < 3; ++plid) diag.push_back(ptslist[tri_diags_[i][plid]]);
    diagonals.push_back(diag);
  }

  double xnormal = 0.0;
  bool found_refplane = false;
  for (std::vector<std::vector<Point*>>::iterator itd = diagonals.begin(); itd != diagonals.end();
       itd++)
  {
    std::vector<Point*> ptl = *itd;
    std::vector<double> RefPlaneTemp = KERNEL::EqnPlaneOfPolygon(ptl);

    scaleEquationOfPlane(RefPlaneTemp);

    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    //---
    // STEP 2: Project all the corner points of the Hex element onto the reference plane.
    // If all these projected points are within the background element, then we consider this as a
    // possible reference plane!
    if (isAllProjectedCornersInsideEle(RefPlaneTemp, points, tol))
    {
      double fac =
          sqrt(pow(RefPlaneTemp[0], 2) + pow(RefPlaneTemp[1], 2) + pow(RefPlaneTemp[2], 2));
      //---
      // STEP 3: Take the reference plane with biggest component in x-direction of the normal
      // vector!
      double xn = fabs(RefPlaneTemp[0]) / fac;

      if (xn > xnormal)
      {
        xnormal = xn;
        RefPlaneEqn = RefPlaneTemp;
        refPtsGmsh_ = ptl;
        found_refplane = true;
      }
    }
  }

  return found_refplane;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the facets of the volumecell                         ager
 *02/16 Sort all the sides based on n_x (the one has more n_x gets on the top) Iterate through all
 *the sides to get the correct reference plane
 *-------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::DirectDivergenceGlobalRefplane::FacetBasedRef(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  if (options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_all &&
      options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_facet)
    return false;

  const plain_facet_set& allfacets = volcell_->Facets();

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh
  // output. First get all the facets of the volumecell and compute the equation of plane for each
  // facet Store them in a data structure which stores the sides based on n_x in descending order
  //---
  std::multimap<double, std::pair<std::vector<double>, std::vector<Point*>>, compareClass>
      facet_data;
  for (plain_facet_set::const_iterator it = allfacets.begin(); it != allfacets.end(); it++)
  {
    std::vector<double> RefPlaneTemp = KERNEL::EqnPlaneOfPolygon((*it)->Points());
    scaleEquationOfPlane(RefPlaneTemp);
    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    facet_data.insert(
        std::make_pair(RefPlaneTemp[0], std::make_pair(RefPlaneTemp, (*it)->Points())));
  }

  double xnormal = 0.0;
  bool found_refplane = false;
  for (std::multimap<double, std::pair<std::vector<double>, std::vector<Point*>>>::iterator it =
           facet_data.begin();
       it != facet_data.end(); it++)
  {
    std::vector<double> RefPlaneTemp = it->second.first;

    scaleEquationOfPlane(RefPlaneTemp);

    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    //---
    // STEP 2: Project all the corner points of the Hex element onto the reference plane.
    // If all these projected points are within the background element, then we consider this as a
    // possible reference plane!

    if (isAllProjectedCornersInsideEle(RefPlaneTemp, points, tol))
    {
      double fac =
          sqrt(pow(RefPlaneTemp[0], 2) + pow(RefPlaneTemp[1], 2) + pow(RefPlaneTemp[2], 2));
      //---
      // STEP 3: Take the reference plane with biggest component in x-direction of the normal
      // vector!
      double xn = fabs(RefPlaneTemp[0]) / fac;
      if (xn > xnormal)
      {
        xnormal = xn;
        RefPlaneEqn = RefPlaneTemp;
        refPtsGmsh_ = it->second.second;
        ;
        found_refplane = true;
      }
    }
  }

  //---
  // Basically using a diagonal reference plane should be enought, otherwise we have to look into
  // that again!
  return found_refplane;
}

/*-------------------------------------------------------------------------------------------------------*
 * Computation of reference plane based on the sides of background element sudhakar 06/15 Sort all
 *the sides based on n_x (the one has more n_x gets on the top) Iterate through all the sides to get
 *the correct reference plane
 *-------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::DirectDivergenceGlobalRefplane::SideBasedRef(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  if (options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_all &&
      options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_side &&
      options_.Direct_Divergence_Refplane() != INPAR::CUT::DirDiv_refplane_diagonal_side)
    return false;

  const std::vector<Side*>& allsides = elem1_->Sides();

  //---
  // STEP 1: Estimate equation of reference plane and store the corresponding points for gmsh
  // output. First get all the sides of the element and compute the equation of plane for each side
  // Store them in a data structure which stores the sides based on n_x in descending order
  //---
  std::multimap<double, std::pair<std::vector<double>, std::vector<Point*>>, compareClass>
      side_data;
  for (std::vector<Side*>::const_iterator it = allsides.begin(); it != allsides.end(); it++)
  {
    const Side* s = *it;
    const std::vector<Node*> nds = s->Nodes();


    for (uint split_quadidx = 0; split_quadidx < 3 * nds.size() - 8; ++split_quadidx)
    {
      std::vector<Point*> ptside;
      if (nds.size() == 4)
      {
        ptside.push_back(nds[side_split_[split_quadidx][0]]->point());
        ptside.push_back(nds[side_split_[split_quadidx][1]]->point());
        ptside.push_back(nds[side_split_[split_quadidx][2]]->point());
      }
      else if (nds.size() == 3)
      {
        for (std::vector<Node*>::const_iterator itn = nds.begin(); itn != nds.end(); itn++)
          ptside.push_back((*itn)->point());
      }
      else
        dserror("Side with another number of nodes than 3 or 4?");

      std::vector<double> RefPlaneTemp = KERNEL::EqnPlaneOfPolygon(ptside);
      scaleEquationOfPlane(RefPlaneTemp);
      if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

      side_data.insert(std::make_pair(RefPlaneTemp[0], std::make_pair(RefPlaneTemp, ptside)));
    }
  }

  double xnormal = 0.0;
  bool found_refplane = false;
  for (std::multimap<double, std::pair<std::vector<double>, std::vector<Point*>>>::iterator it =
           side_data.begin();
       it != side_data.end(); it++)
  {
    std::vector<double> RefPlaneTemp = it->second.first;

    scaleEquationOfPlane(RefPlaneTemp);

    if (fabs(RefPlaneTemp[0]) < REF_PLANE_DIRDIV) continue;

    //---
    // STEP 2: Project all the corner points of the Hex element onto the reference plane.
    // If all these projected points are within the background element, then we consider this as a
    // possible reference plane!
    if (isAllProjectedCornersInsideEle(RefPlaneTemp, points, tol))
    {
      double fac =
          sqrt(pow(RefPlaneTemp[0], 2) + pow(RefPlaneTemp[1], 2) + pow(RefPlaneTemp[2], 2));
      //---
      // STEP 3: Take the reference plane with biggest component in x-direction of the normal
      // vector!
      double xn = fabs(RefPlaneTemp[0]) / fac;
      if (xn > xnormal)
      {
        xnormal = xn;
        RefPlaneEqn = RefPlaneTemp;
        refPtsGmsh_ = it->second.second;
        ;
        found_refplane = true;
      }
    }
  }

  //---
  // Basically using a diagonal reference plane should be enought, otherwise we have to look into
  // that again!
  return found_refplane;
}

/*---------------------------------------------------------------------------------------------------------------*
 * In order to check whether the chosen reference plane is a correct choice,
 * we project all the corner points of the element onto this reference plane sudhakar 06/15 If all
 *these projected points are within the element, then we got the correct ref plane
 *---------------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::DirectDivergenceGlobalRefplane::isAllProjectedCornersInsideEle(
    std::vector<double>& RefPlaneEqn, std::vector<Point*> points, double tol)
{
  for (std::vector<Point*>::iterator it = points.begin(); it != points.end(); it++)
  {
    Point* pt = *it;

    LINALG::Matrix<3, 1> coo;
    pt->Coordinates(coo.A());

    LINALG::Matrix<3, 1> xyz_proj(coo), rst_proj;

    // x-coordinate of corner pt projected in the reference plane
    // y- and z-coordinate remain the same
    xyz_proj(0, 0) =
        (RefPlaneEqn[3] - RefPlaneEqn[1] * coo(1, 0) - RefPlaneEqn[2] * coo(2, 0)) / RefPlaneEqn[0];

    // get the local coordinates of the projected point
    elem1_->LocalCoordinates(xyz_proj, rst_proj);

    // Check whether the local coordinate of the projected point is within the specified limits
    if (std::abs(rst_proj(0, 0)) > 1.0 + tol or std::abs(rst_proj(1, 0)) > 1.0 + tol or
        std::abs(rst_proj(2, 0)) > 1.0 + tol or std::isnan(rst_proj(0, 0)) or
        std::isnan(rst_proj(1, 0)) or std::isnan(rst_proj(2, 0)))
    {
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------------------------------------------*
 * Scale the equation of plane to enable comparison of normals sudhakar 07/15
 *----------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergenceGlobalRefplane::scaleEquationOfPlane(
    std::vector<double>& RefPlaneEqn)
{
  double scale =
      sqrt(pow(RefPlaneEqn[0], 2.0) + pow(RefPlaneEqn[1], 2.0) + pow(RefPlaneEqn[2], 2.0));
  for (unsigned i = 0; i < 4; i++) RefPlaneEqn[i] /= scale;
}

const unsigned GEO::CUT::DirectDivergenceGlobalRefplane::tri_diags_[24][3] = {{1, 6, 7}, {0, 6, 7},
    {0, 1, 7}, {0, 1, 6}, {3, 4, 5}, {2, 4, 5}, {2, 3, 5}, {2, 3, 4}, {6, 3, 0}, {5, 3, 0},
    {5, 6, 0}, {5, 6, 3}, {7, 2, 1}, {4, 2, 1}, {4, 7, 1}, {4, 7, 2}, {4, 6, 2}, {0, 6, 2},
    {0, 4, 2}, {0, 4, 6}, {1, 3, 7}, {5, 3, 7}, {5, 1, 7}, {5, 1, 3}};

const unsigned GEO::CUT::DirectDivergenceGlobalRefplane::side_split_[4][3] = {
    {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
