/*---------------------------------------------------------------------*/
/*!
\file direct_divergence.cpp

\brief Generate main Gauss points when using "DirectDivergence" approach.

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/


#include "direct_divergence.H"
#include "facet_integration.H"
#include "volume_integration.H"

#include "cut_kernel.H"
#include "direct_divergence_refplane.H"

#include "cut_output.H"

#include <Teuchos_TimeMonitor.hpp>


/*-------------------------------------------------------------------------------------------------------------------*
  Create integration points on the facets of the volumecell by triangulating the facets
  A reference facet is identified on which integration weights are set to zero Sudhakar 04/12
*--------------------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::DirectDivergence::VCIntegrationRule(
    std::vector<double>& RefPlaneEqn)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::DirectDivergence::VCIntegrationRule" );

  std::vector<plain_facet_set::const_iterator>
      facetIterator;  // iterators of facets which need to be considered for integration rule
  plain_facet_set::const_iterator IteratorRefFacet;  // iterator for the reference facet
  isRef_ = false;                                    // whether ref plane is falling on facet?

  // get integration facets and reference plane
  Teuchos::RCP<BoundingBox> fbox = Teuchos::rcp(BoundingBox::Create());
  const plain_facet_set& facete = volcell_->Facets();
  // create bounding box around all facets
  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fe = *i;
    const std::vector<Point*>& corn = fe->Points();
    for (std::vector<Point*>::const_iterator p = corn.begin(); p != corn.end(); ++p)
    {
      fbox->AddPoint((*p)->X());
    }
  }
  LINALG::Matrix<3, 2> fvolume = fbox->GetBoundingVolume();
  const double totalVolume = (fvolume(0, 1) - fvolume(0, 0)) * (fvolume(1, 1) - fvolume(1, 0)) *
                             (fvolume(2, 1) - fvolume(2, 0));

  ListFacets(facetIterator, RefPlaneEqn, IteratorRefFacet, isRef_);

  if (isRef_) refFacet_ = *IteratorRefFacet;

  if (facetIterator.size() == 0)
  {
#if EXTENDED_CUT_DEBUG_OUTPUT
    std::cout << "NOTICE: X-component normal is zero on all facets,\
                 Volume of the bounding box surrounding facets is equal to "
              << totalVolume << std::endl;
#endif

    if (totalVolume > MINIMUM_VOLUME_BB_FACETS)
    {
      const plain_facet_set& facete = volcell_->Facets();
      std::cout << "number of facets: " << facete.size() << std::endl;

      for (plain_facet_set::const_iterator f = facete.begin(); f != facete.end(); f++)
        (*f)->Print(std::cout);

      // dump element and facets
      std::ofstream file("facets_x_normal_equal_0_CUTFAIL_DD.pos");
      volcell_->DumpGmsh(file);
      file.close();

      std::string filename1("element_x_normal_equal_0_CUTFAIL_DD.pos");
      std::ofstream file1(filename1.c_str());
      GEO::CUT::OUTPUT::GmshCompleteCutElement(file1, elem1_, false);
      file1.close();

      std::stringstream err_msg;
      err_msg << "x-component normal is zero on all the facets? It should not be. Volume of the "
                 "bounding box around facets is "
              << totalVolume;
      dserror(err_msg.str());
    }
  }

  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp =
      Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints(0));

#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
  std::cout << "Number of facets: " << volcell_->Facets().size() << std::endl;
  std::cout << "Size of facetIterator: " << facetIterator.size() << std::endl;
#endif

  for (unsigned i = 0; i < facetIterator.size(); i++)
  {
    plain_facet_set::const_iterator iter = facetIterator[i];
    Facet* fe = *iter;
    FacetIntegration faee1(fe, elem1_, position_, false, false);

    faee1.DivergenceIntegrationRuleNew(mesh_, cgp);
  }

#if 0  // integrate specified functions using the Gaussian rule generated -- used in postprocessing
  DRT::UTILS::GaussIntegration gi(cgp);
  IntegrateSpecificFuntions( gi, RefPlaneEqn );  //integrate specific functions
#endif

  return cgp;
}

/*-------------------------------------------------------------------------------------------------------*
        Identify the list of facets which need to be triangulated (called as integration facets)
        Get the reference facet that will be used to get internal Gauss rule (sudhakar 04/12) As far
as possible, the reference facet is set on a cut side to reduce the no of Gauss pts Reference facet
is selected so that all internal points are inside the volumecell
*--------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::ListFacets(
    std::vector<plain_facet_set::const_iterator>& facetIterator, std::vector<double>& RefPlaneEqn,
    plain_facet_set::const_iterator& IteratorRefFacet, bool& IsRefFacet)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::DirectDivergence::ListFacets" );

  const plain_facet_set& facete = volcell_->Facets();

  bool RefOnCutSide = false;
  std::vector<std::vector<double>> eqnAllFacets(facete.size());

  // store the iterators of  all warped facets
  std::vector<plain_facet_set::const_iterator> warpFac;

  // check whether all facets of this vc are oriented in a plane
  // if not then some sides are warped
  // we need to generate quadrature rules in global coordinates
  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fe = *i;
    const std::vector<Point*>& corn = fe->CornerPoints();
    bool isPlanar = fe->IsPlanar(mesh_, corn);  // triangulates non-planar facets.

    if (isPlanar == false)  // and !(fe->BelongsToLevelSetSide()) )
    {
      warpFac.push_back(i);
#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
      std::cout << "encountered a WARPED side\n";

      std::cout << "side-Id() " << fe->SideId() << std::endl;

      std::cout << "the side has " << corn.size() << " points" << std::endl;

      for (unsigned j = 0; j < corn.size(); j++) corn[j]->Print(std::cout);

      std::cout << std::endl;
#endif
    }
  }

  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fe = *i;

#ifdef LOCAL
    std::vector<std::vector<double>> cornersLocal;
    fe->CornerPointsLocal(elem1_, cornersLocal, true);
#else
    std::vector<std::vector<double>> cornersLocal = fe->CornerPointsGlobal(elem1_, true);
#endif

    std::vector<double> RefPlaneTemp = KERNEL::EqnPlaneOfPolygon(cornersLocal);
    const int index = std::distance(facete.begin(), i);
    eqnAllFacets[index] = RefPlaneTemp;

    // consider only facet whose x-direction normal componenet is non-zero
    if (fabs(RefPlaneTemp[0]) > TOL_EQN_PLANE)  // This could give issues with non-planar facets?
    {
      TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT::DirectDivergence::ListFacets-tmp1");

#ifdef LOCAL
      if (warpFac.size() > 0)  // if there are warped facets that are not yet processed
      {
        if (i == warpFac[0])
        {
          // reference plane cant be defined over warped facet since the facet itself is not in a
          // plane
          facetIterator.push_back(i);
          warpFac.erase(warpFac.begin());
          continue;
        }
      }

      // store the non-cut facet as reference
      if (isRef_ == false && !fe->OnCutSide() && fabs(RefPlaneTemp[0]) > REF_PLANE_DIRDIV)
      {
        RefPlaneEqn = RefPlaneTemp;
        IteratorRefFacet = i;
        isRef_ = true;
        continue;
      }
      // as far as possible, take cut side as reference
      // because it it possible to delet one or more from facetList
      if (RefOnCutSide == false && fe->OnCutSide() && fabs(RefPlaneTemp[0]) > REF_PLANE_DIRDIV)
      {
        bool addRef = false;

        // if any cut side is on the plane x=C definitely this can be a reference facet
        // this produces all internal points within background element
        if (fabs(RefPlaneTemp[1]) < TOL_EQN_PLANE && fabs(RefPlaneTemp[2]) < TOL_EQN_PLANE)
        {
          addRef = true;
        }

        // for any other cut facet, make sure when it is extended over the whole
        // volumecell, it falls within the background element --> inside internal points
        else
        {
          double x1 = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0;
          double y1 = -1.0, y2 = 1.0, z1 = -1.0, z2 = 1.0;
          x1 = (RefPlaneTemp[3] - RefPlaneTemp[2] * z1 - RefPlaneTemp[1] * y1) / RefPlaneTemp[0];
          x2 = (RefPlaneTemp[3] - RefPlaneTemp[2] * z1 - RefPlaneTemp[1] * y2) / RefPlaneTemp[0];
          x3 = (RefPlaneTemp[3] - RefPlaneTemp[2] * z2 - RefPlaneTemp[1] * y1) / RefPlaneTemp[0];
          x4 = (RefPlaneTemp[3] - RefPlaneTemp[2] * z2 - RefPlaneTemp[1] * y2) / RefPlaneTemp[0];

          // TODO: This is specific to hex8 background element. Extend this for other elements also
          if (fabs(x1) < (1.0 + 1e-8) && fabs(x2) < (1.0 + 1e-8) && fabs(x3) < (1.0 + 1e-8) &&
              fabs(x4) < (1.0 + 1e-8))
          {
            addRef = true;
          }
        }

        if (addRef)
        {
          if (isRef_)  // push already existing reference facet from non-cut side
            facetIterator.push_back(IteratorRefFacet);
          RefPlaneEqn = RefPlaneTemp;
          IteratorRefFacet = i;
          isRef_ = true;
          RefOnCutSide = true;
          continue;
        }
      }
#endif
      facetIterator.push_back(i);  // if not a referece side, push this to facet consideration
    }
  }

#ifdef LOCAL
  // if no reference side found --> no element side and no cut side with x=C
  // this means we create a plane x=0 and assumes this as reference plane
  if (isRef_ == false)
  {
    RefPlaneEqn[0] = 1.0;
    for (unsigned i = 1; i < 4; i++) RefPlaneEqn[i] = 0.0;
  }
#else
  // When we construct integration rule in global coordinate system,
  // we need a reference plane such that when the main Gauss points are projected over
  // this plane, this line of projection must be completely within the background element
  DirectDivergenceGlobalRefplane ddg(elem1_, volcell_, mesh_.GetOptions());
  RefPlaneEqn = ddg.GetReferencePlane();
  refPtsGmsh_ = ddg.GetReferencePointGmsh();
#endif

  // if a1x+a2y+a3z=a4 is the equation of reference plane and
  //    b1x+b2y+b3z=b4 is equation of the considered facet
  // if (a1/a4==b1/b4 && a2/a4==b2/b4 && a3/a4==b3/b4 ) then both reference and the
  //   considered facet are in the same plane, so delete this facet
  if (RefOnCutSide)
  {
    // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::DirectDivergence::ListFacets-tmp2" );

    for (unsigned i = 0; i < facetIterator.size(); i++)
    {
      plain_facet_set::const_iterator iter = facetIterator[i];

      const int index = std::distance(facete.begin(), iter);

      double facetx = eqnAllFacets[index][0];
      double facety = eqnAllFacets[index][1];
      double facetz = eqnAllFacets[index][2];
      double facetRhs = eqnAllFacets[index][3];

      if (fabs(RefPlaneEqn[3]) > TOL_EQN_PLANE &&
          fabs(facetRhs) > TOL_EQN_PLANE)  // planes for which ax+by+cz=d
      {
        if (fabs(RefPlaneEqn[0] / RefPlaneEqn[3] - facetx / facetRhs) < TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[1] / RefPlaneEqn[3] - facety / facetRhs) < TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[2] / RefPlaneEqn[3] - facetz / facetRhs) < TOL_EQN_PLANE)
        {
          facetIterator.erase(facetIterator.begin() + i);
          i--;
        }
      }
      else  // planes for which ax+by+cz=0
      {
        if (fabs(RefPlaneEqn[1] / RefPlaneEqn[0] - facety / facetx) < TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[2] / RefPlaneEqn[0] - facetz / facetx) < TOL_EQN_PLANE &&
            fabs(RefPlaneEqn[3] / RefPlaneEqn[0] - facetRhs / facetx) < TOL_EQN_PLANE)
        {
          facetIterator.erase(facetIterator.begin() + i);
          i--;
        }
      }
    }
  }
}

/*--------------------------------------------------------------------------------------------------------------*
                   Geometry of volumecell and main Gauss pts for visualization sudhakar 04/12
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::DivengenceCellsGMSH(
    const DRT::UTILS::GaussIntegration& gpv, Teuchos::RCP<DRT::UTILS::GaussPoints>& gpmain)
{
#ifdef LOCAL
  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "divergenceCells" << sideno << "_" << this->elem1_->Id() << ".pos";
  std::ofstream file(str.str().c_str());

  volcell_->DumpGmsh(file);

  // Activate this if you think that something is wrong with the shape of vc
  // But it is always in global coordinates
  // volcell_->DumpGmshSolid( file, mesh_ );

#ifdef OUTPUT_GLOBAL_DIVERGENCE_CELLS
  const plain_facet_set& facete = volcell_->Facets();
  file << "View \"EqnOfPlaneNormals \" {\n";
  for (plain_facet_set::const_iterator j = facete.begin(); j != facete.end(); j++)
  {
    // Writes only for GLOBAL Coordinates!
    GEO::CUT::OUTPUT::GmshEqnPlaneNormalDump(file, *j, true);
  }
  file << "};\n";
#endif

  //-----------------------------------------------------------------------
  // write main Gauss points
  file << "Geometry.PointSize=8.0;\n";  // Increase the point size
  file << "View \"Main points \" {\n";
  for (DRT::UTILS::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
  {
    const LINALG::Matrix<3, 1> etaFacet(iquad.Point());
#ifndef OUTPUT_GLOBAL_DIVERGENCE_CELLS
    file << "SP(" << etaFacet(0, 0) << "," << etaFacet(1, 0) << "," << etaFacet(2, 0) << ","
         << "1"
         << "){0.0};" << std::endl;
#else
    LINALG::Matrix<3, 1> etaGlobal;
    elem1_->GlobalCoordinates(etaFacet, etaGlobal);
    file << "SP(" << etaGlobal(0, 0) << "," << etaGlobal(1, 0) << "," << etaGlobal(2, 0) << ","
         << "1"
         << "){0.0};" << std::endl;
#endif
  }
  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,100,0} };\n";  // Changing color to red
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                     // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";                 // Disable legend
  file << "View[PostProcessing.NbViews-1].PointSize = 6.0;";               // increase point size

  //-----------------------------------------------------------------------
  // write internal Gauss points
  /*file<<"Geometry.PointSize=8.0;\n";      // Increase the point size
  file<<"View \"Internal points \" {\n";
  int nu = 0;
  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    DRT::UTILS::GaussIntegration gi = intGRule[nu];
    for ( DRT::UTILS::GaussIntegration::iterator iqu=gi.begin(); iqu!=gi.end(); ++iqu )
    {
      const LINALG::Matrix<3,1> eta( iqu.Point() );
      file<<"SP("<<eta(0,0)<<","<<eta(1,0)<<","<<eta(2,0)<<","<<"1"<<"){0.0};"<<std::endl;
    }
    nu++;
  }
  file<<"};\n";
  file<<"View[PostProcessing.NbViews-1].ColorTable = { {0,0,0} };\n"; // Changing color to red
  file<<"View[PostProcessing.NbViews-1].Light=0;\n";  // Disable the lighting
  file<<"View[PostProcessing.NbViews-1].ShowScale=0;\n";  // Disable legend
  file<<"View[PostProcessing.NbViews-1].PointSize = 4.0;"; // increase point size*/

  //---------------------------------------------------------------------------------------
  // write the coordinates of the reference facet, and change its color for better visualization
  file << "View \"Ref Facet \" {\n";
  // ref plane is on a facet
  if (isRef_)
  {
    std::vector<std::vector<double>> corners;
#ifndef OUTPUT_GLOBAL_DIVERGENCE_CELLS
    refFacet_->CornerPointsLocal(elem1_, corners);
#else
    corners = refFacet_->CornerPointsGlobal(elem1_);
#endif
    //    std::vector<std::vector<double> > corners = refFacet_->CornerPointsGlobal(elem1_);
    for (unsigned i = 0; i < corners.size(); i++)
    {
      const std::vector<double> coords1 = corners[i];
      const std::vector<double> coords2 = corners[(i + 1) % corners.size()];
      file << "SL(" << coords1[0] << "," << coords1[1] << "," << coords1[2] << "," << coords2[0]
           << "," << coords2[1] << "," << coords2[2] << ")"
           << "{0,0};\n";
    }
  }
  else  // ref plane is x=0
  {
    file << "SL(0,-1,-1,0,1,-1)"
         << "{0,0};\n";
    file << "SL(0,1,-1,0,1,1)"
         << "{0,0};\n";
    file << "SL(0,1,1,0,-1,1)"
         << "{0,0};\n";
    file << "SL(0,-1,1,0,-1,-1)"
         << "{0,0};\n";
  }
  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {255,0,0} };\n";  // Changing color to red
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                     // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";                 // Disable legend
  file << "View[PostProcessing.NbViews-1].LineWidth = 4.0;";               // increase line width*/
#else
  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "divergenceCells" << sideno << ".pos";
  std::ofstream file(str.str().c_str());

  GEO::CUT::OUTPUT::GmshCompleteCutElement(file, elem1_);
  volcell_->DumpGmsh(file);

  // Activate this if you doubt that something is wrong with the vc
  // volcell_->DumpGmshSolid( file, mesh_ );

  /*//-----------------------------------------------------------------------
  // Line joining main Gauss points to reference point
  file<<"Geometry.LineWidth=1.0;\n";
  file<<"View \"Connecting lines \" {\n";
  DRT::UTILS::GaussIntegration grule(gpmain);
  for ( DRT::UTILS::GaussIntegration::iterator iquad=grule.begin(); iquad!=grule.end(); ++iquad )
  {
    const LINALG::Matrix<3,1> eta( iquad.Point() );

    double jac = fabs(refpt(0,0)-eta(0,0))*0.5; // jacobian for 1D transformation rule
    if ( jac < JAC_LINE_TOL )
      continue;

    file<<"SL("<<eta(0,0)<<","<<eta(1,0)<<","<<eta(2,0)<<","
        <<refpt(0,0)<<","<<refpt(1,0)<<","<<refpt(2,0)<<"){0.0,0.0};"<<std::endl;
  }
  file<<"};\n";
  file<<"View[PostProcessing.NbViews-1].ColorTable = { {0,0,100} };\n"; // Changing color to red
  file<<"View[PostProcessing.NbViews-1].Light=0;\n";  // Disable the lighting
  file<<"View[PostProcessing.NbViews-1].ShowScale=0;\n";  // Disable legend
  file<<"View[PostProcessing.NbViews-1].PointSize = 6.0;"; // increase point size
  //-----------------------------------------------------------------------*/

  //-----------------------------------------------------------------------
  // Draw diagonal reference side
  file << "Geometry.LineWidth=1.0;\n";
  file << "View \"Diagonal reference side \" {\n";

  for (unsigned itp = 0; itp != refPtsGmsh_.size(); itp++)
  {
    Point* pt1 = refPtsGmsh_[itp];
    Point* pt2 = refPtsGmsh_[(itp + 1) % refPtsGmsh_.size()];
    double co1[3], co2[3];
    pt1->Coordinates(co1);
    pt2->Coordinates(co2);
    file << "SL(" << co1[0] << "," << co1[1] << "," << co1[2] << "," << co2[0] << "," << co2[1]
         << "," << co2[2] << "){0.0,0.0};" << std::endl;
  }

  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,0,100} };\n";  // Changing color to red
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                     // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";                 // Disable legend
  file << "View[PostProcessing.NbViews-1].PointSize = 6.0;";               // increase point size
  //-----------------------------------------------------------------------

  //-----------------------------------------------------------------------
  // write main Gauss points
  file << "Geometry.PointSize=8.0;\n";  // Increase the point size
  file << "View \"Main points \" {\n";
  if (gpv.NumPoints() > 0)
  {
    for (DRT::UTILS::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
    {
      const LINALG::Matrix<3, 1> etaFacet(iquad.Point());
      file << "SP(" << etaFacet(0, 0) << "," << etaFacet(1, 0) << "," << etaFacet(2, 0) << ","
           << "1"
           << "){0.0};" << std::endl;
    }
  }
  file << "};\n";
  file << "View[PostProcessing.NbViews-1].ColorTable = { {0,100,0} };\n";  // Changing color to red
  file << "View[PostProcessing.NbViews-1].Light=0;\n";                     // Disable the lighting
  file << "View[PostProcessing.NbViews-1].ShowScale=0;\n";                 // Disable legend
  file << "View[PostProcessing.NbViews-1].PointSize = 6.0;";               // increase point size
  //-----------------------------------------------------------------------

#endif
}

/*--------------------------------------------------------------------------------------------------------------*
     Compute the volume of the considered cell by integrating 1 using the Gauss rule obtained.
sudhakar 04/12 Then the volume in local coordinates is converted to global coordinate value
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::DebugVolume(const DRT::UTILS::GaussIntegration& gpv, bool& isNeg)
{
  int numint = 0;
  double TotalInteg = 0.0;

  for (DRT::UTILS::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
  {
    const LINALG::Matrix<3, 1> etaFacet(iquad.Point());
    const double weiFacet = iquad.Weight();

    TotalInteg += weiFacet;
    numint++;
  }

  // set the volume of this volumecell
  // the volume from local coordinates is converted in terms of global coordinates
  double volGlobal = 0.0;
#ifdef LOCAL
  if (not elem1_->isShadow())
  {
    switch (elem1_->Shape())
    {
      case DRT::Element::hex8:
      {
        volGlobal =
            elem1_->ScalarFromLocalToGlobal<3, DRT::Element::hex8>(TotalInteg, "LocalToGlobal");
        break;
      }
      case DRT::Element::tet4:
      {
        volGlobal =
            elem1_->ScalarFromLocalToGlobal<3, DRT::Element::tet4>(TotalInteg, "LocalToGlobal");
        break;
      }
      case DRT::Element::wedge6:
      {
        volGlobal =
            elem1_->ScalarFromLocalToGlobal<3, DRT::Element::wedge6>(TotalInteg, "LocalToGlobal");
        break;
      }
      case DRT::Element::pyramid5:
      {
        volGlobal =
            elem1_->ScalarFromLocalToGlobal<3, DRT::Element::pyramid5>(TotalInteg, "LocalToGlobal");
        break;
      }
      default:
        throw std::runtime_error("unsupported element type");
    }
  }

  else
  {
    switch (elem1_->getQuadShape())
    {
      case DRT::Element::hex20:
      {
        volGlobal = elem1_->ScalarFromLocalToGlobal<3, DRT::Element::hex20>(
            TotalInteg, "LocalToGlobal", true);
        break;
      }
      case DRT::Element::hex27:
      {
        volGlobal = elem1_->ScalarFromLocalToGlobal<3, DRT::Element::hex27>(
            TotalInteg, "LocalToGlobal", true);
        break;
      }
      case DRT::Element::tet10:
      {
        volGlobal = elem1_->ScalarFromLocalToGlobal<3, DRT::Element::tet10>(
            TotalInteg, "LocalToGlobal", true);
        break;
      }
      default:
        throw std::runtime_error("unsupported parent Quad element type");
    }
  }
#else
  volGlobal = TotalInteg;
#endif

  if (volGlobal < 0.0 || TotalInteg < 0.0)
  {
    if (fabs(TotalInteg) < REF_VOL_DIRDIV)
    {
      isNeg = true;
      volcell_->SetVolume(0.0);
      std::cout << "----WARNING:::negligible volumecell parent id = "
                << volcell_->ParentElement()->Id() << "---------------" << std::endl;
      std::cout << "volume in local coordinates = " << TotalInteg
                << "\t volume in global coordinates = " << volGlobal << std::endl;

      return;
    }

    mesh_.DebugDump(elem1_, __FILE__, __LINE__);
    std::cout << "volume: " << TotalInteg << std::endl;
    throw std::runtime_error("negative volume predicted by the DirectDivergence integration rule;");
    //    dserror("negative volume predicted by the DirectDivergence integration rule; volume =
    //    %0.20f",TotalInteg);
  }

#ifdef DEBUGCUTLIBRARY  // check the volume with the moment fitting and check the values
#ifdef LOCAL
#ifdef DIRECTDIV_EXTENDED_DEBUG_OUTPUT
  std::cout << "Check volumecell with moment-fitting." << std::endl;
#endif
  VolumeIntegration vi(volcell_, elem1_, volcell_->Position(), 1);
  Epetra_SerialDenseVector volMom = vi.compute_rhs_moment();
  volMom(0) = volcell_->Volume();

  if (fabs(volGlobal - volMom(0)) > 1e-6)
  {
    std::cout << "comparison of volume prediction\n";
    std::cout << std::setprecision(15) << volGlobal << "\t" << volMom(0) << "\n";
    // dserror("volume prediction is wrong");
    throw std::runtime_error("DirectDiv and Momfitting are not producing same volumes!");
    // std::cout << "DirectDiv and Momfitting are not producing same volumes!" << std::endl;
  }
#endif
#endif

  volcell_->SetVolume(volGlobal);
  //  volcell_->SetVolume(TotalInteg);
  if (std::isnan(volGlobal))
  {
    std::cout << "-------------------------------------------------------------\n";
    std::cout << "There are two possible sources of this problem \n";
    std::cout
        << "1. divCells created from facet may fall on a line. Print the main Gauss points from "
           "GEO::CUT::FacetIntegration::DivergenceIntegrationRule(),"
           " if this is the case, all points belong to a particular divCells have NaN weights\n";
    std::cout << "2. GLOBAL::: The reference plane is not correctly chosen. Print the equation of "
                 "reference plane and if the first component "
                 "is close to zero, then the volume is infinity. Check "
                 "GEO::CUT::DirectDivergenceGlobalRefplane::GetReferencePlane() \n";
    throw std::runtime_error("Volume is not a number.");
  }
}

/*--------------------------------------------------------------------------------------------------------------*
         Integrate given polynomials using the gaussian rule generated using directDivergence.
sudhakar 04/12 Can be used for post-processing
*---------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::DirectDivergence::IntegrateSpecificFuntions(const DRT::UTILS::GaussIntegration& gpc)
{
  double TotalInteg = 0.0;

  //#ifdef LOCAL
  DRT::UTILS::GaussIntegration gpv(gpc);
  for (DRT::UTILS::GaussIntegration::iterator iquad = gpv.begin(); iquad != gpv.end(); ++iquad)
  {
    const LINALG::Matrix<3, 1> etaFacet(iquad.Point());
    const double weiFacet = iquad.Weight();
    double xx = etaFacet(0, 0);
    double yy = etaFacet(1, 0);
    double zz = etaFacet(2, 0);
    TotalInteg +=
        (pow(xx, 6) + xx * pow(yy, 4) * zz + xx * xx * yy * yy * zz * zz + pow(zz, 6)) * weiFacet;
  }
  //#else

  //#endif



  std::cout << std::setprecision(20) << "the integral of a predefined function = " << TotalInteg
            << "\n";
}
