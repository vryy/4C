/*----------------------------------------------------------------------*/
/*! \file

\brief Helper methods for cavitation simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_binstrategy/drt_meshfree_multibin.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_gmsh.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

// print out cut situations for analytical void frac computation
//#define GMSHOUT

/*----------------------------------------------------------------------*
 | fluid fraction calculation                              ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateFluidFraction(Teuchos::RCP<const Epetra_Vector> particleradius)
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateFluidFraction");
  Teuchos::RCP<Epetra_FEVector> fluid_fraction =
      Teuchos::rcp(new Epetra_FEVector(*fluiddis_->ElementRowMap()));
  // export element volume to col layout
  Teuchos::RCP<Epetra_Vector> ele_volume_col =
      LINALG::CreateVector(*fluiddis_->ElementColMap(), false);
  LINALG::Export(*ele_volume_, *ele_volume_col);

  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispnp();

  const bool havepbc = BinStrategy()->HavePBC();

  std::set<int> examinedbins;
  int numrownodes = BinStrategy()->BinDiscret()->NodeRowMap()->NumMyElements();
  for (int i = 0; i < numrownodes; ++i)
  {
    DRT::Node* currentparticle = BinStrategy()->BinDiscret()->lRowNode(i);
    DRT::Element** currentbin = currentparticle->Elements();

    const int binId = currentbin[0]->Id();

    if (examinedbins.count(binId) == 1)
    {
      continue;
    }
    else
    {
      examinedbins.insert(binId);
    }

    // search for underlying fluid elements

    // find maximal radius of particles in this bin to get enough neighboring fluid eles
    DRT::Node** particles = currentbin[0]->Nodes();
    double maxradius = -1.0;
    for (int iparticle = 0; iparticle < currentbin[0]->NumNode(); ++iparticle)
    {
      DRT::Node* currparticle = particles[iparticle];
      const double r_p =
          (*particleradius)[BinStrategy()->BinDiscret()->NodeRowMap()->LID(currparticle->Id())];
      maxradius = std::max(r_p, maxradius);
    }
    if (maxradius < 0.0) dserror("maximum radius smaller than zero");

    // get an ijk-range that is large enough
    int ijk[3];
    BinStrategy()->ConvertGidToijk(binId, ijk);

    // scaling factor in order to account for influence of bubble
    int ibinrange[3];
    for (int dim = 0; dim < 3; ++dim)
    {
      ibinrange[dim] = (int)((maxradius * influencescaling_) / BinStrategy()->BinSize()[dim]) + 1;
      if (ibinrange[dim] > 1)
        dserror("not yet tested for such large bubbles -> think about extending ghosting!");
    }

    int ijk_range[] = {ijk[0] - ibinrange[0], ijk[0] + ibinrange[0], ijk[1] - ibinrange[1],
        ijk[1] + ibinrange[1], ijk[2] - ibinrange[2], ijk[2] + ibinrange[2]};

    // variable to store bin ids of surrounding bins
    std::vector<int> binIds;
    binIds.reserve((2 * ibinrange[0] + 1) * (2 * ibinrange[1] + 1) * (2 * ibinrange[2] + 1));

    // get corresponding bin ids in ijk range and fill them into binIds (in gid)
    BinStrategy()->GidsInijkRange(&ijk_range[0], binIds, false);

    // variable to store all fluid elements in neighborhood
    std::set<DRT::Element*> neighboringfluideles;
    for (std::vector<int>::const_iterator i = binIds.begin(); i != binIds.end(); ++i)
    {
      // extract bins from discretization after checking on existence
      const int lid = BinStrategy()->BinDiscret()->ElementColMap()->LID(*i);
      if (lid < 0) continue;
      // extract bins from discretization
      DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(
          BinStrategy()->BinDiscret()->lColElement(lid));

      DRT::Element** currfluideles = currbin->AssociatedEles(bin_fluidcontent_);
      for (int ifluidele = 0; ifluidele < currbin->NumAssociatedEle(bin_fluidcontent_); ++ifluidele)
      {
        neighboringfluideles.insert(currfluideles[ifluidele]);
      }
    }

    // loop over all row particles in this bin and evaluate them
    for (int iparticle = 0; iparticle < currentbin[0]->NumNode(); ++iparticle)
    {
      DRT::Node* currparticle = particles[iparticle];

      // fill particle position
      LINALG::Matrix<3, 1> particleposition(false);
      std::vector<int> lm_b = BinStrategy()->BinDiscret()->Dof(currparticle);
      // convert dof of particle into correct local id in this Epetra_Vector
      int posx = bubblepos->Map().LID(lm_b[0]);
      for (int dim = 0; dim < 3; ++dim)
      {
        // fill position of currentparticle
        particleposition(dim) = (*bubblepos)[posx + dim];
      }

      // get particle radius and influence of bubble
      const int lid = particleradius->Map().LID(currparticle->Id());
      if (lid < 0) dserror("lid to gid %d not on this proc", currparticle->Id());
      const double r_p = (*particleradius)[lid];
      double influence = influencescaling_ * r_p;

      // bubble volume which is assigned to fluid elements
      const double bubblevol = 4.0 / 3.0 * M_PI * std::pow(r_p, 3);

      // variables for periodic boundary condition corrections
      std::vector<bool> pbcdetected;
      std::vector<LINALG::Matrix<3, 1>> pbceleoffset;

      // store all elements in which the bubble might be located fully inside
      std::vector<int> insideeles;

      const int elelid = ele_volume_col->Map().LID((*neighboringfluideles.begin())->Id());
      if (elelid < 0)
        dserror("element id %i is not on this proc", (*neighboringfluideles.begin())->Id());
      const double elevolume = (*ele_volume_col)[elelid];

      // check for very small bubbles whose volume can be directly assigned to a
      // single fluid element
      {
        const double charactelelength = std::pow(elevolume, 1.0 / 3.0);
        // heuristic criterion of 1/20th of the characteristic element length of one fluid element
        // in the neighborhood
        if (influence < 0.05 * charactelelength)
        {
          // loop all surrounding fluid elements and find possible candidates for assigning
          for (std::set<DRT::Element*>::const_iterator ineighbor = neighboringfluideles.begin();
               ineighbor != neighboringfluideles.end(); ++ineighbor)
          {
            DRT::Element* ele = *ineighbor;
            bool pbcdetected_ele;
            static LINALG::Matrix<3, 1> pbceleoffset_ele;
            const bool xaabboverlap = XAABBoverlap(
                ele, influence, particleposition, havepbc, pbcdetected_ele, pbceleoffset_ele);
            if (xaabboverlap == true)
            {
              insideeles.push_back(ele->Id());
              pbcdetected.push_back(pbcdetected_ele);
              pbceleoffset.push_back(pbceleoffset_ele);
            }
          }
          // assign to the correct fluid element
          AssignSmallBubbles(
              bubblevol, particleposition, insideeles, fluid_fraction, pbcdetected, pbceleoffset);

          continue;  // with next particle in this bin
        }
      }

      // process larger bubbles whose volume is distributed among the surrounding fluid elements
      {
        // variable to store not yet normalized volume for each fluid element
        std::map<int, double> volumefraction;
        double vol_influence = 0.0;
        bool tightcut = false;

        // do while volume is not negative which is important for overlapping surfaces/points of
        // bubble and fluid
        do
        {
          // reset variables if necessary
          if (tightcut == true)
          {
            insideeles.clear();
            volumefraction.clear();
            vol_influence = 0.0;
            tightcut = false;
          }

          // loop over surrounding fluid elements
          for (std::set<DRT::Element*>::const_iterator ineighbor = neighboringfluideles.begin();
               ineighbor != neighboringfluideles.end(); ++ineighbor)
          {
            DRT::Element* ele = *ineighbor;
            bool pbcdetected_ele;
            static LINALG::Matrix<3, 1> pbceleoffset_ele;
            const bool xaabboverlap = XAABBoverlap(
                ele, influence, particleposition, havepbc, pbcdetected_ele, pbceleoffset_ele);

            double vol_ele = 0.0;
            if (xaabboverlap == true)
            {
              // store element id and pbc information
              insideeles.push_back(ele->Id());
              pbcdetected.push_back(pbcdetected_ele);
              pbceleoffset.push_back(pbceleoffset_ele);

              switch (void_frac_strategy_)
              {
                case INPAR::CAVITATION::analytical_constpoly:
                case INPAR::CAVITATION::analytical_quadraticpoly:
                case INPAR::CAVITATION::analytical_quarticpoly:
                {
                  DoAnalyticalIntegrationFluidFrac(ele, currparticle->Id(), particleposition,
                      influence, vol_ele, tightcut, pbcdetected_ele, pbceleoffset_ele);

                  // suppress close to -0.0 void values in case of tight cut situations
                  if (vol_ele < 0.0)
                  {
                    if (vol_ele < -elevolume * 1.0e-8 and tightcut == false)
                    {
                      dserror(
                          "volume of influence is negative: %e (ele id: %d, ele volume: %e, bubble "
                          "id: %d)"
                          "-> check cut procedure",
                          vol_ele, elevolume, ele->Id(), currparticle->Id());
                    }
                    vol_ele = 0.0;
                  }
                }
                break;
                case INPAR::CAVITATION::gaussian_integration:
                {
                  DoGaussianIntegrationFluidFrac(
                      ele, particleposition, influence, vol_ele, pbcdetected_ele, pbceleoffset_ele);
                }
                break;
                default:
                  dserror("void fraction calculation strategy does not exist");
                  break;
              }

              // in case of tight cuts, rerun neighboring ele loop with slightly smaller influence
              if (tightcut == true)
              {
                influence *= 0.999;
                break;
              }

              // safety check
              if (vol_ele < 0.0)
                dserror("negative volume occured during void fraction computation");

              // sum and store volume for each fluid element
              vol_influence += vol_ele;
              volumefraction[ele->Id()] = vol_ele;
            }  // end if xaabboverlap
          }    // end loop neighboring eles

        } while (tightcut == true);

        // distribute volumes stored in volumefraction
        if (vol_influence != 0.0)
        {
          // normalize with volume of bubble divided by volume of influence
          const double normalization = bubblevol / vol_influence;
          // assemble void volume of fluid elements
          for (std::map<int, double>::const_iterator iter = volumefraction.begin();
               iter != volumefraction.end(); ++iter)
          {
            const double val = iter->second * normalization;

            // do assembly of void fraction into fluid
            int err = fluid_fraction->SumIntoGlobalValues(1, &(iter->first), &val);
            if (err < 0) dserror("summing into Epetra_FEVector failed");
          }
        }
        else
        {
          if (void_frac_strategy_ != INPAR::CAVITATION::gaussian_integration)
            dserror("analytical integration used but integrated volume is zero");

          AssignSmallBubbles(
              bubblevol, particleposition, insideeles, fluid_fraction, pbcdetected, pbceleoffset);
        }
      }

    }  // loop over particles in one bin
  }    // loop over outer particles

  // call global assemble
  int err = fluid_fraction->GlobalAssemble(Add, false);
  if (err < 0) dserror("global assemble into fluidforces failed");

  // compute fluid fraction in elements (= ele_volume - void fraction)
  fluid_fraction->Update(1.0, *ele_volume_, -1.0);

  // divide element wise fluid volume by element volume
  fluid_fraction->ReciprocalMultiply(1.0, *ele_volume_, *fluid_fraction, 0.0);

  double minfluidfrac = 0.0;
  fluid_fraction->MinValue(&minfluidfrac);
  if (minfluidfrac < 0.05 && fluid_->PhysicalType() == INPAR::FLUID::loma)
    dserror("minimum fluid fraction is low (%f), increase mesh size or reduce particle size",
        minfluidfrac);

  double maxfluidfrac = 0.0;
  fluid_fraction->MaxValue(&maxfluidfrac);
  if (MyRank() == 0)
    std::cout << std::setprecision(16) << "minimum fluid fraction is: " << minfluidfrac
              << " and maximimum fluid fraction is: " << maxfluidfrac << std::endl;

  // apply fluid fraction to fluid on element level for visualization purpose
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_)->SetFluidFraction(fluid_fraction);

  // transform fluid fraction from element level to node level
  // nodal values are stored in pressure dof of vector as expected in time integration
  switch (fluidfrac_reconstr_)
  {
    case INPAR::CAVITATION::fluidfracreco_spr:
      ComputePatchRecoveredFluidFraction(fluid_fraction);
      break;
    case INPAR::CAVITATION::fluidfracreco_l2:
      ComputeL2ProjectedFluidFraction(fluid_fraction);
      break;
    default:
      dserror("desired fluid fraction projection method not available");
      break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | Gaussian integration for void volume calculation        ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::DoGaussianIntegrationFluidFrac(DRT::Element* ele,
    LINALG::Matrix<3, 1>& particleposition, const double influence, double& vol_ele,
    const bool pbcdetected, const LINALG::Matrix<3, 1>& pbceleoffset)
{
  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // get element location vector and ownerships
  std::vector<int> lm_f;
  std::vector<int> lmowner_f;
  std::vector<int> lmstride;
  ele->LocationVector(*fluiddis_, lm_f, lmowner_f, lmstride);

  // Reshape element matrices and vectors and initialize to zero
  elevector1.Size(1);

  // set action in order to calculate fluid fraction in this element
  Teuchos::ParameterList params;
  params.set<int>("action", FLD::calc_volume_gaussint);
  params.set<double>("influence", influence);
  params.set<int>("gp_per_dir", gauss_rule_per_dir_);

  // in case of pbcs: it is faster to modify the particle position than the element coordinates
  LINALG::Matrix<3, 1> modifiedparticlepos(particleposition);
  if (pbcdetected) modifiedparticlepos.Update(-1.0, pbceleoffset, 1.0);
  params.set<LINALG::Matrix<3, 1>>("particlepos", modifiedparticlepos);

  // call the element specific evaluate method (elevec1 = void volume)
  ele->Evaluate(
      params, *fluiddis_, lm_f, elematrix1, elematrix2, elevector1, elevector2, elevector3);

  vol_ele = elevector1[0];

  return;
}


/*----------------------------------------------------------------------*
 | analytic integration for void volume calculation        ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::DoAnalyticalIntegrationFluidFrac(DRT::Element* ele, const int bubbleid,
    const LINALG::Matrix<3, 1>& particleposition, const double influence, double& vol_ele,
    bool& tightcut, const bool pbcdetected, const LINALG::Matrix<3, 1>& pbceleoffset)
{
  // pull out each volume element and ask for surfaces
  std::vector<Teuchos::RCP<DRT::Element>> surfaces = ele->Surfaces();
  const int numsurfacenodes = surfaces[0]->NumNode();

  // center of fluid element
  static LINALG::Matrix<3, 1> centerele;
  centerele.PutScalar(0.0);
  const DRT::Node* nodeele;
  for (int inodeele = 0; inodeele < ele->NumNode(); ++inodeele)
  {
    nodeele = ele->Nodes()[inodeele];
    for (int dim = 0; dim < 3; ++dim)
    {
      centerele(dim) = centerele(dim) + nodeele->X()[dim];
    }
  }
  centerele.Scale(1.0 / ele->NumNode());

  // in case of pbcs: it is faster to modify the particle position than the element coordinates
  LINALG::Matrix<3, 1> modifiedparticlepos(particleposition);
  if (pbcdetected) modifiedparticlepos.Update(-1.0, pbceleoffset, 1.0);

#ifdef GMSHOUT
  PrintBubbleAndFluidEleToGMSH(ele, bubbleid, modifiedparticlepos, influence);
#endif

  std::map<int, LINALG::Matrix<3, 1>> currentpositions;
  std::vector<LINALG::Matrix<3, 1>> normals(ele->NumSurface());
  // unwarp element surfaces if necessary
  UnwarpElement(ele, surfaces, numsurfacenodes, currentpositions, normals);

  // variables for calculating bubblesurfaces
  std::vector<LINALG::Matrix<3, 1>> bubbleX, bubble_X;
  bubbleX.reserve(numsurfacenodes + 2);
  bubble_X.reserve(numsurfacenodes + 2);

  int numsurfevaluated = 0;

  // loop over surfaces of current fluid element and integrate over surfaces
  for (int isurface = 0; isurface < ele->NumSurface(); ++isurface)
  {
    Teuchos::RCP<DRT::Element> currsurf = surfaces[isurface];

    std::vector<LINALG::Matrix<3, 1>> surfacenodes(numsurfacenodes);

    // loop over nodes of current surface
    for (int inode = 0; inode < numsurfacenodes; ++inode)
    {
      const int id = currsurf->NodeIds()[inode];
      surfacenodes[inode].Update(currentpositions[id]);
    }

    // normal of auxiliary plane
    const LINALG::Matrix<3, 1> n(normals[isurface], true);

    // fill corner points of surface (already unwarped)
    std::map<int, LINALG::Matrix<3, 1>> currentpositionssurf;
    for (int inode = 0; inode < numsurfacenodes; ++inode)
    {
      currentpositionssurf[currsurf->NodeIds()[inode]] = surfacenodes[inode];
    }

    // get XAABB of this single surface element
    const LINALG::Matrix<3, 2> xaabb = GEO::getXAABBofPositions(currentpositionssurf);
    bool boundingbox = true;

    // bubble surfaces of the cubic bubble
    const double bubblesurface[] = {modifiedparticlepos(0) + influence,
        modifiedparticlepos(0) - influence, modifiedparticlepos(1) + influence,
        modifiedparticlepos(1) - influence, modifiedparticlepos(2) + influence,
        modifiedparticlepos(2) - influence};

    // test the bounding box touching the bubble influence
    for (int dim = 0; dim < 3; ++dim)
      if (xaabb(dim, 0) > bubblesurface[dim * 2] or xaabb(dim, 1) < bubblesurface[dim * 2 + 1])
      {
        boundingbox = false;
        break;
      }

    if (boundingbox == true)
    {
      // is there divergence in x-direction?
      if (std::abs(n(0)) > 1.0e-5)
      {
        numsurfevaluated += EvaluateSurface(surfacenodes, n, centerele, modifiedparticlepos,
            influence, vol_ele, tightcut, ele->Id(), isurface);
      }
    }

    // store penetration points of +x and -x surfaces of the bubble
    GetPenetrationPointsOfXSurfaces(
        n, surfacenodes, influence, modifiedparticlepos, bubbleX, bubble_X, tightcut, ele->Id());
  }

  // integration over x-oriented bubble surface
  if (bubbleX.size() != 0)
  {
    static LINALG::Matrix<3, 1> n;
    n(0) = 1.0;
    n(1) = 0.0;
    n(2) = 0.0;

    if (bubbleX.size() > 3)
    {
      BuildConvexHull(bubbleX);
    }

    numsurfevaluated += EvaluateSurface(bubbleX, n, modifiedparticlepos, modifiedparticlepos,
        influence, vol_ele, tightcut, ele->Id(), ele->NumSurface());
  }

  // integration over -x-oriented bubble surface
  if (bubble_X.size() != 0)
  {
    static LINALG::Matrix<3, 1> n;
    n(0) = -1.0;
    n(1) = 0.0;
    n(2) = 0.0;

    if (bubble_X.size() > 3)
    {
      BuildConvexHull(bubble_X);
    }
    numsurfevaluated += EvaluateSurface(bubble_X, n, modifiedparticlepos, modifiedparticlepos,
        influence, vol_ele, tightcut, ele->Id(), ele->NumSurface() + 1);
  }

  // the evaluation of a single surface indicates a tight cut situation
  if (numsurfevaluated == 1) tightcut = true;

#ifdef GMSHOUT
  // finish file here
  gmshfilecontent_.close();
#endif

  return;
}


/*----------------------------------------------------------------------*
 | compute integration points for analytic integration     ghamm 08/13  |
 *----------------------------------------------------------------------*/
int CAVITATION::Algorithm::EvaluateSurface(std::vector<LINALG::Matrix<3, 1>>& surfacenodes,
    const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& centerele,
    const LINALG::Matrix<3, 1>& particleposition, const double influence, double& vol_ele,
    bool& tightcut, const int eleid, const int isurface)
{
  //**********************************************************************
  // Evaluating surfaces over the bubble
  // - calculating corner points
  // - calculating penetration points
  // - calculating bubblecorner points
  //**********************************************************************

  // variable to store integration points
  std::vector<LINALG::Matrix<3, 1>> integrationpoints;
  integrationpoints.reserve(12);
  double bubblesurface[] = {particleposition(0) + influence, particleposition(0) - influence,
      particleposition(1) + influence, particleposition(1) - influence,
      particleposition(2) + influence, particleposition(2) - influence};

  double t;
  const double tol = 1.0e-8;

  const int numsurfacenodes = (int)surfacenodes.size();
  static LINALG::Matrix<3, 1> centersurf;
  centersurf.PutScalar(0.0);
  for (int inode = 0; inode < numsurfacenodes; ++inode)
  {
    centersurf.Update(1.0, surfacenodes[inode], 1.0);
  }
  centersurf.Scale(1.0 / numsurfacenodes);

  // n is wrong oriented
  if ((centersurf(0) - centerele(0)) * n(0) + (centersurf(1) - centerele(1)) * n(1) +
          (centersurf(2) - centerele(2)) * n(2) <
      0)
    dserror("you should not show up here");

  static LINALG::Matrix<3, 1> currpoint, difference;
  // loop over lines to find corner points within the influence of the bubble and penetration points
  // which penetrate the bubble surface
  for (int iedges = 0; iedges < numsurfacenodes; ++iedges)
  {
    difference.Update(1.0, surfacenodes[iedges], -1.0, particleposition);

    // catch corner points cp by testing whether NormInf of the current corner point is smaller than
    // the influence radius corner points are edge-points of the fluid surface, which are into the
    // bubble influence area InfNorm due to a cubic bubble
    if (difference.NormInf() <= influence + tol)
    {
      integrationpoints.push_back(surfacenodes[iedges]);
    }

    // catch penetration points pp which penetrate the bubble surface by building a line between two
    // surface points penetration points are at the edge lines of the surface and penetrate the
    // bubble surface and test whether the penetration happens between them with 0<t<1 and the point
    // is in InfNorm influence of the bubble
    for (int ipene = 0; ipene < 6; ++ipene)
    {
      // parameter t for penetrating bubble surface
      if (std::abs(surfacenodes[(iedges + 1) % numsurfacenodes](ipene / 2) -
                   surfacenodes[iedges](ipene / 2)) > tol)
      {
        t = (bubblesurface[ipene] - surfacenodes[iedges](ipene / 2)) /
            (surfacenodes[(iedges + 1) % numsurfacenodes](ipene / 2) -
                surfacenodes[iedges](ipene / 2));
      }
      else
      {
        t = -1.0;
      }

      // current bubble surface is penetrated (first condition)
      if (+tol < t and t < 1 - tol)
      {
        currpoint.Update(
            t, surfacenodes[(iedges + 1) % numsurfacenodes], 1.0 - t, surfacenodes[iedges]);
        difference.Update(1.0, currpoint, -1.0, particleposition);

        // second condition
        // +tol because penetration points are always in the influence distance in NormInf
        if (difference.NormInf() <= influence + tol)
        {
          integrationpoints.push_back(currpoint);
        }
      }
    }
  }

  // get bubble corner points which are in the surface
  // calculates the bubble corner intersection point of the two surfaces
  // A bubble corner point is an intersection point of the fluid surface with the bubble edges
  // calculates the intersection points of the bubble edges with the plane surface (represented with
  // n,centersurf)
  for (int y = 2; y < 4; y++)
    for (int z = 4; z < 6; z++)
    {
      static LINALG::Matrix<3, 1> bubblecorner;

      CalculateBubbleCornerPoint(n, centersurf, y, z, bubblesurface, bubblecorner);
      bool inpoint = CheckPointInSurface(surfacenodes, centersurf, centerele, bubblecorner);
      difference.Update(1.0, bubblecorner, -1.0, particleposition);
      if (difference.NormInf() <= influence + tol and inpoint == true)
      {
        integrationpoints.push_back(bubblecorner);
      }
    }
  // if normal vector in y-direction == 0: line of bubble is parallel to the fluid surface and can
  // be omitted
  if (std::abs(n(1)) > tol)
    for (int x = 0; x < 2; x++)
      for (int z = 4; z < 6; z++)
      {
        static LINALG::Matrix<3, 1> bubblecorner;

        CalculateBubbleCornerPoint(n, centersurf, x, z, bubblesurface, bubblecorner);
        bool inpoint = CheckPointInSurface(surfacenodes, centersurf, centerele, bubblecorner);
        difference.Update(1.0, bubblecorner, -1.0, particleposition);
        if (difference.NormInf() <= influence + tol and inpoint == true)
        {
          integrationpoints.push_back(bubblecorner);
        }
      }
  // if normal vector in z-direction == 0: line of bubble is parallel to the fluid surface and can
  // be omitted
  if (std::abs(n(2)) > tol)
    for (int x = 0; x < 2; x++)
      for (int y = 2; y < 4; y++)
      {
        static LINALG::Matrix<3, 1> bubblecorner;

        CalculateBubbleCornerPoint(n, centersurf, x, y, bubblesurface, bubblecorner);
        bool inpoint = CheckPointInSurface(surfacenodes, centersurf, centerele, bubblecorner);
        difference.Update(1.0, bubblecorner, -1.0, particleposition);
        if (difference.NormInf() <= influence + tol and inpoint == true)
        {
          integrationpoints.push_back(bubblecorner);
        }
      }

  // calculate the ring integral over fluid surfaces
  if (integrationpoints.size() > 2)
  {
    // integration points are ordered for correct ring integral evaluation
    // close points are removed altering integrationpoints.size()
    {
      BuildConvexHull(integrationpoints);

      if (integrationpoints.size() == 2)
      {
        // nothing to add because only two points left which will end up in a zero contribution
        return 0;
      }
      else if (integrationpoints.size() == 1)
      {
        // surface overlap detected
        tightcut = true;
        return 0;
      }
    }

#ifdef GMSHOUT
    PrintIntegrationLinesToGMSH(integrationpoints, isurface);
#endif

    static LINALG::Matrix<3, 1> centerringintgral;
    centerringintgral.PutScalar(0.0);
    const int numintegrationpoints = (int)integrationpoints.size();
    for (int inode = 0; inode < numintegrationpoints; ++inode)
    {
      centerringintgral.Update(1.0, integrationpoints[inode], 1.0);
    }
    centerringintgral.Scale(1.0 / numintegrationpoints);

    for (int iring = 0; iring < numintegrationpoints; ++iring)
    {
      switch (void_frac_strategy_)
      {
        case INPAR::CAVITATION::analytical_constpoly:
        {
          EvaluateTwoPointsConstPoly(n, centerringintgral, integrationpoints[iring],
              integrationpoints[(iring + 1) % numintegrationpoints], vol_ele);
        }
        break;
        case INPAR::CAVITATION::analytical_quadraticpoly:
        {
          EvaluateTwoPointsQuadraticPoly(n, centerringintgral, integrationpoints[iring],
              integrationpoints[(iring + 1) % numintegrationpoints], particleposition, vol_ele,
              influence);
        }
        break;
        case INPAR::CAVITATION::analytical_quarticpoly:
        {
          EvaluateTwoPointsQuarticPoly(n, centerringintgral, integrationpoints[iring],
              integrationpoints[(iring + 1) % numintegrationpoints], particleposition, vol_ele,
              influence);
        }
        break;
        default:
          dserror("volume strategy does not exist");
          break;
      }
    }
  }
  else if (integrationpoints.size() == 2)
  {
    // nothing to add because only two points left which will end up in a zero contribution
    return 0;
  }
  else if (integrationpoints.size() > 0)
  {
    // surface overlap detected
    tightcut = true;
    return 0;
  }

  return 1;
}


/*----------------------------------------------------------------------*
 | find penetration points of bubble surfaces (+x and -x)  ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::GetPenetrationPointsOfXSurfaces(const LINALG::Matrix<3, 1>& n,
    std::vector<LINALG::Matrix<3, 1>>& surfacenodes, const double influence,
    const LINALG::Matrix<3, 1>& particleposition, std::vector<LINALG::Matrix<3, 1>>& bubbleX,
    std::vector<LINALG::Matrix<3, 1>>& bubble_X, bool& tightcut, const int eleid)
{
  // method counts penetration points of one element and stores them
  // bubbleX array of all penetration points
  // bubble_X array of all penetration points

  double t;
  double tol = 1.0e-8;
  double bubblesurface[] = {particleposition(0) + influence, particleposition(0) - influence};
  const int numsurfacenodes = (int)surfacenodes.size();

  // store penetration points of +x and -x-surfaces
  // penetration points are at the edge lines of the surface and penetrate the bubble surface
  // and test whether the penetration happens between them with 0<t<1 and the point is in NormInf
  // influence of the bubble
  for (int iedges = 0; iedges < numsurfacenodes; ++iedges)
  {
    for (int ipene = 0; ipene < 2; ++ipene)
    {
      // parameter t for penetrating bubble surface
      if (std::abs(surfacenodes[(iedges + 1) % numsurfacenodes](0) - surfacenodes[iedges](0)) > tol)
      {
        t = (bubblesurface[ipene] - surfacenodes[iedges](0)) /
            (surfacenodes[(iedges + 1) % numsurfacenodes](0) - surfacenodes[iedges](0));
        if ((-tol * 10.0 < t and t < tol * 10.0) or (1.0 - tol * 10.0 < t and t < 1.0 + tol * 10.0))
        {
          // due to surface overlap, do integration again with smaller influence radius
          tightcut = true;
          return;
        }
      }
      else
      {
        t = -1.0;
      }

      // current bubble surface is penetrated (first condition)
      if (tol < t and t < 1.0 - tol)
      {
        static LINALG::Matrix<3, 1> difference, pp;
        pp.Update(t, surfacenodes[(iedges + 1) % numsurfacenodes], 1.0 - t, surfacenodes[iedges]);

        // surface which  is only described with n(0)=1 does not penetrate
        if (n(0) < 1.0)
        {
          if (ipene == 0)  // +X bubble surface
          {
            if (bubbleX.size() == 0)
            {
              bubbleX.push_back(pp);
            }
            else
            {
              bool alreadyin = false;
              for (size_t i = 0; i < bubbleX.size(); ++i)
              {
                difference.Update(1.0, pp, -1.0, bubbleX[i]);
                if (difference.Norm2() < tol)
                {
                  alreadyin = true;
                  break;
                }
              }
              if (alreadyin == false)
              {
                bubbleX.push_back(pp);
              }
            }
          }
          else if (ipene == 1)  // -X bubble surface
          {
            if (bubble_X.size() == 0)
            {
              bubble_X.push_back(pp);
            }
            else
            {
              bool alreadyin = false;
              for (size_t i = 0; i < bubble_X.size(); ++i)
              {
                difference.Update(1.0, pp, -1.0, bubble_X[i]);
                if (difference.Norm2() < tol)
                {
                  alreadyin = true;
                  break;
                }
              }
              if (alreadyin == false)
              {
                bubble_X.push_back(pp);
              }
            }
          }
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculation of bubble corner points                     ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateBubbleCornerPoint(const LINALG::Matrix<3, 1>& n,
    const LINALG::Matrix<3, 1>& centersurface, int bubblesurface1, int bubblesurface2,
    double* bubblesurface, LINALG::Matrix<3, 1>& bubblecorner)
{
  const double tol = 1.0e-4;

  // calculates the bubble corner intersection point of the two surfaces
  // A bubble corner point is a intersection point of the fluid surface with the bubble edges
  // surface 0 is +x
  // 1 is -x
  // 2 is y
  // 3 is -y
  // 4 is z
  // 5 is -z
  int iset[3] = {0, 0, 0};
  for (int j = 0; j < 6; ++j)
  {
    if (bubblesurface1 == j)
    {
      bubblecorner(j / 2) = bubblesurface[j];
      iset[j / 2] = 1;
    }
    else if (bubblesurface2 == j)
    {
      bubblecorner(j / 2) = bubblesurface[j];
      iset[j / 2] = 1;
    }
  }
  // x direction variable
  if (iset[0] == 0)
  {
    bubblecorner(0) = (n(0) * centersurface(0) + n(1) * centersurface(1) + n(2) * centersurface(2) -
                          n(1) * bubblecorner(1) - n(2) * bubblecorner(2)) /
                      n(0);
  }
  // y direction variable
  if (iset[1] == 0)
  {
    if (std::abs(n(1)) > tol)
    {
      bubblecorner(1) =
          (n(0) * centersurface(0) + n(1) * centersurface(1) + n(2) * centersurface(2) -
              n(0) * bubblecorner(0) - n(2) * bubblecorner(2)) /
          n(1);
    }
    else
    {
      bubblecorner(1) = centersurface(1);
    }
  }
  // z direction variable
  if (iset[2] == 0)
  {
    if (std::abs(n(2)) > tol)
    {
      bubblecorner(2) =
          (n(0) * centersurface(0) + n(1) * centersurface(1) + n(2) * centersurface(2) -
              n(1) * bubblecorner(1) - n(0) * bubblecorner(0)) /
          n(2);
    }
    else
    {
      bubblecorner(2) = centersurface(2);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | check whether a given point is inside of a surface      ghamm 08/13  |
 *----------------------------------------------------------------------*/
bool CAVITATION::Algorithm::CheckPointInSurface(std::vector<LINALG::Matrix<3, 1>>& surfacenodes,
    const LINALG::Matrix<3, 1>& centersurface, const LINALG::Matrix<3, 1>& centerele,
    const LINALG::Matrix<3, 1>& pointtocheck)
{
  // tests whether a point is in the element surface by calculating the cross product and
  // check whether the point is always left or right of the surrounding lines
  bool inpoint = false;
  int leftright = 0;

  int numsurfacenodes = (int)surfacenodes.size();
  for (int iedges = 0; iedges < numsurfacenodes; ++iedges)
  {
    double DY = surfacenodes[(iedges + 1) % numsurfacenodes](1) - surfacenodes[iedges](1);
    double DZ = surfacenodes[(iedges + 1) % numsurfacenodes](2) - surfacenodes[iedges](2);
    double BY = pointtocheck(1) - surfacenodes[iedges](1);
    double BZ = pointtocheck(2) - surfacenodes[iedges](2);
    if ((DY * BZ - DZ * BY) * (centersurface(0) - centerele(0)) > 0)
    {
      // point is left of current line
    }
    else
    {
      // point is right of current line
      leftright += 1;
    }
  }
  // true is returned if the point is always left or right --> inside
  if (leftright == 0 or leftright == numsurfacenodes) inpoint = true;

  return inpoint;
}


/*----------------------------------------------------------------------*
 | build convex hull of points in y-z plane                ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildConvexHull(std::vector<LINALG::Matrix<3, 1>>& surfacenodes)
{
  int np = (int)surfacenodes.size();
  bool out = false;
  std::vector<MORTAR::Vertex> respoly;
  std::vector<MORTAR::Vertex> collconvexhull;

  // temporary storage for transformed points
  Epetra_SerialDenseMatrix transformed(2, np);

  std::vector<int> dummy;
  std::vector<double> coords(3);
  double maxdist = 0.0;
  // transform each convex hull point
  for (int i = 0; i < np; ++i)
  {
    for (int k = 0; k < 3; ++k) coords[k] = surfacenodes[i](k);

    collconvexhull.push_back(MORTAR::Vertex(
        coords, MORTAR::Vertex::lineclip, dummy, NULL, NULL, false, false, NULL, -1.0));
    // x coordinate doesn't matter. projection into y-z layer
    for (int k = 1; k < 3; ++k)
    {
      transformed(k - 1, i) = coords[k];
    }
    // find largest edge length to adapt tolerance
    static LINALG::Matrix<3, 1> diff;
    for (int neighbor = i + 1; neighbor < np; ++neighbor)
    {
      for (int dim = 0; dim < 3; ++dim)
        diff(dim) = surfacenodes[i](dim) - surfacenodes[neighbor](dim);
    }
    maxdist = std::max(maxdist, sqrt(diff(0) * diff(0) + diff(1) * diff(1) + diff(2) * diff(2)));
  }

  // adapt tolerance to account for size of the problem
  double tol = MORTARCLIPTOL * maxdist;

  // sort convex hull points to obtain final ring integral
  int removed = MORTAR::SortConvexHullPoints(out, transformed, collconvexhull, respoly, tol);

  if (removed > 0) np = np - removed;

  // fill in in new order
  surfacenodes.resize(np);
  for (int i = 0; i < np; ++i)
  {
    for (int k = 0; k < 3; ++k)
    {
      surfacenodes[i](k) = respoly[i].Coord()[k];
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate line integral for analytic integration         ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::EvaluateTwoPointsConstPoly(const LINALG::Matrix<3, 1>& n,
    const LINALG::Matrix<3, 1>& centerringintgral, const LINALG::Matrix<3, 1>& point1,
    const LINALG::Matrix<3, 1>& point2, double& vol_ele)
{
  // calculate the line integral of bubble in the fluid element
  // direction of integration automatically corrected
  static LINALG::Matrix<3, 1> m, delta;
  double Y, dY, Z, dZ, a, b, c;

  delta.Update(1.0, point2, -1.0, point1);
  m(0) = 0.0;
  if (n(0) > 0)
  {
    m(1) = delta(2);
    m(2) = -delta(1);
  }
  else
  {
    m(1) = -delta(2);
    m(2) = delta(1);
  }

  // factors for constant function integration
  a = -0.5 * n(1) / n(0);
  b = centerringintgral(0) + (centerringintgral(1) * n(1) + centerringintgral(2) * n(2)) / n(0);
  c = -n(2) / n(0);

  // declare start of integration
  if ((m(1) * ((point1(1) + point2(1)) / 2 - centerringintgral(1)) +
          m(2) * ((point1(2) + point2(2)) / 2 - centerringintgral(2))) > 0)
  {
    // m is correct oriented
    Y = point1(1);
    Z = point1(2);
    dY = delta(1);
    dZ = delta(2);
  }
  else
  {
    // change integration direction;
    Y = point2(1);
    Z = point2(2);
    dY = -delta(1);
    dZ = -delta(2);
  }

  // divergence in Y direction
  const double linecontribution = dZ * (a * (Y * Y + Y * dY + dY * dY / 3.0) + b * (Y + 0.5 * dY) +
                                           c * (Z * (Y + dY * 0.5) + dZ * (Y * 0.5 + dY / 3.0)));
  vol_ele += linecontribution;

  return;
}


/*----------------------------------------------------------------------*
 | evaluate line integral for analytic integration         ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::EvaluateTwoPointsQuadraticPoly(const LINALG::Matrix<3, 1>& n,
    const LINALG::Matrix<3, 1>& centerringintgral, LINALG::Matrix<3, 1>& point1,
    LINALG::Matrix<3, 1>& point2, const LINALG::Matrix<3, 1>& particleposition, double& vol_ele,
    const double influence)
{
  // calculate the line integral of bubble in the fluid element
  // direction of integration automatically corrected
  static LINALG::Matrix<3, 1> m, delta, point1CoTr, point2CoTr, middleofsurface;
  double a, b, c, d, e, f, g, s;

  // coordinate system transformation
  point1CoTr.Update(1.0, point1, -1.0, particleposition);
  point2CoTr.Update(1.0, point2, -1.0, particleposition);
  delta.Update(1.0, point2CoTr, -1.0, point1CoTr);
  middleofsurface.Update(1.0, centerringintgral, -1.0, particleposition);

  const double invn0 = 1.0 / n(0);
  m(0) = 0;
  m(1) = delta(2) * n(0);
  m(2) = delta(1) * -n(0);

  a = -n(1) * invn0;
  b = -n(2) * invn0;
  c = middleofsurface(0) + middleofsurface(1) * n(1) * invn0 + middleofsurface(2) * n(2) * invn0;

  // check direction of integration
  if ((m(1) * (point1CoTr(1) - middleofsurface(1)) + m(2) * (point1CoTr(2) - middleofsurface(2))) >
      0)
  {
    // m is correct oriented -> direction is correct
    d = point1CoTr(1);
    f = point1CoTr(2);
    e = delta(1);
    g = delta(2);
  }
  else
  {
    // change integration direction
    d = point2CoTr(1);
    f = point2CoTr(2);
    e = -delta(1);
    g = -delta(2);
  }

  // divergence in Y direction
  if (g > 1.0e-8 or g < -1.0e-8)
  {
    // precompute pow values
    const double a2 = a * a;
    const double a3 = a2 * a;

    const double b2 = b * b;
    const double b3 = b2 * b;

    const double c2 = c * c;
    const double c3 = c2 * c;

    const double d2 = d * d;
    const double d3 = d2 * d;
    const double d4 = d3 * d;
    const double d5 = d4 * d;
    const double d6 = d5 * d;

    const double e2 = e * e;
    const double e3 = e2 * e;
    const double e4 = e3 * e;
    const double e5 = e4 * e;
    const double e6 = e5 * e;

    const double f2 = f * f;
    const double f3 = f2 * f;
    const double f4 = f3 * f;
    const double f5 = f4 * f;

    const double g2 = g * g;
    const double g3 = g2 * g;
    const double g4 = g3 * g;
    const double g5 = g4 * g;

    const double n2 = influence * influence;
    const double n4 = n2 * n2;
    const double n6 = n4 * n2;

    // add another line integral
    s = 27.0 / (8.0 * n6);  // = (3/(4n^3))^3 * (2n)^3
    const double linecontribution =
        s * g *
        (((-20 * c3 * d3 * f2 - 45 * a * c2 * d4 * f2 - 36 * a2 * c * d5 * f2 - 10 * a3 * d6 * f2 -
              60 * b * c2 * d3 * f3 - 90 * a * b * c * d4 * f3 - 36 * a2 * b * d5 * f3 -
              60 * b2 * c * d3 * f4 - 45 * a * b2 * d4 * f4 - 20 * b3 * d3 * f5 +
              20 * c3 * d3 * n2 + 45 * a * c2 * d4 * n2 + 36 * a2 * c * d5 * n2 +
              10 * a3 * d6 * n2 + 60 * b * c2 * d3 * f * n2 + 90 * a * b * c * d4 * f * n2 +
              36 * a2 * b * d5 * f * n2 + 60 * c3 * d * f2 * n2 + 90 * a * c2 * d2 * f2 * n2 +
              60 * (1 + a2) * c * d3 * f2 * n2 + 60 * b2 * c * d3 * f2 * n2 +
              15 * a * (3 + a2) * d4 * f2 * n2 + 45 * a * b2 * d4 * f2 * n2 +
              180 * b * c2 * d * f3 * n2 + 180 * a * b * c * d2 * f3 * n2 +
              60 * (1 + a2) * b * d3 * f3 * n2 + 20 * b3 * d3 * f3 * n2 +
              180 * b2 * c * d * f4 * n2 + 90 * a * b2 * d2 * f4 * n2 + 60 * b3 * d * f5 * n2 -
              60 * c3 * d * n4 - 90 * a * c2 * d2 * n4 - 60 * (1 + a2) * c * d3 * n4 -
              15 * a * (3 + a2) * d4 * n4 - 180 * b * c2 * d * f * n4 -
              180 * a * b * c * d2 * f * n4 - 60 * (1 + a2) * b * d3 * f * n4 -
              180 * c * d * f2 * n4 - 180 * b2 * c * d * f2 * n4 - 90 * a * d2 * f2 * n4 -
              90 * a * b2 * d2 * f2 * n4 - 180 * b * d * f3 * n4 - 60 * b3 * d * f3 * n4 +
              180 * c * d * n6 + 90 * a * d2 * n6 + 180 * b * d * f * n6) +
            (-30 * c3 * d2 * e * f2 - 90 * a * c2 * d3 * e * f2 - 90 * a2 * c * d4 * e * f2 -
                30 * a3 * d5 * e * f2 - 90 * b * c2 * d2 * e * f3 - 180 * a * b * c * d3 * e * f3 -
                90 * a2 * b * d4 * e * f3 - 90 * b2 * c * d2 * e * f4 - 90 * a * b2 * d3 * e * f4 -
                30 * b3 * d2 * e * f5 - 20 * c3 * d3 * f * g - 45 * a * c2 * d4 * f * g -
                36 * a2 * c * d5 * f * g - 10 * a3 * d6 * f * g - 90 * b * c2 * d3 * f2 * g -
                135 * a * b * c * d4 * f2 * g - 54 * a2 * b * d5 * f2 * g -
                120 * b2 * c * d3 * f3 * g - 90 * a * b2 * d4 * f3 * g - 50 * b3 * d3 * f4 * g +
                30 * c3 * d2 * e * n2 + 90 * a * c2 * d3 * e * n2 + 90 * a2 * c * d4 * e * n2 +
                30 * a3 * d5 * e * n2 + 90 * b * c2 * d2 * e * f * n2 +
                180 * a * b * c * d3 * e * f * n2 + 90 * a2 * b * d4 * e * f * n2 +
                30 * c3 * e * f2 * n2 + 90 * a * c2 * d * e * f2 * n2 + 90 * c * d2 * e * f2 * n2 +
                90 * a2 * c * d2 * e * f2 * n2 + 90 * b2 * c * d2 * e * f2 * n2 +
                90 * a * d3 * e * f2 * n2 + 30 * a3 * d3 * e * f2 * n2 +
                90 * a * b2 * d3 * e * f2 * n2 + 90 * b * c2 * e * f3 * n2 +
                180 * a * b * c * d * e * f3 * n2 + 90 * b * d2 * e * f3 * n2 +
                90 * a2 * b * d2 * e * f3 * n2 + 30 * b3 * d2 * e * f3 * n2 +
                90 * b2 * c * e * f4 * n2 + 90 * a * b2 * d * e * f4 * n2 + 30 * b3 * e * f5 * n2 +
                30 * b * c2 * d3 * g * n2 + 45 * a * b * c * d4 * g * n2 +
                18 * a2 * b * d5 * g * n2 + 60 * c3 * d * f * g * n2 +
                90 * a * c2 * d2 * f * g * n2 + 60 * c * d3 * f * g * n2 +
                60 * a2 * c * d3 * f * g * n2 + 60 * b2 * c * d3 * f * g * n2 +
                45 * a * d4 * f * g * n2 + 15 * a3 * d4 * f * g * n2 +
                45 * a * b2 * d4 * f * g * n2 + 270 * b * c2 * d * f2 * g * n2 +
                270 * a * b * c * d2 * f2 * g * n2 + 90 * b * d3 * f2 * g * n2 +
                90 * a2 * b * d3 * f2 * g * n2 + 30 * b3 * d3 * f2 * g * n2 +
                360 * b2 * c * d * f3 * g * n2 + 180 * a * b2 * d2 * f3 * g * n2 +
                150 * b3 * d * f4 * g * n2 - 30 * c3 * e * n4 - 90 * a * c2 * d * e * n4 -
                90 * c * d2 * e * n4 - 90 * a2 * c * d2 * e * n4 - 90 * a * d3 * e * n4 -
                30 * a3 * d3 * e * n4 - 90 * b * c2 * e * f * n4 -
                180 * a * b * c * d * e * f * n4 - 90 * b * d2 * e * f * n4 -
                90 * a2 * b * d2 * e * f * n4 - 90 * c * e * f2 * n4 - 90 * b2 * c * e * f2 * n4 -
                90 * a * d * e * f2 * n4 - 90 * a * b2 * d * e * f2 * n4 - 90 * b * e * f3 * n4 -
                30 * b3 * e * f3 * n4 - 90 * b * c2 * d * g * n4 - 90 * a * b * c * d2 * g * n4 -
                30 * b * d3 * g * n4 - 30 * a2 * b * d3 * g * n4 - 180 * c * d * f * g * n4 -
                180 * b2 * c * d * f * g * n4 - 90 * a * d2 * f * g * n4 -
                90 * a * b2 * d2 * f * g * n4 - 270 * b * d * f2 * g * n4 -
                90 * b3 * d * f2 * g * n4 + 90 * c * e * n6 + 90 * a * d * e * n6 +
                90 * b * e * f * n6 + 90 * b * d * g * n6) +
            ((-60 * c3 * d * e2 * f2 - 270 * a * c2 * d2 * e2 * f2 - 360 * a2 * c * d3 * e2 * f2 -
                150 * a3 * d4 * e2 * f2 - 180 * b * c2 * d * e2 * f3 -
                540 * a * b * c * d2 * e2 * f3 - 360 * a2 * b * d3 * e2 * f3 -
                180 * b2 * c * d * e2 * f4 - 270 * a * b2 * d2 * e2 * f4 - 60 * b3 * d * e2 * f5 -
                120 * c3 * d2 * e * f * g - 360 * a * c2 * d3 * e * f * g -
                360 * a2 * c * d4 * e * f * g - 120 * a3 * d5 * e * f * g -
                540 * b * c2 * d2 * e * f2 * g - 1080 * a * b * c * d3 * e * f2 * g -
                540 * a2 * b * d4 * e * f2 * g - 720 * b2 * c * d2 * e * f3 * g -
                720 * a * b2 * d3 * e * f3 * g - 300 * b3 * d2 * e * f4 * g - 20 * c3 * d3 * g2 -
                45 * a * c2 * d4 * g2 - 36 * a2 * c * d5 * g2 - 10 * a3 * d6 * g2 -
                180 * b * c2 * d3 * f * g2 - 270 * a * b * c * d4 * f * g2 -
                108 * a2 * b * d5 * f * g2 - 360 * b2 * c * d3 * f2 * g2 -
                270 * a * b2 * d4 * f2 * g2 - 200 * b3 * d3 * f3 * g2 + 60 * c3 * d * e2 * n2 +
                270 * a * c2 * d2 * e2 * n2 + 360 * a2 * c * d3 * e2 * n2 +
                150 * a3 * d4 * e2 * n2 + 180 * b * c2 * d * e2 * f * n2 +
                540 * a * b * c * d2 * e2 * f * n2 + 360 * a2 * b * d3 * e2 * f * n2 +
                90 * a * c2 * e2 * f2 * n2 + 180 * c * d * e2 * f2 * n2 +
                180 * a2 * c * d * e2 * f2 * n2 + 180 * b2 * c * d * e2 * f2 * n2 +
                270 * a * d2 * e2 * f2 * n2 + 90 * a3 * d2 * e2 * f2 * n2 +
                270 * a * b2 * d2 * e2 * f2 * n2 + 180 * a * b * c * e2 * f3 * n2 +
                180 * b * d * e2 * f3 * n2 + 180 * a2 * b * d * e2 * f3 * n2 +
                60 * b3 * d * e2 * f3 * n2 + 90 * a * b2 * e2 * f4 * n2 +
                180 * b * c2 * d2 * e * g * n2 + 360 * a * b * c * d3 * e * g * n2 +
                180 * a2 * b * d4 * e * g * n2 + 120 * c3 * e * f * g * n2 +
                360 * a * c2 * d * e * f * g * n2 + 360 * c * d2 * e * f * g * n2 +
                360 * a2 * c * d2 * e * f * g * n2 + 360 * b2 * c * d2 * e * f * g * n2 +
                360 * a * d3 * e * f * g * n2 + 120 * a3 * d3 * e * f * g * n2 +
                360 * a * b2 * d3 * e * f * g * n2 + 540 * b * c2 * e * f2 * g * n2 +
                1080 * a * b * c * d * e * f2 * g * n2 + 540 * b * d2 * e * f2 * g * n2 +
                540 * a2 * b * d2 * e * f2 * g * n2 + 180 * b3 * d2 * e * f2 * g * n2 +
                720 * b2 * c * e * f3 * g * n2 + 720 * a * b2 * d * e * f3 * g * n2 +
                300 * b3 * e * f4 * g * n2 + 60 * c3 * d * g2 * n2 + 90 * a * c2 * d2 * g2 * n2 +
                60 * c * d3 * g2 * n2 + 60 * a2 * c * d3 * g2 * n2 + 60 * b2 * c * d3 * g2 * n2 +
                45 * a * d4 * g2 * n2 + 15 * a3 * d4 * g2 * n2 + 45 * a * b2 * d4 * g2 * n2 +
                540 * b * c2 * d * f * g2 * n2 + 540 * a * b * c * d2 * f * g2 * n2 +
                180 * b * d3 * f * g2 * n2 + 180 * a2 * b * d3 * f * g2 * n2 +
                60 * b3 * d3 * f * g2 * n2 + 1080 * b2 * c * d * f2 * g2 * n2 +
                540 * a * b2 * d2 * f2 * g2 * n2 + 600 * b3 * d * f3 * g2 * n2 -
                90 * a * c2 * e2 * n4 - 180 * c * d * e2 * n4 - 180 * a2 * c * d * e2 * n4 -
                270 * a * d2 * e2 * n4 - 90 * a3 * d2 * e2 * n4 - 180 * a * b * c * e2 * f * n4 -
                180 * b * d * e2 * f * n4 - 180 * a2 * b * d * e2 * f * n4 - 90 * a * e2 * f2 * n4 -
                90 * a * b2 * e2 * f2 * n4 - 180 * b * c2 * e * g * n4 -
                360 * a * b * c * d * e * g * n4 - 180 * b * d2 * e * g * n4 -
                180 * a2 * b * d2 * e * g * n4 - 360 * c * e * f * g * n4 -
                360 * b2 * c * e * f * g * n4 - 360 * a * d * e * f * g * n4 -
                360 * a * b2 * d * e * f * g * n4 - 540 * b * e * f2 * g * n4 -
                180 * b3 * e * f2 * g * n4 - 180 * c * d * g2 * n4 - 180 * b2 * c * d * g2 * n4 -
                90 * a * d2 * g2 * n4 - 90 * a * b2 * d2 * g2 * n4 - 540 * b * d * f * g2 * n4 -
                180 * b3 * d * f * g2 * n4 + 90 * a * e2 * n6 + 180 * b * e * g * n6)) /
                3 +
            ((-10 * c3 * e3 * f2 - 90 * a * c2 * d * e3 * f2 - 180 * a2 * c * d2 * e3 * f2 -
                100 * a3 * d3 * e3 * f2 - 30 * b * c2 * e3 * f3 - 180 * a * b * c * d * e3 * f3 -
                180 * a2 * b * d2 * e3 * f3 - 30 * b2 * c * e3 * f4 - 90 * a * b2 * d * e3 * f4 -
                10 * b3 * e3 * f5 - 60 * c3 * d * e2 * f * g - 270 * a * c2 * d2 * e2 * f * g -
                360 * a2 * c * d3 * e2 * f * g - 150 * a3 * d4 * e2 * f * g -
                270 * b * c2 * d * e2 * f2 * g - 810 * a * b * c * d2 * e2 * f2 * g -
                540 * a2 * b * d3 * e2 * f2 * g - 360 * b2 * c * d * e2 * f3 * g -
                540 * a * b2 * d2 * e2 * f3 * g - 150 * b3 * d * e2 * f4 * g -
                30 * c3 * d2 * e * g2 - 90 * a * c2 * d3 * e * g2 - 90 * a2 * c * d4 * e * g2 -
                30 * a3 * d5 * e * g2 - 270 * b * c2 * d2 * e * f * g2 -
                540 * a * b * c * d3 * e * f * g2 - 270 * a2 * b * d4 * e * f * g2 -
                540 * b2 * c * d2 * e * f2 * g2 - 540 * a * b2 * d3 * e * f2 * g2 -
                300 * b3 * d2 * e * f3 * g2 - 30 * b * c2 * d3 * g3 - 45 * a * b * c * d4 * g3 -
                18 * a2 * b * d5 * g3 - 120 * b2 * c * d3 * f * g3 - 90 * a * b2 * d4 * f * g3 -
                100 * b3 * d3 * f2 * g3 + 10 * c3 * e3 * n2 + 90 * a * c2 * d * e3 * n2 +
                180 * a2 * c * d2 * e3 * n2 + 100 * a3 * d3 * e3 * n2 + 30 * b * c2 * e3 * f * n2 +
                180 * a * b * c * d * e3 * f * n2 + 180 * a2 * b * d2 * e3 * f * n2 +
                30 * c * e3 * f2 * n2 + 30 * a2 * c * e3 * f2 * n2 + 30 * b2 * c * e3 * f2 * n2 +
                90 * a * d * e3 * f2 * n2 + 30 * a3 * d * e3 * f2 * n2 +
                90 * a * b2 * d * e3 * f2 * n2 + 30 * b * e3 * f3 * n2 +
                30 * a2 * b * e3 * f3 * n2 + 10 * b3 * e3 * f3 * n2 +
                90 * b * c2 * d * e2 * g * n2 + 270 * a * b * c * d2 * e2 * g * n2 +
                180 * a2 * b * d3 * e2 * g * n2 + 90 * a * c2 * e2 * f * g * n2 +
                180 * c * d * e2 * f * g * n2 + 180 * a2 * c * d * e2 * f * g * n2 +
                180 * b2 * c * d * e2 * f * g * n2 + 270 * a * d2 * e2 * f * g * n2 +
                90 * a3 * d2 * e2 * f * g * n2 + 270 * a * b2 * d2 * e2 * f * g * n2 +
                270 * a * b * c * e2 * f2 * g * n2 + 270 * b * d * e2 * f2 * g * n2 +
                270 * a2 * b * d * e2 * f2 * g * n2 + 90 * b3 * d * e2 * f2 * g * n2 +
                180 * a * b2 * e2 * f3 * g * n2 + 30 * c3 * e * g2 * n2 +
                90 * a * c2 * d * e * g2 * n2 + 90 * c * d2 * e * g2 * n2 +
                90 * a2 * c * d2 * e * g2 * n2 + 90 * b2 * c * d2 * e * g2 * n2 +
                90 * a * d3 * e * g2 * n2 + 30 * a3 * d3 * e * g2 * n2 +
                90 * a * b2 * d3 * e * g2 * n2 + 270 * b * c2 * e * f * g2 * n2 +
                540 * a * b * c * d * e * f * g2 * n2 + 270 * b * d2 * e * f * g2 * n2 +
                270 * a2 * b * d2 * e * f * g2 * n2 + 90 * b3 * d2 * e * f * g2 * n2 +
                540 * b2 * c * e * f2 * g2 * n2 + 540 * a * b2 * d * e * f2 * g2 * n2 +
                300 * b3 * e * f3 * g2 * n2 + 90 * b * c2 * d * g3 * n2 +
                90 * a * b * c * d2 * g3 * n2 + 30 * b * d3 * g3 * n2 + 30 * a2 * b * d3 * g3 * n2 +
                10 * b3 * d3 * g3 * n2 + 360 * b2 * c * d * f * g3 * n2 +
                180 * a * b2 * d2 * f * g3 * n2 + 300 * b3 * d * f2 * g3 * n2 - 30 * c * e3 * n4 -
                30 * a2 * c * e3 * n4 - 90 * a * d * e3 * n4 - 30 * a3 * d * e3 * n4 -
                30 * b * e3 * f * n4 - 30 * a2 * b * e3 * f * n4 - 90 * a * b * c * e2 * g * n4 -
                90 * b * d * e2 * g * n4 - 90 * a2 * b * d * e2 * g * n4 -
                90 * a * e2 * f * g * n4 - 90 * a * b2 * e2 * f * g * n4 - 90 * c * e * g2 * n4 -
                90 * b2 * c * e * g2 * n4 - 90 * a * d * e * g2 * n4 -
                90 * a * b2 * d * e * g2 * n4 - 270 * b * e * f * g2 * n4 -
                90 * b3 * e * f * g2 * n4 - 90 * b * d * g3 * n4 - 30 * b3 * d * g3 * n4)) /
                2 +
            (-9 * a * c2 * e4 * f2 - 36 * a2 * c * d * e4 * f2 - 30 * a3 * d2 * e4 * f2 -
                18 * a * b * c * e4 * f3 - 36 * a2 * b * d * e4 * f3 - 9 * a * b2 * e4 * f4 -
                8 * c3 * e3 * f * g - 72 * a * c2 * d * e3 * f * g -
                144 * a2 * c * d2 * e3 * f * g - 80 * a3 * d3 * e3 * f * g -
                36 * b * c2 * e3 * f2 * g - 216 * a * b * c * d * e3 * f2 * g -
                216 * a2 * b * d2 * e3 * f2 * g - 48 * b2 * c * e3 * f3 * g -
                144 * a * b2 * d * e3 * f3 * g - 20 * b3 * e3 * f4 * g - 12 * c3 * d * e2 * g2 -
                54 * a * c2 * d2 * e2 * g2 - 72 * a2 * c * d3 * e2 * g2 - 30 * a3 * d4 * e2 * g2 -
                108 * b * c2 * d * e2 * f * g2 - 324 * a * b * c * d2 * e2 * f * g2 -
                216 * a2 * b * d3 * e2 * f * g2 - 216 * b2 * c * d * e2 * f2 * g2 -
                324 * a * b2 * d2 * e2 * f2 * g2 - 120 * b3 * d * e2 * f3 * g2 -
                36 * b * c2 * d2 * e * g3 - 72 * a * b * c * d3 * e * g3 -
                36 * a2 * b * d4 * e * g3 - 144 * b2 * c * d2 * e * f * g3 -
                144 * a * b2 * d3 * e * f * g3 - 120 * b3 * d2 * e * f2 * g3 -
                12 * b2 * c * d3 * g4 - 9 * a * b2 * d4 * g4 - 20 * b3 * d3 * f * g4 +
                9 * a * c2 * e4 * n2 + 36 * a2 * c * d * e4 * n2 + 30 * a3 * d2 * e4 * n2 +
                18 * a * b * c * e4 * f * n2 + 36 * a2 * b * d * e4 * f * n2 +
                9 * a * e4 * f2 * n2 + 3 * a3 * e4 * f2 * n2 + 9 * a * b2 * e4 * f2 * n2 +
                12 * b * c2 * e3 * g * n2 + 72 * a * b * c * d * e3 * g * n2 +
                72 * a2 * b * d2 * e3 * g * n2 + 24 * c * e3 * f * g * n2 +
                24 * a2 * c * e3 * f * g * n2 + 24 * b2 * c * e3 * f * g * n2 +
                72 * a * d * e3 * f * g * n2 + 24 * a3 * d * e3 * f * g * n2 +
                72 * a * b2 * d * e3 * f * g * n2 + 36 * b * e3 * f2 * g * n2 +
                36 * a2 * b * e3 * f2 * g * n2 + 12 * b3 * e3 * f2 * g * n2 +
                18 * a * c2 * e2 * g2 * n2 + 36 * c * d * e2 * g2 * n2 +
                36 * a2 * c * d * e2 * g2 * n2 + 36 * b2 * c * d * e2 * g2 * n2 +
                54 * a * d2 * e2 * g2 * n2 + 18 * a3 * d2 * e2 * g2 * n2 +
                54 * a * b2 * d2 * e2 * g2 * n2 + 108 * a * b * c * e2 * f * g2 * n2 +
                108 * b * d * e2 * f * g2 * n2 + 108 * a2 * b * d * e2 * f * g2 * n2 +
                36 * b3 * d * e2 * f * g2 * n2 + 108 * a * b2 * e2 * f2 * g2 * n2 +
                36 * b * c2 * e * g3 * n2 + 72 * a * b * c * d * e * g3 * n2 +
                36 * b * d2 * e * g3 * n2 + 36 * a2 * b * d2 * e * g3 * n2 +
                12 * b3 * d2 * e * g3 * n2 + 144 * b2 * c * e * f * g3 * n2 +
                144 * a * b2 * d * e * f * g3 * n2 + 120 * b3 * e * f2 * g3 * n2 +
                36 * b2 * c * d * g4 * n2 + 18 * a * b2 * d2 * g4 * n2 + 60 * b3 * d * f * g4 * n2 -
                9 * a * e4 * n4 - 3 * a3 * e4 * n4 - 12 * b * e3 * g * n4 -
                12 * a2 * b * e3 * g * n4 - 18 * a * e2 * g2 * n4 - 18 * a * b2 * e2 * g2 * n4 -
                36 * b * e * g3 * n4 - 12 * b3 * e * g3 * n4) +
            ((-18 * a2 * c * e5 * f2 - 30 * a3 * d * e5 * f2 - 18 * a2 * b * e5 * f3 -
                45 * a * c2 * e4 * f * g - 180 * a2 * c * d * e4 * f * g -
                150 * a3 * d2 * e4 * f * g - 135 * a * b * c * e4 * f2 * g -
                270 * a2 * b * d * e4 * f2 * g - 90 * a * b2 * e4 * f3 * g - 10 * c3 * e3 * g2 -
                90 * a * c2 * d * e3 * g2 - 180 * a2 * c * d2 * e3 * g2 - 100 * a3 * d3 * e3 * g2 -
                90 * b * c2 * e3 * f * g2 - 540 * a * b * c * d * e3 * f * g2 -
                540 * a2 * b * d2 * e3 * f * g2 - 180 * b2 * c * e3 * f2 * g2 -
                540 * a * b2 * d * e3 * f2 * g2 - 100 * b3 * e3 * f3 * g2 -
                90 * b * c2 * d * e2 * g3 - 270 * a * b * c * d2 * e2 * g3 -
                180 * a2 * b * d3 * e2 * g3 - 360 * b2 * c * d * e2 * f * g3 -
                540 * a * b2 * d2 * e2 * f * g3 - 300 * b3 * d * e2 * f2 * g3 -
                90 * b2 * c * d2 * e * g4 - 90 * a * b2 * d3 * e * g4 - 150 * b3 * d2 * e * f * g4 -
                10 * b3 * d3 * g5 + 18 * a2 * c * e5 * n2 + 30 * a3 * d * e5 * n2 +
                18 * a2 * b * e5 * f * n2 + 45 * a * b * c * e4 * g * n2 +
                90 * a2 * b * d * e4 * g * n2 + 45 * a * e4 * f * g * n2 +
                15 * a3 * e4 * f * g * n2 + 45 * a * b2 * e4 * f * g * n2 + 30 * c * e3 * g2 * n2 +
                30 * a2 * c * e3 * g2 * n2 + 30 * b2 * c * e3 * g2 * n2 +
                90 * a * d * e3 * g2 * n2 + 30 * a3 * d * e3 * g2 * n2 +
                90 * a * b2 * d * e3 * g2 * n2 + 90 * b * e3 * f * g2 * n2 +
                90 * a2 * b * e3 * f * g2 * n2 + 30 * b3 * e3 * f * g2 * n2 +
                90 * a * b * c * e2 * g3 * n2 + 90 * b * d * e2 * g3 * n2 +
                90 * a2 * b * d * e2 * g3 * n2 + 30 * b3 * d * e2 * g3 * n2 +
                180 * a * b2 * e2 * f * g3 * n2 + 90 * b2 * c * e * g4 * n2 +
                90 * a * b2 * d * e * g4 * n2 + 150 * b3 * e * f * g4 * n2 +
                30 * b3 * d * g5 * n2)) /
                3 +
            ((-10 * a3 * e6 * f2 - 72 * a2 * c * e5 * f * g - 120 * a3 * d * e5 * f * g -
                108 * a2 * b * e5 * f2 * g - 45 * a * c2 * e4 * g2 - 180 * a2 * c * d * e4 * g2 -
                150 * a3 * d2 * e4 * g2 - 270 * a * b * c * e4 * f * g2 -
                540 * a2 * b * d * e4 * f * g2 - 270 * a * b2 * e4 * f2 * g2 -
                60 * b * c2 * e3 * g3 - 360 * a * b * c * d * e3 * g3 -
                360 * a2 * b * d2 * e3 * g3 - 240 * b2 * c * e3 * f * g3 -
                720 * a * b2 * d * e3 * f * g3 - 200 * b3 * e3 * f2 * g3 -
                180 * b2 * c * d * e2 * g4 - 270 * a * b2 * d2 * e2 * g4 -
                300 * b3 * d * e2 * f * g4 - 60 * b3 * d2 * e * g5 + 10 * a3 * e6 * n2 +
                36 * a2 * b * e5 * g * n2 + 45 * a * e4 * g2 * n2 + 15 * a3 * e4 * g2 * n2 +
                45 * a * b2 * e4 * g2 * n2 + 60 * b * e3 * g3 * n2 + 60 * a2 * b * e3 * g3 * n2 +
                20 * b3 * e3 * g3 * n2 + 90 * a * b2 * e2 * g4 * n2 + 60 * b3 * e * g5 * n2)) /
                7 -
            (e2 * g *
                (10 * a3 * e4 * f + 18 * a2 * c * e3 * g + 30 * a3 * d * e3 * g +
                    54 * a2 * b * e3 * f * g + 45 * a * b * c * e2 * g2 +
                    90 * a2 * b * d * e2 * g2 + 90 * a * b2 * e2 * f * g2 + 30 * b2 * c * e * g3 +
                    90 * a * b2 * d * e * g3 + 50 * b3 * e * f * g3 + 30 * b3 * d * g4)) /
                4 -
            (e3 * g2 *
                (10 * a3 * e3 + 36 * a2 * b * e2 * g + 45 * a * b2 * e * g2 + 20 * b3 * g3)) /
                9)) /
        180;

    vol_ele += linecontribution;
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate line integral for analytic integration         ghamm 06/16  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::EvaluateTwoPointsQuarticPoly(const LINALG::Matrix<3, 1>& n,
    const LINALG::Matrix<3, 1>& centerringintgral, LINALG::Matrix<3, 1>& point1,
    LINALG::Matrix<3, 1>& point2, const LINALG::Matrix<3, 1>& particleposition, double& vol_ele,
    const double influence)
{
  // calculate the line integral of bubble in the fluid element
  // direction of integration automatically corrected
  static LINALG::Matrix<3, 1> m, delta, point1CoTr, point2CoTr, middleofsurface;
  double a, b, c, d, e, f, g;

  // coordinate system transformation
  point1CoTr.Update(1.0, point1, -1.0, particleposition);
  point2CoTr.Update(1.0, point2, -1.0, particleposition);
  delta.Update(1.0, point2CoTr, -1.0, point1CoTr);
  middleofsurface.Update(1.0, centerringintgral, -1.0, particleposition);

  const double invn0 = 1.0 / n(0);
  m(0) = 0;
  m(1) = delta(2) * n(0);
  m(2) = delta(1) * -n(0);

  a = -n(1) * invn0;
  b = -n(2) * invn0;
  c = middleofsurface(0) + middleofsurface(1) * n(1) * invn0 + middleofsurface(2) * n(2) * invn0;

  // check direction of integration
  if ((m(1) * (point1CoTr(1) - middleofsurface(1)) + m(2) * (point1CoTr(2) - middleofsurface(2))) >
      0)
  {
    // m is correct oriented -> direction is correct
    d = point1CoTr(1);
    f = point1CoTr(2);
    e = delta(1);
    g = delta(2);
  }
  else
  {
    // change integration direction
    d = point2CoTr(1);
    f = point2CoTr(2);
    e = -delta(1);
    g = -delta(2);
  }

  // divergence in Y direction
  if (g > 1.0e-8 or g < -1.0e-8)
  {
    // precompute pow values
    const double a2 = a * a;
    const double a3 = a2 * a;
    const double a4 = a3 * a;
    const double a5 = a4 * a;

    const double b2 = b * b;
    const double b3 = b2 * b;
    const double b4 = b3 * b;
    const double b5 = b4 * b;

    const double c2 = c * c;
    const double c3 = c2 * c;
    const double c4 = c3 * c;
    const double c5 = c4 * c;

    const double d2 = d * d;
    const double d3 = d2 * d;
    const double d4 = d3 * d;
    const double d5 = d4 * d;
    const double d6 = d5 * d;
    const double d7 = d6 * d;
    const double d8 = d7 * d;
    const double d9 = d8 * d;
    const double d10 = d9 * d;

    const double e2 = e * e;
    const double e3 = e2 * e;
    const double e4 = e3 * e;
    const double e5 = e4 * e;
    const double e6 = e5 * e;
    const double e7 = e6 * e;
    const double e8 = e7 * e;
    const double e9 = e8 * e;
    const double e10 = e9 * e;

    const double f2 = f * f;
    const double f3 = f2 * f;
    const double f4 = f3 * f;
    const double f5 = f4 * f;
    const double f6 = f5 * f;
    const double f7 = f6 * f;

    const double g2 = g * g;
    const double g3 = g2 * g;
    const double g4 = g3 * g;
    const double g5 = g4 * g;
    const double g6 = g5 * g;

    const double n2 = influence * influence;
    const double n3 = n2 * influence;
    const double n4 = n3 * influence;
    const double n5 = n4 * influence;
    const double n6 = n5 * influence;
    const double n7 = n6 * influence;
    const double n8 = n7 * influence;
    const double n9 = n8 * influence;
    const double n10 = n9 * influence;
    const double n11 = n10 * influence;
    const double n15 = n11 * n4;

    // add another line integral
    const double s = 8.0 * n3 * 15.0 * 15.0 * 15.0 / (16.0 * 16.0 * 16.0);
    const double linecontribution =
        s * g *
        ((g3 *
             ((252 * c5 * d5 * f + 1050 * a * c4 * d6 * f + 1800 * a2 * c3 * d7 * f +
                  1575 * a3 * c2 * d8 * f + 700 * a4 * c * d9 * f + 126 * a5 * d10 * f +
                  1260 * b * c4 * d5 * f2 + 4200 * a * b * c3 * d6 * f2 +
                  5400 * a2 * b * c2 * d7 * f2 + 3150 * a3 * b * c * d8 * f2 +
                  700 * a4 * b * d9 * f2 + 2520 * b2 * c3 * d5 * f3 + 6300 * a * b2 * c2 * d6 * f3 +
                  5400 * a2 * b2 * c * d7 * f3 + 1575 * a3 * b2 * d8 * f3 +
                  2520 * b3 * c2 * d5 * f4 + 4200 * a * b3 * c * d6 * f4 +
                  1800 * a2 * b3 * d7 * f4 + 1260 * b4 * c * d5 * f5 + 1050 * a * b4 * d6 * f5 +
                  252 * b5 * d5 * f6 - 840 * c5 * d3 * f * n2 - 3150 * a * c4 * d4 * f * n2 -
                  840 * c3 * d5 * f * n2 - 5040 * a2 * c3 * d5 * f * n2 -
                  2100 * a * c2 * d6 * f * n2 - 4200 * a3 * c2 * d6 * f * n2 -
                  1800 * a2 * c * d7 * f * n2 - 1800 * a4 * c * d7 * f * n2 -
                  525 * a3 * d8 * f * n2 - 315 * a5 * d8 * f * n2 - 4200 * b * c4 * d3 * f2 * n2 -
                  12600 * a * b * c3 * d4 * f2 * n2 - 2520 * b * c2 * d5 * f2 * n2 -
                  15120 * a2 * b * c2 * d5 * f2 * n2 - 4200 * a * b * c * d6 * f2 * n2 -
                  8400 * a3 * b * c * d6 * f2 * n2 - 1800 * a2 * b * d7 * f2 * n2 -
                  1800 * a4 * b * d7 * f2 * n2 - 8400 * b2 * c3 * d3 * f3 * n2 -
                  18900 * a * b2 * c2 * d4 * f3 * n2 - 2520 * b2 * c * d5 * f3 * n2 -
                  15120 * a2 * b2 * c * d5 * f3 * n2 - 2100 * a * b2 * d6 * f3 * n2 -
                  4200 * a3 * b2 * d6 * f3 * n2 - 8400 * b3 * c2 * d3 * f4 * n2 -
                  12600 * a * b3 * c * d4 * f4 * n2 - 840 * b3 * d5 * f4 * n2 -
                  5040 * a2 * b3 * d5 * f4 * n2 - 4200 * b4 * c * d3 * f5 * n2 -
                  3150 * a * b4 * d4 * f5 * n2 - 840 * b5 * d3 * f6 * n2 + 1260 * c5 * d * f * n4 +
                  3150 * a * c4 * d2 * f * n4 + 2800 * c3 * d3 * f * n4 +
                  4200 * a2 * c3 * d3 * f * n4 + 6300 * a * c2 * d4 * f * n4 +
                  3150 * a3 * c2 * d4 * f * n4 + 1260 * c * d5 * f * n4 +
                  5040 * a2 * c * d5 * f * n4 + 1260 * a4 * c * d5 * f * n4 +
                  1050 * a * d6 * f * n4 + 1400 * a3 * d6 * f * n4 + 210 * a5 * d6 * f * n4 +
                  6300 * b * c4 * d * f2 * n4 + 12600 * a * b * c3 * d2 * f2 * n4 +
                  8400 * b * c2 * d3 * f2 * n4 + 12600 * a2 * b * c2 * d3 * f2 * n4 +
                  12600 * a * b * c * d4 * f2 * n4 + 6300 * a3 * b * c * d4 * f2 * n4 +
                  1260 * b * d5 * f2 * n4 + 5040 * a2 * b * d5 * f2 * n4 +
                  1260 * a4 * b * d5 * f2 * n4 + 12600 * b2 * c3 * d * f3 * n4 +
                  18900 * a * b2 * c2 * d2 * f3 * n4 + 8400 * b2 * c * d3 * f3 * n4 +
                  12600 * a2 * b2 * c * d3 * f3 * n4 + 6300 * a * b2 * d4 * f3 * n4 +
                  3150 * a3 * b2 * d4 * f3 * n4 + 12600 * b3 * c2 * d * f4 * n4 +
                  12600 * a * b3 * c * d2 * f4 * n4 + 2800 * b3 * d3 * f4 * n4 +
                  4200 * a2 * b3 * d3 * f4 * n4 + 6300 * b4 * c * d * f5 * n4 +
                  3150 * a * b4 * d2 * f5 * n4 + 1260 * b5 * d * f6 * n4 - 4200 * c3 * d * f * n6 -
                  6300 * a * c2 * d2 * f * n6 - 4200 * c * d3 * f * n6 -
                  4200 * a2 * c * d3 * f * n6 - 3150 * a * d4 * f * n6 - 1050 * a3 * d4 * f * n6 -
                  12600 * b * c2 * d * f2 * n6 - 12600 * a * b * c * d2 * f2 * n6 -
                  4200 * b * d3 * f2 * n6 - 4200 * a2 * b * d3 * f2 * n6 -
                  12600 * b2 * c * d * f3 * n6 - 6300 * a * b2 * d2 * f3 * n6 -
                  4200 * b3 * d * f4 * n6 + 6300 * c * d * f * n8 + 3150 * a * d2 * f * n8 +
                  6300 * b * d * f2 * n8) +
                 ((5040 * c5 * d4 * e * f + 25200 * a * c4 * d5 * e * f +
                     50400 * a2 * c3 * d6 * e * f + 50400 * a3 * c2 * d7 * e * f +
                     25200 * a4 * c * d8 * e * f + 5040 * a5 * d9 * e * f +
                     25200 * b * c4 * d4 * e * f2 + 100800 * a * b * c3 * d5 * e * f2 +
                     151200 * a2 * b * c2 * d6 * e * f2 + 100800 * a3 * b * c * d7 * e * f2 +
                     25200 * a4 * b * d8 * e * f2 + 50400 * b2 * c3 * d4 * e * f3 +
                     151200 * a * b2 * c2 * d5 * e * f3 + 151200 * a2 * b2 * c * d6 * e * f3 +
                     50400 * a3 * b2 * d7 * e * f3 + 50400 * b3 * c2 * d4 * e * f4 +
                     100800 * a * b3 * c * d5 * e * f4 + 50400 * a2 * b3 * d6 * e * f4 +
                     25200 * b4 * c * d4 * e * f5 + 25200 * a * b4 * d5 * e * f5 +
                     5040 * b5 * d4 * e * f6 + 252 * c5 * d5 * g + 1050 * a * c4 * d6 * g +
                     1800 * a2 * c3 * d7 * g + 1575 * a3 * c2 * d8 * g + 700 * a4 * c * d9 * g +
                     126 * a5 * d10 * g + 6300 * b * c4 * d5 * f * g +
                     21000 * a * b * c3 * d6 * f * g + 27000 * a2 * b * c2 * d7 * f * g +
                     15750 * a3 * b * c * d8 * f * g + 3500 * a4 * b * d9 * f * g +
                     22680 * b2 * c3 * d5 * f2 * g + 56700 * a * b2 * c2 * d6 * f2 * g +
                     48600 * a2 * b2 * c * d7 * f2 * g + 14175 * a3 * b2 * d8 * f2 * g +
                     32760 * b3 * c2 * d5 * f3 * g + 54600 * a * b3 * c * d6 * f3 * g +
                     23400 * a2 * b3 * d7 * f3 * g + 21420 * b4 * c * d5 * f4 * g +
                     17850 * a * b4 * d6 * f4 * g + 5292 * b5 * d5 * f5 * g -
                     10080 * c5 * d2 * e * f * n2 - 50400 * a * c4 * d3 * e * f * n2 -
                     16800 * c3 * d4 * e * f * n2 - 100800 * a2 * c3 * d4 * e * f * n2 -
                     50400 * a * c2 * d5 * e * f * n2 - 100800 * a3 * c2 * d5 * e * f * n2 -
                     50400 * a2 * c * d6 * e * f * n2 - 50400 * a4 * c * d6 * e * f * n2 -
                     16800 * a3 * d7 * e * f * n2 - 10080 * a5 * d7 * e * f * n2 -
                     50400 * b * c4 * d2 * e * f2 * n2 - 201600 * a * b * c3 * d3 * e * f2 * n2 -
                     50400 * b * c2 * d4 * e * f2 * n2 - 302400 * a2 * b * c2 * d4 * e * f2 * n2 -
                     100800 * a * b * c * d5 * e * f2 * n2 -
                     201600 * a3 * b * c * d5 * e * f2 * n2 - 50400 * a2 * b * d6 * e * f2 * n2 -
                     50400 * a4 * b * d6 * e * f2 * n2 - 100800 * b2 * c3 * d2 * e * f3 * n2 -
                     302400 * a * b2 * c2 * d3 * e * f3 * n2 - 50400 * b2 * c * d4 * e * f3 * n2 -
                     302400 * a2 * b2 * c * d4 * e * f3 * n2 - 50400 * a * b2 * d5 * e * f3 * n2 -
                     100800 * a3 * b2 * d5 * e * f3 * n2 - 100800 * b3 * c2 * d2 * e * f4 * n2 -
                     201600 * a * b3 * c * d3 * e * f4 * n2 - 16800 * b3 * d4 * e * f4 * n2 -
                     100800 * a2 * b3 * d4 * e * f4 * n2 - 50400 * b4 * c * d2 * e * f5 * n2 -
                     50400 * a * b4 * d3 * e * f5 * n2 - 10080 * b5 * d2 * e * f6 * n2 -
                     840 * c5 * d3 * g * n2 - 3150 * a * c4 * d4 * g * n2 - 840 * c3 * d5 * g * n2 -
                     5040 * a2 * c3 * d5 * g * n2 - 2100 * a * c2 * d6 * g * n2 -
                     4200 * a3 * c2 * d6 * g * n2 - 1800 * a2 * c * d7 * g * n2 -
                     1800 * a4 * c * d7 * g * n2 - 525 * a3 * d8 * g * n2 - 315 * a5 * d8 * g * n2 -
                     21000 * b * c4 * d3 * f * g * n2 - 63000 * a * b * c3 * d4 * f * g * n2 -
                     12600 * b * c2 * d5 * f * g * n2 - 75600 * a2 * b * c2 * d5 * f * g * n2 -
                     21000 * a * b * c * d6 * f * g * n2 - 42000 * a3 * b * c * d6 * f * g * n2 -
                     9000 * a2 * b * d7 * f * g * n2 - 9000 * a4 * b * d7 * f * g * n2 -
                     75600 * b2 * c3 * d3 * f2 * g * n2 - 170100 * a * b2 * c2 * d4 * f2 * g * n2 -
                     22680 * b2 * c * d5 * f2 * g * n2 - 136080 * a2 * b2 * c * d5 * f2 * g * n2 -
                     18900 * a * b2 * d6 * f2 * g * n2 - 37800 * a3 * b2 * d6 * f2 * g * n2 -
                     109200 * b3 * c2 * d3 * f3 * g * n2 - 163800 * a * b3 * c * d4 * f3 * g * n2 -
                     10920 * b3 * d5 * f3 * g * n2 - 65520 * a2 * b3 * d5 * f3 * g * n2 -
                     71400 * b4 * c * d3 * f4 * g * n2 - 53550 * a * b4 * d4 * f4 * g * n2 -
                     17640 * b5 * d3 * f5 * g * n2 + 5040 * c5 * e * f * n4 +
                     25200 * a * c4 * d * e * f * n4 + 33600 * c3 * d2 * e * f * n4 +
                     50400 * a2 * c3 * d2 * e * f * n4 + 100800 * a * c2 * d3 * e * f * n4 +
                     50400 * a3 * c2 * d3 * e * f * n4 + 25200 * c * d4 * e * f * n4 +
                     100800 * a2 * c * d4 * e * f * n4 + 25200 * a4 * c * d4 * e * f * n4 +
                     25200 * a * d5 * e * f * n4 + 33600 * a3 * d5 * e * f * n4 +
                     5040 * a5 * d5 * e * f * n4 + 25200 * b * c4 * e * f2 * n4 +
                     100800 * a * b * c3 * d * e * f2 * n4 + 100800 * b * c2 * d2 * e * f2 * n4 +
                     151200 * a2 * b * c2 * d2 * e * f2 * n4 +
                     201600 * a * b * c * d3 * e * f2 * n4 +
                     100800 * a3 * b * c * d3 * e * f2 * n4 + 25200 * b * d4 * e * f2 * n4 +
                     100800 * a2 * b * d4 * e * f2 * n4 + 25200 * a4 * b * d4 * e * f2 * n4 +
                     50400 * b2 * c3 * e * f3 * n4 + 151200 * a * b2 * c2 * d * e * f3 * n4 +
                     100800 * b2 * c * d2 * e * f3 * n4 + 151200 * a2 * b2 * c * d2 * e * f3 * n4 +
                     100800 * a * b2 * d3 * e * f3 * n4 + 50400 * a3 * b2 * d3 * e * f3 * n4 +
                     50400 * b3 * c2 * e * f4 * n4 + 100800 * a * b3 * c * d * e * f4 * n4 +
                     33600 * b3 * d2 * e * f4 * n4 + 50400 * a2 * b3 * d2 * e * f4 * n4 +
                     25200 * b4 * c * e * f5 * n4 + 25200 * a * b4 * d * e * f5 * n4 +
                     5040 * b5 * e * f6 * n4 + 1260 * c5 * d * g * n4 +
                     3150 * a * c4 * d2 * g * n4 + 2800 * c3 * d3 * g * n4 +
                     4200 * a2 * c3 * d3 * g * n4 + 6300 * a * c2 * d4 * g * n4 +
                     3150 * a3 * c2 * d4 * g * n4 + 1260 * c * d5 * g * n4 +
                     5040 * a2 * c * d5 * g * n4 + 1260 * a4 * c * d5 * g * n4 +
                     1050 * a * d6 * g * n4 + 1400 * a3 * d6 * g * n4 + 210 * a5 * d6 * g * n4 +
                     31500 * b * c4 * d * f * g * n4 + 63000 * a * b * c3 * d2 * f * g * n4 +
                     42000 * b * c2 * d3 * f * g * n4 + 63000 * a2 * b * c2 * d3 * f * g * n4 +
                     63000 * a * b * c * d4 * f * g * n4 + 31500 * a3 * b * c * d4 * f * g * n4 +
                     6300 * b * d5 * f * g * n4 + 25200 * a2 * b * d5 * f * g * n4 +
                     6300 * a4 * b * d5 * f * g * n4 + 113400 * b2 * c3 * d * f2 * g * n4 +
                     170100 * a * b2 * c2 * d2 * f2 * g * n4 + 75600 * b2 * c * d3 * f2 * g * n4 +
                     113400 * a2 * b2 * c * d3 * f2 * g * n4 + 56700 * a * b2 * d4 * f2 * g * n4 +
                     28350 * a3 * b2 * d4 * f2 * g * n4 + 163800 * b3 * c2 * d * f3 * g * n4 +
                     163800 * a * b3 * c * d2 * f3 * g * n4 + 36400 * b3 * d3 * f3 * g * n4 +
                     54600 * a2 * b3 * d3 * f3 * g * n4 + 107100 * b4 * c * d * f4 * g * n4 +
                     53550 * a * b4 * d2 * f4 * g * n4 + 26460 * b5 * d * f5 * g * n4 -
                     16800 * c3 * e * f * n6 - 50400 * a * c2 * d * e * f * n6 -
                     50400 * c * d2 * e * f * n6 - 50400 * a2 * c * d2 * e * f * n6 -
                     50400 * a * d3 * e * f * n6 - 16800 * a3 * d3 * e * f * n6 -
                     50400 * b * c2 * e * f2 * n6 - 100800 * a * b * c * d * e * f2 * n6 -
                     50400 * b * d2 * e * f2 * n6 - 50400 * a2 * b * d2 * e * f2 * n6 -
                     50400 * b2 * c * e * f3 * n6 - 50400 * a * b2 * d * e * f3 * n6 -
                     16800 * b3 * e * f4 * n6 - 4200 * c3 * d * g * n6 -
                     6300 * a * c2 * d2 * g * n6 - 4200 * c * d3 * g * n6 -
                     4200 * a2 * c * d3 * g * n6 - 3150 * a * d4 * g * n6 -
                     1050 * a3 * d4 * g * n6 - 63000 * b * c2 * d * f * g * n6 -
                     63000 * a * b * c * d2 * f * g * n6 - 21000 * b * d3 * f * g * n6 -
                     21000 * a2 * b * d3 * f * g * n6 - 113400 * b2 * c * d * f2 * g * n6 -
                     56700 * a * b2 * d2 * f2 * g * n6 - 54600 * b3 * d * f3 * g * n6 +
                     25200 * c * e * f * n8 + 25200 * a * d * e * f * n8 + 25200 * b * e * f2 * n8 +
                     6300 * c * d * g * n8 + 3150 * a * d2 * g * n8 + 31500 * b * d * f * g * n8)) /
                     5 +
                 (5 *
                     (1008 * c5 * d3 * e2 * f + 6300 * a * c4 * d4 * e2 * f +
                         15120 * a2 * c3 * d5 * e2 * f + 17640 * a3 * c2 * d6 * e2 * f +
                         10080 * a4 * c * d7 * e2 * f + 2268 * a5 * d8 * e2 * f +
                         5040 * b * c4 * d3 * e2 * f2 + 25200 * a * b * c3 * d4 * e2 * f2 +
                         45360 * a2 * b * c2 * d5 * e2 * f2 + 35280 * a3 * b * c * d6 * e2 * f2 +
                         10080 * a4 * b * d7 * e2 * f2 + 10080 * b2 * c3 * d3 * e2 * f3 +
                         37800 * a * b2 * c2 * d4 * e2 * f3 + 45360 * a2 * b2 * c * d5 * e2 * f3 +
                         17640 * a3 * b2 * d6 * e2 * f3 + 10080 * b3 * c2 * d3 * e2 * f4 +
                         25200 * a * b3 * c * d4 * e2 * f4 + 15120 * a2 * b3 * d5 * e2 * f4 +
                         5040 * b4 * c * d3 * e2 * f5 + 6300 * a * b4 * d4 * e2 * f5 +
                         1008 * b5 * d3 * e2 * f6 + 126 * c5 * d4 * e * g +
                         630 * a * c4 * d5 * e * g + 1260 * a2 * c3 * d6 * e * g +
                         1260 * a3 * c2 * d7 * e * g + 630 * a4 * c * d8 * e * g +
                         126 * a5 * d9 * e * g + 3150 * b * c4 * d4 * e * f * g +
                         12600 * a * b * c3 * d5 * e * f * g +
                         18900 * a2 * b * c2 * d6 * e * f * g +
                         12600 * a3 * b * c * d7 * e * f * g + 3150 * a4 * b * d8 * e * f * g +
                         11340 * b2 * c3 * d4 * e * f2 * g + 34020 * a * b2 * c2 * d5 * e * f2 * g +
                         34020 * a2 * b2 * c * d6 * e * f2 * g + 11340 * a3 * b2 * d7 * e * f2 * g +
                         16380 * b3 * c2 * d4 * e * f3 * g + 32760 * a * b3 * c * d5 * e * f3 * g +
                         16380 * a2 * b3 * d6 * e * f3 * g + 10710 * b4 * c * d4 * e * f4 * g +
                         10710 * a * b4 * d5 * e * f4 * g + 2646 * b5 * d4 * e * f5 * g +
                         126 * b * c4 * d5 * g2 + 420 * a * b * c3 * d6 * g2 +
                         540 * a2 * b * c2 * d7 * g2 + 315 * a3 * b * c * d8 * g2 +
                         70 * a4 * b * d9 * g2 + 1512 * b2 * c3 * d5 * f * g2 +
                         3780 * a * b2 * c2 * d6 * f * g2 + 3240 * a2 * b2 * c * d7 * f * g2 +
                         945 * a3 * b2 * d8 * f * g2 + 3780 * b3 * c2 * d5 * f2 * g2 +
                         6300 * a * b3 * c * d6 * f2 * g2 + 2700 * a2 * b3 * d7 * f2 * g2 +
                         3528 * b4 * c * d5 * f3 * g2 + 2940 * a * b4 * d6 * f3 * g2 +
                         1134 * b5 * d5 * f4 * g2 - 1008 * c5 * d * e2 * f * n2 -
                         7560 * a * c4 * d2 * e2 * f * n2 - 3360 * c3 * d3 * e2 * f * n2 -
                         20160 * a2 * c3 * d3 * e2 * f * n2 - 12600 * a * c2 * d4 * e2 * f * n2 -
                         25200 * a3 * c2 * d4 * e2 * f * n2 - 15120 * a2 * c * d5 * e2 * f * n2 -
                         15120 * a4 * c * d5 * e2 * f * n2 - 5880 * a3 * d6 * e2 * f * n2 -
                         3528 * a5 * d6 * e2 * f * n2 - 5040 * b * c4 * d * e2 * f2 * n2 -
                         30240 * a * b * c3 * d2 * e2 * f2 * n2 -
                         10080 * b * c2 * d3 * e2 * f2 * n2 -
                         60480 * a2 * b * c2 * d3 * e2 * f2 * n2 -
                         25200 * a * b * c * d4 * e2 * f2 * n2 -
                         50400 * a3 * b * c * d4 * e2 * f2 * n2 -
                         15120 * a2 * b * d5 * e2 * f2 * n2 - 15120 * a4 * b * d5 * e2 * f2 * n2 -
                         10080 * b2 * c3 * d * e2 * f3 * n2 -
                         45360 * a * b2 * c2 * d2 * e2 * f3 * n2 -
                         10080 * b2 * c * d3 * e2 * f3 * n2 -
                         60480 * a2 * b2 * c * d3 * e2 * f3 * n2 -
                         12600 * a * b2 * d4 * e2 * f3 * n2 - 25200 * a3 * b2 * d4 * e2 * f3 * n2 -
                         10080 * b3 * c2 * d * e2 * f4 * n2 -
                         30240 * a * b3 * c * d2 * e2 * f4 * n2 - 3360 * b3 * d3 * e2 * f4 * n2 -
                         20160 * a2 * b3 * d3 * e2 * f4 * n2 - 5040 * b4 * c * d * e2 * f5 * n2 -
                         7560 * a * b4 * d2 * e2 * f5 * n2 - 1008 * b5 * d * e2 * f6 * n2 -
                         252 * c5 * d2 * e * g * n2 - 1260 * a * c4 * d3 * e * g * n2 -
                         420 * c3 * d4 * e * g * n2 - 2520 * a2 * c3 * d4 * e * g * n2 -
                         1260 * a * c2 * d5 * e * g * n2 - 2520 * a3 * c2 * d5 * e * g * n2 -
                         1260 * a2 * c * d6 * e * g * n2 - 1260 * a4 * c * d6 * e * g * n2 -
                         420 * a3 * d7 * e * g * n2 - 252 * a5 * d7 * e * g * n2 -
                         6300 * b * c4 * d2 * e * f * g * n2 -
                         25200 * a * b * c3 * d3 * e * f * g * n2 -
                         6300 * b * c2 * d4 * e * f * g * n2 -
                         37800 * a2 * b * c2 * d4 * e * f * g * n2 -
                         12600 * a * b * c * d5 * e * f * g * n2 -
                         25200 * a3 * b * c * d5 * e * f * g * n2 -
                         6300 * a2 * b * d6 * e * f * g * n2 - 6300 * a4 * b * d6 * e * f * g * n2 -
                         22680 * b2 * c3 * d2 * e * f2 * g * n2 -
                         68040 * a * b2 * c2 * d3 * e * f2 * g * n2 -
                         11340 * b2 * c * d4 * e * f2 * g * n2 -
                         68040 * a2 * b2 * c * d4 * e * f2 * g * n2 -
                         11340 * a * b2 * d5 * e * f2 * g * n2 -
                         22680 * a3 * b2 * d5 * e * f2 * g * n2 -
                         32760 * b3 * c2 * d2 * e * f3 * g * n2 -
                         65520 * a * b3 * c * d3 * e * f3 * g * n2 -
                         5460 * b3 * d4 * e * f3 * g * n2 - 32760 * a2 * b3 * d4 * e * f3 * g * n2 -
                         21420 * b4 * c * d2 * e * f4 * g * n2 -
                         21420 * a * b4 * d3 * e * f4 * g * n2 - 5292 * b5 * d2 * e * f5 * g * n2 -
                         420 * b * c4 * d3 * g2 * n2 - 1260 * a * b * c3 * d4 * g2 * n2 -
                         252 * b * c2 * d5 * g2 * n2 - 1512 * a2 * b * c2 * d5 * g2 * n2 -
                         420 * a * b * c * d6 * g2 * n2 - 840 * a3 * b * c * d6 * g2 * n2 -
                         180 * a2 * b * d7 * g2 * n2 - 180 * a4 * b * d7 * g2 * n2 -
                         5040 * b2 * c3 * d3 * f * g2 * n2 -
                         11340 * a * b2 * c2 * d4 * f * g2 * n2 - 1512 * b2 * c * d5 * f * g2 * n2 -
                         9072 * a2 * b2 * c * d5 * f * g2 * n2 - 1260 * a * b2 * d6 * f * g2 * n2 -
                         2520 * a3 * b2 * d6 * f * g2 * n2 - 12600 * b3 * c2 * d3 * f2 * g2 * n2 -
                         18900 * a * b3 * c * d4 * f2 * g2 * n2 - 1260 * b3 * d5 * f2 * g2 * n2 -
                         7560 * a2 * b3 * d5 * f2 * g2 * n2 - 11760 * b4 * c * d3 * f3 * g2 * n2 -
                         8820 * a * b4 * d4 * f3 * g2 * n2 - 3780 * b5 * d3 * f4 * g2 * n2 +
                         1260 * a * c4 * e2 * f * n4 + 3360 * c3 * d * e2 * f * n4 +
                         5040 * a2 * c3 * d * e2 * f * n4 + 15120 * a * c2 * d2 * e2 * f * n4 +
                         7560 * a3 * c2 * d2 * e2 * f * n4 + 5040 * c * d3 * e2 * f * n4 +
                         20160 * a2 * c * d3 * e2 * f * n4 + 5040 * a4 * c * d3 * e2 * f * n4 +
                         6300 * a * d4 * e2 * f * n4 + 8400 * a3 * d4 * e2 * f * n4 +
                         1260 * a5 * d4 * e2 * f * n4 + 5040 * a * b * c3 * e2 * f2 * n4 +
                         10080 * b * c2 * d * e2 * f2 * n4 +
                         15120 * a2 * b * c2 * d * e2 * f2 * n4 +
                         30240 * a * b * c * d2 * e2 * f2 * n4 +
                         15120 * a3 * b * c * d2 * e2 * f2 * n4 + 5040 * b * d3 * e2 * f2 * n4 +
                         20160 * a2 * b * d3 * e2 * f2 * n4 + 5040 * a4 * b * d3 * e2 * f2 * n4 +
                         7560 * a * b2 * c2 * e2 * f3 * n4 + 10080 * b2 * c * d * e2 * f3 * n4 +
                         15120 * a2 * b2 * c * d * e2 * f3 * n4 +
                         15120 * a * b2 * d2 * e2 * f3 * n4 + 7560 * a3 * b2 * d2 * e2 * f3 * n4 +
                         5040 * a * b3 * c * e2 * f4 * n4 + 3360 * b3 * d * e2 * f4 * n4 +
                         5040 * a2 * b3 * d * e2 * f4 * n4 + 1260 * a * b4 * e2 * f5 * n4 +
                         126 * c5 * e * g * n4 + 630 * a * c4 * d * e * g * n4 +
                         840 * c3 * d2 * e * g * n4 + 1260 * a2 * c3 * d2 * e * g * n4 +
                         2520 * a * c2 * d3 * e * g * n4 + 1260 * a3 * c2 * d3 * e * g * n4 +
                         630 * c * d4 * e * g * n4 + 2520 * a2 * c * d4 * e * g * n4 +
                         630 * a4 * c * d4 * e * g * n4 + 630 * a * d5 * e * g * n4 +
                         840 * a3 * d5 * e * g * n4 + 126 * a5 * d5 * e * g * n4 +
                         3150 * b * c4 * e * f * g * n4 + 12600 * a * b * c3 * d * e * f * g * n4 +
                         12600 * b * c2 * d2 * e * f * g * n4 +
                         18900 * a2 * b * c2 * d2 * e * f * g * n4 +
                         25200 * a * b * c * d3 * e * f * g * n4 +
                         12600 * a3 * b * c * d3 * e * f * g * n4 + 3150 * b * d4 * e * f * g * n4 +
                         12600 * a2 * b * d4 * e * f * g * n4 +
                         3150 * a4 * b * d4 * e * f * g * n4 + 11340 * b2 * c3 * e * f2 * g * n4 +
                         34020 * a * b2 * c2 * d * e * f2 * g * n4 +
                         22680 * b2 * c * d2 * e * f2 * g * n4 +
                         34020 * a2 * b2 * c * d2 * e * f2 * g * n4 +
                         22680 * a * b2 * d3 * e * f2 * g * n4 +
                         11340 * a3 * b2 * d3 * e * f2 * g * n4 +
                         16380 * b3 * c2 * e * f3 * g * n4 +
                         32760 * a * b3 * c * d * e * f3 * g * n4 +
                         10920 * b3 * d2 * e * f3 * g * n4 +
                         16380 * a2 * b3 * d2 * e * f3 * g * n4 + 10710 * b4 * c * e * f4 * g * n4 +
                         10710 * a * b4 * d * e * f4 * g * n4 + 2646 * b5 * e * f5 * g * n4 +
                         630 * b * c4 * d * g2 * n4 + 1260 * a * b * c3 * d2 * g2 * n4 +
                         840 * b * c2 * d3 * g2 * n4 + 1260 * a2 * b * c2 * d3 * g2 * n4 +
                         1260 * a * b * c * d4 * g2 * n4 + 630 * a3 * b * c * d4 * g2 * n4 +
                         126 * b * d5 * g2 * n4 + 504 * a2 * b * d5 * g2 * n4 +
                         126 * a4 * b * d5 * g2 * n4 + 7560 * b2 * c3 * d * f * g2 * n4 +
                         11340 * a * b2 * c2 * d2 * f * g2 * n4 + 5040 * b2 * c * d3 * f * g2 * n4 +
                         7560 * a2 * b2 * c * d3 * f * g2 * n4 + 3780 * a * b2 * d4 * f * g2 * n4 +
                         1890 * a3 * b2 * d4 * f * g2 * n4 + 18900 * b3 * c2 * d * f2 * g2 * n4 +
                         18900 * a * b3 * c * d2 * f2 * g2 * n4 + 4200 * b3 * d3 * f2 * g2 * n4 +
                         6300 * a2 * b3 * d3 * f2 * g2 * n4 + 17640 * b4 * c * d * f3 * g2 * n4 +
                         8820 * a * b4 * d2 * f3 * g2 * n4 + 5670 * b5 * d * f4 * g2 * n4 -
                         2520 * a * c2 * e2 * f * n6 - 5040 * c * d * e2 * f * n6 -
                         5040 * a2 * c * d * e2 * f * n6 - 7560 * a * d2 * e2 * f * n6 -
                         2520 * a3 * d2 * e2 * f * n6 - 5040 * a * b * c * e2 * f2 * n6 -
                         5040 * b * d * e2 * f2 * n6 - 5040 * a2 * b * d * e2 * f2 * n6 -
                         2520 * a * b2 * e2 * f3 * n6 - 420 * c3 * e * g * n6 -
                         1260 * a * c2 * d * e * g * n6 - 1260 * c * d2 * e * g * n6 -
                         1260 * a2 * c * d2 * e * g * n6 - 1260 * a * d3 * e * g * n6 -
                         420 * a3 * d3 * e * g * n6 - 6300 * b * c2 * e * f * g * n6 -
                         12600 * a * b * c * d * e * f * g * n6 - 6300 * b * d2 * e * f * g * n6 -
                         6300 * a2 * b * d2 * e * f * g * n6 - 11340 * b2 * c * e * f2 * g * n6 -
                         11340 * a * b2 * d * e * f2 * g * n6 - 5460 * b3 * e * f3 * g * n6 -
                         1260 * b * c2 * d * g2 * n6 - 1260 * a * b * c * d2 * g2 * n6 -
                         420 * b * d3 * g2 * n6 - 420 * a2 * b * d3 * g2 * n6 -
                         7560 * b2 * c * d * f * g2 * n6 - 3780 * a * b2 * d2 * f * g2 * n6 -
                         6300 * b3 * d * f2 * g2 * n6 + 1260 * a * e2 * f * n8 +
                         630 * c * e * g * n8 + 630 * a * d * e * g * n8 +
                         3150 * b * e * f * g * n8 + 630 * b * d * g2 * n8)) /
                     3 +
                 (5 *
                     (2016 * c5 * d2 * e3 * f + 16800 * a * c4 * d3 * e3 * f +
                         50400 * a2 * c3 * d4 * e3 * f + 70560 * a3 * c2 * d5 * e3 * f +
                         47040 * a4 * c * d6 * e3 * f + 12096 * a5 * d7 * e3 * f +
                         10080 * b * c4 * d2 * e3 * f2 + 67200 * a * b * c3 * d3 * e3 * f2 +
                         151200 * a2 * b * c2 * d4 * e3 * f2 + 141120 * a3 * b * c * d5 * e3 * f2 +
                         47040 * a4 * b * d6 * e3 * f2 + 20160 * b2 * c3 * d2 * e3 * f3 +
                         100800 * a * b2 * c2 * d3 * e3 * f3 + 151200 * a2 * b2 * c * d4 * e3 * f3 +
                         70560 * a3 * b2 * d5 * e3 * f3 + 20160 * b3 * c2 * d2 * e3 * f4 +
                         67200 * a * b3 * c * d3 * e3 * f4 + 50400 * a2 * b3 * d4 * e3 * f4 +
                         10080 * b4 * c * d2 * e3 * f5 + 16800 * a * b4 * d3 * e3 * f5 +
                         2016 * b5 * d2 * e3 * f6 + 504 * c5 * d3 * e2 * g +
                         3150 * a * c4 * d4 * e2 * g + 7560 * a2 * c3 * d5 * e2 * g +
                         8820 * a3 * c2 * d6 * e2 * g + 5040 * a4 * c * d7 * e2 * g +
                         1134 * a5 * d8 * e2 * g + 12600 * b * c4 * d3 * e2 * f * g +
                         63000 * a * b * c3 * d4 * e2 * f * g +
                         113400 * a2 * b * c2 * d5 * e2 * f * g +
                         88200 * a3 * b * c * d6 * e2 * f * g + 25200 * a4 * b * d7 * e2 * f * g +
                         45360 * b2 * c3 * d3 * e2 * f2 * g +
                         170100 * a * b2 * c2 * d4 * e2 * f2 * g +
                         204120 * a2 * b2 * c * d5 * e2 * f2 * g +
                         79380 * a3 * b2 * d6 * e2 * f2 * g + 65520 * b3 * c2 * d3 * e2 * f3 * g +
                         163800 * a * b3 * c * d4 * e2 * f3 * g +
                         98280 * a2 * b3 * d5 * e2 * f3 * g + 42840 * b4 * c * d3 * e2 * f4 * g +
                         53550 * a * b4 * d4 * e2 * f4 * g + 10584 * b5 * d3 * e2 * f5 * g +
                         1260 * b * c4 * d4 * e * g2 + 5040 * a * b * c3 * d5 * e * g2 +
                         7560 * a2 * b * c2 * d6 * e * g2 + 5040 * a3 * b * c * d7 * e * g2 +
                         1260 * a4 * b * d8 * e * g2 + 15120 * b2 * c3 * d4 * e * f * g2 +
                         45360 * a * b2 * c2 * d5 * e * f * g2 +
                         45360 * a2 * b2 * c * d6 * e * f * g2 + 15120 * a3 * b2 * d7 * e * f * g2 +
                         37800 * b3 * c2 * d4 * e * f2 * g2 +
                         75600 * a * b3 * c * d5 * e * f2 * g2 +
                         37800 * a2 * b3 * d6 * e * f2 * g2 + 35280 * b4 * c * d4 * e * f3 * g2 +
                         35280 * a * b4 * d5 * e * f3 * g2 + 11340 * b5 * d4 * e * f4 * g2 +
                         504 * b2 * c3 * d5 * g3 + 1260 * a * b2 * c2 * d6 * g3 +
                         1080 * a2 * b2 * c * d7 * g3 + 315 * a3 * b2 * d8 * g3 +
                         3528 * b3 * c2 * d5 * f * g3 + 5880 * a * b3 * c * d6 * f * g3 +
                         2520 * a2 * b3 * d7 * f * g3 + 5544 * b4 * c * d5 * f2 * g3 +
                         4620 * a * b4 * d6 * f2 * g3 + 2520 * b5 * d5 * f3 * g3 -
                         672 * c5 * e3 * f * n2 - 10080 * a * c4 * d * e3 * f * n2 -
                         6720 * c3 * d2 * e3 * f * n2 - 40320 * a2 * c3 * d2 * e3 * f * n2 -
                         33600 * a * c2 * d3 * e3 * f * n2 - 67200 * a3 * c2 * d3 * e3 * f * n2 -
                         50400 * a2 * c * d4 * e3 * f * n2 - 50400 * a4 * c * d4 * e3 * f * n2 -
                         23520 * a3 * d5 * e3 * f * n2 - 14112 * a5 * d5 * e3 * f * n2 -
                         3360 * b * c4 * e3 * f2 * n2 - 40320 * a * b * c3 * d * e3 * f2 * n2 -
                         20160 * b * c2 * d2 * e3 * f2 * n2 -
                         120960 * a2 * b * c2 * d2 * e3 * f2 * n2 -
                         67200 * a * b * c * d3 * e3 * f2 * n2 -
                         134400 * a3 * b * c * d3 * e3 * f2 * n2 -
                         50400 * a2 * b * d4 * e3 * f2 * n2 - 50400 * a4 * b * d4 * e3 * f2 * n2 -
                         6720 * b2 * c3 * e3 * f3 * n2 - 60480 * a * b2 * c2 * d * e3 * f3 * n2 -
                         20160 * b2 * c * d2 * e3 * f3 * n2 -
                         120960 * a2 * b2 * c * d2 * e3 * f3 * n2 -
                         33600 * a * b2 * d3 * e3 * f3 * n2 - 67200 * a3 * b2 * d3 * e3 * f3 * n2 -
                         6720 * b3 * c2 * e3 * f4 * n2 - 40320 * a * b3 * c * d * e3 * f4 * n2 -
                         6720 * b3 * d2 * e3 * f4 * n2 - 40320 * a2 * b3 * d2 * e3 * f4 * n2 -
                         3360 * b4 * c * e3 * f5 * n2 - 10080 * a * b4 * d * e3 * f5 * n2 -
                         672 * b5 * e3 * f6 * n2 - 504 * c5 * d * e2 * g * n2 -
                         3780 * a * c4 * d2 * e2 * g * n2 - 1680 * c3 * d3 * e2 * g * n2 -
                         10080 * a2 * c3 * d3 * e2 * g * n2 - 6300 * a * c2 * d4 * e2 * g * n2 -
                         12600 * a3 * c2 * d4 * e2 * g * n2 - 7560 * a2 * c * d5 * e2 * g * n2 -
                         7560 * a4 * c * d5 * e2 * g * n2 - 2940 * a3 * d6 * e2 * g * n2 -
                         1764 * a5 * d6 * e2 * g * n2 - 12600 * b * c4 * d * e2 * f * g * n2 -
                         75600 * a * b * c3 * d2 * e2 * f * g * n2 -
                         25200 * b * c2 * d3 * e2 * f * g * n2 -
                         151200 * a2 * b * c2 * d3 * e2 * f * g * n2 -
                         63000 * a * b * c * d4 * e2 * f * g * n2 -
                         126000 * a3 * b * c * d4 * e2 * f * g * n2 -
                         37800 * a2 * b * d5 * e2 * f * g * n2 -
                         37800 * a4 * b * d5 * e2 * f * g * n2 -
                         45360 * b2 * c3 * d * e2 * f2 * g * n2 -
                         204120 * a * b2 * c2 * d2 * e2 * f2 * g * n2 -
                         45360 * b2 * c * d3 * e2 * f2 * g * n2 -
                         272160 * a2 * b2 * c * d3 * e2 * f2 * g * n2 -
                         56700 * a * b2 * d4 * e2 * f2 * g * n2 -
                         113400 * a3 * b2 * d4 * e2 * f2 * g * n2 -
                         65520 * b3 * c2 * d * e2 * f3 * g * n2 -
                         196560 * a * b3 * c * d2 * e2 * f3 * g * n2 -
                         21840 * b3 * d3 * e2 * f3 * g * n2 -
                         131040 * a2 * b3 * d3 * e2 * f3 * g * n2 -
                         42840 * b4 * c * d * e2 * f4 * g * n2 -
                         64260 * a * b4 * d2 * e2 * f4 * g * n2 -
                         10584 * b5 * d * e2 * f5 * g * n2 - 2520 * b * c4 * d2 * e * g2 * n2 -
                         10080 * a * b * c3 * d3 * e * g2 * n2 - 2520 * b * c2 * d4 * e * g2 * n2 -
                         15120 * a2 * b * c2 * d4 * e * g2 * n2 -
                         5040 * a * b * c * d5 * e * g2 * n2 -
                         10080 * a3 * b * c * d5 * e * g2 * n2 - 2520 * a2 * b * d6 * e * g2 * n2 -
                         2520 * a4 * b * d6 * e * g2 * n2 - 30240 * b2 * c3 * d2 * e * f * g2 * n2 -
                         90720 * a * b2 * c2 * d3 * e * f * g2 * n2 -
                         15120 * b2 * c * d4 * e * f * g2 * n2 -
                         90720 * a2 * b2 * c * d4 * e * f * g2 * n2 -
                         15120 * a * b2 * d5 * e * f * g2 * n2 -
                         30240 * a3 * b2 * d5 * e * f * g2 * n2 -
                         75600 * b3 * c2 * d2 * e * f2 * g2 * n2 -
                         151200 * a * b3 * c * d3 * e * f2 * g2 * n2 -
                         12600 * b3 * d4 * e * f2 * g2 * n2 -
                         75600 * a2 * b3 * d4 * e * f2 * g2 * n2 -
                         70560 * b4 * c * d2 * e * f3 * g2 * n2 -
                         70560 * a * b4 * d3 * e * f3 * g2 * n2 -
                         22680 * b5 * d2 * e * f4 * g2 * n2 - 1680 * b2 * c3 * d3 * g3 * n2 -
                         3780 * a * b2 * c2 * d4 * g3 * n2 - 504 * b2 * c * d5 * g3 * n2 -
                         3024 * a2 * b2 * c * d5 * g3 * n2 - 420 * a * b2 * d6 * g3 * n2 -
                         840 * a3 * b2 * d6 * g3 * n2 - 11760 * b3 * c2 * d3 * f * g3 * n2 -
                         17640 * a * b3 * c * d4 * f * g3 * n2 - 1176 * b3 * d5 * f * g3 * n2 -
                         7056 * a2 * b3 * d5 * f * g3 * n2 - 18480 * b4 * c * d3 * f2 * g3 * n2 -
                         13860 * a * b4 * d4 * f2 * g3 * n2 - 8400 * b5 * d3 * f3 * g3 * n2 +
                         2240 * c3 * e3 * f * n4 + 3360 * a2 * c3 * e3 * f * n4 +
                         20160 * a * c2 * d * e3 * f * n4 + 10080 * a3 * c2 * d * e3 * f * n4 +
                         10080 * c * d2 * e3 * f * n4 + 40320 * a2 * c * d2 * e3 * f * n4 +
                         10080 * a4 * c * d2 * e3 * f * n4 + 16800 * a * d3 * e3 * f * n4 +
                         22400 * a3 * d3 * e3 * f * n4 + 3360 * a5 * d3 * e3 * f * n4 +
                         6720 * b * c2 * e3 * f2 * n4 + 10080 * a2 * b * c2 * e3 * f2 * n4 +
                         40320 * a * b * c * d * e3 * f2 * n4 +
                         20160 * a3 * b * c * d * e3 * f2 * n4 + 10080 * b * d2 * e3 * f2 * n4 +
                         40320 * a2 * b * d2 * e3 * f2 * n4 + 10080 * a4 * b * d2 * e3 * f2 * n4 +
                         6720 * b2 * c * e3 * f3 * n4 + 10080 * a2 * b2 * c * e3 * f3 * n4 +
                         20160 * a * b2 * d * e3 * f3 * n4 + 10080 * a3 * b2 * d * e3 * f3 * n4 +
                         2240 * b3 * e3 * f4 * n4 + 3360 * a2 * b3 * e3 * f4 * n4 +
                         630 * a * c4 * e2 * g * n4 + 1680 * c3 * d * e2 * g * n4 +
                         2520 * a2 * c3 * d * e2 * g * n4 + 7560 * a * c2 * d2 * e2 * g * n4 +
                         3780 * a3 * c2 * d2 * e2 * g * n4 + 2520 * c * d3 * e2 * g * n4 +
                         10080 * a2 * c * d3 * e2 * g * n4 + 2520 * a4 * c * d3 * e2 * g * n4 +
                         3150 * a * d4 * e2 * g * n4 + 4200 * a3 * d4 * e2 * g * n4 +
                         630 * a5 * d4 * e2 * g * n4 + 12600 * a * b * c3 * e2 * f * g * n4 +
                         25200 * b * c2 * d * e2 * f * g * n4 +
                         37800 * a2 * b * c2 * d * e2 * f * g * n4 +
                         75600 * a * b * c * d2 * e2 * f * g * n4 +
                         37800 * a3 * b * c * d2 * e2 * f * g * n4 +
                         12600 * b * d3 * e2 * f * g * n4 + 50400 * a2 * b * d3 * e2 * f * g * n4 +
                         12600 * a4 * b * d3 * e2 * f * g * n4 +
                         34020 * a * b2 * c2 * e2 * f2 * g * n4 +
                         45360 * b2 * c * d * e2 * f2 * g * n4 +
                         68040 * a2 * b2 * c * d * e2 * f2 * g * n4 +
                         68040 * a * b2 * d2 * e2 * f2 * g * n4 +
                         34020 * a3 * b2 * d2 * e2 * f2 * g * n4 +
                         32760 * a * b3 * c * e2 * f3 * g * n4 + 21840 * b3 * d * e2 * f3 * g * n4 +
                         32760 * a2 * b3 * d * e2 * f3 * g * n4 +
                         10710 * a * b4 * e2 * f4 * g * n4 + 1260 * b * c4 * e * g2 * n4 +
                         5040 * a * b * c3 * d * e * g2 * n4 + 5040 * b * c2 * d2 * e * g2 * n4 +
                         7560 * a2 * b * c2 * d2 * e * g2 * n4 +
                         10080 * a * b * c * d3 * e * g2 * n4 +
                         5040 * a3 * b * c * d3 * e * g2 * n4 + 1260 * b * d4 * e * g2 * n4 +
                         5040 * a2 * b * d4 * e * g2 * n4 + 1260 * a4 * b * d4 * e * g2 * n4 +
                         15120 * b2 * c3 * e * f * g2 * n4 +
                         45360 * a * b2 * c2 * d * e * f * g2 * n4 +
                         30240 * b2 * c * d2 * e * f * g2 * n4 +
                         45360 * a2 * b2 * c * d2 * e * f * g2 * n4 +
                         30240 * a * b2 * d3 * e * f * g2 * n4 +
                         15120 * a3 * b2 * d3 * e * f * g2 * n4 +
                         37800 * b3 * c2 * e * f2 * g2 * n4 +
                         75600 * a * b3 * c * d * e * f2 * g2 * n4 +
                         25200 * b3 * d2 * e * f2 * g2 * n4 +
                         37800 * a2 * b3 * d2 * e * f2 * g2 * n4 +
                         35280 * b4 * c * e * f3 * g2 * n4 + 35280 * a * b4 * d * e * f3 * g2 * n4 +
                         11340 * b5 * e * f4 * g2 * n4 + 2520 * b2 * c3 * d * g3 * n4 +
                         3780 * a * b2 * c2 * d2 * g3 * n4 + 1680 * b2 * c * d3 * g3 * n4 +
                         2520 * a2 * b2 * c * d3 * g3 * n4 + 1260 * a * b2 * d4 * g3 * n4 +
                         630 * a3 * b2 * d4 * g3 * n4 + 17640 * b3 * c2 * d * f * g3 * n4 +
                         17640 * a * b3 * c * d2 * f * g3 * n4 + 3920 * b3 * d3 * f * g3 * n4 +
                         5880 * a2 * b3 * d3 * f * g3 * n4 + 27720 * b4 * c * d * f2 * g3 * n4 +
                         13860 * a * b4 * d2 * f2 * g3 * n4 + 12600 * b5 * d * f3 * g3 * n4 -
                         3360 * c * e3 * f * n6 - 3360 * a2 * c * e3 * f * n6 -
                         10080 * a * d * e3 * f * n6 - 3360 * a3 * d * e3 * f * n6 -
                         3360 * b * e3 * f2 * n6 - 3360 * a2 * b * e3 * f2 * n6 -
                         1260 * a * c2 * e2 * g * n6 - 2520 * c * d * e2 * g * n6 -
                         2520 * a2 * c * d * e2 * g * n6 - 3780 * a * d2 * e2 * g * n6 -
                         1260 * a3 * d2 * e2 * g * n6 - 12600 * a * b * c * e2 * f * g * n6 -
                         12600 * b * d * e2 * f * g * n6 - 12600 * a2 * b * d * e2 * f * g * n6 -
                         11340 * a * b2 * e2 * f2 * g * n6 - 2520 * b * c2 * e * g2 * n6 -
                         5040 * a * b * c * d * e * g2 * n6 - 2520 * b * d2 * e * g2 * n6 -
                         2520 * a2 * b * d2 * e * g2 * n6 - 15120 * b2 * c * e * f * g2 * n6 -
                         15120 * a * b2 * d * e * f * g2 * n6 - 12600 * b3 * e * f2 * g2 * n6 -
                         2520 * b2 * c * d * g3 * n6 - 1260 * a * b2 * d2 * g3 * n6 -
                         5880 * b3 * d * f * g3 * n6 + 630 * a * e2 * g * n8 +
                         1260 * b * e * g2 * n8)) /
                     7 +
                 5 * (126 * c5 * d * e4 * f + 1575 * a * c4 * d2 * e4 * f +
                         6300 * a2 * c3 * d3 * e4 * f + 11025 * a3 * c2 * d4 * e4 * f +
                         8820 * a4 * c * d5 * e4 * f + 2646 * a5 * d6 * e4 * f +
                         630 * b * c4 * d * e4 * f2 + 6300 * a * b * c3 * d2 * e4 * f2 +
                         18900 * a2 * b * c2 * d3 * e4 * f2 + 22050 * a3 * b * c * d4 * e4 * f2 +
                         8820 * a4 * b * d5 * e4 * f2 + 1260 * b2 * c3 * d * e4 * f3 +
                         9450 * a * b2 * c2 * d2 * e4 * f3 + 18900 * a2 * b2 * c * d3 * e4 * f3 +
                         11025 * a3 * b2 * d4 * e4 * f3 + 1260 * b3 * c2 * d * e4 * f4 +
                         6300 * a * b3 * c * d2 * e4 * f4 + 6300 * a2 * b3 * d3 * e4 * f4 +
                         630 * b4 * c * d * e4 * f5 + 1575 * a * b4 * d2 * e4 * f5 +
                         126 * b5 * d * e4 * f6 + 63 * c5 * d2 * e3 * g +
                         525 * a * c4 * d3 * e3 * g + 1575 * a2 * c3 * d4 * e3 * g +
                         2205 * a3 * c2 * d5 * e3 * g + 1470 * a4 * c * d6 * e3 * g +
                         378 * a5 * d7 * e3 * g + 1575 * b * c4 * d2 * e3 * f * g +
                         10500 * a * b * c3 * d3 * e3 * f * g +
                         23625 * a2 * b * c2 * d4 * e3 * f * g +
                         22050 * a3 * b * c * d5 * e3 * f * g + 7350 * a4 * b * d6 * e3 * f * g +
                         5670 * b2 * c3 * d2 * e3 * f2 * g +
                         28350 * a * b2 * c2 * d3 * e3 * f2 * g +
                         42525 * a2 * b2 * c * d4 * e3 * f2 * g +
                         19845 * a3 * b2 * d5 * e3 * f2 * g + 8190 * b3 * c2 * d2 * e3 * f3 * g +
                         27300 * a * b3 * c * d3 * e3 * f3 * g +
                         20475 * a2 * b3 * d4 * e3 * f3 * g + 5355 * b4 * c * d2 * e3 * f4 * g +
                         8925 * a * b4 * d3 * e3 * f4 * g + 1323 * b5 * d2 * e3 * f5 * g +
                         315 * b * c4 * d3 * e2 * g2 + 1575 * a * b * c3 * d4 * e2 * g2 +
                         2835 * a2 * b * c2 * d5 * e2 * g2 + 2205 * a3 * b * c * d6 * e2 * g2 +
                         630 * a4 * b * d7 * e2 * g2 + 3780 * b2 * c3 * d3 * e2 * f * g2 +
                         14175 * a * b2 * c2 * d4 * e2 * f * g2 +
                         17010 * a2 * b2 * c * d5 * e2 * f * g2 +
                         6615 * a3 * b2 * d6 * e2 * f * g2 + 9450 * b3 * c2 * d3 * e2 * f2 * g2 +
                         23625 * a * b3 * c * d4 * e2 * f2 * g2 +
                         14175 * a2 * b3 * d5 * e2 * f2 * g2 + 8820 * b4 * c * d3 * e2 * f3 * g2 +
                         11025 * a * b4 * d4 * e2 * f3 * g2 + 2835 * b5 * d3 * e2 * f4 * g2 +
                         315 * b2 * c3 * d4 * e * g3 + 945 * a * b2 * c2 * d5 * e * g3 +
                         945 * a2 * b2 * c * d6 * e * g3 + 315 * a3 * b2 * d7 * e * g3 +
                         2205 * b3 * c2 * d4 * e * f * g3 + 4410 * a * b3 * c * d5 * e * f * g3 +
                         2205 * a2 * b3 * d6 * e * f * g3 + 3465 * b4 * c * d4 * e * f2 * g3 +
                         3465 * a * b4 * d5 * e * f2 * g3 + 1575 * b5 * d4 * e * f3 * g3 +
                         63 * b3 * c2 * d5 * g4 + 105 * a * b3 * c * d6 * g4 +
                         45 * a2 * b3 * d7 * g4 + 252 * b4 * c * d5 * f * g4 +
                         210 * a * b4 * d6 * f * g4 + 189 * b5 * d5 * f2 * g4 -
                         315 * a * c4 * e4 * f * n2 - 420 * c3 * d * e4 * f * n2 -
                         2520 * a2 * c3 * d * e4 * f * n2 - 3150 * a * c2 * d2 * e4 * f * n2 -
                         6300 * a3 * c2 * d2 * e4 * f * n2 - 6300 * a2 * c * d3 * e4 * f * n2 -
                         6300 * a4 * c * d3 * e4 * f * n2 - 3675 * a3 * d4 * e4 * f * n2 -
                         2205 * a5 * d4 * e4 * f * n2 - 1260 * a * b * c3 * e4 * f2 * n2 -
                         1260 * b * c2 * d * e4 * f2 * n2 - 7560 * a2 * b * c2 * d * e4 * f2 * n2 -
                         6300 * a * b * c * d2 * e4 * f2 * n2 -
                         12600 * a3 * b * c * d2 * e4 * f2 * n2 -
                         6300 * a2 * b * d3 * e4 * f2 * n2 - 6300 * a4 * b * d3 * e4 * f2 * n2 -
                         1890 * a * b2 * c2 * e4 * f3 * n2 - 1260 * b2 * c * d * e4 * f3 * n2 -
                         7560 * a2 * b2 * c * d * e4 * f3 * n2 - 3150 * a * b2 * d2 * e4 * f3 * n2 -
                         6300 * a3 * b2 * d2 * e4 * f3 * n2 - 1260 * a * b3 * c * e4 * f4 * n2 -
                         420 * b3 * d * e4 * f4 * n2 - 2520 * a2 * b3 * d * e4 * f4 * n2 -
                         315 * a * b4 * e4 * f5 * n2 - 21 * c5 * e3 * g * n2 -
                         315 * a * c4 * d * e3 * g * n2 - 210 * c3 * d2 * e3 * g * n2 -
                         1260 * a2 * c3 * d2 * e3 * g * n2 - 1050 * a * c2 * d3 * e3 * g * n2 -
                         2100 * a3 * c2 * d3 * e3 * g * n2 - 1575 * a2 * c * d4 * e3 * g * n2 -
                         1575 * a4 * c * d4 * e3 * g * n2 - 735 * a3 * d5 * e3 * g * n2 -
                         441 * a5 * d5 * e3 * g * n2 - 525 * b * c4 * e3 * f * g * n2 -
                         6300 * a * b * c3 * d * e3 * f * g * n2 -
                         3150 * b * c2 * d2 * e3 * f * g * n2 -
                         18900 * a2 * b * c2 * d2 * e3 * f * g * n2 -
                         10500 * a * b * c * d3 * e3 * f * g * n2 -
                         21000 * a3 * b * c * d3 * e3 * f * g * n2 -
                         7875 * a2 * b * d4 * e3 * f * g * n2 -
                         7875 * a4 * b * d4 * e3 * f * g * n2 - 1890 * b2 * c3 * e3 * f2 * g * n2 -
                         17010 * a * b2 * c2 * d * e3 * f2 * g * n2 -
                         5670 * b2 * c * d2 * e3 * f2 * g * n2 -
                         34020 * a2 * b2 * c * d2 * e3 * f2 * g * n2 -
                         9450 * a * b2 * d3 * e3 * f2 * g * n2 -
                         18900 * a3 * b2 * d3 * e3 * f2 * g * n2 -
                         2730 * b3 * c2 * e3 * f3 * g * n2 -
                         16380 * a * b3 * c * d * e3 * f3 * g * n2 -
                         2730 * b3 * d2 * e3 * f3 * g * n2 -
                         16380 * a2 * b3 * d2 * e3 * f3 * g * n2 -
                         1785 * b4 * c * e3 * f4 * g * n2 - 5355 * a * b4 * d * e3 * f4 * g * n2 -
                         441 * b5 * e3 * f5 * g * n2 - 315 * b * c4 * d * e2 * g2 * n2 -
                         1890 * a * b * c3 * d2 * e2 * g2 * n2 - 630 * b * c2 * d3 * e2 * g2 * n2 -
                         3780 * a2 * b * c2 * d3 * e2 * g2 * n2 -
                         1575 * a * b * c * d4 * e2 * g2 * n2 -
                         3150 * a3 * b * c * d4 * e2 * g2 * n2 - 945 * a2 * b * d5 * e2 * g2 * n2 -
                         945 * a4 * b * d5 * e2 * g2 * n2 - 3780 * b2 * c3 * d * e2 * f * g2 * n2 -
                         17010 * a * b2 * c2 * d2 * e2 * f * g2 * n2 -
                         3780 * b2 * c * d3 * e2 * f * g2 * n2 -
                         22680 * a2 * b2 * c * d3 * e2 * f * g2 * n2 -
                         4725 * a * b2 * d4 * e2 * f * g2 * n2 -
                         9450 * a3 * b2 * d4 * e2 * f * g2 * n2 -
                         9450 * b3 * c2 * d * e2 * f2 * g2 * n2 -
                         28350 * a * b3 * c * d2 * e2 * f2 * g2 * n2 -
                         3150 * b3 * d3 * e2 * f2 * g2 * n2 -
                         18900 * a2 * b3 * d3 * e2 * f2 * g2 * n2 -
                         8820 * b4 * c * d * e2 * f3 * g2 * n2 -
                         13230 * a * b4 * d2 * e2 * f3 * g2 * n2 -
                         2835 * b5 * d * e2 * f4 * g2 * n2 - 630 * b2 * c3 * d2 * e * g3 * n2 -
                         1890 * a * b2 * c2 * d3 * e * g3 * n2 - 315 * b2 * c * d4 * e * g3 * n2 -
                         1890 * a2 * b2 * c * d4 * e * g3 * n2 - 315 * a * b2 * d5 * e * g3 * n2 -
                         630 * a3 * b2 * d5 * e * g3 * n2 - 4410 * b3 * c2 * d2 * e * f * g3 * n2 -
                         8820 * a * b3 * c * d3 * e * f * g3 * n2 -
                         735 * b3 * d4 * e * f * g3 * n2 - 4410 * a2 * b3 * d4 * e * f * g3 * n2 -
                         6930 * b4 * c * d2 * e * f2 * g3 * n2 -
                         6930 * a * b4 * d3 * e * f2 * g3 * n2 - 3150 * b5 * d2 * e * f3 * g3 * n2 -
                         210 * b3 * c2 * d3 * g4 * n2 - 315 * a * b3 * c * d4 * g4 * n2 -
                         21 * b3 * d5 * g4 * n2 - 126 * a2 * b3 * d5 * g4 * n2 -
                         840 * b4 * c * d3 * f * g4 * n2 - 630 * a * b4 * d4 * f * g4 * n2 -
                         630 * b5 * d3 * f2 * g4 * n2 + 630 * a * c2 * e4 * f * n4 +
                         315 * a3 * c2 * e4 * f * n4 + 630 * c * d * e4 * f * n4 +
                         2520 * a2 * c * d * e4 * f * n4 + 630 * a4 * c * d * e4 * f * n4 +
                         1575 * a * d2 * e4 * f * n4 + 2100 * a3 * d2 * e4 * f * n4 +
                         315 * a5 * d2 * e4 * f * n4 + 1260 * a * b * c * e4 * f2 * n4 +
                         630 * a3 * b * c * e4 * f2 * n4 + 630 * b * d * e4 * f2 * n4 +
                         2520 * a2 * b * d * e4 * f2 * n4 + 630 * a4 * b * d * e4 * f2 * n4 +
                         630 * a * b2 * e4 * f3 * n4 + 315 * a3 * b2 * e4 * f3 * n4 +
                         70 * c3 * e3 * g * n4 + 105 * a2 * c3 * e3 * g * n4 +
                         630 * a * c2 * d * e3 * g * n4 + 315 * a3 * c2 * d * e3 * g * n4 +
                         315 * c * d2 * e3 * g * n4 + 1260 * a2 * c * d2 * e3 * g * n4 +
                         315 * a4 * c * d2 * e3 * g * n4 + 525 * a * d3 * e3 * g * n4 +
                         700 * a3 * d3 * e3 * g * n4 + 105 * a5 * d3 * e3 * g * n4 +
                         1050 * b * c2 * e3 * f * g * n4 + 1575 * a2 * b * c2 * e3 * f * g * n4 +
                         6300 * a * b * c * d * e3 * f * g * n4 +
                         3150 * a3 * b * c * d * e3 * f * g * n4 + 1575 * b * d2 * e3 * f * g * n4 +
                         6300 * a2 * b * d2 * e3 * f * g * n4 +
                         1575 * a4 * b * d2 * e3 * f * g * n4 + 1890 * b2 * c * e3 * f2 * g * n4 +
                         2835 * a2 * b2 * c * e3 * f2 * g * n4 +
                         5670 * a * b2 * d * e3 * f2 * g * n4 +
                         2835 * a3 * b2 * d * e3 * f2 * g * n4 + 910 * b3 * e3 * f3 * g * n4 +
                         1365 * a2 * b3 * e3 * f3 * g * n4 + 315 * a * b * c3 * e2 * g2 * n4 +
                         630 * b * c2 * d * e2 * g2 * n4 + 945 * a2 * b * c2 * d * e2 * g2 * n4 +
                         1890 * a * b * c * d2 * e2 * g2 * n4 +
                         945 * a3 * b * c * d2 * e2 * g2 * n4 + 315 * b * d3 * e2 * g2 * n4 +
                         1260 * a2 * b * d3 * e2 * g2 * n4 + 315 * a4 * b * d3 * e2 * g2 * n4 +
                         2835 * a * b2 * c2 * e2 * f * g2 * n4 +
                         3780 * b2 * c * d * e2 * f * g2 * n4 +
                         5670 * a2 * b2 * c * d * e2 * f * g2 * n4 +
                         5670 * a * b2 * d2 * e2 * f * g2 * n4 +
                         2835 * a3 * b2 * d2 * e2 * f * g2 * n4 +
                         4725 * a * b3 * c * e2 * f2 * g2 * n4 + 3150 * b3 * d * e2 * f2 * g2 * n4 +
                         4725 * a2 * b3 * d * e2 * f2 * g2 * n4 +
                         2205 * a * b4 * e2 * f3 * g2 * n4 + 315 * b2 * c3 * e * g3 * n4 +
                         945 * a * b2 * c2 * d * e * g3 * n4 + 630 * b2 * c * d2 * e * g3 * n4 +
                         945 * a2 * b2 * c * d2 * e * g3 * n4 + 630 * a * b2 * d3 * e * g3 * n4 +
                         315 * a3 * b2 * d3 * e * g3 * n4 + 2205 * b3 * c2 * e * f * g3 * n4 +
                         4410 * a * b3 * c * d * e * f * g3 * n4 +
                         1470 * b3 * d2 * e * f * g3 * n4 + 2205 * a2 * b3 * d2 * e * f * g3 * n4 +
                         3465 * b4 * c * e * f2 * g3 * n4 + 3465 * a * b4 * d * e * f2 * g3 * n4 +
                         1575 * b5 * e * f3 * g3 * n4 + 315 * b3 * c2 * d * g4 * n4 +
                         315 * a * b3 * c * d2 * g4 * n4 + 70 * b3 * d3 * g4 * n4 +
                         105 * a2 * b3 * d3 * g4 * n4 + 1260 * b4 * c * d * f * g4 * n4 +
                         630 * a * b4 * d2 * f * g4 * n4 + 945 * b5 * d * f2 * g4 * n4 -
                         315 * a * e4 * f * n6 - 105 * a3 * e4 * f * n6 - 105 * c * e3 * g * n6 -
                         105 * a2 * c * e3 * g * n6 - 315 * a * d * e3 * g * n6 -
                         105 * a3 * d * e3 * g * n6 - 525 * b * e3 * f * g * n6 -
                         525 * a2 * b * e3 * f * g * n6 - 315 * a * b * c * e2 * g2 * n6 -
                         315 * b * d * e2 * g2 * n6 - 315 * a2 * b * d * e2 * g2 * n6 -
                         945 * a * b2 * e2 * f * g2 * n6 - 315 * b2 * c * e * g3 * n6 -
                         315 * a * b2 * d * e * g3 * n6 - 735 * b3 * e * f * g3 * n6 -
                         105 * b3 * d * g4 * n6) +
                 (14 *
                     (24 * c5 * e5 * f + 600 * a * c4 * d * e5 * f + 3600 * a2 * c3 * d2 * e5 * f +
                         8400 * a3 * c2 * d3 * e5 * f + 8400 * a4 * c * d4 * e5 * f +
                         3024 * a5 * d5 * e5 * f + 120 * b * c4 * e5 * f2 +
                         2400 * a * b * c3 * d * e5 * f2 + 10800 * a2 * b * c2 * d2 * e5 * f2 +
                         16800 * a3 * b * c * d3 * e5 * f2 + 8400 * a4 * b * d4 * e5 * f2 +
                         240 * b2 * c3 * e5 * f3 + 3600 * a * b2 * c2 * d * e5 * f3 +
                         10800 * a2 * b2 * c * d2 * e5 * f3 + 8400 * a3 * b2 * d3 * e5 * f3 +
                         240 * b3 * c2 * e5 * f4 + 2400 * a * b3 * c * d * e5 * f4 +
                         3600 * a2 * b3 * d2 * e5 * f4 + 120 * b4 * c * e5 * f5 +
                         600 * a * b4 * d * e5 * f5 + 24 * b5 * e5 * f6 + 30 * c5 * d * e4 * g +
                         375 * a * c4 * d2 * e4 * g + 1500 * a2 * c3 * d3 * e4 * g +
                         2625 * a3 * c2 * d4 * e4 * g + 2100 * a4 * c * d5 * e4 * g +
                         630 * a5 * d6 * e4 * g + 750 * b * c4 * d * e4 * f * g +
                         7500 * a * b * c3 * d2 * e4 * f * g +
                         22500 * a2 * b * c2 * d3 * e4 * f * g +
                         26250 * a3 * b * c * d4 * e4 * f * g + 10500 * a4 * b * d5 * e4 * f * g +
                         2700 * b2 * c3 * d * e4 * f2 * g + 20250 * a * b2 * c2 * d2 * e4 * f2 * g +
                         40500 * a2 * b2 * c * d3 * e4 * f2 * g +
                         23625 * a3 * b2 * d4 * e4 * f2 * g + 3900 * b3 * c2 * d * e4 * f3 * g +
                         19500 * a * b3 * c * d2 * e4 * f3 * g +
                         19500 * a2 * b3 * d3 * e4 * f3 * g + 2550 * b4 * c * d * e4 * f4 * g +
                         6375 * a * b4 * d2 * e4 * f4 * g + 630 * b5 * d * e4 * f5 * g +
                         300 * b * c4 * d2 * e3 * g2 + 2000 * a * b * c3 * d3 * e3 * g2 +
                         4500 * a2 * b * c2 * d4 * e3 * g2 + 4200 * a3 * b * c * d5 * e3 * g2 +
                         1400 * a4 * b * d6 * e3 * g2 + 3600 * b2 * c3 * d2 * e3 * f * g2 +
                         18000 * a * b2 * c2 * d3 * e3 * f * g2 +
                         27000 * a2 * b2 * c * d4 * e3 * f * g2 +
                         12600 * a3 * b2 * d5 * e3 * f * g2 + 9000 * b3 * c2 * d2 * e3 * f2 * g2 +
                         30000 * a * b3 * c * d3 * e3 * f2 * g2 +
                         22500 * a2 * b3 * d4 * e3 * f2 * g2 + 8400 * b4 * c * d2 * e3 * f3 * g2 +
                         14000 * a * b4 * d3 * e3 * f3 * g2 + 2700 * b5 * d2 * e3 * f4 * g2 +
                         600 * b2 * c3 * d3 * e2 * g3 + 2250 * a * b2 * c2 * d4 * e2 * g3 +
                         2700 * a2 * b2 * c * d5 * e2 * g3 + 1050 * a3 * b2 * d6 * e2 * g3 +
                         4200 * b3 * c2 * d3 * e2 * f * g3 + 10500 * a * b3 * c * d4 * e2 * f * g3 +
                         6300 * a2 * b3 * d5 * e2 * f * g3 + 6600 * b4 * c * d3 * e2 * f2 * g3 +
                         8250 * a * b4 * d4 * e2 * f2 * g3 + 3000 * b5 * d3 * e2 * f3 * g3 +
                         300 * b3 * c2 * d4 * e * g4 + 600 * a * b3 * c * d5 * e * g4 +
                         300 * a2 * b3 * d6 * e * g4 + 1200 * b4 * c * d4 * e * f * g4 +
                         1200 * a * b4 * d5 * e * f * g4 + 900 * b5 * d4 * e * f2 * g4 +
                         30 * b4 * c * d5 * g5 + 25 * a * b4 * d6 * g5 + 54 * b5 * d5 * f * g5 -
                         80 * c3 * e5 * f * n2 - 480 * a2 * c3 * e5 * f * n2 -
                         1200 * a * c2 * d * e5 * f * n2 - 2400 * a3 * c2 * d * e5 * f * n2 -
                         3600 * a2 * c * d2 * e5 * f * n2 - 3600 * a4 * c * d2 * e5 * f * n2 -
                         2800 * a3 * d3 * e5 * f * n2 - 1680 * a5 * d3 * e5 * f * n2 -
                         240 * b * c2 * e5 * f2 * n2 - 1440 * a2 * b * c2 * e5 * f2 * n2 -
                         2400 * a * b * c * d * e5 * f2 * n2 -
                         4800 * a3 * b * c * d * e5 * f2 * n2 - 3600 * a2 * b * d2 * e5 * f2 * n2 -
                         3600 * a4 * b * d2 * e5 * f2 * n2 - 240 * b2 * c * e5 * f3 * n2 -
                         1440 * a2 * b2 * c * e5 * f3 * n2 - 1200 * a * b2 * d * e5 * f3 * n2 -
                         2400 * a3 * b2 * d * e5 * f3 * n2 - 80 * b3 * e5 * f4 * n2 -
                         480 * a2 * b3 * e5 * f4 * n2 - 75 * a * c4 * e4 * g * n2 -
                         100 * c3 * d * e4 * g * n2 - 600 * a2 * c3 * d * e4 * g * n2 -
                         750 * a * c2 * d2 * e4 * g * n2 - 1500 * a3 * c2 * d2 * e4 * g * n2 -
                         1500 * a2 * c * d3 * e4 * g * n2 - 1500 * a4 * c * d3 * e4 * g * n2 -
                         875 * a3 * d4 * e4 * g * n2 - 525 * a5 * d4 * e4 * g * n2 -
                         1500 * a * b * c3 * e4 * f * g * n2 - 1500 * b * c2 * d * e4 * f * g * n2 -
                         9000 * a2 * b * c2 * d * e4 * f * g * n2 -
                         7500 * a * b * c * d2 * e4 * f * g * n2 -
                         15000 * a3 * b * c * d2 * e4 * f * g * n2 -
                         7500 * a2 * b * d3 * e4 * f * g * n2 -
                         7500 * a4 * b * d3 * e4 * f * g * n2 -
                         4050 * a * b2 * c2 * e4 * f2 * g * n2 -
                         2700 * b2 * c * d * e4 * f2 * g * n2 -
                         16200 * a2 * b2 * c * d * e4 * f2 * g * n2 -
                         6750 * a * b2 * d2 * e4 * f2 * g * n2 -
                         13500 * a3 * b2 * d2 * e4 * f2 * g * n2 -
                         3900 * a * b3 * c * e4 * f3 * g * n2 - 1300 * b3 * d * e4 * f3 * g * n2 -
                         7800 * a2 * b3 * d * e4 * f3 * g * n2 - 1275 * a * b4 * e4 * f4 * g * n2 -
                         100 * b * c4 * e3 * g2 * n2 - 1200 * a * b * c3 * d * e3 * g2 * n2 -
                         600 * b * c2 * d2 * e3 * g2 * n2 - 3600 * a2 * b * c2 * d2 * e3 * g2 * n2 -
                         2000 * a * b * c * d3 * e3 * g2 * n2 -
                         4000 * a3 * b * c * d3 * e3 * g2 * n2 - 1500 * a2 * b * d4 * e3 * g2 * n2 -
                         1500 * a4 * b * d4 * e3 * g2 * n2 - 1200 * b2 * c3 * e3 * f * g2 * n2 -
                         10800 * a * b2 * c2 * d * e3 * f * g2 * n2 -
                         3600 * b2 * c * d2 * e3 * f * g2 * n2 -
                         21600 * a2 * b2 * c * d2 * e3 * f * g2 * n2 -
                         6000 * a * b2 * d3 * e3 * f * g2 * n2 -
                         12000 * a3 * b2 * d3 * e3 * f * g2 * n2 -
                         3000 * b3 * c2 * e3 * f2 * g2 * n2 -
                         18000 * a * b3 * c * d * e3 * f2 * g2 * n2 -
                         3000 * b3 * d2 * e3 * f2 * g2 * n2 -
                         18000 * a2 * b3 * d2 * e3 * f2 * g2 * n2 -
                         2800 * b4 * c * e3 * f3 * g2 * n2 - 8400 * a * b4 * d * e3 * f3 * g2 * n2 -
                         900 * b5 * e3 * f4 * g2 * n2 - 600 * b2 * c3 * d * e2 * g3 * n2 -
                         2700 * a * b2 * c2 * d2 * e2 * g3 * n2 - 600 * b2 * c * d3 * e2 * g3 * n2 -
                         3600 * a2 * b2 * c * d3 * e2 * g3 * n2 - 750 * a * b2 * d4 * e2 * g3 * n2 -
                         1500 * a3 * b2 * d4 * e2 * g3 * n2 -
                         4200 * b3 * c2 * d * e2 * f * g3 * n2 -
                         12600 * a * b3 * c * d2 * e2 * f * g3 * n2 -
                         1400 * b3 * d3 * e2 * f * g3 * n2 -
                         8400 * a2 * b3 * d3 * e2 * f * g3 * n2 -
                         6600 * b4 * c * d * e2 * f2 * g3 * n2 -
                         9900 * a * b4 * d2 * e2 * f2 * g3 * n2 -
                         3000 * b5 * d * e2 * f3 * g3 * n2 - 600 * b3 * c2 * d2 * e * g4 * n2 -
                         1200 * a * b3 * c * d3 * e * g4 * n2 - 100 * b3 * d4 * e * g4 * n2 -
                         600 * a2 * b3 * d4 * e * g4 * n2 - 2400 * b4 * c * d2 * e * f * g4 * n2 -
                         2400 * a * b4 * d3 * e * f * g4 * n2 - 1800 * b5 * d2 * e * f2 * g4 * n2 -
                         100 * b4 * c * d3 * g5 * n2 - 75 * a * b4 * d4 * g5 * n2 -
                         180 * b5 * d3 * f * g5 * n2 + 120 * c * e5 * f * n4 +
                         480 * a2 * c * e5 * f * n4 + 120 * a4 * c * e5 * f * n4 +
                         600 * a * d * e5 * f * n4 + 800 * a3 * d * e5 * f * n4 +
                         120 * a5 * d * e5 * f * n4 + 120 * b * e5 * f2 * n4 +
                         480 * a2 * b * e5 * f2 * n4 + 120 * a4 * b * e5 * f2 * n4 +
                         150 * a * c2 * e4 * g * n4 + 75 * a3 * c2 * e4 * g * n4 +
                         150 * c * d * e4 * g * n4 + 600 * a2 * c * d * e4 * g * n4 +
                         150 * a4 * c * d * e4 * g * n4 + 375 * a * d2 * e4 * g * n4 +
                         500 * a3 * d2 * e4 * g * n4 + 75 * a5 * d2 * e4 * g * n4 +
                         1500 * a * b * c * e4 * f * g * n4 + 750 * a3 * b * c * e4 * f * g * n4 +
                         750 * b * d * e4 * f * g * n4 + 3000 * a2 * b * d * e4 * f * g * n4 +
                         750 * a4 * b * d * e4 * f * g * n4 + 1350 * a * b2 * e4 * f2 * g * n4 +
                         675 * a3 * b2 * e4 * f2 * g * n4 + 200 * b * c2 * e3 * g2 * n4 +
                         300 * a2 * b * c2 * e3 * g2 * n4 + 1200 * a * b * c * d * e3 * g2 * n4 +
                         600 * a3 * b * c * d * e3 * g2 * n4 + 300 * b * d2 * e3 * g2 * n4 +
                         1200 * a2 * b * d2 * e3 * g2 * n4 + 300 * a4 * b * d2 * e3 * g2 * n4 +
                         1200 * b2 * c * e3 * f * g2 * n4 + 1800 * a2 * b2 * c * e3 * f * g2 * n4 +
                         3600 * a * b2 * d * e3 * f * g2 * n4 +
                         1800 * a3 * b2 * d * e3 * f * g2 * n4 + 1000 * b3 * e3 * f2 * g2 * n4 +
                         1500 * a2 * b3 * e3 * f2 * g2 * n4 + 450 * a * b2 * c2 * e2 * g3 * n4 +
                         600 * b2 * c * d * e2 * g3 * n4 + 900 * a2 * b2 * c * d * e2 * g3 * n4 +
                         900 * a * b2 * d2 * e2 * g3 * n4 + 450 * a3 * b2 * d2 * e2 * g3 * n4 +
                         2100 * a * b3 * c * e2 * f * g3 * n4 + 1400 * b3 * d * e2 * f * g3 * n4 +
                         2100 * a2 * b3 * d * e2 * f * g3 * n4 + 1650 * a * b4 * e2 * f2 * g3 * n4 +
                         300 * b3 * c2 * e * g4 * n4 + 600 * a * b3 * c * d * e * g4 * n4 +
                         200 * b3 * d2 * e * g4 * n4 + 300 * a2 * b3 * d2 * e * g4 * n4 +
                         1200 * b4 * c * e * f * g4 * n4 + 1200 * a * b4 * d * e * f * g4 * n4 +
                         900 * b5 * e * f2 * g4 * n4 + 150 * b4 * c * d * g5 * n4 +
                         75 * a * b4 * d2 * g5 * n4 + 270 * b5 * d * f * g5 * n4 -
                         75 * a * e4 * g * n6 - 25 * a3 * e4 * g * n6 - 100 * b * e3 * g2 * n6 -
                         100 * a2 * b * e3 * g2 * n6 - 150 * a * b2 * e2 * g3 * n6 -
                         100 * b3 * e * g4 * n6)) /
                     3 +
                 (14 *
                     (150 * a * c4 * e6 * f + 1800 * a2 * c3 * d * e6 * f +
                         6300 * a3 * c2 * d2 * e6 * f + 8400 * a4 * c * d3 * e6 * f +
                         3780 * a5 * d4 * e6 * f + 600 * a * b * c3 * e6 * f2 +
                         5400 * a2 * b * c2 * d * e6 * f2 + 12600 * a3 * b * c * d2 * e6 * f2 +
                         8400 * a4 * b * d3 * e6 * f2 + 900 * a * b2 * c2 * e6 * f3 +
                         5400 * a2 * b2 * c * d * e6 * f3 + 6300 * a3 * b2 * d2 * e6 * f3 +
                         600 * a * b3 * c * e6 * f4 + 1800 * a2 * b3 * d * e6 * f4 +
                         150 * a * b4 * e6 * f5 + 9 * c5 * e5 * g + 225 * a * c4 * d * e5 * g +
                         1350 * a2 * c3 * d2 * e5 * g + 3150 * a3 * c2 * d3 * e5 * g +
                         3150 * a4 * c * d4 * e5 * g + 1134 * a5 * d5 * e5 * g +
                         225 * b * c4 * e5 * f * g + 4500 * a * b * c3 * d * e5 * f * g +
                         20250 * a2 * b * c2 * d2 * e5 * f * g +
                         31500 * a3 * b * c * d3 * e5 * f * g + 15750 * a4 * b * d4 * e5 * f * g +
                         810 * b2 * c3 * e5 * f2 * g + 12150 * a * b2 * c2 * d * e5 * f2 * g +
                         36450 * a2 * b2 * c * d2 * e5 * f2 * g +
                         28350 * a3 * b2 * d3 * e5 * f2 * g + 1170 * b3 * c2 * e5 * f3 * g +
                         11700 * a * b3 * c * d * e5 * f3 * g + 17550 * a2 * b3 * d2 * e5 * f3 * g +
                         765 * b4 * c * e5 * f4 * g + 3825 * a * b4 * d * e5 * f4 * g +
                         189 * b5 * e5 * f5 * g + 225 * b * c4 * d * e4 * g2 +
                         2250 * a * b * c3 * d2 * e4 * g2 + 6750 * a2 * b * c2 * d3 * e4 * g2 +
                         7875 * a3 * b * c * d4 * e4 * g2 + 3150 * a4 * b * d5 * e4 * g2 +
                         2700 * b2 * c3 * d * e4 * f * g2 + 20250 * a * b2 * c2 * d2 * e4 * f * g2 +
                         40500 * a2 * b2 * c * d3 * e4 * f * g2 +
                         23625 * a3 * b2 * d4 * e4 * f * g2 + 6750 * b3 * c2 * d * e4 * f2 * g2 +
                         33750 * a * b3 * c * d2 * e4 * f2 * g2 +
                         33750 * a2 * b3 * d3 * e4 * f2 * g2 + 6300 * b4 * c * d * e4 * f3 * g2 +
                         15750 * a * b4 * d2 * e4 * f3 * g2 + 2025 * b5 * d * e4 * f4 * g2 +
                         900 * b2 * c3 * d2 * e3 * g3 + 4500 * a * b2 * c2 * d3 * e3 * g3 +
                         6750 * a2 * b2 * c * d4 * e3 * g3 + 3150 * a3 * b2 * d5 * e3 * g3 +
                         6300 * b3 * c2 * d2 * e3 * f * g3 + 21000 * a * b3 * c * d3 * e3 * f * g3 +
                         15750 * a2 * b3 * d4 * e3 * f * g3 + 9900 * b4 * c * d2 * e3 * f2 * g3 +
                         16500 * a * b4 * d3 * e3 * f2 * g3 + 4500 * b5 * d2 * e3 * f3 * g3 +
                         900 * b3 * c2 * d3 * e2 * g4 + 2250 * a * b3 * c * d4 * e2 * g4 +
                         1350 * a2 * b3 * d5 * e2 * g4 + 3600 * b4 * c * d3 * e2 * f * g4 +
                         4500 * a * b4 * d4 * e2 * f * g4 + 2700 * b5 * d3 * e2 * f2 * g4 +
                         225 * b4 * c * d4 * e * g5 + 225 * a * b4 * d5 * e * g5 +
                         405 * b5 * d4 * e * f * g5 + 9 * b5 * d5 * g6 -
                         300 * a * c2 * e6 * f * n2 - 600 * a3 * c2 * e6 * f * n2 -
                         1800 * a2 * c * d * e6 * f * n2 - 1800 * a4 * c * d * e6 * f * n2 -
                         2100 * a3 * d2 * e6 * f * n2 - 1260 * a5 * d2 * e6 * f * n2 -
                         600 * a * b * c * e6 * f2 * n2 - 1200 * a3 * b * c * e6 * f2 * n2 -
                         1800 * a2 * b * d * e6 * f2 * n2 - 1800 * a4 * b * d * e6 * f2 * n2 -
                         300 * a * b2 * e6 * f3 * n2 - 600 * a3 * b2 * e6 * f3 * n2 -
                         30 * c3 * e5 * g * n2 - 180 * a2 * c3 * e5 * g * n2 -
                         450 * a * c2 * d * e5 * g * n2 - 900 * a3 * c2 * d * e5 * g * n2 -
                         1350 * a2 * c * d2 * e5 * g * n2 - 1350 * a4 * c * d2 * e5 * g * n2 -
                         1050 * a3 * d3 * e5 * g * n2 - 630 * a5 * d3 * e5 * g * n2 -
                         450 * b * c2 * e5 * f * g * n2 - 2700 * a2 * b * c2 * e5 * f * g * n2 -
                         4500 * a * b * c * d * e5 * f * g * n2 -
                         9000 * a3 * b * c * d * e5 * f * g * n2 -
                         6750 * a2 * b * d2 * e5 * f * g * n2 -
                         6750 * a4 * b * d2 * e5 * f * g * n2 - 810 * b2 * c * e5 * f2 * g * n2 -
                         4860 * a2 * b2 * c * e5 * f2 * g * n2 -
                         4050 * a * b2 * d * e5 * f2 * g * n2 -
                         8100 * a3 * b2 * d * e5 * f2 * g * n2 - 390 * b3 * e5 * f3 * g * n2 -
                         2340 * a2 * b3 * e5 * f3 * g * n2 - 450 * a * b * c3 * e4 * g2 * n2 -
                         450 * b * c2 * d * e4 * g2 * n2 - 2700 * a2 * b * c2 * d * e4 * g2 * n2 -
                         2250 * a * b * c * d2 * e4 * g2 * n2 -
                         4500 * a3 * b * c * d2 * e4 * g2 * n2 - 2250 * a2 * b * d3 * e4 * g2 * n2 -
                         2250 * a4 * b * d3 * e4 * g2 * n2 - 4050 * a * b2 * c2 * e4 * f * g2 * n2 -
                         2700 * b2 * c * d * e4 * f * g2 * n2 -
                         16200 * a2 * b2 * c * d * e4 * f * g2 * n2 -
                         6750 * a * b2 * d2 * e4 * f * g2 * n2 -
                         13500 * a3 * b2 * d2 * e4 * f * g2 * n2 -
                         6750 * a * b3 * c * e4 * f2 * g2 * n2 - 2250 * b3 * d * e4 * f2 * g2 * n2 -
                         13500 * a2 * b3 * d * e4 * f2 * g2 * n2 -
                         3150 * a * b4 * e4 * f3 * g2 * n2 - 300 * b2 * c3 * e3 * g3 * n2 -
                         2700 * a * b2 * c2 * d * e3 * g3 * n2 - 900 * b2 * c * d2 * e3 * g3 * n2 -
                         5400 * a2 * b2 * c * d2 * e3 * g3 * n2 -
                         1500 * a * b2 * d3 * e3 * g3 * n2 - 3000 * a3 * b2 * d3 * e3 * g3 * n2 -
                         2100 * b3 * c2 * e3 * f * g3 * n2 -
                         12600 * a * b3 * c * d * e3 * f * g3 * n2 -
                         2100 * b3 * d2 * e3 * f * g3 * n2 -
                         12600 * a2 * b3 * d2 * e3 * f * g3 * n2 -
                         3300 * b4 * c * e3 * f2 * g3 * n2 - 9900 * a * b4 * d * e3 * f2 * g3 * n2 -
                         1500 * b5 * e3 * f3 * g3 * n2 - 900 * b3 * c2 * d * e2 * g4 * n2 -
                         2700 * a * b3 * c * d2 * e2 * g4 * n2 - 300 * b3 * d3 * e2 * g4 * n2 -
                         1800 * a2 * b3 * d3 * e2 * g4 * n2 - 3600 * b4 * c * d * e2 * f * g4 * n2 -
                         5400 * a * b4 * d2 * e2 * f * g4 * n2 - 2700 * b5 * d * e2 * f2 * g4 * n2 -
                         450 * b4 * c * d2 * e * g5 * n2 - 450 * a * b4 * d3 * e * g5 * n2 -
                         810 * b5 * d2 * e * f * g5 * n2 - 30 * b5 * d3 * g6 * n2 +
                         150 * a * e6 * f * n4 + 200 * a3 * e6 * f * n4 + 30 * a5 * e6 * f * n4 +
                         45 * c * e5 * g * n4 + 180 * a2 * c * e5 * g * n4 +
                         45 * a4 * c * e5 * g * n4 + 225 * a * d * e5 * g * n4 +
                         300 * a3 * d * e5 * g * n4 + 45 * a5 * d * e5 * g * n4 +
                         225 * b * e5 * f * g * n4 + 900 * a2 * b * e5 * f * g * n4 +
                         225 * a4 * b * e5 * f * g * n4 + 450 * a * b * c * e4 * g2 * n4 +
                         225 * a3 * b * c * e4 * g2 * n4 + 225 * b * d * e4 * g2 * n4 +
                         900 * a2 * b * d * e4 * g2 * n4 + 225 * a4 * b * d * e4 * g2 * n4 +
                         1350 * a * b2 * e4 * f * g2 * n4 + 675 * a3 * b2 * e4 * f * g2 * n4 +
                         300 * b2 * c * e3 * g3 * n4 + 450 * a2 * b2 * c * e3 * g3 * n4 +
                         900 * a * b2 * d * e3 * g3 * n4 + 450 * a3 * b2 * d * e3 * g3 * n4 +
                         700 * b3 * e3 * f * g3 * n4 + 1050 * a2 * b3 * e3 * f * g3 * n4 +
                         450 * a * b3 * c * e2 * g4 * n4 + 300 * b3 * d * e2 * g4 * n4 +
                         450 * a2 * b3 * d * e2 * g4 * n4 + 900 * a * b4 * e2 * f * g4 * n4 +
                         225 * b4 * c * e * g5 * n4 + 225 * a * b4 * d * e * g5 * n4 +
                         405 * b5 * e * f * g5 * n4 + 45 * b5 * d * g6 * n4)) /
                     5 +
                 (10 *
                     (720 * a2 * c3 * e7 * f + 5040 * a3 * c2 * d * e7 * f +
                         10080 * a4 * c * d2 * e7 * f + 6048 * a5 * d3 * e7 * f +
                         2160 * a2 * b * c2 * e7 * f2 + 10080 * a3 * b * c * d * e7 * f2 +
                         10080 * a4 * b * d2 * e7 * f2 + 2160 * a2 * b2 * c * e7 * f3 +
                         5040 * a3 * b2 * d * e7 * f3 + 720 * a2 * b3 * e7 * f4 +
                         105 * a * c4 * e6 * g + 1260 * a2 * c3 * d * e6 * g +
                         4410 * a3 * c2 * d2 * e6 * g + 5880 * a4 * c * d3 * e6 * g +
                         2646 * a5 * d4 * e6 * g + 2100 * a * b * c3 * e6 * f * g +
                         18900 * a2 * b * c2 * d * e6 * f * g +
                         44100 * a3 * b * c * d2 * e6 * f * g + 29400 * a4 * b * d3 * e6 * f * g +
                         5670 * a * b2 * c2 * e6 * f2 * g + 34020 * a2 * b2 * c * d * e6 * f2 * g +
                         39690 * a3 * b2 * d2 * e6 * f2 * g + 5460 * a * b3 * c * e6 * f3 * g +
                         16380 * a2 * b3 * d * e6 * f3 * g + 1785 * a * b4 * e6 * f4 * g +
                         126 * b * c4 * e5 * g2 + 2520 * a * b * c3 * d * e5 * g2 +
                         11340 * a2 * b * c2 * d2 * e5 * g2 + 17640 * a3 * b * c * d3 * e5 * g2 +
                         8820 * a4 * b * d4 * e5 * g2 + 1512 * b2 * c3 * e5 * f * g2 +
                         22680 * a * b2 * c2 * d * e5 * f * g2 +
                         68040 * a2 * b2 * c * d2 * e5 * f * g2 +
                         52920 * a3 * b2 * d3 * e5 * f * g2 + 3780 * b3 * c2 * e5 * f2 * g2 +
                         37800 * a * b3 * c * d * e5 * f2 * g2 +
                         56700 * a2 * b3 * d2 * e5 * f2 * g2 + 3528 * b4 * c * e5 * f3 * g2 +
                         17640 * a * b4 * d * e5 * f3 * g2 + 1134 * b5 * e5 * f4 * g2 +
                         1260 * b2 * c3 * d * e4 * g3 + 9450 * a * b2 * c2 * d2 * e4 * g3 +
                         18900 * a2 * b2 * c * d3 * e4 * g3 + 11025 * a3 * b2 * d4 * e4 * g3 +
                         8820 * b3 * c2 * d * e4 * f * g3 + 44100 * a * b3 * c * d2 * e4 * f * g3 +
                         44100 * a2 * b3 * d3 * e4 * f * g3 + 13860 * b4 * c * d * e4 * f2 * g3 +
                         34650 * a * b4 * d2 * e4 * f2 * g3 + 6300 * b5 * d * e4 * f3 * g3 +
                         2520 * b3 * c2 * d2 * e3 * g4 + 8400 * a * b3 * c * d3 * e3 * g4 +
                         6300 * a2 * b3 * d4 * e3 * g4 + 10080 * b4 * c * d2 * e3 * f * g4 +
                         16800 * a * b4 * d3 * e3 * f * g4 + 7560 * b5 * d2 * e3 * f2 * g4 +
                         1260 * b4 * c * d3 * e2 * g5 + 1575 * a * b4 * d4 * e2 * g5 +
                         2268 * b5 * d3 * e2 * f * g5 + 126 * b5 * d4 * e * g6 -
                         720 * a2 * c * e7 * f * n2 - 720 * a4 * c * e7 * f * n2 -
                         1680 * a3 * d * e7 * f * n2 - 1008 * a5 * d * e7 * f * n2 -
                         720 * a2 * b * e7 * f2 * n2 - 720 * a4 * b * e7 * f2 * n2 -
                         210 * a * c2 * e6 * g * n2 - 420 * a3 * c2 * e6 * g * n2 -
                         1260 * a2 * c * d * e6 * g * n2 - 1260 * a4 * c * d * e6 * g * n2 -
                         1470 * a3 * d2 * e6 * g * n2 - 882 * a5 * d2 * e6 * g * n2 -
                         2100 * a * b * c * e6 * f * g * n2 - 4200 * a3 * b * c * e6 * f * g * n2 -
                         6300 * a2 * b * d * e6 * f * g * n2 - 6300 * a4 * b * d * e6 * f * g * n2 -
                         1890 * a * b2 * e6 * f2 * g * n2 - 3780 * a3 * b2 * e6 * f2 * g * n2 -
                         252 * b * c2 * e5 * g2 * n2 - 1512 * a2 * b * c2 * e5 * g2 * n2 -
                         2520 * a * b * c * d * e5 * g2 * n2 -
                         5040 * a3 * b * c * d * e5 * g2 * n2 - 3780 * a2 * b * d2 * e5 * g2 * n2 -
                         3780 * a4 * b * d2 * e5 * g2 * n2 - 1512 * b2 * c * e5 * f * g2 * n2 -
                         9072 * a2 * b2 * c * e5 * f * g2 * n2 -
                         7560 * a * b2 * d * e5 * f * g2 * n2 -
                         15120 * a3 * b2 * d * e5 * f * g2 * n2 - 1260 * b3 * e5 * f2 * g2 * n2 -
                         7560 * a2 * b3 * e5 * f2 * g2 * n2 - 1890 * a * b2 * c2 * e4 * g3 * n2 -
                         1260 * b2 * c * d * e4 * g3 * n2 - 7560 * a2 * b2 * c * d * e4 * g3 * n2 -
                         3150 * a * b2 * d2 * e4 * g3 * n2 - 6300 * a3 * b2 * d2 * e4 * g3 * n2 -
                         8820 * a * b3 * c * e4 * f * g3 * n2 - 2940 * b3 * d * e4 * f * g3 * n2 -
                         17640 * a2 * b3 * d * e4 * f * g3 * n2 -
                         6930 * a * b4 * e4 * f2 * g3 * n2 - 840 * b3 * c2 * e3 * g4 * n2 -
                         5040 * a * b3 * c * d * e3 * g4 * n2 - 840 * b3 * d2 * e3 * g4 * n2 -
                         5040 * a2 * b3 * d2 * e3 * g4 * n2 - 3360 * b4 * c * e3 * f * g4 * n2 -
                         10080 * a * b4 * d * e3 * f * g4 * n2 - 2520 * b5 * e3 * f2 * g4 * n2 -
                         1260 * b4 * c * d * e2 * g5 * n2 - 1890 * a * b4 * d2 * e2 * g5 * n2 -
                         2268 * b5 * d * e2 * f * g5 * n2 - 252 * b5 * d2 * e * g6 * n2 +
                         105 * a * e6 * g * n4 + 140 * a3 * e6 * g * n4 + 21 * a5 * e6 * g * n4 +
                         126 * b * e5 * g2 * n4 + 504 * a2 * b * e5 * g2 * n4 +
                         126 * a4 * b * e5 * g2 * n4 + 630 * a * b2 * e4 * g3 * n4 +
                         315 * a3 * b2 * e4 * g3 * n4 + 280 * b3 * e3 * g4 * n4 +
                         420 * a2 * b3 * e3 * g4 * n4 + 315 * a * b4 * e2 * g5 * n4 +
                         126 * b5 * e * g6 * n4)) /
                     11 +
                 5 * (105 * a3 * c2 * e8 * f + 420 * a4 * c * d * e8 * f + 378 * a5 * d2 * e8 * f +
                         210 * a3 * b * c * e8 * f2 + 420 * a4 * b * d * e8 * f2 +
                         105 * a3 * b2 * e8 * f3 + 30 * a2 * c3 * e7 * g +
                         210 * a3 * c2 * d * e7 * g + 420 * a4 * c * d2 * e7 * g +
                         252 * a5 * d3 * e7 * g + 450 * a2 * b * c2 * e7 * f * g +
                         2100 * a3 * b * c * d * e7 * f * g + 2100 * a4 * b * d2 * e7 * f * g +
                         810 * a2 * b2 * c * e7 * f2 * g + 1890 * a3 * b2 * d * e7 * f2 * g +
                         390 * a2 * b3 * e7 * f3 * g + 70 * a * b * c3 * e6 * g2 +
                         630 * a2 * b * c2 * d * e6 * g2 + 1470 * a3 * b * c * d2 * e6 * g2 +
                         980 * a4 * b * d3 * e6 * g2 + 630 * a * b2 * c2 * e6 * f * g2 +
                         3780 * a2 * b2 * c * d * e6 * f * g2 + 4410 * a3 * b2 * d2 * e6 * f * g2 +
                         1050 * a * b3 * c * e6 * f2 * g2 + 3150 * a2 * b3 * d * e6 * f2 * g2 +
                         490 * a * b4 * e6 * f3 * g2 + 42 * b2 * c3 * e5 * g3 +
                         630 * a * b2 * c2 * d * e5 * g3 + 1890 * a2 * b2 * c * d2 * e5 * g3 +
                         1470 * a3 * b2 * d3 * e5 * g3 + 294 * b3 * c2 * e5 * f * g3 +
                         2940 * a * b3 * c * d * e5 * f * g3 + 4410 * a2 * b3 * d2 * e5 * f * g3 +
                         462 * b4 * c * e5 * f2 * g3 + 2310 * a * b4 * d * e5 * f2 * g3 +
                         210 * b5 * e5 * f3 * g3 + 210 * b3 * c2 * d * e4 * g4 +
                         1050 * a * b3 * c * d2 * e4 * g4 + 1050 * a2 * b3 * d3 * e4 * g4 +
                         840 * b4 * c * d * e4 * f * g4 + 2100 * a * b4 * d2 * e4 * f * g4 +
                         630 * b5 * d * e4 * f2 * g4 + 210 * b4 * c * d2 * e3 * g5 +
                         350 * a * b4 * d3 * e3 * g5 + 378 * b5 * d2 * e3 * f * g5 +
                         42 * b5 * d3 * e2 * g6 - 35 * a3 * e8 * f * n2 - 21 * a5 * e8 * f * n2 -
                         30 * a2 * c * e7 * g * n2 - 30 * a4 * c * e7 * g * n2 -
                         70 * a3 * d * e7 * g * n2 - 42 * a5 * d * e7 * g * n2 -
                         150 * a2 * b * e7 * f * g * n2 - 150 * a4 * b * e7 * f * g * n2 -
                         70 * a * b * c * e6 * g2 * n2 - 140 * a3 * b * c * e6 * g2 * n2 -
                         210 * a2 * b * d * e6 * g2 * n2 - 210 * a4 * b * d * e6 * g2 * n2 -
                         210 * a * b2 * e6 * f * g2 * n2 - 420 * a3 * b2 * e6 * f * g2 * n2 -
                         42 * b2 * c * e5 * g3 * n2 - 252 * a2 * b2 * c * e5 * g3 * n2 -
                         210 * a * b2 * d * e5 * g3 * n2 - 420 * a3 * b2 * d * e5 * g3 * n2 -
                         98 * b3 * e5 * f * g3 * n2 - 588 * a2 * b3 * e5 * f * g3 * n2 -
                         210 * a * b3 * c * e4 * g4 * n2 - 70 * b3 * d * e4 * g4 * n2 -
                         420 * a2 * b3 * d * e4 * g4 * n2 - 420 * a * b4 * e4 * f * g4 * n2 -
                         70 * b4 * c * e3 * g5 * n2 - 210 * a * b4 * d * e3 * g5 * n2 -
                         126 * b5 * e3 * f * g5 * n2 - 42 * b5 * d * e2 * g6 * n2) +
                 (5 * (560 * a4 * c * e9 * f + 1008 * a5 * d * e9 * f + 560 * a4 * b * e9 * f2 +
                          315 * a3 * c2 * e8 * g + 1260 * a4 * c * d * e8 * g +
                          1134 * a5 * d2 * e8 * g + 3150 * a3 * b * c * e8 * f * g +
                          6300 * a4 * b * d * e8 * f * g + 2835 * a3 * b2 * e8 * f2 * g +
                          1080 * a2 * b * c2 * e7 * g2 + 5040 * a3 * b * c * d * e7 * g2 +
                          5040 * a4 * b * d2 * e7 * g2 + 6480 * a2 * b2 * c * e7 * f * g2 +
                          15120 * a3 * b2 * d * e7 * f * g2 + 5400 * a2 * b3 * e7 * f2 * g2 +
                          1260 * a * b2 * c2 * e6 * g3 + 7560 * a2 * b2 * c * d * e6 * g3 +
                          8820 * a3 * b2 * d2 * e6 * g3 + 5880 * a * b3 * c * e6 * f * g3 +
                          17640 * a2 * b3 * d * e6 * f * g3 + 4620 * a * b4 * e6 * f2 * g3 +
                          504 * b3 * c2 * e5 * g4 + 5040 * a * b3 * c * d * e5 * g4 +
                          7560 * a2 * b3 * d2 * e5 * g4 + 2016 * b4 * c * e5 * f * g4 +
                          10080 * a * b4 * d * e5 * f * g4 + 1512 * b5 * e5 * f2 * g4 +
                          1260 * b4 * c * d * e4 * g5 + 3150 * a * b4 * d2 * e4 * g5 +
                          2268 * b5 * d * e4 * f * g5 + 504 * b5 * d2 * e3 * g6 -
                          105 * a3 * e8 * g * n2 - 63 * a5 * e8 * g * n2 -
                          360 * a2 * b * e7 * g2 * n2 - 360 * a4 * b * e7 * g2 * n2 -
                          420 * a * b2 * e6 * g3 * n2 - 840 * a3 * b2 * e6 * g3 * n2 -
                          168 * b3 * e5 * g4 * n2 - 1008 * a2 * b3 * e5 * g4 * n2 -
                          630 * a * b4 * e4 * g5 * n2 - 168 * b5 * e3 * g6 * n2)) /
                     13 +
                 ((252 * a5 * e10 * f + 350 * a4 * c * e9 * g + 630 * a5 * d * e9 * g +
                     1750 * a4 * b * e9 * f * g + 1575 * a3 * b * c * e8 * g2 +
                     3150 * a4 * b * d * e8 * g2 + 4725 * a3 * b2 * e8 * f * g2 +
                     2700 * a2 * b2 * c * e7 * g3 + 6300 * a3 * b2 * d * e7 * g3 +
                     6300 * a2 * b3 * e7 * f * g3 + 2100 * a * b3 * c * e6 * g4 +
                     6300 * a2 * b3 * d * e6 * g4 + 4200 * a * b4 * e6 * f * g4 +
                     630 * b4 * c * e5 * g5 + 3150 * a * b4 * d * e5 * g5 +
                     1134 * b5 * e5 * f * g5 + 630 * b5 * d * e4 * g6)) /
                     7 +
                 (e5 * g *
                     (126 * a5 * e5 + 700 * a4 * b * e4 * g + 1575 * a3 * b2 * e3 * g2 +
                         1800 * a2 * b3 * e2 * g3 + 1050 * a * b4 * e * g4 + 252 * b5 * g5)) /
                     15)) /
                (6300 * n15) +
            ((f2 - n2) *
                ((252 * c5 * d5 * f2 + 1050 * a * c4 * d6 * f2 + 1800 * a2 * c3 * d7 * f2 +
                     1575 * a3 * c2 * d8 * f2 + 700 * a4 * c * d9 * f2 + 126 * a5 * d10 * f2 +
                     1260 * b * c4 * d5 * f3 + 4200 * a * b * c3 * d6 * f3 +
                     5400 * a2 * b * c2 * d7 * f3 + 3150 * a3 * b * c * d8 * f3 +
                     700 * a4 * b * d9 * f3 + 2520 * b2 * c3 * d5 * f4 +
                     6300 * a * b2 * c2 * d6 * f4 + 5400 * a2 * b2 * c * d7 * f4 +
                     1575 * a3 * b2 * d8 * f4 + 2520 * b3 * c2 * d5 * f5 +
                     4200 * a * b3 * c * d6 * f5 + 1800 * a2 * b3 * d7 * f5 +
                     1260 * b4 * c * d5 * f6 + 1050 * a * b4 * d6 * f6 + 252 * b5 * d5 * f7 -
                     252 * c5 * d5 * n2 - 1050 * a * c4 * d6 * n2 - 1800 * a2 * c3 * d7 * n2 -
                     1575 * a3 * c2 * d8 * n2 - 700 * a4 * c * d9 * n2 - 126 * a5 * d10 * n2 -
                     1260 * b * c4 * d5 * f * n2 - 4200 * a * b * c3 * d6 * f * n2 -
                     5400 * a2 * b * c2 * d7 * f * n2 - 3150 * a3 * b * c * d8 * f * n2 -
                     700 * a4 * b * d9 * f * n2 - 840 * c5 * d3 * f2 * n2 -
                     3150 * a * c4 * d4 * f2 * n2 - 840 * (1 + 6 * a2) * c3 * d5 * f2 * n2 -
                     2520 * b2 * c3 * d5 * f2 * n2 - 2100 * a * c2 * d6 * f2 * n2 -
                     4200 * a3 * c2 * d6 * f2 * n2 - 6300 * a * b2 * c2 * d6 * f2 * n2 -
                     1800 * a2 * c * d7 * f2 * n2 - 1800 * a4 * c * d7 * f2 * n2 -
                     5400 * a2 * b2 * c * d7 * f2 * n2 - 525 * a3 * d8 * f2 * n2 -
                     315 * a5 * d8 * f2 * n2 - 1575 * a3 * b2 * d8 * f2 * n2 -
                     4200 * b * c4 * d3 * f3 * n2 - 12600 * a * b * c3 * d4 * f3 * n2 -
                     2520 * b * c2 * d5 * f3 * n2 - 15120 * a2 * b * c2 * d5 * f3 * n2 -
                     2520 * b3 * c2 * d5 * f3 * n2 - 4200 * a * b * c * d6 * f3 * n2 -
                     8400 * a3 * b * c * d6 * f3 * n2 - 4200 * a * b3 * c * d6 * f3 * n2 -
                     1800 * a2 * b * d7 * f3 * n2 - 1800 * a4 * b * d7 * f3 * n2 -
                     1800 * a2 * b3 * d7 * f3 * n2 - 8400 * b2 * c3 * d3 * f4 * n2 -
                     18900 * a * b2 * c2 * d4 * f4 * n2 - 2520 * b2 * c * d5 * f4 * n2 -
                     15120 * a2 * b2 * c * d5 * f4 * n2 - 1260 * b4 * c * d5 * f4 * n2 -
                     2100 * a * b2 * d6 * f4 * n2 - 4200 * a3 * b2 * d6 * f4 * n2 -
                     1050 * a * b4 * d6 * f4 * n2 - 8400 * b3 * c2 * d3 * f5 * n2 -
                     12600 * a * b3 * c * d4 * f5 * n2 - 840 * b3 * d5 * f5 * n2 -
                     5040 * a2 * b3 * d5 * f5 * n2 - 252 * b5 * d5 * f5 * n2 -
                     4200 * b4 * c * d3 * f6 * n2 - 3150 * a * b4 * d4 * f6 * n2 -
                     840 * b5 * d3 * f7 * n2 + 840 * c5 * d3 * n4 + 3150 * a * c4 * d4 * n4 +
                     840 * (1 + 6 * a2) * c3 * d5 * n4 + 2100 * a * c2 * d6 * n4 +
                     4200 * a3 * c2 * d6 * n4 + 1800 * a2 * c * d7 * n4 + 1800 * a4 * c * d7 * n4 +
                     525 * a3 * d8 * n4 + 315 * a5 * d8 * n4 + 4200 * b * c4 * d3 * f * n4 +
                     12600 * a * b * c3 * d4 * f * n4 + 2520 * b * c2 * d5 * f * n4 +
                     15120 * a2 * b * c2 * d5 * f * n4 + 4200 * a * b * c * d6 * f * n4 +
                     8400 * a3 * b * c * d6 * f * n4 + 1800 * a2 * b * d7 * f * n4 +
                     1800 * a4 * b * d7 * f * n4 + 1260 * c5 * d * f2 * n4 +
                     3150 * a * c4 * d2 * f2 * n4 + 1400 * (2 + 3 * a2) * c3 * d3 * f2 * n4 +
                     8400 * b2 * c3 * d3 * f2 * n4 + 6300 * a * c2 * d4 * f2 * n4 +
                     3150 * a3 * c2 * d4 * f2 * n4 + 18900 * a * b2 * c2 * d4 * f2 * n4 +
                     1260 * (1 + 4 * a2 + a4) * c * d5 * f2 * n4 + 2520 * b2 * c * d5 * f2 * n4 +
                     15120 * a2 * b2 * c * d5 * f2 * n4 + 1050 * a * d6 * f2 * n4 +
                     1400 * a3 * d6 * f2 * n4 + 210 * a5 * d6 * f2 * n4 +
                     2100 * a * b2 * d6 * f2 * n4 + 4200 * a3 * b2 * d6 * f2 * n4 +
                     6300 * b * c4 * d * f3 * n4 + 12600 * a * b * c3 * d2 * f3 * n4 +
                     8400 * b * c2 * d3 * f3 * n4 + 12600 * a2 * b * c2 * d3 * f3 * n4 +
                     8400 * b3 * c2 * d3 * f3 * n4 + 6300 * a * (2 + a2) * b * c * d4 * f3 * n4 +
                     12600 * a * b3 * c * d4 * f3 * n4 + 1260 * b * d5 * f3 * n4 +
                     5040 * a2 * b * d5 * f3 * n4 + 1260 * a4 * b * d5 * f3 * n4 +
                     840 * b3 * d5 * f3 * n4 + 5040 * a2 * b3 * d5 * f3 * n4 +
                     12600 * b2 * c3 * d * f4 * n4 + 18900 * a * b2 * c2 * d2 * f4 * n4 +
                     4200 * (2 + 3 * a2) * b2 * c * d3 * f4 * n4 + 4200 * b4 * c * d3 * f4 * n4 +
                     6300 * a * b2 * d4 * f4 * n4 + 3150 * a3 * b2 * d4 * f4 * n4 +
                     3150 * a * b4 * d4 * f4 * n4 + 12600 * b3 * c2 * d * f5 * n4 +
                     12600 * a * b3 * c * d2 * f5 * n4 + 2800 * b3 * d3 * f5 * n4 +
                     4200 * a2 * b3 * d3 * f5 * n4 + 840 * b5 * d3 * f5 * n4 +
                     6300 * b4 * c * d * f6 * n4 + 3150 * a * b4 * d2 * f6 * n4 +
                     1260 * b5 * d * f7 * n4 - 1260 * c5 * d * n6 - 3150 * a * c4 * d2 * n6 -
                     1400 * (2 + 3 * a2) * c3 * d3 * n6 - 6300 * a * c2 * d4 * n6 -
                     3150 * a3 * c2 * d4 * n6 - 1260 * (1 + 4 * a2 + a4) * c * d5 * n6 -
                     1050 * a * d6 * n6 - 1400 * a3 * d6 * n6 - 210 * a5 * d6 * n6 -
                     6300 * b * c4 * d * f * n6 - 12600 * a * b * c3 * d2 * f * n6 -
                     8400 * b * c2 * d3 * f * n6 - 12600 * a2 * b * c2 * d3 * f * n6 -
                     6300 * a * (2 + a2) * b * c * d4 * f * n6 - 1260 * b * d5 * f * n6 -
                     5040 * a2 * b * d5 * f * n6 - 1260 * a4 * b * d5 * f * n6 -
                     4200 * c3 * d * f2 * n6 - 12600 * b2 * c3 * d * f2 * n6 -
                     6300 * a * c2 * d2 * f2 * n6 - 18900 * a * b2 * c2 * d2 * f2 * n6 -
                     4200 * (1 + a2) * c * d3 * f2 * n6 -
                     4200 * (2 + 3 * a2) * b2 * c * d3 * f2 * n6 - 3150 * a * d4 * f2 * n6 -
                     1050 * a3 * d4 * f2 * n6 - 6300 * a * b2 * d4 * f2 * n6 -
                     3150 * a3 * b2 * d4 * f2 * n6 - 12600 * b * c2 * d * f3 * n6 -
                     12600 * b3 * c2 * d * f3 * n6 - 12600 * a * b * c * d2 * f3 * n6 -
                     12600 * a * b3 * c * d2 * f3 * n6 - 4200 * b * d3 * f3 * n6 -
                     4200 * a2 * b * d3 * f3 * n6 - 2800 * b3 * d3 * f3 * n6 -
                     4200 * a2 * b3 * d3 * f3 * n6 - 12600 * b2 * c * d * f4 * n6 -
                     6300 * b4 * c * d * f4 * n6 - 6300 * a * b2 * d2 * f4 * n6 -
                     3150 * a * b4 * d2 * f4 * n6 - 4200 * b3 * d * f5 * n6 -
                     1260 * b5 * d * f5 * n6 + 4200 * c3 * d * n8 + 6300 * a * c2 * d2 * n8 +
                     4200 * (1 + a2) * c * d3 * n8 + 3150 * a * d4 * n8 + 1050 * a3 * d4 * n8 +
                     12600 * b * c2 * d * f * n8 + 12600 * a * b * c * d2 * f * n8 +
                     4200 * b * d3 * f * n8 + 4200 * a2 * b * d3 * f * n8 + 6300 * c * d * f2 * n8 +
                     12600 * b2 * c * d * f2 * n8 + 3150 * a * d2 * f2 * n8 +
                     6300 * a * b2 * d2 * f2 * n8 + 6300 * b * d * f3 * n8 +
                     4200 * b3 * d * f3 * n8 - 6300 * c * d * n10 - 3150 * a * d2 * n10 -
                     6300 * b * d * f * n10) +
                    (630 * c5 * d4 * e * f2 + 3150 * a * c4 * d5 * e * f2 +
                        6300 * a2 * c3 * d6 * e * f2 + 6300 * a3 * c2 * d7 * e * f2 +
                        3150 * a4 * c * d8 * e * f2 + 630 * a5 * d9 * e * f2 +
                        3150 * b * c4 * d4 * e * f3 + 12600 * a * b * c3 * d5 * e * f3 +
                        18900 * a2 * b * c2 * d6 * e * f3 + 12600 * a3 * b * c * d7 * e * f3 +
                        3150 * a4 * b * d8 * e * f3 + 6300 * b2 * c3 * d4 * e * f4 +
                        18900 * a * b2 * c2 * d5 * e * f4 + 18900 * a2 * b2 * c * d6 * e * f4 +
                        6300 * a3 * b2 * d7 * e * f4 + 6300 * b3 * c2 * d4 * e * f5 +
                        12600 * a * b3 * c * d5 * e * f5 + 6300 * a2 * b3 * d6 * e * f5 +
                        3150 * b4 * c * d4 * e * f6 + 3150 * a * b4 * d5 * e * f6 +
                        630 * b5 * d4 * e * f7 + 504 * c5 * d5 * f * g +
                        2100 * a * c4 * d6 * f * g + 3600 * a2 * c3 * d7 * f * g +
                        3150 * a3 * c2 * d8 * f * g + 1400 * a4 * c * d9 * f * g +
                        252 * a5 * d10 * f * g + 3150 * b * c4 * d5 * f2 * g +
                        10500 * a * b * c3 * d6 * f2 * g + 13500 * a2 * b * c2 * d7 * f2 * g +
                        7875 * a3 * b * c * d8 * f2 * g + 1750 * a4 * b * d9 * f2 * g +
                        7560 * b2 * c3 * d5 * f3 * g + 18900 * a * b2 * c2 * d6 * f3 * g +
                        16200 * a2 * b2 * c * d7 * f3 * g + 4725 * a3 * b2 * d8 * f3 * g +
                        8820 * b3 * c2 * d5 * f4 * g + 14700 * a * b3 * c * d6 * f4 * g +
                        6300 * a2 * b3 * d7 * f4 * g + 5040 * b4 * c * d5 * f5 * g +
                        4200 * a * b4 * d6 * f5 * g + 1134 * b5 * d5 * f6 * g -
                        630 * c5 * d4 * e * n2 - 3150 * a * c4 * d5 * e * n2 -
                        6300 * a2 * c3 * d6 * e * n2 - 6300 * a3 * c2 * d7 * e * n2 -
                        3150 * a4 * c * d8 * e * n2 - 630 * a5 * d9 * e * n2 -
                        3150 * b * c4 * d4 * e * f * n2 - 12600 * a * b * c3 * d5 * e * f * n2 -
                        18900 * a2 * b * c2 * d6 * e * f * n2 -
                        12600 * a3 * b * c * d7 * e * f * n2 - 3150 * a4 * b * d8 * e * f * n2 -
                        1260 * c5 * d2 * e * f2 * n2 - 6300 * a * c4 * d3 * e * f2 * n2 -
                        2100 * c3 * d4 * e * f2 * n2 - 12600 * a2 * c3 * d4 * e * f2 * n2 -
                        6300 * b2 * c3 * d4 * e * f2 * n2 - 6300 * a * c2 * d5 * e * f2 * n2 -
                        12600 * a3 * c2 * d5 * e * f2 * n2 -
                        18900 * a * b2 * c2 * d5 * e * f2 * n2 - 6300 * a2 * c * d6 * e * f2 * n2 -
                        6300 * a4 * c * d6 * e * f2 * n2 - 18900 * a2 * b2 * c * d6 * e * f2 * n2 -
                        2100 * a3 * d7 * e * f2 * n2 - 1260 * a5 * d7 * e * f2 * n2 -
                        6300 * a3 * b2 * d7 * e * f2 * n2 - 6300 * b * c4 * d2 * e * f3 * n2 -
                        25200 * a * b * c3 * d3 * e * f3 * n2 - 6300 * b * c2 * d4 * e * f3 * n2 -
                        37800 * a2 * b * c2 * d4 * e * f3 * n2 - 6300 * b3 * c2 * d4 * e * f3 * n2 -
                        12600 * a * b * c * d5 * e * f3 * n2 -
                        25200 * a3 * b * c * d5 * e * f3 * n2 -
                        12600 * a * b3 * c * d5 * e * f3 * n2 - 6300 * a2 * b * d6 * e * f3 * n2 -
                        6300 * a4 * b * d6 * e * f3 * n2 - 6300 * a2 * b3 * d6 * e * f3 * n2 -
                        12600 * b2 * c3 * d2 * e * f4 * n2 -
                        37800 * a * b2 * c2 * d3 * e * f4 * n2 - 6300 * b2 * c * d4 * e * f4 * n2 -
                        37800 * a2 * b2 * c * d4 * e * f4 * n2 - 3150 * b4 * c * d4 * e * f4 * n2 -
                        6300 * a * b2 * d5 * e * f4 * n2 - 12600 * a3 * b2 * d5 * e * f4 * n2 -
                        3150 * a * b4 * d5 * e * f4 * n2 - 12600 * b3 * c2 * d2 * e * f5 * n2 -
                        25200 * a * b3 * c * d3 * e * f5 * n2 - 2100 * b3 * d4 * e * f5 * n2 -
                        12600 * a2 * b3 * d4 * e * f5 * n2 - 630 * b5 * d4 * e * f5 * n2 -
                        6300 * b4 * c * d2 * e * f6 * n2 - 6300 * a * b4 * d3 * e * f6 * n2 -
                        1260 * b5 * d2 * e * f7 * n2 - 630 * b * c4 * d5 * g * n2 -
                        2100 * a * b * c3 * d6 * g * n2 - 2700 * a2 * b * c2 * d7 * g * n2 -
                        1575 * a3 * b * c * d8 * g * n2 - 350 * a4 * b * d9 * g * n2 -
                        1680 * c5 * d3 * f * g * n2 - 6300 * a * c4 * d4 * f * g * n2 -
                        1680 * c3 * d5 * f * g * n2 - 10080 * a2 * c3 * d5 * f * g * n2 -
                        2520 * b2 * c3 * d5 * f * g * n2 - 4200 * a * c2 * d6 * f * g * n2 -
                        8400 * a3 * c2 * d6 * f * g * n2 - 6300 * a * b2 * c2 * d6 * f * g * n2 -
                        3600 * a2 * c * d7 * f * g * n2 - 3600 * a4 * c * d7 * f * g * n2 -
                        5400 * a2 * b2 * c * d7 * f * g * n2 - 1050 * a3 * d8 * f * g * n2 -
                        630 * a5 * d8 * f * g * n2 - 1575 * a3 * b2 * d8 * f * g * n2 -
                        10500 * b * c4 * d3 * f2 * g * n2 - 31500 * a * b * c3 * d4 * f2 * g * n2 -
                        6300 * b * c2 * d5 * f2 * g * n2 - 37800 * a2 * b * c2 * d5 * f2 * g * n2 -
                        3780 * b3 * c2 * d5 * f2 * g * n2 - 10500 * a * b * c * d6 * f2 * g * n2 -
                        21000 * a3 * b * c * d6 * f2 * g * n2 -
                        6300 * a * b3 * c * d6 * f2 * g * n2 - 4500 * a2 * b * d7 * f2 * g * n2 -
                        4500 * a4 * b * d7 * f2 * g * n2 - 2700 * a2 * b3 * d7 * f2 * g * n2 -
                        25200 * b2 * c3 * d3 * f3 * g * n2 -
                        56700 * a * b2 * c2 * d4 * f3 * g * n2 - 7560 * b2 * c * d5 * f3 * g * n2 -
                        45360 * a2 * b2 * c * d5 * f3 * g * n2 - 2520 * b4 * c * d5 * f3 * g * n2 -
                        6300 * a * b2 * d6 * f3 * g * n2 - 12600 * a3 * b2 * d6 * f3 * g * n2 -
                        2100 * a * b4 * d6 * f3 * g * n2 - 29400 * b3 * c2 * d3 * f4 * g * n2 -
                        44100 * a * b3 * c * d4 * f4 * g * n2 - 2940 * b3 * d5 * f4 * g * n2 -
                        17640 * a2 * b3 * d5 * f4 * g * n2 - 630 * b5 * d5 * f4 * g * n2 -
                        16800 * b4 * c * d3 * f5 * g * n2 - 12600 * a * b4 * d4 * f5 * g * n2 -
                        3780 * b5 * d3 * f6 * g * n2 + 1260 * c5 * d2 * e * n4 +
                        6300 * a * c4 * d3 * e * n4 + 2100 * c3 * d4 * e * n4 +
                        12600 * a2 * c3 * d4 * e * n4 + 6300 * a * c2 * d5 * e * n4 +
                        12600 * a3 * c2 * d5 * e * n4 + 6300 * a2 * c * d6 * e * n4 +
                        6300 * a4 * c * d6 * e * n4 + 2100 * a3 * d7 * e * n4 +
                        1260 * a5 * d7 * e * n4 + 6300 * b * c4 * d2 * e * f * n4 +
                        25200 * a * b * c3 * d3 * e * f * n4 + 6300 * b * c2 * d4 * e * f * n4 +
                        37800 * a2 * b * c2 * d4 * e * f * n4 +
                        12600 * a * b * c * d5 * e * f * n4 + 25200 * a3 * b * c * d5 * e * f * n4 +
                        6300 * a2 * b * d6 * e * f * n4 + 6300 * a4 * b * d6 * e * f * n4 +
                        630 * c5 * e * f2 * n4 + 3150 * a * c4 * d * e * f2 * n4 +
                        4200 * c3 * d2 * e * f2 * n4 + 6300 * a2 * c3 * d2 * e * f2 * n4 +
                        12600 * b2 * c3 * d2 * e * f2 * n4 + 12600 * a * c2 * d3 * e * f2 * n4 +
                        6300 * a3 * c2 * d3 * e * f2 * n4 + 37800 * a * b2 * c2 * d3 * e * f2 * n4 +
                        3150 * c * d4 * e * f2 * n4 + 12600 * a2 * c * d4 * e * f2 * n4 +
                        3150 * a4 * c * d4 * e * f2 * n4 + 6300 * b2 * c * d4 * e * f2 * n4 +
                        37800 * a2 * b2 * c * d4 * e * f2 * n4 + 3150 * a * d5 * e * f2 * n4 +
                        4200 * a3 * d5 * e * f2 * n4 + 630 * a5 * d5 * e * f2 * n4 +
                        6300 * a * b2 * d5 * e * f2 * n4 + 12600 * a3 * b2 * d5 * e * f2 * n4 +
                        3150 * b * c4 * e * f3 * n4 + 12600 * a * b * c3 * d * e * f3 * n4 +
                        12600 * b * c2 * d2 * e * f3 * n4 + 18900 * a2 * b * c2 * d2 * e * f3 * n4 +
                        12600 * b3 * c2 * d2 * e * f3 * n4 + 25200 * a * b * c * d3 * e * f3 * n4 +
                        12600 * a3 * b * c * d3 * e * f3 * n4 +
                        25200 * a * b3 * c * d3 * e * f3 * n4 + 3150 * b * d4 * e * f3 * n4 +
                        12600 * a2 * b * d4 * e * f3 * n4 + 3150 * a4 * b * d4 * e * f3 * n4 +
                        2100 * b3 * d4 * e * f3 * n4 + 12600 * a2 * b3 * d4 * e * f3 * n4 +
                        6300 * b2 * c3 * e * f4 * n4 + 18900 * a * b2 * c2 * d * e * f4 * n4 +
                        12600 * b2 * c * d2 * e * f4 * n4 + 18900 * a2 * b2 * c * d2 * e * f4 * n4 +
                        6300 * b4 * c * d2 * e * f4 * n4 + 12600 * a * b2 * d3 * e * f4 * n4 +
                        6300 * a3 * b2 * d3 * e * f4 * n4 + 6300 * a * b4 * d3 * e * f4 * n4 +
                        6300 * b3 * c2 * e * f5 * n4 + 12600 * a * b3 * c * d * e * f5 * n4 +
                        4200 * b3 * d2 * e * f5 * n4 + 6300 * a2 * b3 * d2 * e * f5 * n4 +
                        1260 * b5 * d2 * e * f5 * n4 + 3150 * b4 * c * e * f6 * n4 +
                        3150 * a * b4 * d * e * f6 * n4 + 630 * b5 * e * f7 * n4 +
                        2100 * b * c4 * d3 * g * n4 + 6300 * a * b * c3 * d4 * g * n4 +
                        1260 * b * c2 * d5 * g * n4 + 7560 * a2 * b * c2 * d5 * g * n4 +
                        2100 * a * b * c * d6 * g * n4 + 4200 * a3 * b * c * d6 * g * n4 +
                        900 * a2 * b * d7 * g * n4 + 900 * a4 * b * d7 * g * n4 +
                        2520 * c5 * d * f * g * n4 + 6300 * a * c4 * d2 * f * g * n4 +
                        5600 * c3 * d3 * f * g * n4 + 8400 * a2 * c3 * d3 * f * g * n4 +
                        8400 * b2 * c3 * d3 * f * g * n4 + 12600 * a * c2 * d4 * f * g * n4 +
                        6300 * a3 * c2 * d4 * f * g * n4 + 18900 * a * b2 * c2 * d4 * f * g * n4 +
                        2520 * c * d5 * f * g * n4 + 10080 * a2 * c * d5 * f * g * n4 +
                        2520 * a4 * c * d5 * f * g * n4 + 2520 * b2 * c * d5 * f * g * n4 +
                        15120 * a2 * b2 * c * d5 * f * g * n4 + 2100 * a * d6 * f * g * n4 +
                        2800 * a3 * d6 * f * g * n4 + 420 * a5 * d6 * f * g * n4 +
                        2100 * a * b2 * d6 * f * g * n4 + 4200 * a3 * b2 * d6 * f * g * n4 +
                        15750 * b * c4 * d * f2 * g * n4 + 31500 * a * b * c3 * d2 * f2 * g * n4 +
                        21000 * b * c2 * d3 * f2 * g * n4 + 31500 * a2 * b * c2 * d3 * f2 * g * n4 +
                        12600 * b3 * c2 * d3 * f2 * g * n4 + 31500 * a * b * c * d4 * f2 * g * n4 +
                        15750 * a3 * b * c * d4 * f2 * g * n4 +
                        18900 * a * b3 * c * d4 * f2 * g * n4 + 3150 * b * d5 * f2 * g * n4 +
                        12600 * a2 * b * d5 * f2 * g * n4 + 3150 * a4 * b * d5 * f2 * g * n4 +
                        1260 * b3 * d5 * f2 * g * n4 + 7560 * a2 * b3 * d5 * f2 * g * n4 +
                        37800 * b2 * c3 * d * f3 * g * n4 + 56700 * a * b2 * c2 * d2 * f3 * g * n4 +
                        25200 * b2 * c * d3 * f3 * g * n4 + 37800 * a2 * b2 * c * d3 * f3 * g * n4 +
                        8400 * b4 * c * d3 * f3 * g * n4 + 18900 * a * b2 * d4 * f3 * g * n4 +
                        9450 * a3 * b2 * d4 * f3 * g * n4 + 6300 * a * b4 * d4 * f3 * g * n4 +
                        44100 * b3 * c2 * d * f4 * g * n4 + 44100 * a * b3 * c * d2 * f4 * g * n4 +
                        9800 * b3 * d3 * f4 * g * n4 + 14700 * a2 * b3 * d3 * f4 * g * n4 +
                        2100 * b5 * d3 * f4 * g * n4 + 25200 * b4 * c * d * f5 * g * n4 +
                        12600 * a * b4 * d2 * f5 * g * n4 + 5670 * b5 * d * f6 * g * n4 -
                        630 * c5 * e * n6 - 3150 * a * c4 * d * e * n6 - 4200 * c3 * d2 * e * n6 -
                        6300 * a2 * c3 * d2 * e * n6 - 12600 * a * c2 * d3 * e * n6 -
                        6300 * a3 * c2 * d3 * e * n6 - 3150 * c * d4 * e * n6 -
                        12600 * a2 * c * d4 * e * n6 - 3150 * a4 * c * d4 * e * n6 -
                        3150 * a * d5 * e * n6 - 4200 * a3 * d5 * e * n6 - 630 * a5 * d5 * e * n6 -
                        3150 * b * c4 * e * f * n6 - 12600 * a * b * c3 * d * e * f * n6 -
                        12600 * b * c2 * d2 * e * f * n6 - 18900 * a2 * b * c2 * d2 * e * f * n6 -
                        25200 * a * b * c * d3 * e * f * n6 - 12600 * a3 * b * c * d3 * e * f * n6 -
                        3150 * b * d4 * e * f * n6 - 12600 * a2 * b * d4 * e * f * n6 -
                        3150 * a4 * b * d4 * e * f * n6 - 2100 * c3 * e * f2 * n6 -
                        6300 * b2 * c3 * e * f2 * n6 - 6300 * a * c2 * d * e * f2 * n6 -
                        18900 * a * b2 * c2 * d * e * f2 * n6 - 6300 * c * d2 * e * f2 * n6 -
                        6300 * a2 * c * d2 * e * f2 * n6 - 12600 * b2 * c * d2 * e * f2 * n6 -
                        18900 * a2 * b2 * c * d2 * e * f2 * n6 - 6300 * a * d3 * e * f2 * n6 -
                        2100 * a3 * d3 * e * f2 * n6 - 12600 * a * b2 * d3 * e * f2 * n6 -
                        6300 * a3 * b2 * d3 * e * f2 * n6 - 6300 * b * c2 * e * f3 * n6 -
                        6300 * b3 * c2 * e * f3 * n6 - 12600 * a * b * c * d * e * f3 * n6 -
                        12600 * a * b3 * c * d * e * f3 * n6 - 6300 * b * d2 * e * f3 * n6 -
                        6300 * a2 * b * d2 * e * f3 * n6 - 4200 * b3 * d2 * e * f3 * n6 -
                        6300 * a2 * b3 * d2 * e * f3 * n6 - 6300 * b2 * c * e * f4 * n6 -
                        3150 * b4 * c * e * f4 * n6 - 6300 * a * b2 * d * e * f4 * n6 -
                        3150 * a * b4 * d * e * f4 * n6 - 2100 * b3 * e * f5 * n6 -
                        630 * b5 * e * f5 * n6 - 3150 * b * c4 * d * g * n6 -
                        6300 * a * b * c3 * d2 * g * n6 - 4200 * b * c2 * d3 * g * n6 -
                        6300 * a2 * b * c2 * d3 * g * n6 - 6300 * a * b * c * d4 * g * n6 -
                        3150 * a3 * b * c * d4 * g * n6 - 630 * b * d5 * g * n6 -
                        2520 * a2 * b * d5 * g * n6 - 630 * a4 * b * d5 * g * n6 -
                        8400 * c3 * d * f * g * n6 - 12600 * b2 * c3 * d * f * g * n6 -
                        12600 * a * c2 * d2 * f * g * n6 - 18900 * a * b2 * c2 * d2 * f * g * n6 -
                        8400 * c * d3 * f * g * n6 - 8400 * a2 * c * d3 * f * g * n6 -
                        8400 * b2 * c * d3 * f * g * n6 - 12600 * a2 * b2 * c * d3 * f * g * n6 -
                        6300 * a * d4 * f * g * n6 - 2100 * a3 * d4 * f * g * n6 -
                        6300 * a * b2 * d4 * f * g * n6 - 3150 * a3 * b2 * d4 * f * g * n6 -
                        31500 * b * c2 * d * f2 * g * n6 - 18900 * b3 * c2 * d * f2 * g * n6 -
                        31500 * a * b * c * d2 * f2 * g * n6 -
                        18900 * a * b3 * c * d2 * f2 * g * n6 - 10500 * b * d3 * f2 * g * n6 -
                        10500 * a2 * b * d3 * f2 * g * n6 - 4200 * b3 * d3 * f2 * g * n6 -
                        6300 * a2 * b3 * d3 * f2 * g * n6 - 37800 * b2 * c * d * f3 * g * n6 -
                        12600 * b4 * c * d * f3 * g * n6 - 18900 * a * b2 * d2 * f3 * g * n6 -
                        6300 * a * b4 * d2 * f3 * g * n6 - 14700 * b3 * d * f4 * g * n6 -
                        3150 * b5 * d * f4 * g * n6 + 2100 * c3 * e * n8 +
                        6300 * a * c2 * d * e * n8 + 6300 * c * d2 * e * n8 +
                        6300 * a2 * c * d2 * e * n8 + 6300 * a * d3 * e * n8 +
                        2100 * a3 * d3 * e * n8 + 6300 * b * c2 * e * f * n8 +
                        12600 * a * b * c * d * e * f * n8 + 6300 * b * d2 * e * f * n8 +
                        6300 * a2 * b * d2 * e * f * n8 + 3150 * c * e * f2 * n8 +
                        6300 * b2 * c * e * f2 * n8 + 3150 * a * d * e * f2 * n8 +
                        6300 * a * b2 * d * e * f2 * n8 + 3150 * b * e * f3 * n8 +
                        2100 * b3 * e * f3 * n8 + 6300 * b * c2 * d * g * n8 +
                        6300 * a * b * c * d2 * g * n8 + 2100 * b * d3 * g * n8 +
                        2100 * a2 * b * d3 * g * n8 + 12600 * c * d * f * g * n8 +
                        12600 * b2 * c * d * f * g * n8 + 6300 * a * d2 * f * g * n8 +
                        6300 * a * b2 * d2 * f * g * n8 + 15750 * b * d * f2 * g * n8 +
                        6300 * b3 * d * f2 * g * n8 - 3150 * c * e * n10 - 3150 * a * d * e * n10 -
                        3150 * b * e * f * n10 - 3150 * b * d * g * n10) +
                    (5 *
                        (504 * c5 * d3 * e2 * f2 + 3150 * a * c4 * d4 * e2 * f2 +
                            7560 * a2 * c3 * d5 * e2 * f2 + 8820 * a3 * c2 * d6 * e2 * f2 +
                            5040 * a4 * c * d7 * e2 * f2 + 1134 * a5 * d8 * e2 * f2 +
                            2520 * b * c4 * d3 * e2 * f3 + 12600 * a * b * c3 * d4 * e2 * f3 +
                            22680 * a2 * b * c2 * d5 * e2 * f3 + 17640 * a3 * b * c * d6 * e2 * f3 +
                            5040 * a4 * b * d7 * e2 * f3 + 5040 * b2 * c3 * d3 * e2 * f4 +
                            18900 * a * b2 * c2 * d4 * e2 * f4 +
                            22680 * a2 * b2 * c * d5 * e2 * f4 + 8820 * a3 * b2 * d6 * e2 * f4 +
                            5040 * b3 * c2 * d3 * e2 * f5 + 12600 * a * b3 * c * d4 * e2 * f5 +
                            7560 * a2 * b3 * d5 * e2 * f5 + 2520 * b4 * c * d3 * e2 * f6 +
                            3150 * a * b4 * d4 * e2 * f6 + 504 * b5 * d3 * e2 * f7 +
                            1008 * c5 * d4 * e * f * g + 5040 * a * c4 * d5 * e * f * g +
                            10080 * a2 * c3 * d6 * e * f * g + 10080 * a3 * c2 * d7 * e * f * g +
                            5040 * a4 * c * d8 * e * f * g + 1008 * a5 * d9 * e * f * g +
                            6300 * b * c4 * d4 * e * f2 * g + 25200 * a * b * c3 * d5 * e * f2 * g +
                            37800 * a2 * b * c2 * d6 * e * f2 * g +
                            25200 * a3 * b * c * d7 * e * f2 * g + 6300 * a4 * b * d8 * e * f2 * g +
                            15120 * b2 * c3 * d4 * e * f3 * g +
                            45360 * a * b2 * c2 * d5 * e * f3 * g +
                            45360 * a2 * b2 * c * d6 * e * f3 * g +
                            15120 * a3 * b2 * d7 * e * f3 * g + 17640 * b3 * c2 * d4 * e * f4 * g +
                            35280 * a * b3 * c * d5 * e * f4 * g +
                            17640 * a2 * b3 * d6 * e * f4 * g + 10080 * b4 * c * d4 * e * f5 * g +
                            10080 * a * b4 * d5 * e * f5 * g + 2268 * b5 * d4 * e * f6 * g +
                            1008 * b * c4 * d5 * f * g2 + 3360 * a * b * c3 * d6 * f * g2 +
                            4320 * a2 * b * c2 * d7 * f * g2 + 2520 * a3 * b * c * d8 * f * g2 +
                            560 * a4 * b * d9 * f * g2 + 4536 * b2 * c3 * d5 * f2 * g2 +
                            11340 * a * b2 * c2 * d6 * f2 * g2 + 9720 * a2 * b2 * c * d7 * f2 * g2 +
                            2835 * a3 * b2 * d8 * f2 * g2 + 7560 * b3 * c2 * d5 * f3 * g2 +
                            12600 * a * b3 * c * d6 * f3 * g2 + 5400 * a2 * b3 * d7 * f3 * g2 +
                            5544 * b4 * c * d5 * f4 * g2 + 4620 * a * b4 * d6 * f4 * g2 +
                            1512 * b5 * d5 * f5 * g2 - 504 * c5 * d3 * e2 * n2 -
                            3150 * a * c4 * d4 * e2 * n2 - 7560 * a2 * c3 * d5 * e2 * n2 -
                            8820 * a3 * c2 * d6 * e2 * n2 - 5040 * a4 * c * d7 * e2 * n2 -
                            1134 * a5 * d8 * e2 * n2 - 2520 * b * c4 * d3 * e2 * f * n2 -
                            12600 * a * b * c3 * d4 * e2 * f * n2 -
                            22680 * a2 * b * c2 * d5 * e2 * f * n2 -
                            17640 * a3 * b * c * d6 * e2 * f * n2 -
                            5040 * a4 * b * d7 * e2 * f * n2 - 504 * c5 * d * e2 * f2 * n2 -
                            3780 * a * c4 * d2 * e2 * f2 * n2 - 1680 * c3 * d3 * e2 * f2 * n2 -
                            10080 * a2 * c3 * d3 * e2 * f2 * n2 -
                            5040 * b2 * c3 * d3 * e2 * f2 * n2 - 6300 * a * c2 * d4 * e2 * f2 * n2 -
                            12600 * a3 * c2 * d4 * e2 * f2 * n2 -
                            18900 * a * b2 * c2 * d4 * e2 * f2 * n2 -
                            7560 * a2 * c * d5 * e2 * f2 * n2 - 7560 * a4 * c * d5 * e2 * f2 * n2 -
                            22680 * a2 * b2 * c * d5 * e2 * f2 * n2 -
                            2940 * a3 * d6 * e2 * f2 * n2 - 1764 * a5 * d6 * e2 * f2 * n2 -
                            8820 * a3 * b2 * d6 * e2 * f2 * n2 - 2520 * b * c4 * d * e2 * f3 * n2 -
                            15120 * a * b * c3 * d2 * e2 * f3 * n2 -
                            5040 * b * c2 * d3 * e2 * f3 * n2 -
                            30240 * a2 * b * c2 * d3 * e2 * f3 * n2 -
                            5040 * b3 * c2 * d3 * e2 * f3 * n2 -
                            12600 * a * b * c * d4 * e2 * f3 * n2 -
                            25200 * a3 * b * c * d4 * e2 * f3 * n2 -
                            12600 * a * b3 * c * d4 * e2 * f3 * n2 -
                            7560 * a2 * b * d5 * e2 * f3 * n2 - 7560 * a4 * b * d5 * e2 * f3 * n2 -
                            7560 * a2 * b3 * d5 * e2 * f3 * n2 - 5040 * b2 * c3 * d * e2 * f4 * n2 -
                            22680 * a * b2 * c2 * d2 * e2 * f4 * n2 -
                            5040 * b2 * c * d3 * e2 * f4 * n2 -
                            30240 * a2 * b2 * c * d3 * e2 * f4 * n2 -
                            2520 * b4 * c * d3 * e2 * f4 * n2 - 6300 * a * b2 * d4 * e2 * f4 * n2 -
                            12600 * a3 * b2 * d4 * e2 * f4 * n2 -
                            3150 * a * b4 * d4 * e2 * f4 * n2 - 5040 * b3 * c2 * d * e2 * f5 * n2 -
                            15120 * a * b3 * c * d2 * e2 * f5 * n2 - 1680 * b3 * d3 * e2 * f5 * n2 -
                            10080 * a2 * b3 * d3 * e2 * f5 * n2 - 504 * b5 * d3 * e2 * f5 * n2 -
                            2520 * b4 * c * d * e2 * f6 * n2 - 3780 * a * b4 * d2 * e2 * f6 * n2 -
                            504 * b5 * d * e2 * f7 * n2 - 1260 * b * c4 * d4 * e * g * n2 -
                            5040 * a * b * c3 * d5 * e * g * n2 -
                            7560 * a2 * b * c2 * d6 * e * g * n2 -
                            5040 * a3 * b * c * d7 * e * g * n2 - 1260 * a4 * b * d8 * e * g * n2 -
                            2016 * c5 * d2 * e * f * g * n2 - 10080 * a * c4 * d3 * e * f * g * n2 -
                            3360 * c3 * d4 * e * f * g * n2 -
                            20160 * a2 * c3 * d4 * e * f * g * n2 -
                            5040 * b2 * c3 * d4 * e * f * g * n2 -
                            10080 * a * c2 * d5 * e * f * g * n2 -
                            20160 * a3 * c2 * d5 * e * f * g * n2 -
                            15120 * a * b2 * c2 * d5 * e * f * g * n2 -
                            10080 * a2 * c * d6 * e * f * g * n2 -
                            10080 * a4 * c * d6 * e * f * g * n2 -
                            15120 * a2 * b2 * c * d6 * e * f * g * n2 -
                            3360 * a3 * d7 * e * f * g * n2 - 2016 * a5 * d7 * e * f * g * n2 -
                            5040 * a3 * b2 * d7 * e * f * g * n2 -
                            12600 * b * c4 * d2 * e * f2 * g * n2 -
                            50400 * a * b * c3 * d3 * e * f2 * g * n2 -
                            12600 * b * c2 * d4 * e * f2 * g * n2 -
                            75600 * a2 * b * c2 * d4 * e * f2 * g * n2 -
                            7560 * b3 * c2 * d4 * e * f2 * g * n2 -
                            25200 * a * b * c * d5 * e * f2 * g * n2 -
                            50400 * a3 * b * c * d5 * e * f2 * g * n2 -
                            15120 * a * b3 * c * d5 * e * f2 * g * n2 -
                            12600 * a2 * b * d6 * e * f2 * g * n2 -
                            12600 * a4 * b * d6 * e * f2 * g * n2 -
                            7560 * a2 * b3 * d6 * e * f2 * g * n2 -
                            30240 * b2 * c3 * d2 * e * f3 * g * n2 -
                            90720 * a * b2 * c2 * d3 * e * f3 * g * n2 -
                            15120 * b2 * c * d4 * e * f3 * g * n2 -
                            90720 * a2 * b2 * c * d4 * e * f3 * g * n2 -
                            5040 * b4 * c * d4 * e * f3 * g * n2 -
                            15120 * a * b2 * d5 * e * f3 * g * n2 -
                            30240 * a3 * b2 * d5 * e * f3 * g * n2 -
                            5040 * a * b4 * d5 * e * f3 * g * n2 -
                            35280 * b3 * c2 * d2 * e * f4 * g * n2 -
                            70560 * a * b3 * c * d3 * e * f4 * g * n2 -
                            5880 * b3 * d4 * e * f4 * g * n2 -
                            35280 * a2 * b3 * d4 * e * f4 * g * n2 -
                            1260 * b5 * d4 * e * f4 * g * n2 -
                            20160 * b4 * c * d2 * e * f5 * g * n2 -
                            20160 * a * b4 * d3 * e * f5 * g * n2 -
                            4536 * b5 * d2 * e * f6 * g * n2 - 504 * b2 * c3 * d5 * g2 * n2 -
                            1260 * a * b2 * c2 * d6 * g2 * n2 - 1080 * a2 * b2 * c * d7 * g2 * n2 -
                            315 * a3 * b2 * d8 * g2 * n2 - 3360 * b * c4 * d3 * f * g2 * n2 -
                            10080 * a * b * c3 * d4 * f * g2 * n2 -
                            2016 * b * c2 * d5 * f * g2 * n2 -
                            12096 * a2 * b * c2 * d5 * f * g2 * n2 -
                            1512 * b3 * c2 * d5 * f * g2 * n2 -
                            3360 * a * b * c * d6 * f * g2 * n2 -
                            6720 * a3 * b * c * d6 * f * g2 * n2 -
                            2520 * a * b3 * c * d6 * f * g2 * n2 -
                            1440 * a2 * b * d7 * f * g2 * n2 - 1440 * a4 * b * d7 * f * g2 * n2 -
                            1080 * a2 * b3 * d7 * f * g2 * n2 -
                            15120 * b2 * c3 * d3 * f2 * g2 * n2 -
                            34020 * a * b2 * c2 * d4 * f2 * g2 * n2 -
                            4536 * b2 * c * d5 * f2 * g2 * n2 -
                            27216 * a2 * b2 * c * d5 * f2 * g2 * n2 -
                            1512 * b4 * c * d5 * f2 * g2 * n2 - 3780 * a * b2 * d6 * f2 * g2 * n2 -
                            7560 * a3 * b2 * d6 * f2 * g2 * n2 - 1260 * a * b4 * d6 * f2 * g2 * n2 -
                            25200 * b3 * c2 * d3 * f3 * g2 * n2 -
                            37800 * a * b3 * c * d4 * f3 * g2 * n2 - 2520 * b3 * d5 * f3 * g2 * n2 -
                            15120 * a2 * b3 * d5 * f3 * g2 * n2 - 504 * b5 * d5 * f3 * g2 * n2 -
                            18480 * b4 * c * d3 * f4 * g2 * n2 -
                            13860 * a * b4 * d4 * f4 * g2 * n2 - 5040 * b5 * d3 * f5 * g2 * n2 +
                            504 * c5 * d * e2 * n4 + 3780 * a * c4 * d2 * e2 * n4 +
                            1680 * c3 * d3 * e2 * n4 + 10080 * a2 * c3 * d3 * e2 * n4 +
                            6300 * a * c2 * d4 * e2 * n4 + 12600 * a3 * c2 * d4 * e2 * n4 +
                            7560 * a2 * c * d5 * e2 * n4 + 7560 * a4 * c * d5 * e2 * n4 +
                            2940 * a3 * d6 * e2 * n4 + 1764 * a5 * d6 * e2 * n4 +
                            2520 * b * c4 * d * e2 * f * n4 +
                            15120 * a * b * c3 * d2 * e2 * f * n4 +
                            5040 * b * c2 * d3 * e2 * f * n4 +
                            30240 * a2 * b * c2 * d3 * e2 * f * n4 +
                            12600 * a * b * c * d4 * e2 * f * n4 +
                            25200 * a3 * b * c * d4 * e2 * f * n4 +
                            7560 * a2 * b * d5 * e2 * f * n4 + 7560 * a4 * b * d5 * e2 * f * n4 +
                            630 * a * c4 * e2 * f2 * n4 + 1680 * c3 * d * e2 * f2 * n4 +
                            2520 * a2 * c3 * d * e2 * f2 * n4 + 5040 * b2 * c3 * d * e2 * f2 * n4 +
                            7560 * a * c2 * d2 * e2 * f2 * n4 + 3780 * a3 * c2 * d2 * e2 * f2 * n4 +
                            22680 * a * b2 * c2 * d2 * e2 * f2 * n4 + 2520 * c * d3 * e2 * f2 * n4 +
                            10080 * a2 * c * d3 * e2 * f2 * n4 + 2520 * a4 * c * d3 * e2 * f2 * n4 +
                            5040 * b2 * c * d3 * e2 * f2 * n4 +
                            30240 * a2 * b2 * c * d3 * e2 * f2 * n4 + 3150 * a * d4 * e2 * f2 * n4 +
                            4200 * a3 * d4 * e2 * f2 * n4 + 630 * a5 * d4 * e2 * f2 * n4 +
                            6300 * a * b2 * d4 * e2 * f2 * n4 +
                            12600 * a3 * b2 * d4 * e2 * f2 * n4 + 2520 * a * b * c3 * e2 * f3 * n4 +
                            5040 * b * c2 * d * e2 * f3 * n4 +
                            7560 * a2 * b * c2 * d * e2 * f3 * n4 +
                            5040 * b3 * c2 * d * e2 * f3 * n4 +
                            15120 * a * b * c * d2 * e2 * f3 * n4 +
                            7560 * a3 * b * c * d2 * e2 * f3 * n4 +
                            15120 * a * b3 * c * d2 * e2 * f3 * n4 + 2520 * b * d3 * e2 * f3 * n4 +
                            10080 * a2 * b * d3 * e2 * f3 * n4 + 2520 * a4 * b * d3 * e2 * f3 * n4 +
                            1680 * b3 * d3 * e2 * f3 * n4 + 10080 * a2 * b3 * d3 * e2 * f3 * n4 +
                            3780 * a * b2 * c2 * e2 * f4 * n4 + 5040 * b2 * c * d * e2 * f4 * n4 +
                            7560 * a2 * b2 * c * d * e2 * f4 * n4 +
                            2520 * b4 * c * d * e2 * f4 * n4 + 7560 * a * b2 * d2 * e2 * f4 * n4 +
                            3780 * a3 * b2 * d2 * e2 * f4 * n4 + 3780 * a * b4 * d2 * e2 * f4 * n4 +
                            2520 * a * b3 * c * e2 * f5 * n4 + 1680 * b3 * d * e2 * f5 * n4 +
                            2520 * a2 * b3 * d * e2 * f5 * n4 + 504 * b5 * d * e2 * f5 * n4 +
                            630 * a * b4 * e2 * f6 * n4 + 2520 * b * c4 * d2 * e * g * n4 +
                            10080 * a * b * c3 * d3 * e * g * n4 + 2520 * b * c2 * d4 * e * g * n4 +
                            15120 * a2 * b * c2 * d4 * e * g * n4 +
                            5040 * a * b * c * d5 * e * g * n4 +
                            10080 * a3 * b * c * d5 * e * g * n4 + 2520 * a2 * b * d6 * e * g * n4 +
                            2520 * a4 * b * d6 * e * g * n4 + 1008 * c5 * e * f * g * n4 +
                            5040 * a * c4 * d * e * f * g * n4 + 6720 * c3 * d2 * e * f * g * n4 +
                            10080 * a2 * c3 * d2 * e * f * g * n4 +
                            10080 * b2 * c3 * d2 * e * f * g * n4 +
                            20160 * a * c2 * d3 * e * f * g * n4 +
                            10080 * a3 * c2 * d3 * e * f * g * n4 +
                            30240 * a * b2 * c2 * d3 * e * f * g * n4 +
                            5040 * c * d4 * e * f * g * n4 + 20160 * a2 * c * d4 * e * f * g * n4 +
                            5040 * a4 * c * d4 * e * f * g * n4 +
                            5040 * b2 * c * d4 * e * f * g * n4 +
                            30240 * a2 * b2 * c * d4 * e * f * g * n4 +
                            5040 * a * d5 * e * f * g * n4 + 6720 * a3 * d5 * e * f * g * n4 +
                            1008 * a5 * d5 * e * f * g * n4 + 5040 * a * b2 * d5 * e * f * g * n4 +
                            10080 * a3 * b2 * d5 * e * f * g * n4 +
                            6300 * b * c4 * e * f2 * g * n4 +
                            25200 * a * b * c3 * d * e * f2 * g * n4 +
                            25200 * b * c2 * d2 * e * f2 * g * n4 +
                            37800 * a2 * b * c2 * d2 * e * f2 * g * n4 +
                            15120 * b3 * c2 * d2 * e * f2 * g * n4 +
                            50400 * a * b * c * d3 * e * f2 * g * n4 +
                            25200 * a3 * b * c * d3 * e * f2 * g * n4 +
                            30240 * a * b3 * c * d3 * e * f2 * g * n4 +
                            6300 * b * d4 * e * f2 * g * n4 +
                            25200 * a2 * b * d4 * e * f2 * g * n4 +
                            6300 * a4 * b * d4 * e * f2 * g * n4 +
                            2520 * b3 * d4 * e * f2 * g * n4 +
                            15120 * a2 * b3 * d4 * e * f2 * g * n4 +
                            15120 * b2 * c3 * e * f3 * g * n4 +
                            45360 * a * b2 * c2 * d * e * f3 * g * n4 +
                            30240 * b2 * c * d2 * e * f3 * g * n4 +
                            45360 * a2 * b2 * c * d2 * e * f3 * g * n4 +
                            10080 * b4 * c * d2 * e * f3 * g * n4 +
                            30240 * a * b2 * d3 * e * f3 * g * n4 +
                            15120 * a3 * b2 * d3 * e * f3 * g * n4 +
                            10080 * a * b4 * d3 * e * f3 * g * n4 +
                            17640 * b3 * c2 * e * f4 * g * n4 +
                            35280 * a * b3 * c * d * e * f4 * g * n4 +
                            11760 * b3 * d2 * e * f4 * g * n4 +
                            17640 * a2 * b3 * d2 * e * f4 * g * n4 +
                            2520 * b5 * d2 * e * f4 * g * n4 + 10080 * b4 * c * e * f5 * g * n4 +
                            10080 * a * b4 * d * e * f5 * g * n4 + 2268 * b5 * e * f6 * g * n4 +
                            1680 * b2 * c3 * d3 * g2 * n4 + 3780 * a * b2 * c2 * d4 * g2 * n4 +
                            504 * b2 * c * d5 * g2 * n4 + 3024 * a2 * b2 * c * d5 * g2 * n4 +
                            420 * a * b2 * d6 * g2 * n4 + 840 * a3 * b2 * d6 * g2 * n4 +
                            5040 * b * c4 * d * f * g2 * n4 +
                            10080 * a * b * c3 * d2 * f * g2 * n4 +
                            6720 * b * c2 * d3 * f * g2 * n4 +
                            10080 * a2 * b * c2 * d3 * f * g2 * n4 +
                            5040 * b3 * c2 * d3 * f * g2 * n4 +
                            10080 * a * b * c * d4 * f * g2 * n4 +
                            5040 * a3 * b * c * d4 * f * g2 * n4 +
                            7560 * a * b3 * c * d4 * f * g2 * n4 + 1008 * b * d5 * f * g2 * n4 +
                            4032 * a2 * b * d5 * f * g2 * n4 + 1008 * a4 * b * d5 * f * g2 * n4 +
                            504 * b3 * d5 * f * g2 * n4 + 3024 * a2 * b3 * d5 * f * g2 * n4 +
                            22680 * b2 * c3 * d * f2 * g2 * n4 +
                            34020 * a * b2 * c2 * d2 * f2 * g2 * n4 +
                            15120 * b2 * c * d3 * f2 * g2 * n4 +
                            22680 * a2 * b2 * c * d3 * f2 * g2 * n4 +
                            5040 * b4 * c * d3 * f2 * g2 * n4 + 11340 * a * b2 * d4 * f2 * g2 * n4 +
                            5670 * a3 * b2 * d4 * f2 * g2 * n4 + 3780 * a * b4 * d4 * f2 * g2 * n4 +
                            37800 * b3 * c2 * d * f3 * g2 * n4 +
                            37800 * a * b3 * c * d2 * f3 * g2 * n4 + 8400 * b3 * d3 * f3 * g2 * n4 +
                            12600 * a2 * b3 * d3 * f3 * g2 * n4 + 1680 * b5 * d3 * f3 * g2 * n4 +
                            27720 * b4 * c * d * f4 * g2 * n4 + 13860 * a * b4 * d2 * f4 * g2 * n4 +
                            7560 * b5 * d * f5 * g2 * n4 - 630 * a * c4 * e2 * n6 -
                            1680 * c3 * d * e2 * n6 - 2520 * a2 * c3 * d * e2 * n6 -
                            7560 * a * c2 * d2 * e2 * n6 - 3780 * a3 * c2 * d2 * e2 * n6 -
                            2520 * c * d3 * e2 * n6 - 10080 * a2 * c * d3 * e2 * n6 -
                            2520 * a4 * c * d3 * e2 * n6 - 3150 * a * d4 * e2 * n6 -
                            4200 * a3 * d4 * e2 * n6 - 630 * a5 * d4 * e2 * n6 -
                            2520 * a * b * c3 * e2 * f * n6 - 5040 * b * c2 * d * e2 * f * n6 -
                            7560 * a2 * b * c2 * d * e2 * f * n6 -
                            15120 * a * b * c * d2 * e2 * f * n6 -
                            7560 * a3 * b * c * d2 * e2 * f * n6 - 2520 * b * d3 * e2 * f * n6 -
                            10080 * a2 * b * d3 * e2 * f * n6 - 2520 * a4 * b * d3 * e2 * f * n6 -
                            1260 * a * c2 * e2 * f2 * n6 - 3780 * a * b2 * c2 * e2 * f2 * n6 -
                            2520 * c * d * e2 * f2 * n6 - 2520 * a2 * c * d * e2 * f2 * n6 -
                            5040 * b2 * c * d * e2 * f2 * n6 -
                            7560 * a2 * b2 * c * d * e2 * f2 * n6 - 3780 * a * d2 * e2 * f2 * n6 -
                            1260 * a3 * d2 * e2 * f2 * n6 - 7560 * a * b2 * d2 * e2 * f2 * n6 -
                            3780 * a3 * b2 * d2 * e2 * f2 * n6 - 2520 * a * b * c * e2 * f3 * n6 -
                            2520 * a * b3 * c * e2 * f3 * n6 - 2520 * b * d * e2 * f3 * n6 -
                            2520 * a2 * b * d * e2 * f3 * n6 - 1680 * b3 * d * e2 * f3 * n6 -
                            2520 * a2 * b3 * d * e2 * f3 * n6 - 1260 * a * b2 * e2 * f4 * n6 -
                            630 * a * b4 * e2 * f4 * n6 - 1260 * b * c4 * e * g * n6 -
                            5040 * a * b * c3 * d * e * g * n6 - 5040 * b * c2 * d2 * e * g * n6 -
                            7560 * a2 * b * c2 * d2 * e * g * n6 -
                            10080 * a * b * c * d3 * e * g * n6 -
                            5040 * a3 * b * c * d3 * e * g * n6 - 1260 * b * d4 * e * g * n6 -
                            5040 * a2 * b * d4 * e * g * n6 - 1260 * a4 * b * d4 * e * g * n6 -
                            3360 * c3 * e * f * g * n6 - 5040 * b2 * c3 * e * f * g * n6 -
                            10080 * a * c2 * d * e * f * g * n6 -
                            15120 * a * b2 * c2 * d * e * f * g * n6 -
                            10080 * c * d2 * e * f * g * n6 - 10080 * a2 * c * d2 * e * f * g * n6 -
                            10080 * b2 * c * d2 * e * f * g * n6 -
                            15120 * a2 * b2 * c * d2 * e * f * g * n6 -
                            10080 * a * d3 * e * f * g * n6 - 3360 * a3 * d3 * e * f * g * n6 -
                            10080 * a * b2 * d3 * e * f * g * n6 -
                            5040 * a3 * b2 * d3 * e * f * g * n6 -
                            12600 * b * c2 * e * f2 * g * n6 - 7560 * b3 * c2 * e * f2 * g * n6 -
                            25200 * a * b * c * d * e * f2 * g * n6 -
                            15120 * a * b3 * c * d * e * f2 * g * n6 -
                            12600 * b * d2 * e * f2 * g * n6 -
                            12600 * a2 * b * d2 * e * f2 * g * n6 -
                            5040 * b3 * d2 * e * f2 * g * n6 -
                            7560 * a2 * b3 * d2 * e * f2 * g * n6 -
                            15120 * b2 * c * e * f3 * g * n6 - 5040 * b4 * c * e * f3 * g * n6 -
                            15120 * a * b2 * d * e * f3 * g * n6 -
                            5040 * a * b4 * d * e * f3 * g * n6 - 5880 * b3 * e * f4 * g * n6 -
                            1260 * b5 * e * f4 * g * n6 - 2520 * b2 * c3 * d * g2 * n6 -
                            3780 * a * b2 * c2 * d2 * g2 * n6 - 1680 * b2 * c * d3 * g2 * n6 -
                            2520 * a2 * b2 * c * d3 * g2 * n6 - 1260 * a * b2 * d4 * g2 * n6 -
                            630 * a3 * b2 * d4 * g2 * n6 - 10080 * b * c2 * d * f * g2 * n6 -
                            7560 * b3 * c2 * d * f * g2 * n6 -
                            10080 * a * b * c * d2 * f * g2 * n6 -
                            7560 * a * b3 * c * d2 * f * g2 * n6 - 3360 * b * d3 * f * g2 * n6 -
                            3360 * a2 * b * d3 * f * g2 * n6 - 1680 * b3 * d3 * f * g2 * n6 -
                            2520 * a2 * b3 * d3 * f * g2 * n6 - 22680 * b2 * c * d * f2 * g2 * n6 -
                            7560 * b4 * c * d * f2 * g2 * n6 - 11340 * a * b2 * d2 * f2 * g2 * n6 -
                            3780 * a * b4 * d2 * f2 * g2 * n6 - 12600 * b3 * d * f3 * g2 * n6 -
                            2520 * b5 * d * f3 * g2 * n6 + 1260 * a * c2 * e2 * n8 +
                            2520 * c * d * e2 * n8 + 2520 * a2 * c * d * e2 * n8 +
                            3780 * a * d2 * e2 * n8 + 1260 * a3 * d2 * e2 * n8 +
                            2520 * a * b * c * e2 * f * n8 + 2520 * b * d * e2 * f * n8 +
                            2520 * a2 * b * d * e2 * f * n8 + 630 * a * e2 * f2 * n8 +
                            1260 * a * b2 * e2 * f2 * n8 + 2520 * b * c2 * e * g * n8 +
                            5040 * a * b * c * d * e * g * n8 + 2520 * b * d2 * e * g * n8 +
                            2520 * a2 * b * d2 * e * g * n8 + 5040 * c * e * f * g * n8 +
                            5040 * b2 * c * e * f * g * n8 + 5040 * a * d * e * f * g * n8 +
                            5040 * a * b2 * d * e * f * g * n8 + 6300 * b * e * f2 * g * n8 +
                            2520 * b3 * e * f2 * g * n8 + 2520 * b2 * c * d * g2 * n8 +
                            1260 * a * b2 * d2 * g2 * n8 + 5040 * b * d * f * g2 * n8 +
                            2520 * b3 * d * f * g2 * n8 - 630 * a * e2 * n10 -
                            1260 * b * e * g * n10)) /
                        3 +
                    5 * (126 * c5 * d2 * e3 * f2 + 1050 * a * c4 * d3 * e3 * f2 +
                            3150 * a2 * c3 * d4 * e3 * f2 + 4410 * a3 * c2 * d5 * e3 * f2 +
                            2940 * a4 * c * d6 * e3 * f2 + 756 * a5 * d7 * e3 * f2 +
                            630 * b * c4 * d2 * e3 * f3 + 4200 * a * b * c3 * d3 * e3 * f3 +
                            9450 * a2 * b * c2 * d4 * e3 * f3 + 8820 * a3 * b * c * d5 * e3 * f3 +
                            2940 * a4 * b * d6 * e3 * f3 + 1260 * b2 * c3 * d2 * e3 * f4 +
                            6300 * a * b2 * c2 * d3 * e3 * f4 + 9450 * a2 * b2 * c * d4 * e3 * f4 +
                            4410 * a3 * b2 * d5 * e3 * f4 + 1260 * b3 * c2 * d2 * e3 * f5 +
                            4200 * a * b3 * c * d3 * e3 * f5 + 3150 * a2 * b3 * d4 * e3 * f5 +
                            630 * b4 * c * d2 * e3 * f6 + 1050 * a * b4 * d3 * e3 * f6 +
                            126 * b5 * d2 * e3 * f7 + 504 * c5 * d3 * e2 * f * g +
                            3150 * a * c4 * d4 * e2 * f * g + 7560 * a2 * c3 * d5 * e2 * f * g +
                            8820 * a3 * c2 * d6 * e2 * f * g + 5040 * a4 * c * d7 * e2 * f * g +
                            1134 * a5 * d8 * e2 * f * g + 3150 * b * c4 * d3 * e2 * f2 * g +
                            15750 * a * b * c3 * d4 * e2 * f2 * g +
                            28350 * a2 * b * c2 * d5 * e2 * f2 * g +
                            22050 * a3 * b * c * d6 * e2 * f2 * g +
                            6300 * a4 * b * d7 * e2 * f2 * g + 7560 * b2 * c3 * d3 * e2 * f3 * g +
                            28350 * a * b2 * c2 * d4 * e2 * f3 * g +
                            34020 * a2 * b2 * c * d5 * e2 * f3 * g +
                            13230 * a3 * b2 * d6 * e2 * f3 * g + 8820 * b3 * c2 * d3 * e2 * f4 * g +
                            22050 * a * b3 * c * d4 * e2 * f4 * g +
                            13230 * a2 * b3 * d5 * e2 * f4 * g + 5040 * b4 * c * d3 * e2 * f5 * g +
                            6300 * a * b4 * d4 * e2 * f5 * g + 1134 * b5 * d3 * e2 * f6 * g +
                            1260 * b * c4 * d4 * e * f * g2 + 5040 * a * b * c3 * d5 * e * f * g2 +
                            7560 * a2 * b * c2 * d6 * e * f * g2 +
                            5040 * a3 * b * c * d7 * e * f * g2 + 1260 * a4 * b * d8 * e * f * g2 +
                            5670 * b2 * c3 * d4 * e * f2 * g2 +
                            17010 * a * b2 * c2 * d5 * e * f2 * g2 +
                            17010 * a2 * b2 * c * d6 * e * f2 * g2 +
                            5670 * a3 * b2 * d7 * e * f2 * g2 + 9450 * b3 * c2 * d4 * e * f3 * g2 +
                            18900 * a * b3 * c * d5 * e * f3 * g2 +
                            9450 * a2 * b3 * d6 * e * f3 * g2 + 6930 * b4 * c * d4 * e * f4 * g2 +
                            6930 * a * b4 * d5 * e * f4 * g2 + 1890 * b5 * d4 * e * f5 * g2 +
                            504 * b2 * c3 * d5 * f * g3 + 1260 * a * b2 * c2 * d6 * f * g3 +
                            1080 * a2 * b2 * c * d7 * f * g3 + 315 * a3 * b2 * d8 * f * g3 +
                            1638 * b3 * c2 * d5 * f2 * g3 + 2730 * a * b3 * c * d6 * f2 * g3 +
                            1170 * a2 * b3 * d7 * f2 * g3 + 1764 * b4 * c * d5 * f3 * g3 +
                            1470 * a * b4 * d6 * f3 * g3 + 630 * b5 * d5 * f4 * g3 -
                            126 * c5 * d2 * e3 * n2 - 1050 * a * c4 * d3 * e3 * n2 -
                            3150 * a2 * c3 * d4 * e3 * n2 - 4410 * a3 * c2 * d5 * e3 * n2 -
                            2940 * a4 * c * d6 * e3 * n2 - 756 * a5 * d7 * e3 * n2 -
                            630 * b * c4 * d2 * e3 * f * n2 - 4200 * a * b * c3 * d3 * e3 * f * n2 -
                            9450 * a2 * b * c2 * d4 * e3 * f * n2 -
                            8820 * a3 * b * c * d5 * e3 * f * n2 -
                            2940 * a4 * b * d6 * e3 * f * n2 - 42 * c5 * e3 * f2 * n2 -
                            630 * a * c4 * d * e3 * f2 * n2 - 420 * c3 * d2 * e3 * f2 * n2 -
                            2520 * a2 * c3 * d2 * e3 * f2 * n2 -
                            1260 * b2 * c3 * d2 * e3 * f2 * n2 - 2100 * a * c2 * d3 * e3 * f2 * n2 -
                            4200 * a3 * c2 * d3 * e3 * f2 * n2 -
                            6300 * a * b2 * c2 * d3 * e3 * f2 * n2 -
                            3150 * a2 * c * d4 * e3 * f2 * n2 - 3150 * a4 * c * d4 * e3 * f2 * n2 -
                            9450 * a2 * b2 * c * d4 * e3 * f2 * n2 - 1470 * a3 * d5 * e3 * f2 * n2 -
                            882 * a5 * d5 * e3 * f2 * n2 - 4410 * a3 * b2 * d5 * e3 * f2 * n2 -
                            210 * b * c4 * e3 * f3 * n2 - 2520 * a * b * c3 * d * e3 * f3 * n2 -
                            1260 * b * c2 * d2 * e3 * f3 * n2 -
                            7560 * a2 * b * c2 * d2 * e3 * f3 * n2 -
                            1260 * b3 * c2 * d2 * e3 * f3 * n2 -
                            4200 * a * b * c * d3 * e3 * f3 * n2 -
                            8400 * a3 * b * c * d3 * e3 * f3 * n2 -
                            4200 * a * b3 * c * d3 * e3 * f3 * n2 -
                            3150 * a2 * b * d4 * e3 * f3 * n2 - 3150 * a4 * b * d4 * e3 * f3 * n2 -
                            3150 * a2 * b3 * d4 * e3 * f3 * n2 - 420 * b2 * c3 * e3 * f4 * n2 -
                            3780 * a * b2 * c2 * d * e3 * f4 * n2 -
                            1260 * b2 * c * d2 * e3 * f4 * n2 -
                            7560 * a2 * b2 * c * d2 * e3 * f4 * n2 -
                            630 * b4 * c * d2 * e3 * f4 * n2 - 2100 * a * b2 * d3 * e3 * f4 * n2 -
                            4200 * a3 * b2 * d3 * e3 * f4 * n2 - 1050 * a * b4 * d3 * e3 * f4 * n2 -
                            420 * b3 * c2 * e3 * f5 * n2 - 2520 * a * b3 * c * d * e3 * f5 * n2 -
                            420 * b3 * d2 * e3 * f5 * n2 - 2520 * a2 * b3 * d2 * e3 * f5 * n2 -
                            126 * b5 * d2 * e3 * f5 * n2 - 210 * b4 * c * e3 * f6 * n2 -
                            630 * a * b4 * d * e3 * f6 * n2 - 42 * b5 * e3 * f7 * n2 -
                            630 * b * c4 * d3 * e2 * g * n2 - 3150 * a * b * c3 * d4 * e2 * g * n2 -
                            5670 * a2 * b * c2 * d5 * e2 * g * n2 -
                            4410 * a3 * b * c * d6 * e2 * g * n2 -
                            1260 * a4 * b * d7 * e2 * g * n2 - 504 * c5 * d * e2 * f * g * n2 -
                            3780 * a * c4 * d2 * e2 * f * g * n2 -
                            1680 * c3 * d3 * e2 * f * g * n2 -
                            10080 * a2 * c3 * d3 * e2 * f * g * n2 -
                            2520 * b2 * c3 * d3 * e2 * f * g * n2 -
                            6300 * a * c2 * d4 * e2 * f * g * n2 -
                            12600 * a3 * c2 * d4 * e2 * f * g * n2 -
                            9450 * a * b2 * c2 * d4 * e2 * f * g * n2 -
                            7560 * a2 * c * d5 * e2 * f * g * n2 -
                            7560 * a4 * c * d5 * e2 * f * g * n2 -
                            11340 * a2 * b2 * c * d5 * e2 * f * g * n2 -
                            2940 * a3 * d6 * e2 * f * g * n2 - 1764 * a5 * d6 * e2 * f * g * n2 -
                            4410 * a3 * b2 * d6 * e2 * f * g * n2 -
                            3150 * b * c4 * d * e2 * f2 * g * n2 -
                            18900 * a * b * c3 * d2 * e2 * f2 * g * n2 -
                            6300 * b * c2 * d3 * e2 * f2 * g * n2 -
                            37800 * a2 * b * c2 * d3 * e2 * f2 * g * n2 -
                            3780 * b3 * c2 * d3 * e2 * f2 * g * n2 -
                            15750 * a * b * c * d4 * e2 * f2 * g * n2 -
                            31500 * a3 * b * c * d4 * e2 * f2 * g * n2 -
                            9450 * a * b3 * c * d4 * e2 * f2 * g * n2 -
                            9450 * a2 * b * d5 * e2 * f2 * g * n2 -
                            9450 * a4 * b * d5 * e2 * f2 * g * n2 -
                            5670 * a2 * b3 * d5 * e2 * f2 * g * n2 -
                            7560 * b2 * c3 * d * e2 * f3 * g * n2 -
                            34020 * a * b2 * c2 * d2 * e2 * f3 * g * n2 -
                            7560 * b2 * c * d3 * e2 * f3 * g * n2 -
                            45360 * a2 * b2 * c * d3 * e2 * f3 * g * n2 -
                            2520 * b4 * c * d3 * e2 * f3 * g * n2 -
                            9450 * a * b2 * d4 * e2 * f3 * g * n2 -
                            18900 * a3 * b2 * d4 * e2 * f3 * g * n2 -
                            3150 * a * b4 * d4 * e2 * f3 * g * n2 -
                            8820 * b3 * c2 * d * e2 * f4 * g * n2 -
                            26460 * a * b3 * c * d2 * e2 * f4 * g * n2 -
                            2940 * b3 * d3 * e2 * f4 * g * n2 -
                            17640 * a2 * b3 * d3 * e2 * f4 * g * n2 -
                            630 * b5 * d3 * e2 * f4 * g * n2 -
                            5040 * b4 * c * d * e2 * f5 * g * n2 -
                            7560 * a * b4 * d2 * e2 * f5 * g * n2 -
                            1134 * b5 * d * e2 * f6 * g * n2 - 630 * b2 * c3 * d4 * e * g2 * n2 -
                            1890 * a * b2 * c2 * d5 * e * g2 * n2 -
                            1890 * a2 * b2 * c * d6 * e * g2 * n2 -
                            630 * a3 * b2 * d7 * e * g2 * n2 -
                            2520 * b * c4 * d2 * e * f * g2 * n2 -
                            10080 * a * b * c3 * d3 * e * f * g2 * n2 -
                            2520 * b * c2 * d4 * e * f * g2 * n2 -
                            15120 * a2 * b * c2 * d4 * e * f * g2 * n2 -
                            1890 * b3 * c2 * d4 * e * f * g2 * n2 -
                            5040 * a * b * c * d5 * e * f * g2 * n2 -
                            10080 * a3 * b * c * d5 * e * f * g2 * n2 -
                            3780 * a * b3 * c * d5 * e * f * g2 * n2 -
                            2520 * a2 * b * d6 * e * f * g2 * n2 -
                            2520 * a4 * b * d6 * e * f * g2 * n2 -
                            1890 * a2 * b3 * d6 * e * f * g2 * n2 -
                            11340 * b2 * c3 * d2 * e * f2 * g2 * n2 -
                            34020 * a * b2 * c2 * d3 * e * f2 * g2 * n2 -
                            5670 * b2 * c * d4 * e * f2 * g2 * n2 -
                            34020 * a2 * b2 * c * d4 * e * f2 * g2 * n2 -
                            1890 * b4 * c * d4 * e * f2 * g2 * n2 -
                            5670 * a * b2 * d5 * e * f2 * g2 * n2 -
                            11340 * a3 * b2 * d5 * e * f2 * g2 * n2 -
                            1890 * a * b4 * d5 * e * f2 * g2 * n2 -
                            18900 * b3 * c2 * d2 * e * f3 * g2 * n2 -
                            37800 * a * b3 * c * d3 * e * f3 * g2 * n2 -
                            3150 * b3 * d4 * e * f3 * g2 * n2 -
                            18900 * a2 * b3 * d4 * e * f3 * g2 * n2 -
                            630 * b5 * d4 * e * f3 * g2 * n2 -
                            13860 * b4 * c * d2 * e * f4 * g2 * n2 -
                            13860 * a * b4 * d3 * e * f4 * g2 * n2 -
                            3780 * b5 * d2 * e * f5 * g2 * n2 - 126 * b3 * c2 * d5 * g3 * n2 -
                            210 * a * b3 * c * d6 * g3 * n2 - 90 * a2 * b3 * d7 * g3 * n2 -
                            1680 * b2 * c3 * d3 * f * g3 * n2 -
                            3780 * a * b2 * c2 * d4 * f * g3 * n2 -
                            504 * b2 * c * d5 * f * g3 * n2 -
                            3024 * a2 * b2 * c * d5 * f * g3 * n2 -
                            252 * b4 * c * d5 * f * g3 * n2 - 420 * a * b2 * d6 * f * g3 * n2 -
                            840 * a3 * b2 * d6 * f * g3 * n2 - 210 * a * b4 * d6 * f * g3 * n2 -
                            5460 * b3 * c2 * d3 * f2 * g3 * n2 -
                            8190 * a * b3 * c * d4 * f2 * g3 * n2 - 546 * b3 * d5 * f2 * g3 * n2 -
                            3276 * a2 * b3 * d5 * f2 * g3 * n2 - 126 * b5 * d5 * f2 * g3 * n2 -
                            5880 * b4 * c * d3 * f3 * g3 * n2 - 4410 * a * b4 * d4 * f3 * g3 * n2 -
                            2100 * b5 * d3 * f4 * g3 * n2 + 42 * c5 * e3 * n4 +
                            630 * a * c4 * d * e3 * n4 + 420 * c3 * d2 * e3 * n4 +
                            2520 * a2 * c3 * d2 * e3 * n4 + 2100 * a * c2 * d3 * e3 * n4 +
                            4200 * a3 * c2 * d3 * e3 * n4 + 3150 * a2 * c * d4 * e3 * n4 +
                            3150 * a4 * c * d4 * e3 * n4 + 1470 * a3 * d5 * e3 * n4 +
                            882 * a5 * d5 * e3 * n4 + 210 * b * c4 * e3 * f * n4 +
                            2520 * a * b * c3 * d * e3 * f * n4 + 1260 * b * c2 * d2 * e3 * f * n4 +
                            7560 * a2 * b * c2 * d2 * e3 * f * n4 +
                            4200 * a * b * c * d3 * e3 * f * n4 +
                            8400 * a3 * b * c * d3 * e3 * f * n4 +
                            3150 * a2 * b * d4 * e3 * f * n4 + 3150 * a4 * b * d4 * e3 * f * n4 +
                            140 * c3 * e3 * f2 * n4 + 210 * a2 * c3 * e3 * f2 * n4 +
                            420 * b2 * c3 * e3 * f2 * n4 + 1260 * a * c2 * d * e3 * f2 * n4 +
                            630 * a3 * c2 * d * e3 * f2 * n4 +
                            3780 * a * b2 * c2 * d * e3 * f2 * n4 + 630 * c * d2 * e3 * f2 * n4 +
                            2520 * a2 * c * d2 * e3 * f2 * n4 + 630 * a4 * c * d2 * e3 * f2 * n4 +
                            1260 * b2 * c * d2 * e3 * f2 * n4 +
                            7560 * a2 * b2 * c * d2 * e3 * f2 * n4 + 1050 * a * d3 * e3 * f2 * n4 +
                            1400 * a3 * d3 * e3 * f2 * n4 + 210 * a5 * d3 * e3 * f2 * n4 +
                            2100 * a * b2 * d3 * e3 * f2 * n4 + 4200 * a3 * b2 * d3 * e3 * f2 * n4 +
                            420 * b * c2 * e3 * f3 * n4 + 630 * a2 * b * c2 * e3 * f3 * n4 +
                            420 * b3 * c2 * e3 * f3 * n4 + 2520 * a * b * c * d * e3 * f3 * n4 +
                            1260 * a3 * b * c * d * e3 * f3 * n4 +
                            2520 * a * b3 * c * d * e3 * f3 * n4 + 630 * b * d2 * e3 * f3 * n4 +
                            2520 * a2 * b * d2 * e3 * f3 * n4 + 630 * a4 * b * d2 * e3 * f3 * n4 +
                            420 * b3 * d2 * e3 * f3 * n4 + 2520 * a2 * b3 * d2 * e3 * f3 * n4 +
                            420 * b2 * c * e3 * f4 * n4 + 630 * a2 * b2 * c * e3 * f4 * n4 +
                            210 * b4 * c * e3 * f4 * n4 + 1260 * a * b2 * d * e3 * f4 * n4 +
                            630 * a3 * b2 * d * e3 * f4 * n4 + 630 * a * b4 * d * e3 * f4 * n4 +
                            140 * b3 * e3 * f5 * n4 + 210 * a2 * b3 * e3 * f5 * n4 +
                            42 * b5 * e3 * f5 * n4 + 630 * b * c4 * d * e2 * g * n4 +
                            3780 * a * b * c3 * d2 * e2 * g * n4 +
                            1260 * b * c2 * d3 * e2 * g * n4 +
                            7560 * a2 * b * c2 * d3 * e2 * g * n4 +
                            3150 * a * b * c * d4 * e2 * g * n4 +
                            6300 * a3 * b * c * d4 * e2 * g * n4 +
                            1890 * a2 * b * d5 * e2 * g * n4 + 1890 * a4 * b * d5 * e2 * g * n4 +
                            630 * a * c4 * e2 * f * g * n4 + 1680 * c3 * d * e2 * f * g * n4 +
                            2520 * a2 * c3 * d * e2 * f * g * n4 +
                            2520 * b2 * c3 * d * e2 * f * g * n4 +
                            7560 * a * c2 * d2 * e2 * f * g * n4 +
                            3780 * a3 * c2 * d2 * e2 * f * g * n4 +
                            11340 * a * b2 * c2 * d2 * e2 * f * g * n4 +
                            2520 * c * d3 * e2 * f * g * n4 +
                            10080 * a2 * c * d3 * e2 * f * g * n4 +
                            2520 * a4 * c * d3 * e2 * f * g * n4 +
                            2520 * b2 * c * d3 * e2 * f * g * n4 +
                            15120 * a2 * b2 * c * d3 * e2 * f * g * n4 +
                            3150 * a * d4 * e2 * f * g * n4 + 4200 * a3 * d4 * e2 * f * g * n4 +
                            630 * a5 * d4 * e2 * f * g * n4 + 3150 * a * b2 * d4 * e2 * f * g * n4 +
                            6300 * a3 * b2 * d4 * e2 * f * g * n4 +
                            3150 * a * b * c3 * e2 * f2 * g * n4 +
                            6300 * b * c2 * d * e2 * f2 * g * n4 +
                            9450 * a2 * b * c2 * d * e2 * f2 * g * n4 +
                            3780 * b3 * c2 * d * e2 * f2 * g * n4 +
                            18900 * a * b * c * d2 * e2 * f2 * g * n4 +
                            9450 * a3 * b * c * d2 * e2 * f2 * g * n4 +
                            11340 * a * b3 * c * d2 * e2 * f2 * g * n4 +
                            3150 * b * d3 * e2 * f2 * g * n4 +
                            12600 * a2 * b * d3 * e2 * f2 * g * n4 +
                            3150 * a4 * b * d3 * e2 * f2 * g * n4 +
                            1260 * b3 * d3 * e2 * f2 * g * n4 +
                            7560 * a2 * b3 * d3 * e2 * f2 * g * n4 +
                            5670 * a * b2 * c2 * e2 * f3 * g * n4 +
                            7560 * b2 * c * d * e2 * f3 * g * n4 +
                            11340 * a2 * b2 * c * d * e2 * f3 * g * n4 +
                            2520 * b4 * c * d * e2 * f3 * g * n4 +
                            11340 * a * b2 * d2 * e2 * f3 * g * n4 +
                            5670 * a3 * b2 * d2 * e2 * f3 * g * n4 +
                            3780 * a * b4 * d2 * e2 * f3 * g * n4 +
                            4410 * a * b3 * c * e2 * f4 * g * n4 +
                            2940 * b3 * d * e2 * f4 * g * n4 +
                            4410 * a2 * b3 * d * e2 * f4 * g * n4 +
                            630 * b5 * d * e2 * f4 * g * n4 + 1260 * a * b4 * e2 * f5 * g * n4 +
                            1260 * b2 * c3 * d2 * e * g2 * n4 +
                            3780 * a * b2 * c2 * d3 * e * g2 * n4 +
                            630 * b2 * c * d4 * e * g2 * n4 +
                            3780 * a2 * b2 * c * d4 * e * g2 * n4 +
                            630 * a * b2 * d5 * e * g2 * n4 + 1260 * a3 * b2 * d5 * e * g2 * n4 +
                            1260 * b * c4 * e * f * g2 * n4 +
                            5040 * a * b * c3 * d * e * f * g2 * n4 +
                            5040 * b * c2 * d2 * e * f * g2 * n4 +
                            7560 * a2 * b * c2 * d2 * e * f * g2 * n4 +
                            3780 * b3 * c2 * d2 * e * f * g2 * n4 +
                            10080 * a * b * c * d3 * e * f * g2 * n4 +
                            5040 * a3 * b * c * d3 * e * f * g2 * n4 +
                            7560 * a * b3 * c * d3 * e * f * g2 * n4 +
                            1260 * b * d4 * e * f * g2 * n4 + 5040 * a2 * b * d4 * e * f * g2 * n4 +
                            1260 * a4 * b * d4 * e * f * g2 * n4 + 630 * b3 * d4 * e * f * g2 * n4 +
                            3780 * a2 * b3 * d4 * e * f * g2 * n4 +
                            5670 * b2 * c3 * e * f2 * g2 * n4 +
                            17010 * a * b2 * c2 * d * e * f2 * g2 * n4 +
                            11340 * b2 * c * d2 * e * f2 * g2 * n4 +
                            17010 * a2 * b2 * c * d2 * e * f2 * g2 * n4 +
                            3780 * b4 * c * d2 * e * f2 * g2 * n4 +
                            11340 * a * b2 * d3 * e * f2 * g2 * n4 +
                            5670 * a3 * b2 * d3 * e * f2 * g2 * n4 +
                            3780 * a * b4 * d3 * e * f2 * g2 * n4 +
                            9450 * b3 * c2 * e * f3 * g2 * n4 +
                            18900 * a * b3 * c * d * e * f3 * g2 * n4 +
                            6300 * b3 * d2 * e * f3 * g2 * n4 +
                            9450 * a2 * b3 * d2 * e * f3 * g2 * n4 +
                            1260 * b5 * d2 * e * f3 * g2 * n4 + 6930 * b4 * c * e * f4 * g2 * n4 +
                            6930 * a * b4 * d * e * f4 * g2 * n4 + 1890 * b5 * e * f5 * g2 * n4 +
                            420 * b3 * c2 * d3 * g3 * n4 + 630 * a * b3 * c * d4 * g3 * n4 +
                            42 * b3 * d5 * g3 * n4 + 252 * a2 * b3 * d5 * g3 * n4 +
                            2520 * b2 * c3 * d * f * g3 * n4 +
                            3780 * a * b2 * c2 * d2 * f * g3 * n4 +
                            1680 * b2 * c * d3 * f * g3 * n4 +
                            2520 * a2 * b2 * c * d3 * f * g3 * n4 +
                            840 * b4 * c * d3 * f * g3 * n4 + 1260 * a * b2 * d4 * f * g3 * n4 +
                            630 * a3 * b2 * d4 * f * g3 * n4 + 630 * a * b4 * d4 * f * g3 * n4 +
                            8190 * b3 * c2 * d * f2 * g3 * n4 +
                            8190 * a * b3 * c * d2 * f2 * g3 * n4 + 1820 * b3 * d3 * f2 * g3 * n4 +
                            2730 * a2 * b3 * d3 * f2 * g3 * n4 + 420 * b5 * d3 * f2 * g3 * n4 +
                            8820 * b4 * c * d * f3 * g3 * n4 + 4410 * a * b4 * d2 * f3 * g3 * n4 +
                            3150 * b5 * d * f4 * g3 * n4 - 140 * c3 * e3 * n6 -
                            210 * a2 * c3 * e3 * n6 - 1260 * a * c2 * d * e3 * n6 -
                            630 * a3 * c2 * d * e3 * n6 - 630 * c * d2 * e3 * n6 -
                            2520 * a2 * c * d2 * e3 * n6 - 630 * a4 * c * d2 * e3 * n6 -
                            1050 * a * d3 * e3 * n6 - 1400 * a3 * d3 * e3 * n6 -
                            210 * a5 * d3 * e3 * n6 - 420 * b * c2 * e3 * f * n6 -
                            630 * a2 * b * c2 * e3 * f * n6 - 2520 * a * b * c * d * e3 * f * n6 -
                            1260 * a3 * b * c * d * e3 * f * n6 - 630 * b * d2 * e3 * f * n6 -
                            2520 * a2 * b * d2 * e3 * f * n6 - 630 * a4 * b * d2 * e3 * f * n6 -
                            210 * c * e3 * f2 * n6 - 210 * a2 * c * e3 * f2 * n6 -
                            420 * b2 * c * e3 * f2 * n6 - 630 * a2 * b2 * c * e3 * f2 * n6 -
                            630 * a * d * e3 * f2 * n6 - 210 * a3 * d * e3 * f2 * n6 -
                            1260 * a * b2 * d * e3 * f2 * n6 - 630 * a3 * b2 * d * e3 * f2 * n6 -
                            210 * b * e3 * f3 * n6 - 210 * a2 * b * e3 * f3 * n6 -
                            140 * b3 * e3 * f3 * n6 - 210 * a2 * b3 * e3 * f3 * n6 -
                            630 * a * b * c3 * e2 * g * n6 - 1260 * b * c2 * d * e2 * g * n6 -
                            1890 * a2 * b * c2 * d * e2 * g * n6 -
                            3780 * a * b * c * d2 * e2 * g * n6 -
                            1890 * a3 * b * c * d2 * e2 * g * n6 - 630 * b * d3 * e2 * g * n6 -
                            2520 * a2 * b * d3 * e2 * g * n6 - 630 * a4 * b * d3 * e2 * g * n6 -
                            1260 * a * c2 * e2 * f * g * n6 - 1890 * a * b2 * c2 * e2 * f * g * n6 -
                            2520 * c * d * e2 * f * g * n6 - 2520 * a2 * c * d * e2 * f * g * n6 -
                            2520 * b2 * c * d * e2 * f * g * n6 -
                            3780 * a2 * b2 * c * d * e2 * f * g * n6 -
                            3780 * a * d2 * e2 * f * g * n6 - 1260 * a3 * d2 * e2 * f * g * n6 -
                            3780 * a * b2 * d2 * e2 * f * g * n6 -
                            1890 * a3 * b2 * d2 * e2 * f * g * n6 -
                            3150 * a * b * c * e2 * f2 * g * n6 -
                            1890 * a * b3 * c * e2 * f2 * g * n6 - 3150 * b * d * e2 * f2 * g * n6 -
                            3150 * a2 * b * d * e2 * f2 * g * n6 -
                            1260 * b3 * d * e2 * f2 * g * n6 -
                            1890 * a2 * b3 * d * e2 * f2 * g * n6 -
                            1890 * a * b2 * e2 * f3 * g * n6 - 630 * a * b4 * e2 * f3 * g * n6 -
                            630 * b2 * c3 * e * g2 * n6 - 1890 * a * b2 * c2 * d * e * g2 * n6 -
                            1260 * b2 * c * d2 * e * g2 * n6 -
                            1890 * a2 * b2 * c * d2 * e * g2 * n6 -
                            1260 * a * b2 * d3 * e * g2 * n6 - 630 * a3 * b2 * d3 * e * g2 * n6 -
                            2520 * b * c2 * e * f * g2 * n6 - 1890 * b3 * c2 * e * f * g2 * n6 -
                            5040 * a * b * c * d * e * f * g2 * n6 -
                            3780 * a * b3 * c * d * e * f * g2 * n6 -
                            2520 * b * d2 * e * f * g2 * n6 - 2520 * a2 * b * d2 * e * f * g2 * n6 -
                            1260 * b3 * d2 * e * f * g2 * n6 -
                            1890 * a2 * b3 * d2 * e * f * g2 * n6 -
                            5670 * b2 * c * e * f2 * g2 * n6 - 1890 * b4 * c * e * f2 * g2 * n6 -
                            5670 * a * b2 * d * e * f2 * g2 * n6 -
                            1890 * a * b4 * d * e * f2 * g2 * n6 - 3150 * b3 * e * f3 * g2 * n6 -
                            630 * b5 * e * f3 * g2 * n6 - 630 * b3 * c2 * d * g3 * n6 -
                            630 * a * b3 * c * d2 * g3 * n6 - 140 * b3 * d3 * g3 * n6 -
                            210 * a2 * b3 * d3 * g3 * n6 - 2520 * b2 * c * d * f * g3 * n6 -
                            1260 * b4 * c * d * f * g3 * n6 - 1260 * a * b2 * d2 * f * g3 * n6 -
                            630 * a * b4 * d2 * f * g3 * n6 - 2730 * b3 * d * f2 * g3 * n6 -
                            630 * b5 * d * f2 * g3 * n6 + 210 * c * e3 * n8 +
                            210 * a2 * c * e3 * n8 + 630 * a * d * e3 * n8 +
                            210 * a3 * d * e3 * n8 + 210 * b * e3 * f * n8 +
                            210 * a2 * b * e3 * f * n8 + 630 * a * b * c * e2 * g * n8 +
                            630 * b * d * e2 * g * n8 + 630 * a2 * b * d * e2 * g * n8 +
                            630 * a * e2 * f * g * n8 + 630 * a * b2 * e2 * f * g * n8 +
                            630 * b2 * c * e * g2 * n8 + 630 * a * b2 * d * e * g2 * n8 +
                            1260 * b * e * f * g2 * n8 + 630 * b3 * e * f * g2 * n8 +
                            210 * b3 * d * g3 * n8) +
                    2 * (126 * c5 * d * e4 * f2 + 1575 * a * c4 * d2 * e4 * f2 +
                            6300 * a2 * c3 * d3 * e4 * f2 + 11025 * a3 * c2 * d4 * e4 * f2 +
                            8820 * a4 * c * d5 * e4 * f2 + 2646 * a5 * d6 * e4 * f2 +
                            630 * b * c4 * d * e4 * f3 + 6300 * a * b * c3 * d2 * e4 * f3 +
                            18900 * a2 * b * c2 * d3 * e4 * f3 + 22050 * a3 * b * c * d4 * e4 * f3 +
                            8820 * a4 * b * d5 * e4 * f3 + 1260 * b2 * c3 * d * e4 * f4 +
                            9450 * a * b2 * c2 * d2 * e4 * f4 + 18900 * a2 * b2 * c * d3 * e4 * f4 +
                            11025 * a3 * b2 * d4 * e4 * f4 + 1260 * b3 * c2 * d * e4 * f5 +
                            6300 * a * b3 * c * d2 * e4 * f5 + 6300 * a2 * b3 * d3 * e4 * f5 +
                            630 * b4 * c * d * e4 * f6 + 1575 * a * b4 * d2 * e4 * f6 +
                            126 * b5 * d * e4 * f7 + 1008 * c5 * d2 * e3 * f * g +
                            8400 * a * c4 * d3 * e3 * f * g + 25200 * a2 * c3 * d4 * e3 * f * g +
                            35280 * a3 * c2 * d5 * e3 * f * g + 23520 * a4 * c * d6 * e3 * f * g +
                            6048 * a5 * d7 * e3 * f * g + 6300 * b * c4 * d2 * e3 * f2 * g +
                            42000 * a * b * c3 * d3 * e3 * f2 * g +
                            94500 * a2 * b * c2 * d4 * e3 * f2 * g +
                            88200 * a3 * b * c * d5 * e3 * f2 * g +
                            29400 * a4 * b * d6 * e3 * f2 * g + 15120 * b2 * c3 * d2 * e3 * f3 * g +
                            75600 * a * b2 * c2 * d3 * e3 * f3 * g +
                            113400 * a2 * b2 * c * d4 * e3 * f3 * g +
                            52920 * a3 * b2 * d5 * e3 * f3 * g +
                            17640 * b3 * c2 * d2 * e3 * f4 * g +
                            58800 * a * b3 * c * d3 * e3 * f4 * g +
                            44100 * a2 * b3 * d4 * e3 * f4 * g + 10080 * b4 * c * d2 * e3 * f5 * g +
                            16800 * a * b4 * d3 * e3 * f5 * g + 2268 * b5 * d2 * e3 * f6 * g +
                            5040 * b * c4 * d3 * e2 * f * g2 +
                            25200 * a * b * c3 * d4 * e2 * f * g2 +
                            45360 * a2 * b * c2 * d5 * e2 * f * g2 +
                            35280 * a3 * b * c * d6 * e2 * f * g2 +
                            10080 * a4 * b * d7 * e2 * f * g2 +
                            22680 * b2 * c3 * d3 * e2 * f2 * g2 +
                            85050 * a * b2 * c2 * d4 * e2 * f2 * g2 +
                            102060 * a2 * b2 * c * d5 * e2 * f2 * g2 +
                            39690 * a3 * b2 * d6 * e2 * f2 * g2 +
                            37800 * b3 * c2 * d3 * e2 * f3 * g2 +
                            94500 * a * b3 * c * d4 * e2 * f3 * g2 +
                            56700 * a2 * b3 * d5 * e2 * f3 * g2 +
                            27720 * b4 * c * d3 * e2 * f4 * g2 +
                            34650 * a * b4 * d4 * e2 * f4 * g2 + 7560 * b5 * d3 * e2 * f5 * g2 +
                            5040 * b2 * c3 * d4 * e * f * g3 +
                            15120 * a * b2 * c2 * d5 * e * f * g3 +
                            15120 * a2 * b2 * c * d6 * e * f * g3 +
                            5040 * a3 * b2 * d7 * e * f * g3 + 16380 * b3 * c2 * d4 * e * f2 * g3 +
                            32760 * a * b3 * c * d5 * e * f2 * g3 +
                            16380 * a2 * b3 * d6 * e * f2 * g3 + 17640 * b4 * c * d4 * e * f3 * g3 +
                            17640 * a * b4 * d5 * e * f3 * g3 + 6300 * b5 * d4 * e * f4 * g3 +
                            1008 * b3 * c2 * d5 * f * g4 + 1680 * a * b3 * c * d6 * f * g4 +
                            720 * a2 * b3 * d7 * f * g4 + 2142 * b4 * c * d5 * f2 * g4 +
                            1785 * a * b4 * d6 * f2 * g4 + 1134 * b5 * d5 * f3 * g4 -
                            126 * c5 * d * e4 * n2 - 1575 * a * c4 * d2 * e4 * n2 -
                            6300 * a2 * c3 * d3 * e4 * n2 - 11025 * a3 * c2 * d4 * e4 * n2 -
                            8820 * a4 * c * d5 * e4 * n2 - 2646 * a5 * d6 * e4 * n2 -
                            630 * b * c4 * d * e4 * f * n2 - 6300 * a * b * c3 * d2 * e4 * f * n2 -
                            18900 * a2 * b * c2 * d3 * e4 * f * n2 -
                            22050 * a3 * b * c * d4 * e4 * f * n2 -
                            8820 * a4 * b * d5 * e4 * f * n2 - 315 * a * c4 * e4 * f2 * n2 -
                            420 * c3 * d * e4 * f2 * n2 - 2520 * a2 * c3 * d * e4 * f2 * n2 -
                            1260 * b2 * c3 * d * e4 * f2 * n2 - 3150 * a * c2 * d2 * e4 * f2 * n2 -
                            6300 * a3 * c2 * d2 * e4 * f2 * n2 -
                            9450 * a * b2 * c2 * d2 * e4 * f2 * n2 -
                            6300 * a2 * c * d3 * e4 * f2 * n2 - 6300 * a4 * c * d3 * e4 * f2 * n2 -
                            18900 * a2 * b2 * c * d3 * e4 * f2 * n2 -
                            3675 * a3 * d4 * e4 * f2 * n2 - 2205 * a5 * d4 * e4 * f2 * n2 -
                            11025 * a3 * b2 * d4 * e4 * f2 * n2 - 1260 * a * b * c3 * e4 * f3 * n2 -
                            1260 * b * c2 * d * e4 * f3 * n2 -
                            7560 * a2 * b * c2 * d * e4 * f3 * n2 -
                            1260 * b3 * c2 * d * e4 * f3 * n2 -
                            6300 * a * b * c * d2 * e4 * f3 * n2 -
                            12600 * a3 * b * c * d2 * e4 * f3 * n2 -
                            6300 * a * b3 * c * d2 * e4 * f3 * n2 -
                            6300 * a2 * b * d3 * e4 * f3 * n2 - 6300 * a4 * b * d3 * e4 * f3 * n2 -
                            6300 * a2 * b3 * d3 * e4 * f3 * n2 - 1890 * a * b2 * c2 * e4 * f4 * n2 -
                            1260 * b2 * c * d * e4 * f4 * n2 -
                            7560 * a2 * b2 * c * d * e4 * f4 * n2 -
                            630 * b4 * c * d * e4 * f4 * n2 - 3150 * a * b2 * d2 * e4 * f4 * n2 -
                            6300 * a3 * b2 * d2 * e4 * f4 * n2 - 1575 * a * b4 * d2 * e4 * f4 * n2 -
                            1260 * a * b3 * c * e4 * f5 * n2 - 420 * b3 * d * e4 * f5 * n2 -
                            2520 * a2 * b3 * d * e4 * f5 * n2 - 126 * b5 * d * e4 * f5 * n2 -
                            315 * a * b4 * e4 * f6 * n2 - 1260 * b * c4 * d2 * e3 * g * n2 -
                            8400 * a * b * c3 * d3 * e3 * g * n2 -
                            18900 * a2 * b * c2 * d4 * e3 * g * n2 -
                            17640 * a3 * b * c * d5 * e3 * g * n2 -
                            5880 * a4 * b * d6 * e3 * g * n2 - 336 * c5 * e3 * f * g * n2 -
                            5040 * a * c4 * d * e3 * f * g * n2 - 3360 * c3 * d2 * e3 * f * g * n2 -
                            20160 * a2 * c3 * d2 * e3 * f * g * n2 -
                            5040 * b2 * c3 * d2 * e3 * f * g * n2 -
                            16800 * a * c2 * d3 * e3 * f * g * n2 -
                            33600 * a3 * c2 * d3 * e3 * f * g * n2 -
                            25200 * a * b2 * c2 * d3 * e3 * f * g * n2 -
                            25200 * a2 * c * d4 * e3 * f * g * n2 -
                            25200 * a4 * c * d4 * e3 * f * g * n2 -
                            37800 * a2 * b2 * c * d4 * e3 * f * g * n2 -
                            11760 * a3 * d5 * e3 * f * g * n2 - 7056 * a5 * d5 * e3 * f * g * n2 -
                            17640 * a3 * b2 * d5 * e3 * f * g * n2 -
                            2100 * b * c4 * e3 * f2 * g * n2 -
                            25200 * a * b * c3 * d * e3 * f2 * g * n2 -
                            12600 * b * c2 * d2 * e3 * f2 * g * n2 -
                            75600 * a2 * b * c2 * d2 * e3 * f2 * g * n2 -
                            7560 * b3 * c2 * d2 * e3 * f2 * g * n2 -
                            42000 * a * b * c * d3 * e3 * f2 * g * n2 -
                            84000 * a3 * b * c * d3 * e3 * f2 * g * n2 -
                            25200 * a * b3 * c * d3 * e3 * f2 * g * n2 -
                            31500 * a2 * b * d4 * e3 * f2 * g * n2 -
                            31500 * a4 * b * d4 * e3 * f2 * g * n2 -
                            18900 * a2 * b3 * d4 * e3 * f2 * g * n2 -
                            5040 * b2 * c3 * e3 * f3 * g * n2 -
                            45360 * a * b2 * c2 * d * e3 * f3 * g * n2 -
                            15120 * b2 * c * d2 * e3 * f3 * g * n2 -
                            90720 * a2 * b2 * c * d2 * e3 * f3 * g * n2 -
                            5040 * b4 * c * d2 * e3 * f3 * g * n2 -
                            25200 * a * b2 * d3 * e3 * f3 * g * n2 -
                            50400 * a3 * b2 * d3 * e3 * f3 * g * n2 -
                            8400 * a * b4 * d3 * e3 * f3 * g * n2 -
                            5880 * b3 * c2 * e3 * f4 * g * n2 -
                            35280 * a * b3 * c * d * e3 * f4 * g * n2 -
                            5880 * b3 * d2 * e3 * f4 * g * n2 -
                            35280 * a2 * b3 * d2 * e3 * f4 * g * n2 -
                            1260 * b5 * d2 * e3 * f4 * g * n2 - 3360 * b4 * c * e3 * f5 * g * n2 -
                            10080 * a * b4 * d * e3 * f5 * g * n2 - 756 * b5 * e3 * f6 * g * n2 -
                            2520 * b2 * c3 * d3 * e2 * g2 * n2 -
                            9450 * a * b2 * c2 * d4 * e2 * g2 * n2 -
                            11340 * a2 * b2 * c * d5 * e2 * g2 * n2 -
                            4410 * a3 * b2 * d6 * e2 * g2 * n2 -
                            5040 * b * c4 * d * e2 * f * g2 * n2 -
                            30240 * a * b * c3 * d2 * e2 * f * g2 * n2 -
                            10080 * b * c2 * d3 * e2 * f * g2 * n2 -
                            60480 * a2 * b * c2 * d3 * e2 * f * g2 * n2 -
                            7560 * b3 * c2 * d3 * e2 * f * g2 * n2 -
                            25200 * a * b * c * d4 * e2 * f * g2 * n2 -
                            50400 * a3 * b * c * d4 * e2 * f * g2 * n2 -
                            18900 * a * b3 * c * d4 * e2 * f * g2 * n2 -
                            15120 * a2 * b * d5 * e2 * f * g2 * n2 -
                            15120 * a4 * b * d5 * e2 * f * g2 * n2 -
                            11340 * a2 * b3 * d5 * e2 * f * g2 * n2 -
                            22680 * b2 * c3 * d * e2 * f2 * g2 * n2 -
                            102060 * a * b2 * c2 * d2 * e2 * f2 * g2 * n2 -
                            22680 * b2 * c * d3 * e2 * f2 * g2 * n2 -
                            136080 * a2 * b2 * c * d3 * e2 * f2 * g2 * n2 -
                            7560 * b4 * c * d3 * e2 * f2 * g2 * n2 -
                            28350 * a * b2 * d4 * e2 * f2 * g2 * n2 -
                            56700 * a3 * b2 * d4 * e2 * f2 * g2 * n2 -
                            9450 * a * b4 * d4 * e2 * f2 * g2 * n2 -
                            37800 * b3 * c2 * d * e2 * f3 * g2 * n2 -
                            113400 * a * b3 * c * d2 * e2 * f3 * g2 * n2 -
                            12600 * b3 * d3 * e2 * f3 * g2 * n2 -
                            75600 * a2 * b3 * d3 * e2 * f3 * g2 * n2 -
                            2520 * b5 * d3 * e2 * f3 * g2 * n2 -
                            27720 * b4 * c * d * e2 * f4 * g2 * n2 -
                            41580 * a * b4 * d2 * e2 * f4 * g2 * n2 -
                            7560 * b5 * d * e2 * f5 * g2 * n2 - 1260 * b3 * c2 * d4 * e * g3 * n2 -
                            2520 * a * b3 * c * d5 * e * g3 * n2 -
                            1260 * a2 * b3 * d6 * e * g3 * n2 -
                            10080 * b2 * c3 * d2 * e * f * g3 * n2 -
                            30240 * a * b2 * c2 * d3 * e * f * g3 * n2 -
                            5040 * b2 * c * d4 * e * f * g3 * n2 -
                            30240 * a2 * b2 * c * d4 * e * f * g3 * n2 -
                            2520 * b4 * c * d4 * e * f * g3 * n2 -
                            5040 * a * b2 * d5 * e * f * g3 * n2 -
                            10080 * a3 * b2 * d5 * e * f * g3 * n2 -
                            2520 * a * b4 * d5 * e * f * g3 * n2 -
                            32760 * b3 * c2 * d2 * e * f2 * g3 * n2 -
                            65520 * a * b3 * c * d3 * e * f2 * g3 * n2 -
                            5460 * b3 * d4 * e * f2 * g3 * n2 -
                            32760 * a2 * b3 * d4 * e * f2 * g3 * n2 -
                            1260 * b5 * d4 * e * f2 * g3 * n2 -
                            35280 * b4 * c * d2 * e * f3 * g3 * n2 -
                            35280 * a * b4 * d3 * e * f3 * g3 * n2 -
                            12600 * b5 * d2 * e * f4 * g3 * n2 - 126 * b4 * c * d5 * g4 * n2 -
                            105 * a * b4 * d6 * g4 * n2 - 3360 * b3 * c2 * d3 * f * g4 * n2 -
                            5040 * a * b3 * c * d4 * f * g4 * n2 - 336 * b3 * d5 * f * g4 * n2 -
                            2016 * a2 * b3 * d5 * f * g4 * n2 - 126 * b5 * d5 * f * g4 * n2 -
                            7140 * b4 * c * d3 * f2 * g4 * n2 - 5355 * a * b4 * d4 * f2 * g4 * n2 -
                            3780 * b5 * d3 * f3 * g4 * n2 + 315 * a * c4 * e4 * n4 +
                            420 * c3 * d * e4 * n4 + 2520 * a2 * c3 * d * e4 * n4 +
                            3150 * a * c2 * d2 * e4 * n4 + 6300 * a3 * c2 * d2 * e4 * n4 +
                            6300 * a2 * c * d3 * e4 * n4 + 6300 * a4 * c * d3 * e4 * n4 +
                            3675 * a3 * d4 * e4 * n4 + 2205 * a5 * d4 * e4 * n4 +
                            1260 * a * b * c3 * e4 * f * n4 + 1260 * b * c2 * d * e4 * f * n4 +
                            7560 * a2 * b * c2 * d * e4 * f * n4 +
                            6300 * a * b * c * d2 * e4 * f * n4 +
                            12600 * a3 * b * c * d2 * e4 * f * n4 +
                            6300 * a2 * b * d3 * e4 * f * n4 + 6300 * a4 * b * d3 * e4 * f * n4 +
                            630 * a * c2 * e4 * f2 * n4 + 315 * a3 * c2 * e4 * f2 * n4 +
                            1890 * a * b2 * c2 * e4 * f2 * n4 + 630 * c * d * e4 * f2 * n4 +
                            2520 * a2 * c * d * e4 * f2 * n4 + 630 * a4 * c * d * e4 * f2 * n4 +
                            1260 * b2 * c * d * e4 * f2 * n4 +
                            7560 * a2 * b2 * c * d * e4 * f2 * n4 + 1575 * a * d2 * e4 * f2 * n4 +
                            2100 * a3 * d2 * e4 * f2 * n4 + 315 * a5 * d2 * e4 * f2 * n4 +
                            3150 * a * b2 * d2 * e4 * f2 * n4 + 6300 * a3 * b2 * d2 * e4 * f2 * n4 +
                            1260 * a * b * c * e4 * f3 * n4 + 630 * a3 * b * c * e4 * f3 * n4 +
                            1260 * a * b3 * c * e4 * f3 * n4 + 630 * b * d * e4 * f3 * n4 +
                            2520 * a2 * b * d * e4 * f3 * n4 + 630 * a4 * b * d * e4 * f3 * n4 +
                            420 * b3 * d * e4 * f3 * n4 + 2520 * a2 * b3 * d * e4 * f3 * n4 +
                            630 * a * b2 * e4 * f4 * n4 + 315 * a3 * b2 * e4 * f4 * n4 +
                            315 * a * b4 * e4 * f4 * n4 + 420 * b * c4 * e3 * g * n4 +
                            5040 * a * b * c3 * d * e3 * g * n4 + 2520 * b * c2 * d2 * e3 * g * n4 +
                            15120 * a2 * b * c2 * d2 * e3 * g * n4 +
                            8400 * a * b * c * d3 * e3 * g * n4 +
                            16800 * a3 * b * c * d3 * e3 * g * n4 +
                            6300 * a2 * b * d4 * e3 * g * n4 + 6300 * a4 * b * d4 * e3 * g * n4 +
                            1120 * c3 * e3 * f * g * n4 + 1680 * a2 * c3 * e3 * f * g * n4 +
                            1680 * b2 * c3 * e3 * f * g * n4 +
                            10080 * a * c2 * d * e3 * f * g * n4 +
                            5040 * a3 * c2 * d * e3 * f * g * n4 +
                            15120 * a * b2 * c2 * d * e3 * f * g * n4 +
                            5040 * c * d2 * e3 * f * g * n4 +
                            20160 * a2 * c * d2 * e3 * f * g * n4 +
                            5040 * a4 * c * d2 * e3 * f * g * n4 +
                            5040 * b2 * c * d2 * e3 * f * g * n4 +
                            30240 * a2 * b2 * c * d2 * e3 * f * g * n4 +
                            8400 * a * d3 * e3 * f * g * n4 + 11200 * a3 * d3 * e3 * f * g * n4 +
                            1680 * a5 * d3 * e3 * f * g * n4 +
                            8400 * a * b2 * d3 * e3 * f * g * n4 +
                            16800 * a3 * b2 * d3 * e3 * f * g * n4 +
                            4200 * b * c2 * e3 * f2 * g * n4 +
                            6300 * a2 * b * c2 * e3 * f2 * g * n4 +
                            2520 * b3 * c2 * e3 * f2 * g * n4 +
                            25200 * a * b * c * d * e3 * f2 * g * n4 +
                            12600 * a3 * b * c * d * e3 * f2 * g * n4 +
                            15120 * a * b3 * c * d * e3 * f2 * g * n4 +
                            6300 * b * d2 * e3 * f2 * g * n4 +
                            25200 * a2 * b * d2 * e3 * f2 * g * n4 +
                            6300 * a4 * b * d2 * e3 * f2 * g * n4 +
                            2520 * b3 * d2 * e3 * f2 * g * n4 +
                            15120 * a2 * b3 * d2 * e3 * f2 * g * n4 +
                            5040 * b2 * c * e3 * f3 * g * n4 +
                            7560 * a2 * b2 * c * e3 * f3 * g * n4 +
                            1680 * b4 * c * e3 * f3 * g * n4 +
                            15120 * a * b2 * d * e3 * f3 * g * n4 +
                            7560 * a3 * b2 * d * e3 * f3 * g * n4 +
                            5040 * a * b4 * d * e3 * f3 * g * n4 + 1960 * b3 * e3 * f4 * g * n4 +
                            2940 * a2 * b3 * e3 * f4 * g * n4 + 420 * b5 * e3 * f4 * g * n4 +
                            2520 * b2 * c3 * d * e2 * g2 * n4 +
                            11340 * a * b2 * c2 * d2 * e2 * g2 * n4 +
                            2520 * b2 * c * d3 * e2 * g2 * n4 +
                            15120 * a2 * b2 * c * d3 * e2 * g2 * n4 +
                            3150 * a * b2 * d4 * e2 * g2 * n4 + 6300 * a3 * b2 * d4 * e2 * g2 * n4 +
                            5040 * a * b * c3 * e2 * f * g2 * n4 +
                            10080 * b * c2 * d * e2 * f * g2 * n4 +
                            15120 * a2 * b * c2 * d * e2 * f * g2 * n4 +
                            7560 * b3 * c2 * d * e2 * f * g2 * n4 +
                            30240 * a * b * c * d2 * e2 * f * g2 * n4 +
                            15120 * a3 * b * c * d2 * e2 * f * g2 * n4 +
                            22680 * a * b3 * c * d2 * e2 * f * g2 * n4 +
                            5040 * b * d3 * e2 * f * g2 * n4 +
                            20160 * a2 * b * d3 * e2 * f * g2 * n4 +
                            5040 * a4 * b * d3 * e2 * f * g2 * n4 +
                            2520 * b3 * d3 * e2 * f * g2 * n4 +
                            15120 * a2 * b3 * d3 * e2 * f * g2 * n4 +
                            17010 * a * b2 * c2 * e2 * f2 * g2 * n4 +
                            22680 * b2 * c * d * e2 * f2 * g2 * n4 +
                            34020 * a2 * b2 * c * d * e2 * f2 * g2 * n4 +
                            7560 * b4 * c * d * e2 * f2 * g2 * n4 +
                            34020 * a * b2 * d2 * e2 * f2 * g2 * n4 +
                            17010 * a3 * b2 * d2 * e2 * f2 * g2 * n4 +
                            11340 * a * b4 * d2 * e2 * f2 * g2 * n4 +
                            18900 * a * b3 * c * e2 * f3 * g2 * n4 +
                            12600 * b3 * d * e2 * f3 * g2 * n4 +
                            18900 * a2 * b3 * d * e2 * f3 * g2 * n4 +
                            2520 * b5 * d * e2 * f3 * g2 * n4 + 6930 * a * b4 * e2 * f4 * g2 * n4 +
                            2520 * b3 * c2 * d2 * e * g3 * n4 +
                            5040 * a * b3 * c * d3 * e * g3 * n4 + 420 * b3 * d4 * e * g3 * n4 +
                            2520 * a2 * b3 * d4 * e * g3 * n4 + 5040 * b2 * c3 * e * f * g3 * n4 +
                            15120 * a * b2 * c2 * d * e * f * g3 * n4 +
                            10080 * b2 * c * d2 * e * f * g3 * n4 +
                            15120 * a2 * b2 * c * d2 * e * f * g3 * n4 +
                            5040 * b4 * c * d2 * e * f * g3 * n4 +
                            10080 * a * b2 * d3 * e * f * g3 * n4 +
                            5040 * a3 * b2 * d3 * e * f * g3 * n4 +
                            5040 * a * b4 * d3 * e * f * g3 * n4 +
                            16380 * b3 * c2 * e * f2 * g3 * n4 +
                            32760 * a * b3 * c * d * e * f2 * g3 * n4 +
                            10920 * b3 * d2 * e * f2 * g3 * n4 +
                            16380 * a2 * b3 * d2 * e * f2 * g3 * n4 +
                            2520 * b5 * d2 * e * f2 * g3 * n4 + 17640 * b4 * c * e * f3 * g3 * n4 +
                            17640 * a * b4 * d * e * f3 * g3 * n4 + 6300 * b5 * e * f4 * g3 * n4 +
                            420 * b4 * c * d3 * g4 * n4 + 315 * a * b4 * d4 * g4 * n4 +
                            5040 * b3 * c2 * d * f * g4 * n4 +
                            5040 * a * b3 * c * d2 * f * g4 * n4 + 1120 * b3 * d3 * f * g4 * n4 +
                            1680 * a2 * b3 * d3 * f * g4 * n4 + 420 * b5 * d3 * f * g4 * n4 +
                            10710 * b4 * c * d * f2 * g4 * n4 + 5355 * a * b4 * d2 * f2 * g4 * n4 +
                            5670 * b5 * d * f3 * g4 * n4 - 630 * a * c2 * e4 * n6 -
                            315 * a3 * c2 * e4 * n6 - 630 * c * d * e4 * n6 -
                            2520 * a2 * c * d * e4 * n6 - 630 * a4 * c * d * e4 * n6 -
                            1575 * a * d2 * e4 * n6 - 2100 * a3 * d2 * e4 * n6 -
                            315 * a5 * d2 * e4 * n6 - 1260 * a * b * c * e4 * f * n6 -
                            630 * a3 * b * c * e4 * f * n6 - 630 * b * d * e4 * f * n6 -
                            2520 * a2 * b * d * e4 * f * n6 - 630 * a4 * b * d * e4 * f * n6 -
                            315 * a * e4 * f2 * n6 - 105 * a3 * e4 * f2 * n6 -
                            630 * a * b2 * e4 * f2 * n6 - 315 * a3 * b2 * e4 * f2 * n6 -
                            840 * b * c2 * e3 * g * n6 - 1260 * a2 * b * c2 * e3 * g * n6 -
                            5040 * a * b * c * d * e3 * g * n6 -
                            2520 * a3 * b * c * d * e3 * g * n6 - 1260 * b * d2 * e3 * g * n6 -
                            5040 * a2 * b * d2 * e3 * g * n6 - 1260 * a4 * b * d2 * e3 * g * n6 -
                            1680 * c * e3 * f * g * n6 - 1680 * a2 * c * e3 * f * g * n6 -
                            1680 * b2 * c * e3 * f * g * n6 - 2520 * a2 * b2 * c * e3 * f * g * n6 -
                            5040 * a * d * e3 * f * g * n6 - 1680 * a3 * d * e3 * f * g * n6 -
                            5040 * a * b2 * d * e3 * f * g * n6 -
                            2520 * a3 * b2 * d * e3 * f * g * n6 - 2100 * b * e3 * f2 * g * n6 -
                            2100 * a2 * b * e3 * f2 * g * n6 - 840 * b3 * e3 * f2 * g * n6 -
                            1260 * a2 * b3 * e3 * f2 * g * n6 - 1890 * a * b2 * c2 * e2 * g2 * n6 -
                            2520 * b2 * c * d * e2 * g2 * n6 -
                            3780 * a2 * b2 * c * d * e2 * g2 * n6 -
                            3780 * a * b2 * d2 * e2 * g2 * n6 - 1890 * a3 * b2 * d2 * e2 * g2 * n6 -
                            5040 * a * b * c * e2 * f * g2 * n6 -
                            3780 * a * b3 * c * e2 * f * g2 * n6 - 5040 * b * d * e2 * f * g2 * n6 -
                            5040 * a2 * b * d * e2 * f * g2 * n6 -
                            2520 * b3 * d * e2 * f * g2 * n6 -
                            3780 * a2 * b3 * d * e2 * f * g2 * n6 -
                            5670 * a * b2 * e2 * f2 * g2 * n6 - 1890 * a * b4 * e2 * f2 * g2 * n6 -
                            1260 * b3 * c2 * e * g3 * n6 - 2520 * a * b3 * c * d * e * g3 * n6 -
                            840 * b3 * d2 * e * g3 * n6 - 1260 * a2 * b3 * d2 * e * g3 * n6 -
                            5040 * b2 * c * e * f * g3 * n6 - 2520 * b4 * c * e * f * g3 * n6 -
                            5040 * a * b2 * d * e * f * g3 * n6 -
                            2520 * a * b4 * d * e * f * g3 * n6 - 5460 * b3 * e * f2 * g3 * n6 -
                            1260 * b5 * e * f2 * g3 * n6 - 630 * b4 * c * d * g4 * n6 -
                            315 * a * b4 * d2 * g4 * n6 - 1680 * b3 * d * f * g4 * n6 -
                            630 * b5 * d * f * g4 * n6 + 315 * a * e4 * n8 + 105 * a3 * e4 * n8 +
                            420 * b * e3 * g * n8 + 420 * a2 * b * e3 * g * n8 +
                            630 * a * b2 * e2 * g2 * n8 + 420 * b3 * e * g3 * n8) +
                    14 *
                        (3 * c5 * e5 * f2 + 75 * a * c4 * d * e5 * f2 +
                            450 * a2 * c3 * d2 * e5 * f2 + 1050 * a3 * c2 * d3 * e5 * f2 +
                            1050 * a4 * c * d4 * e5 * f2 + 378 * a5 * d5 * e5 * f2 +
                            15 * b * c4 * e5 * f3 + 300 * a * b * c3 * d * e5 * f3 +
                            1350 * a2 * b * c2 * d2 * e5 * f3 + 2100 * a3 * b * c * d3 * e5 * f3 +
                            1050 * a4 * b * d4 * e5 * f3 + 30 * b2 * c3 * e5 * f4 +
                            450 * a * b2 * c2 * d * e5 * f4 + 1350 * a2 * b2 * c * d2 * e5 * f4 +
                            1050 * a3 * b2 * d3 * e5 * f4 + 30 * b3 * c2 * e5 * f5 +
                            300 * a * b3 * c * d * e5 * f5 + 450 * a2 * b3 * d2 * e5 * f5 +
                            15 * b4 * c * e5 * f6 + 75 * a * b4 * d * e5 * f6 + 3 * b5 * e5 * f7 +
                            60 * c5 * d * e4 * f * g + 750 * a * c4 * d2 * e4 * f * g +
                            3000 * a2 * c3 * d3 * e4 * f * g + 5250 * a3 * c2 * d4 * e4 * f * g +
                            4200 * a4 * c * d5 * e4 * f * g + 1260 * a5 * d6 * e4 * f * g +
                            375 * b * c4 * d * e4 * f2 * g + 3750 * a * b * c3 * d2 * e4 * f2 * g +
                            11250 * a2 * b * c2 * d3 * e4 * f2 * g +
                            13125 * a3 * b * c * d4 * e4 * f2 * g +
                            5250 * a4 * b * d5 * e4 * f2 * g + 900 * b2 * c3 * d * e4 * f3 * g +
                            6750 * a * b2 * c2 * d2 * e4 * f3 * g +
                            13500 * a2 * b2 * c * d3 * e4 * f3 * g +
                            7875 * a3 * b2 * d4 * e4 * f3 * g + 1050 * b3 * c2 * d * e4 * f4 * g +
                            5250 * a * b3 * c * d2 * e4 * f4 * g +
                            5250 * a2 * b3 * d3 * e4 * f4 * g + 600 * b4 * c * d * e4 * f5 * g +
                            1500 * a * b4 * d2 * e4 * f5 * g + 135 * b5 * d * e4 * f6 * g +
                            600 * b * c4 * d2 * e3 * f * g2 + 4000 * a * b * c3 * d3 * e3 * f * g2 +
                            9000 * a2 * b * c2 * d4 * e3 * f * g2 +
                            8400 * a3 * b * c * d5 * e3 * f * g2 +
                            2800 * a4 * b * d6 * e3 * f * g2 + 2700 * b2 * c3 * d2 * e3 * f2 * g2 +
                            13500 * a * b2 * c2 * d3 * e3 * f2 * g2 +
                            20250 * a2 * b2 * c * d4 * e3 * f2 * g2 +
                            9450 * a3 * b2 * d5 * e3 * f2 * g2 +
                            4500 * b3 * c2 * d2 * e3 * f3 * g2 +
                            15000 * a * b3 * c * d3 * e3 * f3 * g2 +
                            11250 * a2 * b3 * d4 * e3 * f3 * g2 +
                            3300 * b4 * c * d2 * e3 * f4 * g2 + 5500 * a * b4 * d3 * e3 * f4 * g2 +
                            900 * b5 * d2 * e3 * f5 * g2 + 1200 * b2 * c3 * d3 * e2 * f * g3 +
                            4500 * a * b2 * c2 * d4 * e2 * f * g3 +
                            5400 * a2 * b2 * c * d5 * e2 * f * g3 +
                            2100 * a3 * b2 * d6 * e2 * f * g3 + 3900 * b3 * c2 * d3 * e2 * f2 * g3 +
                            9750 * a * b3 * c * d4 * e2 * f2 * g3 +
                            5850 * a2 * b3 * d5 * e2 * f2 * g3 + 4200 * b4 * c * d3 * e2 * f3 * g3 +
                            5250 * a * b4 * d4 * e2 * f3 * g3 + 1500 * b5 * d3 * e2 * f4 * g3 +
                            600 * b3 * c2 * d4 * e * f * g4 + 1200 * a * b3 * c * d5 * e * f * g4 +
                            600 * a2 * b3 * d6 * e * f * g4 + 1275 * b4 * c * d4 * e * f2 * g4 +
                            1275 * a * b4 * d5 * e * f2 * g4 + 675 * b5 * d4 * e * f3 * g4 +
                            60 * b4 * c * d5 * f * g5 + 50 * a * b4 * d6 * f * g5 +
                            63 * b5 * d5 * f2 * g5 - 3 * c5 * e5 * n2 - 75 * a * c4 * d * e5 * n2 -
                            450 * a2 * c3 * d2 * e5 * n2 - 1050 * a3 * c2 * d3 * e5 * n2 -
                            1050 * a4 * c * d4 * e5 * n2 - 378 * a5 * d5 * e5 * n2 -
                            15 * b * c4 * e5 * f * n2 - 300 * a * b * c3 * d * e5 * f * n2 -
                            1350 * a2 * b * c2 * d2 * e5 * f * n2 -
                            2100 * a3 * b * c * d3 * e5 * f * n2 -
                            1050 * a4 * b * d4 * e5 * f * n2 - 10 * c3 * e5 * f2 * n2 -
                            60 * a2 * c3 * e5 * f2 * n2 - 30 * b2 * c3 * e5 * f2 * n2 -
                            150 * a * c2 * d * e5 * f2 * n2 - 300 * a3 * c2 * d * e5 * f2 * n2 -
                            450 * a * b2 * c2 * d * e5 * f2 * n2 -
                            450 * a2 * c * d2 * e5 * f2 * n2 - 450 * a4 * c * d2 * e5 * f2 * n2 -
                            1350 * a2 * b2 * c * d2 * e5 * f2 * n2 - 350 * a3 * d3 * e5 * f2 * n2 -
                            210 * a5 * d3 * e5 * f2 * n2 - 1050 * a3 * b2 * d3 * e5 * f2 * n2 -
                            30 * b * c2 * e5 * f3 * n2 - 180 * a2 * b * c2 * e5 * f3 * n2 -
                            30 * b3 * c2 * e5 * f3 * n2 - 300 * a * b * c * d * e5 * f3 * n2 -
                            600 * a3 * b * c * d * e5 * f3 * n2 -
                            300 * a * b3 * c * d * e5 * f3 * n2 - 450 * a2 * b * d2 * e5 * f3 * n2 -
                            450 * a4 * b * d2 * e5 * f3 * n2 - 450 * a2 * b3 * d2 * e5 * f3 * n2 -
                            30 * b2 * c * e5 * f4 * n2 - 180 * a2 * b2 * c * e5 * f4 * n2 -
                            15 * b4 * c * e5 * f4 * n2 - 150 * a * b2 * d * e5 * f4 * n2 -
                            300 * a3 * b2 * d * e5 * f4 * n2 - 75 * a * b4 * d * e5 * f4 * n2 -
                            10 * b3 * e5 * f5 * n2 - 60 * a2 * b3 * e5 * f5 * n2 -
                            3 * b5 * e5 * f5 * n2 - 75 * b * c4 * d * e4 * g * n2 -
                            750 * a * b * c3 * d2 * e4 * g * n2 -
                            2250 * a2 * b * c2 * d3 * e4 * g * n2 -
                            2625 * a3 * b * c * d4 * e4 * g * n2 -
                            1050 * a4 * b * d5 * e4 * g * n2 - 150 * a * c4 * e4 * f * g * n2 -
                            200 * c3 * d * e4 * f * g * n2 - 1200 * a2 * c3 * d * e4 * f * g * n2 -
                            300 * b2 * c3 * d * e4 * f * g * n2 -
                            1500 * a * c2 * d2 * e4 * f * g * n2 -
                            3000 * a3 * c2 * d2 * e4 * f * g * n2 -
                            2250 * a * b2 * c2 * d2 * e4 * f * g * n2 -
                            3000 * a2 * c * d3 * e4 * f * g * n2 -
                            3000 * a4 * c * d3 * e4 * f * g * n2 -
                            4500 * a2 * b2 * c * d3 * e4 * f * g * n2 -
                            1750 * a3 * d4 * e4 * f * g * n2 - 1050 * a5 * d4 * e4 * f * g * n2 -
                            2625 * a3 * b2 * d4 * e4 * f * g * n2 -
                            750 * a * b * c3 * e4 * f2 * g * n2 -
                            750 * b * c2 * d * e4 * f2 * g * n2 -
                            4500 * a2 * b * c2 * d * e4 * f2 * g * n2 -
                            450 * b3 * c2 * d * e4 * f2 * g * n2 -
                            3750 * a * b * c * d2 * e4 * f2 * g * n2 -
                            7500 * a3 * b * c * d2 * e4 * f2 * g * n2 -
                            2250 * a * b3 * c * d2 * e4 * f2 * g * n2 -
                            3750 * a2 * b * d3 * e4 * f2 * g * n2 -
                            3750 * a4 * b * d3 * e4 * f2 * g * n2 -
                            2250 * a2 * b3 * d3 * e4 * f2 * g * n2 -
                            1350 * a * b2 * c2 * e4 * f3 * g * n2 -
                            900 * b2 * c * d * e4 * f3 * g * n2 -
                            5400 * a2 * b2 * c * d * e4 * f3 * g * n2 -
                            300 * b4 * c * d * e4 * f3 * g * n2 -
                            2250 * a * b2 * d2 * e4 * f3 * g * n2 -
                            4500 * a3 * b2 * d2 * e4 * f3 * g * n2 -
                            750 * a * b4 * d2 * e4 * f3 * g * n2 -
                            1050 * a * b3 * c * e4 * f4 * g * n2 - 350 * b3 * d * e4 * f4 * g * n2 -
                            2100 * a2 * b3 * d * e4 * f4 * g * n2 - 75 * b5 * d * e4 * f4 * g * n2 -
                            300 * a * b4 * e4 * f5 * g * n2 - 300 * b2 * c3 * d2 * e3 * g2 * n2 -
                            1500 * a * b2 * c2 * d3 * e3 * g2 * n2 -
                            2250 * a2 * b2 * c * d4 * e3 * g2 * n2 -
                            1050 * a3 * b2 * d5 * e3 * g2 * n2 - 200 * b * c4 * e3 * f * g2 * n2 -
                            2400 * a * b * c3 * d * e3 * f * g2 * n2 -
                            1200 * b * c2 * d2 * e3 * f * g2 * n2 -
                            7200 * a2 * b * c2 * d2 * e3 * f * g2 * n2 -
                            900 * b3 * c2 * d2 * e3 * f * g2 * n2 -
                            4000 * a * b * c * d3 * e3 * f * g2 * n2 -
                            8000 * a3 * b * c * d3 * e3 * f * g2 * n2 -
                            3000 * a * b3 * c * d3 * e3 * f * g2 * n2 -
                            3000 * a2 * b * d4 * e3 * f * g2 * n2 -
                            3000 * a4 * b * d4 * e3 * f * g2 * n2 -
                            2250 * a2 * b3 * d4 * e3 * f * g2 * n2 -
                            900 * b2 * c3 * e3 * f2 * g2 * n2 -
                            8100 * a * b2 * c2 * d * e3 * f2 * g2 * n2 -
                            2700 * b2 * c * d2 * e3 * f2 * g2 * n2 -
                            16200 * a2 * b2 * c * d2 * e3 * f2 * g2 * n2 -
                            900 * b4 * c * d2 * e3 * f2 * g2 * n2 -
                            4500 * a * b2 * d3 * e3 * f2 * g2 * n2 -
                            9000 * a3 * b2 * d3 * e3 * f2 * g2 * n2 -
                            1500 * a * b4 * d3 * e3 * f2 * g2 * n2 -
                            1500 * b3 * c2 * e3 * f3 * g2 * n2 -
                            9000 * a * b3 * c * d * e3 * f3 * g2 * n2 -
                            1500 * b3 * d2 * e3 * f3 * g2 * n2 -
                            9000 * a2 * b3 * d2 * e3 * f3 * g2 * n2 -
                            300 * b5 * d2 * e3 * f3 * g2 * n2 - 1100 * b4 * c * e3 * f4 * g2 * n2 -
                            3300 * a * b4 * d * e3 * f4 * g2 * n2 - 300 * b5 * e3 * f5 * g2 * n2 -
                            300 * b3 * c2 * d3 * e2 * g3 * n2 -
                            750 * a * b3 * c * d4 * e2 * g3 * n2 -
                            450 * a2 * b3 * d5 * e2 * g3 * n2 -
                            1200 * b2 * c3 * d * e2 * f * g3 * n2 -
                            5400 * a * b2 * c2 * d2 * e2 * f * g3 * n2 -
                            1200 * b2 * c * d3 * e2 * f * g3 * n2 -
                            7200 * a2 * b2 * c * d3 * e2 * f * g3 * n2 -
                            600 * b4 * c * d3 * e2 * f * g3 * n2 -
                            1500 * a * b2 * d4 * e2 * f * g3 * n2 -
                            3000 * a3 * b2 * d4 * e2 * f * g3 * n2 -
                            750 * a * b4 * d4 * e2 * f * g3 * n2 -
                            3900 * b3 * c2 * d * e2 * f2 * g3 * n2 -
                            11700 * a * b3 * c * d2 * e2 * f2 * g3 * n2 -
                            1300 * b3 * d3 * e2 * f2 * g3 * n2 -
                            7800 * a2 * b3 * d3 * e2 * f2 * g3 * n2 -
                            300 * b5 * d3 * e2 * f2 * g3 * n2 -
                            4200 * b4 * c * d * e2 * f3 * g3 * n2 -
                            6300 * a * b4 * d2 * e2 * f3 * g3 * n2 -
                            1500 * b5 * d * e2 * f4 * g3 * n2 - 75 * b4 * c * d4 * e * g4 * n2 -
                            75 * a * b4 * d5 * e * g4 * n2 - 1200 * b3 * c2 * d2 * e * f * g4 * n2 -
                            2400 * a * b3 * c * d3 * e * f * g4 * n2 -
                            200 * b3 * d4 * e * f * g4 * n2 -
                            1200 * a2 * b3 * d4 * e * f * g4 * n2 - 75 * b5 * d4 * e * f * g4 * n2 -
                            2550 * b4 * c * d2 * e * f2 * g4 * n2 -
                            2550 * a * b4 * d3 * e * f2 * g4 * n2 -
                            1350 * b5 * d2 * e * f3 * g4 * n2 - 3 * b5 * d5 * g5 * n2 -
                            200 * b4 * c * d3 * f * g5 * n2 - 150 * a * b4 * d4 * f * g5 * n2 -
                            210 * b5 * d3 * f2 * g5 * n2 + 10 * c3 * e5 * n4 +
                            60 * a2 * c3 * e5 * n4 + 150 * a * c2 * d * e5 * n4 +
                            300 * a3 * c2 * d * e5 * n4 + 450 * a2 * c * d2 * e5 * n4 +
                            450 * a4 * c * d2 * e5 * n4 + 350 * a3 * d3 * e5 * n4 +
                            210 * a5 * d3 * e5 * n4 + 30 * b * c2 * e5 * f * n4 +
                            180 * a2 * b * c2 * e5 * f * n4 + 300 * a * b * c * d * e5 * f * n4 +
                            600 * a3 * b * c * d * e5 * f * n4 + 450 * a2 * b * d2 * e5 * f * n4 +
                            450 * a4 * b * d2 * e5 * f * n4 + 15 * c * e5 * f2 * n4 +
                            60 * a2 * c * e5 * f2 * n4 + 15 * a4 * c * e5 * f2 * n4 +
                            30 * b2 * c * e5 * f2 * n4 + 180 * a2 * b2 * c * e5 * f2 * n4 +
                            75 * a * d * e5 * f2 * n4 + 100 * a3 * d * e5 * f2 * n4 +
                            15 * a5 * d * e5 * f2 * n4 + 150 * a * b2 * d * e5 * f2 * n4 +
                            300 * a3 * b2 * d * e5 * f2 * n4 + 15 * b * e5 * f3 * n4 +
                            60 * a2 * b * e5 * f3 * n4 + 15 * a4 * b * e5 * f3 * n4 +
                            10 * b3 * e5 * f3 * n4 + 60 * a2 * b3 * e5 * f3 * n4 +
                            150 * a * b * c3 * e4 * g * n4 + 150 * b * c2 * d * e4 * g * n4 +
                            900 * a2 * b * c2 * d * e4 * g * n4 +
                            750 * a * b * c * d2 * e4 * g * n4 +
                            1500 * a3 * b * c * d2 * e4 * g * n4 + 750 * a2 * b * d3 * e4 * g * n4 +
                            750 * a4 * b * d3 * e4 * g * n4 + 300 * a * c2 * e4 * f * g * n4 +
                            150 * a3 * c2 * e4 * f * g * n4 + 450 * a * b2 * c2 * e4 * f * g * n4 +
                            300 * c * d * e4 * f * g * n4 + 1200 * a2 * c * d * e4 * f * g * n4 +
                            300 * a4 * c * d * e4 * f * g * n4 +
                            300 * b2 * c * d * e4 * f * g * n4 +
                            1800 * a2 * b2 * c * d * e4 * f * g * n4 +
                            750 * a * d2 * e4 * f * g * n4 + 1000 * a3 * d2 * e4 * f * g * n4 +
                            150 * a5 * d2 * e4 * f * g * n4 + 750 * a * b2 * d2 * e4 * f * g * n4 +
                            1500 * a3 * b2 * d2 * e4 * f * g * n4 +
                            750 * a * b * c * e4 * f2 * g * n4 +
                            375 * a3 * b * c * e4 * f2 * g * n4 +
                            450 * a * b3 * c * e4 * f2 * g * n4 + 375 * b * d * e4 * f2 * g * n4 +
                            1500 * a2 * b * d * e4 * f2 * g * n4 +
                            375 * a4 * b * d * e4 * f2 * g * n4 + 150 * b3 * d * e4 * f2 * g * n4 +
                            900 * a2 * b3 * d * e4 * f2 * g * n4 + 450 * a * b2 * e4 * f3 * g * n4 +
                            225 * a3 * b2 * e4 * f3 * g * n4 + 150 * a * b4 * e4 * f3 * g * n4 +
                            100 * b2 * c3 * e3 * g2 * n4 + 900 * a * b2 * c2 * d * e3 * g2 * n4 +
                            300 * b2 * c * d2 * e3 * g2 * n4 +
                            1800 * a2 * b2 * c * d2 * e3 * g2 * n4 +
                            500 * a * b2 * d3 * e3 * g2 * n4 + 1000 * a3 * b2 * d3 * e3 * g2 * n4 +
                            400 * b * c2 * e3 * f * g2 * n4 + 600 * a2 * b * c2 * e3 * f * g2 * n4 +
                            300 * b3 * c2 * e3 * f * g2 * n4 +
                            2400 * a * b * c * d * e3 * f * g2 * n4 +
                            1200 * a3 * b * c * d * e3 * f * g2 * n4 +
                            1800 * a * b3 * c * d * e3 * f * g2 * n4 +
                            600 * b * d2 * e3 * f * g2 * n4 +
                            2400 * a2 * b * d2 * e3 * f * g2 * n4 +
                            600 * a4 * b * d2 * e3 * f * g2 * n4 +
                            300 * b3 * d2 * e3 * f * g2 * n4 +
                            1800 * a2 * b3 * d2 * e3 * f * g2 * n4 +
                            900 * b2 * c * e3 * f2 * g2 * n4 +
                            1350 * a2 * b2 * c * e3 * f2 * g2 * n4 +
                            300 * b4 * c * e3 * f2 * g2 * n4 +
                            2700 * a * b2 * d * e3 * f2 * g2 * n4 +
                            1350 * a3 * b2 * d * e3 * f2 * g2 * n4 +
                            900 * a * b4 * d * e3 * f2 * g2 * n4 + 500 * b3 * e3 * f3 * g2 * n4 +
                            750 * a2 * b3 * e3 * f3 * g2 * n4 + 100 * b5 * e3 * f3 * g2 * n4 +
                            300 * b3 * c2 * d * e2 * g3 * n4 +
                            900 * a * b3 * c * d2 * e2 * g3 * n4 + 100 * b3 * d3 * e2 * g3 * n4 +
                            600 * a2 * b3 * d3 * e2 * g3 * n4 +
                            900 * a * b2 * c2 * e2 * f * g3 * n4 +
                            1200 * b2 * c * d * e2 * f * g3 * n4 +
                            1800 * a2 * b2 * c * d * e2 * f * g3 * n4 +
                            600 * b4 * c * d * e2 * f * g3 * n4 +
                            1800 * a * b2 * d2 * e2 * f * g3 * n4 +
                            900 * a3 * b2 * d2 * e2 * f * g3 * n4 +
                            900 * a * b4 * d2 * e2 * f * g3 * n4 +
                            1950 * a * b3 * c * e2 * f2 * g3 * n4 +
                            1300 * b3 * d * e2 * f2 * g3 * n4 +
                            1950 * a2 * b3 * d * e2 * f2 * g3 * n4 +
                            300 * b5 * d * e2 * f2 * g3 * n4 + 1050 * a * b4 * e2 * f3 * g3 * n4 +
                            150 * b4 * c * d2 * e * g4 * n4 + 150 * a * b4 * d3 * e * g4 * n4 +
                            600 * b3 * c2 * e * f * g4 * n4 +
                            1200 * a * b3 * c * d * e * f * g4 * n4 +
                            400 * b3 * d2 * e * f * g4 * n4 + 600 * a2 * b3 * d2 * e * f * g4 * n4 +
                            150 * b5 * d2 * e * f * g4 * n4 + 1275 * b4 * c * e * f2 * g4 * n4 +
                            1275 * a * b4 * d * e * f2 * g4 * n4 + 675 * b5 * e * f3 * g4 * n4 +
                            10 * b5 * d3 * g5 * n4 + 300 * b4 * c * d * f * g5 * n4 +
                            150 * a * b4 * d2 * f * g5 * n4 + 315 * b5 * d * f2 * g5 * n4 -
                            15 * c * e5 * n6 - 60 * a2 * c * e5 * n6 - 15 * a4 * c * e5 * n6 -
                            75 * a * d * e5 * n6 - 100 * a3 * d * e5 * n6 - 15 * a5 * d * e5 * n6 -
                            15 * b * e5 * f * n6 - 60 * a2 * b * e5 * f * n6 -
                            15 * a4 * b * e5 * f * n6 - 150 * a * b * c * e4 * g * n6 -
                            75 * a3 * b * c * e4 * g * n6 - 75 * b * d * e4 * g * n6 -
                            300 * a2 * b * d * e4 * g * n6 - 75 * a4 * b * d * e4 * g * n6 -
                            150 * a * e4 * f * g * n6 - 50 * a3 * e4 * f * g * n6 -
                            150 * a * b2 * e4 * f * g * n6 - 75 * a3 * b2 * e4 * f * g * n6 -
                            100 * b2 * c * e3 * g2 * n6 - 150 * a2 * b2 * c * e3 * g2 * n6 -
                            300 * a * b2 * d * e3 * g2 * n6 - 150 * a3 * b2 * d * e3 * g2 * n6 -
                            200 * b * e3 * f * g2 * n6 - 200 * a2 * b * e3 * f * g2 * n6 -
                            100 * b3 * e3 * f * g2 * n6 - 150 * a2 * b3 * e3 * f * g2 * n6 -
                            150 * a * b3 * c * e2 * g3 * n6 - 100 * b3 * d * e2 * g3 * n6 -
                            150 * a2 * b3 * d * e2 * g3 * n6 - 300 * a * b2 * e2 * f * g3 * n6 -
                            150 * a * b4 * e2 * f * g3 * n6 - 75 * b4 * c * e * g4 * n6 -
                            75 * a * b4 * d * e * g4 * n6 - 200 * b3 * e * f * g4 * n6 -
                            75 * b5 * e * f * g4 * n6 - 15 * b5 * d * g5 * n6) +
                    2 * (75 * a * c4 * e6 * f2 + 900 * a2 * c3 * d * e6 * f2 +
                            3150 * a3 * c2 * d2 * e6 * f2 + 4200 * a4 * c * d3 * e6 * f2 +
                            1890 * a5 * d4 * e6 * f2 + 300 * a * b * c3 * e6 * f3 +
                            2700 * a2 * b * c2 * d * e6 * f3 + 6300 * a3 * b * c * d2 * e6 * f3 +
                            4200 * a4 * b * d3 * e6 * f3 + 450 * a * b2 * c2 * e6 * f4 +
                            2700 * a2 * b2 * c * d * e6 * f4 + 3150 * a3 * b2 * d2 * e6 * f4 +
                            300 * a * b3 * c * e6 * f5 + 900 * a2 * b3 * d * e6 * f5 +
                            75 * a * b4 * e6 * f6 + 72 * c5 * e5 * f * g +
                            1800 * a * c4 * d * e5 * f * g + 10800 * a2 * c3 * d2 * e5 * f * g +
                            25200 * a3 * c2 * d3 * e5 * f * g + 25200 * a4 * c * d4 * e5 * f * g +
                            9072 * a5 * d5 * e5 * f * g + 450 * b * c4 * e5 * f2 * g +
                            9000 * a * b * c3 * d * e5 * f2 * g +
                            40500 * a2 * b * c2 * d2 * e5 * f2 * g +
                            63000 * a3 * b * c * d3 * e5 * f2 * g +
                            31500 * a4 * b * d4 * e5 * f2 * g + 1080 * b2 * c3 * e5 * f3 * g +
                            16200 * a * b2 * c2 * d * e5 * f3 * g +
                            48600 * a2 * b2 * c * d2 * e5 * f3 * g +
                            37800 * a3 * b2 * d3 * e5 * f3 * g + 1260 * b3 * c2 * e5 * f4 * g +
                            12600 * a * b3 * c * d * e5 * f4 * g +
                            18900 * a2 * b3 * d2 * e5 * f4 * g + 720 * b4 * c * e5 * f5 * g +
                            3600 * a * b4 * d * e5 * f5 * g + 162 * b5 * e5 * f6 * g +
                            1800 * b * c4 * d * e4 * f * g2 +
                            18000 * a * b * c3 * d2 * e4 * f * g2 +
                            54000 * a2 * b * c2 * d3 * e4 * f * g2 +
                            63000 * a3 * b * c * d4 * e4 * f * g2 +
                            25200 * a4 * b * d5 * e4 * f * g2 + 8100 * b2 * c3 * d * e4 * f2 * g2 +
                            60750 * a * b2 * c2 * d2 * e4 * f2 * g2 +
                            121500 * a2 * b2 * c * d3 * e4 * f2 * g2 +
                            70875 * a3 * b2 * d4 * e4 * f2 * g2 +
                            13500 * b3 * c2 * d * e4 * f3 * g2 +
                            67500 * a * b3 * c * d2 * e4 * f3 * g2 +
                            67500 * a2 * b3 * d3 * e4 * f3 * g2 + 9900 * b4 * c * d * e4 * f4 * g2 +
                            24750 * a * b4 * d2 * e4 * f4 * g2 + 2700 * b5 * d * e4 * f5 * g2 +
                            7200 * b2 * c3 * d2 * e3 * f * g3 +
                            36000 * a * b2 * c2 * d3 * e3 * f * g3 +
                            54000 * a2 * b2 * c * d4 * e3 * f * g3 +
                            25200 * a3 * b2 * d5 * e3 * f * g3 +
                            23400 * b3 * c2 * d2 * e3 * f2 * g3 +
                            78000 * a * b3 * c * d3 * e3 * f2 * g3 +
                            58500 * a2 * b3 * d4 * e3 * f2 * g3 +
                            25200 * b4 * c * d2 * e3 * f3 * g3 +
                            42000 * a * b4 * d3 * e3 * f3 * g3 + 9000 * b5 * d2 * e3 * f4 * g3 +
                            7200 * b3 * c2 * d3 * e2 * f * g4 +
                            18000 * a * b3 * c * d4 * e2 * f * g4 +
                            10800 * a2 * b3 * d5 * e2 * f * g4 +
                            15300 * b4 * c * d3 * e2 * f2 * g4 +
                            19125 * a * b4 * d4 * e2 * f2 * g4 + 8100 * b5 * d3 * e2 * f3 * g4 +
                            1800 * b4 * c * d4 * e * f * g5 + 1800 * a * b4 * d5 * e * f * g5 +
                            1890 * b5 * d4 * e * f2 * g5 + 72 * b5 * d5 * f * g6 -
                            75 * a * c4 * e6 * n2 - 900 * a2 * c3 * d * e6 * n2 -
                            3150 * a3 * c2 * d2 * e6 * n2 - 4200 * a4 * c * d3 * e6 * n2 -
                            1890 * a5 * d4 * e6 * n2 - 300 * a * b * c3 * e6 * f * n2 -
                            2700 * a2 * b * c2 * d * e6 * f * n2 -
                            6300 * a3 * b * c * d2 * e6 * f * n2 -
                            4200 * a4 * b * d3 * e6 * f * n2 - 150 * a * c2 * e6 * f2 * n2 -
                            300 * a3 * c2 * e6 * f2 * n2 - 450 * a * b2 * c2 * e6 * f2 * n2 -
                            900 * a2 * c * d * e6 * f2 * n2 - 900 * a4 * c * d * e6 * f2 * n2 -
                            2700 * a2 * b2 * c * d * e6 * f2 * n2 - 1050 * a3 * d2 * e6 * f2 * n2 -
                            630 * a5 * d2 * e6 * f2 * n2 - 3150 * a3 * b2 * d2 * e6 * f2 * n2 -
                            300 * a * b * c * e6 * f3 * n2 - 600 * a3 * b * c * e6 * f3 * n2 -
                            300 * a * b3 * c * e6 * f3 * n2 - 900 * a2 * b * d * e6 * f3 * n2 -
                            900 * a4 * b * d * e6 * f3 * n2 - 900 * a2 * b3 * d * e6 * f3 * n2 -
                            150 * a * b2 * e6 * f4 * n2 - 300 * a3 * b2 * e6 * f4 * n2 -
                            75 * a * b4 * e6 * f4 * n2 - 90 * b * c4 * e5 * g * n2 -
                            1800 * a * b * c3 * d * e5 * g * n2 -
                            8100 * a2 * b * c2 * d2 * e5 * g * n2 -
                            12600 * a3 * b * c * d3 * e5 * g * n2 -
                            6300 * a4 * b * d4 * e5 * g * n2 - 240 * c3 * e5 * f * g * n2 -
                            1440 * a2 * c3 * e5 * f * g * n2 - 360 * b2 * c3 * e5 * f * g * n2 -
                            3600 * a * c2 * d * e5 * f * g * n2 -
                            7200 * a3 * c2 * d * e5 * f * g * n2 -
                            5400 * a * b2 * c2 * d * e5 * f * g * n2 -
                            10800 * a2 * c * d2 * e5 * f * g * n2 -
                            10800 * a4 * c * d2 * e5 * f * g * n2 -
                            16200 * a2 * b2 * c * d2 * e5 * f * g * n2 -
                            8400 * a3 * d3 * e5 * f * g * n2 - 5040 * a5 * d3 * e5 * f * g * n2 -
                            12600 * a3 * b2 * d3 * e5 * f * g * n2 -
                            900 * b * c2 * e5 * f2 * g * n2 -
                            5400 * a2 * b * c2 * e5 * f2 * g * n2 -
                            540 * b3 * c2 * e5 * f2 * g * n2 -
                            9000 * a * b * c * d * e5 * f2 * g * n2 -
                            18000 * a3 * b * c * d * e5 * f2 * g * n2 -
                            5400 * a * b3 * c * d * e5 * f2 * g * n2 -
                            13500 * a2 * b * d2 * e5 * f2 * g * n2 -
                            13500 * a4 * b * d2 * e5 * f2 * g * n2 -
                            8100 * a2 * b3 * d2 * e5 * f2 * g * n2 -
                            1080 * b2 * c * e5 * f3 * g * n2 -
                            6480 * a2 * b2 * c * e5 * f3 * g * n2 -
                            360 * b4 * c * e5 * f3 * g * n2 - 5400 * a * b2 * d * e5 * f3 * g * n2 -
                            10800 * a3 * b2 * d * e5 * f3 * g * n2 -
                            1800 * a * b4 * d * e5 * f3 * g * n2 - 420 * b3 * e5 * f4 * g * n2 -
                            2520 * a2 * b3 * e5 * f4 * g * n2 - 90 * b5 * e5 * f4 * g * n2 -
                            900 * b2 * c3 * d * e4 * g2 * n2 -
                            6750 * a * b2 * c2 * d2 * e4 * g2 * n2 -
                            13500 * a2 * b2 * c * d3 * e4 * g2 * n2 -
                            7875 * a3 * b2 * d4 * e4 * g2 * n2 -
                            3600 * a * b * c3 * e4 * f * g2 * n2 -
                            3600 * b * c2 * d * e4 * f * g2 * n2 -
                            21600 * a2 * b * c2 * d * e4 * f * g2 * n2 -
                            2700 * b3 * c2 * d * e4 * f * g2 * n2 -
                            18000 * a * b * c * d2 * e4 * f * g2 * n2 -
                            36000 * a3 * b * c * d2 * e4 * f * g2 * n2 -
                            13500 * a * b3 * c * d2 * e4 * f * g2 * n2 -
                            18000 * a2 * b * d3 * e4 * f * g2 * n2 -
                            18000 * a4 * b * d3 * e4 * f * g2 * n2 -
                            13500 * a2 * b3 * d3 * e4 * f * g2 * n2 -
                            12150 * a * b2 * c2 * e4 * f2 * g2 * n2 -
                            8100 * b2 * c * d * e4 * f2 * g2 * n2 -
                            48600 * a2 * b2 * c * d * e4 * f2 * g2 * n2 -
                            2700 * b4 * c * d * e4 * f2 * g2 * n2 -
                            20250 * a * b2 * d2 * e4 * f2 * g2 * n2 -
                            40500 * a3 * b2 * d2 * e4 * f2 * g2 * n2 -
                            6750 * a * b4 * d2 * e4 * f2 * g2 * n2 -
                            13500 * a * b3 * c * e4 * f3 * g2 * n2 -
                            4500 * b3 * d * e4 * f3 * g2 * n2 -
                            27000 * a2 * b3 * d * e4 * f3 * g2 * n2 -
                            900 * b5 * d * e4 * f3 * g2 * n2 - 4950 * a * b4 * e4 * f4 * g2 * n2 -
                            1800 * b3 * c2 * d2 * e3 * g3 * n2 -
                            6000 * a * b3 * c * d3 * e3 * g3 * n2 -
                            4500 * a2 * b3 * d4 * e3 * g3 * n2 - 2400 * b2 * c3 * e3 * f * g3 * n2 -
                            21600 * a * b2 * c2 * d * e3 * f * g3 * n2 -
                            7200 * b2 * c * d2 * e3 * f * g3 * n2 -
                            43200 * a2 * b2 * c * d2 * e3 * f * g3 * n2 -
                            3600 * b4 * c * d2 * e3 * f * g3 * n2 -
                            12000 * a * b2 * d3 * e3 * f * g3 * n2 -
                            24000 * a3 * b2 * d3 * e3 * f * g3 * n2 -
                            6000 * a * b4 * d3 * e3 * f * g3 * n2 -
                            7800 * b3 * c2 * e3 * f2 * g3 * n2 -
                            46800 * a * b3 * c * d * e3 * f2 * g3 * n2 -
                            7800 * b3 * d2 * e3 * f2 * g3 * n2 -
                            46800 * a2 * b3 * d2 * e3 * f2 * g3 * n2 -
                            1800 * b5 * d2 * e3 * f2 * g3 * n2 - 8400 * b4 * c * e3 * f3 * g3 * n2 -
                            25200 * a * b4 * d * e3 * f3 * g3 * n2 - 3000 * b5 * e3 * f4 * g3 * n2 -
                            900 * b4 * c * d3 * e2 * g4 * n2 - 1125 * a * b4 * d4 * e2 * g4 * n2 -
                            7200 * b3 * c2 * d * e2 * f * g4 * n2 -
                            21600 * a * b3 * c * d2 * e2 * f * g4 * n2 -
                            2400 * b3 * d3 * e2 * f * g4 * n2 -
                            14400 * a2 * b3 * d3 * e2 * f * g4 * n2 -
                            900 * b5 * d3 * e2 * f * g4 * n2 -
                            15300 * b4 * c * d * e2 * f2 * g4 * n2 -
                            22950 * a * b4 * d2 * e2 * f2 * g4 * n2 -
                            8100 * b5 * d * e2 * f3 * g4 * n2 - 90 * b5 * d4 * e * g5 * n2 -
                            3600 * b4 * c * d2 * e * f * g5 * n2 -
                            3600 * a * b4 * d3 * e * f * g5 * n2 -
                            3780 * b5 * d2 * e * f2 * g5 * n2 - 240 * b5 * d3 * f * g6 * n2 +
                            150 * a * c2 * e6 * n4 + 300 * a3 * c2 * e6 * n4 +
                            900 * a2 * c * d * e6 * n4 + 900 * a4 * c * d * e6 * n4 +
                            1050 * a3 * d2 * e6 * n4 + 630 * a5 * d2 * e6 * n4 +
                            300 * a * b * c * e6 * f * n4 + 600 * a3 * b * c * e6 * f * n4 +
                            900 * a2 * b * d * e6 * f * n4 + 900 * a4 * b * d * e6 * f * n4 +
                            75 * a * e6 * f2 * n4 + 100 * a3 * e6 * f2 * n4 +
                            15 * a5 * e6 * f2 * n4 + 150 * a * b2 * e6 * f2 * n4 +
                            300 * a3 * b2 * e6 * f2 * n4 + 180 * b * c2 * e5 * g * n4 +
                            1080 * a2 * b * c2 * e5 * g * n4 + 1800 * a * b * c * d * e5 * g * n4 +
                            3600 * a3 * b * c * d * e5 * g * n4 + 2700 * a2 * b * d2 * e5 * g * n4 +
                            2700 * a4 * b * d2 * e5 * g * n4 + 360 * c * e5 * f * g * n4 +
                            1440 * a2 * c * e5 * f * g * n4 + 360 * a4 * c * e5 * f * g * n4 +
                            360 * b2 * c * e5 * f * g * n4 + 2160 * a2 * b2 * c * e5 * f * g * n4 +
                            1800 * a * d * e5 * f * g * n4 + 2400 * a3 * d * e5 * f * g * n4 +
                            360 * a5 * d * e5 * f * g * n4 + 1800 * a * b2 * d * e5 * f * g * n4 +
                            3600 * a3 * b2 * d * e5 * f * g * n4 + 450 * b * e5 * f2 * g * n4 +
                            1800 * a2 * b * e5 * f2 * g * n4 + 450 * a4 * b * e5 * f2 * g * n4 +
                            180 * b3 * e5 * f2 * g * n4 + 1080 * a2 * b3 * e5 * f2 * g * n4 +
                            1350 * a * b2 * c2 * e4 * g2 * n4 + 900 * b2 * c * d * e4 * g2 * n4 +
                            5400 * a2 * b2 * c * d * e4 * g2 * n4 +
                            2250 * a * b2 * d2 * e4 * g2 * n4 + 4500 * a3 * b2 * d2 * e4 * g2 * n4 +
                            3600 * a * b * c * e4 * f * g2 * n4 +
                            1800 * a3 * b * c * e4 * f * g2 * n4 +
                            2700 * a * b3 * c * e4 * f * g2 * n4 + 1800 * b * d * e4 * f * g2 * n4 +
                            7200 * a2 * b * d * e4 * f * g2 * n4 +
                            1800 * a4 * b * d * e4 * f * g2 * n4 + 900 * b3 * d * e4 * f * g2 * n4 +
                            5400 * a2 * b3 * d * e4 * f * g2 * n4 +
                            4050 * a * b2 * e4 * f2 * g2 * n4 + 2025 * a3 * b2 * e4 * f2 * g2 * n4 +
                            1350 * a * b4 * e4 * f2 * g2 * n4 + 600 * b3 * c2 * e3 * g3 * n4 +
                            3600 * a * b3 * c * d * e3 * g3 * n4 + 600 * b3 * d2 * e3 * g3 * n4 +
                            3600 * a2 * b3 * d2 * e3 * g3 * n4 + 2400 * b2 * c * e3 * f * g3 * n4 +
                            3600 * a2 * b2 * c * e3 * f * g3 * n4 +
                            1200 * b4 * c * e3 * f * g3 * n4 +
                            7200 * a * b2 * d * e3 * f * g3 * n4 +
                            3600 * a3 * b2 * d * e3 * f * g3 * n4 +
                            3600 * a * b4 * d * e3 * f * g3 * n4 + 2600 * b3 * e3 * f2 * g3 * n4 +
                            3900 * a2 * b3 * e3 * f2 * g3 * n4 + 600 * b5 * e3 * f2 * g3 * n4 +
                            900 * b4 * c * d * e2 * g4 * n4 + 1350 * a * b4 * d2 * e2 * g4 * n4 +
                            3600 * a * b3 * c * e2 * f * g4 * n4 +
                            2400 * b3 * d * e2 * f * g4 * n4 +
                            3600 * a2 * b3 * d * e2 * f * g4 * n4 +
                            900 * b5 * d * e2 * f * g4 * n4 + 3825 * a * b4 * e2 * f2 * g4 * n4 +
                            180 * b5 * d2 * e * g5 * n4 + 1800 * b4 * c * e * f * g5 * n4 +
                            1800 * a * b4 * d * e * f * g5 * n4 + 1890 * b5 * e * f2 * g5 * n4 +
                            360 * b5 * d * f * g6 * n4 - 75 * a * e6 * n6 - 100 * a3 * e6 * n6 -
                            15 * a5 * e6 * n6 - 90 * b * e5 * g * n6 - 360 * a2 * b * e5 * g * n6 -
                            90 * a4 * b * e5 * g * n6 - 450 * a * b2 * e4 * g2 * n6 -
                            225 * a3 * b2 * e4 * g2 * n6 - 200 * b3 * e3 * g3 * n6 -
                            300 * a2 * b3 * e3 * g3 * n6 - 225 * a * b4 * e2 * g4 * n6 -
                            90 * b5 * e * g5 * n6) +
                    5 * (45 * a2 * c3 * e7 * f2 + 315 * a3 * c2 * d * e7 * f2 +
                            630 * a4 * c * d2 * e7 * f2 + 378 * a5 * d3 * e7 * f2 +
                            135 * a2 * b * c2 * e7 * f3 + 630 * a3 * b * c * d * e7 * f3 +
                            630 * a4 * b * d2 * e7 * f3 + 135 * a2 * b2 * c * e7 * f4 +
                            315 * a3 * b2 * d * e7 * f4 + 45 * a2 * b3 * e7 * f5 +
                            105 * a * c4 * e6 * f * g + 1260 * a2 * c3 * d * e6 * f * g +
                            4410 * a3 * c2 * d2 * e6 * f * g + 5880 * a4 * c * d3 * e6 * f * g +
                            2646 * a5 * d4 * e6 * f * g + 525 * a * b * c3 * e6 * f2 * g +
                            4725 * a2 * b * c2 * d * e6 * f2 * g +
                            11025 * a3 * b * c * d2 * e6 * f2 * g +
                            7350 * a4 * b * d3 * e6 * f2 * g + 945 * a * b2 * c2 * e6 * f3 * g +
                            5670 * a2 * b2 * c * d * e6 * f3 * g +
                            6615 * a3 * b2 * d2 * e6 * f3 * g + 735 * a * b3 * c * e6 * f4 * g +
                            2205 * a2 * b3 * d * e6 * f4 * g + 210 * a * b4 * e6 * f5 * g +
                            126 * b * c4 * e5 * f * g2 + 2520 * a * b * c3 * d * e5 * f * g2 +
                            11340 * a2 * b * c2 * d2 * e5 * f * g2 +
                            17640 * a3 * b * c * d3 * e5 * f * g2 +
                            8820 * a4 * b * d4 * e5 * f * g2 + 567 * b2 * c3 * e5 * f2 * g2 +
                            8505 * a * b2 * c2 * d * e5 * f2 * g2 +
                            25515 * a2 * b2 * c * d2 * e5 * f2 * g2 +
                            19845 * a3 * b2 * d3 * e5 * f2 * g2 + 945 * b3 * c2 * e5 * f3 * g2 +
                            9450 * a * b3 * c * d * e5 * f3 * g2 +
                            14175 * a2 * b3 * d2 * e5 * f3 * g2 + 693 * b4 * c * e5 * f4 * g2 +
                            3465 * a * b4 * d * e5 * f4 * g2 + 189 * b5 * e5 * f5 * g2 +
                            1260 * b2 * c3 * d * e4 * f * g3 +
                            9450 * a * b2 * c2 * d2 * e4 * f * g3 +
                            18900 * a2 * b2 * c * d3 * e4 * f * g3 +
                            11025 * a3 * b2 * d4 * e4 * f * g3 + 4095 * b3 * c2 * d * e4 * f2 * g3 +
                            20475 * a * b3 * c * d2 * e4 * f2 * g3 +
                            20475 * a2 * b3 * d3 * e4 * f2 * g3 + 4410 * b4 * c * d * e4 * f3 * g3 +
                            11025 * a * b4 * d2 * e4 * f3 * g3 + 1575 * b5 * d * e4 * f4 * g3 +
                            2520 * b3 * c2 * d2 * e3 * f * g4 +
                            8400 * a * b3 * c * d3 * e3 * f * g4 +
                            6300 * a2 * b3 * d4 * e3 * f * g4 + 5355 * b4 * c * d2 * e3 * f2 * g4 +
                            8925 * a * b4 * d3 * e3 * f2 * g4 + 2835 * b5 * d2 * e3 * f3 * g4 +
                            1260 * b4 * c * d3 * e2 * f * g5 + 1575 * a * b4 * d4 * e2 * f * g5 +
                            1323 * b5 * d3 * e2 * f2 * g5 + 126 * b5 * d4 * e * f * g6 -
                            45 * a2 * c3 * e7 * n2 - 315 * a3 * c2 * d * e7 * n2 -
                            630 * a4 * c * d2 * e7 * n2 - 378 * a5 * d3 * e7 * n2 -
                            135 * a2 * b * c2 * e7 * f * n2 - 630 * a3 * b * c * d * e7 * f * n2 -
                            630 * a4 * b * d2 * e7 * f * n2 - 45 * a2 * c * e7 * f2 * n2 -
                            45 * a4 * c * e7 * f2 * n2 - 135 * a2 * b2 * c * e7 * f2 * n2 -
                            105 * a3 * d * e7 * f2 * n2 - 63 * a5 * d * e7 * f2 * n2 -
                            315 * a3 * b2 * d * e7 * f2 * n2 - 45 * a2 * b * e7 * f3 * n2 -
                            45 * a4 * b * e7 * f3 * n2 - 45 * a2 * b3 * e7 * f3 * n2 -
                            105 * a * b * c3 * e6 * g * n2 - 945 * a2 * b * c2 * d * e6 * g * n2 -
                            2205 * a3 * b * c * d2 * e6 * g * n2 -
                            1470 * a4 * b * d3 * e6 * g * n2 - 210 * a * c2 * e6 * f * g * n2 -
                            420 * a3 * c2 * e6 * f * g * n2 - 315 * a * b2 * c2 * e6 * f * g * n2 -
                            1260 * a2 * c * d * e6 * f * g * n2 -
                            1260 * a4 * c * d * e6 * f * g * n2 -
                            1890 * a2 * b2 * c * d * e6 * f * g * n2 -
                            1470 * a3 * d2 * e6 * f * g * n2 - 882 * a5 * d2 * e6 * f * g * n2 -
                            2205 * a3 * b2 * d2 * e6 * f * g * n2 -
                            525 * a * b * c * e6 * f2 * g * n2 -
                            1050 * a3 * b * c * e6 * f2 * g * n2 -
                            315 * a * b3 * c * e6 * f2 * g * n2 -
                            1575 * a2 * b * d * e6 * f2 * g * n2 -
                            1575 * a4 * b * d * e6 * f2 * g * n2 -
                            945 * a2 * b3 * d * e6 * f2 * g * n2 - 315 * a * b2 * e6 * f3 * g * n2 -
                            630 * a3 * b2 * e6 * f3 * g * n2 - 105 * a * b4 * e6 * f3 * g * n2 -
                            63 * b2 * c3 * e5 * g2 * n2 - 945 * a * b2 * c2 * d * e5 * g2 * n2 -
                            2835 * a2 * b2 * c * d2 * e5 * g2 * n2 -
                            2205 * a3 * b2 * d3 * e5 * g2 * n2 - 252 * b * c2 * e5 * f * g2 * n2 -
                            1512 * a2 * b * c2 * e5 * f * g2 * n2 -
                            189 * b3 * c2 * e5 * f * g2 * n2 -
                            2520 * a * b * c * d * e5 * f * g2 * n2 -
                            5040 * a3 * b * c * d * e5 * f * g2 * n2 -
                            1890 * a * b3 * c * d * e5 * f * g2 * n2 -
                            3780 * a2 * b * d2 * e5 * f * g2 * n2 -
                            3780 * a4 * b * d2 * e5 * f * g2 * n2 -
                            2835 * a2 * b3 * d2 * e5 * f * g2 * n2 -
                            567 * b2 * c * e5 * f2 * g2 * n2 -
                            3402 * a2 * b2 * c * e5 * f2 * g2 * n2 -
                            189 * b4 * c * e5 * f2 * g2 * n2 -
                            2835 * a * b2 * d * e5 * f2 * g2 * n2 -
                            5670 * a3 * b2 * d * e5 * f2 * g2 * n2 -
                            945 * a * b4 * d * e5 * f2 * g2 * n2 - 315 * b3 * e5 * f3 * g2 * n2 -
                            1890 * a2 * b3 * e5 * f3 * g2 * n2 - 63 * b5 * e5 * f3 * g2 * n2 -
                            315 * b3 * c2 * d * e4 * g3 * n2 -
                            1575 * a * b3 * c * d2 * e4 * g3 * n2 -
                            1575 * a2 * b3 * d3 * e4 * g3 * n2 -
                            1890 * a * b2 * c2 * e4 * f * g3 * n2 -
                            1260 * b2 * c * d * e4 * f * g3 * n2 -
                            7560 * a2 * b2 * c * d * e4 * f * g3 * n2 -
                            630 * b4 * c * d * e4 * f * g3 * n2 -
                            3150 * a * b2 * d2 * e4 * f * g3 * n2 -
                            6300 * a3 * b2 * d2 * e4 * f * g3 * n2 -
                            1575 * a * b4 * d2 * e4 * f * g3 * n2 -
                            4095 * a * b3 * c * e4 * f2 * g3 * n2 -
                            1365 * b3 * d * e4 * f2 * g3 * n2 -
                            8190 * a2 * b3 * d * e4 * f2 * g3 * n2 -
                            315 * b5 * d * e4 * f2 * g3 * n2 - 2205 * a * b4 * e4 * f3 * g3 * n2 -
                            315 * b4 * c * d2 * e3 * g4 * n2 - 525 * a * b4 * d3 * e3 * g4 * n2 -
                            840 * b3 * c2 * e3 * f * g4 * n2 -
                            5040 * a * b3 * c * d * e3 * f * g4 * n2 -
                            840 * b3 * d2 * e3 * f * g4 * n2 -
                            5040 * a2 * b3 * d2 * e3 * f * g4 * n2 -
                            315 * b5 * d2 * e3 * f * g4 * n2 - 1785 * b4 * c * e3 * f2 * g4 * n2 -
                            5355 * a * b4 * d * e3 * f2 * g4 * n2 - 945 * b5 * e3 * f3 * g4 * n2 -
                            63 * b5 * d3 * e2 * g5 * n2 - 1260 * b4 * c * d * e2 * f * g5 * n2 -
                            1890 * a * b4 * d2 * e2 * f * g5 * n2 -
                            1323 * b5 * d * e2 * f2 * g5 * n2 - 252 * b5 * d2 * e * f * g6 * n2 +
                            45 * a2 * c * e7 * n4 + 45 * a4 * c * e7 * n4 + 105 * a3 * d * e7 * n4 +
                            63 * a5 * d * e7 * n4 + 45 * a2 * b * e7 * f * n4 +
                            45 * a4 * b * e7 * f * n4 + 105 * a * b * c * e6 * g * n4 +
                            210 * a3 * b * c * e6 * g * n4 + 315 * a2 * b * d * e6 * g * n4 +
                            315 * a4 * b * d * e6 * g * n4 + 105 * a * e6 * f * g * n4 +
                            140 * a3 * e6 * f * g * n4 + 21 * a5 * e6 * f * g * n4 +
                            105 * a * b2 * e6 * f * g * n4 + 210 * a3 * b2 * e6 * f * g * n4 +
                            63 * b2 * c * e5 * g2 * n4 + 378 * a2 * b2 * c * e5 * g2 * n4 +
                            315 * a * b2 * d * e5 * g2 * n4 + 630 * a3 * b2 * d * e5 * g2 * n4 +
                            126 * b * e5 * f * g2 * n4 + 504 * a2 * b * e5 * f * g2 * n4 +
                            126 * a4 * b * e5 * f * g2 * n4 + 63 * b3 * e5 * f * g2 * n4 +
                            378 * a2 * b3 * e5 * f * g2 * n4 + 315 * a * b3 * c * e4 * g3 * n4 +
                            105 * b3 * d * e4 * g3 * n4 + 630 * a2 * b3 * d * e4 * g3 * n4 +
                            630 * a * b2 * e4 * f * g3 * n4 + 315 * a3 * b2 * e4 * f * g3 * n4 +
                            315 * a * b4 * e4 * f * g3 * n4 + 105 * b4 * c * e3 * g4 * n4 +
                            315 * a * b4 * d * e3 * g4 * n4 + 280 * b3 * e3 * f * g4 * n4 +
                            420 * a2 * b3 * e3 * f * g4 * n4 + 105 * b5 * e3 * f * g4 * n4 +
                            63 * b5 * d * e2 * g5 * n4 + 315 * a * b4 * e2 * f * g5 * n4 +
                            126 * b5 * e * f * g6 * n4) +
                    (5 *
                        (105 * a3 * c2 * e8 * f2 + 420 * a4 * c * d * e8 * f2 +
                            378 * a5 * d2 * e8 * f2 + 210 * a3 * b * c * e8 * f3 +
                            420 * a4 * b * d * e8 * f3 + 105 * a3 * b2 * e8 * f4 +
                            480 * a2 * c3 * e7 * f * g + 3360 * a3 * c2 * d * e7 * f * g +
                            6720 * a4 * c * d2 * e7 * f * g + 4032 * a5 * d3 * e7 * f * g +
                            1800 * a2 * b * c2 * e7 * f2 * g + 8400 * a3 * b * c * d * e7 * f2 * g +
                            8400 * a4 * b * d2 * e7 * f2 * g + 2160 * a2 * b2 * c * e7 * f3 * g +
                            5040 * a3 * b2 * d * e7 * f3 * g + 840 * a2 * b3 * e7 * f4 * g +
                            1120 * a * b * c3 * e6 * f * g2 +
                            10080 * a2 * b * c2 * d * e6 * f * g2 +
                            23520 * a3 * b * c * d2 * e6 * f * g2 +
                            15680 * a4 * b * d3 * e6 * f * g2 + 3780 * a * b2 * c2 * e6 * f2 * g2 +
                            22680 * a2 * b2 * c * d * e6 * f2 * g2 +
                            26460 * a3 * b2 * d2 * e6 * f2 * g2 + 4200 * a * b3 * c * e6 * f3 * g2 +
                            12600 * a2 * b3 * d * e6 * f3 * g2 + 1540 * a * b4 * e6 * f4 * g2 +
                            672 * b2 * c3 * e5 * f * g3 + 10080 * a * b2 * c2 * d * e5 * f * g3 +
                            30240 * a2 * b2 * c * d2 * e5 * f * g3 +
                            23520 * a3 * b2 * d3 * e5 * f * g3 + 2184 * b3 * c2 * e5 * f2 * g3 +
                            21840 * a * b3 * c * d * e5 * f2 * g3 +
                            32760 * a2 * b3 * d2 * e5 * f2 * g3 + 2352 * b4 * c * e5 * f3 * g3 +
                            11760 * a * b4 * d * e5 * f3 * g3 + 840 * b5 * e5 * f4 * g3 +
                            3360 * b3 * c2 * d * e4 * f * g4 +
                            16800 * a * b3 * c * d2 * e4 * f * g4 +
                            16800 * a2 * b3 * d3 * e4 * f * g4 + 7140 * b4 * c * d * e4 * f2 * g4 +
                            17850 * a * b4 * d2 * e4 * f2 * g4 + 3780 * b5 * d * e4 * f3 * g4 +
                            3360 * b4 * c * d2 * e3 * f * g5 + 5600 * a * b4 * d3 * e3 * f * g5 +
                            3528 * b5 * d2 * e3 * f2 * g5 + 672 * b5 * d3 * e2 * f * g6 -
                            105 * a3 * c2 * e8 * n2 - 420 * a4 * c * d * e8 * n2 -
                            378 * a5 * d2 * e8 * n2 - 210 * a3 * b * c * e8 * f * n2 -
                            420 * a4 * b * d * e8 * f * n2 - 35 * a3 * e8 * f2 * n2 -
                            21 * a5 * e8 * f2 * n2 - 105 * a3 * b2 * e8 * f2 * n2 -
                            360 * a2 * b * c2 * e7 * g * n2 - 1680 * a3 * b * c * d * e7 * g * n2 -
                            1680 * a4 * b * d2 * e7 * g * n2 - 480 * a2 * c * e7 * f * g * n2 -
                            480 * a4 * c * e7 * f * g * n2 - 720 * a2 * b2 * c * e7 * f * g * n2 -
                            1120 * a3 * d * e7 * f * g * n2 - 672 * a5 * d * e7 * f * g * n2 -
                            1680 * a3 * b2 * d * e7 * f * g * n2 - 600 * a2 * b * e7 * f2 * g * n2 -
                            600 * a4 * b * e7 * f2 * g * n2 - 360 * a2 * b3 * e7 * f2 * g * n2 -
                            420 * a * b2 * c2 * e6 * g2 * n2 -
                            2520 * a2 * b2 * c * d * e6 * g2 * n2 -
                            2940 * a3 * b2 * d2 * e6 * g2 * n2 -
                            1120 * a * b * c * e6 * f * g2 * n2 -
                            2240 * a3 * b * c * e6 * f * g2 * n2 -
                            840 * a * b3 * c * e6 * f * g2 * n2 -
                            3360 * a2 * b * d * e6 * f * g2 * n2 -
                            3360 * a4 * b * d * e6 * f * g2 * n2 -
                            2520 * a2 * b3 * d * e6 * f * g2 * n2 -
                            1260 * a * b2 * e6 * f2 * g2 * n2 - 2520 * a3 * b2 * e6 * f2 * g2 * n2 -
                            420 * a * b4 * e6 * f2 * g2 * n2 - 168 * b3 * c2 * e5 * g3 * n2 -
                            1680 * a * b3 * c * d * e5 * g3 * n2 -
                            2520 * a2 * b3 * d2 * e5 * g3 * n2 - 672 * b2 * c * e5 * f * g3 * n2 -
                            4032 * a2 * b2 * c * e5 * f * g3 * n2 -
                            336 * b4 * c * e5 * f * g3 * n2 - 3360 * a * b2 * d * e5 * f * g3 * n2 -
                            6720 * a3 * b2 * d * e5 * f * g3 * n2 -
                            1680 * a * b4 * d * e5 * f * g3 * n2 - 728 * b3 * e5 * f2 * g3 * n2 -
                            4368 * a2 * b3 * e5 * f2 * g3 * n2 - 168 * b5 * e5 * f2 * g3 * n2 -
                            420 * b4 * c * d * e4 * g4 * n2 - 1050 * a * b4 * d2 * e4 * g4 * n2 -
                            3360 * a * b3 * c * e4 * f * g4 * n2 -
                            1120 * b3 * d * e4 * f * g4 * n2 -
                            6720 * a2 * b3 * d * e4 * f * g4 * n2 -
                            420 * b5 * d * e4 * f * g4 * n2 - 3570 * a * b4 * e4 * f2 * g4 * n2 -
                            168 * b5 * d2 * e3 * g5 * n2 - 1120 * b4 * c * e3 * f * g5 * n2 -
                            3360 * a * b4 * d * e3 * f * g5 * n2 - 1176 * b5 * e3 * f2 * g5 * n2 -
                            672 * b5 * d * e2 * f * g6 * n2 + 35 * a3 * e8 * n4 +
                            21 * a5 * e8 * n4 + 120 * a2 * b * e7 * g * n4 +
                            120 * a4 * b * e7 * g * n4 + 140 * a * b2 * e6 * g2 * n4 +
                            280 * a3 * b2 * e6 * g2 * n4 + 56 * b3 * e5 * g3 * n4 +
                            336 * a2 * b3 * e5 * g3 * n4 + 210 * a * b4 * e4 * g4 * n4 +
                            56 * b5 * e3 * g5 * n4)) /
                        3 +
                    (70 * a4 * c * e9 * f2 + 126 * a5 * d * e9 * f2 + 70 * a4 * b * e9 * f3 +
                        630 * a3 * c2 * e8 * f * g + 2520 * a4 * c * d * e8 * f * g +
                        2268 * a5 * d2 * e8 * f * g + 1575 * a3 * b * c * e8 * f2 * g +
                        3150 * a4 * b * d * e8 * f2 * g + 945 * a3 * b2 * e8 * f3 * g +
                        2160 * a2 * b * c2 * e7 * f * g2 + 10080 * a3 * b * c * d * e7 * f * g2 +
                        10080 * a4 * b * d2 * e7 * f * g2 + 4860 * a2 * b2 * c * e7 * f2 * g2 +
                        11340 * a3 * b2 * d * e7 * f2 * g2 + 2700 * a2 * b3 * e7 * f3 * g2 +
                        2520 * a * b2 * c2 * e6 * f * g3 + 15120 * a2 * b2 * c * d * e6 * f * g3 +
                        17640 * a3 * b2 * d2 * e6 * f * g3 + 5460 * a * b3 * c * e6 * f2 * g3 +
                        16380 * a2 * b3 * d * e6 * f2 * g3 + 2940 * a * b4 * e6 * f3 * g3 +
                        1008 * b3 * c2 * e5 * f * g4 + 10080 * a * b3 * c * d * e5 * f * g4 +
                        15120 * a2 * b3 * d2 * e5 * f * g4 + 2142 * b4 * c * e5 * f2 * g4 +
                        10710 * a * b4 * d * e5 * f2 * g4 + 1134 * b5 * e5 * f3 * g4 +
                        2520 * b4 * c * d * e4 * f * g5 + 6300 * a * b4 * d2 * e4 * f * g5 +
                        2646 * b5 * d * e4 * f2 * g5 + 1008 * b5 * d2 * e3 * f * g6 -
                        70 * a4 * c * e9 * n2 - 126 * a5 * d * e9 * n2 - 70 * a4 * b * e9 * f * n2 -
                        315 * a3 * b * c * e8 * g * n2 - 630 * a4 * b * d * e8 * g * n2 -
                        210 * a3 * e8 * f * g * n2 - 126 * a5 * e8 * f * g * n2 -
                        315 * a3 * b2 * e8 * f * g * n2 - 540 * a2 * b2 * c * e7 * g2 * n2 -
                        1260 * a3 * b2 * d * e7 * g2 * n2 - 720 * a2 * b * e7 * f * g2 * n2 -
                        720 * a4 * b * e7 * f * g2 * n2 - 540 * a2 * b3 * e7 * f * g2 * n2 -
                        420 * a * b3 * c * e6 * g3 * n2 - 1260 * a2 * b3 * d * e6 * g3 * n2 -
                        840 * a * b2 * e6 * f * g3 * n2 - 1680 * a3 * b2 * e6 * f * g3 * n2 -
                        420 * a * b4 * e6 * f * g3 * n2 - 126 * b4 * c * e5 * g4 * n2 -
                        630 * a * b4 * d * e5 * g4 * n2 - 336 * b3 * e5 * f * g4 * n2 -
                        2016 * a2 * b3 * e5 * f * g4 * n2 - 126 * b5 * e5 * f * g4 * n2 -
                        126 * b5 * d * e4 * g5 * n2 - 1260 * a * b4 * e4 * f * g5 * n2 -
                        336 * b5 * e3 * f * g6 * n2) +
                    ((126 * a5 * e10 * f2 + 2800 * a4 * c * e9 * f * g +
                        5040 * a5 * d * e9 * f * g + 3500 * a4 * b * e9 * f2 * g +
                        12600 * a3 * b * c * e8 * f * g2 + 25200 * a4 * b * d * e8 * f * g2 +
                        14175 * a3 * b2 * e8 * f2 * g2 + 21600 * a2 * b2 * c * e7 * f * g3 +
                        50400 * a3 * b2 * d * e7 * f * g3 + 23400 * a2 * b3 * e7 * f2 * g3 +
                        16800 * a * b3 * c * e6 * f * g4 + 50400 * a2 * b3 * d * e6 * f * g4 +
                        17850 * a * b4 * e6 * f2 * g4 + 5040 * b4 * c * e5 * f * g5 +
                        25200 * a * b4 * d * e5 * f * g5 + 5292 * b5 * e5 * f2 * g5 +
                        5040 * b5 * d * e4 * f * g6 - 126 * a5 * e10 * n2 -
                        700 * a4 * b * e9 * g * n2 - 1575 * a3 * b2 * e8 * g2 * n2 -
                        1800 * a2 * b3 * e7 * g3 * n2 - 1050 * a * b4 * e6 * g4 * n2 -
                        252 * b5 * e5 * g5 * n2)) /
                        11 +
                    (e5 * f * g *
                        (126 * a5 * e5 + 700 * a4 * b * e4 * g + 1575 * a3 * b2 * e3 * g2 +
                            1800 * a2 * b3 * e2 * g3 + 1050 * a * b4 * e * g4 + 252 * b5 * g5)) /
                        3)) /
                (6300 * n15) -
            (g2 * (-3 * f2 + n2) *
                (((252 * c5 * d5 + 1050 * a * c4 * d6 + 1800 * a2 * c3 * d7 + 1575 * a3 * c2 * d8 +
                     700 * a4 * c * d9 + 126 * a5 * d10 + 1260 * b * c4 * d5 * f +
                     4200 * a * b * c3 * d6 * f + 5400 * a2 * b * c2 * d7 * f +
                     3150 * a3 * b * c * d8 * f + 700 * a4 * b * d9 * f + 2520 * b2 * c3 * d5 * f2 +
                     6300 * a * b2 * c2 * d6 * f2 + 5400 * a2 * b2 * c * d7 * f2 +
                     1575 * a3 * b2 * d8 * f2 + 2520 * b3 * c2 * d5 * f3 +
                     4200 * a * b3 * c * d6 * f3 + 1800 * a2 * b3 * d7 * f3 +
                     1260 * b4 * c * d5 * f4 + 1050 * a * b4 * d6 * f4 + 252 * b5 * d5 * f5 -
                     840 * c5 * d3 * n2 - 3150 * a * c4 * d4 * n2 - 840 * c3 * d5 * n2 -
                     5040 * a2 * c3 * d5 * n2 - 2100 * a * c2 * d6 * n2 - 4200 * a3 * c2 * d6 * n2 -
                     1800 * a2 * c * d7 * n2 - 1800 * a4 * c * d7 * n2 - 525 * a3 * d8 * n2 -
                     315 * a5 * d8 * n2 - 4200 * b * c4 * d3 * f * n2 -
                     12600 * a * b * c3 * d4 * f * n2 - 2520 * b * c2 * d5 * f * n2 -
                     15120 * a2 * b * c2 * d5 * f * n2 - 4200 * a * b * c * d6 * f * n2 -
                     8400 * a3 * b * c * d6 * f * n2 - 1800 * a2 * b * d7 * f * n2 -
                     1800 * a4 * b * d7 * f * n2 - 8400 * b2 * c3 * d3 * f2 * n2 -
                     18900 * a * b2 * c2 * d4 * f2 * n2 - 2520 * b2 * c * d5 * f2 * n2 -
                     15120 * a2 * b2 * c * d5 * f2 * n2 - 2100 * a * b2 * d6 * f2 * n2 -
                     4200 * a3 * b2 * d6 * f2 * n2 - 8400 * b3 * c2 * d3 * f3 * n2 -
                     12600 * a * b3 * c * d4 * f3 * n2 - 840 * b3 * d5 * f3 * n2 -
                     5040 * a2 * b3 * d5 * f3 * n2 - 4200 * b4 * c * d3 * f4 * n2 -
                     3150 * a * b4 * d4 * f4 * n2 - 840 * b5 * d3 * f5 * n2 + 1260 * c5 * d * n4 +
                     3150 * a * c4 * d2 * n4 + 2800 * c3 * d3 * n4 + 4200 * a2 * c3 * d3 * n4 +
                     6300 * a * c2 * d4 * n4 + 3150 * a3 * c2 * d4 * n4 + 1260 * c * d5 * n4 +
                     5040 * a2 * c * d5 * n4 + 1260 * a4 * c * d5 * n4 + 1050 * a * d6 * n4 +
                     1400 * a3 * d6 * n4 + 210 * a5 * d6 * n4 + 6300 * b * c4 * d * f * n4 +
                     12600 * a * b * c3 * d2 * f * n4 + 8400 * b * c2 * d3 * f * n4 +
                     12600 * a2 * b * c2 * d3 * f * n4 + 12600 * a * b * c * d4 * f * n4 +
                     6300 * a3 * b * c * d4 * f * n4 + 1260 * b * d5 * f * n4 +
                     5040 * a2 * b * d5 * f * n4 + 1260 * a4 * b * d5 * f * n4 +
                     12600 * b2 * c3 * d * f2 * n4 + 18900 * a * b2 * c2 * d2 * f2 * n4 +
                     8400 * b2 * c * d3 * f2 * n4 + 12600 * a2 * b2 * c * d3 * f2 * n4 +
                     6300 * a * b2 * d4 * f2 * n4 + 3150 * a3 * b2 * d4 * f2 * n4 +
                     12600 * b3 * c2 * d * f3 * n4 + 12600 * a * b3 * c * d2 * f3 * n4 +
                     2800 * b3 * d3 * f3 * n4 + 4200 * a2 * b3 * d3 * f3 * n4 +
                     6300 * b4 * c * d * f4 * n4 + 3150 * a * b4 * d2 * f4 * n4 +
                     1260 * b5 * d * f5 * n4 - 4200 * c3 * d * n6 - 6300 * a * c2 * d2 * n6 -
                     4200 * c * d3 * n6 - 4200 * a2 * c * d3 * n6 - 3150 * a * d4 * n6 -
                     1050 * a3 * d4 * n6 - 12600 * b * c2 * d * f * n6 -
                     12600 * a * b * c * d2 * f * n6 - 4200 * b * d3 * f * n6 -
                     4200 * a2 * b * d3 * f * n6 - 12600 * b2 * c * d * f2 * n6 -
                     6300 * a * b2 * d2 * f2 * n6 - 4200 * b3 * d * f3 * n6 + 6300 * c * d * n8 +
                     3150 * a * d2 * n8 + 6300 * b * d * f * n8)) /
                        3 +
                    (5 *
                        (126 * c5 * d4 * e + 630 * a * c4 * d5 * e + 1260 * a2 * c3 * d6 * e +
                            1260 * a3 * c2 * d7 * e + 630 * a4 * c * d8 * e + 126 * a5 * d9 * e +
                            630 * b * c4 * d4 * e * f + 2520 * a * b * c3 * d5 * e * f +
                            3780 * a2 * b * c2 * d6 * e * f + 2520 * a3 * b * c * d7 * e * f +
                            630 * a4 * b * d8 * e * f + 1260 * b2 * c3 * d4 * e * f2 +
                            3780 * a * b2 * c2 * d5 * e * f2 + 3780 * a2 * b2 * c * d6 * e * f2 +
                            1260 * a3 * b2 * d7 * e * f2 + 1260 * b3 * c2 * d4 * e * f3 +
                            2520 * a * b3 * c * d5 * e * f3 + 1260 * a2 * b3 * d6 * e * f3 +
                            630 * b4 * c * d4 * e * f4 + 630 * a * b4 * d5 * e * f4 +
                            126 * b5 * d4 * e * f5 + 126 * b * c4 * d5 * g +
                            420 * a * b * c3 * d6 * g + 540 * a2 * b * c2 * d7 * g +
                            315 * a3 * b * c * d8 * g + 70 * a4 * b * d9 * g +
                            504 * b2 * c3 * d5 * f * g + 1260 * a * b2 * c2 * d6 * f * g +
                            1080 * a2 * b2 * c * d7 * f * g + 315 * a3 * b2 * d8 * f * g +
                            756 * b3 * c2 * d5 * f2 * g + 1260 * a * b3 * c * d6 * f2 * g +
                            540 * a2 * b3 * d7 * f2 * g + 504 * b4 * c * d5 * f3 * g +
                            420 * a * b4 * d6 * f3 * g + 126 * b5 * d5 * f4 * g -
                            252 * c5 * d2 * e * n2 - 1260 * a * c4 * d3 * e * n2 -
                            420 * c3 * d4 * e * n2 - 2520 * a2 * c3 * d4 * e * n2 -
                            1260 * a * c2 * d5 * e * n2 - 2520 * a3 * c2 * d5 * e * n2 -
                            1260 * a2 * c * d6 * e * n2 - 1260 * a4 * c * d6 * e * n2 -
                            420 * a3 * d7 * e * n2 - 252 * a5 * d7 * e * n2 -
                            1260 * b * c4 * d2 * e * f * n2 - 5040 * a * b * c3 * d3 * e * f * n2 -
                            1260 * b * c2 * d4 * e * f * n2 - 7560 * a2 * b * c2 * d4 * e * f * n2 -
                            2520 * a * b * c * d5 * e * f * n2 -
                            5040 * a3 * b * c * d5 * e * f * n2 - 1260 * a2 * b * d6 * e * f * n2 -
                            1260 * a4 * b * d6 * e * f * n2 - 2520 * b2 * c3 * d2 * e * f2 * n2 -
                            7560 * a * b2 * c2 * d3 * e * f2 * n2 -
                            1260 * b2 * c * d4 * e * f2 * n2 -
                            7560 * a2 * b2 * c * d4 * e * f2 * n2 -
                            1260 * a * b2 * d5 * e * f2 * n2 - 2520 * a3 * b2 * d5 * e * f2 * n2 -
                            2520 * b3 * c2 * d2 * e * f3 * n2 -
                            5040 * a * b3 * c * d3 * e * f3 * n2 - 420 * b3 * d4 * e * f3 * n2 -
                            2520 * a2 * b3 * d4 * e * f3 * n2 - 1260 * b4 * c * d2 * e * f4 * n2 -
                            1260 * a * b4 * d3 * e * f4 * n2 - 252 * b5 * d2 * e * f5 * n2 -
                            420 * b * c4 * d3 * g * n2 - 1260 * a * b * c3 * d4 * g * n2 -
                            252 * b * c2 * d5 * g * n2 - 1512 * a2 * b * c2 * d5 * g * n2 -
                            420 * a * b * c * d6 * g * n2 - 840 * a3 * b * c * d6 * g * n2 -
                            180 * a2 * b * d7 * g * n2 - 180 * a4 * b * d7 * g * n2 -
                            1680 * b2 * c3 * d3 * f * g * n2 -
                            3780 * a * b2 * c2 * d4 * f * g * n2 - 504 * b2 * c * d5 * f * g * n2 -
                            3024 * a2 * b2 * c * d5 * f * g * n2 - 420 * a * b2 * d6 * f * g * n2 -
                            840 * a3 * b2 * d6 * f * g * n2 - 2520 * b3 * c2 * d3 * f2 * g * n2 -
                            3780 * a * b3 * c * d4 * f2 * g * n2 - 252 * b3 * d5 * f2 * g * n2 -
                            1512 * a2 * b3 * d5 * f2 * g * n2 - 1680 * b4 * c * d3 * f3 * g * n2 -
                            1260 * a * b4 * d4 * f3 * g * n2 - 420 * b5 * d3 * f4 * g * n2 +
                            126 * c5 * e * n4 + 630 * a * c4 * d * e * n4 + 840 * c3 * d2 * e * n4 +
                            1260 * a2 * c3 * d2 * e * n4 + 2520 * a * c2 * d3 * e * n4 +
                            1260 * a3 * c2 * d3 * e * n4 + 630 * c * d4 * e * n4 +
                            2520 * a2 * c * d4 * e * n4 + 630 * a4 * c * d4 * e * n4 +
                            630 * a * d5 * e * n4 + 840 * a3 * d5 * e * n4 +
                            126 * a5 * d5 * e * n4 + 630 * b * c4 * e * f * n4 +
                            2520 * a * b * c3 * d * e * f * n4 + 2520 * b * c2 * d2 * e * f * n4 +
                            3780 * a2 * b * c2 * d2 * e * f * n4 +
                            5040 * a * b * c * d3 * e * f * n4 +
                            2520 * a3 * b * c * d3 * e * f * n4 + 630 * b * d4 * e * f * n4 +
                            2520 * a2 * b * d4 * e * f * n4 + 630 * a4 * b * d4 * e * f * n4 +
                            1260 * b2 * c3 * e * f2 * n4 + 3780 * a * b2 * c2 * d * e * f2 * n4 +
                            2520 * b2 * c * d2 * e * f2 * n4 +
                            3780 * a2 * b2 * c * d2 * e * f2 * n4 +
                            2520 * a * b2 * d3 * e * f2 * n4 + 1260 * a3 * b2 * d3 * e * f2 * n4 +
                            1260 * b3 * c2 * e * f3 * n4 + 2520 * a * b3 * c * d * e * f3 * n4 +
                            840 * b3 * d2 * e * f3 * n4 + 1260 * a2 * b3 * d2 * e * f3 * n4 +
                            630 * b4 * c * e * f4 * n4 + 630 * a * b4 * d * e * f4 * n4 +
                            126 * b5 * e * f5 * n4 + 630 * b * c4 * d * g * n4 +
                            1260 * a * b * c3 * d2 * g * n4 + 840 * b * c2 * d3 * g * n4 +
                            1260 * a2 * b * c2 * d3 * g * n4 + 1260 * a * b * c * d4 * g * n4 +
                            630 * a3 * b * c * d4 * g * n4 + 126 * b * d5 * g * n4 +
                            504 * a2 * b * d5 * g * n4 + 126 * a4 * b * d5 * g * n4 +
                            2520 * b2 * c3 * d * f * g * n4 + 3780 * a * b2 * c2 * d2 * f * g * n4 +
                            1680 * b2 * c * d3 * f * g * n4 + 2520 * a2 * b2 * c * d3 * f * g * n4 +
                            1260 * a * b2 * d4 * f * g * n4 + 630 * a3 * b2 * d4 * f * g * n4 +
                            3780 * b3 * c2 * d * f2 * g * n4 +
                            3780 * a * b3 * c * d2 * f2 * g * n4 + 840 * b3 * d3 * f2 * g * n4 +
                            1260 * a2 * b3 * d3 * f2 * g * n4 + 2520 * b4 * c * d * f3 * g * n4 +
                            1260 * a * b4 * d2 * f3 * g * n4 + 630 * b5 * d * f4 * g * n4 -
                            420 * c3 * e * n6 - 1260 * a * c2 * d * e * n6 -
                            1260 * c * d2 * e * n6 - 1260 * a2 * c * d2 * e * n6 -
                            1260 * a * d3 * e * n6 - 420 * a3 * d3 * e * n6 -
                            1260 * b * c2 * e * f * n6 - 2520 * a * b * c * d * e * f * n6 -
                            1260 * b * d2 * e * f * n6 - 1260 * a2 * b * d2 * e * f * n6 -
                            1260 * b2 * c * e * f2 * n6 - 1260 * a * b2 * d * e * f2 * n6 -
                            420 * b3 * e * f3 * n6 - 1260 * b * c2 * d * g * n6 -
                            1260 * a * b * c * d2 * g * n6 - 420 * b * d3 * g * n6 -
                            420 * a2 * b * d3 * g * n6 - 2520 * b2 * c * d * f * g * n6 -
                            1260 * a * b2 * d2 * f * g * n6 - 1260 * b3 * d * f2 * g * n6 +
                            630 * c * e * n8 + 630 * a * d * e * n8 + 630 * b * e * f * n8 +
                            630 * b * d * g * n8)) /
                        2 +
                    3 * (168 * c5 * d3 * e2 + 1050 * a * c4 * d4 * e2 + 2520 * a2 * c3 * d5 * e2 +
                            2940 * a3 * c2 * d6 * e2 + 1680 * a4 * c * d7 * e2 +
                            378 * a5 * d8 * e2 + 840 * b * c4 * d3 * e2 * f +
                            4200 * a * b * c3 * d4 * e2 * f + 7560 * a2 * b * c2 * d5 * e2 * f +
                            5880 * a3 * b * c * d6 * e2 * f + 1680 * a4 * b * d7 * e2 * f +
                            1680 * b2 * c3 * d3 * e2 * f2 + 6300 * a * b2 * c2 * d4 * e2 * f2 +
                            7560 * a2 * b2 * c * d5 * e2 * f2 + 2940 * a3 * b2 * d6 * e2 * f2 +
                            1680 * b3 * c2 * d3 * e2 * f3 + 4200 * a * b3 * c * d4 * e2 * f3 +
                            2520 * a2 * b3 * d5 * e2 * f3 + 840 * b4 * c * d3 * e2 * f4 +
                            1050 * a * b4 * d4 * e2 * f4 + 168 * b5 * d3 * e2 * f5 +
                            420 * b * c4 * d4 * e * g + 1680 * a * b * c3 * d5 * e * g +
                            2520 * a2 * b * c2 * d6 * e * g + 1680 * a3 * b * c * d7 * e * g +
                            420 * a4 * b * d8 * e * g + 1680 * b2 * c3 * d4 * e * f * g +
                            5040 * a * b2 * c2 * d5 * e * f * g +
                            5040 * a2 * b2 * c * d6 * e * f * g + 1680 * a3 * b2 * d7 * e * f * g +
                            2520 * b3 * c2 * d4 * e * f2 * g + 5040 * a * b3 * c * d5 * e * f2 * g +
                            2520 * a2 * b3 * d6 * e * f2 * g + 1680 * b4 * c * d4 * e * f3 * g +
                            1680 * a * b4 * d5 * e * f3 * g + 420 * b5 * d4 * e * f4 * g +
                            168 * b2 * c3 * d5 * g2 + 420 * a * b2 * c2 * d6 * g2 +
                            360 * a2 * b2 * c * d7 * g2 + 105 * a3 * b2 * d8 * g2 +
                            504 * b3 * c2 * d5 * f * g2 + 840 * a * b3 * c * d6 * f * g2 +
                            360 * a2 * b3 * d7 * f * g2 + 504 * b4 * c * d5 * f2 * g2 +
                            420 * a * b4 * d6 * f2 * g2 + 168 * b5 * d5 * f3 * g2 -
                            168 * c5 * d * e2 * n2 - 1260 * a * c4 * d2 * e2 * n2 -
                            560 * c3 * d3 * e2 * n2 - 3360 * a2 * c3 * d3 * e2 * n2 -
                            2100 * a * c2 * d4 * e2 * n2 - 4200 * a3 * c2 * d4 * e2 * n2 -
                            2520 * a2 * c * d5 * e2 * n2 - 2520 * a4 * c * d5 * e2 * n2 -
                            980 * a3 * d6 * e2 * n2 - 588 * a5 * d6 * e2 * n2 -
                            840 * b * c4 * d * e2 * f * n2 - 5040 * a * b * c3 * d2 * e2 * f * n2 -
                            1680 * b * c2 * d3 * e2 * f * n2 -
                            10080 * a2 * b * c2 * d3 * e2 * f * n2 -
                            4200 * a * b * c * d4 * e2 * f * n2 -
                            8400 * a3 * b * c * d4 * e2 * f * n2 -
                            2520 * a2 * b * d5 * e2 * f * n2 - 2520 * a4 * b * d5 * e2 * f * n2 -
                            1680 * b2 * c3 * d * e2 * f2 * n2 -
                            7560 * a * b2 * c2 * d2 * e2 * f2 * n2 -
                            1680 * b2 * c * d3 * e2 * f2 * n2 -
                            10080 * a2 * b2 * c * d3 * e2 * f2 * n2 -
                            2100 * a * b2 * d4 * e2 * f2 * n2 - 4200 * a3 * b2 * d4 * e2 * f2 * n2 -
                            1680 * b3 * c2 * d * e2 * f3 * n2 -
                            5040 * a * b3 * c * d2 * e2 * f3 * n2 - 560 * b3 * d3 * e2 * f3 * n2 -
                            3360 * a2 * b3 * d3 * e2 * f3 * n2 - 840 * b4 * c * d * e2 * f4 * n2 -
                            1260 * a * b4 * d2 * e2 * f4 * n2 - 168 * b5 * d * e2 * f5 * n2 -
                            840 * b * c4 * d2 * e * g * n2 - 3360 * a * b * c3 * d3 * e * g * n2 -
                            840 * b * c2 * d4 * e * g * n2 - 5040 * a2 * b * c2 * d4 * e * g * n2 -
                            1680 * a * b * c * d5 * e * g * n2 -
                            3360 * a3 * b * c * d5 * e * g * n2 - 840 * a2 * b * d6 * e * g * n2 -
                            840 * a4 * b * d6 * e * g * n2 - 3360 * b2 * c3 * d2 * e * f * g * n2 -
                            10080 * a * b2 * c2 * d3 * e * f * g * n2 -
                            1680 * b2 * c * d4 * e * f * g * n2 -
                            10080 * a2 * b2 * c * d4 * e * f * g * n2 -
                            1680 * a * b2 * d5 * e * f * g * n2 -
                            3360 * a3 * b2 * d5 * e * f * g * n2 -
                            5040 * b3 * c2 * d2 * e * f2 * g * n2 -
                            10080 * a * b3 * c * d3 * e * f2 * g * n2 -
                            840 * b3 * d4 * e * f2 * g * n2 -
                            5040 * a2 * b3 * d4 * e * f2 * g * n2 -
                            3360 * b4 * c * d2 * e * f3 * g * n2 -
                            3360 * a * b4 * d3 * e * f3 * g * n2 - 840 * b5 * d2 * e * f4 * g * n2 -
                            560 * b2 * c3 * d3 * g2 * n2 - 1260 * a * b2 * c2 * d4 * g2 * n2 -
                            168 * b2 * c * d5 * g2 * n2 - 1008 * a2 * b2 * c * d5 * g2 * n2 -
                            140 * a * b2 * d6 * g2 * n2 - 280 * a3 * b2 * d6 * g2 * n2 -
                            1680 * b3 * c2 * d3 * f * g2 * n2 -
                            2520 * a * b3 * c * d4 * f * g2 * n2 - 168 * b3 * d5 * f * g2 * n2 -
                            1008 * a2 * b3 * d5 * f * g2 * n2 - 1680 * b4 * c * d3 * f2 * g2 * n2 -
                            1260 * a * b4 * d4 * f2 * g2 * n2 - 560 * b5 * d3 * f3 * g2 * n2 +
                            210 * a * c4 * e2 * n4 + 560 * c3 * d * e2 * n4 +
                            840 * a2 * c3 * d * e2 * n4 + 2520 * a * c2 * d2 * e2 * n4 +
                            1260 * a3 * c2 * d2 * e2 * n4 + 840 * c * d3 * e2 * n4 +
                            3360 * a2 * c * d3 * e2 * n4 + 840 * a4 * c * d3 * e2 * n4 +
                            1050 * a * d4 * e2 * n4 + 1400 * a3 * d4 * e2 * n4 +
                            210 * a5 * d4 * e2 * n4 + 840 * a * b * c3 * e2 * f * n4 +
                            1680 * b * c2 * d * e2 * f * n4 + 2520 * a2 * b * c2 * d * e2 * f * n4 +
                            5040 * a * b * c * d2 * e2 * f * n4 +
                            2520 * a3 * b * c * d2 * e2 * f * n4 + 840 * b * d3 * e2 * f * n4 +
                            3360 * a2 * b * d3 * e2 * f * n4 + 840 * a4 * b * d3 * e2 * f * n4 +
                            1260 * a * b2 * c2 * e2 * f2 * n4 + 1680 * b2 * c * d * e2 * f2 * n4 +
                            2520 * a2 * b2 * c * d * e2 * f2 * n4 +
                            2520 * a * b2 * d2 * e2 * f2 * n4 + 1260 * a3 * b2 * d2 * e2 * f2 * n4 +
                            840 * a * b3 * c * e2 * f3 * n4 + 560 * b3 * d * e2 * f3 * n4 +
                            840 * a2 * b3 * d * e2 * f3 * n4 + 210 * a * b4 * e2 * f4 * n4 +
                            420 * b * c4 * e * g * n4 + 1680 * a * b * c3 * d * e * g * n4 +
                            1680 * b * c2 * d2 * e * g * n4 + 2520 * a2 * b * c2 * d2 * e * g * n4 +
                            3360 * a * b * c * d3 * e * g * n4 +
                            1680 * a3 * b * c * d3 * e * g * n4 + 420 * b * d4 * e * g * n4 +
                            1680 * a2 * b * d4 * e * g * n4 + 420 * a4 * b * d4 * e * g * n4 +
                            1680 * b2 * c3 * e * f * g * n4 +
                            5040 * a * b2 * c2 * d * e * f * g * n4 +
                            3360 * b2 * c * d2 * e * f * g * n4 +
                            5040 * a2 * b2 * c * d2 * e * f * g * n4 +
                            3360 * a * b2 * d3 * e * f * g * n4 +
                            1680 * a3 * b2 * d3 * e * f * g * n4 +
                            2520 * b3 * c2 * e * f2 * g * n4 +
                            5040 * a * b3 * c * d * e * f2 * g * n4 +
                            1680 * b3 * d2 * e * f2 * g * n4 +
                            2520 * a2 * b3 * d2 * e * f2 * g * n4 +
                            1680 * b4 * c * e * f3 * g * n4 + 1680 * a * b4 * d * e * f3 * g * n4 +
                            420 * b5 * e * f4 * g * n4 + 840 * b2 * c3 * d * g2 * n4 +
                            1260 * a * b2 * c2 * d2 * g2 * n4 + 560 * b2 * c * d3 * g2 * n4 +
                            840 * a2 * b2 * c * d3 * g2 * n4 + 420 * a * b2 * d4 * g2 * n4 +
                            210 * a3 * b2 * d4 * g2 * n4 + 2520 * b3 * c2 * d * f * g2 * n4 +
                            2520 * a * b3 * c * d2 * f * g2 * n4 + 560 * b3 * d3 * f * g2 * n4 +
                            840 * a2 * b3 * d3 * f * g2 * n4 + 2520 * b4 * c * d * f2 * g2 * n4 +
                            1260 * a * b4 * d2 * f2 * g2 * n4 + 840 * b5 * d * f3 * g2 * n4 -
                            420 * a * c2 * e2 * n6 - 840 * c * d * e2 * n6 -
                            840 * a2 * c * d * e2 * n6 - 1260 * a * d2 * e2 * n6 -
                            420 * a3 * d2 * e2 * n6 - 840 * a * b * c * e2 * f * n6 -
                            840 * b * d * e2 * f * n6 - 840 * a2 * b * d * e2 * f * n6 -
                            420 * a * b2 * e2 * f2 * n6 - 840 * b * c2 * e * g * n6 -
                            1680 * a * b * c * d * e * g * n6 - 840 * b * d2 * e * g * n6 -
                            840 * a2 * b * d2 * e * g * n6 - 1680 * b2 * c * e * f * g * n6 -
                            1680 * a * b2 * d * e * f * g * n6 - 840 * b3 * e * f2 * g * n6 -
                            840 * b2 * c * d * g2 * n6 - 420 * a * b2 * d2 * g2 * n6 -
                            840 * b3 * d * f * g2 * n6 + 210 * a * e2 * n8 + 420 * b * e * g * n8) +
                    (20 *
                        (63 * c5 * d2 * e3 + 525 * a * c4 * d3 * e3 + 1575 * a2 * c3 * d4 * e3 +
                            2205 * a3 * c2 * d5 * e3 + 1470 * a4 * c * d6 * e3 +
                            378 * a5 * d7 * e3 + 315 * b * c4 * d2 * e3 * f +
                            2100 * a * b * c3 * d3 * e3 * f + 4725 * a2 * b * c2 * d4 * e3 * f +
                            4410 * a3 * b * c * d5 * e3 * f + 1470 * a4 * b * d6 * e3 * f +
                            630 * b2 * c3 * d2 * e3 * f2 + 3150 * a * b2 * c2 * d3 * e3 * f2 +
                            4725 * a2 * b2 * c * d4 * e3 * f2 + 2205 * a3 * b2 * d5 * e3 * f2 +
                            630 * b3 * c2 * d2 * e3 * f3 + 2100 * a * b3 * c * d3 * e3 * f3 +
                            1575 * a2 * b3 * d4 * e3 * f3 + 315 * b4 * c * d2 * e3 * f4 +
                            525 * a * b4 * d3 * e3 * f4 + 63 * b5 * d2 * e3 * f5 +
                            315 * b * c4 * d3 * e2 * g + 1575 * a * b * c3 * d4 * e2 * g +
                            2835 * a2 * b * c2 * d5 * e2 * g + 2205 * a3 * b * c * d6 * e2 * g +
                            630 * a4 * b * d7 * e2 * g + 1260 * b2 * c3 * d3 * e2 * f * g +
                            4725 * a * b2 * c2 * d4 * e2 * f * g +
                            5670 * a2 * b2 * c * d5 * e2 * f * g +
                            2205 * a3 * b2 * d6 * e2 * f * g + 1890 * b3 * c2 * d3 * e2 * f2 * g +
                            4725 * a * b3 * c * d4 * e2 * f2 * g +
                            2835 * a2 * b3 * d5 * e2 * f2 * g + 1260 * b4 * c * d3 * e2 * f3 * g +
                            1575 * a * b4 * d4 * e2 * f3 * g + 315 * b5 * d3 * e2 * f4 * g +
                            315 * b2 * c3 * d4 * e * g2 + 945 * a * b2 * c2 * d5 * e * g2 +
                            945 * a2 * b2 * c * d6 * e * g2 + 315 * a3 * b2 * d7 * e * g2 +
                            945 * b3 * c2 * d4 * e * f * g2 + 1890 * a * b3 * c * d5 * e * f * g2 +
                            945 * a2 * b3 * d6 * e * f * g2 + 945 * b4 * c * d4 * e * f2 * g2 +
                            945 * a * b4 * d5 * e * f2 * g2 + 315 * b5 * d4 * e * f3 * g2 +
                            63 * b3 * c2 * d5 * g3 + 105 * a * b3 * c * d6 * g3 +
                            45 * a2 * b3 * d7 * g3 + 126 * b4 * c * d5 * f * g3 +
                            105 * a * b4 * d6 * f * g3 + 63 * b5 * d5 * f2 * g3 -
                            21 * c5 * e3 * n2 - 315 * a * c4 * d * e3 * n2 -
                            210 * c3 * d2 * e3 * n2 - 1260 * a2 * c3 * d2 * e3 * n2 -
                            1050 * a * c2 * d3 * e3 * n2 - 2100 * a3 * c2 * d3 * e3 * n2 -
                            1575 * a2 * c * d4 * e3 * n2 - 1575 * a4 * c * d4 * e3 * n2 -
                            735 * a3 * d5 * e3 * n2 - 441 * a5 * d5 * e3 * n2 -
                            105 * b * c4 * e3 * f * n2 - 1260 * a * b * c3 * d * e3 * f * n2 -
                            630 * b * c2 * d2 * e3 * f * n2 -
                            3780 * a2 * b * c2 * d2 * e3 * f * n2 -
                            2100 * a * b * c * d3 * e3 * f * n2 -
                            4200 * a3 * b * c * d3 * e3 * f * n2 -
                            1575 * a2 * b * d4 * e3 * f * n2 - 1575 * a4 * b * d4 * e3 * f * n2 -
                            210 * b2 * c3 * e3 * f2 * n2 - 1890 * a * b2 * c2 * d * e3 * f2 * n2 -
                            630 * b2 * c * d2 * e3 * f2 * n2 -
                            3780 * a2 * b2 * c * d2 * e3 * f2 * n2 -
                            1050 * a * b2 * d3 * e3 * f2 * n2 - 2100 * a3 * b2 * d3 * e3 * f2 * n2 -
                            210 * b3 * c2 * e3 * f3 * n2 - 1260 * a * b3 * c * d * e3 * f3 * n2 -
                            210 * b3 * d2 * e3 * f3 * n2 - 1260 * a2 * b3 * d2 * e3 * f3 * n2 -
                            105 * b4 * c * e3 * f4 * n2 - 315 * a * b4 * d * e3 * f4 * n2 -
                            21 * b5 * e3 * f5 * n2 - 315 * b * c4 * d * e2 * g * n2 -
                            1890 * a * b * c3 * d2 * e2 * g * n2 - 630 * b * c2 * d3 * e2 * g * n2 -
                            3780 * a2 * b * c2 * d3 * e2 * g * n2 -
                            1575 * a * b * c * d4 * e2 * g * n2 -
                            3150 * a3 * b * c * d4 * e2 * g * n2 - 945 * a2 * b * d5 * e2 * g * n2 -
                            945 * a4 * b * d5 * e2 * g * n2 - 1260 * b2 * c3 * d * e2 * f * g * n2 -
                            5670 * a * b2 * c2 * d2 * e2 * f * g * n2 -
                            1260 * b2 * c * d3 * e2 * f * g * n2 -
                            7560 * a2 * b2 * c * d3 * e2 * f * g * n2 -
                            1575 * a * b2 * d4 * e2 * f * g * n2 -
                            3150 * a3 * b2 * d4 * e2 * f * g * n2 -
                            1890 * b3 * c2 * d * e2 * f2 * g * n2 -
                            5670 * a * b3 * c * d2 * e2 * f2 * g * n2 -
                            630 * b3 * d3 * e2 * f2 * g * n2 -
                            3780 * a2 * b3 * d3 * e2 * f2 * g * n2 -
                            1260 * b4 * c * d * e2 * f3 * g * n2 -
                            1890 * a * b4 * d2 * e2 * f3 * g * n2 -
                            315 * b5 * d * e2 * f4 * g * n2 - 630 * b2 * c3 * d2 * e * g2 * n2 -
                            1890 * a * b2 * c2 * d3 * e * g2 * n2 -
                            315 * b2 * c * d4 * e * g2 * n2 -
                            1890 * a2 * b2 * c * d4 * e * g2 * n2 -
                            315 * a * b2 * d5 * e * g2 * n2 - 630 * a3 * b2 * d5 * e * g2 * n2 -
                            1890 * b3 * c2 * d2 * e * f * g2 * n2 -
                            3780 * a * b3 * c * d3 * e * f * g2 * n2 -
                            315 * b3 * d4 * e * f * g2 * n2 -
                            1890 * a2 * b3 * d4 * e * f * g2 * n2 -
                            1890 * b4 * c * d2 * e * f2 * g2 * n2 -
                            1890 * a * b4 * d3 * e * f2 * g2 * n2 -
                            630 * b5 * d2 * e * f3 * g2 * n2 - 210 * b3 * c2 * d3 * g3 * n2 -
                            315 * a * b3 * c * d4 * g3 * n2 - 21 * b3 * d5 * g3 * n2 -
                            126 * a2 * b3 * d5 * g3 * n2 - 420 * b4 * c * d3 * f * g3 * n2 -
                            315 * a * b4 * d4 * f * g3 * n2 - 210 * b5 * d3 * f2 * g3 * n2 +
                            70 * c3 * e3 * n4 + 105 * a2 * c3 * e3 * n4 +
                            630 * a * c2 * d * e3 * n4 + 315 * a3 * c2 * d * e3 * n4 +
                            315 * c * d2 * e3 * n4 + 1260 * a2 * c * d2 * e3 * n4 +
                            315 * a4 * c * d2 * e3 * n4 + 525 * a * d3 * e3 * n4 +
                            700 * a3 * d3 * e3 * n4 + 105 * a5 * d3 * e3 * n4 +
                            210 * b * c2 * e3 * f * n4 + 315 * a2 * b * c2 * e3 * f * n4 +
                            1260 * a * b * c * d * e3 * f * n4 +
                            630 * a3 * b * c * d * e3 * f * n4 + 315 * b * d2 * e3 * f * n4 +
                            1260 * a2 * b * d2 * e3 * f * n4 + 315 * a4 * b * d2 * e3 * f * n4 +
                            210 * b2 * c * e3 * f2 * n4 + 315 * a2 * b2 * c * e3 * f2 * n4 +
                            630 * a * b2 * d * e3 * f2 * n4 + 315 * a3 * b2 * d * e3 * f2 * n4 +
                            70 * b3 * e3 * f3 * n4 + 105 * a2 * b3 * e3 * f3 * n4 +
                            315 * a * b * c3 * e2 * g * n4 + 630 * b * c2 * d * e2 * g * n4 +
                            945 * a2 * b * c2 * d * e2 * g * n4 +
                            1890 * a * b * c * d2 * e2 * g * n4 +
                            945 * a3 * b * c * d2 * e2 * g * n4 + 315 * b * d3 * e2 * g * n4 +
                            1260 * a2 * b * d3 * e2 * g * n4 + 315 * a4 * b * d3 * e2 * g * n4 +
                            945 * a * b2 * c2 * e2 * f * g * n4 +
                            1260 * b2 * c * d * e2 * f * g * n4 +
                            1890 * a2 * b2 * c * d * e2 * f * g * n4 +
                            1890 * a * b2 * d2 * e2 * f * g * n4 +
                            945 * a3 * b2 * d2 * e2 * f * g * n4 +
                            945 * a * b3 * c * e2 * f2 * g * n4 + 630 * b3 * d * e2 * f2 * g * n4 +
                            945 * a2 * b3 * d * e2 * f2 * g * n4 + 315 * a * b4 * e2 * f3 * g * n4 +
                            315 * b2 * c3 * e * g2 * n4 + 945 * a * b2 * c2 * d * e * g2 * n4 +
                            630 * b2 * c * d2 * e * g2 * n4 + 945 * a2 * b2 * c * d2 * e * g2 * n4 +
                            630 * a * b2 * d3 * e * g2 * n4 + 315 * a3 * b2 * d3 * e * g2 * n4 +
                            945 * b3 * c2 * e * f * g2 * n4 +
                            1890 * a * b3 * c * d * e * f * g2 * n4 +
                            630 * b3 * d2 * e * f * g2 * n4 + 945 * a2 * b3 * d2 * e * f * g2 * n4 +
                            945 * b4 * c * e * f2 * g2 * n4 + 945 * a * b4 * d * e * f2 * g2 * n4 +
                            315 * b5 * e * f3 * g2 * n4 + 315 * b3 * c2 * d * g3 * n4 +
                            315 * a * b3 * c * d2 * g3 * n4 + 70 * b3 * d3 * g3 * n4 +
                            105 * a2 * b3 * d3 * g3 * n4 + 630 * b4 * c * d * f * g3 * n4 +
                            315 * a * b4 * d2 * f * g3 * n4 + 315 * b5 * d * f2 * g3 * n4 -
                            105 * c * e3 * n6 - 105 * a2 * c * e3 * n6 - 315 * a * d * e3 * n6 -
                            105 * a3 * d * e3 * n6 - 105 * b * e3 * f * n6 -
                            105 * a2 * b * e3 * f * n6 - 315 * a * b * c * e2 * g * n6 -
                            315 * b * d * e2 * g * n6 - 315 * a2 * b * d * e2 * g * n6 -
                            315 * a * b2 * e2 * f * g * n6 - 315 * b2 * c * e * g2 * n6 -
                            315 * a * b2 * d * e * g2 * n6 - 315 * b3 * e * f * g2 * n6 -
                            105 * b3 * d * g3 * n6)) /
                        3 +
                    30 *
                        (6 * c5 * d * e4 + 75 * a * c4 * d2 * e4 + 300 * a2 * c3 * d3 * e4 +
                            525 * a3 * c2 * d4 * e4 + 420 * a4 * c * d5 * e4 + 126 * a5 * d6 * e4 +
                            30 * b * c4 * d * e4 * f + 300 * a * b * c3 * d2 * e4 * f +
                            900 * a2 * b * c2 * d3 * e4 * f + 1050 * a3 * b * c * d4 * e4 * f +
                            420 * a4 * b * d5 * e4 * f + 60 * b2 * c3 * d * e4 * f2 +
                            450 * a * b2 * c2 * d2 * e4 * f2 + 900 * a2 * b2 * c * d3 * e4 * f2 +
                            525 * a3 * b2 * d4 * e4 * f2 + 60 * b3 * c2 * d * e4 * f3 +
                            300 * a * b3 * c * d2 * e4 * f3 + 300 * a2 * b3 * d3 * e4 * f3 +
                            30 * b4 * c * d * e4 * f4 + 75 * a * b4 * d2 * e4 * f4 +
                            6 * b5 * d * e4 * f5 + 60 * b * c4 * d2 * e3 * g +
                            400 * a * b * c3 * d3 * e3 * g + 900 * a2 * b * c2 * d4 * e3 * g +
                            840 * a3 * b * c * d5 * e3 * g + 280 * a4 * b * d6 * e3 * g +
                            240 * b2 * c3 * d2 * e3 * f * g + 1200 * a * b2 * c2 * d3 * e3 * f * g +
                            1800 * a2 * b2 * c * d4 * e3 * f * g + 840 * a3 * b2 * d5 * e3 * f * g +
                            360 * b3 * c2 * d2 * e3 * f2 * g +
                            1200 * a * b3 * c * d3 * e3 * f2 * g +
                            900 * a2 * b3 * d4 * e3 * f2 * g + 240 * b4 * c * d2 * e3 * f3 * g +
                            400 * a * b4 * d3 * e3 * f3 * g + 60 * b5 * d2 * e3 * f4 * g +
                            120 * b2 * c3 * d3 * e2 * g2 + 450 * a * b2 * c2 * d4 * e2 * g2 +
                            540 * a2 * b2 * c * d5 * e2 * g2 + 210 * a3 * b2 * d6 * e2 * g2 +
                            360 * b3 * c2 * d3 * e2 * f * g2 + 900 * a * b3 * c * d4 * e2 * f * g2 +
                            540 * a2 * b3 * d5 * e2 * f * g2 + 360 * b4 * c * d3 * e2 * f2 * g2 +
                            450 * a * b4 * d4 * e2 * f2 * g2 + 120 * b5 * d3 * e2 * f3 * g2 +
                            60 * b3 * c2 * d4 * e * g3 + 120 * a * b3 * c * d5 * e * g3 +
                            60 * a2 * b3 * d6 * e * g3 + 120 * b4 * c * d4 * e * f * g3 +
                            120 * a * b4 * d5 * e * f * g3 + 60 * b5 * d4 * e * f2 * g3 +
                            6 * b4 * c * d5 * g4 + 5 * a * b4 * d6 * g4 + 6 * b5 * d5 * f * g4 -
                            15 * a * c4 * e4 * n2 - 20 * c3 * d * e4 * n2 -
                            120 * a2 * c3 * d * e4 * n2 - 150 * a * c2 * d2 * e4 * n2 -
                            300 * a3 * c2 * d2 * e4 * n2 - 300 * a2 * c * d3 * e4 * n2 -
                            300 * a4 * c * d3 * e4 * n2 - 175 * a3 * d4 * e4 * n2 -
                            105 * a5 * d4 * e4 * n2 - 60 * a * b * c3 * e4 * f * n2 -
                            60 * b * c2 * d * e4 * f * n2 - 360 * a2 * b * c2 * d * e4 * f * n2 -
                            300 * a * b * c * d2 * e4 * f * n2 -
                            600 * a3 * b * c * d2 * e4 * f * n2 - 300 * a2 * b * d3 * e4 * f * n2 -
                            300 * a4 * b * d3 * e4 * f * n2 - 90 * a * b2 * c2 * e4 * f2 * n2 -
                            60 * b2 * c * d * e4 * f2 * n2 - 360 * a2 * b2 * c * d * e4 * f2 * n2 -
                            150 * a * b2 * d2 * e4 * f2 * n2 - 300 * a3 * b2 * d2 * e4 * f2 * n2 -
                            60 * a * b3 * c * e4 * f3 * n2 - 20 * b3 * d * e4 * f3 * n2 -
                            120 * a2 * b3 * d * e4 * f3 * n2 - 15 * a * b4 * e4 * f4 * n2 -
                            20 * b * c4 * e3 * g * n2 - 240 * a * b * c3 * d * e3 * g * n2 -
                            120 * b * c2 * d2 * e3 * g * n2 - 720 * a2 * b * c2 * d2 * e3 * g * n2 -
                            400 * a * b * c * d3 * e3 * g * n2 -
                            800 * a3 * b * c * d3 * e3 * g * n2 - 300 * a2 * b * d4 * e3 * g * n2 -
                            300 * a4 * b * d4 * e3 * g * n2 - 80 * b2 * c3 * e3 * f * g * n2 -
                            720 * a * b2 * c2 * d * e3 * f * g * n2 -
                            240 * b2 * c * d2 * e3 * f * g * n2 -
                            1440 * a2 * b2 * c * d2 * e3 * f * g * n2 -
                            400 * a * b2 * d3 * e3 * f * g * n2 -
                            800 * a3 * b2 * d3 * e3 * f * g * n2 -
                            120 * b3 * c2 * e3 * f2 * g * n2 -
                            720 * a * b3 * c * d * e3 * f2 * g * n2 -
                            120 * b3 * d2 * e3 * f2 * g * n2 -
                            720 * a2 * b3 * d2 * e3 * f2 * g * n2 - 80 * b4 * c * e3 * f3 * g * n2 -
                            240 * a * b4 * d * e3 * f3 * g * n2 - 20 * b5 * e3 * f4 * g * n2 -
                            120 * b2 * c3 * d * e2 * g2 * n2 -
                            540 * a * b2 * c2 * d2 * e2 * g2 * n2 -
                            120 * b2 * c * d3 * e2 * g2 * n2 -
                            720 * a2 * b2 * c * d3 * e2 * g2 * n2 -
                            150 * a * b2 * d4 * e2 * g2 * n2 - 300 * a3 * b2 * d4 * e2 * g2 * n2 -
                            360 * b3 * c2 * d * e2 * f * g2 * n2 -
                            1080 * a * b3 * c * d2 * e2 * f * g2 * n2 -
                            120 * b3 * d3 * e2 * f * g2 * n2 -
                            720 * a2 * b3 * d3 * e2 * f * g2 * n2 -
                            360 * b4 * c * d * e2 * f2 * g2 * n2 -
                            540 * a * b4 * d2 * e2 * f2 * g2 * n2 -
                            120 * b5 * d * e2 * f3 * g2 * n2 - 120 * b3 * c2 * d2 * e * g3 * n2 -
                            240 * a * b3 * c * d3 * e * g3 * n2 - 20 * b3 * d4 * e * g3 * n2 -
                            120 * a2 * b3 * d4 * e * g3 * n2 - 240 * b4 * c * d2 * e * f * g3 * n2 -
                            240 * a * b4 * d3 * e * f * g3 * n2 - 120 * b5 * d2 * e * f2 * g3 * n2 -
                            20 * b4 * c * d3 * g4 * n2 - 15 * a * b4 * d4 * g4 * n2 -
                            20 * b5 * d3 * f * g4 * n2 + 30 * a * c2 * e4 * n4 +
                            15 * a3 * c2 * e4 * n4 + 30 * c * d * e4 * n4 +
                            120 * a2 * c * d * e4 * n4 + 30 * a4 * c * d * e4 * n4 +
                            75 * a * d2 * e4 * n4 + 100 * a3 * d2 * e4 * n4 +
                            15 * a5 * d2 * e4 * n4 + 60 * a * b * c * e4 * f * n4 +
                            30 * a3 * b * c * e4 * f * n4 + 30 * b * d * e4 * f * n4 +
                            120 * a2 * b * d * e4 * f * n4 + 30 * a4 * b * d * e4 * f * n4 +
                            30 * a * b2 * e4 * f2 * n4 + 15 * a3 * b2 * e4 * f2 * n4 +
                            40 * b * c2 * e3 * g * n4 + 60 * a2 * b * c2 * e3 * g * n4 +
                            240 * a * b * c * d * e3 * g * n4 + 120 * a3 * b * c * d * e3 * g * n4 +
                            60 * b * d2 * e3 * g * n4 + 240 * a2 * b * d2 * e3 * g * n4 +
                            60 * a4 * b * d2 * e3 * g * n4 + 80 * b2 * c * e3 * f * g * n4 +
                            120 * a2 * b2 * c * e3 * f * g * n4 +
                            240 * a * b2 * d * e3 * f * g * n4 +
                            120 * a3 * b2 * d * e3 * f * g * n4 + 40 * b3 * e3 * f2 * g * n4 +
                            60 * a2 * b3 * e3 * f2 * g * n4 + 90 * a * b2 * c2 * e2 * g2 * n4 +
                            120 * b2 * c * d * e2 * g2 * n4 + 180 * a2 * b2 * c * d * e2 * g2 * n4 +
                            180 * a * b2 * d2 * e2 * g2 * n4 + 90 * a3 * b2 * d2 * e2 * g2 * n4 +
                            180 * a * b3 * c * e2 * f * g2 * n4 + 120 * b3 * d * e2 * f * g2 * n4 +
                            180 * a2 * b3 * d * e2 * f * g2 * n4 + 90 * a * b4 * e2 * f2 * g2 * n4 +
                            60 * b3 * c2 * e * g3 * n4 + 120 * a * b3 * c * d * e * g3 * n4 +
                            40 * b3 * d2 * e * g3 * n4 + 60 * a2 * b3 * d2 * e * g3 * n4 +
                            120 * b4 * c * e * f * g3 * n4 + 120 * a * b4 * d * e * f * g3 * n4 +
                            60 * b5 * e * f2 * g3 * n4 + 30 * b4 * c * d * g4 * n4 +
                            15 * a * b4 * d2 * g4 * n4 + 30 * b5 * d * f * g4 * n4 -
                            15 * a * e4 * n6 - 5 * a3 * e4 * n6 - 20 * b * e3 * g * n6 -
                            20 * a2 * b * e3 * g * n6 - 30 * a * b2 * e2 * g2 * n6 -
                            20 * b3 * e * g3 * n6) +
                    (21 *
                        (3 * c5 * e5 + 75 * a * c4 * d * e5 + 450 * a2 * c3 * d2 * e5 +
                            1050 * a3 * c2 * d3 * e5 + 1050 * a4 * c * d4 * e5 +
                            378 * a5 * d5 * e5 + 15 * b * c4 * e5 * f +
                            300 * a * b * c3 * d * e5 * f + 1350 * a2 * b * c2 * d2 * e5 * f +
                            2100 * a3 * b * c * d3 * e5 * f + 1050 * a4 * b * d4 * e5 * f +
                            30 * b2 * c3 * e5 * f2 + 450 * a * b2 * c2 * d * e5 * f2 +
                            1350 * a2 * b2 * c * d2 * e5 * f2 + 1050 * a3 * b2 * d3 * e5 * f2 +
                            30 * b3 * c2 * e5 * f3 + 300 * a * b3 * c * d * e5 * f3 +
                            450 * a2 * b3 * d2 * e5 * f3 + 15 * b4 * c * e5 * f4 +
                            75 * a * b4 * d * e5 * f4 + 3 * b5 * e5 * f5 +
                            75 * b * c4 * d * e4 * g + 750 * a * b * c3 * d2 * e4 * g +
                            2250 * a2 * b * c2 * d3 * e4 * g + 2625 * a3 * b * c * d4 * e4 * g +
                            1050 * a4 * b * d5 * e4 * g + 300 * b2 * c3 * d * e4 * f * g +
                            2250 * a * b2 * c2 * d2 * e4 * f * g +
                            4500 * a2 * b2 * c * d3 * e4 * f * g +
                            2625 * a3 * b2 * d4 * e4 * f * g + 450 * b3 * c2 * d * e4 * f2 * g +
                            2250 * a * b3 * c * d2 * e4 * f2 * g +
                            2250 * a2 * b3 * d3 * e4 * f2 * g + 300 * b4 * c * d * e4 * f3 * g +
                            750 * a * b4 * d2 * e4 * f3 * g + 75 * b5 * d * e4 * f4 * g +
                            300 * b2 * c3 * d2 * e3 * g2 + 1500 * a * b2 * c2 * d3 * e3 * g2 +
                            2250 * a2 * b2 * c * d4 * e3 * g2 + 1050 * a3 * b2 * d5 * e3 * g2 +
                            900 * b3 * c2 * d2 * e3 * f * g2 +
                            3000 * a * b3 * c * d3 * e3 * f * g2 +
                            2250 * a2 * b3 * d4 * e3 * f * g2 + 900 * b4 * c * d2 * e3 * f2 * g2 +
                            1500 * a * b4 * d3 * e3 * f2 * g2 + 300 * b5 * d2 * e3 * f3 * g2 +
                            300 * b3 * c2 * d3 * e2 * g3 + 750 * a * b3 * c * d4 * e2 * g3 +
                            450 * a2 * b3 * d5 * e2 * g3 + 600 * b4 * c * d3 * e2 * f * g3 +
                            750 * a * b4 * d4 * e2 * f * g3 + 300 * b5 * d3 * e2 * f2 * g3 +
                            75 * b4 * c * d4 * e * g4 + 75 * a * b4 * d5 * e * g4 +
                            75 * b5 * d4 * e * f * g4 + 3 * b5 * d5 * g5 - 10 * c3 * e5 * n2 -
                            60 * a2 * c3 * e5 * n2 - 150 * a * c2 * d * e5 * n2 -
                            300 * a3 * c2 * d * e5 * n2 - 450 * a2 * c * d2 * e5 * n2 -
                            450 * a4 * c * d2 * e5 * n2 - 350 * a3 * d3 * e5 * n2 -
                            210 * a5 * d3 * e5 * n2 - 30 * b * c2 * e5 * f * n2 -
                            180 * a2 * b * c2 * e5 * f * n2 - 300 * a * b * c * d * e5 * f * n2 -
                            600 * a3 * b * c * d * e5 * f * n2 - 450 * a2 * b * d2 * e5 * f * n2 -
                            450 * a4 * b * d2 * e5 * f * n2 - 30 * b2 * c * e5 * f2 * n2 -
                            180 * a2 * b2 * c * e5 * f2 * n2 - 150 * a * b2 * d * e5 * f2 * n2 -
                            300 * a3 * b2 * d * e5 * f2 * n2 - 10 * b3 * e5 * f3 * n2 -
                            60 * a2 * b3 * e5 * f3 * n2 - 150 * a * b * c3 * e4 * g * n2 -
                            150 * b * c2 * d * e4 * g * n2 - 900 * a2 * b * c2 * d * e4 * g * n2 -
                            750 * a * b * c * d2 * e4 * g * n2 -
                            1500 * a3 * b * c * d2 * e4 * g * n2 - 750 * a2 * b * d3 * e4 * g * n2 -
                            750 * a4 * b * d3 * e4 * g * n2 - 450 * a * b2 * c2 * e4 * f * g * n2 -
                            300 * b2 * c * d * e4 * f * g * n2 -
                            1800 * a2 * b2 * c * d * e4 * f * g * n2 -
                            750 * a * b2 * d2 * e4 * f * g * n2 -
                            1500 * a3 * b2 * d2 * e4 * f * g * n2 -
                            450 * a * b3 * c * e4 * f2 * g * n2 - 150 * b3 * d * e4 * f2 * g * n2 -
                            900 * a2 * b3 * d * e4 * f2 * g * n2 - 150 * a * b4 * e4 * f3 * g * n2 -
                            100 * b2 * c3 * e3 * g2 * n2 - 900 * a * b2 * c2 * d * e3 * g2 * n2 -
                            300 * b2 * c * d2 * e3 * g2 * n2 -
                            1800 * a2 * b2 * c * d2 * e3 * g2 * n2 -
                            500 * a * b2 * d3 * e3 * g2 * n2 - 1000 * a3 * b2 * d3 * e3 * g2 * n2 -
                            300 * b3 * c2 * e3 * f * g2 * n2 -
                            1800 * a * b3 * c * d * e3 * f * g2 * n2 -
                            300 * b3 * d2 * e3 * f * g2 * n2 -
                            1800 * a2 * b3 * d2 * e3 * f * g2 * n2 -
                            300 * b4 * c * e3 * f2 * g2 * n2 -
                            900 * a * b4 * d * e3 * f2 * g2 * n2 - 100 * b5 * e3 * f3 * g2 * n2 -
                            300 * b3 * c2 * d * e2 * g3 * n2 -
                            900 * a * b3 * c * d2 * e2 * g3 * n2 - 100 * b3 * d3 * e2 * g3 * n2 -
                            600 * a2 * b3 * d3 * e2 * g3 * n2 -
                            600 * b4 * c * d * e2 * f * g3 * n2 -
                            900 * a * b4 * d2 * e2 * f * g3 * n2 -
                            300 * b5 * d * e2 * f2 * g3 * n2 - 150 * b4 * c * d2 * e * g4 * n2 -
                            150 * a * b4 * d3 * e * g4 * n2 - 150 * b5 * d2 * e * f * g4 * n2 -
                            10 * b5 * d3 * g5 * n2 + 15 * c * e5 * n4 + 60 * a2 * c * e5 * n4 +
                            15 * a4 * c * e5 * n4 + 75 * a * d * e5 * n4 + 100 * a3 * d * e5 * n4 +
                            15 * a5 * d * e5 * n4 + 15 * b * e5 * f * n4 +
                            60 * a2 * b * e5 * f * n4 + 15 * a4 * b * e5 * f * n4 +
                            150 * a * b * c * e4 * g * n4 + 75 * a3 * b * c * e4 * g * n4 +
                            75 * b * d * e4 * g * n4 + 300 * a2 * b * d * e4 * g * n4 +
                            75 * a4 * b * d * e4 * g * n4 + 150 * a * b2 * e4 * f * g * n4 +
                            75 * a3 * b2 * e4 * f * g * n4 + 100 * b2 * c * e3 * g2 * n4 +
                            150 * a2 * b2 * c * e3 * g2 * n4 + 300 * a * b2 * d * e3 * g2 * n4 +
                            150 * a3 * b2 * d * e3 * g2 * n4 + 100 * b3 * e3 * f * g2 * n4 +
                            150 * a2 * b3 * e3 * f * g2 * n4 + 150 * a * b3 * c * e2 * g3 * n4 +
                            100 * b3 * d * e2 * g3 * n4 + 150 * a2 * b3 * d * e2 * g3 * n4 +
                            150 * a * b4 * e2 * f * g3 * n4 + 75 * b4 * c * e * g4 * n4 +
                            75 * a * b4 * d * e * g4 * n4 + 75 * b5 * e * f * g4 * n4 +
                            15 * b5 * d * g5 * n4)) /
                        2 +
                    (70 *
                        (15 * a * c4 * e6 + 180 * a2 * c3 * d * e6 + 630 * a3 * c2 * d2 * e6 +
                            840 * a4 * c * d3 * e6 + 378 * a5 * d4 * e6 + 60 * a * b * c3 * e6 * f +
                            540 * a2 * b * c2 * d * e6 * f + 1260 * a3 * b * c * d2 * e6 * f +
                            840 * a4 * b * d3 * e6 * f + 90 * a * b2 * c2 * e6 * f2 +
                            540 * a2 * b2 * c * d * e6 * f2 + 630 * a3 * b2 * d2 * e6 * f2 +
                            60 * a * b3 * c * e6 * f3 + 180 * a2 * b3 * d * e6 * f3 +
                            15 * a * b4 * e6 * f4 + 18 * b * c4 * e5 * g +
                            360 * a * b * c3 * d * e5 * g + 1620 * a2 * b * c2 * d2 * e5 * g +
                            2520 * a3 * b * c * d3 * e5 * g + 1260 * a4 * b * d4 * e5 * g +
                            72 * b2 * c3 * e5 * f * g + 1080 * a * b2 * c2 * d * e5 * f * g +
                            3240 * a2 * b2 * c * d2 * e5 * f * g +
                            2520 * a3 * b2 * d3 * e5 * f * g + 108 * b3 * c2 * e5 * f2 * g +
                            1080 * a * b3 * c * d * e5 * f2 * g +
                            1620 * a2 * b3 * d2 * e5 * f2 * g + 72 * b4 * c * e5 * f3 * g +
                            360 * a * b4 * d * e5 * f3 * g + 18 * b5 * e5 * f4 * g +
                            180 * b2 * c3 * d * e4 * g2 + 1350 * a * b2 * c2 * d2 * e4 * g2 +
                            2700 * a2 * b2 * c * d3 * e4 * g2 + 1575 * a3 * b2 * d4 * e4 * g2 +
                            540 * b3 * c2 * d * e4 * f * g2 + 2700 * a * b3 * c * d2 * e4 * f * g2 +
                            2700 * a2 * b3 * d3 * e4 * f * g2 + 540 * b4 * c * d * e4 * f2 * g2 +
                            1350 * a * b4 * d2 * e4 * f2 * g2 + 180 * b5 * d * e4 * f3 * g2 +
                            360 * b3 * c2 * d2 * e3 * g3 + 1200 * a * b3 * c * d3 * e3 * g3 +
                            900 * a2 * b3 * d4 * e3 * g3 + 720 * b4 * c * d2 * e3 * f * g3 +
                            1200 * a * b4 * d3 * e3 * f * g3 + 360 * b5 * d2 * e3 * f2 * g3 +
                            180 * b4 * c * d3 * e2 * g4 + 225 * a * b4 * d4 * e2 * g4 +
                            180 * b5 * d3 * e2 * f * g4 + 18 * b5 * d4 * e * g5 -
                            30 * a * c2 * e6 * n2 - 60 * a3 * c2 * e6 * n2 -
                            180 * a2 * c * d * e6 * n2 - 180 * a4 * c * d * e6 * n2 -
                            210 * a3 * d2 * e6 * n2 - 126 * a5 * d2 * e6 * n2 -
                            60 * a * b * c * e6 * f * n2 - 120 * a3 * b * c * e6 * f * n2 -
                            180 * a2 * b * d * e6 * f * n2 - 180 * a4 * b * d * e6 * f * n2 -
                            30 * a * b2 * e6 * f2 * n2 - 60 * a3 * b2 * e6 * f2 * n2 -
                            36 * b * c2 * e5 * g * n2 - 216 * a2 * b * c2 * e5 * g * n2 -
                            360 * a * b * c * d * e5 * g * n2 - 720 * a3 * b * c * d * e5 * g * n2 -
                            540 * a2 * b * d2 * e5 * g * n2 - 540 * a4 * b * d2 * e5 * g * n2 -
                            72 * b2 * c * e5 * f * g * n2 - 432 * a2 * b2 * c * e5 * f * g * n2 -
                            360 * a * b2 * d * e5 * f * g * n2 -
                            720 * a3 * b2 * d * e5 * f * g * n2 - 36 * b3 * e5 * f2 * g * n2 -
                            216 * a2 * b3 * e5 * f2 * g * n2 - 270 * a * b2 * c2 * e4 * g2 * n2 -
                            180 * b2 * c * d * e4 * g2 * n2 -
                            1080 * a2 * b2 * c * d * e4 * g2 * n2 -
                            450 * a * b2 * d2 * e4 * g2 * n2 - 900 * a3 * b2 * d2 * e4 * g2 * n2 -
                            540 * a * b3 * c * e4 * f * g2 * n2 - 180 * b3 * d * e4 * f * g2 * n2 -
                            1080 * a2 * b3 * d * e4 * f * g2 * n2 -
                            270 * a * b4 * e4 * f2 * g2 * n2 - 120 * b3 * c2 * e3 * g3 * n2 -
                            720 * a * b3 * c * d * e3 * g3 * n2 - 120 * b3 * d2 * e3 * g3 * n2 -
                            720 * a2 * b3 * d2 * e3 * g3 * n2 - 240 * b4 * c * e3 * f * g3 * n2 -
                            720 * a * b4 * d * e3 * f * g3 * n2 - 120 * b5 * e3 * f2 * g3 * n2 -
                            180 * b4 * c * d * e2 * g4 * n2 - 270 * a * b4 * d2 * e2 * g4 * n2 -
                            180 * b5 * d * e2 * f * g4 * n2 - 36 * b5 * d2 * e * g5 * n2 +
                            15 * a * e6 * n4 + 20 * a3 * e6 * n4 + 3 * a5 * e6 * n4 +
                            18 * b * e5 * g * n4 + 72 * a2 * b * e5 * g * n4 +
                            18 * a4 * b * e5 * g * n4 + 90 * a * b2 * e4 * g2 * n4 +
                            45 * a3 * b2 * e4 * g2 * n4 + 40 * b3 * e3 * g3 * n4 +
                            60 * a2 * b3 * e3 * g3 * n4 + 45 * a * b4 * e2 * g4 * n4 +
                            18 * b5 * e * g5 * n4)) /
                        9 +
                    12 * (15 * a2 * c3 * e7 + 105 * a3 * c2 * d * e7 + 210 * a4 * c * d2 * e7 +
                             126 * a5 * d3 * e7 + 45 * a2 * b * c2 * e7 * f +
                             210 * a3 * b * c * d * e7 * f + 210 * a4 * b * d2 * e7 * f +
                             45 * a2 * b2 * c * e7 * f2 + 105 * a3 * b2 * d * e7 * f2 +
                             15 * a2 * b3 * e7 * f3 + 35 * a * b * c3 * e6 * g +
                             315 * a2 * b * c2 * d * e6 * g + 735 * a3 * b * c * d2 * e6 * g +
                             490 * a4 * b * d3 * e6 * g + 105 * a * b2 * c2 * e6 * f * g +
                             630 * a2 * b2 * c * d * e6 * f * g + 735 * a3 * b2 * d2 * e6 * f * g +
                             105 * a * b3 * c * e6 * f2 * g + 315 * a2 * b3 * d * e6 * f2 * g +
                             35 * a * b4 * e6 * f3 * g + 21 * b2 * c3 * e5 * g2 +
                             315 * a * b2 * c2 * d * e5 * g2 + 945 * a2 * b2 * c * d2 * e5 * g2 +
                             735 * a3 * b2 * d3 * e5 * g2 + 63 * b3 * c2 * e5 * f * g2 +
                             630 * a * b3 * c * d * e5 * f * g2 + 945 * a2 * b3 * d2 * e5 * f * g2 +
                             63 * b4 * c * e5 * f2 * g2 + 315 * a * b4 * d * e5 * f2 * g2 +
                             21 * b5 * e5 * f3 * g2 + 105 * b3 * c2 * d * e4 * g3 +
                             525 * a * b3 * c * d2 * e4 * g3 + 525 * a2 * b3 * d3 * e4 * g3 +
                             210 * b4 * c * d * e4 * f * g3 + 525 * a * b4 * d2 * e4 * f * g3 +
                             105 * b5 * d * e4 * f2 * g3 + 105 * b4 * c * d2 * e3 * g4 +
                             175 * a * b4 * d3 * e3 * g4 + 105 * b5 * d2 * e3 * f * g4 +
                             21 * b5 * d3 * e2 * g5 - 15 * a2 * c * e7 * n2 -
                             15 * a4 * c * e7 * n2 - 35 * a3 * d * e7 * n2 - 21 * a5 * d * e7 * n2 -
                             15 * a2 * b * e7 * f * n2 - 15 * a4 * b * e7 * f * n2 -
                             35 * a * b * c * e6 * g * n2 - 70 * a3 * b * c * e6 * g * n2 -
                             105 * a2 * b * d * e6 * g * n2 - 105 * a4 * b * d * e6 * g * n2 -
                             35 * a * b2 * e6 * f * g * n2 - 70 * a3 * b2 * e6 * f * g * n2 -
                             21 * b2 * c * e5 * g2 * n2 - 126 * a2 * b2 * c * e5 * g2 * n2 -
                             105 * a * b2 * d * e5 * g2 * n2 - 210 * a3 * b2 * d * e5 * g2 * n2 -
                             21 * b3 * e5 * f * g2 * n2 - 126 * a2 * b3 * e5 * f * g2 * n2 -
                             105 * a * b3 * c * e4 * g3 * n2 - 35 * b3 * d * e4 * g3 * n2 -
                             210 * a2 * b3 * d * e4 * g3 * n2 - 105 * a * b4 * e4 * f * g3 * n2 -
                             35 * b4 * c * e3 * g4 * n2 - 105 * a * b4 * d * e3 * g4 * n2 -
                             35 * b5 * e3 * f * g4 * n2 - 21 * b5 * d * e2 * g5 * n2) +
                    (15 * (105 * a3 * c2 * e8 + 420 * a4 * c * d * e8 + 378 * a5 * d2 * e8 +
                              210 * a3 * b * c * e8 * f + 420 * a4 * b * d * e8 * f +
                              105 * a3 * b2 * e8 * f2 + 360 * a2 * b * c2 * e7 * g +
                              1680 * a3 * b * c * d * e7 * g + 1680 * a4 * b * d2 * e7 * g +
                              720 * a2 * b2 * c * e7 * f * g + 1680 * a3 * b2 * d * e7 * f * g +
                              360 * a2 * b3 * e7 * f2 * g + 420 * a * b2 * c2 * e6 * g2 +
                              2520 * a2 * b2 * c * d * e6 * g2 + 2940 * a3 * b2 * d2 * e6 * g2 +
                              840 * a * b3 * c * e6 * f * g2 + 2520 * a2 * b3 * d * e6 * f * g2 +
                              420 * a * b4 * e6 * f2 * g2 + 168 * b3 * c2 * e5 * g3 +
                              1680 * a * b3 * c * d * e5 * g3 + 2520 * a2 * b3 * d2 * e5 * g3 +
                              336 * b4 * c * e5 * f * g3 + 1680 * a * b4 * d * e5 * f * g3 +
                              168 * b5 * e5 * f2 * g3 + 420 * b4 * c * d * e4 * g4 +
                              1050 * a * b4 * d2 * e4 * g4 + 420 * b5 * d * e4 * f * g4 +
                              168 * b5 * d2 * e3 * g5 - 35 * a3 * e8 * n2 - 21 * a5 * e8 * n2 -
                              120 * a2 * b * e7 * g * n2 - 120 * a4 * b * e7 * g * n2 -
                              140 * a * b2 * e6 * g2 * n2 - 280 * a3 * b2 * e6 * g2 * n2 -
                              56 * b3 * e5 * g3 * n2 - 336 * a2 * b3 * e5 * g3 * n2 -
                              210 * a * b4 * e4 * g4 * n2 - 56 * b5 * e3 * g5 * n2)) /
                        11 +
                    (5 * (70 * a4 * c * e9 + 126 * a5 * d * e9 + 70 * a4 * b * e9 * f +
                             315 * a3 * b * c * e8 * g + 630 * a4 * b * d * e8 * g +
                             315 * a3 * b2 * e8 * f * g + 540 * a2 * b2 * c * e7 * g2 +
                             1260 * a3 * b2 * d * e7 * g2 + 540 * a2 * b3 * e7 * f * g2 +
                             420 * a * b3 * c * e6 * g3 + 1260 * a2 * b3 * d * e6 * g3 +
                             420 * a * b4 * e6 * f * g3 + 126 * b4 * c * e5 * g4 +
                             630 * a * b4 * d * e5 * g4 + 126 * b5 * e5 * f * g4 +
                             126 * b5 * d * e4 * g5)) /
                        6 +
                    (e5 * (126 * a5 * e5 + 700 * a4 * b * e4 * g + 1575 * a3 * b2 * e3 * g2 +
                              1800 * a2 * b3 * e2 * g3 + 1050 * a * b4 * e * g4 + 252 * b5 * g5)) /
                        13)) /
                (3150 * n15));

    vol_ele += linecontribution;
  }

  return;
}


/*----------------------------------------------------------------------*
 | test for overlap of the XAABB of bubble and element     ghamm 11/14  |
 *----------------------------------------------------------------------*/
bool CAVITATION::Algorithm::XAABBoverlap(DRT::Element* ele, const double influence,
    const LINALG::Matrix<3, 1>& particleposition, const bool havepbc, bool& pbcdetected,
    LINALG::Matrix<3, 1>& pbceleoffset)
{
  pbcdetected = false;
  if (havepbc)
  {
    static LINALG::Matrix<3, 1> distance;
    pbceleoffset.PutScalar(0.0);

    for (unsigned idim = 0; idim < 3; ++idim)
    {
      // check distance of fluid element and particle position roughly
      distance(idim) = ele->Nodes()[0]->X()[idim] - particleposition(idim);
      // use heuristic to find out whether fluid element and bubble are connected via a pbc bound
      if (BinStrategy()->HavePBC(idim) and
          (std::abs(distance(idim)) > BinStrategy()->PBCDelta(idim) * 0.5))
      {
        // fluid element and bubble are (hopefully) connected via a pbc bound
        pbcdetected = true;

        // compute necessary offset to account for pbc
        if (distance(idim) > 0.0)
          pbceleoffset(idim) = -BinStrategy()->PBCDelta(idim);
        else
          pbceleoffset(idim) = +BinStrategy()->PBCDelta(idim);
      }
    }
  }

  // get bounding box of current element
  double xaabb[6] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
  for (int inode = 0; inode < ele->NumNode(); ++inode)
  {
    const DRT::Node* node = ele->Nodes()[inode];
    const double* coord = node->X();
    for (size_t i = 0; i < 3; ++i)
    {
      xaabb[i + 0] = std::min(xaabb[i + 0], coord[i]);
      xaabb[i + 3] = std::max(xaabb[i + 3], coord[i]);
    }
  }

  double bubblesurface[6] = {particleposition(0) + influence, particleposition(1) + influence,
      particleposition(2) + influence, particleposition(0) - influence,
      particleposition(1) - influence, particleposition(2) - influence};

  // modify xaabb in case of pbc
  if (pbcdetected)
  {
    for (size_t i = 0; i < 3; ++i)
    {
      xaabb[i] += pbceleoffset(i);
      xaabb[3 + i] += pbceleoffset(i);
    }
  }

  bool overlap = true;
  // test whether the bounding box of the fluid element touches the bubble influence
  for (int dim = 0; dim < 3; ++dim)
  {
    if (xaabb[dim] - GEO::TOL7 > bubblesurface[dim] or
        xaabb[dim + 3] + GEO::TOL7 < bubblesurface[dim + 3])
    {
      overlap = false;
      break;
    }
  }

  return overlap;
}


/*----------------------------------------------------------------------*
 | assign bubble volume to the underlying fluid element    ghamm 11/14  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::AssignSmallBubbles(const double bubblevol,
    const LINALG::Matrix<3, 1>& particleposition, const std::vector<int>& insideeles,
    Teuchos::RCP<Epetra_FEVector> void_volumes, const std::vector<bool>& pbcdetected,
    const std::vector<LINALG::Matrix<3, 1>>& pbceleoffset)
{
  // safety check
  if (insideeles.size() < 1)
    dserror("no underlying fluid element for void fraction computation found.");

  // a very small bubble was processed which cannot be resolved with volume integration
  // just add the contribution to the underlying fluid element
  if (insideeles.size() == 1)
  {
    // in case of pbcs: it is faster to modify the particle position than the element coordinates
    LINALG::Matrix<3, 1> modifiedparticlepos(particleposition);
    if (pbcdetected[0]) modifiedparticlepos.Update(-1.0, pbceleoffset[0], 1.0);

    // do assembly of void fraction into fluid
    int err = void_volumes->SumIntoGlobalValues(1, &insideeles[0], &bubblevol);
    if (err < 0) dserror("summing into Epetra_FEVector failed");

    // safety check, may be moved to debug version only
    {
      DRT::Element* fluidele = fluiddis_->gElement(insideeles[0]);

      static LINALG::Matrix<3, 1> dummy;
      // get coordinates of the particle position in parameter space of the element
      bool insideele = GEO::currentToVolumeElementCoordinates(
          fluidele->Shape(), xyze_cache_[fluidele->LID()], modifiedparticlepos, dummy);

      if (insideele == false)
        dserror(
            "bubble at position x: %f y: %f z: %f is expected to lie inside fluid "
            "element with id: %d",
            particleposition(0), particleposition(1), particleposition(2), fluidele->Id());
    }
  }
  else
  {
    bool insideele = false;
    for (size_t i = 0; i < insideeles.size(); ++i)
    {
      // in case of pbcs: it is faster to modify the particle position than the element coordinates
      LINALG::Matrix<3, 1> modifiedparticlepos(particleposition);
      if (pbcdetected[i]) modifiedparticlepos.Update(-1.0, pbceleoffset[i], 1.0);

      DRT::Element* fluidele = fluiddis_->gElement(insideeles[i]);

      static LINALG::Matrix<3, 1> dummy;
      // get coordinates of the particle position in parameter space of the element
      insideele = GEO::currentToVolumeElementCoordinates(
          fluidele->Shape(), xyze_cache_[fluidele->LID()], modifiedparticlepos, dummy);

      if (insideele == true)
      {
        // do assembly of void fraction into fluid
        int err = void_volumes->SumIntoGlobalValues(1, &insideeles[i], &bubblevol);
        if (err < 0) dserror("summing into Epetra_FEVector failed");
        break;
      }
    }
    if (insideele == false)
      dserror(
          "void fraction for particle (position x: %d y: %d z: %d) could not be assigned "
          "to an underlying fluid element.",
          particleposition(0), particleposition(1), particleposition(2));
  }

  return;
}


/*----------------------------------------------------------------------*
 | get underlying element as well as position in element space          |
 |                                                          ghamm 06/15 |
 *----------------------------------------------------------------------*/
DRT::Element* CAVITATION::Algorithm::GetEleCoordinatesFromPosition(
    const DRT::Node* currparticle,             ///< particle
    const LINALG::Matrix<3, 1>& myposition,    ///< position
    DRT::MESHFREE::MeshfreeMultiBin* currbin,  ///< corresponding bin
    LINALG::Matrix<3, 1>&
        elecoord,  ///< matrix to be filled with particle coordinates in element space
    const bool approxelecoordsinit  ///< bool whether an inital approx. of the ele coords is used
)
{
  // variables to store information about element in which the particle is located
  elecoord.Clear();
  bool insideele = false;

  DRT::Element** fluidelesinbin = currbin->AssociatedEles(bin_fluidcontent_);
  int numfluidelesinbin = currbin->NumAssociatedEle(bin_fluidcontent_);

  std::set<int>::const_iterator eleiter;
  // search for underlying fluid element with fast search if desired
  for (int ele = 0; ele < numfluidelesinbin; ++ele)
  {
    DRT::Element* fluidele = fluidelesinbin[ele];

    // get coordinates of the particle position in parameter space of the element
    insideele = GEO::currentToVolumeElementCoordinates(
        fluidele->Shape(), xyze_cache_[fluidele->LID()], myposition, elecoord, approxelecoordsinit);

    if (insideele == true)
    {
      // leave loop over all fluid eles in bin
      return fluidele;
    }
  }

  // repeat search for underlying fluid element with standard search in case nothing was found
  if (approxelecoordsinit == true)
  {
    for (int ele = 0; ele < numfluidelesinbin; ++ele)
    {
      DRT::Element* fluidele = fluidelesinbin[ele];

      // get coordinates of the particle position in parameter space of the element
      insideele = GEO::currentToVolumeElementCoordinates(
          fluidele->Shape(), xyze_cache_[fluidele->LID()], myposition, elecoord, false);

      if (insideele == true)
      {
        // leave loop over all fluid eles in bin
        return fluidele;
      }
    }
  }

  // repeat search for the underlying fluid element in neighborhood of current bin
  // this should only be necessary on very rare cases

  // gather neighboring bins
  int ijk[3];
  BinStrategy()->ConvertGidToijk(currbin->Id(), ijk);

  // ijk_range contains: i_min   i_max     j_min     j_max    k_min     k_max
  const int ijk_range[] = {ijk[0] - 1, ijk[0] + 1, ijk[1] - 1, ijk[1] + 1, ijk[2] - 1, ijk[2] + 1};
  std::vector<int> binIds;
  binIds.reserve(27);

  // check on existence here
  BinStrategy()->GidsInijkRange(ijk_range, binIds, true);

  for (size_t b = 0; b < binIds.size(); ++b)
  {
    currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(
        BinStrategy()->BinDiscret()->gElement(binIds[b]));

    fluidelesinbin = currbin->AssociatedEles(bin_fluidcontent_);
    numfluidelesinbin = currbin->NumAssociatedEle(bin_fluidcontent_);

    for (int ele = 0; ele < numfluidelesinbin; ++ele)
    {
      DRT::Element* fluidele = fluidelesinbin[ele];

      // get coordinates of the particle position in parameter space of the element
      insideele = GEO::currentToVolumeElementCoordinates(
          fluidele->Shape(), xyze_cache_[fluidele->LID()], myposition, elecoord, false);

      if (insideele == true)
      {
        std::cout << "INFO: underlying fluid element (id: " << fluidele->Id()
                  << ") for currparticle with Id: " << currparticle->Id()
                  << " and position: " << myposition(0) << " " << myposition(1) << " "
                  << myposition(2) << " was found in a neighboring bin: " << currbin->Id()
                  << " on proc " << MyRank() << std::endl;

        // leave loop over all fluid eles in bin
        return fluidele;
      }
    }
  }

  // no underlying fluid element was found
  return NULL;
}


/*----------------------------------------------------------------------*
 | compute velocity at bubble position                     ghamm 12/15  |
 *----------------------------------------------------------------------*/
bool CAVITATION::Algorithm::ComputeVelocityAtBubblePosition(DRT::Node* currparticle,
    LINALG::Matrix<3, 1>& particleposition, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // find out in which fluid element the current particle is located
  if (currparticle->NumElement() != 1)
    dserror("ERROR: A particle is assigned to more than one bin!");
  DRT::Element** currele = currparticle->Elements();
  DRT::MESHFREE::MeshfreeMultiBin* currbin =
      dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);

  static LINALG::Matrix<3, 1> elecoord;
  DRT::Element* targetfluidele = GetEleCoordinatesFromPosition(
      currparticle, particleposition, currbin, elecoord, approxelecoordsinit_);

  if (targetfluidele == NULL)
  {
    std::cout << "WARNING: velocity for bubble (id: " << currparticle->Id()
              << " ) could not be computed at position: " << particleposition(0) << " "
              << particleposition(1) << " " << particleposition(2) << " on proc " << MyRank()
              << std::endl;
    return false;
  }

  // get element location vector and ownerships
  std::vector<int> lm_f;
  std::vector<int> lmowner_f;
  std::vector<int> lmstride;
  targetfluidele->LocationVector(*fluiddis_, lm_f, lmowner_f, lmstride);

  // set action in order to compute velocity -> state with name "vel" expected inside
  Teuchos::ParameterList params;
  params.set<int>("action", FLD::interpolate_velocity_to_given_point);
  params.set<LINALG::Matrix<3, 1>>("elecoords", elecoord);

  // call the element specific evaluate method (elevec1 = fluid vel)
  targetfluidele->Evaluate(params, *fluiddis_, lm_f, elemat1, elemat2, elevec1, elevec2, elevec3);

  // enforce 2D bubble movement for pseudo-2D problem
  if (BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_2Dz) elevec1(2) = 0.0;

  return true;
}


/*----------------------------------------------------------------------*
 | compute pressure at bubble position                     ghamm 06/15  |
 *----------------------------------------------------------------------*/
bool CAVITATION::Algorithm::ComputePressureAtBubblePosition(DRT::Node* currparticle,
    const LINALG::Matrix<3, 1>& particleposition, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // find out in which fluid element the current particle is located
  if (currparticle->NumElement() != 1)
    dserror("ERROR: A particle is assigned to more than one bin!");
  DRT::Element** currele = currparticle->Elements();
  DRT::MESHFREE::MeshfreeMultiBin* currbin =
      dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);

  static LINALG::Matrix<3, 1> elecoord;
  DRT::Element* targetfluidele = GetEleCoordinatesFromPosition(
      currparticle, particleposition, currbin, elecoord, approxelecoordsinit_);

  // fill cache for force computation
  const size_t i = currparticle->LID();
  if (underlyingelecache_.size() > i)
  {
    UnderlyingEle& e = underlyingelecache_[i];
    e.ele = targetfluidele;
    e.elecoord = elecoord;
  }

  if (targetfluidele == NULL)
  {
    std::cout << "WARNING: pressure for bubble (id: " << currparticle->Id()
              << " ) could not be computed at position: " << particleposition(0) << " "
              << particleposition(1) << " " << particleposition(2) << " on proc " << MyRank()
              << std::endl;
    return false;
  }

  switch (targetfluidele->Shape())
  {
    case DRT::Element::hex8:
      ComputePressureAtBubblePositionT<DRT::Element::hex8>(
          targetfluidele, elecoord, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    case DRT::Element::tet4:
      ComputePressureAtBubblePositionT<DRT::Element::tet4>(
          targetfluidele, elecoord, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    default:
      dserror("add desired 3D element type here");
      break;
  }

  return true;
}


/*----------------------------------------------------------------------*
 | compute pressure at bubble position                     ghamm 06/15  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void CAVITATION::Algorithm::ComputePressureAtBubblePositionT(const DRT::Element* targetfluidele,
    const LINALG::Matrix<3, 1>& elecoord, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
#ifdef INLINED_ELE_EVAL
  {
    const int nsd_ = 3;
    static const int numdofpernode_ = nsd_ + 1;
    static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! coordinates of current integration point in reference coordinates
    static LINALG::Matrix<nsd_, 1> xsi_;
    //! node coordinates
    static LINALG::Matrix<nen_, 1> funct_;
    //! pressure at given point
    static LINALG::Matrix<nen_, 1> epre;

    // coordinates of the current point
    xsi_.Update(elecoord);

    // shape functions
    DRT::UTILS::shape_function<distype>(xsi_, funct_);

    //----------------------------------------------------------------------------
    //   Extract pressure from global vectors and compute pressure at point
    //----------------------------------------------------------------------------


    // extract local values of the global vectors
    const std::vector<int>& lm_f = lm_cache_[targetfluidele->LID()];
    std::vector<double> myvel(lm_f.size());

    if (fluiddis_->HasState("vel"))
    {
      // fill the local element vector with the global values
      DRT::UTILS::ExtractMyValues(*fluiddis_->GetState("vel"), myvel, lm_f);

      for (int inode = 0; inode < nen_; ++inode)  // number of nodes
      {
        // fill a scalar field via a pointer
        epre(inode, 0) = myvel[nsd_ + (inode * numdofpernode_)];
      }

      elevec1[0] = funct_.Dot(epre);
    }

    if (fluiddis_->HasState("velnp"))
    {
      // fill the local element vector with the global values
      // fill the local element vector with the global values
      DRT::UTILS::ExtractMyValues(*fluiddis_->GetState("velnp"), myvel, lm_f);

      for (int inode = 0; inode < nen_; ++inode)  // number of nodes
      {
        // fill a scalar field via a pointer
        epre(inode, 0) = myvel[nsd_ + (inode * numdofpernode_)];
      }

      if (elevec1.Length() != 2) dserror("velnp is set, there must be a vel as well");

      elevec1[1] = funct_.Dot(epre);
    }
  }
#else
  {
    // set action in order to compute pressure -> state with name "veln" (and "velnp") expected
    // inside
    Teuchos::ParameterList params;
    params.set<int>("action", FLD::interpolate_pressure_to_given_point);
    params.set<LINALG::Matrix<3, 1>>("elecoords", elecoord);

    // call the element specific evaluate method (elevec1 = fluid press)
    targetfluidele->Evaluate(params, *fluiddis_, lm_cache_[targetfluidele->LID()], elemat1, elemat2,
        elevec1, elevec2, elevec3);
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
 | unwarp fluid element for analytical integration         ghamm 08/16  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::UnwarpElement(const DRT::Element* ele,
    const std::vector<Teuchos::RCP<DRT::Element>>& surfaces, const int numsurfacenodes,
    std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
    std::vector<LINALG::Matrix<3, 1>>& normals)
{
  // find all auxiliary planes for the surfaces
  const int numsurf = ele->NumSurface();
  std::vector<LINALG::Matrix<3, 1>> elecenters(numsurf);
  std::vector<double> ds(numsurf);
  // loop over surfaces of current fluid element and compute auxiliary planes
  for (int isurface = 0; isurface < numsurf; ++isurface)
  {
    Teuchos::RCP<DRT::Element> currsurf = surfaces[isurface];

    std::vector<LINALG::Matrix<3, 1>> surfacenodes(numsurfacenodes);
    static LINALG::Matrix<3, 1> centerauxsurf;
    centerauxsurf.PutScalar(0.0);

    // loop over nodes of current surface
    for (int inode = 0; inode < numsurfacenodes; ++inode)
    {
      const DRT::Node* node = currsurf->Nodes()[inode];
      LINALG::Matrix<3, 1>& position = surfacenodes[inode];

      // fill position with node.X()
      for (int dim = 0; dim < 3; ++dim) position(dim) = node->X()[dim];

      // compute center of auxiliary plane
      centerauxsurf.Update(1.0, position, 1.0);
    }

    centerauxsurf.Scale(1.0 / numsurfacenodes);
    elecenters[isurface] = centerauxsurf;

    // compute auxiliary plane
    static LINALG::Matrix<3, 1> n, r, p;
    r.Update(1.0, surfacenodes[0], -1.0, surfacenodes[2]);
    p.Update(1.0, surfacenodes[1], -1.0, surfacenodes[3 % numsurfacenodes]);

    // calculate normal vector n with cross product of diagonals r,p
    n(0) = r(1) * p(2) - r(2) * p(1);
    n(1) = r(2) * p(0) - r(0) * p(2);
    n(2) = r(0) * p(1) - r(1) * p(0);
    n.Scale(1.0 / n.Norm2());
    normals[isurface] = n;

    // eqn of plane a*nx + b*ny + c*nz + d = 0
    const double d = -(n(0) * centerauxsurf(0) + n(1) * centerauxsurf(1) + n(2) * centerauxsurf(2));
    ds[isurface] = d;
  }

  // store unwarped node positions
  switch (ele->Shape())
  {
    case DRT::Element::hex8:
    {
      // compute intersection points of three planes to unwarp hex8 element

      /// node 0
      static LINALG::Matrix<3, 3> A;
      static LINALG::Matrix<3, 1> b;
      static LINALG::Matrix<3, 1> x;

      A(0, 0) = normals[0](0);
      A(0, 1) = normals[0](1);
      A(0, 2) = normals[0](2);
      A(1, 0) = normals[1](0);
      A(1, 1) = normals[1](1);
      A(1, 2) = normals[1](2);
      A(2, 0) = normals[4](0);
      A(2, 1) = normals[4](1);
      A(2, 2) = normals[4](2);
      b(0) = -ds[0];
      b(1) = -ds[1];
      b(2) = -ds[4];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[0]] = x;

      /// node 1
      A(0, 0) = normals[0](0);
      A(0, 1) = normals[0](1);
      A(0, 2) = normals[0](2);
      A(1, 0) = normals[1](0);
      A(1, 1) = normals[1](1);
      A(1, 2) = normals[1](2);
      A(2, 0) = normals[2](0);
      A(2, 1) = normals[2](1);
      A(2, 2) = normals[2](2);
      b(0) = -ds[0];
      b(1) = -ds[1];
      b(2) = -ds[2];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[1]] = x;

      /// node 2
      x.PutScalar(0.0);
      A(0, 0) = normals[0](0);
      A(0, 1) = normals[0](1);
      A(0, 2) = normals[0](2);
      A(1, 0) = normals[2](0);
      A(1, 1) = normals[2](1);
      A(1, 2) = normals[2](2);
      A(2, 0) = normals[3](0);
      A(2, 1) = normals[3](1);
      A(2, 2) = normals[3](2);
      b(0) = -ds[0];
      b(1) = -ds[2];
      b(2) = -ds[3];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[2]] = x;

      /// node 3
      x.PutScalar(0.0);
      A(0, 0) = normals[0](0);
      A(0, 1) = normals[0](1);
      A(0, 2) = normals[0](2);
      A(1, 0) = normals[3](0);
      A(1, 1) = normals[3](1);
      A(1, 2) = normals[3](2);
      A(2, 0) = normals[4](0);
      A(2, 1) = normals[4](1);
      A(2, 2) = normals[4](2);
      b(0) = -ds[0];
      b(1) = -ds[3];
      b(2) = -ds[4];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[3]] = x;

      /// node 4
      x.PutScalar(0.0);
      A(0, 0) = normals[5](0);
      A(0, 1) = normals[5](1);
      A(0, 2) = normals[5](2);
      A(1, 0) = normals[1](0);
      A(1, 1) = normals[1](1);
      A(1, 2) = normals[1](2);
      A(2, 0) = normals[4](0);
      A(2, 1) = normals[4](1);
      A(2, 2) = normals[4](2);
      b(0) = -ds[5];
      b(1) = -ds[1];
      b(2) = -ds[4];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[4]] = x;

      /// node 5
      x.PutScalar(0.0);
      A(0, 0) = normals[5](0);
      A(0, 1) = normals[5](1);
      A(0, 2) = normals[5](2);
      A(1, 0) = normals[1](0);
      A(1, 1) = normals[1](1);
      A(1, 2) = normals[1](2);
      A(2, 0) = normals[2](0);
      A(2, 1) = normals[2](1);
      A(2, 2) = normals[2](2);
      b(0) = -ds[5];
      b(1) = -ds[1];
      b(2) = -ds[2];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[5]] = x;

      /// node 6
      x.PutScalar(0.0);
      A(0, 0) = normals[5](0);
      A(0, 1) = normals[5](1);
      A(0, 2) = normals[5](2);
      A(1, 0) = normals[2](0);
      A(1, 1) = normals[2](1);
      A(1, 2) = normals[2](2);
      A(2, 0) = normals[3](0);
      A(2, 1) = normals[3](1);
      A(2, 2) = normals[3](2);
      b(0) = -ds[5];
      b(1) = -ds[2];
      b(2) = -ds[3];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[6]] = x;

      /// node 7
      x.PutScalar(0.0);
      A(0, 0) = normals[5](0);
      A(0, 1) = normals[5](1);
      A(0, 2) = normals[5](2);
      A(1, 0) = normals[3](0);
      A(1, 1) = normals[3](1);
      A(1, 2) = normals[3](2);
      A(2, 0) = normals[4](0);
      A(2, 1) = normals[4](1);
      A(2, 2) = normals[4](2);
      b(0) = -ds[5];
      b(1) = -ds[3];
      b(2) = -ds[4];
      A.Invert();
      x.Multiply(A, b);
      currentpositions[ele->NodeIds()[7]] = x;
    }
    break;
    case DRT::Element::tet4:
    {
      // surfaces are already planar
      for (int inode = 0; inode < ele->NumNode(); ++inode)
      {
        const DRT::Node* node = ele->Nodes()[inode];
        static LINALG::Matrix<3, 1> position;

        // fill position with node.X()
        for (int dim = 0; dim < 3; ++dim) position(dim) = node->X()[dim];
        currentpositions[ele->NodeIds()[inode]] = position;
      }
    }
    break;
    default:
      dserror("element type not yet supported");
      break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute superconvergent patch recovery for fluid frac    ghamm 02/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ComputePatchRecoveredFluidFraction(
    Teuchos::RCP<const Epetra_MultiVector> fluidfraction)
{
  Teuchos::ParameterList params;
  // only element center is relevant in this action
  params.set<int>("action", FLD::calc_velgrad_ele_center);

  Teuchos::RCP<Epetra_MultiVector> nodevec;
  if (BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_3D)
    nodevec = DRT::UTILS::ComputeSuperconvergentPatchRecovery<3>(
        fluiddis_, Teuchos::rcp((*fluidfraction)(0), false), "dummy", 1, params);
  else
    nodevec = DRT::UTILS::ComputeSuperconvergentPatchRecovery<2>(
        fluiddis_, Teuchos::rcp((*fluidfraction)(0), false), "dummy", 1, params);

  const int numnode = nodevec->MyLength();
  for (int i = 0; i < numnode; ++i)
  {
    const int nodeid = nodevec->Map().GID(i);

    // fill fluid fraction into pressure dof of dof row map --> always 3D problem!
    const int pressuredof = fluiddis_->Dof(fluiddis_->gNode(nodeid), 3);
    const int lid = fluidfracnp_->Map().LID(pressuredof);

    (*fluidfracnp_)[lid] = (*(*nodevec)(0))[i];
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute nodal fluid fraction based on L2 projection      ghamm 04/14 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ComputeL2ProjectedFluidFraction(
    Teuchos::RCP<const Epetra_MultiVector> fluidfraction)
{
  if (not fluidfraction->Map().SameAs(*fluiddis_->ElementRowMap()))
    dserror("input map is not an element row map");
  Teuchos::RCP<Epetra_Vector> eleveccol =
      Teuchos::rcp(new Epetra_Vector(*fluiddis_->ElementColMap(), true));
  LINALG::Export(*fluidfraction, *eleveccol);

  // handle pbcs if existing
  // build inverse map from slave to master nodes
  std::map<int, int> slavetomastercolnodesmap;
  {
    Teuchos::RCP<std::map<int, std::vector<int>>> allcoupledcolnodes =
        fluiddis_->GetAllPBCCoupledColNodes();

    for (std::map<int, std::vector<int>>::const_iterator masterslavepair =
             allcoupledcolnodes->begin();
         masterslavepair != allcoupledcolnodes->end(); ++masterslavepair)
    {
      // loop slave nodes associated with master
      for (std::vector<int>::const_iterator iter = masterslavepair->second.begin();
           iter != masterslavepair->second.end(); ++iter)
      {
        const int slavegid = *iter;
        slavetomastercolnodesmap[slavegid] = masterslavepair->first;
      }
    }
  }

  // get reduced node row map of fluid field --> will be used for setting up linear system
  const Epetra_Map* fullnoderowmap = fluiddis_->NodeRowMap();
  // remove pbc slave nodes from full noderowmap
  std::vector<int> reducednoderowmap;
  // a little more memory than necessary is possibly reserved here
  reducednoderowmap.reserve(fullnoderowmap->NumMyElements());
  for (int i = 0; i < fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);
    // do not add slave pbc nodes here
    if (slavetomastercolnodesmap.count(nodeid) == 0) reducednoderowmap.push_back(nodeid);
  }

  // build node row map which does not include slave pbc nodes
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(
      -1, (int)reducednoderowmap.size(), &reducednoderowmap[0], 0, fullnoderowmap->Comm()));

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> matrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*noderowmap, 108, false, true));
  // create empty right hand side
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*noderowmap, true);

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // get number of elements
  const int numele = fluiddis_->NumMyColElements();

  // loop column elements
  for (int i = 0; i < numele; ++i)
  {
    DRT::Element* actele = fluiddis_->lColElement(i);
    const int numnode = actele->NumNode();

    // get element location vector for nodes
    lm.resize(numnode);
    lmowner.resize(numnode);

    DRT::Node** nodes = actele->Nodes();
    for (int n = 0; n < numnode; ++n)
    {
      const int nodeid = nodes[n]->Id();

      std::map<int, int>::iterator slavemasterpair = slavetomastercolnodesmap.find(nodeid);
      if (slavemasterpair != slavetomastercolnodesmap.end())
        lm[n] = slavemasterpair->second;
      else
        lm[n] = nodeid;

      // owner of pbc master and slave nodes are identical
      lmowner[n] = nodes[n]->Owner();
    }

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(numnode);
    elematrix1.Shape(numnode, numnode);

    // set action in order to project element fluid fraction to nodal fluid fraction
    Teuchos::ParameterList params;
    params.set<int>("action", FLD::calc_fluidfrac_projection);

    // fill element fluid fraction value into params
    int eid = actele->Id();
    int lid = eleveccol->Map().LID(eid);
    if (lid < 0) dserror("element %i is not on this processor", eid);
    double eleval = (*eleveccol)[lid];
    params.set<double>("elefluidfrac", eleval);

    // call the element specific evaluate method (elevec1 = rhs, elemat1 = mass matrix)
    actele->Evaluate(
        params, *fluiddis_, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);

    // assembling into node maps
    matrix->Assemble(eid, elematrix1, lm, lmowner);
    LINALG::Assemble(*rhs, elevector1, lm, lmowner);
  }  // end element loop

  // finalize the matrix
  matrix->Complete();

  // get solver parameter list of linear solver
  const int solvernumber =
      DRT::Problem::Instance()->CavitationParams().get<int>("FLUIDFRAC_PROJ_SOLVER");
  // solve system
  Teuchos::RCP<Epetra_Vector> nodevec =
      SolveOnNodeBasedVector(solvernumber, noderowmap, matrix, rhs);

  // dofrowmap does not include separate dofs for slave pbc nodes
  for (int i = 0; i < noderowmap->NumMyElements(); ++i)
  {
    const int nodeid = noderowmap->GID(i);

    // fill fluid fraction into pressure dof of dof row map --> always 3D problem!
    const int pressuredof = fluiddis_->Dof(fluiddis_->gNode(nodeid), 3);
    const int lid = fluidfracnp_->Map().LID(pressuredof);

    (*fluidfracnp_)[lid] = (*nodevec)[i];
  }

  return;
}


/*----------------------------------------------------------------------*
 | solver: node based vector and full nullspace             ghamm 04/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CAVITATION::Algorithm::SolveOnNodeBasedVector(const int solvernumber,
    Teuchos::RCP<Epetra_Map> noderowmap, Teuchos::RCP<LINALG::SparseMatrix> matrix,
    Teuchos::RCP<Epetra_Vector> rhs)
{
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(solvernumber);

  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(
      solverparams, fluiddis_->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams, "AZPREC");
  switch (prectyp)
  {
    case INPAR::SOLVER::azprec_ML:
    case INPAR::SOLVER::azprec_MLfluid:
    case INPAR::SOLVER::azprec_MLAPI:
    case INPAR::SOLVER::azprec_MLfluid2:
    case INPAR::SOLVER::azprec_MueLuAMG_sym:
    case INPAR::SOLVER::azprec_MueLuAMG_nonsym:
    {
      Teuchos::ParameterList* preclist_ptr = NULL;
      // switch here between ML and MueLu cases
      if (prectyp == INPAR::SOLVER::azprec_ML or prectyp == INPAR::SOLVER::azprec_MLfluid or
          prectyp == INPAR::SOLVER::azprec_MLAPI or prectyp == INPAR::SOLVER::azprec_MLfluid2)
        preclist_ptr = &((solver->Params()).sublist("ML Parameters"));
      else if (prectyp == INPAR::SOLVER::azprec_MueLuAMG_sym or
               prectyp == INPAR::SOLVER::azprec_MueLuAMG_nonsym)
        preclist_ptr = &((solver->Params()).sublist("MueLu Parameters"));
      else
        dserror("please add correct parameter list");

      Teuchos::ParameterList& preclist = *preclist_ptr;
      preclist.set<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
      // ML would not tolerate this Teuchos::rcp-ptr in its list otherwise
      preclist.set<bool>("ML validate parameter list", false);

      preclist.set("PDE equations", 1);
      preclist.set("null space: dimension", 1);
      preclist.set("null space: type", "pre-computed");
      preclist.set("null space: add default vectors", false);

      // allocate the local length of the rowmap
      const int lrows = noderowmap->NumMyElements();
      Teuchos::RCP<std::vector<double>> ns = Teuchos::rcp(new std::vector<double>(lrows));
      double* nullsp = &((*ns)[0]);

      // compute null space manually
      for (int j = 0; j < lrows; ++j) nullsp[j] = 1.0;

      preclist.set<Teuchos::RCP<std::vector<double>>>("nullspace", ns);
      preclist.set("null space: vectors", nullsp);
    }
    break;
    case INPAR::SOLVER::azprec_ILU:
    case INPAR::SOLVER::azprec_ILUT:
      // do nothing
      break;
    default:
      dserror("You have to choose ML, MueLu or ILU preconditioning");
      break;
  }

  // solution vector based on node row map
  Teuchos::RCP<Epetra_Vector> nodevec = Teuchos::rcp(new Epetra_Vector(*noderowmap));

  // solve
  solver->Solve(matrix->EpetraOperator(), nodevec, rhs, true, true);

  return nodevec;
}


/*----------------------------------------------------------------------*
 | print cut situation for analytical void frac computation ghamm 08/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrintBubbleAndFluidEleToGMSH(DRT::Element* ele, const int bubbleid,
    const LINALG::Matrix<3, 1>& particleposition, const double influence)
{
  std::stringstream ss;
  ss << "bubble_" << bubbleid << "_element_" << ele->Id() << "_step";
  const std::string filename = IO::GMSH::GetFileName(ss.str(), Step(), false, Comm().MyPID());
  gmshfilecontent_.open(filename.c_str());

  // gmsh out of bubble
  {
    gmshfilecontent_ << "View \" "
                     << "bubble"
                     << " \" {\n";
    gmshfilecontent_ << "SH"
                     << "(";

    gmshfilecontent_ << particleposition(0) - influence << "," << particleposition(1) - influence
                     << "," << particleposition(2) - influence << ",";
    gmshfilecontent_ << particleposition(0) + influence << "," << particleposition(1) - influence
                     << "," << particleposition(2) - influence << ",";
    gmshfilecontent_ << particleposition(0) + influence << "," << particleposition(1) + influence
                     << "," << particleposition(2) - influence << ",";
    gmshfilecontent_ << particleposition(0) - influence << "," << particleposition(1) + influence
                     << "," << particleposition(2) - influence << ",";
    gmshfilecontent_ << particleposition(0) - influence << "," << particleposition(1) - influence
                     << "," << particleposition(2) + influence << ",";
    gmshfilecontent_ << particleposition(0) + influence << "," << particleposition(1) - influence
                     << "," << particleposition(2) + influence << ",";
    gmshfilecontent_ << particleposition(0) + influence << "," << particleposition(1) + influence
                     << "," << particleposition(2) + influence << ",";
    gmshfilecontent_ << particleposition(0) - influence << "," << particleposition(1) + influence
                     << "," << particleposition(2) + influence;

    gmshfilecontent_ << "){";
    for (int i = 0; i < ele->NumNode(); ++i)
    {
      if (i != 0) gmshfilecontent_ << ",";
      gmshfilecontent_ << ele->Id();
    }
    gmshfilecontent_ << "};\n";
    // end of current view
    gmshfilecontent_ << "};\n";
  }

  char elementtype;
  // gmsh out of surfaces of fluid ele
  {
    std::vector<Teuchos::RCP<DRT::Element>> surfaces = ele->Surfaces();
    const int numsurfacenodes = surfaces[0]->NumNode();

    gmshfilecontent_ << "View \" "
                     << "fluidele_surfaces"
                     << " \" {\n";
    for (int isurface = 0; isurface < ele->NumSurface(); ++isurface)
    {
      Teuchos::RCP<DRT::Element> surf = surfaces[isurface];

      switch (surf->NumNode())
      {
        case 3:
          elementtype = 'T';
          break;
        case 4:
          elementtype = 'Q';
          break;
        default:
          std::stringstream str;
          str << "unknown element type for " << ele->NumNode() << " nodes!";
          throw std::runtime_error(str.str());
      }

      gmshfilecontent_ << "S" << elementtype << "(";

      // loop over nodes of current surface
      for (int inode = 0; inode < numsurfacenodes; ++inode)
      {
        const DRT::Node* node = surf->Nodes()[inode];
        if (inode != 0) gmshfilecontent_ << ",";
        gmshfilecontent_ << node->X()[0] << "," << node->X()[1] << "," << node->X()[2];
      }

      gmshfilecontent_ << "){";
      for (int i = 0; i < surf->NumNode(); ++i)
      {
        if (i != 0) gmshfilecontent_ << ",";
        gmshfilecontent_ << isurface;
      }
      gmshfilecontent_ << "};\n";
    }
    // end of current view
    gmshfilecontent_ << "};\n";
  }

  return;
}


/*----------------------------------------------------------------------*
 | print integration lines for analytical void frac computation ghamm 08/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrintIntegrationLinesToGMSH(
    const std::vector<LINALG::Matrix<3, 1>>& integrationpoints, const int isurface)
{
  const int numintegrationpoints = (int)integrationpoints.size();
  gmshfilecontent_ << "View \" "
                   << "integrationline_surf_" << isurface << " \" {\n";
  for (int iring = 0; iring < numintegrationpoints; ++iring)
  {
    gmshfilecontent_ << "SL"
                     << "(";

    gmshfilecontent_ << integrationpoints[iring](0) << "," << integrationpoints[iring](1) << ","
                     << integrationpoints[iring](2) << ",";
    gmshfilecontent_ << integrationpoints[(iring + 1) % numintegrationpoints](0) << ","
                     << integrationpoints[(iring + 1) % numintegrationpoints](1) << ","
                     << integrationpoints[(iring + 1) % numintegrationpoints](2);
    gmshfilecontent_ << "){";

    gmshfilecontent_ << iring << "," << iring;
    gmshfilecontent_ << "};\n";
  }
  // end of current view
  gmshfilecontent_ << "};\n";

  return;
}
