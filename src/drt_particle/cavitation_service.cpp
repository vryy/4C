/*----------------------------------------------------------------------*/
/*!
\file cavitation_service.cpp

\brief Helper methods for cavitation simulations

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"


/*----------------------------------------------------------------------*
 | fluid fraction calculation                              ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateFluidFraction(
  Teuchos::RCP<const Epetra_Vector> particleradius
  )
{
  Teuchos::RCP<Epetra_FEVector> fluid_fraction = Teuchos::rcp(new Epetra_FEVector(*fluiddis_->ElementRowMap()));
  // export element volume to col layout
  Teuchos::RCP<Epetra_Vector> ele_volume_col = LINALG::CreateVector(*fluiddis_->ElementColMap(), false);
  LINALG::Export(*ele_volume_, *ele_volume_col);

  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->Dispnp();

  std::set<int> examinedbins;
  int numrownodes = particledis_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    DRT::Node *currentparticle = particledis_->lRowNode(i);
    DRT::Element** currentbin = currentparticle->Elements();

    const int binId=currentbin[0]->Id();

    if(examinedbins.count(binId) == 1)
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
    for(int iparticle=0; iparticle<currentbin[0]->NumNode(); ++iparticle)
    {
      DRT::Node* currparticle = particles[iparticle];
      double r_p = (*particleradius)[ particledis_->NodeRowMap()->LID(currparticle->Id()) ];
      maxradius = std::max(r_p, maxradius);
    }
    if(maxradius < 0.0)
      dserror("maximum radius smaller than zero");

    // get an ijk-range that is large enough
    int ijk[3];
    ConvertGidToijk(binId, ijk);

    // minimal bin size
    double minbin = std::numeric_limits<double>::max();
    for(int dim=0; dim<3; ++dim)
      minbin = std::min(minbin, bin_size_[dim]);

    // scaling factor in order to account for influence of bubble
    double scale = 1.2;
    const int ibinrange = (int)((maxradius*scale)/minbin) + 1;
    int ijk_range[] = {ijk[0]-ibinrange, ijk[0]+ibinrange, ijk[1]-ibinrange, ijk[1]+ibinrange, ijk[2]-ibinrange, ijk[2]+ibinrange};

    if(ibinrange > 3)
      dserror("not yet tested for such large bubbles");

    // variable to store bin ids of surrounding bins
    std::vector<int> binIds;
    binIds.reserve((2*ibinrange+1) * (2*ibinrange+1) * (2*ibinrange+1));

    // get corresponding bin ids in ijk range and fill them into binIds (in gid)
    GidsInijkRange(&ijk_range[0], binIds, false);

    // variable to store all fluid elements in neighborhood
    std::set<DRT::Element*> neighboringfluideles;
    for(std::vector<int>::const_iterator i=binIds.begin(); i!=binIds.end(); ++i)
    {
      // extract bins from discretization after checking on existence
      const int lid = particledis_->ElementColMap()->LID(*i);
      if(lid<0)
        continue;
      // extract bins from discretization
      DRT::MESHFREE::MeshfreeMultiBin* currbin =
          dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>( particledis_->lColElement(lid) );

      DRT::Element** currfluideles = currbin->AssociatedFluidEles();
      for(int ifluidele=0; ifluidele<currbin->NumAssociatedFluidEle(); ++ifluidele)
      {
        neighboringfluideles.insert(currfluideles[ifluidele]);
      }
    }

    // loop over all row particles in this bin and evaluate them
    for(int iparticle=0; iparticle<currentbin[0]->NumNode(); ++iparticle)
    {
      DRT::Node* currparticle = particles[iparticle];

      // fill particle position
      LINALG::Matrix<3,1> particleposition(false);
      std::vector<int> lm_b = particledis_->Dof(currparticle);
      // convert dof of particle into correct local id in this Epetra_Vector
      int posx = bubblepos->Map().LID(lm_b[0]);
      for (int dim=0; dim<3; ++dim)
      {
        // fill position of currentparticle
        particleposition(dim) = (*bubblepos)[posx+dim];
      }

      // get particle radius and influence of bubble
      const int lid = particleradius->Map().LID(currparticle->Id());
      if(lid<0)
        dserror("lid to gid %d not on this proc", currparticle->Id());
      const double r_p = (*particleradius)[lid];
      double influence = scale*r_p;

      // bubble volume which is assigned to fluid elements
      const double bubblevol = 4.0/3.0*M_PI*std::pow(r_p,3);

      // store all elements in which the bubble might be located fully inside
      std::vector<int> insideeles;

      // check for very small bubbles whose volume can be directly assigned to a
      // single fluid element
      {
        const int elelid = ele_volume_col->Map().LID((*neighboringfluideles.begin())->Id());
        if(elelid < 0)
          dserror("element id %i is not on this proc", (*neighboringfluideles.begin())->Id());
        const double charactelelength = std::pow((*ele_volume_col)[elelid], 1.0/3.0);
        // heuristic criterion of 1/20th of the characteristic element length of one fluid element in the neighborhood
        if (influence < 0.05*charactelelength)
        {
        // loop all surrounding fluid elements and find possible candidates for assigning
          for (std::set<DRT::Element*>::const_iterator ineighbor=neighboringfluideles.begin(); ineighbor!=neighboringfluideles.end(); ++ineighbor)
          {
            DRT::Element* ele = *ineighbor;
            const bool xaabboverlap = XAABBoverlap(ele, influence, particleposition);
            if(xaabboverlap == true)
              insideeles.push_back(ele->Id());
          }
          // assign to the correct fluid element
          AssignSmallBubbles(bubblevol, particleposition, insideeles, fluid_fraction);

          continue; // with next particle in this bin
        }
      }

      // process larger bubbles whose volume is distributed among the surrounding fluid elements
      {
        // variable to store not yet normalized volume for each fluid element
        std::map<int, double> volumefraction;
        double vol_influence = 0.0;
        bool surfaceoverlap = false;

        // do while volume is not negative which is important for overlapping surfaces/points of bubble and fluid
        do
        {
          // reset variables if necessary
          if(surfaceoverlap == true)
          {
            insideeles.clear();
            volumefraction.clear();
            vol_influence = 0.0;
            surfaceoverlap = false;
          }

          // loop over surrounding fluid elements
          for (std::set<DRT::Element*>::const_iterator ineighbor=neighboringfluideles.begin(); ineighbor!=neighboringfluideles.end(); ++ineighbor)
          {
            DRT::Element* ele = *ineighbor;

            const bool xaabboverlap = XAABBoverlap(ele, influence, particleposition);

            double vol_ele = 0.0;
            if(xaabboverlap == true)
            {
              // store element id
              insideeles.push_back(ele->Id());
              switch(void_frac_strategy_)
              {
              case INPAR::CAVITATION::analytical_constpoly:
              case INPAR::CAVITATION::analytical_quadraticpoly:
              {
                DoAnalyticalIntegrationFluidFrac(ele, particleposition, influence, vol_ele, surfaceoverlap);

                // in case of polynomial influence, twisted surfaces and tight cut situations,
                // the polynomial can become negative very close to the bubble influence area
                if(vol_ele < 0.0)
                  vol_ele = 0.0;
              }
              break;
              case INPAR::CAVITATION::gaussian_integration:
              {
                DoGaussianIntegrationFluidFrac(ele, particleposition, influence, vol_ele);
              }
              break;
              default:
                dserror("void fraction calculation strategy does not exist");
              break;
              }

              // in case of tight cuts, rerun neighboring ele loop with slightly smaller influence
              if(surfaceoverlap == true)
              {
                influence *= 0.999;
                break;
              }

              // safety check
              if(vol_ele < 0.0)
                dserror("negative volume occured during void fraction computation");

              // sum and store volume for each fluid element
              vol_influence += vol_ele;
              volumefraction[ele->Id()] = vol_ele;
            } // end if xaabboverlap
          } // end loop neighboring eles

        } while(surfaceoverlap == true);

        // distribute volumes stored in volumefraction
        if(vol_influence != 0.0)
        {
          // normalize with volume of bubble divided by volume of influence
          const double normalization = bubblevol / vol_influence;
          // assemble void volume of fluid elements
          for(std::map<int, double>::const_iterator iter=volumefraction.begin(); iter!=volumefraction.end(); ++iter)
          {
            const double val = iter->second * normalization;

            // do assembly of void fraction into fluid
            int err = fluid_fraction->SumIntoGlobalValues(1, &(iter->first), &val);
            if (err<0)
              dserror("summing into Epetra_FEVector failed");
          }
        }
        else
        {
          AssignSmallBubbles(bubblevol, particleposition, insideeles, fluid_fraction);
        }
      }

    } // loop over particles in one bin
  } // loop over outer particles

  // call global assemble
  int err = fluid_fraction->GlobalAssemble(Add, false);
  if (err<0)
    dserror("global assemble into fluidforces failed");

  // compute fluid fraction in elements (= ele_volume - void fraction)
  fluid_fraction->Update(1.0, *ele_volume_, -1.0);

  // divide element wise fluid volume by element volume
  fluid_fraction->ReciprocalMultiply(1.0, *ele_volume_, *fluid_fraction, 0.0);

  double minfluidfrac = 0.0;
  fluid_fraction->MinValue(&minfluidfrac);
  if(minfluidfrac < 0.05 && fluid_->PhysicalType() == INPAR::FLUID::loma)
    dserror("minimum fluid fraction is low (%f), increase mesh size or reduce particle size", minfluidfrac);

  double maxfluidfrac = 0.0;
  fluid_fraction->MaxValue(&maxfluidfrac);
  if(myrank_ == 0)
    std::cout << std::setprecision(16) << "minimum fluid fraction is: " << minfluidfrac << " and maximimum fluid fraction is: " << maxfluidfrac << std::endl;

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
void CAVITATION::Algorithm::DoGaussianIntegrationFluidFrac(
  DRT::Element* ele,
  LINALG::Matrix<3,1>& particleposition,
  const double influence,
  double& vol_ele
  )
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
  ele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

  // Reshape element matrices and vectors and initialize to zero
  elevector1.Size(1);

  // set action in order to calculate fluid fraction in this element
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::calc_volume_gaussint);
  params.set<LINALG::Matrix<3,1> >("particlepos", particleposition);
  params.set<double>("influence", influence);

  params.set<int>("gp_per_dir", gauss_rule_per_dir_);

  // call the element specific evaluate method (elevec1 = void volume)
  ele->Evaluate(params,*fluiddis_,lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

  vol_ele = elevector1[0];

  return;
}


/*----------------------------------------------------------------------*
 | analytic integration for void volume calculation        ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::DoAnalyticalIntegrationFluidFrac(
  DRT::Element* ele,
  LINALG::Matrix<3,1>& particleposition,
  const double influence,
  double& vol_ele,
  bool& surfaceoverlap
  )
{
  // pull out each volume element and ask for surfaces
  std::vector<Teuchos::RCP<DRT::Element> > surfaces = ele->Surfaces();
  int numsurfacenodes = surfaces[0]->NumNode();

  // center of fluid element
  LINALG::Matrix<3,1> centerele(true);
  const DRT::Node* nodeele;
  for(int inodeele=0; inodeele<ele->NumNode(); ++inodeele)
  {
    nodeele = ele->Nodes()[inodeele];
    for(int dim=0; dim<3; ++dim)
    {
      centerele(dim) = centerele(dim)+nodeele->X()[dim];
    }
  }
  centerele.Scale(1.0/ele->NumNode());

  // variables for calculating bubblesurfaces
  std::vector<LINALG::Matrix<3,1> > bubbleX, bubble_X;
  bubbleX.reserve(numsurfacenodes+2);
  bubble_X.reserve(numsurfacenodes+2);

  // loop over surfaces of current fluid element and feed variable elements and currentpositions
  for (int isurface=0; isurface<ele->NumSurface(); ++isurface)
  {
    // fill variables and corner points of surface
    std::map<int, Teuchos::RCP<DRT::Element> > elements;
    std::map<int, LINALG::Matrix<3,1> > currentpositions;

    std::vector<LINALG::Matrix<3,1> > surfacenodes(numsurfacenodes);

    elements[surfaces[isurface]->Id()] = surfaces[isurface];

    // loop over nodes of current surface: surfaces[isurface]
    for(int inode=0; inode<numsurfacenodes; ++inode)
    {
      const DRT::Node* node = surfaces[isurface]->Nodes()[inode];
      LINALG::Matrix<3,1> position;

      // fill position with node.X()
      for (int dim=0; dim<3; ++dim)
      {
        position(dim) = node->X()[dim];
      }
      surfacenodes[inode] = position;
      currentpositions[node->Id()] = position;
    }

    // get XAABB of this single surface element
    LINALG::Matrix<3,2> xaabb = GEO::getXAABBofEles(elements, currentpositions);
    bool boundingbox = true;

    LINALG::Matrix<3,1> n, r, p;
    r.Update(1.0,surfacenodes[0],-1.0,surfacenodes[2]);
    p.Update(1.0,surfacenodes[1],-1.0,surfacenodes[3%numsurfacenodes]);

    //calculate normalvector n with cross product of diagonals r,p
    n(0) = r(1)*p(2)-r(2)*p(1);
    n(1) = r(2)*p(0)-r(0)*p(2);
    n(2) = r(0)*p(1)-r(1)*p(0);
    n.Scale(1.0/n.Norm2());

    // bubblesurfaces of the cubic bubble
    double bubblesurface[] = {particleposition(0)+influence, particleposition(0)-influence,
                              particleposition(1)+influence, particleposition(1)-influence,
                              particleposition(2)+influence, particleposition(2)-influence};

    // test the bounding box touching the bubbleinfluence
    for(int dim=0; dim<3; ++dim)
      if(xaabb(dim,0) > bubblesurface[dim*2] or xaabb(dim,1) < bubblesurface[dim*2+1])
      {
        boundingbox = false;
        break;
      }

    if(boundingbox == true)
    {
      // is there divergence in x-direction?
      if (n(0)<-0.00001 or n(0)>0.00001)
      {
        EvaluateSurface(surfacenodes, n, centerele , particleposition, influence, vol_ele, surfaceoverlap);
      }
    }

    //store penetration points of +x and -x surfaces of the bubble
    GetPenetrationPointsOfXSurfaces(n, surfacenodes, influence, particleposition, bubbleX, bubble_X, surfaceoverlap);
  }

  // integration over x-oriented bubble surface
  if(bubbleX.size() != 0)
  {
    LINALG::Matrix<3,1> n;
    n(0) = 1.0;
    n(1) = 0.0;
    n(2) = 0.0;

    if (bubbleX.size() > 3)
    {
      BuildConvexHull(bubbleX);
    }
    EvaluateSurface(bubbleX, n, particleposition , particleposition, influence, vol_ele, surfaceoverlap);
  }

  // integration over -x-oriented bubble surface
  if(bubble_X.size() != 0)
  {
    LINALG::Matrix<3,1> n;
    n(0) = -1.0;
    n(1) = 0.0;
    n(2) = 0.0;

    if(bubble_X.size() > 3)
    {
      BuildConvexHull(bubble_X);
    }
    EvaluateSurface(bubble_X, n, particleposition , particleposition, influence, vol_ele, surfaceoverlap);
  }

return;
}


/*----------------------------------------------------------------------*
 | compute integration points for analytic integration     ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::EvaluateSurface(
  std::vector<LINALG::Matrix<3,1> >& surfacenodes,
  const LINALG::Matrix<3,1>& n,
  const LINALG::Matrix<3,1>& centerele,
  LINALG::Matrix<3,1>& particleposition,
  const double influence,
  double& vol_ele,
  bool& surfaceoverlap
  )
{
  //**********************************************************************
  // Evaluating surfaces over the bubble
  // - calculating corner points
  // - calculating penetration points
  // - calculating bubblecorner points
  //**********************************************************************

  // variable to store integrationpoints
  std::vector<LINALG::Matrix<3,1> > integrationpoints;
  integrationpoints.reserve(12);
  double bubblesurface[]={particleposition(0)+influence , particleposition(0)-influence ,
                          particleposition(1)+influence , particleposition(1)-influence ,
                          particleposition(2)+influence , particleposition(2)-influence};

  double t;
  double tol=1.0e-8;

  int numsurfacenodes = (int)surfacenodes.size();
  LINALG::Matrix<3,1> centersurf(true);
  for(int inode=0; inode<numsurfacenodes; ++inode)
  {
    centersurf.Update(1.0, surfacenodes[inode], 1.0);
  }
  centersurf.Scale(1.0/numsurfacenodes);

  // n is wrong oriented
  if ((centersurf(0)-centerele(0))*n(0)+(centersurf(1)-centerele(1))*n(1)+(centersurf(2)-centerele(2))*n(2)<0)
    dserror("you should not show up here");

  LINALG::Matrix<3,1> currpoint, difference;
  // loop over lines to find corner points within the influence of the bubble and penetration points which penetrate the bubble surface
  for(int iedges=0; iedges<numsurfacenodes; ++iedges)
  {
    difference.Update(1,surfacenodes[iedges],-1.0,particleposition);

    // catch corner points cp by testing whether NormInf of the current corner point is smaller than the influenceradius
    // corner points are edge-points of the fluid surface, which are into the bubble influence area
    // InfNorm due to a cubic bubble
    if(difference.NormInf() <= influence+tol)
    {
      integrationpoints.push_back(surfacenodes[iedges]);
    }

    // catch penetration points pp which penetrate the bubble surface by building a line between two surfacepoints
    // penetration points are at the edge lines of the surface and penetrate the bubble surface
    // and test whether the penetration happens between them with 0<t<1 and the point is in InfNorm influence of the bubble
    for(int ipene=0; ipene<6; ++ipene)
    {
      // parameter t for penetrating bubble surface
      if((surfacenodes[(iedges+1)%numsurfacenodes](ipene/2)-surfacenodes[iedges](ipene/2))>tol or (surfacenodes[(iedges+1)%numsurfacenodes](ipene/2)-surfacenodes[iedges](ipene/2)<-tol))
        {
          t = (bubblesurface[ipene]-surfacenodes[iedges](ipene/2))/(surfacenodes[(iedges+1)%numsurfacenodes](ipene/2)-surfacenodes[iedges](ipene/2));
        }
      else
        t = -1.0;

      // current bubble surface is penetrated (first condition)
      if(+tol<t and t<1-tol)
      {
        currpoint.Update(t,surfacenodes[(iedges+1)%numsurfacenodes],1.0-t,surfacenodes[iedges]);
        difference.Update(1.0, currpoint, -1.0, particleposition);

        // second condition
        // +tol because penetration points are always in the influence distance in NormInf
        if(difference.NormInf() <= influence+tol)
        {
          integrationpoints.push_back(currpoint);
        }
      }
    }
  }

  // get bubblecorner points which are in the surface
  // calculates the bubblecorner intersection point of the two surfaces
  // A bubblecorner point is an intersection point of the fluid surface with the bubbleedges
  // calculates the intersection points of the bubbleedges with the plane surface (represented with n,centersurf)
  for(int y=2; y<4; y++)
    for(int z=4; z<6; z++)
    {
      LINALG::Matrix<3,1> bubblecorner;

      CalculateBubbleCornerPoint(n, centersurf, y, z, bubblesurface, bubblecorner);
      bool inpoint = CheckPointInSurface(surfacenodes, centersurf, centerele, bubblecorner);
      difference.Update(1,bubblecorner,-1,particleposition);
      if(difference.NormInf()<=influence+tol and inpoint==true)
      {
        integrationpoints.push_back(bubblecorner);
      }
    }
  // if normal vector in y-direction == 0: line of bubble is parallel to the fluid surface and can be omitted
  if(n(1)>tol or n(1)<-tol)
    for(int x=0; x<2; x++)
      for(int z=4; z<6; z++)
      {
        LINALG::Matrix<3,1> bubblecorner;

        CalculateBubbleCornerPoint(n, centersurf, x, z, bubblesurface, bubblecorner);
        bool inpoint = CheckPointInSurface(surfacenodes, centersurf, centerele, bubblecorner);
        difference.Update(1.0, bubblecorner, -1.0, particleposition);
        if(difference.NormInf()<=influence+tol and inpoint==true)
        {
          integrationpoints.push_back(bubblecorner);
        }
      }
  // if normal vector in z-direction == 0: line of bubble is parallel to the fluid surface and can be omitted
  if(n(2)>tol or n(2)<-tol)
    for(int x=0; x<2; x++)
      for(int y=2; y<4; y++)
      {
        LINALG::Matrix<3,1> bubblecorner;

        CalculateBubbleCornerPoint(n, centersurf, x, y, bubblesurface, bubblecorner);
        bool inpoint = CheckPointInSurface(surfacenodes, centersurf, centerele, bubblecorner);
        difference.Update(1.0, bubblecorner,-1.0, particleposition);
        if(difference.NormInf()<=influence+tol and inpoint==true)
        {
          integrationpoints.push_back(bubblecorner);
        }
      }

  // calculate the ring integral over fluid surfaces
  if(integrationpoints.size() > 2)
  {
    // integration points are ordered for correct ring integral evaluation
    // close points are removed altering integrationpoints.size()
    if(integrationpoints.size() > 3)
      BuildConvexHull(integrationpoints);

    int numintegrationpoints = (int)integrationpoints.size();

    LINALG::Matrix<3,1> centerringintgral(true);
    for(int inode=0; inode<numintegrationpoints; ++inode)
    {
      centerringintgral.Update(1.0, integrationpoints[inode], 1.0);
    }
    centerringintgral.Scale(1.0/numintegrationpoints);

    for(int iring=0; iring<numintegrationpoints; ++iring)
    {
      switch(void_frac_strategy_)
      {
      case INPAR::CAVITATION::analytical_constpoly:
      {
        EvaluateTwoPointsConstPoly(n, centerringintgral, integrationpoints[iring], integrationpoints[(iring+1)%numintegrationpoints], vol_ele);
      }
      break;
      case INPAR::CAVITATION::analytical_quadraticpoly:
      {
        EvaluateTwoPointsQuadraticPoly(n, centerringintgral, integrationpoints[iring], integrationpoints[(iring+1)%numintegrationpoints], particleposition, vol_ele, influence);
      }
      break;
      default:
        dserror("volume strategy does not exist");
      break;
      }
    }
  }
  else if(integrationpoints.size() > 0)
  {
    // surface overlap detected
    surfaceoverlap = true;
  }

  return;
}


/*----------------------------------------------------------------------*
 | find penetration points of bubble surfaces (+x and -x)  ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::GetPenetrationPointsOfXSurfaces(
  const LINALG::Matrix<3,1>& n,
  std::vector<LINALG::Matrix<3,1> >& surfacenodes,
  const double influence,
  const LINALG::Matrix<3,1>& particleposition,
  std::vector<LINALG::Matrix<3,1> >& bubbleX,
  std::vector<LINALG::Matrix<3,1> >& bubble_X,
  bool& surfaceoverlap)
{
  // method counts penetration points of one element and stores them
  // bubbleX array of all penetration points
  // bubble_X array of all penetration points

  double t;
  double tol = 1.0e-8;
  double bubblesurface[] = {particleposition(0)+influence , particleposition(0)-influence};
  int numsurfacenodes = (int)surfacenodes.size();

  // store penetration points of +x and -x-surfaces
  // penetration points are at the edge lines of the surface and penetrate the bubble surface
  // and test whether the penetration happens between them with 0<t<1 and the point is in NormInf influence of the bubble
  for(int iedges=0; iedges<numsurfacenodes; ++iedges)
  {
    for(int ipene=0; ipene<2; ++ipene)
    {
      //parameter t for penetrating bubble surface
      if((surfacenodes[(iedges+1)%numsurfacenodes](0)-surfacenodes[iedges](0))>tol or (surfacenodes[(iedges+1)%numsurfacenodes](0)-surfacenodes[iedges](0)<-tol))
        {
          t = (bubblesurface[ipene]-surfacenodes[iedges](0))/(surfacenodes[(iedges+1)%numsurfacenodes](0)-surfacenodes[iedges](0));
          if((-tol*10.0 < t and t < tol*10.0) or (1.0-tol*10.0 < t and t < 1.0+tol*10.0))
          {
            // due to surface overlap, do integration again with smaller influence radius
            surfaceoverlap = true;
          }
        }
      else
        t = -1.0;

      // current bubble surface is penetrated (first condition)
      if(tol<t and t<1.0-tol)
      {
        LINALG::Matrix <3,1> difference, pp;
        pp.Update(t,surfacenodes[(iedges+1)%numsurfacenodes],1-t,surfacenodes[iedges]);

        // surface which  is only described with n(0)=1 does not penetrate
        if(n(0) < 1.0)
        {
          if(ipene == 0) // +X bubble surface
          {
            if(bubbleX.size() == 0)
            {
              bubbleX.push_back(pp);
            }
            else
            {
              bool alreadyin = false;
              for(size_t i=0; i<bubbleX.size(); ++i)
              {
                difference.Update(1.0, pp, -1.0, bubbleX[i]);
                if(difference.Norm2()<tol)
                {
                  alreadyin = true;
                  break;
                }
              }
              if(alreadyin == false)
              {
                bubbleX.push_back(pp);
              }
            }
          }
          else if(ipene == 1) // -X bubble surface
          {
            if(bubble_X.size()==0)
            {
              bubble_X.push_back(pp);
            }
            else
            {
              bool alreadyin = false;
              for(size_t i=0; i<bubble_X.size(); ++i)
              {
                difference.Update(1.0, pp, -1.0, bubble_X[i]);
                if(difference.Norm2()<tol)
                {
                  alreadyin = true;
                  break;
                }
              }
              if(alreadyin == false)
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
void CAVITATION::Algorithm::CalculateBubbleCornerPoint(
  const LINALG::Matrix<3,1>& n,
  const LINALG::Matrix<3,1>& centersurface,
  int bubblesurface1,
  int bubblesurface2,
  double *bubblesurface,
  LINALG::Matrix<3,1>& bubblecorner)
{
  // calculates the bubblecorner intersection point of the two surfaces
  // A bubblecorner point is a intersection point of the fluid surface with the bubbleedges
  // suface 0 is +x
  // 1 is -x
  // 2 is y
  // 3 is -y
  // 4 is z
  // 5 is -z
  int iset[3];
  iset[0] = 0;
  iset[1] = 0;
  iset[2] = 0;
  for(int j=0; j<6; j++)
  {
    if(bubblesurface1 == j)
    {
      bubblecorner(j/2) = bubblesurface[j];
      iset[j/2] = 1;
    }
    else if(bubblesurface2 == j)
    {
      bubblecorner(j/2) = bubblesurface[j];
      iset[j/2] = 1;
    }
  }
  // x direction variable
  if(iset[0] == 0)
  {
    bubblecorner(0) = (n(0)*centersurface(0)+n(1)*centersurface(1)+n(2)*centersurface(2)-n(1)*bubblecorner(1)-n(2)*bubblecorner(2))/n(0);
  }
  // y direction variable
  if(iset[1] == 0)
  {
    if(n(1)>0.0001 or n(1)<-0.0001)
    {
      bubblecorner(1) = (n(0)*centersurface(0)+n(1)*centersurface(1)+n(2)*centersurface(2)-n(0)*bubblecorner(0)-n(2)*bubblecorner(2))/n(1);
    }
    else
    {
      bubblecorner(1) = centersurface(1);
    }
  }
  // z direction variable
  if(iset[2] == 0)
  {
    if(n(2)>0.0001 or n(2)<-0.0001)
    {
      bubblecorner(2) = (n(0)*centersurface(0)+n(1)*centersurface(1)+n(2)*centersurface(2)-n(1)*bubblecorner(1)-n(0)*bubblecorner(0))/n(2);
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
bool CAVITATION::Algorithm::CheckPointInSurface(
  std::vector<LINALG::Matrix<3,1> >& surfacenodes,
  const LINALG::Matrix<3,1>&  centersurface,
  const LINALG::Matrix<3,1>&  centerele,
  const LINALG::Matrix<3,1>&  pointtocheck
  )
{
  // tests whether a point is in the element surface by calculating the cross product and
  // check whether the point is always left or right of the surrounding lines
  bool inpoint = false;
  int leftright = 0;

  int numsurfacenodes = (int)surfacenodes.size();
  for(int iedges=0; iedges<numsurfacenodes; ++iedges)
  {
    double DY = surfacenodes[(iedges+1)%numsurfacenodes](1)-surfacenodes[iedges](1);
    double DZ = surfacenodes[(iedges+1)%numsurfacenodes](2)-surfacenodes[iedges](2);
    double BY = pointtocheck(1)-surfacenodes[iedges](1);
    double BZ = pointtocheck(2)-surfacenodes[iedges](2);
    if((DY*BZ-DZ*BY)*(centersurface(0)-centerele(0))>0)
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
  if(leftright==0 or leftright==numsurfacenodes)
    inpoint = true;

  return inpoint;
}


/*----------------------------------------------------------------------*
 | build convex hull of points in y-z plane                ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildConvexHull(std::vector<LINALG::Matrix<3,1> >& surfacenodes)
{
  int np = (int)surfacenodes.size();
  bool out = false;
  double tol = 1.0e-8;
  std::vector<MORTAR::Vertex> respoly;
  std::vector<MORTAR::Vertex> collconvexhull;

  // temporary storage for transformed points
  Epetra_SerialDenseMatrix transformed(2,np);

  std::vector<int> dummy;
  std::vector<double> coords(3);
  // transform each convex hull point
  for (int i=0; i<np; ++i)
  {
    for(int k=0; k<3; ++k)
      coords[k] = surfacenodes[i](k);

    collconvexhull.push_back(MORTAR::Vertex(coords,MORTAR::Vertex::lineclip,dummy,NULL,NULL,false,false,NULL,-1.0));
    // x coordinate doesn't matter. projection into y-z layer
    for (int k=1; k<3; ++k)
    {
      transformed(k-1,i) = coords[k];
    }
  }

  // sort convex hull points to obtain final ring integral
  int removed = MORTAR::SortConvexHullPoints(out, transformed, collconvexhull, respoly, tol);

  if(removed>0)
    np=np-removed;

  // fill in in new order
  surfacenodes.resize(np);
  for(int i=0; i<np; ++i)
  {
    for(int k=0; k<3; ++k)
    {
      surfacenodes[i](k) = respoly[i].Coord()[k];
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | evaluate line integral for analytic integration         ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::EvaluateTwoPointsConstPoly(
  const LINALG::Matrix<3,1>&  n,
  const LINALG::Matrix<3,1>&  centerringintgral,
  const LINALG::Matrix<3,1>&  point1,
  const LINALG::Matrix<3,1>&  point2,
  double& vol_ele
  )
{
  // calculate the line integral of bubble in the fluid element
  // direction of integration automatically corrected
  LINALG::Matrix<3,1>  m, delta;
  double Y,dY,Z,dZ,a,b,c;

  delta.Update(1.0, point2, -1.0, point1);
  m(0) = 0.0;
  if(n(0) > 0)
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
  a = -0.5*n(1)/n(0);
  b = centerringintgral(0)+(centerringintgral(1)*n(1)+centerringintgral(2)*n(2))/n(0);
  c = -n(2)/n(0);

  // declare start of integration
  if((m(1)*((point1(1)+point2(1))/2-centerringintgral(1))+m(2)*((point1(2)+point2(2))/2-centerringintgral(2))) > 0)
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
  vol_ele += dZ*(a*(Y*Y+Y*dY+dY*dY/3)+b*(Y+0.5*dY)+c*(Z*(Y+dY*0.5)+dZ*(Y*0.5+dY/3)));

  return;
}


/*----------------------------------------------------------------------*
 | evaluate line integral for analytic integration         ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::EvaluateTwoPointsQuadraticPoly(
  const LINALG::Matrix<3,1>& n,
  const LINALG::Matrix<3,1>& centerringintgral,
  LINALG::Matrix<3,1>& point1,
  LINALG::Matrix<3,1>& point2,
  const LINALG::Matrix<3,1>& particleposition,
  double& vol_ele,
  const double influence
  )
{
  // calculate the line integral of bubble in the fluid element
  // direction of integration automatically corrected
  LINALG::Matrix<3,1> m, delta, point1CoTr, point2CoTr, middleofsurface;
  double a,b,c,d,e,f,g,s;

  // coordinate system transformation
  point1CoTr.Update(1.0, point1, -1.0, particleposition);
  point2CoTr.Update(1.0, point2, -1.0, particleposition);
  delta.Update(1.0, point2CoTr, -1.0, point1CoTr);
  middleofsurface.Update(1.0, centerringintgral, -1.0, particleposition);

  m(0) = 0;
  m(1) = delta(2)*n(0);
  m(2) = delta(1)*-n(0);

  a = -n(1)/n(0);
  b = -n(2)/n(0);
  c = middleofsurface(0)+middleofsurface(1)*n(1)/n(0)+middleofsurface(2)*n(2)/n(0);

  // check direction of integration
  if((m(1)*(point1CoTr(1)-middleofsurface(1))+m(2)*(point1CoTr(2)-middleofsurface(2))) > 0)
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
  if(g > 1.0e-8 or g < -1.0e-8 )
  {
    // add another line integral
    s = 27.0/(8.0*pow(influence,6.0));
    vol_ele +=
        s*g*(((-20*pow(c,3)*pow(d,3)*pow(f,2) - 45*a*pow(c,2)*pow(d,4)*pow(f,2) - 36*pow(a,2)*c*pow(d,5)*pow(f,2) - 10*pow(a,3)*pow(d,6)*pow(f,2) - 60*b*pow(c,2)*pow(d,3)*pow(f,3) - 90*a*b*c*pow(d,4)*pow(f,3) - 36*pow(a,2)*b*pow(d,5)*pow(f,3) - 60*pow(b,2)*c*pow(d,3)*pow(f,4) - 45*a*pow(b,2)*pow(d,4)*pow(f,4) - 20*pow(b,3)*pow(d,3)*pow(f,5) + 20*pow(c,3)*pow(d,3)*pow(influence,2) + 45*a*pow(c,2)*pow(d,4)*pow(influence,2) + 36*pow(a,2)*c*pow(d,5)*pow(influence,2) + 10*pow(a,3)*pow(d,6)*pow(influence,2) + 60*b*pow(c,2)*pow(d,3)*f*pow(influence,2) + 90*a*b*c*pow(d,4)*f*pow(influence,2) + 36*pow(a,2)*b*pow(d,5)*f*pow(influence,2) + 60*pow(c,3)*d*pow(f,2)*pow(influence,2) + 90*a*pow(c,2)*pow(d,2)*pow(f,2)*pow(influence,2) + 60*(1 + pow(a,2))*c*pow(d,3)*pow(f,2)* pow(influence,2) + 60*pow(b,2)*c*pow(d,3)*pow(f,2)*pow(influence,2) + 15*a*(3 + pow(a,2))*pow(d,4)* pow(f,2)*pow(influence,2) + 45*a*pow(b,2)*pow(d,4)*pow(f,2)*pow(influence,2) + 180*b*pow(c,2)*d*pow(f,3)*pow(influence,2) + 180*a*b*c*pow(d,2)*pow(f,3)*pow(influence,2) + 60*(1 + pow(a,2))*b*pow(d,3)*pow(f,3)*pow(influence,2) + 20*pow(b,3)*pow(d,3)*pow(f,3)* pow(influence,2) + 180*pow(b,2)*c*d*pow(f,4)*pow(influence,2) + 90*a*pow(b,2)*pow(d,2)*pow(f,4)* pow(influence,2) + 60*pow(b,3)*d*pow(f,5)*pow(influence,2) - 60*pow(c,3)*d*pow(influence,4) - 90*a*pow(c,2)*pow(d,2)*pow(influence,4) - 60*(1 + pow(a,2))*c*pow(d,3)*pow(influence,4) - 15*a*(3 + pow(a,2))*pow(d,4)*pow(influence,4) - 180*b*pow(c,2)*d*f*pow(influence,4) - 180*a*b*c*pow(d,2)*f*pow(influence,4) - 60*(1 + pow(a,2))*b*pow(d,3)*f*pow(influence,4) - 180*c*d*pow(f,2)*pow(influence,4) - 180*pow(b,2)*c*d*pow(f,2)*pow(influence,4) - 90*a*pow(d,2)*pow(f,2)*pow(influence,4) - 90*a*pow(b,2)*pow(d,2)*pow(f,2)*pow(influence,4) - 180*b*d*pow(f,3)*pow(influence,4) - 60*pow(b,3)*d*pow(f,3)*pow(influence,4) + 180*c*d*pow(influence,6) + 90*a*pow(d,2)*pow(influence,6) + 180*b*d*f*pow(influence,6)) + (-30*pow(c,3)*pow(d,2)*e*pow(f,2) - 90*a*pow(c,2)*pow(d,3)*e*pow(f,2) - 90*pow(a,2)*c*pow(d,4)*e*pow(f,2) - 30*pow(a,3)*pow(d,5)*e*pow(f,2) - 90*b*pow(c,2)*pow(d,2)*e*pow(f,3) - 180*a*b*c*pow(d,3)*e*pow(f,3) - 90*pow(a,2)*b*pow(d,4)*e*pow(f,3) - 90*pow(b,2)*c*pow(d,2)*e*pow(f,4) - 90*a*pow(b,2)*pow(d,3)*e*pow(f,4) - 30*pow(b,3)*pow(d,2)*e*pow(f,5) - 20*pow(c,3)*pow(d,3)*f*g - 45*a*pow(c,2)*pow(d,4)*f*g - 36*pow(a,2)*c*pow(d,5)*f*g - 10*pow(a,3)*pow(d,6)*f*g - 90*b*pow(c,2)*pow(d,3)*pow(f,2)*g - 135*a*b*c*pow(d,4)*pow(f,2)*g - 54*pow(a,2)*b*pow(d,5)*pow(f,2)*g - 120*pow(b,2)*c*pow(d,3)*pow(f,3)*g - 90*a*pow(b,2)*pow(d,4)*pow(f,3)*g - 50*pow(b,3)*pow(d,3)*pow(f,4)*g + 30*pow(c,3)*pow(d,2)*e*pow(influence,2) + 90*a*pow(c,2)*pow(d,3)*e*pow(influence,2) + 90*pow(a,2)*c*pow(d,4)*e*pow(influence,2) + 30*pow(a,3)*pow(d,5)*e*pow(influence,2) + 90*b*pow(c,2)*pow(d,2)*e*f*pow(influence,2) + 180*a*b*c*pow(d,3)*e*f*pow(influence,2) + 90*pow(a,2)*b*pow(d,4)*e*f*pow(influence,2) + 30*pow(c,3)*e*pow(f,2)*pow(influence,2) + 90*a*pow(c,2)*d*e*pow(f,2)*pow(influence,2) + 90*c*pow(d,2)*e*pow(f,2)*pow(influence,2) + 90*pow(a,2)*c*pow(d,2)*e*pow(f,2)*pow(influence,2) + 90*pow(b,2)*c*pow(d,2)*e*pow(f,2)* pow(influence,2) + 90*a*pow(d,3)*e*pow(f,2)*pow(influence,2) + 30*pow(a,3)*pow(d,3)*e*pow(f,2)* pow(influence,2) + 90*a*pow(b,2)*pow(d,3)*e*pow(f,2)*pow(influence,2) + 90*b*pow(c,2)*e*pow(f,3)*pow(influence,2) + 180*a*b*c*d*e*pow(f,3)*pow(influence,2) + 90*b*pow(d,2)*e*pow(f,3)*pow(influence,2) + 90*pow(a,2)*b*pow(d,2)*e*pow(f,3)*pow(influence,2) + 30*pow(b,3)*pow(d,2)*e*pow(f,3)*pow(influence,2) + 90*pow(b,2)*c*e*pow(f,4)*pow(influence,2) + 90*a*pow(b,2)*d*e*pow(f,4)*pow(influence,2) + 30*pow(b,3)*e*pow(f,5)*pow(influence,2) + 30*b*pow(c,2)*pow(d,3)*g*pow(influence,2) + 45*a*b*c*pow(d,4)*g*pow(influence,2) + 18*pow(a,2)*b*pow(d,5)*g*pow(influence,2) + 60*pow(c,3)*d*f*g*pow(influence,2) + 90*a*pow(c,2)*pow(d,2)*f*g*pow(influence,2) + 60*c*pow(d,3)*f*g*pow(influence,2) + 60*pow(a,2)*c*pow(d,3)*f*g*pow(influence,2) + 60*pow(b,2)*c*pow(d,3)*f*g*pow(influence,2) + 45*a*pow(d,4)*f*g*pow(influence,2) + 15*pow(a,3)*pow(d,4)*f*g*pow(influence,2) + 45*a*pow(b,2)*pow(d,4)*f*g*pow(influence,2) + 270*b*pow(c,2)*d*pow(f,2)*g*pow(influence,2) + 270*a*b*c*pow(d,2)*pow(f,2)*g*pow(influence,2) + 90*b*pow(d,3)*pow(f,2)*g*pow(influence,2) + 90*pow(a,2)*b*pow(d,3)*pow(f,2)*g*pow(influence,2) + 30*pow(b,3)*pow(d,3)*pow(f,2)*g*pow(influence,2) + 360*pow(b,2)*c*d*pow(f,3)*g*pow(influence,2) + 180*a*pow(b,2)*pow(d,2)*pow(f,3)*g* pow(influence,2) + 150*pow(b,3)*d*pow(f,4)*g*pow(influence,2) - 30*pow(c,3)*e*pow(influence,4) - 90*a*pow(c,2)*d*e*pow(influence,4) - 90*c*pow(d,2)*e*pow(influence,4) - 90*pow(a,2)*c*pow(d,2)*e*pow(influence,4) - 90*a*pow(d,3)*e*pow(influence,4) - 30*pow(a,3)*pow(d,3)*e*pow(influence,4) - 90*b*pow(c,2)*e*f*pow(influence,4) - 180*a*b*c*d*e*f*pow(influence,4) - 90*b*pow(d,2)*e*f*pow(influence,4) - 90*pow(a,2)*b*pow(d,2)*e*f*pow(influence,4) - 90*c*e*pow(f,2)*pow(influence,4) - 90*pow(b,2)*c*e*pow(f,2)*pow(influence,4) - 90*a*d*e*pow(f,2)*pow(influence,4) - 90*a*pow(b,2)*d*e*pow(f,2)*pow(influence,4) - 90*b*e*pow(f,3)*pow(influence,4) - 30*pow(b,3)*e*pow(f,3)*pow(influence,4) - 90*b*pow(c,2)*d*g*pow(influence,4) - 90*a*b*c*pow(d,2)*g*pow(influence,4) - 30*b*pow(d,3)*g*pow(influence,4) - 30*pow(a,2)*b*pow(d,3)*g*pow(influence,4) - 180*c*d*f*g*pow(influence,4) - 180*pow(b,2)*c*d*f*g*pow(influence,4) - 90*a*pow(d,2)*f*g*pow(influence,4) - 90*a*pow(b,2)*pow(d,2)*f*g*pow(influence,4) - 270*b*d*pow(f,2)*g*pow(influence,4) - 90*pow(b,3)*d*pow(f,2)*g*pow(influence,4) + 90*c*e*pow(influence,6) + 90*a*d*e*pow(influence,6) + 90*b*e*f*pow(influence,6) + 90*b*d*g*pow(influence,6)) + ((-60*pow(c,3)*d*pow(e,2)*pow(f,2) - 270*a*pow(c,2)*pow(d,2)*pow(e,2)*pow(f,2) - 360*pow(a,2)*c*pow(d,3)*pow(e,2)*pow(f,2) - 150*pow(a,3)*pow(d,4)*pow(e,2)*pow(f,2) - 180*b*pow(c,2)*d*pow(e,2)*pow(f,3) - 540*a*b*c*pow(d,2)*pow(e,2)*pow(f,3) - 360*pow(a,2)*b*pow(d,3)*pow(e,2)*pow(f,3) - 180*pow(b,2)*c*d*pow(e,2)*pow(f,4) - 270*a*pow(b,2)*pow(d,2)*pow(e,2)*pow(f,4) - 60*pow(b,3)*d*pow(e,2)*pow(f,5) - 120*pow(c,3)*pow(d,2)*e*f*g - 360*a*pow(c,2)*pow(d,3)*e*f*g - 360*pow(a,2)*c*pow(d,4)*e*f*g - 120*pow(a,3)*pow(d,5)*e*f*g - 540*b*pow(c,2)*pow(d,2)*e*pow(f,2)*g - 1080*a*b*c*pow(d,3)*e*pow(f,2)* g - 540*pow(a,2)*b*pow(d,4)*e*pow(f,2)*g - 720*pow(b,2)*c*pow(d,2)*e* pow(f,3)*g - 720*a*pow(b,2)*pow(d,3)*e*pow(f,3)*g - 300*pow(b,3)*pow(d,2)*e*pow(f,4)*g - 20*pow(c,3)*pow(d,3)*pow(g,2) - 45*a*pow(c,2)*pow(d,4)*pow(g,2) - 36*pow(a,2)*c*pow(d,5)*pow(g,2) - 10*pow(a,3)*pow(d,6)*pow(g,2) - 180*b*pow(c,2)*pow(d,3)*f*pow(g,2) - 270*a*b*c*pow(d,4)*f*pow(g,2) - 108*pow(a,2)*b*pow(d,5)*f*pow(g,2) - 360*pow(b,2)*c*pow(d,3)*pow(f,2)*pow(g,2) - 270*a*pow(b,2)*pow(d,4)*pow(f,2)*pow(g,2) - 200*pow(b,3)*pow(d,3)*pow(f,3)*pow(g,2) + 60*pow(c,3)*d*pow(e,2)*pow(influence,2) + 270*a*pow(c,2)*pow(d,2)*pow(e,2)*pow(influence,2) + 360*pow(a,2)*c*pow(d,3)*pow(e,2)*pow(influence,2) + 150*pow(a,3)*pow(d,4)*pow(e,2)*pow(influence,2) + 180*b*pow(c,2)*d*pow(e,2)*f*pow(influence,2) + 540*a*b*c*pow(d,2)*pow(e,2)*f*pow(influence,2) + 360*pow(a,2)*b*pow(d,3)*pow(e,2)*f* pow(influence,2) + 90*a*pow(c,2)*pow(e,2)*pow(f,2)*pow(influence,2) + 180*c*d*pow(e,2)*pow(f,2)* pow(influence,2) + 180*pow(a,2)*c*d*pow(e,2)*pow(f,2)*pow(influence,2) + 180*pow(b,2)*c*d*pow(e,2)*pow(f,2)*pow(influence,2) + 270*a*pow(d,2)*pow(e,2)*pow(f,2)* pow(influence,2) + 90*pow(a,3)*pow(d,2)*pow(e,2)*pow(f,2)*pow(influence,2) + 270*a*pow(b,2)*pow(d,2)*pow(e,2)*pow(f,2)*pow(influence,2) + 180*a*b*c*pow(e,2)*pow(f,3)* pow(influence,2) + 180*b*d*pow(e,2)*pow(f,3)*pow(influence,2) + 180*pow(a,2)*b*d*pow(e,2)* pow(f,3)*pow(influence,2) + 60*pow(b,3)*d*pow(e,2)*pow(f,3)*pow(influence,2) + 90*a*pow(b,2)*pow(e,2)*pow(f,4)*pow(influence,2) + 180*b*pow(c,2)*pow(d,2)*e*g*pow(influence,2) + 360*a*b*c*pow(d,3)*e*g*pow(influence,2) + 180*pow(a,2)*b*pow(d,4)*e*g*pow(influence,2) + 120*pow(c,3)*e*f*g*pow(influence,2) + 360*a*pow(c,2)*d*e*f*g*pow(influence,2) + 360*c*pow(d,2)*e*f*g*pow(influence,2) + 360*pow(a,2)*c*pow(d,2)*e*f*g*pow(influence,2) + 360*pow(b,2)*c*pow(d,2)*e*f*g*pow(influence,2) + 360*a*pow(d,3)*e*f*g*pow(influence,2) + 120*pow(a,3)*pow(d,3)*e*f*g*pow(influence,2) + 360*a*pow(b,2)*pow(d,3)*e*f*g* pow(influence,2) + 540*b*pow(c,2)*e*pow(f,2)*g*pow(influence,2) + 1080*a*b*c*d*e*pow(f,2)*g*pow(influence,2) + 540*b*pow(d,2)*e*pow(f,2)*g* pow(influence,2) + 540*pow(a,2)*b*pow(d,2)*e*pow(f,2)*g*pow(influence,2) + 180*pow(b,3)*pow(d,2)*e*pow(f,2)*g*pow(influence,2) + 720*pow(b,2)*c*e*pow(f,3)*g* pow(influence,2) + 720*a*pow(b,2)*d*e*pow(f,3)*g*pow(influence,2) + 300*pow(b,3)*e*pow(f,4)*g*pow(influence,2) + 60*pow(c,3)*d*pow(g,2)*pow(influence,2) + 90*a*pow(c,2)*pow(d,2)*pow(g,2)*pow(influence,2) + 60*c*pow(d,3)*pow(g,2)*pow(influence,2) + 60*pow(a,2)*c*pow(d,3)*pow(g,2)*pow(influence,2) + 60*pow(b,2)*c*pow(d,3)*pow(g,2)*pow(influence,2) + 45*a*pow(d,4)*pow(g,2)*pow(influence,2) + 15*pow(a,3)*pow(d,4)*pow(g,2)*pow(influence,2) + 45*a*pow(b,2)*pow(d,4)*pow(g,2)*pow(influence,2) + 540*b*pow(c,2)*d*f*pow(g,2)*pow(influence,2) + 540*a*b*c*pow(d,2)*f*pow(g,2)*pow(influence,2) + 180*b*pow(d,3)*f*pow(g,2)*pow(influence,2) + 180*pow(a,2)*b*pow(d,3)*f*pow(g,2)*pow(influence,2) + 60*pow(b,3)*pow(d,3)*f*pow(g,2)* pow(influence,2) + 1080*pow(b,2)*c*d*pow(f,2)*pow(g,2)*pow(influence,2) + 540*a*pow(b,2)*pow(d,2)*pow(f,2)*pow(g,2)*pow(influence,2) + 600*pow(b,3)*d*pow(f,3)*pow(g,2)* pow(influence,2) - 90*a*pow(c,2)*pow(e,2)*pow(influence,4) - 180*c*d*pow(e,2)*pow(influence,4) - 180*pow(a,2)*c*d*pow(e,2)*pow(influence,4) - 270*a*pow(d,2)*pow(e,2)*pow(influence,4) - 90*pow(a,3)*pow(d,2)*pow(e,2)*pow(influence,4) - 180*a*b*c*pow(e,2)*f*pow(influence,4) - 180*b*d*pow(e,2)*f*pow(influence,4) - 180*pow(a,2)*b*d*pow(e,2)*f*pow(influence,4) - 90*a*pow(e,2)*pow(f,2)*pow(influence,4) - 90*a*pow(b,2)*pow(e,2)*pow(f,2)*pow(influence,4) - 180*b*pow(c,2)*e*g*pow(influence,4) - 360*a*b*c*d*e*g*pow(influence,4) - 180*b*pow(d,2)*e*g*pow(influence,4) - 180*pow(a,2)*b*pow(d,2)*e*g*pow(influence,4) - 360*c*e*f*g*pow(influence,4) - 360*pow(b,2)*c*e*f*g*pow(influence,4) - 360*a*d*e*f*g*pow(influence,4) - 360*a*pow(b,2)*d*e*f*g*pow(influence,4) - 540*b*e*pow(f,2)*g*pow(influence,4) - 180*pow(b,3)*e*pow(f,2)*g*pow(influence,4) - 180*c*d*pow(g,2)*pow(influence,4) - 180*pow(b,2)*c*d*pow(g,2)*pow(influence,4) - 90*a*pow(d,2)*pow(g,2)*pow(influence,4) - 90*a*pow(b,2)*pow(d,2)*pow(g,2)*pow(influence,4) - 540*b*d*f*pow(g,2)*pow(influence,4) - 180*pow(b,3)*d*f*pow(g,2)*pow(influence,4) + 90*a*pow(e,2)*pow(influence,6) + 180*b*e*g*pow(influence,6)))/3 + ((-10*pow(c,3)*pow(e,3)*pow(f,2) - 90*a*pow(c,2)*d*pow(e,3)*pow(f,2) - 180*pow(a,2)*c*pow(d,2)*pow(e,3)*pow(f,2) - 100*pow(a,3)*pow(d,3)*pow(e,3)*pow(f,2) - 30*b*pow(c,2)*pow(e,3)*pow(f,3) - 180*a*b*c*d*pow(e,3)*pow(f,3) - 180*pow(a,2)*b*pow(d,2)*pow(e,3)*pow(f,3) - 30*pow(b,2)*c*pow(e,3)*pow(f,4) - 90*a*pow(b,2)*d*pow(e,3)*pow(f,4) - 10*pow(b,3)*pow(e,3)*pow(f,5) - 60*pow(c,3)*d*pow(e,2)*f*g - 270*a*pow(c,2)*pow(d,2)*pow(e,2)*f*g - 360*pow(a,2)*c*pow(d,3)*pow(e,2)*f*g - 150*pow(a,3)*pow(d,4)*pow(e,2)*f*g - 270*b*pow(c,2)*d*pow(e,2)*pow(f,2)*g - 810*a*b*c*pow(d,2)*pow(e,2)*pow(f,2)* g - 540*pow(a,2)*b*pow(d,3)*pow(e,2)*pow(f,2)*g - 360*pow(b,2)*c*d*pow(e,2)*pow(f,3)*g - 540*a*pow(b,2)*pow(d,2)*pow(e,2)*pow(f,3)* g - 150*pow(b,3)*d*pow(e,2)*pow(f,4)*g - 30*pow(c,3)*pow(d,2)*e*pow(g,2) - 90*a*pow(c,2)*pow(d,3)*e*pow(g,2) - 90*pow(a,2)*c*pow(d,4)*e*pow(g,2) - 30*pow(a,3)*pow(d,5)*e*pow(g,2) - 270*b*pow(c,2)*pow(d,2)*e*f*pow(g,2) - 540*a*b*c*pow(d,3)*e*f*pow(g,2) - 270*pow(a,2)*b*pow(d,4)*e*f*pow(g,2) - 540*pow(b,2)*c*pow(d,2)*e*pow(f,2)*pow(g,2) - 540*a*pow(b,2)*pow(d,3)*e*pow(f,2)* pow(g,2) - 300*pow(b,3)*pow(d,2)*e*pow(f,3)*pow(g,2) - 30*b*pow(c,2)*pow(d,3)*pow(g,3) - 45*a*b*c*pow(d,4)*pow(g,3) - 18*pow(a,2)*b*pow(d,5)*pow(g,3) - 120*pow(b,2)*c*pow(d,3)*f*pow(g,3) - 90*a*pow(b,2)*pow(d,4)*f*pow(g,3) - 100*pow(b,3)*pow(d,3)*pow(f,2)*pow(g,3) + 10*pow(c,3)*pow(e,3)*pow(influence,2) + 90*a*pow(c,2)*d*pow(e,3)*pow(influence,2) + 180*pow(a,2)*c*pow(d,2)*pow(e,3)*pow(influence,2) + 100*pow(a,3)*pow(d,3)*pow(e,3)*pow(influence,2) + 30*b*pow(c,2)*pow(e,3)*f*pow(influence,2) + 180*a*b*c*d*pow(e,3)*f*pow(influence,2) + 180*pow(a,2)*b*pow(d,2)*pow(e,3)*f*pow(influence,2) + 30*c*pow(e,3)*pow(f,2)*pow(influence,2) + 30*pow(a,2)*c*pow(e,3)*pow(f,2)*pow(influence,2) + 30*pow(b,2)*c*pow(e,3)*pow(f,2)*pow(influence,2) + 90*a*d*pow(e,3)*pow(f,2)*pow(influence,2) + 30*pow(a,3)*d*pow(e,3)*pow(f,2)*pow(influence,2) + 90*a*pow(b,2)*d*pow(e,3)*pow(f,2)*pow(influence,2) + 30*b*pow(e,3)*pow(f,3)*pow(influence,2) + 30*pow(a,2)*b*pow(e,3)*pow(f,3)*pow(influence,2) + 10*pow(b,3)*pow(e,3)*pow(f,3)*pow(influence,2) + 90*b*pow(c,2)*d*pow(e,2)*g*pow(influence,2) + 270*a*b*c*pow(d,2)*pow(e,2)*g* pow(influence,2) + 180*pow(a,2)*b*pow(d,3)*pow(e,2)*g*pow(influence,2) + 90*a*pow(c,2)*pow(e,2)*f*g*pow(influence,2) + 180*c*d*pow(e,2)*f*g*pow(influence,2) + 180*pow(a,2)*c*d*pow(e,2)*f*g*pow(influence,2) + 180*pow(b,2)*c*d*pow(e,2)*f*g* pow(influence,2) + 270*a*pow(d,2)*pow(e,2)*f*g*pow(influence,2) + 90*pow(a,3)*pow(d,2)*pow(e,2)*f*g*pow(influence,2) + 270*a*pow(b,2)*pow(d,2)*pow(e,2)*f*g* pow(influence,2) + 270*a*b*c*pow(e,2)*pow(f,2)*g*pow(influence,2) + 270*b*d*pow(e,2)*pow(f,2)*g*pow(influence,2) + 270*pow(a,2)*b*d*pow(e,2)*pow(f,2)*g* pow(influence,2) + 90*pow(b,3)*d*pow(e,2)*pow(f,2)*g*pow(influence,2) + 180*a*pow(b,2)*pow(e,2)*pow(f,3)*g*pow(influence,2) + 30*pow(c,3)*e*pow(g,2)*pow(influence,2) + 90*a*pow(c,2)*d*e*pow(g,2)*pow(influence,2) + 90*c*pow(d,2)*e*pow(g,2)*pow(influence,2) + 90*pow(a,2)*c*pow(d,2)*e*pow(g,2)*pow(influence,2) + 90*pow(b,2)*c*pow(d,2)*e*pow(g,2)* pow(influence,2) + 90*a*pow(d,3)*e*pow(g,2)*pow(influence,2) + 30*pow(a,3)*pow(d,3)*e*pow(g,2)* pow(influence,2) + 90*a*pow(b,2)*pow(d,3)*e*pow(g,2)*pow(influence,2) + 270*b*pow(c,2)*e*f*pow(g,2)*pow(influence,2) + 540*a*b*c*d*e*f*pow(g,2)* pow(influence,2) + 270*b*pow(d,2)*e*f*pow(g,2)*pow(influence,2) + 270*pow(a,2)*b*pow(d,2)*e*f*pow(g,2)*pow(influence,2) + 90*pow(b,3)*pow(d,2)*e*f*pow(g,2)* pow(influence,2) + 540*pow(b,2)*c*e*pow(f,2)*pow(g,2)*pow(influence,2) + 540*a*pow(b,2)*d*e*pow(f,2)*pow(g,2)*pow(influence,2) + 300*pow(b,3)*e*pow(f,3)*pow(g,2)* pow(influence,2) + 90*b*pow(c,2)*d*pow(g,3)*pow(influence,2) + 90*a*b*c*pow(d,2)*pow(g,3)* pow(influence,2) + 30*b*pow(d,3)*pow(g,3)*pow(influence,2) + 30*pow(a,2)*b*pow(d,3)*pow(g,3)* pow(influence,2) + 10*pow(b,3)*pow(d,3)*pow(g,3)*pow(influence,2) + 360*pow(b,2)*c*d*f*pow(g,3)* pow(influence,2) + 180*a*pow(b,2)*pow(d,2)*f*pow(g,3)*pow(influence,2) + 300*pow(b,3)*d*pow(f,2)*pow(g,3)*pow(influence,2) - 30*c*pow(e,3)*pow(influence,4) - 30*pow(a,2)*c*pow(e,3)*pow(influence,4) - 90*a*d*pow(e,3)*pow(influence,4) - 30*pow(a,3)*d*pow(e,3)*pow(influence,4) - 30*b*pow(e,3)*f*pow(influence,4) - 30*pow(a,2)*b*pow(e,3)*f*pow(influence,4) - 90*a*b*c*pow(e,2)*g*pow(influence,4) - 90*b*d*pow(e,2)*g*pow(influence,4) - 90*pow(a,2)*b*d*pow(e,2)*g*pow(influence,4) - 90*a*pow(e,2)*f*g*pow(influence,4) - 90*a*pow(b,2)*pow(e,2)*f*g*pow(influence,4) - 90*c*e*pow(g,2)*pow(influence,4) - 90*pow(b,2)*c*e*pow(g,2)*pow(influence,4) - 90*a*d*e*pow(g,2)*pow(influence,4) - 90*a*pow(b,2)*d*e*pow(g,2)*pow(influence,4) - 270*b*e*f*pow(g,2)*pow(influence,4) - 90*pow(b,3)*e*f*pow(g,2)*pow(influence,4) - 90*b*d*pow(g,3)*pow(influence,4) - 30*pow(b,3)*d*pow(g,3)*pow(influence,4)))/2 + (-9*a*pow(c,2)*pow(e,4)*pow(f,2) - 36*pow(a,2)*c*d*pow(e,4)*pow(f,2) - 30*pow(a,3)*pow(d,2)*pow(e,4)*pow(f,2) - 18*a*b*c*pow(e,4)*pow(f,3) - 36*pow(a,2)*b*d*pow(e,4)*pow(f,3) - 9*a*pow(b,2)*pow(e,4)*pow(f,4) - 8*pow(c,3)*pow(e,3)*f*g - 72*a*pow(c,2)*d*pow(e,3)*f*g - 144*pow(a,2)*c*pow(d,2)*pow(e,3)*f*g - 80*pow(a,3)*pow(d,3)*pow(e,3)*f*g - 36*b*pow(c,2)*pow(e,3)*pow(f,2)*g - 216*a*b*c*d*pow(e,3)*pow(f,2)*g - 216*pow(a,2)*b*pow(d,2)*pow(e,3)*pow(f,2)*g - 48*pow(b,2)*c*pow(e,3)*pow(f,3)*g - 144*a*pow(b,2)*d*pow(e,3)*pow(f,3)*g - 20*pow(b,3)*pow(e,3)*pow(f,4)*g - 12*pow(c,3)*d*pow(e,2)*pow(g,2) - 54*a*pow(c,2)*pow(d,2)*pow(e,2)*pow(g,2) - 72*pow(a,2)*c*pow(d,3)*pow(e,2)*pow(g,2) - 30*pow(a,3)*pow(d,4)*pow(e,2)*pow(g,2) - 108*b*pow(c,2)*d*pow(e,2)*f*pow(g,2) - 324*a*b*c*pow(d,2)*pow(e,2)*f* pow(g,2) - 216*pow(a,2)*b*pow(d,3)*pow(e,2)*f*pow(g,2) - 216*pow(b,2)*c*d*pow(e,2)*pow(f,2)*pow(g,2) - 324*a*pow(b,2)*pow(d,2)*pow(e,2)*pow(f,2)* pow(g,2) - 120*pow(b,3)*d*pow(e,2)*pow(f,3)*pow(g,2) - 36*b*pow(c,2)*pow(d,2)*e*pow(g,3) - 72*a*b*c*pow(d,3)*e*pow(g,3) - 36*pow(a,2)*b*pow(d,4)*e*pow(g,3) - 144*pow(b,2)*c*pow(d,2)*e*f*pow(g,3) - 144*a*pow(b,2)*pow(d,3)*e*f*pow(g,3) - 120*pow(b,3)*pow(d,2)*e*pow(f,2)*pow(g,3) - 12*pow(b,2)*c*pow(d,3)*pow(g,4) - 9*a*pow(b,2)*pow(d,4)*pow(g,4) - 20*pow(b,3)*pow(d,3)*f*pow(g,4) + 9*a*pow(c,2)*pow(e,4)*pow(influence,2) + 36*pow(a,2)*c*d*pow(e,4)*pow(influence,2) + 30*pow(a,3)*pow(d,2)*pow(e,4)*pow(influence,2) + 18*a*b*c*pow(e,4)*f*pow(influence,2) + 36*pow(a,2)*b*d*pow(e,4)*f*pow(influence,2) + 9*a*pow(e,4)*pow(f,2)*pow(influence,2) + 3*pow(a,3)*pow(e,4)*pow(f,2)*pow(influence,2) + 9*a*pow(b,2)*pow(e,4)*pow(f,2)*pow(influence,2) + 12*b*pow(c,2)*pow(e,3)*g*pow(influence,2) + 72*a*b*c*d*pow(e,3)*g*pow(influence,2) + 72*pow(a,2)*b*pow(d,2)*pow(e,3)*g*pow(influence,2) + 24*c*pow(e,3)*f*g*pow(influence,2) + 24*pow(a,2)*c*pow(e,3)*f*g*pow(influence,2) + 24*pow(b,2)*c*pow(e,3)*f*g*pow(influence,2) + 72*a*d*pow(e,3)*f*g*pow(influence,2) + 24*pow(a,3)*d*pow(e,3)*f*g*pow(influence,2) + 72*a*pow(b,2)*d*pow(e,3)*f*g*pow(influence,2) + 36*b*pow(e,3)*pow(f,2)*g*pow(influence,2) + 36*pow(a,2)*b*pow(e,3)*pow(f,2)*g*pow(influence,2) + 12*pow(b,3)*pow(e,3)*pow(f,2)*g*pow(influence,2) + 18*a*pow(c,2)*pow(e,2)*pow(g,2)*pow(influence,2) + 36*c*d*pow(e,2)*pow(g,2)*pow(influence,2) + 36*pow(a,2)*c*d*pow(e,2)*pow(g,2)*pow(influence,2) + 36*pow(b,2)*c*d*pow(e,2)*pow(g,2)*pow(influence,2) + 54*a*pow(d,2)*pow(e,2)*pow(g,2)*pow(influence,2) + 18*pow(a,3)*pow(d,2)*pow(e,2)*pow(g,2)*pow(influence,2) + 54*a*pow(b,2)*pow(d,2)*pow(e,2)*pow(g,2)* pow(influence,2) + 108*a*b*c*pow(e,2)*f*pow(g,2)*pow(influence,2) + 108*b*d*pow(e,2)*f*pow(g,2)*pow(influence,2) + 108*pow(a,2)*b*d*pow(e,2)*f*pow(g,2)* pow(influence,2) + 36*pow(b,3)*d*pow(e,2)*f*pow(g,2)*pow(influence,2) + 108*a*pow(b,2)*pow(e,2)*pow(f,2)*pow(g,2)*pow(influence,2) + 36*b*pow(c,2)*e*pow(g,3)*pow(influence,2) + 72*a*b*c*d*e*pow(g,3)*pow(influence,2) + 36*b*pow(d,2)*e*pow(g,3)*pow(influence,2) + 36*pow(a,2)*b*pow(d,2)*e*pow(g,3)*pow(influence,2) + 12*pow(b,3)*pow(d,2)*e*pow(g,3)*pow(influence,2) + 144*pow(b,2)*c*e*f*pow(g,3)*pow(influence,2) + 144*a*pow(b,2)*d*e*f*pow(g,3)* pow(influence,2) + 120*pow(b,3)*e*pow(f,2)*pow(g,3)*pow(influence,2) + 36*pow(b,2)*c*d*pow(g,4)*pow(influence,2) + 18*a*pow(b,2)*pow(d,2)*pow(g,4)*pow(influence,2) + 60*pow(b,3)*d*f*pow(g,4)*pow(influence,2) - 9*a*pow(e,4)*pow(influence,4) - 3*pow(a,3)*pow(e,4)*pow(influence,4) - 12*b*pow(e,3)*g*pow(influence,4) - 12*pow(a,2)*b*pow(e,3)*g*pow(influence,4) - 18*a*pow(e,2)*pow(g,2)*pow(influence,4) - 18*a*pow(b,2)*pow(e,2)*pow(g,2)*pow(influence,4) - 36*b*e*pow(g,3)*pow(influence,4) - 12*pow(b,3)*e*pow(g,3)*pow(influence,4)) + ((-18*pow(a,2)*c*pow(e,5)*pow(f,2) - 30*pow(a,3)*d*pow(e,5)*pow(f,2) - 18*pow(a,2)*b*pow(e,5)*pow(f,3) - 45*a*pow(c,2)*pow(e,4)*f*g - 180*pow(a,2)*c*d*pow(e,4)*f*g - 150*pow(a,3)*pow(d,2)*pow(e,4)*f*g - 135*a*b*c*pow(e,4)*pow(f,2)*g - 270*pow(a,2)*b*d*pow(e,4)*pow(f,2)*g - 90*a*pow(b,2)*pow(e,4)*pow(f,3)*g - 10*pow(c,3)*pow(e,3)*pow(g,2) - 90*a*pow(c,2)*d*pow(e,3)*pow(g,2) - 180*pow(a,2)*c*pow(d,2)*pow(e,3)*pow(g,2) - 100*pow(a,3)*pow(d,3)*pow(e,3)*pow(g,2) - 90*b*pow(c,2)*pow(e,3)*f*pow(g,2) - 540*a*b*c*d*pow(e,3)*f*pow(g,2) - 540*pow(a,2)*b*pow(d,2)*pow(e,3)*f* pow(g,2) - 180*pow(b,2)*c*pow(e,3)*pow(f,2)*pow(g,2) - 540*a*pow(b,2)*d*pow(e,3)*pow(f,2)*pow(g,2) - 100*pow(b,3)*pow(e,3)*pow(f,3)*pow(g,2) - 90*b*pow(c,2)*d*pow(e,2)*pow(g,3) - 270*a*b*c*pow(d,2)*pow(e,2)*pow(g,3) - 180*pow(a,2)*b*pow(d,3)*pow(e,2)*pow(g,3) - 360*pow(b,2)*c*d*pow(e,2)*f*pow(g,3) - 540*a*pow(b,2)*pow(d,2)*pow(e,2)*f*pow(g,3) - 300*pow(b,3)*d*pow(e,2)*pow(f,2)* pow(g,3) - 90*pow(b,2)*c*pow(d,2)*e*pow(g,4) - 90*a*pow(b,2)*pow(d,3)*e* pow(g,4) - 150*pow(b,3)*pow(d,2)*e*f*pow(g,4) - 10*pow(b,3)*pow(d,3)*pow(g,5) + 18*pow(a,2)*c*pow(e,5)*pow(influence,2) + 30*pow(a,3)*d*pow(e,5)*pow(influence,2) + 18*pow(a,2)*b*pow(e,5)*f*pow(influence,2) + 45*a*b*c*pow(e,4)*g*pow(influence,2) + 90*pow(a,2)*b*d*pow(e,4)*g*pow(influence,2) + 45*a*pow(e,4)*f*g*pow(influence,2) + 15*pow(a,3)*pow(e,4)*f*g*pow(influence,2) + 45*a*pow(b,2)*pow(e,4)*f*g*pow(influence,2) + 30*c*pow(e,3)*pow(g,2)*pow(influence,2) + 30*pow(a,2)*c*pow(e,3)*pow(g,2)*pow(influence,2) + 30*pow(b,2)*c*pow(e,3)*pow(g,2)*pow(influence,2) + 90*a*d*pow(e,3)*pow(g,2)*pow(influence,2) + 30*pow(a,3)*d*pow(e,3)*pow(g,2)*pow(influence,2) + 90*a*pow(b,2)*d*pow(e,3)*pow(g,2)*pow(influence,2) + 90*b*pow(e,3)*f*pow(g,2)*pow(influence,2) + 90*pow(a,2)*b*pow(e,3)*f*pow(g,2)*pow(influence,2) + 30*pow(b,3)*pow(e,3)*f*pow(g,2)*pow(influence,2) + 90*a*b*c*pow(e,2)*pow(g,3)*pow(influence,2) + 90*b*d*pow(e,2)*pow(g,3)*pow(influence,2) + 90*pow(a,2)*b*d*pow(e,2)*pow(g,3)*pow(influence,2) + 30*pow(b,3)*d*pow(e,2)*pow(g,3)*pow(influence,2) + 180*a*pow(b,2)*pow(e,2)*f*pow(g,3)* pow(influence,2) + 90*pow(b,2)*c*e*pow(g,4)*pow(influence,2) + 90*a*pow(b,2)*d*e*pow(g,4)* pow(influence,2) + 150*pow(b,3)*e*f*pow(g,4)*pow(influence,2) + 30*pow(b,3)*d*pow(g,5)*pow(influence,2)))/3 + ((-10*pow(a,3)*pow(e,6)*pow(f,2) - 72*pow(a,2)*c*pow(e,5)*f*g - 120*pow(a,3)*d*pow(e,5)*f*g - 108*pow(a,2)*b*pow(e,5)*pow(f,2)*g - 45*a*pow(c,2)*pow(e,4)*pow(g,2) - 180*pow(a,2)*c*d*pow(e,4)*pow(g,2) - 150*pow(a,3)*pow(d,2)*pow(e,4)*pow(g,2) - 270*a*b*c*pow(e,4)*f*pow(g,2) - 540*pow(a,2)*b*d*pow(e,4)*f*pow(g,2) - 270*a*pow(b,2)*pow(e,4)*pow(f,2)*pow(g,2) - 60*b*pow(c,2)*pow(e,3)*pow(g,3) - 360*a*b*c*d*pow(e,3)*pow(g,3) - 360*pow(a,2)*b*pow(d,2)*pow(e,3)*pow(g,3) - 240*pow(b,2)*c*pow(e,3)*f*pow(g,3) - 720*a*pow(b,2)*d*pow(e,3)*f*pow(g,3) - 200*pow(b,3)*pow(e,3)*pow(f,2)*pow(g,3) - 180*pow(b,2)*c*d*pow(e,2)*pow(g,4) - 270*a*pow(b,2)*pow(d,2)*pow(e,2)*pow(g,4) - 300*pow(b,3)*d*pow(e,2)*f*pow(g,4) - 60*pow(b,3)*pow(d,2)*e*pow(g,5) + 10*pow(a,3)*pow(e,6)*pow(influence,2) + 36*pow(a,2)*b*pow(e,5)*g*pow(influence,2) + 45*a*pow(e,4)*pow(g,2)*pow(influence,2) + 15*pow(a,3)*pow(e,4)*pow(g,2)*pow(influence,2) + 45*a*pow(b,2)*pow(e,4)*pow(g,2)*pow(influence,2) + 60*b*pow(e,3)*pow(g,3)*pow(influence,2) + 60*pow(a,2)*b*pow(e,3)*pow(g,3)*pow(influence,2) + 20*pow(b,3)*pow(e,3)*pow(g,3)*pow(influence,2) + 90*a*pow(b,2)*pow(e,2)*pow(g,4)*pow(influence,2) + 60*pow(b,3)*e*pow(g,5)*pow(influence,2)))/ 7 - (pow(e,2)*g*(10*pow(a,3)*pow(e,4)*f + 18*pow(a,2)*c*pow(e,3)*g + 30*pow(a,3)*d*pow(e,3)*g + 54*pow(a,2)*b*pow(e,3)*f*g + 45*a*b*c*pow(e,2)*pow(g,2) + 90*pow(a,2)*b*d*pow(e,2)*pow(g,2) + 90*a*pow(b,2)*pow(e,2)*f*pow(g,2) + 30*pow(b,2)*c*e*pow(g,3) + 90*a*pow(b,2)*d*e*pow(g,3) + 50*pow(b,3)*e*f*pow(g,3) + 30*pow(b,3)*d*pow(g,4)))/4 - (pow(e,3)*pow(g,2)*(10*pow(a,3)*pow(e,3) + 36*pow(a,2)*b*pow(e,2)*g + 45*a*pow(b,2)*e*pow(g,2) + 20*pow(b,3)*pow(g,3)))/9))/180;
  }

  return;
}


/*----------------------------------------------------------------------*
 | test for overlap of the XAABB of bubble and element     ghamm 11/14  |
 *----------------------------------------------------------------------*/
bool CAVITATION::Algorithm::XAABBoverlap(
  DRT::Element* ele,
  const double influence,
  LINALG::Matrix<3,1> particleposition
  )
{
  // get bounding box of current element
  double xaabb[6] = { std::numeric_limits<double>::max(),  std::numeric_limits<double>::max(),
                      std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
  for (int inode=0; inode<ele->NumNode(); ++inode)
  {
    const DRT::Node* node  = ele->Nodes()[inode];
    const double* coord = node->X();
    for (size_t i = 0; i < 3; ++i)
    {
      xaabb[i+0] = std::min(xaabb[i+0], coord[i]);
      xaabb[i+3] = std::max(xaabb[i+3], coord[i]);
    }
  }

  double bubblesurface[6] = {particleposition(0)+influence, particleposition(1)+influence,
                             particleposition(2)+influence, particleposition(0)-influence,
                             particleposition(1)-influence, particleposition(2)-influence};

  bool overlap = true;
  // test whether the bounding box of the fluid element touches the bubble influence
  for(int dim=0; dim<3; ++dim)
  {
    if(xaabb[dim] - GEO::TOL7 > bubblesurface[dim] or xaabb[dim+3] + GEO::TOL7 < bubblesurface[dim+3])
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
void CAVITATION::Algorithm::AssignSmallBubbles(
  const double bubblevol,
  const LINALG::Matrix<3,1> particleposition,
  const std::vector<int> insideeles,
  Teuchos::RCP<Epetra_FEVector> void_volumes
  )
{
  // safety check
  if(insideeles.size() < 1)
    dserror("no underlying fluid element for void fraction computation found.");

  // a very small bubble was processed which cannot be resolved with volume integration
  // just add the contribution to the underlying fluid element
  if(insideeles.size() == 1)
  {
    // do assembly of void fraction into fluid
    int err = void_volumes->SumIntoGlobalValues(1, &insideeles[0], &bubblevol);
    if (err<0)
      dserror("summing into Epetra_FEVector failed");
#ifdef DEBUG
    DRT::Element* fluidele = fluiddis_->gElement(insideeles[0]);
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

    static LINALG::Matrix<3,1> dummy(false);
    // get coordinates of the particle position in parameter space of the element
    bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, dummy);

    if(insideele == false)
      dserror("bubble is expected to lie inside this fluid element");
#endif
  }
  else
  {
    bool insideele = false;
    for(size_t i=0; i<insideeles.size(); ++i)
    {
      DRT::Element* fluidele = fluiddis_->gElement(insideeles[i]);
      const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

      static LINALG::Matrix<3,1> dummy(false);
      // get coordinates of the particle position in parameter space of the element
      insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, dummy);

      if(insideele == true)
      {
        // do assembly of void fraction into fluid
        int err = void_volumes->SumIntoGlobalValues(1, &insideeles[i], &bubblevol);
        if (err<0)
          dserror("summing into Epetra_FEVector failed");
        break;
      }
    }
    if(insideele == false)
      dserror("void fraction could not be assigned to an underlying fluid element.");
  }

  return;
}


/*----------------------------------------------------------------------*
 | get underlying element as well as position in element space          |
 |                                                          ghamm 06/15 |
 *----------------------------------------------------------------------*/
DRT::Element* CAVITATION::Algorithm::GetEleCoordinatesFromPosition(
  LINALG::Matrix<3,1>& myposition,          ///< position
  DRT::MESHFREE::MeshfreeMultiBin* currbin, ///< corresponding bin
  LINALG::Matrix<3,1>& elecoord,            ///< matrix to be filled with particle coordinates in element space
  const bool approxelecoordsinit            ///< bool whether an inital approx. of the ele coords is used
  )
{
  // variables to store information about element in which the particle is located
  DRT::Element* targetfluidele = NULL;
  elecoord.Clear();

  DRT::Element** fluidelesinbin = currbin->AssociatedFluidEles();
  int numfluidelesinbin = currbin->NumAssociatedFluidEle();

  std::set<int>::const_iterator eleiter;
  // search for underlying fluid element with fast search if desired
  for(int ele=0; ele<numfluidelesinbin; ++ele)
  {
    DRT::Element* fluidele = fluidelesinbin[ele];
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

    //get coordinates of the particle position in parameter space of the element
    bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, myposition, elecoord, approxelecoordsinit);

    if(insideele == true)
    {
      targetfluidele = fluidele;
      // leave loop over all fluid eles in bin
      break;
    }
  }

  // repeat search for underlying fluid element with standard search in case nothing was found
  if(targetfluidele == NULL and approxelecoordsinit == true)
  {
    for(int ele=0; ele<numfluidelesinbin; ++ele)
    {
      DRT::Element* fluidele = fluidelesinbin[ele];
      const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

      //get coordinates of the particle position in parameter space of the element
      bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, myposition, elecoord, false);

      if(insideele == true)
      {
        targetfluidele = fluidele;
        // leave loop over all fluid eles in bin
        break;
      }
    }
  }

  return targetfluidele;
}


/*----------------------------------------------------------------------*
 | compute velocity at bubble position                     ghamm 12/15  |
 *----------------------------------------------------------------------*/
bool CAVITATION::Algorithm::ComputeVelocityAtBubblePosition(
  DRT::Node* currparticle,
  LINALG::Matrix<3,1>& particleposition,
  Epetra_SerialDenseMatrix& elemat1,
  Epetra_SerialDenseMatrix& elemat2,
  Epetra_SerialDenseVector& elevec1,
  Epetra_SerialDenseVector& elevec2,
  Epetra_SerialDenseVector& elevec3)
{
  // find out in which fluid element the current particle is located
  if(currparticle->NumElement() != 1)
    dserror("ERROR: A particle is assigned to more than one bin!");
  DRT::Element** currele = currparticle->Elements();
  DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);

  static LINALG::Matrix<3,1> elecoord(false);
  DRT::Element* targetfluidele = GetEleCoordinatesFromPosition(particleposition, currbin, elecoord, approxelecoordsinit_);

  if(targetfluidele == NULL)
  {
    std::cout << "WARNING: velocity for bubble (id: " << currparticle->Id() << " ) could not be computed at position: "
        << particleposition(0) << " " << particleposition(1) << " " << particleposition(2) << " on proc " << myrank_ << std::endl;
    return false;
  }

  // get element location vector and ownerships
  std::vector<int> lm_f;
  std::vector<int> lmowner_f;
  std::vector<int> lmstride;
  targetfluidele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

  // set action in order to compute velocity -> state with name "vel" expected inside
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::interpolate_velocity_to_given_point);
  params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

  // call the element specific evaluate method (elevec1 = fluid vel)
  targetfluidele->Evaluate(params,*fluiddis_,lm_f,elemat1,elemat2,elevec1,elevec2,elevec3);

  return true;
}


/*----------------------------------------------------------------------*
 | compute pressure at bubble position                     ghamm 06/15  |
 *----------------------------------------------------------------------*/
bool CAVITATION::Algorithm::ComputePressureAtBubblePosition(
  DRT::Node* currparticle,
  LINALG::Matrix<3,1>& particleposition,
  Epetra_SerialDenseMatrix& elemat1,
  Epetra_SerialDenseMatrix& elemat2,
  Epetra_SerialDenseVector& elevec1,
  Epetra_SerialDenseVector& elevec2,
  Epetra_SerialDenseVector& elevec3)
{
  // find out in which fluid element the current particle is located
  if(currparticle->NumElement() != 1)
    dserror("ERROR: A particle is assigned to more than one bin!");
  DRT::Element** currele = currparticle->Elements();
  DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);

  static LINALG::Matrix<3,1> elecoord(false);
  DRT::Element* targetfluidele = GetEleCoordinatesFromPosition(particleposition, currbin, elecoord, approxelecoordsinit_);

  // fill cache for force computation
  const size_t i   = currparticle->LID();
  if(underlyingelecache_.size() > i)
  {
    UnderlyingEle& e = underlyingelecache_[i];
    e.ele = targetfluidele;
    e.elecoord = elecoord;
  }

  if(targetfluidele == NULL)
  {
    std::cout << "WARNING: pressure for bubble (id: " << currparticle->Id() << " ) could not be computed at position: "
        << particleposition(0) << " " << particleposition(1) << " " << particleposition(2) << " on proc " << myrank_ << std::endl;
    return false;
  }

  // get element location vector and ownerships
  std::vector<int> lm_f;
  std::vector<int> lmowner_f;
  std::vector<int> lmstride;
  targetfluidele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

  // set action in order to compute pressure -> state with name "veln" (and "velnp") expected inside
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::interpolate_pressure_to_given_point);
  params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

  // call the element specific evaluate method (elevec1 = fluid press)
  targetfluidele->Evaluate(params,*fluiddis_,lm_f,elemat1,elemat2,elevec1,elevec2,elevec3);

  return true;
}


/*----------------------------------------------------------------------*
 | compute superconvergent patch recovery for fluid frac    ghamm 02/16 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ComputePatchRecoveredFluidFraction(
  Teuchos::RCP<const Epetra_MultiVector> fluidfraction
  )
{
  Teuchos::ParameterList params;
  // only element center is relevant in this action
  params.set<int>("action",FLD::calc_velgrad_ele_center);

  Teuchos::RCP<Epetra_MultiVector> nodevec;
  if(particle_dim_ == INPAR::PARTICLE::particle_3D)
    nodevec = DRT::UTILS::ComputeSuperconvergentPatchRecovery<3>(fluiddis_,Teuchos::rcp((*fluidfraction)(0),false), "dummy", 1, params);
  else
    nodevec = DRT::UTILS::ComputeSuperconvergentPatchRecovery<2>(fluiddis_,Teuchos::rcp((*fluidfraction)(0),false), "dummy", 1, params);

  const int numnode = nodevec->MyLength();
  for(int i=0; i<numnode; ++i)
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
  Teuchos::RCP<const Epetra_MultiVector> fluidfraction
  )
{
  if(not fluidfraction->Map().SameAs(*fluiddis_->ElementRowMap()))
    dserror("input map is not an element row map");
  Teuchos::RCP<Epetra_Vector> eleveccol = Teuchos::rcp(new Epetra_Vector(*fluiddis_->ElementColMap(),true));
  LINALG::Export(*fluidfraction, *eleveccol);

  // handle pbcs if existing
  // build inverse map from slave to master nodes
  std::map<int,int> slavetomastercolnodesmap;
  {
    Teuchos::RCP<std::map<int,std::vector<int> > > allcoupledcolnodes = fluiddis_->GetAllPBCCoupledColNodes();

    for(std::map<int,std::vector<int> >::const_iterator masterslavepair = allcoupledcolnodes->begin();
        masterslavepair != allcoupledcolnodes->end() ; ++masterslavepair)
    {
      // loop slave nodes associated with master
      for(std::vector<int>::const_iterator iter=masterslavepair->second.begin(); iter!=masterslavepair->second.end(); ++iter)
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
  for(int i=0; i<fullnoderowmap->NumMyElements(); ++i)
  {
    const int nodeid = fullnoderowmap->GID(i);
    // do not add slave pbc nodes here
    if(slavetomastercolnodesmap.count(nodeid) == 0)
      reducednoderowmap.push_back(nodeid);
  }

  // build node row map which does not include slave pbc nodes
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(-1,(int)reducednoderowmap.size(),&reducednoderowmap[0],0,fullnoderowmap->Comm()));

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> matrix = Teuchos::rcp(new LINALG::SparseMatrix(*noderowmap,108,false,true));
  // create empty right hand side
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*noderowmap,true);

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
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = fluiddis_->lColElement(i);
    const int numnode = actele->NumNode();

    // get element location vector for nodes
    lm.resize(numnode);
    lmowner.resize(numnode);

    DRT::Node** nodes = actele->Nodes();
    for(int n=0; n<numnode; ++n)
    {
      const int nodeid = nodes[n]->Id();

      std::map<int,int>::iterator slavemasterpair = slavetomastercolnodesmap.find(nodeid);
      if(slavemasterpair != slavetomastercolnodesmap.end())
        lm[n] = slavemasterpair->second;
      else
        lm[n] = nodeid;

      // owner of pbc master and slave nodes are identical
      lmowner[n] = nodes[n]->Owner();
    }

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(numnode);
    elematrix1.Shape(numnode,numnode);

    // set action in order to project element fluid fraction to nodal fluid fraction
    Teuchos::ParameterList params;
    params.set<int>("action",FLD::calc_fluidfrac_projection);

    // fill element fluid fraction value into params
    int eid = actele->Id();
    int lid = eleveccol->Map().LID(eid);
    if(lid < 0)
      dserror("element %i is not on this processor", eid);
    double eleval = (*eleveccol)[lid];
    params.set<double >("elefluidfrac", eleval);

    // call the element specific evaluate method (elevec1 = rhs, elemat1 = mass matrix)
    actele->Evaluate(params,*fluiddis_,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);

   // assembling into node maps
   matrix->Assemble(eid,elematrix1,lm,lmowner);
   LINALG::Assemble(*rhs,elevector1,lm,lmowner);
  } //end element loop

  // finalize the matrix
  matrix->Complete();

  // get solver parameter list of linear solver
  const int solvernumber = DRT::Problem::Instance()->CavitationParams().get<int>("FLUIDFRAC_PROJ_SOLVER");
  // solve system
  Teuchos::RCP<Epetra_Vector> nodevec = SolveOnNodeBasedVector(solvernumber, noderowmap, matrix, rhs);

  // dofrowmap does not include separate dofs for slave pbc nodes
  for(int i=0; i<noderowmap->NumMyElements(); ++i)
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
Teuchos::RCP<Epetra_Vector> CAVITATION::Algorithm::SolveOnNodeBasedVector(
  const int solvernumber,
  Teuchos::RCP<Epetra_Map> noderowmap,
  Teuchos::RCP<LINALG::SparseMatrix> matrix,
  Teuchos::RCP<Epetra_Vector> rhs
  )
{
  const Teuchos::ParameterList& solverparams = DRT::Problem::Instance()->SolverParams(solvernumber);

  Teuchos::RCP<LINALG::Solver> solver =
                                  Teuchos::rcp(new LINALG::Solver(solverparams,
                                  fluiddis_->Comm(),
                                  DRT::Problem::Instance()->ErrorFile()->Handle()));

  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams,"AZPREC");
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
    if(prectyp == INPAR::SOLVER::azprec_ML
        or prectyp == INPAR::SOLVER::azprec_MLfluid
        or prectyp == INPAR::SOLVER::azprec_MLAPI
        or prectyp == INPAR::SOLVER::azprec_MLfluid2)
      preclist_ptr = &((solver->Params()).sublist("ML Parameters"));
    else if(prectyp == INPAR::SOLVER::azprec_MueLuAMG_sym
        or prectyp == INPAR::SOLVER::azprec_MueLuAMG_nonsym)
      preclist_ptr = &((solver->Params()).sublist("MueLu Parameters"));
    else
      dserror("please add correct parameter list");

    Teuchos::ParameterList& preclist = *preclist_ptr;
    preclist.set<Teuchos::RCP<std::vector<double> > > ("nullspace",Teuchos::null);
    // ML would not tolerate this Teuchos::rcp-ptr in its list otherwise
    preclist.set<bool>("ML validate parameter list",false);

    preclist.set("PDE equations",1);
    preclist.set("null space: dimension",1);
    preclist.set("null space: type","pre-computed");
    preclist.set("null space: add default vectors",false);

    // allocate the local length of the rowmap
    const int lrows = noderowmap->NumMyElements();
    Teuchos::RCP<std::vector<double> > ns = Teuchos::rcp(new std::vector<double>(lrows));
    double* nullsp = &((*ns)[0]);

    // compute null space manually
    for (int j=0; j<lrows; ++j)
      nullsp[j] = 1.0;

    preclist.set<Teuchos::RCP<std::vector<double> > >("nullspace",ns);
    preclist.set("null space: vectors",nullsp);
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
