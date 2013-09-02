/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.cpp

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "particle_timint_centrdiff.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_calc.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_mortar/mortar_utils.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_cut/cut_position.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_cavitation.H"
#include "../drt_io/io_pstream.H"
#include "../headers/definitions.h"

#include <Epetra_FEVector.h>
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
CAVITATION::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : PARTICLE::Algorithm(comm,params),
  coupalgo_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::CouplingStrategyOverFields>(params,"COUPALGO")),
  void_frac_strategy_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::VoidFractionCalculation>(params,"VOID_FRACTION_CALC")),
  gauss_rule_per_dir_(params.get<int>("NUM_GP_VOID_FRACTION")),
  approxelecoordsinit_((bool)DRT::INPUT::IntegralValue<int>(params,"APPROX_ELECOORDS_INIT")),
  fluiddis_(Teuchos::null),
  fluid_(Teuchos::null),
  ele_volume_(Teuchos::null)
{
  // setup fluid time integrator
  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");
  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(),false));
  fluid_ = fluid->FluidFieldrcp();

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the cavitation algorithm                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // deal with particle inflow
    ParticleInflow();

    // transfer particles into their correct bins
    TransferParticles();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished

}


/*----------------------------------------------------------------------*
 | setup of the system                                      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupSystem()
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Init(bool restarted)
{
  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  CreateBins();

  // gather all fluid coleles in each bin for proper extended ghosting
  std::map<int, std::set<int> > fluideles;
  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcs(fluideles);

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // read out bubble inflow condition and set bubble inflows in corresponding bins
    // assumption: only row bins are available up to here
    BuildBubbleInflowCondition();
  }

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>,PARTICLE::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBins(homelessparticles);

  // ghost bins, particles and fluid elements according to the bins
  SetupGhosting(binrowmap, fluideles);

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // access structural dynamic params list which will be possibly modified while creating the time integrator
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

    // create time integrator based on structural time integration
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> particles =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(), const_cast<Teuchos::ParameterList&>(sdyn), particledis_));
    particles_ = particles->StructureFieldrcp();

    // determine consistent initial acceleration for the particles
    CalculateAndApplyForcesToParticles();
    Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->DetermineMassDampConsistAccel();
  }

  // compute volume of each fluid element and store it
  ele_volume_ = LINALG::CreateVector(*fluiddis_->ElementRowMap(), false);
  int numfluidele = fluiddis_->NumMyRowElements();
  for(int i=0; i<numfluidele; ++i)
  {
    DRT::Element* fluidele = fluiddis_->lRowElement(i);
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));
    double ev = GEO::ElementVolume( fluidele->Shape(), xyze );
    (*ele_volume_)[i] = ev;
  }

  // some output
  IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);
  DRT::UTILS::PrintParallelDistribution(*fluiddis_);

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrepareTimeStep()
{
  PARTICLE::Algorithm::PrepareTimeStep();
  fluid_->PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Integrate()
{
  {
    TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateVoidFraction");
    CalculateVoidFraction();
  }

  // apply forces and solve particle time step
  PARTICLE::Algorithm::Integrate();

  {
    TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::IntegrateFluid");
    fluid_->Solve();
  }

  return;
}


/*----------------------------------------------------------------------*
 | void fraction calculation                               ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateVoidFraction()
{
  Teuchos::RCP<Epetra_FEVector> void_volumes = Teuchos::rcp(new Epetra_FEVector(*fluiddis_->ElementRowMap()));

  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->ExtractDispn();
  Teuchos::RCP<Epetra_Vector> particleradius = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp();

  std::set<int> examinedbins;
  int numrownodes = particledis_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    DRT::Node *currentparticle = particledis_->lRowNode(i);
    DRT::Element** currentbin = currentparticle->Elements();

    int binId=currentbin[0]->Id();

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
    double minbin = bin_size_[0];
    for(int dim=1; dim<3; ++dim)
      minbin = std::min(minbin, bin_size_[dim]);

    // scaling factor in order to account for influence of bubble
    double scale = 1.2;
    int ibinrange = (int)((maxradius*scale)/minbin) + 1;
    int ijk_range[] = {ijk[0]-ibinrange, ijk[0]+ibinrange, ijk[1]-ibinrange, ijk[1]+ibinrange, ijk[2]-ibinrange, ijk[2]+ibinrange};

    if(ibinrange > 3)
      dserror("not yet tested for such large bubbles");

    // variable to store bin ids of surrounding bins
    std::set<int> binIds;

    // get corresponding bin ids in ijk range and fill them into binIds (in gid)
    GidsInijkRange(&ijk_range[0], binIds, true);

    // variable to store all fluid elements in neighborhood
    std::set<DRT::Element*> neighboringfluideles;
    for(std::set<int>::const_iterator i=binIds.begin(); i!=binIds.end(); ++i)
    {
      // extract bins from discretization
      DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>( particledis_->gElement(*i) );
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
      LINALG::Matrix<3,1> particleposition;
      std::vector<int> lm_b = particledis_->Dof(currparticle);
      // convert dof of particle into correct local id in this Epetra_Vector
      int posx = bubblepos->Map().LID(lm_b[0]);
      for (int dim=0; dim<3; ++dim)
      {
        // fill position of currentparticle
        particleposition(dim) = (*bubblepos)[posx+dim];
      }

      // get particle radius and influence of bubble
      double r_p = (*particleradius)[ particledis_->NodeRowMap()->LID(currparticle->Id()) ];
      double influence = scale*r_p;

      double vol_influence;

      // variable to store not yet normalized volume for each fluid element
      std::map<int, double> volumefraction;

      bool surfaceoverlap;

      // do while volume is not negative which is important for overlapping surfaces/points of bubble and fluid
      do
      {
        surfaceoverlap = false;
        vol_influence = 0.0;

        // loop over set neighboring fluideles
        for (std::set<DRT::Element*>::const_iterator ineighbor=neighboringfluideles.begin(); ineighbor!=neighboringfluideles.end(); ++ineighbor)
        {
          DRT::Element* ele = *ineighbor;
          int gid = ele->Id();

          // get bounding box of current element
          LINALG::Matrix<3,2> xaabb;
          for(int inode=0; inode<ele->NumNode(); ++inode)
          {
            const DRT::Node* node = ele->Nodes()[inode];
            LINALG::Matrix<3,1> position;

            // fill with nodal positions
            for (int dim=0; dim<3; ++dim)
            {
              position(dim) = node->X()[dim];
              if(inode==0)
              {
                xaabb(dim,0) = position(dim);
                xaabb(dim,1) = position(dim);
              }
              else
              {
                if(position(dim) < xaabb(dim,0))
                  xaabb(dim,0) = position(dim);
                if(position(dim) > xaabb(dim,1))
                  xaabb(dim,1) = position(dim);
              }
            }
          }


          double bubblesurface[] = {particleposition(0)+influence, particleposition(0)-influence,
                                    particleposition(1)+influence, particleposition(1)-influence,
                                    particleposition(2)+influence, particleposition(2)-influence};

          bool boundingbox = true;
          // test whether the bounding box of the fluid element touches the bubbleinfluence
          for(int dim=0; dim<3; ++dim)
          {
            if(xaabb(dim,0) - GEO::TOL7 > bubblesurface[dim*2] or xaabb(dim,1) + GEO::TOL7 < bubblesurface[dim*2+1])
            {
              boundingbox = false;
              break;
            }
          }

          double vol_ele = 0.0;
          if(boundingbox == true)
          {
            switch(void_frac_strategy_)
            {
            case INPAR::CAVITATION::analytical_constpoly:
            {
              DoAnalyticalIntegration(ele, particleposition, influence, vol_ele, surfaceoverlap);
            }
            break;
            case INPAR::CAVITATION::analytical_quadraticpoly:
            {
              DoAnalyticalIntegration(ele, particleposition, influence, vol_ele, surfaceoverlap);
            }
            break;
            case INPAR::CAVITATION::gaussian_integration:
            {
              DoGaussianIntegration(ele, particleposition, influence, vol_ele);
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

            // in case of polynomial influence, twisted surfaces and tight cut situations,
            // the polynomial can become negative very close to the bubble influence area
            if(vol_ele < 0.0)
              vol_ele = 0.0;

            // sum and store volume for each fluid element
            vol_influence += vol_ele;
            volumefraction[gid] = vol_ele;
          } // end if boundingbox

        } // end loop neighboring eles

      } while(surfaceoverlap == true);

      // normalize with volume of bubble divided by volume of influence
      double normalization = 4.0/3.0*M_PI*pow(r_p,3.0) / vol_influence;
      // assemble void volume of fluid elements
      for(std::map<int, double>::const_iterator iter=volumefraction.begin(); iter!=volumefraction.end(); ++iter)
      {
        double val = iter->second * normalization;

        // do assembly of void fraction into fluid
        int err = void_volumes->SumIntoGlobalValues(1, &(iter->first), &val);
        if (err<0)
          dserror("summing into Epetra_FEVector failed");
      }

    } // loop over particles in one bin
  } // loop over outer particles

  // call global assemble
  int err = void_volumes->GlobalAssemble(Add, false);
  if (err<0)
    dserror("global assemble into fluidforces failed");

  // divide element wise void volume by element volume
  void_volumes->ReciprocalMultiply(1.0, *ele_volume_, *void_volumes, 0.0);

  // apply void fraction to fluid
  fluid_->SetVoidVolume(void_volumes);

  return;
}


/*----------------------------------------------------------------------*
 | Gaussian integration for void fraction calculation      ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::DoGaussianIntegration(
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

  // set action in order to calculate the velocity and material derivative of the velocity
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::void_fraction_gaussian_integration);

  params.set<LINALG::Matrix<3,1> >("particle_pos", particleposition);

  params.set<double>("influence", influence);

  params.set<int>("gp_per_dir", gauss_rule_per_dir_);

  // call the element specific evaluate method (elevec1 = void fraction)
  ele->Evaluate(params,*fluiddis_,lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

  vol_ele = elevector1[0];

  return;
}


/*----------------------------------------------------------------------*
 | analytic integration for void fraction calculation      ghamm 08/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::DoAnalyticalIntegration(
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
    r.Update(1,surfacenodes[0],-1,surfacenodes[2]);
    p.Update(1,surfacenodes[1],-1,surfacenodes[3%numsurfacenodes]);

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
 | calculate fluid forces on particle and apply it         ghamm 01/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateAndApplyForcesToParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::CalculateAndApplyForcesToParticles");
  const int dim = 3;

  fluiddis_->ClearState();
  particledis_->ClearState();

  //TODO 1: Test whether these states are updated enough
  // at the beginning of the time step: veln := velnp(previous step)
  fluiddis_->SetState("veln",fluid_->Veln());
  fluiddis_->SetState("velnm",fluid_->Velnm());

  Teuchos::RCP<const Epetra_Vector> bubblepos = particles_->ExtractDispn();
  Teuchos::RCP<const Epetra_Vector> bubblevel = particles_->ExtractVeln();
  Teuchos::RCP<Epetra_Vector> bubbleradius = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp();

  // vectors to be filled with forces, note: global assemble is needed for fluidforces due to the case with large bins and small fluid eles
  Teuchos::RCP<Epetra_Vector> bubbleforces = LINALG::CreateVector(*particledis_->DofRowMap(),true);
  Teuchos::RCP<Epetra_FEVector> fluidforces = Teuchos::rcp(new Epetra_FEVector(*fluiddis_->DofRowMap()));

  if(fluiddis_->NumMyColElements() <= 0)
    dserror("there is no fluid element to ask for material parameters");
  Teuchos::RCP<MAT::NewtonianFluid> actmat = Teuchos::rcp_dynamic_cast<MAT::NewtonianFluid>(fluiddis_->lColElement(0)->Material());
  if(actmat == Teuchos::null)
    dserror("type cast of fluid material failed");

  double rho_l = actmat->Density();
  double mu_l = actmat->Viscosity();
  double rho_b = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ParticleDensity();

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // only row particles are evaluated
  for(int i=0; i<particledis_->NodeRowMap()->NumMyElements(); ++i)
  {
    DRT::Node* currparticle = particledis_->lRowNode(i);
    // fill particle position
    LINALG::Matrix<3,1> particleposition;
    std::vector<int> lm_b = particledis_->Dof(currparticle);
    int posx = bubblepos->Map().LID(lm_b[0]);
    for (int d=0; d<dim; ++d)
      particleposition(d) = (*bubblepos)[posx+d];


    //--------------------------------------------------------------------
    // 1st step: element coordinates of particle position in fluid element
    //--------------------------------------------------------------------

    // variables to store information about element in which the particle is located
    DRT::Element* targetfluidele = NULL;
    LINALG::Matrix<3,1> elecoord(true);

    // find out in which fluid element the current particle is located
    if(currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
    DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    DRT::Element** fluidelesinbin = currbin->AssociatedFluidEles();
    int numfluidelesinbin = currbin->NumAssociatedFluidEle();

    std::set<int>::const_iterator eleiter;
    // search for underlying fluid element with fast search if desired
    for(int ele=0; ele<numfluidelesinbin; ++ele)
    {
      DRT::Element* fluidele = fluidelesinbin[ele];
      const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

      //get coordinates of the particle position in parameter space of the element
      bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, elecoord, approxelecoordsinit_);

      if(insideele == true)
      {
        targetfluidele = fluidele;
        // leave loop over all fluid eles in bin
        break;
      }
    }

    // repeat search for underlying fluid element with standard search in case nothing was found
    if(targetfluidele == NULL and approxelecoordsinit_ == true)
    {
      for(int ele=0; ele<numfluidelesinbin; ++ele)
      {
        DRT::Element* fluidele = fluidelesinbin[ele];
        const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(fluidele));

        //get coordinates of the particle position in parameter space of the element
        bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, elecoord, false);

        if(insideele == true)
        {
          targetfluidele = fluidele;
          // leave loop over all fluid eles in bin
          break;
        }
      }
    }


    //--------------------------------------------------------------------
    // 2nd step: forces on this bubble are calculated
    //--------------------------------------------------------------------

    if(targetfluidele == NULL)
    {
      std::cout << "INFO: currparticle with Id: " << currparticle->Id() << " and position: " << particleposition(0) << " "
          << particleposition(1) << " " << particleposition(2) << " " << " does not have an underlying fluid element -> no forces calculated" << std::endl;

      std::vector<double> tmpposition(dim);
      for(int d=0; d<dim; ++d)
        tmpposition[d] = particleposition(d);
      int bubbleBinId = ConvertPosToGid(tmpposition);
      std::cout << "particle is in binId: " << bubbleBinId << " while currbin->Id() is " << currbin->Id() <<
          " . The following number of fluid eles is in this bin:" << numfluidelesinbin << std::endl;

      // do not assemble forces for this bubble and continue with next bubble
      continue;
    }

    // get element location vector and ownerships
    std::vector<int> lm_f;
    std::vector<int> lmowner_f;
    std::vector<int> lmstride;
    targetfluidele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(dim);
    elevector2.Size(dim);
    elevector3.Size(dim);

    // set action in order to calculate the velocity and material derivative of the velocity
    Teuchos::ParameterList params;
    params.set<int>("action",FLD::calc_mat_deriv_u_and_rot_u);
    params.set<double>("timestep",Dt());
    params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

    // call the element specific evaluate method (elevec1 = fluid vel u; elevec2 = mat deriv of fluid vel, elevec3 = rot of fluid vel)
    targetfluidele->Evaluate(params,*fluiddis_,lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    // get bubble velocity and acceleration
    std::vector<double> v_bub(lm_b.size());
    DRT::UTILS::ExtractMyValues(*bubblevel,v_bub,lm_b);

    // get bubble radius
    double r_bub = (*bubbleradius)[ particledis_->NodeRowMap()->LID(currparticle->Id()) ];

    // bubble Reynolds number
    LINALG::Matrix<3,1> v_rel(false);
    for (int d=0; d<dim; ++d)
      v_rel(d) = elevector1[d] - v_bub[d];

    double v_relabs = v_rel.Norm2();
    double Re_b = 2.0 * r_bub * v_relabs * rho_l / mu_l;

    bool output = false;
    if(output)
    {
      std::cout << "v_rel: " << v_rel(0) << "  " << v_rel(1) << "  " << v_rel(2) << "  " << std::endl;
      std::cout << "v_relabs: " << v_relabs << std::endl;
      std::cout << "bubble Reynolds number: " << Re_b << std::endl;
    }

    // variable to sum forces for the current bubble under observation
    LINALG::Matrix<3,1> sumforces;
    /*------------------------------------------------------------------*/
    //// 2.1) drag force = 0.5 * c_d * rho_l * Pi * r_b^2 * |u-v| * (u-v) or
    //// Stokes law for very small Re: drag force = 6.0 * Pi * mu_l * r_b * (u-v)
    double coeff1 = 0.0;
    if(Re_b < 0.1)
    {
      coeff1 = 6.0 * M_PI * mu_l * r_bub;
    }
    else
    {
      double c_d = 0.0;
      if(Re_b < 1000.0)
        c_d = 24.0 * (1.0 + 0.15 * pow(Re_b,0.687)) / Re_b;
      else
        c_d = 0.44;

      coeff1 = 0.5 * c_d * rho_l * M_PI * r_bub * r_bub * v_relabs;
    }

    LINALG::Matrix<3,1> dragforce(v_rel);
    dragforce.Scale(coeff1);
    //assemble
    sumforces.Update(1.0, dragforce);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.2) lift force = c_l * rho_l * volume_b * (u-v) x rot_u   with rot_u = nabla x u
    double c_l = 0.5;
    double vol_b = 4.0 / 3.0 * M_PI * r_bub * r_bub* r_bub;
    LINALG::Matrix<3,1> rot_u;
    for (int d=0; d<dim; ++d)
      rot_u(d) = elevector3(d);

    LINALG::Matrix<3,1> liftforce = GEO::computeCrossProduct(v_rel, rot_u);

    double coeff2 = c_l * rho_l * vol_b;
    liftforce.Scale(coeff2);
    //assemble
    sumforces.Update(1.0, liftforce, 1.0);
    // store forces for coupling to fluid
    LINALG::Matrix<3,1> couplingforce(sumforces);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.3) gravity and buoyancy forces = volume_b * rho_bub * g - volume_b * rho_l * ( g - Du/Dt )
    LINALG::Matrix<3,1> Du_Dt;
    for (int d=0; d<dim; ++d)
      Du_Dt(d) = elevector2[d];

    LINALG::Matrix<3,1> grav_buoy_force;
    grav_buoy_force.Update(rho_b, gravity_acc_);
    grav_buoy_force.Update(-rho_l, gravity_acc_, rho_l, Du_Dt, 1.0);
    grav_buoy_force.Scale(vol_b);
    //assemble
    sumforces.Update(1.0, grav_buoy_force, 1.0);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.4) virtual/added mass = c_VM * rho_l * volume_b * ( Du/Dt - Dv/Dt )
    //// Note: implicit treatment of bubble acceleration in added mass, other forces explicit
    //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) + c_VM * rho_l * volume_b * Du/Dt }{ 1 + c_VM * rho_l / rho_b }
    double c_VM = 0.5;
    double coeff3 = c_VM * rho_l * vol_b;
    double coeff4 = 1.0 + c_VM * rho_l / rho_b;

    LINALG::Matrix<3,1> bubbleforce;
    bubbleforce.Update(1.0/coeff4, sumforces, coeff3/coeff4, Du_Dt);
    /*------------------------------------------------------------------*/


    //--------------------------------------------------------------------
    // 3rd step: assemble bubble forces
    //--------------------------------------------------------------------

    // assemble of bubble forces (note: row nodes evaluated)
    Epetra_SerialDenseVector forcecurrbubble(3);
    for(int d=0; d<dim; ++d)
      forcecurrbubble[d] = bubbleforce(d);
    std::vector<int> lmowner_b(lm_b.size(), myrank_);
    LINALG::Assemble(*bubbleforces,forcecurrbubble,lm_b,lmowner_b);

    // coupling forces between fluid and particle only include certain forces
    // calculate added mass force
    LINALG::Matrix<3,1> addedmassforce;
    double m_b = vol_b * rho_b;
    addedmassforce.Update(coeff3, Du_Dt, -coeff3/m_b, bubbleforce);
    //// coupling force = -(dragforce + liftforce + addedmassforce); actio = reactio --> minus sign
    couplingforce.Update(-1.0, addedmassforce, -1.0);

    if(coupalgo_ != INPAR::CAVITATION::OneWay)
    {
      // assemble of fluid forces must be done globally because col entries in the fluid can occur
      // although only row particles are evaluated
      const int numnode = targetfluidele->NumNode();
      Epetra_SerialDenseVector funct(numnode);
      // get shape functions of the element; evaluated at the bubble position --> distribution
      DRT::UTILS::shape_function_3D(funct,elecoord(0,0),elecoord(1,0),elecoord(2,0),targetfluidele->Shape());
      // prepare assembly for fluid forces (pressure degrees do not have to be filled)

      int numdofperfluidele = numnode*(dim+1);
      double val[numdofperfluidele];
      for(int iter=0; iter<numnode; ++iter)
      {
        for(int d=0; d<dim; ++d)
        {
          val[iter*(dim+1) + d] = funct[iter] * couplingforce(d);
        }
      }
      // do assembly of bubble forces on fluid
      int err = fluidforces->SumIntoGlobalValues(numdofperfluidele, &lm_f[0], &val[0]);
      if (err<0)
        dserror("summing into Epetra_FEVector failed");
    }


    //--------------------------------------------------------------------
    // 4th step: output
    //--------------------------------------------------------------------
    if(output)
    {
      // gravity
      LINALG::Matrix<3,1> gravityforce(gravity_acc_);
      gravityforce.Scale(rho_b*vol_b);
      std::cout << "gravity force       : " << gravityforce << std::endl;

      // buoyancy
      double coeff5 = - vol_b * rho_l;
      LINALG::Matrix<3,1> buoyancyforce(true);
      buoyancyforce.Update(coeff5, gravity_acc_);
      std::cout << "buoyancy force      : " << buoyancyforce << std::endl;

      // effective buoyancy / inertia term
      LINALG::Matrix<3,1> effectbuoyancyforce;
      effectbuoyancyforce.Update(-coeff5, Du_Dt);
      std::cout << "effective buoy force: " << effectbuoyancyforce << std::endl;

      // drag, lift and added mass force
      std::cout << "dragforce force     : " << dragforce << std::endl;
      std::cout << "liftforce force     : " << liftforce << std::endl;
      std::cout << "added mass force    : " << addedmassforce << std::endl;

      // sum over all bubble forces
      std::cout << "sum over all forces : " << bubbleforce << std::endl;

      // fluid force
      std::cout << "fluid force         : " << couplingforce << std::endl;
    }

  } // end iparticle

  //--------------------------------------------------------------------
  // 4th step: apply forces to bubbles and fluid field
  //--------------------------------------------------------------------
  particledis_->SetState("particleforces", bubbleforces);

  if(coupalgo_ != INPAR::CAVITATION::OneWay)
  {
    // call global assemble
    int err = fluidforces->GlobalAssemble(Add, false);
    if (err<0)
      dserror("global assemble into fluidforces failed");

    fluid_->ApplyExternalForces(fluidforces);
  }

  return;
}


/*----------------------------------------------------------------------*
 | particles are inserted into domain                      ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ParticleInflow()
{
  std::map<int, std::list<Teuchos::RCP<BubbleSource> > >::const_iterator biniter;

  int timeforinflow = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    // all particles have the same inflow frequency --> it is enough to test one
    // assumption only valid in case of one condition or conditions with identical inflow frequency
    if(biniter->second.size() != 0)
    {
      double inflowtime = 1.0 / biniter->second.front()->inflow_freq_;
      if(Step() % ((int)(inflowtime/Dt())) == 0)
      {
        timeforinflow = 1;
        break;
      }
    }
  }

  int globaltimeforinflow = 0;
  particledis_->Comm().MaxAll(&timeforinflow, &globaltimeforinflow, 1);
  if(globaltimeforinflow == 0)
    return; // no inflow detected


  // initialize bubble id with largest bubble id in use + 1 (on each proc)
  int maxbubbleid = particledis_->NodeRowMap()->MaxAllGID()+1;

  // start filling particles
  int inflowcounter = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      std::vector<double> inflow_position = (*particleiter)->inflow_position_;
      std::set<Teuchos::RCP<DRT::Node>,PARTICLE::Less> homelessparticles;
      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;
      Teuchos::RCP<DRT::Node> newparticle = Teuchos::rcp(new DRT::Node(newbubbleid, &inflow_position[0], myrank_));
      PlaceNodeCorrectly(newparticle, &inflow_position[0], homelessparticles);
      if(homelessparticles.size() != 0)
        dserror("New bubble could not be inserted on this proc! Bubble inflow broken.");
    }
    inflowcounter += (int)biniter->second.size();
  }

  std::cout << "Inflow of " << inflowcounter << " bubbles on proc " << myrank_ << std::endl;

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  particledis_->FillComplete(true, false, true);

  // update of state vectors to the new maps
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();

  // insert data for new bubbles into state vectors
  const Epetra_Map* dofrowmap = particledis_->DofRowMap();
  const Epetra_Map* noderowmap = particledis_->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> disn = particles_->ExtractDispnp();
  Teuchos::RCP<Epetra_Vector> veln = particles_->ExtractVelnp();
  Teuchos::RCP<Epetra_Vector> radiusn = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp();

  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      std::vector<double> inflow_position = (*particleiter)->inflow_position_;
      std::vector<double> inflow_vel = (*particleiter)->inflow_vel_;
      int inflow_vel_curve = (*particleiter)->inflow_vel_curve_;
      double inflow_radius = (*particleiter)->inflow_radius_;
      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;

      double curvefac = 1.0;
      // curves are numbered starting with 1 in the input file
      if(inflow_vel_curve > 0)
        curvefac = DRT::Problem::Instance()->Curve(inflow_vel_curve-1).f(Time());

      DRT::Node* currparticle = particledis_->gNode(newbubbleid);
      // get the first gid of a particle and convert it into a LID
      int lid = dofrowmap->LID(particledis_->Dof(currparticle, 0));
      for(int dim=0; dim<3; ++dim)
      {
        (*disn)[lid+dim] = inflow_position[dim];
        (*veln)[lid+dim] = inflow_vel[dim] * curvefac;
      }
      lid = noderowmap->LID(newbubbleid);
      (*radiusn)[lid] = inflow_radius;
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Update()
{
  PARTICLE::Algorithm::Update();
  fluid_->Update();
  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ReadRestart(int restart)
{
  PARTICLE::Algorithm::ReadRestart(restart);
  fluid_->ReadRestart(restart);
  return;
}


/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CreateBins()
{
  // if not yet specified, get XAABB_ from underlying discretization
  if( XAABB_(2,1) > 0.9e12  and  XAABB_(2,1) < 1.1e12 )
  {
    IO::cout << "XAABB is computed based on the underlying fluid discretization" << IO::endl;
    XAABB_ = GEO::getXAABBofDis(*fluiddis_);
    // local bounding box
    double locmin[3] = {XAABB_(0,0), XAABB_(1,0), XAABB_(2,0)};
    double locmax[3] = {XAABB_(0,1), XAABB_(1,1), XAABB_(2,1)};
    // global bounding box
    double globmin[3];
    double globmax[3];
    // do the necessary communication
    Comm().MinAll(&locmin[0], &globmin[0], 3);
    Comm().MaxAll(&locmax[0], &globmax[0], 3);

    for(int dim=0; dim<3; ++dim)
    {
      XAABB_(dim,0) = globmin[dim];
      XAABB_(dim,1) = globmax[dim];
    }
  }

  // divide global bounding box into bins
  for (int dim = 0; dim < 3; ++dim)
  {
    // std::floor leads to bins that are at least of size cutoff_radius
    bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));
    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
  }

  IO::cout << "Global bounding box size: " << XAABB_;
  IO::cout << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << IO::endl;

  return;
}


/*----------------------------------------------------------------------*
| bins are distributed to the processors                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CAVITATION::Algorithm::DistributeBinsToProcs(std::map<int, std::set<int> >& fluideles)
{
  //--------------------------------------------------------------------
  // 1st step: exploiting bounding box idea for fluid elements and bins
  //--------------------------------------------------------------------
  {
    TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::DistributeBinsToProcs_step1");
    // loop over row fluid elements is enough, extended ghosting always includes standard ghosting
    for (int lid = 0; lid < fluiddis_->NumMyRowElements(); ++lid)
    {
      DRT::Element* fluidele = fluiddis_->lRowElement(lid);
      DRT::Node** fluidnodes = fluidele->Nodes();
      const int numnode = fluidele->NumNode();

      // initialize ijk_range with ijk of first node of fluid element
      int ijk[3];
      {
        const DRT::Node* node = fluidnodes[0];
        const double* coords = node->X();
        ConvertPosToijk(coords, ijk);
      }

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; ++j)
      {
        const DRT::Node* node = fluidnodes[j];
        const double* coords = node->X();
        int ijk[3];
        ConvertPosToijk(coords, ijk);

        for(int dim=0; dim<3; ++dim)
        {
          if(ijk[dim]<ijk_range[dim*2])
            ijk_range[dim*2]=ijk[dim];
          if(ijk[dim]>ijk_range[dim*2+1])
            ijk_range[dim*2+1]=ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range
      std::set<int> binIds;
      GidsInijkRange(&ijk_range[0], binIds, false);

      // assign fluid element to bins
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
        fluideles[*biniter].insert(fluidele->Id());
    }
  }


  //--------------------------------------------------------------------
  // 2nd step: decide which proc will be owner of each bin
  //--------------------------------------------------------------------

  std::vector<int> rowbins;
  {
    TEUCHOS_FUNC_TIME_MONITOR("CAVITATION::Algorithm::DistributeBinsToProcs_step2");
    // NOTE: This part of the setup can be the bottleneck because vectors of all bins
    // are needed on each proc (memory issue!!); std::map could perhaps help when gathering
    // num fluid nodes in each bin, then block wise communication after copying data to vector

    int numbins = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];
    std::vector<int> mynumeles_per_bin(numbins,0);

    std::map<int, std::set<int> >::const_iterator iter;
    for(iter=fluideles.begin(); iter!=fluideles.end(); ++ iter)
    {
      mynumeles_per_bin[iter->first] = iter->second.size();
    }

    // find maximum number of eles in each bin over all procs (init with -1)
    std::vector<int> maxnumeles_per_bin(numbins,-1);
    fluiddis_->Comm().MaxAll(&mynumeles_per_bin[0], &maxnumeles_per_bin[0], numbins);

    // it is possible that several procs have the same number of eles in a bin
    // only proc which has maximum number of eles in a bin writes its rank
    std::vector<int> myrank_per_bin(numbins,-1);
    for(int i=0; i<numbins; ++i)
    {
      if(mynumeles_per_bin[i] == maxnumeles_per_bin[i])
        myrank_per_bin[i] = myrank_;
    }

    mynumeles_per_bin.clear();
    maxnumeles_per_bin.clear();

    // find maximum myrank for each bin over all procs (init with -1)
    std::vector<int> maxmyrank_per_bin(numbins,-1);
    fluiddis_->Comm().MaxAll(&myrank_per_bin[0], &maxmyrank_per_bin[0], numbins);

    // distribute bins to proc with highest rank
    for(int gid=0; gid<numbins; ++gid)
    {
      if(myrank_ == maxmyrank_per_bin[gid])
      {
        Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy", gid, myrank_);
        particledis_->AddElement(bin);
        rowbins.push_back(gid);
      }
    }

    myrank_per_bin.clear();
    maxmyrank_per_bin.clear();
  }

  // return binrowmap (without having called FillComplete on particledis_ so far)
  return Teuchos::rcp(new Epetra_Map(-1,(int)rowbins.size(),&rowbins[0],0,Comm()));
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap, std::map<int, std::set<int> >& fluideles)
{
  //--------------------------------------------------------------------
  // 1st and 2nd step
  //--------------------------------------------------------------------

  PARTICLE::Algorithm::SetupGhosting(binrowmap);


  //--------------------------------------------------------------------
  // 3st step: extend ghosting of underlying fluid discretization according to bin distribution
  //--------------------------------------------------------------------
  std::map<int, std::set<int> > extendedfluidghosting;
  {
    // do communication to gather all elements for extended ghosting
    const int numproc = fluiddis_->Comm().NumProc();

    for (int iproc = 0; iproc < numproc; ++iproc)
    {
      // first: proc i tells all procs how many col bins it has
      int numbin = bincolmap_->NumMyElements();
      fluiddis_->Comm().Broadcast(&numbin, 1, iproc);
      // second: proc i tells all procs which col bins it has
      std::vector<int> binid(numbin,0);
      if(iproc == myrank_)
      {
        int* bincolmap = bincolmap_->MyGlobalElements();
        for (int i=0; i<numbin; ++i)
          binid[i] = bincolmap[i];
      }
      fluiddis_->Comm().Broadcast(&binid[0], numbin, iproc);

      // loop over all own bins and find requested ones
      std::map<int, std::set<int> > sdata;
      std::map<int, std::set<int> > rdata;

      for(int i=0; i<numbin; ++i)
      {
        sdata[binid[i]].insert(fluideles[binid[i]].begin(),fluideles[binid[i]].end());
      }

      LINALG::Gather<int>(sdata, rdata, 1, &iproc, fluiddis_->Comm());

      // proc i has to store the received data
      if(iproc == myrank_)
      {
        extendedfluidghosting = rdata;
      }
    }

    //reduce map of sets to one set and copy to a vector to create fluidcolmap
    std::set<int> redufluideleset;
    std::map<int, std::set<int> >::iterator iter;
    for(iter=extendedfluidghosting.begin(); iter!= extendedfluidghosting.end(); ++iter)
    {
      redufluideleset.insert(iter->second.begin(),iter->second.end());
    }
    std::vector<int> fluidcolgids(redufluideleset.begin(),redufluideleset.end());
    Teuchos::RCP<Epetra_Map> fluidcolmap = Teuchos::rcp(new Epetra_Map(-1,(int)fluidcolgids.size(),&fluidcolgids[0],0,Comm()));

    // create ghosting for fluid eles (each knowing its node ids)
    fluiddis_->ExportColumnElements(*fluidcolmap);

    // create a set of node IDs for each proc (row + ghost)
    std::set<int> nodes;
    for (int lid=0;lid<fluidcolmap->NumMyElements();++lid)
    {
      DRT::Element* ele = fluiddis_->gElement(fluidcolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0; inode<ele->NumNode(); ++inode)
        nodes.insert(nodeids[inode]);
    }

    // copy nodegids to a vector and create nodecolmap
    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,Comm()));

    // create ghosting for nodes
    fluiddis_->ExportColumnNodes(*nodecolmap);

    // do a final fillcomplete to build connectivity
    fluiddis_->FillComplete(true,true,true);

  }

  //--------------------------------------------------------------------
  // 4th step: assign fluid elements to bins
  //--------------------------------------------------------------------
  {
    for(std::map<int, std::set<int> >::const_iterator biniter=extendedfluidghosting.begin(); biniter!=extendedfluidghosting.end(); ++biniter)
    {
      DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->gElement(biniter->first));
      for(std::set<int>::const_iterator fluideleiter=biniter->second.begin(); fluideleiter!=biniter->second.end(); ++fluideleiter)
      {
        int fluideleid = *fluideleiter;
        currbin->AddAssociatedFluidEle(fluideleid, fluiddis_->gElement(fluideleid));
//          cout << "in bin with id:" << currbin->Id() << " is fluid ele with id" << fluideleid << "with pointer" << fluiddis_->gElement(fluideleid) << endl;
      }
    }
  }

#ifdef DEBUG
  // check whether each particle has an underlying fluid element
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for(int i=0; i<fluiddis_->NumMyColNodes(); i++)
  {
    DRT::Node* node = fluiddis_->lColNode(i);
    LINALG::Matrix<3,1> currpos;

    for (int a=0; a<3; a++)
    {
      currpos(a) = node->X()[a];
    }
    currentpositions.insert(std::pair<int,LINALG::Matrix<3,1> >(node->Id(),currpos));
  }
  // start loop over all particles
  for(int k=0; k<particledis_->NumMyColNodes(); k++)
  {
    DRT::Node* particle = particledis_->lColNode(k);
    const double* pos = particle->X();
    LINALG::Matrix<3,1> projpoint;
    for(int dim=0; dim<3; dim++)
      projpoint(dim) = pos[dim];
    bool foundele = false;
    for(int i=0; i<fluiddis_->NumMyColElements(); i++)
    {
      DRT::Element* fluidele = fluiddis_->lColElement(i);

      LINALG::Matrix<3,1> elecoord(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(Teuchos::rcp(fluidele,false), currentpositions));

      //get coordinates of the particle position in parameter space of the element
      foundele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, projpoint, elecoord);

      if(foundele == true)
        break;
    }
    if(foundele == false)
      dserror("particle (Id:%d) was found which does not have fluid support", particle->Id());
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
| build connectivity from fluid elements to bins            ghamm 07/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildElementToBinPointers()
{
  // first call base class to associate potential particle walls
  PARTICLE::Algorithm::BuildElementToBinPointers();

  // loop over column bins and fill fluid elements
  const int numcolbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numcolbin; ++ibin)
  {
    DRT::Element* actele = particledis_->lColElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
    const int numfluidele = actbin->NumAssociatedFluidEle();
    const int* fluideleids = actbin->AssociatedFluidEleIds();
    std::vector<DRT::Element*> fluidelements(numfluidele);
    for(int iele=0; iele<numfluidele; ++iele)
    {
      const int fluideleid = fluideleids[iele];
      fluidelements[iele] = fluiddis_->gElement(fluideleid);
    }
    actbin->BuildFluidElePointers(&fluidelements[0]);
  }

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fluid_->CreateFieldTest());
  PARTICLE::Algorithm::TestResults(comm);
  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Output()
{
  fluid_->Output();
  PARTICLE::Algorithm::Output();
  return;
}


/*----------------------------------------------------------------------*
 | get adjacent bins to corner, where ijk is in 1st octant ghamm 02/13  |
 *----------------------------------------------------------------------*/
std::vector<int> CAVITATION::Algorithm::AdjacentBinstoCorner(int* ijk)
{
  std::vector<int> adjbins;
  adjbins.reserve(8);

  // get all adjacent bins to the current corner, including the bin itself
  for(int i=-1;i<1;i++)
  {
    for(int j=-1;j<1;j++)
    {
      for(int k=-1;k<1;k++)
      {
        int ijk_neighbor[3] = {ijk[0]+i, ijk[1]+j, ijk[2]+k};

        int neighborgid = ConvertijkToGid(&ijk_neighbor[0]);
        if(neighborgid != -1)
        {
          adjbins.push_back(neighborgid);
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return adjbins;
}


/*----------------------------------------------------------------------*
| setup of bubble inflow                                    ghamm 01/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildBubbleInflowCondition()
{
  // build inflow boundary condition
  std::vector<DRT::Condition*> conds;
  particledis_->GetCondition("ParticleInflow", conds);
  // unique bubbleinflow id over all inflow conditions
  int bubbleinflowid = 0;
  for (size_t i=0; i<conds.size(); ++i)
  {
    if(i>0)
      dserror("only taken care of one particle inflow condition so far. "
          "Remedy: bubble_source_ needs to be a vector of the current layout");
    /*
     * inflow condition --> bubble sources
     *
     *  example: num_per_dir = {4, 5, 1}
     *
     *       <-> (dist_x = (vertex2_x-vertex1_x)/(num_per_dir_x-1))
     *
     *   x  x  x  x<-------- vertex2
     *
     *   x  x  x  x
     *
     *   x  x  x  x   ^
     *                | (dist_y = (vertex2_y-vertex1_y)/(num_per_dir_y-1) )
     *   x  x  x  x   ^
     *
     *   x  x  x  x
     *   ^
     *   |
     * vertex1
     *
     */

    // extract data from inflow condition
    const std::vector<double>* vertex1 = conds[i]->Get<std::vector<double> >("vertex1");
    const std::vector<double>* vertex2 = conds[i]->Get<std::vector<double> >("vertex2");
    const std::vector<int>* num_per_dir = conds[i]->Get<std::vector<int> >("num_per_dir");
    const std::vector<double>* inflow_vel = conds[i]->Get<std::vector<double> >("inflow_vel");
    int inflow_vel_curve = conds[i]->GetInt("inflow_vel_curve");
    double inflow_freq = conds[i]->GetDouble("inflow_freq");

    // make sure that a particle material is defined in the dat-file
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
    if (id==-1)
      dserror("Could not find particle material");

    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
    double initial_radius = actmat->initialradius_;

    double inflowtime = 1.0 / inflow_freq;
    if(std::abs(inflowtime/Dt() - (int)(inflowtime/Dt())) > EPS9)
      dserror("1/inflow_freq with inflow_freq = %f cannot be divided by fluid time step %f", inflowtime, Dt());
/* MUST BE ADDED WHEN PARTICLE CONTACT IS CONSIDERED
    double inflow_vel_mag = sqrt((*inflow_vel)[0]*(*inflow_vel)[0] + (*inflow_vel)[1]*(*inflow_vel)[1] + (*inflow_vel)[2]*(*inflow_vel)[2]);
    if(initial_radius/inflow_vel_mag > inflowtime)
      dserror("Overlap for inflowing bubbles expected: initial_radius/inflow_vel_mag = %f s > inflow_freq = %f s", initial_radius/inflow_vel_mag, inflowtime);
*/
    // loop over all bubble inflow positions and fill them into bin when they are on this proc;
    // up to here, only row bins are available
    std::vector<double> source_pos(3);
    for(int z=0; z<(*num_per_dir)[2]; ++z)
    {
      double dist_z = ((*vertex2)[2] - (*vertex1)[2]) / ((((*num_per_dir)[2]-1)!=0) ? ((*num_per_dir)[2]-1) : 1);
      source_pos[2] = (*vertex1)[2] + z * dist_z;
      for(int y=0; y<(*num_per_dir)[1]; ++y)
      {
        double dist_y = ((*vertex2)[1] - (*vertex1)[1]) / ((((*num_per_dir)[1]-1)!=0) ? ((*num_per_dir)[1]-1) : 1);
        source_pos[1] = (*vertex1)[1] + y * dist_y;
        for(int x=0; x<(*num_per_dir)[0]; ++x)
        {
          double dist_x = ((*vertex2)[0] - (*vertex1)[0]) / ((((*num_per_dir)[0]-1)!=0) ? ((*num_per_dir)[0]-1) : 1);
          source_pos[0] = (*vertex1)[0] + x * dist_x;
          // check whether this source position is on this proc
          int binId = ConvertPosToGid(source_pos);
          bool found = particledis_->HaveGlobalElement(binId);
          if(found == true)
          {
            Teuchos::RCP<BubbleSource> bubbleinflow = Teuchos::rcp(new BubbleSource(
                                                                          bubbleinflowid,
                                                                          source_pos,
                                                                          *inflow_vel,
                                                                          inflow_vel_curve,
                                                                          initial_radius,
                                                                          inflow_freq));
            bubble_source_[binId].push_back(bubbleinflow);
#ifdef DEBUG
            if(particledis_->gElement(binId)->Owner() != myrank_)
              dserror("Only row bins should show up here. Either add additional if-case or move ghosting to a later point in time.");
#endif
          }
          bubbleinflowid++;
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | particle source                                         ghamm 02/13  |
 *----------------------------------------------------------------------*/
CAVITATION::BubbleSource::BubbleSource(
  int bubbleinflowid,
  std::vector<double> inflow_position,
  std::vector<double> inflow_vel,
  int inflow_vel_curve,
  double inflow_radius,
  double inflow_freq
  ) :
  inflowid_(bubbleinflowid),
  inflow_position_(inflow_position),
  inflow_vel_(inflow_vel),
  inflow_vel_curve_(inflow_vel_curve),
  inflow_radius_(inflow_radius),
  inflow_freq_(inflow_freq)
{
}
