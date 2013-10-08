/*----------------------------------------------------------------------*/
/*!
\file particle_contact.cpp
\brief Particle collision handling

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_contact.H"
#include "particle_algorithm.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_coordtrafo.H"

#include "../drt_mat/stvenantkirchhoff.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_FEVector.h>


PARTICLE::ParticleCollisionHandler::ParticleCollisionHandler(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm
  ) :
  myrank_(discret->Comm().MyPID()),
  discret_(discret),
  particle_algorithm_(particlealgorithm),
  contact_energy_(0.0),
  g_max_(0.0),
  radiusncol_(Teuchos::null),
  masscol_(Teuchos::null),
  disncol_(Teuchos::null),
  velncol_(Teuchos::null),
  ang_velncol_(Teuchos::null)
{
  // make sure that a particle material is defined in the dat-file
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
  if (id==-1)
    dserror("Could not find particle material");

  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
  // currently all particles have identical density and radius
  double density = actmat->density_;
  nue_ = actmat->poissonratio_;
  young_ = actmat->young_;

  ReadContactParameters(density);

  if(myrank_ == 0)
  {
    FILE *stream = 0;
    stream = fopen("MaxPenetration.txt","w");
    fclose(stream);
  }

  return;
}


/*----------------------------------------------------------------------*
 | set states from time integrator to prepare collisions   ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandler::SetState(
  Teuchos::RCP<Epetra_Vector> radius,
  Teuchos::RCP<Epetra_Vector> mass)
{
  // node based vectors
  radiusncol_ = LINALG::CreateVector(*discret_->NodeColMap(),false);
  LINALG::Export(*radius,*radiusncol_);
  masscol_ = LINALG::CreateVector(*discret_->NodeColMap(),false);
  LINALG::Export(*mass,*masscol_);

  return;
}


/*----------------------------------------------------------------------*
 | compute collisions (inter-particle and particle-wall)   ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandler::ComputeCollisions(
  double dt,
  Teuchos::RCP<Epetra_Vector> &f_contact,
  Teuchos::RCP<Epetra_Vector> &m_contact
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("TimeforContactSearchAndCalculation");

  if(dt>dt_krit_ and myrank_==0)
  {
    std::cout<<"W A R N I N G : time step larger than critical time step!"<<std::endl;
    std::cout<<"W A R N I N G : calculated critical time step: "<<dt_krit_<<" !"<<std::endl;
  }

  contact_energy_ = 0.0;

  // miraculous transformation from row to col layout ...
  disncol_ = discret_->GetState("bubblepos");
  velncol_ = discret_->GetState("bubblevel");
  ang_velncol_ = discret_->GetState("bubbleangvel");

  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn = walldiscret->GetState("walldisn");
  Teuchos::RCP<const Epetra_Vector> wallveln = walldiscret->GetState("wallveln");

  // define vector for contact force
  Teuchos::RCP<Epetra_FEVector> f_structure = Teuchos::rcp(new Epetra_FEVector(*discret_->DofRowMap()));

  // store bins, which have already been examined
  std::set<int> examinedbins;

  // loop over all particles
  for(int i=0; i<discret_->NodeColMap()->NumMyElements(); ++i)
  {
    DRT::Node *currparticle = discret_->lColNode(i);

    if(currparticle->NumElement() != 1)
      dserror("More than one element for this particle");

    DRT::Element** CurrentBin = currparticle->Elements();
    int binId = CurrentBin[0]->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.count(binId) == 1 )
    {
      continue;
    }
    //else: bin is examined for the first time --> new entry in examinedbins_
    else
    {
      examinedbins.insert(binId);
    }

    int ijk[3];
    particle_algorithm_->ConvertGidToijk(binId,ijk);

    // ijk_range contains: i_min   i_max     j_min     j_max    k_min     k_max
    int ijk_range[] = {ijk[0]-1, ijk[0]+1, ijk[1]-1, ijk[1]+1, ijk[2]-1, ijk[2]+1};
    std::set<int> binIds;

    particle_algorithm_->GidsInijkRange(ijk_range,binIds,true);

    // list of all particles in the neighborhood of currparticle
    std::set<DRT::Node*> neighboringparticles;

    // list of walls that border on the CurrentBin
    std::set<DRT::Element*> neighboringwalls;

    FindNeighbors(neighboringparticles, neighboringwalls, binIds);

    DRT::Node **NodesInCurrentBin = CurrentBin[0]->Nodes();
    int numparticle = CurrentBin[0]->NumNode();

    // loop over all particles in CurrentBin
    for(int i=0; i<numparticle; ++i)
    {
      DRT::Node *particle_i = NodesInCurrentBin[i];

      std::vector<int> lm_i;
      lm_i.reserve(3);

      // extract global dof ids and fill into lm_i
      discret_->Dof(particle_i, lm_i);
      int owner_i = particle_i->Owner();

      //position of particle i
      std::vector<double> myposition_i(3);
      DRT::UTILS::ExtractMyValues(*disncol_,myposition_i,lm_i);

      //velocity of particle i
      std::vector<double> myvel_i(3);
      DRT::UTILS::ExtractMyValues(*velncol_,myvel_i,lm_i);

      //angular-velocity of particle i
      std::vector<double> myangvel_i(3);
      DRT::UTILS::ExtractMyValues(*ang_velncol_,myangvel_i,lm_i);

      int lid = discret_->NodeColMap()->LID(particle_i->Id());

      // radius of particle i
      double radius_i = (*radiusncol_)[lid];

      // mass of particle i
      double mass_i = (*masscol_)[lid];

      // evaluate contact with walls first
      std::map<int,WallContactPoint> surfaces;
      std::map<int,WallContactPoint> lines;
      std::map<int,WallContactPoint> nodes;

      std::set<int> unusedIds;

      // check whether there is contact between particle i and neighboring walls
      for(std::set<DRT::Element*>::const_iterator w=neighboringwalls.begin(); w!=neighboringwalls.end();  ++w)
      {
        int numnodes = (*w)->NumNode();
        std::vector<int> lm_wall;
        lm_wall.reserve(numnodes * 3);

        std::vector<int> lmowner;
        std::vector<int> lmstride;
        (*w)->LocationVector(*walldiscret,lm_wall,lmowner,lmstride);

        // nodal displacements
        std::vector<double> nodal_disp(numnodes * 3);
        DRT::UTILS::ExtractMyValues(*walldisn,nodal_disp,lm_wall);

        // get current position of nodes: x = X + u
        std::map<int,LINALG::Matrix<3,1> > nodeCoord;
        DRT::Node** wallnodes = (*w)->Nodes();
        for(int counter=0; counter<numnodes; ++counter)
        {
          LINALG::Matrix<3,1> currpos;
          const double* X = wallnodes[counter]->X();
          currpos(0) = X[0] + nodal_disp[counter*3+0];
          currpos(1) = X[1] + nodal_disp[counter*3+1];
          currpos(2) = X[2] + nodal_disp[counter*3+2];
          nodeCoord[wallnodes[counter]->Id()] = currpos;
        }

        LINALG::Matrix<3,1> nearestPoint;
        LINALG::Matrix<3,1> position_i;

        //transfer entries from myposition_i to position_i
        for(int n=0; n<3; ++n)
        {
          position_i(n)=myposition_i[n];
        }

        //-------find point on wall element with smallest distance to particle_i-------------------
        GEO::ObjectType objecttype = GEO::nearest3DObjectOnElement((*w),nodeCoord,position_i,nearestPoint);
        //-----------------------------------------------------------------------------------------

        LINALG::Matrix<3,1> r_i_wall;
        r_i_wall.Update(1.0, nearestPoint, -1.0, position_i);
        double distance_i_wall = r_i_wall.Norm2();
        double penetration = distance_i_wall-radius_i;

        if(penetration <= 0.0)
        {
          // get pointer to the current object type of closest point
          std::map<int,WallContactPoint> *pointer=0;
          switch(objecttype)
          {
          case GEO::SURFACE_OBJECT:
          {
            pointer = &surfaces;
          }
          break;
          case GEO::LINE_OBJECT:
          {
            pointer = &lines;
          }
          break;
          case GEO::NODE_OBJECT:
          {
            pointer = &nodes;
          }
          break;
          default:
            dserror("unknown object type");
          break;
          }

          // check, whether point has already been detected (e.g. one line element between two surfaces)
          bool insert = true;
          for(std::map<int,WallContactPoint>::const_iterator iter = (*pointer).begin(); iter != (*pointer).end();++iter)
          {
            LINALG::Matrix<3,1> distance_vector;
            distance_vector.Update(1.0, nearestPoint, -1.0, (iter->second).point);
            double distance = distance_vector.Norm2();
            double adaptedtol = GEO::TOL7 * radius_i;

            if (distance < adaptedtol)
            {
              //point has already been detected --> do not insert
              insert = false;
              unusedIds.insert((*w)->Id());
              break;
            }
          }

          // insert contact point with current surface in corresponding map (surf, line, node)
          if(insert)
          {
            WallContactPoint currentContact = { nearestPoint, penetration , nodeCoord ,lm_wall , lmowner };
            (*pointer).insert(std::pair<int,WallContactPoint>((*w)->Id(),currentContact));
          }
        }
        //penetration > 0.0 --> contact impossible
        else
        {
          unusedIds.insert((*w)->Id());
        }
      }

      // find entries of lines and nodes which are within the penetration volume of the current particle
      // hierarchical: surfaces first
      for(std::map<int,WallContactPoint>::const_iterator surfiter = surfaces.begin(); surfiter != surfaces.end(); ++surfiter)
      {
        // within this radius no other contact point can lie: radius = sqrt(r_i^2 - (r_i-|g|)^2)
        double radius_surface = sqrt(pow(radius_i,2.0) - pow(radius_i-fabs((surfiter->second).penetration),2.0));

        for(std::map<int,WallContactPoint>::const_iterator lineiter = lines.begin(); lineiter != lines.end();++lineiter)
        {
          LINALG::Matrix<3,1> distance_vector;
          distance_vector.Update(1.0, (surfiter->second).point, -1.0, (lineiter->second).point);
          double distance = distance_vector.Norm2();
          if(distance <= radius_surface)
            unusedIds.insert(lineiter->first);
        }
        for(std::map<int,WallContactPoint>::const_iterator nodeiter=nodes.begin(); nodeiter!= nodes.end(); ++nodeiter)
        {
          LINALG::Matrix<3,1> distance_vector;
          distance_vector.Update(1.0, (surfiter->second).point, -1.0, (nodeiter->second).point);
          double distance = distance_vector.Norm2();
          if(distance <= radius_surface)
            unusedIds.insert(nodeiter->first);
        }
      }
      // find entries of nodes which are within the penetration volume of the current particle
      // hierarchical: lines next
      for(std::map<int,WallContactPoint>::const_iterator lineiter=lines.begin(); lineiter != lines.end(); ++lineiter)
      {
        // radius = sqrt(r_i^2 - (r_i-|g|)^2)
        double radius_line = sqrt(pow(radius_i,2.0) - pow(radius_i-fabs((lineiter->second).penetration),2.0));

        for(std::map<int,WallContactPoint>::const_iterator nodeiter=nodes.begin(); nodeiter!=nodes.end(); ++nodeiter)
        {
          LINALG::Matrix<3,1> distance_vector;
          distance_vector.Update(1.0, (lineiter->second).point, -1.0, (nodeiter->second).point);
          double distance = distance_vector.Norm2();
          if(distance <= radius_line)
            unusedIds.insert(nodeiter->first);
        }
      }

      // write entries of lines and nodes to surfaces if contact has to be evaluated
      for(std::map<int,WallContactPoint>::const_iterator iter = lines.begin(); iter != lines.end() ;++iter)
      {
        if( !unusedIds.count(iter->first) )
          surfaces.insert(std::pair<int,WallContactPoint>(iter->first,iter->second));
      }
      for(std::map<int,WallContactPoint>::const_iterator iter=nodes.begin(); iter!=nodes.end(); ++iter)
      {
        if( !unusedIds.count(iter->first) )
          surfaces.insert(std::pair<int,WallContactPoint>(iter->first,iter->second));
      }

      // evaluate contact between particle_i and entries of surfaces
      std::map<int, PARTICLE::Collision>& history_wall = static_cast<PARTICLE::ParticleNode*>(particle_i)->Get_history_wall();
      if(history_wall.size() > 3)
        dserror("Contact with more than 3 wall elements. Check whether history is deleted correctly.");

      for(std::map<int,WallContactPoint>::const_iterator iter = surfaces.begin(); iter != surfaces.end();++iter)
      {
        // gid of wall element
        int gid_wall = iter->first;

        //distance-vector
        double r_contact[3] = {0.0};

        // normal-vector
        double normal[3] = {0.0};

        // distance between centre of mass of two particles
        double norm_r_contact = 0.0;

        // velocity v_rel = v_i - v_j
        double v_rel[3] = {0.0};

        // velocity v_rel_tangential
        double v_rel_tangential[3] = {0.0};

        // part of v_rel in normal-direction
        double v_rel_normal = 0.0;

        // penetration
        double g = 0.0;

        // normalised mass
        double m_eff = mass_i;

        // contact force
        double normalcontactforce = 0.0;
        double tangentcontactforce[3] = {0.0};

        //distance-vector and distance--------------------------
        for(int n=0; n<3; ++n)
        {
          //calculate entries of r_contact
          r_contact[n] = (iter->second).point(n) - myposition_i[n];
          //length of r_contact
          norm_r_contact += r_contact[n]*r_contact[n];
        }
        norm_r_contact = sqrt(norm_r_contact);
        //--------------------------------------------------------

        //normal-vector-------------------------------------------
        for(int n=0; n<3; ++n)
        {
          normal[n] = r_contact[n]/norm_r_contact;
        }
        //-------------------------------------------------------

        //penetration--------------------------------------------
        //g = norm_r_contact - radius_i;
        g = (iter->second).penetration;

        if(fabs(g)>g_max_)
          g_max_ = fabs(g);
        //-------------------------------------------------------

        //-------get velocity of contact point-----------------------
        LINALG::Matrix<3,1> vel_nearestPoint(true);
        LINALG::Matrix<2,1> elecoord(true);
        DRT::Element *CurrentEle = walldiscret->gElement(gid_wall);
        const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(CurrentEle, (iter->second).nodalCoordinates));

        // get coordinates of the projection point in parameter space of the element (xi_coordinates)
        GEO::CurrentToSurfaceElementCoordinates(CurrentEle->Shape(), xyze, (iter->second).point, elecoord);

        int numnodes = CurrentEle->NumNode();
        Epetra_SerialDenseVector funct(numnodes);

        // get shape functions of the element evaluated at the projection point
        DRT::UTILS::shape_function_2D(funct,elecoord(0,0),elecoord(1,0),CurrentEle->Shape());

        std::vector<double> nodal_vel(numnodes * 3);
        DRT::UTILS::ExtractMyValues(*wallveln,nodal_vel,(iter->second).lm);
        for(int node=0; node<numnodes; ++node)
        {
          for(int dim=0; dim<3; ++dim)
          {
            vel_nearestPoint(dim) += funct[node] * nodal_vel[node * 3 + dim];
          }
        }
        //-----------------------------------------------------------

        // velocity v_rel = v_i - v_wall and v_rel in normal-direction: v_rel * n
        for(int n=0; n<3; ++n)
        {
          v_rel[n] = myvel_i[n] - vel_nearestPoint(n);
          v_rel_normal += v_rel[n]*normal[n];
        }

        // normal contact force between particle and wall (note: owner_j = -1)
        CalculateNormalContactForce(g, v_rel_normal, mass_i, normalcontactforce, owner_i, -1);

        if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang)
        {
          // velocity v_rel = v_i + omega_i x (r'_i n) - v_wall = (v_i - v_wall) + omega_i x (r'_i n)
          // and velocity v_rel_tangential
          for(int n=0; n<3; ++n)
          {
            v_rel[n] += (radius_i+g) * (myangvel_i[(n+1)%3] * normal[(n+2)%3] - myangvel_i[(n+2)%3] * normal[(n+1)%3]);
            v_rel_tangential[n] = v_rel[n] - v_rel_normal * normal[n];
          }

          // if g < 0 and g_lasttimestep > 0 -> create history variables
          if(!history_wall.count(gid_wall))
          {
             PARTICLE::Collision col;
             // initialize with stick
             col.stick = true;
             // initialize g_t[3]
             for(int n=0; n<3; ++n)
             {
               col.g_t[n] = 0.0;
             }

             // insert new entry
             history_wall.insert(std::pair<int,PARTICLE::Collision>(gid_wall,col));
          }

          // calculation of tangential contact force
          CalculateTangentialContactForce(normalcontactforce, normal, tangentcontactforce,
                      history_wall[gid_wall], v_rel_tangential, m_eff, dt,owner_i, -1);
        }

        // assembly of contact forces and moments
        Epetra_SerialDenseVector val_i(3);
        Epetra_SerialDenseVector m_i(3);
        std::vector<int> lmowner_i(3);

        double r_i = radius_i + g;
        int owner_i = particle_i->Owner();

        for(int n=0; n<3; ++n)
        {
          lmowner_i[n] = owner_i;
          // contact forces
          val_i[n] = normalcontactforce * normal[n] + tangentcontactforce[n];
          // moments: m_i = (r_i * n) x F_t
          m_i[n] = r_i * (normal[(n+1)%3] * tangentcontactforce[(n+2)%3] - normal[(n+2)%3] * tangentcontactforce[(n+1)%3]);
        }

        // do assembly of contact moments
        LINALG::Assemble(*m_contact,m_i,lm_i,lmowner_i);

        // do assembly of contact forces
        LINALG::Assemble(*f_contact,val_i,lm_i,lmowner_i);

        // forces on wall elements
        double nodal_forces[numnodes * 3];
        for(int node=0; node<numnodes; ++node)
        {
          for(int n=0; n<3; ++n)
          {
            nodal_forces[node * 3 + n] = funct[node] *(- val_i[n]);
          }
        }

        // assembly of contact forces on walls
        if(owner_i == myrank_)
        {
          int err = f_structure->SumIntoGlobalValues(numnodes * 3, &((iter->second).lm)[0], &nodal_forces[0]);
          if (err<0)
            dserror("summing into Epetra_FEVector failed");
        }
      } // end for contact points on surfaces

      if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang)
      {
        //delete those entries in history_wall_ which are no longer in contact with particle_i in current time step
        for(std::set<DRT::Element*>::const_iterator w=neighboringwalls.begin(); w != neighboringwalls.end(); ++w)
        {
          int gid_wall = (*w)->Id();
          if( unusedIds.count(gid_wall) and history_wall.count(gid_wall) )
            history_wall.erase(gid_wall);
        }
      }

      // check whether there is contact between particle i and all other particles in the neighborhood except those which
      // have a lower or equal ID than particle i (--> ignoring self-contact)
      std::map<int, PARTICLE::Collision>& history_particle = static_cast<PARTICLE::ParticleNode*>(particle_i)->Get_history_particle();
      if(history_particle.size() > 12)
        dserror("Contact with more than 12 particles particles. Check whether history is deleted correctly.");

      for(std::set<DRT::Node*>::const_iterator j=neighboringparticles.begin(); j!=neighboringparticles.end(); ++j)
      {
        int gid_j = (*j)->Id();
        //evaluate contact only once!
        if(particle_i->Id() >= gid_j)
          continue;

        int owner_j = (*j)->Owner();
        std::vector<int> lm_j;
        lm_j.reserve(3);

        // extract global dof ids
        discret_->Dof((*j), lm_j);

        // position of particle j
        std::vector<double> myposition_j(3);
        DRT::UTILS::ExtractMyValues(*disncol_,myposition_j,lm_j);

        // velocity of particle j
        std::vector<double> myvel_j(3);
        DRT::UTILS::ExtractMyValues(*velncol_,myvel_j,lm_j);

        // angular velocity of particle j
        std::vector<double> myangvel_j(3);
        DRT::UTILS::ExtractMyValues(*ang_velncol_,myangvel_j,lm_j);

        lid = discret_->NodeColMap()->LID(gid_j);

        // radius of particle j
        double radius_j = (*radiusncol_)[lid];

        // mass of particle j
        double mass_j = (*masscol_)[lid];

        // normalized mass
        double m_eff = mass_i * mass_j / (mass_i + mass_j);

        // contact force
        double normalcontactforce = 0.0;
        double tangentcontactforce[3] = {0.0};

        // distance vector and distance between two particles
        double r_contact[3];
        double norm_r_contact = 0.0;
        for(int n=0; n<3; ++n)
        {
          r_contact[n] = myposition_j[n]-myposition_i[n];
          norm_r_contact += r_contact[n]*r_contact[n];
        }
        norm_r_contact = sqrt(norm_r_contact);

        // penetration
        double g = norm_r_contact - radius_i - radius_j;
        // in case of penetration contact forces and moments are calculated
        if(g <= 0.0)
        {
          if(fabs(g)>g_max_)
            g_max_ = fabs(g);

          // normal vector and velocity v_rel = v_i - v_j and part of v_rel in normal-direction: v_rel * n
          // velocity v_rel
          double v_rel[3];
          // normal vector
          double normal[3];
          // part of v_rel in normal-direction
          double v_rel_normal = 0.0;
          for(int n=0; n<3; ++n)
          {
            normal[n] = r_contact[n]/norm_r_contact;
            v_rel[n] = myvel_i[n] - myvel_j[n];
            v_rel_normal += v_rel[n]*normal[n];
          }

          // calculation of normal contact force
          CalculateNormalContactForce(g, v_rel_normal, m_eff, normalcontactforce, owner_i, owner_j);

          // calculation of tangential contact force
          if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang)
          {
            // velocity v_rel = v_i - v_j + omega_i x (r'_i n) + omega_j x (r'_j n)
            // and velocity v_rel_tangential
            double v_rel_tangential[3];
            for(int n=0; n<3; ++n)
            {
              v_rel[n] += (radius_i+g/2.0) * (myangvel_i[(n+1)%3] * normal[(n+2)%3] - myangvel_i[(n+2)%3] * normal[(n+1)%3]) + (radius_j+g/2.0) * (myangvel_j[(n+1)%3] * normal[(n+2)%3] - myangvel_j[(n+2)%3] * normal[(n+1)%3]);
              v_rel_tangential[n] = v_rel[n] - v_rel_normal * normal[n];
            }

            // if history variables does not exist -> create it
            if(!history_particle.count(gid_j))
            {
              PARTICLE::Collision col;
              // initialize with stick
              col.stick = true;
              //initialize g_t[3]
              for(int n=0; n<3; ++n)
              {
                col.g_t[n] = 0.0;
              }

              //insert new entry
              history_particle.insert(std::pair<int,PARTICLE::Collision>(gid_j,col));
            }

            CalculateTangentialContactForce(normalcontactforce, normal, tangentcontactforce,
                      history_particle[gid_j], v_rel_tangential, m_eff, dt, owner_i, owner_j);
          }

          //----------ASSEMBLY---------------------------------------

          Epetra_SerialDenseVector val_i(3);
          Epetra_SerialDenseVector m_i(3);
          std::vector<int> lmowner_i(3);

          Epetra_SerialDenseVector val_j(3);
          Epetra_SerialDenseVector m_j(3);
          std::vector<int> lmowner_j(3);

          double r_i = radius_i + g/2.0;
          double r_j = radius_j + g/2.0;

          // forces
          for(int n=0; n<3; ++n)
          {
            lmowner_i[n] = owner_i;
            lmowner_j[n] = owner_j;

            val_i[n] = normalcontactforce * normal[n] + tangentcontactforce[n];
            // actio = reactio
            val_j[n] = - val_i[n];
          }

          // moments
          for(int n=0; n<3; ++n)
          {
            // m_i = (r_i * n) x F_t
            m_i[n] = r_i * (normal[(n+1)%3] * tangentcontactforce[(n+2)%3] - normal[(n+2)%3] * tangentcontactforce[(n+1)%3]);

            // m_j = r_j/r_i * m_i
            m_j[n] = r_j/r_i * m_i[n];
          }

          // assembly of contact moments
          LINALG::Assemble(*m_contact,m_i,lm_i,lmowner_i);
          LINALG::Assemble(*m_contact,m_j,lm_j,lmowner_j);

          // assembly of contact forces
          LINALG::Assemble(*f_contact,val_i,lm_i,lmowner_i);
          LINALG::Assemble(*f_contact,val_j,lm_j,lmowner_j);

        }
        else // g > 0.0 --> no contact
        {
          // erase entry in history if still existing
          if(history_particle.count(gid_j))
          {
           history_particle.erase(gid_j);
          }
        }

      }
    }
  }

  // call global assemble for particle forces on walls
  int err = f_structure->GlobalAssemble(Add, false);
  if (err<0)
    dserror("global assemble into fluidforces failed");


  radiusncol_ = Teuchos::null;
  masscol_ =Teuchos::null;
  velncol_ = Teuchos::null;
  disncol_ = Teuchos::null;
  ang_velncol_ = Teuchos::null;

  //apply contact-forces
// structure->SetForceInterface(f_structure);

  return contact_energy_;
}


/*----------------------------------------------------------------------*
 | get neighboring particles and wall elements             ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandler::FindNeighbors(
  std::set<DRT::Node*> &neighboringparticles,
  std::set<DRT::Element*> &neighboringwalls,
  std::set<int> &binIds
  )
{
  // loop over all neighboring bins (including the central bin)
  for(std::set<int>::const_iterator bin=binIds.begin(); bin!=binIds.end(); ++bin)
  {
    DRT::Element *neighboringbin = discret_->gElement(*bin);

    // gather neighboring wall elements
    DRT::Element** WallElements = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(neighboringbin)->AssociatedWallEles();
    int NumWallElements = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(neighboringbin)->NumAssociatedWallEle();
    for(int iwall=0;iwall<NumWallElements;iwall++)
      neighboringwalls.insert(WallElements[iwall]);

    // gather neighboring particles
    DRT::Node **NeighboringNodes = neighboringbin->Nodes();
    int numparticle = neighboringbin->NumNode();
    for(int iNodes=0;iNodes<numparticle;iNodes++)
      neighboringparticles.insert(NeighboringNodes[iNodes]);
  }

  return;
}



/*----------------------------------------------------------------------*
 | calculate normal contact force for single contact pair  ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandler::CalculateNormalContactForce(
  double g,
  double v_rel_normal,
  double m_eff,
  double& normalcontactforce,
  int owner_i,
  int owner_j
  )
{
	// damping parameter
	double d = 0.0;

	//--------------------------------------------------------
	// which contact law: for details see Bachelor thesis Niklas Fehn
	// LinSpring = linear spring
	// Hertz = normal force law of Hertz
	// LinSpringDamp = linear spring damper element
	// LeeHerrmann = nonlinear normal force law of Lee and Herrmann
	// KuwabaraKono = nonlinear normal force law of Kuwabara und Kono
	// Tsuji = nonlinear normal force law of Tsuji
	//---------------------------------------------------------

  // damping parameter
	// TODO: for uni-sized particles this can be done once in the beginning -> efficiency
	// TODO: different m_eff for particle-wall and particel-particle needed when uni-size assumption
	if(owner_j < 0) // contact particle-wall
	{
    if(normal_contact_==INPAR::PARTICLE::LinSpringDamp)
    {
      d = 2 * fabs(log(e_wall_)) * sqrt(k_normal_ * m_eff / (pow(log(e_wall_),2.0)+ pow(M_PI,2.0)));
    }
    else
    {
      d = d_normal_wall_;
    }
	}
	else // contact particle-particle
	{
	  if(normal_contact_==INPAR::PARTICLE::LinSpringDamp)
	  {
	    d = 2 * fabs(log(e_)) * sqrt(k_normal_ * m_eff / (pow(log(e_),2.0)+ pow(M_PI,2.0)));
	  }
	  else
	  {
	    d = d_normal_;
	  }
	}

  // contact force
	switch(normal_contact_)
  {
  case INPAR::PARTICLE::LinSpring:
  {
    normalcontactforce = k_normal_ * g;

#ifdef energyoutput
    //monitor E N E R G Y: here: calculate energy of elastic contact
    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 1.0/2.0 * k_normal_ * g * g;
#endif
  }
  break;
  case INPAR::PARTICLE::Hertz:
  {
    normalcontactforce = - k_normal_ * pow(-g,1.5);

#ifdef energyoutput
    //monitor E N E R G Y: here: calculate energy of elastic contact
    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
#endif
  }
  break;
  case INPAR::PARTICLE::LinSpringDamp:
  {
    normalcontactforce = k_normal_ * g - d * v_rel_normal;

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

#ifdef energyoutput
    //monitor E N E R G Y: here: calculate energy of elastic contact
    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 1.0/2.0 * k_normal_ * g * g;
#endif
  }
  break;
  case INPAR::PARTICLE::LeeHerrmann:
  {
    // m_eff = m_i * m_j / ( m_i + m_j)

    normalcontactforce = - k_normal_ * pow(-g,1.5) - m_eff * d * v_rel_normal;

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

#ifdef energyoutput
    //monitor E N E R G Y: here: calculate energy of elastic contact
    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
#endif
  }
  break;
  case INPAR::PARTICLE::KuwabaraKono:
  {
    normalcontactforce = - k_normal_ * pow(-g,1.5) - d * v_rel_normal * pow(-g,0.5);

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

#ifdef energyoutput
    //monitor E N E R G Y: here: calculate energy of elastic contact
    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
#endif
  }
  break;
  case INPAR::PARTICLE::Tsuji:
  {
    normalcontactforce = - k_normal_ * pow(-g,1.5) - d * v_rel_normal * pow(-g,0.25);

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

#ifdef energyoutput
    //monitor E N E R G Y: here: calculate energy of elastic contact
    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
#endif
  }
  break;
  default:
    dserror("specified normal contact law does not exist");
    break;
  }

	return;
}


/*----------------------------------------------------------------------*
 | calculate tangential contact force for single           ghamm 09/13  |
 | contact pair                                                         |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandler::CalculateTangentialContactForce(
  double normalcontactforce,
  double *normal,
  double *tangentcontactforce,
  PARTICLE::Collision &currentColl,
  double *v_rel_tangential,
  double m_eff,
  const double dt,
  int owner_i,
  int owner_j
  )
{
	// damping parameter
	double d = -0.0;
	// stiffness
	double k = -1.0;
  // frictional coefficient
  double mu = -1.0;

	if(owner_j<0) //contact particle-wall
	{
	  // friction
    mu = mu_wall_;
    // stiffness
    k = k_tang_wall_;
    // damping
    if(d_tang_wall_ < 0.0)
    {
      d = 2 * fabs(log(e_wall_)) * sqrt(k_normal_ * m_eff / (pow(log(e_wall_),2.0)+ pow(M_PI,2.0)));
    }
    else
    {
      d = d_tang_wall_;
    }
	}
	else // contact particle-particle
	{
	  // friction
	  mu = mu_;
	  // stiffness
	  k = k_tang_;
	  // damping
	  if(d_tang_ < 0.0)
	  {
	    d = 2 * fabs(log(e_)) * sqrt(k_normal_ * m_eff / (pow(log(e_),2.0)+ pow(M_PI,2.0)));
	  }
	  else
	  {
	    d = d_tang_;
	  }
	}

	// store length of g_t at time n
	double old_length = 0.0;
	double interime = 0.0;
	for(int n=0; n<3; ++n)
	{
		old_length += currentColl.g_t[n] * currentColl.g_t[n];
		interime += normal[n] * currentColl.g_t[n];
	}
	old_length = sqrt(old_length);

	// projection of g_t onto current normal at time n+1
  double new_length = 0.0;
	for(int n=0; n<3; ++n)
	{
		currentColl.g_t[n] += - interime * normal[n];
		new_length += currentColl.g_t[n] * currentColl.g_t[n];
	}
	new_length = sqrt(new_length);

	// ensure that g_t has the same length as before projection
	// if almost no tangential spring elongation, neglect it
	if(new_length > 1.0E-14)
	{
		for(int n=0; n<3; ++n)
		{
			currentColl.g_t[n] = old_length/new_length * currentColl.g_t[n];
		}
	}

	// update of elastic tangential displacement if stick is true
  if(currentColl.stick == true)
  {
    for(int n=0; n<3; ++n)
    {
      currentColl.g_t[n] += v_rel_tangential[n] * dt;
    }
  }

	// calculate tangential test force
  // norm of tangential contact force
  double norm_f_t = 0.0;
	for(int n=0; n<3; ++n)
	{
		tangentcontactforce[n] = - k * currentColl.g_t[n] - d * v_rel_tangential[n];
		norm_f_t += tangentcontactforce[n] * tangentcontactforce[n];
	}
	norm_f_t = sqrt(norm_f_t);

	// Coulomb friction law

	// tangential contact force for "stick" - case----------------------
	if( norm_f_t <= (mu * fabs(normalcontactforce)) )
	{
	  currentColl.stick = true;
		//tangential contact force already calculated
	}
	else //"slip"-case
	{
	  currentColl.stick = false;
		//calculate tangent vector ( unit vector in (test-)tangentcontactforce-direction )
		double tangent[3];
		for(int n=0; n<3; ++n)
		{
			tangent[n] = tangentcontactforce[n] / norm_f_t;
		}

		// calculate tangent contact force and tangential displacements
		for(int n=0; n<3; ++n)
		{
			tangentcontactforce[n] = mu * fabs(normalcontactforce) * tangent[n];
			currentColl.g_t[n] = - 1/k * (tangentcontactforce[n] + d * v_rel_tangential[n]);
		}
	}
	//---------------------------------------------------------------

#ifdef energyoutput
	new_length = 0.0;
	for(int n=0; n<3; ++n)
	{
		new_length += currentColl.g_t[n] * currentColl.g_t[n];
	}
	new_length = sqrt(new_length);

	contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.5 * k * pow(new_length,2.0);
#endif

	return;
}


/*----------------------------------------------------------------------*
 | read initial contact parameters and validate them       ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandler::ReadContactParameters(double density)
{
  //extract input-parameters
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  contact_strategy_ = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ContactStrategy>(particleparams,"CONTACT_STRATEGY");
  normal_contact_ = DRT::INPUT::IntegralValue<INPAR::PARTICLE::NormalContact>(particleparams,"NORMAL_CONTACT_LAW");

  if(contact_strategy_ != INPAR::PARTICLE::None)
  {
    r_min_ = particleparams.get<double>("MIN_RADIUS");
    r_max_ = particleparams.get<double>("MAX_RADIUS");
    v_max_ = particleparams.get<double>("MAX_VELOCITY");
    c_ = particleparams.get<double>("REL_PENETRATION");
    e_ = particleparams.get<double>("COEFF_RESTITUTION");
    e_wall_ = particleparams.get<double>("COEFF_RESTITUTION_WALL");
    mu_wall_ = particleparams.get<double>("FRICT_COEFF_WALL");
    mu_ = particleparams.get<double>("FRICT_COEFF");
    tension_cutoff_ = (DRT::INPUT::IntegralValue<int>(particleparams,"TENSION_CUTOFF") == 1);

    if(r_min_<0.0 or r_max_<0.0 or v_max_<0.0 or c_<0.0)
      dserror("Invalid input parameter (MIN_RADIUS,MAX_RADIUS,MAX_VELOCITY,REL_PENETRATION have to be larger than zero)");

    if((e_<0.0 or e_wall_<0.0) && normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
      dserror("Invalid input parameter COEFF_RESTITUTION for this kind of contact law");

    //critical time step
    double mass_min = 0.0;
    mass_min = density * 4.0/3.0 * M_PI * pow( r_min_ ,3.0 );

    double k_tkrit = 0.0;

    // here the first element of the structural problem is asked for its material and assumed that wall element has same properties
    Teuchos::RCP<DRT::Discretization> structdis = particle_algorithm_->Structure()->Discretization();
    int local_structmatid = -1;
    // in order to allow procs without wall element communication is necessary
    if(structdis->NumMyColElements() != 0)
    {
      Teuchos::RCP<MAT::Material> mat = structdis->lColElement(0)->Material();
      MAT::PAR::Parameter* params = mat->Parameter();
      local_structmatid = params->Id();
    }
    int structmatid = -1;
    structdis->Comm().MaxAll(&local_structmatid, &structmatid, 1);

    Teuchos::RCP<MAT::StVenantKirchhoff> structmat = Teuchos::rcp_dynamic_cast<MAT::StVenantKirchhoff>(MAT::Material::Factory(structmatid));
    if(structmat == Teuchos::null)
      dserror("only stvenantkirchhoff material is supported so far :-(");

    double G_wall = structmat->ShearMod();
    double nue_wall = structmat->PoissonRatio();

    double G= young_ / (2*(1+nue_));

    //kappa - tangential to normal stiffness ratio
    kappa_ = (1-nue_)/(1-0.5*nue_);
    kappa_wall_ = ( (1-nue_)/G + (1-nue_wall)/G_wall ) / ( (1-0.5*nue_)/G + (1-0.5*nue_wall)/G_wall );

    //------------stiffness----------------------------
    switch(normal_contact_)
    {
    case INPAR::PARTICLE::LinSpring:
    case INPAR::PARTICLE::LinSpringDamp:
    {
      //stiffness calculated from relative penetration and some other input parameters (linear spring)
      k_normal_ = 2.0/3.0 * r_max_ * M_PI * density * pow(v_max_,2.0) / pow(c_,2.0);
      //for tangential contact same stiffness is used
      k_tang_ = kappa_ * k_normal_;
      k_tang_wall_ = kappa_wall_ * k_normal_;
      k_tkrit = k_normal_;

      double user_normal_stiffness = particleparams.get<double>("NORMAL_STIFF");
      if(user_normal_stiffness > 0.0)
      {
        //if user specifies normal stiffness, this stiffness will be used as normal and tangential stiffness for simulation
        k_normal_ = user_normal_stiffness;
        //for tangential contact same stiffness is used
        k_tang_ =kappa_ * k_normal_;
        k_tang_wall_ =kappa_wall_ * k_normal_;
        //stiffness used for calculation of critical time step
        k_tkrit = k_normal_;

        std::cout<<"WARNING: stiffness calculated from relative penetration will be overwritten by input NORMAL_STIFF!!!"<<std::endl;
      }
    }
    break;
    case INPAR::PARTICLE::Hertz:
    case INPAR::PARTICLE::LeeHerrmann:
    case INPAR::PARTICLE::KuwabaraKono:
    case INPAR::PARTICLE::Tsuji:
    {
      if(contact_strategy_==INPAR::PARTICLE::NormalAndTang)
        dserror("tangential contact only with linear normal model implemented");

      //stiffness calculated from relative penetration and some other input parameters (Hertz)
      k_normal_ = 10.0/3.0 * M_PI * density * pow(v_max_,2.0) * pow(r_max_,0.5) / pow(2*c_,2.5);
      //stiffness used for calculation of critical time step (linear spring stiffness needed!)
      k_tkrit = 2.0/3.0 * r_max_ * M_PI * density * pow(v_max_,2.0) / pow(c_,2.0);

      double user_normal_stiffness = particleparams.get<double>("NORMAL_STIFF");
      if(user_normal_stiffness > 0.0)
      {
        //if user specifies normal stiffness, this stiffness will be used as normal stiffness for simulation
        k_normal_ = user_normal_stiffness;
        //for tangential contact the user specified (nonlinear) normal stiffness which has to be transformed into a linear normal
        //stiffness with the same relative penetration which is used as (linear) tangential stiffness afterwards
        double value = 2048.0/1875.0 * density * pow(v_max_,2.0) * M_PI * pow(r_max_,3.0) * pow(k_normal_,4.0);
        //stiffness used for calculation of critical time step (linear spring stiffness needed!)
        k_tkrit = pow(value,0.2);

        std::cout<<"WARNING: stiffness calculated from relative penetration will be overwritten by input NORMAL_STIFF!!!"<<std::endl;
      }
    }
    break;
    default:
      dserror("normal contact law does not exist");
    break;
    }
    //---------------------------------------------------------

    //------------------damping--------------------------------
    d_normal_ = -1.0;
    d_tang_ = -1.0;

    if(normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP");
      if(user_normal_damping >= 0.0)
      {
        dserror("Invalid input parameter NORMAL_DAMP for this kind of contact law");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP");
      if(user_tang_damping >= 0.0)
      {
        dserror("Invalid input parameter TANG_DAMP for this kind of contact law");
      }
    }

    if(normal_contact_ == INPAR::PARTICLE::LeeHerrmann || normal_contact_ == INPAR::PARTICLE::KuwabaraKono ||
        normal_contact_ == INPAR::PARTICLE::Tsuji)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP");
      if(user_normal_damping >= 0.0)
      {
        //user has to specify normal damping coefficient
        d_normal_ = user_normal_damping;
      }
      else
      {
        dserror("For this kind of contact law the input parameter NORMAL_DAMP is invalid");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP");
      if(user_tang_damping >= 0.0)
      {
        //user has to specify tangential damping coefficient
        d_tang_ = user_tang_damping;
      }
      else
      {
        if(contact_strategy_==INPAR::PARTICLE::NormalAndTang)
          dserror("For this kind of contact law the input parameter TANG_DAMP is invalid");
      }
    }
    //------------------damping (wall)--------------------------------
    d_normal_wall_ = -1.0;
    d_tang_wall_ = -1.0;

    if(normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP_WALL");
      if(user_normal_damping >= 0.0)
      {
        dserror("Invalid input parameter NORMAL_DAMP_WALL for this kind of contact law");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP_WALL");
      if(user_tang_damping >= 0.0)
      {
        dserror("Invalid input parameter TANG_DAMP_WALL for this kind of contact law");
      }
    }

    if(normal_contact_ == INPAR::PARTICLE::LeeHerrmann || normal_contact_ == INPAR::PARTICLE::KuwabaraKono ||
        normal_contact_ == INPAR::PARTICLE::Tsuji)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP_WALL");
      if(user_normal_damping >= 0.0)
      {
        //user has to specify normal damping coefficient
        d_normal_wall_ = user_normal_damping;
      }
      else
      {
        dserror("For this kind of contact law the input parameter NORMAL_DAMP_WALL is invalid");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP_WALL");
      if(user_tang_damping >= 0.0)
      {
        //user has to specify tangential damping coefficient
        d_tang_wall_ = user_tang_damping;
      }
      else
      {
        if(contact_strategy_==INPAR::PARTICLE::NormalAndTang)
          dserror("For this kind of contact law the input parameter TANG_DAMP_WALL is invalid");
      }

    }
    //---------------------------------------------------------------

    double factor = 0.0;
    double safety = 0.75;

    //initialize factor
    if(contact_strategy_==INPAR::PARTICLE::Normal)
    { factor = 0.34; }
    if(contact_strategy_==INPAR::PARTICLE::NormalAndTang)
    { factor = 0.22; }

    //calculate critical time step
    dt_krit_ = safety * factor * sqrt( mass_min / k_tkrit );

    //check frictional coefficient
    if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang)
    {
      if(mu_<=0.0 or mu_wall_ <=0.0)
       dserror("Friction coefficient invalid");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | assemble energies of particles                          ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandler::EnergyAssemble(double owner_i,double owner_j)
{
  double value = -1.0;

  const int myrank = discret_->Comm().MyPID();

  //contact  with wall
  if(owner_j<0)
  {
    if(owner_i!=myrank)
     value = 0.0;
  }
  //contact between particle_i and particle_j
  else
  {
    if((owner_i==myrank && owner_j!=myrank) or (owner_i!=myrank && owner_j==myrank))
      value = 0.5;
    else if(owner_i!=myrank and owner_j!=myrank)
      value = 0.0;
  }

  return value;
}


/*----------------------------------------------------------------------*
 | print maximum penetration of particles                  ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandler::PrintMaxPenetration(int step, double time)
{
  if(myrank_==0 && step%100==0)
  {
    FILE *stream = 0;

    stream = fopen("MaxPenetration.txt","a");
    fprintf(stream,"%12.10f\t%12.10f\n",time,g_max_);

    fclose(stream);
  }
  return;
}
