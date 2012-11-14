/*!----------------------------------------------------------------------
\file turbulence_statistics_cha_mult_phase.cpp

\brief calculate mean values and fluctuations for turbulent channel
flows.


*----------------------------------------------------------------------*/

#include "turbulence_statistics_bcf.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dofset.H"

/*----------------------------------------------------------------------

  Standard Constructor (public)

  ---------------------------------------------------------------------*/
COMBUST::TurbulenceStatisticsBcf::TurbulenceStatisticsBcf(
  RefCountPtr<DRT::Discretization>   actdis,
  ParameterList&                     params)
  :
  discret_(actdis),
  params_ (params)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
    dserror("Evaluation of turbulence statistics only for 3d channel flow!");

  //----------------------------------------------------------------------
  // switches, control parameters, material parameters

  // type of fluid flow solver: incompressible, Boussinesq approximation, varying density, loma
  physicaltype_ = DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params, "Physical Type");

  // get the plane normal direction from the parameterlist
  {
    string planestring = params_.sublist("TURBULENCE MODEL").get<string>("HOMDIR","not_specified");

    if(planestring == "xz")
    {
      dim_ = 1;
    }
    else if(planestring == "yz")
    {
      dim_ = 0;
    }
    else if(planestring == "xy")
    {
      dim_ = 2;
    }
    else
    {
      dserror("homogeneuous plane for channel flow was specified incorrectly.");
    }
  }

  // this should not happen, but it does not hurt to check anyways
  if(physicaltype_ == INPAR::FLUID::loma)
  {
    dserror("Loma not supported by the multi phase channel flow");
  }

  //--------------------------------------------------------------
  // get the number and the viscosities of the fluids
  {
    // iterate over all fluids and get fluid viscosity from
    // material definition --- for computation of ltau
    visc_.clear();
    matidtoindex_ = Teuchos::rcp(new std::map<int, int>);
    matidtoindex_->clear();
    int index = 0;
    const int nummat = DRT::Problem::Instance()->Materials()->Num();
    for (int id = 1; id-1 < nummat; ++id)
    {
      Teuchos::RCP<const MAT::PAR::Material> mat = DRT::Problem::Instance()->Materials()->ById(id);
      if (mat == Teuchos::null)
	dserror("Could not find material Id %d", id);
      else
      {
	if (mat->Type() == INPAR::MAT::m_fluid)
	{
	  const MAT::PAR::Parameter* matparam = mat->Parameter();
	  const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(matparam);
	  visc_.push_back(actmat->viscosity_);
	  (*matidtoindex_)[actmat->Id()] = index;
	  index++;
	}
      }
    }
    numphase_ = visc_.size();
  }


  // ---------------------------------------------------------------------
  // up to now, there are no records written
  countrecord_ = 0;

  // ---------------------------------------------------------------------
  // compute all planes for sampling

  // available planes of element nodes (polynomial)
  // this is not only used for the nodal sampling, but also for
  // the integral/volume sampling. If an element has 4 nodes
  // in that plane and the remaining nodes are greater than the
  // plane, it belongs to said plane.
  nodeplanes_ = Teuchos::rcp(new vector<double> );

  // allocate array for bounding box
  //
  //	      |  x  |  y  |  z
  //	------+-----+-----+-----
  //	  min |     |	  |
  //	------+-----+-----+-----
  //	  max |     |	  |
  //
  //
  boundingbox_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(2,3));
  for (int row = 0;row<3;++row)
  {
    (*boundingbox_)(0,row) = +10e+19;
    (*boundingbox_)(1,row) = -10e+19;
  }

  {
    // the criterion allows differences in coordinates by 1e-9
    set<double,PlaneSortCriterion> availablecoords;

    // loop nodes, build set of planes accessible on this proc and
    // calculate bounding box
    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = discret_->lRowNode(i);
      availablecoords.insert(node->X()[dim_]);

      for (int row = 0;row<3;++row)
      {
       if ((*boundingbox_)(0,row)>node->X()[row])
       {
         (*boundingbox_)(0,row)=node->X()[row];
       }
       if ((*boundingbox_)(1,row)<node->X()[row])
       {
         (*boundingbox_)(1,row)=node->X()[row];
       }
      }
    }

    // communicate mins
    for (int row = 0;row<3;++row)
    {
      double min;

      discret_->Comm().MinAll(&((*boundingbox_)(0,row)),&min,1);
      (*boundingbox_)(0,row)=min;
    }

    // communicate maxs
    for (int row = 0;row<3;++row)
    {
      double max;

      discret_->Comm().MaxAll(&((*boundingbox_)(1,row)),&max,1);
      (*boundingbox_)(1,row)=max;
    }

    //--------------------------------------------------------------------
    // round robin loop to communicate coordinates to all procs

    {
#ifdef PARALLEL
      int myrank  =discret_->Comm().MyPID();
#endif
      int numprocs=discret_->Comm().NumProc();

      vector<char> sblock;
      vector<char> rblock;


#ifdef PARALLEL
      // create an exporter for point to point comunication
      DRT::Exporter exporter(discret_->Comm());
#endif

      for (int np=0;np<numprocs;++np)
      {
        DRT::PackBuffer data;

       for (set<double,PlaneSortCriterion>::iterator plane=availablecoords.begin();
            plane!=availablecoords.end();
            ++plane)
       {
         DRT::ParObject::AddtoPack(data,*plane);
       }
       data.StartPacking();
       for (set<double,PlaneSortCriterion>::iterator plane=availablecoords.begin();
            plane!=availablecoords.end();
            ++plane)
       {
         DRT::ParObject::AddtoPack(data,*plane);
       }
       swap( sblock, data() );

#ifdef PARALLEL
       MPI_Request request;
       int   tag    =myrank;

       int   frompid=myrank;
       int   topid  =(myrank+1)%numprocs;

       int   length=sblock.size();

       exporter.ISend(frompid,topid,
        &(sblock[0]),sblock.size(),
        tag,request);

       rblock.clear();

       // receive from predecessor
       frompid=(myrank+numprocs-1)%numprocs;
       exporter.ReceiveAny(frompid,tag,rblock,length);

       if(tag!=(myrank+numprocs-1)%numprocs)
       {
	 dserror("received wrong message (ReceiveAny)");
       }

       exporter.Wait(request);

       {
         // for safety
         exporter.Comm().Barrier();
       }
#else
       // dummy communication
       rblock.clear();
       rblock=sblock;
#endif

       // Unpack received block into set of all planes.
       {
          vector<double> coordsvec;

          coordsvec.clear();

          vector<char>::size_type index = 0;
          while (index < rblock.size())
          {
            double onecoord;
            DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
            availablecoords.insert(onecoord);
          }
       }
      }
    }

    //----------------------------------------------------------------------
    // push coordinates of planes in a vector

    {
      for(set<double,PlaneSortCriterion>::iterator coord=availablecoords.begin();
          coord!=availablecoords.end();
          ++coord)
      {
       nodeplanes_->push_back(*coord);
      }
    }
  }

  //---------------------------------------------------------------------
  // determine the number of nodes per plane
  {
    int loccount = 0;
    for (int i=0; i<discret_->NumMyRowNodes(); ++i)
    {
      DRT::Node* node = discret_->lRowNode(i);
      if ((*nodeplanes_)[0] > node->X()[dim_] - 2e-9 and (*nodeplanes_)[0] < node->X()[dim_] + 2e-9)
         loccount++;
    }
    discret_->Comm().SumAll(&loccount,&numnode_,1);
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of in plane mean values
  // which should be one less than the node planes, but
  // for ease of implementation the last plane is simply
  // left at 0.0

  int size = nodeplanes_->size();

  // arrays for integration based averaging
  // --------------------------------------

  sumvol_.resize(numphase_);

  for (size_t i = 0; i < numphase_; i++)
  {
    sumvol_[i] = Teuchos::rcp(new vector<double>);
    sumvol_[i]->resize(size,0.0);
  }

  // arrays for point based averaging
  // --------------------------------

  pointsumnode_.resize(numphase_);

  // first order moments
  pointsumu_.resize(numphase_);
  pointsumv_.resize(numphase_);
  pointsumw_.resize(numphase_);
  pointsump_.resize(numphase_);

  // second order moments
  pointsumsqu_.resize(numphase_);
  pointsumsqv_.resize(numphase_);
  pointsumsqw_.resize(numphase_);
  pointsumsqp_.resize(numphase_);

  pointsumuv_.resize(numphase_);
  pointsumuw_.resize(numphase_);
  pointsumvw_.resize(numphase_);

  for (size_t i = 0; i < numphase_; i++)
  {
    pointsumnode_[i] =  Teuchos::rcp(new vector<double> );
    pointsumnode_[i]->resize(size,0.0);

    // first order moments
    pointsumu_[i] =  Teuchos::rcp(new vector<double> );
    pointsumu_[i]->resize(size,0.0);

    pointsumv_[i] =  Teuchos::rcp(new vector<double> );
    pointsumv_[i]->resize(size,0.0);

    pointsumw_[i] =  Teuchos::rcp(new vector<double> );
    pointsumw_[i]->resize(size,0.0);

    pointsump_[i] =  Teuchos::rcp(new vector<double> );
    pointsump_[i]->resize(size,0.0);

    // second order moments
    pointsumsqu_[i] =  Teuchos::rcp(new vector<double> );
    pointsumsqu_[i]->resize(size,0.0);

    pointsumsqv_[i] =  Teuchos::rcp(new vector<double> );
    pointsumsqv_[i]->resize(size,0.0);

    pointsumsqw_[i] =  Teuchos::rcp(new vector<double> );
    pointsumsqw_[i]->resize(size,0.0);

    pointsumsqp_[i] =  Teuchos::rcp(new vector<double> );
    pointsumsqp_[i]->resize(size,0.0);

    pointsumuv_[i] =  Teuchos::rcp(new vector<double> );
    pointsumuv_[i]->resize(size,0.0);

    pointsumuw_[i] =  Teuchos::rcp(new vector<double> );
    pointsumuw_[i]->resize(size,0.0);

    pointsumvw_[i] =  Teuchos::rcp(new vector<double> );
    pointsumvw_[i]->resize(size,0.0);
  }


  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file(s)

  Teuchos::RefCountPtr<std::ofstream> log;
//  Teuchos::RefCountPtr<std::ofstream> log_res;

  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");

    {
      s.append(".flow_statistics");

      log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
      (*log) << "# Statistics for turbulent incompressible multi phase channel flow\n\n";

      log->flush();
    }
  }

  // clear statistics
  this->ClearStatistics();

  return;
}// TurbulenceStatisticsChaMultPhase::TurbulenceStatisticsChaMultPhase

/*----------------------------------------------------------------------*
   Destructor
 -----------------------------------------------------------------------*/
COMBUST::TurbulenceStatisticsBcf::~TurbulenceStatisticsBcf()
{
  return;
}// TurbulenceStatisticsBcf::~TurbulenceStatisticsBcf()

/*----------------------------------------------------------------------*
       Compute the in-plane mean values of first and second order
       moments for velocities, pressure and Cs are added to global
       'sum' vectors.
 -----------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsBcf::DoTimeSample(
  Teuchos::RefCountPtr<const Epetra_Vector> stdvelnp,
  Teuchos::RefCountPtr<const Epetra_Vector> stdforce,
  Teuchos::RefCountPtr<const DRT::DofSet>   stddofset,
  Teuchos::RefCountPtr<const Epetra_Vector> discretmatchingvelnp,
  Teuchos::RefCountPtr<const Epetra_Vector> phinp)
{
  // we have an additional sample
  numsamp_++;

  //===============================================
  // create all vectors needed for sampling
  //===============================================
  // we cannot keep them as xfem constantly changes
  // the vector distribution
  stddofset_	= stddofset;
  const Epetra_Map* dofrowmap = stddofset_->DofRowMap();

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  toggleu_.resize(numphase_);
  togglev_.resize(numphase_);
  togglew_.resize(numphase_);
  togglep_.resize(numphase_);

  for (size_t i = 0; i < numphase_; i++)
  {
    toggleu_[i]   = LINALG::CreateVector(*dofrowmap,true);
    togglev_[i]   = LINALG::CreateVector(*dofrowmap,true);
    togglew_[i]   = LINALG::CreateVector(*dofrowmap,true);
    togglep_[i]   = LINALG::CreateVector(*dofrowmap,true);
  }

  //----------------------------------------------------------------------
  // create pointer to velnp as member so we do not have to pass it around.
  // the pointer will be reset after the sample has been taken!!!
  stdvelnp_	= stdvelnp;

  pointsquaredvelnp_ = LINALG::CreateVector(*dofrowmap);
  pointsquaredvelnp_->Multiply(1.0,*stdvelnp_,*stdvelnp_,0.0);

  if (discretmatchingvelnp != Teuchos::null)
  {
    fullvelnp_ = discretmatchingvelnp;
  }
  else
  {
    // if a discretmatching velnp was not provided one could try the
    // stdvelnp, but this will most likely fail
    // fullvelnp_ = stdvelnp;
    dserror("The multi phase turbulence statistics object needs a velnp that matches the one of the discretization.");
  }

  if (phinp != Teuchos::null)
  {
    phinp_ = phinp;
  }
  else
  {
    dserror("Hack! The multi phase turbulence statistics object needs a phinp.");
  }


  //==========================================================
  // start the averaging
  //==========================================================

  //----------------------------------------------------------------------
  // loop planes and calculate integral means in each plane in EvalIntMeanValInPlanes
  this->EvaluateIntegralMeanValuesInPlanes();

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane
  this->EvaluatePointwiseMeanValuesInPlanes();

  //----------------------------------------------------------------------
  // Compute forces on top and bottom plate for normalization purposes.
  // We do not differentiate between the fluids in this part.

  for(vector<double>::iterator plane=nodeplanes_->begin();
      plane!=nodeplanes_->end();
      ++plane)
  {
    // only true for top and bottom plane
    if ((*plane-2e-9 < (*nodeplanes_)[0] and *plane+2e-9 > (*nodeplanes_)[0])
        or
        (*plane-2e-9 < (*nodeplanes_)[nodeplanes_->size()-1] and *plane+2e-9 > (*nodeplanes_)[nodeplanes_->size()-1]))
    {
      // toggle vectors are one in the position of a dof in this plane,
      // else 0
      toggleu_[0]->PutScalar(0.0);
      togglev_[0]->PutScalar(0.0);
      togglew_[0]->PutScalar(0.0);

      // activate toggles for in plane dofs
      for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
      {
         DRT::Node* node = discret_->lRowNode(nn);

         // this node belongs to the plane under consideration
         if (node->X()[dim_]<*plane+2e-9 && node->X()[dim_]>*plane-2e-9)
         {
           vector<int> dof = stddofset_->Dof(node);
           double      one = 1.0;

           toggleu_[0]->ReplaceGlobalValues(1,&one,&(dof[0]));
           togglev_[0]->ReplaceGlobalValues(1,&one,&(dof[1]));
           togglew_[0]->ReplaceGlobalValues(1,&one,&(dof[2]));
         }
      }

      // compute forces by dot product
      double inc=0.0;
      {
         double local_inc=0.0;
         for(int rr=0;rr<(*(toggleu_[0])).MyLength();++rr)
         {
           local_inc+=(*(toggleu_[0]))[rr]*(*(toggleu_[0]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);

         if(abs(inc)<1e-9)
         {
           dserror("there are no forced nodes on the boundary\n");
         }

         local_inc=0.0;
         for(int rr=0; rr<stdforce->MyLength(); ++rr)
         {
           local_inc+=(*stdforce)[rr] * (*(toggleu_[0]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         sumforceu_+=inc;

         local_inc=0.0;
         for(int rr=0;rr<stdforce->MyLength();++rr)
         {
           local_inc+=(*stdforce)[rr] * (*(togglev_[0]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         sumforcev_+=inc;


         local_inc=0.0;
         for(int rr=0;rr<stdforce->MyLength();++rr)
         {
           local_inc+=(*stdforce)[rr] * (*(togglew_[0]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         sumforcew_+=inc;
      }
    }
  }

  //=================================================
  // clean up pointers to velnp, etc.
  //=================================================

  stddofset_	     = Teuchos::null;

  toggleu_.clear();
  togglev_.clear();
  togglew_.clear();

  stdvelnp_     = Teuchos::null;
  pointsquaredvelnp_ = Teuchos::null;
  fullvelnp_     = Teuchos::null;

  phinp_     = Teuchos::null;

  return;
}// TurbulenceStatisticsBcf::DoTimeSample


/*----------------------------------------------------------------------*
  Compute in plane means of u,u^2 etc. (integral version)
 -----------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsBcf::EvaluateIntegralMeanValuesInPlanes()
{

  //----------------------------------------------------------------------
  // loop elements and perform integration over homogeneous plane

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  eleparams.set("action","calc_turbulence_statistics");

  // choose what to assemble
  eleparams.set("assemble matrix 1",false);
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",false);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);

  // set parameter list
  eleparams.set("normal direction to homogeneous plane",dim_);
  eleparams.set("coordinate vector for hom. planes",nodeplanes_);

  // set size of vectors
  const int size = sumvol_[0]->size();

  // generate processor local result vectors
  vector<RefCountPtr<vector<double> > > locvol( numphase_);
  vector<RefCountPtr<vector<double> > > globvol(numphase_);

  for (size_t i = 0; i < numphase_; i++)
  {
    locvol[i] = Teuchos::rcp(new vector<double> );
    locvol[i]->resize(size,0.0);

    globvol[i] = Teuchos::rcp(new vector<double>);
    globvol[i]->resize(size,0.0);
  }

  // communicate pointers to the result vectors to the element
  eleparams.set("map materialid to index", matidtoindex_);

  eleparams.set("element volume", &locvol);

  // counts the number of elements in the lowest homogeneous plane
  // (the number is the same for all planes, since we use a structured
  //  cartesian mesh)
  int locprocessedeles=0;

  eleparams.set("count processed elements", &locprocessedeles);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1,converged)", fullvelnp_);

  // call loop over elements
  discret_->Evaluate(eleparams, null, null, null, null, null);
  discret_->ClearState();


  //----------------------------------------------------------------------
  // add contributions from all processors
  // it looks a bit messy due to the vector<RCP<vector<double> > > construct
  for (size_t i = 0; i < numphase_; i++)
  {
    discret_->Comm().SumAll(&((*(locvol[i]))[0]), &((*(globvol[i]))[0]), size);
  }

  //----------------------------------------------------------------------
  // the sums are divided by the layers area to get the area average

  {
    discret_->Comm().SumAll(&locprocessedeles, &numele_, 1);

    for(size_t i = 0; i < sumvol_[0]->size()-1; ++i) // the last plane is empty by definition
    {
      double fullvolplane = 0.0;
      for(size_t j = 0; j < numphase_; j++)
         fullvolplane += (*(globvol[j]))[i];

      if (fullvolplane == 0.0)
         dserror("Apparently we got an empty (vol=0.0) element plane");

      for(size_t j = 0; j < numphase_; j++)
      {
         (*(sumvol_[j]))[i] += (*(globvol[j]))[i] / fullvolplane;
      }
    }
  }

  return;
}// TurbulenceStatisticsBcf::EvaluateIntegralMeanValuesInPlanes()


/*----------------------------------------------------------------------*
  Compute in plane means of u,u^2 etc. (nodal quantities)
  ----------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsBcf::EvaluatePointwiseMeanValuesInPlanes()
{
  int planenum = 0;

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane

  for(vector<double>::iterator plane=nodeplanes_->begin();
      plane!=nodeplanes_->end();
      ++plane)
  {
    // holds the nodes per phase in this plane on this proc
    vector<int> countnodesinplane(numphase_, 0);

    // toggle vectors are one in the position of a dof in this plane,
    // else 0
    for (size_t i = 0; i < numphase_; i++)
    {
      toggleu_[i]->PutScalar(0.0);
      togglev_[i]->PutScalar(0.0);
      togglew_[i]->PutScalar(0.0);
      togglep_[i]->PutScalar(0.0);
    }

    //----------------------------------------------------------------------
    // activate toggles for in plane dofs

    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn) //perhaps we should keep the list nodes per plane...
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the plane under consideration
      if (node->X()[dim_]<*plane+2e-9 && node->X()[dim_]>*plane-2e-9)
      {
         vector<int> dof = stddofset_->Dof(node);
         double      one = 1.0;
         int  index = -1;

         // This is a dirty hack.  The statistics channel should  not know
         // anything about  the combust specific phi.  There ought to be a
         // method to let either the node decide to which phase it belongs
         // or to ask the discretization or similar.
         // Perhaps a node condition?
         const int lid = phinp_->Map().LID(node->Id());
         if (lid < 0)
           dserror("Node %i is not on this proc.", node->Id());
         if ((*phinp_)[lid] < 0 )
           index = 1;
         else
           index = 0;
         // end hack

         toggleu_[index]->ReplaceGlobalValues(1,&one,&(dof[0]));
         togglev_[index]->ReplaceGlobalValues(1,&one,&(dof[1]));
         togglew_[index]->ReplaceGlobalValues(1,&one,&(dof[2]));
         togglep_[index]->ReplaceGlobalValues(1,&one,&(dof[3]));

         // now check whether we have a pbc condition on this node
         vector<DRT::Condition*> mypbc;

         node->GetCondition("SurfacePeriodic",mypbc);

         // yes, we have a pbc
         if (mypbc.size()>0)
         {
           // loop them and check, whether this is a pbc pure master node
           // for all previous conditions
           size_t ntimesmaster = 0;
           for (size_t numcond=0;numcond<mypbc.size();++numcond)
           {
             const string* mymasterslavetoggle
               = mypbc[numcond]->Get<string>("Is slave periodic boundary condition");

             if(*mymasterslavetoggle=="Master")
             {
               ++ntimesmaster;
             } // end is slave?
           } // end loop this conditions

           if(ntimesmaster!=mypbc.size())
           {
             continue;
           }
           // we have a master. Remember this cause we have to extend the patch
           // to the other side...
         }

         (countnodesinplane[index])++;
      }
    }

    vector<int> countnodesinplaneonallprocs(numphase_, 0);

    discret_->Comm().SumAll(&(countnodesinplane[0]), &(countnodesinplaneonallprocs[0]), numphase_);

    for (size_t iphase = 0; iphase < numphase_; ++iphase)
    {
      // if there are no nodes in this plane, there is nothing to do anyways
      if (countnodesinplaneonallprocs[iphase])
      {
         double inc=0.0;
         double local_inc=0.0;
         //----------------------------------------------------------------------
         // sum up the number of nodes per plane
         (*(pointsumnode_[iphase]))[planenum] += countnodesinplaneonallprocs[iphase];

         //----------------------------------------------------------------------
         // compute scalar products from velnp and toggle vec to sum up
         // values in this plane
         local_inc=0.0;
         for(int rr=0; rr<stdvelnp_->MyLength(); ++rr)
         {
           local_inc += (*stdvelnp_)[rr] * (*(toggleu_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumu_[iphase]))[planenum] += inc;

         local_inc=0.0;
         for(int rr=0;rr<stdvelnp_->MyLength();++rr)
         {
           local_inc+=(*stdvelnp_)[rr] * (*(togglev_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumv_[iphase]))[planenum] += inc;

         local_inc=0.0;
         for(int rr=0;rr<stdvelnp_->MyLength();++rr)
         {
           local_inc+=(*stdvelnp_)[rr]*(*(togglew_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumw_[iphase]))[planenum] += inc;

         local_inc=0.0;
         for(int rr=0;rr<stdvelnp_->MyLength();++rr)
         {
           local_inc+=(*stdvelnp_)[rr]*(*(togglep_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsump_[iphase]))[planenum]+=inc;

         //----------------------------------------------------------------------
         // compute scalar products from squaredvelnp and toggle vec to
         // sum up values for second order moments in this plane

         local_inc=0.0;
         for(int rr=0; rr < stdvelnp_->MyLength(); rr += 4)
         {
           local_inc += ((*stdvelnp_)[rr]*(*(toggleu_[iphase]))[rr]) * ((*stdvelnp_)[rr+1]*(*(togglev_[iphase]))[rr+1]);
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumuv_[iphase]))[planenum]+=inc;

         local_inc=0.0;
         for(int rr=0; rr < stdvelnp_->MyLength(); rr += 4)
         {
           local_inc += ((*stdvelnp_)[rr]*(*(toggleu_[iphase]))[rr]) * ((*stdvelnp_)[rr+2]*(*(togglew_[iphase]))[rr+2]);
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumuw_[iphase]))[planenum]+=inc;

         local_inc=0.0;
         for(int rr=0; rr < stdvelnp_->MyLength(); rr += 4)
         {
           local_inc += ((*stdvelnp_)[rr+1]*(*(togglev_[iphase]))[rr+1]) * ((*stdvelnp_)[rr+2]*(*(togglew_[iphase]))[rr+2]);
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumvw_[iphase]))[planenum]+=inc;

         //----------------------------------------------------------------------
         // compute scalar products from squaredvelnp and toggle vec to
         // sum up values for second order moments in this plane
         local_inc=0.0;
         for(int rr=0;rr<pointsquaredvelnp_->MyLength();++rr)
         {
           local_inc+=(*pointsquaredvelnp_)[rr]*(*(toggleu_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumsqu_[iphase]))[planenum] += inc;

         local_inc=0.0;
         for(int rr=0;rr<pointsquaredvelnp_->MyLength();++rr)
         {
           local_inc+=(*pointsquaredvelnp_)[rr]*(*(togglev_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumsqv_[iphase]))[planenum] += inc;

         local_inc=0.0;
         for(int rr=0;rr<pointsquaredvelnp_->MyLength();++rr)
         {
           local_inc+=(*pointsquaredvelnp_)[rr]*(*(togglew_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumsqw_[iphase]))[planenum] += inc;

         local_inc=0.0;
         for(int rr=0;rr<pointsquaredvelnp_->MyLength();++rr)
         {
           local_inc+=(*pointsquaredvelnp_)[rr]*(*(togglep_[iphase]))[rr];
         }
         discret_->Comm().SumAll(&local_inc,&inc,1);
         (*(pointsumsqp_[iphase]))[planenum] += inc;
      }
    } //end for (size_t iphase = 0; iphase < numphase_; ++iphase)
    planenum++;
  }
  return;
}// TurbulenceStatisticsBcf::EvaluatePointwiseMeanValuesInPlanes()



/*----------------------------------------------------------------------*

       Compute a time average of the mean values over all steps
	  since the last output. Dump the result to file.

  ----------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsBcf::TimeAverageMeansAndOutputOfStatistics(int step)
{
  if (numsamp_ == 0)
  {
    dserror("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average

  for(size_t i=0; i<sumvol_[0]->size(); ++i)
  {
    for (size_t j = 0; j < numphase_; j++)
    {
      (*(sumvol_[j]))[i] /= numsamp_;

//(*(sumu_[j])  )[i] /=numsamp_;
//(*(sumv_[j])  )[i] /=numsamp_;
//(*(sumw_[j])  )[i] /=numsamp_;
//(*(sump_[j])  )[i] /=numsamp_;
//
//(*(sumuv_[j]) )[i] /=numsamp_;
//(*(sumuw_[j]) )[i] /=numsamp_;
//(*(sumvw_[j]) )[i] /=numsamp_;
//
//(*(sumsqu_[j]))[i] /=numsamp_;
//(*(sumsqv_[j]))[i] /=numsamp_;
//(*(sumsqw_[j]))[i] /=numsamp_;
//(*(sumsqp_[j]))[i] /=numsamp_;
    }
  }

  for(size_t i=0; i<pointsumu_[0]->size(); ++i)
  {
    for (size_t j = 0; j < numphase_; j++)
    {
      const double numsamp = ((*(pointsumnode_[j]))[i]);
      if(numsamp > 0.0)
      {
         (*(pointsumu_[j]))[i]	/= numsamp;
         (*(pointsumv_[j]))[i]	/= numsamp;
         (*(pointsumw_[j]))[i]	/= numsamp;
         (*(pointsump_[j]))[i]	/= numsamp;

         (*(pointsumsqu_[j]))[i] /= numsamp;
         (*(pointsumsqv_[j]))[i] /= numsamp;
         (*(pointsumsqw_[j]))[i] /= numsamp;
         (*(pointsumsqp_[j]))[i] /= numsamp;

         (*(pointsumuv_[j]))[i]  /= numsamp;
         (*(pointsumuw_[j]))[i]  /= numsamp;
         (*(pointsumvw_[j]))[i]  /= numsamp;
      }
    }
  }

  sumforceu_/=numsamp_;
  sumforcev_/=numsamp_;
  sumforcew_/=numsamp_;

  //----------------------------------------------------------------------
  // evaluate area to calculate u_tau, l_tau (and tau_W)
  double area = 1.0;
  for (int i=0;i<3;i++)
  {
    if(i!=dim_)
    {
      area*=((*boundingbox_)(1,i)-(*boundingbox_)(0,i));
    }
  }
  // there are two Dirichlet boundaries
  area*=2;

  //----------------------------------------------------------------------
  // we expect nonzero forces (tractions) only in flow direction

  // ltau is used to compute y+
  vector<double> ltau(numphase_);
  if	  (sumforceu_>sumforcev_ && sumforceu_>sumforcew_)
  {
    if(abs(sumforceu_)< 1.0e-12)
    {
      dserror("zero force during computation of wall shear stress\n");
    }

    for (size_t i = 0; i < numphase_; i++)
      ltau[i] = visc_[i]/sqrt(sumforceu_/area);
  }
  else if (sumforcev_>sumforceu_ && sumforcev_>sumforcew_)
  {
    for (size_t i = 0; i < numphase_; i++)
      ltau[i] = visc_[i]/sqrt(sumforcev_/area);
  }
  else if (sumforcew_>sumforceu_ && sumforcew_>sumforcev_)
  {
    for (size_t i = 0; i < numphase_; i++)
      ltau[i] = visc_[i]/sqrt(sumforcew_/area);
  }
  else
  {
    dserror("Cannot determine flow direction by traction (seems to be not unique)");
  }

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RefCountPtr<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::app));
    (*log) << "\n\n\n";
    (*log) << "# Statistics record " << countrecord_;
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    for (size_t j = 0; j < numphase_; ++j)
    {
      // get matid by index
      int matid = -1;
      for (std::map<int, int>::const_iterator it = matidtoindex_->begin(); it != matidtoindex_->end(); ++it)
      {
         if ((unsigned)it->second == j)
         {
           matid = it->first;
           break;
         }
      }

      (*log) << "# Material ID  	 : ";
      (*log) << "   " << std::setw(11) << std::setprecision(1) << matid;
      (*log) << &endl;

      (*log) << "# (u_tau)^2 = tau_W/rho : ";
      (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforceu_/area;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforcev_/area;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforcew_/area;
      (*log) << &endl;


      (*log) << "#|-------------------------------------|";
      (*log) << "-------------------------------------------------point";
      (*log) << "wise-------------------------------------------";
      (*log) << "------------------------------------------------------------------------|\n";

      (*log) << "#     y        y+      volfraction";
      (*log) << "   umean         vmean      wmean     pmean";
      (*log) << "   mean u^2      mean v^2      mean w^2";
      (*log) << "      mean u*v      mean u*w     mean v*w    mean p^2    no. samp \n";
      (*log) << std::scientific;
      for(size_t i=0; i<nodeplanes_->size(); ++i)
      {
         // y, y+ and volfraction
         (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*nodeplanes_)[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*nodeplanes_)[i]/ltau[j];
         (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*(sumvol_[j]))[i];

         // pointwise means
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumu_   [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumv_   [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumw_   [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsump_   [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqu_ [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqv_ [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqw_ [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumuv_  [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumuw_  [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumvw_  [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqp_ [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumnode_[j]))[i];
         (*log) << "   \n";
      }
    }
    log->flush();
  } // end myrank 0

  // log was written, so increase counter for records
  countrecord_++;

  return;

}// TurbulenceStatisticsBcf::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*

      Compute a time average of the mean values over all steps
       of the sampling period so far. Dump the result to file.

  ----------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsBcf::DumpStatistics(int step)
{
  if (numsamp_ == 0)
  {
    dserror("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average

  //----------------------------------------------------------------------
  // evaluate area to calculate u_tau, l_tau (and tau_W)
  double area = 1.0;
  for (int i=0;i<3;i++)
  {
    if(i!=dim_)
    {
      area*=((*boundingbox_)(1,i)-(*boundingbox_)(0,i));
    }
  }
  // there are two Dirichlet boundaries
  area*=2;

  //----------------------------------------------------------------------
  // we expect nonzero forces (tractions) only in flow direction

  // ltau is used to compute y+
  vector<double> ltau(numphase_);
  if	  (sumforceu_>sumforcev_ && sumforceu_>sumforcew_)
  {
    for (size_t i = 0; i < numphase_; i++)
      ltau[i] = visc_[i]/sqrt(sumforceu_/(area*numsamp_));
  }
  else if (sumforcev_>sumforceu_ && sumforcev_>sumforcew_)
  {
    for (size_t i = 0; i < numphase_; i++)
      ltau[i] = visc_[i]/sqrt(sumforcev_/(area*numsamp_));
  }
  else if (sumforcew_>sumforceu_ && sumforcew_>sumforcev_)
  {
    for (size_t i = 0; i < numphase_; i++)
      ltau[i] = visc_[i]/sqrt(sumforcew_/(area*numsamp_));
  }
  else
  {
    dserror("Cannot determine flow direction by traction (seems to be not unique)");
  }

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RefCountPtr<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
    (*log) << "# Statistics for turbulent incompressible channel flow (first- and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    for (size_t j = 0; j < numphase_; ++j)
    {
      // get matid by index
      int matid = -1;
      for (std::map<int, int>::const_iterator it = matidtoindex_->begin(); it != matidtoindex_->end(); ++it)
      {
         if ((unsigned)it->second == j)
         {
           matid = it->first;
           break;
         }
      }

      (*log) << "# Material ID  	 : ";
      (*log) << "   " << std::setw(11) << std::setprecision(1) << matid;
      (*log) << &endl;

      (*log) << "# (u_tau)^2 = tau_W/rho : ";
      (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforceu_/(area*numsamp_);
      (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforcev_/(area*numsamp_);
      (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforcew_/(area*numsamp_);
      (*log) << &endl;

      (*log) << "#     y	    y+   volfraction";
      (*log) << "	    umean	  vmean 	wmean	      pmean";
      (*log) << "	 mean u^2      mean v^2      mean w^2	  mean p^2";
      (*log) << "      mean u*v      mean u*w	   mean v*w    no. samp";
      (*log) << "\n";

      (*log) << std::scientific;
      for(size_t i=0; i<nodeplanes_->size(); ++i)
      {
         const double numsamp = ((*(pointsumnode_[j]))[i]);
         if(numsamp > 0.0)
         {
           (*(pointsumu_[j]))[i]   /= numsamp;
           (*(pointsumv_[j]))[i]   /= numsamp;
           (*(pointsumw_[j]))[i]   /= numsamp;
           (*(pointsump_[j]))[i]   /= numsamp;

           (*(pointsumsqu_[j]))[i] /= numsamp;
           (*(pointsumsqv_[j]))[i] /= numsamp;
           (*(pointsumsqw_[j]))[i] /= numsamp;
           (*(pointsumsqp_[j]))[i] /= numsamp;

           (*(pointsumuv_[j]))[i]   /= numsamp;
           (*(pointsumuw_[j]))[i]   /= numsamp;
           (*(pointsumvw_[j]))[i]   /= numsamp;
         }

         (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*nodeplanes_)[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*nodeplanes_)[i]/ltau[j];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(sumvol_[j]))[i]/numsamp_;

         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumu_       [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumv_       [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumw_       [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsump_       [j]))[i];

         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqu_     [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqv_     [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqw_     [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumsqp_     [j]))[i];

         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumuv_       [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumuw_       [j]))[i];
         (*log) << "   " << std::setw(11) << std::setprecision(4) << (*(pointsumvw_       [j]))[i];

         (*log) << "   " << std::setw(11) << std::setprecision(4) << numsamp;

         (*log) << "\n\n";
      }
    }
    log->flush();
  }

  return;

}// TurbulenceStatisticsBcf::DumpStatistics


/*----------------------------------------------------------------------*

		  Reset sums and number of samples to 0

  ----------------------------------------------------------------------*/
void COMBUST::TurbulenceStatisticsBcf::ClearStatistics()
{
  // reset the number of samples
  numsamp_ =0;

  // reset forces (mean values and values at bottom and top wall)
  sumforceu_=0;
  sumforcev_=0;
  sumforcew_=0;

  // reset integral and pointwise averages
  for (size_t j=0; j<numphase_; ++j)
  {
    for(size_t i=0; i<pointsumu_[0]->size(); ++i)
    {
      (*(sumvol_[j]))[i]  =0;

      (*(pointsumnode_[j]))[i]=0;

      (*(pointsumu_[j]))[i]  =0;
      (*(pointsumv_[j]))[i]  =0;
      (*(pointsumw_[j]))[i]  =0;
      (*(pointsump_[j]))[i]  =0;

      (*(pointsumsqu_[j]))[i]=0;
      (*(pointsumsqv_[j]))[i]=0;
      (*(pointsumsqw_[j]))[i]=0;
      (*(pointsumsqp_[j]))[i]=0;

      (*(pointsumuv_[j]))[i] =0;
      (*(pointsumuw_[j]))[i] =0;
      (*(pointsumvw_[j]))[i] =0;
    }
  }

  return;
}// TurbulenceStatisticsBcf::ClearStatistics

