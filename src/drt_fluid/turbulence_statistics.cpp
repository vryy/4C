/*!----------------------------------------------------------------------
\file turbulence_statistics.cpp

\brief calculate mean values and fluctuations for turbulent channel
flows.


*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "turbulence_statistics.H"

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

  <pre>
  o Create vector of homogeneous plane coordinates

  o Allocate 4 distributed toggle vectors and one distributed vector
  for squares
  </pre>

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatistics::TurbulenceStatistics(
  RefCountPtr<DRT::Discretization> actdis,
  ParameterList&                   params)
  :
  discret_(actdis),
  params_ (params)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
  {
    dserror("Evaluation of turbulence statistics only for 3d channel flow!");
  }

  // up to now, there are no records written
  countrecord_ = 0;

  // get the plane normal direction from the parameterlist
  {
    string planestring = params_.sublist("TURBULENCE MODEL").get<string>("CHANNEL_HOMPLANE","xz");

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

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  meanvelnp_    = LINALG::CreateVector(*dofrowmap,true);

  toggleu_      = LINALG::CreateVector(*dofrowmap,true);
  togglev_      = LINALG::CreateVector(*dofrowmap,true);
  togglew_      = LINALG::CreateVector(*dofrowmap,true);
  togglep_      = LINALG::CreateVector(*dofrowmap,true);



  // available planes of element nodes (polynomial)/corners
  // (Nurbs) of elements
  nodeplanes_ = rcp(new vector<double> );


  // available homogeneous (sampling) planes --- there are
  // numsubdivisions layers per element layer between two
  // nodes (Polynomial)/per element layer (Nurbs)
  planecoordinates_ = rcp(new vector<double> );

  const int numsubdivisions=5;


  // try to cast discretisation to nurbs variant
  // this tells you what kind of computation of
  // samples is required
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*actdis));

  // allocate array for bounding box
  //
  //          |  x  |  y  |  z
  //    ------+-----+-----+-----
  //      min |     |     |
  //    ------+-----+-----+-----
  //      max |     |     |
  //
  //
  boundingbox_ = rcp(new Epetra_SerialDenseMatrix(2,3));
  for (int row = 0;row<3;++row)
  {
    (*boundingbox_)(0,row) = +10e+19;
    (*boundingbox_)(1,row) = -10e+19;
  }


  // get fluid viscosity from material definition --- for computation
  // of ltau
  visc_ = mat->m.fluid->viscosity;

  if(nurbsdis == NULL)
  {
    //----------------------------------------------------------------------
    // create set of available homogeneous planes. The normal direction
    // is read from the parameter list
    //----------------------------------------------------------------------
    planecoordinates_ = rcp(new vector<double> );

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
    //--------------------------------------------------------------------
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
	// export set to sendbuffer
	sblock.clear();

	for (set<double,PlaneSortCriterion>::iterator plane=availablecoords.begin();
	     plane!=availablecoords.end();
	     ++plane)
	{
	  DRT::ParObject::AddtoPack(sblock,*plane);
	}
#ifdef PARALLEL
	MPI_Request request;
	int         tag    =myrank;

	int         frompid=myrank;
	int         topid  =(myrank+1)%numprocs;

	int         length=sblock.size();

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

	//--------------------------------------------------
	// Unpack received block into set of all planes.
	{
	  vector<double> coordsvec;

	  coordsvec.clear();

	  int index = 0;
	  while (index < (int)rblock.size())
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
    //----------------------------------------------------------------------
    {
      nodeplanes_ = rcp(new vector<double> );


      for(set<double,PlaneSortCriterion>::iterator coord=availablecoords.begin();
	  coord!=availablecoords.end();
	  ++coord)
      {
	nodeplanes_->push_back(*coord);
      }

      //----------------------------------------------------------------------
      // insert additional sampling planes (to show influence of quadratic
      // shape functions)
      //----------------------------------------------------------------------

      for(unsigned rr =0; rr < nodeplanes_->size()-1; ++rr)
      {
	double delta = ((*nodeplanes_)[rr+1]-(*nodeplanes_)[rr])/((double) numsubdivisions);

	for (int mm =0; mm < numsubdivisions; ++mm)
	{
	  planecoordinates_->push_back((*nodeplanes_)[rr]+delta*mm);
	}
      }
      planecoordinates_->push_back((*nodeplanes_)[(*nodeplanes_).size()-1]);
    }
  }
  else
  {

    // pointwise sampling does not make any sense for Nurbs
    // discretisations since shape functions are not interpolating


    // planecoordinates are determined by the element (cartesian) number
    // in y direction and the number of sampling planes in between
    // and nodeplanes are kept as the corners of elements
    // (to be able to visualise stuff on the element center later on)

    // for nurbs discretisations, all vector sizes are already determined
    // by the knotvector size
    if(dim_!=1)
    {
      dserror("For the nurbs stuff, we require that xz is the hom. plane\n");
    }

    // get nurbs dis' knotvector sizes
    vector<int> n_x_m_x_l(nurbsdis->Return_n_x_m_x_l());

    // get nurbs dis' element numbers
    vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele());

    // get the knotvector itself
    RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

    // resize and initialise to 0
    {
      (*nodeplanes_      ).resize(nele_x_mele_x_lele[1]+1);
      (*planecoordinates_).resize(nele_x_mele_x_lele[1]*(numsubdivisions-1)+1);

      vector<double>::iterator coord;

      for (coord  = (*nodeplanes_).begin();
	   coord != (*nodeplanes_).end()  ;
	   ++coord)
      {
	*coord=0;
      }
      for (coord  = planecoordinates_->begin();
	   coord != planecoordinates_->end()  ;
	   ++coord)
      {
	*coord=0;
      }
    }

    // get element map
    const Epetra_Map* elementmap = nurbsdis->ElementRowMap();

    // loop all available elements
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));
      DRT::Node**   nodes = actele->Nodes();

      // get gid, location in the patch
      int gid = actele->Id();

      vector<int> ele_cart_id=knots->ConvertEleGidToKnotIds(gid);

      // want to loop all control points of the element,
      // so get the number of points
      const int numnp = actele->NumNode();

	// access elements knot span
      std::vector<blitz::Array<double,1> > knots(3);
      (*((*nurbsdis).GetKnotVector())).GetEleKnots(knots,actele->Id());

      // aquire weights from nodes
      blitz::Array<double,1> weights(numnp);

      for (int inode=0; inode<numnp; ++inode)
      {
	DRT::NURBS::ControlPoint* cp
	  =
	  dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

	weights(inode) = cp->W();
      }

      // get shapefunctions, compute all visualisation point positions
      blitz::Array<double,1> nurbs_shape_funct(numnp);

      switch (actele->Shape())
      {
      case DRT::Element::nurbs8:
      case DRT::Element::nurbs27:
      {
	// element local point position
	blitz::Array<double,1> uv(3);

	{
	  // standard

	  //               v
	  //              /
          //  w  7       /   8
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  // 5|/        6/  |
	  //  +---------+   |
	  //  |         |   |
	  //  |         |   +
	  //  |         |  / 4
	  //  |         | /
	  //  |         |/
	  //  +---------+ ----->u
	  // 1           2
	  // use v-coordinate of point 1 and 8
	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }

	  (*nodeplanes_      )[ele_cart_id[1]]                    +=x[1];
	  (*planecoordinates_)[ele_cart_id[1]*(numsubdivisions-1)]+=x[1];

	  for (int isd=0; isd<3; ++isd)
	  {
	    if ((*boundingbox_)(0,isd)>x[isd])
	    {
	      (*boundingbox_)(0,isd)=x[isd];
	    }
	    if ((*boundingbox_)(1,isd)<x[isd])
	    {
	      (*boundingbox_)(1,isd)=x[isd];
	    }
	  }

	  for(int rr=1;rr<numsubdivisions-1;++rr)
	  {
	    uv(1) += 2.0/(numsubdivisions-1);

	    DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						  uv               ,
						  knots            ,
						  weights          ,
						  actele->Shape()  );
	    for (int isd=0; isd<3; ++isd)
	    {
	      double val = 0;
	      for (int inode=0; inode<numnp; ++inode)
	      {
		val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	      }
	      x[isd]=val;
	    }
	    (*planecoordinates_)[ele_cart_id[1]*(numsubdivisions-1)+rr]+=x[1];
	  }


	  // set upper point of element, too (only for last layer)
	  if(ele_cart_id[1]+1 == nele_x_mele_x_lele[1])
	  {
	    // point 8
	    uv(0)=  1.0;
	    uv(1)=  1.0;
	    uv(2)=  1.0;
	    DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						  uv               ,
						  knots            ,
						  weights          ,
						  actele->Shape()  );
	    for (int isd=0; isd<3; ++isd)
	    {
	      double val = 0;
	      for (int inode=0; inode<numnp; ++inode)
	      {
		val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	      }
	      x[isd]=val;
	    }

	    (*nodeplanes_)      [ele_cart_id[1]                   +1]+=x[1];
	    (*planecoordinates_)[(ele_cart_id[1]+1)*(numsubdivisions-1)]+=x[1];

            for (int isd=0; isd<3; ++isd)
            {
              if ((*boundingbox_)(0,isd)>x[isd])
              {
                (*boundingbox_)(0,isd)=x[isd];
              }
              if ((*boundingbox_)(1,isd)<x[isd])
              {
                (*boundingbox_)(1,isd)=x[isd];
              }
            }
	  }

	}
	break;
      }
      default:
	dserror("Unknown element shape for a nurbs element or nurbs type not valid for turbulence calculation\n");
      }
    }


    //----------------------------------------------------------------------
    // add contributions from all processors, normalize
    //----------------------------------------------------------------------
    vector<double> lnodeplanes      (*nodeplanes_      );
    vector<double> lplanecoordinates(*planecoordinates_);

    discret_->Comm().SumAll(&(lnodeplanes[0]      ),&((*nodeplanes_      )[0]),nodeplanes_->size()      );
    discret_->Comm().SumAll(&(lplanecoordinates[0]),&((*planecoordinates_)[0]),planecoordinates_->size());

    {
      (*nodeplanes_      ).resize(nele_x_mele_x_lele[1]+1);
      (*planecoordinates_).resize(nele_x_mele_x_lele[1]*(numsubdivisions-1)+1);

      vector<double>::iterator coord;

      int nelelayer = (nele_x_mele_x_lele[0])*(nele_x_mele_x_lele[2]);

      for (coord  = (*nodeplanes_).begin();
	   coord != (*nodeplanes_).end()  ;
	   ++coord)
      {
	*coord/=(double)(nelelayer);
      }
      for (coord  = planecoordinates_->begin();
	   coord != planecoordinates_->end()  ;
	   ++coord)
      {
	*coord/=(double)(nelelayer);
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

  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of in plane mean values
  //----------------------------------------------------------------------
  int size = planecoordinates_->size();

  //----------------------------------------------------------------------
  // arrays for integration based averaging

  // first order moments
  sumu_ =  rcp(new vector<double> );
  sumu_->resize(size,0.0);

  sumv_ =  rcp(new vector<double> );
  sumv_->resize(size,0.0);

  sumw_ =  rcp(new vector<double> );
  sumw_->resize(size,0.0);

  sump_ =  rcp(new vector<double> );
  sump_->resize(size,0.0);

  // now the second order moments
  sumsqu_ =  rcp(new vector<double> );
  sumsqu_->resize(size,0.0);

  sumsqv_ =  rcp(new vector<double> );
  sumsqv_->resize(size,0.0);

  sumsqw_ =  rcp(new vector<double> );
  sumsqw_->resize(size,0.0);

  sumuv_ =  rcp(new vector<double> );
  sumuv_->resize(size,0.0);

  sumuw_ =  rcp(new vector<double> );
  sumuw_->resize(size,0.0);

  sumvw_ =  rcp(new vector<double> );
  sumvw_->resize(size,0.0);

  sumsqp_ =  rcp(new vector<double> );
  sumsqp_->resize(size,0.0);



  //----------------------------------------------------------------------
  // arrays for point based averaging

  pointsquaredvelnp_  = LINALG::CreateVector(*dofrowmap,true);

  // first order moments
  pointsumu_ =  rcp(new vector<double> );
  pointsumu_->resize(size,0.0);

  pointsumv_ =  rcp(new vector<double> );
  pointsumv_->resize(size,0.0);

  pointsumw_ =  rcp(new vector<double> );
  pointsumw_->resize(size,0.0);

  pointsump_ =  rcp(new vector<double> );
  pointsump_->resize(size,0.0);

  // now the second order moments
  pointsumsqu_ =  rcp(new vector<double> );
  pointsumsqu_->resize(size,0.0);

  pointsumsqv_ =  rcp(new vector<double> );
  pointsumsqv_->resize(size,0.0);

  pointsumsqw_ =  rcp(new vector<double> );
  pointsumsqw_->resize(size,0.0);

  pointsumsqp_ =  rcp(new vector<double> );
  pointsumsqp_->resize(size,0.0);

  //----------------------------------------------------------------------
  // arrays for averaging of Smagorinsky constant etc.
  //
  // means for the Smagorinsky constant
  sumCs_  =  rcp(new vector<double> );
  sumCs_->resize(nodeplanes_->size()-1,0.0);

  incrsumCs_  =  rcp(new vector<double> );
  incrsumCs_->resize(nodeplanes_->size()-1,0.0);

  // means for (Cs*delta)^2
  sumCs_delta_sq_  =  rcp(new vector<double> );
  sumCs_delta_sq_->resize(nodeplanes_->size()-1,0.0);

  incrsumCs_delta_sq_  =  rcp(new vector<double> );
  incrsumCs_delta_sq_->resize(nodeplanes_->size()-1,0.0);

  // means for the effective viscosity
  sumvisceff_  =  rcp(new vector<double> );
  sumvisceff_->resize(nodeplanes_->size()-1,0.0);

  incrsumvisceff_  =  rcp(new vector<double> );
  incrsumvisceff_->resize(nodeplanes_->size()-1,0.0);

  //----------------------------------------------------------------------
  // arrays for averaging of residual, subscales etc.

  // means for comparison of of residual and subscale acceleration
  sumres_    =  rcp(new vector<double> );
  sumres_->resize(3*(nodeplanes_->size()-1),0.0);
  sumres_sq_ =  rcp(new vector<double> );
  sumres_sq_->resize(3*(nodeplanes_->size()-1),0.0);

  sumsacc_   =  rcp(new vector<double> );
  sumsacc_->resize(3*(nodeplanes_->size()-1),0.0);
  sumsacc_sq_=  rcp(new vector<double> );
  sumsacc_sq_->resize(3*(nodeplanes_->size()-1),0.0);

  sumsvelaf_=  rcp(new vector<double> );
  sumsvelaf_->resize(3*(nodeplanes_->size()-1),0.0);
  sumsvelaf_sq_=  rcp(new vector<double> );
  sumsvelaf_sq_->resize(3*(nodeplanes_->size()-1),0.0);

  sumresC_       =  rcp(new vector<double> );
  sumresC_->resize(nodeplanes_->size()-1,0.0);
  sumresC_sq_    =  rcp(new vector<double> );
  sumresC_sq_->resize(nodeplanes_->size()-1,0.0);

  sumspresacc_   =  rcp(new vector<double> );
  sumspresacc_->resize(nodeplanes_->size()-1,0.0);
  sumspresacc_sq_=  rcp(new vector<double> );
  sumspresacc_sq_->resize(nodeplanes_->size()-1,0.0);

  sumspressnp_   =  rcp(new vector<double> );
  sumspressnp_->resize(nodeplanes_->size()-1,0.0);
  sumspressnp_sq_=  rcp(new vector<double> );
  sumspressnp_sq_->resize(nodeplanes_->size()-1,0.0);

  //----------------------------------------------------------------------
  // initialise output
  //----------------------------------------------------------------------
  Teuchos::RefCountPtr<std::ofstream> log;
  Teuchos::RefCountPtr<std::ofstream> log_Cs;
  Teuchos::RefCountPtr<std::ofstream> log_res;

  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistic");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# Flow statistics for turbulent channel flow (first- and second-order moments)\n\n";

    log->flush();

    // additional output for dynamic Smagorinsky model
    if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
        ==
        "CLASSICAL_LES")
    {
      if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Dynamic_Smagorinsky"
         ||
         params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Smagorinsky_with_van_Driest_damping"
        )
      {
        std::string s_smag = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
        s_smag.append(".Cs_statistic");

        log_Cs = Teuchos::rcp(new std::ofstream(s_smag.c_str(),ios::out));
        (*log_Cs) << "# Statistics for turbulent channel flow (Smagorinsky constant)\n\n";
      }
    }

    // output of residuals and subscale quantities
    std::string s_res = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s_res.append(".res_statistic");

    log_res = Teuchos::rcp(new std::ofstream(s_res.c_str(),ios::out));
    (*log_res) << "# Statistics for turbulent channel flow (residuals and subscale quantities)\n";
    (*log_res) << "# All values are first averaged over the integration points in an element \n";
    (*log_res) << "# and after that averaged over a whole element layer in the homogeneous plane\n\n";

  }

  // clear statistics
  this->ClearStatistics();

  return;
}// TurbulenceStatistics::TurbulenceStatistics

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatistics::~TurbulenceStatistics()
{
  return;
}// TurbulenceStatistics::~TurbulenceStatistics()

void FLD::TurbulenceStatistics::DoTimeSample(
  Teuchos::RefCountPtr<Epetra_Vector> velnp,
  Epetra_Vector & force
  )
{
  //----------------------------------------------------------------------
  // we have an additional sample
  //----------------------------------------------------------------------
  numsamp_++;

  //----------------------------------------------------------------------
  // meanvelnp is a refcount copy of velnp
  //----------------------------------------------------------------------
  meanvelnp_->Update(1.0,*velnp,0.0);

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // loop planes and calculate integral means in each plane
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  this->EvaluateIntegralMeanValuesInPlanes();

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // try to cast discretisation to nurbs variant
  // this tells you whether pointwise computation of
  // samples is allowed
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));

  if(nurbsdis == NULL)
  {
    this->EvaluatePointwiseMeanValuesInPlanes();
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // compute forces on top and bottom plate for normalization purposes
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  for(vector<double>::iterator plane=planecoordinates_->begin();
      plane!=planecoordinates_->end();
      ++plane)
  {
    //----------------------------------------------------------------------
    // only true for top and bottom plane
    //----------------------------------------------------------------------
    if ((*plane-2e-9 < (*planecoordinates_)[0]
         &&
         *plane+2e-9 > (*planecoordinates_)[0])
        ||
        (*plane-2e-9 < (*planecoordinates_)[planecoordinates_->size()-1]
         &&
         *plane+2e-9 > (*planecoordinates_)[planecoordinates_->size()-1])
      )
    {
      // toggle vectors are one in the position of a dof in this plane,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);

      //----------------------------------------------------------------------
      // activate toggles for in plane dofs
      //----------------------------------------------------------------------
      for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // this node belongs to the plane under consideration
        if (node->X()[dim_]<*plane+2e-9 && node->X()[dim_]>*plane-2e-9)
        {
          vector<int> dof = discret_->Dof(node);
          double      one = 1.0;

          toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
          togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
          togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        }
      }

      double inc;
      force.Dot(*toggleu_,&inc);
      sumforceu_+=inc;
      force.Dot(*togglev_,&inc);
      sumforcev_+=inc;
      force.Dot(*togglew_,&inc);
      sumforcew_+=inc;
    }
  }

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // add increment of last iteration to the sum of Cs values
  // (statistics for dynamic Smagorinsky model)
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
       ||
       params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Smagorinsky_with_van_Driest_damping"
      )
    {
      for (unsigned rr=0;rr<(*incrsumCs_).size();++rr)
      {
        (*sumCs_         )[rr]+=(*incrsumCs_         )[rr];
        (*sumCs_delta_sq_)[rr]+=(*incrsumCs_delta_sq_)[rr];
        (*sumvisceff_    )[rr]+=(*incrsumvisceff_    )[rr];
      }
    }
  }

  return;
}// TurbulenceStatistics::DoTimeSample

/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatistics::EvaluateIntegralMeanValuesInPlanes()
{

  //----------------------------------------------------------------------
  // loop elements and perform integration over homogeneous plane
  //----------------------------------------------------------------------
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
  eleparams.set("coordinate vector for hom. planes",planecoordinates_);

  // set size of vectors
  int size = sumu_->size();

  // generate processor local result vectors
  RefCountPtr<vector<double> > locsumu =  rcp(new vector<double> );
  locsumu->resize(size,0.0);

  RefCountPtr<vector<double> > locsumv =  rcp(new vector<double> );
  locsumv->resize(size,0.0);

  RefCountPtr<vector<double> > locsumw =  rcp(new vector<double> );
  locsumw->resize(size,0.0);

  RefCountPtr<vector<double> > locsump =  rcp(new vector<double> );
  locsump->resize(size,0.0);

  RefCountPtr<vector<double> > locsumsqu =  rcp(new vector<double> );
  locsumsqu->resize(size,0.0);

  RefCountPtr<vector<double> > locsumsqv =  rcp(new vector<double> );
  locsumsqv->resize(size,0.0);

  RefCountPtr<vector<double> > locsumsqw =  rcp(new vector<double> );
  locsumsqw->resize(size,0.0);

  RefCountPtr<vector<double> > locsumuv  =  rcp(new vector<double> );
  locsumuv->resize(size,0.0);

  RefCountPtr<vector<double> > locsumuw  =  rcp(new vector<double> );
  locsumuw->resize(size,0.0);

  RefCountPtr<vector<double> > locsumvw  =  rcp(new vector<double> );
  locsumvw->resize(size,0.0);

  RefCountPtr<vector<double> > locsumsqp =  rcp(new vector<double> );
  locsumsqp->resize(size,0.0);

  RefCountPtr<vector<double> > globsumu =  rcp(new vector<double> );
  globsumu->resize(size,0.0);

  RefCountPtr<vector<double> > globsumv =  rcp(new vector<double> );
  globsumv->resize(size,0.0);

  RefCountPtr<vector<double> > globsumw =  rcp(new vector<double> );
  globsumw->resize(size,0.0);

  RefCountPtr<vector<double> > globsump =  rcp(new vector<double> );
  globsump->resize(size,0.0);

  RefCountPtr<vector<double> > globsumsqu =  rcp(new vector<double> );
  globsumsqu->resize(size,0.0);

  RefCountPtr<vector<double> > globsumsqv =  rcp(new vector<double> );
  globsumsqv->resize(size,0.0);

  RefCountPtr<vector<double> > globsumsqw =  rcp(new vector<double> );
  globsumsqw->resize(size,0.0);

  RefCountPtr<vector<double> > globsumuv  =  rcp(new vector<double> );
  globsumuv->resize(size,0.0);

  RefCountPtr<vector<double> > globsumuw  =  rcp(new vector<double> );
  globsumuw->resize(size,0.0);

  RefCountPtr<vector<double> > globsumvw  =  rcp(new vector<double> );
  globsumvw->resize(size,0.0);

  RefCountPtr<vector<double> > globsumsqp =  rcp(new vector<double> );
  globsumsqp->resize(size,0.0);

  // communicate pointers to the result vectors to the element
  eleparams.set("mean velocity u"     ,locsumu);
  eleparams.set("mean velocity v"     ,locsumv);
  eleparams.set("mean velocity w"     ,locsumw);
  eleparams.set("mean pressure p"     ,locsump);

  eleparams.set("mean value u^2",locsumsqu);
  eleparams.set("mean value v^2",locsumsqv);
  eleparams.set("mean value w^2",locsumsqw);
  eleparams.set("mean value uv" ,locsumuv );
  eleparams.set("mean value uw" ,locsumuw );
  eleparams.set("mean value vw" ,locsumvw );
  eleparams.set("mean value p^2",locsumsqp);

  // counts the number of elements in the lowest homogeneous plane
  // (the number is the same for all planes, since we use a structured
  //  cartesian mesh)
  int locprocessedeles=0;

  eleparams.set("count processed elements",&locprocessedeles);


  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p (n+1,converged)"    ,meanvelnp_);

  // call loop over elements
  discret_->Evaluate(eleparams,null,null,null,null,null);
  discret_->ClearState();


  //----------------------------------------------------------------------
  // add contributions from all processors
  //----------------------------------------------------------------------

  discret_->Comm().SumAll(&((*locsumu)[0]),&((*globsumu)[0]),size);
  discret_->Comm().SumAll(&((*locsumv)[0]),&((*globsumv)[0]),size);
  discret_->Comm().SumAll(&((*locsumw)[0]),&((*globsumw)[0]),size);
  discret_->Comm().SumAll(&((*locsump)[0]),&((*globsump)[0]),size);

  discret_->Comm().SumAll(&((*locsumsqu)[0]),&((*globsumsqu)[0]),size);
  discret_->Comm().SumAll(&((*locsumsqv)[0]),&((*globsumsqv)[0]),size);
  discret_->Comm().SumAll(&((*locsumsqw)[0]),&((*globsumsqw)[0]),size);
  discret_->Comm().SumAll(&((*locsumuv)[0]) ,&((*globsumuv)[0]) ,size);
  discret_->Comm().SumAll(&((*locsumuw)[0]) ,&((*globsumuw)[0]) ,size);
  discret_->Comm().SumAll(&((*locsumvw)[0]) ,&((*globsumvw)[0]) ,size);
  discret_->Comm().SumAll(&((*locsumsqp)[0]),&((*globsumsqp)[0]),size);


  //----------------------------------------------------------------------
  // the sums are divided by the number of elements to get the area
  // average
  //----------------------------------------------------------------------
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));

  if(nurbsdis == NULL)
  {
    discret_->Comm().SumAll(&locprocessedeles,&numele_,1);
  }
  else
  {
    // get nurbs dis' element numbers
    vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele());

    numele_ = nele_x_mele_x_lele[0]*nele_x_mele_x_lele[2];
  }


  for(unsigned i=0; i<planecoordinates_->size(); ++i)
  {
    (*sumu_)[i]  +=(*globsumu)[i];
    (*sumv_)[i]  +=(*globsumv)[i];
    (*sumw_)[i]  +=(*globsumw)[i];
    (*sump_)[i]  +=(*globsump)[i];

    (*sumsqu_)[i]+=(*globsumsqu)[i];
    (*sumsqv_)[i]+=(*globsumsqv)[i];
    (*sumsqw_)[i]+=(*globsumsqw)[i];
    (*sumuv_)[i] +=(*globsumuv) [i];
    (*sumuw_)[i] +=(*globsumuw) [i];
    (*sumvw_)[i] +=(*globsumvw) [i];
    (*sumsqp_)[i]+=(*globsumsqp)[i];
  }

  return;

}// TurbulenceStatistics::EvaluateIntegralMeanValuesInPlanes()


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatistics::EvaluatePointwiseMeanValuesInPlanes()
{
  int planenum = 0;

  //----------------------------------------------------------------------
  // pointwise multiplication to get squared values
  //----------------------------------------------------------------------
  pointsquaredvelnp_->Multiply(1.0,*meanvelnp_,*meanvelnp_,0.0);


  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane
  //----------------------------------------------------------------------
  for(vector<double>::iterator plane=planecoordinates_->begin();
      plane!=planecoordinates_->end();
      ++plane)
  {

    // toggle vectors are one in the position of a dof in this plane,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes in plane (required to calc. in plane mean)
    int countnodesinplane=0;

    //----------------------------------------------------------------------
    // activate toggles for in plane dofs
    //----------------------------------------------------------------------
    for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    {
      DRT::Node* node = discret_->lRowNode(nn);

      // this node belongs to the plane under consideration
      if (node->X()[dim_]<*plane+2e-9 && node->X()[dim_]>*plane-2e-9)
      {
        vector<int> dof = discret_->Dof(node);
        double      one = 1.0;

        toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
        togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
        togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
        togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));

        // now check whether we have a pbc condition on this node
        vector<DRT::Condition*> mypbc;

        node->GetCondition("SurfacePeriodic",mypbc);

        // yes, we have a pbc
        if (mypbc.size()>0)
        {
          // loop them and check, whether this is a pbc pure master node
          // for all previous conditions
          unsigned ntimesmaster = 0;
          for (unsigned numcond=0;numcond<mypbc.size();++numcond)
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
        countnodesinplane++;
      }
    }

    int countnodesinplaneonallprocs=0;

    discret_->Comm().SumAll(&countnodesinplane,&countnodesinplaneonallprocs,1);

    if (countnodesinplaneonallprocs)
    {
      //----------------------------------------------------------------------
      // compute scalar products from velnp and toggle vec to sum up
      // values in this plane
      //----------------------------------------------------------------------
      double inc;
      meanvelnp_->Dot(*toggleu_,&inc);
      (*pointsumu_)[planenum]+=inc/countnodesinplaneonallprocs;
      meanvelnp_->Dot(*togglev_,&inc);
      (*pointsumv_)[planenum]+=inc/countnodesinplaneonallprocs;
      meanvelnp_->Dot(*togglew_,&inc);
      (*pointsumw_)[planenum]+=inc/countnodesinplaneonallprocs;
      meanvelnp_->Dot(*togglep_,&inc);
      (*pointsump_)[planenum]+=inc/countnodesinplaneonallprocs;

      //----------------------------------------------------------------------
      // compute scalar products from squaredvelnp and toggle vec to
      // sum up values for second order moments in this plane
      //----------------------------------------------------------------------
      pointsquaredvelnp_->Dot(*toggleu_,&inc);
      (*pointsumsqu_)[planenum]+=inc/countnodesinplaneonallprocs;
      pointsquaredvelnp_->Dot(*togglev_,&inc);
      (*pointsumsqv_)[planenum]+=inc/countnodesinplaneonallprocs;
      pointsquaredvelnp_->Dot(*togglew_,&inc);
      (*pointsumsqw_)[planenum]+=inc/countnodesinplaneonallprocs;
      pointsquaredvelnp_->Dot(*togglep_,&inc);
      (*pointsumsqp_)[planenum]+=inc/countnodesinplaneonallprocs;
    }
    planenum++;
  }

  return;
}// TurbulenceStatistics::EvaluatePointwiseMeanValuesInPlanes()


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatistics::TimeAverageMeansAndOutputOfStatistics(int step)
{
  if (numsamp_ == 0)
  {
    dserror("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average
  int aux = numele_*numsamp_;

  for(unsigned i=0; i<planecoordinates_->size(); ++i)
  {

    (*sumu_)[i]   /=aux;
    (*sumv_)[i]   /=aux;
    (*sumw_)[i]   /=aux;
    (*sump_)[i]   /=aux;

    (*sumuv_)[i]  /=aux;
    (*sumuw_)[i]  /=aux;
    (*sumvw_)[i]  /=aux;

    (*sumsqu_)[i] /=aux;
    (*sumsqv_)[i] /=aux;
    (*sumsqw_)[i] /=aux;
    (*sumsqp_)[i] /=aux;
  }

  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));

  if(nurbsdis == NULL)
  {
    for(unsigned i=0; i<planecoordinates_->size(); ++i)
    {
      // the pointwise values have already been normalised by
      // "countnodesinplaneonallprocs", so we just divide by
      // the number of time samples
      (*pointsumu_)[i]   /=numsamp_;
      (*pointsumv_)[i]   /=numsamp_;
      (*pointsumw_)[i]   /=numsamp_;
      (*pointsump_)[i]   /=numsamp_;

      (*pointsumsqu_)[i] /=numsamp_;
      (*pointsumsqv_)[i] /=numsamp_;
      (*pointsumsqw_)[i] /=numsamp_;
      (*pointsumsqp_)[i] /=numsamp_;
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
  int flowdirection =0;

  // ltau is used to compute y+
  double ltau = 0;
  if      (sumforceu_>sumforcev_ && sumforceu_>sumforcew_)
  {
    flowdirection=0;
    ltau = visc_/sqrt(sumforceu_/area);
  }
  else if (sumforcev_>sumforceu_ && sumforcev_>sumforcew_)
  {
    flowdirection=1;
    ltau = visc_/sqrt(sumforcev_/area);
  }
  else if (sumforcew_>sumforceu_ && sumforcew_>sumforcev_)
  {
    flowdirection=2;
    ltau = visc_/sqrt(sumforcew_/area);
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
    s.append(".flow_statistic");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
    (*log) << "\n\n\n";
    (*log) << "# Statistics record " << countrecord_;
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    (*log) << "# (u_tau)^2 = tau_W/rho : ";
    (*log) << "   " << setw(11) << setprecision(4) << sumforceu_/area;
    (*log) << "   " << setw(11) << setprecision(4) << sumforcev_/area;
    (*log) << "   " << setw(11) << setprecision(4) << sumforcew_/area;
    (*log) << &endl;


    (*log) << "#|-------------------";
    (*log) << "----------------------------------------------------------";
    (*log) << "--integration based-------------------------";
    (*log) << "----------------------------------------------------------|";
    (*log) << "-------------------------------------------------point";
    (*log) << "wise---------------------------------------";
    (*log) << "------------|\n";

    (*log) << "#     y            y+";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "        mean u^2      mean v^2      mean w^2";
    (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2";
    (*log) << "       umean         vmean         wmean         pmean";
    (*log) << "        mean u^2      mean v^2      mean w^2";
    (*log) << "       mean p^2 \n";
    (*log) << scientific;
    for(unsigned i=0; i<planecoordinates_->size(); ++i)
    {
      (*log) <<  " "  << setw(11) << setprecision(4) << (*planecoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << (*planecoordinates_)[i]/ltau;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumu_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << (*sumv_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << (*sumw_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << (*sump_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << ((*sumsqu_)[i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*sumsqv_)[i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*sumsqw_)[i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*sumuv_) [i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*sumuw_) [i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*sumvw_) [i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*sumsqp_)[i]);
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumu_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumv_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumw_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsump_)   [i];
      (*log) << "   " << setw(11) << setprecision(4) << ((*pointsumsqu_)[i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*pointsumsqv_)[i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*pointsumsqw_)[i]);
      (*log) << "   " << setw(11) << setprecision(4) << ((*pointsumsqp_)[i]);
      (*log) << "   \n";
    }
    log->flush();
  }

  if (discret_->Comm().MyPID()==0)
  {

    // additional output for dynamic Smagorinsky model
    if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
        ==
        "CLASSICAL_LES")
    {
      if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Dynamic_Smagorinsky"
         ||
         params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Smagorinsky_with_van_Driest_damping"
        )
      {
        // get the outfile
        Teuchos::RefCountPtr<std::ofstream> log_Cs;

        std::string s_smag = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
        s_smag.append(".Cs_statistic");

        log_Cs = Teuchos::rcp(new std::ofstream(s_smag.c_str(),ios::app));

        (*log_Cs) << "\n\n\n";
        (*log_Cs) << "# Statistics record " << countrecord_;
        (*log_Cs) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";


        (*log_Cs) << "#     y      ";
        (*log_Cs) << "     Cs     ";
        (*log_Cs) << "   (Cs*hk)^2 ";
        (*log_Cs) << "    visceff  ";
        (*log_Cs) << &endl;
        (*log_Cs) << scientific;
        for (unsigned rr=0;rr<sumCs_->size();++rr)
        {
          (*log_Cs) << setw(11) << setprecision(4) << 0.5*((*nodeplanes_)[rr+1]+(*nodeplanes_)[rr]) << "  " ;
          (*log_Cs) << setw(11) << setprecision(4) << ((*sumCs_)[rr])/(numele_*numsamp_) << "  ";
          (*log_Cs) << setw(11) << setprecision(4) << ((*sumCs_delta_sq_)[rr])/(numele_*numsamp_)<< "  " ;
          (*log_Cs) << setw(11) << setprecision(4) << ((*sumvisceff_)[rr])/(numele_*numsamp_) << &endl;
        }
        log_Cs->flush();
      }
    }
  }

  if (discret_->Comm().MyPID()==0)
  {
    Teuchos::RefCountPtr<std::ofstream> log_res;

    // output of residuals and subscale quantities
    std::string s_res = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s_res.append(".res_statistic");

    log_res = Teuchos::rcp(new std::ofstream(s_res.c_str(),ios::app));

    (*log_res) << "\n\n\n";
    (*log_res) << "# Statistics record " << countrecord_;
    (*log_res) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";
    (*log_res) << "#       y    ";
    (*log_res) << "    res_x  ";
    (*log_res) << "      res_y  ";
    (*log_res) << "      res_z  ";
    (*log_res) << "     sacc_x  ";
    (*log_res) << "     sacc_y  ";
    (*log_res) << "     sacc_z  ";
    (*log_res) << "     svel_x  ";
    (*log_res) << "     svel_y  ";
    (*log_res) << "     svel_z  ";
    (*log_res) << "      resC   ";
    (*log_res) << "     spacc   ";
    (*log_res) << "    spresnp   ";

    (*log_res) << "   res_sq_x  ";
    (*log_res) << "   res_sq_y  ";
    (*log_res) << "   res_sq_z  ";
    (*log_res) << "   sacc_sq_x ";
    (*log_res) << "   sacc_sq_y ";
    (*log_res) << "   sacc_sq_z ";
    (*log_res) << "   svel_sq_x ";
    (*log_res) << "   svel_sq_y ";
    (*log_res) << "   svel_sq_z ";
    (*log_res) << "    resC_sq  ";
    (*log_res) << "   spacc_sq  ";
    (*log_res) << "  spvelnp_sq "  <<&endl;

    (*log_res) << scientific;
    for (unsigned rr=0;rr<nodeplanes_->size()-1;++rr)
    {
      (*log_res)  << setw(11) << setprecision(4) << 0.5*((*nodeplanes_)[rr+1]+(*nodeplanes_)[rr]) << "  " ;

      (*log_res)  << setw(11) << setprecision(4) << (*sumres_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumsvelaf_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsvelaf_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsvelaf_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumresC_)[rr]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumspresacc_)[rr]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumspressnp_)[rr]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumres_sq_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_sq_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_sq_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_sq_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_sq_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_sq_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumsvelaf_sq_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsvelaf_sq_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsvelaf_sq_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumresC_sq_)[rr]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumspresacc_sq_)[rr]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumspressnp_sq_)[rr]/(numele_*numsamp_) << "  ";


      (*log_res)  << &endl;
    }
    log_res->flush();
  }


  // log was written, so increase counter for records
  countrecord_++;

  return;

}// TurbulenceStatistics::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatistics::DumpStatistics(int step)
{
  if (numsamp_ == 0)
  {
    dserror("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average
  int aux = numele_*numsamp_;

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
  int flowdirection =0;

  // ltau is used to compute y+
  double ltau = 0;
  if      (sumforceu_>sumforcev_ && sumforceu_>sumforcew_)
  {
    flowdirection=0;
    ltau = visc_/sqrt(sumforceu_/(area*numsamp_));
  }
  else if (sumforcev_>sumforceu_ && sumforcev_>sumforcew_)
  {
    flowdirection=1;
    ltau = visc_/sqrt(sumforcev_/(area*numsamp_));
  }
  else if (sumforcew_>sumforceu_ && sumforcew_>sumforcev_)
  {
    flowdirection=2;
    ltau = visc_/sqrt(sumforcew_/(area*numsamp_));
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
    s.append(".flow_statistic");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# Flow statistics for turbulent flow in a channel (first- and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    (*log) << "# (u_tau)^2 = tau_W/rho : ";
    (*log) << "   " << setw(11) << setprecision(4) << sumforceu_/(area*numsamp_);
    (*log) << "   " << setw(11) << setprecision(4) << sumforcev_/(area*numsamp_);
    (*log) << "   " << setw(11) << setprecision(4) << sumforcew_/(area*numsamp_);
    (*log) << &endl;

    (*log) << "#     y            y+";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "        mean u^2      mean v^2      mean w^2";
    (*log) << "      mean u*v      mean u*w      mean v*w        Varp   \n";

    (*log) << scientific;
    for(unsigned i=0; i<planecoordinates_->size(); ++i)
    {
      (*log) <<  " "  << setw(11) << setprecision(4) << (*planecoordinates_)[i];
      (*log) << "   " << setw(11) << setprecision(4) << (*planecoordinates_)[i]/ltau;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumu_)  [i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumv_)  [i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumw_)  [i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sump_)  [i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumsqu_)[i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumsqv_)[i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumsqw_)[i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumuv_) [i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumuw_) [i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumvw_) [i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*sumsqp_)[i]/aux;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumu_)  [i]/numsamp_;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumv_)  [i]/numsamp_;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumw_)  [i]/numsamp_;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsump_)  [i]/numsamp_;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumsqu_)[i]/numsamp_;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumsqv_)[i]/numsamp_;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumsqw_)[i]/numsamp_;
      (*log) << "   " << setw(11) << setprecision(4) << (*pointsumsqp_)[i]/numsamp_;
      (*log) << "   \n";
    }
    log->flush();
  }

  if (discret_->Comm().MyPID()==0)
  {

    // additional output for dynamic Smagorinsky model
    if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
        ==
        "CLASSICAL_LES")
    {
      if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Dynamic_Smagorinsky"
         ||
         params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
         ==
         "Smagorinsky_with_van_Driest_damping"
        )
      {
        // get the outfile
        Teuchos::RefCountPtr<std::ofstream> log_Cs;

        std::string s_smag = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
        s_smag.append(".Cs_statistic");

        log_Cs = Teuchos::rcp(new std::ofstream(s_smag.c_str(),ios::out));
        (*log_Cs) << "# Smagorinsky parameter statistics for turbulent flow in a channel";
        (*log_Cs) << "\n\n\n";
        (*log_Cs) << "# Statistics record ";
        (*log_Cs) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";


        (*log_Cs) << "#     y      ";
        (*log_Cs) << "     Cs     ";
        (*log_Cs) << "   (Cs*hk)^2 ";
        (*log_Cs) << "    visceff  ";
        (*log_Cs) << &endl;
        (*log_Cs) << scientific;
        for (unsigned rr=0;rr<sumCs_->size();++rr)
        {
          (*log_Cs) << setw(11) << setprecision(4) << 0.5*((*nodeplanes_)[rr+1]+(*nodeplanes_)[rr]) << "  " ;
          (*log_Cs) << setw(11) << setprecision(4) << ((*sumCs_)[rr])/(numele_*numsamp_) << "  ";
          (*log_Cs) << setw(11) << setprecision(4) << ((*sumCs_delta_sq_)[rr])/(numele_*numsamp_)<< "  " ;
          (*log_Cs) << setw(11) << setprecision(4) << ((*sumvisceff_)[rr])/(numele_*numsamp_) << &endl;
        }
        log_Cs->flush();
      }
    }
  }

  if (discret_->Comm().MyPID()==0)
  {
    Teuchos::RefCountPtr<std::ofstream> log_res;

    // output of residuals and subscale quantities
    std::string s_res = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s_res.append(".res_statistic");

    log_res = Teuchos::rcp(new std::ofstream(s_res.c_str(),ios::out));
    (*log_res) << "# Residual statistics for turbulent flow in a channel";
    (*log_res) << "\n\n\n";
    (*log_res) << "# Statistics record ";
    (*log_res) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";
    (*log_res) << "#       y    ";
    (*log_res) << "    res_x  ";
    (*log_res) << "      res_y  ";
    (*log_res) << "      res_z  ";
    (*log_res) << "     sacc_x  ";
    (*log_res) << "     sacc_y  ";
    (*log_res) << "     sacc_z   ";
    (*log_res) << "   res_sq_x  ";
    (*log_res) << "   res_sq_y  ";
    (*log_res) << "   res_sq_z  ";
    (*log_res) << "   sacc_sq_x ";
    (*log_res) << "   sacc_sq_y ";
    (*log_res) << "   sacc_sq_z "<<&endl;

    (*log_res) << scientific;
    for (unsigned rr=0;rr<nodeplanes_->size()-1;++rr)
    {
      (*log_res)  << setw(11) << setprecision(4) << 0.5*((*nodeplanes_)[rr+1]+(*nodeplanes_)[rr]) << "  " ;
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumres_sq_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_sq_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumres_sq_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_sq_)[3*rr  ]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_sq_)[3*rr+1]/(numele_*numsamp_) << "  ";
      (*log_res)  << setw(11) << setprecision(4) << (*sumsacc_sq_)[3*rr+2]/(numele_*numsamp_) << "  ";

      (*log_res)  << &endl;
    }
    log_res->flush();
  }

  return;

}// TurbulenceStatistics::DumpStatistics


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatistics::ClearStatistics()
{
  // reset the number of samples
  numsamp_ =0;

  // reset forces
  sumforceu_=0;
  sumforcev_=0;
  sumforcew_=0;

  // reset integral and pointwise averages
  for(unsigned i=0; i<planecoordinates_->size(); ++i)
  {
    (*sumu_)[i]  =0;
    (*sumv_)[i]  =0;
    (*sumw_)[i]  =0;
    (*sump_)[i]  =0;

    (*sumuv_ )[i]=0;
    (*sumuw_ )[i]=0;
    (*sumvw_ )[i]=0;
    (*sumsqu_)[i]=0;
    (*sumsqv_)[i]=0;
    (*sumsqw_)[i]=0;
    (*sumsqp_)[i]=0;

    (*pointsumu_)[i]  =0;
    (*pointsumv_)[i]  =0;
    (*pointsumw_)[i]  =0;
    (*pointsump_)[i]  =0;

    (*pointsumsqu_)[i]=0;
    (*pointsumsqv_)[i]=0;
    (*pointsumsqw_)[i]=0;
    (*pointsumsqp_)[i]=0;

  }

  meanvelnp_->PutScalar(0.0);

  // reset smapling for dynamic Smagorinsky model
  if (params_.sublist("TURBULENCE MODEL").get<string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES")
      ==
      "CLASSICAL_LES")
  {
    if(params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Dynamic_Smagorinsky"
       ||
       params_.sublist("TURBULENCE MODEL").get<string>("PHYSICAL_MODEL","no_model")
       ==
       "Smagorinsky_with_van_Driest_damping"
      )
    {
      for (unsigned rr=0;rr<sumCs_->size();++rr)
      {
        // reset value
        (*sumCs_)         [rr]=0;
        (*sumCs_delta_sq_)[rr]=0;
        (*sumvisceff_)    [rr]=0;
      }
    }
  }

  // reset residuals and subscale averages
  for (unsigned rr=0;rr<sumres_->size()/3;++rr)
  {
    (*sumres_      )[3*rr  ]=0;
    (*sumres_      )[3*rr+1]=0;
    (*sumres_      )[3*rr+2]=0;

    (*sumsacc_     )[3*rr  ]=0;
    (*sumsacc_     )[3*rr+1]=0;
    (*sumsacc_     )[3*rr+2]=0;

    (*sumsvelaf_   )[3*rr  ]=0;
    (*sumsvelaf_   )[3*rr+1]=0;
    (*sumsvelaf_   )[3*rr+2]=0;

    (*sumres_sq_   )[3*rr  ]=0;
    (*sumres_sq_   )[3*rr+1]=0;
    (*sumres_sq_   )[3*rr+2]=0;

    (*sumsacc_sq_  )[3*rr  ]=0;
    (*sumsacc_sq_  )[3*rr+1]=0;
    (*sumsacc_sq_  )[3*rr+2]=0;

    (*sumsvelaf_sq_)[3*rr  ]=0;
    (*sumsvelaf_sq_)[3*rr+1]=0;
    (*sumsvelaf_sq_)[3*rr+2]=0;
  }
  for (unsigned rr=0;rr<sumresC_->size();++rr)
  {
    (*sumresC_       )[rr]=0;
    (*sumspresacc_   )[rr]=0;
    (*sumspressnp_   )[rr]=0;

    (*sumresC_sq_    )[rr]=0;
    (*sumspresacc_sq_)[rr]=0;
    (*sumspressnp_sq_)[rr]=0;
  }

  return;
}// TurbulenceStatistics::ClearStatistics



#endif /* CCADISCRET       */
