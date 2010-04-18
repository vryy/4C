/*!----------------------------------------------------------------------
  \file turbulence_statistics_ccy.cpp

\brief Compute (time and space) averaged values for turbulent flows
       around a rotating cylinder and write them to files.


*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "turbulence_statistics_ccy.H"

/*----------------------------------------------------------------------

                  Standard Constructor (public)

  ---------------------------------------------------------------------*/
FLD::TurbulenceStatisticsCcy::TurbulenceStatisticsCcy(
  RefCountPtr<DRT::Discretization> actdis             ,
  bool                             alefluid           ,
  RefCountPtr<Epetra_Vector>       dispnp             ,
  ParameterList&                   params
  )
  :
  discret_            (actdis             ),
  alefluid_           (alefluid           ),
  dispnp_             (dispnp             ),
  params_             (params             )
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
  {
    dserror("Evaluation of turbulence statistics only for 3d flows!");
  }

  //----------------------------------------------------------------------
  // allocate some vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  meanvelnp_     = LINALG::CreateVector(*dofrowmap,true);

  //----------------------------------------------------------------------
  // switches, control parameters, material parameters

  // get the plane normal direction from the parameterlist
  {
    string planestring = params_.sublist("TURBULENCE MODEL").get<string>("HOMDIR","not_specified");

    if(planestring == "z")
    {
      dim_ = 2;
    }
    else
    {
      dserror("homogeneuous direction for this flow was specified incorrectly. (need z)");
    }
  }

  // ---------------------------------------------------------------------
  // up to now, there are no records written
  countrecord_ = 0;

  // ---------------------------------------------------------------------
  // compute all planes for sampling

  // available shells of element corners (Nurbs) of elements
  nodeshells_ = rcp(new vector<double> );

  // available homogeneous (sampling) shells --- there are
  // numsubdivisions layers per element layer
  shellcoordinates_ = rcp(new vector<double> );

  const int numsubdivisions=5;

  // try to cast discretisation to nurbs variant
  // this tells you what kind of computation of
  // samples is required
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*actdis));

  if(nurbsdis==NULL)
  {
    dserror("Need Nurbs mesh for turbulent flows around a circular cylinder\n");
  }
  else
  {

    // real pointwise control point sampling does not make any sense
    // for Nurbs discretisations since shape functions are not interpolating

    // radial shellcoordinates are determined by the element
    // (cartesian) number in the second knotspan direction and
    // the number of sampling shells in between are added

    // for nurbs discretisations, all vector sizes are already determined
    // by the knotvector size

    // get nurbs dis' knotvector sizes
    vector<int> n_x_m_x_l(nurbsdis->Return_n_x_m_x_l(0));

    // get nurbs dis' element numbers
    vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(0));

    // get the knotvector itself
    RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

    // resize and initialise to 0
    {
      (*nodeshells_      ).resize(nele_x_mele_x_lele[1]+1);
      (*shellcoordinates_).resize(nele_x_mele_x_lele[1]*(numsubdivisions-1)+1);

      vector<double>::iterator coord;

      for (coord  = (*nodeshells_).begin();
	   coord != (*nodeshells_).end()  ;
	   ++coord)
      {
	*coord=0;
      }
      for (coord  = shellcoordinates_->begin();
	   coord != shellcoordinates_->end()  ;
	   ++coord)
      {
	*coord=0;
      }
    }

    // count numbers of nodes in homogeneous shells
    int nodeshellsize=nodeshells_->size();

    vector<int> nodeshells_numnodes      (nodeshellsize,0);
    vector<int> shellcoordinates_numnodes((*shellcoordinates_).size(),0);

    // get element map
    const Epetra_Map* elementmap = nurbsdis->ElementRowMap();

    // loop all available elements
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      // get element pointer
      DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));

      // want to loop all control points of the element,
      // so get the number of points
      const int numnp = actele->NumNode();

      // get the elements control points/nodes
      DRT::Node**   nodes = actele->Nodes();

      // aquire weights from nodes
      Epetra_SerialDenseVector weights(numnp);

      for (int inode=0; inode<numnp; ++inode)
      {
	DRT::NURBS::ControlPoint* cp
	  =
	  dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

	weights(inode) = cp->W();
      }

      // get gid, location in the patch
      int gid = actele->Id();

      int patchid=0;

      vector<int> ele_cart_id(3);
      knots->ConvertEleGidToKnotIds(gid,patchid,ele_cart_id);

      // access elements knot span
      std::vector<Epetra_SerialDenseVector> knots(3);
      bool zero_size=(*((*nurbsdis).GetKnotVector())).GetEleKnots(knots,actele->Id());

      // zero sized elements have to be skipped
      if(zero_size)
      {
        continue;
      }

      // get shapefunctions, compute all visualisation point positions
      Epetra_SerialDenseVector nurbs_shape_funct(numnp);

      switch (actele->Shape())
      {
      case DRT::Element::nurbs8:
      case DRT::Element::nurbs27:
      {
	// element local point position
	Epetra_SerialDenseVector uv(3);

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
	  // use r-coordinate of point 2 and 4
	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
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

	  (*nodeshells_      )[ele_cart_id[1]]                    +=sqrt(x[0]*x[0]+x[1]*x[1]);
	  (*shellcoordinates_)[ele_cart_id[1]*(numsubdivisions-1)]+=sqrt(x[0]*x[0]+x[1]*x[1]);

          nodeshells_numnodes      [ele_cart_id[1]                    ]+=1;
          shellcoordinates_numnodes[ele_cart_id[1]*(numsubdivisions-1)]+=1;


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
	    (*shellcoordinates_)[ele_cart_id[1]*(numsubdivisions-1)+rr]+=sqrt(x[0]*x[0]+x[1]*x[1]);
            ++(shellcoordinates_numnodes[ele_cart_id[1]*(numsubdivisions-1)+rr]);
	  }


	  // set upper point of element, too (only for last layer)
	  if(ele_cart_id[1]+1 == nele_x_mele_x_lele[1])
	  {
	    // point 8
	    uv(0)=  1.0;
	    uv(1)=  1.0;
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

	    (*nodeshells_)      [ele_cart_id[1]                      +1]+=sqrt(x[0]*x[0]+x[1]*x[1]);
	    (*shellcoordinates_)[(ele_cart_id[1]+1)*(numsubdivisions-1)]+=sqrt(x[0]*x[0]+x[1]*x[1]);
            ++(nodeshells_numnodes      [ele_cart_id[1]                      +1]);
            ++(shellcoordinates_numnodes[(ele_cart_id[1]+1)*(numsubdivisions-1)]);
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

    std::vector<double> lnodeplanes      (*nodeshells_      );
    std::vector<double> lplanecoordinates(*shellcoordinates_);

    vector<int> lnodeshells_numnodes      (nodeshells_numnodes);
    vector<int> lshellcoordinates_numnodes(shellcoordinates_numnodes);

    discret_->Comm().SumAll(&(lnodeplanes[0]      ),&((*nodeshells_      )[0]),nodeshells_->size()      );
    discret_->Comm().SumAll(&(lplanecoordinates[0]),&((*shellcoordinates_)[0]),shellcoordinates_->size());

    discret_->Comm().SumAll(&(lnodeshells_numnodes      [0]),&(nodeshells_numnodes      [0]),nodeshells_numnodes.size());
    discret_->Comm().SumAll(&(lshellcoordinates_numnodes[0]),&(shellcoordinates_numnodes[0]),shellcoordinates_numnodes.size());

    {
      (*nodeshells_      ).resize(nele_x_mele_x_lele[1]+1);
      (*shellcoordinates_).resize(nele_x_mele_x_lele[1]*(numsubdivisions-1)+1);

      for (unsigned rr=0;rr<(*nodeshells_).size();++rr)
      {
        if(fabs(nodeshells_numnodes[rr])<1e-9)
        {
          dserror("zero nodes in shell layer %d\n",rr);
        }

        (*nodeshells_)[rr]/=(double)(nodeshells_numnodes[rr]);
      }


      for (unsigned rr=0;rr<(*shellcoordinates_).size();++rr)
      {
        if(fabs(shellcoordinates_numnodes[rr])<1e-9)
        {
          dserror("zero nodes in sampling shell layer %d\n",rr);
        }

        (*shellcoordinates_)[rr]/=(double)(shellcoordinates_numnodes[rr]);
      }
    }
  }

  //----------------------------------------------------------------------
  // sort shellcoordinates and nodeshells
  {
    set<double,PlaneSortCriterion> shellset;

    vector<double>::iterator coord;
    set<double,PlaneSortCriterion>::iterator shell;

    {
      for (coord  = (*nodeshells_).begin();
           coord != (*nodeshells_).end()  ;
           ++coord)
      {
        shellset.insert(*coord);
      }

      int rr=0;
      for (shell  = shellset.begin();
           shell != shellset.end()  ;
           ++shell)
      {
        (*nodeshells_)[rr]=*shell;
        ++rr;
      }

    }

    shellset.clear();

    {

      for (coord  = (*shellcoordinates_).begin();
           coord != (*shellcoordinates_).end()  ;
           ++coord)
      {
        shellset.insert(*coord);
      }

      int rr=0;
      for (shell  = shellset.begin();
           shell != shellset.end()  ;
           ++shell)
      {
        (*shellcoordinates_)[rr]=*shell;
        ++rr;
      }
    }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of in plane mean values

  int size = shellcoordinates_->size();

  // arrays for point based averaging
  // --------------------------------

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
  pointsumuu_ =  rcp(new vector<double> );
  pointsumuu_->resize(size,0.0);

  pointsumvv_ =  rcp(new vector<double> );
  pointsumvv_->resize(size,0.0);

  pointsumww_ =  rcp(new vector<double> );
  pointsumww_->resize(size,0.0);

  pointsumpp_ =  rcp(new vector<double> );
  pointsumpp_->resize(size,0.0);

  pointsumuv_  =  rcp(new vector<double> );
  pointsumuv_->resize(size,0.0);

  pointsumuw_  =  rcp(new vector<double> );
  pointsumuw_->resize(size,0.0);

  pointsumvw_  =  rcp(new vector<double> );
  pointsumvw_->resize(size,0.0);

  //----------------------------------------------------------------------
  // initialise output

  Teuchos::RefCountPtr<std::ofstream> log;

  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::out));
    (*log) << "# Statistics for turbulent incompressible flow in a rotating cylinder (first- and second-order moments)\n\n";

    log->flush();
  }

  // clear statistics
  this->ClearStatistics();

  return;
}// TurbulenceStatisticsCcy::TurbulenceStatisticsCcy

/*----------------------------------------------------------------------*

                           Destructor

 -----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsCcy::~TurbulenceStatisticsCcy()
{
  return;
}// TurbulenceStatisticsCcy::~TurbulenceStatisticsCcy()

/*----------------------------------------------------------------------*

       Compute the in-plane mean values of first and second order
       moments for velocities, pressure and Cs are added to global
                            'sum' vectors.

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::DoTimeSample(
  Teuchos::RefCountPtr<Epetra_Vector> velnp
  )
{
  // we have an additional sample
  numsamp_++;

  // meanvelnp is a refcount copy of velnp

  meanvelnp_->Update(1.0,*velnp,0.0);

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane
  this->EvaluatePointwiseMeanValuesInPlanes();

  return;
}// TurbulenceStatisticsCcy::DoTimeSample


/*----------------------------------------------------------------------*

          Compute in plane means of u,u^2 etc. (nodal quantities)

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::EvaluatePointwiseMeanValuesInPlanes()
{

  const int numsubdivisions=5;

  //----------------------------------------------------------------------
  // sort shellcoordinates and nodeshells
  map<double,int   ,PlaneSortCriterion> countpoints;
  map<double,double,PlaneSortCriterion> meanu;
  map<double,double,PlaneSortCriterion> meanv;
  map<double,double,PlaneSortCriterion> meanw;
  map<double,double,PlaneSortCriterion> meanp;
  map<double,double,PlaneSortCriterion> meanuu;
  map<double,double,PlaneSortCriterion> meanvv;
  map<double,double,PlaneSortCriterion> meanww;
  map<double,double,PlaneSortCriterion> meanpp;
  map<double,double,PlaneSortCriterion> meanuv;
  map<double,double,PlaneSortCriterion> meanuw;
  map<double,double,PlaneSortCriterion> meanvw;

  for (vector<double>::iterator coord  = (*shellcoordinates_).begin();
       coord != (*shellcoordinates_).end()  ;
       ++coord)
  {
    double r=*coord;

    meanu.insert(pair<double,double>(r,0.0));
    meanv.insert(pair<double,double>(r,0.0));
    meanw.insert(pair<double,double>(r,0.0));
    meanp.insert(pair<double,double>(r,0.0));

    meanuu.insert(pair<double,double>(r,0.0));
    meanvv.insert(pair<double,double>(r,0.0));
    meanww.insert(pair<double,double>(r,0.0));
    meanpp.insert(pair<double,double>(r,0.0));
    meanuv.insert(pair<double,double>(r,0.0));
    meanuw.insert(pair<double,double>(r,0.0));
    meanvw.insert(pair<double,double>(r,0.0));

    countpoints.insert(pair<double,int>(r,0));
  }

  // try to cast discretisation to nurbs variant
  // this tells you what kind of computation of
  // samples is required
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));

  nurbsdis->SetState("velnp",meanvelnp_);

  // get nurbs dis' knotvector sizes
  vector<int> n_x_m_x_l(nurbsdis->Return_n_x_m_x_l(0));

  // get nurbs dis' element numbers
  vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(0));

  // get the knotvector itself
  RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

  // get element map
  const Epetra_Map* elementmap = nurbsdis->ElementRowMap();

  // loop all available elements
  for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
  {
    // get element pointer
    DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));

    // want to loop all control points of the element,
    // so get the number of points
    const int numnp = actele->NumNode();

    // get the elements control points/nodes
    DRT::Node**   nodes = actele->Nodes();

    // aquire weights from nodes
    Epetra_SerialDenseVector weights(numnp);

    for (int inode=0; inode<numnp; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

	weights(inode) = cp->W();
    }
    // get gid, location in the patch
    int gid = actele->Id();

    int patchid=0;

    vector<int> ele_cart_id(3);
    knots->ConvertEleGidToKnotIds(gid,patchid,ele_cart_id);

    // access elements knot span
    std::vector<Epetra_SerialDenseVector> knots(3);
    bool zero_size=(*((*nurbsdis).GetKnotVector())).GetEleKnots(knots,actele->Id());

    // zero sized elements have to be skipped
    if(zero_size)
    {
      continue;
    }

    // get shapefunctions, compute all visualisation point positions
    Epetra_SerialDenseVector nurbs_shape_funct(numnp);

    // extract local values from the global vectors
    vector<int> lm;
    vector<int> lmowner;

    actele->LocationVector(*nurbsdis,lm,lmowner);

    // extract local values from global vector
    vector<double> myvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*(nurbsdis->GetState("velnp")),myvelnp,lm);

    // create Matrix objects
    LINALG::Matrix<3,27> evelnp;
    LINALG::Matrix<27,1> eprenp;

    // insert velocity  into element array
    for (int i=0;i<27;++i)
    {
      const int fi=4*i;

      evelnp(0,i) = myvelnp[  fi];
      evelnp(1,i) = myvelnp[1+fi];
      evelnp(2,i) = myvelnp[2+fi];

      eprenp(  i) = myvelnp[3+fi];
    }

    switch (actele->Shape())
    {
    case DRT::Element::nurbs27:
    {
      LINALG::Matrix<3,1> vel;

      // element local point position
      Epetra_SerialDenseVector uv(3);

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
        // use r-coordinate of point 1 and 8
        // temporary x vector
        std::vector<double> x(3);

        // point 2
        uv(0)=  1.0;
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

        const double r=sqrt(x[0]*x[0]+x[1]*x[1]);

        {
          double val=0;
          for (int inode=0; inode<numnp; ++inode)
          {
            val+=nurbs_shape_funct(inode)*evelnp(0,inode);
          }
          vel(0)=val;

          val=0;
          for (int inode=0; inode<numnp; ++inode)
          {
            val+=nurbs_shape_funct(inode)*evelnp(1,inode);
          }
          vel(1)=val;

          val=0;
          for (int inode=0; inode<numnp; ++inode)
          {
            val+=nurbs_shape_funct(inode)*evelnp(2,inode);
          }
          vel(2)=val;

          val=0;
          for (int inode=0; inode<numnp; ++inode)
          {
            val+=nurbs_shape_funct(inode)*eprenp(inode);
          }
          meanp [r]+=val;
          meanpp[r]+=val*val;

          map<double,int,PlaneSortCriterion>::iterator shell=countpoints.find(r);
          if(shell==countpoints.end())
          {
            dserror("radial coordinate %12.5e was not map\n",r);
          }
          else
          {
            shell->second+=1;
          }

          double uphi=1.0/r*(x[0]*vel(1)-x[1]*vel(0));
          double ur  =1.0/r*(x[0]*vel(0)+x[1]*vel(1));

          meanu[r]+=uphi;
          meanv[r]+=ur;
          meanw[r]+=vel(2);

          meanuu[r]+=uphi*uphi;
          meanvv[r]+=ur*ur;
          meanww[r]+=vel(2)*vel(2);

          meanuv[r]+=uphi*ur;
          meanuw[r]+=uphi*vel(2);
          meanvw[r]+=ur*vel(2);
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

          const double r=sqrt(x[0]*x[0]+x[1]*x[1]);

          {
            double val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*evelnp(0,inode);
            }
            vel(0)=val;

            val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*evelnp(1,inode);
            }
            vel(1)=val;

            val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*evelnp(2,inode);
            }
            vel(2)=val;

            val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*eprenp(inode);
            }
            meanp[r]+=val;
            meanpp[r]+=val*val;

            map<double,int,PlaneSortCriterion>::iterator shell=countpoints.find(r);
            if(shell==countpoints.end())
            {
              dserror("radial coordinate %12.5e was not map\n",r);
            }
            else
            {
              shell->second+=1;
            }

            double uphi=1.0/r*(x[0]*vel(1)-x[1]*vel(0));
            double ur  =1.0/r*(x[0]*vel(0)+x[1]*vel(1));

            meanu[r]+=uphi;
            meanv[r]+=ur;
            meanw[r]+=vel(2);

            meanuu[r]+=uphi*uphi;
            meanvv[r]+=ur*ur;
            meanww[r]+=vel(2)*vel(2);

            meanuv[r]+=uphi*ur;
            meanuw[r]+=uphi*vel(2);
            meanvw[r]+=ur*vel(2);
          }
        }

        // set upper point of element, too (only for last layer)
        if(ele_cart_id[1]+1 == nele_x_mele_x_lele[1])
        {
          // point 4
          uv(0)=  1.0;
          uv(1)=  1.0;
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


          const double r=sqrt(x[0]*x[0]+x[1]*x[1]);

          {
            double val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*evelnp(0,inode);
            }
            vel(0)=val;

            val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*evelnp(1,inode);
            }
            vel(1)=val;

            val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*evelnp(2,inode);
            }
            vel(2)=val;

            val=0;
            for (int inode=0; inode<numnp; ++inode)
            {
              val+=nurbs_shape_funct(inode)*eprenp(inode);
            }
            meanp[r]+=val;
            meanpp[r]+=val*val;

            map<double,int,PlaneSortCriterion>::iterator shell=countpoints.find(r);
            if(shell==countpoints.end())
            {
              dserror("radial coordinate %12.5e was not map\n",r);
            }
            else
            {
              shell->second+=1;
            }

            double uphi=1.0/r*(x[0]*vel(1)-x[1]*vel(0));
            double ur  =1.0/r*(x[0]*vel(0)+x[1]*vel(1));

            meanu[r]+=uphi;
            meanv[r]+=ur;
            meanw[r]+=vel(2);

            meanuu[r]+=uphi*uphi;
            meanvv[r]+=ur*ur;
            meanww[r]+=vel(2)*vel(2);

            meanuv[r]+=uphi*ur;
            meanuw[r]+=uphi*vel(2);
            meanvw[r]+=ur*vel(2);
          }
        }
      }
      break;
    }
    default:
      dserror("Unknown element shape for a nurbs element or nurbs type not valid for turbulence calculation\n");
    }
  } // end element loop

  // communicate results among processors
  int size=countpoints.size();
  int rr;

  // collect number of samples
  vector<int> lpointcount;
  vector<int> pointcount(size);

  for (map<double,int,PlaneSortCriterion>::iterator shell  = countpoints.begin();
       shell != countpoints.end();
       ++shell)
  {
    lpointcount.push_back(shell->second);
  }
  discret_->Comm().SumAll(&(lpointcount[0]),&(pointcount[0]),size);

  // collect number of samples
  vector<double> lmeanu;
  vector<double> lmeanv;
  vector<double> lmeanw;
  vector<double> lmeanp;

  vector<double> lmeanuu;
  vector<double> lmeanvv;
  vector<double> lmeanww;
  vector<double> lmeanpp;

  vector<double> lmeanuv;
  vector<double> lmeanuw;
  vector<double> lmeanvw;

  vector<double> gmeanu(size);
  vector<double> gmeanv(size);
  vector<double> gmeanw(size);
  vector<double> gmeanp(size);

  vector<double> gmeanuu(size);
  vector<double> gmeanvv(size);
  vector<double> gmeanww(size);
  vector<double> gmeanpp(size);

  vector<double> gmeanuv(size);
  vector<double> gmeanuw(size);
  vector<double> gmeanvw(size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanu.begin();
       shell != meanu.end();
       ++shell)
  {
    if(fabs(pointcount[rr])<1e-6)
    {
      dserror("zero pointcount during computation of averages, layer %d\n",rr);
    }

    lmeanu.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanu[0]),&((gmeanu)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanv.begin();
       shell != meanv.end();
       ++shell)
  {
    lmeanv.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanv[0]),&((gmeanv)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanw.begin();
       shell != meanw.end();
       ++shell)
  {
    lmeanw.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanw[0]),&((gmeanw)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanp.begin();
       shell != meanp.end();
       ++shell)
  {
    lmeanp.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanp[0]),&((gmeanp)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanuu.begin();
       shell != meanuu.end();
       ++shell)
  {

    lmeanuu.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanuu[0]),&((gmeanuu)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanvv.begin();
       shell != meanvv.end();
       ++shell)
  {
    lmeanvv.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanvv[0]),&((gmeanvv)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanww.begin();
       shell != meanww.end();
       ++shell)
  {
    lmeanww.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanww[0]),&((gmeanww)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanpp.begin();
       shell != meanpp.end();
       ++shell)
  {
    lmeanpp.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanpp[0]),&((gmeanpp)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanuv.begin();
       shell != meanuv.end();
       ++shell)
  {

    lmeanuv.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanuv[0]),&((gmeanuv)[0]),size);

  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanuw.begin();
       shell != meanuw.end();
       ++shell)
  {
    lmeanuw.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanuw[0]),&((gmeanuw)[0]),size);


  rr=0;
  for (map<double,double,PlaneSortCriterion>::iterator shell  = meanvw.begin();
       shell != meanvw.end();
       ++shell)
  {
    lmeanvw.push_back(shell->second/pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanvw[0]),&(gmeanvw[0]),size);

  for(int mm=0;mm<size;++mm)
  {
    (*pointsumu_ )[mm]+=gmeanu [mm];
    (*pointsumv_ )[mm]+=gmeanv [mm];
    (*pointsumw_ )[mm]+=gmeanw [mm];
    (*pointsump_ )[mm]+=gmeanp [mm];

    (*pointsumuu_)[mm]+=gmeanuu[mm];
    (*pointsumvv_)[mm]+=gmeanvv[mm];
    (*pointsumww_)[mm]+=gmeanww[mm];
    (*pointsumpp_)[mm]+=gmeanpp[mm];

    (*pointsumuv_)[mm]+=gmeanuv[mm];
    (*pointsumuw_)[mm]+=gmeanuw[mm];
    (*pointsumvw_)[mm]+=gmeanvw[mm];
  }

  return;
}// TurbulenceStatisticsCcy::EvaluatePointwiseMeanValuesInPlanes()


/*----------------------------------------------------------------------*

       Compute a time average of the mean values over all steps
          since the last output. Dump the result to file.

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::TimeAverageMeansAndOutputOfStatistics(int step)
{
  if (numsamp_ == 0)
  {
    dserror("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RefCountPtr<std::ofstream> log;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s = params_.sublist("TURBULENCE MODEL").get<string>("statistics outfile");
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),ios::app));
    (*log) << "\n\n\n";
    (*log) << "# Statistics record " << countrecord_;
    (*log) << " (Steps " << step-numsamp_+1 << "--" << step <<")\n";

    (*log) << "#|-------------------";
    (*log) << "---------------------------------------------------";
    (*log) << "--point based (interpolated)-------------------------";
    (*log) << "----------------------------------------------------------|";
    (*log) << "\n";

    (*log) << "#      y       ";
    (*log) << "    u_theta    ";
    (*log) << "      u_r      ";
    (*log) << "      u_z      ";
    (*log) << "       p       ";
    (*log) << "u_theta*u_theta";
    (*log) << "    u_r*u_r    ";
    (*log) << "    u_z*u_z    ";
    (*log) << "      p*p      ";
    (*log) << "  u_theta*u_r  ";
    (*log) << "  u_theta*u_z  ";
    (*log) << "    u_r*u_z    ";
    (*log) << "\n";
    (*log) << scientific;
    for(unsigned i=0; i<shellcoordinates_->size(); ++i)
    {
      // y and y+
      (*log) <<  " "  << setw(11) << setprecision(4) << (*shellcoordinates_)[i];

      // pointwise means
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumu_ )[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumv_ )[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumw_ )[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsump_ )[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumuu_)[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumvv_)[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumww_)[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumpp_)[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumuv_)[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumuw_)[i]/numsamp_;
      (*log) << "    " << setw(11) << setprecision(4) << (*pointsumvw_)[i]/numsamp_;
      (*log) << "\n";
    }
    log->flush();
  } // end myrank 0


  // log was written, so increase counter for records
  countrecord_++;

  return;

}// TurbulenceStatisticsCcy::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*

                  Reset sums and number of samples to 0

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::ClearStatistics()
{
  // reset the number of samples
  numsamp_ =0;

  // reset integral and pointwise averages
  for(unsigned i=0; i<shellcoordinates_->size(); ++i)
  {
    (*pointsumu_)[i]  =0;
    (*pointsumv_)[i]  =0;
    (*pointsumw_)[i]  =0;
    (*pointsump_)[i]  =0;

    (*pointsumuu_)[i]=0;
    (*pointsumvv_)[i]=0;
    (*pointsumww_)[i]=0;
    (*pointsumpp_)[i]=0;

    (*pointsumuv_)[i]=0;
    (*pointsumuw_)[i]=0;
    (*pointsumvw_)[i]=0;
  }

  meanvelnp_->PutScalar(0.0);

  return;
}// TurbulenceStatisticsCcy::ClearStatistics


#endif /* CCADISCRET       */
