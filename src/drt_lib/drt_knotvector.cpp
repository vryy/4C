/*!----------------------------------------------------------------------
\file drt_knotvector.cpp

<pre>
Maintainer: Peter Gamnitzer
            gammi@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_knotvector.H"
#include <blitz/array.h>


/*----------------------------------------------------------------------*
 |  empty ctor (public)                                      gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::Knotvector() :
  ParObject              (         ),
  dim_                   (0        ),
  npatches_              (0        ),
  filled_                (false    ),
  degree_                (0        ),
  n_x_m_x_l_             (0        ),
  nele_x_mele_x_lele_    (0        ),
  interpolation_         (0        ),
  offsets_               (0        ),
  knot_values_           (0        )
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::Knotvector(
  int                  dim       ,
  int                  npatches
  ) :
  ParObject              (        ),
  dim_                   (dim     ),
  npatches_              (npatches),
  filled_                (false   ),
  degree_                (npatches),
  n_x_m_x_l_             (npatches),
  nele_x_mele_x_lele_    (npatches),
  interpolation_         (npatches),
  offsets_               (npatches),
  knot_values_           (npatches)
{
  // check if there are any patches
  if(npatches_<0)
  {
    dserror("we need at least one patch\n");
  }

  // resize degrees
  
  // loop patches, resize to dimension
  for(int rr=0;rr<npatches_;++rr)
  {
    (degree_[rr]).resize(dim_);
  }

  // resize n_x_m_x_l, 
  
  // loop patches, resize to dimension
  for(int rr=0;rr<npatches_;++rr)
  {
    (n_x_m_x_l_         [rr]).resize(dim_);
    (nele_x_mele_x_lele_[rr]).resize(dim_);
  }

  // initialise interpolation on all patches
  for(int rr=0;rr<npatches_;++rr)
  {
    (interpolation_[rr]).resize(dim_);
    for(int mm=0;mm<dim_;++mm)
    {
      (interpolation_[rr])[mm]=knotvector_is_not_defined;
    }
  }  

  // provide knot vectors for all patches and dimensions
  for(int rr=0;rr<npatches_;++rr)
  {
    (knot_values_[rr]).resize(dim_);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::~Knotvector()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy ctor (public)                                       gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::Knotvector(const DRT::NURBS::Knotvector & old)
:
  ParObject     (old               ),
  dim_          (old.dim_          ),
  npatches_     (old.npatches_     ),
  filled_       (old.filled_       ),
  degree_       (old.degree_       ),
  n_x_m_x_l_    (old.n_x_m_x_l_    ),
  interpolation_(old.interpolation_),
  offsets_      (old.offsets_      ),
  knot_values_  (old.npatches_     )
{
  // deep copy knot vectors

  for(int np=0;np<npatches_;++np)
  {
    (knot_values_[np]).resize(dim_);
    for(int rr=0;rr<dim_;++rr)
    {
      ((knot_values_[np])[rr]) = Teuchos::rcp(new vector<double>);
      *((knot_values_[np])[rr])= *((old.knot_values_[np])[rr]);
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | convert an element gid to its corresponding triple knot index        |
 |                                                  (public) gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::ConvertEleGidToKnotIds(
  const int      gid        ,
  int         &  npatch     ,
  vector<int> &  loc_cart_id)
{

  if((int)loc_cart_id.size()!= dim_)
  {
    dserror("size vector not of appropriate size (%d,%d)\n",(int)loc_cart_id.size(),dim_);
  }

  if(filled_==false)
  {
    dserror("cannot convert ele ids when filled is false\n");
  }

  // gid is at least in patch 0 (or higher)
  npatch=0;

  for(int np=1;np<npatches_;++np)
  {
    // if this is true, gid is in this patch (or higher)
    if(gid>=offsets_[np])
    {
      npatch++;
    }
    else
    {
      break;
    }
  }

  // reduce gid by patchoffset to get patch local id
  const int locid = gid-offsets_[npatch];

  if(dim_==3)
  {
    // locid = num_u+num_v*nele+num_w*nele*mele     (3d)
    //         |              |       |       |
    //         +--------------+       +-------+
    //            inthislayer          uv_layer
    int uv_layer   = 
      (nele_x_mele_x_lele_[npatch])[0]
      *
      (nele_x_mele_x_lele_[npatch])[1]; 

    // compute num_w
    loc_cart_id[2]   = locid/uv_layer;

    // see above
    int inthislayer= locid%uv_layer;

    // compute num_v and num_u
    loc_cart_id[1]   = inthislayer/(nele_x_mele_x_lele_[npatch])[0];
    loc_cart_id[0]   = inthislayer%(nele_x_mele_x_lele_[npatch])[0];
  }
  else if (dim_==2)
  {
    // locid = num_u+num_v*nele                     (2d)

    // compute num_v and num_u
    loc_cart_id[0]   = locid%(nele_x_mele_x_lele_[npatch])[0];
    loc_cart_id[1]   = locid/(nele_x_mele_x_lele_[npatch])[0];
  }
  else
  {
    dserror("dim_ not available\n");
  }

  return;
}


/*----------------------------------------------------------------------*
 | get element knot vectors to a given element id   (public) gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::GetEleKnots(
  vector<blitz::Array<double,1> > & eleknots,
  int                               gid
  )
{  
  //------------------------------------------------
  // determine the segments knot values 
  // --- remember, nurbs are a cartesian thing, 
  // that means there is a matching between control 
  // point ids and knot ids ....

  // this is the number of knots associated with  
  // this specific element: 
  //
  //         +----------------+
  //         |                |
  //         | (2*degree_u+2) |
  //         | (2*degree_v+2) |
  //         | (2*degree_w+2) |
  //         |                |
  //         +----------------+

 
  if(filled_==false)
  {
    dserror("cannot get ele knots when filled is false\n");
  }

  eleknots.resize(dim_);

  // get base indices and patch number from element gid
  vector<int> cartids(dim_);
  int         npatch;

  ConvertEleGidToKnotIds(gid,npatch,cartids);

  // use them to aquire the required knots
  for(int rr=0;rr<dim_;++rr)
  {
    (eleknots[rr]).resize(2*(degree_[npatch])[rr]+2);
    
    for(int mm=0;mm<2*(degree_[npatch])[rr]+2;++mm)
    {
      (eleknots[rr])(mm)=(*((knot_values_[npatch])[rr]))[cartids[rr]+mm];
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | set knots in one direction                       (public) gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::SetKnots(
  const int                     & direction       , 
  const int                     & npatch          , 
  const int                     & degree          , 
  const int                     & numknots        , 
  const std::string             & knotvectortype  ,
  Teuchos::RCP<vector<double> >   directions_knots)
{

  // filled is false now since new add new knots
  filled_=false;

  if(direction<0 || direction>dim_-1)
  {
    dserror("direction has to in[0...dim_]\n");
  }

  if(npatch<0 || npatch>npatches_-1)
  {
    dserror("patchnumber is invalid\n");
  }

  // set the type
  if(knotvectortype=="Interpolated")
  {
    (interpolation_[npatch])[direction]=knotvector_is_interpolating;
  }
  else if (knotvectortype=="Periodic")
  {
    (interpolation_[npatch])[direction]=knotvector_is_periodic;
  }
  else
  {
    dserror("unknown knotvector-type '%s'\n",knotvectortype.c_str());
  }

  // set the degree of the added knotvector
  (degree_   [npatch])[direction]=degree;

  // set the size of the added knotvector
  (n_x_m_x_l_[npatch])[direction]=numknots;

  // set the actual values
  (knot_values_[npatch])[direction]=directions_knots;

  return;
}

/*----------------------------------------------------------------------*
 | finish                                           (public) gammi 05/08|
 *----------------------------------------------------------------------*/ 
void DRT::NURBS::Knotvector::FinishKnots()
{
  //--------------------------------------------------
  // plausibility checks

  // check if there are any patches
  if(npatches_<0)
  {
    dserror("we need at least one patch\n");
  }

  // check degrees
  if((int)degree_.size()!=npatches_)
  {
    dserror("each patch needs its own degree information\n");
  }
  else
  {
    // loop patches, check dimensions
    for(int rr=0;rr<npatches_;++rr)
    {
      if((int)(degree_[rr]).size()!=dim_)
      {
	dserror("size mismatch: degree\n");
      }
    }
  }

  // check n_x_m_x_ls
  if((int)n_x_m_x_l_.size()!=npatches_)
  {
    dserror("each patch needs its own n_x_m_x_l information\n");
  }
  else
  {
    // loop patches, check dimensions
    for(int rr=0;rr<npatches_;++rr)
    {
      if((int)(n_x_m_x_l_[rr]).size()!=dim_)
      {
	dserror("size mismatch: n_x_m_x_l\n");
      }
    }
  }

  // do we have a knotvector for each dimension 
  // and each patch?
  if((int)knot_values_.size()!=npatches_)
  {
    dserror("each patch needs its own knotvector\n");
  }

  // loop patches
  for(int np=0;np<npatches_;++np)
  {
    if ((int)(knot_values_[np]).size()!=dim_)
    {
      dserror("knotvector of patch has to be of size dim\n");
    }

    for(int rr=0;rr<dim_;++rr)
    {
      // is the knotvector of this dimension nonempty?
      if((knot_values_[np])[rr]==Teuchos::null)
      {
	dserror("no knotvector available in this direction\n");
      }

      // has it the correct size?
      if((int)(*((knot_values_[np])[rr])).size()
	 !=
	 (n_x_m_x_l_[np])[rr])
      {
	dserror("knotvector size mismatch to n_x_m_x_l_ %d!=%d\n",
		(*((knot_values_[np])[rr])).size(),
		(n_x_m_x_l_[np])[rr]);
      }
    
      // is interpolation/periodicity assigned correctly?
      if((interpolation_[np])[rr]==knotvector_is_not_defined)
      {
	dserror("undefined knotvector type\n");
      }
      else if((interpolation_[np])[rr]==knotvector_is_interpolating)
      {
	// for interpolating knot vectors, the first and last
	// knots have to be repeated degree+1 times
	double firstval = (*((knot_values_[np])[rr]))[                     0];
	double lastval  = (*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr]-1];

	for(int mm=1;mm<(degree_[np])[rr]+1;++mm)
	{
	  double db =
	    abs((*((knot_values_[np])[rr]))[                       mm]-firstval);
	  double de = 
	    abs((*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr]-1-mm]-lastval );
	  
	  if(de>1e-9||db>1e-9)
	  {
	    dserror("need multiple knots at the beginning and end of an interpolated knotvector\n");
	  }
	}
      }
      else if((interpolation_[np])[rr]==knotvector_is_periodic)
      {
	// for periodic knot vectors, distances between the 
	// degree+1 first and last nodes have to be equal
	for(int mm=1;mm<(degree_[np])[rr]+1;++mm)
	{
	  double db = 
	    (*((knot_values_[np])[rr]))[mm  ]
	    -
	    (*((knot_values_[np])[rr]))[mm-1];
	  double de = 
	    (*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr]  -mm]
	    -
	    (*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr]-1-mm];
	  
	  if(abs(de-db)>1e-9)
	  {
	    dserror("periodic knotvector doesn't obey periodicity\n");
	  }
	}
      }
    } // loop dimensions
  } // end loop patches

  //--------------------------------------------------
  // generate offset arrays for element to patch
  // mapping and size of element arrays of patches

  // get patches element distribution
  for(int rr=0;rr<npatches_;++rr)
  {
    for(int mm=0;mm<dim_;++mm)
    {
      (nele_x_mele_x_lele_[rr])[mm]=(n_x_m_x_l_[rr])[mm]-2*(degree_[rr])[mm]-1;
    }
  }

  // get element ordering among patches
  offsets_[0]=0;
  for(int rr=1;rr<npatches_;++rr)
  {
    int nele_inpatch=1;
    for(int mm=0;mm<dim_;++mm)
    {
      nele_inpatch*=(nele_x_mele_x_lele_[rr-1])[mm];
    }
    offsets_[rr]=offsets_[rr-1]+nele_inpatch;
  }

  //--------------------------------------------------
  // the knotvector is OK
  filled_=true;

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 05/08 |
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add number of patches
  AddtoPack(data,npatches_);

  // add dimension
  AddtoPack(data,dim_);

  // add degree vector  
  for(int np=0;np<npatches_;++np)
  {
    AddtoPack(data,degree_[np]);
  }

  // add knotvector size
  for(int np=0;np<npatches_;++np)
  {
    AddtoPack(data,n_x_m_x_l_[np]);
  }

  // add element numbers in all cartesian
  // directions
  for(int np=0;np<npatches_;++np)
  {
    AddtoPack(data,nele_x_mele_x_lele_[np]);
  }

  // add Knotvector types
  for(int np=0;np<npatches_;++np)
  {
    for(int rr=0;rr<dim_;++rr)
    {
      AddtoPack(data,(interpolation_[np])[rr]);
    }
  }

  // add patch offsets
  AddtoPack(data,offsets_);

  // add Knotvector coordinates itself
  for(int np=0;np<npatches_;++np)
  {
    for(int rr=0;rr<dim_;++rr)
    {
      AddtoPack(data,(*((knot_values_[np])[rr])));
    }
  }

  return;
}



/*----------------------------------------------------------------------*
 |  Unpack Knotvectors data                                    (public) |
 |                                                          gammi 05/08 |
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::Unpack(const vector<char>& data)
{
  int position = 0;

  filled_=false;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract number of patches
  ExtractfromPack(position,data,npatches_);

  // extract dimension
  ExtractfromPack(position,data,dim_);

  // resize all vectors 
  degree_            .resize(npatches_);
  n_x_m_x_l_         .resize(npatches_);
  nele_x_mele_x_lele_.resize(npatches_);
  interpolation_     .resize(npatches_);
  knot_values_       .resize(npatches_);

  for(int np=0;np<npatches_;++np)
  {
    (degree_            [np]).resize(dim_);
    (n_x_m_x_l_         [np]).resize(dim_);
    (nele_x_mele_x_lele_[np]).resize(dim_);
    (interpolation_     [np]).resize(dim_);
    (knot_values_       [np]).resize(dim_);
  }
  
  // extract degree vector 
  for(int np=0;np<npatches_;++np)
  {
    ExtractfromPack(position,data,degree_[np]);
  }

  // extract knotvector size
  for(int np=0;np<npatches_;++np)
  {
    ExtractfromPack(position,data,n_x_m_x_l_[np]);
  }
    
  // extract element numbers in all cartesian
  // directions
  for(int np=0;np<npatches_;++np)
  {
    ExtractfromPack(position,data,nele_x_mele_x_lele_[np]);
  }

  // extract knotvector types
  for(int np=0;np<npatches_;++np)
  {
    for(int rr=0;rr<dim_;++rr)
    {
      ExtractfromPack(position,data,((interpolation_[np])[rr]));
    }
  }

  // extract patch offsets
  ExtractfromPack(position,data,offsets_);
  
  // extract knotvector coordinates itself
  for(int np=0;np<npatches_;++np)
  {
    for(int rr=0;rr<dim_;++rr)
    {
      (knot_values_[np])[rr]
	=
	Teuchos::rcp(new vector<double>((n_x_m_x_l_[np])[rr]));
    
      ExtractfromPack(position,data,(*((knot_values_[np])[rr])));
    }
  }

  return;
}

#endif
