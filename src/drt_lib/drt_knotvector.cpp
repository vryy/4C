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
  filled_                (false    ),
  degree_                (0        ),
  n_x_m_x_l_             (0        ),
  nele_x_mele_x_lele_    (0        ),
  interpolation_         (0        ),
  knot_values_           (0        )
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::Knotvector(
  int         dim       ,
  vector<int> degree    ,
  vector<int> n_x_m_x_l
  ) :
  ParObject              (         ),
  dim_                   (dim      ),
  filled_                (false    ),
  degree_                (degree   ),
  n_x_m_x_l_             (n_x_m_x_l),
  nele_x_mele_x_lele_    (n_x_m_x_l),
  interpolation_         (dim_     ),
  knot_values_           (dim_     )
{
  // check dimensions
  if((int)degree_.size()!=dim_)
  {
    dserror("size mismatch: degree\n");
  }
  if((int)n_x_m_x_l_.size()!=dim_)
  {
    dserror("size mismatch: n_x_m_x_l_\n");
  }

  // initialise interpolation
  for(int rr=0;rr<dim_;++rr)
  {
    interpolation_[rr]=knotvector_is_not_defined;
  }
  
  // get element distribution
  for(int rr=0;rr<dim_;++rr)
  {
    nele_x_mele_x_lele_[rr]-=2*degree_[rr]+1;
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
  filled_       (old.filled_       ),
  degree_       (old.degree_       ),
  n_x_m_x_l_    (old.n_x_m_x_l_    ),
  interpolation_(old.interpolation_),
  knot_values_  (old.dim_          )
{
  // deep copy knot vectors
  for(int rr=0;rr<dim_;++rr)
  {
    knot_values_[rr]   = Teuchos::rcp(new vector<double>);
    *(knot_values_[rr])= *(old.knot_values_[rr]);
  }
}


/*----------------------------------------------------------------------*
 | convert an element gid to its corresponding triple knot index        |
 |                                                  (public) gammi 05/08|
 *----------------------------------------------------------------------*/
vector<int> DRT::NURBS::Knotvector::ConvertEleGidToKnotIds(int gid)
{
  
  vector<int> knotindex(dim_);



  if(dim_==3)
  {
    // gid = num_u+num_v*nele+num_w*nele*mele     (3d)
    //       |              |       |       |
    //       +--------------+       +-------+
    //          inthislayer          uv_layer
    int uv_layer   = nele_x_mele_x_lele_[0]*nele_x_mele_x_lele_[1]; 

    // compute num_w
    knotindex[2]   = gid/uv_layer;

    // see above
    int inthislayer= gid%uv_layer;

    // compute num_v and num_u
    knotindex[1]   = inthislayer/nele_x_mele_x_lele_[0];
    knotindex[0]   = inthislayer%nele_x_mele_x_lele_[0];
  }
  else if (dim_==2)
  {
    // gid = num_u+num_v*nele                     (2d)

    // compute num_v and num_u
    knotindex[0]   = gid%nele_x_mele_x_lele_[0];
    knotindex[1]   = gid/nele_x_mele_x_lele_[0];
  }
  else
  {
    dserror("dim_ not available\n");
  }

  return(knotindex);
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

  // get base indices from element gid
  vector<int> cartids(ConvertEleGidToKnotIds(gid));

  // use them to aquire the required knots
  for(int rr=0;rr<dim_;++rr)
  {
    (eleknots[rr]).resize(2*degree_[rr]+2);
    
    for(int mm=0;mm<2*degree_[rr]+2;++mm)
    {
      (eleknots[rr])(mm)=(*(knot_values_[rr]))[cartids[rr]+mm];
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | set knots in one direction                       (public) gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::SetKnots(
  const int                     & direction       , 
  const std::string             & knotvectortype  ,
  Teuchos::RCP<vector<double> > & directions_knots)
{
  // filled is false now since new add new knots
  filled_=false;

  if(direction<0 || direction>dim_-1)
  {
    dserror("direction has to in[0...dim_]\n");
  }

  // set the type
  if(knotvectortype=="Interpolated")
  {
    interpolation_[direction]=knotvector_is_interpolating;
  }
  else if (knotvectortype=="Periodic")
  {
    interpolation_[direction]=knotvector_is_periodic;
  }
  else
  {
    dserror("unknown knotvector-type '%s'\n",knotvectortype.c_str());
  }

  // set the actual values
  knot_values_[direction]=directions_knots;

  return;
}

/*----------------------------------------------------------------------*
 | finish                                           (public) gammi 05/08|
 *----------------------------------------------------------------------*/ 
void DRT::NURBS::Knotvector::FinishKnots()
{
  // do we have a knotvector for each dimension?
  if ((int)knot_values_.size()!=dim_)
  {
    dserror("knotvector has to be of size dim\n");
  }

  for(int rr=0;rr<dim_;++rr)
  {
    // is the knotvector of this dimension nonempty?
    if(knot_values_[rr]==Teuchos::null)
    {
      dserror("no knotvector available in this direction\n");
    }

    // has it the correct size?
    if((int)(*knot_values_[rr]).size()!=n_x_m_x_l_[rr])
    {
      dserror("knotvector size mismatch to n_x_m_x_l_ %d!=%d\n",(*knot_values_[rr]).size(),n_x_m_x_l_[rr]);
    }
    
    // is interpolation/periodicity assigned correctly?
    if(interpolation_[rr]==knotvector_is_not_defined)
    {
      dserror("undefined knotvector type\n");
    }
    else if(interpolation_[rr]==knotvector_is_interpolating)
    {
      // for interpolating knot vectors, the first and last
      // knots have to be repeated degree+1 times
      double firstval = (*knot_values_[rr])[               0];
      double lastval  = (*knot_values_[rr])[n_x_m_x_l_[rr]-1];
      for(int mm=1;mm<degree_[rr]+1;++mm)
      {
	double db = abs((*knot_values_[rr])[mm]-firstval);
	double de = abs((*knot_values_[rr])[n_x_m_x_l_[rr]-1-mm]-lastval);
	
	if(de>1e-9||db>1e-9)
	{
	  dserror("need multiple knots at the beginning and end of an interpolated knotvector\n");
      	}
      }
    }
    else if(interpolation_[rr]==knotvector_is_periodic)
    {
      // for periodic knot vectors, distances between the 
      // degree+1 first and last nodes have to be equal
      for(int mm=1;mm<degree_[rr]+1;++mm)
      {
	double db = (*knot_values_[rr])[mm  ]
	            -
	            (*knot_values_[rr])[mm-1];
	double de = (*knot_values_[rr])[n_x_m_x_l_[rr]  -mm]
		    -
		    (*knot_values_[rr])[n_x_m_x_l_[rr]-1-mm];
	
	if(abs(de-db)>1e-9)
	{
	  dserror("periodic knotvector doesn't obey periodicity\n");
      	}
      }
    }
  }

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

  // add dimension
  AddtoPack(data,dim_);

  // add degree vector  
  AddtoPack(data,degree_);

  // add knotvector size
  AddtoPack(data,n_x_m_x_l_);
    
  // add element numbers in all cartesian
  // directions
  AddtoPack(data,nele_x_mele_x_lele_);

  // add Knotvector types
  for(int rr=0;rr<dim_;++rr)
  {
    AddtoPack(data,interpolation_[rr]);
  }

  // add Knotvector coordinates itself
  for(int rr=0;rr<dim_;++rr)
  {
    AddtoPack(data,(*(knot_values_[rr])));
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

  // extract dimension
  ExtractfromPack(position,data,dim_);

  // resize all vectors 
  degree_            .resize(dim_);
  n_x_m_x_l_         .resize(dim_);
  nele_x_mele_x_lele_.resize(dim_);
  interpolation_     .resize(dim_);
  knot_values_       .resize(dim_);
  
  // extract degree vector  
  ExtractfromPack(position,data,degree_);

  // extract knotvector size
  ExtractfromPack(position,data,n_x_m_x_l_);
    
  // extract element numbers in all cartesian
  // directions
  ExtractfromPack(position,data,nele_x_mele_x_lele_);

  // extract knotvector types
  for(int rr=0;rr<dim_;++rr)
  {
    ExtractfromPack(position,data,(interpolation_[rr]));
  }
  
  // extract knotvector coordinates itself
  for(int rr=0;rr<dim_;++rr)
  {
    knot_values_[rr]=Teuchos::rcp(new vector<double>(n_x_m_x_l_[rr]));
    
    ExtractfromPack(position,data,(*(knot_values_[rr])));
  }

  return;
}

#endif
