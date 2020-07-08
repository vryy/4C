/*---------------------------------------------------------------------*/
/*! \file

\brief knot vectors for nurbs problems (isogeometric analysis)

\level 1

*/
/*---------------------------------------------------------------------*/

#include "drt_knotvector.H"


DRT::NURBS::KnotvectorObjectType DRT::NURBS::KnotvectorObjectType::instance_;


DRT::ParObject* DRT::NURBS::KnotvectorObjectType::Create(const std::vector<char>& data)
{
  return NULL;
}


/*----------------------------------------------------------------------*
 |  empty ctor (public)                                      gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::Knotvector()
    : ParObject(),
      dim_(0),
      npatches_(0),
      filled_(false),
      degree_(0),
      n_x_m_x_l_(0),
      nele_x_mele_x_lele_(0),
      interpolation_(0),
      offsets_(0),
      knot_values_(0)
{
  return;
}  // DRT::NURBS::Knotvector::Knotvector (empty)

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::Knotvector(const int dim, const int npatches)
    : ParObject(),
      dim_(dim),
      npatches_(npatches),
      filled_(false),
      degree_(npatches),
      n_x_m_x_l_(npatches),
      nele_x_mele_x_lele_(npatches),
      interpolation_(npatches),
      offsets_(npatches),
      knot_values_(npatches)
{
  // check if there are any patches
  if (npatches_ < 0)
  {
    dserror("A knot vector needs at least one patch.");
  }

  // resize degrees

  // loop over all patches and resize to dimension
  for (int rr = 0; rr < npatches_; ++rr)
  {
    (degree_[rr]).resize(dim_);
    (n_x_m_x_l_[rr]).resize(dim_);
    (nele_x_mele_x_lele_[rr]).resize(dim_);
    (interpolation_[rr]).resize(dim_);
    (knot_values_[rr]).resize(dim_);

    // initialize interpolation type
    for (int mm = 0; mm < dim_; ++mm) (interpolation_[rr])[mm] = knotvector_is_not_defined;
  }

  return;
}  // DRT::NURBS::Knotvector::Knotvector (standard)

/*----------------------------------------------------------------------*
 |  copy ctor (public)                                       gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::Knotvector::Knotvector(const DRT::NURBS::Knotvector& old)
    : ParObject(old),
      dim_(old.dim_),
      npatches_(old.npatches_),
      filled_(false),
      degree_(old.degree_),
      n_x_m_x_l_(old.n_x_m_x_l_),
      nele_x_mele_x_lele_(old.nele_x_mele_x_lele_),
      interpolation_(old.interpolation_),
      offsets_(old.offsets_),
      knot_values_(old.npatches_)
{
  // deep copy knot vectors

  for (int np = 0; np < npatches_; ++np)
  {
    (knot_values_[np]).resize(dim_);
    for (int rr = 0; rr < dim_; ++rr)
    {
      ((knot_values_[np])[rr]) = Teuchos::rcp(new std::vector<double>);
      *((knot_values_[np])[rr]) = *((old.knot_values_[np])[rr]);
    }
  }
  return;
}  // DRT::NURBS::Knotvector::Knotvector (copy)


/*----------------------------------------------------------------------*
 | convert an element gid to its corresponding triple knot index        |
 |                                                  (public) gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::ConvertEleGidToKnotIds(
    const int gid, int& npatch, std::vector<int>& loc_cart_id) const
{
  if ((int)loc_cart_id.size() != dim_)
  {
    dserror("size vector not of appropriate size (%d,%d)\n", (int)loc_cart_id.size(), dim_);
  }

  if (filled_ == false)
  {
    dserror("cannot convert ele ids when filled is false\n");
  }

  // get number of patch containing the node
  npatch = ReturnPatchId(gid);

  // reduce gid by patchoffset to get patch local id
  const int locid = gid - offsets_[npatch];

  if (dim_ == 3)
  {
    // locid = num_u+num_v*nele+num_w*nele*mele     (3d)
    //         |              |       |       |
    //         +--------------+       +-------+
    //            inthislayer          uv_layer
    int uv_layer = (nele_x_mele_x_lele_[npatch])[0] * (nele_x_mele_x_lele_[npatch])[1];

    // compute num_w
    loc_cart_id[2] = locid / uv_layer;

    // see above
    int inthislayer = locid % uv_layer;

    // compute num_v and num_u
    loc_cart_id[1] = inthislayer / (nele_x_mele_x_lele_[npatch])[0];
    loc_cart_id[0] = inthislayer % (nele_x_mele_x_lele_[npatch])[0];
  }
  else if (dim_ == 2)
  {
    // locid = num_u+num_v*nele                     (2d)

    // compute num_v and num_u
    loc_cart_id[0] = locid % (nele_x_mele_x_lele_[npatch])[0];
    loc_cart_id[1] = locid / (nele_x_mele_x_lele_[npatch])[0];
  }
  else
  {
    dserror("dim_ not available\n");
  }

  return;
}  // ConvertEleGidToKnotIds


/*----------------------------------------------------------------------*
 | get element knot vectors to a given element id   (public) gammi 05/08|
 *----------------------------------------------------------------------*/
bool DRT::NURBS::Knotvector::GetEleKnots(
    std::vector<Epetra_SerialDenseVector>& eleknots, int gid) const
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


  if (filled_ == false)
  {
    dserror("cannot get ele knots when filled is false\n");
  }

  eleknots.resize(dim_);

  // get base indices and patch number from element gid
  std::vector<int> cartids(dim_);
  int npatch;

  // ToDo The current NURBS implementation relies totally on a set of consecutively
  // increasing GIDs starting at 0. This should be fixed if mixed discretizations
  // shall be considered. Otherwise the elements belonging to the NURBS
  // discretization must be a closed set placed always at the beginning
  // of an element list which seems to be a unnecessary and non-intuitive restriction.
  // hiermeier 11/17

  //  gid = gid - 3773;

  ConvertEleGidToKnotIds(gid, npatch, cartids);

  // use them to aquire the required knots
  for (int rr = 0; rr < dim_; ++rr)
  {
    (eleknots[rr]).Size(2 * (degree_[npatch])[rr] + 2);

    for (int mm = 0; mm < 2 * (degree_[npatch])[rr] + 2; ++mm)
    {
      (eleknots[rr])(mm) = (*((knot_values_[npatch])[rr]))[cartids[rr] + mm];
    }
  }

  // check local knotspan for multiple knots indicating zero
  // sized elements
  bool zero_size = false;

  for (int rr = 0; rr < dim_; ++rr)
  {
    const double size =
        (eleknots[rr])((degree_[npatch])[rr]) - (eleknots[rr])((degree_[npatch])[rr] + 1);

    if (fabs(size) < 1e-12)
    {
      zero_size = true;
    }
  }

  return (zero_size);
}  // DRT::NURBS::Knotvector::GetEleKnots


/*----------------------------------------------------------------------*
 | extract element surface's knot vectors out of the knot vector of the |
 | parent element. On the fly, get orientation of normal vector.        |
 |                                                  (public) gammi 05/09|
 *----------------------------------------------------------------------*/
bool DRT::NURBS::Knotvector::GetBoundaryEleAndParentKnots(
    std::vector<Epetra_SerialDenseVector>& eleknots,
    std::vector<Epetra_SerialDenseVector>& surfknots, double& normalfac, int pgid,
    const int surfaceid) const
{
  // get parent element local knotspan to extract the surface's knotspan
  // from
  //
  // immediately, check for multiple knots indicating zero sized elements
  bool zero_size = GetEleKnots(eleknots, pgid);

  // locate patch containing parent element
  const int np = ReturnPatchId(pgid);

  if (dim_ == 3)
  {
    // check for nurbs27 elements to get numbering
    if ((degree_[np])[0] == 2 && (degree_[np])[1] == 2 && (degree_[np])[2] == 2)
    {
      switch (surfaceid)
      {
        case 0:
        {
          // t=-1
          /*
                  parent               surface

                   s|                    s|
                    |                     |
                +---+---+             +---+---+
               6|  7|  8|      r     6|  7|  8|      r
                +   +-- +  -----      +   +-- +  -----
               3|  4   5|            3|  4   5|
                +---+---+             +---+---+
               0   1   2             0   1   2
          */
          surfknots[0].Size(eleknots[0].Length());
          surfknots[1].Size(eleknots[1].Length());
          surfknots[0] = eleknots[0];
          surfknots[1] = eleknots[1];

          normalfac = -1.0;
          break;
        }
        case 1:
        {
          // t=+1
          /*
                  parent               surface

                   s|                    s|
                    |                     |
                +---+---+             +---+---+
              24| 25| 26|      r     6|  7|  8|      r
                +   +-- +  -----      +   +-- +  -----
              21| 22  23|            3|  4   5|
                +---+---+             +---+---+
              18  19  20             0   1   2
          */
          surfknots[0].Size(eleknots[0].Length());
          surfknots[1].Size(eleknots[1].Length());
          surfknots[0] = eleknots[0];
          surfknots[1] = eleknots[1];

          normalfac = 1.0;
          break;
        }
        case 2:
        {
          // s=-1
          /*
                  parent               surface

                   t|                    s|
                    |                     |
                +---+---+             +---+---+
              18| 19| 20|      r     6|  7|  8|      r
                +   +-- +  -----      +   +-- +  -----
               9| 10  11|            3|  4   5|
                +---+---+             +---+---+
               0   1   2             0   1   2
          */
          surfknots[0].Size(eleknots[0].Length());
          surfknots[1].Size(eleknots[2].Length());
          surfknots[0] = eleknots[0];
          surfknots[1] = eleknots[2];

          normalfac = 1.0;
          break;
        }
        case 3:
        {
          // s=+1
          /*
                  parent               surface

                   t|                    s|
                    |                     |
                +---+---+             +---+---+
              24| 25| 26|    r       6|  7|  8|      r
                +   +-- + ----        +   +-- +  -----
              15| 16  17|            3|  4   5|
                +---+---+             +---+---+
               6   7   8             0   1   2
          */
          surfknots[0].Size(eleknots[0].Length());
          surfknots[1].Size(eleknots[2].Length());
          surfknots[0] = eleknots[0];
          surfknots[1] = eleknots[2];

          normalfac = -1.0;
          break;
        }
        case 4:
        {
          // r=+1
          /*
                  parent               surface

                   t|                    s|
                    |                     |
                +---+---+             +---+---+
              20| 23| 26|      s     6|  7|  8|      r
                +   +-- +  -----      +   +-- +  -----
              11| 14  17|            3|  4   5|
                +---+---+             +---+---+
               2   5   8             0   1   2
          */
          surfknots[0].Size(eleknots[1].Length());
          surfknots[1].Size(eleknots[2].Length());
          surfknots[0] = eleknots[1];
          surfknots[1] = eleknots[2];

          normalfac = 1.0;
          break;
        }
        case 5:
        {
          // r=-1
          /*
                  parent               surface

                   t|                    s|
                    |                     |
                +---+---+             +---+---+
              18| 21| 24|      s     6|  7|  8|      r
                +   +-- +  -----      +   +-- +  -----
               9| 12  15|            3|  4   5|
                +---+---+             +---+---+
               0   3   6             0   1   2
          */
          surfknots[0].Size(eleknots[1].Length());
          surfknots[1].Size(eleknots[2].Length());
          surfknots[0] = eleknots[1];
          surfknots[1] = eleknots[2];

          normalfac = -1.0;
          break;
        }
        default:
          dserror("invalid number of surfaces, unable to determine intpoint in parent");
      }
    }  // if(degree_[0]==2 && degree_[1]==2 && degree_[2]==2)
    else
    {
      dserror("surface knot vector extraction not available for this degree\n");
    }
  }  // if(dim_==3)
  else if (dim_ == 2)
  {
    // check for nurbs9 elements to apply numbering
    if ((degree_[np])[0] == 2 && (degree_[np])[1] == 2)
    {
      switch (surfaceid)
      {
        case 0:
        {
          // s=-1
          /*

                  parent                line

                               r                     r
                +---+---+  -----      +---+---+ ------
               0   1   2             0   1   2

          */
          surfknots[0].Size(eleknots[0].Length());
          surfknots[0] = eleknots[0];

          normalfac = 1.0;
          break;
        }
        case 1:
        {
          // r=+1
          /*
                  parent               surface

                   s|                        r|
                    |                         |
                    +                         +
                   8|                        2|
                    +                         +
                   5|                        1|
                    +                         +
                    2                         0
          */
          surfknots[0].Size(eleknots[1].Length());
          surfknots[0] = eleknots[1];

          normalfac = 1.0;
          break;
        }
        case 2:
        {
          // s=+1
          /*

                  parent                line

                               r                           r
                +---+---+  -----             +---+---+ -----
               6   7   8                    0   1   2

          */
          surfknots[0].Size(eleknots[0].Length());
          surfknots[0] = eleknots[0];

          normalfac = -1.0;
          break;
        }
        case 3:
        {
          // r=-1
          /*
                  parent               surface

                   s|                        r|
                    |                         |
                    +                         +
                   6|                        2|
                    +                         +
                   3|                        1|
                    +                         +
                    0                         0
          */

          surfknots[0].Size(eleknots[1].Length());
          surfknots[0] = eleknots[1];

          normalfac = -1.0;
          break;
        }
        default:
          dserror("invalid number of surfaces, unable to determine intpoint in parent");
      }
    }  // if(degree_[0]==2 && degree_[1]==2)
    else
    {
      dserror("line knot vector extraction not available for this degree\n");
    }
  }  // if(dim_==2)
  else
  {
    dserror("surface knot vector extraction only in 2d and 3d\n");
  }

  return (zero_size);
}  // DRT::NURBS::Knotvector::GetBoundaryEleAndParentKnots


/*----------------------------------------------------------------------*
 | set knots in one direction                       (public) gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::SetKnots(const int& direction, const int& npatch, const int& degree,
    const int& numknots, const std::string& knotvectortype,
    Teuchos::RCP<std::vector<double>> directions_knots)
{
  // filled is false now since new add new knots
  filled_ = false;

  if (direction < 0 || direction > dim_ - 1)
  {
    dserror("direction has to in[0...dim_]\n");
  }

  if (npatch < 0 || npatch > npatches_ - 1)
  {
    dserror("patchnumber is invalid\n");
  }

  // set the type
  if (knotvectortype == "Interpolated")
  {
    (interpolation_[npatch])[direction] = knotvector_is_interpolating;
  }
  else if (knotvectortype == "Periodic")
  {
    (interpolation_[npatch])[direction] = knotvector_is_periodic;
  }
  else
  {
    dserror("unknown knotvector-type '%s'\n", knotvectortype.c_str());
  }

  // set the degree of the added knotvector
  (degree_[npatch])[direction] = degree;

  // set the size of the added knotvector
  (n_x_m_x_l_[npatch])[direction] = numknots;

  // set the actual values
  (knot_values_[npatch])[direction] = directions_knots;

  return;
}  // DRT::NURBS::Knotvector::SetKnots

/*----------------------------------------------------------------------*
 | finish                                           (public) gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::FinishKnots(const int smallest_gid_in_dis)
{
  //--------------------------------------------------
  // plausibility checks

  // empty knot vector --- do not finish or set
  if (npatches_ == 0)
  {
    return;
  }

  // check if there are any patches
  if (npatches_ < 0)
  {
    dserror("There are no patches. We need at least one patch.");
  }

  // check degrees
  if ((int)degree_.size() != npatches_)
  {
    dserror("Each patch needs its own degree information.");
  }
  else
  {
    // loop patches, check dimensions
    for (int rr = 0; rr < npatches_; ++rr)
    {
      if ((int)(degree_[rr]).size() != dim_)
      {
        dserror("size mismatch: degree\n");
      }
    }
  }

  // check n_x_m_x_ls
  if ((int)n_x_m_x_l_.size() != npatches_)
  {
    dserror("each patch needs its own n_x_m_x_l information\n");
  }
  else
  {
    // loop patches, check dimensions
    for (int rr = 0; rr < npatches_; ++rr)
    {
      if ((int)(n_x_m_x_l_[rr]).size() != dim_)
      {
        dserror("size mismatch: n_x_m_x_l\n");
      }
    }
  }

  // do we have a knotvector for each dimension
  // and each patch?
  if ((int)knot_values_.size() != npatches_)
  {
    dserror("each patch needs its own knotvector\n");
  }

  // loop patches
  for (int np = 0; np < npatches_; ++np)
  {
    if ((int)(knot_values_[np]).size() != dim_)
    {
      dserror("knotvector of patch has to be of size dim\n");
    }

    for (int rr = 0; rr < dim_; ++rr)
    {
      // is the knotvector of this dimension nonempty?
      if ((knot_values_[np])[rr] == Teuchos::null)
      {
        dserror("no knotvector available in this direction\n");
      }

      // has it the correct size?
      if ((int)(*((knot_values_[np])[rr])).size() != (n_x_m_x_l_[np])[rr])
      {
        dserror("knotvector size mismatch to n_x_m_x_l_ %d!=%d\n",
            (*((knot_values_[np])[rr])).size(), (n_x_m_x_l_[np])[rr]);
      }

      // is interpolation/periodicity assigned correctly?
      if ((interpolation_[np])[rr] == knotvector_is_not_defined)
      {
        dserror("undefined knotvector type\n");
      }
      else if ((interpolation_[np])[rr] == knotvector_is_interpolating)
      {
        // for interpolating knot vectors, the first and last
        // knots have to be repeated degree+1 times
        double firstval = (*((knot_values_[np])[rr]))[0];
        double lastval = (*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr] - 1];

        for (int mm = 1; mm < (degree_[np])[rr] + 1; ++mm)
        {
          double db = abs((*((knot_values_[np])[rr]))[mm] - firstval);
          double de = abs((*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr] - 1 - mm] - lastval);

          if (de > 1e-9 || db > 1e-9)
          {
            dserror("need multiple knots at the beginning and end of an interpolated knotvector\n");
          }
        }
      }
      else if ((interpolation_[np])[rr] == knotvector_is_periodic)
      {
        // for periodic knot vectors, distances between the
        // degree+1 first and last nodes have to be equal
        for (int mm = 1; mm < (degree_[np])[rr] + 1; ++mm)
        {
          double db = (*((knot_values_[np])[rr]))[mm] - (*((knot_values_[np])[rr]))[mm - 1];
          double de = (*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr] - mm] -
                      (*((knot_values_[np])[rr]))[(n_x_m_x_l_[np])[rr] - 1 - mm];
          if (abs(de - db) > 1e-9)
          {
            dserror("periodic knotvector doesn't obey periodicity\n");
          }
        }
      }
    }  // loop dimensions
  }    // end loop patches

  //--------------------------------------------------
  // generate offset arrays for element to patch
  // mapping and size of element arrays of patches

  // get patches element distribution
  for (int rr = 0; rr < npatches_; ++rr)
  {
    for (int mm = 0; mm < dim_; ++mm)
    {
      (nele_x_mele_x_lele_[rr])[mm] = (n_x_m_x_l_[rr])[mm] - 2 * (degree_[rr])[mm] - 1;
    }
  }

  // get element ordering among patches
  offsets_[0] = smallest_gid_in_dis;
  for (int rr = 1; rr < npatches_; ++rr)
  {
    int nele_inpatch = 1;
    for (int mm = 0; mm < dim_; ++mm)
    {
      nele_inpatch *= (nele_x_mele_x_lele_[rr - 1])[mm];
    }
    offsets_[rr] = offsets_[rr - 1] + nele_inpatch;
  }

  //--------------------------------------------------
  // the knotvector is OK
  filled_ = true;

  return;
}  // DRT::NURBS::Knotvector::FinishKnots

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 05/08 |
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::Pack(DRT::PackBuffer& data) const
{
  // we don't need the PackBuffer for the knotvector (at the moment)
  // DRT::PackBuffer::SizeMarker sm( data );
  // sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add number of patches
  AddtoPack(data, npatches_);

  // add dimension
  AddtoPack(data, dim_);

  // add degree vector
  for (int np = 0; np < npatches_; ++np)
  {
    AddtoPack(data, degree_[np]);
  }

  // add knotvector size
  for (int np = 0; np < npatches_; ++np)
  {
    AddtoPack(data, n_x_m_x_l_[np]);
  }

  // add element numbers in all cartesian
  // directions
  for (int np = 0; np < npatches_; ++np)
  {
    AddtoPack(data, nele_x_mele_x_lele_[np]);
  }

  // add Knotvector types
  for (int np = 0; np < npatches_; ++np)
  {
    for (int rr = 0; rr < dim_; ++rr)
    {
      AddtoPack(data, (interpolation_[np])[rr]);
    }
  }

  // add patch offsets
  AddtoPack(data, offsets_);

  // add Knotvector coordinates itself
  for (int np = 0; np < npatches_; ++np)
  {
    for (int rr = 0; rr < dim_; ++rr)
    {
      AddtoPack(data, (*((knot_values_[np])[rr])));
    }
  }

  return;
}  // DRT::NURBS::Knotvector::Pack

/*----------------------------------------------------------------------*
 |  Unpack Knotvectors data                                    (public) |
 |                                                          gammi 05/08 |
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  filled_ = false;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract number of patches
  ExtractfromPack(position, data, npatches_);

  // extract dimension
  ExtractfromPack(position, data, dim_);

  // resize all vectors
  degree_.resize(npatches_);
  n_x_m_x_l_.resize(npatches_);
  nele_x_mele_x_lele_.resize(npatches_);
  interpolation_.resize(npatches_);
  knot_values_.resize(npatches_);

  for (int np = 0; np < npatches_; ++np)
  {
    (degree_[np]).resize(dim_);
    (n_x_m_x_l_[np]).resize(dim_);
    (nele_x_mele_x_lele_[np]).resize(dim_);
    (interpolation_[np]).resize(dim_);
    (knot_values_[np]).resize(dim_);
  }

  // extract degree vector
  for (int np = 0; np < npatches_; ++np)
  {
    ExtractfromPack(position, data, degree_[np]);
  }

  // extract knotvector size
  for (int np = 0; np < npatches_; ++np)
  {
    ExtractfromPack(position, data, n_x_m_x_l_[np]);
  }

  // extract element numbers in all cartesian
  // directions
  for (int np = 0; np < npatches_; ++np)
  {
    ExtractfromPack(position, data, nele_x_mele_x_lele_[np]);
  }

  // extract knotvector types
  for (int np = 0; np < npatches_; ++np)
  {
    for (int rr = 0; rr < dim_; ++rr)
    {
      (interpolation_[np])[rr] =
          static_cast<DRT::NURBS::Knotvector::KnotvectorType>(ExtractInt(position, data));
    }
  }

  // extract patch offsets
  ExtractfromPack(position, data, offsets_);

  // extract knotvector coordinates itself
  for (int np = 0; np < npatches_; ++np)
  {
    for (int rr = 0; rr < dim_; ++rr)
    {
      (knot_values_[np])[rr] = Teuchos::rcp(new std::vector<double>((n_x_m_x_l_[np])[rr]));

      ExtractfromPack(position, data, (*((knot_values_[np])[rr])));
    }
  }

  return;
}  // DRT::NURBS::Knotvector::Unpack


/*----------------------------------------------------------------------*
 |  Return number of zero sized elements in knotspan of this patch      |
 |  (public)                                                gammi 05/08 |
 *----------------------------------------------------------------------*/
std::vector<int> DRT::NURBS::Knotvector::Return_n_zerosize_ele(const int npatch)
{
  if (!filled_)
  {
    dserror("can't access data. knotvector not completed\n");
  }

  std::vector<int> num_zero_sized(dim_);

  switch (dim_)
  {
    case 1:
    {
      for (int rr = (degree_[npatch])[0]; rr < (n_x_m_x_l_[npatch])[0] - (degree_[npatch])[0] - 1;
           ++rr)
      {
        double size = 1.0;

        size *= (*(knot_values_[npatch])[0])[rr + 1] - (*(knot_values_[npatch])[0])[rr];

        if (fabs(size) < 1e-12)
        {
          ++(num_zero_sized[0]);
          break;
        }
      }
      break;
    }
    case 2:
    {
      double size = 0.0;
      for (int rr = (degree_[npatch])[0]; rr < (n_x_m_x_l_[npatch])[0] - (degree_[npatch])[0] - 1;
           ++rr)
      {
        size = (*(knot_values_[npatch])[0])[rr + 1] - (*(knot_values_[npatch])[0])[rr];
        if (fabs(size) < 1e-12)
        {
          ++(num_zero_sized[0]);
        }
      }
      for (int mm = (degree_[npatch])[1]; mm < (n_x_m_x_l_[npatch])[1] - (degree_[npatch])[1] - 1;
           ++mm)
      {
        size = (*(knot_values_[npatch])[1])[mm + 1] - (*(knot_values_[npatch])[1])[mm];

        if (fabs(size) < 1e-12)
        {
          ++(num_zero_sized[1]);
          break;
        }
      }
      break;
    }
    case 3:
    {
      double size = 0.0;
      for (int rr = (degree_[npatch])[0]; rr < (n_x_m_x_l_[npatch])[0] - (degree_[npatch])[0] - 1;
           ++rr)
      {
        size = (*(knot_values_[npatch])[0])[rr + 1] - (*(knot_values_[npatch])[0])[rr];
        if (fabs(size) < 1e-12)
        {
          ++(num_zero_sized[0]);
        }
      }
      for (int mm = (degree_[npatch])[1]; mm < (n_x_m_x_l_[npatch])[1] - (degree_[npatch])[1] - 1;
           ++mm)
      {
        size = (*(knot_values_[npatch])[1])[mm + 1] - (*(knot_values_[npatch])[1])[mm];

        if (fabs(size) < 1e-12)
        {
          ++(num_zero_sized[1]);
          break;
        }
      }
      for (int kk = (degree_[npatch])[2]; kk < (n_x_m_x_l_[npatch])[2] - (degree_[npatch])[2] - 1;
           ++kk)
      {
        size = (*(knot_values_[npatch])[2])[kk + 1] - (*(knot_values_[npatch])[2])[kk];

        if (fabs(size) < 1e-12)
        {
          ++(num_zero_sized[2]);
          break;
        }
      }
      break;
    }
    default:
      dserror("implemented only for 1,2 and 3 dimensions\n");
  }

  return (num_zero_sized);
}  // DRT::NURBS::Knotvector::Return_n_zerosize_ele(const int npatch)



/*----------------------------------------------------------------------*
 |  Return the global id of the next nonzero sized element in the       |
 |  knotspan                  (public)                      gammi 04/09 |
 *----------------------------------------------------------------------*/
int DRT::NURBS::Knotvector::Return_next_nonzero_ele_gid(const int zero_ele_gid)
{
  std::vector<int> zero_ele_cart_id(dim_);
  std::vector<int> nonzero_ele_cart_id(dim_);

  int npatch = -1;

  ConvertEleGidToKnotIds(zero_ele_gid, npatch, zero_ele_cart_id);

  std::vector<int> count(dim_, -1);

  for (int dir = 0; dir < dim_; ++dir)
  {
    int location;
    location = zero_ele_cart_id[dir] + degree_[npatch][dir];

    double size = 0.0;
    while (fabs(size) < 1e-12)
    {
      size =
          (*(knot_values_[npatch])[dir])[location + 1] - (*(knot_values_[npatch])[dir])[location];
      ++location;
      ++count[dir];
    }
  }

  for (int dir = 0; dir < dim_; ++dir)
  {
    nonzero_ele_cart_id[dir] = zero_ele_cart_id[dir] + count[dir];
  }
  int nextnonzero_gid = ConvertEleKnotIdsToGid(npatch, nonzero_ele_cart_id);

  return (nextnonzero_gid);
}  // DRT::NURBS::Knotvector::Return_next_nonzero_ele_gid

/*----------------------------------------------------------------------*
 | convert an element local id + patch number to its corresponding gid  |
 |                                                  (public) gammi 04/09|
 *----------------------------------------------------------------------*/
int DRT::NURBS::Knotvector::ConvertEleKnotIdsToGid(
    const int& npatch, const std::vector<int>& loc_cart_id)
{
  if (!filled_)
  {
    dserror("can't access data. knotvector not completed\n");
  }

  int gid = -1;

  switch (dim_)
  {
    case 1:
    {
      //   gid = num_u+num_v*nele                     (1d)
      //         |              |
      //         +--------------+
      //            inthislayer
      gid = loc_cart_id[0];
      gid += offsets_[npatch];

      break;
    }
    case 2:
    {
      //   gid = num_u+num_v*nele                     (2d)
      //         |              |
      //         +--------------+
      //            inthislayer
      gid = loc_cart_id[1] * nele_x_mele_x_lele_[npatch][0];
      gid += loc_cart_id[0];
      gid += offsets_[npatch];

      break;
    }
    case 3:
    {
      //   gid = num_u+num_v*nele+num_w*nele*mele     (3d)
      //         |              |       |       |
      //         +--------------+       +-------+
      //            inthislayer          uv_layer

      gid = loc_cart_id[2] * nele_x_mele_x_lele_[npatch][0] * nele_x_mele_x_lele_[npatch][1];
      gid += loc_cart_id[1] * nele_x_mele_x_lele_[npatch][0];
      gid += loc_cart_id[0];
      gid += offsets_[npatch];

      break;
    }
    default:
      dserror("implemented only for 1,2 and 3 dimensions\n");
  }

  return (gid);
}  // DRT::NURBS::Knotvector::ConvertEleKnotIdsToGid


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::NURBS::Knotvector::Print(std::ostream& os) const
{
  os << "\nPrinting a Knotvector: " << std::endl;
  for (int patch = 0; patch < npatches_; ++patch)
  {
    os << "patch " << patch << ":\n";
    for (int direction = 0; direction < dim_; ++direction)
    {
      os << "  direction " << direction << ": ";
      for (std::size_t i = 0; i < knot_values_[patch][direction]->size(); ++i)
        os << (*(knot_values_[patch][direction]))[i] << ", ";
      os << "\n";
    }
    os << std::endl;
  }
}
