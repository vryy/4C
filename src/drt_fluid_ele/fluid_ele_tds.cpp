/*!----------------------------------------------------------------------
\file fluid_ele_tds.cpp
\brief
\level 3
<pre>
\maintainer Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#include "fluid_ele_tds.H"


FLD::TDSEleDataType FLD::TDSEleDataType::instance_;


/*----------------------------------------------------------------------*
 |  ctor (public)                                              gjb 12/12|
 *----------------------------------------------------------------------*/
FLD::TDSEleData::TDSEleData()
{
    saccn_ .Shape(0,0);
    svelnp_.Shape(0,0);
    sveln_ .Shape(0,0);

    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         gjb 12/12|
 *----------------------------------------------------------------------*/
/*
FLD::TDSEleData::TDSEleData()(const FLD::TDSEleData& old) :
saccn_      (old.saccn_      ),
svelnp_     (old.svelnp_     ),
sveln_      (old.sveln_      )
{
    return;
}
*/

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // history variables
  AddtoPack(data,saccn_.M());
  AddtoPack(data,saccn_.N());

  int size = saccn_.M()*saccn_.N()*sizeof(double);

  AddtoPack(data,saccn_ .A(),size);
  AddtoPack(data,svelnp_.A(),size);
  AddtoPack(data,sveln_ .A(),size);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");

  // history variables (subgrid-scale velocities, accelerations and pressure)
  {
    int firstdim;
    int secondim;

    ExtractfromPack(position,data,firstdim);
    ExtractfromPack(position,data,secondim);


    saccn_ .Shape(firstdim,secondim);
    svelnp_.Shape(firstdim,secondim);
    sveln_ .Shape(firstdim,secondim);


    int size = firstdim*secondim*sizeof(double);

    ExtractfromPack(position,data,&(saccn_ .A()[0]),size);
    ExtractfromPack(position,data,&(svelnp_.A()[0]),size);
    ExtractfromPack(position,data,&(sveln_ .A()[0]),size);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 12/12 |
 *----------------------------------------------------------------------*/
FLD::TDSEleData::~TDSEleData()
{
  return;
}

/*----------------------------------------------------------------------*
 |  activate time dependent subgrid scales (public)           gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::ActivateTDS(
    int nquad,
    int nsd,
    double** saccn,
    double** sveln,
    double** svelnp
 )
{
  if(saccn_.M() != nsd
      ||
      saccn_.N() != nquad)
  {
    saccn_ .Shape(nsd,nquad);
    memset(saccn_.A() ,0,nsd*nquad*sizeof(double));

    sveln_ .Shape(nsd,nquad);
    memset(sveln_.A() ,0,nsd*nquad*sizeof(double));

    svelnp_.Shape(nsd,nquad);
    memset(svelnp_.A(),0,nsd*nquad*sizeof(double));
  }
  // for what's this exactly?
  if ( saccn !=NULL ) *saccn = saccn_.A();
  if ( sveln !=NULL ) *sveln = sveln_.A();
  if ( svelnp!=NULL ) *svelnp = svelnp_.A();
}



/*----------------------------------------------------------------------*
 |  activate time dependent subgrid-scales (public)           gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::UpdateSvelnpInOneDirection(
    const double  fac1,
    const double  fac2,
    const double  fac3,
    const double  resM ,
    const double  alphaF,
    const int     dim,
    const int     iquad,
    double&       svelaf
    )
{

    /*
        ~n+1           ~n           ~ n            n+1
        u    =  fac1 * u  + fac2 * acc  -fac3 * res
         (i)

    */

  svelnp_(dim,iquad)=
    fac1*sveln_(dim,iquad)
    +
    fac2*saccn_(dim,iquad)
    -
    fac3*resM;

  /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

  */
  svelaf=
    alphaF      *svelnp_(dim,iquad)
    +
    (1.0-alphaF)*sveln_ (dim,iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  activate time dependend subgrid scales (public)           gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::UpdateSvelnpInOneDirection(
    const double  fac1,
    const double  fac2,
    const double  fac3,
    const double  resM,
    const double  alphaF,
    const int     dim,
    const int     iquad,
    double&       svelnp,
    double&       svelaf
    )
    {
      UpdateSvelnpInOneDirection(fac1,
                                 fac2,
                                 fac3,
                                 resM,
                                 alphaF,
                                 dim ,
                                 iquad,
                                 svelaf);

      svelnp=svelnp_(dim,iquad);

      return;
    }

/*----------------------------------------------------------------------*
 |  update time dependent subgrid scales (public)             gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::Update(const double dt, const double gamma)
{

// the old subgrid-scale acceleration for the next timestep is calculated
// on the fly, not stored on the element
/*
                 ~n+1   ~n
         ~ n+1     u    - u     ~ n   / 1.0-gamma \
        acc    =   --------- - acc * |  ---------  |
                   gamma*dt           \   gamma   /

         ~ n       ~ n+1   / 1.0-gamma \
        acc    =    acc * |  ---------  |
 */

// variable in space dimensions
  const int nsd = saccn_.M(); // does this always hold?
for(int rr=0;rr<nsd;++rr)
{
  for(int mm=0;mm<svelnp_.N();++mm)
  {
    saccn_(rr,mm) =
        (svelnp_(rr,mm)-sveln_(rr,mm))/(gamma*dt)
        -
        saccn_(rr,mm)*(1.0-gamma)/gamma;
  }
}

// most recent subgrid-scale velocity becomes the old subscale velocity
// for the next timestep
//
//  ~n   ~n+1
//  u <- u
//
// variable in space dimensions
for(int rr=0;rr<nsd;++rr)
{
  for(int mm=0;mm<svelnp_.N();++mm)
  {
    sveln_(rr,mm)=svelnp_(rr,mm);
  }
}

return;
}

