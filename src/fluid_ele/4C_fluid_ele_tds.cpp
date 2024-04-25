/*-----------------------------------------------------------*/
/*! \file

\brief time-dependent subgrid scale functionality


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_fluid_ele_tds.hpp"

FOUR_C_NAMESPACE_OPEN


FLD::TDSEleDataType FLD::TDSEleDataType::instance_;


/*----------------------------------------------------------------------*
 |  ctor (public)                                              gjb 12/12|
 *----------------------------------------------------------------------*/
FLD::TDSEleData::TDSEleData()
{
  saccn_.shape(0, 0);
  svelnp_.shape(0, 0);
  sveln_.shape(0, 0);

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
void FLD::TDSEleData::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // history variables
  AddtoPack(data, saccn_.numRows());
  AddtoPack(data, saccn_.numCols());

  int size = saccn_.numRows() * saccn_.numCols() * sizeof(double);

  AddtoPack(data, saccn_.values(), size);
  AddtoPack(data, svelnp_.values(), size);
  AddtoPack(data, sveln_.values(), size);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // history variables (subgrid-scale velocities, accelerations and pressure)
  {
    int firstdim;
    int secondim;

    ExtractfromPack(position, data, firstdim);
    ExtractfromPack(position, data, secondim);


    saccn_.shape(firstdim, secondim);
    svelnp_.shape(firstdim, secondim);
    sveln_.shape(firstdim, secondim);


    int size = firstdim * secondim * sizeof(double);

    ExtractfromPack(position, data, saccn_.values(), size);
    ExtractfromPack(position, data, svelnp_.values(), size);
    ExtractfromPack(position, data, sveln_.values(), size);
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  activate time dependent subgrid scales (public)           gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::ActivateTDS(
    int nquad, int nsd, double** saccn, double** sveln, double** svelnp)
{
  if (saccn_.numRows() != nsd || saccn_.numCols() != nquad)
  {
    saccn_.shape(nsd, nquad);
    sveln_.shape(nsd, nquad);
    svelnp_.shape(nsd, nquad);
  }
  // for what's this exactly?
  if (saccn != nullptr) *saccn = saccn_.values();
  if (sveln != nullptr) *sveln = sveln_.values();
  if (svelnp != nullptr) *svelnp = svelnp_.values();
}



/*----------------------------------------------------------------------*
 |  activate time dependent subgrid-scales (public)           gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::UpdateSvelnpInOneDirection(const double fac1, const double fac2,
    const double fac3, const double resM, const double alphaF, const int dim, const int iquad,
    double& svelaf)
{
  /*
      ~n+1           ~n           ~ n            n+1
      u    =  fac1 * u  + fac2 * acc  -fac3 * res
       (i)

  */

  svelnp_(dim, iquad) = fac1 * sveln_(dim, iquad) + fac2 * saccn_(dim, iquad) - fac3 * resM;

  /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

  */
  svelaf = alphaF * svelnp_(dim, iquad) + (1.0 - alphaF) * sveln_(dim, iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  activate time dependend subgrid scales (public)           gjb 12/12 |
 *----------------------------------------------------------------------*/
void FLD::TDSEleData::UpdateSvelnpInOneDirection(const double fac1, const double fac2,
    const double fac3, const double resM, const double alphaF, const int dim, const int iquad,
    double& svelnp, double& svelaf)
{
  UpdateSvelnpInOneDirection(fac1, fac2, fac3, resM, alphaF, dim, iquad, svelaf);

  svelnp = svelnp_(dim, iquad);

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
  const int nsd = saccn_.numRows();  // does this always hold?
  for (int rr = 0; rr < nsd; ++rr)
  {
    for (int mm = 0; mm < svelnp_.numCols(); ++mm)
    {
      saccn_(rr, mm) = (svelnp_(rr, mm) - sveln_(rr, mm)) / (gamma * dt) -
                       saccn_(rr, mm) * (1.0 - gamma) / gamma;
    }
  }

  // most recent subgrid-scale velocity becomes the old subscale velocity
  // for the next timestep
  //
  //  ~n   ~n+1
  //  u <- u
  //
  // variable in space dimensions
  for (int rr = 0; rr < nsd; ++rr)
  {
    for (int mm = 0; mm < svelnp_.numCols(); ++mm)
    {
      sveln_(rr, mm) = svelnp_(rr, mm);
    }
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
