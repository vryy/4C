
#ifdef CCADISCRET

#include <vector>
#include "newtonianfluid.H"

extern struct _MATERIAL *mat;


MAT::NewtonianFluid::NewtonianFluid()
  : matdata_(NULL)
{
}


MAT::NewtonianFluid::NewtonianFluid(MATERIAL* matdata)
  : matdata_(matdata)
{
}


void MAT::NewtonianFluid::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;
  AddtoPack(data,matdata);
}


void MAT::NewtonianFluid::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}

#endif
