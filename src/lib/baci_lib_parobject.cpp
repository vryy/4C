/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for handling of parallel data exchange

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_lib_parobject.H"


void DRT::ParObject::AddtoPack(PackBuffer& data, const ParObject& obj) { obj.Pack(data); }

void DRT::ParObject::AddtoPack(PackBuffer& data, const ParObject* obj) { obj->Pack(data); }

void DRT::ParObject::AddtoPack(PackBuffer& data, const CORE::LINALG::SerialDenseMatrix& stuff)
{
  int m = stuff.numRows();
  int n = stuff.numCols();
  AddtoPack(data, m);
  AddtoPack(data, n);
  double* A = stuff.values();
  AddtoPack(data, A, n * m * sizeof(double));
}

void DRT::ParObject::AddtoPack(PackBuffer& data, const CORE::LINALG::SerialDenseVector& stuff)
{
  int m = stuff.length();
  AddtoPack(data, m);
  double* A = stuff.values();
  AddtoPack(data, A, m * sizeof(double));
}

void DRT::ParObject::AddtoPack(PackBuffer& data, const std::string& stuff)
{
  int numele = stuff.size();
  AddtoPack(data, numele);
  AddtoPack(data, stuff.data(), numele * sizeof(char));
}

void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position,
    const std::vector<char>& data, CORE::LINALG::SerialDenseMatrix& stuff)
{
  int m = 0;
  ExtractfromPack(position, data, m);
  int n = 0;
  ExtractfromPack(position, data, n);
  stuff.reshape(m, n);
  double* a = stuff.values();
  if (m * n > 0) ExtractfromPack(position, data, a, n * m * sizeof(double));
}

void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position,
    const std::vector<char>& data, CORE::LINALG::SerialDenseVector& stuff)
{
  int m = 0;
  ExtractfromPack(position, data, m);
  stuff.resize(m);
  double* a = stuff.values();
  if (m > 0) ExtractfromPack(position, data, a, m * sizeof(double));
}

void DRT::ParObject::ExtractfromPack(
    std::vector<char>::size_type& position, const std::vector<char>& data, std::string& stuff)
{
  int dim = 0;
  ExtractfromPack(position, data, dim);
  stuff.resize(dim);
  int size = dim * sizeof(char);
  ExtractfromPack(position, data, stuff.data(), size);
}
