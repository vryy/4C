/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for handling of parallel data exchange

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_comm_parobject.hpp"

FOUR_C_NAMESPACE_OPEN


void Core::Communication::ParObject::AddtoPack(PackBuffer& data, const ParObject& obj)
{
  obj.Pack(data);
}

void Core::Communication::ParObject::AddtoPack(PackBuffer& data, const ParObject* obj)
{
  obj->Pack(data);
}

void Core::Communication::ParObject::AddtoPack(
    PackBuffer& data, const Core::LinAlg::SerialDenseMatrix& stuff)
{
  int m = stuff.numRows();
  int n = stuff.numCols();
  AddtoPack(data, m);
  AddtoPack(data, n);
  double* A = stuff.values();
  AddtoPack(data, A, n * m * sizeof(double));
}

void Core::Communication::ParObject::AddtoPack(
    PackBuffer& data, const Core::LinAlg::SerialDenseVector& stuff)
{
  int m = stuff.length();
  AddtoPack(data, m);
  double* A = stuff.values();
  AddtoPack(data, A, m * sizeof(double));
}

void Core::Communication::ParObject::AddtoPack(PackBuffer& data, const std::string& stuff)
{
  int numele = stuff.size();
  AddtoPack(data, numele);
  AddtoPack(data, stuff.data(), numele * sizeof(char));
}

void Core::Communication::ParObject::ExtractfromPack(std::vector<char>::size_type& position,
    const std::vector<char>& data, Core::LinAlg::SerialDenseMatrix& stuff)
{
  int m = 0;
  ExtractfromPack(position, data, m);
  int n = 0;
  ExtractfromPack(position, data, n);
  stuff.reshape(m, n);
  double* a = stuff.values();
  if (m * n > 0) ExtractfromPack(position, data, a, n * m * sizeof(double));
}

void Core::Communication::ParObject::ExtractfromPack(std::vector<char>::size_type& position,
    const std::vector<char>& data, Core::LinAlg::SerialDenseVector& stuff)
{
  int m = 0;
  ExtractfromPack(position, data, m);
  stuff.resize(m);
  double* a = stuff.values();
  if (m > 0) ExtractfromPack(position, data, a, m * sizeof(double));
}

void Core::Communication::ParObject::ExtractfromPack(
    std::vector<char>::size_type& position, const std::vector<char>& data, std::string& stuff)
{
  int dim = 0;
  ExtractfromPack(position, data, dim);
  stuff.resize(dim);
  int size = dim * sizeof(char);
  ExtractfromPack(position, data, stuff.data(), size);
}

int Core::Communication::ExtractAndAssertId(std::vector<char>::size_type& position,
    const std::vector<char>& data, const int desired_type_id)
{
  int type_id = 0;
  Core::Communication::ParObject::ExtractfromPack(position, data, type_id);

  std::string error_message = "Wrong instance type data. The extracted type id is " +
                              std::to_string(type_id) + ", while the desired type id is " +
                              std::to_string(desired_type_id);
  FOUR_C_ASSERT(type_id == desired_type_id, error_message.c_str());

  return type_id;
}
FOUR_C_NAMESPACE_CLOSE
