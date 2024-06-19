/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for handling of parallel data exchange

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_comm_parobject.hpp"

FOUR_C_NAMESPACE_OPEN


void Core::Communication::ParObject::add_to_pack(PackBuffer& data, const ParObject& obj)
{
  obj.pack(data);
}

void Core::Communication::ParObject::add_to_pack(PackBuffer& data, const ParObject* obj)
{
  obj->pack(data);
}

void Core::Communication::ParObject::add_to_pack(
    PackBuffer& data, const Core::LinAlg::SerialDenseMatrix& stuff)
{
  int m = stuff.numRows();
  int n = stuff.numCols();
  add_to_pack(data, m);
  add_to_pack(data, n);
  double* A = stuff.values();
  add_to_pack(data, A, n * m * sizeof(double));
}

void Core::Communication::ParObject::add_to_pack(
    PackBuffer& data, const Core::LinAlg::SerialDenseVector& stuff)
{
  int m = stuff.length();
  add_to_pack(data, m);
  double* A = stuff.values();
  add_to_pack(data, A, m * sizeof(double));
}

void Core::Communication::ParObject::add_to_pack(PackBuffer& data, const std::string& stuff)
{
  int numele = stuff.size();
  add_to_pack(data, numele);
  add_to_pack(data, stuff.data(), numele * sizeof(char));
}

void Core::Communication::ParObject::extract_from_pack(std::vector<char>::size_type& position,
    const std::vector<char>& data, Core::LinAlg::SerialDenseMatrix& stuff)
{
  int m = 0;
  extract_from_pack(position, data, m);
  int n = 0;
  extract_from_pack(position, data, n);
  stuff.reshape(m, n);
  double* a = stuff.values();
  if (m * n > 0) extract_from_pack(position, data, a, n * m * sizeof(double));
}

void Core::Communication::ParObject::extract_from_pack(std::vector<char>::size_type& position,
    const std::vector<char>& data, Core::LinAlg::SerialDenseVector& stuff)
{
  int m = 0;
  extract_from_pack(position, data, m);
  stuff.resize(m);
  double* a = stuff.values();
  if (m > 0) extract_from_pack(position, data, a, m * sizeof(double));
}

void Core::Communication::ParObject::extract_from_pack(
    std::vector<char>::size_type& position, const std::vector<char>& data, std::string& stuff)
{
  int dim = 0;
  extract_from_pack(position, data, dim);
  stuff.resize(dim);
  int size = dim * sizeof(char);
  extract_from_pack(position, data, stuff.data(), size);
}

int Core::Communication::ExtractAndAssertId(std::vector<char>::size_type& position,
    const std::vector<char>& data, const int desired_type_id)
{
  int type_id = 0;
  Core::Communication::ParObject::extract_from_pack(position, data, type_id);

  std::string error_message = "Wrong instance type data. The extracted type id is " +
                              std::to_string(type_id) + ", while the desired type id is " +
                              std::to_string(desired_type_id);
  FOUR_C_ASSERT(type_id == desired_type_id, error_message.c_str());

  return type_id;
}
FOUR_C_NAMESPACE_CLOSE
