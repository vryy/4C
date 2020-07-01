/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for handling of parallel data exchange

\level 0


*/
/*---------------------------------------------------------------------*/

#include "drt_parobject.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ParObject::ParObject() { return; }

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ParObject::ParObject(const DRT::ParObject& old) { return; }

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ParObject::~ParObject() { return; }

void DRT::ParObject::AddtoPack(PackBuffer& data, const ParObject& obj) { obj.Pack(data); }

void DRT::ParObject::AddtoPack(PackBuffer& data, const ParObject* obj) { obj->Pack(data); }

/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseMatrix specialization                   (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const Epetra_SerialDenseMatrix& stuff)
{
  int m = stuff.M();
  int n = stuff.N();
  AddtoPack(data, m);
  AddtoPack(data, n);
  double* A = stuff.A();
  AddtoPack(data, A, n * m * sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseVector specialization                   (public) |
 |                                                     TK & MAF  05/08  |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const Epetra_SerialDenseVector& stuff)
{
  int m = stuff.Length();
  AddtoPack(data, m);
  double* A = stuff.Values();
  AddtoPack(data, A, m * sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a LINALG::SerialDenseMatrix specialization                  (public) |
 |                                                          henke 12/09 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const LINALG::SerialDenseMatrix& stuff)
{
  int m = stuff.M();
  int n = stuff.N();
  AddtoPack(data, m);
  AddtoPack(data, n);
  double* A = stuff.A();
  AddtoPack(data, A, n * m * sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a string specialization                                     (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::string& stuff)
{
  int numele = stuff.size();
  AddtoPack(data, numele);
  AddtoPack(data, &stuff[0], numele * sizeof(char));
  return;
}

/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseMatrix specialization                   (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position,
    const std::vector<char>& data, Epetra_SerialDenseMatrix& stuff)
{
  int m = 0;
  ExtractfromPack(position, data, m);
  int n = 0;
  ExtractfromPack(position, data, n);
  stuff.Reshape(m, n);
  double* a = stuff.A();
  if (m * n > 0) ExtractfromPack(position, data, a, n * m * sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseVector specialization                   (public) |
 |                                                     TK & MAF  05/08  |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position,
    const std::vector<char>& data, Epetra_SerialDenseVector& stuff)
{
  int m = 0;
  ExtractfromPack(position, data, m);
  stuff.Resize(m);
  double* a = stuff.Values();
  if (m > 0) ExtractfromPack(position, data, a, m * sizeof(double));
  return;
}

/*----------------------------------------------------------------------*
 | a LINALG::SerialDenseMatrix specialization                  (public) |
 |                                                          henke 12/09 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position,
    const std::vector<char>& data, LINALG::SerialDenseMatrix& stuff)
{
  int m = 0;
  ExtractfromPack(position, data, m);
  int n = 0;
  ExtractfromPack(position, data, n);
  stuff.Reshape(m, n);
  double* a = stuff.A();
  if (m * n > 0) ExtractfromPack(position, data, a, n * m * sizeof(double));
  return;
}

/*----------------------------------------------------------------------*
 | a string specialization                                     (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(
    std::vector<char>::size_type& position, const std::vector<char>& data, std::string& stuff)
{
  int dim = 0;
  ExtractfromPack(position, data, dim);
  stuff.resize(dim);
  int size = dim * sizeof(char);
  ExtractfromPack(position, data, &stuff[0], size);
  return;
}
