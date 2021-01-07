/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of free service functions for elements

\level 3
*/
/*----------------------------------------------------------------------*/

#include "so_element_service.H"
#include <Epetra_SerialDenseMatrix.h>
#include "so_hex8.H"

template <class T>
void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(Epetra_MultiVector& global_data,
    const T& nodal_data, const DRT::ELEMENTS::So_base* ele, bool nodal_average)
{
  for (int i = 0; i < nodal_data.M(); ++i)
  {
    const int lid = global_data.Map().LID(ele->NodeIds()[i]);
    if (lid >= 0)  // rownode
    {
      const double invmyadjele = (nodal_average) ? 1.0 / ele->Nodes()[i]->NumElement() : 1.0;
      for (int j = 0; j < nodal_data.N(); ++j)
      {
        (*(global_data(j)))[lid] += nodal_data(i, j) * invmyadjele;
      }
    }
  }
}

template <class T>
void DRT::ELEMENTS::AssembleAveragedElementValues(
    Epetra_MultiVector& global_data, const T& gp_data, const DRT::ELEMENTS::So_base* ele)
{
  const Epetra_BlockMap& elemap = global_data.Map();
  int lid = elemap.LID(ele->Id());
  if (lid != -1)
  {
    for (int i = 0; i < gp_data.N(); ++i)
    {
      double& s = (*(global_data(i)))[lid];  // resolve pointer for faster access
      s = 0.;
      for (unsigned j = 0; j < gp_data.M(); ++j)
      {
        s += gp_data(j, i);
      }
      s *= 1.0 / gp_data.M();
    }
  }
}

void DRT::ELEMENTS::AssembleGaussPointValues(
    std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data,
    const Epetra_SerialDenseMatrix& gp_data, const DRT::ELEMENTS::So_base* ele)
{
  for (int gp = 0; gp < gp_data.M(); ++gp)
  {
    const Epetra_BlockMap& elemap = global_data[gp]->Map();
    int lid = elemap.LID(ele->Id());
    if (lid != -1)
    {
      for (int i = 0; i < gp_data.N(); ++i)
      {
        (*((*global_data[gp])(i)))[lid] += gp_data(gp, i);
      }
    }
  }
}

void DRT::ELEMENTS::AssembleNodalElementCount(
    Epetra_IntVector& global_count, const DRT::ELEMENTS::So_base* ele)
{
  for (int n = 0; n < ele->NumNode(); ++n)
  {
    const int lid = global_count.Map().LID(ele->NodeIds()[n]);

    if (lid != -1)
    {
      global_count[lid] += 1;
    }
  }
}


// explicit template instantiations
template void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(Epetra_MultiVector&,
    const LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D>&, const DRT::ELEMENTS::So_base*, bool);
template void DRT::ELEMENTS::AssembleExtrapolatedNodalValues(
    Epetra_MultiVector&, const LINALG::SerialDenseMatrix&, const DRT::ELEMENTS::So_base*, bool);

template void DRT::ELEMENTS::AssembleAveragedElementValues(Epetra_MultiVector&,
    const LINALG::Matrix<NUMNOD_SOH8, MAT::NUM_STRESS_3D>&, const DRT::ELEMENTS::So_base*);
template void DRT::ELEMENTS::AssembleAveragedElementValues(
    Epetra_MultiVector&, const LINALG::SerialDenseMatrix&, const DRT::ELEMENTS::So_base*);