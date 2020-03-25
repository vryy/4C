/*----------------------------------------------------------------------*/
/*! \file

\brief

\level 3

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/
#include "anisotropy_utils.H"
#include "../drt_lib/drt_pack_buffer.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_matelast/elast_aniso_structuraltensor_strategy.H"
#include "../drt_lib/drt_inputreader.H"

void MAT::ReadAnisotropyFiber(
    DRT::INPUT::LineDefinition* linedef, std::string specifier, LINALG::Matrix<3, 1>& fiber_vector)
{
  std::vector<double> fiber;
  linedef->ExtractDoubleVector(std::move(specifier), fiber);
  double f1norm = 0.;
  // normalization
  for (std::vector<double>::size_type i = 0; i < 3; ++i)
  {
    f1norm += fiber[i] * fiber[i];
  }
  f1norm = sqrt(f1norm);

  if (f1norm < 1e-9)
  {
    dserror("The given fiber is not a vector but zero.");
  }

  // fill final fiber vector
  for (std::vector<double>::size_type i = 0; i < 3; ++i)
  {
    fiber_vector(i) = fiber[i] / f1norm;
  }
}

template <typename T>
void MAT::ComputeStructuralTensors(std::vector<std::vector<LINALG::Matrix<3, 1>>>& fibers,
    std::vector<std::vector<T>>& structural_tensor,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& strategy)
{
  // Need to compute the stuctural tensors
  if (Teuchos::is_null(strategy))
  {
    dserror("Structural tensor strategy is null!");
  }

  structural_tensor.resize(fibers.size());
  for (std::vector<std::vector<LINALG::Matrix<3, 1>>>::size_type gp = 0; gp < fibers.size(); ++gp)
  {
    structural_tensor[gp].resize(fibers[gp].size());
    for (std::vector<LINALG::Matrix<3, 1>>::size_type i = 0; i < fibers[gp].size(); ++i)
    {
      T A(false);
      strategy->SetupStructuralTensor(fibers[gp][i], A);

      structural_tensor[gp][i].Update(A);
    }
  }
}

template <typename T>
void MAT::PackFiberVector(DRT::PackBuffer& buffer, const std::vector<std::vector<T>>& vct)
{
  DRT::ParObject::AddtoPack(buffer, static_cast<int>(vct.size()));

  for (const auto& list : vct)
  {
    DRT::ParObject::AddtoPack(buffer, list);
  }
}

template <typename T>
void MAT::UnpackFiberVector(std::vector<char>::size_type& position, const std::vector<char>& data,
    std::vector<std::vector<T>>& vct)
{
  vct.clear();
  int numfibs = DRT::ParObject::ExtractInt(position, data);
  for (int i = 0; i < numfibs; ++i)
  {
    std::vector<T> mat(0);
    DRT::ParObject::ExtractfromPack(position, data, mat);
    vct.emplace_back(mat);
  }
}

/*----------------------------------------------------------------------------*/
// explicit instantiation of template functions
template void MAT::ComputeStructuralTensors(std::vector<std::vector<LINALG::Matrix<3, 1>>>& fibers,
    std::vector<std::vector<LINALG::Matrix<6, 1>>>& structural_tensor,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& strategy);
template void MAT::ComputeStructuralTensors(std::vector<std::vector<LINALG::Matrix<3, 1>>>& fibers,
    std::vector<std::vector<LINALG::Matrix<3, 3>>>& structural_tensor,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& strategy);
template void MAT::PackFiberVector(
    DRT::PackBuffer& buffer, const std::vector<std::vector<LINALG::Matrix<3, 3>>>& vct);
template void MAT::UnpackFiberVector(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::vector<LINALG::Matrix<3, 3>>>& vct);
template void MAT::PackFiberVector(
    DRT::PackBuffer& buffer, const std::vector<std::vector<LINALG::Matrix<3, 1>>>& vct);
template void MAT::UnpackFiberVector(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::vector<LINALG::Matrix<3, 1>>>& vct);
template void MAT::PackFiberVector(
    DRT::PackBuffer& buffer, const std::vector<std::vector<LINALG::Matrix<6, 1>>>& vct);
template void MAT::UnpackFiberVector(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::vector<LINALG::Matrix<6, 1>>>& vct);
