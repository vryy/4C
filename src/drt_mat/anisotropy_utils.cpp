/*----------------------------------------------------------------------*/
/*! \file

\brief

\level 3

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
  f1norm = std::sqrt(f1norm);

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

template <typename T, unsigned int numfib>
void MAT::ComputeStructuralTensors(std::vector<std::array<LINALG::Matrix<3, 1>, numfib>>& fibers,
    std::vector<std::array<T, numfib>>& structural_tensor,
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
    for (std::vector<LINALG::Matrix<3, 1>>::size_type i = 0; i < numfib; ++i)
    {
      T A(false);
      strategy->SetupStructuralTensor(fibers[gp].at(i), A);

      structural_tensor[gp].at(i).Update(A);
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

template <typename T, unsigned int numfib>
void MAT::PackFiberArray(DRT::PackBuffer& buffer, const std::vector<std::array<T, numfib>>& vct)
{
  DRT::ParObject::AddtoPack(buffer, static_cast<int>(vct.size()));

  for (const auto& list : vct)
  {
    for (const auto& fiber : list)
    {
      DRT::ParObject::AddtoPack<T::Rows(), T::Cols()>(buffer, fiber);
    }
  }
}

template <typename T>
void MAT::UnpackFiberVector(std::vector<char>::size_type& position, const std::vector<char>& data,
    std::vector<std::vector<T>>& vct)
{
  vct.clear();
  int numgps = DRT::ParObject::ExtractInt(position, data);
  for (int i = 0; i < numgps; ++i)
  {
    std::vector<T> mat(0);
    DRT::ParObject::ExtractfromPack(position, data, mat);
    vct.emplace_back(mat);
  }
}

template <typename T, unsigned int numfib>
void MAT::UnpackFiberArray(std::vector<char>::size_type& position, const std::vector<char>& data,
    std::vector<std::array<T, numfib>>& vct)
{
  vct.clear();
  int numgps = DRT::ParObject::ExtractInt(position, data);
  for (int i = 0; i < numgps; ++i)
  {
    std::array<T, numfib> mat;
    for (unsigned int j = 0; j < numfib; ++j)
    {
      DRT::ParObject::ExtractfromPack<T::Rows(), T::Cols()>(position, data, mat.at(j));
    }
    vct.emplace_back(mat);
  }
}

/*----------------------------------------------------------------------------*/
// explicit instantiation of template functions
template void MAT::ComputeStructuralTensors<LINALG::Matrix<6, 1>, 1u>(
    std::vector<std::array<LINALG::Matrix<3, 1>, 1>>& fibers,
    std::vector<std::array<LINALG::Matrix<6, 1>, 1>>& structural_tensor,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& strategy);
template void MAT::ComputeStructuralTensors<LINALG::Matrix<3, 3>, 1u>(
    std::vector<std::array<LINALG::Matrix<3, 1>, 1>>& fibers,
    std::vector<std::array<LINALG::Matrix<3, 3>, 1>>& structural_tensor,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& strategy);
template void MAT::ComputeStructuralTensors<LINALG::Matrix<6, 1>, 2u>(
    std::vector<std::array<LINALG::Matrix<3, 1>, 2>>& fibers,
    std::vector<std::array<LINALG::Matrix<6, 1>, 2>>& structural_tensor,
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& strategy);
template void MAT::ComputeStructuralTensors<LINALG::Matrix<3, 3>, 2u>(
    std::vector<std::array<LINALG::Matrix<3, 1>, 2>>& fibers,
    std::vector<std::array<LINALG::Matrix<3, 3>, 2>>& structural_tensor,
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

template void MAT::PackFiberArray<LINALG::Matrix<3, 1>, 1u>(
    DRT::PackBuffer& buffer, const std::vector<std::array<LINALG::Matrix<3, 1>, 1>>& vct);
template void MAT::PackFiberArray<LINALG::Matrix<3, 3>, 1u>(
    DRT::PackBuffer& buffer, const std::vector<std::array<LINALG::Matrix<3, 3>, 1>>& vct);
template void MAT::PackFiberArray<LINALG::Matrix<6, 1>, 1u>(
    DRT::PackBuffer& buffer, const std::vector<std::array<LINALG::Matrix<6, 1>, 1>>& vct);

template void MAT::UnpackFiberArray<LINALG::Matrix<3, 1>, 1>(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::array<LINALG::Matrix<3, 1>, 1>>& vct);
template void MAT::UnpackFiberArray<LINALG::Matrix<3, 3>, 1>(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::array<LINALG::Matrix<3, 3>, 1>>& vct);
template void MAT::UnpackFiberArray<LINALG::Matrix<6, 1>, 1>(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::array<LINALG::Matrix<6, 1>, 1>>& vct);

template void MAT::PackFiberArray<LINALG::Matrix<3, 1>, 2u>(
    DRT::PackBuffer& buffer, const std::vector<std::array<LINALG::Matrix<3, 1>, 2>>& vct);
template void MAT::PackFiberArray<LINALG::Matrix<3, 3>, 2u>(
    DRT::PackBuffer& buffer, const std::vector<std::array<LINALG::Matrix<3, 3>, 2>>& vct);
template void MAT::PackFiberArray<LINALG::Matrix<6, 1>, 2u>(
    DRT::PackBuffer& buffer, const std::vector<std::array<LINALG::Matrix<6, 1>, 2>>& vct);

template void MAT::UnpackFiberArray<LINALG::Matrix<3, 1>, 2>(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::array<LINALG::Matrix<3, 1>, 2>>& vct);
template void MAT::UnpackFiberArray<LINALG::Matrix<3, 3>, 2>(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::array<LINALG::Matrix<3, 3>, 2>>& vct);
template void MAT::UnpackFiberArray<LINALG::Matrix<6, 1>, 2>(std::vector<char>::size_type& position,
    const std::vector<char>& data, std::vector<std::array<LINALG::Matrix<6, 1>, 2>>& vct);
