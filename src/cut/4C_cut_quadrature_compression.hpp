/*----------------------------------------------------------------------*/
/*! \file
\brief see paper by Sudhakar


\level 2
 */
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_CUT_QUADRATURE_COMPRESSION_HPP
#define FOUR_C_CUT_QUADRATURE_COMPRESSION_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_gausspoints.hpp"

namespace Teuchos
{
  template <typename OrdinalType, typename ScalarType>
  class SerialDenseMatrix;

  template <typename OrdinalType, typename ScalarType>
  class SerialDenseVector;

  template <class T>
  class RCP;
}  // namespace Teuchos

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class VolumeCell;


    class QuadratureCompression
    {
     public:
      QuadratureCompression();

      bool perform_compression_of_quadrature(
          CORE::FE::GaussPointsComposite& gin, CORE::GEO::CUT::VolumeCell* vc);

      // CORE::FE::GaussIntegration get_compressed_quadrature(){ return *gout_; }
      Teuchos::RCP<CORE::FE::GaussPoints> get_compressed_quadrature() { return gout_; }

     private:
      void FormMatrixSystem(CORE::FE::GaussPointsComposite& gin,
          Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs);

      void Teuchos_GELS(Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol);

      /*!
      \brief Solve the under/over determined system by performing QR decomposition, achived using
      Teuchos framework
       */
      void qr_decomposition_teuchos(Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol);

      /*!
      \brief Solve the under/over determined system by performing QR decomposition, achived using
      LAPACK
       */
      void qr_decomposition_lapack(Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol);

      bool compress_leja_points(CORE::FE::GaussPointsComposite& gin,
          Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol);

      void compute_and_print_error(CORE::FE::GaussPointsComposite& gin,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol, std::vector<int>& work, int& na);

      void GetPivotalRows(
          Teuchos::RCP<CORE::LINALG::IntSerialDenseVector>& work_temp, std::vector<int>& work);

      bool is_this_value_already_in_dense_vector(
          int& input, std::vector<int>& vec, int upper_range, int& index);

      Teuchos::RCP<CORE::FE::GaussPoints> form_new_quadrature_rule(
          CORE::FE::GaussPointsComposite& gin, Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol,
          std::vector<int>& work, int& na);
      int GetCorrectIndex(int& input, std::vector<int>& vec, int upper_range);

      void write_compressed_quadrature_gmsh(
          CORE::FE::GaussPointsComposite& gin, CORE::GEO::CUT::VolumeCell* vc);

      void integrate_predefined_polynomials(CORE::FE::GaussPointsComposite& gin);

      Teuchos::RCP<CORE::FE::GaussPoints> gout_;
    };

  }  // namespace CUT
}  // namespace CORE::GEO
FOUR_C_NAMESPACE_CLOSE

#endif
