/*----------------------------------------------------------------------*/
/*! \file
\brief see paper by Sudhakar


\level 2
 */
/*----------------------------------------------------------------------*/


#ifndef BACI_CUT_QUADRATURE_COMPRESSION_H
#define BACI_CUT_QUADRATURE_COMPRESSION_H

#include "baci_config.H"

#include "baci_discretization_fem_general_utils_gausspoints.H"

namespace Teuchos
{
  template <typename OrdinalType, typename ScalarType>
  class SerialDenseMatrix;

  template <typename OrdinalType, typename ScalarType>
  class SerialDenseVector;

  template <class T>
  class RCP;
}  // namespace Teuchos

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class VolumeCell;


    class QuadratureCompression
    {
     public:
      QuadratureCompression();

      bool PerformCompressionOfQuadrature(
          CORE::FE::GaussPointsComposite& gin, CORE::GEO::CUT::VolumeCell* vc);

      // CORE::FE::GaussIntegration GetCompressedQuadrature(){ return *gout_; }
      Teuchos::RCP<CORE::FE::GaussPoints> GetCompressedQuadrature() { return gout_; }

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
      void QR_decomposition_Teuchos(Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol);

      /*!
      \brief Solve the under/over determined system by performing QR decomposition, achived using
      LAPACK
       */
      void QR_decomposition_LAPACK(Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol);

      bool Compress_Leja_points(CORE::FE::GaussPointsComposite& gin,
          Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& mat,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol);

      void ComputeAndPrintError(CORE::FE::GaussPointsComposite& gin,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& rhs,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol, std::vector<int>& work, int& na);

      void GetPivotalRows(
          Teuchos::RCP<CORE::LINALG::IntSerialDenseVector>& work_temp, std::vector<int>& work);

      bool isThisValueAlreadyInDenseVector(
          int& input, std::vector<int>& vec, int upper_range, int& index);

      Teuchos::RCP<CORE::FE::GaussPoints> FormNewQuadratureRule(CORE::FE::GaussPointsComposite& gin,
          Teuchos::RCP<CORE::LINALG::SerialDenseVector>& sol, std::vector<int>& work, int& na);
      int GetCorrectIndex(int& input, std::vector<int>& vec, int upper_range);

      void WriteCompressedQuadratureGMSH(
          CORE::FE::GaussPointsComposite& gin, CORE::GEO::CUT::VolumeCell* vc);

      void IntegratePredefinedPolynomials(CORE::FE::GaussPointsComposite& gin);

      Teuchos::RCP<CORE::FE::GaussPoints> gout_;
    };

  }  // namespace CUT
}  // namespace CORE::GEO
BACI_NAMESPACE_CLOSE

#endif  // CUT_QUADRATURE_COMPRESSION_H
