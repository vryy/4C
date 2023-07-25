/*----------------------------------------------------------------------*/
/*! \file
\brief see paper by Sudhakar


\level 2
 */
/*----------------------------------------------------------------------*/


#ifndef BACI_CUT_QUADRATURE_COMPRESSION_H
#define BACI_CUT_QUADRATURE_COMPRESSION_H

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


namespace CORE::GEO
{
  namespace CUT
  {
    class VolumeCell;


    class QuadratureCompression
    {
     public:
      QuadratureCompression();

      ~QuadratureCompression(){};

      bool PerformCompressionOfQuadrature(
          CORE::DRT::UTILS::GaussPointsComposite& gin, CORE::GEO::CUT::VolumeCell* vc);

      // CORE::DRT::UTILS::GaussIntegration GetCompressedQuadrature(){ return *gout_; }
      Teuchos::RCP<CORE::DRT::UTILS::GaussPoints> GetCompressedQuadrature() { return gout_; }

     private:
      void FormMatrixSystem(CORE::DRT::UTILS::GaussPointsComposite& gin,
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>& mat,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs);

      void Teuchos_GELS(Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>& mat,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol);

      /*!
      \brief Solve the under/over determined system by performing QR decomposition, achived using
      Teuchos framework
       */
      void QR_decomposition_Teuchos(Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>& mat,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol);

      /*!
      \brief Solve the under/over determined system by performing QR decomposition, achived using
      LAPACK
       */
      void QR_decomposition_LAPACK(Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>& mat,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol);

      bool Compress_Leja_points(CORE::DRT::UTILS::GaussPointsComposite& gin,
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>& mat,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol);

      void ComputeAndPrintError(CORE::DRT::UTILS::GaussPointsComposite& gin,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol, std::vector<int>& work,
          int& na);

      void GetPivotalRows(
          Teuchos::RCP<Teuchos::SerialDenseVector<int, int>>& work_temp, std::vector<int>& work);

      bool isThisValueAlreadyInDenseVector(
          int& input, std::vector<int>& vec, int upper_range, int& index);

      Teuchos::RCP<CORE::DRT::UTILS::GaussPoints> FormNewQuadratureRule(
          CORE::DRT::UTILS::GaussPointsComposite& gin,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol, std::vector<int>& work,
          int& na);
      int GetCorrectIndex(int& input, std::vector<int>& vec, int upper_range);

      void WriteCompressedQuadratureGMSH(
          CORE::DRT::UTILS::GaussPointsComposite& gin, CORE::GEO::CUT::VolumeCell* vc);

      void IntegratePredefinedPolynomials(CORE::DRT::UTILS::GaussPointsComposite& gin);

      Teuchos::RCP<CORE::DRT::UTILS::GaussPoints> gout_;
    };

  }  // namespace CUT
}  // namespace CORE::GEO
#endif  // CUT_QUADRATURE_COMPRESSION_H
