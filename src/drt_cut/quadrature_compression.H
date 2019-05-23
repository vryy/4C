/*----------------------------------------------------------------------*/
/*!
\file quadrature_compression.H
\brief see paper by Sudhakar

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

\level 2
 */
/*----------------------------------------------------------------------*/


#ifndef QUADRATURE_COMPRESSION_H_
#define QUADRATURE_COMPRESSION_H_

#include "../drt_fem_general/drt_utils_gausspoints.H"

namespace Teuchos
{
  template <typename OrdinalType, typename ScalarType>
  class SerialDenseMatrix;

  template <typename OrdinalType, typename ScalarType>
  class SerialDenseVector;

  template <class T>
  class RCP;
}  // namespace Teuchos


namespace GEO
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
          DRT::UTILS::GaussPointsComposite& gin, GEO::CUT::VolumeCell* vc);

      // DRT::UTILS::GaussIntegration GetCompressedQuadrature(){ return *gout_; }
      Teuchos::RCP<DRT::UTILS::GaussPoints> GetCompressedQuadrature() { return gout_; }

     private:
      void FormMatrixSystem(DRT::UTILS::GaussPointsComposite& gin,
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

      bool Compress_Leja_points(DRT::UTILS::GaussPointsComposite& gin,
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int, double>>& mat,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol);

      void ComputeAndPrintError(DRT::UTILS::GaussPointsComposite& gin,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& rhs,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol, std::vector<int>& work,
          int& na);

      void GetPivotalRows(
          Teuchos::RCP<Teuchos::SerialDenseVector<int, int>>& work_temp, std::vector<int>& work);

      bool isThisValueAlreadyInDenseVector(
          int& input, std::vector<int>& vec, int upper_range, int& index);

      Teuchos::RCP<DRT::UTILS::GaussPoints> FormNewQuadratureRule(
          DRT::UTILS::GaussPointsComposite& gin,
          Teuchos::RCP<Teuchos::SerialDenseVector<int, double>>& sol, std::vector<int>& work,
          int& na);
      int GetCorrectIndex(int& input, std::vector<int>& vec, int upper_range);

      void WriteCompressedQuadratureGMSH(
          DRT::UTILS::GaussPointsComposite& gin, GEO::CUT::VolumeCell* vc);

      void IntegratePredefinedPolynomials(DRT::UTILS::GaussPointsComposite& gin);

      Teuchos::RCP<DRT::UTILS::GaussPoints> gout_;
    };

  }  // namespace CUT
}  // namespace GEO
#endif /* QUADRATURE_COMPRESSION_H_ */
