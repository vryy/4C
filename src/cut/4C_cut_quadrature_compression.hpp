// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_CUT_QUADRATURE_COMPRESSION_HPP
#define FOUR_C_CUT_QUADRATURE_COMPRESSION_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class SerialDenseMatrix;
}

namespace Cut
{
  class VolumeCell;


  class QuadratureCompression
  {
   public:
    QuadratureCompression();

    bool perform_compression_of_quadrature(
        Core::FE::GaussPointsComposite& gin, Cut::VolumeCell* vc);

    // Core::FE::GaussIntegration get_compressed_quadrature(){ return *gout_; }
    std::shared_ptr<Core::FE::GaussPoints> get_compressed_quadrature() { return gout_; }

   private:
    void form_matrix_system(Core::FE::GaussPointsComposite& gin,
        Core::LinAlg::SerialDenseMatrix& mat, Core::LinAlg::SerialDenseVector& rhs);

    void teuchos_gels(std::shared_ptr<Core::LinAlg::SerialDenseMatrix>& mat,
        std::shared_ptr<Core::LinAlg::SerialDenseVector>& rhs,
        Core::LinAlg::SerialDenseVector& sol);

    /*!
    \brief Solve the under/over determined system by performing QR decomposition, achived using
    Teuchos framework
     */
    void qr_decomposition_teuchos(std::shared_ptr<Core::LinAlg::SerialDenseMatrix>& mat,
        std::shared_ptr<Core::LinAlg::SerialDenseVector>& rhs,
        std::shared_ptr<Core::LinAlg::SerialDenseVector>& sol);

    /*!
    \brief Solve the under/over determined system by performing QR decomposition, achived using
    LAPACK
     */
    void qr_decomposition_lapack(std::shared_ptr<Core::LinAlg::SerialDenseMatrix>& mat,
        std::shared_ptr<Core::LinAlg::SerialDenseVector>& rhs,
        Core::LinAlg::SerialDenseVector& sol);

    bool compress_leja_points(Core::FE::GaussPointsComposite& gin,
        std::shared_ptr<Core::LinAlg::SerialDenseMatrix>& mat,
        std::shared_ptr<Core::LinAlg::SerialDenseVector>& rhs,
        std::shared_ptr<Core::LinAlg::SerialDenseVector>& sol);

    void compute_and_print_error(Core::FE::GaussPointsComposite& gin,
        Core::LinAlg::SerialDenseVector& rhs, Core::LinAlg::SerialDenseVector& sol,
        std::vector<int>& work, int& na);

    void get_pivotal_rows(Core::LinAlg::IntSerialDenseVector& work_temp, std::vector<int>& work);

    bool is_this_value_already_in_dense_vector(
        int& input, std::vector<int>& vec, int upper_range, int& index);

    std::shared_ptr<Core::FE::GaussPoints> form_new_quadrature_rule(
        Core::FE::GaussPointsComposite& gin, Core::LinAlg::SerialDenseVector& sol,
        std::vector<int>& work, int& na);
    int get_correct_index(int& input, std::vector<int>& vec, int upper_range);

    void write_compressed_quadrature_gmsh(Core::FE::GaussPointsComposite& gin, Cut::VolumeCell* vc);

    void integrate_predefined_polynomials(Core::FE::GaussPointsComposite& gin);

    std::shared_ptr<Core::FE::GaussPoints> gout_;
  };

}  // namespace Cut

FOUR_C_NAMESPACE_CLOSE

#endif
