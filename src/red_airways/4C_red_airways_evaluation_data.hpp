// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_RED_AIRWAYS_EVALUATION_DATA_HPP
#define FOUR_C_RED_AIRWAYS_EVALUATION_DATA_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <memory>


FOUR_C_NAMESPACE_OPEN

namespace Discret::ReducedLung
{
  /**
   * Store all global vectors that are required to evaluate different kind of reduced airway
   * elements.
   */
  struct EvaluationData
  {
   protected:
    EvaluationData() = default;

   public:
    std::shared_ptr<Core::LinAlg::Vector<double>> acinar_vnp_strain;
    std::shared_ptr<Core::LinAlg::Vector<double>> acinar_vnp;
    std::shared_ptr<Core::LinAlg::Vector<double>> acinar_vn;
    std::shared_ptr<Core::LinAlg::Vector<double>> acinar_v;

    std::shared_ptr<Core::LinAlg::Vector<double>> qin_nm;
    std::shared_ptr<Core::LinAlg::Vector<double>> qin_n;
    std::shared_ptr<Core::LinAlg::Vector<double>> qin_np;

    std::shared_ptr<Core::LinAlg::Vector<double>> x_n;
    std::shared_ptr<Core::LinAlg::Vector<double>> x_np;
    std::shared_ptr<Core::LinAlg::Vector<double>> open;

    std::shared_ptr<Core::LinAlg::Vector<double>> p_extn;
    std::shared_ptr<Core::LinAlg::Vector<double>> p_extnp;
    std::shared_ptr<Core::LinAlg::Vector<double>> airway_acinus_dep;
    bool compute_awacinter{};

    std::shared_ptr<Core::LinAlg::Vector<double>> qout_np;
    std::shared_ptr<Core::LinAlg::Vector<double>> qout_n;
    std::shared_ptr<Core::LinAlg::Vector<double>> qout_nm;

    std::shared_ptr<Core::LinAlg::Vector<double>> p0np;
    std::shared_ptr<Core::LinAlg::Vector<double>> p0n;
    std::shared_ptr<Core::LinAlg::Vector<double>> p0nm;
    std::shared_ptr<Core::LinAlg::Vector<double>> elemArea0;

    std::shared_ptr<Core::LinAlg::Vector<double>> acini_e_volume;
    std::shared_ptr<Core::LinAlg::Vector<double>> elemVolume;
    std::shared_ptr<Core::LinAlg::Vector<double>> elemVolumen;
    std::shared_ptr<Core::LinAlg::Vector<double>> elemVolumenp;

    std::shared_ptr<Core::LinAlg::Vector<double>> generations;

    std::shared_ptr<Core::LinAlg::Vector<double>> elemRadiusnp;

    std::shared_ptr<Core::LinAlg::Vector<double>> cfl;

    std::shared_ptr<Core::LinAlg::Vector<double>> bcval;
    std::shared_ptr<Core::LinAlg::Vector<double>> dbctog;

    std::shared_ptr<Core::LinAlg::Vector<double>> acini_bc;

    double lungVolume_np{};
    double lungVolume_n{};
    double lungVolume_nm{};

    double time{};
    double dt{};

    static EvaluationData& get()
    {
      static EvaluationData evaluation_data;
      return evaluation_data;
    }

    EvaluationData(const EvaluationData&) = delete;
    EvaluationData(EvaluationData&&) = delete;
    EvaluationData& operator=(const EvaluationData&) = delete;
    EvaluationData& operator=(EvaluationData&&) = delete;
  };
}  // namespace Discret::ReducedLung

FOUR_C_NAMESPACE_CLOSE

#endif
