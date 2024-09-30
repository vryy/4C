/*---------------------------------------------------------------------*/
/*! \file
\brief Data for airway elements
\level 3
*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_RED_AIRWAYS_EVALUATION_DATA_HPP
#define FOUR_C_RED_AIRWAYS_EVALUATION_DATA_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <Teuchos_RCP.hpp>


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
    Teuchos::RCP<Core::LinAlg::Vector<double>> acinar_vnp_strain;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acinar_vnp;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acinar_vn;
    Teuchos::RCP<Core::LinAlg::Vector<double>> acinar_v;

    Teuchos::RCP<Core::LinAlg::Vector<double>> qin_nm;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qin_n;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qin_np;

    Teuchos::RCP<Core::LinAlg::Vector<double>> x_n;
    Teuchos::RCP<Core::LinAlg::Vector<double>> x_np;
    Teuchos::RCP<Core::LinAlg::Vector<double>> open;

    Teuchos::RCP<Core::LinAlg::Vector<double>> p_extn;
    Teuchos::RCP<Core::LinAlg::Vector<double>> p_extnp;
    Teuchos::RCP<Core::LinAlg::Vector<double>> airway_acinus_dep;
    bool compute_awacinter{};

    Teuchos::RCP<Core::LinAlg::Vector<double>> qout_np;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qout_n;
    Teuchos::RCP<Core::LinAlg::Vector<double>> qout_nm;

    Teuchos::RCP<Core::LinAlg::Vector<double>> p0np;
    Teuchos::RCP<Core::LinAlg::Vector<double>> p0n;
    Teuchos::RCP<Core::LinAlg::Vector<double>> p0nm;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemArea0;

    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_e_volume;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemVolume;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemVolumen;
    Teuchos::RCP<Core::LinAlg::Vector<double>> elemVolumenp;

    Teuchos::RCP<Core::LinAlg::Vector<double>> generations;

    bool solveScatra{};

    Teuchos::RCP<Core::LinAlg::Vector<double>> junVolMix_Corrector;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatran;
    Teuchos::RCP<Core::LinAlg::Vector<double>> scatranp;

    Teuchos::RCP<Core::LinAlg::Vector<double>> e1scatran;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e2scatran;

    Teuchos::RCP<Core::LinAlg::Vector<double>> e1scatranp;
    Teuchos::RCP<Core::LinAlg::Vector<double>> e2scatranp;

    Teuchos::RCP<Core::LinAlg::Vector<double>> dscatranp;

    Teuchos::RCP<Core::LinAlg::Vector<double>> elemRadiusnp;

    Teuchos::RCP<Core::LinAlg::Vector<double>> cfl;

    Teuchos::RCP<Core::LinAlg::Vector<double>> po2;

    Teuchos::RCP<Core::LinAlg::Vector<double>> bcval;
    Teuchos::RCP<Core::LinAlg::Vector<double>> dbctog;

    Teuchos::RCP<Core::LinAlg::Vector<double>> acini_bc;

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
