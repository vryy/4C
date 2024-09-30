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
    Teuchos::RCP<Core::LinAlg::Vector> acinar_vnp_strain;
    Teuchos::RCP<Core::LinAlg::Vector> acinar_vnp;
    Teuchos::RCP<Core::LinAlg::Vector> acinar_vn;
    Teuchos::RCP<Core::LinAlg::Vector> acinar_v;

    Teuchos::RCP<Core::LinAlg::Vector> qin_nm;
    Teuchos::RCP<Core::LinAlg::Vector> qin_n;
    Teuchos::RCP<Core::LinAlg::Vector> qin_np;

    Teuchos::RCP<Core::LinAlg::Vector> x_n;
    Teuchos::RCP<Core::LinAlg::Vector> x_np;
    Teuchos::RCP<Core::LinAlg::Vector> open;

    Teuchos::RCP<Core::LinAlg::Vector> p_extn;
    Teuchos::RCP<Core::LinAlg::Vector> p_extnp;
    Teuchos::RCP<Core::LinAlg::Vector> airway_acinus_dep;
    bool compute_awacinter{};

    Teuchos::RCP<Core::LinAlg::Vector> qout_np;
    Teuchos::RCP<Core::LinAlg::Vector> qout_n;
    Teuchos::RCP<Core::LinAlg::Vector> qout_nm;

    Teuchos::RCP<Core::LinAlg::Vector> p0np;
    Teuchos::RCP<Core::LinAlg::Vector> p0n;
    Teuchos::RCP<Core::LinAlg::Vector> p0nm;
    Teuchos::RCP<Core::LinAlg::Vector> elemArea0;

    Teuchos::RCP<Core::LinAlg::Vector> acini_e_volume;
    Teuchos::RCP<Core::LinAlg::Vector> elemVolume;
    Teuchos::RCP<Core::LinAlg::Vector> elemVolumen;
    Teuchos::RCP<Core::LinAlg::Vector> elemVolumenp;

    Teuchos::RCP<Core::LinAlg::Vector> generations;

    bool solveScatra{};

    Teuchos::RCP<Core::LinAlg::Vector> junVolMix_Corrector;
    Teuchos::RCP<Core::LinAlg::Vector> scatran;
    Teuchos::RCP<Core::LinAlg::Vector> scatranp;

    Teuchos::RCP<Core::LinAlg::Vector> e1scatran;
    Teuchos::RCP<Core::LinAlg::Vector> e2scatran;

    Teuchos::RCP<Core::LinAlg::Vector> e1scatranp;
    Teuchos::RCP<Core::LinAlg::Vector> e2scatranp;

    Teuchos::RCP<Core::LinAlg::Vector> dscatranp;

    Teuchos::RCP<Core::LinAlg::Vector> elemRadiusnp;

    Teuchos::RCP<Core::LinAlg::Vector> cfl;

    Teuchos::RCP<Core::LinAlg::Vector> po2;

    Teuchos::RCP<Core::LinAlg::Vector> bcval;
    Teuchos::RCP<Core::LinAlg::Vector> dbctog;

    Teuchos::RCP<Core::LinAlg::Vector> acini_bc;

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
