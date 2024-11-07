// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_so3_hex8.hpp"
#include "4C_so3_surface.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void Discret::Elements::SoHex8::soh8_homog(Teuchos::ParameterList& params)
{
  if (Global::Problem::instance(0)->get_communicators()->sub_comm()->MyPID() == owner())
  {
    double homogdens = 0.;
    const static std::vector<double> weights = soh8_weights();

    for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
    {
      const double density = material()->density(gp);
      homogdens += detJ_[gp] * weights[gp] * density;
    }

    double homogdensity = params.get<double>("homogdens", 0.0);
    params.set("homogdens", homogdensity + homogdens);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Set EAS internal variables on the microscale (public)       lw 04/08|
 *----------------------------------------------------------------------*/
// the microscale internal EAS data have to be saved separately for every
// macroscopic Gauss point and set before the determination of microscale
// stiffness etc.

void Discret::Elements::SoHex8::soh8_set_eas_multi(Teuchos::ParameterList& params)
{
  if (eastype_ != soh8_easnone)
  {
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldalpha =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldalpha", nullptr);
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldfeas =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldfeas", nullptr);
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldKaainv =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldKaainv", nullptr);
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldKda =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldKda", nullptr);

    if (oldalpha == nullptr || oldfeas == nullptr || oldKaainv == nullptr || oldKda == nullptr)
      FOUR_C_THROW("Cannot get EAS internal data from parameter list for multi-scale problems");

    easdata_.alpha = *(*oldalpha)[id()];
    easdata_.feas = *(*oldfeas)[id()];
    easdata_.invKaa = *(*oldKaainv)[id()];
    easdata_.Kda = *(*oldKda)[id()];
  }
}


/*----------------------------------------------------------------------*
 |  Initialize EAS internal variables on the microscale         lw 03/08|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex8::soh8_eas_init_multi(Teuchos::ParameterList& params)
{
  if (eastype_ != soh8_easnone)
  {
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> lastalpha =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "lastalpha", nullptr);
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldalpha =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldalpha", nullptr);
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldfeas =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldfeas", nullptr);
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldKaainv =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldKaainv", nullptr);
    std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>> oldKda =
        params
            .get<std::shared_ptr<std::map<int, std::shared_ptr<Core::LinAlg::SerialDenseMatrix>>>>(
                "oldKda", nullptr);

    (*lastalpha)[id()] = std::make_shared<Core::LinAlg::SerialDenseMatrix>(neas_, 1);
    (*oldalpha)[id()] = std::make_shared<Core::LinAlg::SerialDenseMatrix>(neas_, 1);
    (*oldfeas)[id()] = std::make_shared<Core::LinAlg::SerialDenseMatrix>(neas_, 1);
    (*oldKaainv)[id()] = std::make_shared<Core::LinAlg::SerialDenseMatrix>(neas_, neas_);
    (*oldKda)[id()] = std::make_shared<Core::LinAlg::SerialDenseMatrix>(neas_, NUMDOF_SOH8);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                              lw 05/08|
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex8::soh8_read_restart_multi()
{
  std::shared_ptr<Core::Mat::Material> mat = material();

  if (mat->material_type() == Core::Materials::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<Mat::MicroMaterial*>(mat.get());
    int eleID = id();
    bool eleowner = false;
    if (Global::Problem::instance()->get_dis("structure")->get_comm().MyPID() == owner())
      eleowner = true;

    for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp) micro->read_restart(gp, eleID, eleowner);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
