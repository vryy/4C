/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet4 Element
\level 2
*----------------------------------------------------------------------*/

#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_so3_tet4.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void Discret::ELEMENTS::SoTet4::sotet4_homog(Teuchos::ParameterList& params)
{
  if (Global::Problem::instance(0)->get_communicators()->sub_comm()->MyPID() == owner())
  {
    double homogdens = 0.;
    const static std::vector<double> weights = so_tet4_1gp_weights();

    for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
    {
      const double density = material()->density(gp);
      homogdens += V_ * weights[gp] * density;
    }

    double homogdensity = params.get<double>("homogdens", 0.0);
    params.set("homogdens", homogdensity + homogdens);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                              lw 05/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoTet4::sotet4_read_restart_multi()
{
  Teuchos::RCP<Core::Mat::Material> mat = material();

  if (mat->material_type() == Core::Materials::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<Mat::MicroMaterial*>(mat.get());
    int eleID = id();
    bool eleowner = false;
    if (Global::Problem::instance()->get_dis("structure")->get_comm().MyPID() == owner())
      eleowner = true;

    for (int gp = 0; gp < NUMGPT_SOTET4; ++gp) micro->read_restart(gp, eleID, eleowner);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
