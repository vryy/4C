/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet4 Element
\level 2
*----------------------------------------------------------------------*/

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
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
  if (Global::Problem::Instance(0)->GetCommunicators()->SubComm()->MyPID() == Owner())
  {
    double homogdens = 0.;
    const static std::vector<double> weights = so_tet4_1gp_weights();

    for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
    {
      const double density = Material()->Density(gp);
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
  Teuchos::RCP<Core::Mat::Material> mat = Material();

  if (mat->MaterialType() == Core::Materials::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<Mat::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (Global::Problem::Instance()->GetDis("structure")->Comm().MyPID() == Owner())
      eleowner = true;

    for (int gp = 0; gp < NUMGPT_SOTET4; ++gp) micro->read_restart(gp, eleID, eleowner);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
