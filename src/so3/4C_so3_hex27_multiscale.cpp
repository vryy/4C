/*----------------------------------------------------------------------*/
/*! \file
\brief multiscale functionality for tri-quadratic displacement based solid element
\level 2


*----------------------------------------------------------------------*/


#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_so3_hex27.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                                |
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void Discret::ELEMENTS::SoHex27::soh27_homog(Teuchos::ParameterList& params)
{
  if (Global::Problem::Instance(0)->GetCommunicators()->SubComm()->MyPID() == Owner())
  {
    double homogdens = 0.;
    const static std::vector<double> weights = soh27_weights();

    for (int gp = 0; gp < NUMGPT_SOH27; ++gp)
    {
      const double density = Material()->Density(gp);
      homogdens += detJ_[gp] * weights[gp] * density;
    }

    double homogdensity = params.get<double>("homogdens", 0.0);
    params.set("homogdens", homogdensity + homogdens);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                                      |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex27::soh27_read_restart_multi()
{
  Teuchos::RCP<Core::Mat::Material> mat = Material();

  if (mat->MaterialType() == Core::Materials::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<Mat::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (Global::Problem::Instance()->GetDis("structure")->Comm().MyPID() == Owner())
      eleowner = true;

    for (int gp = 0; gp < NUMGPT_SOH27; ++gp) micro->read_restart(gp, eleID, eleowner);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
