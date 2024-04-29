/*----------------------------------------------------------------------*/
/*! \file
\brief multiscale variant of 3D quadratic serendipity element
\level 2


*----------------------------------------------------------------------*/


#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_so3_hex20.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                                |
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::SoHex20::soh20_homog(Teuchos::ParameterList& params)
{
  if (GLOBAL::Problem::Instance(0)->GetCommunicators()->SubComm()->MyPID() == Owner())
  {
    double homogdens = 0.;
    const static std::vector<double> weights = soh20_weights();

    for (int gp = 0; gp < NUMGPT_SOH20; ++gp)
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
void DRT::ELEMENTS::SoHex20::soh20_read_restart_multi()
{
  Teuchos::RCP<CORE::MAT::Material> mat = Material();

  if (mat->MaterialType() == CORE::Materials::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<MAT::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (GLOBAL::Problem::Instance()->GetDis("structure")->Comm().MyPID() == Owner())
      eleowner = true;

    for (int gp = 0; gp < NUMGPT_SOH20; ++gp) micro->ReadRestart(gp, eleID, eleowner);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
