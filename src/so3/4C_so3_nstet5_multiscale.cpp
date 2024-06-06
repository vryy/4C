/*----------------------------------------------------------------------*/
/*! \file
\brief multiscale variant of NStet5 element

\level 3

*----------------------------------------------------------------------*/

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_so3_nstet5.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void Discret::ELEMENTS::NStet5::nstet5_homog(Teuchos::ParameterList& params)
{
  if (Global::Problem::Instance(0)->GetCommunicators()->SubComm()->MyPID() == Owner())
  {
    const double density = Material()->Density(0);

    double homogdens = V_ * density;

    double homogdensity = params.get<double>("homogdens", 0.0);
    params.set("homogdens", homogdensity + homogdens);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                              lw 05/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::NStet5::nstet5_read_restart_multi()
{
  const int gp = 0;  // there is only one Gauss point

  Teuchos::RCP<Core::Mat::Material> mat = Material();

  if (mat->MaterialType() == Core::Materials::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<Mat::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (Global::Problem::Instance()->GetDis("structure")->Comm().MyPID() == Owner())
      eleowner = true;

    micro->read_restart(gp, eleID, eleowner);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
