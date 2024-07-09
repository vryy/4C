/*----------------------------------------------------------------------*/
/*! \file
\brief multiscale variant of NStet5 element

\level 3

*----------------------------------------------------------------------*/

#include "4C_comm_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
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
  if (Global::Problem::instance(0)->get_communicators()->sub_comm()->MyPID() == owner())
  {
    const double density = material()->density(0);

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

  Teuchos::RCP<Core::Mat::Material> mat = material();

  if (mat->material_type() == Core::Materials::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<Mat::MicroMaterial*>(mat.get());
    int eleID = id();
    bool eleowner = false;
    if (Global::Problem::instance()->get_dis("structure")->get_comm().MyPID() == owner())
      eleowner = true;

    micro->read_restart(gp, eleID, eleowner);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
