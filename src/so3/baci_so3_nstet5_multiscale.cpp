/*----------------------------------------------------------------------*/
/*! \file
\brief multiscale variant of NStet5 element

\level 3

*----------------------------------------------------------------------*/

#include "baci_so3_nstet5.H"
#include "baci_mat_micromaterial.H"
#include "baci_lib_globalproblem.H"
#include "baci_comm_utils.H"
#include "baci_lib_discret.H"



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::NStet5::nstet5_homog(Teuchos::ParameterList& params)
{
  if (DRT::Problem::Instance(0)->GetCommunicators()->SubComm()->MyPID() == Owner())
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
void DRT::ELEMENTS::NStet5::nstet5_read_restart_multi()
{
  const int gp = 0;  // there is only one Gauss point

  Teuchos::RCP<MAT::Material> mat = Material();

  if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
  {
    auto* micro = dynamic_cast<MAT::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (DRT::Problem::Instance()->GetDis("structure")->Comm().MyPID() == Owner()) eleowner = true;

    micro->ReadRestart(gp, eleID, eleowner);
  }

  return;
}
