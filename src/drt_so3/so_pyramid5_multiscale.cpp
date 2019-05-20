/*!----------------------------------------------------------------------
\file so_pyramid5_multiscale.cpp
\brief multiscale functionality for pyramid shaped solid element
\level 2

<pre>
\maintainer Christoph Meier
</pre>

*----------------------------------------------------------------------*/


#include "so_pyramid5.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_discret.H"



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                                |
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_pyramid5::sop5_homog(Teuchos::ParameterList& params)
{
  if (DRT::Problem::Instance(0)->GetNPGroup()->SubComm()->MyPID() == Owner())
  {
    double homogdens = 0.;
    const static std::vector<double> weights = sop5_weights();

    for (int gp = 0; gp < NUMGPT_SOP5; ++gp)
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
void DRT::ELEMENTS::So_pyramid5::sop5_read_restart_multi()
{
  Teuchos::RCP<MAT::Material> mat = Material();

  if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
  {
    MAT::MicroMaterial* micro = static_cast<MAT::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (DRT::Problem::Instance()->GetDis("structure")->Comm().MyPID() == Owner()) eleowner = true;

    for (int gp = 0; gp < NUMGPT_SOP5; ++gp) micro->ReadRestart(gp, eleID, eleowner);
  }

  return;
}
