/*!----------------------------------------------------------------------
\file so_hex27_multiscale.cpp
\brief

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#include "so_hex27.H"
#include "../drt_lib/drt_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/micromaterial.H"
#include "../drt_lib/drt_globalproblem.H"

using namespace std; // cout etc.

extern struct _GENPROB     genprob;



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                                |
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_hex27::soh27_homog(ParameterList&  params)
{
  double homogdens = 0.;
  const static std::vector<double> weights = soh27_weights();
  const double density = Material()->Density();

  for (int gp=0; gp<NUMGPT_SOH27; ++gp)
  {
    homogdens += detJ_[gp] * weights[gp] * density;
  }

  double homogdensity = params.get<double>("homogdens", 0.0);
  params.set("homogdens", homogdensity+homogdens);

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex27::soh27_read_restart_multi()
{
  RefCountPtr<MAT::Material> mat = Material();

  if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
  {
    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm().MyPID()==Owner()) eleowner = true;

    for (int gp=0; gp<NUMGPT_SOH27; ++gp)
      micro->ReadRestart(gp, eleID, eleowner);
  }

  return;
}


#endif
#endif
