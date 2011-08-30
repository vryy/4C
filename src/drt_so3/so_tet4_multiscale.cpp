/*!----------------------------------------------------------------------
\file so_tet4_multiscale.cpp
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
#include "so_tet4.H"
#include "../drt_mat/micromaterial.H"

using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_tet4::sotet4_homog(ParameterList&  params)
{
  double homogdens = 0.;
  const static vector<double> weights = so_tet4_1gp_weights();
  const double density = Material()->Density();

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    homogdens += V_ * weights[gp] * density;
  }

  double homogdensity = params.get<double>("homogdens", 0.0);
  params.set("homogdens", homogdensity+homogdens);

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                              lw 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::sotet4_read_restart_multi()
{
  const int ele_ID = Id();
  RefCountPtr<MAT::Material> mat = Material();

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

    micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_readrestart");
  }
  return;
}


/*----------------------------------------------------------------------*
 |  New result files on the microscale                          lw 08/11|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::sotet4_multi_invana_init()
{
  const int ele_ID = Id();
  RefCountPtr<MAT::Material> mat = Material();

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

    micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_invana_init");
  }
  return;
}

#endif
#endif
