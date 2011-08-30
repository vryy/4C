/*!----------------------------------------------------------------------
\file so_ptet_multiscale.cpp
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
#include "so_ptet.H"
#include "../drt_mat/micromaterial.H"

using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::Ptet::ptet_homog(ParameterList&  params)
{
  const double density = Material()->Density();

  double homogdens = V_ * density;

  double homogdensity = params.get<double>("homogdens", 0.0);
  params.set("homogdens", homogdensity+homogdens);

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                              lw 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::ptet_read_restart_multi()
{
  const int ele_ID = Id();
  const int gp = 0; // there is only one Gauss point

  RefCountPtr<MAT::Material> mat = Material();
  MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
  micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_readrestart");

  return;
}


/*----------------------------------------------------------------------*
 |  New result files on the microscale                          lw 08/11|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::ptet_multi_invana_init()
{
  const int ele_ID = Id();
  const int gp = 0; // there is only one Gauss point

  RefCountPtr<MAT::Material> mat = Material();
  MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
  micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_invana_init");

  return;
}

#endif
#endif
