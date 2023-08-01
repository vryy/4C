/*---------------------------------------------------------------------*/
/*! \file
\brief Data for airway elements
\level 3
*/
/*---------------------------------------------------------------------*/

#ifndef BACI_RED_AIRWAYS_ELEM_PARAMS_H
#define BACI_RED_AIRWAYS_ELEM_PARAMS_H

namespace DRT::REDAIRWAYS
{
  struct ElemParams
  {
    double qout_np;
    double qout_n;
    double qout_nm;
    double qin_np;
    double qin_n;
    double qin_nm;
    double volnp;
    double voln;
    double acin_vnp;
    double acin_vn;
    double lungVolume_np;
    double lungVolume_n;
    double lungVolume_nm;
    double x_np;
    double x_n;
    double p_extn;
    double p_extnp;
    double open;
  };
}  // namespace DRT::REDAIRWAYS

#endif
