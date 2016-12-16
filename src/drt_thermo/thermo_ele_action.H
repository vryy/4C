/*----------------------------------------------------------------------*/
/*!
 \file thermo_ele_action.H

 \brief enums defining action of the thermo elements and related helper functions

\level 1

\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
 *----------------------------------------------------------------------*/


#ifndef THERMO_ELE_ACTION_H_
#define THERMO_ELE_ACTION_H_

namespace THR
{

/*--------------------------------------------------------------------------
 | enum that provides all possible thermo actions
 *--------------------------------------------------------------------------*/
enum Action {
    none,
    calc_thermo_fint,
    calc_thermo_fintcapa,
    calc_thermo_finttang,
    calc_thermo_heatflux,
    postproc_thermo_heatflux,
    integrate_shape_functions,
    calc_thermo_update_istep,
    calc_thermo_reset_istep,
    calc_thermo_energy,
    calc_thermo_coupltang,
    calc_thermo_fintcond,
    calc_thermo_fext,
    calc_thermo_error
}; // enum Action

/*--------------------------------------------------------------------------
 | enum that provides all possible thermo actions on a boundary
 *--------------------------------------------------------------------------*/
enum BoundaryAction {
  ba_none,
  calc_thermo_fextconvection,
  calc_thermo_fextconvection_coupltang,
  calc_normal_vectors,
  ba_integrate_shape_functions
};

/*!
 * \brief translate to string for screen output
 */
inline std::string ActionToString(const Action action)
{
    std::string s = "";
    switch (action)
    {
    case none:                      s = "none";  break;
    case calc_thermo_fint:          s = "calc_thermo_fint";  break;
    case calc_thermo_fintcapa:      s = "calc_thermo_fintcapa";  break;
    case calc_thermo_finttang:      s = "calc_thermo_finttang";  break;
    case calc_thermo_heatflux:      s = "calc_thermo_heatflux";  break;
    case postproc_thermo_heatflux:  s = "postproc_thermo_heatflux";  break;
    case integrate_shape_functions: s = "integrate_shape_functions";  break;
    case calc_thermo_update_istep:  s = "calc_thermo_update_istep";  break;
    case calc_thermo_reset_istep:   s = "calc_thermo_reset_istep";  break;
    case calc_thermo_energy:        s = "calc_thermo_energy";  break;
    case calc_thermo_coupltang:     s = "calc_thermo_coupltang";  break;
    default:
        std::cout << action << std::endl;
        dserror("no string for this action defined!");
        break;
    };
    return s;
}

inline std::string BoundaryActionToString(const BoundaryAction baction)
{
    std::string s = "";
    switch (baction)
    {
    case ba_none:                                s = "ba_none";  break;
    case calc_thermo_fextconvection:             s = "calc_thermo_fextconvection";  break;
    case calc_thermo_fextconvection_coupltang:   s = "calc_thermo_fextconvection_coupltang";  break;
    case calc_normal_vectors:                    s = "calc_normal_vectors";  break;
    case ba_integrate_shape_functions:           s = "ba_integrate_shape_functions";  break;
    default:
        std::cout << baction << std::endl;
        dserror("no string for this boundary action defined!");
        break;
    };
    return s;
}


} //namespace THR


#endif /* THERMO_ELE_ACTION_H_ */
