/*!----------------------------------------------------------------------*
\file combust_dyn.H

\brief control of dynamic combustion analysis

\level 2

<pre>
\maintainer Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

#ifndef COMBUST_DYN_H
#define COMBUST_DYN_H


/*!

  \brief control routine for dynamic combustion analysis

  \detailed

  The approach to simulation of premixed combustion followed by
  the Emmy-Noether-Group at LNM consists of solving a coupled
  two-field problem. One field is the fluid and the second field
  is the level-set function indicating the location of the flame
  interface. The level-set function is a scalar transport equation
  of convection (and diffusion)-type. The scalar is called "G" and
  therefore the level-set function is often called the "G-function".

  Implementing the interaction between both fields with the FEM,
  results in solving both (non-linear) systems of equations separately
  and coupling them via exchanging information at the flame interface.
  Two discretizations are created in BACI, one of type "fluid" and
  another one of type "condif".

  "combust_dyn" is the upper most control routine for the dynamic
  combustion analysis in the DRT::COMBUST module.
  It controls the following tasks:

    o create the G-function discretization from the fluid discretization
    o call a combustion algorithm to solve the multi-field problem
    o check the results

  remark: in this function, the notation is changed from "scatra"-
  style to "gfunc"-style . That is, even though we are still dealing
  with a scalar transport field using scatra elements and so on,
  for convenience the "scarta" discretization is from now on named
  after the G-funciton:

  disnumcdf -> disnumgff
  condifdis -> gfuncdis
  ...

  \author henke
  \date 06/08
*/

void combust_dyn();


#endif /* #ifndef COMBUST_DYN_H */
