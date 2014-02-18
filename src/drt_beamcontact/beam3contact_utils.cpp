/*!----------------------------------------------------------------------
\file beam3contact_utils.cpp
\brief A set of utility functions for beam contact

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Christoh Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam3contact_utils.H"


  //Cast of FAD to double
  double BEAMCONTACT::CastToDouble(FAD a)
  {
    return a.val();
  }

  //Cast of double to double
  double BEAMCONTACT::CastToDouble(double a)
  {
    return a;
  }

  //Calculate Norm of a scalar FAD or double quantity
  double BEAMCONTACT::Norm(double a)
  {
    return sqrt(a*a);
  }

  //Calculate Norm of a scalar FAD or double quantity
  FAD BEAMCONTACT::Norm(FAD a)
  {
    return pow(a*a,0.5);
  }
