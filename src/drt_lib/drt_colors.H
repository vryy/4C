/*----------------------------------------------------------------------*/
/*!
\file drt_colors.H

\brief Color escape sequences for those that prefer colorful output.

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_COLORS_H
#define DRT_COLORS_H

#ifdef COLOROUTPUT
#define BLACK               "[30m"
#define BLACK_LIGHT         "[30;1m"
#define RED                 "[31m"
#define RED_LIGHT           "[31;1m"
#define GREEN               "[32m"
#define GREEN_LIGHT         "[32;1m"
#define YELLOW              "[33m"
#define YELLOW_LIGHT        "[33;1m"
#define BLUE                "[34m"
#define BLUE_LIGHT          "[34;1m"
#define MAGENTA             "[35m"
#define MAGENTA_LIGHT       "[35;1m"
#define BLUE2               "[36m"
#define BLUE2_LIGHT         "[36;1m"
#define GRAY                "[37m"
#define GRAY_LIGHT          "[37;1m"
#define END_COLOR           "[m"
#else
#define BLACK               ""
#define BLACK_LIGHT         ""
#define RED                 ""
#define RED_LIGHT           ""
#define GREEN               ""
#define GREEN_LIGHT         ""
#define YELLOW              ""
#define YELLOW_LIGHT        ""
#define BLUE                ""
#define BLUE_LIGHT          ""
#define MAGENTA             ""
#define MAGENTA_LIGHT       ""
#define BLUE2               ""
#define BLUE2_LIGHT         ""
#define GRAY                ""
#define GRAY_LIGHT          ""
#define END_COLOR           ""
#endif

#endif
