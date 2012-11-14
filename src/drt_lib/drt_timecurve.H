/*----------------------------------------------------------------------*/
/*!
\file drt_timecurve.H

\brief Managing and evaluating of time curves

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
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_TIMECURVE_H
#define DRT_TIMECURVE_H

#include <iostream>
#include <vector>
#include <Teuchos_RCP.hpp>

#include "drt_parser.H"


namespace DRT {

namespace INPUT {

  class Lines;
  class DatFileReader;

}

namespace UTILS {

  class TimeSlice;

  /*====================================================================*/
  /// time curve consisting of a bunch of time slices
  /*!
    A time curve consisting of one or more (adjacent) time
    slices. Each time slice can be of a different type, so you can
    have a total mess of a curve. But honestly, the slice class you
    will normally want to use is ExprTimeSlice. The other ones are
    provided to stay compatible with old code only.

    All time curves are managed by TimeCurveManager.
   */
  class TimeCurve
  {
  public:

    /// evaluate time curve at given time
    /*!
      If time t is less than the starting time of the first slice,
      this starting time is used instead. If time t is larger than the
      ending time of the last slice, this ending time is used.
     */
    double f(double t);

    /// evaluate time curve and its derivatives
    std::vector<double> FctDer(const double t, const unsigned deg);

    /// return the end time of the last time slice of this curve
    double end();

    /// add a new time slice at the and of this curve
    void AddSlice(Teuchos::RCP<TimeSlice> slice);

    void Print(std::ostream& out) const;


  private:

    /// all our time slices ordered by time
    std::vector<Teuchos::RCP<TimeSlice> > slices_;


  };



  /// class that manages all defined time curves
  class TimeCurveManager
  {
  public:

    TimeCurveManager() {}

    /// Definition of valid input section
    Teuchos::RCP<DRT::INPUT::Lines> ValidTimeCurveLines();

    /// read the baci input file and set up all curves
    void ReadInput(DRT::INPUT::DatFileReader& reader);

    /// return a specified time curve
    TimeCurve& Curve(int num);

  private:

    /// no copying
    TimeCurveManager(const TimeCurveManager&);

    /// no copying
    TimeCurveManager& operator=(const TimeCurveManager&);

    /// all curves defined
    std::vector<Teuchos::RCP<TimeCurve> > curves_;
  };

}  // end namespace UTILS
}  // end namespace DRT


void PrintTimeCurveDatHeader();
/// output operator
std::ostream& operator<<(std::ostream& out, const DRT::UTILS::TimeCurve& curve);

#endif
