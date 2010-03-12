/*!----------------------------------------------------------------------
\file patspec.cpp

\brief A collection of methods to modify patient specific geometries

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "patspec.H"

/*----------------------------------------------------------------------*
 |                                                             gee 03/10|
 *----------------------------------------------------------------------*/
void PATSPEC::PatientSpecificGeometryComputation(DRT::Discretization& dis)
{
  if (!dis.Comm().MyPID())
  {
    cout << "____________________________________________________________\n";
    cout << "Entering patient specific structural preprocessing (PATSPEC)\n";
    cout << "************************************************************\n";
  }
  
  return;
}


#endif
