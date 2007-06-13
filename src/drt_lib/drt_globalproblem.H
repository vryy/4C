/*----------------------------------------------------------------------*/
/*!
\file drt_globalproblem.H

\brief global list of problems

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef DRT_GLOBALPROBLEM_H
#define DRT_GLOBALPROBLEM_H

#include <vector>
#include <Teuchos_RefCountPtr.hpp>

#include "drt_discret.H"
#include "drt_inputreader.H"

#ifdef PARALLEL
#include <mpi.h>
#endif
extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

namespace DRT
{

/*----------------------------------------------------------------------*/
/// Global problem instance that keeps the discretizations
/*!
  The point of the Problem class is to store the discretizations in a
  structure that parallels the old ccarat field variables. (The fields
  are needed never the less, the old input routines depend on them.)
  So it is the responsibilty of the input functions to setup the
  ccarat fields and the Problem class in the same way.

  This class behaves like a singleton, so there is always one instance
  available. But you can have more than one instance of Problem. In
  normal situations this will not be needed. So don't bother. Just
  call the static Instance() function to get the global instance and
  access your discretizations.

  In the special case that you want to read more that one input file,
  however, you will need to handle the fields from each file
  separately. This is something you cannot do with ccarat's field
  variable so in that case you will have to rely on a second Problem
  instance.

  \author u.kue
  \date 06/07
 */
/*----------------------------------------------------------------------*/
class Problem
{
public:

  /// return an instance of this class
  static RefCountPtr<Problem> Instance(int num=0);

  /// return number of problem instances
  static unsigned NumInstances() { return instances_.size(); }

  /// find a particular discretization
  RefCountPtr<Discretization> Dis(int fieldnum, int disnum) const;

  /// add a discretization to a specific field
  void AddDis(int fieldnum, RefCountPtr<Discretization> dis);

  /// set given discretization to specific field at specific position
  void SetDis(int fieldnum, int disnum, RefCountPtr<Discretization> dis);

  /// add material vector
  void AddMaterial(const _MATERIAL& m);

  /// activate material vector of problem instance as global mat
  void ActivateMaterial();

  /// tell number of known fields
  unsigned NumFields() const { return discretizations_.size(); }

  /// tell number of discretizations in a field
  unsigned NumDis(int fieldnum) const { return discretizations_[fieldnum].size(); }

private:

  /// private default constructor to disallow creation of instances
  Problem() {}

  /// no copying
  Problem(const Problem&);

  /// no copying
  Problem& operator=(const Problem&);

  /// the single instance
  static vector<RefCountPtr<Problem> > instances_;

  /// discretizations of this problem
  vector<vector<RefCountPtr<Discretization> > > discretizations_;

  /// material vector
  vector<_MATERIAL> material_;
};

}

#endif
#endif
