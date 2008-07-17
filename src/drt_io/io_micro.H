/*----------------------------------------------------------------------*/
/*!
 * \file io_micro.H
\brief output context of one discretization

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#if 0

#ifndef IO_DRT_MICRO_H
#define IO_DRT_MICRO_H

#ifdef CCADISCRET

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_BlockMap.h>

#include "io.H"

#include <map>
#include <string>

#include "../drt_lib/drt_discret.H"
//#include "io.h"
#include "io_hdf.H"

#ifdef BINIO
#include <H5LT.h>
#endif

extern "C" {
#include "../pss_full/pss_table.h"
}

namespace IO
{

/*!
  \brief The output context of a discretization at microscale

  Create an object of this class for every discretization those mesh
  and results you want to write. Data are written in parallel to
  processor local files. The first process additionally maintains the
  (plain text) control file that glues all result files together.

  \author m.kue
  \date 02/07
*/
class MicroDiscretizationWriter : public DiscretizationWriter
{
public:

  /// Writer constructor always needs a discretization
  MicroDiscretizationWriter(RefCountPtr<DRT::Discretization> dis, int probnum, int ele, int gp);

  /// cleanup, close hdf5 files
  ~MicroDiscretizationWriter();

  /// get name of restart file
  string GetRestartFileName() { return name_restart_; };


private:

  /// name of restart file if present
  string name_restart_;
};

class MicroDiscretizationReader : public DiscretizationReader
{
public:

  MicroDiscretizationReader(Teuchos::RCP<DRT::Discretization> dis, int step, string name);

  ~MicroDiscretizationReader();

private:

  MAP microfile_;
};

}

#endif
#endif
#endif
