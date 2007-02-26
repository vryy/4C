/*----------------------------------------------------------------------*/
/*!
\file io_drt.H

\class OUT_DRT

\brief output context of one discretization

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef IO_DRT_H
#define IO_DRT_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "../discret/drt_discret.H"
#include "Teuchos_RefCountPtr.hpp"
#include "io.h"

#ifdef BINIO
#include <H5LT.h>
#endif

extern "C" {
#include "../headers/standardtypes.h"
}

#include "../discret/dstrc.H" /* needs standardtypes.h */


/*!
  \class OUT_DRT

  \brief The output context of a discretization

  Create an object of this class for every discretization those mesh
  and results you want to write. Data are written in parallel to
  processor local files. The first process additionally maintains the
  (plain text) control file that glues all result files together.

  \author m.kue
  \date 02/07
*/
class OUT_DRT
{
public:

  // @{ \name Constructors and destructors

  OUT_DRT(RefCountPtr<DRT::Discretization> dis);
  ~OUT_DRT();

  //@}

  // @{ \name Output methods

  //! write result header to control file
  // You will want to call this once each time step _before_ the
  // result data is written.
  // \param step : current time step
  // \param time : current absolute time
  void new_step(int step, double time);

  //! write a result vector
  //
  // There will be an entry in the current result step in the control
  // file that points to this vector
  //
  // \param name : control file entry name
  // \param vec  : the result data vector
  void write_vector(string name, RefCountPtr<Epetra_Vector> vec);

  //! write new "field" group to control file including node and
  // element chunks
  void write_mesh(int step, double time);

  //@}

private:

  //! open new mesh file
  void create_mesh_file(int step);

  //! open new result file
  void create_result_file(int step);

  //! find position of my discretization in (old) ccarat arrays
  void findposition(int& field_pos, unsigned int& dis);

  //! my discretization
  RefCountPtr<DRT::Discretization> dis_;

  unsigned int disnum_;
  int field_pos_;
  int step_;
  double time_;

#ifdef BINIO
  hid_t meshfile_;
  hid_t resultfile_;
  string meshfilename_;
  string resultfilename_;
  hid_t meshgroup_;
  hid_t resultgroup_;
#endif

  int resultfile_changed_;
  int meshfile_changed_;
};

#endif
#endif
#endif
