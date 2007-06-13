/*----------------------------------------------------------------------*/
/*!
\file post_drt_common.H

\brief drt binary filter library

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

The drt binary filter library contains classes to open and read binary
files. Based on these different filter applications can be written.

Note: The new filter design uses the discretization, node and element
classes. So postprocessing calculations are easily possible. Element
specific output operations can be done by the elements themselfes.

The new filters can be used in parallel.

*/
/*----------------------------------------------------------------------*/

#ifndef POST_DRT_COMMON_H
#define POST_DRT_COMMON_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <vector>
#include <string>

#include "../io/hdf_reader.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "../drt_lib/drt_discret.H"
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

extern "C" {
#include "../headers/standardtypes.h"
#include "../pss_full/pss_table.h"
#include "../pss_full/pss_table_iter.h"
}

#include "../drt_lib/dstrc.H" /* needs standardtypes.h */

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>


#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

class PostProblem;
class PostField;
class PostResult;
class PostChunk;

using namespace std;
using namespace Teuchos;

// external declarations not found otherwise
extern "C"
{
void distribute_drt_grids();
}


/*!
  \class PostProblem

  \brief Foundation class of the filter

  Each filter needs one instance of PostProblem. This sets up the
  filter environment including MPI in a parallel setting. Some rather
  obscure global ccarat variables are handled here as well.

  \author m.kue
  \date 02/07
*/
class PostProblem
{
public:
  // @{ \name Constructors and destructors

  //! setup filter and read command line
  //
  // \param CLP (in)  : preconfigured command line processor
  // \param argc (in) : number of command line arguments
  // \param argv (in) : command line arguments
  PostProblem(Teuchos::CommandLineProcessor& CLP, int argc, char** argv);
  ~PostProblem();

  //@}

  // @{ \name Access methods

  //! returns the basename of the output files
  string basename() const { return basename_; }

  //! sets the basename of the output files to name
  void set_basename(string name) { basename_ = name; }

  //! returns the number of result groups that are stored in the control file
  int num_results() { return static_cast<int>(result_group_.size()); }

  //! returns a pointer to a vector of all result groups
  vector<MAP*>* result_groups() { return &result_group_; }

  //! returns the number of discretizations
  int num_discr();

  //! returns the number of dimensions the problem has
  int num_dim() { return ndim_; }

  //! returns a pointer to the num-th discretization
  PostField* get_discretization(int num);

  int num_output_procs();

  //! returns the first step in the range, that the user specified
  int start() { return start_; }

  //! returns the last step in the range, that the user specified
  int end() { return end_; }

  //! returns the step parameter the user specified
  int step() { return step_; }

  RefCountPtr<Epetra_Comm> comm();

  //! returns the directory name in which the control files and data files are
  string input_dir() { return input_dir_; }

  //! returns type of the simulation
  PROBLEM_TYP Problemtype() { return problemtype_; }

  //@}

private:

  //! initializes all the data a filter needs.
  //
  // This function is called by the Constructor.
  void setup_filter(const char* output_name);

  //! reads the mesh files
  //
  // calls 'getfield()' for each 'field'-entry in the mesh file
  // (currently it reads only the fields with step ==  0).
  // This function is called by the Constructor.
  void read_meshes();

  //! creates and returns a PostField insance from a field MAP.
  //
  // Keeps track of global field variable
  PostField getfield(MAP* field_info);

  //! type of the problem we post process
  PROBLEM_TYP problemtype_;

  //! number of dimensions
  int ndim_;

  //! output file base name (possibly including directory)
  string basename_;

  //! memory version of control file
  /*! If a restarted run is post processed this is a merged control
    file constructed from all control files available. */
  MAP control_table_;

  //! the meshes of our problem (fields are discretizations here!)
  vector<PostField> fields_;

  /* start, stop and step numbers. a python like slice. */
  /* We don't have to read each result. This is set by command line
   * arguments. */
  int start_;
  int end_;
  int step_;

  //! MAP pointers into the control file map.
  // We don't own these pointers, so no reference counting.
  vector<MAP*> result_group_;

  string input_dir_;
  RefCountPtr<Epetra_Comm> comm_;
};


/*!
  \class PostField

  \brief The representation of one mesh.

  Problems with fixed mesh topologies will want to create one instance
  of this class for each mesh of the problem.

  Problems with changing meshes can be done, but will need some more
  work.

  \author m.kue
  \date 02/07
*/
class PostField
{
public:
  // @{ \name Constructors and destructors

  PostField(string name, RefCountPtr<Epetra_Comm> comm, PostProblem* problem,
            int field_pos, string field_name, FIELDTYP type, int disnum,
            int numnd, int numele);
  PostField(RefCountPtr<DRT::Discretization> dis, PostProblem* problem,
            int field_pos, string field_name, FIELDTYP type, int disnum,
            int numnd, int numele);

  ~PostField();

  //@}

  // @{ \name Access methods

  //! returns the discretization of this field
  RefCountPtr<DRT::Discretization> discretization() { return dis_; }

  //! returns the number of output processors
  int num_output_procs() { return num_output_procs_; }

  //! set the number of output processors
  void set_num_output_procs(int num) { num_output_procs_ = num; }

  //! returns the number of global Dof-Ids
  int global_id_num() { return dis_->DofRowMap()->NumGlobalElements(); }

  //! returns the number of the discretization stored in this field
  int disnum() { return disnum_; }

  //! returns the PostProbem instance
  PostProblem* problem() { return problem_; }

  //! returns the index of this field
  int field_pos() { return field_pos_; }

  //! returns the name of the field
  string name() { return field_name_; }

  //! returns the type of this field
  FIELDTYP type() { return type_; }

  //! returns the number of nodes
  int num_nodes() { return numnd_; }

  //! returns the number of elements
  int num_elements() { return numele_; }

  //@}

private:
  RefCountPtr<DRT::Discretization> dis_;
  PostProblem* problem_;
  int num_output_procs_;
  int field_pos_;
  string field_name_;
  FIELDTYP type_;
  int disnum_;
  int numnd_;
  int numele_;
};


/*!
  \class PostResult

  \brief Representation of the results of one time step.

  Most filters will loop all time steps from the control file.

  \author m.kue
  \date 02/07
*/
class PostResult
{
public:
  // @{ \name Constructors and destructors
  PostResult(PostField* field);
  ~PostResult();

  //@}

  // @{ \name Reader methods

  //! goes to the next result block in the slice
  //
  // Automatically opens new result files.
  //
  // \returns 1 if a new result block has been found, otherwise returns 0
  int next_result();

  //! reads the data of the result vector 'name' from the current result block
  //
  // \param name (in): result vector name inside the control file
  //
  // \returns vector with results read
  RefCountPtr<Epetra_Vector> read_result(string name);

  //@}

  // @{ \name Access methods

  //! returns the field in which this result exists
  PostField* field() { return field_; }

  //! returns the MAP all the result groups are read from
  MAP* group() { return group_; }

  //! returns time of this result
  double time() { return map_read_real(group_, "time"); }

  //! returns step number of this result
  int step() { return map_read_int(group_, "step"); }

  //@}
private:

  //! closes all the currently open result files
  void close_result_files();

  //! opens result files.
  //
  // The name is taken from the "result_file" entry in the block
  // 'field_info'
  //
  // \param field_info (in) : result block map for control file
  void open_result_files(MAP* field_info);

  //! Tell whether a given result group belongs to this result.
  int match_field_result(MAP* result_group);

  PostField* field_;
  int pos_;
  MAP* group_;
  HDFReader file_;
};


inline int PostProblem::num_discr()
{
  return static_cast<int>(fields_.size());
}


#endif
#endif
#endif
