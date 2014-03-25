/*----------------------------------------------------------------------*/
/*!
\file post_drt_vtu.H

\brief VTU filter

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*/
/*----------------------------------------------------------------------*/
#ifndef POST_DRT_VTU_H
#define POST_DRT_VTU_H


#include "../drt_lib/drt_dserror.H"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>


// forward declarations
class PostField;
class PostResult;
namespace DRT
{
  class Discretization;
  class Node;
}

//! defines different result types (node-based, element-based, Gauss-Points?)
enum ResultType
{
  no_restype,          // unknown result type
  dofbased,            // result values based on dofrowmap
  nodebased,           // result values based on noderowmap
  elementbased,        // result values based on elementrowmap
  elementdof,          // result values based on elementdofrowmap (unused)
  max_restype          // end marker. must be the last entry
};

/*
 \brief Base class for VTU output generation

 \author kronbichler
 \date 03/14
*/
class VtuWriter
{
public:
  //! constructor
  VtuWriter(PostField* field,
            const std::string filename);

  //! destructor
  virtual ~VtuWriter()
  {}

  //! write the whole thing
  void WriteFiles();

  //! write the whole thing
  void WriteFilesChangingGeom()
  {
    CurrentlyNotImplemented();
  }

protected:

  //! look for problem dependent result entries and write them
  virtual void WriteAllResults(PostField* field) = 0;

  /*!
   \brief write one time step of a result

   The results are taken from a reconstructed
   Epetra_Vector. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part.
   */
  void WriteResult(
    const std::string groupname, ///< name of the result group in the control file
    const std::string name, ///< name of the result to be written
    const ResultType restype, ///< type of the result to be written (nodal-/element-based)
    const int numdf, ///< number of dofs per node to this result
    const int from=0, ///< start position of values in nodes
    const bool fillzeros=false ///< zeros are filled to vtu file when no data is available
    );

  /*!
   \brief write the geometry of one time step
   */
  void WriteGeo();

  /*!
   \brief write prologue of the VTU file
   */
  void WriteVtuHeader ();

  /*!
   \brief write epilogue of the VTU file
   */
  void WriteVtuFooter ();

  /// write all element based results
  /*!
    \note This method sees only those element results that are defined in
    the first result group.
   */
  void WriteElementResults(PostField* field)
  {
    if (myrank_ == 0)
      std::cout << "WriteElementResults currently not implemented" << std::endl;
  }

  /*!
   \brief writes debug output as long as this filter is incomplete (TODO MK: should go away at some point)
   */
  void CurrentlyNotImplemented () const
  {
    dserror("Functionality currently not implemented");
  }

  //! Current post field
  PostField* field_;

  //! Current file name
  const std::string filename_;

  //! Rank of current processor; only processor 0 does the output
  const int myrank_;

  //! Output stream for current time step
  std::ofstream currentout_;

  //! Time value for the current time step
  double time_;

  //! Value of the output step for searching the data base
  int timestep_;

  //! Current cycle step (e.g. in nonlinear iteration, not used yet)
  int cycle_;
};

#endif
