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
#include "../post_drt_common/post_drt_common.H"
#include "../post_drt_common/post_writer_base.H"

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


/*
 \brief Base class for VTU output generation

 \author kronbichler
 \date 03/14
*/
class VtuWriter : public PostWriterBase
{
public:
  //! constructor. Initializes the writer to a certain field.
  VtuWriter (PostField* field,
             const std::string &name);

  //! destructor
  virtual ~VtuWriter()
  {}

  //! write the whole thing
  virtual void WriteFiles(PostFilterBase &filter);

  //! write the whole thing
  virtual void WriteFilesChangingGeom(PostFilterBase &filter)
  {
    CurrentlyNotImplemented();
  }

protected:

  /*!
   \brief write one time step of a result

   The results are taken from a reconstructed
   Epetra_Vector. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part.
   */
  virtual void WriteResult(
    const std::string groupname, ///< name of the result group in the control file
    const std::string name, ///< name of the result to be written
    const ResultType restype, ///< type of the result to be written (nodal-/element-based)
    const int numdf, ///< number of dofs per node to this result
    const int from=0, ///< start position of values in nodes
    const bool fillzeros=false ///< zeros are filled to vtu file when no data is available
    );


  virtual void WriteResultOneTimeStep(
            PostResult& result,   ///< result group in the control file
            const std::string groupname, ///< name of the result group in the control file
            const std::string name, ///< name of the result to be written
            const ResultType restype, ///< type of the result to be written (nodal-/element-based)
            const int numdf, ///< number of dofs per node to this result
            bool firststep, ///< bool whether this is the first time step
            bool laststep,  ///< bool whether this is the last time step
            const int from=0 ///< start position of values in nodes
            )
  {
    dserror("Not yet implemented");
  }

  /*!
   \brief write a particular variable to file

   Write results. Some variables need interaction with the post filter,
   e.g. structural stresses that do some element computations before output.
   To allow for a generic interface, the calling site needs to supply a
   class derived from SpecialFieldInterface that knows which function to call.

   \author kronbichler
   \date 04/14
   */
  virtual void WriteSpecialField (
      SpecialFieldInterface &special,
      PostResult& result,   ///< result group in the control file
      const ResultType  restype,
      const std::string &groupname,
      const std::vector<std::string> &fieldnames,
      const std::string &outinfo);

  void WriteDofResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_Vector> &data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf,
    const int from,
    const bool fillzeros) const;

  virtual void WriteNodalResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf) const;

  virtual void WriteElementResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf,
    const int from) const;

  /*!
   \brief write the geometry of one time step
   */
  void WriteGeo();

  /*!
   \brief write the solution data collected in the given vector
   */
  void WriteSolutionVector (const std::vector<double>& solution,
                            const int ncomponents,
                            const std::string &name,
                            std::ofstream& file) const;

  /*!
   \brief write prologue of the VTU file
   */
  void WriteVtuHeader ();

  /*!
   \brief write epilogue of the VTU file
   */
  void WriteVtuFooter (const std::string &filenamebase);


  /*!
   \brief writes debug output as long as this filter is incomplete (TODO MK: should go away at some point)
   */
  void CurrentlyNotImplemented () const
  {
    dserror("Functionality currently not implemented");
  }

  //! Total number of MPI processors (all processors output their row elements)
  int numproc_;

  //! Output stream for current time step
  std::ofstream currentout_;

  //! Output stream for master file of current time step (only proc 0)
  std::ofstream currentmasterout_;

  //! Time value for the current time step
  double time_;

  //! Value of the output step for searching the data base
  int timestep_;

  //! Current cycle step (e.g. in nonlinear iteration, not used yet)
  int cycle_;
};

#endif
