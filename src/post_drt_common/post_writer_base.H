/*----------------------------------------------------------------------*/
/*! \file

\brief contains base class for a generic output filter (ensight and vtk are derived from this class)

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/
/*----------------------------------------------------------------------*/

#ifndef POST_DRT_WRITER_BASE_H_
#define POST_DRT_WRITER_BASE_H_

#include "post_filter_base.H"

#include <string>
#include <fstream>

class PostField;
class Epetra_MultiVector;


//! Special writer class that is used to invoke a particular method in the output writer
//! e.g.
struct SpecialFieldInterface
{
  virtual ~SpecialFieldInterface() {}

  virtual std::vector<int> NumDfMap() = 0;

  virtual void operator()(std::vector<Teuchos::RCP<std::ofstream>>& files, PostResult& result,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::vector<std::string>& names) = 0;
};


//! Base class for various output writers that use generic interfaces (Ensight, VTU)
class PostWriterBase
{
 public:
  //! constructor. initializes the writer to a certain field
  PostWriterBase(PostField* field, const std::string& filename);

  //! destructor
  virtual ~PostWriterBase() {}

  //! return the field specified at construction
  PostField* GetField()
  {
    dsassert(field_ != NULL, "No field has been set");
    return field_;
  }

  const PostField* GetField() const
  {
    dsassert(field_ != NULL, "No field has been set");
    return field_;
  }

  virtual void WriteFiles(PostFilterBase& filter) = 0;

  //! write the whole thing for changing geometries
  virtual void WriteFilesChangingGeom(PostFilterBase& filter) = 0;

  /*!
   \brief write all time steps of a result

   Write results. The results are taken from a reconstructed
   Epetra_Vector. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part.
   */
  virtual void WriteResult(
      const std::string groupname,  ///< name of the result group in the control file
      const std::string name,       ///< name of the result to be written
      const ResultType restype,     ///< type of the result to be written (nodal-/element-based)
      const int numdf,              ///< number of dofs per node to this result
      const int from = 0,           ///< start position of values in nodes
      const bool fillzeros = false  ///< zeros are filled to ensight file when no data is available
      ) = 0;

  /*!
   \brief write all time steps of a result in one time step

   Write results. The results are taken from a reconstructed
   Epetra_Vector. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part. Currently file continuation is not supported
   because Paraview is not able to load it due to some weird wild card
   issue.
   */
  virtual void WriteResultOneTimeStep(PostResult& result,  ///< result group in the control file
      const std::string groupname,  ///< name of the result group in the control file
      const std::string name,       ///< name of the result to be written
      const ResultType restype,     ///< type of the result to be written (nodal-/element-based)
      const int numdf,              ///< number of dofs per node to this result
      bool firststep,               ///< bool whether this is the first time step
      bool laststep,                ///< bool whether this is the last time step
      const int from = 0            ///< start position of values in nodes
      ) = 0;

  /*!
   \brief write a particular variable to file

   Write results. Some variables need interaction with the post filter,
   e.g. structural stresses that do some element computations before output.
   To allow for a generic interface, the calling site needs to supply a
   class derived from SpecialFieldInterface that knows which function to call.
   */
  virtual void WriteSpecialField(SpecialFieldInterface& special,
      PostResult& result,  ///< result group in the control file
      const ResultType restype, const std::string& groupname,
      const std::vector<std::string>& fieldnames, const std::string& outinfo) = 0;

  /*!
   \brief Write one step of a nodal result
   */
  virtual void WriteNodalResultStep(std::ofstream& file,
      const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf) = 0;

  /*!
   \brief Write one step of an element result
   */
  virtual void WriteElementResultStep(std::ofstream& file,
      const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from) = 0;

 protected:
  PostField* field_;
  std::string filename_;
  unsigned int myrank_;   ///< global processor id
  unsigned int numproc_;  ///< number of processors
};

#endif /* POST_DRT_WRITER_BASE_H_ */
