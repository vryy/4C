/*----------------------------------------------------------------------*/
/*! \file

 \brief contains base class for a generic output filter (ensight and vtk are derived from this
 class)


 \level 1
 */

#ifndef FOUR_C_POST_FILTER_BASE_HPP
#define FOUR_C_POST_FILTER_BASE_HPP

#include "baci_config.hpp"

#include "baci_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
class PostField;
class PostResult;
class PostWriterBase;


//! defines different result types (node-based, element-based, Gauss-Points?)
enum ResultType
{
  no_restype,    ///< unknown result type
  dofbased,      ///< result values based on dofrowmap
  nodebased,     ///< result values based on noderowmap
  elementbased,  ///< result values based on elementrowmap
  elementdof,    ///< result values based on elementdofrowmap (unused)
  max_restype    ///< end marker. must be the last entry
};


//! Base class for various field writers
class PostFilterBase
{
 public:
  /*
  //! Constructor. Initializes the writer to the current field
   */
  PostFilterBase(PostField* field, const std::string& name);

  //! destructor
  virtual ~PostFilterBase() = default;
  //! Writes the whole thing by invoking writer->WriteFiles()
  void WriteFiles();

  //! get the underlying writer object
  PostWriterBase& GetWriter()
  {
    FOUR_C_ASSERT(writer_ != Teuchos::null, "Not initialized");
    return *writer_;
  }

  //! look for problem dependent result entries and write them
  virtual void WriteAllResults(PostField* field) = 0;

  //! look for problem dependent result entries and write them for one time step
  virtual void WriteAllResultsOneTimeStep(PostResult& result, bool firststep, bool laststep)
  {
    FOUR_C_THROW("not implemented");
  }

  /// write all element based results
  /*!
    \note This method sees only those element results that are defined in
    the first result group.
   */
  void WriteElementResults(PostField* field) { WriteAnyResults(field, "element", elementbased); }

  /// write all node based results
  /*!
    \note This method sees only those node results that are defined in
    the first result group.
   */
  void WriteNodeResults(PostField* field) { WriteAnyResults(field, "node", nodebased); }

  /// write all dof based results
  /*!
    \note This method sees only those dof results that are defined in
    the first result group.
   */
  void WriteDofResults(PostField* field) { WriteAnyResults(field, "dof", dofbased); }

  /// write all results of the given type
  /*!
    \note This method sees only those results that are defined in the first
    result group (that is the first time step).
   */
  void WriteAnyResults(PostField* field, const char* type, const ResultType restype);

 protected:
  //! The actual writer object
  Teuchos::RCP<PostWriterBase> writer_;
};

FOUR_C_NAMESPACE_CLOSE

#endif
