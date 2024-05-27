/*----------------------------------------------------------------------*/
/*! \file

\brief VTK filter


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_POST_VTK_WRITER_HPP
#define FOUR_C_POST_VTK_WRITER_HPP


#include "4C_config.hpp"

#include "4C_io_vtk_writer_base.hpp"  // LIBB64
#include "4C_post_writer_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>
#include <zlib.h>

#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN
class PostField;
class PostResult;


namespace DRT
{
  class Discretization;
  class Node;
}  // namespace DRT


/*
 \brief Base class for VTK output generation

 \author kronbichler
 \date 03/14
*/
class PostVtkWriter : public PostWriterBase
{
 public:
  //! constructor. Initializes the writer to a certain field.
  PostVtkWriter(PostField* field, const std::string& name);

  //! write the whole thing
  void WriteFiles(PostFilterBase& filter) override;

 protected:
  //! Return the opening xml tag for this writer type
  virtual const std::string& writer_opening_tag() const = 0;

  //! Return the parallel opening xml tag for this writer type
  virtual const std::string& writer_p_opening_tag() const = 0;

  //! Return a vector of parallel piece tags for each file
  virtual const std::vector<std::string>& writer_p_piece_tags() const = 0;

  //! Give every writer a chance to do preparations before writing
  virtual void writer_prep_timestep() = 0;

  //! Return the parallel file suffix including the dot for this file type
  virtual const std::string& writer_p_suffix() const = 0;

  //! Return the string of this writer type
  virtual const std::string& writer_string() const = 0;

  //! Return the file suffix including the dot for this file type
  virtual const std::string& writer_suffix() const = 0;

  //! Write a single result step
  virtual void write_dof_result_step(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf, const int from,
      const bool fillzeros) = 0;

  //! Write a single result step
  void write_nodal_result_step(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf) override = 0;

  //! Write a single result step
  void write_element_result_step(std::ofstream& file, const Teuchos::RCP<Epetra_MultiVector>& data,
      std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
      const std::string& groupname, const std::string& name, const int numdf,
      const int from) override = 0;

  //! write the geometry of one time step
  virtual void write_geo() = 0;

  /*!
   \brief write one time step of a result

   The results are taken from a reconstructed
   Epetra_Vector. In many cases this vector will contain just one
   variable (displacements) and thus is easy to write as a whole. At
   other times, however, there is more than one result (velocity,
   pressure) and we want to write just one part of it. So we have to
   specify which part.
   */
  void WriteResult(const std::string groupname,  ///< name of the result group in the control file
      const std::string name,                    ///< name of the result to be written
      const ResultType restype,     ///< type of the result to be written (nodal-/element-based)
      const int numdf,              ///< number of dofs per node to this result
      const int from = 0,           ///< start position of values in nodes
      const bool fillzeros = false  ///< zeros are filled to vtk file when no data is available
      ) override;


  void write_result_one_time_step(PostResult& result,  ///< result group in the control file
      const std::string groupname,  ///< name of the result group in the control file
      const std::string name,       ///< name of the result to be written
      const ResultType restype,     ///< type of the result to be written (nodal-/element-based)
      const int numdf,              ///< number of dofs per node to this result
      bool firststep,               ///< bool whether this is the first time step
      bool laststep,                ///< bool whether this is the last time step
      const int from = 0            ///< start position of values in nodes
      ) override
  {
    FOUR_C_THROW("Not yet implemented");
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
  void WriteSpecialField(SpecialFieldInterface& special,
      PostResult& result,  ///< result group in the control file
      const ResultType restype, const std::string& groupname,
      const std::vector<std::string>& fieldnames, const std::string& outinfo) override;


  //! write the solution data collected in the given vector
  virtual void write_solution_vector(const std::vector<double>& solution, const int ncomponents,
      const std::string& name, std::ofstream& file) const;

  //! write prologue of the VTK file
  virtual void write_vtk_header();

  //! write epilogue of the VTK file
  virtual void write_vtk_footer();

  //! write a master file for VTK that collects links to all time steps
  virtual void write_vtk_master_file(const std::vector<std::pair<double, std::string>>& filenames,
      const std::string& dirname) const;

  //! writes debug output as long as this filter is incomplete (TODO MK: should go away at some
  //! point)
  virtual void currently_not_implemented() const
  {
    FOUR_C_THROW("Functionality currently not implemented");
  }

  enum Phase
  {
    INIT,
    POINTS,
    CELLS,
    FINAL
  };

  Phase currentPhase_;

  //! Output stream for current time step
  std::ofstream currentout_;

  //! Output stream for master file of current time step (only proc 0)
  std::ofstream currentmasterout_;

  //! Part of the filename w/ timestep w/out processor id
  std::string filenamebase_;

  unsigned int ntdigits_;

  unsigned int npdigits_;

  //! Time value for the current time step
  double time_;

  //! Value of the output step for searching the data base
  int timestep_;

  //! Current cycle step (e.g. in nonlinear iteration, not used yet)
  int cycle_;

  //! toggle between text and binary output
  bool write_binary_output_;
};

FOUR_C_NAMESPACE_CLOSE

#endif
