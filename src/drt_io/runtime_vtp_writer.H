/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output in vtk/vtp format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef IO_RUNTIME_VTP_WRITER_H
#define IO_RUNTIME_VTP_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include <Teuchos_RCP.hpp>

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
class VtpWriter;

namespace DRT
{
  class Discretization;
}

/*-----------------------------------------------------------------------------------------------*/
/* namespace */

/*!
 * \brief This object allows to write visualization output
 *        - in vtk/vtp format
 *        - at runtime
 *        - in parallel
 *        - binary-encoded
 *
 * \author grill
 * \date 04/17
 */
class RuntimeVtpWriter
{
 public:
  /// Constructor
  RuntimeVtpWriter();

  /// Destructor
  virtual ~RuntimeVtpWriter() {}

  /** \brief reset current simulation time and time step number
   *
   *  \author grill
   *  \date 04/17 */
  void Initialize(unsigned int myrank, unsigned int num_processors,
      unsigned int max_number_timesteps_to_be_written,
      const std::string& path_existing_output_directory, const std::string& simulation_name,
      const std::string& geometry_name, const std::string& restart_name, double restart_time,
      bool write_binary_output);

  /** \brief reset current simulation time and time step number and geometry name
   *
   *  \author grill
   *  \date 04/17 */
  void SetupForNewTimeStepAndGeometry(
      double time, unsigned int timestep, const std::string& geometryname);

  /** \brief reset geometry data for poly data set all at once
   *
   *  \author grill
   *  \date 04/17 */
  void ResetGeometry(const std::vector<double>& point_coordinates);

  /** \brief get mutable point coordinate vector for poly data set
   *
   *         this allows you to append geometry data in multiple steps
   *
   *  \author grill
   *  \date 04/17 */
  std::vector<double>& GetMutablePointCoordinateVector();


  /** \brief append visualization data vector with n components per point
   *
   *  \author grill
   *  \date 04/17 */
  // Todo template <typename T>, double, int,
  void AppendVisualizationPointDataVector(const std::vector<double>& datavalues,
      const unsigned int num_components_per_point, const std::string& dataname);


  /** \brief write all required VTK files to filesystem
   *
   *  \author grill
   *  \date 04/17 */
  void WriteFiles();

  /** \brief write a VTK collection file summarizing all written files of this object to filesystem
   *
   *  \author grill
   *  \date 04/17 */
  void WriteCollectionFileOfAllWrittenFiles(const std::string& collectionfilename);


 private:
  /** \brief write all stored data vectors to VTK files
   *
   *  \author grill
   *  \date 04/17 */
  void WriteAllPointDataVectorsToFiles();


 private:
  //! the actual vtp writer object
  Teuchos::RCP<VtpWriter> vtp_writer_;


  //! @name member variables required for the definition of a poly data set
  //! @{

  //! point coordinates with three components each
  std::vector<double> point_coordinates_;

  //! @}


  //! key: result name, entry: (solution data vector, num_components)
  std::map<std::string, std::pair<std::vector<double>, unsigned int>> point_data_vectors_;
};

#endif
