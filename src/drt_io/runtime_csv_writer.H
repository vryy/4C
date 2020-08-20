/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief Write output in csv format at runtime in serial
\level 2
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef IO_RUNTIME_CSV_WRITER_H
#define IO_RUNTIME_CSV_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include <string>
#include <vector>
#include <map>

/*!
 * \brief Write output in csv format at runtime in serial
 *
 * This object allows to write output of data in csv file format at runtime in serial, i.e., on
 * processor 0 only.
 *
 * \note In the case of a restart, the content prior to the restart step is copied from the output
 * file the restart is read from.
 *
 * \author Sebastian Fuchs \date 08/2020
 */
class RuntimeCsvWriter
{
 public:
  //! constructor
  explicit RuntimeCsvWriter(unsigned int myrank);

  //! destructor
  ~RuntimeCsvWriter() = default;

  //! init csv writer
  void Init(const std::string& outputname);

  //! register data vector
  void RegisterDataVector(
      const std::string& dataname, const unsigned int numcomponents, const unsigned int precision);

  //! setup csv writer
  void Setup();

  //! reset current time and time step number
  void ResetTimeAndTimeStep(double time, unsigned int timestep)
  {
    time_ = time;
    timestep_ = timestep;
  }

  //! append data vector
  void AppendDataVector(const std::string& dataname, const std::vector<double>& datavalues);

  //! write csv file to filesystem
  void WriteFile();

 private:
  //! write file header
  void WriteFileHeader();

  //! global processor id of this processor
  const unsigned int myrank_;

  //! current time
  double time_;

  //! current time step
  unsigned int timestep_;

  //! output name
  std::string outputname_;

  //! full path to output file
  std::string fullpathoutputfile_;

  //! key: result name, entry: (data vector, precision)
  std::map<std::string, std::pair<std::vector<double>, unsigned int>> data_vectors_;

  //! flag to check if csv writer is initialized
  bool isinit_ = false;

  //! flag to check if csv writer is setup
  bool issetup_ = false;
};

#endif
