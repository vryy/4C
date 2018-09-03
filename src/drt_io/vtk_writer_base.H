/*----------------------------------------------------------------------*/
/*!
\file vtk_writer_base.H

\brief VTK writer base class

\maintainer Martin Kronbichler

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef IO_VTK_WRITER_BASE_H
#define IO_VTK_WRITER_BASE_H


#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include <stdint.h>
#include <zlib.h>


#include "../drt_lib/drt_dserror.H"



namespace LIBB64
{
  /**
   * Do a base64 encoding of the given data.
   *
   * The function allocates memory as necessary and returns a pointer to
   * it. The calling function must release this memory again.
   */
  char* encode_block(const char* data, const int data_size);

  template <typename T>
  void writeCompressedBlock(const std::vector<T>& data, std::ostream& out)
  {
    uLongf compressed_data_length = compressBound(data.size() * sizeof(T));
    char* compressed_data = new char[compressed_data_length];
    int err = compress2((Bytef*)compressed_data, &compressed_data_length, (const Bytef*)&data[0],
        data.size() * sizeof(T), Z_BEST_COMPRESSION);
    if (err != Z_OK) dserror("zlib compression failed");

    // now encode the compression header
    const uint32_t compression_header[4] = {1, /* number of blocks */
        (uint32_t)(data.size() * sizeof(T)),   /* size of block */
        (uint32_t)(data.size() * sizeof(T)),   /* size of last block */
        (uint32_t)compressed_data_length};     /* list of compressed sizes of blocks */

    char* encoded_header =
        encode_block((char*)&compression_header[0], 4 * sizeof(compression_header[0]));
    out << encoded_header;
    delete[] encoded_header;

    // next do the compressed
    // data encoding in base64
    char* encoded_data = encode_block(compressed_data, compressed_data_length);
    delete[] compressed_data;

    out << encoded_data;
    out << std::endl;
    delete[] encoded_data;
  }

  /**
   \brief Helper function to determine output file string from time step number
   */
  std::string int2string(const unsigned int i, const unsigned int digits);

  /**
   \brief Helper function to determine number of digits for a given integer value
   */
  inline unsigned int ndigits(unsigned int number)
  {
    // start numbering from 0, so need count digits based on number one less
    if (number > 1) number -= 1;
    unsigned int digits = 0;
    while (number > 0)
    {
      digits++;
      number /= 10;
    }
    return digits;
  }

}  // namespace LIBB64


/*
 \brief Base class for VTK output generation

 \author kronbichler
 \date 03/14
*/
class VtkWriterBase
{
 public:
  //! constructor
  VtkWriterBase();

  //! destructor
  virtual ~VtkWriterBase() {}

  /** \brief initialize object with required input data
   *
   *  \author grill
   *  \date 03/17 */
  void Initialize(unsigned int myrank, unsigned int num_processors,
      unsigned int max_number_timesteps_to_be_written,
      const std::string& path_existing_working_directory,
      const std::string& name_new_vtk_subdirectory, const std::string& geometry_name,
      const std::string& restart_name, double restart_time, bool write_binary_output);

  /** \brief set class variable storing working directory and create it if not existing
   *
   *  \author grill
   *  \date 03/17 */
  void SetAndCreateVtkWorkingDirectory(const std::string& path_existing_working_directory,
      const std::string& name_new_vtk_subdirectory);

  /** \brief reset current simulation time and time step number
   *
   *  \author grill
   *  \date 03/17 */
  void ResetTimeAndTimeStep(double time, unsigned int timestepnumber);

  /** \brief reset geometry name
   *
   *  \author grill
   *  \date 04/17 */
  void ResetGeometryName(const std::string& geometryname);

  /** \brief initialize the required file streams for processor individual file and master file
   *
   *  \author grill
   *  \date 03/17 */
  void InitializeVtkFileStreamsForNewGeometryAndOrTimeStep();

  //! write prologue of all required vtk files
  void WriteVtkHeadersAndFieldData();

  //! write epilogue of all required vtk files
  void WriteVtkFooters();

  //! write a VTK collection file that summarizes paths to all written files (e.g. for all time
  //! steps) note: this only includes files written during the 'lifetime' of this writer object
  void WriteVtkCollectionFileForAllWrittenMasterFiles(const std::string& collectionfilename) const;


 protected:
  //! write a data vector as DataArray to corresponding vtk files
  // Todo template <typename T>
  void WriteDataArray(
      const std::vector<double>& data, const int num_components, const std::string& name);

  //! generate the part of the filename that expresses the processor ID
  const std::string& GetPartOfFileNameIndicatingProcessorId(unsigned int processor_id) const;


  //! Return the opening xml tag for this writer type
  virtual const std::string& WriterOpeningTag() const = 0;

  //! Return the parallel opening xml tag for this writer type
  virtual const std::string& WriterPOpeningTag() const = 0;

  //! Return a vector of parallel piece tags for each file
  virtual const std::vector<std::string>& WriterPPieceTags() const = 0;

  //! Return the parallel file suffix including the dot for this file type
  virtual const std::string& WriterPSuffix() const = 0;

  //! Return the string of this writer type
  virtual const std::string& WriterString() const = 0;

  //! Return the file suffix including the dot for this file type
  virtual const std::string& WriterSuffix() const = 0;


  //! throw error if file stream not ready to write into
  void ThrowErrorIfInvalidFileStream(const std::ofstream& filestream) const;

 private:
  //! write prologue of the VTK master file (handled by proc 0)
  void WriteVtkHeaderMasterFile(const std::string& byteorder);

  //! write prologue of the VTK file on this processor
  void WriteVtkHeaderThisProcessor(const std::string& byteorder);


  // if previously set, write time and cycle as 'field data'
  void WriteTimeOrCycleAsFieldData();


  //! write the required information about the DataArray to master file
  // Todo template <typename T>
  void WriteDataArrayMasterFile(const int num_components, const std::string& name);

  //! write the data array to this processor's file
  // Todo template <typename T>
  void WriteDataArrayThisProcessor(
      const std::vector<double>& data, const int num_components, const std::string& name);


  //! write epilogue of of the VTK master file (handled by proc 0)
  void WriteVtkFooterMasterFile();

  //! write epilogue of the VTK file on this processor
  void WriteVtkFooterThisProcessor();


  //! initialize the individual vtk file stream on each processor
  void InitializeVtkFileStreamThisProcessor();

  //! initialize the vtk 'master file' stream (handled by proc 0)
  void InitializeVtkMasterFileStream();

  //! append current master file name and time to collection file content
  void AppendMasterFileAndTimeToCollectionFileMidSectionContent(const std::string& master_file_name,
      const std::string& master_file_directory_name, double time);

  //! write a VTK collection file that summarizes paths to all the given files
  void WriteVtkCollectionFileForGivenListOfMasterFiles(const std::string& collectionfilename,
      const std::vector<std::pair<double, std::string>>& masterfiles_time_and_name) const;

  /** write a VTK collection file that summarizes paths to all the given files
   *
   *  \author eichinger
   *  \date 05/17 */
  void CreateRestartedInitialCollectionFileMidSection(
      const std::string& geometryname, const std::string& restartfilename, double restart_time);

  //! construct full path and name for VTK collection file stream from given file name
  std::string GetVtkCollectionFileFullPathAndName(const std::string& collectionfilename) const;

  //! write the header into a given VTK collection file stream
  void WriteHeaderIntoGivenVtkCollectionFileStream(std::ofstream& collectionfilestream) const;

  //! write the footer into a given VTK collection file stream
  void WriteFooterIntoGivenVtkCollectionFileStream(std::ofstream& collectionfilestream) const;

  //! determine the name of the subdirectory by extracting the part after the last '/'
  //! in the full path of the vtk working directory
  std::string DetermineVtkSubdirectoryNameFromFullVtkWorkingPath() const;

  //! write given master file and time value into given collection file stream
  void WriteMasterFileAndTimeValueIntoGivenVtkCollectionFileStream(
      std::ofstream& collectionfilestream, const std::string& master_file_name,
      const std::string& master_file_directory_name, double time) const;

 protected:
  //! subsequent phases in the process of writing data to VTK files:
  //  header -> point data -> cell data -> footer
  enum Phase
  {
    INIT,
    POINTS,
    CELLS,
    FINAL,
    VAGUE
  };

  // current phase
  Phase currentPhase_;


  //! Output stream for current processor-specific file
  std::ofstream currentout_;

  //! Output stream for current master file (only proc 0)
  std::ofstream currentmasterout_;

  //! stream for midsection of the vtk collection file ('.pvd')
  //! containing [time value and filename] of all yet written master files
  std::stringstream collection_file_midsection_cumulated_content_;


  //! full path of the working directory for all vtk files to be written into
  std::string working_directory_full_path_;

  //! Part of the current filename w/ timestep w/out processor id
  //  processor-specific files and master file share this part of the filename
  std::string filename_base_;


  //! number of digits (i.e. field width) which is reserved for time step number in file names
  unsigned int num_timestep_digits_;

  //! number of digits (i.e. field width) which is reserved for number of processors in file names
  unsigned int num_processor_digits_;


  //! name of the geometry
  std::string geometry_name_;

  //! Time value for the current time step
  double time_;

  //! output step
  int timestep_;

  //! restart step
  int restart_;

  //! Current cycle step (e.g. in nonlinear iteration, not used yet)
  int cycle_;


  //! toggle between ascii and binary output
  bool write_binary_output_;


  //! global processor id of this processor
  unsigned int myrank_;

  //! number of involved processors
  unsigned int numproc_;
};

#endif
