/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief Write output in csv format at runtime in serial
\level 2
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BACI_IO_RUNTIME_CSV_WRITER_HPP
#define BACI_IO_RUNTIME_CSV_WRITER_HPP

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include "baci_config.hpp"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

BACI_NAMESPACE_OPEN

namespace IO
{
  class OutputControl;
  class RuntimeCsvWriterImpl;

  /*!
   * \brief Write output in csv format at runtime in serial
   *
   * Output is only written on proc 0
   *
   * \note In the case of a restart, the content prior to the restart step is copied from the output
   * file the restart is read from.
   *
   */
  class RuntimeCsvWriter
  {
   public:
    RuntimeCsvWriter(int myrank, const IO::OutputControl& output_control, std::string outputname);

    ~RuntimeCsvWriter();

    //! Register name of column @p dataname for a future call to AppendDataVector(). The column can
    //! be split into @p numcomponents columns which are indicated by a post fixed identifier
    //! (":<number of component>"). The numerical precision is @p precision.
    void RegisterDataVector(const std::string& dataname, unsigned int numcomponents, int precision);

    //! set current time and time step number
    void ResetTimeAndTimeStep(double time, unsigned int timestep);

    //! append data @p datavalues to column @p dataname. In case of sub columns @p datavalues must
    //! have as many entries as the number of subcolumns
    void AppendDataVector(const std::string& dataname, const std::vector<double>& datavalues);

    //! write one line to csv file. Data must have been passed via AppendDataVector()
    void WriteCollectedDataToFile();

    //! write @p data to file at @p time and @p timestep
    void WriteDataToFile(double time, unsigned int timestep,
        const std::map<std::string, std::vector<double>>& data) const;

   private:
    std::unique_ptr<RuntimeCsvWriterImpl> implementation_;
  };
}  // namespace IO

BACI_NAMESPACE_CLOSE

#endif