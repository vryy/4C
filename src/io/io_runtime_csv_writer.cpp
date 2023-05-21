/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief Write output in csv format at runtime in serial
\level 2
*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "io_runtime_csv_writer.H"

#include "utils_exceptions.H"
#include "lib_globalproblem.H"

#include "io.H"
#include "io_control.H"

#include <iostream>
#include <utility>

/*-----------------------------------------------------------------------------------------------*/
namespace IO
{
  //! Interface for the implementation of RuntimeCsvWriter based on an inheritance graph. The pure
  //! virtual interface class has two derived classes: One for proc 0 that does all the writing and
  //! one for all other procs that do nothing.
  class RuntimeCsvWriterImpl
  {
   public:
    virtual ~RuntimeCsvWriterImpl() = default;

    //! register data vector
    virtual void RegisterDataVector(
        const std::string& dataname, unsigned int numcomponents, int precision) = 0;

    //! reset current time and time step number
    virtual void ResetTimeAndTimeStep(double time, unsigned int timestep) = 0;

    //! append data vector
    virtual void AppendDataVector(
        const std::string& dataname, const std::vector<double>& datavalues) = 0;

    //! write csv file to filesystem
    virtual void WriteFile() = 0;
  };

  //! This class does the writing to the file. It is created on proc 0.
  class RuntimeCsvWriterProc0 : public RuntimeCsvWriterImpl
  {
   public:
    explicit RuntimeCsvWriterProc0(std::string outputname);

    void RegisterDataVector(
        const std::string& dataname, unsigned int numcomponents, int precision) override;

    void ResetTimeAndTimeStep(const double time, const unsigned int timestep) override
    {
      time_ = time;
      timestep_ = timestep;
    }

    void AppendDataVector(
        const std::string& dataname, const std::vector<double>& datavalues) override;

    void WriteFile() override;

   private:
    void WriteFileHeader();

    //! key: result name, entry: (data vector, precision)
    std::map<std::string, std::pair<std::vector<double>, int>> data_vectors_;

    //! full path to output file
    std::string fullpathoutputfile_;

    bool header_line_written_;

    //! output name
    std::string outputname_;

    //! current time
    double time_;

    //! current time step
    unsigned int timestep_;
  };

  //! This class does intentionally nothing. It is created on all other procs except for proc 0.
  class RuntimeCsvWriterOtherProcs : public RuntimeCsvWriterImpl
  {
   public:
    void RegisterDataVector(
        const std::string& dataname, unsigned int numcomponents, int precision) override
    {
    }

    void ResetTimeAndTimeStep(const double time, const unsigned int timestep) override {}

    void AppendDataVector(
        const std::string& dataname, const std::vector<double>& datavalues) override
    {
    }

    void WriteFile() override {}
  };

  void RuntimeCsvWriter::RegisterDataVector(
      const std::string& dataname, unsigned int numcomponents, int precision)
  {
    implementation_->RegisterDataVector(dataname, numcomponents, precision);
  }

  void RuntimeCsvWriter::ResetTimeAndTimeStep(double time, unsigned int timestep)
  {
    implementation_->ResetTimeAndTimeStep(time, timestep);
  }

  void RuntimeCsvWriter::AppendDataVector(
      const std::string& dataname, const std::vector<double>& datavalues)
  {
    implementation_->AppendDataVector(dataname, datavalues);
  }

  void RuntimeCsvWriter::WriteFile() { implementation_->WriteFile(); }

  RuntimeCsvWriter::RuntimeCsvWriter(int myrank, std::string outputname)
  {
    if (myrank == 0)
      implementation_ = std::make_unique<RuntimeCsvWriterProc0>(std::move(outputname));
    else
      implementation_ = std::make_unique<RuntimeCsvWriterOtherProcs>();
  }

  RuntimeCsvWriter::~RuntimeCsvWriter() = default;

  RuntimeCsvWriterProc0::RuntimeCsvWriterProc0(std::string outputname)
      : header_line_written_(false), outputname_(std::move(outputname)), time_(0.0), timestep_(-1)
  {
    // determine full path to output prefix
    const std::string fullpathoutputprefix =
        (DRT::Problem::Instance()->OutputControlFile()->FileName());

    // set full path to output file
    fullpathoutputfile_ = fullpathoutputprefix + "-" + outputname_ + ".csv";

    // clear content
    {
      std::ofstream outputfile(fullpathoutputfile_, std::ios_base::out | std::ios_base::trunc);
      outputfile.close();
    }

    // in case of restart copy content of restart file to output file prior to restart time step
    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // determine full path to restart prefix
      const std::string fullpathrestartprefix =
          DRT::Problem::Instance()->OutputControlFile()->RestartName();

      // set full path to restart file
      const std::string fullpathrestartfile = fullpathrestartprefix + "-" + outputname_ + ".csv";

      std::ifstream restartfile(fullpathrestartfile, std::ios_base::out);

      // check if file was found
      if (not restartfile)
        dserror("restart file '%s' could not be found", fullpathrestartfile.c_str());

      std::stringstream sectionpriorrestart;

      // loop over lines of restart file
      std::string line;
      while (std::getline(restartfile, line))
      {
        std::istringstream sline(line);

        // get first word of line
        std::string word;
        std::getline(sline, word, ',');

        // skip first line of restart file
        if (word == "step") continue;

        // get time step of current line
        const int timestep = std::stoi(word);

        // append current line
        if (timestep <= restart) sectionpriorrestart << line << "\n";
      }

      restartfile.close();

      // write to file
      std::ofstream outputfile(fullpathoutputfile_, std::ios_base::app);
      outputfile << sectionpriorrestart.str();
      outputfile.close();
      header_line_written_ = true;
    }
  }

  void RuntimeCsvWriterProc0::RegisterDataVector(
      const std::string& dataname, const unsigned int numcomponents, const int precision)
  {
    data_vectors_[dataname] = {std::vector<double>(numcomponents, 0.0), precision};
  }

  void RuntimeCsvWriterProc0::AppendDataVector(
      const std::string& dataname, const std::vector<double>& datavalues)
  {
    if (not data_vectors_.count(dataname))
      dserror("data vector '%s' not registered!", dataname.c_str());

    if ((data_vectors_[dataname].first).size() != datavalues.size())
      dserror("size of data vector '%s' changed!", dataname.c_str());

    data_vectors_[dataname].first = datavalues;
  }

  void RuntimeCsvWriterProc0::WriteFile()
  {
    if (data_vectors_.empty()) dserror("no data vectors registered!");

    if (!header_line_written_) WriteFileHeader();

    std::ofstream outputfile(fullpathoutputfile_, std::ios_base::app);
    outputfile << timestep_ << "," << time_;

    for (const auto& data : data_vectors_)
    {
      const std::vector<double>& data_vector = (data.second).first;
      const int precision = (data.second).second;

      outputfile << std::setprecision(precision) << std::scientific;

      for (const auto& data_vector_component : data_vector)
        outputfile << "," << data_vector_component;
    }

    outputfile << "\n";
    outputfile.close();
  }

  void RuntimeCsvWriterProc0::WriteFileHeader()
  {
    std::ofstream outputfile(fullpathoutputfile_, std::ios_base::app);
    outputfile << "step,time";

    for (const auto& data : data_vectors_)
    {
      const std::string& dataname = data.first;
      const unsigned int numcomponents = ((data.second).first).size();

      if (numcomponents == 1)
        outputfile << "," << dataname;
      else
        for (unsigned int i = 0; i < numcomponents; ++i) outputfile << "," << dataname << ":" << i;
    }

    outputfile << "\n";
    outputfile.close();

    header_line_written_ = true;
  }
}  // namespace IO
