/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief Write output in csv format at runtime in serial
\level 2
*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "runtime_csv_writer.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include <iostream>
#include <fstream>
#include <sstream>

/*-----------------------------------------------------------------------------------------------*/

RuntimeCsvWriter::RuntimeCsvWriter(unsigned int myrank) : myrank_(myrank), time_(0.0), timestep_(-1)
{
  // empty constructor
}

void RuntimeCsvWriter::Init(const std::string& outputname)
{
  if (myrank_ != 0) return;

  // set output name
  outputname_ = outputname;

  // determine full path to output prefix
  const std::string fullpathoutputprefix(DRT::Problem::Instance()->OutputControlFile()->FileName());

  // set full path to output file
  fullpathoutputfile_ = fullpathoutputprefix + "-" + outputname_ + ".csv";

  // clear content
  std::ofstream outputfile(fullpathoutputfile_, std::ios_base::out | std::ios_base::trunc);
  outputfile.close();

  isinit_ = true;
}

void RuntimeCsvWriter::RegisterDataVector(
    const std::string& dataname, const unsigned int numcomponents, const unsigned int precision)
{
  if (myrank_ != 0) return;

  if (not isinit_) dserror("csv writer not initialized!");
  if (issetup_) dserror("register data vector prior to calling setup csv writer!");

  data_vectors_[dataname] = std::make_pair(std::vector<double>(numcomponents, 0.0), precision);
}

void RuntimeCsvWriter::Setup()
{
  if (myrank_ != 0) return;

  if (not isinit_) dserror("csv writer not initialized!");
  if (data_vectors_.size() == 0) dserror("no data vectors registered!");

  // write file header
  WriteFileHeader();

  // restart step
  const int restart = DRT::Problem::Instance()->Restart();

  // in case of restart copy content of restart file to output file prior to restart time step
  if (restart)
  {
    // determine full path to restart prefix
    const std::string fullpathrestartprefix =
        DRT::Problem::Instance()->OutputControlFile()->RestartName();

    // set full path to restart file
    const std::string fullpathrestartfile = fullpathrestartprefix + "-" + outputname_ + ".csv";

    std::ifstream restartfile(fullpathrestartfile, std::ios_base::out);

    // check if file was found
    if (not restartfile) dserror("restart file '%s' could not be found", fullpathrestartfile);

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
  }

  issetup_ = true;
}

void RuntimeCsvWriter::AppendDataVector(
    const std::string& dataname, const std::vector<double>& datavalues)
{
  if (myrank_ != 0) return;

  if (not isinit_) dserror("csv writer not initialized!");
  if (not issetup_) dserror("csv writer not setup!");

#ifdef DEBUG
  if (not data_vectors_.count(dataname))
    dserror("data vector '%s' not registered!", dataname.c_str());

  if ((data_vectors_[dataname].first).size() != datavalues.size())
    dserror("size of data vector '%s' changed!", dataname.c_str());
#endif

  data_vectors_[dataname].first = datavalues;
}

void RuntimeCsvWriter::WriteFile() const
{
#ifdef DEBUG
  if (myrank_ != 0)
  {
    if (isinit_) dserror("csv writer should only be initialized on proc 0!");
    if (issetup_) dserror("csv writer should only be setup on proc 0!");
  }
#endif

  if (myrank_ != 0) return;

  std::ofstream outputfile(fullpathoutputfile_, std::ios_base::app);
  outputfile << timestep_ << "," << time_;

  for (auto& data : data_vectors_)
  {
    const std::vector<double>& data_vector = (data.second).first;
    const unsigned int precision = (data.second).second;

    outputfile << std::setprecision(precision) << std::scientific;

    for (auto data_vector_component : data_vector) outputfile << "," << data_vector_component;
  }

  outputfile << "\n";
  outputfile.close();
}

void RuntimeCsvWriter::WriteFileHeader() const
{
  if (myrank_ != 0) return;

  std::ofstream outputfile(fullpathoutputfile_, std::ios_base::app);
  outputfile << "step,time";

  for (auto& data : data_vectors_)
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
}
