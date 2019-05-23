/*!----------------------------------------------------------------------

\brief basic helper classes for image reconstruction

\level 2

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
*----------------------------------------------------------------------*/

#include "pat_utils.H"
#include "pat_imagereconstruction.H"

#include "Epetra_CrsMatrix.h"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_timestepping/timintmstep.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "acou_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_io/io_control.H"

bool TIMEWINDOW = false;  // true;
double TIMEWINDOW_BEGIN = 6.0;
double TIMEWINDOW_END = 22.5;

/*----------------------------------------------------------------------*/
ACOU::PATMonitorManager::PATMonitorManager(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<Teuchos::ParameterList> params, int illusetup)
{
  discret_ = discret;

  // determine dimensionality of the problem
  ndim_ = DRT::Problem::Instance()->NDim();

  // get tomograph type
  tomotype_ = DRT::INPUT::IntegralValue<INPAR::ACOU::TomographType>(
      params->sublist("PA IMAGE RECONSTRUCTION"), "TOMOGRAPHTYPE");

  // read monitor file
  std::string monitorfilename;
  if (illusetup == 0)
    monitorfilename = params->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("MONITORFILE");
  else if (illusetup == 1)
    monitorfilename = "partial_illu_secondsetup.monitor";  //"partillu_coarse_fw_2ndillu.monitor";
  else if (illusetup == 2)
    monitorfilename = "partial_illu_thirdsetup.monitor";  //"partillu_coarse_fw_3rdillu.monitor";
  else
    dserror("PATMonitorManager only for three monitor files implemented");

  if (monitorfilename == "none.monitor")
    dserror("No monitor file provided but required for optoacoustic image reconstruction");

  // insert path to monitor file if necessary
  if (monitorfilename[0] != '/')
  {
    std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
    std::string::size_type pos = filename.rfind('/');
    if (pos != std::string::npos)
    {
      std::string path = filename.substr(0, pos + 1);
      monitorfilename.insert(monitorfilename.begin(), path.begin(), path.end());
    }
  }

  // open monitor file and read it
  FILE* file = fopen(monitorfilename.c_str(), "rb");
  if (file == NULL) dserror("Could not open monitor file %s", monitorfilename.c_str());

  char buffer[150000];
  DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
  char* foundit = NULL;

  // read steps
  foundit = strstr(buffer, "steps");
  foundit += strlen("steps");
  unsigned int nsteps = strtol(foundit, &foundit, 10);

  dtacou_ = params->get<double>("TIMESTEP");
  if (params->get<double>("MAXTIME") / dtacou_ < params->get<int>("NUMSTEP"))
    nsteps_ = params->get<double>("MAXTIME") / dtacou_ + 1;
  else
    nsteps_ = params->get<int>("NUMSTEP") + 1;

  times_.resize(nsteps_);
  std::vector<double> timesteps(nsteps);

  // read mics
  foundit = strstr(buffer, "mics");
  foundit += strlen("mics");
  ngmics_ = nmics_ = strtol(foundit, &foundit, 10);
  positions_.resize(nmics_ * ndim_);
  cellidsset_ = false;

  for (unsigned int i = 0; i < nmics_; ++i)
  {
    DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
    foundit = buffer;
    for (unsigned int j = 0; j < ndim_;
         ++j)  // if it is a two dimensional problem, the third coordinate is discarded
      positions_[i * ndim_ + j] = strtod(foundit, &foundit);
  }

  // initialize other relevant quantities
  simulatedvalues_.clear();
  simulatedvalues_.resize(nmics_ * nsteps_);
  measurementvalues_.resize(nmics_ * nsteps_);

  // read comment lines
  foundit = buffer;
  DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
  while (strstr(buffer, "#"))
    DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);

  // read in the values for each node
  unsigned int count = 0;
  Epetra_SerialDenseVector mcurve(nmics_ * nsteps);
  for (unsigned int i = 0; i < nsteps; ++i)
  {
    // read the time step
    timesteps[i] = strtod(foundit, &foundit);
    for (unsigned int j = 0; j < nmics_; ++j) mcurve[count++] = strtod(foundit, &foundit);
    DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
    foundit = buffer;
  }
  if (count != nmics_ * nsteps) dserror("Number of measured pressure values wrong on input");

  // close input file
  fclose(file);

  // security check
  eps_ = dtacou_ / 1000.0;
  if (timesteps[0] > eps_) dserror("monitor file has to start at 0 but starts at %e", times_[0]);
  maxtime_ = std::min(params->get<int>("NUMSTEP") * dtacou_, params->get<double>("MAXTIME"));
  if (timesteps[nsteps - 1] < maxtime_ - eps_)
    dserror(
        "you only provide monitor values until %f but you want to simulate until %f. nsteps %f, "
        "nsteps x nmics %f",
        timesteps[nsteps - 1], maxtime_, nsteps, nmics_ * nsteps);

  // interpolate the read-in monitor values to the baci times
  for (unsigned int t = 0; t < nsteps_; ++t) times_[t] = double(t) * dtacou_;

  if (dtacou_ < timesteps[1] + eps_ &&
      dtacou_ > timesteps[1] - eps_)  // no interpolation necessary, just copy data
  {
    for (unsigned int i = 0; i < nmics_ * nsteps_; ++i) measurementvalues_[i] = mcurve[i];
  }
  else  // time steps not equally sized
  {
    // first line is simple because it is t=0
    for (unsigned int i = 0; i < nmics_; ++i) measurementvalues_[i] = mcurve[i];

    // all following times
    double dtmon = timesteps[1];
    for (unsigned int s = 1; s < nsteps_; ++s)
    {
      double sought_time = times_[s];

      int index = sought_time / dtmon;  // index of the smaller value
      double timequotient =
          (sought_time - timesteps[index]) / (timesteps[index + 1] - timesteps[index]);
      for (unsigned int i = 0; i < nmics_; ++i)
        measurementvalues_[s * nmics_ + i] =
            mcurve[index * nmics_ + i] +
            timequotient * (mcurve[(index + 1) * nmics_ + i] - mcurve[index * nmics_ + i]);
    }
  }

  // read quantities concerning impulse response
  // check if we need the impulse response and if so, transform it to the used time step
  double dtimpresp = params->sublist("PA IMAGE RECONSTRUCTION").get<double>("IMPULSERESPONSE_DT");
  if (dtimpresp != 0.0)
    consider_impulse_response_ = true;
  else
    consider_impulse_response_ = false;
  if (consider_impulse_response_)
  {
    // std::cout<<"consider_impulse_response_ "<<consider_impulse_response_<<std::endl;
    // get file name in which the impulse response is held
    std::string impulseresponsefilename =
        params->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("IMPULSERESPONSE");

    // check if file is given
    if (impulseresponsefilename == "none.impresp")
      dserror(
          "if you set IMPULSERESPONSE_DT != 0.0 you have to provide an impulse response in a file");

    // insert path to file if necessary
    if (impulseresponsefilename[0] != '/')
    {
      std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos != std::string::npos)
      {
        std::string path = filename.substr(0, pos + 1);
        impulseresponsefilename.insert(impulseresponsefilename.begin(), path.begin(), path.end());
      }
    }

    // open file
    FILE* file = fopen(impulseresponsefilename.c_str(), "rb");
    if (file == NULL)
      dserror("Could not open impulse response file %s", impulseresponsefilename.c_str());

    // read file
    char buffer[150000];
    DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, impulseresponsefilename);
    char* foundit = NULL;
    char* test = NULL;
    foundit = buffer;

    std::vector<double> impulseresponse;
    int num_imprespvals = 0;
    double norm = 0.0;
    do
    {
      impulseresponse.push_back(strtod(foundit, &foundit));
      test = fgets(buffer, 150000, file);
      foundit = buffer;
      num_imprespvals++;
      norm += impulseresponse[num_imprespvals - 1];
    } while (test != NULL);

    // in case time steps are the same just copy the impulse response
    if (dtimpresp == dtacou_)
    {
      impulse_response_.resize(num_imprespvals);
      for (int i = 0; i < num_imprespvals; ++i)
        impulse_response_[i] = impulseresponse[i] / std::abs(norm);
    }
    else  // otherwise interpolate the given impulse response to the dtacou_ timestep
    {
      int num_baciimprespvals = dtimpresp * num_imprespvals / dtacou_;
      impulse_response_.resize(num_baciimprespvals);
      double fac = double(num_imprespvals) / std::abs(norm) / double(num_baciimprespvals);

      for (int i = 0; i < num_baciimprespvals; ++i)
      {
        double actualt = i * dtacou_;           // we need values for this time point
        int impresindex = actualt / dtimpresp;  // corresponds to this index
        impulse_response_[i] =
            fac * (impulseresponse[impresindex] +
                      (impulseresponse[impresindex + 1] - impulseresponse[impresindex]) *
                          (actualt - impresindex * dtimpresp) / (dtimpresp));
      }
    }
  }  // read impulse response end
}

/*----------------------------------------------------------------------*/
void ACOU::PATMonitorManager::SetCellIds(
    std::vector<int> cellidsin, std::vector<bool> rowelesin, double fullfacemeasure)
{
  cellidsset_ = true;

  if (nmics_ != cellidsin.size())
    dserror("the number of cell ids should be the number of detectors");

  fullfacemeasure_ = fullfacemeasure;

  std::vector<bool> owned(nmics_, false);

  // reduce the number of stored values by checking the discretization on this processor and the
  // nodes in the monitor line condition this is necessary because otherwise the objective function
  // is evaluated wrongly if the simulation is run on several processors
  int countowned = 0;
  int countrowelesin = 0;
  for (unsigned int i = 0; i < cellidsin.size(); ++i)
  {
    if (cellidsin[i] >= 0)
    {
      owned[i] = true;
      countowned++;
    }
    if (rowelesin[i]) countrowelesin++;
  }

  // recalculate local quantities
  std::vector<double> measurementvaluesshort(countowned * nsteps_);
  std::vector<double> simulatedvaluesshort(countowned * nsteps_);
  std::vector<unsigned int> cellidsshort(countowned);
  std::vector<double> positionsshort(countowned * ndim_);
  std::vector<bool> cellrowflagshort(countowned);
  std::vector<double> facemeasuresshort(countowned);
  int count = 0;
  for (unsigned int i = 0; i < cellidsin.size(); ++i)
  {
    if (owned[i])
    {
      cellidsshort[count] = cellidsin[i];
      for (unsigned int s = 0; s < nsteps_; ++s)
      {
        measurementvaluesshort[s * countowned + count] = measurementvalues_[s * nmics_ + i];
        simulatedvaluesshort[s * countowned + count] = simulatedvalues_[s * nmics_ + i];
      }

      for (unsigned int d = 0; d < ndim_; ++d)
        positionsshort[count * ndim_ + d] = positions_[i * ndim_ + d];
      cellrowflagshort[count] = rowelesin[i];
      count++;
    }
  }

  // reset everything
  nmics_ = countowned;
  cellrowflag_.clear();
  cellrowflag_ = cellrowflagshort;
  cellids_.clear();
  cellids_ = cellidsshort;
  measurementvalues_.clear();
  measurementvalues_ = measurementvaluesshort;
  simulatedvalues_.clear();
  simulatedvalues_ = simulatedvaluesshort;
  positions_.clear();
  positions_ = positionsshort;

  sourcevalues_.clear();
  sourcevalues_ = measurementvaluesshort;  // in the first run we set it like this (then it works
                                           // for reduction simulations)

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PATMonitorManager::EvaluateMonitorValues(
    double time, const std::vector<double>& point, double& value)
{
  if (nmics_ == 0)  // in computations with many processors this function might be called even
                    // though it is not relevant
    return;

  time = maxtime_ - time;  // recalculate for adjoint evaluation

  // find closest microphone
  switch (tomotype_)
  {
    case INPAR::ACOU::pat_circle:  // every tomograph determines distances different, e.g. for this
                                   // one the z-component does not matter
    {
      if (ndim_ == 2)
      {
        double smallestdistance = std::numeric_limits<double>::max();
        double secondsmallestdistance = std::numeric_limits<double>::max();
        unsigned int smallestdistindex = -1;
        unsigned int secondsmallestdistindex = -1;
        for (unsigned int m = 0; m < nmics_; ++m)
        {
          double currentdistance = std::sqrt(
              (positions_[m * ndim_ + 0] - point[0]) * (positions_[m * ndim_ + 0] - point[0]) +
              (positions_[m * ndim_ + 1] - point[1]) * (positions_[m * ndim_ + 1] - point[1]));

          if (currentdistance < smallestdistance)
          {
            smallestdistance = currentdistance;
            smallestdistindex = m;
          }
        }
        for (unsigned int m = 0; m < nmics_; ++m)
        {
          double currentdistance = std::sqrt(
              (positions_[m * ndim_ + 0] - point[0]) * (positions_[m * ndim_ + 0] - point[0]) +
              (positions_[m * ndim_ + 1] - point[1]) * (positions_[m * ndim_ + 1] - point[1]));
          if (currentdistance < secondsmallestdistance &&
              currentdistance > smallestdistance - eps_ && m != smallestdistindex)
          {
            secondsmallestdistance = currentdistance;
            secondsmallestdistindex = m;
          }
        }

        // interpolate the values in time and space
        // first: evaluate in time
        int index = std::round(time / dtacou_);  // find time index

        double value_time_smallestdistance = sourcevalues_[smallestdistindex + index * nmics_];
        if (nmics_ > 1)
        {
          double value_time_secondsmallestdistance =
              sourcevalues_[secondsmallestdistindex + index * nmics_];

          value_time_smallestdistance *= ngmics_ / fullfacemeasure_;
          value_time_secondsmallestdistance *= ngmics_ / fullfacemeasure_;

          // second: interpolate in space (this is tomograph specific!! this one does not care about
          // z-component!)
          value = secondsmallestdistance / (smallestdistance + secondsmallestdistance) *
                      value_time_smallestdistance +
                  smallestdistance / (smallestdistance + secondsmallestdistance) *
                      value_time_secondsmallestdistance;
        }
        else
          value = value_time_smallestdistance * ngmics_ / fullfacemeasure_;
      }
      else if (ndim_ == 3)
      {
        std::vector<double> alldistances(nmics_);

        for (unsigned int m = 0; m < nmics_; ++m)
        {
          alldistances[m] = std::sqrt(
              (positions_[m * ndim_ + 0] - point[0]) * (positions_[m * ndim_ + 0] - point[0]) +
              (positions_[m * ndim_ + 1] - point[1]) * (positions_[m * ndim_ + 1] - point[1]) +
              (positions_[m * ndim_ + 2] - point[2]) * (positions_[m * ndim_ + 2] - point[2]));
        }

        std::vector<unsigned int> indices(4);
        double temp = std::numeric_limits<double>::max();
        for (unsigned int m = 0; m < nmics_; ++m)
        {
          double currentdistance = alldistances[m];
          if (currentdistance < temp)
          {
            temp = currentdistance;
            indices[0] = m;
          }
        }
        temp = std::numeric_limits<double>::max();
        for (unsigned int m = 0; m < nmics_; ++m)
        {
          double currentdistance = alldistances[m];
          if (currentdistance < temp && m != indices[0])
          {
            temp = currentdistance;
            indices[1] = m;
          }
        }
        temp = std::numeric_limits<double>::max();
        for (unsigned int m = 0; m < nmics_; ++m)
        {
          double currentdistance = alldistances[m];
          if (currentdistance < temp && m != indices[0] && m != indices[1])
          {
            temp = currentdistance;
            indices[2] = m;
          }
        }
        temp = std::numeric_limits<double>::max();
        for (unsigned int m = 0; m < nmics_; ++m)
        {
          double currentdistance = alldistances[m];
          if (currentdistance < temp && m != indices[0] && m != indices[1] && m != indices[2])
          {
            temp = currentdistance;
            indices[3] = m;
          }
        }
        std::vector<double> distances(4);
        for (unsigned int i = 0; i < 4; ++i) distances[i] = alldistances[indices[i]];
        double sumdistances = distances[0] + distances[1] + distances[2] + distances[3];

        // interpolate the values in time and space
        // first: evaluate in time
        int index = std::round(time / dtacou_);  // find time index

        std::vector<double> values_time(4);
        for (unsigned int i = 0; i < 4; ++i)
          values_time[i] = sourcevalues_[indices[i] + index * nmics_] * ngmics_ / fullfacemeasure_;
        if (nmics_ > 1)
        {
          // second: interpolate in space (this is tomograph specific!! this one does not care about
          // z-component!)
          value = (sumdistances - distances[0]) * values_time[0];
          for (unsigned int i = 1; i < 4; ++i)
            value += (sumdistances - distances[i]) * values_time[i];
          value /= (3.0 * sumdistances);
        }
        else
          value = sourcevalues_[indices[0] + index * nmics_] * ngmics_ / fullfacemeasure_;
      }
      else
        dserror("MonitorManager only implemented for 2D and 3D");
    }
    break;
    default:
      dserror("other tomograph types except circular not implemented");
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PATMonitorManager::StoreForwardValues(
    const double time, const std::vector<double> values)
{
  int index = std::round(time / dtacou_);

  for (unsigned int v = 0; v < values.size(); ++v) simulatedvalues_[nmics_ * index + v] = values[v];
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PATMonitorManager::ConvolveSignals()
{
  if (consider_impulse_response_)
  {
    // convolve forward signals with the impulse response
    std::vector<double> simulatedvalues_convolved(simulatedvalues_.size(), 0.0);
    for (unsigned m = 0; m < nmics_; ++m)
      for (unsigned s = 0; s < nsteps_; ++s)
        for (unsigned i = 0; i < impulse_response_.size() && i <= s; ++i)
          simulatedvalues_convolved[nmics_ * s + m] +=
              simulatedvalues_[nmics_ * (s - i) + m] * impulse_response_[i];

    // overwrite
    simulatedvalues_ = simulatedvalues_convolved;

    for (unsigned int i = 0; i < sourcevalues_.size(); ++i) sourcevalues_[i] = 0.0;

    // adjoint convolution for source values
    for (unsigned m = 0; m < nmics_; ++m)
      for (unsigned s = 0; s < nsteps_; ++s)
        for (unsigned i = 0; i < impulse_response_.size() && i + s < nsteps_; ++i)
          sourcevalues_[nmics_ * s + m] +=
              (simulatedvalues_[nmics_ * (s + i) + m] - measurementvalues_[nmics_ * (s + i) + m]) *
              impulse_response_[i];
  }
  else
  {
    // no convolution, only evaluation of source as difference
    if (TIMEWINDOW)
    {
      for (unsigned m = 0; m < nmics_; ++m)
        for (unsigned s = 0; s < nsteps_; ++s)
        {
          unsigned i = nmics_ * s + m;
          if (s * dtacou_ >= TIMEWINDOW_BEGIN && s * dtacou_ <= TIMEWINDOW_END)
            sourcevalues_[i] = simulatedvalues_[i] - measurementvalues_[i];
          else
          {
            sourcevalues_[i] = 0.0;
            simulatedvalues_[i] = 0.0;
            measurementvalues_[i] = 0.0;
          }
        }
    }
    else
    {
      for (unsigned i = 0; i < simulatedvalues_.size(); ++i)
      {
        // if((simulatedvalues_[i]-measurementvalues_[i])>5e-3)
        //  std::cout<<"i "<<i<<" simu "<<simulatedvalues_[i]<<" meas "<<measurementvalues_[i]<<"
        //  mic "<<i%nmics_<<" step "<<int(i/nmics_)<<" nmics "<<nmics_<<" at
        //  "<<positions_[(i%nmics_)*ndim_]<<" "<<positions_[(i%nmics_)*ndim_+1]<<"
        //  "<<positions_[(i%nmics_)*ndim_+2]<<std::endl;
        sourcevalues_[i] = simulatedvalues_[i] - measurementvalues_[i];
      }
    }
  }

  return;
}

void ACOU::PATMonitorManager::WriteMonitorFileInvana(std::string filename)
{
  int myrank = discret_->Comm().MyPID();
  std::vector<int> l_numrowdetecperproc(discret_->Comm().NumProc(), -1);
  std::vector<int> g_numrowdetecperproc(discret_->Comm().NumProc(), -1);

  // write header
  FILE* fp = NULL;
  if (myrank == 0)
  {
    fp = fopen(filename.c_str(), "w");
    if (fp == NULL) dserror("Couldn't open file.");
  }

  int l_numrowdetec = 0;
  for (unsigned int i = 0; i < cellrowflag_.size(); ++i)
    if (cellrowflag_[i]) l_numrowdetec++;

  l_numrowdetecperproc[myrank] = l_numrowdetec;
  discret_->Comm().MaxAll(
      &l_numrowdetecperproc[0], &g_numrowdetecperproc[0], l_numrowdetecperproc.size());

  int g_numrowdetec = g_numrowdetecperproc[0];
  for (int p = 1; p < discret_->Comm().NumProc(); ++p) g_numrowdetec += g_numrowdetecperproc[p];

  if (myrank == 0)
  {
    fprintf(fp, "steps %d ", nsteps_);
    fprintf(fp, "mics %d\n", g_numrowdetec);

    for (unsigned int i = 0; i < nmics_; ++i)
      if (cellrowflag_[i])
      {
        if (ndim_ == 2)
          fprintf(fp, "%e %e %e\n", positions_[i * ndim_], positions_[i * ndim_ + 1], 0.0);
        else if (ndim_ == 3)
          fprintf(fp, "%e %e %e\n", positions_[i * ndim_], positions_[i * ndim_ + 1],
              positions_[i * ndim_ + 2]);
      }
  }

  for (int p = 1; p < discret_->Comm().NumProc(); ++p)
  {
    std::vector<double> positionstemp(g_numrowdetec * ndim_);
    int count = 0;
    if (p == myrank)
    {
      for (unsigned int i = 0; i < cellrowflag_.size(); ++i)
      {
        if (cellrowflag_[i])
        {
          for (unsigned int d = 0; d < ndim_; ++d)
            positionstemp[count * ndim_ + d] = positions_[i * ndim_ + d];
          count++;
        }
      }
    }

    discret_->Comm().Broadcast(&positionstemp[0], g_numrowdetec * ndim_, p);
    if (myrank == 0)
    {
      for (int i = 0; i < g_numrowdetecperproc[0]; ++i)
      {
        if (ndim_ == 2)
          fprintf(fp, "%e %e %e\n", positionstemp[i * ndim_], positionstemp[i * ndim_ + 1], 0.0);
        else if (ndim_ == 3)
          fprintf(fp, "%e %e %e\n", positionstemp[i * ndim_], positionstemp[i * ndim_ + 1],
              positionstemp[i * ndim_ + 2]);
      }
    }
  }

  if (myrank == 0)
  {
    fprintf(fp, "#\n");
    fprintf(fp, "#\n");
    fprintf(fp, "#\n");
  }

  // get all measurement values from the other procs and then write all at once
  std::vector<double> l_allmeasurementsonproczero(g_numrowdetec * nsteps_, 0.0);
  if (myrank == 0)
  {
    int count = 0;
    for (unsigned int i = 0; i < nmics_; ++i)
    {
      if (cellrowflag_[i])
      {
        for (unsigned int t = 0; t < nsteps_; ++t)
          l_allmeasurementsonproczero[count + g_numrowdetec * t] = simulatedvalues_[i + nmics_ * t];
        count++;
      }
    }
  }

  for (int p = 1; p < discret_->Comm().NumProc(); ++p)
  {
    int offset = g_numrowdetecperproc[0];
    if (p == myrank)
    {
      for (int i = 1; i < p; ++i) offset += g_numrowdetecperproc[p];

      int count = 0;
      for (unsigned int i = 0; i < nmics_; ++i)
      {
        if (cellrowflag_[i])
        {
          for (unsigned int t = 0; t < nsteps_; ++t)
            l_allmeasurementsonproczero[offset + count + g_numrowdetec * t] =
                simulatedvalues_[i + nmics_ * t];
          count++;
        }
      }
    }
  }
  std::vector<double> g_allmeasurementsonproczero(g_numrowdetec * nsteps_, 0.0);
  discret_->Comm().SumAll(
      &l_allmeasurementsonproczero[0], &g_allmeasurementsonproczero[0], g_numrowdetec * nsteps_);
  if (myrank == 0)
  {
    for (unsigned int t = 0; t < nsteps_; ++t)
    {
      fprintf(fp, "%e ", times_[t]);
      for (int d = 0; d < g_numrowdetec; ++d)
        fprintf(fp, "%e ", g_allmeasurementsonproczero[t * g_numrowdetec + d]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PATMonitorManager::EvaluateError()
{
  double error = 0.0;
  for (unsigned int m = 0; m < sourcevalues_.size(); ++m)
    if (cellrowflag_[int(m % nmics_)]) error += sourcevalues_[m] * sourcevalues_[m];

  //  for(int t=0; t<times_.size(); ++t)
  //  {
  //    std::cout<<times_[t]<<" ";
  //      for(int d=0; d<nmics_; ++d)
  //        std::cout<<simulatedvalues_[t*nmics_+d]<<" ";
  //      std::cout<<std::endl;
  //
  //  }
  //
  //
  //  for(int t=0; t<times_.size(); ++t)
  //  {
  //    std::cout<<times_[t]<<" ";
  //      for(int d=0; d<nmics_; ++d)
  //        std::cout<<measurementvalues_[t*nmics_+d]<<" ";
  //      std::cout<<std::endl;
  //
  //  }

  double gerror = 0.0;

  discret_->Comm().SumAll(&error, &gerror, 1);

  return 0.5 * gerror;
}

/*----------------------------------------------------------------------*/
ACOU::PATSearchDirection::PATSearchDirection(INPAR::ACOU::OptimizationType opti) { opti_ = opti; }

/*----------------------------------------------------------------------*/
void ACOU::PATSearchDirection::Setup(const Epetra_Map* map, const Epetra_Map* uniquemap)
{
  if (opti_ == INPAR::ACOU::inv_lbfgs)
  {
    actsize_ = 0;
    ssize_ = 10;

    sstore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, map, true));
    ystore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, map, true));

    map_ = Teuchos::rcp(new Epetra_Map(*map));
    uniquemap_ = Teuchos::rcp(new Epetra_Map(*uniquemap));

    oldparams_ = Teuchos::rcp(new Epetra_Vector(*map_, false));
    oldgradient_ = Teuchos::rcp(new Epetra_Vector(*map_, false));
  }
  // nothing to do for steepest descent

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ACOU::PATSearchDirection::ComputeDirection(
    Teuchos::RCP<Epetra_Vector> gradient, Teuchos::RCP<Epetra_Vector> params, int iter)
{
  Teuchos::RCP<Epetra_Vector> direction = Teuchos::rcp(new Epetra_Vector(gradient->Map(), false));
  if (opti_ == INPAR::ACOU::inv_lbfgs)
  {
    if (iter == 0)
    {
      oldgradient_->Update(1.0, *gradient, 0.0);
      oldparams_->Update(1.0, *params, 0.0);
      direction->Update(-1.0, *gradient, 0.0);
    }
    else
    {
      // store vectors
      if (iter <= ssize_)
      {
        actsize_ += 1;

        sstore_->Resize(-actsize_ + 1, 0, map_.get(), false);
        ystore_->Resize(-actsize_ + 1, 0, map_.get(), false);

        Epetra_Vector s(*map_, false);
        s.Update(1.0, *params, -1.0, *oldparams_, 0.0);
        sstore_->UpdateSteps(s);

        s.Update(1.0, *gradient, -1.0, *oldgradient_, false);
        ystore_->UpdateSteps(s);
      }

      // compute direction
      direction->Update(1.0, *gradient, 0.0);
      std::vector<double> alpha;

      // loop steps
      for (int i = 0; i > -actsize_; i--)
      {
        double aa = 0.0;
        double bb = 0.0;

        ((*ystore_)(i))->Dot(*(*sstore_)(i), &aa);
        ((*sstore_)(i))->Dot((*direction), &bb);

        alpha.push_back(1 / aa * bb);

        direction->Update(-1.0 * alpha.back(), *(*ystore_)(i), 1.0);
      }

      for (int i = -actsize_ + 1; i <= 0; i++)
      {
        double aa = 0.0;
        double bb = 0.0;
        double beta = 0.0;

        ((*ystore_)(i))->Dot(*(*sstore_)(i), &aa);
        ((*ystore_)(i))->Dot((*direction), &bb);

        beta = 1 / aa * bb;
        double alphac = alpha.back();
        alpha.pop_back();

        direction->Update(alphac - beta, *(*sstore_)(i), 1.0);
      }

      // minimization not maximization
      direction->Scale(-1.0);
    }
  }
  else
  {
    if (0)
      direction->Update(-1000.0, *gradient, 0.0);
    else
    {
      double maxval = 0.0;
      gradient->MaxValue(&maxval);
      double minval = 0.0;
      gradient->MinValue(&minval);

      if (abs(minval) > maxval && (minval > 1.0e-10 || minval < -1.0e-10))
        direction->Update(-0.1 / abs(minval), *gradient, 0.0);
      else if (abs(maxval) > 1.0e-10)
        direction->Update(-0.1 / abs(maxval), *gradient, 0.0);
      else
      {
        double mean = 0.0;
        gradient->MeanValue(&mean);
        direction->Update(-0.5 / abs(mean), *gradient, 0.0);
      }
    }
  }
  return direction;
}

/*----------------------------------------------------------------------*/
ACOU::PATLineSearch::PATLineSearch(Teuchos::RCP<PatImageReconstruction> imagereconstruction)
{
  itermax_ = imagereconstruction->acouparams_->sublist("PA IMAGE RECONSTRUCTION")
                 .get<int>("INV_LS_MAX_RUN");

  c1_ = 1.0e-12;
  c2_ = 0.9;

  alpha_max_ = 10.0;  // initial step length ususally allows a change of 0.1

  imagereconstruction_ = imagereconstruction;
  myrank_ = imagereconstruction_->myrank_;
}

/*----------------------------------------------------------------------*/
void ACOU::PATLineSearch::Init(double J0, Teuchos::RCP<Epetra_Vector> gradient,
    Teuchos::RCP<Epetra_Vector> direction, Teuchos::RCP<Epetra_Vector> state,
    const Epetra_Map* uniquemap)
{
  // initialize objective function value
  J_0_ = J0;
  J_i_ = J0;
  J_im1_ = J0;

  // initialize vectorial quantities
  dir_ = direction;
  step_ = Teuchos::rcp(new Epetra_Vector(dir_->Map()));
  state_ = Teuchos::rcp(new Epetra_Vector(*state));

  // initialize unique map
  uniquemap_ = Teuchos::rcp(new Epetra_Map(*uniquemap));

  // initialize norm of derivative of line search function
  imagereconstruction_->CalculateGradDirNorm(*dir_, *uniquemap_, &normgradphi_0_);
  normgradphi_i_ = 0.0;

  // set step lengths
  alpha_i_ = 1.0;
  alpha_im1_ = 0.0;
  alpha_x_ = 0.0;

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PATLineSearch::Run()
{
  if (!myrank_)
    std::cout << "*************** RUN LINE SEARCH: J_0 " << J_0_ << " normgradphi_0 "
              << normgradphi_0_ << std::endl;

  for (int i = 0; i < itermax_; ++i)
  {
    if (!myrank_)
      std::cout << "*************** line search iteration " << i << " of maximal " << itermax_
                << " line search iterations, alpha_i_ " << alpha_i_ << ", J_0_ " << J_0_
                << ", J_i_ " << J_i_ << ", normgradphi_0_ " << normgradphi_0_ << ", normgradphi_i_ "
                << normgradphi_i_ << std::endl;

    // update parameters
    step_->Update(1.0, *state_, 0.0);
    step_->Update(alpha_i_, *dir_, 1.0);

    double maxval = -1.0;
    step_->MaxValue(&maxval);
    while (maxval < 0.0)
    {
      alpha_i_ /= 2.0;
      step_->Update(1.0, *state_, 0.0);
      step_->Update(alpha_i_, *dir_, 1.0);
      step_->MaxValue(&maxval);
      if (!myrank_)
        std::cout << "warning: had to reduce line search step length to " << alpha_i_
                  << ", all values were negative" << std::endl;
    }

    imagereconstruction_->ReplaceParams(step_);

    // solve forward problem
    imagereconstruction_->SolveStandardScatra();
    imagereconstruction_->SolveStandardAcou();

    // evaluate objective function
    J_i_ = imagereconstruction_->EvalulateObjectiveFunction();

    // check first condition
    if (J_i_ > J_0_ + c1_ * alpha_i_ * normgradphi_0_ || (J_i_ >= J_im1_ && i > 0))
    {
      if (!myrank_)
        std::cout << "*************** line search condition 1 met, J_i " << J_i_ << ", J_0 " << J_0_
                  << std::endl;
      alpha_x_ = Zoom(alpha_im1_, alpha_i_, J_im1_, false);
      break;
    }
    else if (!myrank_)
      std::cout << "*************** line search condition 1 NOT met, J_i " << J_i_ << ", J_0 "
                << J_0_ << std::endl;

    // solve adjoint problem
    imagereconstruction_->SolveAdjointAcou();
    imagereconstruction_->SolveAdjointScatra();

    // calculate gradient
    imagereconstruction_->EvaluateGradient();
    imagereconstruction_->CalculateGradDirNorm(*dir_, *uniquemap_, &normgradphi_i_);

    // check second condition
    if (std::abs(normgradphi_i_) <= -c2_ * normgradphi_0_)
    {
      alpha_x_ = alpha_i_;
      if (!myrank_)
        std::cout << "*************** line search condition 2 met, |\\/phi|_i " << normgradphi_i_
                  << ", |\\/phi|_0 " << normgradphi_0_ << std::endl;
      break;
    }
    else if (!myrank_)
      std::cout << "*************** line search condition 2 NOT met, |\\/phi|_i " << normgradphi_i_
                << ", |\\/phi|_0 " << normgradphi_0_ << std::endl;

    // check third condition
    if (normgradphi_i_ >= 0)
    {
      if (!myrank_) std::cout << "*************** line search condition 3 met" << std::endl;
      alpha_x_ = Zoom(alpha_i_, alpha_im1_, J_im1_, true);
      break;
    }
    else if (!myrank_)
      std::cout << "*************** line search condition 3 not met" << std::endl;

    // update alphas
    alpha_im1_ = alpha_i_;
    alpha_i_ *= 2.0;  // = PredictStepLength();

    if (alpha_i_ > alpha_max_ && J_0_ > J_i_)
    {
      alpha_x_ = alpha_im1_;
      break;
    }
    else if (alpha_i_ > alpha_max_ && J_0_ <= J_i_)
    {
      alpha_x_ = 0.0;
      break;
    }
  }

  if (J_i_ < J_0_)  // if(alpha_x_ != 0.0)
  {
    if (!myrank_) std::cout << "*************** line search succeeded" << std::endl;
    return true;
  }
  else
  {
    if (!myrank_) std::cout << "*************** line search failed" << std::endl;

    // replace params by original params
    imagereconstruction_->ReplaceParams(state_);

    return false;
  }
}

/*----------------------------------------------------------------------*/
double ACOU::PATLineSearch::PredictStepLength()
{
  // step length must become longer when we are here!
  double alpha = 0.0;

  // calculate the coefficients of the quartic polynomial
  // double d = J_0_;
  double c = normgradphi_0_;
  double b = -1.0 / alpha_i_ / alpha_i_ *
             (alpha_i_ * normgradphi_i_ + 2.0 * normgradphi_0_ * alpha_i_ - 3.0 * J_i_ + 3 * J_0_);
  double a = 1.0 / alpha_i_ / alpha_i_ / alpha_i_ *
             (J_i_ - J_0_ - normgradphi_0_ * alpha_i_ - b * alpha_i_ * alpha_i_);

  // calculate the minima of the quartic polynomial
  double radi = b * b / 9.0 / a / a - c / 3.0 / a;
  if (radi > 0.0)
  {
    double alpha1 = -b / 3.0 / a + sqrt(radi);
    double alpha2 = -b / 3.0 / a - sqrt(radi);
    // check if the results suit us
    if (alpha1 > alpha_i_ && alpha1 < 10.0 * alpha_i_)
      alpha = alpha1;
    else if (alpha2 > alpha_i_ && alpha2 < 10.0 * alpha_i_)
      alpha = alpha2;
    else if (alpha1 > 10.0 * alpha_i_ && alpha2 > 10.0 * alpha_i_)
      alpha = 10.0 * alpha_i_;
    else if (alpha1 < alpha_i_ && alpha2 < alpha_i_)
      alpha = 2.0 * alpha_i_;
    else
      alpha = 2.0 * alpha_i_;
  }
  else
  {
    // quadratic interpolation
    alpha = -normgradphi_0_ / 2.0 / (J_i_ - J_0_ - normgradphi_0_ * alpha_i_) * alpha_i_;
    if (alpha < alpha_i_)
      alpha = 2.0 * alpha_i_;
    else if (alpha > 10.0 * alpha_i_)
      alpha = 10.0 * alpha_i_;
  }
  return alpha;
}

/*----------------------------------------------------------------------*/
double ACOU::PATLineSearch::Zoom(double alpha_lo, double alpha_hi, double J_alpha_lo, bool turn)
{
  double alpha_j = 0.0;
  double J_j = 0.0;
  double normgradphi_j = 0.0;

  // zoom iterations
  for (int j = 0; j < itermax_; ++j)
  {
    // update alpha
    if (turn)
      alpha_j = alpha_hi + (alpha_lo - alpha_hi) / 3.0;
    else
      alpha_j = (alpha_lo + alpha_hi) / 2.0;

    // output for user
    if (!myrank_)
      std::cout << "*************** zoom iteration " << j << ": alpha_lo " << alpha_lo
                << " alpha_hi " << alpha_hi << " alpha_j " << alpha_j << std::endl;

    // update parameters
    step_->Update(1.0, *state_, 0.0);
    step_->Update(alpha_j, *dir_, 1.0);
    imagereconstruction_->ReplaceParams(step_);

    // solve forward problem
    imagereconstruction_->SolveStandardScatra();
    imagereconstruction_->SolveStandardAcou();

    // evaluate objective function
    J_j = imagereconstruction_->EvalulateObjectiveFunction();

    // output for user
    if (!myrank_)
      std::cout << "J_j " << J_j << " J_0_ " << J_0_ << " J_0_+c1_... "
                << J_0_ + c1_ * alpha_j * normgradphi_0_ << " J_alpha_lo " << J_alpha_lo
                << std::endl;

    if (J_j > J_0_ + c1_ * alpha_j * normgradphi_0_ || J_j >= J_alpha_lo)
    {
      if (!myrank_)
        std::cout << "*************** zoom condition 1 met"
                  << ", J_j " << J_j << ", J_0 " << J_0_ << std::endl;
      alpha_hi = alpha_j;
    }
    else
    {
      // solve adjoint problem
      imagereconstruction_->SolveAdjointAcou();
      imagereconstruction_->SolveAdjointScatra();

      // calculate gradient
      imagereconstruction_->EvaluateGradient();
      imagereconstruction_->CalculateGradDirNorm(*dir_, *uniquemap_, &normgradphi_j);

      // check second condition
      if (std::abs(normgradphi_j) <= -c2_ * normgradphi_0_)
      {
        if (!myrank_)
          std::cout << "*************** zoom condition 2 met, |\\/phi|_j " << normgradphi_j
                    << ", |\\/phi|_0 " << normgradphi_0_ << std::endl;
        J_i_ = J_j;
        return alpha_j;
      }
      else
      {
        if (!myrank_)
          std::cout << "*************** zoom condition 2 NOT met, |\\/phi|_j " << normgradphi_j
                    << ", |\\/phi|_0 " << normgradphi_0_ << std::endl;
      }

      // check third condition
      if (normgradphi_j * (alpha_hi - alpha_lo) >= 0)
      {
        if (!myrank_) std::cout << "*************** zoom condition 3 met" << std::endl;
        alpha_hi = alpha_lo;
      }
      else if (!myrank_)
        std::cout << "*************** zoom condition 3 NOT met" << std::endl;

      alpha_lo = alpha_j;
    }
  }

  if (J_j < J_0_)  // line search does not fulfill both wolfe conditions but i would not call it
                   // fail anyhow
  {
    J_i_ = J_j;
    return alpha_j;
  }
  else
    return 0.0;
}


/*----------------------------------------------------------------------*/
ACOU::PATRegula::PATRegula(INPAR::ACOU::RegulaType regulatype, double tikhweight, double tvdweight,
    double tvdeps, Teuchos::RCP<DRT::Discretization> discret)
    : type_(regulatype),
      tikhweight_(tikhweight),
      tvdweight_(tvdweight),
      tvdeps_(tvdeps),
      tvdmatrix_(Teuchos::null),
      discret_(discret)
{
  if (type_ == INPAR::ACOU::pat_regula_tvd || type_ == INPAR::ACOU::pat_regula_tikhtvd)
  {
    // create matrix
    tvdmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *discret_->ElementRowMap(), 6, false));

    // fill matrix
    FillTvdMatrix();
  }
}

/*----------------------------------------------------------------------*/
void ACOU::PATRegula::Evaluate(Teuchos::RCP<Epetra_Vector> params, double* val)
{
  if (type_ == INPAR::ACOU::pat_regula_tikh || type_ == INPAR::ACOU::pat_regula_tikhtvd)
  {
    double norm2 = 0.0;
    params->Norm2(&norm2);
    *val += 0.5 * tikhweight_ * norm2 * norm2;
  }
  if (type_ == INPAR::ACOU::pat_regula_tvd || type_ == INPAR::ACOU::pat_regula_tikhtvd)
  {
    // export parameters to column map
    Epetra_Vector paramscol(tvdmatrix_->ColMap(), false);
    LINALG::Export(*params, paramscol);

    double functvalue = 0.0;
    for (int i = 0; i < params->MyLength(); i++)
    {
      // get weights of neighbouring parameters
      int lenindices = tvdmatrix_->NumMyEntries(i);
      int numindex;
      std::vector<int> indices(lenindices, 0);
      std::vector<double> weights(lenindices, 0);
      tvdmatrix_->ExtractMyRowCopy(i, lenindices, numindex, &weights[0], &indices[0]);

      double rowsum = 0.0;
      double rowval = (*params)[i];
      for (int j = 0; j < lenindices; j++)
      {
        if (indices[j] != i)  // skip substracting from itself
        {
          double colval = (paramscol)[indices[j]];
          rowsum += weights[j] * (colval - rowval) * (colval - rowval);
        }
      }
      // sum over all the rows
      functvalue += sqrt(rowsum + tvdeps_);
    }

    // collect contributions from all procs
    double gfunctvalue = 0.0;
    discret_->Comm().SumAll(&functvalue, &gfunctvalue, 1);

    *val += tvdweight_ * gfunctvalue;
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PATRegula::EvaluateGradient(
    Teuchos::RCP<Epetra_Vector> params, Teuchos::RCP<Epetra_Vector> gradient)
{
  if (type_ == INPAR::ACOU::pat_regula_tikh || type_ == INPAR::ACOU::pat_regula_tikhtvd)
  {
    gradient->Update(1.0, *params, 1.0);
  }
  if (type_ == INPAR::ACOU::pat_regula_tvd || type_ == INPAR::ACOU::pat_regula_tikhtvd)
  {
    // export parameters to column map
    Epetra_Vector paramscol(tvdmatrix_->ColMap(), false);
    LINALG::Export(*params, paramscol);

    Epetra_Vector gradientcol(tvdmatrix_->ColMap(), true);
    for (int i = 0; i < params->MyLength(); i++)
    {
      // get weights of neighbouring parameters
      int lenindices = tvdmatrix_->NumMyEntries(i);
      int numindex;
      std::vector<int> indices(lenindices, 0);
      std::vector<double> weights(lenindices, 0);
      tvdmatrix_->ExtractMyRowCopy(i, lenindices, numindex, &weights[0], &indices[0]);

      // value of the parameter in this row
      double rowval = (*params)[i];

      // denominator
      double denom = 0.0;
      for (int j = 0; j < lenindices; j++)
      {
        if (indices[j] != i)  // skip substracting from itself
        {
          double colval = (paramscol)[indices[j]];
          denom += weights[j] * (colval - rowval) * (colval - rowval);
        }
      }
      denom = sqrt(denom + tvdeps_);

      // nominator i
      double nomi = 0.0;
      for (int j = 0; j < lenindices; j++)
      {
        if (indices[j] != i)  // skip substracting from itself
        {
          double colval = (paramscol)[indices[j]];
          nomi += weights[j] * (colval - rowval);
        }
      }
      nomi = nomi * (-1);

      // insert contributions into dof i of the gradient
      gradientcol.SumIntoMyValue(i, 0, nomi / denom);

      // nominator j
      double nomj = 0.0;
      for (int j = 0; j < lenindices; j++)
      {
        if (indices[j] != i)  // skip substracting from itself
        {
          double colval = (paramscol)[indices[j]];
          nomj = weights[j] * (colval - rowval);
          // insert contributions into dof j of the gradient
          gradientcol.SumIntoMyValue(indices[j], 0, nomj / denom);
        }
      }
    }

    // bring back gradient to the unique layout he came here with.
    // we have to add up off proc components so we cannot use
    // the LINALG::Export since it only provides CombineMode::Insert
    Epetra_MultiVector tmp(gradient->Map(), gradient->NumVectors(), true);
    Epetra_Export exporter(gradientcol.Map(), tmp.Map());
    int err = tmp.Export(gradientcol, exporter, Add);
    if (err) dserror("Export using exporter returned err=%d", err);

    gradient->Update(tvdweight_, tmp, 1.0);
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PATRegula::FillTvdMatrix()
{
  std::map<std::vector<int>, std::vector<int>>
      facemap;                                    // map faces to corresponding elements/parameters
  std::map<std::vector<int>, double> faceweight;  // map of faces to its area

  // loop all elements
  for (int e = 0; e < discret_->NumMyRowElements(); ++e)
  {
    // get the element
    DRT::Element* ele = discret_->lRowElement(e);
    int pgid = ele->Id();

    // decide whether we deal with a 2D or 3D discretization
    unsigned int nele = 0;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    std::vector<std::vector<int>> faces;
    if (ele->NumSurface() > 1)  // 2D boundary element and 3D parent element
    {
      nele = ele->NumSurface();
      faces = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
    }
    else if (ele->NumSurface() == 1)  // 1D boundary element and 2D parent element
    {
      nele = ele->NumLine();
      faces = DRT::UTILS::getEleNodeNumberingLines(distype);
    }
    else
      dserror("creating internal faces for 1D elements (would be points) not implemented");

    // safety check
    if (nele != faces.size()) dserror("number of surfaces or lines does not match!");

    // get the surface/line elements for area computation
    std::vector<Teuchos::RCP<DRT::Element>> surfs;
    if (ele->NumSurface() > 1)
      surfs = ele->Surfaces();
    else if (ele->NumSurface() == 1)
      surfs = ele->Lines();
    else
      dserror("creating internal faces for 1D elements (would be points) not implemented");

    // get nodes of each of this element's face
    for (unsigned int iele = 0; iele < nele; iele++)
    {
      // allocate node vectors
      unsigned int nnode = faces[iele].size();
      std::vector<int> nodeids(nnode);

      // get connectivity info
      for (unsigned int inode = 0; inode < nnode; inode++)
        nodeids[inode] = ele->NodeIds()[faces[iele][inode]];

      // get the area of this face
      Teuchos::ParameterList p;
      // SetAction(p);
      if (discret_->Name() == "scatra")
        p.set<int>("action", SCATRA::bd_integrate_shape_functions);
      else if (discret_->Name() == "acou")
        p.set<int>("action", ACOU::bd_integrate);
      else
        dserror("what happenend?");

      p.set("area", 0.0);
      DRT::Element::LocationArray la(discret_->NumDofSets());
      surfs[iele]->LocationVector(*discret_, la, false);
      // initialize element vectors
      int ndof = ele->NumNode() * 3;
      Epetra_SerialDenseMatrix elematrix1(ndof, ndof, false);
      Epetra_SerialDenseMatrix elematrix2(ndof, ndof, false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);
      surfs[iele]->Evaluate(
          p, *discret_, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      double area = p.get("area", -1.0);
      if (area < 0.0) dserror("area computation of surface failed");

      // sort the nodes in faces to make them to be used for keys in the facemap.
      // they are not unique (in a parallel layout sense) though since a face can
      // be on multiple processors
      std::sort(nodeids.begin(), nodeids.end());

      // store parameter global id corresponding to this face
      facemap[nodeids].push_back(pgid);
      faceweight[nodeids] = area;
    }
  }  // loop local elements

  /*------------------------------------------------------------------- */
  // STEP 2: gather for each face on each proc the same set of corresponding
  // elements; ntargetprocs is equal to the total number of processors to make
  // data redundant on all procs
  // the face weight dont need to be communicated since they should be the
  // same on every proc having a face

  const int numprocs = discret_->Comm().NumProc();
  int allproc[numprocs];
  for (int i = 0; i < numprocs; ++i) allproc[i] = i;

  for (int i = 0; i < discret_->Comm().NumProc(); i++)
  {
    // by looping the procs one by one we need to make sure that every proc
    // does the same number of loops regardless of its own lenght of keys
    // in the facemap. so the actual number of loops is defined by the current
    // proc and distributed in localmapsize
    int localmapsize;
    if (discret_->Comm().MyPID() == i) localmapsize = facemap.size();
    discret_->Comm().Broadcast(&localmapsize, 1, i);

    // now iterate as often as there are faces on proc i
    std::map<std::vector<int>, std::vector<int>>::iterator face_it(facemap.begin());
    for (int j = 0; j < localmapsize; j++)
    {
      // get length of current face-key on all procs
      int keylength = 0;
      if (discret_->Comm().MyPID() == i) keylength = face_it->first.size();
      discret_->Comm().Broadcast(&(keylength), 1, i);

      // send current face-key to all procs
      std::vector<int> facekey(keylength, 0);
      if (discret_->Comm().MyPID() == i) facekey = face_it->first;
      discret_->Comm().Broadcast(&(facekey[0]), keylength, i);

      // check whether one of the other procs also has this key and write IDs of
      // procs who own this face in "sowningprocs" and distribute this knowledge
      // among all procs to rowningprocs
      std::map<std::vector<int>, std::vector<int>>::iterator face_abroad(facemap.find(facekey));
      std::vector<int> sowningprocs;
      std::vector<int> rowningprocs;
      if (face_abroad != facemap.end() && discret_->Comm().MyPID() != i)
        sowningprocs.push_back(discret_->Comm().MyPID());
      LINALG::Gather(sowningprocs, rowningprocs, numprocs, allproc, discret_->Comm());

      // now bring parameters corresponding to this face on the other procs to proc i
      // (they are send to all procs but only proc i stores them in the map with
      // the current face-key)
      std::vector<int> sparams;
      std::vector<int> rparams;
      if (std::find(rowningprocs.begin(), rowningprocs.end(), discret_->Comm().MyPID()) !=
          rowningprocs.end())
        sparams = facemap[facekey];
      LINALG::Gather(sparams, rparams, numprocs, allproc, discret_->Comm());

      // store additional elements on proc i
      if (discret_->Comm().MyPID() == i)
      {
        for (int iele = 0; iele < (int)rparams.size(); iele++)
        {
          // depending on which proc comes first it might be that parameters
          // are added redundantly to a face which leads to undesired summation
          // upon inserting weights into the matrix. So only add if not existent yet
          if (std::find(facemap[facekey].begin(), facemap[facekey].end(), rparams[iele]) ==
              facemap[facekey].end())
            facemap[facekey].push_back(rparams[iele]);
        }
      }

      // increase face pointer on this proc
      if (discret_->Comm().MyPID() == i) face_it++;

    }  // faces on each proc
  }    // procs

  /*------------------------------------------------------------------- */
  // STEP 3: sort elements into graph according to neighbor information
  // in facemap
  std::map<std::vector<int>, std::vector<int>>::iterator faces;
  for (faces = facemap.begin(); faces != facemap.end(); faces++)
  {
    // all parameter ids connected to this face
    std::vector<int> parameters = faces->second;

    // weight corresponding to these parameters
    std::vector<double> weights(parameters.size(), faceweight[faces->first]);

    for (int iele = 0; iele < (int)parameters.size(); iele++)
    {
      int globalrow = parameters[iele];
      if (discret_->ElementRowMap()->LID(globalrow) >= 0)
      {
        // like this the diagonal entries are inserted redundantly and summed up
        // after FillComplete() is called; they are more or less useless anyways
        tvdmatrix_->InsertGlobalValues(globalrow, parameters.size(), &weights[0], &parameters[0]);
      }
    }
  }

  // finalize the graph ...
  tvdmatrix_->FillComplete();
  tvdmatrix_->OptimizeStorage();

  // delete diagonal values
  Teuchos::RCP<Epetra_Vector> diagonal =
      Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap(), true));
  tvdmatrix_->ReplaceDiagonalValues(*diagonal);

  return;
}
