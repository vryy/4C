/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_timemonitor.cpp

\brief Fast time monitor. E. g. to measure the element evaluation times.

\maintainer Michael Hiermeier

\date May 23, 2018

\level 3
*/
/*----------------------------------------------------------------------------*/

#include "contact_aug_timemonitor.H"
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_Time.hpp>
#include <Epetra_Comm.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename enum_class>
CONTACT::AUG::TimeMonitor<enum_class>::TimeMonitor()
{
  static_assert(std::is_same<unsigned, typename std::underlying_type<enum_class>::type>::value,
      "The template ENUM_CLASS must use UNSIGNED INT as underlying type!");

  timings_.resize(static_cast<unsigned>(enum_class::MAX_TIME_ID), std::make_pair(-1.0, 0.0));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename enum_class>
void CONTACT::AUG::TimeMonitor<enum_class>::TimeMonitor::reset()
{
  std::fill(timings_.begin(), timings_.end(), std::make_pair(-1.0, 0.0));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename enum_class>
void CONTACT::AUG::TimeMonitor<enum_class>::start(const enum_class id)
{
  timings_[static_cast<unsigned>(id)].first = Teuchos::Time::wallTime();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename enum_class>
void CONTACT::AUG::TimeMonitor<enum_class>::stop(const enum_class id)
{
  std::pair<double, double>& begin_time = timings_[static_cast<unsigned>(id)];
  if (begin_time.first == -1.0) dserror("Call start() first!");

  double& accumulated_time = begin_time.second;
  last_incr_ = Teuchos::Time::wallTime() - begin_time.first;
  accumulated_time += last_incr_;

  begin_time.first = -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename enum_class>
double CONTACT::AUG::TimeMonitor<enum_class>::getMyTotalTime() const
{
  double my_total_time = 0.0;
  for (auto& t : timings_) my_total_time += t.second;

  return my_total_time;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename enum_class>
void CONTACT::AUG::TimeMonitor<enum_class>::write(std::ostream& os)
{
  int mypid = 0;
  if (comm_) mypid = comm_->MyPID();

  if (mypid == 0)
  {
    os << std::string(100, '=') << std::endl;
    os << "CONTACT::AUG::TimeMonitor - Final Overview:\n";
  }
  for (unsigned i = 0; i < static_cast<unsigned>(enum_class::MAX_TIME_ID); ++i)
  {
    double gtime = 0.0;
    if (comm_)
      comm_->SumAll(&timings_[i].second, &gtime, 1);
    else
      gtime = timings_[i].second;

    if (gtime == 0.0) continue;

    std::string name = TimeID2Str(static_cast<enum_class>(i));
    if (mypid == 0)
    {
      os << std::string(100, '-') << "\n";
      os << "TOTAL - " << std::left << std::setw(72) << name << ": " << std::scientific
         << std::setprecision(5) << gtime << " [sec.]\n";
    }

    if (comm_)
    {
      const int numproc = comm_->NumProc();
      std::vector<double> lproc_timings(numproc, 0.0);
      lproc_timings[mypid] = timings_[i].second;

      std::vector<double> gproc_timings(numproc, 0.0);
      comm_->SumAll(lproc_timings.data(), gproc_timings.data(), numproc);

      if (mypid == 0)
      {
        for (int p = 0; p < numproc; ++p)
        {
          os << "proc #" << std::right << std::setw(3) << p << " - " << std::left << std::setw(68)
             << name << ": " << std::scientific << std::setprecision(5) << gproc_timings[p]
             << " [sec.]\n"
             << std::left;
        }
      }
    }
  }
  if (mypid == 0) os << std::string(100, '=') << std::endl;

  // wait till everything is done
  if (comm_) comm_->Barrier();
}

/*----------------------------------------------------------------------------*/
template class CONTACT::AUG::TimeMonitor<CONTACT::AUG::TimeID>;
template class CONTACT::AUG::TimeMonitor<CONTACT::AUG::GlobalTimeID>;
