/*---------------------------------------------------------------------*/
/*! \file

\brief printing error messages

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_config.H"

#include "baci_utils_exceptions.H"

#ifdef BACI_WITH_BACKTRACE
#define BOOST_STACKTRACE_USE_BACKTRACE
#endif

#include <boost/stacktrace.hpp>

#include <cstring>
#include <iostream>
#include <sstream>

BACI_NAMESPACE_OPEN

namespace CORE
{
  namespace INTERNAL
  {
    void ErrorHelper::operator()(const std::string& format, ...) const
    {
      this->operator()(format.c_str());
    }

    void ErrorHelper::operator()(const char* format, ...) const
    {
      int initialized;
      MPI_Initialized(&initialized);
      int myrank = 0;
      if (initialized > 0)
      {
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      }

      va_list arglist;
      va_start(arglist, format);
      constexpr int buffer_size = 8096;
      std::array<char, buffer_size> buffer;
      const int written_size = std::vsnprintf(buffer.data(), buffer_size, format, arglist);
      va_end(arglist);

      std::string formatted_msg(buffer.data(), std::min(buffer_size, written_size));

      std::stringstream compound_message;
      compound_message << "PROC " << myrank << " ERROR in " << file << ", line " << line_number
                       << ":\n";
      compound_message << formatted_msg;
      compound_message << "\n------------------\n";

      throw CORE::Exception(compound_message.str());
    }
  }  // namespace INTERNAL

  const char* Exception::what() const noexcept
  {
    what_message_ = message + to_string(*stacktrace);
    return what_message_.c_str();
  }

  // This number tells the stack trace class to skip a certain number of frames that are introduced
  // by our exception framework itself. Users are not interested in these calls, and provding them
  // exposes unnecessary details.
  constexpr std::size_t skip_frames = 2;

  Exception::Exception(std::string message)
      : message(std::move(message)),
        stacktrace(std::make_unique<boost::stacktrace::stacktrace>(skip_frames, /*max_depth=*/-1))
  {
  }

  // Default the destructor here to facilitate an incomplete stacktrace type as member variable.
  Exception::~Exception() = default;
}  // namespace CORE

BACI_NAMESPACE_CLOSE
