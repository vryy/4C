# check for advanced features that give rise to extra preprocessor flags

CHECK_CXX_SOURCE_RUNS(
  "
  #include <cxxabi.h>
  #include <stdio.h>
  #include <execinfo.h>
  #include <stdlib.h>
  #include <string>
  #include <stdexcept>
  template <int dummy>
  bool f()
  {
  int nptrs;
  void *buffer[100];
  char **strings;
  nptrs   = backtrace(buffer, 100);
  strings = backtrace_symbols(buffer, nptrs);
  std::string entry (strings[0]);
  std::string applicationname = entry.substr(0,entry.find('('));
  std::string filename;
  const std::size_t address0 = entry.find(\" [\");
  const std::size_t address1 = entry.find(\"]\", address0);
  std::string progname = \"addr2line \" + entry.substr(address0+2,address1-2-address0) + \" -e \" + applicationname;
  FILE *stream = popen( progname.c_str(), \"r\");
  if (stream != 0) {
  char path[4096];
  if (fgets(path, sizeof(path)-1, stream) != 0) {
  std::string pathname (path);
  if (pathname.find('?') == std::string::npos)
  filename = \"  (\" +
  pathname.substr(pathname.find_last_of('/')+1,
  pathname.length()-pathname.find_last_of('/')-2) + ')';
  }
  }
  pclose(stream);
  const std::size_t start = entry.find('('), end = entry.find('+');
  std::string functionname = entry.substr(start+1,end-start-1);
  int status;
  char *p = abi::__cxa_demangle(functionname.c_str(), 0, 0, &status);
  bool completed = false;
  if (status == 0)
  {
  std::string demangled(p);
  if (demangled.find(\"bool f<1>()\") != std::string::npos)
  completed = true;
  }
  free(p);
  free(strings);
  return !completed;
  }
  int main(int argc, char **argv)
  {
  return f<1>();
  }
  "
  ADVANCED_STACKTR)

if(ADVANCED_STACKTR)
  add_definitions("-DENABLE_ADVANCED_STACKTR")
endif(ADVANCED_STACKTR)
