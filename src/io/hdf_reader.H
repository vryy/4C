// -*- C++ -*-
#ifndef HDF_READER_H
#define HDF_READER_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <H5LT.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "../discret/drt_discret.H"
#include "Teuchos_RefCountPtr.hpp"

extern "C" {
#include "../headers/standardtypes.h"
}
#include "../discret/dstrc.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

class HDFReader
{
public:
  HDFReader(string dir, int num_proc=1);
  HDFReader();
  ~HDFReader();
  void open(string basename);
  RefCountPtr<std::vector<char> > read_element_data(int step, int new_proc_num=0,
                                                    int my_id=0);
  RefCountPtr<std::vector<char> > read_node_data(int step, int new_proc_num=0,
                                                 int my_id=0);
  RefCountPtr<Epetra_Vector>
  read_result_data(string id_path, string value_path, int new_proc_num, int my_id,
                   RefCountPtr<Epetra_Comm> Comm, int gId_num);
  void close();
  void set_num_output_procs(int num);
private:
  RefCountPtr<std::vector<char> > read_char_data(std::string path,
                                                 int start=0, int end=-1);
  RefCountPtr<std::vector<int> > read_int_data(std::string path,
                                               int start=0, int end=-1);
  RefCountPtr<std::vector<double> > read_double_data(std::string path,
                                                     int start=0, int end=-1);
  void calculate_range(int new_proc_num, int my_id, int& start, int& end);
  std::vector<std::string> filenames_;
  std::vector<hid_t> files_;
  string input_dir_;
  int num_output_proc_;
};


#endif
#endif
#endif // HDF_READER_H
