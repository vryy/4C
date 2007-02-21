#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "hdf_reader.H"

using namespace std;

/*----------------------------------------------------------------------*
 * The Constructor of the HDFReader (num_proc defaults to 1)
 *----------------------------------------------------------------------*/
HDFReader::HDFReader(string dir, int num_proc):
  filenames_(0),
  files_(0),
  input_dir_(dir),
  num_output_proc_(num_proc)
{
  DSTraceHelper("HDFReader::HDFReader");
  if (num_proc < 1)
  {
    dserror("Invalid value for num_output_proc: %i, in HDFReader",num_proc);
  }
  files_.resize(num_proc,-1);
  filenames_.resize(num_proc);
}

/*----------------------------------------------------------------------*
 * Another Constructor without any arguments
 *----------------------------------------------------------------------*/
HDFReader::HDFReader():
  filenames_(),
  files_(),
  num_output_proc_(1)
{
  DSTraceHelper("HDFReader::HDFReader");
  files_.resize(1,-1);
  filenames_.resize(1);
}


/*----------------------------------------------------------------------*
 * The Destructor
 *----------------------------------------------------------------------*/
HDFReader::~HDFReader()
{
  DSTraceHelper("HDFReader::~HDFReader");
  close();
}

/*----------------------------------------------------------------------*
 * With num_output_proc_ == 1 this function opens the result data file
 * with name basename. When num_output_proc_ > 1 it opens the result
 * files of all processors, by appending .p<proc_num> to the basename.
 *----------------------------------------------------------------------*/
void HDFReader::open(string basename)
{
  DSTraceHelper("HDFReader::open");
  for (int i = 0; i < num_output_proc_; ++i)
  {
    if (files_[i] != -1)
    {
      herr_t status = H5Fclose(files_[i]);
      if (status < 0)
        dserror("Failed to close HDF-file %s", filenames_[i].c_str());
    }
    ostringstream buf;
    if (num_output_proc_>1)
    {
      buf << input_dir_
          << basename
          << ".p" << i;
    }
    else
    {
      buf << input_dir_
          << basename;
    }
    filenames_[i] = buf.str();
    files_[i] = H5Fopen(buf.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (files_[i] < 0)
      dserror("Failed to open HDF-file %s", filenames_[i].c_str());
  }
}
/*----------------------------------------------------------------------*
 * reads the packed element data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files
 *----------------------------------------------------------------------*/
RefCountPtr<vector<char> > HDFReader::read_element_data(int step, int new_proc_num,
                                                        int my_id)
{
  DSTraceHelper("HDFReader::read_element_data");
  if (files_[0] == -1)
    dserror("Tried to read data without opening any file");
  ostringstream path;
  path << "/step" << step << "/elements";
  int start,end;
  if (new_proc_num == 0 && my_id == 0)
  {
    start = 0;
    end = num_output_proc_;
  }
  else
  {
    calculate_range(new_proc_num,my_id,start,end);
  }
  return read_char_data(path.str(),start,end);
}

/*----------------------------------------------------------------------*
 * reads the packed node data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files
 *----------------------------------------------------------------------*/
RefCountPtr<vector<char> > HDFReader::read_node_data(int step, int new_proc_num,
                                                     int my_id)
{
  DSTraceHelper("HDFReader::read_node_data");
  if (files_[0] == -1)
    dserror("Tried to read data without opening any file");
  ostringstream path;
  path << "/step" << step << "/nodes";
  int start,end;
  if (new_proc_num == 0 && my_id == 0)
  {
    start = 0;
    end = num_output_proc_;
  }
  else
  {
    calculate_range(new_proc_num,my_id,start,end);
  }
  RefCountPtr<vector<char> > d = read_char_data(path.str(),start,end);
  return d;
}

/*----------------------------------------------------------------------*
 * reads the dataset 'path' in all the files in the range [start,end)
 * and returns all the data in one vector. The data is assumed to by
 * of type char (private)
 *----------------------------------------------------------------------*/
RefCountPtr<vector<char> >
HDFReader::read_char_data(string path, int start, int end)
{
  DSTraceHelper("HDFReader::read_char_data");
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  RefCountPtr<vector<char> > data = rcp(new vector<char>);
  for (int i = start; i < end; ++i)
  {
    DSTraceHelper("HDFReader::read_char_data");
    hid_t dataset = H5Dopen(files_[i],path.c_str());
    if (dataset < 0)
      dserror("Failed to open dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      dserror("Failed to get dataspace from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    hsize_t dim,maxdim;
    hsize_t res = H5Sget_simple_extent_dims(dataspace,&dim,&maxdim);
    if (res < 0)
      dserror("Failed to get size from dataspace in HDF-file %s",
              filenames_[i].c_str());
    data->resize(offset+dim);
    herr_t status = H5Dread(dataset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
                            H5P_DEFAULT,&((*data)[offset]));
    offset += dim;
    if (status < 0)
      dserror("Failed to read data from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    status = H5Sclose(dataspace);
    if (status < 0)
      dserror("Failed to close node dataspace",
              path.c_str(),filenames_[i].c_str());
    status = H5Dclose(dataset);
    if (status < 0)
      dserror("Failed to close node dataset",
              path.c_str(),filenames_[i].c_str());
  }
  return data;
}

/*----------------------------------------------------------------------*
 * reads the dataset 'path' in all the files in the range [start,end)
 * and returns all the data in one vector<int> (private)
 *----------------------------------------------------------------------*/
RefCountPtr<vector<int> >
HDFReader::read_int_data(string path, int start, int end)
{
  DSTraceHelper("HDFReader::read_int_data");
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  RefCountPtr<vector<int> > data = rcp(new vector<int>);
  for (int i = start; i < end; ++i)
  {
    DSTraceHelper("HDFReader::read_int_data");
    hid_t dataset = H5Dopen(files_[i],path.c_str());
    if (dataset < 0)
      dserror("Failed to open dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      dserror("Failed to get dataspace from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    hsize_t dim,maxdim;
    hsize_t res = H5Sget_simple_extent_dims(dataspace,&dim,&maxdim);
    if (res < 0)
      dserror("Failed to get size from dataspace in HDF-file %s",
              filenames_[i].c_str());
    data->resize(offset+dim);
    herr_t status = H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
                            H5P_DEFAULT,&((*data)[offset]));
    offset += dim;
    if (status < 0)
      dserror("Failed to read data from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    status = H5Sclose(dataspace);
    if (status < 0)
      dserror("Failed to close node dataspace",
              path.c_str(),filenames_[i].c_str());
    status = H5Dclose(dataset);
    if (status < 0)
      dserror("Failed to close node dataset",
              path.c_str(),filenames_[i].c_str());
  }
  return data;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<vector<double> >
HDFReader::read_double_data(string path, int start, int end)
{
  DSTraceHelper("HDFReader::read_double_data");
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  RefCountPtr<vector<double> > data = rcp(new vector<double>);
  for (int i = start; i < end; ++i)
  {
    DSTraceHelper("HDFReader::read_double_data");
    hid_t dataset = H5Dopen(files_[i],path.c_str());
    if (dataset < 0)
      dserror("Failed to open dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      dserror("Failed to get dataspace from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    hsize_t dim,maxdim;
    hsize_t res = H5Sget_simple_extent_dims(dataspace,&dim,&maxdim);
    if (res < 0)
      dserror("Failed to get size from dataspace in HDF-file %s",
              filenames_[i].c_str());
    data->resize(offset+dim);
    herr_t status = H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
                            H5P_DEFAULT,&((*data)[offset]));
    offset += dim;
    if (status < 0)
      dserror("Failed to read data from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    status = H5Sclose(dataspace);
    if (status < 0)
      dserror("Failed to close node dataspace",
              path.c_str(),filenames_[i].c_str());
    status = H5Dclose(dataset);
    if (status < 0)
      dserror("Failed to close node dataset",
              path.c_str(),filenames_[i].c_str());
  }
  return data;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector>
HDFReader::read_result_data(string id_path, string value_path, int new_proc_num,
                            int my_id, RefCountPtr<Epetra_Comm> Comm, int gId_num)
{
  DSTraceHelper("HDFReader::read_result_data");
  if (files_[0] == -1)
    dserror("Tried to read data without opening any file");
  int start, end;
  calculate_range(new_proc_num,my_id,start,end);
  RefCountPtr<vector<int> > ids;
  ids = read_int_data(id_path,start,end);
  Epetra_Map map = Epetra_Map(gId_num,static_cast<int>(ids->size()),
                              &((*ids)[0]),0,*Comm);
  RefCountPtr<Epetra_Vector> res = rcp(new Epetra_Vector(map,false));
  RefCountPtr<vector<double> > values;
  values = read_double_data(value_path,start,end);
  copy(values->begin(),values->end(),res->Values());
  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void HDFReader::close()
{
  DSTraceHelper("HDFReader::close");
  for (int i = 0; i < num_output_proc_; ++i)
  {
    if (files_[i] != -1)
    {
      herr_t status = H5Fclose(files_[i]);
      if (status < 0)
        dserror("Failed to close HDF-file %s", filenames_[i].c_str());
      files_[i] = -1;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void HDFReader::calculate_range(int new_proc_num, int my_id, int& start, int& end)
{
  int mod = num_output_proc_ % new_proc_num;
  if (my_id < mod)
  {
    start = (num_output_proc_ / new_proc_num + 1)*my_id;
    end = (num_output_proc_ / new_proc_num + 1)*(my_id+1);
  }
  else
  {
    start = (num_output_proc_ / new_proc_num + 1)*mod +
            (num_output_proc_ / new_proc_num)*(my_id - mod);
    end = (num_output_proc_ / new_proc_num + 1)*mod +
            (num_output_proc_ / new_proc_num)*(my_id - mod + 1);;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void HDFReader::set_num_output_procs(int num)
{
  DSTraceHelper("HDFReader::set_num_output_procs");
  close();
  num_output_proc_ = num;
  files_.resize(num,-1);
  filenames_.resize(num);
}

#endif
#endif
