#ifdef CCADISCRET
#ifdef BINIO

#include <iostream>
#include "hdf_reader.H"

using namespace std;

/*----------------------------------------------------------------------*
 * The Constructor of the HDFReader (num_proc defaults to 1)
 *----------------------------------------------------------------------*/
IO::HDFReader::HDFReader(string dir):
  filenames_(0),
  files_(0),
  input_dir_(dir),
  num_output_proc_(0)
{
}

/*----------------------------------------------------------------------*
 * The Destructor
 *----------------------------------------------------------------------*/
IO::HDFReader::~HDFReader()
{
  Close();
}

/*----------------------------------------------------------------------*
 * With num_output_proc_ == 1 this function opens the result data file
 * with name basename. When num_output_proc_ > 1 it opens the result
 * files of all processors, by appending .p<proc_num> to the basename.
 *----------------------------------------------------------------------*/
void IO::HDFReader::Open(string basename,int num_output_procs)
{
  Close();
  num_output_proc_ = num_output_procs;
  for (int i = 0; i < num_output_proc_; ++i)
  {
    ostringstream buf;
    buf << input_dir_ << basename;
    if (num_output_proc_>1)
    {
      buf << ".p" << i;
    }
    filenames_.push_back(buf.str());
    files_.push_back(H5Fopen(buf.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    if (files_[i] < 0)
      dserror("Failed to open HDF-file %s", filenames_[i].c_str());
  }
}
/*----------------------------------------------------------------------*
 * reads the packed element data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files
 *----------------------------------------------------------------------*/
RefCountPtr<std::vector<char> > IO::HDFReader::ReadElementData(int step, int new_proc_num, int my_id) const
{
  if (files_.size()==0 || files_[0] == -1)
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
    CalculateRange(new_proc_num,my_id,start,end);
  }
  return ReadCharData(path.str(),start,end);
}

#if 1
/*----------------------------------------------------------------------*
 * reads the packed bc data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files                                          gammi 05/07
 *----------------------------------------------------------------------*/
RefCountPtr<std::vector<char> > IO::HDFReader::ReadCondition(
        const int step,
        const int new_proc_num,
        const int my_id,
        const string condname) const
{
  if (files_.size()==0 || files_[0] == -1)
    dserror("Tried to read data without opening any file");

  ostringstream path;
  path << "/step" << step << "/" << condname;

  // only one proc (PROC 0) wrote this conditions to the mesh file
  const int start =0;
  const int end   =1;

  /* Save old error handler */
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto(NULL, NULL);

  /* Probe. Likely to fail, but that's okay */
  const hid_t dataset = H5Dopen(files_[0],(path.str()).c_str());

  /* Restore previous error handler */
  H5Eset_auto(old_func, old_client_data);

  RefCountPtr<std::vector<char> > block = rcp(new std::vector<char>());
  if (dataset > -1)
  {
      block = ReadCharData(path.str(),start,end);
  }

  return block;
}
#endif

/*----------------------------------------------------------------------*
 * reads the packed node data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files
 *----------------------------------------------------------------------*/
RefCountPtr<std::vector<char> > IO::HDFReader::ReadNodeData(
        int step,
        int new_proc_num,
        int my_id) const
{
  if (files_.size()==0 || files_[0] == -1)
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
    CalculateRange(new_proc_num,my_id,start,end);
  }
  RefCountPtr<vector<char> > d = ReadCharData(path.str(),start,end);
  return d;
}

/*----------------------------------------------------------------------*
 * reads the dataset 'path' in all the files in the range [start,end)
 * and returns all the data in one vector. The data is assumed to by
 * of type char (private)
 *----------------------------------------------------------------------*/
RefCountPtr<std::vector<char> >
IO::HDFReader::ReadCharData(string path, int start, int end) const
{
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  RefCountPtr<vector<char> > data = rcp(new vector<char>);
  for (int i = start; i < end; ++i)
  {
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
RefCountPtr<std::vector<int> >
IO::HDFReader::ReadIntData(string path, int start, int end) const
{
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  RefCountPtr<vector<int> > data = rcp(new vector<int>);
  for (int i = start; i < end; ++i)
  {
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
      dserror("Failed to close node dataspace %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    status = H5Dclose(dataset);
    if (status < 0)
      dserror("Failed to close node dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
  }
  return data;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<std::vector<double> >
IO::HDFReader::ReadDoubleData(string path, int start, int end, std::vector<int>& lengths) const
{
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  RefCountPtr<vector<double> > data = rcp(new vector<double>);
  for (int i = start; i < end; ++i)
  {
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
    lengths.push_back(dim);
    offset += dim;
    if (status < 0)
      dserror("Failed to read data from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    status = H5Sclose(dataspace);
    if (status < 0)
      dserror("Failed to close node dataspace %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    status = H5Dclose(dataset);
    if (status < 0)
      dserror("Failed to close node dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
  }
  return data;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<Epetra_MultiVector>
IO::HDFReader::ReadResultData(string id_path, string value_path, int columns, const Epetra_Comm& Comm) const
{
  int new_proc_num = Comm.NumProc();
  int my_id = Comm.MyPID();

  if (files_.size()==0 || files_[0] == -1)
    dserror("Tried to read data without opening any file");
  int start, end;
  CalculateRange(new_proc_num,my_id,start,end);

  RefCountPtr<vector<int> > ids = ReadIntData(id_path,start,end);
  Epetra_Map map(-1,static_cast<int>(ids->size()), &((*ids)[0]),0,Comm);

  RefCountPtr<Epetra_MultiVector> res;
  if (columns==1)
    res = rcp(new Epetra_Vector(map,false));
  else
    res = rcp(new Epetra_MultiVector(map,columns,false));

  std::vector<int> lengths;
  RCP<std::vector<double> > values = ReadDoubleData(value_path,start,end,lengths);

  if (static_cast<int>(values->size()) != res->MyLength()*res->NumVectors())
    dserror("vector value size mismatch: %d != %d",values->size(),res->MyLength()*res->NumVectors());

  // Rearrange multi vectors that are read with fewer processors than written.
  int offset = 0;
  for (int i = start; i < end; ++i)
  {
    int l = lengths[i-start];
    for (int c=0; c<columns; ++c)
    {
      copy(&(*values)[offset+ c   *l/columns],
           &(*values)[offset+(c+1)*l/columns],
           &res->Values()[c*res->MyLength()+offset/columns]);
    }
    offset += l;
  }

  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::HDFReader::Close()
{
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
  filenames_.resize(0);
  files_.resize(0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::HDFReader::CalculateRange(int new_proc_num, int my_id, int& start, int& end) const
{
  const int mod = num_output_proc_ % new_proc_num;
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
          (num_output_proc_ / new_proc_num)*(my_id - mod + 1);
  }
}

#endif
#endif
