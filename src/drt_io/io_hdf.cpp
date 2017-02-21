/*----------------------------------------------------------------------*/
/*!
\file io_hdf.cpp

\brief Helpers to read HDF5 based output.

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/


#include <iostream>
#include "io_hdf.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 * The Constructor of the HDFReader (num_proc defaults to 1)
 *----------------------------------------------------------------------*/
IO::HDFReader::HDFReader(std::string dir):
  filenames_(0),
  files_(0),
  input_dir_(dir),
  num_output_proc_(0)
{
  // inhibit delayed closure, throws error if file contents still in use
  H5Plist_ = H5Pcreate(H5P_FILE_ACCESS);
  herr_t status = H5Pset_fclose_degree( H5Plist_, H5F_CLOSE_WEAK);
  if (status < 0)
    dserror("Failed to set file access list");
}

/*----------------------------------------------------------------------*
 * The Destructor
 *----------------------------------------------------------------------*/
IO::HDFReader::~HDFReader()
{
  Close();
  herr_t status = H5Pclose(H5Plist_);
  if (status < 0)
    dserror("Failed to close file access list");
}

/*----------------------------------------------------------------------*
 * With num_output_proc_ == 1 this function opens the result data file
 * with name basename. When num_output_proc_ > 1 it opens the result
 * files of all processors, by appending .p<proc_num> to the basename.
 *----------------------------------------------------------------------*/
void IO::HDFReader::Open(std::string basename,int num_output_procs,int new_proc_num,int my_id)
{
  int start;
  int end;
  num_output_proc_ = num_output_procs;
  CalculateRange(new_proc_num, my_id, start, end);
  Close();
  for (int i = 0; i < num_output_proc_; ++i)
  {
    std::ostringstream buf;
    buf << input_dir_ << basename;
    if (num_output_proc_>1)
    {
      buf << ".p" << i;
    }
    if (i>=start and i<end)
    {
      filenames_.push_back(buf.str());
      files_.push_back(H5Fopen(buf.str().c_str(), H5F_ACC_RDONLY, H5Plist_));
      if (files_[i] < 0)
        dserror("Failed to open HDF-file %s", filenames_[i].c_str());
    }
    else
    {
      filenames_.push_back("");
      files_.push_back(-1);
    }
  }
}
/*----------------------------------------------------------------------*
 * reads the packed element data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> > IO::HDFReader::ReadElementData(int step, int new_proc_num, int my_id) const
{
  if (files_.size()==0)
    dserror("Tried to read data without opening any file");
  std::ostringstream path;
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


/*----------------------------------------------------------------------*
 * reads the packed bc data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files                                          gammi 05/07
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> > IO::HDFReader::ReadCondition(
        const int step,
        const int new_proc_num,
        const int my_id,
        const std::string condname) const
{
  if (files_.size()==0)
    dserror("Tried to read data without opening any file");

  std::ostringstream path;
  path << "/step" << step << "/" << condname;

  // only one proc (PROC 0) wrote this conditions to the mesh file
  const int start =0;
  const int end   =1;

  Teuchos::RCP<std::vector<char> > block;
  block = ReadCharData(path.str(),start,end);

  return block;
}


/*----------------------------------------------------------------------*
 * reads the packed knotvector data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files                                          gammi 05/08
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> > IO::HDFReader::ReadKnotvector(
  const int step) const
{
  if (files_.size()==0)
    dserror("Tried to read data without opening any file");

  std::ostringstream path;
  path << "/step" << step << "/" << "knotvector";

  // only one proc (PROC 0) wrote the knotvector to the mesh file
  const int start =0;
  const int end   =1;

  Teuchos::RCP<std::vector<char> > block;
  block = ReadCharData(path.str(),start,end);

  return block;
}


/*----------------------------------------------------------------------*
 * reads the packed node data from the mesh files
 * Note: this function should only be called when the HDFReader opened
 * the mesh files
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> > IO::HDFReader::ReadNodeData(
        int step,
        int new_proc_num,
        int my_id) const
{
  if (files_.size()==0)
    dserror("Tried to read data without opening any file");
  std::ostringstream path;
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
  Teuchos::RCP<std::vector<char> > d = ReadCharData(path.str(),start,end);
  return d;
}

/*----------------------------------------------------------------------*
 * reads the dataset 'path' in all the files in the range [start,end)
 * and returns all the data in one vector. The data is assumed to by
 * of type char (private)
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> >
IO::HDFReader::ReadCharData(std::string path, int start, int end) const
{
  if (end == -1)
    end = num_output_proc_;
  hsize_t offset = 0;
  Teuchos::RCP<std::vector<char> > data = Teuchos::rcp(new std::vector<char>);
  for (int i = start; i < end; ++i)
  {
    const char* cpath = path.c_str();
    hid_t dataset = H5Dopen(files_[i],cpath);
    if (dataset < 0)
    {
      dserror("Failed to open dataset %s in HDF-file %s",
              cpath,filenames_[i].c_str());
    }
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      dserror("Failed to get dataspace from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    int rank = H5Sget_simple_extent_ndims(dataspace);
    switch (rank)
    {
    case 0:
      break;
    case 1:
    {
      hsize_t dim,maxdim;
      int res = H5Sget_simple_extent_dims(dataspace,&dim,&maxdim);
      if (res < 0)
        dserror("Failed to get size from dataspace in HDF-file %s",
                filenames_[i].c_str());
      data->resize(offset+dim);
      herr_t status = H5Dread(dataset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
                              H5P_DEFAULT,&((*data)[offset]));
      offset += dim;
      if (status < 0)
        printf("Failed to read data from dataset %s in HDF-file %s. "
            "This can be tolerated in case you have procs without row elements!",
                path.c_str(),filenames_[i].c_str());
      break;
    }
    default:
      dserror("HDF5 rank=%d unsupported", rank);
      break;
    }

    herr_t status = H5Sclose(dataspace);
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
 * and returns all the data in one std::vector<int> (private)
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<int> >
IO::HDFReader::ReadIntData(std::string path, int start, int end) const
{
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  Teuchos::RCP<std::vector<int> > data = Teuchos::rcp(new std::vector<int>);
  for (int i = start; i < end; ++i)
  {
    hid_t dataset = H5Dopen(files_[i],path.c_str());
    if (dataset < 0)
    {
      dserror("Failed to open dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    }
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
      dserror("Failed to get dataspace from dataset %s in HDF-file %s",
              path.c_str(),filenames_[i].c_str());
    int rank = H5Sget_simple_extent_ndims(dataspace);
    switch (rank)
    {
    case 0:
      break;
    case 1:
    {
      hsize_t dim,maxdim;
      int res = H5Sget_simple_extent_dims(dataspace,&dim,&maxdim);
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
      break;
    }
    default:
      dserror("HDF5 rank=%d unsupported", rank);
      break;
    }

    herr_t status = H5Sclose(dataspace);
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
Teuchos::RCP<std::vector<double> >
IO::HDFReader::ReadDoubleData(std::string path, int start, int end, std::vector<int>& lengths) const
{
  if (end == -1)
    end = num_output_proc_;
  int offset = 0;
  Teuchos::RCP<std::vector<double> > data = Teuchos::rcp(new std::vector<double>);
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
    int rank = H5Sget_simple_extent_ndims(dataspace);
    switch (rank)
    {
    case 0:
      lengths.push_back(0);
      break;
    case 1:
    {
      hsize_t dim,maxdim;
      int res = H5Sget_simple_extent_dims(dataspace,&dim,&maxdim);
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
      break;
    }
    default:
      dserror("HDF5 rank=%d unsupported", rank);
      break;
    }

    herr_t status = H5Sclose(dataspace);
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
Teuchos::RCP<Epetra_MultiVector>
IO::HDFReader::ReadResultData(std::string id_path, std::string value_path, int columns, const Epetra_Comm& Comm) const
{
  int new_proc_num = Comm.NumProc();
  int my_id = Comm.MyPID();

  if (files_.size()==0)
    dserror("Tried to read data without opening any file");
  int start, end;
  CalculateRange(new_proc_num,my_id,start,end);

  Teuchos::RCP<std::vector<int> > ids = ReadIntData(id_path,start,end);
  Epetra_Map map(-1,static_cast<int>(ids->size()), &((*ids)[0]),0,Comm);

  Teuchos::RCP<Epetra_MultiVector> res;
  if (columns==1)
    res = Teuchos::rcp(new Epetra_Vector(map,false));
  else
    res = Teuchos::rcp(new Epetra_MultiVector(map,columns,false));

  std::vector<int> lengths;
  Teuchos::RCP<std::vector<double> > values = ReadDoubleData(value_path,start,end,lengths);

  if (static_cast<int>(values->size()) != res->MyLength()*res->NumVectors())
    dserror("vector value size mismatch: %d != %d",values->size(),res->MyLength()*res->NumVectors());

  // Rearrange multi vectors that are read with fewer processors than written.
  int offset = 0;
  for (int i = start; i < end; ++i)
  {
    int l = lengths[i-start];
    for (int c=0; c<columns; ++c)
    {
      std::copy(&(*values)[offset+ c   *l/columns],
           &(*values)[offset+(c+1)*l/columns],
           &res->Values()[c*res->MyLength()+offset/columns]);
    }
    offset += l;
  }

  return res;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> >
IO::HDFReader::ReadResultDataVecChar(std::string id_path, std::string value_path, int columns, const Epetra_Comm& Comm, Teuchos::RCP<Epetra_Map>& elemap) const
{
  if (columns!=1)
    dserror("got multivector, std::vector<char> expected");

  int new_proc_num = Comm.NumProc();
  int my_id = Comm.MyPID();

  if (files_.size()==0)
    dserror("Tried to read data without opening any file");
  int start, end;
  CalculateRange(new_proc_num,my_id,start,end);

  Teuchos::RCP<std::vector<int> > ids = ReadIntData(id_path,start,end);
  //cout << "size of ids:" << (*ids).size() << endl;
  Epetra_Map map(-1,static_cast<int>(ids->size()), &((*ids)[0]),0,Comm);
  elemap = Teuchos::rcp(new Epetra_Map(map));

  Teuchos::RCP<std::vector<char> > res = ReadCharData(value_path,start,end);
  return res;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<char> >
IO::HDFReader::ReadCharVector(std::string value_path, const Epetra_Comm& Comm) const
{
  int new_proc_num = Comm.NumProc();
  int my_id = Comm.MyPID();

  if ( files_.size() == 0 )
    dserror("Tried to read data without opening any file");
  int start, end;
  CalculateRange( new_proc_num, my_id, start, end );

  Teuchos::RCP<std::vector<char> > res = ReadCharData( value_path, start ,end );
  return res;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::HDFReader::Close()
{
  for (int i = 0; i < num_output_proc_ and i<static_cast<int>(files_.size()); ++i)
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


