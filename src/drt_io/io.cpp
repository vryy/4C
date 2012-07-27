/*----------------------------------------------------------------------*/
/*!
 * \file io.cpp
\brief output context of one discretization

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#include "io.H"
#include "io_control.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../pss_full/pss_cpp.h" // access to legacy parser module

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationReader::DiscretizationReader(Teuchos::RCP<DRT::Discretization> dis,
                                               Teuchos::RCP<IO::InputControl> input,
                                               int step)
  : dis_(dis),
    input_(input)
{
#ifdef BINIO
  FindResultGroup(step, input_->ControlFile());
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationReader::DiscretizationReader(Teuchos::RCP<DRT::Discretization> dis,
                                               int step)
  : dis_(dis),
    input_(DRT::Problem::Instance()->InputControlFile())
{
#ifdef BINIO
  FindResultGroup(step, input_->ControlFile());
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::ReadVector(Teuchos::RCP<Epetra_Vector> vec, std::string name)
{
#ifdef BINIO
  MAP* result = map_read_map(restart_step_, name.c_str());
  int columns;
  if (map_find_int(result,"columns",&columns))
  {
    if (columns != 1)
      dserror("got multivector with name '%s', vector expected", name.c_str());
  }
  ReadMultiVector(vec, name);
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::ReadSerialDenseMatrix(Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > mapdata,
                                                     std::string name)
{
#ifdef BINIO
  MAP* result = map_read_map(restart_step_, name.c_str());
  std::string id_path = map_read_string(result, "ids");
  std::string value_path = map_read_string(result, "values");
  int columns = map_find_int(result,"columns",&columns);
  if (not map_find_int(result,"columns",&columns))
  {
    columns = 1;
  }
  if (columns != 1)
    dserror("got multivector with name '%s', vector<char> expected", name.c_str());

  Teuchos::RCP<Epetra_Map> elemap;
  Teuchos::RCP<std::vector<char> > data = reader_->ReadResultDataVecChar(id_path, value_path, columns,
                                                                        dis_->Comm(), elemap);

  std::vector<char>::size_type position=0;
  for (int i=0;i<elemap->NumMyElements();++i)
  {
    Teuchos::RCP<Epetra_SerialDenseMatrix> matrix = Teuchos::rcp(new Epetra_SerialDenseMatrix);
    DRT::ParObject::ExtractfromPack(position, *data, *matrix);
    (*mapdata)[elemap->GID(i)]=matrix;
  }
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::ReadMultiVector(Teuchos::RCP<Epetra_MultiVector> vec, std::string name)
{
#ifdef BINIO
  MAP* result = map_read_map(restart_step_, name.c_str());
  std::string id_path = map_read_string(result, "ids");
  std::string value_path = map_read_string(result, "values");
  int columns;
  if (not map_find_int(result,"columns",&columns))
  {
    columns = 1;
  }
  Teuchos::RCP<Epetra_MultiVector> nv = reader_->ReadResultData(id_path, value_path, columns, dis_->Comm());
  LINALG::Export(*nv, *vec);
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::ReadMesh(int step)
{
#ifdef BINIO
  FindMeshGroup(step, input_->ControlFile());

  Teuchos::RCP<std::vector<char> > nodedata =
    meshreader_->ReadNodeData(step,dis_->Comm().NumProc(),dis_->Comm().MyPID());

  Teuchos::RCP<std::vector<char> > elementdata =
    meshreader_->ReadElementData(step,dis_->Comm().NumProc(),dis_->Comm().MyPID());

  // before we unpack nodes/elements we store a copy of the nodal row/col map
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(*dis_->NodeRowMap()));
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(*dis_->NodeColMap()));

  // unpack nodes and elements and redistirbuted to current layout

  // take care --- we are just adding elements to the discretisation
  // that means depending on the current distribution and the
  // distribution of the data read we might increase the
  // number of elements in dis_
  // the call to redistribute deletes the unnecessary elements,
  // so everything should be OK
  dis_->UnPackMyNodes(nodedata);
  dis_->UnPackMyElements(elementdata);
  dis_->Redistribute(*noderowmap,*nodecolmap);

#endif
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::ReadRedundantDoubleVector( Teuchos::RCP<std::vector<double> >& doublevec,
							  const std::string name)
{
#ifdef BINIO
  int length;

  if (dis_->Comm().MyPID() == 0)
  {
    // only proc0 reads the vector entities
    // was muß aber jetzt hier wirklich geschehen???

    MAP* result = map_read_map(restart_step_, name.c_str());
    std::string value_path = map_read_string(result, "values");

    doublevec = reader_->ReadDoubleVector(value_path);

    length = doublevec->size();
  }

  // communicate the length of the vector to come
  dis_->Comm().Broadcast(&length,1,0);

  // make vector having the correct length on all procs
  doublevec->resize(length);

  // now distribute information to all procs
  dis_->Comm().Broadcast ( &((*doublevec)[0]), length, 0 );
#endif
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IO::DiscretizationReader::ReadInt(std::string name)
{
  return map_read_int(restart_step_, name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double IO::DiscretizationReader::ReadDouble(std::string name)
{
  return map_read_real(restart_step_, name.c_str());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::FindGroup(int step,
                                         MAP* file,
                                         const char* caption,
                                         const char* filestring,
                                         MAP*& result_info,
                                         MAP*& file_info)
{
#ifdef BINIO

  SYMBOL *symbol;

  /*
   * Iterate all symbols under the name "result" and get the one that
   * matches the given step. Note that this iteration starts from the
   * last result group and goes backward. */

  std::string name = dis_->Name();

  symbol = map_find_symbol(file, caption);
  while (symbol != NULL)
  {
    if (symbol_is_map(symbol))
    {
      MAP* map;
      symbol_get_map(symbol, &map);
      if (map_has_string(map, "field", name.c_str()) and
          map_has_int(map, "step", step))
      {
        result_info = map;
        break;
      }
    }
    symbol = symbol->next;
  }
  if (symbol == NULL)
  {
    dserror("No restart entry for discretization '%s' step %d in symbol table. "
            "Control file corrupt?",
            name.c_str(),step);
  }

  /*--------------------------------------------------------------------*/
  /* open file to read */

  /* We have a symbol and its map that corresponds to the step we are
   * interessted in. Now we need to continue our search to find the
   * step that defines the output file used for our step. */

  while (symbol != NULL)
  {
    if (symbol_is_map(symbol))
    {
      MAP* map;
      symbol_get_map(symbol, &map);
      if (map_has_string(map, "field", name.c_str()))
      {
        /*
         * If one of these files is here the other one has to be
         * here, too. If it's not, it's a bug in the input. */
        if (map_symbol_count(map, filestring) > 0)
        {
          file_info = map;
          break;
        }
      }
    }
    symbol = symbol->next;
  }

  /* No restart files defined? */
  if (symbol == NULL)
  {
    dserror("no restart file definitions found in control file");
  }

#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::FindResultGroup(int step, MAP* file)
{
#ifdef BINIO

  MAP *result_info = NULL;
  MAP *file_info = NULL;

  FindGroup(step,file,"result","result_file",result_info,file_info);
  reader_ = OpenFiles("result_file",file_info);

  restart_step_ = result_info;

#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::FindMeshGroup(int step, MAP* file)
{
#ifdef BINIO

  MAP *result_info = NULL;
  MAP *file_info = NULL;

  FindGroup(step,file,"field","mesh_file",result_info,file_info);
  meshreader_ = OpenFiles("mesh_file",file_info);

  // We do not need result_info as we are interested in the mesh files only.

#endif
}


#ifdef BINIO

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<IO::HDFReader>
IO::DiscretizationReader::OpenFiles(const char* filestring,
                                    MAP* result_step)
{
  int numoutputproc;
  if (!map_find_int(result_step,"num_output_proc",&numoutputproc))
  {
    numoutputproc = 1;
  }

  const std::string name = input_->FileName();

  std::string dirname;
  const std::string::size_type pos = name.find_last_of('/');
  if (pos==std::string::npos)
  {
    dirname = "";
  }
  else
  {
    dirname = name.substr(0,pos+1);
  }

  const std::string filename = map_read_string(result_step, filestring);

  const Epetra_Comm& comm = dis_->Comm();
  Teuchos::RCP<HDFReader> reader = Teuchos::rcp(new HDFReader(dirname));
  reader->Open(filename,numoutputproc,comm.NumProc(),comm.MyPID());
  return reader;
}

#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationWriter::DiscretizationWriter(Teuchos::RCP<DRT::Discretization> dis,
                                               Teuchos::RCP<OutputControl> output)
  :
  dis_(dis),
  step_(0),
  time_(0),
#ifdef BINIO
  meshfile_(-1),
  resultfile_(-1),
  meshfilename_(),
  resultfilename_(),
  meshgroup_(-1),
  resultgroup_(-1),
#endif
  resultfile_changed_(-1),
  meshfile_changed_(-1),
  output_(output)
{
#ifdef BINIO
  Check();
#else
  cerr << "compiled without BINIO: no output will be written\n";
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationWriter::DiscretizationWriter(Teuchos::RCP<DRT::Discretization> dis)
  :
  dis_(dis),
  step_(0),
  time_(0),
#ifdef BINIO
  meshfile_(-1),
  resultfile_(-1),
  meshfilename_(),
  resultfilename_(),
  meshgroup_(-1),
  resultgroup_(-1),
#endif
  resultfile_changed_(-1),
  meshfile_changed_(-1),
  output_(DRT::Problem::Instance()->OutputControlFile())
{
#ifdef BINIO
  Check();
#else
  cerr << "compiled without BINIO: no output will be written\n";
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationWriter::~DiscretizationWriter()
{
#ifdef BINIO
  if (meshfile_ != -1)
  {
    const herr_t status = H5Fclose(meshfile_);
    if (status < 0)
    {
      dserror("Failed to close HDF file %s", meshfilename_.c_str());
    }
  }
  if (resultfile_ != -1)
  {
    // apparently H5Fclose(resultfile_); does not close the file, if there are still
    // some open groups or datasets, so close them first

   	// Get the number of open groups
    int num_og = H5Fget_obj_count(resultfile_, H5F_OBJ_GROUP);

    // get vector to store ids of open groups
    std::vector<hid_t> oid_list(num_og, -1);

    herr_t status = H5Fget_obj_ids(resultfile_,H5F_OBJ_GROUP,num_og,&(oid_list[0]));
    if (status < 0)
       dserror("Failed to get id's of open groups in resultfile");

    // loop over open groups
    for( int i = 0; i < num_og; i++ )
    {
      const herr_t status_g = H5Gclose(oid_list[i]);
          if (status_g < 0)
            dserror("Failed to close HDF-group in resultfile");
    }
    //now close the result file
    const herr_t status_c = H5Fclose(resultfile_);
    if (status_c < 0)
    {
      dserror("Failed to close HDF file %s", resultfilename_.c_str());
    }
  }
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::Check()
{
  // The discretization name will appear in the control file and the post
  // filters expect certain names (for certain types of fields). Thus we
  // restrict the names that can be given to a discretization.
  //
  // If you want to create a dummy discretization choose the name 'none'.
  /*
  const char* names[] = FIELDNAMES;
  for (int i=0; names[i]!=NULL; ++i)
    if (dis_->Name()==names[i])
      return;
  dserror("illegal discretization name '%s'",dis_->Name().c_str());
  */

  // We do not restrict the discretization names anymore. The postfilters
  // do not care about the actual field names.          gjb 04/12
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::CreateMeshFile(const int step)
{
#ifdef BINIO

  std::ostringstream meshname;

  meshname << output_->FileName()
           << ".mesh."
           << dis_->Name()
           << ".s" << step
    ;
  meshfilename_ = meshname.str();
  if (dis_->Comm().NumProc() > 1) {
    meshname << ".p"
             << dis_->Comm().MyPID();
  }

  if (meshfile_ != -1)
  {
    const herr_t status = H5Fclose(meshfile_);
    if (status < 0)
    {
      dserror("Failed to close HDF file %s", meshfilename_.c_str());
    }
  }

  meshfile_ = H5Fcreate(meshname.str().c_str(),
                        H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (meshfile_ < 0)
    dserror("Failed to open file %s", meshname.str().c_str());
  meshfile_changed_ = step;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::CreateResultFile(const int step)
{
#ifdef BINIO

  std::ostringstream resultname;
  resultname << output_->FileName()
             << ".result."
             << dis_->Name()
             << ".s" << step
    ;
  resultfilename_ = resultname.str();
  if (dis_->Comm().NumProc() > 1) {
    resultname << ".p"
               << dis_->Comm().MyPID();
  }
  if (resultfile_ != -1)
  {
    herr_t status = H5Fclose(resultfile_);
    if (status < 0)
    {
      dserror("Failed to close HDF file %s", resultfilename_.c_str());
    }
  }

  // we will never refer to maps stored in other files
  mapcache_.clear();
  mapstack_.clear();

  resultfile_ = H5Fcreate(resultname.str().c_str(),
                          H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (resultfile_ < 0)
    dserror("Failed to open file %s", resultname.str().c_str());
  resultfile_changed_ = step;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::NewResultFile(int numb_run)
{
  resultfile_changed_ = -1;
  meshfile_changed_ = -1;
  output_->NewResultFile(numb_run);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::NewResultFile(std::string name_appendix, int numb_run)
{
  resultfile_changed_ = -1;
  meshfile_changed_ = -1;
  output_->NewResultFile(name_appendix, numb_run);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::OverwriteResultFile()
{
  resultfile_changed_ = -1;
  meshfile_changed_ = -1;
  output_->OverwriteResultFile();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::NewStep(const int step, const double time)
{
#ifdef BINIO

  bool write_file = false;

  step_ = step;
  time_ = time;
  std::ostringstream groupname;
  groupname << "step" << step_;

  if (resultgroup_ != -1)
  {
    const herr_t status = H5Gclose(resultgroup_);
    if (status < 0)
    {
      dserror("Failed to close HDF group in file %s",resultfilename_.c_str());
    }
  }

  if (step_ - resultfile_changed_ >= output_->FileSteps() or
      resultfile_changed_ == -1)
  {
    CreateResultFile(step_);
    write_file = true;
  }

  resultgroup_ = H5Gcreate(resultfile_,groupname.str().c_str(),0);
  if (resultgroup_ < 0)
    dserror("Failed to write HDF-group in resultfile");

  if (dis_->Comm().MyPID() == 0)
  {
    output_->ControlFile()
      << "result:\n"
      << "    field = \"" << dis_->Name() << "\"\n"
      << "    time = " << time << "\n"
      << "    step = " << step << "\n\n";

    if (write_file)
    {
      if (dis_->Comm().NumProc() > 1)
      {
        output_->ControlFile()
          << "    num_output_proc = " << dis_->Comm().NumProc() << "\n";
      }
      std::string filename;
      const std::string::size_type pos = resultfilename_.find_last_of('/');
      if (pos==std::string::npos)
        filename = resultfilename_;
      else
        filename = resultfilename_.substr(pos+1);
      output_->ControlFile()
        << "    result_file = \"" << filename << "\"\n\n";
    }
    output_->ControlFile() << std::flush;
  }
  const herr_t status = H5Fflush(resultgroup_,H5F_SCOPE_LOCAL);
  if (status < 0)
  {
    dserror("Failed to flush HDF file %s", resultfilename_.c_str());
  }
#endif
}

/*----------------------------------------------------------------------*/
/*write double to control file                                  tk 04/08*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteDouble(const std::string name, const double value)
{
#ifdef BINIO

  if (dis_->Comm().MyPID() == 0)
  {
    // using a local stringstream we make sure that we do not change
    // the output formatting of control file permanently
    std::stringstream s;
     s<< "    " << name << " = " << std::scientific << std::setprecision(16)
      << value << "\n\n" << std::flush;
    output_->ControlFile()<<s.str() << std::flush;
  }

#endif
}

/*----------------------------------------------------------------------*/
/*write int to control file                                     tk 04/08*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteInt(const std::string name, const int value)
{
#ifdef BINIO

  if (dis_->Comm().MyPID() == 0)
  {
    output_->ControlFile()
      << "    " << name << " = " << value << "\n\n" << std::flush;
  }
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteVector(const std::string name,
                                           Teuchos::RCP<Epetra_MultiVector> vec,
                                           IO::DiscretizationWriter::VectorType vt)
{
#ifdef BINIO

  std::string valuename = name + ".values";
  double* data = vec->Values();
  const hsize_t size = vec->MyLength() * vec->NumVectors();
  int length = 1;
  if (size==0)
    length = 0;
  const herr_t make_status = H5LTmake_dataset_double(resultgroup_,valuename.c_str(),length,&size,data);
  if (make_status < 0)
    dserror("Failed to create dataset in HDF-resultfile. status=%d", make_status);

  std::string idname;

  // We maintain a map cache to avoid rewriting the same map all the
  // time. The idea is that a map is never modified once it is
  // constructed. Thus the internal data class can be used to find
  // identical maps easily. This will not find all identical maps, but
  // all maps with the same data pointer are guaranteed to be
  // identical.

  std::ostringstream groupname;
  groupname << "/step"
            << step_
            << "/"
    ;

  valuename = groupname.str()+valuename;

  const Epetra_BlockMapData* mapdata = vec->Map().DataPtr();
  std::map<const Epetra_BlockMapData*, std::string>::const_iterator m = mapcache_.find(mapdata);
  if (m!=mapcache_.end())
  {
    // the map has been written already, just link to it again
    idname = m->second;
  }
  else
  {
    const hsize_t mapsize = vec->MyLength();
    idname = name + ".ids";
    int* ids = vec->Map().MyGlobalElements();
    int length = 1;
    if (size==0)
      length = 0;
    const herr_t make_status = H5LTmake_dataset_int(resultgroup_,idname.c_str(),length,&mapsize,ids);
    if (make_status < 0)
      dserror("Failed to create dataset in HDF-resultfile");

    idname = groupname.str()+idname;

    // remember where we put the map
    mapcache_[mapdata] = idname;

    // Make a copy of the map. This is a RCP copy internally. We just make
    // sure here the map stays alive as long as we keep our cache. Otherwise
    // subtle errors could occur.
    mapstack_.push_back(vec->Map());
  }

  if (dis_->Comm().MyPID() == 0)
  {
    std::string vectortype;
    switch (vt)
    {
    case dofvector:
      vectortype = "dof";
      break;
    case nodevector:
      vectortype = "node";
      break;
    case elementvector:
      vectortype = "element";
      break;
    default:
      dserror("unknown vector type %d", vt);
      break;
    }
    output_->ControlFile()
      << "    " << name << ":\n"
      << "        type = \"" << vectortype << "\"\n"
      << "        columns = " << vec->NumVectors() << "\n"
      << "        values = \"" << valuename.c_str() << "\"\n"
      << "        ids = \"" << idname.c_str() << "\"\n\n"  // different names + other informations?
      << std::flush;
  }
  const herr_t flush_status = H5Fflush(resultgroup_,H5F_SCOPE_LOCAL);
  if (flush_status < 0)
  {
    dserror("Failed to flush HDF file %s", resultfilename_.c_str());
  }

#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteVector(const std::string name,
                                           const std::vector<char>& vec,
                                           const Epetra_Map& elemap,
                                           IO::DiscretizationWriter::VectorType vt)
{
#ifdef BINIO

  std::string valuename = name + ".values";
  const hsize_t size = vec.size();
  const char* data = &vec[0];
  const herr_t make_status = H5LTmake_dataset_char(resultgroup_,valuename.c_str(),1,&size,data);
  if (make_status < 0)
    dserror("Failed to create dataset in HDF-resultfile");

  std::string idname;

  // We maintain a map cache to avoid rewriting the same map all the
  // time. The idea is that a map is never modified once it is
  // constructed. Thus the internal data class can be used to find
  // identical maps easily. This will not find all identical maps, but
  // all maps with the same data pointer are guaranteed to be
  // identical.

  std::ostringstream groupname;
  groupname << "/step"
            << step_
            << "/"
    ;

  valuename = groupname.str()+valuename;

  const Epetra_BlockMapData* mapdata = elemap.DataPtr();
  std::map<const Epetra_BlockMapData*, std::string>::const_iterator m = mapcache_.find(mapdata);
  if (m!=mapcache_.end())
  {
    // the map has been written already, just link to it again
    idname = m->second;
  }
  else
  {
    const hsize_t mapsize = elemap.NumMyElements();
    idname = name + ".ids";
    int* ids = elemap.MyGlobalElements();
    const herr_t make_status = H5LTmake_dataset_int(resultgroup_,idname.c_str(),1,&mapsize,ids);
    if (make_status < 0)
      dserror("Failed to create dataset in HDF-resultfile");

    idname = groupname.str()+idname;

    // remember where we put the map
    mapcache_[mapdata] = idname;

    // Make a copy of the map. This is a RCP copy internally. We just make
    // sure here the map stays alive as long as we keep our cache. Otherwise
    // subtle errors could occur.
    mapstack_.push_back(elemap);
  }

  if (dis_->Comm().MyPID() == 0)
  {
    std::string vectortype;
    switch (vt)
    {
    case dofvector:
      vectortype = "dof";
      break;
    case nodevector:
      vectortype = "node";
      break;
    case elementvector:
      vectortype = "element";
      break;
    default:
      dserror("unknown vector type %d", vt);
      break;
    }
    output_->ControlFile()
      << "    " << name << ":\n"
      << "        type = \"" << vectortype << "\"\n"
      << "        columns = 1\n"
      << "        values = \"" << valuename << "\"\n"
      << "        ids = \"" << idname << "\"\n\n" // different names + other informations?
      << std::flush;
  }
  const herr_t flush_status = H5Fflush(resultgroup_,H5F_SCOPE_LOCAL);
  if (flush_status < 0)
  {
    dserror("Failed to flush HDF file %s", resultfilename_.c_str());
  }
#endif
}


/*----------------------------------------------------------------------*
 *                                                          a.ger 11/07 *
 *----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteCondition(const std::string condname) const
{
#ifdef BINIO
  // put condition into block
  Teuchos::RCP<std::vector<char> > block = dis_->PackCondition(condname);

  // write block to file. Note: Block can be empty, if the condition is not found,
  // which means it is not used -> so no dserror() here
  if(!block->empty())
  {
    if (dis_->Comm().MyPID() == 0)
    {
      output_->ControlFile() << "    condition = \"" << condname << "\"\n";
    }

    hsize_t dim = static_cast<hsize_t>(block->size());
    const herr_t status = H5LTmake_dataset_char(
            meshgroup_,
            condname.c_str(),
            1,
            &dim,
            &((*block)[0]));
    if (status < 0)
      dserror("Failed to create dataset in HDF-meshfile");
  }
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteMesh(const int step, const double time)
{
#ifdef BINIO


  bool write_file = false;

  if (step - meshfile_changed_ >= output_->FileSteps() or
      meshfile_changed_ == -1)
  {
    CreateMeshFile(step);
    write_file = true;
  }
  std::ostringstream name;
  name << "step" << step;
  meshgroup_ = H5Gcreate(meshfile_,name.str().c_str(),0);
  if (meshgroup_ < 0)
    dserror("Failed to write group in HDF-meshfile");

  // procs without row elements must not write data
  Teuchos::RCP<std::vector<char> > elementdata = dis_->PackMyElements();
  hsize_t dim = static_cast<hsize_t>(elementdata->size());
  if (elementdata->size() != 0)
  {
    const herr_t element_status = H5LTmake_dataset_char(meshgroup_,"elements",1,&dim,&((*elementdata)[0]));
    if (element_status < 0)
      dserror("Failed to create dataset in HDF-meshfile for group none");
  }
  else
  {
    const herr_t element_status = H5LTmake_dataset_char(meshgroup_,"elements",0,&dim,&((*elementdata)[0]));
    if (element_status < 0)
      dserror("Failed to create dataset in HDF-meshfile");
  }

  // procs without row nodes must not write data
  Teuchos::RCP<std::vector<char> > nodedata = dis_->PackMyNodes();
  dim = static_cast<hsize_t>(nodedata->size());
  if (nodedata->size() != 0)
  {
    const herr_t node_status = H5LTmake_dataset_char(meshgroup_,"nodes",1,&dim,&((*nodedata)[0]));
    if (node_status < 0)
      dserror("Failed to create dataset in HDF-meshfile");
  }
  else
  {
    const herr_t node_status = H5LTmake_dataset_char(meshgroup_,"nodes",0,&dim,&((*nodedata)[0]));
    if (node_status < 0)
      dserror("Failed to create dataset in HDF-meshfile");
  }

  // ... write other mesh informations
  if (dis_->Comm().MyPID() == 0)
  {
    output_->ControlFile()
      << "field:\n"
      << "    field = \"" << dis_->Name() << "\"\n"
      << "    time = " << time << "\n"
      << "    step = " << step << "\n\n"
      << "    num_nd = " << dis_->NumGlobalNodes() << "\n"
      << "    num_ele = " << dis_->NumGlobalElements() << "\n"
      << "    num_dof = " << dis_->DofRowMap(0)->NumGlobalElements() << "\n\n"
      ;

    WriteCondition("SurfacePeriodic");
    WriteCondition("LinePeriodic");

    // knotvectors for nurbs-discretisation
    WriteKnotvector();

    if (write_file)
    {
      if (dis_->Comm().NumProc() > 1)
      {
        output_->ControlFile()
          << "    num_output_proc = " << dis_->Comm().NumProc() << "\n";
      }
      std::string filename;
      std::string::size_type pos = meshfilename_.find_last_of('/');
      if (pos==std::string::npos)
        filename = meshfilename_;
      else
        filename = meshfilename_.substr(pos+1);
      output_->ControlFile()
        << "    mesh_file = \"" << filename << "\"\n\n";
    }
    output_->ControlFile() << std::flush;
  }
  const herr_t flush_status = H5Fflush(meshgroup_,H5F_SCOPE_LOCAL);
  if (flush_status < 0)
  {
    dserror("Failed to flush HDF file %s", meshfilename_.c_str());
  }
  const herr_t close_status = H5Gclose(meshgroup_);
  if (close_status < 0)
  {
    dserror("Failed to close HDF group in file %s", meshfilename_.c_str());
  }
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteMesh(const int step, const double time, std::string name_base_file )
{
#ifdef BINIO

  // ... write other mesh informations
  if (dis_->Comm().MyPID() == 0)
  {
    output_->ControlFile()
      << "field:\n"
      << "    field = \"" << dis_->Name() << "\"\n"
      << "    time = " << time << "\n"
      << "    step = " << step << "\n\n"
      << "    num_nd = " << dis_->NumGlobalNodes() << "\n"
      << "    num_ele = " << dis_->NumGlobalElements() << "\n"
      << "    num_dof = " << dis_->DofRowMap(0)->NumGlobalElements() << "\n\n"
      ;

    //WriteCondition("SurfacePeriodic");
    //WriteCondition("LinePeriodic");

    // knotvectors for nurbs-discretisation
    //WriteKnotvector();
    // create name for meshfile as in createmeshfile which is not called here
    std::ostringstream meshname;

    meshname << name_base_file
               << ".mesh."
               << dis_->Name()
               << ".s" << step
        ;
      meshfilename_ = meshname.str();

      if (dis_->Comm().NumProc() > 1)
      {
        output_->ControlFile()
          << "    num_output_proc = " << dis_->Comm().NumProc() << "\n";
      }
      std::string filename;
      std::string::size_type pos = meshfilename_.find_last_of('/');
      if (pos==std::string::npos)
        filename = meshfilename_;
      else
        filename = meshfilename_.substr(pos+1);
      output_->ControlFile()
        << "    mesh_file = \"" << filename << "\"\n\n";
       // << "    mesh_file = \"" << name_base_file << "\"\n\n";

    output_->ControlFile() << std::flush;
  }

#endif
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteElementData()
{
#ifdef BINIO
  std::map<std::string,int>::const_iterator fool;
  std::map<std::string,int> names;   // contains name and dimension of data

  // loop all elements and build map of data names and dimensions
  const Epetra_Map* elerowmap = dis_->ElementRowMap();
  for (int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    // get names and dimensions from every element
    dis_->lRowElement(i)->VisNames(names);
  }
  
  // By applying GatherAll we get the combined map including all elemental values
  // which where found by VisNames
  const Epetra_Comm& comm = dis_->Comm();
  LINALG::GatherAll(names, comm);

  // make sure there's no name with a dimension of less than 1
  for (fool = names.begin(); fool!= names.end(); ++fool)
    if (fool->second<1) dserror("Dimension of data must be at least 1");

  // loop all names aquired form the elements and fill data vectors
  std::vector<double> eledata(0);
  for (fool = names.begin(); fool!= names.end(); ++fool)
  {
    const int dimension = fool->second;
    eledata.resize(dimension);
    for (int i=0; i<dimension; ++i) eledata[i] = 0.0;

    // MultiVector stuff from the elements is put in
    Epetra_MultiVector sysdata(*elerowmap,dimension,true);

    for (int i=0; i<elerowmap->NumMyElements(); ++i)
    {
      // get data for a given name from element & put in sysdata
      dis_->lRowElement(i)->VisData(fool->first,eledata);
      if ((int)eledata.size() != dimension)
        dserror("element manipulated size of visualization data");
      for (int j=0; j<dimension; ++j) (*sysdata(j))[i] = eledata[j];
    }

    WriteVector(fool->first, Teuchos::rcp(&sysdata,false), IO::DiscretizationWriter::elementvector);

  } // for (fool = names.begin(); fool!= names.end(); ++fool)

#endif
}


/*----------------------------------------------------------------------*
 *                                                          gammi 05/08 *
 *----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteKnotvector() const
{
#ifdef BINIO
  // try a dynamic cast of the discretisation to a nurbs discretisation
  DRT::NURBS::NurbsDiscretization* nurbsdis
    =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis_));

  if(nurbsdis!=NULL)
  {
    // get knotvector from nurbsdis
    Teuchos::RCP<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

    // put knotvector into block
    DRT::PackBuffer block;
    knots->Pack(block);
    block.StartPacking();
    knots->Pack(block);

    // write block to file
    if(!block().empty())
    {
      hsize_t dim = static_cast<hsize_t>(block().size());
      const herr_t status = H5LTmake_dataset_char(
        meshgroup_,
        "knotvector",
        1,
        &dim,
        &(block()[0]));
      if (status < 0)
        dserror("Failed to create dataset in HDF-meshfile");
    }
    else
    {
      dserror("block empty --- couldn't write knots\n");
    }
  }
#endif

  return;
}

/*----------------------------------------------------------------------*/
/* write a stl vector of doubles from proc0                             */
/*----------------------------------------------------------------------*/
  void IO::DiscretizationWriter::WriteRedundantDoubleVector(const std::string name,
							    Teuchos::RCP<std::vector<double> > doublevec)
{
#ifdef BINIO
  if (dis_->Comm().MyPID() == 0)
  {
    // only proc0 writes the vector entities to the binary data
    // an appropriate name has to be provided
    std::string valuename = name + ".values";
    const hsize_t size = doublevec->size();
    int length = 1;
    if (size==0)
      length = 0;
    const herr_t make_status = H5LTmake_dataset_double(resultgroup_,valuename.c_str(),length,&size,&((*doublevec)[0]));
    if (make_status < 0)
      dserror("Failed to create dataset in HDF-resultfile. status=%d", make_status);

    // do I need the following naming stuff?
    std::ostringstream groupname;

    groupname << "/step"
	      << step_
	      << "/"
      ;

    valuename = groupname.str()+valuename;

    // a comment is also added to the control file
    output_->ControlFile()
      << "    " << name << ":\n"
      << "        values = \"" << valuename.c_str() << "\"\n\n"
      << std::flush;

    const herr_t flush_status = H5Fflush(resultgroup_,H5F_SCOPE_LOCAL);
    if (flush_status < 0)
    {
      dserror("Failed to flush HDF file %s", resultfilename_.c_str());
    }
  } // endif proc0

#endif
}
