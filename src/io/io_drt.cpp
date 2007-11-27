/*----------------------------------------------------------------------*/
/*!
\file io_drt.cpp

\brief output context of one discretization

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <iostream>
#include <sstream>
#include <string>

#include "io_drt.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_globalproblem.H"

using namespace std;

#ifdef BINIO

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for output.

  This structure needs to be initialized at startup. The whole output
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern BIN_OUT_MAIN bin_out_main;

/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for input.

  This structure needs to be initialized at startup. The whole input
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern BIN_IN_MAIN bin_in_main;

/*----------------------------------------------------------------------*/
/*!
  \brief All fields names.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern CHAR* fieldnames[];

#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
static void FindPosition(RefCountPtr<DRT::Discretization> dis, int& field_pos, int& disnum, int probnum=0)
{
#ifdef BINIO

  for (field_pos=0; field_pos<genprob.numfld; ++field_pos)
  {
    for (disnum=0; disnum<field[field_pos].ndis; ++disnum)
    {
      if (DRT::Problem::Instance(probnum)->Dis(field_pos,disnum).get() == dis.get())
      {
        return;
      }
    }
  }
  // no field found
  dserror("unregistered field object");
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationReader::DiscretizationReader(RefCountPtr<DRT::Discretization> dis, int step)
  : dis_(dis)
{
  restart_step_ = FindResultGroup(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::ReadVector(RefCountPtr<Epetra_Vector> vec, string name)
{
#ifdef BINIO
  MAP* result = map_read_map(restart_step_, const_cast<char*>(name.c_str()));
  string id_path = map_read_string(result, "ids");
  string value_path = map_read_string(result, "values");
  RefCountPtr<Epetra_Vector> nv = reader_->ReadResultData(id_path, value_path, dis_->Comm());
  LINALG::Export(*nv, *vec);
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::ReadMesh(int step)
{
#ifdef BINIO
  FindMeshGroup(step);

  RefCountPtr<vector<char> > nodedata =
    meshreader_->ReadNodeData(step,dis_->Comm().NumProc(),dis_->Comm().MyPID());

  RefCountPtr<vector<char> > elementdata =
   meshreader_->ReadElementData(step,dis_->Comm().NumProc(),dis_->Comm().MyPID());

  // before we unpack nodes/elements we store a copy of the nodal row/col map
  RefCountPtr<Epetra_Map> noderowmap = rcp(new Epetra_Map(*dis_->NodeRowMap()));
  RefCountPtr<Epetra_Map> nodecolmap = rcp(new Epetra_Map(*dis_->NodeColMap()));

  // unpack nodes and elements and redistirbuted to current layout
  dis_->UnPackMyNodes(nodedata);
  dis_->UnPackMyElements(elementdata);
  dis_->Redistribute(*noderowmap,*nodecolmap);
  int err = dis_->FillComplete();
  if (err) dserror("FillComplete() returned %d",err);

#endif
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IO::DiscretizationReader::ReadInt(string name)
{
  return map_read_int(restart_step_, const_cast<char*>(name.c_str()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double IO::DiscretizationReader::ReadDouble(string name)
{
  return map_read_real(restart_step_, const_cast<char*>(name.c_str()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAP *IO::DiscretizationReader::FindResultGroup(int step)
{
#ifdef BINIO

  MAP *result_info = NULL;
  SYMBOL *symbol;
  int field_pos;
  int disnum;

  FindPosition(dis_, field_pos, disnum);

  FIELD *actfield = &(field[field_pos]);

  /*
   * Iterate all symbols under the name "result" and get the one that
   * matches the given step. Note that this iteration starts from the
   * last result group and goes backward. */

  symbol = map_find_symbol(&(bin_in_main.table), "result");
  while (symbol != NULL)
  {
    if (symbol_is_map(symbol))
    {
      MAP* map;
      symbol_get_map(symbol, &map);
      if (map_has_string(map, "field", fieldnames[actfield->fieldtyp]) &&
          map_has_int(map, "field_pos", field_pos) &&
          map_has_int(map, "discretization", disnum) &&
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
    dserror("No restart entry for step %d in symbol table. Control file corrupt?",step);
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
      if (map_has_string(map, "field", fieldnames[actfield->fieldtyp]) &&
          map_has_int(map, "field_pos", field_pos) &&
          map_has_int(map, "discretization", disnum))
      {
        /*
         * If one of these files is here the other one has to be
         * here, too. If it's not, it's a bug in the input. */
        if (map_symbol_count(map, "result_file") > 0)
        {
          OpenDataFiles(map);
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

  return result_info;
#else
  return NULL;
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAP *IO::DiscretizationReader::FindMeshGroup(int step)
{
#ifdef BINIO

  MAP *result_info = NULL;
  SYMBOL *symbol;
  int field_pos;
  int disnum;

  FindPosition(dis_, field_pos, disnum);

  FIELD *actfield = &(field[field_pos]);

  /*
   * Iterate all symbols under the name "result" and get the one that
   * matches the given step. Note that this iteration starts from the
   * last result group and goes backward. */

  symbol = map_find_symbol(&(bin_in_main.table), "field");
  while (symbol != NULL)
  {
    if (symbol_is_map(symbol))
    {
      MAP* map;
      symbol_get_map(symbol, &map);
      if (map_has_string(map, "field", fieldnames[actfield->fieldtyp]) &&
          map_has_int(map, "field_pos", field_pos) &&
          map_has_int(map, "discretization", disnum) &&
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
    dserror("No restart entry for step %d in symbol table. Control file corrupt?",step);
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
      if (map_has_string(map, "field", fieldnames[actfield->fieldtyp]) &&
          map_has_int(map, "field_pos", field_pos) &&
          map_has_int(map, "discretization", disnum))
      {
        /*
         * If one of these files is here the other one has to be
         * here, too. If it's not, it's a bug in the input. */
        if (map_symbol_count(map, "mesh_file") > 0)
        {
          OpenMeshFiles(map);
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

  return result_info;
#else
  return NULL;
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::OpenDataFiles(MAP* result_step)
{
#ifdef BINIO
  int numoutputproc;
  if (!map_find_int(result_step,"num_output_proc",&numoutputproc))
  {
    numoutputproc = 1;
  }

  string name = bin_out_main.name;

  string dirname;
  string::size_type pos = name.find_last_of('/');
  if (pos==string::npos)
  {
    dirname = "";
  }
  else
  {
    dirname = name.substr(0,pos+1);
  }

  string filename;
  filename = map_read_string(result_step, "result_file");

  reader_ = rcp(new HDFReader(dirname));
  reader_->Open(filename,numoutputproc);
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationReader::OpenMeshFiles(MAP* result_step)
{
#ifdef BINIO
  int numoutputproc;
  if (!map_find_int(result_step,"num_output_proc",&numoutputproc))
  {
    numoutputproc = 1;
  }

  string name = bin_out_main.name;

  string dirname;
  string::size_type pos = name.find_last_of('/');
  if (pos==string::npos)
  {
    dirname = "";
  }
  else
  {
    dirname = name.substr(0,pos+1);
  }

  string filename;
  filename = map_read_string(result_step, "mesh_file");

  meshreader_ = rcp(new HDFReader(dirname));
  meshreader_->Open(filename,numoutputproc);
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationWriter::DiscretizationWriter(RefCountPtr<DRT::Discretization> dis, int probnum):
  dis_(dis),
  disnum_(0),
  field_pos_(0),
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
  meshfile_changed_(-1)

{
#ifdef BINIO
  cfname_ = bin_out_main.name;
  cf_ = bin_out_main.control_file;
  steps_per_file_ = bin_out_main.steps_per_file;
#endif

  FindPosition(dis_, field_pos_, disnum_, probnum);
#ifndef BINIO
  cerr << "compiled without BINIO: no output will be written\n";
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IO::DiscretizationWriter::~DiscretizationWriter()
{
#ifdef BINIO
  herr_t status;
  if (meshfile_ != -1)
  {
    status = H5Fclose(meshfile_);
    if (status < 0)
    {
      dserror("Failed to close HDF file %s", meshfilename_.c_str());
    }
  }
  if (resultfile_ != -1)
  {
    status = H5Fclose(resultfile_);
    if (status < 0)
    {
      dserror("Failed to close HDF file %s", resultfilename_.c_str());
    }
  }
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::CreateMeshFile(int step)
{
#ifdef BINIO

  ostringstream meshname;

  meshname << cfname_
           << ".mesh."
           << fieldnames[field[field_pos_].fieldtyp]
           << ".f" << field_pos_
           << ".d" << disnum_
           << ".s" << step
    ;
  meshfilename_ = meshname.str();
  if (dis_->Comm().NumProc() > 1) {
    meshname << ".p"
             << dis_->Comm().MyPID();
  }

  if (meshfile_ != -1)
  {
    herr_t status = H5Fclose(meshfile_);
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
void IO::DiscretizationWriter::CreateResultFile(int step)
{
#ifdef BINIO

  ostringstream resultname;
  resultname << cfname_
             << ".result."
             << fieldnames[field[field_pos_].fieldtyp]
             << ".f" << field_pos_
             << ".d" << disnum_
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

  resultfile_ = H5Fcreate(resultname.str().c_str(),
                          H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (resultfile_ < 0)
    dserror("Failed to open file %s", resultname.str().c_str());
  resultfile_changed_ = step;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::NewStep(int step, double time)
{
#ifdef BINIO

  herr_t status;
  bool write_file = false;

  step_ = step;
  time_ = time;
  ostringstream groupname;
  groupname << "step" << step_;

  if (resultgroup_ != -1)
  {
    status = H5Gclose(resultgroup_);
    if (status < 0)
    {
      dserror("Failed to close HDF group in file %s",resultfilename_.c_str());
    }
  }

  if (step_ - resultfile_changed_ >= steps_per_file_
      || resultfile_changed_ == -1)
  {
    CreateResultFile(step_);
    write_file = true;
  }

  resultgroup_ = H5Gcreate(resultfile_,groupname.str().c_str(),0);
  if (resultgroup_ < 0)
    dserror("Failed to write HDF-group in resultfile");

  if (dis_->Comm().MyPID() == 0)
  {
    fprintf(cf_,
            "result:\n"
            "    field = \"%s\"\n"
            "    field_pos = %i\n"
            "    discretization = %i\n"
            "    time = %f\n"
            "    step = %i\n\n",
            fieldnames[field[field_pos_].fieldtyp],
            field_pos_,
            disnum_,
            time,
            step
      );
    if (write_file)
    {
      if (dis_->Comm().NumProc() > 1)
      {
        fprintf(cf_,
                "    num_output_proc = %d\n",
                dis_->Comm().NumProc());
      }
      string filename;
      string::size_type pos = resultfilename_.find_last_of('/');
      if (pos==string::npos)
        filename = resultfilename_;
      else
        filename = resultfilename_.substr(pos+1);
      fprintf(cf_,
              "    result_file = \"%s\"\n\n",
              filename.c_str()
        );
    }
    fflush(cf_);
  }
  status = H5Fflush(resultgroup_,H5F_SCOPE_LOCAL);
  if (status < 0)
  {
    dserror("Failed to flush HDF file %s", resultfilename_.c_str());
  }
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteVector(string name, RefCountPtr<Epetra_Vector> vec)
{
#ifdef BINIO

  herr_t status;
  string valuename = name + ".values";
  double* data = vec->Values();
  hsize_t size = vec->MyLength();
  status = H5LTmake_dataset_double(resultgroup_,valuename.c_str(),1,&size,data);
  if (status < 0)
    dserror("Failed to create dataset in HDF-resultfile");

  string idname;

  // We maintain a map cache to avoid rewriting the same map all the
  // time. The idea is that a map is never modified once it is
  // constructed. Thus the internal data class can be used to find
  // identical maps easily. This will not find all identical maps, but
  // all maps with the same data pointer are guaranteed to be
  // identical.

  ostringstream groupname;
  groupname << "/step"
            << step_
            << "/"
    ;

  valuename = groupname.str()+valuename;

  const Epetra_BlockMapData* mapdata = vec->Map().DataPtr();
  std::map<const Epetra_BlockMapData*, std::string>::iterator m = mapcache_.find(mapdata);
  if (m!=mapcache_.end())
  {
    // the map has been written already, just link to it again
    idname = m->second;
  }
  else
  {
    idname = name + ".ids";
    int* ids = vec->Map().MyGlobalElements();
    status = H5LTmake_dataset_int(resultgroup_,idname.c_str(),1,&size,ids);
    if (status < 0)
      dserror("Failed to create dataset in HDF-resultfile");

    idname = groupname.str()+idname;

    // remember where we put the map
    mapcache_[mapdata] = idname;
  }

  if (dis_->Comm().MyPID() == 0)
  {
    fprintf(cf_,
            "    %s:\n"
            "        values = \"%s\"\n"
            "        ids = \"%s\"\n\n",  // different names + other informations?
            name.c_str(),
            valuename.c_str(),
            idname.c_str()
      );
    fflush(cf_);
  }
  status = H5Fflush(resultgroup_,H5F_SCOPE_LOCAL);
  if (status < 0)
  {
    dserror("Failed to flush HDF file %s", resultfilename_.c_str());
  }
#endif
}

// write a specific condition
void IO::DiscretizationWriter::WriteCondition(const string condname)
{
  //--------------------------------------------------
  // pack boundary condition
  RefCountPtr<vector<char> > block = rcp(new vector<char>);

  // get boundary conditions
  vector<DRT::Condition*> cond;
  dis_->GetCondition(condname,cond);

  for (vector<DRT::Condition*>::const_iterator i = cond.begin();
       i!=cond.end();
       ++i)
  {
    vector<char> conddata;
    (*i)->Pack(conddata);
    DRT::ParObject::AddtoPack(*block,conddata);
  }

  //--------------------------------------------------
  // write block to file
  if(!block->empty())
  {
    hsize_t dim[] = {static_cast<hsize_t>(block->size())};
    const herr_t status = H5LTmake_dataset_char(
            meshgroup_,
            condname.c_str(),
            1,
            dim,
            &((*block)[0]));
    if (status < 0)
      dserror("Failed to create dataset in HDF-meshfile");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IO::DiscretizationWriter::WriteMesh(int step, double time)
{
#ifdef BINIO

  bool write_file = false;

  if (step - meshfile_changed_ >= steps_per_file_
    || meshfile_changed_ == -1)
  {
    CreateMeshFile(step);
    write_file = true;
  }
  ostringstream name;
  name << "step" << step;
  meshgroup_ = H5Gcreate(meshfile_,name.str().c_str(),0);
  if (meshgroup_ < 0)
    dserror("Failed to write group in HDF-meshfile");

  RefCountPtr<vector<char> > elementdata = dis_->PackMyElements();
  if (elementdata->size()==0)
    dserror("no element data no proc %d. Too few elements?", dis_->Comm().MyPID());
  hsize_t dim[] = {static_cast<hsize_t>(elementdata->size())};
  const herr_t element_status = H5LTmake_dataset_char(meshgroup_,"elements",1,dim,&((*elementdata)[0]));
  if (element_status < 0)
    dserror("Failed to create dataset in HDF-meshfile");

  RefCountPtr<vector<char> > nodedata = dis_->PackMyNodes();
  dim[0] = static_cast<hsize_t>(nodedata->size());
  const herr_t node_status = H5LTmake_dataset_char(meshgroup_,"nodes",1,dim,&((*nodedata)[0]));
  if (node_status < 0)
    dserror("Failed to create dataset in HDF-meshfile");

  // ... write other mesh informations

  if (dis_->Comm().MyPID() == 0)
  {
#if 1
    WriteCondition("SurfacePeriodic");
    WriteCondition("LinePeriodic");
      
//    {
//      //--------------------------------------------------
//      // pack all periodic boundary conditions
//      RefCountPtr<vector<char> > pbcblock = rcp(new vector<char>);
//      // get all periodic boundary conditions
//      vector<DRT::Condition*> percond;
//
//      // first the surfaces
//      dis_->GetCondition("SurfacePeriodic",percond);
//
//      for (vector<DRT::Condition*>::iterator i=percond.begin();
//           i!=percond.end();
//           ++i)
//      {
//        vector<char> perconddata;
//        (*i)->Pack(perconddata);
//        DRT::ParObject::AddtoPack(*pbcblock,perconddata);
//      }
//
//      // then the lines
//      dis_->GetCondition("LinePeriodic",percond);
//
//      for (vector<DRT::Condition*>::iterator i=percond.begin();
//           i!=percond.end();
//           ++i)
//      {
//        vector<char> perconddata;
//        (*i)->Pack(perconddata);
//        DRT::ParObject::AddtoPack(*pbcblock,perconddata);
//      }
//
//
//      //--------------------------------------------------
//      // write pbcblock to file
//      if(!pbcblock->empty())
//      {
//        dim[0] = static_cast<hsize_t>(pbcblock->size());
//        const herr_t bc_status = H5LTmake_dataset_char(meshgroup_,
//                                       "periodicbc",
//                                       1,
//                                       dim,
//                                       &((*pbcblock)[0]));
//        if (bc_status < 0)
//          dserror("Failed to create dataset in HDF-meshfile");
//      }
//    }
#endif
    fprintf(cf_,
            "field:\n"
            "    field = \"%s\"\n"
            "    field_pos = %i\n"
            "    discretization = %i\n"
            "    dis_name = \"%s\"\n\n"
            "    step = %i\n"
            "    time = %f\n\n"
            "    num_nd = %i\n"
            "    num_ele = %i\n"
            "    num_dof = %i\n\n",
            fieldnames[field[field_pos_].fieldtyp],
            field_pos_,
            disnum_,
            dis_->Name().c_str(),
            step,
            time,
            dis_->NumGlobalNodes(),
            dis_->NumGlobalElements(),
            dis_->DofRowMap()->NumGlobalElements()  // is that right ???

      );
    if (write_file)
    {
      if (dis_->Comm().NumProc() > 1)
      {
        fprintf(cf_,
                "    num_output_proc = %d\n",
                dis_->Comm().NumProc());
      }
      string filename;
      string::size_type pos = meshfilename_.find_last_of('/');
      if (pos==string::npos)
        filename = meshfilename_;
      else
        filename = meshfilename_.substr(pos+1);
      fprintf(cf_,
              "    mesh_file = \"%s\"\n\n",
              filename.c_str());
    }
    fflush(cf_);
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

#endif
