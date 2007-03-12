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
#ifdef TRILINOS_PACKAGE

#include <iostream>
#include <sstream>
#include <string>

#include "io_drt.H"

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
  \brief All fields names.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern CHAR* fieldnames[];

#endif

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DiscretizationWriter::DiscretizationWriter(RefCountPtr<DRT::Discretization> dis):
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
  FindPosition(field_pos_,disnum_);
#ifndef BINIO
  cerr << "compiled without BINIO: no output will be written\n";
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DiscretizationWriter::~DiscretizationWriter()
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
void DiscretizationWriter::CreateMeshFile(int step)
{
#ifdef BINIO

  ostringstream meshname;

  meshname << bin_out_main.name
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
void DiscretizationWriter::CreateResultFile(int step)
{
#ifdef BINIO

  ostringstream resultname;
  resultname << bin_out_main.name
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

  resultfile_ = H5Fcreate(resultname.str().c_str(),
                          H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if (resultfile_ < 0)
    dserror("Failed to open file %s", resultname.str().c_str());
  resultfile_changed_ = step;
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DiscretizationWriter::FindPosition(int& field_pos, unsigned int& disnum)
{
#ifdef BINIO

  for (field_pos=0;field_pos<genprob.numfld; ++field_pos)
  {
    vector<RefCountPtr<DRT::Discretization> >* discretizations =
      static_cast<vector<RefCountPtr<DRT::Discretization> >*>
      (field[field_pos].ccadis);
    for (disnum=0;disnum<discretizations->size();++disnum)
    {
      if ((*discretizations)[disnum] == dis_) {
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
void DiscretizationWriter::NewStep(int step, double time)
{
#ifdef BINIO

  herr_t status;
  bool write_file = false;

  step_ = step;
  time_ = time;
  ostringstream groupname;
  groupname << "step" << step_;

  if (resultgroup_ != -1) {
    status = H5Gclose(resultgroup_);
    if (status < 0)
    {
      dserror("Failed to close HDF group in file %s",resultfilename_.c_str());
    }
  }

  if (step_ - resultfile_changed_ >= bin_out_main.steps_per_file
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
    fprintf(bin_out_main.control_file,
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
    if (write_file) {
      if (dis_->Comm().NumProc() > 1)
        fprintf(bin_out_main.control_file,
                "    num_output_proc = %d\n",
                dis_->Comm().NumProc());
      fprintf(bin_out_main.control_file,
              "    result_file = \"%s\"\n\n",
              resultfilename_.c_str()
        );
    }
    fflush(bin_out_main.control_file);
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
void DiscretizationWriter::WriteVector(string name, RefCountPtr<Epetra_Vector> vec)
{
#ifdef BINIO

  herr_t status;
  string valuename = name;
  valuename.append(".values");
  double* data = vec->Values();
  hsize_t size = vec->MyLength();
  status = H5LTmake_dataset_double(resultgroup_,valuename.c_str(),1,&size,data);
  if (status < 0)
    dserror("Failed to create dataset in HDF-resultfile");
  string idname = name;
  idname.append(".ids");
  int* ids = vec->Map().MyGlobalElements();
  status = H5LTmake_dataset_int(resultgroup_,idname.c_str(),1,&size,ids);
  if (status < 0)
    dserror("Failed to create dataset in HDF-resultfile");

  if (dis_->Comm().MyPID() == 0)
  {
    ostringstream groupname;
    groupname << "/step"
              << step_
              << "/"
      ;

    valuename = groupname.str()+valuename;
    idname = groupname.str()+idname;
    fprintf(bin_out_main.control_file,
            "    %s:\n"
            "        values = \"%s\"\n"
            "        ids = \"%s\"\n\n",  // different names + other informations?
            name.c_str(),
            valuename.c_str(),
            idname.c_str()
      );
    fflush(bin_out_main.control_file);
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
void DiscretizationWriter::WriteMesh(int step, double time)
{
#ifdef BINIO

  herr_t status;
  bool write_file = false;

  if (step - meshfile_changed_ >= bin_out_main.steps_per_file
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
  hsize_t dim[] = {static_cast<hsize_t>(elementdata->size())};
  status = H5LTmake_dataset_char(meshgroup_,"elements",1,dim,&((*elementdata)[0]));
  if (status < 0)
    dserror("Failed to create dataset in HDF-meshfile");

  RefCountPtr<vector<char> > nodedata = dis_->PackMyNodes();
  dim[0] = static_cast<hsize_t>(nodedata->size());
  status = H5LTmake_dataset_char(meshgroup_,"nodes",1,dim,&((*nodedata)[0]));
  if (status < 0)
    dserror("Failed to create dataset in HDF-meshfile");

  // ... write other mesh informations

  if (dis_->Comm().MyPID() == 0)
  {
    fprintf(bin_out_main.control_file,
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
        fprintf(bin_out_main.control_file,
                "    num_output_proc = %d\n",
                dis_->Comm().NumProc());
      }
      fprintf(bin_out_main.control_file,
              "    mesh_file = \"%s\"\n\n",
              meshfilename_.c_str());
    }
    fflush(bin_out_main.control_file);
  }
  status = H5Fflush(meshgroup_,H5F_SCOPE_LOCAL);
  if (status < 0)
  {
    dserror("Failed to flush HDF file %s", meshfilename_.c_str());
  }
  status = H5Gclose(meshgroup_);
  if (status < 0)
  {
    dserror("Failed to close HDF group in file %s", meshfilename_.c_str());
  }
#endif
}

#endif
#endif
