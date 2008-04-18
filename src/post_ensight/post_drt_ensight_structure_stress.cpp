/*!
  \file post_drt_ensight_structure_stress.cpp

  \brief postprocessing of structural stresses

  <pre>
  Maintainer: Lena Wiechert
              wiechert@lnm.mw.tum.de
              http://www.lnm.mw.tum.de/Members/wiechert
              089-289-15303
  </pre>

*/

#ifdef CCADISCRET

#include "post_drt_ensight_writer.H"
#include <string>
#include "post_drt_ensight_single_field_writers.H"
#include "../drt_lib/drt_utils.H"

using namespace std;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::PostStress(const string groupname, const string stresstype)
{
  PostResult result = PostResult(field_);
  result.next_result();

  if (!map_has_map(result.group(), groupname.c_str()))
    return;

  //--------------------------------------------------------------------
  // calculation and output of nodal stresses in xyz-reference frame
  //--------------------------------------------------------------------

  if (stresstype == "ndxyz")
  {
    WriteNodalStress(groupname, result);
  }

  //-------------------------------------------------------------------------
  // calculation and output of element center stresses in xyz-reference frame
  //-------------------------------------------------------------------------

  else if (stresstype == "cxyz")
  {
    WriteElementCenterStress(groupname, result);
  }

  //-----------------------------------------------------------------------------------
  // calculation and output of nodal and element center stresses in xyz-reference frame
  //-----------------------------------------------------------------------------------

  else if (stresstype == "cxyz_ndxyz")
  {
    WriteNodalStress(groupname, result);

    // reset result for postprocessing and output of element center stresses
    PostResult resultelestress = PostResult(field_);
    resultelestress.next_result();
    WriteElementCenterStress(groupname,resultelestress);
  }

  else if (stresstype == "nd123")
  {
    WriteNodalEigenStress(groupname, result);
  }

  else if (stresstype == "c123")
  {
    WriteElementCenterEigenStress(groupname, result);
  }

  else if (stresstype == "c123_nd123")
  {
    WriteNodalEigenStress(groupname, result);

    // reset result for postprocessing and output of element center stresses
    PostResult resultelestress = PostResult(field_);
    resultelestress.next_result();
    WriteElementCenterEigenStress(groupname,resultelestress);
  }

  else
  {
    dserror("Unknown stress/strain type");
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalStress(const string groupname,
                                              PostResult& result)
{
  string name;
  string out;

  if (groupname=="gauss_2PK_stresses_xyz")
  {
    name="nodal_2PK_stresses_xyz";
    out="2nd Piola-Kirchhoff stresses";
  }
  else if (groupname=="gauss_cauchy_stresses_xyz")
  {
    name="nodal_cauchy_stresses_xyz";
    out="Cauchy stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name="nodal_GL_strains_xyz";
    out="Green-Lagrange strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const string filename = filename_ + "_"+ field_->name() + "."+ name;
  ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_==0)
    cout<<"writing node-based " << out << endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "node";

  WriteNodalStressStep(file,result,resultfilepos,groupname,name,6);
  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    stepsize = ((int) file.tellp())-startfilepos;
    if (stepsize <= 0) dserror("found invalid step size for result file");
  }
  else
    stepsize = 1; //use dummy value on other procs

  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
    }

    WriteNodalStressStep(file,result,resultfilepos,groupname,name,6);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = 6;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalStressStep(ofstream& file,
                                                  PostResult& result,
                                                  map<string, vector<ofstream::pos_type> >& resultfilepos,
                                                  const string groupname,
                                                  const string name,
                                                  const int numdf) const
{
  //--------------------------------------------------------------------
  // calculate nodal stresses from gauss point stresses
  //--------------------------------------------------------------------

  const RefCountPtr<std::map<int,RefCountPtr<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();

  ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","ndxyz");
  p.set("gpstressmap", data);
  RefCountPtr<Epetra_Vector> normal_stresses = LINALG::CreateVector(*(dis->DofRowMap()),true);
  RefCountPtr<Epetra_Vector> shear_stresses = LINALG::CreateVector(*(dis->DofRowMap()),true);
  dis->Evaluate(p,null,null,normal_stresses,shear_stresses,null);

  const Epetra_Map* nodemap = dis->NodeRowMap();
  RefCountPtr<Epetra_MultiVector> nodal_stresses = rcp(new Epetra_MultiVector(*nodemap, 6));

  const int numnodes = dis->NumMyRowNodes();

  for (int i=0;i<numnodes;++i)
  {
    (*((*nodal_stresses)(0)))[i] = (*normal_stresses)[3*i];
    (*((*nodal_stresses)(1)))[i] = (*normal_stresses)[3*i+1];
    (*((*nodal_stresses)(2)))[i] = (*normal_stresses)[3*i+2];
    (*((*nodal_stresses)(3)))[i] = (*shear_stresses)[3*i];
    (*((*nodal_stresses)(4)))[i] = (*shear_stresses)[3*i+1];
    (*((*nodal_stresses)(5)))[i] = (*shear_stresses)[3*i+2];
  }

  const Epetra_BlockMap& datamap = nodal_stresses->Map();

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_MultiVector> data_proc0 = rcp(new Epetra_MultiVector(*proc0map_,numdf));
  Epetra_Import proc0dofimporter(*proc0map_,datamap);
  int err = data_proc0->Import(*nodal_stresses,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);


  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------

  vector<ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, "coordinates");

  //--------------------------------------------------------------------
  // write results
  //--------------------------------------------------------------------

  const int finalnumnode = proc0map_->NumGlobalElements();

  if (myrank_==0) // ensures pointer dofgids is valid
  {
    for (int idf=0; idf<numdf; ++idf)
    {
      for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0map_
      {
        Write(file, static_cast<float>((*((*data_proc0)(idf)))[inode]));
      }
    }
  } // if (myrank_==0)

  Write(file, "END TIME STEP");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterStress(const string groupname,
                                                      PostResult& result)
{
  string name;
  string out;

  if (groupname=="gauss_2PK_stresses_xyz")
  {
    name="element_2PK_stresses_xyz";
    out="2nd Piola-Kirchhoff stresses";
  }
  else if (groupname=="gauss_cauchy_stresses_xyz")
  {
    name="element_cauchy_stresses_xyz";
    out="Cauchy stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name="element_GL_strains_xyz";
    out="Green-Lagrange strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const string filename = filename_ + "_"+ field_->name() + "."+ name;
  ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_==0)
    cout<<"writing element-based center " << out << endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "element";

  WriteElementCenterStressStep(file,result,resultfilepos,groupname,name,6);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    stepsize = ((int) file.tellp())-startfilepos;
    if (stepsize <= 0) dserror("found invalid step size for result file");
  }
  else
    stepsize = 1; //use dummy value on other procs

  while (result.next_result())
  {
    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
    }
    WriteElementCenterStressStep(file,result,resultfilepos,groupname,name,6);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = 6;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterStressStep(ofstream& file,
                                                          PostResult& result,
                                                          map<string, vector<ofstream::pos_type> >& resultfilepos,
                                                          const string groupname,
                                                          const string name,
                                                          const int numdf) const
{
  //--------------------------------------------------------------------
  // calculate element center stresses from gauss point stresses
  //--------------------------------------------------------------------

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const RefCountPtr<std::map<int,RefCountPtr<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);
  ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", data);
  RCP<Epetra_MultiVector> elestress = rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),6));
  p.set("elestress",elestress);
  dis->Evaluate(p,null,null,null,null,null);
  if (elestress==null)
  {
    dserror("vector containing element center stresses/strains not available");
  }

  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------

  vector<ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);

  const Epetra_BlockMap& datamap = elestress->Map();

  // do stupid conversion into Epetra map
  RefCountPtr<Epetra_Map> epetradatamap;
  epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));

  RefCountPtr<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);
  // sort proc0datamap so that we can loop it and get nodes in ascending order.
  std::vector<int> sortmap;
  sortmap.reserve(proc0datamap->NumMyElements());
  sortmap.assign(proc0datamap->MyGlobalElements(), proc0datamap->MyGlobalElements()+proc0datamap->NumMyElements());
  std::sort(sortmap.begin(), sortmap.end());
  proc0datamap = Teuchos::rcp(new Epetra_Map(-1, sortmap.size(), &sortmap[0], 0, proc0datamap->Comm()));

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_MultiVector> data_proc0 = rcp(new Epetra_MultiVector(*proc0datamap,numdf));
  Epetra_Import proc0dofimporter(*proc0datamap,datamap);
  int err = data_proc0->Import(*elestress,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  //--------------------------------------------------------------------
  // specify the element type
  //--------------------------------------------------------------------
  // loop over the different element types present
  if (myrank_==0)
  {
    if (eleGidPerDisType_.empty()==true) dserror("no element types available");
  }
  EleGidPerDisType::const_iterator iter;
  for (iter=eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const string ensighteleString = GetEnsightString(iter->first);
    const int numelepertype = (iter->second).size();
    vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file, ensighteleString);

    //------------------------------------------------------------------
    // write results
    //------------------------------------------------------------------

    if (myrank_==0) // ensures pointer dofgids is valid
    {
      for (int idf=0; idf<numdf; ++idf)
      {
        for (int iele=0; iele<numelepertype; iele++) // inode == lid of node because we use proc0map_
        {
          // extract element global id
          const int gid = actelegids[iele];
          // get the dof local id w.r.t. the final datamap
          int lid = proc0datamap->LID(gid);
          if (lid > -1)
          {
            Write(file, static_cast<float>((*((*data_proc0)(idf)))[lid]));
          }
        }
      }
    } // if (myrank_==0)
  }
  Write(file, "END TIME STEP");
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalEigenStress(const string groupname,
                                                   PostResult& result)
{
  vector<string> name(6);
  string out;

  if (groupname=="gauss_2PK_stresses_xyz")
  {
    name[0]="nodal_2PK_stresses_eigenval1";
    name[1]="nodal_2PK_stresses_eigenval2";
    name[2]="nodal_2PK_stresses_eigenval3";
    name[3]="nodal_2PK_stresses_eigenvec1";
    name[4]="nodal_2PK_stresses_eigenvec2";
    name[5]="nodal_2PK_stresses_eigenvec3";
    out="principal 2nd Piola-Kirchhoff stresses";
  }
  else if (groupname=="gauss_cauchy_stresses_xyz")
  {
    name[0]="nodal_cauchy_stresses_eigenval1";
    name[1]="nodal_cauchy_stresses_eigenval2";
    name[2]="nodal_cauchy_stresses_eigenval3";
    name[3]="nodal_cauchy_stresses_eigenvec1";
    name[4]="nodal_cauchy_stresses_eigenvec2";
    name[5]="nodal_cauchy_stresses_eigenvec3";
    out="principal Cauchy stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name[0]="nodal_GL_strains_eigenval1";
    name[1]="nodal_GL_strains_eigenval2";
    name[2]="nodal_GL_strains_eigenval3";
    name[3]="nodal_GL_strains_eigenvec1";
    name[4]="nodal_GL_strains_eigenvec2";
    name[5]="nodal_GL_strains_eigenvec3";
    out="principal Green-Lagrange strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  vector<bool> multiple_files(6);
  for (int i=0;i<6;++i)
  {
    multiple_files[i] = false;
  }

  // open file
  vector<string> filenames(6);
  for (int i=0;i<6;++i)
  {
    filenames[i] = filename_ + "_"+ field_->name() + "."+ name[i];
  }

  ofstream file0;
  ofstream file1;
  ofstream file2;
  ofstream file3;
  ofstream file4;
  ofstream file5;
  vector<int> startfilepos(6);
  for (int i=0;i<6;++i)
    startfilepos[i] = 0;
  if (myrank_==0)
  {
    file0.open(filenames[0].c_str());
    startfilepos[0] = file0.tellp(); // file position should be zero, but we stay flexible
    file1.open(filenames[1].c_str());
    startfilepos[1] = file1.tellp(); // file position should be zero, but we stay flexible
    file2.open(filenames[2].c_str());
    startfilepos[2] = file2.tellp(); // file position should be zero, but we stay flexible
    file3.open(filenames[3].c_str());
    startfilepos[3] = file3.tellp(); // file position should be zero, but we stay flexible
    file4.open(filenames[4].c_str());
    startfilepos[4] = file4.tellp(); // file position should be zero, but we stay flexible
    file5.open(filenames[5].c_str());
    startfilepos[5] = file5.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  vector<int> stepsize(6);
  for (int i=0;i<6;++i)
  {
    stepsize[i]=0;
  }

  if (myrank_==0)
    cout << "writing node-based " << out << endl;

  // store information for later case file creation
  for (int i=0;i<6;++i)
  {
    variableresulttypemap_[name[i]] = "node";
  }

  WriteNodalEigenStressStep(file0,file1,file2,file3,file4,file5,result,resultfilepos,groupname,name,6);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    stepsize[0] = ((int) file0.tellp())-startfilepos[0];
    if (stepsize[0] <= 0) dserror("found invalid step size for result file");
    stepsize[1] = ((int) file1.tellp())-startfilepos[1];
    if (stepsize[1] <= 0) dserror("found invalid step size for result file");
    stepsize[2] = ((int) file2.tellp())-startfilepos[2];
    if (stepsize[2] <= 0) dserror("found invalid step size for result file");
    stepsize[3] = ((int) file3.tellp())-startfilepos[3];
    if (stepsize[3] <= 0) dserror("found invalid step size for result file");
    stepsize[4] = ((int) file4.tellp())-startfilepos[4];
    if (stepsize[4] <= 0) dserror("found invalid step size for result file");
    stepsize[5] = ((int) file5.tellp())-startfilepos[5];
    if (stepsize[5] <= 0) dserror("found invalid step size for result file");
  }
  else
  {
    for (int i=0;i<6;++i)
    {
      stepsize[i] = 1; //use dummy value on other procs
    }
  }

  while (result.next_result())
  {
    const int indexsize0 = 80+2*sizeof(int)+(file0.tellp()/stepsize[0]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file0.tellp())+stepsize[0]+indexsize0>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[0];
      FileSwitcher(file0,mf,filesetmap_,resultfilepos,stepsize[0],name[0],filenames[0]);
    }
    const int indexsize1 = 80+2*sizeof(int)+(file1.tellp()/stepsize[1]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file1.tellp())+stepsize[1]+indexsize1>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[1];
      FileSwitcher(file1,mf,filesetmap_,resultfilepos,stepsize[1],name[1],filenames[1]);
    }
    const int indexsize2 = 80+2*sizeof(int)+(file2.tellp()/stepsize[2]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file2.tellp())+stepsize[2]+indexsize2>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[2];
      FileSwitcher(file2,mf,filesetmap_,resultfilepos,stepsize[2],name[2],filenames[2]);
    }
    const int indexsize3 = 80+2*sizeof(int)+(file3.tellp()/stepsize[3]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file3.tellp())+stepsize[3]+indexsize3>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[3];
      FileSwitcher(file3,mf,filesetmap_,resultfilepos,stepsize[3],name[3],filenames[3]);
    }
    const int indexsize4 = 80+2*sizeof(int)+(file4.tellp()/stepsize[4]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file4.tellp())+stepsize[4]+indexsize4>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[4];
      FileSwitcher(file4,mf,filesetmap_,resultfilepos,stepsize[4],name[4],filenames[4]);
    }
    const int indexsize5 = 80+2*sizeof(int)+(file5.tellp()/stepsize[5]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file5.tellp())+stepsize[5]+indexsize5>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[5];
      FileSwitcher(file5,mf,filesetmap_,resultfilepos,stepsize[5],name[5],filenames[5]);
    }

    WriteNodalEigenStressStep(file0,file1,file2,file3,file4,file5,result,resultfilepos,groupname,name,6);
  }
  // store information for later case file creation

  filesetmap_[name[0]].push_back(file0.tellp()/stepsize[0]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[0]] = 1;
  variablefilenamemap_[name[0]] = filenames[0];
  filesetmap_[name[1]].push_back(file1.tellp()/stepsize[1]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[1]] = 1;
  variablefilenamemap_[name[1]] = filenames[1];
  filesetmap_[name[2]].push_back(file2.tellp()/stepsize[2]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[2]] = 1;
  variablefilenamemap_[name[2]] = filenames[2];
  filesetmap_[name[3]].push_back(file3.tellp()/stepsize[3]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[3]] = 3;
  variablefilenamemap_[name[3]] = filenames[3];
  filesetmap_[name[4]].push_back(file4.tellp()/stepsize[4]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[4]] = 3;
  variablefilenamemap_[name[4]] = filenames[4];
  filesetmap_[name[5]].push_back(file5.tellp()/stepsize[5]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[5]] = 3;
  variablefilenamemap_[name[5]] = filenames[5];

  // store solution times vector for later case file creation
  for (int i=0;i<6;++i)
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name[i]] = restimes;
  }

  //append index table
  WriteIndexTable(file0, resultfilepos[name[0]]);
  WriteIndexTable(file1, resultfilepos[name[1]]);
  WriteIndexTable(file2, resultfilepos[name[2]]);
  WriteIndexTable(file3, resultfilepos[name[3]]);
  WriteIndexTable(file4, resultfilepos[name[4]]);
  WriteIndexTable(file5, resultfilepos[name[5]]);

  for (int i=0;i<6;++i) resultfilepos[name[i]].clear();

  if (file0.is_open())
    file0.close();
  if (file1.is_open())
    file1.close();
  if (file2.is_open())
    file2.close();
  if (file3.is_open())
    file3.close();
  if (file4.is_open())
    file4.close();
  if (file5.is_open())
    file5.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalEigenStressStep(ofstream& file0,
                                                       ofstream& file1,
                                                       ofstream& file2,
                                                       ofstream& file3,
                                                       ofstream& file4,
                                                       ofstream& file5,
                                                       PostResult& result,
                                                       map<string, vector<ofstream::pos_type> >& resultfilepos,
                                                       const string groupname,
                                                       vector<string> name,
                                                       const int numdf) const
{
  //--------------------------------------------------------------------
  // calculate nodal stresses from gauss point stresses
  //--------------------------------------------------------------------

  const RefCountPtr<std::map<int,RefCountPtr<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();

  ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","ndxyz");
  p.set("gpstressmap", data);
  RefCountPtr<Epetra_Vector> normal_stresses = LINALG::CreateVector(*(dis->DofRowMap()),true);
  RefCountPtr<Epetra_Vector> shear_stresses = LINALG::CreateVector(*(dis->DofRowMap()),true);
  dis->Evaluate(p,null,null,normal_stresses,shear_stresses,null);

  const Epetra_BlockMap& datamap = normal_stresses->Map();

  // do stupid conversion into Epetra map
  RefCountPtr<Epetra_Map> epetradatamap;
  epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));

  RefCountPtr<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_Vector> normal_data_proc0 = rcp(new Epetra_Vector(*proc0datamap));
  RefCountPtr<Epetra_Vector> shear_data_proc0 = rcp(new Epetra_Vector(*proc0datamap));
  Epetra_Import proc0dofimporter(*proc0datamap,*epetradatamap);
  int err1 = normal_data_proc0->Import(*normal_stresses,proc0dofimporter,Insert);
  if (err1>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err1);
  int err2 = shear_data_proc0->Import(*shear_stresses,proc0dofimporter,Insert);
  if (err2>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err2);

  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------

  vector<ofstream::pos_type>& filepos0 = resultfilepos[name[0]];
  Write(file0, "BEGIN TIME STEP");
  filepos0.push_back(file0.tellp());
  Write(file0, "description");
  Write(file0, "part");
  Write(file0, field_->field_pos()+1);
  Write(file0, "coordinates");
  vector<ofstream::pos_type>& filepos1 = resultfilepos[name[1]];
  Write(file1, "BEGIN TIME STEP");
  filepos1.push_back(file1.tellp());
  Write(file1, "description");
  Write(file1, "part");
  Write(file1, field_->field_pos()+1);
  Write(file1, "coordinates");
  vector<ofstream::pos_type>& filepos2 = resultfilepos[name[2]];
  Write(file2, "BEGIN TIME STEP");
  filepos2.push_back(file2.tellp());
  Write(file2, "description");
  Write(file2, "part");
  Write(file2, field_->field_pos()+1);
  Write(file2, "coordinates");
  vector<ofstream::pos_type>& filepos3 = resultfilepos[name[3]];
  Write(file3, "BEGIN TIME STEP");
  filepos3.push_back(file3.tellp());
  Write(file3, "description");
  Write(file3, "part");
  Write(file3, field_->field_pos()+1);
  Write(file3, "coordinates");
  vector<ofstream::pos_type>& filepos4 = resultfilepos[name[4]];
  Write(file4, "BEGIN TIME STEP");
  filepos4.push_back(file4.tellp());
  Write(file4, "description");
  Write(file4, "part");
  Write(file4, field_->field_pos()+1);
  Write(file4, "coordinates");
  vector<ofstream::pos_type>& filepos5 = resultfilepos[name[5]];
  Write(file5, "BEGIN TIME STEP");
  filepos5.push_back(file5.tellp());
  Write(file5, "description");
  Write(file5, "part");
  Write(file5, field_->field_pos()+1);
  Write(file5, "coordinates");


  //--------------------------------------------------------------------
  // write results
  //--------------------------------------------------------------------

  const int finalnumnode = proc0map_->NumGlobalElements();

  if (myrank_==0) // ensures pointer dofgids is valid
  {
    vector<Epetra_SerialDenseMatrix> eigenvec(finalnumnode, Epetra_SerialDenseMatrix(3,3));
    vector<Epetra_SerialDenseVector> eigenval(finalnumnode, Epetra_SerialDenseVector(3));

    for (int i=0;i<finalnumnode;++i)
    {
      (eigenvec[i])(0,0) = (*normal_data_proc0)[3*i];
      (eigenvec[i])(0,1) = (*shear_data_proc0)[3*i];
      (eigenvec[i])(0,2) = (*shear_data_proc0)[3*i+2];
      (eigenvec[i])(1,0) = (eigenvec[i])(0,1);
      (eigenvec[i])(1,1) = (*normal_data_proc0)[3*i+1];
      (eigenvec[i])(1,2) = (*shear_data_proc0)[3*i+1];
      (eigenvec[i])(2,0) = (eigenvec[i])(0,2);
      (eigenvec[i])(2,1) = (eigenvec[i])(1,2);
      (eigenvec[i])(2,2) = (*normal_data_proc0)[3*i+2];

      LINALG::SymmetricEigen((eigenvec[i]), eigenval[i], 3, 'V');   // option 'V' enables calculation of eigenvectors
    }

    for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0datamap
    {
      Write(file0, static_cast<float>((eigenval[inode])[0]));
      Write(file1, static_cast<float>((eigenval[inode])[1]));
      Write(file2, static_cast<float>((eigenval[inode])[2]));
    }

    for (int idf=0; idf<3; ++idf)
    {
      for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0datamap
      {
        Write(file3, static_cast<float>((eigenvec[inode])(idf,0)));
        Write(file4, static_cast<float>((eigenvec[inode])(idf,1)));
        Write(file5, static_cast<float>((eigenvec[inode])(idf,2)));
      }
    }
  } // if (myrank_==0)

  Write(file0, "END TIME STEP");
  Write(file1, "END TIME STEP");
  Write(file2, "END TIME STEP");
  Write(file3, "END TIME STEP");
  Write(file4, "END TIME STEP");
  Write(file5, "END TIME STEP");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterEigenStress(const string groupname,
                                                           PostResult& result)
{
  vector<string> name(6);
  string out;

  if (groupname=="gauss_2PK_stresses_xyz")
  {
    name[0]="element_2PK_stresses_eigenval1";
    name[1]="element_2PK_stresses_eigenval2";
    name[2]="element_2PK_stresses_eigenval3";
    name[3]="element_2PK_stresses_eigenvec1";
    name[4]="element_2PK_stresses_eigenvec2";
    name[5]="element_2PK_stresses_eigenvec3";
    out="principal 2nd Piola-Kirchhoff stresses";
  }
  else if (groupname=="gauss_cauchy_stresses_xyz")
  {
    name[0]="element_cauchy_stresses_eigenval1";
    name[1]="element_cauchy_stresses_eigenval2";
    name[2]="element_cauchy_stresses_eigenval3";
    name[3]="element_cauchy_stresses_eigenvec1";
    name[4]="element_cauchy_stresses_eigenvec2";
    name[5]="element_cauchy_stresses_eigenvec3";
    out="principal Cauchy stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name[0]="element_GL_strains_eigenval1";
    name[1]="element_GL_strains_eigenval2";
    name[2]="element_GL_strains_eigenval3";
    name[3]="element_GL_strains_eigenvec1";
    name[4]="element_GL_strains_eigenvec2";
    name[5]="element_GL_strains_eigenvec3";
    out="principal Green-Lagrange strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  vector<bool> multiple_files(6);
  for (int i=0;i<6;++i)
  {
    multiple_files[i] = false;
  }

  // open file
  vector<string> filenames(6);
  for (int i=0;i<6;++i)
  {
    filenames[i] = filename_ + "_"+ field_->name() + "."+ name[i];
  }

  ofstream file0;
  ofstream file1;
  ofstream file2;
  ofstream file3;
  ofstream file4;
  ofstream file5;
  vector<int> startfilepos(6);
  for (int i=0;i<6;++i)
    startfilepos[i] = 0;
  if (myrank_==0)
  {
    file0.open(filenames[0].c_str());
    startfilepos[0] = file0.tellp(); // file position should be zero, but we stay flexible
    file1.open(filenames[1].c_str());
    startfilepos[1] = file1.tellp(); // file position should be zero, but we stay flexible
    file2.open(filenames[2].c_str());
    startfilepos[2] = file2.tellp(); // file position should be zero, but we stay flexible
    file3.open(filenames[3].c_str());
    startfilepos[3] = file3.tellp(); // file position should be zero, but we stay flexible
    file4.open(filenames[4].c_str());
    startfilepos[4] = file4.tellp(); // file position should be zero, but we stay flexible
    file5.open(filenames[5].c_str());
    startfilepos[5] = file5.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  vector<int> stepsize(6);
  for (int i=0;i<6;++i)
  {
    stepsize[i]=0;
  }

  if (myrank_==0)
    cout << "writing element-based center " << out << endl;

  // store information for later case file creation
  for (int i=0;i<6;++i)
  {
    variableresulttypemap_[name[i]] = "element";
  }

  WriteElementCenterEigenStressStep(file0,file1,file2,file3,file4,file5,result,resultfilepos,groupname,name,6);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    stepsize[0] = ((int) file0.tellp())-startfilepos[0];
    if (stepsize[0] <= 0) dserror("found invalid step size for result file");
    stepsize[1] = ((int) file1.tellp())-startfilepos[1];
    if (stepsize[1] <= 0) dserror("found invalid step size for result file");
    stepsize[2] = ((int) file2.tellp())-startfilepos[2];
    if (stepsize[2] <= 0) dserror("found invalid step size for result file");
    stepsize[3] = ((int) file3.tellp())-startfilepos[3];
    if (stepsize[3] <= 0) dserror("found invalid step size for result file");
    stepsize[4] = ((int) file4.tellp())-startfilepos[4];
    if (stepsize[4] <= 0) dserror("found invalid step size for result file");
    stepsize[5] = ((int) file5.tellp())-startfilepos[5];
    if (stepsize[5] <= 0) dserror("found invalid step size for result file");
  }
  else
  {
    for (int i=0;i<6;++i)
    {
      stepsize[i] = 1; //use dummy value on other procs
    }
  }

  while (result.next_result())
  {
    const int indexsize0 = 80+2*sizeof(int)+(file0.tellp()/stepsize[0]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file0.tellp())+stepsize[0]+indexsize0>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[0];
      FileSwitcher(file0,mf,filesetmap_,resultfilepos,stepsize[0],name[0],filenames[0]);
    }
    const int indexsize1 = 80+2*sizeof(int)+(file1.tellp()/stepsize[1]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file1.tellp())+stepsize[1]+indexsize1>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[1];
      FileSwitcher(file1,mf,filesetmap_,resultfilepos,stepsize[1],name[1],filenames[1]);
    }
    const int indexsize2 = 80+2*sizeof(int)+(file2.tellp()/stepsize[2]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file2.tellp())+stepsize[2]+indexsize2>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[2];
      FileSwitcher(file2,mf,filesetmap_,resultfilepos,stepsize[2],name[2],filenames[2]);
    }
    const int indexsize3 = 80+2*sizeof(int)+(file3.tellp()/stepsize[3]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file3.tellp())+stepsize[3]+indexsize3>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[3];
      FileSwitcher(file3,mf,filesetmap_,resultfilepos,stepsize[3],name[3],filenames[3]);
    }
    const int indexsize4 = 80+2*sizeof(int)+(file4.tellp()/stepsize[4]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file4.tellp())+stepsize[4]+indexsize4>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[4];
      FileSwitcher(file4,mf,filesetmap_,resultfilepos,stepsize[4],name[4],filenames[4]);
    }
    const int indexsize5 = 80+2*sizeof(int)+(file5.tellp()/stepsize[5]+2)*sizeof(long);
    if (static_cast<long unsigned int>(file5.tellp())+stepsize[5]+indexsize5>= FILE_SIZE_LIMIT_)
    {
      bool mf = multiple_files[5];
      FileSwitcher(file5,mf,filesetmap_,resultfilepos,stepsize[5],name[5],filenames[5]);
    }

    WriteElementCenterEigenStressStep(file0,file1,file2,file3,file4,file5,result,resultfilepos,groupname,name,6);
  }
  // store information for later case file creation

  filesetmap_[name[0]].push_back(file0.tellp()/stepsize[0]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[0]] = 1;
  variablefilenamemap_[name[0]] = filenames[0];
  filesetmap_[name[1]].push_back(file1.tellp()/stepsize[1]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[1]] = 1;
  variablefilenamemap_[name[1]] = filenames[1];
  filesetmap_[name[2]].push_back(file2.tellp()/stepsize[2]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[2]] = 1;
  variablefilenamemap_[name[2]] = filenames[2];
  filesetmap_[name[3]].push_back(file3.tellp()/stepsize[3]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[3]] = 3;
  variablefilenamemap_[name[3]] = filenames[3];
  filesetmap_[name[4]].push_back(file4.tellp()/stepsize[4]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[4]] = 3;
  variablefilenamemap_[name[4]] = filenames[4];
  filesetmap_[name[5]].push_back(file5.tellp()/stepsize[5]);// has to be done BEFORE writing the index table
  variablenumdfmap_[name[5]] = 3;
  variablefilenamemap_[name[5]] = filenames[5];

  // store solution times vector for later case file creation
  for (int i=0;i<6;++i)
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name[i]] = restimes;
  }

  //append index table
  WriteIndexTable(file0, resultfilepos[name[0]]);
  WriteIndexTable(file1, resultfilepos[name[1]]);
  WriteIndexTable(file2, resultfilepos[name[2]]);
  WriteIndexTable(file3, resultfilepos[name[3]]);
  WriteIndexTable(file4, resultfilepos[name[4]]);
  WriteIndexTable(file5, resultfilepos[name[5]]);

  for (int i=0;i<6;++i) resultfilepos[name[i]].clear();

  if (file0.is_open())
    file0.close();
  if (file1.is_open())
    file1.close();
  if (file2.is_open())
    file2.close();
  if (file3.is_open())
    file3.close();
  if (file4.is_open())
    file4.close();
  if (file5.is_open())
    file5.close();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterEigenStressStep(ofstream& file0,
                                                               ofstream& file1,
                                                               ofstream& file2,
                                                               ofstream& file3,
                                                               ofstream& file4,
                                                               ofstream& file5,
                                                               PostResult& result,
                                                               map<string, vector<ofstream::pos_type> >& resultfilepos,
                                                               const string groupname,
                                                               vector<string> name,
                                                               const int numdf) const
{
  //--------------------------------------------------------------------
  // calculate nodal stresses from gauss point stresses
  //--------------------------------------------------------------------

  const RefCountPtr<std::map<int,RefCountPtr<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();

  ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", data);
  RCP<Epetra_MultiVector> elestress = rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),6));
  p.set("elestress",elestress);
  dis->Evaluate(p,null,null,null,null,null);
  if (elestress==null)
  {
    dserror("vector containing element center stresses/strains not available");
  }

  const Epetra_BlockMap& datamap = elestress->Map();

  // do stupid conversion into Epetra map
  RefCountPtr<Epetra_Map> epetradatamap;
  epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));

  RefCountPtr<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);
  // sort proc0datamap so that we can loop it and get nodes in ascending order.
  std::vector<int> sortmap;
  sortmap.reserve(proc0datamap->NumMyElements());
  sortmap.assign(proc0datamap->MyGlobalElements(), proc0datamap->MyGlobalElements()+proc0datamap->NumMyElements());
  std::sort(sortmap.begin(), sortmap.end());
  proc0datamap = Teuchos::rcp(new Epetra_Map(-1, sortmap.size(), &sortmap[0], 0, proc0datamap->Comm()));

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_MultiVector> data_proc0 = rcp(new Epetra_MultiVector(*proc0datamap,numdf));
  Epetra_Import proc0dofimporter(*proc0datamap,datamap);
  int err = data_proc0->Import(*elestress,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  //--------------------------------------------------------------------
  // write some key words
  //--------------------------------------------------------------------

  vector<ofstream::pos_type>& filepos0 = resultfilepos[name[0]];
  Write(file0, "BEGIN TIME STEP");
  filepos0.push_back(file0.tellp());
  Write(file0, "description");
  Write(file0, "part");
  Write(file0, field_->field_pos()+1);
  vector<ofstream::pos_type>& filepos1 = resultfilepos[name[1]];
  Write(file1, "BEGIN TIME STEP");
  filepos1.push_back(file1.tellp());
  Write(file1, "description");
  Write(file1, "part");
  Write(file1, field_->field_pos()+1);
  vector<ofstream::pos_type>& filepos2 = resultfilepos[name[2]];
  Write(file2, "BEGIN TIME STEP");
  filepos2.push_back(file2.tellp());
  Write(file2, "description");
  Write(file2, "part");
  Write(file2, field_->field_pos()+1);
  vector<ofstream::pos_type>& filepos3 = resultfilepos[name[3]];
  Write(file3, "BEGIN TIME STEP");
  filepos3.push_back(file3.tellp());
  Write(file3, "description");
  Write(file3, "part");
  Write(file3, field_->field_pos()+1);
  vector<ofstream::pos_type>& filepos4 = resultfilepos[name[4]];
  Write(file4, "BEGIN TIME STEP");
  filepos4.push_back(file4.tellp());
  Write(file4, "description");
  Write(file4, "part");
  Write(file4, field_->field_pos()+1);
  vector<ofstream::pos_type>& filepos5 = resultfilepos[name[5]];
  Write(file5, "BEGIN TIME STEP");
  filepos5.push_back(file5.tellp());
  Write(file5, "description");
  Write(file5, "part");
  Write(file5, field_->field_pos()+1);

  //--------------------------------------------------------------------
  // specify the element type
  //--------------------------------------------------------------------
  // loop over the different element types present
  if (myrank_==0)
  {
    if (eleGidPerDisType_.empty()==true) dserror("no element types available");
  }
  EleGidPerDisType::const_iterator iter;
  for (iter=eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const string ensighteleString = GetEnsightString(iter->first);
    const int numelepertype = (iter->second).size();
    vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file0, ensighteleString);
    Write(file1, ensighteleString);
    Write(file2, ensighteleString);
    Write(file3, ensighteleString);
    Write(file4, ensighteleString);
    Write(file5, ensighteleString);

    //--------------------------------------------------------------------
    // write results
    //--------------------------------------------------------------------

    if (myrank_==0) // ensures pointer dofgids is valid
    {
      vector<Epetra_SerialDenseMatrix> eigenvec(numelepertype, Epetra_SerialDenseMatrix(3,3));
      vector<Epetra_SerialDenseVector> eigenval(numelepertype, Epetra_SerialDenseVector(3));

      for (int i=0;i<numelepertype;++i)
      {
        // extract element global id
        const int gid = actelegids[i];
        // get the dof local id w.r.t. the final datamap
        int lid = proc0datamap->LID(gid);

        (eigenvec[i])(0,0) = (*(*data_proc0)(0))[lid];
        (eigenvec[i])(0,1) = (*(*data_proc0)(3))[lid];
        (eigenvec[i])(0,2) = (*(*data_proc0)(5))[lid];
        (eigenvec[i])(1,0) = (eigenvec[i])(0,1);
        (eigenvec[i])(1,1) = (*(*data_proc0)(1))[lid];
        (eigenvec[i])(1,2) = (*(*data_proc0)(4))[lid];
        (eigenvec[i])(2,0) = (eigenvec[i])(0,2);
        (eigenvec[i])(2,1) = (eigenvec[i])(1,2);
        (eigenvec[i])(2,2) = (*(*data_proc0)(2))[lid];

        LINALG::SymmetricEigen((eigenvec[i]), eigenval[i], 3, 'V');   // option 'V' enables calculation of eigenvectors
      }

      for (int iele=0; iele<numelepertype; iele++)
      {
        Write(file0, static_cast<float>((eigenval[iele])[0]));
        Write(file1, static_cast<float>((eigenval[iele])[1]));
        Write(file2, static_cast<float>((eigenval[iele])[2]));
      }

      for (int idf=0; idf<3; ++idf)
      {
        for (int iele=0; iele<numelepertype; iele++)
        {
          Write(file3, static_cast<float>((eigenvec[iele])(idf,0)));
          Write(file4, static_cast<float>((eigenvec[iele])(idf,1)));
          Write(file5, static_cast<float>((eigenvec[iele])(idf,2)));
        }
      }
    } // if (myrank_==0)

    Write(file0, "END TIME STEP");
    Write(file1, "END TIME STEP");
    Write(file2, "END TIME STEP");
    Write(file3, "END TIME STEP");
    Write(file4, "END TIME STEP");
    Write(file5, "END TIME STEP");
  }

  return;
}

#endif
