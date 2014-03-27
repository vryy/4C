/*!
  \file post_drt_ensight_structure_stress.cpp

  \brief postprocessing of structural stresses

  <pre>
  Maintainer: Lena Wiechert
              yoshihara@lnm.mw.tum.de
              http://www.lnm.mw.tum.de/Members/wiechert
              089-289-15303
  </pre>

*/



#include "post_drt_ensight_writer.H"
#include "../post_drt_common/post_drt_common.H"
#include <string>
#include "post_drt_ensight_single_field_writers.H"
#include "../linalg/linalg_utils.H"
#include "../pss_full/pss_cpp.h"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::PostStress(const std::string groupname, const std::string stresstype)
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
void StructureEnsightWriter::WriteNodalStress(const std::string groupname,
                                              PostResult& result)
{
  std::string name;
  std::string out;

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
  else if (groupname=="gauss_2PK_coupling_stresses_xyz")
  {
    name="nodal_2PK_coupling_stresses_xyz";
    out="2nd Piola-Kirchhoff coupling stresses";
  }
  else if (groupname=="gauss_cauchy_coupling_stresses_xyz")
  {
    name="nodal_cauchy_coupling_stresses_xyz";
    out="Cauchy coupling stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name="nodal_GL_strains_xyz";
    out="Green-Lagrange strains";
  }
  else if (groupname=="gauss_EA_strains_xyz")
  {
    name="nodal_EA_strains_xyz";
    out="Euler-Almansi strains";
  }
  else if (groupname=="gauss_LOG_strains_xyz")
  {
    name="nodal_LOG_strains_xyz";
    out="Logarithmic strains";
  }
  else if (groupname=="gauss_pl_GL_strains_xyz")
  {
    name="nodal_pl_GL_strains_xyz";
    out="Plastic Green-Lagrange strains";
  }
  else if (groupname=="gauss_pl_EA_strains_xyz")
  {
    name="nodal_pl_EA_strains_xyz";
    out="Plastic Euler-Almansi strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const std::string filename = filename_ + "_"+ field_->name() + "."+ name;
  std::ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_==0)
    std::cout<<"writing node-based " << out << std::endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "node";

  WriteNodalStressStep(file,result,resultfilepos,groupname,name);
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

    WriteNodalStressStep(file,result,resultfilepos,groupname,name);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = 6;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
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
void StructureEnsightWriter::WriteNodalStressStep(std::ofstream& file,
                                                  PostResult& result,
                                                  std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                                  const std::string groupname,
                                                  const std::string name) const
{
  //--------------------------------------------------------------------
  // calculate nodal stresses from gauss point stresses
  //--------------------------------------------------------------------

  const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* noderowmap = dis->NodeRowMap();

  Teuchos::ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","ndxyz");
  p.set("gpstressmap", data);
  Epetra_MultiVector* tmp = new Epetra_MultiVector(*noderowmap,6,true);
  Teuchos::RCP<Epetra_MultiVector> nodal_stress = Teuchos::rcp(tmp);
  p.set("poststress",nodal_stress);
  dis->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (nodal_stress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }

  EnsightWriter::WriteNodalResultStep(file, nodal_stress, resultfilepos, groupname, name, 6);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterStress(const std::string groupname,
                                                      PostResult& result)
{
  std::string name;
  std::string out;

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
  else if (groupname=="gauss_2PK_coupling_stresses_xyz")
  {
    name="element_2PK_coupling_stresses_xyz";
    out="2nd Piola-Kirchhoff coupling stresses";
  }
  else if (groupname=="gauss_cauchy_coupling_stresses_xyz")
  {
    name="element_cauchy_coupling_stresses_xyz";
    out="Cauchy coupling stresses";
  }
  else if (groupname=="gauss_GL_strains_xyz")
  {
    name="element_GL_strains_xyz";
    out="Green-Lagrange strains";
  }
  else if (groupname=="gauss_EA_strains_xyz")
  {
    name="element_EA_strains_xyz";
    out="Euler-Almansi strains";
  }
  else if (groupname=="gauss_LOG_strains_xyz")
  {
    name="element_LOG_strains_xyz";
    out="Logarithmic strains";
  }
  else if (groupname=="gauss_pl_GL_strains_xyz")
  {
    name="element_pl_GL_strains_xyz";
    out="Plastic Green-Lagrange strains";
  }
  else if (groupname=="gauss_pl_EA_strains_xyz")
  {
    name="element_pl_EA_strains_xyz";
    out="Plastic Euler-Almansi strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  bool multiple_files = false;

  // open file
  const std::string filename = filename_ + "_"+ field_->name() + "."+ name;
  std::ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  std::map<std::string, std::vector<std::fstream::pos_type> > resultfilepos;
  int stepsize = 0;

  if (myrank_==0)
    std::cout<<"writing element-based center " << out << std::endl;

  // store information for later case file creation
  variableresulttypemap_[name] = "element";

  WriteElementCenterStressStep(file,result,resultfilepos,groupname,name);

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
    WriteElementCenterStressStep(file,result,resultfilepos,groupname,name);
  }
  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = 6;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
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
void StructureEnsightWriter::WriteElementCenterStressStep(std::ofstream& file,
                                                          PostResult& result,
                                                          std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                                          const std::string groupname,
                                                          const std::string name) const
{
  //--------------------------------------------------------------------
  // calculate element center stresses from gauss point stresses
  //--------------------------------------------------------------------

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);
  Teuchos::ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", data);
  Teuchos::RCP<Epetra_MultiVector> elestress = Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),6));
  p.set("poststress",elestress);
  dis->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (elestress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }
  EnsightWriter::WriteElementResultStep(file, elestress, resultfilepos, groupname, name, 6, 0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalEigenStress(const std::string groupname,
                                                   PostResult& result)
{
	int numfiles = 6;
	std::vector<std::string> name(numfiles);
	std::string out;

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
  else if (groupname=="gauss_2PK_coupling_stresses_xyz")
  {
    name[0]="nodal_2PK_coupling_stresses_eigenval1";
    name[1]="nodal_2PK_coupling_stresses_eigenval2";
    name[2]="nodal_2PK_coupling_stresses_eigenval3";
    name[3]="nodal_2PK_coupling_stresses_eigenvec1";
    name[4]="nodal_2PK_coupling_stresses_eigenvec2";
    name[5]="nodal_2PK_coupling_stresses_eigenvec3";
    out="principal 2nd Piola-Kirchhoff coupling stresses";
  }
  else if (groupname=="gauss_cauchy_coupling_stresses_xyz")
  {
    name[0]="nodal_cauchy_coupling_stresses_eigenval1";
    name[1]="nodal_cauchy_coupling_stresses_eigenval2";
    name[2]="nodal_cauchy_coupling_stresses_eigenval3";
    name[3]="nodal_cauchy_coupling_stresses_eigenvec1";
    name[4]="nodal_cauchy_coupling_stresses_eigenvec2";
    name[5]="nodal_cauchy_coupling_stresses_eigenvec3";
    out="principal Cauchy coupling stresses";
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
  else if (groupname=="gauss_EA_strains_xyz")
  {
    name[0]="nodal_EA_strains_eigenval1";
    name[1]="nodal_EA_strains_eigenval2";
    name[2]="nodal_EA_strains_eigenval3";
    name[3]="nodal_EA_strains_eigenvec1";
    name[4]="nodal_EA_strains_eigenvec2";
    name[5]="nodal_EA_strains_eigenvec3";
    out="principal Euler-Almansi strains";
  }
  else if (groupname=="gauss_LOG_strains_xyz")
  {
    name[0]="nodal_LOG_strains_eigenval1";
    name[1]="nodal_LOG_strains_eigenval2";
    name[2]="nodal_LOG_strains_eigenval3";
    name[3]="nodal_LOG_strains_eigenvec1";
    name[4]="nodal_LOG_strains_eigenvec2";
    name[5]="nodal_LOG_strains_eigenvec3";
    out="principal Logarithmic strains";
  }
  else if (groupname=="gauss_pl_GL_strains_xyz")
  {
    name[0]="nodal_pl_GL_strains_eigenval1";
    name[1]="nodal_pl_GL_strains_eigenval2";
    name[2]="nodal_pl_GL_strains_eigenval3";
    name[3]="nodal_pl_GL_strains_eigenvec1";
    name[4]="nodal_pl_GL_strains_eigenvec2";
    name[5]="nodal_pl_GL_strains_eigenvec3";
    out="principal plastic Green-Lagrange strains";
  }
  else if (groupname=="gauss_pl_EA_strains_xyz")
  {
    name[0]="nodal_pl_EA_strains_eigenval1";
    name[1]="nodal_pl_EA_strains_eigenval2";
    name[2]="nodal_pl_EA_strains_eigenval3";
    name[3]="nodal_pl_EA_strains_eigenvec1";
    name[4]="nodal_pl_EA_strains_eigenvec2";
    name[5]="nodal_pl_EA_strains_eigenvec3";
    out="principal plastic Euler-Almansi strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  std::vector<bool> multiple_files(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    multiple_files[i] = false;
  }

  // open file
  std::vector<std::string> filenames(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    filenames[i] = filename_ + "_"+ field_->name() + "."+ name[i];
  }

  std::vector<Teuchos::RCP<std::ofstream> > files(numfiles);
  std::vector<int> startfilepos(numfiles);
  for (int i=0;i<numfiles;++i)
    startfilepos[i] = 0;
  for (int i=0;i<numfiles;++i)
  {
    files[i] = Teuchos::rcp(new std::ofstream);

    if (myrank_==0)
    {
      files[i]->open(filenames[i].c_str());
      startfilepos[i] = files[i]->tellp(); // file position should be zero, but we stay flexible
    }
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  std::vector<int> stepsize(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    stepsize[i]=0;
  }

  if (myrank_==0)
    std::cout << "writing node-based " << out << std::endl;

  // store information for later case file creation
  for (int i=0;i<numfiles;++i)
  {
    variableresulttypemap_[name[i]] = "node";
  }

  WriteNodalEigenStressStep(files,result,resultfilepos,groupname,name);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    for (int i=0;i<numfiles;++i)
    {
      stepsize[i] = ((int) files[i]->tellp())-startfilepos[i];
      if (stepsize[i] <= 0) dserror("found invalid step size for result file");
    }
  }
  else
  {
    for (int i=0;i<numfiles;++i)
    {
      stepsize[i] = 1; //use dummy value on other procs
    }
  }

  while (result.next_result())
  {
    for (int i=0;i<numfiles;++i)
    {
      const int indexsize = 80+2*sizeof(int)+(files[i]->tellp()/stepsize[i]+2)*sizeof(long);
      if (static_cast<long unsigned int>(files[i]->tellp())+stepsize[i]+indexsize>= FILE_SIZE_LIMIT_)
      {
        bool mf = multiple_files[i];
        FileSwitcher(*(files[i]),mf,filesetmap_,resultfilepos,stepsize[i],name[i],filenames[i]);
      }
    }

    WriteNodalEigenStressStep(files,result,resultfilepos,groupname,name);
  }
  // store information for later case file creation

  if (numfiles==6)
  {
    for (int i=0;i<numfiles;++i)
    {
      filesetmap_[name[i]].push_back(files[i]->tellp()/stepsize[i]);// has to be done BEFORE writing the index table
      if (i<3) variablenumdfmap_[name[i]] = 1;
      else     variablenumdfmap_[name[i]] = 3;
      variablefilenamemap_[name[i]] = filenames[i];
    }
  }
  else
  {
    for (int i=0;i<numfiles;++i)
    {
      filesetmap_[name[i]].push_back(files[i]->tellp()/stepsize[i]);// has to be done BEFORE writing the index table
      if (i<2) variablenumdfmap_[name[i]] = 1;
      else     variablenumdfmap_[name[i]] = 3;
      variablefilenamemap_[name[i]] = filenames[i];
    }
  }

  // store solution times vector for later case file creation
  for (int i=0;i<numfiles;++i)
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name[i]] = restimes;
  }

  //append index table
  for (int i=0;i<numfiles;++i)
  {
    WriteIndexTable(*(files[i]), resultfilepos[name[i]]);
    resultfilepos[name[i]].clear();
    if (files[i]->is_open()) files[i]->close();
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteNodalEigenStressStep(std::vector<Teuchos::RCP<std::ofstream> > files,
                                                       PostResult& result,
                                                       std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                                       const std::string groupname,
                                                       std::vector<std::string> name)
{
  //--------------------------------------------------------------------
  // calculate nodal stresses from gauss point stresses
  //--------------------------------------------------------------------

  const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* noderowmap = dis->NodeRowMap();

  Teuchos::ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","ndxyz");
  p.set("gpstressmap", data);
  Teuchos::RCP<Epetra_MultiVector> nodal_stress = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,6,true));
  p.set("poststress",nodal_stress);
  dis->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (nodal_stress==Teuchos::null)
  {
    dserror("vector containing nodal stresses/strains not available");
  }

  // Epetra_MultiVector with eigenvalues (3) and eigenvectors (9 components) in each row (=node)
  std::vector<Teuchos::RCP<Epetra_MultiVector> > nodal_eigen_val_vec(6);
  for (int i=0; i<3; ++i)
    nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,1));
  for (int i=3; i<6; ++i)
    nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,3));

  const int numnodes = dis->NumMyRowNodes();
  bool threedim = true;
  if (field_->problem()->num_dim()==2) threedim = false;

  // the three-dimensional case
  if (threedim)
  {
    for (int i=0;i<numnodes;++i)
    {
      Epetra_SerialDenseMatrix eigenvec(3,3);
      Epetra_SerialDenseVector eigenval(3);

      eigenvec(0,0) = (*((*nodal_stress)(0)))[i];
      eigenvec(0,1) = (*((*nodal_stress)(3)))[i];
      eigenvec(0,2) = (*((*nodal_stress)(5)))[i];
      eigenvec(1,0) = eigenvec(0,1);
      eigenvec(1,1) = (*((*nodal_stress)(1)))[i];
      eigenvec(1,2) = (*((*nodal_stress)(4)))[i];
      eigenvec(2,0) = eigenvec(0,2);
      eigenvec(2,1) = eigenvec(1,2);
      eigenvec(2,2) = (*((*nodal_stress)(2)))[i];

      LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

      for (int d=0; d<3; ++d)
      {
        (*((*nodal_eigen_val_vec[d])(0)))[i] = eigenval(d);
        for (int e=0; e<3; ++e)
          (*((*nodal_eigen_val_vec[d+3])(e)))[i] = eigenvec(e,d);
      }
    }
  }
  // the two-dimensional case
  else
  {
    for (int i=0;i<numnodes;++i)
    {
      Epetra_SerialDenseMatrix eigenvec(2,2);
      Epetra_SerialDenseVector eigenval(2);

      eigenvec(0,0) = (*((*nodal_stress)(0)))[i];
      eigenvec(0,1) = (*((*nodal_stress)(3)))[i];
      eigenvec(1,0) = eigenvec(0,1);
      eigenvec(1,1) = (*((*nodal_stress)(1)))[i];

      LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

      (*((*nodal_eigen_val_vec[0])(0)))[i] = eigenval(0);
      (*((*nodal_eigen_val_vec[1])(0)))[i] = eigenval(1);
      (*((*nodal_eigen_val_vec[2])(0)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[3])(0)))[i] = eigenvec(0,0);
      (*((*nodal_eigen_val_vec[3])(1)))[i] = eigenvec(1,0);
      (*((*nodal_eigen_val_vec[3])(2)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[4])(0)))[i] = eigenvec(0,1);
      (*((*nodal_eigen_val_vec[4])(1)))[i] = eigenvec(1,1);
      (*((*nodal_eigen_val_vec[4])(2)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[5])(0)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[5])(1)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[5])(2)))[i] = 0.0;
    }
  }

  for (int i=0; i<3; ++i)
    EnsightWriter::WriteNodalResultStep(*files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 1);
  for (int i=3; i<6; ++i)
    EnsightWriter::WriteNodalResultStep(*files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 3);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterEigenStress(const std::string groupname,
                                                           PostResult& result)
{
  int numfiles = 6;
  std::vector<std::string> name(numfiles);
  std::string out;

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
  else if (groupname=="gauss_2PK_coupling_stresses_xyz")
  {
    name[0]="element_2PK_coupling_stresses_eigenval1";
    name[1]="element_2PK_coupling_stresses_eigenval2";
    name[2]="element_2PK_coupling_stresses_eigenval3";
    name[3]="element_2PK_coupling_stresses_eigenvec1";
    name[4]="element_2PK_coupling_stresses_eigenvec2";
    name[5]="element_2PK_coupling_stresses_eigenvec3";
    out="principal 2nd Piola-Kirchhoff coupling stresses";
  }
  else if (groupname=="gauss_cauchy_coupling_stresses_xyz")
  {
    name[0]="element_cauchy_coupling_stresses_eigenval1";
    name[1]="element_cauchy_coupling_stresses_eigenval2";
    name[2]="element_cauchy_coupling_stresses_eigenval3";
    name[3]="element_cauchy_coupling_stresses_eigenvec1";
    name[4]="element_cauchy_coupling_stresses_eigenvec2";
    name[5]="element_cauchy_coupling_stresses_eigenvec3";
    out="principal Cauchy coupling stresses";
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
  else if (groupname=="gauss_EA_strains_xyz")
  {
    name[0]="element_EA_strains_eigenval1";
    name[1]="element_EA_strains_eigenval2";
    name[2]="element_EA_strains_eigenval3";
    name[3]="element_EA_strains_eigenvec1";
    name[4]="element_EA_strains_eigenvec2";
    name[5]="element_EA_strains_eigenvec3";
    out="principal Euler-Almansi strains";
  }
  else if (groupname=="gauss_LOG_strains_xyz")
  {
    name[0]="element_LOG_strains_eigenval1";
    name[1]="element_LOG_strains_eigenval2";
    name[2]="element_LOG_strains_eigenval3";
    name[3]="element_LOG_strains_eigenvec1";
    name[4]="element_LOG_strains_eigenvec2";
    name[5]="element_LOG_strains_eigenvec3";
    out="principal Logarithmic strains";
  }
  else if (groupname=="gauss_pl_GL_strains_xyz")
  {
    name[0]="element_pl_GL_strains_eigenval1";
    name[1]="element_pl_GL_strains_eigenval2";
    name[2]="element_pl_GL_strains_eigenval3";
    name[3]="element_pl_GL_strains_eigenvec1";
    name[4]="element_pl_GL_strains_eigenvec2";
    name[5]="element_pl_GL_strains_eigenvec3";
    out="principal plastic Green-Lagrange strains";
  }
  else if (groupname=="gauss_pl_EA_strains_xyz")
  {
    name[0]="element_pl_EA_strains_eigenval1";
    name[1]="element_pl_EA_strains_eigenval2";
    name[2]="element_pl_EA_strains_eigenval3";
    name[3]="element_pl_EA_strains_eigenvec1";
    name[4]="element_pl_EA_strains_eigenvec2";
    name[5]="element_pl_EA_strains_eigenvec3";
    out="principal plastic Euler-Almansi strains";
  }
  else
  {
    dserror("trying to write something that is not a stress or a strain");
    exit(1);
  }

  // new for file continuation
  std::vector<bool> multiple_files(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    multiple_files[i] = false;
  }

  // open file
  std::vector<std::string> filenames(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    filenames[i] = filename_ + "_"+ field_->name() + "."+ name[i];
  }

  std::vector<Teuchos::RCP<std::ofstream> > files(numfiles);
  std::vector<int> startfilepos(numfiles);
  for (int i=0;i<numfiles;++i)
    startfilepos[i] = 0;

  for (int i=0;i<numfiles;++i)
  {
    files[i] = Teuchos::rcp(new std::ofstream);

    if (myrank_==0)
    {
      files[i]->open(filenames[i].c_str());
      startfilepos[i] = files[i]->tellp(); // file position should be zero, but we stay flexible
    }
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  std::vector<int> stepsize(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    stepsize[i]=0;
  }

  if (myrank_==0)
    std::cout << "writing element-based center " << out << std::endl;

  // store information for later case file creation
  for (int i=0;i<numfiles;++i)
  {
    variableresulttypemap_[name[i]] = "element";
  }

  WriteElementCenterEigenStressStep(files,result,resultfilepos,groupname,name);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    for (int i=0;i<numfiles;++i)
    {
      stepsize[i] = ((int) files[i]->tellp())-startfilepos[i];
      if (stepsize[i] <= 0) dserror("found invalid step size for result file");
    }
  }
  else
  {
    for (int i=0;i<numfiles;++i)
    {
      stepsize[i] = 1; //use dummy value on other procs
    }
  }

  while (result.next_result())
  {
    for (int i=0;i<numfiles;++i)
    {
      const int indexsize = 80+2*sizeof(int)+(files[i]->tellp()/stepsize[i]+2)*sizeof(long);
      if (static_cast<long unsigned int>(files[i]->tellp())+stepsize[i]+indexsize>= FILE_SIZE_LIMIT_)
      {
        bool mf = multiple_files[i];
        FileSwitcher(*(files[i]),mf,filesetmap_,resultfilepos,stepsize[i],name[i],filenames[i]);
      }
    }

    WriteElementCenterEigenStressStep(files,result,resultfilepos,groupname,name);
  }
  // store information for later case file creation

  if (numfiles==6)
  {
    for (int i=0;i<numfiles;++i)
    {
      filesetmap_[name[i]].push_back(files[i]->tellp()/stepsize[i]);// has to be done BEFORE writing the index table
      if (i<3) variablenumdfmap_[name[i]] = 1;
      else     variablenumdfmap_[name[i]] = 3;
      variablefilenamemap_[name[i]] = filenames[i];
    }
  }
  else
  {
    for (int i=0;i<numfiles;++i)
    {
      filesetmap_[name[i]].push_back(files[i]->tellp()/stepsize[i]);// has to be done BEFORE writing the index table
      if (i<2) variablenumdfmap_[name[i]] = 1;
      else     variablenumdfmap_[name[i]] = 3;
      variablefilenamemap_[name[i]] = filenames[i];
    }
  }

  // store solution times vector for later case file creation
  for (int i=0;i<numfiles;++i)
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name[i]] = restimes;
  }

  //append index table
  for (int i=0;i<numfiles;++i)
  {
    WriteIndexTable(*(files[i]), resultfilepos[name[i]]);
    resultfilepos[name[i]].clear();
    if (files[i]->is_open()) files[i]->close();
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteElementCenterEigenStressStep(std::vector<Teuchos::RCP<std::ofstream> > files,
                                                               PostResult& result,
                                                               std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                                               const std::string groupname,
                                                               std::vector<std::string> name)
{
  //--------------------------------------------------------------------
  // calculate nodal stresses from gauss point stresses
  //--------------------------------------------------------------------

  const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > data =
    result.read_result_serialdensematrix(groupname);

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  Teuchos::ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","cxyz");
  p.set("gpstressmap", data);
  Teuchos::RCP<Epetra_MultiVector> elestress = Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),6));
  p.set("poststress",elestress);
  dis->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (elestress==Teuchos::null)
  {
    dserror("vector containing element center stresses/strains not available");
  }

  std::vector<Teuchos::RCP<Epetra_MultiVector> > nodal_eigen_val_vec(6);
  for (int i=0; i<3; ++i)
    nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),1));
  for (int i=3; i<6; ++i)
    nodal_eigen_val_vec[i] = Teuchos::rcp(new Epetra_MultiVector(*(dis->ElementRowMap()),3));

  const int numnodes = dis->NumMyRowNodes();
  bool threedim = true;
  if (field_->problem()->num_dim()==2) threedim = false;

  // the three-dimensional case
  if (threedim)
  {
    for (int i=0;i<dis->NumMyRowElements();++i)
    {
      Epetra_SerialDenseMatrix eigenvec(3,3);
      Epetra_SerialDenseVector eigenval(3);

      eigenvec(0,0) = (*((*elestress)(0)))[i];
      eigenvec(0,1) = (*((*elestress)(3)))[i];
      eigenvec(0,2) = (*((*elestress)(5)))[i];
      eigenvec(1,0) = eigenvec(0,1);
      eigenvec(1,1) = (*((*elestress)(1)))[i];
      eigenvec(1,2) = (*((*elestress)(4)))[i];
      eigenvec(2,0) = eigenvec(0,2);
      eigenvec(2,1) = eigenvec(1,2);
      eigenvec(2,2) = (*((*elestress)(2)))[i];

      LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

      for (int d=0; d<3; ++d)
      {
        (*((*nodal_eigen_val_vec[d])(0)))[i] = eigenval(d);
        for (int e=0; e<3; ++e)
          (*((*nodal_eigen_val_vec[d+3])(e)))[i] = eigenvec(e,d);
      }
    }
  }
  // the two-dimensional case
  else
  {
    for (int i=0;i<numnodes;++i)
    {
      Epetra_SerialDenseMatrix eigenvec(2,2);
      Epetra_SerialDenseVector eigenval(2);

      eigenvec(0,0) = (*((*elestress)(0)))[i];
      eigenvec(0,1) = (*((*elestress)(3)))[i];
      eigenvec(1,0) = eigenvec(0,1);
      eigenvec(1,1) = (*((*elestress)(1)))[i];

      LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

      (*((*nodal_eigen_val_vec[0])(0)))[i] = eigenval(0);
      (*((*nodal_eigen_val_vec[1])(0)))[i] = eigenval(1);
      (*((*nodal_eigen_val_vec[2])(0)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[3])(0)))[i] = eigenvec(0,0);
      (*((*nodal_eigen_val_vec[3])(1)))[i] = eigenvec(1,0);
      (*((*nodal_eigen_val_vec[3])(2)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[4])(0)))[i] = eigenvec(0,1);
      (*((*nodal_eigen_val_vec[4])(1)))[i] = eigenvec(1,1);
      (*((*nodal_eigen_val_vec[4])(2)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[5])(0)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[5])(1)))[i] = 0.0;
      (*((*nodal_eigen_val_vec[5])(2)))[i] = 0.0;
    }
  }

  for (int i=0; i<3; ++i)
    EnsightWriter::WriteElementResultStep(*files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 1, 0);
  for (int i=3; i<6; ++i)
    EnsightWriter::WriteElementResultStep(*files[i], nodal_eigen_val_vec[i], resultfilepos, groupname, name[i], 3, 0);
}

