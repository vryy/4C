/*----------------------------------------------------------------------*/
/*!
\file post_drt_gid.cpp

\brief GiD filter

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
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>
#include <string>

#include "../post_drt_common/post_drt_common.H"
extern "C" {
#include "../headers/standardtypes.h"
}

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_node.H"

using namespace std;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class EnsightWriter
{
public:

  EnsightWriter(PostField* field, string filename);
  virtual ~EnsightWriter() {}

  //! write the whole thing
  void WriteFiles();

protected:

  //! look for problem dependent result entries and write them
  virtual void WriteResults(PostField* field) = 0;

  /*!
    \brief write all time steps of a result

    Write nodal results. The results are taken from a reconstructed
    Epetra_Vector. In many cases this vector will contain just one
    variable (displacements) and thus is easy to write as a whole. At
    other times, however, there is more than one result (velocity,
    pressure) and we want to write just one part of it. So we have to
    specify which part.

    \param groupname  (in): name of the result group in the control file
    \param name       (in): name of the result to be written
    \param numdf      (in): number of dofs per node to this result
    \param from       (in): start position of values in nodes

    \author u.kue
    \date 03/07
  */
  void WriteResults(string groupname, string name, int numdf, int from=0);

private:

  // sizeof(int)==4 required!
  //void Write(int i) { geofile_.write(reinterpret_cast<const char*>(&i),sizeof(int)); }
  //void Write(unsigned i) { geofile_.write(reinterpret_cast<const char*>(&i),sizeof(int)); }

  template <class T> void Write(ofstream& os,T i) { os.write(reinterpret_cast<const char*>(&i),sizeof(T)); }
  void Write(ofstream& os,string s) { WriteString(os,s); }
  void Write(ofstream& os,const char* s) { WriteString(os,s); }

  void WriteString(ofstream& stream, string str);
  void WriteCoordinates();
  void WriteCells();
  void WriteResult(ofstream& file, PostResult& result, string groupname, string name, int numdf, int from);
  void WriteIndexTable(ofstream& file, const vector<ofstream::pos_type>& filepos);

  PostField* field_;
  string filename_;
  ofstream casefile_;
  ofstream geofile_;
  vector<ofstream::pos_type> fileposition_;
  vector<double> time_;
  map<string, vector<ofstream::pos_type> > resultfilepos_;
  vector<vector<int> > filesets_;

  static const unsigned FILE_SIZE_LIMIT = 0x7fffffff; // 2GB
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EnsightWriter::EnsightWriter(PostField* field, string filename)
  : field_(field), filename_(filename)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteFiles()
{
  string geofilename = filename_ + "_" + field_->name() + ".geo";
  geofile_.open(geofilename.c_str());
  if (!geofile_)
    dserror("failed to open file: %s", geofilename.c_str());

  Write(geofile_,"C Binary");

  // loop all results
  PostResult result = PostResult(field_);
  if (result.next_result())
  {
    time_.push_back(result.time());

    Write(geofile_,"BEGIN TIME STEP");
    fileposition_.push_back(geofile_.tellp());
    Write(geofile_,field_->name() + " geometry");
    Write(geofile_,"Comment");
    Write(geofile_,"node id given");
    Write(geofile_,"element id off");

    // part + partnumber + comment
    // careful! field_->field_pos() returns the position of the ccarat
    // field, ignoring the discretizations. So if there are many
    // discretizations in one field, we have to do something different...
    Write(geofile_,"part");
    Write(geofile_,field_->field_pos()+1);
    Write(geofile_,field_->name() + " field");

    Write(geofile_,"coordinates");
    Write(geofile_,field_->num_nodes());

    // write the grid information
    WriteCoordinates();
    WriteCells();

    Write(geofile_,"END TIME STEP");
  }
  while (result.next_result())
  {
    time_.push_back(result.time());
  }

  // append index table
  WriteIndexTable(geofile_,fileposition_);

  // now do the case file

  string casefilename = filename_ + "_" + field_->name() + ".case";
  casefile_.open(casefilename.c_str());
  casefile_ << "# created using post_drt_ensight\n"
            << "FORMAT\n\n"
            << "type:\tensight gold\n\n"
            << "GEOMETRY\n\n"
            << "model:\t2\t2\t" << geofilename << "\n\n"
            << "VARIABLE\n\n";

  // whatever result we need
  WriteResults(field_);

  casefile_ << "\nTIME\n"
            << "time set:\t\t1\n"
            << "number of steps:\t" << time_.size() << "\n"
            << "time values: ";
  for (unsigned i=0; i<time_.size(); ++i)
  {
    casefile_ << time_[i] << " ";
    if (i%8==0 && i!=0)
      casefile_ << "\n";
  }
  casefile_ << "\n\ntime set:\t\t2\n"
            << "number of steps:\t1\n"
            << "time values: " << time_[0] << "\n";
  casefile_ << "\n"
            << "FILE\n"
            << "file set:\t\t1\n"
            << "number of steps:\t" << time_.size() << "\n\n"
            << "file set:\t\t2\n"
            << "number of steps:\t1\n"
    ;
  for (unsigned int i = 0; i < filesets_.size(); ++i)
  {
    casefile_ << "\nfile set:\t\t" << 3+i << "\n";
    for (unsigned int j = 0; j < filesets_[i].size(); ++j)
    {
      casefile_ << "filename index:\t" << 1+j << "\n"
                << "number of steps:\t" << filesets_[i][j] << "\n";
    }
  }

  casefile_.close();
  geofile_.close();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteCoordinates()
{
  RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* nodemap = dis->NodeRowMap();

  if (nodemap->NumMyElements()!=nodemap->NumGlobalElements())
    dserror("filter cannot be run in parallel");

  int numnp = nodemap->NumMyElements();
  // write node ids
  for (int i=0; i<numnp; ++i)
  {
    Write(geofile_,nodemap->GID(i)+1);
  }

  // ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n
  for (int j=0; j<3; ++j)
  {
    for (int i=0; i<numnp; i++)
    {
      int gid = nodemap->GID(i);
      DRT::Node* actnode = dis->gNode(gid);
      Write(geofile_,static_cast<float>(actnode->X()[j]));
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteCells()
{
  RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* nodemap = dis->NodeRowMap();
  if (nodemap->NumMyElements()!=nodemap->NumGlobalElements())
    dserror("filter cannot be run in parallel");

  const Epetra_Map* elementmap = dis->ElementRowMap();
  if (elementmap->NumMyElements()!=elementmap->NumGlobalElements())
    dserror("filter cannot be run in parallel");

  int numele = elementmap->NumMyElements();

  // ==================================================================
  // We expect all elements in a mesh to be of the same type (shape
  // and everything)
  DRT::Element* actele = dis->gElement(elementmap->GID(0));

  switch (actele->Shape())
  {
  case DRT::Element::hex8:
    Write(geofile_,"hexa8");
    Write(geofile_,numele);
    break;
  case DRT::Element::hex20:
    Write(geofile_,"hexa20");
    Write(geofile_,numele);
    break;
  case DRT::Element::hex27:
    Write(geofile_,"hexa8");
    Write(geofile_,8*numele);
    break;
  case DRT::Element::tet4:
    Write(geofile_,"tetra4");
    Write(geofile_,numele);
    break;
  case DRT::Element::quad4:
    Write(geofile_,"quad4");
    Write(geofile_,numele);
    break;
  case DRT::Element::quad8:
    Write(geofile_,"quad8");
    Write(geofile_,numele);
    break;
  case DRT::Element::quad9:
    Write(geofile_,"quad4");
    Write(geofile_,4*numele);
    break;
  case DRT::Element::tri3:
    Write(geofile_,"tria3");
    Write(geofile_,numele);
    break;
  default:
    dserror("element type : %d", actele->Shape());
    break;
  }

  // ==================================================================

  for (int i=0; i<numele; ++i)
  {
    actele = dis->gElement(elementmap->GID(i));

    switch (actele->Shape())
    {
    // standard case with direct support
    case DRT::Element::hex20:
    case DRT::Element::hex8:
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::tet4:
    case DRT::Element::tri3:
    {
      int numnp = actele->NumNode();
      DRT::Node** nodes = actele->Nodes();
      for (int j=0; j<numnp; ++j)
      {
        Write(geofile_,nodemap->LID(nodes[j]->Id())+1);
      }
      break;
    }

    // special cases

    case DRT::Element::hex27:
    {
      DRT::Node** nodes = actele->Nodes();

      /*------------sub element 1*/
      Write(geofile_,nodemap->LID(nodes[ 0]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 8]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[20]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[11]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[12]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[21]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[24]->Id())+1);

      /*------------sub element 2*/
      Write(geofile_,nodemap->LID(nodes[ 8]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 1]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 9]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[20]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[21]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[13]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[22]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);

      /*------------sub element 3*/
      Write(geofile_,nodemap->LID(nodes[20]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 9]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 2]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[10]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[22]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[14]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[23]->Id())+1);

      /*------------sub element 4*/
      Write(geofile_,nodemap->LID(nodes[11]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[20]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[10]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 3]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[24]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[23]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[15]->Id())+1);

      /*------------sub element 5*/
      Write(geofile_,nodemap->LID(nodes[12]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[21]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[24]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 4]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[16]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[25]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[19]->Id())+1);

      /*------------sub element 6*/
      Write(geofile_,nodemap->LID(nodes[21]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[13]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[22]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[16]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 5]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[17]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[25]->Id())+1);

      /*------------sub element 7*/
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[22]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[14]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[23]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[25]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[17]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 6]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[18]->Id())+1);

      /*------------sub element 8*/
      Write(geofile_,nodemap->LID(nodes[24]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[26]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[23]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[15]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[19]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[25]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[18]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[ 7]->Id())+1);
      break;
    }

    case DRT::Element::quad9:
    {
      DRT::Node** nodes = actele->Nodes();

      /*------------sub element 1*/
      Write(geofile_,nodemap->LID(nodes[0]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[4]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[8]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[7]->Id())+1);

      /*------------sub element 2*/
      Write(geofile_,nodemap->LID(nodes[4]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[1]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[5]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[8]->Id())+1);

      /*------------sub element 3*/
      Write(geofile_,nodemap->LID(nodes[8]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[5]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[2]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[6]->Id())+1);

      /*------------sub element 4*/
      Write(geofile_,nodemap->LID(nodes[7]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[8]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[6]->Id())+1);
      Write(geofile_,nodemap->LID(nodes[3]->Id())+1);
      break;
    }

    default:
      dserror("ups");
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteResults(string groupname, string name, int numdf, int from)
{
  PostResult result = PostResult(field_);
  result.next_result();
  if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
    return;

  // new for file continuation
  unsigned int fileset = 1;
  const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
  if (nodemap->NumMyElements()!=nodemap->NumGlobalElements())
      dserror("filter cannot be run in parallel");
  int numnp = nodemap->NumMyElements();
  int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);
  int indexsize = 0;

  string filename = filename_ + "_" + field_->name() + "." + name;
  ofstream file(filename.c_str());

  WriteResult(file, result, groupname, name, numdf, from);
  while (result.next_result())
  {
    indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize >=
        FILE_SIZE_LIMIT) {
      if (fileset == 1) {
        fileset = filesets_.size()+3;
        vector<int> numsteps;
        numsteps.push_back(file.tellp()/stepsize);
        filesets_.push_back(numsteps);
        // append index table
        WriteIndexTable(file,resultfilepos_[name]);
        resultfilepos_[name].clear();
        file.close();
        rename(filename.c_str(), (filename+"001").c_str());
        file.open((filename+"002").c_str());
      }
      else {
        filesets_[fileset-3].push_back(file.tellp()/stepsize);
        ostringstream newfilename;
        newfilename << filename;
        newfilename.width(3);
        newfilename.fill('0');
        newfilename << filesets_[fileset-3].size()+1;
        // append index table
        WriteIndexTable(file,resultfilepos_[name]);
        resultfilepos_[name].clear();
        file.close();
        file.open(newfilename.str().c_str());
      }
    }
    WriteResult(file, result, groupname, name, numdf, from);
  }

  // append index table
  WriteIndexTable(file,resultfilepos_[name]);
  resultfilepos_[name].clear();

  if (fileset != 1)
  {
    filesets_[fileset-3].push_back(file.tellp()/stepsize);
    filename += "***";
  }

  if (numdf==6)
  {
    casefile_ << "tensor symm per node:\t1\t" << fileset << "\t" << name << "\t"
              << filename << "\n";
  }
  else if (numdf==3 or numdf==2)
  {
    casefile_ << "vector per node:\t1\t" << fileset << "\t" << name << "\t"
              << filename << "\n";
  }
  else if (numdf==1)
  {
    casefile_ << "scalar per node:\t1\t" << fileset << "\t" << name << "\t"
              << filename << "\n";
  }
}


/*----------------------------------------------------------------------*/
//! Write nodal values. Each node has to have the same number of dofs.
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteResult(ofstream& file, PostResult& result,
                                string groupname, string name, int numdf, int from)
{
  vector<ofstream::pos_type>& filepos = resultfilepos_[name];
  Write(file,"BEGIN TIME STEP");
  filepos.push_back(file.tellp());

  Write(file,"description");
  Write(file,"part");
  Write(file,field_->field_pos()+1);
  Write(file,"coordinates");

  RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* nodemap = dis->NodeRowMap();

  if (nodemap->NumMyElements()!=nodemap->NumGlobalElements())
    dserror("filter cannot be run in parallel");

  RefCountPtr<Epetra_Vector> data = result.read_result(groupname);
  const Epetra_BlockMap& datamap = data->Map();

  int numnp = nodemap->NumMyElements();
  for (int j=0; j<numdf; ++j)
  {
    for (int i=0; i<numnp; i++)
    {
      DRT::Node* n = dis->lRowNode(i);
      int lid = datamap.LID(dis->Dof(n,from+j));
      if (lid > -1)
      {
        Write(file,static_cast<float>((*data)[lid]));
      }
      else
      {
        // Assume we have to write a value here.
        Write<float>(file,0.);
      }
    }
  }

  // 2 component vectors in a 3d problem require a row of zeros.
  // do we really need this?
  if (numdf==2)
  {
    for (int i=0; i<numnp; i++)
    {
      Write<float>(file,0.);
    }
  }

  Write(file,"END TIME STEP");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteIndexTable(ofstream& file, const vector<ofstream::pos_type>& filepos)
{
  ofstream::pos_type lastpos = file.tellp();
  unsigned steps = filepos.size();
  Write(file,steps);
  for (unsigned i=0; i<steps; ++i)
  {
    Write<long>(file,filepos[i]);
  }
  Write(file,0);
  Write<long>(file,lastpos);
  Write(file,"FILE_INDEX");
}


/*----------------------------------------------------------------------*/
/*!
  \brief write strings of exactly 80 chars
 */
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteString(ofstream& stream, string str)
{
  // we need to write 80 bytes per string
  vector<char> s(str.begin(),str.end());
  while (s.size()<80)
    s.push_back('\0');
  stream.write(&s[0],80);
}



/*----------------------------------------------------------------------*/
/*
  \brief Writer for structural problems
 */
/*----------------------------------------------------------------------*/
class StructureEnsightWriter : public EnsightWriter
{
public:
  StructureEnsightWriter(PostField* field, string filename)
    : EnsightWriter(field,filename) {}

protected:

  virtual void WriteResults(PostField* field);
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteResults(PostField* field)
{
  EnsightWriter::WriteResults("displacement","displacement",field->problem()->num_dim());
  EnsightWriter::WriteResults("velocity","velocity",field->problem()->num_dim());
  EnsightWriter::WriteResults("acceleration","acceleration",field->problem()->num_dim());
}


/*----------------------------------------------------------------------*/
/*
  \brief Writer for fluid problems
 */
/*----------------------------------------------------------------------*/
class FluidEnsightWriter : public EnsightWriter
{
public:
  FluidEnsightWriter(PostField* field, string filename)
    : EnsightWriter(field,filename) {}

protected:

  virtual void WriteResults(PostField* field);
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidEnsightWriter::WriteResults(PostField* field)
{
  EnsightWriter::WriteResults("velnp","velocity",field->problem()->num_dim());
  EnsightWriter::WriteResults("velnp","pressure",1,field->problem()->num_dim());
  EnsightWriter::WriteResults("residual","residual",field->problem()->num_dim());
  EnsightWriter::WriteResults("dispnp","displacement",field->problem()->num_dim());
}


/*----------------------------------------------------------------------*/
/*
  \brief Writer for ale problems
 */
/*----------------------------------------------------------------------*/
class AleEnsightWriter : public EnsightWriter
{
public:
  AleEnsightWriter(PostField* field, string filename)
    : EnsightWriter(field,filename) {}

protected:

  virtual void WriteResults(PostField* field);
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleEnsightWriter::WriteResults(PostField* field)
{
  EnsightWriter::WriteResults("dispnp","displacement",field->problem()->num_dim());
}


/*----------------------------------------------------------------------*/
/*!
  \brief filter main routine

  Write binary ensight format.

  The ens_checker that is part of ensight can be used to verify the
  files generated here.

  \author u.kue
  \date 03/07
 */
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  Teuchos::CommandLineProcessor My_CLP;
  My_CLP.setDocString(
    "Post DRT ensight Filter\n"
    );

  PostProblem problem = PostProblem(My_CLP,argc,argv);

#if 0
  for (int i = 0; i<problem.num_discr(); ++i)
  {
    PostField* field = problem.get_discretization(i);
    StructureEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
  }
#endif

  // each problem type is different and writes different results
  switch (problem.Problemtype())
  {
  case prb_fsi:
  {
    string basename = problem.basename();
//     PostField* structfield = problem.get_discretization(0);
//     StructureEnsightWriter structwriter(structfield, basename);
//     structwriter.WriteFiles();

    PostField* fluidfield = problem.get_discretization(1);
    FluidEnsightWriter fluidwriter(fluidfield, basename);
    fluidwriter.WriteFiles();
    break;
  }
  case prb_structure:
  {
    PostField* field = problem.get_discretization(0);
    StructureEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
    break;
  }
  case prb_fluid:
  {
    PostField* field = problem.get_discretization(0);
    FluidEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
    break;
  }
  case prb_ale:
  {
    PostField* field = problem.get_discretization(0);
    AleEnsightWriter writer(field, problem.basename());
    writer.WriteFiles();
    break;
  }
  default:
    dserror("problem type %d not yet supported", problem.Problemtype());
  }

  return 0;
}


#endif
#endif
