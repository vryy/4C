/*----------------------------------------------------------------------*/
/*!
\file post_drt_monitor.cpp

\brief monitoring filter for one data

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/foerster
            089 - 289-15262
</pre>

*/
/*----------------------------------------------------------------------*/

/*!
\addtogroup Monitoring
*//*! @{ (documentation module open)*/
#ifdef CCADISCRET

#include <string>
#include <Teuchos_CommandLineProcessor.hpp>

#include "../post_drt_common/post_drt_common.H"
#include "../drt_lib/drt_discret.H"



/*----------------------------------------------------------------------*/
/*!
\brief pure virtual class to do monitoring
*/
/*----------------------------------------------------------------------*/
class MonWriter
{
public:
  //! constructor
  MonWriter(){}
  //! destructor
  virtual ~MonWriter()
  {
  }

  //! write something
  virtual void WriteMonFile(PostProblem& problem, string& infieldtype, int node);

protected:

  virtual PostField* GetFieldPtr(PostProblem& problem) = 0;

  virtual void CheckInfieldType(string& infieldtype) = 0;

  virtual void FieldError(int node) = 0;

  virtual void WriteHeader(ofstream& outfile) = 0;

  virtual void WriteTableHead(ofstream& outfile, int dim) = 0;

  virtual void WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim) = 0;

private:
  // undesired copy constructor
  MonWriter(const MonWriter& old);
  // undesired = operator
  MonWriter& operator= (const MonWriter& old);

}; // end of class MonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MonWriter::WriteMonFile(PostProblem& problem, string& infieldtype, int node)
{
  // create my output file
  string filename = problem.outname() + ".mon";
  ofstream outfile;
  outfile.open(filename.c_str());

  // int numdis = problem.num_discr();

  // get pointer to discretisation of actual field
  PostField* field = GetFieldPtr(problem);

  CheckInfieldType(infieldtype);

  // pointer (rcp) to actual discretisation
  RCP< DRT::Discretization > mydiscrete = field->discretization();

  // test, if this node belongs to me
  bool ismynode = mydiscrete->HaveGlobalNode(node);
    if (!ismynode) // if this node does not belong to this field ( or proc, but we should be seriell)
    FieldError(node);

  // pointer to my actual node
  const DRT::Node* mynode = mydiscrete->gNode(node);

  const Epetra_Map* mydofrowmap = mydiscrete->DofRowMap();

  // space dimension of the problem
  int dim = problem.num_dim();

  // global nodal dof numbers
  std::vector<int> gdof = mydiscrete->Dof(mynode);
  std::vector<int> ldof(dim+1); // local nodal dof numbers

  for(int i=0; i < dim+1; ++i)
  {
    ldof[i] = mydofrowmap->LID(gdof[i]);
  }

  // write header
  WriteHeader(outfile);
  outfile << node << "\n";
  outfile << "# control information: nodal coordinates   ";
  outfile << "x = " << mynode->X()[0] << "    ";
  outfile << "y = " << mynode->X()[1] << "    ";
  outfile << "z = " << mynode->X()[2] << "\n";
  outfile << "#\n";

  WriteTableHead(outfile,dim);

  // get actual results of total problem
  PostResult result = PostResult(field);

  // this is a loop over all time steps that should be written
  // writing step size is considered
  while(result.next_result())
    WriteResult(outfile,result,ldof,dim);

  outfile.close();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/ 
class FieldMonWriter : public MonWriter
{
public:
  //! constructor
  FieldMonWriter(){}
  //! destructor
  virtual ~FieldMonWriter()
  {}

protected:

  virtual PostField* GetFieldPtr(PostProblem& problem);

private:
}; // end of class FieldMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FieldMonWriter::GetFieldPtr(PostProblem& problem)
{
  return problem.get_discretization(0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/ 
class FluidMonWriter : public FieldMonWriter
{
public:
  //! constructor
  FluidMonWriter(){}
  //! destructor
  virtual ~FluidMonWriter()
  {}

protected:

  virtual void CheckInfieldType(string& infieldtype);

  void FieldError(int node);

  virtual void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  virtual void WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim);

private:
}; // end of class FluidMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidMonWriter::CheckInfieldType(string& infieldtype)
{
  if (infieldtype != "fluid")
    cout << "\nPure fluid problem, field option other than fluid has been ignored!\n\n";
}

void FluidMonWriter::FieldError(int node)
{
  dserror("Node %i does not belong to fluid field!",node);
}

void FluidMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# fluid problem, writing nodal data of node ";
}

void FluidMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "# step   time     u_x      u_y      p\n";
    break;
  case 3:
   outfile << "# step   time     u_x      u_y      u_z      p\n";  
   break;
  default: 
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

void FluidMonWriter::WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim)
{
  // get actual result vector
  RCP< Epetra_Vector > resvec = result.read_result("velnp");

  // do output
  outfile << "   " << result.step() << "    " << result.time() << "  ";
  for (int i=0; i<dim+1; ++i)
  {
    outfile << (*resvec)[ldof[i]] << "   ";
  }
  outfile << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class StructMonWriter : public FieldMonWriter
{
public:
  //! constructor
  StructMonWriter(){}
  //! destructor
  virtual ~StructMonWriter()
  {}

protected:

  virtual void CheckInfieldType(string& infieldtype);

  void FieldError(int node);

  virtual void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  void WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim);

private:
}; // end of class StructMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructMonWriter::CheckInfieldType(string& infieldtype)
{
  if (infieldtype != "structure")
    cout << "\nPure structural problem, field option other than structure has been ignored!\n\n";
}

void StructMonWriter::FieldError(int node)
{;
  dserror("Node %i does not belong to structure field!", node);
}

void StructMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# structure problem, writing nodal data of node ";
}


void StructMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "# step   time     d_x      d_y      v_x      v_y      a_x      a_y\n";
    break;
  case 3:
   outfile << "# step   time     d_x      d_y      d_z      v_x      v_y      v_z      a_x      a_y      a_z\n";  
   break;
  default: 
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

void StructMonWriter::WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim)
{
  // get actual result vector for displacement
  RCP< Epetra_Vector > resvec = result.read_result("displacement");

  // do output of time step data
  outfile << "   " << result.step() << "    " << result.time() << "  ";
  // output displacement
  for (int i=0; i<dim; ++i)
  {
    outfile << (*resvec)[ldof[i]] << "   ";
  }
  // get actual result vector for velocity
  resvec = result.read_result("velocity");
  // output velocity
  for (int i=0; i<dim; ++i)
  {
    outfile << (*resvec)[ldof[i]] << "   ";
  }
  // get actual result vector for acceleration
  resvec = result.read_result("acceleration");
  // output acceleration
  for (int i=0; i<dim; ++i)
  {
    outfile << (*resvec)[ldof[i]] << "   ";
  }
  outfile << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class AleMonWriter : public FieldMonWriter
{
public:
  //! constructor
  AleMonWriter(){}
  //! destructor
  virtual ~AleMonWriter()
  {}

protected:

  virtual void CheckInfieldType(string& infieldtype);

  void FieldError(int node);

  void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  void WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim);

private:
}; // end of class AleMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleMonWriter::CheckInfieldType(string& infieldtype)
{
  if (infieldtype != "ale")
    cout << "\nPure ALE problem, field option other than ale has been ignored!\n\n";
}

void AleMonWriter::FieldError(int node)
{
  dserror("Node %i does not belong to ALE field!", node);
}

void AleMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# ALE problem, writing nodal data of node ";
}


void AleMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "# step   time     d_x      d_y\n";
    break;
  case 3:
   outfile << "# step   time     d_x      d_y      d_z\n";  
   break;
  default: 
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

void AleMonWriter::WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim)
{
  // get actual result vector for displacement
  RCP< Epetra_Vector > resvec = result.read_result("displacement");

  // do output of time step data
  outfile << "   " << result.step() << "    " << result.time() << "  ";
  // output displacement
  for (int i=0; i<dim; ++i)
  {
    outfile << (*resvec)[ldof[i]] << "   ";
  }
  outfile << "\n";
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiFluidMonWriter : public FluidMonWriter
{
public:
  //! constructor
  FsiFluidMonWriter(){}
  //! destructor
  virtual ~FsiFluidMonWriter()
  {}

protected:

  void CheckInfieldType(string& infieldtype){};

  PostField* GetFieldPtr(PostProblem& problem);

  void WriteHeader(ofstream& outfile);

  void WriteTableHead(ofstream& outfile, int dim);

  void WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim);

private:
}; // end of class FsiFluidMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FsiFluidMonWriter::GetFieldPtr(PostProblem& problem)
{
  // get pointer to discretisation of actual field
  PostField* myfield = problem.get_discretization(1);
  if (myfield->type() != fluid)
    dserror("Fieldtype of field 1 is not fluid.");
  return myfield;
}

void FsiFluidMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# FSI problem, writing nodal data of fluid node ";
}

void FsiFluidMonWriter::WriteTableHead(ofstream& outfile, int dim)
{
  switch (dim)
  {
  case 2:
    outfile << "# step   time     d_x      d_y      u_x      u_y      p\n";
    break;
  case 3:
   outfile << "# step   time     d_x      d_y      d_z     u_x      u_y      u_z      p\n";  
   break;
  default: 
    dserror("Number of dimensions in space differs from 2 and 3!");
  }
}

void FsiFluidMonWriter::WriteResult(ofstream& outfile, PostResult& result, std::vector<int>& ldof, int dim)
{
  // get actual result vector for displacement
  RCP< Epetra_Vector > resvec = result.read_result("dispnp");

  // do output of general time step data
  outfile << "   " << result.step() << "    " << result.time() << "  ";
  // do output of displacement
  for (int i=0; i<dim; ++i)
  {
    outfile << (*resvec)[ldof[i]] << "   ";
  }
  // get actual result vector for velocity
  resvec = result.read_result("velnp");
  // do output for velocity and pressure
  for (int i=0; i<dim+1; ++i)
  {
    outfile << (*resvec)[ldof[i]] << "   ";
  }
  outfile << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiStructMonWriter : public StructMonWriter
{
public:
  //! constructor
  FsiStructMonWriter(){}
  //! destructor
  virtual ~FsiStructMonWriter()
  {}

protected:

  void CheckInfieldType(string& infieldtype){};

  PostField* GetFieldPtr(PostProblem& problem);

  void WriteHeader(ofstream& outfile);

private:
}; // end of class FsiStructMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FsiStructMonWriter::GetFieldPtr(PostProblem& problem)
{
  // get pointer to discretisation of actual field
  PostField* myfield = problem.get_discretization(0);
  if (myfield->type() != structure)
    dserror("Fieldtype of field 1 is not structure.");
  return myfield;
}

void FsiStructMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# FSI problem, writing nodal data of structure node ";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
class FsiAleMonWriter : public AleMonWriter
{
public:
  //! constructor
  FsiAleMonWriter(){}
  //! destructor
  virtual ~FsiAleMonWriter()
  {}

protected:

  void CheckInfieldType(string& infieldtype){};

  PostField* GetFieldPtr(PostProblem& problem);

  void WriteHeader(ofstream& outfile);

private:
}; // end of class FsiAleMonWriter


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PostField* FsiAleMonWriter::GetFieldPtr(PostProblem& problem)
{
  // get pointer to discretisation of actual field
  PostField* myfield = problem.get_discretization(1);
  if (myfield->type() != fluid)
    dserror("Fieldtype of field 1 is not fluid.");
  return myfield;
}

void FsiAleMonWriter::WriteHeader(ofstream& outfile)
{
  outfile << "# FSI problem, writing nodal data of ALE node ";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*!
 \brief filter main routine for monitoring filter

 Write ASCII file of one nodes history

 Note: Works in seriell version only! Requires to read one instance of the discretisation!!

 \author chfoe
 \date 11/07
 */
int main(int argc, char** argv)
{
  // command line processor to deal with arguments
  Teuchos::CommandLineProcessor my_comlinproc;
  my_comlinproc.setDocString("Post DRT monitoring filter\n\nwrite nodal result data of specified node into outfile.mon");

  bool required=true;
  /* Set an additional integer command line option which is the global node Id 
     of the node you're interested in */
  int node = 0;
  my_comlinproc.setOption("node", &node, "Global node number",required);
  /* Set a std::string command line option */
  std::string infieldtype = "fluid";
  my_comlinproc.setOption("field", &infieldtype, "Field to which output node belongs (fluid, structure, ale)");

  // my post processing problem itself
  PostProblem problem(my_comlinproc,argc,argv);

    
  switch (problem.Problemtype())
  {
    case prb_fsi:
    {
      if(infieldtype == "fluid")
      {
        FsiFluidMonWriter mymonwriter;
        mymonwriter.WriteMonFile(problem,infieldtype,node);
      }
      else if(infieldtype == "structure")
      {
        FsiStructMonWriter mymonwriter;
        mymonwriter.WriteMonFile(problem,infieldtype,node);
      }
      else if(infieldtype == "ale")
      {
        dserror("There is no ALE output. Displacements of fluid nodes can be printed.");
        FsiAleMonWriter mymonwriter;
        mymonwriter.WriteMonFile(problem,infieldtype,node);
      }
      else
        dserror("handling for monitoring of this fieldtype not yet implemented");
      break;
    }
    case prb_structure:
    {
      StructMonWriter mymonwriter;
      mymonwriter.WriteMonFile(problem,infieldtype,node);
      break;
    }
    case prb_fluid:
    {
      FluidMonWriter mymonwriter;
      mymonwriter.WriteMonFile(problem,infieldtype,node);
      break;
    }
    case prb_ale:
    {
      AleMonWriter mymonwriter;
      mymonwriter.WriteMonFile(problem,infieldtype,node);
      break;
    }
    default:
        dserror("problem type %d not yet supported", problem.Problemtype());
  }

  return 0;
}
#endif
/*! @} (documentation module close)*/
