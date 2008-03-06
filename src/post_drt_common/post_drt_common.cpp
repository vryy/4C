/*----------------------------------------------------------------------*/
/*!
\file post_drt_common.cpp

\brief drt binary filter library

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <algorithm>
#include <stack>
#include "post_drt_common.H"

#include "../io/hdf_reader.H"
#include <sstream>

#include "../drt_lib/drt_globalproblem.H"
#include <EpetraExt_Transpose_CrsGraph.h>

#if 1
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fluid/drt_periodicbc.H"
#endif

extern "C" {

  void create_communicators();

}

/*----------------------------------------------------------------------*
 * Some functions and global variables, most of them are only
 * important for linking.
 *----------------------------------------------------------------------*/



/* There are some global variables in ccarat that are needed by the
 * service functions. We need to specify them here and set them up
 * properly. */
struct _FILES   allfiles;
struct _PAR     par;
struct _FIELD  *field;
struct _GENPROB genprob;
struct _MATERIAL *mat;
struct _IO_FLAGS ioflags;

// not actually used, but referenced in par_assignmesh.c
struct _PARTITION  *partition;
struct _SOLVAR  *solv;

char* fieldnames[] = FIELDNAMES;



/*----------------------------------------------------------------------*
 * The main part of this file. All the functions of the three classes
 * PostProblem, PostField and PostResult are defined here.
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 * the Constructor of PostProblem
 *----------------------------------------------------------------------*/
PostProblem::PostProblem(Teuchos::CommandLineProcessor& CLP,
                         int argc, char** argv)
  : start_(0),end_(-1),step_(1)
{
#ifdef PARALLEL
  MPI_Init(&argc,&argv);
#endif

#ifdef DEBUG
  dsinit();
#endif

  string file = "xxx";
  string output;

  CLP.throwExceptions(false);
  CLP.setOption("start",&start_,"first time step to read");
  CLP.setOption("end",&end_,"last time step to read");
  CLP.setOption("step",&step_,"number of time steps to jump");
  CLP.setOption("file",&file,"control file to open");
  CLP.setOption("output",&output,"output file name [defaults to control file name]");
  CLP.setOption("stresstype",&stresstype_,"stress output type [cxyz or ndxyz]");

  CommandLineProcessor::EParseCommandLineReturn
    parseReturn = CLP.parse(argc,argv);

  if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED)
  {
    exit(0);
  }
  if (parseReturn != CommandLineProcessor::PARSE_SUCCESSFUL)
  {
    exit(1);
  }

  if (file=="")
  {
    CLP.printHelpMessage(argv[0],cout);
    exit(1);
  }

  if (file.length()<=8 or file.substr(file.length()-8,8)!=".control")
  {
    file += ".control";
  }

  if (output=="")
  {
    output = file.substr(0,file.length()-8);
  }

  if (stresstype_=="")
  {
    stresstype_ = "none";
  }

  result_group_ = vector<MAP*>();
  setup_filter(file, output);

  ndim_ = map_read_int(&control_table_, "ndim");
  dsassert((ndim_ == 2) || (ndim_ == 3), "illegal dimension");

  char* problem_names[] = PROBLEMNAMES;
  char* type = map_read_string(&control_table_, "problem_type");
  int i;
  for (i=0; problem_names[i] != NULL; ++i)
  {
    if (strcmp(type, problem_names[i])==0)
    {
      problemtype_ = static_cast<PROBLEM_TYP>(i);
      break;
    }
  }
  if (problem_names[i] == NULL)
  {
    dserror("unknown problem type '%s'", type);
  }


  /*--------------------------------------------------------------------*/
  /* collect all result groups */
  SYMBOL* symbol = map_find_symbol(&control_table_, "result");
  while (symbol != NULL)
  {
    if (!symbol_is_map(symbol))
    {
      dserror("failed to get result group");
    }
    result_group_.push_back(symbol_map(symbol));
    symbol = symbol->next;
  }
  genprob.numfld = 0;

  read_meshes();
}

/*----------------------------------------------------------------------*
 * the Destructor
 *----------------------------------------------------------------------*/
PostProblem::~PostProblem()
{
  destroy_map(&control_table_);
#ifdef PARALLEL
  MPI_Finalize();
#endif
}


/*----------------------------------------------------------------------*
 * returns a pointer to the num-th discretization
 *----------------------------------------------------------------------*/
PostField* PostProblem::get_discretization(int num)
{
  for (unsigned i=0; i<fields_.size(); ++i)
  {
    if (fields_[i].field_pos()==num)
      return &fields_[i];
  }
  dserror("no field with position %d", num);
  return NULL;
}



/*----------------------------------------------------------------------*
 * returns the Epetra Communicator object
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Comm> PostProblem::comm()
{
  return comm_;
}


/*----------------------------------------------------------------------*
 * initializes all the data a filter needs. This function is called by
 * the Constructor.  (private)
 *----------------------------------------------------------------------*/
void PostProblem::setup_filter(string control_file_name, string output_name)
{
  int length;
  MAP* table;
  MAP temp_table;

#ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &par.myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &par.nprocs);
  Epetra_MpiComm* com = new Epetra_MpiComm(MPI_COMM_WORLD);
  comm_ = rcp(com);
#else
  par.myrank=0;
  par.nprocs=1;
  Epetra_SerialComm* com = new Epetra_SerialComm();
  comm_ = rcp(com);
#endif

  /* The warning system is not set up. It's rather stupid anyway. */

  /* We need to open the error output file. The other ones are not
   * important. */
  basename_ = control_file_name.substr(0,control_file_name.length()-8);
  outname_ = output_name;
  length = output_name.length();
  strcpy(allfiles.outputfile_name, outname_.c_str());
  strcpy(allfiles.outputfile_name+length-8, ".post.log");
  allfiles.out_err = fopen(allfiles.outputfile_name, "w");

  parse_control_file(&control_table_, control_file_name.c_str());

  /*
   * Now that we've read the control file given by the user we have to
   * take care of any previous (restarted) control files. These files
   * build a chain. So as long as a previous file exists we have to
   * open it and read any groups of results with smaller step numbers
   * than the results we've already read. */

  /* The general idea is to merge the different control files in
   * memory. If one step is written several times the last version is
   * used. */

  table = &control_table_;

  /* copy directory information */
  string::size_type separator = basename_.rfind('/', string::npos);
  if (separator != string::npos)
  {
    input_dir_ = basename_.substr(0,separator+1);
  }
  else
  {
    input_dir_ = "";
  }

  while (map_symbol_count(table, "restarted_run") > 0)
  {
    FILE* f;
    SYMBOL* first_result;
    SYMBOL* previous_results;
    INT first_step;
    SYMBOL dummy_symbol;
    INT counter;

    /* copy directory information */
    control_file_name = input_dir_;

    /* copy file name */
    control_file_name += map_read_string(table, "restarted_run");
    control_file_name += ".control";

    /* test open to see if it exists */
    f = fopen(control_file_name.c_str(), "rb");
    if (f == NULL)
    {
      printf("Restarted control file '%s' does not exist. Skip previous results.\n",
             control_file_name.c_str());
      break;
    }
    fclose(f);

    /* copy all the result steps that are previous to this file */
    /* We assume that the results are ordered! */

    /*------------------------------------------------------------------*/
    /* find the first result in the current table */
    first_result = map_find_symbol(&control_table_, "result");
    if (first_result == NULL)
    {
      dserror("no result sections in control file '%s'\n", control_file_name.c_str());
    }
    while (first_result->next != NULL)
    {
      first_result = first_result->next;
    }
    first_step = map_read_int(symbol_map(first_result), "step");


    /*------------------------------------------------------------------*/
    /* done with this control file */
    if (table != &control_table_)
    {
      destroy_map(table);
    }

    /* The first time we reach this place we had just used the main
     * control table. But from now on we are interessted in the
     * previous control files we read. */
    table = &temp_table;

    /* read the previous control file */
    parse_control_file(table, control_file_name.c_str());
    printf("read restarted control file: %s\n", control_file_name.c_str());

    /* find the previous results */

    counter = 0;

    /*
     * the dummy_symbol is a hack that allows us to treat all results
     * in the list the same way (use the same code). Without it we'd
     * need special conditions for the first entry. */
    previous_results = &dummy_symbol;
    previous_results->next = map_find_symbol(table, "result");
    while (previous_results->next != NULL)
    {
      SYMBOL* result;
      INT step;
      result = previous_results->next;
      step = map_read_int(symbol_map(result), "step");

      if (step < first_step)
      {
        /* found it */
        /* Now we simply switch all previous results to our main
         * map. The assumption is a perfect ordering */
        map_prepend_symbols(&control_table_, "result", result,
                            map_symbol_count(table, "result") - counter);
        previous_results->next = NULL;

        /*
         * In case all results go to the main map we have to disconnect
         * them explicitly. */
        if (previous_results == &dummy_symbol)
        {
          map_disconnect_symbols(table, "result");
        }
        break;
      }

      /* Not found yet. Go up one result. */
      previous_results = previous_results->next;
      counter += 1;
    }
  }

}

/*----------------------------------------------------------------------*
 * reads the mesh files and calls 'getfield()' for each 'field'-entry
 * in the mesh file (currently it reads only the fields with step ==
 * 0). This function is called by the Constructor. (private)
 *----------------------------------------------------------------------*/
void PostProblem::read_meshes()
{
  SYMBOL* mesh = map_find_symbol(&control_table_,"field");
  if (mesh == NULL)
    dserror("No field found.");

  // We have to reverse the traversal of meshes we get from the control file
  // in order to get the same dof numbers in all discretizations as we had
  // during our calculation.
  // The order inside the control file is important!
  // Discretizations have to be FillComplete()ed in the same order as during
  // the calculation!
  std::stack<MAP*> meshstack;

  while (mesh != NULL)
  {
    // only those fields with a mesh file entry are readable here
    // (each control file is bound to include at least one of those)
    if (map_find_symbol(symbol_map(mesh), "mesh_file")!=NULL)
    {
      meshstack.push(symbol_map(mesh));
    }
    mesh = mesh->next;
  }
  mesh = NULL;

  while (not meshstack.empty())
  {
    MAP* meshmap = meshstack.top();
    meshstack.pop();

    int field_pos = map_read_int(meshmap, "field_pos");
    int disnum = map_read_int(meshmap, "discretization");

    bool havefield = false;
    for (unsigned i=0; i<fields_.size(); ++i)
    {
      if (fields_[i].field_pos()==field_pos and
          fields_[i].disnum()==disnum)
      {
        havefield = true;
        break;
      }
    }

    // only read a field that has not yet been read
    // for now we do not care at which step this field was defined
    // if we want to support changing meshes one day, we'll need to
    // change this code...
    if (not havefield)
    {
      int step;
      if (!map_find_int(meshmap,"step",&step))
        dserror("No step information in field.");

      PostField currfield = getfield(meshmap);

      int num_output_procs;
      if (!map_find_int(meshmap,"num_output_proc",&num_output_procs))
      {
        num_output_procs = 1;
      }
      currfield.set_num_output_procs(num_output_procs);
      char* fn;
      if (!map_find_string(meshmap,"mesh_file",&fn))
        dserror("No meshfile name for discretization %s.", currfield.discretization()->Name().c_str());
      string filename = fn;
      IO::HDFReader reader = IO::HDFReader(input_dir_);
      reader.Open(filename,num_output_procs);

      RefCountPtr<vector<char> > node_data =
        reader.ReadNodeData(step, comm_->NumProc(), comm_->MyPID());
      currfield.discretization()->UnPackMyNodes(node_data);

      RefCountPtr<vector<char> > element_data =
        reader.ReadElementData(step, comm_->NumProc(), comm_->MyPID());
      currfield.discretization()->UnPackMyElements(element_data);

      // read periodic boundary conditions if available
      RefCountPtr<vector<char> > cond_pbcsline =
        reader.ReadCondition(step, comm_->NumProc(), comm_->MyPID(), "LinePeriodic");
      currfield.discretization()->UnPackCondition(cond_pbcsline, "LinePeriodic");
      RefCountPtr<vector<char> > cond_pbcssurf =
        reader.ReadCondition(step, comm_->NumProc(), comm_->MyPID(), "SurfacePeriodic");
      currfield.discretization()->UnPackCondition(cond_pbcssurf, "SurfacePeriodic");

      // read XFEMCoupling boundary conditions if available
      RefCountPtr<vector<char> > cond_xfem =
        reader.ReadCondition(step, comm_->NumProc(), comm_->MyPID(), "XFEMCoupling");
      currfield.discretization()->UnPackCondition(cond_xfem, "XFEMCoupling");

      // before we can call discretization functions (in parallel) we
      // need field based communicators.
      create_communicators();

      // setup of parallel layout: create ghosting of already distributed nodes+elems
#ifdef PARALLEL
      setup_ghosting(currfield.discretization());
#endif

      //distribute_drt_grids();
      currfield.discretization()->FillComplete();

      // -------------------------------------------------------------------
      // connect degrees of freedom for periodic boundary conditions
      // -------------------------------------------------------------------
      if(!cond_pbcssurf->empty())
      {
        PeriodicBoundaryConditions::PeriodicBoundaryConditions pbc(currfield.discretization());
        pbc.UpdateDofsForPeriodicBoundaryConditions();
      }

      fields_.push_back(currfield);
    }
  }
}

/*----------------------------------------------------------------------*
 * creates and returns a PostField insance from a field MAP. keeps
 * track of global field variable (private)
 *----------------------------------------------------------------------*/
PostField PostProblem::getfield(MAP* field_info)
{
  int field_pos = map_read_int(field_info, "field_pos");
  int disnum = map_read_int(field_info, "discretization");
  char* field_name = map_read_string(field_info, "field");
  char* dis_name = map_read_string(field_info, "dis_name");
  int numnd = map_read_int(field_info, "num_nd");
  int numele = map_read_int(field_info, "num_ele");
  int type;

  char* fieldnames[] = FIELDNAMES;
  int i;
  for (i=0; fieldnames[i]!=NULL; ++i)
  {
    if (strcmp(fieldnames[i], field_name)==0)
    {
      type = i;
      break;
    }
  }
  if (fieldnames[i]==NULL)
  {
    type = i;
    dserror("unknown field type '%d'", i);
  }

  if (field_pos >= genprob.numfld)
  {
    if (field)
    {
      field = (FIELD*)CCAREALLOC(field,(field_pos+1)*sizeof(FIELD));
    }
    else
    {
      field = (FIELD*)CCAMALLOC((field_pos+1)*sizeof(FIELD));
    }

    for (i=genprob.numfld; i<field_pos+1; ++i)
    {
      field[i].dis = NULL;
      field[i].ndis = 0;
    }
    genprob.numfld = field_pos+1;
  }

  field[field_pos].fieldtyp = (FIELDTYP)type;

  if (field[field_pos].ndis <= disnum)
  {
    field[field_pos].ndis = disnum+1;
  }

  RefCountPtr<DRT::Discretization> dis = rcp(new DRT::Discretization(dis_name,comm_));
  DRT::Problem::Instance()->SetDis(field_pos, disnum, dis);
  return PostField(dis, this, field_pos, field_name,
                   (FIELDTYP)type, disnum,numnd,numele);
}


/*-----------------------------------------------------------------------*
 * set up the parallel discretization layout(ghosting!) (private) gjb 11/07
 *-----------------------------------------------------------------------*/
void PostProblem::setup_ghosting(RefCountPtr<DRT::Discretization> dis)
{
    // this section is strongly oriented on what is done during the
	// usual BACI setup phase.
    // reference: src/drt_lib/drt_inputreader.cpp
	// ToDo: make PostProblem::setup_ghosting a method of the dicretization class itself

	int numnode=dis->NumMyColNodes();
    vector<int> nids(numnode);         // vector for global node ids

    // we have to know about the global node ids on each processor.
    // since the current discretization is NOT(!!!) Filled() yet, the following was
    // the only solution to create the ghosting for ALREADY DISTRIBUTED nodes/elements.
    int nodecount=0;
    int ngid=0;
    while(nodecount < numnode)
    {
    	if (dis->HaveGlobalNode(ngid)) // do we have this global node id on this processor?
    	{ nids[nodecount]= ngid;
    	  nodecount+=1;
    	}
  	    ngid+=1;
    }

    // now create a preliminary node row map
    RefCountPtr<Epetra_Map> rownodes = rcp(new Epetra_Map(-1,nids.size(),&nids[0],0,*comm_));
    nids.clear();

    // construct graphs
    RefCountPtr<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*rownodes,81,false));
    RefCountPtr<Epetra_CrsGraph> finalgraph = rcp(new Epetra_CrsGraph(Copy,*rownodes,81,false));

#if 0
    //cout<<dis->NumMyColNodes()<<endl;
    //rownodes->Print(cout);
    dis->Print(cout);
#endif

    // loop over all elements located on this processor (no ghosting existent)
    list<vector<int> > elementnodes;

    int numele = dis->NumMyColElements();
    int elecount=0;
    int elegid=0;
    while(elecount < numele)
    {
	  if (dis->HaveGlobalElement(elegid)) // do we have this global element id on this processor?
	  {
	    elecount+=1;
	    // get the node ids of this element...
	    const int  numnode = dis->gElement(elegid)->NumNode();
	    const int* nodeids1 = dis->gElement(elegid)->NodeIds();
	    // ... and store them locally
	    elementnodes.push_back(vector<int>(nodeids1, nodeids1+numnode));
	  }
	elegid+=1;
    } // while

#if 0
      cout << "\n\nelementnodes: size=" << elementnodes.size() << endl;
      for (list<vector<int> >::iterator i=elementnodes.begin();
           i!=elementnodes.end();
           ++i)
      {
        copy(i->begin(), i->end(), ostream_iterator<int>(cout, " "));
        cout << endl;
      }
#endif
#if 0
      // storage set for refused row numbers
      set<int> refusedrowgid;
#endif

    // initial fill-up of the node dependency graph
    for (list<vector<int> >::iterator i=elementnodes.begin();
           i!=elementnodes.end();
           ++i)
      {
        // get the node ids of this element
        int  numelnodes = static_cast<int>(i->size());
        int* nodeids = &(*i)[0];

        // loop nodes and add this topology to the row in the graph of every node
        for (int i=0; i<numelnodes; ++i)
        {
          // ensure first that node with nodeids[i] belongs to this processor
          if (rownodes->MyGID(nodeids[i]))
          {
        	int err = graph->InsertGlobalIndices(nodeids[i],numelnodes,nodeids);
            if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
            // do the same for finalgraph
        	err = finalgraph->InsertGlobalIndices(nodeids[i],numelnodes,nodeids);
            if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
          }
          else
          {  // this case does NOT happen in the normal BACI inputreader;
        	 // however, here we loose information, which has to be fetched from the transposed graph
#if 0
              cout<<"Proc "<<par.myrank<<" refused InsertGlobelIndices with nodeids[i]="<<nodeids[i]<<endl;
              // store refused row number for later use?
        	  refusedrowgid.insert(nodeids[i]);
#endif
          }
        }
      } // end for-loop over elementnodes

    elementnodes.clear();


      // finalize construction of initial graph
      // FillComplete() is necessary for the following transposition)
      int err = graph->FillComplete(*rownodes,*rownodes);
      if (err) dserror("graph->FillComplete returned %d",err);

      // create transformation object
      RefCountPtr<EpetraExt::CrsGraph_Transpose::CrsGraph_Transpose> graphtransposer = rcp(new EpetraExt::CrsGraph_Transpose());
      // create graph object
      RefCountPtr<Epetra_CrsGraph> tgraph = rcp(new Epetra_CrsGraph(Copy,*rownodes,81,false));
      Epetra_CrsGraph& new_graph = ((*tgraph));
      // finally do the transposition
      new_graph=(*graphtransposer)(*graph);
      // free memory of the graph transposer object
      graphtransposer=null;

#if 0
      if (!(refusedrowgid.empty()))
      {
    	  set<int>::iterator j;
    	  for (j = refusedrowgid.begin(); j != refusedrowgid.end(); ++j)
    	  	{
             cout<<"Proc: "<<par.myrank<<"  Iterator = "<<*j<<endl;
    	  	}
      }
#endif

    // loop over my rows of transposed graph and insert dependencies to final graph
    for (int lid = 0; lid < tgraph->NumMyRows(); ++lid)
    {
     int gid = graph->GRID(lid); // get global id
     int maxrowlength=tgraph->NumGlobalIndices(gid);
     int numelnodes=0;
     vector<int> nodeentries(maxrowlength,0);
     int* nodeids = &(nodeentries)[0];
     // what is the Global row index??
     //cout<<"Proc: "<<par.myrank<<"  Iterator = "<<j <<"  gid: "<<gid<<endl;
     if (gid > (-1))
          {
        	 // int err = tgraph->ExtractGlobalRowView(gid,numelnodes,nodeids);
        	  int err = tgraph->ExtractGlobalRowCopy(gid,maxrowlength,numelnodes,nodeids);
              if (err<0) dserror("error while extracting global row copy from transposed graph: %d",err);
#if 0
              for (int k=0;k<numelnodes;++k)
               {
        	    cout<<"nodeids["<<k<<"] = "<<nodeids[k]<<endl;
               }
#endif
          // insert node gids into final graph
          err = finalgraph->InsertGlobalIndices(gid,numelnodes,nodeids);
          if (err<0) dserror("finalgraph->InsertGlobalIndices returned %d",err);
          }
    }

    // finalize construction of final graph
    err = finalgraph->FillComplete(*rownodes,*rownodes);
    if (err) dserror("graph->FillComplete returned %d",err);

    // no partition of the graph using metis here, since we want to keep the currently
    // existing parallel distribution of nodes and elements

    // replace rownodes, colnodes with row and column maps from the graph
    // do stupid conversion from Epetra_BlockMap to Epetra_Map
    const Epetra_BlockMap& brow = finalgraph->RowMap();
    const Epetra_BlockMap& bcol = finalgraph->ColMap();
    rownodes = rcp(new Epetra_Map(brow.NumGlobalElements(),
                                   brow.NumMyElements(),
                                   brow.MyGlobalElements(),
                                   0,
                                   *comm_));
    RefCountPtr<Epetra_Map> colnodes;
    colnodes = rcp(new Epetra_Map(bcol.NumGlobalElements(),
                                   bcol.NumMyElements(),
                                   bcol.MyGlobalElements(),
                                   0,
                                   *comm_));

#if 0
    //rownodes->Print(cout);
    //colnodes->Print(cout);
    //graph->Print(cout);
    //tgraph->Print(cout);
    finalgraph->Print(cout);
#endif

    // clean up
    graph = null;
    finalgraph = null;
    tgraph = null;

    // distribute ghost nodes resolving the node dependencies given by the final graph
    dis->ExportColumnNodes(*colnodes);

    // now we do the element stuff
    RefCountPtr<Epetra_Map> elerowmap;
    RefCountPtr<Epetra_Map> elecolmap;

    // now we have all elements in a linear map roweles
    // build resonable maps for elements from the
    // already valid and final node maps
    // note that nothing is actually redistributed in here
    dis->BuildElementRowColumn(*rownodes,*colnodes,elerowmap,elecolmap);

    // we can now export elements to resonable row element distribution
    dis->ExportRowElements(*elerowmap);

    // export to the column map / create ghosting of elements
    dis->ExportColumnElements(*elecolmap);

#if 0
    dis->Print(cout);
#endif

    return;
}


/*----------------------------------------------------------------------*
 * The Constructor of PostField
 *----------------------------------------------------------------------*/
PostField::PostField(string name, RefCountPtr<Epetra_Comm> comm,
                     PostProblem* problem, int field_pos, string field_name,
                     FIELDTYP type, int disnum, const int numnd, const int numele):
  dis_(rcp(new DRT::Discretization(name,comm))),
  problem_(problem),
  field_pos_(field_pos),
  field_name_(field_name),
  type_(type),
  disnum_(disnum),
  numnd_(numnd),
  numele_(numele)
{
}

/*----------------------------------------------------------------------*
 * Another Constructor of PostField. This one takes the discretization
 * directly
 *----------------------------------------------------------------------*/
PostField::PostField(RefCountPtr<DRT::Discretization> dis, PostProblem* problem,
                     int field_pos, string field_name, FIELDTYP type,
                     int disnum, const int numnd, const int numele)
  : dis_(dis),
    problem_(problem),
    field_pos_(field_pos),
    field_name_(field_name),
    type_(type),
    disnum_(disnum),
    numnd_(numnd),
    numele_(numele)
{
}

/*----------------------------------------------------------------------*
 * The Destructor
 *----------------------------------------------------------------------*/
PostField::~PostField()
{
}





/*----------------------------------------------------------------------*
 * The Constructor of PostResult
 *----------------------------------------------------------------------*/
PostResult::PostResult(PostField* field):
  field_(field),
  pos_(-1),
  group_(NULL),
  file_((field->problem()->input_dir()))
{
}


/*----------------------------------------------------------------------*
 * The Destructor of PostResult
 *----------------------------------------------------------------------*/
PostResult::~PostResult()
{
  close_result_files();
}

/*----------------------------------------------------------------------*
 * get timesteps when the solution is written
 *----------------------------------------------------------------------*/
vector<double> PostResult::get_result_times(const string& fieldname)
{
    vector<double> times; // timesteps when the solution is written

    if (this->next_result())
        times.push_back(this->time());
    else
        dserror("no solution found in field '%s'", fieldname.c_str());

    while (this->next_result())
        times.push_back(this->time());
    return times;
}

/*----------------------------------------------------------------------*
 * get timesteps when the specific solution vector >name< is written
 *                                                               gjb02/08
 *----------------------------------------------------------------------*/
vector<double> PostResult::get_result_times(
        const string& fieldname,
        const string& groupname)
{
    vector<double> times; // timesteps when the solution is written

    if (this->next_result(groupname))
        times.push_back(this->time());
    else
        dserror("no solution found in field '%s'", fieldname.c_str());

    while (this->next_result(groupname))
        times.push_back(this->time());
    return times;
}

/*----------------------------------------------------------------------*
 * loads the next result block and opens new result files if there are
 * any. Returns 1 when a new result block has been found, otherwise
 * returns 0
 *----------------------------------------------------------------------*/
int PostResult::next_result()
{
  PostProblem* problem = field_->problem();
  int ret = 0;

  for (int i=pos_+1; i<problem->num_results(); ++i)
  {
    INT step;
    MAP* map;
    map = (*problem->result_groups())[problem->num_results()-1-i];

    if (match_field_result(map))
    {

      /*
       * Open the new files if there are any.
       *
       * If one of these files is here the other one has to be
       * here, too. If it's not, it's a bug in the input. */
      if ((map_symbol_count(map, "result_file") > 0))
      {
	close_result_files();
        open_result_files(map);
      }

      /*
       * We use the real step numbers here. That is a user has to give
       * the real numbers, too. Maybe that's the best way to handle
       * it. */
      /* In case of FSI everything else hurts even more. */
      step = map_read_int(map, "step");

      /* we are only interessted if the result matches the slice */
      if ((step >= problem->start()) &&
          ((step <= problem->end()) || (problem->end() == -1)) &&
          ((step - problem->start()) % problem->step() == 0))
      {
        pos_ = i;
        group_ = map;
        ret = 1;
        break;
      }
    }
  }
  return ret;
}


/*----------------------------------------------------------------------*
 * loads the next result block that contains written result values
 * specified by a given groupname. Returns 1 when a new result block has
 * been found, otherwise returns 0                              gjb 02/08
 *----------------------------------------------------------------------*/
int PostResult::next_result(const string& groupname)
{
    int ret = next_result();
    // go on, until the specified result is contained or end of time slice reached
    while((!map_has_map(group_, groupname.c_str())) && (ret == 1))
    {
        ret = next_result();
    }
    return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether a given result group belongs to this result.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
int PostResult::match_field_result(MAP* result_group) const
{
  return (strcmp(map_read_string(result_group, "field"),
                 fieldnames[field_->type()]) == 0) &&
    (map_read_int(result_group, "field_pos") == field_->field_pos()) &&
    (map_read_int(result_group, "discretization") == field_->disnum());
}

/*----------------------------------------------------------------------*
 * closes all the currently open result files
 *----------------------------------------------------------------------*/
void PostResult::close_result_files()
{
  file_.Close();
}

/*----------------------------------------------------------------------*
 * opens result files. The name is taken from the "result_file" entry
 * in the block 'field_info'
 *----------------------------------------------------------------------*/
void PostResult::open_result_files(MAP* field_info)
{
  int num_output_procs;
  if (!map_find_int(field_info,"num_output_proc",&num_output_procs))
  {
    num_output_procs = 1;
  }
  string basename = map_read_string(field_info,"result_file");
  //field_->problem()->set_basename(basename);
  file_.Open(basename,num_output_procs);
}

/*----------------------------------------------------------------------*
 * reads the data of the result vector 'name' from the current result
 * block and returns it as an Epetra Vector.
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> PostResult::read_result(const string name)
{
  MAP* result = map_read_map(group_, name.c_str());
  int columns;
  if (map_find_int(result,"columns",&columns))
  {
    if (columns != 1)
      dserror("got multivector with name '%s', vector expected", name.c_str());
  }
  return Teuchos::rcp_dynamic_cast<Epetra_Vector>(read_multi_result(name));
}

/*----------------------------------------------------------------------*
 * reads the data of the result vector 'name' from the current result
 * block and returns it as an std::vector<char>. the corresponding
 * elemap is returned, too.
 *----------------------------------------------------------------------*/
RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > >
PostResult::read_result_serialdensematrix(const string name)
{
  RefCountPtr<Epetra_Comm> comm = field_->problem()->comm();
  MAP* result = map_read_map(group_, name.c_str());
  string id_path = map_read_string(result, "ids");
  string value_path = map_read_string(result, "values");
  int columns = map_find_int(result,"columns",&columns);
  if (not map_find_int(result,"columns",&columns))
  {
    columns = 1;
  }
  if (columns != 1)
    dserror("got multivector with name '%s', vector<char> expected", name.c_str());

  RefCountPtr<Epetra_Map> elemap;
  RefCountPtr<std::vector<char> > data = file_.ReadResultDataVecChar(id_path, value_path, columns,
                                                                     *comm, elemap);

  RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > mapdata = rcp(new std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> >);
  int position=0;
//   cout << "elemap:\n" << *elemap << endl;
//   cout << "myelenum: " << elemap->NumMyElements() << endl;
  for (int i=0;i<elemap->NumMyElements();++i)
  {
    RefCountPtr<Epetra_SerialDenseMatrix> gpstress = rcp(new Epetra_SerialDenseMatrix);
    DRT::ParObject::ExtractfromPack(position, *data, *gpstress);
    (*mapdata)[elemap->GID(i)]=gpstress;
  }

  return mapdata;
}

/*----------------------------------------------------------------------*
 * reads the data of the result vector 'name' from the current result
 * block and returns it as an Epetra Vector.
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_MultiVector> PostResult::read_multi_result(const string name)
{
  RefCountPtr<Epetra_Comm> comm = field_->problem()->comm();
  MAP* result = map_read_map(group_, name.c_str());
  string id_path = map_read_string(result, "ids");
  string value_path = map_read_string(result, "values");
  int columns;
  if (not map_find_int(result,"columns",&columns))
  {
    columns = 1;
  }
  return file_.ReadResultData(id_path, value_path, columns, *comm);
}

#endif
