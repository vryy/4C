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
#ifdef TRILINOS_PACKAGE

#include <algorithm>
#include "post_drt_common.H"

#include "../io/hdf_reader.H"
#include <sstream>


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

// not actually used, but referenced in par_assignmesh.c
struct _PARTITION  *partition;
struct _SOLVAR  *solv;

char* fieldnames[] = FIELDNAMES;





/*----------------------------------------------------------------------*
 * The main part of this file. All the functions of the three classes
 * PostProbem, PostField and PostResult are defined here.
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

  /* We need to take two steps back. dsinit is too close to ccarat. */
  dstrc_exit();
  dstrc_exit();

  DSTraceHelper dst("PostProblem::PostProblem");
#endif

  string file;

  CLP.throwExceptions(false);
  CLP.setOption("start",&start_,"first time step to read");
  CLP.setOption("end",&end_,"last time step to read");
  CLP.setOption("step",&step_,"number of time steps to jump");
  CLP.setOption("file",&file,"control file to open");

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

  result_group_ = vector<MAP*>();
  setup_filter(file.c_str());

  ndim_ = map_read_int(&control_table_, "ndim");
  dsassert((ndim_ == 2) || (ndim_ == 3), "illegal dimension");

  // in the old version the problem type was checked here.


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
  DSTraceHelper dst("PostProblem::~PostProblem");
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
  DSTraceHelper dst("PostProblem::get_discretization");
  if (num >= static_cast<int>(fields_.size()))
    dserror("index for field vector out of range.(index=%i,vector_size=%i)",
            num,fields_.size());
  return &fields_[num];
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
void PostProblem::setup_filter(const char* output_name)
{
  DSTraceHelper dst("PostProblem::setup_filter");
  int length;
  string control_file_name;
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
  length = strlen(output_name);
  if ((length > 8) && (strcmp(output_name+length-8, ".control")==0))
  {
    control_file_name = output_name;
    basename_ = control_file_name.substr(0,length-8);
    strcpy(allfiles.outputfile_name, output_name);
    strcpy(allfiles.outputfile_name+length-8, ".post.log");
  }
  else
  {
    basename_ = output_name;
    control_file_name = basename_ + ".control";
    strcpy(allfiles.outputfile_name, output_name);
    strcpy(allfiles.outputfile_name+length, ".post.log");
  }
  allfiles.out_err = fopen(allfiles.outputfile_name, "w");

  parse_control_file(&control_table_, (char*)control_file_name.c_str());

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

  while (map_symbol_count(table, "restarted_run") > 0)
  {
    INT pos;
    CHAR* separator;
    FILE* f;
    SYMBOL* first_result;
    SYMBOL* previous_results;
    INT first_step;
    SYMBOL dummy_symbol;
    INT counter;

    /* copy directory information */
    separator = rindex(output_name, '/');
    if (separator != NULL)
    {
      pos = separator-output_name+1;
      control_file_name = output_name;
      control_file_name = control_file_name.substr(0,pos);
    }
    else
    {
      pos = 0;
      control_file_name = "../";
    }
    input_dir_ = control_file_name;

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
    parse_control_file(table, (char*)control_file_name.c_str());
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
  DSTraceHelper dst("PostProblem::read_mesh");
  SYMBOL* mesh = map_find_symbol(&control_table_,"field");
  if (mesh == NULL)
    dserror("No field found.");

  while (mesh != NULL)
  {
    int step;
    if (!map_find_int(symbol_map(mesh),"step",&step))
      dserror("No step information in field.");
    if (step == 0)
    {
      PostField currfield = getfield(symbol_map(mesh));

      int num_output_procs;
      if (!map_find_int(symbol_map(mesh),"num_output_proc",&num_output_procs))
      {
        num_output_procs = 1;
      }
      currfield.set_num_output_procs(num_output_procs);
      char* fn;
      if (!map_find_string(symbol_map(mesh),"mesh_file",&fn))
        dserror("No meshfile name for discretization %s.", currfield.discretization()->Name().c_str());
      string filename = fn;
      HDFReader reader = HDFReader(input_dir_,num_output_procs);
      reader.open(filename);

      RefCountPtr<vector<char> > node_data =
        reader.read_node_data(step, comm_->NumProc(), comm_->MyPID());
      currfield.discretization()->UnPackMyNodes(node_data);

      RefCountPtr<vector<char> > element_data =
        reader.read_element_data(step, comm_->NumProc(), comm_->MyPID());
      currfield.discretization()->UnPackMyElements(element_data);

      // before we can call discretization functions (in parallel) we
      // need field based communicators.
      create_communicators();

      distribute_drt_grids();
      currfield.discretization()->FillComplete();
      fields_.push_back(currfield);
    }
    mesh = mesh->next;
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
      vector<RefCountPtr<DRT::Discretization> >* discretization =
        new vector<RefCountPtr<DRT::Discretization> >();

      field[i].ccadis = (void*)discretization;
      field[i].dis = NULL;
      field[i].ndis = 0;
    }
    genprob.numfld = field_pos+1;
  }

  field[field_pos].fieldtyp = (FIELDTYP)type;

  vector<RefCountPtr<DRT::Discretization> >* discretization =
    (vector<RefCountPtr<DRT::Discretization> >*)field[field_pos].ccadis;
  if (static_cast<int>(discretization->size()) <= disnum)
  {
    field[field_pos].ndis = disnum+1;
    discretization->resize(disnum+1);
  }
  (*discretization)[disnum] = rcp(new DRT::Discretization(dis_name,comm_));
  return PostField((*discretization)[disnum], this, field_pos, field_name,
                   (FIELDTYP)type, disnum,numnd,numele);
}




/*----------------------------------------------------------------------*
 * The Constructor of PostField
 *----------------------------------------------------------------------*/
PostField::PostField(string name, RefCountPtr<Epetra_Comm> comm,
                     PostProblem* problem, int field_pos, string field_name,
                     FIELDTYP type, int disnum, int numnd, int numele):
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
                     int disnum, int numnd, int numele)
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
  file_((field->problem()->input_dir()),field->num_output_procs())
{
  DSTraceHelper("PostResult::PostResult");
  int num_output_proc = field_->num_output_procs();
  file_ = HDFReader((field_->problem()->input_dir()), num_output_proc);
}


/*----------------------------------------------------------------------*
 * The Destructor of PostResult
 *----------------------------------------------------------------------*/
PostResult::~PostResult()
{
  close_result_files();
}


/*----------------------------------------------------------------------*
 * loads the next result block and opens new result files if there are
 * any. Returns 1 when a new result block has been found, otherwise
 * returns 0
 *----------------------------------------------------------------------*/
int PostResult::next_result()
{
  DSTraceHelper("PostResult::next_result");
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



/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether a given result group belongs to this result.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
int PostResult::match_field_result(MAP* result_group)
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
  DSTraceHelper("PostResult::close_result_files");
  file_.close();
}

/*----------------------------------------------------------------------*
 * opens result files. The name is taken from the "result_file" entry
 * in the block 'field_info'
 *----------------------------------------------------------------------*/
void PostResult::open_result_files(MAP* field_info)
{
  DSTraceHelper("PostResult::open_result_files");
  string basename = map_read_string(field_info,"result_file");
  field_->problem()->set_basename(basename);
  file_.open(basename);
}

/*----------------------------------------------------------------------*
 * reads the data of the result vector 'name' from the current result
 * block and returns it as an Epetra Vector.
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Vector> PostResult::read_result(string name)
{
  RefCountPtr<Epetra_Comm> comm = field_->problem()->comm();
  int gId_num = field_->global_id_num();
  MAP* result = map_read_map(group_, const_cast<char*>(name.c_str()));
  string id_path = map_read_string(result, "ids");
  string value_path = map_read_string(result, "values");
  return file_.read_result_data(id_path.c_str(), value_path.c_str(), comm->NumProc(),
                                comm->MyPID(), comm, gId_num);
}

#endif
#endif
