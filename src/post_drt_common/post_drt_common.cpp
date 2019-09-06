/*----------------------------------------------------------------------*/
/*! \file

\brief drt binary filter library

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/



#include "post_drt_common.H"

#include <stack>
#include "../drt_lib/drt_exporter.H"
#include <EpetraExt_Transpose_CrsGraph.h>
#include "../drt_lib/drt_parobjectregister.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_periodicbc.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_inpar/inpar_problemtype.H"

#include "../pss_full/pss_cpp.h"
extern "C"
{
#include "../pss_full/pss_table_iter.h"
}

/*----------------------------------------------------------------------*
 * The main part of this file. All the functions of the three classes
 * PostProblem, PostField and PostResult are defined here.
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 * the Constructor of PostProblem
 *----------------------------------------------------------------------*/
PostProblem::PostProblem(Teuchos::CommandLineProcessor& CLP, int argc, char** argv)
    : start_(0), end_(-1), step_(1), mortar_(false)
{
#ifdef PARALLEL
  MPI_Init(&argc, &argv);
#endif

  std::string file = "xxx";
  std::string output;
  filter_ = "ensight";
  struct_vel_acc_ = "no";
  struct_mat_disp_ = "no";
  struct_rot_ = "no";
  std::string mortar_string = "no";

  int printparobjecttypes = 0;

  CLP.throwExceptions(false);
  CLP.setOption("filter", &filter_, "filter to run [ensight, gid, vtu, vtu_node_based, vti]");
  CLP.setOption("start", &start_, "first time step to read");
  CLP.setOption("end", &end_, "last time step to read");
  CLP.setOption("step", &step_, "number of time steps to jump");
  CLP.setOption("file", &file, "control file to open");
  CLP.setOption("output", &output, "output file name [defaults to control file name]");
  CLP.setOption("stresstype", &stresstype_,
      "stress output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption("stress", &stresstype_,
      "stress output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption("straintype", &straintype_,
      "strain output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption("strain", &straintype_,
      "strain output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption("mortar", &mortar_string, "Do post-processing of mortar interfaces [yes]");
  CLP.setOption("optquantitytype", &optquantitytype_,
      "optional quantity output type [cxyz, ndxyz, cxyz_ndxyz]");
  CLP.setOption(
      "optquantity", &optquantitytype_, "optional quantity output type [cxyz, ndxyz, cxyz_ndxyz]");
  CLP.setOption("heatfluxtype", &heatfluxtype_,
      "heatflux output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption("heatflux", &heatfluxtype_,
      "heatflux output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption("tempgradtype", &tempgradtype_,
      "tempgrad output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption("tempgrad", &tempgradtype_,
      "tempgrad output type [cxyz, ndxyz, cxyz_ndxyz, c123, nd123, c123_nd123]");
  CLP.setOption(
      "structvelacc", &struct_vel_acc_, "structural velocity and acceleration output [yes]");
  CLP.setOption("rotation", &struct_rot_, "structural rotation matrix R [yes]");
  CLP.setOption("structmatdisp", &struct_mat_disp_, "material displacement output output [yes]");
  CLP.setOption("printparobjecttypes", &printparobjecttypes,
      "print names of parobject types (registration hack)");
  CLP.setOption("outputtype", &outputtype_,
      "binary (bin) or ascii (ascii) output, option works for vtu filter only");
  Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = CLP.parse(argc, argv);

  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
  {
    exit(1);
  }

  if (printparobjecttypes)
  {
    // hack so that the parobject types are registered
    PrintParObjectList();
    exit(0);
  }

  if (file == "")
  {
    CLP.printHelpMessage(argv[0], std::cout);
    exit(1);
  }

  if (file.length() <= 8 or file.substr(file.length() - 8, 8) != ".control")
  {
    file += ".control";
  }

  if (output == "")
  {
    output = file.substr(0, file.length() - 8);
  }

  if (stresstype_ == "")
  {
    stresstype_ = "none";
  }

  if (straintype_ == "")
  {
    straintype_ = "none";
  }

  if (optquantitytype_ == "")
  {
    optquantitytype_ = "none";
  }

  if (heatfluxtype_ == "")
  {
    heatfluxtype_ = "none";
  }

  if (tempgradtype_ == "")
  {
    tempgradtype_ = "none";
  }

  if (outputtype_ == "")
  {
    outputtype_ = "bin";
  }

  if (mortar_string == "yes")
  {
    mortar_ = true;
  }

  result_group_ = std::vector<MAP*>();
  setup_filter(file, output);

  ndim_ = map_read_int(&control_table_, "ndim");
  dsassert((ndim_ == 1) || (ndim_ == 2) || (ndim_ == 3), "illegal dimension");

  const char* type = map_read_string(&control_table_, "problem_type");
  const std::string probtype(type);
  problemtype_ = INPAR::PROBLEMTYPE::StringToProblemType(probtype);

  spatial_approx_ = map_read_string(&control_table_, "spatial_approximation");

  if (spatial_approx_ != "Nurbs" and spatial_approx_ != "Polynomial" and
      spatial_approx_ != "Meshfree" and spatial_approx_ != "HDG")
  {
    dserror("unknown type of spatial approximation '%s'", spatial_approx_.c_str());
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
PostField* PostProblem::get_discretization(const int num)
{
  if (num >= static_cast<int>(fields_.size()))
  {
    std::cout << "You asked for discretization " << num
              << " (counting from zero), but there are only " << fields_.size()
              << " discretization(s)!";
    dserror("This is a bug!");
  }
  return &fields_[num];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int PostProblem::field_pos(const PostField* field) const
{
  for (std::vector<PostField>::const_iterator i = fields_.begin(); i != fields_.end(); ++i)
  {
    if (&*i == field)
    {
      return field - &fields_[0];
    }
  }
  dserror("field not in list");
  return -1;
}


/*----------------------------------------------------------------------*
 * returns the Epetra Communicator object
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Comm> PostProblem::comm() { return comm_; }


/*----------------------------------------------------------------------*
 * initializes all the data a filter needs. This function is called by
 * the Constructor.  (private)
 *----------------------------------------------------------------------*/
void PostProblem::setup_filter(std::string control_file_name, std::string output_name)
{
  MAP temp_table;

  comm_ = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));

  /* The warning system is not set up. It's rather stupid anyway. */

  basename_ = control_file_name.substr(0, control_file_name.length() - 8);
  outname_ = output_name;

  parse_control_file(&control_table_, control_file_name.c_str(), MPI_COMM_WORLD);

  /*
   * Now that we've read the control file given by the user we have to
   * take care of any previous (restarted) control files. These files
   * build a chain. So as long as a previous file exists we have to
   * open it and read any groups of results with smaller step numbers
   * than the results we've already read. */

  /* The general idea is to merge the different control files in
   * memory. If one step is written several times the last version is
   * used. */

  MAP* table = &control_table_;

  /* copy directory information */
  std::string::size_type separator = basename_.rfind('/', std::string::npos);
  if (separator != std::string::npos)
  {
    input_dir_ = basename_.substr(0, separator + 1);
  }
  else
  {
    input_dir_ = "";
  }

  while (map_symbol_count(table, "restarted_run") > 0)
  {
    /* copy directory information */
    control_file_name = input_dir_;

    /* copy file name */
    control_file_name += map_read_string(table, "restarted_run");
    control_file_name += ".control";

    /* test open to see if it exists */
    FILE* f = fopen(control_file_name.c_str(), "rb");
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
    SYMBOL* first_result = map_find_symbol(&control_table_, "result");
    if (first_result == NULL)
    {
      dserror("no result sections in control file '%s'\n", control_file_name.c_str());
    }
    while (first_result->next != NULL)
    {
      first_result = first_result->next;
    }
    const INT first_step = map_read_int(symbol_map(first_result), "step");


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
    parse_control_file(table, control_file_name.c_str(), MPI_COMM_WORLD);
    if (comm_->MyPID() == 0) printf("read restarted control file: %s\n", control_file_name.c_str());

    /* find the previous results */
    {
      INT counter = 0;

      /*
       * the dummy_symbol is a hack that allows us to treat all results
       * in the list the same way (use the same code). Without it we'd
       * need special conditions for the first entry. */
      SYMBOL dummy_symbol;
      SYMBOL* previous_results = &dummy_symbol;
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
          map_prepend_symbols(
              &control_table_, "result", result, map_symbol_count(table, "result") - counter);
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

    /* find the previous mesh files */
    {
      INT counter = 0;

      /*
       * the dummy_symbol is a hack that allows us to treat all results
       * in the list the same way (use the same code). Without it we'd
       * need special conditions for the first entry. */
      SYMBOL dummy_symbol;
      SYMBOL* previous_fields = &dummy_symbol;
      previous_fields->next = map_find_symbol(table, "field");
      while (previous_fields->next != NULL)
      {
        SYMBOL* field;
        INT step;
        field = previous_fields->next;
        step = map_read_int(symbol_map(field), "step");

        if (step < first_step)
        {
          /* found it */
          /* Now we simply switch all previous fields to our main
           * map. The assumption is a perfect ordering */
          map_prepend_symbols(
              &control_table_, "field", field, map_symbol_count(table, "field") - counter);
          previous_fields->next = NULL;

          /*
           * In case all fields go to the main map we have to disconnect
           * them explicitly. */
          if (previous_fields == &dummy_symbol)
          {
            map_disconnect_symbols(table, "field");
          }
          break;
        }

        /* Not found yet. Go up one field. */
        previous_fields = previous_fields->next;
        counter += 1;
      }
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
  SYMBOL* mesh = map_find_symbol(&control_table_, "field");
  if (mesh == NULL) dserror("No field found.");

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
    if (map_find_symbol(symbol_map(mesh), "mesh_file") != NULL)
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

    std::string name = map_read_string(meshmap, "field");

    bool havefield = false;
    for (unsigned i = 0; i < fields_.size(); ++i)
    {
      if (fields_[i].name() == name)
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
      if (!map_find_int(meshmap, "step", &step)) dserror("No step information in field.");

      PostField currfield = getfield(meshmap);

      int num_output_procs;
      if (!map_find_int(meshmap, "num_output_proc", &num_output_procs))
      {
        num_output_procs = 1;
      }
      currfield.set_num_output_procs(num_output_procs);
      char* fn;
      if (!map_find_string(meshmap, "mesh_file", &fn))
        dserror(
            "No meshfile name for discretization %s.", currfield.discretization()->Name().c_str());
      std::string filename = fn;
      IO::HDFReader reader = IO::HDFReader(input_dir_);
      reader.Open(filename, num_output_procs, comm_->NumProc(), comm_->MyPID());

      if (currfield.num_nodes() != 0)
      {
        Teuchos::RCP<std::vector<char>> node_data =
            reader.ReadNodeData(step, comm_->NumProc(), comm_->MyPID());
        currfield.discretization()->UnPackMyNodes(node_data);
      }

      if (currfield.num_elements() != 0)
      {
        Teuchos::RCP<std::vector<char>> element_data =
            reader.ReadElementData(step, comm_->NumProc(), comm_->MyPID());
        currfield.discretization()->UnPackMyElements(element_data);
      }

      Teuchos::RCP<std::vector<char>> cond_pbcsline;
      Teuchos::RCP<std::vector<char>> cond_pbcssurf;

      for (SYMBOL* condition = map_find_symbol(meshmap, "condition"); condition != NULL;
           condition = condition->next)
      {
        char* condname;
        if (not symbol_get_string(condition, &condname)) dserror("condition name expected");

        // read periodic boundary conditions if available
        if (std::string(condname) == "LinePeriodic")
        {
          DRT::Exporter exporter(*comm_);

          cond_pbcsline = Teuchos::rcp(new std::vector<char>());

          if (comm_->MyPID() == 0)
          {
            cond_pbcsline =
                reader.ReadCondition(step, comm_->NumProc(), comm_->MyPID(), "LinePeriodic");

#ifdef PARALLEL
            // distribute condition to all procs
            if (comm_->NumProc() > 1)
            {
              MPI_Request request;
              int tag = -1;
              int frompid = 0;
              int topid = -1;

              for (int np = 1; np < comm_->NumProc(); ++np)
              {
                tag = np;
                topid = np;

                exporter.ISend(
                    frompid, topid, &((*cond_pbcsline)[0]), (*cond_pbcsline).size(), tag, request);
              }
            }
#endif
          }
          else
          {
#ifdef PARALLEL
            int length = -1;
            int frompid = 0;
            int mypid = comm_->MyPID();

            std::vector<char> rblock;

            exporter.ReceiveAny(frompid, mypid, rblock, length);

            *cond_pbcsline = rblock;
#endif
          }

          currfield.discretization()->UnPackCondition(cond_pbcsline, "LinePeriodic");
        }
        else if (std::string(condname) == "SurfacePeriodic")
        {
          DRT::Exporter exporter(*comm_);

          cond_pbcssurf = Teuchos::rcp(new std::vector<char>());

          if (comm_->MyPID() == 0)
          {
            cond_pbcssurf =
                reader.ReadCondition(step, comm_->NumProc(), comm_->MyPID(), "SurfacePeriodic");

            // distribute condition to all procs
            if (comm_->NumProc() > 1)
            {
#ifdef PARALLEL
              MPI_Request request;
              int tag = -1;
              int frompid = 0;
              int topid = -1;

              for (int np = 1; np < comm_->NumProc(); ++np)
              {
                tag = np;
                topid = np;

                exporter.ISend(
                    frompid, topid, &((*cond_pbcssurf)[0]), (*cond_pbcssurf).size(), tag, request);
              }
#endif
            }
          }
          else
          {
#ifdef PARALLEL
            int length = -1;
            int frompid = 0;
            int mypid = comm_->MyPID();

            std::vector<char> rblock;

            exporter.ReceiveAny(frompid, mypid, rblock, length);

            *cond_pbcssurf = rblock;
#endif
          }
          currfield.discretization()->UnPackCondition(cond_pbcssurf, "SurfacePeriodic");
        }
        else
        {
          if (comm_->MyPID() == 0)
            printf("condition name '%s' not supported, continue anyway", condname);
          //          dserror("condition name '%s' not supported", condname);
        }
      }

      if ((currfield.problem()->Problemtype() == prb_particle_old or
              currfield.problem()->Problemtype() == prb_cavitation) and
          currfield.discretization()->Name() == "particle")
      {
        Teuchos::RCP<DRT::Discretization> particledis = currfield.discretization();
        // replace dof set to have an IndependentDofSet for particles
        Teuchos::RCP<DRT::IndependentDofSet> independentdofset =
            Teuchos::rcp(new DRT::IndependentDofSet(true));
        particledis->ReplaceDofSet(independentdofset);
        if (comm_->MyPID() == 0)
        {
          // find maximum number of possible particles during simulation
          int maxnodeid = get_max_nodeid("particle");
          // insert maxnodeid+1 particles into discret (on proc 0) placed at the origin
          double coords[] = {0.0, 0.0, 0.0};
          for (int id = 0; id <= maxnodeid; ++id)
          {
            // add a rigid shpere element for each particle (element id euqal to node id)
            Teuchos::RCP<DRT::Element> ele = DRT::UTILS::Factory("RIGIDSPHERE", "dummy", id, 0);
            particledis->AddElement(ele);

            Teuchos::RCP<DRT::Node> particle = Teuchos::rcp(new DRT::Node(id, coords, 0));
            particledis->AddNode(particle);
            ele->SetNodeIds(1, &id);
          }
          currfield.set_num_nodes(maxnodeid + 1);
        }
      }

      // read knot vectors for nurbs discretisations
      if (spatial_approx_ == "Nurbs")
      {
        // try a dynamic cast of the discretisation to a nurbs discretisation
        DRT::NURBS::NurbsDiscretization* nurbsdis =
            dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*currfield.discretization()));

        if (nurbsdis == NULL)
          dserror("Discretization %s is not a NurbsDiscretization",
              currfield.discretization()->Name().c_str());

        Teuchos::RCP<std::vector<char>> packed_knots;
        if (comm_->MyPID() == 0)
          packed_knots = reader.ReadKnotvector(step);
        else
          packed_knots = Teuchos::rcp(new std::vector<char>());

        // distribute knots to all procs
        if (comm_->NumProc() > 1)
        {
          DRT::Exporter exporter(nurbsdis->Comm());

          if (comm_->MyPID() == 0)
          {
#ifdef PARALLEL
            MPI_Request request;
            int tag = -1;
            int frompid = 0;
            int topid = -1;

            for (int np = 1; np < comm_->NumProc(); ++np)
            {
              tag = np;
              topid = np;

              exporter.ISend(
                  frompid, topid, &((*packed_knots)[0]), (*packed_knots).size(), tag, request);
            }
#endif
          }
          else
          {
#ifdef PARALLEL
            int length = -1;
            int frompid = 0;
            int mypid = comm_->MyPID();

            std::vector<char> rblock;

            exporter.ReceiveAny(frompid, mypid, rblock, length);

            *packed_knots = rblock;
#endif
          }
        }

        Teuchos::RCP<DRT::NURBS::Knotvector> knots = Teuchos::rcp(new DRT::NURBS::Knotvector());

        knots->Unpack(*packed_knots);

        if (nurbsdis == NULL)
        {
          dserror("expected a nurbs discretisation for spatial approx. Nurbs\n");
        }

        if (nurbsdis->Comm().NumProc() != 1)
          nurbsdis->SetupGhostingWrongNameDoNotUse(false, false, false);
        else
          nurbsdis->FillComplete(false, false, false);


        if (!(nurbsdis->Filled()))
        {
          dserror("nurbsdis was not fc\n");
        }

        int smallest_gid_in_dis = nurbsdis->ElementRowMap()->MinAllGID();

        knots->FinishKnots(smallest_gid_in_dis);

        nurbsdis->SetKnotVector(knots);

        // do initialisation
        currfield.discretization()->FillComplete();
      }
      else
      {
        // setup of parallel layout: create ghosting of already distributed nodes+elems
        if (currfield.discretization()->Comm().NumProc() != 1)
        {
          currfield.discretization()->SetupGhostingWrongNameDoNotUse();
        }
        else
          currfield.discretization()->FillComplete();
      }

      // -------------------------------------------------------------------
      // connect degrees of freedom for periodic boundary conditions
      // -------------------------------------------------------------------
      // parallel execution?!
      if ((cond_pbcssurf != Teuchos::null and not cond_pbcssurf->empty()) or
          (cond_pbcsline != Teuchos::null and not cond_pbcsline->empty()))
      {
        PeriodicBoundaryConditions pbc(currfield.discretization());
        pbc.UpdateDofsForPeriodicBoundaryConditions();
      }

      fields_.push_back(currfield);
    }
  }
}

/*----------------------------------------------------------------------*
 * in case of changing geometries, the mesh has to be read in each step
 *----------------------------------------------------------------------*/
void PostProblem::re_read_mesh(int fieldpos, std::string fieldname, int outputstep)
{
  SYMBOL* mesh = map_find_symbol(&control_table_, "field");
  if (mesh == NULL) dserror("No field found.");

  // We have to reverse the traversal of meshes we get from the control file
  // in order to get the same dof numbers in all discretizations as we had
  // during our calculation.
  // The order inside the control file is important!
  // Discretizations have to be FillComplete()ed in the same order as during
  // the calculation!

  MAP* meshmap = symbol_map(mesh);
  bool found = false;

  // search for the desired time step
  while (found == false and mesh != NULL)
  {
    // only those fields with a mesh file entry are readable here
    if (map_find_symbol(meshmap, "mesh_file") != NULL)
    {
      int step;
      if (!map_find_int(meshmap, "step", &step)) dserror("No step information in field.");

      std::string name = map_read_string(meshmap, "field");

      if (step == outputstep and not name.compare(fieldname))
      {
        found = true;
      }
      else
      {
        mesh = mesh->next;
        meshmap = symbol_map(mesh);
      }
    }
  }
  mesh = NULL;
  if (found == false) dserror("mesh for desired time step (%d) not found", outputstep);

  {
    int step;
    if (!map_find_int(meshmap, "step", &step)) dserror("No step information in field.");

    PostField currfield = getfield(meshmap);

    int num_output_procs;
    if (!map_find_int(meshmap, "num_output_proc", &num_output_procs))
    {
      num_output_procs = 1;
    }
    currfield.set_num_output_procs(num_output_procs);
    char* fn;
    if (!map_find_string(meshmap, "mesh_file", &fn))
      dserror(
          "No meshfile name for discretization %s.", currfield.discretization()->Name().c_str());
    std::string filename = fn;
    IO::HDFReader reader = IO::HDFReader(input_dir_);
    reader.Open(filename, num_output_procs, comm_->NumProc(), comm_->MyPID());

    Teuchos::RCP<std::vector<char>> node_data =
        reader.ReadNodeData(step, comm_->NumProc(), comm_->MyPID());
    currfield.discretization()->UnPackMyNodes(node_data);

    if (currfield.num_elements() != 0)
    {
      Teuchos::RCP<std::vector<char>> element_data =
          reader.ReadElementData(step, comm_->NumProc(), comm_->MyPID());
      currfield.discretization()->UnPackMyElements(element_data);
    }

    Teuchos::RCP<std::vector<char>> cond_pbcsline;
    Teuchos::RCP<std::vector<char>> cond_pbcssurf;

    for (SYMBOL* condition = map_find_symbol(meshmap, "condition"); condition != NULL;
         condition = condition->next)
    {
      char* condname;
      if (not symbol_get_string(condition, &condname)) dserror("condition name expected");

      // read periodic boundary conditions if available
      if (std::string(condname) == "LinePeriodic")
      {
        dserror("periodic boundary conditions with changing geometries have never been tested");

        DRT::Exporter exporter(*comm_);

        cond_pbcsline = Teuchos::rcp(new std::vector<char>());

        if (comm_->MyPID() == 0)
        {
          cond_pbcsline =
              reader.ReadCondition(step, comm_->NumProc(), comm_->MyPID(), "LinePeriodic");

#ifdef PARALLEL
          // distribute condition to all procs
          if (comm_->NumProc() > 1)
          {
            MPI_Request request;
            int tag = -1;
            int frompid = 0;
            int topid = -1;

            for (int np = 1; np < comm_->NumProc(); ++np)
            {
              tag = np;
              topid = np;

              exporter.ISend(
                  frompid, topid, &((*cond_pbcsline)[0]), (*cond_pbcsline).size(), tag, request);
            }
          }
#endif
        }
        else
        {
#ifdef PARALLEL
          int length = -1;
          int frompid = 0;
          int mypid = comm_->MyPID();

          std::vector<char> rblock;

          exporter.ReceiveAny(frompid, mypid, rblock, length);

          *cond_pbcsline = rblock;
#endif
        }

        currfield.discretization()->UnPackCondition(cond_pbcsline, "LinePeriodic");
      }
      else if (std::string(condname) == "SurfacePeriodic")
      {
        DRT::Exporter exporter(*comm_);

        cond_pbcssurf = Teuchos::rcp(new std::vector<char>());

        if (comm_->MyPID() == 0)
        {
          cond_pbcssurf =
              reader.ReadCondition(step, comm_->NumProc(), comm_->MyPID(), "SurfacePeriodic");

          // distribute condition to all procs
          if (comm_->NumProc() > 1)
          {
#ifdef PARALLEL
            MPI_Request request;
            int tag = -1;
            int frompid = 0;
            int topid = -1;

            for (int np = 1; np < comm_->NumProc(); ++np)
            {
              tag = np;
              topid = np;

              exporter.ISend(
                  frompid, topid, &((*cond_pbcssurf)[0]), (*cond_pbcssurf).size(), tag, request);
            }
#endif
          }
        }
        else
        {
#ifdef PARALLEL
          int length = -1;
          int frompid = 0;
          int mypid = comm_->MyPID();

          std::vector<char> rblock;

          exporter.ReceiveAny(frompid, mypid, rblock, length);

          *cond_pbcssurf = rblock;
#endif
        }
        currfield.discretization()->UnPackCondition(cond_pbcssurf, "SurfacePeriodic");
      }
      else
      {
        if (comm_->MyPID() == 0)
          printf("condition name '%s' not supported, continue anyway", condname);
        //          dserror("condition name '%s' not supported", condname);
      }
    }

    // read knot vectors for nurbs discretisations
    if (spatial_approx_ == "Nurbs")
    {
      // try a dynamic cast of the discretisation to a nurbs discretisation
      DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*currfield.discretization()));

      if (nurbsdis == NULL)
        dserror("Discretization %s is not a NurbsDiscretization",
            currfield.discretization()->Name().c_str());

      Teuchos::RCP<std::vector<char>> packed_knots;
      if (comm_->MyPID() == 0)
        packed_knots = reader.ReadKnotvector(step);
      else
        packed_knots = Teuchos::rcp(new std::vector<char>());

      // distribute knots to all procs
      if (comm_->NumProc() > 1)
      {
        DRT::Exporter exporter(nurbsdis->Comm());

        if (comm_->MyPID() == 0)
        {
#ifdef PARALLEL
          MPI_Request request;
          int tag = -1;
          int frompid = 0;
          int topid = -1;

          for (int np = 1; np < comm_->NumProc(); ++np)
          {
            tag = np;
            topid = np;

            exporter.ISend(
                frompid, topid, &((*packed_knots)[0]), (*packed_knots).size(), tag, request);
          }
#endif
        }
        else
        {
#ifdef PARALLEL
          int length = -1;
          int frompid = 0;
          int mypid = comm_->MyPID();

          std::vector<char> rblock;

          exporter.ReceiveAny(frompid, mypid, rblock, length);

          *packed_knots = rblock;
#endif
        }
      }

      Teuchos::RCP<DRT::NURBS::Knotvector> knots = Teuchos::rcp(new DRT::NURBS::Knotvector());

      knots->Unpack(*packed_knots);

      if (nurbsdis == NULL)
      {
        dserror("expected a nurbs discretisation for spatial approx. Nurbs\n");
      }

      if (nurbsdis->Comm().NumProc() != 1)
        nurbsdis->SetupGhostingWrongNameDoNotUse(false, false, false);
      else
        nurbsdis->FillComplete(false, false, false);


      if (!(nurbsdis->Filled()))
      {
        dserror("nurbsdis was not fc\n");
      }

      int smallest_gid_in_dis = nurbsdis->ElementRowMap()->MinAllGID();

      knots->FinishKnots(smallest_gid_in_dis);

      nurbsdis->SetKnotVector(knots);

      // do initialisation
      currfield.discretization()->FillComplete();
    }
    else
    {
      // setup of parallel layout: create ghosting of already distributed nodes+elems
      if (currfield.discretization()->Comm().NumProc() != 1)
      {
        currfield.discretization()->SetupGhostingWrongNameDoNotUse();
      }
      else
        currfield.discretization()->FillComplete();
    }

    // -------------------------------------------------------------------
    // connect degrees of freedom for periodic boundary conditions
    // -------------------------------------------------------------------
    // parallel execution?!
    if ((cond_pbcssurf != Teuchos::null and not cond_pbcssurf->empty()) or
        (cond_pbcsline != Teuchos::null and not cond_pbcsline->empty()))
    {
      PeriodicBoundaryConditions pbc(currfield.discretization());
      pbc.UpdateDofsForPeriodicBoundaryConditions();
    }

    // store newly created postfield and replace old one
    fields_[fieldpos] = currfield;
  }
}

/*----------------------------------------------------------------------*
 * creates and returns a PostField instance from a field MAP. (private)
 *----------------------------------------------------------------------*/
PostField PostProblem::getfield(MAP* field_info)
{
  const char* field_name = map_read_string(field_info, "field");
  const int numnd = map_read_int(field_info, "num_nd");
  const int numele = map_read_int(field_info, "num_ele");

  Teuchos::RCP<DRT::Discretization> dis;

  if (spatial_approx_ == "Polynomial" or spatial_approx_ == "Meshfree" or spatial_approx_ == "HDG")
  {
    dis = Teuchos::rcp(new DRT::Discretization(field_name, comm_));
  }
  else if (spatial_approx_ == "Nurbs")
  {
    dis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(field_name, comm_));
  }
  else
  {
    dserror("Unknown type of spatial approximation\n");
  }

  return PostField(dis, this, field_name, numnd, numele);
}

/*--------------------------------------------------------------------------*
 * loop all fields and get maximum node id for given field name ghamm 03/13
 *--------------------------------------------------------------------------*/
int PostProblem::get_max_nodeid(const std::string& fieldname)
{
  SYMBOL* mesh = map_find_symbol(&control_table_, "field");

  if (mesh == NULL) dserror("No field found.");

  int maxnodeid = -1;
  while (mesh != NULL)
  {
    MAP* meshmap = symbol_map(mesh);
    mesh = mesh->next;

    std::string name = map_read_string(meshmap, "field");

    if (name == fieldname)
    {
      const int stepmaxnodeid = map_read_int(meshmap, "max_nodeid");
      if (stepmaxnodeid > maxnodeid) maxnodeid = stepmaxnodeid;
    }
  }

  return maxnodeid;
}


/*----------------------------------------------------------------------*
 * Constructor of PostField.
 *----------------------------------------------------------------------*/
PostField::PostField(Teuchos::RCP<DRT::Discretization> dis, PostProblem* problem,
    std::string field_name, const int numnd, const int numele)
    : dis_(dis), problem_(problem), field_name_(field_name), numnd_(numnd), numele_(numele)
{
}

/*----------------------------------------------------------------------*
 * The Destructor
 *----------------------------------------------------------------------*/
PostField::~PostField() {}


/*----------------------------------------------------------------------*
 * The Constructor of PostResult
 *----------------------------------------------------------------------*/
PostResult::PostResult(PostField* field)
    : field_(field), pos_(-1), group_(NULL), file_((field->problem()->input_dir()))
{
}


/*----------------------------------------------------------------------*
 * The Destructor of PostResult
 *----------------------------------------------------------------------*/
PostResult::~PostResult() { close_result_files(); }

/*----------------------------------------------------------------------*
 * get timesteps when the solution is written
 *----------------------------------------------------------------------*/
std::vector<double> PostResult::get_result_times(const std::string& fieldname)
{
  std::vector<double> times;  // timesteps when the solution is written

#if 0
    if (this->next_result())
        times.push_back(this->time());
    else
        dserror("no solution found in field '%s'", fieldname.c_str());
#endif

  while (this->next_result()) times.push_back(this->time());

  if (times.size() == 0)
  {
    dserror(
        "PostResult::get_result_times(fieldname='%s'):\n  no solution steps found in specified "
        "timestep range! Check --start, --end, --step parameters.",
        fieldname.c_str());
  }

  return times;
}

/*----------------------------------------------------------------------*
 * get timesteps when the specific solution vector >name< is written
 *                                                               gjb02/08
 *----------------------------------------------------------------------*/
std::vector<double> PostResult::get_result_times(
    const std::string& fieldname, const std::string& groupname)
{
  std::vector<double> times;  // timesteps when the solution is written

  if (this->next_result(groupname))
    times.push_back(this->time());
  else
    dserror("no solution found in field '%s'", fieldname.c_str());

  while (this->next_result(groupname)) times.push_back(this->time());

  if (times.size() == 0)
  {
    dserror(
        "PostResult::get_result_times(fieldname='%s', groupname='%s'):\n  no solution steps found "
        "in specified timestep range! Check --start, --end, --step parameters.",
        fieldname.c_str(), groupname.c_str());
  }

  return times;
}

/*----------------------------------------------------------------------*
 * get times and steps when the solution is written
 *----------------------------------------------------------------------*/
void PostResult::get_result_timesandsteps(
    const std::string& fieldname, std::vector<double>& times, std::vector<int>& steps)
{
  while (this->next_result())
  {
    times.push_back(this->time());
    steps.push_back(this->step());
  }

  if (times.size() == 0 or steps.size() == 0)
  {
    dserror(
        "PostResult::get_result_timesandsteps(fieldname='%s'):\n  no solution steps found in "
        "specified range! Check --start, --end, --step parameters.",
        fieldname.c_str());
  }

  return;
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

  for (int i = pos_ + 1; i < problem->num_results(); ++i)
  {
    MAP* map = (*problem->result_groups())[problem->num_results() - 1 - i];

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
      const INT step = map_read_int(map, "step");

      /* we are only interessted if the result matches the slice */
      if ((step >= problem->start()) && ((step <= problem->end()) || (problem->end() == -1)) &&
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
int PostResult::next_result(const std::string& groupname)
{
  int ret = next_result();
  // go on, until the specified result is contained or end of time slice reached
  while ((!map_has_map(group_, groupname.c_str())) && (ret == 1))
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
  return field_->name() == map_read_string(result_group, "field");
}

/*----------------------------------------------------------------------*
 * closes all the currently open result files
 *----------------------------------------------------------------------*/
void PostResult::close_result_files() { file_.Close(); }

/*----------------------------------------------------------------------*
 * opens result files. The name is taken from the "result_file" entry
 * in the block 'field_info'
 *----------------------------------------------------------------------*/
void PostResult::open_result_files(MAP* field_info)
{
  int num_output_procs;
  if (!map_find_int(field_info, "num_output_proc", &num_output_procs))
  {
    num_output_procs = 1;
  }
  const std::string basename = map_read_string(field_info, "result_file");
  // field_->problem()->set_basename(basename);
  Epetra_Comm& comm = *field_->problem()->comm();
  file_.Open(basename, num_output_procs, comm.NumProc(), comm.MyPID());
}

/*----------------------------------------------------------------------*
 * reads the data of the result vector 'name' from the current result
 * block and returns it as an Epetra Vector.
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> PostResult::read_result(const std::string name)
{
  MAP* result = map_read_map(group_, name.c_str());
  int columns;
  if (map_find_int(result, "columns", &columns))
  {
    if (columns != 1) dserror("got multivector with name '%s', vector expected", name.c_str());
  }
  return Teuchos::rcp_dynamic_cast<Epetra_Vector>(read_multi_result(name));
}

/*----------------------------------------------------------------------*
 * reads the data of the result vector 'name' from the current result
 * block and returns it as an std::vector<char>. the corresponding
 * elemap is returned, too.
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>>
PostResult::read_result_serialdensematrix(const std::string name)
{
  Teuchos::RCP<Epetra_Comm> comm = field_->problem()->comm();
  MAP* result = map_read_map(group_, name.c_str());
  std::string id_path = map_read_string(result, "ids");
  std::string value_path = map_read_string(result, "values");
  int columns = map_find_int(result, "columns", &columns);
  if (not map_find_int(result, "columns", &columns))
  {
    columns = 1;
  }
  if (columns != 1)
    dserror("got multivector with name '%s', std::vector<char> expected", name.c_str());

  Teuchos::RCP<Epetra_Map> elemap;
  Teuchos::RCP<std::vector<char>> data =
      file_.ReadResultDataVecChar(id_path, value_path, columns, *comm, elemap);

  Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> mapdata =
      Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>);
  std::vector<char>::size_type position = 0;
  //   std::cout << "elemap:\n" << *elemap << std::endl;
  //   std::cout << "myelenum: " << elemap->NumMyElements() << std::endl;
  for (int i = 0; i < elemap->NumMyElements(); ++i)
  {
    Teuchos::RCP<Epetra_SerialDenseMatrix> gpstress = Teuchos::rcp(new Epetra_SerialDenseMatrix);
    DRT::ParObject::ExtractfromPack(position, *data, *gpstress);
    (*mapdata)[elemap->GID(i)] = gpstress;
  }

  const Epetra_Map& elecolmap = *field_->discretization()->ElementColMap();
  DRT::Exporter ex(*elemap, elecolmap, *comm);
  ex.Export(*mapdata);

  return mapdata;
}

/*----------------------------------------------------------------------*
 * reads the data of the result vector 'name' from the current result
 * block and returns it as an Epetra Vector.
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> PostResult::read_multi_result(const std::string name)
{
  const Teuchos::RCP<Epetra_Comm> comm = field_->problem()->comm();
  MAP* result = map_read_map(group_, name.c_str());
  const std::string id_path = map_read_string(result, "ids");
  const std::string value_path = map_read_string(result, "values");
  int columns;
  if (not map_find_int(result, "columns", &columns))
  {
    columns = 1;
  }
  return file_.ReadResultData(id_path, value_path, columns, *comm);
}

//! returns time of this result
double PostResult::time() const { return map_read_real(group_, "time"); }

//! returns step number of this result
int PostResult::step() const { return map_read_int(group_, "step"); }



//! returns the number of global Dof-Ids
int PostField::global_id_num() const { return dis_->DofRowMap()->NumGlobalElements(); }
