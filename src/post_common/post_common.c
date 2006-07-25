/*!
\file
\brief Code that is common to all filters.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Filters like this one are special inhabitants of the ccarat
world. They are always single processor applications yet they share
some code with ccarat and are closely linked to ccarat internals.

The general idea is that we cannot load the whole result data into
memory at once.

All filters need to read the ccarat output. Thus there are some common
functions.

\author u.kue
\date 10/04

*/

#include <assert.h>
#include <strings.h>
#include <stdarg.h>
#include "../headers/standardtypes.h"
#include "post_common.h"


/* There are some global variables in ccarat that are needed by the
 * service functions. We need to specify them here and set them up
 * properly. */
struct _FILES           allfiles;
struct _PAR     par;
static INT      num_para=6;
static INT      parameter[]={'l', 's','u', 'g', 'w', 'o'};
static INT      first=1;

#ifdef DEBUG
struct _CCA_TRACE         trace;
#endif


/*----------------------------------------------------------------------*/
/*!
  \brief All fields names.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
CHAR* fieldnames[] = FIELDNAMES;

/*----------------------------------------------------------------------*/
/*!
  \brief All dis type names.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
CHAR* distypenames[] = DISTYPENAMES;

/*----------------------------------------------------------------------*/
/*!
  \brief All element type names.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
CHAR* elementnames[] = ELEMENTNAMES;


/*----------------------------------------------------------------------*/
/*!
  \brief This is ccarat setup in a hurry.

  Initialize ds tracer. Open log and control file. Read control file.

  \param output_name   (i) control file name as given at the command line
  \param control_table (o) map for the control file's contents
  \param basename      (o) the control file's name

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void setup_filter(CHAR* output_name, MAP* control_table, CHAR* basename)
{
  INT length;
  CHAR control_file_name[100];
  MAP* table;
  MAP temp_table;

  par.myrank = 0;
  par.nprocs = 1;

  /* Do we want to have tracing? We could. */
#ifdef DEBUG
  dstrc_enter("setup_filter");
#endif

  /* The warning system is not set up. It's rather stupid anyway. */

  /* We need to open the error output file. The other ones are not
   * important. */
  length = strlen(output_name);
  if ((length > 8) && (strcmp(output_name+length-8, ".control")==0)) {
    /* dsassert isn't working yet. */
    assert(length-8+13 < 100);
    strcpy(allfiles.outputfile_name, output_name);
    strcpy(allfiles.outputfile_name+length-8, ".post.log");
    strcpy(control_file_name, output_name);
    strcpy(basename, output_name);
    basename[length-8] = '\0';
  }
  else {
    /* dsassert isn't working yet. */
    assert(length+13 < 100);
    strcpy(allfiles.outputfile_name, output_name);
    strcpy(allfiles.outputfile_name+length, ".post.log");
    strcpy(control_file_name, output_name);
    strcpy(control_file_name+length, ".control");
    strcpy(basename, output_name);
  }
  allfiles.out_err = fopen(allfiles.outputfile_name, "w");

  parse_control_file(control_table, control_file_name);

  /*
   * Now that we've read the control file given by the user we have to
   * take care of any previous (restarted) control files. These files
   * build a chain. So as long as a previous file exists we have to
   * open it and read any groups of results with smaller step numbers
   * than the results we've already read. */

  /* The general idea is to merge the different control files in
   * memory. If one step is written several times the last version is
   * used. */

  table = control_table;

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
      strncpy(control_file_name, output_name, pos);
    }
    else
    {
      pos = 0;
    }

    /* copy file name */
    strcpy(&control_file_name[pos], map_read_string(table, "restarted_run"));
    strcat(control_file_name, ".control");

    /* test open to see if it exists */
    f = fopen(control_file_name, "rb");
    if (f == NULL)
    {
      printf("Restarted control file '%s' does not exist. Skip previous results.\n",
             control_file_name);
      break;
    }
    fclose(f);

    /* copy all the result steps that are previous to this file */
    /* We assume that the results are ordered! */

    /*------------------------------------------------------------------*/
    /* find the first result in the current table */
    first_result = map_find_symbol(control_table, "result");
    if (first_result == NULL)
    {
      dserror("no result sections in control file '%s'\n", control_file_name);
    }
    while (first_result->next != NULL)
    {
      first_result = first_result->next;
    }
    first_step = map_read_int(symbol_map(first_result), "step");


    /*------------------------------------------------------------------*/
    /* done with this control file */
    if (table != control_table)
    {
      destroy_map(table);
    }

    /* The first time we reach this place we had just used the main
     * control table. But from now on we are interessted in the
     * previous control files we read. */
    table = &temp_table;

    /* read the previous control file */
    parse_control_file(table, control_file_name);
    printf("read restarted control file: %s\n", control_file_name);

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
        map_prepend_symbols(control_table, "result", result,
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

#ifdef DEBUG
  dstrc_exit();
#endif
}
DOUBLE linear_interpolation(INT x, INT y, INT z,  ELEMENT* element, INT k)
{
  DOUBLE result;

  result=-0.125*(1-x)*(1-y)*(1-z)*(2+x+y+z) * (element->node[0]->x[k]);
  result-=0.125*(1+x)*(1-y)*(1-z)*(2-x+y+z) * (element->node[1]->x[k]);
  result-=0.125*(1+x)*(1+y)*(1-z)*(2-x-y+z) * (element->node[2]->x[k]);
  result-=0.125*(1-x)*(1+y)*(1-z)*(2+x-y+z) * (element->node[3]->x[k]);
  result-=0.125*(1-x)*(1-y)*(1+z)*(2+x+y-z) * (element->node[4]->x[k]);
  result-=0.125*(1+x)*(1-y)*(1+z)*(2-x+y-z) * (element->node[5]->x[k]);
  result-=0.125*(1+x)*(1+y)*(1+z)*(2-x-y-z) * (element->node[6]->x[k]);
  result-=0.125*(1-x)*(1+y)*(1+z)*(2+x-y-z) * (element->node[7]->x[k]);
  result+=0.25*(1-x*x)*(1-y)*(1-z) * (element->node[8]->x[k]);
  result+=0.25*(1+x)*(1-y*y)*(1-z) * (element->node[9]->x[k]);
  result+=0.25*(1-x*x)*(1+y)*(1-z) * (element->node[10]->x[k]);
  result+=0.25*(1-x)*(1-y*y)*(1-z) * (element->node[11]->x[k]);
  result+=0.25*(1-x)*(1-y)*(1-z*z) * (element->node[12]->x[k]);
  result+=0.25*(1+x)*(1-y)*(1-z*z) * (element->node[13]->x[k]);
  result+=0.25*(1+x)*(1+y)*(1-z*z) * (element->node[14]->x[k]);
  result+=0.25*(1-x)*(1+y)*(1-z*z) * (element->node[15]->x[k]);
  result+=0.25*(1-x*x)*(1-y)*(1+z) * (element->node[16]->x[k]);
  result+=0.25*(1+x)*(1-y*y)*(1+z) * (element->node[17]->x[k]);
  result+=0.25*(1-x*x)*(1+y)*(1+z) * (element->node[18]->x[k]);
  result+=0.25*(1-x)*(1-y*y)*(1+z) * (element->node[19]->x[k]);

  return result;
}
/*------------------------------------------------------------*/
/* This function tests, if a generated node is yet existing.
 * If true, it updates the existing node,
 * if false, it creates a new node */
INT test_node(INT numnp,INT numnp_old, POST_DISCRETIZATION* discret, INT i, INT j, DOUBLE* x_temp, INT* nodelist)
{ INT k, m, n;
  DOUBLE f;

  if (first!=1)
  {
    for (m=numnp_old; m<numnp; m++)
    {
      if (nodelist[m-numnp_old]!=1)
      {
        n=0;
        for (k=0; k<3; k++)
        { f = fabs(discret->node[m].x[k]-x_temp[k]);
          if (f<EPS9)
            n++;
        }
        if (n==3)
        { for (k=0; k<3; k++)
          { f=discret->node[m].x[k]+x_temp[k];
            discret->node[m].x[k] = f/2;
          }
          discret->element[i].node[j] = &discret->node[m];
          discret->node[m].numele++;
          nodelist[m-numnp_old]=1;
          goto end;
        }
      }
    }
  }
  discret->node[numnp].Id_loc = numnp;
  discret->node[numnp].Id = numnp;
  discret->node[numnp].proc = 0;
  discret->node[numnp].numele = 1;
  for (k=0; k<3; k++)
    discret->node[numnp].x[k] = x_temp[k];
  discret->element[i].node[j] = &(discret->node[numnp]);
  numnp++;
  first=0;

  end:
  return numnp;
}

/*----------------------------------------------------------------------*/
/*!
  \brief A Hack.

  This a yet another hack. This function is called by dserror and
  closes all open files --- in ccarat. The filters are not that
  critical. Thus we do nothing here. We just have to have this
  function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void io_emergency_close_files()
{
}


/*----------------------------------------------------------------------*/
/*!
  \brief The log level of the filter

  Specifies what messages are to appear on the screen. Any messages
  are written to the log file as long as the log level is greater than
  zero.

  The value can be set via the command line option -l.

  \author u.kue
  \date 01/05
*/
/*----------------------------------------------------------------------*/
static INT loglevel = 2;


/*----------------------------------------------------------------------*/
/*!
  \brief Filter log output

  Writes the message to the log file. If the given level is small
  enough it prints the message on screen as well.

  \param level (i) level of this message
  \param msg   (i) any text followed by printf like arguments

  \author u.kue
  \date 01/05
*/
/*----------------------------------------------------------------------*/
void post_log(INT level, CHAR* msg, ...)
{
  va_list ap;

  if (loglevel > 0)
  {
    va_start(ap, msg);

    vfprintf(allfiles.out_err,msg,ap);
    if (level <= loglevel)
      vprintf(msg,ap);

    va_end(ap);
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Init a translation table.

  The purpose of these tables is to find the internal enum value to an
  external number read from the data file. We need to translate each
  enum value we read. Otherwise the number written would be an
  implementation detail, thus we could not change our implementation
  once we had some binary files.

  \param table (o) table to be initialized
  \param group (i) control file group that defines the external values
  \param names (i) list of names ordered by internal number, NULL terminated

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void init_translation_table(TRANSLATION_TABLE* table,
                            MAP* group,
                            CHAR** names)
{
  INT i;
  INT count;

#ifdef DEBUG
  dstrc_enter("init_translation_table");
#endif

  /*
   * It seems to be more convenient to require the list of names to be
   * null terminated and not to ask for it's length in advance. Thus
   * we have to count the names. But we can easily to this because
   * this function is never called in the inner loop. We've got the
   * time. ;) */
  count = 0;
  for (count=0; names[count] != NULL; ++count)
  {}

  table->group = group;
  table->table = CCACALLOC(count, sizeof(INT));

  /* set everything to -1. */
  memset(table->table, 0xff, count*sizeof(INT));

  table->length = count;
  for (i=0; i<count; ++i)
  {
    if (map_symbol_count(group, names[i]) > 0)
    {
      INT num;
      num = map_read_int(group, names[i]);
      /*if ((num < 0) || (num >= count))*/
      if (num < 0)
      {
        dserror("illegal external number for name '%s': %d", names[i], num);
      }
      if (num >= table->length)
      {
        table->table = CCAREALLOC(table->table, (num+1)*sizeof(INT));
        memset(table->table+table->length, 0xff, (num-table->length)*sizeof(INT));
        table->length = num+1;
      }

      /* extern -> intern translation */
      table->table[num] = i;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clear the memory occupied by the translation table.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void destroy_translation_table(TRANSLATION_TABLE* table)
{
#ifdef DEBUG
  dstrc_enter("destroy_translation_table");
#endif

  CCAFREE(table->table);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Open data files.

  \param result     (i/o) result those files are to be opened
  \param field_info (i)   group that contains the filename definitions
  \param prefix     (i)   file name prefix

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
static void open_data_files(RESULT_DATA* result, MAP* field_info, CHAR* prefix)
{
  CHAR* filename;
  CHAR var_name[100];
  PROBLEM_DATA* problem;

#ifdef DEBUG
  dstrc_enter("open_data_files");
#endif

  problem = result->field->problem;

  sprintf(var_name, "%s_value_file", prefix);
  filename = map_read_string(field_info, var_name);

  /* It's misleading to look in the current directory by default. */

  /*result->value_file = fopen(filename, "rb");
    if (result->value_file == NULL)*/
  {
    CHAR buf[100];
    strcpy(buf, problem->input_dir);
    strcat(buf, filename);

    /* windows asks for a binary flag here... */
    result->value_file = fopen(buf, "rb");
    if (result->value_file == NULL)
    {
      dserror("failed to open file '%s'", filename);
    }
#ifdef LOWMEM
    post_log(2, "open file: '%s'\n", buf);
#endif
  }
  sprintf(var_name, "%s_size_file", prefix);

  filename = map_read_string(field_info, var_name);
  /*result->size_file = fopen(filename, "rb");
    if (result->size_file == NULL)*/
  {
    CHAR buf[100];
    strcpy(buf, problem->input_dir);
    strcat(buf, filename);

    /* windows asks for a binary flag here... */
    result->size_file = fopen(buf, "rb");
    if (result->size_file == NULL)
    {
      dserror("failed to open file '%s'", filename);
    }
#ifdef LOWMEM
    post_log(2, "open file: '%s'\n", buf);
#endif
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract one discretization's data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_field_data(PROBLEM_DATA* problem, FIELD_DATA* field, MAP* field_info)
{
  INT i;
  CHAR* type;

#ifdef DEBUG
  dstrc_enter("init_field_data");
#endif

  field->problem = problem;

  field->group = field_info;
  field->field_pos = map_read_int(field_info, "field_pos");
  field->disnum = map_read_int(field_info, "discretization");
  field->numele = map_read_int(field_info, "numele");
  field->numnp = map_read_int(field_info, "numnp");
  field->numdf = map_read_int(field_info, "numdof");

  field->name = map_read_string(field_info, "field");
  for (i=0; fieldnames[i]!=NULL; ++i)
  {
    if (strcmp(fieldnames[i], field->name)==0)
    {
      field->type = i;
      break;
    }
  }
  if (fieldnames[i]==NULL)
  {
    dserror("unknown field type '%s'", type);
  }

  /*--------------------------------------------------------------------*/
  /* Open the data files. */

  /* The fake variables. */
  field->head.pos = -1;
  field->head.field = field;
  field->head.group = field_info;

  open_data_files(&(field->head), field_info, "mesh");

  /*--------------------------------------------------------------------*/
  /* setup chunk structures */

  init_chunk_data(&(field->head), &(field->ele_param), "ele_param");
  init_chunk_data(&(field->head), &(field->mesh), "mesh");
  init_chunk_data(&(field->head), &(field->coords), "coords");

  /*--------------------------------------------------------------------*/
  /* special problems demand special attention. */

#ifdef D_SHELL8
  if (map_has_string(field_info, "shell8_problem", "yes"))
  {
    field->is_shell8_problem = 1;
  }
  else
  {
    field->is_shell8_problem = 0;
  }
#endif

#ifdef D_SHELL9
  if (map_has_string(field_info, "shell9_problem", "yes"))
  {
    /* This is a shell9 problem. There is guaranteed to be just one
     * type of element. The element_type flags are ignored. */

    field->is_shell9_problem = 1;

    field->s9_smooth_results = map_has_string(field_info, "shell9_smoothed", "yes");
    field->s9_layers = map_read_int(field_info, "shell9_layers");
  }
  else {
    field->is_shell9_problem = 0;
  }
#endif
#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief Print a short message how to use the filter and exit.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
static void usage(CHAR* progname)
{
  CHAR *index;

  index = rindex(progname, '/');

  if (strcmp(&index[1], "post_visual3")==0)
  {
    printf("usage: %s [options] control-file\n", progname);
    printf("\n"
         " options:\n"
         "    -s beg:end[:step]        read from beg to end-1 every step\n"
         "    -l [level]               set log level (0==none, infty==everything)\n\n"
         "    -w                       set white background for movie creation\n"
         "    -g                       set grey colour scale\n"
         "    -u                       make an unsteady grid problem steady\n\n");
    exit(1);
  }
  if (strcmp(&index[1], "post_file_manager")==0)
  {
    printf("usage: %s [options] control-file new-filename(optional)\n", progname);
    printf("\n"
         " options:\n"
         "    -s beg:end[:step]        read from beg to end-1 every step\n"
         "    -l [level]               set log level (0==none, infty==everything)\n"
         "    -o                       overwrite all existing files\n"
         "\n");
    exit(1);
  }

printf("usage: %s [options] control-file\n", progname);
printf("\n"
         "  options:\n"
         "    -s beg:end[:step]        read from beg to end-1 every step\n"
         "    -l [level]               set log level (0==none, infty==everything)\n\n");
exit(1);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the problem's data.

  \param problem       (o) problem object to be initialized
  \param argc          (i) number of command line arguments
  \param argv          (i) command line arguments

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_problem_data(PROBLEM_DATA* problem, INT argc, CHAR** argv)
{
  CHAR* problem_names[] = PROBLEMNAMES;
  CHAR* type;
  INT i, j, check;
  SYMBOL* symbol;
  CHAR* separator;
  MAP* control_table;
  CHAR* filename;
  INT lastarg=0;

#ifdef DEBUG
  dsinit();

  /* We need to take two steps back. dsinit is too close to ccarat. */
  dstrc_exit();
  dstrc_exit();

  dstrc_enter("init_problem_data");
#endif

  /*--------------------------------------------------------------------*/
  /* default values */

  problem->start = 0;
  problem->end = -1;
  problem->step = 1;

  /*--------------------------------------------------------------------*/
  /* process command line arguments */

  if (argc < 2)
  {
    usage(argv[0]);
  }

  if ((argc == 2) && (strcmp(argv[1], "-h")==0))
  {
    usage(argv[0]);
  }

  for (i=1; i<argc; ++i)
  {
    CHAR* arg;

    arg = argv[i];

    if (arg[0] == '-')
    {
      switch (arg[1])
      {
      case 's':                 /* slices */
      {
        if (arg[2] != '\0')
        {
          arg = &(arg[2]);
        }
        else
        {
          i += 1;
          if (i == argc-1)
          {
            /* dserror is not initialized yet */
            printf("%s: option '-s' must be followed by a slice like this: 'beg:end[:step]'\n", argv[0]);
            exit(1);
          }
          arg = argv[i];
        }

        /* simple parsing, only limited error checking */
        problem->start = atoi(arg);
        arg = strstr(arg, ":");
        if (arg == NULL)
        {
          /* dserror is not initialized yet */
          printf("%s: option '-s' must be followed by a slice like this: 'beg:end[:step]'\n", argv[0]);
          exit(1);
        }

        /* we support things like 'beg::step' and 'beg:' */
        if ((arg[1] != ':') && (arg[1] != '\0'))
        {
          problem->end = atoi(arg+1);
        }
        arg = strstr(arg+1, ":");
        if (arg != NULL)
        {
          problem->step = atoi(arg+1);
          if (problem->step < 1)
          {
            /* dserror is not initialized yet */
            printf("%s: step must be greater than zero\n", argv[0]);
            exit(1);
          }
        }
        break;
      }
      case 'h':
        usage(argv[0]);
        break;
      case 'l':
        if (arg[2] != '\0')
        {
          arg = &(arg[2]);
        }
        else
        {
          if (i+1 == argc-1)
          {
            /* no explizit number */
            loglevel = 3;
            break;
          }
          arg = argv[i+1];
          if (arg[0] >= '0' && arg[0] <= '9')
          {
            i += 1;
          }
          else {
            /* no explizit number */
            loglevel = 3;
            break;
          }
        }
        loglevel = atoi(arg);
        break;

        default:
          check=0;
          for (j=0;j<num_para;j++)
          {
            if (arg[1]==parameter[j]) check=1;
          }
          if (check==0)
          {
            printf("unsupported option '%s'", arg);
            usage(argv[0]);
          }
          else break;
      }
      if (i==argc-1) usage(argv[0]);
      if (argv[i+1][0]!='-') lastarg=i;

    }
  }

  /*--------------------------------------------------------------------*/
  /* setup fake ccarat environment and read control file */

  control_table = &(problem->control_table);

  setup_filter(argv[lastarg+1], control_table, problem->basename);

  dsassert(map_has_string(control_table, "version", "0.2"),
           "expect version 0.2 control file");

  /* Debug output */
  /*map_print(stdout, &control_table, 0);*/

  /*--------------------------------------------------------------------*/
  /* read general information */

  problem->ndim = map_read_int(control_table, "ndim");
  dsassert((problem->ndim == 2) || (problem->ndim == 3), "illegal dimension");

  type = map_read_string(control_table, "problem_type");
  for (i=0; problem_names[i] != NULL; ++i)
  {
    if (strcmp(type, problem_names[i])==0)
    {
      problem->type = i;
      break;
    }
  }
  if (problem_names[i] == NULL)
  {
    dserror("unknown problem type '%s'", type);
  }

  /*--------------------------------------------------------------------*/
  /* Find the input directory by looking at the control file
   * name. We look in the current directory for data files and if that
   * fails we look in the input directory. */

  /* This is the unix version. Different input directories are not
   * supported on windows. */
  filename = problem->basename;
  separator = rindex(filename, '/');
  if (separator == NULL)
  {
    problem->input_dir[0] = '\0';
  }
  else
  {
    INT n;

    /* 'separator-filename' gives the number of chars before the
     * separator. But we want to copy the separator as well. */
    n = separator-filename+1;
    dsassert(n < 100, "file name overflow");
    strncpy(problem->input_dir, filename, n);
    problem->input_dir[n] = '\0';
  }

  /*--------------------------------------------------------------------*/
  /* get the meaning of the elements' chunks */

  setup_element_variables_map(control_table);

  /*--------------------------------------------------------------------*/
  /* collect all result groups */

  problem->num_results = map_symbol_count(control_table, "result");
  if (problem->num_results == 0)
  {
    dserror("no results found");
  }
  problem->result_group = (MAP**)CCACALLOC(problem->num_results, sizeof(MAP*));

  /* find the first result group */
  symbol = map_find_symbol(control_table, "result");

  /* We rely on the fact that groups are linked in reverse order. */
  /* That is results are written ordered by time step. */
  for (i=problem->num_results-1; i>=0; --i)
  {
    if (!symbol_is_map(symbol))
    {
      dserror("failed to get result group");
    }

    problem->result_group[i] = symbol_map(symbol);

    symbol = symbol->next;
  }

  /*--------------------------------------------------------------------*/
  /* setup all fields */

  problem->num_discr = map_symbol_count(control_table, "field");
  if (problem->num_discr==0)
  {
    dserror("no field group found");
  }
  problem->discr = (FIELD_DATA*)CCACALLOC(problem->num_discr, sizeof(FIELD_DATA));

  problem->numele = 0;
  problem->numnp = 0;

  /* find the first field (the last one that has been written) */
  symbol = map_find_symbol(control_table, "field");

  /* read all fields headers, open the data files */
  for (i=0; i<problem->num_discr; ++i)
  {
    if (!symbol_is_map(symbol))
    {
      dserror("failed to get field group");
    }

    init_field_data(problem, &(problem->discr[i]), symbol_map(symbol));

    problem->numele += problem->discr[i].numele;
    problem->numnp  += problem->discr[i].numnp;

    symbol = symbol->next;
  }

  /*--------------------------------------------------------------------*/
  /* setup the translation tables */

  init_translation_table(&(problem->element_type),
                         map_read_map(control_table, "element_names"),
                         elementnames);

  init_translation_table(&(problem->distype),
                         map_read_map(control_table, "distype_names"),
                         distypenames);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether a given result group belongs to this field.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
INT match_field_result(FIELD_DATA *field, MAP *result_group)
{
  return (strcmp(map_read_string(result_group, "field"),
                 fieldnames[field->type]) == 0) &&
    (map_read_int(result_group, "field_pos") == field->field_pos) &&
    (map_read_int(result_group, "discretization") == field->disnum);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Initialize the result data.

  You need to call \a next_result to get to the first result of this
  discretization.

  \author u.kue
  \date 11/04
  \sa next_result
*/
/*----------------------------------------------------------------------*/
void init_result_data(FIELD_DATA* field, RESULT_DATA* result)
{
#ifdef DEBUG
  dstrc_enter("init_result_data");
#endif

  result->field = field;
  result->pos = -1;
  result->value_file = NULL;
  result->size_file = NULL;

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Cleanup result data.

  Please note that there must not be any chunk data on this result
  after this function has been called.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void destroy_result_data(RESULT_DATA* result)
{
#ifdef DEBUG
  dstrc_enter("destroy_result_data");
#endif

  if (result->value_file != NULL)
  {
    fclose(result->value_file);
    fclose(result->size_file);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Go to the next result of this discretization.

  \param result (i/o) on input the current result data

  \return zero if there is no further result. non-zero otherwise.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
INT next_result(RESULT_DATA* result)
{
  PROBLEM_DATA* problem;
  INT i;
  INT ret = 0;

#ifdef DEBUG
  dstrc_enter("next_result");
#endif

  problem = result->field->problem;

  for (i=result->pos+1; i<problem->num_results; ++i)
  {
    INT step;
    MAP* map;
    map = problem->result_group[i];

    if (match_field_result(result->field, map))
    {

      /*
       * Open the new files if there are any.
       *
       * If one of these files is here the other one has to be
       * here, too. If it's not, it's a bug in the input. */
      if ((map_symbol_count(map, "result_value_file") > 0) ||
          (map_symbol_count(map, "result_size_file") > 0))
      {
        if (result->value_file != NULL)
        {
          fclose(result->value_file);
          fclose(result->size_file);
        }
        open_data_files(result, map, "result");
      }

      /*
       * We use the real step numbers here. That is a user has to give
       * the real numbers, too. Maybe that's the best way to handle
       * it. */
      /* In case of FSI everything else hurts even more. */
      step = map_read_int(map, "step");

      /* we are only interessted if the result matches the slice */
      if ((step >= problem->start) &&
          ((step <= problem->end) || (problem->end == -1)) &&
          ((step - problem->start) % problem->step == 0))
      {
        result->pos = i;
        result->group = map;
        ret = 1;
        break;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*======================================================================*/
/*======================================================================*/


  /* on little endian machines we have to convert */
  /* We have 8 byte doubles and 4 byte integer by definition. Nothing
   * else. */

#ifdef IS_LITTLE_ENDIAN

  /* very specific swap macro */
#define SWAP_CHAR(c1,c2) { CHAR t; t=c1; c1=c2; c2=t; }


/*----------------------------------------------------------------------*/
/*!
  \brief Convert the read big endian values to little endian.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
static void byteswap_doubles(DOUBLE* data, INT length)
{
  INT i;

  for (i=0; i<length; ++i)
  {
    CHAR* ptr;

    ptr = (CHAR*)&(data[i]);
    SWAP_CHAR(ptr[0], ptr[7]);
    SWAP_CHAR(ptr[1], ptr[6]);
    SWAP_CHAR(ptr[2], ptr[5]);
    SWAP_CHAR(ptr[3], ptr[4]);
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Convert the read big endian values to little endian.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
static void byteswap_ints(INT* data, INT length)
{
  INT i;
  for (i=0; i<length; ++i)
  {
    CHAR* ptr;

    ptr = (CHAR*)&(data[i]);
    SWAP_CHAR(ptr[0], ptr[3]);
    SWAP_CHAR(ptr[1], ptr[2]);
  }
}

#undef SWAP_CHAR

#else

/* noops to make life easier further down */
#define byteswap_doubles(data, length)
#define byteswap_ints(data, length)

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief Set up the chunk structure to iterate the chunk's entries.

  \param result (i) on input the current result data
  \param chunk  (o) the chunk to be initialized
  \param name   (i) the name of the group in the control file

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void init_chunk_data(RESULT_DATA* result, CHUNK_DATA* chunk, CHAR* name)
{
#ifndef LOWMEM
  CHAR* type;
#endif

#ifdef DEBUG
  dstrc_enter("init_chunk_data");
#endif

  chunk->result = result;
  chunk->group = map_read_map(result->group, name);

  chunk->value_entry_length = map_read_int(chunk->group, "value_entry_length");
  chunk->value_offset = map_read_int(chunk->group, "value_offset");

  chunk->size_entry_length = map_read_int(chunk->group, "size_entry_length");
  chunk->size_offset = map_read_int(chunk->group, "size_offset");

#ifdef LOWMEM

  /* Low memory! We read one entry only. This way we have to reread
   * some entries many times. */

  if (chunk->value_entry_length > 0)
  {
    chunk->value_buf = (DOUBLE*)CCACALLOC(chunk->value_entry_length, sizeof(DOUBLE));
  }
  else
  {
    chunk->value_buf = NULL;
  }

  if (chunk->size_entry_length > 0)
  {
    chunk->size_buf = (INT*)CCACALLOC(chunk->size_entry_length, sizeof(INT));
  }
  else
  {
    chunk->size_buf = NULL;
  }

#else

  /* More memory (smaller problem size). We read the whole chunk at
   * once. This is supposed to be fast. */

  type = map_read_string(chunk->group, "type");

  if (strcmp(type, "element") == 0)
  {
    INT length = chunk->value_entry_length*result->field->numele;

    if (length > 0)
    {
      chunk->value_data = (DOUBLE*)CCACALLOC(length, sizeof(DOUBLE));
      fseek(chunk->result->value_file, chunk->value_offset, SEEK_SET);
      if (fread(chunk->value_data, sizeof(DOUBLE),
                length,
                chunk->result->value_file) != length)
      {
        dserror("failed to read value file of field '%s'", result->field->name);
      }
      byteswap_doubles(chunk->value_data, length);
    }
    else
    {
      chunk->value_data = NULL;
    }

    length = chunk->size_entry_length*result->field->numele;

    if (length > 0)
    {
      chunk->size_data = (INT*)CCACALLOC(length, sizeof(INT));
      fseek(chunk->result->size_file, chunk->size_offset, SEEK_SET);
      if (fread(chunk->size_data, sizeof(INT),
                length,
                chunk->result->size_file) != length)
      {
        dserror("failed to read size file of field '%s'", result->field->name);
      }
      byteswap_ints(chunk->size_data, length);
    }
    else
    {
      chunk->size_data = NULL;
    }
  }
  else if (strcmp(type, "node") == 0)
  {
    INT length = chunk->value_entry_length*result->field->numnp;

    if (length > 0)
    {
      INT length_read;
      chunk->value_data = (DOUBLE*)CCACALLOC(length, sizeof(DOUBLE));
      fseek(chunk->result->value_file, chunk->value_offset, SEEK_SET);
      length_read = fread(chunk->value_data, sizeof(DOUBLE),
                          length,
                          chunk->result->value_file);
      if (length_read != length)
      {
        dserror("failed to read value file of field '%s'", result->field->name);
      }
      byteswap_doubles(chunk->value_data, length);
    }
    else
    {
      chunk->value_data = NULL;
    }

    length = chunk->size_entry_length*result->field->numnp;

    if (length > 0)
    {
      INT length_read;
      chunk->size_data = (INT*)CCACALLOC(length, sizeof(INT));
      fseek(chunk->result->size_file, chunk->size_offset, SEEK_SET);
      length_read = fread(chunk->size_data, sizeof(INT),
                          length,
                          chunk->result->size_file);
      if (length_read != length)
      {
        dserror("failed to read size file of field '%s'", result->field->name);
      }
      byteswap_ints(chunk->size_data, length);
    }
    else
    {
      chunk->size_data = NULL;
    }
  }
  else
  {
    dserror("chunk type '%s' not supported", type);
  }

#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Free the chunk data.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void destroy_chunk_data(CHUNK_DATA* chunk)
{
#ifdef DEBUG
  dstrc_enter("destroy_chunk_data");
#endif

#ifdef LOWMEM

  if (chunk->value_buf != NULL)
    CCAFREE(chunk->value_buf);

  if (chunk->size_buf != NULL)
    CCAFREE(chunk->size_buf);

#else

  if (chunk->value_data != NULL)
    CCAFREE(chunk->value_data);

  if (chunk->size_data != NULL)
    CCAFREE(chunk->size_data);

#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read one size entry form the file and store it to this
  chunk_data's internal buffer.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void chunk_read_size_entry(CHUNK_DATA* chunk, INT id)
{
#ifdef DEBUG
  dstrc_enter("chunk_read_size_entry");
#endif

  dsassert(chunk->size_entry_length > 0, "cannot read empty entry");

#ifdef LOWMEM

  fseek(chunk->result->size_file,
        chunk->size_offset + chunk->size_entry_length*id*sizeof(INT),
        SEEK_SET);
  if (fread(chunk->size_buf, sizeof(INT),
            chunk->size_entry_length,
            chunk->result->size_file) != chunk->size_entry_length)
  {
    dserror("failed to read size file of field '%s'", chunk->result->field->name);
  }
  byteswap_ints(chunk->size_buf, chunk->size_entry_length);

#else

  chunk->size_buf = &(chunk->size_data[chunk->size_entry_length*id]);

#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read one value entry form the file and store it to this
  chunk_data's internal buffer.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void chunk_read_value_entry(CHUNK_DATA* chunk, INT id)
{
#ifdef DEBUG
  dstrc_enter("chunk_read_value_entry");
#endif

  dsassert(chunk->value_entry_length > 0, "cannot read empty entry");

#ifdef LOWMEM

  fseek(chunk->result->value_file,
        chunk->value_offset + chunk->value_entry_length*id*sizeof(DOUBLE),
        SEEK_SET);
  if (fread(chunk->value_buf, sizeof(DOUBLE),
            chunk->value_entry_length,
            chunk->result->value_file) != chunk->value_entry_length)
  {
    dserror("failed to read value file of field '%s'", chunk->result->field->name);
  }
  byteswap_doubles(chunk->value_buf, chunk->value_entry_length);

#else

  chunk->value_buf = &(chunk->value_data[chunk->value_entry_length*id]);

#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read the element parameters common to all elements.

  This is just a service to make life easier and yet to do extensive
  error checking.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void get_element_params(FIELD_DATA* field,
                        INT i,
                        INT* Id, INT* el_type, INT* dis, INT* numnp)
{
#ifdef DEBUG
  dstrc_enter("get_element_params");
#endif

  chunk_read_size_entry(&(field->ele_param), i);

  *Id      = field->ele_param.size_buf[element_variables.ep_size_Id];

  *el_type = field->ele_param.size_buf[element_variables.ep_size_eltyp];
  if ((*el_type < 0) || (*el_type >= field->problem->element_type.length))
  {
    dserror("element type %d exceeds range", *el_type);
  }

  /* translate to internal value */
  *el_type = field->problem->element_type.table[*el_type];

  *dis     = field->ele_param.size_buf[element_variables.ep_size_distyp];
  if ((*dis < 0) || (*dis >= field->problem->distype.length))
  {
    dserror("element dis %d exceeds range", *dis);
  }

  *numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*!
  \brief Set up a (fake) discretization.

  Create the node and element arrays, read node coordinates and mesh
  connectivity.

  \param discret       (o) Uninitialized discretization object
  \param problem       (i) problem data
  \param field         (i) general field data
  \param redef_hex20   (i) make hex27 out of hex20 elements

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void init_post_discretization(POST_DISCRETIZATION* discret,
                              PROBLEM_DATA* problem,
                              FIELD_DATA* field,
                              INT redef_hex20)
{
  INT i;
  INT j;
  INT *offset;
  INT numnp;
  INT numnp_old;
  INT* nodelist=NULL;

#ifdef DEBUG
  dstrc_enter("init_post_discretization");
#endif

  discret->field = field;

  discret->node = (NODE*)CCACALLOC(field->numnp, sizeof(NODE));
  discret->element = (ELEMENT*)CCACALLOC(field->numele, sizeof(ELEMENT));
  numnp_old=field->numnp;

  /*--------------------------------------------------------------------*/
  /* read the node coordinates */
  for (i=0; i<field->numnp; ++i)
  {
    discret->node[i].Id_loc = i;
    discret->node[i].proc = 0;
    discret->node[i].numele = 0;

    chunk_read_size_entry(&(field->coords), i);
    discret->node[i].Id = field->coords.size_buf[element_variables.ep_size_Id];

    chunk_read_value_entry(&(field->coords), i);
    for (j=0; j<field->coords.value_entry_length; ++j)
      discret->node[i].x[j] = field->coords.value_buf[j];
  }

  /*--------------------------------------------------------------------*/
  /* read the mesh */
  chunk_read_size_entry(&(field->ele_param), 0);
  first=1;
  for (i=0; i<field->numele; ++i)
  {
    INT Id;
    INT el_type;
    INT numnp;
    INT distype;
    INT j;

    chunk_read_size_entry(&(field->ele_param), i);

    Id      = field->ele_param.size_buf[element_variables.ep_size_Id];
    el_type = field->ele_param.size_buf[element_variables.ep_size_eltyp];
    distype = field->ele_param.size_buf[element_variables.ep_size_distyp];
    numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

    /* external >= internal */
    el_type = field->problem->element_type.table[el_type];
    distype = field->problem->distype.table[distype];

    chunk_read_size_entry(&(field->mesh), i);

    discret->element[i].Id = Id;
    discret->element[i].Id_loc = i;
    discret->element[i].proc = 0;

    if (redef_hex20 && distype==hex20)
    {
      DOUBLE x[3];
      NODE* tempnode[20];
      int k;
      /* hex20-elements are considered as hex27-elements
       * missing nodes will be generated
       * coordinates of the new nodes are interpolated immediately,
       * solutions in function lin_interopol */
      if (i==0)
      {
        discret->node = (NODE*)CCAREALLOC(discret->node, 2*(field->numnp*sizeof(NODE)));
        nodelist=(INT*)CCACALLOC(field->numnp, sizeof(INT));
        for (k=0;k<field->numnp;k++)
        {
          nodelist[k]=0;
        }
      }

      discret->element[i].numnp = 27;
      discret->element[i].node = (NODE**)CCACALLOC(27, sizeof(NODE*));
      discret->element[i].eltyp = el_type;
      discret->element[i].distyp = distype;

      for (j=0; j<20; ++j)
      {
        discret->element[i].node[j] = &(discret->node[field->mesh.size_buf[j]]);
        discret->element[i].node[j]->numele++;
      }
/* The node map of hierarchical hex20-elements have to be
 * redefined, so that the nodes are on the same place like nodes with the
 * same number in a hex27-element */

/* node map of hex27 elements                       */
/* this are the node maps used in GiD               */

/*  z y           7              18             6   */
/*  |/            o--------------o--------------o   */
/*  o-x          /:             /              /|   */
/*              / :            /              / |   */
/*             /  :           /              /  |   */
/*          19/   :        25/           17 /   |   */
/*           o--------------o--------------o    |   */
/*          /     :        /              /|    |   */
/*         /    15o       /    23o       / |  14o   */
/*        /       :      /              /  |   /|   */
/*      4/        :  16 /             5/   |  / |   */
/*      o--------------o--------------o    | /  |   */
/*      |         :    |   26         |    |/   |   */
/*      |  24o    :    |    o         |  22o    |   */
/*      |         :    |       10     |   /|    |   */
/*      |        3o....|.........o....|../.|....o   */
/*      |        .     |              | /  |   / 2  */
/*      |       .    21|            13|/   |  /     */
/*   12 o--------------o--------------o    | /      */
/*      |     .        |              |    |/       */
/*      |  11o         | 20o          |    o        */
/*      |   .          |              |   / 9       */
/*      |  .           |              |  /          */
/*      | .            |              | /           */
/*      |.             |              |/            */
/*      o--------------o--------------o             */
/*      0              8              1             */

/* node map of hierarchical hex20 elements          */
/*                3              11             0   */
/*  z y           o--------------o--------------o   */
/*  |/           /:                            /|   */
/*  o-x         / :                           / |   */
/*             /  :                          /  |   */
/*          10/   :                       8 /   |   */
/*           o    :                        o    |   */
/*          /     :                       /     |   */
/*         /    19o                      /    16o   */
/*        /       :                     /       |   */
/*      2/        :    9              1/        |   */
/*      o--------------o--------------o         |   */
/*      |         :                   |         |   */
/*      |         :                   |         |   */
/*      |         :            15     |         |   */
/*      |        7o..............o....|.........o   */
/*      |        .                    |        / 4  */
/*      |       .                   17|       /     */
/*   18 o      .                      o      /      */
/*      |     .                       |     /       */
/*      |  14o                        |    o        */
/*      |   .                         |   / 12      */
/*      |  .                          |  /          */
/*      | .                           | /           */
/*      |.                            |/            */
/*      o--------------o--------------o             */
/*      6             13              5             */

      for (j=12; j<20; j++)
        tempnode[j] = discret->element[i].node[j];
/*       if (distype == h_hex20) */
/*       { discret->element[i].node[12]=tempnode[16]; */
/*         discret->element[i].node[13]=tempnode[17]; */
/*         discret->element[i].node[14]=tempnode[18]; */
/*         discret->element[i].node[15]=tempnode[19]; */
/*         discret->element[i].node[16]=tempnode[12]; */
/*         discret->element[i].node[17]=tempnode[13]; */
/*         discret->element[i].node[18]=tempnode[14]; */
/*         discret->element[i].node[19]=tempnode[15]; */
/*       } */

/* node map of hex20-elements                          */
/*                 7              18             6     */
/*                o--------------o--------------o      */
/*  z y          /:                            /|      */
/*  |/          / :                           / |      */
/*  o-x        /  :                          /  |      */
/*          19/   :                      17 /   |      */
/*           o    :                        o    |      */
/*          /     :                       /     |      */
/*         /    15o                      /    14o      */
/*        /       :                     /       |      */
/*      4/        :                   5/        |      */
/*      o--------------o--------------o         |      */
/*      |         :     16            |         |      */
/*      |         :                   |         |      */
/*      |         :            10     |         |      */
/*      |        3o..............o....|.........o      */
/*      |        .                    |        / 2     */
/*      |       .                   13|       /        */
/*   12 o      .                      o      /         */
/*      |     .                       |     /          */
/*      |  11o                        |    o           */
/*      |   .                         |   / 9          */
/*      |  .                          |  /             */
/*      | .                           | /              */
/*      |.                            |/               */
/*      o--------------o--------------o                */
/*      0              8              1                */

      /* linear interpolation of coordinates */
      for (k=0; k<3; k++)
        x[k]=linear_interpolation( 0,  0, -1, &discret->element[i], k);
      field->numnp=test_node(field->numnp,numnp_old, discret, i, 20, x, nodelist);

      for (k=0; k<3; k++)
        x[k]=linear_interpolation( 0, -1,  0, &discret->element[i], k);
      field->numnp=test_node(field->numnp,numnp_old, discret, i, 21, x, nodelist);

      for (k=0; k<3; k++)
        x[k]=linear_interpolation( 1,  0,  0, &discret->element[i], k);
      field->numnp=test_node(field->numnp,numnp_old, discret, i, 22, x, nodelist);

      for (k=0; k<3; k++)
        x[k]=linear_interpolation( 0,  1,  0, &discret->element[i], k);
      field->numnp=test_node(field->numnp, numnp_old,discret, i, 23, x, nodelist);

      for (k=0; k<3; k++)
        x[k]=linear_interpolation(-1,  0,  0, &discret->element[i], k);
      field->numnp=test_node(field->numnp, numnp_old,discret, i, 24, x, nodelist);

      for (k=0; k<3; k++)
        x[k]=linear_interpolation( 0,  0,  1, &discret->element[i], k);
      field->numnp=test_node(field->numnp, numnp_old,discret, i, 25, x, nodelist);

      for (k=0; k<3; k++)
        x[k]=linear_interpolation( 0,  0,  0, &discret->element[i], k);
      field->numnp=test_node(field->numnp, numnp_old,discret, i, 26, x, nodelist);

      if (i==field->numele-1)
        discret->node=(NODE*)CCAREALLOC(discret->node, (field->numnp)*sizeof(NODE));
    }
    else
    {
      discret->element[i].node = (NODE**)CCACALLOC(numnp, sizeof(NODE*));
      discret->element[i].numnp = numnp;
      discret->element[i].eltyp = el_type;
      discret->element[i].distyp = distype;
      for (j=0; j<numnp; ++j)
      {
        discret->element[i].node[j] = &(discret->node[field->mesh.size_buf[j]]);
        discret->element[i].node[j]->numele++;
      }
    }
  }

  if (nodelist!=NULL)
  {
    CCAFREE(nodelist);
  }

  offset=(INT*)CCACALLOC(field->numnp, sizeof(INT));

  for (i=0;i<field->numnp; ++i)
  {
    discret->node[i].element=(ELEMENT**)CCACALLOC(discret->node[i].numele, sizeof(ELEMENT*));
    offset[i]=0;
  }
  for (i=0;i<field->numele;++i)
  {
    numnp = discret->element[i].numnp;
    for (j=0;j<numnp;++j)
    {
      discret->element[i].node[j]->element[offset[discret->element[i].node[j]->Id_loc]]=&discret->element[i];
      offset[discret->element[i].node[j]->Id_loc]++;
    }
  }
  CCAFREE(offset);

#ifdef DEBUG
  dstrc_exit();
#endif
}


void destroy_post_discretization(POST_DISCRETIZATION* discret)
{
  INT i;

#ifdef DEBUG
  dstrc_enter("destroy_post_discretization");
#endif

  for (i=0; i<discret->field->numnp; ++i)
  {
    CCAFREE(discret->node[i].element);
  }

  for (i=0; i<discret->field->numele; ++i)
  {
    CCAFREE(discret->element[i].node);
  }

  CCAFREE(discret->node);
  CCAFREE(discret->element);

#ifdef DEBUG
  dstrc_exit();
#endif
}

void init_post_design(POST_DESIGN* design, FIELD_DATA* field)
{
  INT ndvol;
  INT ndsurf;
  INT ndline;
  INT ndnode;
  INT i;


#ifdef DEBUG
  dstrc_enter("init_post_design");
#endif

  design->field = field;
  design->nnode=field->numnp;

  map_find_int(field->group, "ndvol", &ndvol);
  design->volume=(VOLUME*)CCACALLOC(ndvol, sizeof(VOLUME));
  design->ndvol=ndvol;

  map_find_int(field->group, "ndsurf", &ndsurf);
  design->surface=(SURFACE*)CCACALLOC(ndsurf, sizeof(SURFACE));
  design->ndsurf=ndsurf;

  map_find_int(field->group, "ndline", &ndline);
  design->line=(LINE*)CCACALLOC(ndline, sizeof(LINE));
  design->ndline=ndline;

  map_find_int(field->group, "ndnode", &ndnode);
  design->dnode=(DESIGNNODE*)CCACALLOC(ndline, sizeof(DESIGNNODE));
  design->ndnode=ndnode;

  /*initialisation of design variables */
  for (i=0;i<design->ndvol;i++)
  {
    design->volume[i].id=i;
    design->volume[i].ndsurf=0;
    design->volume[i].nnode=0;
    design->volume[i].surface=(SURFACE**)CCACALLOC(design->ndsurf, sizeof(SURFACE*));
    design->volume[i].node=(NODE**)CCACALLOC(design->nnode, sizeof(NODE*));
  }

  for (i=0;i<design->ndsurf;i++)
  {
    design->surface[i].id=i;
    design->surface[i].ndline=0;
    design->surface[i].nnode=0;
    design->surface[i].line=(LINE**)CCACALLOC(design->ndline, sizeof(LINE*));
    design->surface[i].node=(NODE**)CCACALLOC(design->nnode, sizeof(NODE*));
  }

  for (i=0;i<design->ndline;i++)
  {
    design->line[i].id=i;
    design->line[i].nnode=0;
    design->line[i].ndnode=0;
    design->line[i].dnode=(DESIGNNODE**)CCACALLOC(design->ndnode, sizeof(DESIGNNODE*));
    design->line[i].node=(NODE**)CCACALLOC(design->nnode, sizeof(NODE*));
  }

  for (i=0;i<design->ndnode;i++)
  {
    design->dnode[i].id=i;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


#ifdef D_FSI

/*----------------------------------------------------------------------*/
/*!
  \brief Find the connection between ale and fluid nodes.

  \param problem       (i) The problem data
  \param struct_field  (i) struct field data; might be NULL
  \param fluid_field   (i) fluid field data
  \param ale_field     (i) ale field data
  \param _fluid_struct_connect  (o) per fluid node array that gives
                                    the corresponding structure nodes
  \param _fluid_ale_connect     (o) gives the corresponding ale nodes

  \warning The connection array are allocated here but must be freed
  by the caller.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void post_find_fsi_coupling(PROBLEM_DATA *problem,
                            FIELD_DATA *struct_field,
                            FIELD_DATA *fluid_field,
                            FIELD_DATA *ale_field,
                            INT **_fluid_struct_connect,
                            INT **_fluid_ale_connect)
{
  INT i;

  INT *fluid_struct_connect = NULL;
  INT *fluid_ale_connect;

#ifdef DEBUG
  dstrc_enter("post_find_fsi_coupling");
#endif

  /* read the fluid node coordinates one by one and search for
   * matching struct and ale nodes. */

  /* We have to allocate flags for all fluid nodes. */

  fluid_ale_connect = (INT*)CCACALLOC(fluid_field->numnp, sizeof(INT));
  if (struct_field != NULL) {
    fluid_struct_connect = (INT*)CCACALLOC(fluid_field->numnp, sizeof(INT));
  }

  /* This is a quadratic loop. If it turns out to be too slow one
   * could implement some quad- or octtree algorithm. */
  for (i=0; i<fluid_field->numnp; ++i) {
    INT n_ale;

    chunk_read_value_entry(&(fluid_field->coords), i);

    /* search the structure nodes */
    if (struct_field != NULL) {
      INT n_struct;

      /* no corresponding struct node by default */
      fluid_struct_connect[i] = -1;

      /* search the struct node */
      for (n_struct=0; n_struct<struct_field->numnp; ++n_struct) {
        INT k;
        DOUBLE diff = 0;

        chunk_read_value_entry(&(struct_field->coords), n_struct);

        /* quadratic error norm */
        for (k=0; k<problem->ndim; ++k) {
          DOUBLE d = fluid_field->coords.value_buf[k] - struct_field->coords.value_buf[k];
          diff += d*d;
        }

        /*
         * If the difference is very small we've found the corresponding
         * node. The tolerance here might be too big for very fine
         * meshes. I don't know. (ccarat uses the same tolerance but
         * applies it to the sqare root. I don't want to waste the time
         * to calculate it...)
         *
         * In fluid_struct_connect we store the local indices, that is
         * no real ids. */
        if (diff < 1e-10) {
          fluid_struct_connect[i] = n_struct;
          break;
        }
      }
    }

    /* search the ale nodes */

    /* no corresponding ale node by default */
    fluid_ale_connect[i] = -1;

    /* search the ale node */
    for (n_ale=0; n_ale<ale_field->numnp; ++n_ale) {
      INT k;
      DOUBLE diff = 0;

      chunk_read_value_entry(&(ale_field->coords), n_ale);

      /* quadratic error norm */
      for (k=0; k<problem->ndim; ++k) {
        DOUBLE d = fluid_field->coords.value_buf[k] - ale_field->coords.value_buf[k];
        diff += d*d;
      }

      /*
       * If the difference is very small we've found the corresponding
       * node. The tolerance here might be too big for very fine
       * meshes. I don't know. (ccarat uses the same tolerance but
       * applies it to the sqare root. I don't want to waste the time
       * to calculate it...)
       *
       * In fluid_ale_connect we store the local indices, that is no
       * real ids. */
      if (diff < 1e-10) {
        fluid_ale_connect[i] = n_ale;
        break;
      }
    }
  }

  *_fluid_struct_connect = fluid_struct_connect;
  *_fluid_ale_connect = fluid_ale_connect;

#ifdef DEBUG
  dstrc_exit();
#endif
}
#endif

