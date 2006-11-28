/*!
\file
\brief compare two binary result files
<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

\author u.kue
\date 11/06

*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "../post_common/post_common.h"
#include "post_cmp.h"
#include "../pss_full/pss_set.h"

#ifdef IS_LITTLE_ENDIAN

  /* very specific swap macro */
#define SWAP_CHAR(c1,c2) { CHAR t; t=c1; c1=c2; c2=t; }

#else

#define SWAP_CHAR(c1,c2)

#endif


int main(int argc, char** argv)
{
  PROBLEM_DATA problem1;
  PROBLEM_DATA problem2;
  char** singleargs;
  INT i;
  FILE* cf;
  CHAR buf[1024];
  MAP_ITERATOR iterator;

  /* Here we want to read two control files. So reconstruct the arg
   * vector. */

  singleargs = (char**)CCAMALLOC((argc-1)*sizeof(char*));
  for (i=0; i<argc-2; ++i)
    singleargs[i] = argv[i];

  singleargs[argc-2] = argv[argc-2];
  init_problem_data(&problem1, argc-1, singleargs);

  singleargs[argc-2] = argv[argc-1];
  init_problem_data(&problem2, argc-1, singleargs);

  if (problem1.type != problem2.type)
  {
    if (((problem1.type == prb_fluid   ) && (problem2.type == prb_fluid_pm)) ||
        ((problem1.type == prb_fluid_pm) && (problem2.type == prb_fluid   )))
    {
      /* fine */
    }
    else
      dserror("problem types differ");
  }
  
  if (problem1.ndim != problem2.ndim)
    dserror("dimensions differ");
  
  if (problem1.num_discr != problem2.num_discr)
    dserror("number of discretizations differ");

  /* setup difference control file */
  sprintf(buf, "%s-%s.control", argv[argc-2], argv[argc-1]);
  cf = fopen(buf, "w");

  /* Write the preamble and the mesh from the first input
   * problem. This assumes that this problem is in the current
   * directory. This way we avoid to duplicate the mesh files. */

  init_map_iterator(&iterator, &problem1.control_table);
  while (next_map_node(&iterator))
  {
    MAP_NODE* node;
    node = iterator_get_node(&iterator);
    if (strcmp(node->key,"result")!=0 &&
        strcmp(node->key,"restart")!=0)
    {
      printf("write '%s'\n", node->key);
      symbol_print(cf, node->key, node->symbol, 0);
    }
  }
    
  for (i=0; i<problem1.num_discr; ++i)
  {
    INT resnum = 0;
    FIELD_DATA* field1;
    FIELD_DATA* field2;
    RESULT_DATA result1;
    RESULT_DATA result2;

    static FILE* vf[10];
    static FILE* sf[10];
    
    field1 = &(problem1.discr[i]);
    field2 = &(problem2.discr[i]);

    if ((field1->numele != field2->numele) ||
        (field1->numnp != field2->numnp))
      dserror("discretizations do not match");

    init_result_data(field1, &result1);
    init_result_data(field2, &result2);
    while (next_result(&result1) && next_result(&result2))
    {
      DOUBLE time1;
      DOUBLE time2;

      time1 = map_read_real(result1.group, "time");
      time2 = map_read_real(result2.group, "time");
      if (fabs(time1-time2) > 1e-6)
        dserror("time difference: t1=%f, t2=%f", time1, time2);

      if (i>=10)
        dserror("too many discetizations");

      fprintf(cf, "result:\n"
              "    field = \"%s\"\n"
              "    field_pos = %d\n"
              "    discretization = %d\n"
              "    time = %f\n"
              "    step = %d\n"
              "\n", 
              map_read_string(result1.group, "field"),
              map_read_int(result1.group, "field_pos"),
              map_read_int(result1.group, "discretization"),
              map_read_real(result1.group, "time"), 
              map_read_int(result1.group, "step"));
        
      if (vf[i]==NULL)
      {
        sprintf(buf, "%s-%s.result.%s.f%d.d%d.s0.values",
                argv[argc-2], argv[argc-1],
                map_read_string(result1.group, "field"),
                map_read_int(result1.group, "field_pos"),
                map_read_int(result1.group, "discretization"));
        vf[i] = fopen(buf, "w");
        fprintf(cf, "    result_value_file = \"%s\"\n", buf);
          
        sprintf(buf, "%s-%s.result.%s.f%d.d%d.s0.sizes",
                argv[argc-2], argv[argc-1],
                map_read_string(result1.group, "field"),
                map_read_int(result1.group, "field_pos"),
                map_read_int(result1.group, "discretization"));
        sf[i] = fopen(buf, "w");
        fprintf(cf, "    result_size_file = \"%s\"\n", buf);
        fprintf(cf, "\n");
      }

      switch (problem1.type)
      {
      case prb_fluid:
      case prb_fluid_pm:
      {
        DOUBLE L1[3] = {0, 0, 0};
        DOUBLE L2[3] = {0, 0, 0};
        DOUBLE Linf[3] = {0, 0, 0};
        INT j;
        INT k;
        CHUNK_DATA chunk1;
        CHUNK_DATA chunk2;

        init_chunk_data(&result1, &chunk1, "velocity");
        init_chunk_data(&result2, &chunk2, "velocity");

        if ((chunk1.value_entry_length != chunk2.value_entry_length) ||
            (chunk1.value_entry_length > 3))
          dserror("illegal value entry length");

        fprintf(cf, "    velocity:\n"
                "        size_entry_length = 0\n"
                "        size_offset = 0\n"
                "        type = \"node\"\n"
                "        value_entry_length = %d\n"
                "        value_offset = %ld\n"
                "\n",
                chunk1.value_entry_length, ftell(vf[i]));
        
        for (j=0; j<field1->numnp; ++j)
        {
          chunk_read_value_entry(&chunk1, j);
          chunk_read_value_entry(&chunk2, j);

          for (k=0; k<chunk1.value_entry_length; ++k)
          {
            CHAR* ptr;
            DOUBLE diff = chunk1.value_buf[k] - chunk2.value_buf[k];
            L1[k] += FABS(diff);
            L2[k] += diff*diff;
            Linf[k] = MAX(Linf[k], FABS(diff));

            ptr = (CHAR*)&diff;
            SWAP_CHAR(ptr[0], ptr[7]);
            SWAP_CHAR(ptr[1], ptr[6]);
            SWAP_CHAR(ptr[2], ptr[5]);
            SWAP_CHAR(ptr[3], ptr[4]);
            fwrite(ptr, sizeof(DOUBLE), 1, vf[i]);
          }
        }

        if (resnum==0)
          printf("\t\t|L1|\t\t|L2|\t\t|Linf|\n");
        printf("%f:\n", time1);
        for (k=0; k<chunk1.value_entry_length; ++k)
        {
          char dir;
          dir = 'x'+k;
          L2[k] = sqrt(L2[k]);
          printf("\t%c:\t%f\t%f\t%f\n", dir, L1[k], L2[k], Linf[k]);
        }

        destroy_chunk_data(&chunk1);
        destroy_chunk_data(&chunk2);
        break;
      }
      default:
      {
        char *names[] = PROBLEMNAMES;
        dserror("problem type '%s' not supported", names[problem1.type]);
      }
      }
      
      resnum += 1;
    }

    destroy_result_data(&result1);
    destroy_result_data(&result2);
  }

  fclose(cf);
  
  return 0;
}
