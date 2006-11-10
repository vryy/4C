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


int main(int argc, char** argv)
{
  PROBLEM_DATA problem1;
  PROBLEM_DATA problem2;
  char** singleargs;
  INT i;

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
    dserror("problem types differ");
  
  if (problem1.ndim != problem2.ndim)
    dserror("dimensions differ");
  
  if (problem1.num_discr != problem2.num_discr)
    dserror("number of discretizations differ");

  for (i=0; i<problem1.num_discr; ++i)
  {
    INT resnum = 0;
    FIELD_DATA* field1;
    FIELD_DATA* field2;
    RESULT_DATA result1;
    RESULT_DATA result2;
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

      switch (problem1.type)
      {
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
        
        for (j=0; j<field1->numnp; ++j)
        {
          chunk_read_value_entry(&chunk1, j);
          chunk_read_value_entry(&chunk2, j);

          for (k=0; k<chunk1.value_entry_length; ++k)
          {
            DOUBLE diff = chunk1.value_buf[k] - chunk2.value_buf[k];
            L1[k] += FABS(diff);
            L2[k] += diff*diff;
            Linf[k] = MAX(Linf[k], FABS(diff));
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
  }
  
  return 0;
}
