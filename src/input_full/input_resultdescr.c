/*!
\file
\brief Read the predefined (expected) results from the input file.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

\author uk
\date 06/04

*/
#ifdef RESULTTEST

#include "../headers/standardtypes.h"

extern struct _FILES allfiles;
extern struct _GENPROB genprob;
extern struct _RESULTDESCR* resultdescr;
extern struct _PAR par;


/*!
 * \brief Read all the result description lines.
 *
 * Read all the result description lines. Such a line looks like:
 *
 * <pre>
 * fieldname DIS x NODE y POSITION sol(row,col) NAME name VALUE value TOLERANCE tolerance
 * </pre>
 *
 * or
 *
 * <pre>
 * fieldname DIS x ELEMENT y POSITION somename(i,...,k) NAME name VALUE value TOLERANCE tolerance
 * </pre>
 *
 * or
 *
 * <pre>
 * fieldname SPECIAL x
 * </pre>
 */
void inp_resultdescr()
{
  INT ierr;
  INT i;

#ifdef DEBUG
  dstrc_enter("inp_resultdescr");
#endif

  genprob.numresults = 0;

  if (frfind("--RESULT DESCRIPTION")==1) {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0) {
      genprob.numresults++;
      frread();
    }

    resultdescr = (RESULTDESCR*)CCACALLOC(genprob.numresults, sizeof(RESULTDESCR));

    if (frfind("--RESULT DESCRIPTION")==1) {
      frread();
      for (i = 0; strncmp(allfiles.actplace,"------",6)!=0; ++i) {

        frchk("STRUCTURE",&ierr);
        if (ierr==1) {
          resultdescr[i].field = structure;
        }
        else {
          frchk("FLUID",&ierr);
          if (ierr==1) {
            resultdescr[i].field = fluid;
          }
          else {
            frchk("ALE",&ierr);
            if (ierr==1) {
              resultdescr[i].field = ale;
            }
            else {
              frchk("THERMAL",&ierr);
              if (ierr==1) {
                resultdescr[i].field = thermal;
              }
              else {
                dserror("Unknown field type");
              }
            }
          }
        }

        frint("DIS", &(resultdescr[i].dis), &ierr);
        if (ierr==1) {
          resultdescr[i].dis--;
        }
        else {
          /* abuse the dis field for an extra parameter */
          frint("SPECIAL", &(resultdescr[i].dis), &ierr);
          if (ierr==1) {
            resultdescr[i].node = -1;
            resultdescr[i].element = -1;
            goto end;
          }
          dserror("Failed to read discretisation number");
        }

        frint("NODE", &(resultdescr[i].node), &ierr);
        if (ierr==1) {
          resultdescr[i].node--;
          resultdescr[i].element = -1;
        }
        else {
          resultdescr[i].node = -1;
          frint("ELEMENT", &(resultdescr[i].element), &ierr);
          if (ierr==1) {
            resultdescr[i].element--;
          }
          else {
            dserror("Failed to read node or element number");
          }
        }

        frchar("POSITION", resultdescr[i].position, &ierr);
        if (ierr!=1) {
          dserror("Failed to read position specification");
        }

        frchar("NAME", resultdescr[i].name, &ierr);
        if (ierr!=1) {
          dserror("Failed to read name");
        }

        frdouble("VALUE", &(resultdescr[i].value), &ierr);
        if (ierr!=1) {
          dserror("Failed to read expected value");
        }

        frdouble("TOLERANCE", &(resultdescr[i].tolerance), &ierr);
        if (ierr!=1) {
          dserror("Failed to read tolerance");
        }

        frread();
      }
    }
    else {
      dserror("Unreliable input. Panic.");
    }
  }
  else
  {
    if (par.myrank==0)
    {
      printf("\n" MAGENTA_LIGHT "No result test section found. Skip testing." END_COLOR "\n");
    }
  }

end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif
