/* NOT USED ANY MORE !!!!!!!!!!!!!!!!!!!! 
   YOU CAN FIND THE ACUTAL fluid.h in ../headers 
*/   
#ifdef D_FLUID2
#include "../fluid2/fluid2.h"
typedef union _FLUID_DATA
{
   struct _F2_DATA  f2data;
} FLUID_DATA;
#else
typedef union _FLUID_DATA
{
   int          dummy;
} FLUID_DATA;
#endif
