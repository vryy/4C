#define MAX(x,y)  ((x) < (y) ? (y) : (x))
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#define ABS(x)    ((x) <  0  ? (-x) : (x))
#define FABS(x)   ((x) < 0.0 ? (-(x)) : (x))
#define RAD       (atan(1.0)/45.0)
#define PI        (asin(1.0)*2.0)
#define MAXNOD    (27)
typedef struct _NODE
{
     int                        Id;             
} NODE;
typedef struct _ELEMENT
{  
     int                        Id;                 
} ELEMENT;
