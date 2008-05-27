#ifdef CCADISCRET
#include "drt_utils_bspline.H"

//--------------------------------------------------
// Constructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(
  const int                     degree,
  const blitz::Array<double,1>  local_knotvector
  )
:
  myknotvector_(local_knotvector),
  degree_      (degree)
{
  return;
}

//--------------------------------------------------
// Copy constructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::BsplinePolynomial(
  const BsplinePolynomial& old
  )
  :
  degree_      (old.degree_)
{
  myknotvector_=old.myknotvector_;
  return;
}

//--------------------------------------------------
// Destructor
//--------------------------------------------------
DRT::NURBS::UTILS::BsplinePolynomial::~BsplinePolynomial()
{
  return;
}

//--------------------------------------------------
// Compute ldofid's Bspline value at point x
//--------------------------------------------------
void DRT::NURBS::UTILS::BsplinePolynomial::EvaluateBspline (
  double       & bspline_value,
  const double   x            ,
  const int      ldofid         
  )
{

  //                        ^
  //             ****       ^        +-----------+  
  //            *    *      ^        | ldofid==0 | 
  //           *      *     ^        +-----------+  
  //         **        **   ^                   
  //      ***            ***^   
  //  +***---+-----+-----+--***+-----+-----+-----+
  //                        ^
  //                        ^                    
  //                   **** ^        +-----------+             
  //                  *    *^        | ldofid==1 |             
  //                 *      *        +-----------+             
  //               **       ^**                    
  //            ***         ^  ***                
  //  +-----+***---+-----+-----+--***+-----+-----+
  //                        ^
  //                        ^
  //  +-----------+         ^****
  //  | ldofid==2 |         *    *
  //  +-----------+        *^     *
  //                     ** ^      **    
  //                  ***   ^        ***    
  //  +-----+-----+***---+-----+-----+--***+-----+
  //                        ^
  //                        ^   
  //  +-----------+         ^      ****
  //  | ldofid==3 |         ^     *    *
  //  +-----------+         ^    *      *
  //                        ^  **        **    
  //                        ***            ***    
  //  +-----+-----+-----+***---+-----+-----+--***+
  //                        ^
  //                        ^
  //                        x
    
#ifdef DEBUG
  // consistency vhecks
  if(0>ldofid || ldofid>degree_)
  {
    this->Throwerror("ldofid not in allowed (element dof) range [0;degree]\n");
  }

  // The node support looks like this:
  //
  //  |<----degree----->|     |<----degree----->|
  //  |                 |     |                 |
  //  |                 |     |                 |
  //
  //  +-----+-----+-----+-----+-----+-----+-----+
  //
  //                    |     |
  //                    |     |
  //                    |<--->|
  //
  //              element containing x
  //
  // check that we really have got something like this:

  if(myknotvector_.extent(blitz::firstDim)!=2*degree_+2)
  {
    std::string errorstring;
    
    errorstring.append("node support size is not 2*degree_+2\n");
    errorstring.append("Cannot compute bspline values\n");
    this->Throwerror(errorstring);
  }

  if(x<myknotvector_(degree_  )-1e-9
     ||
     x>myknotvector_(degree_+1)+1e-9)
  {
    this->Throwerror("Point not in evaluation interval\n");
  }
#endif
    
  // define the vector of values at x of all initial polynomials 
  vector<double> bspline(degree_+1);
    
  // The nonzero initial bspline polynomial and the intervals
  // that define the compact support of the bspline number lid
  // of given degree:
  // 
  //
  //  ldofid = 0:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //  +-----+-----+-----+--x--+
  //  
  //  all other initial bsplines are zero at x
  //  (empty intervals indicate the compact support of
  //   the other bspline polynomials which contribute
  //   to the bspline value associated with the given 
  //   dof lid. The union of all intervals defines the 
  //   support of the computed bspline of degree 3)
  //
  //
  //  ldofid = 1:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //        +-----+-----+--x--+-----+            
  //  
  //
  //  ldofid = 2:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //              +-----+--x--+-----+-----+       
  //  
  //
  //  ldofid = 3:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //                    +--x--+-----+-----+-----+
  //  
  //                    |     |
  //                    |     |
  //                    |<--->|
  //
  //              element containing x
  //

  for(int rr=0;rr<degree_+1;rr++)
  {
    bspline[rr]=0;
  }
  bspline[degree_-ldofid]=1;

  //        |        |        |        |
  //        | rr==0  | rr==1  | rr==2  |
  //        |        |        |        |
  //
  //      N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
  //        |       /|       /|       / 
  //        |      / |      / |      / 
  //        |     /  |     /  |     /  
  //        |    /   |    /   |    /   
  //        |   /    |   /    |   /    
  //        |  /     |  /     |  /     
  //        | /      | /      | /      
  //      N(0,1)  N(1,1)   N(2,1)             p == 1
  //        |       /|       / 
  //        |      / |      / 
  //        |     /  |     /  
  //        |    /   |    /   
  //        |   /    |   /    
  //        |  /     |  /     
  //        | /      | /      
  //      N(0,2)  N(1,2)                      p == 2
  //        |       / 
  //        |      / 
  //        |     /  
  //        |    /   
  //        |   /    
  //        |  /     
  //        | /      
  //      N(0,3)                              p == 3
  //
  //
  //
  // memory is reused, i.e. in the end, N(0,3) is contained
  // in bspline[0]
  //

  // loop all rows in the upper table
  for(int p=0;p<degree_;++p)
  {	
    // do computation of bspline values of specified degree,
    // corresponding to one row in the scheme above
    for(int rr=0;rr<degree_-p;++rr)
    {
      // id of first bspline function of this combination
      int  first = ldofid+rr;

      double fact1;
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(first+p+1)-myknotvector_(first))<10e-9)
      {
	fact1=0;
      }
      else
      {
	fact1 =(x-myknotvector_(first));
	fact1/=(myknotvector_(first+p+1)-myknotvector_(first));
      }

      double fact2;
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(first+p+2)-myknotvector_(first+1))<10e-9)
      {
	fact2=0;
      }
      else
      {
	fact2 =(myknotvector_(first+p+2)-x);
	fact2/=(myknotvector_(first+p+2)-myknotvector_(first+1));
      }
      // do the actual bspline recursion --- memory is reused!
      bspline[rr]=
	fact1*bspline[rr]
	+
	fact2*bspline[rr+1];
    }
  }

  // set the output
  bspline_value=bspline[0];
    
  return;
}


//--------------------------------------------------
// Compute ldofid's Bspline value at point x
// In addiditon, compute its first derivative
//--------------------------------------------------
void DRT::NURBS::UTILS::BsplinePolynomial::EvaluateBsplineAndDeriv(
  double       & bsplineval,
  double       & bsplineder,
  const double   x         ,
  const int      ldofid         
  )
{

  //                        ^
  //             ****       ^        +-----------+  
  //            *    *      ^        | ldofid==0 | 
  //           *      *     ^        +-----------+  
  //         **        **   ^                   
  //      ***            ***^   
  //  +***---+-----+-----+--***+-----+-----+-----+
  //                        ^
  //                        ^                    
  //                   **** ^        +-----------+             
  //                  *    *^        | ldofid==1 |             
  //                 *      *        +-----------+             
  //               **       ^**                    
  //            ***         ^  ***                
  //  +-----+***---+-----+-----+--***+-----+-----+
  //                        ^
  //                        ^
  //  +-----------+         ^****
  //  | ldofid==2 |         *    *
  //  +-----------+        *^     *
  //                     ** ^      **    
  //                  ***   ^        ***    
  //  +-----+-----+***---+-----+-----+--***+-----+
  //                        ^
  //                        ^   
  //  +-----------+         ^      ****
  //  | ldofid==3 |         ^     *    *
  //  +-----------+         ^    *      *
  //                        ^  **        **    
  //                        ***            ***    
  //  +-----+-----+-----+***---+-----+-----+--***+
  //                        ^
  //                        ^
  //                        x

#ifdef DEBUG    
  // consistency vhecks
  if(0>ldofid || ldofid>degree_)
  {
    this->Throwerror("lodfid not in allowed (element dof) range [0;degree]\n");
  }

  // The node support looks like this:
  //
  //  |<----degree----->|     |<----degree----->|
  //  |                 |     |                 |
  //  |                 |     |                 |
  //
  //  +-----+-----+-----+-----+-----+-----+-----+
  //
  //                    |     |
  //                    |     |
  //                    |<--->|
  //
  //              element containing x
  //
  // check that we really have got something like this:

  if(myknotvector_.extent(blitz::firstDim)!=2*degree_+2)
  {
    std::string errorstring;
    
    errorstring.append("node support size is not 2*degree_+2\n");
    errorstring.append("Cannot compute bspline values\n");
    this->Throwerror(errorstring);
  }

  if(x<myknotvector_(degree_  )-1e-9
     ||
     x>myknotvector_(degree_+1)+1e-9)
  {
    this->Throwerror("Point not in Evaluation interval\n");
  }
#endif

  // define the vector of values at x of all initial polynomials 
  vector<double> bspline(degree_+1);
    
  // The nonzero initial bspline polynomial and the intervals
  // that define the compact support of the bspline number lid
  // of given degree:
  // 
  //
  //  ldofid = 0:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //  +-----+-----+-----+--x--+
  //  
  //  all other initial bsplines are zero at x
  //  (empty intervals indicate the compact support of
  //   the other bspline polynomials which contribute
  //   to the bspline value associated with the given 
  //   dof lid. The union of all intervals defines the 
  //   support of the computed bspline of degree 3)
  //
  //
  //  ldofid = 1:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //        +-----+-----+--x--+-----+            
  //  
  //
  //  ldofid = 2:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //              +-----+--x--+-----+-----+       
  //  
  //
  //  ldofid = 3:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //                    +--x--+-----+-----+-----+
  //  
  //                    |     |
  //                    |     |
  //                    |<--->|
  //
  //              element containing x
  //

  // initial values for recursion
  //
  //           +-
  //           | 
  //    0      |  1   for x contained in interval i
  //   N (x) = |
  //    i      |  0   otherwise 
  //           |
  //           +-
  //

  for(int rr=0;rr<degree_+1;rr++)
  {
    bspline[rr]=0;
  }
  bspline[degree_-ldofid]=1;

  //        |        |        |        |
  //        | rr==0  | rr==1  | rr==2  |
  //        |        |        |        |
  //
  //      N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
  //        |       /|       /|       / 
  //        |      / |      / |      / 
  //        |     /  |     /  |     /  
  //        |    /   |    /   |    /   
  //        |   /    |   /    |   /    
  //        |  /     |  /     |  /     
  //        | /      | /      | /      
  //      N(0,1)  N(1,1)   N(2,1)             p == 1
  //        |       /|       / 
  //        |      / |      / 
  //        |     /  |     /  
  //        |    /   |    /   
  //        |   /    |   /    
  //        |  /     |  /     
  //        | /      | /      
  //      N(0,2)  N(1,2)                      p == 2
  //
  //
  //
  // ----------------------------------------------------
  //
  //
  //
  //      N(0,2)  N(1,2)                      p == 2   +--
  //        |\      /|                                 |
  //        | \    / |                                 | branch
  //        |  \  /  |                                 |
  //        |   \/   |                                 | for first
  //        |   /\   |                                 |
  //        |  /  \  |                                 | derivatives
  //        | /    \ |                                 |
  //     val(0,3) der(0,3)                    p == 3   +--
  //
  //
  //
  // memory is reused on the first level. For the last
  // step, we have additional memory tobe able to access
  // N(0,3) twice
  //

  // loop all rows in the upper table up to the last 
  // but one. Both arguments are still required to compute
  // the derivatives, so do not throw them away 
  // (or overwrite)
  for(int p=0;p<degree_-1;++p)
  {	
    // do computation of bspline values of specified degree,
    // corresponding to one row in the scheme above
    for(int rr=0;rr<degree_-p;++rr)
    {
      // id of first bspline function of this combination
      int  i = ldofid+rr;

      // recursion for the computation of the basis 
      // function
      //
      //          x - x                x     - x            
      //  p            i     p-1        i+p+1         p-1     
      // N (x) = -------- * N   (x) + ------------ * N   (x)  
      //  i      x   - x     i        x     - x       i+1     
      //          i+p   i              i+p+1   i+1          
      //
      //        |        |           |            |
      //        +--------+           +------------+
      //           fact1                  fact2
      //

      double fact1;
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(i+p+1)-myknotvector_(i))<10e-9)
      {
	fact1=0;
      }
      else
      {
	fact1 =(x-myknotvector_(i));
	fact1/=(myknotvector_(i+p+1)-myknotvector_(i));
      }

      double fact2;
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(i+p+2)-myknotvector_(i+1))<10e-9)
      {
	fact2=0;
      }
      else
      {
	fact2 =(myknotvector_(i+p+2)-x);
	fact2/=(myknotvector_(i+p+2)-myknotvector_(i+1));
      }
      // do the actual bspline recursion --- memory is reused!
      bspline[rr]=
	fact1*bspline[rr]
	+
	fact2*bspline[rr+1];
    }
  }
  
  //---------------------------------------------------
  // do computation of bspline value in the last level 
  // corresponding to one row in the scheme above
    
  double fact1;
  // the first part of the if statement allows to 
  // enforce interpolation using multiple nodes
  // the second part is a part of the bspline recursion
  if(fabs(myknotvector_(ldofid+degree_)-myknotvector_(ldofid))<10e-9)
  {
    fact1=0;
  }
  else
  {
    fact1 =(x-myknotvector_(ldofid));
    fact1/=(myknotvector_(ldofid+degree_)-myknotvector_(ldofid));
  }
  
  double fact2;
  // the first part of the if statement allows to 
  // enforce interpolation using multiple nodes
  // the second part is a part of the bspline recursion
  if(fabs(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1))<10e-9)
  {
    fact2=0;
  }
  else
  {
    fact2 =(myknotvector_(ldofid+degree_+1)-x);
    fact2/=(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1));
  }
  // do the actual bspline recursion --- memory is reused!
  bsplineval=
    fact1*bspline[0]
    +
    fact2*bspline[1];

  //---------------------------------------------------
  // do computation of bspline derivatives from the 
  // last level corresponding to one row in the scheme 
  // above

  if(degree_>0)
  {
    //
    //                            
    //   p          p        p-1           p         p-1     
    // N'  (x) = -------- * N   (x) - ----------- * N   (x)  
    //   i       x   - x     i        x     - x      i+1     
    //            i+p   i              i+p+1   i+1          
    //
    //          |        |           |            |
    //          +--------+           +------------+
    //             fact1                  fact2
    
    // the first part of the if statement allows to 
    // enforce interpolation using multiple nodes
    // the second part computes fact1 in the equation
    // above
    if(fabs(myknotvector_(ldofid+degree_)-myknotvector_(ldofid))<10e-9)
    {
      fact1=0;
    }
    else
    {
      fact1 =degree_;
      fact1/=(myknotvector_(ldofid+degree_)-myknotvector_(ldofid));
    }
    
    // the first part of the if statement allows to 
    // enforce interpolation using multiple nodes
    // the second part is a part of the bspline recursion
    // see above
    if(fabs(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1))<10e-9)
    {
      fact2=0;
    }
    else
    {
      fact2 =degree_;
      fact2/=(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1));
    }
    
    // compute the actual bspline derivative formula
    bsplineder=
      fact1*bspline[0]
      -
      fact2*bspline[1];
  }
  else
  {
    // piecewise constants get 0 derivatives 
    bsplineder=0;
  }

  return;
}

//--------------------------------------------------
// Compute ldofid's Bspline value at point x
// In addiditon, compute its first and second
// derivative
//--------------------------------------------------
void DRT::NURBS::UTILS::BsplinePolynomial::EvaluateBsplineFirstAndSecondDeriv(
  double &             bsplineval  ,
  double &             bsplineder  ,
  double &             bsplineder2 ,
  const double         x           ,
  const int            ldofid         
  )
{


  /*
 
            x - x                x     - x            
    p            i     p-1        i+p+1         p-1     
   N (x) = -------- * N   (x) + ------------ * N   (x)  
    i      x   - x     i        x     - x       i+1     
            i+p   i              i+p+1   i+1          

     p          p        p-1           p         p-1     
   N'  (x) = -------- * N   (x) - ----------- * N   (x)  
     i       x   - x     i        x     - x      i+1     
              i+p   i              i+p+1   i+1          

  -------------------------------------------------------

     p          p        p-1           p         p-1     
   N'' (x) = -------- * N'  (x) - ----------- * N'  (x)  
     i       x   - x     i        x     - x      i+1     
              i+p   i              i+p+1   i+1          

     p-1         p-1       p-2          p-1        p-2     
   N'  (x) = ---------- * N   (x) - ----------- * N   (x)  
     i       x     - x     i        x   - x        i+1     
              i+p-1   i              i+p   i+1          


  */

  //                        ^
  //             ****       ^        +-----------+  
  //            *    *      ^        | ldofid==0 | 
  //           *      *     ^        +-----------+  
  //         **        **   ^                   
  //      ***            ***^   
  //  +***---+-----+-----+--***+-----+-----+-----+
  //                        ^
  //                        ^                    
  //                   **** ^        +-----------+             
  //                  *    *^        | ldofid==1 |             
  //                 *      *        +-----------+             
  //               **       ^**                    
  //            ***         ^  ***                
  //  +-----+***---+-----+-----+--***+-----+-----+
  //                        ^
  //                        ^
  //  +-----------+         ^****
  //  | ldofid==2 |         *    *
  //  +-----------+        *^     *
  //                     ** ^      **    
  //                  ***   ^        ***    
  //  +-----+-----+***---+-----+-----+--***+-----+
  //                        ^
  //                        ^   
  //  +-----------+         ^      ****
  //  | ldofid==3 |         ^     *    *
  //  +-----------+         ^    *      *
  //                        ^  **        **    
  //                        ***            ***    
  //  +-----+-----+-----+***---+-----+-----+--***+
  //                        ^
  //                        ^
  //                        x
    
#ifdef DEBUG
  // consistency vhecks
  if(0>ldofid || ldofid>degree_)
  {
    this->Throwerror("lodfid not in allowed ldofid (element dof) range [0;degree]\n");
  }

  // The node support looks like this:
  //
  //  |<----degree----->|     |<----degree----->|
  //  |                 |     |                 |
  //  |                 |     |                 |
  //
  //  +-----+-----+-----+-----+-----+-----+-----+
  //
  //                    |     |
  //                    |     |
  //                    |<--->|
  //
  //              element containing x
  //
  // check that we really have got something like this:

  if(myknotvector_.extent(blitz::firstDim)!=2*degree_+2)
  {
    std::string errorstring;
    
    errorstring.append("node support size is not 2*degree_+2\n");
    errorstring.append("Cannot compute bspline values\n");
    this->Throwerror(errorstring);
  }

  if(x<myknotvector_(degree_  )-1e-9
     ||
     x>myknotvector_(degree_+1)+1e-9)
  {
    this->Throwerror("Point not in evaluation interval\n");
  }
#endif
    
  // define the vector of values at x of all initial polynomials 
  vector<double> bspline(degree_+1);
    
  // The nonzero initial bspline polynomial and the intervals
  // that define the compact support of the bspline number lid
  // of given degree:
  // 
  //
  //  ldofid = 0:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //  +-----+-----+-----+--x--+
  //  
  //  all other initial bsplines are zero at x
  //  (empty intervals indicate the compact support of
  //   the other bspline polynomials which contribute
  //   to the bspline value associated with the given 
  //   dof lid. The union of all intervals defines the 
  //   support of the computed bspline of degree 3)
  //
  //
  //  ldofid = 1:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //        +-----+-----+--x--+-----+            
  //  
  //
  //  ldofid = 2:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //              +-----+--x--+-----+-----+       
  //  
  //
  //  ldofid = 3:
  //                    +-----+
  //                    |     |
  //                    |     |
  //                    |     |
  //                    +--x--+-----+-----+-----+
  //  
  //                    |     |
  //                    |     |
  //                    |<--->|
  //
  //              element containing x
  //

  // initial values for recursion
  //
  //           +-
  //           | 
  //    0      |  1   for x contained in interval i
  //   N (x) = |
  //    i      |  0   otherwise 
  //           |
  //           +-
  //

  for(int rr=0;rr<degree_+1;rr++)
  {
    bspline[rr]=0;
  }
  bspline[degree_-ldofid]=1;

  //        |        |        |        |
  //        | rr==0  | rr==1  | rr==2  |
  //        |        |        |        |
  //
  //      N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
  //        |       /|       /|       / 
  //        |      / |      / |      / 
  //        |     /  |     /  |     /  
  //        |    /   |    /   |    /   
  //        |   /    |   /    |   /    
  //        |  /     |  /     |  /     
  //        | /      | /      | /      
  //      N(0,1)  N(1,1)   N(2,1)             p == 1
  //
  //
  //
  // ====================================================
  //
  //
  //
  //      N(0,1)  N(1,1)   N(2,1)             p == 1   +--	    
  //        |       /|       / 			       |	    
  //        |      / |      / 			       | branch	    
  //        |     /  |     /  			       |	    
  //        |    /   |    /   			       | for second  
  //        |   /    |   /    			       |	    
  //        |  /     |  /     			       | derivatives
  //        | /      | /      			       |	    
  //      N(0,2)  N(1,2)                      p == 2   +--          
  //
  //
  //
  // ====================================================
  //
  //
  //
  //      N(0,2)  N(1,2)                      p == 2   +--
  //        |\      /|                                 |
  //        | \    / |                                 | branch
  //        |  \  /  |                                 |
  //        |   \/   |                                 | for first
  //        |   /\   |                                 |
  //        |  /  \  |                                 | derivatives
  //        | /    \ |                                 |
  //     val(0,3) der(0,3)                    p == 3   +--          
  //
  //
  //
  // memory is reused on the first level. For the last
  // step, we have additional memory tobe able to access
  // N(0,3) twice
  //

  // temps 
  double fact1;
  double fact2;

  // loop all rows in the upper table up to the last 
  // but one. Both arguments are still required to compute
  // the derivatives, so do not throw them away 
  // (or overwrite)
  for(int p=0;p<degree_-2;++p)
  {	
    // do computation of bspline values of specified degree,
    // corresponding to one row in the scheme above
    for(int rr=0;rr<degree_-p;++rr)
    {
      // id of first bspline function of this combination
      int  i = ldofid+rr;

      // recursion for the computation of the basis 
      // function
      //
      //          x - x                x     - x            
      //  p            i     p-1        i+p+1         p-1     
      // N (x) = -------- * N   (x) + ------------ * N   (x)  
      //  i      x   - x     i        x     - x       i+1     
      //          i+p   i              i+p+1   i+1          
      //
      //        |        |           |            |
      //        +--------+           +------------+
      //           fact1                  fact2
      //

      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(i+p+1)-myknotvector_(i))<10e-9)
      {
	fact1=0;
      }
      else
      {
	fact1 =(x-myknotvector_(i));
	fact1/=(myknotvector_(i+p+1)-myknotvector_(i));
      }

      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(i+p+2)-myknotvector_(i+1))<10e-9)
      {
	fact2=0;
      }
      else
      {
	fact2 =(myknotvector_(i+p+2)-x);
	fact2/=(myknotvector_(i+p+2)-myknotvector_(i+1));
      }
      // do the actual bspline recursion --- memory is reused!
      bspline[rr]=
	fact1*bspline[rr]
	+
	fact2*bspline[rr+1];
    }
  }

  //
  // ====================================================
  //

  //---------------------------------------------------
  // do computation of both bspline derivatives  
  // from the p-1 level 

  blitz::Array<double,1> pmo_deriv(2);

  if(degree_>1)
  {
    for(int rr=0;rr<2;++rr)
    {
      // id of first bspline function of this combination
      int  i = ldofid+rr;

      //
      //                            
      //
      //     p-1         p-1       p-2          p-1        p-2     
      //   N'  (x) = ---------- * N   (x) - ----------- * N   (x)  
      //     i       x     - x     i        x   - x        i+1     
      //              i+p-1   i              i+p   i+1          
      //             |        |           |            |
      //             +--------+           +------------+
      //                fact1                  fact2
      
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part computes fact1 in the equation
      // above
      if(fabs(myknotvector_(i+degree_-1)-myknotvector_(i))<10e-9)
      {
	fact1=0;
      }
      else
      {
	fact1 =degree_-1;
	fact1/=(myknotvector_(i+degree_-1)-myknotvector_(i));
      }
      
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      // see above
      if(fabs(myknotvector_(i+degree_)-myknotvector_(i+1))<10e-9)
      {
	fact2=0;
      }
      else
      {
	fact2 =degree_-1;
	fact2/=(myknotvector_(i+degree_)-myknotvector_(i+1));
      }
      
      // compute the actual bspline derivative formula
      pmo_deriv(rr)=
	fact1*bspline[rr]
	-
	fact2*bspline[rr+1];
    }

    // do computation of bspline values of degree-1,
    
    for(int rr=0;rr<2;++rr)
    {
      // id of first bspline function of this combination
      int  i = ldofid+rr;
      
      // recursion for the computation of the basis 
      // function
      //
      //          x - x                x     - x            
      //  p            i     p-1        i+p+1         p-1     
      // N (x) = -------- * N   (x) + ------------ * N   (x)  
      //  i      x   - x     i        x     - x       i+1     
      //          i+p   i              i+p+1   i+1          
      //
      //        |        |           |            |
      //        +--------+           +------------+
      //           fact1                  fact2
      //
      
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(i+degree_-1)-myknotvector_(i))<10e-9)
      {
	fact1=0;
      }
      else
      {
	fact1 =(x-myknotvector_(i));
	fact1/=(myknotvector_(i+degree_-1)-myknotvector_(i));
      }
      
      // the first part of the if statement allows to 
      // enforce interpolation using multiple nodes
      // the second part is a part of the bspline recursion
      if(fabs(myknotvector_(i+degree_)-myknotvector_(i+1))<10e-9)
      {
	fact2=0;
      }
      else
      {
	fact2 =(myknotvector_(i+degree_)-x);
	fact2/=(myknotvector_(i+degree_)-myknotvector_(i+1));
      }
      // do the actual bspline recursion --- memory is reused!
      bspline[rr]=
	fact1*bspline[rr]
	+
	fact2*bspline[rr+1];
    }
  }
  else
  {
    // piecewise constants get 0 derivatives 
    pmo_deriv=0;
  }


  //---------------------------------------------------
  // evaluate the second derivatives

  //   p          p        p-1           p         p-1     
  // N'' (x) = -------- * N'  (x) - ----------- * N'  (x)  
  //   i       x   - x     i        x     - x      i+1     
  //            i+p   i              i+p+1   i+1          
  //
  //          |        |           |            |
  //          +--------+           +------------+
  //             fact1                  fact2
  
  // the first part of the if statement allows to 
  // enforce interpolation using multiple nodes
  // the second part computes fact1 in the equation
  // above
  if(fabs(myknotvector_(ldofid+degree_)-myknotvector_(ldofid))<10e-9)
  {
    fact1=0;
  }
  else
  {
    fact1 =degree_;
    fact1/=(myknotvector_(ldofid+degree_)-myknotvector_(ldofid));
  }
    
  // the first part of the if statement allows to 
  // enforce interpolation using multiple nodes
  // the second part is a part of the bspline recursion
  // see above
  if(fabs(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1))<10e-9)
  {
    fact2=0;
  }
  else
  {
    fact2 =degree_;
    fact2/=(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1));
  }
    
  // compute the actual bspline derivative formula
  bsplineder2=
    fact1*pmo_deriv(0)
    -
    fact2*pmo_deriv(1);

  //
  // ====================================================
  //

  //---------------------------------------------------
  // do computation of bspline value in the last level 
    
  // the first part of the if statement allows to 
  // enforce interpolation using multiple nodes
  // the second part is a part of the bspline recursion
  if(fabs(myknotvector_(ldofid+degree_)-myknotvector_(ldofid))<10e-9)
  {
    fact1=0;
  }
  else
  {
    fact1 =(x-myknotvector_(ldofid));
    fact1/=(myknotvector_(ldofid+degree_)-myknotvector_(ldofid));
  }
  
  // the first part of the if statement allows to 
  // enforce interpolation using multiple nodes
  // the second part is a part of the bspline recursion
  if(fabs(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1))<10e-9)
  {
    fact2=0;
  }
  else
  {
    fact2 =(myknotvector_(ldofid+degree_+1)-x);
    fact2/=(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1));
  }
  // do the actual bspline recursion --- memory is reused!
  bsplineval=
    fact1*bspline[0]
    +
    fact2*bspline[1];

  //---------------------------------------------------
  // do computation of bspline derivatives from the 
  // last level corresponding to one row in the scheme 
  // above

  if(degree_>0)
  {
    //
    //                            
    //   p          p        p-1           p         p-1     
    // N'  (x) = -------- * N   (x) - ----------- * N   (x)  
    //   i       x   - x     i        x     - x      i+1     
    //            i+p   i              i+p+1   i+1          
    //
    //          |        |           |            |
    //          +--------+           +------------+
    //             fact1                  fact2
    
    // the first part of the if statement allows to 
    // enforce interpolation using multiple nodes
    // the second part computes fact1 in the equation
    // above
    if(fabs(myknotvector_(ldofid+degree_)-myknotvector_(ldofid))<10e-9)
    {
      fact1=0;
    }
    else
    {
      fact1 =degree_;
      fact1/=(myknotvector_(ldofid+degree_)-myknotvector_(ldofid));
    }
    
    // the first part of the if statement allows to 
    // enforce interpolation using multiple nodes
    // the second part is a part of the bspline recursion
    // see above
    if(fabs(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1))<10e-9)
    {
      fact2=0;
    }
    else
    {
      fact2 =degree_;
      fact2/=(myknotvector_(ldofid+degree_+1)-myknotvector_(ldofid+1));
    }
    
    // compute the actual bspline derivative formula
    bsplineder=
      fact1*bspline[0]
      -
      fact2*bspline[1];
  }
  else
  {
    // piecewise constants get 0 derivatives 
    bsplineder=0;
  }
  
  return;
}

#endif
