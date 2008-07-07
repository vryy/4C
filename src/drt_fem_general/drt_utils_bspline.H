#ifdef CCADISCRET

#ifndef BSPLINE_H
#define BSPLINE_H

#include<math.h>
#include<string>
#include<vector>
#include<iostream>
#include<fstream>

#include <blitz/array.h>
#include "../drt_lib/drt_dserror.H"

using namespace std;


namespace DRT
{
namespace NURBS
{
namespace UTILS
{

/*!
\file bspline.H
\class BsplinePolynomial
\brief A class which allows to evaluate all bspline
       polynomials which are nonzero in the 'central'
       knotvector interval
*/
class BsplinePolynomial
{
public:
  //--------------------------------------------------
  /*! \brief constructor  
  
      \param degree           int (i)
                              degree of polynomial, for
                              consistency checks
  
      \param local_knotvector blitz::Array<double,1> (i)
                              knot range that contains 
                              the compact support of 
                              all bsplines which are  
                              not zero in x 
  
   
      The knot support looks like this:

      \verbatim
 
      |<----degree----->|     |<----degree----->|
      |                 |     |                 |
      |                 |     |                 |
  
      +-----+-----+-----+-----+-----+-----+-----+
  
                        |     |
                        |     |
                        |<--->|
    
                  element containing x
		  
      \endverbatim

  */
  //--------------------------------------------------
  BsplinePolynomial (
    const int                     degree,
    const blitz::Array<double,1>  local_knotvector
    );
  //--------------------------------------------------
  //! \brief destructor
  //--------------------------------------------------
  virtual ~BsplinePolynomial();
  //--------------------------------------------------
  //! \brief copy constructor
  //! 
  //! \param old const BsplinePolynomial (i)  
  //!            bspline polynomial to copy
  //--------------------------------------------------
  BsplinePolynomial(const BsplinePolynomial& old);

  //--------------------------------------------------
  /*! \brief compute a bspline value at position x
  
   \param bspline_value double (o)
                        the bspline value computed
  
   \param  x            double (i) 
                        position to evaluate bspline
  
   \param  ldofid       int    (i) 
                        which bspline to evaluate,
                        ldofid is the local dof id,
                        i.e. the number of the 
                        bspline we want to evaluate.
                        0 is the leftmost nonzero 
                        bspline, degree the right 
                        most. All of them have 
                        degree+1 nonzero intervals

  
   example for the use of ldofid: (degree is 3) 

   \verbatim
                          ^                    
               ****       ^        +-----------+  
              *    *      ^        | ldofid==0 | 
             *      *     ^        +-----------+  
           **        **   ^                    
        ***            ***^                    
    +***---+-----+-----+--***+-----+-----+-----+
                          ^
                          ^                    
                     **** ^        +-----------+             
                    *    *^        | ldofid==1 |             
                   *      *        +-----------+             
                 **       ^**                    
              ***         ^  ***                
    +-----+***---+-----+-----+--***+-----+-----+
                          ^
                          ^
    +-----------+         ^****
    | ldofid==2 |         *    *
    +-----------+        *^     *
                       ** ^      **    
                    ***   ^        ***    
    +-----+-----+***---+-----+-----+--***+-----+
                          ^
                          ^   
    +-----------+         ^      ****
    | ldofid==3 |         ^     *    *
    +-----------+         ^    *      *
                          ^  **        **       
                          ***            ***
    +-----+-----+-----+***---+-----+-----+--***+
                          ^
                          ^
                          x

   \endverbatim


   Recursion formula used:
  
   \verbatim
 
            x - x                x     - x            
    p            i     p-1        i+p+1         p-1     
   N (x) = -------- * N   (x) + ------------ * N   (x)  
    i      x   - x     i        x     - x       i+1     
            i+p   i              i+p+1   i+1          

   \endverbatim

   The recursion is started from the following functions
   (depending on ldofid). 

   \verbatim
  
             +-
             | 
      0      |  1   for x contained in interval i
     N (x) = |
      i      |  0   otherwise 
             |
             +-
  
   \endverbatim

   In the example above (degree 3), we have:
   
   \verbatim
   
  
    ldofid = 0:
                      +-----+
                      |     |
                      |     |
                      |     |
    +-----+-----+-----+--x--+
    

  
    ldofid = 1:
                      +-----+
                      |     |
                      |     |
                      |     |
          +-----+-----+--x--+-----+            
    
  
    ldofid = 2:
                      +-----+
                      |     |
                      |     |
                      |     |
                +-----+--x--+-----+-----+       
    
  
    ldofid = 3:
                      +-----+
                      |     |
                      |     |
                      |     |
                      +--x--+-----+-----+-----+
    
                      |     |
                      |     |
                      |<--->|
  
                element containing x
  
   \endverbatim



   the recursion is done in Aitken-Neville-style:

   \verbatim
          |        |        |        |
          | rr==0  | rr==1  | rr==2  |
          |        |        |        |
  
        N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
          |       /|       /|       / 
          |      / |      / |      / 
          |     /  |     /  |     /  
          |    /   |    /   |    /   
          |   /    |   /    |   /    
          |  /     |  /     |  /     
          | /      | /      | /      
        N(0,1)  N(1,1)   N(2,1)             p == 1
          |       /|       / 
          |      / |      / 
          |     /  |     /  
          |    /   |    /   
          |   /    |   /    
          |  /     |  /     
          | /      | /      
        N(0,2)  N(1,2)                      p == 2
          |       / 
          |      / 
          |     /  
          |    /   
          |   /    
          |  /     
          | /      
        N(0,3)                              p == 3
   
   \endverbatim
  
   memory is reused, i.e. in the end, N(0,3) is 
   contained in bspline[0]

  */
  //--------------------------------------------------
  void EvaluateBspline (
    double       & bspline_value,
    const double   x            ,
    const int      ldofid         
    );

  //--------------------------------------------------
  /*! \brief Compute ldofid's Bspline value at point x
             In addiditon, compute its first derivative
  
   \param bsplineval   double   (o)
                       bspline value

   \param bsplineder   double   (o) 
                       first derivative of bspline
  
   \param x            double   (i) 
                       position to evaluate bspline
  
   \param ldofid       int      (i) 
                       which bspline to evaluate,
                       ldofid is the local dof id,
                       i.e. the number of the 
                       bspline we want to evaluate.
                       0 is the leftmost nonzero 
                       bspline, degree the right 
                       most. All of them have 
                       degree+1 nonzero intervals

   Recursion formulas used:
  
   \verbatim
 
            x - x                x     - x            
    p            i     p-1        i+p+1         p-1     
   N (x) = -------- * N   (x) + ------------ * N   (x)  
    i      x   - x     i        x     - x       i+1     
            i+p   i              i+p+1   i+1          

     p          p        p-1           p         p-1     
   N'  (x) = -------- * N   (x) - ----------- * N   (x)  
     i       x   - x     i        x     - x      i+1     
              i+p   i              i+p+1   i+1          

   \endverbatim

   mind that the recursion scheme is changed
   (we have to descent into a branch for the 
    first derivatives):

   \verbatim

          |        |        |        |
          | rr==0  | rr==1  | rr==2  |
          |        |        |        |
  
        N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
          |       /|       /|       / 
          |      / |      / |      / 
          |     /  |     /  |     /  
          |    /   |    /   |    /   
          |   /    |   /    |   /    
          |  /     |  /     |  /     
          | /      | /      | /      
        N(0,1)  N(1,1)   N(2,1)             p == 1
          |       /|       / 
          |      / |      / 
          |     /  |     /  
          |    /   |    /   
          |   /    |   /    
          |  /     |  /     
          | /      | /      
        N(0,2)  N(1,2)                      p == 2
  
  
  
   ----------------------------------------------------
  
  
  
        N(0,2)  N(1,2)                      p == 2   +--
          |\      /|                                 |
          | \    / |                                 | branch
          |  \  /  |                                 |
          |   \/   |                                 | for first
          |   /\   |                                 |
          |  /  \  |                                 | derivatives
          | /    \ |                                 |
       val(0,3) der(0,3)                    p == 3   +--
  
  
  
   memory is reused on the first level. For the last
   step, we have additional memory to be able to access
   N(0,3) twice
  

   \endverbatim

  */
  void EvaluateBsplineAndDeriv(
    double       & bsplineval,
    double       & bsplineder,
    const double   x         ,
    const int      ldofid         
    );

  //--------------------------------------------------
  /*! \brief Compute ldofid's Bspline value at point x
             In addiditon, compute its first and
	     second derivative
  
  \param bsplineval  double  (o)
                     bspline value

  \param bsplineder  double  (o)
                     first derivative of bspline

  \param bsplineder2 double  (o)
                     second derivative of bspline
  
  \param x           double  (i)
                     position to evaluate bspline
  
  \param ldofid      int     (i)
                     which bspline to evaluate,
	 	     ldofid is the local dof id,
	 	     i.e. the number of the 
		     bspline we want to evaluate.
		     0 is the leftmost nonzero 
		     bspline, degree the right 
		     most. All of them have 
		     degree+1 nonzero intervals

   (Recursion) formulas:
  
   \verbatim

   Recursion for Basis functions:
 
            x - x                x     - x            
    p            i     p-1        i+p+1         p-1     
   N (x) = -------- * N   (x) + ------------ * N   (x)  
    i      x   - x     i        x     - x       i+1     
            i+p   i              i+p+1   i+1          

   Computation of first derivatives from basis functions
   of the previous level:

     p          p        p-1           p         p-1     
   N'  (x) = -------- * N   (x) - ----------- * N   (x)  
     i       x   - x     i        x     - x      i+1     
              i+p   i              i+p+1   i+1          

   Computation of second derivatives:

  -------------------------------------------------------

     p          p        p-1           p         p-1     
   N'' (x) = -------- * N'  (x) - ----------- * N'  (x)  
     i       x   - x     i        x     - x      i+1     
              i+p   i              i+p+1   i+1          

   with

     p-1         p-1       p-2          p-1        p-2     
   N'  (x) = ---------- * N   (x) - ----------- * N   (x)  
     i       x     - x     i        x   - x        i+1     
              i+p-1   i              i+p   i+1          

   \endverbatim
   

   \verbatim


          |        |        |        |
          | rr==0  | rr==1  | rr==2  |
          |        |        |        |
  
        N(0,0)  N(1,0)   N(2,0)   N(3,0)    p == 0
          |       /|       /|       / 
          |      / |      / |      / 
          |     /  |     /  |     /  
          |    /   |    /   |    /   
          |   /    |   /    |   /    
          |  /     |  /     |  /     
          | /      | /      | /      
        N(0,1)  N(1,1)   N(2,1)             p == 1
  
  
  
   ====================================================
  
  
  
        N(0,1)  N(1,1)   N(2,1)             p == 1   +--	    
          |       /|       / 			       |	    
          |      / |      / 			       | branch	    
          |     /  |     /  			       |	    
          |    /   |    /   			       | for second  
          |   /    |   /    			       |	    
          |  /     |  /     			       | derivatives
          | /      | /      			       |	    
        N(0,2)  N(1,2)                      p == 2   +--          
  
  
  
   ====================================================
  
  
  
        N(0,2)  N(1,2)                      p == 2   +--
          |\      /|                                 |
          | \    / |                                 | branch
          |  \  /  |                                 |
          |   \/   |                                 | for first
          |   /\   |                                 |
          |  /  \  |                                 | derivatives
          | /    \ |                                 |
       val(0,3) der(0,3)                    p == 3   +--          

   \endverbatim

  */
  void EvaluateBsplineFirstAndSecondDeriv(
    double &             bsplineval ,
    double &             bsplineder ,
    double &             bsplineder2,
    const double         x          ,
    const int            ldofid
    );

  //--------------------------------------------------
  //! \brief Print some information on the bspline
  //--------------------------------------------------
  void   PrintBspline()
    {
      cout << "------------------\n";
      cout << "Bspline knotvector\n";
      cout << myknotvector_;
      cout << &endl;
      cout << "Bspline degree: " << degree_ <<"\n";
      cout << &endl;
      cout << "Allows to evaluate bsplines at positions ";
      cout << "in interval [";
      cout << myknotvector_(degree_  );
      cout << ",";
      cout << myknotvector_(degree_+1);
      cout << "]\n";

      return;
    };
  
private:

  //! the part of the knotvector where bspline 
  //! polynomials of the center-interval are non-zero
  blitz::Array<double,1> myknotvector_;

  //! a working array for the construction of bsplines
  vector<double>         bspline_;

  //! the degree of the bspline-polynomial under 
  //! consideration
  int                    degree_;

  //! temporary doubles
  double                 fact_[4];

  //! both bspline derivatives from the p-1 level
  double                 pmo_deriv_[2];

  //--------------------------------------------------
  //! \brief throw an errormessage
  //--------------------------------------------------
  void Throwerror(const string errormessage)
    {
      // give some information on bspline
      PrintBspline();

      // and the throw the error and exit with a
      // sigsegv
      dserror(errormessage.c_str());
      return;
    };

};

} // namespace UTILS 
} // namespace NURBS 
} // namespace DRT

#endif
#endif
