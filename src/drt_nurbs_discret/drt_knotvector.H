/*!----------------------------------------------------------------------
\file drt_knotvector.H

\brief knot vectors for nurbs problems (isogeometric analysis)

       the class is containing the data structures + some 
       service functions (do we have interpolation, is it an 
       open knot vector or periodic, access methods etc)

       ParObject is implemented to be able to write the knotvector
       to disc for io

<pre>
Maintainer: Peter Gamnitzer
            gammi@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef KNOTVECTOR_H
#define KNOTVECTOR_H

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_dserror.H"

#include "Teuchos_RefCountPtr.hpp"


using namespace std;


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

namespace NURBS
{
  
/*!
\brief A class to manage a nurbs knotvector

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class Knotvector  : public DRT::ParObject
{
public:
  //! @name Construction, destruction and copying

  /*! \brief standard constructor

      \param  dim       (i)
              defining the dimension of the knotspan
      \param  npatches  (i)
              defining the number of patches

   */
  Knotvector(int                  dim       ,
	     int                  npatches  );

  
  /*!
    \brief empty constructor
   */
  Knotvector();


  //! \brief destructor
  virtual ~Knotvector();

  /*!
      \brief copy constructor
      
      \param old const Knotvector (i)         
  */
  Knotvector(const DRT::NURBS::Knotvector & old);

  //@}

  //! @name Access methods for shapefunction 
  //!       evaluation

  /*!
      \brief convert an element gid to its 
             corresponding triple knot index
      
      \param gid (i)
             given a global element id, this
	     method will return the base indices
	     of the elements 'local' knotspans

      \param npatch (o)	     
             The number of the patch this element 
             belongs to.

      \param loc_cart_id (o)
             patch local cartesian base indices 
             for knotvectors, elements and control points

      \note This method heavily relies on the 
            cartesian structure of the knotspan.
            We assume a element/control point numbering
	    as follows:
	    
      \verbatim
      
            patch local cartesian element index: 
	    
	             (num_u,num_v,num_w)

	    gid = patchoffset+num_u+num_v*nele                 (2d)

	    gid = patchoffset+num_u+num_v*nele+num_w*nele*mele (3d)
	    
      \endverbatim

      \note Here, all elements in a patch are assumed to be
            numbered consecutively and patchoffset is the nuber
            of all elements from previous patches.
	        
	    Example: Element gid=7, linear element, 2D, npatch 1
	    
                     domain size: nele=6,mele=2	    

                     The numbers in the picture are control
		     point gids

      \verbatim

             v ^
               |
                  +----------------+
              24  | 25    26    27 |  28     29    30    31
                  |                |       
              16  | 17    18    19 |  20     21    22    23
                  |                |       
               8  |  9    10    11 |  12     13    14    15
                  +----------------+       
               0     1     2     3     4      5     6     7  ->
                                                              u
         will return (1,1)

      \endverbatim

      The numbering of the control points, is assumed to be 
      accordingly.

  */
  void ConvertEleGidToKnotIds(
    const int      gid        , 
    int         &  npatch     ,
    vector<int> &  loc_cart_id);

  /*!
    \brief get element knot vectors to a given element id

           This method will be called before any 
           shapefunction evaluation

    \param eleknots (o)
           the element local knotvector to the 
	   given a global element id 
      
    \param gid (i)
           given a global element id 

    \note  This method heavily relies on the 
           cartesian structure of the knotspan
           and the consecutive order of elements/control 
           points in the patches.
	    
	   Example: Element gid=7, linear element, 2D
                    The numbers in the picture are control
		    point gids (see above)

      \verbatim

      knots_[0]

      ||----|----|----|----|----|----|----||


      knots_[1]

      ||----|----|----||

      base index (1,1) will allow to access the local knotspan

      eleknots[0]

            |----|----|----|----|           


      eleknots[1]

       |----|----|----||


      \endverbatim

  */
  void GetEleKnots(
    vector<Epetra_SerialDenseVector> & eleknots,
    int                                gid
    );

  //@}

  //! @name Insert methods

  /*!
    \brief set knots for a patch in one direction

    \param direction (i)
           the direction to which the knotvector
	   corresponds
    \param npatch (i)
           the number of the patch this knotvector
	   belongs to
    \param degree (i)
           the degree of the bspline polynomial 
	   belonging to this knot vector
    \param numknots (i)
           the number knots that will be added here
    \param knotvectortype (i)
           specifies whether we add a periodic or
	   interpolating knot vector
    \param directions_knots (i)
           the knotvector to be inserted

  */
  void SetKnots(
    const int                     & direction        , 
    const int                     & npatch           , 
    const int                     & degree          , 
    const int                     & numknots        , 
    const std::string             & knotvectortype   ,
    Teuchos::RCP<vector<double> >   directions_knots);

  //@}

  //! @name Checks

  /*!
    \brief finish
    
    counting the number of knots added, doing
    consistency checks for size and periodicity

    Calculate offset arrays for patches for 
    access methods.
    
    An unfinished knotvector can not be accessed!
  */
  void FinishKnots();

  //@}

  //! @name Pack/Unpack for io --- implementation of ParObject's
  //!       virtual classes 

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined in the 
  ParObject header file
  */
  virtual int UniqueParObjectId() const { return ParObject_Knotvector;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this class

  */
  virtual void Unpack(const vector<char>& data);

  //@}

  //! @name variable access methods

  /*!
  \brief Return the degree of the nurbs patch

      \param npatch (i)	     
             The number of the patch

  */
  virtual vector<int> ReturnDegree(const int npatch)
    {
      if(!filled_)
      {
	dserror("can't access data. knotvector not completed\n");
      }
      return(degree_[npatch]);
    };

  /*!
  \brief Return the knot sizes of the nurbs patch

      \param npatch (i)	     
             The number of the patch

  */
  virtual vector<int> Return_n_x_m_x_l(const int npatch)
    {
      if(!filled_)
      {
	dserror("can't access data. knotvector not completed\n");
      }
      return(n_x_m_x_l_[npatch]);
    };

  /*!
  \brief Return the element sizes of the nurbs patch

      \param npatch (i)	     
             The number of the patch

  */
  virtual vector<int> Return_nele_x_mele_x_lele(const int npatch)
    {
      if(!filled_)
      {
	dserror("can't access data. knotvector not completed\n");
      }
      return(nele_x_mele_x_lele_[npatch]);
    };


  /*!
  \brief Return the element offsets of nurbs patches

  \return offsets array

  */
  virtual vector<int> ReturnOffsets()
    {
      if(!filled_)
      {
	dserror("can't access data. knotvector not completed\n");
      }
      return(offsets_);
    };


  /*!
  \brief Return the number of patches

  \return number of patches

  */
  virtual int ReturnNP()
    {
      if(!filled_)
      {
	dserror("can't access data. knotvector not completed\n");
      }
      return(npatches_);
    };

  //@}

private:
  //! Knotvector types
  enum KnotvectorType
  {
    knotvector_is_interpolating,
    knotvector_is_periodic     ,
    knotvector_is_not_defined     
  };

  //! dimension 
  int                                    dim_;

  //! number of patches
  int                                    npatches_;

  //! indicates that knots are ready for access
  bool                                   filled_;

  //! degree of bspline-polynomials defined on this knotvector
  vector<vector<int> >                   degree_;

  //! number of knots in each direction
  vector<vector<int> >                   n_x_m_x_l_;

  //! number of elements in each direction
  vector<vector<int> >                   nele_x_mele_x_lele_ ;

  //! are the component closed or open knotvectors?
  vector<vector<KnotvectorType> >        interpolation_;

  //! number of elements in each direction
  vector<int>                            offsets_;

  //! the actual values
  vector<vector<Teuchos::RCP<vector<double> > > >
                                         knot_values_;

};

} // namespace NURBS
  
} // namespace DRT

#endif // KNOTVECTOR_H
#endif // CCADISCRET
