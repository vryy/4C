/*!----------------------------------------------------------------------
\file drt_contact_projector.H

\class CONTACT::Projector

\brief A class to perform projections of nodes onto opposing CElements

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_PROJECTOR_H
#define DRT_CONTACT_PROJECTOR_H

#include "drt_cnode.H"
#include "drt_celement.H"

using namespace std;
using namespace Teuchos;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\class Projector

\brief A class to perform projections of nodes onto opposing CElements


\author popp (popp@lnm.mw.tum.de)
*/

class Projector
{
public:
	
	//! @name Constructors and destructors
	
	/*!
	\brief Standard constructor
	
	Constructs an instance of this class.<br>
  Note that this is \b not a collective call as projections are
  performed in parallel by individual processes.
  
  \param twoD (in): true if problem is 2D, false if 3D
	*/
	explicit Projector(bool twoD);

	/*!
	\brief Destructor

	*/
	virtual ~Projector() {}
	
	//@}
	
	//! @name Access methods
	
	/*!
	\brief Return whether this instance was created for 2D or 3D projections
	
	*/
	inline bool IsTwoDimensional() {return twoD_;}
	
	//@}
	
	//! @name 2D and 3D projection methods
	
	/*!
	\brief Project a node onto a CElement along the node's normal
	
	Used to project a slave side node onto a CElement of the master side
  
  This method will compute the coordinates of a projection of a Node in
  the local coordinate system of a CElement. The projection point will
  not neccesarily fall inside the CElement. However, if the projection
  point is far outside the segment's boundaries, problems with the
  internal nonlinear iteration might occur and a warning is issued when
  convergence can not be achieved in a limited number of iterations.
  
  \param node (in): Slave node to project
  \param ele (in) : Master Celement to project on
  \param xi (out) : Local coordinates of projection on CElement 
  */
  bool Project_NodalNormal(CONTACT::CNode& node,
  												 CONTACT::CElement& ele,
  												 double xi[]);
  
  /*!
  \brief Project a node onto a Celement along the interpolated
         outward normal field of the Celement
    
  Used to project a master side node onto a Celement of the slave side
    
  This method will compute the coordinates of a projection of a Node in
  the local coordinate system of a Celement. The projection point will
  not neccesarily fall inside the Celement. However, if the projection
  point is far outside the segment's boundaries, problems with the
  internal nonlinear iteration might occur and a warning is issued when
  convergence can not be achieved in a limited number of iterations.
    
  \param node (in): Master Node to project
  \param ele (in) : Slave Celement to project on
  \param xi (out) : Local coordinates of projection on Celement
  */
  bool Project_ElementNormal(CONTACT::CNode& node,
  													 CONTACT::CElement& ele,
  													 double xi[]);
 
  /*!
  \brief Project a Gauss point onto a CElement along GP normal
  	
  Used to project a slave side GP onto a CElement of the master side
    
   This method will compute the coordinates of a projection of a Gauss
   point in the local coordinate system of a CElement.
    
   \param gpele (in): Slave CElement containing GP to project
   \param gpeta (in): Local coordinates of GP on gpele
   \param ele (in)  : Master Celement to project on
   \param xi (out)  : Local coordinates of projection on CElement 
   */
   bool Project_GaussPoint(CONTACT::CElement& gpele,
  		 										 const double* gpeta,
    											 CONTACT::CElement& ele,
    											 double xi[]);
    
  /*!
  \brief Evaluate F for nodal normal projection

  */
  double Evaluate_F_NodalNormal(CONTACT::CNode& node,
  															CONTACT::CElement& ele,
  															const double* eta);  
  
  /*!
  \brief Evaluate gradF for nodal normal projection

  */
  double Evaluate_gradF_NodalNormal(CONTACT::CNode& node,
  																	CONTACT::CElement& ele,
  																	const double* eta);

  /*!
  \brief Evaluate F for element normal projection

  */
  double Evaluate_F_ElementNormal(CONTACT::CNode& node,
  																CONTACT::CElement& ele,
  																const double* eta);  
  
  /*!
  \brief Evaluate gradF for element normal projection

  */
  double Evaluate_gradF_ElementNormal(CONTACT::CNode& node,
  																		CONTACT::CElement& ele,
  																		const double* eta);
  
  /*!
  \brief Evaluate F for Gauss point projection

  */
  double Evaluate_F_GaussPoint(const double* gpx,
  														 const double* gpn,
   														 CONTACT::CElement& ele,
   														 const double* eta);  
    
  /*!
  \brief Evaluate gradF for Gauss point projection

  */
  double Evaluate_gradF_GaussPoint(const double* gpn,
			 														 CONTACT::CElement& ele,
			 														 const double* eta);  
  
  //@}
  		
protected:
		
		// don't want = operator and cctor
		Projector operator = (const Projector& old);
		Projector(const Projector& old);
	
		bool				twoD_;					// problem dimension, true if 2D

}; // class Projector
}	// namespace CONTACT


#endif  // #ifndef DRT_CONTACT_PROJECTOR_H
#endif  // #ifdef CCADISCRET
