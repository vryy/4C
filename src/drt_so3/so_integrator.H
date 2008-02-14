/*!----------------------------------------------------------------------
\file so_tet10.H

\class DRT::ELEMENTS::So_tet10

\brief A C++ version of the 10-node tet solid element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by: Alexander Volf
			alexander.volf@mytum.de     
</pre>

*----------------------------------------------------------------------*/
#if defined D_SOH8 || D_SOTET
#ifdef CCADISCRET

#ifndef SO_INTEGRATOR_H
#define SO_INTEGRATOR_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"

#include "../drt_mat/material.H"

#define NUMNOD_SOTET10      10      // number of nodes
#define NODDOF_SOTET10       3      // number of dofs per node
#define NUMDOF_SOTET10      30      // total dofs per element
#define NUMSTR_SOTET10       6      // number of strains/stresses 	/****/
#define NUMGPT_SOTET10       4      // total gauss points per element	/****/
#define NUMDIM_SOTET10       3      // number of dimensions 			/****/
#define NUMCOORD_SOTET10	 4 		// number of shape function cooordinates (ksi1-ksi4) 
#define NUMNOD_SOTET10_FACE  6		// number of nodes on a TET10 face (which is a TRI6) 
#define NUMGPT_SOTET10_FACE	 3		// number of GP    on a TET10 face (which is a TRI6) 

// forward declarations
struct _MATERIAL;
class SoDisp;
/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

namespace ELEMENTS
{
// forward declarations	
	class SoDisp;
/*!
\class DRT::ELEMENTS::So_integrator

\brief A general class agregating all necessery data for discrete integration

<h3>Purpose</h3>
The class offers a standard interface for discrete integration, mainly the agregated
shape functions and shape function derivatives, it is a purely general class and 
does NOT own a constructor. Constructors exist only for derived classes. It has however  
a virtual destructor, that should be the same for all derived classes. Derived classes are not
supposed to add extra attributes, they may however add extra methods.

pointer to a vector<> of shape functions stored in Epetra_Vectors			
evaluated at each gp, for each node. It may look like this:



\f[
   \text{shapefct[gp]} =
   \begin{bmatrix}
   		N_1[gp] & N_2[gp] & \hdots & N_X[gp]
   \end{bmatrix}
\f]

pointer to a vector<> of shape function derivatives stored in
Epetra_Matrices evaluated at each gp, for \f$\xi_1-\xi_{4}\f$ it may look like this:	
													
													
\f[
	\text{deriv\_gp} = 
	\begin{bmatrix}								    
		N_{1,xsi1} &  N_{1,xsi2} &  N_{1,xsi3} &  N_{1,xsi4} \\
		N_{2,xsi1} &  N_{2,xsi2} &  N_{2,xsi3} &  N_{2,xsi4} \\ 
		\vdots & \vdots & \vdots & \vdots \\
		N_{X,xsi1} &  N_{X,xsi2} & N_{X,xsi3} & N_{X,xsi4}
	\end{bmatrix}
\f]

NOTE: X stands for the number of NODES (and shape functions) in a element


\author a.vlf
*/

 /*---------------------------------------------------------------------------------*
  | class Discrete_integrator                                           vlf 09/07   |
  |                                                                                 |
  | a general class agregating all necessery data for discrete integration	        |
  |                                                                                 |
  |  pointer to vector<> of shape function stored in Epetra_Vectors	                |
  |  evaluated at each gp , for each node  it looks like this:                      |
  |  shapefct[gp] = [N_1(gp)  N_2(gp)  -- N_X(gp)]                                  |
  |                                                                                 |
  |  pointer to vector<> of shape function derivatives stored in                    |
  |	 Epetra_Matrices evaluated at each gp, for xsi1-4                               |
  |  it looks like this:														    |
  |                 [  N1_,xsi1   N1_,xsi2    N1_,xsi3   N1_,xsi4 ]	                |
  |  deriv_gp[gp] = [  N2_,xsi1   N2_,xsi2    N2_,xsi3   N2_,xsi4 ] (at a given gp)	|
  | 		        [     |	         |           |         |      ]	                |
  |                 [ NX_,xsi1    NX_,xsi2    NX_,xsi3   NX_,xsi4 ]   				|
  |  where X stands for the number of NODES (and shape functions) in a element		|
  |                                                                                 |
  *---------------------------------------------------------------------------------*/ 
  class So_integrator
  {
	public: 
		vector<Epetra_SerialDenseVector> shapefct_gp;  
      	vector<Epetra_SerialDenseMatrix> deriv_gp;
      	Epetra_SerialDenseVector weights;
      	int num_gp;
      	int num_nodes;
      	int num_coords;
      	virtual ~So_integrator();
  };
 /*---------------------------------------------------------------------------------*
  | class Tet10_integrator_4point                                       vlf 09/07   |
  |                                                                                 |
  | a special class for quadratic 10-node tetrahedra using a 4 point quadrature   	|
  |	it is exaxt for linear and quadratic polynomials                                |
  |                                                                                 |
  *---------------------------------------------------------------------------------*/
  class Integrator_tet10_4point : public So_integrator
  {
  	public:
  		Integrator_tet10_4point(void);
  };
 
 /*---------------------------------------------------------------------------------*
  | class Tet_integrator_14point                                        vlf 09/07   |
  |                                                                                 |
  | a special class for quadratic 10-node tetrahedra using a 14 point quadrature   	|
  |	it is exaxt for 0-4th order polynomials                                         |
  |                                                                                 |
  *---------------------------------------------------------------------------------*/
  
  class Integrator_tet10_14point : public So_integrator
  {
  	public:
  		Integrator_tet10_14point(void);
  };
  
  /*---------------------------------------------------------------------------------*
  | class Tet_integrator_1point											vlf 09/07   |
  |																					|
  | a special class for quadratic 4-node tetrahedra using a 4 point quadrature   	|
  |	it is exaxt for linear and quadratic polynomials								|
  |																					|
  *---------------------------------------------------------------------------------*/
  class Integrator_tet4_1point : public So_integrator
  {
  	public:
  		Integrator_tet4_1point(void);
  };

  class Integrator_tri3_1point : public So_integrator
  {
  	public:
  		Integrator_tri3_1point(void);
  };
  
  class Integrator_tet4_4point : public So_integrator
  {
  	public:
  		Integrator_tet4_4point(void);
  };
} // namespace ELEMENTS
} // namespace DRT

long double det_volf(Epetra_SerialDenseMatrix& in_matrix);


#endif  // #ifndef SO_INTEGRATOR
#endif  // #ifdef  CCADISCRET
#endif  // #if D_SOTET
