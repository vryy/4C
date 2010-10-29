/*!----------------------------------------------------------------------
\file mortar_analytical.H

\brief A set of analytical solutions for convergence analysis of mortar methods

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_ANALYTICAL_H
#define MORTAR_ANALYTICAL_H

#include "../linalg/linalg_fixedsizematrix.H"

//------------------------------------------------------------------------
// available analytical solutions
//------------------------------------------------------------------------
// 0 = no output of error norms
// 1 = 3D beam bending example
// 2 = 3D pressurized hollow sphere example
//------------------------------------------------------------------------
#define MORTARANALYTICALSOL 0

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms
*/
namespace MORTAR
{
  /*!
  \brief Analytical solutions for 3D elasticity problems

  \param pos (in)       : position where analytical solution is evaluated
  \param uanalyt (out)  : analytical displacement solution
  \param epsanalyt (out): analytical strain solution
  */
  inline void AnalyticalSolutions3D(const LINALG::Matrix<3,1>& pos,
  		                       LINALG::Matrix<3,1>& uanalyt,
  		                       LINALG::Matrix<6,1>& epsanalyt,
  		                       LINALG::Matrix<3,3>& derivanalyt)
	{
		//----------------------------------------------------------------------
		// 3D beam bending solution
		//----------------------------------------------------------------------
		if (MORTARANALYTICALSOL==1)
		{
			// displacements
			uanalyt(0,0) = (pos(0,0)*pos(1,0)) / 10.0;
			uanalyt(1,0) = (-pos(0,0)*pos(0,0)-0.3*pos(1,0)*pos(1,0)+0.3*pos(2,0)*pos(2,0)) / 20.0;
			uanalyt(2,0) = (-0.3*pos(1,0)*pos(2,0)) / 10.0;

			// strains
			epsanalyt(0,0) =  pos(1,0) / 10.0;
			epsanalyt(1,0) = -0.3*pos(1,0) / 10.0;
			epsanalyt(2,0) = -0.3*pos(1,0) / 10.0;
			epsanalyt(3,0) =  0.0;
			epsanalyt(4,0) =  0.0;
			epsanalyt(5,0) =  0.0;

			// displacement derivatives
			derivanalyt(0,0) =  pos(1,0) / 10.0;
			derivanalyt(0,1) =  pos(0,0) / 10.0;
			derivanalyt(0,2) =  0.0;
			derivanalyt(1,0) = -pos(0,0) / 10.0;
			derivanalyt(1,1) = -0.3*pos(1,0) / 10.0;
			derivanalyt(1,2) =  0.3*pos(2,0) / 10.0;
			derivanalyt(2,0) =  0.0;
			derivanalyt(2,1) = -0.3*pos(2,0) / 10.0;
			derivanalyt(2,2) = -0.3*pos(1,0) / 10.0;
		}

		//----------------------------------------------------------------------
		// 3D pressurized sphere solution
		//----------------------------------------------------------------------
		else if (MORTARANALYTICALSOL==2)
		{
			dserror("ERROR: 3D pressurized sphere solution not yet implemented");
		}

		//----------------------------------------------------------------------
		// Other cases
		//----------------------------------------------------------------------
		else
		{
			dserror("ERROR: Other 3D analytical solutions not yet implemented");
		}

		return;
	}

} // namespace MORTAR

#endif  // #ifndef MORTAR_ANALYTICAL_H
#endif  // #ifdef CCADISCRET
