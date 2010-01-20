/*!----------------------------------------------------------------------
\file contact_coupling2d.cpp
\brief A class for mortar coupling of ONE slave element and ONE master
       element of a contact interface in 2D.

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

#include "contact_coupling2d.H"
#include "contact_integrator.H"
#include "../drt_mortar/mortar_defines.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::CoCoupling2d::CoCoupling2d(DRT::Discretization& idiscret, int dim,
                                    MORTAR::MortarElement& sele,
                                    MORTAR::MortarElement& mele) :
MORTAR::Coupling2d(idiscret,dim,sele,mele)
{
  // empty constructor
  
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
CONTACT::CoCoupling2d::CoCoupling2d(const MORTAR::MortarInterface::ShapeFcnType shapefcn,
                                    DRT::Discretization& idiscret, int dim,
                                    MORTAR::MortarElement& sele,
                                    MORTAR::MortarElement& mele) :
MORTAR::Coupling2d(shapefcn,idiscret,dim,sele,mele)
{
  // empty constructor
  
  return;
}

/*----------------------------------------------------------------------*
 |  Integrate slave / master overlap (public)                 popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoCoupling2d::IntegrateOverlap(vector<double>& xiproj)
{
  // explicitely defined shapefunction type needed
  if( shapefcn_ == MORTAR::MortarInterface::Undefined)
    dserror("ERROR: IntegrateOverlap called without specific shape function defined!");
  
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj entries integrate the Mortar   */
  /* matrices D and M and the weighted gap function g~ on the overlap   */
  /* of the current sl / ma pair.                                       */
  /**********************************************************************/

  //local working copies of input variables
  double sxia = xiproj[0];
  double sxib = xiproj[1];
  double mxia = xiproj[2];
  double mxib = xiproj[3];

  // create a CONTACT integrator instance with correct NumGP and Dim
  CONTACT::CoIntegrator integrator(shapefcn_,sele_.Shape());

  // do the overlap integration (integrate and linearize both M and gap)
  int nrow = sele_.NumNode();
  int ncol = mele_.NumNode();
  RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),nrow*Dim()));
  RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),ncol*Dim()));
  RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));
  integrator.IntegrateDerivSegment2D(sele_,sxia,sxib,mele_,mxia,mxib,dseg,mseg,gseg);

  // do the two assemblies into the slave nodes
#ifdef MORTARONELOOP
  integrator.AssembleD(Comm(),sele_,*dseg);
#endif // #ifdef MORTARONELOOP
  integrator.AssembleM(Comm(),sele_,mele_,*mseg);
  
  // also do assembly of weighted gap vector
  integrator.AssembleG(Comm(),sele_,*gseg);

  return true;
}

#endif //#ifdef CCADISCRET
