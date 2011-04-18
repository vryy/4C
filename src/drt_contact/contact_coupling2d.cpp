/*!----------------------------------------------------------------------
\file contact_coupling2d.cpp
\brief Classes for mortar contact coupling in 2D.

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
#include "../drt_mortar/mortar_element.H"
#include "../drt_lib/drt_globalproblem.H"


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
CONTACT::CoCoupling2d::CoCoupling2d(const INPAR::MORTAR::ShapeFcn shapefcn,
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
bool CONTACT::CoCoupling2d::IntegrateOverlap()
{
  // explicitely defined shapefunction type needed
  if( shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateOverlap called without specific shape function defined!");
  
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj_ entries integrate the Mortar  */
  /* matrices D and M and the weighted gap function g~ on the overlap   */
  /* of the current sl / ma pair.                                       */
  /**********************************************************************/

  // no integration if no overlap
  if (!overlap_) return false;

  //local working copies of input variables
  double sxia = xiproj_[0];
  double sxib = xiproj_[1];
  double mxia = xiproj_[2];
  double mxib = xiproj_[3];

  // create a CONTACT integrator instance with correct NumGP and Dim
  CONTACT::CoIntegrator integrator(shapefcn_,SlaveElement().Shape());

  // do the overlap integration (integrate and linearize both M and gap and wear)
  int nrow = SlaveElement().NumNode();
  int ncol = MasterElement().NumNode();
  RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),nrow*Dim()));
  RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),ncol*Dim()));
  RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));
  RCP<Epetra_SerialDenseVector> wseg = Teuchos::null;
  if((DRT::Problem::Instance()->MeshtyingAndContactParams()).get<double>("WEARCOEFF")>0.0)
    wseg = rcp(new Epetra_SerialDenseVector(nrow));

  integrator.IntegrateDerivSegment2D(SlaveElement(),sxia,sxib,MasterElement(),mxia,mxib,dseg,mseg,gseg,wseg);

  // do the two assemblies into the slave nodes
  integrator.AssembleD(Comm(),SlaveElement(),*dseg);
  integrator.AssembleM(Comm(),SlaveElement(),MasterElement(),*mseg);
  
  // also do assembly of weighted gap vector
  integrator.AssembleG(Comm(),SlaveElement(),*gseg);
  
  // assemble wear
  if((DRT::Problem::Instance()->MeshtyingAndContactParams()).get<double>("WEARCOEFF")>0.0)
    integrator.AssembleWear(Comm(),SlaveElement(),*wseg);

  return true;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::CoCoupling2dManager::CoCoupling2dManager(DRT::Discretization& idiscret, int dim,
                                                  MORTAR::MortarElement* sele,
                                                  vector<MORTAR::MortarElement*> mele) :
shapefcn_(INPAR::MORTAR::shape_undefined),
idiscret_(idiscret),
dim_(dim),
sele_(sele),
mele_(mele)
{
  // evaluate coupling
  EvaluateCoupling();

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
CONTACT::CoCoupling2dManager::CoCoupling2dManager(const INPAR::MORTAR::ShapeFcn shapefcn,
                                                 DRT::Discretization& idiscret, int dim,
                                                 MORTAR::MortarElement* sele,
                                                 vector<MORTAR::MortarElement*> mele) :
shapefcn_(shapefcn),
idiscret_(idiscret),
dim_(dim),
sele_(sele),
mele_(mele)
{
  // evaluate coupling
  EvaluateCoupling();

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate coupling pairs                                   popp 03/09|
 *----------------------------------------------------------------------*/
bool CONTACT::CoCoupling2dManager::EvaluateCoupling()
{
  // loop over all master elements associated with this slave element
  for (int m=0;m<(int)MasterElements().size();++m)
  {
    // create CoCoupling2d object and push back
    Coupling().push_back(rcp(new CoCoupling2d(shapefcn_,idiscret_,dim_,SlaveElement(),MasterElement(m))));

    // project the element pair
    Coupling()[m]->Project();

    // check for element overlap
    Coupling()[m]->DetectOverlap();

    // integrate the element overlap
    Coupling()[m]->IntegrateOverlap();
  }

  return true;
}

#endif //#ifdef CCADISCRET
