/*!----------------------------------------------------------------------
\file drt_contact_coupling2d.cpp
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

#include "drt_contact_coupling2d.H"
#include "drt_contact_projector.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
CONTACT::Coupling2d::Coupling2d(DRT::Discretization& idiscret, int dim,
                            CONTACT::CElement& sele, CONTACT::CElement& mele,
                            Epetra_SerialDenseMatrix& csegs) :
idiscret_(idiscret),
dim_(dim),
sele_(sele),
mele_(mele),
contactsegs_(csegs)
{
  // *********************************************************************
  // the two-dimensional case
  // *********************************************************************
  // prepare overlap integration
  vector<bool> hasproj(4);
  vector<double> xiproj(4);
  bool overlap = false;

  // project the element pair
  Project(hasproj,xiproj);

  // check for element overlap
  overlap = DetectOverlap(hasproj,xiproj);

  // integrate the element overlap
  if (overlap) IntegrateOverlap(xiproj);

  return;
}


/*----------------------------------------------------------------------*
 |  Project slave / master element pair (public)              popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling2d::Project(vector<bool>& hasproj,
                                vector<double>& xiproj)
{
  // initialize projection status
  hasproj[0] = false;   // slave 0 end node
  hasproj[1] = false;   // slave 1 end node
  hasproj[2] = false;   // master 0 end node
  hasproj[3] = false;   // master 1 end node

  // get slave and master element nodes
  DRT::Node** mysnodes = sele_.Nodes();
  if (!mysnodes)
    dserror("ERROR: IntegrateOverlap: Null pointer for mysnodes!");
  DRT::Node** mymnodes = mele_.Nodes();
  if (!mymnodes)
    dserror("ERROR: IntegrateOverlap: Null pointer for mymnodes!");

  // create a projector instance of problem dimension Dim()
  CONTACT::Projector projector(Dim());

  // project slave nodes onto master element
  for (int i=0;i<sele_.NumNode();++i)
  {
    CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(mysnodes[i]);
    double xi[2] = {0.0, 0.0};
    projector.ProjectNodalNormal(*snode,mele_,xi);

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!
    if ((-1.0-CONTACTPROJTOL<=xi[0]) && (xi[0]<=1.0+CONTACTPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i<2)
      {
        hasproj[i]=true;
        xiproj[i]=xi[0];
      }
      // nevertheless we need the inner node projection status later (weighted gap)
      snode->HasProj()=true;
    }
  }

  // project master nodes onto slave element
  for (int i=0;i<2;++i)
  {
    CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(mymnodes[i]);
    double xi[2] = {0.0, 0.0};
    projector.ProjectElementNormal(*mnode,sele_,xi);

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!!!
    if ((-1.0-CONTACTPROJTOL<=xi[0]) && (xi[0]<=1.0+CONTACTPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i<2)
      {
        hasproj[i+2]=true;
        xiproj[i+2]=xi[0];
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Detect overlap of slave / master pair (public)            popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling2d::DetectOverlap(vector<bool>& hasproj,
                                      vector<double>& xiproj)
{
  /**********************************************************************/
  /* OVERLAP CASES                                                      */
  /* Depending on mxi and sxi overlap will be decided!                  */
  /* Even for 3noded CElements only the two end nodes matter in 2D!     */
  /* There are several cases how the 2 elements can overlap. Handle all */
  /* of them, including the ones that they don't overlap at all!        */
  /**********************************************************************/

  // For the non-overlapping cases, the possibility of an identical local
  // node numbering direction for both sides is taken into account!!
  // (this can happen, when elements far from each other are projected,
  // which actually should be impossible due to the search radius
  // condition in the potential contact pair search above!
  // But you never know...)

  // For the overlapping cases, it is a prerequisite that the two local
  // node numbering directions are opposite!!
  // (this is the case, when the elements are sufficiently near each other,
  // which is ensured by only processing nodes that fulfill the
  // search radius condition above!)

  // CAUTION: The bool output variable in this method is a REAL output
  // variable, determining whether there is an overlap or not!

  // initialize local working variables
  bool overlap = false;
  double sxia = 0.0;
  double sxib = 0.0;
  double mxia = 0.0;
  double mxib = 0.0;

  // local working copies of input variables
  bool s0hasproj = hasproj[0];
  bool s1hasproj = hasproj[1];
  bool m0hasproj = hasproj[2];
  bool m1hasproj = hasproj[3];

  vector<double> sprojxi(2);
  sprojxi[0] = xiproj[0];
  sprojxi[1] = xiproj[1];

  vector<double> mprojxi(2);
  mprojxi[0] = xiproj[2];
  mprojxi[1] = xiproj[3];


  /* CASE 1 (NO OVERLAP):
     no feasible projection found for any of the 4 outer element nodes  */

  if (!s0hasproj && !s1hasproj && !m0hasproj && !m1hasproj)
  {
    //do nothing
  }

  /* CASES 2-5 (NO OVERLAP):
     feasible projection found only for 1 of the 4 outer element nodes
     (this can happen due to the necessary projection tolerance!!!)     */

  else if  (s0hasproj && !s1hasproj && !m0hasproj && !m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=sprojxi[0]) && (sprojxi[0]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  else if  (!s0hasproj && s1hasproj && !m0hasproj && !m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=sprojxi[1]) && (sprojxi[1]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  else if  (!s0hasproj && !s1hasproj && m0hasproj && !m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=mprojxi[0]) && (mprojxi[0]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  else if  (!s0hasproj && !s1hasproj && !m0hasproj && m1hasproj)
  {
    if ((-1.0+CONTACTPROJTOL<=mprojxi[1]) && (mprojxi[1]<=1.0-CONTACTPROJTOL))
    {
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", sele_.Id(), mele_.Id());
    }
  }

  /* CASE 6 (OVERLAP):
     feasible projection found for all 4 outer element nodes
     (this can happen due to the necessary projection tolerance!!!)     */

  else if (s0hasproj && s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;

    // internal case 1 for global CASE 6
    // (equivalent to global CASE 7, slave fully projects onto master)
    if ((sprojxi[0]<1.0) && (sprojxi[1]>-1.0))
    {
      sxia = -1.0;
      sxib = 1.0;
      mxia = sprojxi[1];      // local node numbering always anti-clockwise!!!
      mxib = sprojxi[0];
      //cout << "Problem solved with internal case 1!" << endl;
    }

    // internal case 2 for global CASE 6
    // (equivalent to global CASE 8, master fully projects onto slave)
    else if ((mprojxi[0]<1.0) && (mprojxi[1]>-1.0))
    {
      mxia = -1.0;
      mxib = 1.0;
      sxia = mprojxi[1];      // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
      //cout << "Problem solved with internal case 2!" << endl;
    }

    // internal case 3 for global CASE 6
    // (equivalent to global CASE 9, both nodes no. 0 project successfully)
    else if ((sprojxi[0]<1.0+CONTACTPROJLIM) && (mprojxi[0]<1.0+CONTACTPROJLIM))
    {
      sxia = -1.0;
      sxib = mprojxi[0];      // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
      //cout << "Problem solved with internal case 3!" << endl;
    }

    // internal case 4 for global CASE 6
    // (equivalent to global CASE 10, both nodes no. 1 project successfully)
    else if ((sprojxi[1]>-1.0-CONTACTPROJLIM) && (mprojxi[1]>-1.0-CONTACTPROJLIM))
    {
      sxia = mprojxi[1];
      sxib = 1.0;            // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
      //cout << "Problem solved with internal case 4!" << endl;
    }

    // unknown internal case for global CASE 6
    else
    {
      cout << "CONTACT::Interface::IntegrateOverlap "<< endl << "has detected '4 projections'-case for Sl./Ma. pair "
           << sele_.Id() << "/" << mele_.Id() << endl;
      cout << "SElement Node IDs: " << (sele_.Nodes()[0])->Id() << " " << (sele_.Nodes()[1])->Id() << endl;
      cout << "MElement Node IDs: " << (mele_.Nodes()[0])->Id() << " " << (mele_.Nodes()[1])->Id() << endl;
      cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Unknown overlap case found in global case 6!");
    }
  }

  /* CASES 7-8 (OVERLAP):
     feasible projections found for both nodes of one element, this
      means one of the two elements is projecting fully onto the other!  */

  else if (s0hasproj && s1hasproj && !m0hasproj && !m1hasproj)
  {
    overlap = true;
    sxia = -1.0;
    sxib = 1.0;
    mxia = sprojxi[1];      // local node numbering always anti-clockwise!!!
    mxib = sprojxi[0];
  }

  else if (!s0hasproj && !s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;
    mxia = -1.0;
    mxib = 1.0;
    sxia = mprojxi[1];      // local node numbering always anti-clockwise!!!
    sxib = mprojxi[0];
  }

  /* CASES 9-10 (OVERLAP):
     feasible projections found for one node of each element, due to
     node numbering only identical local node ID pairs possible!        */

  else if (s0hasproj && !s1hasproj && m0hasproj && !m1hasproj)
  {
    // do the two elements really have an overlap?
    if ((sprojxi[0]>-1.0+CONTACTPROJLIM) && (mprojxi[0]>-1.0+CONTACTPROJLIM))
    {
      overlap = true;
      sxia = -1.0;
      sxib = mprojxi[0];      // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
    }
  }

  else if (!s0hasproj && s1hasproj && !m0hasproj && m1hasproj)
  {
    // do the two elements really have an overlap?
    if ((sprojxi[1]<1.0-CONTACTPROJLIM) && (mprojxi[1]<1.0-CONTACTPROJLIM))
    {
      overlap = true;
      sxia = mprojxi[1];
      sxib = 1.0;            // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
    }
  }

  /* CASES 11-14 (OVERLAP):
     feasible projections found for 3 out of the total 4 nodes,
     this can either lead to cases 7/8 or 9/10!                         */
  else if (s0hasproj && s1hasproj && m0hasproj && !m1hasproj)
  {
    overlap = true;
    // equivalent to global case 7
    if (mprojxi[0]>1.0)
    {
      sxia = -1.0;
      sxib = 1.0;
      mxia = sprojxi[1];    // local node numbering always anti-clockwise!!!
      mxib = sprojxi[0];
    }
    // equivalent to global case 9
    else
    {
      sxia = -1.0;
      sxib = mprojxi[0];    // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
    }
  }

  else if (s0hasproj && s1hasproj && !m0hasproj && m1hasproj)
  {
    overlap = true;
    // equivalent to global case 7
    if (mprojxi[1]<-1.0)
    {
      sxia = -1.0;
      sxib = 1.0;
      mxia = sprojxi[1];  // local node numbering always anti-clockwise!!!
      mxib = sprojxi[0];
    }
    // equivalent to global case 10
    else
    {
      sxia = mprojxi[1];
      sxib = 1.0;          // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
    }
  }

  else if (s0hasproj && !s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;
    // equivalent to global case 8
    if (sprojxi[0]>1.0)
    {
      mxia = -1.0;
      mxib = 1.0;
      sxia = mprojxi[1];      // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
    }
    // equivalent to global case 9
    else
    {
      sxia = -1.0;
      sxib = mprojxi[0];    // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
    }
  }

  else if (!s0hasproj && s1hasproj && m0hasproj && m1hasproj)
  {
    overlap = true;
    // equivalent to global case 8
    if (sprojxi[1]<-1.0)
    {
      mxia = -1.0;
      mxib = 1.0;
      sxia = mprojxi[1];  // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
    }
    // equivalent to global case 10
    else
    {
      sxia = mprojxi[1];
      sxib = 1.0;          // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
    }
  }

  /* CASE DEFAULT: unknown overlap case                                  */
  else
  {
    cout << "SElement: " << sele_.NodeIds()[0] << " " << sele_.NodeIds()[1] << endl;
    cout << "MElement: " << mele_.NodeIds()[0] << " " << mele_.NodeIds()[1] << endl;
    cout << "s0: " << s0hasproj << " s1: " << s1hasproj << endl;
    cout << "m0: " << m0hasproj << " m1: " << m1hasproj << endl;
    dserror("ERROR: IntegrateOverlap: Unknown overlap case found!");
  }

  // check for 1:1 node projections and for infeasible limits
  if ((sxia<-1.0) || (sxib>1.0) || (mxia<-1.0) || (mxib>1.0))
  {
    if (abs(sxia+1.0)<CONTACTPROJLIM) sxia=-1.0;
    if (abs(sxib-1.0)<CONTACTPROJLIM) sxib= 1.0;
    if (abs(mxia+1.0)<CONTACTPROJLIM) mxia=-1.0;
    if (abs(mxib-1.0)<CONTACTPROJLIM) mxib= 1.0;

    if ((sxia<-1.0) || (sxib>1.0) || (mxia<-1.0) || (mxib>1.0))
    {
      cout << "Slave: " << sxia << " " << sxib << endl;
      cout << "Master: " << mxia << " " << mxib << endl;
      dserror("ERROR: IntegrateOverlap: Determined infeasible limits!");
    }
  }

  // update integration limits in xiproj
  xiproj[0]=sxia;
  xiproj[1]=sxib;
  xiproj[2]=mxia;
  xiproj[3]=mxib;

  // prepare gmsh visualization
#ifdef DEBUG
  if (overlap)
  {
    double sxialoc[2] = {sxia, 0.0};
    double sxibloc[2] = {sxib, 0.0};
    double mxialoc[2] = {mxia, 0.0};
    double mxibloc[2] = {mxib, 0.0};

    double sxiaglob[3] = {0.0, 0.0, 0.0};
    double sxibglob[3] = {0.0, 0.0, 0.0};
    double mxiaglob[3] = {0.0, 0.0, 0.0};
    double mxibglob[3] = {0.0, 0.0, 0.0};

    sele_.LocalToGlobal(sxialoc,sxiaglob,0);
    sele_.LocalToGlobal(sxibloc,sxibglob,0);
    mele_.LocalToGlobal(mxialoc,mxiaglob,0);
    mele_.LocalToGlobal(mxibloc,mxibglob,0);

    Epetra_SerialDenseMatrix& segs = CSegs();
    segs.Reshape(segs.M()+1,12);
    segs(segs.M()-1,0)  = sxiaglob[0];
    segs(segs.M()-1,1)  = sxiaglob[1];
    segs(segs.M()-1,2)  = sxiaglob[2];
    segs(segs.M()-1,3)  = mxibglob[0];
    segs(segs.M()-1,4)  = mxibglob[1];
    segs(segs.M()-1,5)  = mxibglob[2];
    segs(segs.M()-1,6)  = mxiaglob[0];
    segs(segs.M()-1,7)  = mxiaglob[1];
    segs(segs.M()-1,8)  = mxiaglob[2];
    segs(segs.M()-1,9)  = sxibglob[0];
    segs(segs.M()-1,10) = sxibglob[1];
    segs(segs.M()-1,11) = sxibglob[2];
  }
#endif // #ifdef DEBUG

  return overlap;
}

/*----------------------------------------------------------------------*
 |  Integrate slave / master overlap (public)                 popp 04/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Coupling2d::IntegrateOverlap(vector<double>& xiproj)
{
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj entries integrate the Mortar   */
  /* matrix M and the weighted gap function g~ on the overlap of the    */
  /* current slave / master CElement pair                               */
  /**********************************************************************/

  //local working copies of input variables
  double sxia = xiproj[0];
  double sxib = xiproj[1];
  double mxia = xiproj[2];
  double mxib = xiproj[3];

  // create an integrator instance with correct NumGP and Dim
  CONTACT::Integrator integrator(sele_.Shape());

  // do the overlap integration (integrate and linearize both M and gap)
  int nrow = sele_.NumNode();
  int ncol = mele_.NumNode();
  RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),nrow*Dim()));
  RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*Dim(),ncol*Dim()));
  RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));
  integrator.IntegrateDerivSegment2D(sele_,sxia,sxib,mele_,mxia,mxib,dseg,mseg,gseg);

  // do the two assemblies into the slave nodes
#ifdef CONTACTONEMORTARLOOP
  integrator.AssembleD(Comm(),sele_,*dseg);
#endif // #ifdef CONTACTONEMORTARLOOP
  integrator.AssembleM(Comm(),sele_,mele_,*mseg);
  integrator.AssembleG(Comm(),sele_,*gseg);

  /*----------------------------------------------------------------------
  // check for the modification of the M matrix for curved interfaces
  // (based on the paper by M. Puso / B. Wohlmuth, IJNME, 2005)
  // (note: we assume that the modification is not useful for mortar CONTACT,
  // but only for mortar MESH TYING as published!)
  //----------------------------------------------------------------------
  bool modification = false;

  // conditions for modification
  // (1) linear shape functions for the slave elements
  // (2) 2D problems (3D unknown / unpublished ??)
  // (3) curved interface, n1!=n2
  // (4) use of dual shape functions for LM (always true in our case !!)
  if (sele_.Shape()==DRT::Element::line2)
  {
    if (integrator.Dim()==2)
    {
      CNode* snode0 = static_cast<CNode*>(sele_.Nodes()[0]);
      CNode* snode1 = static_cast<CNode*>(sele_.Nodes()[1]);

      const double* n0 = snode0->n();
      const double* n1 = snode1->n();
      double delta = (n0[0]-n1[0])*(n0[0]-n1[0]) + (n0[1]-n1[1])*(n0[1]-n1[1]);

      if (delta>1.0e-8)
        modification = true;
    }
  }

  // integrate and assemble the modification, if necessary

  if (modification)
  {
    RCP<Epetra_SerialDenseMatrix> mmodseg = integrator.IntegrateMmod2D(sele_,sxia,sxib,mele_,mxia,mxib);
    integrator.AssembleMmod(Comm(),sele_,mele_,*mmodseg);
  }
  //--------------------------------------------------------------------*/

  return true;
}

#endif //#ifdef CCADISCRET
