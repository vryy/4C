/*!----------------------------------------------------------------------
\file mortar_coupling2d.cpp
\brief Classes for mortar coupling in 2D.

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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include "mortar_coupling2d.H"
#include "mortar_node.H"
#include "mortar_element.H"
#include "mortar_projector.H"
#include "mortar_integrator.H"
#include "mortar_defines.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Coupling2d::Coupling2d(DRT::Discretization& idiscret,
                               int dim, bool quad, INPAR::MORTAR::LagMultQuad lmtype,
                               MORTAR::MortarElement& sele,
                               MORTAR::MortarElement& mele) :
shapefcn_(INPAR::MORTAR::shape_undefined),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
lmtype_(lmtype),
sele_(sele),
mele_(mele),
overlap_(false)
{
  // initialize variables
  hasproj_.resize(4);
  xiproj_.resize(4);

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 06/09|
 *----------------------------------------------------------------------*/
MORTAR::Coupling2d::Coupling2d(const INPAR::MORTAR::ShapeFcn shapefcn,
                               DRT::Discretization& idiscret,
                               int dim, bool quad, INPAR::MORTAR::LagMultQuad lmtype,
                               MORTAR::MortarElement& sele,
                               MORTAR::MortarElement& mele) :
shapefcn_(shapefcn),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
lmtype_(lmtype),
sele_(sele),
mele_(mele),
overlap_(false)
{
  // initialize variables
  hasproj_.resize(4);
  xiproj_.resize(4);
  
  return;
}

/*----------------------------------------------------------------------*
 |  get communicator  (public)                                popp 06/09|
 *----------------------------------------------------------------------*/
const Epetra_Comm& MORTAR::Coupling2d::Comm() const
{
  return idiscret_.Comm();
}

/*----------------------------------------------------------------------*
 |  Rough check if elements are near (with normals)           popp 11/10|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling2d::RoughCheckOrient()
{
  // we first need the master element center
  double loccenter[2] = {0.0, 0.0};

  // compute the unit normal vector at the slave element center
  double nsc[3] = {0.0, 0.0, 0.0};
  SlaveElement().ComputeUnitNormalAtXi(loccenter,nsc);

  // compute the unit normal vector at the master element center
  double nmc[3] = {0.0, 0.0, 0.0};
  MasterElement().ComputeUnitNormalAtXi(loccenter,nmc);

  // check orientation of the two normals
  double dot = nsc[0]*nmc[0]+nsc[1]*nmc[1]+nsc[2]*nmc[2];
  if (dot < -1.0e-12) return true;
  else                return false;
}

/*----------------------------------------------------------------------*
 |  Project slave / master element pair (public)              popp 04/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling2d::Project()
{
  // initialize projection status
  hasproj_[0] = false;   // slave 0 end node
  hasproj_[1] = false;   // slave 1 end node
  hasproj_[2] = false;   // master 0 end node
  hasproj_[3] = false;   // master 1 end node

  // rough check of orientation of element centers
  // if slave and master element center normals form an
  // angle > 90° the pair will not be considered further
  bool orient = RoughCheckOrient();
  if (!orient) return false;

  // get slave and master element nodes
  DRT::Node** mysnodes = SlaveElement().Nodes();
  if (!mysnodes) dserror("ERROR: IntegrateOverlap: Null pointer for mysnodes!");
  DRT::Node** mymnodes = MasterElement().Nodes();
  if (!mymnodes) dserror("ERROR: IntegrateOverlap: Null pointer for mymnodes!");

  // create a projector instance of problem dimension Dim()
  MORTAR::MortarProjector projector(Dim());

  // project slave nodes onto master element
  for (int i=0;i<SlaveElement().NumNode();++i)
  {
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(mysnodes[i]);
    double xi[2] = {0.0, 0.0};
    projector.ProjectNodalNormal(*snode,MasterElement(),xi);

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!
    if ((-1.0-MORTARPROJTOL<=xi[0]) && (xi[0]<=1.0+MORTARPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i<2)
      {
        hasproj_[i]=true;
        xiproj_[i]=xi[0];
      }
      // nevertheless we need the inner node projection status later (weighted gap)
      snode->HasProj()=true;
    }
  }

  // project master nodes onto slave element
  for (int i=0;i<2;++i)
  {
    MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mymnodes[i]);
    double xi[2] = {0.0, 0.0};
    projector.ProjectElementNormal(*mnode,SlaveElement(),xi);

    // save projection if it is feasible
    // we need an expanded feasible domain in order to check pathological
    // cases due to round-off error and iteration tolerances later!!!
    if ((-1.0-MORTARPROJTOL<=xi[0]) && (xi[0]<=1.0+MORTARPROJTOL))
    {
      // for element overlap only the outer nodes are of interest
      if (i<2)
      {
        hasproj_[i+2]=true;
        xiproj_[i+2]=xi[0];
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Detect overlap of slave / master pair (public)            popp 04/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling2d::DetectOverlap()
{
  /**********************************************************************/
  /* OVERLAP CASES                                                      */
  /* Depending on mxi and sxi overlap will be decided!                  */
  /* Even for 3noded elements only the two end nodes matter in 2D!      */
  /* There are several cases how the 2 elements can overlap. Handle all */
  /* of them, including the ones that they don't overlap at all!        */
  /**********************************************************************/

  // For the non-overlapping cases, the possibility of an identical local
  // node numbering direction for both sides is taken into account!!
  // (this can happen, when elements far from each other are projected,
  // which actually should be impossible due to the search radius
  // condition in the potential coupling pair search above!
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
  bool s0hasproj = hasproj_[0];
  bool s1hasproj = hasproj_[1];
  bool m0hasproj = hasproj_[2];
  bool m1hasproj = hasproj_[3];

  std::vector<double> sprojxi(2);
  sprojxi[0] = xiproj_[0];
  sprojxi[1] = xiproj_[1];

  std::vector<double> mprojxi(2);
  mprojxi[0] = xiproj_[2];
  mprojxi[1] = xiproj_[3];


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
    if ((-1.0+MORTARPROJTOL<=sprojxi[0]) && (sprojxi[0]<=1.0-MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (SlaveElement().Nodes()[0])->Id() << " " << (SlaveElement().Nodes()[1])->Id() << endl;
      std::cout << "MElement Node IDs: " << (MasterElement().Nodes()[0])->Id() << " " << (MasterElement().Nodes()[1])->Id() << endl;
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", SlaveElement().Id(), MasterElement().Id());
    }
  }

  else if  (!s0hasproj && s1hasproj && !m0hasproj && !m1hasproj)
  {
    if ((-1.0+MORTARPROJTOL<=sprojxi[1]) && (sprojxi[1]<=1.0-MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (SlaveElement().Nodes()[0])->Id() << " " << (SlaveElement().Nodes()[1])->Id() << endl;
      std::cout << "MElement Node IDs: " << (MasterElement().Nodes()[0])->Id() << " " << (MasterElement().Nodes()[1])->Id() << endl;
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", SlaveElement().Id(), MasterElement().Id());
    }
  }

  else if  (!s0hasproj && !s1hasproj && m0hasproj && !m1hasproj)
  {
    if ((-1.0+MORTARPROJTOL<=mprojxi[0]) && (mprojxi[0]<=1.0-MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (SlaveElement().Nodes()[0])->Id() << " " << (SlaveElement().Nodes()[1])->Id() << endl;
      std::cout << "MElement Node IDs: " << (MasterElement().Nodes()[0])->Id() << " " << (MasterElement().Nodes()[1])->Id() << endl;
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", SlaveElement().Id(), MasterElement().Id());
    }
  }

  else if  (!s0hasproj && !s1hasproj && !m0hasproj && m1hasproj)
  {
    if ((-1.0+MORTARPROJTOL<=mprojxi[1]) && (mprojxi[1]<=1.0-MORTARPROJTOL))
    {
      std::cout << "SElement Node IDs: " << (SlaveElement().Nodes()[0])->Id() << " " << (SlaveElement().Nodes()[1])->Id() << endl;
      std::cout << "MElement Node IDs: " << (MasterElement().Nodes()[0])->Id() << " " << (MasterElement().Nodes()[1])->Id() << endl;
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: IntegrateOverlap: Significant overlap ignored S%i M%i!", SlaveElement().Id(), MasterElement().Id());
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
      //std::cout << "Problem solved with internal case 1!" << endl;
    }

    // internal case 2 for global CASE 6
    // (equivalent to global CASE 8, master fully projects onto slave)
    else if ((mprojxi[0]<1.0) && (mprojxi[1]>-1.0))
    {
      mxia = -1.0;
      mxib = 1.0;
      sxia = mprojxi[1];      // local node numbering always anti-clockwise!!!
      sxib = mprojxi[0];
      //std::cout << "Problem solved with internal case 2!" << endl;
    }

    // internal case 3 for global CASE 6
    // (equivalent to global CASE 9, both nodes no. 0 project successfully)
    else if ((sprojxi[0]<1.0+MORTARPROJLIM) && (mprojxi[0]<1.0+MORTARPROJLIM))
    {
      sxia = -1.0;
      sxib = mprojxi[0];      // local node numbering always anti-clockwise!!!
      mxia = -1.0;
      mxib = sprojxi[0];
      //std::cout << "Problem solved with internal case 3!" << endl;
    }

    // internal case 4 for global CASE 6
    // (equivalent to global CASE 10, both nodes no. 1 project successfully)
    else if ((sprojxi[1]>-1.0-MORTARPROJLIM) && (mprojxi[1]>-1.0-MORTARPROJLIM))
    {
      sxia = mprojxi[1];
      sxib = 1.0;            // local node numbering always anti-clockwise!!!
      mxia = sprojxi[1];
      mxib = 1.0;
      //std::cout << "Problem solved with internal case 4!" << endl;
    }

    // unknown internal case for global CASE 6
    else
    {
      std::cout << "MORTAR::Coupling2d::DetectOverlap "<< endl << "has detected '4 projections'-case for Sl./Ma. pair "
                << SlaveElement().Id() << "/" << MasterElement().Id() << endl;
      std::cout << "SElement Node IDs: " << (SlaveElement().Nodes()[0])->Id() << " " << (SlaveElement().Nodes()[1])->Id() << endl;
      std::cout << "MElement Node IDs: " << (MasterElement().Nodes()[0])->Id() << " " << (MasterElement().Nodes()[1])->Id() << endl;
      std::cout << "SPROJXI_0: " << sprojxi[0] << " SPROJXI_1: " << sprojxi[1] << endl;
      std::cout << "MPROJXI_0: " << mprojxi[0] << " MPROJXI_1: " << mprojxi[1] << endl;
      dserror("ERROR: DetectOverlap: Unknown overlap case found in global case 6!");
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
    if ((sprojxi[0]>-1.0+MORTARPROJLIM) && (mprojxi[0]>-1.0+MORTARPROJLIM))
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
    if ((sprojxi[1]<1.0-MORTARPROJLIM) && (mprojxi[1]<1.0-MORTARPROJLIM))
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
    std::cout << "SElement: " << SlaveElement().NodeIds()[0] << " " << SlaveElement().NodeIds()[1] << endl;
    std::cout << "MElement: " << MasterElement().NodeIds()[0] << " " << MasterElement().NodeIds()[1] << endl;
    std::cout << "s0: " << s0hasproj << " s1: " << s1hasproj << endl;
    std::cout << "m0: " << m0hasproj << " m1: " << m1hasproj << endl;
    dserror("ERROR: IntegrateOverlap: Unknown overlap case found!");
  }

  // check for 1:1 node projections and for infeasible limits
  if ((sxia<-1.0) || (sxib>1.0) || (mxia<-1.0) || (mxib>1.0))
  {
    if (abs(sxia+1.0)<MORTARPROJLIM) sxia=-1.0;
    if (abs(sxib-1.0)<MORTARPROJLIM) sxib= 1.0;
    if (abs(mxia+1.0)<MORTARPROJLIM) mxia=-1.0;
    if (abs(mxib-1.0)<MORTARPROJLIM) mxib= 1.0;

    if ((sxia<-1.0) || (sxib>1.0) || (mxia<-1.0) || (mxib>1.0))
    {
      std::cout << "Slave: " << sxia << " " << sxib << endl;
      std::cout << "Master: " << mxia << " " << mxib << endl;
      dserror("ERROR: IntegrateOverlap: Determined infeasible limits!");
    }
  }

  // update integration limits in xiproj_
  xiproj_[0]=sxia;
  xiproj_[1]=sxib;
  xiproj_[2]=mxia;
  xiproj_[3]=mxib;

  // store overlap information
  overlap_ = overlap;

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate slave / master overlap (public)                 popp 04/08|
 *----------------------------------------------------------------------*/
bool MORTAR::Coupling2d::IntegrateOverlap()
{
  // explicitely defined shapefunction type needed
  if (shapefcn_ == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateOverlap called without specific shape function defined!");
  
  /**********************************************************************/
  /* INTEGRATION                                                        */
  /* Depending on overlap and the xiproj_ entries integrate the Mortar  */
  /* matrices D and M on the overlap of the current sl / ma pair.       */
  /**********************************************************************/

  // no integration if no overlap
  if (!overlap_) return false;

  // set segmentation status of all slave nodes
  // (hassegment_ of a slave node is true if ANY segment/cell
  // is integrated that contributes to this slave node)
  int nnodes = SlaveElement().NumNode();
  DRT::Node** mynodes = SlaveElement().Nodes();
  if (!mynodes) dserror("ERROR: Null pointer!");
  for (int k=0;k<nnodes;++k)
  {
    MORTAR::MortarNode* mycnode = static_cast<MORTAR::MortarNode*> (mynodes[k]);
    if (!mycnode) dserror("ERROR: Null pointer!");
    mycnode->HasSegment()=true;
  }
  
  //local working copies of input variables
  double sxia = xiproj_[0];
  double sxib = xiproj_[1];
  double mxia = xiproj_[2];
  double mxib = xiproj_[3];

  // create an integrator instance with correct NumGP and Dim
  MORTAR::MortarIntegrator integrator(shapefcn_,SlaveElement().Shape());

  // *******************************************************************
  // different options for mortar integration
  // *******************************************************************
  // (1) no quadratic element(s) involved -> linear LM interpolation
  // (2) quadratic element(s) involved -> quadratic LM interpolation
  // (3) quadratic element(s) involved -> linear LM interpolation
  // (4) quadratic element(s) involved -> piecew. linear LM interpolation
  // *******************************************************************
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

  // *******************************************************************
  // cases (1), (2) and (3)
  // *******************************************************************
  if (!Quad() ||
      (Quad() && lmtype==INPAR::MORTAR::lagmult_quad_quad) ||
      (Quad() && lmtype==INPAR::MORTAR::lagmult_lin_lin))
  {
    // do the overlap integration (integrate and linearize both M and gap)
    int nrow = SlaveElement().NumNode();
    int ncol = MasterElement().NumNode();
    int ndof = static_cast<MORTAR::MortarNode*>(SlaveElement().Nodes()[0])->NumDof();

    // mortar matrix dimensions depend on the actual number of slave / master DOFs
    // per node, which is NOT always identical to the problem dimension
    // (e.g. fluid meshtying -> 4 DOFs per node in a 3D problem)
    // -> thus the following check has been commented out (popp 04/2011)
    // if (ndof != Dim()) dserror("ERROR: Problem dimension and dofs per node not identical");

    Teuchos::RCP<Epetra_SerialDenseMatrix> dseg = Teuchos::rcp(new Epetra_SerialDenseMatrix(nrow*ndof,nrow*ndof));
    Teuchos::RCP<Epetra_SerialDenseMatrix> mseg = Teuchos::rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
    Teuchos::RCP<Epetra_SerialDenseVector> gseg = Teuchos::rcp(new Epetra_SerialDenseVector(nrow));
    integrator.IntegrateDerivSegment2D(SlaveElement(),sxia,sxib,MasterElement(),mxia,mxib,lmtype,dseg,mseg,gseg);

    // do the two assemblies into the slave nodes
    integrator.AssembleD(Comm(),SlaveElement(),*dseg);
    integrator.AssembleM(Comm(),SlaveElement(),MasterElement(),*mseg);
  }

  // *******************************************************************
  // case (4)
  // *******************************************************************
  else if (Quad() && lmtype==INPAR::MORTAR::lagmult_pwlin_pwlin)
  {
     dserror("ERROR: Piecewise linear LM interpolation not (yet?) implemented in 2D");
  }

  // *******************************************************************
  // other cases
  // *******************************************************************
  else
  {
    dserror("ERROR: IntegrateOverlap: Invalid case for 2D mortar coupling LM interpolation");
  }

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
  if (SlaveElement().Shape()==DRT::Element::line2)
  {
    if (integrator.Dim()==2)
    {
      MortarNode* snode0 = static_cast<MortarNode*>(SlaveElement().Nodes()[0]);
      MortarNode* snode1 = static_cast<MortarNode*>(SlaveElement().Nodes()[1]);

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
    Teuchos::RCP<Epetra_SerialDenseMatrix> mmodseg = integrator.IntegrateMmod2D(SlaveElement(),sxia,sxib,MasterElement(),mxia,mxib);
    integrator.AssembleMmod(Comm(),SlaveElement(),MasterElement(),*mmodseg);
  }
  //--------------------------------------------------------------------*/

  return true;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Coupling2dManager::Coupling2dManager(DRT::Discretization& idiscret,
                                             int dim, bool quad,
                                             INPAR::MORTAR::LagMultQuad lmtype,
                                             MORTAR::MortarElement* sele,
                                             std::vector<MORTAR::MortarElement*> mele) :
shapefcn_(INPAR::MORTAR::shape_undefined),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
lmtype_(lmtype),
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
MORTAR::Coupling2dManager::Coupling2dManager(const INPAR::MORTAR::ShapeFcn shapefcn,
                                             DRT::Discretization& idiscret,
                                             int dim, bool quad,
                                             INPAR::MORTAR::LagMultQuad lmtype,
                                             MORTAR::MortarElement* sele,
                                             std::vector<MORTAR::MortarElement*> mele) :
shapefcn_(shapefcn),
idiscret_(idiscret),
dim_(dim),
quad_(quad),
lmtype_(lmtype),
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
bool MORTAR::Coupling2dManager::EvaluateCoupling()
{
  // loop over all master elements associated with this slave element
  for (int m=0;m<(int)MasterElements().size();++m)
  {
    // create Coupling2d object and push back
    Coupling().push_back(Teuchos::rcp(new Coupling2d(
      shapefcn_,idiscret_,dim_,quad_,lmtype_,SlaveElement(),MasterElement(m))));

    // project the element pair
    Coupling()[m]->Project();

    // check for element overlap
    Coupling()[m]->DetectOverlap();

    // integrate the element overlap
    Coupling()[m]->IntegrateOverlap();
  }

  return true;
}

