/*!----------------------------------------------------------------------
\file mortar_integrator.cpp
\brief A class to perform integrations of Mortar matrices on the overlap
        of two MortarElements in 1D and 2D

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

#include "mortar_integrator.H"
#include "mortar_node.H"
#include "mortar_element.H"
#include "mortar_projector.H"
#include "mortar_coupling3d_classes.H"
#include "mortar_defines.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 07/09|
 *----------------------------------------------------------------------*/
MORTAR::MortarIntegrator::MortarIntegrator(Teuchos::ParameterList& params,
                               DRT::Element::DiscretizationType eletype) :
imortar_(params)
{
  InitializeGP(eletype);
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points                                   popp 06/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::InitializeGP(DRT::Element::DiscretizationType eletype)
{
  //**********************************************************************
  // Create integration points according to eletype!
  // Note that our standard Gauss rules are:
  // 5 points: for integrals on 1D lines            (1,2,3,4,5,6,7,8,9,10)
  //           -> degree of precision: 9
  // 7 points: for integrals on 2D triangles        (1,3,6,7,12,37,64)
  //           -> degree of precision: 5
  // 9 points: for integrals on 2D quadrilaterals   (1,4,9,16,25,36,49,64)
  //           -> degree of precision: 5
  //**********************************************************************

  // get input
  const Teuchos::ParameterList& smortar  = DRT::Problem::Instance()->MortarCouplingParams();
  int numgp = smortar.get<int>("NUMGP_PER_DIM");

  INPAR::MORTAR::IntType integrationtype =
    DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(smortar,"INTTYPE");

  switch(eletype)
  {
  case DRT::Element::line2:
  case DRT::Element::line3:
  {
    // set default value for Gauss rule first
    dim_=2;
    DRT::UTILS::GaussRule1D mygaussrule = DRT::UTILS::intrule_line_5point;

    // GP switch if non-zero value provided by user
    if(integrationtype==INPAR::MORTAR::inttype_elements ||integrationtype==INPAR::MORTAR::inttype_elements_BS)
    {
      if (numgp>0)
      {
        switch(numgp)
        {
          case 1:
          {
            mygaussrule = DRT::UTILS::intrule_line_1point;
            break;
          }
          case 2:
          {
            mygaussrule = DRT::UTILS::intrule_line_2point;
            break;
          }
          case 3:
          {
            mygaussrule = DRT::UTILS::intrule_line_3point;
            break;
          }
          case 4:
          {
            mygaussrule = DRT::UTILS::intrule_line_4point;
            break;
          }
          case 5:
          {
            mygaussrule = DRT::UTILS::intrule_line_5point;
            break;
          }
          case 6:
          {
            mygaussrule = DRT::UTILS::intrule_line_6point;
            break;
          }
          case 7:
          {
            mygaussrule = DRT::UTILS::intrule_line_7point;
            break;
          }
          case 8:
          {
            mygaussrule = DRT::UTILS::intrule_line_8point;
            break;
          }
          case 9:
          {
            mygaussrule = DRT::UTILS::intrule_line_9point;
            break;
          }
          case 10:
          {
            mygaussrule = DRT::UTILS::intrule_line_10point;
            break;
          }
          case 16:
          {
            mygaussrule = DRT::UTILS::intrule_line_16point;
            break;
          }
          case 20:
          {
            mygaussrule = DRT::UTILS::intrule_line_20point;
            break;
          }
          case 32:
          {
            mygaussrule = DRT::UTILS::intrule_line_32point;
            break;
          }
          default:
          {
            dserror("Requested GP-Number is not implemented!");
            break;
          }
        }
      }
    }

    const DRT::UTILS::IntegrationPoints1D intpoints(mygaussrule);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),2);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=0.0;
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    // set default value for Gauss rule first
    dim_=3;
    DRT::UTILS::GaussRule2D mygaussrule=DRT::UTILS::intrule_tri_7point;

    // GP switch if non-zero value provided by user
    if(integrationtype==INPAR::MORTAR::inttype_elements)// ||integrationtype==INPAR::MORTAR::inttype_elements_BS)
    {
      if (numgp>0)
      {
        switch(numgp)
        {
          case 1:
          {
            mygaussrule = DRT::UTILS::intrule_tri_1point;
            break;
          }
          case 2:
          {
            mygaussrule = DRT::UTILS::intrule_tri_6point;
            break;
          }
          case 3:
          {
            mygaussrule = DRT::UTILS::intrule_tri_7point;
            break;
          }
          case 4:
          {
            mygaussrule = DRT::UTILS::intrule_tri_12point;
            break;
          }
          case 5:
          {
            mygaussrule = DRT::UTILS::intrule_tri_12point;
            break;
          }
          case 6:
          {
            mygaussrule = DRT::UTILS::intrule_tri_37point;
            break;
          }
          case 7:
          {
            mygaussrule = DRT::UTILS::intrule_tri_37point;
            break;
          }
          case 8:
          {
            mygaussrule = DRT::UTILS::intrule_tri_64point;
            break;
          }
          case 9:
          {
            mygaussrule = DRT::UTILS::intrule_tri_64point;
            break;
          }
          case 10:
          {
            mygaussrule = DRT::UTILS::intrule_tri_64point;
            break;
          }
          case 20:
          {
            mygaussrule = DRT::UTILS::intrule_tri_64point;
            break;
          }
          default:
          {
            dserror("Requested GP-Number is not implemented!");
            break;
          }
        }
      }
    }

    const DRT::UTILS::IntegrationPoints2D intpoints(mygaussrule);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),2);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=intpoints.qxg[i][1];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    // set default value for Gauss rule first
    dim_=3;
    DRT::UTILS::GaussRule2D mygaussrule=DRT::UTILS::intrule_quad_9point;

    // GP switch if non-zero value provided by user
    if(integrationtype==INPAR::MORTAR::inttype_elements ||integrationtype==INPAR::MORTAR::inttype_elements_BS)
    {
      if (numgp>0)
      {
        switch(numgp)
        {
          case 2:
          {
            mygaussrule = DRT::UTILS::intrule_quad_4point;
            break;
          }
          case 3:
          {
            mygaussrule = DRT::UTILS::intrule_quad_9point;
            break;
          }
          case 4:
          {
            mygaussrule = DRT::UTILS::intrule_quad_16point;
            break;
          }
          case 5:
          {
            mygaussrule = DRT::UTILS::intrule_quad_25point;
            break;
          }
          case 6:
          {
            mygaussrule = DRT::UTILS::intrule_quad_36point;
            break;
          }
          case 7:
          {
            mygaussrule = DRT::UTILS::intrule_quad_49point;
            break;
          }
          case 8:
          {
            mygaussrule = DRT::UTILS::intrule_quad_64point;
            break;
          }
          case 9:
          {
            mygaussrule = DRT::UTILS::intrule_quad_81point;
            break;
          }
          case 10:
          {
            mygaussrule = DRT::UTILS::intrule_quad_100point;
            break;
          }
          case 16:
          {
            mygaussrule = DRT::UTILS::intrule_quad_256point;
            break;
          }
          case 20:
          {
            mygaussrule = DRT::UTILS::intrule_quad_400point;
            break;
          }
          case 32:
          {
            mygaussrule = DRT::UTILS::intrule_quad_1024point;
            break;
          }
          default:
          {
            dserror("Requested GP-Number is not implemented!");
            break;
          }
        }
      }
    }

    const DRT::UTILS::IntegrationPoints2D intpoints(mygaussrule);
    ngp_ = intpoints.nquad;
    coords_.Reshape(nGP(),2);
    weights_.resize(nGP());
    for (int i=0;i<nGP();++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=intpoints.qxg[i][1];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  default:
  {
    dserror("ERROR: MortarIntegrator: This contact element type is not implemented!");
    break;
  }
  } // switch(eletype)
  
  return;
}

/*--------------------------------------------------------------------------------------*
 | Integrate without segmentation --> more GP required                       farah 01/13|
 | Integration over the entire Slave-Element: no mapping sxi->eta                       |
 | required                                                                             |
 *--------------------------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::EleBased_Integration(
       Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
       Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
       MORTAR::MortarElement& sele,
       std::vector<MORTAR::MortarElement*> meles,
       bool *boundary_ele)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // number of nodes (slave, master)
  int nrow    = sele.NumNode();
  int ndof    = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();
  int mndof   = static_cast<MORTAR::MortarNode*>(meles[0]->Nodes()[0])->NumDof();
  int nodemaster  = meles[0]->NumNode();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(nodemaster);
  LINALG::SerialDenseMatrix mderiv(nodemaster,1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,1);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  if (lmtype == INPAR::MORTAR::lagmult_lin_lin && sele.Shape() == DRT::Element::line3)
  {
    bound = false; // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
  }

  double sxia=-1;
  double sxib=1;
  double sxi[2]={0.0 , 0.0};

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp) //loop to the end //Number of the GP
  {
    bool is_on_mele=false;
    bool kink_projection=false;
    // coordinates and weight of the GP
    double eta[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (linlm)
      sele.EvaluateShapeLagMultLin(ShapeFcn(),sxi,lmval,lmderiv,nrow);//nrow
    else
      sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

    sele.EvaluateShape(sxi,sval,sderiv,nrow);

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;


    //********************************************************************
    //  loop over all involved masterelements
    //********************************************************************
    for (int nummaster=0;nummaster<(int)meles.size();++nummaster)
    {
      // get Master element nodes themselves
      DRT::Node** mnodes = meles[nummaster]->Nodes();
      if(!mnodes) dserror("ERROR: EleBased_Integration: Null pointer!");

      // project Gauss point onto master element
      double mxi[2] = {0.0, 0.0};
      sxi[0]=eta[0];

      MORTAR::MortarProjector projector(2); //2D Projector
      projector.ProjectGaussPoint(sele,eta,*meles[nummaster],mxi);

      // evaluate trace space shape functions (on both elements)
      meles[nummaster]->EvaluateShape(mxi,mval,mderiv,meles[nummaster]->NumNode());

      // check GP projection
      if ((mxi[0]>=-1) && (mxi[0]<=1) && (kink_projection==false))
      {
        kink_projection=true;
        is_on_mele=true;
        //***************************
        // standard shape functions
        //***************************
        if (ShapeFcn() == INPAR::MORTAR::shape_standard)
        {
          // loop over all mseg matrix entries
          // !!! nrow represents the slave Lagrange multipliers !!!
          // !!! ncol represents the master dofs                !!!
          // (this DOES matter here for mseg, as it might
          // sometimes be rectangular, not quadratic!)
          for (int j=0; j<nrow*ndof; ++j)
          {
            if (meles.size()>0)
            {
              // integrate mseg
              for (int k=0; k<nodemaster*mndof; ++k)
              {
                int jindex = (int)(j/ndof);
                int kindex = (int)(k/mndof);

                  DRT::Node** mnodes = meles[nummaster]->Nodes();
                MORTAR::MortarNode* mymnode = static_cast<MORTAR::MortarNode*> (mnodes[kindex]);
                if (!mymnode) dserror("ERROR: Null pointer!");

                // multiply the two shape functions
                double prod = lmval[jindex]*mval[kindex];

                // isolate the mseg and dseg entries to be filled
                // (both the main diagonal and every other secondary diagonal)
                // and add current Gauss point's contribution to mseg and dseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*mndof)))
                {
                (*mseg)(j, k+nummaster*nodemaster*mndof) += prod*dxdsxi*dsxideta*wgt;
                }
              }
            }
          }
        }//std shape

        //***********************
        // dual shape functions
        //***********************
        else if (ShapeFcn() == INPAR::MORTAR::shape_dual)
        {
          //dserror("ERROR: noch nicht vervollständigt");
          // loop over all mseg matrix entries
          // nrow represents the slave Lagrange multipliers !!!
          // ncol represents the master dofs !!!
          // (this DOES matter here for mseg, as it might
          // sometimes be rectangular, not quadratic!)
          for (int j=0;j<nrow*ndof;++j)
          {
            // for dual shape functions we can make use
            // of the row summing lemma: D_jj = Sum(k) M_jk
            // hence, they can be combined into one single loop

            // integrate mseg and dseg (no boundary modification)
            for (int k=0;k<nodemaster*mndof;++k)
            {
              int jindex = (int)(j/ndof);
              int kindex = (int)(k/ndof);

              // multiply the two shape functions
              double prod = lmval[jindex]*mval[kindex];

              // isolate the mseg and dseg entries to be filled
              // (both the main diagonal and every other secondary diagonal for mseg)
              // (only the main diagonal for dseg)
              // and add current Gauss point's contribution to mseg
              if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              {
                (*mseg)(j, k+nummaster*nodemaster*mndof) += prod*dxdsxi*dsxideta*wgt;
                if(!bound)
                  {
                    (*dseg)(j,j) += prod*dxdsxi*dsxideta*wgt;
                  }
              }
            }
            // integrate dseg (boundary modification)
            if (bound)
            {
              MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[(int)(j/ndof)]);
              if (!mymrtrnode) dserror("ERROR: EleBased_integration: Null pointer!");
              bool j_boundnode = mymrtrnode->IsOnBound();

              for (int k=0;k<nrow*ndof;++k)
              {
                MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[(int)(k/ndof)]);
                if (!mymrtrnode2) dserror("ERROR: EleBased_Integration: Null pointer!");
                bool k_boundnode = mymrtrnode2->IsOnBound();

                int jindex = (int)(j/ndof);
                int kindex = (int)(k/ndof);

                // do not assemble off-diagonal terms if j,k are both non-boundary nodes
                if (!j_boundnode && !k_boundnode && (jindex!=kindex)) continue;

                // multiply the two shape functions
                double prod = lmval[jindex]*sval[kindex];

                // isolate the dseg entries to be filled
                // (both the main diagonal and every other secondary diagonal)
                // and add current Gauss point's contribution to dseg
                if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                  (*dseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
              }
            }
          } //j-loop
       } //end - if dualshapefunction
      } // end - if Pojection on MasterElement
    }//loop-end over all involved master elements

    // strong discontinuity --> Boundary Segmentation
    if(is_on_mele==false)
    {
      *boundary_ele=true;
    }

    //**********************************************************
    //Dseg
    //**********************************************************
    if (ShapeFcn() == INPAR::MORTAR::shape_standard)
    {
      if(is_on_mele==true)
      {
        //integrate dseg
        for (int j=0; j<nrow*ndof; ++j)
        {
          for (int k=0; k<nrow*ndof; ++k)
          {
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);

            // multiply the two shape functions
            double prod = lmval[jindex]*sval[kindex];

            // isolate the mseg and dseg entries to be filled
            // (both the main diagonal and every other secondary diagonal)
            // and add current Gauss point's contribution to mseg and dseg
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            {
              (*dseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
            }
          }
        }
      }
    }
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    {
      if (is_on_mele==true)
      {
        for (int j=0;j<nrow*ndof;++j)
        {
        // integrate dseg (boundary modification)
        if (bound)
        {
          MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[(int)(j/ndof)]);
          if (!mymrtrnode) dserror("ERROR: EleBased_integration: Null pointer!");
          bool j_boundnode = mymrtrnode->IsOnBound();

          for (int k=0;k<nrow*ndof;++k)
          {
          MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[(int)(k/ndof)]);
          if (!mymrtrnode2) dserror("ERROR: EleBased_Integration: Null pointer!");
          bool k_boundnode = mymrtrnode2->IsOnBound();

          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // do not assemble off-diagonal terms if j,k are both non-boundary nodes
          if (!j_boundnode && !k_boundnode && (jindex!=kindex)) continue;

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
          }
        }
        }
      }
    }

    }//loop-end over all gp

    //consistency check
    bool zero=true;
    for (int m=0;m<(int)((meles.size())*4);++m)
    {
      for (int n=0;n<(int)4;++n)
      {
        if((*mseg)(n,m)!=0)
        {
          zero=false;
        }
      }
    }
    if (zero==true)
    {
      for (int i=0;i<(int)4;++i)
      {
        for (int j=0; j<(int)4;++j)
        {
          (*dseg)(i,j)=0.0;
        }
      }
    }

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize D on a slave element (2D / 3D)    popp 02/09|
 |  This method integrates the element D matrix and stores it in dseg.  |
 |  Moreover, derivatives LinD are built and stored directly into the   |
 |  adajcent nodes. (Thus this method combines EVERYTHING before done   |
 |  aperarately in IntegrateD and DerivD!)                              |
 |  ********** modified version: responds to ShapeFcn() ***   popp 05/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::IntegrateDerivSlave2D3D(
      MORTAR::MortarElement& sele, double* sxia, double* sxib,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg)
{
  //**********************************************************************
  dserror("ERROR: IntegrateDerivSlave2D3D method is outdated!");
  //**********************************************************************

  // explicitely defined shapefunction type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSlave2D3D called without specific shape function defined!");
    
  //check input data
  if (!sele.IsSlave())
    dserror("ERROR: IntegrateDerivSlave2D3D called on a non-slave MortarElement!");
  if ((sxia[0]<-1.0) || (sxia[1]<-1.0) || (sxib[0]>1.0) || (sxib[1]>1.0))
    dserror("ERROR: IntegrateDerivSlave2D3D called with infeasible slave limits!");

  // number of nodes (slave)
  int nrow = sele.NumNode();
  int ncol = sele.NumNode();
  int ndof = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty objects for shape fct. evaluation
  LINALG::SerialDenseVector val(nrow);
  LINALG::SerialDenseMatrix deriv(nrow,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), 0.0};
    if (Dim()==3) eta[1] = Coordinate(gp,1);
    double wgt = Weight(gp);

    // evaluate trace space and Lagrange multiplier shape functions
    sele.EvaluateShape(eta,val,deriv,nrow);
    sele.EvaluateShapeLagMult(ShapeFcn(),eta,lmval,lmderiv,nrow);

    // evaluate the Jacobian det
    double dxdsxi = sele.Jacobian(eta);

    // compute element D matrix ******************************************
    // loop over all dtemp matrix entries
    // nrow represents the Lagrange multipliers !!!
    // ncol represents the dofs !!!
    for (int j=0;j<nrow*ndof;++j)
    {
      for (int k=0;k<ncol*ndof;++k)
      {
        int jindex = (int)(j/ndof);
        int kindex = (int)(k/ndof);

        // multiply the two shape functions
        double prod = lmval[jindex]*val[kindex];

        // isolate the dseg entries to be filled
        // (both the main diagonal and every other secondary diagonal)
        // and add current Gauss point's contribution to dseg
        if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          (*dseg)(j,k) += prod*dxdsxi*wgt;
      }
    }
    // compute element D matrix ******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 1D slave / master overlap (2D)  popp 02/09|
 |  This method integrates the overlap M matrix and weighted gap g~     |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM, IntegrateG, DerivM and DerivG!)                         |
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::IntegrateDerivSegment2D(
     MORTAR::MortarElement& sele, double& sxia, double& sxib,
     MORTAR::MortarElement& mele, double& mxia, double& mxib,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{ 
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

  // explicitely defined shapefunction type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivSegment2D called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateAndDerivSegment called on a wrong type of MortarElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
    dserror("ERROR: IntegrateAndDerivSegment called with infeasible master limits!");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,1);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateAndDerivSegment: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }
  
  // decide whether linear LM are used for quadratic FE here
  bool linlm = false;
  if (lmtype == INPAR::MORTAR::lagmult_lin_lin && sele.Shape() == DRT::Element::line3)
  {
    bound = false; // crosspoints and linear LM NOT at the same time!!!!
    linlm = true;
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    double sxi[2] = {0.0, 0.0};
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // project Gauss point onto master element
    double mxi[2] = {0.0, 0.0};
    MORTAR::MortarProjector projector(2);
    projector.ProjectGaussPoint(sele,sxi,mele,mxi);

    // simple version (no Gauss point projection)
    // double test[2] = {0.0, 0.0};
    // test[0] = 0.5*(1-eta[0])*mxib + 0.5*(1+eta[0])*mxia;
    // mxi[0]=test[0];

    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      dserror("ERROR: IntegrateAndDerivSegment: Gauss point projection failed! mxi=%d",mxi[0]);
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (linlm)
      sele.EvaluateShapeLagMultLin(ShapeFcn(),sxi,lmval,lmderiv,nrow);
    else
      sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);
    
    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;
    
    // compute segment D/M matrix ****************************************
    // standard shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_standard)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to mseg and dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
          }
        }

        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal)
          // and add current Gauss point's contribution to mseg and dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*dseg)(j, k) += prod*dxdsxi*dsxideta*wgt;
          }
        }
      } // nrow*ndof loop
    }
    
    // dual shape functions
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    { 
      // loop over all mseg matrix entries
      // nrow represents the slave Lagrange multipliers !!!
      // ncol represents the master dofs !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0;j<nrow*ndof;++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(k) M_jk
        // hence, they can be combined into one single loop

        // integrate mseg and dseg (no boundary modification)
        for (int k=0;k<ncol*ndof;++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg and dseg entries to be filled
          // (both the main diagonal and every other secondary diagonal for mseg)
          // (only the main diagonal for dseg)
          // and add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
            if(!bound) (*dseg)(j,j) += prod*dxdsxi*dsxideta*wgt;
          }
        }

        // integrate dseg (boundary modification)
        if (bound)
        {
          MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[(int)(j/ndof)]);
          if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
          bool j_boundnode = mymrtrnode->IsOnBound();
          
          for (int k=0;k<nrow*ndof;++k)
          {
            MORTAR::MortarNode* mymrtrnode2 = static_cast<MORTAR::MortarNode*>(mynodes[(int)(k/ndof)]);
            if (!mymrtrnode2) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
            bool k_boundnode = mymrtrnode2->IsOnBound();
            
            int jindex = (int)(j/ndof);
            int kindex = (int)(k/ndof);
            
            // do not assemble off-diagonal terms if j,k are both non-boundary nodes
            if (!j_boundnode && !k_boundnode && (jindex!=kindex)) continue;
              
            // multiply the two shape functions 
            double prod = lmval[jindex]*sval[kindex];
    
            // isolate the dseg entries to be filled
            // (both the main diagonal and every other secondary diagonal)
            // and add current Gauss point's contribution to dseg
            if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              (*dseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
          }
        }
      }
    } // ShapeFcn() switch
    // compute segment D/M matrix ****************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate Mmod on slave / master overlap (2D)             popp 01/08|
 |  This method integrates the modification to the Mortar matrix M      |
 |  for curved interface (Paper by Puso/Wohlmuth) from given local      |
 |  coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib                          |
 |  Output is an Epetra_SerialDenseMatrix holding the int. values       |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_SerialDenseMatrix> MORTAR::MortarIntegrator::IntegrateMmod2D(
                                       MORTAR::MortarElement& sele,
                                       double& sxia, double& sxib,
                                       MORTAR::MortarElement& mele,
                                       double& mxia, double& mxib)
{
  //**********************************************************************
  dserror("ERROR: IntegrateMmod2D method is outdated!");
  //**********************************************************************

  //check for problem dimension
  if (Dim()!=2) dserror("ERROR: 2D integration method called for non-2D problem");

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateMmod2D called on a wrong type of MortarElement pair!");
  if ((sxia<-1.0) || (sxib>1.0))
    dserror("ERROR: IntegrateMmod2D called with infeasible slave limits!");
  if ((mxia<-1.0) || (mxib>1.0))
      dserror("ERROR: IntegrateMmod2D called with infeasible master limits!");

  // create empty mmodseg object and wrap it with RCP
  int nrow  = sele.NumNode();
  int nrowdof = Dim();
  int ncol  = mele.NumNode();
  int ncoldof = Dim();
  Teuchos::RCP<Epetra_SerialDenseMatrix> mmodseg = Teuchos::rcp(new Epetra_SerialDenseMatrix(nrow*nrowdof,ncol*ncoldof));

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,1);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,1);

  // loop over all Gauss points for integration
  for (int gp=0;gp<nGP();++gp)
  {
    double eta[2] = {Coordinate(gp,0), 0.0};
    double wgt = Weight(gp);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // coordinate transformation sxi->eta (slave MortarElement->Overlap)
    sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;

    // project Gauss point onto master element
    MORTAR::MortarProjector projector(2);
    projector.ProjectGaussPoint(sele,sxi,mele,mxi);

    // check GP projection
    if ((mxi[0]<mxia) || (mxi[0]>mxib))
    {
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
      std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      dserror("ERROR: IntegrateMmod2D: Gauss point projection failed!");
    }

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // build the delta function of slave side shape functions
    double deltasval = sval[0]-sval[1];

    // evaluate the two slave side Jacobians
    double dxdsxi = sele.Jacobian(sxi);
    double dsxideta = -0.5*sxia + 0.5*sxib;

    /* loop over all mmodseg matrix entries
       nrow represents the slave Lagrange multipliers !!!
       ncol represents the master dofs !!!
       (this DOES matter here for mmodseg, as it might
       sometimes be rectangular, not quadratic!)              */
    for (int j=0;j<nrow*nrowdof;++j)
    {
      for (int k=0;k<ncol*ncoldof;++k)
      {
        // multiply the two shape functions
        int mindex = (int)(k/ncoldof);
        double prod = 0.5*deltasval*mval[mindex];
        // add current Gauss point's contribution to mmodseg
        (*mmodseg)(j,k) += prod*dxdsxi*dsxideta*wgt;
      }
    }
  } // for (int gp=0;gp<nGP();++gp)

  // prepare computation of purely geometric part of Mmod entries
  MortarNode* snode0 = static_cast<MortarNode*>(sele.Nodes()[0]);
  MortarNode* snode1 = static_cast<MortarNode*>(sele.Nodes()[1]);

  // normals
  double n[2][2];
  n[0][0] = snode0->MoData().n()[0];
  n[0][1] = snode0->MoData().n()[1];
  n[1][0] = snode1->MoData().n()[0];
  n[1][1] = snode1->MoData().n()[1];

  // scalar product n1 * n2
  double n1n2 = 0.0;
  for (int i=0;i<2;++i)
    n1n2+=n[0][i]*n[1][i];

  // vector product n1 x n2
  double n1xn2 = n[0][0]*n[1][1] - n[0][1]*n[1][0];

  // // multiply geometric part onto Mmod
  for (int i=0;i<ncol;++i)
  {
    (*mmodseg)(0,0+i*ncoldof) *=  (1.0-n1n2);
    (*mmodseg)(1,0+i*ncoldof) *=  n1xn2;
    (*mmodseg)(0,1+i*ncoldof) *= -n1xn2;
    (*mmodseg)(1,1+i*ncoldof) *=  (1.0-n1n2);

    (*mmodseg)(2,0+i*ncoldof) *=  (n1n2-1.0);
    (*mmodseg)(3,0+i*ncoldof) *=  n1xn2;
    (*mmodseg)(2,1+i*ncoldof) *= -n1xn2;
    (*mmodseg)(3,1+i*ncoldof) *=  (n1n2-1.0);
  }

  return mmodseg;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 02/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::IntegrateDerivCell3D(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     Teuchos::RCP<MORTAR::IntCell> cell,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType dt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3D called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3D called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);
  
  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // note that the third component of sxi is necessary!
    // (although it will always be 0.0 of course)
    double tempsxi[3] = {0.0, 0.0, 0.0};
    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};
    double projalpha = 0.0;

    // get Gauss point in slave element coordinates
    cell->LocalToGlobal(eta,tempsxi,0);
    sxi[0] = tempsxi[0];
    sxi[1] = tempsxi[1];

    // project Gauss point onto master element
    MORTAR::MortarProjector projector(3);
    projector.ProjectGaussPoint3D(sele,sxi,mele,mxi,projalpha);

    // check GP projection
    double tol = 0.01;
    if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
        std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
        std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // evaluate Lagrange multiplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the two Jacobians (int. cell and slave element)
    double jaccell = cell->Jacobian(eta);
    double jacslave = sele.Jacobian(sxi);

    // compute cell D/M matrix *******************************************
    // standard shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_standard)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg  
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*mseg)(j, k) += prod*jaccell*jacslave*wgt;
        }
        
        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              (*dseg)(j, k) += prod*jaccell*jacslave*wgt;
        }
      } 
    }
    
    // dual shape functions
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(k) M_jk
        // hence, they can be combined into one single loop
        
        // integrate mseg and dseg
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);
  
          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg  
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*jaccell*jacslave*wgt;
            (*dseg)(j, j) += prod*jaccell*jacslave*wgt;
          }
        }
      }
    }
    // compute cell D/M matrix *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize without segmentation             farah 01/13|
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::IntegrateDerivCell3D_EleBased(
     MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg,
     bool *boundary_ele)
{
  // explicitely defined shapefunction type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");

  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType dt = meles[0]->Shape();

  // check input data
  for (int test=0; test<(int)meles.size();++test)
  {
    if ((!sele.IsSlave()) || (meles[test]->IsSlave()))
      dserror("ERROR: IntegrateDerivCell3D called on a wrong type of MortarElement pair!");
  }

  int msize   = meles.size();
  int nrow    = sele.NumNode();
  int nmnode  = meles[0]->NumNode();
  int ndof    = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(nmnode);
  LINALG::SerialDenseMatrix mderiv(nmnode,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // note that the third component of sxi is necessary!
    // (although it will always be 0.0 of course)
    //double tempsxi[3] = {0.0, 0.0, 0.0};
    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};
    double projalpha = 0.0;

    sxi[0]=eta[0];
    sxi[1]=eta[1];

    //check for Boundary Segmentation
    bool projactable_gp=false;

    //*******************************************************************
    // loop over meles
    //*******************************************************************
    for(int nummaster=0;nummaster<msize;++nummaster)
    {
      // project Gauss point onto master element
      MORTAR::MortarProjector projector(3);
      projector.ProjectGaussPoint3D(sele,sxi,*meles[nummaster],mxi,projalpha);

      bool is_on_mele=true;

      // check GP projection
      double tol = 0.00;
      if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
      {
        if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
        {
//          std::cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
//          std::cout << "Slave ID: " << sele.Id() << " Master ID: " << meles[nummaster]->Id() << std::endl;
//          std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
//          std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
//          std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;

          is_on_mele=false;
        }
      }
      else
      {
        if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
        {
//          std::cout << "\n***Warning: IntegrateDerivCell3D: Gauss point projection outside!";
//          std::cout << "Slave ID: " << sele.Id() << " Master ID: " << meles[nummaster]->Id() << std::endl;
//          std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
//          std::cout << "Gauss point: " << sxi[0] << " " << sxi[1] << std::endl;
//          std::cout << "Projection: " << mxi[0] << " " << mxi[1] << std::endl;

          is_on_mele=false;
        }
      }

      if (is_on_mele==true)
      {
        projactable_gp=true;

        // evaluate Lagrange multiplier shape functions (on slave element)
        sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

        // evaluate trace space shape functions (on both elements)
        sele.EvaluateShape(sxi,sval,sderiv,nrow);
        meles[nummaster]->EvaluateShape(mxi,mval,mderiv,nmnode);

        // evaluate the two Jacobians (int. cell and slave element)
        //double jaccell = cell->Jacobian(eta);
        double jacslave = sele.Jacobian(sxi);

        // compute cell D/M matrix *******************************************
        // standard shape functions
        if (ShapeFcn() == INPAR::MORTAR::shape_standard)
        {
          // loop over all mseg matrix entries
          // !!! nrow represents the slave Lagrange multipliers !!!
          // !!! ncol represents the master dofs                !!!
          // (this DOES matter here for mseg, as it might
          // sometimes be rectangular, not quadratic!)
          for (int j=0; j<nrow*ndof; ++j)
          {
            // integrate mseg
            for (int k=0; k<nmnode*ndof; ++k)
            {
              int jindex = (int)(j/ndof);
              int kindex = (int)(k/ndof);

              // multiply the two shape functions
              double prod = lmval[jindex]*mval[kindex];

              // isolate the mseg entries to be filled and
              // add current Gauss point's contribution to mseg
              if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                (*mseg)(j, k+nummaster*nmnode*ndof)+= prod*jacslave*wgt; //*jaccell
            }

            // integrate dseg
            for (int k=0; k<nrow*ndof; ++k)
            {
              int jindex = (int)(j/ndof);
              int kindex = (int)(k/ndof);

              // multiply the two shape functions
              double prod = lmval[jindex]*sval[kindex];

              // isolate the mseg entries to be filled and
              // add current Gauss point's contribution to mseg
              if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
                  (*dseg)(j, k) += prod*jacslave*wgt;//*jaccell
            }
          }
        }

        // dual shape functions
        else if (ShapeFcn() == INPAR::MORTAR::shape_dual)
        {
          // loop over all mseg matrix entries
          // !!! nrow represents the slave Lagrange multipliers !!!
          // !!! ncol represents the master dofs                !!!
          // (this DOES matter here for mseg, as it might
          // sometimes be rectangular, not quadratic!)
          for (int j=0; j<nrow*ndof; ++j)
          {
            // for dual shape functions we can make use
            // of the row summing lemma: D_jj = Sum(k) M_jk
            // hence, they can be combined into one single loop

            // integrate mseg and dseg
            for (int k=0; k<nmnode*ndof; ++k)
            {
              int jindex = (int)(j/ndof);
              int kindex = (int)(k/ndof);

              // multiply the two shape functions
              double prod = lmval[jindex]*mval[kindex];

              // isolate the mseg entries to be filled and
              // add current Gauss point's contribution to mseg
              if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
              {
                (*mseg)(j, k+nummaster*nmnode*ndof) += prod*jacslave*wgt; //*jaccell
                (*dseg)(j, j) += prod*jacslave*wgt; //*jaccell
              }
            }
          }
        }
      }//is_on_mele==true
    }//loop over meles

    // strong discontinuity --> Boundary Segmentation
    if(projactable_gp==false)
    {
      *boundary_ele=true;
    }

  }

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix and weighted gap g~        |
 |  and stores it in mseg and gseg respectively. Moreover, derivatives  |
 |  LinM and Ling are built and stored directly into the adjacent nodes.|
 |  (Thus this method combines EVERYTHING before done separately in     |
 |  IntegrateM3D, IntegrateG3D, DerivM3D and DerivG3D!)                 |
 |  This is the auxiliary plane coupling version!!!                     |
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::IntegrateDerivCell3DAuxPlane(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // explicitely defined shapefunction type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without specific shape function defined!");
    
  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master element
  DRT::Element::DiscretizationType sdt = sele.Shape();
  DRT::Element::DiscretizationType mdt = mele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlane called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int ndof = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta,globgp,0);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // project Gauss point onto slave element
    // project Gauss point onto master element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector projector(3);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,sele,sxi,sprojalpha);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,mele,mxi,mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt==DRT::Element::quad4 || sdt==DRT::Element::quad8 || sdt==DRT::Element::quad9)
    {
      if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // evaluate Lagrange mutliplier shape functions (on slave element)
    sele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmval,lmderiv,nrow);

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(sxi,sval,sderiv,nrow);
    mele.EvaluateShape(mxi,mval,mderiv,ncol);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // compute cell D/M matrix *******************************************
    // standard shape functions
    if (ShapeFcn() == INPAR::MORTAR::shape_standard)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*mseg)(j, k) += prod*jac*wgt;
        }

        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j, k) += prod*jac*wgt;
        }
      }
    }
    
    // dual shape functions
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual)
    {
      // loop over all mseg matrix entries
      // !!! nrow represents the slave Lagrange multipliers !!!
      // !!! ncol represents the master dofs                !!!
      // (this DOES matter here for mseg, as it might
      // sometimes be rectangular, not quadratic!)
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(k) M_jk
        // hence, they can be combined into one single loop

        // integrate mseg and dseg
        for (int k=0; k<ncol*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[kindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
          {
            (*mseg)(j, k) += prod*jac*wgt;
            (*dseg)(j, j) += prod*jac*wgt;
          }
        }
      } // nrow*ndof loop
    } // ShapeFcn() switch
    // compute cell D/M matrix *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 01/11|
 |  This method integrates the cell M matrix (and possibly D matrix)    |
 |  and stores it in mseg and dseg respectively.                        |
 |  This is the QUADRATIC coupling version!!!                           |
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::IntegrateDerivCell3DQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     Teuchos::RCP<MORTAR::IntCell> cell,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();

  // explicitely defined shapefunction type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DQuad called without specific shape function defined!");

  /*std::cout << std::endl;
  std::cout << "Slave type: " << sele.Shape() << std::endl;
  std::cout << "SlaveElement Nodes:";
  for (int k=0;k<sele.NumNode();++k) std::cout << " " << sele.NodeIds()[k];
  std::cout << "\nMaster type: " << mele.Shape() << std::endl;
  std::cout << "MasterElement Nodes:";
  for (int k=0;k<mele.NumNode();++k) std::cout << " " << mele.NodeIds()[k];
  std::cout << std::endl;
  std::cout << "SlaveSub type: " << sintele.Shape() << std::endl;
  std::cout << "SlaveSubElement Nodes:";
  for (int k=0;k<sintele.NumNode();++k) std::cout << " " << sintele.NodeIds()[k];
  std::cout << "\nMasterSub type: " << mintele.Shape() << std::endl;
  std::cout << "MasterSubElement Nodes:";
  for (int k=0;k<mintele.NumNode();++k) std::cout << " " << mintele.NodeIds()[k];  */

  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of master IntElement
  DRT::Element::DiscretizationType mdt = mintele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DQuad called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DQuad called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int nintrow = sintele.NumNode();
  int ndof = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);
  LINALG::SerialDenseVector lmintval(nintrow);
  LINALG::SerialDenseMatrix lmintderiv(nintrow,2,true);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");

  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether displacement shape fct. modification has to be considered or not
  // this is the case for dual quadratic and linear Lagrange multipliers on quad9/quad8/tri6 elements
  bool dualquad3d = false;
  if ( (ShapeFcn() == INPAR::MORTAR::shape_dual) &&
       (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin) &&
       (sele.Shape() == DRT::Element::quad9 || sele.Shape() == DRT::Element::quad8 || sele.Shape() == DRT::Element::tri6) )
  {
    dualquad3d = true;
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // note that the third component of sxi is necessary!
    // (although it will always be 0.0 of course)
    double tempsxi[3] = {0.0, 0.0, 0.0};
    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};
    double projalpha = 0.0;

    // get Gauss point in slave integration element coordinates
    cell->LocalToGlobal(eta,tempsxi,0);
    sxi[0] = tempsxi[0];
    sxi[1] = tempsxi[1];

    // project Gauss point onto master integration element
    MORTAR::MortarProjector projector(3);
    projector.ProjectGaussPoint3D(sintele,sxi,mintele,mxi,projalpha);

    // check GP projection (MASTER)
    double tol = 0.01;
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DQuad: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DQuad: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // map Gauss point back to slave element (affine map)
    // map Gauss point back to master element (affine map)
    double psxi[2] = {0.0, 0.0};
    double pmxi[2] = {0.0, 0.0};
    sintele.MapToParent(sxi,psxi);
    mintele.MapToParent(mxi,pmxi);

    //std::cout << "SInt-GP:    " << sxi[0] << " " << sxi[1] << std::endl;
    //std::cout << "MInt-GP:    " << mxi[0] << " " << mxi[1] << std::endl;
    //std::cout << "SParent-GP: " << psxi[0] << " " << psxi[1] << std::endl;
    //std::cout << "MParent-GP: " << pmxi[0] << " " << pmxi[1] << std::endl;

    // evaluate Lagrange multiplier shape functions (on slave element)
    if (bound)
      sele.EvaluateShapeLagMultLin(ShapeFcn(),psxi,lmval,lmderiv,nrow);
    else
    {
      sele.EvaluateShapeLagMult(ShapeFcn(),psxi,lmval,lmderiv,nrow);
      sintele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmintval,lmintderiv,nintrow);
    }

    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(psxi,sval,sderiv,nrow,dualquad3d);
    mele.EvaluateShape(pmxi,mval,mderiv,ncol);

    // evaluate the two Jacobians (int. cell and slave element)
    double jaccell = cell->Jacobian(eta);
    double jacslave = sintele.Jacobian(sxi);
    double jac = jaccell * jacslave;

    // compute cell D/M matrix *******************************************
    // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
    if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
        (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin))
    {
      // compute all mseg and dseg matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }

        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }

    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    else if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
             lmtype == INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      // compute all mseg (and dseg) matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nintrow*ndof; ++j)
      {
        // integrate mseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);

          // multiply the two shape functions
          double prod = lmintval[jindex]*mval[lindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }

        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmintval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }

    // CASE 4: Dual LM shape functions and quadratic interpolation
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      // compute all mseg (and dseg) matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(l) M_jl
        // hence, they can be combined into one single loop

        // integrate mseg and dseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];

          // isolate the mseg/dseg entries to be filled and
          // add current Gauss point's contribution to mseg/dseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
          {
            (*mseg)(j,l) += prod*jac*wgt;
            (*dseg)(j,j) += prod*jac*wgt;
          }
        }
      }
    }

    // CASE 5: Dual LM shape functions and linear interpolation
    // (here, we must NOT ignore the small off-diagonal terms for accurate convergence)
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual &&
             lmtype == INPAR::MORTAR::lagmult_lin_lin)
    {
      // compute all mseg and dseg matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }

        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }

    // INVALID CASES
    else
    {
      dserror("ERROR: Invalid integration case for 3D quadratic mortar!");
    }
    // compute cell D/M matrix *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Integrate and linearize a 2D slave / master cell (3D)     popp 03/09|
 |  This method integrates the cell M matrix (and possibly D matrix)    |
 |  and stores it in mseg and dseg respectively.                        |
 |  This is the QUADRATIC auxiliary plane coupling version!!!           |
 *----------------------------------------------------------------------*/
void MORTAR::MortarIntegrator::IntegrateDerivCell3DAuxPlaneQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg)
{
  // get LMtype
  INPAR::MORTAR::LagMultQuad lmtype = LagMultQuad();
  
  // explicitely defined shapefunction type needed
  if (ShapeFcn() == INPAR::MORTAR::shape_undefined)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad called without specific shape function defined!");
  
  /*std::cout << std::endl;
  std::cout << "Slave type: " << sele.Shape() << std::endl;
  std::cout << "SlaveElement Nodes:";
  for (int k=0;k<sele.NumNode();++k) std::cout << " " << sele.NodeIds()[k];
  std::cout << "\nMaster type: " << mele.Shape() << std::endl;
  std::cout << "MasterElement Nodes:";
  for (int k=0;k<mele.NumNode();++k) std::cout << " " << mele.NodeIds()[k];
  std::cout << std::endl;
  std::cout << "SlaveSub type: " << sintele.Shape() << std::endl;
  std::cout << "SlaveSubElement Nodes:";
  for (int k=0;k<sintele.NumNode();++k) std::cout << " " << sintele.NodeIds()[k];
  std::cout << "\nMasterSub type: " << mintele.Shape() << std::endl;
  std::cout << "MasterSubElement Nodes:";
  for (int k=0;k<mintele.NumNode();++k) std::cout << " " << mintele.NodeIds()[k];  */

  //check for problem dimension
  if (Dim()!=3) dserror("ERROR: 3D integration method called for non-3D problem");

  // discretization type of slave and master IntElement
  DRT::Element::DiscretizationType sdt = sintele.Shape();
  DRT::Element::DiscretizationType mdt = mintele.Shape();

  // check input data
  if ((!sele.IsSlave()) || (mele.IsSlave()))
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad called on a wrong type of MortarElement pair!");
  if (cell==Teuchos::null)
    dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad called without integration cell");

  // number of nodes (slave, master)
  int nrow = sele.NumNode();
  int ncol = mele.NumNode();
  int nintrow = sintele.NumNode();
  int ndof = static_cast<MORTAR::MortarNode*>(sele.Nodes()[0])->NumDof();

  // create empty vectors for shape fct. evaluation
  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,2,true);
  LINALG::SerialDenseVector mval(ncol);
  LINALG::SerialDenseMatrix mderiv(ncol,2,true);
  LINALG::SerialDenseVector lmval(nrow);
  LINALG::SerialDenseMatrix lmderiv(nrow,2,true);
  LINALG::SerialDenseVector lmintval(nintrow);
  LINALG::SerialDenseMatrix lmintderiv(nintrow,2,true);

  // get slave element nodes themselves
  DRT::Node** mynodes = sele.Nodes();
  if(!mynodes) dserror("ERROR: IntegrateDerivCell3DAuxPlaneQuad: Null pointer!");
  
  // decide whether boundary modification has to be considered or not
  // this is element-specific (is there a boundary node in this element?)
  bool bound = false;
  for (int k=0;k<nrow;++k)
  {
    MORTAR::MortarNode* mymrtrnode = static_cast<MORTAR::MortarNode*>(mynodes[k]);
    if (!mymrtrnode) dserror("ERROR: IntegrateDerivSegment2D: Null pointer!");
    bound += mymrtrnode->IsOnBound();
  }

  // decide whether displacement shape fct. modification has to be considered or not
  // this is the case for dual quadratic and linear Lagrange multipliers on quad9/quad8/tri6 elements
  bool dualquad3d = false;
  if ( (ShapeFcn() == INPAR::MORTAR::shape_dual) &&
       (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin) &&
       (sele.Shape() == DRT::Element::quad9 || sele.Shape() == DRT::Element::quad8 || sele.Shape() == DRT::Element::tri6) )
  {
    dualquad3d = true;
  }

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<nGP();++gp)
  {
    // coordinates and weight
    double eta[2] = {Coordinate(gp,0), Coordinate(gp,1)};
    double wgt = Weight(gp);

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta,globgp,0);

    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};

    // project Gauss point onto slave integration element
    // project Gauss point onto master integration element
    double sprojalpha = 0.0;
    double mprojalpha = 0.0;
    MORTAR::MortarProjector projector(3);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,sintele,sxi,sprojalpha);
    projector.ProjectGaussPointAuxn3D(globgp,auxn,mintele,mxi,mprojalpha);

    // check GP projection (SLAVE)
    double tol = 0.01;
    if (sdt==DRT::Element::quad4 || sdt==DRT::Element::quad8 || sdt==DRT::Element::quad9)
    {
      if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }
    else
    {
      if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Slave Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      }
    }

    // check GP projection (MASTER)
    if (mdt==DRT::Element::quad4 || mdt==DRT::Element::quad8 || mdt==DRT::Element::quad9)
    {
      if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }
    else
    {
      if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
      {
        std::cout << "\n***Warning: IntegrateDerivCell3DAuxPlane: Master Gauss point projection outside!";
        std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
        std::cout << "GP local: " << eta[0] << " " << eta[1] << std::endl;
        std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      }
    }

    // map Gauss point back to slave element (affine map)
    // map Gauss point back to master element (affine map)
    double psxi[2] = {0.0, 0.0};
    double pmxi[2] = {0.0, 0.0};
    sintele.MapToParent(sxi,psxi);
    mintele.MapToParent(mxi,pmxi);

    //std::cout << "SInt-GP:    " << sxi[0] << " " << sxi[1] << std::endl;
    //std::cout << "MInt-GP:    " << mxi[0] << " " << mxi[1] << std::endl;
    //std::cout << "SParent-GP: " << psxi[0] << " " << psxi[1] << std::endl;
    //std::cout << "MParent-GP: " << pmxi[0] << " " << pmxi[1] << std::endl;

    // evaluate Lagrange multiplier shape functions (on slave element)  
    if (bound)
      sele.EvaluateShapeLagMultLin(ShapeFcn(),psxi,lmval,lmderiv,nrow);
    else
    {
      sele.EvaluateShapeLagMult(ShapeFcn(),psxi,lmval,lmderiv,nrow);
      sintele.EvaluateShapeLagMult(ShapeFcn(),sxi,lmintval,lmintderiv,nintrow);
    }
    
    // evaluate trace space shape functions (on both elements)
    sele.EvaluateShape(psxi,sval,sderiv,nrow,dualquad3d);
    mele.EvaluateShape(pmxi,mval,mderiv,ncol);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // compute cell D/M matrix *******************************************
    // CASE 1/2: Standard LM shape functions and quadratic or linear interpolation
    if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
        (lmtype == INPAR::MORTAR::lagmult_quad_quad || lmtype == INPAR::MORTAR::lagmult_lin_lin))
    {
      // compute all mseg and dseg matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);
  
          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }
  
        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }
    
    // CASE 3: Standard LM shape functions and piecewise linear interpolation
    else if (ShapeFcn() == INPAR::MORTAR::shape_standard &&
             lmtype == INPAR::MORTAR::lagmult_pwlin_pwlin)
    {
      // compute all mseg (and dseg) matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nintrow*ndof; ++j)
      {
        // integrate mseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);
  
          // multiply the two shape functions
          double prod = lmintval[jindex]*mval[lindex];
  
          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }
  
        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmintval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }
    
    // CASE 4: Dual LM shape functions and quadratic interpolation
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual &&
             lmtype == INPAR::MORTAR::lagmult_quad_quad)
    {
      // compute all mseg (and dseg) matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // for dual shape functions we can make use
        // of the row summing lemma: D_jj = Sum(l) M_jl
        // hence, they can be combined into one single loop

        // integrate mseg and dseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];

          // isolate the mseg/dseg entries to be filled and
          // add current Gauss point's contribution to mseg/dseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
          {
            (*mseg)(j,l) += prod*jac*wgt;
            (*dseg)(j,j) += prod*jac*wgt;
          }
        }
      }
    }

    // CASE 5: Dual LM shape functions and linear interpolation
    // (here, we must NOT ignore the small off-diagonal terms for accurate convergence)
    else if (ShapeFcn() == INPAR::MORTAR::shape_dual &&
             lmtype == INPAR::MORTAR::lagmult_lin_lin)
    {
      // compute all mseg and dseg matrix entries
      // loop over Lagrange multiplier dofs j
      for (int j=0; j<nrow*ndof; ++j)
      {
        // integrate mseg
        for (int l=0; l<ncol*ndof; ++l)
        {
          int jindex = (int)(j/ndof);
          int lindex = (int)(l/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*mval[lindex];

          // isolate the mseg entries to be filled and
          // add current Gauss point's contribution to mseg
          if ((j==l) || ((j-jindex*ndof)==(l-lindex*ndof)))
            (*mseg)(j,l) += prod*jac*wgt;
        }

        // integrate dseg
        for (int k=0; k<nrow*ndof; ++k)
        {
          int jindex = (int)(j/ndof);
          int kindex = (int)(k/ndof);

          // multiply the two shape functions
          double prod = lmval[jindex]*sval[kindex];

          // isolate the dseg entries to be filled and
          // add current Gauss point's contribution to dseg
          if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
            (*dseg)(j,k) += prod*jac*wgt;
        }
      }
    }
    
    // INVALID CASES
    else
    {
      dserror("ERROR: Invalid integration case for 3D quadratic mortar!");
    }
    // compute cell D/M matrix *******************************************
  }
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave             |
 |  element to the D map of the adjacent slave nodes.                   |
 *----------------------------------------------------------------------*/
bool MORTAR::MortarIntegrator::AssembleD(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                          Epetra_SerialDenseMatrix& dseg)
{
  //problem dimension via numdof of an arbitrary slave node
  DRT::Node** testnodes = sele.Nodes();
  MORTAR::MortarNode* tnode = static_cast<MORTAR::MortarNode*>(testnodes[0]);
  int dim = tnode->NumDof();

  bool nozero=false;
  for (int i=0;i<dim*sele.NumNode();++i)
  {
    for (int j=0;j<dim*sele.NumNode();++j)
    {
      if(dseg(i,j)!=0.0)
      {
        nozero=true;
        break;
      }
    }
  }

  if(nozero==true)
  {
    // get adjacent nodes to assemble to
    DRT::Node** snodes = sele.Nodes();
    if (!snodes)
      dserror("ERROR: AssembleD: Null pointer for snodes!");

    // loop over all slave nodes
    for (int slave=0;slave<sele.NumNode();++slave)
    {
      MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
      int sndof = snode->NumDof();

      // only process slave node rows that belong to this proc
      if (snode->Owner() != comm.MyPID())
        continue;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (snode->IsOnBound())
        continue;

      // loop over all dofs of the slave node
      for (int sdof=0;sdof<sndof;++sdof)
      {
        // loop over all slave nodes again ("master nodes")
        for (int master=0;master<sele.NumNode();++master)
        {
          MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(snodes[master]);
          const int* mdofs = mnode->Dofs();
          int mndof = mnode->NumDof();

          // loop over all dofs of the slave node again ("master dofs")
          for (int mdof=0;mdof<mndof;++mdof)
          {
            int col = mdofs[mdof];
            double val = dseg(slave*sndof+sdof,master*mndof+mdof);

            // BOUNDARY NODE MODIFICATION **********************************
            // The boundary nodes have already been re-defined  as being
            // master nodes, so all we have to do here is to shift the terms
            // associated with a boundary node from D to the respective
            // place in M! Mind the minus sign, which is needed because M
            // will later enter the contact forces with another minus sign.
            // *************************************************************
            if (mnode->IsOnBound())
            {
              double minusval = -val;
              snode->AddMValue(sdof,col,minusval);
            }
            else
            {
              snode->AddDValue(sdof,col,val);
            }
          }
        }
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution (2D / 3D)                         popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool MORTAR::MortarIntegrator::AssembleD(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                         MORTAR::IntElement& sintele,
                                         Epetra_SerialDenseMatrix& dseg)
{
  // get adjacent int nodes to assemble to
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("ERROR: AssembleD: Null pointer for sintnodes!");
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleD: Null pointer for snodes!");

  // loop over all slave int nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    MORTAR::MortarNode* sintnode = static_cast<MORTAR::MortarNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all slave nodes ("master nodes")
      for (int master=0;master<sele.NumNode();++master)
      {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(snodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the slave node ("master dofs")
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = dseg(slave*sintndof+sdof,master*mndof+mdof);

          // assembly
          sintnode->AddDValue(sdof,col,val);
        }
      }
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution (2D / 3D)                         popp 01/08|
 |  This method assembles the contrubution of a 1D/2D slave and master  |
 |  overlap pair to the M map of the adjacent slave nodes.              |
 *----------------------------------------------------------------------*/
bool MORTAR::MortarIntegrator::AssembleM(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                         MORTAR::MortarElement& mele,
                                         Epetra_SerialDenseMatrix& mseg)
{
  // get adjacent slave nodes and master nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleM: Null pointer for snodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleM: Null pointer for mnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
    int sndof = snode->NumDof();

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mseg(slave*sndof+sdof,master*mndof+mdof);
          snode->AddMValue(sdof,col,val);
        }
      }
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution (2D / 3D)                         popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool MORTAR::MortarIntegrator::AssembleM(const Epetra_Comm& comm,
                                         MORTAR::IntElement& sintele,
                                         MORTAR::MortarElement& mele,
                                         Epetra_SerialDenseMatrix& mseg)
{
  // get adjacent slave int nodes and master nodes
  DRT::Node** sintnodes = sintele.Nodes();
  if (!sintnodes) dserror("ERROR: AssembleM: Null pointer for sintnodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes) dserror("ERROR: AssembleM: Null pointer for mnodes!");

  // loop over all slave int nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    MORTAR::MortarNode* sintnode = static_cast<MORTAR::MortarNode*>(sintnodes[slave]);
    int sintndof = sintnode->NumDof();

    // only process slave int node rows that belong to this proc
    if (sintnode->Owner() != comm.MyPID()) continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sintndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mseg(slave*sintndof+sdof,master*mndof+mdof);
          sintnode->AddMValue(sdof,col,val);
        }
      }
    }
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution for ele. b. int. (2D / 3D)       farah 01/13|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool MORTAR::MortarIntegrator::AssembleM_EleBased(const Epetra_Comm& comm,
                                         MORTAR::MortarElement& sele,
                                         std::vector<MORTAR::MortarElement*> meles,
                                         Epetra_SerialDenseMatrix& mseg)
{
  //problem dimension via numdof of an arbitrary slave node
  DRT::Node** testnodes = sele.Nodes();
  MORTAR::MortarNode* tnode = static_cast<MORTAR::MortarNode*>(testnodes[0]);
  int dim = tnode->NumDof();

  for ( int nummasterele=0;nummasterele<(int)meles.size();++nummasterele )
  {
    //if all entries zero-valued --> no assembly
    bool nozero=false;
    for ( int i=0;i<dim*meles[nummasterele]->NumNode();++i  )
    {
      for (int j=0;j<dim*sele.NumNode();++j)
      {
        if(mseg(j,i+nummasterele*meles[nummasterele]->NumNode()*dim) != 0.0 )
        {
          nozero=true;
          break;
        }
      }
    }

    if (nozero==true)
    {
      // get adjacent slave nodes and master nodes
      DRT::Node** snodes = sele.Nodes();
      if (!snodes)
      dserror("ERROR: AssembleM: Null pointer for snodes!");
      DRT::Node** mnodes = meles[nummasterele]->Nodes();
      if (!mnodes)
      dserror("ERROR: AssembleM: Null pointer for mnodes!");

      // loop over all slave nodes
      for (int slave=0;slave<sele.NumNode();++slave)
      {
      MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
      int sndof = snode->NumDof();

      // only process slave node rows that belong to this proc
      if (snode->Owner() != comm.MyPID())
        continue;

      // do not process slave side boundary nodes
      // (their row entries would be zero anyway!)
      if (snode->IsOnBound())
        continue;

      // loop over all dofs of the slave node
      for (int sdof=0;sdof<sndof;++sdof)
      {
        // loop over all master nodes
        for (int master=0;master<meles[nummasterele]->NumNode();++master)
        {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs(); //global dofs
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];

          double val = mseg(slave*sndof+sdof,master*mndof+mdof+nummasterele*mndof*meles[nummasterele]->NumNode());
          snode->AddMValue(sdof,col,val);
        }
        }
      }
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble Mmod contribution (2D)                            popp 01/08|
 |  This method assembles the contribution of a 1D slave / master        |
 |  overlap pair to the Mmod map of the adjacent slave nodes.            |
 *----------------------------------------------------------------------*/
bool MORTAR::MortarIntegrator::AssembleMmod(const Epetra_Comm& comm,
                                            MORTAR::MortarElement& sele,
                                            MORTAR::MortarElement& mele,
                                            Epetra_SerialDenseMatrix& mmodseg)
{
  //**********************************************************************
  dserror("ERROR: AssembleMmod method is outdated!");
  //**********************************************************************

  // get adjacent slave nodes and master nodes
  DRT::Node** snodes = sele.Nodes();
  if (!snodes)
    dserror("ERROR: AssembleMmod: Null pointer for snodes!");
  DRT::Node** mnodes = mele.Nodes();
  if (!mnodes)
    dserror("ERROR: AssembleMmod: Null pointer for mnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    MORTAR::MortarNode* snode = static_cast<MORTAR::MortarNode*>(snodes[slave]);
    //const int* sdofs = snode->Dofs();
    int sndof = snode->NumDof();

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID())
      continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound())
      continue;

    // loop over all dofs of the slave node
    for (int sdof=0;sdof<sndof;++sdof)
    {
      // loop over all master nodes
      for (int master=0;master<mele.NumNode();++master)
      {
        MORTAR::MortarNode* mnode = static_cast<MORTAR::MortarNode*>(mnodes[master]);
        const int* mdofs = mnode->Dofs();
        int mndof = mnode->NumDof();

        // loop over all dofs of the master node
        for (int mdof=0;mdof<mndof;++mdof)
        {
          int col = mdofs[mdof];
          double val = mmodseg(slave*sndof+sdof,master*mndof+mdof);
          snode->AddMmodValue(sdof,col,val);
        }
      }
    }
    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << std::endl;
    std::map<int, double> nodemap0 = (snode->GetMmod())[0];
    std::map<int, double> nodemap1 = (snode->GetMmod())[1];
    typedef std::map<int,double>::const_iterator CI;

    std::cout << "Row dof id: " << sdofs[0] << std::endl;;
    for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
      std::cout << p->first << '\t' << p->second << std::endl;

    std::cout << "Row dof id: " << sdofs[1] << std::endl;
    for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
      std::cout << p->first << '\t' << p->second << std::endl;
#endif // #ifdef DEBUG
     */
  }

  return true;
}

