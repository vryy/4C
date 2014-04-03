/*!----------------------------------------------------------------------
\file volmortar_integrator.cpp

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "volmortar_integrator.H"
#include "volmortar_shape.H"
#include "volmortar_defines.H"
#include "volmortar_cell.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_mortar/mortar_calc_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::VolMortarIntegrator()
{
  InitializeGP();
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points                                  farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::InitializeGP(bool integrateele, int domain)
{
  DRT::Element::DiscretizationType intshape = DRT::Element::dis_none;

  if(integrateele)
  {
    if(domain==0)
      intshape=distypeS;
    else if(domain==1)
      intshape=distypeM;
    else
      dserror("integration domain not specified!");
  }
  else
  {
    if(ndim_==2)
      intshape=DRT::Element::tri3;
    else if(ndim_==3)
      intshape=DRT::Element::tet4;
    else
      dserror("wrong dimension!");
  }


  //*******************************
  // choose Gauss rule according
  //*******************************
  switch(intshape)
  {
  case DRT::Element::tri3:
  {
    DRT::UTILS::GaussRule2D mygaussrule=DRT::UTILS::intrule_tri_7point;

    const DRT::UTILS::IntegrationPoints2D intpoints(mygaussrule);
    ngp_ = intpoints.nquad;
    coords_.Reshape(ngp_,2);
    weights_.resize(ngp_);
    for (int i=0;i<ngp_;++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=intpoints.qxg[i][1];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  case DRT::Element::tet4:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_tet_24point;

    const DRT::UTILS::IntegrationPoints3D intpoints(mygaussrule);
    ngp_ = intpoints.nquad;
    coords_.Reshape(ngp_,3);
    weights_.resize(ngp_);
    for (int i=0;i<ngp_;++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=intpoints.qxg[i][1];
      coords_(i,2)=intpoints.qxg[i][2];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  case DRT::Element::hex8:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_hex_27point;

    const DRT::UTILS::IntegrationPoints3D intpoints(mygaussrule);
    ngp_ = intpoints.nquad;
    coords_.Reshape(ngp_,3);
    weights_.resize(ngp_);
    for (int i=0;i<ngp_;++i)
    {
      coords_(i,0)=intpoints.qxg[i][0];
      coords_(i,1)=intpoints.qxg[i][1];
      coords_(i,2)=intpoints.qxg[i][2];
      weights_[i]=intpoints.qwgt[i];
    }
    break;
  }
  default:
  {
    dserror("ERROR: VolMortarIntegrator: This element type is not implemented!");
    break;
  }
  } // switch(eletype)

  return;
}

/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::IntegrateCells2D(
     DRT::Element& sele,
     DRT::Element& mele,
     Teuchos::RCP<MORTAR::IntCell> cell,
     LINALG::SparseMatrix& dmatrix,
     LINALG::SparseMatrix& mmatrix,
     Teuchos::RCP<const DRT::Discretization> slavedis,
     Teuchos::RCP<const DRT::Discretization> masterdis)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1>             sval;
  LINALG::Matrix<nm_,1>             mval;
  LINALG::Matrix<ns_,1>             lmval;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<ngp_;++gp)
  {
//    // coordinates and weight
    double eta[2] = {coords_(gp,0), coords_(gp,1)};
    double wgt = weights_[gp];

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta,globgp,0);

    // map gp into slave and master para space
    double sxi[2] = {0.0, 0.0};
    double mxi[2] = {0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(sele,globgp,sxi);
    MORTAR::UTILS::GlobalToLocal<distypeM>(mele,globgp,mxi);

    // Check parameter space mapping
    bool proj = CheckMapping2D(sele,mele,sxi,mxi);
    if (proj==false)
      dserror("Mapping failed!");

    // evaluate trace space shape functions (on both elements)
    UTILS::volmortar_shape_function_2D(sval, sxi[0],sxi[1],distypeS);
    UTILS::volmortar_shape_function_2D(mval, mxi[0],mxi[1],distypeM);

    // evaluate Lagrange mutliplier shape functions (on slave element)
    //UTILS::volmortar_shape_function_2D(lmval, sxi[0],sxi[1],distypeS);
    UTILS::volmortar_dualshape_function_2D(lmval,sele, sxi[0],sxi[1],distypeS);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // compute segment D/M matrix ****************************************
    // standard shape functions
    if (false)//(shapefcn_ == INPAR::MORTAR::shape_standard)
    {
      for (int j=0; j<ns_; ++j)
      {
        DRT::Node* cnode = sele.Nodes()[j];
        int nsdof=slavedis->NumDof(1,cnode);

        //loop over slave dofs
        for (int jdof=0;jdof<nsdof;++jdof)
        {
          int row = slavedis->Dof(1,cnode,jdof);

          ////////////////////////////////////////
          // integrate M
          for (int k=0; k<nm_; ++k)
          {
            DRT::Node* mnode = mele.Nodes()[k];
            int nmdof=masterdis->NumDof(0,mnode);

            for (int kdof=0;kdof<nmdof;++kdof)
            {
              int col = masterdis->Dof(0,mnode,kdof);

              // multiply the two shape functions
              double prod = lmval(j)*mval(k)*jac*wgt;

              // dof to dof
              if (jdof==kdof)
              {
                if(abs(prod)>VOLMORTARINTTOL) mmatrix.Assemble(prod, row, col);
              }
            }
          }

          ////////////////////////////////////////
          // integrate D
          for (int k=0; k<ns_; ++k)
          {
            DRT::Node* snode = sele.Nodes()[k];
            int nddof=slavedis->NumDof(1,snode);

            for (int kdof=0;kdof<nddof;++kdof)
            {
              int col = slavedis->Dof(1,snode,kdof);

              // multiply the two shape functions
              double prod = lmval(j)*sval(k)*jac*wgt;

              // dof to dof
              if (jdof==kdof)
              {
                if(abs(prod)>VOLMORTARINTTOL) dmatrix.Assemble(prod, row, col);
              }
            }
          }
        }
      }
    }
    else  // DUAL
    {
      for (int j=0;j<ns_;++j)
      {
        DRT::Node* cnode = sele.Nodes()[j];
        int nsdof=slavedis->NumDof(1,cnode);

        //loop over slave dofs
        for (int jdof=0;jdof<nsdof;++jdof)
        {
          int row = slavedis->Dof(1,cnode,jdof);

          ////////////////////////////////////////////////////////////////
          // integrate M and D
          for (int k=0; k<nm_; ++k)
          {
            DRT::Node* mnode = mele.Nodes()[k];
            int nmdof=masterdis->NumDof(0,mnode);

            for (int kdof=0;kdof<nmdof;++kdof)
            {
              int col = masterdis->Dof(0,mnode,kdof);

              // multiply the two shape functions
              double prod = lmval(j)*mval(k)*jac*wgt;

              // dof to dof
              if (jdof==kdof)
              {
                if (abs(prod)>VOLMORTARINTTOL) mmatrix.Assemble(prod, row, col);
                if (abs(prod)>VOLMORTARINTTOL) dmatrix.Assemble(prod, row, row);
              }
            }
          }
        ////////////////////////////////////////////////////////////////
        }
      }
    }
  }//end gp loop

  return;
}

/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::IntegrateCells3D(
    DRT::Element& Aele,
    DRT::Element& Bele,
    Teuchos::RCP<VOLMORTAR::Cell> cell,
    LINALG::SparseMatrix& dmatrix_A,
    LINALG::SparseMatrix& mmatrix_A,
    LINALG::SparseMatrix& dmatrix_B,
    LINALG::SparseMatrix& mmatrix_B,
    Teuchos::RCP<const DRT::Discretization> Adis,
    Teuchos::RCP<const DRT::Discretization> Bdis)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1>             sval_A;
  LINALG::Matrix<nm_,1>             mval_A;
  LINALG::Matrix<ns_,1>             lmval_A;

  LINALG::Matrix<nm_,1>             lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<ngp_;++gp)
  {
    // coordinates and weight
    double eta[3] = {coords_(gp,0), coords_(gp,1), coords_(gp,2)};
    double wgt = weights_[gp];

    // evaluate the integration cell Jacobian
    double jac = cell->Vol();

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta,globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    double Bxi[3] = {0.0, 0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(Aele,globgp,Axi);
    MORTAR::UTILS::GlobalToLocal<distypeM>(Bele,globgp,Bxi);

    // Check parameter space mapping
    CheckMapping3D(Aele,Bele,Axi,Bxi);
//    std::cout << "GP= " << globgp[0] << "  " << globgp[1] << "  " << globgp[2] << std::endl;
//    cell->Print();
//    bool proj = CheckMapping3D(Aele,Bele,Axi,Bxi);
//    if (proj==false)
//      jac=0.0;

    //---------------------------------------------------
    //---------------------------------------------------
    //---------------------------------------------------
    //---------------------------------------------------

    // evaluate trace space shape functions (on both elements)
    UTILS::volmortar_shape_function_3D(sval_A, Axi[0],Axi[1], Axi[2],distypeS);
    UTILS::volmortar_shape_function_3D(mval_A, Bxi[0],Bxi[1], Bxi[2],distypeM);

    // evaluate Lagrange mutliplier shape functions (on slave element)
    UTILS::volmortar_dualshape_function_3D(lmval_A,Aele, Axi[0],Axi[1],Axi[2],distypeS);
    // ---
    UTILS::volmortar_dualshape_function_3D(lmval_B,Bele, Bxi[0],Bxi[1],Bxi[2],distypeM);


    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j=0;j<ns_;++j)
    {
      DRT::Node* cnode = Aele.Nodes()[j];
      int nsdof=Adis->NumDof(1,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Adis->Dof(1,cnode,jdof);

        // integrate M and D
        for (int k=0; k<nm_; ++k)
        {
          DRT::Node* mnode = Bele.Nodes()[k];
          int nmdof=Bdis->NumDof(0,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Bdis->Dof(0,mnode,kdof);

            // multiply the two shape functions
            double prod = lmval_A(j)*mval_A(k)*jac*wgt;

            // dof to dof
            if (jdof==kdof)
            {
              if (abs(prod)>VOLMORTARINTTOL) mmatrix_A.Assemble(prod, row, col);
              if (abs(prod)>VOLMORTARINTTOL) dmatrix_A.Assemble(prod, row, row);
            }
          }
        }
      }
    }

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j=0;j<nm_;++j)
    {
      DRT::Node* cnode = Bele.Nodes()[j];
      int nsdof=Bdis->NumDof(1,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Bdis->Dof(1,cnode,jdof);

        // integrate M and D
        for (int k=0; k<ns_; ++k)
        {
          DRT::Node* mnode = Aele.Nodes()[k];
          int nmdof=Adis->NumDof(0,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Adis->Dof(0,mnode,kdof);

            // multiply the two shape functions
            double prod = lmval_B(j)*sval_A(k)*jac*wgt;

            // dof to dof
            if (jdof==kdof)
            {
              if (abs(prod)>VOLMORTARINTTOL) mmatrix_B.Assemble(prod, row, col);
              if (abs(prod)>VOLMORTARINTTOL) dmatrix_B.Assemble(prod, row, row);
            }
          }
        }
      }
    }

  }//end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::IntegrateEle3D(
    int domain,
     DRT::Element& Aele,
     DRT::Element& Bele,
     LINALG::SparseMatrix& dmatrix_A,
     LINALG::SparseMatrix& mmatrix_A,
     LINALG::SparseMatrix& dmatrix_B,
     LINALG::SparseMatrix& mmatrix_B,
     Teuchos::RCP<const DRT::Discretization> Adis,
     Teuchos::RCP<const DRT::Discretization> Bdis)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1>             sval_A;
  LINALG::Matrix<nm_,1>             mval_A;
  LINALG::Matrix<ns_,1>             lmval_A;

  LINALG::Matrix<nm_,1>             lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<ngp_;++gp)
  {
//    // coordinates and weight
    double eta[3] = {coords_(gp,0), coords_(gp,1), coords_(gp,2)};
    double wgt = weights_[gp];
    double jac = 0.0;
    double globgp[3] = {0.0, 0.0, 0.0};

    if(domain==0)
    {
      // evaluate the integration cell Jacobian
      jac = UTILS::Jacobian<distypeS>(eta,Aele);

      // get global Gauss point coordinates
      UTILS::LocalToGlobal<distypeS>(Aele,eta,globgp);
    }
    else if(domain==1)
    {
      // evaluate the integration cell Jacobian
      jac = UTILS::Jacobian<distypeM>(eta,Bele);

      // get global Gauss point coordinates
      UTILS::LocalToGlobal<distypeM>(Bele,eta,globgp);
    }
    else
      dserror("wrong domain for integration!");


    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    double Bxi[3] = {0.0, 0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(Aele,globgp,Axi);
    MORTAR::UTILS::GlobalToLocal<distypeM>(Bele,globgp,Bxi);

    // Check parameter space mapping
    CheckMapping3D(Aele,Bele,Axi,Bxi);
//    std::cout << "==================================" << std::endl;
//    if(domain==0)
//    {
//      for(int i=0;i<Aele.NumNode();++i)
//      {
//        std::cout << "create vertex " << Aele.Nodes()[i]->X()[0] <<"  "<< Aele.Nodes()[i]->X()[1] <<"  "<< Aele.Nodes()[i]->X()[2] <<std::endl;
//      }
//    }
//    if(domain==1)
//    {
//      for(int i=0;i<Bele.NumNode();++i)
//      {
//        std::cout << "create vertex " << Bele.Nodes()[i]->X()[0] <<"  "<< Bele.Nodes()[i]->X()[1] <<"  "<< Bele.Nodes()[i]->X()[2] <<std::endl;
//      }
//    }

//    bool proj = CheckMapping3D(Aele,Bele,Axi,Bxi);
//    if (proj==false)
//      jac=0.0;

    //---------------------------------------------------
    //---------------------------------------------------
    //---------------------------------------------------
    //---------------------------------------------------

    // evaluate trace space shape functions (on both elements)
    UTILS::volmortar_shape_function_3D(sval_A, Axi[0],Axi[1], Axi[2],distypeS);
    UTILS::volmortar_shape_function_3D(mval_A, Bxi[0],Bxi[1], Bxi[2],distypeM);

    // evaluate Lagrange mutliplier shape functions (on slave element)
    UTILS::volmortar_dualshape_function_3D(lmval_A,Aele, Axi[0],Axi[1],Axi[2],distypeS);
    // ---
    UTILS::volmortar_dualshape_function_3D(lmval_B,Bele, Bxi[0],Bxi[1],Bxi[2],distypeM);


    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j=0;j<ns_;++j)
    {
      DRT::Node* cnode = Aele.Nodes()[j];
      int nsdof=Adis->NumDof(1,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Adis->Dof(1,cnode,jdof);

        // integrate M and D
        for (int k=0; k<nm_; ++k)
        {
          DRT::Node* mnode = Bele.Nodes()[k];
          int nmdof=Bdis->NumDof(0,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Bdis->Dof(0,mnode,kdof);

            // multiply the two shape functions
            double prod = lmval_A(j)*mval_A(k)*jac*wgt;

            // dof to dof
            if (jdof==kdof)
            {
              if (abs(prod)>VOLMORTARINTTOL) mmatrix_A.Assemble(prod, row, col);
              if (abs(prod)>VOLMORTARINTTOL) dmatrix_A.Assemble(prod, row, row);
            }
          }
        }
      }
    }

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j=0;j<nm_;++j)
    {
      DRT::Node* cnode = Bele.Nodes()[j];
      int nsdof=Bdis->NumDof(1,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Bdis->Dof(1,cnode,jdof);

        // integrate M and D
        for (int k=0; k<ns_; ++k)
        {
          DRT::Node* mnode = Aele.Nodes()[k];
          int nmdof=Adis->NumDof(0,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Adis->Dof(0,mnode,kdof);

            // multiply the two shape functions
            double prod = lmval_B(j)*sval_A(k)*jac*wgt;

            // dof to dof
            if (jdof==kdof)
            {
              if (abs(prod)>VOLMORTARINTTOL) mmatrix_B.Assemble(prod, row, col);
              if (abs(prod)>VOLMORTARINTTOL) dmatrix_B.Assemble(prod, row, row);
            }
          }
        }
      }
    }

  }//end gp loop

  return;
}

/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
bool VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::CheckMapping2D(DRT::Element& sele,
                                                                     DRT::Element& mele,
                                                                     double* sxi, double* mxi)
{
  // check GP projection (SLAVE)
  double tol = 0.01;
  if (distypeS==DRT::Element::quad4 || distypeS==DRT::Element::quad8 || distypeS==DRT::Element::quad9)
  {
    if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      return false;
    }
  }
  else
  {
    if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      return false;
    }
  }

  // check GP projection (MASTER)
  if (distypeM==DRT::Element::quad4 || distypeM==DRT::Element::quad8 || distypeM==DRT::Element::quad9)
  {
    if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      return false;
    }
  }
  else
  {
    if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      return false;
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
bool VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::CheckMapping3D(DRT::Element& sele,
                                                                       DRT::Element& mele,
                                                                       double* sxi, double* mxi)
{
  // check GP projection (SLAVE)
  double tol = 1e-5;
  if (distypeS==DRT::Element::hex8 || distypeS==DRT::Element::hex20 || distypeS==DRT::Element::hex27)
  {
    if (sxi[0]<-1.0-tol || sxi[1]<-1.0-tol || sxi[2]<-1.0-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[2]>1.0+tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] << std::endl;

      for(int i=0;i<sele.NumNode();++i)
      {
        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      }
      std::cout << "------------" << std::endl;
      for(int i=0;i<mele.NumNode();++i)
      {
        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      }

      return false;
    }
  }
  else if(distypeS==DRT::Element::tet4 || distypeS==DRT::Element::tet10)
  {
    if(sxi[0]<0.0-tol || sxi[1]<0.0-tol || sxi[2]<0.0-tol || (sxi[0]+sxi[1]+sxi[2])>1.0+tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Master GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] << std::endl;
      for(int i=0;i<sele.NumNode();++i)
      {
        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      }
      std::cout << "------------" << std::endl;
      for(int i=0;i<mele.NumNode();++i)
      {
        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      }
      return false;
    }
  }
  else
    dserror("Wrong element type!");

  // check GP projection (MASTER)
  if (distypeM==DRT::Element::hex8 || distypeM==DRT::Element::hex20 || distypeM==DRT::Element::hex27)
  {
    if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[2]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[2]>1.0+tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] << std::endl;
      for(int i=0;i<sele.NumNode();++i)
      {
        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      }
      std::cout << "------------" << std::endl;
      for(int i=0;i<mele.NumNode();++i)
      {
        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      }
      return false;
    }
  }
  else if(distypeM==DRT::Element::tet4 || distypeM==DRT::Element::tet10)
  {
    if(mxi[0]<0.0-tol || mxi[1]<0.0-tol || mxi[2]<0.0-tol || (mxi[0]+mxi[1]+mxi[2])>1.0+tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] << std::endl;
      for(int i=0;i<sele.NumNode();++i)
      {
        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
      }
      std::cout << "------------" << std::endl;
      for(int i=0;i<mele.NumNode();++i)
      {
        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
      }
      return false;
    }
  }
  else
    dserror("Wrong element type!");

  return true;
}
/*----------------------------------------------------------------------*
 |  possible slave/master elment pairs                       farah 01/14|
 *----------------------------------------------------------------------*/
//slave quad4
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::quad4,DRT::Element::quad4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::quad4,DRT::Element::tri3>;

//slave tri3
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tri3,DRT::Element::quad4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tri3,DRT::Element::tri3>;

//slave hex8
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::tet4>;

//slave tet4
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::tet4>;
