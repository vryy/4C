/*!----------------------------------------------------------------------
\file volmortar_integrator.cpp

\brief integration routines for the volmortar framework

\level 1

<pre>
\maintainer Philipp Farah
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
#include "../drt_cut/cut_volumecell.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 02/15|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS>
VOLMORTAR::VolMortarIntegratorEleBased<distypeS>::VolMortarIntegratorEleBased(Teuchos::ParameterList& params)
{
  // get type of quadratic modification
  dualquad_ = DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::DualQuad>(params,"DUALQUAD");
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points for ele-based integration        farah 02/15|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS>
void VOLMORTAR::VolMortarIntegratorEleBased<distypeS>::InitializeGP()
{
  // init shape of integration domain
  DRT::Element::DiscretizationType intshape = distypeS;

  //*******************************
  // choose Gauss rule accordingly
  //*******************************
  switch(intshape)
  {
  //*******************************
  //               2D
  //*******************************
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
  case DRT::Element::tri6:
  {
    DRT::UTILS::GaussRule2D mygaussrule=DRT::UTILS::intrule_tri_12point;

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
  case DRT::Element::quad4:
  {
    DRT::UTILS::GaussRule2D mygaussrule=DRT::UTILS::intrule_quad_64point;

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
  case DRT::Element::quad8:
  {
    DRT::UTILS::GaussRule2D mygaussrule=DRT::UTILS::intrule_quad_64point;

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
  case DRT::Element::quad9:
  {
    DRT::UTILS::GaussRule2D mygaussrule=DRT::UTILS::intrule_quad_64point;

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
  //*******************************
  //               3D
  //*******************************
  case DRT::Element::tet4:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_tet_45point;

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
  case DRT::Element::tet10:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_tet_45point;

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
  case DRT::Element::hex20:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_hex_125point;

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
  case DRT::Element::hex27:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_hex_125point;

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
  case DRT::Element::pyramid5:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_pyramid_8point;

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
  //*******************************
  //            Default
  //*******************************
  default:
  {
    dserror("ERROR: VolMortarIntegrator: This element type is not implemented!");
    break;
  }
  } // switch(eletype)

  return;
}

/*----------------------------------------------------------------------*
 |  Initialize gauss points for ele-based integration        farah 02/15|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS>
void VOLMORTAR::VolMortarIntegratorEleBased<distypeS>::IntegrateEleBased3D(
    DRT::Element& sele,
    std::vector<int>& foundeles,
    LINALG::SparseMatrix& D,
    LINALG::SparseMatrix& M,
    Teuchos::RCP<const DRT::Discretization> Adis,
    Teuchos::RCP<const DRT::Discretization> Bdis,
    int dofseta,
    int dofsetb,
    const Teuchos::RCP<const Epetra_Map>& PAB_dofrowmap,
    const Teuchos::RCP<const Epetra_Map>& PAB_dofcolmap)
{
  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<ngp_;++gp)
  {
    // coordinates and weight
    double eta[3]    = {coords_(gp,0), coords_(gp,1), 0.0};

    if(ndim_==3)
      eta[2] = coords_(gp,2);

    double wgt       = weights_[gp];
    double jac       = 0.0;
    double globgp[3] = {0.0, 0.0, 0.0};

    // quantities for eval. outside gp
    double gpdist    = 1.0e12;
    int    gpid      = 0;
    double AuxXi[3]  = {0.0, 0.0, 0.0};

    // evaluate the integration cell Jacobian
    jac = UTILS::Jacobian<distypeS>(eta,sele);

    // get global Gauss point coordinates
    UTILS::LocalToGlobal<distypeS>(sele,eta,globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(sele,globgp,Axi);

    // loop over beles
    for(int found=0;found<(int)foundeles.size();++found)
    {
      //get master element
      DRT::Element* Bele = Bdis->gElement(foundeles[found]);
      DRT::Element::DiscretizationType shape = Bele->Shape();

      bool proj = false;

      switch (shape)
      {
      //************************************************
      //                    2D
      //************************************************
      case DRT::Element::tri3:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::tri3>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::tri6:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::tri6>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::quad4:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::quad4>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::quad8:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::quad8>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::quad9:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::quad9>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);
        break;
      }
      //************************************************
      //                    3D
      //************************************************
      case DRT::Element::hex8:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::hex8>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::hex20:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::hex20>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::hex27:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::hex27>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::tet4:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::tet4>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::tet10:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::tet10>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      case DRT::Element::pyramid5:
      {
        proj=VolMortarEleBasedGP<distypeS,DRT::Element::pyramid5>(
            sele,
            Bele,
            foundeles,
            found,
            gpid,
            jac,
            wgt,
            gpdist,
            Axi,
            AuxXi,
            globgp,
            dualquad_,
            D,
            M,
            Adis,
            Bdis,
            dofseta,
            dofsetb,
            PAB_dofrowmap,
            PAB_dofcolmap);

        break;
      }
      default:
      {
        dserror("ERROR: unknown shape!");
        break;
      }
      }

      // if gp evaluated break ele loop
      if(proj == true)
        break;
      else
        continue;
    }//beles
  }//end gp loop
  return;
}


/*----------------------------------------------------------------------*
 |  possible elements for ele-based integration              farah 02/15|
 *----------------------------------------------------------------------*/
template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::quad4>;
template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::quad8>;
template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::quad9>;

template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::tri3>;
template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::tri6>;

template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::hex20>;
template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::hex27>;

template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::tet4>;
template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::tet10>;

template class VOLMORTAR::VolMortarIntegratorEleBased<DRT::Element::pyramid5>;


/*----------------------------------------------------------------------*
 |  gp evaluation                                            farah 02/15|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
bool VOLMORTAR::VolMortarEleBasedGP(
    DRT::Element& sele,
    DRT::Element* mele,
    std::vector<int>& foundeles,
    int& found, int& gpid,
    double& jac,
    double& wgt,
    double& gpdist,
    double* Axi,
    double* AuxXi,
    double* globgp,
    INPAR::VOLMORTAR::DualQuad& dq,
    LINALG::SparseMatrix& D,
    LINALG::SparseMatrix& M,
    Teuchos::RCP<const DRT::Discretization> Adis,
    Teuchos::RCP<const DRT::Discretization> Bdis,
    int dofseta,
    int dofsetb,
    const Teuchos::RCP<const Epetra_Map>& PAB_dofrowmap,
    const Teuchos::RCP<const Epetra_Map>& PAB_dofcolmap)
{
  //! ns_: number of slave element nodes
  static const int ns_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement;

  //! nm_: number of master element nodes
  static const int nm_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;

  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1> sval_A;
  LINALG::Matrix<nm_,1> mval_A;
  LINALG::Matrix<ns_,1> lmval_A;

  double Bxi[3] = {0.0, 0.0, 0.0};

  bool converged = true;
  MORTAR::UTILS::GlobalToLocal<distypeM>(*mele,globgp,Bxi, converged);
  if(!converged and found!=((int)foundeles.size()-1))
    return false;

  // save distance of gp
  double l = sqrt(Bxi[0]*Bxi[0] + Bxi[1]*Bxi[1] + Bxi[2]*Bxi[2]);
  if(l<gpdist)
  {
    gpdist = l;
    gpid   = foundeles[found];
    AuxXi[0]  = Bxi[0];
    AuxXi[1]  = Bxi[1];
    AuxXi[2]  = Bxi[2];
  }

  // Check parameter space mapping
  bool proj = CheckMapping<distypeS,distypeM>(sele,*mele,Axi,Bxi);

  // if gp outside continue or eval nearest gp
  if(!proj and (found!=((int)foundeles.size()-1)))
    return false;
  else if(!proj and found==((int)foundeles.size()-1))
  {
    Bxi[0] = AuxXi[0];
    Bxi[1] = AuxXi[1];
    Bxi[2] = AuxXi[2];
    mele = Bdis->gElement(gpid);
  }

  // for "master" side
  UTILS::shape_function<distypeS>(sval_A,Axi,dq);
  UTILS::shape_function<distypeM>(mval_A,Bxi);

  // evaluate Lagrange multiplier shape functions (on slave element)
  UTILS::dual_shape_function<distypeS>(lmval_A,Axi,sele,dq);

  // compute cell D/M matrix ****************************************
  // dual shape functions
  for (int j=0;j<ns_;++j)
  {
    DRT::Node* cnode = sele.Nodes()[j];
    if (cnode->Owner() != Adis->Comm().MyPID())
      continue;

    const int nsdof=Adis->NumDof(dofseta,cnode);

//    if (true)//(shapefcn_ == INPAR::MORTAR::shape_standard)
//    {
//      for (int j=0; j<ns_; ++j)
//      {
//        DRT::Node* cnode = sele.Nodes()[j];
//        int nsdof=Adis->NumDof(dofseta,cnode);
//
//        //loop over slave dofs
//        for (int jdof=0;jdof<nsdof;++jdof)
//        {
//          int row = Adis->Dof(dofseta,cnode,jdof);
//
//          ////////////////////////////////////////
//          // integrate M
//          for (int k=0; k<nm_; ++k)
//          {
//            DRT::Node* mnode = mele->Nodes()[k];
//            int nmdof=Bdis->NumDof(dofsetb,mnode);
//
//            for (int kdof=0;kdof<nmdof;++kdof)
//            {
//              int col = Bdis->Dof(dofsetb,mnode,kdof);
//
//              // multiply the two shape functions
//              double prod = sval_A(j)*mval_A(k)*jac*wgt;
//
//              // dof to dof
//              if (jdof==kdof)
//              {
//                if(abs(prod)>VOLMORTARINTTOL) M.Assemble(prod, row, col);
//              }
//            }
//          }
//
//          ////////////////////////////////////////
//          // integrate D
//          for (int k=0; k<ns_; ++k)
//          {
//            DRT::Node* snode = sele.Nodes()[k];
//            int nddof=Adis->NumDof(dofseta,snode);
//
//            for (int kdof=0;kdof<nddof;++kdof)
//            {
//              int col = Adis->Dof(dofseta,snode,kdof);
//
//              // multiply the two shape functions
//              double prod = sval_A(j)*sval_A(k)*jac*wgt;
//
//              // dof to dof
//              if (jdof==kdof)
//              {
//                if(abs(prod)>VOLMORTARINTTOL) D.Assemble(prod, row, col);
//              }
//            }
//          }
//        }
//      }
//    }

    //loop over slave dofs
    for (int jdof=0;jdof<nsdof;++jdof)
    {
      const int row = Adis->Dof(dofseta,cnode,jdof);

      if(not PAB_dofrowmap->MyGID(row))
        continue;

      // integrate D
      const double prod2 = lmval_A(j)*sval_A(j)*jac*wgt;
      if (abs(prod2)>VOLMORTARINTTOL)
        D.Assemble(prod2, row, row);

      // integrate M
      for (int k=0; k<nm_; ++k)
      {
        DRT::Node* mnode = mele->Nodes()[k];
        const int col = Bdis->Dof(dofsetb,mnode,jdof);

        if(not PAB_dofcolmap->MyGID(col))
          continue;

        // multiply the two shape functions
        const double prod = lmval_A(j)*mval_A(k)*jac*wgt;

        if (abs(prod)>VOLMORTARINTTOL)
          M.Assemble(prod, row, col);
      }
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  possible slave/master element pairs                      farah 02/15|
 *----------------------------------------------------------------------*/
////slave quad4
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::quad4,DRT::Element::quad4>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::quad4,DRT::Element::tri3>;
//
////slave tri3
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tri3,DRT::Element::quad4>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tri3,DRT::Element::tri3>;
//
////slave hex8
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex8,DRT::Element::tet4>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex8,DRT::Element::tet10>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex8,DRT::Element::hex8>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex8,DRT::Element::hex27>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex8,DRT::Element::hex20>;
//
////slave hex20
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex20,DRT::Element::tet4>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex20,DRT::Element::tet10>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex20,DRT::Element::hex8>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex20,DRT::Element::hex27>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex20,DRT::Element::hex20>;
//
////slave hex27
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex27,DRT::Element::tet4>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex27,DRT::Element::tet10>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex27,DRT::Element::hex8>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex27,DRT::Element::hex27>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::hex27,DRT::Element::hex20>;
//
////slave tet4
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet4,DRT::Element::tet4>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet4,DRT::Element::tet10>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet4,DRT::Element::hex8>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet4,DRT::Element::hex27>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet4,DRT::Element::hex20>;
//
////slave tet10
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet10,DRT::Element::tet4>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet10,DRT::Element::tet10>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet10,DRT::Element::hex8>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet10,DRT::Element::hex27>;
//template class VOLMORTAR::VolMortarEleBasedGP<DRT::Element::tet10,DRT::Element::hex20>;


/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::VolMortarIntegrator(Teuchos::ParameterList& params)
{
  // get type of quadratic modification
  dualquad_ = DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::DualQuad>(params,"DUALQUAD");

  // define gp rule
  InitializeGP();
}


/*----------------------------------------------------------------------*
 |  Initialize gauss points                                  farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::InitializeGP(bool integrateele, int domain, DRT::Element::DiscretizationType shape)
{
  // init shape of integration domain
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
      intshape=shape;
    else
      dserror("wrong dimension!");
  }

  //*******************************
  // choose Gauss rule accordingly
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
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_tet_45point;

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
  case DRT::Element::tet10:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_tet_45point;

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
  case DRT::Element::hex20:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_hex_125point;

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
  case DRT::Element::hex27:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_hex_125point;

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
  case DRT::Element::pyramid5:
  {
    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_pyramid_8point;

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
     Teuchos::RCP<const DRT::Discretization> masterdis,
     int sdofset,
     int mdofset)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1> sval;
  LINALG::Matrix<nm_,1> mval;
  LINALG::Matrix<ns_,1> lmval;

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
    double sxi[3] = {0.0, 0.0 , 0.0};
    double mxi[3] = {0.0, 0.0 , 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(sele,globgp,sxi);
    MORTAR::UTILS::GlobalToLocal<distypeM>(mele,globgp,mxi);

    // Check parameter space mapping
    bool proj = CheckMapping2D(sele,mele,sxi,mxi);
    if (proj==false)
      dserror("ERROR: Mapping failed!");

    // evaluate trace space shape functions (on both elements)
    UTILS::shape_function<distypeS>(sval,sxi);
    UTILS::shape_function<distypeM>(mval,mxi);

    // evaluate Lagrange mutliplier shape functions (on slave element)
    //UTILS::volmortar_shape_function_2D(lmval, sxi[0],sxi[1],distypeS);
    UTILS::dual_shape_function<distypeS>(lmval,sxi,sele);

    // evaluate the integration cell Jacobian
    double jac = cell->Jacobian(eta);

    // compute segment D/M matrix ****************************************
    // standard shape functions
    if (false)//(shapefcn_ == INPAR::MORTAR::shape_standard)
    {
      for (int j=0; j<ns_; ++j)
      {
        DRT::Node* cnode = sele.Nodes()[j];
        int nsdof=slavedis->NumDof(sdofset,cnode);

        //loop over slave dofs
        for (int jdof=0;jdof<nsdof;++jdof)
        {
          int row = slavedis->Dof(sdofset,cnode,jdof);

          ////////////////////////////////////////
          // integrate M
          for (int k=0; k<nm_; ++k)
          {
            DRT::Node* mnode = mele.Nodes()[k];
            int nmdof=masterdis->NumDof(mdofset,mnode);

            for (int kdof=0;kdof<nmdof;++kdof)
            {
              int col = masterdis->Dof(mdofset,mnode,kdof);

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
            int nddof=slavedis->NumDof(sdofset,snode);

            for (int kdof=0;kdof<nddof;++kdof)
            {
              int col = slavedis->Dof(sdofset,snode,kdof);

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

        if (cnode->Owner() != slavedis->Comm().MyPID())
          continue;

        int nsdof=slavedis->NumDof(sdofset,cnode);

        //loop over slave dofs
        for (int jdof=0;jdof<nsdof;++jdof)
        {
          int row = slavedis->Dof(sdofset,cnode,jdof);

          ////////////////////////////////////////////////////////////////
          // integrate M and D
          for (int k=0; k<nm_; ++k)
          {
            DRT::Node* mnode = mele.Nodes()[k];
            int nmdof=masterdis->NumDof(mdofset,mnode);

            for (int kdof=0;kdof<nmdof;++kdof)
            {
              int col = masterdis->Dof(mdofset,mnode,kdof);

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
    Teuchos::RCP<const DRT::Discretization> Bdis,
    int sdofset_A,
    int mdofset_A,
    int sdofset_B,
    int mdofset_B)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1> sval_A;
  LINALG::Matrix<nm_,1> mval_A;
  LINALG::Matrix<ns_,1> lmval_A;
  LINALG::Matrix<nm_,1> lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<ngp_;++gp)
  {
    // coordinates and weight
    double eta[3] = {coords_(gp,0), coords_(gp,1), coords_(gp,2)};
    double wgt = weights_[gp];

    // get global Gauss point coordinates
    double globgp[3] = {0.0, 0.0, 0.0};
    cell->LocalToGlobal(eta,globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    double Bxi[3] = {0.0, 0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(Aele,globgp,Axi);
    MORTAR::UTILS::GlobalToLocal<distypeM>(Bele,globgp,Bxi);

    // evaluate the integration cell Jacobian
    double jac=0.0;
    if (cell->Shape()==DRT::Element::tet4)
      jac      = cell->Vol();
    else if(cell->Shape()==DRT::Element::hex8)
      jac = cell->CalcJac(eta);
    else
      dserror("used shape not supported in volmortar integrator!");

    // Check parameter space mapping
    //std::cout << "globgp " << globgp[0] <<"  "<< globgp[1] <<"  "<< globgp[2] <<std::endl;

    bool check = CheckMapping3D(Aele,Bele,Axi,Bxi);
    if(!check)
      continue;

    // evaluate trace space shape functions (on both elements)
    UTILS::shape_function<distypeS>(sval_A, Axi);
    UTILS::shape_function<distypeM>(mval_A, Bxi);

    // evaluate Lagrange multiplier shape functions (on slave element)
    UTILS::dual_shape_function<distypeS>(lmval_A,Axi,Aele,dualquad_);
    UTILS::dual_shape_function<distypeM>(lmval_B,Bxi,Bele,dualquad_);

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j=0;j<ns_;++j)
    {
      DRT::Node* cnode = Aele.Nodes()[j];
      int nsdof=Adis->NumDof(sdofset_A,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Adis->Dof(sdofset_A,cnode,jdof);

        // integrate M and D
        for (int k=0; k<nm_; ++k)
        {
          DRT::Node* mnode = Bele.Nodes()[k];
          int nmdof=Bdis->NumDof(mdofset_A,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Bdis->Dof(mdofset_A,mnode,kdof);

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
      int nsdof=Bdis->NumDof(sdofset_B,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Bdis->Dof(sdofset_B,cnode,jdof);

        // integrate M and D
        for (int k=0; k<ns_; ++k)
        {
          DRT::Node* mnode = Aele.Nodes()[k];
          int nmdof=Adis->NumDof(mdofset_B,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Adis->Dof(mdofset_B,mnode,kdof);

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
 |  Compute D/M entries for Volumetric Mortar                farah 04/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::IntegrateCells3D_DirectDiveregence(
    DRT::Element& Aele,
    DRT::Element& Bele,
    GEO::CUT::VolumeCell& vc,
    Teuchos::RCP<DRT::UTILS::GaussPoints> intpoints,
    bool switched_conf,
    LINALG::SparseMatrix& dmatrix_A,
    LINALG::SparseMatrix& mmatrix_A,
    LINALG::SparseMatrix& dmatrix_B,
    LINALG::SparseMatrix& mmatrix_B,
    Teuchos::RCP<const DRT::Discretization> Adis,
    Teuchos::RCP<const DRT::Discretization> Bdis,
    int sdofset_A,
    int mdofset_A,
    int sdofset_B,
    int mdofset_B)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1> sval_A;
  LINALG::Matrix<nm_,1> mval_A;
  LINALG::Matrix<ns_,1> lmval_A;
  LINALG::Matrix<nm_,1> lmval_B;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<intpoints->NumPoints();++gp)
  {
    double weight_out = intpoints->Weight(gp);
    // coordinates and weight
    double eta[3] = {intpoints->Point(gp)[0], intpoints->Point(gp)[1], intpoints->Point(gp)[2]};

    double globgp[3] = {0.0, 0.0, 0.0};

    if(switched_conf)
      UTILS::LocalToGlobal<distypeS>(Aele,eta,globgp);
    else
      UTILS::LocalToGlobal<distypeM>(Bele,eta,globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    double Bxi[3] = {0.0, 0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(Aele,globgp,Axi);
    MORTAR::UTILS::GlobalToLocal<distypeM>(Bele,globgp,Bxi);

//      std::cout << "-------------------------------------" << std::endl;
//      std::cout << "globgp= " << globgp[0] << "  " << globgp[1] << "  " << globgp[2] << std::endl;
//      std::cout << "eta= " << eta[0] << "  " << eta[1] << "  " << eta[2] << std::endl;
//      std::cout << "Axi= " << Axi[0] << "  " << Axi[1] << "  " << Axi[2] << std::endl;
//      std::cout << "Bxi= " << Bxi[0] << "  " << Bxi[1] << "  " << Bxi[2] << std::endl;

    // evaluate the integration cell Jacobian
    double jac = 0.0;

    if(switched_conf)
      jac    = UTILS::Jacobian<distypeS>(Axi,Aele);
    else
      jac    = UTILS::Jacobian<distypeM>(Bxi,Bele);

    // Check parameter space mapping
    //CheckMapping3D(Aele,Bele,Axi,Bxi);

    // evaluate trace space shape functions (on both elements)
    UTILS::shape_function<distypeS>(sval_A, Axi);
    UTILS::shape_function<distypeM>(mval_A, Bxi);

    // evaluate Lagrange multiplier shape functions (on slave element)
    UTILS::dual_shape_function<distypeS>(lmval_A, Axi,Aele,dualquad_);
    UTILS::dual_shape_function<distypeM>(lmval_B, Bxi,Bele,dualquad_);

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j=0;j<ns_;++j)
    {
      DRT::Node* cnode = Aele.Nodes()[j];
      int nsdof=Adis->NumDof(sdofset_A,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Adis->Dof(sdofset_A,cnode,jdof);

        // integrate M and D
        for (int k=0; k<nm_; ++k)
        {
          DRT::Node* mnode = Bele.Nodes()[k];
          int nmdof=Bdis->NumDof(mdofset_A,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Bdis->Dof(mdofset_A,mnode,kdof);

            // multiply the two shape functions
            double prod = lmval_A(j)*mval_A(k)*jac*weight_out;
//              std::cout << "PROD1 = " << prod  << " row= " << row << "  col= " << col << "  j= " << j<<  "  nsdof= " << nsdof<< std::endl;
//              cnode->Print(cout);
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
      int nsdof=Bdis->NumDof(sdofset_B,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Bdis->Dof(sdofset_B,cnode,jdof);

        // integrate M and D
        for (int k=0; k<ns_; ++k)
        {
          DRT::Node* mnode = Aele.Nodes()[k];
          int nmdof=Adis->NumDof(mdofset_B,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Adis->Dof(sdofset_B,mnode,kdof);

            // multiply the two shape functions
            double prod = lmval_B(j)*sval_A(k)*jac*weight_out;
//              std::cout << "PROD2 = " << prod  << " row= " << row << "  col= " << col<< std::endl;
//              cnode->Print(cout);

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
 |  Compute D/M entries for Volumetric Mortar                farah 04/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::IntegrateEleBased3D_ADis(DRT::Element& Aele,
    std::vector<int>& foundeles,
    LINALG::SparseMatrix& dmatrix_A,
    LINALG::SparseMatrix& mmatrix_A,
    Teuchos::RCP<const DRT::Discretization> Adis,
    Teuchos::RCP<const DRT::Discretization> Bdis,
    int dofsetA,
    int dofsetB)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1> sval_A;
  LINALG::Matrix<nm_,1> mval_A;
  LINALG::Matrix<ns_,1> lmval_A;

  //**********************************************************************
  // loop over all Gauss points for integration
  //**********************************************************************
  for (int gp=0;gp<ngp_;++gp)
  {
    // coordinates and weight
    double eta[3]    = {coords_(gp,0), coords_(gp,1), coords_(gp,2)};
    double wgt       = weights_[gp];
    double jac       = 0.0;
    double globgp[3] = {0.0, 0.0, 0.0};

    // quantities for eval. outside gp
    double gpdist    = 1.0e12;
    int    gpid      = 0;
    double AuxXi[3]  = {0.0, 0.0, 0.0};

    // evaluate the integration cell Jacobian
    jac = UTILS::Jacobian<distypeS>(eta,Aele);

    // get global Gauss point coordinates
    UTILS::LocalToGlobal<distypeS>(Aele,eta,globgp);

    // map gp into A and B para space
    double Axi[3] = {0.0, 0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeS>(Aele,globgp,Axi);

    // loop over beles
    for(int found=0;found<(int)foundeles.size();++found)
    {
      //get master element
      DRT::Element* Bele = Bdis->gElement(foundeles[found]);
      double Bxi[3] = {0.0, 0.0, 0.0};

      bool converged = true;
      MORTAR::UTILS::GlobalToLocal<distypeM>(*Bele,globgp,Bxi, converged);
      if(!converged and found!=((int)foundeles.size()-1))
        continue;

      // save distance of gp
      double l = sqrt(Bxi[0]*Bxi[0] + Bxi[1]*Bxi[1] + Bxi[2]*Bxi[2]);
      if(l<gpdist)
      {
        gpdist = l;
        gpid   = foundeles[found];
        AuxXi[0]  = Bxi[0];
        AuxXi[1]  = Bxi[1];
        AuxXi[2]  = Bxi[2];
      }

      // Check parameter space mapping
      bool proj = CheckMapping3D(Aele,*Bele,Axi,Bxi);

      // if gp outside continue or eval nearest gp
      if(!proj and (found!=((int)foundeles.size()-1)))
        continue;
      else if(!proj and found==((int)foundeles.size()-1))
      {
        Bxi[0] = AuxXi[0];
        Bxi[1] = AuxXi[1];
        Bxi[2] = AuxXi[2];
        Bele = Bdis->gElement(gpid);
      }

      // for "master" side
      UTILS::shape_function<distypeS>(sval_A,Axi,dualquad_);
      UTILS::shape_function<distypeM>(mval_A,Bxi);

      // evaluate Lagrange multiplier shape functions (on slave element)
      UTILS::dual_shape_function<distypeS>(lmval_A,Axi,Aele,dualquad_);

      // compute cell D/M matrix ****************************************
      // dual shape functions
      for (int j=0;j<ns_;++j)
      {
        DRT::Node* cnode = Aele.Nodes()[j];
        if (cnode->Owner() != Adis->Comm().MyPID())
          continue;

        int nsdof=Adis->NumDof(dofsetA,cnode);

        //loop over slave dofs
        for (int jdof=0;jdof<nsdof;++jdof)
        {
          int row = Adis->Dof(dofsetA,cnode,jdof);

          // integrate D
          double prod2 = lmval_A(j)*sval_A(j)*jac*wgt;
          if (abs(prod2)>VOLMORTARINTTOL)
            dmatrix_A.Assemble(prod2, row, row);

          // integrate M
          for (int k=0; k<nm_; ++k)
          {
            DRT::Node* mnode = Bele->Nodes()[k];
            int col = Bdis->Dof(dofsetB,mnode,jdof);

            // multiply the two shape functions
            double prod = lmval_A(j)*mval_A(k)*jac*wgt;

            if (abs(prod)>VOLMORTARINTTOL)
              mmatrix_A.Assemble(prod, row, col);
          }
        }
      }

      break;
    }//beles
  }//end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 04/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void VOLMORTAR::VolMortarIntegrator<distypeS,distypeM>::IntegrateEleBased3D_BDis(DRT::Element& Bele,
    std::vector<int>& foundeles,
    LINALG::SparseMatrix& dmatrix_B,
    LINALG::SparseMatrix& mmatrix_B,
    Teuchos::RCP<const DRT::Discretization> Adis,
    Teuchos::RCP<const DRT::Discretization> Bdis,
    int dofsetA,
    int dofsetB)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1> mval_A;
  LINALG::Matrix<nm_,1> sval_B;
  LINALG::Matrix<nm_,1> lmval_B;

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

    // quantities for eval. outside gp
    double gpdist    = 1.0e12;
    int    gpid      = 0;
    double AuxXi[3]  = {0.0, 0.0, 0.0};

    // evaluate the integration cell Jacobian
    jac = UTILS::Jacobian<distypeM>(eta,Bele);

    // get global Gauss point coordinates
    UTILS::LocalToGlobal<distypeM>(Bele,eta,globgp);

    // map gp into A and B para space
    double Bxi[3] = {0.0, 0.0, 0.0};
    MORTAR::UTILS::GlobalToLocal<distypeM>(Bele,globgp,Bxi);

      // loop over beles
    for(int found=0;found<(int)foundeles.size();++found)
    {
        //get master element
        DRT::Element* Aele = Adis->gElement(foundeles[found]);
      double Axi[3] = {0.0, 0.0, 0.0};

      bool converged = true;
      MORTAR::UTILS::GlobalToLocal<distypeS>(*Aele,globgp,Axi,converged);
      if(!converged and found!=((int)foundeles.size()-1))
        continue;

      // save distance of gp
      double l = sqrt(Axi[0]*Axi[0] + Axi[1]*Axi[1] + Axi[2]*Axi[2]);
      if(l<gpdist)
      {
        gpdist = l;
        gpid=foundeles[found];
        AuxXi[0]  = Axi[0];
        AuxXi[1]  = Axi[1];
        AuxXi[2]  = Axi[2];
      }

      // Check parameter space mapping
      bool proj = CheckMapping3D(*Aele,Bele,Axi,Bxi);

      // if gp outside continue or eval nearest gp
      if(!proj and (found!=((int)foundeles.size()-1)))
        continue;
      else if(!proj and found==((int)foundeles.size()-1))
      {
        Axi[0] = AuxXi[0];
        Axi[1] = AuxXi[1];
        Axi[2] = AuxXi[2];
        Aele = Adis->gElement(gpid);
      }

      // evaluate trace space shape functions (on both elements)
      UTILS::shape_function<distypeM>(sval_B,Bxi,dualquad_);
      UTILS::shape_function<distypeS>(mval_A,Axi);

      // evaluate Lagrange multiplier shape functions (on slave element)
      UTILS::dual_shape_function<distypeM>(lmval_B,Bxi,Bele,dualquad_);
      // compute cell D/M matrix ****************************************
      // dual shape functions
      for (int j=0;j<nm_;++j)
      {
        DRT::Node* cnode = Bele.Nodes()[j];
        if (cnode->Owner() != Bdis->Comm().MyPID())
          continue;

        int nsdof=Bdis->NumDof(dofsetB,cnode);

        //loop over slave dofs
        for (int jdof=0;jdof<nsdof;++jdof)
        {
          int row = Bdis->Dof(dofsetB,cnode,jdof);

          // integrate D
          double prod2 = lmval_B(j)*sval_B(j)*jac*wgt;
          if (abs(prod2)>VOLMORTARINTTOL)
            dmatrix_B.Assemble(prod2, row, row);

          // integrate M
          for (int k=0; k<ns_; ++k)
          {
            DRT::Node* mnode = Aele->Nodes()[k];
            int col = Adis->Dof(dofsetA,mnode,jdof);

            // multiply the two shape functions
            double prod = lmval_B(j)*mval_A(k)*jac*wgt;

            if (abs(prod)>VOLMORTARINTTOL)
              mmatrix_B.Assemble(prod, row, col);
          }
        }
      }

      break;
    }//beles
  }//end gp loop

  return;
}


/*----------------------------------------------------------------------*
 |  Compute D/M entries for Volumetric Mortar                farah 01/14|
 |  This function is for element-wise integration when an               |
 |  element is completely located within an other element               |
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
     Teuchos::RCP<const DRT::Discretization> Bdis,
     int sdofset_A,
     int mdofset_A,
     int sdofset_B,
     int mdofset_B)
{
  // create empty vectors for shape fct. evaluation
  LINALG::Matrix<ns_,1> sval_A;
  LINALG::Matrix<nm_,1> mval_A;
  LINALG::Matrix<ns_,1> lmval_A;
  LINALG::Matrix<nm_,1> lmval_B;

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

    // evaluate trace space shape functions (on both elements)
    UTILS::shape_function<distypeS>(sval_A,Axi);
    UTILS::shape_function<distypeM>(mval_A,Bxi);

    // evaluate Lagrange multiplier shape functions (on slave element)
    UTILS::dual_shape_function<distypeS>(lmval_A,Axi,Aele,dualquad_);
    UTILS::dual_shape_function<distypeM>(lmval_B,Bxi,Bele,dualquad_);

    // compute cell D/M matrix ****************************************
    // dual shape functions
    for (int j=0;j<ns_;++j)
    {
      DRT::Node* cnode = Aele.Nodes()[j];
      int nsdof=Adis->NumDof(sdofset_A,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Adis->Dof(sdofset_A,cnode,jdof);

        // integrate M and D
        for (int k=0; k<nm_; ++k)
        {
          DRT::Node* mnode = Bele.Nodes()[k];
          int nmdof=Bdis->NumDof(mdofset_A,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Bdis->Dof(mdofset_A,mnode,kdof);

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
      int nsdof=Bdis->NumDof(sdofset_B,cnode);

      //loop over slave dofs
      for (int jdof=0;jdof<nsdof;++jdof)
      {
        int row = Bdis->Dof(sdofset_B,cnode,jdof);

        // integrate M and D
        for (int k=0; k<ns_; ++k)
        {
          DRT::Node* mnode = Aele.Nodes()[k];
          int nmdof=Adis->NumDof(mdofset_B,mnode);

          for (int kdof=0;kdof<nmdof;++kdof)
          {
            int col = Adis->Dof(mdofset_B,mnode,kdof);

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
  const double tol = 1e-10;
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
  else if (distypeS==DRT::Element::tri3 || distypeS==DRT::Element::tri6)
  {
    if (sxi[0]<-tol || sxi[1]<-tol || sxi[0]>1.0+tol || sxi[1]>1.0+tol || sxi[0]+sxi[1]>1.0+2*tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << std::endl;
      return false;
    }
  }
  else
    dserror("Wrong element type!");

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
  else if (distypeS==DRT::Element::tri3 || distypeS==DRT::Element::tri6)
  {
    if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
    {
      std::cout << "\n***Warning: Gauss point projection outside!";
      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << std::endl;
      return false;
    }
  }
  else
    dserror("Wrong element type!");

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
//      std::cout << "\n***Warning: Gauss point projection outside!";
//      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
//      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] << std::endl;
//
//      for(int i=0;i<sele.NumNode();++i)
//      {
//        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
//      }
//      std::cout << "------------" << std::endl;
//      for(int i=0;i<mele.NumNode();++i)
//      {
//        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
//      }

      return false;
    }
  }
  else if(distypeS==DRT::Element::tet4 || distypeS==DRT::Element::tet10)
  {
    if(sxi[0]<0.0-tol || sxi[1]<0.0-tol || sxi[2]<0.0-tol || (sxi[0]+sxi[1]+sxi[2])>1.0+tol)
    {
//      std::cout << "\n***Warning: Gauss point projection outside!";
//      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
//      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] << std::endl;
//      for(int i=0;i<sele.NumNode();++i)
//      {
//        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
//      }
//      std::cout << "------------" << std::endl;
//      for(int i=0;i<mele.NumNode();++i)
//      {
//        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
//      }
      return false;
    }
  }
  else if(distypeS==DRT::Element::pyramid5)
  {
    if (sxi[2]<0.0-tol || -sxi[0]+sxi[2]>1.0+tol || sxi[0]+sxi[2]>1.0+tol || -sxi[1]+sxi[2]>1.0+tol || sxi[1]+sxi[2]>1.0+tol)
    {
//      std::cout << "\n***Warning: Gauss point projection outside!";
//      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
//      std::cout << "Slave GP projection: " << sxi[0] << " " << sxi[1] << " " << sxi[2] << std::endl;
//      for(int i=0;i<sele.NumNode();++i)
//      {
//        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
//      }
//      std::cout << "------------" << std::endl;
//      for(int i=0;i<mele.NumNode();++i)
//      {
//        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
//      }
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
//      std::cout << "\n***Warning: Gauss point projection outside!";
//      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
//      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] << std::endl;
//      for(int i=0;i<sele.NumNode();++i)
//      {
//        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
//      }
//      std::cout << "------------" << std::endl;
//      for(int i=0;i<mele.NumNode();++i)
//      {
//        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
//      }
      return false;
    }
  }
  else if(distypeM==DRT::Element::tet4 || distypeM==DRT::Element::tet10)
  {
    if(mxi[0]<0.0-tol || mxi[1]<0.0-tol || mxi[2]<0.0-tol || (mxi[0]+mxi[1]+mxi[2])>1.0+tol)
    {
//      std::cout << "\n***Warning: Gauss point projection outside!";
//      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
//      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] << std::endl;
//      for(int i=0;i<sele.NumNode();++i)
//      {
//        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
//      }
//      std::cout << "------------" << std::endl;
//      for(int i=0;i<mele.NumNode();++i)
//      {
//        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
//      }
      return false;
    }
  }
  else if(distypeM==DRT::Element::pyramid5)
  {
    if (mxi[2]<0.0-tol || -mxi[0]+mxi[2]>1.0+tol || mxi[0]+mxi[2]>1.0+tol || -mxi[1]+mxi[2]>1.0+tol || mxi[1]+mxi[2]>1.0+tol)
    {
//      std::cout << "\n***Warning: Gauss point projection outside!";
//      std::cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << std::endl;
//      std::cout << "Master GP projection: " << mxi[0] << " " << mxi[1] << " " << mxi[2] << std::endl;
//      for(int i=0;i<sele.NumNode();++i)
//      {
//        std::cout << "create vertex " << sele.Nodes()[i]->X()[0] <<"  "<< sele.Nodes()[i]->X()[1] <<"  "<< sele.Nodes()[i]->X()[2] <<std::endl;
//      }
//      std::cout << "------------" << std::endl;
//      for(int i=0;i<mele.NumNode();++i)
//      {
//        std::cout << "create vertex " << mele.Nodes()[i]->X()[0] <<"  "<< mele.Nodes()[i]->X()[1] <<"  "<< mele.Nodes()[i]->X()[2] <<std::endl;
//      }
      return false;
    }
  }
  else
    dserror("Wrong element type!");

  return true;
}


/*----------------------------------------------------------------------*
 |  possible slave/master element pairs                       farah 01/14|
 *----------------------------------------------------------------------*/
//slave quad4
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::quad4,DRT::Element::quad4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::quad4,DRT::Element::tri3>;

//slave tri3
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tri3,DRT::Element::quad4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tri3,DRT::Element::tri3>;

//slave hex8
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::tet4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::tet10>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::hex27>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::hex20>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex8,DRT::Element::pyramid5>;

//slave hex20
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex20,DRT::Element::tet4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex20,DRT::Element::tet10>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex20,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex20,DRT::Element::hex27>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex20,DRT::Element::hex20>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex20,DRT::Element::pyramid5>;

//slave hex27
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex27,DRT::Element::tet4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex27,DRT::Element::tet10>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex27,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex27,DRT::Element::hex27>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex27,DRT::Element::hex20>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::hex27,DRT::Element::pyramid5>;

//slave tet4
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::tet4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::tet10>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::hex27>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::hex20>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet4,DRT::Element::pyramid5>;

//slave tet10
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet10,DRT::Element::tet4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet10,DRT::Element::tet10>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet10,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet10,DRT::Element::hex27>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet10,DRT::Element::hex20>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::tet10,DRT::Element::pyramid5>;

//slave pyramid 5
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::pyramid5,DRT::Element::tet4>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::pyramid5,DRT::Element::tet10>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::pyramid5,DRT::Element::hex8>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::pyramid5,DRT::Element::hex27>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::pyramid5,DRT::Element::hex20>;
template class VOLMORTAR::VolMortarIntegrator<DRT::Element::pyramid5,DRT::Element::pyramid5>;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 06/14|
 *----------------------------------------------------------------------*/
VOLMORTAR::ConsInterpolator::ConsInterpolator()
{
  // empty
}


/*----------------------------------------------------------------------*
 |  interpolate (public)                                     farah 06/14|
 *----------------------------------------------------------------------*/
void VOLMORTAR::ConsInterpolator::Interpolate(
    DRT::Node* node,
    LINALG::SparseMatrix& pmatrix,
    Teuchos::RCP<const DRT::Discretization> nodediscret,
    Teuchos::RCP<const DRT::Discretization> elediscret,
    std::vector<int>& foundeles,
    std::pair<int,int>& dofset,
    const Teuchos::RCP<const Epetra_Map>& P_dofrowmap,
    const Teuchos::RCP<const Epetra_Map>& P_dofcolmap)
{
  // check ownership
  if (node->Owner() != nodediscret->Comm().MyPID())
    return;

  // map gp into A and B para space
  double nodepos[3] = {node->X()[0], node->X()[1], node->X()[2]};
  double dist       = 1.0e12;
  int    eleid      = 0;
  double AuxXi[3]   = {0.0, 0.0, 0.0};

  // element loop (brute force)
  for(int found=0;found<(int)foundeles.size();++found)
  {
    bool proj = false;

    //get master element
    DRT::Element* ele = elediscret->gElement(foundeles[found]);

    switch (ele->Shape())
    {
    // 2D --------------------------------------------
    case DRT::Element::tri3:
    {
      proj = ConsInterpolatorEval<DRT::Element::tri3>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::tri6:
    {
      proj = ConsInterpolatorEval<DRT::Element::tri6>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::quad4:
    {
      proj = ConsInterpolatorEval<DRT::Element::quad4>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::quad8:
    {
      proj = ConsInterpolatorEval<DRT::Element::quad8>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::quad9:
    {
      proj = ConsInterpolatorEval<DRT::Element::quad9>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }

    // 3D --------------------------------------------
    case DRT::Element::hex8:
    {
      proj = ConsInterpolatorEval<DRT::Element::hex8>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::hex20:
    {
      proj = ConsInterpolatorEval<DRT::Element::hex20>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::hex27:
    {
      proj = ConsInterpolatorEval<DRT::Element::hex27>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::tet4:
    {
      proj = ConsInterpolatorEval<DRT::Element::tet4>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::tet10:
    {
      proj = ConsInterpolatorEval<DRT::Element::tet10>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    case DRT::Element::pyramid5:
    {
      proj = ConsInterpolatorEval<DRT::Element::pyramid5>(
          node,
          ele,
          pmatrix,
          nodediscret,
          elediscret,
          foundeles,
          found,
          eleid,
          dist,
          AuxXi,
          nodepos,
          dofset,
          P_dofrowmap,
          P_dofcolmap);
      break;
    }
    default:
    {
      dserror("ERROR: unknown shape!");
      break;
    }
    } // end switch

    // if node evaluated break ele loop
    if(proj == true)
      break;
    else
      continue;

    break;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  node evaluation                                          farah 02/15|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool VOLMORTAR::ConsInterpolatorEval(
    DRT::Node* node,
    DRT::Element* ele,
    LINALG::SparseMatrix& pmatrix,
    Teuchos::RCP<const DRT::Discretization> nodediscret,
    Teuchos::RCP<const DRT::Discretization> elediscret,
    std::vector<int>& foundeles,
    int& found,
    int& eleid,
    double& dist,
    double* AuxXi,
    double* nodepos,
    std::pair<int,int>& dofset,
    const Teuchos::RCP<const Epetra_Map>& P_dofrowmap,
    const Teuchos::RCP<const Epetra_Map>& P_dofcolmap)
{
  //! ns_: number of slave element nodes
  static const int n_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  double xi[3]   = {0.0, 0.0, 0.0};
  bool converged = true;

  MORTAR::UTILS::GlobalToLocal<distype>(*ele,nodepos,xi,converged);

  // no convergence of local newton?
  if(!converged and found!=((int)foundeles.size()-1))
    return false;

  // save distance of gp
  double l = sqrt(xi[0]*xi[0] + xi[1]*xi[1] + xi[2]*xi[2]);
  if(l<dist)
  {
    dist = l;
    eleid = foundeles[found];
    AuxXi[0] = xi[0];
    AuxXi[1] = xi[1];
    AuxXi[2] = xi[2];
  }

  // Check parameter space mapping
  bool proj = CheckMapping<distype>(*ele,xi);

  // if node outside --> continue or eval nearest gp
  if(!proj and found!=((int)foundeles.size()-1))
  {
    return false;
  }
  else if(!proj and found==((int)foundeles.size()-1))
  {
    xi[0] = AuxXi[0];
    xi[1] = AuxXi[1];
    xi[2] = AuxXi[2];
    ele = elediscret->gElement(eleid);
  }

  // get values
  LINALG::Matrix<n_,1>  val;
  UTILS::shape_function<distype>(val,xi);

  int nsdof=nodediscret->NumDof(dofset.first,node);

  //loop over slave dofs
  for (int jdof=0;jdof<nsdof;++jdof)
  {
    const int row = nodediscret->Dof(dofset.first,node,jdof);

    if(not P_dofrowmap->MyGID(row))
      continue;

    for (int k=0; k<ele->NumNode(); ++k)
    {
      DRT::Node* bnode = ele->Nodes()[k];
      const int col = elediscret->Dof(dofset.second,bnode,jdof);

      if(not P_dofcolmap->MyGID(col))
        continue;

      const double val2=val(k);
      //if (abs(val2)>VOLMORTARINTTOL)
        pmatrix.Assemble(val2, row, col);
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  possible elements for interpolation                      farah 06/14|
 *----------------------------------------------------------------------*/
//template class VOLMORTAR::ConsInterpolator<DRT::Element::quad4>;
//template class VOLMORTAR::ConsInterpolator<DRT::Element::quad8>;
//template class VOLMORTAR::ConsInterpolator<DRT::Element::quad9>;
//
//template class VOLMORTAR::ConsInterpolator<DRT::Element::tri3>;
//template class VOLMORTAR::ConsInterpolator<DRT::Element::tri6>;
//
//template class VOLMORTAR::ConsInterpolator<DRT::Element::hex8>;
//template class VOLMORTAR::ConsInterpolator<DRT::Element::hex20>;
//template class VOLMORTAR::ConsInterpolator<DRT::Element::hex27>;
//
//template class VOLMORTAR::ConsInterpolator<DRT::Element::tet4>;
//template class VOLMORTAR::ConsInterpolator<DRT::Element::tet10>;

