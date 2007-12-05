/*----------------------------------------------------------------------*/
/*!
\file xfluid3_stationary.cpp

\brief Internal implementation of XFluid3 element -- stationary formulation

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xfluid3_stationary.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"

#include <Epetra_SerialDenseSolver.h>

using namespace XFEM::PHYSICS;

DRT::Elements::XFluid3Stationary::XFluid3Stationary(
		const int iel,
		const int numparamvelx,
		const int numparamvely,
		const int numparamvelz,
		const int numparampres)
  : iel_(iel),
    xyze_(3,iel_,blitz::ColumnMajorArray<2>()),
    edeadng_(3,iel_,blitz::ColumnMajorArray<2>()),
    funct_(iel_),
    deriv_(3,iel_,blitz::ColumnMajorArray<2>()),
    deriv2_(6,iel_,blitz::ColumnMajorArray<2>()),
    xjm_(3,3,blitz::ColumnMajorArray<2>()),
    xji_(3,3,blitz::ColumnMajorArray<2>()),
    vderxy_(3,3,blitz::ColumnMajorArray<2>()),
    pderxy_(3),
    vderxy2_(3,6,blitz::ColumnMajorArray<2>()),
    derxy_(3,iel_,blitz::ColumnMajorArray<2>()),
    derxy2_(6,iel_,blitz::ColumnMajorArray<2>()),
    bodyforce_(3),
    velino_(3),
    velint_(3),
    gradp_(3),
    tau_(3),
    enr_viscs2_(3,3,numparamvelx,blitz::ColumnMajorArray<3>()),
    enr_conv_c_(numparamvelx),
    enr_conv_g_(numparamvelx),
    enr_conv_r_(3,3,numparamvelx,blitz::ColumnMajorArray<3>()),
    rhsint_(3),
    conv_old_(3),
    visc_old_(3),
    res_old_(3),
    xder2_(6,3,blitz::ColumnMajorArray<2>()),
    enr_funct_(numparamvelx),
    enr_derxy_(3,numparamvelx,blitz::ColumnMajorArray<2>()),
    enr_derxy2_(6,numparamvelx,blitz::ColumnMajorArray<2>())
{
}


inline void integrateMatrix(
        const XFEM::ElementDofManager&    dofman,
        blitz::Array<double,2>&           estif,
        const XFEM::PHYSICS::Field&       testfield,
        const blitz::Array<double,1>&     testshape,
        const double&                     fac,
        const XFEM::PHYSICS::Field&       trialfield,
        const blitz::Array<double,1>&     trialshape
        )
{

    const int numparamtest  = dofman.NumDofPerField(testfield);
    const int numparamtrial = dofman.NumDofPerField(trialfield);

    const vector<int> testdof  = dofman.LocalDofPosPerField(testfield);
    const vector<int> trialdof = dofman.LocalDofPosPerField(trialfield);

    for (int ui=0; ui<numparamtrial; ++ui)
    {
        const int trialpos = trialdof[ui];

        for (int vi=0; vi<numparamtest; ++vi)
        {
            const int testpos = testdof[vi];

            estif(testpos, trialpos) += fac*testshape(vi)*trialshape(ui) ;
        }
    }
}

inline void integrateVector(
        const XFEM::ElementDofManager&    dofman,
        blitz::Array<double,1>&           eforce,
        const XFEM::PHYSICS::Field&       testfield,
        const blitz::Array<double,1>&     testshape,
        const double&                     fac
        )
{

    const int numparamtest  = dofman.NumDofPerField(testfield);

    const vector<int> testdof  = dofman.LocalDofPosPerField(testfield);

    for (int vi=0; vi<numparamtest; ++vi)
    {
        const int testpos = testdof[vi];

        eforce(testpos) += fac*testshape(vi);
    }
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (private)                  gjb 11/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Stationary::Sysmat(XFluid3* ele,
                                       const RCP<XFEM::InterfaceHandle>  ih,
                                       const blitz::Array<double,2>&     evelnp,
                                       const blitz::Array<double,1>&     eprenp,
                                       blitz::Array<double,2>&           estif,
                                       blitz::Array<double,1>&           eforce,
                                       struct _MATERIAL*       material,
                                       double                  pseudotime,
                                       bool                    newton ,
                                       bool                    pstab  ,
                                       bool                    supg   ,
                                       bool                    vstab  ,
                                       bool                    cstab  
  )
{

  // set element data
  const DRT::Element::DiscretizationType distype = ele->Shape();

  //! information about domain integration cells
  const XFEM::DomainIntCells   domainIntCells   = ih->domainIntCells(ele->Id(),distype);
  //! information about boundary integration cells
  const XFEM::BoundaryIntCells boundaryIntCells = ih->boundaryIntCells(ele->Id());

  // get node coordinates
  DRT::Node** const nodes = ele->Nodes();
  for (int inode=0; inode<iel_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
  }

  // dead load in element nodes
  ////////////////////////////////////////////////////BodyForce(ele,pseudotime);

  // get viscosity
  // check here, if we really have a fluid !!
  dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
  const double visc = material->m.fluid->viscosity;

  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
//   blitz::fourthIndex l;   // Placeholder for the fourth index

  blitz::Range _  = blitz::Range::all();
//   blitz::Range ux = blitz::Range(0, 4*iel_-4, 4);
//   blitz::Range uy = blitz::Range(1, 4*iel_-3, 4);
//   blitz::Range uz = blitz::Range(2, 4*iel_-2, 4);
//   blitz::Range p  = blitz::Range(3, 4*iel_-1, 4);

  // TODO: put in global dofmanager and check whether the global dofmanager gives the same results as the stored local one. only then we can be sure that we made no mistake
  
  const XFEM::ElementDofManager& dofman = ele->eleDofManager_;

  // stabilization parameter
  // This has to be done before anything else is calculated because
  // we use the same arrays internally.
  CalTauStationary(ele,evelnp,distype,visc);

  // flag for higher order elements
  const bool higher_order_ele = ele->isHigherOrderElement(distype);

  // integrate of integration cell
  dsassert(domainIntCells.empty() == false, "this is a bug!");

  DRT::Utils::GaussRule3D gaussrule;
  bool standard_integration = true;
  if (domainIntCells.size() == 1){
      gaussrule = ele->gaussrule_;
      standard_integration = true;
  }
  else {
      gaussrule = DRT::Utils::intrule_tet_11point;
      standard_integration = false;
  }

  for (XFEM::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell < domainIntCells.end(); ++cell)
  {

  // gaussian points
  const DRT::Utils::IntegrationPoints3D intpoints(gaussrule);

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point in cell coordinates \eta
    const double cell_e0 = intpoints.qxg[iquad][0];
    const double cell_e1 = intpoints.qxg[iquad][1];
    const double cell_e2 = intpoints.qxg[iquad][2];

    const vector<double> e = cell->modifyGaussRule3D(standard_integration,cell_e0,cell_e1,cell_e2);

    // coordinates of the current integration point in element coordinates \xi
    const double e1 = e[0];
    const double e2 = e[1];
    const double e3 = e[2];
    const double detcell = e[3];

    // position of the gausspoint in physical coordinates
    const blitz::Array<double,1> gauss_pos(blitz::sum(funct_(j)*xyze_(i,j),j));

    // shape functions and their derivatives
    DRT::Utils::shape_function_3D(funct_,e1,e2,e3,distype);
    DRT::Utils::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

    // get Jacobian matrix and determinant
    // actually compute its transpose....
    /*
      +-            -+ T      +-            -+
      | dx   dx   dx |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dr   dr   dr |
      |              |        |              |
      | dy   dy   dy |        | dx   dy   dz |
      | --   --   -- |   =    | --   --   -- |
      | dr   ds   dt |        | ds   ds   ds |
      |              |        |              |
      | dz   dz   dz |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dt   dt   dt |
      +-            -+        +-            -+
    */
    xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
    const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                       xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                       xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                       xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                       xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                       xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
    const double fac = intpoints.qwgt[iquad]*det*detcell;

    if (det < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %lf", ele->Id(), det);
    }

    // inverse of jacobian
    xji_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det;
    xji_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det;
    xji_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det;
    xji_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det;
    xji_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det;
    xji_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det;
    xji_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det;
    xji_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det;
    xji_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det;

    // compute global derivates
    derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

    // compute second global derivative
    if (higher_order_ele)
    {
      DRT::Utils::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
      gder2(ele);
    }
    else
    {
      derxy2_  = 0.;
    }


    // enrich shape functions and derivatives
    // simplest case: jump or void enrichment
    // -> no chain rule, since enrichment function derivative is zero

    int dofcounter = 0;
    for (int inode=0; inode<iel_; inode++)
    {
      const int gid = nodes[inode]->Id();
      const blitz::Array<double,1> nodalpos(xyze_(_,inode));

      const std::set<XFEM::FieldEnr>  enrfieldset = dofman.FieldEnrSetPerNode(gid);
      for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
          if (enrfield->getField() == XFEM::PHYSICS::Velx)
          {
              const XFEM::Enrichment enr = enrfield->getEnrichment();
              const double enrval = dofman.enrValue(enr,gauss_pos,nodalpos);
              enr_funct_(dofcounter) = funct_(inode) * enrval;
              enr_derxy_(_,dofcounter) = derxy_(_,inode) * enrval;
              // compute second global derivative
              if (higher_order_ele)
              {
                enr_derxy2_(_,dofcounter) = derxy2_(_,inode) * enrval;
              }
              else
              {
                enr_derxy2_(_,dofcounter) = 0.;
              }
              dofcounter += 1;
          }
      }
    }

    // get velocities (n+g,i) at integration point
    velint_ = blitz::sum(enr_funct_(j)*evelnp(i,j),j);

    // get velocity (np,i) derivatives at integration point
    vderxy_ = blitz::sum(enr_derxy_(j,k)*evelnp(i,k),k);

    // calculate 2nd velocity derivatives at integration point
    if (higher_order_ele)
    {
      vderxy2_ = blitz::sum(enr_derxy2_(j,k)*evelnp(i,k),k);
    }
    else
    {
      vderxy2_ = 0.;
    }

    // get pressure gradients
    gradp_ = blitz::sum(enr_derxy_(i,j)*eprenp(j),j);

    double press = blitz::sum(enr_funct_*eprenp);

    // get bodyforce in gausspoint
    bodyforce_ = 0.0;
    //////////////////////////////////////////bodyforce_ = blitz::sum(enr_edeadng_(i,j)*enr_funct_(j),j);

    // perform integration for entire matrix and rhs

    // stabilisation parameter
    const double tau_M  = tau_(0)*fac;
    const double tau_Mp = tau_(1)*fac;
    const double tau_C  = tau_(2)*fac;

    /*------------------------- evaluate rhs vector at integration point ---*/
    //   rhsint_ = histvec_(i) + bodyforce_(i);
    // histvec is always zero in stationary case (!):
    rhsint_ = bodyforce_(i);    

    /*----------------- get numerical representation of single operators ---*/

    /* Convective term  u_old * grad u_old: */
    conv_old_ = blitz::sum(vderxy_(i, j)*velint_(j), j);

    /* Viscous term  div epsilon(u_old) */
    visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,3) + vderxy2_(0,2) + vderxy2_(2,4));
    visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,3) + vderxy2_(1,2) + vderxy2_(2,5));
    visc_old_(2) = vderxy2_(2,2) + 0.5 * (vderxy2_(2,0) + vderxy2_(0,4) + vderxy2_(2,1) + vderxy2_(1,5));

    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    enr_conv_c_ = blitz::sum(enr_derxy_(j,i)*velint_(j), j);

    /*--- convective grid part u_G * grad (funct) -----------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    enr_conv_g_ = 0.0;


    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /* /                                     \
       |  u_old_x,x   u_old_x,y   u_old x,z  |
       |                                     |
       |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
       |                                     |
       |  u_old_z,x   u_old_z,y   u_old_z,z  |
       \                                     /
       with  N .. form function matrix                                   */
    enr_conv_r_ = vderxy_(i, j)*enr_funct_(k);

    /*--- viscous term  - grad * epsilon(u): ----------------------------*/
    /*   /                                                \
         |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
       1 |                                                |
     - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
       2 |                                                |
         |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
         \                                                /

         with N_x .. x-line of N
         N_y .. y-line of N                                             */

    enr_viscs2_(0,0,_) = - 0.5 * (2.0 * enr_derxy2_(0,_) + enr_derxy2_(1,_) + enr_derxy2_(2,_));
    enr_viscs2_(0,1,_) = - 0.5 *  enr_derxy2_(3,_);
    enr_viscs2_(0,2,_) = - 0.5 *  enr_derxy2_(4,_);
    enr_viscs2_(1,0,_) = - 0.5 *  enr_derxy2_(3,_);
    enr_viscs2_(1,1,_) = - 0.5 * (enr_derxy2_(0,_) + 2.0 * enr_derxy2_(1,_) + enr_derxy2_(2,_));
    enr_viscs2_(1,2,_) = - 0.5 *  enr_derxy2_(5,_);
    enr_viscs2_(2,0,_) = - 0.5 *  enr_derxy2_(4,_);
    enr_viscs2_(2,1,_) = - 0.5 *  enr_derxy2_(5,_);
    enr_viscs2_(2,2,_) = - 0.5 * (enr_derxy2_(0,_) + enr_derxy2_(1,_) + 2.0 * enr_derxy2_(2,_));

    /* pressure gradient term derxy, funct without or with integration   *
     * by parts, respectively                                            */

    /*--------------------------------- now build single stiffness terms ---*/
    {
      // evaluate residual once for all stabilisation right hand sides
      res_old_ = -rhsint_+(conv_old_+gradp_-2*visc*visc_old_);

      //----------------------------------------------------------------------
      //                            GALERKIN PART

      /* convection, convective part */
      /*
                   /                       \
                  |  / n+1       \          |
                  | | u   o nabla | Du , v  |
                  |  \ (i)       /          |
                   \                       /
      */
      integrateMatrix(dofman, estif, Velx, enr_funct_(_), fac, Velx, enr_conv_c_(_));
      integrateMatrix(dofman, estif, Vely, enr_funct_(_), fac, Vely, enr_conv_c_(_));
      integrateMatrix(dofman, estif, Velz, enr_funct_(_), fac, Velz, enr_conv_c_(_));

      /* Viskositaetsterm */
      /*
                    /                        \
                   |       /  \         / \   |
                   |  eps | Du | , eps | v |  |
                   |       \  /         \ /   |
                    \                        /
      */
      integrateMatrix(dofman, estif, Velx, enr_derxy_(0,_), 2.0*visc*fac, Velx, enr_derxy_(0,_));
      integrateMatrix(dofman, estif, Velx, enr_derxy_(1,_),     visc*fac, Velx, enr_derxy_(1,_));
      integrateMatrix(dofman, estif, Velx, enr_derxy_(2,_),     visc*fac, Velx, enr_derxy_(2,_));

      integrateMatrix(dofman, estif, Velx, enr_derxy_(1,_),     visc*fac, Vely, enr_derxy_(0,_));
      integrateMatrix(dofman, estif, Velx, enr_derxy_(2,_),     visc*fac, Velz, enr_derxy_(0,_));
      integrateMatrix(dofman, estif, Vely, enr_derxy_(0,_),     visc*fac, Velx, enr_derxy_(1,_));

      integrateMatrix(dofman, estif, Vely, enr_derxy_(0,_),     visc*fac, Vely, enr_derxy_(0,_));
      integrateMatrix(dofman, estif, Vely, enr_derxy_(1,_), 2.0*visc*fac, Vely, enr_derxy_(1,_));
      integrateMatrix(dofman, estif, Vely, enr_derxy_(2,_),     visc*fac, Vely, enr_derxy_(2,_));

      integrateMatrix(dofman, estif, Vely, enr_derxy_(2,_),     visc*fac, Velz, enr_derxy_(1,_));
      integrateMatrix(dofman, estif, Velz, enr_derxy_(0,_),     visc*fac, Velx, enr_derxy_(2,_));
      integrateMatrix(dofman, estif, Velz, enr_derxy_(1,_),     visc*fac, Vely, enr_derxy_(2,_));

      integrateMatrix(dofman, estif, Velz, enr_derxy_(0,_),     visc*fac, Velz, enr_derxy_(0,_));
      integrateMatrix(dofman, estif, Velz, enr_derxy_(1,_),     visc*fac, Velz, enr_derxy_(1,_));
      integrateMatrix(dofman, estif, Velz, enr_derxy_(2,_), 2.0*visc*fac, Velz, enr_derxy_(2,_));

      /* Druckterm */
      /*
                      /                \
                     |                  |
                     |  Dp , nabla o v  |
                     |                  |
                      \                /
      */
      integrateMatrix(dofman, estif, Velx, enr_derxy_(0,_), -fac, Pres, enr_funct_(_));
      integrateMatrix(dofman, estif, Vely, enr_derxy_(1,_), -fac, Pres, enr_funct_(_));
      integrateMatrix(dofman, estif, Velz, enr_derxy_(2,_), -fac, Pres, enr_funct_(_));

      /* Divergenzfreiheit */
      /*
                     /                \
                    |                  |
                    | nabla o Du  , q  |
                    |                  |
                     \                /
      */
      integrateMatrix(dofman, estif, Pres, enr_funct_(_), fac, Velx, enr_derxy_(0,_));
      integrateMatrix(dofman, estif, Pres, enr_funct_(_), fac, Vely, enr_derxy_(1,_));
      integrateMatrix(dofman, estif, Pres, enr_funct_(_), fac, Velz, enr_derxy_(2,_));


      if (newton)
      {
          /*  convection, reactive part
                 /                         \
                |  /          \   n+1       |
                | | Du o nabla | u     , v  |
                |  \          /   (i)       |
                 \                         /
          */
          integrateMatrix(dofman, estif, Velx, enr_funct_(_), fac, Velx, enr_conv_r_(0, 0, _));
          integrateMatrix(dofman, estif, Velx, enr_funct_(_), fac, Vely, enr_conv_r_(0, 1, _));
          integrateMatrix(dofman, estif, Velx, enr_funct_(_), fac, Velz, enr_conv_r_(0, 2, _));

          integrateMatrix(dofman, estif, Vely, enr_funct_(_), fac, Velx, enr_conv_r_(1, 0, _));
          integrateMatrix(dofman, estif, Vely, enr_funct_(_), fac, Vely, enr_conv_r_(1, 1, _));
          integrateMatrix(dofman, estif, Vely, enr_funct_(_), fac, Velz, enr_conv_r_(1, 2, _));

          integrateMatrix(dofman, estif, Velz, enr_funct_(_), fac, Velx, enr_conv_r_(2, 0, _));
          integrateMatrix(dofman, estif, Velz, enr_funct_(_), fac, Vely, enr_conv_r_(2, 1, _));
          integrateMatrix(dofman, estif, Velz, enr_funct_(_), fac, Velz, enr_conv_r_(2, 2, _));
      }

      /* convection */
      integrateVector(dofman, eforce, Velx, enr_conv_r_(0, 0, _), -fac*velint_(0));
      integrateVector(dofman, eforce, Velx, enr_conv_r_(0, 1, _), -fac*velint_(1));
      integrateVector(dofman, eforce, Velx, enr_conv_r_(0, 2, _), -fac*velint_(2));

      integrateVector(dofman, eforce, Vely, enr_conv_r_(1, 0, _), -fac*velint_(0));
      integrateVector(dofman, eforce, Vely, enr_conv_r_(1, 1, _), -fac*velint_(1));
      integrateVector(dofman, eforce, Vely, enr_conv_r_(1, 2, _), -fac*velint_(2));

      integrateVector(dofman, eforce, Velz, enr_conv_r_(2, 0, _), -fac*velint_(0));
      integrateVector(dofman, eforce, Velz, enr_conv_r_(2, 1, _), -fac*velint_(1));
      integrateVector(dofman, eforce, Velz, enr_conv_r_(2, 2, _), -fac*velint_(2));

      /* pressure */
      integrateVector(dofman, eforce, Velx, enr_derxy_(0, _), press*fac);
      integrateVector(dofman, eforce, Vely, enr_derxy_(1, _), press*fac);
      integrateVector(dofman, eforce, Velz, enr_derxy_(2, _), press*fac);

      /* viscosity */
      integrateVector(dofman, eforce, Velx, enr_derxy_(0, _), -2.0*visc*fac*vderxy_(0, 0));
      integrateVector(dofman, eforce, Velx, enr_derxy_(1, _),     -visc*fac*vderxy_(0, 1));
      integrateVector(dofman, eforce, Velx, enr_derxy_(1, _),     -visc*fac*vderxy_(1, 0));
      integrateVector(dofman, eforce, Velx, enr_derxy_(2, _),     -visc*fac*vderxy_(0, 2));
      integrateVector(dofman, eforce, Velx, enr_derxy_(2, _),     -visc*fac*vderxy_(2, 0));

      integrateVector(dofman, eforce, Vely, enr_derxy_(0, _),     -visc*fac*vderxy_(0, 1));
      integrateVector(dofman, eforce, Vely, enr_derxy_(0, _),     -visc*fac*vderxy_(1, 0));
      integrateVector(dofman, eforce, Vely, enr_derxy_(1, _), -2.0*visc*fac*vderxy_(1, 1));
      integrateVector(dofman, eforce, Vely, enr_derxy_(2, _),     -visc*fac*vderxy_(1, 2));
      integrateVector(dofman, eforce, Vely, enr_derxy_(2, _),     -visc*fac*vderxy_(2, 1));

      integrateVector(dofman, eforce, Velz, enr_derxy_(0, _),     -visc*fac*vderxy_(0, 2));
      integrateVector(dofman, eforce, Velz, enr_derxy_(0, _),     -visc*fac*vderxy_(2, 0));
      integrateVector(dofman, eforce, Velz, enr_derxy_(1, _),     -visc*fac*vderxy_(1, 2));
      integrateVector(dofman, eforce, Velz, enr_derxy_(1, _),     -visc*fac*vderxy_(2, 1));
      integrateVector(dofman, eforce, Velz, enr_derxy_(2, _), -2.0*visc*fac*vderxy_(2, 2));

      // source term of the right hand side
      integrateVector(dofman, eforce, Velx, enr_funct_(_), fac*rhsint_(0));
      integrateVector(dofman, eforce, Vely, enr_funct_(_), fac*rhsint_(1));
      integrateVector(dofman, eforce, Velz, enr_funct_(_), fac*rhsint_(2));

      // continuity equation
      integrateVector(dofman, eforce, Pres, enr_conv_r_(0, 0, _), -fac);
      integrateVector(dofman, eforce, Pres, enr_conv_r_(1, 1, _), -fac);
      integrateVector(dofman, eforce, Pres, enr_conv_r_(2, 2, _), -fac);

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if(pstab)
      {


        /* pressure stabilisation: convection, convective part */
        /*
                  /                            \
                 |  / n+1       \               |
                 | | u   o nabla | Du , nabla q |
                 |  \ (i)       /               |
                  \                            /
        */
        integrateMatrix(dofman, estif, Pres, enr_derxy_(0,_), tau_Mp, Velx, enr_conv_c_(_));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(1,_), tau_Mp, Vely, enr_conv_c_(_));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(2,_), tau_Mp, Velz, enr_conv_c_(_));

        /* pressure stabilisation: viscosity (-L_visc_u) */
        /*
                 /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                 \                              /
        */
        integrateMatrix(dofman, estif, Pres, enr_derxy_(0,_), 2.0*visc*tau_Mp, Velx, enr_viscs2_(0, 0, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(1,_), 2.0*visc*tau_Mp, Velx, enr_viscs2_(0, 1, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(2,_), 2.0*visc*tau_Mp, Velx, enr_viscs2_(0, 2, _));

        integrateMatrix(dofman, estif, Pres, enr_derxy_(0,_), 2.0*visc*tau_Mp, Vely, enr_viscs2_(0, 1, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(1,_), 2.0*visc*tau_Mp, Vely, enr_viscs2_(1, 1, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(2,_), 2.0*visc*tau_Mp, Vely, enr_viscs2_(1, 2, _));

        integrateMatrix(dofman, estif, Pres, enr_derxy_(0,_), 2.0*visc*tau_Mp, Velz, enr_viscs2_(0, 2, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(1,_), 2.0*visc*tau_Mp, Velz, enr_viscs2_(1, 2, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(2,_), 2.0*visc*tau_Mp, Velz, enr_viscs2_(2, 2, _));

        /* pressure stabilisation: pressure( L_pres_p) */
        /*
                    /                    \
                   |                      |
                   |  nabla Dp , nabla q  |
                   |                      |
                    \                    /
        */
        integrateMatrix(dofman, estif, Pres, enr_derxy_(0,_), tau_Mp, Pres, enr_derxy_(0, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(1,_), tau_Mp, Pres, enr_derxy_(1, _));
        integrateMatrix(dofman, estif, Pres, enr_derxy_(2,_), tau_Mp, Pres, enr_derxy_(2, _));

        if (newton)
        {
          dserror("not translated to XFEM integrator yet!");
          for (int ui=0; ui<iel_; ++ui)
          {
            const int UDF = ui*4;
            const int VDF = ui*4+1;
            const int WDF = ui*4+2;
            for (int vi=0; vi<iel_; ++vi)
            {
              const int dPDF = vi*4+3;
              /*  pressure stabilisation: convection, reactive part

                  /                             \
                 |  /          \   n+1           |
                 | | Du o nabla | u     , grad q |
                 |  \          /   (i)           |
                  \                             /

              */
              estif(dPDF, UDF) += tau_Mp*(enr_derxy_(0, vi)*enr_conv_r_(0, 0, ui)
                                          +
                                          enr_derxy_(1, vi)*enr_conv_r_(1, 0, ui)
                                          +
                                          enr_derxy_(2, vi)*enr_conv_r_(2, 0, ui)) ;
              estif(dPDF, VDF) += tau_Mp*(enr_derxy_(0, vi)*enr_conv_r_(0, 1, ui)
                                          +
                                          enr_derxy_(1, vi)*enr_conv_r_(1, 1, ui)
                                          +
                                          enr_derxy_(2, vi)*enr_conv_r_(2, 1, ui)) ;
              estif(dPDF, WDF) += tau_Mp*(enr_derxy_(0, vi)*enr_conv_r_(0, 2, ui)
                                          +
                                          enr_derxy_(1, vi)*enr_conv_r_(1, 2, ui)
                                          +
                                          enr_derxy_(2, vi)*enr_conv_r_(2, 2, ui)) ;

            } // vi
          } // ui
        } // if newton

        // pressure stabilisation
        integrateVector(dofman, eforce, Pres, enr_derxy_(0, _), -tau_Mp*res_old_(0));
        integrateVector(dofman, eforce, Pres, enr_derxy_(1, _), -tau_Mp*res_old_(1));
        integrateVector(dofman, eforce, Pres, enr_derxy_(2, _), -tau_Mp*res_old_(2));
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg)
      {
          /* supg stabilisation: convective part ( L_conv_u) */
          /*
               /                                           \
              |    / n+1        \        / n+1        \     |
              |   | u    o nabla | Du , | u    o nabla | v  |
              |    \ (i)        /        \ (i)        /     |
               \                                           /
          */
          integrateMatrix(dofman, estif, Velx, enr_conv_c_(_), tau_M, Velx, enr_conv_c_(_));
          integrateMatrix(dofman, estif, Vely, enr_conv_c_(_), tau_M, Vely, enr_conv_c_(_));
          integrateMatrix(dofman, estif, Velz, enr_conv_c_(_), tau_M, Velz, enr_conv_c_(_));
          
          /* supg stabilisation: pressure part  ( L_pres_p) */
          /*
                    /                              \
                   |              / n+1       \     |
                   |  nabla Dp , | u   o nabla | v  |
                   |              \ (i)       /     |
                    \                              /
          */
          integrateMatrix(dofman, estif, Velx, enr_conv_c_(_), tau_M, Pres, enr_derxy_(0, _));
          integrateMatrix(dofman, estif, Vely, enr_conv_c_(_), tau_M, Pres, enr_derxy_(1, _));
          integrateMatrix(dofman, estif, Velz, enr_conv_c_(_), tau_M, Pres, enr_derxy_(2, _));
          
          /* supg stabilisation: viscous part  (-L_visc_u) */
          /*
                /                                        \
               |               /  \    / n+1        \     |
               |  nabla o eps | Du |, | u    o nabla | v  |
               |               \  /    \ (i)        /     |
                \                                        /
          */
          integrateMatrix(dofman, estif, Velx, enr_conv_c_(_), 2.0*visc*tau_M, Velx, enr_viscs2_(0, 0, _));
          integrateMatrix(dofman, estif, Velx, enr_conv_c_(_), 2.0*visc*tau_M, Vely, enr_viscs2_(0, 1, _));
          integrateMatrix(dofman, estif, Velx, enr_conv_c_(_), 2.0*visc*tau_M, Velz, enr_viscs2_(0, 2, _));
          
          integrateMatrix(dofman, estif, Vely, enr_conv_c_(_), 2.0*visc*tau_M, Velx, enr_viscs2_(0, 1, _));
          integrateMatrix(dofman, estif, Vely, enr_conv_c_(_), 2.0*visc*tau_M, Vely, enr_viscs2_(1, 1, _));
          integrateMatrix(dofman, estif, Vely, enr_conv_c_(_), 2.0*visc*tau_M, Velz, enr_viscs2_(1, 2, _));
          
          integrateMatrix(dofman, estif, Velz, enr_conv_c_(_), 2.0*visc*tau_M, Velx, enr_viscs2_(0, 2, _));
          integrateMatrix(dofman, estif, Velz, enr_conv_c_(_), 2.0*visc*tau_M, Vely, enr_viscs2_(1, 2, _));
          integrateMatrix(dofman, estif, Velz, enr_conv_c_(_), 2.0*visc*tau_M, Velz, enr_viscs2_(2, 2, _));

        if (newton)
        {
          for (int ui=0; ui<iel_; ++ui)
          {
            const int UDF = ui*4;
            const int VDF = ui*4+1;
            const int WDF = ui*4+2;
            for (int vi=0; vi<iel_; ++vi)
            {
              const int dUDF = vi*4;
              const int dVDF = vi*4+1;
              const int dWDF = vi*4+2;
              /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
              /*
                       /                                           \
                      |    / n+1        \   n+1    /          \     |
                      |   | u    o nabla | u    , | Du o nabla | v  |
                      |    \ (i)        /   (i)    \          /     |
                       \                                           /

                       /                                           \
                      |    /          \   n+1    / n+1        \     |
                      |   | Du o nabla | u    , | u    o nabla | v  |
                      |    \          /   (i)    \ (i)        /     |
                       \                                           /
              */
              estif(dUDF, UDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(0, 0, ui)
                                         +
                                         velint_(0)*enr_derxy_(0, vi)*enr_conv_r_(0, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(0, vi)*enr_conv_r_(0, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(0, vi)*enr_conv_r_(0, 2, ui)) ;
              estif(dUDF, VDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(0, 1, ui)
                                         +
                                         velint_(0)*enr_derxy_(1, vi)*enr_conv_r_(0, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(1, vi)*enr_conv_r_(0, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(1, vi)*enr_conv_r_(0, 2, ui)) ;
              estif(dUDF, WDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(0, 2, ui)
                                         +
                                         velint_(0)*enr_derxy_(2, vi)*enr_conv_r_(0, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(2, vi)*enr_conv_r_(0, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(2, vi)*enr_conv_r_(0, 2, ui)) ;
              estif(dVDF, UDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(1, 0, ui)
                                         +
                                         velint_(0)*enr_derxy_(0, vi)*enr_conv_r_(1, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(0, vi)*enr_conv_r_(1, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(0, vi)*enr_conv_r_(1, 2, ui)) ;
              estif(dVDF, VDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(1, 1, ui)
                                         +
                                         velint_(0)*enr_derxy_(1, vi)*enr_conv_r_(1, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(1, vi)*enr_conv_r_(1, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(1, vi)*enr_conv_r_(1, 2, ui)) ;
              estif(dVDF, WDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(1, 2, ui)
                                         +
                                         velint_(0)*enr_derxy_(2, vi)*enr_conv_r_(1, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(2, vi)*enr_conv_r_(1, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(2, vi)*enr_conv_r_(1, 2, ui)) ;
              estif(dWDF, UDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(2, 0, ui)
                                         +
                                         velint_(0)*enr_derxy_(0, vi)*enr_conv_r_(2, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(0, vi)*enr_conv_r_(2, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(0, vi)*enr_conv_r_(2, 2, ui)) ;
              estif(dWDF, VDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(2, 1, ui)
                                         +
                                         velint_(0)*enr_derxy_(1, vi)*enr_conv_r_(2, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(1, vi)*enr_conv_r_(2, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(1, vi)*enr_conv_r_(2, 2, ui)) ;
              estif(dWDF, WDF) += tau_M*(enr_conv_c_(vi)*enr_conv_r_(2, 2, ui)
                                         +
                                         velint_(0)*enr_derxy_(2, vi)*enr_conv_r_(2, 0, ui)
                                         +
                                         velint_(1)*enr_derxy_(2, vi)*enr_conv_r_(2, 1, ui)
                                         +
                                         velint_(2)*enr_derxy_(2, vi)*enr_conv_r_(2, 2, ui)) ;


              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                            /                               \
                           |         n+1    /          \     |
                           |  nabla p    , | Du o nabla | v  |
                           |         (i)    \          /     |
                            \                               /
              */
              estif(dUDF, UDF) += tau_M*enr_funct_(ui)*gradp_(0)*enr_derxy_(0, vi) ;
              estif(dUDF, VDF) += tau_M*enr_funct_(ui)*gradp_(0)*enr_derxy_(1, vi) ;
              estif(dUDF, WDF) += tau_M*enr_funct_(ui)*gradp_(0)*enr_derxy_(2, vi) ;
              estif(dVDF, UDF) += tau_M*enr_funct_(ui)*gradp_(1)*enr_derxy_(0, vi) ;
              estif(dVDF, VDF) += tau_M*enr_funct_(ui)*gradp_(1)*enr_derxy_(1, vi) ;
              estif(dVDF, WDF) += tau_M*enr_funct_(ui)*gradp_(1)*enr_derxy_(2, vi) ;
              estif(dWDF, UDF) += tau_M*enr_funct_(ui)*gradp_(2)*enr_derxy_(0, vi) ;
              estif(dWDF, VDF) += tau_M*enr_funct_(ui)*gradp_(2)*enr_derxy_(1, vi) ;
              estif(dWDF, WDF) += tau_M*enr_funct_(ui)*gradp_(2)*enr_derxy_(2, vi) ;


              /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
              /*
                      /                                         \
                     |               / n+1 \    /          \     |
                     |  nabla o eps | u     |, | Du o nabla | v  |
                     |               \ (i) /    \          /     |
                      \                                         /
              */
              estif(dUDF, UDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(0)*enr_derxy_(0, vi) ;
              estif(dUDF, VDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(0)*enr_derxy_(1, vi) ;
              estif(dUDF, WDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(0)*enr_derxy_(2, vi) ;
              estif(dVDF, UDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(1)*enr_derxy_(0, vi) ;
              estif(dVDF, VDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(1)*enr_derxy_(1, vi) ;
              estif(dVDF, WDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(1)*enr_derxy_(2, vi) ;
              estif(dWDF, UDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(2)*enr_derxy_(0, vi) ;
              estif(dWDF, VDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(2)*enr_derxy_(1, vi) ;
              estif(dWDF, WDF) += -2.0*visc*tau_M*enr_funct_(ui)*visc_old_(2)*enr_derxy_(2, vi) ;


              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                          /                             \
                         |              /          \     |
                         |  rhsint   , | Du o nabla | v  |
                         |              \          /     |
                          \                             /

              */
              estif(dUDF, UDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(0, vi)*rhsint_(0)) ;
              estif(dUDF, VDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(1, vi)*rhsint_(0)) ;
              estif(dUDF, WDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(2, vi)*rhsint_(0)) ;
              estif(dVDF, UDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(0, vi)*rhsint_(1)) ;
              estif(dVDF, VDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(1, vi)*rhsint_(1)) ;
              estif(dVDF, WDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(2, vi)*rhsint_(1)) ;
              estif(dWDF, UDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(0, vi)*rhsint_(2)) ;
              estif(dWDF, VDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(1, vi)*rhsint_(2)) ;
              estif(dWDF, WDF) += -(tau_M*enr_funct_(ui)*enr_derxy_(2, vi)*rhsint_(2)) ;

            } // vi
          } // ui
        } // if newton

        // supg stabilisation
        integrateVector(dofman, eforce, Velx, enr_conv_c_(_), -tau_M*res_old_(0));
        integrateVector(dofman, eforce, Vely, enr_conv_c_(_), -tau_M*res_old_(1));
        integrateVector(dofman, eforce, Velz, enr_conv_c_(_), -tau_M*res_old_(2));
      }


      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART


//  !!!!
//  viscous part of stabilisation is switched off!
//  vstab is set to false within fluid3_evaluate.cpp
//!!!!

      if(vstab)
      {
        for (int ui=0; ui<iel_; ++ui)
        {
          const int UDF = ui*4;
          const int VDF = ui*4+1;
          const int WDF = ui*4+2;
          const int PDF = ui*4+3;
          for (int vi=0; vi<iel_; ++vi)
          {
            const int dUDF = vi*4;
            const int dVDF = vi*4+1;
            const int dWDF = vi*4+2;
            /* viscous stabilisation, inertia part */
            /*
                        /                  \
                       |                    |
                       |  Du , div eps (v)  |
                       |                    |
                        \                  /
            */
            estif(dUDF, UDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(0, 0, vi) ;
            estif(dUDF, VDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(0, 1, vi) ;
            estif(dUDF, WDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(0, 2, vi) ;
            estif(dVDF, UDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(0, 1, vi) ;
            estif(dVDF, VDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(1, 1, vi) ;
            estif(dVDF, WDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(1, 2, vi) ;
            estif(dWDF, UDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(0, 2, vi) ;
            estif(dWDF, VDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(1, 2, vi) ;
            estif(dWDF, WDF) += 2.0*visc*tau_Mp*enr_funct_(ui)*enr_viscs2_(2, 2, vi) ;

            /* viscous stabilisation, convective part */
            /*
                 /                                \
                |  / n+1       \                   |
                | | u   o nabla | Du , div eps (v) |
                |  \ (i)       /                   |
                 \                                /
            */
            estif(dUDF, UDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(0, 0, vi) ;
            estif(dUDF, VDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(0, 1, vi) ;
            estif(dUDF, WDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(0, 2, vi) ;
            estif(dVDF, UDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(0, 1, vi) ;
            estif(dVDF, VDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(1, 1, vi) ;
            estif(dVDF, WDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(1, 2, vi) ;
            estif(dWDF, UDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(0, 2, vi) ;
            estif(dWDF, VDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(1, 2, vi) ;
            estif(dWDF, WDF) += 2.0*visc*tau_Mp*enr_conv_c_(ui)*enr_viscs2_(2, 2, vi) ;


            /* viscous stabilisation, pressure part ( L_pres_p) */
            /*
                     /                        \
                    |                          |
                    |  nabla Dp , div eps (v)  |
                    |                          |
                     \                        /
            */
            estif(dUDF, PDF) += 2.0*visc*tau_Mp*(enr_derxy_(0, ui)*enr_viscs2_(0, 0, vi)
                                                 +
                                                 enr_derxy_(1, ui)*enr_viscs2_(0, 1, vi)
                                                 +
                                                 enr_derxy_(2, ui)*enr_viscs2_(0, 2, vi)) ;
            estif(dVDF, PDF) += 2.0*visc*tau_Mp*(enr_derxy_(0, ui)*enr_viscs2_(0, 1, vi)
                                                 +
                                                 enr_derxy_(1, ui)*enr_viscs2_(1, 1, vi)
                                                 +
                                                 enr_derxy_(2, ui)*enr_viscs2_(1, 2, vi)) ;
            estif(dWDF, PDF) += 2.0*visc*tau_Mp*(enr_derxy_(0, ui)*enr_viscs2_(0, 2, vi)
                                                 +
                                                 enr_derxy_(1, ui)*enr_viscs2_(1, 2, vi)
                                                 +
                                                 enr_derxy_(2, ui)*enr_viscs2_(2, 2, vi)) ;

            /* viscous stabilisation, viscous part (-L_visc_u) */
            /*
               /                                 \
              |               /  \                |
              |  nabla o eps | Du | , div eps (v) |
              |               \  /                |
               \                                 /
            */
            estif(dUDF, UDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 0, ui)*enr_viscs2_(0, 0, vi)
                                                        +
                                                        enr_viscs2_(0, 1, ui)*enr_viscs2_(0, 1, vi)
                                                        +
                                                        enr_viscs2_(0, 2, ui)*enr_viscs2_(0, 2, vi)) ;
            estif(dUDF, VDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 0, vi)*enr_viscs2_(0, 1, ui)
                                                        +
                                                        enr_viscs2_(0, 1, vi)*enr_viscs2_(1, 1, ui)
                                                        +
                                                        enr_viscs2_(0, 2, vi)*enr_viscs2_(1, 2, ui)) ;
            estif(dUDF, WDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 0, vi)*enr_viscs2_(0, 2, ui)
                                                        +
                                                        enr_viscs2_(0, 1, vi)*enr_viscs2_(1, 2, ui)
                                                        +
                                                        enr_viscs2_(0, 2, vi)*enr_viscs2_(2, 2, ui)) ;
            estif(dVDF, UDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 0, ui)*enr_viscs2_(0, 1, vi)
                                                        +
                                                        enr_viscs2_(0, 1, ui)*enr_viscs2_(1, 1, vi)
                                                        +
                                                        enr_viscs2_(0, 2, ui)*enr_viscs2_(1, 2, vi)) ;
            estif(dVDF, VDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 1, ui)*enr_viscs2_(0, 1, vi)
                                                        +
                                                        enr_viscs2_(1, 1, ui)*enr_viscs2_(1, 1, vi)
                                                        +
                                                        enr_viscs2_(1, 2, ui)*enr_viscs2_(1, 2, vi)) ;
            estif(dVDF, WDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 1, vi)*enr_viscs2_(0, 2, ui)
                                                        +
                                                        enr_viscs2_(1, 1, vi)*enr_viscs2_(1, 2, ui)
                                                        +
                                                        enr_viscs2_(1, 2, vi)*enr_viscs2_(2, 2, ui)) ;
            estif(dWDF, UDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 0, ui)*enr_viscs2_(0, 2, vi)
                                                        +
                                                        enr_viscs2_(0, 1, ui)*enr_viscs2_(1, 2, vi)
                                                        +
                                                        enr_viscs2_(0, 2, ui)*enr_viscs2_(2, 2, vi)) ;
            estif(dWDF, VDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 1, ui)*enr_viscs2_(0, 2, vi)
                                                        +
                                                        enr_viscs2_(1, 1, ui)*enr_viscs2_(1, 2, vi)
                                                        +
                                                        enr_viscs2_(1, 2, ui)*enr_viscs2_(2, 2, vi)) ;
            estif(dWDF, WDF) += 4.0*(visc*visc)*tau_Mp*(enr_viscs2_(0, 2, ui)*enr_viscs2_(0, 2, vi)
                                                        +
                                                        enr_viscs2_(1, 2, ui)*enr_viscs2_(1, 2, vi)
                                                        +
                                                        enr_viscs2_(2, 2, ui)*enr_viscs2_(2, 2, vi)) ;
          } // vi
        } // ui

        if (newton)
        {
          for (int ui=0; ui<iel_; ++ui)
          {
            const int UDF = ui*4;
            const int VDF = ui*4+1;
            const int WDF = ui*4+2;
            for (int vi=0; vi<iel_; ++vi)
            {
              const int dUDF = vi*4;
              const int dVDF = vi*4+1;
              const int dWDF = vi*4+2;
              /* viscous stabilisation, reactive part of convection */
              /*
                   /                                 \
                  |  /          \   n+1               |
                  | | Du o nabla | u    , div eps (v) |
                  |  \          /   (i)               |
                   \                                 /
              */
              estif(dUDF, UDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 0, vi)*enr_conv_r_(0, 0, ui)
                                                   +
                                                   enr_viscs2_(0, 1, vi)*enr_conv_r_(1, 0, ui)
                                                   +
                                                   enr_viscs2_(0, 2, vi)*enr_conv_r_(2, 0, ui)) ;
              estif(dUDF, VDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 0, vi)*enr_conv_r_(0, 1, ui)
                                                   +
                                                   enr_viscs2_(0, 1, vi)*enr_conv_r_(1, 1, ui)
                                                   +
                                                   enr_viscs2_(0, 2, vi)*enr_conv_r_(2, 1, ui)) ;
              estif(dUDF, WDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 0, vi)*enr_conv_r_(0, 2, ui)
                                                   +
                                                   enr_viscs2_(0, 1, vi)*enr_conv_r_(1, 2, ui)
                                                   +
                                                   enr_viscs2_(0, 2, vi)*enr_conv_r_(2, 2, ui)) ;
              estif(dVDF, UDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 1, vi)*enr_conv_r_(0, 0, ui)
                                                   +
                                                   enr_viscs2_(1, 1, vi)*enr_conv_r_(1, 0, ui)
                                                   +
                                                   enr_viscs2_(1, 2, vi)*enr_conv_r_(2, 0, ui)) ;
              estif(dVDF, VDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 1, vi)*enr_conv_r_(0, 1, ui)
                                                   +
                                                   enr_viscs2_(1, 1, vi)*enr_conv_r_(1, 1, ui)
                                                   +
                                                   enr_viscs2_(1, 2, vi)*enr_conv_r_(2, 1, ui)) ;
              estif(dVDF, WDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 1, vi)*enr_conv_r_(0, 2, ui)
                                                   +
                                                   enr_viscs2_(1, 1, vi)*enr_conv_r_(1, 2, ui)
                                                   +
                                                   enr_viscs2_(1, 2, vi)*enr_conv_r_(2, 2, ui)) ;
              estif(dWDF, UDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 2, vi)*enr_conv_r_(0, 0, ui)
                                                   +
                                                   enr_viscs2_(1, 2, vi)*enr_conv_r_(1, 0, ui)
                                                   +
                                                   enr_viscs2_(2, 2, vi)*enr_conv_r_(2, 0, ui)) ;
              estif(dWDF, VDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 2, vi)*enr_conv_r_(0, 1, ui)
                                                   +
                                                   enr_viscs2_(1, 2, vi)*enr_conv_r_(1, 1, ui)
                                                   +
                                                   enr_viscs2_(2, 2, vi)*enr_conv_r_(2, 1, ui)) ;
              estif(dWDF, WDF) += 2.0*visc*tau_Mp*(enr_viscs2_(0, 2, vi)*enr_conv_r_(0, 2, ui)
                                                   +
                                                   enr_viscs2_(1, 2, vi)*enr_conv_r_(1, 2, ui)
                                                   +
                                                   enr_viscs2_(2, 2, vi)*enr_conv_r_(2, 2, ui)) ;
            } // vi
          } // ui
        } // if newton

        /* viscous stabilisation */
        integrateVector(dofman, eforce, Velx, enr_viscs2_(0, 0, _), -2.0*visc*tau_Mp*res_old_(0));
        integrateVector(dofman, eforce, Velx, enr_viscs2_(0, 1, _), -2.0*visc*tau_Mp*res_old_(1));
        integrateVector(dofman, eforce, Velx, enr_viscs2_(0, 2, _), -2.0*visc*tau_Mp*res_old_(2));
        
        integrateVector(dofman, eforce, Vely, enr_viscs2_(0, 1, _), -2.0*visc*tau_Mp*res_old_(0));
        integrateVector(dofman, eforce, Vely, enr_viscs2_(1, 1, _), -2.0*visc*tau_Mp*res_old_(1));
        integrateVector(dofman, eforce, Vely, enr_viscs2_(1, 2, _), -2.0*visc*tau_Mp*res_old_(2));
                
        integrateVector(dofman, eforce, Velz, enr_viscs2_(0, 2, _), -2.0*visc*tau_Mp*res_old_(0));
        integrateVector(dofman, eforce, Velz, enr_viscs2_(1, 2, _), -2.0*visc*tau_Mp*res_old_(1));
        integrateVector(dofman, eforce, Velz, enr_viscs2_(2, 2, _), -2.0*visc*tau_Mp*res_old_(2));
        
      } // endif vstab

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if(cstab)
      {
        /* continuity stabilisation */
        /*
                   /                        \
                  |                          |
                  | nabla o Du  , nabla o v  |
                  |                          |
                   \                        /
        */
        integrateMatrix(dofman, estif, Velx, enr_derxy_(0, _), tau_C, Velx, enr_derxy_(0, _));
        integrateMatrix(dofman, estif, Velx, enr_derxy_(0, _), tau_C, Vely, enr_derxy_(1, _));
        integrateMatrix(dofman, estif, Velx, enr_derxy_(0, _), tau_C, Velz, enr_derxy_(2, _));
        
        integrateMatrix(dofman, estif, Vely, enr_derxy_(1, _), tau_C, Velx, enr_derxy_(0, _));
        integrateMatrix(dofman, estif, Vely, enr_derxy_(1, _), tau_C, Vely, enr_derxy_(1, _));
        integrateMatrix(dofman, estif, Vely, enr_derxy_(1, _), tau_C, Velz, enr_derxy_(2, _));
        
        integrateMatrix(dofman, estif, Velz, enr_derxy_(2, _), tau_C, Velx, enr_derxy_(0, _));
        integrateMatrix(dofman, estif, Velz, enr_derxy_(2, _), tau_C, Vely, enr_derxy_(1, _));
        integrateMatrix(dofman, estif, Velz, enr_derxy_(2, _), tau_C, Velz, enr_derxy_(2, _));
        
        const double tau_C_divunp=tau_C*(vderxy_(0, 0)+vderxy_(1, 1)+vderxy_(2, 2));

        integrateVector(dofman, eforce, Velx, enr_derxy_(0, _), -tau_C_divunp);
        integrateVector(dofman, eforce, Vely, enr_derxy_(1, _), -tau_C_divunp);
        integrateVector(dofman, eforce, Velz, enr_derxy_(2, _), -tau_C_divunp);

      } // endif cstab
    }
  }
  } // end loop over integration cells
  return;
}



//
// calculate stabilization parameter
//
void DRT::Elements::XFluid3Stationary::CalTauStationary(
  XFluid3* ele,
  const blitz::Array<double,2>&           evelnp,
  const DRT::Element::DiscretizationType  distype,
  const double                            visc
  )
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // use one point gauss rule to calculate tau at element center
  DRT::Utils::GaussRule3D integrationrule_stabili=DRT::Utils::intrule3D_undefined;
  switch (distype)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    integrationrule_stabili = DRT::Utils::intrule_hex_1point;
    break;
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    integrationrule_stabili = DRT::Utils::intrule_tet_1point;
    break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    integrationrule_stabili = DRT::Utils::intrule_wedge_1point;
    break;
  case DRT::Element::pyramid5:
    integrationrule_stabili = DRT::Utils::intrule_pyramid_1point;
    break;
  default:
    dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::Utils::IntegrationPoints3D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];
  const double e3    = intpoints.qxg[0][2];
  const double wquad = intpoints.qwgt[0];

  DRT::Utils::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::Utils::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
  case DRT::Element::tet4:
  case DRT::Element::pyramid5:
  case DRT::Element::hex8:
  case DRT::Element::wedge6:
    mk = 0.333333333333333333333;
    break;
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge15:
    mk = 0.083333333333333333333;
    break;
  default:
    dserror("type unknown!\n");
  }

  // get velocities at element center
  velint_ = blitz::sum(funct_(j)*evelnp(i,j),j);

  // get Jacobian matrix and determinant
  xjm_ = blitz::sum(deriv_(i,k)*xyze_(j,k),k);
  const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                     xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                     xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                     xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                     xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                     xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
  const double vol = wquad*det;

  // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol/PI),(1.0/3.0))/sqrt(3.0);

  // inverse of jacobian
  xji_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det;
  xji_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det;
  xji_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det;
  xji_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det;
  xji_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det;
  xji_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det;
  xji_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det;
  xji_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det;
  xji_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det;

  // compute global derivates
  derxy_ = blitz::sum(xji_(i,k)*deriv_(k,j),k);

  // get velocity norm
  const double vel_norm = sqrt(blitz::sum(velint_*velint_));

  // normed velocity at element centre
  if (vel_norm>=1e-6)
  {
    velino_ = velint_/vel_norm;
  }
  else
  {
    velino_ = 0.;
    velino_(0) = 1;
  }

  // get streamlength
  const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
  const double strle = 2.0/val;

  // calculate tau
  // stabilization parameters for stationary case

  // compute tau_Mu
  const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
  const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
  tau_(0) = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);

  // compute tau_Mp
  const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
  const double xi_tau_mp = DMAX(re_tau_mp,1.0);
  tau_(1) = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);

  // compute tau_C
  const double xi_tau_c = DMIN(re_tau_mp, 1.0);
  tau_(2) = 0.5*vel_norm*hk*xi_tau_c;
  
}



/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Stationary::BodyForce(XFluid3* ele, const double pseudotime)
{
  vector<DRT::Condition*> myneumcond;
  DRT::Node** nodes = ele->Nodes();
  dserror("not adapted to xfem (, yet)!!!");
  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for (int inode=0;inode<iel_;inode++)
  {
    nodes[inode]->GetCondition("VolumeNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one VolumeNeumann cond on one node");
    }
    if (myneumcond.size()==1)
    {
      nodecount++;
    }
  }

  if (nodecount == iel_)
  {
    // find out whether we will use a (pseudo-)time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a (pseudo-)timecurve
    {
      // factor for the intermediate step
      if(pseudotime >= 0.0)
      {
        curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(pseudotime);
      }
      else
      {
	// do not compute an "alternative" curvefac here since a negative pseudotime value
	// indicates an error.
        dserror("Negative pseudotime value in body force calculation: time = %f",pseudotime);
        //curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a (pseudo-)timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for (int jnode=0; jnode<iel_; jnode++)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<3;isd++)
      {
        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_ = 0.;
  }
}


/*----------------------------------------------------------------------*
 |  calculate second global derivatives w.r.t. x,y,z at point r,s,t
 |                                            (private)      gammi 07/07
 |
 | From the six equations
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dr^2     dr | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ------ = -- | --*-- + --*-- + --*-- |
 |  ds^2     ds | ds dx   ds dy   ds dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dt^2     dt | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dr     ds | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | dt dr     dt | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dt     ds | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 | the matrix (jacobian-bar matrix) system
 |
 | +-                                                                                         -+   +-    -+
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dr/          \dr/           \dr/             dr dr           dr dr           dr dr     |   | dx^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \ds/          \ds/           \ds/             ds ds           ds ds           ds ds     |   | dy^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dt/          \dt/           \dt/             dt dt           dt dt           dt dt     |   | dz^2 |
 | |                                                                                           | * |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr ds         dr ds          dr ds        dr ds   ds dr   dr ds   ds dr  dr ds   ds dr  |   | dxdy |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr dt         dr dt          dr dt        dr dt   dt dr   dr dt   dt dr  dr dt   dt dr  |   | dxdz |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dt ds         dt ds          dt ds        dt ds   ds dt   dt ds   ds dt  dt ds   ds dt  |   | dydz |
 | +-                                                                                         -+   +-    -+
 |
 |                  +-    -+     +-                           -+
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dr^2 |     | dr^2 dx   dr^2 dy   dr^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | ds^2 |     | ds^2 dx   ds^2 dy   ds^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dt^2 |     | dt^2 dx   dt^2 dy   dt^2 dz |
 |              =   |      |  -  |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drds |     | drds dx   drds dy   drds dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drdt |     | drdt dx   drdt dy   drdt dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2z dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dtds |     | dtds dx   dtds dy   dtds dz |
 |                  +-    -+     +-                           -+
 |
 |
 | is derived. This is solved for the unknown global derivatives.
 |
 |
 |             jacobian_bar * derxy2 = deriv2 - xder2 * derxy
 |                                              |           |
 |                                              +-----------+
 |                                              'chainrulerhs'
 |                                     |                    |
 |                                     +--------------------+
 |                                          'chainrulerhs'
 |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3Stationary::gder2(XFluid3* ele)
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // initialize and zero out everything
  static Epetra_SerialDenseMatrix bm(6,6);

  // calculate elements of jacobian_bar matrix
  bm(0,0) = xjm_(0,0)*xjm_(0,0);
  bm(1,0) = xjm_(1,0)*xjm_(1,0);
  bm(2,0) = xjm_(2,0)*xjm_(2,0);
  bm(3,0) = xjm_(0,0)*xjm_(1,0);
  bm(4,0) = xjm_(0,0)*xjm_(2,0);
  bm(5,0) = xjm_(2,0)*xjm_(1,0);

  bm(0,1) = xjm_(0,1)*xjm_(0,1);
  bm(1,1) = xjm_(1,1)*xjm_(1,1);
  bm(2,1) = xjm_(2,1)*xjm_(2,1);
  bm(3,1) = xjm_(0,1)*xjm_(1,1);
  bm(4,1) = xjm_(0,1)*xjm_(2,1);
  bm(5,1) = xjm_(2,1)*xjm_(1,1);

  bm(0,2) = xjm_(0,2)*xjm_(0,2);
  bm(1,2) = xjm_(1,2)*xjm_(1,2);
  bm(2,2) = xjm_(2,2)*xjm_(2,2);
  bm(3,2) = xjm_(0,2)*xjm_(1,2);
  bm(4,2) = xjm_(0,2)*xjm_(2,2);
  bm(5,2) = xjm_(2,2)*xjm_(1,2);

  bm(0,3) = 2.*xjm_(0,0)*xjm_(0,1);
  bm(1,3) = 2.*xjm_(1,0)*xjm_(1,1);
  bm(2,3) = 2.*xjm_(2,0)*xjm_(2,1);
  bm(3,3) = xjm_(0,0)*xjm_(1,1)+xjm_(1,0)*xjm_(0,1);
  bm(4,3) = xjm_(0,0)*xjm_(2,1)+xjm_(2,0)*xjm_(0,1);
  bm(5,3) = xjm_(1,0)*xjm_(2,1)+xjm_(2,0)*xjm_(1,1);

  bm(0,4) = 2.*xjm_(0,0)*xjm_(0,2);
  bm(1,4) = 2.*xjm_(1,0)*xjm_(1,2);
  bm(2,4) = 2.*xjm_(2,0)*xjm_(2,2);
  bm(3,4) = xjm_(0,0)*xjm_(1,2)+xjm_(1,0)*xjm_(0,2);
  bm(4,4) = xjm_(0,0)*xjm_(2,2)+xjm_(2,0)*xjm_(0,2);
  bm(5,4) = xjm_(1,0)*xjm_(2,2)+xjm_(2,0)*xjm_(1,2);

  bm(0,5) = 2.*xjm_(0,1)*xjm_(0,2);
  bm(1,5) = 2.*xjm_(1,1)*xjm_(1,2);
  bm(2,5) = 2.*xjm_(2,1)*xjm_(2,2);
  bm(3,5) = xjm_(0,1)*xjm_(1,2)+xjm_(1,1)*xjm_(0,2);
  bm(4,5) = xjm_(0,1)*xjm_(2,2)+xjm_(2,1)*xjm_(0,2);
  bm(5,5) = xjm_(1,1)*xjm_(2,2)+xjm_(2,1)*xjm_(1,2);

  /*------------------ determine 2nd derivatives of coord.-functions */

  /*
  |
  |         0 1 2              0...iel-1
  |        +-+-+-+             +-+-+-+-+        0 1 2
  |        | | | | 0           | | | | | 0     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | 0
  |        | | | | 1           | | | | | 1   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 2           | | | | | 2     +-+-+-+
  |        +-+-+-+       =     +-+-+-+-+       | | | | .
  |        | | | | 3           | | | | | 3     +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 4           | | | | | 4   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 5           | | | | | 5     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | iel-1
  |		     	      	     	       +-+-+-+
  |
  |        xder2               deriv2          xyze^T
  |
  |
  |                                     +-                  -+
  |  	   	    	    	        | d^2x   d^2y   d^2z |
  |  	   	    	    	        | ----   ----   ---- |
  | 	   	   	   	        | dr^2   dr^2   dr^2 |
  | 	   	   	   	        |                    |
  | 	   	   	   	        | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  | 	   	   	   	        | ds^2   ds^2   ds^2 |
  | 	   	   	   	        |                    |
  | 	   	   	   	        | d^2x   d^2y   d^2z |
  | 	   	   	   	        | ----   ----   ---- |
  | 	   	   	   	        | dt^2   dt^2   dt^2 |
  |               yields    xder2  =    |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drds   drds   drds |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drdt   drdt   drdt |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dsdt   dsdt   dsdt |
  | 	   	   	   	        +-                  -+
  |
  |
  */

  xder2_ = blitz::sum(deriv2_(i,k)*xyze_(j,k),k);

  /*
  |        0...iel-1             0 1 2
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 0          | | | | 0
  |        +-+-+-+-+            +-+-+-+            0...iel-1
  |        | | | | | 1          | | | | 1         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 0
  |        | | | | | 2          | | | | 2         +-+-+-+-+
  |        +-+-+-+-+       =    +-+-+-+       *   | | | | | 1 * (-1)
  |        | | | | | 3          | | | | 3         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 2
  |        | | | | | 4          | | | | 4         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 5          | | | | 5          derxy
  |        +-+-+-+-+            +-+-+-+
  |
  |       chainrulerhs          xder2
  */

  derxy2_ = -blitz::sum(xder2_(i,k)*derxy_(k,j),k);

  /*
  |        0...iel-1            0...iel-1         0...iel-1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 0          | | | | | 0       | | | | | 0
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 1          | | | | | 1       | | | | | 1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 2          | | | | | 2       | | | | | 2
  |        +-+-+-+-+       =    +-+-+-+-+    +    +-+-+-+-+
  |        | | | | | 3          | | | | | 3       | | | | | 3
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 4          | | | | | 4       | | | | | 4
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 5          | | | | | 5       | | | | | 5
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |
  |       chainrulerhs         chainrulerhs        deriv2
  */

  derxy2_ += deriv2_;

  /* make LR decomposition and solve system for all right hand sides
   * (i.e. the components of chainrulerhs)
  |
  |          0  1  2  3  4  5         i        i
  | 	   +--+--+--+--+--+--+       +-+      +-+
  | 	   |  |  |  |  |  |  | 0     | | 0    | | 0
  | 	   +--+--+--+--+--+--+       +-+      +-+
  | 	   |  |  |  |  |  |  | 1     | | 1    | | 1
  | 	   +--+--+--+--+--+--+       +-+      +-+
  | 	   |  |  |  |  |  |  | 2     | | 2    | | 2
  | 	   +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
  |        |  |  |  |  |  |  | 3     | | 3    | | 3
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 4     | | 4    | | 4
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 5     | | 5    | | 5
  |        +--+--+--+--+--+--+       +-+      +-+
  |                                   |        |
  |                                   |        |
  |                                   derxy2[i]|
  |		                               |
  |		                               chainrulerhs[i]
  |
  |	  yields
  |
  |                      0...iel-1
  |                      +-+-+-+-+
  |                      | | | | | 0 = drdr
  |                      +-+-+-+-+
  |                      | | | | | 1 = dsds
  |                      +-+-+-+-+
  |                      | | | | | 2 = dtdt
  |            derxy2 =  +-+-+-+-+
  |                      | | | | | 3 = drds
  |                      +-+-+-+-+
  |                      | | | | | 4 = drdt
  |                      +-+-+-+-+
  |                      | | | | | 5 = dsdt
  |    	          	 +-+-+-+-+
  */

  Epetra_SerialDenseMatrix ederxy2(View,derxy2_.data(),6,6,iel_);

  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(bm);

  // No need for a separate rhs. We assemble the rhs to the solution
  // vector. The solver will destroy the rhs and return the solution.
  solver.SetVectors(ederxy2,ederxy2);
  solver.Solve();

  return;
}


#endif
#endif
