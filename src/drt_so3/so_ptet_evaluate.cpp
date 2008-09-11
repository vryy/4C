/*!----------------------------------------------------------------------*
\file so_ptet_evaluate.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_ptet.H"
#include "so_integrator.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/hyperpolyconvex.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::InitElement(DRT::ELEMENTS::PtetRegister* myregister)
{
  myregister_ = myregister;

  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,NUMDIM_PTET> xrefe;
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,NUMDIM_PTET+1> J;
  {
    // compute element volume
    DRT::Node** nodes = Nodes();
    for (int i=0; i<NUMNOD_PTET; ++i)
    {
      const double* x = nodes[i]->X();
      J(i,0) = 1.0;
      J(i,1) = xrefe(i,0) = x[0];
      J(i,2) = xrefe(i,1) = x[1];
      J(i,3) = xrefe(i,2) = x[2];
    }
    V_ = J.Determinant()/6.0;
    if (V_==0.0)     dserror("Element volume is zero");
    else if (V_<0.0) dserror("Element volume is negative");
  }
  // compute derivatives of shape functions w.r.t. to material coords.
  {
    // one gauss point at 0.25/0.25/0.25
    // gauss point weight is 1.0, so skip it
    const double gploc = 0.25;
    LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,1> funct;
    ShapeFunction(funct,gploc,gploc,gploc,gploc);
    LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,NUMCOORD_PTET> deriv;
    ShapeFunctionDerivatives(deriv);
    LINALG::FixedSizeSerialDenseMatrix<NUMCOORD_PTET-1,NUMCOORD_PTET> tmp;
    tmp.MultiplyTN(xrefe,deriv);
    for (int i=0; i<4; i++) J(0,i)=1;
    for (int row=0;row<3;row++)
      for (int col=0;col<4;col++)
        J(row+1,col)=tmp(row,col);

    LINALG::FixedSizeSerialDenseMatrix<NUMCOORD_PTET,NUMDIM_PTET> Iaug(true); // initialize to zero
    Iaug(1,0)=1;
    Iaug(2,1)=1;
    Iaug(3,2)=1;
    LINALG::FixedSizeSerialDenseMatrix<NUMCOORD_PTET,NUMDIM_PTET> partials;
    LINALG::FixedSizeSerialDenseSolver<NUMCOORD_PTET,NUMCOORD_PTET,NUMDIM_PTET> solver;
    solver.SetMatrix(J);
    solver.SetVectors(partials,Iaug);
    solver.FactorWithEquilibration(true);
    int err  = solver.Factor();
    int err2 = solver.Solve();
    if (err || err2) dserror("Inversion of Jacobian failed");
    /* structure of nxyz_:
    **             [   dN_1     dN_1     dN_1   ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **    nxyz_ =  [     |        |        |    ]
    **             [                            ]
    **             [   dN_4     dN_4     dN_4   ]
    **             [  -------  -------  ------- ]
    **             [    dX       dY       dZ    ]
    */
    nxyz_.Multiply(deriv,partials);
  }


  return;
} // DRT::ELEMENTS::Ptet::InitElement



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Ptet::Evaluate(ParameterList& params,
                                  DRT::Discretization&      discretization,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseMatrix& elemat1_epetra,
                                  Epetra_SerialDenseMatrix& elemat2_epetra,
                                  Epetra_SerialDenseVector& elevec1_epetra,
                                  Epetra_SerialDenseVector& elevec2_epetra,
                                  Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,NUMDOF_PTET> elemat1(elemat1_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,NUMDOF_PTET> elemat2(elemat2_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,          1> elevec1(elevec1_epetra.A(),true);
  LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,          1> elevec2(elevec2_epetra.A(),true);
  // start with "none"
  DRT::ELEMENTS::Ptet::ActionType act = Ptet::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                act = Ptet::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                act = Ptet::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")           act = Ptet::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")            act = Ptet::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")            act = Ptet::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")           act = Ptet::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                  act = Ptet::calc_struct_stress;
  else if (action=="postprocess_stress")                  act = Ptet::postprocess_stress;
  else if (action=="calc_struct_eleload")                 act = Ptet::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")                 act = Ptet::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")            act = Ptet::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")          act = Ptet::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")             act = Ptet::calc_struct_reset_istep;
  else dserror("Unknown type of action for Ptet");

  // what should the element do
  switch(act)
  {
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      ptetnlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1);
      if (act==calc_struct_nlnstifflmass) ptetlumpmass(&elemat2);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,NUMDOF_PTET>* elemat1ptr = NULL;
      if (elemat1.IsInitialized()) elemat1ptr = &elemat1;
      ptetnlnstiffmass(lm,mydisp,myres,elemat1ptr,NULL,&elevec1);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // The element does not calc stresses, it just takes the nodal
      // stresses/strains that have been computed by the register class
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (stressdata==null) dserror("Cannot get stress 'data'");
      if (straindata==null) dserror("Cannot get strain 'data'");
      LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,NUMSTR_PTET> stress;
      LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,NUMSTR_PTET> strain;
      map<int,vector<double> >& nodestress = myregister_->nodestress_;
      map<int,vector<double> >& nodestrain = myregister_->nodestrain_;
      for (int i=0; i<NumNode(); ++i)
      {
        int gid = Nodes()[i]->Id();
        map<int,vector<double> >::iterator foolstress = nodestress.find(gid);
        map<int,vector<double> >::iterator foolstrain = nodestrain.find(gid);
        if (foolstress==nodestress.end() || foolstrain==nodestrain.end())
          dserror("Cannot find marching nodal stresses/strains");
        vector<double>& nstress = foolstress->second;
        vector<double>& nstrain = foolstrain->second;
        for (int j=0; j<6; ++j)
        {
          stress(i,j) = nstress[j];
          strain(i,j) = nstrain[j];
        }
      }
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,NUMSTR_PTET> gpstress(((*gpstressmap)[gid])->A(),true);
      if (stresstype=="ndxyz")
      {
        for (int i=0;i<NUMNOD_PTET;++i)
        {
          elevec1(3*i)=gpstress(i,0);
          elevec1(3*i+1)=gpstress(i,1);
          elevec1(3*i+2)=gpstress(i,2);
        }
        for (int i=0;i<NUMNOD_PTET;++i)
        {
          elevec2(3*i)=gpstress(i,3);
          elevec2(3*i+1)=gpstress(i,4);
          elevec2(3*i+2)=gpstress(i,5);
        }
      }
      else if (stresstype=="cxyz" || stresstype=="cxyz_ndxyz")
        dserror("The Ptet does not do element stresses, nodal only (ndxyz), because its a nodal tet!");
      else
        dserror("unknown type of stress/strain output on element level");
    }
    break;

    case calc_struct_eleload:
      dserror("this class is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_update_imrlike:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_reset_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,NUMDOF_PTET> myemat(true);
      ptetnlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not implemented");
    break;

    // linear stiffness
    case calc_struct_linstiff:
      dserror("action calc_struct_linstiff currently not supported");
    break;

    default:
      dserror("Unknown type of action for Ptet");
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::ptetnlnstiffmass(
      vector<int>&              lm,                                                // location matrix
      vector<double>&           disp,                                              // current displacements
      vector<double>&           residual,                                          // current residuum
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,NUMDOF_PTET>* stiffmatrix,    // element stiffness matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,NUMDOF_PTET>* massmatrix,     // element mass matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,          1>* force)          // stress output options
{
  //--------------------------------------------------- geometry update
  if (!FisNew_) DeformationGradient(disp);
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_PTET,NUMDIM_PTET>& defgrd = F_;
  // reset the bool indicating that the stored deformation gradient is 'fresh'
  FisNew_ = false;

  //--------------------------- Right Cauchy-Green tensor C = = F^T * F
  LINALG::FixedSizeSerialDenseMatrix<NUMDIM_PTET,NUMDIM_PTET> cauchygreen;
  cauchygreen.MultiplyTN(defgrd,defgrd);

  // --Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_PTET,1> glstrain;
  glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
  glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
  glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
  glstrain(3) = cauchygreen(0,1);
  glstrain(4) = cauchygreen(1,2);
  glstrain(5) = cauchygreen(2,0);

  //------------------------------------ B=operator ( same as in hex8 !)
  /*
  ** B = F : N,xyz
  **
  **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
  **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
  **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
  ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
  **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
  **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
  **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
  **      [                                                         ]
  **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
  **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
  **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
  **      [                                                         ]
  **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
  **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
  **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
  */

  // 6x12 n_stresses * number degrees of freedom per element
  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_PTET,NUMDOF_PTET> bop;
  for (int i=0; i<NUMNOD_PTET; i++)
  {
    bop(0,NODDOF_PTET*i+0) = defgrd(0,0)*nxyz_(i,0);
    bop(0,NODDOF_PTET*i+1) = defgrd(1,0)*nxyz_(i,0);
    bop(0,NODDOF_PTET*i+2) = defgrd(2,0)*nxyz_(i,0);
    bop(1,NODDOF_PTET*i+0) = defgrd(0,1)*nxyz_(i,1);
    bop(1,NODDOF_PTET*i+1) = defgrd(1,1)*nxyz_(i,1);
    bop(1,NODDOF_PTET*i+2) = defgrd(2,1)*nxyz_(i,1);
    bop(2,NODDOF_PTET*i+0) = defgrd(0,2)*nxyz_(i,2);
    bop(2,NODDOF_PTET*i+1) = defgrd(1,2)*nxyz_(i,2);
    bop(2,NODDOF_PTET*i+2) = defgrd(2,2)*nxyz_(i,2);
    /* ~~~ */
    bop(3,NODDOF_PTET*i+0) = defgrd(0,0)*nxyz_(i,1) + defgrd(0,1)*nxyz_(i,0);
    bop(3,NODDOF_PTET*i+1) = defgrd(1,0)*nxyz_(i,1) + defgrd(1,1)*nxyz_(i,0);
    bop(3,NODDOF_PTET*i+2) = defgrd(2,0)*nxyz_(i,1) + defgrd(2,1)*nxyz_(i,0);
    bop(4,NODDOF_PTET*i+0) = defgrd(0,1)*nxyz_(i,2) + defgrd(0,2)*nxyz_(i,1);
    bop(4,NODDOF_PTET*i+1) = defgrd(1,1)*nxyz_(i,2) + defgrd(1,2)*nxyz_(i,1);
    bop(4,NODDOF_PTET*i+2) = defgrd(2,1)*nxyz_(i,2) + defgrd(2,2)*nxyz_(i,1);
    bop(5,NODDOF_PTET*i+0) = defgrd(0,2)*nxyz_(i,0) + defgrd(0,0)*nxyz_(i,2);
    bop(5,NODDOF_PTET*i+1) = defgrd(1,2)*nxyz_(i,0) + defgrd(1,0)*nxyz_(i,2);
    bop(5,NODDOF_PTET*i+2) = defgrd(2,2)*nxyz_(i,0) + defgrd(2,0)*nxyz_(i,2);
  }

  //------------------------------------------------- call material law
#if 1 // dev stab on cauchy stresses
  Epetra_SerialDenseMatrix cmat_epetra(NUMSTR_PTET,NUMSTR_PTET);
  Epetra_SerialDenseVector stress_epetra(NUMSTR_PTET);
  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_PTET,NUMSTR_PTET> cmat(cmat_epetra.A(),true); // Views
  LINALG::FixedSizeSerialDenseMatrix<NUMSTR_PTET,1> stress(stress_epetra.A(),true);
  double density = -999.99;
  {
    // do deviatoric F, C, E
    const double J = defgrd.Determinant();
    LINALG::FixedSizeSerialDenseMatrix<NUMDIM_PTET,NUMDIM_PTET> Cbar(cauchygreen);
    Cbar.Scale(pow(J,-2./3.));
    LINALG::SerialDenseVector glstrainbar(NUMSTR_PTET);
    glstrainbar(0) = 0.5 * (Cbar(0,0) - 1.0);
    glstrainbar(1) = 0.5 * (Cbar(1,1) - 1.0);
    glstrainbar(2) = 0.5 * (Cbar(2,2) - 1.0);
    glstrainbar(3) = Cbar(0,1);
    glstrainbar(4) = Cbar(1,2);
    glstrainbar(5) = Cbar(2,0);
    Epetra_SerialDenseMatrix Fbar(Copy,defgrd.A(),NUMDIM_PTET,NUMDIM_PTET,NUMDIM_PTET);
    Fbar.Scale(pow(J,-1./3.));

    SelectMaterial(stress_epetra,cmat_epetra,density,glstrainbar,Fbar,0);

    // define stuff we need to do the split
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_PTET,NUMSTR_PTET> cmatdev; // set to zero?
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_PTET,1> stressdev;

    // do just the deviatoric components
    PtetRegister::DevStressTangent(stressdev,cmatdev,cmat,stress,cauchygreen);
    stress = stressdev;
    cmat = cmatdev;
    stress.Scale(ALPHA_PTET);
    cmat.Scale(ALPHA_PTET);
  }
#endif

#if 0 // original puso tet
  Epetra_SerialDenseMatrix cmat(NUMSTR_PTET,NUMSTR_PTET);
  Epetra_SerialDenseVector stress(NUMSTR_PTET);
  double density = -999.99;
  SelectMaterial(stress,cmat,density,glstrain,defgrd,0);
  stress.Scale(ALPHA_PTET);
  cmat.Scale(ALPHA_PTET);
#endif

  if (force)
  {
    // integrate internal force vector f = f + (B^T . sigma) * V_
    force->MultiplyTN(V_,bop,stress,1.0);
  }
  if (stiffmatrix)
  {
    // integrate elastic stiffness matrix
    // keu = keu + (B^T . C . B) * V_
    LINALG::FixedSizeSerialDenseMatrix<NUMSTR_PTET,NUMDOF_PTET> cb;
    cb.Multiply(cmat,bop);
    stiffmatrix->MultiplyTN(V_,bop,cb,1.0);
    // integrate `geometric' stiffness matrix and add to keu
    double sBL[3];
    const double V = Volume();
    for (int i=0; i<NUMNOD_PTET; ++i)
    {
      sBL[0] = V*(stress(0) * nxyz_(i,0) + stress(3) * nxyz_(i,1) + stress(5) * nxyz_(i,2));
      sBL[1] = V*(stress(3) * nxyz_(i,0) + stress(1) * nxyz_(i,1) + stress(4) * nxyz_(i,2));
      sBL[2] = V*(stress(5) * nxyz_(i,0) + stress(4) * nxyz_(i,1) + stress(2) * nxyz_(i,2));
      for (int j=0; j<NUMNOD_PTET; ++j)
      {
        double BsB = 0.0;
        for (int dim=0; dim<NUMDIM_PTET; ++dim)
          BsB += nxyz_(j,dim) * sBL[dim];
        (*stiffmatrix)(NUMDIM_PTET*i+0,NUMDIM_PTET*j+0) += BsB;
        (*stiffmatrix)(NUMDIM_PTET*i+1,NUMDIM_PTET*j+1) += BsB;
        (*stiffmatrix)(NUMDIM_PTET*i+2,NUMDIM_PTET*j+2) += BsB;
      }
    }

  } // if (force && stiffmatrix)

  //------------------------------------------ do mass matrix if desired
  if (massmatrix)
  {
    // for mass matrix use a 4 gauss points integration:
    // ( 1 gauss point is not enough!)
    const double alpha  = (5.0 + 3.0*sqrt(5.0))/20.0;
    const double beta   = (5.0 - sqrt(5.0))/20.0;
    const double weight = 0.25;
    const double V = Volume();
    double xsi[4][4];
    xsi[0][0] = alpha;   xsi[0][1] = beta ;   xsi[0][2] = beta ;   xsi[0][3] = beta ;
    xsi[1][0] = beta ;   xsi[1][1] = alpha;   xsi[1][2] = beta ;   xsi[1][3] = beta ;
    xsi[2][0] = beta ;   xsi[2][1] = beta ;   xsi[2][2] = alpha;   xsi[2][3] = beta ;
    xsi[3][0] = beta ;   xsi[3][1] = beta ;   xsi[3][2] = beta ;   xsi[3][3] = alpha;
    for (int gp=0; gp<4; ++gp)
    {
      LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,1> funct;
      ShapeFunction(funct,xsi[gp][0],xsi[gp][1],xsi[gp][2],xsi[gp][3]);
      const double f = density * V * weight;
      for (int i=0; i<NUMNOD_PTET; ++i)
        for (int j=0; j<NUMNOD_PTET; ++j)
        {
          const double fac = funct(i) * funct(j) * f;
          (*massmatrix)(NUMDIM_PTET*i+0,NUMDIM_PTET*j+0) += fac;
          (*massmatrix)(NUMDIM_PTET*i+1,NUMDIM_PTET*j+1) += fac;
          (*massmatrix)(NUMDIM_PTET*i+2,NUMDIM_PTET*j+2) += fac;
        }
    }
  } // if (massmatrix)

  return;
} // DRT::ELEMENTS::Ptet::ptetnlnstiffmass


/*----------------------------------------------------------------------*
 |  lump mass matrix                                         bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::ptetlumpmass(LINALG::FixedSizeSerialDenseMatrix<NUMDOF_PTET,NUMDOF_PTET>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 |  compute deformation gradient (protected)                   gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::DeformationGradient(vector<double>& disp)
{
  LINALG::FixedSizeSerialDenseMatrix<NUMNOD_PTET,NUMDIM_PTET> xdisp;
  for (int i=0; i<NUMNOD_PTET; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_PTET+0];
    xdisp(i,1) = disp[i*NODDOF_PTET+1];
    xdisp(i,2) = disp[i*NODDOF_PTET+2];
  }
  F_.MultiplyTN(xdisp,nxyz_);
  F_(0,0)+=1;
  F_(1,1)+=1;
  F_(2,2)+=1;
  FisNew_ = true;
  return;
}


/*----------------------------------------------------------------------*
 | material laws for Ptet (protected)                          gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::SelectMaterial(
                                Epetra_SerialDenseVector& stress,
                                Epetra_SerialDenseMatrix& cmat,
                                double& density,
                                const Epetra_SerialDenseVector& glstrain,
                                const Epetra_SerialDenseMatrix& defgrd,
                                int gp)
{
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast<MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(&glstrain,&cmat,&stress);
      density = stvk->Density();
    }
    break;
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast<MAT::NeoHooke*>(mat.get());
      neo->Evaluate(&glstrain,&cmat,&stress);
      density = neo->Density();
    }
    break;
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast<MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(&glstrain,&cmat,&stress);
      density = aaa->Density();
    }
    break;
    default:
      dserror("Illegal type %d of material for element Ptet tet4", mat->MaterialType());
    break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // DRT::ELEMENTS::Ptet::SelectMaterial


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Ptet::EvaluateNeumann(ParameterList& params,
                                         DRT::Discretization&      discretization,
                                         DRT::Condition&           condition,
                                         vector<int>&              lm,
                                         Epetra_SerialDenseVector& elevec1)
{
  dserror("DRT::ELEMENTS::Ptet::EvaluateNeumann not implemented");
  return -1;
} // DRT::ELEMENTS::Ptet::EvaluateNeumann







#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
