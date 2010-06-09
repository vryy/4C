/*!----------------------------------------------------------------------*
\file so_nstet_evaluate.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_nstet.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::InitElement()
{
  LINALG::Matrix<NUMNOD_NSTET,NUMDIM_NSTET> xrefe;
  LINALG::Matrix<NUMNOD_NSTET,NUMDIM_NSTET+1> J;
  {
    // compute element volume
    DRT::Node** nodes = Nodes();
    for (int i=0; i<NUMNOD_NSTET; ++i)
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
    LINALG::Matrix<NUMNOD_NSTET,1> funct;
    ShapeFunction(funct,gploc,gploc,gploc,gploc);
    LINALG::Matrix<NUMNOD_NSTET,NUMCOORD_NSTET> deriv;
    ShapeFunctionDerivatives(deriv);
    LINALG::Matrix<NUMCOORD_NSTET-1,NUMCOORD_NSTET> tmp;
    tmp.MultiplyTN(xrefe,deriv);
    for (int i=0; i<4; i++) J(0,i)=1;
    for (int row=0;row<3;row++)
      for (int col=0;col<4;col++)
        J(row+1,col)=tmp(row,col);

    LINALG::Matrix<NUMCOORD_NSTET,NUMDIM_NSTET> Iaug(true); // initialize to zero
    Iaug(1,0)=1;
    Iaug(2,1)=1;
    Iaug(3,2)=1;
    LINALG::Matrix<NUMCOORD_NSTET,NUMDIM_NSTET> partials;
    LINALG::FixedSizeSerialDenseSolver<NUMCOORD_NSTET,NUMCOORD_NSTET,NUMDIM_NSTET> solver;
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
} // DRT::ELEMENTS::NStet::InitElement



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NStet::Evaluate(ParameterList& params,
                                  DRT::Discretization&      discretization,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseMatrix& elemat1_epetra,
                                  Epetra_SerialDenseMatrix& elemat2_epetra,
                                  Epetra_SerialDenseVector& elevec1_epetra,
                                  Epetra_SerialDenseVector& elevec2_epetra,
                                  Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_NSTET,          1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_NSTET,          1> elevec2(elevec2_epetra.A(),true);
  // start with "none"
  DRT::ELEMENTS::NStet::ActionType act = NStet::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                act = NStet::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                act = NStet::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")           act = NStet::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")            act = NStet::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")            act = NStet::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")           act = NStet::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                  act = NStet::calc_struct_stress;
  else if (action=="postprocess_stress")                  act = NStet::postprocess_stress;
  else if (action=="calc_struct_eleload")                 act = NStet::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")                 act = NStet::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")            act = NStet::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")          act = NStet::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")             act = NStet::calc_struct_reset_istep;
  else dserror("Unknown type of action for NStet");

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
      nstetnlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1);
      if (act==calc_struct_nlnstifflmass) nstetlumpmass(&elemat2);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET>* elemat1ptr = NULL;
      if (elemat1.IsInitialized()) elemat1ptr = &elemat1;
      nstetnlnstiffmass(lm,mydisp,myres,elemat1ptr,NULL,&elevec1);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
#if 1 // original nodal strain stresses
    {
      // The element does not calc stresses, it just takes the nodal
      // stresses/strains that have been computed by the register class
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (stressdata==null) dserror("Cannot get stress 'data'");
      if (straindata==null) dserror("Cannot get strain 'data'");
      LINALG::Matrix<NUMNOD_NSTET,NUMSTR_NSTET> stress;
      LINALG::Matrix<NUMNOD_NSTET,NUMSTR_NSTET> strain;
      map<int,vector<double> >& nodestress = ElementType().nodestress_;
      map<int,vector<double> >& nodestrain = ElementType().nodestrain_;
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
#else // classical tet4 stresses
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==null) dserror("Cannot get 'stress' data");
      if (straindata==null) dserror("Cannot get 'strain' data");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET> stress(true); // set to zero
      LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET> strain(true);
      INPAR::STR::StressType iostress = params.get<INPAR::STR::StressType>("iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain = params.get<INPAR::STR::StrainType>("iostrain", INPAR::STR::strain_none);
      calcstress(lm,mydisp,&stress,&strain,iostress,iostrain);
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
#endif
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
#if 1 // nodal strain style stresses
    {
      const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMNOD_NSTET,NUMSTR_NSTET> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz")
      {
        for (int i=0;i<NUMNOD_NSTET;++i)
        {
          elevec1(3*i)=gpstress(i,0);
          elevec1(3*i+1)=gpstress(i,1);
          elevec1(3*i+2)=gpstress(i,2);
        }
        for (int i=0;i<NUMNOD_NSTET;++i)
        {
          elevec2(3*i)=gpstress(i,3);
          elevec2(3*i+1)=gpstress(i,4);
          elevec2(3*i+2)=gpstress(i,5);
        }
      }
      else if (stresstype=="cxyz" || stresstype=="cxyz_ndxyz")
        dserror("The NStet does not do element stresses, nodal only (ndxyz), because its a nodal tet!");
      else
        dserror("unknown type of stress/strain output on element level");
    }
#else
    {
      const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        so_tet4_expol(gpstress, elevec1, elevec2);
      }
      else if (stresstype=="cxyz")
      {
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < NUMSTR_NSTET; ++i)
          {
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_NSTET; ++j)
              s += gpstress(j,i);
            s /= NUMGPT_NSTET;
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        so_tet4_expol(gpstress, elevec1, elevec2);

        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_NSTET; ++i)
          {
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_NSTET; ++j)
              s += gpstress(j,i);
            s /= NUMGPT_NSTET;
          }
        }
      }
      else {
        dserror("unknown type of stress/strain output on element level");
      }
    }
#endif
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
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET> myemat(true);
      nstetnlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1);
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
      dserror("Unknown type of action for NStet");
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::nstetnlnstiffmass(
      vector<int>&              lm,                              // location matrix
      vector<double>&           disp,                            // current displacements
      vector<double>&           residual,                        // current residuum
      LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET>* stiffmatrix,    // element stiffness matrix
      LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET>* massmatrix,     // element mass matrix
      LINALG::Matrix<NUMDOF_NSTET,          1>* force)           // stress output options
{
  //--------------------------------------------------- geometry update
  if (!FisNew_) DeformationGradient(disp);
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd = F_;

  // reset the bool indicating that the stored deformation gradient is 'fresh'
  FisNew_ = false;

  //--------------------------- Right Cauchy-Green tensor C = = F^T * F
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> cauchygreen;
  cauchygreen.MultiplyTN(defgrd,defgrd);

  // --Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  LINALG::Matrix<NUMSTR_NSTET,1> glstrain;
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
  LINALG::Matrix<NUMSTR_NSTET,NUMDOF_NSTET> bop;
  for (int i=0; i<NUMNOD_NSTET; i++)
  {
    bop(0,NODDOF_NSTET*i+0) = defgrd(0,0)*nxyz_(i,0);
    bop(0,NODDOF_NSTET*i+1) = defgrd(1,0)*nxyz_(i,0);
    bop(0,NODDOF_NSTET*i+2) = defgrd(2,0)*nxyz_(i,0);
    bop(1,NODDOF_NSTET*i+0) = defgrd(0,1)*nxyz_(i,1);
    bop(1,NODDOF_NSTET*i+1) = defgrd(1,1)*nxyz_(i,1);
    bop(1,NODDOF_NSTET*i+2) = defgrd(2,1)*nxyz_(i,1);
    bop(2,NODDOF_NSTET*i+0) = defgrd(0,2)*nxyz_(i,2);
    bop(2,NODDOF_NSTET*i+1) = defgrd(1,2)*nxyz_(i,2);
    bop(2,NODDOF_NSTET*i+2) = defgrd(2,2)*nxyz_(i,2);
    /* ~~~ */
    bop(3,NODDOF_NSTET*i+0) = defgrd(0,0)*nxyz_(i,1) + defgrd(0,1)*nxyz_(i,0);
    bop(3,NODDOF_NSTET*i+1) = defgrd(1,0)*nxyz_(i,1) + defgrd(1,1)*nxyz_(i,0);
    bop(3,NODDOF_NSTET*i+2) = defgrd(2,0)*nxyz_(i,1) + defgrd(2,1)*nxyz_(i,0);
    bop(4,NODDOF_NSTET*i+0) = defgrd(0,1)*nxyz_(i,2) + defgrd(0,2)*nxyz_(i,1);
    bop(4,NODDOF_NSTET*i+1) = defgrd(1,1)*nxyz_(i,2) + defgrd(1,2)*nxyz_(i,1);
    bop(4,NODDOF_NSTET*i+2) = defgrd(2,1)*nxyz_(i,2) + defgrd(2,2)*nxyz_(i,1);
    bop(5,NODDOF_NSTET*i+0) = defgrd(0,2)*nxyz_(i,0) + defgrd(0,0)*nxyz_(i,2);
    bop(5,NODDOF_NSTET*i+1) = defgrd(1,2)*nxyz_(i,0) + defgrd(1,0)*nxyz_(i,2);
    bop(5,NODDOF_NSTET*i+2) = defgrd(2,2)*nxyz_(i,0) + defgrd(2,0)*nxyz_(i,2);
  }

  //------------------------------------------------- call material law
  //------------------------------------------------------ stabilization
  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET> cmat(true);
  LINALG::Matrix<NUMSTR_NSTET,1> stress(true);
  double density = -999.99;
  switch(stabtype_)
  {
    case DRT::ELEMENTS::so_nstet4_voldev:
    {
      VolDevStab(defgrd,cauchygreen,stress,cmat,density);
      if (force) VolDevStabLinear(defgrd,force);
    }
    break;
    case DRT::ELEMENTS::so_nstet4_dev:
    {
      DevStab(defgrd,cauchygreen,stress,cmat,density);
    }
    break;
    case DRT::ELEMENTS::so_nstet4_puso:
    {
      SelectMaterial(stress,cmat,density,glstrain,defgrd,0);
      stress.Scale(ALPHA_NSTET);
      cmat.Scale(ALPHA_NSTET);
    }
    break;
    case DRT::ELEMENTS::so_nstet4_stab_none:
    break;
    default:
      dserror("Unknown type of stabilization");
    break;
  }

  //-------------------------------------------------------------------------------
  if (force)
  {
    // integrate internal force vector f = f + (B^T . sigma) * V_
    force->MultiplyTN(V_,bop,stress,1.0);
  }
  if (stiffmatrix)
  {
    // integrate elastic stiffness matrix
    // keu = keu + (B^T . C . B) * V_
    LINALG::Matrix<NUMSTR_NSTET,NUMDOF_NSTET> cb;
    cb.Multiply(cmat,bop);
    stiffmatrix->MultiplyTN(V_,bop,cb,1.0);
    // integrate `geometric' stiffness matrix and add to keu
    double sBL[3];
    const double V = Volume();
    for (int i=0; i<NUMNOD_NSTET; ++i)
    {
      sBL[0] = V*(stress(0) * nxyz_(i,0) + stress(3) * nxyz_(i,1) + stress(5) * nxyz_(i,2));
      sBL[1] = V*(stress(3) * nxyz_(i,0) + stress(1) * nxyz_(i,1) + stress(4) * nxyz_(i,2));
      sBL[2] = V*(stress(5) * nxyz_(i,0) + stress(4) * nxyz_(i,1) + stress(2) * nxyz_(i,2));
      for (int j=0; j<NUMNOD_NSTET; ++j)
      {
        double BsB = 0.0;
        for (int dim=0; dim<NUMDIM_NSTET; ++dim)
          BsB += nxyz_(j,dim) * sBL[dim];
        (*stiffmatrix)(NUMDIM_NSTET*i+0,NUMDIM_NSTET*j+0) += BsB;
        (*stiffmatrix)(NUMDIM_NSTET*i+1,NUMDIM_NSTET*j+1) += BsB;
        (*stiffmatrix)(NUMDIM_NSTET*i+2,NUMDIM_NSTET*j+2) += BsB;
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
    const double f = density * V * weight;
    for (int gp=0; gp<4; ++gp)
    {
      LINALG::Matrix<NUMNOD_NSTET,1> funct;
      ShapeFunction(funct,xsi[gp][0],xsi[gp][1],xsi[gp][2],xsi[gp][3]);
      for (int i=0; i<NUMNOD_NSTET; ++i)
        for (int j=0; j<NUMNOD_NSTET; ++j)
        {
          const double fac = funct(i) * funct(j) * f;
          (*massmatrix)(NUMDIM_NSTET*i+0,NUMDIM_NSTET*j+0) += fac;
          (*massmatrix)(NUMDIM_NSTET*i+1,NUMDIM_NSTET*j+1) += fac;
          (*massmatrix)(NUMDIM_NSTET*i+2,NUMDIM_NSTET*j+2) += fac;
        }
    }
  } // if (massmatrix)

  return;
} // DRT::ELEMENTS::NStet::nstetnlnstiffmass


/*----------------------------------------------------------------------*
 |  lump mass matrix                                         bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::nstetlumpmass(LINALG::Matrix<NUMDOF_NSTET,NUMDOF_NSTET>* emass)
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
void DRT::ELEMENTS::NStet::DeformationGradient(vector<double>& disp)
{
  LINALG::Matrix<NUMNOD_NSTET,NUMDIM_NSTET> xdisp;
  for (int i=0; i<NUMNOD_NSTET; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_NSTET+0];
    xdisp(i,1) = disp[i*NODDOF_NSTET+1];
    xdisp(i,2) = disp[i*NODDOF_NSTET+2];
  }
  F_.MultiplyTN(xdisp,nxyz_);
  F_(0,0)+=1;
  F_(1,1)+=1;
  F_(2,2)+=1;
  FisNew_ = true;
  return;
}

#if 0
/*----------------------------------------------------------------------*
 |  compute deformed volume (protected)                        gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::DeformedVolume(vector<double>& disp)
{
  LINALG::Matrix<NUMNOD_NSTET,NUMDIM_NSTET+1> J;
  for (int i=0; i<NUMNOD_NSTET; ++i)
  {
    J(i,0) = 1.0;
    J(i,1) = Nodes()[i]->X()[0] + disp[i*NODDOF_NSTET+0];
    J(i,2) = Nodes()[i]->X()[1] + disp[i*NODDOF_NSTET+1];
    J(i,3) = Nodes()[i]->X()[2] + disp[i*NODDOF_NSTET+2];
  }
  v_ = J.Determinant()/6.0;
  if (v_==0.0)     dserror("Element spatial volume is zero");
  else if (v_<0.0) dserror("Element spatial volume is negative");
  visNew_ = true;
  return;
}
#endif

/*----------------------------------------------------------------------*
 | material laws for NStet (protected)                          gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::SelectMaterial(
                      LINALG::Matrix<6,1>& stress,
                      LINALG::Matrix<6,6>& cmat,
                      double& density,
                      LINALG::Matrix<6,1>& glstrain,
                      LINALG::Matrix<3,3>& defgrd,
                      int gp)
{
  Epetra_SerialDenseVector stress_e(View,stress.A(),stress.Rows());
  Epetra_SerialDenseMatrix cmat_e(View,cmat.A(),cmat.Rows(),cmat.Rows(),cmat.Columns());
  const Epetra_SerialDenseVector glstrain_e(View,glstrain.A(),glstrain.Rows());
  //Epetra_SerialDenseMatrix defgrd_e(View,defgrd.A(),defgrd.Rows(),defgrd.Rows(),defgrd.Columns());


  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast<MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(&glstrain_e,&cmat_e,&stress_e);
      density = stvk->Density();
    }
    break;
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast<MAT::NeoHooke*>(mat.get());
      neo->Evaluate(&glstrain_e,&cmat_e,&stress_e);
      density = neo->Density();
    }
    break;
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast<MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(&glstrain_e,&cmat_e,&stress_e);
      density = aaa->Density();
    }
    break;
    case INPAR::MAT::m_lung_ogden: /* lung tissue material with Ogden for volumetric part */
    {
      MAT::LungOgden* lungog = static_cast <MAT::LungOgden*>(mat.get());
      lungog->Evaluate(&glstrain,&cmat,&stress);
      density = lungog->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_penalty: /* lung tissue material with penalty function for incompressibility constraint */
    {
      MAT::LungPenalty* lungpen = static_cast <MAT::LungPenalty*>(mat.get());

      lungpen->Evaluate(&glstrain,&cmat,&stress);

      density = lungpen->Density();
      return;
      break;
    }
    default:
      dserror("Illegal type %d of material for element NStet tet4", mat->MaterialType());
    break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // DRT::ELEMENTS::NStet::SelectMaterial



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NStet::EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1)
{
  dserror("DRT::ELEMENTS::NStet::EvaluateNeumann not implemented");
  return -1;
} // DRT::ELEMENTS::NStet::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  evaluate the element stress                                gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::calcstress(vector<int>&              lm,
                          vector<double>&           disp,
                          LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>* stresses,
                          LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>* elestrain,
                          const INPAR::STR::StressType     iostress,
                          const INPAR::STR::StrainType     iostrain)
{
  //--------------------------------------------------- geometry update
  if (!FisNew_) DeformationGradient(disp);
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd = F_;
  // reset the bool indicating that the stored deformation gradient is 'fresh'
  FisNew_ = false;

  //--------------------------- Right Cauchy-Green tensor C = = F^T * F
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> cauchygreen;
  cauchygreen.MultiplyTN(defgrd,defgrd);

  // --Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  LINALG::Matrix<NUMSTR_NSTET,1> glstrain;
  glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
  glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
  glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
  glstrain(3) = cauchygreen(0,1);
  glstrain(4) = cauchygreen(1,2);
  glstrain(5) = cauchygreen(2,0);

  switch (iostrain)
  {
  case INPAR::STR::strain_gl:
  {
    if (elestrain == NULL) dserror("no strain data available");
    for (int i = 0; i < 3; ++i)
      (*elestrain)(0,i) = glstrain(i);
    for (int i = 3; i < 6; ++i)
      (*elestrain)(0,i) = 0.5 * glstrain(i);
  }
  break;
  case INPAR::STR::strain_ea:
  {
    if (elestrain == NULL) dserror("no strain data available");

    // rewriting Green-Lagrange strains in matrix format
    LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> gl;
    gl(0,0) = glstrain(0);
    gl(0,1) = 0.5*glstrain(3);
    gl(0,2) = 0.5*glstrain(5);
    gl(1,0) = gl(0,1);
    gl(1,1) = glstrain(1);
    gl(1,2) = 0.5*glstrain(4);
    gl(2,0) = gl(0,2);
    gl(2,1) = gl(1,2);
    gl(2,2) = glstrain(2);

    // inverse of deformation gradient
    //Epetra_SerialDenseMatrix invdefgrd(defgrd); // make a copy here otherwise defgrd is destroyed!
    //LINALG::NonsymInverse3x3(invdefgrd);
    LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> invdefgrd;
    invdefgrd.Invert(defgrd);

    LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> temp;
    LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> euler_almansi;
    temp.Multiply(gl,invdefgrd);
    euler_almansi.MultiplyTN(invdefgrd,temp);

    (*elestrain)(0,0) = euler_almansi(0,0);
    (*elestrain)(0,1) = euler_almansi(1,1);
    (*elestrain)(0,2) = euler_almansi(2,2);
    (*elestrain)(0,3) = euler_almansi(0,1);
    (*elestrain)(0,4) = euler_almansi(1,2);
    (*elestrain)(0,5) = euler_almansi(0,2);
  }
  break;
  case INPAR::STR::strain_none:
    break;
  default:
    dserror("requested strain option not available");
  }

  //------------------------------------------------- call material law
  //------------------------------------------------------ stabilization
  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET> cmat(true); // Views
  LINALG::Matrix<NUMSTR_NSTET,1> stress(true);
  double density = -999.99;
  SelectMaterial(stress,cmat,density,glstrain,defgrd,0);

  switch (iostress)
  {
  case INPAR::STR::stress_2pk:
  {
    if (stresses == NULL) dserror("no stress data available");
    for (int i = 0; i < NUMSTR_NSTET; ++i)
      (*stresses)(0,i) = stress(i);
  }
  break;
  case INPAR::STR::stress_cauchy:
  {
    if (stresses == NULL) dserror("no stress data available");
    double detF = defgrd.Determinant();

    LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> pkstress;
    pkstress(0,0) = stress(0);
    pkstress(0,1) = stress(3);
    pkstress(0,2) = stress(5);
    pkstress(1,0) = pkstress(0,1);
    pkstress(1,1) = stress(1);
    pkstress(1,2) = stress(4);
    pkstress(2,0) = pkstress(0,2);
    pkstress(2,1) = pkstress(1,2);
    pkstress(2,2) = stress(2);

    LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> temp;
    LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> cauchystress;
    temp.Multiply(1.0/detF,defgrd,pkstress,0.);
    cauchystress.MultiplyNT(temp,defgrd);

    (*stresses)(0,0) = cauchystress(0,0);
    (*stresses)(0,1) = cauchystress(1,1);
    (*stresses)(0,2) = cauchystress(2,2);
    (*stresses)(0,3) = cauchystress(0,1);
    (*stresses)(0,4) = cauchystress(1,2);
    (*stresses)(0,5) = cauchystress(0,2);
  }
  break;
  case INPAR::STR::stress_none:
    break;
  default:
    dserror("requested stress type not available");
  }




  return;
} // DRT::ELEMENTS::NStet::calcstress


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::so_tet4_expol
(
    LINALG::Matrix<NUMGPT_NSTET,NUMSTR_NSTET>& stresses,
    LINALG::Matrix<NUMDOF_NSTET,1>& elevec1,
    LINALG::Matrix<NUMDOF_NSTET,1>& elevec2
)
{
  static LINALG::Matrix<NUMNOD_NSTET, NUMGPT_NSTET> expol;
  static bool isfilled;

  if (isfilled==false)
  {
    expol(0,0)=1.0;
    expol(1,0)=1.0;
    expol(2,0)=1.0;
    expol(3,0)=1.0;

    isfilled=true;
  }

  LINALG::Matrix<NUMNOD_NSTET,NUMSTR_NSTET> nodalstresses;
  nodalstresses.Multiply(expol,stresses);

  for (int i=0;i<NUMNOD_NSTET;++i)
  {
     elevec1(3*i)=nodalstresses(i,0);
     elevec1(3*i+1)=nodalstresses(i,1);
     elevec1(3*i+2)=nodalstresses(i,2);
   }
   for (int i=0;i<NUMNOD_NSTET;++i)
   {
     elevec2(3*i)=nodalstresses(i,3);
     elevec2(3*i+1)=nodalstresses(i,4);
     elevec2(3*i+2)=nodalstresses(i,5);
   }
  return;
}

/*----------------------------------------------------------------------*
 |  volumetric stabilization                                   gee 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::VolStabilization(
                    map<int,DRT::ELEMENTS::NStet*>&          adjele,
                    map<int,DRT::Node*>&                     adjnode,
                    vector<int>&                             adjlm,
                    vector<int>&                             adjlmowner,
                    vector<map<int,DRT::ELEMENTS::NStet*> >& nodaladjele,
                    vector<map<int,DRT::Node*> >&            nodaladjnode,
                    LINALG::SerialDenseMatrix&               stiff,
                    LINALG::SerialDenseVector&               force)
{

  //---------------------------------------------------------------------
  // build averaged deformation gradient and volume for each node
  const int                    ndofperpatch=(int)adjlm.size();
  if (isnan(ndofperpatch)) dserror("Found nan in ele %d",Id());
  vector<double>               J(NumNode());
  vector<LINALG::Matrix<3,3> > Fnode(NumNode());
  vector<double>               Vnode(NumNode(),0.0);
  const double                 Jele = F_.Determinant();
  for (int i=0; i<NumNode(); ++i)
  {
    Fnode[i].PutScalar(0.0);
    map<int,DRT::ELEMENTS::NStet*>::iterator fool;
    for (fool = nodaladjele[i].begin(); fool != nodaladjele[i].end(); ++fool)
    {
      const double Velequart = fool->second->Volume()/NUMNOD_NSTET;
      Vnode[i] += Velequart;
      Fnode[i].Update(Velequart,fool->second->F_,1.0);
    }
    Fnode[i].Scale(1.0/Vnode[i]);
    J[i] = Fnode[i].Determinant();
    if (J[i]<1.0e-8) dserror("Nodal detF almost zero or negative. detF %15.10e",J[i]);
  }

  //-----------------------------------------------------------------------
  // do positioning map node -> dofsinpatch
  map<int,int>  node_pos;
  map<int,DRT::Node*>::iterator pnode;
  int count=0;
  for (pnode=adjnode.begin(); pnode != adjnode.end(); ++pnode)
  {
    node_pos[pnode->first] = count;
    count++;
  }


  //----------------------------------------------------------------------
  // build nodal quantities strains,stresses, pressure, Grad(\delta u) and assemble them
  LINALG::Matrix<4,1>      p;                           // nodal pressure (linear within element)
  Epetra_SerialDenseMatrix G(NumNode(),ndofperpatch);   // G^I = ( F^{-T}:Grad(u) )^I
  for (int node=0; node<NumNode(); ++node)
  {
    LINALG::Matrix<6,6> cmatnode(true);    // nodal material tangent
    LINALG::Matrix<6,1> stressnode(true);  // nodal PKII stress
    LINALG::Matrix<6,1> gl;
    LINALG::Matrix<3,3> cg;
    cg.MultiplyTN(Fnode[node],Fnode[node]);
    gl(0) = 0.5 * (cg(0,0) - 1.0);
    gl(1) = 0.5 * (cg(1,1) - 1.0);
    gl(2) = 0.5 * (cg(2,2) - 1.0);
    gl(3) = cg(0,1);
    gl(4) = cg(1,2);
    gl(5) = cg(2,0);
    double density = 0.0;
    RCP<MAT::Material> mat = nodaladjele[node].begin()->second->Material();
    // maybe do this with volumetric strain only????
    ElementType().SelectMaterial(mat,stressnode,cmatnode,density,gl,Fnode[node],0);
    // compute p at node I: p = -1/3 * J^-1 * S:C
    p(node) = 0.0;
    p(node) += cg(0,0) * stressnode(0);
    p(node) += cg(0,1) * stressnode(3);
    p(node) += cg(0,2) * stressnode(5);
    p(node) += cg(1,0) * stressnode(3);
    p(node) += cg(1,1) * stressnode(1);
    p(node) += cg(1,2) * stressnode(4);
    p(node) += cg(2,0) * stressnode(5);
    p(node) += cg(2,1) * stressnode(4);
    p(node) += cg(2,2) * stressnode(2);
    p(node) *= ((-1.0/3.0)/J[node]);


    //----------------------------------------------------------------------
    // build nodal averaged Grad(\delta u) operator and assemble it
    map<int,DRT::ELEMENTS::NStet*>:: iterator fool;
    for (fool = nodaladjele[node].begin(); fool != nodaladjele[node].end(); ++fool)
    {
      DRT::ELEMENTS::NStet* ele = fool->second;
      const double VeletoVnodal = (ele->Volume()/4.0)/Vnode[node];
      // element discrete Grad(u) operator
      LINALG::Matrix<NUMNOD_NSTET,NUMDIM_NSTET>& nxyz = ele->nxyz_;
      LINALG::Matrix<9,12> Bl(true);
      for (int i=0; i<ele->NumNode(); ++i)
      {
        Bl(0,NODDOF_NSTET*i+0) = nxyz(i,0);
        Bl(1,NODDOF_NSTET*i+1) = nxyz(i,1);
        Bl(2,NODDOF_NSTET*i+2) = nxyz(i,2);

        Bl(3,NODDOF_NSTET*i+0) = nxyz(i,1);
        Bl(4,NODDOF_NSTET*i+1) = nxyz(i,2);
        Bl(5,NODDOF_NSTET*i+2) = nxyz(i,0);

        Bl(6,NODDOF_NSTET*i+0) = nxyz(i,2);
        Bl(7,NODDOF_NSTET*i+1) = nxyz(i,0);
        Bl(8,NODDOF_NSTET*i+2) = nxyz(i,1);
      }

      // element transposed inverse deformation gradient F^{-T}
      LINALG::Matrix<9,1> FinvT;
      {
        LINALG::Matrix<3,3> tmp;
        tmp.Invert(ele->F_);
        FinvT(0) = tmp(0,0);
        FinvT(1) = tmp(1,1);
        FinvT(2) = tmp(2,2);
        FinvT(3) = tmp(1,0);
        FinvT(4) = tmp(2,1);
        FinvT(5) = tmp(0,2);
        FinvT(6) = tmp(2,0);
        FinvT(7) = tmp(0,1);
        FinvT(8) = tmp(1,2);
      }
      // Gele = (Vele/4 / Vnodal) * F^-T : Grad(u)
      LINALG::Matrix<1,12> Gele;
      Gele.MultiplyTN(VeletoVnodal,FinvT,Bl);
      //----------------------------------------------------------------------
      // Assemble element contribution Gele into nodal patch
      for (int i=0; i<ele->NumNode(); ++i)
      {
        const int pos = node_pos[ele->Nodes()[i]->Id()];
        G(node,NODDOF_NSTET*pos+0) += Gele(i*NODDOF_NSTET+0);
        G(node,NODDOF_NSTET*pos+1) += Gele(i*NODDOF_NSTET+1);
        G(node,NODDOF_NSTET*pos+2) += Gele(i*NODDOF_NSTET+2);
      }
    } // for (fool = nodaladjele[node].begin(); fool != nodaladjele[node].end(); ++fool)
  } // for (int node=0; node<NumNode(); ++node)

  //----------------------------------------------------------------------
  // compute matrices M E for projection
  LINALG::Matrix<NUMNOD_NSTET,NUMNOD_NSTET> M(true);
  LINALG::Matrix<1,NUMNOD_NSTET>            E(true);
  {
    // 4 Gauss point integration
    const double alpha  = (5.0 + 3.0*sqrt(5.0))/20.0;
    const double beta   = (5.0 - sqrt(5.0))/20.0;
    const double weight = 0.25 * Volume();
    double xsi[4][4];
    xsi[0][0] = alpha;   xsi[0][1] = beta ;   xsi[0][2] = beta ;   xsi[0][3] = beta ;
    xsi[1][0] = beta ;   xsi[1][1] = alpha;   xsi[1][2] = beta ;   xsi[1][3] = beta ;
    xsi[2][0] = beta ;   xsi[2][1] = beta ;   xsi[2][2] = alpha;   xsi[2][3] = beta ;
    xsi[3][0] = beta ;   xsi[3][1] = beta ;   xsi[3][2] = beta ;   xsi[3][3] = alpha;
    for (int gp=0; gp<4; ++gp)
    {
      LINALG::Matrix<NUMNOD_NSTET,1> funct;
      ShapeFunction(funct,xsi[gp][0],xsi[gp][1],xsi[gp][2],xsi[gp][3]);
      for (int i=0; i<NUMNOD_NSTET; ++i)
      {
        E(i) += funct(i) * weight;
        for (int j=0; j<NUMNOD_NSTET; ++j)
          M(i,j) += funct(i) * funct(j) * weight;
      }
    }
  }


  //----------------------------------------------------------------------
  LINALG::Matrix<NUMNOD_NSTET,1>            phat;
  LINALG::Matrix<NUMNOD_NSTET,NUMNOD_NSTET> MmEDE;
  MmEDE.MultiplyTN(E,E);
  MmEDE.Update(1.0,M,-1.0/Volume());
  phat.Multiply(MmEDE,p);
  //printf("Id %d p %15.10e %15.10e %15.10e %15.10e \n",Id(),p(0),p(1),p(2),p(3));
  //printf("phat %15.10e %15.10e %15.10e %15.10e \n\n",phat(0),phat(1),phat(2),phat(3));
  //fflush(stdout);

  //----------------------------------------------------------------------
  // Multiply nodal Grad \delta u operator with projection and pressure and element volume
  Epetra_SerialDenseVector ephat(View,phat.A(),NUMNOD_NSTET);
  force.Multiply('T','N',-Jele*Volume(),G,ephat,0.0);

#if 0
  //---------------------------------------------------------------------
  // build projected pressure
  // pbar = D^-1 E p = 1/V * E p
  // (note that pbar is just the average of the nodal values ;-) )
  LINALG::Matrix<1,1> pbar;
  pbar.Multiply(E,p);
  pbar.Scale(1.0/Volume());
#endif

#if 0
  //----------------------------------------------------------------------
  // compute the element constant pressure for comparison
  {
    LINALG::Matrix<3,3>& F = F_;
    LINALG::Matrix<6,6> cmat(true);    // material tangent
    LINALG::Matrix<6,1> stress(true);  // PKII stress
    LINALG::Matrix<6,1> gl;
    LINALG::Matrix<3,3> cg;
    cg.MultiplyTN(F,F);
    gl(0) = 0.5 * (cg(0,0) - 1.0);
    gl(1) = 0.5 * (cg(1,1) - 1.0);
    gl(2) = 0.5 * (cg(2,2) - 1.0);
    gl(3) = cg(0,1);
    gl(4) = cg(1,2);
    gl(5) = cg(2,0);
    double density = 0.0;
    SelectMaterial(stress,cmat,density,gl,F,0);
    double pele = 0.0;
    pele += cg(0,0) * stress(0);
    pele += cg(0,1) * stress(3);
    pele += cg(0,2) * stress(5);
    pele += cg(1,0) * stress(3);
    pele += cg(1,1) * stress(1);
    pele += cg(1,2) * stress(4);
    pele += cg(2,0) * stress(5);
    pele += cg(2,1) * stress(4);
    pele += cg(2,2) * stress(2);

    pele *= ((-1.0/3.0)/F.Determinant());

    printf("ele %10d pbar %15.10e pele %15.10e\n",Id(),pbar(0),pele);
    fflush(stdout);
  }
#endif





  return;
} // DRT::ELEMENTS::NStet::VolStabilization







#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
