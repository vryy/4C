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

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
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
void DRT::ELEMENTS::Ptet::InitElement()
{
  LINALG::SerialDenseMatrix xrefe(NUMNOD_PTET,NUMDIM_PTET);
  LINALG::SerialDenseMatrix J(NUMNOD_PTET,NUMDIM_PTET+1);
  {
    // compute element volume
    for (int i=0; i<NUMNOD_PTET; ++i)
    {
      J(i,0) = 1.0;
      J(i,1) = xrefe(i,0) = Nodes()[i]->X()[0];
      J(i,2) = xrefe(i,1) = Nodes()[i]->X()[1];
      J(i,3) = xrefe(i,2) = Nodes()[i]->X()[2];
    }
    V_ = LINALG::DeterminantLU(J)/6.0;
    if (V_==0.0)     dserror("Element volume is zero");
    else if (V_<0.0) dserror("Element volume is negative");
  }
  // compute derivatives of shape functions w.r.t. to material coords.
  {
    // one gauss point at 0.25/0.25/0.25
    // gauss point weight is 1.0, so skip it
    const double gploc = 0.25;
    LINALG::SerialDenseVector funct(NUMNOD_PTET);
    ShapeFunction(funct,gploc,gploc,gploc,gploc);
    LINALG::SerialDenseMatrix deriv(NUMNOD_PTET,NUMCOORD_PTET);
    ShapeFunctionDerivatives(deriv);
    LINALG::SerialDenseMatrix tmp(NUMCOORD_PTET-1,NUMCOORD_PTET);
    tmp.Multiply('T','N',1.0,xrefe,deriv,0.0);
    for (int i=0; i<4; i++) J(0,i)=1;
    for (int row=0;row<3;row++)
      for (int col=0;col<4;col++)
        J(row+1,col)=tmp(row,col);
    Epetra_SerialDenseMatrix Iaug(NUMCOORD_PTET,NUMDIM_PTET);
    Iaug(1,0)=1;
    Iaug(2,1)=1;
    Iaug(3,2)=1;
    Epetra_SerialDenseMatrix partials(NUMCOORD_PTET,NUMDIM_PTET);
    Epetra_SerialDenseSolver solver;
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
    nxyz_.LightShape(NUMNOD_PTET,NUMDIM_PTET);
    nxyz_.Multiply('N','N',1.0,deriv,partials,0.0);
  }
  
  
  return;
} // DRT::ELEMENTS::Ptet::InitElement



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Ptet::Evaluate(ParameterList& params,
                                  DRT::Discretization&      discretization,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseMatrix& elemat1,
                                  Epetra_SerialDenseMatrix& elemat2,
                                  Epetra_SerialDenseVector& elevec1,
                                  Epetra_SerialDenseVector& elevec2,
                                  Epetra_SerialDenseVector& elevec3)
{
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
  else if (action=="calc_struct_stress")                  act = Ptet::calc_struct_stress;
  else if (action=="postprocess_stress")                  act = Ptet::postprocess_stress;
  else if (action=="calc_struct_eleload")                 act = Ptet::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")                 act = Ptet::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")            act = Ptet::calc_struct_update_istep;
  else if (action=="calc_struct_update_genalpha_imrlike") act = Ptet::calc_struct_update_genalpha_imrlike;
  else dserror("Unknown type of action for Ptet");

  // what should the element do
  switch(act) 
  {
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass: 
    {
      //cout << "calc_struct_nlnstiffmass\n";
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      ptetnlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,false);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: 
    {
      //cout << "calc_struct_nlnstiff\n";
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      ptetnlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,false);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress: 
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==null) dserror("Cannot get stress 'data'");
      if (straindata==null) dserror("Cannot get strain 'data'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      Epetra_SerialDenseMatrix stress(NUMGPT_PTET,NUMSTR_PTET);
      Epetra_SerialDenseMatrix strain(NUMGPT_PTET,NUMSTR_PTET);
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      if (iostrain!="euler_almansi") 
        ptetnlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,cauchy);
      else
       dserror("requested option euler_almansi not yet implemented for ptet");
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
      const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      RCP<Epetra_SerialDenseMatrix> gpstress = (*gpstressmap)[gid];

      if (stresstype=="ndxyz") 
      {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_PTET,NUMSTR_PTET);
        so_tet4_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_PTET);

        for (int i=0;i<NUMNOD_PTET;++i)
        {
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_PTET;++i)
        {
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_PTET;++i)
        {
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
      }
      else if (stresstype=="cxyz") 
      {
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_PTET; ++i) 
          {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_PTET; ++j) 
            {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_PTET * (*gpstress)(j,i);
            }
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") 
      {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_PTET,NUMSTR_PTET);
        so_tet4_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_PTET);

        for (int i=0;i<NUMNOD_PTET;++i)
        {
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_PTET;++i)
        {
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_PTET;++i)
        {
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_PTET; ++i) 
          {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_PTET; ++j) 
            {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_PTET * (*gpstress)(j,i);
            }
          }
        }
      }
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

    case calc_struct_update_genalpha_imrlike: 
    {
      ;// there is nothing to do here at the moment
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
      dserror("Case 'calc_struct_internalforce' not implemented");
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
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      Epetra_SerialDenseMatrix* elestress,      // stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      // strains at GP
      const bool                cauchy)         // stress output options
{
  //--------------------------------------------------- geometry update
  if (!FisNew_) DeformationGradient(disp);
  Epetra_SerialDenseMatrix& defgrd = F_;
  // reset the bool indicating that the stored deformation gradient is 'fresh'
  FisNew_ = false;
  
  //--------------------------- Right Cauchy-Green tensor C = = F^T * F
  LINALG::SerialDenseMatrix cauchygreen(NUMDIM_PTET,NUMDIM_PTET);
  cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);
  
  // --Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  LINALG::SerialDenseVector glstrain(NUMSTR_PTET);
  glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
  glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
  glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
  glstrain(3) = cauchygreen(0,1);
  glstrain(4) = cauchygreen(1,2);
  glstrain(5) = cauchygreen(2,0);
  
  //----------------------- output strains if desired for visualization
  if (elestrain)
  {
    for (int i = 0; i < 3; ++i) (*elestrain)(0,i) = glstrain(i);
    for (int i = 3; i < 6; ++i) (*elestrain)(0,i) = 0.5 * glstrain(i);
  }
  
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
  LINALG::SerialDenseMatrix bop(NUMSTR_PTET,NUMDOF_PTET);
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
#if 1
  const double third = 1./3.;
  Epetra_SerialDenseMatrix Idev(NUMSTR_PTET,NUMSTR_PTET);
  Idev(0,0) =  2.0;  Idev(0,1) = -1.0;  Idev(0,2) = -1.0;
  Idev(1,0) = -1.0;  Idev(1,1) =  2.0;  Idev(1,2) = -1.0;
  Idev(2,0) = -1.0;  Idev(2,1) = -1.0;  Idev(2,2) =  2.0;
  Idev(3,3) = 1.0;
  Idev(4,4) = 1.0;
  Idev(5,5) = 1.0;
  Idev.Scale(third);
  
  // keep only deviatoric part of strains
  LINALG::SerialDenseVector gltmp(glstrain);
  Idev.Multiply(false,gltmp,glstrain);
  

  Epetra_SerialDenseMatrix cmat(NUMSTR_PTET,NUMSTR_PTET);
  Epetra_SerialDenseVector stress(NUMSTR_PTET);
  double density = -999.99;
  SelectMaterial(stress,cmat,density,glstrain,defgrd,0);
  stress.Scale(ALPHA_PTET);
  cmat.Scale(ALPHA_PTET);

  // keep only the deviatoric part of stresses
  Epetra_SerialDenseVector stresstmp(stress);
  Idev.Multiply(false,stresstmp,stress);

  LINALG::SerialDenseMatrix tmp(NUMSTR_PTET,NUMSTR_PTET);
  tmp.Multiply('N','N',1.0,cmat,Idev,0.0);
  cmat.Multiply('N','N',1.0,Idev,tmp,0.0);

#else
  Epetra_SerialDenseMatrix cmat(NUMSTR_PTET,NUMSTR_PTET);
  Epetra_SerialDenseVector stress(NUMSTR_PTET);
  double density = -999.99;
  SelectMaterial(stress,cmat,density,glstrain,defgrd,0);
  stress.Scale(ALPHA_PTET);
  cmat.Scale(ALPHA_PTET);
#endif
   
  //------------------------- return stresses desired for visualization
  if (elestress)
  {
    // return 2nd Piola-Kirchhoff stresses
    if (!cauchy)
      for (int i = 0; i < NUMSTR_PTET; ++i) (*elestress)(0,i) = (1./ALPHA_PTET)*stress(i);
    // Cauchy stresses
    else
    {
        const double detF = defgrd(0,0)*defgrd(1,1)*defgrd(2,2) +
                            defgrd(0,1)*defgrd(1,2)*defgrd(2,0) +
                            defgrd(0,2)*defgrd(1,0)*defgrd(2,1) -
                            defgrd(0,2)*defgrd(1,1)*defgrd(2,0) -
                            defgrd(0,0)*defgrd(1,2)*defgrd(2,1) -
                            defgrd(0,1)*defgrd(1,0)*defgrd(2,2);
        LINALG::SerialDenseMatrix pkstress(NUMDIM_PTET,NUMDIM_PTET);
        pkstress(0,0) = stress(0);
        pkstress(0,1) = stress(3);
        pkstress(0,2) = stress(5);
        pkstress(1,0) = pkstress(0,1);
        pkstress(1,1) = stress(1);
        pkstress(1,2) = stress(4);
        pkstress(2,0) = pkstress(0,2);
        pkstress(2,1) = pkstress(1,2);
        pkstress(2,2) = stress(2);
        pkstress.Scale(1./ALPHA_PTET);

        LINALG::SerialDenseMatrix temp(NUMDIM_PTET,NUMDIM_PTET);
        LINALG::SerialDenseMatrix cauchystress(NUMDIM_PTET,NUMDIM_PTET);
        temp.Multiply('N','N',1.0/detF,defgrd,pkstress,0.);
        cauchystress.Multiply('N','T',1.0,temp,defgrd,0.);
        (*elestress)(0,0) = cauchystress(0,0);
        (*elestress)(0,1) = cauchystress(1,1);
        (*elestress)(0,2) = cauchystress(2,2);
        (*elestress)(0,3) = cauchystress(0,1);
        (*elestress)(0,4) = cauchystress(1,2);
        (*elestress)(0,5) = cauchystress(0,2);
      
    }
  }
  
  if (force && stiffmatrix)
  {
    // integrate internal force vector f = f + (B^T . sigma) * V_ 
    (*force).Multiply('T','N',V_,bop,stress,1.0);

    // integrate elastic stiffness matrix
    // keu = keu + (B^T . C . B) * V_
    LINALG::SerialDenseMatrix cb(NUMSTR_PTET,NUMDOF_PTET);
    cb.Multiply('N','N',1.0,cmat,bop,0.0); 
    stiffmatrix->Multiply('T','N',V_,bop,cb,1.0);
    // integrate `geometric' stiffness matrix and add to keu
    double sBL[3];
    const double V = Volume();
    for (int i=0; i<NUMNOD_PTET; ++i)
    {
      sBL[0] = V*stress(0) * nxyz_(i,0) + V*stress(3) * nxyz_(i,1) + V*stress(5) * nxyz_(i,2);
      sBL[1] = V*stress(3) * nxyz_(i,0) + V*stress(1) * nxyz_(i,1) + V*stress(4) * nxyz_(i,2);
      sBL[2] = V*stress(5) * nxyz_(i,0) + V*stress(4) * nxyz_(i,1) + V*stress(2) * nxyz_(i,2);
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
      LINALG::SerialDenseVector funct(NUMNOD_PTET);
      ShapeFunction(funct,xsi[gp][0],xsi[gp][1],xsi[gp][2],xsi[gp][3]);
      for (int i=0; i<NUMNOD_PTET; ++i)
        for (int j=0; j<NUMNOD_PTET; ++j)
        {
          const double fac = funct(i) * funct(j) * density * V * weight;
          (*massmatrix)(NUMDIM_PTET*i+0,NUMDIM_PTET*j+0) += fac;
          (*massmatrix)(NUMDIM_PTET*i+1,NUMDIM_PTET*j+1) += fac;
          (*massmatrix)(NUMDIM_PTET*i+2,NUMDIM_PTET*j+2) += fac;
        }
    }
  } // if (massmatrix)

  return;
} // DRT::ELEMENTS::Ptet::ptetnlnstiffmass




/*----------------------------------------------------------------------*
 |  compute deformation gradient (protected)                   gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ptet::DeformationGradient(vector<double>& disp)
{
  LINALG::SerialDenseMatrix xdisp(NUMNOD_PTET,NUMDIM_PTET);
  F_.Shape(NUMDIM_PTET,NUMDIM_PTET);
  for (int i=0; i<NUMNOD_PTET; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_PTET+0];
    xdisp(i,1) = disp[i*NODDOF_PTET+1];
    xdisp(i,2) = disp[i*NODDOF_PTET+2];
  }
  F_.Multiply('T','N',1.0,xdisp,nxyz_,0.0);
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
