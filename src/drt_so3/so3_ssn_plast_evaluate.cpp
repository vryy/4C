/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast_evaluate.cpp
\brief

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>

*/


/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast.H"

#include "so3_ssn_plast_fwd.hpp"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/plasticelasthyper.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/material_service.H"
#include "../drt_inpar/inpar_structure.H"

// headers of supported hyperelastic-materials

/*----------------------------------------------------------------------*
 | evaluate the element (public)                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::Evaluate(
  Teuchos::ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra
  )
{
  // start with "none"
  ActionType act = none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                        act =  calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                        act =  calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")                   act =  calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")                    act =  calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")                    act =  calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")                   act =  calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                          act =  calc_struct_stress;
  else if (action=="calc_struct_update_istep")                    act =  calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")                     act =  calc_struct_reset_istep;
  else if (action=="postprocess_stress")                          act =  postprocess_stress;
  else
    dserror("Unknown type of action for So3_Plast: %s",action.c_str());

  // what should the element do
  switch(act)
  {
  //============================================================================
  // linear stiffness
  case calc_struct_linstiff:
  case calc_struct_linstiffmass:
  {
    dserror("linear kinematics version out dated");
    break;
  }

  //============================================================================
  // nonlinear stiffness
  case calc_struct_internalforce:
  {
    // internal force vector
    LINALG::Matrix<numdofperelement_,1> elevec1(elevec1_epetra.A(),true);
    // elemat1+2, elevec2+3 are not used anyway

    // need current displacement and residual/incremental displacements
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0,"displacement");
    Teuchos::RCP<const Epetra_Vector> res
      = discretization.GetState(0,"residual displacement");
    if ( (disp == Teuchos::null) or (res==Teuchos::null) )
      dserror("Cannot get state vectors 'displacement' and/or residual");
    std::vector<double> mydisp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
    std::vector<double> myres(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);
    // create a dummy element matrix to apply linearised EAS-stuff onto
    LINALG::Matrix<numdofperelement_,numdofperelement_> myemat(true);

    // default: geometrically non-linear analysis with Total Lagrangean approach
    nln_stiffmass(la[0].lm_,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
        INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());


    break;
  }

  //============================================================================
  // nonlinear stiffness
  case calc_struct_nlnstiff:
  {
    // stiffness
    LINALG::Matrix<numdofperelement_,numdofperelement_> elemat1(elemat1_epetra.A(),true);
    LINALG::Matrix<numdofperelement_,numdofperelement_>* matptr = NULL;
    if (elemat1.IsInitialized()) matptr = &elemat1;
    LINALG::Matrix<numdofperelement_,1> elevec1(elevec1_epetra.A(),true);

    // need current displacement and residual forces
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
    if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
    std::vector<double> mydisp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
    std::vector<double> myres(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);

    // default: geometrically non-linear analysis with Total Lagrangean approach
    nln_stiffmass(la[0].lm_,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
        INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());

    break;
  }  // calc_struct_nlnstiff

  //============================================================================
  // (non)linear stiffness, mass matrix and internal force vector
  case calc_struct_nlnstiffmass:
  case calc_struct_nlnstifflmass:
  {
    // need current displacement and residual forces
     Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
     Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
     if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
     std::vector<double> mydisp(la[0].lm_.size());
     DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
     std::vector<double> myres(la[0].lm_.size());
     DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);
     // stiffness
     LINALG::Matrix<numdofperelement_,numdofperelement_> elemat1(elemat1_epetra.A(),true);
     // mass
     LINALG::Matrix<numdofperelement_,numdofperelement_> elemat2(elemat2_epetra.A(),true);
     // internal force
     LINALG::Matrix<numdofperelement_,1> elevec1(elevec1_epetra.A(),true);

     // default: geometrically non-linear analysis with Total Lagrangean approach
     nln_stiffmass(la[0].lm_,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
         INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());

     if(act==calc_struct_nlnstifflmass)
       // lump mass matrix
       // we assume #elemat2 is a square matrix
       for (int c=0; c<elemat2_epetra.N(); ++c)  // parse columns
       {
         double d = 0.0;
         for (int r=0; r<elemat2_epetra.M(); ++r)  // parse rows
         {
           d += elemat2(r,c);  // accumulate row entries
           elemat2(r,c) = 0.0;
         }
         elemat2(c,c) = d;  // apply sum of row entries on diagonal
       }
    break;
  }  // calc_struct_nlnstiff(l)mass

  case calc_struct_stress:
  {
    // elemat1+2,elevec1-3 are not used anyway

    // nothing to do for ghost elements
    if (discretization.Comm().MyPID() == Owner())
    {
      Teuchos::RCP<const Epetra_Vector> disp
        = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> res
        = discretization.GetState(0,"residual displacement");
      if ( (disp == Teuchos::null) or (res == Teuchos::null) )
        dserror("Cannot get state vectors 'displacement'");
      Teuchos::RCP<std::vector<char> > stressdata = params.get<Teuchos::RCP<std::vector<char> > >("stress", Teuchos::null);
      Teuchos::RCP<std::vector<char> > straindata = params.get<Teuchos::RCP<std::vector<char> > >("strain",Teuchos:: null);
      if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
      if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");

      std::vector<double> mydisp((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
      std::vector<double> myres((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);

      LINALG::Matrix<numgpt_post,numstr_> stress;
      LINALG::Matrix<numgpt_post,numstr_> strain;
      INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(la[0].lm_,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,
          iostress,iostrain,discretization.Comm().MyPID());
      {
        DRT::PackBuffer data;
        this->AddtoPack(data, stress);
        data.StartPacking();
        this->AddtoPack(data, stress);
        std::copy(data().begin(),data().end(),std::back_inserter(*stressdata));
      }
      {
        DRT::PackBuffer data;
        this->AddtoPack(data, strain);
        data.StartPacking();
        this->AddtoPack(data, strain);
        std::copy(data().begin(),data().end(),std::back_inserter(*straindata));
      }
    }

    break;
  }


  //============================================================================
  // required for predictor TangDis --> can be helpful in compressible case!
  case calc_struct_reset_istep:
  {
    for (int i=0; i<numgpt_; i++)
    {
      (*KbbInv_)[i].Scale(0.);
      (*Kbd_)[i].Scale(0.);
      (*fbeta_)[i].Scale(0.);
    }
    break;
  }

  //============================================================================
  case calc_struct_update_istep:
  {
    // update plastic deformation
    // default: geometrically non-linear analysis with Total Lagrangean approach
    UpdatePlasticDeformation_nln(plspintype_);
    break;
  }  // calc_struct_update_istep

  // note that in the following, quantities are always referred to as
  // "stresses" etc. although they might also apply to strains
  // (depending on what this routine is called for from the post filter)
  case postprocess_stress:
  {
    const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
      params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
    if (gpstressmap==Teuchos::null)
      dserror("no gp stress/strain map available for postprocessing");
    std::string stresstype = params.get<std::string>("stresstype","ndxyz");
    int gid = Id();
    LINALG::Matrix<numgpt_post,numstr_> gpstress(((*gpstressmap)[gid])->A(),true);

    Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
    if (poststress==Teuchos::null)
      dserror("No element stress/strain vector available");

   if (stresstype=="ndxyz")
    {
      if (distype==DRT::Element::hex8)
        soh8_expol(gpstress, *poststress);
      else
        dserror("only element centered stress for so3_plast");
    }
    else if (stresstype=="cxyz")
    {
      const Epetra_BlockMap& elemap = poststress->Map();
      int lid = elemap.LID(Id());
      if (lid!=-1)
      {
        for (int i = 0; i < numstr_; ++i)
        {
          double& s = (*((*poststress)(i)))[lid]; // resolve pointer for faster access
          s = 0.;
          for (int j = 0; j < numgpt_post; ++j)
          {
            s += gpstress(j,i);
          }
          s *= 1.0/numgpt_post;
        }
      }
    }
    else
    {
      dserror("unknown type of stress/strain output on element level");
    }
    break;
  }

// calc_struct_stifftemp

  //============================================================================
  default:
    dserror("Unknown type of action for So3_Plast");
    break;
  } // action

  return 0;
}  // Evaluate()





/*----------------------------------------------------------------------*
 | calculate the nonlinear B-operator                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::CalculateBop(
  LINALG::Matrix<numstr_,numdofperelement_>* bop,
  LINALG::Matrix<nsd_,nsd_>* defgrd,
  LINALG::Matrix<nsd_,nen_>* N_XYZ
  )
{
  // lump mass matrix
  if (bop != NULL)
  {
    /* non-linear B-operator (may so be called, meaning of B-operator is not so
    **  sharp in the non-linear realm) *
    **   B = F^{i,T} . B_L *
    ** with linear B-operator B_L =  N_XYZ (6x24) = (3x8)
    **
    **   B    =   F  . N_XYZ
    ** (6x24)   (3x3) (3x8)
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
    for (int i=0; i<nen_; ++i)
    {
      (*bop)(0,numdofpernode_*i+0) = (*defgrd)(0,0)*(*N_XYZ)(0,i);
      (*bop)(0,numdofpernode_*i+1) = (*defgrd)(1,0)*(*N_XYZ)(0,i);
      (*bop)(0,numdofpernode_*i+2) = (*defgrd)(2,0)*(*N_XYZ)(0,i);
      (*bop)(1,numdofpernode_*i+0) = (*defgrd)(0,1)*(*N_XYZ)(1,i);
      (*bop)(1,numdofpernode_*i+1) = (*defgrd)(1,1)*(*N_XYZ)(1,i);
      (*bop)(1,numdofpernode_*i+2) = (*defgrd)(2,1)*(*N_XYZ)(1,i);
      (*bop)(2,numdofpernode_*i+0) = (*defgrd)(0,2)*(*N_XYZ)(2,i);
      (*bop)(2,numdofpernode_*i+1) = (*defgrd)(1,2)*(*N_XYZ)(2,i);
      (*bop)(2,numdofpernode_*i+2) = (*defgrd)(2,2)*(*N_XYZ)(2,i);
      /* ~~~ */
      (*bop)(3,numdofpernode_*i+0) = (*defgrd)(0,0)*(*N_XYZ)(1,i) + (*defgrd)(0,1)*(*N_XYZ)(0,i);
      (*bop)(3,numdofpernode_*i+1) = (*defgrd)(1,0)*(*N_XYZ)(1,i) + (*defgrd)(1,1)*(*N_XYZ)(0,i);
      (*bop)(3,numdofpernode_*i+2) = (*defgrd)(2,0)*(*N_XYZ)(1,i) + (*defgrd)(2,1)*(*N_XYZ)(0,i);
      (*bop)(4,numdofpernode_*i+0) = (*defgrd)(0,1)*(*N_XYZ)(2,i) + (*defgrd)(0,2)*(*N_XYZ)(1,i);
      (*bop)(4,numdofpernode_*i+1) = (*defgrd)(1,1)*(*N_XYZ)(2,i) + (*defgrd)(1,2)*(*N_XYZ)(1,i);
      (*bop)(4,numdofpernode_*i+2) = (*defgrd)(2,1)*(*N_XYZ)(2,i) + (*defgrd)(2,2)*(*N_XYZ)(1,i);
      (*bop)(5,numdofpernode_*i+0) = (*defgrd)(0,2)*(*N_XYZ)(0,i) + (*defgrd)(0,0)*(*N_XYZ)(2,i);
      (*bop)(5,numdofpernode_*i+1) = (*defgrd)(1,2)*(*N_XYZ)(0,i) + (*defgrd)(1,0)*(*N_XYZ)(2,i);
      (*bop)(5,numdofpernode_*i+2) = (*defgrd)(2,2)*(*N_XYZ)(0,i) + (*defgrd)(2,0)*(*N_XYZ)(2,i);
    }
  }
}  // CalculateBop()


/*----------------------------------------------------------------------*
 | calculate the linear B-operator                          seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::CalculateBoplin(
  LINALG::Matrix<numstr_,numdofperelement_>* boplin,
  LINALG::Matrix<nsd_,nen_>* N_XYZ
  )
{
  // lump mass matrix
  if (boplin != NULL)
  {
    // linear B-operator B = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    for (int i=0; i<nen_; ++i)
    {
      (*boplin)(0,numdofpernode_*i+0) = (*N_XYZ)(0,i);
      (*boplin)(0,numdofpernode_*i+1) = 0.0;
      (*boplin)(0,numdofpernode_*i+2) = 0.0;
      (*boplin)(1,numdofpernode_*i+0) = 0.0;
      (*boplin)(1,numdofpernode_*i+1) = (*N_XYZ)(1,i);
      (*boplin)(1,numdofpernode_*i+2) = 0.0;
      (*boplin)(2,numdofpernode_*i+0) = 0.0;
      (*boplin)(2,numdofpernode_*i+1) = 0.0;
      (*boplin)(2,numdofpernode_*i+2) = (*N_XYZ)(2,i);
      /* ~~~ */
      (*boplin)(3,numdofpernode_*i+0) = (*N_XYZ)(1,i);
      (*boplin)(3,numdofpernode_*i+1) = (*N_XYZ)(0,i);
      (*boplin)(3,numdofpernode_*i+2) = 0.0;
      (*boplin)(4,numdofpernode_*i+0) = 0.0;
      (*boplin)(4,numdofpernode_*i+1) = (*N_XYZ)(2,i);
      (*boplin)(4,numdofpernode_*i+2) = (*N_XYZ)(1,i);
      (*boplin)(5,numdofpernode_*i+0) = (*N_XYZ)(2,i);
      (*boplin)(5,numdofpernode_*i+1) = 0.0;
      (*boplin)(5,numdofpernode_*i+2) = (*N_XYZ)(0,i);
    }
  }
}  // CalculateBoplin()

/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)  seitz 04/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::EvaluateNeumann(Teuchos::ParameterList&   params,
                                            DRT::Discretization&      discretization,
                                            DRT::Condition&           condition,
                                            std::vector<int>&         lm,
                                            Epetra_SerialDenseVector& elevec1,
                                            Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
  // **

  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int> >("funct");
  LINALG::Matrix<nsd_,1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim=0; dim<nsd_; dim++)
      if ((*funct)[dim] > 0)
        havefunct = havefunct or true;

  // update element geometry
  LINALG::Matrix<nen_,nsd_> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<nen_; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<numgpt_; ++gp) {

    // shape functions (shapefunct) and their first derivatives (deriv)
    LINALG::Matrix<nen_,1> shapefunct;
    DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
    LINALG::Matrix<nsd_,nen_> deriv;
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    // compute the Jacobian matrix
    LINALG::Matrix<nsd_,nsd_> jac;
    jac.Multiply(deriv,xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct) {
      for (int dim=0; dim<nsd_; dim++) {
        xrefegp(dim) = 0.0;
        for (int nodid=0; nodid<nen_; ++nodid)
          xrefegp(dim) += shapefunct(nodid) * xrefe(nodid,dim);
      }
    }

    // integration factor
    const double fac = wgt_[gp] * curvefac * detJ;
    // distribute/add over element load vector
    for(int dim=0; dim<nsd_; dim++) {
      // function evaluation
      const int functnum = (funct) ? (*funct)[dim] : -1;
      const double functfac
        = (functnum>0)
        ? DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,xrefegp.A(),time,NULL)
        : 1.0;
      const double dim_fac = (*onoff)[dim] * (*val)[dim] * fac * functfac;
      for (int nodid=0; nodid<nen_; ++nodid) {
        elevec1[nodid*nsd_+dim] += shapefunct(nodid) * dim_fac;
      }
    }

  }/* ==================================================== end of Loop over GP */

  return 0;
}



/*----------------------------------------------------------------------*
 | initialise Jacobian                                      seitz 07/13 |
 | is called once in Initialize() in so3_ssn_plast_eletypes.cpp         |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::InitJacobianMapping()
{
  // get the material coordinates
  LINALG::Matrix<nen_,nsd_> xrefe;
  for (int i=0; i<nen_; ++i)
  {
    Node** nodes = Nodes();
    if (!nodes) dserror("Nodes() returned null pointer");
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);
  xsi_.resize(numgpt_);

  // initialise the derivatives of the shape functions
  LINALG::Matrix<nsd_,nen_> deriv;

  // coordinates of the current integration point (xsi_)
  for (int gp=0; gp<numgpt_; ++gp)
  {
    // first derivatives of shape functions (deriv)
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    // compute Jacobian matrix and determinant
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
    // derivatives of coordinates w.r.t material coordinates xjm_ = dx/ds
    invJ_[gp].Multiply(deriv,xrefe);
    // xij_ = ds/dx
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] < 1.0E-16)
      dserror("ZERO OR NEGATIVE JACOBIAN DETERMINANT: %f",detJ_[gp]);
  }  // end gp loop

  return;
}  // InitJacobianMapping()}


/*----------------------------------------------------------------------*
 | internal force, stiffness and mass for f-bar elements    seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::nln_stiffmass(
    std::vector<int>&              lm,             // location matrix
    std::vector<double>&           disp,           // current displacements
    std::vector<double>&           residual,       // current residual displ
    LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix, // element stiffness matrix
    LINALG::Matrix<numdofperelement_,numdofperelement_>* massmatrix,  // element mass matrix
    LINALG::Matrix<numdofperelement_,1>* force,                 // element internal force vector
    LINALG::Matrix<numgpt_post,numstr_>* elestress,   // stresses at GP
    LINALG::Matrix<numgpt_post,numstr_>* elestrain,   // strains at GP
    Teuchos::ParameterList&        params,         // algorithmic parameters e.g. time
    const INPAR::STR::StressType   iostress,  // stress output option
    const INPAR::STR::StrainType   iostrain,  // strain output option
    const int MyPID  // processor id
    )
{
  // update element geometry
  LINALG::Matrix<nen_,3> xrefe;  // X, material coord. of element
  LINALG::Matrix<nen_,3> xcurr;  // x, current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<nen_; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*numdofpernode_+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*numdofpernode_+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*numdofpernode_+2];
  }

  // we need the (residual) displacement -- current increment of displacement
  LINALG::Matrix<numdofperelement_,1> res_d;
  for (int i = 0; i<numdofperelement_; ++i)
  {
    res_d(i) = residual[i];
  }

  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  // by N_XYZ = J^-1 * N_rst
  LINALG::Matrix<nsd_,nen_> N_XYZ;
  // build deformation gradient wrt to material configuration
  LINALG::Matrix<nsd_,nsd_> defgrd(false);
  // shape functions and their first derivatives
  LINALG::Matrix<nen_,1> shapefunct;
  LINALG::Matrix<nsd_,nen_> deriv;

  // ---------------------- deformation gradient at centroid of element
  double detF_0 = -1.0;
  LINALG::Matrix<nsd_,nsd_> defgrd_0(false);
  LINALG::Matrix<nsd_,nsd_> invdefgrd_0(false);
  LINALG::Matrix<3,nen_> N_XYZ_0(false);
  if(fbar_)
  {
    //element coordinate derivatives at centroid
    LINALG::Matrix<3,nen_> N_rst_0(false);
    DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0.0, 0.0, 0.0, DRT::Element::hex8);

    //inverse jacobian matrix at centroid
    LINALG::Matrix<3,3> invJ_0(false);
    invJ_0.Multiply(N_rst_0,xrefe);
    invJ_0.Invert();
    //material derivatives at centroid
    N_XYZ_0.Multiply(invJ_0,N_rst_0);

    //deformation gradient and its determinant at centroid
    defgrd_0.MultiplyTT(xcurr,N_XYZ_0);
    invdefgrd_0.Invert(defgrd_0);
    detF_0 = defgrd_0.Determinant();
  }

  // 3x3 LINALG matrix for temporary stuff
  LINALG::Matrix<nsd_,nsd_> tmp1(false);
  LINALG::Matrix<nsd_,nsd_> tmp2(false);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // get references from parameter list
  double& lp_inc = params.get<double>("Lp_increment_square");
  double& lp_res = params.get<double>("Lp_residual_square");
  double& eas_inc= params.get<double>("EAS_increment_square");
  double& eas_res= params.get<double>("EAS_residual_square");
  int& num_active_gp = params.get<int>("number_active_plastic_gp");
  INPAR::STR::PredEnum pred = INPAR::STR::pred_vague;
  if (params.isParameter("predict_type"))
    pred = params.get<INPAR::STR::PredEnum>("predict_type");
  bool tang_pred = false;
  if (params.isParameter("eval_tang_pred"))
    tang_pred = params.get<bool>("eval_tang_pred");

  /* evaluation of EAS variables (which are constant for the following):
  ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
  ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
  ** -> T0^{-T}
  */
  std::vector<Epetra_SerialDenseMatrix>* M_GP = NULL;   // EAS matrix M at all GPs
  LINALG::SerialDenseMatrix M;      // EAS matrix M at current GP
  double detJ0;                     // detJ(origin)
  Epetra_SerialDenseVector delta_alpha_eas(neas_); // incremental change of the alphas in this iteration
  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  LINALG::Matrix<numstr_,numstr_> T0invT;  // trafo matrix

  if (eastype_!=soh8p_easnone)
  {
    EasSetup(&M_GP,detJ0,T0invT,xrefe);

    /* end of EAS Update ******************/
    // add Kda . res_d to feas
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    delta_alpha_eas.Size(neas_);
    if (stiffmatrix!=NULL)
    {
      // constant predictor
      if (pred == INPAR::STR::pred_constdis)
        alpha_eas_=alpha_eas_last_timestep_;

      // tangential predictor
      else if (pred == INPAR::STR::pred_tangdis)
        for (int i=0; i<neas_; i++)
          (*alpha_eas_)(i) = 0.
          + (*alpha_eas_last_timestep_)(i)
          + (*alpha_eas_delta_over_last_timestep_)(i)
          ;

      // do usual recovery
      else if (pred == INPAR::STR::pred_vague)
      {
        switch(eastype_)
        {
        case soh8p_easmild:
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,numdofperelement_,1>(1.0, feas_->A(), 1.0, Kad_->A(), res_d.A());
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,soh8p_easmild,1>(0.0,delta_alpha_eas,-1.0,*KaaInv_,*feas_);
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,soh8p_easmild,1>(1.0,*alpha_eas_,-1.0,*KaaInv_,*feas_);
          break;
        case soh8p_easfull:
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,numdofperelement_,1>(1.0, feas_->A(), 1.0, Kad_->A(), res_d.A());
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,soh8p_easfull,1>(0.0,delta_alpha_eas,-1.0,*KaaInv_,*feas_);
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,soh8p_easfull,1>(1.0,*alpha_eas_,-1.0,*KaaInv_,*feas_);
          break;
        case soh8p_easnone: break;
        default: dserror("Don't know what to do with EAS type %d", eastype_); break;
        }
      }
    }
    eas_inc += pow(delta_alpha_eas.Norm2(),2.);

    // reset EAS matrices
    KaaInv_->Shape(neas_,neas_);
    Kad_->Shape(neas_,numdofperelement_);
    feas_->Size(neas_);
    /* end of EAS Update ******************/
  }

  // EAS matrix block
  Epetra_SerialDenseMatrix Kda(numdofperelement_,neas_);
  // temporary Epetra matrix for this and that
  Epetra_SerialDenseMatrix tmp;
  // RHS of EAS equation without condensed plasticity
  // (to calculate the residual of this equation)
  Epetra_SerialDenseVector feas_uncondensed(neas_);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<numgpt_; ++gp)
  {
    // shape functions (shapefunct) and their first derivatives (deriv)
    DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);

    /* get the inverse of the Jacobian matrix which looks like:
     **            [ x_,r  y_,r  z_,r ]^-1
     **     J^-1 = [ x_,s  y_,s  z_,s ]
     **            [ x_,t  y_,t  z_,t ]
     */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // (material) deformation gradient
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ);

    // calcualte total rcg
    // no F-bar modification here
    LINALG::Matrix<3,3> total_cauchygreen(false);
    total_cauchygreen.MultiplyTN(defgrd,defgrd);
    // total Cauchy green in voigt notation
    LINALG::Matrix<numstr_,1> RCG;
    for (int i=0; i<3; i++) RCG(i)=total_cauchygreen(i,i);
    RCG(3)=total_cauchygreen(0,1)*2.;
    RCG(4)=total_cauchygreen(1,2)*2.;
    RCG(5)=total_cauchygreen(0,2)*2.;

    // deformation gradient consistent with (potentially EAS-modified) GL strains
    // without eas this is equal to the regular defgrd.
    LINALG::Matrix<3,3> defgrd_mod(defgrd);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8p_easnone)
    {
      // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
      LINALG::Matrix<numstr_,1> total_glstrain(false);
      total_glstrain(0) = 0.5 * (total_cauchygreen(0,0) - 1.0);
      total_glstrain(1) = 0.5 * (total_cauchygreen(1,1) - 1.0);
      total_glstrain(2) = 0.5 * (total_cauchygreen(2,2) - 1.0);
      total_glstrain(3) = total_cauchygreen(0,1);
      total_glstrain(4) = total_cauchygreen(1,2);
      total_glstrain(5) = total_cauchygreen(2,0);
      M.LightShape(numstr_,neas_);
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      //Epetra_SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      switch(eastype_)
      {
      case soh8p_easfull:
        LINALG::DENSEFUNCTIONS::multiply<double,numstr_,numstr_,soh8p_easfull>(M.A(), detJ0/detJ, T0invT.A(), (M_GP->at(gp)).A());
        LINALG::DENSEFUNCTIONS::multiply<double,numstr_,soh8p_easfull,1>(1.0,total_glstrain.A(),1.0,M.A(),alpha_eas_->A());
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double,numstr_,numstr_,soh8p_easmild>(M.A(), detJ0/detJ, T0invT.A(), (M_GP->at(gp)).A());
        LINALG::DENSEFUNCTIONS::multiply<double,numstr_,soh8p_easmild,1>(1.0,total_glstrain.A(),1.0,M.A(),alpha_eas_->A());
        break;
      case soh8p_easnone: break;
      default: dserror("Don't know what to do with EAS type %d", eastype_); break;
      }

      // calculate deformation gradient consistent with modified GL strain tensor
      CalcConsistentDefgrd(defgrd,total_glstrain,defgrd_mod);
    } // ------------------------------------------------------------------ EAS

    // calculate nonlinear B-operator
    LINALG::Matrix<numstr_,numdofperelement_> bop(false);
    CalculateBop(&bop,&defgrd,&N_XYZ);

    // strain output *********************************
    if (eastype_!=soh8p_easnone && numgpt_!=8)
    {
      // no stress output currently
    }
    else
    {
      switch (iostrain)
      {
      case INPAR::STR::strain_gl:
      {
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        LINALG::Matrix<numstr_,1> total_glstrain(false);
        total_glstrain(0) = 0.5 * (total_cauchygreen(0,0) - 1.0);
        total_glstrain(1) = 0.5 * (total_cauchygreen(1,1) - 1.0);
        total_glstrain(2) = 0.5 * (total_cauchygreen(2,2) - 1.0);
        total_glstrain(3) = total_cauchygreen(0,1);
        total_glstrain(4) = total_cauchygreen(1,2);
        total_glstrain(5) = total_cauchygreen(2,0);

        if (elestrain == NULL) dserror("strain data not available");
        for (int i = 0; i < 3; ++i)
          (*elestrain)(gp,i) = total_glstrain(i);
        for (int i = 3; i < 6; ++i)
          (*elestrain)(gp,i) = 0.5 * total_glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == NULL) dserror("strain data not available");

        // inverse of deformation gradient
        LINALG::Matrix<3,3> invdefgrd;
        invdefgrd.Invert(defgrd);

        LINALG::Matrix<3,3> total_euler_almansi(true);
        for (int i=0; i<3; i++) total_euler_almansi(i,i) =1.;
        tmp1.MultiplyTN(invdefgrd,invdefgrd);
        total_euler_almansi.MultiplyTN(-1.,invdefgrd,tmp1,1.);
        total_euler_almansi.Scale(0.5);

        (*elestrain)(gp,0) = total_euler_almansi(0,0);
        (*elestrain)(gp,1) = total_euler_almansi(1,1);
        (*elestrain)(gp,2) = total_euler_almansi(2,2);
        (*elestrain)(gp,3) = total_euler_almansi(0,1);
        (*elestrain)(gp,4) = total_euler_almansi(1,2);
        (*elestrain)(gp,5) = total_euler_almansi(0,2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
      {
        dserror("requested strain type not available");
        break;
      }
      }
    }
    // end of strain output **************************

    // variables needed for F-bar
    double detF = -1.;
    LINALG::Matrix<nsd_,nsd_> invdefgrd(false);
    double f_bar_factor=1.;

    // calculate modified deformation gradient
    if (fbar_)
    {
      // inverse and determinant
      detF=invdefgrd.Invert(defgrd);

    // check for element distortion
    if (detF<=0. || detF_0<=0.)
      dserror("element distortion too large");

    // modify deformation gradient
    // ATTENTION: defgrd_mod now contains the F-bar deformation gradient
    // The compatible displacement-based deformation gradient is still
    // available in defgrd.
    f_bar_factor=pow(detF_0/detF,1/3.);
    defgrd_mod.Scale(f_bar_factor);
    }

    // plastic flow increment
    LINALG::Matrix<nsd_,nsd_> deltaLp;

    // recover plastic variables
    if (HavePlasticSpin())
    {
      if (eastype_!=soh8p_easnone)
        RecoverPlasticity<plspin>(res_d,pred,gp,MyPID,deltaLp,lp_inc,(stiffmatrix!=NULL),&delta_alpha_eas);
      else
        RecoverPlasticity<plspin>(res_d,pred,gp,MyPID,deltaLp,lp_inc,(stiffmatrix!=NULL));
    }
    else
    {
      if (eastype_!=soh8p_easnone)
        RecoverPlasticity<zerospin>(res_d,pred,gp,MyPID,deltaLp,lp_inc,(stiffmatrix!=NULL),&delta_alpha_eas);
      else
        RecoverPlasticity<zerospin>(res_d,pred,gp,MyPID,deltaLp,lp_inc,(stiffmatrix!=NULL));
    }

    // material call *********************************************
    LINALG::Matrix<numstr_,1> pk2;
    LINALG::Matrix<numstr_,numstr_> cmat;
    plmat->EvaluateElast(&defgrd_mod,&deltaLp,params,&pk2,&cmat,gp);
    // material call *********************************************

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < numstr_; ++i)
        (*elestress)(gp,i) = pk2(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF = defgrd.Determinant();

      LINALG::Matrix<3,3> pkstress;
      pkstress(0,0) = pk2(0);
      pkstress(0,1) = pk2(3);
      pkstress(0,2) = pk2(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = pk2(1);
      pkstress(1,2) = pk2(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = pk2(2);

      LINALG::Matrix<3,3> cauchystress;
      tmp1.Multiply(1.0/detF,defgrd,pkstress);
      cauchystress.MultiplyNT(tmp1,defgrd);

      (*elestress)(gp,0) = cauchystress(0,0);
      (*elestress)(gp,1) = cauchystress(1,1);
      (*elestress)(gp,2) = cauchystress(2,2);
      (*elestress)(gp,3) = cauchystress(0,1);
      (*elestress)(gp,4) = cauchystress(1,2);
      (*elestress)(gp,5) = cauchystress(0,2);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
    {
      dserror("requested stress type not available");
      break;
    }
    }

    // integrate usual internal force and stiffness matrix
    double detJ_w = detJ*wgt_[gp];
    // integrate elastic internal force vector **************************
    // update internal force vector
    if (force != NULL)
    {
      if (fbar_)
        force->MultiplyTN(detJ_w/f_bar_factor, bop, pk2, 1.0);
      else
        force->MultiplyTN(detJ_w, bop, pk2, 1.0);
    }

    // additional f-bar derivatives
    LINALG::Matrix<numdofperelement_,1> htensor(true);

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<numstr_,numdofperelement_> cb;
      cb.Multiply(cmat,bop);
      if (fbar_)
        stiffmatrix->MultiplyTN(detJ_w*f_bar_factor,bop,cb,1.0);
      else
        stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<numstr_,1> sfac(pk2); // auxiliary integrated stress
      if (fbar_)
        sfac.Scale(detJ_w/f_bar_factor); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      else
        sfac.Scale(detJ_w); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<nen_; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
                             + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
                             + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
                             + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<nen_; ++jnod)
        {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<3; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************

      // integrate additional fbar matrix**************************************
      if (fbar_)
      {
        for(int n=0;n<numdofperelement_;n++)
          for(int i=0;i<3;i++)
            htensor(n) += invdefgrd_0(i,n%3)*N_XYZ_0(i,n/3)-invdefgrd(i,n%3)*N_XYZ(i,n/3);

        LINALG::Matrix<numstr_,1> ccg;
        ccg.Multiply(cmat,RCG);

        LINALG::Matrix<numdofperelement_,1> bopccg(false); // auxiliary integrated stress
        bopccg.MultiplyTN(detJ_w*f_bar_factor/3.0,bop,ccg);

        LINALG::Matrix<numdofperelement_,1> bops(false); // auxiliary integrated stress
        bops.MultiplyTN(-detJ_w/f_bar_factor/3.0,bop,pk2);
        stiffmatrix->MultiplyNT(1.,bops,htensor,1.);
        stiffmatrix->MultiplyNT(1.,bopccg,htensor,1.);
      }
      // end of integrate additional fbar matrix*****************************

      // EAS technology: integrate matrices --------------------------------- EAS
      if (!tang_pred)
        if (eastype_ != soh8p_easnone)
        {
          // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
          // integrate Kda: Kad += (M^T . cmat . B) * detJ * w(gp)
          // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
          LINALG::SerialDenseMatrix cM(numstr_, neas_); // temporary c . M
          switch(eastype_)
          {
          case soh8p_easfull:
            LINALG::DENSEFUNCTIONS::multiply<double,numstr_,numstr_,soh8p_easfull>(cM.A(), cmat.A(), M.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,soh8p_easfull>(1.0, *KaaInv_, detJ_w, M, cM);
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,numdofperelement_>(1.0, Kad_->A(), detJ_w, M.A(), cb.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,numdofperelement_,numstr_,soh8p_easfull>(1.0, Kda.A(), detJ_w,cb.A(), M.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,1>(1.0, feas_->A(), detJ_w, M.A(), pk2.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,1>(1.0, feas_uncondensed.A(), detJ_w, M.A(), pk2.A());
            break;
          case soh8p_easmild:
            LINALG::DENSEFUNCTIONS::multiply<double,numstr_,numstr_,soh8p_easmild>(cM.A(), cmat.A(), M.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,soh8p_easmild>(1.0, *KaaInv_, detJ_w, M, cM);
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,numdofperelement_>(1.0, Kad_->A(), detJ_w, M.A(), cb.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,numdofperelement_,numstr_,soh8p_easmild>(1.0, Kda.A(), detJ_w,cb.A(), M.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,1>(1.0, feas_->A(), detJ_w, M.A(), pk2.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,1>(1.0, feas_uncondensed.A(), detJ_w, M.A(), pk2.A());
            break;
          case soh8p_easnone: break;
          default: dserror("Don't know what to do with EAS type %d", eastype_); break;
          }
        } // ---------------------------------------------------------------- EAS
    } // end of stiffness matrix

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      double density = Material()->Density();
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<nen_; ++inod)
      {
        ifactor = shapefunct(inod) * factor;
        for (int jnod=0; jnod<nen_; ++jnod)
        {
          massfactor = shapefunct(inod) * ifactor;     // intermediate factor
          (*massmatrix)(3*inod+0,3*jnod+0) += massfactor;
          (*massmatrix)(3*inod+1,3*jnod+1) += massfactor;
          (*massmatrix)(3*inod+2,3*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

    // plastic modifications
    if ( (stiffmatrix!=NULL || force!=NULL) && !tang_pred)
    {
      if (HavePlasticSpin())
      {
        if (fbar_)
          CondensePlasticity<plspin>(defgrd_mod,deltaLp,bop,&RCG,MyPID,detJ_w,gp,params,force,stiffmatrix,num_active_gp,lp_res,NULL,NULL,&f_bar_factor,&htensor);
        else if (eastype_!=soh8p_easnone)
          CondensePlasticity<plspin>(defgrd_mod,deltaLp,bop,NULL,MyPID,detJ_w,gp,params,force,stiffmatrix,num_active_gp,lp_res,&M,&Kda);
        else
          CondensePlasticity<plspin>(defgrd_mod,deltaLp,bop,NULL,MyPID,detJ_w,gp,params,force,stiffmatrix,num_active_gp,lp_res);
      }
      else
      {
        if (fbar_)
          CondensePlasticity<zerospin>(defgrd_mod,deltaLp,bop,&RCG,MyPID,detJ_w,gp,params,force,stiffmatrix,num_active_gp,lp_res,NULL,NULL,&f_bar_factor,&htensor);
        else if (eastype_!=soh8p_easnone)
          CondensePlasticity<zerospin>(defgrd_mod,deltaLp,bop,NULL,MyPID,detJ_w,gp,params,force,stiffmatrix,num_active_gp,lp_res,&M,&Kda);
        else
          CondensePlasticity<zerospin>(defgrd_mod,deltaLp,bop,NULL,MyPID,detJ_w,gp,params,force,stiffmatrix,num_active_gp,lp_res);
      }
    }// plastic modifications
  } // gp loop

  // Static condensation EAS --> stiff ********************************
  if (stiffmatrix != NULL && !tang_pred && eastype_!=soh8p_easnone)
  {
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(*KaaInv_);
    solve_for_inverseKaa.Invert();

    tmp.Shape(numdofperelement_,neas_);
    switch (eastype_)
    {
    case soh8p_easfull:
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easfull,soh8p_easfull>
        (0.,tmp.A(),1.,Kda.A(),KaaInv_->A());
      if (stiffmatrix!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easfull,numdofperelement_>
          (1.,stiffmatrix->A(),-1.,tmp.A(),Kad_->A());
      if (force!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easfull,1>
          (1.,force->A(),-1.,tmp.A(),feas_->A());
      break;
    case soh8p_easmild:
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easmild,soh8p_easmild>
        (0.,tmp.A(),1.,Kda.A(),KaaInv_->A());
      if (stiffmatrix!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easmild,numdofperelement_>
          (1.,stiffmatrix->A(),-1.,tmp.A(),Kad_->A());
      if (force!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easmild,1>
          (1.,force->A(),-1.,tmp.A(),feas_->A());
      break;
    case soh8p_easnone:
      break;
    default: dserror("Don't know what to do with EAS type %d", eastype_); break;
    }
    eas_res += pow(feas_uncondensed.Norm2(),2.);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  condense plastic degrees of freedom                     seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
template<int spintype>
void DRT::ELEMENTS::So3_Plast<distype>::CondensePlasticity(
    const LINALG::Matrix<nsd_,nsd_>& defgrd,
    const LINALG::Matrix<nsd_,nsd_>& deltaLp,
    const LINALG::Matrix<numstr_,numdofperelement_>& bop,
    const LINALG::Matrix<numstr_,1>* RCG,
    const int MyPID,
    const double detJ_w,
    const int gp,
    Teuchos::ParameterList& params,
    LINALG::Matrix<numdofperelement_,1>* force,
    LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix,
    int& num_active_gp,
    double& lp_res,
    Epetra_SerialDenseMatrix* M,
    Epetra_SerialDenseMatrix* Kda,
    const double* f_bar_factor,
    const LINALG::Matrix<numdofperelement_,1>* htensor
    )
{
  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // temporary Epetra matrix for matrix-matrix-matrix products
  Epetra_SerialDenseMatrix tmp;

  // second material call ****************************************************
  LINALG::Matrix<numstr_,spintype+1> dpk2ddp;
  LINALG::Matrix<spintype+1,1> ncp;
  LINALG::Matrix<spintype+1,spintype+1> dncpddp;
  LINALG::Matrix<spintype+1,numstr_> dncpdc;
  bool active=false;
  bool elast=false;
  bool as_converged=true;
  plmat->EvaluatePlast(&defgrd,&deltaLp,params,&dpk2ddp,
      &ncp,&dncpdc,&dncpddp,&active,&elast,&as_converged,gp);
  if (MyPID==Owner()) lp_res+=pow(ncp.Norm2(),2.);
  if (active && Owner()==MyPID) ++num_active_gp;
  if (as_converged==false) params.set("unconverged_active_set",true);
  // *************************************************************************

  // Simple matrix do delete the linear dependent row in Voigt-notation.
  // Since all entries are deviatoric, we drop the third entry on the main diagonal.
  LINALG::Matrix<spintype+1,spintype> voigt_red(true);
  switch(spintype)
  {
  case zerospin:
    voigt_red(0,0)=voigt_red(1,1)=voigt_red(3,2)=voigt_red(4,3)=voigt_red(5,4)=1.;
    break;
  case plspin:
    voigt_red(0,0)=voigt_red(1,1)=voigt_red(3,2)=voigt_red(4,3)=voigt_red(5,4)=voigt_red(6,5)=voigt_red(7,6)=voigt_red(8,7)=1.;
    break;
  default:
    dserror("Don't know what to do with plastic spin type %d", spintype);
    break;
  }

  // We have to adapt the derivatives to account for the deviatoric structure
  // of delta D^p.
  LINALG::Matrix<spintype+1,spintype> dDpdbeta;
  switch(spintype)
  {
  case zerospin:
    dDpdbeta(0,0)=dDpdbeta(1,1)=dDpdbeta(5,4)=dDpdbeta(4,3)=dDpdbeta(3,2)=1.;
    dDpdbeta(2,0)=dDpdbeta(2,1)=-1.;
    break;
  case plspin:
    dDpdbeta(0,0)=dDpdbeta(1,1)=+1.;
    dDpdbeta(2,0)=dDpdbeta(2,1)=-1.;
    dDpdbeta(3,2)=1.;
    dDpdbeta(4,3)=1.;
    dDpdbeta(5,4)=1.;
    dDpdbeta(6,2)=1.;
    dDpdbeta(7,3)=1.;
    dDpdbeta(8,4)=1.;

    dDpdbeta(3,5)=1.;
    dDpdbeta(4,6)=1.;
    dDpdbeta(5,7)=1.;
    dDpdbeta(6,5)=-1.;
    dDpdbeta(7,6)=-1.;
    dDpdbeta(8,7)=-1.;
    break;
  default:
    dserror("Don't know what to do with plastic spin type %d", spintype);
    break;
  }

  // derivative of internal force w.r.t. beta
  LINALG::Matrix<numdofperelement_,spintype> kdbeta;
  // derivative of pk2 w.r.t. beta
  LINALG::Matrix<numstr_,spintype> dpk2db;
  dpk2db.Multiply(dpk2ddp,dDpdbeta);
  if (fbar_)
    kdbeta.MultiplyTN(detJ_w/(*f_bar_factor),bop,dpk2db);
  else
    kdbeta.MultiplyTN(detJ_w,bop,dpk2db);

  // EAS matrix block
  Epetra_SerialDenseMatrix Kab(neas_,spintype);
  switch(eastype_)
  {
  case soh8p_easnone:
    // do nothing
    break;
  case soh8p_easmild:
    LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,spintype>(1.,Kab.A(),detJ_w,M->A(),dpk2db.A());
    break;
  case soh8p_easfull:
    LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,spintype>(1.,Kab.A(),detJ_w,M->A(),dpk2db.A());
    break;
  case soh8p_eassosh8:
    LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_eassosh8,numstr_,spintype>(1.,Kab.A(),detJ_w,M->A(),dpk2db.A());
    break;
  default:
    dserror("Don't know what to do with EAS type %d", eastype_);
    break;
  }

  // gauss points that require a full condensation
  if (!elast)
  {
    // apply chain rule for kbb block
    LINALG::Matrix<spintype+1,spintype> dNCPdb;
    dNCPdb.Multiply(dncpddp,dDpdbeta);
    LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,spintype+1,spintype>
      (0.,(*KbbInv_)[gp].A(),1.,voigt_red.A(),dNCPdb.A());

    // apply chain rule for kbd block
    LINALG::Matrix<spintype+1,numdofperelement_> dNCPdd;
    if (fbar_)
    {
      if (RCG==NULL) dserror("CondensePlasticity(...) requires RCG in case of FBAR");
      LINALG::Matrix<spintype+1,1> tmp61;
      tmp61.Multiply(.5,dncpdc,(*RCG));
      dNCPdd.MultiplyNT((*f_bar_factor)*(*f_bar_factor)*2./3.,tmp61,*htensor,0.);
      dNCPdd.Multiply((*f_bar_factor)*(*f_bar_factor),dncpdc,bop,1.);
    }
    else
      dNCPdd.Multiply(dncpdc,bop);
    LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,spintype+1,numdofperelement_>
      (0.,(*Kbd_)[gp].A(),1.,voigt_red.A(),dNCPdd.A());

    // EAS block kba
    if (eastype_!=soh8p_easnone)
    {
      Kba_->at(gp).Shape(spintype,neas_);
      tmp.Shape(spintype+1,neas_);
      switch(eastype_)
      {
      case soh8p_easnone:
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double,spintype+1,numstr_,soh8p_easmild>(0.,tmp.A(),1.,dncpdc.A(),M->A());
        LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,spintype+1,soh8p_easmild>(0.,Kba_->at(gp).A(),1.,voigt_red.A(),tmp.A());
       break;
      case soh8p_easfull:
        LINALG::DENSEFUNCTIONS::multiply<double,spintype+1,numstr_,soh8p_easfull>(0.,tmp.A(),1.,dncpdc.A(),M->A());
        LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,numstr_,soh8p_easfull>(0.,Kba_->at(gp).A(),1.,voigt_red.A(),tmp.A());
        break;
      case soh8p_eassosh8:
        LINALG::DENSEFUNCTIONS::multiply<double,spintype+1,numstr_,soh8p_eassosh8>(0.,tmp.A(),1.,dncpdc.A(),M->A());
        LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,numstr_,soh8p_eassosh8>(0.,Kba_->at(gp).A(),1.,voigt_red.A(),tmp.A());
        break;
      default: dserror("Don't know what to do with EAS type %d", eastype_); break;
      }
    }

    // residual
    LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,spintype+1,1>
      (0.,(*fbeta_)[gp].A(),1.,voigt_red.A(),ncp.A());

    // **************************************************************
    // static condensation of inner variables
    // **************************************************************
    //inverse matrix block [k_beta beta]_ij
    Epetra_SerialDenseSolver solve_for_kbbinv;
    solve_for_kbbinv.SetMatrix((*KbbInv_)[gp]);
    int err = solve_for_kbbinv.Invert();
    if (err != 0)
      dserror("Inversion of Kbb failed");

    // temporary  Kdb.Kbb^-1
    LINALG::Matrix<numdofperelement_,spintype> KdbKbb;
    LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,spintype>
      (0.,KdbKbb.A(),1.,kdbeta.A(),(*KbbInv_)[gp].A());

    // "plastic displacement stiffness"
    // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
    if (stiffmatrix!=NULL)
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,numdofperelement_>
        (1.,stiffmatrix->A(),-1.,KdbKbb.A(),(*Kbd_)[gp].A());

    // "plastic internal force"
    // plFint = [K_db.K_bb^-1].f_b
    if (force!=NULL)
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,1>
        (1.,force->A(),-1.,KdbKbb.A(),(*fbeta_)[gp].A());

    if (eastype_!=soh8p_easnone)
    {
      // condense plasticity into EAS matrix blocks
      tmp.Shape(neas_,spintype);
      switch (eastype_)
      {
      case soh8p_easfull:
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,soh8p_easfull>
          (1.,Kda->A(),-1.,KdbKbb.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,spintype>
          (0.,tmp.A(),1.,Kab.A(),KbbInv_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,numdofperelement_>
          (1.,Kad_->A(),-1.,tmp.A(),Kbd_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,soh8p_easfull>
          (1.,KaaInv_->A(),-1.,tmp.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,1>
          (1.,feas_->A(),-1.,tmp.A(),fbeta_->at(gp).A());
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,soh8p_easmild>
          (1.,Kda->A(),-1.,KdbKbb.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,spintype>
          (0.,tmp.A(),1.,Kab.A(),KbbInv_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,numdofperelement_>
          (1.,Kad_->A(),-1.,tmp.A(),Kbd_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,soh8p_easmild>
          (1.,KaaInv_->A(),-1.,tmp.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,1>
          (1.,feas_->A(),-1.,tmp.A(),fbeta_->at(gp).A());
        break;
      case soh8p_eassosh8:
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,soh8p_eassosh8>
          (1.,Kda->A(),-1.,KdbKbb.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,spintype>
          (0.,tmp.A(),1.,Kab.A(),KbbInv_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,numdofperelement_>
          (1.,Kad_->A(),-1.,tmp.A(),Kbd_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,soh8p_eassosh8>
          (1.,KaaInv_->A(),-1.,tmp.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,1>
          (1.,feas_->A(),-1.,tmp.A(),fbeta_->at(gp).A());
        break;
      case soh8p_easnone:
        // do nothing
        break;
      default: dserror("Don't know what to do with EAS type %d", eastype_); break;
      }
    }
    // **************************************************************
    // static condensation of inner variables
    // **************************************************************
  }

  // simple condensation as kbb is diagonal
  // This destinction is not necessary. However, we use the known diagonal
  // structure to avoid the inversion of the 5x5 matrix kbb.
  else
  {
    if (dDp_last_iter_->at(gp).NormInf()>0.)
    {
      if (force!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,1>
          (1.,force->A(),-1.,kdbeta.A(),(*dDp_last_iter_)[gp].A());

      switch (eastype_)
      {
      case soh8p_easnone:
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,1>
          (1.,feas_->A(),-1.,Kab.A(),(*dDp_last_iter_)[gp].A());
        break;
      case soh8p_easfull:
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,1>
          (1.,feas_->A(),-1.,Kab.A(),(*dDp_last_iter_)[gp].A());
        break;
      case soh8p_eassosh8:
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,1>
          (1.,feas_->A(),-1.,Kab.A(),(*dDp_last_iter_)[gp].A());
        break;
      default: dserror("Don't know what to do with EAS type %d", eastype_); break;
      }
    }

    (*KbbInv_)[gp].Scale(0.);
    (*fbeta_)[gp].Scale(0.);
    (*Kbd_)[gp].Scale(0.);
    if (eastype_!=soh8p_easnone)
      Kba_->at(gp).Shape(spintype,neas_);
    (*dDp_last_iter_)[gp].Scale(0.);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  condense plastic degrees of freedom                     seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
template<int spintype>
void DRT::ELEMENTS::So3_Plast<distype>::RecoverPlasticity(
    const LINALG::Matrix<numdofperelement_,1>& res_d,
    const INPAR::STR::PredEnum& pred,
    const int gp,
    const int MyPID,
    LINALG::Matrix<3,3>& deltaLp,
    double& lp_inc,
    bool recover,
    Epetra_SerialDenseVector* delta_alpha_eas
    )
{
  // recover condensed variables from last iteration step *********
  if (recover)
  {
  // temporary Epetra matrix
    Epetra_SerialDenseVector tmp_v(spintype);
    Epetra_SerialDenseMatrix tmp_m(spintype,numdofperelement_);
    switch(pred)
    {
    // constant predictor
    case INPAR::STR::pred_constdis:
      (*dDp_last_iter_)[gp].Scale(0.);
      break;

    // tangential predictor
    case INPAR::STR::pred_tangdis:
      // do nothing
      break;

    // do usual recovery
    case INPAR::STR::pred_vague:
      // first part
      LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,1>
        (0.,tmp_v.A(),1.,(*KbbInv_)[gp].A(),(*fbeta_)[gp].A());

      // second part
      LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,numdofperelement_>
        (0.,tmp_m.A(),1.,(*KbbInv_)[gp].A(),(*Kbd_)[gp].A());
      LINALG::DENSEFUNCTIONS::multiply<double,spintype,numdofperelement_,1>
        (1.,tmp_v.A(),1.,tmp_m.A(),res_d.A());

      // EAS part
      if (eastype_!=soh8p_easnone)
      {
        tmp_m.Shape(spintype,neas_);
        switch (eastype_)
        {
        case soh8p_easmild:
          LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,soh8p_easmild>(0.,tmp_m.A(),1.,KbbInv_->at(gp).A(),Kba_->at(gp).A());
          LINALG::DENSEFUNCTIONS::multiply<double,spintype,soh8p_easmild,1>(1.,tmp_v.A(),1.,tmp_m.A(),delta_alpha_eas->A());
          break;
        case soh8p_easfull:
          LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,soh8p_easfull>(0.,tmp_m.A(),1.,KbbInv_->at(gp).A(),Kba_->at(gp).A());
          LINALG::DENSEFUNCTIONS::multiply<double,spintype,soh8p_easfull,1>(1.,tmp_v.A(),1.,tmp_m.A(),delta_alpha_eas->A());
          break;
        case soh8p_eassosh8:
          LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,soh8p_eassosh8>(0.,tmp_m.A(),1.,KbbInv_->at(gp).A(),Kba_->at(gp).A());
          LINALG::DENSEFUNCTIONS::multiply<double,spintype,soh8p_eassosh8,1>(1.,tmp_v.A(),1.,tmp_m.A(),delta_alpha_eas->A());
          break;
        case soh8p_easnone:
          break;
        default: dserror("Don't know what to do with EAS type %d", eastype_); break;
        }
      }
      LINALG::DENSEFUNCTIONS::update<double,spintype,1>(1.,(*dDp_last_iter_)[gp],-1.,tmp_v);

      if (MyPID==Owner())
      {
        if (spintype==zerospin)
        {
          lp_inc +=tmp_v(0)*tmp_v(0)
                    +tmp_v(1)*tmp_v(1)
                    +(-tmp_v(0)-tmp_v(1))*(-tmp_v(0)-tmp_v(1))
                    +tmp_v(2)*tmp_v(2)*2.
                    +tmp_v(3)*tmp_v(3)*2.
                    +tmp_v(4)*tmp_v(4)*2.;
        }
        else if (spintype==plspin)
        {
          lp_inc +=tmp_v(0)*tmp_v(0)
                    +tmp_v(1)*tmp_v(1)
                    +(-tmp_v(0)-tmp_v(1))*(-tmp_v(0)-tmp_v(1));
          for (int i=2; i<8; i++)
            lp_inc +=tmp_v(i)*tmp_v(i)*2.;
        }
        else
          dserror("Don't know what to do with Spin type %d", spintype);
      }
      break;

    // unknown predictor type
    default:
      dserror("semi-smooth Newton plasticity algorithm doesn't know "
          "what to do with predictor type %i",pred);
      break;
    }
  }
  // end of recover **********************************************

  // current plastic flow increment
  deltaLp(0,0) = (*dDp_last_iter_)[gp](0);
  deltaLp(1,1) = (*dDp_last_iter_)[gp](1);
  deltaLp(2,2) = -1.0*((*dDp_last_iter_)[gp](0)+(*dDp_last_iter_)[gp](1));
  deltaLp(0,1) = (*dDp_last_iter_)[gp](2);
  deltaLp(1,0) = (*dDp_last_iter_)[gp](2);
  deltaLp(1,2) = (*dDp_last_iter_)[gp](3);
  deltaLp(2,1) = (*dDp_last_iter_)[gp](3);
  deltaLp(0,2) = (*dDp_last_iter_)[gp](4);
  deltaLp(2,0) = (*dDp_last_iter_)[gp](4);
  if (spintype==plspin)
  {
    deltaLp(0,1) += (*dDp_last_iter_)[gp](5);
    deltaLp(1,0) -= (*dDp_last_iter_)[gp](5);
    deltaLp(1,2) += (*dDp_last_iter_)[gp](6);
    deltaLp(2,1) -= (*dDp_last_iter_)[gp](6);
    deltaLp(0,2) += (*dDp_last_iter_)[gp](7);
    deltaLp(2,0) -= (*dDp_last_iter_)[gp](7);
  }
  return;
}



/*----------------------------------------------------------------------*
 |  update plastic deformation for nonlinear kinematics     seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::UpdatePlasticDeformation_nln(PlSpinType spintype)
{
  // loop over all Gauss points
  for (int gp=0; gp<numgpt_; gp++)
  {
    LINALG::Matrix<3,3> deltaLp;
    deltaLp(0,0) = (*dDp_last_iter_)[gp](0);
    deltaLp(1,1) = (*dDp_last_iter_)[gp](1);
    deltaLp(2,2) = -1.0*((*dDp_last_iter_)[gp](0)+(*dDp_last_iter_)[gp](1));
    deltaLp(0,1) = (*dDp_last_iter_)[gp](2);
    deltaLp(1,0) = (*dDp_last_iter_)[gp](2);
    deltaLp(1,2) = (*dDp_last_iter_)[gp](3);
    deltaLp(2,1) = (*dDp_last_iter_)[gp](3);
    deltaLp(0,2) = (*dDp_last_iter_)[gp](4);
    deltaLp(2,0) = (*dDp_last_iter_)[gp](4);
    if (spintype==plspin)
    {
      deltaLp(0,1) += (*dDp_last_iter_)[gp](5);
      deltaLp(1,0) -= (*dDp_last_iter_)[gp](5);
      deltaLp(1,2) += (*dDp_last_iter_)[gp](6);
      deltaLp(2,1) -= (*dDp_last_iter_)[gp](6);
      deltaLp(0,2) += (*dDp_last_iter_)[gp](7);
      deltaLp(2,0) -= (*dDp_last_iter_)[gp](7);
    }
    static_cast<MAT::PlasticElastHyper*>(Material().get())->UpdateGP(gp,&deltaLp);

    (*KbbInv_)[gp].Scale(0.);
    (*Kbd_)[gp].Scale(0.);
    (*fbeta_)[gp].Scale(0.);
  }

  if (eastype_!=soh8p_easnone)
  {
    for (int i=0; i<neas_; i++)
    {
      (*alpha_eas_delta_over_last_timestep_)(i) = (*alpha_eas_)(i)-(*alpha_eas_last_timestep_)(i);
      (*alpha_eas_last_timestep_)(i) = (*alpha_eas_)(i);
    }
    Kad_->Shape(neas_,numdofperelement_);
    KaaInv_->Shape(neas_,neas_);
    feas_->Size(neas_);
  }

  return;
}

