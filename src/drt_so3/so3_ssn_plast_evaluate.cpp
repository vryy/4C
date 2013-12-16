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

// headers of supported hyperelastic-materials

/*----------------------------------------------------------------------*
 | evaluate the element (public)                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Evaluate(
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

    // fbar element?
    DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * eleFBAR
      = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* >(this);

    // default: geometrically non-linear analysis with Total Lagrangean approach
    if (kintype_ == geo_nonlinear)
    {
      if (HaveHillPlasticity())
      {
        if (!eleFBAR)
          nln_stiffmass_hill(la[0].lm_,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
        else
          nln_stiffmassHill_fbar(la[0].lm_,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
      }
      else
      {
        if (!eleFBAR)
          nln_stiffmass(la[0].lm_,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
        else
          nln_stiffmass_fbar(la[0].lm_,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
      }
    }
    // geometric geo_linear
    else if (kintype_ == geo_linear)
    {
      lin_stiffmass(la[0].lm_,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
          INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
      if (HaveHillPlasticity())
        dserror("hill plasticity only for non-linear kinematics");
    }
    else
      dserror("unknown kinematics type");


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

    // fbar element?
    DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * eleFBAR
      = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* >(this);


    // default: geometrically non-linear analysis with Total Lagrangean approach
    if (kintype_ == geo_nonlinear)
    {
      if (HaveHillPlasticity())
      {
        if (!eleFBAR)
          nln_stiffmass_hill(la[0].lm_,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
        else
          nln_stiffmassHill_fbar(la[0].lm_,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
      }
      else
      {
        if (!eleFBAR)
          nln_stiffmass(la[0].lm_,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
        else
          nln_stiffmass_fbar(la[0].lm_,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
              INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
      }
    }
    // geometric geo_linear
    else if (kintype_ == geo_linear)
    {
      lin_stiffmass(la[0].lm_,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
          INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
      if (HaveHillPlasticity())
        dserror("hill plasticity only for non-linear kinematics");
    }
    else
      dserror("unknown kinematics type");

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

     // fbar element?
     DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * eleFBAR
       = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* >(this);

     // default: geometrically non-linear analysis with Total Lagrangean approach
     if (kintype_ == geo_nonlinear)
     {
       if (HaveHillPlasticity())
       {
         if (!eleFBAR)
           nln_stiffmass_hill(la[0].lm_,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
               INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
         else
           nln_stiffmassHill_fbar(la[0].lm_,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
               INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
       }
       else
       {
         if (!eleFBAR)
           nln_stiffmass(la[0].lm_,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
               INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
         else
           nln_stiffmass_fbar(la[0].lm_,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
               INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
       }
     }
     // geometric geo_linear
     else if (kintype_ == geo_linear)
     {
       lin_stiffmass(la[0].lm_,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
           INPAR::STR::stress_none,INPAR::STR::strain_none,discretization.Comm().MyPID());
       if (HaveHillPlasticity())
         dserror("hill plasticity only for non-linear kinematics");
     }
     else
       dserror("unknown kinematics type");

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
    if (discretization.Comm().MyPID() == so3_ele::Owner())
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

      // fbar element?
      DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * eleFBAR
        = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* >(this);


      // default: geometrically non-linear analysis with Total Lagrangean approach
      if (kintype_ == geo_nonlinear)
      {
        if (HaveHillPlasticity())
        {
          if (!eleFBAR)
            nln_stiffmass_hill(la[0].lm_,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,
                iostress,iostrain,discretization.Comm().MyPID());
          else
            nln_stiffmassHill_fbar(la[0].lm_,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,
                iostress,iostrain,discretization.Comm().MyPID());
        }
        else
        {
          if (!eleFBAR)
            nln_stiffmass(la[0].lm_,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,
                iostress,iostrain,discretization.Comm().MyPID());
          else
            nln_stiffmass_fbar(la[0].lm_,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,
                iostress,iostrain,discretization.Comm().MyPID());
        }
      }
      // geometric geo_linear
      else if (kintype_ == geo_linear)
      {
        lin_stiffmass(la[0].lm_,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,
            iostress,iostrain,discretization.Comm().MyPID());
        if (HaveHillPlasticity())
          dserror("hill plasticity only for non-linear kinematics");
      }

      else
        dserror("unknown kinematics type");
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
      if (HaveHillPlasticity())
      {
        mDLp_last_iter_->at(i)=mDLp_last_timestep_->at(i);
        KbbInvHill_->at(i).Clear();
        KbdHill_->at(i).Clear();
        fbetaHill_->at(i).Clear();
      }
      else
      {
        DalphaK_last_iter_->at(i)=DalphaK_last_timestep_->at(i);
        KbbInv_->at(i).Clear();
        Kbd_->at(i).Clear();
        fbeta_->at(i).Clear();
      }
    }
    break;
  }

  //============================================================================
  case calc_struct_update_istep:
  {
    // update plastic deformation
    // default: geometrically non-linear analysis with Total Lagrangean approach
    if (kintype_ == geo_nonlinear)
    {
      if (HaveHillPlasticity())
        UpdatePlasticDeformationHill_nln();
      else
        UpdatePlasticDeformation_nln();
    }
    // geometric geo_linear
    else if (kintype_ == geo_linear)
      UpdatePlasticDeformation_lin();
    else
      dserror("unknown kinematics type");

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
    LINALG::Matrix<numgpt_post,nsd_> gpstress(((*gpstressmap)[gid])->A(),true);

    Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
    if (poststress==Teuchos::null)
      dserror("No element stress/strain vector available");

    if (stresstype=="ndxyz")
    {
      dserror("only element centered stress for so3_plast");
//      // extrapolate stresses/strains at Gauss points to nodes
//      soh8_expol(gpstress, *poststress);
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
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::CalculateBop(
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
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::CalculateBoplin(
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
 | initialise Jacobian                                      seitz 07/13 |
 | is called once in Initialize() in so3_ssn_plast_eletypes.cpp         |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::InitJacobianMapping()
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
    // get the coordinates of Gauss points, here use intrepid
    const double* gpcoord = intpoints_.Point(gp);
    for (int idim=0; idim<nsd_; idim++)
    {
      xsi_[gp](idim) = gpcoord[idim];
    }
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
}  // InitJacobianMapping()

/*----------------------------------------------------------------------*
 | internal force, nonlinear stiffness, and mass            seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::nln_stiffmass(
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
  // update element geometry 3D:
  LINALG::Matrix<nen_,nsd_> xrefe;  // X, material coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;  // x, current  coord. of element

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

  // 3x3 LINALG matrix for temporary stuff
  LINALG::Matrix<nsd_,nsd_> tmp1(false);
  LINALG::Matrix<nsd_,nsd_> tmp2(false);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // get plastic material parameters
  double kinhard=plmat->Kinhard();
  double isohard=plmat->Isohard();
  double expisohard=plmat->Expisohard();
  double infyield=plmat->Infyield();
  double inityield=plmat->Inityield();
  double cpl=plparams_->get<double>("SEMI_SMOOTH_CPL");
  double stab_s=plparams_->get<double>("STABILIZATION_S");

  // converged active set
  bool converged_active_set=true;

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

    // calcualte total rcg/gl/ea for output
    LINALG::Matrix<nsd_,nsd_> total_cauchygreen(false);
    total_cauchygreen.MultiplyTN(defgrd,defgrd);
    // total Cauchy green in voigt notation
    LINALG::Matrix<6,1> RCG;
    for (int i=0; i<3; i++) RCG(i)=total_cauchygreen(i,i);
    RCG(3)=total_cauchygreen(0,1)*2.;
    RCG(4)=total_cauchygreen(1,2)*2.;
    RCG(5)=total_cauchygreen(0,2)*2.;

    // calculate nonlinear B-operator
    LINALG::Matrix<numstr_,numdofperelement_> bop(false);
    CalculateBop(&bop,&defgrd,&N_XYZ);

    // strain output *********************************
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
    }// end of strain output **************************

    // recover condensed variables from last iteration step *********
    // first part
    LINALG::Matrix<5,1> tmp51;
    tmp51.Multiply(KbbInv_->at(gp),fbeta_->at(gp));

    // second part
    LINALG::Matrix<5,numdofperelement_> tmp524;
    tmp524.Multiply(KbbInv_->at(gp),Kbd_->at(gp));
    tmp51.Multiply(1.,tmp524,res_d,1.);
    DalphaK_last_iter_->at(gp).Update(-1.,tmp51,1.);
    double inrement_norm_sqare=tmp51(0)*tmp51(0)
                              +tmp51(1)*tmp51(1)
                              +(-tmp51(0)-tmp51(1))*(-tmp51(0)-tmp51(1))
                              +tmp51(2)*tmp51(2)*2.
                              +tmp51(3)*tmp51(3)*2.
                              +tmp51(4)*tmp51(4)*2.;
    params.get<double>("Lp_increment_square")+=inrement_norm_sqare;
    // end of recover *********************************************

    // inverse plastic deformation gradient
    LINALG::Matrix<nsd_,nsd_> InvPlasticDefgrd(false);
    // inverse plastic deformation gradient at last time step
    LINALG::Matrix<nsd_,nsd_> InvPlasticDefgrdLast = last_plastic_defgrd_inverse_->at(gp);
    // current plastic flow increment
    LINALG::Matrix<nsd_,nsd_> DeltaAlphaK(false);
    DeltaAlphaK(0,0) = DalphaK_last_iter_->at(gp)(0);
    DeltaAlphaK(1,1) = DalphaK_last_iter_->at(gp)(1);
    DeltaAlphaK(2,2) = -1.0*(DalphaK_last_iter_->at(gp)(0)+DalphaK_last_iter_->at(gp)(1));
    DeltaAlphaK(0,1) = DalphaK_last_iter_->at(gp)(2);
    DeltaAlphaK(1,0) = DalphaK_last_iter_->at(gp)(2);
    DeltaAlphaK(1,2) = DalphaK_last_iter_->at(gp)(3);
    DeltaAlphaK(2,1) = DalphaK_last_iter_->at(gp)(3);
    DeltaAlphaK(0,2) = DalphaK_last_iter_->at(gp)(4);
    DeltaAlphaK(2,0) = DalphaK_last_iter_->at(gp)(4);

    tmp1=DeltaAlphaK;
    MatrixExponential3x3(tmp1);
    InvPlasticDefgrd.Multiply(InvPlasticDefgrdLast,tmp1);

    // material call
    LINALG::Matrix<numstr_,1> pk2_stress(true);
    LINALG::Matrix<6,6> Cmat_ABCD(true);
    LINALG::Matrix<6,9> dpk2dfpinv(true);
    LINALG::Matrix<3,3> mandelstress(true);
    LINALG::Matrix<6,6> dmdc(true);
    LINALG::Matrix<6,9> dmdfpinv(true);
    params.set<int>("gp",gp);
    params.set<int>("eleID",Id());
    plmat->Evaluate(&defgrd,&InvPlasticDefgrd,params,&pk2_stress,&Cmat_ABCD,&dpk2dfpinv,&mandelstress,&dmdc,&dmdfpinv);

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < numstr_; ++i)
        (*elestress)(gp,i) = pk2_stress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF = defgrd.Determinant();

      LINALG::Matrix<3,3> pkstress;
      pkstress(0,0) = pk2_stress(0);
      pkstress(0,1) = pk2_stress(3);
      pkstress(0,2) = pk2_stress(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = pk2_stress(1);
      pkstress(1,2) = pk2_stress(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = pk2_stress(2);

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

    // equivalent stress eta
    LINALG::Matrix<nsd_,nsd_> eta(mandelstress);
    for (int i=0; i<nsd_; i++)
      eta(i,i) -= 1./3.*(mandelstress(0,0) + mandelstress(1,1) + mandelstress(2,2));
    eta.Update(2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
    eta.Update(2./3.*kinhard,DeltaAlphaK,1.);
    LINALG::Matrix<5,1> eta_vec(false);
    eta_vec(0) = eta(0,0);
    eta_vec(1) = eta(1,1);
    eta_vec(2) = 0.5*(eta(0,1)+eta(1,0));
    eta_vec(3) = 0.5*(eta(2,1)+eta(1,2));
    eta_vec(4) = 0.5*(eta(0,2)+eta(2,0));

    // eta_trial
    LINALG::Matrix<nsd_,nsd_> eta_trial(eta);
    eta_trial.Update(-1.*cpl,DeltaAlphaK,1.);
    LINALG::Matrix<5,1> eta_trial_vec(false);
    eta_trial_vec(0) = eta_trial(0,0);
    eta_trial_vec(1) = eta_trial(1,1);
    eta_trial_vec(2) = 0.5*(eta_trial(0,1)+eta_trial(1,0));
    eta_trial_vec(3) = 0.5*(eta_trial(2,1)+eta_trial(1,2));
    eta_trial_vec(4) = 0.5*(eta_trial(0,2)+eta_trial(2,0));

    // absolute values
    double absetatrial=eta_trial.Norm2();
    double abseta=eta.Norm2();
    double Dissipation=0.;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        Dissipation-=eta(i,j)*DeltaAlphaK(i,j);

    double deltaAlphaI=0.;
    if (Dissipation>0. && abseta!=0.)
      deltaAlphaI=Dissipation/abseta*sqrt(2./3.);

    // current yield stress equivalent (yield stress scaled by sqrt(2/3))
    double Ypl=0.;
    if (Dissipation>0.)
      Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+ deltaAlphaI)))
          + isohard*(last_alpha_isotropic_->at(gp)(0,0)+ deltaAlphaI) +inityield);
    else
      Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0))))
          + isohard*(last_alpha_isotropic_->at(gp)(0,0)) +inityield);

    // check activity state
    // inactive
    if (Ypl<absetatrial)
    {
      if (activity_state_->at(gp)==false) // gp switches state
        if (abs(Ypl-absetatrial)>AS_CONVERGENCE_TOL*inityield
            || DalphaK_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
          converged_active_set = false;
      activity_state_->at(gp) = true;
    }
    // active
    else
    {
      if (activity_state_->at(gp)==true) // gp switches state
        if (abs(Ypl-absetatrial)>AS_CONVERGENCE_TOL*inityield
            || DalphaK_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
        converged_active_set = false;
      activity_state_->at(gp) = false;
    }

    // integrate usual internal force and stiffness matrix
    double detJ_w = detJ*intpoints_.Weight(gp);
    // integrate elastic internal force vector **************************
    // update internal force vector
    if (force != NULL)
    {
      force->MultiplyTN(detJ_w, bop, pk2_stress, 1.0);
    }

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,numdofperelement_> cb;
      cb.Multiply(Cmat_ABCD,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(pk2_stress); // auxiliary integrated stress
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
          for (int idim=0; idim<nsd_; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************
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
          (*massmatrix)(nsd_*inod+0,nsd_*jnod+0) += massfactor;
          (*massmatrix)(nsd_*inod+1,nsd_*jnod+1) += massfactor;
          (*massmatrix)(nsd_*inod+2,nsd_*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

    // plastic modifications
    if (
        (stiffmatrix!=NULL || force!=NULL)
        &&
        (activity_state_->at(gp)==true || DalphaK_last_iter_->at(gp).NormInf()!=0.)
       )
    {
      // variables needed for condensation and calculated seperately for active and inactive Gauss points
      LINALG::Matrix<5,numdofperelement_> kbetad(true);
      LINALG::Matrix<5,5> kbetabeta(true);
      LINALG::Matrix<5,1> force_beta(true);
      LINALG::Matrix<numdofperelement_,5> kdbeta;

      // 5x5 identitiy matrix
      LINALG::Matrix<5,5> id5(true);
      for (int i=0; i<5; i++)
        id5(i,i)=1.;

        // damping parameter apl
        double apl=1.;
        if (Ypl/abseta<1.) apl=Ypl/abseta;

        // derivative of the matrix exponential
        // for ease of notation, we don't do the derivative w.r.t. DeltaLp
        // but -DeltaLp=DeltaAlphaK
        LINALG::Matrix<6,6> Dexp(false);
        MatrixExponentialDerivativeSym3x3(DeltaAlphaK,Dexp);

        // derivative of inverse plastic deformation gradient w.r.t. flow increment
        LINALG::Matrix<9,5> DFpiDbeta(true);
        for (int a=0; a<3; a++)
          for (int A=0; A<3; A++)
            for (int j=0; j<5; j++)
              for (int b=0; b<3; b++)
              {
                if (j==0 || j==1)
                  DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) +=
                      (Dexp(VOIGT3X3SYM_[a][b],j)-Dexp(VOIGT3X3SYM_[a][b],2))*InvPlasticDefgrdLast(A,b);
                else
                  DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) +=
                      Dexp(VOIGT3X3SYM_[a][b],j+1)*InvPlasticDefgrdLast(A,b)*2.;
              }

        // **************************************************************
        // stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
        // **************************************************************
        // tensor dSdbeta
        // in index notation this contains
        //                   d S_AB
        // dSdbeta_ABj = --------------
        //                  d beta_j
        LINALG::Matrix<6,5> dSdbeta;
        dSdbeta.Multiply(dpk2dfpinv,DFpiDbeta);

        // Calculate stiffness matrix [k^e_{d,beta}]_ij (i=1..numdof; j=1..5)
        kdbeta.MultiplyTN(detJ_w,bop,dSdbeta);
        // **************************************************************
        // end of stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
        // **************************************************************

        if (activity_state_->at(gp)==true || Dissipation>0. )
        {
        // derivative of Mandel stress w.r.t. beta
        //                      d bar Sigma_ab
        // dSigmaDbeta_abj = --------------------
        //                         d beta_j
        LINALG::Matrix<6,5> dsigmadbeta;
        dsigmadbeta.Multiply(dmdfpinv,DFpiDbeta);

        //                    d eta_ab
        // detadbeta_abj = ---------------
        //                    d beta_j
        // in Voigt notation
        // as eta is traceless, there are only 5 components
        LINALG::Matrix<5,5> detadbeta;
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
          {
            // diagonal entries
            if (i==0 || i==1)
              detadbeta(i,j) = dsigmadbeta(i,j)
              - 1./3.* ( dsigmadbeta(0,j) + dsigmadbeta(1,j) + dsigmadbeta(2,j));
            else
              detadbeta(i,j) = dsigmadbeta(i+1,j);
          }
        detadbeta.Update(2./3.*kinhard,id5,1.);

        // calculate derivative detadd
        //             d eta_ab
        // detadd = --------------
        //              d d_i
        LINALG::Matrix<5,numdofperelement_> detadd;

        // derivative of Mandel stress tensor
        LINALG::Matrix<6,numdofperelement_> dSigmadd(true);
        dSigmadd.Multiply(dmdc,bop);

        for (int i=0; i<5; i++)
          for (int j=0; j<numdofperelement_; j++)
          {
            // diagonal entries
            if (i==0 || i==1)
              detadd(i,j) = dSigmadd(i,j) - 1./3. * (dSigmadd(0,j) + dSigmadd(1,j) + dSigmadd(2,j));
            else
              detadd(i,j) = dSigmadd(i+1,j);
          }

        LINALG::Matrix<5,1> dYpl_dbeta(true);
        LINALG::Matrix<numdofperelement_,1> dYpl_dd(true);
        if (Dissipation>0.)
        {
          for (int j=0;j<2;j++)
          {
            for (int i=0; i<numdofperelement_;i++)
            {
              dYpl_dd(i) += -Dissipation/abseta/abseta*(2.*eta_vec(j) + eta_vec((j+1)%2))/abseta*detadd(j,i);
              dYpl_dd(i) -= (2.*DalphaK_last_iter_->at(gp)(j) + DalphaK_last_iter_->at(gp)((j+1)%2))*detadd(j,i)/abseta;
            }
            for (int i=0; i<5;i++)
            {
              dYpl_dbeta(i) += -Dissipation/abseta/abseta*(2.*eta_vec(j) + eta_vec((j+1)%2))/abseta*detadbeta(j,i);
              dYpl_dbeta(i) -= (2.*DalphaK_last_iter_->at(gp)(j)+DalphaK_last_iter_->at(gp)((j+1)%2))*detadbeta(j,i)/abseta;
            }
            dYpl_dbeta(j) -= (2.*eta_vec(j) + eta_vec((j+1)%2))/abseta;
          }
          for (int j=2; j<5; j++)
          {
            for (int i=0; i<numdofperelement_; i++)
            {
              dYpl_dd(i) += -Dissipation/abseta/abseta*2.*eta_vec(j)*detadd(j,i)/abseta;
              dYpl_dd(i) -= 2.*DalphaK_last_iter_->at(gp)(j)*detadd(j,i)/abseta;
            }
            for (int i=0; i<5; i++)
            {
              dYpl_dbeta(i) += -Dissipation/abseta/abseta*2.*eta_vec(j)*detadbeta(j,i)/abseta;
              dYpl_dbeta(i) -= 2.*DalphaK_last_iter_->at(gp)(j)*detadbeta(j,i)/abseta;
            }
            dYpl_dbeta(j) -= 2.*eta_vec(j)/abseta;
          }
          dYpl_dd.Scale(2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI))));
          dYpl_dbeta.Scale(2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI))));
        }

        if (activity_state_->at(gp)==true)
        {
          // communicate number of active plastic gauss points back to time integration
          // don't sum up for ghost elements
          if (MyPID == so3_ele::Owner()) params.get<int>("number_active_plastic_gp")++;

          // **************************************************************
          // stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
          // linearization of the complementarity function
          // **************************************************************

          LINALG::Matrix<5,1> dabs_eta_trial_dbeta(true);
          for (int j=0;j<2; j++)
            for (int i=0; i<5; i++)
              dabs_eta_trial_dbeta(i) += (2.*eta_trial_vec(j)+ eta_trial_vec((j+1)%2))*(detadbeta(j,i)-id5(i,j)*cpl)/absetatrial;
          for (int j=2; j<5; j++)
            for (int i=0; i<5; i++)
              dabs_eta_trial_dbeta(i) += (2.*eta_trial_vec(j))*(detadbeta(j,i)-id5(i,j)*cpl)/absetatrial;

          // build kbb from the previous linearizations
          KbbInv_->at(gp).Clear();
          KbbInv_->at(gp).Update(1.-Ypl/absetatrial,detadbeta,1.);
          KbbInv_->at(gp).Update(Ypl/absetatrial*cpl,id5,1.);
          KbbInv_->at(gp).MultiplyNT(-1./absetatrial,eta_trial_vec,dYpl_dbeta,1.);
          KbbInv_->at(gp).MultiplyNT((1-stab_s)*Ypl/pow(absetatrial,2.),eta_trial_vec,dabs_eta_trial_dbeta,1.);
          KbbInv_->at(gp).MultiplyNT(stab_s*apl/absetatrial,eta_vec,dabs_eta_trial_dbeta,1.);

          // **************************************************************
          // end of stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
          // **************************************************************

          // **************************************************************
          // stiffness matrix [k^e_{beta d}]_ij (i=1..5; j=1..numdof)
          // linearization of the complementarity function
          // **************************************************************

          LINALG::Matrix<numdofperelement_,1> dabs_eta_trial_dd(true);
          for (int j=0;j<2; j++)
            for (int i=0; i<numdofperelement_; i++)
              dabs_eta_trial_dd(i) += (2.*eta_trial_vec(j)+ eta_trial_vec((j+1)%2))*(detadd(j,i))/absetatrial;
          for (int j=2; j<5; j++)
            for (int i=0; i<numdofperelement_; i++)
              dabs_eta_trial_dd(i) += (2.*eta_trial_vec(j))*(detadd(j,i))/absetatrial;

          Kbd_->at(gp).Clear();
          Kbd_->at(gp).Update(1.-1.*Ypl/absetatrial,detadd,1.);
          Kbd_->at(gp).MultiplyNT(stab_s*apl/absetatrial,eta_vec,dabs_eta_trial_dd,1.);
          Kbd_->at(gp).MultiplyNT((1.-stab_s)*Ypl/absetatrial/absetatrial,eta_trial_vec,dabs_eta_trial_dd,1.);
          Kbd_->at(gp).MultiplyNT(-1./absetatrial,eta_trial_vec,dYpl_dd,1.);
          // **************************************************************
          // end of stiffness matrix [k^e_{beta d}]_ij (i=1..5; j=1..numdof)
          // **************************************************************

          // **************************************************************
          // right hand side term for complementarity function f^int_i (i=1..5)
          // **************************************************************
          fbeta_->at(gp).Update(eta_vec);
          fbeta_->at(gp).Update(-1.*Ypl/absetatrial,eta_trial_vec,1.);
          // **************************************************************
          // end of right hand side term for complementarity function f^int_i (i=1..5)
          // **************************************************************
        } // acrtive Gauss points

        // inactive Gauss point with plastic history within this load/time step
        else if (activity_state_->at(gp)==false && Dissipation>0.)
        {
          // the complementarity function is
          // C^pl = - Ypl^s * cplparam_ * delta alpha^k

          // Complementarity function independent from displacements
          Kbd_->at(gp).Clear();
          Kbd_->at(gp).MultiplyNT(stab_s*cpl/Ypl,DalphaK_last_iter_->at(gp),dYpl_dd,1.);

          KbbInv_->at(gp).Update(cpl,id5,0.);
          KbbInv_->at(gp).MultiplyNT(stab_s/Ypl * cpl,DalphaK_last_iter_->at(gp),dYpl_dbeta,1.);

          // right hand side term
          fbeta_->at(gp).Update(cpl,DalphaK_last_iter_->at(gp),0.);
        }

        // **************************************************************
        // static condensation of inner variables
        // **************************************************************
        //inverse matrix block [k_beta beta]_ij
        Epetra_SerialDenseMatrix Kbetabeta_epetra(5,5);
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            Kbetabeta_epetra(i,j) = KbbInv_->at(gp)(i,j);
        // we need the inverse of K_beta beta
        Epetra_SerialDenseSolver solve_for_inverseKbb;
        solve_for_inverseKbb.SetMatrix(Kbetabeta_epetra);
        solve_for_inverseKbb.Invert();
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            KbbInv_->at(gp)(i,j) = Kbetabeta_epetra(i,j);

        LINALG::Matrix<numdofperelement_,5> KdbKbb; // temporary  Kdb.Kbb^-1
        KdbKbb.Multiply(kdbeta,KbbInv_->at(gp));

        // "plastic displacement stiffness"
        // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
        if (stiffmatrix!=NULL) stiffmatrix->Multiply(-1.,KdbKbb,Kbd_->at(gp),1.);

        // "plastic internal force"
        // plFint = [K_db.K_bb^-1].f_b
        if (force!=NULL) force->Multiply(-1.,KdbKbb,fbeta_->at(gp),1.);
      }
        else if (activity_state_->at(gp)==false && Dissipation<0.)
        {
          // C^{pl} = cpl*DeltaAlphaK
          // independent of displacements
          Kbd_->at(gp).Clear();
          // We can state the inverse right away
          KbbInv_->at(gp).Update(id5);
          // right hand side
          fbeta_->at(gp).Update(DalphaK_last_iter_->at(gp));

          // condensation to internal force vector
          if (force!=NULL)
            force->Multiply(-1.,kdbeta,fbeta_->at(gp),1.);
        }
    } // modification for plastic Gauss points
    // reset for elastic GP with no plastic flow increment in this time step
    else
    {
      Kbd_->at(gp).Clear();
      fbeta_->at(gp).Clear();
      KbbInv_->at(gp).Clear();
    }

    // square of the residual L2 norm
    double residual_norm_sqare=fbeta_->at(gp)(0)*fbeta_->at(gp)(0)
                              +fbeta_->at(gp)(1)*fbeta_->at(gp)(1)
                              +(-fbeta_->at(gp)(0)-fbeta_->at(gp)(1))*(-fbeta_->at(gp)(0)-fbeta_->at(gp)(1))
                              +fbeta_->at(gp)(2)*fbeta_->at(gp)(2)*2.
                              +fbeta_->at(gp)(3)*fbeta_->at(gp)(3)*2.
                              +fbeta_->at(gp)(4)*fbeta_->at(gp)(4)*2.;
        params.get<double>("Lp_residual_square")+=residual_norm_sqare;

  } // gp loop

  // communicate unconverged active set to time integration
  if (converged_active_set==false)
    params.set("unconverged_active_set",true);

  return;
    }

/*----------------------------------------------------------------------*
 | internal force, nonlinear stiffness, and mass            seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::nln_stiffmass_hill(
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
  // update element geometry 3D:
  LINALG::Matrix<nen_,nsd_> xrefe;  // X, material coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;  // x, current  coord. of element

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

  // 3x3 LINALG matrix for temporary stuff
  LINALG::Matrix<nsd_,nsd_> tmp1(false);
  LINALG::Matrix<nsd_,nsd_> tmp2(false);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // get plastic material parameters
  double kinhard=plmat->Kinhard();
  double isohard=plmat->Isohard();
  double expisohard=plmat->Expisohard();
  double infyield=plmat->Infyield();
  double inityield=plmat->Inityield();
  double cpl=plparams_->get<double>("SEMI_SMOOTH_CPL");
  double stab_s=plparams_->get<double>("STABILIZATION_S");
  LINALG::Matrix<5,5> PlAniso(plmat->PlAniso());
  LINALG::Matrix<5,5> InvPlAniso(plmat->InvPlAniso());
  double PlSpinEta=-plmat->PlSpinEta();

  // converged active set
  bool converged_active_set=true;

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

    // calcualte total rcg/gl/ea for output
    LINALG::Matrix<nsd_,nsd_> total_cauchygreen(false);
    total_cauchygreen.MultiplyTN(defgrd,defgrd);
    // total Cauchy green in voigt notation
    LINALG::Matrix<6,1> RCG;
    for (int i=0; i<3; i++) RCG(i)=total_cauchygreen(i,i);
    RCG(3)=total_cauchygreen(0,1)*2.;
    RCG(4)=total_cauchygreen(1,2)*2.;
    RCG(5)=total_cauchygreen(0,2)*2.;

    // calculate nonlinear B-operator
    LINALG::Matrix<numstr_,numdofperelement_> bop(false);
    CalculateBop(&bop,&defgrd,&N_XYZ);

    // strain output *********************************
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
    }// end of strain output **************************

    // recover condensed variables from last iteration step *********
    // first part
    LINALG::Matrix<8,1> tmp81;
    tmp81.Multiply(KbbInvHill_->at(gp),fbetaHill_->at(gp));

    // second part
    LINALG::Matrix<8,numdofperelement_> tmp824;
    tmp824.Multiply(KbbInvHill_->at(gp),KbdHill_->at(gp));
    tmp81.Multiply(1.,tmp824,res_d,1.);
    mDLp_last_iter_->at(gp).Update(-1.,tmp81,1.);
    double increment_norm_sqare=tmp81(0)*tmp81(0)
                              +tmp81(1)*tmp81(1)
                              +(-tmp81(0)-tmp81(1))*(-tmp81(0)-tmp81(1));
    for (int i=2; i<8; i++) increment_norm_sqare+=tmp81(i)*tmp81(i)*2.;
    params.get<double>("Lp_increment_square")+=increment_norm_sqare;
    // end of recover *********************************************

    // inverse plastic deformation gradient
    LINALG::Matrix<nsd_,nsd_> InvPlasticDefgrd(false);
    // inverse plastic deformation gradient at last time step
    LINALG::Matrix<nsd_,nsd_> InvPlasticDefgrdLast = last_plastic_defgrd_inverse_->at(gp);

    // current plastic flow increment =-Delta Lp
    LINALG::Matrix<nsd_,nsd_> mDLp(false);
    // current kinematic hardening increment
    LINALG::Matrix<nsd_,nsd_> DeltaAlphaK(false);
    LINALG::Matrix<5,1> DeltaAlphaK_vec(false);
    for (int i=0; i<5; i++) DeltaAlphaK_vec(i) = mDLp_last_iter_->at(gp)(i);
    // -Dp
    DeltaAlphaK(0,0) = mDLp_last_iter_->at(gp)(0);
    DeltaAlphaK(1,1) = mDLp_last_iter_->at(gp)(1);
    DeltaAlphaK(2,2) = -1.0*(mDLp_last_iter_->at(gp)(0)+mDLp_last_iter_->at(gp)(1));
    DeltaAlphaK(0,1) = mDLp_last_iter_->at(gp)(2);
    DeltaAlphaK(1,0) = mDLp_last_iter_->at(gp)(2);
    DeltaAlphaK(1,2) = mDLp_last_iter_->at(gp)(3);
    DeltaAlphaK(2,1) = mDLp_last_iter_->at(gp)(3);
    DeltaAlphaK(0,2) = mDLp_last_iter_->at(gp)(4);
    DeltaAlphaK(2,0) = mDLp_last_iter_->at(gp)(4);
    mDLp.Update(DeltaAlphaK);
    // Wp
    LINALG::Matrix<3,3> Wp(true);
    Wp(0,1) -= mDLp_last_iter_->at(gp)(5);
    Wp(1,0) += mDLp_last_iter_->at(gp)(5);
    Wp(1,2) -= mDLp_last_iter_->at(gp)(6);
    Wp(2,1) += mDLp_last_iter_->at(gp)(6);
    Wp(0,2) -= mDLp_last_iter_->at(gp)(7);
    Wp(2,0) += mDLp_last_iter_->at(gp)(7);
    mDLp.Update(-1.,Wp,1.);

    // compute matrix exponential
    tmp1=mDLp;
    MatrixExponential3x3(tmp1);
    InvPlasticDefgrd.Multiply(InvPlasticDefgrdLast,tmp1);

    // material call
    LINALG::Matrix<numstr_,1> pk2_stress(true);
    LINALG::Matrix<6,6> Cmat_ABCD(true);
    LINALG::Matrix<6,9> dpk2dfpinv(true);
    LINALG::Matrix<3,3> mandelstress(true);
    LINALG::Matrix<6,6> dmdc(true);
    LINALG::Matrix<6,9> dmdfpinv(true);
    params.set<int>("gp",gp);
    params.set<int>("eleID",Id());
    plmat->Evaluate(&defgrd,&InvPlasticDefgrd,params,&pk2_stress,&Cmat_ABCD,&dpk2dfpinv,&mandelstress,&dmdc,&dmdfpinv);

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < numstr_; ++i)
        (*elestress)(gp,i) = pk2_stress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF = defgrd.Determinant();

      LINALG::Matrix<3,3> pkstress;
      pkstress(0,0) = pk2_stress(0);
      pkstress(0,1) = pk2_stress(3);
      pkstress(0,2) = pk2_stress(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = pk2_stress(1);
      pkstress(1,2) = pk2_stress(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = pk2_stress(2);

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

    // equivalent stress eta
    LINALG::Matrix<nsd_,nsd_> eta(mandelstress);
    for (int i=0; i<nsd_; i++)
      eta(i,i) -= 1./3.*(mandelstress(0,0) + mandelstress(1,1) + mandelstress(2,2));
    eta.Update(2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
    eta.Update(2./3.*kinhard,DeltaAlphaK,1.);
    LINALG::Matrix<5,1> eta_vec(false);
    eta_vec(0) = eta(0,0);
    eta_vec(1) = eta(1,1);
    eta_vec(2) = 0.5*(eta(0,1)+eta(1,0));
    eta_vec(3) = 0.5*(eta(2,1)+eta(1,2));
    eta_vec(4) = 0.5*(eta(0,2)+eta(2,0));
    LINALG::Matrix<5,1> Aeta(false);
    Aeta.Multiply(PlAniso,eta_vec);

    // eta_trial
    LINALG::Matrix<5,1> eta_trial_vec(eta_vec);
    eta_trial_vec.Multiply(-cpl,InvPlAniso,DeltaAlphaK_vec,1.);
    LINALG::Matrix<5,1> AetaTr(false);
    AetaTr.Multiply(PlAniso,eta_trial_vec);

    // Matrix norms
    double AnormEtatrial=0.; // sqrt( eta_tr : A : eta_tr )
    LINALG::Matrix<1,1> tmp11;
    LINALG::Matrix<5,1> tmp51;
    tmp1(0,0) = 2./3.*AetaTr(0)-1./3.*AetaTr(1);
    tmp1(1,1) =-1./3.*AetaTr(0)+2./3.*AetaTr(1);
    tmp1(2,2) = -tmp1(0,0)-tmp1(1,1);
    tmp1(0,1) = AetaTr(2)/2.;
    tmp1(1,0) = AetaTr(2)/2.;
    tmp1(1,2) = AetaTr(3)/2.;
    tmp1(2,1) = AetaTr(3)/2.;
    tmp1(0,2) = AetaTr(4)/2.;
    tmp1(2,0) = AetaTr(4)/2.;
    tmp11.MultiplyTN(eta_trial_vec,AetaTr);
    AnormEtatrial=sqrt(tmp11(0,0));

    tmp1(0,0) = 2./3.*Aeta(0)-1./3.*Aeta(1);
    tmp1(1,1) =-1./3.*Aeta(0)+2./3.*Aeta(1);
    tmp1(2,2) = -tmp1(0,0)-tmp1(1,1);
    tmp1(0,1) = Aeta(2)/2.;
    tmp1(1,0) = Aeta(2)/2.;
    tmp1(1,2) = Aeta(3)/2.;
    tmp1(2,1) = Aeta(3)/2.;
    tmp1(0,2) = Aeta(4)/2.;
    tmp1(2,0) = Aeta(4)/2.;
    double EucNormAeta = tmp1.Norm2();
    tmp11.MultiplyTN(eta_vec,Aeta);
    double AnormEta=sqrt(tmp11(0,0));

    double DpAe=0.;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        DpAe+=DeltaAlphaK(i,j)*tmp1(i,j);

    if (DpAe<0.)
      deltaAlphaI_->at(gp) = -sqrt(2./3.)*DpAe*AnormEta/EucNormAeta/EucNormAeta;

    double Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_->at(gp))))
        + isohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_->at(gp)) +inityield);

      // check activity state
      // inactive
      if (Ypl<AnormEtatrial)
      {
        if (activity_state_->at(gp)==false) // gp switches state
          if (abs(Ypl-AnormEtatrial)>AS_CONVERGENCE_TOL*inityield
              || mDLp_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
            converged_active_set = false;
        activity_state_->at(gp) = true;
      }
      // active
      else
      {
        if (activity_state_->at(gp)==true) // gp switches state
          if (abs(Ypl-AnormEtatrial)>AS_CONVERGENCE_TOL*inityield
              || mDLp_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
            converged_active_set = false;
        activity_state_->at(gp) = false;
      }

      // integrate usual internal force and stiffness matrix
      double detJ_w = detJ*intpoints_.Weight(gp);
      // integrate elastic internal force vector **************************
      // update internal force vector
      if (force != NULL)
      {
        force->MultiplyTN(detJ_w, bop, pk2_stress, 1.0);
      }

      // update stiffness matrix
      if (stiffmatrix != NULL)
      {
        // integrate `elastic' and `initial-displacement' stiffness matrix
        // keu = keu + (B^T . C . B) * detJ * w(gp)
        LINALG::Matrix<6,numdofperelement_> cb;
        cb.Multiply(Cmat_ABCD,bop);
        stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

        // integrate `geometric' stiffness matrix and add to keu *****************
        LINALG::Matrix<6,1> sfac(pk2_stress); // auxiliary integrated stress
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
            for (int idim=0; idim<nsd_; ++idim)
              bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
            (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
            (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
            (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
          }
        } // end of integrate `geometric' stiffness******************************
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
            (*massmatrix)(nsd_*inod+0,nsd_*jnod+0) += massfactor;
            (*massmatrix)(nsd_*inod+1,nsd_*jnod+1) += massfactor;
            (*massmatrix)(nsd_*inod+2,nsd_*jnod+2) += massfactor;
          }
        }
      } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++


      // plastic modifications
      if (
          (stiffmatrix!=NULL || force!=NULL)
          &&
          (activity_state_->at(gp)==true || mDLp_last_iter_->at(gp).NormInf()!=0.)
         )
      {
          // variables needed for condensation and calculated seperately for active and inactive Gauss points
          LINALG::Matrix<8,numdofperelement_> kbetad(true);
          LINALG::Matrix<8,8> kbetabeta(true);
          LINALG::Matrix<8,1> force_beta(true);
          LINALG::Matrix<numdofperelement_,8> kdbeta;

          // derivative of symmetric complementarity function w.r.t. displacements
          LINALG::Matrix<5,numdofperelement_> DcplSymDd(true);
          // derivative of symmetric complementarity function
          LINALG::Matrix<5,8> DcplSymDbeta(true);
          // complementarity function
          LINALG::Matrix<5,1> CplSym(true);
          // plastic spin equation
          LINALG::Matrix<3,1> SpEq(true);
          // derivative of spin equation w.r.t. displacements
          LINALG::Matrix<3,numdofperelement_> dSpEqDd(true);
          // derivative of spin equation w.r.t. beta
          LINALG::Matrix<3,8> dSpEqDbeta(true);

          // delta alpha_K d beta
          LINALG::Matrix<5,8> DalphaKdbeta(true);
          for (int i=0; i<5; i++)
            DalphaKdbeta(i,i)=1.;
          // damping parameter apl
          double apl=1.;
          if (Ypl/AnormEta<1.) apl=Ypl/AnormEta;

          // derivative of the matrix exponential
          LINALG::Matrix<9,9> Dexp(false);
          MatrixExponentialDerivative3x3(mDLp,Dexp);

          // derivative of inverse plastic deformation gradient w.r.t. flow increment
          LINALG::Matrix<9,8> DFpiDbeta(true);
          for (int a=0; a<3; a++)
            for (int A=0; A<3; A++)
              for (int j=0; j<8; j++)
                for (int b=0; b<3; b++)
                {
                  // main diagonal entries
                  if (j==0 || j==1)
                    DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) +=
                    (Dexp(VOIGT3X3NONSYM_[b][a],j)-Dexp(VOIGT3X3NONSYM_[b][a],2))*InvPlasticDefgrdLast(A,b);
                  if (j==2 || j==3 || j==4)
                    DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) += (Dexp(VOIGT3X3NONSYM_[b][a],j+1)+Dexp(VOIGT3X3NONSYM_[b][a],j+4))*InvPlasticDefgrdLast(A,b);
                  if (j==5 ||j==6 || j==7)
                    DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) += (Dexp(VOIGT3X3NONSYM_[b][a],j-2)-Dexp(VOIGT3X3NONSYM_[b][a],j+1))*InvPlasticDefgrdLast(A,b);
                }

          // **************************************************************
          // stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
          // **************************************************************
          // tensor dSdbeta
          // in index notation this contains
          //                   d S_AB
          // dSdbeta_ABj = --------------
          //                  d beta_j
          LINALG::Matrix<6,8> dSdbeta;
          dSdbeta.Multiply(dpk2dfpinv,DFpiDbeta);

          // Calculate stiffness matrix [k^e_{d,beta}]_ij (i=1..numdof; j=1..5)
          kdbeta.MultiplyTN(detJ_w,bop,dSdbeta);
          // **************************************************************
          // end of stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
          // **************************************************************

          // symmetric NCP stuff
          if (activity_state_->at(gp)==true || DpAe<0. )
          {
          // derivative of Mandel stress w.r.t. beta
          //                      d bar Sigma_ab
          // dSigmaDbeta_abj = --------------------
          //                         d beta_j
          LINALG::Matrix<6,8> dsigmadbeta;
          dsigmadbeta.Multiply(dmdfpinv,DFpiDbeta);

          //                    d eta_ab
          // detadbeta_abj = ---------------
          //                    d beta_j
          // in Voigt notation
          // as eta is traceless, there are only 5 components
          LINALG::Matrix<5,8> detadbeta;
          LINALG::Matrix<5,8> detatrialdbeta;
          for (int i=0; i<5; i++)
            for (int j=0; j<8; j++)
            {
              // diagonal entries
              if (i==0 || i==1)
              {
                detadbeta(i,j) = dsigmadbeta(i,j)
                                         - 1./3.* ( dsigmadbeta(0,j) + dsigmadbeta(1,j) + dsigmadbeta(2,j)) + (2./3.*kinhard*(i==j));
              }
              else
              {
                detadbeta(i,j) = dsigmadbeta(i+1,j) + (2./3.*kinhard*(i==j));
              }
              detatrialdbeta(i,j) = detadbeta(i,j);
              if (j<5)
                detatrialdbeta(i,j) -= cpl*InvPlAniso(i,j);
            }

          // calculate derivative detadd
          //             d eta_ab
          // detadd = --------------
          //              d d_i
          LINALG::Matrix<5,numdofperelement_> detadd;

          // derivative of Mandel stress tensor
          LINALG::Matrix<6,numdofperelement_> dSigmadd(true);
          dSigmadd.Multiply(dmdc,bop);

          for (int i=0; i<5; i++)
            for (int j=0; j<numdofperelement_; j++)
            {
              // diagonal entries
              if (i==0 || i==1)
                detadd(i,j) = dSigmadd(i,j) - 1./3. * (dSigmadd(0,j) + dSigmadd(1,j) + dSigmadd(2,j));
              else
                detadd(i,j) = dSigmadd(i+1,j);
            }

           LINALG::Matrix<8,1> dYpl_dbeta(true);
           LINALG::Matrix<numdofperelement_,1> dYpl_dd(true);
           if (DpAe<0.)
           {
             LINALG::Matrix<5,8> dAetaDbeta;
             dAetaDbeta.Multiply(PlAniso,detadbeta);
             LINALG::Matrix<5,numdofperelement_> dAetaDd;
             dAetaDd.Multiply(PlAniso,detadd);
             LINALG::Matrix<8,1> dAnormEtaDbeta(true);
              dAnormEtaDbeta.MultiplyTN(1./AnormEta,detadbeta,Aeta);
              LINALG::Matrix<numdofperelement_,1> dAnormEtaDd(true);
               dAnormEtaDd.MultiplyTN(1./AnormEta,detadd,Aeta);
             for (int i=0; i<8; i++)
             {
               dYpl_dbeta(i) += 4./3.*Aeta(0)*dAetaDbeta(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dbeta(i) -= 2./3.*Aeta(1)*dAetaDbeta(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dbeta(i) -= 2./3.*Aeta(0)*dAetaDbeta(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dbeta(i) += 4./3.*Aeta(1)*dAetaDbeta(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dbeta(i) += Aeta(2)*dAetaDbeta(2,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dbeta(i) += Aeta(3)*dAetaDbeta(3,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dbeta(i) += Aeta(4)*dAetaDbeta(4,i) *AnormEta*DpAe/pow(EucNormAeta,4.);

               dYpl_dbeta(i) -= DeltaAlphaK_vec(0)*dAetaDbeta(0,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(i) -= DeltaAlphaK_vec(2)*dAetaDbeta(2,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(i) -= DeltaAlphaK_vec(3)*dAetaDbeta(3,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(i) -= DeltaAlphaK_vec(4)*dAetaDbeta(4,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(i) -= DeltaAlphaK_vec(1)*dAetaDbeta(1,i)*AnormEta/EucNormAeta/EucNormAeta;
             }

             dYpl_dbeta(0) -= Aeta(0)*AnormEta/EucNormAeta/EucNormAeta;
             dYpl_dbeta(1) -= Aeta(1)*AnormEta/EucNormAeta/EucNormAeta;
             dYpl_dbeta(2) -= Aeta(2)*AnormEta/EucNormAeta/EucNormAeta;
             dYpl_dbeta(3) -= Aeta(3)*AnormEta/EucNormAeta/EucNormAeta;
             dYpl_dbeta(4) -= Aeta(4)*AnormEta/EucNormAeta/EucNormAeta;

             dYpl_dbeta.Update(-DpAe/EucNormAeta/EucNormAeta,dAnormEtaDbeta,1.);

             for (int i=0; i<numdofperelement_; i++)
             {
               dYpl_dd(i) += 4./3.*Aeta(0)*dAetaDd(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dd(i) -= 2./3.*Aeta(1)*dAetaDd(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dd(i) -= 2./3.*Aeta(0)*dAetaDd(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dd(i) += 4./3.*Aeta(1)*dAetaDd(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dd(i) += Aeta(2)*dAetaDd(2,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dd(i) += Aeta(3)*dAetaDd(3,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
               dYpl_dd(i) += Aeta(4)*dAetaDd(4,i) *AnormEta*DpAe/pow(EucNormAeta,4.);

               dYpl_dd(i) -= DeltaAlphaK_vec(0)*dAetaDd(0,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dd(i) -= DeltaAlphaK_vec(2)*dAetaDd(2,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dd(i) -= DeltaAlphaK_vec(3)*dAetaDd(3,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dd(i) -= DeltaAlphaK_vec(4)*dAetaDd(4,i)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dd(i) -= DeltaAlphaK_vec(1)*dAetaDd(1,i)*AnormEta/EucNormAeta/EucNormAeta;
             }

             dYpl_dd.Update(-DpAe/EucNormAeta/EucNormAeta,dAnormEtaDd,1.);

             dYpl_dbeta.Scale(2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI_->at(gp)))));
             dYpl_dd.Scale(   2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI_->at(gp)))));

           }

          if (activity_state_->at(gp)==true)
          {
            // communicate number of active plastic gauss points back to time integration
            // don't sum up for ghost elements
            if (MyPID == so3_ele::Owner()) params.get<int>("number_active_plastic_gp")++;


            CplSym.Update(1.,eta_vec,0.);
            CplSym.Update(-Ypl/AnormEtatrial,eta_trial_vec,1.);

            // **************************************************************
            // stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
            // linearization of the complementarity function
            // **************************************************************

            LINALG::Matrix<5,1> Aetatr;
            Aetatr.Multiply(PlAniso,eta_trial_vec);
            LINALG::Matrix<8,1> dAnormEtaTrialDbeta(true);
            dAnormEtaTrialDbeta.MultiplyTN(1./AnormEtatrial,detatrialdbeta,Aetatr);

            LINALG::Matrix<5,5> etasAetatr;
            LINALG::Matrix<5,1> etas;
            etas.Update(1.-stab_s*Ypl/AnormEtatrial,eta_trial_vec);
            etas.Update(stab_s*apl,eta_vec,1.);
            etasAetatr.MultiplyNT(etas,AetaTr);

            DcplSymDbeta.Update(1.-Ypl/AnormEtatrial,detadbeta,1.);
            DcplSymDbeta.Multiply(cpl*Ypl/AnormEtatrial,InvPlAniso,DalphaKdbeta,1.);
            DcplSymDbeta.MultiplyNT(-1./AnormEtatrial,eta_trial_vec,dYpl_dbeta,1.);
            DcplSymDbeta.MultiplyNT((1.-stab_s)*Ypl/AnormEtatrial/AnormEtatrial,eta_trial_vec,dAnormEtaTrialDbeta,1.);
            DcplSymDbeta.MultiplyNT(stab_s*apl/AnormEtatrial,eta_vec,dAnormEtaTrialDbeta,1.);

            // **************************************************************
            // end of stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
            // **************************************************************

            // **************************************************************
            // stiffness matrix [k^e_{beta d}]_ij (i=1..5; j=1..numdof)
            // linearization of the complementarity function
            // **************************************************************
            LINALG::Matrix<numdofperelement_,1> dAnormEtatrialDd(true);
            dAnormEtatrialDd.MultiplyTN(1./AnormEtatrial,detadd,Aetatr,1.);
            LINALG::Matrix<numdofperelement_,1> dAnormEtaDd(true);
            dAnormEtaDd.MultiplyTN(1./AnormEta,detadd,Aeta,1.);

            DcplSymDd.Update(1.-Ypl/AnormEtatrial,detadd,1.);
            DcplSymDd.MultiplyNT((1.-stab_s)*Ypl/AnormEtatrial/AnormEtatrial,eta_trial_vec,dAnormEtatrialDd,1.);
            DcplSymDd.MultiplyNT(stab_s*apl/AnormEtatrial,eta_vec,dAnormEtatrialDd,1.);
            DcplSymDd.MultiplyNT(1./AnormEtatrial,eta_trial_vec,dYpl_dd,1.);

          } // active Gauss points

          else // inactive Gauss points
          {
            LINALG::Matrix<5,1> AiDalphak;
            AiDalphak.Multiply(InvPlAniso,DeltaAlphaK_vec);

            DcplSymDd.Clear();
            DcplSymDd.MultiplyNT(stab_s*cpl/Ypl,AiDalphak,dYpl_dd,1.);

            DcplSymDbeta.Multiply(cpl,InvPlAniso,DalphaKdbeta);
            DcplSymDbeta.MultiplyNT(cpl*stab_s/Ypl,AiDalphak,dYpl_dbeta,1.);

            CplSym.Multiply(cpl,InvPlAniso,DeltaAlphaK_vec);
          }
          // symmetric NCP done

          // plastic spin stuff
          // Wp-eta_sp/inityield* (-eta*delta alpha_k + delta aplpha_k * eta) = 0
          // <==> Wp = eta_sp/inityield * (eta*Dp - Dp*eta)
          // This equation is apriori skew-symmetric so it reduces to 3 independent
          // scalar-valued equations.
          SpEq(0) = -mDLp_last_iter_->at(gp)(5) -PlSpinEta/inityield *
              ( DeltaAlphaK_vec(0)*mandelstress(0,1)
                  +DeltaAlphaK_vec(2)*mandelstress(1,1)
                  +DeltaAlphaK_vec(4)*mandelstress(1,2)
                  -DeltaAlphaK_vec(2)*mandelstress(0,0)
                  -DeltaAlphaK_vec(1)*mandelstress(0,1)
                  -DeltaAlphaK_vec(3)*mandelstress(0,2) );
          SpEq(1) = -mDLp_last_iter_->at(gp)(6) -PlSpinEta/inityield *
              ( DeltaAlphaK_vec(2)*mandelstress(0,2)
                  +DeltaAlphaK_vec(1)*mandelstress(1,2)*2.
                  +DeltaAlphaK_vec(3)*mandelstress(2,2)
                  -DeltaAlphaK_vec(4)*mandelstress(0,1)
                  -DeltaAlphaK_vec(3)*mandelstress(1,1)
                  +DeltaAlphaK_vec(0)*mandelstress(1,2) );
          SpEq(2) = -mDLp_last_iter_->at(gp)(7) -PlSpinEta/inityield *
              ( DeltaAlphaK_vec(0)*mandelstress(0,2)*2.
                  +DeltaAlphaK_vec(2)*mandelstress(1,2)
                  +DeltaAlphaK_vec(4)*mandelstress(2,2)
                  -DeltaAlphaK_vec(4)*mandelstress(0,0)
                  -DeltaAlphaK_vec(3)*mandelstress(0,1)
                  +DeltaAlphaK_vec(1)*mandelstress(0,2) );

          // derivative of SpEq w.r.t. beta
          for (int i=0; i<8; i++)
          {
            dSpEqDbeta(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(0)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[1][1],i);
            dSpEqDbeta(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[1][2],i);
            dSpEqDbeta(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[0][0],i);
            dSpEqDbeta(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(1)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[0][2],i);

            dSpEqDbeta(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[0][2],i);
            dSpEqDbeta(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dsigmadbeta(VOIGT3X3SYM_[1][2],i)*2.;
            dSpEqDbeta(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[2][2],i);
            dSpEqDbeta(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[1][1],i);
            dSpEqDbeta(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dsigmadbeta(VOIGT3X3SYM_[1][2],i);

            dSpEqDbeta(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dsigmadbeta(VOIGT3X3SYM_[0][2],i)*2.;
            dSpEqDbeta(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[1][2],i);
            dSpEqDbeta(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[2][2],i);
            dSpEqDbeta(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[0][0],i);
            dSpEqDbeta(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dsigmadbeta(VOIGT3X3SYM_[0][2],i);
          }
          dSpEqDbeta(0,5) -=1.;
          dSpEqDbeta(0,0) -=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(0,2) -=PlSpinEta/inityield*mandelstress(1,1);
          dSpEqDbeta(0,4) -=PlSpinEta/inityield*mandelstress(1,2);
          dSpEqDbeta(0,2) +=PlSpinEta/inityield*mandelstress(0,0);
          dSpEqDbeta(0,1) +=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(0,3) +=PlSpinEta/inityield*mandelstress(0,2);

          dSpEqDbeta(1,6) -=1.;
          dSpEqDbeta(1,2) -=PlSpinEta/inityield*mandelstress(0,2);
          dSpEqDbeta(1,1) -=PlSpinEta/inityield*mandelstress(1,2)*2.;
          dSpEqDbeta(1,3) -=PlSpinEta/inityield*mandelstress(2,2);
          dSpEqDbeta(1,4) +=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(1,3) +=PlSpinEta/inityield*mandelstress(1,1);
          dSpEqDbeta(1,0) +=PlSpinEta/inityield*mandelstress(1,2);

          dSpEqDbeta(2,7) -=1.;
          dSpEqDbeta(2,0) -=PlSpinEta/inityield*mandelstress(0,2)*2.;
          dSpEqDbeta(2,2) -=PlSpinEta/inityield*mandelstress(1,2);
          dSpEqDbeta(2,4) -=PlSpinEta/inityield*mandelstress(2,2);
          dSpEqDbeta(2,4) +=PlSpinEta/inityield*mandelstress(0,0);
          dSpEqDbeta(2,3) +=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(2,1) -=PlSpinEta/inityield*mandelstress(0,2);

          // derivative of SpEq w.r.t. displacements
          for (int i=0; i<numdofperelement_; i++)
          {
            dSpEqDd(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(0)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[1][1],i);
            dSpEqDd(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[1][2],i);
            dSpEqDd(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[0][0],i);
            dSpEqDd(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(1)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[0][2],i);

            dSpEqDd(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[0][2],i);
            dSpEqDd(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dSigmadd(VOIGT3X3SYM_[1][2],i)*2.;
            dSpEqDd(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[2][2],i);
            dSpEqDd(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[1][1],i);
            dSpEqDd(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dSigmadd(VOIGT3X3SYM_[1][2],i);

            dSpEqDd(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dSigmadd(VOIGT3X3SYM_[0][2],i)*2.;
            dSpEqDd(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[1][2],i);
            dSpEqDd(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[2][2],i);
            dSpEqDd(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[0][0],i);
            dSpEqDd(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dSigmadd(VOIGT3X3SYM_[0][2],i);
          }

//          // FD check beta
//          double epsilon=1.e-10;
//          for (int i=0; i<8; i++)
//          {
//            mDLp_last_iter_->at(gp)(i) += epsilon;
//            LINALG::Matrix<3,3> mDLp_new(false);
//            LINALG::Matrix<3,3> DeltaAlphaK_new(false);
//            LINALG::Matrix<5,1> DeltaAlphaK_vec_new(false);
//            for (int iii=0; iii<5; iii++) DeltaAlphaK_vec_new(iii) = mDLp_last_iter_->at(gp)(iii);
//            // -Dp
//            DeltaAlphaK_new(0,0) = mDLp_last_iter_->at(gp)(0);
//            DeltaAlphaK_new(1,1) = mDLp_last_iter_->at(gp)(1);
//            DeltaAlphaK_new(2,2) = -1.0*(mDLp_last_iter_->at(gp)(0)+mDLp_last_iter_->at(gp)(1));
//            DeltaAlphaK_new(0,1) = mDLp_last_iter_->at(gp)(2);
//            DeltaAlphaK_new(1,0) = mDLp_last_iter_->at(gp)(2);
//            DeltaAlphaK_new(1,2) = mDLp_last_iter_->at(gp)(3);
//            DeltaAlphaK_new(2,1) = mDLp_last_iter_->at(gp)(3);
//            DeltaAlphaK_new(0,2) = mDLp_last_iter_->at(gp)(4);
//            DeltaAlphaK_new(2,0) = mDLp_last_iter_->at(gp)(4);
//            mDLp_new.Update(DeltaAlphaK_new);
//            // -Wp
//            mDLp_new(0,1) += mDLp_last_iter_->at(gp)(5);
//            mDLp_new(1,0) -= mDLp_last_iter_->at(gp)(5);
//            mDLp_new(1,2) += mDLp_last_iter_->at(gp)(6);
//            mDLp_new(2,1) -= mDLp_last_iter_->at(gp)(6);
//            mDLp_new(0,2) += mDLp_last_iter_->at(gp)(7);
//            mDLp_new(2,0) -= mDLp_last_iter_->at(gp)(7);
//            double absDeltaAlphaK_new=DeltaAlphaK_new.Norm2();
//            tmp1=mDLp_new;
//            MatrixExponential3x3(tmp1);
//
//            LINALG::Matrix<3,3> InvPlasticDefgrd_new;
//            InvPlasticDefgrd_new.Multiply(InvPlasticDefgrdLast,tmp1);
////            for (int i=0; i<3; i++)
////              for (int j=0; j<3; j++)
////              {
////                LINALG::Matrix<3,3> InvPlasticDefgrd_new(InvPlasticDefgrd);
////                InvPlasticDefgrd_new(i,j) +=epsilon;
//
//            // material call
//            LINALG::Matrix<numstr_,1> pk2_stress_new(true);
//            LINALG::Matrix<6,6> Cmat_ABCD_new(true);
//            LINALG::Matrix<6,9> dpk2dfpinv_new(true);
//            LINALG::Matrix<3,3> mandelstress_new(true);
//            LINALG::Matrix<6,6> dmdc_new(true);
//            LINALG::Matrix<6,9> dmdfpinv_new(true);
//            params.set<int>("gp",gp);
//            params.set<int>("eleID",Id());
//            plmat->Evaluate(&defgrd,&InvPlasticDefgrd_new,params,&pk2_stress_new,&Cmat_ABCD_new,&dpk2dfpinv_new,&mandelstress_new,&dmdc_new,&dmdfpinv_new);
//
//            // equivalent stress eta
//            LINALG::Matrix<nsd_,nsd_> eta_new(mandelstress_new);
//            for (int ii=0; ii<nsd_; ii++)
//              eta_new(ii,ii) -= 1./3.*(mandelstress_new(0,0) + mandelstress_new(1,1) + mandelstress_new(2,2));
//            eta_new.Update(2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
//            eta_new.Update(2./3.*kinhard,DeltaAlphaK_new,1.);
//            LINALG::Matrix<5,1> eta_vec_new(false);
//            eta_vec_new(0) = eta_new(0,0);
//            eta_vec_new(1) = eta_new(1,1);
//            eta_vec_new(2) = 0.5*(eta_new(0,1)+eta_new(1,0));
//            eta_vec_new(3) = 0.5*(eta_new(2,1)+eta_new(1,2));
//            eta_vec_new(4) = 0.5*(eta_new(0,2)+eta_new(2,0));
//
//            LINALG::Matrix<5,1> Aeta_new;
//            Aeta_new.Multiply(PlAniso,eta_vec_new);
//            tmp1(0,0) = 2./3.*Aeta_new(0)-1./3.*Aeta_new(1);
//            tmp1(1,1) =-1./3.*Aeta_new(0)+2./3.*Aeta_new(1);
//            tmp1(2,2) = -tmp1(0,0)-tmp1(1,1);
//            tmp1(0,1) = Aeta_new(2)/2.;
//            tmp1(1,0) = Aeta_new(2)/2.;
//            tmp1(1,2) = Aeta_new(3)/2.;
//            tmp1(2,1) = Aeta_new(3)/2.;
//            tmp1(0,2) = Aeta_new(4)/2.;
//            tmp1(2,0) = Aeta_new(4)/2.;
//            double EucNormAeta_new = tmp1.Norm2();
//            double DpAe_new=0.;
//            for (int ii=0; ii<3; ii++)
//              for (int jj=0; jj<3; jj++)
//                DpAe_new+=DeltaAlphaK_new(ii,jj)*tmp1(ii,jj);
//
//            // eta_trial
//            LINALG::Matrix<5,1> eta_trial_vec_new(eta_vec_new);
//            eta_trial_vec_new.Multiply(-cpl,InvPlAniso,DeltaAlphaK_vec,1.);
//
//            double AnormEtatrial_new=0.; // sqrt( eta_tr : A : eta_tr )
//            LINALG::Matrix<1,1> tmp11;
//            LINALG::Matrix<5,1> tmp51;
//            tmp51.Multiply(PlAniso,eta_trial_vec_new);
//            tmp11.MultiplyTN(eta_trial_vec_new,tmp51);
//            AnormEtatrial_new=sqrt(tmp11(0,0));
//            double AnormEta_new=0.; // sqrt( eta_tr : A : eta_tr )
//            tmp51.Multiply(PlAniso,eta_vec_new);
//            tmp11.MultiplyTN(eta_vec_new,tmp51);
//            AnormEta_new=sqrt(tmp11(0,0));
//
//            double deltaAlphaI_new=0.;
//            if (DpAe<0.)
//              deltaAlphaI_new = -sqrt(2./3.)*DpAe_new*AnormEta_new/EucNormAeta_new/EucNormAeta_new;
//
//            double Ypl_new = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_new)))
//                + isohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_new) +inityield);
//
//
////            // eta_trial
////            LINALG::Matrix<5,1> eta_trial_vec_new(eta_vec_new);
////            eta_trial_vec_new.Multiply(-cpl,InvPlAniso,DeltaAlphaK_vec_new,1.);
////
////            double AnormEtatrial_new=0.; // sqrt( eta_tr : A : eta_tr )
////            LINALG::Matrix<1,1> tmp11;
////            LINALG::Matrix<5,1> tmp51;
////            tmp51.Multiply(PlAniso,eta_trial_vec_new);
////            tmp11.MultiplyTN(eta_trial_vec_new,tmp51);
////            AnormEtatrial_new=sqrt(tmp11(0,0));
////
////            double Ypl_new=0.;
////            if (Dissipation>0.)
////              Ypl_new = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+ sqrt(2./3.)*absDeltaAlphaK_new)))
////                  + isohard*(last_alpha_isotropic_->at(gp)(0,0)+ sqrt(2./3.)*absDeltaAlphaK_new) +inityield);
////            else
////              Ypl_new = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0))))
////                  + isohard*(last_alpha_isotropic_->at(gp)(0,0)) +inityield);
//
//
//            LINALG::Matrix<5,1> CplSym_new(true);
//            CplSym_new.Update(1.,eta_vec_new,0.);
//            CplSym_new.Update(-Ypl_new/AnormEtatrial_new,eta_trial_vec_new,1.);
//
//            LINALG::Matrix<numdofperelement_,1> FintGp_new;
//            FintGp_new.MultiplyTN(detJ_w, bop, pk2_stress_new, 0.0);
////
//            LINALG::Matrix<3,1> SpEq_new;
//            SpEq_new(0) = -mDLp_last_iter_->at(gp)(5) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec_new(0)*eta_vec_new(2)
//                    +DeltaAlphaK_vec_new(2)*eta_vec_new(1)
//                    +DeltaAlphaK_vec_new(4)*eta_vec_new(3)
//                    -DeltaAlphaK_vec_new(2)*eta_vec_new(0)
//                    -DeltaAlphaK_vec_new(1)*eta_vec_new(2)
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(4) );
//            SpEq_new(1) = -mDLp_last_iter_->at(gp)(6) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec_new(2)*eta_vec_new(4)
//                    +DeltaAlphaK_vec_new(1)*eta_vec_new(3)*2.
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(0)
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(1)*2.
//                    -DeltaAlphaK_vec_new(4)*eta_vec_new(2)
//                    +DeltaAlphaK_vec_new(0)*eta_vec_new(3) );
//            SpEq_new(2) = -mDLp_last_iter_->at(gp)(7) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec_new(0)*eta_vec_new(4)*2.
//                    +DeltaAlphaK_vec_new(2)*eta_vec_new(3)
//                    -DeltaAlphaK_vec_new(4)*eta_vec_new(0)*2.
//                    -DeltaAlphaK_vec_new(4)*eta_vec_new(1)
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(2)
//                    +DeltaAlphaK_vec_new(1)*eta_vec_new(4) );
//
//            std::cout << std::scientific;
//            std::cout.precision(8);
//            for (int a=0; a<5; a++)
////              for (int A=0; A<3; A++)
//              {
//                double ref=DcplSymDbeta(a,i);
//                double FD=(CplSym_new(a)-CplSym(a))/epsilon;
//                double relE=abs((FD-ref)/ref);
//                std::cout << "Aai:"<<i<< "\tFD: " << FD << "\tref: " << ref
//                    << "\trelE: " << relE << std::endl;
//              }
//            mDLp_last_iter_->at(gp)(i) -= epsilon;
//          }
//              std::cout << std::endl;
////dserror("stop");

//          // FD check d
//          for (int FDdisp=0; FDdisp<numdofperelement_; FDdisp++)
//          {
//            double epsilon = 1.e-10;
//            disp.at(FDdisp) += epsilon;
//            // calculate stuff with disturbed displacement
//            for (int i=0; i<nen_; ++i)
//            {
//              const double* x = nodes[i]->X();
//              xrefe(i,0) = x[0];
//              xrefe(i,1) = x[1];
//              xrefe(i,2) = x[2];
//
//              xcurr(i,0) = xrefe(i,0) + disp[i*numdofpernode_+0];
//              xcurr(i,1) = xrefe(i,1) + disp[i*numdofpernode_+1];
//              xcurr(i,2) = xrefe(i,2) + disp[i*numdofpernode_+2];
//            }
//            // compute derivatives N_XYZ at gp w.r.t. material coordinates
//            // by N_XYZ = J^-1 * N_rst
//            N_XYZ.Multiply(invJ_[gp],deriv);
//            detJ = detJ_[gp];
//
//            // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
//            LINALG::Matrix<3,3> defgrd_new;
//            defgrd_new.MultiplyTT(xcurr,N_XYZ);
//
//            // material call
//            LINALG::Matrix<numstr_,1> pk2_stress_new(true);
//            LINALG::Matrix<6,6> Cmat_ABCD_new(true);
//            LINALG::Matrix<6,9> dpk2dfpinv_new(true);
//            LINALG::Matrix<3,3> mandelstress_new(true);
//            LINALG::Matrix<6,6> dmdc_new(true);
//            LINALG::Matrix<6,9> dmdfpinv_new(true);
//            params.set<int>("gp",gp);
//            params.set<int>("eleID",Id());
//            plmat->Evaluate(&defgrd_new,&InvPlasticDefgrd,params,&pk2_stress_new,&Cmat_ABCD_new,&dpk2dfpinv_new,&mandelstress_new,&dmdc_new,&dmdfpinv_new);
//
//            // equivalent stress eta
//            LINALG::Matrix<nsd_,nsd_> eta_new(mandelstress_new);
//            for (int ii=0; ii<nsd_; ii++)
//              eta_new(ii,ii) -= 1./3.*(mandelstress_new(0,0) + mandelstress_new(1,1) + mandelstress_new(2,2));
//            eta_new.Update(2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
//            eta_new.Update(2./3.*kinhard,DeltaAlphaK,1.);
//            LINALG::Matrix<5,1> eta_vec_new(false);
//            eta_vec_new(0) = eta_new(0,0);
//            eta_vec_new(1) = eta_new(1,1);
//            eta_vec_new(2) = 0.5*(eta_new(0,1)+eta_new(1,0));
//            eta_vec_new(3) = 0.5*(eta_new(2,1)+eta_new(1,2));
//            eta_vec_new(4) = 0.5*(eta_new(0,2)+eta_new(2,0));
//
//
//            LINALG::Matrix<5,1> Aeta_new;
//            Aeta_new.Multiply(PlAniso,eta_vec_new);
//            tmp1(0,0) = 2./3.*Aeta_new(0)-1./3.*Aeta_new(1);
//            tmp1(1,1) =-1./3.*Aeta_new(0)+2./3.*Aeta_new(1);
//            tmp1(2,2) = -tmp1(0,0)-tmp1(1,1);
//            tmp1(0,1) = Aeta_new(2)/2.;
//            tmp1(1,0) = Aeta_new(2)/2.;
//            tmp1(1,2) = Aeta_new(3)/2.;
//            tmp1(2,1) = Aeta_new(3)/2.;
//            tmp1(0,2) = Aeta_new(4)/2.;
//            tmp1(2,0) = Aeta_new(4)/2.;
//            double EucNormAeta_new = tmp1.Norm2();
//            double DpAe_new=0.;
//            for (int ii=0; ii<3; ii++)
//              for (int jj=0; jj<3; jj++)
//                DpAe_new+=DeltaAlphaK(ii,jj)*tmp1(ii,jj);
//
//            // eta_trial
//            LINALG::Matrix<5,1> eta_trial_vec_new(eta_vec_new);
//            eta_trial_vec_new.Multiply(-cpl,InvPlAniso,DeltaAlphaK_vec,1.);
//
//            double AnormEtatrial_new=0.; // sqrt( eta_tr : A : eta_tr )
//            LINALG::Matrix<1,1> tmp11;
//            LINALG::Matrix<5,1> tmp51;
//            tmp51.Multiply(PlAniso,eta_trial_vec_new);
//            tmp11.MultiplyTN(eta_trial_vec_new,tmp51);
//            AnormEtatrial_new=sqrt(tmp11(0,0));
//            double AnormEta_new=0.; // sqrt( eta_tr : A : eta_tr )
//            tmp51.Multiply(PlAniso,eta_vec_new);
//            tmp11.MultiplyTN(eta_vec_new,tmp51);
//            AnormEta_new=sqrt(tmp11(0,0));
//
//            double deltaAlphaI_new=0.;
//            if (DpAe<0.)
//              deltaAlphaI_new = -sqrt(2./3.)*DpAe_new*AnormEta_new/EucNormAeta_new/EucNormAeta_new;
//
//            double Ypl_new = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_new)))
//                + isohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_new) +inityield);
//
//
//            LINALG::Matrix<5,1> CplSym_new(true);
//            CplSym_new.Update(1.-Ypl/AnormEtatrial_new,Aeta_new,1.);
//            CplSym_new.Update(Ypl/AnormEtatrial_new*cpl,DeltaAlphaK_vec,1.);
//            CplSym_new.Scale(pow(AnormEtatrial_new,stab_s));
//
//            LINALG::Matrix<3,1> SpEq_new;
//            SpEq_new(0) = -mDLp_last_iter_->at(gp)(5) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec(0)*eta_vec_new(2)
//                    +DeltaAlphaK_vec(2)*eta_vec_new(1)
//                    +DeltaAlphaK_vec(4)*eta_vec_new(3)
//                    -DeltaAlphaK_vec(2)*eta_vec_new(0)
//                    -DeltaAlphaK_vec(1)*eta_vec_new(2)
//                    -DeltaAlphaK_vec(3)*eta_vec_new(4) );
//            SpEq_new(1) = -mDLp_last_iter_->at(gp)(6) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec(2)*eta_vec_new(4)
//                    +DeltaAlphaK_vec(1)*eta_vec_new(3)*2.
//                    -DeltaAlphaK_vec(3)*eta_vec_new(0)
//                    -DeltaAlphaK_vec(3)*eta_vec_new(1)*2.
//                    -DeltaAlphaK_vec(4)*eta_vec_new(2)
//                    +DeltaAlphaK_vec(0)*eta_vec_new(3) );
//            SpEq_new(2) = -mDLp_last_iter_->at(gp)(7) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec(0)*eta_vec_new(4)*2.
//                    +DeltaAlphaK_vec(2)*eta_vec_new(3)
//                    -DeltaAlphaK_vec(4)*eta_vec_new(0)*2.
//                    -DeltaAlphaK_vec(4)*eta_vec_new(1)
//                    -DeltaAlphaK_vec(3)*eta_vec_new(2)
//                    +DeltaAlphaK_vec(1)*eta_vec_new(4) );
//
//            std::cout << std::scientific;
//            std::cout.precision(8);
////            for (int a=0; a<5; a++)
//              //                for (int A=0; A<3; A++)
//              {
//                double ref=dYpl_dd(FDdisp);
//                double FD=(Ypl_new-Ypl)/epsilon;
//                double relE=abs((FD-ref)/ref);
//                std::cout << "ai:"<<FDdisp << "\tFD: " << FD << "\tref: " << ref
//                    << "\trelE: " << relE << std::endl;
//              }
//            // undo FD
//            disp.at(FDdisp) -= epsilon;
//            for (int i=0; i<nen_; ++i)
//            {
//              const double* x = nodes[i]->X();
//              xrefe(i,0) = x[0];
//              xrefe(i,1) = x[1];
//              xrefe(i,2) = x[2];
//
//              xcurr(i,0) = xrefe(i,0) + disp[i*numdofpernode_+0];
//              xcurr(i,1) = xrefe(i,1) + disp[i*numdofpernode_+1];
//              xcurr(i,2) = xrefe(i,2) + disp[i*numdofpernode_+2];
//            }
//          }
//          std::cout << std::endl;

          // build matrices kbb and kbd
          for (int i=0; i<5; i++)
            for (int j=0; j<numdofperelement_; j++)
              kbetad(i,j) = DcplSymDd(i,j);
          for (int i=0; i<3; i++)
            for (int j=0; j<numdofperelement_; j++)
              kbetad(i+5,j) = dSpEqDd(i,j);

          for (int i=0; i<5; i++)
            for (int j=0; j<8; j++)
              kbetabeta(i,j) = DcplSymDbeta(i,j);
          for (int i=0; i<3; i++)
            for (int j=0; j<8; j++)
              kbetabeta(i+5,j) = dSpEqDbeta(i,j);

          for (int i=0; i<5; i++)
            force_beta(i) = CplSym(i);
          for (int i=0; i<3; i++)
            force_beta(i+5) = SpEq(i);

        // **************************************************************
        // static condensation of inner variables
        // **************************************************************
        //inverse matrix block [k_beta beta]_ij
        LINALG::Matrix<8,8> InvKbetabeta;
        Epetra_SerialDenseMatrix Kbetabeta_epetra(8,8);
        for (int i=0; i<8; i++)
          for (int j=0; j<8; j++)
            Kbetabeta_epetra(i,j) = kbetabeta(i,j);
        // we need the inverse of K_beta beta
        Epetra_SerialDenseSolver solve_for_inverseKbb;
        solve_for_inverseKbb.SetMatrix(Kbetabeta_epetra);
        solve_for_inverseKbb.Invert();
        for (int i=0; i<8; i++)
          for (int j=0; j<8; j++)
            InvKbetabeta(i,j) = Kbetabeta_epetra(i,j);
        // store for recover step
        KbbInvHill_->at(gp) = InvKbetabeta;
        fbetaHill_->at(gp) = force_beta;
        KbdHill_->at(gp) = kbetad;

        LINALG::Matrix<numdofperelement_,8> KdbKbb; // temporary  Kdb.Kbb^-1
        KdbKbb.Multiply(kdbeta,InvKbetabeta);

        // "plastic displacement stiffness"
        // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
        if (stiffmatrix!=NULL) stiffmatrix->Multiply(-1.,KdbKbb,kbetad,1.);
        // "plastic internal force"
        // plFint = [K_db.K_bb^-1].f_b
        if (force!=NULL) force->Multiply(-1.,KdbKbb,force_beta,1.);

        LINALG::Matrix<numdofperelement_,1> force_update;
        force_update.Multiply(-1.,KdbKbb,force_beta,0.);
        LINALG::Matrix<numdofperelement_,numdofperelement_> stiff_update;
        stiff_update.Multiply(-1.,KdbKbb,kbetad,0.);
        }
        else if (mDLp_last_iter_->at(gp).NormInf()!=0.)
        {
          KbbInvHill_->at(gp).Clear();
          for (int i=0; i<8; i++) KbbInvHill_->at(gp)(i,i)=1.;
          KbdHill_->at(gp).Clear();
          fbetaHill_->at(gp).Update(mDLp_last_iter_->at(gp));

          // condensation to internal force vector
          if (force!=NULL)
            force->Multiply(-1.,kdbeta,fbetaHill_->at(gp),1.);
        }
      } // modification for plastic Gauss points
      else
      {
        KbbInvHill_->at(gp).Clear();
        fbetaHill_ ->at(gp).Clear();
        KbdHill_   ->at(gp).Clear();
      }

      // square of the residual L2 norm
      double residual_norm_sqare=fbetaHill_->at(gp)(0)*fbetaHill_->at(gp)(0)
                                +fbetaHill_->at(gp)(1)*fbetaHill_->at(gp)(1)
                                +(-fbetaHill_->at(gp)(0)-fbetaHill_->at(gp)(1))*(-fbetaHill_->at(gp)(0)-fbetaHill_->at(gp)(1));
      for (int i=2; i<8; i++) residual_norm_sqare+=fbetaHill_->at(gp)(i)*fbetaHill_->at(gp)(i)*2.;
      params.get<double>("Lp_residual_square")+=residual_norm_sqare;

  } // gp loop

  // communicate unconverged active set to time integration
  if (converged_active_set==false)
    params.set("unconverged_active_set",true);

  return;
}

/*----------------------------------------------------------------------*
 | internal force, stiffness and mass for f-bar elements    seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::nln_stiffmass_fbar(
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
  // check for f-bar element
  DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * eleFBAR
  = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* >(this);
  if (!eleFBAR)
    dserror("invalid call for nln_stiffmass_fbar from non-fbar element");

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
  LINALG::Matrix<3,nen_> N_XYZ;
  // build deformation gradient wrt to material configuration
  LINALG::Matrix<3,3> defgrd(false);
  // shape functions and their first derivatives
  LINALG::Matrix<nen_,1> shapefunct;
  LINALG::Matrix<3,nen_> deriv;

  // ---------------------- deformation gradient at centroid of element
  double detF_0 = -1.0;
  LINALG::Matrix<3,3> invdefgrd_0(false);
  LINALG::Matrix<3,nen_> N_XYZ_0(false);
  //element coordinate derivatives at centroid
  LINALG::Matrix<3,nen_> N_rst_0(false);
  DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, DRT::Element::hex8);

  //inverse jacobian matrix at centroid
  LINALG::Matrix<3,3> invJ_0(false);
  invJ_0.Multiply(N_rst_0,xrefe);
  invJ_0.Invert();
  //material derivatives at centroid
  N_XYZ_0.Multiply(invJ_0,N_rst_0);

  //deformation gradient and its determinant at centroid
  LINALG::Matrix<3,3> defgrd_0(false);
  defgrd_0.MultiplyTT(xcurr,N_XYZ_0);
  invdefgrd_0.Invert(defgrd_0);
  detF_0 = defgrd_0.Determinant();

  // 3x3 LINALG matrix for temporary stuff
  LINALG::Matrix<3,3> tmp1(false);
  LINALG::Matrix<3,3> tmp2(false);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // get plastic material parameters
  double kinhard=plmat->Kinhard();
  double isohard=plmat->Isohard();
  double expisohard=plmat->Expisohard();
  double infyield=plmat->Infyield();
  double inityield=plmat->Inityield();
  double cpl=plparams_->get<double>("SEMI_SMOOTH_CPL");
  double stab_s=plparams_->get<double>("STABILIZATION_S");

  // converged active set
  bool converged_active_set=true;

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
    LINALG::Matrix<3,3> invdefgrd;
    double detF=invdefgrd.Invert(defgrd);

    // check for element distortion
    if (detF<=0. || detF_0<=0.)
      dserror("element distortion too large");

    // F_bar deformation gradient =(detF_0/detF)^1/3*F
    LINALG::Matrix<3,3> defgrd_bar(defgrd);
    double f_bar_factor=pow(detF_0/detF,1/3.);
    defgrd_bar.Scale(f_bar_factor);

    // calcualte total rcg/gl/ea for output
    LINALG::Matrix<3,3> total_cauchygreen(false);
    total_cauchygreen.MultiplyTN(defgrd,defgrd);
    // total Cauchy green in voigt notation
    LINALG::Matrix<6,1> RCG;
    for (int i=0; i<3; i++) RCG(i)=total_cauchygreen(i,i);
    RCG(3)=total_cauchygreen(0,1)*2.;
    RCG(4)=total_cauchygreen(1,2)*2.;
    RCG(5)=total_cauchygreen(0,2)*2.;

    // calculate nonlinear B-operator
    LINALG::Matrix<numstr_,numdofperelement_> bop(false);
    CalculateBop(&bop,&defgrd,&N_XYZ);

    // strain output *********************************
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
    // end of strain output **************************

    // recover condensed variables from last iteration step *********
    // first part
    LINALG::Matrix<5,1> tmp51;
    tmp51.Multiply(KbbInv_->at(gp),fbeta_->at(gp));

    // second part
    LINALG::Matrix<5,numdofperelement_> tmp524;
    tmp524.Multiply(KbbInv_->at(gp),Kbd_->at(gp));
    tmp51.Multiply(1.,tmp524,res_d,1.);
    DalphaK_last_iter_->at(gp).Update(-1.,tmp51,1.);
    double increment_norm_sqare=tmp51(0)*tmp51(0)
                              +tmp51(1)*tmp51(1)
                              +(-tmp51(0)-tmp51(1))*(-tmp51(0)-tmp51(1))
                              +tmp51(2)*tmp51(2)*2.
                              +tmp51(3)*tmp51(3)*2.
                              +tmp51(4)*tmp51(4)*2.;
    params.get<double>("Lp_increment_square")+=increment_norm_sqare;
    // end of recover **********************************************

    // current plastic flow increment
    LINALG::Matrix<3,3> DeltaAlphaK(false);
    DeltaAlphaK(0,0) = DalphaK_last_iter_->at(gp)(0);
    DeltaAlphaK(1,1) = DalphaK_last_iter_->at(gp)(1);
    DeltaAlphaK(2,2) = -1.0*(DalphaK_last_iter_->at(gp)(0)+DalphaK_last_iter_->at(gp)(1));
    DeltaAlphaK(0,1) = DalphaK_last_iter_->at(gp)(2);
    DeltaAlphaK(1,0) = DalphaK_last_iter_->at(gp)(2);
    DeltaAlphaK(1,2) = DalphaK_last_iter_->at(gp)(3);
    DeltaAlphaK(2,1) = DalphaK_last_iter_->at(gp)(3);
    DeltaAlphaK(0,2) = DalphaK_last_iter_->at(gp)(4);
    DeltaAlphaK(2,0) = DalphaK_last_iter_->at(gp)(4);

    // inverse plastic deformation gradient
    LINALG::Matrix<3,3> InvPlasticDefgrd(false);
    // inverse plastic deformation gradient at last time step
    LINALG::Matrix<3,3> InvPlasticDefgrdLast = last_plastic_defgrd_inverse_->at(gp);
    tmp1=DeltaAlphaK;
    MatrixExponential3x3(tmp1);
    InvPlasticDefgrd.Multiply(InvPlasticDefgrdLast,tmp1);

    // material call
    LINALG::Matrix<numstr_,1> pk2_stress_bar(true);
    LINALG::Matrix<6,6> Cmat_ABCD_bar(true);
    LINALG::Matrix<6,9> dpk2dfpinv(true);
    LINALG::Matrix<3,3> mandelstress(true);
    LINALG::Matrix<6,6> dmdc(true);
    LINALG::Matrix<6,9> dmdfpinv(true);
    params.set<int>("gp",gp);
    params.set<int>("eleID",Id());
    plmat->Evaluate(&defgrd_bar,&InvPlasticDefgrd,params,&pk2_stress_bar,&Cmat_ABCD_bar,&dpk2dfpinv,&mandelstress,&dmdc,&dmdfpinv);

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < numstr_; ++i)
        (*elestress)(gp,i) = pk2_stress_bar(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF = defgrd.Determinant();

      LINALG::Matrix<3,3> pkstress;
      pkstress(0,0) = pk2_stress_bar(0);
      pkstress(0,1) = pk2_stress_bar(3);
      pkstress(0,2) = pk2_stress_bar(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = pk2_stress_bar(1);
      pkstress(1,2) = pk2_stress_bar(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = pk2_stress_bar(2);

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
    // equivalent stress eta
    LINALG::Matrix<3,3> eta_bar(mandelstress);
    for (int i=0; i<3; i++)
      eta_bar(i,i) -= 1./3.*(mandelstress(0,0) + mandelstress(1,1) + mandelstress(2,2));
    eta_bar.Update(2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
    eta_bar.Update(2./3.*kinhard,DeltaAlphaK,1.);
    LINALG::Matrix<5,1> eta_bar_vec(false);
    eta_bar_vec(0) = eta_bar(0,0);
    eta_bar_vec(1) = eta_bar(1,1);
    eta_bar_vec(2) = 0.5*(eta_bar(0,1)+eta_bar(1,0));
    eta_bar_vec(3) = 0.5*(eta_bar(2,1)+eta_bar(1,2));
    eta_bar_vec(4) = 0.5*(eta_bar(0,2)+eta_bar(2,0));

    // eta_trial
    LINALG::Matrix<3,3> eta_bar_trial(eta_bar);
    eta_bar_trial.Update(-1.*cpl,DeltaAlphaK,1.);
    LINALG::Matrix<5,1> eta_bar_trial_vec(false);
    eta_bar_trial_vec(0) = eta_bar_trial(0,0);
    eta_bar_trial_vec(1) = eta_bar_trial(1,1);
    eta_bar_trial_vec(2) = 0.5*(eta_bar_trial(0,1)+eta_bar_trial(1,0));
    eta_bar_trial_vec(3) = 0.5*(eta_bar_trial(2,1)+eta_bar_trial(1,2));
    eta_bar_trial_vec(4) = 0.5*(eta_bar_trial(0,2)+eta_bar_trial(2,0));

    // || eta_trial || = sqrt(eta_trial_ij * eta_trial_ij)
    double absetatrial_bar=eta_bar_trial.Norm2();
    double abseta_bar=eta_bar.Norm2();
    double Dissipation=0.;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        Dissipation-=eta_bar(i,j)*DeltaAlphaK(i,j);

    double deltaAlphaI=0.;
    if (Dissipation>0. && abseta_bar!=0.)
      deltaAlphaI=Dissipation/abseta_bar*sqrt(2./3.);

    // current yield stress equivalent (yield stress scaled by sqrt(2/3))
    double Ypl=0.;
    if (Dissipation>0.)
      Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+ deltaAlphaI)))
          + isohard*(last_alpha_isotropic_->at(gp)(0,0)+ deltaAlphaI) +inityield);
    else
      Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0))))
          + isohard*(last_alpha_isotropic_->at(gp)(0,0)) +inityield);

    // check activity state
    // inactive
    if (Ypl<absetatrial_bar)
    {
      if (activity_state_->at(gp)==false) // gp switches state
      {
        if (abs(Ypl-absetatrial_bar)>AS_CONVERGENCE_TOL*inityield
            || DalphaK_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
          converged_active_set = false;
      }
      activity_state_->at(gp) = true;
    }
    // active
    else
    {
      if (activity_state_->at(gp)==true) // gp switches state
      {
        if (abs(Ypl-absetatrial_bar)>AS_CONVERGENCE_TOL*inityield
            || DalphaK_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
          converged_active_set = false;
      }
      activity_state_->at(gp) = false;
    }

    // integrate usual internal force and stiffness matrix
    double detJ_w = detJ*intpoints_.Weight(gp);
    // integrate elastic internal force vector **************************
    // update internal force vector
    if (force != NULL)
    {
      force->MultiplyTN(detJ_w/f_bar_factor, bop, pk2_stress_bar, 1.0);
    }

    // additional f-bar derivatives
    LINALG::Matrix<numdofperelement_,1> htensor(true);

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      for(int n=0;n<numdofperelement_;n++)
        for(int i=0;i<3;i++)
          htensor(n) += invdefgrd_0(i,n%3)*N_XYZ_0(i,n/3)-invdefgrd(i,n%3)*N_XYZ(i,n/3);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,numdofperelement_> cb;
      cb.Multiply(Cmat_ABCD_bar,bop);
      stiffmatrix->MultiplyTN(detJ_w*f_bar_factor,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(pk2_stress_bar); // auxiliary integrated stress
      sfac.Scale(detJ_w/f_bar_factor); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
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
      LINALG::Matrix<numstr_,1> ccg;
      ccg.Multiply(Cmat_ABCD_bar,RCG);

      LINALG::Matrix<numdofperelement_,1> bopccg(false); // auxiliary integrated stress
      bopccg.MultiplyTN(detJ_w*f_bar_factor/3.0,bop,ccg);

      LINALG::Matrix<numdofperelement_,1> bops(false); // auxiliary integrated stress
      bops.MultiplyTN(-detJ_w/f_bar_factor/3.0,bop,pk2_stress_bar);
      stiffmatrix->MultiplyNT(1.,bops,htensor,1.);
      stiffmatrix->MultiplyNT(1.,bopccg,htensor,1.);
      // end of integrate additional fbar matrix*****************************
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
    if (
        (stiffmatrix!=NULL || force!=NULL)
        &&
        (activity_state_->at(gp)==true || DalphaK_last_iter_->at(gp).NormInf()!=0.)
       )
    {
      // variables needed for condensation and calculated seperately for active and inactive Gauss points
      LINALG::Matrix<numdofperelement_,5> kdbeta;

      // 5x5 identitiy matrix
      LINALG::Matrix<5,5> id5(true);
      for (int i=0; i<5; i++)
        id5(i,i)=1.;

        // damping parameter apl
        double apl=1.;
        if (Ypl/abseta_bar<1.) apl=Ypl/abseta_bar;

        // derivative of the matrix exponential
        // for ease of notation, we don't do the derivative w.r.t. DeltaLp
        // but -DeltaLp=DeltaAlphaK
        LINALG::Matrix<6,6> Dexp(false);
        MatrixExponentialDerivativeSym3x3(DeltaAlphaK,Dexp);

        LINALG::Matrix<9,5> DFpiDbeta(true);
        for (int a=0; a<3; a++)
          for (int A=0; A<3; A++)
            for (int j=0; j<5; j++)
              for (int b=0; b<3; b++)
              {
                if (j==0 || j==1)
                  DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) +=
                      (Dexp(VOIGT3X3SYM_[a][b],j)-Dexp(VOIGT3X3SYM_[a][b],2))*InvPlasticDefgrdLast(A,b);
                else
                  DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) +=
                      Dexp(VOIGT3X3SYM_[a][b],j+1)*InvPlasticDefgrdLast(A,b)*2.;
              }

        // **************************************************************
        // stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
        // **************************************************************
        // tensor dSdbeta
        // in index notation this contains
        //                   d S_AB
        // dSdbeta_ABj = --------------
        //                  d beta_j
        LINALG::Matrix<6,5> dSdbeta;
        dSdbeta.Multiply(dpk2dfpinv,DFpiDbeta);

        // Calculate stiffness matrix [k^e_{d,beta}]_ij (i=1..numdof; j=1..5)
        kdbeta.MultiplyTN(detJ_w/f_bar_factor,bop,dSdbeta);
        // **************************************************************
        // end of stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
        // **************************************************************

        if (activity_state_->at(gp)==true || Dissipation>0.)
        {

        // derivative of Mandel stress w.r.t. beta
        //                      d bar Sigma_ab
        // dSigmaDbeta_abj = --------------------
        //                         d beta_j
        LINALG::Matrix<6,5> dsigmadbeta;
        dsigmadbeta.Multiply(dmdfpinv,DFpiDbeta);

        //                    d eta_ab
        // detadbeta_abj = ---------------
        //                    d beta_j
        // in Voigt notation
        // as eta is traceless, there are only 5 components
        LINALG::Matrix<5,5> detadbeta;
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
          {
            // diagonal entries
            if (i==0 || i==1)
              detadbeta(i,j) = dsigmadbeta(i,j)
              - 1./3.* ( dsigmadbeta(0,j) + dsigmadbeta(1,j) + dsigmadbeta(2,j));
            else
              detadbeta(i,j) = dsigmadbeta(i+1,j);
          }
        detadbeta.Update(2./3.*kinhard,id5,1.);

        // calculate derivative detadd
        //             d eta_ab
        // detadd = --------------
        //              d d_i
        LINALG::Matrix<5,numdofperelement_> detadd;

        // derivative of Mandel stress tensor
        LINALG::Matrix<6,numdofperelement_> dSigmadd(true);
        dSigmadd.Multiply(dmdc,bop);
        LINALG::Matrix<6,1> tmp61;
        tmp61.Multiply(dmdc,RCG);
        dSigmadd.MultiplyNT(1./3.,tmp61,htensor,1.);
        dSigmadd.Scale(f_bar_factor*f_bar_factor);

        for (int i=0; i<5; i++)
          for (int j=0; j<numdofperelement_; j++)
          {
            // diagonal entries
            if (i==0 || i==1)
              detadd(i,j) = dSigmadd(i,j) - 1./3. * (dSigmadd(0,j) + dSigmadd(1,j) + dSigmadd(2,j));
            else
              detadd(i,j) = dSigmadd(i+1,j);
          }

        LINALG::Matrix<5,1> dYpl_dbeta(true);
        LINALG::Matrix<numdofperelement_,1> dYpl_dd(true);
        if (Dissipation>0.)
        {
          for (int j=0;j<2;j++)
          {
            for (int i=0; i<numdofperelement_;i++)
            {
              dYpl_dd(i) += -Dissipation/abseta_bar/abseta_bar*(2.*eta_bar_vec(j) + eta_bar_vec((j+1)%2))/abseta_bar*detadd(j,i);
              dYpl_dd(i) -= (2.*DalphaK_last_iter_->at(gp)(j) + DalphaK_last_iter_->at(gp)((j+1)%2))*detadd(j,i)/abseta_bar;
            }
            for (int i=0; i<5;i++)
            {
              dYpl_dbeta(i) += -Dissipation/abseta_bar/abseta_bar*(2.*eta_bar_vec(j) + eta_bar_vec((j+1)%2))/abseta_bar*detadbeta(j,i);
              dYpl_dbeta(i) -= (2.*DalphaK_last_iter_->at(gp)(j)+DalphaK_last_iter_->at(gp)((j+1)%2))*detadbeta(j,i)/abseta_bar;
            }
            dYpl_dbeta(j) -= (2.*eta_bar_vec(j) + eta_bar_vec((j+1)%2))/abseta_bar;
          }
          for (int j=2; j<5; j++)
          {
            for (int i=0; i<numdofperelement_; i++)
            {
              dYpl_dd(i) += -Dissipation/abseta_bar/abseta_bar*2.*eta_bar_vec(j)*detadd(j,i)/abseta_bar;
              dYpl_dd(i) -= 2.*DalphaK_last_iter_->at(gp)(j)*detadd(j,i)/abseta_bar;
            }
            for (int i=0; i<5; i++)
            {
              dYpl_dbeta(i) += -Dissipation/abseta_bar/abseta_bar*2.*eta_bar_vec(j)*detadbeta(j,i)/abseta_bar;
              dYpl_dbeta(i) -= 2.*DalphaK_last_iter_->at(gp)(j)*detadbeta(j,i)/abseta_bar;
            }
            dYpl_dbeta(j) -= 2.*eta_bar_vec(j)/abseta_bar;
          }
          dYpl_dd.Scale(2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI))));
          dYpl_dbeta.Scale(2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI))));
        }

        if (activity_state_->at(gp)==true)
        {
          // communicate number of active plastic gauss points back to time integration
          // don't sum up for ghost elements
          if (MyPID == so3_ele::Owner())
          {
            params.get<int>("number_active_plastic_gp")++;
          }
          
          // **************************************************************
          // stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
          // linearization of the complementarity function
          // **************************************************************

          LINALG::Matrix<5,1> dabs_eta_trial_dbeta(true);
          for (int j=0;j<2; j++)
          {
            for (int i=0; i<5; i++)
              dabs_eta_trial_dbeta(i) += (2.*eta_bar_trial_vec(j)+ eta_bar_trial_vec((j+1)%2))*(detadbeta(j,i)-id5(i,j)*cpl)/absetatrial_bar;
          }
          for (int j=2; j<5; j++)
          {
            for (int i=0; i<5; i++)
              dabs_eta_trial_dbeta(i) += (2.*eta_bar_trial_vec(j))*(detadbeta(j,i)-id5(i,j)*cpl)/absetatrial_bar;
          }

          // build kbb from all previous linearizations
          KbbInv_->at(gp).Clear();
          KbbInv_->at(gp).Update(1.-Ypl/absetatrial_bar,detadbeta,1.);
          KbbInv_->at(gp).Update(Ypl/absetatrial_bar*cpl,id5,1.);
          KbbInv_->at(gp).MultiplyNT(-1./absetatrial_bar,eta_bar_trial_vec,dYpl_dbeta,1.);
          KbbInv_->at(gp).MultiplyNT((1-stab_s)*Ypl/absetatrial_bar/absetatrial_bar,eta_bar_trial_vec,dabs_eta_trial_dbeta,1.);
          KbbInv_->at(gp).MultiplyNT(stab_s*apl/absetatrial_bar,eta_bar_vec,dabs_eta_trial_dbeta,1.);

          // **************************************************************
          // end of stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
          // **************************************************************

          // **************************************************************
          // stiffness matrix [k^e_{beta d}]_ij (i=1..5; j=1..numdof)
          // linearization of the complementarity function
          // **************************************************************

          LINALG::Matrix<numdofperelement_,1> dabs_eta_trial_dd(true);
          for (int j=0;j<2; j++)
            for (int i=0; i<numdofperelement_; i++)
              dabs_eta_trial_dd(i) += (2.*eta_bar_trial_vec(j)+ eta_bar_trial_vec((j+1)%2))*(detadd(j,i))/absetatrial_bar;
          for (int j=2; j<5; j++)
            for (int i=0; i<numdofperelement_; i++)
              dabs_eta_trial_dd(i) += (2.*eta_bar_trial_vec(j))*(detadd(j,i))/absetatrial_bar;

          Kbd_->at(gp).Clear();
          Kbd_->at(gp).Update(1.-1.*Ypl/absetatrial_bar,detadd,1.);
          Kbd_->at(gp).MultiplyNT(stab_s*apl/absetatrial_bar,eta_bar_vec,dabs_eta_trial_dd,1.);
          Kbd_->at(gp).MultiplyNT((1.-stab_s)*Ypl/absetatrial_bar/absetatrial_bar,eta_bar_trial_vec,dabs_eta_trial_dd,1.);
          Kbd_->at(gp).MultiplyNT(-1./absetatrial_bar,eta_bar_trial_vec,dYpl_dd,1.);

          // **************************************************************
          // end of stiffness matrix [k^e_{beta d}]_ij (i=1..5; j=1..numdof)
          // **************************************************************


          // **************************************************************
          // right hand side term for complementarity function f^int_i (i=1..5)
          // **************************************************************
          fbeta_->at(gp).Update(eta_bar_vec);
          fbeta_->at(gp).Update(-1.*Ypl/absetatrial_bar,eta_bar_trial_vec,1.);

          // **************************************************************
          // end of right hand side term for complementarity function f^int_i (i=1..5)
          // **************************************************************
        } // acrtive Gauss points

        // inactive Gauss point with plastic history within this load/time step
        else if (activity_state_->at(gp)==false && Dissipation>0.)
        {
          // the complementarity function is
          // C^pl = - Ypl^s * cplparam_ * delta alpha^k

          // Complementarity function independent from displacements
          Kbd_->at(gp).Clear();
          Kbd_->at(gp).MultiplyNT(stab_s*cpl/Ypl,DalphaK_last_iter_->at(gp),dYpl_dd,1.);

          KbbInv_->at(gp).Update(cpl,id5,0.);
          KbbInv_->at(gp).MultiplyNT(stab_s/Ypl * cpl,DalphaK_last_iter_->at(gp),dYpl_dbeta,1.);

          // right hand side term
          fbeta_->at(gp).Update(cpl,DalphaK_last_iter_->at(gp),0.);
        }

        // **************************************************************
        // static condensation of inner variables
        // **************************************************************
        //inverse matrix block [k_beta beta]_ij
        Epetra_SerialDenseMatrix Kbetabeta_epetra(5,5);
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            Kbetabeta_epetra(i,j) = KbbInv_->at(gp)(i,j);
        // we need the inverse of K_beta beta
        Epetra_SerialDenseSolver solve_for_inverseKbb;
        solve_for_inverseKbb.SetMatrix(Kbetabeta_epetra);
        solve_for_inverseKbb.Invert();
        for (int i=0; i<5; i++)
          for (int j=0; j<5; j++)
            KbbInv_->at(gp)(i,j) = Kbetabeta_epetra(i,j);

        LINALG::Matrix<numdofperelement_,5> KdbKbb; // temporary  Kdb.Kbb^-1
        KdbKbb.Multiply(kdbeta,KbbInv_->at(gp));

        // "plastic displacement stiffness"
        // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
        if (stiffmatrix!=NULL) stiffmatrix->Multiply(-1.,KdbKbb,Kbd_->at(gp),1.);

        // "plastic internal force"
        // plFint = [K_db.K_bb^-1].f_b
        if (force!=NULL) force->Multiply(-1.,KdbKbb,fbeta_->at(gp),1.);
      }
      else if (activity_state_->at(gp)==false && Dissipation<0.)
      {
        // C^{pl} = cpl*DeltaAlphaK
        // independent of displacements
        Kbd_->at(gp).Clear();
        // We can state the inverse right away
        KbbInv_->at(gp).Update(1./cpl,id5);
        // right hand side
        fbeta_->at(gp).Update(cpl,DalphaK_last_iter_->at(gp));

        // condensation to internal force vector
        if (force!=NULL)
          force->Multiply(-1./cpl,kdbeta,fbeta_->at(gp),1.);
      }
    } // modification for plastic Gauss points
    // reset for elastic GP with no plastic flow increment in this time step
    else
    {
      Kbd_->at(gp).Clear();
      fbeta_->at(gp).Clear();
      KbbInv_->at(gp).Clear();
    }

    // square of the residual L2 norm
    double residual_norm_sqare=fbeta_->at(gp)(0)*fbeta_->at(gp)(0)
                              +fbeta_->at(gp)(1)*fbeta_->at(gp)(1)
                              +(-fbeta_->at(gp)(0)-fbeta_->at(gp)(1))*(-fbeta_->at(gp)(0)-fbeta_->at(gp)(1))
                              +fbeta_->at(gp)(2)*fbeta_->at(gp)(2)*2.
                              +fbeta_->at(gp)(3)*fbeta_->at(gp)(3)*2.
                              +fbeta_->at(gp)(4)*fbeta_->at(gp)(4)*2.;
        params.get<double>("Lp_residual_square")+=residual_norm_sqare;

  } // gp loop

  // communicate unconverged active set to time integration
  if (converged_active_set==false)
    params.set("unconverged_active_set",true);

  return;
}

/*----------------------------------------------------------------------*
 | internal force, stiffness and mass for f-bar elements    seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::nln_stiffmassHill_fbar(
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
  // check for f-bar element
  DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8> * eleFBAR
    = dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>* >(this);
  if (!eleFBAR)
    dserror("invalid call for nln_stiffmass_fbar from non-fbar element");

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
  LINALG::Matrix<3,nen_> N_XYZ;
  // build deformation gradient wrt to material configuration
  LINALG::Matrix<3,3> defgrd(false);
  // shape functions and their first derivatives
  LINALG::Matrix<nen_,1> shapefunct;
  LINALG::Matrix<3,nen_> deriv;

  // ---------------------- deformation gradient at centroid of element
  double detF_0 = -1.0;
  LINALG::Matrix<3,3> invdefgrd_0(false);
  LINALG::Matrix<3,nen_> N_XYZ_0(false);
  //element coordinate derivatives at centroid
  LINALG::Matrix<3,nen_> N_rst_0(false);
  DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, DRT::Element::hex8);

  //inverse jacobian matrix at centroid
  LINALG::Matrix<3,3> invJ_0(false);
  invJ_0.Multiply(N_rst_0,xrefe);
  invJ_0.Invert();
  //material derivatives at centroid
  N_XYZ_0.Multiply(invJ_0,N_rst_0);

  //deformation gradient and its determinant at centroid
  LINALG::Matrix<3,3> defgrd_0(false);
  defgrd_0.MultiplyTT(xcurr,N_XYZ_0);
  invdefgrd_0.Invert(defgrd_0);
  detF_0 = defgrd_0.Determinant();

  // 3x3 LINALG matrix for temporary stuff
  LINALG::Matrix<3,3> tmp1(false);
  LINALG::Matrix<3,3> tmp2(false);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // get plastic material parameters
  double kinhard=plmat->Kinhard();
  double isohard=plmat->Isohard();
  double expisohard=plmat->Expisohard();
  double infyield=plmat->Infyield();
  double inityield=plmat->Inityield();
  double cpl=plparams_->get<double>("SEMI_SMOOTH_CPL");
  double stab_s=plparams_->get<double>("STABILIZATION_S");
  LINALG::Matrix<5,5> PlAniso(plmat->PlAniso());
  LINALG::Matrix<5,5> InvPlAniso(plmat->InvPlAniso());
  double PlSpinEta=-plmat->PlSpinEta();

  // converged active set
  bool converged_active_set=true;

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
    LINALG::Matrix<3,3> invdefgrd;
    double detF=invdefgrd.Invert(defgrd);

    // check for element distortion
    if (detF<=0. || detF_0<=0.)
      dserror("element distortion too large");


    // F_bar deformation gradient =(detF_0/detF)^1/3*F
    LINALG::Matrix<3,3> defgrd_bar(defgrd);
    double f_bar_factor=pow(detF_0/detF,1/3.);
    defgrd_bar.Scale(f_bar_factor);

    // calcualte total rcg/gl/ea for output
    LINALG::Matrix<3,3> total_cauchygreen(false);
    total_cauchygreen.MultiplyTN(defgrd,defgrd);
    // total Cauchy green in voigt notation
    LINALG::Matrix<6,1> RCG;
    for (int i=0; i<3; i++) RCG(i)=total_cauchygreen(i,i);
    RCG(3)=total_cauchygreen(0,1)*2.;
    RCG(4)=total_cauchygreen(1,2)*2.;
    RCG(5)=total_cauchygreen(0,2)*2.;

    // calculate nonlinear B-operator
    LINALG::Matrix<numstr_,numdofperelement_> bop(false);
    CalculateBop(&bop,&defgrd,&N_XYZ);

    // strain output *********************************
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
    // end of strain output **************************

    // recover condensed variables from last iteration step *********
    // first part
    LINALG::Matrix<8,1> tmp81;
    tmp81.Multiply(KbbInvHill_->at(gp),fbetaHill_->at(gp));

    // second part
    LINALG::Matrix<8,numdofperelement_> tmp824;
    tmp824.Multiply(KbbInvHill_->at(gp),KbdHill_->at(gp));
    tmp81.Multiply(1.,tmp824,res_d,1.);
    mDLp_last_iter_->at(gp).Update(-1.,tmp81,1.);
    double increment_norm_sqare=tmp81(0)*tmp81(0)
                              +tmp81(1)*tmp81(1)
                              +(-tmp81(0)-tmp81(1))*(-tmp81(0)-tmp81(1));
    for (int i=2; i<8; i++) increment_norm_sqare+=tmp81(i)*tmp81(i)*2.;
    params.get<double>("Lp_increment_square")+=increment_norm_sqare;
    // end of recover *********************************************

     // current plastic flow increment =-Delta Lp
     LINALG::Matrix<nsd_,nsd_> mDLp(false);
     // current kinematic hardening increment
     LINALG::Matrix<nsd_,nsd_> DeltaAlphaK(false);
     LINALG::Matrix<5,1> DeltaAlphaK_vec(false);
     for (int i=0; i<5; i++) DeltaAlphaK_vec(i) = mDLp_last_iter_->at(gp)(i);
     // -Dp
     DeltaAlphaK(0,0) = mDLp_last_iter_->at(gp)(0);
     DeltaAlphaK(1,1) = mDLp_last_iter_->at(gp)(1);
     DeltaAlphaK(2,2) = -1.0*(mDLp_last_iter_->at(gp)(0)+mDLp_last_iter_->at(gp)(1));
     DeltaAlphaK(0,1) = mDLp_last_iter_->at(gp)(2);
     DeltaAlphaK(1,0) = mDLp_last_iter_->at(gp)(2);
     DeltaAlphaK(1,2) = mDLp_last_iter_->at(gp)(3);
     DeltaAlphaK(2,1) = mDLp_last_iter_->at(gp)(3);
     DeltaAlphaK(0,2) = mDLp_last_iter_->at(gp)(4);
     DeltaAlphaK(2,0) = mDLp_last_iter_->at(gp)(4);
     mDLp.Update(DeltaAlphaK);
     // Wp
     LINALG::Matrix<3,3> Wp(true);
     Wp(0,1) -= mDLp_last_iter_->at(gp)(5);
     Wp(1,0) += mDLp_last_iter_->at(gp)(5);
     Wp(1,2) -= mDLp_last_iter_->at(gp)(6);
     Wp(2,1) += mDLp_last_iter_->at(gp)(6);
     Wp(0,2) -= mDLp_last_iter_->at(gp)(7);
     Wp(2,0) += mDLp_last_iter_->at(gp)(7);
     mDLp.Update(-1.,Wp,1.);

     // inverse plastic deformation gradient
      LINALG::Matrix<nsd_,nsd_> InvPlasticDefgrd(false);
      // inverse plastic deformation gradient at last time step
      LINALG::Matrix<nsd_,nsd_> InvPlasticDefgrdLast = last_plastic_defgrd_inverse_->at(gp);

     // compute matrix exponential
     tmp1=mDLp;
     MatrixExponential3x3(tmp1);
     InvPlasticDefgrd.Multiply(InvPlasticDefgrdLast,tmp1);

     // material call
     LINALG::Matrix<numstr_,1> pk2_stress_bar(true);
     LINALG::Matrix<6,6> Cmat_ABCD_bar(true);
     LINALG::Matrix<6,9> dpk2dfpinv(true);
     LINALG::Matrix<3,3> mandelstress(true);
     LINALG::Matrix<6,6> dmdc(true);
     LINALG::Matrix<6,9> dmdfpinv(true);
     params.set<int>("gp",gp);
     params.set<int>("eleID",Id());
     plmat->Evaluate(&defgrd_bar,&InvPlasticDefgrd,params,&pk2_stress_bar,&Cmat_ABCD_bar,&dpk2dfpinv,&mandelstress,&dmdc,&dmdfpinv);

     // return gp stresses
     switch (iostress)
     {
     case INPAR::STR::stress_2pk:
     {
       if (elestress == NULL) dserror("stress data not available");
       for (int i = 0; i < numstr_; ++i)
           (*elestress)(gp,i) = pk2_stress_bar(i);
     }
     break;
     case INPAR::STR::stress_cauchy:
     {
       if (elestress == NULL) dserror("stress data not available");
       const double detF = defgrd.Determinant();

       LINALG::Matrix<3,3> pkstress;
       pkstress(0,0) = pk2_stress_bar(0);
       pkstress(0,1) = pk2_stress_bar(3);
       pkstress(0,2) = pk2_stress_bar(5);
       pkstress(1,0) = pkstress(0,1);
       pkstress(1,1) = pk2_stress_bar(1);
       pkstress(1,2) = pk2_stress_bar(4);
       pkstress(2,0) = pkstress(0,2);
       pkstress(2,1) = pkstress(1,2);
       pkstress(2,2) = pk2_stress_bar(2);

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

     // equivalent stress eta
     LINALG::Matrix<nsd_,nsd_> eta(mandelstress);
     for (int i=0; i<nsd_; i++)
       eta(i,i) -= 1./3.*(mandelstress(0,0) + mandelstress(1,1) + mandelstress(2,2));
     eta.Update(2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
     eta.Update(2./3.*kinhard,DeltaAlphaK,1.);
     LINALG::Matrix<5,1> eta_vec(false);
     eta_vec(0) = eta(0,0);
     eta_vec(1) = eta(1,1);
     eta_vec(2) = 0.5*(eta(0,1)+eta(1,0));
     eta_vec(3) = 0.5*(eta(2,1)+eta(1,2));
     eta_vec(4) = 0.5*(eta(0,2)+eta(2,0));
     LINALG::Matrix<5,1> Aeta(false);
     Aeta.Multiply(PlAniso,eta_vec);

     // eta_trial
     LINALG::Matrix<5,1> eta_trial_vec(eta_vec);
     eta_trial_vec.Multiply(-cpl,InvPlAniso,DeltaAlphaK_vec,1.);
     LINALG::Matrix<5,1> AetaTr(false);
     AetaTr.Multiply(PlAniso,eta_trial_vec);

     // Matrix norms
     double AnormEtatrial=0.; // sqrt( eta_tr : A : eta_tr )
     LINALG::Matrix<1,1> tmp11;
     LINALG::Matrix<5,1> tmp51;
     tmp1(0,0) = 2./3.*AetaTr(0)-1./3.*AetaTr(1);
     tmp1(1,1) =-1./3.*AetaTr(0)+2./3.*AetaTr(1);
     tmp1(2,2) = -tmp1(0,0)-tmp1(1,1);
     tmp1(0,1) = AetaTr(2)/2.;
     tmp1(1,0) = AetaTr(2)/2.;
     tmp1(1,2) = AetaTr(3)/2.;
     tmp1(2,1) = AetaTr(3)/2.;
     tmp1(0,2) = AetaTr(4)/2.;
     tmp1(2,0) = AetaTr(4)/2.;
     tmp11.MultiplyTN(eta_trial_vec,AetaTr);
     AnormEtatrial=sqrt(tmp11(0,0));

     tmp1(0,0) = 2./3.*Aeta(0)-1./3.*Aeta(1);
     tmp1(1,1) =-1./3.*Aeta(0)+2./3.*Aeta(1);
     tmp1(2,2) = -tmp1(0,0)-tmp1(1,1);
     tmp1(0,1) = Aeta(2)/2.;
     tmp1(1,0) = Aeta(2)/2.;
     tmp1(1,2) = Aeta(3)/2.;
     tmp1(2,1) = Aeta(3)/2.;
     tmp1(0,2) = Aeta(4)/2.;
     tmp1(2,0) = Aeta(4)/2.;
     double EucNormAeta = tmp1.Norm2();
     tmp11.MultiplyTN(eta_vec,Aeta);
     double AnormEta=sqrt(tmp11(0,0));

     double DpAe=0.;
     for (int i=0; i<3; i++)
       for (int j=0; j<3; j++)
         DpAe+=DeltaAlphaK(i,j)*tmp1(i,j);

     if (DpAe<0.)
       deltaAlphaI_->at(gp) = -sqrt(2./3.)*DpAe*AnormEta/EucNormAeta/EucNormAeta;

     double Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_->at(gp))))
         + isohard*(last_alpha_isotropic_->at(gp)(0,0)+  deltaAlphaI_->at(gp)) +inityield);

     // check activity state
     // inactive
     if (Ypl<AnormEtatrial)
     {
       if (activity_state_->at(gp)==false) // gp switches state
         if (abs(Ypl-AnormEtatrial)>AS_CONVERGENCE_TOL*inityield
             || mDLp_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
           converged_active_set = false;
       activity_state_->at(gp) = true;
     }
     // active
     else
     {
       if (activity_state_->at(gp)==true) // gp switches state
         if (abs(Ypl-AnormEtatrial)>AS_CONVERGENCE_TOL*inityield
             || mDLp_last_iter_->at(gp).NormInf()>AS_CONVERGENCE_TOL*inityield/cpl)
           converged_active_set = false;
       activity_state_->at(gp) = false;
     }

      // integrate usual internal force and stiffness matrix
      double detJ_w = detJ*intpoints_.Weight(gp);
      // integrate elastic internal force vector **************************
      // update internal force vector
      if (force != NULL)
      {
        force->MultiplyTN(detJ_w/f_bar_factor, bop, pk2_stress_bar, 1.0);
      }

      // additional f-bar derivatives
      LINALG::Matrix<numdofperelement_,1> htensor(true);

      // update stiffness matrix
      if (stiffmatrix != NULL)
      {
        for(int n=0;n<numdofperelement_;n++)
          for(int i=0;i<3;i++)
            htensor(n) += invdefgrd_0(i,n%3)*N_XYZ_0(i,n/3)-invdefgrd(i,n%3)*N_XYZ(i,n/3);

        // integrate `elastic' and `initial-displacement' stiffness matrix
        // keu = keu + (B^T . C . B) * detJ * w(gp)
        LINALG::Matrix<6,numdofperelement_> cb;
        cb.Multiply(Cmat_ABCD_bar,bop);
        stiffmatrix->MultiplyTN(detJ_w*f_bar_factor,bop,cb,1.0);

        // integrate `geometric' stiffness matrix and add to keu *****************
        LINALG::Matrix<6,1> sfac(pk2_stress_bar); // auxiliary integrated stress
        sfac.Scale(detJ_w/f_bar_factor); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
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
         LINALG::Matrix<numstr_,1> ccg;
         ccg.Multiply(Cmat_ABCD_bar,RCG);

         LINALG::Matrix<numdofperelement_,1> bopccg(false); // auxiliary integrated stress
         bopccg.MultiplyTN(detJ_w*f_bar_factor/3.0,bop,ccg);

         LINALG::Matrix<numdofperelement_,1> bops(false); // auxiliary integrated stress
         bops.MultiplyTN(-detJ_w/f_bar_factor/3.0,bop,pk2_stress_bar);
         stiffmatrix->MultiplyNT(1.,bops,htensor,1.);
         stiffmatrix->MultiplyNT(1.,bopccg,htensor,1.);
         // end of integrate additional fbar matrix*****************************
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
      if (
          (stiffmatrix!=NULL || force!=NULL)
          &&
          (activity_state_->at(gp)==true || mDLp_last_iter_->at(gp).NormInf()!=0.)
         )
      {
        // variables needed for condensation and calculated seperately for active and inactive Gauss points
        LINALG::Matrix<8,numdofperelement_> kbetad(true);
        LINALG::Matrix<8,8> kbetabeta(true);
        LINALG::Matrix<8,1> force_beta(true);
        LINALG::Matrix<numdofperelement_,8> kdbeta;

        // derivative of symmetric complementarity function w.r.t. displacements
        LINALG::Matrix<5,numdofperelement_> DcplSymDd(true);
        // derivative of symmetric complementarity function
        LINALG::Matrix<5,8> DcplSymDbeta(true);
        // complementarity function
        LINALG::Matrix<5,1> CplSym(true);
        // plastic spin equation
        LINALG::Matrix<3,1> SpEq(true);
        // derivative of spin equation w.r.t. displacements
        LINALG::Matrix<3,numdofperelement_> dSpEqDd(true);
        // derivative of spin equation w.r.t. beta
        LINALG::Matrix<3,8> dSpEqDbeta(true);

        // delta alpha_K d beta
        LINALG::Matrix<5,8> DalphaKdbeta(true);
        for (int i=0; i<5; i++)
          DalphaKdbeta(i,i)=1.;

          // damping parameter apl
          double apl=1.;
          if (Ypl/AnormEta<1.) apl=Ypl/AnormEta;

          // derivative of the matrix exponential
          LINALG::Matrix<9,9> Dexp(false);
          MatrixExponentialDerivative3x3(mDLp,Dexp);

          // derivative of inverse plastic deformation gradient w.r.t. flow increment
          LINALG::Matrix<9,8> DFpiDbeta(true);
          for (int a=0; a<3; a++)
            for (int A=0; A<3; A++)
              for (int j=0; j<8; j++)
                for (int b=0; b<3; b++)
                {
                  // main diagonal entries
                  if (j==0 || j==1)
                    DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) +=
                    (Dexp(VOIGT3X3NONSYM_[b][a],j)-Dexp(VOIGT3X3NONSYM_[b][a],2))*InvPlasticDefgrdLast(A,b);
                  if (j==2 || j==3 || j==4)
                    DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) += (Dexp(VOIGT3X3NONSYM_[b][a],j+1)+Dexp(VOIGT3X3NONSYM_[b][a],j+4))*InvPlasticDefgrdLast(A,b);
                  if (j==5 ||j==6 || j==7)
                    DFpiDbeta(VOIGT3X3NONSYM_[A][a],j) += (Dexp(VOIGT3X3NONSYM_[b][a],j-2)-Dexp(VOIGT3X3NONSYM_[b][a],j+1))*InvPlasticDefgrdLast(A,b);
                }

          // **************************************************************
          // stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
          // **************************************************************
          // tensor dSdbeta
          // in index notation this contains
          //                   d S_AB
          // dSdbeta_ABj = --------------
          //                  d beta_j
          LINALG::Matrix<6,8> dSdbeta;
          dSdbeta.Multiply(dpk2dfpinv,DFpiDbeta);

          // Calculate stiffness matrix [k^e_{d,beta}]_ij (i=1..numdof; j=1..5)
          kdbeta.MultiplyTN(detJ_w/f_bar_factor,bop,dSdbeta);
          // **************************************************************
          // end of stiffness matrix [k^e_{d beta}]_ij (i=1..numdof; j=1..5)
          // **************************************************************

          // Due to stabilization, we have to treat inactive nodes as well
          if (activity_state_->at(gp)==true || DpAe>0.)
          {
          // calculate derivative detadd
           //             d eta_ab
           // detadd = --------------
           //              d d_i
           LINALG::Matrix<5,numdofperelement_> detadd;

           // derivative of Mandel stress tensor
           LINALG::Matrix<6,numdofperelement_> dSigmadd(true);
           dSigmadd.Multiply(dmdc,bop);
           LINALG::Matrix<6,1> tmp61;
           tmp61.Multiply(dmdc,RCG);
           dSigmadd.MultiplyNT(1./3.,tmp61,htensor,1.);
           dSigmadd.Scale(f_bar_factor*f_bar_factor);

           for (int i=0; i<5; i++)
             for (int j=0; j<numdofperelement_; j++)
             {
               // diagonal entries
               if (i==0 || i==1)
                 detadd(i,j) = dSigmadd(i,j) - 1./3. * (dSigmadd(0,j) + dSigmadd(1,j) + dSigmadd(2,j));
               else
                 detadd(i,j) = dSigmadd(i+1,j);
             }

            // derivative of Mandel stress w.r.t. beta
            //                      d bar Sigma_ab
            // dSigmaDbeta_abj = --------------------
            //                         d beta_j
            LINALG::Matrix<6,8> dsigmadbeta;
            dsigmadbeta.Multiply(dmdfpinv,DFpiDbeta);

            //                    d eta_ab
            // detadbeta_abj = ---------------
            //                    d beta_j
            // in Voigt notation
            // as eta is traceless, there are only 5 components
            LINALG::Matrix<5,8> detadbeta;
            LINALG::Matrix<5,8> detatrialdbeta;
            for (int i=0; i<5; i++)
              for (int j=0; j<8; j++)
              {
                // diagonal entries
                if (i==0 || i==1)
                {
                  detadbeta(i,j) = dsigmadbeta(i,j)
                                           - 1./3.* ( dsigmadbeta(0,j) + dsigmadbeta(1,j) + dsigmadbeta(2,j)) + (2./3.*kinhard*(i==j));
                }
                else
                {
                  detadbeta(i,j) = dsigmadbeta(i+1,j) + (2./3.*kinhard*(i==j));
                }
                detatrialdbeta(i,j) = detadbeta(i,j);
                if (j<5)
                  detatrialdbeta(i,j) -= cpl*InvPlAniso(i,j);
              }

            LINALG::Matrix<8,1> dAnormEtaDbeta(true);
             dAnormEtaDbeta.MultiplyTN(1./AnormEta,detadbeta,Aeta);


             LINALG::Matrix<8,1> dYpl_dbeta(true);
             LINALG::Matrix<numdofperelement_,1> dYpl_dd(true);
             if (DpAe<0.)
             {
               LINALG::Matrix<5,8> dAetaDbeta;
               dAetaDbeta.Multiply(PlAniso,detadbeta);
               LINALG::Matrix<5,numdofperelement_> dAetaDd;
               dAetaDd.Multiply(PlAniso,detadd);
               LINALG::Matrix<8,1> dAnormEtaDbeta(true);
                dAnormEtaDbeta.MultiplyTN(1./AnormEta,detadbeta,Aeta);
                LINALG::Matrix<numdofperelement_,1> dAnormEtaDd(true);
                 dAnormEtaDd.MultiplyTN(1./AnormEta,detadd,Aeta);
               for (int i=0; i<8; i++)
               {
                 dYpl_dbeta(i) += 4./3.*Aeta(0)*dAetaDbeta(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dbeta(i) -= 2./3.*Aeta(1)*dAetaDbeta(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dbeta(i) -= 2./3.*Aeta(0)*dAetaDbeta(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dbeta(i) += 4./3.*Aeta(1)*dAetaDbeta(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dbeta(i) += Aeta(2)*dAetaDbeta(2,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dbeta(i) += Aeta(3)*dAetaDbeta(3,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dbeta(i) += Aeta(4)*dAetaDbeta(4,i) *AnormEta*DpAe/pow(EucNormAeta,4.);

                 dYpl_dbeta(i) -= DeltaAlphaK_vec(0)*dAetaDbeta(0,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dbeta(i) -= DeltaAlphaK_vec(2)*dAetaDbeta(2,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dbeta(i) -= DeltaAlphaK_vec(3)*dAetaDbeta(3,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dbeta(i) -= DeltaAlphaK_vec(4)*dAetaDbeta(4,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dbeta(i) -= DeltaAlphaK_vec(1)*dAetaDbeta(1,i)*AnormEta/EucNormAeta/EucNormAeta;
               }

               dYpl_dbeta(0) -= Aeta(0)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(1) -= Aeta(1)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(2) -= Aeta(2)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(3) -= Aeta(3)*AnormEta/EucNormAeta/EucNormAeta;
               dYpl_dbeta(4) -= Aeta(4)*AnormEta/EucNormAeta/EucNormAeta;

               dYpl_dbeta.Update(-DpAe/EucNormAeta/EucNormAeta,dAnormEtaDbeta,1.);

               for (int i=0; i<numdofperelement_; i++)
               {
                 dYpl_dd(i) += 4./3.*Aeta(0)*dAetaDd(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dd(i) -= 2./3.*Aeta(1)*dAetaDd(0,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dd(i) -= 2./3.*Aeta(0)*dAetaDd(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dd(i) += 4./3.*Aeta(1)*dAetaDd(1,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dd(i) += Aeta(2)*dAetaDd(2,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dd(i) += Aeta(3)*dAetaDd(3,i) *AnormEta*DpAe/pow(EucNormAeta,4.);
                 dYpl_dd(i) += Aeta(4)*dAetaDd(4,i) *AnormEta*DpAe/pow(EucNormAeta,4.);

                 dYpl_dd(i) -= DeltaAlphaK_vec(0)*dAetaDd(0,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dd(i) -= DeltaAlphaK_vec(2)*dAetaDd(2,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dd(i) -= DeltaAlphaK_vec(3)*dAetaDd(3,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dd(i) -= DeltaAlphaK_vec(4)*dAetaDd(4,i)*AnormEta/EucNormAeta/EucNormAeta;
                 dYpl_dd(i) -= DeltaAlphaK_vec(1)*dAetaDd(1,i)*AnormEta/EucNormAeta/EucNormAeta;
               }

               dYpl_dd.Update(-DpAe/EucNormAeta/EucNormAeta,dAnormEtaDd,1.);

               dYpl_dbeta.Scale(2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI_->at(gp)))));
               dYpl_dd.Scale(   2./3.*(isohard+(infyield-inityield)*expisohard*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+deltaAlphaI_->at(gp)))));

             }

          if (activity_state_->at(gp)==true)
          {
            // communicate number of active plastic gauss points back to time integration
            // don't sum up for ghost elements
            if (MyPID == so3_ele::Owner()) params.get<int>("number_active_plastic_gp")++;

            // Symmetric complementarity function
            CplSym.Update(1.,eta_vec,0.);
            CplSym.Update(-Ypl/AnormEtatrial,eta_trial_vec,1.);

            // **************************************************************
            // stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
            // linearization of the complementarity function
            // **************************************************************

            LINALG::Matrix<5,1> Aetatr;
            Aetatr.Multiply(PlAniso,eta_trial_vec);
            LINALG::Matrix<8,1> dAnormEtaTrialDbeta(true);
            dAnormEtaTrialDbeta.MultiplyTN(1./AnormEtatrial,detatrialdbeta,Aetatr);

            DcplSymDbeta.Update(1.-Ypl/AnormEtatrial,detadbeta,1.);
            DcplSymDbeta.Multiply(cpl*Ypl/AnormEtatrial,InvPlAniso,DalphaKdbeta,1.);
            DcplSymDbeta.MultiplyNT(-1./AnormEtatrial,eta_trial_vec,dYpl_dbeta,1.);
            DcplSymDbeta.MultiplyNT((1.-stab_s)*Ypl/AnormEtatrial/AnormEtatrial,eta_trial_vec,dAnormEtaTrialDbeta,1.);
            DcplSymDbeta.MultiplyNT(stab_s*apl/AnormEtatrial,eta_vec,dAnormEtaTrialDbeta,1.);

            // **************************************************************
            // end of stiffness matrix [k^e_{beta beta}]_ij (i=1..5; j=1..5)
            // **************************************************************

            // **************************************************************
            // stiffness matrix [k^e_{beta d}]_ij (i=1..5; j=1..numdof)
            // linearization of the complementarity function
            // **************************************************************
            LINALG::Matrix<numdofperelement_,1> dAnormEtatrialDd(true);
            dAnormEtatrialDd.MultiplyTN(1./AnormEtatrial,detadd,Aetatr,1.);
            LINALG::Matrix<numdofperelement_,1> dAnormEtaDd(true);
            dAnormEtaDd.MultiplyTN(1./AnormEta,detadd,Aeta,1.);

            DcplSymDd.Update(1.-Ypl/AnormEtatrial,detadd,1.);
            DcplSymDd.MultiplyNT((1.-stab_s)*Ypl/AnormEtatrial/AnormEtatrial,eta_trial_vec,dAnormEtatrialDd,1.);
            DcplSymDd.MultiplyNT(stab_s*apl/AnormEtatrial,eta_vec,dAnormEtatrialDd,1.);
            DcplSymDd.MultiplyNT(1./AnormEtatrial,eta_trial_vec,dYpl_dd,1.);
            // **************************************************************
            // end of stiffness matrix [k^e_{beta d}]_ij (i=1..5; j=1..numdof)
            // **************************************************************

          } // acrtive Gauss points

          // inactive Gauss point with plastic history within this load/time step
          else
          {
            LINALG::Matrix<5,1> AiDalphak;
            AiDalphak.Multiply(InvPlAniso,DeltaAlphaK_vec);

            DcplSymDd.Clear();
            DcplSymDd.MultiplyNT(stab_s*cpl/Ypl,AiDalphak,dYpl_dd,1.);

            DcplSymDbeta.Multiply(cpl,InvPlAniso,DalphaKdbeta,1.);
            DcplSymDbeta.MultiplyNT(cpl*stab_s/Ypl,AiDalphak,dYpl_dbeta,1.);

            CplSym.Multiply(cpl,InvPlAniso,DeltaAlphaK_vec);
          }
          // symmetric NCP done

          // plastic spin stuff
          // Wp-eta_sp/inityield* (-eta*delta alpha_k + delta aplpha_k * eta) = 0
          // <==> Wp = eta_sp/inityield * (eta*Dp - Dp*eta)
          // This equation is apriori skew-symmetric so it reduces to 3 independent
          // scalar-valued equations.
          SpEq(0) = -mDLp_last_iter_->at(gp)(5) -PlSpinEta/inityield *
              ( DeltaAlphaK_vec(0)*mandelstress(0,1)
                  +DeltaAlphaK_vec(2)*mandelstress(1,1)
                  +DeltaAlphaK_vec(4)*mandelstress(1,2)
                  -DeltaAlphaK_vec(2)*mandelstress(0,0)
                  -DeltaAlphaK_vec(1)*mandelstress(0,1)
                  -DeltaAlphaK_vec(3)*mandelstress(0,2) );
          SpEq(1) = -mDLp_last_iter_->at(gp)(6) -PlSpinEta/inityield *
              ( DeltaAlphaK_vec(2)*mandelstress(0,2)
                  +DeltaAlphaK_vec(1)*mandelstress(1,2)*2.
                  +DeltaAlphaK_vec(3)*mandelstress(2,2)
                  -DeltaAlphaK_vec(4)*mandelstress(0,1)
                  -DeltaAlphaK_vec(3)*mandelstress(1,1)
                  +DeltaAlphaK_vec(0)*mandelstress(1,2) );
          SpEq(2) = -mDLp_last_iter_->at(gp)(7) -PlSpinEta/inityield *
              ( DeltaAlphaK_vec(0)*mandelstress(0,2)*2.
                  +DeltaAlphaK_vec(2)*mandelstress(1,2)
                  +DeltaAlphaK_vec(4)*mandelstress(2,2)
                  -DeltaAlphaK_vec(4)*mandelstress(0,0)
                  -DeltaAlphaK_vec(3)*mandelstress(0,1)
                  +DeltaAlphaK_vec(1)*mandelstress(0,2) );

          // derivative of SpEq w.r.t. beta
          for (int i=0; i<8; i++)
          {
            dSpEqDbeta(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(0)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[1][1],i);
            dSpEqDbeta(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[1][2],i);
            dSpEqDbeta(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[0][0],i);
            dSpEqDbeta(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(1)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[0][2],i);

            dSpEqDbeta(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[0][2],i);
            dSpEqDbeta(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dsigmadbeta(VOIGT3X3SYM_[1][2],i)*2.;
            dSpEqDbeta(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[2][2],i);
            dSpEqDbeta(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[1][1],i);
            dSpEqDbeta(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dsigmadbeta(VOIGT3X3SYM_[1][2],i);

            dSpEqDbeta(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dsigmadbeta(VOIGT3X3SYM_[0][2],i)*2.;
            dSpEqDbeta(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dsigmadbeta(VOIGT3X3SYM_[1][2],i);
            dSpEqDbeta(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[2][2],i);
            dSpEqDbeta(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dsigmadbeta(VOIGT3X3SYM_[0][0],i);
            dSpEqDbeta(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dsigmadbeta(VOIGT3X3SYM_[0][1],i);
            dSpEqDbeta(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dsigmadbeta(VOIGT3X3SYM_[0][2],i);
          }

          dSpEqDbeta(0,0) -=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(0,2) -=PlSpinEta/inityield*mandelstress(1,1);
          dSpEqDbeta(0,4) -=PlSpinEta/inityield*mandelstress(1,2);
          dSpEqDbeta(0,2) +=PlSpinEta/inityield*mandelstress(0,0);
          dSpEqDbeta(0,1) +=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(0,3) +=PlSpinEta/inityield*mandelstress(0,2);

          dSpEqDbeta(1,2) -=PlSpinEta/inityield*mandelstress(0,2);
          dSpEqDbeta(1,1) -=PlSpinEta/inityield*mandelstress(1,2)*2.;
          dSpEqDbeta(1,3) -=PlSpinEta/inityield*mandelstress(2,2);
          dSpEqDbeta(1,4) +=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(1,3) +=PlSpinEta/inityield*mandelstress(1,1);
          dSpEqDbeta(1,0) +=PlSpinEta/inityield*mandelstress(1,2);

          dSpEqDbeta(2,0) -=PlSpinEta/inityield*mandelstress(0,2)*2.;
          dSpEqDbeta(2,2) -=PlSpinEta/inityield*mandelstress(1,2);
          dSpEqDbeta(2,4) -=PlSpinEta/inityield*mandelstress(2,2);
          dSpEqDbeta(2,4) +=PlSpinEta/inityield*mandelstress(0,0);
          dSpEqDbeta(2,3) +=PlSpinEta/inityield*mandelstress(0,1);
          dSpEqDbeta(2,1) -=PlSpinEta/inityield*mandelstress(0,2);

          dSpEqDbeta(0,5) -=1.;
          dSpEqDbeta(1,6) -=1.;
          dSpEqDbeta(2,7) -=1.;

          // derivative of SpEq w.r.t. displacements
          for (int i=0; i<numdofperelement_; i++)
          {
            dSpEqDd(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(0)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[1][1],i);
            dSpEqDd(0,i) -= PlSpinEta/inityield*DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[1][2],i);
            dSpEqDd(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[0][0],i);
            dSpEqDd(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(1)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(0,i) += PlSpinEta/inityield*DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[0][2],i);

            dSpEqDd(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[0][2],i);
            dSpEqDd(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dSigmadd(VOIGT3X3SYM_[1][2],i)*2.;
            dSpEqDd(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[2][2],i);
            dSpEqDd(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(1,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[1][1],i);
            dSpEqDd(1,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dSigmadd(VOIGT3X3SYM_[1][2],i);

            dSpEqDd(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(0)*dSigmadd(VOIGT3X3SYM_[0][2],i)*2.;
            dSpEqDd(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(2)*dSigmadd(VOIGT3X3SYM_[1][2],i);
            dSpEqDd(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[2][2],i);
            dSpEqDd(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(4)*dSigmadd(VOIGT3X3SYM_[0][0],i);
            dSpEqDd(2,i) += PlSpinEta/inityield *DeltaAlphaK_vec(3)*dSigmadd(VOIGT3X3SYM_[0][1],i);
            dSpEqDd(2,i) -= PlSpinEta/inityield *DeltaAlphaK_vec(1)*dSigmadd(VOIGT3X3SYM_[0][2],i);
          }

//          // FD check beta
//          double epsilon=1.e-10;
//          for (int i=0; i<8; i++)
//          {
//            mDLp_last_iter_->at(gp)(i) += epsilon;
//            LINALG::Matrix<3,3> mDLp_new(false);
//            LINALG::Matrix<3,3> DeltaAlphaK_new(false);
//            LINALG::Matrix<5,1> DeltaAlphaK_vec_new(false);
//            for (int iii=0; iii<5; iii++) DeltaAlphaK_vec_new(iii) = mDLp_last_iter_->at(gp)(iii);
//            // -Dp
//            DeltaAlphaK_new(0,0) = mDLp_last_iter_->at(gp)(0);
//            DeltaAlphaK_new(1,1) = mDLp_last_iter_->at(gp)(1);
//            DeltaAlphaK_new(2,2) = -1.0*(mDLp_last_iter_->at(gp)(0)+mDLp_last_iter_->at(gp)(1));
//            DeltaAlphaK_new(0,1) = mDLp_last_iter_->at(gp)(2);
//            DeltaAlphaK_new(1,0) = mDLp_last_iter_->at(gp)(2);
//            DeltaAlphaK_new(1,2) = mDLp_last_iter_->at(gp)(3);
//            DeltaAlphaK_new(2,1) = mDLp_last_iter_->at(gp)(3);
//            DeltaAlphaK_new(0,2) = mDLp_last_iter_->at(gp)(4);
//            DeltaAlphaK_new(2,0) = mDLp_last_iter_->at(gp)(4);
//            mDLp_new.Update(DeltaAlphaK_new);
//            // -Wp
//            mDLp_new(0,1) += mDLp_last_iter_->at(gp)(5);
//            mDLp_new(1,0) -= mDLp_last_iter_->at(gp)(5);
//            mDLp_new(1,2) += mDLp_last_iter_->at(gp)(6);
//            mDLp_new(2,1) -= mDLp_last_iter_->at(gp)(6);
//            mDLp_new(0,2) += mDLp_last_iter_->at(gp)(7);
//            mDLp_new(2,0) -= mDLp_last_iter_->at(gp)(7);
//            double absDeltaAlphaK_new=DeltaAlphaK_new.Norm2();
//            tmp1=mDLp_new;
//            MatrixExponential3x3(tmp1);
//
//            LINALG::Matrix<3,3> InvPlasticDefgrd_new;
//            InvPlasticDefgrd_new.Multiply(InvPlasticDefgrdLast,tmp1);
//
////            for (int i=0; i<3; i++)
////              for (int j=0; j<3; j++)
////              {
////                LINALG::Matrix<3,3> InvPlasticDefgrd_new(InvPlasticDefgrd);
////                InvPlasticDefgrd_new(i,j) +=epsilon;
//
//            // material call
//            LINALG::Matrix<numstr_,1> pk2_stress_new(true);
//            LINALG::Matrix<6,6> Cmat_ABCD_new(true);
//            LINALG::Matrix<6,9> dpk2dfpinv_new(true);
//            LINALG::Matrix<3,3> mandelstress_new(true);
//            LINALG::Matrix<6,6> dmdc_new(true);
//            LINALG::Matrix<6,9> dmdfpinv_new(true);
//            params.set<int>("gp",gp);
//            params.set<int>("eleID",Id());
//            plmat->Evaluate(&defgrd_bar,&InvPlasticDefgrd_new,params,&pk2_stress_new,&Cmat_ABCD_new,&dpk2dfpinv_new,&mandelstress_new,&dmdc_new,&dmdfpinv_new);
//
//            // equivalent stress eta
//            LINALG::Matrix<nsd_,nsd_> eta_new(mandelstress_new);
//            for (int ii=0; ii<nsd_; ii++)
//              eta_new(ii,ii) -= 1./3.*(mandelstress_new(0,0) + mandelstress_new(1,1) + mandelstress_new(2,2));
//            eta_new.Update(2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
//            eta_new.Update(2./3.*kinhard,DeltaAlphaK_new,1.);
//            LINALG::Matrix<5,1> eta_vec_new(false);
//            eta_vec_new(0) = eta_new(0,0);
//            eta_vec_new(1) = eta_new(1,1);
//            eta_vec_new(2) = 0.5*(eta_new(0,1)+eta_new(1,0));
//            eta_vec_new(3) = 0.5*(eta_new(2,1)+eta_new(1,2));
//            eta_vec_new(4) = 0.5*(eta_new(0,2)+eta_new(2,0));
//
//            LINALG::Matrix<5,1> Aeta_new;
//            Aeta_new.Multiply(PlAniso,eta_vec_new);
//
//            // eta_trial
//            LINALG::Matrix<5,1> eta_trial_vec_new(eta_vec_new);
//            eta_trial_vec_new.Multiply(-cpl,InvPlAniso,DeltaAlphaK_vec_new,1.);
//
//            double AnormEtatrial_new=0.; // sqrt( eta_tr : A : eta_tr )
//            LINALG::Matrix<1,1> tmp11;
//            LINALG::Matrix<5,1> tmp51;
//            tmp51.Multiply(PlAniso,eta_trial_vec_new);
//            tmp11.MultiplyTN(eta_trial_vec_new,tmp51);
//            AnormEtatrial_new=sqrt(tmp11(0,0));
//            double AnormEta_new=0.; // sqrt( eta_tr : A : eta_tr )
//            tmp51.Multiply(PlAniso,eta_vec_new);
//            tmp11.MultiplyTN(eta_vec_new,tmp51);
//            AnormEta_new=sqrt(tmp11(0,0));
//
//
//            LINALG::Matrix<5,1> Aetatrial_new;
//            Aetatrial_new.Multiply(PlAniso,eta_trial_vec_new);
//            tmp1.Clear();
//            tmp1(0,0) = 2./3.*Aetatrial_new(0) -1./3.*Aetatrial_new(1);
//            tmp1(1,1) =-1./3.*Aetatrial_new(0) +2./3.*Aetatrial_new(1);
//            tmp1(2,2) =-1./3.*Aetatrial_new(0) -1./3.*Aetatrial_new(1);
//            tmp1(0,1) = Aetatrial_new(2)/2.;
//            tmp1(1,0) = Aetatrial_new(2)/2.;
//            tmp1(1,2) = Aetatrial_new(3)/2.;
//            tmp1(2,1) = Aetatrial_new(3)/2.;
//            tmp1(0,2) = Aetatrial_new(4)/2.;
//            tmp1(2,0) = Aetatrial_new(4)/2.;
//            double EucNormAetaTrial_new = tmp1.Norm2();
//
//
//            tmp1.Clear();
//            tmp1(0,0) = 2./3.*Aeta_new(0)-1./3.*Aeta_new(1);
//            tmp1(1,1) =-1./3.*Aeta_new(0)+2./3.*Aeta_new(1);
//            tmp1(2,2) =-1./3.*Aeta_new(0)-1./3.*Aeta_new(1);
//            tmp1(0,1) = Aeta_new(2)/2.;
//            tmp1(1,0) = Aeta_new(2)/2.;
//            tmp1(1,2) = Aeta_new(3)/2.;
//            tmp1(2,1) = Aeta_new(3)/2.;
//            tmp1(0,2) = Aeta_new(4)/2.;
//            tmp1(2,0) = Aeta_new(4)/2.;
//            double EucNormAeta_new = tmp1.Norm2();
//            double DpAeta_new=0.;
//            for (int ii=0; ii<3; ii++)
//              for (int j=0; j<3; j++)
//                DpAeta_new += DeltaAlphaK_new(ii,j)*tmp1(ii,j);
//            double deltaAlphaI_new=0.;
//            if (DpAe<0.)
//              deltaAlphaI_new=-sqrt(2./3.)*AnormEta_new*DpAeta_new/EucNormAeta_new/EucNormAeta_new;
//
//
////            double deltaAlphaI_new=0.;
////            if (Dissipation>0.)
////              deltaAlphaI_new = sqrt(2./3.)*absDeltaAlphaK_new/EucNormAetaTrial_new*AnormEtatrial_new;
//
//
//            double Ypl_new=0.;
//            if (DpAe<0.)
//              Ypl_new = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+ deltaAlphaI_new)))
//                  + isohard*(last_alpha_isotropic_->at(gp)(0,0)+ deltaAlphaI_new) +inityield);
//            else
//              Ypl_new = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0))))
//                  + isohard*(last_alpha_isotropic_->at(gp)(0,0)) +inityield);
//
//
//            LINALG::Matrix<5,1> CplSym_new(true);
//            CplSym_new.Update(1.-Ypl_new/AnormEtatrial_new,Aeta_new,1.);
//            CplSym_new.Update(Ypl_new/AnormEtatrial_new*cpl,DeltaAlphaK_vec_new,1.);
//            //            CplSym_new.Scale(pow(AnormEtatrial_new,stab_s));
//            //
//            LINALG::Matrix<numdofperelement_,1> FintGp_new;
//            FintGp_new.MultiplyTN(detJ_w, bop, pk2_stress_new, 0.0);
//            //
//            LINALG::Matrix<3,1> SpEq_new;
//            SpEq_new(0) = -mDLp_last_iter_->at(gp)(5) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec_new(0)*eta_vec_new(2)
//                    +DeltaAlphaK_vec_new(2)*eta_vec_new(1)
//                    +DeltaAlphaK_vec_new(4)*eta_vec_new(3)
//                    -DeltaAlphaK_vec_new(2)*eta_vec_new(0)
//                    -DeltaAlphaK_vec_new(1)*eta_vec_new(2)
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(4) );
//            SpEq_new(1) = -mDLp_last_iter_->at(gp)(6) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec_new(2)*eta_vec_new(4)
//                    +DeltaAlphaK_vec_new(1)*eta_vec_new(3)*2.
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(0)
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(1)*2.
//                    -DeltaAlphaK_vec_new(4)*eta_vec_new(2)
//                    +DeltaAlphaK_vec_new(0)*eta_vec_new(3) );
//            SpEq_new(2) = -mDLp_last_iter_->at(gp)(7) -PlSpinEta/inityield *
//                ( DeltaAlphaK_vec_new(0)*eta_vec_new(4)*2.
//                    +DeltaAlphaK_vec_new(2)*eta_vec_new(3)
//                    -DeltaAlphaK_vec_new(4)*eta_vec_new(0)*2.
//                    -DeltaAlphaK_vec_new(4)*eta_vec_new(1)
//                    -DeltaAlphaK_vec_new(3)*eta_vec_new(2)
//                    +DeltaAlphaK_vec_new(1)*eta_vec_new(4) );
//
//            std::cout << std::scientific;
//            std::cout.precision(8);
//            for (int a=0; a<3; a++)
//              for (int A=0; A<3; A++)
//              {
//                double ref=DFpiDbeta(VOIGT3X3NONSYM_[A][a],i);
//                double FD=(InvPlasticDefgrd_new(A,a)-InvPlasticDefgrd(A,a))/epsilon;
//                double relE=0.; if (ref!=0.) relE=abs((FD-ref)/ref);
//                std::cout << "aAi:"<<a<<A<<i<< "\tFD: " << FD << "\tref: " << ref
//                    << "\trelE: " << relE << std::endl;
//              }
//              mDLp_last_iter_->at(gp)(i) -= epsilon;
//          }
//          std::cout << std::endl;
////          dserror("stop");

          // build matrices kbb and kbd
          for (int i=0; i<5; i++)
            for (int j=0; j<numdofperelement_; j++)
              kbetad(i,j) = DcplSymDd(i,j);
          for (int i=0; i<3; i++)
            for (int j=0; j<numdofperelement_; j++)
              kbetad(i+5,j) = dSpEqDd(i,j);

          for (int i=0; i<5; i++)
            for (int j=0; j<8; j++)
              kbetabeta(i,j) = DcplSymDbeta(i,j);
          for (int i=0; i<3; i++)
            for (int j=0; j<8; j++)
              kbetabeta(i+5,j) = dSpEqDbeta(i,j);

          for (int i=0; i<5; i++)
            force_beta(i) = CplSym(i);
          for (int i=0; i<3; i++)
            force_beta(i+5) = SpEq(i);

          // **************************************************************
          // static condensation of inner variables
          // **************************************************************
          //inverse matrix block [k_beta beta]_ij
          LINALG::Matrix<8,8> InvKbetabeta;
          Epetra_SerialDenseMatrix Kbetabeta_epetra(8,8);
          for (int i=0; i<8; i++)
            for (int j=0; j<8; j++)
              Kbetabeta_epetra(i,j) = kbetabeta(i,j);
          // we need the inverse of K_beta beta
          Epetra_SerialDenseSolver solve_for_inverseKbb;
          solve_for_inverseKbb.SetMatrix(Kbetabeta_epetra);
          solve_for_inverseKbb.Invert();
          for (int i=0; i<8; i++)
            for (int j=0; j<8; j++)
              InvKbetabeta(i,j) = Kbetabeta_epetra(i,j);
          // store for recover step
          KbbInvHill_->at(gp) = InvKbetabeta;
          fbetaHill_->at(gp) = force_beta;
          KbdHill_->at(gp) = kbetad;

          LINALG::Matrix<numdofperelement_,8> KdbKbb; // temporary  Kdb.Kbb^-1
          KdbKbb.Multiply(kdbeta,InvKbetabeta);

          // "plastic displacement stiffness"
          // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
          if (stiffmatrix!=NULL) stiffmatrix->Multiply(-1.,KdbKbb,kbetad,1.);
          // "plastic internal force"
          // plFint = [K_db.K_bb^-1].f_b
          if (force!=NULL) force->Multiply(-1.,KdbKbb,force_beta,1.);
          }
          else if (mDLp_last_iter_->at(gp).NormInf()!=0.)
          {
            KbbInvHill_->at(gp).Clear();
            for (int i=0; i<8; i++) KbbInvHill_->at(gp)(i,i)=1.;
            KbdHill_->at(gp).Clear();
            fbetaHill_->at(gp).Update(mDLp_last_iter_->at(gp));

            // condensation to internal force vector
            if (force!=NULL)
              force->Multiply(-1.,kdbeta,fbetaHill_->at(gp),1.);
          }
        } // modification for plastic Gauss points
      else
      {
        KbbInvHill_->at(gp).Clear();
        fbetaHill_ ->at(gp).Clear();
        KbdHill_   ->at(gp).Clear();
      }

      // square of the residual L2 norm
      double residual_norm_sqare=fbetaHill_->at(gp)(0)*fbetaHill_->at(gp)(0)
                                +fbetaHill_->at(gp)(1)*fbetaHill_->at(gp)(1)
                                +(-fbetaHill_->at(gp)(0)-fbetaHill_->at(gp)(1))*(-fbetaHill_->at(gp)(0)-fbetaHill_->at(gp)(1));
      for (int i=2; i<8; i++) residual_norm_sqare+=fbetaHill_->at(gp)(i)*fbetaHill_->at(gp)(i)*2.;
      params.get<double>("Lp_residual_square")+=residual_norm_sqare;

    } // gp loop

    // communicate unconverged active set to time integration
    if (converged_active_set==false)
      params.set("unconverged_active_set",true);

    return;
}

/*----------------------------------------------------------------------*
 | internal force, linear stiffness and mass                seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::lin_stiffmass(
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
  // no longer compatible with PlasticHyperElast
  dserror("DRT::ELEMENTS::So3_Plast<so3_ele,distype>::lin_stiffmass is out dated");

//  // update element geometry hex8, 3D: (8x3)
//  LINALG::Matrix<nen_,nsd_> xrefe;  // X, material coord. of element
//  LINALG::Matrix<nen_,nsd_> xcurr;  // x, current  coord. of element
//
//  DRT::Node** nodes = Nodes();
//  for (int i=0; i<nen_; ++i)
//  {
//    const double* x = nodes[i]->X();
//    xrefe(i,0) = x[0];
//    xrefe(i,1) = x[1];
//    xrefe(i,2) = x[2];
//
//    xcurr(i,0) = xrefe(i,0) + disp[i*numdofpernode_+0];
//    xcurr(i,1) = xrefe(i,1) + disp[i*numdofpernode_+1];
//    xcurr(i,2) = xrefe(i,2) + disp[i*numdofpernode_+2];
//  }
//
//  // we need the (residual) displacement -- current increment of displacement
//  LINALG::Matrix<numdofperelement_,1> res_d;
//  LINALG::Matrix<numdofperelement_,1> nodaldisp;
//  for (int i = 0; i<numdofperelement_; ++i)
//  {
//    nodaldisp(i) = disp[i];
//    res_d(i) = residual[i];
//  }
//
//  // compute derivatives N_XYZ at gp w.r.t. material coordinates
//  // by N_XYZ = J^-1 * N_rst
//  LINALG::Matrix<nsd_,nen_> N_XYZ;
//  // shape functions and their first derivatives
//  LINALG::Matrix<nen_,1> shapefunct;
//  LINALG::Matrix<nsd_,nen_> deriv;
//
//  // 3x3 LINALG matrix for temporary stuff
//  LINALG::Matrix<nsd_,nsd_> tmp1(false);
//  LINALG::Matrix<nsd_,nsd_> tmp2(false);
//
//  // get plastic material paramters
//  MAT::PlasticSemiSmooth* plmat = dynamic_cast<MAT::PlasticSemiSmooth*>(plasticmat_.access_private_ptr());
//  double kinhard=plmat->KinHard();
//  double isohard=plmat->IsoHard();
//  double expisohard=plmat->ExpIsoHard();
//  double infyield=plmat->InfYield();
//  double inityield=plmat->InitYield();
//
//  // converged active set
//  bool converged_active_set=true;
//
//  /* =========================================================================*/
//  /* ================================================= Loop over Gauss Points */
//  /* =========================================================================*/
//  for (int gp=0; gp<numgpt_; ++gp)
//  {
//    // shape functions (shapefunct) and their first derivatives (deriv)
//    DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
//    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp],deriv);
//
//    /* get the inverse of the Jacobian matrix which looks like:
//     **            [ x_,r  y_,r  z_,r ]^-1
//     **     J^-1 = [ x_,s  y_,s  z_,s ]
//     **            [ x_,t  y_,t  z_,t ]
//     */
//    // compute derivatives N_XYZ at gp w.r.t. material coordinates
//    // by N_XYZ = J^-1 * N_rst
//    N_XYZ.Multiply(invJ_[gp],deriv); // (6.21)
//    double detJ = detJ_[gp]; // (6.22)
//
//    // calculate the linear B-operator B_L = N_XYZ
//    LINALG::Matrix<numstr_,numdofperelement_> boplin;
//    CalculateBoplin(&boplin,&N_XYZ);
//
//
//    // approximate linearised strain tensor using common naming of strain vector
//    // linstrain={E11,E22,E33,2*E12,2*E23,2*E31}
//    LINALG::Matrix<numstr_,1> total_linstrain;
//    // E = epsilon_GL == epsilon_1
//    // build the linearised strain epsilon = B . d
//    total_linstrain.Multiply(boplin,nodaldisp);
//
//    // linear strain in matrix notation
//    LINALG::Matrix<nsd_,nsd_> total_linstrain_matrix;
//    total_linstrain_matrix(0,0) = total_linstrain(0);
//    total_linstrain_matrix(0,1) = 0.5*total_linstrain(3);
//    total_linstrain_matrix(0,2) = 0.5*total_linstrain(5);
//    total_linstrain_matrix(1,0) = total_linstrain_matrix(0,1);
//    total_linstrain_matrix(1,1) = total_linstrain(1);
//    total_linstrain_matrix(1,2) = 0.5*total_linstrain(4);
//    total_linstrain_matrix(2,0) = total_linstrain_matrix(0,2);
//    total_linstrain_matrix(2,1) = total_linstrain_matrix(1,2);
//    total_linstrain_matrix(2,2) = total_linstrain(2);
//
//    // return gp strains (only in case of stress/strain output)
//    switch (iostrain)
//    {
//    // in the linear realm all these strain measures are equal
//    case INPAR::STR::strain_gl:
//    case INPAR::STR::strain_ea:
//    {
//      if (elestrain == NULL) dserror("strain data not available");
//      for (int i = 0; i < 3; ++i)
//        (*elestrain)(gp,i) = total_linstrain(i);
//      for (int i = 3; i < 6; ++i)
//        (*elestrain)(gp,i) = 0.5 * total_linstrain(i);
//    }
//    break;
//    case INPAR::STR::strain_none:
//      break;
//    default:
//    {
//      dserror("requested strain type not available");
//      break;
//    }
//    }
//
//    // recover condensed variables from last iteration step *********
//    if (res_d.NormInf()!=0.)
//    {
//      // first part
//      LINALG::Matrix<5,1> tmp51;
//      tmp51.Multiply(KbbInv_->at(gp),fbeta_->at(gp));
//      DalphaK_last_iter_->at(gp).Update(-1.,tmp51,1.);
//
//      // second part
//      LINALG::Matrix<5,numdofperelement_> tmp524;
//      tmp524.Multiply(KbbInv_->at(gp),Kbd_->at(gp));
//      tmp51.Multiply(tmp524,res_d);
//      DalphaK_last_iter_->at(gp).Update(-1.,tmp51,1.);
//    }// end of recover **********************************************
//
//    // current plastic flow increment
//    LINALG::Matrix<nsd_,nsd_> DeltaAlphaK(false);
//    DeltaAlphaK(0,0) = DalphaK_last_iter_->at(gp)(0);
//    DeltaAlphaK(1,1) = DalphaK_last_iter_->at(gp)(1);
//    DeltaAlphaK(2,2) = -1.0*(DalphaK_last_iter_->at(gp)(0)+DalphaK_last_iter_->at(gp)(1));
//    DeltaAlphaK(0,1) = DalphaK_last_iter_->at(gp)(2);
//    DeltaAlphaK(1,0) = DalphaK_last_iter_->at(gp)(2);
//    DeltaAlphaK(1,2) = DalphaK_last_iter_->at(gp)(3);
//    DeltaAlphaK(2,1) = DalphaK_last_iter_->at(gp)(3);
//    DeltaAlphaK(0,2) = DalphaK_last_iter_->at(gp)(4);
//    DeltaAlphaK(2,0) = DalphaK_last_iter_->at(gp)(4);
//
//    // additive elastic/plastic split
//    LINALG::Matrix<nsd_,nsd_> elastic_linstrain_matrix(total_linstrain_matrix);
//    elastic_linstrain_matrix.Update(-1.,last_alpha_kinematic_->at(gp),1.);
//    elastic_linstrain_matrix.Update(-1.,DeltaAlphaK,1.);
//
//    // elastic strain in Voigt notation
//    LINALG::Matrix<numstr_,1> elastic_linstrain;
//    elastic_linstrain(0) = elastic_linstrain_matrix(0,0);
//    elastic_linstrain(1) = elastic_linstrain_matrix(1,1);
//    elastic_linstrain(2) = elastic_linstrain_matrix(2,2);
//    elastic_linstrain(3) = 2.*elastic_linstrain_matrix(0,1);
//    elastic_linstrain(4) = 2.*elastic_linstrain_matrix(1,2);
//    elastic_linstrain(5) = 2.*elastic_linstrain_matrix(0,2);
//
//    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
//     ** For semi-smooth plasticity only isotropic materials are valid.
//     ** Only SVK material for linear kinematics!
//     */
//    if (Material()->MaterialType()!=INPAR::MAT::m_stvenant)
//      dserror("so3_Plast works only with svk material for linear kinematics, sorry");
//    LINALG::Matrix<numstr_,numstr_> cmat(true);
//    LINALG::Matrix<numstr_,1> stress(true);
//    params.set<int>("gp",gp);
//    params.set<int>("eleID",Id());
//    EvaluateELasticMaterial(NULL,&elastic_linstrain,params,&stress,&cmat);
//    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc
//
//    // return gp stresses
//    switch (iostress)
//    {
//    // linear kinematics --> equal stress
//    case INPAR::STR::stress_2pk:
//    case INPAR::STR::stress_cauchy:
//    {
//      if (elestress == NULL) dserror("stress data not available");
//      for (int i = 0; i < numstr_; ++i)
//        (*elestress)(gp,i) = stress(i);
//    }
//    break;
//    case INPAR::STR::stress_none:
//      break;
//    default:
//    {
//      dserror("requested stress type not available");
//      break;
//    }
//    }
//
//    // stress in matrix form
//    LINALG::Matrix<3,3> stress_matrix;
//    stress_matrix(0,0) = stress(0);
//    stress_matrix(1,1) = stress(1);
//    stress_matrix(2,2) = stress(2);
//    stress_matrix(0,1) = stress(3);
//    stress_matrix(1,0) = stress_matrix(0,1);
//    stress_matrix(1,2) = stress(4);
//    stress_matrix(2,1) = stress_matrix(1,2);
//    stress_matrix(0,2) = stress(5);
//    stress_matrix(2,0) = stress_matrix(0,2);
//
//    // eta
//    LINALG::Matrix<3,3> eta(stress_matrix);
//    for (int i=0; i<3; i++)
//      eta(i,i) -= 1./3.*(stress_matrix(0,0) + stress_matrix(1,1) + stress_matrix(2,2));
//    eta.Update(-2./3.*kinhard,last_alpha_kinematic_->at(gp),1.);
//    eta.Update(-2./3.*kinhard,DeltaAlphaK,1.);
//    LINALG::Matrix<5,1> eta_vec(false);
//    eta_vec(0) = eta(0,0);
//    eta_vec(1) = eta(1,1);
//    eta_vec(2) = 0.5*(eta(0,1)+eta(1,0));
//    eta_vec(3) = 0.5*(eta(2,1)+eta(1,2));
//    eta_vec(4) = 0.5*(eta(0,2)+eta(2,0));
//
//    // eta_trial
//    LINALG::Matrix<3,3> eta_trial(eta);
//    eta_trial.Update(cpl_,DeltaAlphaK,1.0);
//    LINALG::Matrix<5,1> eta_trial_vec(false);
//    eta_trial_vec(0) = eta_trial(0,0);
//    eta_trial_vec(1) = eta_trial(1,1);
//    eta_trial_vec(2) = 0.5*(eta_trial(0,1)+eta_trial(1,0));
//    eta_trial_vec(3) = 0.5*(eta_trial(2,1)+eta_trial(1,2));
//    eta_trial_vec(4) = 0.5*(eta_trial(0,2)+eta_trial(2,0));
//
//    // absolute values
//    double absetatrial=eta_trial.Norm2();
//    double absDeltaAlphaK=DeltaAlphaK.Norm2();
//    double abseta=eta.Norm2();
//    double Dissipation=0.;
//    for (int i=0; i<3; i++)
//      for (int j=0; j<3; j++)
//        Dissipation-=eta(i,j)*DeltaAlphaK(i,j);
//
//    // current yield stress equivalent (yield stress scaled by sqrt(2/3))
//    double Ypl=0.;
//    if (Dissipation>0.)
//      Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+ sqrt(2./3.)*absDeltaAlphaK)))
//          + isohard*(last_alpha_isotropic_->at(gp)(0,0)+ sqrt(2./3.)*absDeltaAlphaK) +inityield);
//    else
//      Ypl = sqrt(2./3.) * ((infyield - inityield)*(1.-exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0))))
//          + isohard*(last_alpha_isotropic_->at(gp)(0,0)) +inityield);
//
//    // check activity state
//    // inactive
//    if (Ypl<absetatrial)
//    {
//      if (activity_state_->at(gp)==false) // gp switches state
//        converged_active_set = false;
//      activity_state_->at(gp) = true;
//    }
//    // active
//    else
//    {
//      if (activity_state_->at(gp)==true) // gp switches state
//        converged_active_set = false;
//      activity_state_->at(gp) = false;
//    }
//
//    // integrate elastic internal force vector **************************
//    double detJ_w = detJ*intpoints_.Weight(gp);
//
//    // update internal force vector
//    if (force != NULL)
//    {
//      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
//      force->MultiplyTN(detJ_w, boplin, stress, 1.0);
//    }
//    // end integrate elastic internal force vector **********************
//
//    // update/integrate `elastic' and `initial-displacement' stiffness matrix
//    if (stiffmatrix != NULL)
//    {
//      // keu = keu + (B^T . C . B) * detJ * w(gp)
//      LINALG::Matrix<6,numdofperelement_> cb;
//      cb.Multiply(cmat,boplin);
//      stiffmatrix->MultiplyTN(detJ_w,boplin,cb,1.0);
//    }
//
//    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
//    {
//      double density = Material()->Density();
//      // integrate consistent mass matrix
//      const double factor = detJ_w * density;
//      double ifactor, massfactor;
//      for (int inod=0; inod<nen_; ++inod)
//      {
//        ifactor = shapefunct(inod) * factor;
//        for (int jnod=0; jnod<nen_; ++jnod)
//        {
//          massfactor = shapefunct(inod) * ifactor;     // intermediate factor
//          (*massmatrix)(nsd_*inod+0,nsd_*jnod+0) += massfactor;
//          (*massmatrix)(nsd_*inod+1,nsd_*jnod+1) += massfactor;
//          (*massmatrix)(nsd_*inod+2,nsd_*jnod+2) += massfactor;
//        }
//      }
//
//    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
//
//    // plastic modifications
//    if (stiffmatrix!=NULL || force!=NULL)
//    {
//      // variables needed for condensation and calculated seperately for active and inactive Gauss points
//      LINALG::Matrix<5,numdofperelement_> kbetad(true);
//      LINALG::Matrix<5,5> kbetabeta(true);
//      LINALG::Matrix<5,1> force_beta(true);
//      LINALG::Matrix<numdofperelement_,5> kdbeta;
//
//      // 5x5 identitiy matrix
//      LINALG::Matrix<5,5> id5(true);
//      for (int i=0; i<5; i++)
//        id5(i,i)=1.;
//
//      // Due to stabilization, we have to treat inactive nodes as well
//      if (activity_state_->at(gp)==true || Dissipation>0.)
//      {
//        // damping parameter apl
//        double apl=1.;
//        if (Ypl/abseta<1.) apl=Ypl/abseta;
//
//        // matrix to eliminate the 6th stress variable as the stress is traceless
//        LINALG::Matrix<6,5> DDalphakDbetaVM(true);
//        DDalphakDbetaVM(0,0) = 1.;  // d Dalphak_0,0 d beta_0
//        DDalphakDbetaVM(1,1) = 1.;  // d Dalphak_1,1 d beta_1
//        DDalphakDbetaVM(2,0) = -1.; // d Dalphak_2,2 d beta_0
//        DDalphakDbetaVM(2,1) = -1.; // d Dalphak_2,2 d beta_1
//        DDalphakDbetaVM(3,2) = 2.;  // d Dalphak_0,1 d beta_2
//        DDalphakDbetaVM(4,3) = 2.;  // d Dalphak_1,2 d beta_3
//        DDalphakDbetaVM(5,4) = 2.;  // d Dalphak_0,2 d beta_4
//
//        // [K_d beta]_ij ***************************************************
//        LINALG::Matrix<numstr_,5> dsigmadbeta;
//        dsigmadbeta.Multiply(-1.,cmat,DDalphakDbetaVM);
//        kdbeta.MultiplyTN(detJ_w,boplin,dsigmadbeta);
//
//        if (activity_state_->at(gp)==true)
//        {
//          // communicate number of active plastic gauss points back to time integration
//          // don't sum up for ghost elements
//          if (MyPID == so3_ele::Owner()) params.get<int>("number_active_plastic_gp")++;
//
//          // [K_beta d]_ij ***************************************************
//
//          LINALG::Matrix<numstr_,numdofperelement_> dSigmadd;
//          dSigmadd.Multiply(cmat,boplin);
//
//          LINALG::Matrix<5,numdofperelement_> detadd(false);
//          for (int i=0; i<5; i++)
//            for (int j=0; j<numdofperelement_; j++)
//            {
//              // diagonal entries
//              if (i==0 || i==1)
//                detadd(i,j) = dSigmadd(VOIGT3X3SYM_[i][i],j)
//                - 1./3.*(  dSigmadd(0,j)
//                    +dSigmadd(1,j)
//                    +dSigmadd(2,j)  );
//              if (i==2)
//                detadd(i,j) = dSigmadd(3,j);
//              if (i==3)
//                detadd(i,j) = dSigmadd(4,j);
//              if (i==4)
//                detadd(i,j) = dSigmadd(5,j);
//            }
//
//          LINALG::Matrix<numdofperelement_,1> dabs_eta_trial_dd(true);
//          for (int j=0;j<2; j++)
//            for (int i=0; i<numdofperelement_; i++)
//              dabs_eta_trial_dd(i) += (2.*eta_trial_vec(j)+ eta_trial_vec((j+1)%2))*(detadd(j,i))/absetatrial;
//          for (int j=2; j<5; j++)
//            for (int i=0; i<numdofperelement_; i++)
//              dabs_eta_trial_dd(i) += (2.*eta_trial_vec(j))*(detadd(j,i))/absetatrial;
//
//          kbetad.Clear();
//          kbetad.Update(1.-1.*Ypl/absetatrial,detadd,1.);
//          kbetad.MultiplyNT(stab_s_*apl/absetatrial,eta_vec,dabs_eta_trial_dd,1.);
//          kbetad.MultiplyNT((1.-stab_s_)*Ypl/absetatrial/absetatrial,eta_trial_vec,dabs_eta_trial_dd,1.);
//
//          // store for recover step
//          Kbd_->at(gp) = kbetad;
//
//          // [K_ beta beta]_ij *************************************************
//
//          //                    d eta_ab
//          // detadbeta_abj = ---------------
//          //                    d beta_j
//          // in Voigt notation
//          // as eta is traceless, there are only 5 components
//          LINALG::Matrix<5,5> detadbeta;
//          for (int i=0; i<5; i++)
//            for (int j=0; j<5; j++)
//            {
//              // diagonal entries
//              if (i==0 || i==1)
//                detadbeta(i,j) = dsigmadbeta(i,j)
//                - 1./3.*(  dsigmadbeta(0,j)
//                    +dsigmadbeta(1,j)
//                    +dsigmadbeta(2,j)  )
//                    - 2./3.*kinhard * (i==j);
//              // off-diagonal entries
//              if (i==2)
//                detadbeta(i,j) = dsigmadbeta(3,j) - 2./3.*kinhard * (i==j);
//              if (i==3)
//                detadbeta(i,j) = dsigmadbeta(4,j) - 2./3.*kinhard * (i==j);
//              if (i==4)
//                detadbeta(i,j) = dsigmadbeta(5,j) - 2./3.*kinhard * (i==j);
//            }
//
//          LINALG::Matrix<5,1> dYpl_dbeta(true);
//          LINALG::Matrix<5,1> dabs_eta_trial_dbeta(true);
//          for (int j=0;j<2; j++)
//          {
//            if (Dissipation>0.) dYpl_dbeta(j) += 2./3. * (2.*DalphaK_last_iter_->at(gp)(j) + DalphaK_last_iter_->at(gp)((j+1)%2))/absDeltaAlphaK;
//            for (int i=0; i<5; i++)
//              dabs_eta_trial_dbeta(i) += (2.*eta_trial_vec(j)+ eta_trial_vec((j+1)%2))*(detadbeta(j,i)+id5(i,j)*cpl_)/absetatrial;
//          }
//          for (int j=2; j<5; j++)
//          {
//            if (Dissipation>0.) dYpl_dbeta(j) += 2./3. * 2.*DalphaK_last_iter_->at(gp)(j)/absDeltaAlphaK;
//            for (int i=0; i<5; i++)
//              dabs_eta_trial_dbeta(i) += (2.*eta_trial_vec(j))*(detadbeta(j,i)+id5(i,j)*cpl_)/absetatrial;
//          }
//          dYpl_dbeta.Scale(isohard + expisohard*(infyield-inityield)*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+ sqrt(2./3.)*absDeltaAlphaK)));
//
//          // build kbb from all previous linearizations
//          kbetabeta.Clear();
//          kbetabeta.Update(1.-Ypl/absetatrial,detadbeta,1.);
//          kbetabeta.Update(-Ypl/absetatrial*cpl_,id5,1.);
//          kbetabeta.MultiplyNT(-1./absetatrial,eta_trial_vec,dYpl_dbeta,1.);
//          kbetabeta.MultiplyNT((1.-stab_s_)*Ypl/pow(absetatrial,2.),eta_trial_vec,dabs_eta_trial_dbeta,1.);
//          kbetabeta.MultiplyNT(stab_s_*apl/absetatrial,eta_vec,dabs_eta_trial_dbeta,1.);
//
//          // right hand side term
//          force_beta.Update(eta_vec);
//          force_beta.Update(-1.*Ypl/absetatrial,eta_trial_vec,1.);
//          // store for recover step
//          fbeta_->at(gp) = force_beta;
//        } // active gp
//
//        // inactive Gauss point with plastic history within this load/time step
//        else if (activity_state_->at(gp)==false && Dissipation>0.)
//        {
//          // the complementarity function is
//          // C^pl = - Ypl^s * cplparam_ * delta alpha^k
//
//          // Complementarity function independent from displacements
//          kbetad.Clear();
//          Kbd_->at(gp).Clear();
//
//          LINALG::Matrix<5,1> dYpl_dbeta(true);
//          for (int j=0;j<2; j++)
//            dYpl_dbeta(j) += 2./3. * (2.*DalphaK_last_iter_->at(gp)(j) + DalphaK_last_iter_->at(gp)((j+1)%2))/absDeltaAlphaK;
//          for (int j=2; j<5; j++)
//            dYpl_dbeta(j) += 2./3. * 2.*DalphaK_last_iter_->at(gp)(j)/absDeltaAlphaK;
//          dYpl_dbeta.Scale(isohard + expisohard*(infyield-inityield)*exp(-expisohard*(last_alpha_isotropic_->at(gp)(0,0)+ sqrt(2./3.)*absDeltaAlphaK)));
//
//          kbetabeta.Update(-1.*cpl_,id5,0.);
//          kbetabeta.MultiplyNT(-1.*stab_s_/Ypl * cpl_,DalphaK_last_iter_->at(gp),dYpl_dbeta,1.);
//
//          // right hand side term
//          force_beta.Update(-1.*cpl_,DalphaK_last_iter_->at(gp),0.);
//          fbeta_->at(gp).Update(-1.*cpl_,DalphaK_last_iter_->at(gp),0.);
//        }
//      } // active or dissipation >0.
//
//      else
//      {
//        // Complementarity function independent from displacements
//        kbetad.Clear();
//        Kbd_->at(gp).Clear();
//        kdbeta.Clear();
//        kbetabeta.Update(-cpl_,id5);
//
//        // right hand side term
//        force_beta.Update(-cpl_,DalphaK_last_iter_->at(gp),0.);
//        fbeta_->at(gp).Update(-cpl_,DalphaK_last_iter_->at(gp),0.);
//      }
//
//      // **************************************************************
//      // static condensation of inner variables
//      // **************************************************************
//      //inverse matrix block [k_beta beta]_ij
//      LINALG::Matrix<5,5> InvKbetabeta;
//      Epetra_SerialDenseMatrix Kbetabeta_epetra(5,5);
//      for (int i=0; i<5; i++)
//        for (int j=0; j<5; j++)
//          Kbetabeta_epetra(i,j) = kbetabeta(i,j);
//      // we need the inverse of K_beta beta
//      Epetra_SerialDenseSolver solve_for_inverseKbb;
//      solve_for_inverseKbb.SetMatrix(Kbetabeta_epetra);
//      solve_for_inverseKbb.Invert();
//      for (int i=0; i<5; i++)
//        for (int j=0; j<5; j++)
//          InvKbetabeta(i,j) = Kbetabeta_epetra(i,j);
//      // store for recover step
//      KbbInv_->at(gp) = InvKbetabeta;
//
//      LINALG::Matrix<numdofperelement_,5> KdbKbb; // temporary  Kdb.Kbb^-1
//      KdbKbb.Multiply(kdbeta,InvKbetabeta);
//
//      // "plastic displacement stiffness"
//      // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
//      if (stiffmatrix!=NULL) stiffmatrix->Multiply(-1.,KdbKbb,kbetad,1.);
//
//      // "plastic internal force"
//      // plFint = [K_db.K_bb^-1].f_b
//      if (force!=NULL) force->Multiply(-1.,KdbKbb,force_beta,1.);
//
//    } // if (stiffmatrix!=NULL || force!=NULL)
//  } // gp loop
//
//  // communicate unconverged active set to time integration
//  if (converged_active_set==false)
//    params.set("unconverged_active_set",true);

  return;
}

/*----------------------------------------------------------------------*
 |  update plastic deformation for nonlinear kinematics     seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::UpdatePlasticDeformation_nln()
{
  // loop over all Gauss points
  for (int gp=0; gp<numgpt_; gp++)
  {
    if(activity_state_->at(gp)==true)
    {
    LINALG::Matrix<3,3> Dalphak;
    Dalphak(0,0) = DalphaK_last_iter_->at(gp)(0);
    Dalphak(1,1) = DalphaK_last_iter_->at(gp)(1);
    Dalphak(2,2) = -1.0*(DalphaK_last_iter_->at(gp)(0)+DalphaK_last_iter_->at(gp)(1));
    Dalphak(0,1) = DalphaK_last_iter_->at(gp)(2);
    Dalphak(1,0) = DalphaK_last_iter_->at(gp)(2);
    Dalphak(1,2) = DalphaK_last_iter_->at(gp)(3);
    Dalphak(2,1) = DalphaK_last_iter_->at(gp)(3);
    Dalphak(0,2) = DalphaK_last_iter_->at(gp)(4);
    Dalphak(2,0) = DalphaK_last_iter_->at(gp)(4);
      double absDalphak=0.;
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          absDalphak += Dalphak(i,j) * Dalphak(i,j);
      absDalphak = sqrt(absDalphak);

      //cout << "at gp: " << gp << ": delta alpha isotropic: " << absDalphak << endl;
      // evolution equation for isotropic hardening variable
      // alpha^i,n+1 = alpha^i,n + || delta alpha^k ||
      last_alpha_isotropic_->at(gp)(0) = last_alpha_isotropic_->at(gp)(0) + sqrt(2./3.)*absDalphak;

      // evolution equation for kinematic hardening variable
      // alpha^k,n+1 = alpha^k,n + delta alpha^k
      last_alpha_kinematic_->at(gp).Update(1.,Dalphak,1.);

      // evolution equation for plastic deformation gradient
      // F^p,n+1 = exp(-delta alpha^k) * Fp,n
      LINALG::Matrix<3,3> tmp(Dalphak);
      tmp.Scale(-1.);
      MatrixExponential3x3(tmp);
      LINALG::Matrix<3,3> FpLast(last_plastic_defgrd_inverse_->at(gp));
      FpLast.Invert();
      last_plastic_defgrd_inverse_->at(gp).Multiply(tmp,FpLast);
      last_plastic_defgrd_inverse_->at(gp).Invert();
    }

    DalphaK_last_timestep_->at(gp)=DalphaK_last_iter_->at(gp);
    KbbInv_->at(gp).Clear();
    Kbd_->at(gp).Clear();
    fbeta_->at(gp).Clear();
  }

  return;
}

/*----------------------------------------------------------------------*
 |  update plastic deformation for nonlinear kinematics     seitz 09/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::UpdatePlasticDeformationHill_nln()
{
  // loop over all Gauss points
  for (int gp=0; gp<numgpt_; gp++)
  {
    if(activity_state_->at(gp)==true)
    {
      LINALG::Matrix<nsd_,nsd_> DeltaAlphaK(false);
      LINALG::Matrix<nsd_,nsd_> mDLp;
      // -Dp
      DeltaAlphaK(0,0) = mDLp_last_iter_->at(gp)(0);
      DeltaAlphaK(1,1) = mDLp_last_iter_->at(gp)(1);
      DeltaAlphaK(2,2) = -1.0*(mDLp_last_iter_->at(gp)(0)+mDLp_last_iter_->at(gp)(1));
      DeltaAlphaK(0,1) = mDLp_last_iter_->at(gp)(2);
      DeltaAlphaK(1,0) = mDLp_last_iter_->at(gp)(2);
      DeltaAlphaK(1,2) = mDLp_last_iter_->at(gp)(3);
      DeltaAlphaK(2,1) = mDLp_last_iter_->at(gp)(3);
      DeltaAlphaK(0,2) = mDLp_last_iter_->at(gp)(4);
      DeltaAlphaK(2,0) = mDLp_last_iter_->at(gp)(4);
      mDLp.Update(DeltaAlphaK);
      // Wp
      LINALG::Matrix<3,3> Wp(true);
      Wp(0,1) -= mDLp_last_iter_->at(gp)(5);
      Wp(1,0) += mDLp_last_iter_->at(gp)(5);
      Wp(1,2) -= mDLp_last_iter_->at(gp)(6);
      Wp(2,1) += mDLp_last_iter_->at(gp)(6);
      Wp(0,2) -= mDLp_last_iter_->at(gp)(7);
      Wp(2,0) += mDLp_last_iter_->at(gp)(7);
      mDLp.Update(-1.,Wp,1.);

      //cout << "at gp: " << gp << ": delta alpha isotropic: " << absDalphak << endl;
      // evolution equation for isotropic hardening variable
      // alpha^i,n+1 = alpha^i,n + || delta alpha^k ||
      last_alpha_isotropic_->at(gp)(0) = last_alpha_isotropic_->at(gp)(0) + deltaAlphaI_->at(gp);

      // evolution equation for kinematic hardening variable
      // alpha^k,n+1 = alpha^k,n + delta alpha^k
      last_alpha_kinematic_->at(gp).Update(1.,DeltaAlphaK,1.);

      // evolution equation for plastic deformation gradient
      // F^p,n+1 = exp(-delta alpha^k) * Fp,n
      LINALG::Matrix<3,3> tmp(mDLp);
      tmp.Scale(-1.);
      MatrixExponential3x3(tmp);
      LINALG::Matrix<3,3> FpLast(last_plastic_defgrd_inverse_->at(gp));
      FpLast.Invert();
      last_plastic_defgrd_inverse_->at(gp).Multiply(tmp,FpLast);
      last_plastic_defgrd_inverse_->at(gp).Invert();
    }
    mDLp_last_timestep_->at(gp)=mDLp_last_iter_->at(gp);
    KbbInvHill_->at(gp).Clear();
    KbdHill_->at(gp).Clear();
    fbetaHill_->at(gp).Clear();
    deltaAlphaI_->at(gp)=0.;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  update plastic deformation for linear kinematics        seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::UpdatePlasticDeformation_lin()
{
  // loop over all Gauss points
  for (int gp=0; gp<numgpt_; gp++)
  {
    if (activity_state_->at(gp)==true)
    {
      LINALG::Matrix<3,3> Dalphak;
      Dalphak(0,0) = DalphaK_last_iter_->at(gp)(0);
      Dalphak(1,1) = DalphaK_last_iter_->at(gp)(1);
      Dalphak(2,2) = -1.0*(DalphaK_last_iter_->at(gp)(0)+DalphaK_last_iter_->at(gp)(1));
      Dalphak(0,1) = DalphaK_last_iter_->at(gp)(2);
      Dalphak(1,0) = DalphaK_last_iter_->at(gp)(2);
      Dalphak(1,2) = DalphaK_last_iter_->at(gp)(3);
      Dalphak(2,1) = DalphaK_last_iter_->at(gp)(3);
      Dalphak(0,2) = DalphaK_last_iter_->at(gp)(4);
      Dalphak(2,0) = DalphaK_last_iter_->at(gp)(4);

      double absDalphak=0.;
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          absDalphak += Dalphak(i,j) * Dalphak(i,j);
      absDalphak = sqrt(absDalphak);

      // evolution equation for isotropic hardening variable
      // alpha^i,n+1 = alpha^i,n + || delta alpha^k ||
      last_alpha_isotropic_->at(gp)(0) = last_alpha_isotropic_->at(gp)(0) + absDalphak;

      // evolution equation for kinematic hardening variable
      // alpha^k,n+1 = alpha^k,n + delta alpha^k
      last_alpha_kinematic_->at(gp).Update(1.,Dalphak,1.);
    }
    DalphaK_last_iter_->at(gp).Clear();
  }// gauss point loop
  return;
}


/*----------------------------------------------------------------------*
 |  matrix exponential                                      seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::MatrixExponential3x3( LINALG::Matrix<3,3>& MatrixInOut )
{
  double Norm=MatrixInOut.Norm2();
  // direct calculation for zero-matrix
  if (Norm==0.)
  {
    MatrixInOut.Clear();
    for (int i=0; i<3; i++)
      MatrixInOut(i,i)=1.;
    return;
  }

  // Calculation of matrix exponential via power series. This is usually
  // faster than by polar decomposition for matrices are close to zero.
  // For small plastic increments this is the case
  LINALG::Matrix<3,3> In(MatrixInOut);
  int n=0;
  int facn=1;
  MatrixInOut.Clear();
  for (int i=0; i<3; i++)
    MatrixInOut(i,i)=1.;
  LINALG::Matrix<3,3> tmp(MatrixInOut);
  LINALG::Matrix<3,3> tmp2(MatrixInOut);
  while (n<50 && tmp.Norm2()/facn>1.e-16)
  {
    n++;
    facn*=n;
    tmp.Multiply(tmp2,In);
    tmp2=tmp;
    MatrixInOut.Update(1./facn,tmp,1.);
  }
  if (n==50) dserror("matrix exponential unconverged in %i steps",n);

  return;
}

/*---------------------------------------------------------------------------*
 |  matrix exponential derivative of a symmetric matrix          seitz 07/13 |
 *---------------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::MatrixExponentialDerivativeSym3x3(const LINALG::Matrix<3,3> MatrixIn, LINALG::Matrix<6,6>& MatrixExpDeriv)
{
  double Norm=MatrixIn.Norm2();

  LINALG::Matrix<6,6> id4sharp(true);
  for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;

  // direct calculation for zero-matrix
  if (Norm==0.)
  {
    MatrixExpDeriv = id4sharp;
    return;
  }

  if(Norm<0.) // fixme
  {
    // see Souza-Neto: Computational Methods for plasticity, Box B.2.
    int nmax=0;
    int nIter=0;
    int nfac=1;
    LINALG::Matrix<3,3> tmp1;
    LINALG::Matrix<3,3> tmp2(true);
    for (int i=0; i<3; i++) tmp2(i,i)=1.;

    // all needed powers of X
    std::vector<LINALG::Matrix<3,3> > Xn;
    Xn.resize(0);
    Xn.push_back(tmp2);

    // all needed factorials
    std::vector<int> fac;
    fac.resize(0);
    fac.push_back(nfac);

    // compute nmax and Xn
    while (nIter<50 && tmp2.Norm2()/nfac>1.e-16)
    {
      nIter++;
      nfac *= nIter;
      fac.push_back(nfac);
      tmp1.Multiply(tmp2,MatrixIn);
      Xn.push_back(tmp1);
      tmp2=tmp1;
    }
    if (nIter==50) dserror("matrix exponential unconverged in %i steps",nIter);
    nmax=nIter;

    // compose derivative of matrix exponential (symmetric Voigt-notation)
    MatrixExpDeriv.Clear();
    for (int n=1; n<=nmax; n++)
    {
      for (int m=1; m<=n/2; m++)
        AddToSymMatrixExponentialDeriv(1./fac.at(n),Xn.at(m-1),Xn.at(n-m),MatrixExpDeriv);
      if (n%2==1)
        AddToSymMatrixExponentialDeriv(0.5/fac.at(n),Xn.at((n-1)/2),Xn.at((n-1)/2),MatrixExpDeriv);
    }
  }
  else
  {
    double EWtolerance=1.e-12;

    LINALG::Matrix<3,3> EV(MatrixIn);
    LINALG::Matrix<3,3> EW;
    LINALG::SYEV(EV,EW,EV);

    LINALG::Matrix<3,1> vec1;
    LINALG::Matrix<3,1> vec2;
    LINALG::Matrix<3,3> tmp1;
    LINALG::Matrix<3,3> tmp2;

    MatrixExpDeriv.Clear();
    // souza eq. (A.52)
    // note: EW stored in ascending order

    //  d X^2 / d X  =  1/2 * (  delta_jk X_lj + delta_il X_kj
    //                         + delta_jl X_ik + delta_kj X_il )
    //
    // y_i = log(x_i)
    // dy_i / dx_j = delta_ij 1/x_i

    LINALG::Matrix<3,3> id2(true);
    for (int i=0; i<3; i++)
      id2(i,i) =1.0 ;
    //  // --------------------------------- switch by number of equal eigenvalues

    if (abs(EW(0,0)-EW(1,1))<EWtolerance && abs(EW(1,1)-EW(2,2))<EWtolerance ) // ------------------ x_a == x_b == x_c
    {
      // calculate derivative
      MatrixExpDeriv = id4sharp;
      MatrixExpDeriv.Scale(exp(EW(0,0)));
    }

    else if ( ( abs(EW(0,0)-EW(1,1))<EWtolerance && abs(EW(1,1)-EW(2,2))>EWtolerance ) ||
        ( abs(EW(0,0)-EW(1,1))>EWtolerance && abs(EW(1,1)-EW(2,2))<EWtolerance )  ) // ---- x_a != x_b == x_c or x_a == x_b != x_c
    {
      // factors
      double s1=0.0;
      double s2=0.0;
      double s3=0.0;
      double s4=0.0;
      double s5=0.0;
      double s6=0.0;

      int a=0;
      int c=0;

      // switch which two EW are equal
      if ( abs(EW(0,0)-EW(1,1))<EWtolerance && abs(EW(1,1)-EW(2,2))>EWtolerance ) // ----------------------- x_a == x_b != x_c
      {
        a=2;
        c=0;
      }
      else if ( abs(EW(0,0)-EW(1,1))>EWtolerance && abs(EW(1,1)-EW(2,2))<EWtolerance) // ------------------ x_a != x_b == x_c
      {
        a=0;
        c=2;
      }
      else
        dserror("you should not be here");

      // in souza eq. (A.53):
      s1 = ( exp(EW(a,a)) - exp(EW(c,c)) ) / ( pow( EW(a,a) - EW(c,c),2.0 ) )  -  exp(EW(c,c)) / (EW(a,a)-EW(c,c));
      s2 = 2.0 * EW(c,c) * (exp(EW(a,a))-exp(EW(c,c)))/(pow(EW(a,a)-EW(c,c),2.0)) - (EW(a,a)+EW(c,c))/(EW(a,a)-EW(c,c)) * exp(EW(c,c));
      s3 = 2.0 * (exp(EW(a,a))-exp(EW(c,c)))/(pow(EW(a,a)-EW(c,c),3.0)) - (exp(EW(a,a)) + exp(EW(c,c)))/(pow(EW(a,a)-EW(c,c),2.0));
      s4 = EW(c,c)*s3;
      s5 = s4;
      s6 = EW(c,c)*EW(c,c) * s3;

      // calculate derivative
      MAT::AddToCmatDerivTensorSquare(MatrixExpDeriv,s1,MatrixIn,1.);
      MatrixExpDeriv.Update(-s2,id4sharp,1.);
      MAT::ElastSymTensorMultiply(MatrixExpDeriv,-1.*s3,MatrixIn,MatrixIn,1.);
      MAT::ElastSymTensorMultiply(MatrixExpDeriv,s4,MatrixIn,id2,1.);
      MAT::ElastSymTensorMultiply(MatrixExpDeriv,s5,id2,MatrixIn,1.);
      MAT::ElastSymTensorMultiply(MatrixExpDeriv,s6,id2,id2,1.);
    }

    else if ( abs(EW(0,0)-EW(1,1))>EWtolerance && abs(EW(1,1)-EW(2,2))>EWtolerance ) // ----------------- x_a != x_b != x_c
    {
      for (int a=0; a<3; a++) // loop over all eigenvalues
      {
        int b = (a+1)%3;
        int c = (a+2)%3;

        LINALG::Matrix<3,1> ea;
        LINALG::Matrix<3,1> eb;
        LINALG::Matrix<3,1> ec;
        for (int i=0; i<3; i++)
        {
          ea(i) = EV(i,a);
          eb(i) = EV(i,b);
          ec(i) = EV(i,c);
        }
        LINALG::Matrix<3,3> Ea;
        Ea.MultiplyNT(ea,ea);
        LINALG::Matrix<3,3> Eb;
        Eb.MultiplyNT(eb,eb);
        LINALG::Matrix<3,3> Ec;
        Ec.MultiplyNT(ec,ec);

        double fac = exp(EW(a,a)) / ( (EW(a,a)-EW(b,b)) * (EW(a,a)-EW(c,c)) );

        // + d X^2 / d X
        MAT::AddToCmatDerivTensorSquare(MatrixExpDeriv,fac,MatrixIn,1.);

        // - (x_b + x_c) I_s
        MatrixExpDeriv.Update(-1.*(EW(b,b)+EW(c,c))*fac,id4sharp,1.);

        // - [(x_a - x_b) + (x_a - x_c)] E_a \dyad E_a
        MAT::ElastSymTensorMultiply(MatrixExpDeriv,-1.*fac * ( (EW(a,a)-EW(b,b)) + (EW(a,a)-EW(c,c)) ),Ea,Ea,1.);


        // - (x_b - x_c) (E_b \dyad E_b)
        MAT::ElastSymTensorMultiply(MatrixExpDeriv,-1.*fac * (EW(b,b) - EW(c,c)),Eb,Eb,1.);

        // + (x_b - x_c) (E_c \dyad E_c)
        MAT::ElastSymTensorMultiply(MatrixExpDeriv,fac * (EW(b,b) - EW(c,c)),Ec,Ec,1.);

        // dy / dx_a E_a \dyad E_a
        MAT::ElastSymTensorMultiply(MatrixExpDeriv,exp(EW(a,a)),Ea,Ea,1.);
      } // end loop over all eigenvalues

    }

    else dserror("you should not be here.");
  }
  return;
}

/*---------------------------------------------------------------------------*
 |  matrix exponential derivative of a symmetric matrix          seitz 09/13 |
 *---------------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::MatrixExponentialDerivative3x3(const LINALG::Matrix<3,3> MatrixIn, LINALG::Matrix<9,9>& MatrixExpDeriv)
{
  // see Souza-Neto: Computational Methods for plasticity, Box B.2.
  int nmax=0;
  int nIter=0;
  int nfac=1;
  LINALG::Matrix<3,3> tmp1;
  LINALG::Matrix<3,3> tmp2(true);
  for (int i=0; i<3; i++) tmp2(i,i)=1.;

  // all needed powers of X
  std::vector<LINALG::Matrix<3,3> > Xn;
  Xn.resize(0);
  Xn.push_back(tmp2);

  // all needed factorials
  std::vector<int> fac;
  fac.resize(0);
  fac.push_back(nfac);

  // compute nmax and Xn
  while (nIter<50 && tmp2.Norm2()/nfac>1.e-16)
  {
    nIter++;
    nfac *= nIter;
    fac.push_back(nfac);
    tmp1.Multiply(tmp2,MatrixIn);
    Xn.push_back(tmp1);
    tmp2=tmp1;
  }
  if (nIter==50) dserror("matrix exponential unconverged in %i steps",nIter);
  nmax=nIter;

  // compose derivative of matrix exponential (non-symmetric Voigt-notation)
  MatrixExpDeriv.Clear();
  for (int n=1; n<=nmax; n++)
    for (int m=1; m<=n; m++)
      AddToMatrixExponentialDeriv(1./fac.at(n),Xn.at(m-1),Xn.at(n-m),MatrixExpDeriv);

  return;
}

/*---------------------------------------------------------------------------*
 |  add terms for matrix exponential derivative of a symmetric matrix        |
 | via power series                                              seitz 08/13 |
 *---------------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::AddToSymMatrixExponentialDeriv(const double fac,
    const LINALG::Matrix<3,3> A,const LINALG::Matrix<3,3> B, LINALG::Matrix<6,6>& Dexp)
{
  Dexp(0,0) += 2. * fac * A(0,0) * B(0,0);
  Dexp(0,3) += fac * ( A(0,0) * B(0,1) + A(0,1) * B(0,0));
  Dexp(0,5) += fac * ( A(0,0) * B(0,2) + A(0,2) * B(0,0));
  Dexp(0,1) += 2. * fac * A(0,1) * B(0,1);
  Dexp(0,4) += fac * ( A(0,1) * B(0,2) + A(0,2) * B(0,1));
  Dexp(0,2) += 2. * fac * A(0,2) * B(0,2);

  Dexp(3,0) += fac * (A(0,0) * B(0,1) + A(0,1) * B(0,0));
  Dexp(3,3) += 0.5 * fac * (A(0,0) * B(1,1) + 2. * A(0,1) * B(0,1) + A(1,1) * B(0,0));
  Dexp(3,5) += 0.5 * fac * (A(0,0) * B(1,2) + A(0,2) * B(0,1) + A(1,2) * B(0,0) + A(0,1) * B(0,2));
  Dexp(3,1) +=  fac * ( A(0,1) * B(1,1) +  A(1,1) * B(0,1));
  Dexp(3,4) += 0.5 * fac * (A(0,1) * B(1,2) + A(0,2) * B(1,1) + A(1,2) * B(0,1) + A(1,1) * B(0,2));
  Dexp(3,2) += fac * ( A(0,2) * B(1,2) + A(1,2) * B(0,2));

  Dexp(5,0) += fac * (A(0,0) * B(0,2) + A(0,2) * B(0,0));
  Dexp(5,3) += 0.5 * fac * (A(0,0) * B(1,2) + A(0,2) * B(0,1) + A(1,2) * B(0,0) + A(0,1) * B(0,2));
  Dexp(5,5) += 0.5 * fac * (A(0,0) * B(2,2) + 2. * A(0,2) * B(0,2) + A(2,2) * B(0,0));
  Dexp(5,1) += fac * (A(0,1) * B(1,2) + A(1,2) * B(0,1));
  Dexp(5,4) += 0.5 * fac * (A(0,1) * B(2,2) + A(0,2) * B(1,2) + A(2,2) * B(0,1) + A(1,2) * B(0,2));
  Dexp(5,2) += fac * ( A(0,2) * B(2,2) +  A(2,2) * B(0,2));

  Dexp(1,0) += 2. * fac * A(0,1) * B(0,1);
  Dexp(1,3) +=  fac * ( A(0,1) * B(1,1) + A(1,1) * B(0,1));
  Dexp(1,5) +=  fac * ( A(0,1) * B(1,2) +  A(1,2) * B(0,1));
  Dexp(1,1) += 2. * fac * A(1,1) * B(1,1);
  Dexp(1,4) +=  fac * ( A(1,1) * B(1,2) + A(1,2) * B(1,1));
  Dexp(1,2) += 2. * fac * A(1,2) * B(1,2);

  Dexp(4,0) +=  fac * ( A(0,1) * B(0,2) +  A(0,2) * B(0,1));
  Dexp(4,3) += 0.5 * fac * (A(0,1) * B(1,2) + A(0,2) * B(1,1) + A(1,2) * B(0,1) + A(1,1) * B(0,2));
  Dexp(4,5) += 0.5 * fac * (A(0,1) * B(2,2) + A(0,2) * B(1,2) + A(2,2) * B(0,1) + A(1,2) * B(0,2));
  Dexp(4,1) +=  fac * ( A(1,1) * B(1,2) + A(1,2) * B(1,1));
  Dexp(4,4) += 0.5 * fac * (A(1,1) * B(2,2) + 2. * A(1,2) * B(1,2) + A(2,2) * B(1,1));
  Dexp(4,2) += fac * ( A(1,2) * B(2,2) + A(2,2) * B(1,2));

  Dexp(2,0) += 2. * fac * A(0,2) * B(0,2);
  Dexp(2,3) += fac * (A(0,2) * B(1,2) +  A(1,2) * B(0,2));
  Dexp(2,5) += fac * ( A(0,2) * B(2,2) +  A(2,2) * B(0,2));
  Dexp(2,1) += 2. * fac * A(1,2) * B(1,2);
  Dexp(2,4) += fac * ( A(1,2) * B(2,2) + A(2,2) * B(1,2));
  Dexp(2,2) += 2. * fac * A(2,2) * B(2,2);

  return;
}

/*---------------------------------------------------------------------------*
 |  add terms for matrix exponential derivative of a symmetric matrix        |
 | via power series                                              seitz 09/13 |
 *---------------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::AddToMatrixExponentialDeriv(const double fac,
    const LINALG::Matrix<3,3> A,const LINALG::Matrix<3,3> B, LINALG::Matrix<9,9>& Dexp)
{
  Dexp(0,0) += fac * A(0,0) * B(0,0);
  Dexp(0,3) += fac * A(0,0) * B(1,0);
  Dexp(0,5) += fac * A(0,0) * B(2,0);
  Dexp(0,6) += fac * A(0,1) * B(0,0);
  Dexp(0,1) += fac * A(0,1) * B(1,0);
  Dexp(0,4) += fac * A(0,1) * B(2,0);
  Dexp(0,8) += fac * A(0,2) * B(0,0);
  Dexp(0,7) += fac * A(0,2) * B(1,0);
  Dexp(0,2) += fac * A(0,2) * B(2,0);

  Dexp(3,0) += fac * A(0,0) * B(0,1);
  Dexp(3,3) += fac * A(0,0) * B(1,1);
  Dexp(3,5) += fac * A(0,0) * B(2,1);
  Dexp(3,6) += fac * A(0,1) * B(0,1);
  Dexp(3,1) += fac * A(0,1) * B(1,1);
  Dexp(3,4) += fac * A(0,1) * B(2,1);
  Dexp(3,8) += fac * A(0,2) * B(0,1);
  Dexp(3,7) += fac * A(0,2) * B(1,1);
  Dexp(3,2) += fac * A(0,2) * B(2,1);

  Dexp(5,0) += fac * A(0,0) * B(0,2);
  Dexp(5,3) += fac * A(0,0) * B(1,2);
  Dexp(5,5) += fac * A(0,0) * B(2,2);
  Dexp(5,6) += fac * A(0,1) * B(0,2);
  Dexp(5,1) += fac * A(0,1) * B(1,2);
  Dexp(5,4) += fac * A(0,1) * B(2,2);
  Dexp(5,8) += fac * A(0,2) * B(0,2);
  Dexp(5,7) += fac * A(0,2) * B(1,2);
  Dexp(5,2) += fac * A(0,2) * B(2,2);

  Dexp(6,0) += fac * A(1,0) * B(0,0);
  Dexp(6,3) += fac * A(1,0) * B(1,0);
  Dexp(6,5) += fac * A(1,0) * B(2,0);
  Dexp(6,6) += fac * A(1,1) * B(0,0);
  Dexp(6,1) += fac * A(1,1) * B(1,0);
  Dexp(6,4) += fac * A(1,1) * B(2,0);
  Dexp(6,8) += fac * A(1,2) * B(0,0);
  Dexp(6,7) += fac * A(1,2) * B(1,0);
  Dexp(6,2) += fac * A(1,2) * B(2,0);

  Dexp(1,0) += fac * A(1,0) * B(0,1);
  Dexp(1,3) += fac * A(1,0) * B(1,1);
  Dexp(1,5) += fac * A(1,0) * B(2,1);
  Dexp(1,6) += fac * A(1,1) * B(0,1);
  Dexp(1,1) += fac * A(1,1) * B(1,1);
  Dexp(1,4) += fac * A(1,1) * B(2,1);
  Dexp(1,8) += fac * A(1,2) * B(0,1);
  Dexp(1,7) += fac * A(1,2) * B(1,1);
  Dexp(1,2) += fac * A(1,2) * B(2,1);

  Dexp(4,0) += fac * A(1,0) * B(0,2);
  Dexp(4,3) += fac * A(1,0) * B(1,2);
  Dexp(4,5) += fac * A(1,0) * B(2,2);
  Dexp(4,6) += fac * A(1,1) * B(0,2);
  Dexp(4,1) += fac * A(1,1) * B(1,2);
  Dexp(4,4) += fac * A(1,1) * B(2,2);
  Dexp(4,8) += fac * A(1,2) * B(0,2);
  Dexp(4,7) += fac * A(1,2) * B(1,2);
  Dexp(4,2) += fac * A(1,2) * B(2,2);

  Dexp(8,0) += fac * A(2,0) * B(0,0);
  Dexp(8,3) += fac * A(2,0) * B(1,0);
  Dexp(8,5) += fac * A(2,0) * B(2,0);
  Dexp(8,6) += fac * A(2,1) * B(0,0);
  Dexp(8,1) += fac * A(2,1) * B(1,0);
  Dexp(8,4) += fac * A(2,1) * B(2,0);
  Dexp(8,8) += fac * A(2,2) * B(0,0);
  Dexp(8,7) += fac * A(2,2) * B(1,0);
  Dexp(8,2) += fac * A(2,2) * B(2,0);

  Dexp(7,0) += fac * A(2,0) * B(0,1);
  Dexp(7,3) += fac * A(2,0) * B(1,1);
  Dexp(7,5) += fac * A(2,0) * B(2,1);
  Dexp(7,6) += fac * A(2,1) * B(0,1);
  Dexp(7,1) += fac * A(2,1) * B(1,1);
  Dexp(7,4) += fac * A(2,1) * B(2,1);
  Dexp(7,8) += fac * A(2,2) * B(0,1);
  Dexp(7,7) += fac * A(2,2) * B(1,1);
  Dexp(7,2) += fac * A(2,2) * B(2,1);

  Dexp(2,0) += fac * A(2,0) * B(0,2);
  Dexp(2,3) += fac * A(2,0) * B(1,2);
  Dexp(2,5) += fac * A(2,0) * B(2,2);
  Dexp(2,6) += fac * A(2,1) * B(0,2);
  Dexp(2,1) += fac * A(2,1) * B(1,2);
  Dexp(2,4) += fac * A(2,1) * B(2,2);
  Dexp(2,8) += fac * A(2,2) * B(0,2);
  Dexp(2,7) += fac * A(2,2) * B(1,2);
  Dexp(2,2) += fac * A(2,2) * B(2,2);
}
