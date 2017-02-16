/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast_evaluate.cpp
\maintainer Alexander Seitz
\brief
\level 2
\maintainer Alexander Seitz


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
  dsassert(kintype_==INPAR::STR::kinem_nonlinearTotLag,"only geometricallly nonlinear formluation for plasticity!");

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
  else if (action == "calc_struct_stifftemp")                     act =  calc_struct_stifftemp;
  else if (action=="calc_struct_energy")                          act =  calc_struct_energy;
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

    // initialise the vectors
    // Evaluate() is called the first time in StructureBaseAlgorithm: at this
    // stage the coupling field is not yet known. Pass coupling vectors filled
    // with zeros
    // the size of the vectors is the length of the location vector/nsd_
    // velocities for TSI problem
    std::vector<double> myvel(0);
    std::vector<double> mytempnp(0);
    std::vector<double> mytempres(0);
    if (tsi_)
    {
      if (discretization.HasState(0,"velocity"))
      {
        // get the velocities
        Teuchos::RCP<const Epetra_Vector> vel
        = discretization.GetState(0,"velocity");
         if (vel == Teuchos::null)
           dserror("Cannot get state vectors 'velocity'");
        // extract the velocities
        myvel.resize((la[0].lm_).size());
        DRT::UTILS::ExtractMyValues(*vel,myvel,la[0].lm_);
      }
      if (discretization.HasState(1,"temperature"))
      {
        Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1,"temperature");
        if (tempnp==Teuchos::null)
          dserror("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the dimension of the problem
        const int numdofpernode_thr = discretization.NumDof(1,Nodes()[0]);
        if (la[1].Size() != nen_*numdofpernode_thr)
          dserror("Location vector length for temperature does not match!");
        // extract the current temperatures
        mytempnp.resize( ( (la[0].lm_).size() )/nsd_, 0.0 );
        DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
      }
      if (discretization.HasState(1,"residual temperature"))
      {
        Teuchos::RCP<const Epetra_Vector> tempres = discretization.GetState(1,"residual temperature");
        if (tempres==Teuchos::null)
          dserror("Cannot get state vector 'tempres'");

        // the temperature field has only one dof per node, disregarded by the dimension of the problem
        const int numdofpernode_thr = discretization.NumDof(1,Nodes()[0]);
        if (la[1].Size() != nen_*numdofpernode_thr)
          dserror("Location vector length for temperature does not match!");
        // extract the current temperatures
        mytempres.resize( ( (la[0].lm_).size() )/nsd_, 0.0 );
        DRT::UTILS::ExtractMyValues(*tempres,mytempres,la[1].lm_);
      }
    }

    // default: geometrically non-linear analysis with Total Lagrangean approach
    nln_stiffmass(la[0].lm_,mydisp,myres,myvel,mytempnp,mytempres,&myemat,NULL,&elevec1,NULL,NULL,NULL,params,
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
    LINALG::Matrix<numdofperelement_,1> elevec3(elevec3_epetra.A(),true);

    // need current displacement and residual forces
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
    if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
    std::vector<double> mydisp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);
    std::vector<double> myres(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*res,myres,la[0].lm_);

    // initialise the vectors
    // Evaluate() is called the first time in StructureBaseAlgorithm: at this
    // stage the coupling field is not yet known. Pass coupling vectors filled
    // with zeros
    // the size of the vectors is the length of the location vector/nsd_
    // velocities for TSI problem
    std::vector<double> myvel(0);
    std::vector<double> mytempnp(0);
    std::vector<double> mytempres(0);
    if (tsi_)
    {
      if (discretization.HasState(0,"velocity"))
      {
        // get the velocities
        Teuchos::RCP<const Epetra_Vector> vel
        = discretization.GetState(0,"velocity");
         if (vel == Teuchos::null)
           dserror("Cannot get state vectors 'velocity'");
        // extract the velocities
        myvel.resize((la[0].lm_).size());
        DRT::UTILS::ExtractMyValues(*vel,myvel,la[0].lm_);
      }

      if (discretization.HasState(1,"temperature"))
      {
        Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1,"temperature");
        if (tempnp==Teuchos::null)
          dserror("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the dimension of the problem
        const int numdofpernode_thr = discretization.NumDof(1,Nodes()[0]);
        if (la[1].Size() != nen_*numdofpernode_thr)
          dserror("Location vector length for temperature does not match!");
        // extract the current temperatures
        mytempnp.resize( ( (la[0].lm_).size() )/nsd_, 0.0 );
        DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
      }
      if (discretization.HasState(1,"residual temperature"))
      {
        Teuchos::RCP<const Epetra_Vector> tempres = discretization.GetState(1,"residual temperature");
        if (tempres==Teuchos::null)
          dserror("Cannot get state vector 'tempres'");

        // the temperature field has only one dof per node, disregarded by the dimension of the problem
        const int numdofpernode_thr = discretization.NumDof(1,Nodes()[0]);
        if (la[1].Size() != nen_*numdofpernode_thr)
          dserror("Location vector length for temperature does not match!");
        // extract the current temperatures
        mytempres.resize( ( (la[0].lm_).size() )/nsd_, 0.0 );
        DRT::UTILS::ExtractMyValues(*tempres,mytempres,la[1].lm_);
      }
    }

    // default: geometrically non-linear analysis with Total Lagrangean approach
    nln_stiffmass(la[0].lm_,mydisp,myres,myvel,mytempnp,mytempres,matptr,NULL,&elevec1,&elevec3,NULL,NULL,params,
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
     // internal force without condensation
     LINALG::Matrix<numdofperelement_,1> elevec3(elevec3_epetra.A(),true);

     // initialise the vectors
     // Evaluate() is called the first time in StructureBaseAlgorithm: at this
     // stage the coupling field is not yet known. Pass coupling vectors filled
     // with zeros
     // the size of the vectors is the length of the location vector/nsd_
     // velocities for TSI problem
     std::vector<double> myvel(0);
     std::vector<double> mytempnp(0);
     std::vector<double> mytempres(0);
     if (tsi_)
     {
       if (discretization.HasState(0,"velocity"))
       {
         // get the velocities
         Teuchos::RCP<const Epetra_Vector> vel
         = discretization.GetState(0,"velocity");
          if (vel == Teuchos::null)
            dserror("Cannot get state vectors 'velocity'");
         // extract the velocities
         myvel.resize((la[0].lm_).size());
         DRT::UTILS::ExtractMyValues(*vel,myvel,la[0].lm_);
       }
       if (discretization.HasState(1,"temperature"))
       {
         Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1,"temperature");
         if (tempnp==Teuchos::null)
           dserror("Cannot get state vector 'tempnp'");

         // the temperature field has only one dof per node, disregarded by the dimension of the problem
         const int numdofpernode_thr = discretization.NumDof(1,Nodes()[0]);
         if (la[1].Size() != nen_*numdofpernode_thr)
           dserror("Location vector length for temperature does not match!");
         // extract the current temperatures
         mytempnp.resize( ( (la[0].lm_).size() )/nsd_, 0.0 );
         DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
       }
       if (discretization.HasState(1,"residual temperature"))
       {
         Teuchos::RCP<const Epetra_Vector> tempres = discretization.GetState(1,"residual temperature");
         if (tempres==Teuchos::null)
           dserror("Cannot get state vector 'tempres'");

         // the temperature field has only one dof per node, disregarded by the dimension of the problem
         const int numdofpernode_thr = discretization.NumDof(1,Nodes()[0]);
         if (la[1].Size() != nen_*numdofpernode_thr)
           dserror("Location vector length for temperature does not match!");
         // extract the current temperatures
         mytempres.resize( ( (la[0].lm_).size() )/nsd_, 0.0 );
         DRT::UTILS::ExtractMyValues(*tempres,mytempres,la[1].lm_);
       }
     }

     // default: geometrically non-linear analysis with Total Lagrangean approach
     nln_stiffmass(la[0].lm_,mydisp,myres,myvel,mytempnp,mytempres,&elemat1,&elemat2,&elevec1,&elevec3,NULL,NULL,params,
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

      // initialise the vectors
      // Evaluate() is called the first time in StructureBaseAlgorithm: at this
      // stage the coupling field is not yet known. Pass coupling vectors filled
      // with zeros
      // the size of the vectors is the length of the location vector/nsd_
      // velocities for TSI problem
      std::vector<double> myvel(0);
      std::vector<double> mytempnp(0);
      std::vector<double> mytempres(0);
      if (tsi_)
      {
        if (discretization.HasState(0,"velocity"))
        {
          // get the velocities
          Teuchos::RCP<const Epetra_Vector> vel
          = discretization.GetState(0,"velocity");
           if (vel == Teuchos::null)
             dserror("Cannot get state vectors 'velocity'");
          // extract the velocities
          myvel.resize((la[0].lm_).size());
          DRT::UTILS::ExtractMyValues(*vel,myvel,la[0].lm_);
        }
        if (discretization.HasState(1,"temperature"))
        {
          Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(1,"temperature");
          if (tempnp==Teuchos::null)
            dserror("Cannot get state vector 'tempnp'");

          // the temperature field has only one dof per node, disregarded by the dimension of the problem
          const int numdofpernode_thr = discretization.NumDof(1,Nodes()[0]);
          if (la[1].Size() != nen_*numdofpernode_thr)
            dserror("Location vector length for temperature does not match!");
          // extract the current temperatures
          mytempnp.resize( ( (la[0].lm_).size() )/nsd_, 0.0 );
          DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
        }
      }

      LINALG::Matrix<numgpt_post,numstr_> stress;
      LINALG::Matrix<numgpt_post,numstr_> strain;
      INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

      // default: geometrically non-linear analysis with Total Lagrangean approach
      nln_stiffmass(la[0].lm_,mydisp,myres,myvel,mytempnp,mytempres,NULL,NULL,NULL,NULL,&stress,&strain,params,
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
      KbbInv_[i].Scale(0.);
      Kbd_   [i].Scale(0.);
      fbeta_ [i].Scale(0.);
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

  //============================================================================
  case calc_struct_stifftemp: // here we want to build the K_dT block for monolithic TSI
  {
    // stiffness
    LINALG::Matrix<numdofperelement_,nen_> k_dT(elemat1_epetra.A(),true);

    // calculate matrix block
    nln_kdT_tsi(&k_dT,params);
  }
  break;

  case calc_struct_energy:
  {
    // need current displacement
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0,"displacement");
    std::vector<double> mydisp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);

    elevec1_epetra(0) = CalcIntEnergy(mydisp,params);
  }
  break;

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

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < nsd_)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = nsd_; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // find out whether we will use time curves and get the factors
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  std::vector<double> curvefacs(nsd_, 1.0);
  for (int i=0; i < nsd_; ++i)
  {
    const int curvenum = (curve) ? (*curve)[i] : -1;
    if (curvenum>=0 && usetime)
      curvefacs[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }


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
    const double fac = wgt_[gp] * detJ;
    // distribute/add over element load vector
    for(int dim=0; dim<nsd_; dim++) {
      // function evaluation
      const int functnum = (funct) ? (*funct)[dim] : -1;
      const double functfac
        = (functnum>0)
        ? DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,xrefegp.A(),time,NULL)
        : 1.0;
      const double dim_fac = (*onoff)[dim] * (*val)[dim] * fac * curvefacs[dim] * functfac;
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
    std::vector<double>& vel,                      // current velocities
    std::vector<double>&           temp,           // current temperatures
    std::vector<double>&           temp_res,       // current temperatures
    LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix, // element stiffness matrix
    LINALG::Matrix<numdofperelement_,numdofperelement_>* massmatrix,  // element mass matrix
    LINALG::Matrix<numdofperelement_,1>* force,                 // element internal force vector
    LINALG::Matrix<numdofperelement_,1>* force_str,          // structure force
    LINALG::Matrix<numgpt_post,numstr_>* elestress,   // stresses at GP
    LINALG::Matrix<numgpt_post,numstr_>* elestrain,   // strains at GP
    Teuchos::ParameterList&        params,         // algorithmic parameters e.g. time
    const INPAR::STR::StressType   iostress,  // stress output option
    const INPAR::STR::StrainType   iostrain,  // strain output option
    const int MyPID  // processor id
    )
{
  if (data_==Teuchos::null)
    if (params.isParameter("PlastSsnData"))
    data_=params.get<Teuchos::RCP<UTILS::PlastSsnData> >("PlastSsnData");

  // do the evaluation of tsi terms
  const bool eval_tsi = (temp.size()!=0);

  // update element geometry
  LINALG::Matrix<nen_,3> xrefe(false);      // X, material coord. of element
  LINALG::Matrix<nen_,3> xcurr(false);      // x, current  coord. of element
  LINALG::Matrix<nen_,3> xcurrrate(false);  // x, rate of current  coord. of element
  LINALG::Matrix<nen_,1> etemp(false);      // vector of the current element temperatures
  LINALG::Matrix<nen_,1> res_T(true);      // vector of the current element residual temperatures

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
    if (eval_tsi)
    {
      etemp(i)=temp[i];
      if (neas_!=0)
        if (temp_res.size()!=0)
          res_T(i)=temp_res[i];

      if (vel.size()!=0)
      {
        xcurrrate(i,0) = vel[i*numdofpernode_+0];
        xcurrrate(i,1) = vel[i*numdofpernode_+1];
        xcurrrate(i,2) = vel[i*numdofpernode_+2];
      }
    }
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
    LINALG::Matrix<nsd_,nen_> N_rst_0(false);
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
  {
    if (tsi_==true)
      dserror("TSI with so3Plast elements only with PlasticElastHyper material");
  }
  // get references from parameter list
  double& lp_inc = data_->pl_inc_;
  double& lp_res = data_->pl_res_;
  double& eas_inc= data_->eas_inc_;
  double& eas_res= data_->eas_res_;
  int& num_active_gp = data_->num_active_;
  INPAR::STR::PredEnum pred = data_->pred_type_;
  bool no_condensation = data_->no_pl_condensation_;

  // do not recover condensed variables if it is a TSI predictor step
  bool no_recovery=data_->no_recovery_;
  // time integration factor (for TSI)
  double theta=-1.;
  // time step size (for TSI)
  double dt=-1.;
  if (eval_tsi && (stiffmatrix!=NULL || force!=NULL))
  {
    // get time integration data
    theta=data_->scale_timint_;
    dt=data_->dt_;
    if (theta<=0 || dt<=0)
      dserror("time integration parameters not provided in element for TSI problem");
  }
  // check if we need to split the residuals (for Newton line search)
  // if true an additional global vector is assembled containing
  // the internal forces without the condensed EAS entries and the norm
  // of the EAS residual is calculated
  bool split_res = data_->split_res_;

  /* evaluation of EAS variables (which are constant for the following):
  ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
  ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
  ** -> T0^{-T}
  */
  std::vector<Epetra_SerialDenseMatrix>* M_GP = NULL;   // EAS matrix M at all GPs
  LINALG::SerialDenseMatrix M;      // EAS matrix M at current GP
  double detJ0;                     // detJ(origin)
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
    if (stiffmatrix!=NULL && !no_recovery)
    {
      // this is a line search step, i.e. the direction of the eas increments
      // has been calculated by a Newton step and now it is only scaled
      if (data_->ls_)
      {
        double alpha_ls=data_->alpha_ls_;
        // undo step
        alpha_eas_inc_->Scale(-1.);
        alpha_eas_->operator +=(*alpha_eas_inc_);
        // scale increment
        alpha_eas_inc_->Scale(-1.*alpha_ls);
        // add reduced increment
        alpha_eas_->operator +=(*alpha_eas_inc_);
      }
      else
      {
        // constant predictor
        if (pred == INPAR::STR::pred_constdis)
          *alpha_eas_=*alpha_eas_last_timestep_;

        // tangential predictor
        else if (pred == INPAR::STR::pred_tangdis || pred == INPAR::STR::pred_constvel || pred == INPAR::STR::pred_constacc)
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
            if (KaT_!=Teuchos::null) LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,nen_,1>(1.,feas_->A(),1.,KaT_->A(),res_T.A());
            LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,soh8p_easmild,1>(0.0,*alpha_eas_inc_,-1.0,*KaaInv_,*feas_);
            LINALG::DENSEFUNCTIONS::update<double,soh8p_easmild,1>(1.0,*alpha_eas_,1.0,*alpha_eas_inc_);
            break;
          case soh8p_easfull:
            LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,numdofperelement_,1>(1.0, feas_->A(), 1.0, Kad_->A(), res_d.A());
            if (KaT_!=Teuchos::null) LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,nen_,1>(1.,feas_->A(),1.,KaT_->A(),res_T.A());
            LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,soh8p_easfull,1>(0.0,*alpha_eas_inc_,-1.0,*KaaInv_,*feas_);
            LINALG::DENSEFUNCTIONS::update<double,soh8p_easfull,1>(1.0,*alpha_eas_,1.0,*alpha_eas_inc_);
            break;
          case soh8p_easnone: break;
          default: dserror("Don't know what to do with EAS type %d", eastype_); break;
          }
        }
      }
    }
    eas_inc += pow(alpha_eas_inc_->Norm2(),2.);

    // reset EAS matrices
    KaaInv_->Shape(neas_,neas_);
    Kad_->Shape(neas_,numdofperelement_);
    if (KaT_!=Teuchos::null) KaT_->Shape(neas_,nen_);
    feas_->Size(neas_);
    /* end of EAS Update ******************/
  }

  // EAS matrix block
  Epetra_SerialDenseMatrix Kda(numdofperelement_,neas_);
  std::vector<Epetra_SerialDenseVector> dHda(0);
  if (eastype_!=soh8p_easnone && eval_tsi)
    dHda.resize(numgpt_,Epetra_SerialDenseVector(neas_));
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
    LINALG::Matrix<numstr_,1> RCG;// total Cauchy green in strain-like voigt notation
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
    // inverse and determinant
    detF=invdefgrd.Invert(defgrd);

    // calculate modified deformation gradient
    if (fbar_)
    {
      // inverse and determinant
      detF=invdefgrd.Invert(defgrd);

    // check for element distortion
    if (detF<=0. || detF_0<=0.)
    {
      // check, if errors are tolerated or should throw a dserror
      bool error_tol=false;
      if (params.isParameter("tolerate_errors"))
        error_tol=params.get<bool>("tolerate_errors");
      if (error_tol)
      {
        params.set<bool>("eval_error",true);
        stiffmatrix->Clear();
        force->Clear();
        return;
      }
      else
        dserror("element distortion too large");
    }
    // modify deformation gradient
    // ATTENTION: defgrd_mod now contains the F-bar deformation gradient
    // The compatible displacement-based deformation gradient is still
    // available in defgrd.
    f_bar_factor=pow(detF_0/detF,1/3.);
    defgrd_mod.Scale(f_bar_factor);
    }

    // Gauss point temperature
    double gp_temp=-1.e12;
    if (eval_tsi) gp_temp = etemp.Dot(shapefunct);

    // plastic flow increment
    LINALG::Matrix<nsd_,nsd_> deltaLp;

    if (plmat!=NULL)
    {
      // recover plastic variables
      if (HavePlasticSpin())
          RecoverPlasticity<plspin>(res_d,pred,gp,MyPID,gp_temp,params,deltaLp,lp_inc,(stiffmatrix!=NULL && !no_recovery));
      else
          RecoverPlasticity<zerospin>(res_d,pred,gp,MyPID,gp_temp,params,deltaLp,lp_inc,(stiffmatrix!=NULL && !no_recovery));
    }
    // material call *********************************************
    LINALG::Matrix<numstr_,1> pk2;
    LINALG::Matrix<numstr_,numstr_> cmat;
    if (plmat!=NULL)
    {
      plmat->EvaluateElast(&defgrd_mod,&deltaLp,params,&pk2,&cmat,gp,Id());
      if (eval_tsi) plmat->EvaluateThermalStress(&defgrd_mod,gp_temp,params,&pk2,&cmat,gp,Id());
    }
    else
    {
      LINALG::Matrix<numstr_,1> total_glstrain(false);
      total_glstrain(0) = 0.5 * (total_cauchygreen(0,0) - 1.0);
      total_glstrain(1) = 0.5 * (total_cauchygreen(1,1) - 1.0);
      total_glstrain(2) = 0.5 * (total_cauchygreen(2,2) - 1.0);
      total_glstrain(3) = total_cauchygreen(0,1);
      total_glstrain(4) = total_cauchygreen(1,2);
      total_glstrain(5) = total_cauchygreen(2,0);
      params.set<int>("gp",gp);
      SolidMaterial()->Evaluate(&defgrd_mod,&total_glstrain,params,&pk2,&cmat,Id());
    }
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
    if (force_str != NULL && split_res)
    {
      if (fbar_)
        force_str->MultiplyTN(detJ_w/f_bar_factor, bop, pk2, 1.0);
      else
        force_str->MultiplyTN(detJ_w, bop, pk2, 1.0);
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
      if (!no_condensation)
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
          massfactor = shapefunct(jnod) * ifactor;     // intermediate factor
          (*massmatrix)(3*inod+0,3*jnod+0) += massfactor;
          (*massmatrix)(3*inod+1,3*jnod+1) += massfactor;
          (*massmatrix)(3*inod+2,3*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

    if (eval_tsi && (stiffmatrix!=NULL || force!=NULL) && plmat!=NULL)
    {
      // volumetric part of K_dT = Gough-Joule effect**********************************
      // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccc
      // get the thermal material tangent
      LINALG::Matrix<numstr_,1> cTvol(true);
      LINALG::Matrix<numstr_,numstr_> dcTvoldE;
      plmat->EvaluateCTvol(&defgrd_mod,params,&cTvol,&dcTvoldE,gp,Id());
      // end of call material law ccccccccccccccccccccccccccccccccccccccccccccc
      if (fbar_)
        (*dFintdT_)[gp].MultiplyTN(detJ_w/(f_bar_factor),bop,cTvol,0.);
      else
        (*dFintdT_)[gp].MultiplyTN(detJ_w,bop,cTvol,0.);
      if (eastype_!=soh8p_easnone)
      {
        LINALG::Matrix<numstr_,nen_> cTm;
        cTm.MultiplyNT(cTvol,shapefunct);
        switch(eastype_)
        {
        case soh8p_easfull:
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,nen_>
            (1.,KaT_->A(),detJ_w,M.A(),cTm.A());
          break;
        case soh8p_easmild:
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,nen_>
            (1.,KaT_->A(),detJ_w,M.A(),cTm.A());
          break;
        case soh8p_easnone: break;
        default: dserror("Don't know what to do with EAS type %d", eastype_); break;
        }
      }

      // elastic heating ******************************************************
      plmat->HepDiss(gp)=0.;
      plmat->dHepDT(gp) =0.;
      plmat->dHepDissDd(gp).Size(numdofperelement_);
      plmat->dHepDissDd(gp).Scale(0.);
      if (eastype_==soh8p_easnone)
      {
        if (fbar_)
        {
          LINALG::Matrix<3,3> defgrd_rate_0;
          defgrd_rate_0.MultiplyTT(xcurrrate,N_XYZ_0);

          double he_fac;
          double he_fac_deriv;
          plmat->EvaluateGoughJoule(detF_0,Id(),he_fac,he_fac_deriv);

          double fiddfdot =0.;
          for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
              fiddfdot+=invdefgrd_0(j,i)*defgrd_rate_0(i,j);

          double j_dot = 0.;
          for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
              j_dot += detF_0 * invdefgrd_0(j,i)*defgrd_rate_0(i,j);

          double He=he_fac *gp_temp * j_dot;

          plmat->HepDiss(gp)=He;

          // derivative of elastic heating w.r.t. temperature *******************
          plmat->dHepDT(gp) =he_fac * j_dot;

          LINALG::Matrix<numdofperelement_,1> deriv_jdot_d(true);
          LINALG::Matrix<numdofperelement_,1> deriv_j_d(true);

          for (int i=0;i<3;++i)
            for (int n=0;n<numdofperelement_;++n)
              deriv_j_d(n) += detF_0 * invdefgrd_0(i,n%3) * N_XYZ_0(i,n/3);

          double timefac_d=1./(theta*dt);
          LINALG::Matrix<3,3> tmp;
          tmp.Multiply(invdefgrd_0,defgrd_rate_0);
          LINALG::Matrix<3,3> tmp2;
          tmp2.Multiply(tmp,invdefgrd_0);

          deriv_jdot_d.Update(fiddfdot,deriv_j_d,1.);
          for (int i=0;i<3;++i)
            for (int n=0;n<numdofperelement_;++n)
              deriv_jdot_d(n) +=
                                 detF_0 * invdefgrd_0(i,n%3) * N_XYZ_0(i,n/3) * timefac_d
                                -detF_0 * tmp2       (i,n%3) * N_XYZ_0(i,n/3)
                                ;

          LINALG::Matrix<numdofperelement_,1> dHedd(true);
          dHedd.Update(gp_temp*he_fac,deriv_jdot_d,1.);
          dHedd.Update(he_fac_deriv*gp_temp*j_dot,deriv_j_d,1.);

          LINALG::DENSEFUNCTIONS::update<double,numdofperelement_,1>(plmat->dHepDissDd(gp).A(),dHedd.A());
        }

        else
        {
          LINALG::Matrix<3,3> defgrd_rate;
          defgrd_rate.MultiplyTT(xcurrrate,N_XYZ);

          double he_fac;
          double he_fac_deriv;
          plmat->EvaluateGoughJoule(detF,Id(),he_fac,he_fac_deriv);

          double fiddfdot =0.;
          for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
              fiddfdot+=invdefgrd(j,i)*defgrd_rate(i,j);

          double j_dot = 0.;
          for (int i=0;i<3;++i)
            for (int j=0;j<3;++j)
              j_dot += detF * invdefgrd(j,i)*defgrd_rate(i,j);

          double He=he_fac *gp_temp * j_dot;

          plmat->HepDiss(gp)=He;

          // derivative of elastic heating w.r.t. temperature *******************
          plmat->dHepDT(gp) =he_fac * j_dot;

          LINALG::Matrix<numdofperelement_,1> deriv_jdot_d(true);
          LINALG::Matrix<numdofperelement_,1> deriv_j_d(true);

          for (int i=0;i<3;++i)
            for (int n=0;n<numdofperelement_;++n)
              deriv_j_d(n) += detF * invdefgrd(i,n%3) * N_XYZ(i,n/3);

          double timefac_d=1./(theta*dt);
          LINALG::Matrix<3,3> tmp;
          tmp.Multiply(invdefgrd,defgrd_rate);
          LINALG::Matrix<3,3> tmp2;
          tmp2.Multiply(tmp,invdefgrd);

          deriv_jdot_d.Update(fiddfdot,deriv_j_d,1.);
          for (int i=0;i<3;++i)
            for (int n=0;n<numdofperelement_;++n)
              deriv_jdot_d(n) +=
                                 detF * invdefgrd(i,n%3) * N_XYZ(i,n/3) * timefac_d
                                -detF * tmp2     (i,n%3) * N_XYZ(i,n/3)
                                ;

          LINALG::Matrix<numdofperelement_,1> dHedd(true);
          dHedd.Update(gp_temp*he_fac,deriv_jdot_d,1.);
          dHedd.Update(he_fac_deriv*gp_temp*j_dot,deriv_j_d,1.);

          LINALG::DENSEFUNCTIONS::update<double,numdofperelement_,1>(plmat->dHepDissDd(gp).A(),dHedd.A());
        }
      }
      else
      {
        LINALG::Matrix<3,3> defgrd_rate;
        defgrd_rate.MultiplyTT(xcurrrate,N_XYZ);
        double timefac_d=1./(theta*dt);
        // Gough-Joule effect
        plmat->HepDiss(gp)=0.;
        plmat->dHepDT(gp) =0.;
        plmat->dHepDissDd(gp).Size(numdofperelement_);
        // Like this it should be easier to do EAS as well
        LINALG::Matrix<3,3> RCGrate;
        RCGrate.MultiplyTN(defgrd_rate,defgrd);
        RCGrate.MultiplyTN(1.,defgrd,defgrd_rate,1.);
        LINALG::Matrix<6,1> RCGrateVec;
        for (int i=0; i<3; ++i) RCGrateVec(i,0) = RCGrate(i,i);
        RCGrateVec(3,0) = 2.*RCGrate(0,1);
        RCGrateVec(4,0) = 2.*RCGrate(1,2);
        RCGrateVec(5,0) = 2.*RCGrate(0,2);

        // enhance the deformation rate
        if (eastype_!=soh8p_easnone)
        {
          Epetra_SerialDenseVector alpha_dot(neas_);
          switch(eastype_)
          {
          case soh8p_easmild:
            // calculate EAS-rate
            LINALG::DENSEFUNCTIONS::update<double,soh8p_easmild,1>(0.,alpha_dot,1.,*alpha_eas_);
            LINALG::DENSEFUNCTIONS::update<double,soh8p_easmild,1>(1.,alpha_dot,-1.,*alpha_eas_last_timestep_);
            alpha_dot.Scale(timefac_d);
            // enhance the strain rate
            // factor 2 because we deal with RCGrate and not GLrate
            LINALG::DENSEFUNCTIONS::multiply<double,numstr_,soh8p_easmild,1>
              (1.,RCGrateVec.A(),2.,M.A(),alpha_dot.A());
            break;
          case soh8p_easfull:
            // calculate EAS-rate
            LINALG::DENSEFUNCTIONS::update<double,soh8p_easfull,1>(0.,alpha_dot,1.,*alpha_eas_);
            LINALG::DENSEFUNCTIONS::update<double,soh8p_easfull,1>(1.,alpha_dot,-1.,*alpha_eas_last_timestep_);
            alpha_dot.Scale(timefac_d);
            // enhance the strain rate
            // factor 2 because we deal with RCGrate and not GLrate
            LINALG::DENSEFUNCTIONS::multiply<double,numstr_,soh8p_easfull,1>
              (1.,RCGrateVec.A(),2.,M.A(),alpha_dot.A());
            break;
          case soh8p_easnone: break;
          default: dserror("Don't know what to do with EAS type %d", eastype_); break;
          }
        }// enhance the deformation rate

        // heating ************************************************************
        double He=.5*gp_temp*cTvol.Dot(RCGrateVec);

        plmat->HepDiss(gp)=He;

        // derivative of elastic heating w.r.t. temperature *******************
        plmat->dHepDT(gp) =.5*cTvol.Dot(RCGrateVec);

        // derivative of elastic heating w.r.t. displacement ******************
        LINALG::Matrix<numdofperelement_,1> dHedd(true);
        LINALG::Matrix<6,nen_*nsd_> boprate(false);  // (6x24)
        CalculateBop(&boprate, &defgrd_rate, &N_XYZ);

        LINALG::Matrix<6,1> tmp61;
        tmp61.MultiplyTN(dcTvoldE,RCGrateVec);
        dHedd.MultiplyTN(.5*gp_temp,bop,tmp61,1.);

        dHedd.MultiplyTN(timefac_d*gp_temp,bop,cTvol,1.);
        dHedd.MultiplyTN(gp_temp,boprate,cTvol,1.);

        // derivative of elastic heating w.r.t. EAS alphas *******************
        if (eastype_!=soh8p_easnone)
        {
          switch(eastype_)
          {
          case soh8p_easmild:
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,1>
              (0.,dHda[gp].A(),.5*gp_temp,M.A(),tmp61.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,1>
              (1.,dHda[gp].A(),gp_temp*timefac_d,M.A(),cTvol.A());
            break;
          case soh8p_easfull:
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,1>
              (0.,dHda[gp].A(),.5*gp_temp,M.A(),tmp61.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,1>
              (1.,dHda[gp].A(),gp_temp*timefac_d,M.A(),cTvol.A());
            break;
          case soh8p_easnone: break;
          default: dserror("Don't know what to do with EAS type %d", eastype_); break;
          }
        }

        plmat->dHepDissDd(gp).Scale(0.);
        LINALG::DENSEFUNCTIONS::update<double,numdofperelement_,1>(plmat->dHepDissDd(gp).A(),dHedd.A());
      }
    } // if (eval_tsi)

    // plastic modifications
    if ( (stiffmatrix!=NULL || force!=NULL) && !no_condensation && plmat!=NULL)
    {
      if (HavePlasticSpin())
      {
        if (fbar_)
          CondensePlasticity<plspin>(defgrd_mod,deltaLp,bop,&N_XYZ,&RCG,MyPID,detJ_w,
              gp,gp_temp,params,force,stiffmatrix,num_active_gp,lp_res,
              NULL,NULL,NULL,&f_bar_factor,&htensor);
        else if (eastype_!=soh8p_easnone)
          CondensePlasticity<plspin>(defgrd_mod,deltaLp,bop,&N_XYZ,NULL,MyPID,detJ_w,
              gp,gp_temp,params,force,stiffmatrix,num_active_gp,lp_res,&M,&Kda,&dHda);
        else
          CondensePlasticity<plspin>(defgrd_mod,deltaLp,bop,&N_XYZ,NULL,MyPID,detJ_w,
              gp,gp_temp,params,force,stiffmatrix,num_active_gp,lp_res);
      }
      else
      {
        if (fbar_)
          CondensePlasticity<zerospin>(defgrd_mod,deltaLp,bop,&N_XYZ,&RCG,MyPID,detJ_w,
              gp,gp_temp,params,force,stiffmatrix,num_active_gp,lp_res,
              NULL,NULL,NULL,&f_bar_factor,&htensor);
        else if (eastype_!=soh8p_easnone)
          CondensePlasticity<zerospin>(defgrd_mod,deltaLp,bop,&N_XYZ,NULL,MyPID,detJ_w,
              gp,gp_temp,params,force,stiffmatrix,num_active_gp,lp_res,&M,&Kda,&dHda);
        else
          CondensePlasticity<zerospin>(defgrd_mod,deltaLp,bop,&N_XYZ,NULL,MyPID,detJ_w,
              gp,gp_temp,params,force,stiffmatrix,num_active_gp,lp_res);
      }
    }// plastic modifications
  } // gp loop

  // Static condensation EAS --> stiff ********************************
  if (stiffmatrix != NULL && !no_condensation && eastype_!=soh8p_easnone)
  {
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(*KaaInv_);
    solve_for_inverseKaa.Invert();

    Epetra_SerialDenseMatrix kdakaai(numdofperelement_,neas_);
    switch (eastype_)
    {
    case soh8p_easfull:
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easfull,soh8p_easfull>
        (0.,kdakaai.A(),1.,Kda.A(),KaaInv_->A());
      if (stiffmatrix!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easfull,numdofperelement_>
          (1.,stiffmatrix->A(),-1.,kdakaai.A(),Kad_->A());
      if (force!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easfull,1>
          (1.,force->A(),-1.,kdakaai.A(),feas_->A());
      break;
    case soh8p_easmild:
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easmild,soh8p_easmild>
        (0.,kdakaai.A(),1.,Kda.A(),KaaInv_->A());
      if (stiffmatrix!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easmild,numdofperelement_>
          (1.,stiffmatrix->A(),-1.,kdakaai.A(),Kad_->A());
      if (force!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easmild,1>
          (1.,force->A(),-1.,kdakaai.A(),feas_->A());
      break;
    case soh8p_easnone:
      break;
    default: dserror("Don't know what to do with EAS type %d", eastype_); break;
    }

    // TSI with EAS
    if (eval_tsi)
    {
      Epetra_SerialDenseVector dHdaKaai(neas_);
      switch (eastype_)
      {
      case soh8p_easfull:
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easfull,nen_>
          (0.,KdT_eas_->A(),-1.,kdakaai.A(),KaT_->A());
        for (int gp=0; gp<numgpt_; ++gp)
        {
          LINALG::DENSEFUNCTIONS::multiply<double,1,soh8p_easfull,soh8p_easfull>
            (0.,dHdaKaai.A(),1.,dHda.at(gp).A(),KaaInv_->A());
          LINALG::DENSEFUNCTIONS::multiplyTN<double,numdofperelement_,soh8p_easfull,1>
            (1.,plmat->dHepDissDd(gp).A(),-1.,Kad_->A(),dHdaKaai.A());
          LINALG::DENSEFUNCTIONS::multiply<double,1,soh8p_easfull,1>
            (1.,&(plmat->HepDiss(gp)),-1.,dHdaKaai.A(),feas_->A());
          LINALG::DENSEFUNCTIONS::multiplyTN<double,nen_,soh8p_easfull,1>
            (0.,plmat->dHepDTeas()->at(gp).A(),-1.,KaT_->A(),dHdaKaai.A());
        }
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,soh8p_easmild,nen_>
          (0.,KdT_eas_->A(),-1.,kdakaai.A(),KaT_->A());
        for (int gp=0; gp<numgpt_; ++gp)
        {
          LINALG::DENSEFUNCTIONS::multiply<double,1,soh8p_easmild,soh8p_easmild>
            (0.,dHdaKaai.A(),1.,dHda.at(gp).A(),KaaInv_->A());
          LINALG::DENSEFUNCTIONS::multiplyTN<double,numdofperelement_,soh8p_easmild,1>
            (1.,plmat->dHepDissDd(gp).A(),-1.,Kad_->A(),dHdaKaai.A());
          LINALG::DENSEFUNCTIONS::multiply<double,1,soh8p_easmild,1>
            (1.,&(plmat->HepDiss(gp)),-1.,dHdaKaai.A(),feas_->A());
          LINALG::DENSEFUNCTIONS::multiplyTN<double,nen_,soh8p_easmild,1>
            (0.,plmat->dHepDTeas()->at(gp).A(),-1.,KaT_->A(),dHdaKaai.A());
        }
        break;
      case soh8p_easnone:
        break;
      default: dserror("Don't know what to do with EAS type %d", eastype_); break;
      }
    }
    eas_res += pow(feas_uncondensed.Norm2(),2.);
  }

  // rhs norm of eas equations
  if (eastype_!=soh8p_easnone)
    if (split_res)
      // only add for row-map elements
      if (params.get<int>("MyPID")==Owner())
        params.get<double>("cond_rhs_norm") += pow(feas_uncondensed.Norm2(),2.);
  return;
}

/*----------------------------------------------------------------------*
 | coupling term k_dT in monolithic TSI                     seitz 06/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::nln_kdT_tsi(
    LINALG::Matrix<numdofperelement_,nen_>* k_dT,
    Teuchos::ParameterList& params
)
{
  if (k_dT==NULL)
    return;

  // shape functions
  LINALG::Matrix<nen_,1> shapefunct(false);

  for (int gp=0; gp<numgpt_; gp++)
  {
    // shape functions
    DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
    // update linear coupling matrix K_dT
    k_dT->MultiplyNT(1.,(*dFintdT_)[gp],shapefunct,1.);
  }

  // EAS part
  if (eastype_!=soh8p_easnone)
  {
    if (KdT_eas_==Teuchos::null)
      dserror("for TSI with EAS the block KdT_eas_ should be acessible here");
    k_dT->Update(1.,*KdT_eas_,1.);
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
    const LINALG::Matrix<nsd_,nen_>* N_XYZ,
    const LINALG::Matrix<numstr_,1>* RCG,
    const int MyPID,
    const double detJ_w,
    const int gp,
    const double temp,
    Teuchos::ParameterList& params,
    LINALG::Matrix<numdofperelement_,1>* force,
    LINALG::Matrix<numdofperelement_,numdofperelement_>* stiffmatrix,
    int& num_active_gp,
    double& lp_res,
    Epetra_SerialDenseMatrix* M,
    Epetra_SerialDenseMatrix* Kda,
    std::vector<Epetra_SerialDenseVector>* dHda,
    const double* f_bar_factor,
    const LINALG::Matrix<numdofperelement_,1>* htensor
    )
{
  bool eval_tsi=tsi_ && (temp!=-1.e12);

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  // temporary Epetra matrix for matrix-matrix-matrix products
  Epetra_SerialDenseMatrix tmp;

  // Nitsche contact
  LINALG::Matrix<numstr_,1>* cauchy_ptr=NULL;
  LINALG::Matrix<numstr_,spintype+1> d_cauchy_ddp;
  LINALG::Matrix<numstr_,spintype+1>* d_cauchy_ddp_ptr=NULL;
  LINALG::Matrix<numstr_,numstr_> d_cauchy_dC;
  LINALG::Matrix<numstr_,numstr_>* d_cauchy_dC_ptr=NULL;
  LINALG::Matrix<numstr_,9> d_cauchy_dF;
  LINALG::Matrix<numstr_,9>* d_cauchy_dF_ptr=NULL;
  if (is_nitsche_contact_)
  {
    cauchy_ptr=&(cauchy_.at(gp));
    d_cauchy_ddp_ptr=&d_cauchy_ddp;
    d_cauchy_dC_ptr=&d_cauchy_dC;
    d_cauchy_dF_ptr=&d_cauchy_dF;
    if (eval_tsi)
      dserror("no Nitsche thermo contact yet");
  }

  // second material call ****************************************************
  LINALG::Matrix<numstr_,spintype+1> dpk2ddp;
  LINALG::Matrix<spintype+1,1> ncp;
  LINALG::Matrix<spintype+1,spintype+1> dncpddp;
  LINALG::Matrix<spintype+1,numstr_> dncpdc;
  LINALG::Matrix<spintype+1,1> dncpdT;
  LINALG::Matrix<numstr_,1> dHdC;
  LINALG::Matrix<spintype+1,1> dHdLp;
  bool active=false;
  bool elast=false;
  bool as_converged=true;
  if (!eval_tsi)
    plmat->EvaluatePlast(&defgrd,&deltaLp,0,params,&dpk2ddp,
        &ncp,&dncpdc,&dncpddp,&active,&elast,&as_converged,gp,NULL,NULL,NULL,data_->dt_,Id(),
        cauchy_ptr,d_cauchy_ddp_ptr,d_cauchy_dC_ptr,d_cauchy_dF_ptr);
  else
    plmat->EvaluatePlast(&defgrd,&deltaLp,temp,params,&dpk2ddp,
        &ncp,&dncpdc,&dncpddp,&active,&elast,&as_converged,gp,&dncpdT,&dHdC,&dHdLp,data_->dt_,Id());

  if (MyPID==Owner()) lp_res+=pow(ncp.Norm2(),2.);
  if (active && Owner()==MyPID) ++num_active_gp;
  if (as_converged==false) data_->unconverged_active_set_=true;
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

  LINALG::Matrix<numstr_,spintype> d_cauchy_db;
  if (is_nitsche_contact_)
  {
    if (N_XYZ==NULL)
      dserror("shape derivative not provided");

    LINALG::Matrix<9,numdofperelement_> dFdd;
    for (int k=0;k<nen_;++k)
    {
      for (int d=0;d<3;++d)
        dFdd(d,k*nsd_+d)=(*N_XYZ)(d,k);
      dFdd(3,k*nsd_+0)=(*N_XYZ)(1,k);
      dFdd(4,k*nsd_+1)=(*N_XYZ)(2,k);
      dFdd(5,k*nsd_+0)=(*N_XYZ)(2,k);

      dFdd(6,k*nsd_+1)=(*N_XYZ)(0,k);
      dFdd(7,k*nsd_+2)=(*N_XYZ)(1,k);
      dFdd(8,k*nsd_+2)=(*N_XYZ)(0,k);
    }
    cauchy_deriv_.at(gp).Clear();
    if (fbar_)
    {
      if (RCG==NULL) dserror("CondensePlasticity(...) requires RCG in case of FBAR");
      LINALG::Matrix<6,1> tmp61;
      tmp61.Multiply(.5,d_cauchy_dC,(*RCG));
      cauchy_deriv_.at(gp).MultiplyNT((*f_bar_factor)*(*f_bar_factor)*2./3.,tmp61,*htensor,1.);
      cauchy_deriv_.at(gp).Multiply((*f_bar_factor)*(*f_bar_factor),d_cauchy_dC,bop,1.);

      LINALG::Matrix<9,1> f; for(int i=0;i<3;++i) f(i)=defgrd(i,i);
      f(3)=defgrd(0,1); f(4)=defgrd(1,2); f(5)=defgrd(0,2);
      f(6)=defgrd(1,0); f(7)=defgrd(2,1); f(8)=defgrd(2,0);

      tmp61.Multiply(1.,d_cauchy_dF,f,0.);
      cauchy_deriv_.at(gp).MultiplyNT(*f_bar_factor/3.,tmp61,*htensor,1.);
      cauchy_deriv_.at(gp).Multiply(*f_bar_factor,d_cauchy_dF,dFdd,1.);
    }
    else
    {
      cauchy_deriv_.at(gp).Multiply(1.,d_cauchy_dC,bop,1.);
      cauchy_deriv_.at(gp).Multiply(1.,d_cauchy_dF,dFdd,1.);
    }
    d_cauchy_db.Multiply(d_cauchy_ddp,dDpdbeta);
  }

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
      (0.,KbbInv_[gp].A(),1.,voigt_red.A(),dNCPdb.A());

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
      (0.,Kbd_[gp].A(),1.,voigt_red.A(),dNCPdd.A());

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
      (0.,fbeta_[gp].A(),1.,voigt_red.A(),ncp.A());

    // **************************************************************
    // static condensation of inner variables
    // **************************************************************
    //inverse matrix block [k_beta beta]_ij
    Epetra_SerialDenseSolver solve_for_kbbinv;
    solve_for_kbbinv.SetMatrix(KbbInv_[gp]);
    int err = solve_for_kbbinv.Invert();
    if (err != 0)
    {
      // check, if errors are tolerated or should throw a dserror
      bool error_tol=false;
      if (params.isParameter("tolerate_errors"))
        error_tol=params.get<bool>("tolerate_errors");
      if (error_tol)
      {
        params.set<bool>("eval_error",true);
        return;
      }
      else
        dserror("Inversion of Kbb failed");
    }
    // temporary  Kdb.Kbb^-1
    LINALG::Matrix<numdofperelement_,spintype> KdbKbb;
    LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,spintype>
      (0.,KdbKbb.A(),1.,kdbeta.A(),KbbInv_[gp].A());

    // "plastic displacement stiffness"
    // plstiff = [k_d beta] * [k_beta beta]^-1 * [k_beta d]
    if (stiffmatrix!=NULL)
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,numdofperelement_>
        (1.,stiffmatrix->A(),-1.,KdbKbb.A(),Kbd_[gp].A());

    // "plastic internal force"
    // plFint = [K_db.K_bb^-1].f_b
    if (force!=NULL)
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,1>
        (1.,force->A(),-1.,KdbKbb.A(),fbeta_[gp].A());

    // TSI
    if (eval_tsi)
    {
      // thermal derivative
      (*KbT_).at(gp).Size(spintype);
      LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,spintype+1,1>
      (0.,(*KbT_)[gp].A(),1.,voigt_red.A(),dncpdT.A());

      // condense to K_dT
      LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,1>
        (1.,(*dFintdT_)[gp].A(),-1.,KdbKbb.A(),(*KbT_)[gp].A());

      // Plastic heating and dissipation dC
      LINALG::Matrix<numdofperelement_,1> dHepDissDd;
      if (fbar_)
      {
        if (RCG==NULL) dserror("CondensePlasticity(...) requires RCG in case of FBAR");
        LINALG::Matrix<1,1> tmp11;
        tmp11.MultiplyTN(.5,dHdC,(*RCG));
        dHepDissDd.Multiply((*f_bar_factor)*(*f_bar_factor)*2./3.,*htensor,tmp11,0.);
        dHepDissDd.MultiplyTN((*f_bar_factor)*(*f_bar_factor),bop,dHdC,1.);
      }
      else
        dHepDissDd.MultiplyTN(bop,dHdC);

      // store in material
      LINALG::DENSEFUNCTIONS::update<double,numdofperelement_,1>(1.,plmat->dHepDissDd(gp).A(),1.,dHepDissDd.A());

      // Plastic heating and dissipation dbeta
        LINALG::Matrix<spintype,1> dHepDissDbeta;
      if (spintype!=zerospin)
        dserror("no TSI with plastic spin yet");
      else
        LINALG::DENSEFUNCTIONS::multiplyTN<double,zerospin,zerospin+1,1>
          (0.,dHepDissDbeta.A(),1.,dDpdbeta.A(),dHdLp.A());

      // condense the heating terms
      LINALG::Matrix<spintype,1> dHdbKbbi;
      LINALG::DENSEFUNCTIONS::multiplyTN<double,spintype,spintype,1>
        (0.,dHdbKbbi.A(),1.,KbbInv_[gp].A(),dHepDissDbeta.A());
      LINALG::DENSEFUNCTIONS::multiplyTN<double,1,spintype,1>
        (1.,&(plmat->HepDiss(gp)),-1.,dHdbKbbi.A(),fbeta_[gp].A());
      LINALG::DENSEFUNCTIONS::multiplyTN<double,numdofperelement_,spintype,1>
        (1.,plmat->dHepDissDd(gp).A(),-1.,Kbd_[gp].A(),dHdbKbbi.A());
      LINALG::DENSEFUNCTIONS::multiplyTN<double,1,spintype,1>
        (1.,&(plmat->dHepDT(gp)),-1.,(*KbT_)[gp].A(),dHdbKbbi.A());

      // TSI with EAS
      if (eastype_!=soh8p_easnone)
      {
        // error checks
        if (dHda==NULL)
          dserror("dHda is NULL pointer");
        else if ((int)dHda->size()!=numgpt_)
          dserror("dHda has wrong size");

        switch (eastype_)
        {
        case soh8p_easmild:
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,numstr_,1>
            (1.,dHda->at(gp).A(),1.,M->A(),dHdC.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easmild,spintype,1>
            (1.,dHda->at(gp).A(),-1.,Kba_->at(gp).A(),dHdbKbbi.A());

          break;
        case soh8p_easfull:
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,numstr_,1>
            (1.,dHda->at(gp).A(),1.,M->A(),dHdC.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_easfull,spintype,1>
            (1.,dHda->at(gp).A(),-1.,Kba_->at(gp).A(),dHdbKbbi.A());
          break;
        case soh8p_eassosh8:
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_eassosh8,numstr_,1>
            (1.,dHda->at(gp).A(),1.,M->A(),dHdC.A());
          LINALG::DENSEFUNCTIONS::multiplyTN<double,soh8p_eassosh8,spintype,1>
            (1.,dHda->at(gp).A(),-1.,Kba_->at(gp).A(),dHdbKbbi.A());
          break;
        case soh8p_easnone: break;
        default: dserror("Don't know what to do with EAS type %d", eastype_); break;
        }
      }
    } // TSI

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
          (0.,tmp.A(),1.,Kab.A(),KbbInv_[gp].A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,numdofperelement_>
          (1.,Kad_->A(),-1.,tmp.A(),Kbd_[gp].A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,soh8p_easfull>
          (1.,KaaInv_->A(),-1.,tmp.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,1>
          (1.,feas_->A(),-1.,tmp.A(),fbeta_[gp].A());
        if (eval_tsi)
        {
          Epetra_SerialDenseMatrix kbTm(spintype,nen_);
          LINALG::Matrix<nen_,1> shapefunct;
          DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
          LINALG::DENSEFUNCTIONS::multiplyNT<double,spintype,1,nen_>
            (0.,kbTm.A(),1.,KbT_->at(gp).A(),shapefunct.A());
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,nen_>
            (1.,KaT_->A(),-1.,tmp.A(),kbTm.A());
        }
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,soh8p_easmild>
          (1.,Kda->A(),-1.,KdbKbb.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,spintype>
          (0.,tmp.A(),1.,Kab.A(),KbbInv_[gp].A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,numdofperelement_>
          (1.,Kad_->A(),-1.,tmp.A(),Kbd_[gp].A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,soh8p_easmild>
          (1.,KaaInv_->A(),-1.,tmp.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,1>
          (1.,feas_->A(),-1.,tmp.A(),fbeta_[gp].A());
        if (eval_tsi)
        {
          Epetra_SerialDenseMatrix kbTm(spintype,nen_);
          LINALG::Matrix<nen_,1> shapefunct;
          DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
          LINALG::DENSEFUNCTIONS::multiplyNT<double,spintype,1,nen_>
            (0.,kbTm.A(),1.,KbT_->at(gp).A(),shapefunct.A());
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,nen_>
            (1.,KaT_->A(),-1.,tmp.A(),kbTm.A());
        }
        break;
      case soh8p_eassosh8:
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,soh8p_eassosh8>
          (1.,Kda->A(),-1.,KdbKbb.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,spintype>
          (0.,tmp.A(),1.,Kab.A(),KbbInv_[gp].A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,numdofperelement_>
          (1.,Kad_->A(),-1.,tmp.A(),Kbd_[gp].A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,soh8p_eassosh8>
          (1.,KaaInv_->A(),-1.,tmp.A(),Kba_->at(gp).A());
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,1>
          (1.,feas_->A(),-1.,tmp.A(),fbeta_[gp].A());
        if (eval_tsi)
        {
          Epetra_SerialDenseMatrix kbTm(spintype,nen_);
          LINALG::Matrix<nen_,1> shapefunct;
          DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct);
          LINALG::DENSEFUNCTIONS::multiplyNT<double,spintype,1,nen_>
            (0.,kbTm.A(),1.,KbT_->at(gp).A(),shapefunct.A());
          LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,nen_>
            (1.,KaT_->A(),-1.,tmp.A(),kbTm.A());
        }
        break;
      case soh8p_easnone:
        // do nothing
        break;
      default: dserror("Don't know what to do with EAS type %d", eastype_); break;
      }
    }

    if (is_nitsche_contact_)
    {
      LINALG::Matrix<numstr_,spintype> tmp;
      tmp.Multiply(d_cauchy_db,LINALG::Matrix<spintype,spintype>(KbbInv_.at(gp).A(),true));
      cauchy_.at(gp).Multiply(-1.,tmp,LINALG::Matrix<spintype,1>(fbeta_.at(gp).A(),true),1.);
      cauchy_deriv_.at(gp).Multiply(-1.,tmp,LINALG::Matrix<spintype,numdofperelement_>(Kbd_.at(gp).A(),true),1.);
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
    if (dDp_last_iter_[gp].NormInf()>0.)
    {
      if (force!=NULL)
        LINALG::DENSEFUNCTIONS::multiply<double,numdofperelement_,spintype,1>
          (1.,force->A(),-1.,kdbeta.A(),dDp_last_iter_[gp].A());

      switch (eastype_)
      {
      case soh8p_easnone:
        break;
      case soh8p_easmild:
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easmild,spintype,1>
          (1.,feas_->A(),-1.,Kab.A(),dDp_last_iter_[gp].A());
        break;
      case soh8p_easfull:
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_easfull,spintype,1>
          (1.,feas_->A(),-1.,Kab.A(),dDp_last_iter_[gp].A());
        break;
      case soh8p_eassosh8:
        LINALG::DENSEFUNCTIONS::multiply<double,soh8p_eassosh8,spintype,1>
          (1.,feas_->A(),-1.,Kab.A(),dDp_last_iter_[gp].A());
        break;
      default: dserror("Don't know what to do with EAS type %d", eastype_); break;
      }
    }

    KbbInv_[gp].Scale(0.);
    for (int i=0; i<KbbInv_[gp].RowDim(); ++i)
      KbbInv_[gp](i,i)=1./(plmat->cpl());
    fbeta_[gp].Scale(0.);
    fbeta_[gp]=dDp_last_iter_[gp];
    fbeta_[gp].Scale(plmat->cpl());
    Kbd_[gp].Scale(0.);
    if (eastype_!=soh8p_easnone)
      Kba_->at(gp).Shape(spintype,neas_);
    if (KbT_!=Teuchos::null)
      KbT_->at(gp).Scale(0.);
  }

  // rhs norm of eas equations
  if (data_->split_res_)
    // only add for row-map elements
    if (data_->myPID==Owner())
      params.get<double>("cond_rhs_norm") += detJ_w*detJ_w*pow(fbeta_[gp].Norm2(),2.);

  return;
}

/*----------------------------------------------------------------------*
 |  recover plastic degrees of freedom                      seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
template<int spintype>
void DRT::ELEMENTS::So3_Plast<distype>::RecoverPlasticity(
    const LINALG::Matrix<numdofperelement_,1>& res_d,
    const INPAR::STR::PredEnum& pred,
    const int gp,
    const int MyPID,
    const double temp,
    const Teuchos::ParameterList& params,
    LINALG::Matrix<3,3>& deltaLp,
    double& lp_inc,
    bool recover
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
      dDp_last_iter_[gp].Scale(0.);
      break;

    // tangential predictor
    case INPAR::STR::pred_tangdis:
    case INPAR::STR::pred_constvel:
    case INPAR::STR::pred_constacc:
      // do nothing
      break;

    // do usual recovery
    case INPAR::STR::pred_vague:
      // this is a line search step, i.e. the direction of the eas increments
      // has been calculated by a Newton step and now it is only scaled
      if (params.isParameter("alpha_ls"))
      {
        double alpha_ls=params.get<double>("alpha_ls");
        // undo step
        dDp_inc_[gp].Scale(-1.);
        dDp_last_iter_[gp] += dDp_inc_[gp];
        // scale increment
        dDp_inc_[gp].Scale(-1.*alpha_ls);
        // add reduced increment
        dDp_last_iter_[gp] += dDp_inc_[gp];
      }
      // do recovery
      else
      {
        // first part
        LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,1>
        (0.,dDp_inc_[gp].A(),-1.,KbbInv_[gp].A(),fbeta_[gp].A());

        // second part
        LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,numdofperelement_>
        (0.,tmp_m.A(),1.,KbbInv_[gp].A(),Kbd_[gp].A());
        LINALG::DENSEFUNCTIONS::multiply<double,spintype,numdofperelement_,1>
        (1.,dDp_inc_[gp].A(),-1.,tmp_m.A(),res_d.A());

        // thermal part
        if (KbT_!=Teuchos::null)
          if (temp!=-1.e12) // this is the default value provided if no temperatures are available
          {
            LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,1>
            (1.,dDp_inc_[gp].A(),-1.*(temp-(*temp_last_)[gp]),KbbInv_[gp].A(),(*KbT_)[gp].A());
            (*temp_last_)[gp]=temp;
          }

        // EAS part
        if (eastype_!=soh8p_easnone)
        {
          tmp_m.Shape(spintype,neas_);
          switch (eastype_)
          {
          case soh8p_easmild:
            LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,soh8p_easmild>(0.,tmp_m.A(),1.,KbbInv_[gp].A(),Kba_->at(gp).A());
            LINALG::DENSEFUNCTIONS::multiply<double,spintype,soh8p_easmild,1>(1.,dDp_inc_[gp].A(),-1.,tmp_m.A(),alpha_eas_inc_->A());
            break;
          case soh8p_easfull:
            LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,soh8p_easfull>(0.,tmp_m.A(),1.,KbbInv_[gp].A(),Kba_->at(gp).A());
            LINALG::DENSEFUNCTIONS::multiply<double,spintype,soh8p_easfull,1>(1.,dDp_inc_[gp].A(),-1.,tmp_m.A(),alpha_eas_inc_->A());
            break;
          case soh8p_eassosh8:
            LINALG::DENSEFUNCTIONS::multiply<double,spintype,spintype,soh8p_eassosh8>(0.,tmp_m.A(),1.,KbbInv_[gp].A(),Kba_->at(gp).A());
            LINALG::DENSEFUNCTIONS::multiply<double,spintype,soh8p_eassosh8,1>(1.,dDp_inc_[gp].A(),-1.,tmp_m.A(),alpha_eas_inc_->A());
            break;
          case soh8p_easnone:
            break;
          default: dserror("Don't know what to do with EAS type %d", eastype_); break;
          }
        }// EAS part

        LINALG::DENSEFUNCTIONS::update<double,spintype,1>(1.,dDp_last_iter_[gp],1.,dDp_inc_[gp]);
      }
      if (MyPID==Owner())
        lp_inc += pow(dDp_inc_[gp].Norm2(),2.);
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
  deltaLp(0,0) = dDp_last_iter_[gp](0);
  deltaLp(1,1) = dDp_last_iter_[gp](1);
  deltaLp(2,2) = -1.0*(dDp_last_iter_[gp](0)+dDp_last_iter_[gp](1));
  deltaLp(0,1) = dDp_last_iter_[gp](2);
  deltaLp(1,0) = dDp_last_iter_[gp](2);
  deltaLp(1,2) = dDp_last_iter_[gp](3);
  deltaLp(2,1) = dDp_last_iter_[gp](3);
  deltaLp(0,2) = dDp_last_iter_[gp](4);
  deltaLp(2,0) = dDp_last_iter_[gp](4);
  if (spintype==plspin)
  {
    deltaLp(0,1) += dDp_last_iter_[gp](5);
    deltaLp(1,0) -= dDp_last_iter_[gp](5);
    deltaLp(1,2) += dDp_last_iter_[gp](6);
    deltaLp(2,1) -= dDp_last_iter_[gp](6);
    deltaLp(0,2) += dDp_last_iter_[gp](7);
    deltaLp(2,0) -= dDp_last_iter_[gp](7);
  }
  return;
}



/*----------------------------------------------------------------------*
 |  update plastic deformation for nonlinear kinematics     seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::UpdatePlasticDeformation_nln(PlSpinType spintype)
{
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
  {
    // loop over all Gauss points
    for (int gp=0; gp<numgpt_; gp++)
    {
      LINALG::Matrix<3,3> deltaLp;
      deltaLp(0,0) = dDp_last_iter_[gp](0);
      deltaLp(1,1) = dDp_last_iter_[gp](1);
      deltaLp(2,2) = -1.0*(dDp_last_iter_[gp](0)+dDp_last_iter_[gp](1));
      deltaLp(0,1) = dDp_last_iter_[gp](2);
      deltaLp(1,0) = dDp_last_iter_[gp](2);
      deltaLp(1,2) = dDp_last_iter_[gp](3);
      deltaLp(2,1) = dDp_last_iter_[gp](3);
      deltaLp(0,2) = dDp_last_iter_[gp](4);
      deltaLp(2,0) = dDp_last_iter_[gp](4);
      if (spintype==plspin)
      {
        deltaLp(0,1) += dDp_last_iter_[gp](5);
        deltaLp(1,0) -= dDp_last_iter_[gp](5);
        deltaLp(1,2) += dDp_last_iter_[gp](6);
        deltaLp(2,1) -= dDp_last_iter_[gp](6);
        deltaLp(0,2) += dDp_last_iter_[gp](7);
        deltaLp(2,0) -= dDp_last_iter_[gp](7);
      }
      static_cast<MAT::PlasticElastHyper*>(Material().get())->UpdateGP(gp,&deltaLp);

      KbbInv_[gp].Scale(0.);
      Kbd_[gp].Scale(0.);
      fbeta_[gp].Scale(0.);
      if (tsi_)
        (*KbT_)[gp].Scale(0.);
    }
  }
  else
  {
    SolidMaterial()->Update();
  }

  if (eastype_!=soh8p_easnone)
  {
    for (int i=0; i<neas_; i++)
    {
      (*alpha_eas_delta_over_last_timestep_)(i) = (*alpha_eas_)(i)-(*alpha_eas_last_timestep_)(i);
      (*alpha_eas_last_timestep_)(i) = (*alpha_eas_)(i);
    }
    Kad_->Scale(0.);
    KaaInv_->Scale(0.);
    feas_->Scale(0.);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calculate internal energy of the element (private)                  |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::So3_Plast<distype>::CalcIntEnergy(
  std::vector<double>&   disp          , // current displacements
  Teuchos::ParameterList&         params        ) // strain output option
{
  double energy=0.;

  // update element geometry
  LINALG::Matrix<nen_,3> xrefe(false);      // X, material coord. of element
  LINALG::Matrix<nen_,3> xcurr(false);      // x, current  coord. of element

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
    LINALG::Matrix<nsd_,nen_> N_rst_0(false);
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

  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("elastic strain energy in so3plast elements only for plastic material");

  /* evaluation of EAS variables (which are constant for the following):
  ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
  ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
  ** -> T0^{-T}
  */
  std::vector<Epetra_SerialDenseMatrix>* M_GP = NULL;   // EAS matrix M at all GPs
  LINALG::SerialDenseMatrix M;      // EAS matrix M at current GP
  double detJ0;                     // detJ(origin)
  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  LINALG::Matrix<numstr_,numstr_> T0invT;  // trafo matrix

  if (eastype_!=soh8p_easnone)
    EasSetup(&M_GP,detJ0,T0invT,xrefe);

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
    LINALG::Matrix<numstr_,1> RCG;// total Cauchy green in strain-like voigt notation
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

    if (fbar_)
    {
      double detF=defgrd.Determinant();
      double f_bar_factor=std::pow(detF_0/detF,1/3.);
      defgrd_mod.Scale(f_bar_factor);
    }

    double psi = plmat->StrainEnergy(defgrd_mod,gp,Id());

    double detJ_w = detJ*wgt_[gp];
    energy += detJ_w*psi;

  } // gp loop

  return energy;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::GetCauchyAtXi(
        const LINALG::Matrix<3,1>& xi,
        const std::vector<double>& disp,
        const LINALG::Matrix<3,1>& n,
        const LINALG::Matrix<3,1>& t,
        double& sigma_nt,
        Epetra_SerialDenseMatrix* dsntdd,
        Epetra_SerialDenseMatrix* d2sntdd2,
        Epetra_SerialDenseMatrix* d2sntDdDn,
        Epetra_SerialDenseMatrix* d2sntDdDt,
        Epetra_SerialDenseMatrix* d2sntDdDpxi,
        LINALG::Matrix<3,1>* dsntdn,
        LINALG::Matrix<3,1>* dsntdt,
        LINALG::Matrix<3,1>* dsntdpxi)
{
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    GetCauchyAtXiPlast(xi,disp,n,t,sigma_nt,dsntdd,d2sntdd2,d2sntDdDn,d2sntDdDt,
            d2sntDdDpxi,dsntdn,dsntdt,dsntdpxi);
  else
    GetCauchyAtXiElast(xi,disp,n,t,sigma_nt,dsntdd,d2sntdd2,d2sntDdDn,d2sntDdDt,
            d2sntDdDpxi,dsntdn,dsntdt,dsntdpxi);
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::GetCauchyAtXiElast(
        const LINALG::Matrix<3,1>& xi,
        const std::vector<double>& disp,
        const LINALG::Matrix<3,1>& n,
        const LINALG::Matrix<3,1>& t,
        double& sigma_nt,
        Epetra_SerialDenseMatrix* dsntdd,
        Epetra_SerialDenseMatrix* d2sntdd2,
        Epetra_SerialDenseMatrix* d2sntDdDn,
        Epetra_SerialDenseMatrix* d2sntDdDt,
        Epetra_SerialDenseMatrix* d2sntDdDpxi,
        LINALG::Matrix<3,1>* dsntdn,
        LINALG::Matrix<3,1>* dsntdt,
        LINALG::Matrix<3,1>* dsntdpxi)
{
  if (fbar_ || eastype_!=soh8p_easnone)
    dserror("cauchy stress not available for fbar or eas elements");

  sigma_nt=0.;

  LINALG::Matrix<nen_,nsd_> xrefe;  // reference coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();

  for (int i=0; i<nen_; ++i)
  {
    const double* x = nodes[i]->X();
    for (int d =0; d<nsd_; ++d)
    {
      xrefe(i,d) = x[d];
      xcurr(i,d) = xrefe(i,d) + disp[i*nsd_+d];
    }
  }

  LINALG::Matrix<nsd_,nen_> deriv;
  DRT::UTILS::shape_function_deriv1<distype>(xi,deriv);

  LINALG::Matrix<nsd_,nen_> N_XYZ;

  LINALG::Matrix<nsd_,nsd_> invJ;
  invJ.Multiply(deriv,xrefe);
  invJ.Invert();
  N_XYZ.Multiply(invJ,deriv);
  LINALG::Matrix<nsd_,nsd_> defgrd;
  defgrd.MultiplyTT(xcurr,N_XYZ);

  LINALG::Matrix<nsd_,nsd_> b;
  b.MultiplyNT(defgrd,defgrd);

  LINALG::Matrix<nsd_,nen_> F_N_XYZ;
  F_N_XYZ.Multiply(defgrd,N_XYZ);

  LINALG::Matrix<MAT::NUM_STRESS_3D,numdofperelement_> bop; // here: linearization of b=FF^T !!!
  if (dsntdd || d2sntDdDn || d2sntDdDt || d2sntDdDpxi)

  {
    for (int i=0; i<nen_; ++i)
    {
      for (int x=0;x<nsd_;++x)
        bop(x,nsd_*i+x) += F_N_XYZ(x,i); // defgrd(x,g)*N_XYZ(g,i);
      {
        bop(3,nsd_*i+0) +=  F_N_XYZ(1,i); //defgrd(1,g)*N_XYZ(g,i);
        bop(3,nsd_*i+1) +=  F_N_XYZ(0,i); //defgrd(0,g)*N_XYZ(g,i);
        bop(4,nsd_*i+2) +=  F_N_XYZ(1,i); //defgrd(1,g)*N_XYZ(g,i);
        bop(4,nsd_*i+1) +=  F_N_XYZ(2,i); //defgrd(2,g)*N_XYZ(g,i);
        bop(5,nsd_*i+0) +=  F_N_XYZ(2,i); //defgrd(2,g)*N_XYZ(g,i);
        bop(5,nsd_*i+2) +=  F_N_XYZ(0,i); //defgrd(0,g)*N_XYZ(g,i);
      }
    }
    bop.Scale(2.);
  }

  LINALG::Matrix<6,1> dsntdb;
  LINALG::Matrix<6,6> d2sntdb2;
  LINALG::Matrix<6,nsd_> d2sntDbDn;
  LINALG::Matrix<6,nsd_> d2sntDbDt;
  SolidMaterial()->EvaluateCauchy(b,n,t,sigma_nt,&dsntdb,&d2sntdb2,&d2sntDbDn,&d2sntDbDt,dsntdn,dsntdt,0);

  if (dsntdd)
  {
    dsntdd     ->Reshape(numdofperelement_,1);
    LINALG::Matrix<numdofperelement_,1> dsntdd_m(dsntdd->A(),true);
    dsntdd_m.MultiplyTN(bop,dsntdb);
  }

  if (d2sntDdDn)
  {
    d2sntDdDn  ->Reshape(numdofperelement_,nsd_);
    LINALG::Matrix<numdofperelement_,nsd_> d2sntDdDn_m(d2sntDdDn->A(),true);
    d2sntDdDn_m.MultiplyTN(bop,d2sntDbDn);
  }

  if (d2sntDdDt)
  {
    d2sntDdDt  ->Reshape(numdofperelement_,nsd_);
    LINALG::Matrix<numdofperelement_,nsd_> d2sntDdDt_m(d2sntDdDt->A(),true);
    d2sntDdDt_m.MultiplyTN(bop,d2sntDbDt);
  }

  if (d2sntdd2)
  {
    d2sntdd2   ->Reshape(numdofperelement_,numdofperelement_);
    LINALG::Matrix<numdofperelement_,numdofperelement_> d2sntdd2_linalg(d2sntdd2->A(),true);
    d2sntdd2_linalg.Clear();
    LINALG::Matrix<6,numdofperelement_> d2sntdb2Bop;
    d2sntdb2Bop.Multiply(d2sntdb2,bop);
    d2sntdd2_linalg.MultiplyTN(1.,bop,d2sntdb2Bop,1.);

    LINALG::Matrix<nen_,nen_> NN;
    NN.MultiplyTN(N_XYZ,N_XYZ);

    for (int i=0; i<nen_; ++i)
    {
      for (int k=0;k<nen_;++k)
        // for (int g=0;g<3;++g)
      {
        for (int x=0;x<nsd_;++x)
          d2sntdd2_linalg(nsd_*i+x,nsd_*k+x) += dsntdb(x)*2.*NN(i,k); // *N_XYZ(g,i)*N_XYZ(g,k);

        d2sntdd2_linalg(nsd_*i+0,nsd_*k+1)   += dsntdb(3)*2.*NN(i,k); // *N_XYZ(g,i)*N_XYZ(g,k);
        d2sntdd2_linalg(nsd_*i+1,nsd_*k+0)   += dsntdb(3)*2.*NN(i,k); // *N_XYZ(g,i)*N_XYZ(g,k);

        d2sntdd2_linalg(nsd_*i+1,nsd_*k+2)   += dsntdb(4)*2.*NN(i,k); // *N_XYZ(g,i)*N_XYZ(g,k);
        d2sntdd2_linalg(nsd_*i+2,nsd_*k+1)   += dsntdb(4)*2.*NN(i,k); // *N_XYZ(g,i)*N_XYZ(g,k);

        d2sntdd2_linalg(nsd_*i+0,nsd_*k+2)   += dsntdb(5)*2.*NN(i,k); // *N_XYZ(g,i)*N_XYZ(g,k);
        d2sntdd2_linalg(nsd_*i+2,nsd_*k+0)   += dsntdb(5)*2.*NN(i,k); // *N_XYZ(g,i)*N_XYZ(g,k);
      }
    }
  }


  if (d2sntDdDpxi)
  {
    d2sntDdDpxi->Reshape(numdofperelement_,nsd_);

    LINALG::Matrix < DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2,
    nen_ > deriv2;
    DRT::UTILS::shape_function_deriv2<distype>(xi,deriv2);

    LINALG::Matrix<nen_,nsd_> xXF(xcurr);
    xXF.MultiplyNT(-1.,xrefe,defgrd,1.);

    LINALG::Matrix<nsd_,DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2> xXFsec;
    xXFsec.MultiplyTT(1.,xXF,deriv2,0.);

    LINALG::Matrix<nsd_,nsd_> jift;
    jift.MultiplyTT(invJ,defgrd);

    LINALG::Matrix<MAT::NUM_STRESS_3D,nsd_> dbdxi(true);

    int VOIGT3X3SYM_[3][3] = {{0,3,5},{3,1,4},{5,4,2}};

    for (int i=0;i<nsd_;++i)
      for (int j=0;j<nsd_;++j)
      {
        dbdxi(VOIGT3X3SYM_[i][j],0)+=
            xXFsec(i,0)*jift(0,j)
            +xXFsec(i,3)*jift(1,j)
            +xXFsec(i,4)*jift(2,j)
            +xXFsec(j,0)*jift(0,i)
            +xXFsec(j,3)*jift(1,i)
            +xXFsec(j,4)*jift(2,i);

        dbdxi(VOIGT3X3SYM_[i][j],1)+=
            xXFsec(i,3)*jift(0,j)
            +xXFsec(i,1)*jift(1,j)
            +xXFsec(i,5)*jift(2,j)
            +xXFsec(j,3)*jift(0,i)
            +xXFsec(j,1)*jift(1,i)
            +xXFsec(j,5)*jift(2,i);

        dbdxi(VOIGT3X3SYM_[i][j],2)+=
            xXFsec(i,4)*jift(0,j)
            +xXFsec(i,5)*jift(1,j)
            +xXFsec(i,2)*jift(2,j)
            +xXFsec(j,4)*jift(0,i)
            +xXFsec(j,5)*jift(1,i)
            +xXFsec(j,2)*jift(2,i);
      }

    dsntdpxi->MultiplyTN(dbdxi,dsntdb);

    LINALG::Matrix<numdofperelement_,nsd_> d2sntDdDpxi_m(d2sntDdDpxi->A(),true);
    d2sntDdDpxi_m.Clear();

    LINALG::Matrix<6,numdofperelement_> d2sntdb2Bop;
    d2sntdb2Bop.Multiply(d2sntdb2,bop);
    d2sntDdDpxi_m.MultiplyTN(d2sntdb2Bop,dbdxi);

    LINALG::Matrix<nsd_,nen_> invJ_N_XYZ;
    invJ_N_XYZ.MultiplyTN(invJ,N_XYZ);
    LINALG::Matrix <DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2,nsd_ > Xsec;
    Xsec.Multiply(deriv2,xrefe);
    LINALG::Matrix<nen_,6> N_XYZ_Xsec;
    N_XYZ_Xsec.MultiplyTT(N_XYZ,Xsec);
    for (int i=0; i<nen_; ++i)
    {
      {
        int x=0;
        for (x=0;x<nsd_;++x)
        {
          d2sntDdDpxi_m(nsd_*i+x,0)+=dsntdb(x)*2.*(
              invJ_N_XYZ(0,i)*xXFsec    (x,0)
              +invJ_N_XYZ(1,i)*xXFsec    (x,3)
              +invJ_N_XYZ(2,i)*xXFsec    (x,4)
              +jift      (0,x)*deriv2    (0,i)
              +jift      (1,x)*deriv2    (3,i)
              +jift      (2,x)*deriv2    (4,i)
              -jift      (0,x)*N_XYZ_Xsec(i,0)
              -jift      (1,x)*N_XYZ_Xsec(i,3)
              -jift      (2,x)*N_XYZ_Xsec(i,4));

          d2sntDdDpxi_m(nsd_*i+x,1)+=dsntdb(x)*2.*(
              invJ_N_XYZ(0,i)*xXFsec    (x,3)
              +invJ_N_XYZ(1,i)*xXFsec    (x,1)
              +invJ_N_XYZ(2,i)*xXFsec    (x,5)
              +jift      (0,x)*deriv2    (3,i)
              +jift      (1,x)*deriv2    (1,i)
              +jift      (2,x)*deriv2    (5,i)
              -jift      (0,x)*N_XYZ_Xsec(i,3)
              -jift      (1,x)*N_XYZ_Xsec(i,1)
              -jift      (2,x)*N_XYZ_Xsec(i,5));

          d2sntDdDpxi_m(nsd_*i+x,2)+=dsntdb(x)*2.*(
              invJ_N_XYZ(0,i)*xXFsec    (x,4)
              +invJ_N_XYZ(1,i)*xXFsec    (x,5)
              +invJ_N_XYZ(2,i)*xXFsec    (x,2)
              +jift      (0,x)*deriv2    (4,i)
              +jift      (1,x)*deriv2    (5,i)
              +jift      (2,x)*deriv2    (2,i)
              -jift      (0,x)*N_XYZ_Xsec(i,4)
              -jift      (1,x)*N_XYZ_Xsec(i,5)
              -jift      (2,x)*N_XYZ_Xsec(i,2));
        }

        x=1;
        int y=0;
        d2sntDdDpxi_m(nsd_*i+y,0)+=dsntdb(3)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,0)
            +invJ_N_XYZ(1,i)*xXFsec    (x,3)
            +invJ_N_XYZ(2,i)*xXFsec    (x,4)
            +jift      (0,x)*deriv2    (0,i)
            +jift      (1,x)*deriv2    (3,i)
            +jift      (2,x)*deriv2    (4,i)
            -jift      (0,x)*N_XYZ_Xsec(i,0)
            -jift      (1,x)*N_XYZ_Xsec(i,3)
            -jift      (2,x)*N_XYZ_Xsec(i,4));
        d2sntDdDpxi_m(nsd_*i+y,1)+=dsntdb(3)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,3)
            +invJ_N_XYZ(1,i)*xXFsec    (x,1)
            +invJ_N_XYZ(2,i)*xXFsec    (x,5)
            +jift      (0,x)*deriv2    (3,i)
            +jift      (1,x)*deriv2    (1,i)
            +jift      (2,x)*deriv2    (5,i)
            -jift      (0,x)*N_XYZ_Xsec(i,3)
            -jift      (1,x)*N_XYZ_Xsec(i,1)
            -jift      (2,x)*N_XYZ_Xsec(i,5));
        d2sntDdDpxi_m(nsd_*i+y,2)+=dsntdb(3)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,4)
            +invJ_N_XYZ(1,i)*xXFsec    (x,5)
            +invJ_N_XYZ(2,i)*xXFsec    (x,2)
            +jift      (0,x)*deriv2    (4,i)
            +jift      (1,x)*deriv2    (5,i)
            +jift      (2,x)*deriv2    (2,i)
            -jift      (0,x)*N_XYZ_Xsec(i,4)
            -jift      (1,x)*N_XYZ_Xsec(i,5)
            -jift      (2,x)*N_XYZ_Xsec(i,2));
        x=0;y=1;
        d2sntDdDpxi_m(nsd_*i+y,0)+=dsntdb(3)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,0)
            +invJ_N_XYZ(1,i)*xXFsec    (x,3)
            +invJ_N_XYZ(2,i)*xXFsec    (x,4)
            +jift      (0,x)*deriv2    (0,i)
            +jift      (1,x)*deriv2    (3,i)
            +jift      (2,x)*deriv2    (4,i)
            -jift      (0,x)*N_XYZ_Xsec(i,0)
            -jift      (1,x)*N_XYZ_Xsec(i,3)
            -jift      (2,x)*N_XYZ_Xsec(i,4));
        d2sntDdDpxi_m(nsd_*i+y,1)+=dsntdb(3)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,3)
            +invJ_N_XYZ(1,i)*xXFsec    (x,1)
            +invJ_N_XYZ(2,i)*xXFsec    (x,5)
            +jift      (0,x)*deriv2    (3,i)
            +jift      (1,x)*deriv2    (1,i)
            +jift      (2,x)*deriv2    (5,i)
            -jift      (0,x)*N_XYZ_Xsec(i,3)
            -jift      (1,x)*N_XYZ_Xsec(i,1)
            -jift      (2,x)*N_XYZ_Xsec(i,5));
        d2sntDdDpxi_m(nsd_*i+y,2)+=dsntdb(3)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,4)
            +invJ_N_XYZ(1,i)*xXFsec    (x,5)
            +invJ_N_XYZ(2,i)*xXFsec    (x,2)
            +jift      (0,x)*deriv2    (4,i)
            +jift      (1,x)*deriv2    (5,i)
            +jift      (2,x)*deriv2    (2,i)
            -jift      (0,x)*N_XYZ_Xsec(i,4)
            -jift      (1,x)*N_XYZ_Xsec(i,5)
            -jift      (2,x)*N_XYZ_Xsec(i,2));

        x=1;y=2;
        d2sntDdDpxi_m(nsd_*i+y,0)+=dsntdb(4)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,0)
            +invJ_N_XYZ(1,i)*xXFsec    (x,3)
            +invJ_N_XYZ(2,i)*xXFsec    (x,4)
            +jift      (0,x)*deriv2    (0,i)
            +jift      (1,x)*deriv2    (3,i)
            +jift      (2,x)*deriv2    (4,i)
            -jift      (0,x)*N_XYZ_Xsec(i,0)
            -jift      (1,x)*N_XYZ_Xsec(i,3)
            -jift      (2,x)*N_XYZ_Xsec(i,4));
        d2sntDdDpxi_m(nsd_*i+y,1)+=dsntdb(4)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,3)
            +invJ_N_XYZ(1,i)*xXFsec    (x,1)
            +invJ_N_XYZ(2,i)*xXFsec    (x,5)
            +jift      (0,x)*deriv2    (3,i)
            +jift      (1,x)*deriv2    (1,i)
            +jift      (2,x)*deriv2    (5,i)
            -jift      (0,x)*N_XYZ_Xsec(i,3)
            -jift      (1,x)*N_XYZ_Xsec(i,1)
            -jift      (2,x)*N_XYZ_Xsec(i,5));
        d2sntDdDpxi_m(nsd_*i+y,2)+=dsntdb(4)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,4)
            +invJ_N_XYZ(1,i)*xXFsec    (x,5)
            +invJ_N_XYZ(2,i)*xXFsec    (x,2)
            +jift      (0,x)*deriv2    (4,i)
            +jift      (1,x)*deriv2    (5,i)
            +jift      (2,x)*deriv2    (2,i)
            -jift      (0,x)*N_XYZ_Xsec(i,4)
            -jift      (1,x)*N_XYZ_Xsec(i,5)
            -jift      (2,x)*N_XYZ_Xsec(i,2));
        x=2;y=1;
        d2sntDdDpxi_m(nsd_*i+y,0)+=dsntdb(4)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,0)
            +invJ_N_XYZ(1,i)*xXFsec    (x,3)
            +invJ_N_XYZ(2,i)*xXFsec    (x,4)
            +jift      (0,x)*deriv2    (0,i)
            +jift      (1,x)*deriv2    (3,i)
            +jift      (2,x)*deriv2    (4,i)
            -jift      (0,x)*N_XYZ_Xsec(i,0)
            -jift      (1,x)*N_XYZ_Xsec(i,3)
            -jift      (2,x)*N_XYZ_Xsec(i,4));
        d2sntDdDpxi_m(nsd_*i+y,1)+=dsntdb(4)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,3)
            +invJ_N_XYZ(1,i)*xXFsec    (x,1)
            +invJ_N_XYZ(2,i)*xXFsec    (x,5)
            +jift      (0,x)*deriv2    (3,i)
            +jift      (1,x)*deriv2    (1,i)
            +jift      (2,x)*deriv2    (5,i)
            -jift      (0,x)*N_XYZ_Xsec(i,3)
            -jift      (1,x)*N_XYZ_Xsec(i,1)
            -jift      (2,x)*N_XYZ_Xsec(i,5));
        d2sntDdDpxi_m(nsd_*i+y,2)+=dsntdb(4)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,4)
            +invJ_N_XYZ(1,i)*xXFsec    (x,5)
            +invJ_N_XYZ(2,i)*xXFsec    (x,2)
            +jift      (0,x)*deriv2    (4,i)
            +jift      (1,x)*deriv2    (5,i)
            +jift      (2,x)*deriv2    (2,i)
            -jift      (0,x)*N_XYZ_Xsec(i,4)
            -jift      (1,x)*N_XYZ_Xsec(i,5)
            -jift      (2,x)*N_XYZ_Xsec(i,2));

        x=0;y=2;
        d2sntDdDpxi_m(nsd_*i+y,0)+=dsntdb(5)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,0)
            +invJ_N_XYZ(1,i)*xXFsec    (x,3)
            +invJ_N_XYZ(2,i)*xXFsec    (x,4)
            +jift      (0,x)*deriv2    (0,i)
            +jift      (1,x)*deriv2    (3,i)
            +jift      (2,x)*deriv2    (4,i)
            -jift      (0,x)*N_XYZ_Xsec(i,0)
            -jift      (1,x)*N_XYZ_Xsec(i,3)
            -jift      (2,x)*N_XYZ_Xsec(i,4));
        d2sntDdDpxi_m(nsd_*i+y,1)+=dsntdb(5)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,3)
            +invJ_N_XYZ(1,i)*xXFsec    (x,1)
            +invJ_N_XYZ(2,i)*xXFsec    (x,5)
            +jift      (0,x)*deriv2    (3,i)
            +jift      (1,x)*deriv2    (1,i)
            +jift      (2,x)*deriv2    (5,i)
            -jift      (0,x)*N_XYZ_Xsec(i,3)
            -jift      (1,x)*N_XYZ_Xsec(i,1)
            -jift      (2,x)*N_XYZ_Xsec(i,5));
        d2sntDdDpxi_m(nsd_*i+y,2)+=dsntdb(5)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,4)
            +invJ_N_XYZ(1,i)*xXFsec    (x,5)
            +invJ_N_XYZ(2,i)*xXFsec    (x,2)
            +jift      (0,x)*deriv2    (4,i)
            +jift      (1,x)*deriv2    (5,i)
            +jift      (2,x)*deriv2    (2,i)
            -jift      (0,x)*N_XYZ_Xsec(i,4)
            -jift      (1,x)*N_XYZ_Xsec(i,5)
            -jift      (2,x)*N_XYZ_Xsec(i,2));
        x=2;y=0;
        d2sntDdDpxi_m(nsd_*i+y,0)+=dsntdb(5)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,0)
            +invJ_N_XYZ(1,i)*xXFsec    (x,3)
            +invJ_N_XYZ(2,i)*xXFsec    (x,4)
            +jift      (0,x)*deriv2    (0,i)
            +jift      (1,x)*deriv2    (3,i)
            +jift      (2,x)*deriv2    (4,i)
            -jift      (0,x)*N_XYZ_Xsec(i,0)
            -jift      (1,x)*N_XYZ_Xsec(i,3)
            -jift      (2,x)*N_XYZ_Xsec(i,4));
        d2sntDdDpxi_m(nsd_*i+y,1)+=dsntdb(5)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,3)
            +invJ_N_XYZ(1,i)*xXFsec    (x,1)
            +invJ_N_XYZ(2,i)*xXFsec    (x,5)
            +jift      (0,x)*deriv2    (3,i)
            +jift      (1,x)*deriv2    (1,i)
            +jift      (2,x)*deriv2    (5,i)
            -jift      (0,x)*N_XYZ_Xsec(i,3)
            -jift      (1,x)*N_XYZ_Xsec(i,1)
            -jift      (2,x)*N_XYZ_Xsec(i,5));
        d2sntDdDpxi_m(nsd_*i+y,2)+=dsntdb(5)*2.*(
            invJ_N_XYZ(0,i)*xXFsec    (x,4)
            +invJ_N_XYZ(1,i)*xXFsec    (x,5)
            +invJ_N_XYZ(2,i)*xXFsec    (x,2)
            +jift      (0,x)*deriv2    (4,i)
            +jift      (1,x)*deriv2    (5,i)
            +jift      (2,x)*deriv2    (2,i)
            -jift      (0,x)*N_XYZ_Xsec(i,4)
            -jift      (1,x)*N_XYZ_Xsec(i,5)
            -jift      (2,x)*N_XYZ_Xsec(i,2));
      }
    }
  }
  return;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::GetCauchyAtXiPlast(
        const LINALG::Matrix<3,1>& xi,
        const std::vector<double>& disp,
        const LINALG::Matrix<3,1>& n,
        const LINALG::Matrix<3,1>& t,
        double& sigma_nt,
        Epetra_SerialDenseMatrix* dsntdd,
        Epetra_SerialDenseMatrix* d2sntdd2,
        Epetra_SerialDenseMatrix* d2sntDdDn,
        Epetra_SerialDenseMatrix* d2sntDdDt,
        Epetra_SerialDenseMatrix* d2sntDdDpxi,
        LINALG::Matrix<3,1>* dsntdn,
        LINALG::Matrix<3,1>* dsntdt,
        LINALG::Matrix<3,1>* dsntdpxi)
{
  if (distype!=DRT::Element::hex8 || numgpt_!=8)
    dserror("only for hex8 with 8 gp");
  if (Material()->MaterialType()!=INPAR::MAT::m_plelasthyper)
    dserror("only PlasticElastHyper materials here");
  if ((int)cauchy_.size()!=numgpt_ || (int)cauchy_deriv_.size()!=numgpt_)
    dserror("have you evaluated the cauchy stress???");

  sigma_nt=0.;
  if (dsntdpxi) dsntdpxi->Clear();

  LINALG::Matrix<3,3>nttn;
  nttn.MultiplyNT(.5,n,t,1.);
  nttn.MultiplyNT(.5,t,n,1.);
  LINALG::Matrix<6,1> nttn_v;
  for (int i=0;i<3;++i) nttn_v(i)=nttn(i,i);
  nttn_v(3)=nttn(0,1)+nttn(1,0);
  nttn_v(4)=nttn(2,1)+nttn(1,2);
  nttn_v(5)=nttn(0,2)+nttn(2,0);

  LINALG::Matrix<3,1> xi_expol(xi);
  xi_expol.Scale(sqrt(3.));

  LINALG::Matrix<nen_,1> shapefunct;
  LINALG::Matrix<nsd_,nen_> deriv;
  DRT::UTILS::shape_function<distype>(xi_expol,shapefunct);
  DRT::UTILS::shape_function_deriv1<distype>(xi_expol,deriv);

  LINALG::Matrix<numstr_,1> cauchy_expol;
  LINALG::Matrix<numstr_,numdofperelement_> cauchy_deriv_expol;

  LINALG::Matrix<6,1> tmp61;
  for (int gp=0;gp<numgpt_;++gp)
  {
    cauchy_expol.Update(shapefunct(gp),cauchy_.at(gp),1.);
    cauchy_deriv_expol.Update(shapefunct(gp),cauchy_deriv_.at(gp),1.);
    if (dsntdpxi)
      for (int d=0;d<nsd_;++d)
        (*dsntdpxi)(d)+=cauchy_.at(gp).Dot(nttn_v)*deriv(d,gp)*sqrt(3.);
  }

  sigma_nt=cauchy_expol.Dot(nttn_v);

  if (dsntdd)
  {
    dsntdd->Reshape(numdofperelement_,1);
    LINALG::Matrix<numdofperelement_,1>(dsntdd->A(),true).
        MultiplyTN(cauchy_deriv_expol,nttn_v);
  }
  if (d2sntdd2)
  {
    d2sntdd2   ->Reshape(numdofperelement_,numdofperelement_);
    d2sntdd2   ->Scale(0.);
  }

  LINALG::Matrix<numstr_,nsd_> d_nttn_v_dn,d_nttn_v_dt;
  for (int i=0;i<nsd_;++i)
    for (int j=0;j<nsd_;++j)
      for (int a=0;a<nsd_;++a)
      {
        d_nttn_v_dn(VOIGT3X3SYM_[i][j],a)+=.5*((i==a)*t(j)+(j==a)*t(i));
        d_nttn_v_dt(VOIGT3X3SYM_[i][j],a)+=.5*((i==a)*n(j)+(j==a)*n(i));
      }
  if (dsntdn)
    dsntdn->MultiplyTN(d_nttn_v_dn,cauchy_expol);
  if (dsntdt)
    dsntdt->MultiplyTN(d_nttn_v_dt,cauchy_expol);

  if (d2sntDdDn)
  {
    d2sntDdDn  ->Reshape(numdofperelement_,nsd_);
    LINALG::Matrix<numdofperelement_,nsd_>(d2sntDdDn->A(),true).
        MultiplyTN(cauchy_deriv_expol,d_nttn_v_dn);
  }

  if (d2sntDdDt)
  {
    d2sntDdDt  ->Reshape(numdofperelement_,nsd_);
    LINALG::Matrix<numdofperelement_,nsd_>(d2sntDdDt->A(),true).
        MultiplyTN(cauchy_deriv_expol,d_nttn_v_dt);
  }
  if (d2sntDdDpxi)
  {
    d2sntDdDpxi   ->Reshape(numdofperelement_,nsd_);
    d2sntDdDpxi   ->Scale(0.);
  }
}
