/*----------------------------------------------------------------------*/
/*!
\file acinus_impl.cpp

\brief Internal implementation of RedAcinus element. Methods implemented here
       are called by acinus_evaluate.cpp by DRT::ELEMENTS::RedAcinus::Evaluate()
       with the corresponding action.

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/




#include "acinus_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/maxwell_0d_acinus.H"
#include "../drt_mat/air_0d_O2_saturation.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinusImplInterface* DRT::ELEMENTS::RedAcinusImplInterface::Impl(DRT::ELEMENTS::RedAcinus* red_acinus)
{
  switch (red_acinus->Shape())
  {
  case DRT::Element::line2:
  {
    static AcinusImpl<DRT::Element::line2>* acinus;
    if (acinus==NULL)
    {
      acinus = new AcinusImpl<DRT::Element::line2>;
    }
    return acinus;
  }
  default:
    dserror("shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->NumNode());
    break;
  }
  return NULL;
}


/*----------------------------------------------------------------------*
  | constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcinusImpl<distype>::AcinusImpl()
{

}


/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcinusImpl<distype>::Evaluate(
  RedAcinus*                 ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  Teuchos::RCP<MAT::Material> mat)
{

  const int elemVecdim = elevec1_epetra.Length () ;
  std::vector<int>::iterator it_vcr;

  //Get control parameters for time integration
  // get time-step size
  const double dt = params.get<double>("time step size");
  // get time
  const double time = params.get<double>("total time");

  //Get control parameters for stabilization and higher-order elements (currently unused)
  // flag for higher order elements
  // bool higher_order_ele = ele->isHigherOrderElement(distype);

  //Get all general state vectors: flow, pressure,
  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");
  Teuchos::RCP<const Epetra_Vector> pn   = discretization.GetState("pn");
  Teuchos::RCP<const Epetra_Vector> pnm  = discretization.GetState("pnm");

  Teuchos::RCP<const Epetra_Vector> ial  = discretization.GetState("intr_ac_link");

  Teuchos::RCP<Epetra_Vector> acinar_vnp  = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  Teuchos::RCP<Epetra_Vector> acinar_vn   = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vn");


  Teuchos::RCP<Epetra_Vector> qin_nm  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_nm");
  Teuchos::RCP<Epetra_Vector> qin_n   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_n");
  Teuchos::RCP<Epetra_Vector> qin_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  Teuchos::RCP<Epetra_Vector> qout_np = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");
  Teuchos::RCP<Epetra_Vector> qout_n  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_n");
  Teuchos::RCP<Epetra_Vector> qout_nm = params.get<Teuchos::RCP<Epetra_Vector> >("qout_nm");


  if (pnp==Teuchos::null || pn==Teuchos::null || pnm==Teuchos::null )
    dserror("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  //Extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  //Extract local values from the global vectors
  std::vector<double> mypn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,mypn,lm);

  //Extract local values from the global vectors
  std::vector<double> mypnm(lm.size());
  DRT::UTILS::ExtractMyValues(*pnm,mypnm,lm);

  //Extract local values from the global vectors
  std::vector<double> myial(lm.size());
  DRT::UTILS::ExtractMyValues(*ial,myial,lm);

  //Create objects for element arrays
  Epetra_SerialDenseVector epnp(elemVecdim);
  Epetra_SerialDenseVector epn (elemVecdim);
  Epetra_SerialDenseVector epnm(elemVecdim);
  for (int i=0;i<elemVecdim;++i)
  {
    //Split area and volumetric flow rate, insert into element arrays
    epnp(i)   = mypnp[i];
    epn(i)    = mypn[i];
    epnm(i)   = mypnm[i];
  }

  double e_acin_e_vnp;
  double e_acin_e_vn;

  for (int i=0;i<elemVecdim;++i)
  {
    //Split area and volumetric flow rate, insert into element arrays
    e_acin_e_vnp = (*acinar_vnp)[ele->LID()];
    e_acin_e_vn  = (*acinar_vn )[ele->LID()];
  }

  //Get the volumetric flow rate from the previous time step
  Teuchos::ParameterList elem_params;
  elem_params.set<double>("qout_np",(*qout_np)[ele->LID()]);
  elem_params.set<double>("qout_n" ,(*qout_n )[ele->LID()]);
  elem_params.set<double>("qout_nm",(*qout_nm)[ele->LID()]);
  elem_params.set<double>("qin_np" ,(*qin_np )[ele->LID()]);
  elem_params.set<double>("qin_n"  ,(*qin_n  )[ele->LID()]);
  elem_params.set<double>("qin_nm" ,(*qin_nm )[ele->LID()]);

  elem_params.set<double>("acin_vnp" ,e_acin_e_vnp);
  elem_params.set<double>("acin_vn"  ,e_acin_e_vn );

  elem_params.set<double>("lungVolume_np",params.get<double>("lungVolume_np"));
  elem_params.set<double>("lungVolume_n",params.get<double>("lungVolume_n"));
  elem_params.set<double>("lungVolume_nm",params.get<double>("lungVolume_nm"));

  //Call routine for calculating element matrix and right hand side
  Sysmat(ele,
         epnp,
         epn,
         epnm,
         elemat1_epetra,
         elevec1_epetra,
         mat,
         elem_params,
         time,
         dt);

  double Ao = 0.0;
  ele->getParams("Area",Ao);

  //Put zeros on second line of matrix and rhs in case of interacinar linker
  if (myial[1] > 0.0)
  {
    elemat1_epetra(1,0)=0.0;
    elemat1_epetra(1,1)=0.0;
    elevec1_epetra(1)  =0.0;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::Initial(
  RedAcinus*                             ele,
  Teuchos::ParameterList&                params,
  DRT::Discretization&                   discretization,
  std::vector<int>&                      lm,
  Teuchos::RCP<const MAT::Material>      material)
{

  const int   myrank  = discretization.Comm().MyPID();

  Teuchos::RCP<Epetra_Vector> p0np    = params.get<Teuchos::RCP<Epetra_Vector> >("p0np");
  Teuchos::RCP<Epetra_Vector> p0n     = params.get<Teuchos::RCP<Epetra_Vector> >("p0n");
  Teuchos::RCP<Epetra_Vector> p0nm    = params.get<Teuchos::RCP<Epetra_Vector> >("p0nm");

  Teuchos::RCP<Epetra_Vector> generations   = params.get<Teuchos::RCP<Epetra_Vector> >("generations");
  Teuchos::RCP<Epetra_Vector> a_bc          = params.get<Teuchos::RCP<Epetra_Vector> >("acini_bc");

  //  Teuchos::RCP<Epetra_Vector> a_volume      = params.get<Teuchos::RCP<Epetra_Vector> >("acini_volume");
  Teuchos::RCP<Epetra_Vector> a_e_volume    = params.get<Teuchos::RCP<Epetra_Vector> >("acini_e_volume");

  std::vector<int> lmstride;
  Teuchos::RCP<std::vector<int> > lmowner = Teuchos::rcp(new std::vector<int>);
  ele->LocationVector(discretization,lm,*lmowner,lmstride);

  //Initialize the pressure vectors
  if(myrank == (*lmowner)[0])
  {
    int    gid = lm[0];
    double val = 0.0;
    p0np->ReplaceGlobalValues(1,&val,&gid);
    p0n ->ReplaceGlobalValues(1,&val,&gid);
    p0nm->ReplaceGlobalValues(1,&val,&gid);
  }

  //Find the volume of an acinus element
  {
    int    gid2 = ele->Id();
    double acin_vol = 0.0;
    ele->getParams("AcinusVolume",acin_vol);
    a_e_volume->ReplaceGlobalValues(1,&acin_vol,&gid2);
  }

  //Get the generation numbers
  for (int i = 0; i<2; i++)
  {
    if(ele->Nodes()[i]->GetCondition("RedAirwayEvalLungVolCond"))
    {
      // find the acinus condition
      int    gid = ele->Id();
      double val = 1.0;
      a_bc->ReplaceGlobalValues(1,&val,&gid);
    }
  }
  {
    int    gid = ele->Id();
    int    generation = -1;
    double val = double(generation);
    generations->ReplaceGlobalValues(1,&val,&gid);
  }

  bool solveScatra = params.get<bool>("solveScatra");
  Teuchos::RCP<Epetra_Vector>  junVolMix_Corrector;
  Teuchos::RCP<Epetra_Vector>  scatranp;
  Teuchos::RCP<Epetra_Vector>  e1scatranp;
  Teuchos::RCP<Epetra_Vector>  e2scatranp;

  if (solveScatra)
  {
    junVolMix_Corrector  = params.get<Teuchos::RCP<Epetra_Vector> >("junVolMix_Corrector");
    scatranp             = params.get<Teuchos::RCP<Epetra_Vector> >("scatranp");
    e1scatranp           = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
    e2scatranp           = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");
  }

  if (solveScatra)
  {
    double A=0.0;
    double V=0.0;
    ele->getParams("AcinusVolume",V);
    ele->getParams("Area",A);
    int    gid = lm[1];
    junVolMix_Corrector->ReplaceGlobalValues(1,&A,&gid);

    for (int sci=0;sci<iel;sci++)
    {
      int sgid = lm[sci];
      int    esgid = ele->Id();
      // -------------------------------------------------------------
      //
      // -------------------------------------------------------------
      if (ele->Nodes()[sci]->GetCondition("RedAirwayScatraAirCond"))
      {
        double intSat = ele->Nodes()[sci]->GetCondition("RedAirwayScatraAirCond")->GetDouble("INITIAL_CONCENTRATION");
        int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
        // check if O2 properties material exists
        if (id==-1)
        {
          dserror("A material defining O2 properties in air could not be found");
          exit(1);
        }
        const MAT::PAR::Parameter* smat = DRT::Problem::Instance()->Materials()->ParameterById(id);
        const MAT::PAR::Air_0d_O2_saturation* actmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(smat);

        // get atmospheric pressure
        double patm = actmat->atmospheric_p_;
        // get number of O2 moles per unit volume of O2
        double nO2perVO2 = actmat->nO2_per_VO2_;

        // calculate the PO2 at nodes
        double pO2 = intSat*patm;

        // calculate VO2
        double vO2 = V*(pO2/patm);
        // evaluate initial concentration
        double intConc = nO2perVO2*vO2/V;

        scatranp->ReplaceGlobalValues(1,&intConc,&sgid);
        e1scatranp->ReplaceGlobalValues(1,&intConc,&esgid);
        e2scatranp->ReplaceGlobalValues(1,&intConc,&esgid);
      }
      else
      {
        dserror("0D Acinus scatra must be predefined as \"air\" only");
        exit(1);
      }
    }
  }
}//AcinusImpl::Initial


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::Sysmat(
  RedAcinus*                               ele,
  Epetra_SerialDenseVector&                epnp,
  Epetra_SerialDenseVector&                epn,
  Epetra_SerialDenseVector&                epnm,
  Epetra_SerialDenseMatrix&                sysmat,
  Epetra_SerialDenseVector&                rhs,
  Teuchos::RCP<const MAT::Material>        material,
  Teuchos::ParameterList &                          params,
  double                                   time,
  double                                   dt)
{

  // Decide which acinus material should be used
  if((material->MaterialType() == INPAR::MAT::m_0d_maxwell_acinus_neohookean) ||
      (material->MaterialType() == INPAR::MAT::m_0d_maxwell_acinus_exponential) ||
      (material->MaterialType() == INPAR::MAT::m_0d_maxwell_acinus_doubleexponential) ||
      (material->MaterialType() == INPAR::MAT::m_0d_maxwell_acinus_ogden))
  {
    double VolAcinus;
    ele->getParams("AcinusVolume",VolAcinus);
    double volAlvDuct;
    ele->getParams("AlveolarDuctVolume",volAlvDuct);
    const double NumOfAcini = double(floor(VolAcinus/volAlvDuct));

    const Teuchos::RCP<MAT::Maxwell_0d_acinus> acinus_mat = Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus>(ele->Material());

    //Evaluate material law for acinus
    acinus_mat->Evaluate(epnp,
                         epn,
                         epnm,
                         sysmat,
                         rhs,
                         params,
                         NumOfAcini,
                         volAlvDuct,
                         time,
                         dt);
  }
  else
  {
    dserror("Material law is not a valid reduced dimensional lung acinus material.");
    exit(1);
  }
}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::EvaluateTerminalBC(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Epetra_SerialDenseVector&    rhs,
  Teuchos::RCP<MAT::Material>   material)
{
  const int   myrank  = discretization.Comm().MyPID();

  //Get total time
  const double time = params.get<double>("total time");

  //The number of nodes
  const int numnode = lm.size();
  std::vector<int>::iterator it_vcr;

  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (pnp==Teuchos::null)
    dserror("Cannot get state vectors 'pnp'");

  //Extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  //Create objects for element arrays
  Epetra_SerialDenseVector epnp(numnode);

  //Get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    //Split area and volumetric flow rate, insert into element arrays
    epnp(i)    = mypnp[i];
  }

  /**
   * Resolve the BCs
   **/
  for(int i = 0; i<ele->NumNode(); i++)
  {
    if (ele->Nodes()[i]->Owner()== myrank)
    {
      if(ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond") || ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond") || ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond"))
      {
        std::string Bc;
        double BCin = 0.0;
        if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
        {
          DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond");
          //Get the type of prescribed bc
          Bc = *(condition->Get<std::string>("boundarycond"));

          const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");
          const  std::vector<int>*    curve  = condition->Get<std::vector<int>    >("curve");
          const std::vector<int>*     functions = condition->Get<std::vector<int> >("funct");

          // Read in the value of the applied BC
          // Get factor of first CURVE
          double curvefac = 1.0;
          if((*curve)[0]>=0)
          {
            curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
            BCin = (*vals)[0]*curvefac;
          }
          else
          {
            dserror("no boundary condition defined!");
            exit(1);
          }

          // Get factor of FUNCT
          int functnum = -1;
          if (functions) functnum = (*functions)[0];
          else functnum = -1;

          double functionfac = 0.0;
          if(functnum>0)
          {
            functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,(ele->Nodes()[i])->X(),time,NULL);
          }

          // Get factor of second CURVE
          int curve2num = -1;
          double curve2fac = 1.0;
          if (curve) curve2num = (*curve)[1];
          if (curve2num>=0 )
            curve2fac = DRT::Problem::Instance()->Curve(curve2num).f(time);

          // Add first_CURVE + FUNCTION * second_CURVE
          BCin += functionfac*curve2fac;

          //Get the local id of the node to whom the bc is prescribed
          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i]->Id(),discretization.Comm().MyPID());
            exit(1);
          }
        }
        /**
         * For Art_redD_3D_CouplingCond bc
         **/
        else if (ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond"))
        {
          const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond");

          Teuchos::RCP<Teuchos::ParameterList> CoupledTo3DParams  =
            params.get<Teuchos::RCP<Teuchos::ParameterList > >("coupling with 3D fluid params");
          // -----------------------------------------------------------------
          // If the parameter list is empty, then something is wrong!
          // -----------------------------------------------------------------
          if (CoupledTo3DParams.get()==NULL)
          {
            dserror("Cannot prescribe a boundary condition from 3D to reduced D, if the parameters passed don't exist");
            exit(1);
          }

          // -----------------------------------------------------------------
          // Read in Condition type
          // -----------------------------------------------------------------
          //        Type = *(condition->Get<std::string>("CouplingType"));
          // -----------------------------------------------------------------
          // Read in coupling variable rescribed by the 3D simulation
          //
          //     In this case a map called map3D has the following form:
          //     +-----------------------------------------------------------+
          //     |           std::map< string               ,  double        >    |
          //     |     +------------------------------------------------+    |
          //     |     |  ID  | coupling variable name | variable value |    |
          //     |     +------------------------------------------------+    |
          //     |     |  1   |   flow1                |     0.12116    |    |
          //     |     +------+------------------------+----------------+    |
          //     |     |  2   |   pressure2            |    10.23400    |    |
          //     |     +------+------------------------+----------------+    |
          //     |     .  .   .   ....                 .     .......    .    |
          //     |     +------+------------------------+----------------+    |
          //     |     |  N   |   variableN            |    value(N)    |    |
          //     |     +------+------------------------+----------------+    |
          //     +-----------------------------------------------------------+
          // -----------------------------------------------------------------

          int ID = condition->GetInt("ConditionID");
          Teuchos::RCP<std::map<std::string,double> > map3D;
          map3D   = CoupledTo3DParams->get<Teuchos::RCP<std::map<std::string,double > > >("3D map of values");

          // find the applied boundary variable
          std::stringstream stringID;
          stringID<< "_"<<ID;
          for (std::map<std::string,double>::iterator itr = map3D->begin(); itr!=map3D->end(); itr++)
          {
            std::string VariableWithId = itr->first;
            size_t found;
            found= VariableWithId.rfind(stringID.str());
            if (found!=std::string::npos)
            {
              Bc   = std::string(VariableWithId,0,found);
              BCin = itr->second;
              break;
            }
          }
        }
        /**
         * For RedAcinusVentilatorCond bc
         **/
        else if (ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond"))
        {
          DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAcinusVentilatorCond");
          // Get the type of prescribed bc
          Bc  = *(condition->Get<std::string>("phase1"));

          double period  = condition->GetDouble("period");
          double period1 = condition->GetDouble("phase1_period");

          unsigned int phase_number = 0;

          if (fmod(time,period) > period1)
          {
            phase_number = 1;
            Bc = *(condition->Get<std::string>("phase2"));
          }

          const  std::vector<int>*    curve  = condition->Get<std::vector<int> >("curve");
          double curvefac = 1.0;
          const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");

          // Read in the value of the applied BC
          if((*curve)[phase_number]>=0)
          {
            curvefac = DRT::Problem::Instance()->Curve((*curve)[phase_number]).f(time);
            BCin = (*vals)[phase_number]*curvefac;
          }
          else
          {
            dserror("no boundary condition defined!");
            exit(1);
          }

          //Get the local id of the node to whom the bc is prescribed
          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i]->Id(),discretization.Comm().MyPID());
            exit(1);
          }
        }
        else
        {

        }

        /**
         * For pressure or VolumeDependentPleuralPressure bc
         **/
        if (Bc == "pressure" || Bc == "VolumeDependentPleuralPressure")
        {
          Teuchos::RCP<Epetra_Vector> bcval  = params.get<Teuchos::RCP<Epetra_Vector> >("bcval");
          Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector> >("dbctog");

          if (bcval==Teuchos::null||dbctog==Teuchos::null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }

          if (Bc == "VolumeDependentPleuralPressure")
          {
            DRT::Condition * pplCond = ele->Nodes()[i]->GetCondition("RedAirwayVolDependentPleuralPressureCond");
            double Pp_np = 0.0;
            if (pplCond)
            {
              const  std::vector<int>*    curve  = pplCond->Get<std::vector<int>    >("curve");
              double curvefac = 1.0;
              const  std::vector<double>* vals   = pplCond->Get<std::vector<double> >("val");

              //Read in the value of the applied BC
              if((*curve)[0]>=0)
              {
                curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
              }

              //Get parameters for VolumeDependentPleuralPressure condition
              std::string ppl_Type = *(pplCond->Get<std::string>("TYPE"));
              double ap = pplCond->GetDouble("P_PLEURAL_0");
              double bp = pplCond->GetDouble("P_PLEURAL_LIN");
              double cp = pplCond->GetDouble("P_PLEURAL_NONLIN");
              double dp = pplCond->GetDouble("TAU");
              double RV    = pplCond->GetDouble("RV");
              double TLC   = pplCond->GetDouble("TLC");

              //Safety check: in case of polynomial TLC is not used
              if (((ppl_Type == "Linear_Polynomial") or (ppl_Type == "Nonlinear_Polynomial"))
                and (TLC != 0.0))
              {
                dserror("TLC is not used for the following type of VolumeDependentPleuralPressure BC: %s.\n Set TLC = 0.0",ppl_Type.c_str());
              }

              //Safety check: in case of Ogden TLC, P_PLEURAL_0, and P_PLEURAL_LIN
              if ((ppl_Type == "Nonlinear_Ogden")
                and ((TLC != 0.0) or (ap != 0.0) or (bp != 0.0) or (dp == 0.0)))
              {
                dserror("Parameters are not set correctly for Nonlinear_Ogden. Only P_PLEURAL_NONLIN, TAU and RV are used. Set all others to zero. TAU is not allowed to be zero.");
              }

              if (ppl_Type == "Linear_Polynomial")
              {
                const double lungVolumenp = params.get<double>("lungVolume_n");
                Pp_np = ap + bp*(lungVolumenp-RV) + cp*pow((lungVolumenp-RV),dp);
              }
              else if (ppl_Type == "Linear_Exponential")
              {
                const double lungVolumenp = params.get<double>("lungVolume_n");
                const double TLCnp= (lungVolumenp-RV)/(TLC-RV);
                Pp_np = ap +  bp*TLCnp+ cp*exp(dp*TLCnp);
              }
              else if (ppl_Type == "Linear_Ogden")
              {
                const double lungVolumenp = params.get<double>("lungVolume_n");
                Pp_np = RV / lungVolumenp * cp / dp *(1-pow(RV/lungVolumenp,dp));
              }
              else if (ppl_Type == "Nonlinear_Polynomial")
              {
                const double lungVolumenp = params.get<double>("lungVolume_np");
                Pp_np = ap + bp*(lungVolumenp-RV) + cp*pow((lungVolumenp-RV),dp);
              }
              else if (ppl_Type == "Nonlinear_Exponential")
              {
                const double lungVolumenp = params.get<double>("lungVolume_np");
                const double TLCnp= (lungVolumenp-RV)/(TLC-RV);
                Pp_np = ap +  bp*TLCnp+ cp*exp(dp*TLCnp);
              }
              else if (ppl_Type == "Nonlinear_Ogden")
              {
                const double lungVolumenp = params.get<double>("lungVolume_np");
                Pp_np = RV / lungVolumenp * cp / dp *(1-pow(RV/lungVolumenp,dp));
              }
              else
              {
                dserror("Unknown volume pleural pressure type: %s",ppl_Type.c_str());
              }
              Pp_np *= curvefac*((*vals)[0]);
            }
            else
            {
              dserror("No volume dependent pleural pressure condition was defined");
            }
            BCin += Pp_np;
          }

          //Set pressure at node i
          int    gid;
          double val;

          gid = lm[i];
          val = BCin;
          bcval->ReplaceGlobalValues(1,&val,&gid);

          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1,&val,&gid);

        }
        /**
         * For flow bc
         **/
        else if (Bc == "flow")
        {
          // ----------------------------------------------------------
          // Since a node might belong to multiple elements then the
          // flow might be added to the rhs multiple time.
          // To fix this the flow is divided by the number of elements
          // (which is the number of branches). Thus the sum of the
          // final added values is the actual prescribed flow.
          // ----------------------------------------------------------
          int numOfElems = (ele->Nodes()[i])->NumElement();
          BCin /= double(numOfElems);
          rhs(i) += -BCin + rhs(i);
        }
        else
        {
          dserror("prescribed [%s] is not defined for reduced acinuss",Bc.c_str());
          exit(1);
        }

      }
      /**
       * If the node is a terminal node, but no b.c is prescribed to it
       * then a zero output pressure is assumed
       **/
      else
      {
        if (ele->Nodes()[i]->NumElement() == 1)
        {
          //Get the local id of the node to whome the bc is prescribed
          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i],discretization.Comm().MyPID());
            exit(1);
          }

          Teuchos::RCP<Epetra_Vector> bcval  = params.get<Teuchos::RCP<Epetra_Vector> >("bcval");
          Teuchos::RCP<Epetra_Vector> dbctog = params.get<Teuchos::RCP<Epetra_Vector> >("dbctog");

          if (bcval==Teuchos::null||dbctog==Teuchos::null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }

          //Set pressure=0.0 at node i
          int    gid;
          double val;

          gid = lm[i];
          val = 0.0;
          bcval->ReplaceGlobalValues(1,&val,&gid);

          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1,&val,&gid);
        }
      } // END of if there is no BC but the node still is at the terminal

    } // END of if node is available on this processor
  } // End of node i has a condition
}


/*----------------------------------------------------------------------*
 |  Calculate flowrate at current iteration step and the    ismail 01/10|
 |  correponding acinus volume via dV = 0.5*(qnp+qn)*dt                 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::CalcFlowRates(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)

{
  const int elemVecdim = lm.size() ;

  //Get control parameters for time integration
  //Get time-step size
  const double dt = params.get<double>("time step size");
  //Get time
  const double time = params.get<double>("total time");

  //Get control parameters for stabilization and higher-order elements
  // flag for higher order elements
  //  bool higher_order_ele = ele->isHigherOrderElement(distype);

  //Get all general state vectors: flow, pressure,
  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");
  Teuchos::RCP<const Epetra_Vector> pn   = discretization.GetState("pn");
  Teuchos::RCP<const Epetra_Vector> pnm  = discretization.GetState("pnm");

  Teuchos::RCP<Epetra_Vector> qin_nm  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_nm");
  Teuchos::RCP<Epetra_Vector> qin_n   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_n");
  Teuchos::RCP<Epetra_Vector> qin_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  Teuchos::RCP<Epetra_Vector> qout_np = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");
  Teuchos::RCP<Epetra_Vector> qout_n  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_n");
  Teuchos::RCP<Epetra_Vector> qout_nm = params.get<Teuchos::RCP<Epetra_Vector> >("qout_nm");

  Teuchos::RCP<Epetra_Vector> acinar_vn          = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vn");
  Teuchos::RCP<Epetra_Vector> acinar_vnp         = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  Teuchos::RCP<Epetra_Vector> a_volume_strain_np = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp_strain");


  if (pnp==Teuchos::null || pn==Teuchos::null || pnm==Teuchos::null )
    dserror("Cannot get state vectors 'pnp', 'pn', and/or 'pnm''");

  //Extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  //Extract local values from the global vectors
  std::vector<double> mypn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,mypn,lm);

  //Extract local values from the global vectors
  std::vector<double> mypnm(lm.size());
  DRT::UTILS::ExtractMyValues(*pnm,mypnm,lm);

  //Create objects for element arrays
  Epetra_SerialDenseVector epnp(elemVecdim);
  Epetra_SerialDenseVector epn (elemVecdim);
  Epetra_SerialDenseVector epnm(elemVecdim);
  for (int i=0;i<elemVecdim;++i)
  {
    //Split area and volumetric flow rate, insert into element arrays
    epnp(i)   = mypnp[i];
    epn(i)    = mypn[i];
    epnm(i)   = mypnm[i];
  }

  double e_acin_vnp = 0.0;
  double e_acin_vn = 0.0;

  for (int i=0;i<elemVecdim;++i)
  {
    //Split area and volumetric flow rate, insert into element arrays
    e_acin_vnp = (*acinar_vnp)[ele->LID()];
    e_acin_vn  = (*acinar_vn )[ele->LID()];
  }

  //Get the volumetric flow rate from the previous time step
  Teuchos::ParameterList elem_params;
  elem_params.set<double>("qout_np",(*qout_np)[ele->LID()]);
  elem_params.set<double>("qout_n" ,(*qout_n )[ele->LID()]);
  elem_params.set<double>("qout_nm",(*qout_nm)[ele->LID()]);
  elem_params.set<double>("qin_np" ,(*qin_np )[ele->LID()]);
  elem_params.set<double>("qin_n"  ,(*qin_n  )[ele->LID()]);
  elem_params.set<double>("qin_nm" ,(*qin_nm )[ele->LID()]);

  elem_params.set<double>("acin_vnp" ,e_acin_vnp);
  elem_params.set<double>("acin_vn"  ,e_acin_vn );

  Epetra_SerialDenseMatrix sysmat (elemVecdim, elemVecdim,true);
  Epetra_SerialDenseVector  rhs (elemVecdim);

  //Call routine for calculating element matrix and right hand side
  Sysmat(ele,
         epnp,
         epn,
         epnm,
         sysmat,
         rhs,
         material,
         elem_params,
         time,
         dt);

  double qn = (*qin_n  )[ele->LID()];
  double qnp= -1.0*(sysmat(0,0)*epnp(0) + sysmat(0,1)*epnp(1) - rhs(0));

  int gid = ele->Id();

  qin_np  -> ReplaceGlobalValues(1,&qnp,&gid);
  qout_np -> ReplaceGlobalValues(1,&qnp,&gid);

  //Calculate the new volume of the acinus due to the incoming flow; 0.5*(qnp+qn)*dt
  {
    double acinus_volume = e_acin_vn;
    acinus_volume +=  0.5*(qnp+qn)*dt;
    acinar_vnp-> ReplaceGlobalValues(1,& acinus_volume,&gid);

    //Calculate correponding acinar strain
    double vo = 0.0;
    ele->getParams("AcinusVolume",vo);
    double avs_np = (acinus_volume - vo)/vo;
    a_volume_strain_np -> ReplaceGlobalValues(1,&avs_np,&gid);
  }
}


/*----------------------------------------------------------------------*
 |  Calculate element volume                                ismail 07/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::CalcElemVolume(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)

{
  //Get all essential vector variables
  Teuchos::RCP<Epetra_Vector> acinar_vnp   = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  Teuchos::RCP<Epetra_Vector> elemVolumenp = params.get<Teuchos::RCP<Epetra_Vector> >("elemVolumenp");

  //Get acinus size
  double evolnp = (*elemVolumenp)[ele->LID()];

  //Get element global ID
  int gid = ele->Id();

  //Update elem
  elemVolumenp->ReplaceGlobalValues(1,& evolnp,&gid);

}

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::GetCoupledValues(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>   material)
{
  const int   myrank  = discretization.Comm().MyPID();

  //The number of nodes
  const int numnode = lm.size();
  std::vector<int>::iterator it_vcr;

  Teuchos::RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (pnp==Teuchos::null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epnp(numnode);

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)    = mypnp[i];
  }

  // ---------------------------------------------------------------------------------
  // Resolve the BCs
  // ---------------------------------------------------------------------------------
  for(int i = 0; i<ele->NumNode(); i++)
  {
    if (ele->Nodes()[i]->Owner()== myrank)
    {
      if(ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond"))
      {

          const DRT::Condition *condition = ele->Nodes()[i]->GetCondition("Art_redD_3D_CouplingCond");
          Teuchos::RCP<Teuchos::ParameterList> CoupledTo3DParams  =
            params.get<Teuchos::RCP<Teuchos::ParameterList > >("coupling with 3D fluid params");
          // -----------------------------------------------------------------
          // If the parameter list is empty, then something is wrong!
          // -----------------------------------------------------------------
          if (CoupledTo3DParams.get()==NULL)
          {
            dserror("Cannot prescribe a boundary condition from 3D to reduced D, if the parameters passed don't exist");
            exit(1);
          }


        // -----------------------------------------------------------------
        // Compute the variable solved by the reduced D simulation to be
        // passed to the 3D simulation
        //
        //     In this case a map called map1D has the following form:
        //     +-----------------------------------------------------------+
        //     |              std::map< string            ,  double        > >  |
        //     |     +------------------------------------------------+    |
        //     |     |  ID  | coupling variable name | variable value |    |
        //     |     +------------------------------------------------+    |
        //     |     |  1   |   flow1                |     xxxxxxx    |    |
        //     |     +------+------------------------+----------------+    |
        //     |     |  2   |   pressure2            |     xxxxxxx    |    |
        //     |     +------+------------------------+----------------+    |
        //     |     .  .   .   ....                 .     .......    .    |
        //     |     +------+------------------------+----------------+    |
        //     |     |  N   |   variable(N)          | trash value(N) |    |
        //     |     +------+------------------------+----------------+    |
        //     +-----------------------------------------------------------+
        // -----------------------------------------------------------------

        int ID = condition->GetInt("ConditionID");
        Teuchos::RCP<std::map<std::string,double> >  map1D;
        map1D   = CoupledTo3DParams->get<Teuchos::RCP<std::map<std::string,double> > >("reducedD map of values");

        std::string returnedBC = *(condition->Get<std::string>("ReturnedVariable"));

        double BC3d = 0.0;
        if (returnedBC  == "flow")
        {
          // MUST BE DONE
        }
        else if (returnedBC == "pressure")
        {
          BC3d     = epnp(i);
        }
        else
        {
          std::string str = (*condition->Get<std::string>("ReturnedVariable"));
          dserror("%s, is an unimplimented type of coupling",str.c_str());
          exit(1);
        }
        std::stringstream returnedBCwithId;
        returnedBCwithId << returnedBC <<"_" << ID;

        //        std::cout<<"Return ["<<returnedBC<<"] form 1D problem to 3D SURFACE of ID["<<ID<<"]: "<<BC3d<<std::endl;

        // -----------------------------------------------------------------
        // Check whether the coupling wrapper has already initialized this
        // map else wise we will have problems with parallelization, that's
        // because of the preassumption that the map is filled and sorted
        // Thus we can use parallel addition
        // -----------------------------------------------------------------

        std::map<std::string,double>::iterator itrMap1D;
        itrMap1D = map1D->find(returnedBCwithId.str());
        if (itrMap1D == map1D->end())
        {
          dserror("The 3D map for (1D - 3D coupling) has no variable (%s) for ID [%d]",returnedBC.c_str(),ID );
          exit(1);
        }

        // update the 1D map
        (*map1D)[returnedBCwithId.str()] = BC3d;
      }
    } // END of if node is available on this processor
  } // End of node i has a condition
}


/*----------------------------------------------------------------------*
 |  calculate the amount of fluid mixing inside a           ismail 02/13|
 |  junction                                                            |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::GetJunctionVolumeMix(RedAcinus*                   ele,
                                                              Teuchos::ParameterList&      params,
                                                              DRT::Discretization&         discretization,
                                                              Epetra_SerialDenseVector&    volumeMix_np,
                                                              std::vector<int>&            lm,
                                                              Teuchos::RCP<MAT::Material>           material)
{
  // get flow rate out at time step n+1
  Teuchos::RCP<Epetra_Vector> qout_np = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");
  // get acinar volumes at time step n+1
  Teuchos::RCP<Epetra_Vector> acinar_vnp         = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");

  // get the element qout
  double q_out = (*qout_np)[ele->LID()];

  // if transport is flowing into the acinus
  if( q_out >= 0.0)
  {
    ele->getParams("Area",volumeMix_np(1));
  }
  // els if transport is flowing out of the acinus
  else
  {
    ele->getParams("Area",volumeMix_np(0));
    ele->getParams("Area",volumeMix_np(1));
  }

  // extra treatment if an acinus is not connected to anything else
  for (int i=0;i<iel;i++)
  {
    if( ele->Nodes()[i]->NumElement()==1)
      ele->getParams("Area",volumeMix_np(i));
  }
}


/*----------------------------------------------------------------------*
 |  calculate the scalar transport                          ismail 02/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::SolveScatra(RedAcinus*                   ele,
                                                     Teuchos::ParameterList&      params,
                                                     DRT::Discretization&         discretization,
                                                     Epetra_SerialDenseVector&    scatranp,
                                                     Epetra_SerialDenseVector&    volumeMix_np,
                                                     std::vector<int>&            lm,
                                                     Teuchos::RCP<MAT::Material>           material)
{
  const int   myrank  = discretization.Comm().MyPID();
  Teuchos::RCP<Epetra_Vector> qin_np   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  Teuchos::RCP<Epetra_Vector> qout_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");

  Teuchos::RCP<Epetra_Vector> e1scatran = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatran");
  Teuchos::RCP<Epetra_Vector> e2scatran = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatran");

  Teuchos::RCP<Epetra_Vector> e1scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
  Teuchos::RCP<Epetra_Vector> e2scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");

  Teuchos::RCP<const Epetra_Vector> volumeMix = discretization.GetState("junctionVolumeInMix");

  Teuchos::RCP<Epetra_Vector> acinar_vn          = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vn");
  Teuchos::RCP<Epetra_Vector> acinar_vnp         = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");

  double volumenp = (*acinar_vnp)[ele->LID()];
  double volumen = (*acinar_vn )[ele->LID()];

  // extract local values from the global vectors
  std::vector<double> myvolmix(lm.size());
  DRT::UTILS::ExtractMyValues(*volumeMix,myvolmix,lm);
  // get area
  double area = myvolmix[1];

  // get the elements Qin and Qout
  double q_out = (*qout_np)[ele->LID()];
  double q_in  = (*qin_np)[ele->LID()];
  double e1s   = (*e1scatran)[ele->LID()];
  double e2s   = (*e2scatran)[ele->LID()];

  // get time step size
  // const double dt = params.get<double>("time step size");

  // get time
  const double time = params.get<double>("total time");

  //--------------------------------------------------------------------
  // get element length
  //--------------------------------------------------------------------

  // evaluate velocity at nodes (1) and (2)
  double vel1 = q_in /area;
  double vel2 = q_out/area;

  LINALG::Matrix<2,1> velv;
  velv(0,0)=vel1;
  velv(1,0)=vel2;
  // get average velocity
  double vel = vel2;

  //--------------------------------------------------------------------
  // if vel>=0 then node(2) is analytically evaluated;
  //  ---> node(1) is either prescribed or comes from the junction
  // if vel< 0 then node(1) is analytically evaluated;
  //  ---> node(2) is either prescribed or comes from the junction
  //--------------------------------------------------------------------

  if (vel >=0.0)
  {
    // extrapolate the analytical solution
    int gid = ele->Id();
    double scnp = 0.0;
    scnp = (e2s*volumen + e1s*(volumenp-volumen))/(volumenp);
    e2scatranp->ReplaceGlobalValues(1,&scnp,&gid);
  }
  else
  {
    // extrapolate the analytical solution
    double scnp = 0.0;
    int gid = ele->Id();
    scnp = (e2s*volumen + e2s*(volumenp-volumen))/(volumenp);
    {
      e2scatranp->ReplaceGlobalValues(1,&scnp,&gid);
      e1scatranp->ReplaceGlobalValues(1,&scnp,&gid);
    }
  }

  for (int i=0;i<2;i++)
  {
    if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedScatraCond") &&     myrank == ele->Nodes()[i]->Owner())
    {
      double scnp =0.0;
      DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedScatraCond");
      // Get the type of prescribed bc

      const  std::vector<int>*    curve  = condition->Get<std::vector<int>    >("curve");
      double curvefac = 1.0;
      const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");

      // -----------------------------------------------------------------
      // Read in the value of the applied BC
      // -----------------------------------------------------------------
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      if (curvenum>=0 )
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

      scnp = (*vals)[0]*curvefac;

      const std::vector<int>*    functions = condition->Get<std::vector<int> >("funct");
      int functnum = -1;
      if (functions) functnum = (*functions)[0];
      else functnum = -1;

      double functionfac = 0.0;
      if(functnum>0)
      {
        functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,(ele->Nodes()[i])->X(),time,NULL);
      }
      scnp += functionfac;

      // ----------------------------------------------------
      // convert O2 saturation to O2 concentration
      // ----------------------------------------------------
      // get O2 properties in air
      int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
      // check if O2 properties material exists
      if (id==-1)
      {
        dserror("A material defining O2 properties in air could not be found");
        exit(1);
      }
      const MAT::PAR::Parameter* smat = DRT::Problem::Instance()->Materials()->ParameterById(id);
      const MAT::PAR::Air_0d_O2_saturation* actmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(smat);

      // get atmospheric pressure
      double patm = actmat->atmospheric_p_;
      // get number of O2 moles per unit volume of O2
      double nO2perVO2 = actmat->nO2_per_VO2_;
      // calculate the PO2 at nodes
      double pO2 = scnp*patm;
      // calculate VO2
      double vO2 = volumenp*(pO2/patm);
      // evaluate initial concentration
      scnp = nO2perVO2*vO2/volumenp;
      //------------
      if(i==0)
      {
        int    gid = ele->Id();
        double val = scnp;
        if (vel<0.0)
          val = (*e1scatranp)[ele->LID()];
//        if (ele->Owner()==myrank)
        {
          e1scatranp->ReplaceGlobalValues(1,&val,&gid);
        }
        scatranp(0) = val*area;
      }
      else
      {
        int    gid = ele->Id();
        double val = scnp;
        if (vel>=0.0)
          val = (*e2scatranp)[ele->LID()];
//        if (ele->Owner()==myrank)
        {
          e2scatranp->ReplaceGlobalValues(1,&val,&gid);
        }
        scatranp(1) = val*area;
      }
    }
  }


  {
    scatranp(1) = (*e2scatranp)[ele->LID()]*area;
  }
  if (vel <  0.0)
  {
    scatranp(0) = (*e1scatranp)[ele->LID()]*area;
  }
}//SolveScatra



/*----------------------------------------------------------------------*
 |  calculate the scalar transport                          ismail 02/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::SolveScatraBifurcations(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  Epetra_SerialDenseVector&    scatranp,
  Epetra_SerialDenseVector&    volumeMix_np,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{

  Teuchos::RCP<Epetra_Vector> qin_np   = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");
  Teuchos::RCP<Epetra_Vector> qout_np  = params.get<Teuchos::RCP<Epetra_Vector> >("qout_np");

  Teuchos::RCP<const Epetra_Vector> scatran  = discretization.GetState("scatranp");

  Teuchos::RCP<Epetra_Vector> e1scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
  Teuchos::RCP<Epetra_Vector> e2scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");

  Teuchos::RCP<const Epetra_Vector> volumeMix = discretization.GetState("junctionVolumeInMix");

  // extract local values from the global vectors
  std::vector<double> myvolmix(lm.size());
  DRT::UTILS::ExtractMyValues(*volumeMix,myvolmix,lm);
  // get area
  double area = myvolmix[1];

  // get the elements Qin and Qout
  double q_out = (*qout_np)[ele->LID()];
  double q_in  = (*qin_np)[ele->LID()];

  // extract local values from the global vectors
  std::vector<double> myscatran(lm.size());
  DRT::UTILS::ExtractMyValues(*scatran,myscatran,lm);

  // evaluate velocity at nodes (1) and (2)
  double vel1 = q_in /area;
  double vel2 = q_out/area;

  LINALG::Matrix<2,1> velv;
  velv(0,0)=vel1;
  velv(1,0)=vel2;
  // get average velocity
  double vel = 0.5*(vel1+vel2);

  //--------------------------------------------------------------------
  // if vel>=0 then node(2) is analytically evaluated;
  //  ---> node(1) is either prescribed or comes from the junction
  // if vel< 0 then node(1) is analytically evaluated;
  //  ---> node(2) is either prescribed or comes from the junction
  //--------------------------------------------------------------------
  if (vel >=0.0)
  {
    // extrapolate the analytical solution
    double scnp = myscatran[0];
    int gid = ele->Id();
    e1scatranp->ReplaceGlobalValues(1,&scnp,&gid);
  }
  else
  {
    // extrapolate the analytical solution
    double scnp = myscatran[1];
    int gid = ele->Id();
    e2scatranp->ReplaceGlobalValues(1,&scnp,&gid);
  }
}//SolveScatraBifurcations


/*----------------------------------------------------------------------*
 |  calculate the scalar transport                          ismail 02/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::UpdateScatra(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{

  const int   myrank  = discretization.Comm().MyPID();

  Teuchos::RCP<const Epetra_Vector> dscatranp  = discretization.GetState("dscatranp");
  Teuchos::RCP<Epetra_Vector> dscatranp_m = params.get<Teuchos::RCP<Epetra_Vector> >("dscatranp");
  Teuchos::RCP<Epetra_Vector> qin_np     = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");

  // get flowrate
  double qin = (*qin_np)[ele->LID()];

  // extract local values from the global vectors
  std::vector<double> mydscatra(lm.size());
  DRT::UTILS::ExtractMyValues(*dscatranp,mydscatra,lm);

  //--------------------------------------------------------------------
  // if vel>=0 then node(2) is analytically evaluated;
  //  ---> node(1) is either prescribed or comes from the junction
  // if vel< 0 then node(1) is analytically evaluated;
  //  ---> node(2) is either prescribed or comes from the junction
  //--------------------------------------------------------------------
  if (qin<0.0)
  {
    int gid = lm[1];
    double val = mydscatra[1];
    if (myrank == ele->Nodes()[1]->Owner())
    {
      dscatranp_m->ReplaceGlobalValues(1,&val,&gid);
    }
  }
}//UpdateScatra


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::UpdateElem12Scatra(
  RedAcinus*                  ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
  Teuchos::RCP<const Epetra_Vector> scatranp = discretization.GetState("scatranp");
  Teuchos::RCP<const Epetra_Vector> dscatranp = discretization.GetState("dscatranp");
  Teuchos::RCP<const Epetra_Vector> volumeMix = discretization.GetState("junctionVolumeInMix");

  Teuchos::RCP<Epetra_Vector> qin_np     = params.get<Teuchos::RCP<Epetra_Vector> >("qin_np");
  Teuchos::RCP<Epetra_Vector> e1scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e1scatranp");
  Teuchos::RCP<Epetra_Vector> e2scatranp = params.get<Teuchos::RCP<Epetra_Vector> >("e2scatranp");

  // extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*scatranp,myscatranp,lm);

  // extract local values from the global vectors
  std::vector<double> mydscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*dscatranp,mydscatranp,lm);

  // extract local values from the global vectors
  std::vector<double> myvolmix(lm.size());
  DRT::UTILS::ExtractMyValues(*volumeMix,myvolmix,lm);

  // get flowrate
  double qin = (*qin_np)[ele->LID()];
  // Get the average concentration

  // ---------------------------------------------------------------------
  // element scatra must be updated only at the capillary nodes.
  // ---------------------------------------------------------------------
  //  double e2s = (*e2scatranp)[ele->LID()] + mydscatranp[1]*myvolmix[1];
  double e2s = myscatranp[1];

  int gid = ele->Id();
  e2scatranp->ReplaceGlobalValues(1,&e2s,&gid);
  if (qin<0.0)
  {
    e1scatranp->ReplaceGlobalValues(1,&e2s,&gid);
  }
}



/*----------------------------------------------------------------------*
 |  calculate PO2 from concentration                        ismail 06/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::EvalPO2FromScatra(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
  const int   myrank  = discretization.Comm().MyPID();

  // get Po2 vector
  Teuchos::RCP<Epetra_Vector> po2        = params.get<Teuchos::RCP<Epetra_Vector> >("PO2");
  // get Po2 vector
  Teuchos::RCP<const Epetra_Vector> scatran = discretization.GetState("scatranp");
  // get acinar volume
  Teuchos::RCP<Epetra_Vector> acinar_vnp = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_vnp");
  // -------------------------------------------------------------------
  // extract scatra values
  // -------------------------------------------------------------------
  // extract local values from the global vectors
  std::vector<double> myscatran(lm.size());
  DRT::UTILS::ExtractMyValues(*scatran,myscatran,lm);

  // -------------------------------------------------------------------
  // find out if the material type is Air or Blood
  // -------------------------------------------------------------------
  std::string fluidType = "none";
  // if RedAirwayScatraAirCond then material type is air
  if (ele->Nodes()[0]->GetCondition("RedAirwayScatraAirCond")!=NULL
    && ele->Nodes()[1]->GetCondition("RedAirwayScatraAirCond")!=NULL)
  {
    fluidType = "air";
  }
  else
  {
    dserror("A scalar transport element must be defined either as \"air\"");
    exit(1);
  }

  // define a empty pO2 vector
  double pO2 = 0.0;

  // -------------------------------------------------------------------
  // Get O2 properties in air
  // -------------------------------------------------------------------
  if (fluidType == "air")
  {
    // -----------------------------------------------------------------
    // Get O2 properties in air
    // -----------------------------------------------------------------

    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_0d_o2_air_saturation);
    // check if O2 properties material exists
    if (id==-1)
    {
      dserror("A material defining O2 properties in air could not be found");
      exit(1);
    }
    const MAT::PAR::Parameter* smat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::Air_0d_O2_saturation* actmat = static_cast<const MAT::PAR::Air_0d_O2_saturation*>(smat);

    // get atmospheric pressure
    double patm = actmat->atmospheric_p_;
    // get number of O2 moles per unit volume of O2
    double nO2perVO2 = actmat->nO2_per_VO2_;

    // -----------------------------------------------------------------
    // Calculate Vo2 in air
    // -----------------------------------------------------------------
    // get airway volume
    double vAir = (*acinar_vnp)[ele->LID()];
    // calculate the VO2 at nodes
    double vO2 = (vAir*myscatran[lm.size()-1])/nO2perVO2;
    // calculate PO2 at nodes
    pO2 = patm*vO2/vAir;
  }
  else
  {
    dserror("A scalar transport element must be defined either as \"air\" or \"blood\"");
    exit(1);
  }

  // -------------------------------------------------------------------
  // Set element pO2 to PO2 vector
  // -------------------------------------------------------------------
  int    gid = lm[lm.size()-1];
  double val = pO2;
  if(myrank == ele->Nodes()[lm.size()-1]->Owner())
  {
    po2->ReplaceGlobalValues(1,&val,&gid);
  }

}// EvalPO2FromScatra


/*----------------------------------------------------------------------*
 |  calculate essential nodal values                        ismail 06/13|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcinusImpl<distype>::EvalNodalEssentialValues(
  RedAcinus*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  Epetra_SerialDenseVector&    nodal_surface,
  Epetra_SerialDenseVector&    nodal_volume,
  Epetra_SerialDenseVector&    nodal_avg_scatra,
  std::vector<int>&            lm,
  Teuchos::RCP<MAT::Material>           material)
{
  //Get all general state vectors: flow, pressure,
  Teuchos::RCP<Epetra_Vector> acinar_e_v  = params.get<Teuchos::RCP<Epetra_Vector> >("acinar_v");
  Teuchos::RCP<const Epetra_Vector> scatranp = discretization.GetState("scatranp");

  //Extract scatra values
  //Extract local values from the global vectors
  std::vector<double> myscatranp(lm.size());
  DRT::UTILS::ExtractMyValues(*scatranp,myscatranp,lm);

  //Find the volume of an acinus
  //Get the current acinar volume
  double volAcinus = (*acinar_e_v)[ele->LID()];
  //Set nodal volume
  nodal_volume[1] = volAcinus;

  //Find the average scalar transport concentration
  //Set nodal flowrate
  nodal_avg_scatra[0] = myscatranp[1];
  nodal_avg_scatra[1] = myscatranp[1];

  //Find the total gas exchange surface inside an acinus
  // get the initial volume of an acinus
  double volAcinus0;
  ele->getParams("AcinusVolume",volAcinus0);
  // get the initial volume of an alveolar duct
  double volAlvDuct0;
  ele->getParams("AlveolarDuctVolume",volAlvDuct0);
  // find the  number of alveolar duct
  const double numOfAlvDucts = double(floor(volAcinus0/volAlvDuct0));
  // define the number of alveoli per alveolar duct
  const double nAlveoliPerAlveolarDuct = 36.0;
  // define the number of alveoli per duct
  const double nAlveoliPerDuct = 4.0;
  // find the volume of one alveolus
  const double volAlveolus = volAcinus/numOfAlvDucts/nAlveoliPerAlveolarDuct;
  // find the surface of one alveolus
  const double surfAlveolus = (6.0+12.0*sqrt(3.0))*pow(volAlveolus/(8.0*sqrt(2.0)),2.0/3.0);
  // get the length of an edge of an alveolus
  const double al = pow(volAlveolus/(8.0*sqrt(2.0)),1.0/3.0)/3.0;
  // find the surface of an alveolar duct
  const double surfAlveolarDuct = (nAlveoliPerAlveolarDuct - 2.0*nAlveoliPerDuct)*surfAlveolus + 6.0*(al*al);
  // find the surface of an acinus
  const double surfAcinus = surfAlveolarDuct*numOfAlvDucts;
  // set nodal surface area
  nodal_surface[1] = surfAcinus;
}
