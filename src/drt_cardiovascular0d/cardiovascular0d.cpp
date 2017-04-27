/*!----------------------------------------------------------------------
\file cardiovascular0d.cpp

\brief Monolithic coupling of 3D structural dynamics and 0D cardiovascular flow models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>
*----------------------------------------------------------------------*/

#include <iostream>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_so3/so_surface.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "cardiovascular0d.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                              mhv 10/13|
 *----------------------------------------------------------------------*/
UTILS::Cardiovascular0D::Cardiovascular0D(Teuchos::RCP<DRT::Discretization> discr,
    const std::string& conditionname,
    std::vector<int>& curID):
    actdisc_(discr),
    cardiovascular0dcond_(0),
    cardiovascular0dstructcoupcond_(0),
    cardiovascular0dtype_(none),
    atrium_model_(DRT::INPUT::IntegralValue<INPAR::CARDIOVASCULAR0D::Cardvasc0DAtriumModel>(DRT::Problem::Instance()->Cardiovascular0DStructuralParams().sublist("SYS-PUL CIRCULATION PARAMETERS"),"ATRIUM_MODEL")),
    ventricle_model_(DRT::INPUT::IntegralValue<INPAR::CARDIOVASCULAR0D::Cardvasc0DVentricleModel>(DRT::Problem::Instance()->Cardiovascular0DStructuralParams().sublist("SYS-PUL CIRCULATION PARAMETERS"),"VENTRICLE_MODEL")),
    respiratory_model_(DRT::INPUT::IntegralValue<INPAR::CARDIOVASCULAR0D::Cardvasc0DRespiratoryModel>(DRT::Problem::Instance()->Cardiovascular0DStructuralParams().sublist("RESPIRATORY PARAMETERS"),"RESPIRATORY_MODEL")),
    gaussrule_(DRT::UTILS::intrule2D_undefined)
{
  actdisc_->GetCondition(conditionname,cardiovascular0dcond_);
  if (cardiovascular0dcond_.size())
  {
    cardiovascular0dtype_=GetCardiovascular0DType(conditionname);
    std::vector<int> curcoupID;
    for (unsigned int i=0; i<cardiovascular0dcond_.size();i++)
    {
      //std::vector<int> curID(i);
      curID.push_back(cardiovascular0dcond_[i]->GetInt("id"));
    }

    std::vector<DRT::Condition*> surfneumcond;
    //std::vector<int> tmp;
    Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
    if (structdis == Teuchos::null) dserror("no structure discretization available");

    // first get all Neumann conditions on structure
    structdis->GetCondition("SurfaceNeumannCardiovascular0D",surfneumcond);
    unsigned int numneumcond = surfneumcond.size();
    if (numneumcond == 0) dserror("no Neumann conditions on structure");

    // now filter those Neumann conditions that are due to the coupling
    //std::vector<DRT::Condition*> coupcond;
    for (unsigned int k = 0; k < numneumcond; ++k)
    {
      DRT::Condition* actcond = surfneumcond[k];
      if (actcond->Type() == DRT::Condition::Cardiovascular0DStructureCoupling)
        cardiovascular0dstructcoupcond_.push_back(actcond);
    }
    unsigned int numcond = cardiovascular0dstructcoupcond_.size();

    if (numcond == 0) dserror("no coupling conditions found");

    //if (cardiovascular0dcond_.size() != cardiovascular0dstructcoupcond_.size()) dserror("Coupling conditions do not match cardiovascular0d conditions!");

    std::vector<int> wkID(cardiovascular0dcond_.size());
    for (unsigned int i=0; i<cardiovascular0dcond_.size(); i++)
    {
      wkID[i]=(cardiovascular0dcond_[i])->GetInt("id");
    }

    // safety checks for closed-loop vascular model
    if(cardiovascular0dtype_ == cardvasc0d_syspulcirculation or cardiovascular0dtype_ == cardvascrespir0d_syspulperiphcirculation)
    {
      std::vector<const std::string*> condtype(cardiovascular0dcond_.size());
      for (unsigned int i=0; i<cardiovascular0dcond_.size(); i++)
      {
        condtype[i] = cardiovascular0dcond_[i]->Get<std::string>("type");

        if (atrium_model_ == INPAR::CARDIOVASCULAR0D::atr_elastance_0d or atrium_model_ == INPAR::CARDIOVASCULAR0D::atr_prescribed)
        {
          if (*condtype[i] == "atrium_left" or *condtype[i] == "atrium_right")
            dserror("Set ATRIUM_MODEL to '3D' if you want to couple the 0D vascular system to a 3D atrial structure!");
        }
      }

      switch (atrium_model_)
      {
        case INPAR::CARDIOVASCULAR0D::atr_elastance_0d:
        case INPAR::CARDIOVASCULAR0D::atr_prescribed:
        {
          if (ventricle_model_ == INPAR::CARDIOVASCULAR0D::ventr_structure_3d)
          {
            if (cardiovascular0dcond_.size() == 2)
            {

              if (*condtype[0] != "ventricle_left" and *condtype[1] != "ventricle_left")
                dserror("No left/right ventricle type of condition specified!");
              if (*condtype[0] != "ventricle_right" and *condtype[1] != "ventricle_right")
                dserror("No left/right ventricle type of condition specified!");
            }
            else
              dserror("You need 2 conditions (left + right ventricle)!");
          }
          if (ventricle_model_ == INPAR::CARDIOVASCULAR0D::ventr_elastance_0d or ventricle_model_ == INPAR::CARDIOVASCULAR0D::ventr_prescribed)
          {
            if (cardiovascular0dcond_.size() == 1)
            {
              if (*condtype[0] != "dummy")
                dserror("Only specify 1 dummy condition!");
            }
            else
              dserror("You're only allowed to specify 1 (dummy) condition!");
          }
        }
        break;
        case INPAR::CARDIOVASCULAR0D::atr_structure_3d:
        {
          if (ventricle_model_ == INPAR::CARDIOVASCULAR0D::ventr_elastance_0d)
            dserror("You cannot use 3D atria with 0D ventricles!");

          if (cardiovascular0dcond_.size() == 4)
          {
            if (*condtype[0] != "atrium_left" and *condtype[1] != "atrium_left" and *condtype[2] != "atrium_left" and *condtype[3] != "atrium_left")
              dserror("ATRIUM_MODEL is set to '3D' but you don't have a left/right atrium type of condition specified!");
            if (*condtype[0] != "atrium_right" and *condtype[1] != "atrium_right" and *condtype[2] != "atrium_right" and *condtype[3] != "atrium_right")
              dserror("ATRIUM_MODEL is set to '3D' but you don't have a left/right atrium type of condition specified!");

            if (*condtype[0] != "ventricle_left" and *condtype[1] != "ventricle_left" and *condtype[2] != "ventricle_left" and *condtype[3] != "ventricle_left")
              dserror("No left/right ventricle type of condition specified!");
            if (*condtype[0] != "ventricle_right" and *condtype[1] != "ventricle_right" and *condtype[2] != "ventricle_right" and *condtype[3] != "ventricle_right")
              dserror("No left/right ventricle type of condition specified!");
          }
          else
            dserror("You need 4 conditions (left + right ventricle, and left + right atrium)!");
        }
        break;
        default:
          dserror("Unknown ATRIUM_MODEL!");
        break;
      } // end of switch
    } // end if (cardiovascular0dtype_ == cardvasc0d_syspulcirculation)

    std::vector<int> coupcondID(cardiovascular0dstructcoupcond_.size());
    //set Neumann line to condition
    for (unsigned int i=0; i<cardiovascular0dstructcoupcond_.size(); i++)
    {
      coupcondID[i]=(cardiovascular0dstructcoupcond_[i])->GetInt("coupling_id");

      cardiovascular0dstructcoupcond_[i]->Add("type","neum_orthopressure");
      std::vector<int> onoff(6,0);
      onoff[0] = 1;
      cardiovascular0dstructcoupcond_[i]->Add("onoff",onoff);
      std::vector<double> val(6,0.0);
      cardiovascular0dstructcoupcond_[i]->Add("val",val);

    }

    if(cardiovascular0dcond_.size() != cardiovascular0dstructcoupcond_.size())
      dserror("Number of cardiovascular 0D conditions has to be equal to number of cardiovascular 0d structure coupling conditions!");

    if (*std::min_element(wkID.begin(),wkID.end()) != 0)
      dserror("Start your ID numbering from 0 on!");
    if (*std::min_element(coupcondID.begin(),coupcondID.end()) != 0)
      dserror("Start your coupling_id numbering from 0 on!");

    if (*std::min_element(wkID.begin(),wkID.end()) != *std::min_element(coupcondID.begin(),coupcondID.end()))
      dserror("Min cardiovascular0d id not equal to min cardiovascular0d structure coupling id!");
    if (*std::max_element(wkID.begin(),wkID.end()) != *std::max_element(coupcondID.begin(),coupcondID.end()))
      dserror("Max cardiovascular0d id not equal to max cardiovascular0d structure coupling id!");

    if (*std::max_element(wkID.begin(),wkID.end()) != static_cast<int>(cardiovascular0dcond_.size())-1)
      dserror("Max ID should be the number of conditions minus 1!");
    if (*std::max_element(coupcondID.begin(),coupcondID.end()) != static_cast<int>(cardiovascular0dstructcoupcond_.size())-1)
      dserror("Max coupling_id should be the number of conditions minus 1!");

  }
  else
  {
    cardiovascular0dtype_=none;
  }
}

/*-----------------------------------------------------------------------*
|(private)                                                      mhv 10/13|
 *-----------------------------------------------------------------------*/
UTILS::Cardiovascular0D::Cardiovascular0DType UTILS::Cardiovascular0D::GetCardiovascular0DType(const std::string& name)
{
  if (name=="Cardiovascular0D4ElementWindkesselStructureCond")
    return cardvasc0d_4elementwindkessel;
  else if (name=="Cardiovascular0DArterialProxDistStructureCond")
    return cardvasc0d_arterialproxdist;
  else if (name=="Cardiovascular0DSysPulCirculationStructureCond")
    return cardvasc0d_syspulcirculation;
  else if (name=="CardiovascularRespiratory0DSysPulPeriphCirculationStructureCond")
    return cardvascrespir0d_syspulperiphcirculation;
  return none;
}

/*------------------------------------------------------------------------*
|(public)                                                      mhv 10/13  |
|Initialization routine computes ref base values and activates conditions |
 *------------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::Initialize(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2)
{
  dserror("Overridden by derived class!");
  return;
}


/*-----------------------------------------------------------------------*
|(public)                                                       mhv 10/13|
|Evaluate Cardiovascular0D functions, choose the right action based on type    |
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::Evaluate(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat1,
    Teuchos::RCP<LINALG::SparseOperator> sysmat2,
    Teuchos::RCP<LINALG::SparseOperator> sysmat3,
    Teuchos::RCP<Epetra_Vector>    sysvec1,
    Teuchos::RCP<Epetra_Vector>    sysvec2,
    Teuchos::RCP<Epetra_Vector>    sysvec3,
    Teuchos::RCP<Epetra_Vector>    sysvec4,
    Teuchos::RCP<Epetra_Vector>    sysvec5
    )
{
  dserror("Overridden by derived class!");
  return;
}


void UTILS::Cardiovascular0D::EvaluateDStructDp(
    Teuchos::ParameterList&        params,
    Teuchos::RCP<LINALG::SparseOperator> sysmat
    )
{
  // get structural time-integrator dependent values
  double sc_strtimint = params.get("scale_timint",1.0);

  int numdof_per_cond = 1;

  // choose action
  switch (cardiovascular0dtype_)
  {
  case cardvasc0d_4elementwindkessel:
    numdof_per_cond = 3;
    break;
  case cardvasc0d_arterialproxdist:
    numdof_per_cond = 4;
    break;
  case cardvasc0d_syspulcirculation:
    break;
  case cardvascrespir0d_syspulperiphcirculation:
    break;
  case none:
    return;
  default:
    dserror("Unknown Cardiovascular0D type to be evaluated in Cardiovascular0D class!");
    break;
  }

  const int offsetID = params.get<int>("OffsetID");

  std::vector<int> gindex_syspulcirculation(16);
  gindex_syspulcirculation[0] = offsetID;
  for (int j = 1; j < 16; j++) gindex_syspulcirculation[j] = gindex_syspulcirculation[0]+j;

  std::vector<int> gindex_syspulperiphcirculation(34);
  gindex_syspulperiphcirculation[0] = offsetID;
  for (int j = 1; j < 34; j++) gindex_syspulperiphcirculation[j] = gindex_syspulperiphcirculation[0]+j;

  // loop over cardiovascular0d structure coupling conditions
  /* here we do tge loop to assemble the offdiagonal stiffness block dfext/dcvdof (0,1 block)
  this is the derivative of the orthopressure Neumann load (external load vector fext) w.r.t. the pressure*/
  for (unsigned int i = 0; i < cardiovascular0dstructcoupcond_.size(); ++i)
  {
    DRT::Condition& coupcond = *(cardiovascular0dstructcoupcond_[i]);

    int coupcondID=coupcond.GetInt("coupling_id");
    params.set("coupling_id",coupcondID);

    Teuchos::RCP<const Epetra_Vector> disp = params.get<Teuchos::RCP<const Epetra_Vector> >("new disp");
    actdisc_->SetState("displacement",disp);

    // global and local ID of this bc in the redundant vectors
    std::vector<int> gindex(numdof_per_cond);
    gindex[0] = numdof_per_cond*coupcondID + offsetID;
    for (int j = 1; j < numdof_per_cond; j++) gindex[j] = gindex[0]+j;

    std::map<int,Teuchos::RCP<DRT::Element> >& geom = coupcond.Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      curr->second->LocationVector(*actdisc_,lm,lmowner,lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      Epetra_SerialDenseVector elevector;
      elevector.Size(eledim);

      DRT::Element* element = curr->second.get();
      int numnode=element->NumNode();

      // allocate vector for shape functions and matrix for derivatives
      LINALG::SerialDenseVector funct(numnode);
      LINALG::SerialDenseMatrix deriv(2,numnode);
      LINALG::SerialDenseMatrix xc;

      xc.LightShape(numnode,3);

      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
      Teuchos::RCP<const Epetra_Vector> curdispl = actdisc_->GetState("displacement");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*curdispl,mydisp,lm);

      for (int j=0; j<numnode; ++j)
      {
        xc(j,0) = element->Nodes()[j]->X()[0] + mydisp[j*3+0];
        xc(j,1) = element->Nodes()[j]->X()[1] + mydisp[j*3+1];
        xc(j,2) = element->Nodes()[j]->X()[2] + mydisp[j*3+2];
      }

      /*----------------------------------------------------------------------*
      |               start loop over integration points                     |
      *----------------------------------------------------------------------*/
      DRT::Element::DiscretizationType shape = element->Shape();
      // type of gaussian integration
      switch(shape)
      {
      case DRT::Element::tri3:
        gaussrule_ = DRT::UTILS::intrule_tri_3point;
      break;
      case DRT::Element::tri6:
        gaussrule_ = DRT::UTILS::intrule_tri_6point;
      break;
      case DRT::Element::quad4:
        gaussrule_ = DRT::UTILS::intrule_quad_4point;
      break;
      case DRT::Element::quad8:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::quad9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      case DRT::Element::nurbs9:
        gaussrule_ = DRT::UTILS::intrule_quad_9point;
      break;
      default:
          dserror("shape type unknown!\n");
          break;
      }

      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        // set gausspoints from integration rule
        Epetra_SerialDenseVector e(2);
        e(0) = intpoints.qxg[gp][0];
        e(1) = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element

        DRT::UTILS::shape_function_2D(funct,e(0),e(1),shape);
        DRT::UTILS::shape_function_2D_deriv1(deriv,e(0),e(1),shape);

        //stuff to get spatial Neumann
        const int numdim = 3;
        LINALG::SerialDenseMatrix gp_coord(1,numdim);

        std::vector<double> normal(3);

        // note that the length of this normal is the area dA
        // compute dXYZ / drs
        LINALG::SerialDenseMatrix dxyzdrs(2,3);
        dxyzdrs.Multiply('N','N',1.0,deriv,xc,0.0);

        normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
        normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
        normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

        const double fac = intpoints.qwgt[gp];
        for (int node=0; node < numnode; ++node)
          for(int dim=0 ; dim<3; dim++)
            elevector[node*3+dim] += funct[node] * normal[dim] * fac;

      }

      int eid = curr->second->Id();

      // assemble the offdiagonal stiffness block (0,1 block) arising from dR_struct/dcvdof
      // assemble to rectangular matrix. The col corresponds to the Cardiovascular0D ID.
      std::vector<int> colvec(1);

      // choose action
      switch (cardiovascular0dtype_)
      {
      case cardvasc0d_4elementwindkessel:
        colvec[0]=gindex[0];
        break;
      case cardvasc0d_arterialproxdist:
        colvec[0]=gindex[0];
        break;
      case cardvasc0d_syspulcirculation:
      {
        for (unsigned int j = 0; j < cardiovascular0dcond_.size(); ++j)
        {
          DRT::Condition& cond = *(cardiovascular0dcond_[j]);
          int id_cardvasc0d = cond.GetInt("id");
          if (coupcondID == id_cardvasc0d)
          {
            // get the type of the corresponding cardiovascular0D condition
            const std::string* conditiontype = cardiovascular0dcond_[j]->Get<std::string>("type");
            if (*conditiontype == "ventricle_left") colvec[0]=gindex_syspulcirculation[3];
            if (*conditiontype == "ventricle_right") colvec[0]=gindex_syspulcirculation[11];
            if (*conditiontype == "atrium_left") colvec[0]=gindex_syspulcirculation[0];
            if (*conditiontype == "atrium_right") colvec[0]=gindex_syspulcirculation[8];
            if (*conditiontype == "dummy") colvec[0]=gindex_syspulcirculation[0];
          }
        }
      }
        break;
      case cardvascrespir0d_syspulperiphcirculation:
      {
        for (unsigned int j = 0; j < cardiovascular0dcond_.size(); ++j)
        {
          DRT::Condition& cond = *(cardiovascular0dcond_[j]);
          int id_cardvasc0d = cond.GetInt("id");
          if (coupcondID == id_cardvasc0d)
          {
            // get the type of the corresponding cardiovascular0D condition
            const std::string* conditiontype = cardiovascular0dcond_[j]->Get<std::string>("type");
            if (*conditiontype == "ventricle_left") colvec[0]=gindex_syspulperiphcirculation[3];
            if (*conditiontype == "ventricle_right") colvec[0]=gindex_syspulperiphcirculation[27];
            if (*conditiontype == "atrium_left") colvec[0]=gindex_syspulperiphcirculation[0];
            if (*conditiontype == "atrium_right") colvec[0]=gindex_syspulperiphcirculation[24];
            if (*conditiontype == "dummy") colvec[0]=gindex_syspulperiphcirculation[0];
          }
        }
      }
        break;
      case none:
        return;
      default:
        dserror("Unknown Cardiovascular0D type to be evaluated in Cardiovascular0D class!");
        break;
      }

      elevector.Scale(sc_strtimint);
      sysmat->Assemble(eid,lmstride,elevector,lm,lmowner,colvec);
    }

  }

  return;
}


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void UTILS::Cardiovascular0D::SetState
(
    const std::string& state,  ///< name of state to set
    Teuchos::RCP<Epetra_Vector> V  ///< values to set
)
{
  actdisc_->SetState(state,V);
}

