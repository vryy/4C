/*!----------------------------------------------------------------------
\file boxfilter.cpp
//
\brief Filter routines for dynamic Smagorinsky and dynamic Vreman model

Documentation see header.

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/


#include "boxfilter.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/newtonianfluid.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     krank 09/13|
 *----------------------------------------------------------------------*/
FLD::Boxfilter::Boxfilter(
  Teuchos::RCP<DRT::Discretization>     actdis,
  Teuchos::ParameterList&      params)
  :
  // call constructor for "nontrivial" objects
  discret_            (actdis             ),
  params_             (params             ),
  physicaltype_       (DRT::INPUT::get<INPAR::FLUID::PhysicalType>(params_, "Physical Type")),
  //  available control settings
  apply_dynamic_smagorinsky_(false),
  vreman_dynamic_(false),
  apply_box_filter_(false),
  loma_(false),
  incomp_(false),
  velocity_(false),
  reynoldsstress_(false),
  modeled_subgrid_stress_(false),
  expression_(false),
  strainrate_(false),
  alphaij_(false),
  alpha2_(false),
  finescale_velocity_(false),
  densvelocity_(false),
  densstrainrate_(false),
  density_(false),
  phi_(false),
  phi2_(false),
  phiexpression_(false),
  alphaijsc_(false)
{
  Teuchos::ParameterList *  modelparams =&(params_.sublist("TURBULENCE MODEL"));

  if (modelparams->get<std::string>("TURBULENCE_APPROACH","DNS_OR_RESVMM_LES") == "CLASSICAL_LES")
  {
    if(modelparams->get<std::string>("PHYSICAL_MODEL","no_model") == "Dynamic_Smagorinsky")
      apply_dynamic_smagorinsky_ = true;

    if(modelparams->get<std::string>("PHYSICAL_MODEL","no_model") == "Multifractal_Subgrid_Scales")
      apply_box_filter_ = true;

    if (physicaltype_ == INPAR::FLUID::loma)
      loma_ = true;

    if(physicaltype_ == INPAR::FLUID::incompressible)
      incomp_=true;

    if(modelparams->get<std::string>("PHYSICAL_MODEL","no_model") == "Dynamic_Vreman")
      vreman_dynamic_ = true;
  }

  dynsmag_loma_on_ = (loma_ and apply_dynamic_smagorinsky_);

  if (apply_dynamic_smagorinsky_)
  {
    velocity_ = true;
    reynoldsstress_ = true;
    modeled_subgrid_stress_ = true;
  }
  if (apply_box_filter_)
  {
    velocity_ = true;
    reynoldsstress_ = true;
    finescale_velocity_ = true;
  }
  if (dynsmag_loma_on_)
  {
    densvelocity_ = true;
    densstrainrate_ = true;
    density_ = true;
    velocity_ = true;
    reynoldsstress_ = true;
    modeled_subgrid_stress_ = true;
  }
  if (loma_ and vreman_dynamic_)
    dserror("Dynamic Vreman model not implemented for loma!");

  return;
}


/*----------------------------------------------------------------------*
 | Destructor (public)                                                  |
 |                                                           gammi 09/08|
 *----------------------------------------------------------------------*/
FLD::Boxfilter::~Boxfilter()
{
  return;
}


/*----------------------------------------------------------------------*
 | add some scatra specific parameters                  rasthofer 08/12 |
 * ---------------------------------------------------------------------*/
void FLD::Boxfilter::AddScatra(
  Teuchos::RCP<DRT::Discretization>     scatradis
  )
{
  scatradiscret_ = scatradis;

  return;
}

void FLD::Boxfilter::InitializeVreman()
{
  strainrate_ = true;
  expression_ = true;
  alphaij_ = true;
  alpha2_ = true;

  return;
}

void FLD::Boxfilter::InitializeVremanScatra(
  Teuchos::RCP<DRT::Discretization>     scatradis
  )
{
  scatradiscret_ = scatradis;

  phi_=true;
  phi2_=true;
  phiexpression_=true;
  alphaijsc_=true;
  return;
}


/*---------------------------------------------------------------------*
 | Perform box filter operation                                        |
 *---------------------------------------------------------------------*/
void FLD::Boxfilter::ApplyFilter(
  const Teuchos::RCP<const Epetra_Vector> velocity,
  const Teuchos::RCP<const Epetra_Vector> scalar,
  const double                            thermpress,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{

  // perform filtering depending on the LES model
  ApplyBoxFilter(velocity,scalar,thermpress,dirichtoggle);

  return;
}

void FLD::Boxfilter::ApplyFilterScatra(
  const Teuchos::RCP<const Epetra_Vector>      scalar,
  const double                                 thermpress,
  const Teuchos::RCP<const Epetra_Vector>      dirichtoggle,
  const int                                    ndsvel
  )
{

  // perform filtering depending on the LES model
  ApplyBoxFilterScatra(scalar,thermpress,dirichtoggle,ndsvel);

  return;
}

/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                            rasthofer |
 *----------------------------------------------------------------------*/
void FLD::Boxfilter::ApplyBoxFilter(
  const Teuchos::RCP<const Epetra_Vector> velocity,
  const Teuchos::RCP<const Epetra_Vector> scalar,
  const double                            thermpress,
  const Teuchos::RCP<const Epetra_Vector> dirichtoggle
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("ApplyFilterForDynamicComputationOfCs");

  // LES turbulence modeling is only valid for 3 dimensions
  const int numdim =3;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList filterparams;
  // action for elements
  filterparams.set<int>("action",FLD::calc_fluid_box_filter);
  filterparams.set("thermpress",thermpress);

  // set state vector to pass distributed vector to the element
  discret_->ClearState();
  discret_->SetState("u and p (trial)",velocity);
  discret_->SetState("T (trial)",scalar);

  // dummies
  Epetra_SerialDenseMatrix emat1;
  Epetra_SerialDenseMatrix emat2;
  Epetra_SerialDenseVector evec1;
  Epetra_SerialDenseVector evec2;
  Epetra_SerialDenseVector evec3;

  // ---------------------------------------------------------------
  // get a vector layout from the discretization to construct
  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  // alloc an additional vector to store/add up the patch volume
  Teuchos::RCP<Epetra_Vector> patchvol     = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));

  // free mem and reallocate to zero out vecs
  if(velocity_)
    filtered_vel_                   = Teuchos::null;
  if(reynoldsstress_)
    filtered_reynoldsstress_        = Teuchos::null;
  if(modeled_subgrid_stress_)
    filtered_modeled_subgrid_stress_ = Teuchos::null;
  if(densvelocity_)
    filtered_dens_vel_ = Teuchos::null;
  if(density_)
    filtered_dens_ = Teuchos::null;
  if(densstrainrate_)
      filtered_dens_strainrate_ = Teuchos::null;
  if (finescale_velocity_)
    fs_vel_ = Teuchos::null;
  if (strainrate_)
    filtered_strainrate_=Teuchos::null;
  if (expression_)
    filtered_expression_ = Teuchos::null;
  if (alphaij_)
    filtered_alphaij_ = Teuchos::null;
  if (alpha2_)
    filtered_alpha2_ =Teuchos::null;

  if(velocity_)
    filtered_vel_                   = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));
  if(reynoldsstress_)
    filtered_reynoldsstress_        = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  if(modeled_subgrid_stress_)
    filtered_modeled_subgrid_stress_= Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  if(densvelocity_)
    filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));
  if(density_)
    filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  if(densstrainrate_)
    filtered_dens_strainrate_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  if (strainrate_)
    filtered_strainrate_= Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  if (expression_)
    filtered_expression_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  if (alphaij_)
    filtered_alphaij_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));
  if (alpha2_)
    filtered_alpha2_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));

  if (finescale_velocity_)
    fs_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim       ,true));

  // ---------------------------------------------------------------
  // do the integration of the (not normalized) box filter function
  // on the element

  // loop all elements on this proc (including ghosted ones)
  for (int nele=0;nele<discret_->NumMyColElements();++nele)
  {
    // get the element
    DRT::Element* ele = discret_->lColElement(nele);

    // provide vectors for filtered quantities
      Teuchos::RCP<std::vector<double> > vel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
      Teuchos::RCP<std::vector<std::vector<double> > > reynoldsstress_hat = Teuchos::rcp(new std::vector<std::vector<double> >);
      Teuchos::RCP<std::vector<std::vector<double> > > modeled_subgrid_stress = Teuchos::rcp(new std::vector<std::vector<double> >);
    // set to dimensions
    if(reynoldsstress_)
      (*reynoldsstress_hat).resize(numdim);
    if(modeled_subgrid_stress_)
      (*modeled_subgrid_stress).resize(numdim);
    for(int rr=0;rr<numdim;rr++)
    {
      if(reynoldsstress_)
        ((*reynoldsstress_hat)[rr]).resize(numdim);
      if(modeled_subgrid_stress_)
        ((*modeled_subgrid_stress)[rr]).resize(numdim);
    }
    // initialize with zeros
    for(int rr=0;rr<numdim;rr++)
    {
      for(int ss=0;ss<numdim;ss++)
      {
        if(reynoldsstress_)
          (*reynoldsstress_hat)[rr][ss] = 0.0;
        if(modeled_subgrid_stress_)
          (*modeled_subgrid_stress)[rr][ss] = 0.0;
      }
    }
    Teuchos::RCP<std::vector<double> > densvel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    // and set them in parameter list
    filterparams.set<Teuchos::RCP<std::vector<double> > >("densvel_hat",densvel_hat);

    filterparams.set<Teuchos::RCP<std::vector<double> > >("vel_hat",vel_hat);
    filterparams.set<Teuchos::RCP<std::vector<std::vector<double> > > >("reynoldsstress_hat",reynoldsstress_hat);
    filterparams.set<Teuchos::RCP<std::vector<std::vector<double> > > >("modeled_subgrid_stress",modeled_subgrid_stress);

    //Vreman_initialization
    Teuchos::RCP<std::vector<std::vector<double> > > strainrate_hat = Teuchos::rcp(new std::vector<std::vector<double> >);
    Teuchos::RCP<std::vector<std::vector<double> > > alphaij_hat = Teuchos::rcp(new std::vector<std::vector<double> >);
    if(strainrate_)
      (*strainrate_hat).resize(numdim);
    if(alphaij_)
      (*alphaij_hat).resize(numdim);
    for(int rr=0;rr<numdim;rr++)
    {
      if(strainrate_)
        ((*strainrate_hat)[rr]).resize(numdim);
      if(alphaij_)
        ((*alphaij_hat)[rr]).resize(numdim);
    }
    // initialize with zeros
    for(int rr=0;rr<numdim;rr++)
    {
      for(int ss=0;ss<numdim;ss++)
      {
        if(strainrate_)
          (*strainrate_hat)[rr][ss] = 0.0;
        if(alphaij_)
          (*alphaij_hat)[rr][ss] = 0.0;
      }
    }

    //if(strainrate_)
      filterparams.set<Teuchos::RCP<std::vector<std::vector<double> > > >("strainrate_hat",strainrate_hat);
    //if(alphaij_)
      filterparams.set<Teuchos::RCP<std::vector<std::vector<double> > > >("alphaij_hat",alphaij_hat);

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*discret_,lm,lmowner,lmstride);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(filterparams,
                            *discret_,
                            lm,
                            emat1,emat2,
                            evec1,evec2,evec2);
    if (err) dserror("Proc %d: Element %d returned err=%d",
                     discret_->Comm().MyPID(),ele->Id(),err);

    // get contribution to patch volume of this element. Add it up.
    double volume_contribution = filterparams.get<double>("volume_contribution");
    //if(density_)
    double dens_hat = filterparams.get<double>("dens_hat");
    //if(densstrainrate_)
    double dens_strainrate_hat = filterparams.get<double>("dens_strainrate_hat");

    double expression_hat = filterparams.get<double>("expression_hat");
    double alpha2_hat = filterparams.get<double>("alpha2_hat");

    // loop all nodes of this element, add values to the global vectors
    DRT::Node** elenodes=ele->Nodes();
    for(int nn=0;nn<ele->NumNode();++nn)
    {
      DRT::Node* node = (elenodes[nn]);

      // we are interested only in  row nodes
      if(node->Owner() == discret_->Comm().MyPID())
      {

        // now assemble the computed values into the global vector
        int    id = (node->Id());

        patchvol->SumIntoGlobalValues(1,&volume_contribution,&id);


        if(density_)
          filtered_dens_->SumIntoGlobalValues(1,&dens_hat,&id);
        if(densstrainrate_)
          filtered_dens_strainrate_->SumIntoGlobalValues(1,&dens_strainrate_hat,&id);
        if(expression_)
          filtered_expression_->SumIntoGlobalValues(1,&expression_hat,&id);
        if(alpha2_)
          filtered_alpha2_->SumIntoGlobalValues(1,&alpha2_hat,&id);

        for (int idim =0;idim<numdim;++idim)
        {
          if(velocity_)
          {
            double val = (*vel_hat)[idim];
            ((*filtered_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
          }
          if (densvelocity_)
          {
            double val = (*densvel_hat)[idim];
            ((*filtered_dens_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
          }


          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;
            if (reynoldsstress_)
            {
              double val = (*reynoldsstress_hat)[idim][jdim];
              ((*filtered_reynoldsstress_ )       (ij))->SumIntoGlobalValues(1,&val,&id);
            }
            if (modeled_subgrid_stress_)
            {
              double val = (*modeled_subgrid_stress)[idim][jdim];
              ((*filtered_modeled_subgrid_stress_)(ij))->SumIntoGlobalValues(1,&val,&id);
            }
            if (strainrate_)
            {
              double val = (*strainrate_hat)[idim][jdim];
              ((*filtered_strainrate_)(ij))->SumIntoGlobalValues(1,&val,&id);
            }
            if (alphaij_)
            {
              double val = (*alphaij_hat)[idim][jdim];
              ((*filtered_alphaij_)(ij))->SumIntoGlobalValues(1,&val,&id);
            }


          }
        }

      }
    }
  } // end elementloop

  // ---------------------------------------------------------------
  // send add values from masters and slaves
  {
    std::map<int, std::vector<int> >::iterator masternode;

    double val;
    //if(velocity_)
    std::vector<double> vel_val(3);
    std::vector<std::vector<double> > reystress_val;
    if (reynoldsstress_)
    {
      reystress_val.resize(3);
      for(int rr=0;rr<3;rr++)
        (reystress_val[rr]).resize(3);
    }
    std::vector<std::vector<double> > modeled_subgrid_stress_val;
    if (modeled_subgrid_stress_)
    {
      modeled_subgrid_stress_val.resize(3);
      for(int rr=0;rr<3;rr++)
        (modeled_subgrid_stress_val[rr]).resize(3);
    }
    std::vector<std::vector<double> > strainrate_val;
    if (strainrate_)
    {
      strainrate_val.resize(3);
      for(int rr=0;rr<3;rr++)
        (strainrate_val[rr]).resize(3);
    }
    std::vector<std::vector<double> > alphaij_val;
    if (alphaij_)
    {
      alphaij_val.resize(3);
      for(int rr=0;rr<3;rr++)
        (alphaij_val[rr]).resize(3);
    }
    // loma specific quantities
    std::vector<double> dens_vel_val(3);
    double dens_val;
    double dens_strainrate_val;
    double expression_val;
    double alpha2_val;

    Teuchos::RCP<std::map<int,std::vector<int> > >  pbcmapmastertoslave = discret_->GetAllPBCCoupledColNodes();
    // loop all master nodes on this proc
    for(masternode =pbcmapmastertoslave->begin();
        masternode!=pbcmapmastertoslave->end();
        ++masternode)
    {
      // loop only owned nodes
      if ((discret_->gNode(masternode->first))->Owner() != discret_->Comm().MyPID())
        continue;

      // add all slave values to master value
      std::vector<int>::iterator slavenode;

      int lid = noderowmap->LID(masternode->first);
      if (lid < 0) dserror("nodelid < 0 ?");

      val = (*patchvol)[lid];

      if(density_)
        dens_val = (*filtered_dens_)[lid];
      if(densstrainrate_)
        dens_strainrate_val = (*filtered_dens_strainrate_)[lid];
      if(expression_)
        expression_val =(*filtered_expression_)[lid];
      if(alpha2_)
        alpha2_val =(*filtered_alpha2_)[lid];

      for (int idim =0;idim<numdim;++idim)
      {
        if(velocity_)
          vel_val[idim]=((*((*filtered_vel_)(idim)))[lid]);

        if (densvelocity_)
          dens_vel_val[idim] = ((*((*filtered_dens_vel_)(idim)))[lid]);

        for (int jdim =0;jdim<numdim;++jdim)
        {
          const int ij = numdim*idim+jdim;
          if (reynoldsstress_)
            reystress_val             [idim][jdim] = (*((*filtered_reynoldsstress_         ) (ij)))[lid];
          if (modeled_subgrid_stress_)
            modeled_subgrid_stress_val[idim][jdim] = (*((*filtered_modeled_subgrid_stress_ ) (ij)))[lid];
          if (strainrate_)
            strainrate_val            [idim][jdim] = (*((*filtered_strainrate_             ) (ij)))[lid];
          if (alphaij_)
            alphaij_val               [idim][jdim] = (*((*filtered_alphaij_                ) (ij)))[lid];
        }
      }

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        lid = noderowmap->LID(*slavenode);
        val += (*patchvol)[lid];

        if(density_)
          dens_val += (*filtered_dens_)[lid];
        if(densstrainrate_)
          dens_strainrate_val += (*filtered_dens_strainrate_)[lid];
        if(expression_)
          expression_val += (*filtered_expression_)[lid];
        if(alpha2_)
          alpha2_val += (*filtered_alpha2_)[lid];

        for (int idim =0;idim<numdim;++idim)
        {
          if(velocity_)
            vel_val[idim] += ((*((*filtered_vel_)(idim)))[lid]);

          if (densvelocity_)
            dens_vel_val[idim] += ((*((*filtered_dens_vel_)(idim)))[lid]);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;
            if (reynoldsstress_)
              reystress_val             [idim][jdim] += (*((*filtered_reynoldsstress_         ) (ij)))[lid];
            if (modeled_subgrid_stress_)
              modeled_subgrid_stress_val[idim][jdim] += (*((*filtered_modeled_subgrid_stress_ ) (ij)))[lid];
            if (strainrate_)
              strainrate_val            [idim][jdim] += (*((*filtered_strainrate_             ) (ij)))[lid];
            if (alphaij_)
              alphaij_val               [idim][jdim] += (*((*filtered_alphaij_                ) (ij)))[lid];
          } // end loop jdim
        } // end loop idim
      }  // end loop slaves

      // replace value by sum
      lid = noderowmap->LID(masternode->first);
      int error = patchvol->ReplaceMyValues(1,&val,&lid);
      if (error != 0) dserror("dof not on proc");

      {
        int err = 0;
        if(density_)
          err += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
        if(densstrainrate_)
          err += filtered_dens_strainrate_->ReplaceMyValues(1,&dens_strainrate_val,&lid);
        if(expression_)
          err += filtered_expression_->ReplaceMyValues(1,&expression_val,&lid);
        if(alpha2_)
          err += filtered_alpha2_->ReplaceMyValues(1,&alpha2_val,&lid);
        if (err != 0) dserror("dof not on proc");
      }

      for (int idim =0;idim<numdim;++idim)
      {
        int err = 0;
        if(velocity_)
          err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);

        if (densvelocity_)
          err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);

        for (int jdim =0;jdim<numdim;++jdim)
        {
          const int ij = numdim*idim+jdim;
          if (reynoldsstress_)
            err += ((*filtered_reynoldsstress_        )(ij))->ReplaceMyValues(1,&(reystress_val             [idim][jdim]),&lid);
          if (modeled_subgrid_stress_)
            err += ((*filtered_modeled_subgrid_stress_)(ij))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
          if (strainrate_)
            err += ((*filtered_strainrate_            )(ij))->ReplaceMyValues(1,&(strainrate_val            [idim][jdim]),&lid);
          if (alphaij_)
            err += ((*filtered_alphaij_               )(ij))->ReplaceMyValues(1,&(alphaij_val               [idim][jdim]),&lid);
        } // end loop jdim
        if (err != 0) dserror("dof not on proc");
      } // end loop idim

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        int err = 0;
        lid = noderowmap->LID(*slavenode);
        err += patchvol->ReplaceMyValues(1,&val,&lid);

        {
          int err = 0;
          if(density_)
            err += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
          if(densstrainrate_)
            err += filtered_dens_strainrate_->ReplaceMyValues(1,&dens_strainrate_val,&lid);
          if(expression_)
            err += filtered_expression_->ReplaceMyValues(1,&expression_val,&lid);
          if(alpha2_)
            err += filtered_alpha2_->ReplaceMyValues(1,&alpha2_val,&lid);
          if (err != 0) dserror("dof not on proc");
        }

        for (int idim =0;idim<numdim;++idim)
        {
          if(velocity_)
            err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);

          if (densvelocity_)
            err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;
            if (reynoldsstress_)
              err += ((*filtered_reynoldsstress_        )(ij))->ReplaceMyValues(1,&(reystress_val             [idim][jdim]),&lid);
            if (modeled_subgrid_stress_)
              err += ((*filtered_modeled_subgrid_stress_)(ij))->ReplaceMyValues(1,&(modeled_subgrid_stress_val[idim][jdim]),&lid);
            if (strainrate_)
              err += ((*filtered_strainrate_            )(ij))->ReplaceMyValues(1,&(strainrate_val            [idim][jdim]),&lid);
            if (alphaij_)
              err += ((*filtered_alphaij_               )(ij))->ReplaceMyValues(1,&(alphaij_val               [idim][jdim]),&lid);
          } // end loop jdim
        } // end loop idim
        if (err != 0) dserror("dof not on proc");
      } // end loop slaves
    } // end loop masters
  }

  // ---------------------------------------------------------------
  // replace values at dirichlet nodes

  {
    // get a rowmap for the dofs
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode       = discret_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = discret_->Dof(lnode);

      // check whether the node is on a wall, i.e. all velocity dofs
      // are Dirichlet constrained
      int is_dirichlet_node = 0;
      int is_no_slip_node = 0;
      for (int index=0;index<numdim;++index)
      {
        int gid = nodedofset[index];
        int lid = dofrowmap->LID(gid);

        if ((*dirichtoggle)[lid]==1) //this is a dirichlet node
        {
          is_dirichlet_node++;
          double vel_i = (*velocity)[lid];
          if (abs(vel_i) < 1e-12)
          {
            is_no_slip_node++;
          }
        }
      }

      // this node is on a dirichlet boundary
      if (is_dirichlet_node == numdim)
      {
        int err = 0;

        // determine volume
        double thisvol = (*patchvol)[lnodeid];

        // determine density
        double dens = 1.0;
        if (physicaltype_ == INPAR::FLUID::incompressible) // this is important to have here,
        {                                                                                 //  since, for the pure box filter application,
           // get fluid viscosity from material definition                                //  we do not want to multiply the reynolds stress by density
          int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
          if (id==-1)
            dserror("Could not find Newtonian fluid material");
          else
          {
            const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
            const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
            // we need the kinematic viscosity here
            dens = actmat->density_;
          }
        }
        if(density_)
          dens = (*filtered_dens_)[lnodeid]/thisvol;


        // set density (only required for loma)
        // set value to mean value
        // we already divide by the corresponding volume of all contributing elements,
        // since we set the volume to 1.0 in the next step in order not to modify the dirichlet values
        if(density_)
          err += filtered_dens_->ReplaceMyValues(1,&dens,&lnodeid);

        // this node is on a wall
        if (is_no_slip_node == numdim)
        {
          // Peter style
          double val = 0.0;
          if(densstrainrate_)
            err += filtered_dens_strainrate_->ReplaceMyValues(1,&val,&lnodeid);
          if(expression_)
            err += filtered_expression_->ReplaceMyValues(1,&val,&lnodeid);
          if(alpha2_)
            err += filtered_alpha2_->ReplaceMyValues(1,&val,&lnodeid);

        }
        else
        {
          if(densstrainrate_)
          {
            double val = (*filtered_dens_strainrate_)[lnodeid]/thisvol;
            err += filtered_dens_strainrate_->ReplaceMyValues(1,&val,&lnodeid);
          }
          if(expression_)
          {
            double val = (*filtered_expression_)[lnodeid]/thisvol;
            err += filtered_expression_->ReplaceMyValues(1,&val,&lnodeid);
          }
          if(alpha2_)
          {
            double val = (*filtered_alpha2_)[lnodeid]/thisvol;
            err += filtered_alpha2_->ReplaceMyValues(1,&val,&lnodeid);
          }
        }



        for (int idim =0;idim<numdim;++idim)
        {
          int gid_i = nodedofset[idim];
          int lid_i = dofrowmap->LID(gid_i);
          double valvel_i = (*velocity)[lid_i];
          if (velocity_)
          {
            err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&valvel_i,&lnodeid);
          }
          //dens*reynoldsstress not in parameter list until now?
          if (densvelocity_)//=loma
          {
            // note: for incompressible flow, this vector is rebuild in calculation of Lij and Mij
            double valdensvel_i = dens*valvel_i;
            err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&valdensvel_i,&lnodeid);
          }

          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;

            if (reynoldsstress_){
              int gid_j = nodedofset[jdim];
              int lid_j = dofrowmap->LID(gid_j);

              double valvel_j = (*velocity)[lid_j];
              double valvel_ij= dens * valvel_i * valvel_j;
              // remember: density = 1.0 for pure box filter application
              err += ((*filtered_reynoldsstress_         ) (ij))->ReplaceMyValues(1,&valvel_ij,&lnodeid);
            }

            if (is_no_slip_node == numdim)
            {
              // set value to zero (original Peter style)
              double val = 0.0;
              if (modeled_subgrid_stress_)
                err += ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              // remark: setting the modeled stresses equal to zero improves the estimated friction Reynolds number!
              if (strainrate_)
                err += ((*filtered_strainrate_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              if (alphaij_)
                err += ((*filtered_alphaij_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
            }
          else
            {
              // set value to mean value
              // we already divide by the corresponding volume of all contributing elements,
              // since we set the volume to 1.0 in the next step in order not to modify the dirichlet values
              if (modeled_subgrid_stress_)
              {
                double val = ((*((*filtered_modeled_subgrid_stress_ ) (ij)))[lnodeid])/thisvol;
                err += ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              }
              if (strainrate_)
              {
                double val = ((*((*filtered_strainrate_ ) (ij)))[lnodeid])/thisvol;
                err += ((*filtered_strainrate_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              }
              if (alphaij_)
              {
                double val = ((*((*filtered_alphaij_ ) (ij)))[lnodeid])/thisvol;
                err += ((*filtered_alphaij_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
              }
            }
          } // end loop jdim
        } // end loop idim

        double volval = 1.0;
        err += patchvol->ReplaceMyValues(1,&volval,&lnodeid);
        if (err!=0) dserror("dof/node not on proc");
      }// is dirichlet node
    } // end loop all nodes
  }


  // ---------------------------------------------------------------
  // scale vectors by element patch sizes --- this corresponds to
  // the normalization of the box filter function

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<discret_->NumMyRowNodes();++lnodeid)
  {
    double thisvol = (*patchvol)[lnodeid];

    int err = 0;
    double val = 0.0;
    if (density_)
    {
      val = (*filtered_dens_)[lnodeid]/thisvol;
      err += filtered_dens_->ReplaceMyValues(1,&val,&lnodeid);
    }
    if (densstrainrate_)
    {
      val = (*filtered_dens_strainrate_)[lnodeid]/thisvol;
      err += filtered_dens_strainrate_->ReplaceMyValues(1,&val,&lnodeid);
    }
    if (expression_)
    {
      val = (*filtered_expression_)[lnodeid]/thisvol;
      err += filtered_expression_->ReplaceMyValues(1,&val,&lnodeid);
    }
    if (alpha2_)
    {
      val = (*filtered_alpha2_)[lnodeid]/thisvol;
      err += filtered_alpha2_->ReplaceMyValues(1,&val,&lnodeid);
    }

    for (int idim =0;idim<3;++idim)
    {
      if (velocity_)
      {
        val = ((*((*filtered_vel_)(idim)))[lnodeid])/thisvol;
        err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      }

      if (densvelocity_)
      {
        val = ((*((*filtered_dens_vel_)(idim)))[lnodeid])/thisvol;
        err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      }

      for (int jdim =0;jdim<3;++jdim)
      {
        const int ij = numdim*idim+jdim;

        if (reynoldsstress_)
        {
          val = ((*((*filtered_reynoldsstress_ ) (ij)))[lnodeid])/thisvol;
          err += ((*filtered_reynoldsstress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
        }
        if (modeled_subgrid_stress_)
        {
          val = ((*((*filtered_modeled_subgrid_stress_ ) (ij)))[lnodeid])/thisvol;
          err += ((*filtered_modeled_subgrid_stress_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
        }
        if (strainrate_)
        {
          val = ((*((*filtered_strainrate_ ) (ij)))[lnodeid])/thisvol;
          err += ((*filtered_strainrate_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
        }
        if (alphaij_)
        {
          val = ((*((*filtered_alphaij_ ) (ij)))[lnodeid])/thisvol;
          err += ((*filtered_alphaij_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
        }
      } // end loop jdim
      if (err!=0) dserror("dof not on proc");
    } // end loop idim
  } // end loop nodes

  // clean up
  discret_->ClearState();

  //calculate fine scale velocities
  if (finescale_velocity_)
  {
    // fine scale veocity requires filtered velocity
    if (not velocity_)
      dserror("filtered velocity is required in the box filter to calculate the fine scale velocity");
    // loop all elements on this proc
    for (int nid=0;nid<discret_->NumMyRowNodes();++nid)
    {
      // get the node
      DRT::Node* node = discret_->lRowNode(nid);
      // get global ids of all dofs of the node
      std::vector<int> dofs= discret_->Dof(node);
      //we only loop over all velocity dofs
      for(int d=0;d<discret_->NumDof(node)-1;++d)
      {
        // get global id of the dof
        int gid = dofs[d];
        // get local dof id corresponding to the global id
        int lid = discret_->DofRowMap()->LID(gid);
        // filtered velocity and all scale velocity
        double filteredvel = (*((*filtered_vel_)(d)))[nid];
        double vel = (*velocity)[lid];
        // calculate fine scale velocity
        double val = vel - filteredvel;
        // calculate fine scale velocity
        int err = ((*fs_vel_)(d))->ReplaceMyValues(1,&val,&nid);
        if (err!=0) dserror("dof not on proc");
      }
    }
  }

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = discret_->NodeColMap();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  if (velocity_)
    col_filtered_vel_                    = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  if (reynoldsstress_)
    col_filtered_reynoldsstress_         = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (modeled_subgrid_stress_)
    col_filtered_modeled_subgrid_stress_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (finescale_velocity_)
    col_fs_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  if (densvelocity_)
    col_filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  if (density_)
    col_filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  if(densstrainrate_)
    col_filtered_dens_strainrate_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  if (strainrate_)
    col_filtered_strainrate_         = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (alphaij_)
    col_filtered_alphaij_         = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,9,true));
  if (expression_)
      col_filtered_expression_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  if (alpha2_)
      col_filtered_alpha2_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));

  // export filtered vectors in rowmap to columnmap format
  if (velocity_)
    LINALG::Export(*filtered_vel_                   ,*col_filtered_vel_                   );
  if (reynoldsstress_)
    LINALG::Export(*filtered_reynoldsstress_        ,*col_filtered_reynoldsstress_        );
  if (modeled_subgrid_stress_)
    LINALG::Export(*filtered_modeled_subgrid_stress_,*col_filtered_modeled_subgrid_stress_);
  if (finescale_velocity_)
    LINALG::Export(*fs_vel_                   ,*col_fs_vel_                   );
  if (densvelocity_)
    LINALG::Export(*filtered_dens_vel_,*col_filtered_dens_vel_);
  if (density_)
    LINALG::Export(*filtered_dens_,*col_filtered_dens_);
  if(densstrainrate_)
    LINALG::Export(*filtered_dens_strainrate_,*col_filtered_dens_strainrate_);
  if (strainrate_)
      LINALG::Export(*filtered_strainrate_        ,*col_filtered_strainrate_        );
  if (alphaij_)
      LINALG::Export(*filtered_alphaij_        ,*col_filtered_alphaij_        );
  if (expression_)
    LINALG::Export(*filtered_expression_,*col_filtered_expression_);
  if (alpha2_)
    LINALG::Export(*filtered_alpha2_,*col_filtered_alpha2_);
return;
}




/*----------------------------------------------------------------------*
 | perform box filtering                                      (private) |
 |                                                      rasthofer 08/12 |
 *----------------------------------------------------------------------*/
void FLD::Boxfilter::ApplyBoxFilterScatra(
  const Teuchos::RCP<const Epetra_Vector>      scalar,
  const double                                 thermpress,
  const Teuchos::RCP<const Epetra_Vector>      dirichtoggle,
  const int                                    ndsvel
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("ApplyFilterForDynamicComputationOfPrt");
  if (apply_box_filter_ == true) dserror("not yet considered");
  // LES turbulence modeling is only valid for 3 dimensions
  const int numdim =3;

  // generate a parameterlist for communication and control
  Teuchos::ParameterList filterparams;
  // action for elements
  filterparams.set<int>("action",SCATRA::calc_scatra_box_filter);

  // add number of dofset associated with velocity related dofs to parameter list
  filterparams.set<int>("ndsvel",ndsvel);

  filterparams.set("thermpress",thermpress);

  // set state vector to pass distributed vector to the element
  scatradiscret_->ClearState();
  scatradiscret_->SetState("scalar",scalar);

  // dummies
  Epetra_SerialDenseMatrix emat1;
  Epetra_SerialDenseMatrix emat2;
  Epetra_SerialDenseVector evec1;
  Epetra_SerialDenseVector evec2;
  Epetra_SerialDenseVector evec3;

  // ---------------------------------------------------------------
  // get a vector layout from the discretization to construct
  const Epetra_Map* noderowmap = scatradiscret_->NodeRowMap();

  // alloc an additional vector to store/add up the patch volume
  Teuchos::RCP<Epetra_Vector> patchvol     = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));

  // free mem and reallocate to zero out vecs
  filtered_dens_vel_temp_ = Teuchos::null;
  filtered_dens_rateofstrain_temp_ = Teuchos::null;
  filtered_vel_ = Teuchos::null;
  filtered_dens_vel_ = Teuchos::null;
  filtered_temp_ = Teuchos::null;
  filtered_dens_temp_ = Teuchos::null;
  filtered_dens_ = Teuchos::null;
  if (phi_)
    filtered_phi_ = Teuchos::null;
  if (phi2_)
    filtered_phi2_ = Teuchos::null;
  if (phiexpression_)
    filtered_phiexpression_ = Teuchos::null;
  if (alphaijsc_)
    filtered_alphaijsc_ =Teuchos::null;

  filtered_dens_vel_temp_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_dens_rateofstrain_temp_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  filtered_temp_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  filtered_dens_temp_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  if (phi_)
    filtered_phi_=Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim,true));
  if (phi2_)
    filtered_phi2_=Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  if (phiexpression_)
    filtered_phiexpression_=Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  if (alphaijsc_)
    filtered_alphaijsc_ =Teuchos::rcp(new Epetra_MultiVector(*noderowmap,numdim*numdim,true));

  // ---------------------------------------------------------------
  // do the integration of the (not normalized) box filter function
  // on the element

  // loop all elements on this proc (including ghosted ones)
  for (int nele=0;nele<scatradiscret_->NumMyColElements();++nele)
  {
    // get the element
    DRT::Element* ele = scatradiscret_->lColElement(nele);

    // provide vectors for filtered quantities //declaration necessary even if not used
    Teuchos::RCP<std::vector<double> > vel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    Teuchos::RCP<std::vector<double> > densvel_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    Teuchos::RCP<std::vector<double> > densveltemp_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    Teuchos::RCP<std::vector<double> > densstraintemp_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    Teuchos::RCP<std::vector<double> > phi_hat = Teuchos::rcp(new std::vector<double> ((numdim),0.0));
    Teuchos::RCP<std::vector<std::vector<double> > > alphaijsc_hat = Teuchos::rcp(new std::vector<std::vector<double> >);
    if(alphaijsc_)
    {
      (*alphaijsc_hat).resize(numdim);
      for(int rr=0;rr<numdim;rr++)
          ((*alphaijsc_hat)[rr]).resize(numdim);
      // initialize with zeros
      for(int rr=0;rr<numdim;rr++)
      {
        for(int ss=0;ss<numdim;ss++)
            (*alphaijsc_hat)[rr][ss] = 0.0;
      }
    }
    // and set them in parameter list
    filterparams.set<Teuchos::RCP<std::vector<double> > >("vel_hat",vel_hat);
    filterparams.set<Teuchos::RCP<std::vector<double> > >("densvel_hat",densvel_hat);
    filterparams.set<Teuchos::RCP<std::vector<double> > >("densveltemp_hat",densveltemp_hat);
    filterparams.set<Teuchos::RCP<std::vector<double> > >("densstraintemp_hat",densstraintemp_hat);
    filterparams.set<Teuchos::RCP<std::vector<double> > >("phi_hat",phi_hat);
    filterparams.set<Teuchos::RCP<std::vector<std::vector<double> > > >("alphaijsc_hat",alphaijsc_hat);

    // initialize variables for filtered scalar quantities
    double dens_hat = 0.0;
    double temp_hat = 0.0;
    double dens_temp_hat = 0.0;
    double phi2_hat=0.0;
    double phiexpression_hat=0.0;

    // initialize volume contribution
    double volume_contribution = 0.0;

    // get element location vector, dirichlet flags and ownerships
    DRT::Element::LocationArray la(scatradiscret_->NumDofSets());
    ele->LocationVector(*scatradiscret_,la,false);

    // call the element evaluate method to integrate functions
    // against heaviside function element
    int err = ele->Evaluate(filterparams,
                            *scatradiscret_,
                            la,
                            emat1,emat2,
                            evec1,evec2,evec2);
    if (err) dserror("Proc %d: Element %d returned err=%d",
                     scatradiscret_->Comm().MyPID(),ele->Id(),err);

    // get contribution to patch volume of this element. Add it up.
    //double volume_contribution = filterparams.get<double>("volume_contribution");
    volume_contribution = filterparams.get<double>("volume_contribution");

    // filtered scalar quantities
    dens_hat = filterparams.get<double>("dens_hat");
    temp_hat = filterparams.get<double>("temp_hat");
    dens_temp_hat = filterparams.get<double>("dens_temp_hat");
    if (phi2_)
      phi2_hat=filterparams.get<double>("phi2_hat");
    if (phiexpression_)
      phiexpression_hat=filterparams.get<double>("phiexpression_hat");
    // loop all nodes of this element, add values to the global vectors
    DRT::Node** elenodes=ele->Nodes();
    for(int nn=0;nn<ele->NumNode();++nn)
    {
      DRT::Node* node = (elenodes[nn]);

      // we are interested only in  row nodes
      if(node->Owner() == scatradiscret_->Comm().MyPID())
      {

        // now assemble the computed values into the global vector
        int    id = (node->Id());

        patchvol->SumIntoGlobalValues(1,&volume_contribution,&id);
        filtered_dens_->SumIntoGlobalValues(1,&dens_hat,&id);
        filtered_temp_->SumIntoGlobalValues(1,&temp_hat,&id);
        filtered_dens_temp_->SumIntoGlobalValues(1,&dens_temp_hat,&id);
        if(phi2_)
          filtered_phi2_->SumIntoGlobalValues(1,&phi2_hat,&id);
        if (phiexpression_)
          filtered_phiexpression_->SumIntoGlobalValues(1,&phiexpression_hat,&id);
        for (int idim =0;idim<numdim;++idim)
        {
           double val = (*vel_hat)[idim];
          ((*filtered_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
          val = (*densveltemp_hat)[idim];
          ((*filtered_dens_vel_temp_)(idim))->SumIntoGlobalValues(1,&val,&id);
          val = (*densstraintemp_hat)[idim];
          ((*filtered_dens_rateofstrain_temp_)(idim))->SumIntoGlobalValues(1,&val,&id);
          val = (*densvel_hat)[idim];
          ((*filtered_dens_vel_)(idim))->SumIntoGlobalValues(1,&val,&id);
          if(phi_)
          {
            val = (*phi_hat)[idim];
            ((*filtered_phi_)(idim))->SumIntoGlobalValues(1,&val,&id);
          }
          if (alphaijsc_)
          {
            for (int jdim =0;jdim<numdim;++jdim)
            {
              const int ij = numdim*idim+jdim;
              double val = (*alphaijsc_hat)[idim][jdim];
              ((*filtered_alphaijsc_)(ij))->SumIntoGlobalValues(1,&val,&id);
            }
          }
        }
      }
    }
  } // end elementloop

  // ---------------------------------------------------------------
  // send add values from masters and slaves
  {
    std::map<int, std::vector<int> >::iterator masternode;
    double val = 0.0;
    std::vector<double> vel_val(3);
    std::vector<double> dens_vel_val(3);
    std::vector<double> dens_vel_temp_val(3);
    std::vector<double> dens_strain_temp_val(3);
    std::vector<double> phi_val(3);
    std::vector<std::vector<double> > alphaijsc_val;
    if (alphaijsc_)
    {
      alphaijsc_val.resize(3);
      for(int rr=0;rr<3;rr++)
        (alphaijsc_val[rr]).resize(3);
    }
    double temp_val = 0.0;
    double dens_val = 0.0;
    double dens_temp_val = 0.0;
    double phi2_val=0.0;
    double phiexpression_val=0.0;

    // loop all master nodes on this proc
    Teuchos::RCP<std::map<int,std::vector<int> > >  pbcmapmastertoslave = scatradiscret_->GetAllPBCCoupledColNodes();
    // loop all master nodes on this proc
    for(masternode =pbcmapmastertoslave->begin();
        masternode!=pbcmapmastertoslave->end();
        ++masternode)
    {
      // loop only owned nodes
      if ((scatradiscret_->gNode(masternode->first))->Owner() != scatradiscret_->Comm().MyPID())
        continue;

      // add all slave values to master value
      std::vector<int>::iterator slavenode;

      int lid = noderowmap->LID(masternode->first);
      if (lid < 0) dserror("nodelid < 0 ?");

      val = (*patchvol)[lid];

      dens_val = (*filtered_dens_)[lid];
      dens_temp_val = (*filtered_dens_temp_)[lid];
      temp_val = (*filtered_temp_)[lid];
      if (phi2_)
        phi2_val = (*filtered_phi2_)[lid];
      if (phiexpression_)
        phiexpression_val=(*filtered_phiexpression_)[lid];

      for (int idim =0;idim<numdim;++idim)
      {
        vel_val[idim] = ((*((*filtered_vel_)(idim)))[lid]);
        dens_vel_val[idim] = ((*((*filtered_dens_vel_)(idim)))[lid]);
        dens_vel_temp_val[idim] = ((*((*filtered_dens_vel_temp_)(idim)))[lid]);
        dens_strain_temp_val[idim] = ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lid]);
        if (phi_)
          phi_val[idim] = ((*((*filtered_phi_)(idim)))[lid]);
        if (alphaijsc_)
        {
          for (int jdim =0;jdim<numdim;++jdim)
          {
            const int ij = numdim*idim+jdim;
            alphaijsc_val               [idim][jdim] = (*((*filtered_alphaijsc_                ) (ij)))[lid];
          }
        }
      }

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        lid = noderowmap->LID(*slavenode);
        val += (*patchvol)[lid];

        dens_val += (*filtered_dens_)[lid];
        dens_temp_val += (*filtered_dens_temp_)[lid];
        temp_val += (*filtered_temp_)[lid];
        if (phi2_)
          phi2_val += (*filtered_phi2_)[lid];
        if (phiexpression_)
          phiexpression_val += (*filtered_phiexpression_)[lid];

        for (int idim =0;idim<numdim;++idim)
        {
          vel_val[idim] +=((*((*filtered_vel_)(idim)))[lid]);
          dens_vel_val[idim] += ((*((*filtered_dens_vel_)(idim)))[lid]);
          dens_vel_temp_val[idim] += ((*((*filtered_dens_vel_temp_)(idim)))[lid]);
          dens_strain_temp_val[idim] += ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lid]);
          if (phi_)
            phi_val[idim] += ((*((*filtered_phi_)(idim)))[lid]);
          if (alphaijsc_)
          {
            for (int jdim =0;jdim<numdim;++jdim)
            {
              const int ij = numdim*idim+jdim;
              alphaijsc_val [idim][jdim] += (*((*filtered_alphaijsc_ ) (ij)))[lid];
            }
          }

        }
      }  // end loop slaves

      // replace value by sum
      lid = noderowmap->LID(masternode->first);
      int error = patchvol->ReplaceMyValues(1,&val,&lid);
      if (error != 0) dserror("dof not on proc");

      int e = 0;
      e += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
      e += filtered_dens_temp_->ReplaceMyValues(1,&dens_temp_val,&lid);
      e += filtered_temp_->ReplaceMyValues(1,&temp_val,&lid);
      if (phi2_)
        e += filtered_phi2_->ReplaceMyValues(1,&phi2_val,&lid);
      if (phiexpression_)
        e += filtered_phiexpression_->ReplaceMyValues(1,&phiexpression_val,&lid);
      if (e != 0) dserror("dof not on proc");

      for (int idim =0;idim<numdim;++idim)
      {
        int err = 0;
        err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);
        err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);
        err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&(dens_vel_temp_val[idim]),&lid);
        err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&(dens_strain_temp_val[idim]),&lid);
        if (phi_)
          err += ((*filtered_phi_)(idim))->ReplaceMyValues(1,&(phi_val[idim]),&lid);
        if (alphaijsc_)
        {
          for (int jdim =0;jdim<numdim;++jdim)
            {
              const int ij = numdim*idim+jdim;
                err += ((*filtered_alphaijsc_   )(ij))->ReplaceMyValues(1,&(alphaijsc_val  [idim][jdim]),&lid);
            } // end loop jdim
        }
        if (err != 0) dserror("dof not on proc");
      }

      // loop all this masters slaves
      for(slavenode=(masternode->second).begin();slavenode!=(masternode->second).end();++slavenode)
      {
        int err = 0;
        lid = noderowmap->LID(*slavenode);
        err += patchvol->ReplaceMyValues(1,&val,&lid);

        err += filtered_dens_->ReplaceMyValues(1,&dens_val,&lid);
        err += filtered_dens_temp_->ReplaceMyValues(1,&dens_temp_val,&lid);
        err += filtered_temp_->ReplaceMyValues(1,&temp_val,&lid);
        if (phi2_)
          err += filtered_phi2_->ReplaceMyValues(1,&phi2_val,&lid);
        if(phiexpression_)
          err += filtered_phiexpression_->ReplaceMyValues(1,&phiexpression_val,&lid);

        for (int idim =0;idim<numdim;++idim)
        {
          err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&(vel_val[idim]),&lid);
          err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&(dens_vel_val[idim]),&lid);
          err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&(dens_vel_temp_val[idim]),&lid);
          err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&(dens_strain_temp_val[idim]),&lid);
          if (phi_)
            err += ((*filtered_phi_)(idim))->ReplaceMyValues(1,&(phi_val[idim]),&lid);
          if (alphaijsc_)
          {
            for (int jdim =0;jdim<numdim;++jdim)
              {
                const int ij = numdim*idim+jdim;
                  err += ((*filtered_alphaijsc_   )(ij))->ReplaceMyValues(1,&(alphaijsc_val  [idim][jdim]),&lid);
              } // end loop jdim
          }
        }

        if (err != 0) dserror("dof not on proc");
      } // end loop slaves
    } // end loop masters
  }

  // ---------------------------------------------------------------
  // extract convective velocity from scatra discretization
  Teuchos::RCP<const Epetra_Vector> convel = scatradiscret_->GetState(ndsvel,"convective velocity field");
  if (convel == Teuchos::null)
    dserror("Cannot extract convective velocity field");

  // replace values at dirichlet nodes
  {
    // get a rowmap for the dofs
    const Epetra_Map* dofrowmap = scatradiscret_->DofRowMap();

    // as we want to identify nodes at walls,
    // we have to be sure that fluid and scatra are still matching
    if (not scatradiscret_->NodeRowMap()->SameAs(*(discret_->NodeRowMap())))
      dserror("Fluid and ScaTra noderowmaps are NOT identical.");

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<scatradiscret_->NumMyRowNodes();++lnodeid)
    {
      // get the processor local node
      DRT::Node*  lnode = scatradiscret_->lRowNode(lnodeid);
      std::vector<int> nodedofs = scatradiscret_->Dof(ndsvel,lnode);
      // get the corresponding porcessor local fluid node
      DRT::Node*  fluidlnode = discret_->lRowNode(lnodeid);

      // do we have a dirichlet boundary conditions in the fluid
      std::vector<DRT::Condition*> dbccond;
      fluidlnode->GetCondition("Dirichlet",dbccond);

      // yes, we have a dirichlet boundary condition
      if (dbccond.size()>0)
      {
#if DEBUG
        if ((lnode->X()[0]!=fluidlnode->X()[0]) or
            (lnode->X()[1]!=fluidlnode->X()[1]) or
            (lnode->X()[2]!=fluidlnode->X()[2]))
          dserror("Nodes do not match.");
#endif
        // we only want to modify nodes at the wall, as the model should vanish there
        // check, whether we have a no-slip node
        int no_slip_node = 0;
        for (int idim=0; idim<numdim; idim++)
        {
          // get global and local dof IDs
          const int gid = nodedofs[idim];
          const int lid = convel->Map().LID(gid);
          if (lid < 0) dserror("Local ID not found in map for given global ID!");

          double vel_i = (*convel)[lid];
          if (abs(vel_i) < 1e-12)
            no_slip_node++;
        }

        // yes, we have a no-slip node
        if (no_slip_node == numdim)
        {
          // do we also have a temperature dirichlet boundary condition
          // get the set of temperature degrees of freedom associated with the node
          std::vector<int> nodedofset = scatradiscret_->Dof(0,lnode);
          if (nodedofset.size()>1)
            dserror("Dynamic Smagorinsky or dynamic Vreman currently only implemented for one scalar field!");

          // check whether the dofs are Dirichlet constrained
          bool is_dirichlet_node = false;
          int gid = nodedofset[0];
          int lid = dofrowmap->LID(gid);

          //this is a dirichlet node
          if ((*dirichtoggle)[lid]==1)
            is_dirichlet_node = true;

          //get volume
          double thisvol = (*patchvol)[lnodeid];
          // and density
          double dens = (*filtered_dens_)[lnodeid]/thisvol;
          int err = 0;
          err += filtered_dens_->ReplaceMyValues(1,&dens,&lnodeid);

          double temp = 0.0;
          if (is_dirichlet_node)
          {
            temp  = (*scalar)[lid];
            err += filtered_temp_->ReplaceMyValues(1,&temp,&lnodeid);
            double val = dens*temp;
            err += filtered_dens_temp_->ReplaceMyValues(1,&val,&lnodeid);
            if (phi2_)
            {
              double val=0.0;
              err += filtered_phi2_->ReplaceMyValues(1,&val,&lnodeid);
            }
            if (phiexpression_)
            {
              double val=0.0;
              err += filtered_phiexpression_->ReplaceMyValues(1,&val,&lnodeid);
            }
          }
          else
          {
            temp = (*filtered_temp_)[lnodeid]/thisvol;
            err += filtered_temp_->ReplaceMyValues(1,&temp,&lnodeid);
            double val = (*filtered_dens_temp_)[lnodeid]/thisvol;
            err += filtered_dens_temp_->ReplaceMyValues(1,&val,&lnodeid);
          }

          for (int idim=0; idim<numdim; idim++)
          {
            // get global and local dof IDs
            const int gid = nodedofs[idim];
            const int lid = convel->Map().LID(gid);
            if (lid < 0) dserror("Local ID not found in map for given global ID!");

            double valvel_i = (*convel)[lid];
            err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&valvel_i,&lnodeid);

            double valdensvel_i = dens*valvel_i;
            err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&valdensvel_i,&lnodeid);

            double dvtval = dens*temp*valvel_i;
            err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&dvtval,&lnodeid);

            // Peter style
            double drtval = 0.0;
            err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&drtval,&lnodeid);

            if (phi_)
            {
              double drtval = 0.0;
              err += ((*filtered_phi_)(idim))->ReplaceMyValues(1,&drtval,&lnodeid);
            }
            if (alphaijsc_)
            {
              for (int jdim =0;jdim<numdim;++jdim)
              {
                const int ij = numdim*idim+jdim;
                if (no_slip_node == numdim)
                {
                  // set value to zero (original Peter style)
                  double val = 0.0;
                  err += ((*filtered_alphaijsc_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
                }
              else
                {
                  // set value to mean value
                  // we already divide by the corresponding volume of all contributing elements,
                  // since we set the volume to 1.0 in the next step in order not to modify the dirichlet values
                  double val = ((*((*filtered_alphaijsc_ ) (ij)))[lnodeid])/thisvol;
                  err += ((*filtered_alphaijsc_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
                }
              } // end loop jdim
            }

            // alternative: see comment in ApplyBoxFilter() for velocity field
            //double drtval = ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lnodeid])/thisvol;
            //err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&drtval,&lnodeid);
          }

          double volval = 1.0;
          err += patchvol->ReplaceMyValues(1,&volval,&lnodeid);
          if (err!=0) dserror("dof/node not on proc");
        }
      }
    }
  }

  // ---------------------------------------------------------------
  // scale vectors by element patch sizes --- this corresponds to
  // the normalization of the box filter function

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<scatradiscret_->NumMyRowNodes();++lnodeid)
  {
    double thisvol = (*patchvol)[lnodeid];

    int err = 0;
    double val = 0.0;

    val = (*filtered_temp_)[lnodeid]/thisvol;
    err += filtered_temp_->ReplaceMyValues(1,&val,&lnodeid);
    val = (*filtered_dens_)[lnodeid]/thisvol;
    err += filtered_dens_->ReplaceMyValues(1,&val,&lnodeid);
    val = (*filtered_dens_temp_)[lnodeid]/thisvol;
    err += filtered_dens_temp_->ReplaceMyValues(1,&val,&lnodeid);
    if (phi2_)
    {
      val = (*filtered_phi2_)[lnodeid]/thisvol;
      err += filtered_phi2_->ReplaceMyValues(1,&val,&lnodeid);
    }
    if (phiexpression_)
    {
      val = (*filtered_phiexpression_)[lnodeid]/thisvol;
      err += filtered_phiexpression_->ReplaceMyValues(1,&val,&lnodeid);
    }
    for (int idim =0;idim<3;++idim)
    {
      val = ((*((*filtered_vel_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      val = ((*((*filtered_dens_vel_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_dens_vel_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      val = ((*((*filtered_dens_vel_temp_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_dens_vel_temp_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      val = ((*((*filtered_dens_rateofstrain_temp_)(idim)))[lnodeid])/thisvol;
      err += ((*filtered_dens_rateofstrain_temp_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      if (phi_)
      {
        val = ((*((*filtered_phi_)(idim)))[lnodeid])/thisvol;
        err += ((*filtered_phi_)(idim))->ReplaceMyValues(1,&val,&lnodeid);
      }
      if (alphaijsc_)
      {
        for (int jdim =0;jdim<3;++jdim)
        {
          const int ij = numdim*idim+jdim;
          val = ((*((*filtered_alphaijsc_ ) (ij)))[lnodeid])/thisvol;
          err += ((*filtered_alphaijsc_ ) (ij))->ReplaceMyValues(1,&val,&lnodeid);
        } // end loop jdim
      }
    } // end loop idim
    if (err!=0) dserror("dof not on proc");
  } // end loop nodes

  // clean up
  scatradiscret_->ClearState();

  // ----------------------------------------------------------
  // the communication part: Export from row to column map

  // get the column map in order to communicate the result to all ghosted nodes
  const Epetra_Map* nodecolmap = scatradiscret_->NodeColMap();

  // allocate distributed vectors in col map format to have the filtered
  // quantities available on ghosted nodes
  col_filtered_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_dens_vel_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_dens_vel_temp_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_dens_rateofstrain_temp_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  col_filtered_temp_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  col_filtered_dens_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  col_filtered_dens_temp_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  if (phi_)
    col_filtered_phi_ = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,3,true));
  if(phi2_)
    col_filtered_phi2_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  if (phiexpression_)
    col_filtered_phiexpression_ = Teuchos::rcp(new Epetra_Vector(*nodecolmap,true));
  if (alphaijsc_)
      col_filtered_alphaijsc_         = Teuchos::rcp(new Epetra_MultiVector(*nodecolmap,9,true));

  // export filtered vectors in rowmap to columnmap format
  LINALG::Export(*filtered_vel_,*col_filtered_vel_);
  LINALG::Export(*filtered_dens_vel_,*col_filtered_dens_vel_);
  LINALG::Export(*filtered_dens_vel_temp_,*col_filtered_dens_vel_temp_);
  LINALG::Export(*filtered_dens_rateofstrain_temp_,*col_filtered_dens_rateofstrain_temp_);
  LINALG::Export(*filtered_temp_,*col_filtered_temp_);
  LINALG::Export(*filtered_dens_,*col_filtered_dens_);
  LINALG::Export(*filtered_dens_temp_,*col_filtered_dens_temp_);
  if (phi_)
    LINALG::Export(*filtered_phi_,*col_filtered_phi_);
  if (phi2_)
    LINALG::Export(*filtered_phi2_,*col_filtered_phi2_);
  if (phiexpression_)
    LINALG::Export(*filtered_phiexpression_,*col_filtered_phiexpression_);
  if (alphaijsc_)
      LINALG::Export(*filtered_alphaijsc_        ,*col_filtered_alphaijsc_        );

  return;
}



