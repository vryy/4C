/*!----------------------------------------------------------------------
\file patspec.cpp

\brief Methods for spring and dashpot constraints / boundary conditions:

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/


#include "springdashpot.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <iostream>
#include "../drt_io/io_pstream.H" // has to go before io.H
#include "../drt_io/io.H"
#include <Epetra_Time.h>

#include "../drt_adapter/adapter_coupling_nonlin_mortar.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*
 |                                                         pfaller Apr15|
 *----------------------------------------------------------------------*/
UTILS::SpringDashpot::SpringDashpot(
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<DRT::Condition> cond):
    actdisc_(dis),
    spring_(cond),
    stiff_tens_(spring_->GetDouble("SPRING_STIFF_TENS")),
    stiff_comp_(spring_->GetDouble("SPRING_STIFF_COMP")),
    offset_(spring_->GetDouble("SPRING_OFFSET")),
    viscosity_(spring_->GetDouble("DASHPOT_VISCOSITY")),
    coupling_(spring_->GetInt("coupling id")),
    nodes_(spring_->Nodes()),
    area_(),
    gap0_(),
    gap_(),
    gapdt_(),
    dgap_(),
    normals_(),
    dnormals_()
{
  // set type of this spring
  SetSpringType();

  if (springtype_ != cursurfnormal && coupling_>=0)
    dserror("Coupling of spring dashpot to reference surface only possible for DIRECTION cursurfnormal.");

  if (springtype_ == cursurfnormal && coupling_==-1)
    dserror("Coupling id necessary for DIRECTION cursurfnormal.");

  if (springtype_ == cursurfnormal && viscosity_ != 0)
    dserror("Viscous damping not yet implemented for DIRECTION cursurfnormal.");

  // get geometry
  std::map<int,Teuchos::RCP<DRT::Element> >& geom = spring_->Geometry();

  // calculate nodal area
  if (!actdisc_->Comm().MyPID()) IO::cout << "Computing area for spring dashpot condition...\n";
  GetArea(geom);

  // get normal vectors if necessary
  if (springtype_ == cursurfnormal)
    InitializeCurSurfNormal();
  else if (springtype_ == refsurfnormal)
    GetRefNormals(geom);

  return;
}

/*----------------------------------------------------------------------*
 |                                                             mhv 01/14|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpot::Evaluate(
    Teuchos::RCP<LINALG::SparseOperator> stiff,
    Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> velo,
    Teuchos::ParameterList parlist)
{
  if (disp==Teuchos::null) dserror("Cannot find displacement state in discretization");

  // reset last newton step (common for all spring types)
  gap_.clear();
  gapdt_.clear();
  springstress_.clear();

  // get time integrator properties
  const double gamma = parlist.get("scale_gamma",0.0);
  const double beta = parlist.get("scale_beta",1.0);
  const double ts_size = parlist.get("time_step_size",1.0);
  const double time_scale = gamma/(beta*ts_size);

  if (springtype_ == cursurfnormal)
  {
    GetCurNormals(disp);
    stiff->UnComplete(); // sparsity pattern might change
  }

  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j=0; j<(int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d",gid);

      // get nodal values
      const double nodalarea = area_[gid]; // nodal area
      const std::vector<double> normal = normals_[gid]; // normalized nodal normal

      const int numdof = actdisc_->NumDof(0,node);
      assert (numdof==3);
      std::vector<int> dofs = actdisc_->Dof(0,node);

      // initialize. calculation of normals and displacements differs for each spring variant
      std::vector<double> u(numdof); // displacement vector of condition nodes
      std::vector<double> v(numdof); // velocity vector of condition nodes
      if (springtype_ == refsurfnormal ||
          springtype_ == all)
      {
        for (int k=0; k<numdof; ++k)
        {
          u[k] = (*disp)[disp->Map().LID(dofs[k])]; // extract nodal displacement
          v[k] = (*velo)[velo->Map().LID(dofs[k])]; // extract nodal velocity
        }
      }

      // projection of displacement/velocity onto nodal reference normal
      if (springtype_ == refsurfnormal)
      {
        double gap = 0.; // projection of displacement vector to refsurfnormal (u \cdot N) = gap function
        double gapdt = 0.; // projection of velocity vector to refsurfnormal (v \cdot N) = gap function
        for (int k=0; k<numdof; ++k)
        {
          gap -= u[k]*normal[k]; // project displacement on normal
          gapdt -= v[k]*normal[k]; // project velocity on normal
        }

        // insert values into class variables for postprocessing
        gap_.insert(std::pair<int, double>(gid, gap));
        gapdt_.insert(std::pair<int, double>(gid, gapdt));
      }

      // assemble into residual and stiffness matrix
      // spring dashpot acts in every surface dof direction
      if (springtype_ == all)
      {
        for (int k=0; k<numdof; ++k)
        {
          if (stiff_tens_ != stiff_comp_)
            dserror("SPRING_STIFF_TENS != SPRING_STIFF_COMP: Different spring moduli for tension and compression not supported "
                "when specifying 'all' as DIRECTION (no ref surface normal information is calculated for that case)! "
                "Only possible for DIRECTION 'refsurfnormal' or 'cursurfnormal'.");

          const double val  = nodalarea*(stiff_tens_*(u[k]-offset_) + viscosity_*v[k]);
          const double dval = nodalarea*(stiff_tens_ + viscosity_*gamma/(beta*ts_size));

          const int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
          if (err) dserror("SumIntoGlobalValues failed!");
          stiff->Assemble(dval,dofs[k],dofs[k]);
        }
      }
      // spring dashpot acts in ref/cur normal direction
      else if (springtype_ == refsurfnormal ||
               springtype_ == cursurfnormal)
      {
        // spring displacement
        const double gap = gap_[gid];
        const double gapdt = gapdt_[gid];

        // select spring stiffnes
        double springstiff = 0.0;
        if (gap > 0)
          springstiff = stiff_tens_; // gap positive: tensile spring
        if (gap <= 0)
          springstiff = stiff_comp_; // gap negative: compressive spring

        // assemble
        std::vector<double> temp(numdof, 0.0);
        for (int k=0; k<numdof; ++k)
        {
          // force
          const double val = - nodalarea*(springstiff*(gap-offset_) + viscosity_*gapdt)*normal[k];
          const int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
          if (err) dserror("SumIntoGlobalValues failed!");

          // store negative value of internal force for output (=reaction force)
          temp[k] = - val;

          // stiffness
          if (springtype_ == refsurfnormal)
          {
            // linearize displacement and velocity
            for (int m=0; m<numdof; ++m)
            {
              const double dval = nodalarea*(springstiff + viscosity_*time_scale)*normal[k]*normal[m];
              stiff->Assemble(dval,dofs[k],dofs[m]);
            }
          }
          else if (springtype_ == cursurfnormal)
          {
            std::map<int,double> dgap = dgap_[gid];
            std::vector<GEN::pairedvector<int,double> > dnormal = dnormals_[gid];

            // check if projection exists
            if(!dnormal.empty() && !dgap.empty())
            {
              // linearize gap
              for(std::map<int, double>::iterator i = dgap.begin(); i != dgap.end(); ++i)
              {
                const double dval = -nodalarea*springstiff*(i->second)*normal[k];
                stiff->Assemble(dval,dofs[k],i->first);
              }

              // linearize normal
              for(GEN::pairedvector<int,double>::iterator i = dnormal[k].begin(); i != dnormal[k].end(); ++i)
              {
                const double dval = -nodalarea*springstiff*(gap-offset_)*(i->second);
                stiff->Assemble(dval,dofs[k],i->first);
              }
            }
          }
        }
        // add to output
        springstress_.insert(std::pair<int, std::vector<double> >(gid, temp));
      }
    } //node owned by processor
  } //loop over nodes

  if (springtype_ == cursurfnormal)
    stiff->Complete(); // sparsity pattern might have changed

  return;
}

/*----------------------------------------------------------------------*
 |                                                         pfaller Jan14|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpot::OutputGapNormal(
    Teuchos::RCP<Epetra_Vector> &gap,
    Teuchos::RCP<Epetra_MultiVector> &normals,
    Teuchos::RCP<Epetra_MultiVector> &stress)
{
  // export gap function
  for(std::map<int, double>::iterator i = gap_.begin(); i != gap_.end(); ++i)
  {
    // global id -> local id
    const int lid = gap->Map().LID(i->first);
    // local id on processor
    if (lid>=0)
      (*gap)[lid] += i->second;
  }

  // export normal
  for(std::map<int, std::vector<double> >::iterator i = normals_.begin(); i != normals_.end(); ++i)
  {
    // global id -> local id
    const int lid = normals->Map().LID(i->first);
    // local id on processor
    if (lid>=0)
    {
      // copy all components of normal vector
      (*(*normals)(0))[lid] += (i->second)[0];
      (*(*normals)(1))[lid] += (i->second)[1];
      (*(*normals)(2))[lid] += (i->second)[2];
    }
  }

  // export spring stress
  for(std::map<int, std::vector<double> >::iterator i = springstress_.begin(); i != springstress_.end(); ++i)
  {
    // global id -> local id
    const int lid = stress->Map().LID(i->first);
    // local id on processor
    if (lid>=0)
    {
      // copy all components of normal vector
      (*(*stress)(0))[lid] += (i->second)[0];
      (*(*stress)(1))[lid] += (i->second)[1];
      (*(*stress)(2))[lid] += (i->second)[2];
    }
  }

  return;
}

/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::InitializeCurSurfNormal()
{
  // create MORTAR interface
  mortar_ = Teuchos::rcp(new ADAPTER::CouplingNonLinMortar());

  // create CONTACT elements at interface for normal and gap calculation
  mortar_->SetupSpringDashpot(actdisc_, actdisc_, spring_, coupling_, actdisc_->Comm());

  // create temp vectors for gap initialization
  std::map<int,std::map<int,double> > tmpdgap_;
  std::map<int, std::vector<double> > tmpnormals_;
  std::map<int,std::vector<GEN::pairedvector<int,double> > > tmpdnormals_;

  // empty displacement vector
  Teuchos::RCP<Epetra_Vector> disp;
  disp = LINALG::CreateVector(*(actdisc_->DofRowMap()), true);

  // initialize gap in reference configuration
  mortar_->Interface()->EvaluateDistances(disp,tmpnormals_,tmpdnormals_,gap0_,tmpdgap_);

  return;
}

/*-----------------------------------------------------------------------*
|(private) adapted from mhv 01/14                           pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::GetArea(const std::map<int,Teuchos::RCP<DRT::Element> >& geom)
{
  std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();

    Teuchos::ParameterList eparams;

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    element->LocationVector(*(actdisc_),lm,lmowner,lmstride);
    Epetra_SerialDenseMatrix dummat(0,0);
    Epetra_SerialDenseVector dumvec(0);
    Epetra_SerialDenseVector elevector;
    const int eledim = (int)lm.size();
    elevector.Size(eledim);

    eparams.set("action","calc_struct_area");
    eparams.set("area",0.0);
    element->Evaluate(eparams,*(actdisc_),lm,dummat,dummat,dumvec,dumvec,dumvec);

    DRT::Element::DiscretizationType shape = element->Shape();

    double a = eparams.get("area",-1.0);

    // loop over all nodes of the element that share the area
    // do only contribute to my own row nodes
    double apernode = 0.;
    for (int i=0; i<element->NumNode(); ++i)
    {
      /* here we have to take care to assemble the right stiffness to the nodes!!! (mhv 05/2014):
          we do some sort of "manual" gauss integration here since we have to pay attention to assemble
          the correct stiffness in case of quadratic surface elements*/

      switch(shape)
      {
      case DRT::Element::tri3:
        apernode = a / element->NumNode();
        break;
      case DRT::Element::tri6:
      {
        //integration of shape functions over parameter element surface
        double int_N_cornernode = 0.;
        double int_N_edgemidnode = 1./6.;

        int numcornernode = 3;
        int numedgemidnode = 3;

        double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;

        //corner nodes
        if (i==0) apernode = int_N_cornernode * a / weight;
        if (i==1) apernode = int_N_cornernode * a / weight;
        if (i==2) apernode = int_N_cornernode * a / weight;
        //edge mid nodes
        if (i==3) apernode = int_N_edgemidnode * a / weight;
        if (i==4) apernode = int_N_edgemidnode * a / weight;
        if (i==5) apernode = int_N_edgemidnode * a / weight;
      }
      break;
      case DRT::Element::quad4:
        apernode = a / element->NumNode();
        break;
      case DRT::Element::quad8:
      {
        //integration of shape functions over parameter element surface
        double int_N_cornernode = -1./3.;
        double int_N_edgemidnode = 4./3.;

        int numcornernode = 4;
        int numedgemidnode = 4;

        double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;

        //corner nodes
        if (i==0) apernode = int_N_cornernode * a / weight;
        if (i==1) apernode = int_N_cornernode * a / weight;
        if (i==2) apernode = int_N_cornernode * a / weight;
        if (i==3) apernode = int_N_cornernode * a / weight;
        //edge mid nodes
        if (i==4) apernode = int_N_edgemidnode * a / weight;
        if (i==5) apernode = int_N_edgemidnode * a / weight;
        if (i==6) apernode = int_N_edgemidnode * a / weight;
        if (i==7) apernode = int_N_edgemidnode * a / weight;
      }
      break;
      case DRT::Element::quad9:
      {
        //integration of shape functions over parameter element surface
        double int_N_cornernode = 1./9.;
        double int_N_edgemidnode = 4./9.;
        double int_N_centermidnode = 16./9.;

        int numcornernode = 4;
        int numedgemidnode = 4;
        int numcentermidnode = 1;

        double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode + numcentermidnode * int_N_centermidnode;

        //corner nodes
        if (i==0) apernode = int_N_cornernode * a / weight;
        if (i==1) apernode = int_N_cornernode * a / weight;
        if (i==2) apernode = int_N_cornernode * a / weight;
        if (i==3) apernode = int_N_cornernode * a / weight;
        //edge mid nodes
        if (i==4) apernode = int_N_edgemidnode * a / weight;
        if (i==5) apernode = int_N_edgemidnode * a / weight;
        if (i==6) apernode = int_N_edgemidnode * a / weight;
        if (i==7) apernode = int_N_edgemidnode * a / weight;
        //center mid node
        if (i==8) apernode = int_N_centermidnode * a / weight;
      }
      break;
      case DRT::Element::nurbs9:
        dserror("Not yet implemented for Nurbs! To do: Apply the correct weighting of the area per node!");
        break;
      default:
        dserror("shape type unknown!\n");
        break;
      }

      const int gid = element->Nodes()[i]->Id();
      if (!actdisc_->NodeRowMap()->MyGID(gid)) continue;

      // store area in map (gid, area). erase old value before adding new one
      const double newarea = area_[gid] + apernode;
      area_.erase(gid);
      area_.insert(std::pair<int, double>(gid, newarea));
    }
  } // for (ele=geom.begin(); ele != geom.end(); ++ele)

  return;
}

/*-----------------------------------------------------------------------*
|(private) adapted from mhv 01/14                           pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::GetRefNormals(const std::map<int,Teuchos::RCP<DRT::Element> >& geom)
{
  // a vector for all row dofs to hold normals interpolated to the nodes
  Epetra_Vector refnodalnormals(*actdisc_->DofRowMap(),true);

  std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    element->LocationVector(*(actdisc_),lm,lmowner,lmstride);
    Epetra_SerialDenseMatrix dummat(0,0);
    Epetra_SerialDenseVector dumvec(0);
    Epetra_SerialDenseVector elevector;
    const int eledim = (int)lm.size();
    elevector.Size(eledim);

    Teuchos::ParameterList eparams2;
    eparams2.set("action","calc_ref_nodal_normals");
    element->Evaluate(eparams2,*(actdisc_),lm,dummat,dummat,elevector,dumvec,dumvec);
    LINALG::Assemble(refnodalnormals,elevector,lm,lmowner);
  }

  const std::vector<int>& nds = *nodes_;
  for (int j=0; j<(int)nds.size(); ++j)
  {
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];

      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d",gid);

      int numdof = actdisc_->NumDof(node);
      std::vector<int> dofs = actdisc_->Dof(node);

      assert (numdof==3);

      std::vector<double> temp(numdof, 0.0);
      double temp_norm_sq = 0.;
      for (int k=0; k<numdof; ++k)
      {
        temp[k] = refnodalnormals[refnodalnormals.Map().LID(dofs[k])];
        temp_norm_sq += temp[k]*temp[k];
      }

      // calculate vector length
      const double temp_norm = sqrt(temp_norm_sq);

      if (temp_norm == 0.)
        dserror("Nodal normal has length 0.");

      // normalize vector
      for (int k=0; k<numdof; ++k)
        temp[k] /= temp_norm;

      // insert to map
      normals_.insert(std::pair<int, std::vector<double> >(gid, temp));
    }
  }

  return;
}


/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::GetCurNormals(Teuchos::RCP<Epetra_Vector> disp)
{
  // reset last newton step (only for curnormal)
  dgap_.clear();
  normals_.clear();
  dnormals_.clear();

  // temp nodal gap
  std::map<int, double> tmpgap;

  // calculate normals and gap using CONTACT elements
  mortar_->Interface()->EvaluateDistances(disp,normals_,dnormals_,tmpgap,dgap_);

  // subtract reference gap from current gap (gap in reference configuration is zero everywhere)
  for(std::map<int, double>::iterator i = tmpgap.begin(); i != tmpgap.end(); ++i)
  {
    std::map<int, double>::iterator j = gap0_.find(i->first);
    if(j == gap0_.end())
      gap_[i->first] = i->second;
//      dserror("The maps of reference gap and current gap are inconsistent.");
    else
      gap_[i->first] = i->second - j->second;
  }

  return;
}

/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::SetSpringType()
{
  // get spring direction from condition
  const std::string* dir = spring_->Get<std::string>("direction");

  if (*dir=="all")
    springtype_ = all;
  else if (*dir=="refsurfnormal")
    springtype_ = refsurfnormal;
  else if (*dir=="cursurfnormal")
    springtype_ = cursurfnormal;
  else
    dserror("Invalid direction option! Choose DIRECTION all, DIRECTION refsurfnormal or DIRECTION cursurfnormal!");
}
