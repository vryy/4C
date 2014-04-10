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



/*----------------------------------------------------------------------*
 |                                                             mhv 01/14|
 *----------------------------------------------------------------------*/
void SPRINGDASHPOT::SpringDashpot(Teuchos::RCP<DRT::Discretization> dis)
{

  //------------test discretization for presence of spring dashpot tissue condition
  std::vector<DRT::Condition*> springdashpotcond;

  // this is a copy of the discretization passed in to hand it to the element level.
  DRT::Discretization& actdis = *(dis);
  dis->GetCondition("SpringDashpot",springdashpotcond);
  if ((int)springdashpotcond.size())
  {
    //params->set("havespringdashpot", true); // so the condition is evaluated in STR::ApplyForceStiffSpringDashpot
    if (!dis->Comm().MyPID()) IO::cout << "Computing area for spring dashpot condition...\n";

    // loop over all spring dashpot conditions
    for (int cond=0; cond<(int)springdashpotcond.size(); ++cond)
    {

      const std::string* dir = springdashpotcond[cond]->Get<std::string>("direction");

      // a vector for all row nodes to hold element area contributions
      Epetra_Vector nodalarea(*dis->NodeRowMap(),true);
      // a vector for all row dofs to hold normals interpolated to the nodes
      Epetra_Vector refnodalnormals(*dis->DofRowMap(),true);

      std::map<int,Teuchos::RCP<DRT::Element> >& geom = springdashpotcond[cond]->Geometry();
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator ele;
      for (ele=geom.begin(); ele != geom.end(); ++ele)
      {
        DRT::Element* element = ele->second.get();

        Teuchos::ParameterList eparams;
        Teuchos::ParameterList eparams2;
        eparams.set("action","calc_struct_area");
        eparams.set("area",0.0);
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        element->LocationVector(actdis,lm,lmowner,lmstride);
        const int eledim = (int)lm.size();
        Epetra_SerialDenseMatrix dummat(0,0);
        Epetra_SerialDenseVector dumvec(0);
        Epetra_SerialDenseVector elevector;
        elevector.Size(eledim);
        element->Evaluate(eparams,actdis,lm,dummat,dummat,dumvec,dumvec,dumvec);

        // when refsurfnormal direction is chosen, call the respective element evaluation routine and assemble
        if (*dir == "direction_refsurfnormal")
        {
          eparams2.set("action","calc_ref_nodal_normals");
          element->Evaluate(eparams2,actdis,lm,dummat,dummat,elevector,dumvec,dumvec);
          LINALG::Assemble(refnodalnormals,elevector,lm,lmowner);
        }

        double a = eparams.get("area",-1.0);

        // loop over all nodes of the element that share the area
        // do only contribute to my own row nodes
        const double apernode = a / element->NumNode();
        for (int i=0; i<element->NumNode(); ++i)
        {
          int gid = element->Nodes()[i]->Id();
          if (!dis->NodeRowMap()->MyGID(gid)) continue;
          nodalarea[dis->NodeRowMap()->LID(gid)] += apernode;
        }
      } // for (ele=geom.begin(); ele != geom.end(); ++ele)

      // now we have the area per node, put it in a vector that is equal to the nodes vector
      // consider only my row nodes
      const std::vector<int>* nodes = springdashpotcond[cond]->Nodes();

      std::vector<double> apern(nodes->size(),0.0);
      for (int i=0; i<(int)nodes->size(); ++i)
      {
        int gid = (*nodes)[i];
        if (!nodalarea.Map().MyGID(gid)) continue;
        apern[i] = nodalarea[nodalarea.Map().LID(gid)];
      }
      // set vector to the condition
      (*springdashpotcond[cond]).Add("areapernode",apern);


      if (*dir == "direction_refsurfnormal")
      {
        std::vector<double> refndnorms(3 * nodes->size(),0.0);

        const std::vector<int>& nds = *nodes;
        for (int j=0; j<(int)nds.size(); ++j)
        {
          if (dis->NodeRowMap()->MyGID(nds[j]))
          {
            int gid = nds[j];

            DRT::Node* node = dis->gNode(gid);
            if (!node) dserror("Cannot find global node %d",gid);

            int numdof = dis->NumDof(node);
            std::vector<int> dofs = dis->Dof(node);

            assert (numdof==3);

            for (int k=0; k<numdof; ++k)
            {
              //if (!refnodalnormals.Map().MyGID(dofs[k])) continue;
              refndnorms[numdof*j+k] = refnodalnormals[refnodalnormals.Map().LID(dofs[k])];
            }
          }
        }
        // set vector to the condition
        (*springdashpotcond[cond]).Add("refnodalnormals",refndnorms);
      }

    } // for (int cond=0; cond<(int)springdashpotcond.size(); ++cond)

  }// if ((int)springdashpotcond.size())
  //-------------------------------------------------------------------------


  return;
}


/*----------------------------------------------------------------------*
 |                                                             mhv 01/14|
 *----------------------------------------------------------------------*/
void SPRINGDASHPOT::EvaluateSpringDashpot(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<LINALG::SparseOperator> stiff,
    Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> velo,
    Teuchos::ParameterList parlist)
{

  if (disp==Teuchos::null) dserror("Cannot find displacement state in discretization");

  double gamma = parlist.get("scale_gamma",0.0);
  double beta = parlist.get("scale_beta",1.0);
  double ts_size = parlist.get("time_step_size",1.0);

  std::vector<DRT::Condition*> springdashpotcond;
  discret->GetCondition("SpringDashpot",springdashpotcond);


  for (int i=0; i<(int)springdashpotcond.size(); ++i)
  {
    const std::vector<int>* nodes = springdashpotcond[i]->Nodes();
    double springstiff = springdashpotcond[i]->GetDouble("spring_stiff");
    double springoffset = springdashpotcond[i]->GetDouble("spring_offset");
    double dashpotvisc = springdashpotcond[i]->GetDouble("dashpot_visc");
    const std::string* dir = springdashpotcond[i]->Get<std::string>("direction");

    const std::vector<double>* areapernode = springdashpotcond[i]->Get< std::vector<double> >("areapernode");
    const std::vector<double>* refnodalnormals = springdashpotcond[i]->Get< std::vector<double> >("refnodalnormals");


    const std::vector<int>& nds = *nodes;
    for (int j=0; j<(int)nds.size(); ++j)
    {

      if (discret->NodeRowMap()->MyGID(nds[j]))
      {
        int gid = nds[j];
        double nodalarea = (*areapernode)[j];

        DRT::Node* node = discret->gNode(gid);

        if (!node) dserror("Cannot find global node %d",gid);


        int numdof = discret->NumDof(node);
        std::vector<int> dofs = discret->Dof(node);

        assert (numdof==3);

        // displacement vector of condition nodes
        std::vector<double> u(numdof);
        for (int k=0; k<numdof; ++k)
        {
          u[k] = (*disp)[disp->Map().LID(dofs[k])];
        }

        // velocity vector of condition nodes
        std::vector<double> v(numdof);
        for (int k=0; k<numdof; ++k)
        {
          v[k] = (*velo)[velo->Map().LID(dofs[k])];
        }

        // assemble into residual and stiffness matrix for case that spring / dashpot acts in every surface dof direction
        if (*dir == "direction_all")
        {
          for (int k=0; k<numdof; ++k)
          {
            double val = nodalarea*(springstiff*(u[k]-springoffset) + dashpotvisc*v[k]);
            int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
            if (err) dserror("SumIntoGlobalValues failed!");
            stiff->Assemble(nodalarea*(springstiff + dashpotvisc*gamma/(beta*ts_size)),dofs[k],dofs[k]);
          }
        }
        // assemble into residual and stiffness matrix for case that spring / dashpot acts in reference surface normal direction
        else if (*dir == "direction_refsurfnormal")
        {
          // extract averaged nodal ref normal and compute its absolute value
          std::vector<double> unitrefnormal(numdof);
          double temp_ref = 0.;
          for (int k=0; k<numdof; ++k)
          {
            unitrefnormal[k] = (*refnodalnormals)[numdof*j+k];
            temp_ref += unitrefnormal[k]*unitrefnormal[k];
          }

          double unitrefnormalabsval = sqrt(temp_ref);

          for(int k=0; k<numdof; ++k)
          {
            unitrefnormal[k] /= unitrefnormalabsval;
          }

          //dyadic product of ref normal with itself (N \otimes N)
          Epetra_SerialDenseMatrix N_x_N(numdof,numdof);
          for (int l=0; l<numdof; ++l)
            for (int m=0; m<numdof; ++m)
              N_x_N(l,m)=unitrefnormal[l]*unitrefnormal[m];

          for (int k=0; k<numdof; ++k)
          {
            for (int m=0; m<numdof; ++m)
            {
              double val = nodalarea*N_x_N(k,m)*(springstiff*(u[m]-springoffset) + dashpotvisc*v[m]);
              int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
              if (err) dserror("SumIntoGlobalValues failed!");
              stiff->Assemble(nodalarea*(springstiff + dashpotvisc*gamma/(beta*ts_size))*N_x_N(k,m),dofs[k],dofs[m]);
            }
          }
        }
        else dserror("Invalid direction option! Choose direction_all or direction_refsurfnormal!");

      } //node owned by processor

    } //loop over nodes

  } //loop over conditions

  return;
}



