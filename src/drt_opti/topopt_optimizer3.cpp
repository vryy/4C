/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer3.cpp

\brief optimizer of the topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::EvaluateObjective3()
{
  TEUCHOS_FUNC_TIME_MONITOR(" evaluate objective value");

  obj_value_ = 0.0;

//  DRT::Node* node;
//  double value = 0.0;
//  const double dt = fluidParams_->get<double>("time step size");
//
//  Teuchos::RCP<Epetra_Vector> fluidvel = rcp(new Epetra_Vector(*discret_->DofRowMap(),false));
//  Teuchos::RCP<Epetra_Vector> adjointvel = rcp(new Epetra_Vector(*discret_->DofRowMap(),false));
//
//  const int dim = DRT::Problem::Instance()->NDim();
//  std::vector<double> nodalfluidvel(dim);
//  std::vector<double> nodaladjointvel(dim);
//
//  double dissipation_fac = 0.0;
//  if (fluidAdjoint3Parameter_->ObjDissipationTerm()) dissipation_fac = fluidAdjoint3Parameter_->ObjDissipationFac();
//
//  for (int inode=0;inode<discret_->NumMyRowNodes();inode++)
//  {
//    node = discret_->lRowNode(inode);
//
//    vector<int> lm = discret_->Dof(node);
//    lm.pop_back(); // delete pressure dof
//
//    value = 0.0;
//
//    if (fluidAdjoint3Parameter_->IsStationary())
//    {
//      fluidvel = vel_->find(1)->second;
//      adjointvel = adjointvel_->find(1)->second;
//
//      DRT::UTILS::ExtractMyValues(*fluidvel,nodalfluidvel,lm);
//      DRT::UTILS::ExtractMyValues(*adjointvel,nodaladjointvel,lm);
//
//      for (int idim=0;idim<dim;idim++)
//        value +=
//            dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
//            -nodalfluidvel[idim]*nodaladjointvel[idim]; // adjoint part
//    }
//    else
//    {
//      for (size_t timestep=0;timestep<=num_timesteps_;timestep++)
//      {
//        fluidvel = vel_->find(timestep)->second;
//        adjointvel = adjointvel_->find(timestep)->second;
//
//        DRT::UTILS::ExtractMyValues(*fluidvel,nodalfluidvel,lm);
//        DRT::UTILS::ExtractMyValues(*adjointvel,nodaladjointvel,lm);
//
//        if (timestep!=0 && timestep!=num_timesteps_) // default case, weight 1
//        {
//          for (int idim=0;idim<dim;idim++)
//            value +=
//                dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
//                -nodalfluidvel[idim]*nodaladjointvel[idim]; // adjoint part
//        }
//        else // first and last time step, weight 0.5
//        {
//          for (int idim=0;idim<dim;idim++)
//            value += 0.5*(
//                dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
//                -nodalfluidvel[idim]*nodaladjointvel[idim]); // adjoint part
//        }
//      }
//
//      value *= dt; // scale with time step size
//    }
//
//    int err = obj_grad_->SumIntoMyValue(inode,0,value);
//    if (err)
//      dserror("error while adding value to gradient of objective");
//  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TOPOPT::Optimizer::EvaluateGradient3()
{
  TEUCHOS_FUNC_TIME_MONITOR(" evaluate objective gradient");
  cout << "WARNING: THIS IS CURRENTLY NOT A WEAK GRADIENT, BUT ITS STRONG REPRESENTATION" << endl;
  obj_grad_->PutScalar(0.0);

  DRT::Node* node;
  double value = 0.0;
  const double dt = fluidParams_->get<double>("time step size");

  Teuchos::RCP<Epetra_Vector> fluidvel = rcp(new Epetra_Vector(*discret_->DofRowMap(),false));
  Teuchos::RCP<Epetra_Vector> adjointvel = rcp(new Epetra_Vector(*discret_->DofRowMap(),false));

  const int dim = DRT::Problem::Instance()->NDim();
  std::vector<double> nodalfluidvel(dim);
  std::vector<double> nodaladjointvel(dim);

  double dissipation_fac = 0.0;
  if (fluidParams_->get<bool>("OBJECTIVE_DISSIPATION")) dissipation_fac = fluidParams_->get<double>("DISSIPATION_FAC");

  for (int inode=0;inode<discret_->NumMyRowNodes();inode++)
  {
    node = discret_->lRowNode(inode);

    vector<int> lm = discret_->Dof(node);
    lm.pop_back(); // delete pressure dof

    value = 0.0;

    if (fluidParams_->get<int>("time int algo")==INPAR::FLUID::timeint_stationary)
    {
      fluidvel = vel_->find(1)->second;
      adjointvel = adjointvel_->find(1)->second;

      DRT::UTILS::ExtractMyValues(*fluidvel,nodalfluidvel,lm);
      DRT::UTILS::ExtractMyValues(*adjointvel,nodaladjointvel,lm);

      for (int idim=0;idim<dim;idim++)
        value +=
            dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
            -nodalfluidvel[idim]*nodaladjointvel[idim]; // adjoint part
    }
    else
    {
      for (int timestep=0;timestep<=fluidParams_->get<int>("max number timesteps");timestep++)
      {
        fluidvel = vel_->find(timestep)->second;
        adjointvel = adjointvel_->find(timestep)->second;

        DRT::UTILS::ExtractMyValues(*fluidvel,nodalfluidvel,lm);
        DRT::UTILS::ExtractMyValues(*adjointvel,nodaladjointvel,lm);
        if (inode==6) cout << "initial value is " << value << endl;
        if (timestep!=0 && timestep!=fluidParams_->get<int>("max number timesteps")) // default case, weight 1
        {
          for (int idim=0;idim<dim;idim++)
          {
            value +=
                dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
                -nodalfluidvel[idim]*nodaladjointvel[idim]; // adjoint part
          }
        }
        else // first and last time step, weight 0.5
        {
          for (int idim=0;idim<dim;idim++)
          {
            value += 0.5*(
                dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
                -nodalfluidvel[idim]*nodaladjointvel[idim]); // adjoint part
          }
        }
      }
      value *= dt; // scale with time step size
    }

    int err = obj_grad_->SumIntoMyValue(inode,0,value);
    if (err)
      dserror("error while adding value to gradient of objective");
  }
}



