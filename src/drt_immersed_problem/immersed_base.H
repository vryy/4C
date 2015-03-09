/*!----------------------------------------------------------------------
\file immersed_base.cpp

\brief base class for all immersed algorithms

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_BASE_H
#define IMMERSED_BASE_H

#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_element.H" // todo remove when moving interpolatetoimmersedintpoint to drt_lib
#include "../drt_so3/so_hex8.H"
#include "../drt_fluid_ele/fluid_ele_immersed.H"
#include "../drt_immersed_problem/immersed_node.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_calc_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_lib/drt_exporter.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

namespace DRT{
  class Discretization;
}

namespace IMMERSED{

class ImmersedBase
{
public:

  /// setup base class for all immersed algorithms
  explicit ImmersedBase();

  std::vector<int> DetermineImmersionDomain(Teuchos::RCP<DRT::Discretization> backgrounddis,
                                            Teuchos::RCP<DRT::Discretization> immerseddis,
                                            DRT::AssembleStrategy* strategy,
                                            bool firstcall=false
                                            );

  std::vector<int> DetermineImmersionBoundaryDomain(Teuchos::RCP<DRT::Discretization> backgrounddis,
                                                    Teuchos::RCP<DRT::Discretization> immerseddis,
                                                    const std::string& condname,
                                                    bool gpversion = false,
                                                    bool firstcall=false,
                                                    const Teuchos::RCP<std::map<int,Teuchos::RCP<std::map<int,Teuchos::RCP<std::vector<double> > > > > >& gpmap = Teuchos::null
                                                    );

  void CreateGhosting(const Teuchos::RCP<DRT::Discretization> distobeghosted);

  void CreateVolumeCondition(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, DRT::Condition::ConditionType condtype, std::string condname);

  void UpdateVolumeCondition(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, DRT::Condition::ConditionType condtype, std::string condname);

  void DetermineDirichletDofsFromNodes(const Teuchos::RCP<DRT::Discretization>& discret, const std::vector<int>& nodeids, std::vector<int>& vectofill);

  void EvaluateWithInternalCommunication(Teuchos::RCP<DRT::Discretization> dis, DRT::AssembleStrategy* strategy, std::map<int,std::set<int> >& elementstoeval, Teuchos::RCP<GEO::SearchTree> structsearchtree, std::map<int,LINALG::Matrix<3,1> >& currpositions_struct);

  /*----------------------------------------------------------------------*
   * considering the need of communication inside the loop over all       |
   * conditioned elements                                     rauch 05/14 |
   *----------------------------------------------------------------------*/
  void EvaluateInterpolationCondition
  (
      Teuchos::RCP<DRT::Discretization> evaldis,
      Teuchos::ParameterList& params,
      DRT::AssembleStrategy & strategy,
      const std::string& condstring,
      const int condid
  );

};// class ImmersedBase

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
  template < DRT::Element::DiscretizationType sourcedistype,
             DRT::Element::DiscretizationType targetdistype >
  int InterpolateToImmersedIntPoint(const Teuchos::RCP<DRT::Discretization> sourcedis,
                                                            const Teuchos::RCP<DRT::Discretization> targetdis,
                                                            DRT::Element& targetele,
                                                            const std::vector<double>& targetxi,
                                                            const std::vector<double>& targetedisp,
                                                            int action,
                                                            std::vector<double>& targetdata
                                                            )
  {
    //////////////////////////////////////////////////////////////////////////////////
    ///////
    ///////  SOURCE  -->  TARGET
    ///////
    ///////  maps source values to target discretization.
    ///////
    ///////  comments on parallelism:
    ///////
    ///////  every proc needs to enter this routine to make communication work.
    ///////  in the extreme scenario where a proc owns no target elements at all
    ///////  no evaluate will be called and this proc doesn't enter this routine
    ///////  which means in turn, that a Barrier() would never be reached. in
    ///////  this case the attempt to communicate would definitely fail.
    ///////  so make sure, that each proc owns at least one target element and/or
    ///////  it is also to be assured by external measures, that each proc
    ///////  reenters this routine until the last element on every proc is
    ///////  evaluated (e.g. dummy call of Evaluate(...)).
    ///////
    /////////////////////////////////////////////////////////////////////////////////

    // error flag
    int err = 0;

    // get communicator
    const Epetra_Comm& comm = sourcedis->Comm();

    // get the global problem
    DRT::Problem* problem = DRT::Problem::Instance();

    // dimension of global problem
    const int globdim = problem->NDim();

    // pointer to background element
    DRT::Element* sourceele;

    // declare vector of global current coordinates
    double x [globdim];

    // dimension of source discretization
    static const int source_dim = DRT::UTILS::DisTypeToDim<sourcedistype>::dim;
    static const int target_dim = DRT::UTILS::DisTypeToDim<targetdistype>::dim;

    // declarations and definitions for round robin loop
    const int numproc  = comm.NumProc();
    const int myrank   = comm.MyPID();                     // myrank
    const int torank   = (myrank + 1) % numproc;           // sends to
    const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    DRT::Exporter exporter(comm);

    // get current global coordinates of the given point xi of the target dis
    MORTAR::UTILS::LocalToCurrentGlobal<targetdistype>(targetele,globdim,&targetxi[0],targetedisp,&x[0]);
    std::vector<double> xvec(globdim);
    xvec[0]=x[0];
    xvec[1]=x[1];
    xvec[2]=x[2];

    // vector to fill by sourcedis->Evaluate (needs to be resized in calc class)
    Teuchos::RCP<Epetra_SerialDenseVector> vectofill = Teuchos::rcp(new Epetra_SerialDenseVector(6));
    (*vectofill)(0)=-1234.0;

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Epetra_SerialDenseVector> xi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(source_dim));
    Teuchos::RCP<Epetra_SerialDenseVector> targetxi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(target_dim));

    Teuchos::RCP<Epetra_SerialDenseVector> normal_at_targetpoint = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    std::vector<double> normal_vec(3);// normal_vec[0]=0.0;normal_vec[1]=0.0; normal_vec[2]=0.0;

    Epetra_SerialDenseMatrix dummy1;
    Epetra_SerialDenseMatrix dummy2;
    Epetra_SerialDenseVector dummy3;

    Teuchos::ParameterList params;

    params.set<std::string>("action", "calc_cur_normal_at_point");
    for(int i=0; i<target_dim;++i)
      (*targetxi_dense)(i) = targetxi[i];

    DRT::Element::LocationArray targetla(1);
    targetele.LocationVector(*targetdis,targetla,false);

    targetele.Evaluate(params,*targetdis,targetla,dummy1,dummy2,*normal_at_targetpoint,*targetxi_dense,dummy3);
    normal_at_targetpoint->Scale(1.0/(normal_at_targetpoint->Norm2()));
    normal_vec[0] = (*normal_at_targetpoint)(0);
    normal_vec[1] = (*normal_at_targetpoint)(1);
    normal_vec[2] = (*normal_at_targetpoint)(2);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////
    //////     round robin loop
    //////
    ////// in this loop every processor sends the current coordinates of the given target
    ////// point xi to his neighbors and checks whether the recieved coordinates lie within
    ////// one of the rank's searchbox elements.
    //////
    ////// if so, this neighbor interpolates the necessary quantities to this point.
    ////// the results, along with the information whether the point has already been
    ////// matched and with the id of the rank who requested the data originally is sent
    ////// around.
    //////
    ////// when the loop ends, all information are back on the proc who requested the data
    ////// originally and for every point the data should be stored in rdata.
    //////
    //////
    ////// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    //////
    //////////////////////////////////////////////////////////////////////////////////////////////
    // given point already matched?
    int matched = 0;
    // owner of given point
    int owner = myrank;
    // length of vectofill
    int datalength = (int)vectofill->Length();

    // get possible elements being intersected by immersed structure
    DRT::Condition* searchbox = sourcedis->GetCondition("ImmersedSearchbox");
    std::map<int,Teuchos::RCP<DRT::Element> >& searchboxgeom = searchbox->Geometry();

    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      //////////////////////////////////////////////////////////////////////
      // loop over all background elements
      //
      // determine background element in which the given point lies.
      //
      // if the current proc can already match the given point xi (xvec)
      //
      /////////////////////////////////////////////////////////////////////
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
      LINALG::Matrix<source_dim,1> xi(true);
      if(!matched)
      {
        // every proc that has searchboxgeom elements
        for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
        {
          bool converged = false;

          // get parameter space coords xi in source element
          MORTAR::UTILS::GlobalToLocal<sourcedistype>(*(curr->second),&xvec[0],&xi(0),converged);

//  # ifdef DEBUG
//          if(!converged)
//          {
//            dserror("MAPPING FROM SOURCE POINT TO TARGET DID NOT CONVERGE");
//          }
//  # endif
          /////////////////////////////////////////////////
          ////
          ////     MATCH
          ////
          /////////////////////////////////////////////////
          double tol = 1e-12;
          if ((abs(xi(0))-1.0)<tol and (abs(xi(1))-1.0)<tol and (abs(xi(2))-1.0)<tol)
          { // given point lies in element curr

            bool validsource = true;
            double scalarproduct= 0.0;
            Teuchos::RCP<Epetra_SerialDenseVector> vector;

            sourceele = curr->second.getRawPtr();

            // parameter list
            params.set<int>("action", action);
            // fill locationarray
            DRT::Element::LocationArray sourcela(1);
            sourceele->LocationVector(*sourcedis,sourcela,false);

            (*xi_dense)(0)=xi(0);
            (*xi_dense)(1)=xi(1);
            (*xi_dense)(2)=xi(2);

            // does targetpoint lie on an element edge?
            // only the sourceele located outward from the boundary should interpolate
            if( (abs(xi(0))<(1.0+tol) and abs(xi(0))>(1.0-tol)) or (abs(xi(1))<(1.0+tol) and abs(xi(1))>(1.0-tol)) or (abs(xi(2))<(1.0+tol) and abs(xi(2))>(1.0-tol)))
            {
//              // calc sourcele diagonal length
//              const double* X0 = sourceele->Nodes()[0]->X();
//              const double* X6 = sourceele->Nodes()[1]->X();
//              double srcelediaglength = sqrt(pow(X0[0]-X6[0],2)+pow(X0[1]-X6[1],2)+pow(X0[2]-X6[2],2));

              // get nodal coords of source ele
              vector = Teuchos::rcp(new Epetra_SerialDenseVector(3));
              std::vector<double> distances(sourceele->NumNode());
              for(int sourcenode=0; sourcenode<sourceele->NumNode(); sourcenode++)
              {
                const double* X;
                X = sourceele->Nodes()[sourcenode]->X();

                distances[sourcenode] = sqrt(pow(xvec[0]-X[0],2)+pow(xvec[1]-X[1],2)+pow(xvec[2]-X[2],2));

              }

              // get max. distance and corresponding index of src node
              int maxdistanceindex = -1;
              double maxvalue = 0.0;
              for(int i=0; i<(int)distances.size();++i)
              {
                if(distances[i]>maxvalue)
                {
                  maxvalue = distances[i];
                  maxdistanceindex = i;
                }
              }

              const double* X;
              X = sourceele->Nodes()[maxdistanceindex]->X();
              // build vector directed from targetele point (xvec) to sourceel node (X)
              // xvec should not lie at same ele surface as sourcenode
              (*vector)(0) = X[0]-xvec[0];
              (*vector)(1) = X[1]-xvec[1];
              (*vector)(2) = X[2]-xvec[2];
              vector->Scale(vector->Norm2());

              // build scalar product between normal and vector
              scalarproduct = (*vector)(0) * normal_vec[0] + (*vector)(1) * normal_vec[1] + (*vector)(2) * normal_vec[2];

              if(scalarproduct > 0.0)
                validsource = true;
              else
                validsource = false;

#ifdef DEBUG
              if(abs(scalarproduct) < 1e-13 and validsource)
              {
                std::cout<<"scalarproduct = "<<scalarproduct<<std::endl;
                dserror("normal vec and vec from point to node perpendicular.\n"
                        "no robust determination of outward lying src ele possible");
              }
              if(abs(scalarproduct) < -1e-12 and validsource)
              {
                std::cout<<"scalarproduct = "<<scalarproduct<<std::endl;
                dserror("normal vec and vec from point to node have no sharp angle.\n"
                        "no robust determination of outward laying src ele possible");
              }
#endif

              //count++;
            } // if on element edge


            // evaluate and perform action handed in to this function
            params.set<int>("action", action);
            if(validsource)
            {
//              std::cout<<"\n ---------- \n PROC "<<comm.MyPID()<<" ; sourcele : "<<sourceele->Id()<<std::endl;
//              std::cout<<"POINT: ["<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<"]"<<std::endl;
//              std::cout<<"normal: ["<<normal_vec[0]<<" "<<normal_vec[1]<<" "<<normal_vec[2]<<"]"<<std::endl;
//              std::cout<<"vector: "<<*vector<<std::endl;
//              std::cout<<"scalarproduct: "<<scalarproduct<<std::endl;
              sourceele->Evaluate(params,*sourcedis,sourcela,dummy1,dummy2,*vectofill ,*xi_dense,dummy3);
              matched = 1;

              Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersed>(curr->second)->SetIsImmersedBoundary(1);

              //std::cout<<"Tgt. Ele: "<<targetele.Id()<<" lies in src. ele: "<<curr->second->Id()<<" at  pos. ["<<std::setprecision(5)<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<"]"<<std::endl;

              DRT::Node** nodes = curr->second->Nodes();
              for(int node=0;node<curr->second->NumNode();++node)
                static_cast<IMMERSED::ImmersedNode* >(nodes[node])->SetIsImmersedBoundary(1);

              break; // break loop over all source eles
            }
            else
              matched = 0;

          }// if match
        }// loop over all seachbox elements
      }// enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////

      if(numproc>1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)vectofill->Length());

            // point coordinate
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            if(matched==1)
            {
              for(int dim=0;dim<(int)vectofill->Length();++dim)
              {
                data.AddtoPack((*vectofill)(dim));
              }
            }
          }

          // normal vector
          for(int dim=0;dim<globdim;++dim)
          {
            data.AddtoPack(normal_vec[dim]);
          }

          std::swap(sdata, data());

        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from,tag,rdata,length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched    = DRT::ParObject::ExtractInt(position,rdata);
        owner      = DRT::ParObject::ExtractInt(position,rdata);
        datalength = DRT::ParObject::ExtractInt(position,rdata);

        for(int i=0;i<globdim;++i)
          xvec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        if(matched == 1)
        {
          for(int dim=0;dim<datalength;++dim)
            (*vectofill)(dim) = DRT::ParObject::ExtractDouble(position,rdata);
        }

        for(int i=0;i<globdim;++i)
          normal_vec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      } // only if num procs > 1

      if(curr==(searchboxgeom.end()--) and irobin == (numproc-1) and matched == 0)
        dserror("could not match given point on any proc. Element(0)=%f",(*vectofill)(0));

    } // end for irobin

#ifdef DEBUG
    if(owner != comm.MyPID())
      dserror("after round robin every proc should have recieved his original data");
#endif

    // now every proc should store the requested quantities in vectofill for his targetxi
    // time to return that stuff
    for(int dim=0;dim<datalength;++dim)
      targetdata[dim]=(*vectofill)(dim);


    return err;
  };


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  template < DRT::Element::DiscretizationType sourcedistype,
  DRT::Element::DiscretizationType targetdistype >
  int InterpolateToBackgrdIntPoint(const Teuchos::RCP<DRT::Discretization> sourcedis,
      const Teuchos::RCP<DRT::Discretization> targetdis,
      DRT::Element& targetele,
      const std::vector<double>& targetxi,
      const std::vector<double>& targetedisp,
      const std::string action,
      std::vector<double>& targetdata,
      bool& match
  )
  {
    //////////////////////////////////////////////////////////////////////////////////
    ///////
    ///////  SOURCE  -->  TARGET
    ///////  $$ generalized version of InterpolateToImmersedIntPoint; search in discretization $$
    ///////  $$ instead of geometry; gets string action instead of int;                        $$
    ///////
    ///////  maps source values to target discretization.
    ///////
    ///////  comments on parallelism:
    ///////
    ///////  every proc needs to enter this routine to make communication work.
    ///////  in the extreme scenario where a proc owns no target elements at all
    ///////  no evaluate will be called and this proc doesn't enter this routine
    ///////  which means in turn, that a Barrier() would never be reached. in
    ///////  this case the attempt to communicate would definitely fail.
    ///////  so make sure, that each proc owns at least one target element and/or
    ///////  it is also to be assured by external measures, that each proc
    ///////  reenters this routine until the last element on every proc is
    ///////  evaluated (e.g. dummy call of Evaluate(...)).
    ///////
    /////////////////////////////////////////////////////////////////////////////////

    // error flag
    int err = 0;

    // get communicator
    const Epetra_Comm& comm = sourcedis->Comm();
    //std::cout<<"PROC "<<comm.MyPID()<<std::endl;
    // get the global problem
    DRT::Problem* problem = DRT::Problem::Instance();

    // dimension of global problem
    const int globdim = problem->NDim();

    // pointer to background element
    DRT::Element* sourceele;

    // declare vector of global current coordinates
    double x [globdim];

    // dimension of source discretization
    static const int source_dim = DRT::UTILS::DisTypeToDim<sourcedistype>::dim;

    // declarations and definitions for round robin loop
    const int numproc  = comm.NumProc();
    const int myrank   = comm.MyPID();                     // myrank
    const int torank   = (myrank + 1) % numproc;           // sends to
    const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    DRT::Exporter exporter(comm);

    Teuchos::RCP<const Epetra_Vector> dispnp = sourcedis->GetState("displacement");
    DRT::Element::LocationArray la(sourcedis->NumDofSets());

    // get current global coordinates of the given point xi of the target dis
    MORTAR::UTILS::LocalToCurrentGlobal<targetdistype>(targetele,globdim,&targetxi[0],targetedisp,&x[0]);
    std::vector<double> xvec(globdim);
    xvec[0]=x[0];
    xvec[1]=x[1];
    xvec[2]=x[2];

    // vector to fill by sourcedis->Evaluate (needs to be resized in calc class)
    Teuchos::RCP<Epetra_SerialDenseVector> vectofill = Teuchos::rcp(new Epetra_SerialDenseVector(4));
    (*vectofill)(0)=-1234.0;
    if(match) // was already matched, add no contribution
    {(*vectofill)(0)=0.0;
    (*vectofill)(1)=0.0;
    (*vectofill)(2)=0.0;}

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Epetra_SerialDenseVector> xi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(source_dim));

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////
    //////     round robin loop
    //////
    ////// in this loop every processor sends the current coordinates of the given target
    ////// point xi to his neighbors and checks whether the recieved coordinates lie within
    ////// one of the rank's searchbox elements.
    //////
    ////// if so, this neighbor interpolates the necessary quantities to this point.
    ////// the results, along with the information whether the point has already been
    ////// matched and with the id of the rank who requested the data originally is sent
    ////// around.
    //////
    ////// when the loop ends, all information are back on the proc who requested the data
    ////// originally and for every point the data should be stored in rdata.
    //////
    //////
    ////// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    //////
    //////////////////////////////////////////////////////////////////////////////////////////////
    // given point already matched?
    int matched = 0;
    if(match == true)
      matched = 1;

    // owner of given point
    int owner = myrank;
    // length of vectofill
    int datalength = (int)vectofill->Length();

    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      //////////////////////////////////////////////////////////////////////
      // loop over all source elements
      //
      // determine source element in which the given point lies.
      //
      // if the current proc can already match the given point xi (xvec)
      //
      /////////////////////////////////////////////////////////////////////
      const Epetra_Map* source_elecolmap = sourcedis->ElementColMap();
      int mynumcolele = source_elecolmap->NumMyElements();
      LINALG::Matrix<source_dim,1> xi(true);
      if(matched == 0)
      {
        for (int colele=0; colele < mynumcolele; ++colele)
        {
          bool converged = false;
          sourceele = sourcedis->gElement(source_elecolmap->GID(colele));

          // get parameter space coords xi in source element
          sourceele->LocationVector(*sourcedis,la,false);
          std::vector<double> mydispnp(la[0].lm_.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,la[0].lm_);

          std::vector<double> targeteleposition1(3);
          targeteleposition1[0]=(targetele.Nodes())[1]->X()[0];
          targeteleposition1[1]=(targetele.Nodes())[1]->X()[1];
          targeteleposition1[2]=(targetele.Nodes())[1]->X()[2];

          std::vector<double> targeteleposition7(3);
          targeteleposition7[0]=(targetele.Nodes())[7]->X()[0];
          targeteleposition7[1]=(targetele.Nodes())[7]->X()[1];
          targeteleposition7[2]=(targetele.Nodes())[7]->X()[2];

          double diagonaltargetelelength = sqrt(pow(targeteleposition1[0]-targeteleposition7[0],2)+pow(targeteleposition1[1]-targeteleposition7[1],2)+pow(targeteleposition1[2]-targeteleposition7[2],2));

          std::vector<double> sourceeleposition1(3);
          sourceeleposition1[0]=(sourceele->Nodes())[1]->X()[0] + mydispnp[1*3+0];
          sourceeleposition1[1]=(sourceele->Nodes())[1]->X()[1] + mydispnp[1*3+1];
          sourceeleposition1[2]=(sourceele->Nodes())[1]->X()[2] + mydispnp[1*3+2];

          // distance between given point of target dis and sourceelenode 1
          double dist = sqrt(pow(xvec[0]-sourceeleposition1[0],2)+pow(xvec[1]-sourceeleposition1[1],2)+pow(xvec[2]-sourceeleposition1[2],2));

          // only try to match given target point and sourceele if sourceele is in vicinity of target point
          if(dist < 1.6*diagonaltargetelelength)
          {
            MORTAR::UTILS::GlobalToCurrentLocal<sourcedistype>(*sourceele,&mydispnp[0],&xvec[0],&xi(0),converged);
            if(converged==false)
            {
              xi(0) = 2.0; xi(1) = 2.0; xi(2) = 2.0;
            }
          }
          else
          {
            xi(0) = 2.0; xi(1) = 2.0; xi(2) = 2.0;
          }

# ifdef DEBUG
          if(!converged)
          { std::vector<double> sourceeleposition1(3);
          sourceeleposition1[0]=(sourceele->Nodes())[1]->X()[0] + mydispnp[1*3+0];
          sourceeleposition1[1]=(sourceele->Nodes())[1]->X()[1] + mydispnp[1*3+1];;
          sourceeleposition1[2]=(sourceele->Nodes())[1]->X()[2] + mydispnp[1*3+2];;

          std::vector<double> sourceeleposition7(3);
          sourceeleposition7[0]=(sourceele->Nodes())[7]->X()[0];
          sourceeleposition7[1]=(sourceele->Nodes())[7]->X()[1];
          sourceeleposition7[2]=(sourceele->Nodes())[7]->X()[2];

          double diagonalsourcelength = sqrt(pow(sourceeleposition1[0]-sourceeleposition7[0],2)+pow(sourceeleposition1[1]-sourceeleposition7[1],2)+pow(sourceeleposition1[2]-sourceeleposition7[2],2));
          double distance = sqrt(pow(xvec[0]-sourceeleposition7[0],2)+pow(xvec[1]-sourceeleposition7[1],2)+pow(xvec[2]-sourceeleposition7[2],2));

          //              std::cout<<"UNCONVERGED MAPPING : (maybe not critical)"<<std::endl;
          //              std::cout<<"Global coordinate of target point : "<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<std::endl;
          //              std::cout<<"Global coordinate of source node1 : "<<sourceeleposition1[0]<<" "<<sourceeleposition1[1]<<" "<<sourceeleposition1[2]<<std::endl;
          //              std::cout<<"Diagonal length of source element : "<<diagonalsourcelength<<std::endl;
          //              std::cout<<"Distance between target point and source node 7 : "<<distance<<std::endl;

          if(distance < 5*diagonalsourcelength)
            dserror("MAPPING FROM SOURCE POINT TO TARGET DID NOT CONVERGE");
          else
          {
            xi(0)=-1234.0; xi(1)=-1234.0; xi(2)=-1234.0;
          }
          }
# endif
          /////////////////////////////////////////////////
          ////
          ////     MATCH
          ////
          /////////////////////////////////////////////////
          double tol = 1e-11;
          if ((abs(xi(0))-1.0)<tol and (abs(xi(1))-1.0)<tol and (abs(xi(2))-1.0)<tol)
          { // given point lies in element curr
            matched = 1;
            match = true;

            // parameter list
            Teuchos::ParameterList params;
            params.set<std::string>("action", action);
            // fill locationarray
            DRT::Element::LocationArray la(1);
            sourceele->LocationVector(*sourcedis,la,false);

            (*xi_dense)(0)=xi(0);
            (*xi_dense)(1)=xi(1);
            (*xi_dense)(2)=xi(2);

            Epetra_SerialDenseMatrix dummy1;
            Epetra_SerialDenseMatrix dummy2;
            Epetra_SerialDenseVector dummy3;

            // evaluate and perform action handed in to this function
            sourceele->Evaluate(params,*sourcedis,la,dummy1,dummy2,*vectofill ,*xi_dense,dummy3);

            break; // break loop over all source elements
          }// if match

          if(colele==(mynumcolele-1) and matched == 0 and irobin == (numproc-1)) // this point lies not within discretization at all
          {(*vectofill)(0)=0.0;(*vectofill)(1)=0.0;(*vectofill)(2)=0.0;(*vectofill)(3)=0.0;}
        }// loop over all seachbox elements
      }// enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////

      if(numproc>1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)vectofill->Length());

            // point coordinate
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            if(matched==1)
            {
              for(int dim=0;dim<(int)vectofill->Length();++dim)
              {
                data.AddtoPack((*vectofill)(dim));
              }
            }
          }
          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from,tag,rdata,length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched    = DRT::ParObject::ExtractInt(position,rdata);
        owner      = DRT::ParObject::ExtractInt(position,rdata);
        datalength = DRT::ParObject::ExtractInt(position,rdata);

        for(int i=0;i<globdim;++i)
          xvec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        if(matched == 1)
        {
          for(int dim=0;dim<datalength;++dim)
            (*vectofill)(dim) = DRT::ParObject::ExtractDouble(position,rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      } // only if numproc > 1

    } // end for irobin

#ifdef DEBUG
    if(owner != comm.MyPID())
      dserror("after round robin every proc should have recieved his original data");
#endif

    // now every proc should store the requested quantities in vectofill for his targetxi
    // time to return that stuff
    for(int dim=0;dim<datalength;++dim)
      targetdata[dim]=(*vectofill)(dim);


    return err;
  };


  /*----------------------------------------------------------------------*/
  /* Analogous to InterpolateToImmersedIntPoint() but no evaluate is      */
  /* called.                                                              */
  /*----------------------------------------------------------------------*/
  template < DRT::Element::DiscretizationType sourcedistype,
  DRT::Element::DiscretizationType targetdistype >
  int DetermineBackgrdElementsWithImmersedBoundary(const Teuchos::RCP<DRT::Discretization> sourcedis,
      const Teuchos::RCP<DRT::Discretization> targetdis,
      DRT::Element& targetele,
      const std::vector<double>& targetxi,
      const std::vector<double>& targetedisp
  )
  {
    // error flag
    int err = 0;

    // get communicator
    const Epetra_Comm& comm = sourcedis->Comm();

    // get the global problem
    DRT::Problem* problem = DRT::Problem::Instance();

    // dimension of global problem
    const int globdim = problem->NDim();

    // pointer to background element
    DRT::Element* sourceele;

    // declare vector of global current coordinates
    double x [globdim];

    // dimension of source discretization
    static const int source_dim = DRT::UTILS::DisTypeToDim<sourcedistype>::dim;
    static const int target_dim = DRT::UTILS::DisTypeToDim<targetdistype>::dim;

    // declarations and definitions for round robin loop
    const int numproc  = comm.NumProc();
    const int myrank   = comm.MyPID();                     // myrank
    const int torank   = (myrank + 1) % numproc;           // sends to
    const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    DRT::Exporter exporter(comm);

    // get current global coordinates of the given point xi of the target dis
    MORTAR::UTILS::LocalToCurrentGlobal<targetdistype>(targetele,globdim,&targetxi[0],targetedisp,&x[0]);
    std::vector<double> xvec(globdim);
    xvec[0]=x[0];
    xvec[1]=x[1];
    xvec[2]=x[2];

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Epetra_SerialDenseVector> xi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(source_dim));
    Teuchos::RCP<Epetra_SerialDenseVector> targetxi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(target_dim));

    Teuchos::RCP<Epetra_SerialDenseVector> normal_at_targetpoint = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    std::vector<double> normal_vec(3);// normal_vec[0]=0.0;normal_vec[1]=0.0; normal_vec[2]=0.0;

    Epetra_SerialDenseMatrix dummy1;
    Epetra_SerialDenseMatrix dummy2;
    Epetra_SerialDenseVector dummy3;

    Teuchos::ParameterList params;

    params.set<std::string>("action", "calc_cur_normal_at_point");
    for(int i=0; i<target_dim;++i)
      (*targetxi_dense)(i) = targetxi[i];

    DRT::Element::LocationArray targetla(1);
    targetele.LocationVector(*targetdis,targetla,false);

    targetele.Evaluate(params,*targetdis,targetla,dummy1,dummy2,*normal_at_targetpoint,*targetxi_dense,dummy3);
    normal_at_targetpoint->Scale(1.0/(normal_at_targetpoint->Norm2()));
    normal_vec[0] = (*normal_at_targetpoint)(0);
    normal_vec[1] = (*normal_at_targetpoint)(1);
    normal_vec[2] = (*normal_at_targetpoint)(2);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////
    //////     round robin loop
    //////
    //////////////////////////////////////////////////////////////////////////////////////////////
    // given point already matched?
    int matched = 0;
    // owner of given point
    int owner = myrank;

    // get possible elements being intersected by immersed structure
    DRT::Condition* searchbox = sourcedis->GetCondition("ImmersedSearchbox");
    std::map<int,Teuchos::RCP<DRT::Element> >& searchboxgeom = searchbox->Geometry();

    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      //////////////////////////////////////////////////////////////////////
      // loop over all background elements
      //
      // determine background element in which the given point lies.
      //
      // if the current proc can already match the given point xi (xvec)
      //
      /////////////////////////////////////////////////////////////////////
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
      LINALG::Matrix<source_dim,1> xi(true);
      if(!matched)
      {
        // every proc that has searchboxgeom elements
        for (curr=searchboxgeom.begin(); curr!=searchboxgeom.end(); ++curr)
        {
          bool converged = false;

          // get parameter space coords xi in source element
          MORTAR::UTILS::GlobalToLocal<sourcedistype>(*(curr->second),&xvec[0],&xi(0),converged);

//# ifdef DEBUG
//          if(!converged)
//          {
//            dserror("MAPPING FROM SOURCE POINT TO TARGET DID NOT CONVERGE");
//          }
//# endif
          /////////////////////////////////////////////////
          ////
          ////     MATCH
          ////
          /////////////////////////////////////////////////
          double tol = 1e-12;
          if ((abs(xi(0))-1.0)<tol and (abs(xi(1))-1.0)<tol and (abs(xi(2))-1.0)<tol)
          { // given point lies in element curr

            bool validsource = true;
            double scalarproduct= 0.0;
            Teuchos::RCP<Epetra_SerialDenseVector> vector;

            sourceele = curr->second.getRawPtr();

            // fill locationarray
            DRT::Element::LocationArray sourcela(1);
            sourceele->LocationVector(*sourcedis,sourcela,false);

            (*xi_dense)(0)=xi(0);
            (*xi_dense)(1)=xi(1);
            (*xi_dense)(2)=xi(2);

            // does targetpoint lie on an element edge?
            // only the sourceele located outward from the boundary should interpolate
            if( (abs(xi(0))<(1.0+tol) and abs(xi(0))>(1.0-tol)) or (abs(xi(1))<(1.0+tol) and abs(xi(1))>(1.0-tol)) or (abs(xi(2))<(1.0+tol) and abs(xi(2))>(1.0-tol)))
            {
              //              // calc sourcele diagonal length
              //              const double* X0 = sourceele->Nodes()[0]->X();
              //              const double* X6 = sourceele->Nodes()[1]->X();
              //              double srcelediaglength = sqrt(pow(X0[0]-X6[0],2)+pow(X0[1]-X6[1],2)+pow(X0[2]-X6[2],2));

              // get nodal coords of source ele
              vector = Teuchos::rcp(new Epetra_SerialDenseVector(3));
              std::vector<double> distances(sourceele->NumNode());
              for(int sourcenode=0; sourcenode<sourceele->NumNode(); sourcenode++)
              {
                const double* X;
                X = sourceele->Nodes()[sourcenode]->X();

                distances[sourcenode] = sqrt(pow(xvec[0]-X[0],2)+pow(xvec[1]-X[1],2)+pow(xvec[2]-X[2],2));

              }

              // get max. distance and corresponding index of src node
              int maxdistanceindex = -1;
              double maxvalue = 0.0;
              for(int i=0; i<(int)distances.size();++i)
              {
                if(distances[i]>maxvalue)
                {
                  maxvalue = distances[i];
                  maxdistanceindex = i;
                }
              }

              const double* X;
              X = sourceele->Nodes()[maxdistanceindex]->X();
              // build vector directed from targetele point (xvec) to sourceel node (X)
              // xvec should not lie at same ele surface as sourcenode
              (*vector)(0) = X[0]-xvec[0];
              (*vector)(1) = X[1]-xvec[1];
              (*vector)(2) = X[2]-xvec[2];
              vector->Scale(vector->Norm2());

              // build scalar product between normal and vector
              scalarproduct = (*vector)(0) * normal_vec[0] + (*vector)(1) * normal_vec[1] + (*vector)(2) * normal_vec[2];

              if(scalarproduct > 0.0)
                validsource = true;
              else
                validsource = false;

#ifdef DEBUG
              if(abs(scalarproduct) < 1e-13 and validsource)
              {
                std::cout<<"scalarproduct = "<<scalarproduct<<std::endl;
                dserror("normal vec and vec from point to node perpendicular.\n"
                        "no robust determination of outward lying src ele possible");
              }
              if(abs(scalarproduct) < -1e-12 and validsource)
              {
                std::cout<<"scalarproduct = "<<scalarproduct<<std::endl;
                dserror("normal vec and vec from point to node have no sharp angle.\n"
                        "no robust determination of outward laying src ele possible");
              }
#endif

              //count++;
            } // if on element edge

            if(validsource)
            {
              matched = 1;

              Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::FluidImmersed>(curr->second)->SetIsImmersedBoundary(1);

              DRT::Node** nodes = curr->second->Nodes();
              for(int node=0;node<curr->second->NumNode();++node)
                static_cast<IMMERSED::ImmersedNode* >(nodes[node])->SetIsImmersedBoundary(1);

              break; // break loop over all source eles
            }
            else
              matched = 0;

          }// if match
        }// loop over all seachbox elements
      }// enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////


      if(numproc>1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);

            // point coordinate
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            // normal vector
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(normal_vec[dim]);
            }

            std::swap(sdata, data());
          }
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from,tag,rdata,length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched    = DRT::ParObject::ExtractInt(position,rdata);
        owner      = DRT::ParObject::ExtractInt(position,rdata);

        for(int i=0;i<globdim;++i)
          xvec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        for(int i=0;i<globdim;++i)
          normal_vec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      } // only if parallel simulation (numproc > 1)
    } // end for irobin

#ifdef DEBUG
    if(owner != comm.MyPID())
      dserror("after round robin every proc should have recieved his original data");
#endif

    return err;
  };


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  template < DRT::Element::DiscretizationType sourcedistype,
  DRT::Element::DiscretizationType targetdistype >
  int InterpolateToBackgrdIntPointFAST(
      std::map<int,std::set<int> >& curr_subset_of_structdis,
      const Teuchos::RCP<DRT::Discretization> sourcedis,
      const Teuchos::RCP<DRT::Discretization> targetdis,
      DRT::Element& targetele,
      const std::vector<double>& targetxi,
      const std::vector<double>& targetedisp,
      const std::string action,
      std::vector<double>& targetdata,
      bool& match,
      bool doCommunication = true
  )
  {
    //////////////////////////////////////////////////////////////////////////////////
    ///////
    /////// FAST version of InterpolateToBackgrdIntPoint
    ///////
    /////// employing search tree for immersed discretization.
    ///////
    ///////
    /////////////////////////////////////////////////////////////////////////////////

    // error flag
    int err = 0;

    // get communicator
    const Epetra_Comm& comm = sourcedis->Comm();

    // get the global problem
    DRT::Problem* problem = DRT::Problem::Instance();

    // dimension of global problem
    const int globdim = problem->NDim();

    // pointer to background element
    DRT::Element* sourceele;

    // declare vector of global current coordinates
    double x [globdim];

    // dimension of source discretization
    static const int source_dim = DRT::UTILS::DisTypeToDim<sourcedistype>::dim;

    // declarations and definitions for round robin loop
          int numproc  = comm.NumProc();
    const int myrank   = comm.MyPID();                     // myrank
    const int torank   = (myrank + 1) % numproc;           // sends to
    const int fromrank = (myrank + numproc - 1) % numproc; // recieves from

    if(numproc == 1)
      doCommunication=false;

    if(doCommunication==false)
      numproc=1;

    Teuchos::RCP<const Epetra_Vector> dispnp = sourcedis->GetState("displacement");
    DRT::Element::LocationArray la(sourcedis->NumDofSets());

    // get current global coordinates of the given point xi of the target dis
    MORTAR::UTILS::LocalToCurrentGlobal<targetdistype>(targetele,globdim,&targetxi[0],targetedisp,&x[0]);
    std::vector<double> xvec(globdim);
    xvec[0]=x[0];
    xvec[1]=x[1];
    xvec[2]=x[2];

    // vector to be filled in sourcedis->Evaluate (needs to be resized in calc class)
    Teuchos::RCP<Epetra_SerialDenseVector> vectofill = Teuchos::rcp(new Epetra_SerialDenseVector(4));
    (*vectofill)(0)=-1234.0;
    if(match) // was already matched, add no contribution
    {(*vectofill)(0)=0.0;
     (*vectofill)(1)=0.0;
     (*vectofill)(2)=0.0;
     (*vectofill)(3)=0.0; }

    // save parameter space coordinate in serial dense vector
    Teuchos::RCP<Epetra_SerialDenseVector> xi_dense = Teuchos::rcp(new Epetra_SerialDenseVector(source_dim));

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////
    //////     round robin loop
    //////
    ////// in this loop every processor sends the current coordinates of the given target
    ////// point xi to his neighbors and checks whether the recieved coordinates lie within
    ////// one of the rank's searchbox elements.
    //////
    ////// if so, this neighbor interpolates the necessary quantities to this point.
    ////// the results, along with the information whether the point has already been
    ////// matched and with the id of the rank who requested the data originally is sent
    ////// around.
    //////
    ////// when the loop ends, all information are back on the proc who requested the data
    ////// originally and for every point the data should be stored in rdata.
    //////
    //////
    ////// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    //////
    //////////////////////////////////////////////////////////////////////////////////////////////
    // given point already matched?
    int matched = 0;
    if(match == true)
      matched = 1;

    // owner of given point
    int owner = myrank;
    // length of vectofill
    int datalength = (int)vectofill->Length();

    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      //////////////////////////////////////////////////////////////////////
      // loop over all source elements
      //
      // determine source element in which the given point lies.
      //
      // if the current proc can already match the given point xi (xvec)
      //
      /////////////////////////////////////////////////////////////////////
      LINALG::Matrix<source_dim,1> xi(true);
      if(matched == 0)
      {
        for(std::map<int, std::set<int> >::const_iterator closele = curr_subset_of_structdis.begin(); closele != curr_subset_of_structdis.end(); closele++)
        {
          for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
          {
            //std::cout<<"tgtele "<<targetele.Id()<<": "<<*eleIter<<std::endl;
            bool converged = false;
            sourceele = sourcedis->gElement(*eleIter);

            // get parameter space coords xi in source element
            sourceele->LocationVector(*sourcedis,la,false);
            std::vector<double> mydispnp(la[0].lm_.size());
            DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,la[0].lm_);

            std::vector<double> targeteleposition1(3);
            targeteleposition1[0]=(targetele.Nodes())[1]->X()[0];
            targeteleposition1[1]=(targetele.Nodes())[1]->X()[1];
            targeteleposition1[2]=(targetele.Nodes())[1]->X()[2];

            std::vector<double> targeteleposition7(3);
            targeteleposition7[0]=(targetele.Nodes())[7]->X()[0];
            targeteleposition7[1]=(targetele.Nodes())[7]->X()[1];
            targeteleposition7[2]=(targetele.Nodes())[7]->X()[2];

            double diagonaltargetelelength = sqrt(pow(targeteleposition1[0]-targeteleposition7[0],2)+pow(targeteleposition1[1]-targeteleposition7[1],2)+pow(targeteleposition1[2]-targeteleposition7[2],2));

            std::vector<double> sourceeleposition1(3);
            sourceeleposition1[0]=(sourceele->Nodes())[1]->X()[0] + mydispnp[1*3+0];
            sourceeleposition1[1]=(sourceele->Nodes())[1]->X()[1] + mydispnp[1*3+1];
            sourceeleposition1[2]=(sourceele->Nodes())[1]->X()[2] + mydispnp[1*3+2];

            // distance between given point of target dis and sourceelenode 1
            double dist = sqrt(pow(xvec[0]-sourceeleposition1[0],2)+pow(xvec[1]-sourceeleposition1[1],2)+pow(xvec[2]-sourceeleposition1[2],2));

            // only try to match given target point and sourceele if sourceele is in vicinity of target point
            if(dist < 1.6*diagonaltargetelelength)
            {
              MORTAR::UTILS::GlobalToCurrentLocal<sourcedistype>(*sourceele,&mydispnp[0],&xvec[0],&xi(0),converged);
              if(converged==false)
              {
                xi(0) = 2.0; xi(1) = 2.0; xi(2) = 2.0;
              }
            }
            else
            {
              xi(0) = 2.0; xi(1) = 2.0; xi(2) = 2.0;
            }

//# ifdef DEBUG
//            if(!converged and dist < 1.6*diagonaltargetelelength)
//            { std::vector<double> sourceeleposition1(3);
//            sourceeleposition1[0]=(sourceele->Nodes())[1]->X()[0] + mydispnp[1*3+0];
//            sourceeleposition1[1]=(sourceele->Nodes())[1]->X()[1] + mydispnp[1*3+1];;
//            sourceeleposition1[2]=(sourceele->Nodes())[1]->X()[2] + mydispnp[1*3+2];;
//
//            std::vector<double> sourceeleposition7(3);
//            sourceeleposition7[0]=(sourceele->Nodes())[7]->X()[0];
//            sourceeleposition7[1]=(sourceele->Nodes())[7]->X()[1];
//            sourceeleposition7[2]=(sourceele->Nodes())[7]->X()[2];
//
//            //double diagonalsourcelength = sqrt(pow(sourceeleposition1[0]-sourceeleposition7[0],2)+pow(sourceeleposition1[1]-sourceeleposition7[1],2)+pow(sourceeleposition1[2]-sourceeleposition7[2],2));
//            //double distance = sqrt(pow(xvec[0]-sourceeleposition7[0],2)+pow(xvec[1]-sourceeleposition7[1],2)+pow(xvec[2]-sourceeleposition7[2],2));
//
//            //              std::cout<<"UNCONVERGED MAPPING : (maybe not critical)"<<std::endl;
//            //              std::cout<<"Global coordinate of target point : "<<xvec[0]<<" "<<xvec[1]<<" "<<xvec[2]<<std::endl;
//            //              std::cout<<"Global coordinate of source node1 : "<<sourceeleposition1[0]<<" "<<sourceeleposition1[1]<<" "<<sourceeleposition1[2]<<std::endl;
//            //              std::cout<<"Diagonal length of source element : "<<diagonalsourcelength<<std::endl;
//            //              std::cout<<"Distance between target point and source node 7 : "<<distance<<std::endl;
//
//            dserror("MAPPING FROM SOURCE POINT TO TARGET DID NOT CONVERGE");
//            xi(0)=-1234.0; xi(1)=-1234.0; xi(2)=-1234.0;
//
//            }
//# endif
            /////////////////////////////////////////////////
            ////
            ////     MATCH
            ////
            /////////////////////////////////////////////////
            double tol = 1e-11;
            if ((abs(xi(0))-1.0)<tol and (abs(xi(1))-1.0)<tol and (abs(xi(2))-1.0)<tol)
            { // given point lies in element curr
              matched = 1;
              match = true;

              // parameter list
              Teuchos::ParameterList params;
              params.set<std::string>("action", action);
              // fill locationarray
              DRT::Element::LocationArray la(1);
              sourceele->LocationVector(*sourcedis,la,false);

              (*xi_dense)(0)=xi(0);
              (*xi_dense)(1)=xi(1);
              (*xi_dense)(2)=xi(2);

              Epetra_SerialDenseMatrix dummy1;
              Epetra_SerialDenseMatrix dummy2;
              Epetra_SerialDenseVector dummy3;

              // evaluate and perform action handed in to this function
              sourceele->Evaluate(params,*sourcedis,la,dummy1,dummy2,*vectofill ,*xi_dense,dummy3);

              break; // break loop over all source elements
            }// if match

            if(closele == (curr_subset_of_structdis.end()--) and eleIter == ((closele->second).end()--) and matched == 0 and irobin == (numproc-1)) // this point lies not within discretization at all
            {(*vectofill)(0)=0.0;(*vectofill)(1)=0.0;(*vectofill)(2)=0.0;(*vectofill)(3)=0.0;}
          }// loop over all seachbox elements
        }// loop over closeele set
      }// enter loop only if point not yet matched
      /////////////////////////////////////////////////////////////////////

      if(doCommunication)
      {
        DRT::Exporter exporter(comm);

        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(matched);
            data.AddtoPack(owner);
            data.AddtoPack((int)vectofill->Length());

            // point coordinate
            for(int dim=0;dim<globdim;++dim)
            {
              data.AddtoPack(xvec[dim]);
            }

            if(matched==1)
            {
              for(int dim=0;dim<(int)vectofill->Length();++dim)
              {
                data.AddtoPack((*vectofill)(dim));
              }
            }
          }
          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from,tag,rdata,length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        matched    = DRT::ParObject::ExtractInt(position,rdata);
        owner      = DRT::ParObject::ExtractInt(position,rdata);
        datalength = DRT::ParObject::ExtractInt(position,rdata);

        for(int i=0;i<globdim;++i)
          xvec[i] = DRT::ParObject::ExtractDouble(position,rdata);

        if(matched == 1)
        {
          for(int dim=0;dim<datalength;++dim)
            (*vectofill)(dim) = DRT::ParObject::ExtractDouble(position,rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        comm.Barrier();
      } // only if numproc > 1

    } // end for irobin

#ifdef DEBUG
    if(owner != comm.MyPID())
      dserror("after round robin every proc should have recieved his original data");
#endif

    // now every proc should store the requested quantities in vectofill for his targetxi
    // time to return that stuff
    for(int dim=0;dim<datalength;++dim)
      targetdata[dim]=(*vectofill)(dim);


    return err;
  };


}// namsepace IMMERSED

#endif

