/*!----------------------------------------------------------------------
\file turbulence_statistics_bfda.cpp

\level 2
<pre>

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
<pre>
*----------------------------------------------------------------------*/
/*--created in a Semesterarbeit by Sebastian Buechner 2014 ------------*/
/*---------------------------------------------------------------------*/
/*                                                                     */
/*  ___________          ______________                                */
/* |           \ ______|              |                                */
/* |            _______O              |     ----------> z direction    */
/* |___________/       |______________|                                */
/*                                                                     */
/* (Origin of coordinate system is marked by O)                        */
/*                                                                     */
/*                                                                     */
/*---------------------------------------------------------------------*/
/*
*/



#include "turbulence_statistics_bfda.H"





//#define COMBINE_SAMPLES

/*----------------------------------------------------------------------*/
/*!
  \brief Standard Constructor (public)

  <pre>
  o Create vector "zcoordinates" with coordinates along z-axis
  o Create matrix "rcoordinates_" with coordinates of radial coordinates of evaluation planes
    columns are evaluation planes corresponding to the positions in "posEvaluation_"
  </pre>

*/
/*----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsBfda::TurbulenceStatisticsBfda(
  Teuchos::RCP<DRT::Discretization> actdis,
  Teuchos::ParameterList&           params,
  const std::string&                statistics_outfilename):
  discret_      (actdis),
  params_       (params),
  statistics_outfilename_(statistics_outfilename)
{
  if (discret_->Comm().MyPID()==0)
  {
    std::cout << "This is the turbulence statistics manager of a benchmark of blood flowing through a channel:" << std::endl;
  }

  /*---------------------------------------------------------------------*/
  /*---------------------------------------------------------------------*/
  // Parameters:
  // z-coordinates for evaluation
  posEvaluation_.clear();
  posEvaluation_.push_back(  -0.088  );
  posEvaluation_.push_back(  -0.064  );
  posEvaluation_.push_back(  -0.048  );
  posEvaluation_.push_back(  -0.02   );
  posEvaluation_.push_back(  -0.008  );
  posEvaluation_.push_back(   0.0    );
  posEvaluation_.push_back(   0.008  );
  posEvaluation_.push_back(   0.016  );
  posEvaluation_.push_back(   0.024  );
  posEvaluation_.push_back(   0.032  );
  posEvaluation_.push_back(   0.06   );
  posEvaluation_.push_back(   0.08   );
  /*---------------------------------------------------------------------*/
  /*---------------------------------------------------------------------*/

  // number of evaluation planes
  int numPosEvaluation = posEvaluation_.size();

  // Vector for evaluation planes of actual mesh (nodes not exactly on evaluation planes
  actPosEvaluation_.clear();

  // Resize matrix for radial coordinates in a first step (final size is later set)
  rcoordinates_.reshape(1,numPosEvaluation);


  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim!=3)
    dserror("Evaluation of turbulence statistics only for 3d flow problems!");

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  togglew_      = LINALG::CreateVector(*dofrowmap,true);
  togglep_      = LINALG::CreateVector(*dofrowmap,true);


  //----------------------------------------------------------------------
  // create sets of coordinates
  //----------------------------------------------------------------------
  // the criterion allows differences in coordinates by 1e-9
  std::set<double,LineSortCriterion> zavcoords;
  std::set<double,LineSortCriterion> ravcoords;


  // loop nodes and build sets of lines in z-direction
  // accessible on this proc
  for (int i=0; i<discret_->NumMyRowNodes(); ++i)
  {
    // create node object
    DRT::Node* node = discret_->lRowNode(i);

    // Is the actual node on z-axis?
    // => get z-coordinates of nodes on z-axis
    if (node->X()[1]<2e-9 && node->X()[1]>-2e-9
        &&node->X()[0]<2e-9 && node->X()[0]>-2e-9)
      //then insert z-coordinate to vector "zavcoords"
      zavcoords.insert(node->X()[2]);
  }

  //--------------------------------------------------------------------
  // round robin loop to communicate coordinates to all procs
  //--------------------------------------------------------------------
//  {
#ifdef PARALLEL
    int myrank  =discret_->Comm().MyPID();
#endif
    int numprocs=discret_->Comm().NumProc();
//
    std::vector<char> sblock;
    std::vector<char> rblock;
//
#ifdef PARALLEL
    // create an exporter for point to point communication
    DRT::Exporter exporter(discret_->Comm());
#endif


    // first, communicate coordinates in x1-direction
    for (int np=0;np<numprocs;++np)
    {
            DRT::PackBuffer data;

            for (std::set<double,LineSortCriterion>::iterator zline=zavcoords.begin();
                 zline!=zavcoords.end();
                 ++zline)
            {
              DRT::ParObject::AddtoPack(data,*zline);
            }
            data.StartPacking();
            for (std::set<double,LineSortCriterion>::iterator zline=zavcoords.begin();
                 zline!=zavcoords.end();
                 ++zline)
            {
              DRT::ParObject::AddtoPack(data,*zline);
            }
            std::swap( sblock, data() );

      #ifdef PARALLEL
            MPI_Request request;
            int         tag    =myrank;

            int         frompid=myrank;
            int         topid  =(myrank+1)%numprocs;

            int         length=sblock.size();

            exporter.ISend(frompid,topid,
                           &(sblock[0]),sblock.size(),
                           tag,request);

            rblock.clear();

            // receive from predecessor
            frompid=(myrank+numprocs-1)%numprocs;
            exporter.ReceiveAny(frompid,tag,rblock,length);

            if(tag!=(myrank+numprocs-1)%numprocs)
            {
              dserror("received wrong message (ReceiveAny)");
            }

            exporter.Wait(request);

            {
              // for safety
              exporter.Comm().Barrier();
            }
      #else
            // dummy communication
            rblock.clear();
            rblock=sblock;
      #endif

            //--------------------------------------------------
            // Unpack received block into set of all planes.
            {
              std::vector<double> coordsvec;

              coordsvec.clear();

              std::vector<char>::size_type index = 0;
              while (index < rblock.size())
              {
                double onecoord;
                DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
                zavcoords.insert(onecoord);
              }
            }
    }





  //----------------------------------------------------------------------
  // push coordinates in vectors (copy zavcoords to a new vector zcoordinates)
  //----------------------------------------------------------------------
//  {
    zcoordinates_ = Teuchos::rcp(new std::vector<double> );

    for(std::set<double,LineSortCriterion>::iterator coord1=zavcoords.begin();
        coord1!=zavcoords.end();
        ++coord1)
    {
      zcoordinates_->push_back(*coord1);
    }
//  }



  //----------------------------------------------------------------------
  // number of coordinates in z-direction
  //----------------------------------------------------------------------
  numzcoor_ = zcoordinates_->size();



  //----------------------------------------------------------------------
  // Error Message if no nodes are on z-axis (nothing for evaluation)
  //----------------------------------------------------------------------
  if (discret_->Comm().MyPID()==0)
  {
    if (numzcoor_==0)
    {
      std::cout << std::endl << "WARNING!" << std::endl;
      std::cout << "No nodes are on the centerline! Please check if center of cells is on the z-axis so that there are no nodes on it!" << std::endl << std::endl;
    }
    else
      std::cout << "Number of nodes on z-axis: " << numzcoor_ << std::endl;
  }









  // loop nodes and build sets of lines in radial direction
  // not exact position of evaluation is evaluated because there could be no nodes
  // instead the next higher z-coordinate is used if the exact position is not available
  // accessible on this proc

    // Is the actual node on y-axis of an evaluation plane?
    // => get radial coordinates of nodes on evaluation planes
    // Loop through all evaluation planes
    for (int actPosEval=0; actPosEval<numPosEvaluation; actPosEval++)
    {
      // Reset variable for radial coordinates on each proc
      ravcoords.clear();


      int actRadNode = 0;

      // get exact or next z coordinate next to the evaluation plane where a node exists
      // therefore loop through the entries in zcoordinates_ until an entry is equal or bigger
      // than the desired position of the evaluation plane
      double actZ = 0;
      for (std::set<double,LineSortCriterion>::iterator temp_ptr_actZ=zavcoords.begin(); temp_ptr_actZ!=zavcoords.end(); ++temp_ptr_actZ)
      {
        if(*temp_ptr_actZ>=posEvaluation_[actPosEval])
        {
          actZ = *temp_ptr_actZ;
          actPosEvaluation_.push_back(*temp_ptr_actZ);
          break;
        }
      }

      // Loop to get all radial coordinates of this evaluation plane
      for (int i=0; i<discret_->NumMyRowNodes(); ++i)
      {
        // create node object
        DRT::Node* node = discret_->lRowNode(i);

        // Get all nodes on evaluation plane at this z position with y=0 and x>2e-9
        if (node->X()[2]<actZ+2e-9 && node->X()[2]>actZ-2e-9
            && node->X()[1]<2e-9 && node->X()[1]>-2e-9
            && node->X()[0]>2e-9)
          {
          //insert radial coordinate (here x) to matrix
          //rcoordinates_(actRadNode, actPosEval)=1.0; //(*node).X()[0]; ///////////////////////////////////////////////////////////////////////////////////////////
          ravcoords.insert((*node).X()[0]);
          actRadNode++;
          }
      }


      int countActRadNodeOnAllProcs = 0;
      discret_->Comm().SumAll(&actRadNode,&countActRadNodeOnAllProcs,1);


    #ifdef PARALLEL
        int myrank  =discret_->Comm().MyPID();
    #endif
        int numprocs=discret_->Comm().NumProc();

        std::vector<char> sblock;
        std::vector<char> rblock;

    #ifdef PARALLEL
        // create an exporter for point to point communication
        DRT::Exporter exporter(discret_->Comm());
    #endif


    // first, communicate coordinates in x1-direction
    for (int np=0;np<numprocs;++np)
    {
              DRT::PackBuffer data;

              for (std::set<double,LineSortCriterion>::iterator zline=ravcoords.begin();
                   zline!=ravcoords.end();
                   ++zline)
              {
                DRT::ParObject::AddtoPack(data,*zline);
              }
              data.StartPacking();
              for (std::set<double,LineSortCriterion>::iterator zline=ravcoords.begin();
                   zline!=ravcoords.end();
                   ++zline)
              {
                DRT::ParObject::AddtoPack(data,*zline);
              }
              std::swap( sblock, data() );

        #ifdef PARALLEL
              MPI_Request request;
              int         tag    =myrank;

              int         frompid=myrank;
              int         topid  =(myrank+1)%numprocs;

              int         length=sblock.size();

              exporter.ISend(frompid,topid,
                             &(sblock[0]),sblock.size(),
                             tag,request);

              rblock.clear();

              // receive from predecessor
              frompid=(myrank+numprocs-1)%numprocs;
              exporter.ReceiveAny(frompid,tag,rblock,length);

              if(tag!=(myrank+numprocs-1)%numprocs)
              {
                dserror("received wrong message (ReceiveAny)");
              }

              exporter.Wait(request);

              {
                // for safety
                exporter.Comm().Barrier();
              }
        #else
              // dummy communication
              rblock.clear();
              rblock=sblock;
        #endif

              //--------------------------------------------------
              // Unpack received block into set of all planes.
              {
                std::vector<double> coordsvec;

                coordsvec.clear();
//                ravcoords.clear();  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                std::vector<char>::size_type index = 0;
                while (index < rblock.size())
                {
                  double onecoord;
                  DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
                  ravcoords.insert(onecoord);
                }
              }
    }


    // Is actual number of radial coordinates higher than current maximum (= number of rows of matrix)
    // +1 at actRadNode because counting starts at 0 but numRows starts counting at 1
    if((countActRadNodeOnAllProcs+1)>rcoordinates_.numRows()) rcoordinates_.reshape(countActRadNodeOnAllProcs+1,numPosEvaluation);


      //----------------------------------------------------------------------
      // push ravcoords in matrix rcoordinates_
      //----------------------------------------------------------------------
        int count_coord1=1;
        for(std::set<double,LineSortCriterion>::iterator coord1=ravcoords.begin();
            coord1!=ravcoords.end();
            ++coord1)
        {
          rcoordinates_(count_coord1,actPosEval)=*coord1;
          count_coord1++;
        }
    }// End of Loop for evaluation planes

  // Get maximum number of radial points
  numrstatlocations_=rcoordinates_.numRows();

  //----------------------------------------------------------------------
  // Error Message if no nodes are in radial direction in a specific evaluation plane (nothing for evaluation)
  //----------------------------------------------------------------------
  if (discret_->Comm().MyPID()==0)
  {
      if (numrstatlocations_==0)
      {
        std::cout << std::endl << "WARNING!" << std::endl;
        std::cout << "No nodes are on the evaluation planes! Please check if center of cells is on the z-axis so that there are no nodes on it!" << std::endl << std::endl;
      }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of mean values
  //----------------------------------------------------------------------
    zsumw_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
    zsumw_->Reshape(1,numzcoor_);

    zsump_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
    zsump_->Reshape(1,numzcoor_);

    rsumw_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
    rsumw_->Reshape(numrstatlocations_,numPosEvaluation);

    rsump_ =  Teuchos::rcp(new Epetra_SerialDenseMatrix);
    rsump_->Reshape(numrstatlocations_,numPosEvaluation);


  // set number of samples to zero
    numsamp_ = 0;

  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file

  Teuchos::RCP<std::ofstream> log;
  Teuchos::RCP<std::ofstream> log_2;

  if (discret_->Comm().MyPID()==0)
  {
    std::string s(statistics_outfilename_);;
    std::string s2(statistics_outfilename_);

//    s.append(".flow_statistics");
    s.append(".flow_statistics_along_z");
    s2.append(".flow_statistics_evaluation_planes");

    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow of blood flowing through a channel\n\n";

    log_2 = Teuchos::rcp(new std::ofstream(s2.c_str(),std::ios::out));
    (*log_2) << "# Statistics for turbulent incompressible flow of blood flowing through a channel\n\n";

    log->flush();
    log_2->flush();
  }

  return;
}// TurbulenceStatisticsBfda::TurbulenceStatisticsBfda


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsBfda::~TurbulenceStatisticsBfda()
{
  return;
}// TurbulenceStatisticsBfda::~TurbulenceStatisticsBfda()


void FLD::TurbulenceStatisticsBfda::DoTimeSample(
Teuchos::RCP<Epetra_Vector> velnp
)
{
  //----------------------------------------------------------------------
  // increase sample counter
  //----------------------------------------------------------------------
  numsamp_++;

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // sampling of velocity/pressure values along z axis
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  int znodnum = -1;

  // Big loop along z-axis
  for (std::vector<double>::iterator zline=zcoordinates_->begin();
       zline!=zcoordinates_->end();
       ++zline)
  {
    znodnum++;

      // toggle vectors are one in the position of a dof of this node, else 0
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // count the number of nodes in x3-direction contributing to this nodal value
      int countnodes=0;

      //write a 1.0 at the position of the actual node of the processor in the toggle vectors
      for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
      {
        DRT::Node* node = discret_->lRowNode(nn);

        // If node is on z-axis then get toggle vector for pressure and velocity at actual node               // this is the wall node
        if ((node->X()[2]<(*zline+2e-9) and node->X()[2]>(*zline-2e-9))   and   (node->X()[0]<(2e-9) and node->X()[0]>(-2e-9))   and   (node->X()[1]<(2e-9) and node->X()[1]>(-2e-9)))
        {
          // Creation of vector dof that stores the values that variable "discret" points to
          // dof = [position of u     position of v     position of w     position of p]
          // position in terms of position of all degrees of freedom: u1 v1 w1 p1 u2 v2 w2 p2 u3 ...
          std::vector<int> dof = discret_->Dof(node);

          double           one = 1.0;
          togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));
          togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));


          // increase counter for evaluated nodes
          countnodes++;
        }
      }

      // Count nodes on all procs (should be 1)
      int countnodesonallprocs=1;
      discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);
      if (countnodesonallprocs)
      {
        //----------------------------------------------------------------------
        // get values for velocity derivative and pressure
        //----------------------------------------------------------------------
        double w;
        velnp->Dot(*togglew_,&w);
        double p;
        velnp->Dot(*togglep_,&p);

        //----------------------------------------------------------------------
        // add spatial mean values to statistical sample
        //----------------------------------------------------------------------
        int rnodnum = 0;
        (*zsumw_)(rnodnum,znodnum)+=w;
        (*zsump_)(rnodnum,znodnum)+=p;
      }
  }// End of big loop along z-axis
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------


  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // Get mean velocity from evaluation planes
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  // actual radial coordinate
  double actR=0.0;

  for (unsigned actPosEval=0; actPosEval<actPosEvaluation_.size(); actPosEval++)
  {
    for (int rnodnum=0;
        rnodnum<numrstatlocations_;
        rnodnum++)
    {
        // Get position
        actR = rcoordinates_(rnodnum, actPosEval);

        // Less nodes at the beginning -> only calculate until the last node and not until the matrix is finished (0 else)
        if((rcoordinates_(rnodnum, actPosEval)==0.0) and (rnodnum!=0))
        {
          (*rsumw_)(rnodnum,actPosEval)=-1.0;
          (*rsump_)(rnodnum,actPosEval)=-1.0;
          continue;
        }

        togglew_->PutScalar(0.0);
        togglep_->PutScalar(0.0);

        int countnodes=0;

        // get mean value of 4 nodes on a circle on the axis
        for (int circ=0; circ<4; circ++)
        {
          // Calculate coordinates of point
          double actX=sin(circ/2.0*3.141592)*actR;
          double actY=cos(circ/2.0*3.141592)*actR;

          // Because tolerances should not be added set nearly zero values to zero
          if (actX < 1e-8 and actX > -1e-8) actX=0.0;
          if (actY < 1e-8 and actY > -1e-8) actY=0.0;

          // Put 1.0 in toggle vectors if node on proc is desired node
          for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
          {
            DRT::Node* node = discret_->lRowNode(nn);

            // If node is on desired position then get toggle vector for pressure and velocity at actual node
            if ((node->X()[1]<(actY+2e-9) and node->X()[1]>(actY-2e-9))   and   (node->X()[0]<(actX+2e-9) and node->X()[0]>(actX-2e-9))
                and   (node->X()[2]<(actPosEvaluation_[actPosEval]+2e-9) and node->X()[2]>(actPosEvaluation_[actPosEval]-2e-9)))
            {
              std::vector<int> dof = discret_->Dof(node);
              double one = 1.0;
              togglep_->ReplaceGlobalValues(1,&one,&(dof[3]));
              togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
              countnodes++;
            }
          }
        }

        // Sum nodes on all procs (should be 4 nodes at a specific radius per evaluation plane except of r=0 where it is 1)
        int countnodesonallprocs=4;
        discret_->Comm().SumAll(&countnodes,&countnodesonallprocs,1);

        // At r=0 the loop above wrote 4 times a 1.0 at the same position of the toggle vector so this is a special case
        if(rcoordinates_(rnodnum, actPosEval)>-1e-9 and rcoordinates_(rnodnum, actPosEval)<1e-9) countnodesonallprocs=1;

        if (countnodesonallprocs)
        {
          // Get sum of values of 4 nodes at specific radius and evaluation plane
          double w;
          velnp->Dot(*togglew_,&w);
          double p;
          velnp->Dot(*togglep_,&p);

          // Calculate mean value
          double wsm=w/countnodesonallprocs;
          double psm=p/countnodesonallprocs;

          // Add mean value to matrix
          (*rsumw_)(rnodnum,actPosEval)+=wsm;
          (*rsump_)(rnodnum,actPosEval)+=psm;
        }
     }// End of big loop along radial coordinate
  }//End of loop through all evaluation planes
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------

  return;
}// TurbulenceStatisticsBfda::DoTimeSample


/*----------------------------------------------------------------------*
 *
 *----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsBfda::DumpStatistics(int step)
{
  //----------------------------------------------------------------------------------------------
  // file *.flow_statistics_along_z
  //----------------------------------------------------------------------------------------------
  Teuchos::RCP<std::ofstream> log;
  Teuchos::RCP<std::ofstream> log_2;
  if (discret_->Comm().MyPID()==0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics_along_z");    // s.append(".flow_statistics");
    log = Teuchos::rcp(new std::ofstream(s.c_str(),std::ios::out));

    // Output of mean velocity and pressure along z-axis
    (*log) << "Output file of FDA blood flow benchmark – evaluation of z-axis                \n";
    (*log) << "------------------------------------------------------------------------------\n";
    (*log) << "                                                                              \n";
    (*log) << "Number of samples: " << numsamp_ << "                                         \n";
    (*log) << "\n";
    (*log) << "Positions_of_nodes_on_z-axis Mean_values_of_velocity_w_in_z-direction Mean_values_of_pressure_p_in_z-direction \n";

    for (unsigned i=0; i<zcoordinates_->size(); ++i)
    {

        (*log) <<  " "  << std::setw(11) << std::setprecision(4) << (*zcoordinates_)[i];
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*zsumw_)(0,i);
        (*log) << "   " << std::setw(11) << std::setprecision(4) << (*zsump_)(0,i);
        (*log) << "\n";
    }


    //----------------------------------------------------------------------------------------------
    // file *.flow_statistics_evaluation_planes
    //----------------------------------------------------------------------------------------------
    std::string s2(statistics_outfilename_);
    s2.append(".flow_statistics_evaluation_planes");
    log_2 = Teuchos::rcp(new std::ofstream(s2.c_str(),std::ios::out));

    // Output of mean velocity and pressure of evaluation planes
    (*log_2) << "Output file of FDA blood flow benchmark – evaluation of evaluation planes     \n";
    (*log_2) << "------------------------------------------------------------------------------\n";
    (*log_2) << "                                                                              \n";
    (*log_2) << "Number_of_samples: " << numsamp_ << "                                         \n";
    (*log_2) << "\n\n";
    (*log_2) << "Radial_coordinate ";
    for(unsigned actEvalPlane=0; actEvalPlane<posEvaluation_.size();actEvalPlane++)
    (*log_2) << "Evaluation_plane z-coordinate_of_evaluation_plane Radial_coordinate  Mean_values_of_velocity_w_in_radial_direction Mean_values_of_pressure_p_in_radial_direction ";

    (*log_2) << "\n";

    for (int i=0; i<rcoordinates_.numRows(); ++i)
    {
        for (int j=0; j<rcoordinates_.numCols(); ++j)
        {
          (*log_2) <<  " "  << std::setw(11) << std::setprecision(4) << j+1;
          (*log_2) <<  " "  << std::setw(11) << std::setprecision(4) << actPosEvaluation_[j];
          (*log_2) <<  " "  << std::setw(11) << std::setprecision(4) << rcoordinates_[j][i];
          (*log_2) << "   " << std::setw(11) << std::setprecision(4) << (*rsumw_)(i,j);
          (*log_2) << "   " << std::setw(11) << std::setprecision(4) << (*rsump_)(i,j);
        }
        (*log_2) << "\n";
    }
    log->flush();
    log_2->flush();
  }
  return;
}// TurbulenceStatisticsBfda::DumpStatistics
