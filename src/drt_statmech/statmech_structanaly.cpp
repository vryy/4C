/*!----------------------------------------------------------------------
\file statmech_structanaly.cpp
\brief structural analysis and output for StatMech network structures

\maintainer Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276

*----------------------------------------------------------------------*/

#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_truss3/truss3.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"


//MEASURETIME activates measurement of computation time for certain parts of the code
//#define MEASURETIME


/*------------------------------------------------------------------------------*
 | output nodal displacements                                      mueller 5/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputNodalDisplacements(const Epetra_Vector&      disrow,
                                                         const std::ostringstream& filename)
{
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");

    std::stringstream dispnode;
    // retrieve translational node displacements
    for(int i=0; i<discret_->DofColMap()->NumMyElements(); i=i+6)
      dispnode << discol[i]<<" "<< discol[i+1]<<" "<< discol[i+2]<<std::endl;

    fputs(dispnode.str().c_str(), fp);
    fclose(fp);
  }
  return;
}

/*------------------------------------------------------------------------------*
 | output nodal positions                                          mueller 9/13 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputNodalPositions(const Epetra_Vector&      disrow,
                                                     const std::ostringstream& filename)
{
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");

    std::stringstream posnode;
    // retrieve translational node displacements
    for(int i=0; i<discret_->NodeColMap()->NumMyElements(); i++)
    {
      DRT::Node* colnode = discret_->lColNode(i);
      std::vector<int> dofnode = discret_->Dof(colnode);
      posnode <<std::scientific<<std::setprecision(8) << colnode->X()[0]+discol[dofnode[0]]<<"\t"
                                                      << colnode->X()[1]+discol[dofnode[1]]<<"\t"
                                                      << colnode->X()[2]+discol[dofnode[2]]<<std::endl;
    }
    fputs(posnode.str().c_str(), fp);
    fclose(fp);
  }
  return;
}

/*----------------------------------------------------------------------*
 | output for density-density-correlation-function(public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrOutput(const Epetra_Vector&      disrow,
                                             const std::ostringstream& filename,
                                             const int&                istep,
                                             const double&             dt)
{
  /*Output:
   * (1) structure number and characteristic length (radius, thickness)
   * (2) internal energy
   * (3) histogram of inter-crosslinker distance of crosslinker elements -> Density-Density-Correlation
   *     in three columns (length of one column: numbins) due to periodic continuation of the boundary volume
   * (4) histograms of spherical coordinates (azimuth angle phi, polar angle theta/ cos(theta)
   * (5) radial density distribution
   */
  if(!discret_->Comm().MyPID())
    std::cout<<"\n\n====================== Analysis of structural polymorphism ======================"<<std::endl;

  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2) || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
  // storage vector for shifted crosslinker LIDs(crosslinkermap)
  LINALG::Matrix<3,1> boxcenter;
  LINALG::Matrix<3,1> centershift;
  std::vector<int> crosslinkerentries;

#ifdef MEASURETIME
  double t0 = Teuchos::Time::wallTime();
#endif

  // determine the new center point of the periodic volume
  DDCorrShift(&boxcenter, &centershift, &crosslinkerentries);

#ifdef MEASURETIME
  double t1 = Teuchos::Time::wallTime();
#endif

  if(!discret_->Comm().MyPID())
    std::cout<<crosslinkerentries.size()<<" crosslinker elements found...\n"<<std::endl;

  // Determine current network structure
  LINALG::Matrix<3,1> cog;
  DDCorrCurrentStructure(disrow, &cog, &centershift, &crosslinkerentries, istep, filename);
  // store center of gravity for later use in GmshOutput()
  cog_ = cog;

#ifdef MEASURETIME
  double t2 = Teuchos::Time::wallTime();
#endif

  // Compute internal energy
  std::vector<double> internalenergy;
  internalenergy.clear();
  const Teuchos::RCP<Epetra_Vector> disp = Teuchos::rcp(new Epetra_Vector(disrow));
  ComputeInternalEnergy(disp, internalenergy,dt, filename);

#ifdef MEASURETIME
  double t3 = Teuchos::Time::wallTime();
#endif

  //calculcate the distance of crosslinker elements to all crosslinker elements (crosslinkermap)
  // testwise, base centershift upon calculated cog_ (rise in adequacy?)
  LINALG::Matrix<3,1> newcentershift;
  for(int i=0; i<(int)cog.M(); i++)
    newcentershift(i) = cog(i) - periodlength_->at(i)/2.0;

  // write the numbers of free, one-bonded, and two-bonded crosslink molecules
  CrosslinkCount(filename);

#ifdef MEASURETIME
  double t4 = Teuchos::Time::wallTime();
#endif

  // Orientation correlation function of filaments
  OrientationCorrelation(disrow, istep);

#ifdef MEASURETIME
  double t5 = Teuchos::Time::wallTime();
#endif

  // Compute the the network mesh size in dependency to the radial distance to a given COG
  ComputeLocalMeshSize(disrow, newcentershift, istep);

#ifdef MEASURETIME
  double t6 = Teuchos::Time::wallTime();
#endif

  // MultiVector because result vector will be of length 3*ddcorrrowmap_->MyLength()
  Epetra_MultiVector crosslinksperbinrow(*ddcorrrowmap_,9 , true);
  Epetra_MultiVector crosslinksperbinrotrow(*ddcorrrowmap_,3 , true);
  DDCorrFunction(crosslinksperbinrow, crosslinksperbinrotrow, &newcentershift);

#ifdef MEASURETIME
  double t7 = Teuchos::Time::wallTime();
#endif

  // calculation of filament element orientation in spherical coordinates, sorted into histogram
  Epetra_Vector phibinsrow(*ddcorrrowmap_, true);
  Epetra_Vector thetabinsrow(*ddcorrrowmap_, true);
  Epetra_Vector costhetabinsrow(*ddcorrrowmap_, true);
  SphericalCoordsDistribution(disrow, phibinsrow, thetabinsrow, costhetabinsrow);

#ifdef MEASURETIME
  double t8 = Teuchos::Time::wallTime();
#endif

  Epetra_Vector radialdistancesrow(*ddcorrrowmap_, true);
  RadialDensityDistribution(radialdistancesrow, centershift);

#ifdef MEASURETIME
  double t9 = Teuchos::Time::wallTime();
#endif

  // Import
  Epetra_MultiVector crosslinksperbincol(*ddcorrcolmap_,crosslinksperbinrow.NumVectors() , true);
  Epetra_MultiVector crosslinksperbinrotcol(*ddcorrcolmap_,crosslinksperbinrotrow.NumVectors() , true);
  Epetra_Vector phibinscol(*ddcorrcolmap_, true);
  Epetra_Vector thetabinscol(*ddcorrcolmap_, true);
  Epetra_Vector costhetabinscol(*ddcorrcolmap_, true);
  Epetra_Vector radialdistancescol(*ddcorrcolmap_, true);
  Epetra_Import importer(*ddcorrcolmap_, *ddcorrrowmap_);
  crosslinksperbincol.Import(crosslinksperbinrow, importer, Insert);
  crosslinksperbinrotcol.Import(crosslinksperbinrotrow, importer, Insert);
  phibinscol.Import(phibinsrow, importer, Insert);
  thetabinscol.Import(thetabinsrow, importer, Insert);
  costhetabinscol.Import(costhetabinsrow, importer, Insert);
  radialdistancescol.Import(radialdistancesrow, importer, Insert);

  // Add the processor-specific data up
  std::vector<std::vector<int> > crosslinksperbin(numbins, std::vector<int>(crosslinksperbincol.NumVectors(),0));
  std::vector<std::vector<int> > crosslinksperbinrot(numbins, std::vector<int>(crosslinksperbinrotcol.NumVectors(),0));
  std::vector<int> phibins(numbins, 0);
  std::vector<int> thetabins(numbins, 0);
  std::vector<int> costhetabins(numbins, 0);
  std::vector<int> radialdistbins(numbins, 0);
  //int total = 0;
  for(int i=0; i<numbins; i++)
    for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
    {
      for(int col=0; col<crosslinksperbincol.NumVectors(); col++)
        crosslinksperbin[i][col] += (int)crosslinksperbincol[col][pid*numbins+i];
      for(int col=0; col<crosslinksperbinrotcol.NumVectors(); col++)
        crosslinksperbinrot[i][col] += (int)crosslinksperbinrotcol[col][pid*numbins+i];
      phibins[i] += (int)phibinscol[pid*numbins+i];
      thetabins[i] += (int)thetabinscol[pid*numbins+i];
      costhetabins[i] += (int)costhetabinscol[pid*numbins+i];
      radialdistbins[i] += (int)radialdistancescol[pid*numbins+i];
    }

  // write data to file
  if(discret_->Comm().MyPID()==0)
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream histogram;
    // first part of output
    for(int i=0; i<numbins; i++)
    {
      histogram<<i+1<<"    ";
      for(int j=0; j<(int)crosslinksperbin[i].size(); j++)
        histogram<<crosslinksperbin[i][j]<<"    ";
      for(int j=0; j<(int)crosslinksperbinrot[i].size(); j++)
        histogram<<crosslinksperbinrot[i][j]<<"    ";
      histogram<<"  "<<phibins[i]<<"    "<<thetabins[i]<<"    "<<costhetabins[i]<<"    "<<radialdistbins[i]<<std::endl;
      //std::cout<<"  "<<phibins[i]<<"    "<<thetabins[i]<<"    "<<costhetabins[i]<<"    "<<radialdistbins[i]<<std::endl;
    }

    fputs(histogram.str().c_str(), fp);
    fclose(fp);
  }
  if(!discret_->Comm().MyPID())
  {
#ifdef MEASURETIME
    std::cout<<"\n=================Time  Measurement================"<<std::endl;
    std::cout<<"StatMechOutput::DDCorrOutput"<<std::endl;
    std::cout<<"DDCorrShift                 :\t"<<std::setprecision(4)<<t1-t0<<"\ts"<<std::endl;
    std::cout<<"DDCorrCurrentStructure      :\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<std::endl;
    std::cout<<"ComputeInternalEnergy       :\t"<<std::setprecision(4)<<t3-t2<<"\ts"<<std::endl;
    std::cout<<"CrosslinkCount              :\t"<<std::setprecision(4)<<t4-t3<<"\ts"<<std::endl;
    std::cout<<"OrientationCorrelation      :\t"<<std::setprecision(4)<<t5-t4<<"\ts"<<std::endl;
    std::cout<<"ComputeLocalMeshSize        :\t"<<std::setprecision(4)<<t6-t5<<"\ts"<<std::endl;
    std::cout<<"DDCorrFunction              :\t"<<std::setprecision(4)<<t7-t6<<"\ts"<<std::endl;
    std::cout<<"SphericalCoordsDistribution :\t"<<std::setprecision(4)<<t8-t7<<"\ts"<<std::endl;
    std::cout<<"RadialDensityDistribution   :\t"<<std::setprecision(4)<<t9-t8<<"\ts"<<std::endl;
    std::cout<<"Communication               :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t9<<"\ts"<<std::endl;
    std::cout<<"=================================================="<<std::endl;
    std::cout<<"total time                  :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t0<<"\ts"<<std::endl;
#endif
    std::cout<<"================================================================================="<<std::endl;
  }
}//STATMECH::StatMechManager::DDCorrOutput()

/*------------------------------------------------------------------------------*
 | Selects raster point with the smallest average distance to all crosslinker   |
 | elements, makes it the new center of the boundary box and shifts crosslinker |
 | positions.                                                                   |
 |                                                        (public) mueller 11/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrShift(LINALG::Matrix<3,1>* boxcenter, LINALG::Matrix<3,1>* centershift, std::vector<int>* crosslinkerentries)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2) || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  double periodlength = periodlength_->at(0);
  int numrasterpoints = statmechparams_.get<int>("NUMRASTERPOINTS", 3);
  // smallest average distance among all the average distances between the rasterpoints and all crosslinker elements
  // (init with 2*pl,so that it is definitely overwritten by the first "real" value)
  double  smallestdistance = 2.0*periodlength_->at(0);

  //store crosslinker element position within crosslinkerbond to crosslinkerentries
  for(int i=0; i<crosslinkerbond_->MyLength(); i++)
    if((*crosslinkerbond_)[0][i]>-0.9 && (*crosslinkerbond_)[1][i]>-0.9)
      crosslinkerentries->push_back(i);

  int numcrossele = (int)crosslinkerentries->size();

  // determine the new center of the box
  if(numcrossele>0)
  {
    for(int i=0; i<numrasterpoints; i++)
      for(int j=0; j<numrasterpoints; j++)
        for(int k=0; k<numrasterpoints; k++)
        {
          double averagedistance = 0.0;
          LINALG::Matrix<3,1> currentrasterpoint;
          LINALG::Matrix<3,1> currentcentershift;

          // calculate current raster point
          currentrasterpoint(0) = i*periodlength/(numrasterpoints-1);
          currentrasterpoint(1) = j*periodlength/(numrasterpoints-1);
          currentrasterpoint(2) = k*periodlength/(numrasterpoints-1);

          // calculate the center shift (difference vector between regular center and new center of the boundary box)
          for(int l=0; l<(int)currentrasterpoint.M(); l++)
            currentcentershift(l) = currentrasterpoint(l)-periodlength/2.0;

          // calculate average distance of crosslinker elements to raster point
          for(int l=0; l<numcrossele; l++)
          {
            // get the crosslinker position in question and shift it according to new boundary box center
            LINALG::Matrix<3,1> distance;
            for(int m=0; m<(int)distance.M(); m++)
            {
              distance(m) = (*crosslinkerpositions_)[m][(*crosslinkerentries)[l]];
              if (distance(m) > periodlength+currentcentershift(m))
                distance(m) -= periodlength;
              if (distance(m) < 0.0+currentcentershift(m))
                distance(m) += periodlength;
            }
            distance -= currentrasterpoint;
            averagedistance += distance.Norm2();
          }
          averagedistance /= (double)crosslinkerentries->size();

          if(averagedistance<smallestdistance)
          {
            smallestdistance = averagedistance;
            *boxcenter = currentrasterpoint;
            *centershift = currentcentershift;
          }
        }
  }
  else
  {
    boxcenter->PutScalar(periodlength/2.0);
    centershift->PutScalar(0.0);
  }
  //if(!discret_->Comm().MyPID())
    //std::cout<<"Box Center(2): "<<(*boxcenter)[0]<<", "<<(*boxcenter)[1]<<", "<<(*boxcenter)[2]<<std::endl;
  return;
}//STATMECH::StatMechManager::DDCorrShift()

/*------------------------------------------------------------------------------*
 | Determine current network structure and output network type as single        |
 | characteristic number. Also, output filament orientations.                   |
 |                                                        (public) mueller 11/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrCurrentStructure(const Epetra_Vector& disrow,
                                             LINALG::Matrix<3,1>* cog,
                                             LINALG::Matrix<3,1>* centershift,
                                             std::vector<int>* crosslinkerentries,
                                             const int& istep,
                                             const std::ostringstream& filename,
                                             bool filorientoutput)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  double periodlength = periodlength_->at(0);
  // number of crosslinker elements
  int numcrossele = (int)crosslinkerentries->size();
  std::vector<int> crosslinksinvolume(3,0);

  // get column map displacements
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  /// calculate center of gravity with respect to shiftedpositions for bound crosslinkers
  std::vector<LINALG::Matrix<3,1> > shiftedpositions;
  cog->Clear();
  // shift positions according to new center point (raster point)
  for(int i=0; i<numcrossele; i++)
  {
    LINALG::Matrix<3,1> currposition;
    for(int j=0; j<(int)currposition.M(); j++)
    {
      currposition(j) = (*crosslinkerpositions_)[j][(*crosslinkerentries)[i]];
      if (currposition(j) > periodlength+(*centershift)(j))
        currposition(j) -= periodlength;
      if (currposition(j) < 0.0+(*centershift)(j))
        currposition(j) += periodlength;
    }
    (*cog) += currposition;
    if(discret_->Comm().MyPID()==0)
      shiftedpositions.push_back(currposition);
  }
  if(numcrossele != 0)
    cog->Scale(1.0/(double)numcrossele);

  // zero out for new run
  trafo_->Zero();

  // calculations done by Proc 0 only
  if(discret_->Comm().MyPID()==0)
  {
    // number indicating structure type
    int structurenumber = 0;
    // indices for layer plane vectors
    int dir1 = -1, dir2 = -1;
    double rlink = statmechparams_.get<double>("R_LINK", 1.0);
    // number of nested intervals
    int maxexponent = (int)ceil(log(periodlength/rlink)/log(2.0))*2;
    // vector for test volumes (V[0]-sphere, V[1]-cylinder, V[2]-layer/homogenous network)
    std::vector<double> volumes(3,pow(10*periodlength, 3.0));
    // characteristic lengths of "volumes"
    std::vector<double> characlength(3,10*periodlength);
    // crosslinker fraction included in test volume
    std::vector<double> crossfraction(3,0.0);
    // number of iterations until crosslinker fraction lies within the given threshold fraction +/- tolerance
    std::vector<int> niter(3,0);

    // iterated vectors
    // bundle
    LINALG::Matrix<3,1> cylvec;
    // layer
    std::vector<LINALG::Matrix<3,1> > layervectors;
    // cluster-layer vectors (layer vectors determined when actually, a cluster phase is detected)
    std::vector<LINALG::Matrix<3,1> > clusterlayervecs;

    // get the intersections of the axis of a cylinder with the (two) cube faces
    std::vector<LINALG::Matrix<3,1> > intersections;
    // coordinates of intersection points of layer-type volume with the cube edges
    std::vector<LINALG::Matrix<3,1> > interseccoords;

    if(numcrossele>0)
    {
  /// calculate normed vectors and output filament element orientations
      // normed vectors for structural analysis (projections of base vectors e1, e2, e3 onto structure)
      std::vector<LINALG::Matrix<3,1> > normedvectors;
      for(int i=0; i<3; i++)
      {
        LINALG::Matrix<3,1> normedi;
        normedi.Clear();
        normedvectors.push_back(normedi);
      }

  /// determine normed vectors as well as output filament element vectors
      std::ostringstream orientfilename;
      orientfilename << "./FilamentOrientations_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
      FilamentOrientations(discol, &normedvectors, orientfilename, filorientoutput);

      // select the vector best fitting the axis of the bundle cylinder: vector with the greatest length -> smallest
      // enclosed angle with the cylinder axis
      int startiterindex = 0;
      double veclength = 0.0;
      std::cout<<"Normed vectors:"<<std::endl;
      for(int i=0; i<(int)normedvectors.size(); i++)
      {
        if(normedvectors[i].Norm2() > veclength)
        {
          veclength = normedvectors[i].Norm2();
          startiterindex = i;
        }
        // Scale normed vectors to unit length
        normedvectors[i].Scale(1.0/normedvectors[i].Norm2());
        std::cout<<i<<": "<<normedvectors[i](0)<<" "<<normedvectors[i](1)<<" "<<normedvectors[i](2)<<std::endl;
      }

  /// determine network structure
      // threshold fraction of crosslinkers
      double pthresh = 0.9;

      // calculate smallest possible test volumes that fulfill the requirement /numcrossele >= pthresh
      for(int i=0; i<(int)volumes.size(); i++)
      {
        switch(i)
        {
          // spherical volume
          case 0:
          {
            bool leaveloop = false;
            // tolerance
            double lowertol = 0.05;
            double uppertol = 0.01;
            // initial search radius
            double radius = periodlength/2.0;
            // fraction of crosslinks within test volume
            double pr = 0.0;
            int exponent = 1;

            // determine the two normed vectors with the largest cross product 2-Norm (closest to being perpendicular)
            // this is done in order to capture the smooth transformation from cluster to layer. More precisely, we
            // calculate the layer triad. The information we gain from this procedure is useful once we observe a
            // filament orientation distribution which is no longer homogeneous.
            double crossprodnorm = 0.0;
            for(int j=0; j<3; j++)
              for(int k=0; k<3; k++)
                if(k>j)
                {
                  // cross product
                  LINALG::Matrix<3,1> crossvec;
                  crossvec(0) = normedvectors[j](1)*normedvectors[k](2) - normedvectors[j](2)*normedvectors[k](1);
                  crossvec(1) = normedvectors[j](2)*normedvectors[k](0) - normedvectors[j](0)*normedvectors[k](2);
                  crossvec(2) = normedvectors[j](0)*normedvectors[k](1) - normedvectors[j](1)*normedvectors[k](0);
                  if(crossvec.Norm2()>crossprodnorm)
                  {
                    crossprodnorm = crossvec.Norm2();
                    dir1=j;
                    dir2=k;
                  }
                }

            // store initial vectors to be iterated
            clusterlayervecs.push_back(normedvectors[dir1]);
            clusterlayervecs.push_back(normedvectors[dir2]);

            // iterate both vectors in order to obtain projections into the layer plane
            const int maxiterations = 25;
            std::cout<<"\nVector iteration:"<<std::endl;
            std::cout<<"Cluster1: ";
            DDCorrIterateVector(discol, &clusterlayervecs[0], maxiterations);
            std::cout<<"Cluster2: ";
            DDCorrIterateVector(discol, &clusterlayervecs[1], maxiterations);

            // cross product n_1 x n_2, plane normal
            LINALG::Matrix<3,1> normal;
            normal(0) = clusterlayervecs[0](1)*clusterlayervecs[1](2) - clusterlayervecs[0](2)*clusterlayervecs[1](1);
            normal(1) = clusterlayervecs[0](2)*clusterlayervecs[1](0) - clusterlayervecs[0](0)*clusterlayervecs[1](2);
            normal(2) = clusterlayervecs[0](0)*clusterlayervecs[1](1) - clusterlayervecs[0](1)*clusterlayervecs[1](0);
            normal.Scale(1.0/normal.Norm2());
            clusterlayervecs.push_back(normal);

            // loop as long as pr has not yet reached pthresh
            // flag indicating convergence (set to false, if max no. of iterations is reached)
            bool converged = true;
            while(!leaveloop)
            {
              int rcount = 0;
              // loop over crosslinker elements
              for(int j=0; j<numcrossele; j++)
              {
                // get distance of crosslinker element to center of gravity
                LINALG::Matrix<3,1> dist = shiftedpositions[j];
                dist -= *cog;
                if(dist.Norm2()<=radius)
                  rcount++;
              }
              pr = double(rcount)/double(numcrossele);

              exponent++;
              niter[0]++;
              // new radius
              if(exponent<=maxexponent)
              {
                if((pr<pthresh-lowertol || pr>pthresh+uppertol))
                {
                  // determine "growth direction"
                  double sign;
                  if(pr<pthresh)
                    sign = 1.0;
                  else
                    sign = -1.0;
                  radius += sign*periodlength/pow(2.0,(double)exponent);
                }
                else
                {
                  crosslinksinvolume[i] = rcount;
                  crossfraction[i] = pr;
                  leaveloop = true;
                }
              }
              else
              {
                converged = false;
                crosslinksinvolume[i] = rcount;
                crossfraction[i] = pr;
                leaveloop = true;
              }
            }
            // store characteristic length and test sphere volume
            if(converged)
            {
              // store cluster diameter
              characlength[i] = 2.0*radius;
              volumes[i] = 4/3 * M_PI * pow(radius, 3.0);
            }
          }
          break;
          // cylindrical volume
          case 1:
          {
            bool leaveloop = false;
            // tolerance
            double lowertol = 0.02;
            double uppertol = 0.02;
            double radius = periodlength/2.0;
            double cyllength = 0.0;
            double pr = 0.0;
            int exponent = 1;

            // iterate to obtain fitting normed direction
            LINALG::Matrix<3,1> normj = normedvectors[startiterindex];
            const int maxiterations = 25;
            std::cout<<"Bundle:   ";
            DDCorrIterateVector(discol, &normj, maxiterations);
            // for output
            cylvec = normj;

            // cube face boundaries of jk-surface of cubical volume
            LINALG::Matrix<3,2> surfaceboundaries;
            for(int j=0; j<(int)surfaceboundaries.M(); j++)
            {
              surfaceboundaries(j,0) = (*centershift)(j);
              surfaceboundaries(j,1) = (*centershift)(j)+periodlength;
            }
            for(int j=0; j<(int)surfaceboundaries.N(); j++)
              for(int k=0; k<3; k++)
                for(int l=0; l<3; l++)
                  if(l>k)
                    for(int m=0; m<3; m++)
                      if(m!=k && m!=l)
                      {
                        LINALG::Matrix<3,1> currentintersection;
                        // known intersection component
                        currentintersection(m) = surfaceboundaries(m,j);
                        // get line parameter
                        double lambdaline = (currentintersection(m)-(*cog)(m))/normj(m);
                        currentintersection(k) = (*cog)(k)+lambdaline*normj(k);
                        currentintersection(l) = (*cog)(l)+lambdaline*normj(l);
                        // check if intersection lies on volume boundary
                        if(currentintersection(k)<=surfaceboundaries(k,1) && currentintersection(k)>=surfaceboundaries(k,0) &&
                           currentintersection(l)<=surfaceboundaries(l,1) && currentintersection(l)>=surfaceboundaries(l,0))
                          intersections.push_back(currentintersection);
                      }
            LINALG::Matrix<3,1> deltaisecs = intersections[1];
            deltaisecs -= intersections[0];
            cyllength = deltaisecs.Norm2();

            // calculate fraction of crosslinkers within cylinder
            // flag indicating convergence (set to false, if max no. of iterations is reached)
            bool converged = true;
            while(!leaveloop)
            {
              int rcount = 0;
              // loop over crosslinker elements
              for(int j=0; j<numcrossele; j++)
              {
                // get distance of crosslinker element to normed1 through center of gravity
                // intersection line-plane
                LINALG::Matrix<3,1> crosstocog = shiftedpositions[j];
                crosstocog -= *cog;

                double numerator = crosstocog.Dot(normj);
                double denominator = (normj).Dot(normj);
                double lambdaisec = numerator/denominator;
                // intersection and distance of crosslinker to intersection
                LINALG::Matrix<3,1> isecpt = normj;
                isecpt.Scale(lambdaisec);
                isecpt += (*cog);
                LINALG::Matrix<3,1> deltaiseccog = shiftedpositions[j];
                deltaiseccog -= isecpt;
                double distance = deltaiseccog.Norm2();

                if(distance<=radius)
                  rcount++;
              }
              pr = double(rcount)/double(numcrossele);

              exponent++;
              niter[1]++;
              // new radius
              if(exponent<=maxexponent)
              {
                if((pr<pthresh-lowertol || pr>pthresh+uppertol))
                {
                  // determine "growth direction"
                  double sign;
                  if(pr<pthresh)
                    sign = 1.0;
                  else
                    sign = -1.0;
                  radius += sign*periodlength/pow(2.0,(double)exponent);
                }
                else
                {
                  crosslinksinvolume[i] = rcount;
                  crossfraction[i] = pr;
                  leaveloop = true;
                }
              }
              else
              {
                converged = false;
                crosslinksinvolume[i] = rcount;
                crossfraction[i] = pr;
                leaveloop = true;
              }
            }
            if(converged)
            {
              // store bundle diameter
              characlength[i] = 2.0*radius;
              volumes[i] = M_PI*radius*radius*cyllength;
            }
          }
          break;
          // cuboid layer volume
          case 2:
          {
            bool leaveloop = false;
            // tolerance
            double lowertol = 0.02;
            double uppertol = 0.02;
            double thickness = periodlength/2.0;
            double pr = 0.0;
            int exponent = 1;

            //determine the two normed vectors with the largest cross product 2-Norm (closest to being perpendicular)
            double crossprodnorm = 0.0;
            for(int j=0; j<3; j++)
              for(int k=0; k<3; k++)
                if(k>j)
                {
                  // cross product
                  LINALG::Matrix<3,1> crossvec;
                  crossvec(0) = normedvectors[j](1)*normedvectors[k](2) - normedvectors[j](2)*normedvectors[k](1);
                  crossvec(1) = normedvectors[j](2)*normedvectors[k](0) - normedvectors[j](0)*normedvectors[k](2);
                  crossvec(2) = normedvectors[j](0)*normedvectors[k](1) - normedvectors[j](1)*normedvectors[k](0);
                  if(crossvec.Norm2()>crossprodnorm)
                  {
                    crossprodnorm = crossvec.Norm2();
                    dir1=j;
                    dir2=k;
                  }
                }

            // store initial vectors to be iterated
            layervectors.push_back(normedvectors[dir1]);
            layervectors.push_back(normedvectors[dir2]);

            // iterate both vectors in order to obtain projections into the layer plane
            const int maxiterations = 25;
            std::cout<<"Layer1:  ";
            DDCorrIterateVector(discol, &layervectors[0], maxiterations);
            std::cout<<"Layer2:  ";
            DDCorrIterateVector(discol, &layervectors[1], maxiterations);

            // cross product n_1 x n_2, plane normal
            LINALG::Matrix<3,1> normal;
            normal(0) = layervectors[0](1)*layervectors[1](2) - layervectors[0](2)*layervectors[1](1);
            normal(1) = layervectors[0](2)*layervectors[1](0) - layervectors[0](0)*layervectors[1](2);
            normal(2) = layervectors[0](0)*layervectors[1](1) - layervectors[0](1)*layervectors[1](0);
            normal.Scale(1.0/normal.Norm2());
            layervectors.push_back(normal);

            // crosslinkerpositions which are considered to be within the cuboid
            std::vector<LINALG::Matrix<3,1> > crosslinkswithinvolume;

            // flag indicating convergence (set to false, if max no. of iterations is reached)
            bool converged = true;
            while(!leaveloop)
            {
              for(int j=0; j<numcrossele; j++)
              {
                // given, that cog E plane with normal vector "normal"
                // constant in Hessian normal form
                double d = normal.Dot((*cog));
                // distance of crosslinker element to plane
                double pn = normal.Dot(shiftedpositions[j]);
                double disttoplane = fabs(pn-d);

                if(disttoplane <= thickness)
                  crosslinkswithinvolume.push_back(shiftedpositions[j]);
              }
              pr = double(crosslinkswithinvolume.size())/double(numcrossele);

              exponent++;
              niter[2]++;

              if(exponent<=maxexponent)
              {
                if((pr<pthresh-lowertol || pr>pthresh+uppertol))
                {
                  double sign;
                  if(pr<pthresh)
                    sign = 1.0;
                  else
                    sign = -1.0;
                  thickness += sign*periodlength/pow(2.0,(double)exponent);
                  crosslinkswithinvolume.clear();
                }
                else
                {
                  crosslinksinvolume[i] = (int)crosslinkswithinvolume.size();
                  crossfraction[i] = pr;
                  leaveloop = true;
                }
              }
              else
              {
                converged = false;
                crosslinksinvolume[i] = (int)crosslinkswithinvolume.size();
                crossfraction[i] = pr;
                leaveloop = true;
              }
            }
            // calculation of the volume
            // according to number of intersection points (only if crosslinker fraction calculation converged)
            if(converged)
            {
              // cube face boundaries of jk-surface of cubical volume
              LINALG::Matrix<3,2> surfaceboundaries;
              for(int j=0; j<(int)surfaceboundaries.M(); j++)
              {
                surfaceboundaries(j,0) = (*centershift)(j);
                surfaceboundaries(j,1) = (*centershift)(j)+periodlength;
              }

              // first step: find the intersection points layer x volume egdes.
              // At first, we assume a layer delimited by the volume boundaries)
              for(int surf=0; surf<(int)surfaceboundaries.N(); surf++) // (two) planes perpendicular to l-direction
              {
                for(int j=0; j<3; j++) // spatial component j
                  for(int k=0; k<3; k++) // spatial component k
                    if(k>j)// above diagonal
                      for(int l=0; l<3; l++) // spatial component l
                        if(l!=j && l!=k)
                        {
                          LINALG::Matrix<3,1> coords;
                          coords(l) = surfaceboundaries(l,surf);
                          for(int edge=0; edge<(int)surfaceboundaries.N(); edge++)
                            for(int m=0; m<(int)surfaceboundaries.M(); m++)
                              if(l!=m)
                              {
                                coords(m) = surfaceboundaries(m,edge);
                                for(int n=0; n<(int)surfaceboundaries.M(); n++)
                                  if(n!=m && n!=l)
                                  {
                                    coords(n) = (((*cog)(l)-coords(l))*normal(l) - (coords(m)-(*cog)(m))*normal(m))/normal(n) + (*cog)(n);
                                    double lowerbound = surfaceboundaries(n,0);
                                    double upperbound = surfaceboundaries(n,1);
                                    if((coords(n)>= lowerbound || fabs(coords(n)-lowerbound)<1e-6) && (coords(n)<= upperbound || fabs(coords(n)-upperbound)<1e-6))
                                    {
                                      bool redundant = false;
                                      // check for redundant entries
                                      if((int)interseccoords.size()>0)
                                      {
                                        LINALG::Matrix<3,1> check = coords;
                                        for(int p=0; p<(int)interseccoords.size(); p++)
                                        {
                                          check -= interseccoords[p];
                                          if(check.Norm2()<1e-4)
                                            redundant = true;
                                          check = coords;
                                        }
                                      }
                                      if(!redundant)
                                        interseccoords.push_back(coords);
                                    }
                                  }
                              }
                        }
              }

              //second step: Approximate the base area of the volume/prism
              // calculate and store unit direction vectors corners-COG
              std::vector<LINALG::Matrix<3,1> > itercoords = interseccoords;
              std::vector<LINALG::Matrix<3,1> > dirvecs;
              std::vector<double> initdistances;
              //std::cout<<"initdist = ";
              for(int j=0; j<(int)interseccoords.size(); j++)
              {
                // directional vectors pointing from center of gravity towards the corners
                LINALG::Matrix<3,1> dirvec = interseccoords[j];
                dirvec -= *cog;
                initdistances.push_back(dirvec.Norm2());
                //std::cout<<dirvec.Norm2()<<" ";
                dirvec.Scale(1.0/dirvec.Norm2());
                dirvecs.push_back(dirvec);
              }
              //std::cout<<std::endl;

              // calculate crosslinker projection positions in layer plane
              std::vector<LINALG::Matrix<3,1> > crossprojections;
              for(int j=0; j<(int)crosslinkswithinvolume.size(); j++)
              {
                // calculate projection of crosslinker onto layer plane
                LINALG::Matrix<3,1> projection = crosslinkswithinvolume[j];
                LINALG::Matrix<3,1> diffvec = normal;
                // line paramater
                double mu = (interseccoords[0].Dot(normal) - crosslinkswithinvolume[j].Dot(normal)) / (normal.Dot(normal));
                diffvec.Scale(mu);
                projection += diffvec;
                crossprojections.push_back(projection);
              }

              // calculate edge directions of the layer poygon
              std::vector<LINALG::Matrix<3,1> > edgedirs;
              std::vector<int> neworder;
              for(int j=0; j<(int)interseccoords.size(); j++)
                for(int k=0; k<(int)interseccoords.size(); k++)
                  if(k>j)
                    for(int l=0; l<(int)interseccoords[j].M(); l++)
                      if(fabs((interseccoords[j])(l)-(interseccoords[k])(l))<1e-7)
                      {
                        LINALG::Matrix<3,1> jdir = interseccoords[k];
                        jdir -= interseccoords[j];
                        jdir.Scale(1.0/jdir.Norm2());
                        edgedirs.push_back(jdir);
                        neworder.push_back(j);
                      }

              // iterate as long as layer volume contains 90-95% of the crosslinks
              leaveloop = false;
              int numiter = 30;
              int iter = 0;
              int rcount = 0;
              double uplim = 0.01;
              double lowlim = 0.03;
              double threshold = 0.95;
              // start with 1.0 (100% of crosslinkers within layer volume)
              double crossfrac = 1.0;
              while(!leaveloop)
              {
                if(iter<numiter)
                {
                  if(crossfrac<threshold-lowlim || crossfrac>threshold+uplim)
                  {
                    //calculate new set of layer corner coordinates
                    for(int j=0; j<(int)itercoords.size(); j++)
                    {
                      LINALG::Matrix<3,1> deltaj = dirvecs[j];
                      deltaj.Scale(initdistances[j]/pow(1.4,(double)(iter+1)));
                      // depending on crossfrac, choose in which direction the triangle grows:
                      if(crossfrac>threshold)
                        itercoords[j] -= deltaj;
                      if(crossfrac<threshold)
                        itercoords[j] += deltaj;
                    }

                    // check if crosslinkers lie within the volume with the new base area
                    // number of crosslinkers
                    rcount = 0;
                    for(int j=0; j<(int)crosslinkswithinvolume.size(); j++)
                    {
                      // distance of crosslinker projection to center of gravity
                      LINALG::Matrix<3,1> crosstocog = crossprojections[j];
                      crosstocog -= *cog;
                      double dcrosstocog = crosstocog.Norm2();
                      crosstocog.Scale(1.0/dcrosstocog);

                      // calculate shortest positive distance for intersections of the line (cog->crosslinker) with the layer polygon from the center of gravity
                      double dclosestisec = 10*periodlength;
                      for(int k=0; k<(int)edgedirs.size(); k++)
                      {
                        // line paramete (=length)
                        double numerator = edgedirs[k](1)*(itercoords[neworder[k]](0)-(*cog)(0)) - edgedirs[k](0)*(itercoords[neworder[k]](1)-(*cog)(1));
                        double denominator = crosstocog(0)*edgedirs[k](1) - crosstocog(1)*edgedirs[k](0);
                        // lambda either positive (i.e. intersection lies in the direction of the crosslinker) or negative
                        double lambda = numerator/denominator;
                        // in case of parallel directional vectors, lambda takes big value (does not matter anyway as this case is not important for the coming calculations)
                        if(fabs(denominator) < 1e-8)
                          lambda = 1e8;
                        // save smallest positive lambda. It marks the closest intersection of the (cog->crosslinker)-line with an edge
                        if(lambda>0.0 && lambda<dclosestisec)
                          dclosestisec = lambda;
                      }
                      // if the distance to the closest intersection is bigger than the distance between crosslinker and center of gravity,
                      // the linker in question lies within the iterated volume
                      if(dclosestisec>=dcrosstocog)
                        rcount++;
                    }
                    //new crosslinker fraction in volume
                    crossfrac = double(rcount)/double(crosslinkswithinvolume.size());
                    //std::cout<<"i="<<iter<<", crossfrac = "<<crossfrac<<", rcount = "<<rcount<<"/"<<crosslinkswithinvolume.size()<<std::endl;
                    iter++;
                  }
                  else
                  {
                    // set itercoords as new interseccoords
                    interseccoords = itercoords;
                    crosslinksinvolume[i] = rcount;
                    crossfraction[i] *= pr;
                    std::cout<<"-> adjusted volume after "<<iter<<" iterations with "<<rcount<<"/"<<crosslinkswithinvolume.size()<<"( p="<<crossfrac<<" )"<<std::endl;
                    leaveloop = true;
                  }
                }
                else
                {
                  // set itercoords as new interseccoords
                  interseccoords = itercoords;
                  crosslinksinvolume[i] = rcount;
                  crossfraction[i] *= pr;
                  std::cout<<"-> adjusted volume after maxiter = "<<iter<<" iterations with "<<rcount<<"/"<<crosslinkswithinvolume.size()<<"( p="<<crossfrac<<" )"<<std::endl;
                  leaveloop = true;
                }
              }

              // calculation of layer volume
              switch((int)interseccoords.size())
              {
                // triangle
                case 3:
                {
                  // area of the triangle via cosine rule
                  LINALG::Matrix<3,1> c = interseccoords[0];
                  c -= interseccoords[1];
                  LINALG::Matrix<3,1> a = interseccoords[1];
                  a -= interseccoords[2];
                  LINALG::Matrix<3,1> b = interseccoords[0];
                  b -= interseccoords[2];
                  double cl = c.Norm2();
                  double al = a.Norm2();
                  double bl = b.Norm2();
                  double alpha = acos((cl*cl+bl*bl-al*al)/(2.0*bl*cl));
                  double h = bl * sin(alpha);

                  // volume of layer (factor 0.5 missing, since "real" thickness is thickness*2.0)
                  volumes[i] = cl*h*thickness;
                }
                break;
                // square/rectangle/trapezoid
                case 4:
                {
                  // edges
                  LINALG::Matrix<3,1> a = interseccoords[1];
                  a -= interseccoords[0];
                  LINALG::Matrix<3,1> c = interseccoords[3];
                  c -= interseccoords[2];
                  LINALG::Matrix<3,1> d = interseccoords[2];
                  d -= interseccoords[0];
                  // diagonal
                  LINALG::Matrix<3,1> f = interseccoords[2];
                  f -= interseccoords[1];
                  double al = a.Norm2();
                  double cl = c.Norm2();
                  double dl = d.Norm2();
                  double fl = f.Norm2();
                  double alpha = acos((al*al+dl*dl-fl*fl)/(2.0*al*dl));
                  double h = dl * sin(alpha);

                  volumes[i] = (al+cl) * h * thickness;
                  //std::cout<<"layer volume_rec = "<<volumes[i]<<std::endl;
                }
                break;
                // hexahedron
                case 6:
                {
                  double hexvolume = 0.0;
                  // indices mapping correct order of hexagon edges
                  for(int j=0; j<(int)interseccoords.size(); j++)
                    for(int k=0; k<(int)interseccoords.size(); k++)
                      if(k<j)
                        for(int l=0; l<3; l++)
                          if(fabs(interseccoords[i](l)-interseccoords[j](l))<1e-7)// components identical
                          {
                            // get edge of j-th triangle within hexagon
                            LINALG::Matrix<3,1> a = interseccoords[k];
                            a -= *cog;
                            LINALG::Matrix<3,1> b = interseccoords[j];
                            b -= *cog;
                            LINALG::Matrix<3,1> c = interseccoords[k];
                            c -= interseccoords[j];

                            double al = a.Norm2();
                            double bl = b.Norm2();
                            double cl = c.Norm2();
                            double alpha = acos((cl*cl+bl*bl-al*al)/(2.0*bl*cl));
                            double h = bl * sin(alpha);

                            hexvolume += cl*h*thickness;
                          }
                  volumes[i] = hexvolume;
                }
                break;
              }
              characlength[i] = 2.0*thickness;
            }
          }
          break;
        }
      }
    }
    // smallest volume
    double minimalvol = 9e99;
    int minimum = 0;
    for(int i=0; i<(int)volumes.size(); i++)
      if(volumes[i]<minimalvol)
      {
        minimalvol = volumes[i];
        structurenumber = i;
        minimum = i;
      }
    if(structurenumber==0 && characlength[0]>=periodlength)
      structurenumber = 3;

    std::cout<<"\nVolumes: "<<std::endl;
    for(int i=0; i<(int)volumes.size(); i++)
      std::cout<<std::fixed<<std::setprecision(6)<<"V("<<i<<"): "<<volumes[i]<<"  l_c: "<<characlength[i]<<"  p_cross: "<<crossfraction[i]<<" crosslinks: "<<crosslinksinvolume[i]<<"/"<<numcrossele<<"  niter: "<<niter[i]<<std::endl;

  /// std::cout and return network structure
    // write to output files
    // append structure number to DDCorr output
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");
    std::stringstream structuretype;
    structuretype<<structurenumber<<"    "<<characlength[minimum];
    for(int j=0; j<15; j++)
      structuretype<<"    "<<0.0;
    structuretype<<std::endl;
    fputs(structuretype.str().c_str(), fp);
    fclose(fp);

    // clear the vector first
    testvolumepos_.clear();
    if(numcrossele>0)
    {
      switch(structurenumber)
      {
        // cluster
        case 0:
        {
          std::cout<<"\nNetwork structure: Cluster"<<std::endl;
          characlength_ = characlength[structurenumber];
          structuretype_ = structurenumber;

          // calculate the trafo_ matrix (as if we had a layer)

          // adjust second plane direction so that we get an orthonormal basis
          clusterlayervecs[1](0) = clusterlayervecs[0](1)*clusterlayervecs[2](2) - clusterlayervecs[0](2)*clusterlayervecs[2](1);
          clusterlayervecs[1](1) = clusterlayervecs[0](2)*clusterlayervecs[2](0) - clusterlayervecs[0](0)*clusterlayervecs[2](2);
          clusterlayervecs[1](2) = clusterlayervecs[0](0)*clusterlayervecs[2](1) - clusterlayervecs[0](1)*clusterlayervecs[2](0);
          clusterlayervecs[1].Scale(1/clusterlayervecs[1].Norm2()); // hm, not necessary

          // build the base
          for(int i=0; i<trafo_->M(); i++)
            for(int j=0; j<trafo_->N(); j++)
              (*trafo_)(i,j) = clusterlayervecs[i](j);
        }
        break;
        // bundle
        case 1:
        {
          std::cout<<"\nNetwork structure: Bundle"<<std::endl;
          std::cout<<"axis vector: "<<cylvec(0)<<" "<<cylvec(1)<<" "<<cylvec(2)<<std::endl;
          structuretype_ = structurenumber;
          characlength_ = characlength[structurenumber];
          for(int i=0; i<(int)intersections.size(); i++)
            testvolumepos_.push_back(intersections[i]);

          // save trafo matrix for later use in DDCorrFunction()
          // second direction
          LINALG::Matrix<3,1> raddir1;
          raddir1(0) = 0.0;
          raddir1(1) = 1.0;
          raddir1(2) = cylvec(1)/cylvec(2);
          raddir1.Scale(1.0/raddir1.Norm2());
          // third direction
          LINALG::Matrix<3,1> raddir2;
          raddir2(0) = raddir1(1)*cylvec(2) - raddir1(2)*cylvec(1);
          raddir2(1) = raddir1(2)*cylvec(0) - raddir1(0)*cylvec(2);
          raddir2(2) = raddir1(0)*cylvec(1) - raddir1(1)*cylvec(0);
          raddir2.Scale(1.0/raddir2.Norm2());

          // build the base
          for(int j=0; j<trafo_->N(); j++)
          {
            (*trafo_)(0,j) = cylvec(j);
            (*trafo_)(1,j) = raddir1(j);
            (*trafo_)(2,j)= raddir2(j);
          }
        }
        break;
        // layer
        case 2:
        {
          std::cout<<"\nNetwork structure: Layer ( ";
          switch((int)(interseccoords.size()))
          {
            case 3: std::cout<<"triangular shape )"; break;
            case 4: std::cout<<"rectangular shape )"; break;
            case 6: std::cout<<"haxagonal shape )"; break;
          }
          std::cout<<std::endl;
          for(int i=0; i<(int)layervectors.size(); i++)
            std::cout<<"layer vector "<<i+1<<": "<<layervectors[i](0)<<" "<<layervectors[i](1)<<" "<<layervectors[i](2)<<std::endl;
          structuretype_ = structurenumber;
          characlength_ = characlength[structurenumber];
          for(int i=0; i<(int)interseccoords.size(); i++)
            testvolumepos_.push_back(interseccoords[i]);

          // calculate second plane vector (which is now exactly orthogonal to vec1)
          layervectors[1](0) = layervectors[0](1)*layervectors[2](2) - layervectors[0](2)*layervectors[2](1);
          layervectors[1](1) = layervectors[0](2)*layervectors[2](0) - layervectors[0](0)*layervectors[2](2);
          layervectors[1](2) = layervectors[0](0)*layervectors[2](1) - layervectors[0](1)*layervectors[2](0);
          layervectors[1].Scale(1/layervectors[1].Norm2()); // hm, not necessary

          // build the base
          for(int i=0; i<trafo_->M(); i++)
            for(int j=0; j<trafo_->N(); j++)
              (*trafo_)(i,j) = layervectors[i](j);
        }
        break;
        // homogeneous
        case 3:
        {
          std::cout<<"\nNetwork structure: Homogeneous network"<<std::endl;
          structuretype_ = structurenumber;
          characlength_ = characlength[0];

          // save trafo matrix for later use in DDCorrFunction()
          for(int i=0; i<trafo_->M(); i++)
            (*trafo_)(i,i) = 1.0;
        }
        break;
      }
    }
    else
    {
      std::cout<<"\nNetwork structure: Homogeneous network"<<std::endl;
      structuretype_ = structurenumber;
      characlength_ = characlength[0];

      // save trafo matrix for later use in DDCorrFunction()
      for(int i=0; i<trafo_->M(); i++)
        (*trafo_)(i,i) = 1.0;
    }
  }

  // Communicate trafo_ to other procs
  std::vector<double> localtrafo(9,0.0);
  std::vector<double> globaltrafo(9,0.0);
  for(int i=0; i<trafo_->M(); i++)
    for(int j=0; j<trafo_->N(); j++)
    {
      localtrafo.at(3*i+j) = (*trafo_)(i,j);
      discret_->Comm().SumAll(&localtrafo[3*i+j], &globaltrafo[3*i+j], 1);
      (*trafo_)(i,j) = globaltrafo.at(3*i+j);
    }
  //std::cout<<*trafo_<<std::endl;
}//DDCorrCurrentStructure()

/*------------------------------------------------------------------------------*                                                 |
 | density-density correlation function                   (private) mueller 01/11|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrIterateVector(const Epetra_Vector& discol, LINALG::Matrix<3,1>* vectorj, const int& maxiterations)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  // get filament number conditions
  std::vector<DRT::Condition*> filaments(0);
  discret_->GetCondition("FilamentNumber", filaments);
  bool vectorconverged = false;
  int iteration = 0;
  double periodlength = periodlength_->at(0);
  double tolangle = M_PI/180.0; // 1°

  while(!vectorconverged)
  {
    if(iteration<maxiterations)
    {
      LINALG::Matrix<3,1> vectorjp;
      vectorjp.Clear();
      for(int fil=0; fil<(int)filaments.size(); fil++)
      {
        // get next filament
        DRT::Condition* currfilament = filaments[fil];
        for(int node=1; node<(int)currfilament->Nodes()->size(); node++)
        {
          // obtain column map LIDs
          int gid0 = currfilament->Nodes()->at(node-1);
          int gid1 = currfilament->Nodes()->at(node);
          int nodelid0 = discret_->NodeColMap()->LID(gid0);
          int nodelid1 = discret_->NodeColMap()->LID(gid1);
          DRT::Node* node0 = discret_->lColNode(nodelid0);
          DRT::Node* node1 = discret_->lColNode(nodelid1);

          // calculate directional vector between nodes
          LINALG::Matrix<3, 1> dirvec;
          for(int dof=0; dof<3; dof++)
          {
            int dofgid0 = discret_->Dof(node0)[dof];
            int dofgid1 = discret_->Dof(node1)[dof];
            double poscomponent0 = node0->X()[dof]+discol[discret_->DofColMap()->LID(dofgid0)];
            double poscomponent1 = node1->X()[dof]+discol[discret_->DofColMap()->LID(dofgid1)];
            // check for periodic boundary shift and correct accordingly
            if (fabs(poscomponent1 - periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
              poscomponent1 -= periodlength;
            else if (fabs(poscomponent1 + periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
              poscomponent1 += periodlength;

            dirvec(dof) = poscomponent1-poscomponent0;
          }
          // normed vector
          dirvec.Scale(1.0/dirvec.Norm2());

          if(acos(dirvec.Dot((*vectorj)))>(M_PI/2.0))
            dirvec.Scale(-1.0);
          vectorjp += dirvec;
        }
      }
      vectorjp.Scale(1.0/vectorjp.Norm2());
      // check angle between old n_j and n_(j+1)
      double vecvecangle = acos(vectorjp.Dot((*vectorj)));
      if(vecvecangle < tolangle)
      {
        std::cout<<" vector converged after "<<iteration+1<<" iteration(s) with angle "<<vecvecangle/M_PI*180.0<<" deg"<<std::endl;
        vectorconverged = true;
      }
      *vectorj = vectorjp;
    }
    else
    {
      std::cout<<"...Vector did not converge after "<<maxiterations<<" iterations. Continuing..."<<std::endl;
      vectorconverged = true;
    }
    iteration++;
  }
}//DDCorrIterateVector()
/*------------------------------------------------------------------------------*                                                 |
 | density-density correlation function                  (private) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrFunction(Epetra_MultiVector& crosslinksperbinrow, Epetra_MultiVector& crosslinksperbinrotrow, LINALG::Matrix<3,1>* centershift)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
  double periodlength = periodlength_->at(0);


/// preliminary operations to set workframe
  // exporter and importer
  Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
  Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);
  // tranfer map format of crosslinkerbond_ and crosslinkerpositions_
  Epetra_MultiVector crosslinkerbondtrans(*transfermap_,2,true);
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_,3,true);
  if(discret_->Comm().MyPID()!=0)
  {
    crosslinkerbond_->PutScalar(0.0);
    crosslinkerpositions_->PutScalar(0.0);
  }
  // distribute information to processors
  crosslinkerbondtrans.Export(*crosslinkerbond_, crosslinkexporter, Add);
  crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);

  // reproduce redundancy prior to Export (not needed in this method, but in the course of DDCorroutput())
  crosslinkerbond_->Import(crosslinkerbondtrans, crosslinkimporter, Insert);
  crosslinkerpositions_->Import(crosslinkerpositionstrans, crosslinkimporter, Insert);

/// preparations for calculation of inter-crosslink distances between central box and surrounding box crosslinker elements
  /* rotate crosslinker position into new (material) coordinate system, then construct
   * new surrounding box crosslinker position with periodic continuations according to the new
   * coordinate directions. This will in some cases lead to a cutoff of crosslinkers which come to lie
   * outside of the rotated volume. These crosslinkers are neglected.
   * The reasons why we are allowed do this, lie with the structure of the network as well as the object of observation.
   * When considering the worst case, we lose 19.6% of the volume (rotation around two axis, 45° each). Still, this does
   * not matter in our case. Since we want to study stationary phases like the cluster or the layer phase,
   * we only consider time intervals where the phase in question does not change its shape anymore. In addition,
   * we try to center the structure as well as possible within the boundary volume by shifting the center of the
   * periodic box. Clusters are not affected by the rotation and cutoff since they constitute a structural singularity,
   * so to speak, and therefore have no contace to the volume boundaries. They are not affected by periodic BCs and are
   * strongly localized density peaks. Layers tend to span over at least one direction of the periodic boundary volume.
   * If this stage is reached, the cutoff will not affect the layer, i.e. its crosslinkers at all, since there is no material
   * in the corners of the boundary volume anymore.
   * Ergo, we wait until this state is reached and then start to evaluate (done via Matlab, not directly implemented).
   * For homogeneous networks, of course, things are different. However, we do not need to worry about it since we have
   * equally distributed filament directions and hence no need to rotate neither the coordinate system nor the entire volume.
   */
/// calculate  shifted as well as shifted and rotated crosslinker positions of the center box
  Epetra_MultiVector centerboxpostrans(*transfermap_,3);
  Epetra_MultiVector centerboxrotpostrans(*transfermap_,3);
  for(int i=0; i<transfermap_->NumMyElements(); i++)
  {
    LINALG::SerialDenseMatrix crosspos(3,1);
    for(int j=0; j<crosspos.M(); j++)
    {
      // 1. original position
      crosspos(j,0) = crosslinkerpositionstrans[j][i];
      // 2. shift, so that crosspos comes to lie within not yet rotated new center box
      if (crosspos(j,0) > periodlength+(*centershift)(j))
        crosspos(j,0) -= periodlength;
      if (crosspos(j,0) < 0.0+(*centershift)(j))
        crosspos(j,0) += periodlength;
      // 2.5 store standard shifted center box position
      centerboxpostrans[j][i] = crosspos(j,0);
      // 3. translate origin
      crosspos(j,0) -= cog_(j);
    }
    // 4. transform the crosslink molecule position into material coordinates
    LINALG::SerialDenseMatrix crossposrot(3,1);
    trafo_->Multiply(false, crosspos, crossposrot);
    for(int j=0; j<centerboxrotpostrans.NumVectors(); j++)
      centerboxrotpostrans[j][i] = crossposrot(j,0);
  }

/// calculate shifted and rotated crosslinker positions for the surrounding boxes (including the center box)
  Epetra_MultiVector positionstrans(*transfermap_, 3*27, true);
  Epetra_MultiVector rotatedpositionstrans(*transfermap_, 3*27, true);
  std::vector<std::vector<int> > boxindices;
  int boxnumber=0;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        // box indices
        std::vector<int> ijk;
        ijk.push_back(i);
        ijk.push_back(j);
        ijk.push_back(k);
        boxindices.push_back(ijk);
        // obtain standard crosslinker position and shift it according to current box with respect to box (i=1;j=1;k=1)
        for(int l=0; l<transfermap_->NumMyElements(); l++)
          for(int m=0; m<centerboxpostrans.NumVectors(); m++)
          {
            // shift standard crosslinker positions to box with indices (i,j,k)
            positionstrans[boxnumber*3+m][l] = centerboxpostrans[m][l] + (ijk[m]-1)*periodlength;
            // shift the rotated crosslinker positions to box with indices (i,j,k)
            rotatedpositionstrans[boxnumber*3+m][l] = centerboxrotpostrans[m][l] + (ijk[m]-1)*periodlength;
          }
        boxnumber++;
      }
  // communicate positions for redundancy (crosslinkermap format)
  Epetra_MultiVector positions(*crosslinkermap_, 3*27, true);
  Epetra_MultiVector rotatedpositions(*crosslinkermap_, 3*27, true);
  positions.Import(positionstrans, crosslinkimporter, Insert);
  rotatedpositions.Import(rotatedpositionstrans, crosslinkimporter, Insert);

  // retrieve center of gravity (global coordinates)
  Epetra_SerialDenseMatrix cog(3,1);
  for(int i=0; i<cog.M(); i++)
    cog(i,0) = cog_(i);

/// sort inter-crosslink distances into respective bins
  //parallel brute force from here on
  /* 1. loop over boxes
   * 2. loop over the crosslinkermap format crosslinkers of the boxes
   * 3. loop over the transfermap format crosslinkers of the actual (center box)
   * 4. if(): only calculate deltaxij if the following if both indices i and j stand for crosslinker elements or periodic
   *          images of the those elements
   */
  for(int boxnum=0; boxnum<boxnumber; boxnum++)
    for(int j=0; j<crosslinkermap_->NumMyElements(); j++)
      for(int i=0; i<transfermap_->NumMyElements(); i++)
        if(crosslinkerbondtrans[0][i]>-0.9 && crosslinkerbondtrans[1][i]>-0.9 && (*crosslinkerbond_)[0][j]>-0.9 && (*crosslinkerbond_)[1][j]>-0.9)
        {
          // sort standard (GLOBAL coordinates) inter-crosslink distances into bins
          for(int m=0; m<centerboxpostrans.NumVectors(); m++)
          {
            // centerbox crosslinker component m
            double cboxposm = centerboxpostrans[m][i];
            // surr. box crosslinker component m
            double surrboxposm = positions[boxnum*3+m][j];

            double distm = fabs(surrboxposm-cboxposm);
            int absbin = (int)floor(distm/periodlength*(double)numbins);
            if(absbin==3*numbins)
              absbin--;
            int thebin = absbin%numbins;
            int thecol = (int)floor((double)absbin/(double)numbins)+3*m;
            crosslinksperbinrow[thecol][thebin] += 1.0;
          }

          // sort inter-crosslink distances of ROTATED system into bins
          for(int m=0; m<centerboxrotpostrans.NumVectors(); m++)
          {
            //centerbox crosslinker component m
            double cboxposrotm = centerboxrotpostrans[m][i];
            // surrounding box crosslinker position component m
            double surrboxposrotm = rotatedpositions[boxnum*3+m][j];
            int surrboxindex = boxindices[boxnum][m];
            // check if both crosslinkers lie within the shifted and rotated boundary volume
            double deltacbox = fabs(cboxposrotm-(*centershift)(m));
            double deltasbox = fabs(surrboxposrotm-((*centershift)(m)+(surrboxindex-1)*periodlength));
            if(deltacbox<=periodlength/2.0 && deltasbox<=periodlength/2.0) // inside the rotated volume
            {
              // determine bin for rotated fixed system coordinates
              double distm = fabs(surrboxposrotm-cboxposrotm);
              int thebin = (int)floor(distm/periodlength*(double)numbins);
              int thecol = m;
              // only distances [0;H[
              if(thebin<numbins)
                crosslinksperbinrotrow[thecol][thebin] += 1.0;
            }
          }
        }
  return;
}//STATMECH::StatMechManager::DDCorrFunction()

/*------------------------------------------------------------------------------*                                                 |
 | simply counts the number of free, one-bonded, and two-bonded crosslinkers    |
 |                                                        (public) mueller 4/11 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkCount(const std::ostringstream& filename)
{
  if(discret_->Comm().MyPID()==0)
  {
    int free = 0;
    int onebond = 0;
    int twobond = 0;

    for(int i=0; i<crosslinkerbond_->MyLength(); i++)
    {
      // count number of node IDs in crosslinkerbond_ entry i
      int numofnodes = 0;
      for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
        if((*crosslinkerbond_)[j][i]>-0.9)
          numofnodes++;
      // increment according to numofnodes
      switch(numofnodes)
      {
        // free
        case 0:
          free++;
        break;
        // crosslink with one bond
        case 1:
          onebond++;
        break;
        // crosslink with two bonds
        case 2:
          twobond++;
        break;
      }
    }

    // write to file
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream ccount;
    ccount<<free<<"    "<<onebond<<"    "<<twobond<<"    ";
    for(int i=0; i<13; i++)
      ccount<<"    0";
    ccount<<std::endl;
    fputs(ccount.str().c_str(), fp);
    fclose(fp);
  }

  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | output the coverage of crosslinker binding sites (nodes) and the distribution|
 | of bound crosslinkers                                 (private) mueller 5/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkCoverageOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, bool coverageonly)
{
  if(coverageonly)
  {
    Teuchos::RCP<Epetra_Vector> bspotstatustrans = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_));
    CommunicateVector(bspotstatustrans, bspotstatus_,true,false,false,true);

    // check for occupied binding spots
    int partialsum = 0;
    for(int i=0; i<bspotstatustrans->MyLength(); i++)
    {
      if((*filamentnumber_)[bspotcolmap_->LID(bspotrowmap_->GID(i))]==0 && (*bspotstatustrans)[i]>-0.1)
        partialsum++;
      else if((int)(*filamentnumber_)[bspotcolmap_->LID(bspotrowmap_->GID(i))]>0)
        break;
    }
    // sum up processor-specific values and communicate
    int sum = 0;
    discret_->Comm().SumAll(&partialsum,&sum,1);

    if(!discret_->Comm().MyPID())
    {
      // file pointer
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "a");
      std::stringstream coverage;

      coverage << sum << " "<< std::setprecision(5) <<(double)(sum)/(double)(numbspots_) <<std::endl;
      // print to file and close
      fputs(coverage.str().c_str(), fp);
      fclose(fp);
    }
  }
  else
  {
    /* Consider the horizontal filament of a loom setup. Go along this filament and count the
     * number of occupied binding spots. Also, store the spatial distribution of occupied spots
     * for analysis of cluster size, etc.*/
    Epetra_Vector discol(*discret_->DofColMap(),true);
    LINALG::Export(disrow,discol);
    std::map<int, LINALG::Matrix<3, 1> > currentpositions;
    std::map<int, LINALG::Matrix<3, 1> > currentrotations;
    GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);

    if(!discret_->Comm().MyPID())
    {
      // file pointer
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "a");
      std::stringstream coverage;

      for(int i=0; i<bspotstatus_->MyLength(); i++)
      {
        // consider only filament 0 (horizontal filament)
        if((int)(*filamentnumber_)[i]==0 && (*bspotstatus_)[i]>-0.1)
          // store both singly and doubly bound linkers
          coverage<<bspotcolmap_->GID(i)<<"  "<<(int)(*numbond_)[(int)(*bspotstatus_)[i]]<<"  "<< std::scientific << std::setprecision(15) <<(currentpositions.find(i)->second)(0)<<std::endl;
        else if((int)(*filamentnumber_)[i]>0)
          break;
      }
      // print to file and close
      fputs(coverage.str().c_str(), fp);
      fclose(fp);
    }
  }
  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | linker spot counter and check of interfilament orient. (public) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OrientationCorrelation(const Epetra_Vector& disrow, const int &istep)
{
  /* what this does:
   * 1) orientation correlation function for all proximal binding spot pairs (regardless of orientation)
   * 2) count the number of overall possible binding spots for crosslinkers (proximal and correctly oriented)
   * 3) order parameter
   */
  if (DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT"))
  {
    if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
      dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

    Epetra_Vector discol(*discret_->DofColMap(), true);
    LINALG::Export(disrow, discol);

    Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::null;
    if(statmechparams_.get<double>("ILINK",0.0)>0.0)
      bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    GetBindingSpotPositions(discol,bspotpositions, bspotrotations);

    // NODAL TRIAD UPDATE
    Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::null;
    // beam elements only
      bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
      GetBindingSpotTriads(bspotrotations, bspottriadscol);

    // distance and orientation checks
    int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
    // max. distance between two binding spots
    double maxdist = periodlength_->at(0)*sqrt(3.0);
    // max. angle
    double maxangle = M_PI/2.0;
    // minimal and maximal linker search radii
    double rmin = statmechparams_.get<double>("R_LINK",0.0)-statmechparams_.get<double>("DeltaR_LINK",0.0);
    double rmax = statmechparams_.get<double>("R_LINK",0.0)+statmechparams_.get<double>("DeltaR_LINK",0.0);
    // number of overall independent combinations
    int numnodes = discret_->NodeColMap()->NumMyElements();
    int numcombinations = (numnodes*numnodes-numnodes)/2;
    // combinations on each processor
    int combinationsperproc = (int)floor((double)numcombinations/(double)discret_->Comm().NumProc());
    int remainder = numcombinations%combinationsperproc;

    int bindingspots = 0;
    int combicount = 0;
    Teuchos::RCP<Epetra_Vector> anglesrow = Teuchos::rcp(new Epetra_Vector(*ddcorrrowmap_, true));
    // including correlation on the same filament
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsincrow = Teuchos::rcp(new Epetra_MultiVector(*ddcorrrowmap_,2, true));
    // excluding correlation on identical filament
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsexcrow = Teuchos::rcp(new Epetra_MultiVector(*ddcorrrowmap_,2, true));
    // loop over crosslinkermap_ (column map, same for all procs: maps all crosslink molecules)
    for(int mypid=0; mypid<discret_->Comm().NumProc(); mypid++)
    {
      bool quitloop = false;
      if(mypid==discret_->Comm().MyPID())
      {
        bool continueloop = false;
        int appendix = 0;
        if(mypid==discret_->Comm().NumProc()-1)
          appendix = remainder;

        for(int i=0; i<bspotcolmap_->NumMyElements(); i++)
        {
          for(int j=0; j<bspotcolmap_->NumMyElements(); j++)
          {
            // start adding crosslink from here
            if(i==(*startindex_)[2*mypid] && j==(*startindex_)[2*mypid+1])
              continueloop = true;
            // only entries above main diagonal and within limits of designated number of crosslink molecules per processor
            if(j>i && continueloop)
            {
              if(combicount<combinationsperproc+appendix)
              {
                combicount++;

                Epetra_SerialDenseMatrix LID(2,1);
                LID(0,0) = i;
                LID(1,0) = j;

                LINALG::Matrix<3,1> distance;
                for(int k=0; k<(int)distance.M(); k++)
                  distance(k) = (*bspotpositions)[k][i]-(*bspotpositions)[k][j];

                // current distance bin
                int currdistbin = (int)floor(distance.Norm2()/maxdist*numbins);
                // reduce bin if distance == maxdist
                if(currdistbin==numbins)
                  currdistbin--;

                // direction between currently considered two nodes
                LINALG::Matrix<3,1> direction(distance);
                direction.Scale(1.0/direction.Norm2());
                Teuchos::RCP<double> phifil = Teuchos::rcp(new double(0.0));
                bool orientation = CheckOrientation(direction,discol,bspottriadscol,LID,phifil);

                // increment count for that bin
                (*orderparameterbinsincrow)[0][currdistbin] += 1.0;
                // order parameter
                (*orderparameterbinsincrow)[1][currdistbin] += (3.0*cos((*phifil))*cos((*phifil))-1)/2.0;
                if((*filamentnumber_)[i]!=(*filamentnumber_)[j])
                {
                  (*orderparameterbinsexcrow)[0][currdistbin] += 1.0;
                  (*orderparameterbinsexcrow)[1][currdistbin] += (3.0*cos((*phifil))*cos((*phifil))-1)/2.0;
                }
                // proximity check
                if(distance.Norm2()>rmin && distance.Norm2()<rmax)
                {
                  // if angular constraints are met, increase binding spot count
                  if(orientation)
                    bindingspots++;

                  // determine the bin
                  int curranglebin = (int)floor((*phifil)/maxangle*numbins);
                  // in case the distance is exactly periodlength*sqrt(3)
                  if(curranglebin==numbins)
                    curranglebin--;
                  (*anglesrow)[curranglebin] += 1.0;
                }
              }
              else
              {
                quitloop = true;
                break;
              }
            }
          }
          if(quitloop)
            break;
        }
        if(quitloop)
          break;
      }
      else
        continue;
    }
    // Export
    // add up binding spots
    int bspotsglob = 0;
    discret_->Comm().SumAll(&bindingspots,&bspotsglob,1);

    Teuchos::RCP<Epetra_Vector> anglescol = Teuchos::rcp(new Epetra_Vector(*ddcorrcolmap_, true));
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsinccol = Teuchos::rcp(new Epetra_MultiVector(*ddcorrcolmap_,2, true));
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsexccol = Teuchos::rcp(new Epetra_MultiVector(*ddcorrcolmap_,2, true));
    CommunicateVector(anglesrow,anglescol, false, true,false);
    CommunicateMultiVector(orderparameterbinsincrow, orderparameterbinsinccol,false,true,false);
    CommunicateMultiVector(orderparameterbinsexcrow, orderparameterbinsexccol,false,true,false);

    // Add the processor-specific data up
    std::vector<int> angles(numbins, 0);
    std::vector<std::vector<double> > orderparameterinc(numbins, std::vector<double>(2,0.0));
    std::vector<std::vector<double> > orderparameterexc(numbins, std::vector<double>(2,0.0));
    for(int i=0; i<numbins; i++)
      for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
      {
        angles[i] += (int)(*anglescol)[pid*numbins+i];
        orderparameterinc[i][0] += (*orderparameterbinsinccol)[0][pid*numbins+i];
        orderparameterinc[i][1] += (*orderparameterbinsinccol)[1][pid*numbins+i];
        orderparameterexc[i][0] += (*orderparameterbinsexccol)[0][pid*numbins+i];
        orderparameterexc[i][1] += (*orderparameterbinsexccol)[1][pid*numbins+i];
      }

    // average values
    for(int i=0; i<numbins; i++)
      if(orderparameterinc[i][0]>0.0) // i.e. >0
      {
        orderparameterinc[i][1] /= orderparameterinc[i][0];
        orderparameterexc[i][1] /= orderparameterexc[i][0];
      }
      else
      {
        orderparameterinc[i][1] = -99.0;
        orderparameterexc[i][1] = -99.0;
      }

    // write data to file
    if(!discret_->Comm().MyPID())
    {
      std::ostringstream orientfilename;
      orientfilename << outputrootpath_ << "/StatMechOutput/LinkerSpotsOrCorr_"<<std::setw(6) << std::setfill('0') << istep <<".dat";

      FILE* fp = NULL;
      fp = fopen(orientfilename.str().c_str(), "w");
      std::stringstream histogram;

      histogram<<bspotsglob<<"    "<<-99<<"    "<<-99<<std::endl;
      for(int i=0; i<numbins; i++)
        histogram<<i+1<<"    "<<angles[i]<<"    "<<std::setprecision(12)<<orderparameterinc[i][1]<<"    "<<orderparameterinc[i][0]<<"    "<<orderparameterexc[i][1]<<"    "<<orderparameterexc[i][0]<<std::endl;
      //write content into file and close it
      fputs(histogram.str().c_str(), fp);
      fclose(fp);
    }
  }
}//NumLinkerSpotsAndOrientation()

/*------------------------------------------------------------------------------*                                                 |
 | computes the mesh size of the network dep. on distance to cog                |
 |                                                        (public) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ComputeLocalMeshSize(const Epetra_Vector& disrow, LINALG::Matrix<3,1>& centershift, const int &istep)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  double periodlength = periodlength_->at(0);
  double maxdist = periodlength/2.0*sqrt(3.0);
  int numbins = statmechparams_.get<int>("HISTOGRAMBINS",1);
  // center of gravity
  LINALG::Matrix<3,1> cog;
  for(int i=0; i<(int)cog.M(); i++)
    cog(i) = centershift(i) + periodlength/2.0;

  // 1. set up of a row map vector containing the global node positions
  std::vector<LINALG::Matrix<3,1> > xrelnodes;
  xrelnodes.clear();
  for (int i=0; i<discret_->NumMyRowNodes(); i++)
  {
    //get pointer at a node
    const DRT::Node* node = discret_->lRowNode(i);
    //get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = discret_->Dof(node);
    // global position and shift according to new centershift
    LINALG::Matrix<3, 1> xglob;

    // shift according to given center (cog)
    for(int j=0; j<(int)xglob.M(); j++)
    {
      // get the global node position
      xglob(j) = node->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode[j])];
      // shift the j-th component according to new center
      if (xglob(j) > periodlength+centershift(j))
        xglob(j) -= periodlength;
      if (xglob(j) < centershift(j))
        xglob(j) += periodlength;
    }
    xglob -= cog;
    xrelnodes.push_back(xglob);
  }

  Epetra_Vector fillengthrow(*ddcorrrowmap_, true);
  // calculate sqrt(DV/DL)
  for(int i=1; i<discret_->NumMyRowNodes(); i++)
  {
    // column map LID of the two nodes in question
    int gid0 = discret_->NodeRowMap()->GID(i-1);
    int gid1 = discret_->NodeRowMap()->GID(i);
    int collid0 = discret_->NodeColMap()->LID(gid0);
    int collid1 = discret_->NodeColMap()->LID(gid1);

    // make sure both nodes lie on the same filament
    if((*filamentnumber_)[collid0] == (*filamentnumber_)[collid1])
    {
      // calculate node bins
      int index0 = i-1;
      int index1 = i;
      int bin0 = (int)floor(xrelnodes[i-1].Norm2()/maxdist*numbins);
      int bin1 = (int)floor(xrelnodes[i].Norm2()/maxdist*numbins);
      if(bin0 == numbins)
        bin0--;
      if(bin1 == numbins)
        bin1--;

      // case: the two nodes lie within different bins: add length segments binwise
      if(bin0 != bin1)
      {
        // switch in order to follow convention
        if(bin1<bin0)
        {
          index0 = i;
          index1 = i-1;
          bin0 = (int)floor(xrelnodes[i].Norm2()/maxdist*numbins);
          bin1 = (int)floor(xrelnodes[i-1].Norm2()/maxdist*numbins);
        }

        // check, if element is broken and unshift in order to obtain correct directional vector
        // calculate directional vector between nodes
        LINALG::Matrix<3,1> unshift = xrelnodes[index1];

        for(int j=0; j<(int)xrelnodes[index0].M(); j++)
        {
          // check for periodic boundary shift and correct accordingly
          if (fabs(xrelnodes[index1](j) - periodlength - xrelnodes[index0](j)) < fabs(xrelnodes[index1](j) - xrelnodes[index0](j)))
            unshift(j) -= periodlength;
          else if (fabs(xrelnodes[index1](j) + periodlength - xrelnodes[index0](j)) < fabs(xrelnodes[index1](j) - xrelnodes[index0](j)))
            unshift(j) += periodlength;
        }


        LINALG::Matrix<3,1> dirvec = unshift;
        dirvec -= xrelnodes[index0];
        // directional unit vector
        dirvec.Scale(1.0/dirvec.Norm2());

        // number of bin intersections
        int numisecs = abs(bin1-bin0);
        // starting point
        LINALG::Matrix<3,1> xstartj = xrelnodes[index0];

        /*-----------------------------------
         * short example:
         *    bin:   16    17    18    19
         *         |  o--|-----|-----|--o  |
         * abslim: 16    17    18    19    20
         * binlim:       0     1     2
          -----------------------------------*/
        //std::cout<<"numisecs = "<<numisecs<<std::endl;
        for(int j=0; j<numisecs; j++)
        {
          // all segments except last segment
          if(j<=numisecs-1)
          {
            int currbin = bin0+j;
            // j-th bin limit
            double nextbinlimit = (double)(currbin+1)/(double)numbins*maxdist;

            // line parameter mu for intersection of line with spherical shell:
            // note: we take the solution that is (may be?) >=0 since we chose the directional vector accordingly

            //std::cout<<"  j="<<j<<", xstartj: "<<xstartj(0)<<" "<<xstartj(1)<<" "<<xstartj(2)<<std::endl;
            double a = xstartj.Dot(dirvec);
            double b = sqrt(a*a-xstartj.Norm2()*xstartj.Norm2()+nextbinlimit*nextbinlimit);
            double mu = -a + b;

            //std::cout<<"  mu("<<j<<") = "<<mu<<std::endl;

            // intersection of the line with the spherical shell
            LINALG::Matrix<3,1> isection = dirvec;
            isection.Scale(mu);
            isection += xstartj;
            // calculate element (filament) segment added to the respective bin
            LINALG::Matrix<3,1> segment = isection;
            segment -= xstartj;
            fillengthrow[currbin] += segment.Norm2();
            // store current intersection in case there is more than one
            xstartj = isection;
            //if(numisecs>0)
            //  std::cout<<"    j="<<j+1<<", previs: "<<xstartj(0)<<" "<<xstartj(1)<<" "<<xstartj(2)<<std::endl;
          }
          // last segment
          if(j==numisecs-1)
          {
            //std::cout<<"  xfinale: "<<xstartj(0)<<" "<<xstartj(1)<<" "<<xstartj(2)<<std::endl;
            LINALG::Matrix<3,1> segment = xrelnodes[index1];
            segment -= xstartj;
            fillengthrow[bin1] += segment.Norm2();
            xstartj.Clear();
          }
        }
      }
      else // just add the entire element length
      {
        LINALG::Matrix<3,1> elength = xrelnodes[i];
        elength -= xrelnodes[i-1];
        fillengthrow[bin0] += elength.Norm2();
      }
    }
  }

  Epetra_Vector fillengthcol(*ddcorrcolmap_,true);
  Epetra_Import importer(*ddcorrcolmap_, *ddcorrrowmap_);
  fillengthcol.Import(fillengthrow,importer,Insert);
  // Add the processor-specific data up
  std::vector<double> fillength(numbins, 0.0);
  for(int i=0; i<numbins; i++)
    for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
      fillength[i] += fillengthcol[pid*numbins+i];

  // Proc 0 section
  if(discret_->Comm().MyPID()==0)
  {
    std::ostringstream filename;
    filename << outputrootpath_ << "/StatMechOutput/LocalMeshSize_"<<std::setw(6) << std::setfill('0') << istep <<".dat";

    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");
    std::stringstream histogram;

    for(int i=0; i<numbins; i++)
      histogram<<i+1<<"    "<<std::setprecision(12)<<fillength[i]<<std::endl;

    //write content into file and close it
    fputs(histogram.str().c_str(), fp);
    fclose(fp);
  }
}

/*------------------------------------------------------------------------------*
 | Output of linker positions and bond status             (public) mueller 03/14|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputLinkerPositions(std::ostringstream& filename)
{
  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    std::stringstream filecontent;
    fp = fopen(filename.str().c_str(), "a");

    for(int i=0; i<numbond_->MyLength(); i++)
      filecontent<<std::scientific<<std::setprecision(8)<<(*crosslinkerpositions_)[0][i]<<"\t"<<(*crosslinkerpositions_)[1][i]<<"\t"<<(*crosslinkerpositions_)[2][i]<<"\t"<<(*numbondconv_) [i]<<"\t"<<(*numbond_)[i]<<std::endl;
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }
  return;
}

/*------------------------------------------------------------------------------*
 |Output of change of inclusive angles                                          |
 |between filaments                                     (public) mukherjee 03/15|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputDeltaTheta(std::ostringstream& filename)
{

    FILE* fp = NULL;
    std::stringstream filecontent;
    fp = fopen(filename.str().c_str(), "a");
    for(int i=0; i<discret_->ElementRowMap()->NumMyElements(); i++)
    {
        DRT::Element* ele = discret_->lRowElement(i);
        const DRT::ElementType &eot = ele->ElementType();
        if(eot==DRT::ELEMENTS::Truss3Type::Instance())
        {
          DRT::ELEMENTS::Truss3* crossele = dynamic_cast<DRT::ELEMENTS::Truss3*> (ele);
          double DeltaTheta= crossele->GetDeltaTheta();
          filecontent<<std::scientific<<std::setprecision(8)<<DeltaTheta<<std::endl;
        }
        else if(eot==DRT::ELEMENTS::Beam3Type::Instance())
        {
          DRT::ELEMENTS::Beam3* crossele = dynamic_cast<DRT::ELEMENTS::Beam3*> (ele);
          if(crossele==NULL)
            return;
          LINALG::Matrix<3,1> Tcurr1_alt(true);
          LINALG::Matrix<3,1> Tcurr2_alt(true);
          LINALG::Matrix<3,1> Tref1_alt(true);
          LINALG::Matrix<3,1> Tref2_alt(true);
          crossele->TcurrBeam3r(Tcurr1_alt,Tcurr2_alt);
          crossele->TrefBeam3r(Tref1_alt,Tref2_alt);
          double Phi_ref_alt= GetTheta(Tref1_alt,Tref2_alt);
          double Phi_curr_alt= GetTheta(Tcurr1_alt,Tcurr2_alt);
          double DeltaTheta= abs(Phi_ref_alt-Phi_curr_alt);
          filecontent<<std::scientific<<std::setprecision(8)<<DeltaTheta<<std::endl;
        }
    }
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);

  return;
}

/*------------------------------------------------------------------------------*
 | Output of linker length for spring linkers           (public) mukherjee 03/15|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputLinkerLength(std::ostringstream& filename)
{
    FILE* fp = NULL;
    std::stringstream filecontent;
    fp = fopen(filename.str().c_str(), "a");
    for(int i=0; i<discret_->ElementRowMap()->NumMyElements(); i++)
    {
        DRT::Element* ele = discret_->lRowElement(i);
        const DRT::ElementType &eot = ele->ElementType();
        if(eot==DRT::ELEMENTS::Truss3Type::Instance()) // i.e. it's a spring linker
        {
          DRT::ELEMENTS::Truss3* crossele = dynamic_cast<DRT::ELEMENTS::Truss3*> (ele);
          double lcurrCrosslink= crossele->Lcurr();
          filecontent<<std::scientific<<std::setprecision(8)<<lcurrCrosslink<<std::endl;
        }
        if(eot==DRT::ELEMENTS::Beam3Type::Instance()) // i.e. it's a beam linker
        {
          DRT::ELEMENTS::Beam3* crossele = dynamic_cast<DRT::ELEMENTS::Beam3*> (ele);
          double lcurrCrosslink= crossele->Lcurr();
          filecontent<<std::scientific<<std::setprecision(8)<<lcurrCrosslink<<std::endl;
        }
    }
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);

  return;
}

/*------------------------------------------------------------------------------*
 | Output of relative motion between linker nodes         (public) mueller 05/13|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputSlidingMotion(const Epetra_Vector& disrow, std::ostringstream& filename)
{
  // initialize on first call
  if(linkernodepairs_==Teuchos::null)
  {
    linkernodepairs_ = Teuchos::rcp(new std::vector<std::vector<int> >);
    for(int i=0; i<numbond_->MyLength(); i++)
    {
      if((*numbond_)[i]>1.9)
      {
        std::vector<int> pair(2, -1);
        pair[0] = (int)(*crosslinkerbond_)[0][i];
        pair[1] = (int)(*crosslinkerbond_)[1][i];
        linkernodepairs_->push_back(pair);
      }
    }
  }

  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);
  std::map<int, LINALG::Matrix<3, 1> > currentpositions;
  std::map<int, LINALG::Matrix<3, 1> > currentrotations;
  GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);

  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    std::stringstream filecontent;
    fp = fopen(filename.str().c_str(), "a");

    for(int i=0; i<(int)linkernodepairs_->size(); i++)
    {
      int bspotlid0 = bspotcolmap_->LID((*linkernodepairs_)[i][0]);
      int bspotlid1 = bspotcolmap_->LID((*linkernodepairs_)[i][1]);
      std::map<int, LINALG::Matrix<3,1> >::const_iterator posbspot0 = currentpositions.find(bspotlid0);
      std::map<int, LINALG::Matrix<3,1> >::const_iterator posbspot1 = currentpositions.find(bspotlid1);

      // write binding spot position
      filecontent<<std::setprecision(8)<<(posbspot0->second)(0)<<" "<<(posbspot0->second)(1)<<" "<<(posbspot0->second)(2)<<" ";
      filecontent<<std::setprecision(8)<<(posbspot1->second)(0)<<" "<<(posbspot1->second)(1)<<" "<<(posbspot1->second)(2)<<" ";
      // write binding spot status
      filecontent<<(*bspotstatus_)[bspotlid0]<<" "<<(*bspotstatus_)[bspotlid1]<<" ";
      // write linker status of the linker attached to above's binding spots
      if((*bspotstatus_)[bspotlid0]>-0.9)
        filecontent<<(*numbond_)[(*bspotstatus_)[bspotlid0]];
      if((*bspotstatus_)[bspotlid1]>-0.9)
        filecontent<<" "<<(*numbond_)[(*bspotstatus_)[bspotlid1]];
      filecontent<<std::endl;
    }
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }
  return;
}

/*------------------------------------------------------------------------------*
 | Structure COG & inertia tensor output                  (public) mueller 05/13|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::StructureCOGInertiaTensorOutput(const int&           istep,
                                                                const double&        time,
                                                                const Epetra_Vector& disrow,
                                                                std::ostringstream&  filename)
{
  // export row displacement to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);
  std::map<int, LINALG::Matrix<3, 1> > currentpositions;
  std::map<int, LINALG::Matrix<3, 1> > currentrotations;
  GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);

  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    std::stringstream filecontent;
    fp = fopen(filename.str().c_str(), "a");

    filecontent << std::scientific << std::setprecision(15) << istep<<"  "<<time<<"  "<<discret_->NumMyColNodes()<<std::endl;

    // calculate center of gravity of the structure
    LINALG::Matrix<3,1> COG(true);
    for(std::map< int,LINALG::Matrix<3,1> >::const_iterator it = currentpositions.begin(); it!=currentpositions.end(); it++)
      COG += it->second;
    COG.Scale(1.0/(double)discret_->NumMyColNodes());
//    std::cout<<"COG =\n"<<COG<<std::endl;
//    filecontent << std::scientific << std::setprecision(15) << COG(0)<<"  "<<COG(1)<<"  "<<COG(2)<<std::endl;

    // calculate relative position to COG
    for(int i=0; i<discret_->NumMyColNodes(); i++)
    {
      std::map< int,LINALG::Matrix<3,1> >::iterator posi = currentpositions.find(i);
      posi->second -= COG;
    }
    // calculate inertial tensor
    LINALG::Matrix<3,3> I_ij(true);
    LINALG::Matrix<3,3> delta(true);
    delta(0,0) = 1.0;
    delta(1,1) = 1.0;
    delta(2,2) = 1.0;
    for(std::map< int,LINALG::Matrix<3,1> >::const_iterator it = currentpositions.begin(); it!=currentpositions.end(); it++)
      for(int i=0; i<(int)I_ij.M(); i++)
        for(int j=0; j<(int)I_ij.N(); j++)
          I_ij(i,j) += (it->second).Dot(it->second)*delta(i,j) - (it->second)(i)*(it->second)(j);
    // Scale in order to avoid large numbers
    I_ij.Scale(1.0/double(discret_->NumMyColNodes()));

    for(int i=0; i<(int)I_ij.M(); i++)
    {
      for(int j=0; j<(int)I_ij.N(); j++)
        filecontent << std::scientific << std::setprecision(15) << I_ij(i,j)<<"  ";
      filecontent<<std::endl;
    }


    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }
}

/*------------------------------------------------------------------------------*                                                 |
| distribution of spherical coordinates                  (public) mueller 12/10|
*------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::SphericalCoordsDistribution(const Epetra_Vector& disrow,
                                                  Epetra_Vector& phibinsrow,
                                                  Epetra_Vector& thetabinsrow,
                                                  Epetra_Vector& costhetabinsrow)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
  double periodlength = periodlength_->at(0);

  for(int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* element = discret_->lRowElement(i);
    // consider filament elements only
    if(element->Id()<basisnodes_)
    {
      int gid0 = element->Nodes()[0]->Id();
      int gid1 = element->Nodes()[1]->Id();
      int lid0 = discret_->NodeRowMap()->LID(gid0);
      int lid1 = discret_->NodeRowMap()->LID(gid1);
      DRT::Node* node0 = discret_->lRowNode(lid0);
      DRT::Node* node1 = discret_->lRowNode(lid1);

      // calculate directional vector between nodes
      double dirlength = 0.0;
      Epetra_SerialDenseMatrix dirvec(3,1);
      for(int dof=0; dof<dirvec.M(); dof++)
      {
        int dofgid0 = discret_->Dof(node0)[dof];
        int dofgid1 = discret_->Dof(node1)[dof];
        double poscomponent0 = node0->X()[dof]+disrow[discret_->DofRowMap()->LID(dofgid0)];
        double poscomponent1 = node1->X()[dof]+disrow[discret_->DofRowMap()->LID(dofgid1)];
        // check for periodic boundary shift and correct accordingly
        if (fabs(poscomponent1 - periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
          poscomponent1 -= periodlength;
        else if (fabs(poscomponent1 + periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
          poscomponent1 += periodlength;

        dirvec(dof,0) = poscomponent1-poscomponent0;
        dirlength += dirvec(dof,0)*dirvec(dof,0);
      }
      dirlength = sqrt(dirlength);
      // normed directional vector
      dirvec.Scale(1.0/dirlength);

      Epetra_SerialDenseMatrix dirvecrot(3,1,true);
      trafo_->Multiply(false,dirvec,dirvecrot);

      // transform into spherical coordinates (phi E [-pi;pi], theta E [0; pi]) and sort into appropriate bin
      double phi = atan2(dirvecrot(1,0),dirvecrot(0,0)) + M_PI;
      double theta = acos(dirvecrot(2,0));

      int phibin = (int)floor(phi/(2*M_PI)*numbins);
      int thetabin = (int)floor(theta/M_PI*numbins);
      int costhetabin = (int)floor((cos(theta)+1.0)/2.0*numbins);
      if(phibin == numbins)
        phibin--;
      if(thetabin == numbins)
        thetabin--;
      if(costhetabin == numbins)
        costhetabin--;
      if(phibin<0 || thetabin<0)
        dserror("bin smaller zero");
      phibinsrow[phibin] += 1.0;
      thetabinsrow[thetabin] += 1.0;
      costhetabinsrow[costhetabin] += 1.0;
    }
  }
  return;
}//SphericalCoordsDistribution()

/*------------------------------------------------------------------------------*
 | radial crosslinker density distribution                (public) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::RadialDensityDistribution(Epetra_Vector& radialdistancesrow, LINALG::Matrix<3,1>& centershift)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  // simpler version taking into account only the original boundary volume
  double periodlength = periodlength_->at(0);
  double maxdistance = periodlength/2.0*sqrt(3.0);
  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);

  LINALG::Matrix<3,1> cog;
  for(int i=0; i<(int)cog.M(); i++)
    cog(i) = centershift(i) + periodlength/2.0;

  // Export to transfer map format
  Epetra_MultiVector crosslinkerbondtrans(*transfermap_, 2, true);
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, 3, true);
  Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
  Epetra_Export crosslinkimporter(*crosslinkermap_, *transfermap_);
  // note: it seems to be necessary to clear all vectors other than on Proc 0.
  // otherwise, i.e. using method "Insert" on Export, from time to time (no clear pattern has emerged so far)
  // incorrect data is written to transfer vectors. Odd!
  if(discret_->Comm().MyPID()!=0)
  {
    crosslinkerbond_->PutScalar(0.0);
    crosslinkerpositions_->PutScalar(0.0);
  }
  //Export
  crosslinkerbondtrans.Export(*crosslinkerbond_, crosslinkexporter, Add);
  crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);
  // Reimport to create pre-Export status on all procs
  crosslinkerbond_->Import(crosslinkerbondtrans, crosslinkimporter, Insert);
  crosslinkerpositions_->Import(crosslinkerpositionstrans, crosslinkimporter, Insert);

  // calculate distance of crosslinkers to center of gravity
  for(int i=0; i<crosslinkerbondtrans.MyLength(); i++)
    if(crosslinkerbondtrans[0][i]>-0.9 && crosslinkerbondtrans[1][i]>-0.9)
    {
      LINALG::Matrix<3,1> distance;
      // shift coordinates according to new boundary box center
      for(int j=0; j<crosslinkerpositionstrans.NumVectors(); j++)
      {
        distance(j) = crosslinkerpositionstrans[j][i];
        if (distance(j) > periodlength+centershift(j))
          distance(j) -= periodlength;
        if (distance(j) < centershift(j))
          distance(j) += periodlength;
      }
      distance -= cog;
      // determine histogram bin for current crosslinker element
      int currbin = (int)floor(distance.Norm2()/maxdistance * numbins);
      if(currbin==numbins)
        currbin--;
      radialdistancesrow[currbin] += 1.0;
    }
  return;
}//RadialDensityDistribution()

/*----------------------------------------------------------------------*
 | filament orientations and output               (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::FilamentOrientations(const Epetra_Vector& discol, std::vector<LINALG::Matrix<3,1> >* normedvectors, const std::ostringstream& filename, bool fileoutput)
{
  /* Output of filament element orientations (Proc 0 only):
   * format: filamentnumber    d_x  d_y  d_z
   */

  if(discret_->Comm().MyPID()==0)
  {
    if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2)  || periodlength_->at(0) <= 0.0)
      dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

    double periodlength = periodlength_->at(0);

    FILE* fp = NULL;
    if(fileoutput)
      fp = fopen(filename.str().c_str(), "w");
    std::stringstream fileleorientation;

    // get filament number conditions
    std::vector<DRT::Condition*> filaments(0);
    discret_->GetCondition("FilamentNumber", filaments);

    for(int fil=0; fil<(int)filaments.size(); fil++)
    {
      // get next filament
      DRT::Condition* currfilament = filaments[fil];
      for(int node=1; node<(int)currfilament->Nodes()->size(); node++)
      {
        // obtain column map LIDs
        int gid0 = currfilament->Nodes()->at(node-1);
        int gid1 = currfilament->Nodes()->at(node);
        int nodelid0 = discret_->NodeColMap()->LID(gid0);
        int nodelid1 = discret_->NodeColMap()->LID(gid1);
        DRT::Node* node0 = discret_->lColNode(nodelid0);
        DRT::Node* node1 = discret_->lColNode(nodelid1);

        // calculate directional vector between nodes
        LINALG::Matrix<3, 1> dirvec;
        for(int dof=0; dof<3; dof++)
        {
          int dofgid0 = discret_->Dof(node0)[dof];
          int dofgid1 = discret_->Dof(node1)[dof];
          double poscomponent0 = node0->X()[dof]+discol[discret_->DofColMap()->LID(dofgid0)];
          double poscomponent1 = node1->X()[dof]+discol[discret_->DofColMap()->LID(dofgid1)];
          // check for periodic boundary shift and correct accordingly
          if (fabs(poscomponent1 - periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
            poscomponent1 -= periodlength;
          else if (fabs(poscomponent1 + periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
            poscomponent1 += periodlength;

          dirvec(dof) = poscomponent1-poscomponent0;
        }
        // normed vector
        dirvec.Scale(1/dirvec.Norm2());

        // add element directional vectors up
        for(int i=0; i<(int)normedvectors->size(); i++)
        {
          // base vector
          LINALG::Matrix<3,1> ei;
          LINALG::Matrix<3,1> vi = dirvec;

          ei.Clear();
          ei(i) = 1.0;

          if(acos(vi.Dot(ei))>(M_PI/2.0))
            vi.Scale(-1.0);
          (*normedvectors)[i] += vi;
        }

        // write normed directional vector to stream
        if(fileoutput)
          fileleorientation<<fil<<"    "<<std::setprecision(12)<<dirvec(0)<<" "<<dirvec(1)<<" "<<dirvec(2)<<std::endl;
      }
    }
    if(fileoutput)
    {
      fputs(fileleorientation.str().c_str(), fp);
      fclose(fp);
    }
  }
}// StatMechManager:FilamentOrientations()

/*----------------------------------------------------------------------*
 | output of free filament length segments       (private) mueller 10/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputFreeFilamentLength(const Epetra_Vector&      disrow,
                                                         const std::ostringstream& filename)
{
  Epetra_Vector discol(*(discret_->DofColMap()),true);
  LINALG::Export(disrow,discol);
  Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_, 3, true));
  GetBindingSpotPositions(discol, bspotpositions, Teuchos::null);

  // Proc 0 only
  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");

    std::stringstream freefillength;

    double freelength = 0.0;

    std::vector<double> endtoend(6,0.0);
    // initialize with first bspot position of 1st filament
    for(int i=0; i<bspotpositions->NumVectors(); i++)
      endtoend[i] = (*bspotpositions)[i][0];
    //double summedlength = 0.0;
    for(int n=1; n<bspotpositions->MyLength(); n++)
    {
      // standard linkers: bspot is equal to node
      int bspotn = n-1;
      int bspotnp = n;
      // interpolated binding sites require additional measures in order to acquire the required node IDs
      if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
      {
        // take the first node (0) associated with the binding spots n-1 and n. Due to the design of the interpolated binding sites,
        // it is impossible that nodes 0 and 1 lie on different filaments.
        bspotn = discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][n-1]);
        bspotnp = discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][n]);
      }
      // binding spots have to be on the same filament
      if((*filamentnumber_)[bspotn] == (*filamentnumber_)[bspotnp])
      {
        // first binding spot pair of a) a filament or b) after a doubly bound crosslinker
        if(freelength>0.0 && crosslinkermap_->LID((int)(*bspotstatus_)[bspotn])>-1) // = crosslinker attached
        {
          if((*crosslink2element_)[(*bspotstatus_)[bspotn]]>1.9) //  = crosslinker element
          {
            // calculate end-to-end distance
            std::vector<double> tmppos(3,0.0);
            for(int i=0; i<bspotpositions->NumVectors(); i++)
            {
              tmppos[i] = (*bspotpositions)[i][bspotn];
              endtoend[i+3] = (*bspotpositions)[i][bspotn];
            }
            UnshiftPositions(endtoend,2,false);
            double e2elength = sqrt((endtoend[3]-endtoend[0])*(endtoend[3]-endtoend[0]) +
                                    (endtoend[4]-endtoend[1])*(endtoend[4]-endtoend[1]) +
                                    (endtoend[5]-endtoend[2])*(endtoend[5]-endtoend[2]));
            // prepare for next end to end distance measurement
            for(int i=0; i<bspotpositions->NumVectors(); i++)
              endtoend[i] = (*bspotpositions)[i][bspotn];
            freefillength<<std::setprecision(6)<<freelength<<"\t"<<e2elength<<std::endl;
            //summedlength += freelength;
            // after every crosslinker element, a new free filament segment begins. Hence reset.
            freelength = 0.0;
          }
        }

        std::vector<double> bspotpos(6,0.0);
        for(int i=0; i<bspotpositions->NumVectors(); i++)
        {
          bspotpos[i] = (*bspotpositions)[i][bspotn];
          bspotpos[i+3] = (*bspotpositions)[i][bspotnp];
        }
        UnshiftPositions(bspotpos, 2, false);
        double addedlength = sqrt((bspotpos[3]-bspotpos[0])*(bspotpos[3]-bspotpos[0]) +
                                  (bspotpos[4]-bspotpos[1])*(bspotpos[4]-bspotpos[1]) +
                                  (bspotpos[5]-bspotpos[2])*(bspotpos[5]-bspotpos[2]));
        freelength += addedlength;
      }
      else // write whatever is left as free filament length before hopping to the nex filament
      {
        for(int i=0; i<bspotpositions->NumVectors(); i++)
          endtoend[i+3] = (*bspotpositions)[i][bspotn];
        UnshiftPositions(endtoend,2,false);
        double e2elength = sqrt((endtoend[3]-endtoend[0])*(endtoend[3]-endtoend[0]) +
                                (endtoend[4]-endtoend[1])*(endtoend[4]-endtoend[1]) +
                                (endtoend[5]-endtoend[2])*(endtoend[5]-endtoend[2]));
        // prepare for next filament
        for(int i=0; i<bspotpositions->NumVectors(); i++)
          endtoend[i] = (*bspotpositions)[i][bspotnp];
        //summedlength += freelength;
        //summedlength = 0.0;
        freefillength<<std::setprecision(6)<<freelength<<"\t"<<e2elength<<std::endl;
        // reset in order to start fresh with the next filament
        freelength = 0.0;
      }
    }
    fputs(freefillength.str().c_str(), fp);
    fclose(fp);
  }

  return;
}
