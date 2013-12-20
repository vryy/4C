/*!

\file turbulence_hit_initial_field.cpp

\brief routines to initialize homogeneous isotropic turbulence simulations

References are

<pre>

</pre>

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*/

#include <complex>
#include <math.h>

#ifdef HAVE_FFTW
#include"fftw3.h"
#endif

#include "turbulence_hit_initial_field.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

namespace FLD
{

/*--------------------------------------------------------------*
 | constructor                                  rasthofer 04/13 |
 *--------------------------------------------------------------*/
HomIsoTurbInitialField::HomIsoTurbInitialField(
        FluidImplicitTimeInt& timeint,
        const INPAR::FLUID::InitialField initfield)
        :
        discret_(timeint.discret_),
        velnp_(timeint.velnp_),
        veln_(timeint.veln_),
        type_(initfield)
{
  // determine number of modes
  // number of modes equal to number of elements in one spatial direction
  // this does not yield the correct value
  // nummodes_ = (int) pow((double) discret_->NumGlobalElements(),1.0/3.0);
  switch (discret_->NumGlobalElements())
  {
    case 512:
    {
      nummodes_ = 8;
      break;
    }
    case 32768:
    {
      nummodes_ = 32;
      break;
    }
    case 110592:
    {
      nummodes_ = 48;
      break;
    }
    case 262144:
    {
      nummodes_ = 64;
      break;
    }
    default:
    {
      dserror("Set problem size! %i", discret_->NumGlobalElements());
      break;
    }
  }

  //-------------------------------------------------
  // create set of node coordinates
  //-------------------------------------------------

  // the criterion allows differences in coordinates by 1e-9
  std::set<double,LineSortCriterion> coords;
  // loop all nodes and store x1-coordinate
  for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
  {
    DRT::Node* node = discret_->lRowNode(inode);
    if ((node->X()[1]<2e-9 && node->X()[1]>-2e-9) and (node->X()[2]<2e-9 && node->X()[2]>-2e-9))
    coords.insert(node->X()[0]);
  }
  // communicate coordinates to all procs via round Robin loop
  {
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

    // communicate coordinates
    for (int np=0;np<numprocs;++np)
    {
      DRT::PackBuffer data;

      for (std::set<double,LineSortCriterion>::iterator x1line = coords.begin();
           x1line != coords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(data,*x1line);
      }
      data.StartPacking();
      for (std::set<double,LineSortCriterion>::iterator x1line = coords.begin();
           x1line != coords.end();
           ++x1line)
      {
        DRT::ParObject::AddtoPack(data,*x1line);
      }
      std::swap( sblock, data() );

#ifdef PARALLEL
      MPI_Request request;
      int tag = myrank;

      int frompid = myrank;
      int topid  = (myrank+1)%numprocs;

      int length=sblock.size();

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

      // unpack received block into set of all coordinates
      {
        std::vector<double> coordsvec;

        coordsvec.clear();

        std::vector<char>::size_type index = 0;
        while (index < rblock.size())
        {
          double onecoord;
          DRT::ParObject::ExtractfromPack(index,rblock,onecoord);
          coords.insert(onecoord);
        }
      }
    }
  }
  // push coordinates in vectors
  {
    coordinates_ = Teuchos::rcp(new std::vector<double> );

    for(std::set<double,LineSortCriterion>::iterator coord1 = coords.begin();
        coord1 != coords.end();
        ++coord1)
    {
      coordinates_->push_back(*coord1);
    }
  }

  //-------------------------------------------------
  // set-up experimental data
  //-------------------------------------------------
  // non-dimensionalize and store experimental data

  if (type_ == INPAR::FLUID::initialfield_hit_comte_bellot_corrsin)
    // of Comte-Bellot-Corrsin experiment
    PrepareExparimentalData();

  return;
}


// dummy function to check the functionality of the fftw library
#if 0
void HomIsoTurbInitialField::CalculateInitialField()
{
  Teuchos::RCP<Teuchos::Array <std::complex<double> > > u1_hat = Teuchos::rcp( new Teuchos::Array<std::complex<double> >(nummodes_*nummodes_*(nummodes_/2+1)));
  Teuchos::RCP<Teuchos::Array <double> > u1 = Teuchos::rcp( new Teuchos::Array<double>(nummodes_*nummodes_*nummodes_));

  for (int i = 0; i < nummodes_; i++)
  {
    for (int j = 0; j< nummodes_; j++)
    {
      for (int l = 0; l< nummodes_; l++){
      int pos = l+ nummodes_ * j + nummodes_ * nummodes_ * i;
      double myi = (double) i;
      double myj = (double) j;
      double myl = (double) l;
      (*u1)[pos] = (myi+1) * (myj+1) * (myl+1);
      }
    }
  }

  // output on screen
  std::cout << "u1_hat " << u1->size() << std::endl;
  for (int i=0; i< u1->size(); i++)
    std::cout << (*u1)[i] << std::endl;

  // set-up
  fftw_plan fft = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_,
                                       &((*u1)[0]),
                                       (reinterpret_cast<fftw_complex*>(&((*u1_hat)[0]))),
                                       FFTW_ESTIMATE);
  // fft
  fftw_execute(fft);

  for (int i=0; i< u1_hat->size(); i++)
      (*u1_hat)[i] /= nummodes_*nummodes_*nummodes_;

  std::cout << "u1_hat " << u1_hat->size() << std::endl;
  for (int i=0; i< u1_hat->size(); i++)
    std::cout << "real  " << real((*u1_hat)[i]) << "    imag  " << imag((*u1_hat)[i]) << std::endl;

  // backward
  Teuchos::RCP<Teuchos::Array <double> > u1_back = Teuchos::rcp( new Teuchos::Array<double>(nummodes_*nummodes_*nummodes_));
  fft = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_, //instead of (nummodes_/2+1): this is important to have here
                             (reinterpret_cast<fftw_complex*>(&((*u1_hat)[0]))),
                             &((*u1_back)[0]),
                             FFTW_ESTIMATE);
  fftw_execute(fft);
  // free memory
  fftw_destroy_plan(fft);

  // output on screen
//  std::cout << "u1_hat " << std::endl;
//  for (int i=0; i< u1_hat->size(); i++)
//    std::cout << "real  " << real((*u1_hat)[i]) << "    imag  " << imag((*u1_hat)[i]) << std::endl;
  for (int i=0; i< u1_back->size(); i++)
    std::cout << (*u1_back)[i] << std::endl;

  return;
}
#endif

#if 1
/*--------------------------------------------------------------*
 | calculate initial field using fft            rasthofer 04/13 |
 *--------------------------------------------------------------*/
void HomIsoTurbInitialField::CalculateInitialField()
{
#if HAVE_FFTW

  // set and initialize working arrays
  Teuchos::RCP<Teuchos::Array <std::complex<double> > > u1_hat = Teuchos::rcp( new Teuchos::Array<std::complex<double> >(nummodes_*nummodes_*(nummodes_/2+1)));
  Teuchos::RCP<Teuchos::Array <std::complex<double> > > u2_hat = Teuchos::rcp( new Teuchos::Array<std::complex<double> >(nummodes_*nummodes_*(nummodes_/2+1)));
  Teuchos::RCP<Teuchos::Array <std::complex<double> > > u3_hat = Teuchos::rcp( new Teuchos::Array<std::complex<double> >(nummodes_*nummodes_*(nummodes_/2+1)));

  Teuchos::RCP<Teuchos::Array <double> > u1 = Teuchos::rcp( new Teuchos::Array<double>(nummodes_*nummodes_*nummodes_));
  Teuchos::RCP<Teuchos::Array <double> > u2 = Teuchos::rcp( new Teuchos::Array<double>(nummodes_*nummodes_*nummodes_));
  Teuchos::RCP<Teuchos::Array <double> > u3 = Teuchos::rcp( new Teuchos::Array<double>(nummodes_*nummodes_*nummodes_));

  //-------------------------------------------------
  // construction of initial field in spectral space
  //-------------------------------------------------

  // as given in the respective literature, the Fourier coefficients are
  // evaluated in the following intervals
  // k_1: [-nummodes_/2,(nummodes_/2-1)]
  // k_2: [-nummodes_/2,(nummodes_/2-1)]
  // k_3: [-nummodes_/2,0]
  for (int k_1 = (-nummodes_/2); k_1 <= (nummodes_/2-1); k_1++)
  {
    for (int k_2 = (-nummodes_/2); k_2 <= (nummodes_/2-1); k_2++)
    {
      for (int k_3 = (-nummodes_/2); k_3 <= 0; k_3++)
      {
        //get position in u1_hat
        const int pos = (k_3+nummodes_/2) + (nummodes_/2+1) * ((k_2+nummodes_/2) + nummodes_ * (k_1+nummodes_/2));

        if (k_1 == (-nummodes_/2) or k_2 == (-nummodes_/2) or k_3 == (-nummodes_/2))
        {
          // odd-ball wave numbers are set to zero to ensure that solution is a real function
          ((*u1_hat)[pos]).real(0.0);
          // this is important to have here
          ((*u1_hat)[pos]).imag(0.0);
          // remaining analogously
          ((*u2_hat)[pos]).real(0.0);
          ((*u2_hat)[pos]).imag(0.0);
          ((*u3_hat)[pos]).real(0.0);
          ((*u3_hat)[pos]).imag(0.0);
        }
        else if (k_1 == 0 and k_2 == 0 and k_3 == 0)
        {
          // likewise set to zero since there will not be any conjugate complex
          ((*u1_hat)[pos]).real(0.0);
          // this is important to have here
          ((*u1_hat)[pos]).imag(0.0);
          // remaining analogously
          ((*u2_hat)[pos]).real(0.0);
          ((*u2_hat)[pos]).imag(0.0);
          ((*u3_hat)[pos]).real(0.0);
          ((*u3_hat)[pos]).imag(0.0);
        }
        else
        {
          bool calculate = true;
          // check if conjugate complex has already been set
          int pos_conj = -999;
          if (k_3 == 0)
          {
            pos_conj = (-k_3+nummodes_/2) + (nummodes_/2+1) * ((-k_2+nummodes_/2) + nummodes_ * (-k_1+nummodes_/2));

            if (pos_conj < pos)
              calculate = false;
          }

          if (calculate)
          {
            const double k_12 = sqrt(k_1*k_1 + k_2*k_2);
            const double k = sqrt(k_1*k_1 + k_2*k_2 + k_3*k_3);

            double random_theta1 = 0.0;
            double random_theta2 = 0.0;
            double random_phi = 0.0;

            // random numbers are created by one processor and
            // then send to the other processors
            // this ensures that all processors construct the same
            // initial field, which is important to get a matching
            // velocity field in physical space
            if (discret_->Comm().MyPID() == 0)
            {
              DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();
              // set range [0;1] (default: [-1;1])
//              random->SetRandRange(0.0,1.0);
//              random_theta1 = random->Uni();
//              random_theta2 = random->Uni();
//              random_phi = random->Uni();
              random_theta1 = 0.5*random->Uni()+0.5;
              random_theta2 = 0.5*random->Uni()+0.5;
              random_phi = 0.5*random->Uni()+0.5;
            }
            discret_->Comm().Broadcast(&random_theta1,1,0);
            discret_->Comm().Broadcast(&random_theta2,1,0);
            discret_->Comm().Broadcast(&random_phi,1,0);

            // estimate energy at wave number from energy spectrum
            double energy = 0.0;
            if (type_ == INPAR::FLUID::initialfield_hit_comte_bellot_corrsin)
              energy = InterpolateEnergyFromSpectrum(k);
            else
              energy = CalculateEnergyFromSpectrum(k);

            if (energy < 0.0)
            {
              std::cout << "k " << k << std::endl;
              std::cout << "k1  " << k_1 << std::endl;
              std::cout << "k2  " << k_2 << std::endl;
              std::cout << "k3  " << k_3 << std::endl;
              dserror("Negative energy!");
            }

            // remark on the literature:
            // Collis 2002: sqrt(energy/(2*PI*k))
            const double fac = sqrt(energy/(2*PI*k*k));
            // Rogallo 1981: sqrt(energy/(4*PI*k*k))
            // the missing factor 1/2 of Collis version compared to Rogallo version
            // is related to the definition of E from u_i
            // here, we have E = 1/2 * u_i *u_i (see statistics manager)

            // real part, imaginary part
            std::complex<double> alpha(fac * cos(2*PI*random_theta1) * cos(2*PI*random_phi),fac * sin(2*PI*random_theta1) * cos(2*PI*random_phi));
            std::complex<double> beta(fac * cos(2*PI*random_theta2) * sin(2*PI*random_phi),fac * sin(2*PI*random_theta2) * sin(2*PI*random_phi));

            // construct velocity from alpha, beta and k
            if (k_12 > 1.0e-9)
            {
              (*u1_hat)[pos] = (alpha * k * ((double)k_2) + beta * ((double)k_1) * ((double)k_3)) / (k * k_12);
              (*u2_hat)[pos] = (beta * ((double)k_2) * ((double)k_3) - alpha * k * ((double)k_1)) / (k * k_12);
              if (k_3 == 0)
                (*u3_hat)[pos] = -beta * k_12 / k;
              else
                (*u3_hat)[pos] = -(((double)k_1) * ((*u1_hat)[pos]) + ((double)k_2) * ((*u2_hat)[pos])) / ((double)k_3);
            }
            else
            {
              (*u1_hat)[pos] = alpha;
              (*u2_hat)[pos] = beta;
              (*u3_hat)[pos] = -(((double)k_1) * ((*u1_hat)[pos]) + ((double)k_2) * ((*u2_hat)[pos])) / ((double)k_3);
            }
          }
          else
          {
            (*u1_hat)[pos] = conj((*u1_hat)[pos_conj]);
            (*u2_hat)[pos] = conj((*u2_hat)[pos_conj]);
            (*u3_hat)[pos] = conj((*u3_hat)[pos_conj]);
          }
        }
      }
    }
  }

  // transfer to FFTW structure
  // FFTW assumes wave numbers in the following intervals
  // k_1: [0,(nummodes_-1)]
  // k_2: [0,(nummodes_-1)]
  // k_3: [0,nummodes_/2]
  // using peridocity and conjugate symmetry allows for setting
  // the Fourier coefficients in the required interval
  Teuchos::RCP<Teuchos::Array <std::complex<double> > > u1_hat_fftw = Teuchos::rcp( new Teuchos::Array<std::complex<double> >(nummodes_*nummodes_*(nummodes_/2+1)));
  Teuchos::RCP<Teuchos::Array <std::complex<double> > > u2_hat_fftw = Teuchos::rcp( new Teuchos::Array<std::complex<double> >(nummodes_*nummodes_*(nummodes_/2+1)));
  Teuchos::RCP<Teuchos::Array <std::complex<double> > > u3_hat_fftw = Teuchos::rcp( new Teuchos::Array<std::complex<double> >(nummodes_*nummodes_*(nummodes_/2+1)));

  for (int fftw_k_1 = 0; fftw_k_1 <= (nummodes_-1); fftw_k_1++)
  {
    for (int fftw_k_2 = 0; fftw_k_2 <= (nummodes_-1); fftw_k_2++)
    {
      for (int fftw_k_3 = 0; fftw_k_3 <= (nummodes_/2); fftw_k_3++)
      {
        int k_1 = -999;
        int k_2 = -999;
        int k_3 = -999;
        bool conjugate = false;

        if ((fftw_k_1 >= 0 and fftw_k_1 <= (nummodes_/2-1)) and
            (fftw_k_2 >= 0 and fftw_k_2 <= (nummodes_/2-1)) and
            fftw_k_3 == 0)
        {
          // wave number vector is part of construction domain
          // and this value is taken
          k_1 = fftw_k_1;
          k_2 = fftw_k_2;
          k_3 = fftw_k_3;
        }
        else
        {
          // see whether the negative wave vector is in the construction domain
          k_1 = -fftw_k_1;
          k_2 = -fftw_k_2;
          k_3 = -fftw_k_3;
          if ((k_1 >= (-nummodes_/2) and k_1 <= (nummodes_/2-1)) and
              (k_2 >= (-nummodes_/2) and k_2 <= (nummodes_/2-1)) and
              (k_3 >= (-nummodes_/2) and k_3 <= 0))
          {
            // negative wave vector is in the construction domain
            // conjugate complex have to be used here
            conjugate = true;
          }
          else
          {
            // if negative wave vector is not in the construction domain
            // we have to shift it into the domain using the periodicity of the
            // wave number field
            // -k_3 always lies within the construction domain!
            if (k_1 < (-nummodes_/2))
               k_1 += nummodes_;
            if (k_2 < (-nummodes_/2))
              k_2 += nummodes_;

            conjugate = true;
          }
        }

        // get position in u1_hat_fftw
        const int pos = fftw_k_3 + (nummodes_/2+1) * (fftw_k_2 + nummodes_ * fftw_k_1);
        // and in u1_hat
        const int pos_cond = (k_3+nummodes_/2) + (nummodes_/2+1) * ((k_2+nummodes_/2) + nummodes_ * (k_1+nummodes_/2));

        // set value
        if (not conjugate)
        {
          (*u1_hat_fftw)[pos] = (*u1_hat)[pos_cond];
          (*u2_hat_fftw)[pos] = (*u2_hat)[pos_cond];
          (*u3_hat_fftw)[pos] = (*u3_hat)[pos_cond];
        }
        else
        {
          (*u1_hat_fftw)[pos] = conj((*u1_hat)[pos_cond]);
          (*u2_hat_fftw)[pos] = conj((*u2_hat)[pos_cond]);
          (*u3_hat_fftw)[pos] = conj((*u3_hat)[pos_cond]);
        }
      }
    }
  }

  //----------------------------------------
  // fast Fourier transformation using FFTW
  //----------------------------------------

#ifdef HAVE_FFTW
  // set-up
  fftw_plan fft = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
                                       (reinterpret_cast<fftw_complex*>(&((*u1_hat_fftw)[0]))),
                                       &((*u1)[0]), FFTW_ESTIMATE);
  // fft
  fftw_execute(fft);
  // free memory
  fftw_destroy_plan(fft);

  // similar for the remaining two directions
  fftw_plan fft_2 = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
                                         (reinterpret_cast<fftw_complex*>(&((*u2_hat_fftw)[0]))),
                                         &((*u2)[0]), FFTW_ESTIMATE);
  fftw_execute(fft_2);
  // free memory
  fftw_destroy_plan(fft_2);
  fftw_plan fft_3 = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
                                         (reinterpret_cast<fftw_complex*>(&((*u3_hat_fftw)[0]))),
                                         &((*u3)[0]), FFTW_ESTIMATE);
  fftw_execute(fft_3);
  // free memory
  fftw_destroy_plan(fft_3);
  fftw_cleanup();
#else
    dserror("FFTW required for HIT!");
#endif

  //----------------------------------------
  // set velocity field
  //----------------------------------------

  for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
  {
    // get node
    DRT::Node* node = discret_->lRowNode(inode);

    // get coordinates
    LINALG::Matrix<3,1> xyz(true);
    for (int idim = 0; idim < 3; idim++)
        xyz(idim,0) = node->X()[idim];

//    std::cout << "coords " << xyz << std::endl;

    // get global ids of all dofs of the node
    std::vector<int> dofs = discret_->Dof(node);

    // determine position
    std::vector<int> loc(3);
    for (int idim = 0; idim < 3; idim++)
    {
      for (std::size_t rr = 0; rr < coordinates_->size(); rr++)
      {
        if ((xyz(idim,0) <= ((*coordinates_)[rr]+2e-9)) and (xyz(idim,0) >= ((*coordinates_)[rr]-2e-9)))
        {
          // due to periodic boundary conditions,
          // the value at the last node is equal to the one at the first node
          // using this strategy, no special care is required for slave nodes
          if ((int) rr < nummodes_)
           loc[idim] = rr;
          else
           loc[idim] = 0;

          break;
        }
      }
    }

    // get position in transferred velocity vectors u_1, u_2 and u_3
    const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

    // get local dof id corresponding to the global id
    int lid = discret_->DofRowMap()->LID(dofs[0]);
    // set value
    int err = velnp_->ReplaceMyValues(1,&((*u1)[pos]),&lid);
    // analogous for remaining directions
    lid = discret_->DofRowMap()->LID(dofs[1]);
    err = velnp_->ReplaceMyValues(1,&((*u2)[pos]),&lid);
    lid = discret_->DofRowMap()->LID(dofs[2]);
    err = velnp_->ReplaceMyValues(1,&((*u3)[pos]),&lid);
    if (err >0)
      dserror("Could not set initial field!");
  }

  // initialize veln_ as well
  veln_->Update(1.0,*velnp_ ,0.0);

  return;
#else
  dserror("FFTW required");
#endif
}
#endif


/*--------------------------------------------------------------*
 | set and non-dimensionalize experimental data rasthofer 04/13 |
 *--------------------------------------------------------------*/
void HomIsoTurbInitialField::PrepareExparimentalData()
{
  //----------------------------------------
  // set-up wave numbers
  //----------------------------------------

  // wave number given from experiment [cm-1]
  // have to be transferred to [m-1] and non-dimensionalized
  k_exp_.resize(19);
  k_exp_[0] = 0.20;
  k_exp_[1] = 0.25;
  k_exp_[2] = 0.30;
  k_exp_[3] = 0.40;
  k_exp_[4] = 0.50;
  k_exp_[5] = 0.70;
  k_exp_[6] = 1.00;
  k_exp_[7] = 1.50;
  k_exp_[8] = 2.00;
  k_exp_[9] = 2.50;
  k_exp_[10] = 3.00;
  k_exp_[11] = 4.00;
  k_exp_[12] = 6.00;
  k_exp_[13] = 8.00;
  k_exp_[14] = 10.00;
  k_exp_[15] = 12.50;
  k_exp_[16] = 15.00;
  k_exp_[17] = 17.50;
  k_exp_[18] = 20.00;

  // non-dimensionalize wave number (according to Collis 2002)
  // grid size of experiment
  const double M = 0.0508;
  // domain length
  const double L = 10.0 * M;
  // reference length
  const double L_ref = L / (2.0 * PI);

  for (std::size_t rr = 0; rr < k_exp_.size(); rr ++)
      k_exp_[rr] *= (L_ref/0.01);

  //----------------------------------------
  // set-up energy
  //----------------------------------------

  // energy spectrum given from experiment [cm3/s2]
  // have to be transferred to [m3/s2] and non-dimensionalized
  E_exp_.resize(19);
  E_exp_[0] = 129.0;
  E_exp_[1] = 230.0;
  E_exp_[2] = 322.0;
  E_exp_[3] = 435.0;
  E_exp_[4] = 457.0;
  E_exp_[5] = 380.0;
  E_exp_[6] = 270.0;
  E_exp_[7] = 168.0;
  E_exp_[8] = 120.0;
  E_exp_[9] = 89.0;
  E_exp_[10] = 70.3;
  E_exp_[11] = 47.0;
  E_exp_[12] = 24.7;
  E_exp_[13] = 12.6;
  E_exp_[14] = 7.42;
  E_exp_[15] = 3.96;
  E_exp_[16] = 2.33;
  E_exp_[17] = 1.34;
  E_exp_[18] = 0.80;

  // non-dimensionalize energy spectrum
  // inlet velocity of experiment
  const double U_0 = 10.0;
  // reference time
  const double t_ref = 64.0 * M / U_0;

  for (std::size_t rr = 0; rr < E_exp_.size(); rr ++)
      E_exp_[rr] *= ((0.01*0.01*0.01)*(t_ref*t_ref)/(L_ref*L_ref*L_ref));

  return;
}


/*--------------------------------------------------------------*
 | get energy for given wave number             rasthofer 04/13 |
 *--------------------------------------------------------------*/
double HomIsoTurbInitialField::InterpolateEnergyFromSpectrum(double k)
{
  double energy = 0.0;

  // the smallest value for k=1, the next smaller 1.41
  // the following required extrapolation for k<1.6 yields
  // negative values for k=1, which are not physical
  // therefore, energy is set to zero for this case

  // determine position of k
  int position = -1;
  for (std::size_t rr = 0; rr < k_exp_.size(); rr ++)
  {
    if (k < k_exp_[rr])
    {
      position = rr;
      break;
    }
  }

  if (position == -1)
    dserror("Could not determine wave number!");

  if (position > 0)
    // interpolate energy
    energy = E_exp_[position] + (E_exp_[position-1] - E_exp_[position]) / (k_exp_[position] - k_exp_[position-1]) * (k_exp_[position] - k);
  else
    // extrapolate energy
    energy = E_exp_[position+1] - (E_exp_[position+1] - E_exp_[position]) / (k_exp_[position+1] - k_exp_[position]) * (k_exp_[position+1] - k);

  // see above
  if (k == 1)
      energy = 0.0;

  return energy;
}


/*--------------------------------------------------------------*
 | get energy for given wave number             rasthofer 05/13 |
 *--------------------------------------------------------------*/
double HomIsoTurbInitialField::CalculateEnergyFromSpectrum(double k)
{
  // remark: k>0 here!

  double energy = 0.0;

  if (type_ == INPAR::FLUID::initialfield_forced_hit_simple_algebraic_spectrum)
  {
    // initial spectrum as used in Hickel et al. 2006
    energy = 0.5 * pow(k,-5.0/3.0);
  }
  else if (type_ == INPAR::FLUID::initialfield_passive_hit_const_input)
  {
    if (k <=2)
      energy = 0.1 * 1.0;
    else
      energy = 0.1 * pow(2.0,5.0/3.0) * pow(k,-5.0/3.0);
  }
  else if (type_ == INPAR::FLUID::initialfield_forced_hit_numeric_spectrum)
  {
    // initial spectrum as used in Bazilevs et al. 2007 (from Langford & Moser 1999)
    std::vector<double> k_vec(48);
    k_vec[0] = 1.43551;
    k_vec[1] = 1.59008;
    k_vec[2] = 1.74499;
    k_vec[3] = 1.93294;
    k_vec[4] = 2.10165;
    k_vec[5] = 2.3064;
    k_vec[6] = 2.50765;
    k_vec[7] = 2.70109;
    k_vec[8] = 2.90939;
    k_vec[9] = 3.1926;
    k_vec[10] = 3.47102;
    k_vec[11] = 3.77363;
    k_vec[12] = 4.14107;
    k_vec[13] = 4.58665;
    k_vec[14] = 5.22421;
    k_vec[15] = 5.78648;
    k_vec[16] = 6.40926;
    k_vec[17] = 7.0994;
    k_vec[18] = 7.93698;
    k_vec[19] = 8.87338;
    k_vec[20] = 9.92049;
    k_vec[21] = 11.0909;
    k_vec[22] = 12.6323;
    k_vec[23] = 14.2543;
    k_vec[24] = 16.085;
    k_vec[25] = 17.8161;
    k_vec[26] = 19.7332;
    k_vec[27] = 21.6541;
    k_vec[28] = 24.2077;
    k_vec[29] = 26.8124;
    k_vec[30] = 29.4218;
    k_vec[31] = 32.8899;
    k_vec[32] = 36.4263;
    k_vec[33] = 40.343;
    k_vec[34] = 44.6798;
    k_vec[35] = 48.5728;
    k_vec[36] = 52.805;
    k_vec[37] = 55.826;
    k_vec[38] = 60.131;
    k_vec[39] = 64.1623;
    k_vec[40] = 69.1053;
    k_vec[41] = 73.7417;
    k_vec[42] = 78.6891;
    k_vec[43] = 82.4224;
    k_vec[44] = 87.1378;
    k_vec[45] = 92.1229;
    k_vec[46] = 98.3035;
    k_vec[47] = 103.935;

    std::vector<double> e_vec(48);
    e_vec[0] = 14.2756;
    e_vec[1] = 12.8539;
    e_vec[2] = 11.5745;
    e_vec[3] = 10.7899;
    e_vec[4] = 10.0599;
    e_vec[5] = 9.05865;
    e_vec[6] = 8.15762;
    e_vec[7] = 6.85399;
    e_vec[8] = 5.56221;
    e_vec[9] = 4.51327;
    e_vec[10] = 3.79176;
    e_vec[11] = 3.07691;
    e_vec[12] = 2.58484;
    e_vec[13] = 2.09723;
    e_vec[14] = 1.82356;
    e_vec[15] = 1.53182;
    e_vec[16] = 1.28676;
    e_vec[17] = 1.15861;
    e_vec[18] = 0.973181;
    e_vec[19] = 0.817432;
    e_vec[20] = 0.710861;
    e_vec[21] = 0.597094;
    e_vec[22] = 0.501464;
    e_vec[23] = 0.406811;
    e_vec[24] = 0.34168;
    e_vec[25] = 0.287018;
    e_vec[26] = 0.232874;
    e_vec[27] = 0.188958;
    e_vec[28] = 0.148072;
    e_vec[29] = 0.12014;
    e_vec[30] = 0.0941574;
    e_vec[31] = 0.0688353;
    e_vec[32] = 0.0503267;
    e_vec[33] = 0.0367947;
    e_vec[34] = 0.0259835;
    e_vec[35] = 0.0196708;
    e_vec[36] = 0.0148917;
    e_vec[37] = 0.0112761;
    e_vec[38] = 0.00915091;
    e_vec[39] = 0.00646395;
    e_vec[40] = 0.0047269;
    e_vec[41] = 0.00357899;
    e_vec[42] = 0.00270985;
    e_vec[43] = 0.00212455;
    e_vec[44] = 0.00160872;
    e_vec[45] = 0.00121814;
    e_vec[46] = 0.000922317;
    e_vec[47] = 0.000775034;

    // determine position of k
    int position = -1;
    for (std::size_t rr = 0; rr < k_vec.size(); rr ++)
    {
      if (k < k_vec[rr])
      {
        position = rr;
        break;
      }
    }

    if (position == -1)
      dserror("Could not determine wave number!");

    if (position > 0)
      // interpolate energy
      energy = e_vec[position] + (e_vec[position-1] - e_vec[position]) / (k_vec[position] - k_vec[position-1]) * (k_vec[position] - k);
    else
      // extrapolate energy
      energy = e_vec[position+1] - (e_vec[position+1] - e_vec[position]) / (k_vec[position+1] - k_vec[position]) * (k_vec[position+1] - k);
  }

  return energy;
}


};
