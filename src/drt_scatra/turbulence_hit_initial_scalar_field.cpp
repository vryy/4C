/*----------------------------------------------------------------------*/
/*! \file


\brief routines to initialize homogeneous isotropic turbulence simulations with passive scalar
transport

\level 2

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#include <complex>
#include <math.h>

#ifdef HAVE_FFTW
#include "fftw3.h"
#endif

#include "turbulence_hit_initial_scalar_field.H"

#include "scatra_timint_implicit.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

namespace SCATRA
{
  /*--------------------------------------------------------------*
   | constructor                                  rasthofer 04/13 |
   *--------------------------------------------------------------*/
  HomIsoTurbInitialScalarField::HomIsoTurbInitialScalarField(
      ScaTraTimIntImpl& timeint, const INPAR::SCATRA::InitialField initfield)
      : discret_(timeint.discret_), phinp_(timeint.phinp_), phin_(timeint.phin_), type_(initfield)
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
      case 1728:
      {
        nummodes_ = 12;
        break;
      }
      case 4096:
      {
        nummodes_ = 16;
        break;
      }
      case 13824:
      {
        nummodes_ = 24;
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
    std::set<double, LineSortCriterion> coords;
    // loop all nodes and store x1-coordinate
    for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
    {
      DRT::Node* node = discret_->lRowNode(inode);
      if ((node->X()[1] < 2e-9 && node->X()[1] > -2e-9) and
          (node->X()[2] < 2e-9 && node->X()[2] > -2e-9))
        coords.insert(node->X()[0]);
    }

    // communicate coordinates to all procs via round Robin loop
    {
#ifdef PARALLEL
      int myrank = discret_->Comm().MyPID();
#endif
      int numprocs = discret_->Comm().NumProc();

      std::vector<char> sblock;
      std::vector<char> rblock;

#ifdef PARALLEL
      // create an exporter for point to point communication
      DRT::Exporter exporter(discret_->Comm());
#endif

      // communicate coordinates
      for (int np = 0; np < numprocs; ++np)
      {
        DRT::PackBuffer data;

        for (std::set<double, LineSortCriterion>::iterator x1line = coords.begin();
             x1line != coords.end(); ++x1line)
        {
          DRT::ParObject::AddtoPack(data, *x1line);
        }
        data.StartPacking();
        for (std::set<double, LineSortCriterion>::iterator x1line = coords.begin();
             x1line != coords.end(); ++x1line)
        {
          DRT::ParObject::AddtoPack(data, *x1line);
        }
        std::swap(sblock, data());

#ifdef PARALLEL
        MPI_Request request;
        int tag = myrank;

        int frompid = myrank;
        int topid = (myrank + 1) % numprocs;

        int length = sblock.size();

        exporter.ISend(frompid, topid, &(sblock[0]), sblock.size(), tag, request);

        rblock.clear();

        // receive from predecessor
        frompid = (myrank + numprocs - 1) % numprocs;
        exporter.ReceiveAny(frompid, tag, rblock, length);

        if (tag != (myrank + numprocs - 1) % numprocs)
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
        rblock = sblock;
#endif

        // unpack received block into set of all coordinates
        {
          std::vector<double> coordsvec;

          coordsvec.clear();

          std::vector<char>::size_type index = 0;
          while (index < rblock.size())
          {
            double onecoord;
            DRT::ParObject::ExtractfromPack(index, rblock, onecoord);
            coords.insert(onecoord);
          }
        }
      }
    }
    // push coordinates in vectors
    {
      coordinates_ = Teuchos::rcp(new std::vector<double>);

      for (std::set<double, LineSortCriterion>::iterator coord1 = coords.begin();
           coord1 != coords.end(); ++coord1)
      {
        coordinates_->push_back(*coord1);
      }
    }

    return;
  }


  /*--------------------------------------------------------------*
   | calculate initial field using fft            rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void HomIsoTurbInitialScalarField::CalculateInitialField()
  {
#ifdef HAVE_FFTW
    // set and initialize working arrays
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> phi_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

    Teuchos::RCP<Teuchos::Array<double>> phi =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    //-------------------------------------------------
    // construction of initial field in spectral space
    //-------------------------------------------------

    // as given in the respective literature, the Fourier coefficients are
    // evaluated in the following intervals
    // k_1: [-nummodes_/2,(nummodes_/2-1)]
    // k_2: [-nummodes_/2,(nummodes_/2-1)]
    // k_3: [-nummodes_/2,0]
    for (int k_1 = (-nummodes_ / 2); k_1 <= (nummodes_ / 2 - 1); k_1++)
    {
      for (int k_2 = (-nummodes_ / 2); k_2 <= (nummodes_ / 2 - 1); k_2++)
      {
        for (int k_3 = (-nummodes_ / 2); k_3 <= 0; k_3++)
        {
          // get position in phi_hat
          const int pos =
              (k_3 + nummodes_ / 2) +
              (nummodes_ / 2 + 1) * ((k_2 + nummodes_ / 2) + nummodes_ * (k_1 + nummodes_ / 2));

          if (k_1 == (-nummodes_ / 2) or k_2 == (-nummodes_ / 2) or k_3 == (-nummodes_ / 2))
          {
            // odd-ball wave numbers are set to zero to ensure that solution is real function
            ((*phi_hat)[pos]).real(0.0);
            // this is important to have here
            ((*phi_hat)[pos]).imag(0.0);
          }
          else if (k_1 == 0 and k_2 == 0 and k_3 == 0)
          {
            // likewise set to zero since there will not be any conjugate complex
            ((*phi_hat)[pos]).real(0.0);
            // this is important to have here
            ((*phi_hat)[pos]).imag(0.0);
          }
          else
          {
            bool calculate = true;
            // check if conjugate complex has already been set
            int pos_conj = -999;
            if (k_3 == 0)
            {
              pos_conj = (-k_3 + nummodes_ / 2) +
                         (nummodes_ / 2 + 1) *
                             ((-k_2 + nummodes_ / 2) + nummodes_ * (-k_1 + nummodes_ / 2));

              if (pos_conj < pos) calculate = false;
            }

            if (calculate)
            {
              const double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);

              double random_theta = 0.0;

              // random numbers are created by one processor and
              // then send to the other processors
              // this ensures that all processors construct the same
              // initial field, which is important to get a matching
              // scalar field in physical space
              if (discret_->Comm().MyPID() == 0)
              {
                DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();
                // set range [0;1] (default: [-1;1])
                //              random->SetRandRange(0.0,1.0);
                //              random_theta = random->Uni();
                // use this version to get random field different from fluid
                random_theta = 0.5 * random->Uni() + 0.5;
              }
              discret_->Comm().Broadcast(&random_theta, 1, 0);

              // estimate energy at wave number from energy spectrum
              const double energy = CalculateEnergyFromSpectrum(k);

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
              const double fac = sqrt(energy / (2 * PI * k * k));
              // Rogallo 1981: sqrt(energy/(4*PI*k*k))
              // the missing factor 1/2 of Collis version compared to Rogallo version
              // is related to the definition of E from phi
              // here, we have E = 1/2 * phi * phi (see statistics manager)

              // real part, imaginary part
              std::complex<double> alpha(
                  fac * cos(2 * PI * random_theta), fac * sin(2 * PI * random_theta));
              (*phi_hat)[pos] = alpha;
            }
            else
            {
              (*phi_hat)[pos] = conj((*phi_hat)[pos_conj]);
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
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> phi_hat_fftw = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

    for (int fftw_k_1 = 0; fftw_k_1 <= (nummodes_ - 1); fftw_k_1++)
    {
      for (int fftw_k_2 = 0; fftw_k_2 <= (nummodes_ - 1); fftw_k_2++)
      {
        for (int fftw_k_3 = 0; fftw_k_3 <= (nummodes_ / 2); fftw_k_3++)
        {
          int k_1 = -999;
          int k_2 = -999;
          int k_3 = -999;
          bool conjugate = false;

          if ((fftw_k_1 >= 0 and fftw_k_1 <= (nummodes_ / 2 - 1)) and
              (fftw_k_2 >= 0 and fftw_k_2 <= (nummodes_ / 2 - 1)) and fftw_k_3 == 0)
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
            if ((k_1 >= (-nummodes_ / 2) and k_1 <= (nummodes_ / 2 - 1)) and
                (k_2 >= (-nummodes_ / 2) and k_2 <= (nummodes_ / 2 - 1)) and
                (k_3 >= (-nummodes_ / 2) and k_3 <= 0))
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
              if (k_1 < (-nummodes_ / 2)) k_1 += nummodes_;
              if (k_2 < (-nummodes_ / 2)) k_2 += nummodes_;

              conjugate = true;
            }
          }

          // get position in phi_hat_fftw
          const int pos = fftw_k_3 + (nummodes_ / 2 + 1) * (fftw_k_2 + nummodes_ * fftw_k_1);
          // and in phi_hat
          const int pos_cond =
              (k_3 + nummodes_ / 2) +
              (nummodes_ / 2 + 1) * ((k_2 + nummodes_ / 2) + nummodes_ * (k_1 + nummodes_ / 2));

          // set value
          if (not conjugate)
          {
            (*phi_hat_fftw)[pos] = (*phi_hat)[pos_cond];
          }
          else
          {
            (*phi_hat_fftw)[pos] = conj((*phi_hat)[pos_cond]);
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
        (reinterpret_cast<fftw_complex*>(&((*phi_hat_fftw)[0]))), &((*phi)[0]), FFTW_ESTIMATE);
    // fft
    fftw_execute(fft);
    // free memory
    fftw_destroy_plan(fft);
    fftw_cleanup();
#else
    dserror("FFTW required for HIT!");
#endif

    //----------------------------------------
    // set scalar field
    //----------------------------------------

    for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
    {
      // get node
      DRT::Node* node = discret_->lRowNode(inode);

      // get coordinates
      LINALG::Matrix<3, 1> xyz(true);
      for (int idim = 0; idim < 3; idim++) xyz(idim, 0) = node->X()[idim];

      // get global ids of all dofs of the node
      std::vector<int> dofs = discret_->Dof(0, node);
      if (dofs.size() > 1) dserror("Only one dof per node for homogeneous isotropic turbulence!");

      // determine position
      std::vector<int> loc(3);
      for (int idim = 0; idim < 3; idim++)
      {
        for (std::size_t rr = 0; rr < coordinates_->size(); rr++)
        {
          if ((xyz(idim, 0) <= ((*coordinates_)[rr] + 2e-9)) and
              (xyz(idim, 0) >= ((*coordinates_)[rr] - 2e-9)))
          {
            // due to periodic boundary conditions,
            // the value at the last node is equal to the one at the first node
            // using this strategy, no special care is required for slave nodes
            if ((int)rr < nummodes_)
              loc[idim] = rr;
            else
              loc[idim] = 0;

            break;
          }
        }
      }

      // get position in transferred phi vector
      const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

      // get local dof id corresponding to the global id
      int lid = discret_->DofRowMap()->LID(dofs[0]);
      // set value
      int err = phinp_->ReplaceMyValues(1, &((*phi)[pos]), &lid);

      if (err > 0) dserror("Could not set initial field!");
    }

    // initialize phin_ as well
    phin_->Update(1.0, *phinp_, 0.0);

    return;
#else
    dserror("FFTW required");
#endif
  }


  /*--------------------------------------------------------------*
   | get energy for given wave number             rasthofer 05/13 |
   *--------------------------------------------------------------*/
  double HomIsoTurbInitialScalarField::CalculateEnergyFromSpectrum(double k)
  {
    // remark: k > 0 here
    double energy = 0.0;

    if (type_ == INPAR::SCATRA::initialfield_forced_hit_low_Sc)
    {
      if (k > 2.0)
        energy = 0.1 * pow(2.0, 5.0 / 3.0) * pow(k, -5.0 / 3.0);
      else
        energy = 0.1 * 1.0;
    }
    else if (type_ == INPAR::SCATRA::initialfield_forced_hit_high_Sc)
    {
      if (k > 2.0)
        energy = 0.1 * 2.0 * pow(k, -1.0);
      else
        energy = 0.1 * 1.0;
    }
    else
      dserror("Unkown initial field!");

    return energy;
  }


};  // namespace SCATRA
