/*----------------------------------------------------------------------*/
/*!

\brief routines to calculate forcing for homogeneous isotropic turbulence simulations

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#include <complex>

#ifdef HAVE_FFTW
#include "fftw3.h"
#endif

#include "turbulence_hit_forcing.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_timint_genalpha.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid/fluid_xwall.H"
#include "../drt_fluid/fluid_timint_hdg.H"

#define USE_TRAGET_SPECTRUM
//#define TIME_UPDATE_FORCING_SPECTRUM

namespace FLD
{
  /*--------------------------------------------------------------*
   | constructor                                  rasthofer 04/13 |
   *--------------------------------------------------------------*/
  HomIsoTurbForcing::HomIsoTurbForcing(FluidImplicitTimeInt& timeint)
      : ForcingInterface(),
        forcing_type_(DRT::INPUT::IntegralValue<INPAR::FLUID::ForcingType>(
            timeint.params_->sublist("TURBULENCE MODEL"), "FORCING_TYPE")),
        discret_(timeint.discret_),
        forcing_(timeint.forcing_),
        velnp_(timeint.velnp_),
        velaf_(timeint.velaf_),
        threshold_wavenumber_(
            timeint.params_->sublist("TURBULENCE MODEL").get<double>("THRESHOLD_WAVENUMBER", 0)),
        is_genalpha_(false),
        num_force_steps_(
            timeint.params_->sublist("TURBULENCE MODEL").get<int>("FORCING_TIME_STEPS", 0)),
        dt_(timeint.dta_),
        Pin_(timeint.params_->sublist("TURBULENCE MODEL").get<double>("POWER_INPUT", 0)),
        E_kf_(0.0)
  {
    // set gen-alpha
    TimIntGenAlpha* timeint_genalpha = dynamic_cast<TimIntGenAlpha*>(&timeint);

    if (timeint_genalpha != NULL) is_genalpha_ = true;

    if (timeint.special_flow_ == "forced_homogeneous_isotropic_turbulence" or
        timeint.special_flow_ == "scatra_forced_homogeneous_isotropic_turbulence")
      flow_type_ = forced_homogeneous_isotropic_turbulence;
    else
      flow_type_ = decaying_homogeneous_isotropic_turbulence;

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

    //-------------------------------------------------
    // create set of wave numbers
    //-------------------------------------------------

    // push wave numbers in vector
    {
      wavenumbers_ = Teuchos::rcp(new std::vector<double>);

      wavenumbers_->resize((std::size_t)nummodes_);
      for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++) (*wavenumbers_)[rr] = rr;
    }

    // set size of energy-spectrum vector
    energyspectrum_n_ = Teuchos::rcp(new std::vector<double>);
    energyspectrum_n_->resize(wavenumbers_->size());
    energyspectrum_np_ = Teuchos::rcp(new std::vector<double>);
    energyspectrum_np_->resize(wavenumbers_->size());
    // and initialize with zeros, just to be sure
    for (std::size_t rr = 0; rr < energyspectrum_n_->size(); rr++)
    {
      (*energyspectrum_n_)[rr] = 0.0;
      (*energyspectrum_np_)[rr] = 0.0;
    }

    // forcing factor: factor to multiply Fourier coefficients of velocity
    force_fac_ =
        Teuchos::rcp(new Epetra_SerialDenseVector(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

    return;
  }


  /*--------------------------------------------------------------*
   | initialize energy spectrum by initial field  rasthofer 05/13 |
   *--------------------------------------------------------------*/
  void HomIsoTurbForcing::SetInitialSpectrum(INPAR::FLUID::InitialField init_field_type)
  {
    if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum)
    {
#ifdef USE_TRAGET_SPECTRUM
      if (init_field_type == INPAR::FLUID::initfield_forced_hit_simple_algebraic_spectrum)
      {
        for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
        {
          if ((*wavenumbers_)[rr] > 0.0)
            (*energyspectrum_n_)[rr] = 0.5 * pow((*wavenumbers_)[rr], -5.0 / 3.0);
          else
            (*energyspectrum_n_)[rr] = 0.0;
        }
      }
      else if (init_field_type == INPAR::FLUID::initfield_passive_hit_const_input)
      {
        (*energyspectrum_n_)[0] = 0.0;
        for (std::size_t rr = 1; rr < wavenumbers_->size(); rr++)
        {
          if ((*wavenumbers_)[rr] > 0.0 and (*wavenumbers_)[rr] <= 2.0)
            (*energyspectrum_n_)[rr] = 0.1 * 1.0;
          else
            (*energyspectrum_n_)[rr] =
                0.1 * pow(2.0, 5.0 / 3.0) * pow((*wavenumbers_)[rr], -5.0 / 3.0);
        }
      }
      else if (init_field_type == INPAR::FLUID::initfield_hit_comte_bellot_corrsin)
      {
        //----------------------------------------
        // set-up wave numbers
        //----------------------------------------

        // wave number given from experiment [cm-1]
        // have to be transferred to [m-1] and non-dimensionalized
        std::vector<double> k_exp(19);
        k_exp[0] = 0.20;
        k_exp[1] = 0.25;
        k_exp[2] = 0.30;
        k_exp[3] = 0.40;
        k_exp[4] = 0.50;
        k_exp[5] = 0.70;
        k_exp[6] = 1.00;
        k_exp[7] = 1.50;
        k_exp[8] = 2.00;
        k_exp[9] = 2.50;
        k_exp[10] = 3.00;
        k_exp[11] = 4.00;
        k_exp[12] = 6.00;
        k_exp[13] = 8.00;
        k_exp[14] = 10.00;
        k_exp[15] = 12.50;
        k_exp[16] = 15.00;
        k_exp[17] = 17.50;
        k_exp[18] = 20.00;

        // non-dimensionalize wave number (according to Collis 2002)
        // grid size of experiment
        const double M = 0.0508;
        // domain length
        const double L = 10.0 * M;
        // reference length
        const double L_ref = L / (2.0 * PI);

        for (std::size_t rr = 0; rr < k_exp.size(); rr++) k_exp[rr] *= (L_ref / 0.01);

        //----------------------------------------
        // set-up energy
        //----------------------------------------

        // energy spectrum given from experiment [cm3/s2]
        // have to be transferred to [m3/s2] and non-dimensionalized
        std::vector<double> E_exp(19);
        E_exp[0] = 129.0;
        E_exp[1] = 230.0;
        E_exp[2] = 322.0;
        E_exp[3] = 435.0;
        E_exp[4] = 457.0;
        E_exp[5] = 380.0;
        E_exp[6] = 270.0;
        E_exp[7] = 168.0;
        E_exp[8] = 120.0;
        E_exp[9] = 89.0;
        E_exp[10] = 70.3;
        E_exp[11] = 47.0;
        E_exp[12] = 24.7;
        E_exp[13] = 12.6;
        E_exp[14] = 7.42;
        E_exp[15] = 3.96;
        E_exp[16] = 2.33;
        E_exp[17] = 1.34;
        E_exp[18] = 0.80;

        // non-dimensionalize energy spectrum
        // inlet velocity of experiment
        const double U_0 = 10.0;
        // reference time
        const double t_ref = 64.0 * M / U_0;

        for (std::size_t rr = 0; rr < E_exp.size(); rr++)
          E_exp[rr] *= ((0.01 * 0.01 * 0.01) * (t_ref * t_ref) / (L_ref * L_ref * L_ref));

        // set energy spectrum at k=0 to 0
        (*energyspectrum_n_)[0] = 0.0;
        for (std::size_t rr = 1; rr < wavenumbers_->size(); rr++)
        {
          // the smallest value for k=1, the next smaller 1.41
          // the following required extrapolation for k<1.6 yields
          // negative values for k=1, which are not physical
          // therefore, energy is set to zero for this case -> only initial field
          // here extrapolation from log-plot is used -> 0.006

          // determine position of k
          int position = -1;
          for (std::size_t i = 0; i < k_exp.size(); i++)
          {
            if ((*wavenumbers_)[rr] < k_exp[i])
            {
              position = i;
              break;
            }
          }

          if (position == -1) dserror("Could not determine wave number!");

          if (position > 0)
            // interpolate energy
            (*energyspectrum_n_)[rr] =
                E_exp[position] + (E_exp[position - 1] - E_exp[position]) /
                                      (k_exp[position] - k_exp[position - 1]) *
                                      (k_exp[position] - (*wavenumbers_)[rr]);
          else
            // extrapolate energy
            (*energyspectrum_n_)[rr] =
                E_exp[position + 1] - (E_exp[position + 1] - E_exp[position]) /
                                          (k_exp[position + 1] - k_exp[position]) *
                                          (k_exp[position + 1] - (*wavenumbers_)[rr]);

          // see above
          if ((*wavenumbers_)[rr] == 1) (*energyspectrum_n_)[rr] = 0.006;
        }
      }
      else
        dserror("Other initial spectra than simple algebraic spectrum not yet implemented!");


#else
      CalculateForcing(0);
      TimeUpdateForcing();
#endif
    }

    return;
  }


  /*--------------------------------------------------------------*
   | activate calculation of forcing              rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void HomIsoTurbForcing::ActivateForcing(const bool activate)
  {
    activate_ = activate;
    return;
  }


  /*--------------------------------------------------------------*
   | calculate volume force                       rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void HomIsoTurbForcing::CalculateForcing(const int step)
  {
#if HAVE_FFTW
    // check if forcing is selected
    if (flow_type_ == forced_homogeneous_isotropic_turbulence or
        (flow_type_ == decaying_homogeneous_isotropic_turbulence and step <= num_force_steps_))
    {
      //-------------------------------------------------------------------------------
      // calculate Fourier transformation of velocity
      //-------------------------------------------------------------------------------

      // set and initialize working arrays
      Teuchos::RCP<Teuchos::Array<std::complex<double>>> u1_hat = Teuchos::rcp(
          new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
      Teuchos::RCP<Teuchos::Array<std::complex<double>>> u2_hat = Teuchos::rcp(
          new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
      Teuchos::RCP<Teuchos::Array<std::complex<double>>> u3_hat = Teuchos::rcp(
          new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

      Teuchos::RCP<Teuchos::Array<double>> local_u1 =
          Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
      Teuchos::RCP<Teuchos::Array<double>> local_u2 =
          Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
      Teuchos::RCP<Teuchos::Array<double>> local_u3 =
          Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

      Teuchos::RCP<Teuchos::Array<double>> global_u1 =
          Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
      Teuchos::RCP<Teuchos::Array<double>> global_u2 =
          Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
      Teuchos::RCP<Teuchos::Array<double>> global_u3 =
          Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

      //-----------------------------------
      // prepare Fourier transformation
      //-----------------------------------

      // set solution in local vectors for velocity

      for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
      {
        // get node
        DRT::Node* node = discret_->lRowNode(inode);

        // get coordinates
        LINALG::Matrix<3, 1> xyz(true);
        for (int idim = 0; idim < 3; idim++) xyz(idim, 0) = node->X()[idim];

        // get global ids of all dofs of the node
        std::vector<int> dofs = discret_->Dof(node);

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

        // get position in velocity vectors local_u_1, local_u_2 and local_u_3
        const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

        // get local dof id corresponding to the global id
        int lid = discret_->DofRowMap()->LID(dofs[0]);
        // set value
        if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum or
            (forcing_type_ == INPAR::FLUID::fixed_power_input and (not is_genalpha_)))
        {
          (*local_u1)[pos] = (*velnp_)[lid];
          // analogously for remaining directions
          lid = discret_->DofRowMap()->LID(dofs[1]);
          (*local_u2)[pos] = (*velnp_)[lid];
          lid = discret_->DofRowMap()->LID(dofs[2]);
          (*local_u3)[pos] = (*velnp_)[lid];
        }
        else
        {
          (*local_u1)[pos] = (*velaf_)[lid];
          // analogously for remaining directions
          lid = discret_->DofRowMap()->LID(dofs[1]);
          (*local_u2)[pos] = (*velaf_)[lid];
          lid = discret_->DofRowMap()->LID(dofs[2]);
          (*local_u3)[pos] = (*velaf_)[lid];
        }
      }

      // get values form all processors
      // number of nodes without slave nodes
      const int countallnodes = nummodes_ * nummodes_ * nummodes_;
      discret_->Comm().SumAll(&((*local_u1)[0]), &((*global_u1)[0]), countallnodes);

      discret_->Comm().SumAll(&((*local_u2)[0]), &((*global_u2)[0]), countallnodes);

      discret_->Comm().SumAll(&((*local_u3)[0]), &((*global_u3)[0]), countallnodes);

      //----------------------------------------
      // fast Fourier transformation using FFTW
      //----------------------------------------

      // note: this is not very efficient, since each
      // processor does the fft and there is no communication

#ifdef HAVE_FFTW
      // set-up
      fftw_plan fft = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u1)[0]),
          (reinterpret_cast<fftw_complex*>(&((*u1_hat)[0]))), FFTW_ESTIMATE);
      // fft
      fftw_execute(fft);
      // free memory
      fftw_destroy_plan(fft);
      // analogously for remaining directions
      fftw_plan fft_2 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u2)[0]),
          (reinterpret_cast<fftw_complex*>(&((*u2_hat)[0]))), FFTW_ESTIMATE);
      fftw_execute(fft_2);
      // free memory
      fftw_destroy_plan(fft_2);
      fftw_plan fft_3 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u3)[0]),
          (reinterpret_cast<fftw_complex*>(&((*u3_hat)[0]))), FFTW_ESTIMATE);
      fftw_execute(fft_3);
      // free memory
      fftw_destroy_plan(fft_3);
      fftw_cleanup();
#else
      dserror("FFTW required for HIT!");
#endif

      // scale solution (not done in the fftw routine)
      for (int i = 0; i < u1_hat->size(); i++)
      {
        (*u1_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
        (*u2_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
        (*u3_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
      }

      //----------------------------------------
      // compute energy spectrum
      //----------------------------------------

      // transfer from FFTW structure to intervals around zero
      // FFTW assumes wave numbers in the following intervals
      // k_1: [0,(nummodes_-1)]
      // k_2: [0,(nummodes_-1)]
      // k_3: [0,nummodes_/2]
      // here, we would like to have
      // k_1: [-nummodes_/2,(nummodes_/2-1)]
      // k_2: [-nummodes_/2,(nummodes_/2-1)]
      // k_3: [-nummodes_/2,0]
      // using peridocity and conjugate symmetry allows for setting
      // the Fourier coefficients in the required interval

      // reset energy spectrum at time n+1/n+af
      for (std::size_t rr = 0; rr < energyspectrum_np_->size(); rr++)
        (*energyspectrum_np_)[rr] = 0.0;

      // reset (just to be sure)
      E_kf_ = 0.0;

      // the complete number of modes is required here
      // hence, we have k_3: [-nummodes_/2,(nummodes_/2-1)]
#if 0  // outdated, to be removed
    for (int fftw_k_1 = 0; fftw_k_1 <= (nummodes_-1); fftw_k_1++)
    {
      for (int fftw_k_2 = 0; fftw_k_2 <= (nummodes_-1); fftw_k_2++)
      {
        for (int fftw_k_3 = 0; fftw_k_3 <= (nummodes_-1); fftw_k_3++)
        {

          int pos_fftw_k_1 = -999;
          int pos_fftw_k_2 = -999;
          int pos_fftw_k_3 = -999;

          int k_1 = -999;
          int k_2 = -999;
          int k_3 = -999;

          if ((fftw_k_1 >= 0 and fftw_k_1 <= (nummodes_/2-1)) and
              (fftw_k_2 >= 0 and fftw_k_2 <= (nummodes_/2-1)) and
              (fftw_k_3 >= 0 and fftw_k_3 <= 0 )) //and fftw_k_3 <= (nummodes_/2-1)))
          {
            // wave number vector is part of construction domain
            // and this value is taken
            k_1 = fftw_k_1;
            k_2 = fftw_k_2;
            k_3 = fftw_k_3;

            pos_fftw_k_1 = fftw_k_1;
            pos_fftw_k_2 = fftw_k_2;
            pos_fftw_k_3 = fftw_k_3;
          }
          else
          {
            // check if current wave vector is in the domain of fftw
            // if not shift to domain using periodicity and conjugate symmetry
            if (fftw_k_3 > (nummodes_/2))
            {
              // find corresponding value in fftw-range first
              int intermediate_fftw_3 = fftw_k_3 - nummodes_;
              // since intermediate_fftw_3 is still outside of domain,
              // get position of conjugate complex
              intermediate_fftw_3 *= -1;
              int intermediate_fftw_1 = -fftw_k_1;
              int intermediate_fftw_2 = -fftw_k_2;
              // check if this position is inside intervals
              if ((intermediate_fftw_1 >= 0 and intermediate_fftw_1 <= (nummodes_-1))
                   and (intermediate_fftw_2 >= 0 and intermediate_fftw_2 <= (nummodes_-1))
                   and (intermediate_fftw_3 >= 0 and intermediate_fftw_3 <= (nummodes_/2)))
              {
                pos_fftw_k_1 = intermediate_fftw_1;
                pos_fftw_k_2 = intermediate_fftw_2;
                pos_fftw_k_3 = intermediate_fftw_3;
              }
              else
              {
                // shift intermediate_fftw_1 and intermediate_fftw_2 into fftw interval
                // intermediate_fftw_3 always lies within the fftw interval
                if (intermediate_fftw_1 < 0)
                  intermediate_fftw_1 += nummodes_;
                if (intermediate_fftw_2 < 0)
                  intermediate_fftw_2 += nummodes_;

                pos_fftw_k_1 = intermediate_fftw_1;
                pos_fftw_k_2 = intermediate_fftw_2;
                pos_fftw_k_3 = intermediate_fftw_3;
              }
            }
            else
            {
              pos_fftw_k_1 = fftw_k_1;
              pos_fftw_k_2 = fftw_k_2;
              pos_fftw_k_3 = fftw_k_3;
            }

            // see whether the negative wave vector is in the construction domain
            k_1 = -pos_fftw_k_1;
            k_2 = -pos_fftw_k_2;
            k_3 = -pos_fftw_k_3;
            if ((k_1 >= (-nummodes_/2) and k_1 <= (nummodes_/2-1)) and
                (k_2 >= (-nummodes_/2) and k_2 <= (nummodes_/2-1)) and
                (k_3 >= (-nummodes_/2) and k_3 <= 0 )) //and k_3 <= (nummodes_/2-1)))
            {

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
            }
          }
#else
      for (int k_1 = (-nummodes_ / 2); k_1 <= (nummodes_ / 2 - 1); k_1++)
      {
        for (int k_2 = (-nummodes_ / 2); k_2 <= (nummodes_ / 2 - 1); k_2++)
        {
          for (int k_3 = (-nummodes_ / 2); k_3 <= (nummodes_ / 2 - 1); k_3++)
          {
            // initialize position in FFTW vectors
            int pos_fftw_k_1 = -999;
            int pos_fftw_k_2 = -999;
            int pos_fftw_k_3 = -999;

            // check if current wave vector lies within the fftw domain
            if ((k_1 >= 0 and k_1 <= (nummodes_ / 2 - 1)) and
                (k_2 >= 0 and k_2 <= (nummodes_ / 2 - 1)) and
                (k_3 >= 0 and k_3 <= (nummodes_ / 2 - 1)))
            {
              pos_fftw_k_1 = k_1;
              pos_fftw_k_2 = k_2;
              pos_fftw_k_3 = k_3;
            }
            else
            {
              // if k_3 is < 0, we have to take the conjugate
              // to get into the FFTW domain
              if (k_3 < 0)
              {
                int k_conj_1 = -k_1;
                int k_conj_2 = -k_2;
                int k_conj_3 = -k_3;

                // check if conjugate wave vector lies within the fftw domain
                // this has to be fulfilled for k_3 but not for k_1 and k_2
                if ((k_conj_1 >= 0 and k_conj_1 <= (nummodes_ - 1)) and
                    (k_conj_2 >= 0 and k_conj_2 <= (nummodes_ - 1)) and
                    (k_conj_3 >= 0 and k_conj_3 <= (nummodes_ - 1)))
                {
                  pos_fftw_k_1 = k_conj_1;
                  pos_fftw_k_2 = k_conj_2;
                  pos_fftw_k_3 = k_conj_3;
                }
                else
                {
                  if (not(k_conj_3 >= 0 and k_conj_3 <= (nummodes_ - 1)))
                    dserror("k_3 in fftw domain expected!");

                  // shift k_1 and k_2 into fftw domain
                  if (k_conj_1 < 0) k_conj_1 += nummodes_;
                  if (k_conj_2 < 0) k_conj_2 += nummodes_;

                  if ((k_conj_1 >= 0 and k_conj_1 <= (nummodes_ - 1)) and
                      (k_conj_2 >= 0 and k_conj_2 <= (nummodes_ - 1)) and
                      (k_conj_3 >= 0 and k_conj_3 <= (nummodes_ - 1)))
                  {
                    pos_fftw_k_1 = k_conj_1;
                    pos_fftw_k_2 = k_conj_2;
                    pos_fftw_k_3 = k_conj_3;
                  }
                  else
                    dserror("Position in fftw domain expected!");
                }
              }
              else
              {
                int k_shift_1 = k_1;
                int k_shift_2 = k_2;
                int k_shift_3 = k_3;

                if (not(k_shift_3 >= 0 and k_shift_3 <= (nummodes_ - 1)))
                  dserror("k_3 in fftw domain expected!");

                // shift k_1 and k_2 into fftw domain
                if (k_shift_1 < 0) k_shift_1 += nummodes_;
                if (k_shift_2 < 0) k_shift_2 += nummodes_;

                if ((k_shift_1 >= 0 and k_shift_1 <= (nummodes_ - 1)) and
                    (k_shift_2 >= 0 and k_shift_2 <= (nummodes_ - 1)) and
                    (k_shift_3 >= 0 and k_shift_3 <= (nummodes_ - 1)))
                {
                  pos_fftw_k_1 = k_shift_1;
                  pos_fftw_k_2 = k_shift_2;
                  pos_fftw_k_3 = k_shift_3;
                }
                else
                  dserror("Position in fftw domain expected!");
              }
            }
#endif

      // get position in u1_hat
      const int pos =
          pos_fftw_k_3 + (nummodes_ / 2 + 1) * (pos_fftw_k_2 + nummodes_ * pos_fftw_k_1);

      // calculate energy
      // E = 1/2 * u_i * conj(u_i)
      // u_i * conj(u_i) = real(u_i)^2 + imag(u_i)^2
      // const std::complex<double> energy = 0.5 * ((*u1_hat)[pos] * conj((*u1_hat)[pos])
      //                                          + (*u2_hat)[pos] * conj((*u2_hat)[pos])
      //                                          + (*u3_hat)[pos] * conj((*u3_hat)[pos]));
      // instead
      const double energy =
          0.5 * (norm((*u1_hat)[pos]) + norm((*u2_hat)[pos]) + norm((*u3_hat)[pos]));

      if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum)
      {
        // get wave number
        const double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);

        // insert into sampling vector
        // find position via k
        for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
        {
          if (k > ((*wavenumbers_)[rr] - 0.5) and k <= ((*wavenumbers_)[rr] + 0.5))
          {
            (*energyspectrum_np_)[rr] += energy;
          }
        }
      }
      else if (forcing_type_ == INPAR::FLUID::fixed_power_input)
      {
        if ((abs(k_1) < threshold_wavenumber_ and abs(k_2) < threshold_wavenumber_ and
                abs(k_3) < threshold_wavenumber_) and
            ((k_1 * k_1 + k_2 * k_2 + k_3 * k_3) > 0))
          E_kf_ += energy;
      }
      else
        dserror("Unknown forcing type!");
    }
  }
}

//--------------------------------------------------------
// forcing factor from energy spectrum
//--------------------------------------------------------

// reset to zero (just to be sure)
for (int rr = 0; rr < (nummodes_ * nummodes_ * (nummodes_ / 2 + 1)); rr++) (*force_fac_)(rr) = 0.0;

for (int fftw_k_1 = 0; fftw_k_1 <= (nummodes_ - 1); fftw_k_1++)
{
  for (int fftw_k_2 = 0; fftw_k_2 <= (nummodes_ - 1); fftw_k_2++)
  {
    for (int fftw_k_3 = 0; fftw_k_3 <= (nummodes_ / 2 + 1); fftw_k_3++)
    {
      int k_1 = -999;
      int k_2 = -999;
      int k_3 = -999;

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
        }
        else
        {
          // if negative wave vector is not in the construction domain
          // we have to shift it into the domain using the periodicity of the
          // wave number field
          // -k_3 always lies within the construction domain!
          if (k_1 < (-nummodes_ / 2)) k_1 += nummodes_;
          if (k_2 < (-nummodes_ / 2)) k_2 += nummodes_;
        }
      }

      if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum)
      {
        // compute linear compensation factor from energy spectrum

        // following Hickel 2007
        //         (1/(2*E))* dE/dt if k <= k_t
        // C(k) =
        //         0                else
        // threshold wave number k_t
        // E intermediate energy spectrum, i.e. solution without forcing

        // get wave number
        double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);
        if (k <= threshold_wavenumber_ and k > 0)
        {
          int rr_k = 0;
          for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
          {
            if ((*wavenumbers_)[rr] > k)
            {
              rr_k = rr;
              break;
            }
          }

          // interpolate spectrum to wave number
          double E_n = Interpolate(k, (*wavenumbers_)[rr_k - 1], (*wavenumbers_)[rr_k],
              (*energyspectrum_n_)[rr_k - 1], (*energyspectrum_n_)[rr_k]);
          double E_np = Interpolate(k, (*wavenumbers_)[rr_k - 1], (*wavenumbers_)[rr_k],
              (*energyspectrum_np_)[rr_k - 1], (*energyspectrum_np_)[rr_k]);

          // get position in fac-vector
          const int pos = fftw_k_3 + (nummodes_ / 2 + 1) * (fftw_k_2 + nummodes_ * fftw_k_1);

          // calculate C
          (*force_fac_)(pos) = (1.0 / (2.0 * E_np)) * (E_np - E_n) / dt_;
        }
      }
      else if (forcing_type_ == INPAR::FLUID::fixed_power_input)
      {
        // compute power input

        // following Bazilevs et al. 2007
        //
        // fac = Pin/(2*E_kf)
        // Pin: fixed power input
        // E_kf: energy of wave numbers with |k_i|<k_f

        // note: we have "-" here, since we take -fac*u in the function below
        const double fac = -Pin_ / (2.0 * E_kf_);

        if ((abs(k_1) < threshold_wavenumber_ and abs(k_2) < threshold_wavenumber_ and
                abs(k_3) < threshold_wavenumber_) and
            ((k_1 * k_1 + k_2 * k_2 + k_3 * k_3) > 0))
        {
          // get position in fac-vector
          const int pos = fftw_k_3 + (nummodes_ / 2 + 1) * (fftw_k_2 + nummodes_ * fftw_k_1);

          // set factor
          (*force_fac_)(pos) = fac;
        }
      }
      else
        dserror("Unknown forcing type!");
    }
  }
}

// reset E_kf_
E_kf_ = 0.0;
}

return;
#else
    dserror("FFTW required");
#endif
}


/*--------------------------------------------------------------*
 | get forcing                                   rasthofer 04/13 |
 *--------------------------------------------------------------*/
void HomIsoTurbForcing::UpdateForcing(const int step)
{
#if HAVE_FFTW
  // check if forcing is selected
  if (activate_ and
      (flow_type_ == forced_homogeneous_isotropic_turbulence or
          (flow_type_ == decaying_homogeneous_isotropic_turbulence and step <= num_force_steps_)))
  {
    //-------------------------------------------------------------------------------
    // calculate Fourier transformation of velocity
    //-------------------------------------------------------------------------------

    // set and initialize working arrays
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u1_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u2_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u3_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

    Teuchos::RCP<Teuchos::Array<double>> local_u1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> local_u2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> local_u3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    Teuchos::RCP<Teuchos::Array<double>> global_u1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> global_u2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> global_u3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    //-----------------------------------
    // prepare Fourier transformation
    //-----------------------------------

    // set solution in local vectors for velocity

    for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
    {
      // get node
      DRT::Node* node = discret_->lRowNode(inode);

      // get coordinates
      LINALG::Matrix<3, 1> xyz(true);
      for (int idim = 0; idim < 3; idim++) xyz(idim, 0) = node->X()[idim];

      // get global ids of all dofs of the node
      std::vector<int> dofs = discret_->Dof(node);

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

      // get position in velocity vectors local_u_1, local_u_2 and local_u_3
      const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

      // get local dof id corresponding to the global id
      int lid = discret_->DofRowMap()->LID(dofs[0]);
      // set value
      if (not is_genalpha_)
      {
        (*local_u1)[pos] = (*velnp_)[lid];
        // analogously for remaining directions
        lid = discret_->DofRowMap()->LID(dofs[1]);
        (*local_u2)[pos] = (*velnp_)[lid];
        lid = discret_->DofRowMap()->LID(dofs[2]);
        (*local_u3)[pos] = (*velnp_)[lid];
      }
      else
      {
        (*local_u1)[pos] = (*velaf_)[lid];
        // analogously for remaining directions
        lid = discret_->DofRowMap()->LID(dofs[1]);
        (*local_u2)[pos] = (*velaf_)[lid];
        lid = discret_->DofRowMap()->LID(dofs[2]);
        (*local_u3)[pos] = (*velaf_)[lid];
      }
    }

    // get values form all processors
    // number of nodes without slave nodes
    const int countallnodes = nummodes_ * nummodes_ * nummodes_;
    discret_->Comm().SumAll(&((*local_u1)[0]), &((*global_u1)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u2)[0]), &((*global_u2)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u3)[0]), &((*global_u3)[0]), countallnodes);

    //----------------------------------------
    // fast Fourier transformation using FFTW
    //----------------------------------------

    // note: this is not very efficient, since each
    // processor does the fft and there is no communication

#ifdef HAVE_FFTW
    // set-up
    fftw_plan fft = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u1)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u1_hat)[0]))), FFTW_ESTIMATE);
    // fft
    fftw_execute(fft);
    // free memory
    fftw_destroy_plan(fft);
    // analogously for remaining directions
    fftw_plan fft_2 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u2)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u2_hat)[0]))), FFTW_ESTIMATE);
    fftw_execute(fft_2);
    // free memory
    fftw_destroy_plan(fft_2);
    fftw_plan fft_3 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u3)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u3_hat)[0]))), FFTW_ESTIMATE);
    fftw_execute(fft_3);
    // free memory
    fftw_destroy_plan(fft_3);
    fftw_cleanup();
#else
      dserror("FFTW required for HIT!");
#endif

    // scale solution (not done in the fftw routine)
    for (int i = 0; i < u1_hat->size(); i++)
    {
      (*u1_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
      (*u2_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
      (*u3_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
    }

    //----------------------------------------
    // set forcing vector
    //----------------------------------------

    // Fourier coefficients of forcing
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> f1_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> f2_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> f3_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    // where f_hat = -C(k) * u_hat according to Hickel 2007
    // C denotes a linear compensation factor
    // or where C denotes the dissipation dependent factor from
    // Bazilevs et al.

    Teuchos::RCP<Teuchos::Array<double>> f1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> f2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> f3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    for (int rr = 0; rr < (nummodes_ * nummodes_ * (nummodes_ / 2 + 1)); rr++)
    {
      (*f1_hat)[rr] = -((*force_fac_)(rr)) * ((*u1_hat)[rr]);
      (*f2_hat)[rr] = -((*force_fac_)(rr)) * ((*u2_hat)[rr]);
      (*f3_hat)[rr] = -((*force_fac_)(rr)) * ((*u3_hat)[rr]);
    }

    //----------------------------------------
    // fast Fourier transformation using FFTW
    //----------------------------------------

#ifdef HAVE_FFTW
    // setup
    fftw_plan fft_back = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
        (reinterpret_cast<fftw_complex*>(&((*f1_hat)[0]))), &((*f1)[0]), FFTW_ESTIMATE);
    // fft
    fftw_execute(fft_back);
    // free memory
    fftw_destroy_plan(fft_back);

    // similar for the remaining two directions
    fftw_plan fft_back_2 = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
        (reinterpret_cast<fftw_complex*>(&((*f2_hat)[0]))), &((*f2)[0]), FFTW_ESTIMATE);
    fftw_execute(fft_back_2);
    // free memory
    fftw_destroy_plan(fft_back_2);
    fftw_plan fft_back_3 = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
        (reinterpret_cast<fftw_complex*>(&((*f3_hat)[0]))), &((*f3)[0]), FFTW_ESTIMATE);
    fftw_execute(fft_back_3);
    // free memory
    fftw_destroy_plan(fft_back_3);
    fftw_cleanup();
#else
      dserror("FFTW required for HIT!");
#endif

    //----------------------------------------
    // set force vector
    //----------------------------------------

    for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
    {
      // get node
      DRT::Node* node = discret_->lRowNode(inode);

      // get coordinates
      LINALG::Matrix<3, 1> xyz(true);
      for (int idim = 0; idim < 3; idim++) xyz(idim, 0) = node->X()[idim];

      // get global ids of all dofs of the node
      std::vector<int> dofs = discret_->Dof(node);

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

      // get position in transferred velocity vectors u_1, u_2 and u_3
      const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

      // get local dof id corresponding to the global id
      int lid = discret_->DofRowMap()->LID(dofs[0]);
      // set value
      int err = forcing_->ReplaceMyValues(1, &((*f1)[pos]), &lid);
      // analogous for remaining directions
      lid = discret_->DofRowMap()->LID(dofs[1]);
      err = forcing_->ReplaceMyValues(1, &((*f2)[pos]), &lid);
      lid = discret_->DofRowMap()->LID(dofs[2]);
      err = forcing_->ReplaceMyValues(1, &((*f3)[pos]), &lid);
      if (err > 0) dserror("Could not set forcing!");
    }
  }
  else
    // set force to zero
    forcing_->PutScalar(0.0);

  return;
#else
    dserror("FFTW required");
#endif
}


/*--------------------------------------------------------------*
 | time update of energy spectrum               rasthofer 04/13 |
 *--------------------------------------------------------------*/
void HomIsoTurbForcing::TimeUpdateForcing()
{
#ifdef TIME_UPDATE_FORCING_SPECTRUM
  if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum)
  {
    // compute E^n+1 from forced solution
    CalculateForcing(0);
    // update energy spectrum
    for (std::size_t rr = 0; rr < energyspectrum_n_->size(); rr++)
    {
      (*energyspectrum_n_)[rr] = (*energyspectrum_np_)[rr];
      (*energyspectrum_np_)[rr] = 0.0;
    }
  }
#endif

  // reset to zero
  E_kf_ = 0.0;

  for (int rr = 0; rr < (nummodes_ * nummodes_ * (nummodes_ / 2 + 1)); rr++)
    (*force_fac_)(rr) = 0.0;

  forcing_->PutScalar(0.0);

  return;
}


/*--------------------------------------------------------------*
 | constructor                                         bk 03/15 |
 *--------------------------------------------------------------*/
HomIsoTurbForcingHDG::HomIsoTurbForcingHDG(FluidImplicitTimeInt& timeint)
    : HomIsoTurbForcing(timeint)
{
  // here we are using the interior velocity
  TimIntHDG* hdgfluid = dynamic_cast<TimIntHDG*>(&timeint);
  if (hdgfluid == NULL) dserror("this should be a hdg time integer");

  // we want to use the interior velocity here
  velnp_ = hdgfluid->ReturnIntVelnp();

  //-----------------------------------
  // determine number of modes
  //-----------------------------------

  // number of modes equal to 5 times number of elements in one spatial direction
  nummodes_ *= 5;

  //-------------------------------------------------
  // create set of node coordinates
  //-------------------------------------------------

  // push coordinates in vector
  {
    Teuchos::RCP<std::vector<double>> copycoordinates = Teuchos::rcp(new std::vector<double>);

    for (std::vector<double>::iterator coord1 = coordinates_->begin();
         coord1 != coordinates_->end(); ++coord1)
    {
      copycoordinates->push_back(*coord1);
    }

    coordinates_->clear();

    double elesize = abs(copycoordinates->at(1) - copycoordinates->at(0));
    // use 5 sampling locations in each element in each direction
    const double localcoords[5] = {0.9, 0.7, 0.5, 0.3, 0.1};
    for (std::vector<double>::iterator coord1 = copycoordinates->begin();
         coord1 != copycoordinates->end(); ++coord1)
    {
      if (coord1 != copycoordinates->begin())
        for (int i = 0; i < 5; i++) coordinates_->push_back(*coord1 - elesize * localcoords[i]);
    }
  }


  //-------------------------------------------------
  // create set of wave numbers
  //-------------------------------------------------

  // push wave numbers in vector
  {
    wavenumbers_ = Teuchos::rcp(new std::vector<double>);

    wavenumbers_->resize((std::size_t)nummodes_);
    for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++) (*wavenumbers_)[rr] = rr;
  }

  // set size of energy-spectrum vector
  energyspectrum_n_->resize(wavenumbers_->size());
  energyspectrum_np_->resize(wavenumbers_->size());
  // and initialize with zeros, just to be sure
  for (std::size_t rr = 0; rr < energyspectrum_n_->size(); rr++)
  {
    (*energyspectrum_n_)[rr] = 0.0;
    (*energyspectrum_np_)[rr] = 0.0;
  }

  // forcing factor: factor to multiply Fourier coefficients of velocity
  force_fac_ =
      Teuchos::rcp(new Epetra_SerialDenseVector(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

  return;
}

/*--------------------------------------------------------------*
 | initialize energy spectrum by initial field         bk 04/15 |
 *--------------------------------------------------------------*/
void HomIsoTurbForcingHDG::SetInitialSpectrum(INPAR::FLUID::InitialField init_field_type)
{
#ifdef USE_TRAGET_SPECTRUM
  HomIsoTurbForcing::SetInitialSpectrum(init_field_type);
#else
    dserror("only USE_TARGET_SPECTRUM implemented for HDG");
#endif
  return;
}

/*--------------------------------------------------------------*
 | calculate volume force                              bk 03/15 |
 *--------------------------------------------------------------*/
void HomIsoTurbForcingHDG::CalculateForcing(const int step)
{
#if HAVE_FFTW
  // check if forcing is selected
  if (flow_type_ == forced_homogeneous_isotropic_turbulence or
      (flow_type_ == decaying_homogeneous_isotropic_turbulence and step <= num_force_steps_))
  {
    //-------------------------------------------------------------------------------------------------
    // calculate Fourier transformation of velocity
    //-------------------------------------------------------------------------------------------------

    // set and initialize working arrays
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u1_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u2_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u3_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

    Teuchos::RCP<Teuchos::Array<double>> local_u1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> local_u2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> local_u3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    Teuchos::RCP<Teuchos::Array<double>> global_u1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> global_u2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> global_u3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    //-----------------------------------
    // prepare Fourier transformation
    //-----------------------------------

    // set solution in local vectors for velocity

    // procedure: first, go down to the element and get velocity at 5x5x5 positions
    // also get coordinates of these positions
    // the insert the values according to their coordinates in local_u1, etc.
    // then everything can be transformed



    //  //new for HDG: go down to element and get 5x5x5 values
    //  call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action", FLD::interpolate_hdg_for_hit);

    if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum or
        (forcing_type_ == INPAR::FLUID::fixed_power_input and (not is_genalpha_)))
    {
      discret_->SetState(1, "intvelnp", velnp_);
    }
    else
      dserror("it seems like you need velaf_ here, which is not implemented for hit and hdg yet");

    std::vector<int> dummy;
    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;

    for (int el = 0; el < discret_->NumMyRowElements(); ++el)
    {
      Epetra_SerialDenseVector interpolVec;
      DRT::Element* ele = discret_->lRowElement(el);

      interpolVec.Resize(5 * 5 * 5 * 6);  // 5*5*5 points: velx, vely, velz, x, y, z

      ele->Evaluate(params, *discret_, dummy, dummyMat, dummyMat, interpolVec, dummyVec, dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i = 0; i < 5 * 5 * 5; ++i)
      {
        // get coordinates
        LINALG::Matrix<3, 1> xyz(true);
        for (int d = 0; d < 3; ++d) xyz(d) = interpolVec(i * 6 + d + 3);
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
                dserror("I think that this should not happen");

              break;
            }
          }
        }

        // get position in velocity vectors local_u_1, local_u_2 and local_u_3
        const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

        // set value
        (*local_u1)[pos] = interpolVec(i * 6 + 0);
        (*local_u2)[pos] = interpolVec(i * 6 + 1);
        (*local_u3)[pos] = interpolVec(i * 6 + 2);
      }
    }

    // get values form all processors
    // number of nodes without slave nodes
    const int countallnodes = nummodes_ * nummodes_ * nummodes_;
    discret_->Comm().SumAll(&((*local_u1)[0]), &((*global_u1)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u2)[0]), &((*global_u2)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u3)[0]), &((*global_u3)[0]), countallnodes);

    //----------------------------------------
    // fast Fourier transformation using FFTW
    //----------------------------------------

    // note: this is not very efficient, since each
    // processor does the fft and there is no communication

#ifdef HAVE_FFTW
    // set-up
    fftw_plan fft = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u1)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u1_hat)[0]))), FFTW_ESTIMATE);
    // fft
    fftw_execute(fft);
    // free memory
    fftw_destroy_plan(fft);

    // analogously for remaining directions
    fftw_plan fft_2 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u2)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u2_hat)[0]))), FFTW_ESTIMATE);
    fftw_execute(fft_2);
    // free memory
    fftw_destroy_plan(fft_2);
    fftw_plan fft_3 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u3)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u3_hat)[0]))), FFTW_ESTIMATE);
    fftw_execute(fft_3);
    // free memory
    fftw_destroy_plan(fft_3);
    fftw_cleanup();
#else
      dserror("FFTW required for HIT!");
#endif

    // scale solution (not done in the fftw routine)
    for (int i = 0; i < u1_hat->size(); i++)
    {
      (*u1_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
      (*u2_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
      (*u3_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
    }

    //----------------------------------------
    // compute energy spectrum
    //----------------------------------------

    // transfer from FFTW structure to intervals around zero
    // FFTW assumes wave numbers in the following intervals
    // k_1: [0,(nummodes_-1)]
    // k_2: [0,(nummodes_-1)]
    // k_3: [0,nummodes_/2]
    // here, we would like to have
    // k_1: [-nummodes_/2,(nummodes_/2-1)]
    // k_2: [-nummodes_/2,(nummodes_/2-1)]
    // k_3: [-nummodes_/2,0]
    // using peridocity and conjugate symmetry allows for setting
    // the Fourier coefficients in the required interval

    // the complete number of modes is required here
    // hence, we have k_3: [-nummodes_/2,(nummodes_/2-1)]

    // reset energy spectrum at time n+1/n+af
    for (std::size_t rr = 0; rr < energyspectrum_np_->size(); rr++) (*energyspectrum_np_)[rr] = 0.0;

    // reset (just to be sure)
    E_kf_ = 0.0;

    for (int k_1 = (-nummodes_ / 2); k_1 <= (nummodes_ / 2 - 1); k_1++)
    {
      for (int k_2 = (-nummodes_ / 2); k_2 <= (nummodes_ / 2 - 1); k_2++)
      {
        for (int k_3 = (-nummodes_ / 2); k_3 <= (nummodes_ / 2 - 1); k_3++)
        {
          // initialize position in FFTW vectors
          int pos_fftw_k_1 = -999;
          int pos_fftw_k_2 = -999;
          int pos_fftw_k_3 = -999;

          // check if current wave vector lies within the fftw domain
          if ((k_1 >= 0 and k_1 <= (nummodes_ / 2 - 1)) and
              (k_2 >= 0 and k_2 <= (nummodes_ / 2 - 1)) and
              (k_3 >= 0 and k_3 <= (nummodes_ / 2 - 1)))
          {
            pos_fftw_k_1 = k_1;
            pos_fftw_k_2 = k_2;
            pos_fftw_k_3 = k_3;
          }
          else
          {
            // if k_3 is < 0, we have to take the conjugate
            // to get into the FFTW domain
            if (k_3 < 0)
            {
              int k_conj_1 = -k_1;
              int k_conj_2 = -k_2;
              int k_conj_3 = -k_3;

              // check if conjugate wave vector lies within the fftw domain
              // this has to be fulfilled for k_3 but not for k_1 and k_2
              if ((k_conj_1 >= 0 and k_conj_1 <= (nummodes_ - 1)) and
                  (k_conj_2 >= 0 and k_conj_2 <= (nummodes_ - 1)) and
                  (k_conj_3 >= 0 and k_conj_3 <= (nummodes_ - 1)))
              {
                pos_fftw_k_1 = k_conj_1;
                pos_fftw_k_2 = k_conj_2;
                pos_fftw_k_3 = k_conj_3;
              }
              else
              {
                if (not(k_conj_3 >= 0 and k_conj_3 <= (nummodes_ - 1)))
                  dserror("k_3 in fftw domain expected!");

                // shift k_1 and k_2 into fftw domain
                if (k_conj_1 < 0) k_conj_1 += nummodes_;
                if (k_conj_2 < 0) k_conj_2 += nummodes_;

                if ((k_conj_1 >= 0 and k_conj_1 <= (nummodes_ - 1)) and
                    (k_conj_2 >= 0 and k_conj_2 <= (nummodes_ - 1)) and
                    (k_conj_3 >= 0 and k_conj_3 <= (nummodes_ - 1)))
                {
                  pos_fftw_k_1 = k_conj_1;
                  pos_fftw_k_2 = k_conj_2;
                  pos_fftw_k_3 = k_conj_3;
                }
                else
                  dserror("Position in fftw domain expected!");
              }
            }
            else
            {
              int k_shift_1 = k_1;
              int k_shift_2 = k_2;
              int k_shift_3 = k_3;

              if (not(k_shift_3 >= 0 and k_shift_3 <= (nummodes_ - 1)))
                dserror("k_3 in fftw domain expected!");

              // shift k_1 and k_2 into fftw domain
              if (k_shift_1 < 0) k_shift_1 += nummodes_;
              if (k_shift_2 < 0) k_shift_2 += nummodes_;

              if ((k_shift_1 >= 0 and k_shift_1 <= (nummodes_ - 1)) and
                  (k_shift_2 >= 0 and k_shift_2 <= (nummodes_ - 1)) and
                  (k_shift_3 >= 0 and k_shift_3 <= (nummodes_ - 1)))
              {
                pos_fftw_k_1 = k_shift_1;
                pos_fftw_k_2 = k_shift_2;
                pos_fftw_k_3 = k_shift_3;
              }
              else
                dserror("Position in fftw domain expected!");
            }
          }


          // get position in u1_hat
          const int pos =
              pos_fftw_k_3 + (nummodes_ / 2 + 1) * (pos_fftw_k_2 + nummodes_ * pos_fftw_k_1);

          // calculate energy
          // E = 1/2 * u_i * conj(u_i)
          // u_i * conj(u_i) = real(u_i)^2 + imag(u_i)^2
          // const std::complex<double> energy = 0.5 * ((*u1_hat)[pos] * conj((*u1_hat)[pos])
          //                                          + (*u2_hat)[pos] * conj((*u2_hat)[pos])
          //                                          + (*u3_hat)[pos] * conj((*u3_hat)[pos]));
          // instead
          const double energy =
              0.5 * (norm((*u1_hat)[pos]) + norm((*u2_hat)[pos]) + norm((*u3_hat)[pos]));

          if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum)
          {
            // get wave number
            const double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);

            // insert into sampling vector
            // find position via k
            for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
            {
              if (k > ((*wavenumbers_)[rr] - 0.5) and k <= ((*wavenumbers_)[rr] + 0.5))
              {
                (*energyspectrum_np_)[rr] += energy;
              }
            }
          }
          else if (forcing_type_ == INPAR::FLUID::fixed_power_input)
          {
            if ((abs(k_1) < threshold_wavenumber_ and abs(k_2) < threshold_wavenumber_ and
                    abs(k_3) < threshold_wavenumber_) and
                ((k_1 * k_1 + k_2 * k_2 + k_3 * k_3) > 0))
              E_kf_ += energy;
          }
          else
            dserror("Unknown forcing type!");
        }
      }
    }

    //--------------------------------------------------------
    // forcing factor from energy spectrum
    //--------------------------------------------------------

    // reset to zero (just to be sure)
    for (int rr = 0; rr < (nummodes_ * nummodes_ * (nummodes_ / 2 + 1)); rr++)
      (*force_fac_)(rr) = 0.0;

    for (int fftw_k_1 = 0; fftw_k_1 <= (nummodes_ - 1); fftw_k_1++)
    {
      for (int fftw_k_2 = 0; fftw_k_2 <= (nummodes_ - 1); fftw_k_2++)
      {
        for (int fftw_k_3 = 0; fftw_k_3 <= (nummodes_ / 2 + 1); fftw_k_3++)
        {
          int k_1 = -999;
          int k_2 = -999;
          int k_3 = -999;

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
            }
            else
            {
              // if negative wave vector is not in the construction domain
              // we have to shift it into the domain using the periodicity of the
              // wave number field
              // -k_3 always lies within the construction domain!
              if (k_1 < (-nummodes_ / 2)) k_1 += nummodes_;
              if (k_2 < (-nummodes_ / 2)) k_2 += nummodes_;
            }
          }

          if (forcing_type_ == INPAR::FLUID::linear_compensation_from_intermediate_spectrum)
          {
            // compute linear compensation factor from energy spectrum

            // following Hickel 2007
            //         (1/(2*E))* dE/dt if k <= k_t
            // C(k) =
            //         0                else
            // threshold wave number k_t
            // E intermediate energy spectrum, i.e. solution without forcing

            // get wave number
            double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);
            if (k <= threshold_wavenumber_ and k > 0)
            {
              int rr_k = 0;
              for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
              {
                if ((*wavenumbers_)[rr] > k)
                {
                  rr_k = rr;
                  break;
                }
              }

              // interpolate spectrum to wave number
              double E_n = Interpolate(k, (*wavenumbers_)[rr_k - 1], (*wavenumbers_)[rr_k],
                  (*energyspectrum_n_)[rr_k - 1], (*energyspectrum_n_)[rr_k]);
              double E_np = Interpolate(k, (*wavenumbers_)[rr_k - 1], (*wavenumbers_)[rr_k],
                  (*energyspectrum_np_)[rr_k - 1], (*energyspectrum_np_)[rr_k]);

              // get position in fac-vector
              const int pos = fftw_k_3 + (nummodes_ / 2 + 1) * (fftw_k_2 + nummodes_ * fftw_k_1);

              // calculate C
              (*force_fac_)(pos) = (1.0 / (2.0 * E_np)) * (E_np - E_n) / dt_;
            }
          }
          else if (forcing_type_ == INPAR::FLUID::fixed_power_input)
          {
            // compute power input

            // following Bazilevs et al. 2007
            //
            // fac = Pin/(2*E_kf)
            // Pin: fixed power input
            // E_kf: energy of wave numbers with |k_i|<k_f

            // note: we have "-" here, since we take -fac*u in the function below
            const double fac = -Pin_ / (2.0 * E_kf_);

            if ((abs(k_1) < threshold_wavenumber_ and abs(k_2) < threshold_wavenumber_ and
                    abs(k_3) < threshold_wavenumber_) and
                ((k_1 * k_1 + k_2 * k_2 + k_3 * k_3) > 0))
            {
              // get position in fac-vector
              const int pos = fftw_k_3 + (nummodes_ / 2 + 1) * (fftw_k_2 + nummodes_ * fftw_k_1);

              // set factor
              (*force_fac_)(pos) = fac;
            }
          }
          else
            dserror("Unknown forcing type!");
        }
      }
    }

    // reset E_kf_
    E_kf_ = 0.0;
  }
  discret_->ClearState(true);
  return;
#else
    dserror("FFTW required");
#endif
}

/*--------------------------------------------------------------*
 | get forcing                                         bk 03/15 |
 *--------------------------------------------------------------*/
void HomIsoTurbForcingHDG::UpdateForcing(const int step)
{
#if HAVE_FFTW
  // check if forcing is selected
  if (activate_ and
      (flow_type_ == forced_homogeneous_isotropic_turbulence or
          (flow_type_ == decaying_homogeneous_isotropic_turbulence and step <= num_force_steps_)))
  {
    //-------------------------------------------------------------------------------------------------
    // calculate Fourier transformation of velocity
    //-------------------------------------------------------------------------------------------------

    // set and initialize working arrays
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u1_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u2_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> u3_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));

    Teuchos::RCP<Teuchos::Array<double>> local_u1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> local_u2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> local_u3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    Teuchos::RCP<Teuchos::Array<double>> global_u1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> global_u2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> global_u3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    //-----------------------------------
    // prepare Fourier transformation
    //-----------------------------------

    // set solution in local vectors for velocity

    // procedure: first, go down to the element and get velocity at 5x5x5 positions
    // also get coordinates of these positions
    // the insert the values according to their coordinates in local_u1, etc.
    // then everything can be transformed



    //  //new for HDG: go down to element and get 5x5x5 values
    //  call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action", FLD::interpolate_hdg_for_hit);

    discret_->SetState(1, "intvelnp", velnp_);

    std::vector<int> dummy;
    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;


    for (int el = 0; el < discret_->NumMyRowElements(); ++el)
    {
      Epetra_SerialDenseVector interpolVec;
      DRT::Element* ele = discret_->lRowElement(el);

      interpolVec.Resize(5 * 5 * 5 * 6);  // 5*5*5 points: velx, vely, velz, x, y, z

      ele->Evaluate(params, *discret_, dummy, dummyMat, dummyMat, interpolVec, dummyVec, dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i = 0; i < 5 * 5 * 5; ++i)
      {
        // get coordinates
        LINALG::Matrix<3, 1> xyz(true);
        for (int d = 0; d < 3; ++d) xyz(d) = interpolVec(i * 6 + d + 3);
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
                dserror("I think that this should not happen");

              break;
            }
          }
        }

        // get position in velocity vectors local_u_1, local_u_2 and local_u_3
        const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

        // set value
        (*local_u1)[pos] = interpolVec(i * 6 + 0);
        (*local_u2)[pos] = interpolVec(i * 6 + 1);
        (*local_u3)[pos] = interpolVec(i * 6 + 2);
      }
    }

    // get values form all processors
    // number of nodes without slave nodes
    const int countallnodes = nummodes_ * nummodes_ * nummodes_;
    discret_->Comm().SumAll(&((*local_u1)[0]), &((*global_u1)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u2)[0]), &((*global_u2)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u3)[0]), &((*global_u3)[0]), countallnodes);

    //----------------------------------------
    // fast Fourier transformation using FFTW
    //----------------------------------------

    // note: this is not very efficient, since each
    // processor does the fft and there is no communication

#ifdef HAVE_FFTW
    // set-up
    fftw_plan fft = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u1)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u1_hat)[0]))), FFTW_ESTIMATE);
    // fft
    fftw_execute(fft);
    // free memory
    fftw_destroy_plan(fft);

    // analogously for remaining directions
    fftw_plan fft_2 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u2)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u2_hat)[0]))), FFTW_ESTIMATE);
    fftw_execute(fft_2);
    // free memory
    fftw_destroy_plan(fft_2);
    fftw_plan fft_3 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_u3)[0]),
        (reinterpret_cast<fftw_complex*>(&((*u3_hat)[0]))), FFTW_ESTIMATE);
    fftw_execute(fft_3);
    // free memory
    fftw_destroy_plan(fft_3);
    fftw_cleanup();
#else
      dserror("FFTW required for HIT!");
#endif

    // scale solution (not done in the fftw routine)
    for (int i = 0; i < u1_hat->size(); i++)
    {
      (*u1_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
      (*u2_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
      (*u3_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
    }

    //----------------------------------------
    // set forcing vector
    //----------------------------------------

    // Fourier coefficients of forcing
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> f1_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> f2_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<std::complex<double>>> f3_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    // where f_hat = -C(k) * u_hat according to Hickel 2007
    // C denotes a linear compensation factor
    // or where C denotes the dissipation dependent factor from
    // Bazilevs et al.

    Teuchos::RCP<Teuchos::Array<double>> f1 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> f2 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> f3 =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));

    for (int rr = 0; rr < (nummodes_ * nummodes_ * (nummodes_ / 2 + 1)); rr++)
    {
      (*f1_hat)[rr] = -((*force_fac_)(rr)) * ((*u1_hat)[rr]);
      (*f2_hat)[rr] = -((*force_fac_)(rr)) * ((*u2_hat)[rr]);
      (*f3_hat)[rr] = -((*force_fac_)(rr)) * ((*u3_hat)[rr]);
    }

    //----------------------------------------
    // fast Fourier transformation using FFTW
    //----------------------------------------

#ifdef HAVE_FFTW
    // setup
    fftw_plan fft_back = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
        (reinterpret_cast<fftw_complex*>(&((*f1_hat)[0]))), &((*f1)[0]), FFTW_ESTIMATE);
    // fft
    fftw_execute(fft_back);
    // free memory
    fftw_destroy_plan(fft_back);

    // similar for the remaining two directions
    fftw_plan fft_back_2 = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
        (reinterpret_cast<fftw_complex*>(&((*f2_hat)[0]))), &((*f2)[0]), FFTW_ESTIMATE);
    fftw_execute(fft_back_2);
    // free memory
    fftw_destroy_plan(fft_back_2);
    fftw_plan fft_back_3 = fftw_plan_dft_c2r_3d(nummodes_, nummodes_, nummodes_,
        (reinterpret_cast<fftw_complex*>(&((*f3_hat)[0]))), &((*f3)[0]), FFTW_ESTIMATE);
    fftw_execute(fft_back_3);
    // free memory
    fftw_destroy_plan(fft_back_3);
    fftw_cleanup();
#else
      dserror("FFTW required for HIT!");
#endif

    //----------------------------------------
    // set force vector
    //----------------------------------------

    // for 1st evaluate
    // use same params as above
    // this is a dummy, forcing_ should be zero is written in the first components of interpolVec
    discret_->ClearState(true);

    discret_->SetState(1, "intvelnp", velnp_);

    // this is the real value
    discret_->SetState(1, "forcing", forcing_);

    // for 2nd evaluate
    const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
    Epetra_SerialDenseVector elevec1, elevec3;
    Epetra_SerialDenseMatrix elemat1, elemat2;
    Teuchos::ParameterList initParams;
    initParams.set<int>("action", FLD::project_hdg_force_on_dof_vec_for_hit);

    // loop over all elements on the processor
    DRT::Element::LocationArray la(2);
    for (int el = 0; el < discret_->NumMyRowElements(); ++el)
    {
      // 1st evaluate
      DRT::Element* ele = discret_->lRowElement(el);

      std::vector<int> dummy;
      Epetra_SerialDenseMatrix dummyMat;
      Epetra_SerialDenseVector dummyVec;
      Epetra_SerialDenseVector interpolVec;
      interpolVec.Resize(5 * 5 * 5 * 6);  // 5*5*5 points: velx, vely, velz, x, y, z

      ele->Evaluate(params, *discret_, dummy, dummyMat, dummyMat, interpolVec, dummyVec, dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i = 0; i < 5 * 5 * 5; ++i)
      {
        // get coordinates
        LINALG::Matrix<3, 1> xyz(true);
        for (int d = 0; d < 3; ++d) xyz(d) = interpolVec(i * 6 + d + 3);
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
                dserror("I think that this should not happen");

              break;
            }
          }
        }

        // get position in velocity vectors local_u_1, local_u_2 and local_u_3
        const int pos = loc[2] + nummodes_ * (loc[1] + nummodes_ * loc[0]);

        // set value
        interpolVec(i * 6 + 0) = (*f1)[pos];
        interpolVec(i * 6 + 1) = (*f2)[pos];
        interpolVec(i * 6 + 2) = (*f3)[pos];
      }

      // 2nd evaluate
      ele->LocationVector(*discret_, la, false);
      if (elevec1.M() != discret_->NumDof(1, ele)) elevec1.Shape(discret_->NumDof(1, ele), 1);

      ele->Evaluate(
          initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, interpolVec, elevec3);

      if (ele->Owner() == discret_->Comm().MyPID())
      {
        std::vector<int> localDofs = discret_->Dof(1, ele);
        dsassert(localDofs.size() == static_cast<std::size_t>(elevec1.M()), "Internal error");
        for (unsigned int i = 0; i < localDofs.size(); ++i)
          localDofs[i] = intdofrowmap->LID(localDofs[i]);
        forcing_->ReplaceMyValues(localDofs.size(), elevec1.A(), &localDofs[0]);
      }
    }
  }
  else
    // set force to zero
    forcing_->PutScalar(0.0);
  discret_->ClearState(true);

  return;
#else
    dserror("FFTW required");
#endif
  return;
}

/*--------------------------------------------------------------*
 | constructor                                  bk        12/14 |
 *--------------------------------------------------------------*/
PeriodicHillForcing::PeriodicHillForcing(FluidImplicitTimeInt& timeint)
    : ForcingInterface(),
      discret_(timeint.discret_),
      forcing_(timeint.forcing_),
      velnp_(timeint.velnp_),
      velaf_(timeint.velaf_),
      myxwall_(timeint.xwall_),
      oldforce_(0.0),
      oldflow_(49.46),
      idealmassflow_(49.46),
      length_(252.0),
      step_(1),
      count_(0),
      sum_(0.0)
{
  if (discret_->Comm().MyPID() == 0)
    std::cout << "\nforcing for periodic hill such that a mass flow of " << idealmassflow_
              << " is achieved" << std::endl;
  std::vector<DRT::Condition*> bodycond;
  discret_->GetCondition("VolumeNeumann", bodycond);
  const std::vector<double>* val = bodycond[0]->Get<std::vector<double>>("val");
  oldforce_ = val->at(0);
}

/*--------------------------------------------------------------*
 | time update of periodic hill forcing                bk 12/14 |
 *--------------------------------------------------------------*/
void PeriodicHillForcing::UpdateForcing(const int step)
{
  step_ = step;
  return;
}

/*--------------------------------------------------------------*
 | time update of periodic hill forcing                bk 12/14 |
 *--------------------------------------------------------------*/
void PeriodicHillForcing::TimeUpdateForcing()
{
  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", FLD::calc_mass_flow_periodic_hill);

  if (myxwall_ != Teuchos::null) myxwall_->SetXWallParams(eleparams);

  eleparams.set<double>("length", length_);

  discret_->SetState("velnp", velnp_);

  const Epetra_Map* elementrowmap = discret_->ElementRowMap();
  Teuchos::RCP<Epetra_MultiVector> massflvec =
      Teuchos::rcp(new Epetra_MultiVector(*elementrowmap, 1, true));

  // optional: elementwise defined div u may be written to standard output file (not implemented
  // yet)
  discret_->EvaluateScalars(eleparams, massflvec);

  discret_->ClearState();

  Teuchos::RCP<Epetra_MultiVector> massflvecneg =
      Teuchos::rcp(new Epetra_MultiVector(*elementrowmap, 1, true));

  // take into account negative mass flux at the inflow
  for (int i = 0; i < discret_->ElementRowMap()->NumMyElements(); ++i)
  {
    double locflow = (*((*massflvec)(0)))[i];
    if (locflow < -1.0e-9)
    {
      (*((*massflvec)(0)))[i] = 0.0;
      (*((*massflvecneg)(0)))[i] = locflow;
    }
  }

  double massflowpos = 0.0;
  massflvec->Norm1(&massflowpos);
  double massflowneg = 0.0;
  massflvecneg->Norm1(&massflowneg);
  double massflow = massflowpos - massflowneg;

  double dm = massflow - oldflow_;
  double dgoalm = idealmassflow_ - massflow;

  double newforce = 0.0;

  // the inital value does not have any meaning, so assume 0 for first step
  if (step_ < 2) dm = 0.0;

  // new estimated force
  // first contribution makes the system want to get back to the ideal value (spring)
  // the second contribution is a penalty on quick changes (damping)
  // the constants are empirical
  newforce = 500.0 * dgoalm - 30000.0 * dm + oldforce_;

  // now insert values in vector
  forcing_->PutScalar(0.0);

  for (int i = 0; i < discret_->NodeRowMap()->NumMyElements(); ++i)
  {
    int gid = discret_->NodeRowMap()->GID(i);
    DRT::Node* node = discret_->gNode(gid);
    if (!node) dserror("Cannot find node");

    int firstgdofid = discret_->Dof(0, node, 0);

    int firstldofid = forcing_->Map().LID(firstgdofid);

    int err = forcing_->ReplaceMyValue(firstldofid, 0, newforce);
    if (err != 0) dserror("something went wrong during replacemyvalue");
  }

  oldforce_ = newforce;
  oldflow_ = massflow;

  // some statistical data
  count_++;
  sum_ += newforce;

  // provide some information about the current condition
  if (discret_->Comm().MyPID() == 0)
    std::cout << "current mass flux:  " << oldflow_ << "/" << 49.46 << "  force:  " << oldforce_
              << "/" << sum_ / (double)count_ << std::endl;
  return;
}

}  // namespace FLD
