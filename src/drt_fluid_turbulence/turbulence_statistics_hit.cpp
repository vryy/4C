/*----------------------------------------------------------------------*/
/*!

\brief routines for homogeneous isotropic turbulence

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#include <complex>

#ifdef HAVE_FFTW
#include "fftw3.h"
#endif

#include "turbulence_statistics_hit.H"

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/standardtypes_cpp.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"

namespace FLD
{
  /*--------------------------------------------------------------*
   | constructor                                  rasthofer 04/13 |
   *--------------------------------------------------------------*/
  TurbulenceStatisticsHit::TurbulenceStatisticsHit(Teuchos::RCP<DRT::Discretization> actdis,
      Teuchos::ParameterList& params, const std::string& statistics_outfilename, const bool forced)
      : discret_(actdis), params_(params), statistics_outfilename_(statistics_outfilename)
  {
    // set type
    if (forced == true)
      type_ = forced_homogeneous_isotropic_turbulence;
    else
      type_ = decaying_homogeneous_isotropic_turbulence;

    //-----------------------------------
    // output to screen
    //-----------------------------------

    if (discret_->Comm().MyPID() == 0)
    {
      std::cout << "This is the turbulence statistics manager of\n";
      std::cout << "homogeneous isotropic turbulence:\n\n";
      if (type_ == forced_homogeneous_isotropic_turbulence)
        std::cout << "   FORCED HOMOGENEOUS ISOTROPIC TURBULENCE\n";
      else
        std::cout << "   DECAYING HOMOGENEOUS ISOTROPIC TURBULENCE\n";
      std::cout << std::endl;
    }

    //-----------------------------------
    // determine number of modes
    //-----------------------------------

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

    // push coordinates in vector
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
    energyspectrum_ = Teuchos::rcp(new std::vector<double>);
    energyspectrum_->resize(wavenumbers_->size());
    // and initialize with zeros, just to be sure
    for (std::size_t rr = 0; rr < energyspectrum_->size(); rr++) (*energyspectrum_)[rr] = 0.0;

    // set size of dissipation-spectrum vector
    dissipationspectrum_ = Teuchos::rcp(new std::vector<double>);
    dissipationspectrum_->resize(wavenumbers_->size());
    // and initialize with zeros, just to be sure
    for (std::size_t rr = 0; rr < dissipationspectrum_->size(); rr++)
      (*dissipationspectrum_)[rr] = 0.0;

    // set size of scalar-variance-spectrum vector
    scalarvariancespectrum_ = Teuchos::rcp(new std::vector<double>);
    scalarvariancespectrum_->resize(wavenumbers_->size());
    // and initialize with zeros, just to be sure
    for (std::size_t rr = 0; rr < scalarvariancespectrum_->size(); rr++)
      (*scalarvariancespectrum_)[rr] = 0.0;

    //-------------------------------------------------
    // initialize remaining variables
    //-------------------------------------------------

    // sum over velocity vector
    sumvel_ = Teuchos::rcp(new std::vector<double>);
    sumvel_->resize(3);

    // sum over squares of velocity vector components
    sumvelvel_ = Teuchos::rcp(new std::vector<double>);
    sumvelvel_->resize(3);

    // allocate some (toggle) vectors
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    toggleu_ = LINALG::CreateVector(*dofrowmap, true);
    togglev_ = LINALG::CreateVector(*dofrowmap, true);
    togglew_ = LINALG::CreateVector(*dofrowmap, true);

    // set number of samples to zero
    numsamp_ = 0;

    // time-step size
    dt_ = params_.get<double>("time step size");

    // get fluid viscosity from material definition
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
    if (id == -1)
      dserror("Could not find Newtonian fluid material");
    else
    {
      const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
      const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
      // we need the kinematic viscosity here
      double dens = actmat->density_;
      visc_ = actmat->viscosity_ / dens;
    }

    //-------------------------------------------------
    // set step for output of energy spectrum
    //-------------------------------------------------
    // decaying homogeneous isotropic turbulence

    if (type_ == decaying_homogeneous_isotropic_turbulence)
    {
      // get number of forcing steps
      int num_forcing_steps = params_.sublist("TURBULENCE MODEL").get<int>("FORCING_TIME_STEPS", 0);

      // determine output steps for energy spectrum of decaying case
      // experimental data are available for
      // t*U_0/M = 42 corresponds to 0+(number of forcing steps) in simulation
      // t*U_0/M = 98 corresponds to 56+(number of forcing steps) in simulation
      // t*U_0/M = 171 corresponds to 129+(number of forcing steps) in simulation
      std::vector<double> times_exp(2);
      times_exp[0] = 56.0;
      times_exp[1] = 129.0;
      // using
      // grid size of experiment M = 0.0508;
      // domain length  L = 10.0 * M;
      // reference length L_ref = L / (2.0 * PI);
      // inlet velocity of experiment U_0 = 10.0;
      // reference time t_ref = 64.0 * M / U_0;
      times_exp[0] /= 64.0;
      times_exp[1] /= 64.0;

      // set steps in vector
      outsteps_ = Teuchos::rcp(new std::vector<int>);
      if (num_forcing_steps != 0) outsteps_->push_back(num_forcing_steps);
      for (std::size_t rr = 0; rr < times_exp.size(); rr++)
      {
        outsteps_->push_back(((int)(times_exp[rr] / dt_)) + num_forcing_steps);
      }
    }
    else
      outsteps_ = Teuchos::null;

    //-------------------------------------------------------------------------
    // initialize output and initially open respective statistics output file
    //-------------------------------------------------------------------------

    Teuchos::RCP<std::ofstream> log_1;
    Teuchos::RCP<std::ofstream> log_2;

    if (discret_->Comm().MyPID() == 0)
    {
      std::string s(statistics_outfilename_);
      s.append(".energy_spectra");

      log_1 = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
      (*log_1) << "# Energy and dissipation spectra for incompressible homogeneous isotropic "
                  "turbulence\n\n\n\n";

      log_1->flush();

      s = statistics_outfilename_;
      s.append(".kinetic_energy");

      log_2 = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
      (*log_2) << "# Evolution of kinetic energy for incompressible homogeneous isotropic "
                  "turbulence\n\n\n";

      log_2->flush();
    }

    //-------------------------------------------------------------------------
    // write resolved turbulent kinetic energy for given discretization
    //-------------------------------------------------------------------------
    if (type_ == decaying_homogeneous_isotropic_turbulence)
      CalculateResolvedEnergyDecayingTurbulence();

    return;
  }


  /*--------------------------------------------------------------*
   | deconstructor                                rasthofer 04/13 |
   *--------------------------------------------------------------*/
  TurbulenceStatisticsHit::~TurbulenceStatisticsHit() { return; }


  /*--------------------------------------------------------------*
   | do sampling                                  rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp)
  {
#if HAVE_FFTW
    //-------------------------------------------------------------------------------------------------
    // calculate energy spectrum via Fourier transformation
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
      (*local_u1)[pos] = (*velnp)[lid];
      // analogously for remaining directions
      lid = discret_->DofRowMap()->LID(dofs[1]);
      (*local_u2)[pos] = (*velnp)[lid];
      lid = discret_->DofRowMap()->LID(dofs[2]);
      (*local_u3)[pos] = (*velnp)[lid];
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

          // get wave number
          const double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);

          // calculate energy
          // E = 1/2 * u_i * conj(u_i)
          // u_i * conj(u_i) = real(u_i)^2 + imag(u_i)^2
          // const std::complex<double> energy = 0.5 * ((*u1_hat)[pos] * conj((*u1_hat)[pos])
          //                                          + (*u2_hat)[pos] * conj((*u2_hat)[pos])
          //                                          + (*u3_hat)[pos] * conj((*u3_hat)[pos]));
          // instead
          const double energy =
              0.5 * (norm((*u1_hat)[pos]) + norm((*u2_hat)[pos]) + norm((*u3_hat)[pos]));

          // insert into sampling vector
          // find position via k
          for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
          {
            if (k > ((*wavenumbers_)[rr] - 0.5) and k <= ((*wavenumbers_)[rr] + 0.5))
            {
              (*energyspectrum_)[rr] += energy;
              // also compute the dissipation spectrum (not yet carefully validated)
              (*dissipationspectrum_)[rr] +=
                  (((*wavenumbers_)[rr]) * ((*wavenumbers_)[rr]) * energy);
            }
          }
        }
      }
    }

    //-------------------------------------------------------------------------------------------------
    // calculate means in physical space
    //-------------------------------------------------------------------------------------------------

    //----------------------------------
    // initialize toggle vectors
    //----------------------------------

    // toggle vectors are one in the position of a dof of this node,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      // get node
      DRT::Node* node = discret_->lRowNode(nn);

      // get global dof ids
      std::vector<int> dof = discret_->Dof(node);
      double one = 1.0;

      // set one in respective position
      toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
      togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
      togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
    }

    // compute squared values of velocity
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    Teuchos::RCP<Epetra_Vector> squaredvelnp = LINALG::CreateVector(*dofrowmap, true);
    squaredvelnp->Multiply(1.0, *velnp, *velnp, 0.0);

    //----------------------------------
    // get values for velocity
    //----------------------------------

    // velocity components
    double u;
    double v;
    double w;
    velnp->Dot(*toggleu_, &u);
    velnp->Dot(*togglev_, &v);
    velnp->Dot(*togglew_, &w);

    // square of velocity components
    double uu;
    double vv;
    double ww;
    squaredvelnp->Dot(*toggleu_, &uu);
    squaredvelnp->Dot(*togglev_, &vv);
    squaredvelnp->Dot(*togglew_, &ww);

    //-------------------------------------------------
    // add spatial mean values to statistical sample
    //-------------------------------------------------

    (*sumvel_)[0] = u / countallnodes;
    (*sumvel_)[1] = v / countallnodes;
    (*sumvel_)[2] = w / countallnodes;

    (*sumvelvel_)[0] = uu / countallnodes;
    (*sumvelvel_)[1] = vv / countallnodes;
    (*sumvelvel_)[2] = ww / countallnodes;

    //----------------------------------------------------------------------
    // increase sample counter
    //----------------------------------------------------------------------

    // for forced case only, since there is not any statistic-stationary state
    // for the decaying case (merely averaging in space)
    if (type_ == forced_homogeneous_isotropic_turbulence) numsamp_++;

    return;
#else
    dserror("FFTW required");
#endif
  }


  /*--------------------------------------------------------------*
   | do sampling                                  rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::DoScatraTimeSample(
      Teuchos::RCP<Epetra_Vector> velnp, Teuchos::RCP<Epetra_Vector> phinp)
  {
#if HAVE_FFTW
    //-------------------------------------------------------------------------------------------------
    // calculate energy spectrum via Fourier transformation
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

    Teuchos::RCP<Teuchos::Array<std::complex<double>>> phi_hat = Teuchos::rcp(
        new Teuchos::Array<std::complex<double>>(nummodes_ * nummodes_ * (nummodes_ / 2 + 1)));
    Teuchos::RCP<Teuchos::Array<double>> local_phi =
        Teuchos::rcp(new Teuchos::Array<double>(nummodes_ * nummodes_ * nummodes_));
    Teuchos::RCP<Teuchos::Array<double>> global_phi =
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
      (*local_u1)[pos] = (*velnp)[lid];
      // analogously for remaining directions
      lid = discret_->DofRowMap()->LID(dofs[1]);
      (*local_u2)[pos] = (*velnp)[lid];
      lid = discret_->DofRowMap()->LID(dofs[2]);
      (*local_u3)[pos] = (*velnp)[lid];
    }

    // set also solution of scalar field

    for (int inode = 0; inode < scatradiscret_->NumMyRowNodes(); inode++)
    {
      // get node
      DRT::Node* node = scatradiscret_->lRowNode(inode);

      // get coordinates
      LINALG::Matrix<3, 1> xyz(true);
      for (int idim = 0; idim < 3; idim++) xyz(idim, 0) = node->X()[idim];

      // get global ids of all dofs of the node
      std::vector<int> dofs = scatradiscret_->Dof(0, node);
      if (dofs.size() > 1) dserror("Only one scatra dof per node expected!");

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
      int lid = scatradiscret_->DofRowMap()->LID(dofs[0]);
      // set value
      (*local_phi)[pos] = (*phinp)[lid];
    }

    // get values form all processors
    // number of nodes without slave nodes
    const int countallnodes = nummodes_ * nummodes_ * nummodes_;
    discret_->Comm().SumAll(&((*local_u1)[0]), &((*global_u1)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u2)[0]), &((*global_u2)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_u3)[0]), &((*global_u3)[0]), countallnodes);

    discret_->Comm().SumAll(&((*local_phi)[0]), &((*global_phi)[0]), countallnodes);

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

    // as well as phi
    fftw_plan fft_4 = fftw_plan_dft_r2c_3d(nummodes_, nummodes_, nummodes_, &((*global_phi)[0]),
        (reinterpret_cast<fftw_complex*>(&((*phi_hat)[0]))), FFTW_ESTIMATE);
    fftw_execute(fft_4);
    // free memory
    fftw_destroy_plan(fft_4);
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
      (*phi_hat)[i] /= nummodes_ * nummodes_ * nummodes_;
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

          // get wave number
          const double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);

          // calculate energy
          // E = 1/2 * u_i * conj(u_i)
          // u_i * conj(u_i) = real(u_i)^2 + imag(u_i)^2
          // const std::complex<double> energy = 0.5 * ((*u1_hat)[pos] * conj((*u1_hat)[pos])
          //                                          + (*u2_hat)[pos] * conj((*u2_hat)[pos])
          //                                          + (*u3_hat)[pos] * conj((*u3_hat)[pos]));
          // instead
          const double energy =
              0.5 * (norm((*u1_hat)[pos]) + norm((*u2_hat)[pos]) + norm((*u3_hat)[pos]));
          const double variance = 0.5 * norm((*phi_hat)[pos]);

          // insert into sampling vector
          // find position via k
          for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
          {
            if (k > ((*wavenumbers_)[rr] - 0.5) and k <= ((*wavenumbers_)[rr] + 0.5))
            {
              (*energyspectrum_)[rr] += energy;
              (*scalarvariancespectrum_)[rr] += variance;
              // also compute the dissipation spectrum (not yet carefully validated)
              (*dissipationspectrum_)[rr] += ((*wavenumbers_)[rr] * (*wavenumbers_)[rr] * energy);
            }
          }
        }
      }
    }

    //-------------------------------------------------------------------------------------------------
    // calculate means in physical space
    //-------------------------------------------------------------------------------------------------

    // add scalar field field here if required

    //----------------------------------
    // initialize toggle vectors
    //----------------------------------

    // toggle vectors are one in the position of a dof of this node,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);

    for (int nn = 0; nn < discret_->NumMyRowNodes(); ++nn)
    {
      // get node
      DRT::Node* node = discret_->lRowNode(nn);

      // get global dof ids
      std::vector<int> dof = discret_->Dof(node);
      double one = 1.0;

      // set one in respective position
      toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
      togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
      togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
    }

    // compute squared values of velocity
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    Teuchos::RCP<Epetra_Vector> squaredvelnp = LINALG::CreateVector(*dofrowmap, true);
    squaredvelnp->Multiply(1.0, *velnp, *velnp, 0.0);

    //----------------------------------
    // get values for velocity
    //----------------------------------

    // velocity components
    double u;
    double v;
    double w;
    velnp->Dot(*toggleu_, &u);
    velnp->Dot(*togglev_, &v);
    velnp->Dot(*togglew_, &w);

    // square of velocity components
    double uu;
    double vv;
    double ww;
    squaredvelnp->Dot(*toggleu_, &uu);
    squaredvelnp->Dot(*togglev_, &vv);
    squaredvelnp->Dot(*togglew_, &ww);

    //-------------------------------------------------
    // add spatial mean values to statistical sample
    //-------------------------------------------------

    (*sumvel_)[0] = u / countallnodes;
    (*sumvel_)[1] = v / countallnodes;
    (*sumvel_)[2] = w / countallnodes;

    (*sumvelvel_)[0] = uu / countallnodes;
    (*sumvelvel_)[1] = vv / countallnodes;
    (*sumvelvel_)[2] = ww / countallnodes;

    //----------------------------------------------------------------------
    // increase sample counter
    //----------------------------------------------------------------------

    numsamp_++;

    return;
#else
    dserror("FFTW required");
#endif
  }


  /*--------------------------------------------------------------*
   | evaluation of dissipation rate and rbvmm-related quantities  |
   |                                              rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::EvaluateResiduals(
      std::map<std::string, Teuchos::RCP<Epetra_Vector>> statevecs)
  {
    dserror("EvaluateResiduals() not yet implemented for hit!");
    return;
  }


  /*--------------------------------------------------------------*
   | dump statistics to file                      rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::DumpStatistics(int step, bool multiple_records)
  {
    //------------------------------
    // compute remaining quantities
    //------------------------------
    // decaying homogeneous isotropic turbulence

    // resolved turbulent kinetic energy from energy spectrum per unit mass
    double q_E = 0.0;
    // resolved turbulent kinetic energy from velocity fluctuations per unit mass
    double q_u = 0.0;
    // mean-flow kinetic energy  per unit mass
    double mke = 0.0;
    // velocity fluctuations
    double u_prime = 0.0;
    // Taylor scale
    double lambda = 0.0;
    // dissipation
    double diss = 0.0;
    // Taylor scale Reynolds number
    double Re_lambda = 0.0;

    if (type_ == decaying_homogeneous_isotropic_turbulence)
    {
      // start for k=1 not k=0, see Diss Hickel
      for (std::size_t rr = 1; rr < energyspectrum_->size(); rr++)
      {
        // build sum up to cut-off wave number
        if ((*wavenumbers_)[rr] <= ((((double)nummodes_) / 2) - 1)) q_E += (*energyspectrum_)[rr];
      }

      for (int rr = 0; rr < 3; rr++)
        q_u += 0.5 * ((*sumvelvel_)[rr] - (*sumvel_)[rr] * (*sumvel_)[rr]);

      for (int rr = 0; rr < 3; rr++) mke += 0.5 * (*sumvel_)[rr] * (*sumvel_)[rr];
    }

    if (type_ == forced_homogeneous_isotropic_turbulence)
    {
      // start for k=1 not k=0, see Diss Hickel
      for (std::size_t rr = 1; rr < energyspectrum_->size(); rr++)
      {
        if ((*wavenumbers_)[rr] <= ((((double)nummodes_) / 2) - 1))
        {
          q_E += (((*energyspectrum_)[rr]) / numsamp_);
          // is this the same? yes, it is!
          // diss += ((*dissipationspectrum_)[rr]/numsamp_);
          diss +=
              ((*wavenumbers_)[rr] * (*wavenumbers_)[rr] * (((*energyspectrum_)[rr]) / numsamp_));
        }
      }

      u_prime = sqrt(2.0 / 3.0 * q_E);
      lambda = sqrt(5.0 * q_E / diss);
      Re_lambda = lambda * u_prime / visc_;
    }

    //------------------------------
    // write results to file
    //------------------------------

    if (discret_->Comm().MyPID() == 0)
    {
      Teuchos::RCP<std::ofstream> log_k;

      std::string s_k(statistics_outfilename_);
      s_k.append(".energy_spectra");

      if (step == 0 and type_ == decaying_homogeneous_isotropic_turbulence)
      {
        log_k = Teuchos::rcp(new std::ofstream(s_k.c_str(), std::ios::app));

        (*log_k) << "# Energy spectrum of initial field (non-dimensionalized form)\n";
        (*log_k) << "#     k              E\n";
        (*log_k) << std::scientific;
        for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
        {
          (*log_k) << " " << std::setw(11) << std::setprecision(4) << (*wavenumbers_)[rr];
          (*log_k) << "     " << std::setw(11) << std::setprecision(4) << (*energyspectrum_)[rr];
          (*log_k) << "\n";

          if ((*wavenumbers_)[rr] >= (3 - 2e-9) and (*wavenumbers_)[rr] <= (3 + 2e-9))
            std::cout << "k " << (*wavenumbers_)[rr] << " energy  " << (*energyspectrum_)[rr]
                      << std::endl;
        }
        (*log_k) << "\n\n\n";
      }
      else
      {
        if (type_ == forced_homogeneous_isotropic_turbulence)
        {
          if (not multiple_records)
          {
            log_k = Teuchos::rcp(new std::ofstream(s_k.c_str(), std::ios::out));
            (*log_k) << "# Energy and dissipation spectra for incompressible homogeneous isotropic "
                        "turbulence\n\n\n";
          }
          else
            log_k = Teuchos::rcp(new std::ofstream(s_k.c_str(), std::ios::app));

          (*log_k) << "# Statistics record ";
          (*log_k) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";
          (*log_k) << "# tke = " << q_E << "    Taylor scale = " << lambda
                   << "    Re_lambda = " << Re_lambda << "\n";
          (*log_k) << std::scientific;
          (*log_k) << "#     k              E              D\n";
          for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
          {
            (*log_k) << " " << std::setw(11) << std::setprecision(4) << (*wavenumbers_)[rr];
            (*log_k) << "     " << std::setw(11) << std::setprecision(4)
                     << (*energyspectrum_)[rr] / numsamp_;
            (*log_k) << "     " << std::setw(11) << std::setprecision(4)
                     << 2 * visc_ * (*dissipationspectrum_)[rr] / numsamp_;
            (*log_k) << "\n";
          }
          (*log_k) << "\n\n\n";
        }
        else
        {
          log_k = Teuchos::rcp(new std::ofstream(s_k.c_str(), std::ios::app));

          bool print = false;
          for (std::size_t rr = 0; rr < outsteps_->size(); rr++)
          {
            if (step == (*outsteps_)[rr])
            {
              print = true;
              break;
            }
          }

          if (print)
          {
            (*log_k) << "# Energy spectrum at time: " << step * dt_ << " \n";
            (*log_k) << "#     k              E              D\n";
            (*log_k) << std::scientific;
            for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
            {
              (*log_k) << " " << std::setw(11) << std::setprecision(4) << (*wavenumbers_)[rr];
              (*log_k) << "     " << std::setw(11) << std::setprecision(4)
                       << (*energyspectrum_)[rr];
              (*log_k) << "     " << std::setw(11) << std::setprecision(4)
                       << 2 * visc_ * (*dissipationspectrum_)[rr];
              (*log_k) << "\n";
            }
            (*log_k) << "\n\n\n";
          }
        }
      }

      log_k->flush();

      if (type_ == decaying_homogeneous_isotropic_turbulence)
      {
        Teuchos::RCP<std::ofstream> log_t;

        std::string s_t(statistics_outfilename_);
        s_t.append(".kinetic_energy");

        log_t = Teuchos::rcp(new std::ofstream(s_t.c_str(), std::ios::app));

        if (step == 0) (*log_t) << "#     t               q(E)          q(u'u')          MKE\n";

        (*log_t) << std::scientific;
        (*log_t) << " " << std::setw(11) << std::setprecision(4) << step * dt_;
        (*log_t) << "     " << std::setw(11) << std::setprecision(4) << q_E;
        (*log_t) << "     " << std::setw(11) << std::setprecision(4) << q_u;
        (*log_t) << "     " << std::setw(11) << std::setprecision(4) << mke;

        (*log_t) << "\n";
        log_t->flush();
      }
    }

    return;
  }


  /*--------------------------------------------------------------*
   | dump statistics to file                      rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::DumpScatraStatistics(int step, bool multiple_records)
  {
    //------------------------------
    // compute remaining quantities
    //------------------------------
    // decaying homogeneous isotropic turbulence

    // resolved turbulent kinetic energy from energy spectrum per unit mass
    double q_E = 0.0;
    // velocity fluctuations
    double u_prime = 0.0;
    // Taylor scale
    double lambda = 0.0;
    // dissipation
    double diss = 0.0;
    // Taylor scale Reynolds number
    double Re_lambda = 0.0;

    if (type_ == forced_homogeneous_isotropic_turbulence)
    {
      for (std::size_t rr = 1; rr < energyspectrum_->size(); rr++)
      {
        if ((*wavenumbers_)[rr] <= ((((double)nummodes_) / 2) - 1))
        {
          q_E += ((*energyspectrum_)[rr] / numsamp_);
          // is this the same? yes, it is!
          // diss += ((*dissipationspectrum_)[rr]/numsamp_);
          diss += ((*wavenumbers_)[rr] * (*wavenumbers_)[rr] * ((*energyspectrum_)[rr] / numsamp_));
        }
      }

      u_prime = sqrt(2.0 / 3.0 * q_E);
      lambda = sqrt(5.0 * q_E / diss);
      Re_lambda = lambda * u_prime / visc_;
    }

    //------------------------------
    // write results to file
    //------------------------------

    if (discret_->Comm().MyPID() == 0)
    {
      Teuchos::RCP<std::ofstream> log_k;

      std::string s_k = statistics_outfilename_;
      s_k.append(".energy_spectra");

      //    {
      //      if (type_ == forced_homogeneous_isotropic_turbulence)
      //      {
      if (not multiple_records)
      {
        log_k = Teuchos::rcp(new std::ofstream(s_k.c_str(), std::ios::out));
        (*log_k) << "# Energy and dissipation spectra for incompressible homogeneous isotropic "
                    "turbulence\n\n\n";
      }
      else
        log_k = Teuchos::rcp(new std::ofstream(s_k.c_str(), std::ios::app));

      (*log_k) << "# Statistics record ";
      (*log_k) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";
      (*log_k) << "# tke = " << q_E << "    Taylor scale = " << lambda
               << "    Re_lambda = " << Re_lambda << "\n";
      (*log_k) << std::scientific;
      (*log_k) << "#     k              E              E_phi\n";
      for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
      {
        (*log_k) << " " << std::setw(11) << std::setprecision(4) << (*wavenumbers_)[rr];
        (*log_k) << "     " << std::setw(11) << std::setprecision(4)
                 << (*energyspectrum_)[rr] / numsamp_;
        (*log_k) << "     " << std::setw(11) << std::setprecision(4)
                 << (*scalarvariancespectrum_)[rr] / numsamp_;
        //          (*log_k) << "     " << std::setw(11) << std::setprecision(4) <<
        //          2*visc_*(*dissipationspectrum_)[rr]/numsamp_;
        (*log_k) << "\n";
      }
      (*log_k) << "\n\n\n";
      //      }
      //    }

      log_k->flush();
    }

    return;
  }


  /*--------------------------------------------------------------*
   | reset statistics to zero                     rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::ClearStatistics()
  {
    for (std::size_t rr = 0; rr < energyspectrum_->size(); rr++)
    {
      (*energyspectrum_)[rr] = 0.0;
      (*dissipationspectrum_)[rr] = 0.0;
      (*scalarvariancespectrum_)[rr] = 0.0;
    }

    for (std::size_t rr = 0; rr < sumvel_->size(); rr++)
    {
      (*sumvel_)[rr] = 0.0;
      (*sumvelvel_)[rr] = 0.0;
    }

    numsamp_ = 0;

    return;
  }


  /*--------------------------------------------------------------*
   | reset statistics to zero                     rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::ClearScatraStatistics()
  {
    for (std::size_t rr = 0; rr < energyspectrum_->size(); rr++)
    {
      (*energyspectrum_)[rr] = 0.0;
      (*dissipationspectrum_)[rr] = 0.0;
      (*scalarvariancespectrum_)[rr] = 0.0;
    }

    for (std::size_t rr = 0; rr < sumvel_->size(); rr++)
    {
      (*sumvel_)[rr] = 0.0;
      (*sumvelvel_)[rr] = 0.0;
    }

    numsamp_ = 0;

    return;
  }


  /*--------------------------------------------------------------*
   | calculate the resolved energy for the given discretization   |
   |                                              rasthofer 04/13 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHit::CalculateResolvedEnergyDecayingTurbulence()
  {
    //----------------------------------------
    // set reference values
    //----------------------------------------

    // non-dimensionalize wave number (according to Collis 2002)
    // grid size of experiment
    const double M = 0.0508;
    // domain length
    const double L = 10.0 * M;
    // reference length
    const double L_ref = L / (2.0 * PI);
    // non-dimensionalize energy spectrum
    // inlet velocity of experiment
    const double U_0 = 10.0;
    // reference time
    const double t_ref = 64.0 * M / U_0;

    //----------------------------------------
    // set-up wave numbers
    //----------------------------------------

    // wave number given from experiment [cm-1]
    // have to be transferred to [m-1] and non-dimensionalized
    std::vector<double> k_exp(20);
    k_exp[0] = 0.15;
    k_exp[1] = 0.20;
    k_exp[2] = 0.25;
    k_exp[3] = 0.30;
    k_exp[4] = 0.40;
    k_exp[5] = 0.50;
    k_exp[6] = 0.70;
    k_exp[7] = 1.00;
    k_exp[8] = 1.50;
    k_exp[9] = 2.00;
    k_exp[10] = 2.50;
    k_exp[11] = 3.00;
    k_exp[12] = 4.00;
    k_exp[13] = 6.00;
    k_exp[14] = 8.00;
    k_exp[15] = 10.00;
    k_exp[16] = 12.50;
    k_exp[17] = 15.00;
    k_exp[18] = 17.50;
    k_exp[19] = 20.00;

    for (std::size_t rr = 0; rr < k_exp.size(); rr++) k_exp[rr] *= (L_ref / 0.01);

    //----------------------------------------
    // set-up energy
    //----------------------------------------

    // energy spectrum given from experiment [cm3/s2]
    // have to be transferred to [m3/s2] and non-dimensionalized
    // first measured location
    std::vector<double> E_exp_1(20);
    E_exp_1[0] = -1.0;
    E_exp_1[1] = 129.0;
    E_exp_1[2] = 230.0;
    E_exp_1[3] = 322.0;
    E_exp_1[4] = 435.0;
    E_exp_1[5] = 457.0;
    E_exp_1[6] = 380.0;
    E_exp_1[7] = 270.0;
    E_exp_1[8] = 168.0;
    E_exp_1[9] = 120.0;
    E_exp_1[10] = 89.0;
    E_exp_1[11] = 70.3;
    E_exp_1[12] = 47.0;
    E_exp_1[13] = 24.7;
    E_exp_1[14] = 12.6;
    E_exp_1[15] = 7.42;
    E_exp_1[16] = 3.96;
    E_exp_1[17] = 2.33;
    E_exp_1[18] = 1.34;
    E_exp_1[19] = 0.80;

    // second measured location
    std::vector<double> E_exp_2(20);
    E_exp_2[0] = -1.0;
    E_exp_2[1] = 106.0;
    E_exp_2[2] = 196.0;
    E_exp_2[3] = 195.0;
    E_exp_2[4] = 202.0;
    E_exp_2[5] = 168.0;
    E_exp_2[6] = 127.0;
    E_exp_2[7] = 79.2;
    E_exp_2[8] = 47.8;
    E_exp_2[9] = 34.6;
    E_exp_2[10] = 28.6;
    E_exp_2[11] = 23.1;
    E_exp_2[12] = 14.3;
    E_exp_2[13] = 5.95;
    E_exp_2[14] = 2.23;
    E_exp_2[15] = 0.900;
    E_exp_2[16] = 0.363;
    E_exp_2[17] = 0.162;
    E_exp_2[18] = 0.0660;
    E_exp_2[19] = 0.0330;

    // third measured location
    std::vector<double> E_exp_3(20);
    E_exp_3[0] = 49.7;
    E_exp_3[1] = 92.0;
    E_exp_3[2] = 120;
    E_exp_3[3] = 125;
    E_exp_3[4] = 98.0;
    E_exp_3[5] = 81.5;
    E_exp_3[6] = 60.2;
    E_exp_3[7] = 39.4;
    E_exp_3[8] = 24.1;
    E_exp_3[9] = 16.5;
    E_exp_3[10] = 12.5;
    E_exp_3[11] = 9.12;
    E_exp_3[12] = 5.62;
    E_exp_3[13] = 1.69;
    E_exp_3[14] = 0.520;
    E_exp_3[15] = 0.161;
    E_exp_3[16] = 0.0520;
    E_exp_3[17] = 0.0141;
    E_exp_3[18] = -1.0;
    E_exp_3[19] = -1.0;

    for (std::size_t rr = 0; rr < E_exp_1.size(); rr++)
    {
      E_exp_1[rr] *= ((0.01 * 0.01 * 0.01) * (t_ref * t_ref) / (L_ref * L_ref * L_ref));
      E_exp_2[rr] *= ((0.01 * 0.01 * 0.01) * (t_ref * t_ref) / (L_ref * L_ref * L_ref));
      E_exp_3[rr] *= ((0.01 * 0.01 * 0.01) * (t_ref * t_ref) / (L_ref * L_ref * L_ref));
    }

    //----------------------------------------
    // get cut-off wave number
    //----------------------------------------

    // determine cut-off wave number
    // via number of elements in each spatial direction
    // k_c = Pi/h, where h=2Pi/nele (domain (2Pi)^3 assumed)
    // note: nele = nummodes_
    // cut-off wave number
    const double k_c = ((double)nummodes_) / 2.0 - 1;

    //---------------------------------------------------------------
    // calculate resolved and subgrid-scale turbulent kinetic energy
    //---------------------------------------------------------------

    // variables for total turbulent kinetic energy at
    // the three measure locations
    double q_1 = 0.0;
    double q_2 = 0.0;
    double q_3 = 0.0;

    // variables for resolved turbulent kinetic energy at
    // the three measure locations
    double q_1_r = 0.0;
    double q_2_r = 0.0;
    double q_3_r = 0.0;

    // integration of energy spectrum via trapezoidal rule
    int rr_k_max = 0.0;
#if 0
  for (std::size_t rr = 1; rr < E_exp_1.size(); rr++)
  {
    std::cout << "rr " << rr << std::endl;
    if (E_exp_1[rr-1] > 0.0)
    {
      std::cout << "1 " << 0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_1[rr] + E_exp_1[rr-1]) << std::endl;
      std::cout << "k_exp[rr] " << k_exp[rr] << std::endl;
      std::cout << "k_exp[rr-1] " << k_exp[rr-1] << std::endl;
      std::cout << "E_exp_1[rr] " << E_exp_1[rr] << std::endl;
      std::cout << "E_exp_1[rr-1] " << E_exp_1[rr-1] << std::endl;
      q_1 += (0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_1[rr] + E_exp_1[rr-1]));
      std::cout << q_1 << std::endl;
      if (k_exp[rr] <= k_c)
      {
        std::cout << "rr for k_c " << rr << "   " << k_exp[rr] << std::endl;
        q_1_r += 0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_1[rr] + E_exp_1[rr-1]);
        rr_k_max = rr;
      }
    }
    if (E_exp_2[rr-1] > 0.0)
    {
      q_2 += 0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_2[rr] + E_exp_2[rr-1]);
      if (k_exp[rr] <= (k_c))
        q_2_r += 0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_2[rr] + E_exp_2[rr-1]);
    }
    if (E_exp_2[rr] > 0.0)
    {
      q_3 += 0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_3[rr] + E_exp_3[rr-1]);
      if (k_exp[rr] <= (k_c))
        q_3_r += 0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_3[rr] + E_exp_3[rr-1]);
    }
  }
#endif
    for (std::size_t rr = 1; rr < E_exp_1.size(); rr++)
    {
      if (E_exp_1[rr - 1] > 0.0)
      {
        q_1 += IntegrateTrapezoidalRule(k_exp[rr - 1], k_exp[rr], E_exp_1[rr - 1], E_exp_1[rr]);
        if (k_exp[rr] <= k_c)
        {
          q_1_r += IntegrateTrapezoidalRule(k_exp[rr - 1], k_exp[rr], E_exp_1[rr - 1], E_exp_1[rr]);
          rr_k_max = rr;
        }
      }
      if (E_exp_2[rr - 1] > 0.0)
      {
        q_2 += IntegrateTrapezoidalRule(k_exp[rr - 1], k_exp[rr], E_exp_2[rr - 1], E_exp_2[rr]);
        if (k_exp[rr] <= k_c)
          q_2_r += IntegrateTrapezoidalRule(k_exp[rr - 1], k_exp[rr], E_exp_2[rr - 1], E_exp_2[rr]);
      }
      if (E_exp_2[rr] > 0.0)
      {
        q_3 += IntegrateTrapezoidalRule(k_exp[rr - 1], k_exp[rr], E_exp_3[rr - 1], E_exp_3[rr]);
        if (k_exp[rr] <= k_c)
          q_3_r += IntegrateTrapezoidalRule(k_exp[rr - 1], k_exp[rr], E_exp_3[rr - 1], E_exp_3[rr]);
      }
    }
    // integrate resolved turbulent kinetic energy from k_max to k_c,
    // which is missing above
    if (rr_k_max < (int)(E_exp_1.size() - 1))
    {
      const double E1_k_c = Interpolate(
          k_c, k_exp[rr_k_max], k_exp[rr_k_max + 1], E_exp_1[rr_k_max], E_exp_1[rr_k_max + 1]);
      q_1_r += IntegrateTrapezoidalRule(k_exp[rr_k_max], k_c, E_exp_1[rr_k_max], E1_k_c);

      const double E2_k_c = Interpolate(
          k_c, k_exp[rr_k_max], k_exp[rr_k_max + 1], E_exp_2[rr_k_max], E_exp_2[rr_k_max + 1]);
      q_2_r += IntegrateTrapezoidalRule(k_exp[rr_k_max], k_c, E_exp_2[rr_k_max], E2_k_c);

      if (E_exp_3[rr_k_max + 1] > 0.0 and E_exp_3[rr_k_max] > 0.0)
      {
        const double E3_k_c = Interpolate(
            k_c, k_exp[rr_k_max], k_exp[rr_k_max + 1], E_exp_3[rr_k_max], E_exp_3[rr_k_max + 1]);
        q_3_r += IntegrateTrapezoidalRule(k_exp[rr_k_max], k_c, E_exp_3[rr_k_max], E3_k_c);
      }
    }
#if 0
  // integration of energy spectrum via Simpson rule
  for (std::size_t rr = 2; rr < E_exp_1.size(); rr++)
  {
    std::cout << "rr " << rr << std::endl;
    if (E_exp_1[rr-1] > 0.0 and E_exp_1[rr-2] > 0.0)
    {
      //std::cout << "1 " << 0.5 * (k_exp[rr] - k_exp[rr-1]) * (E_exp_1[rr] + E_exp_1[rr-1]) << std::endl;
      std::cout << "k_exp[rr] " << k_exp[rr] << std::endl;
      std::cout << "k_exp[rr-1] " << k_exp[rr-1] << std::endl;
      std::cout << "k_exp[rr-2] " << k_exp[rr-2] << std::endl;
      std::cout << "E_exp_1[rr] " << E_exp_1[rr] << std::endl;
      std::cout << "E_exp_1[rr-1] " << E_exp_1[rr-1] << std::endl;
      std::cout << "E_exp_1[rr-2] " << E_exp_1[rr-2] << std::endl;
      q_1 += (1.0/6.0 * (k_exp[rr] - k_exp[rr-2]) * (E_exp_1[rr] + 4 * E_exp_1[rr-1] + E_exp_1[rr-2]));
      std::cout << q_1 << std::endl;
      if (k_exp[rr] <= k_c)
        q_1_r += (1.0/6.0 * (k_exp[rr] - k_exp[rr-2]) * (E_exp_1[rr] + 4 * E_exp_1[rr-1] + E_exp_1[rr-2]));
      rr++;
    }
  }
#endif

    // subgrid turbulent kinetic energy at
    // the three measure locations
    const double q_1_sgs = q_1 - q_1_r;
    const double q_2_sgs = q_2 - q_2_r;
    const double q_3_sgs = q_3 - q_3_r;

    //------------------------------
    // write results to file
    //------------------------------

    Teuchos::RCP<std::ofstream> log;
    if (discret_->Comm().MyPID() == 0)
    {
      std::string s = statistics_outfilename_;
      s.append(".kinetic_energy");

      log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::app));
      (*log)
          << "# Turbulent kinetic energy from experimental spectrum (non-dimensionalized form)\n";
      (*log) << "#     t           q (total)      q (resolved)    q (subgrid)\n";
      (*log) << std::scientific;
      (*log) << " " << std::setw(11) << std::setprecision(4) << 42.0 * M / U_0 / t_ref;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_1;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_1_r;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_1_sgs;
      (*log) << "\n";
      (*log) << " " << std::setw(11) << std::setprecision(4) << 98.0 * M / U_0 / t_ref;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_2;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_2_r;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_2_sgs;
      (*log) << "\n";
      (*log) << " " << std::setw(11) << std::setprecision(4) << 171.0 * M / U_0 / t_ref;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_3;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_3_r;
      (*log) << "     " << std::setw(11) << std::setprecision(4) << q_3_sgs;
      (*log) << "\n\n\n";
      log->flush();
    }

    return;
  }

  /*--------------------------------------------------------------*
   | constructor                                         bk 03/15 |
   *--------------------------------------------------------------*/
  TurbulenceStatisticsHitHDG::TurbulenceStatisticsHitHDG(Teuchos::RCP<DRT::Discretization> actdis,
      Teuchos::ParameterList& params, const std::string& statistics_outfilename, const bool forced)
      : TurbulenceStatisticsHit(actdis, params, statistics_outfilename, forced)
  {
    //-----------------------------------
    // output to screen
    //-----------------------------------

    if (discret_->Comm().MyPID() == 0)
    {
      std::cout << "This is the turbulence statistics manager for HDG\n" << std::endl;
    }

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
    energyspectrum_ = Teuchos::rcp(new std::vector<double>);
    energyspectrum_->resize(wavenumbers_->size());
    // and initialize with zeros, just to be sure
    for (std::size_t rr = 0; rr < energyspectrum_->size(); rr++) (*energyspectrum_)[rr] = 0.0;

    // set size of dissipation-spectrum vector
    dissipationspectrum_ = Teuchos::rcp(new std::vector<double>);
    dissipationspectrum_->resize(wavenumbers_->size());
    // and initialize with zeros, just to be sure
    for (std::size_t rr = 0; rr < dissipationspectrum_->size(); rr++)
      (*dissipationspectrum_)[rr] = 0.0;

    // set size of scalar-variance-spectrum vector
    scalarvariancespectrum_ = Teuchos::rcp(new std::vector<double>);
    scalarvariancespectrum_->resize(wavenumbers_->size());
    // and initialize with zeros, just to be sure
    for (std::size_t rr = 0; rr < scalarvariancespectrum_->size(); rr++)
      (*scalarvariancespectrum_)[rr] = 0.0;


    //-------------------------------------------------------------------------
    // initialize output and initially open respective statistics output file
    //-------------------------------------------------------------------------

    Teuchos::RCP<std::ofstream> log_1;
    Teuchos::RCP<std::ofstream> log_2;

    if (discret_->Comm().MyPID() == 0)
    {
      std::string s(statistics_outfilename_);
      s.append(".energy_spectra");

      log_1 = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::app));
      (*log_1) << "# Using 5 points in every element for FFT \n\n\n\n";

      log_1->flush();

      s = statistics_outfilename_;
      s.append(".kinetic_energy");

      log_2 = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::app));
      (*log_2) << "# Use HDG-trace variables only \n\n\n";

      log_2->flush();
    }

    return;
  }

  /*--------------------------------------------------------------*
   | deconstructor                                       bk 03/15 |
   *--------------------------------------------------------------*/
  TurbulenceStatisticsHitHDG::~TurbulenceStatisticsHitHDG() { return; }

  /*--------------------------------------------------------------*
   | do sampling                                         bk 03/15 |
   *--------------------------------------------------------------*/
  void TurbulenceStatisticsHitHDG::DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp)
  {
#if HAVE_FFTW
    //-------------------------------------------------------------------------------------------------
    // calculate energy spectrum via Fourier transformation
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

    discret_->SetState(1, "intvelnp", velnp);

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

          // get wave number
          const double k = sqrt(k_1 * k_1 + k_2 * k_2 + k_3 * k_3);

          // calculate energy
          // E = 1/2 * u_i * conj(u_i)
          // u_i * conj(u_i) = real(u_i)^2 + imag(u_i)^2
          // const std::complex<double> energy = 0.5 * ((*u1_hat)[pos] * conj((*u1_hat)[pos])
          //                                          + (*u2_hat)[pos] * conj((*u2_hat)[pos])
          //                                          + (*u3_hat)[pos] * conj((*u3_hat)[pos]));
          // instead
          const double energy =
              0.5 * (norm((*u1_hat)[pos]) + norm((*u2_hat)[pos]) + norm((*u3_hat)[pos]));

          // insert into sampling vector
          // find position via k
          for (std::size_t rr = 0; rr < wavenumbers_->size(); rr++)
          {
            if (k > ((*wavenumbers_)[rr] - 0.5) and k <= ((*wavenumbers_)[rr] + 0.5))
            {
              (*energyspectrum_)[rr] += energy;
              // also compute the dissipation spectrum (not yet carefully validated)
              (*dissipationspectrum_)[rr] +=
                  (((*wavenumbers_)[rr]) * ((*wavenumbers_)[rr]) * energy);
            }
          }
        }
      }
    }

    //-------------------------------------------------------------------------------------------------
    // calculate means in physical space
    //-------------------------------------------------------------------------------------------------

    //----------------------------------
    // initialize toggle vectors
    //----------------------------------

    // toggle vectors are one in the position of a dof of this node,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);

    // this is currently not supported and has to be adapted to HDG
    //  for (int nn=0; nn<discret_->NumMyRowNodes(); ++nn)
    //  {
    //    // get node
    //    DRT::Node* node = discret_->lRowNode(nn);
    //
    //    // get global dof ids
    //    std::vector<int> dof = discret_->Dof(node);
    //    double one = 1.0;
    //
    //    // set one in respective position
    //    toggleu_->ReplaceGlobalValues(1,&one,&(dof[0]));
    //    togglev_->ReplaceGlobalValues(1,&one,&(dof[1]));
    //    togglew_->ReplaceGlobalValues(1,&one,&(dof[2]));
    //  }
    //
    //  // compute squared values of velocity
    //  const Epetra_Map* dofrowmap = discret_->DofRowMap();
    //  Teuchos::RCP<Epetra_Vector> squaredvelnp = LINALG::CreateVector(*dofrowmap,true);
    //  squaredvelnp->Multiply(1.0,*velnp,*velnp,0.0);
    //
    //  //----------------------------------
    //  // get values for velocity
    //  //----------------------------------
    //
    //  // velocity components
    //  double u;
    //  double v;
    //  double w;
    //  velnp->Dot(*toggleu_,&u);
    //  velnp->Dot(*togglev_,&v);
    //  velnp->Dot(*togglew_,&w);
    //
    //  // square of velocity components
    //  double uu;
    //  double vv;
    //  double ww;
    //  squaredvelnp->Dot(*toggleu_,&uu);
    //  squaredvelnp->Dot(*togglev_,&vv);
    //  squaredvelnp->Dot(*togglew_,&ww);
    //
    //  //-------------------------------------------------
    //  // add spatial mean values to statistical sample
    //  //-------------------------------------------------
    //
    //  (*sumvel_)[0] = u/countallnodes;
    //  (*sumvel_)[1] = v/countallnodes;
    //  (*sumvel_)[2] = w/countallnodes;
    //
    //  (*sumvelvel_)[0] = uu/countallnodes;
    //  (*sumvelvel_)[1] = vv/countallnodes;
    //  (*sumvelvel_)[2] = ww/countallnodes;


    //----------------------------------------------------------------------
    // increase sample counter
    //----------------------------------------------------------------------

    // for forced case only, since there is not any statistic-stationary state
    // for the decaying case (merely averaging in space)
    if (type_ == forced_homogeneous_isotropic_turbulence) numsamp_++;
    discret_->ClearState(true);
    return;
#else
    dserror("FFTW required");
#endif
  }


}  // end namespace FLD
