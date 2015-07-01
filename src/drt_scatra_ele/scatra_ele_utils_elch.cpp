/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_utils_elch.cpp

\brief utility class supporting element evaluation for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_utils_elch.H"

#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_elch.H"
#include "../headers/definitions.h"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 07/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleUtilsElch<distype>* DRT::ELEMENTS::ScaTraEleUtilsElch<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    bool create
    )
{
  static ScaTraEleUtilsElch<distype>* instance;

  if(create)
  {
    if(instance == NULL)
      instance = new ScaTraEleUtilsElch<distype>(numdofpernode,numscal);
  }

  else if(instance != NULL)
  {
    delete instance;
    instance = NULL;
  }

  return instance;
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElch<distype>::Done()
{
  // delete singleton
  Instance(0,0,false);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 07/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleUtilsElch<distype>::ScaTraEleUtilsElch(const int numdofpernode,const int numscal) :
  numdofpernode_(numdofpernode),
  numscal_(numscal)
{
  return;
}


/*---------------------------------------------------------------------------------------------------------*
 | evaluation of electrochemistry kinetics at integration point on domain or boundary element   fang 07/15 |
 *---------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleUtilsElch<distype>::EvaluateElchKineticsAtIntegrationPoint(
    const DRT::Element*                           ele,        ///< current element
    Epetra_SerialDenseMatrix&                     emat,       ///< element matrix
    Epetra_SerialDenseVector&                     erhs,       ///< element right-hand side vector
    const std::vector<LINALG::Matrix<nen_,1> >&   conreact,   ///< concentration values of reactive species at element nodes
    const LINALG::Matrix<nen_,1>&                 pot,        ///< el. potential values at element nodes
    const LINALG::Matrix<nen_,1>&                 phihist,    ///< element history vector for potential at electrode
    const double                                  timefac,    ///< time factor
    const double                                  fac,        ///< Gauss integration factor
    const LINALG::Matrix<nen_,1>&                 funct,      ///< shape functions at int. point
    const Teuchos::RCP<const MAT::Material>&      material,   ///< material
    const Teuchos::RCP<DRT::Condition>&           cond,       ///< condition
    const int                                     nume,       ///< number of transferred electrons
    const std::vector<int>&                       stoich,     ///< stoichiometry of the reaction
    const double                                  valence_k,  ///< valence of the single reactant
    const int                                     kinetics,   ///< desired electrode kinetics model
    const double                                  pot0,       ///< actual electrode potential on metal side
    const double                                  frt,        ///< factor F/RT
    const double                                  fns,        ///< factor fns = s_k / (nume * faraday * (-1))
    const double                                  scalar,     ///< scaling factor for element matrix and right-hand side vector contributions
    const int                                     k           ///< index of evaluated scalar
) const
{
  // for pre-multiplication of i0 with 1/(F z_k)
  const double faraday = INPAR::ELCH::faraday_const;    // unit of F: C/mol or mC/mmol or µC/µmol

  // concentration of active species at integration point
  std::vector<double> conint(numscal_,0.0);

  // concentration is evaluated at all GP since some reaction models depend on all concentrations
  for(int kk=0;kk<numscal_;++kk)
    conint[kk] = funct.Dot(conreact[kk]);

  // el. potential at integration point
  const double potint = funct.Dot(pot);

  // history of potential phi on electrode boundary at integration point
  const double phihistint = funct.Dot(phihist);

  // electrode potential difference (epd) at integration point
  const double epd = (pot0 - potint);

  // concentration-dependent Butler-Volmer law(s)
  switch(kinetics)
  {
  case INPAR::SCATRA::butler_volmer:
  case INPAR::SCATRA::butler_volmer_yang1997:
  {
    // read model-specific parameters
    const double alphaa = cond->GetDouble("alpha_a");
    const double alphac = cond->GetDouble("alpha_c");
    const double dlcap = cond->GetDouble("dl_spec_cap");
    double pot0hist = 0.0;
    if(dlcap!=0.0)
      pot0hist = cond->GetDouble("pot0hist");
    double       i0 = cond->GetDouble("i0");
    if (i0 < -EPS14) dserror("i0 is negative, \n"
                             "a positive definition is necessary due to generalized reaction models: %f",i0);
    // add time factor
    i0*=timefac;
    const double gamma = cond->GetDouble("gamma");
    const double refcon = cond->GetDouble("refcon");
    if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

    if(valence_k!=nume)
      dserror("Kinetic model Butler-Volmer: The number of transferred electrons need to be  \n "
          "the same as the charge number of the reacting species %i", k);

    // open circuit potential is assumed to be zero
    const double ocp = 0.0;

    // overpotential based on opencircuit potential
    const double eta = epd - ocp;

#if 0
    // print all parameters read from the current condition
    std::cout<<"kinetic model  = "<<*kinetics<<std::endl;
    std::cout<<"react. species = "<<speciesid<<std::endl;
    std::cout<<"pot0(mod.)     = "<<pot0<<std::endl;
    std::cout<<"curvenum       = "<<curvenum<<std::endl;
    std::cout<<"alpha_a        = "<<alphaa<<std::endl;
    std::cout<<"alpha_c        = "<<alphac<<std::endl;
    std::cout<<"i0(mod.)       = "<<i0<<std::endl;
    std::cout<<"gamma          = "<<gamma<<std::endl;
    std::cout<<"refcon         = "<<refcon<<std::endl;
    std::cout<<"F/RT           = "<<frt<<std::endl<<std::endl;
    std::cout<<"time factor    = "<<timefac<<std::endl;
    std::cout<<"alpha_F        = "<<alphaF<<std::endl;
#endif

#ifdef DEBUG
    // some safety checks/ user warnings
    if ((alphaa*frt*eta) > 100.0)
      std::cout<<"WARNING: Exp(alpha_a...) in Butler-Volmer law is near overflow!"
      <<exp(alphaa*frt*eta)<<std::endl;
    if (((-alphac)*frt*eta) > 100.0)
      std::cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
      <<exp((-alphac)*frt*eta)<<std::endl;
#endif
    double pow_conint_gamma_k = 0.0;
    if ((conint[k]/refcon) < EPS13)
    {
      pow_conint_gamma_k = std::pow(EPS13,gamma);
#ifdef DEBUG
      std::cout<<"WARNING: Rel. Conc. in Butler-Volmer formula is zero/negative: "<<(conint[k]/refcon)<<std::endl;
      std::cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<std::endl;
#endif
    }
    else
      pow_conint_gamma_k = std::pow(conint[k]/refcon,gamma);

    if (kinetics==INPAR::SCATRA::butler_volmer)
    {
      // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
      const double expterm = exp(alphaa*frt*eta)-exp((-alphac)*frt*eta);

      double concterm = 0.0;
      if (conint[k] > EPS13)
        concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
      else
        concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

      for (int vi=0; vi<nen_; ++vi)
      {
        const double fac_fns_i0_funct_vi = scalar*fac*fns*i0*funct(vi);

        // ------matrix: d(R_k)/dx = d(theta*dt*(-1)*(w_k,j_k))/dx
        for (int ui=0; ui<nen_; ++ui)
        {
          emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*concterm*funct(ui)*expterm;
          emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*pow_conint_gamma_k*(((-alphaa)*frt*exp(alphaa*frt*eta))+((-alphac)*frt*exp((-alphac)*frt*eta)))*funct(ui);
        }

        // -----right-hand-side: -R_k = -theta*dt*(-1)*(w_k,j_k)
        erhs[vi*numdofpernode_+k] -= -fac_fns_i0_funct_vi*pow_conint_gamma_k*expterm;
      }

      if(dlcap!=0.0)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          //TODO: Do we need epsilon here
          //add terms of double layer capacitance current density
          for (int ui=0; ui<nen_; ++ui)
            emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += fac*funct(vi)*funct(ui)*dlcap/(nume*faraday);

          // -----right-hand-side: -R_k = -(theta*dt*(-1)*(w_k,j_k)
          erhs[vi*numdofpernode_+k] += fac*funct(vi)*dlcap/(nume*faraday)*(phihistint-pot0hist-potint+pot0);
        }
      }
    } // end if(kinetics=="Butler-Volmer")

    else if (kinetics==INPAR::SCATRA::butler_volmer_yang1997)
    {
      if(dlcap!=0.0)
        dserror("double layer charging is not implemented for Butler-Volmer-Yang electrode kinetics");

      // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
      double concterm = 0.0;
      if ((conint[k]/refcon) > EPS13)
        concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
      else
        concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

      for (int vi=0; vi<nen_; ++vi)
      {
        const double fac_fns_i0_funct_vi = scalar*fac*fns*i0*funct(vi);
        // ------matrix: d(R_k)/dx = d(theta*dt*(-1)*(w_k,j_k))/dx
        for (int ui=0; ui<nen_; ++ui)
        {
          emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*funct(ui)*(-(concterm*exp((-alphac)*frt*eta)));
          emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*(((-alphaa)*frt*exp(alphaa*frt*eta))+(pow_conint_gamma_k*(-alphac)*frt*exp((-alphac)*frt*eta)))*funct(ui);
        }

        // -----right-hand-side: -R_k = -theta*dt*(-1)*(w_k,j_k)
        erhs[vi*numdofpernode_+k] -= -fac_fns_i0_funct_vi*(exp(alphaa*frt*eta)-(pow_conint_gamma_k*exp((-alphac)*frt*eta)));
      }
    } // if (kinetics=="Butler-Volmer-Yang1997")
    else
      dserror("You should not be here!! Two options: Butler-Volmer-Yang1997 and Butler-Volmer-Yang1997 ");

    break;
  }

  // Tafel law (see phd-thesis Georg Bauer, pp.25):
  // implementation of cathodic path: i_n = i_0 * (-exp(-alpha * frt* eta)
  // -> cathodic reaction path: i_0 > 0 and alpha > 0
  // -> anodic reacton path:    i_0 < 0 and alpha < 0
  case INPAR::SCATRA::tafel:
  {
    // read model-specific parameter
    const double alpha = cond->GetDouble("alpha");
    double       i0 = cond->GetDouble("i0");
    i0*=timefac;
    const double gamma = cond->GetDouble("gamma");
    const double refcon = cond->GetDouble("refcon");
    if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);
    const double dlcap = cond->GetDouble("dl_spec_cap");
    if(dlcap!=0.0) dserror("double layer charging is not implemented for Tafel electrode kinetics");

    if(valence_k!=nume)
      dserror("Kinetic model Butler-Volmer: The number of transferred electrons need to be  \n "
          "the same as the charge number of the reacting species %i", k);
    // opencircuit potential is assumed to be zero
    const double ocp = 0.0;
    // overpotential based on opencircuit potential
    const double eta = epd - ocp;

    // concentration-dependent Tafel law
    double pow_conint_gamma_k(0.0);

#ifdef DEBUG
    // some safety checks/ user warnings
    if (((-alpha)*frt*eta) > 100.0)
      std::cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
        <<exp((-alpha)*frt*eta)<<std::endl;
#endif
    if ((conint[k]/refcon) < EPS13)
    {
      pow_conint_gamma_k = std::pow(EPS13,gamma);
#ifdef DEBUG
      std::cout<<"WARNING: Rel. Conc. in Tafel formula is zero/negative: "<<(conint[k]/refcon)<<std::endl;
      std::cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<std::endl;
#endif
    }
    else
      pow_conint_gamma_k = std::pow(conint[k]/refcon,gamma);

    const double expterm = -exp((-alpha)*frt*eta);

    double concterm = 0.0;
    // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
    if (conint[k] > EPS13)
      concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
    else
      concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

    for (int vi=0; vi<nen_; ++vi)
    {
      const double fac_fns_i0_funct_vi = scalar*fac*fns*i0*funct(vi);
      // ---------------------matrix
      for (int ui=0; ui<nen_; ++ui)
      {
        emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*concterm*funct(ui)*expterm;
        emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*pow_conint_gamma_k*(-alpha)*frt*exp((-alpha)*frt*eta)*funct(ui); // do not forget the (-1) from differentiation of eta!
      }
      // ------------right-hand-side
      erhs[vi*numdofpernode_+k] -= -fac_fns_i0_funct_vi*pow_conint_gamma_k*expterm;
    }

    break;
  }

  // linear law:  i_n = frt*i_0*((alphaa+alpha_c)*(V_M - phi)) -> changed 10/13
  // previously implemented: i_n = i_0*(alphaa*frt*(V_M - phi) + 1.0)
  //                         -> linearization in respect to anodic branch!!
  //                         this is not the classical verion of a linear electrode kinetics law
  case INPAR::SCATRA::linear:
  {
    // read model-specific parameter
    const double alphaa = cond->GetDouble("alpha");
    double       i0 = cond->GetDouble("i0");
    const double dlcap = cond->GetDouble("dl_spec_cap");
    double pot0hist = 0.0;
    if(dlcap!=0.0)
      pot0hist = cond->GetDouble("pot0hist");
    i0*=timefac;
    if (i0 < -EPS14) dserror("i0 is negative, \n"
                             "a positive definition is necessary due to generalized reaction models: %f",i0);
    const double gamma = cond->GetDouble("gamma");
    const double refcon = cond->GetDouble("refcon");
    if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

    if(valence_k!=nume)
      dserror("Kinetic model Butler-Volmer: The number of transferred electrons need to be  \n "
          "the same as the charge number of the reacting species %i", k);
    // opencircuit potential is assumed to be zero
    const double ocp = 0.0;
    // overpotential based on opencircuit potential
    const double eta = epd - ocp;

    double pow_conint_gamma_k = 0.0;
    if ((conint[k]/refcon) < EPS13)
    {
      pow_conint_gamma_k = std::pow(EPS13,gamma);
#ifdef DEBUG
      std::cout<<"WARNING: Rel. Conc. in Tafel formula is zero/negative: "<<(conint[k]/refcon)<<std::endl;
      std::cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<std::endl;
#endif
    }
    else
      pow_conint_gamma_k = std::pow(conint[k]/refcon,gamma);
    const double linearfunct = (alphaa*frt*eta);
    // note: gamma==0 deactivates concentration dependency
    double concterm = 0.0;
    if (conint[k] > EPS13)
      concterm = gamma*pow(conint[k],(gamma-1.0))/pow(refcon,gamma);
    else
      dserror("Better stop here!");

    for (int vi=0; vi<nen_; ++vi)
    {
      const double fac_fns_i0_funct_vi = scalar*fac*fns*i0*funct(vi);
      const int fvi = vi*numdofpernode_+k;
      // ---------------------matrix
      for (int ui=0; ui<nen_; ++ui)
      {
        emat(fvi,ui*numdofpernode_+k) += -fac_fns_i0_funct_vi*concterm*funct(ui)*linearfunct;
        emat(fvi,ui*numdofpernode_+numscal_) += -fac_fns_i0_funct_vi*pow_conint_gamma_k*(-alphaa)*frt*funct(ui);
      }
      // ------------right-hand-side
      erhs[fvi] -= -fac_fns_i0_funct_vi*pow_conint_gamma_k*linearfunct;
    }

    if(dlcap!=0.0)
    {
      for (int vi=0; vi<nen_; ++vi)
      {//add terms of double layer capacitance current density
        for (int ui=0; ui<nen_; ++ui)
          emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += fac*funct(vi)*funct(ui)*dlcap/(nume*faraday);

        // -----right-hand-side: -R_k = -(theta*dt*(-1)*(w_k,j_k)
        erhs[vi*numdofpernode_+k] += fac*funct(vi)*dlcap/(nume*faraday)*(phihistint-pot0hist-potint+pot0);
      }
    }

    break;
  }

  case INPAR::SCATRA::butler_volmer_newman:
  {
    // "Electrochemical systems"
    // Newman and Thomas-Alyea, 2004
    // General stoichiometry: pp. 212-213, e.q. 8.26
    // consideration of a elementary step of the form:
    // Sum_i s_i M_i ->  ne-
    // n is one if charge transfer is involved, multiple electron transfers "being unlikely in
    // an elementary step

    const double k_a = cond->GetDouble("k_a");
    const double k_c = cond->GetDouble("k_c");
    const double beta = cond->GetDouble("beta");
    const double dlcap = cond->GetDouble("dl_spec_cap");
    if(dlcap!=0.0) dserror("double layer charging is not implemented for Butler-Volmer-Newman electrode kinetics");

    //reaction order of the cathodic and anodic reactants of ionic species k
    std::vector<int> q(numscal_,0);
    std::vector<int> p(numscal_,0);

    for(int ii=0; ii<numscal_; ii++)
    {
      //according to the convention: anodic reactant is positiv
      if(stoich[ii] > 0)
      {
        q[ii] = 0;
        p[ii] = stoich[ii];
      }
      //according to the convention: cathodic reactant is negative
      else
      {
        q[ii]= -stoich[ii];
        p[ii] = 0;
      }
    }

#ifdef DEBUG
    // some safety checks/ user warnings
    if (((1-beta)*frt*epd) > 100.0)
      std::cout<<"WARNING: Exp((1-beta)...) in Butler-Volmer law is near overflow!"
      <<exp((1-beta)*frt*epd)<<std::endl;
    if (((-beta)*frt*epd) > 100.0)
      std::cout<<"WARNING: Exp(-beta...) in Butler-Volmer law is near overflow!"
      <<exp((-beta)*frt*epd)<<std::endl;
#endif

    double pow_conint_p(1.0);      //product over i (c_i)^(p_i)
    double pow_conint_q(1.0);      //product over i (c_i)^(q_i)
    std::vector<double> pow_conint_p_derivative(numscal_,1.0);  //pow_conint_p derivated after conint[nspec]
    std::vector<double> pow_conint_q_derivative(numscal_,1.0); //pow_conint_q derivated after conint[nspec]

    //concentration term (product of cathodic and anodic species)
    for(int kk=0; kk<numscal_; ++kk)
    {
      if ((conint[kk]) < EPS13) // 1.0E-16)
      {
        pow_conint_p *= std::pow(EPS13,p[kk]);
        pow_conint_q *= std::pow(EPS13,q[kk]);
#ifdef DEBUG
        std::cout<<"WARNING: Rel. Conc. of species" <<k<<" in Butler-Volmer formula is zero/negative: "<<(conint[k])<<std::endl;
        std::cout<<"-> Replacement value: pow(1.0E-16,p[ispec]) = "<< pow(EPS13,p[k]) << " pow(1.0E-13,q[k]) = "<< pow(EPS13,q[k]) <<std::endl;
#endif
      }
      else
      {
        pow_conint_p *= std::pow((conint[kk]),p[kk]);
        pow_conint_q *= std::pow((conint[kk]),q[kk]);
      }
    }

    //derivation of concentration term  with respect to ionic species kk
    for(int kk=0; kk<numscal_; ++kk)
    {
      pow_conint_p_derivative[kk] = pow_conint_p*p[kk]/conint[kk];
      pow_conint_q_derivative[kk] = pow_conint_q*q[kk]/conint[kk];
    }

    // loop over reacting species; determines the line of the matrix
    const double expterma = exp((1-beta)*nume*frt*epd);
    const double exptermc = exp((-beta)*nume*frt*epd);

    for (int vi=0; vi<nen_; ++vi)
    {
      // see Wittmann, Erweiterte Reaktionsmodelle für die numerische Simulation von
      // elektrochemischen Systemen, p.20, equ. 3.4
      const double fac_fns_funct_vi = scalar*faraday*nume*fac*fns*funct(vi);
      for (int ui=0; ui<nen_; ++ui)
      {
        //loop over the columns of the matrix, makes sure that the linearisation w.r.t the first concentration is added to the first column
        for(int kk=0; kk<numscal_; ++kk)
        {
          emat(vi*numdofpernode_+k,ui*numdofpernode_+kk) += -fac_fns_funct_vi*((k_a*expterma*pow_conint_p_derivative[kk]) - (k_c*exptermc*pow_conint_q_derivative[kk]))*funct(ui)*timefac;
        }
        //linearisation w.r.t the potential
        emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += -fac_fns_funct_vi*((-k_a*(1-beta)*nume*frt*expterma*pow_conint_p) - (k_c*beta*nume*frt*exptermc*pow_conint_q))*funct(ui)*timefac;
      }
      // ------------right-hand-side
      erhs[vi*numdofpernode_+k] -= -(fac_fns_funct_vi*((k_a*expterma*pow_conint_p)-(k_c*exptermc*pow_conint_q)))*timefac;
    }

    break;
  }

  case INPAR::SCATRA::butler_volmer_bard:
  {
    // "Electrochemical Methods Fundamentals and Applications"
    // Bard and Faulkner, 2001, pp. 94 ff; pp. 99 eq. 3.4.10
    // reaction model for a one-step, one-electron process (elementar step)
    // O + e -> R
    const double e0 = cond->GetDouble("e0");
    const double k0 = cond->GetDouble("k0");
    const double beta = cond->GetDouble("beta");
    const double c_c0 = cond->GetDouble("c_c0");
    const double c_a0 = cond->GetDouble("c_a0");
    const double dlcap = cond->GetDouble("dl_spec_cap");
    if(dlcap!=0.0) dserror("double layer charging is not implemented for Butler-Volmer-Bard electrode kinetics");

    if(nume!=1)
      dserror("electron != 1; \n "
          "this Butler-Volmer-equation (Bard/Faulkner) works for elementary steps (one electron) only!");

    // only one reactant and product are supported by the basic model
    // only stoichiometry of 1
    {
      int check1 = 0;
      int check2 = 0;
      for(int kk=0; kk<numscal_;kk++)
      {
        if(abs(stoich[kk])>1)
          dserror("Stoichiometry is larger than 1!! \n"
                  "This is not supported by the reaction model based on Bard");

        check1 += abs(stoich[kk]);
        check2 += stoich[kk];
      }
      if (check1>2 or check1==0)
        dserror("More than one reactant or product defined!! \n"
                "This is not supported by the reaction model based on Bard");

      // At the moment it is not checked if two products or reactants are defined
    }

    // equilibrium potential (equilpot):
    // defined in Bard, 2001, p.98, eq. 3.4.3
    const double equilpot = e0 + (log(c_c0/c_a0))/(frt*nume);
    // overpotential based on equilibrium potential
    const double eta_equilpot = epd - equilpot;

    // negative sign: we look at electon flow
    const double i0 = k0*pow(c_c0,1-beta)*pow(c_a0,beta)*nume*faraday;

    // reactant or product not a species in the electrolyte
    // -> concentration = 1.0
    double conctermc = 1.0;
    double concterma = 1.0;
    double conctermc_der = 1.0;
    double concterma_der = 1.0;
    //species id of the anodic and cathodic reactant
    int anodic = 0;
    int cathodic = 0;
    bool checkc = 0;
    bool checka = 0;

    // concentration terms for anodic and cathodic reaction
    // only one reactant and product are supported by the basic model
    // only stoichiometry of 1
    for(int kk=0; kk<numscal_;kk++)
    {
      if(stoich[kk]==1)
      {
        concterma = conint[kk]/c_a0;
        concterma_der = 1.0/c_a0;
        anodic = kk;
        checka = true;
      }
      else if(stoich[kk]==-1)
      {
        conctermc = conint[kk]/c_c0;
        conctermc_der = 1.0/c_c0;
        cathodic = kk;
        checkc = true;
      }
    }

#ifdef DEBUG
    // some safety checks/ user warnings
    if (((1-beta)*(frt*nume)*eta_equilpot) > 100.0)
      std::cout<<"WARNING: Exp((1-beta)...) in Butler-Volmer law is near overflow!"
      <<exp((1-beta)*(frt*nume)*eta_equilpot)<<std::endl;
    if (((-beta)*(frt*nume)*eta_equilpot) > 100.0)
      std::cout<<"WARNING: Exp(-beta...) in Butler-Volmer law is near overflow!"
      <<exp((-beta)*(frt*nume)*eta_equilpot)<<std::endl;
#endif

    const double expterma = exp((1-beta) * (frt*nume) * eta_equilpot);
    const double exptermc = exp((-beta) * (frt*nume) * eta_equilpot);

    for (int vi=0; vi<nen_; ++vi)
    {
      const double fac_i0_funct_vi = scalar*fac*fns*i0*funct(vi);
      // ---------------------matrix
      for (int ui=0; ui<nen_; ++ui)
      {
        //derivation wrt concentration
        if(checkc == true)
          emat(vi*numdofpernode_+k,ui*numdofpernode_+cathodic) += fac_i0_funct_vi*conctermc_der*exptermc*funct(ui)*timefac;
        if(checka == true)
          emat(vi*numdofpernode_+k,ui*numdofpernode_+anodic)   += -fac_i0_funct_vi*concterma_der*expterma*funct(ui)*timefac;
        //derivation wrt potential
        emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_)
          += -fac_i0_funct_vi*(-concterma*(1-beta)*frt*nume*expterma-conctermc*beta*nume*frt*exptermc)*funct(ui)*timefac;
      }
      // ------------right-hand-side
      erhs[vi*numdofpernode_+k] -= -fac_i0_funct_vi*(concterma*expterma - conctermc*exptermc)*timefac;
    }

    break;
  }

  case INPAR::SCATRA::nernst:
    break;

  default:
  {
    dserror("Kinetic model not implemented!");
    break;
  }
  }

  return;
} // DRT::ELEMENTS::ScaTraEleUtilsElch<distype>::EvaluateElchKineticsAtIntegrationPoint


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleUtilsElch<DRT::Element::nurbs27>;
