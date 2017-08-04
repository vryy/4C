/*----------------------------------------------------------------------*/
/*!
\file scatra_multiscale.H

\brief auxiliary material for macro-scale elements in multi-scale simulations of scalar transport problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_MULTISCALE_H
#define SCATRA_MULTISCALE_H

#include <Teuchos_RCP.hpp>

namespace MAT
{
  // forward declaration
  class ScatraMultiScaleGP;

  namespace PAR
  {
    class Material;

    //! material parameters
    class ScatraMultiScale
    {
      public:

        //! constructor
        ScatraMultiScale(
            Teuchos::RCP<MAT::PAR::Material> matdata
            );

        //! destructor
        virtual ~ScatraMultiScale(){return;};

        //! return name of micro-scale input file
        const std::string MicroInputFileName() const {return microfile_;};

        //! return number of micro-scale discretization
        int MicroDisNum() const {return microdisnum_;};

        //! return specific micro-scale surface area
        double A_s() const {return A_s_;};

      protected:

        //! @name material parameters
        //@{
        //! name of micro-scale input file
        const std::string microfile_;

        //! number of micro-scale discretization
        const int microdisnum_;

        //! specific micro-scale surface area
        const double A_s_;
        //@}
    };  // class MAT::PAR::ScatraMultiScale
  } // namespace PAR

  /*----------------------------------------------------------------------*/
  //! material wrapper
  class ScatraMultiScale
  {
    public:

      //! destructor
      virtual ~ScatraMultiScale(){return;};

      //! initialize multi-scale scalar transport material
      void Initialize(
          const int   ele_id,   //!< macro-scale element ID
          const int   gp_id     //!< macro-scale Gauss point ID
          );

      //! prepare time step on micro scale
      void PrepareTimeStep(
          const int                    gp_id,        //!< macro-scale Gauss point ID
          const std::vector<double>&   phinp_macro   //!< macro-scale state variables
          ) const;

      //! evaluate multi-scale scalar transport material
      void Evaluate(
          const int                    gp_id,           //!< macro-scale Gauss point ID
          const std::vector<double>&   phinp_macro,     //!< macro-scale state variables
          double&                      q_micro,         //!< micro-scale flux
          std::vector<double>&         dq_dphi_micro,   //!< derivatives of micro-scale flux w.r.t. macro-scale state variables
          const bool                   solve = true     //!< flag indicating whether micro-scale problem should be solved
          ) const;

      //! evaluate mean concentration on micro scale
      double EvaluateMeanConcentration(
          const int   gp_id   //!< macro-scale Gauss point ID
          ) const;

      //! update multi-scale scalar transport material
      void Update(
          const int   gp_id   //!< macro-scale Gauss point ID
          ) const;

      //! create output on micro scale
      void Output(
          const int   gp_id   //!< macro-scale Gauss point ID
          ) const;

      //! read restart on micro scale
      void ReadRestart(
          const int   gp_id   //!< macro-scale Gauss point ID
          ) const;

      //! return name of micro-scale input file
      const std::string MicroInputFileName() const {return Params()->MicroInputFileName();};

      //! return number of micro-scale discretization
      int MicroDisNum() const {return Params()->MicroDisNum();};

      //! return specific micro-scale surface area
      double A_s() const {return Params()->A_s();};

    protected:

      //! construct empty material
      ScatraMultiScale();

    private:

      //! material parameters
      virtual const MAT::PAR::ScatraMultiScale* Params() const = 0;

      //! map between Gauss point ID and Gauss point submaterial
      std::map<int,Teuchos::RCP<ScatraMultiScaleGP> > matgp_;
  }; // material wrapper
} // namespace MAT
#endif
