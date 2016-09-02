/*----------------------------------------------------------------------*/
/*!
\file drt_function.H

\brief Managing and evaluating of spatial functions

<pre>
\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_FUNCTION_H
#define DRT_FUNCTION_H


#include <iostream>
#include <vector>
#include <complex>
#include <Teuchos_RCP.hpp>

#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_timecurve.H"
#include "drt_dserror.H"
#include "drt_parser.H"
//#include "drt_inputreader.H"

namespace MAT {
class Material;
}

namespace Sacado
{
  namespace Fad
  {
    template <typename> class DFad;
  }
}

namespace DRT {

  class Discretization;
  class Condition;
  class Problem;


namespace INPUT {

  class Lines;
  class DatFileReader;

}

namespace UTILS {


  /*!

  \brief abstract interface for spatial functions.

  Spatial functions can be defined either scalar-valued or vector-valued.

  There exist function definitions using explicitly defined spatial functions
  for test problems (beltrami, kim-moin etc).

  Spatial function based on parsed expressions are provided as well.

  */
  class Function
  {
  public:
    virtual ~Function() {}

    /*!

    \brief evaluate function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated

    */
    virtual double Evaluate(int index, const double* x, double t, DRT::Discretization* dis) = 0;

    /*!

    \brief evaluate derivatives function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated
    \param t     (i) Absolute time in which the
                     function will be evaluated
    \param dis   (i) discretization

    */
    virtual std::vector<std::vector<double> > FctDer(int index, const double* x, const double t, DRT::Discretization* dis)
      {
        return std::vector<std::vector<double> >(0);
      };

    /*!

    \brief Return the number of components of this spatial function
    (1 for scalar functions, dim for vector-valued functions)

    \return number of components

    */
    virtual int NumberComponents()
      {
        return(-1);
      };

  };


  /// class that manages all defined spatial functions
  class FunctionManager
  {
  public:

    FunctionManager() {}

    Teuchos::RCP<DRT::INPUT::Lines> ValidFunctionLines();

    /// read the baci input file and set up all functions
    void ReadInput(DRT::INPUT::DatFileReader& reader);

    /// return a specified function
    Function& Funct(int num);

  private:

    /// no copying
    FunctionManager(const FunctionManager&);

    /// no copying
    FunctionManager& operator=(const FunctionManager&);

    /// all functions defined. The functions may be vector valued!
    std::vector<Teuchos::RCP<Function> > functions_;
  };


  /// special implementation for beltrami flow (velocity, pressure)
  class BeltramiUP : public Function
  {
  public:

    BeltramiUP(int mat_id);

    BeltramiUP( Teuchos::RCP<MAT::Material> & mat );


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    virtual int NumberComponents()
    {
      return(4);
    };

  private:
    double density_;
    double kinviscosity_;
  };



  /// special implementation beltrami flow (gradient of velocity)
  class BeltramiGradU : public Function
  {
  public:

    BeltramiGradU(int mat_id );

    BeltramiGradU( Teuchos::RCP<MAT::Material> & mat );


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents()
    {
      return(9);
    };

  private:
    double kinviscosity_;
  };



  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (velocity, pressure)
  class KimMoinUP : public Function
  {
  public:

    KimMoinUP(int mat_id, bool is_stationary);

    KimMoinUP( Teuchos::RCP<MAT::Material> & mat, bool is_stationary );


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,v,w,p)

     */
    virtual int NumberComponents()
    {
      return(4);
    };

  private:
    double density_;
    double kinviscosity_;
    bool   is_stationary_;
  };



  /// special implementation for 2d (implemented for 3D) stationary kim-moin flow (gradient of velocity)
  class KimMoinGradU : public Function
  {
  public:

    KimMoinGradU(int mat_id, bool is_stationary );

    KimMoinGradU( Teuchos::RCP<MAT::Material> & mat, bool is_stationary );


    /*!

      \brief evaluate function at given position in space

      \param index (i) index defines the function-component which will
                       be evaluated
      \param x     (i) The point in space in which the function will be
                       evaluated

     */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

      \brief Return the number of components of this spatial function
      (This is a vector-valued function)

      \return number of components (u,x , u,y , u,z , v,x , v,y , v,z , w,x , w,y , w,z )

     */
    virtual int NumberComponents()
    {
      return(9);
    };

  private:
    double kinviscosity_;
    bool   is_stationary_;
  };

  /// spatial function based on parsed expression
  class ExprFunction : public Function
  {
  public:
    /*!

    \brief construct spatial function from expression with given origin

    \note  Upon construction, the object defines a spatial function
           returning the same function value for every spatial dimension.
           If a vector-valued spatial function is required, further
           expressions can be added via the AddExpr function. In this
           case, the spatial function defined upon construction will be
           the first component of the vector-valued spatial function.

    \param buf (i) (c-string) expression to be parsed during evaluation
    \param x   (i) x-coordinate of the origin of the coordinate system
    \param y   (i) y-coordinate of the origin of the coordinate system
    \param z   (i) z-coordinate of the origin of the coordinate system

    */
    ExprFunction(char* buf, double x, double y, double z);

    /*!

    \brief Default constructor creating empty object. Expressions are
           added with add function

    */
    ExprFunction();

    /*!

    \brief clean up parse tree

    */
    virtual ~ExprFunction();


    /*!

    \brief evaluate function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated

    */
    virtual double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief evaluate derivatives function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated
    \param t     (i) Absolute time in which the
                     function will be evaluated
    \param dis   (i) discretization

    */
    virtual std::vector<std::vector<double> > FctDer(int index, const double* x, const double t, DRT::Discretization* dis);

    /*!

    \brief add expression to an existing ExprFunction in order to extend
           it to a vector-valued spatial function.

           Every call to AddExpr adds one more component to the
           vector-valued function.

    \param buf (i) (c-string) expression to be parsed during evaluation of this component
    \param x   (i) x-coordinate of the origin of the coordinate system of this component
    \param y   (i) y-coordinate of the origin of the coordinate system of this component
    \param z   (i) z-coordinate of the origin of the coordinate system of this component

    */
    virtual void AddExpr(std::string buf, double x, double y, double z);
    /*!

    \brief Return the number of components of this spatial function
    (1 for scalar functions, dim for vector-valued functions)

    \return number of components

    */
    virtual int NumberComponents()
    {
      return(expr_.size());
    };

  protected:

    /*
      for scalar spatial functions returning the same value for all
      dimensions:

     -----------------------------------------+-------+-------+-------
      spatial function for dimension          |   0   |   1   |   2
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component |         x_[0]
     -----------------------------------------+-----------------------
      origin of spatial function, y-component |         y_[0]
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component |         z_[0]
     -----------------------------------------+-------+-------+-------


      for vector-valued spatial functions, returning separate values
      for all dimensions:

     -----------------------------------------+-------+-------+-------
      spatial function for dimension          |   0   |   1   |   2
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component | x_[0] | x_[1] | x_[2]
     -----------------------------------------+-------+-------+-------
      origin of spatial function, y-component | y_[0] | y_[1] | y_[2]
     -----------------------------------------+-------+-------+-------
      origin of spatial function, x-component | z_[0] | z_[1] | z_[2]
     -----------------------------------------+-------+-------+-------

    */

    int dim_; //! problem dimension determining evaluated components

    std::vector<double> x_; //! origin(s) of spatial function, x-component
    std::vector<double> y_; //! origin(s) of spatial function, y-component
    std::vector<double> z_; //! origin(s) of spatial function, z-component

    std::vector<Teuchos::RCP<DRT::PARSER::Parser<double> > > expr_; //! expression syntax tree(s)
    std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > > > > exprd_;
  };

  /// spatial function based on parsed expression
  class VariableExprFunction : public ExprFunction
  {
  public:
    /*!

    \brief Default constructor creating empty object. Expressions are
           added with add function

    */
    VariableExprFunction();

    /*!

    \brief clean up parse tree

    */
    virtual ~VariableExprFunction();

    /*!

    \brief add expression to an existing ExprFunction in order to extend
           it to a vector-valued spatial function.

           Every call to AddExpr adds one more component to the
           vector-valued function.

    \param buf (i) (c-string) expression to be parsed during evaluation of this component
    \param x   (i) x-coordinate of the origin of the coordinate system of this component
    \param y   (i) y-coordinate of the origin of the coordinate system of this component
    \param z   (i) z-coordinate of the origin of the coordinate system of this component

    */
    void AddExpr(std::string buf, double x, double y, double z);

    /*!

    \brief evaluate function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable

    */
    double Evaluate(int index, const std::vector<std::pair<std::string,double> >& variables);

    /*!

    \brief evaluate derivatives function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated
    \param t     (i) Absolute time in which the
                     function will be evaluated
    \param dis   (i) discretization

    */
    std::vector<std::vector<double> > FctDer(int index, const double* x, const double t, DRT::Discretization* dis);

    /*!

    \brief  at given position in space for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable

    */
    std::vector<std::vector<double> > FctDer(int index, const std::vector<std::pair<std::string,double> >& variables);

    /*!
    \brief  check if a given name is a variable

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param varname  (i)  Variable name to be checked

    */
    bool IsVariable(int index, const std::string& varname) const;

    /*!
    \brief  add a variable to the parser object

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param varname  (i)  Variable name

    \param varvalue (i)  Value of the variable

    */
    void AddVariable(int index, const std::string& varname, double varvalue);

    /*!
    \brief  parse the underlying expression
    */
    void ParseExpressions();

  private:

    bool isparsed_; //! flag indicating whether expression was parsed
  };

  /// special implementation for 3d Beltrami flow
  class BeltramiFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents()
    {
      return(4);
    };

  };


  /// special implementation for 2d Kim-Moin flow
  class KimMoinFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (this is a vector-valued functions)

    \return number of components (u,v,p)

    */
    virtual int NumberComponents()
    {
      return(3);
    };

  };


  /// special implementation for 2d Bochev test case (velocity and pressure)
  class BochevUPFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,p)

    */
    virtual int NumberComponents()
    {
      return(3);
    };

  };


  /// special implementation for 2d Bochev test case (rhs function)
  class BochevRHSFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (f1,f2)

    */
    virtual int NumberComponents()
    {
      return(2);
    };

  };



  /// special implementation for beltrami flow (rhs)
  class BeltramiRHS : public Function
  {
  public:


    BeltramiRHS(int mat_id, bool is_stokes);

    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w)

    */
    virtual int NumberComponents()
    {
      return(3);
    };

  private:
    double kinviscosity_;
    bool is_stokes_;

  };



  /// special implementation for 2d(3D) stationary kim-moin flow (rhs) for pure stokes equation
  class KimMoinRHS : public Function
  {
  public:


    KimMoinRHS(int mat_id, bool is_stationary, bool is_stokes);

    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w)

    */
    virtual int NumberComponents()
    {
      return(3);
    };

  private:
    double kinviscosity_;
    bool is_stationary_;
    bool is_stokes_;

  };


  /// special implementation for (randomly) disturbed 3d turbulent
  /// boundary-layer profile
  /// (currently fixed for low-Mach-number flow through a backward-facing step,
  ///  but may be easily manipulated to fit other profiles in other geometries)
  class TurbBouLayerFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents()
    {
      return(4);
    };

  };


  /// special implementation for (randomly) disturbed 3d turbulent
  /// boundary-layer profile
  /// (incompressible flow over backward-facing step,
  ///  corresponding to geometry of DNS by Le, Moin and Kim)
  class TurbBouLayerFunctionBFS : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents()
    {
      return(4);
    };

  };


  /// special implementation for (randomly) disturbed 3d turbulent boundary-layer profile
  /// (incompressible flow in the ORACLES test rig)
  class TurbBouLayerFunctionORACLES : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

  };


  /// special implementation for Womersley blood flow
  class WomersleyFunction : public Function
  {
  public:

    // ctor
    WomersleyFunction(bool locsys, int e, int mat, int curve, bool fsi);


    // evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);
    // Bessel Functions of the first kind and order 0 or 1 for a complex argument
    std::complex<double> BesselJ01(std::complex<double> z, bool order);
    // perform a discrete fourier transformation on given array
    void DFT(std::vector<double> *data, std::vector< std::complex<double> > *resdata, const int N);

  private:

    bool                   isinit_;
    // switch for use of local coordinate systems
    bool                   locsys_;
    int                    locsysid_;
    // current edge node radius
    double                  radius_;
    // a time curve
    DRT::UTILS::TimeCurve& tc_;
    // number of the material in the input file
    int                    mat_;
    double                 viscosity_;
    // FSI switch
    bool                   fsi_;
    // toggle coordinate transformation of edge node (once per time step)
    bool                    dotrafo_;
    // store t_(n-1) for comparison with t_n
    double                 tnminus1_;

    // further variables
    // number of harmonics that are used in the synthesis of the timecurve
    int                     noharm_;
    // time curve frequency
    double                 fbase_;
           // time curve value of the previous time step (needed in current version to circumvent division by 0)
          double                 tcprevious_;
    // imaginary number i
    std::complex<double>   i_;
    // storage vector for velocity@1s for profile transition 0<t<1
    std::vector<double>     vtemp_;
    // storage vector for Fourier Transform output
    std::vector< std::complex<double> > fouphyscurve_;

    // exist after init phase if locsys_==true
    // to do: move this stuff to separate class in locsys.H/.cpp
    // with functionality to transform spatial vectors from/to local system
    Condition*             locsyscond_;
    std::vector<double>    normal_;
    std::vector<double>    tangent1_;
    std::vector<double>    tangent2_;
    std::vector<double>    origin_;

    // exist after init phase if dirich_==true
    // for edge nodes (polar coordinates)
    // location of smallest modulus of phi, phi<0
    int                    iminminus_;
    // location of smallest modulus of phi, phi>0
    int                    iminplus_;
    // location of largest modulus of phi, phi<0
    int                    imaxminus_;
    // location of largest modulus of phi, phi>0
    int                    imaxplus_;
    // vector with inflow surface ids
    std::vector<int>       surfnodeids_;
    // vector with edge node ids
    std::vector<int>       nodeids_;
    // phase vector
    std::vector<double>    phi_;
    // distances between center of cross section and nodes
    std::vector<double>    noderadius_;
  };


  /// special implementation for stationary 2d Jeffery-Hamel flow
  class JefferyHamelFlowFunction : public Function
  {
  public:
    /*!

    \brief evaluate function at given position in space

    \param index (i) index defines the function-component which will
                     be evaluated
    \param x     (i) The point in space in which the function will be
                     evaluated

    */
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    /*!

    \brief Return the number of components of this spatial function
    (This is a vector-valued function)

    \return number of components (u,v,w,p)

    */
    virtual int NumberComponents()
    {
      return(4);
    };

    /*!

    compute radial flow u as a function of the angle alpha

    \return radial velocity u

    \note this function is static such that we can get the radial velocity
          without creating the Function object

    */
    static double RadialVelocity(
      const double& theta ///< angle between 0 and PI/4 (range is checked in debug mode)
      );

  };


  /// special implementation for controlled rotations
  class ControlledRotationFunction : public Function
  {
  public:

    /// ctor
    ControlledRotationFunction(std::string fileName, std::string type, double origin_x, double origin_y, double origin_z);

    /// evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

  private:
    // Condition type: STRUCTURE=1, FLUID=2
    int type_;

    // Origin, about which the rotation shall be performed
    LINALG::Matrix<3,1> origin_;

    // Time of previous time step (at t-deltaT)
    double timeOld_;

    // Number of maneuver cells (variables)
    const int NUMMANEUVERCELLS_;

    // Number of maneuvers
    int numManeuvers_;

    // Double Vector containing maneuver information (t, omegaDot_x_B, omegaDot_y_B, omegaDot_z_B)
    std::vector<double> maneuvers_;

    // Previous angular acceleration (at t-deltaT)
    LINALG::Matrix<3,1> omegaDotOld_B_;

    // Current angular rate (at t)
    LINALG::Matrix<3,1> omega_B_;

    // Satellite's current attitude trafo matrix from B- to I-system (at t)
    LINALG::Matrix<3,3> satAtt_dcm_IB_;

    // Satellite's current attitude quaternion from B- to I-system (at t)
    LINALG::Matrix<4,1> satAtt_q_IB_;
  };

  /// special implementation for acceleration profiles
  class AccelerationProfileFunction : public Function
  {
  public:

    /// ctor
    AccelerationProfileFunction(std::string fileName);

    /// evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

  private:
    // Time of previous time step (at t-deltaT)
    double timeOld_;

    // Number of cells (variables) (time + 3-dim acc = 4)
    const int NUMACCELERATIONCELLS_;

    // Number of acceleration rows
    int numAccelerations_;

    // Double Vector containing acceleration information (t, acc_x_B, acc_y_B, acc_z_B)
    std::vector<double> accelerations_;

    // Current acceleration (at t)
    LINALG::Matrix<3,1> acc_B_;
  };

  /// special implementation for ramping to a specified value
  class RampToValueFunction : public Function
  {
  public:

    /// ctor
    RampToValueFunction(double value, double startTime, double duration, std::string type);

    /// evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

  private:
    // Time of previous time step (at t-deltaT)
    double timeOld_;

    // Value to which is ramped
    double value_;

    // Time to start the ramping
    double startTime_;

    // Duration of the ramping
    double duration_;

    // Ramping type (1: ZERORAMPCONST, 2: ZERORAMPZERO)
    int type_;

    // Current value (at t)
    double valueCurrent_;
  };

  /// special implementation for the node normals of a specified geometry
  class NodeNormalFunction : public Function
  {
  public:
    /// ctor
    NodeNormalFunction(std::string type, std::vector<double>* origin, double radius, double cylinderHeight, std::vector<double>* orientation, double CassiniA);

    /// Evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

    // Determine node normal of a sphere
    double NodeNormalSphere(int index, const double* xp, std::vector<double> origin);

    // Determine node normal of a Cassini volume
    double NodeNormalCassini(int index, const double* xp, std::vector<double> origin);

  private:
    // Problem dimension (2D or 3D, i.e. 2 or 3)
    int dim_;

    // Geometry type
    int type_;

    // Origin of the geometry
    std::vector<double> origin_;

    // Radius
    double radius_;

    // Cylinder height, divided by 2
    double cylinderHeightHalf_;

    // Orientation of the geometry (symmetry axis)
    std::vector<double> orientation_;

    // Cassini value a
    double CassiniA_;

    // New reference origin, if nodes lie in the positive or negative spherical part
    std::vector<double> originRefSpherePos_;
    std::vector<double> originRefSphereNeg_;
  };

  /// special implementation for the rotation vector used in locsys conditions
  class RotationVectorForNormalSystemFunction : public Function
  {
  public:

    /// ctor
    RotationVectorForNormalSystemFunction(int geoFunct);

    /// evaluate function at given position in space
    double Evaluate(int index, const double* x, double t, DRT::Discretization* dis);

  private:
    // Problem dimension (2D or 3D, i.e. 2 or 3)
    int dim_;

    // Function that calculates the node normal for the specified geometry
    int geoFunct_;
  };

} // namespace UTILS

} // namespace DRT


void PrintFunctionDatHeader();


#endif
