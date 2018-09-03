/*!----------------------------------------------------------------------
\file acou_pml.cpp
\brief PML definition for HDG acoustics

<pre>
\level 3

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
 *----------------------------------------------------------------------*/

#include "acou_pml.H"

#ifdef HAVE_DEAL_II


#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/lac/parallel_vector.h>

#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

#include <fstream>
#include <iostream>
#include <iomanip>


namespace ACOU
{
  using namespace dealii;


  /*
   * This function evaluates the scalar profil-function  of the pml definition
   * depending on the distance from the interface and the width of the layer.
   * For distance < 0 the return value is in general zero.
   *
   * Input:
   *  distance:   the distance xi = n(x-x_0)
   *  width:      the thickness of the layer
   *  fct_choice: wich function should be used
   *              0 = constant
   *              1 = linear
   *              2 = quadratic
   *              3 = cubic
   *              4 = quartic
   *              5 = hyperbolic
   *              6 = shifted hyperbolic
   *
   * Output:
   *  value sigma of the attenuation function
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  inline double InterfaceDataPML<dim, Number>::sigma_fct(
      double distance, double width, int fct_choice) const
  {
    double return_value = 0.0;
    double argument = width > 1e-3 ? distance / width : 1e3 * distance;
    if (distance >= 0) switch (fct_choice)
      {
        case 0:
          return_value = 1.0;
          break;
        case 1:
          return_value = argument;
          break;
        case 2:
          return_value = argument * argument;
          break;
        case 3:
          return_value = argument * argument * argument;
          break;
        case 4:
          return_value = argument * argument * argument * argument;
          break;
        case 5:
          return_value = 1.0 / (kappa - argument);
          break;
        case 6:
          return_value = (1.0 / (kappa - argument)) - (1.0 / kappa);
          break;
        default:
          Assert(false, ExcMessage("attenuation function not valid"));
      }
    else
      return_value = 0.0;

    // scale with the attenuation value
    return_value *= attenuation;

    return return_value;
  }


  /*
   * This function evaluates integral of the scalar profil-function  of the pml definition
   * depending on the distance from the interface and the width of the layer.
   * For distance < 0 the return value is in general zero.
   *
   * Input:
   *  distance:   the distance xi = n(x-x_0)
   *  width:      the thickness of the layer
   *  fct_choice: wich function should be used
   *              0 = constant
   *              1 = linear
   *              2 = quadratic
   *              3 = cubic
   *              4 = quartic
   *              5 = hyperbolic
   *              6 = shifted hyperbolic
   *
   * Output:
   *  integral value of the attenuation function
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  inline double InterfaceDataPML<dim, Number>::sigma_fct_integral(
      double distance, double width, int fct_choice) const
  {
    double return_value = 0.0;
    double argument = width > 1e-3 ? distance / width : 1e3 * distance;
    if (distance >= 0) switch (fct_choice)
      {
        case 0:
          return_value = argument * width;
          break;
        case 1:
          return_value = 0.50 * width * argument * argument;
          break;
        case 2:
          return_value = (width / 3.0) * argument * argument * argument;
          break;
        case 3:
          return_value = 0.25 * width * argument * argument * argument * argument;
          break;
        case 4:
          return_value = 0.20 * width * argument * argument * argument * argument * argument;
          break;
        case 5:
          return_value = width * std::log(kappa / (kappa - argument));
          break;
        case 6:
          return_value = width * (std::log(kappa / (kappa - argument)) - (distance / kappa));
          break;
        default:
          Assert(false, ExcMessage("attenuation function not valid"));
      }
    else
      return_value = 0.0;

    // scale with the attenuation value
    return_value *= attenuation;

    return return_value;
  }

  /*
   * The standard constructor of the InterfaceDataPML class initialize all
   * variables with standard values.
   * The value of kappa will be set to kappa=1.25.
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  InterfaceDataPML<dim, Number>::InterfaceDataPML()
      : width(0), R_0(0), kappa(1.25), is_used(false), geometry(hyper_cube), attenuation(0.0)
  {
    normal_vector = 0;
    x_0 = Point<dim>();
  }


  /*
   * This constructor of the InterfaceDataPML class initialize all
   * variables with standard values, exepting the center (\vec{x_0}) and the radius.
   * This variables will be set with the input values of the
   * constructor.
   *
   * Input:
   *  center: \vec{x_0} of the pml definition
   *  radius: the radius of a circular or sperical pml
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  InterfaceDataPML<dim, Number>::InterfaceDataPML(const Point<dim> center, const double radius)
      : width(0), R_0(radius), kappa(1.25), is_used(false), geometry(hyper_cube), attenuation(0.0)
  {
    normal_vector = 0;
    x_0 = center;
  }

  /*
   * This function returns the normal vector of the layer and
   * works only for rectilinear layer properly.
   *
   * Output:
   *  normal vector \vec{n} of the layer
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  Tensor<1, dim, Number> InterfaceDataPML<dim, Number>::get_normal_vector() const
  {
    return normal_vector;
  }

  /*
   * This function returns the normal vector of the layer
   * depending on the coordinates of the quadrature point and
   * works for rectilinear and circular layer,
   *
   * Input:
   *  point: the coordinates of the quadrature point
   *
   * Output:
   *  normal vector \vec{n} of the layer
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  Tensor<1, dim, Number> InterfaceDataPML<dim, Number>::get_normal_vector(
      const Point<dim> point) const
  {
    Tensor<1, dim, value_type> return_vector;
    return_vector = 0;
    switch (geometry)
    {
      case hyper_cube:
      {
        return_vector = normal_vector;
        break;
      }
      case hyper_sphere_radial:
      case hyper_sphere_tangential:
      case hyper_sphere_radial_inv:
      case hyper_sphere_tangential_inv:
      {
        Tensor<1, dim, value_type> vec_diff;
        for (unsigned int i = 0; i < dim; ++i) vec_diff[i] = point(i) - x_0(i);

        double norm = vec_diff.norm();
        if (norm > 1e-14) return_vector = vec_diff / norm;
        break;
      }
      case hyper_cylinder_radial:
      case hyper_cylinder_tangential:
      {
        Tensor<1, dim, value_type> vec_diff;
        for (unsigned int i = 0; i < 2; ++i) vec_diff[i] = point(i) - x_0(i);

        double norm = vec_diff.norm();
        if (norm > 1e-14) return_vector = vec_diff / norm;
        break;
      }
      default:
        Assert(false, ExcNotImplemented());
    }

    return return_vector;
  }


  /*
   * This function returns the distance of the input-point
   * according to the formula xi = n(x-x_0).
   *
   * Input:
   *  point: the coordinates of the quadrature point
   *
   * Output:
   *  distance \xi = \vec{n}(\vec{x}-\vec{x_0})
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  double InterfaceDataPML<dim, Number>::get_distance(const Point<dim> point) const
  {
    double tmp = 0.0;
    Tensor<1, dim, value_type> vec_diff;

    for (unsigned int i = 0; i < dim; ++i) vec_diff[i] = point(i) - x_0(i);

    switch (geometry)
    {
      case hyper_cube:
      {
        tmp = vec_diff * normal_vector;
        break;
      }
      case hyper_sphere_radial:
      case hyper_sphere_tangential:
      case hyper_sphere_radial_inv:
      case hyper_sphere_tangential_inv:
      {
        tmp = vec_diff * get_normal_vector(point);
        break;
      }
      case hyper_cylinder_radial:
      case hyper_cylinder_tangential:
      {
        Tensor<1, dim, value_type> vec_tmp = get_normal_vector(point);
        tmp = 0.0;
        for (unsigned int i = 0; i < 2; ++i) tmp += vec_diff[i] * vec_tmp[i];
        break;
      }
      default:
        Assert(false, ExcNotImplemented());
    }

    return (tmp < 0 ? 0.0 : tmp);
  }



  /*
   * This function returns the sigma value at the quadrature point.
   *
   * Input:
   *  point: the coordinates of the quadrature point
   *
   * Output:
   *  sigma: attenuation value
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  double InterfaceDataPML<dim, Number>::get_sigma(const Point<dim> point) const
  {
    double distance = get_distance(point);
    double sigma = 0.0;
    if (distance > 1e-14)
    {
      switch (geometry)
      {
        case hyper_cube:
        {
          sigma = sigma_fct(distance, width, attenuation_fct);
          break;
        }
        case hyper_sphere_radial:
        case hyper_cylinder_radial:
        {
          sigma = sigma_fct(distance - R_0, width, attenuation_fct);
          break;
        }
        case hyper_sphere_tangential:
        case hyper_cylinder_tangential:
        {
          sigma = sigma_fct_integral(distance - R_0, width, attenuation_fct) / distance;
          break;
        }
        case hyper_sphere_radial_inv:
        {
          if (distance > R_0)
            sigma = 0.0;
          else if (distance > 1e-2)
          {
            sigma = sigma_fct(R_0 - distance, width, attenuation_fct);
          }
          else
            sigma = 0.0;
          break;
        }
        case hyper_sphere_tangential_inv:
        {
          if (distance > R_0)
            sigma = 0.0;
          else if (distance > 1e-2)
          {
            sigma = sigma_fct_integral(R_0, width, attenuation_fct) / distance;
            sigma -= sigma_fct_integral(R_0 - distance, width, attenuation_fct) / distance;
          }
          else
            sigma = 0.0;
          break;
        }
        default:
        {
          Assert(false, ExcNotImplemented());
        }
      }
    }
    return sigma;
  }


  /*
   * This function returns the tensor \vec{n}\otimes\vec{n} of a
   * rectilinear layer
   *
   * Input:
   *  point: the coordinates of the quadrature point
   *
   * Output:
   *  tensor: \vec{n}\otimes\vec{n}
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  Tensor<2, dim, Number> InterfaceDataPML<dim, Number>::get_tensor(const Point<dim> point) const
  {
    Tensor<2, dim, value_type> tensor;
    Tensor<1, dim, value_type> normal_vector;

    tensor = 0;
    normal_vector = get_normal_vector(point);
    switch (geometry)
    {
      case hyper_cube:
      case hyper_sphere_radial:
      case hyper_cylinder_radial:
      case hyper_sphere_radial_inv:
      {
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j) tensor[i][j] = normal_vector[i] * normal_vector[j];
        break;
      }
      case hyper_sphere_tangential:
      case hyper_cylinder_tangential:
      case hyper_sphere_tangential_inv:
      {
        for (unsigned int i = 0; i < dim; ++i)
        {
          tensor[i][i] = 1;
          for (unsigned int j = 0; j < dim; ++j)
          {
            tensor[i][j] -= normal_vector[i] * normal_vector[j];
          }
        }
        break;
      }
      default:
      {
        Assert(false, ExcNotImplemented());
      }
    }
    return tensor;
  }


  /*
   * This function indicates, if the layer is at the quadrature point
   * active. This will be checked with the distance \xi=\vec{n}(\vec{x}-\vec{x_0}) > 0 ?
   *
   * Input:
   *  point: the coordinates of the quadrature point
   *
   * Output:
   *  flag: false = layer is not active
   *        true  = layer is active
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  bool InterfaceDataPML<dim, Number>::is_active(const Point<dim> point) const
  {
    double tmp = 0.0;
    Tensor<1, dim, value_type> vec_diff;
    for (unsigned int i = 0; i < dim; ++i) vec_diff[i] = point(i) - x_0(i);

    switch (geometry)
    {
      case hyper_cube:
      {
        tmp = vec_diff * normal_vector;
        break;
      }
      case hyper_sphere_radial:
      case hyper_cylinder_radial:
      case hyper_sphere_tangential:
      case hyper_cylinder_tangential:
      {
        tmp = (vec_diff * get_normal_vector(point)) - R_0;
        break;
      }
      case hyper_sphere_radial_inv:
      case hyper_sphere_tangential_inv:
      {
        tmp = R_0 - (vec_diff * get_normal_vector(point));
        break;
      }
      default:
      {
        Assert(false, ExcMessage("false geometry"));
      }
    }

    // if(point(0)<-0.5&&point(1)>0.5)
    //  tmp=1.0;

    return (tmp <= 0 ? false : true);
  }


  /*
   * This function returns a enum of the layer geometry.
   *
   * Input:
   *  none
   *
   * Output:
   *  geometry: which geometry is set
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  enum OuterGeometry InterfaceDataPML<dim, Number>::get_geometry() const
  {
    return geometry;
  }



  /*
   * Standard Constructor of the class AttenuationPML.
   * A quadratic profil function is set by default.
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  AttenuationPML<dim, Number>::AttenuationPML() : attenuation_function(quadratic), n_layer(0)
  {
  }


  /*
   * This function returns the number of layer which are
   * globally available. The circular and spherical pml
   * is split into a radial and tangential component.
   *
   * Input:
   *  none
   *
   * Output:
   *  number of global layer
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  unsigned int AttenuationPML<dim, Number>::get_n_layer()
  {
    return n_layer;
  }


  /*
   * This function indicates, if the layer is at the quadrature point
   * active. This will be checked with the distance \xi=\vec{n}(\vec{x}-\vec{x_0}) > 0 ?
   *
   * Input:
   *  point: the coordinates of the quadrature point
   *
   * Output:
   *  flag: false = layer is not active
   *        true  = layer is active
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  bool AttenuationPML<dim, Number>::is_layer_active(unsigned int layer, Point<dim> point)
  {
    if (layer > n_layer) dserror("layer number is not valid!");

    return interface_pml[layer].is_active(point);
  }


  /*
   * This function read a pml definition from a file.
   * The input for the pml definition should look something like this:
   * PML X0  1  0 0 N_VEC  1  0 0  WIDTH 1 SIGMA 10  FKT CONSTANT
   * PML X0  0  1 0 N_VEC  0  1 0  WIDTH 1 SIGMA 10  FKT LINEAR
   * PML X0 -1  0 0 N_VEC -1  0 0  WIDTH 1 SIGMA 10  FKT QUADRATIC
   * PML X0  0 -1 0 N_VEC  0 -1 0  WIDTH 1 SIGMA 10  FKT CUBIC
   * // COMMENT ...
   * PML X0  1  0 0 N_VEC  1  0 0  WIDTH 1 SIGMA 10  FKT QUARTIC
   * PML X0  0  1 0 R_SPH  1       WIDTH 1 SIGMA 10  FKT HYPBERBOLIC
   * PML X0  0  1 0 R_CYL  1       WIDTH 1 SIGMA 10  FKT SHIFTED_HYPBERBOLIC
   * PML X0  0  0 0 R_INV  1       WIDTH 1 SIGMA 10  FKT QUADRATIC
   *
   * Input:
   *  filename: filename of the file which contain the pml definition
   *
   * Output:
   *  none
   *
   * Author: Matthias Kufner
   * Date: 04.04.2016
   */
  template <int dim, typename Number>
  void AttenuationPML<dim, Number>::read_pml_definition(std::string filename)
  {
    if (filename != std::string("none"))
    {
      // insert path to file if necessary
      if (filename[0] != '/')
      {
        std::string pathfilename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
        std::string::size_type pos = pathfilename.rfind('/');
        if (pos != std::string::npos)
        {
          std::string path = pathfilename.substr(0, pos + 1);
          filename.insert(filename.begin(), path.begin(), path.end());
        }
      }


      // read the input file
      std::ifstream inputfile(filename.c_str());
      if (!inputfile)
      {
        std::cout
            << (std::string("file \"") + filename + std::string("\" could not be opened")).c_str()
            << std::endl;
        dserror("");
        return;
      }

      // read the pml data from the input file
      std::string word = std::string("");
      int error_flag = 0;
      n_layer = 0;
      InterfaceDataPML<dim, Number> interface_tmp;
      while ((!inputfile.eof()) && (error_flag == 0))
      {
        word = std::string("");
        // key word to define a pml
        inputfile >> word;
        if (word == std::string("PML"))
        {
          // read the reference point
          inputfile >> word;
          if (word == std::string("X0"))
          {
            Point<dim> point_tmp;
            double tmp;
            for (unsigned int i = 0; i < 3; ++i)
            {
              if (i < dim)
                inputfile >> point_tmp(i);
              else
                inputfile >> tmp;
            }

            interface_tmp.set_x0(point_tmp);
          }
          else
            error_flag |= 1;

          // read the normal vector for a rectilinear layer or the radius
          // for a circular or spherical layer
          inputfile >> word;
          if ((word == std::string("N_VEC")) || (word == std::string("R_SPH")) ||
              (word == std::string("R_CYL")) || (word == std::string("R_INV")))
          {
            if (word == std::string("N_VEC"))
            {
              interface_tmp.set_geometry(hyper_cube);

              Tensor<1, dim, value_type> normal_vector_tmp;
              double tmp;
              for (unsigned int i = 0; i < 3; ++i)
              {
                if (i < dim)
                  inputfile >> normal_vector_tmp[i];
                else
                  inputfile >> tmp;
              }

              double norm = normal_vector_tmp.norm();
              if (norm > 1e-14)
                interface_tmp.set_normal_vector(normal_vector_tmp / norm);
              else
                dserror("norm of normal vector must be greater than zero!");
            }
            else if (word == std::string("R_SPH"))
            {
              interface_tmp.set_geometry(hyper_sphere);
              double tmp;
              inputfile >> tmp;
              interface_tmp.set_R0(tmp);
            }
            else if (word == std::string("R_CYL"))
            {
              interface_tmp.set_geometry(hyper_cylinder);
              double tmp;
              inputfile >> tmp;
              interface_tmp.set_R0(tmp);
            }
            else if (word == std::string("R_INV"))
            {
              interface_tmp.set_geometry(hyper_sphere_inv);
              double tmp;
              inputfile >> tmp;
              interface_tmp.set_R0(tmp);
            }
          }
          else
            error_flag |= 1;

          // read the width of the layer
          inputfile >> word;
          if (word == std::string("WIDTH"))
          {
            double tmp;
            inputfile >> tmp;
            interface_tmp.set_width(tmp);
          }
          else
            error_flag |= 1;

          // read the factor for the attenuation
          inputfile >> word;
          if (word == std::string("SIGMA"))
          {
            double tmp;
            inputfile >> tmp;
            interface_tmp.set_attenuation(tmp);
          }
          else
            error_flag |= 1;

          // read the choice of the attenuation function
          inputfile >> word;
          if (word == std::string("FKT"))
          {
            inputfile >> word;
            if (word == std::string("CONSTANT"))
              interface_tmp.set_attenuation_fct(constant);
            else if (word == std::string("LINEAR"))
              interface_tmp.set_attenuation_fct(linear);
            else if (word == std::string("QUADRATIC"))
              interface_tmp.set_attenuation_fct(quadratic);
            else if (word == std::string("CUBIC"))
              interface_tmp.set_attenuation_fct(cubic);
            else if (word == std::string("QUARTIC"))
              interface_tmp.set_attenuation_fct(quartic);
            else if (word == std::string("HYPBERBOLIC"))
              interface_tmp.set_attenuation_fct(hyperbolic);
            else if (word == std::string("SHIFTED_HYPERBOLIC"))
              interface_tmp.set_attenuation_fct(shifted_hyperbolic);
            else
              error_flag |= 1;
          }
          else
            error_flag |= 1;

          // save the layer
          switch (interface_tmp.get_geometry())
          {
            case hyper_cube:
            {
              interface_pml.push_back(InterfaceDataPML<dim, Number>());
              interface_pml[n_layer] = interface_tmp;
              n_layer++;
              break;
            }
            case hyper_sphere:
            {
              // first layer in radial direction
              interface_pml.push_back(InterfaceDataPML<dim, Number>());
              interface_pml[n_layer] = interface_tmp;
              interface_pml[n_layer].set_geometry(hyper_sphere_radial);
              // second layer in tangential direction
              n_layer++;
              interface_pml.push_back(InterfaceDataPML<dim, Number>());
              interface_pml[n_layer] = interface_tmp;
              interface_pml[n_layer].set_geometry(hyper_sphere_tangential);
              n_layer++;
              break;
            }
            case hyper_cylinder:
            {
              // first layer in radial direction
              interface_pml.push_back(InterfaceDataPML<dim, Number>());
              interface_pml[n_layer] = interface_tmp;
              interface_pml[n_layer].set_geometry(hyper_cylinder_radial);
              // second layer in tangential direction
              n_layer++;
              interface_pml.push_back(InterfaceDataPML<dim, Number>());
              interface_pml[n_layer] = interface_tmp;
              interface_pml[n_layer].set_geometry(hyper_cylinder_tangential);
              n_layer++;
              break;
            }
            case hyper_sphere_inv:
            {
              // first layer in radial direction
              interface_pml.push_back(InterfaceDataPML<dim, Number>());
              interface_pml[n_layer] = interface_tmp;
              interface_pml[n_layer].set_geometry(hyper_sphere_radial_inv);
              // second layer in tangential direction
              n_layer++;
              interface_pml.push_back(InterfaceDataPML<dim, Number>());
              interface_pml[n_layer] = interface_tmp;
              interface_pml[n_layer].set_geometry(hyper_sphere_tangential_inv);
              n_layer++;
              break;
            }
            default:
            {
              Assert(false, ExcNotImplemented());
            }
          }
        }
        else
        {
          if (word.substr(0, 2) == std::string("//"))
          {
            char line[1024];
            inputfile.getline(line, 1024);
          }
          else if (word != std::string(""))
            error_flag |= 1;
        }
      }
      // output an example for the definition of the pml data
      if ((!inputfile.eof()) || (error_flag != 0))
      {
        std::cout << (std::string("something went wrong reading from file ") + filename).c_str();
        std::cout << std::endl;
        std::cout << "the input for the pml definition should look something like this: "
                  << std::endl;
        std::cout << "PML X0  1  0 0 N_VEC  1  0 0  WIDTH 1 SIGMA 10  FKT CONSTANT \n"
                  << "PML X0  0  1 0 N_VEC  0  1 0  WIDTH 1 SIGMA 10  FKT LINEAR \n"
                  << "PML X0 -1  0 0 N_VEC -1  0 0  WIDTH 1 SIGMA 10  FKT QUADRATIC \n"
                  << "PML X0  0 -1 0 N_VEC  0 -1 0  WIDTH 1 SIGMA 10  FKT CUBIC \n"
                  << "// COMMENT ... \n"
                  << "PML X0  1  0 0 N_VEC  1  0 0  WIDTH 1 SIGMA 10  FKT QUARTIC \n"
                  << "PML X0  0  1 0 R_SPH  1       WIDTH 1 SIGMA 10  FKT HYPBERBOLIC \n"
                  << "PML X0  0  1 0 R_CYL  1       WIDTH 1 SIGMA 10  FKT SHIFTED_HYPERBOLIC \n"
                  << "PML X0  0  0 0 R_INV  1       WIDTH 1 SIGMA 10  FKT QUADRATIC \n"
                  << std::endl;
        dserror("you have to use this notation");
        return;
      }

      // reserve memory for the orthogonality map
      for (unsigned int i = 0; i < n_layer; ++i)
        orthogonality.push_back(std::vector<int>(n_layer, 1));

      // write the orthogonality map
      for (unsigned int i = 0; i < n_layer; ++i)
      {
        orthogonality[i][i] = 0;
        for (unsigned int j = i + 1; j < n_layer; ++j)
        {
          if ((interface_pml[i].get_geometry() == hyper_cube) &&
              (interface_pml[j].get_geometry() == hyper_cube))
          {
            double scalar_product =
                interface_pml[i].get_normal_vector() * interface_pml[j].get_normal_vector();
            if (std::abs(scalar_product) < 1e-14) orthogonality[i][j] = 0;
          }
          else if ((interface_pml[i].get_geometry() == hyper_sphere_radial) &&
                   (interface_pml[j].get_geometry() == hyper_sphere_tangential))
          {
            orthogonality[i][j] = 0;
          }
          else if ((interface_pml[i].get_geometry() == hyper_cylinder_radial) &&
                   (interface_pml[j].get_geometry() == hyper_cylinder_tangential))
          {
            orthogonality[i][j] = 0;
          }
          else if ((interface_pml[i].get_geometry() == hyper_sphere_radial_inv) &&
                   (interface_pml[j].get_geometry() == hyper_sphere_tangential_inv))
          {
            orthogonality[i][j] = 0;
          }
          else if (((interface_pml[i].get_geometry() == hyper_cube) &&
                       ((interface_pml[j].get_geometry() & hyper_cylinder) != 0)) ||
                   ((interface_pml[j].get_geometry() == hyper_cube) &&
                       ((interface_pml[i].get_geometry() & hyper_cylinder) != 0)))
          {
            int flag = 0;
            for (unsigned int i = 0; i < 2; ++i)
              if (std::abs(interface_pml[i].get_normal_vector()[i]) > 1e-14) flag |= 1;
            if (flag)
              orthogonality[i][j] = 1;
            else
              orthogonality[i][j] = 0;
          }
        }
      }

      // mirrow the data, because the matrix is symmetric
      for (unsigned int i = 0; i < n_layer; ++i)
        for (unsigned int j = i + 1; j < n_layer; ++j) orthogonality[j][i] = orthogonality[i][j];

      // debug -> write out the orthogonality map
      /*std::cout << "orthogonality: " << std::endl;
      for (unsigned int i = 0; i < n_layer; ++i)
      {
        for (unsigned int j = 0; j < n_layer; ++j)
          std::cout << orthogonality[i][j] << "  ";
        std::cout << std::endl;
      }*/

      inputfile.close();
    }
    else
      dserror("no pml definition file provided");

    return;
  }

  /*
   * This function calculates the eigenvalues \mu_j, the matrix A and
   * the tensors N\otimesN depending on the position of the quadrature point.
   * The layer_reference holds the indices of the layer which are involved
   * for the calculation.
   *
   * Input:
   *  layer_reference: this vector holds the indices of the layer, which are involved in the element
   *  vec_array_no: index of the VectorizedArray variables
   *  position: coordinates of the quadrature point
   *
   * Output:
   *  sigma_values: this vector holds the values of the scalar profil function sigma the
   *                position of the quadrature point
   *  eigenvalues: this vector holds the eigenvalues of the tensor A
   *  Matrix_A: this tensor includes the tensor function A
   *  eigen_tensors: this vector holds the tensors N\otimesP
   *
   *  Author: Matthias Kufner
   *  Date: 04.04.2016
   */
  template <int dim, typename Number>
  void AttenuationPML<dim, Number>::get_matrix(const std::vector<int> &layer_reference,
      const unsigned int vec_array_no, const Point<dim> &position,
      Tensor<1, dim, VectorizedArray<value_type>> &sigma_values,
      Tensor<1, dim, VectorizedArray<value_type>> &eigenvalues,
      Tensor<2, dim, VectorizedArray<value_type>> &Matrix_A,
      Tensor<3, dim, VectorizedArray<value_type>> &eigen_tensors) const
  {
    unsigned int n_layer = layer_reference.size();
    Tensor<2, dim, value_type> tensor_tmp;

    // initialize the variables
    for (unsigned int n = 0; n < dim; ++n)
    {
      sigma_values[n][vec_array_no] = 0.0;
      eigenvalues[n][vec_array_no] = 0.0;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
        {
          Matrix_A[i][j][vec_array_no] = 0.0;
          eigen_tensors[n][i][j][vec_array_no] = 0.0;
        }
    }

    // check if all involved layer are orthogonal
    unsigned int ortho_flag = 0;
    for (unsigned int i = 0; i < n_layer; ++i)
      for (unsigned int j = i + 1; j < n_layer; ++j)
        if (orthogonality[layer_reference[i]][layer_reference[j]]) ortho_flag |= 1;

    // ortho_flag = 1; // debug: use only the spectral decomposition
    // if all layer are orthogonal, then use the analytic solution
    if (!ortho_flag)
    {
      for (unsigned int n = 0; n < n_layer; ++n)
      {
        sigma_values[n][vec_array_no] = interface_pml[layer_reference[n]].get_sigma(position);
        tensor_tmp = interface_pml[layer_reference[n]].get_tensor(position);
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
          {
            eigen_tensors[n][i][j][vec_array_no] = tensor_tmp[i][j];
            Matrix_A[i][j][vec_array_no] += sigma_values[n][vec_array_no] * tensor_tmp[i][j];
          }
      }
      for (unsigned int i = 0; i < dim; ++i)
        eigenvalues[i][vec_array_no] = sigma_values[i][vec_array_no];
    }
    // if the layer are not orthogonal then use the spectral decomposition
    else if (dim == 2)
    {
      std::vector<value_type> sigma_tmp(n_layer, 0.0);
      for (unsigned int n = 0; n < n_layer; ++n)
      {
        sigma_tmp[n] = interface_pml[layer_reference[n]].get_sigma(position);
        tensor_tmp = interface_pml[layer_reference[n]].get_tensor(position);
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            Matrix_A[i][j][vec_array_no] += sigma_tmp[n] * tensor_tmp[i][j];
      }

      // calculate the eigenvalues and eigenvectors with the mohr's circle
      value_type tmp, phi, R, x_M;
      tmp = Matrix_A[1][1][vec_array_no] - Matrix_A[0][0][vec_array_no];
      phi = -0.5 * std::atan2(Matrix_A[0][1][vec_array_no], 0.5 * tmp);
      tmp *= 0.25 * tmp;
      tmp += Matrix_A[0][1][vec_array_no] * Matrix_A[0][1][vec_array_no];
      R = std::sqrt(tmp);
      x_M = 0.5 * (Matrix_A[1][1][vec_array_no] + Matrix_A[0][0][vec_array_no]);

      // eigenvalues
      eigenvalues[0][vec_array_no] = x_M + R;
      eigenvalues[1][vec_array_no] = x_M - R;

      // eigenvectors
      Tensor<2, dim, value_type> eigenvectors;
      eigenvectors[0][0] = std::sin(phi);
      eigenvectors[1][0] = std::cos(phi);
      eigenvectors[0][1] = -eigenvectors[1][0];
      eigenvectors[1][1] = +eigenvectors[0][0];

      for (unsigned int n = 0; n < dim; ++n)
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            eigen_tensors[n][i][j][vec_array_no] =
                eigenvalues[n][vec_array_no] * eigenvectors[n][i] * eigenvectors[n][j];


      for (unsigned int i = 0; i < dim; ++i) sigma_values[i][vec_array_no] = 1.0;
    }
    // the spectral decomposition is still implemented and tested in 2d
    else
    {
      std::cout << "not orthogonal" << std::endl;
      dserror("spectral decomposition only implemented for 2D");
    }
  }



  template class AttenuationPML<2, double>;
  template class AttenuationPML<3, double>;
  template class AttenuationPML<2, float>;
  template class AttenuationPML<3, float>;

}  // end namespace ACOU
#endif  // HAVE_DEAL_II
