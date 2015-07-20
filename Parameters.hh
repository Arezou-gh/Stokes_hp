#ifndef _PARAMETERS_HH_
#define _PARAMETERS_HH_

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

// Exemple 1 is the first exemple of
// Exemple 2 is the second exemple of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
// Exemple 3 is the first exemple of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
enum EXEMPLE{exemple_1,exemple_2,exemple_3};
enum QUADRATURE{gauss_legendre,gauss_lobatto};

class Parameters
{
  public :
    Parameters(const std::string &input_filename);

    bool get_verbose() const;

    EXEMPLE get_exemple() const;

    QUADRATURE get_quadrature() const;

    unsigned int get_dim() const;

    unsigned int get_max_degree() const;

    unsigned int get_max_n_cycles() const;

    double get_theta() const;

    double get_tolerance() const;

  private :
    void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);

    ParameterHandler prm;
    bool verbose;
    EXEMPLE exemple;
    QUADRATURE quadrature;
    unsigned int dim;
    unsigned int max_degree;
    unsigned int max_n_cycles;
    double theta;
    double tolerance;
};


inline
bool Parameters::get_verbose() const
{
  return verbose;
}


inline 
EXEMPLE Parameters::get_exemple() const
{
  return exemple;
}

inline
QUADRATURE Parameters::get_quadrature() const
{
  return quadrature;
}


inline
unsigned int Parameters::get_dim() const
{
  return dim;
}


inline
unsigned int Parameters::get_max_degree() const
{
  return max_degree;
}


inline
unsigned int Parameters::get_max_n_cycles() const
{
  return max_n_cycles;
}


inline
double Parameters::get_theta() const
{
  return theta;
}


inline
double Parameters::get_tolerance() const
{
  return tolerance;
}

#endif
