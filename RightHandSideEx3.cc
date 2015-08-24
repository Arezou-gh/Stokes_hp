#include "RightHandSideEx3.hh"

#include <cmath>

using std::exp;
using std::cos;
using std::sin;
using std::pow;


template <int dim>
double RightHandSideEx3<dim>::value(const Point<dim> & p, const unsigned int component) const
{
  const double x(p[0]);
  const double y(p[1]);
  double val(0.);

  if (component==0)
    val = 4.*y*(sin(pow(x,2)+pow(y,2)) + 2.*pow(x,2)*cos(pow(x,2)+pow(y,2))) - 
      20.*x*exp(-10.*(pow(x,2)+pow(y,2)));
  if (component==1)
    val = -4.*x*(sin(pow(x,2)+pow(y,2)) + 2.*pow(y,2)*cos(pow(x,2)+pow(y,2))) -
      20.*y*exp(-10.*(pow(x,2)+pow(y,2)));

  return val;
}


template <int dim>
void
RightHandSideEx3<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = this->value(p,c);
}


//Explicit initialization
template class RightHandSideEx3<2>;
template class RightHandSideEx3<3>;