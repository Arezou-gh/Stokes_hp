#ifndef _RIGHTHANDSIDEEx3_HH_
#define _RIGHTHANDSIDEEx3_HH_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

template <int dim>
class RightHandSideEx3 : public Function<dim>
{
  public:
    RightHandSideEx3() : Function<dim>(dim+1) {}
    virtual double value (const Point<dim> &p, const unsigned int component) const;
    virtual void vector_value (const Point<dim> &p, Vector<double> &value) const;
};

#endif
