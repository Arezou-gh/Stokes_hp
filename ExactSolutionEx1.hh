#ifndef _EXACTSOLUTIONEX1_HH_
#define _EXACTSOLUTIONEX1_HH_

#include <vector>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

template <int dim>
class ExactSolutionEx1 : public Function<dim>
{
  public:
	ExactSolutionEx1 () : Function<dim>(dim+1) {}
    virtual void vector_value (const Point<dim> &p,
        Vector<double>   &value) const;
    virtual void vector_gradient (const Point<dim> &p,
        std::vector< Tensor< 1, dim > > &gradients) const;
};

#endif
