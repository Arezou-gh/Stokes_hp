#include "Parameters.hh"

#include <fstream>
#include <string>

Parameters::Parameters(const std::string &input_filename)
{
  ParameterHandler prm;

  // Declate the parameters
  declare_parameters(prm);

  // Read the file
  std::ifstream parameters_file(input_filename.c_str());
  AssertThrow(parameters_file,ExcMessage("Input parameters file not found."));
  const bool success = prm.read_input(parameters_file);
  AssertThrow(success,ExcMessage("Invalid input parameters file."));
  parameters_file.close();

  // Parse the parameters
  parse_parameters(prm);
}


void Parameters::declare_parameters(ParameterHandler &prm)
{
  prm.declare_entry("Verbose", "false", Patterns::Bool(), "Verbose output");
  prm.declare_entry("Exemple", "1", Patterns::Integer(1,3), "Exemple to run");
  prm.declare_entry("Quadrature", "GaussLegendre", Patterns::Selection(
        "GaussLegendre|GaussLobatto"), "Type of quadrature to use");
  prm.declare_entry("Max degree", "10", Patterns::Integer(2), "Maximum degree used");
  prm.declare_entry("Max n cycles", "10", Patterns::Integer(0), "Maximum degree used");
  prm.declare_entry("Theta", "0.5", Patterns::Double(0.,1.), "Refinement parameters");
  prm.declare_entry("Tolerance", "1e-12", Patterns::Double(0.), "Solver tolerance");
}


void Parameters::parse_parameters(ParameterHandler &prm)
{
  verbose = prm.get_bool("Verbose");

  unsigned int ex_number = prm.get_integer("Exemple");
  switch (ex_number)
  {
    case 1: 
      {
        exemple = exemple_1;
        dim = 2;
        break;
      }
    case 2:
      {
        exemple = exemple_2;
        dim = 2;
        break;
      }
    case 3:
      {
        exemple = exemple_3;
        dim = 2;
        break;
      }
    default:
      {
        AssertThrow(false,ExcMessage("Unknow exemple"));
      }
  }

  std::string input = prm.get("Quadrature");
  if (input.compare("GaussLegendre")==0)
    quadrature = gauss_legendre;
  else
    quadrature = gauss_lobatto;

  max_degree = prm.get_integer("Max degree");

  max_n_cycles = prm.get_integer("Max n cycles");
  
  theta = prm.get_double("Theta");
  
  tolerance = prm.get_double("Tolerance");
}