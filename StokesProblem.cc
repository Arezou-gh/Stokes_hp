
//find the best (optimum) theta
// give the error-ErrEstimator plot for the optimum theta
// give the error-theta graph, for different theta's
// give the h- , p- and hp- error-dofs plot for the optimum theta which you found.
//test the code with the new example in Houston's paper

#include "StokesProblem.hh"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include "SchurComplement.hh"
//#include "RightHandSideLocal.hh"

/*......................................................................................*/
// constructor

template <int dim>
StokesProblem<dim>::StokesProblem(int exemple): dof_handler(triangulation),  max_degree (10), Tolerance (1e-16)

{
	if (exemple==1)
		exact_solution = new ExactSolutionEx1<dim>();
	else
		exact_solution = new ExactSolutionEx2<dim>();
	for (unsigned int degree=1; degree<=max_degree; ++degree)
	{

/*
		// GaussLobatto polynomials and quadratures:


		fe_collection.push_back (FESystem<dim>(FE_Q<dim> (QGaussLobatto<1> (degree + 2)), dim,
				FE_Q<dim> (QGaussLobatto<1> (degree + 1)), 1));

		quadrature_collection.push_back(QGaussLobatto<dim> (degree+3));
		face_quadrature_collection.push_back ( QGaussLobatto<dim-1> (degree+3));

		quadrature_collection_Err.push_back( QGaussLobatto<dim> (degree+4));
		face_quadrature_collection_Err.push_back ( QGaussLobatto<dim-1> (degree+4));
*/

		fe_collection.push_back (FESystem<dim>(FE_Q<dim> (degree + 2), dim,
				FE_Q<dim> (degree+1), 1));

		quadrature_collection.push_back(QGauss<dim> (degree+3));
		face_quadrature_collection.push_back (QGauss<dim-1> (degree+3));

		quadrature_collection_Err.push_back(QGauss<dim> (degree+4));
		face_quadrature_collection_Err.push_back (QGauss<dim-1> (degree+4));
	}
	fe_collection.push_back (FESystem<dim>(FE_Nothing<dim>(), dim,
			FE_Nothing<dim>(), 1));
	quadrature_collection.push_back(QGauss<dim>(1));
	face_quadrature_collection.push_back (QGauss<dim-1>(1));
}

/*.....................................................................................*/
template <int dim>
StokesProblem <dim>::~StokesProblem()
{
 dof_handler.clear();
 delete exact_solution;
}

/*......................................................................................*/
template <int dim>
bool StokesProblem <dim>::decreasing (const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &i, const std::pair<double,typename hp::DoFHandler<dim>::active_cell_iterator > &j)
{
 return ((i.first) > (j.first));
}

/*......................................................................................*/
// Generate mesh


template <int dim>
void StokesProblem <dim>::generate_mesh(){

/*
 GridGenerator::hyper_cube (triangulation, -1, 1);
 triangulation.refine_global (1);

 std::ofstream out ("grid-hyper_cube.eps");
 GridOut grid_out;
 grid_out.write_eps (triangulation, out);
 */


 std::vector<Point<dim> > vertices (8);

 vertices [0]=Point<dim> (-1,-1);
 vertices [1]=Point<dim> (0,-1);
 vertices [2]=Point<dim> (-1,0);
 vertices [3]=Point<dim> (0,0);
 vertices [4]=Point<dim> (1,0);
 vertices [5]=Point<dim> (-1,1);
 vertices [6]=Point<dim> (0,1);
 vertices [7]=Point<dim> (1,1);

 const unsigned int n_cells=3;
 std::vector<CellData<dim> > cell(n_cells);
 cell[0].vertices[0]=0;
 cell[0].vertices[1]=1;
 cell[0].vertices[2]=2;
 cell[0].vertices[3]=3;

 cell[1].vertices[0]=2;
 cell[1].vertices[1]=3;
 cell[1].vertices[2]=5;
 cell[1].vertices[3]=6;

 cell[2].vertices[0]=3;
 cell[2].vertices[1]=4;
 cell[2].vertices[2]=6;
 cell[2].vertices[3]=7;

 triangulation.create_triangulation(vertices,cell,SubCellData());
 triangulation.refine_global (1);
 std::ofstream out ("grid-L-Shape.eps");
 GridOut grid_out;
 grid_out.write_eps (triangulation, out);

 // std::cout<<"Number of active cells: "<< triangulation.n_active_cells() << std::endl;
 // std::cout<<"Total number of cells: " << triangulation.n_cells() << std::endl ;
}
//.....................................................................................
template <int dim>
void StokesProblem <dim>::set_global_active_fe_indices (hp::DoFHandler<dim> &dof_handler)
{
 typename hp::DoFHandler<dim>::active_cell_iterator cell= dof_handler.begin_active(), end_cell = dof_handler.end();
 for (; cell!=end_cell; ++cell)

   cell->set_active_fe_index (0);
}
/*......................................................................................*/
// setup system()

template <int dim>
void StokesProblem <dim>::setup_system(){

	system_matrix.clear();

	// set_global_active_fe_indices(dof_handler);

	dof_handler.distribute_dofs (fe_collection);

	DoFRenumbering::Cuthill_McKee (dof_handler);

	std::vector<unsigned int> block_component (dim+1, 0);
	block_component[dim]=1;
	DoFRenumbering::component_wise(dof_handler, block_component);


	//const double evaluate_pressure=  2.0 * std::exp (-1) * std::sin (-1) - (2.0 * (1.0 - std::exp(1.0)) * (std::cos (1.0) - 1.0)) / 3.0; ;
	// std::cout<< "evaluate_pressure: " << evaluate_pressure;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//


	{
		constraints.clear ();
		FEValuesExtractors::Vector velocities(0);
		DoFTools::make_hanging_node_constraints (dof_handler, constraints);
		VectorTools::interpolate_boundary_values (dof_handler,0,(* exact_solution),constraints , fe_collection.component_mask(velocities));


		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
		// Since with Dirichlet velocity Bdry, pressure will be defined up to a constant, in order to make the
		//solution to be unique, we need to add an additional constraint for Pressure
		// We choose for example the first cell of triangulation and do as follow:
		//
		typename hp::DoFHandler<dim>::active_cell_iterator
		first_cell = dof_handler.begin_active();
		//  std::cout<< "vertex _0  of the first_cell   " << first_cell->vertex(0) << std::endl;
		std::vector<types::global_dof_index> local_dof_indices (first_cell->get_fe().dofs_per_cell);
		first_cell->get_dof_indices(local_dof_indices);
		Point<dim> Pnt_in_ref = first_cell->get_fe().unit_support_point(first_cell->get_fe().component_to_system_index(dim,0));
		MappingQ1<dim> mapping;
		Point<dim> Pnt_in_real = mapping.transform_unit_to_real_cell(first_cell,Pnt_in_ref);

		types::global_dof_index first_pressure_dof = local_dof_indices[first_cell->get_fe().component_to_system_index(dim,0)];
		//std::cout<< "first_pressure_dof" << first_pressure_dof << std::endl;

		// component_to_system_index: "Compute the shape function for the given vector component and index."
		constraints.add_line (first_pressure_dof);
		Vector<double> values(3);
		 exact_solution->vector_value(Pnt_in_real,values);
		// std::cout<< "Pnt_in_real :" << Pnt_in_real<<" "<<values[2]<< std::endl;
		constraints.set_inhomogeneity (first_pressure_dof,values[dim]);

		//  constraints.set_homog (local_dof_indices[fe.comp_to_system_i(dim,0)],
		//  	       exact_pressure(dof_handler.begin_active()->vertex(0)));

		//****************/

	}
	constraints.close();
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	std::vector<types::global_dof_index> dofs_per_block (2);
	DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
	const unsigned int n_u = dofs_per_block[0],
			n_p = dofs_per_block[1];
	// std::cout << "   Number of active cells: "
	//        << triangulation.n_active_cells()
	//        << std::endl

	std::cout   << "   Number of degrees of freedom: "
			<< dof_handler.n_dofs()
			<< " (" << n_u << '+' << n_p << ')'
			<< std::endl;
	{
		BlockCompressedSimpleSparsityPattern csp (2,2);
		csp.block(0,0).reinit (n_u, n_u);
		csp.block(1,0).reinit (n_p, n_u);
		csp.block(0,1).reinit (n_u, n_p);
		csp.block(1,1).reinit (n_p, n_p);
		csp.collect_sizes();
		DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
		sparsity_pattern.copy_from (csp);
	}
	system_matrix.reinit (sparsity_pattern);
	solution.reinit (2);
	solution.block(0).reinit (n_u);
	solution.block(1).reinit (n_p);
	solution.collect_sizes ();
	system_rhs.reinit (2);
	system_rhs.block(0).reinit (n_u);
	system_rhs.block(1).reinit (n_p);
	system_rhs.collect_sizes ();
}
/*......................................................................................*/
// assemble system

template <int dim>
void StokesProblem <dim>::assemble_system () {
	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients);

	FullMatrix<double> local_matrix;
	Vector<double> local_rhs;
	std::vector<types::global_dof_index> local_dof_indices;

	//const RightHandSide<dim> rhs_function;
	std::vector<Vector<double> >  rhs_values;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);


	//std::vector<SymmetricTensor<2,dim> > symgrad_phi_u;
	std::vector<Tensor<2,dim> > grad_phi_u;
	std::vector<double> div_phi_u;
	std::vector<Tensor<1,dim> > phi_u;
	std::vector<double> phi_p;

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
		local_matrix.reinit (dofs_per_cell, dofs_per_cell);
		local_rhs.reinit (dofs_per_cell);
		local_matrix=0;
		local_rhs=0;

		hp_fe_values.reinit (cell);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;

		rhs_values.resize(n_q_points, Vector<double>(dim+1));
		rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

		//symgrad_phi_u.resize(dofs_per_cell);
		grad_phi_u.resize(dofs_per_cell);
		div_phi_u.resize(dofs_per_cell);
		phi_u.resize (dofs_per_cell);
		phi_p.resize(dofs_per_cell);

		for (unsigned int q=0; q<n_q_points; ++q)
		{
			for (unsigned int k=0; k<dofs_per_cell; ++k)
			{
				//symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
				grad_phi_u[k] = fe_values[velocities].gradient (k, q);
				div_phi_u[k] = fe_values[velocities].divergence (k, q);
				phi_u[k] = fe_values[velocities].value (k, q);
				phi_p[k] = fe_values[pressure].value (k, q);
			}

			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				for (unsigned int j=0; j<dofs_per_cell; ++j)
				{
					// assemble the local system, like step-22
					//local_matrix(i,j) += (symgrad_phi_u[i] * symgrad_phi_u[j]
					// - div_phi_u[i] * phi_p[j]
					// - phi_p[i] * div_phi_u[j])
					//   * fe_values.JxW(q);



// This is my try to assemble the system matrix using double contract for grad_phi_u....
       		  //and also this is the corresponding system matrix to solve with direct solver!
					local_matrix(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
						- div_phi_u[i] * phi_p[j]
					- phi_p[i] * div_phi_u[j])* JxW_values[q];


/*
					// system matrix corresponding to the iterative solver
					local_matrix(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
							- div_phi_u[i] * phi_p[j]
							                       - phi_p[i] * div_phi_u[j]+ phi_p[i] * phi_p[j])* JxW_values[q];

*/
				} // end of loop over 'j'

				local_rhs(i) += (phi_u[i][0] * rhs_values[q](0) + phi_u[i][1] * rhs_values [q](1)) * JxW_values[q];
			} // end of loop 'i'


		} // end of loop 'q'

		//local system to global system
		local_dof_indices.resize (dofs_per_cell);
		cell->get_dof_indices (local_dof_indices);
		constraints.distribute_local_to_global (local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
	} // end of iteration for cells

}


/*......................................................................................*/
// Solve

// pressure mass matrix to precondition the Schur complement
//  Later, when solving, we then precondition the Schur complement with M^{1}B^{T} where the matrix
//A is related to the Laplace operator
// Thus, solving with A is a lot more complicated: the matrix is badly conditioned and we know that we
//need many iterations unless we have a very good preconditioner

//in 2d, we use the ultimate preconditioner, namely a direct sparse LU decomposition of the matrix.
//This is implemented using the SparseDirectUMFPACK class that uses the UMFPACK direct solver to compute the decomposition.



template <int dim>
void StokesProblem <dim>::solve ()
{
	SparseDirectUMFPACK A_inverse;
	A_inverse.initialize (system_matrix,
			SparseDirectUMFPACK::AdditionalData());
	A_inverse.vmult (solution, system_rhs);


	constraints.distribute (solution);
	//note: for square domain even without normalization of zero pressure mean value, the error orders was attained
	 //basically we did this normalization for L-shaped domains

	solution.block (1).add (-1.0 * pressure_mean_value ());
	constraints.distribute (solution);	
	
}



/*
template <int dim>
void StokesProblem <dim>::solve ()
{
 SparseDirectUMFPACK A_inverse;
 A_inverse.initialize (system_matrix.block(0,0),
     SparseDirectUMFPACK::AdditionalData());
 Vector<double> tmp (solution.block(0).size());
 {
   Vector<double> schur_rhs (solution.block(1).size());
   A_inverse.vmult (tmp, system_rhs.block(0));
   system_matrix.block(1,0).vmult (schur_rhs, tmp);
   schur_rhs -= system_rhs.block(1);

   SchurComplement schur_complement (system_matrix, A_inverse);
   SolverControl solver_control (solution.block(1).size(),
   1e-6*schur_rhs.l2_norm());

  // SolverCG<>    cg (solver_control);
   SolverGMRES<>     gmres (solver_control);
   SparseDirectUMFPACK preconditioner;
   preconditioner.initialize (system_matrix.block(1,1),
	SparseDirectUMFPACK::AdditionalData());
   gmres.solve(schur_complement, solution.block(1), schur_rhs,
			preconditioner);
  // cg.solve (schur_complement, solution.block(1), schur_rhs, preconditioner);
   
   //cout<<" residuals of each step " << solver_control.enable_history_data() << endl;
   constraints.distribute (solution);
   std::cout << "   "
                << solver_control.last_step()
                << " GMRES iterations for Stokes subsystem."
                << std::endl; 
   //  std::cout << "  "
   //<< solver_control.last_step()
   //<< " outer CG Schur complement iterations for pressure"
   //<< std::endl;
 }
 system_matrix.block(0,1).vmult (tmp, solution.block(1));
 tmp *= -1.0;
 tmp += system_rhs.block(0);
 A_inverse.vmult (solution.block(0), tmp);
 constraints.distribute (solution);
 
 // Normalized the solution with keeping pressure's mean value=0.
 solution.block (1).add (-1.0 * pressure_mean_value ());
 constraints.distribute (solution);
}
*/

/*......................................................................................*/
template <int dim>
double StokesProblem <dim>::pressure_mean_value () const
{

	// get pressure such that satisfies mean value property:
	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_JxW_values);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<double> values;
	double domain_mean_val_p=0;
	double measure_domain=0;
	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		hp_fe_values.reinit (cell);

		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;
		values.resize(n_q_points);
		fe_values[pressure].get_function_values(solution, values);
		for (unsigned int q=0; q<n_q_points; ++q)
		{
			domain_mean_val_p += values[q]*JxW_values[q];
			measure_domain += JxW_values[q];
		}//q
	}//cell
	// 3 here is the area corresponding to our 3 cells.
	// return domain_mean_val_p/3.0;
	return domain_mean_val_p / measure_domain;
}

/*......................................................................................*/

template <int dim>
double StokesProblem <dim>::exact_pressure_mean_value () const
{

	// get pressure such that satisfies mean value property:
	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values | update_quadrature_points|update_JxW_values);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<Vector<double> > values;
	double domain_mean_val_p=0;
	double measure_domain=0;
	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	//ExactSolution<dim> exact_solution;
	for (; cell!=endc; ++cell)
	{
		hp_fe_values.reinit (cell);

		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;
		values.resize(n_q_points,Vector<double>(dim+1));
		exact_solution->vector_value_list(fe_values.get_quadrature_points(), values);
		for (unsigned int q=0; q<n_q_points; ++q)
		{
			domain_mean_val_p += values[q][dim]*JxW_values[q];
			measure_domain += JxW_values[q];
		}//q
	}//cell
	// 3 here is the area corresponding to our 3 cells.
	// return domain_mean_val_p/3.0;
	return domain_mean_val_p / measure_domain;
}
/*......................................................................................*/
// compute_error

template <int dim>
void StokesProblem <dim>::compute_error (Vector<double> &error_per_cell, Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err, Vector<double> & Vec_Velocity_Err)
{
	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection_Err, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<double> values;
	std::vector<Tensor<2,dim> > gradients;
	std::vector<Tensor<1,dim> > velocity_values;
	
	std::vector<std::vector<Tensor<1,dim> > > exact_solution_gradients;
	std::vector<Vector<double> > exact_solution_values;

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();

	const double mean_exact_pressure = exact_pressure_mean_value ();
	//std::cout << "*** " << mean_exact_pressure << std::endl;

	unsigned int cell_index=0;
	for (; cell!=endc; ++cell,++cell_index)
	{
		double subtract_p=0.;
		double grad_u_vals=0.;
		double u_vals=0.;
		hp_fe_values.reinit (cell);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const std::vector<Point<dim> >& quadrature_points = fe_values.get_quadrature_points();
		const unsigned int n_q_points =fe_values.n_quadrature_points;

		velocity_values.resize(n_q_points);
		gradients.resize(n_q_points);
		values.resize(n_q_points);
		exact_solution_gradients.resize(n_q_points , std::vector<Tensor<1,dim> > (dim+1));
		exact_solution_values.resize(n_q_points, Vector<double> (dim+1));
 
		fe_values[velocities].get_function_values(solution, velocity_values);
		fe_values[velocities].get_function_gradients(solution, gradients);
		fe_values[pressure].get_function_values(solution, values);

		exact_solution->vector_gradient_list(quadrature_points, exact_solution_gradients);
		exact_solution->vector_value_list(quadrature_points, exact_solution_values);
 
		double diff_laplace_u_grad_p=0;

		for (unsigned int q=0; q<n_q_points; ++q)
		{
			//std::cout<<values[q]<<" ";
			values[q] -= exact_solution_values[q](dim);
			
			// std::cout<< exact_solution_values[q](dim)<<" "<<values[q]<<std::endl;
			subtract_p +=values[q]*values[q]* JxW_values[q];

			for (unsigned int i=0; i<dim; ++i)
			{
				velocity_values[q][i]-=exact_solution_values[q](i);
				gradients[q][i]-=exact_solution_gradients[q][i];
			}

			grad_u_vals += gradients[q].norm_square() * JxW_values[q];
			u_vals += velocity_values[q].norm_square() * JxW_values[q];

		} // q


		error_per_cell(cell_index) = sqrt (subtract_p + grad_u_vals);

		Vect_Pressure_Err(cell_index)=sqrt(subtract_p);
		Vect_grad_Velocity_Err(cell_index)=sqrt(grad_u_vals);
		Vec_Velocity_Err(cell_index)=sqrt(u_vals);

	}// cell
	//  std::cout<< "Vector of Compute Error per Cell: " << error_per_cell<< std::endl ;
	//double L1_norm=error_per_cell.l1_norm();
	//std::cout<< "L1_norm of ERROR is: "<< L1_norm << std::endl;

	std::cout<< std::endl;
	double L2_norm_grad_velocity_Err= Vect_grad_Velocity_Err.l2_norm();
	std::cout<< "L2_norm_grad_velocity_Err : "<< L2_norm_grad_velocity_Err << std::endl;
	std::cout<< std::endl;
	double L2_norm_velocity_Err= Vec_Velocity_Err.l2_norm();
	std::cout<< "L2_norm_velocity_Err : "<< L2_norm_velocity_Err << std::endl;
	std::cout<< std::endl;
	std::cout<< std::endl;
	double L2_norm_pressure_Err=Vect_Pressure_Err.l2_norm();
	std::cout<< "L2_norm_pressure_Err : "<< L2_norm_pressure_Err << std::endl;
	std::cout<< std::endl;
	std::cout<< std::endl;
	double L2_norm_total_Err= sqrt (std::pow (L2_norm_grad_velocity_Err,2)+ std::pow (L2_norm_pressure_Err,2));
	std::cout<< "L2_norm of Tottal_ERROR is : "<< L2_norm_total_Err << std::endl;
	std::cout<< std::endl;

}
/*......................................................................................*/
// compute_estimator

template <int dim>
void StokesProblem <dim>::estimate (Vector<double> &est_per_cell)
{

	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);
	hp::FEFaceValues<dim> hp_fe_face_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
	hp::FEFaceValues<dim> hp_neighbor_face_values(fe_collection, face_quadrature_collection, update_gradients);
	hp::FESubfaceValues<dim> hp_subface_values(fe_collection, face_quadrature_collection, update_JxW_values|update_gradients|update_normal_vectors);
	hp::FESubfaceValues<dim> hp_neighbor_subface_values(fe_collection, face_quadrature_collection, update_gradients);

	std::vector<Tensor<1,dim> > gradients_p;
	std::vector<double> divergences;
	std::vector<Tensor<1,dim> > laplacians;

	std::vector<Tensor<2,dim> > gradients;
	std::vector<Tensor<2,dim> > neighbor_gradients;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	const RightHandSide<dim> rhs_function;
	std::vector<Vector<double> >  rhs_values;

	//est_per_cell.reinit (triangulation.n_active_cells());
	Vector<double> res_est_per_cell(triangulation.n_active_cells());
	Vector<double> Jump_est_per_cell(triangulation.n_active_cells());

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	unsigned int cell_index=0;

	for (; cell!=endc; ++cell,++cell_index)
	{
		hp_fe_values.reinit (cell);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		const unsigned int n_q_points = fe_values.n_quadrature_points;

		rhs_values.resize(n_q_points, Vector<double>(dim+1));
		rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

		divergences.resize(n_q_points);
		gradients_p.resize(n_q_points) ;
		laplacians.resize(n_q_points);

		fe_values[pressure].get_function_gradients(solution, gradients_p);
		fe_values[velocities].get_function_divergences(solution, divergences);
		fe_values[velocities].get_function_laplacians(solution, laplacians);


		double term2=0;//divergence term in estimator definition
		double term1=0;// for the residual term in estimator definition
		for (unsigned int q=0; q<n_q_points; ++q){
			term2 += (divergences[q])*(divergences[q])*JxW_values[q];

			for (unsigned int i=0; i<2; ++i)
				gradients_p[q][i]-= (rhs_values[q](i)+ laplacians[q][i]);

			term1+= contract(gradients_p[q],gradients_p[q])*JxW_values[q];
		}// q
		res_est_per_cell(cell_index)= pow((cell->diameter())/(cell->get_fe().degree), 2.0 ) * (term1) + term2;

		// ........................................... compute jump_est_per_cell..............................................................
		double term3=0;//jumpped part of the estimator
		for (unsigned int face_number=0; face_number<GeometryInfo<2>::faces_per_cell; ++face_number)

			if ((cell->face(face_number)->at_boundary()==false)	&& (cell->face(face_number)->has_children() == false)
					&& (cell->face(face_number)->level() == cell->level()))
			{
				const unsigned int q_index = std::max (cell->active_fe_index(),
						cell->neighbor(face_number)->active_fe_index());

				hp_fe_face_values.reinit       (cell,                        face_number,                             q_index);
				hp_neighbor_face_values.reinit (cell->neighbor(face_number), cell->neighbor_of_neighbor(face_number), q_index);

				const FEFaceValues<2> &neighbor_face_values =hp_neighbor_face_values.get_present_fe_values ();
				const FEFaceValues<2> &fe_face_values = hp_fe_face_values.get_present_fe_values ();

				const std::vector<double>& JxW_values = fe_face_values.get_JxW_values ();

				const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

				gradients.resize(n_face_q_points);
				neighbor_gradients.resize(n_face_q_points);

				neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
				fe_face_values[velocities].get_function_gradients(solution, gradients);

				std::vector<Tensor<1,dim> > jump_per_face;
				jump_per_face.resize(n_face_q_points);
				double jump_val=0;

				for (unsigned int q=0; q<n_face_q_points; ++q)
				{
					for (unsigned int i=0; i<2; ++i){
						for (unsigned int j=0; j<2; ++j){
							jump_per_face[q][i] = (gradients[q][i][j]-neighbor_gradients[q][i][j]) *(fe_face_values.normal_vector(q)[j]);
						}
					}
					jump_val += contract(jump_per_face[q],jump_per_face[q])*JxW_values[q];
				}//q_per_face
				term3 +=(cell->face(face_number)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;
			} //if level

			// else if the neighbor has children

			else if ( (cell->face(face_number)->at_boundary()==false) && (cell->face(face_number)->has_children() == true))
			{
				for (unsigned int subface=0;
						subface< cell->face(face_number)->n_children(); ++subface)
				{
					const unsigned int q_index = std::max(cell->neighbor_child_on_subface (face_number, subface)->active_fe_index(), cell->active_fe_index());

					hp_neighbor_face_values.reinit (cell->neighbor_child_on_subface (face_number, subface), cell->neighbor_of_neighbor(face_number), q_index);
					hp_subface_values.reinit (cell,face_number, subface, q_index);

					const FEFaceValues<2> &neighbor_face_values  = hp_neighbor_face_values.get_present_fe_values ();
					const FESubfaceValues<2> &fe_subface_values = hp_subface_values.get_present_fe_values ();


					const std::vector<double>& JxW_values = fe_subface_values.get_JxW_values ();

					const unsigned int n_subface_q_points = fe_subface_values.n_quadrature_points;//?

							gradients.resize(n_subface_q_points);
					neighbor_gradients.resize(n_subface_q_points);

					neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
					fe_subface_values[velocities].get_function_gradients(solution, gradients);

					std::vector<Tensor<1,dim> > jump_per_subface;
					jump_per_subface.resize(n_subface_q_points);//?

							double jump_val=0;
							for (unsigned int q=0; q<n_subface_q_points; ++q)
							{
								for (unsigned int i=0; i<2; ++i){
									for (unsigned int j=0; j<2; ++j){
										jump_per_subface[q][j] += (gradients[q][i][j]- neighbor_gradients[q][i][j])*(fe_subface_values.normal_vector(q)[j]);
									}//j
								}// i
								jump_val += contract(jump_per_subface[q],jump_per_subface[q])*(JxW_values[q]);
							}//q_per_subface
							term3 +=(cell->face(face_number)->child(subface)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;
				}// subface
			}//else if

		// if the neighbor is coarser

			else if ( (cell->face(face_number)->at_boundary()==false) && (cell->neighbor_is_coarser(face_number)))

			{
				const unsigned int q_index = std::max(cell->active_fe_index(),cell->neighbor(face_number)->active_fe_index());
				hp_fe_face_values.reinit(cell, face_number,q_index);
				hp_neighbor_subface_values.reinit(cell->neighbor(face_number),cell->neighbor_of_coarser_neighbor(face_number).first, cell->neighbor_of_coarser_neighbor(face_number).second,q_index);

				const FEFaceValues<dim> &fe_face_values  = hp_fe_face_values.get_present_fe_values ();
				const FESubfaceValues<dim> &neighbor_subface_values = hp_neighbor_subface_values.get_present_fe_values ();


				const std::vector<double>& JxW_values = fe_face_values.get_JxW_values ();

				const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

				gradients.resize(n_face_q_points);
				neighbor_gradients.resize(n_face_q_points);

				neighbor_subface_values[velocities].get_function_gradients(solution, neighbor_gradients);
				fe_face_values[velocities].get_function_gradients(solution, gradients);

				std::vector<Tensor<1,dim> > jump_per_face;
				jump_per_face.resize(n_face_q_points);
				double jump_val=0;
				for (unsigned int q=0;
						q<n_face_q_points; ++q)
				{
					for (unsigned int i=0; i<2; ++i){
						for (unsigned int j=0; j<2; ++j){
							jump_per_face[q][i] += (gradients[q][i][j]- neighbor_gradients[q][i][j])*(fe_face_values.normal_vector(q)[j]);
						}//j
					}// i
					jump_val += contract(jump_per_face[q],jump_per_face[q])*JxW_values[q];
				}//q_per_face
				term3 +=(cell->face(face_number)->diameter())/(2.0 * cell->get_fe().degree)*jump_val;

			} // else if coarse neighbor

		Jump_est_per_cell(cell_index) = term3;
		est_per_cell(cell_index)=sqrt(Jump_est_per_cell(cell_index)+res_est_per_cell(cell_index));
	}//cell
}// func.estimate ()

/*......................................................................................*/
// get_layers_of_patch_around_cell

template <int dim>
std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> StokesProblem <dim>::get_patch_around_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell)
{

 std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch;
 std::set<typename hp::DoFHandler<dim>::active_cell_iterator> cells_done;

 patch.push_back (cell);
 cells_done.insert(cell);
 //  i counter for the number of patch layers ... n_layers=1 here (1 level of patch around cell)
 for (unsigned int i=0; i<1; ++i)
   {
     const unsigned int patch_size = patch.size();
     for (unsigned int j=0; j<patch_size; ++j)
       {
         for (unsigned int face_number=0; face_number< GeometryInfo<dim>::faces_per_cell; ++face_number)
           {
             if (patch[j]->face(face_number)->at_boundary()==false)
               {
                 if (patch[j]->face(face_number)->has_children() == false)
                   {
                     typename hp::DoFHandler<dim>::active_cell_iterator celll = patch[j]->neighbor(face_number);
                     if (cells_done.count(celll)==0)
                       {
                         patch.push_back(celll);
                         cells_done.insert(celll);
                       }
                   }
                 else
                   for (unsigned int subface=0; subface< patch[j]->face(face_number)->n_children(); ++subface)
                     {
                       typename hp::DoFHandler<dim>::active_cell_iterator child_cell = patch[j]->neighbor_child_on_subface (face_number, subface);
                       if (cells_done.count(child_cell)==0)
                         {
                           patch.push_back(child_cell);
                           cells_done.insert(child_cell);
                         }
                     }// for subface
               } // if at_boundary()==false
           }// face_number
       } // for j
   } // for i
 return patch;
}

/*......................................................................................*/
//get_cells_at_coarsest_common_level

template <int dim>
std::vector<typename hp::DoFHandler<dim>::cell_iterator> StokesProblem <dim>::get_cells_at_coarsest_common_level (const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch)
{
 Assert (patch.size() > 0, ExcMessage("vector containing patch cells should not be an empty vector!"));
 //Assert (patch.size() > 0, ExcInternalError());
 unsigned int min_level = static_cast<unsigned int> (patch[0]->level());
 unsigned int max_level = static_cast<unsigned int> (patch[0]->level());
 for (unsigned int i=0; i<patch.size();++i)
   {
     min_level = std::min (min_level, static_cast<unsigned int> (patch[i]->level()) );
     max_level = std::max (max_level, static_cast<unsigned int> (patch[i]->level()) );
   }

 std::set<typename hp::DoFHandler<dim>::cell_iterator>  uniform_cells;

 typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::const_iterator  patch_c;

 for (patch_c=patch.begin(); patch_c!=patch.end () ; ++patch_c){
	  if (static_cast<unsigned int>((*patch_c)->level()) == min_level)
		  uniform_cells.insert (*patch_c);
	  else
	  {
		  typename hp::DoFHandler<dim>::cell_iterator parent = *patch_c;

		  while (static_cast<unsigned int> (parent->level()) > min_level)
			  parent = parent-> parent();
		  uniform_cells.insert (parent);
	  }
 }

 return std::vector<typename hp::DoFHandler<dim>::cell_iterator> (uniform_cells.begin(), uniform_cells.end());

}	//get_cells_at_coarsest_common_level

/*......................................................................................*/
//build_triangulation_from_patch
template <int dim>
void StokesProblem <dim>::build_triangulation_from_patch (const std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>  &patch,
		Triangulation<dim> &local_triangulation, unsigned int &level_h_refine, unsigned int &level_p_refine, std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> & patch_to_global_tria_map )
		{
	std::vector<typename hp::DoFHandler<dim>::cell_iterator> uniform_cells = get_cells_at_coarsest_common_level (patch);

	level_h_refine=static_cast<unsigned int> (patch[0]->level());
	level_p_refine=static_cast<unsigned int> (patch[0]->active_fe_index());

	local_triangulation.clear();
	std::vector<Point<dim> > vertices;
	const unsigned int n_uniform_cells=uniform_cells.size();
	std::vector<CellData<dim> > cells(n_uniform_cells);
	unsigned int k=0;// for enumerating cells
	unsigned int i=0;// for enumerating vertices

	typename std::vector<typename hp::DoFHandler<dim>::cell_iterator>::const_iterator uniform_c;

	for (uniform_c=uniform_cells.begin(); uniform_c!=uniform_cells.end(); ++uniform_c)
	{
		bool repeat_vertex;
		for (unsigned int j=0;  j< GeometryInfo<dim>::vertices_per_cell; ++j)
		{
			Point<dim> position=(*uniform_c)->vertex (j);
			repeat_vertex=false;

			for (unsigned int m=0; m<i; ++m)
			{

				if (position == vertices[m]){
					repeat_vertex=true;
					cells[k].vertices[j]=m ;
					break;
				}
			}
			if (repeat_vertex==false)
			{
				vertices.push_back(position);
				cells[k].vertices[j]=i;
				i=i+1;
			}

		}//for vertices_per_cell
		k=k+1;
	}//uniform_c

	local_triangulation.create_triangulation(vertices,cells,SubCellData());

	//std::cout<< "local_triangulation.create_triangulation   . . . " << std::endl;
	Assert (local_triangulation.n_active_cells() == uniform_cells.size(), ExcInternalError());

	local_triangulation.clear_user_flags ();
	unsigned int index=0;
	std::map<typename Triangulation<dim>::cell_iterator, typename hp::DoFHandler<dim>::cell_iterator> patch_to_global_tria_map_tmp;
	for (typename Triangulation<dim>::cell_iterator coarse_cell_t = local_triangulation.begin(); coarse_cell_t != local_triangulation.end(); ++coarse_cell_t, ++index)

	{
		patch_to_global_tria_map_tmp.insert (std::make_pair(coarse_cell_t, uniform_cells[index]));
		AssertThrow (( std::fabs(coarse_cell_t->center()(0) - uniform_cells[index]->center()(0))<1e-16 && std::fabs(coarse_cell_t->center()(1) - uniform_cells[index]->center()(1)) <1e-16)  , ExcInternalError());
	}
	//std::cout<< "Before refinement_necessary " << std::endl;
	bool refinement_necessary;
	do
	{
		refinement_necessary = false;
		for (typename Triangulation<dim>::active_cell_iterator cell_tt = local_triangulation.begin_active(); cell_tt != local_triangulation.end(); ++cell_tt)

		{
			// if (patch_to_global_tria_map.count(cell_tt)==0)
			if (patch_to_global_tria_map_tmp[cell_tt]->has_children())
			{
				// std::cout<< "cell_tt -> set_refine_flag()" << std::endl;
				cell_tt -> set_refine_flag();
				refinement_necessary = true;

			}
			else for (unsigned int i=0; i<patch.size(); ++i){
				if (patch_to_global_tria_map_tmp[cell_tt]==patch[i]){
					// this flag shows that this cell is on the patch
					//  std::cout<< " shows that   patch_to_global_tria_map_tmp[cell_tt]   is on the patch " <<std::endl;
					cell_tt->set_user_flag();
					break;}
			}
		}//for tria_patch.begin...

		if (refinement_necessary)
		{
			local_triangulation.execute_coarsening_and_refinement ();

			for (typename Triangulation<dim>::cell_iterator cell_ttt = local_triangulation.begin(); cell_ttt != local_triangulation.end(); ++cell_ttt)
			{
				if (cell_ttt-> has_children())
				{
					//std::cout<< "cell_ttt-> has_children  : " << cell_ttt-> has_children() << std::endl;
					//std::cout<< " cell_ttt ->n_children()  " << cell_ttt ->n_children() << std::endl;


					// Note: Since the cell got children, then it should not be in the map anymore...children may be added into the map, instead


					// these children may not yet be in the map
					for (unsigned int c=0; c< cell_ttt ->n_children(); ++c)
					{

						if (patch_to_global_tria_map_tmp.find(cell_ttt->child(c)) == patch_to_global_tria_map_tmp.end())
						{
							//std::cout<< " making pair ... "<<std::endl;
							patch_to_global_tria_map_tmp.insert (std::make_pair(cell_ttt ->child(c), patch_to_global_tria_map_tmp[cell_ttt]->child(c)));
							//std::cout<< "before assert ...       pairing the cells with the same centers    ... " << std::endl;
							AssertThrow( (std::fabs (cell_ttt ->child(c)->center()(0)- patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(0)) < 1e-16 && std::fabs (cell_ttt ->child(c)->center()(1)-patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(1)) < 1e-16) , ExcInternalError());
							// std::cout<< "after assert ... " << std::endl;
						}
					}
					//patch_to_global_tria_map_tmp.erase(cell_ttt);

				}
			}
			// how do we know that here they pair the cells with the same center?
		}
		//  std::cout<<"the condition of DO loop is still satisfied... refinement_necessary " << refinement_necessary <<std::endl;

	}

	while (refinement_necessary);

	//std::cout<< "After refinement_necessary " << std::endl;
	typename std::map<typename Triangulation<dim>::cell_iterator, typename hp::DoFHandler<dim>::cell_iterator>::iterator map_tmp_it = patch_to_global_tria_map_tmp.begin(),
			map_tmp_end = patch_to_global_tria_map_tmp.end();

	for (; map_tmp_it!=map_tmp_end; ++map_tmp_it)
		patch_to_global_tria_map[map_tmp_it->first] = map_tmp_it->second;

		} // build_triangulation
/*......................................................................................................................................*/
//set_active_fe_indices for cells on and out of each patch


// mark cells in the "local_triangulation" that exist in the patch: go
// through all cells in the "local_triangulation" and see whether the
// corresponding cell in the global triangulation is part of
// the 'patch' list of cells

template <int dim>
void StokesProblem <dim>::set_active_fe_indices (hp::DoFHandler<dim> &local_dof_handler,
std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> & patch_to_global_tria_map)
{
 typename hp::DoFHandler<dim>::active_cell_iterator patch_c= local_dof_handler.begin_active(), end_patch_c = local_dof_handler.end();
 for (; patch_c!=end_patch_c; ++patch_c){

     if (patch_c->user_flag_set()==true)
       {

          typename hp::DoFHandler<dim>::active_cell_iterator global_cell = patch_to_global_tria_map[patch_c];


          patch_c->set_active_fe_index (  global_cell->active_fe_index() );

         //patch_c->set_active_fe_index (0);
       }
     else if (patch_c->user_flag_set()==false)
       {
       //typename hp::DoFHandler<dim>::active_cell_iterator global_cell = patch_to_global_tria_map[patch_c];

   	  // which assigns FE_Nothing for the cells out of patch
         patch_c->set_active_fe_index (max_degree);
       }
     else
       Assert (false, ExcNotImplemented());
 }

}

//.................................................................................................................................
template <int dim>
void StokesProblem <dim>:: patch_output (unsigned int patch_number, const unsigned int cycle, hp::DoFHandler<dim> &local_dof_handler, BlockVector<double> &local_solu)
{

 std::vector<std::string> solution_names (dim, "patch_velocity");
 solution_names.push_back ("patch_pressure");

 std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
 (dim, DataComponentInterpretation::component_is_part_of_vector);
 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

 DataOut<dim,hp::DoFHandler<dim> > patch_data_out;

 patch_data_out.attach_dof_handler (local_dof_handler);


 patch_data_out.add_data_vector (local_solu, solution_names, DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, data_component_interpretation);
 patch_data_out.build_patches ();

 std::string filename = "patch_solution-" +
     Utilities::int_to_string (cycle, 2) +
     +"-"+Utilities::int_to_string (patch_number, 2) +".vtu";
 std::ofstream output (filename.c_str());
 patch_data_out.write_vtu (output);
}


/*......................................................................................................................................*/
//Compute h_convergence_estimator   &   h_workload_number  for each patch around cell

template <int dim>
void StokesProblem <dim>::h_patch_conv_load_no ( const unsigned int cycle , double &h_convergence_est_per_cell,
		unsigned int &h_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		 unsigned int & patch_number)
		{

	Triangulation<dim> local_triangulation;
	unsigned int level_h_refine;
	unsigned int level_p_refine;
	std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> patch_to_global_tria_map;

	std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell(cell);
	build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine, patch_to_global_tria_map);
	hp::DoFHandler<dim> local_dof_handler(local_triangulation);

	//std::cout<< "patch_number is :" << " " << patch_number << std:: endl ;

	set_active_fe_indices (local_dof_handler,patch_to_global_tria_map );
	local_dof_handler.distribute_dofs (fe_collection);

 DoFRenumbering::Cuthill_McKee (local_dof_handler);
 //std::cout<< "number of dofs on the current patch " << local_dof_handler.n_dofs() << std::endl ;

 h_convergence_est_per_cell=0.;
 double h_solu_norm_per_patch=0.;

 ConstraintMatrix constraints_patch;
 BlockSparsityPattern sparsity_pattern_patch;

 std::vector<unsigned int> block_component_patch (dim+1, 0);
 block_component_patch[dim]=1;
 DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

 std::vector<types::global_dof_index> dofs_per_block_patch (2);
 DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);

//......................................................................................................
 BlockVector<double> local_solu (dofs_per_block_patch);

 // Here we are trying to project the values of the global vector "solution" into vector "local_solu" which is
 //solution over patch cells corresponding to cell "K".

 typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
 for (; patch_cl !=end_patch_cl; ++patch_cl)
 {

	  const unsigned int   dofs_per_cl = patch_cl->get_fe().dofs_per_cell;
	  // we check if the corresponding finite element for this cell is not 'FE_Nothing!' and it takes usual finite element.
	  if (dofs_per_cl!=0)
	  {

		  Vector<double> local_solution_values(dofs_per_cl);
		  typename hp::DoFHandler<dim>::active_cell_iterator global_cell = patch_to_global_tria_map[patch_cl];

		  global_cell-> get_dof_values (solution,local_solution_values);
		  patch_cl->set_dof_values (local_solution_values, local_solu);
	  }
 }


 /*
 std::cout<<"the local projected solution is: "<< std::endl;
 local_solu.print(std::cout);
 std::cout<<"the global solution is: "<< std::endl;
 solution.print(std::cout);

  */

//......................  Solution Transfer .......................  h_refinement of patch cells ................................. //

 bool need_to_refine = false;

 typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
 for (; patch_cc!=end_patch_cc; ++patch_cc)
 {
	  if (static_cast<unsigned int> (patch_cc->level()) <  (level_h_refine+1) )
	  {
		  need_to_refine = true;
		  patch_cc->set_refine_flag();
	  }
 }//patch_cc

 if (need_to_refine == true)
 {
	  // user flags will be overwritten by
	  // execute_coarsening_and_refinement. save their values into
	  // the material_id, since that one not only survives
	  // refinement but is also inherited to the children
	  for (typename hp::DoFHandler<dim>::cell_iterator cell_ttt = local_dof_handler.begin(); cell_ttt != local_dof_handler.end(); ++cell_ttt)
	    if (cell_ttt->user_flag_set())
	      cell_ttt->set_material_id (1);
	    else
	      cell_ttt->set_material_id (0);

	  local_triangulation.prepare_coarsening_and_refinement ();
	  SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim>> solution_transfer(local_dof_handler);
	  solution_transfer.prepare_for_pure_refinement();


         //local_triangulation.execute_refinement ();
	  local_triangulation.execute_coarsening_and_refinement ();

	  // get user flags back out of the material_id field
	  for (typename hp::DoFHandler<dim>::cell_iterator cell_ttt = local_dof_handler.begin(); cell_ttt != local_dof_handler.end(); ++cell_ttt)
	    if (cell_ttt->material_id() == 1)
	      cell_ttt->set_user_flag();
	    else
	      cell_ttt->clear_user_flag();
// I believe here I need this function to be un comment! despite of h-refinement

	 // set_active_fe_indices (local_dof_handler,patch_to_global_tria_map );
	  local_dof_handler.distribute_dofs (fe_collection);

	 // std::vector<unsigned int> block_component_patch (dim+1, 0);
	 // block_component_patch[dim]=1;
	  DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

	  // take a copy of the tmp vector
	  // BlockVector<double>   local_solu(temp);

	  DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
	  // resize the vector temp to the correct size
	  BlockVector<double> temp (dofs_per_block_patch);
	  solution_transfer.refine_interpolate(local_solu , temp);
	  local_solu = temp;
 }


//......................  setup_h_patch_system and  patch_rhs .............. //

 unsigned int local_system_size = local_dof_handler. n_dofs();
 h_workload_num = local_dof_handler. n_dofs();

 //patch_output (patch_number , cycle, local_dof_handler, local_solu);

 {
   constraints_patch.clear ();
   FEValuesExtractors::Vector velocities(0);
   DoFTools::make_hanging_node_constraints(local_dof_handler, constraints_patch);
/*
   VectorTools::interpolate_boundary_values (local_dof_handler,
                                              0,
                                              ZeroFunction<dim>(dim+1),
                                              constraints_patch,
                                              fe_collection.component_mask(velocities));

*/

   //......................... Zero_Bdry_Condition_on_Patch .......................................//
   {
   	typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
   	for (; patch_cl !=end_patch_cl; ++patch_cl)
   	{
   		std::vector<types::global_dof_index> local_face_dof_indices ((patch_cl->get_fe()).dofs_per_face);
   		if (patch_cl->user_flag_set() == true)
   		{
   			for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
   			{
   				bool face_is_on_patch_Bdry = false;
   				if ((patch_cl->face(f)->at_boundary()) ||  ((patch_cl->neighbor(f)->has_children() == false)
   						&&
   						(patch_cl->neighbor(f)->user_flag_set() == false)))
   					face_is_on_patch_Bdry = true;
   				else if ((patch_cl->face(f)->at_boundary()) || (patch_cl->neighbor(f)->has_children() == true))
   				{
   					for (unsigned int sf=0; sf< patch_cl->face(f)->n_children(); ++sf)
   						if (patch_cl->neighbor_child_on_subface (f, sf) -> user_flag_set() == false)
   						{
   							face_is_on_patch_Bdry = true;
   							break;
   						}
   				}// else if

   				if (face_is_on_patch_Bdry)
   				{
   					patch_cl->face(f)->get_dof_indices (local_face_dof_indices, patch_cl->active_fe_index());
   					for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
   						// system_to_component_index: "Compute vector component and index of this shape function within the shape functions
   						// corresponding to this component from the index of a shape function within this finite element"
   						if ((patch_cl->get_fe()).face_system_to_component_index(i).first < dim)
						  {
   							constraints_patch.add_line (local_face_dof_indices[i]);
						  }
   				}
   			}// face f
   		}// if user_flag_true

   	}// patch_cl
   }

 }
 constraints_patch.close();

 DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
 const unsigned int n_u=dofs_per_block_patch[0], n_p=dofs_per_block_patch[1];

 {
   BlockCompressedSetSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
   DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
   sparsity_pattern_patch.copy_from(csp);
 }

 BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
 BlockVector<double> patch_solution (dofs_per_block_patch);
 BlockVector<double> patch_rhs (dofs_per_block_patch);

 // .........................................  assemble  patch_system  and patch_rhs .............................. //

 hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

 FullMatrix<double> local_matrix_patch;
 Vector<double> local_rhs_patch;
 Vector<double> local_rhs1;
 Vector<double> local_rhs2;
 std::vector<types::global_dof_index> local_dof_indices;

 std::vector<Vector<double> >  rhs_values;
 const RightHandSide<dim> rhs_function;

 const FEValuesExtractors::Vector velocities (0);
 const FEValuesExtractors::Scalar pressure (dim);

 //std::vector<SymmetricTensor<2,dim> > symgrad_phi_u;
 std::vector<Tensor<2,dim> > grad_phi_u;
 std::vector<double> div_phi_u;
 std::vector<Tensor<1,dim> > phi_u;
 std::vector<double> phi_p;

 std::vector<Tensor<1,dim> > gradients_p;
 std::vector<double> divergences;
 std::vector<Tensor<1,dim> > laplacians;

 std::vector<double> values;
 std::vector<Tensor<2,dim> > gradients;


 typename hp::DoFHandler<dim>::active_cell_iterator patch_cll= local_dof_handler.begin_active(), end_patch_cll = local_dof_handler.end();
 for (; patch_cll!=end_patch_cll; ++patch_cll){

     const unsigned int   dofs_per_cell = patch_cll->get_fe().dofs_per_cell;
     if (dofs_per_cell!=0) {
         local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
         local_rhs_patch.reinit (dofs_per_cell);
         local_rhs1.reinit (dofs_per_cell);
         local_rhs2.reinit (dofs_per_cell);

         hp_fe_values.reinit (patch_cll);
         const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
         const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
         const unsigned int n_q_points = fe_values.n_quadrature_points;

         rhs_values.resize(n_q_points, Vector<double>(dim+1));
         rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

         grad_phi_u.resize(dofs_per_cell);
         //symgrad_phi_u.resize(dofs_per_cell);
         div_phi_u.resize(dofs_per_cell);
         phi_u.resize (dofs_per_cell);
         phi_p.resize(dofs_per_cell);

         divergences.resize(n_q_points);
         gradients_p.resize(n_q_points) ;
         laplacians.resize(n_q_points);

         fe_values[pressure].get_function_gradients(local_solu, gradients_p);
         fe_values[velocities].get_function_divergences(local_solu, divergences);
         fe_values[velocities].get_function_laplacians(local_solu, laplacians);


         for (unsigned int q=0; q<n_q_points; ++q)
           {
             for (unsigned int k=0; k<dofs_per_cell; ++k)
               {
                 grad_phi_u[k] = fe_values[velocities].gradient (k, q);
                 //symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                 phi_u[k] = fe_values[velocities].value (k, q);
                 phi_p[k] = fe_values[pressure].value (k, q);
               }


             for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
           	  for (unsigned int j=0; j<dofs_per_cell; ++j)
           	  {

           	  local_matrix_patch(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])  + (phi_p[i] * phi_p[j]))* JxW_values[q];

           		//  local_matrix_patch(i,j) += ((symgrad_phi_u[i] * symgrad_phi_u[j])
           				 // + (phi_p[i] * phi_p[j]))
           				//  * JxW_values[q];

           	  } // end of loop over 'j'


                 local_rhs1(i)+= ((rhs_values[q](0)+ laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+
					  (rhs_values[q](1)+ laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];
                 local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];
                 local_rhs_patch(i)= local_rhs1(i)+local_rhs2(i);
               }// i
           }//q

         //for (unsigned int k=0; k<phi_p.size(); ++k){
       	//  std::cout<< "phi_p [" << k  << "] =" <<phi_p[k]<< std::endl;
       	//  std::cout<< " symgrad_phi_u[ " << k << "]= " <<  symgrad_phi_u[k] << std::endl ;
        // }

         local_dof_indices.resize (dofs_per_cell);
         patch_cll->get_dof_indices (local_dof_indices);
         constraints_patch.distribute_local_to_global (local_matrix_patch, local_rhs_patch, local_dof_indices, patch_system, patch_rhs);
     }
 }// for patch_cll

/*
 for (unsigned int l=0; l< dofs_per_block_patch[0]+dofs_per_block_patch[1]; ++l)
	  for (unsigned int m=0; m< dofs_per_block_patch[0]+dofs_per_block_patch[1] ; ++m)
		  std::cout<<  "patch system : [ " << l<< " , "  << m << "]= "  << patch_system.el (l,m) << "............." <<std::endl;
*/

//  ................................................  symmetry check  ...........................................................
 {
   BlockVector<double> temp1 (dofs_per_block_patch);
   BlockVector<double> temp2 (dofs_per_block_patch);
   BlockVector<double> unit_vector (dofs_per_block_patch);
   for (unsigned int i=0; i<local_dof_handler.n_dofs(); ++i)
     {
       for (unsigned int j=0; j<local_dof_handler.n_dofs(); ++j)
         {
           if (i==j)
             unit_vector[j]=1;
           else
             unit_vector[j]=0;
         } // for j

       patch_system.vmult (temp1,unit_vector );
       patch_system.Tvmult (temp2, unit_vector);

       for (unsigned int k=0; k<local_dof_handler.n_dofs(); ++k)
         {
           if(std::fabs(temp1[k]-temp2[k]) > 1e-8)
             //if (temp1[k]!=temp2[k])
             {
               std::cout<< "symmetricity check  "<< std::endl;
               std::cout<< "comparison between row and column number: " << i << std::endl;
               std::cout<< "vector temp1 [ "<< k << " ]=" << temp1[k] << std::endl;
               std::cout<< "vector temp2 [ "<< k << " ]=" << temp2[k] << std::endl;
             }
         }// for k

     }// for i
 }

 // patch_system.block(0,1).print(std::cout);


/*
// direct solver
 SparseDirectUMFPACK A_inverse_stiffness;
 A_inverse_stiffness.initialize (patch_system.block(0,0),
		 SparseDirectUMFPACK::AdditionalData());
 A_inverse_stiffness.vmult (patch_solution.block(0), patch_rhs.block(0));
 //constraints_patch.distribute (patch_solution.block(0));

 SparseDirectUMFPACK A_inverse_mass;
 A_inverse_mass.initialize (patch_system.block(1,1),
		 SparseDirectUMFPACK::AdditionalData());
 A_inverse_mass.vmult (patch_solution.block(1), patch_rhs.block(1));
 //constraints_patch.distribute (patch_solution.block(1));

 constraints_patch.distribute (patch_solution);
*/




 // iterative solver
 SolverControl           solver_control_stiffness (patch_rhs.block(0).size(),1e-8*patch_rhs.block(0).l2_norm());
 SolverCG<>              cg_stiff (solver_control_stiffness);

 PreconditionSSOR<> preconditioner_stiffness;
 preconditioner_stiffness.initialize(patch_system.block(0,0), 1.2);
// std::cout<<"before cg_stiff"<<std::endl;
 cg_stiff.solve (patch_system.block(0,0), patch_solution.block(0), patch_rhs.block(0),
		  preconditioner_stiffness);
// std::cout<<"after cg_stiff"<<std::endl;


 SolverControl           solver_control_mass (patch_rhs.block(1).size(),1e-8*patch_rhs.block(1).l2_norm());
 SolverCG<>              cg_mass (solver_control_mass);

 PreconditionSSOR<> preconditioner_mass;
 preconditioner_mass.initialize(patch_system.block(1,1), 1.2);
// std::cout<<"before cg_mass"<<std::endl;
 cg_mass.solve (patch_system.block(1,1), patch_solution.block(1), patch_rhs.block(1),
		  preconditioner_mass);
// std::cout<<"after cg_mass"<<std::endl;
 constraints_patch.distribute (patch_solution);




//.......................................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//

 double pressure_val=0;
 double grad_u_val=0;
 typename hp::DoFHandler<dim>::active_cell_iterator patch_cel= local_dof_handler.begin_active(), end_patch_cel = local_dof_handler.end();
 for (; patch_cel!=end_patch_cel; ++patch_cel)
 {
	  const unsigned int   dofs_per_cel = patch_cel->get_fe().dofs_per_cell;
	  if (dofs_per_cel!=0)
	  {
		  hp_fe_values.reinit (patch_cel);
		  const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
		  const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
		  const unsigned int n_q_points = fe_values.n_quadrature_points;

		  gradients.resize(n_q_points);
		  values.resize(n_q_points);

		  fe_values[velocities].get_function_gradients(patch_solution, gradients);
		  fe_values[pressure].get_function_values(patch_solution, values);


		  for (unsigned int q=0; q<n_q_points; ++q)
		  {
			  pressure_val +=values[q]*values[q]* JxW_values[q];

			  for (unsigned int i=0; i<dim; ++i)

				  grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];
			  //grad_u_val +=double_contract(gradients[q],gradients[q])* JxW_values[q];
		  } // q
		  h_solu_norm_per_patch +=pressure_val + grad_u_val;
	  }

 }// cells on patch
 h_convergence_est_per_cell =sqrt(h_solu_norm_per_patch);
// std::cout<< "Improvement in the solutions after refinement in h- : "<< h_convergence_est_per_cell << std::endl;

		}  // function h_patch_con_load_num


/*......................................................................................................................................*/



//Compute p_convergence_estimator   &   p_workload_number  for each patch around cell

template <int dim>
void StokesProblem <dim>::p_patch_conv_load_no ( const unsigned int cycle , double &p_convergence_est_per_cell,
		unsigned int &p_workload_num, const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		unsigned int & patch_number)

		{

	Triangulation<dim> local_triangulation;
	unsigned int level_h_refine;
	unsigned int level_p_refine;
	std::map<typename Triangulation<dim>::active_cell_iterator, typename hp::DoFHandler<dim>::active_cell_iterator> patch_to_global_tria_map;

	std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch = get_patch_around_cell(cell);

	build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine, patch_to_global_tria_map);
	hp::DoFHandler<dim> local_dof_handler(local_triangulation);

	//std::cout<< "patch_number is :" << " " << patch_number << std:: endl ;

	set_active_fe_indices (local_dof_handler,patch_to_global_tria_map);
	local_dof_handler.distribute_dofs (fe_collection);

	DoFRenumbering::Cuthill_McKee (local_dof_handler);
	//std::cout<< "number of dofs on the current patch " << local_dof_handler.n_dofs() << std::endl ;

	p_convergence_est_per_cell=0.;
	double p_solu_norm_per_patch=0.;

	ConstraintMatrix constraints_patch;
	BlockSparsityPattern sparsity_pattern_patch;

	std::vector<unsigned int> block_component_patch (dim+1, 0);
	block_component_patch[dim]=1;
	DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

	std::vector<types::global_dof_index> dofs_per_block_patch (2);
	DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);

	BlockVector<double> local_solu (dofs_per_block_patch);

	// Here we are trying to project the values of the global vector "solution" into vector "local_solu" which is
	//solution over patch cells corresponding to cell "K".

	typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
	for (; patch_cl !=end_patch_cl; ++patch_cl)
	{

		const unsigned int   dofs_per_cl = patch_cl->get_fe().dofs_per_cell;
		// we check if the corresponding finite element for this cell is not 'FE_Nothing!' and it takes usual finite element.
		if (dofs_per_cl!=0)
		{

			Vector<double> local_solution_values(dofs_per_cl);
			typename hp::DoFHandler<dim>::active_cell_iterator global_cl = patch_to_global_tria_map[patch_cl];

			global_cl-> get_dof_values (solution,local_solution_values);
			patch_cl->set_dof_values (local_solution_values, local_solu);
		}
	}

	//......................  Solution Transfer .......................  p_refinement of patch cells ................................. //
	bool need_to_p_refine = false;


	typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
	for (; patch_cc!=end_patch_cc; ++patch_cc)
	{

	typename hp::DoFHandler<dim>::active_cell_iterator global_cell = patch_to_global_tria_map[patch_cc];


		if ( ( global_cell->active_fe_index() +1) < (fe_collection.size()-1) && global_cell->active_fe_index() < (level_p_refine+1))
			need_to_p_refine = true;

	}//for patch_c


	if (need_to_p_refine == true)
	{
		// user flags will be overwritten by
		// execute_coarsening_and_refinement. save their values into
		// the material_id, since that one not only survives
		// refinement but is also inherited to the children
		for (typename hp::DoFHandler<dim>::cell_iterator cell_ttt = local_dof_handler.begin(); cell_ttt != local_dof_handler.end(); ++cell_ttt)
			if (cell_ttt->user_flag_set())
				cell_ttt->set_material_id (1);
			else
				cell_ttt->set_material_id (0);

		 local_triangulation.prepare_coarsening_and_refinement ();


		SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim>> solution_transfer(local_dof_handler);
		solution_transfer.prepare_for_pure_refinement();
		// solution_transfer.prepare_for_coarsening_and_refinement(local_solu);


		typename hp::DoFHandler<dim>::active_cell_iterator patch_cc= local_dof_handler.begin_active(), end_patch_cc = local_dof_handler.end();
		for (; patch_cc!=end_patch_cc; ++patch_cc)
		{

		        // since here the fe_collection.size()=7 (i.e., 6 indices in total), and we also know the last index hold for the fe_nothing FE, therefore we will
		        // set_active_fe_index for the cells up to index 4. (the reason is that for example we cannot exceed the polynomial degree for the index 5...it is
		        // already the last index before the fe_nothing (index=6))
		        // It is also worth to see how they did p-refinement in step-27.
		        // if  (cell->active_fe_index()+1 < fe_collection.size()))



		       typename hp::DoFHandler<dim>::active_cell_iterator global_cc = patch_to_global_tria_map[patch_cc];

                       patch_cc->set_active_fe_index (global_cc->active_fe_index() +1 );

		}

		// get user flags back out of the material_id field
		for (typename hp::DoFHandler<dim>::cell_iterator cell_ttt = local_dof_handler.begin(); cell_ttt != local_dof_handler.end(); ++cell_ttt)
			if (cell_ttt->material_id() == 1)
				cell_ttt->set_user_flag();
			else
				cell_ttt->clear_user_flag();


               set_active_fe_indices (local_dof_handler,patch_to_global_tria_map );
	        local_dof_handler.distribute_dofs (fe_collection);

		DoFRenumbering::component_wise (local_dof_handler, block_component_patch);

		DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
		// resize the vector temp to the correct size
		BlockVector<double> temp (dofs_per_block_patch);
               solution_transfer.refine_interpolate(local_solu , temp);
		//solution_transfer.interpolate(local_solu , temp);
		local_solu = temp;

		//  we use  the clear() function, for deleting all stored data in SolutionTransfer and reinitializing it .
		//solution_transfer.clear();

	}

	//......................  setup_p_patch_system and  patch_rhs .............. //

	unsigned int local_system_size = local_dof_handler. n_dofs();
	p_workload_num = local_dof_handler. n_dofs();

	//patch_output (patch_number , cycle, local_dof_handler, local_solu);

	{
		constraints_patch.clear ();

		FEValuesExtractors::Vector velocities(0);
		DoFTools::make_hanging_node_constraints(local_dof_handler, constraints_patch);

		//......................... Zero_Bdry_Condition_on_Patch .......................................//
		{
			typename hp::DoFHandler<dim>::active_cell_iterator patch_cl= local_dof_handler.begin_active(), end_patch_cl = local_dof_handler.end();
			for (; patch_cl !=end_patch_cl; ++patch_cl)
			{
				std::vector<types::global_dof_index> local_face_dof_indices ((patch_cl->get_fe()).dofs_per_face);
				if (patch_cl->user_flag_set() == true)
				{
					for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
					{
						bool face_is_on_patch_Bdry = false;
						if ((patch_cl->face(f)->at_boundary()) ||  ((patch_cl->neighbor(f)->has_children() == false)
								&&
								(patch_cl->neighbor(f)->user_flag_set() == false)))
							face_is_on_patch_Bdry = true;
						else if ((patch_cl->face(f)->at_boundary()) || (patch_cl->neighbor(f)->has_children() == true))
						{
							for (unsigned int sf=0; sf< patch_cl->face(f)->n_children(); ++sf)
								if (patch_cl->neighbor_child_on_subface (f, sf) -> user_flag_set() == false)
								{
									face_is_on_patch_Bdry = true;
									break;
								}
						}// else if

						if (face_is_on_patch_Bdry)
						{
							patch_cl->face(f)->get_dof_indices (local_face_dof_indices, patch_cl->active_fe_index());
							for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
								// system_to_component_index: "Compute vector component and index of this shape function within the shape functions
								// corresponding to this component from the index of a shape function within this finite element"
								if ((patch_cl->get_fe()).face_system_to_component_index(i).first < dim)
								{
									constraints_patch.add_line (local_face_dof_indices[i]);
								}
						}
					}// face f
				}// if user_flag_true

			}// patch_cl
		}


	}
	constraints_patch.close();


       //constraints_patch.distribute(local_solu );
	DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
	const unsigned int n_u=dofs_per_block_patch[0], n_p=dofs_per_block_patch[1];

	{
		BlockCompressedSetSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
		DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
		sparsity_pattern_patch.copy_from(csp);
	}

	BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
	BlockVector<double> patch_solution (dofs_per_block_patch);
	BlockVector<double> patch_rhs (dofs_per_block_patch);

	// .........................................  assemble  patch_system  and patch_rhs .............................. //

	hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_quadrature_points|update_JxW_values|update_gradients|update_hessians);

	FullMatrix<double> local_matrix_patch;
	Vector<double> local_rhs_patch;
	Vector<double> local_rhs1;
	Vector<double> local_rhs2;
	std::vector<types::global_dof_index> local_dof_indices;

	std::vector<Vector<double> >  rhs_values;
	const RightHandSide<dim> rhs_function;

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	//std::vector<SymmetricTensor<2,dim> > symgrad_phi_u;
	std::vector<Tensor<2,dim> > grad_phi_u;
	std::vector<double> div_phi_u;
	std::vector<Tensor<1,dim> > phi_u;
	std::vector<double> phi_p;

	std::vector<Tensor<1,dim> > gradients_p;
	std::vector<double> divergences;
	std::vector<Tensor<1,dim> > laplacians;

	std::vector<double> values;
	std::vector<Tensor<2,dim> > gradients;


	typename hp::DoFHandler<dim>::active_cell_iterator patch_cll= local_dof_handler.begin_active(), end_patch_cll = local_dof_handler.end();
	for (; patch_cll!=end_patch_cll; ++patch_cll){

		const unsigned int   dofs_per_cell = patch_cll->get_fe().dofs_per_cell;
		if (dofs_per_cell!=0) {
			local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
			local_rhs_patch.reinit (dofs_per_cell);
			local_rhs1.reinit (dofs_per_cell);
			local_rhs2.reinit (dofs_per_cell);

			hp_fe_values.reinit (patch_cll);
			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
			const unsigned int n_q_points = fe_values.n_quadrature_points;

			rhs_values.resize(n_q_points, Vector<double>(dim+1));
			rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

			grad_phi_u.resize(dofs_per_cell);
			//symgrad_phi_u.resize(dofs_per_cell);
			div_phi_u.resize(dofs_per_cell);
			phi_u.resize (dofs_per_cell);
			phi_p.resize(dofs_per_cell);

			divergences.resize(n_q_points);
			gradients_p.resize(n_q_points) ;
			laplacians.resize(n_q_points);

			fe_values[pressure].get_function_gradients(local_solu, gradients_p);
			fe_values[velocities].get_function_divergences(local_solu, divergences);
			fe_values[velocities].get_function_laplacians(local_solu, laplacians);


			for (unsigned int q=0; q<n_q_points; ++q)
			{
				for (unsigned int k=0; k<dofs_per_cell; ++k)
				{
					grad_phi_u[k] = fe_values[velocities].gradient (k, q);
					//symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
					phi_u[k] = fe_values[velocities].value (k, q);
					phi_p[k] = fe_values[pressure].value (k, q);
				}


				for (unsigned int i=0; i<dofs_per_cell; ++i)
				{
					for (unsigned int j=0; j<dofs_per_cell; ++j)
					{
						local_matrix_patch(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])  + (phi_p[i] * phi_p[j]))* JxW_values[q];

						//local_matrix_patch(i,j) += ((symgrad_phi_u[i] * symgrad_phi_u[j])
						//+ (phi_p[i] * phi_p[j]))
						//* JxW_values[q];
					} // end of loop over 'j'


					local_rhs1(i)+= ((rhs_values[q](0)+laplacians[q][0]-gradients_p[q][0])*(phi_u[i][0])+
							(rhs_values[q](1)+laplacians[q][1]-gradients_p[q][1])*(phi_u[i][1]))* JxW_values[q];
					local_rhs2(i)+= (phi_p[i]*divergences[q])*JxW_values[q];
					local_rhs_patch(i)= local_rhs1(i)+local_rhs2(i);
				}// i
			}//q

			local_dof_indices.resize (dofs_per_cell);
			patch_cll->get_dof_indices (local_dof_indices);
			constraints_patch.distribute_local_to_global (local_matrix_patch, local_rhs_patch, local_dof_indices, patch_system, patch_rhs);
		}
	}// for patch_cll


	/*
	                   for (unsigned int l=0; l< dofs_per_block_patch[0]+dofs_per_block_patch[1]; ++l)
	                 	  for (unsigned int m=0; m< dofs_per_block_patch[0]+dofs_per_block_patch[1] ; ++m)
	                 		  std::cout<<  "patch system : [ " << l<< " , "  << m << "]= "  << patch_system.el (l,m) << "............." <<std::endl;
	 */
	//  ................................................  symmetry check  ...........................................................//

	{
		BlockVector<double> temp1 (dofs_per_block_patch);
		BlockVector<double> temp2 (dofs_per_block_patch);
		BlockVector<double> unit_vector (dofs_per_block_patch);
		for (unsigned int i=0; i<local_dof_handler.n_dofs(); ++i)
		{
			for (unsigned int j=0; j<local_dof_handler.n_dofs(); ++j)
			{
				if (i==j)
					unit_vector[j]=1;
				else
					unit_vector[j]=0;
			} // for j

			patch_system.vmult (temp1,unit_vector );
			patch_system.Tvmult (temp2, unit_vector);

			for (unsigned int k=0; k<local_dof_handler.n_dofs(); ++k)
			{
				if(std::fabs(temp1[k]-temp2[k]) > 1e-8)
					//if (temp1[k]!=temp2[k])
				{
					std::cout<< "symmetricity check  "<< std::endl;
					std::cout<< "comparison between row and column number: " << i << std::endl;
					std::cout<< "vector temp1 [ "<< k << " ]=" << temp1[k] << std::endl;
					std::cout<< "vector temp2 [ "<< k << " ]=" << temp2[k] << std::endl;
				}
			}// for k

		}// for i
	}

	// patch_system.block(0,1).print(std::cout);


	/*
	// direct solver
	 SparseDirectUMFPACK A_inverse_stiffness;
	 A_inverse_stiffness.initialize (patch_system.block(0,0),
			 SparseDirectUMFPACK::AdditionalData());
	 A_inverse_stiffness.vmult (patch_solution.block(0), patch_rhs.block(0));
	 //constraints_patch.distribute (patch_solution.block(0));

	 SparseDirectUMFPACK A_inverse_mass;
	 A_inverse_mass.initialize (patch_system.block(1,1),
			 SparseDirectUMFPACK::AdditionalData());
	 A_inverse_mass.vmult (patch_solution.block(1), patch_rhs.block(1));
	 //constraints_patch.distribute (patch_solution.block(1));

	 constraints_patch.distribute (patch_solution);
*/


	// iterative solver
	SolverControl           solver_control_stiffness (patch_rhs.block(0).size(),1e-8*patch_rhs.block(0).l2_norm());
	SolverCG<>              cg_stiff (solver_control_stiffness);

	PreconditionSSOR<> preconditioner_stiffness;
	preconditioner_stiffness.initialize(patch_system.block(0,0), 1.2);
	// std::cout<<"before cg_stiff"<<std::endl;
	cg_stiff.solve (patch_system.block(0,0), patch_solution.block(0), patch_rhs.block(0),
			preconditioner_stiffness);
	// std::cout<<"after cg_stiff"<<std::endl;


	SolverControl           solver_control_mass (patch_rhs.block(1).size(),1e-8*patch_rhs.block(1).l2_norm());
	SolverCG<>              cg_mass (solver_control_mass);

	PreconditionSSOR<> preconditioner_mass;
	preconditioner_mass.initialize(patch_system.block(1,1), 1.2);
	// std::cout<<"before cg_mass"<<std::endl;
	cg_mass.solve (patch_system.block(1,1), patch_solution.block(1), patch_rhs.block(1),
			preconditioner_mass);
	// std::cout<<"after cg_mass"<<std::endl;
	constraints_patch.distribute (patch_solution);


	//.......................................  get the L2 norm of the gradient of velocity solution and pressure value  .....................//

	double pressure_val=0;
	double grad_u_val=0;
	typename hp::DoFHandler<dim>::active_cell_iterator patch_cel= local_dof_handler.begin_active(), end_patch_cel = local_dof_handler.end();
	for (; patch_cel!=end_patch_cel; ++patch_cel)
	{
		const unsigned int   dofs_per_cel = patch_cel->get_fe().dofs_per_cell;
		if (dofs_per_cel!=0)
		{
			hp_fe_values.reinit (patch_cel);
			const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
			const std::vector<double>& JxW_values = fe_values.get_JxW_values ();
			const unsigned int n_q_points = fe_values.n_quadrature_points;

			gradients.resize(n_q_points);
			values.resize(n_q_points);

			fe_values[velocities].get_function_gradients(patch_solution, gradients);
			fe_values[pressure].get_function_values(patch_solution, values);


			for (unsigned int q=0; q<n_q_points; ++q)
			{
				pressure_val +=values[q]*values[q]* JxW_values[q];

				for (unsigned int i=0; i<dim; ++i)

					grad_u_val += contract(gradients[q][i],gradients[q][i])* JxW_values[q];
				//grad_u_val +=double_contract(gradients[q],gradients[q])* JxW_values[q];
			} // q
			p_solu_norm_per_patch +=pressure_val + grad_u_val;
		}

	}// cells on patch
	p_convergence_est_per_cell =sqrt(p_solu_norm_per_patch);
	//std::cout<< "Improvement in the solutions after refinement in h- : "<< p_convergence_est_per_cell << std::endl;
		}

/*..............................................   marking_cells   .....................................................*/
template <int dim>
void StokesProblem <dim>:: marking_cells (const unsigned int cycle, Vector<float> & marked_cells, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set,
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map, Vector<double> & h_Conv_Est, Vector<double> &p_Conv_Est, Vector<double> &hp_Conv_Est)
		{
	// based on Dorfler's paper, this \theta should be chosen as a number close to minimum of convergence estimator! \theta \in (0, min(convergence estimator)_K)

	std::vector<std::pair<double, typename hp::DoFHandler<dim>::active_cell_iterator> > to_be_sorted;
	std::vector<std::pair<typename hp::DoFHandler<dim>::active_cell_iterator,bool > > to_be_sorted_with_refine_info;

	Vector<double> est_per_cell (triangulation.n_active_cells());
	estimate(est_per_cell);

	//  this vector "convergence_est_per_cell" will be finalized after checking out which h- or p- refinement are going to be chosen for each cell
	Vector<double> convergence_est_per_cell (triangulation.n_active_cells());
	Vector<double> hp_Conv_Est2 (triangulation.n_active_cells());

	h_Conv_Est.reinit(triangulation.n_active_cells());
	p_Conv_Est.reinit(triangulation.n_active_cells());
	hp_Conv_Est.reinit(triangulation.n_active_cells());

	unsigned int cell_index=0;
	unsigned int patch_number=0;

	typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell , ++cell_index, ++patch_number)
	{
		//if(cell_index==1)
		//{
		double indicator_per_cell =0.0;

		double h_convergence_est_per_cell;
		unsigned int h_workload_num;
		h_patch_conv_load_no (cycle ,h_convergence_est_per_cell,h_workload_num, cell, patch_number);

		h_Conv_Est(cell_index)=h_convergence_est_per_cell;
		h_convergence_est_per_cell = h_convergence_est_per_cell /est_per_cell(cell_index);

		double p_convergence_est_per_cell;
		unsigned int p_workload_num;
		p_patch_conv_load_no (cycle ,p_convergence_est_per_cell,p_workload_num, cell, patch_number);

		p_Conv_Est(cell_index)=p_convergence_est_per_cell;
		p_convergence_est_per_cell = p_convergence_est_per_cell /est_per_cell(cell_index);

		double h_ratio= h_convergence_est_per_cell /  h_workload_num ;
		//double h_ratio= h_convergence_est_per_cell;
		double p_ratio= p_convergence_est_per_cell /  p_workload_num ;
		//double p_ratio= p_convergence_est_per_cell;

		if (h_ratio > p_ratio)
		{
			convergence_est_per_cell(cell_index)=h_convergence_est_per_cell;
			indicator_per_cell= convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
			hp_Conv_Est(cell_index)=indicator_per_cell;
			hp_Conv_Est2(cell_index)=h_Conv_Est(cell_index);
			p_ref_map[cell] = false;

			std::cout<< "H-refinement_marking ...  =>  p_ref_map[cell] = " << p_ref_map[cell] << std::endl;
			p_ref_map.insert (std::make_pair(cell, false));

		}

		else
		{

			convergence_est_per_cell(cell_index)=p_convergence_est_per_cell;
			indicator_per_cell=convergence_est_per_cell(cell_index)*est_per_cell(cell_index);
			hp_Conv_Est(cell_index)=indicator_per_cell;
			hp_Conv_Est2(cell_index)=p_Conv_Est(cell_index);
			p_ref_map[cell] = true;

			std::cout<< "P-refinement_marking ...  =>  p_ref_map[cell] = "  << p_ref_map[cell] << std::endl;
			p_ref_map.insert (std::make_pair(cell, true));

		}  //else

		to_be_sorted.push_back(std::make_pair(indicator_per_cell,cell));
		//to_be_sorted.push_back(std::make_pair(convergence_est_per_cell(cell_index),cell));

		//  }   // if index==1

	}// cell
	double min_conv_estimator = convergence_est_per_cell(0);
	//double min_conv_estimator = hp_Conv_Est(0);

	unsigned int index=0;

	typename hp::DoFHandler<dim>::active_cell_iterator
	celll = dof_handler.begin_active(),
	endcl = dof_handler.end();
	for (; celll!=endcl; ++celll , ++index)
	{
		min_conv_estimator = std::min ( min_conv_estimator, convergence_est_per_cell(index) );
		//min_conv_estimator = std::min ( min_conv_estimator, hp_Conv_Est(index) );
		std::cout<<std::endl;
		std::cout<< "convergence_est_per_cell [" << index << " ] ="<<  convergence_est_per_cell(index) << std::endl;
		std::cout<< "est_per_cell [" << index << " ] ="<< est_per_cell(index) << std::endl;
		std::cout<< "  hp_Conv_Est(index) [" << index << " ] ="<<   hp_Conv_Est(index) << " ?= " <<  " hp 2 _Conv_Est(index)[" << index << " ] ="<<   hp_Conv_Est2(index) << std::endl;
		std::cout<< "refinement_marking ...  =>  p-ref==1,  h-ref==0 : "  << p_ref_map[celll] << std::endl;
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<< "min_conv_estimator = " << min_conv_estimator << std::endl;
    //double theta= min_conv_estimator-0.0001;
	double theta= 0.5;
	std::cout<< "theta = " << theta << std::endl;
	std::cout<<std::endl;
	std::sort (to_be_sorted.begin(), to_be_sorted.end(), std_cxx1x::bind(&StokesProblem<dim>::decreasing,this,std_cxx1x::_1,std_cxx1x::_2));

	//for (unsigned int i=0; i< to_be_sorted. size(); ++i)
	//  std::cout<< to_be_sorted[i].first << std::endl;

	double L2_norm=est_per_cell.l2_norm();
	double sum=0;
	for (unsigned int i=0; i< to_be_sorted. size(); ++i)
	{
		to_be_sorted_with_refine_info.push_back(std::make_pair(to_be_sorted[i].second, p_ref_map[to_be_sorted[i].second]));
		typename hp::DoFHandler<dim>::active_cell_iterator  cell_sort=to_be_sorted[i].second;
		sum+= (to_be_sorted[i].first)*(to_be_sorted[i].first);
		//  std::cout<< "SUM :  " << sum << " &" << " (theta*(L2_norm))^2 " << (theta*(L2_norm))*(theta*(L2_norm)) << std::endl;

		candidate_cell_set.push_back (cell_sort);
		// std::cout<< "p_ref_map[cell_sort] :  " << p_ref_map[cell_sort] << std::endl;
		// std::cout<< "to_be_sorted_with_refine_info[i].second : " << to_be_sorted_with_refine_info[i].second << std::endl;

		if (sum >= (theta*(L2_norm))*(theta*(L2_norm)))

			break;
	}
	unsigned int n= candidate_cell_set.size();
	std::cout<< std::endl;
	std::cout<< "number of candidate_cell_set "<< n << std::endl;
	std::cout<< std::endl;
	//..............................................................................................................................

	marked_cells =0.;
	unsigned int i=0;
	typename hp::DoFHandler<dim>::active_cell_iterator  cel = dof_handler.begin_active(),
			endcel = dof_handler.end();
	for (; cel!=endcel; ++cel,++i)
	{
		typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::iterator  mark_candidate;
		for (mark_candidate=candidate_cell_set.begin(); mark_candidate!=candidate_cell_set.end(); ++ mark_candidate)
		{
			if (cel == (*mark_candidate))
				marked_cells(i)=1;

		}
	}

		}// marking_cells ()

//.................................................................................................................................
//Output_result
template <int dim>
/*
void StokesProblem <dim>::output_results (const unsigned int cycle , Vector<float> & marked_cells , Vector<double> &est_per_cell , Vector<double> &error_per_cell, Vector<double> &Vect_Pressure_Err, Vector<double> &Vect_grad_Velocity_Err ,
   		Vector<double> & h_Conv_Est, Vector<double> &p_Conv_Est, Vector<double> &hp_Conv_Est )
 */  		
   		

 // for uniform h-refinemnet
void StokesProblem <dim>::output_results (const unsigned int cycle , Vector<double> & est_per_cell , Vector<double> & error_per_cell, Vector<double> & Vect_Pressure_Err, Vector<double> & Vect_grad_Velocity_Err, Vector<double> & Vec_Velocity_Err )

{

 Vector<float> fe_degrees (triangulation.n_active_cells());
 {
   typename hp::DoFHandler<dim>::active_cell_iterator
   cell = dof_handler.begin_active(),
   endc = dof_handler.end();
   for (unsigned int index=0; cell!=endc; ++cell, ++index)
     fe_degrees(index)= fe_collection[cell->active_fe_index()].degree;
 }

 std::vector<std::string> solution_names (dim, "velocity");

 //std::vector<std::string> solution_names;
 solution_names.push_back ("x_velocity");
 solution_names.push_back ("y_velocity");
 solution_names.push_back ("pressure");

 //std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
 //(dim, DataComponentInterpretation::component_is_part_of_vector);

 std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation
 (dim, DataComponentInterpretation::component_is_part_of_vector);

 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);
 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);
 data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);


 DataOut<dim,hp::DoFHandler<dim> > data_out;

 data_out.attach_dof_handler (dof_handler);
 

 data_out.add_data_vector (solution, solution_names,DataOut<dim,hp::DoFHandler<dim> >::type_dof_data,data_component_interpretation);


 //data_out.add_data_vector (marked_cells, "marked_cells");
 data_out.add_data_vector (fe_degrees, "fe_degree");
 data_out.add_data_vector (est_per_cell, "Error_Estimator");
 data_out.add_data_vector (error_per_cell, "Error");
 data_out.add_data_vector (Vect_Pressure_Err, "Pressure_Error");
 data_out.add_data_vector (Vect_grad_Velocity_Err, "Grad_Velocity_Error");
 data_out.add_data_vector (Vec_Velocity_Err, "Vec_Velocity_Err");
 
/*
 data_out.add_data_vector (h_Conv_Est, "h_refine_Conv_Est");
 data_out.add_data_vector (p_Conv_Est, "p_refine_Conv_Est");
 data_out.add_data_vector (hp_Conv_Est, "hp_refine_Conv_Est");
*/

 data_out.build_patches ();
 std::string filename = "solution-" +
		  Utilities::int_to_string (cycle, 2) +".vtu";
 std::ofstream output (filename.c_str());
 data_out.write_vtu (output);
}

/*..............................................................................................................*/
//refine_in_h_p()

template <int dim>
void StokesProblem <dim>:: refine_in_h_p (const unsigned int cycle, std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> &candidate_cell_set,
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > &p_ref_map )
		{
	bool need_to_h_refine=false;
	unsigned int candidate_INDEX=0;

	typename std::vector<typename hp::DoFHandler<dim>::active_cell_iterator>::iterator  cell_candidate;
	for (cell_candidate=candidate_cell_set.begin(); cell_candidate!=candidate_cell_set.end(); ++ cell_candidate, ++candidate_INDEX)
	{
		std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch_cells = get_patch_around_cell (* cell_candidate);
		unsigned int level_p_refine = static_cast<unsigned int>((*cell_candidate)->active_fe_index());
		unsigned int level_h_refine = static_cast<unsigned int>((*cell_candidate)->level());

		std::cout<< "  p_ref_map[*cell_candidate] :  " << p_ref_map[*cell_candidate] << std::endl;

		if (p_ref_map[*cell_candidate]==false)
		{
			need_to_h_refine = true;
			(*cell_candidate)->set_refine_flag();
			std::cout<< "h-refinement flags" <<std::endl;

		}//if

		else if(p_ref_map[*cell_candidate]==true)
		{

			if ( ((*cell_candidate)-> active_fe_index()+ 1) <  (fe_collection.size()-1)  )
			{
				(*cell_candidate)->set_active_fe_index ((*cell_candidate)->active_fe_index() + 1);

				std::cout<< "p-refinement flags" <<std::endl;
			}

		}

	}



	if ( need_to_h_refine==true)
		triangulation.execute_coarsening_and_refinement();
	//std::cout<< "adaptive h-refinement is done" << std::endl;



bool cell_changed=false;
do
{
cell_changed=false;
unsigned int count_cell=0;
typename hp::DoFHandler<dim>::active_cell_iterator
 cell = dof_handler.begin_active(),
 endc = dof_handler.end();
 for (; cell!=endc; ++cell,++count_cell)
   {
	// std::cout<<std::endl;
        std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> patch_cells = get_patch_around_cell (cell);
      //  std::cout<< "count_cell 1= "<< count_cell <<std::endl;
      //  unsigned int level_p_refine = static_cast<unsigned int>(cell->active_fe_index());
       // unsigned int level_h_refine = static_cast<unsigned int>(cell->level());
      //  std::cout<< "patch_cells.size() "<< patch_cells.size() <<std::endl;

     for (unsigned int i=1; i<patch_cells.size(); ++i)
       {

         if (patch_cells[i]->active_fe_index()+1 < (patch_cells[0]->active_fe_index()))
          {
        	// std::cout<< "count_cell 2= "<< count_cell <<std::endl;
        	// std::cout<<std::endl;
        	// std::cout <<  " patch_cells[i]->active_fe_index()= " << patch_cells[i]->active_fe_index()  << "    " <<  "patch_cells[0]->active_fe_index()= " << patch_cells[0]->active_fe_index()-1 <<std::endl;
        	// std::cout<<std::endl;
             patch_cells[i]->set_active_fe_index(patch_cells[0]->active_fe_index()-1);
             cell_changed=true;
          }


          else if (patch_cells[i]->active_fe_index() > (patch_cells[0]->active_fe_index()+1))
                {
        	     // std::cout<< "count_cell 3= "<< count_cell <<std::endl;
                  patch_cells[0]-> set_active_fe_index (patch_cells[i]->active_fe_index()-1);
                  cell_changed=true;

                }

          }// patch_cells

}// for cell

}//do
while (cell_changed==true);


		}//hp_refine


/*......................................................................................................................................................*/
template <int dim>
void StokesProblem <dim>::run()
{
	for (unsigned int cycle=0; cycle<150; ++cycle)
	{

		std::cout<< std::endl;

		std::cout<< "-----------------------------------------------------------" << std::endl;
		std::cout<< std::endl;
		std::cout << "Cycle " << cycle << ':' << std::endl;
		if (cycle == 0)
		{
			generate_mesh();
			set_global_active_fe_indices(dof_handler);

		}
		std::cout<<"Number of active cells: "<< triangulation.n_active_cells() << std::endl;

		std::cout<<"Total number of cells: " << triangulation.n_cells() << std::endl ;
		setup_system ();
		assemble_system();
		solve ();

		Vector<double> error_per_cell (triangulation.n_active_cells());
		Vector<double> Vect_Pressure_Err(triangulation.n_active_cells());
		Vector<double> Vect_grad_Velocity_Err(triangulation.n_active_cells());
		Vector<double> Vect_Velocity_Err(triangulation.n_active_cells());
		
		compute_error  (error_per_cell, Vect_Pressure_Err, Vect_grad_Velocity_Err, Vect_Velocity_Err);
		Vector<double> est_per_cell (triangulation.n_active_cells());
		estimate(est_per_cell);


		// for uniform refinement:
		double L2_norm_est= est_per_cell.l2_norm();
		std::cout<< "L2_norm of ERROR Estimate is: "<< L2_norm_est << std::endl;

		//  output_results (cycle , est_per_cell , error_per_cell, Vect_Pressure_Err, Vect_grad_Velocity_Err );    
		triangulation.refine_global (1);

		//  std::cout<< "Vector of Error Estimate: "<< est_per_cell << std::endl;
	
		

/*		
// for adaptive h- and p- refinment:
		double L2_norm_est= est_per_cell.l2_norm();
		std::cout<< "L2_norm of ERROR Estimate is: "<< L2_norm_est << std::endl;
		std::cout<<std::endl;
		Vector<float> marked_cells(triangulation.n_active_cells());
		std::vector<typename hp::DoFHandler<dim>::active_cell_iterator> candidate_cell_set;
		std::map<typename hp::DoFHandler<dim>::active_cell_iterator, bool > p_ref_map;

		Vector<double> h_Conv_Est;
		Vector<double> p_Conv_Est;
		Vector<double> hp_Conv_Est;

		marking_cells(cycle,  marked_cells, candidate_cell_set, p_ref_map,  h_Conv_Est, p_Conv_Est , hp_Conv_Est);
		output_results(cycle, marked_cells, est_per_cell, error_per_cell, Vect_Pressure_Err, Vect_grad_Velocity_Err,  h_Conv_Est, p_Conv_Est , hp_Conv_Est);

		refine_in_h_p(cycle,  candidate_cell_set, p_ref_map);

		if (L2_norm_est < Tolerance)
			break;
*/		 


	}// cycle

}//run

//Explicit initialization

template class StokesProblem<2>;

