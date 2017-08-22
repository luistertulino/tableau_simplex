#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

int num_variables, num_restrictions, total_width;
int min_max; // Indicates if the problem is of minimization (0) or maximization (1)
double z = 0; // The inicial value of the objective function
bool there_is_bigger = false; // Indicates if there is some bigger than or equal restriction
int count_bigger = 0; //Count the number of >= restrictions
int iterations = 1;

string separator = "----------------------------------------------------------------------";

void insert_problem_data(double objective_function[], double restrictions_table[], double restrictions[], int less_or_bigger[]){

	cout << "\nInsert -1 if it is a minimization problem or 1 if it is a maximization problem: "; cin >> min_max;

	cout << "\nInsert the coefficients of the objective function, separated by space:\n";
	{
		for (int i = 0; i < num_variables; i++)
		{
			cin >> objective_function[i]; objective_function[i] *= min_max;
		}
		for (int i = num_variables; i < total_width; ++i)
		{
			objective_function[i] = 0;
		}
	}

	cout << "\nInsert the restrictions.\n" <<
			 "Type the coefficients separated by space, the type of restriction (<= or >=), and, in the end, type the restriction.\n" <<
			 "If a variable does not appear in a restriction, use 0 as coefficient.\n\n";
	{
		for (int i = 0; i < num_restrictions; i++)
		{	

			// reads the coefficients of the restriction
			for (int j = 0; j < num_variables; j++)
			{
				cin >> restrictions_table[i*total_width+j];
			}

			string type_restriction;
			cin >> type_restriction;

			if (type_restriction.compare("<=") == 0)
			{
				less_or_bigger[i] = 1;
			}
			else if (type_restriction.compare("=") == 0)
			{
				less_or_bigger[i] = 0;
				there_is_bigger = true;
				count_bigger++;
			}
			else {
				less_or_bigger[i] = -1;
				there_is_bigger = true;
				count_bigger++;
			}

			cin >> restrictions[i];

			for (int j = num_variables; j < total_width; j++)
			{
				restrictions_table[i*total_width+j] = (j - num_variables == i) ? less_or_bigger[i] : 0;
				/* 
					If the current collumn has the same number of the current restriction, then this slots is the slack variable
					of the restriction, so it must receive 1 or -1
				*/
			}
		}
	}
}

void print_tableau(double objective_function[], double restrictions_table[], double restrictions[], int base_variables[]){
	cout << "\n";

	cout << "z = " << z << "  "; cout << "z_j - c_j = ";
	for (int i = 0; i < total_width; ++i)
	{
		cout << objective_function[i] << " ";
	}
	cout << "\n\n";

	for (int i = 0; i < num_restrictions; i++)
	{
		cout << "x_" << base_variables[i] << " " << restrictions[i] << "               ";

		for (int j = 0; j < total_width; j++)
		{
			cout << restrictions_table[i*total_width+j] << " ";
		}
		cout << "\n";
	}

	cout << separator;
}

void first_phase_print_tableau(double z_function[], double w_function[], double restrictions_table[], double restrictions[], int base_variables[], int new_total_width, double w){
	cout << "\n";

	cout << "w = " << w << "  "; cout << "w_j - c_j = ";
	for (int i = 0; i < new_total_width; ++i)
	{
		cout << w_function[i] << " ";
	}
	cout << "\n";

	cout << "z = " << z << "  "; cout << "z_j - c_j = ";
	for (int i = 0; i < total_width; ++i)
	{
		cout << z_function[i] << " ";
	}
	cout << "\n\n";

	for (int i = 0; i < num_restrictions; i++)
	{
		cout << "x_" << base_variables[i] << " " << restrictions[i] << "               ";

		for (int j = 0; j < new_total_width; j++)
		{
			cout << restrictions_table[i*new_total_width+j] << " ";
		}
		cout << "\n";
	}

	cout << separator;
}

void find_biggest_zj_cj(double &biggest_zj_cj, int &index_biggest_zj_cj, double objective_function[], bool not_in_base[]){

	int i = 0;

	while(!not_in_base[i]) { i++; }

	biggest_zj_cj = objective_function[i];
	index_biggest_zj_cj = i;
	/* Here, we are defining the initial biggest zj-cj as the first j that are not in the base, 
	so that we don't need the set a huge value as a initial value */

	for ( ; i < total_width; i++)
	{
		if (!not_in_base[i]) continue; // This variable is in the base, so continue searching
		else {
			if (objective_function[i] > biggest_zj_cj)
			{
				// This strictly bigger guarantees the application of Bland's rule
				biggest_zj_cj = objective_function[i];
				index_biggest_zj_cj = i;
			}
		}
	}

}

int find_index_lowest_ratio(int base_variables[], double restrictions[], double restrictions_table[], int index_biggest_zj_cj, int &index_row_pivot){

	int index_lowest_ratio;
	double lowest_ratio;
	//bool found = false;

	int i = 0;
	while(i < num_restrictions && restrictions_table[i*total_width+index_biggest_zj_cj] <= 0){
		i++;
	}

	if (i == num_restrictions)
	{
		// In this case, the solution is unlimited.
		return -1;
	}

	/* 
		Note the difference between the row of the pivot (index_row_pivot = i;) 
		and the variable selected to leave the base (index_lowest_ratio = base_variables[i];)
	*/

	index_lowest_ratio = base_variables[i];
	lowest_ratio = restrictions[i] / restrictions_table[i*total_width+index_biggest_zj_cj];
	index_row_pivot = i;

	for (i = i+1; i < num_restrictions; i++)
	{
		if (restrictions_table[i*total_width+index_biggest_zj_cj] <= 0) { continue; }
		else if (restrictions[i] / restrictions_table[i*total_width+index_biggest_zj_cj] < lowest_ratio)
		{
			lowest_ratio = restrictions[i] / restrictions_table[i*total_width+index_biggest_zj_cj];
			index_lowest_ratio = base_variables[i];
			index_row_pivot = i;
		}
		else if (restrictions[i] / restrictions_table[i*total_width+index_biggest_zj_cj] == lowest_ratio && base_variables[i] < index_lowest_ratio)
		{
			lowest_ratio = restrictions[i] / restrictions_table[i*total_width+index_biggest_zj_cj];
			index_lowest_ratio = base_variables[i];
			index_row_pivot = i;
		}
		else { continue; }
	}

	return index_lowest_ratio;
	
}

int first_phase(double objective_function[], double restrictions_table[], double restrictions[], int base_variables[], int less_or_bigger[]){

	int aux_total_width = total_width + count_bigger;

	double first_phase_objective_function[aux_total_width];
	double first_phase_restrictions_table[num_restrictions*aux_total_width];
	bool not_in_base[aux_total_width];

	int w = 0; // The value of first phase function
	int artificial_variables = 0; // Number of artificial variables used

	// The initialization of the objective function of the 1s phase
	{
		for (int i = 0; i < total_width; i++)
		{
			first_phase_objective_function[i] = 0;
		}
		for (int i = total_width; i < aux_total_width; i++)
		{
			first_phase_objective_function[i] = -1;
		}
	

		for (int i = 0; i < aux_total_width; i++)
		{
			not_in_base[i] = true;
		}		

		// Replace the slack variables by artificial variables

		for (int i = 0; i < num_restrictions; i++)
		{
			int x_i = base_variables[i];

			if (less_or_bigger[i] == 1)
			{
				// This is a <= restriction: no need for artificial varible
				// And this variable will be in the base
				not_in_base[x_i] = false;
			}
			else{
				// This is a >= restriction: we need to "create" a artificial variable
				not_in_base[total_width+artificial_variables] = false;
				base_variables[i] = total_width + artificial_variables;
				// The original base variable will not be in the base
				artificial_variables++;

				w += restrictions[i];
			}
		}

		// We generate the new tableau and update the wj - cj
		for (int i = 0; i < num_restrictions; i++)
		{
			for (int j = 0; j < total_width; j++)
			{
				// These are the values of the original tableau
				first_phase_restrictions_table[i*aux_total_width+j] = restrictions_table[i*total_width+j];

				if (less_or_bigger[i] == -1 || less_or_bigger[i] == 0)
				{
					first_phase_objective_function[j] += first_phase_restrictions_table[i*aux_total_width+j];
				}
			}

			for (int j = total_width; j < aux_total_width; j++)
			{
				/* These are the values according to the artificial variables
				   If the current column equals to the base variable of the current row,
				   this is where the artificial variable has 1 as coefficent
				*/
				if (base_variables[i] == j)
				{
					first_phase_restrictions_table[i*aux_total_width+j] = 1.0;
				}
				else{
					first_phase_restrictions_table[i*aux_total_width+j] = 0.0;
				}

				if (less_or_bigger[i] == -1 || less_or_bigger[i] == 0)
				{
					first_phase_objective_function[j] += first_phase_restrictions_table[i*aux_total_width+j];
				}
			}
		}
	}

	cout << "\n" << "-----------------------------------FASE ONE-----------------------------------" << "\n\n";

	cout << "\nINITIAL TABLEAU:\n\n";

	first_phase_print_tableau(objective_function, first_phase_objective_function, first_phase_restrictions_table, restrictions, base_variables, aux_total_width, w);

	{
		char a;
		cout << "\nType 'c' to continue.\n";
		cin >> a;
		while(a != 'c'){
			cout << "Wrong letter! Type 'c' to continue.\n";
			cin >> a;
		}
	}	

	double biggest_wj_cj;
	int index_biggest_wj_cj;

	swap(total_width, aux_total_width);
	find_biggest_zj_cj(biggest_wj_cj, index_biggest_wj_cj, first_phase_objective_function, not_in_base);
	// In the 1st phase, the correct width is aux_total_width, so we need to swap and re-swap (to not lose the original total_width)
	swap(total_width, aux_total_width);

	while(biggest_wj_cj > 0){

		cout << "\nIteration " << iterations << endl;

		swap(total_width, aux_total_width);

		int index_row_pivot; // The row where the pivot is
		int index_lowest_ratio = find_index_lowest_ratio(base_variables, restrictions, first_phase_restrictions_table, index_biggest_wj_cj, index_row_pivot);
		// The varible with the lowest ratio b_i/y_ik, so that that varible is gonna leave the base

		swap(total_width, aux_total_width);

		// Note the difference between the row of the pivot and the variable selected to leave the base

		if (index_lowest_ratio != -1)
		{
			double b_i = restrictions[index_row_pivot];
			double y_ik = first_phase_restrictions_table[index_row_pivot*aux_total_width+index_biggest_wj_cj];
			cout << "Variable to enter the base: x_" << index_biggest_wj_cj << ". "
				 << "Value of z_j - c_j = " << biggest_wj_cj << "\n"
			 	 << "Variable to leave the base: x_" << index_lowest_ratio << ". "
			 	 << "Value of the ratio b_i/y_ik = " << b_i << "/" << y_ik << " = " << b_i/y_ik << "\n";


			// Here, we are replacing the basic variable with the lowest ratio for the variable with the biggest zj-cj
			not_in_base[index_biggest_wj_cj] = false;
			not_in_base[index_lowest_ratio] = true;
			base_variables[index_row_pivot] = index_biggest_wj_cj;

			if (index_lowest_ratio >= total_width)
			{
				artificial_variables--;
			}
			
			// Now, we divide de pivot's row by the value of the pivot
			double pivot = first_phase_restrictions_table[index_row_pivot*aux_total_width + index_biggest_wj_cj];
			for (int j = 0; j < aux_total_width; j++)
			{			
				first_phase_restrictions_table[index_row_pivot*aux_total_width + j] /= pivot;
			}
			restrictions[index_row_pivot] /= pivot;

			// Now, we update the restrictions' table
			for (int i = 0; i < num_restrictions; i++)
			{
				if (i == index_row_pivot) continue;

				// For each row of the table
				double value_to_canonize = first_phase_restrictions_table[i*aux_total_width+index_biggest_wj_cj];

				for (int j = 0; j < aux_total_width; j++)
				{
					// For each column of the row i
					
					double value_in_pivot_row = first_phase_restrictions_table[index_row_pivot*aux_total_width+j];
					first_phase_restrictions_table[i*aux_total_width+j] += (-1) * value_to_canonize * value_in_pivot_row;
				}

				// We also update the restrictions
				restrictions[i] += (-1) * value_to_canonize * restrictions[index_row_pivot];
			}

			//Now, we update the array of wj-cj
			double value_to_canonize_wj_cj = first_phase_objective_function[index_biggest_wj_cj];
			for (int j = 0; j < aux_total_width; j++)
			{
				first_phase_objective_function[j] += (-1) * value_to_canonize_wj_cj * first_phase_restrictions_table[index_row_pivot*aux_total_width+j];
			}

			// We update the value of objective function
			w += (-1) * value_to_canonize_wj_cj * restrictions[index_row_pivot];

			double value_to_canonize_zj_cj = objective_function[index_biggest_wj_cj];
			for (int j = 0; j < total_width; j++)
			{
				objective_function[j] += (-1) * value_to_canonize_zj_cj * first_phase_restrictions_table[index_row_pivot*aux_total_width+j];
			}
			
			z += (-1) * value_to_canonize_zj_cj * restrictions[index_row_pivot];

		}
		else{
			cout << "\n\n" << "This problem is unlimited!\n\n";
			return -1;
		}

		iterations++;

		cout << "\nNew tableau:\n";
		first_phase_print_tableau(objective_function, first_phase_objective_function, first_phase_restrictions_table, restrictions, base_variables, aux_total_width, w);

		swap(total_width, aux_total_width);
		find_biggest_zj_cj(biggest_wj_cj, index_biggest_wj_cj, first_phase_objective_function, not_in_base);
		// In the 1st phase, the correct width is aux_total_width, so we need to swap and re-swap (to not lose the original total_width)
		swap(total_width, aux_total_width);

		{
			char a;
			cout << "\nType 'c' to continue.\n";
			cin >> a;
			while(a != 'c'){
				cout << "Wrong letter! Type 'c' to continue.\n";
				cin >> a;
			}
		}

	}

	for (int i = 0; i < num_restrictions; i++)
	{
		for (int j = 0; j < total_width; j++)
		{
			restrictions_table[i*total_width+j] = first_phase_restrictions_table[i*aux_total_width+j];
		}
	}

	cout << "\n" << "-------------------------------END OF FASE ONE--------------------------------" << "\n\n";

	return artificial_variables;
}

void pivoting(double objective_function[], double restrictions_table[], double restrictions[], int less_or_bigger[]){

	int base_variables[num_restrictions]; // Variables indexes that composes the base
	bool not_in_base[total_width]; /* Array that indicates if the variables are or are not out of the base
										false -> is in the base
										true -> is NOT in the base. In this case, such variable is a candidate to enter the base
	
									*/

	{
		for (int i = 0; i < num_variables; i++)
		{
			not_in_base[i] = true;
		}
		for (int i = num_variables; i < total_width; i++)
		{
			not_in_base[i] = false;
		}

		for (int i = 0; i <= num_restrictions; i++)
		{
			base_variables[i] = i + num_variables;
			// The initial basic variable are the slack variables
		}
	}

	if (there_is_bigger)
	{
		int a_v = first_phase(objective_function, restrictions_table, restrictions, base_variables, less_or_bigger);

		if (a_v == 0)
		{
			for (int j = 0; j < total_width; j++)
			{
				not_in_base[j] = true;
			}
			for (int i = 0; i < num_restrictions; i++)
			{
				not_in_base[base_variables[i]] = false;
			}

			cout << "\n" << "-----------------------------------FASE TWO-----------------------------------" << "\n\n";
		}
		else{
			// There is some artificial variable in the base
			cout << "\n" << "There is " << a_v << " artificial variables in the base." << "\n" << "The problem is infeasible!" << "\n";
			return;
		}
	}
	else{
		cout << "\nINITIAL TABLEAU:\n\n";
	}

	print_tableau(objective_function, restrictions_table, restrictions, base_variables);

	double biggest_zj_cj;
	int index_biggest_zj_cj;

	find_biggest_zj_cj(biggest_zj_cj, index_biggest_zj_cj, objective_function, not_in_base);

	char a;
	cout << "\nType 'c' to continue.\n";
	cin >> a;
	while(a != 'c'){
		cout << "Wrong letter! Type 'c' to continue.\n";
		cin >> a;
	}

	while(biggest_zj_cj > 0){
		// If the biggest zj-cj <= 0, the stopping criterion is met

		cout << "\nIteration " << iterations << endl;

		int index_row_pivot; // The row where the pivot is
		int index_lowest_ratio = find_index_lowest_ratio(base_variables, restrictions, restrictions_table, index_biggest_zj_cj, index_row_pivot);
		// The varible with the lowest ratio b_i/y_ik, so that that varible is gonna leave the base

		// Note the difference between the row of the pivot and the variable selected to leave the base

		if (index_lowest_ratio != -1)
		{
			double b_i = restrictions[index_row_pivot];
			double y_ik = restrictions_table[index_row_pivot*total_width+index_biggest_zj_cj];
			cout << "Variable to enter the base: x_" << index_biggest_zj_cj << ". "
				 << "Value of z_j - c_j = " << biggest_zj_cj << "\n"
			 	 << "Variable to leave the base: x_" << index_lowest_ratio << ". "
			 	 << "Value of the ratio b_i/y_ik = " << b_i << "/" << y_ik << " = " << b_i/y_ik << "\n";


			// Here, we are replacing the basic variable with the lowest ratio for the variable with the biggest zj-cj
			not_in_base[index_biggest_zj_cj] = false;
			not_in_base[index_lowest_ratio] = true;
			base_variables[index_row_pivot] = index_biggest_zj_cj;
			
			// Now, we divide de pivot's row by the value of the pivot
			double pivot = restrictions_table[index_row_pivot*total_width + index_biggest_zj_cj];
			for (int j = 0; j < total_width; j++)
			{			
				restrictions_table[index_row_pivot*total_width + j] /= pivot;
			}
			restrictions[index_row_pivot] /= pivot;

			// Now, we update the restrictions' table
			for (int i = 0; i < num_restrictions; i++)
			{
				if (i == index_row_pivot) continue;

				// For each row of the table
				double value_to_canonize = restrictions_table[i*total_width+index_biggest_zj_cj];

				//cout << "\nvalue_to_canonize: " << value_to_canonize << "\n";

				for (int j = 0; j < total_width; j++)
				{
					// For each column of the row i
					
					double value_in_pivot_row = restrictions_table[index_row_pivot*total_width+j];
					restrictions_table[i*total_width+j] += (-1) * value_to_canonize * value_in_pivot_row;
				}

				// We also update the restrictions
				restrictions[i] += (-1) * value_to_canonize * restrictions[index_row_pivot];
			}

			// Now, we update the array of zj-cj
			double value_to_canonize_zj_cj = objective_function[index_biggest_zj_cj];
			for (int j = 0; j < total_width; j++)
			{
				objective_function[j] += (-1) * value_to_canonize_zj_cj * restrictions_table[index_row_pivot*total_width+j];
			}

			// And, finally, we update the value of objective function
			z += (-1) * value_to_canonize_zj_cj * restrictions[index_row_pivot];

		}
		else{
			cout << "\n\n" << "This problem is unlimited!\n\n";
			return;
		}

		iterations++;

		cout << "\nNew tableau:\n";
		print_tableau(objective_function, restrictions_table, restrictions, base_variables);

		find_biggest_zj_cj(biggest_zj_cj, index_biggest_zj_cj, objective_function, not_in_base);

		//cout << "biggest_zj_cj: " << biggest_zj_cj << endl;
		cout << "\nType 'c' to continue.\n";
		cin >> a;
		while(a != 'c'){
			cout << "Wrong letter! Type 'c' to continue.\n";
			cin >> a;
		}

	}


}

int main()
{

	cout << "\nInsert the number of variables of the linear program: "; cin >> num_variables;
	cout << "\nInsert the number of restrictions: "; cin >> num_restrictions;
	total_width = num_variables + num_restrictions;
	/*
		The complete "matrix" is composed of the natural coefficients and the slack variables, that's why we use this sum
	*/

	double objective_function[total_width], restrictions_table[num_restrictions*total_width], restrictions[num_restrictions];
	int less_or_bigger[num_restrictions];

	insert_problem_data(objective_function, restrictions_table, restrictions, less_or_bigger);

	pivoting(objective_function, restrictions_table, restrictions, less_or_bigger);

	return 0;
}