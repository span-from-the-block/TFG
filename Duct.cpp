#include <cstdint>
#include <algorithm>
#include <tuple>
#include <atomic>
#include <thread>
#include <cmath>
#include <iostream>
#include <cstdint>
#include <sstream>
#include <random>
#include <vector>
#include "Duct.h"
//#include "cyclic_barrier.h"

using namespace std;

struct Coordinate
{
	public:
		int x;
		int y;
		int z;

		Coordinate(int x_ = 0, int y_ = 0, int z_ = 0) : x(x_), y(y_), z(z_) {}
		bool operator == (const Coordinate& c) const { return x == c.x && y == c.y && z == c.z; }
		bool operator != (const Coordinate& c) const { return !(x == c.x && y == c.y && z == c.z); }
};

static const Coordinate null_coordinate(0, 0, 0);

thread_local auto rng = default_random_engine{};
thread_local std::random_device rd;
thread_local std::mt19937 re(rd());


thread_local int Duct::total_population = 0;
thread_local int Duct::dysfunctional_population = 0;
thread_local int Duct::total_reproductions = 0;
thread_local int Duct::dysfunctional_reproductions = 0;

thread_local int Duct::stem_counter = 0;
thread_local int Duct::bipotent_counter = 0;
thread_local int Duct::proglumepi_counter = 0;
thread_local int Duct::lumepi_counter = 0;
thread_local int Duct::mut_stem_counter = 0;
thread_local int Duct::mut_bipotent_counter = 0;
thread_local int Duct::mut_proglumepi_counter = 0;
thread_local int Duct::mut_lumepi_counter = 0;

thread_local int Duct::thread_index = 0;
thread_local vector<int> Duct::index_array;
thread_local Coordinate  Duct::boundaries = null_coordinate;

//Not a member function!
float Duct::nextInt(int begin, int end)
{
	std::uniform_int_distribution<int> random(begin, end);
	return random(re);
}

int  Duct::get_generation_number() { return generation_number; }
int  Duct::get_length()            { return length; }
int  Duct::get_radius()            { return radius; }
int  Duct::get_diameter()          { return diameter; }
int  Duct::get_capacity()          { return capacity; }
bool Duct::get_HGP()               { return HGP; }
int  Duct::is_DCIS_achieved()      { return DCIS_achieved ? DCIS_generation : 0; }

Duct::Duct(int l, int r, bool h, int n) :
	length(l),
	radius(r),
	HGP(h),
	diameter(r*2),
	n_threads(n),
	local_initialized(n,false),
	duct_type         (l, vector<vector<uint8_t>>(r*2, vector<uint8_t>(r*2, 0))),
	duct_housekeeping (l, vector<vector<uint8_t>>(r*2, vector<uint8_t>(r*2, 0))),
	duct_protooncogene(l, vector<vector<uint8_t>>(r*2, vector<uint8_t>(r*2, 0))),
	duct_supressor    (l, vector<vector<uint8_t>>(r*2, vector<uint8_t>(r*2, 0))),
	duct_apoptosis    (l, vector<vector<uint8_t>>(r*2, vector<uint8_t>(r*2, 0))),
	thread_locks(n-1)
{
	for(int z=0 ; z<l ; z++)
	{
		for(int x=0 ; x<=(r/2) ; x++)
		{
			for(int y=x ; y<=r ; y++)
			{
				if (get_point_radius(x,y) == r)
				{
					duct_type.at(z).at(x).at(y)         = BASAL;
					duct_type.at(z).at(x).at(2*r-y)     = BASAL;
					duct_type.at(z).at(2*r-x).at(y)     = BASAL;
					duct_type.at(z).at(2*r-x).at(2*r-y) = BASAL;
					// simetría diagonal
					duct_type.at(z).at(y).at(x)         = BASAL;
					duct_type.at(z).at(y).at(2*r-x)     = BASAL;
					duct_type.at(z).at(2*r-y).at(x)     = BASAL;
					duct_type.at(z).at(2*r-y).at(2*r-x) = BASAL;
				}
			}
		}
	}

	generation_number = 0;
	initialized 	  = false;
	DCIS_achieved     = false;
	DCIS_generation   = 0;

	capacity = 0;
	stem_capacity = 0;
	for(int i=0 ; i<diameter ; i++)
	{
		for(int j=0 ; j<diameter ; j++)
		{
			if(get_point_radius(i,j) < radius)
				capacity++;
			if(get_point_radius(i,j) == radius-1 || get_point_radius(i,j) == radius-2)
				stem_capacity++;
		}
	}
	capacity *= length;
	stem_capacity *= length;

	for(int i=0 ; i<length ; i++)
	{
		index_array.push_back(i);
	}
	if(n_threads > 1)
	{
		thread_barrier = new CyclicBarrier(n_threads + 1); // Hilos que lanzamos mas éste, el principal
	}

	reproductions_history.push_back(0);
	dysf_reproductions_history.push_back(0);
}

/* Constructor para crear hilos (begin y end inclusive en rango) */
void Duct::run(int begin, int end, int index)
{
	thread_index = index;
	boundaries = Coordinate(begin, end, 0);
	// [begin, end] ambos inclusive
	for(int z=begin ; z<=end ; z++)
	{
		index_array.push_back(z);
	}

	for (int i=0; i<total_generations; i++)
	{
		thread_barrier->await();
		next_generation();
		update_global_counters();
	}
}

void Duct::next_generation()
{
	total_population = 0;
	dysfunctional_population = 0;
	total_reproductions = 0;
	dysfunctional_reproductions = 0;

	stem_counter = 0;
	bipotent_counter = 0;
	proglumepi_counter = 0;
	lumepi_counter = 0;
	mut_stem_counter = 0;
	mut_bipotent_counter = 0;
	mut_proglumepi_counter = 0;
	mut_lumepi_counter = 0;

	vector<int> index_array_aux(index_array);

	if(!initialized)
	{
		local_initialized[thread_index] = true;
		shuffle(begin(index_array_aux), end(index_array_aux), rng);
	}

	for(uint32_t z_i = 0, z ; z_i < index_array.size() ; z_i++) {
		z = index_array_aux[z_i];

		bool left_boundary_flag = false;
		bool right_boundary_flag = false;

		if(n_threads > 1)
		{
			// Barrera izquierda
			if(boundaries.x != 0 && z <= boundaries.x+1)
			{
				left_boundary_flag = true;
				thread_locks[thread_index-1].lock();
			}
			// Barrera derecha
			if(boundaries.y != length-1 && z >= boundaries.y -1)
			{
				right_boundary_flag = true;
				thread_locks[thread_index].lock();
			}
		}

		for(int y=0 ; y<diameter ; y++)
		{
			for(int x=0 ; x<diameter ; x++)
			{
				Coordinate current_cell_coord = Coordinate(x,y,z);

				if(initialized)
					normal_routine(current_cell_coord);
				else
					init_routine(current_cell_coord);
			}
		}

		if(n_threads > 1)
		{
			if(left_boundary_flag)
				thread_locks[thread_index-1].unlock();
			if(right_boundary_flag)
				thread_locks[thread_index].unlock();
		}
	}
}

void Duct::execute(int nGens)
{
	if(generation_number == 0)
	{
		init_routine_stem();
	}

	total_generations = nGens;

	// SECUENCIAL
	if(n_threads == 1)
	{
		for (int i=0; i<total_generations; i++)
		{
			if(i%1000 == 0)
			{
			for(int z=0 ; z<length ; z++)
			{
				for(int y=0 ; y<diameter ; y++)
				{
					for(int x=0 ; x<diameter ; x++)
					{
						if(cell_is_null(Coordinate(x,y,z)))
						{
							cout << ".";
							continue;
						}
						switch(get_cell_type(Coordinate(x,y,z)))
						{
							case BASAL:
								cout << "O"; break;
							case PROG_STEM:
								cout << "S"; break;
							case PROG_BIPOTENT:
								cout << "B"; break;
							case PROG_MYOEPI:
								cout << "M"; break;
							case MYOEPI:
								cout << "m"; break;
							case PROG_LUMEPI:
								cout << "L"; break;
							case LUMEPI:
								cout << "l"; break;
						}
					}
					cout << endl;
				}
			}
			cout << endl << i << endl;
			}

			initialized = local_initialized[0];
			next_generation();

			generation_number++;
			total_history.push_back(total_population);
			dysfunctional_history.push_back(dysfunctional_population);
			reproductions_history.push_back(
				total_reproductions + reproductions_history[reproductions_history.size()-1]);
			dysf_reproductions_history.push_back(
				dysfunctional_reproductions + dysf_reproductions_history[dysf_reproductions_history.size()-1]);

			stem_mean = (float)mut_stem_counter/(float)stem_counter;
			bipotent_mean = (float)mut_bipotent_counter/(float)bipotent_counter;
			proglumepi_mean = (float)mut_proglumepi_counter/(float)proglumepi_counter;
			lumepi_mean = (float)mut_lumepi_counter/(float)lumepi_counter;

		}
	}
	// PARALELO
	else
	{
		int step = length/n_threads;
		for(int i=0 ; i<n_threads ; i++)
		{
			int begin =     i*step;
			int end   = (i+1)*step-1;
			if ((i+1) == n_threads)
				end = length-1;
			thread_array.push_back(thread(Duct::run, this, begin, end, i));
		}
		
		for (int i=0; i<total_generations ; i++)
		{
			if(i%1000 == 0)
			{
			for(int z=0 ; z<length ; z++)
			{
				for(int y=0 ; y<diameter ; y++)
				{
					for(int x=0 ; x<diameter ; x++)
					{
						if(cell_is_null(Coordinate(x,y,z)))
						{
							cout << ".";
							continue;
						}
						switch(get_cell_type(Coordinate(x,y,z)))
						{
							case BASAL:
								cout << "O"; break;
							case PROG_STEM:
								cout << "S"; break;
							case PROG_BIPOTENT:
								cout << "B"; break;
							case PROG_MYOEPI:
								cout << "M"; break;
							case MYOEPI:
								cout << "m"; break;
							case PROG_LUMEPI:
								cout << "L"; break;
							case LUMEPI:
								cout << "l"; break;
						}
					}
					cout << endl;
				}
			}
			cout << endl << i << endl;
			}
			
			if(!initialized)
			{
				bool init_flag = true;
				for (int t=0; t<n_threads; t++)
				{
					init_flag = init_flag && local_initialized[t];
				}
				initialized = init_flag;
			}
			thread_barrier->await();
		}

		//for (auto t : thread_array)
		for(int t=0 ; t<n_threads ; t++)
		{
			thread_array.at(t).join();
		}
	}
}

int Duct::get_point_radius(int x, int y)
{ return (int)(sqrt((x-radius)*(x-radius)
				   +(y-radius)*(y-radius))+1); }

void Duct::init_routine_stem()
{
	if(radius == 1)
		return;

	// 1. Draw the BASAL membrane by designating outer boundary points in a cylindrical orientation

	/* Hecho en constructor de SliceSec */

	// 2. Place STEM CELLS within the basal membrane
	//    (less than 5% of the population - 200 for a 24,000 cell pop. => [ ~ 0.84% ])

	// [ALT]: Con nuestro modelo de radio, las células que forman las dos capas iniciales son 20,800 contra las
	//		  24,000 del paper.
	float stem_cells_percentage = 1;  // Porcentaje del total de células que queremos que sean células madre
	int n_stem_cells = (int)((((float)stem_capacity) / 100) * stem_cells_percentage);

	if(n_stem_cells == 0)
		n_stem_cells = 1;

	int chosen_x, chosen_y, chosen_z;
	// Elegimos posición de células madre por polling de puntos aleatorios
	for(int counter = 1 ; counter <= n_stem_cells ; counter++)
	{
		// Elegir SliceSec aleatorio
		chosen_z = nextInt(0, length-1);
		// Recorrer slice, elegir punto aleatorio hasta encontrar uno de radio r-1
		do {
			chosen_x = nextInt(0, diameter-1);
			chosen_y = nextInt(0, diameter-1);
		} while (get_point_radius(chosen_x, chosen_y) != radius-1 ||
				 cell_is_null(Coordinate(chosen_x, chosen_y, chosen_z)) == false);

		Coordinate new_coordinate = Coordinate(chosen_x, chosen_y, chosen_z);
		set_new_stem_cell(new_coordinate);
	}
	total_population += n_stem_cells;
}


void Duct::init_routine(const Coordinate& current_cell_coord)
{
	// Si la coordenada no contiene ninguna célula o una célula basal, pasamos a la siguiente
	if(cell_is_null(current_cell_coord) == true || get_cell_type(current_cell_coord) == BASAL)
		return;

	uint8_t current_cell_type = get_cell_type(current_cell_coord);
	bool end_flag = true; // Indica si la inicialización del ducto ha acabado
	bool reproduction_flag = false;

	if(!initialized)
	{
		// 4. Allow PROGENITOR cells to REPRODUCE into VACANT NEIGHBORING lattice points WITHOUT CELLS
		//    according to progenitor hierarchy
		if (current_cell_type == PROG_STEM     ||
			current_cell_type == PROG_BIPOTENT ||
			current_cell_type == PROG_MYOEPI   ||
			current_cell_type == PROG_LUMEPI)
		{
			if (reproduce(current_cell_coord, current_cell_type) == true)
			{
				end_flag = false;
				reproduction_flag = true;
			}
		}

		// 5. Allow NON-STEM cells to MIGRATE 1 point to a VACANT NEIGHBORING lattice point consistent
		//    with the bi-layer ductal structure
		if(current_cell_type != PROG_STEM)
		{
			if (migrate(current_cell_coord, current_cell_type) == true)
				end_flag = false;
		}

		if(end_flag == false)
			local_initialized.at(thread_index) = end_flag;
	}

	this->update_local_counters(current_cell_coord, false, false, reproduction_flag, false);
}



void Duct::update_local_counters(const Coordinate& current_cell_coord,
								 bool dysfunctional_flag,
								 bool cancerous_flag,
								 bool reproduction_flag,
								 bool cancerous_reproduction_flag)
{
	total_population++;
	if(dysfunctional_flag)
		dysfunctional_population++;
	if (reproduction_flag)
	{
		total_reproductions++;
		if (cancerous_reproduction_flag)
		{
			dysfunctional_reproductions++;
		}
	}
	// Para Graph1
	if(!DCIS_achieved && cancerous_flag)
	{
		DCIS_achieved = true;
		DCIS_generation = generation_number;
	}
	update_cell_counters(current_cell_coord);
}

/**
 * @return true si la coordenada c se encuentra dentro del array ducto
 *		   (no necesariamente dentro de la membrana basal)
 */
bool Duct::coordinate_is_inbounds(const Coordinate& c)
{
	return (c.x >= 0) && (c.x < get_diameter()) &&
		   (c.y >= 0) && (c.y < get_diameter()) &&
		   (c.z >= 0) && (c.z < get_length());
}

void Duct::get_neighbors(const Coordinate& c, vector<Coordinate>& neighbors)
{
	int x = c.x;
	int y = c.y;
	int z = c.z;
	// Las 6 coordenadas "principales" (vecindad directa)
	neighbors.at(0) = (Coordinate(x,   y,   z-1));
	neighbors.at(1) = (Coordinate(x,   y,   z+1));
	neighbors.at(2) = (Coordinate(x,   y-1, z));
	neighbors.at(3) = (Coordinate(x-1, y,   z));
	neighbors.at(4) = (Coordinate(x+1, y,   z));
	neighbors.at(5) = (Coordinate(x,   y+1, z));
	// Vecinos diagonales
	neighbors.at(6) = (Coordinate(x-1, y-1, z-1)); // diagonal
	neighbors.at(7) = (Coordinate(x,   y-1, z-1)); // diagonal
	neighbors.at(8) = (Coordinate(x+1, y-1, z-1)); // diagonal
	neighbors.at(9) = (Coordinate(x-1, y,   z-1)); // diagonal
	neighbors.at(10) = (Coordinate(x+1, y,   z-1)); // diagonal
	neighbors.at(11) = (Coordinate(x-1, y+1, z-1)); // diagonal
	neighbors.at(12) = (Coordinate(x,   y+1, z-1)); // diagonal
	neighbors.at(13) = (Coordinate(x+1, y+1, z-1)); // diagonal
	neighbors.at(14) = (Coordinate(x-1, y-1, z+1)); // diagonal
	neighbors.at(15) = (Coordinate(x,   y-1, z+1)); // diagonal
	neighbors.at(16) = (Coordinate(x+1, y-1, z+1)); // diagonal
	neighbors.at(17) = (Coordinate(x-1, y,   z+1)); // diagonal
	neighbors.at(18) = (Coordinate(x+1, y,   z+1)); // diagonal
	neighbors.at(19) = (Coordinate(x-1, y+1, z+1)); // diagonal
	neighbors.at(20) = (Coordinate(x,   y+1, z+1)); // diagonal
	neighbors.at(21) = (Coordinate(x+1, y+1, z+1)); // diagonal
	neighbors.at(22) = (Coordinate(x-1, y-1, z));   // diagonal
	neighbors.at(23) = (Coordinate(x+1, y-1, z));   // diagonal
	neighbors.at(24) = (Coordinate(x-1, y+1, z));   // diagonal
	neighbors.at(25) = (Coordinate(x+1, y+1, z));   // diagonal
}

/**
 * @return Devuelve una coordenada aleatoria que cumple las siguientes condiciones:
 *			- Es vecina de cell_coord
 *			- Se encuentra dentro de la membrana basal
 *			- Está vacía (no está ocupada por una célula)
 *		   Si tal coordenada no existe, devuelve null
 */
Coordinate Duct::get_vacant_neighbor(const Coordinate& cell_coord)
{
	vector<Coordinate> neighbors(26);
	get_neighbors(cell_coord, neighbors);

	shuffle(begin(neighbors), end(neighbors), rng);

	for(uint32_t i=0 ; i < neighbors.size() ; i++)
	{
		Coordinate aux_coord = neighbors.at(i);
		if(coordinate_is_inbounds(aux_coord) == true &&
		   cell_is_null(aux_coord) == true &&
		   get_point_radius(aux_coord.x, aux_coord.y) < radius)
		{
			return neighbors.at(i);
		}
	}
	return null_coordinate;
}

/**
 * @param cell_coord Coordenada a la que se va a migrar
 * @param old_coord  Coordenada desde la que se migra
 * @return true si cell_coord tiene alguna célula vecina
 */
bool Duct::has_adjadcent_neighbor(const Coordinate& cell_coord, const Coordinate& old_coord)
{
	vector<Coordinate> neighbors(26);
	get_neighbors(cell_coord, neighbors);
	for(uint32_t i=0 ; i < neighbors.size() ; i++)
	{
		Coordinate aux_coord = neighbors.at(i);
		if(coordinate_is_inbounds(aux_coord) == true &&
		   cell_is_null(aux_coord) == false /*&&
		   !aux_coord.equals(old_coord)*/)
		{
			return true;
		}
	}
	return false;
}

bool Duct::reproduce(const Coordinate& current_cell_coord, uint8_t current_cell_type)
{
	vector<Coordinate> neighbors(26);
	get_neighbors(current_cell_coord, neighbors);
	shuffle(begin(neighbors), end(neighbors), rng);


	/*
	INIT:	Allow progenitor cells to reproduce into vacant neighboring lattice points
			without cells according to progenitor hierarchy
	NORMAL: if an adjacent neighboring point is vacant and does not break the bi-layer
			ductal structure (Fig. 3C), the cell reproduces
	*/

	/* [ALT]: (a) Permitir que una célula mioepitelial se reproduzca en el sitio que sea, y luego
				  se corrija ésto en la migración, o bien
			  (b) Permitir que sólo se genere una célula en un sitio que le 'corresponde' [actual]
	*/

	int max_radius = radius-1;
	int min_radius = radius-2;

	if(current_cell_type == PROG_LUMEPI || current_cell_type == LUMEPI)
		max_radius--;
	else if(current_cell_type == PROG_MYOEPI || current_cell_type == MYOEPI)
		min_radius++;

	Coordinate new_cell_coord;
	bool end_flag = false;
	for(uint32_t i=0 ; i<neighbors.size() ; i++)
	{
		new_cell_coord = neighbors.at(i);
		if(coordinate_is_inbounds(new_cell_coord) == true &&
		   cell_is_null(new_cell_coord) == true &&
		   get_point_radius(new_cell_coord.x, new_cell_coord.y) <= max_radius &&
		   get_point_radius(new_cell_coord.x, new_cell_coord.y) >= min_radius)
		{
			end_flag = true;
			break;
		}
	}

	if(end_flag == false)
	{
		return false;
	}

	set_new_cell_reproduce(current_cell_coord, new_cell_coord, current_cell_type);

	return true;
}

bool Duct::cancerous_reproduce(const Coordinate& current_cell_coord, uint8_t current_cell_type)
{
	Coordinate new_cell_coord;
	Coordinate pushed_position = null_coordinate;

	vector<Coordinate> neighbors(26);
	get_neighbors(current_cell_coord, neighbors);
	shuffle(begin(neighbors), end(neighbors), rng);

	bool end_flag = false;
	for (uint32_t i=0 ; i<neighbors.size() ; i++)
	{
		new_cell_coord = neighbors.at(i);

		if (coordinate_is_inbounds(new_cell_coord) == true)
		{
			if (cell_is_null(new_cell_coord) == true &&
				get_point_radius(new_cell_coord.x, new_cell_coord.y) < radius &&
				has_adjadcent_neighbor(new_cell_coord, current_cell_coord))
			{
				end_flag = true;
				break;
			}
			// [ALT] Sólo empujamos células LUMEPI (?)
			if ((pushed_position = get_vacant_neighbor(new_cell_coord)) != null_coordinate)
			{
				end_flag = true;
				break;
			}
		}
	}

	if(end_flag == false)
	{
		return false;
	}

	if(pushed_position != null_coordinate)
	{
		push_cell(new_cell_coord, pushed_position);
	}

	set_new_cell_reproduce(current_cell_coord, new_cell_coord, current_cell_type);

	return true;
}

bool Duct::hormonal_reproduce(const Coordinate& current_cell_coord, uint8_t current_cell_type)
{
	vector<Coordinate> neighbors(26);
	get_neighbors(current_cell_coord, neighbors);
	shuffle(begin(neighbors), end(neighbors), rng);

	Coordinate new_cell_coord = null_coordinate;

	bool end_flag = false;
	Coordinate pushed_position = null_coordinate;
	for (uint32_t i=0 ; i<neighbors.size() ; i++)
	{
		new_cell_coord = neighbors.at(i);
		if (coordinate_is_inbounds(new_cell_coord) != false &&
			get_point_radius(new_cell_coord.x, new_cell_coord.y) < radius &&
			cell_is_null(new_cell_coord) == false &&
			(pushed_position = get_vacant_neighbor(new_cell_coord)) != null_coordinate
		   )
		{
			end_flag = true;
			break;
		}
	}

	if(end_flag == false)
	{
		return false;
	}

	push_cell(new_cell_coord, pushed_position);

	set_new_cell_reproduce(current_cell_coord, new_cell_coord, current_cell_type);

	return true;
}

bool Duct::migrate(const Coordinate& current_cell_coord, uint8_t current_cell_type)
{
	// Cortocircuitar cuanto antes
	if(get_point_radius(current_cell_coord.x, current_cell_coord.y) < radius - 3)
	{
		return false;
	}

	vector<Coordinate> neighbors(26);
	get_neighbors(current_cell_coord, neighbors);
	shuffle(begin(neighbors), end(neighbors), rng);

	int max_radius = radius-1;
	int min_radius = radius-2;

	if(current_cell_type == PROG_LUMEPI || current_cell_type == LUMEPI)
		max_radius--;
	else if(current_cell_type == PROG_MYOEPI ||current_cell_type == MYOEPI)
		min_radius++;

	Coordinate new_cell_coord = null_coordinate;
	bool end_flag = false;
	for (uint32_t i=0 ; i<neighbors.size() ; i++)
	{
		new_cell_coord = neighbors.at(i);
		if (coordinate_is_inbounds(new_cell_coord) == true &&
			cell_is_null(new_cell_coord) == true &&
			get_point_radius(new_cell_coord.x, new_cell_coord.y) <= max_radius &&
			get_point_radius(new_cell_coord.x, new_cell_coord.y) >= min_radius)

		{
			end_flag = true;
			break;
		}

	}

	if(end_flag == false)
	{
		return false;
	}

	push_cell(current_cell_coord, new_cell_coord);

	return true;
}

void Duct::update_global_counters()
{
	// Actualizar contadores
	total_population_global += total_population;
	dysfunctional_population_global += dysfunctional_population;
	total_reproductions_global += total_reproductions;
	dysfunctional_reproductions_global += dysfunctional_reproductions;
	stem_counter_global += stem_counter;
	bipotent_counter_global += bipotent_counter;
	proglumepi_counter_global += proglumepi_counter;
	lumepi_counter_global += lumepi_counter;
	mut_stem_counter_global += mut_stem_counter;
	mut_bipotent_counter_global += mut_bipotent_counter;
	mut_proglumepi_counter_global += mut_proglumepi_counter;
	mut_lumepi_counter_global += mut_lumepi_counter;

	if(++thread_counter == n_threads)
	{
		total_history.push_back(total_population_global.load());
		dysfunctional_history.push_back(dysfunctional_population_global.load());
		reproductions_history.push_back(
			total_reproductions_global.load() + reproductions_history[reproductions_history.size()-1]);
		dysf_reproductions_history.push_back(
			dysfunctional_reproductions_global.load() + dysf_reproductions_history[dysf_reproductions_history.size()-1]);

		stem_mean = (float)mut_stem_counter_global.load()/(float)stem_counter_global.load();
		bipotent_mean = (float)mut_bipotent_counter_global.load()/(float)bipotent_counter_global.load();
		proglumepi_mean = (float)mut_proglumepi_counter_global.load()/(float)proglumepi_counter_global.load();
		lumepi_mean = (float)mut_lumepi_counter_global.load()/(float)lumepi_counter_global.load();

		total_population_global = 0;
		dysfunctional_population_global = 0;
		total_reproductions_global = 0;
		dysfunctional_reproductions_global = 0;
		stem_counter_global = 0;
		bipotent_counter_global = 0;
		proglumepi_counter_global = 0;
		lumepi_counter_global = 0;
		mut_stem_counter_global = 0;
		mut_bipotent_counter_global = 0;
		mut_proglumepi_counter_global = 0;
		mut_lumepi_counter_global = 0;

		thread_counter = 0;
		generation_number++;
	}
}

vector<vector<vector<uint8_t> > > Duct::get_duct() { return duct_type; }
vector<vector<uint8_t> > Duct::get_slice(int z)    { return duct_type.at(z); }
uint8_t Duct::get_cell_type(const Coordinate& co)          { return    		 duct_type.at(co.z).at(co.y).at(co.x); }
uint8_t Duct::get_gene_housekeeping(const Coordinate& co)  { return  duct_housekeeping.at(co.z).at(co.y).at(co.x); }
uint8_t Duct::get_gene_protooncogene(const Coordinate& co) { return duct_protooncogene.at(co.z).at(co.y).at(co.x); }
uint8_t Duct::get_gene_supressor(const Coordinate& co)     { return     duct_supressor.at(co.z).at(co.y).at(co.x); }
uint8_t Duct::get_gene_apoptosis(const Coordinate& co)     { return     duct_apoptosis.at(co.z).at(co.y).at(co.x); }


void Duct::set_cell (const Coordinate& co, uint8_t t, uint8_t hk, uint8_t po, uint8_t sp, uint8_t ap)
{
	         duct_type.at(co.z).at(co.y).at(co.x) = t;
	 duct_housekeeping.at(co.z).at(co.y).at(co.x) = hk;
	duct_protooncogene.at(co.z).at(co.y).at(co.x) = po;
	    duct_supressor.at(co.z).at(co.y).at(co.x) = sp;
	    duct_apoptosis.at(co.z).at(co.y).at(co.x) = ap;
}

void Duct::erase_cell(const Coordinate& co)
{
	         duct_type.at(co.z).at(co.y).at(co.x) = 0;
	 duct_housekeeping.at(co.z).at(co.y).at(co.x) = 0;
	duct_protooncogene.at(co.z).at(co.y).at(co.x) = 0;
	    duct_supressor.at(co.z).at(co.y).at(co.x) = 0;
	    duct_apoptosis.at(co.z).at(co.y).at(co.x) = 0;
}

void Duct::set_cell_type(const Coordinate& co, uint8_t type)
	{ duct_type.at(co.z).at(co.y).at(co.x) = type; }

// Muta 10% del genoma
uint8_t Duct::mutate_HGP()
{
	uint8_t counter = (uint8_t)(((((float)gene_length) / 100) * 10));
	return counter;
}

// Mutation Rate representado como por mil, no por cien
uint8_t Duct::mutate(uint8_t gene, int mutation_rate)
{
	for(int i=0 ; i < (gene_length - gene) ; i++)
	{
		if(nextInt(1, 1000) <= mutation_rate && gene < 32)
			gene++;
	}
	return gene;
}

bool Duct::mark_cell_for_death(const Coordinate& co)
{
	return housekeeping_is_dysfunctional(co)
		   || (protooncogene_is_dysfunctional(co) && !supressor_is_dysfunctional(co) && !apoptosis_is_dysfunctional(co))
		   || (!apoptosis_is_dysfunctional(co) &&
	           get_point_radius(co.x, co.y) < radius-2);
}

// Método que usamos para ver si la célula es cancerígena según el paper
bool Duct::is_considered_cancerous(const Coordinate& co)
{
	return //housekeeping_is_dysfunctional(co) ||
		   protooncogene_is_dysfunctional(co)  ||
		   //supressor_is_dysfunctional(co)    ||
		   apoptosis_is_dysfunctional(co);
}

int Duct::get_mutations(const Coordinate& co)
{
	return get_gene_housekeeping(co)  +
		   get_gene_protooncogene(co) +
		   get_gene_supressor(co)     +
		   get_gene_apoptosis(co);
}

bool Duct::is_dysfunctional(const Coordinate& co)
{
	return get_gene_housekeeping(co)  >= 16 ||
		   get_gene_protooncogene(co) >= 16 ||
		   get_gene_supressor(co)	  >= 16  ||
		   get_gene_apoptosis(co)	  >= 16;
}

bool Duct::housekeeping_is_dysfunctional(const Coordinate& co)	 { return  duct_housekeeping.at(co.z).at(co.y).at(co.x) >= 16; }
bool Duct::protooncogene_is_dysfunctional(const Coordinate& co) { return duct_protooncogene.at(co.z).at(co.y).at(co.x) >= 16; }
bool Duct::supressor_is_dysfunctional(const Coordinate& co)     { return     duct_supressor.at(co.z).at(co.y).at(co.x) >= 16; }
bool Duct::apoptosis_is_dysfunctional(const Coordinate& co)     { return     duct_apoptosis.at(co.z).at(co.y).at(co.x) >= 16; }

bool Duct::cell_is_null(const Coordinate& co) { return get_cell_type(co) == 0; }

void Duct::set_new_cell_reproduce(const Coordinate& current_cell_coord, const Coordinate& new_cell_coord,
										  uint8_t current_cell_type)
{

	uint8_t new_cell_type = 0;
	switch (current_cell_type)
	{
		// [OPT] Creo que los switch se compilan como una cadena de if-else;
		// por lo que es mejor poner la opción más común primero
		case PROG_LUMEPI:
			new_cell_type = LUMEPI;
			break;
		case PROG_BIPOTENT:
			// [ALT] Decidimos qué tipo de célula va a generar una bipotente en base a
			// la posición en que se va a generar
			if (get_point_radius(new_cell_coord.x, new_cell_coord.y) == get_radius()-1)
				new_cell_type = PROG_MYOEPI;
			else
				new_cell_type = PROG_LUMEPI;
			break;
		case PROG_STEM:
			new_cell_type = PROG_BIPOTENT;
			break;
		case PROG_MYOEPI:
			new_cell_type = MYOEPI;
			break;
		default:
			break;
	}

	// Mutamos tanto madre como hija
	uint8_t parent_hk, parent_po, parent_sp, parent_ap;
	uint8_t new_cell_hk, new_cell_po, new_cell_sp, new_cell_ap;

	new_cell_hk = parent_hk = get_gene_housekeeping(current_cell_coord);
	new_cell_po = parent_po = get_gene_protooncogene(current_cell_coord);
	new_cell_sp = parent_sp = get_gene_supressor(current_cell_coord);
	new_cell_ap = parent_ap = get_gene_apoptosis(current_cell_coord);

	parent_hk   = mutate(parent_hk, mutation_rate);
	parent_po   = mutate(parent_po, mutation_rate);
	parent_sp   = mutate(parent_sp, mutation_rate);
	parent_ap   = mutate(parent_ap, mutation_rate);
	new_cell_hk = mutate(new_cell_hk, mutation_rate);
	new_cell_po = mutate(new_cell_po, mutation_rate);
	new_cell_sp = mutate(new_cell_sp, mutation_rate);
	new_cell_ap = mutate(new_cell_ap, mutation_rate);
	set_cell(current_cell_coord, current_cell_type, parent_hk, parent_po, parent_sp, parent_ap);
	set_cell(new_cell_coord, new_cell_type, new_cell_hk, new_cell_po, new_cell_sp, new_cell_ap);
}

void Duct::set_new_stem_cell(const Coordinate& new_coordinate)
{
	uint8_t new_stem_cell_type = PROG_STEM;
	uint8_t new_stem_cell_housekeeping  = 0;
	uint8_t new_stem_cell_protooncogene = 0;
	uint8_t new_stem_cell_supressor     = 0;
	uint8_t new_stem_cell_apoptosis     = 0;

	// 3. If HGP, mutate all genes in all stem cells by 10%
	if(HGP)
	{
		// [ALT]: Mutación determinista (mutar 10% del genoma) o estocástica (mutar con 10% de probabilidad) ?
		new_stem_cell_housekeeping  = mutate_HGP();
		new_stem_cell_protooncogene = mutate_HGP();
		new_stem_cell_supressor     = mutate_HGP();
		new_stem_cell_apoptosis     = mutate_HGP();
	}
	set_cell(new_coordinate, new_stem_cell_type,
			 new_stem_cell_housekeeping,
			 new_stem_cell_protooncogene,
			 new_stem_cell_supressor,
			 new_stem_cell_apoptosis);
}

void Duct::push_cell(const Coordinate& new_cell_coord, const Coordinate& pushed_position)
{
	set_cell(pushed_position, get_cell_type(new_cell_coord),
							  get_gene_housekeeping(new_cell_coord),
							  get_gene_protooncogene(new_cell_coord),
							  get_gene_supressor(new_cell_coord),
							  get_gene_apoptosis(new_cell_coord));
	erase_cell(new_cell_coord);
}

void Duct::normal_routine(const Coordinate& current_cell_coord)
{
	if(cell_is_null(current_cell_coord) == true
	   || get_cell_type(current_cell_coord) == BASAL)
	{
		return;
	}

	uint8_t current_cell_type = get_cell_type(current_cell_coord);

	// 1. If housekeeping or apoptosis genes are dysfunctional (50% mutated), mark the cell for death
	// 2. If apoptosis genes are functional and cell is not adjacent to the basal membrane
	//    or a myoepithelial cell, mark the cell for death
	// 3. If dead, remove from lattice
	if(mark_cell_for_death(current_cell_coord))
	{
		erase_cell(current_cell_coord);
		return;
	}

	bool reproduction_flag = false;
	bool cancerous_reproduction_flag = false;

	// 4. Progenitor cells may reproduce uint8_to neighboring positions according to
	//	  progenitor hierarchy for three reasons:
	if (current_cell_type == PROG_STEM     ||
		current_cell_type == PROG_BIPOTENT ||
		current_cell_type == PROG_MYOEPI   ||
		current_cell_type == PROG_LUMEPI)
	{
		// Cancerous reproduction—if the cell has a dysfunctional (50%) protooncogene
		// and tumor suppressor gene, cell reproduces uint8_to any vacant neighboring lattice
		// point or by pushing a neighboring cell to a vacant lattice point and reproducing
		if(protooncogene_is_dysfunctional(current_cell_coord) &&
		   supressor_is_dysfunctional(current_cell_coord))
		{
			if (cancerous_reproduce(current_cell_coord, current_cell_type))
			{
				reproduction_flag = true;
				if (is_considered_cancerous(current_cell_coord))
					cancerous_reproduction_flag = true;
			}
		}
		else
		{
			// Hormonal reproduction—if cell is stochastically chosen (1% chance), cell
			// pushes a adjcance neighboring cell to a vacant lattice point and reproduces
			// [ALT]: No permitimos rep. hormonal de MYOEPI para que no haya MYOEPIs en la capa LUMEPI
			if(current_cell_type != PROG_MYOEPI && nextInt(1, 100) == 1)
			{
				if (hormonal_reproduce(current_cell_coord, current_cell_type))
				{
					reproduction_flag = true;
					if (is_considered_cancerous(current_cell_coord))
						cancerous_reproduction_flag = true;
				}
			}
			else
			{
				// Normal reproduction—if an adjacent neighboring point is vacant and does
				// not break the bi-layer ductal structure (Fig. 3C), the cell reproduces
				if (reproduce(current_cell_coord, current_cell_type))
				{
					reproduction_flag = true;
					if (is_considered_cancerous(current_cell_coord))
						cancerous_reproduction_flag = true;
				}
			}
		}
	}

	// 5. Allow non-stem cells to migrate one point to a vacant neighboring lattice point
	//	  consistent with the bi-layer ductal structure (Fig. 3C)
	if(current_cell_type != PROG_STEM)
	{
		migrate(current_cell_coord, current_cell_type);
	}

	update_local_counters(current_cell_coord,
                                is_dysfunctional(current_cell_coord),
                                is_considered_cancerous(current_cell_coord),
                                reproduction_flag,
                                cancerous_reproduction_flag);
}

void Duct::update_cell_counters(const Coordinate& co)
{
    uint8_t current_cell_type = get_cell_type(co);
	switch(current_cell_type)
	{
		case PROG_STEM:
			stem_counter++;
			mut_stem_counter += get_mutations(co);
			break;
		case PROG_BIPOTENT:
			bipotent_counter++;
			mut_bipotent_counter += get_mutations(co);
			break;
		case PROG_LUMEPI:
			proglumepi_counter++;
			mut_proglumepi_counter += get_mutations(co);
			break;
		case LUMEPI:
			lumepi_counter++;
			mut_lumepi_counter += get_mutations(co);
			break;
	}
}

Duct::~Duct()
{
	if(n_threads > 1)
		delete thread_barrier;
}
