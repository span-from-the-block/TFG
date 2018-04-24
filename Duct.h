#ifndef DUCT_H_
#define DUCT_H_

#include <cstdint>
#include <tuple>
#include <atomic>
#include <thread>
#include <cstdint>
#include <random>
#include <vector>
#include "cyclic_barrier.h"

using namespace std;

//             x    y    z
//typedef tuple<int, int, int> Coordinate;
struct Coordinate;

class Duct {

	public: 
		static const uint8_t BASAL 		   = 1;
		static const uint8_t PROG_STEM	   = 2;
		static const uint8_t PROG_BIPOTENT = 3;
		static const uint8_t PROG_MYOEPI   = 4;
		static const uint8_t PROG_LUMEPI   = 5;
		static const uint8_t MYOEPI        = 6;
		static const uint8_t LUMEPI        = 7;

		static const uint8_t gene_length   = 32;
		static const int     mutation_rate = 5; // Representado como mutation_rate/1000

		vector<vector<vector<uint8_t> > > duct_type;
		vector<vector<vector<uint8_t> > > duct_housekeeping;
		vector<vector<vector<uint8_t> > > duct_protooncogene;
		vector<vector<vector<uint8_t> > > duct_supressor;
		vector<vector<vector<uint8_t> > > duct_apoptosis;

		int	length;
		int	radius;
		bool HGP; 		   // Indica presencia o no de HGP
		int n_threads; // nº de hilos
		int	diameter;
		int	capacity;	   // Capacidad total del ducto
		int	stem_capacity; // Capacidad de las dos capas principales del ducto

		int	generation_number;
		bool initialized; 	  // Indica si se ha inicializado el ducto
		vector<bool> local_initialized;

		// Contadores (para GUI)
		// Contadores de población y reproducción 
		// Contadores locales
		static thread_local int total_population;
		static thread_local int dysfunctional_population;
		static thread_local int total_reproductions;
		static thread_local int dysfunctional_reproductions;
		// Contadores globales
		atomic<int> total_population_global;
		atomic<int> dysfunctional_population_global;
		atomic<int> total_reproductions_global;
		atomic<int> dysfunctional_reproductions_global;

		vector<int> total_history;
		vector<int> dysfunctional_history;
		vector<int> reproductions_history;
		vector<int> dysf_reproductions_history;

		// Contadores de mutación
		// Contadores locales
		static thread_local int stem_counter; 
		static thread_local int bipotent_counter;
		static thread_local int proglumepi_counter;
		static thread_local int lumepi_counter;
		static thread_local int mut_stem_counter;
		static thread_local int mut_bipotent_counter;
		static thread_local int mut_proglumepi_counter;
		static thread_local int mut_lumepi_counter;
	    // Contadores globales
	    atomic<int> stem_counter_global; 
		atomic<int> bipotent_counter_global;
		atomic<int> proglumepi_counter_global;
		atomic<int> lumepi_counter_global;
		atomic<int> mut_stem_counter_global;
		atomic<int> mut_bipotent_counter_global;
		atomic<int> mut_proglumepi_counter_global;
		atomic<int> mut_lumepi_counter_global;
	    float stem_mean;
		float bipotent_mean;
		float proglumepi_mean;
		float lumepi_mean;

		bool DCIS_achieved;   // Indica si se ha alcanzado DCIS o no
		int  DCIS_generation; // Generación en que se generó DCIS

		int   	   total_generations;

		CyclicBarrier* thread_barrier;
		vector<mutex>  thread_locks;
		atomic<int>	   thread_counter;
		vector<thread> thread_array;

	    static thread_local int         thread_index;
    	static thread_local vector<int> index_array;
    	static thread_local Coordinate  boundaries;

		int  get_generation_number();
		int  get_length();
		int  get_radius();
		int  get_diameter();
		int  get_capacity();
		bool get_HGP();
		int  is_DCIS_achieved();

		Duct(int, int, bool, int); 
		void run(int, int, int);
		void execute(int);
		void next_generation();
		void init_routine_stem();
		void init_routine(const Coordinate&);
		void update_local_counters(const Coordinate&, bool, bool, bool, bool);
		
		
		int get_point_radius(int, int);
		bool coordinate_is_inbounds(const Coordinate&);
		void get_neighbors(const Coordinate&, vector<Coordinate>&);
		Coordinate get_vacant_neighbor(const Coordinate&);
		bool has_adjadcent_neighbor(const Coordinate&, const Coordinate&);
		bool reproduce(const Coordinate&, uint8_t);
		bool cancerous_reproduce(const Coordinate&, uint8_t);
		bool hormonal_reproduce(const Coordinate&, uint8_t);
		bool migrate(const Coordinate&, uint8_t);
		void update_global_counters();
		float nextInt(int, int);

		void set_Duct_Type();
		vector<vector<vector<uint8_t> > > get_duct();
		vector<vector<uint8_t> > get_slice(int);
		uint8_t get_cell_type(const Coordinate&);
		uint8_t get_gene_housekeeping(const Coordinate&);
		uint8_t get_gene_protooncogene(const Coordinate&);
		uint8_t get_gene_supressor(const Coordinate&);
		uint8_t get_gene_apoptosis(const Coordinate&);


		void set_cell(const Coordinate&, uint8_t, uint8_t, uint8_t, uint8_t, uint8_t);
		void erase_cell(const Coordinate&);
		void set_cell_type(const Coordinate&, uint8_t);
		uint8_t mutate_HGP();
		uint8_t mutate(uint8_t, int);
		bool mark_cell_for_death(const Coordinate&);
		bool is_considered_cancerous(const Coordinate&);
		int get_mutations(const Coordinate&);
		bool is_dysfunctional(const Coordinate&);
		bool housekeeping_is_dysfunctional(const Coordinate&);
		bool protooncogene_is_dysfunctional(const Coordinate&);
		bool supressor_is_dysfunctional(const Coordinate&);
		bool apoptosis_is_dysfunctional(const Coordinate&);

		bool cell_is_null(const Coordinate&);

		void set_new_cell_reproduce(const Coordinate&, const Coordinate&, uint8_t);
		void set_new_stem_cell(const Coordinate&);
		void push_cell(const Coordinate&, const Coordinate&);
		void normal_routine(const Coordinate&);
		void update_cell_counters(const Coordinate&);
		~Duct();
};

#endif