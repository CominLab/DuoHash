/*
 * MultiSpacedQmer.h
 *
 *       Author: enrico
 *	Modified by: Leonardo
 */


#ifndef INCLUDE_SPACED_MULTISPACEDQMER_H_
	#define INCLUDE_SPACED_MULTISPACEDQMER_H_

	#include "Spaced/SpacedQmer_Multi.h"
	
	
	// PreviousShiftMulti contains all the information needed to retrive 
	// some position from a hash that has already been computed
	struct PreviousShiftMulti {
		// Ones of the current seed that do not find a correspondence with the previous hash
		Position one_to_remove;
		
		// Ones of the current seed that find a correspondene with the previous hash
		Position one_to_keep;

		// one_exit represents by how much the previous hash has to be shifted in order to recreate the 
		// required overlapping. It can be negative, in that case it represent a shift in the opposite direction.
		int one_exit = 0;
		
		// shift_min represents the offset between the current position in the sequence and the one relative 
		// to the previous hash
		int shift_min = 0;
		
		// seed_num is the seed index, in the vector of seeds, from whose hashing we have to recover the positions
		int seed_num = 0;
		
		// mask needed to remove the unwantend values (one_to_remove) from the partial current hash that can be 
		// computed thanks to the information of the previous hash in order to be able to merge together more than
		// one partial hash 
		uint64_t mask = 0;
	};

	// The current hash is computed starting from a group of previous hashes that can be stored in a vector
	typedef std::vector<PreviousShiftMulti> V_PreviusShiftMulti;


	// In the struct groupPrevious are stored both the info about all the previous hashes used and which 
	// positions must be recomputed 
	struct groupPrevious
	{
		V_PreviusShiftMulti prev;

		// not_covered contains the positions that must be computed anew from the sequence because it is not possible to 
		// recover them from the previous hashes
		Position not_covered;
	};

	// V_groupPrevious contains all the information to manage both the standard case (in which the 
	// "best" group (found thanks our greedy approach) is all available) and the transient case 
	// (in which we must take in account that the availability of previous hashes is limited)
	typedef std::vector<groupPrevious> V_groupPrevious;


	struct SeedInfo
	{
		V_groupPrevious group_previous;
		// contains all the position in the seed in which a one is present. It is independent form the shift w.r.t. the sequence so it is stored only once for seed.  
		Position pos_ones;
	};

	// Contains all the computed information for each seed, for standard and each step of transient and of all the 
	// previous hashes for each case [seedindex][standard/transient_step].prev[previous_shift_index] 
	typedef std::vector<SeedInfo> MultiSeedInfo;


	// Struct that it is passed to the hashing function in case of hashing row ordered. In this case we need the additional information about
	// the initial and the final transient parts.
	struct MultiSeedInfoRow
	{
		MultiSeedInfo info;
		size_t transient1;
		size_t transient2;
	};


	// class that computes the information required for efficiently compute the hashing of DNA sequences by 
	// considering simultaneusly more than one spaced seed, all having the same weight and the same length
	// It also returns these information as a reference.
	class MultiSpacedQmer {
		public:
			MultiSpacedQmer() {}
			MultiSpacedQmer(const std::vector<SpacedQmer> &spaced_qmers);

			inline size_t GetLength() const
			{
				return this->spaced_qmers.size();
			}
			
			// gives the information required for computing the hashing by columns
			inline const MultiSeedInfo& Get_multi_seed_info_col() const
			{
				return this->multi_seed_info_col;
			}
			
			// gives the information required for computing the hashing by rows
			inline const MultiSeedInfoRow& Get_multi_seed_info_row() const
			{
				return this->multi_seed_info_row;
			}
				
			void reset(std::vector<SpacedQmer> spaced_qmers);

		private:
			// structure for each of single seeds.
			std::vector<SpacedQmer> spaced_qmers;
			std::vector<std::string> spaced_qs; // Vettore di stringhe rappresentante lo spaced seeds.
			std::vector<Position> v_pos_ones; // Vettore di vettori rappresentante le posizioni degli 1 in ogni seed
			std::vector<Position> v_pos_pos_ones; // Vettore di vettori rappresentanti le posizioni degli 1 in ogni seed rispetto al loro ordine.

			// final vector with all the needed information inside
			MultiSeedInfo multi_seed_info_col; // Per l'hashing column-wise
			MultiSeedInfoRow multi_seed_info_row; // Per l'hashing row-wise

			// compute info for each seed for every case (standard and transient steps)
			void SetMultiSeedInfoCol();
			// compute info for position i for each seed
			void ProcessMultiSeedCol(size_t i);

			// compute info for each seed for every case (standard or transient steps)
			void SetMultiSeedInfoRow();
			// compute info for position i for each seed
			void ProcessMultiSeedRow(MultiSeedInfo& multi_seed_info_row, int i);
	};


	// needed to visualize the information about a previous shift
	void print_shift_multi(PreviousShiftMulti s);

#endif /* INCLUDE_SPACED_MULTISPACEDQMER_H_ */
