/*
 * MultiHashFunction.h
 *
 *	Modified by: Leonardo
 */

#ifndef INCLUDE_HASH_MULTIHASHFUNCTION_H_
	#define INCLUDE_HASH_MULTIHASHFUNCTION_H_

	#include "Hash/HashFunction.h"
	#include <omp.h>




	namespace _hash_
	{
		inline static void GetHashes_speedup_multi_previous_Rotated(	const std::string& s_Str,
																		const SpacedQmer_Multi& spaced_qmers,
																		Hash_V_V& hashes)
		{
			auto get_hash = [&](size_t curr_spaced, size_t curr_idx_hash, const PreviusShift_Ext& curr_shift)
			{
				Hash& curr_hash = hashes[curr_spaced][curr_idx_hash];

				if(curr_shift.current_sp_ptr->GetWeight() < curr_shift.GetSize())
					GetHash(s_Str, curr_idx_hash, *curr_shift.current_sp_ptr, curr_hash);
				else
				{
					const Hash& prev_hash = hashes[curr_shift.prev_qmer][curr_idx_hash - curr_shift.shift_min];
					compute_hash_for_speedup_previous(	s_Str, curr_shift,
														curr_shift.current_sp_ptr->GetPosOne(), curr_shift.prev_sp_ptr->GetPosOne(),
														prev_hash, curr_hash, curr_idx_hash);
				}
			};


			hashes.resize(spaced_qmers.size());


			if (s_Str.size() < spaced_qmers[0].GetQ()) // Almeno 1 hash
			{
				for (size_t i = 0; i < spaced_qmers.size(); i++)
					hashes[i].clear();

				return;
			}


			const size_t n_hashes = s_Str.size() - spaced_qmers[0].GetQ() + 1;
			for (size_t i = 0; i < spaced_qmers.size(); i++)
				hashes[i].resize(n_hashes);


			const std::vector<PreviusShift_Ext_V>& v_shift_min = spaced_qmers.getShiftMinRotated();
			const PreviusShift_Ext_V& zero_prev = v_shift_min[0];
			for (size_t s = 0; s < spaced_qmers.size(); s++)
			{
				const PreviusShift_Ext& zero = zero_prev[s];
				if (!zero.isCorrectSpacedPrevious())
					GetHash(s_Str, 0, *zero.current_sp_ptr, hashes[s][0]); // Il primo è da computare a parte
				else
					get_hash(s, 0, zero);
			}

			size_t i = 1;
			for (; i < std::min(v_shift_min.size(), n_hashes); i++)
			{
				const PreviusShift_Ext_V& prev_i = v_shift_min[i];
				for (size_t ss = 0; ss < spaced_qmers.size(); ss++)
					get_hash(ss, i, prev_i[ss]);
			}

			const PreviusShift_Ext_V& prev_i = v_shift_min.back();
			for (; i < n_hashes; i++)
				for (size_t ss = 0; ss < spaced_qmers.size(); ss++)
					get_hash(ss, i, prev_i[ss]);
		}




		// ------------------------------------------------- //
		//                                                   //
		//            I S S H   M U L T I S E E D            //
		//                                                   //
		// ------------------------------------------------- //
		inline static void compute_hash_with_ISSH_multi_v1(	const std::string& s_Str, const std::vector<V_V_PreviusShift>& VV_shifts,
															const std::vector<Position>& v_pos_one,
															Hash_V_V& hashes, size_t& curr_idx_hash)
		{
			// At this point the transitory part is done, in each position it is possible to retrive all the 
			// positions, but the last one.
			// We know the last position and therefore the corresponding value inside the DNA sequence.
			const Position& pos_not_covered_yet = VV_shifts[0][0][0].one_to_change;
			const size_t& i_to_change = pos_not_covered_yet[0];
			
			for (size_t j = 0; j < VV_shifts.size(); ++j)
			{
				Hash& curr_hash = hashes[j][curr_idx_hash];
				curr_hash.encoding = 0;
				
				
				for (const PreviousShift& curr_sp_shift : VV_shifts[j][0])
				{
					uint64_t partial = hashes[j][curr_idx_hash - curr_sp_shift.shift_min].encoding;
					partial >>= (2 * curr_sp_shift.one_exit);
					partial &= curr_sp_shift.mask;
					curr_hash.encoding |= partial;
				}


				// identify the position of the unknown character inside the DNA sequence
				const unsigned char& tmp_char = s_Str[curr_idx_hash + v_pos_one[0][i_to_change]];
				curr_hash.encoding |= fEncodingTable[(tmp_char << 5) + i_to_change];
			}
		}


		inline static void GetHashes_with_ISSH_multi_v1(    const std::string& s_Str, const std::vector<SpacedQmer>& spaced_qmers,
															const std::vector<V_V_PreviusShift>& VV_shifts,
															const std::vector<Position>& v_pos_one,
															const size_t& max_transient_length,
															Hash_V_V& hashes)
		{
			hashes.resize(spaced_qmers.size());


			if (s_Str.size() < spaced_qmers[0].GetQ()) // Almeno 1 hash
			{
				for (size_t i = 0; i < spaced_qmers.size(); i++)
					hashes[i].clear();

				return;
			}


			const size_t n_hashes = s_Str.size() - spaced_qmers[0].GetQ() + 1;
			for (size_t i = 0; i < spaced_qmers.size(); i++)
				hashes[i].resize(n_hashes);
			

			//const size_t partial_limit = std::min(max_transient_length, n_hashes);
			size_t i = 1;
			
			// Per ogni spaced-seed, inizializza la struttura del risultato e calcola la parte del transiente
			for (size_t j = 0; j < spaced_qmers.size(); j++)
			{
				// Carico il vettore dei gruppi di shift precedenti che mi serve per
				// gestire sia il caso a regime sia tutto il transitorio.
				const V_V_PreviusShift& curr_seed = VV_shifts[j];
				
				// Per lo 0-esimo hash, ovviamente, non posso recuperare niente: lo calcolo posizione per posizione.
				GetHash(s_Str, 0, spaced_qmers[j].GetPosOne(), hashes[j][0]);
				
				i = 1;
				const size_t curr_max = curr_seed.size();
				for (; i < std::min(max_transient_length, n_hashes); i++)
				{
					// Se il vettore di shift è vuoto è sicuro che non posso recuperare nessuna posizione:
					// per evitare controlli e passaggi inutili si utilizza la funzione posizione per posizione.
					if (i < curr_max && curr_seed[i].size() == 0)
						GetHash(s_Str, i, spaced_qmers[j].GetPosOne(), hashes[j][i]);
					else
					{
						const size_t shift_group = (i < curr_max) ? i : 0;
						compute_hash_with_ISSH(s_Str, curr_seed[shift_group], v_pos_one[j], hashes[j], i);
					}
				}
			}

			// once the transient part is computed the standard computation can start for all the remaining columns
			// compute the hashes at position i for all the seeds 
			for (; i < n_hashes; i++)
				compute_hash_with_ISSH_multi_v1(s_Str, VV_shifts, v_pos_one, hashes, i);
		}




		inline static void compute_hash_with_ISSH_multi_col(	const std::string& s_Str, 
																const MultiSeedInfo& VV_shifts,
																Hash_V_V& hashes, size_t& curr_idx_hash)
		{
			// find what is the current situation, if the computation is standard or it is a transient step
			// the seed length is considered by 1+last position in which a one is present in the seed
			// since n_hashes = sequence_length - seed_length + 1 we can compute it as
			
			int shift_group = (curr_idx_hash < VV_shifts[0].group_previous.size() - 1) ? curr_idx_hash + 1 : 0;
			
			// for each of the seed in the group	
			for (size_t j = 0; j < VV_shifts.size(); j++)
			{
				// Get the reference to the current hash to compute
				Hash& curr_hash = hashes[j][curr_idx_hash];
				curr_hash.encoding = 0;

				const groupPrevious& curr_group = VV_shifts[j].group_previous[shift_group];
				const V_PreviusShiftMulti& curr_group_prev = curr_group.prev;
				
				// for faster performance the first hash is computed by calling GetHash instead of recovering all the positions later
				if (VV_shifts[j].group_previous[shift_group].prev.size() == 0)
					GetHash(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, hashes[j][curr_idx_hash]);
				else
				{
					for (const PreviousShiftMulti& curr_sp_shift : curr_group_prev)
					{
						uint64_t partial = hashes[curr_sp_shift.seed_num][curr_idx_hash - curr_sp_shift.shift_min].encoding;
						partial = (curr_sp_shift.one_exit >= 0) ? partial >> (2 * curr_sp_shift.one_exit) : partial << (-2 * curr_sp_shift.one_exit);
						partial &= curr_sp_shift.mask;
						curr_hash.encoding |= partial;
					}

					// compute from the sequence all the missing positions
					const Position& pos_not_covered_yet = curr_group.not_covered;
					const Position& pos_one = VV_shifts[j].pos_ones;
					for (const size_t& i_to_change : pos_not_covered_yet)
					{
						const unsigned char& tmp_char = s_Str[curr_idx_hash + pos_one[i_to_change]];
						curr_hash.encoding |= fEncodingTable[(tmp_char << 5) + i_to_change];
					}
				}
			}
		}


		inline static void GetHashes_with_ISSH_multi_col(const std::string& s_Str, const MultiSeedInfo& VV_shifts, Hash_V_V& hashes)
		{
			hashes.resize(VV_shifts.size());


			if (s_Str.size() + 1 < VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size() - 1]) // Almeno 1 hash
			{
				for (size_t i = 0; i < VV_shifts.size(); i++)
					hashes[i].clear();

				return;
			}


			const size_t n_hashes = s_Str.size() - VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size() - 1];
			for (size_t i = 0; i < VV_shifts.size(); i++)
				hashes[i].resize(n_hashes);

			
			// compute a column at a time
			for (size_t i = 0; i < n_hashes; i++)
				compute_hash_with_ISSH_multi_col(s_Str, VV_shifts, hashes, i);
		}




		inline static void compute_hash_with_ISSH_multi_col_parallel(	const std::string& s_Str, 
																		const MultiSeedInfo& VV_shifts, 
																		Hash_V_V& hashes, size_t& curr_idx_hash, size_t& offset)
		{
			// find what is the current situation, if the computation is standard or it is a transient step
			int shift_group = (offset < VV_shifts[0].group_previous.size() - 1) ? offset + 1 : 0;
			
			for (size_t j = 0; j < VV_shifts.size(); j++)
			{
				// Get the reference to the current hash to compute
				Hash& curr_hash = hashes[j][curr_idx_hash];
				curr_hash.encoding = 0;

				const groupPrevious& curr_group = VV_shifts[j].group_previous[shift_group];
				const V_PreviusShiftMulti& curr_group_prev = curr_group.prev;
				
				// for faster performance the first hash is computed by calling GetHash instead of recovering all the positions later
				if (VV_shifts[j].group_previous[shift_group].prev.size() == 0)
					GetHash(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, hashes[j][curr_idx_hash]);
				else
				{	
					for (const PreviousShiftMulti& curr_sp_shift : curr_group_prev)
					{
						uint64_t partial = hashes[curr_sp_shift.seed_num][curr_idx_hash - curr_sp_shift.shift_min].encoding;
						partial = (curr_sp_shift.one_exit >= 0) ? partial >> (2 * curr_sp_shift.one_exit) : partial << (-2 * curr_sp_shift.one_exit);
						partial &= curr_sp_shift.mask;
						curr_hash.encoding |= partial;
					}

					// compute from the seqence all the missing positions
					const Position& pos_not_covered_yet = curr_group.not_covered;
					const Position& pos_one = VV_shifts[j].pos_ones;
					for (const size_t& i_to_change : pos_not_covered_yet)
					{
						const unsigned char& tmp_char = s_Str[curr_idx_hash + pos_one[i_to_change]];
						curr_hash.encoding |= fEncodingTable[(tmp_char << 5) + i_to_change];
					}
				}
			}
		}


		inline static void GetHashes_with_ISSH_multi_col_parallel(const std::string& s_Str, const MultiSeedInfo& VV_shifts, Hash_V_V& hashes)
		{
			hashes.resize(VV_shifts.size());


			if (s_Str.size() + 1 < VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size() - 1]) // Almeno 1 hash
			{
				for (size_t i = 0; i < VV_shifts.size(); i++)
					hashes[i].clear();

				return;
			}


			const size_t n_hashes = s_Str.size() - VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size() - 1];	
			for (size_t i = 0; i < VV_shifts.size(); i++)
				hashes[i].resize(n_hashes);
			

			// compute a column at a time
			#pragma omp parallel
			{
				size_t offset = 0;

				#pragma omp for schedule(static)
				for (size_t i = 0; i < n_hashes; i++)
				{
					compute_hash_with_ISSH_multi_col_parallel(s_Str, VV_shifts, hashes, i, offset);
					offset++;
				}
			}
		}




		inline static void GetHashes_with_ISSH_multi_row(const std::string& s_Str, const MultiSeedInfoRow& rowInfo, Hash_V_V& hashes)
		{
			const MultiSeedInfo & VV_shifts = rowInfo.info;

			
			hashes.resize(VV_shifts.size());


			if (s_Str.size() + 1 < VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size() - 1]) // Almeno 1 hash
			{
				for (size_t i = 0; i < VV_shifts.size(); i++)
					hashes[i].clear();

				return;
			}


			const size_t n_hashes = s_Str.size() - VV_shifts[0].pos_ones[VV_shifts[0].pos_ones.size() - 1];
			for (size_t i = 0; i < VV_shifts.size(); i++)
				hashes[i].resize(n_hashes);


			size_t transient1 = rowInfo.transient1;
			size_t transient2 = rowInfo.transient2;
			if (n_hashes >= transient1 + transient2)
			{
				// initialize the matrix of hashes, one for each seed
				for (size_t j = 0; j < VV_shifts.size(); j++)
				{
					const V_groupPrevious & V_shift = VV_shifts[j].group_previous;
					
					// compute all the hashes for the current spaced seed
					for (size_t curr_idx_hash = 0; curr_idx_hash < n_hashes; curr_idx_hash++)
					{
						int shift_group = (curr_idx_hash >= transient1 && curr_idx_hash < n_hashes - transient2) ? 0 : ((curr_idx_hash < transient1) ? (curr_idx_hash + 1) : (V_shift.size() - (n_hashes - curr_idx_hash)));
						const groupPrevious& curr_group = V_shift[shift_group];
						const V_PreviusShiftMulti& curr_group_prev = curr_group.prev;
					
						// Get the reference to the current hash to compute
						Hash& curr_hash = hashes[j][curr_idx_hash];
						curr_hash.encoding = 0;
						
						// for faster performance the first hash is computed by calling GetHash instead of recovering all the positions later
						if (VV_shifts[j].group_previous[shift_group].prev.size() == 0)
							GetHash(s_Str, curr_idx_hash, VV_shifts[j].pos_ones, hashes[j][curr_idx_hash]);
						else
						{	
							// use previously computed hashes to recover the position
							for(const PreviousShiftMulti& curr_sp_shift : curr_group_prev)
							{	
								uint64_t partial = hashes[curr_sp_shift.seed_num][curr_idx_hash - curr_sp_shift.shift_min].encoding;
								partial = (curr_sp_shift.one_exit >= 0) ? partial >> (2 * curr_sp_shift.one_exit) : partial << (-2 * curr_sp_shift.one_exit);
								partial &= curr_sp_shift.mask;
								curr_hash.encoding |= partial;
							}

							// compute from the sequence all the missing positions
							const Position& pos_not_covered_yet = curr_group.not_covered;
							const Position& pos_one = VV_shifts[j].pos_ones;
							for(const size_t& i_to_change : pos_not_covered_yet)
							{
								const unsigned char& tmp_char = s_Str[curr_idx_hash + pos_one[i_to_change]];
								curr_hash.encoding |= fEncodingTable[(tmp_char << 5) + i_to_change];
							}
						}
					}
				}
			}
			else
				for(size_t j = 0; j < VV_shifts.size(); j++)
				{
					const Position& pos_one = VV_shifts[j].pos_ones;
					for(size_t curr_idx_hash = 0; curr_idx_hash < (size_t)n_hashes; curr_idx_hash++)
						GetHash(s_Str, curr_idx_hash, pos_one, hashes[j][curr_idx_hash]);
				}
		}
	}


#endif // INCLUDE_HASH_MULTIHASHFUNCTION_H_