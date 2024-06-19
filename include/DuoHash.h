/*
 * DuoHash.h
 *
 *  Created on: 18/jun/2024
 *      Author: Leonardo
 */

#ifndef INCLUDE_DUOHASH_H_
	#define INCLUDE_DUOHASH_H_

	#include "Hash/MultiHashFunction.h"

	#include <fstream>
	#include <string>
	#include <vector>

	
	
	
	bool loadFile(std::string fileName, std::vector<std::string>& lines)
	{
		std::cerr << "Load file... " << std::flush;

		std::ifstream fileStream(fileName);
		if (!fileStream.is_open())
		{
			std::cerr << "Unable to open file " << fileName << ".\n" << std::flush;
			return false;
		}

		std::string line;
		while (std::getline(fileStream, line))
			lines.push_back(line);

		fileStream.close();

		std::cerr << "Complete\n" << std::flush;
		return true;
	}




	// Classe per la gestione dei seed singoli
	class DuoHash
	{
		public:
			using PostProcessingFunction = void (*)(Hash_V& hashes, const size_t& k);


			DuoHash() {}
			DuoHash(const std::vector<std::string>& sequences, const SpacedQmer& spaced): sequences(sequences), spaced(spaced), k(spaced.GetWeight()) {}
			virtual ~DuoHash() {}


			void init(const std::vector<std::string>& sequences, const SpacedQmer& spaced)
			{
				this->sequences = sequences;
				this->spaced = spaced;
				this->k = this->spaced.GetWeight();
			}



			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//     G E T   E N C O D I N G   S I N G L E   S E E D
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			inline void GetEncoding_naive(std::vector<Hash_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());

				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_naive(this->sequences[seq], this->spaced, hashes[seq]);
			}
			inline void GetEncoding_naive(std::vector<Hash_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());

				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_naive(this->sequences[seq], this->spaced, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_FSH(std::vector<Hash_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_speedup_previous(this->sequences[seq], this->spaced, hashes[seq]);
			}
			inline void GetEncoding_FSH(std::vector<Hash_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_speedup_previous(this->sequences[seq], this->spaced, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_ISSH(std::vector<Hash_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_with_ISSH(this->sequences[seq], this->spaced, hashes[seq]);
			}
			inline void GetEncoding_ISSH(std::vector<Hash_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_with_ISSH(this->sequences[seq], this->spaced, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void PrintFASTA(const std::vector<Hash_V>& v_hashes)
			{
				for (const Hash_V& hashes : v_hashes)
					for (const Hash& hash : hashes)
						std::cout << ">\n" << hash.spacedKmer << "\n";
			}
			inline void PrintFASTA(const std::vector<Hash_V>& v_hashes, std::string fileName)
			{
				std::ofstream outfile = std::ofstream(fileName + ".fa");

				for (const Hash_V& hashes : v_hashes)
					for (const Hash& hash : hashes)
						outfile << ">\n" << hash.spacedKmer << "\n";
			}



			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//                M I S C E L L A N E O U S
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			inline const std::vector<std::string>& getSequences() const
			{
				return this->sequences;
			}

			inline const size_t getReads_avg() const
			{
				size_t read_avg = 0;
				for (size_t i = 0; i < this->sequences.size(); i++)
					read_avg += this->sequences[i].size();

				return read_avg / this->sequences.size();
			}

			inline const size_t getRead_min() const
			{
				size_t read_min = this->sequences[0].size();
				for (size_t i = 0; i < this->sequences.size(); i++)
					read_min = std::min(read_min, this->sequences[i].size());

				return read_min;
			}

			inline const size_t getRead_max() const
			{
				size_t read_max = this->sequences[0].size();
				for (size_t i = 0; i < this->sequences.size(); i++)
					read_max = std::max(read_max, this->sequences[i].size());

				return read_max;
			}

			inline const SpacedQmer& getSpacedQmer() const
			{
				return this->spaced;
			}

		private:
			std::vector<std::string> sequences;
			SpacedQmer spaced;
			size_t k;
	};


	// Classe per la gestione dei seed multipli
	class DuoHash_multi
	{
		public:
			using PostProcessingFunction = void (*)(Hash_V_V& hashes, const size_t& k);


			DuoHash_multi() {}
			DuoHash_multi(const std::vector<std::string>& sequences, const std::vector<SpacedQmer>& multi_spaced): sequences(sequences), multi_spaced(multi_spaced), k(multi_spaced[0].GetWeight())
			{
				this->spaced_qmers.init(this->multi_spaced);

				this->VV_shifts.resize(this->multi_spaced.size());
				this->v_pos_one.resize(this->multi_spaced.size());

				this->max_transient_length = 0;
				for(size_t j = 0; j < this->multi_spaced.size(); j++)
				{	
					this->VV_shifts[j] = this->multi_spaced[j].GetMultipleShifts();
					this->v_pos_one[j] = this->multi_spaced[j].GetPosOne();
					this->max_transient_length = std::max(this->VV_shifts[j].size(), this->max_transient_length);
				}

				MultiSpacedQmer MultiSeed(multi_spaced);
				this->infoCol = MultiSeed.Get_multi_seed_info_col();
				this->infoRow = MultiSeed.Get_multi_seed_info_row();
			}
			virtual ~DuoHash_multi() {}


			void init(const std::vector<std::string>& sequences, const std::vector<SpacedQmer>& multi_spaced)
			{
				this->sequences = sequences;
				this->multi_spaced = multi_spaced;
				this->k = this->multi_spaced[0].GetWeight();

				
				this->spaced_qmers.init(this->multi_spaced);

				this->VV_shifts.resize(this->multi_spaced.size());
				this->v_pos_one.resize(this->multi_spaced.size());

				this->max_transient_length = 0;
				for(size_t j = 0; j < this->multi_spaced.size(); j++)
				{	
					this->VV_shifts[j] = this->multi_spaced[j].GetMultipleShifts();
					this->v_pos_one[j] = this->multi_spaced[j].GetPosOne();
					this->max_transient_length = std::max(this->VV_shifts[j].size(), this->max_transient_length);
				}

				MultiSpacedQmer MultiSeed(multi_spaced);
				this->infoCol = MultiSeed.Get_multi_seed_info_col();
				this->infoRow = MultiSeed.Get_multi_seed_info_row();
			}



			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//      G E T   E N C O D I N G   M U L T I   S E E D
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			inline void GetEncoding_naive(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					hashes[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						_hash_::GetHashes_naive(this->sequences[seq], this->multi_spaced[ss], hashes[seq][ss]);
				}
			}
			inline void GetEncoding_naive(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					hashes[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						_hash_::GetHashes_naive(this->sequences[seq], this->multi_spaced[ss], hashes[seq][ss]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_FSH(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					hashes[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						_hash_::GetHashes_speedup_previous(this->sequences[seq], this->multi_spaced[ss], hashes[seq][ss]);
				}
			}
			inline void GetEncoding_FSH(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					hashes[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						_hash_::GetHashes_speedup_previous(this->sequences[seq], this->multi_spaced[ss], hashes[seq][ss]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_ISSH(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					hashes[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						_hash_::GetHashes_with_ISSH(this->sequences[seq], this->multi_spaced[ss], hashes[seq][ss]);
				}
			}
			inline void GetEncoding_ISSH(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					hashes[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						_hash_::GetHashes_with_ISSH(this->sequences[seq], this->multi_spaced[ss], hashes[seq][ss]);
					postProcessing(hashes[seq], this->k);
				}
			}


			inline void GetEncoding_FSH_multi(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_speedup_multi_previous_Rotated(this->sequences[seq], this->spaced_qmers, hashes[seq]);
			}
			inline void GetEncoding_FSH_multi(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_speedup_multi_previous_Rotated(this->sequences[seq], this->spaced_qmers, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_MISSH_v1(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_with_ISSH_multi_v1(this->sequences[seq], this->multi_spaced, this->VV_shifts, this->v_pos_one, this->max_transient_length, hashes[seq]);
			}
			inline void GetEncoding_MISSH_v1(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_with_ISSH_multi_v1(this->sequences[seq], this->multi_spaced, this->VV_shifts, this->v_pos_one, this->max_transient_length, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_MISSH_col(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_with_ISSH_multi_col(this->sequences[seq], this->infoCol, hashes[seq]);
			}
			inline void GetEncoding_MISSH_col(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_with_ISSH_multi_col(this->sequences[seq], this->infoCol, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_MISSH_col_parallel(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_with_ISSH_multi_col_parallel(this->sequences[seq], this->infoCol, hashes[seq]);
			}
			inline void GetEncoding_MISSH_col_parallel(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_with_ISSH_multi_col_parallel(this->sequences[seq], this->infoCol, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void GetEncoding_MISSH_row(std::vector<Hash_V_V>& hashes)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					_hash_::GetHashes_with_ISSH_multi_row(this->sequences[seq], this->infoRow, hashes[seq]);
			}
			inline void GetEncoding_MISSH_row(std::vector<Hash_V_V>& hashes, PostProcessingFunction postProcessing)
			{
				hashes.clear();
				hashes.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					_hash_::GetHashes_with_ISSH_multi_row(this->sequences[seq], this->infoRow, hashes[seq]);
					postProcessing(hashes[seq], this->k);
				}
			}

			inline void PrintFASTA(const std::vector<Hash_V_V>& v_hashes_v)
			{
				for (size_t ss = 0; ss < multi_spaced.size(); ss++)
					for (size_t seq = 0; seq < sequences.size(); seq++)
						for (size_t j = 0; j < sequences[seq].size() - multi_spaced[ss].GetQ() + 1; j++)
							std::cout << ">\n" << v_hashes_v[seq][ss][j].spacedKmer << "\n";
			}
			inline void PrintFASTA(const std::vector<Hash_V_V>& v_hashes_v, std::string fileName)
			{
				for (size_t ss = 0; ss < multi_spaced.size(); ss++)
				{
					std::ofstream outfile = std::ofstream(fileName + "_" + std::to_string(ss) + "_.fa");

					for (size_t seq = 0; seq < sequences.size(); seq++)
						for (size_t j = 0; j < sequences[seq].size() - multi_spaced[ss].GetQ() + 1; j++)
							outfile << ">\n" << v_hashes_v[seq][ss][j].spacedKmer << "\n";
				}
			}



			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//                M I S C E L L A N E O U S
			// +++++++++++++++++++++++++++++++++++++++++++++++++++++++
			inline const std::vector<std::string>& getSequences() const
			{
				return this->sequences;
			}

			inline const size_t getReads_avg() const
			{
				size_t read_avg = 0;
				for (size_t i = 0; i < this->sequences.size(); i++)
					read_avg += this->sequences[i].size();

				return read_avg / this->sequences.size();
			}

			inline const size_t getRead_min() const
			{
				size_t read_min = this->sequences[0].size();
				for (size_t i = 0; i < this->sequences.size(); i++)
					read_min = std::min(read_min, this->sequences[i].size());

				return read_min;
			}

			inline const size_t getRead_max() const
			{
				size_t read_max = this->sequences[0].size();
				for (size_t i = 0; i < this->sequences.size(); i++)
					read_max = std::max(read_max, this->sequences[i].size());

				return read_max;
			}

			inline const std::vector<SpacedQmer>& getSpacedQmers() const
			{
				return this->multi_spaced;
			}

		private:
			std::vector<std::string> sequences;
			std::vector<SpacedQmer> multi_spaced;
			size_t k;

			SpacedQmer_Multi spaced_qmers;

			std::vector<V_V_PreviusShift> VV_shifts;
			std::vector<Position> v_pos_one;
			size_t max_transient_length;

			MultiSeedInfo infoCol;
			MultiSeedInfoRow infoRow;
	};

#endif /* INCLUDE_DUOHASH_H_ */
