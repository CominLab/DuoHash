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
			inline void GetEncoding_naive(std::vector<Encoding_V>& v_encodings)
			{
				v_encodings.clear();
				v_encodings.resize(this->sequences.size());

				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_naive(this->sequences[seq], this->spaced, v_encodings[seq]);
			}
			inline void GetEncoding_naive(std::vector<Encoding_V>& v_encodings, std::vector<Hashing_V>& v_hashings)
			{
				v_encodings.clear();
				v_hashings.clear();
				v_encodings.resize(this->sequences.size());
				v_hashings.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_naive(this->sequences[seq], this->spaced, v_encodings[seq]);
					getHashes(v_encodings[seq], this->k, v_hashings[seq]);
				}
			}
			inline void GetEncoding_naive(std::vector<Encoding_V>& v_encodings, std::vector<SpacedKmer_V>& v_spacedKmers)
			{
				v_encodings.clear();
				v_spacedKmers.clear();
				v_encodings.resize(this->sequences.size());
				v_spacedKmers.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_naive(this->sequences[seq], this->spaced, v_encodings[seq]);
					getSpacedKmers(v_encodings[seq], this->k, v_spacedKmers[seq]);
				}
			}
			inline void GetEncoding_naive(std::vector<Encoding_V>& v_encodings, std::vector<Hashing_V>& v_hashings, std::vector<SpacedKmer_V>& v_spacedKmers)
			{
				v_encodings.clear();
				v_hashings.clear();
				v_spacedKmers.clear();
				v_encodings.resize(this->sequences.size());
				v_hashings.resize(this->sequences.size());
				v_spacedKmers.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_naive(this->sequences[seq], this->spaced, v_encodings[seq]);
					getBoth(v_encodings[seq], this->k, v_hashings[seq], v_spacedKmers[seq]);
				}
			}


			inline void GetEncoding_FSH(std::vector<Encoding_V>& v_encodings)
			{
				v_encodings.clear();
				v_encodings.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_speedup_previous(this->sequences[seq], this->spaced, v_encodings[seq]);
			}
			inline void GetEncoding_FSH(std::vector<Encoding_V>& v_encodings, std::vector<Hashing_V>& v_hashings)
			{
				v_encodings.clear();
				v_hashings.clear();
				v_encodings.resize(this->sequences.size());
				v_hashings.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_speedup_previous(this->sequences[seq], this->spaced, v_encodings[seq]);
					getHashes(v_encodings[seq], this->k, v_hashings[seq]);
				}
			}
			inline void GetEncoding_FSH(std::vector<Encoding_V>& v_encodings, std::vector<SpacedKmer_V>& v_spacedKmers)
			{
				v_encodings.clear();
				v_spacedKmers.clear();
				v_encodings.resize(this->sequences.size());
				v_spacedKmers.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_speedup_previous(this->sequences[seq], this->spaced, v_encodings[seq]);
					getSpacedKmers(v_encodings[seq], this->k, v_spacedKmers[seq]);
				}
			}
			inline void GetEncoding_FSH(std::vector<Encoding_V>& v_encodings, std::vector<Hashing_V>& v_hashings, std::vector<SpacedKmer_V>& v_spacedKmers)
			{
				v_encodings.clear();
				v_hashings.clear();
				v_spacedKmers.clear();
				v_encodings.resize(this->sequences.size());
				v_hashings.resize(this->sequences.size());
				v_spacedKmers.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_speedup_previous(this->sequences[seq], this->spaced, v_encodings[seq]);
					getBoth(v_encodings[seq], this->k, v_hashings[seq], v_spacedKmers[seq]);
				}
			}


			inline void GetEncoding_ISSH(std::vector<Encoding_V>& v_encodings)
			{
				v_encodings.clear();
				v_encodings.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_with_ISSH(this->sequences[seq], this->spaced, v_encodings[seq]);
			}
			inline void GetEncoding_ISSH(std::vector<Encoding_V>& v_encodings, std::vector<Hashing_V>& v_hashings)
			{
				v_encodings.clear();
				v_hashings.clear();
				v_encodings.resize(this->sequences.size());
				v_hashings.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH(this->sequences[seq], this->spaced, v_encodings[seq]);
					getHashes(v_encodings[seq], this->k, v_hashings[seq]);
				}
			}
			inline void GetEncoding_ISSH(std::vector<Encoding_V>& v_encodings, std::vector<SpacedKmer_V>& v_spacedKmers)
			{
				v_encodings.clear();
				v_spacedKmers.clear();
				v_encodings.resize(this->sequences.size());
				v_spacedKmers.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH(this->sequences[seq], this->spaced, v_encodings[seq]);
					getSpacedKmers(v_encodings[seq], this->k, v_spacedKmers[seq]);
				}
			}
			inline void GetEncoding_ISSH(std::vector<Encoding_V>& v_encodings, std::vector<Hashing_V>& v_hashings, std::vector<SpacedKmer_V>& v_spacedKmers)
			{
				v_encodings.clear();
				v_hashings.clear();
				v_spacedKmers.clear();
				v_encodings.resize(this->sequences.size());
				v_hashings.resize(this->sequences.size());
				v_spacedKmers.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH(this->sequences[seq], this->spaced, v_encodings[seq]);
					getBoth(v_encodings[seq], this->k, v_hashings[seq], v_spacedKmers[seq]);
				}
			}


			inline void PrintFASTA(const std::vector<SpacedKmer_V>& v_spacedKmers)
			{
				for (const SpacedKmer_V& spacedKmers : v_spacedKmers)
					for (const SpacedKmer& spacedKmer : spacedKmers)
						std::cout << ">\n" << spacedKmer.spacedKmer << "\n";
			}
			inline void PrintFASTA(const std::vector<SpacedKmer_V>& v_spacedKmers, std::string fileName)
			{
				std::ofstream outfile = std::ofstream(fileName + ".fa");

				for (const SpacedKmer_V& spacedKmers : v_spacedKmers)
					for (const SpacedKmer& spacedKmer : spacedKmers)
						outfile << ">\n" << spacedKmer.spacedKmer << "\n";
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
			inline void GetEncoding_naive(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_naive(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
				}
			}
			inline void GetEncoding_naive(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_naive(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_naive(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_naive(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_naive(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_naive(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void GetEncoding_FSH(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_speedup_previous(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
				}
			}
			inline void GetEncoding_FSH(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_speedup_previous(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_FSH(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_speedup_previous(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_FSH(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_speedup_previous(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void GetEncoding_ISSH(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_with_ISSH(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
				}
			}
			inline void GetEncoding_ISSH(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_with_ISSH(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_ISSH(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_with_ISSH(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_ISSH(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					v_encodings_v[seq].resize(this->multi_spaced.size());
					for(size_t ss = 0; ss < this->multi_spaced.size(); ss++)
						GetHashes_with_ISSH(this->sequences[seq], this->multi_spaced[ss], v_encodings_v[seq][ss]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void GetEncoding_FSH_multi(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_speedup_multi_previous_Rotated(this->sequences[seq], this->spaced_qmers, v_encodings_v[seq]);
			}
			inline void GetEncoding_FSH_multi(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_speedup_multi_previous_Rotated(this->sequences[seq], this->spaced_qmers, v_encodings_v[seq]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_FSH_multi(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_speedup_multi_previous_Rotated(this->sequences[seq], this->spaced_qmers, v_encodings_v[seq]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_FSH_multi(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_speedup_multi_previous_Rotated(this->sequences[seq], this->spaced_qmers, v_encodings_v[seq]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void GetEncoding_MISSH_v1(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_with_ISSH_multi_v1(this->sequences[seq], this->multi_spaced, this->VV_shifts, this->v_pos_one, this->max_transient_length, v_encodings_v[seq]);
			}
			inline void GetEncoding_MISSH_v1(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_v1(this->sequences[seq], this->multi_spaced, this->VV_shifts, this->v_pos_one, this->max_transient_length, v_encodings_v[seq]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_v1(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_v1(this->sequences[seq], this->multi_spaced, this->VV_shifts, this->v_pos_one, this->max_transient_length, v_encodings_v[seq]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_v1(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_v1(this->sequences[seq], this->multi_spaced, this->VV_shifts, this->v_pos_one, this->max_transient_length, v_encodings_v[seq]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void GetEncoding_MISSH_col(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_with_ISSH_multi_col(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
			}
			inline void GetEncoding_MISSH_col(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_col(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_col(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_col(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_col(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_col(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void GetEncoding_MISSH_col_parallel(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_with_ISSH_multi_col_parallel(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
			}
			inline void GetEncoding_MISSH_col_parallel(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_col_parallel(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_col_parallel(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_col_parallel(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_col_parallel(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_col_parallel(this->sequences[seq], this->infoCol, v_encodings_v[seq]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void GetEncoding_MISSH_row(std::vector<Encoding_V_V>& v_encodings_v)
			{
				v_encodings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
					GetHashes_with_ISSH_multi_row(this->sequences[seq], this->infoRow, v_encodings_v[seq]);
			}
			inline void GetEncoding_MISSH_row(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_row(this->sequences[seq], this->infoRow, v_encodings_v[seq]);
					getHashes(v_encodings_v[seq], this->k, v_hashings_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_row(std::vector<Encoding_V_V>& v_encodings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_row(this->sequences[seq], this->infoRow, v_encodings_v[seq]);
					getSpacedKmers(v_encodings_v[seq], this->k, v_spacedKmers_v[seq]);
				}
			}
			inline void GetEncoding_MISSH_row(std::vector<Encoding_V_V>& v_encodings_v, std::vector<Hashing_V_V>& v_hashings_v, std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				v_encodings_v.clear();
				v_hashings_v.clear();
				v_spacedKmers_v.clear();
				v_encodings_v.resize(this->sequences.size());
				v_hashings_v.resize(this->sequences.size());
				v_spacedKmers_v.resize(this->sequences.size());
				
				for (size_t seq = 0; seq < this->sequences.size(); seq++)
				{
					GetHashes_with_ISSH_multi_row(this->sequences[seq], this->infoRow, v_encodings_v[seq]);
					getBoth(v_encodings_v[seq], this->k, v_hashings_v[seq], v_spacedKmers_v[seq]);
				}
			}


			inline void PrintFASTA(const std::vector<SpacedKmer_V_V>& v_spacedKmers_v)
			{
				for (size_t ss = 0; ss < multi_spaced.size(); ss++)
					for (size_t seq = 0; seq < sequences.size(); seq++)
						for (size_t j = 0; j < sequences[seq].size() - multi_spaced[ss].GetQ() + 1; j++)
							std::cout << ">\n" << v_spacedKmers_v[seq][ss][j].spacedKmer << "\n";
			}
			inline void PrintFASTA(const std::vector<SpacedKmer_V_V>& v_spacedKmers_v, std::string fileName)
			{
				for (size_t ss = 0; ss < multi_spaced.size(); ss++)
				{
					std::ofstream outfile = std::ofstream(fileName + "_" + std::to_string(ss) + "_.fa");

					for (size_t seq = 0; seq < sequences.size(); seq++)
						for (size_t j = 0; j < sequences[seq].size() - multi_spaced[ss].GetQ() + 1; j++)
							outfile << ">\n" << v_spacedKmers_v[seq][ss][j].spacedKmer << "\n";
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
