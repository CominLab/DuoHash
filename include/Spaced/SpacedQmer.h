/*
 * SpacedQmer.h
 *
 *   Created on: 20/lug/2016
 *       Author: samuele
 *	Modified by: enrico
 *	Modified by: Leonardo
 */

#ifndef INCLUDE_SPACED_SPACEDQMER_H_
	#define INCLUDE_SPACED_SPACEDQMER_H_

	#include <algorithm>
	#include <bitset>
	#include <cmath>
	#include <cstdint>
	#include <iostream>
	#include <limits>
	#include <string>
	#include <vector>


	//For unit
	struct Pos_Ones
	{
		size_t index_one = std::numeric_limits<size_t>::max();
		size_t n_one = 0;
		size_t pos_start = 0;
		size_t n_one_before = 0; //numero di 1 che precedono il primo 1 dello unit
	};
	typedef std::vector<Pos_Ones> V_Pos_Ones;

	//for previous
	//position is a vector of index referring to the positions inside the seed
	typedef std::vector<size_t> Position;
	typedef std::vector<Position> V_Position;

	// PreviousShift è la struttura che mi contiene tutte le informazioni
	// riguardanti uno shift precedente a quello che sto considerando dello spaced seed
	// mi serve per calcolare più velocemente lo speedup previous
	struct PreviousShift
	{
		// vettori utilizzati principalmente da FSH
		Position one_to_change;
		Position one_to_remove;

		Position one_to_keep;
		
		// one_exit è l'indice del primo uno sovrapponibile del vettore di uno dello
		// spaced seed di quanto devo traslare l'hash.
		size_t one_exit = 0;
		
		// shift_min rappresenta lo shift dell'attuale previousShift rispetto a quello "principale"
		size_t shift_min = 0;
		
		// maschera necessaria al calcolo tramite ISSH
		uint64_t mask = 0;
		
		inline size_t GetSize() const
		{
			return this->one_to_change.size() + this->one_to_remove.size() + this->one_exit + (this->one_to_remove.empty() ? 0 : 2);
		}
	};

	// associato ad uno spacedqmer ho un vettore delle possibili traslazioni dello stesso seed
	typedef std::vector<PreviousShift> V_PreviusShift;
	typedef std::vector<V_PreviusShift> V_V_PreviusShift;




	// funzioni di supporto che stampano qualche risultato
	void printp(Position p);
	void print_shift(PreviousShift s);

	// funzioni di supporto che permettono di gestire con più chiarezza il vettore
	// delle posizioni che non sono ancora state recuperate
	void deleteElement(Position& pos, size_t index);
	bool isContained(Position pointer, Position pos, size_t index);




	class SpacedQmer {
		public:
			SpacedQmer();
			SpacedQmer(std::string spaced_qmer, size_t numprev);

			inline size_t GetWeight() const
			{
				return this->pos_one.size();
			}

			inline size_t GetQ() const
			{
				return this->spaced_q.length();
			}

			inline bool isOne(size_t index) const
			{
				return this->spaced_q[index] == '1';
			}

			inline const Position& GetPosOne() const
			{
				return this->pos_one;
			}

			inline const V_PreviusShift& GetShiftMinChange() const 
			{
				return this->shift_min_change;
			}

			inline const V_V_PreviusShift& GetMultipleShifts() const
			{
				return this->multiple_shifts;
			}
			
			inline const V_V_PreviusShift* GetMultipleShiftsPointer() const
			{
				return &this->multiple_shifts;
			}

			inline const std::string& toString() const
			{
				return this->spaced_q;
			}

			void reset(std::string spaced_qmer, size_t numprev);

		private:
			// the actual string of ones and zeros, the original spaced seed
			std::string spaced_q;

			//pos one is a vector of index corresponding to ones in the seed
			Position pos_one;

			// vettore di posizioni nel vettore pos_one dove inizialmente devo ancora recupereare
			// alcuna corrispondenza precedente.
			Position pos_pos_one;

			//it contains all the data of all the possible overlapping shifted seeds
			V_PreviusShift shift_min_change;

			// è un vettore di gruppi di hash precedenti che voglio riutilizzare, il vettore che
			// mi serve per calcolare gli hash gestisce sia regime sia transitorio
			V_V_PreviusShift multiple_shifts;

			// numero di hash precedenti che utilizzo a regime per recuperare le posizioni
			// conterrà il valore passato da terminale
			size_t num_prev = 0;

			void SaveIndexOne();
			void GetShiftMax(V_PreviusShift& shift_max);
			void SetMultipleShifts(size_t index);
			void SetAllMultipleShift();
			void SetBitMasks();
	};

	// funzioni di supporto che stampano qualche risultato
	void printp(Position p);
	void print_shift(PreviousShift s);

	// funzioni di supporto che permettono di gestire con più chiarezza il vettore
	// delle posizioni che non sono ancora state recuperate
	void deleteElement(Position& pos, size_t index);
	bool isContained(Position pointer, Position pos, size_t index);

#endif /* INCLUDE_SPACED_SPACEDQMER_H_ */
