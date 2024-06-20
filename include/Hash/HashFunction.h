/*
 * HashFunction.h
 *
 *	Modified by: Leonardo
 */

#ifndef INCLUDE_HASH_HASHFUNCTION_H_
	#define INCLUDE_HASH_HASHFUNCTION_H_

	#include "Hash/HashType.h"




	// Hash per spaced qmer con tutti 1
	inline static void GetHash(const std::string& s_Str, const size_t& startQmer, const size_t& length, Encoding& encoding)
	{
		encoding = 0;
		for(size_t j = 0; j < length; j++)
		{
			const unsigned char& tmp_char = s_Str[startQmer + j];
			encoding |= fEncodingTable[(tmp_char << 5) + j];
		}
	}

	// Hash per spaced qmer con *
	inline static void GetHash(const std::string& s_Str, const size_t& startQmer, const Position& pos_one, Encoding& encoding)
	{
		encoding = 0;
		for(size_t j = 0; j < pos_one.size(); j++)
		{
			const unsigned char& tmp_char = s_Str[startQmer + pos_one[j]];
			encoding |= fEncodingTable[(tmp_char << 5) + j];
		}
	}

	inline static void GetHash(const std::string& s_Str, const size_t& startQmer, const SpacedQmer& spaced_qmer, Encoding& encoding)
	{
		GetHash(s_Str, startQmer, spaced_qmer.GetPosOne(), encoding);
	}
	




	// ----------------------------------- //
	//   G E T   H A S H E S   N A I V E   //
	// ----------------------------------- //
	inline static void GetHashes_naive(const std::string& s_Str, const SpacedQmer& spaced_qmer, Encoding_V& encodings)
	{
		encodings.clear();

		if (s_Str.size() < spaced_qmer.GetQ()) // Almeno 1 hash
			return;

		const size_t n_encodings = s_Str.size() - spaced_qmer.GetQ() + 1;
		encodings.resize(n_encodings);

		const Position& pos_one = spaced_qmer.GetPosOne();
		for (size_t pos = 0; pos < n_encodings; pos++)
			GetHash(s_Str, pos, pos_one, encodings[pos]);
	}




	// --------------------------------------- //
	//   G E T   H A S H E S   S P E E D U P   //
	// --------------------------------------- //
	inline static void _compute_encoding_for_speedup_previous(	const std::string& s_Str, const PreviousShift& curr_sp_shift,
																const Position& pos_one_current, const Position& pos_one_prev,
																const Encoding& prev_encoding, Encoding& curr_encoding, size_t& idx_curr_encoding)
	{
		curr_encoding = prev_encoding >> (2 * curr_sp_shift.one_exit);
		

		// Rimuovo gli 1 da eliminare
		uint64_t mask = (2 * pos_one_current.size() != 64) ? ((uint64_t)0b01 << (2 * pos_one_current.size())) - 1 : (uint64_t)(-1);
		for (size_t pos : curr_sp_shift.one_to_remove)
			mask ^= (uint64_t)0b11 << (2 * pos);
		curr_encoding &= mask;
		


		// Aggiorna posizioni da cambiare su hash
		for (const size_t& i_to_change : curr_sp_shift.one_to_change)
		{
			const unsigned char& tmp_char = s_Str[idx_curr_encoding + pos_one_current[i_to_change]];
			curr_encoding |= fEncodingTable[(tmp_char << 5) + i_to_change];
		}


		// Aggiunge i valori all'estremità destra del seed, quelli che non potrei calcolare con sovrapposizioni
		if (pos_one_current.size() == pos_one_prev.size())
		{
			for (size_t j = pos_one_current.size() - curr_sp_shift.one_exit; j < pos_one_current.size(); j++)
			{
				const unsigned char& tmp_char = s_Str[idx_curr_encoding + pos_one_current[j]];
				curr_encoding |= fEncodingTable[(tmp_char << 5) + j];
			}
		}
	}


	inline static void GetHashes_speedup_previous(const std::string& s_Str, const SpacedQmer& spaced_qmer, Encoding_V& encodings)
	{
		auto _get_encoding = [&](size_t curr_idx_encoding, const PreviousShift& curr_shift)
		{
			Encoding& curr_encoding = encodings[curr_idx_encoding];

			if (spaced_qmer.GetWeight() < curr_shift.GetSize())
				GetHash(s_Str, curr_idx_encoding, spaced_qmer.GetPosOne(), curr_encoding);
			else
			{
				const Encoding& prev_encoding = encodings[curr_idx_encoding - curr_shift.shift_min];
				_compute_encoding_for_speedup_previous(	s_Str, curr_shift,
														spaced_qmer.GetPosOne(), spaced_qmer.GetPosOne(),
														prev_encoding, curr_encoding, curr_idx_encoding);
			}
		};



		encodings.clear();

		if (s_Str.size() < spaced_qmer.GetQ()) // Almeno 1 hash
			return;

		const size_t n_encodings = s_Str.size() - spaced_qmer.GetQ() + 1;
		encodings.resize(n_encodings);

		const V_PreviusShift& shift = spaced_qmer.GetShiftMinChange();

		// Devo computare il primo hash a parte
		GetHash(s_Str, 0, spaced_qmer.GetPosOne(), encodings[0]);

		size_t i = 1;
		// Per tutte le posizioni che contemplano gli shift nel primo pezzo di sequenza
		for (; i < std::min(shift.size(), n_encodings); i++)
			_get_encoding(i, shift[i]);

		const PreviousShift& last_shift = shift.back();
		for (; i < n_encodings; i++)
			_get_encoding(i, last_shift);
	}




	// --------------------------------- //
	//   G E T   H A S H E S   I S S H   //
	// --------------------------------- //
	inline static void _compute_encoding_with_ISSH(	const std::string& s_Str, const V_PreviusShift& shifts,
													const Position& pos_one,
													Encoding_V& encodings, size_t& curr_idx_encoding)
	{
		Encoding& curr_encoding = encodings[curr_idx_encoding];
		curr_encoding = 0;

		// Vado a prendere le posizioni che non posso recuperare tramite gli hash precedenti
		// e che andrò a calcolare tramite la funzione di hash alla fine.
		// Questo vettore sicuramente conterrà l'ultima posizione.
		// È stato inserito dentro a shift[0].one_to_change perché era un vettore di posizioni
		// presente in FSH e che non veniva utilizzato in ISSH.
		const Position& pos_not_covered_yet = shifts[0].one_to_change;
		
		
		// Faccio un ciclo per ciascuno degli hash del gruppo da cui devo recuperare le posizioni.
		for (const PreviousShift& curr_sp_shift : shifts)
		{
			// In partial salvo l'hash precedente da cui voglio recuperare posizioni.
			// Faccio lo shift degli hash in modo da avere la giusta sovrapposizione.
			// Facendo l'and con la maschera cancello i bit degli hash che sono sbagliati e mantengo quelli giusti.
			// Infine aggiungo all'hash corrente
			Encoding partial = encodings[curr_idx_encoding - curr_sp_shift.shift_min];
			partial >>= (2 * curr_sp_shift.one_exit);
			partial &= curr_sp_shift.mask;
			curr_encoding |= partial;
		}

		// A questo punto si devono inserire, calcolati uno ad uno con la funzione di codifica,
		// tutti i valori contenuti nelle posizioni di pos_not_covered_yet.
		// Dentro questo vettore ho anche la posizione dell'ultimo valore dell'hash, che va calcolato in ogni caso.
		for (const size_t& i_to_change : pos_not_covered_yet)
		{
			const unsigned char& tmp_char = s_Str[curr_idx_encoding + pos_one[i_to_change]];
			curr_encoding |= fEncodingTable[(tmp_char << 5) + i_to_change];
		}
	}


	inline static void GetHashes_with_ISSH(const std::string& s_Str, const SpacedQmer& spaced_qmer, Encoding_V& encodings)
	{
		encodings.clear();

		if (s_Str.size() < spaced_qmer.GetQ()) // Almeno 1 hash
			return;

		const size_t n_encodings = s_Str.size() - spaced_qmer.GetQ() + 1;
		encodings.resize(n_encodings);

		// Carico il vettore dei gruppi di shift precedenti che mi serve per
		// gestire sia il caso a regime sia tutto il transitorio.
		const V_V_PreviusShift& V_shifts = spaced_qmer.GetMultipleShifts();

		const Position& pos_one = spaced_qmer.GetPosOne();

		// Per lo 0-esimo hash ovviamente non posso recuperare niente:
		// lo calcolo posizione per posizione.
		GetHash(s_Str, 0, pos_one, encodings[0]);

		// Il transitorio mi basta farlo dal primo alla lunghezza del vettore di vettori
		// di shift che voglio andare a utilizzare che gestisce il transitorio dall'1
		// fino alla fine. Ad ogni passo del trasitorio chiamo la stessa funzione
		// passandogli il gruppo di hash precedenti corrispondente

		size_t i = 1;
		// Transitorio
		for (; i < std::min(V_shifts.size(), n_encodings); i++)
		{
			// Se il vettore di shift è vuoto è sicuro che non posso recuperare nessuna posizione:
			// per evitare controlli e passaggi inutili si utilizza la funzione posizione per posizione.
			if (V_shifts[i].size() == 0)
				GetHash(s_Str, i, pos_one, encodings[i]);
			else
				_compute_encoding_with_ISSH(s_Str, V_shifts[i], pos_one, encodings, i);
		}

		// Regime
		// Nel caso a regime si passa alla funzione compute sempre lo stesso gruppo di hash
		const V_PreviusShift& first_shifts = V_shifts[0];
		for (; i < n_encodings; i++)
		{
			// In questo caso si suppone che a regime siano recuperabili tutte le posizioni,
			// quindi non ha senso vedere se il vettore di posizioni precedenti è vuoto.
			_compute_encoding_with_ISSH(s_Str, first_shifts, pos_one, encodings, i);
		}
	}




	// ----------------------------------------- //
	//   G E T   H A S H E S   O N E   P A S S   //
	// ----------------------------------------- //
	inline static void GetHashes_one_pass(const std::string& s_Str, const SpacedQmer& spaced_qmer, Encoding_V& encodings)
	{
		encodings.clear();

		const size_t SEED_Q = spaced_qmer.GetQ();
		if (s_Str.size() < SEED_Q) // Almeno 1 hash
			return;

		const size_t n_encodings = s_Str.size() - SEED_Q + 1;
		encodings.resize(n_encodings);

		
		
		const Position& pos_one = spaced_qmer.GetPosOne();
		const size_t POS_ONE_SIZE = pos_one.size();
		
		// Calcolo gli hash da 0 a spaced_qmer.GetQ() e da nHashes - spaced_qmer.GetQ() a nHashes
		for (size_t pos = 0; pos < SEED_Q; pos++)
		{
			const size_t n_encodings_PLUS_POS = n_encodings + pos;

			const unsigned char& char2add_testa = s_Str[pos];
			const size_t offset_testa = char2add_testa << 5;

			const unsigned char& char2add_coda = s_Str[n_encodings_PLUS_POS];
			const size_t offset_coda = char2add_coda << 5;

			for (size_t j = 0; j < POS_ONE_SIZE; j++)
			{
				const size_t POS_ONE_J = pos_one[j];
				if (POS_ONE_J > pos)
				{
					// Sto calcolando gli hash di coda, da nHashes - spaced_qmer.GetQ() a nHashes
					encodings[n_encodings_PLUS_POS - POS_ONE_J] |= fEncodingTable[offset_coda + j];
				}
				else 
				{
					// Sto calcolando gli hash di testa, da 0 a spaced_qmer.GetQ()
					encodings[pos - POS_ONE_J] |= fEncodingTable[offset_testa + j];
				}
			}
		}

		// Calcolo gli hash da spaced_qmer.GetQ() a nHashes - spaced_qmer.GetQ()
		for (size_t pos = SEED_Q; pos < n_encodings; pos++)
		{
			const unsigned char& char2add = s_Str[pos];
			const size_t offset = char2add << 5;

			for (size_t j = 0; j < POS_ONE_SIZE; j++)
			{
				encodings[pos - pos_one[j]] |= fEncodingTable[offset + j];
			}
		}
	}




#endif // INCLUDE_HASH_HASHFUNCTION_H_