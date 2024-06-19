/*
 * SpacedQmer_Multi.h
 *
 *	Modified by: Leonardo
 */


#ifndef INCLUDE_SPACED_SPACEDQMER_MULTI_H_
	#define INCLUDE_SPACED_SPACEDQMER_MULTI_H_

	#include "Spaced/SpacedQmer.h"


	//For multi unit
	struct MapUnit
	{
		std::vector<size_t> n_one;
		std::vector<V_Pos_Ones> v_v_pos;
	};

	//for previous
	struct PreviusShift_Ext : public PreviousShift 
	{
		size_t prev_qmer = std::numeric_limits<size_t>::max();
		const SpacedQmer* current_sp_ptr;
		const SpacedQmer* prev_sp_ptr;
		
		bool isCorrectSpacedPrevious() const
		{
			return this->prev_qmer != std::numeric_limits<size_t>::max();
		}
	};
	typedef std::vector<PreviusShift_Ext> PreviusShift_Ext_V;


	inline static void GetTableShiftFirstOnSecond(	const std::vector<SpacedQmer>& v_spaced, 
													std::vector<std::vector<PreviusShift_Ext_V>>& table_shift_change)
	{
		table_shift_change.clear();
		table_shift_change.resize(v_spaced.size());
		for (size_t i = 0; i < table_shift_change.size(); ++i)
			table_shift_change[i].resize(v_spaced.size());

		// calcola posizioni da salvare con i vari shift e poi salva
		for (size_t ss = 0; ss < v_spaced.size(); ++ss)
		{
			for (size_t s = 0; s < v_spaced.size(); ++s) // devo far shiftare su tutti i seed
			{
				// reset shift di ss rispetto a s
				table_shift_change[ss][s].resize(v_spaced[s].toString().size()); // deve esser lungo almeno quanto la lunghezza di s
				for (size_t i = 0; i < table_shift_change[ss][s].size(); ++i)
					table_shift_change[ss][s][i].prev_qmer = s; // qmer di riferimento per calcolare l'hash dello spaced ss

				// calcola posizioni da salvare con i vari shift poi salva
				size_t init = 0;
				for (size_t i = 0; i < v_spaced[s].toString().size(); ++i) // per tutti gli shift possibili (calcola anche lo 0 anche se non necessario per il momento)
				{
					// Cerco indice di partenza su vettore di posizione degli uno del spaced che resta fermo, quindi s
					bool find = false;
					for (size_t j = init; j < v_spaced[s].GetPosOne().size(); ++j)
					{
						if (v_spaced[s].GetPosOne()[j] >= i)
						{
							init = j;
							find = true;
							break;
						}
					}
					
					if (!find)
						init = v_spaced[s].GetPosOne().size(); // Serve per saltare prossimo ciclo senza altri controlli

					// per tutte le posizioni di uno da analizzare
					size_t j;
					for (j = init; j < v_spaced[s].GetPosOne().size() && (j - init) < v_spaced[ss].GetPosOne().size(); ++j)
					{
						// confronta indici incolonnati tenendo conto dello shift (del vettore con indice s),
						// verifica se diversi, se si c'è operazione di rimozione e inserimento da fare in quel punto
						// j-init parte da 0, quindi 0-1-2-3... etc. quindi sto confrontando il primo pezzo
						if (v_spaced[ss].GetPosOne()[j - init] != v_spaced[s].GetPosOne()[j] - i)
						{
							table_shift_change[ss][s][i].one_to_remove.push_back(j - init);
							table_shift_change[ss][s][i].one_to_change.push_back(j - init);
						}
					}

					////////////////////////////////////////////////////////////////////
					// rimuovi dopo l'ultimo 1 analizzato del vettore riferimento (s)
					// aggiungi dopo ultimo 1 su spaced da calcolare (ss)

					//	TODO:versione numero 1 seed diversa
					//	for(size_t h = j; h < v_spaced[s].GetPosOne().size(); ++h)//rimanenti 1 di s da eliminare
					//		shift_firstid_on_secondid[ss][s][i].one_to_remove.push_back(h);
					//	for(size_t h = j-init; h < v_spaced[ss].GetPosOne().size(); ++h)//rimanenti 1 di ss da inserire
					//		shift_firstid_on_secondid[ss][s][i].one_to_change.push_back(h);

					////////////////////////////////////////////////////////////////////

					// salva il numero di shift di 1 da effettuare
					table_shift_change[ss][s][i].one_exit = init; // shift di uno da effettuare su s
				}
			}
		}
	}


	inline static void GetShiftMinChange(	const std::vector<SpacedQmer>& v_spaced,
											std::vector<PreviusShift_Ext_V>& v_shift_min)
	{
		// Inizializza matrice per salvare tutti i possibili shift tutti contro tutti
		std::vector<std::vector<PreviusShift_Ext_V>> shift_firstid_on_secondid;
		GetTableShiftFirstOnSecond(v_spaced, shift_firstid_on_secondid);

		size_t max_size = 0;
		for (size_t ss = 0; ss < v_spaced.size(); ++ss)
			for (size_t s = 0; s < v_spaced.size(); ++s)
				if (max_size < shift_firstid_on_secondid[ss][s].size())
					max_size = shift_firstid_on_secondid[ss][s].size();

		// Ora ho tutti gli shift possibili tra tutti. Devo scegliere quali effettuare ad ogni passo,
		// tenendo conto delle rimozioni, che posson esser diverse dagli inserimenti da fare.
		v_shift_min.clear();
		v_shift_min.resize(v_spaced.size());
		for (size_t ss = 0; ss < v_spaced.size(); ++ss)
		{
			// v_shift_min[ss].resize(1); // almeno un elemento, i==0 deve esistere vuoto, shift parte da 1
			for (size_t i = 0; i < max_size; ++i)
			{
				size_t limit_evaluation_previous = (i == 0) ? ss : v_spaced.size();
				for (size_t s = 0; s < limit_evaluation_previous; ++s)
				{
					if (i < shift_firstid_on_secondid[ss][s].size())
					{
						PreviusShift_Ext current_shift_min = shift_firstid_on_secondid[ss][s][i];
						size_t current_size = current_shift_min.GetSize();
						current_shift_min.shift_min = i; // Guess

						if (v_shift_min[ss].size() <= i) // o non c'è niente salvato
						{
							if (i > 0) // lo esegue solo 1 volta al primo seed che valuta, poi so già quale è il migliore poiché viene salvato
							{
								const PreviusShift_Ext& saved_prev_shift_min = v_shift_min[ss][i - 1];
								size_t saved_size = saved_prev_shift_min.GetSize();

								bool size_saved_better = saved_size < current_size;
								if (size_saved_better && saved_prev_shift_min.isCorrectSpacedPrevious())
								{
									current_shift_min = saved_prev_shift_min;
									current_size = current_shift_min.GetSize();
								}
							}
							if (!(i == 0 && current_size > v_spaced[ss].GetWeight()))
								v_shift_min[ss].push_back(current_shift_min);
						}
						else
						{
							const PreviusShift_Ext& saved_shift_min = v_shift_min[ss][i];
							size_t saved_size = saved_shift_min.GetSize();

							bool size_saved_better = saved_size < current_size; // < preferisci il nuovo corrente, con <= preferisco il salvato
							if (!size_saved_better || !saved_shift_min.isCorrectSpacedPrevious())
								v_shift_min[ss][i] = current_shift_min;
						}
					}
				}
				if (i == 0 && v_shift_min[ss].empty())
					v_shift_min[ss].push_back(PreviusShift_Ext());
			}
		}

		// TODO: da salvare anche gli shift migliori sul finale se i seed han dimensioni diverse,
		// cioè per ogni seed e per ogni shift, salvare il migliore tra tutti i seed testati per quel dato shift.
		// questo, insieme al transitorio iniziale, permette di capire da chi e meglio computare cosa.

		// Settare ptr
		for (size_t ss = 0; ss < v_shift_min.size(); ++ss)
		{
			for (size_t i = 0; i < v_shift_min[ss].size(); ++i)
			{
				v_shift_min[ss][i].current_sp_ptr = &v_spaced[ss];
				if (v_shift_min[ss][i].isCorrectSpacedPrevious())
					v_shift_min[ss][i].prev_sp_ptr = &v_spaced[v_shift_min[ss][i].prev_qmer];
			}
		}
	}

	inline static void GetShiftMinChangeRotated(	const std::vector<SpacedQmer>& v_spaced,
													std::vector<PreviusShift_Ext_V>& v_shift_min)
	{
		//Inizializza matrice per salvare tutti i possibili shift tutti contro tutti
		std::vector<PreviusShift_Ext_V> shift_firstid_on_secondid;
		GetShiftMinChange(v_spaced, shift_firstid_on_secondid);
		size_t size_max = max_element(	shift_firstid_on_secondid.begin(), 
										shift_firstid_on_secondid.end(),
										[](PreviusShift_Ext_V& a, PreviusShift_Ext_V& b)
										{
											return a.size() < b.size();
										})->size();
		v_shift_min.clear();
		v_shift_min.resize(size_max, PreviusShift_Ext_V(v_spaced.size()));
		for (size_t ss = 0; ss < shift_firstid_on_secondid.size(); ++ss)
		{
			for (size_t i = 0; i < size_max; ++i)
			{
				if (i < shift_firstid_on_secondid[ss].size())//elemento c'è, quindi copia
					v_shift_min[i][ss] = shift_firstid_on_secondid[ss][i];
				else
					v_shift_min[i][ss] = shift_firstid_on_secondid[ss].back();
			}
		}
	}

	class SpacedQmer_Multi {
		public:
			SpacedQmer_Multi() {}
			SpacedQmer_Multi(const std::vector<SpacedQmer>& v_spaced)
			{
				this->init(v_spaced);
			}

			inline SpacedQmer& operator[](size_t i)
			{
				return this->v_spaced[i];
			}

			inline const SpacedQmer& operator[](size_t i) const
			{
				return this->v_spaced[i];
			}

			inline size_t size() const
			{
				return this->v_spaced.size();
			}

			inline const MapUnit& getMapUnit() const
			{
				return this->map_unit;
			}

			inline const std::vector<PreviusShift_Ext_V>& getShiftMin() const
			{
				return this->v_shift_min;
			}

			inline const std::vector<PreviusShift_Ext_V>& getShiftMinRotated() const
			{
				return this->v_shift_min_rotated;
			}

			inline void init(const std::vector<SpacedQmer>& v_spaced)
			{
				this->v_spaced = v_spaced;
				GetShiftMinChange(this->v_spaced, this->v_shift_min);
				GetShiftMinChangeRotated(this->v_spaced, this->v_shift_min_rotated);
			}

		private:
			std::vector<SpacedQmer> v_spaced;
			MapUnit map_unit;
			std::vector<PreviusShift_Ext_V> v_shift_min;
			std::vector<PreviusShift_Ext_V> v_shift_min_rotated;
	};

#endif /* INCLUDE_SPACED_SPACEDQMER_MULTI_H_ */
