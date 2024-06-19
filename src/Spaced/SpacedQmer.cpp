/*
 * SpacedQmer.cpp
 *
 *	Modified by: Leonardo
 */

#include "Spaced/SpacedQmer.h"


SpacedQmer::SpacedQmer()
{
	this->reset("", 0);
}

SpacedQmer::SpacedQmer(std::string spaced_qmer, size_t numprev)
{
	this->reset(spaced_qmer, numprev);
}

void SpacedQmer::reset(std::string spaced_qmer, size_t numprev)
{
	this->num_prev = numprev;
	this->spaced_q = spaced_qmer;
	this->SaveIndexOne();
	
	// chiamo la funzione che fa i calcoli per l'utilizzo di un singolo hash precedente
	this->GetShiftMax(this->shift_min_change);
	
	// chiamo la funzione che fa i calcoli per trovare i gruppi di shift prencedenti
	// per calcolare l'hash
	this->SetAllMultipleShift();
}

void SpacedQmer::SaveIndexOne()
{
	this->pos_one.clear();
	
	size_t k = 0;
	for (size_t i = 0; i < this->spaced_q.size(); ++i)
		if (this->isOne(i))
		{
			this->pos_one.push_back(i);
			this->pos_pos_one.push_back(k);
			k++;
		}
}

void SpacedQmer::GetShiftMax(V_PreviusShift& shift_max)
{
	shift_max.clear();
	shift_max.resize(this->spaced_q.size());
	
	size_t init = 0;

	bool find;
	for (size_t i = 1; i < this->spaced_q.size(); ++i) // per tutti gli shift possibili a parte il primo
	{
		find = false;
		for (size_t j = init; j < this->pos_one.size(); ++j)
		{
			if (this->pos_one[j] >= i)
			{
				init = j;
				find = true;
				break;
			}
		}
		
		if (!find)
			init = this->pos_one.size(); // Serve per saltare prossimo ciclo senza altri controlli
		
		for (size_t j = init; j < this->pos_one.size(); ++j) // per tutte le posizioni del secondo vettore traslato
		{
			if (this->pos_one[j - init] != this->pos_one[j] - i)
			{
				shift_max[i].one_to_remove.push_back(j - init);
				shift_max[i].one_to_change.push_back(j - init);
			}
			else
				shift_max[i].one_to_keep.push_back(j - init);
		}

		//		il prossimo passaggio non è necessario in quanto i bit si spostano tutti di 2*init quindi son già annullati,
		//		basta ricordarsi di settarli nella funzione di hash
		//		for(size_t j = this->pos_one.size()-init; j < this->pos_one.size(); ++j)//rimanenti da inserire tutti
		//			this->shift_min_change[i].one_to_change.push_back(j);
		
		shift_max[i].one_exit = init;
		shift_max[i].shift_min = i; // Guess

		const PreviousShift& prev_shift_min = shift_max[shift_max[i - 1].shift_min];
		size_t size_previus = prev_shift_min.GetSize();
		size_t size_current = shift_max[i].GetSize();
		if (i > 1 && size_previus < size_current)
			shift_max[i] = shift_max[i - 1];
	}
}

void SpacedQmer::SetMultipleShifts(size_t index)
{
	// preparo il vettore di multiple_shifts[index] che contiene le informazioni
	// sugli hash precedenti da riutilizzare
	this->multiple_shifts[index].clear();

	// Conterrà la posizioni relative al gruppo corrente che non è stato possibile
	// recuperare
	Position pos_not_covered_yet = this->pos_pos_one;
	// Se index è zero si sta calcolando il gruppo migliore a regime, altrimenti
	// si sta calcolando il transitorio: ci sarà un numero inferiore di shift
	// precedenti da poter utilizzare.
	// Devo gestire due possibili limitazioni: o il transitorio o quella limitata
	// dal parametro -num.
	size_t num_max_previous;
	if (this->num_prev == 0)
		num_max_previous = (index > this->spaced_q.size() || index == 0) ? this->spaced_q.size() : index;
	else
		num_max_previous = (index > this->spaced_q.size() || index == 0) ? this->num_prev : std::min(index, this->num_prev);

	// si continua a cercare di recuperare posizioni da hash precedenti finché non
	// si raggiunge il numero massimo di previous impostato o finisco le posizioni
	// recuperabili sicuramente una posizione (l'ultima) va calcolata tramite la
	// funzione di codifica.

	for (size_t k = 0; k < num_max_previous && pos_not_covered_yet.size() > 1; k++)
	{
		// è stato deciso di usare il vettore di shift sempre tutto riempito e di
		// usare push_back per inserire nuovi valori in modo da avere sempre la
		// lunghezza minima. Trovato il migliore curr_best viene memorizzato.
		PreviousShift curr_best;
		
		// max_num_shifts è diverso da num_max_previous perché nel caso di -num x
		// non è detto che si voglia utilizzare gli shift più vicini a quello che si
		// deve calcolare, ma deve poter cercare gli x hash precedenti migliori su
		// tutti i precedenti. L'unica limitazione per max_num_shifts è dovuta al
		// transitorio
		size_t max_num_shifts = (index != 0) ? index : this->spaced_q.size();

		// Si effettua un ciclo su tutti i possibili shift del seed che
		// rappresentano gli hash disponibili per identificare quello che permette
		// di recuperare più posizioni
		for(size_t i = 1; i <= max_num_shifts; ++i)// per tutti gli shift possibili a parte il primo
		{
			// ciclo su tutto il vettore di posizioni di uno, con cardinalità pari al peso
			// dello spaced seed che effettua il calcolo delle sovrapposizioni per ogni
			// possibile punto di attacco. Questo calcolo potrebbe essere effettuato
			// in modo più intelligente ricordando le posizioni corrispondenti allo
			// stesso seed e allo stesso punto di attacco e facendo in modo di non
			// fare calcoli in più.
			// Anche l'hash precedente che si considera ha la possibilità a sua
			// volta di essere traslato rispetto all'hash che sto calcolando rispetto
			// al punto di attacco.
			
			// init è l'offset della traslazione dell'hash da cui voglio provare a
			// recuperare le posizioni.

			// init parte da uno: per lo zero non avrei mai sovrapposizioni perché
			// significherebbe cercare di recuperare qualche posizione sovrapponendo
			// a partire dal primo valore di entrambi l'hash da calcolare ad uno precedente
			// faccio finire il ciclo quando ho finito tutti gli uno a sinistra del
			// primo uno del seed che sto calcolando. Il numero di uno devo calcolarlo
			// con un ciclo for

			size_t num_one_before_shift = 0;
			for (; num_one_before_shift < this->pos_one.size() && this->pos_one[num_one_before_shift] < i; num_one_before_shift++);

			for (size_t init = 1; init <= num_one_before_shift; ++init)
			{
				// ho bisogno di una struttura di supporto temporanea per confrontarla con
				// quella già salvata in curr_best
				PreviousShift temp;
			
				for (size_t j = init; j < this->pos_one.size(); ++j) //per tutte le posizioni del secondo vettore traslato
				{
					// confronta indici incolonnati tenendo conto dello shift, verifica se
					//  diversi. Se lo sono c'è un'operazione da fare in quel punto.
					if ((j - init) < this->pos_one.size())
					{
						if (this->pos_one[j] < i || this->pos_one[j - init] != this->pos_one[j] - i || !isContained(pos_not_covered_yet, this->pos_one, this->pos_one[j - init]))
						{
							// questi due non servono nel caso ISSH in realtà. one_to_remove
							// viene utilizzato solamente per generare le maschere, mentre
							// one_to_change viene utilizzato solo nel primo hash da
							// utilizzare di ciascun gruppo per identificare le posizioni
							// da ricalcolare rispetto a tutto il gruppo.
							temp.one_to_remove.push_back(j - init);
							temp.one_to_change.push_back(j - init);
						}
						else
							temp.one_to_keep.push_back(j - init);
					}
				}

				// salvo del temp corrente a quanti hash fa mi sto riferendo e a che
				// traslazione di questo
				temp.one_exit = init;
				temp.shift_min = i;
				
				// controllo se la posizione ultima che ho trovato è migliore di
				// quella trovata precedentemente
				if (temp.one_to_keep.size() > curr_best.one_to_keep.size())
					curr_best = temp;
			}
		}

		for (size_t j = 0; j < curr_best.one_to_keep.size(); j++)
			deleteElement(pos_not_covered_yet, curr_best.one_to_keep[j]);
		
		// inserisco in questo gruppo di shift lo shift migliore appena calcolato solo se
		// permette di recuperare almeno una posizione
		if (curr_best.one_to_keep.size() != 0)
			this->multiple_shifts[index].push_back(curr_best);
	}
	
	// Necessario perché passare alle funzioni che calcolano l'hash un ulteriore vettore
	// di posizioni da cambiare aumenta di molto il tempo totale. In questo modo
	// riesco a sfruttare il fatto di avere il vettore di posizioni one_to_change
	// inutilizzato perché la maschera la faccio con one_to_remove che è uguale.
	// Inoltre non mi serve che siano aggiornati tutti gli one_to_change, ma ne
	// basta uno per ogni gruppetto di shifts precedenti che vado a considerare,
	// si cambia il primo (all'indice 0) che c'è in tutti i casi
	if(this->multiple_shifts[index].size() > 0)
		this->multiple_shifts[index][0].one_to_change = pos_not_covered_yet;
}

void SpacedQmer::SetAllMultipleShift()
{
	// calcolo il gruppo di shift precedenti migliore a regime
	this->multiple_shifts.resize(1);
	this->SetMultipleShifts(0);
	
	// calcolo la lunghezza del transitorio: vedo quale tra gli shift a regime è
	// il più distante: fino a quella distanza ho il transitorio da gestire.
	// Per ogni passo del transitorio devo calcolare il suo migliore gruppo di
	// traslazioni
	size_t furthest_pos = this->multiple_shifts[0][0].shift_min;
	for (size_t i = 1; i < this->multiple_shifts[0].size(); i++)
		furthest_pos = std::max(furthest_pos, this->multiple_shifts[0][i].shift_min);
	
	// aumento la lunghezza del vettore che contiene i gruppi di shift precedente
	// in modo da fargli contenere il numero di shift del trasitorio più quello a
	// regime
	this->multiple_shifts.resize(furthest_pos);

	// calcolo tutti i gruppi di shift per ogni passo del transitorio
	for(size_t k = 1; k < furthest_pos; k++)
		this->SetMultipleShifts(k);

	// Per ciascun gruppo vado a calcolare le maschere per cancellare i valori
	// che non interessano
	this->SetBitMasks();

	//Queste stampe servono a capire come si comporta la funzione di calcolo dell'hash
	//dato il seed corrente. A regime, se non è stato limitato il numero di hash da
	//utilizzare si deve calcolare solamente l'ultima posizione.
	//Per attivare le stampe basta cambiare false in true.
	if (false)
	{
		std::cout << "Seed: " << this->spaced_q << "\nha un transitiorio lungo " << this->multiple_shifts.size() << "\n";
		for (size_t i = 1; i <  this->multiple_shifts.size(); i++)
		{
			if (this->multiple_shifts[i].size() == 0)
				std::cout << "\nGruppo di shift "<< i << " ricalcola tutte le posizioni\n";
			else
			{
				std::cout << "\nGruppo di shift " << i << " utilizza " << this->multiple_shifts[i].size() << " hash\n";
				std::cout << "Numero posizioni da calcolare: " << this->multiple_shifts[i][0].one_to_change.size() << "\n";
				// Decommentando questo si vedono i dettagli sui singoli hash
				/* for(size_t j = 0; j < this->multiple_shifts[i].size(); j++)
				{
					print_shift(multiple_shifts[i][j]);
				} */
			}
		}
		std::cout << "\nGruppo di shift a regime utilizza " << this->multiple_shifts[0].size() << " hash\n";
		std::cout << "Numero posizioni da calcolare: " << this->multiple_shifts[0][0].one_to_change.size() << "\n\n\n";
		// Come sopra
		/* for(size_t j = 0; j < this->multiple_shifts[0].size(); j++)
		{
			print_shift(this->multiple_shifts[0][j]);
		} */
	}
}

void SpacedQmer::SetBitMasks()
{
	for(size_t k = 0; k < this->multiple_shifts.size(); k++)
		for(size_t i = 0; i < this->multiple_shifts[k].size(); i++)
			for (size_t pos : this->multiple_shifts[k][i].one_to_keep)
				this->multiple_shifts[k][i].mask |= (uint64_t)0b11 << (2 * pos);
}



// funzioni di supporto che stampano qualche risultato
void printp(Position p)
{
	for (size_t i = 0; i < p.size(); ++i)
		std::cout << p[i] << " ";
}

void print_shift(PreviousShift s)
{
	std::cout << "\none_to change= ";
	printp(s.one_to_change);
	// vector containing index at which value the hash has to be remove
	std::cout << "\none_to_remove= ";
	printp(s.one_to_remove);

	std::cout << "\none_to_keep= ";
	printp(s.one_to_keep);

	std::cout << "\none_exit= " << s.one_exit;
	std::cout << "\nshift_min= " << s.shift_min;
	std::cout << "\nsize totale= " << s.GetSize() << "\n";
	if (s.mask != 0)
		std::cout << "Maschera calcolata: " << std::bitset<42>(s.mask);
	std::cout << "\n\n";
}

// funzioni di supporto che permettono di gestire con più chiarezza il vettore
// delle posizioni che non sono ancora state recuperate
void deleteElement(Position& pos, size_t index)
{
	// trova l'indice dell'elemento da rimuovere
	auto it = find(pos.begin(), pos.end(), index);
	if (it != pos.end())
		pos.erase(it);
}

bool isContained(Position pointer, Position pos, size_t index)
{
	for (size_t i = 0; i < pointer.size(); i++)
		if (pos[pointer[i]] == index)
			return true;
	
	return false;
}