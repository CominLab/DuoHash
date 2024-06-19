// Compile with  g++-13 -std=c++0x -O3 -fopenmp -I../include -L../build -lDuoHash -o main main.cpp

#include <DuoHash.h>
#include <chrono>
#include <iomanip>


int main()
{
	std::string dataset = "./TestInputFile/paired.txt";
	std::vector<std::string> sequences;
	if (!loadFile(dataset, sequences))
		std::exit(1);

	std::string seedset = "./Seeds/W22L31.fna";
	std::vector<SpacedQmer> multi_spaced;
	{
		std::vector<std::string> tmp;
		if (!loadFile(seedset, tmp))
			std::exit(1);

		for (auto seed : tmp)
			multi_spaced.push_back(SpacedQmer(seed, 0));
	}




	std::cerr << "\nTest with single seed\n";
	DuoHash test(sequences, multi_spaced[0]);
	std::vector<Hash_V> hashes;

	// GetEncoding_ISSH  /  get only forward & reverse hashing
	{
		std::cerr << "GetEncoding_ISSH... ";
		
		auto t1 = std::chrono::steady_clock::now();
		test.GetEncoding_ISSH(hashes, getHashes);
		auto t2 = std::chrono::steady_clock::now();
		
		auto t = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		std::cerr << "OK. Executed in " << (double)t.count() / 1000 << "ms\n";
	}

	// PrintFASTA
	{
		std::cerr << "PrintFASTA... ";

		auto t1 = std::chrono::steady_clock::now();
		test.PrintFASTA(hashes, "single");
		auto t2 = std::chrono::steady_clock::now();

		auto t = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		std::cerr << "OK. Executed in " << (double)t.count() / 1000 << "ms\n";
	}


	


	std::cerr << "\nTest with multiple seed\n";
	DuoHash_multi test_multi(sequences, multi_spaced);
	std::vector<Hash_V_V> hashes_multi;

	// GetEncoding_MISSH_v1
	{
		std::cerr << "GetEncoding_MISSH_v1... ";

		auto t1 = std::chrono::steady_clock::now();
		test_multi.GetEncoding_MISSH_v1(hashes_multi, getBoth);
		auto t2 = std::chrono::steady_clock::now();
		
		auto t = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		std::cerr << "OK. Executed in " << (double)t.count() / 1000 << "ms\n";
	}

	// PrintFASTA
	{
		std::cerr << "PrintFASTA... ";

		auto t1 = std::chrono::steady_clock::now();
		test_multi.PrintFASTA(hashes_multi, "multi");
		auto t2 = std::chrono::steady_clock::now();

		auto t = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
		std::cerr << "OK. Executed in " << (double)t.count() / 1000 << "ms\n";
	}
}