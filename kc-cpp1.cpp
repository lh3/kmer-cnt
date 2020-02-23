#include <string>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <cstdlib>
#include "ketopt.h"

typedef std::unordered_map<std::string, int> counter_t;

static void count_seq(counter_t &h, int k, std::string &seq)
{
	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	for (size_t i = 0; i <= seq.length() - k; ++i) {
		auto kmer_for = seq.substr(i, k);
		if (kmer_for.find('N') < kmer_for.length()) continue;
		auto kmer_rev = kmer_for;
		std::reverse(kmer_rev.begin(), kmer_rev.end()); // reverse
		for (size_t j = 0; j < kmer_rev.length(); ++j) { // complement
			if (kmer_rev[j] == 'A') kmer_rev[j] = 'T';
			else if (kmer_rev[j] == 'T') kmer_rev[j] = 'A';
			else if (kmer_rev[j] == 'C') kmer_rev[j] = 'G';
			else if (kmer_rev[j] == 'G') kmer_rev[j] = 'C';
		}
		++h[kmer_for < kmer_rev? kmer_for : kmer_rev];
	}
}

static void count_file(counter_t &h, const char *fn, int k)
{
	std::ifstream file(fn);
	std::string line, seq;
	while (getline(file, line)) {
		if (line[0] == '>') {
			if (seq.length() > 0) 
				count_seq(h, k, seq);
			seq = "";
		} else seq += line;
	}
	if (seq.length() > 0)
		count_seq(h, k, seq);
	file.close();
}

static void print_hist(const counter_t &h)
{
	uint64_t cnt[256];
	for (int i = 0; i < 256; ++i) cnt[i] = 0;
	for (auto &k : h) {
		int c = k.second < 255? k.second : 255;
		++cnt[c];
	}
	for (int i = 1; i < 256; ++i)
		std::cout << i << '\t' << cnt[i] << std::endl;
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, k = 31;
	while ((c = ketopt(&o, argc, argv, 1, "k:", 0)) >= 0) {
		if (c == 'k') k = ::atoi(o.arg);
	}
	if (argc - o.ind == 0) {
		std::cout << "Usage: kc-cpp1 [-k " << k << "] <in.fa>" << std::endl;
		return 1;
	}
	counter_t h;
	count_file(h, argv[o.ind], k);
	print_hist(h);
	return 0;
}
