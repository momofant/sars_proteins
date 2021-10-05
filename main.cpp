/*
 * main.cpp
 * Count number of potential proteins in genomic sequence
 * Read sequence from file passed as 1st parameter. Minimal protein length (in aminoacids) passes as 2nd parameter.
 * Sequence is treated as protein if starts with "atg" codone and ends with either "taa", "tag", "tga" in the same reading frame (distance % 3 == 0)
 *
 *  Created on: 21 сент. 2021 г.
 *      Author: anna
 */
#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<vector>

std::string readFasta(const char* fileName)
{
	std::ifstream t(fileName);
	std::stringstream buffer;
	buffer << t.rdbuf();
	return buffer.str();
}

bool isStartCodone(const std::string& sequence, int start) const
{
	return sequence.compare(start, 3, "atg") == 0;
}

bool isStopCodone(const std::string& sequence, int start) const
{
	return sequence.compare(start, 3, "tag") == 0 || sequence.compare(start, 3, "taa") == 0 || sequence.compare(start, 3, "tga") == 0;
}

int getNextStop(const std::string& sequence, int start, int len)
{
	int i = start;
	while(i < len)
	{
		if(isStopCodone(sequence, i))
			return i;
		else
			i+=3;
	}
	return 0;
}

void checkAndReadInputs(int argc, char* argv[], std::string* fasta, int* minProteinSize)
{
	if(argc < 2)
	{
		std::cerr << "You must provide at least one argument\n";
		exit(0);
	}
	*fasta = readFasta(argv[1]);
	*minProteinSize = 3;
	if(argc > 2)
	{
		std::string arg = argv[2];
		try
		{
			std::size_t pos;
			*minProteinSize = 3 * std::stoi(arg, &pos);
			if (pos < arg.size())
			{
				std::cerr << "Trailing characters after number: " << arg << '\n';
			}
		}
		catch (std::invalid_argument const &ex)
		{
			std::cerr << "Invalid number: " << arg << '\n';
		}
		catch (std::out_of_range const &ex)
		{
			std::cerr << "Number out of range: " << arg << '\n';
		}
	}
}

void prepareSequence(std::string* fastaPtr)
{
	(*fastaPtr).erase(std::remove((*fastaPtr).begin(), (*fastaPtr).end(), '\n'), (*fastaPtr).end());
	std::transform((*fastaPtr).begin(), (*fastaPtr).end(), (*fastaPtr).begin(), ::tolower);
}

int main (int argc, char* argv[])
{
	std::string fasta;
	int minProteinSize;
	checkAndReadInputs(argc, argv, &fasta, &minProteinSize);
	prepareSequence(&fasta);
	int i = 0;
	int maxCheckChar = fasta.size() - 2;
	std::vector<int> proteinStart, proteinEnd;
	while( i < maxCheckChar )
	{
		if(isStartCodone(fasta, i))
		{
			int j = getNextStop(fasta, i+3, maxCheckChar);
			if(j==0) //end of string, stop not found
			{
				break;
			}
			else
			{
				if(j - i + 1 > minProteinSize)
				{
					proteinStart.push_back(i);
					proteinEnd.push_back(j);
				}
				++i;
			}
		}
		else
			++i;
	}

	//Output the result
	auto ie = begin(proteinEnd);
	int cnt = 1;
	std::cout << "Number of proteins found: " << proteinStart.size() << "\n";
	for (auto it = begin (proteinStart); it != end (proteinStart); ++it)
	{
	    std::cout << "Protein " << cnt++ << ": " << "(pos " << *it << " len "<< (*ie - *it + 1)/3 << ").." << fasta.substr(*it, *ie - *it + 3) << "..\n";
	    ++ie;
	}
	return 0;
}
