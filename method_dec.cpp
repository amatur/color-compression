//version: jan 27, 2023, 10:40
#include<cmph.h> //#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sdsl/bit_vectors.hpp>
#include<unordered_map>
#include<sstream>
#include<string>
#include <queue>
#include <map>
#include <climits> // for u_int32_t_BIT
#include <iterator>
#include <algorithm>
#include <iostream>
#include<sstream>
#include "BooPHF.h"
using namespace std;
using namespace sdsl;

#include <unordered_map>

class OutputFile{
	public:
		string filename;
		std::ofstream fs;
	OutputFile(){

	}
	OutputFile(string filename){
		this->filename = filename;
		fs.open (filename.c_str(),  std::fstream::out );
	}
	void init(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::out);
	}
	void write(string towrite){
		fs << towrite; // <<endl;
	}
	~OutputFile(){
		fs.close();
	}
};

class InputFile{
	public:
	string filename;
	std::fstream fs;
	InputFile(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::in);
	}
	InputFile(){
	}
	void init(const std::string filename){
		this->filename=filename;
		this->fs.open(this->filename, fstream::in);
	}
	void rewind(){
		this->fs.close();
		this->fs.open(this->filename, fstream::in);
	}

	~InputFile(){
		fs.close();
	}
};


typedef std::vector<bool> HuffCode;
typedef std::map<u_int32_t, HuffCode> HuffCodeMap;

namespace Huffman{
	/// @brief source rosetta code
	class INode
	{
		public:
			const int f;
			virtual ~INode() {}
		protected:
			INode(int f) : f(f) {}
	};
	class InternalNode : public INode
	{
	public:
		INode *const left;
		INode *const right;

		InternalNode(INode* c0, INode* c1) : INode(c0->f + c1->f), left(c0), right(c1) {}
		~InternalNode()
		{
			delete left;
			delete right;
		}
	};
	class LeafNode : public INode
	{
	public:
		u_int32_t c;

		LeafNode(int f, u_int32_t c) : INode(f), c(c) {}
	};

	struct NodeCmp
	{
		bool operator()(const INode* lhs, const INode* rhs) const { return lhs->f > rhs->f; }
	};

	INode* BuildTree(u_int32_t* frequencies, u_int32_t UniqueSymbols)
	{
		std::priority_queue<INode*, std::vector<INode*>, NodeCmp> trees;
	
		for (u_int32_t i = 0; i < UniqueSymbols; ++i)
		{
			if(frequencies[i] != 0)
				trees.push(new LeafNode(frequencies[i], (u_int32_t)i));
		}
		while (trees.size() > 1)
		{
			INode* childR = trees.top();
			trees.pop();

			INode* childL = trees.top();
			trees.pop();

			INode* parent = new InternalNode(childR, childL);
			trees.push(parent);
		}
		return trees.top();
	}

	void GenerateCodes(const INode* node, const HuffCode& prefix, HuffCodeMap& outCodes)
	{
		if (const LeafNode* lf = dynamic_cast<const LeafNode*>(node))
		{
			outCodes[lf->c] = prefix;
		}
		else if (const InternalNode* in = dynamic_cast<const InternalNode*>(node))
		{
			HuffCode leftPrefix = prefix;
			leftPrefix.push_back(false);
			GenerateCodes(in->left, leftPrefix, outCodes);

			HuffCode rightPrefix = prefix;
			rightPrefix.push_back(true);
			GenerateCodes(in->right, rightPrefix, outCodes);
		}
	}

	u_int32_t HuffDecode(const INode* root, string s)
	{
		string ans = "";
		u_int32_t ansint=0;
		const INode* curr = root;
		for (int i = 0; i < s.size(); i++) {
			if(  const LeafNode* lf  = dynamic_cast<const LeafNode*>(curr) ){
					ansint = (u_int32_t)(lf->c);
					return ansint;
					curr = root;
			}else if(const InternalNode* internal   = dynamic_cast<const InternalNode*>(curr) ){
				if (s[i] == '0')
					curr = internal->left;
				else
					curr = internal->right;

				if(const LeafNode* lf  = dynamic_cast<const LeafNode*>(curr)){
					ansint = (u_int32_t)(lf->c);
					cout<< ansint;
				}
			}else{
				exit(2);
			}
		}
		// cout<<ans<<endl;
		return ansint;
	}
}

using namespace Huffman;


class COLESS_Decompress{
public:
	uint64_t convert_binary_string_to_uint(string& str, int start, int end, int block_sz2){ //convert_binary_string_to_uint
		uint64_t res = 0;
		int block_sz = end - start + 1;
// 		assert(block_sz==block_sz2);
        int i =0;
		for (int64_t j = end; j >= start; j--) {
			if (str[j]=='1') {
				res |= 1 << i;
			}
			i+=1;
		}
    	return res;
	}

	uint64_t read_uint(string& str, uint64_t& b_it, int block_sz){ //convert_binary_string_to_uint
		uint64_t res = 0;
		//int block_sz = end - start + 1;
		uint64_t end = block_sz + b_it - 1;
// 		assert(block_sz==block_sz2);
        uint64_t i = 0;
        uint64_t j = end;

		while(true){
			if (str[j]=='1') {
				res |= 1 << i;
			}
			i+=1;
            if(j!=b_it){
                j--;
            }else{
                break;
            }
		}
		b_it += block_sz;
    	return res;
	}

	COLESS_Decompress(long num_kmers, int M, int C,  string sdsl_file=""){
		int lm = ceil(log2(M));
		int lc = ceil(log2(C));

		if(1==0){

    	    rrr_vector<256> rrr_map;
	        rrr_vector<256> rrr_main;
			// rrr_vector rrr_map = rrr_vector<256>();
			// rrr_main =  rrr_vector<256>();
			load_from_file(rrr_map, "rrr_map");
			load_from_file(rrr_map, "rrr_main");
			std::ofstream out("str_bv_mapping.txt");
			out << rrr_map;
			out.close();
			stringstream ss_rrr_map;
			ss_rrr_map << rrr_map;
			string str_map = ss_rrr_map.str();
		}

		InputFile file_bb_map("bb_map");
		string str_map;
		getline(file_bb_map.fs, str_map);
		file_bb_map.fs.close();
		

		uint64_t b_it =  0;
		for(int i = 0; i < M; i++){
			uint64_t colclass = read_uint(str_map, b_it, C);
			cout<<colclass<<endl;
		}

        //decompress bb_main
        //decompress local hash table : bug is in local hash table


		
		//stringstream ss_cc_map;
		
		// rrr_vector<256> cc_map = rrr_vector<256>();
		// load_from_file(cc_map, index_file);
		// cout<<cc_map<<endl;
		// ss_cc_map<<cc_map;
	}
};


int main (int argc, char* argv[]){
	vector<string> args(argv + 1, argv + argc);
    string dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname; //string tmp_dir;
    int M, C;
	long num_kmers=0;
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: tool -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -t <tmp-dir>" << endl;
            return 0;
        } else if (*i == "-i") {
            dedup_bitmatrix_fname = *++i;
        } else if (*i == "-d") {
            dup_bitmatrix_fname = *++i;
        }else if (*i == "-c") {
            C = std::stoi(*++i);
        }else if (*i == "-m") {
            M = std::stoi(*++i);
        }else if (*i == "-k") {
            num_kmers = std::stol(*++i);
        }else if (*i == "-s") {
            spss_boundary_fname = *++i;
		}
		// else if (*i == "-t") {
        //     tmp_dir  = *++i;
		// }
    }

	COLESS_Decompress cdec(num_kmers, M, C);
	return EXIT_SUCCESS;
}

