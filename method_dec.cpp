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



namespace CMPH{
	cmph_t *hash_cmph = NULL;
	void create_table(string filename, int nelem ){
		FILE * keys_fd = fopen(filename.c_str(), "r");
		
		if (keys_fd == NULL) 
		{
		fprintf(stderr, "File not found\n");
		exit(1);
		}	
		// Source of keys
		cmph_io_adapter_t *source = cmph_io_nlfile_adapter(keys_fd);
	
		cmph_config_t *config = cmph_config_new(source);
		cmph_config_set_algo(config, CMPH_CHM);
		hash_cmph = cmph_new(config);
		cmph_config_destroy(config);
		
		cmph_io_nlfile_adapter_destroy(source);   
		fclose(keys_fd);
	}

	unsigned int lookup(string str){	
		const char *key = str.c_str(); 
		//Find key
		unsigned int id = cmph_search(hash_cmph, key, (cmph_uint32)strlen(key));
		// fprintf(stderr, "Id:%u\n", id);
		//Destroy hash
		//cmph_destroy(hash);
		return id;
	}

	void mphf_destroy(){
		cmph_destroy(hash_cmph);
	}
}
using namespace CMPH;


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
    int max_run = 16;
    int lmaxrun = 4;
    string* global_table;
    int C;
    int M;
    int lm, lc;
    OutputFile dec_ess_color;

    COLESS_Decompress(long num_kmers, int M, int C)
    {
        this->C = C;
        this->M = M;
        lm = ceil(log2(M));
        lc = ceil(log2(C));
        global_table = new string[M];
        dec_ess_color.init("dec_ess_color");
    }

    uint64_t convert_binary_string_to_uint(string &str, int start, int end, int block_sz2)
    { // convert_binary_string_to_uint
        uint64_t res = 0;
        int block_sz = end - start + 1;
        // 		assert(block_sz==block_sz2);
        int i = 0;
        for (int64_t j = end; j >= start; j--)
        {
            if (str[j] == '1')
            {
                res |= 1 << i;
            }
            i += 1;
        }
        return res;
    }

    char read_one_bit(string& str, uint64_t& b_it){ //convert_binary_string_to_uint
        return str[b_it++];
    }

    int read_number_encoded_in_unary_zero(string& str, uint64_t& b_it){ //convert_binary_string_to_uint
        int length = 0;
        while(str[b_it++]=='0'){
            length+=1;
        }
        return length;
    }
    int read_number_encoded_in_unary_one(string& str, uint64_t& b_it){ //convert_binary_string_to_uint
        int length = 0;
        while(str[b_it++]=='1'){
            length+=1;
        }
        return length;
    }

    
    string read_color_vector(string& str, uint64_t& b_it){
        string col_vec = str.substr(b_it, C);
        b_it+=C;
        return col_vec;
    }

    void flip_bit(string& s, int pos){
        if(s[pos] == '1')   s[pos]='0';
        if(s[pos] == '0')   s[pos]='1';
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

    void run()
    {
        OutputFile color_global("color_global");
       
        if (1 == 0)
        {
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

        uint64_t b_it = 0;
        for (int i = 0; i < M; i++)
        {
            string col_vector = read_color_vector(str_map, b_it);
            color_global.fs << col_vector << endl;
            global_table[i] = col_vector;
        }
        color_global.fs.close();
        create_table(color_global.filename, M);

        b_it = 0;
        // decompress bb_main: only logm

        vector<int> differ_run;

        string last_col_vector = "";
        while (b_it < str_map.length())
        {
            char c = read_one_bit(str_map, b_it);
            if (c == '0')
            {
                if (differ_run.size())
                {
                    for (int d : differ_run)
                    {
                        flip_bit(last_col_vector, d);
                        dec_ess_color.fs << last_col_vector << endl;
                    }
                    differ_run.clear();
                }         
                uint64_t col_class = read_uint(str_map, b_it, lm);
                cout<<col_class<<" "<<"col_class"<<endl;
                last_col_vector = global_table[col_class];
                dec_ess_color.fs << last_col_vector << endl;
            }
            if (c == '1')
            {
                char c2 = read_one_bit(str_map, b_it);
                if (c2 == '1')
                { // run
                    if (differ_run.size() != 0)
                    {
                        for (int d : differ_run)
                        {
                            flip_bit(last_col_vector, d);
                            dec_ess_color.fs << last_col_vector << endl;
                        }
                        differ_run.clear();
                    }
                    int q = read_number_encoded_in_unary_one(str_map, b_it);
                    assert(read_one_bit(str_map, b_it) == '0');
                    int rem = read_uint(str_map, b_it, lmaxrun);
                    int skip = q * max_run + rem;
                    while (skip)
                    {
                        dec_ess_color.fs << last_col_vector << endl;
                        skip--;
                    }
                }
                else
                { // lc 10
                    int differing_bit = read_uint(str_map, b_it, lc);
                    differ_run.push_back(differing_bit);
                }
            }
        }
    }

    // decompress local hash table : bug is in local hash table
    // stringstream ss_cc_map;
    //  rrr_vector<256> cc_map = rrr_vector<256>();
    //  load_from_file(cc_map, index_file);
    //  cout<<cc_map<<endl;
    //  ss_cc_map<<cc_map;
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
    cdec.run();
	return EXIT_SUCCESS;
}

