//hello
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
using namespace std;
using namespace sdsl;


//sort -T=~/s/tmp/ export TMPDIR=/tmp
//position uint64_t


	// uint64_t write_block(int block_sz, uint8_t category, int value, vector<int>& positions){
	// 	uint64_t b_it = 0;
	// 	write_number_at_loc(positions, num, block_size, b_it);

	// 	//write category
	// 	if(category==0){ //either log M or log U
	// 		write_number_at_loc(positions, category, 1, b_it);
	// 	}else if(category==1){ //log C 
	// 		write_number_at_loc(positions, category, 2, b_it);
	// 	}if(category==2){ // RUN of 1
	// 		write_number_at_loc(positions, category, 2, b_it);
			
	// 	}
	// 	return b_it; // at the end b_it equals size of vector
	// }
class OutputFile{
	public:
		string filename;
		std::ofstream fs;
	OutputFile(){

	}
	OutputFile(string filename){
		this->filename = filename;
		fs.open (filename,  std::fstream::out );
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
class LogFile : public OutputFile	//derived class
{
	public:
		void log(string param_name, string param_value, string delim=":")
		{
			fs << param_name << delim << param_value << endl;
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




namespace CMPH{
    void create_cmph_table(string key_filename, cmph_t * hash){
		cmph_io_adapter_t *source;
		FILE * keys_fd; 

        keys_fd = fopen(key_filename.c_str(), "r"); //Open file with newline separated list of keys
        hash = NULL;
        if (keys_fd == NULL) 
        {
          fprintf(stderr, ("File "+key_filename+" not found\n").c_str());
          exit(1);
        }	
        // Source of keys
        source = cmph_io_nlfile_adapter(keys_fd);
        cmph_config_t *config = cmph_config_new(source);
        cmph_config_set_algo(config, CMPH_BDZ);
        hash = cmph_new(config);
        //cmph_config_destroy(config);
		//fclose(keys_fd);
		//cmph_io_nlfile_adapter_destroy(source);   
    }

    unsigned int lookup(cmph_t *hash, string key_str){ //Find key
       const char *key = key_str.c_str();
	   cout<<key<<endl;
       unsigned int id = cmph_search(hash, key, (cmph_uint32)strlen(key));
       return id;
    }

    // ~CMPH(){
    //   //Destroy hash
    //   cmph_destroy(hash);
    //   cmph_io_nlfile_adapter_destroy(source);   
    //   
    // }
};
using namespace Huffman;
using namespace CMPH;

class COLESS{
public:
	InputFile dedup_bitmatrix_file, spss_boundary_file, dup_bitmatrix_file, tmp_dir;
	LogFile logfile_main;
	long num_kmers;
	int M;
	int C;
	cmph_t* cmp_hash_ptr;
	const int max_run = 16;
	vector<uint64_t> positions;
	HuffCodeMap huff_code_map;
	enum category { CATEGORY_RUN=3, CATEGORY_COLCLASS=0, CATEGORY_COLVEC=2 };
	int lm = 0;
	int lc = 0;

	COLESS(long num_kmers, int M, int C, string dedup_bitmatrix_fname, string dup_bitmatrix_fname, string spss_boundary_fname){
		dedup_bitmatrix_file.init(dedup_bitmatrix_fname);
		spss_boundary_file.init(spss_boundary_fname);
		dup_bitmatrix_file.init(dup_bitmatrix_fname);
		this->num_kmers = num_kmers;
		this->M = M;
		this->C = C;
		this->lm = ceil(log2(M));
		this->lc = ceil(log2(C));
		logfile_main.init("log_coless");
	}


	unsigned int lookup(string key_str){ //Find key
       const char *key = key_str.c_str();
       unsigned int id = cmph_search(cmp_hash_ptr, key, (cmph_uint32)strlen(key));
       return id;
    }

	int hammingDistance (uint64_t x, uint64_t y) {
		uint64_t res = x ^ y;
		return __builtin_popcountll (res) ;
	}

	float get_average(vector<uint64_t> v){
		if(v.size()==0){
			return 0;
		}
		uint64_t summ = 0;
		for (uint64_t e:  v){
			summ+=e;
		}
		return summ/1.0/v.size();
	}

	float get_average(vector<int> v){
		if(v.size()==0){
			return 0;
		}
		uint64_t summ = 0;
		for (uint64_t e:  v){
			summ+=e;
		}
		return summ/1.0/v.size();
	}
	void write_number_at_loc_advanced_by_block_sz(vector<uint64_t> & positions, uint64_t num, uint64_t loc_advanced_by_block_sz){ //requires loc_advanced_by_block_sz += block_size; 
		int64_t j=0;
		while(num!=0)
		{
			if(num%2 == 1){
				positions.push_back(loc_advanced_by_block_sz-1+j); //b[loc_advanced_by_block_sz-1+j] = num%2;
			}
			j--;
			num /= 2;
			
		}
	}

	void write_number_at_loc(vector<uint64_t> & positions, uint64_t num, int block_size, uint64_t& b_it ){
		write_number_at_loc_advanced_by_block_sz(positions, num, b_it+block_size);
		b_it += block_size; //successfully written and place on next bit; if size is 2, (0,1) written, now val is 2.
	}

	void write_binary_string_at_loc(vector<uint64_t> & positions, string binarystring, uint64_t& b_it){
		for (size_t i = 0; i< binarystring.length(); i++) {
			if (binarystring[i]=='1'){
				positions.push_back(b_it+i);
			}
		}
		b_it += binarystring.length();
	}


	bit_vector store_as_sdsl(vector<uint64_t>& positions, uint64_t bv_size, string filename){
		//bit_vector bv = bit_vector(bv_size, 0);
		bit_vector bv(bv_size, 0);
		for (uint64_t p: positions){
			bv[p] = 1;
		}
		rrr_vector<256> rrr_bv(bv);
		//cout << "rrr_MB_bv_mapping="<<size_in_bytes(rrr_bv_mapping)/1024.0/1024.0 << endl;
		store_to_file(rrr_bv, filename);	//"rrr_bv_mapping.sdsl"
		return bv;
	}

	void store_as_binarystring(vector<uint64_t>& positions, uint64_t bv_size, string filename){
		OutputFile binarystring_file(filename);
		sort(positions.begin(), positions.end());
		uint64_t bvi = 0;
		for (uint64_t k = 0; k<bv_size; k++){
			if(bvi < positions.size()){
				if(positions[bvi]==k){
					binarystring_file.write("1");
					bvi++;
				}else{
					binarystring_file.write("0");
				}
			}else{
				binarystring_file.write("0");
			}	
		}
	}

	void store_global_color_class_table(){
		vector<uint64_t> positions;  //wasteful
		uint64_t b_it = 0;

		uint64_t* array_hi = new uint64_t[M];	// maintaing upto C/2 bits
		uint64_t* array_lo = new uint64_t[M];	// maintaing upto C/2 bits

		for(int x=0; x<M; x++){
			string bv_line;
			getline(dedup_bitmatrix_file.fs, bv_line);
			unsigned int idx = lookup(bv_line);		// returns an if in range (0 to M-1)
			array_hi[idx] = std::stoull(bv_line.substr(0,std::min(64,int(C))), nullptr, 2) ;

			array_lo[idx]=0;
			if(C > 64){
				string ss=bv_line.substr(64,(C-64));
				array_lo[idx]=std::stoull(ss, nullptr, 2);
			}
		}
		cout << "Expected_MB_bv_mapping="<<(C*M)/8.0/1024.0/1024.0 << endl;
		dedup_bitmatrix_file.fs.close();

		//	write_number_at_loc(positions, array_hi[idx], min(64, C), b_it ); //array_hi[x] higher uint64_t
		//		write_number_at_loc(positions, array_lo[idx], C-64, b_it ); //array_lo[x] lower uint64_t
		//
		store_as_binarystring(positions, b_it, "bb_map" );
		cout << "expected_MB_bv_mapping="<<(C*M)/8.0/1024.0/1024.0 << endl;
		cout << "rrr_MB_bv_mapping="<<size_in_bytes(store_as_sdsl(positions, b_it, "rrr_bv_mapping.sdsl" ))/1024.0/1024.0 << endl;
	}


	void method1_pass0(bool skip = false){ //load_huffman_table();//int M -> variable string   //writehuffman[in
		// void get_freq_count(){ // scan through all the color vectors to get freq count, global and local
		// // table of size M 
		// }

		double t_begin,t_end; struct timeval timet;
		printf("Construct a MPHF with  %lli elements  \n",M);
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
		create_cmph_table(dedup_bitmatrix_file.filename, cmp_hash_ptr);
		 
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
		double elapsed = t_end - t_begin;
		printf("CMPH constructed perfect hash for %llu keys in %.2fs\n", M,elapsed);

		//if(!skip){
			gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
			OutputFile cmp_keys("cmp_keys");  // get frequency count
			for (uint64_t i=0; i < num_kmers; i+=1){
				string bv_line;
				getline (dup_bitmatrix_file.fs,bv_line);
				cmp_keys.fs<< lookup(bv_line)<<endl;
			}
			gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
			printf("CMPH lookup for %llu keys in %.2fs\n", num_kmers, M,t_end - t_begin);
		//}
		
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
		system("cat cmp_keys | sort -n | uniq -c | rev | cut -f 2 -d\" \" | rev > frqeuency_sorted");
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
		printf("Sorting and getting frequencies for %llu keys in %.2fs\n", num_kmers, M,t_end - t_begin);	

		//
		InputFile infile_freq("frqeuency_sorted");
		string line;
		// Build frequency table
		u_int32_t *frequencies = new u_int32_t[M]; // M -> no. of unique symbols
		std::fill_n(frequencies, M, 0);
		u_int32_t i = 0;
		while(getline(infile_freq.fs, line)){
			stringstream ss(line);
			u_int32_t a ;
			ss >> a; 
			frequencies[i++]= a;
		}		
		INode* root = BuildTree(frequencies, M);
        GenerateCodes(root, HuffCode(), huff_code_map); // huff_code_map is filled: uint32t colclassid-> vector bool
		delete frequencies;
		delete root;
	}

	void method1_pass1(bool skip_pass = false){
		if(skip_pass) return;

		store_global_color_class_table();
		// bit vector values
		uint64_t b_it=0;
		vector<uint64_t> positions; // positions for main vector

		vector<char> spss_boundary; 
		for (uint64_t i=0; i < num_kmers; i+=1){ //load spss_boundary vector in memory from disk
			string spss_line;
			getline (spss_boundary_file.fs,spss_line); 
			spss_boundary.push_back(spss_line[0]); //this kmer starts a simplitig
		}

		//per simplitig values
		set<uint32_t> local_col_classes_uniq; //get the bool HuffCodeMap[M-1]
		int use_local_hash = 0;
		int use_local_hash_nonrun = 0;
		int use_local_hash_huff = 0;
		int use_local_hash_huff_nonrun = 0;
		uint64_t skip=0;
		int case_run = 0;
		int case_lm = 0;
		int case_nonrun = 0;
		int case_dlc = 0;
		uint64_t sum_length_huff = 0;
		uint64_t sum_length_huff_nonrun = 0;
		uint64_t sum_length_huff_uniq = 0;	
		uint64_t sum_length_huff_uniq_nonrun = 0;
		
		//per kmer values
		uint64_t num_kmer_in_simplitig = 0;
		uint64_t curr_bv_hi = 0;
		uint64_t curr_bv_lo = 0;
		uint64_t prev_bv_hi = 0;
		uint64_t prev_bv_lo = 0;

		//pass 1: collect if local or not 
		for (uint64_t i=0; i < num_kmers; i+=1){ 
			//load the color vector of current k-mer from disk to "curr_bv_hi/lo"
			string bv_line;
			getline (dup_bitmatrix_file.fs,bv_line); // bv line = color vector C bits
			curr_bv_hi = std::stoull(bv_line.substr(0,std::min(64, C)), nullptr, 2);
			curr_bv_lo = 0;
			if(C > 63){
				curr_bv_lo = std::stoull(bv_line.substr(64,bv_line.length()-64), nullptr, 2);
			} 

			//per kmer task
			num_kmer_in_simplitig+=1;  //start of simplitig id: num_kmer_in_simplitig
			
			unsigned int curr_kmer_cc_id = lookup(bv_line); //uint64_t num = bphf->lookup(curr_bv);
				
			if(spss_boundary[i]=='0'){ // non-start
				int hd_hi = hammingDistance(prev_bv_hi, curr_bv_hi);
				int hd_lo = hammingDistance(prev_bv_lo, curr_bv_lo);
				int hd= hd_hi+hd_lo;
				if(hd==0){
					skip+=1;	
					case_run+=1;	
				}else{
					if(skip!=0){ 	//not skipped, write lm
						write_number_at_loc(positions, CATEGORY_RUN, 2, b_it);
						write_number_at_loc(positions, skip, 1+floor(log2(skip)), b_it);
					}
					skip=0;


					sum_length_huff_nonrun += huff_code_map[curr_kmer_cc_id].size();
					if(hd*lc < lm){
						case_dlc += 1;
						cout<<"hd: "<<hd<<endl;
					}else{
						case_lm += 1;
						sum_length_huff += huff_code_map[curr_kmer_cc_id].size();
					}
					
				}
			}else{	//start of simplitig, so logm_case
				case_lm+=1;
				// if(skip!=0){ 	//not skipped, write lm
				// 	write_number_at_loc(positions, skip, 1+floor(log2(skip)), b_it);
				// }
				skip=0;
				
				sum_length_huff += huff_code_map[curr_kmer_cc_id].size();
				sum_length_huff_nonrun += huff_code_map[curr_kmer_cc_id].size();

				local_col_classes_uniq.insert(curr_kmer_cc_id);
				// write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
			}

			if(spss_boundary[(i+1)%num_kmers]=='1'){	// end k-mer of simplitig
				// decide what to do
				int l = local_col_classes_uniq.size(); //case_lm
				int ll = ceil(log2(l));

				for(uint32_t uniq_col_class_id: local_col_classes_uniq){
					sum_length_huff_uniq += huff_code_map[uniq_col_class_id].size();
					sum_length_huff_uniq_nonrun += huff_code_map[uniq_col_class_id].size();
				}

				use_local_hash = ( (ll - lm ) * case_lm + lm * (1+l) < 0 ) ;  //ll*case_lm + (lm + l*lm) ::: lm * case_lm 
				use_local_hash_nonrun = ( (ll - lm ) * case_nonrun + lm * (1+l) < 0 ) ;  //ll*case_lm + (lm + l*lm) ::: lm * case_lm 
				use_local_hash_huff = ( (ll*case_lm - sum_length_huff + lm + sum_length_huff_uniq) < 0);
				use_local_hash_huff_nonrun = ( ll*case_nonrun - sum_length_huff_nonrun + lm + sum_length_huff_uniq_nonrun < 0 );

				logfile_main.fs<<use_local_hash<<" "<<use_local_hash_nonrun<<" "<<use_local_hash_huff<<" "<<use_local_hash_huff_nonrun<<" "<<num_kmer_in_simplitig<<endl;
				
				//re-init for new simplitig
				//vector<uint32_t>().swap(local_col_classes_uniq);//
				local_col_classes_uniq.clear();
				num_kmer_in_simplitig = 0;

				
				use_local_hash = use_local_hash_nonrun = use_local_hash_huff = use_local_hash_huff_nonrun = 0;
				skip=0;
				case_run = case_lm = case_nonrun = case_dlc = 0;
				sum_length_huff = sum_length_huff_nonrun = sum_length_huff_uniq = sum_length_huff_uniq_nonrun =  0;

			}
		}
	}
};

	// void method1_pass2(){
	// 	//try this way
	// 	write_number_at_loc(positions, CATEGORY_COLCLASS, 1, b_it);
	// 	write_binary_string_at_loc(positions,  huff(M), b_it_ptr);
	// 	for (uint64_t i=0; i < num_kmers; i+=1){ // For each k-mer in union kmer set
	// 		string bv_line;
	// 		getline (dup_bitmatrix_file.fs,bv_line); // bv line = color vector C bits
	// 		curr_bv = bv_line;
	// 		num_kmer_in_simplitig+=1;  //start of simplitig id: num_kmer_in_simplitig
			

	// 		curr_bv_hi = std::stoull(bv_line.substr(0,std::min(64, C)), nullptr, 2);
	// 		curr_bv_lo = 0;
	// 		if(C > 63){
	// 			curr_bv_lo = std::stoull(bv_line.substr(64,bv_line.length()-64), nullptr, 2);
	// 		} 
			
	// 		if(i!=0){ // non-start
	// 			int hd_hi = hammingDistance(prev_bv_hi, curr_bv_hi);
	// 			int hd_lo = hammingDistance(prev_bv_lo, curr_bv_lo);
	// 			int hd= hd_hi+hd_lo;
	// 			if(hd*lc < lm){
	// 				skip+=1;
	// 				uint64_t A_hi = prev_bv_hi;
	// 				uint64_t B_hi = curr_bv_hi;
	// 				uint64_t A_lo = prev_bv_lo;
	// 				uint64_t B_lo = curr_bv_lo;

	// 				if(hd!=0){
	// 					//category 1
	// 					if(hd<4){
	// 						for (int i_bit=0; i_bit < lc; i_bit+=1){
	// 							if ((( A_hi >>  i_bit) & 1) != (( B_hi >>  i_bit) & 1)){ // check if the bit at the 'i'th position is different
	// 								//cout<<"different: "<<i_bit<<endl;
	// 								write_number_at_loc(positions, i_bit, lc, b_it);
	// 							}
	// 						}
	// 						for (int i_bit=0; i_bit < lc; i_bit+=1){
	// 							if ((( A_lo >>  i_bit) & 1) != (( B_lo >>  i_bit) & 1)){
	// 								//cout<<"different: "<<i_bit<<endl;
	// 								write_number_at_loc(positions, i_bit, lc, b_it);
	// 							}
	// 						}
	// 					}
	// 				}
					
	// 			}else{
	// 				num_logm_case+=1;

	// 				if(skip!=0){ 	//not skipped, write lm
	// 					write_number_at_loc(positions, skip, 1+floor(log2(skip)), b_it);
	// 				}
	// 				skip=0;
	// 				unsigned int col_class_id = cmp.lookup(bv_line); //uint64_t num = bphf->lookup(curr_bv);
	// 				local_col_classes_uniq.insert(col_class_id);
	// 				local_col_classes_nonuniq.push_back(col_class_id);
	// 				write_number_at_loc(positions, col_class_id, lm, b_it);
	// 			}
				
	// 		}else{
	// 			num_logm_case+=1;
	// 			num_logm_case_global+=1;

	// 			b_it += lm;
	// 			unsigned int col_class_id = cmp.lookup(bv_line);//uint64_t num = bphf->lookup(curr_bv);
	// 			local_col_classes_uniq.insert(col_class_id);
	// 			local_col_classes_nonuniq.push_back(col_class_id);
	// 			write_number_at_loc(positions, col_class_id, lm, b_it);
	// 		} 
	// 		prev_bv_hi = curr_bv_hi;
	// 		prev_bv_lo = curr_bv_lo;

	
	// 		if (spss_boundary[(i+1)%num_kmers]=='1') {   //|| i==num_kmers-1	//potential_bug
	// 			if(local_col_classes_uniq.size()!=0)
	// 				o_caseuniqm.fs << local_col_classes_uniq.size()<<" "<<local_col_classes_nonuniq.size()<<" "<<local_col_classes_uniq.size()*lm - local_col_classes_nonuniq.size()*(lm - ceil(log2(local_col_classes_uniq.size())))<<" "<<get_average(local_runs_of_0)<<" "<<local_runs_of_0.size()<<" "<<num_kmer_in_simplitig<<endl;
	// 			else
	// 				o_caseuniqm.fs << local_col_classes_uniq.size()<<" "<<local_col_classes_nonuniq.size()<<" "<<0<<" "<<get_average(local_runs_of_0)<<" "<<local_runs_of_0.size()<<" "<<num_kmer_in_simplitig<<endl;

	// 			local_col_classes_uniq.clear();
	// 			local_col_classes_nonuniq.clear();
	// 		}
	// 	}

	// 	uint64_t num_bits = b_it;
	// 	cout<<num_bits<<" bits in method1"<<endl;
	// 	bit_vector b = bit_vector(num_bits, 0);
	// 	cout<<"init success"<<endl;

	// 	for (uint64_t p : positions)
	// 	{
	// 		b[p]=1;
	// 	}

	// 	cout<<"Num1s_vector_main="<<positions.size()<<endl;
	// 	cout<<"before_rrr_MB_vector_main="<<size_in_bytes(b)/1024.0/1024.0<<endl;
	// 	rrr_vector<256> rb(b);
	// 	float rb_bytes = size_in_bytes(rb);
	// 	float totsize= rb_bytes
	// 	cout << "rrr size, without MPHF (MB)= "<< totsize/1024.0/1024.0  << endl;
	// 	//cout << "rrr size, all (MB) = "<< (size_in_bytes(rrr_bv_mapping)+totsize)/1024.0/1024.0 << endl;

	// 	logfile_main.log("rrr_MB_vector_main",rb_bytes/1024.0/1024.0);
	// 	store_to_file(rb, "rrrbv_1.sdsl");
	// }


// namespace BitManip
// {
// 	uint64_t bitsToShort(char* bits, int unit=64) {
// 		int j = unit;
// 		uint64_t res = 0;
// 		for (int j = unit; j > 0; j--) {
// 			int i = unit-j;
// 			if (bits[j]=='1') {
// 				res |= 1 << i;
// 			}
// 		}
//     	return res;
// 	}
// 	bool* shortToBits(short value) {
// 		bool* bits = new bool[64];
// 		int count = 0;
// 		while(value) {
// 			if (value&1)
// 				bits[count] = 1;
// 			else
// 				bits[count] = 0;
// 			value>>=1;
// 			count++;
// 		}
// 		return bits;
// 	}

// 	void write_bits_from_binary_string(std::ostream & output, std::string const & input)
// 	{
// 		unsigned char c = 0;
// 		int bits = 0;

// 		for (auto i = input.begin(); i != input.end(); ++i) {
// 			if (*i == '0' || *i == '1') {
// 				c = (c << 2);
// 				if (*i == '1') {
// 					++c;
// 				}

// 				if (++bits == 8) {
// 					output << c;
// 					c = 0;
// 					bits = 0;
// 				}
// 			}
// 		}

// 		if (bits > 0) {
// 			while (bits < 8) {
// 				c <<= 2;
// 				++bits;
// 			}
// 			output << c;
// 		}
// 	}
// }


int main (int argc, char* argv[]){
	vector<string> args(argv + 1, argv + argc);
    string dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname;
	//string tmp_dir;
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

	COLESS coless(num_kmers, M, C, dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname);
	
	coless.method1_pass0(true);
	coless.method1_pass1();

	return EXIT_SUCCESS;
}

