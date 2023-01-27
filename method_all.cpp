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



class Hashtable {
    std::unordered_map<uint64_t, uint64_t> htmap; // m_to_l
	uint64_t curr_id = 0;

public:
	Hashtable(){
		curr_id = 0;
	}

    uint64_t put_and_getid(uint64_t key) {
		if(htmap.count(key) > 0){ // present
			return htmap[key];
		}  else {	// absent
			htmap[key] = curr_id;
			curr_id+=1;
			return curr_id-1; 
		}
    }

    // const void *get(int key) {
    //         return htmap[key];
    // }

	bool exists(int key){
		return htmap.count(key) > 0;
	}

	void clear(){
		htmap.clear();
		curr_id = 0;
	}

};

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


namespace BPHF{

    typedef boomphf::SingleHashFunctor<uint64_t>  hasher_t;
    typedef boomphf::mphf<  uint64_t, hasher_t  > boophf_t;
		 boophf_t * bphf; 
    // void construct_bphf_table( int *& data, int nelem, boophf_t * &bphf ){
    //     int nthreads = 8;
    //     double t_begin,t_end; struct timeval timet;
    //     printf("Construct a BooPHF with  %lli elements  \n",nelem);
    //     gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
    //     auto data_iterator = boomphf::range(static_cast<const int*>(data), static_cast<const int*>(data+nelem));
    //     double gammaFactor = 7.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
    //     bphf = new boomphf::mphf<int,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
    //     gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
    //     printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,t_end - t_begin);
    // }

	void create_table(string filename, int nelem ){
		InputFile infile(filename);
		uint64_t* data = (uint64_t * ) calloc(nelem,sizeof(uint64_t));
		string bv_line;
		int i = 0;
		while (getline(infile.fs,bv_line )){
			data[i++] = std::stoull(bv_line, nullptr, 2) ;
		}
		int nthreads = 8;
        double t_begin,t_end; struct timeval timet;
        printf("Construct a BooPHF with  %lli elements  \n",nelem);
        gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
        auto data_iterator = boomphf::range(static_cast<const uint64_t*>(data), static_cast<const uint64_t*>(data+nelem));
        double gammaFactor = 7.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
        bphf = new boomphf::mphf<uint64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
        gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
        printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,t_end - t_begin);
	}

	unsigned int lookup(string str){	
		return bphf->lookup(std::stoull(str, nullptr, 2));
	}

	void mphf_destroy(){
		delete bphf;
	}
}   
//using namespace BPHF; 

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

// class CMPH{
//     public:
//       cmph_t *hash;
//       cmph_io_adapter_t *source;
//       FILE * keys_fd; 

// 	// CMPH(){
// 	// }
//     CMPH(string key_filename){
//         keys_fd = fopen(key_filename.c_str(), "r"); //Open file with newline separated list of keys
//         hash = NULL;
//         if (keys_fd == NULL) 
//         {
//           fprintf(stderr, ("File "+key_filename+" not found\n").c_str());
//           exit(1);
//         }	
//         // Source of keys
//         source = cmph_io_nlfile_adapter(keys_fd);
//         cmph_config_t *config = cmph_config_new(source);
//         cmph_config_set_algo(config, CMPH_FCH);
//         hash = cmph_new(config);
//         // cmph_config_destroy(config);
//     }

//     unsigned int lookup(string key_str){ //Find key
//        const char *key = key_str.c_str();
//        unsigned int id = cmph_search(hash, key, (cmph_uint32)strlen(key));
//        return id;
//     }

//     // ~CMPH(){
//     //   //Destroy hash
//     //   cmph_destroy(hash);
//     //   cmph_io_nlfile_adapter_destroy(source);   
//     //   fclose(keys_fd);
//     // }
// };






class COLESS{
public:
	InputFile dedup_bitmatrix_file, spss_boundary_file, dup_bitmatrix_file, tmp_dir;
	LogFile logfile_main;
	LogFile debug1;
	LogFile debug2;
	long num_kmers;
	int M;
	int C;
	const int max_run = 16;
	vector<uint64_t> positions;
	HuffCodeMap huff_code_map;
	uint64_t CATEGORY_RUN=(uint64_t) 3;
	uint64_t CATEGORY_COLCLASS=(uint64_t) 0;
	uint64_t CATEGORY_COLVEC=(uint64_t) 2;
	int lm = 0;
	int lc = 0;

	int* per_simplitig_l;
	vector<char> spss_boundary; 

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
		debug1.init("debug1");
		debug2.init("debug2");

	}



	~COLESS(){
		mphf_destroy();
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
	void write_number_at_loc_advanced_by_block_sz(vector<uint64_t> & positions, uint64_t num, uint64_t loc_advanced_by_block_sz, uint64_t block_sz){ //requires loc_advanced_by_block_sz += block_size; 
		uint64_t j=0;
		uint64_t begin = loc_advanced_by_block_sz;

		while(num!=0)
		{
			if(num%2 == 1){
				positions.push_back(loc_advanced_by_block_sz-1-j); //b[loc_advanced_by_block_sz-1+j] = num%2;
			}
			j++;
			num /= 2;
			
		}

		debug1.fs<<-j<<" "<<block_sz<<endl;
		assert (-j<=block_sz);

	}

	void write_one(vector<uint64_t> & positions, uint64_t& b_it ){
		positions.push_back(b_it);
		b_it+=1;
	}

	void write_zero(vector<uint64_t> & positions, uint64_t& b_it ){
		b_it+=1;
	}
	
	void write_number_at_loc(vector<uint64_t> & positions, uint64_t num, uint64_t block_size, uint64_t& b_it ){
		write_number_at_loc_advanced_by_block_sz(positions, num, b_it+block_size, block_size);
		b_it += block_size; //successfully written and place on next bit; if size is 2, (0,1) written, now val is 2.
	}

	void write_unary_one_at_loc(vector<uint64_t> & positions, uint64_t unary_num, uint64_t& b_it ){
		for(uint64_t i = 0; i<unary_num; i++ ){
			positions.push_back(b_it);
			b_it+=1;
		}
	}
	void write_unary_zero_at_loc(vector<uint64_t> & positions, uint64_t unary_num, uint64_t& b_it ){
		b_it+=unary_num;
	}

	void write_binary_string_at_loc(vector<uint64_t> & positions, string binarystring, uint64_t& b_it){
		for (size_t i = 0; i< binarystring.length(); i++) {
			if (binarystring[i]=='1'){
				positions.push_back(b_it+i);
				
			}
		}
		b_it += binarystring.length();
	}

	void write_binary_vector_at_loc(vector<uint64_t> & positions, vector<bool> binary_vector, uint64_t& b_it){
		for (size_t i = 0; i< binary_vector.size(); i++) {
			if (binary_vector[i]== 1){
				positions.push_back(b_it+i);
			}
		}
		b_it += binary_vector.size();
	}

	bit_vector store_as_sdsl(vector<uint64_t>& positions, uint64_t bv_size, string filename){
		
		//bit_vector bv = bit_vector(bv_size, 0);
		bit_vector bv(bv_size, 0);
		for (uint64_t p: positions){
			bv[p] = 1;
		}
		if(filename=="rrr_main"){
			debug2.fs<<bv;
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

		LogFile log_num_color_in_class;
		log_num_color_in_class.init("log_num_color_in_class"); 
		for(int x=0; x<M; x++){
			string bv_line;
			getline(dedup_bitmatrix_file.fs, bv_line);
			unsigned int idx = lookup(bv_line);		// returns an if in range (0 to M-1) 
			assert(idx < M);

			array_hi[idx] = std::stoull(bv_line.substr(0,std::min(64,int(C))), nullptr, 2) ; 
			cout<<array_hi[idx]<<endl;
			write_number_at_loc(positions, array_hi[idx], min(64, C), b_it ); //array_hi[x] higher uint64_t

			array_lo[idx]=0;
			if(C > 64){
				string ss=bv_line.substr(64,(C-64));
				array_lo[idx]=std::stoull(ss, nullptr, 2);
				write_number_at_loc(positions, array_lo[idx], C-64, b_it ); //array_lo[x] lower uint64_t
			}

			int num_ones_in_color = __builtin_popcountll(array_hi[idx]) + __builtin_popcountll(array_lo[idx]) ;
			log_num_color_in_class.fs << num_ones_in_color <<endl;
		}
		cout << "Expected_MB_bv_mapping="<<(C*M)/8.0/1024.0/1024.0 << endl;
		dedup_bitmatrix_file.fs.close();

 		
		store_as_binarystring(positions, b_it, "bb_map" );
		store_as_sdsl(positions, b_it, "rrr_map" );

		cout << "expected_MB_bv_mapping="<<(C*M)/8.0/1024.0/1024.0 << endl;
		cout << "rrr_MB_bv_mapping="<<size_in_bytes(store_as_sdsl(positions, b_it, "rrr_bv_mapping.sdsl" ))/1024.0/1024.0 << endl;
		delete array_hi;
		delete array_lo;
	}


	void method1_pass0(bool skip_pass = false){ //load_huffman_table();//int M -> variable string   //writehuffman[in
		// void get_freq_count(){ // scan through all the color vectors to get freq count, global and local
		// // table of size M 
		// }

		double t_begin,t_end; struct timeval timet;
		printf("Construct a MPHF with  %lli elements  \n",M);
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
		create_table(dedup_bitmatrix_file.filename, M );
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
		double elapsed = t_end - t_begin;
		printf("CMPH constructed perfect hash for %llu keys in %.2fs\n", M,elapsed);

		//if(!skip_pass){
			gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
			OutputFile cmp_keys("cmp_keys");  // get frequency count
			for (uint64_t i=0; i < num_kmers; i+=1){
				string bv_line;
				getline (dup_bitmatrix_file.fs,bv_line);
				cmp_keys.fs<<lookup(bv_line)<<endl;
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

	void method1_pass1(bool skip_pass = false){ //decide whether to use local hash table, can skip
		vector<uint64_t> positions_local_table;
		uint64_t b_it_local_table = 0;
		
		dup_bitmatrix_file.rewind();
		store_global_color_class_table();
		// bit vector values
		uint64_t b_it=0;
		vector<uint64_t> positions; // positions for main vector

		for (uint64_t i=0; i < num_kmers; i+=1){ //load spss_boundary vector in memory from disk
			string spss_line;
			getline (spss_boundary_file.fs,spss_line); 
			spss_boundary.push_back(spss_line[0]); //this kmer starts a simplitig
		}
		per_simplitig_l = new int[spss_boundary.size()];

		//per simplitig values
		set<uint32_t> local_col_classes_uniq; //get the bool HuffCodeMap[M-1]
		set<uint32_t> local_col_classes_uniq_nonrun; //get the bool HuffCodeMap[M-1]

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

		OutputFile all_ls("all_ls");
		InputFile cmp_keys("cmp_keys");


		int simplitig_it = 0;
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
			
			//unsigned int curr_kmer_cc_id = lookup(bv_line); //uint64_t num = bphf->lookup(curr_bv);
			string curr_kmer_cc_id_str;
			getline(cmp_keys.fs, curr_kmer_cc_id_str);
			unsigned int curr_kmer_cc_id = std::stoull(curr_kmer_cc_id_str, nullptr, 10); 
			
				
			if(spss_boundary[i]=='0'){ // non-start
				int hd_hi = hammingDistance(prev_bv_hi, curr_bv_hi);
				int hd_lo = hammingDistance(prev_bv_lo, curr_bv_lo);
				int hd= hd_hi+hd_lo;
				if(hd==0){	//CAT=RUN
					skip+=1;	
					case_run+=1;	
				}else{ //CAT=NRUN
					case_nonrun += 1;
					if(skip!=0){ 	//not skipped, write lm
						//write_number_at_loc(positions, CATEGORY_RUN, 2, b_it);
						//write_number_at_loc(positions, skip, 1+floor(log2(skip)), b_it);
					}
					skip=0;

					sum_length_huff_nonrun += huff_code_map[curr_kmer_cc_id].size();
					local_col_classes_uniq_nonrun.insert(curr_kmer_cc_id);
					if(hd*lc < lm){ //CAT=LC
						case_dlc += 1;
						
						//cout<<"hd: "<<hd<<endl;
					}else{ //CAT=LM
						case_lm += 1;
						sum_length_huff += huff_code_map[curr_kmer_cc_id].size();
						local_col_classes_uniq.insert(curr_kmer_cc_id);
					}
					
				}
			}else{	//start of simplitig, so CAT=LM
				case_lm+=1;
				case_nonrun +=1;
				// if(skip!=0){ 	//not skipped, write lm
				// 	write_number_at_loc(positions, skip, 1+floor(log2(skip)), b_it);
				// }
				skip=0;
				
				sum_length_huff += huff_code_map[curr_kmer_cc_id].size();
				sum_length_huff_nonrun += huff_code_map[curr_kmer_cc_id].size();

				local_col_classes_uniq.insert(curr_kmer_cc_id);
				local_col_classes_uniq_nonrun.insert(curr_kmer_cc_id);

				// write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
			}

			if(spss_boundary[(i+1)%num_kmers]=='1'){	// end k-mer of simplitig
				// decide what to do
				int l = local_col_classes_uniq.size(); //case_lm

				per_simplitig_l[simplitig_it] = l;

				int ll = ceil(log2(l)*1.0);
				int l_nrun = local_col_classes_uniq_nonrun.size(); //case_lm
				int ll_nrun = ceil(log2(l_nrun)*1.0);


				case_nonrun = case_dlc + case_lm;
				assert(sum_length_huff_uniq==0);
				assert(sum_length_huff_uniq_nonrun==0);

				for(uint32_t uniq_col_class_id: local_col_classes_uniq){
					sum_length_huff_uniq += huff_code_map[uniq_col_class_id].size();
				}

				for(uint32_t uniq_col_class_id: local_col_classes_uniq_nonrun){
					sum_length_huff_uniq_nonrun += huff_code_map[uniq_col_class_id].size();
				}

				
				write_number_at_loc(positions_local_table, 1, 1, b_it_local_table); //num, bsize
				write_number_at_loc(positions_local_table, l, lm, b_it_local_table);
			


				for(uint32_t uniq_col_class_id: local_col_classes_uniq){
					//write_binary_vector_at_loc(positions_local_table, huff_code_map[uniq_col_class_id], b_it_local_table);
				}

				



				all_ls.fs << l <<" "<<ll<<" "<<l_nrun<<" "<<ll_nrun<<endl;  
				all_ls.fs << "P "<<ll*case_lm<<" "<<lm*case_lm<<" "<<sum_length_huff<<" "<<lm + sum_length_huff_uniq<<" "<<lm * (1+l)<<endl;
				use_local_hash = ( (ll - lm ) * case_lm + lm * (1+l)  ) ;  //ll*case_lm + (lm + l*lm) ::: lm * case_lm 
				use_local_hash_nonrun = ( (ll_nrun - lm ) * case_nonrun + lm * (1+l_nrun) ) ;  //ll*case_lm + (lm + l*lm) ::: lm * case_lm 
				use_local_hash_huff =  (ll*case_lm - sum_length_huff + lm + sum_length_huff_uniq) ;
				use_local_hash_huff_nonrun = ( ll_nrun*case_nonrun - sum_length_huff_nonrun + lm + sum_length_huff_uniq_nonrun  );

				//logfile_main.fs<<use_local_hash<<" "<<use_local_hash_nonrun<<" "<<use_local_hash_huff<<" "<<use_local_hash_huff_nonrun<<" "<<num_kmer_in_simplitig<<endl;
				logfile_main.fs<<ll*case_lm<<" "<<lm*case_lm<<" "<<sum_length_huff<<" "<<ll_nrun*case_lm<<" "<<lm*case_nonrun<<" "<<sum_length_huff_nonrun<<" "<<use_local_hash<<" "<<use_local_hash_nonrun<<" "<<use_local_hash_huff<<" "<<use_local_hash_huff_nonrun<<" s "<<case_run<<" c "<<case_dlc<<" m "<<case_lm<<" "<<num_kmer_in_simplitig<<" "<<sum_length_huff_uniq<<" "<<sum_length_huff_uniq_nonrun<<endl;
				
				//re-init for new simplitig
				//vector<uint32_t>().swap(local_col_classes_uniq);//
				local_col_classes_uniq.clear();
				local_col_classes_uniq_nonrun.clear();
				num_kmer_in_simplitig = 0;

				
				use_local_hash = use_local_hash_nonrun = use_local_hash_huff = use_local_hash_huff_nonrun = 0;
				skip=0;
				case_run = case_lm = case_nonrun = case_dlc = 0;
				sum_length_huff = sum_length_huff_nonrun = sum_length_huff_uniq = sum_length_huff_uniq_nonrun =  0;
			
				simplitig_it+=1;
			}
			prev_bv_hi = curr_bv_hi;
			prev_bv_lo = curr_bv_lo;
		}
		store_as_binarystring(positions_local_table, b_it_local_table, "bb_local_table");
		store_as_sdsl(positions_local_table, b_it_local_table, "rrr_local_table");
		
		
	}

	void method1_pass2(){
		vector<uint64_t> positions;
		uint64_t b_it = 0;
		dup_bitmatrix_file.rewind();

		uint64_t curr_bv_hi =  0;
		uint64_t curr_bv_lo = 0;
		uint64_t prev_bv_hi = 0;
		uint64_t prev_bv_lo = 0;
		uint64_t skip = 0;

		InputFile cmp_keys("cmp_keys");
		int simplitig_it = 0;
		int l = per_simplitig_l[0];
		int ll = ceil(log2(l));
		int lm_or_ll = ll;
		Hashtable local_ht;
		for (uint64_t i=0; i < num_kmers; i+=1){ 
			
			//load the color vector of current k-mer from disk to "curr_bv_hi/lo"
			string bv_line;
			getline (dup_bitmatrix_file.fs,bv_line); // bv line = color vector C bits
			curr_bv_hi = std::stoull(bv_line.substr(0,std::min(64, C)), nullptr, 2);
			curr_bv_lo = 0;
			if(C > 63){
				curr_bv_lo = std::stoull(bv_line.substr(64,bv_line.length()-64), nullptr, 2);
			} 

			//unsigned int curr_kmer_cc_id = lookup(bv_line); //uint64_t num = bphf->lookup(curr_bv);
			string curr_kmer_cc_id_str;
			getline(cmp_keys.fs, curr_kmer_cc_id_str);
			unsigned int curr_kmer_cc_id = std::stoull(curr_kmer_cc_id_str, nullptr, 10); 
			
			if(spss_boundary[i]=='0'){ // non-start
				int hd_hi = hammingDistance(prev_bv_hi, curr_bv_hi);
				int hd_lo = hammingDistance(prev_bv_lo, curr_bv_lo);
				int hd= hd_hi+hd_lo;
				int lmaxrun = ceil(log2(max_run));
				
				if(hd==0){	//CATEGORY=RUN
					skip+=1;	
					//case_run+=1;	
				}else{ //CATEGORY=NOT_RUN
					//case_nonrun += 1;
					if(skip!=0){ 	//not skipped, run break, write lm
						// write_number_at_loc(positions, CATEGORY_RUN, (uint64_t) 2, b_it);
						// write_unary_one_at_loc(positions, (uint64_t) skip, b_it);

						int q = floor(skip/max_run);
						int rem = skip % max_run;
						assert(skip == q*max_run + rem); //skip = q*max_run + rem
						// write_number_at_loc(positions, CATEGORY_RUN, (uint64_t) 2, b_it);
						// write_unary_zero_at_loc(positions, (uint64_t) q, b_it);
						// write_one(positions, b_it);
						// write_number_at_loc(positions, (uint64_t) rem, (uint64_t) lmaxrun, b_it);
					}
					skip=0;

					//if(hd*(lc + 1) < huff_code_map[curr_kmer_cc_id].size()){ //CATEGORY=LC
					if(hd*(lc + 1) < lm){ //CATEGORY=LC

						//case_dlc += 1;
						//write_number_at_loc(positions, CATEGORY_COLVEC, 2, b_it);
						for (int i_bit=0; i_bit < lc; i_bit+=1){
							if ((( prev_bv_hi >>  i_bit) & 1) != (( curr_bv_hi >>  i_bit) & 1)){ 
								//write_number_at_loc(positions, i_bit, lc, b_it); // i_bit is the different bit loc
							}
						}
						for (int i_bit=0; i_bit < lc; i_bit+=1){
							if ((( prev_bv_lo >>  i_bit) & 1) != (( curr_bv_lo >>  i_bit) & 1)){
								//write_number_at_loc(positions, i_bit, lc, b_it);	//i_bit is the different bit loc
							}
						}
					}else{ //CATEGORY=LM
						//case_lm += 1;
						write_number_at_loc(positions, CATEGORY_COLCLASS, 1, b_it);
						write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
						assert(curr_kmer_cc_id<M && curr_kmer_cc_id>0);

						//write_number_at_loc(positions, local_ht.put_and_getid(curr_kmer_cc_id), ll, b_it);
						//write_binary_vector_at_loc(positions, huff_code_map[curr_kmer_cc_id], b_it);
					}	
				}
			}else{	//start of simplitig, so CAT=LM
				l = per_simplitig_l[simplitig_it];
				ll = ceil(log2(l));
				lm_or_ll = ll;

				//case_lm+=1;
				//case_nonrun +=1;
				
				write_number_at_loc(positions, CATEGORY_COLCLASS, 1, b_it);
				//write_number_at_loc(positions, local_ht.put_and_getid(curr_kmer_cc_id), ll, b_it);
				write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
				assert(curr_kmer_cc_id<M && curr_kmer_cc_id>0);
				//write_binary_vector_at_loc(positions, huff_code_map[curr_kmer_cc_id], b_it);



				skip=0;
			}

			if(spss_boundary[(i+1)%num_kmers]=='1'){	// end k-mer of simplitig
				local_ht.clear();
				skip=0;
				simplitig_it+=1;
			}

			prev_bv_hi=curr_bv_hi;
			prev_bv_lo=curr_bv_lo;

		}

		OutputFile positions_out("positions_out");
		for(uint64_t tt: positions){
			positions_out.fs<<tt<<endl;
		}
		cout<<"b_it_size: "<<b_it<<endl;
		store_as_sdsl(positions, b_it, "rrr_main");
		store_as_binarystring(positions, b_it, "bb_main");
	}
};

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
	
	coless.method1_pass0();
	coless.method1_pass1();
	// coless.method1_pass2();

	//COLESS_Decompress cdec(num_kmers, M, C);


	return EXIT_SUCCESS;
}

