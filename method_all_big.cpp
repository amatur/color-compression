//version: jan 30: not working, pause
#include<cmph.h> //#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <vector>
#include <stack>

#include <deque>

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

const bool DEBUG_MODE = false;

namespace TimeMeasure
{
	double t_begin,t_end; struct timeval timet;
	void time_start(){
		gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	}
	void time_end(string msg){
		gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
		cout<<msg<<" time = ";
		printf("%.2fs\n",t_end - t_begin);
	}
} using namespace TimeMeasure;



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
	void close(){
		fs.close();
	}
	~OutputFile(){
		fs.close();
	}
};
class DebugFile : public OutputFile	//derived class
{
	public:
		DebugFile(string filename){
			if(!DEBUG_MODE){
					this->filename = filename;
					fs.open (filename.c_str(),  std::fstream::out );
			}

		}
		DebugFile(){

		}
		void init(const std::string filename)
		{
			if(!DEBUG_MODE){
				this->filename = filename;
				this->fs.open(this->filename, fstream::out);
			}
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

	void close(){
		fs.close();
	}
	~InputFile(){
		fs.close();
	}
};



class Hashtable {

	public:
    std::unordered_map<uint32_t, uint32_t> htmap; // m_to_l
	uint32_t curr_id = 0;
	Hashtable(){
		curr_id = 0;
	}

    uint64_t put_and_getid(uint32_t key) {
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

	bool exists(uint32_t key){
		return htmap.count(key) > 0;
	}

	void clear(){
		htmap.clear();
		curr_id = 0;
	}

	vector<uint32_t> get_array(){
		vector<uint32_t> array(curr_id, 0);
		for (auto x : htmap){
			array[x.second] =  x.first  ;
			//cout<<x.first<<"->"<<x.second<<endl;
		}
		return array;
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
        //double t_begin,t_end; struct timeval timet;
        printf("Construct a BooPHF with  %lli elements  \n",nelem);
        //gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
        auto data_iterator = boomphf::range(static_cast<const uint64_t*>(data), static_cast<const uint64_t*>(data+nelem));
        double gammaFactor = 7.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
        bphf = new boomphf::mphf<uint64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
        
		//gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
        //printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,t_end - t_begin);
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

	

	// u_int32_t HuffDecode(const INode* root, string s, int& loc)
	// {
	// 	string ans = "";
	// 	u_int32_t ansint=0;
	// 	const INode* curr = root;
	// 	for (int i = 0; i < s.size(); i++) {
	// 		if(  const LeafNode* lf  = dynamic_cast<const LeafNode*>(curr) ){
	// 				ansint = (u_int32_t)(lf->c);
	// 				return ansint;
	// 				curr = root;
	// 		}else if(const InternalNode* internal   = dynamic_cast<const InternalNode*>(curr) ){
	// 			if (s[i] == '0')
	// 				curr = internal->left;
	// 			else
	// 				curr = internal->right;

	// 			if(const LeafNode* lf  = dynamic_cast<const LeafNode*>(curr)){
	// 				ansint = (u_int32_t)(lf->c);
	// 				cout<< ansint;
	// 			}
	// 		}else{
	// 			exit(2);
	// 		}
	// 	}
	// 	// cout<<ans<<endl;
	// 	return ansint;
	// }
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
	InputFile dedup_bitmatrix_file, spss_boundary_file, dup_bitmatrix_file;
	
	
	DebugFile logfile_main;
	DebugFile debug1;
	DebugFile debug2;
	DebugFile all_ls;
	long num_kmers;
	int M;
	int C;
	int max_run = 16;
	vector<uint64_t> positions;
	HuffCodeMap huff_code_map;
	uint64_t CATEGORY_RUN=(uint64_t) 3;
	uint64_t CATEGORY_COLCLASS=(uint64_t) 0;
	uint64_t CATEGORY_COLVEC=(uint64_t) 2;
	int lm = 0;
	int lc = 0;
	string* global_table;
	int* per_simplitig_l;
	bool* per_simplitig_use_local;
	
	//per_simplitig_d
	//per simplitig use_local_hash
	vector<char> spss_boundary; 

	//run param
	int d_class_diff = 1; //0,1,2

	bool USE_LOCAL_TABLE = true;
    bool USE_HUFFMAN = true;
	bool ALWAYS_LOCAL_OR_GLOBAL = false;

	COLESS(long num_kmers, int M, int C, string dedup_bitmatrix_fname, string dup_bitmatrix_fname, string spss_boundary_fname, int max_run){
		dedup_bitmatrix_file.init(dedup_bitmatrix_fname);
		spss_boundary_file.init(spss_boundary_fname);
		dup_bitmatrix_file.init(dup_bitmatrix_fname);
		this->max_run = max_run;
		this->num_kmers = num_kmers;
		this->M = M;
		this->C = C;
		this->lm = ceil(log2(M));
		this->lc = ceil(log2(C));
		logfile_main.init("log_coless");
		debug1.init("debug1");
		debug2.init("debug2");

		all_ls.init("all_ls");

	}



	~COLESS(){
		mphf_destroy();
		delete per_simplitig_l;
		delete per_simplitig_use_local;
	}

	// template <typename T> void dump_to_disk(T& vec, uint64_t last_written_pos, fstream fs)
	// {
		
	// 	fs << 
	// }

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
		stack<uint64_t> qpositions;
		if(num!=0){
			if(block_sz==0)
				cout<<"must be non zero block size"<<endl;
		}
		while(num!=0)
		{
			if(num%2 == 1){
				//positions.push_back(loc_advanced_by_block_sz-1-j); //b[loc_advanced_by_block_sz-1+j] = num%2;
				qpositions.push(loc_advanced_by_block_sz-1-j);
			}
			j++;
			num /= 2;
			
		}
		while(!qpositions.empty()){
			positions.push_back(qpositions.top());
			qpositions.pop();
		}

		if(DEBUG_MODE) debug1.fs<<-j<<" "<<block_sz<<endl;
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
		uint64_t lastp  =0;
		for (uint64_t p: positions){
			bv[p] = 1;
		}
		if(filename=="rrr_main"){
			if(DEBUG_MODE) debug2.fs<<bv;
		}
		rrr_vector<256> rrr_bv(bv);
		//cout << "rrr_MB_bv_mapping="<<size_in_bytes(rrr_bv_mapping)/1024.0/1024.0 << endl;
		store_to_file(rrr_bv, filename);	//"rrr_bv_mapping.sdsl"

		return bv;
	}

	void store_as_binarystring(vector<uint64_t>& positions, uint64_t bv_size, string filename){
		OutputFile binarystring_file(filename);
		//sort(positions.begin(), positions.end());
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

		//LogFile log_num_color_in_class;
		//log_num_color_in_class.init("log_num_color_in_class"); 
		
		global_table = new string[M];
		for(int x=0; x<M; x++){
			string bv_line;
			getline(dedup_bitmatrix_file.fs, bv_line);
			unsigned int idx = lookup(bv_line);		// returns an if in range (0 to M-1) 
			assert(idx < M);
			global_table[idx] = bv_line;
			assert(x==idx);

			array_lo[idx] = std::stoull(bv_line.substr(0,std::min(64,int(C))), nullptr, 2) ; 
			write_number_at_loc(positions, array_lo[idx], min(64, C), b_it ); //array_hi[x] higher uint64_t
			array_hi[idx]=0;
			if(C > 64){
				string ss=bv_line.substr(64,(C-64));
				array_hi[idx]=std::stoull(ss, nullptr, 2);
				write_number_at_loc(positions, array_hi[idx], C-64, b_it ); //array_lo[x] lower uint64_t
			}

			// if(DEBUG_MODE) {
			// 	int num_ones_in_color = __builtin_popcountll(array_hi[idx]) + __builtin_popcountll(array_lo[idx]) ;
			// 	log_num_color_in_class.fs << num_ones_in_color <<endl;
			// }

		}
		dedup_bitmatrix_file.fs.close();

		store_as_binarystring(positions, b_it, "bb_map" );
		store_as_sdsl(positions, b_it, "rrr_map" );

		cout << "expected_MB_bv_mapping="<<(C*M)/8.0/1024.0/1024.0 << endl;
		cout << "rrr_MB_bv_mapping="<<size_in_bytes(store_as_sdsl(positions, b_it, "rrr_bv_mapping.sdsl" ))/1024.0/1024.0 << endl;
		delete array_hi;
		delete array_lo;
	}


	void method1_pass0(){ //load_huffman_table();//int M -> variable string   //writehuffman[in
				
		time_start();
		create_table(dedup_bitmatrix_file.filename, M );
		time_end("CMPH constructed perfect hash for "+to_string(M)+" keys.");

		time_start();
		OutputFile cmp_keys("cmp_keys");  // get frequency count
		
		for (uint64_t i=0; i < num_kmers; i+=1){
			string bv_line;
			getline (dup_bitmatrix_file.fs,bv_line);
			cmp_keys.fs<<lookup(bv_line)<<endl;
		}
		time_end("CMPH lookup for "+to_string(num_kmers)+"keys.");
		cmp_keys.close();

		time_start();
		system("cat cmp_keys | sort -n | uniq -c | rev | cut -f 2 -d\" \" | rev > frequency_sorted");
		time_end("Sorting and getting freq for "+to_string(num_kmers)+" keys.");
		
		time_start();
		InputFile infile_freq("frequency_sorted");
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
		time_end("Read freq for "+to_string(M)+" values.");
		infile_freq.close();

		time_start();
		INode* root = BuildTree(frequencies, M);
        GenerateCodes(root, HuffCode(), huff_code_map); // huff_code_map is filled: uint32t colclassid-> vector bool
		delete frequencies;
		delete root;
		time_end("Build huffman tree on " +to_string(M)+" values.");

	}

	void method1_pass1(bool skip_pass = false){ //decide whether to use local hash table, can skip
		dup_bitmatrix_file.rewind();
		time_start();
		store_global_color_class_table();
		time_end("Written global table for "+to_string(M)+" values.");
		
		uint64_t b_it=0;
		vector<uint64_t> positions; // positions for main vector

		for (uint64_t i=0; i < num_kmers; i+=1){ //load spss_boundary vector in memory from disk
			string spss_line;
			getline (spss_boundary_file.fs,spss_line); 
			spss_boundary.push_back(spss_line[0]); //this kmer starts a simplitig
		}
		per_simplitig_l = new int[spss_boundary.size()];
		per_simplitig_use_local =  new bool[spss_boundary.size()];

		//per simplitig values		
		Hashtable local_hash_table;
		int use_local_hash_nonrun = 0;
		int use_local_hash_huff_nonrun = 0;
		uint64_t sum_length_huff_nonrun = 0;
		uint64_t sum_length_huff_uniq_nonrun = 0;
		uint64_t num_kmer_in_simplitig = 0;
		//
		uint64_t skip=0;
		int case_run = 0;
		int case_lm = 0;
		int case_nonrun = 0;
		int case_dlc = 0;
		//		
		vector<uint64_t> positions_local_table;
		uint64_t b_it_local_table = 0;

		//per kmer values
		uint64_t curr_bv_hi = 0;
		uint64_t curr_bv_lo = 0;
		uint64_t prev_bv_hi = 0;
		uint64_t prev_bv_lo = 0;
	
		//InputFile cmp_keys("cmp_keys");
		int simplitig_it = 0;
		for (uint64_t i=0; i < num_kmers; i+=1){ 
			string bv_line;
			getline (dup_bitmatrix_file.fs,bv_line); // bv line = color vector C bits
			curr_bv_lo = std::stoull(bv_line.substr(0,std::min(64, C)), nullptr, 2);
			curr_bv_hi = 0;
			if(C >= 64){
				curr_bv_hi = std::stoull(bv_line.substr(64,bv_line.length()-64), nullptr, 2);
			} 

			//per kmer task
			num_kmer_in_simplitig+=1;  //start of simplitig id: num_kmer_in_simplitig
			
			unsigned int curr_kmer_cc_id = lookup(bv_line); //uint64_t num = bphf->lookup(curr_bv);
			// string curr_kmer_cc_id_str;
			// getline(cmp_keys.fs, curr_kmer_cc_id_str);
			// unsigned int curr_kmer_cc_id = std::stoull(curr_kmer_cc_id_str, nullptr, 10); 
			
			if(spss_boundary[i]=='0'){ // non-start
				int hd_hi = hammingDistance(prev_bv_hi, curr_bv_hi);
				int hd_lo = hammingDistance(prev_bv_lo, curr_bv_lo);
				int hd= hd_hi+hd_lo;
				if(hd==0){	//CAT=RUN
					skip+=1;	
					case_run+=1;	
				}else{ //CAT=NRUN
					case_nonrun += 1;
					skip=0;
					local_hash_table.put_and_getid(curr_kmer_cc_id);
					if(hd*(1+lc) < lm && hd == 1){ //CAT=LC
						case_dlc += 1;
					}else{ //CAT=LM
						case_lm += 1;
					}
				}
			}else{	//start of simplitig, so CAT=LM
				case_lm+=1;
				case_nonrun +=1;
				skip=0;
				sum_length_huff_nonrun += huff_code_map[curr_kmer_cc_id].size();
				local_hash_table.put_and_getid(curr_kmer_cc_id);
			}

			if(spss_boundary[(i+1)%num_kmers]=='1'){	// end k-mer of simplitig
				int l = local_hash_table.curr_id; 
				int ll = ceil(log2(l)*1.0);
				per_simplitig_l[simplitig_it] = l;

				case_nonrun = case_dlc + case_lm;
				
				if(USE_LOCAL_TABLE){
					if(!ALWAYS_LOCAL_OR_GLOBAL){
						write_number_at_loc(positions_local_table, 1, 1, b_it_local_table); //if always use local table, skip
					}
					
					write_number_at_loc(positions_local_table, l, lm, b_it_local_table);

					vector<uint32_t> local_ht_arr = local_hash_table.get_array();
					for(uint32_t i = 0 ; i< local_hash_table.curr_id; i++){
						uint32_t uniq_col_class_id = local_ht_arr[i];
						sum_length_huff_uniq_nonrun += huff_code_map[uniq_col_class_id].size();
						write_binary_vector_at_loc(positions_local_table, huff_code_map[uniq_col_class_id], b_it_local_table);
					}
					local_ht_arr.clear();
				}

				if(DEBUG_MODE) all_ls.fs << l <<" "<<ll<<endl;  
				use_local_hash_nonrun = ( (ll - lm ) * case_nonrun + lm * (1+l) ) ;  //ll*case_lm + (lm + l*lm) ::: lm * case_lm 
				use_local_hash_huff_nonrun = ( ll*case_nonrun - sum_length_huff_nonrun + lm + sum_length_huff_uniq_nonrun  );

				if(use_local_hash_huff_nonrun < 0){
					per_simplitig_use_local[simplitig_it] = true;
				}else{
					per_simplitig_use_local[simplitig_it] = false;
				}



				//logfile_main.fs<<use_local_hash<<" "<<use_local_hash_nonrun<<" "<<use_local_hash_huff<<" "<<use_local_hash_huff_nonrun<<" "<<num_kmer_in_simplitig<<endl;
				if(DEBUG_MODE) logfile_main.fs<<num_kmer_in_simplitig<<" "<< l << " " << l*lm <<" "<<sum_length_huff_uniq_nonrun<<endl;
				
				//re-init for new simplitig
				local_hash_table.clear();
				num_kmer_in_simplitig = 0;

				use_local_hash_nonrun =  use_local_hash_huff_nonrun = 0;
				skip=0;
				case_run = case_lm = case_nonrun = case_dlc = 0;
				sum_length_huff_nonrun = sum_length_huff_uniq_nonrun =  0;
			
				simplitig_it+=1;
			}
			prev_bv_hi = curr_bv_hi;
			prev_bv_lo = curr_bv_lo;
		}
		store_as_binarystring(positions_local_table, b_it_local_table, "bb_local_table");
		store_as_sdsl(positions_local_table, b_it_local_table, "rrr_local_table");
	}

	void method1_pass2()
	{
		vector<uint64_t> positions;
		uint64_t b_it = 0;
		dup_bitmatrix_file.rewind();
		DebugFile cases_smc("cases_smc");
		uint64_t curr_bv_hi = 0;
		uint64_t curr_bv_lo = 0;
		uint64_t prev_bv_hi = 0;
		uint64_t prev_bv_lo = 0;
		uint64_t skip = 0;

		//InputFile cmp_keys("cmp_keys");
		int simplitig_it = 0;
		int l = per_simplitig_l[0];
		int ll = ceil(log2(l));
		int lm_or_ll;
		if (USE_LOCAL_TABLE)
		{
			lm_or_ll = ll;
		}
		else
		{
			lm_or_ll = lm;
		}
		Hashtable local_ht;
		int lmaxrun = ceil(log2(max_run));
		for (uint64_t i = 0; i < num_kmers; i += 1)
		{
			l = per_simplitig_l[simplitig_it];
			ll = ceil(log2(l));
			if(DEBUG_MODE) all_ls.fs << l << endl;

			// load the color vector of current k-mer from disk to "curr_bv_hi/lo"
			string bv_line;
			getline(dup_bitmatrix_file.fs, bv_line); // bv line = color vector C bits

			curr_bv_lo = std::stoull(bv_line.substr(0, std::min(64, C)), nullptr, 2);
			curr_bv_hi = 0;
			if (C > 64 - 1)
			{
				curr_bv_hi = std::stoull(bv_line.substr(64, bv_line.length() - 64), nullptr, 2);
			}

			unsigned int curr_kmer_cc_id = lookup(bv_line); //uint64_t num = bphf->lookup(curr_bv);
			// string curr_kmer_cc_id_str;
			// getline(cmp_keys.fs, curr_kmer_cc_id_str);
			// unsigned int curr_kmer_cc_id = std::stoull(curr_kmer_cc_id_str, nullptr, 10);

			if (spss_boundary[i] == '0')
			{ // non-start
				int hd_hi = hammingDistance(prev_bv_hi, curr_bv_hi);
				int hd_lo = hammingDistance(prev_bv_lo, curr_bv_lo);
				int hd = hd_hi + hd_lo;

				if (hd == 0)
				{ // CATEGORY=RUN
					skip += 1;
					// case_run+=1;
					if(DEBUG_MODE) cases_smc.fs << "r" << endl;
				}
				else
				{ // CATEGORY=NOT_RUN
					// case_nonrun += 1;
					if (skip != 0)
					{ // not skipped, run break, write lm
						// paul method
						{
							int q, rem;
							q = floor(skip / max_run);
							rem = skip % max_run;
							assert(skip == q * max_run + rem); // skip = q*max_run + rem
							write_number_at_loc(positions, CATEGORY_RUN, (uint64_t)2, b_it);
							write_unary_one_at_loc(positions, (uint64_t)q, b_it);
							write_zero(positions, b_it);
							write_number_at_loc(positions, (uint64_t)rem, (uint64_t)lmaxrun, b_it);
						}

						// my method
						{
							// write_number_at_loc(positions, CATEGORY_RUN, (uint64_t) 2, b_it);
							// write_unary_one_at_loc(positions, (uint64_t) ceil(log2(skip)), b_it);
							// write_zero(positions, b_it);
							// write_number_at_loc(positions, (uint64_t) skip, (uint64_t) ceil(log2(skip)), b_it);
						}
					}
					skip = 0;

					if (hd * (lc + 1) < lm_or_ll && hd == 1)
					{ // CATEGORY=LC
						// if(hd*(lc + 1) < huff_code_map[curr_kmer_cc_id].size() && hd==1 ){ //CATEGORY=LC
						// if(hd*(lc + 1) < lm && hd==1){ //CATEGORY=LC
						if(DEBUG_MODE) cases_smc.fs << "d" << endl;

						// case_dlc += 1;

						for (int i_bit = 0; i_bit < 64 && i_bit < C; i_bit += 1)
						{
							if (((prev_bv_lo >> i_bit) & 1) != ((curr_bv_lo >> i_bit) & 1))
							{
								write_number_at_loc(positions, CATEGORY_COLVEC, 2, b_it);
								write_number_at_loc(positions, i_bit, lc, b_it); // i_bit is the different bit loc
							}
						}
						for (int i_bit = 64; i_bit < C; i_bit += 1)
						{
							int actual_i_bit = i_bit - 64;
							if (((prev_bv_hi >> actual_i_bit) & 1) != ((curr_bv_hi >> actual_i_bit) & 1))
							{
								write_number_at_loc(positions, CATEGORY_COLVEC, 2, b_it);
								write_number_at_loc(positions, i_bit, lc, b_it); // i_bit is the different bit loc
							}
						}
					}
					else
					{ // CATEGORY=LM
						if(DEBUG_MODE) cases_smc.fs << "l" << endl;

						// case_lm += 1;

						write_number_at_loc(positions, CATEGORY_COLCLASS, 1, b_it);
						if (USE_LOCAL_TABLE)
						{
							uint64_t localid = local_ht.put_and_getid(curr_kmer_cc_id);
							if (ll == 0 && localid == 1)
							{
								cout << "trouble" << endl;
							}
							write_number_at_loc(positions, localid, ll, b_it);
						}
						else
						{
							if (USE_HUFFMAN)
							{
								write_binary_vector_at_loc(positions, huff_code_map[curr_kmer_cc_id], b_it);
							}
							else
							{
								write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
							}
							// assert(curr_kmer_cc_id<M && curr_kmer_cc_id>0);
						}
					}
				}
			}
			else
			{ // start of simplitig, so CAT=LM
				if(DEBUG_MODE) cases_smc.fs << "l" << endl;

				l = per_simplitig_l[simplitig_it];
				ll = ceil(log2(l));
				lm_or_ll = ll;

				// case_lm+=1;
				// case_nonrun +=1;

				write_number_at_loc(positions, CATEGORY_COLCLASS, 1, b_it);
				if (USE_LOCAL_TABLE)
				{
					uint64_t localid = local_ht.put_and_getid(curr_kmer_cc_id);
					if (ll == 0)
					{
						assert(localid == 0);
					}
					write_number_at_loc(positions, localid, ll, b_it);
				}
				else
				{
					if (USE_HUFFMAN==true)
						write_binary_vector_at_loc(positions, huff_code_map[curr_kmer_cc_id], b_it);
					else
						write_number_at_loc(positions, curr_kmer_cc_id, lm, b_it);
					// assert(curr_kmer_cc_id<M && curr_kmer_cc_id>0);
				}
			}

			if (spss_boundary[(i + 1) % num_kmers] == '1')
			{ // end k-mer of simplitig
				local_ht.clear();
				simplitig_it += 1;
				if (USE_LOCAL_TABLE){

					lm_or_ll = ll;
				}
				if (skip != 0)
				{ // not skipped, run break, write lm
					int q, rem;
					q = floor(skip / max_run);
					rem = skip % max_run;
					assert(skip == q * max_run + rem); // skip = q*max_run + rem
					// paul method
					write_number_at_loc(positions, CATEGORY_RUN, (uint64_t)2, b_it);
					write_unary_one_at_loc(positions, (uint64_t)q, b_it);
					write_zero(positions, b_it);
					write_number_at_loc(positions, (uint64_t)rem, (uint64_t)lmaxrun, b_it);
					// my method //100001
				}
				skip = 0;
			}
			prev_bv_hi = curr_bv_hi;
			prev_bv_lo = curr_bv_lo;
		}

		DebugFile positions_out("positions_out");
		for (uint64_t tt : positions)
		{
			if(DEBUG_MODE) positions_out.fs << tt << endl;
		}
		cout << "b_it_size: " << b_it << endl;
		store_as_sdsl(positions, b_it, "rrr_main");
		store_as_binarystring(positions, b_it, "bb_main");
	}
};

int main (int argc, char* argv[]){
	vector<string> args(argv + 1, argv + argc);
    string dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname;
	//string tmp_dir;
    int M, C;
	int max_run = 16;
	long num_kmers=0;
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: tool -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers> -s <spss-bound> -x <max-run>" << endl;
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
		}else if (*i == "-x") {
            max_run = std::stoi(*++i);
		}
		// else if (*i == "-t") {
        //     tmp_dir  = *++i;
		// }
    }

	COLESS coless(num_kmers, M, C, dedup_bitmatrix_fname, dup_bitmatrix_fname, spss_boundary_fname, max_run);
	
	coless.method1_pass0();
	time_start();
	coless.method1_pass1();
	time_end("pass1.");

	time_start();
	coless.method1_pass2();
	time_end("pass2.");


	//COLESS_Decompress cdec(num_kmers, M, C);


	return EXIT_SUCCESS;
}
