#include "BooPHF.h"

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

#include<sstream>

using namespace std;
using namespace sdsl;

namespace BitManip
{
	uint64_t bitsToShort(char* bits, int unit=64) {
		int j = unit;
		uint64_t res = 0;
		for (int j = unit; j > 0; j--) {
			int i = unit-j;
			if (bits[j]=='1') {
				res |= 1 << i;
			}
		}
    	return res;
	}
	bool* shortToBits(short value) {
		bool* bits = new bool[64];
		int count = 0;
		while(value) {
			if (value&1)
				bits[count] = 1;
			else
				bits[count] = 0;
			value>>=1;
			count++;
		}
		return bits;
	}
}

// 	void load_mph_arr(string sdsl_file){
// 			//load bv
// 			bit_vector b ;
// 			size_t ones = rank_support_v<1>(&b)(b.size());
// 			bit_vector::select_1_type b_sel(&b);

// 			for (size_t i=1; i <= ones; ++i)
// 				int x = b_sel(i);

// 			}
// 			bool
// 	// uint8_t* uint64_to_arr(){
// 	// 	uint64_t number = 23425432542254234532;
// 	// 	uint8_t* result= new uint8_t;
// 	// 	for(int i = 0; i < 8; i++) {
// 	// 		result[i] = uint8_t((number >> 8*(7 - i)) & 0xFF);
// 	// 	}
// 	// 	return 
// 	// }
// }
// class Decompression
// {
// 	public:
// 	rrr_vector<256> rrr_bv;
// 	//rrr_vector<256> rrr_cc_map;
// 	stringstream ss_cc_map;

// 	Decompression(string index_file="rrr_bv_mapping.sdsl", string bit_file="rrrbv.sdsl", int M, int C){
// 		rrr_bv = rrr_vector<256>();
// 		load_from_file(rrr_bv, bit_file);
// 		cout<<rrr_bv<<endl;

// 		std::ofstream out("str_bv_mapping.txt");
// 		out << rrr_bv;
// 		out.close();
		
// 		rrr_vector<256> cc_map = rrr_vector<256>();
// 		load_from_file(cc_map, index_file);
// 		cout<<cc_map<<endl;
// 		ss_cc_map<<cc_map;
// 	}

// 	void getArray(){
// 		char ch;
// 		fstream fin("str_bv_mapping.txt", fstream::in);
// 		while (fin >> noskipws >> ch) {
// 			cout << ch; // Or whatever
// 			if(ch=='1'){
// 				char* bits = new char[lm];
// 				for (int c=0; c<lm; c++){
// 					fin >> noskipws >> bits[c];
// 				}
// 				uint64_t numberar = BitManip::bitsToShort(bits, lm);
// 				//outfile<<getColorVector(numberar)<<endl;
// 			}else{

// 			}
// 		}
// 	}
	
// };


class OutputFile{
  public:
  string filename;
  std::ofstream fs;

  OutputFile(string filename){
    this->filename = filename;
    fs.open (filename,  std::fstream::out );
  }
  void write(string towrite){
      fs << towrite;
      // <<endl;
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
};


typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;

int main (int argc, char* argv[]){
	vector<string> args(argv + 1, argv + argc);
    string in_bitmatrix_fname, dup_bitmatrix_fname;
    int M, num_colors;
	long num_kmers=0;
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: foomatic -i <DE-DUP-bitmatrix> -d <dup-bitmatrix> -c <num-colors> -m <M> -k <num-kmers>" << endl;
            return 0;
        } else if (*i == "-i") {
            in_bitmatrix_fname = *++i;
        } else if (*i == "-d") {
            dup_bitmatrix_fname = *++i;
        }else if (*i == "-c") {
            num_colors = std::stoi(*++i);
        }else if (*i == "-m") {
            M = std::stoi(*++i);
        }else if (*i == "-k") {
            num_kmers = std::stol(*++i);
        }
    }

    InputFile dedup_bitmatrix_file(in_bitmatrix_fname);
    string bv_line;

	InputFile dup_bitmatrix_file(dup_bitmatrix_fname);
    

	//PARAMETERS
	u_int64_t nelem = M;
	uint nthreads = 1;

	uint64_t ii, jj;
	u_int64_t *data;

	uint64_t rab = 0;
    data = (u_int64_t * ) calloc(nelem+rab,sizeof(u_int64_t));
	
	for (u_int64_t i = 0; i < nelem+rab; i++){
        getline (dedup_bitmatrix_file.fs,bv_line);
		data[i] = std::stoi(bv_line, nullptr, 2);
		
	}
	cout<<data[0]<<"->this0"<<endl;
	cout<<data[1]<<"->this1"<<endl;
	cout<<data[2]<<"->this2"<<endl;
	
	boophf_t * bphf = NULL;
	double t_begin,t_end; struct timeval timet;
	printf("Construct a BooPHF with  %lli elements  \n",nelem);
	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	
	// mphf takes as input a c++ range. A simple array of keys can be wrapped with boomphf::range
	// but could be from a user defined iterator (enabling keys to be read from a file or from some complex non-contiguous structure)
	auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(data), static_cast<const u_int64_t*>(data+nelem));
	
	double gammaFactor = 1.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query

	//build the mphf
	bphf = new boomphf::mphf<u_int64_t,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
	
	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
	double elapsed = t_end - t_begin;
	
	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,elapsed);
	printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);
	
	//query mphf like this
	uint64_t  idx = bphf->lookup(data[0]);
	printf(" example query  %lli ----->  %llu \n",data[0],idx);

	uint64_t *array;
	array=(uint64_t *) malloc(M*sizeof(uint64_t));

	for(int x=0; x<M; x++){
		array[bphf->lookup(data[x])] = data[x];
	}

	int block_size_here = num_colors;
	vector<int> ppositions;
	int pos_in_bv=0;
	bit_vector bv_mapping = bit_vector(num_colors*M, 0);
	for(int x=0; x<M; x++){
		uint64_t pnum = array[x];
		int j=0;
		pos_in_bv+=block_size_here;
		while(pnum!=0)
		{
			if(pnum%2 == 1){
				bv_mapping[pos_in_bv-1+j]=1;
			}
			j--;
			pnum /= 2;
		}
	}

	cout<<"size of mapping"<<endl;
	cout << "expected size ="<<(num_colors*M)/8.0 << endl;
    cout << "actual size = "<<size_in_bytes(bv_mapping) << endl;

    rrr_vector<256> rrr_bv_mapping(bv_mapping);
	cout << "rrr size = "<<size_in_bytes(rrr_bv_mapping) << endl;

	store_to_file(rrr_bv_mapping, "rrr_bv_mapping.sdsl");


	free(data);
	OutputFile opt0("bv_opt0.txt");

	u_int64_t curr_bv, prev_bv;
    int lm = ceil(log2(M));
	int lc = ceil(log2(num_colors));
    int block_size = lm;


	vector<u_int64_t> positions;
	int skipped = num_kmers/2;
    int num_bits_predicted = (lm+1)*(num_kmers-skipped)+skipped;
    // bit_vector b = bit_vector((lm+1)*(num_kmers-skipped)+skipped, 0);
    // cout<<"num_bits="<<(lm+1)*(num_kmers-skipped)+skipped<<endl;
    u_int64_t b_it=0;
    for (u_int64_t i=0; i < num_kmers; i+=1){
		getline (dup_bitmatrix_file.fs,bv_line);
		curr_bv = std::stoi(bv_line, nullptr, 2);
        
		
        if(i!=0){
            if(curr_bv == prev_bv){
                //skipped //b[b_it++]=0;
                b_it++;
            }else{
                //b[b_it]=1;
				positions.push_back(b_it);
                b_it += block_size+1;
                uint64_t num = bphf->lookup(curr_bv);

				int64_t j=0;
                 while(num!=0)
                {
                    //b[b_it-1+j] = num%2;
					if(num%2 == 1){
						positions.push_back(b_it-1+j);
					}
                    j--;
                    num /= 2;
                }
            }
        }else{
			
            //b[b_it]=1;
			positions.push_back(b_it);
                b_it += block_size+1;
                uint64_t num = bphf->lookup(curr_bv);

				int64_t j=0;
                 while(num!=0)
                {
                    //b[b_it-1+j] = num%2;
					if(num%2 == 1){
						positions.push_back(b_it-1+j);
					}
                    j--;
                    num /= 2;
                }
        } 
		prev_bv=curr_bv;
    }

	u_int64_t num_bits = b_it;
	cout<<"predicted "<<num_bits_predicted<<" actual "<<num_bits<<endl;
    //cout<<b<<endl;   

	cout<<num_bits<<"NUM"<<endl;
	bit_vector b = bit_vector(num_bits, 0);
	cout<<"init success"<<endl;
	for (u_int64_t p : positions)
	{
		//cout<<p<<endl;
		b[p]=1;
	}

	cout<<"size of matrix"<<endl;
    cout << "expected size ="<< num_bits/8.0 << endl;
    cout << "actual size = "<< size_in_bytes(b) << endl;

    rrr_vector<256> rrrb(b);
	cout << "rrr size = "<< size_in_bytes(rrrb) << endl;

    opt0.fs.close();
	
	
	store_to_file(rrrb, "rrrbv.sdsl");
	

// 	rrr_vector<256> rrr_load = rrr_vector<256>();
// 	load_from_file(rrr_load, "rrrbv.sdsl");
// 	cout<<rrr_load<<endl;


// rrr_vector<256> rrr_load2 = rrr_vector<256>();
// 	load_from_file(rrr_load2, "rrr_bv_mapping.sdsl");
// 	cout<<rrr_load2<<endl;

	
//	Decompress::start();
	//freeing up
	//free(data);
	delete bphf;
	return EXIT_SUCCESS;
}

