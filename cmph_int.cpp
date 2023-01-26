#include "BooPHF.h"
#include <cmph.h>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include <stdio.h>
#include <string.h>
#include <string.h>

  // Create minimal perfect hash function from in-memory vector
using namespace std;

    
namespace BPHF{
    typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
    typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;
}   
using namespace BPHF; 

int main(){
    int nthreads = 8;
    int nelem = 10;
    int *data = (int * ) calloc(nelem,sizeof(int));
	for (u_int64_t i = 0; i < nelem; i++){
       data[i] = i*100;
	}


    boophf_t * bphf = NULL;
	double t_begin,t_end; struct timeval timet;
	printf("Construct a BooPHF with  %lli elements  \n",nelem);
	gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
	
	// mphf takes as input a c++ range. A simple array of keys can be wrapped with boomphf::range
	// but could be from a user defined iterator (enabling keys to be read from a file or from some complex non-contiguous structure)
	auto data_iterator = boomphf::range(static_cast<const int*>(data), static_cast<const int*>(data+nelem));
	double gammaFactor = 4.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query

	//build the mphf
	bphf = new boomphf::mphf<int,hasher_t>(nelem,data_iterator,nthreads,gammaFactor);
	
	gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);	
	printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", nelem,t_end - t_begin;);
	//printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/nelem);
	
	//query mphf like this
	int  idx = bphf->lookup(data[0]);
	printf(" example query  %lli ----->  %llu \n",data[0],idx);
    return 0;
}

// namespace CMPH{

//     void create_table_from_vector(unsigned int nkeys, cmph_t *hash2, char** vv){
//         cmph_io_adapter_t *source = cmph_io_vector_adapter((char **)vv, nkeys);
//         //Create minimal perfect hash function using the brz algorithm.
//         cmph_config_t *config = cmph_config_new(source);
//         cmph_config_set_algo(config, CMPH_CHM);
//         hash2 = cmph_new(config);
//         cmph_config_destroy(config);
//         cmph_io_vector_adapter_destroy(source);   
// 	}

//     unsigned int lookup(string str, cmph_t* hash_cmph){	
// 		const char *key = str.c_str(); 
// 		//Find key
// 		unsigned int id = cmph_search(hash_cmph, key, (cmph_uint32)strlen(key));
// 		// fprintf(stderr, "Id:%u\n", id);
// 		//Destroy hash
// 		//cmph_destroy(hash);
// 		return id;
// 	}
// }
// using namespace CMPH;

//   int main(int argc, char **argv)
//   { 
  
//       // Creating a filled vector
      
//       char **vv ;
//       vv = new char*;
//       vv[0]="1";
//       vv[1]="10";
//       vv[2]="20";
//       unsigned int nkeys = 3;

//       cmph_t *hash_cmph_new = NULL;
//       create_table_from_vector(nkeys, hash_cmph_new, vv);
      
//       unsigned int i = 0;
//       while (i < nkeys) {
//           char *key = vv[i];
//           unsigned int id = cmph_search(hash_cmph_new, key, (cmph_uint32)strlen(key));
//           fprintf(stderr, "key:%s -- hash:%u\n", key, id);
//           i++;
//       }
  
//       //Destroy hash
//       cmph_destroy(hash_cmph_new); 
//       return 0;
//   }