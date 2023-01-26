#include <cmph.h>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include <stdio.h>
  #include <string.h>
  #include <string.h>
  // Create minimal perfect hash function from in-memory vector
using namespace std;

namespace CMPH{

    void create_table_from_vector(unsigned int nkeys, cmph_t *hash2, char** vv){
        cmph_io_adapter_t *source = cmph_io_vector_adapter((char **)vv, nkeys);
        //Create minimal perfect hash function using the brz algorithm.
        cmph_config_t *config = cmph_config_new(source);
        cmph_config_set_algo(config, CMPH_CHM);
        cmph_config_set_mphf_fd(config, mphf_fd);
        hash2 = cmph_new(config);
        cmph_config_destroy(config);
        cmph_io_vector_adapter_destroy(source);   
	}

    unsigned int lookup(string str, cmph* hash_cmph){	
		const char *key = str.c_str(); 
		//Find key
		unsigned int id = cmph_search(hash_cmph, key, (cmph_uint32)strlen(key));
		// fprintf(stderr, "Id:%u\n", id);
		//Destroy hash
		//cmph_destroy(hash);
		return id;
	}
}
using namespace CMPH;

  int main(int argc, char **argv)
  { 
  
      // Creating a filled vector
      
      char **vv ;
      vv = new char**;
      vv[0]="1";
      vv[1]="10";
      vv[2]="20";
      unsigned int nkeys = 3;

      cmph_t *hash_cmph_new = NULL;
      create_table_from_vector(nkeys, hash_cmph_new, vv);
      
      unsigned int i = 0;
      while (i < nkeys) {
          char *key = vv[i];
          unsigned int id = cmph_search(hash_cmph_new, key, (cmph_uint32)strlen(key));
          fprintf(stderr, "key:%s -- hash:%u\n", key, id);
          i++;
      }
  
      //Destroy hash
      cmph_destroy(hash_cmph_new); 
      return 0;
  }