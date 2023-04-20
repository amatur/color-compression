
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include<vector>

using namespace std;

bool little_endian = true;

inline bool is_system_little_endian()
{
    if ( htonl(47) == 47 ) {
        return false;
    // Big endian
    } else {
        return true;//higher byte higher address
    // Little endian.
    }

    const int value { 0x01 };
    const void * address = static_cast<const void *>(&value);
    const unsigned char * least_significant_address = static_cast<const unsigned char *>(address);
    return (*least_significant_address == 0x01);
}

namespace BinaryIO
{
	void write_binary_from_pos(uint8_t* bitvector, vector<uint64_t>& positions, uint64_t &b_it){
        
		int num_blocks_req = ceil(b_it/8);

		//little endian (lo=0,hi=num_blocks_req-1)
        if(little_endian){
            for(int i = 0; i<num_blocks_req; i++){
                bitvector[i] = 0x00;
            }
        }
		
	}
	

	// void to_kmer_binary_file(const std::string& outfile, vector<uint64_t> his, vector<uint64_t> los)
	// {
	// 	ofstream out(outfile, ios::out | ios::binary);

	// 	uint8_t array[8];

	// 	// k value
	// 	array[0] =  this->k        & 0xFF;
	// 	array[1] = (this->k >>  8) & 0xFF;
	// 	array[2] = (this->k >> 16) & 0xFF;
	// 	array[3] = (this->k >> 24) & 0xFF;
	// 	array[4] = (this->k >> 32) & 0xFF;
	// 	array[5] = (this->k >> 40) & 0xFF;
	// 	array[6] = (this->k >> 48) & 0xFF;
	// 	array[7] = (this->k >> 56) & 0xFF;
	// 	out.write((char *)array, 8);

	// 	// number of kmers
	// 	array[0] =  this->kmers.size()        & 0xFF;
	// 	array[1] = (this->kmers.size() >>  8) & 0xFF;
	// 	array[2] = (this->kmers.size() >> 16) & 0xFF;
	// 	array[3] = (this->kmers.size() >> 24) & 0xFF;
	// 	array[4] = (this->kmers.size() >> 32) & 0xFF;
	// 	array[5] = (this->kmers.size() >> 40) & 0xFF;
	// 	array[6] = (this->kmers.size() >> 48) & 0xFF;
	// 	array[7] = (this->kmers.size() >> 56) & 0xFF;
	// 	out.write((char *)array, 8);

	// 	for (const uint64_t kmer : this->kmers)
	// 	{
	// 		array[0] =  kmer        & 0xFF;
	// 		array[1] = (kmer >>  8) & 0xFF;
	// 		array[2] = (kmer >> 16) & 0xFF;
	// 		array[3] = (kmer >> 24) & 0xFF;
	// 		array[4] = (kmer >> 32) & 0xFF;
	// 		array[5] = (kmer >> 40) & 0xFF;
	// 		array[6] = (kmer >> 48) & 0xFF;
	// 		array[7] = (kmer >> 56) & 0xFF;
	// 		out.write((char *)array, 8);
	// 	}

	// 	out.close();
	// }
	
} 
using namespace BinaryIO;

int main (int argc, char* argv[]){
    little_endian = is_system_little_endian();
    cout<<little_endian<<endl;
}