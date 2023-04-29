
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <cmath>
#include<fstream>
#include<vector>
#include <bitset>
#include <climits>
#include <cstring>
#include <iostream>
 

using namespace std;

namespace BinaryIO
{	
// 	uint64_t little_endian_value =
//    (uint64_t)buf[0] + ((uint64_t)buf[1] << 8) + ((uint64_t)buf[2] << 16) + ... + ((uint64_t)buf[7] << 56);

// uint64_t big_endian_value =
//    (uint64_t)buf[7] + ((uint64_t)buf[6] << 8) + ((uint64_t)buf[5] << 16) + ... + ((uint64_t)buf[0] << 56);
	void serialise_64bit( char* dest, uint64_t n)
	{
		dest[0] = (n >> 56) & 0xff;
		dest[1] = (n >> 48) & 0xff;
		dest[2] = (n >> 40) & 0xff;
		dest[3] = (n >> 32) & 0xff;
		dest[4] = (n >> 24) & 0xff;
		dest[5] = (n >> 16) & 0xff;
		dest[6] = (n >>  8) & 0xff;
		dest[7] = (n >>  0) & 0xff;
	}
	uint64_t deserialise_64bit( char buf[8])
	{
		//uint64_t little_endian_value = (uint64_t)buf[0] + ((uint64_t)buf[1] << 8) + ((uint64_t)buf[2] << 16) + ((uint64_t)buf[3] << 24) + ((uint64_t)buf[4] << 32) + ((uint64_t)buf[5] << 40) + ((uint64_t)buf[6] << 48) + ((uint64_t)buf[7] << 56);
		uint64_t big_endian_value = (uint64_t)buf[7] + ((uint64_t)buf[6] << 8) + ((uint64_t)buf[5] << 16) + ((uint64_t)buf[4] << 24) + ((uint64_t)buf[3] << 32) + ((uint64_t)buf[2] << 40) + ((uint64_t)buf[1] << 48) + ((uint64_t)buf[0] << 56);
		cout<<big_endian_value<<endl;
		return big_endian_value;
	}
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

		const int value = { 0x01 };
		const void * address = static_cast<const void *>(&value);
		const unsigned char * least_significant_address = static_cast<const unsigned char *>(address);
		return (*least_significant_address == 0x01);
	}
	void write_binary_bv_from_pos_vector(vector<uint64_t>& positions, uint64_t &b_it, string filename){
		ofstream fout(filename);
        //b_it indicates size of vector: if 8 length, then b_it = 8, bv 0 to 7 should have values
		uint64_t total_num_blocks_req = ceil(b_it/8);
		char* blocks = new char[total_num_blocks_req];
		// for(uint64_t i = 0; i< b_it; i++){
		// 	blocks[i] = 0;
		// }
		memset(blocks, 0, total_num_blocks_req*sizeof(blocks[0]));
		//vector<uint8_t> blocks(total_num_blocks_req, 0);
		for(uint64_t p : positions){
			//uint64_t block_id = floor(p/8);
			//uint8_t block_offset = p%8;
			blocks[(int)floor(p/8)] |= (1<<(p%8));
		}
		//fout.write(reinterpret_cast<uint64_t *>(&b_it), sizeof(b_it));
		char char64_b_it[8];
		serialise_64bit(char64_b_it, b_it);
		fout.write((char *)char64_b_it, 8);

		uint64_t i = 0;
		while(i < total_num_blocks_req){
			fout.write(&blocks[i++], 1);
		}

		delete[] blocks;
		//little endian (lo=0,hi=num_blocks_req-1)
        // if(little_endian){
        //     for(int i = 0; i<num_blocks_req; i++){
        //         bitvector[i] = 0x00;
        //     }
        // }
		fout.close();
	}
	vector<char> read_binary_bv_into_char_array(string filename){
		ifstream infile(filename);
		char uintbuffer[8];
		cout<<"s:"<<sizeof(uintbuffer)<<endl;
		infile.read ((char*)&uintbuffer, sizeof(uintbuffer));
		uint64_t b_it = deserialise_64bit(uintbuffer);
		cout<<"b_it:"<<b_it<<endl;
		vector<char> bv(b_it, '0');
		const int buffer_size = 3;
		uint8_t buffer[buffer_size];
		uint64_t pos = 0;

		int total_num_blocks = ceil(b_it/8);
		while(total_num_blocks > 0){
			//cout<<"read"<<total_num_blocks<<endl;
			infile.read((char *)&buffer,sizeof(buffer));
			total_num_blocks-=buffer_size;
			for(int j = 0; j<buffer_size; j++){
				uint8_t i = 1;
				while (true) {
					// Unset current bit and set the next bit in 'i'
					
					if(i & buffer[j]){
						bv[pos] = '1';
					}
					
			
					// increment position
					pos++;
					if(i==0x80){
						break;
					}
					i = i << 1;
				}
			}
		}
		cout<<"writing bv: "<<endl;
		for(uint64_t i = 0; i< b_it; i++){
			cout<<bv[i];
		}
		infile.close();
		return bv;
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
	//ifstream test.txt
	vector<uint64_t> pos;
	uint64_t b_it = 65536;
	pos.push_back(0);
	pos.push_back(1);
	pos.push_back(2);
	pos.push_back(3);
	pos.push_back(4);
	pos.push_back(5);
	pos.push_back(6);

	pos.push_back(7);
		pos.push_back(9);
	//pos.push_back(70);
	//pos.push_back(20);
	pos.push_back(65534);
	pos.push_back(65535);
	write_binary_bv_from_pos_vector(pos, b_it, "example.bin");
	vector<char> spss_boundary = read_binary_bv_into_char_array("example.bin");
	//cout<<spss_boundary.size()<<endl;

}