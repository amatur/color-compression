
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
namespace BinaryIO
{	
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
		uint64_t big_endian_value = (uint64_t)buf[7] + ((uint64_t)buf[6] << 8) + ((uint64_t)buf[5] << 16) + ((uint64_t)buf[4] << 24) + ((uint64_t)buf[3] << 32) + ((uint64_t)buf[2] << 40) + ((uint64_t)buf[1] << 48) + ((uint64_t)buf[0] << 56);
		//uint64_t little_endian_value = (uint64_t)buf[0] + ((uint64_t)buf[1] << 8) + ((uint64_t)buf[2] << 16) + ((uint64_t)buf[3] << 24) + ((uint64_t)buf[4] << 32) + ((uint64_t)buf[5] << 40) + ((uint64_t)buf[6] << 48) + ((uint64_t)buf[7] << 56);
		//cout<<big_endian_value<<endl;
		return big_endian_value;
	}
	

	/** Write bitvector to disk given positions
	 * @param b_it b_it indicates size of vector: if 8 length, then b_it = 8, bv 0 to 7 should have values
	 **/
	void write_binary_bv_from_pos_vector(vector<uint64_t>& positions, uint64_t &b_it, string filename){
		ofstream fout(filename);
		uint64_t total_num_blocks_req = ceil(b_it/8);
		char* blocks = new char[total_num_blocks_req];
		memset(blocks, 0, total_num_blocks_req*sizeof(blocks[0]));
		for(uint64_t p : positions){
			blocks[(int)floor(p/8)] |= (1<<(p%8)); //uint64_t block_id = floor(p/8); //uint8_t block_offset = p%8;
		}

		char char64_b_it[8];
		serialise_64bit(char64_b_it, b_it);
		fout.write((char *)char64_b_it, 8);

		uint64_t i = 0;
		while(i < total_num_blocks_req){
			fout.write(&blocks[i++], 1);
		}

		delete[] blocks;
		fout.close();
	}
	vector<char> read_binary_bv_into_char_array(string filename){
		ifstream infile(filename);
		char uintbuffer[8];
		infile.read ((char*)&uintbuffer, sizeof(uintbuffer));
		uint64_t b_it = deserialise_64bit(uintbuffer);
		vector<char> bv(b_it, '0');
		const int buffer_size = 3;
		uint8_t buffer[buffer_size];
		uint64_t pos = 0;

		int total_num_blocks = ceil(b_it/8);
		while(total_num_blocks > 0){
			infile.read((char *)&buffer,sizeof(buffer));
			total_num_blocks-=buffer_size;
			for(int j = 0; j<buffer_size; j++){
				uint8_t i = 1;
				while (true) {
					if(i & buffer[j]){
						bv[pos] = '1';
					}
					pos++;
					if(i==0x80){
						break;
					}
					i = i << 1; // Unset current bit and set the next bit in 'i'
				}
			}
		}
		// cout<<"writing bv: "<<endl;
		// for(uint64_t i = 0; i< b_it; i++){
		// 	cout<<bv[i];
		// }
		infile.close();
		return bv;
	}
} 
using namespace BinaryIO;

int main (int argc, char* argv[]){
	vector<uint64_t> pos;
	uint64_t b_it = 65536;

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