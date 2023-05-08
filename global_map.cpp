// C++ program for the above approach
  
#include <bits/stdc++.h>
#include<map>
#include<utility>
#include<vector>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;  

#define SUPPORTED_COLOR 128

    void flip_bit(string& s, int pos){
        if(s[C-pos-1] == '1')  {
            s[C-pos-1]='0';
        } else{
            s[C-pos-1]='1';
        }
    }
class ColorBitVector{
    int MAX_COLORS=128;
    int C;
    uint64_t bvs[MAX_COLORS/64];

    void load_from_string(string bv){

    }
    string convert_to_string(){

    }
    int get_hd(ColorBitVector& a, ColorBitVector& b){
        a[1],b[1]+ a[0],b[0]
    }
    vector<int> get_differing_bits(ColorBitVector& a, ColorBitVector& b){
        vector<int> differing_bits;
        uint64_t prev_bv_lo = a[0];
        uint64_t curr_bv_lo = b[0];

        uint64_t prev_bv_hi = a[1];
        uint64_t curr_bv_hi = b[1];
        for (int i_bit = 0; i_bit < 64 && i_bit < C; i_bit += 1)
        {
            if (((prev_bv_lo >> i_bit) & 1) != ((curr_bv_lo >> i_bit) & 1))
            {
                if(C<64){
                    differing_bits.push_back(i_bit); //0 to 64
                }else{
                    differing_bits.push_back(C-64+i_bit); //36 to 99
                }
            }
        }
        for (int i_bit = 64; i_bit < C; i_bit += 1)
        {
            int actual_i_bit = i_bit - 64;  // 0 to C-65, 35
            if (((prev_bv_hi >> actual_i_bit) & 1) != ((curr_bv_hi >> actual_i_bit) & 1))
            {
                differing_bits.push_back(actual_i_bit);
            }
        }
        return differing_bits;
    }
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
	void convert_binary_bv_into_string_file(string filename, string outfilename){
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
		infile.close();
		return bv;
	}
} 
using namespace BinaryIO;

// DSU data structure
// path compression + rank by union
  void write_number_at_loc_advanced_by_block_sz(vector<uint64_t> & positions, uint64_t num, uint64_t loc_advanced_by_block_sz, uint64_t block_sz){ //requires loc_advanced_by_block_sz += block_size; 
		int64_t j=0;
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

		// if(DEBUG_MODE) debug1.fs<<-j<<" "<<block_sz<<endl;
		if (j > block_sz){
			cout<<"error in block"<<endl;
		}

	}

	
	void write_number_at_loc(vector<uint64_t> & positions, uint64_t num, uint64_t block_size, uint64_t& b_it ){
		write_number_at_loc_advanced_by_block_sz(positions, num, b_it+block_size, block_size);
		b_it += block_size; //successfully written and place on next bit; if size is 2, (0,1) written, now val is 2.
    }

    
class DSU {
    int64_t* parent;
    int64_t* rank;
  
public:
    DSU(int64_t n)
    {
        parent = new int64_t[n];
        rank = new int64_t[n];
  
        for (int64_t i = 0; i < n; i++) {
            parent[i] = -1;
            rank[i] = 1;
        }
    }
  
    // Find function
    int64_t find(int64_t i)
    {
        if (parent[i] == -1)
            return i;
  
        return parent[i] = find(parent[i]);
    }
  
    // Union function
    void unite(int64_t x, int64_t y)
    {
        int64_t s1 = find(x);
        int64_t s2 = find(y);
  
        if (s1 != s2) {
            if (rank[s1] < rank[s2]) {
                parent[s1] = s2;
            }
            else if (rank[s1] > rank[s2]) {
                parent[s2] = s1;
            }
            else {
                parent[s2] = s1;
                rank[s1] += 1;
            }
        }
    }
};
  
class Graph {
    vector<vector<int64_t> > edgelist;
    int64_t V;
    std::map<std::pair<int64_t, int64_t>, int64_t> posmap;
  
public:
    Graph(int64_t V) { this->V = V; }

  
    void addEdge(int64_t x, int64_t y, int w)
    {
        int64_t tempx = x;
        int64_t tempy = y;
        if (y<x){
            tempx = y;
            tempy = x;
        }
        
        if(posmap.count(make_pair(tempx,tempy)) > 0){
            if(w<edgelist[posmap[make_pair(tempx,tempy)]][0]){
                edgelist[posmap[make_pair(tempx,tempy)]][0] = w;
            }
        }else{
            posmap[make_pair(tempx,tempy)] = edgelist.size();
            edgelist.push_back({ w, tempx, tempy });
        }
        
        

    }
  
    void kruskals_mst()
    {
        // 1. Sort all edges
        sort(edgelist.begin(), edgelist.end());

        // set< vector<int> > already_seen;


        // Initialize the DSU
        DSU s(V);
        int64_t ans = 0;
        cout << "Following are the edges in the "
                "constructed MST"
             << endl;
        for (auto edge : edgelist) {
            int64_t w = edge[0];
            int64_t x = edge[1];
            int64_t y = edge[2];

  
            // Take this edge in MST if it does
            // not forms a cycle

            //  auto it = already_seen.find({x,y});
            //  auto it2 = already_seen.find({x,y});  

            
            //if ( it == already_seen.end() && it2 == already_seen.end() ) {  
                if (s.find(x) != s.find(y)) {
                    s.unite(x, y);
                    ans += w;
                    cout << x << " -- " << y << " == " << w << endl;
                    // already_seen.insert({x,y});
                    // already_seen.insert({y,x});

                }   
            //}   

           
        }
  
        cout << "Minimum Cost Spanning Tree: " << ans;
    }
};
  

int hammingDistance (uint64_t x, uint64_t y) {
    uint64_t res = x ^ y;
    return __builtin_popcountll (res) ;
}


double get_rrr_bv_MB(vector<uint64_t>& positions, uint64_t bv_size){
    bit_vector bv(bv_size, 0);
    uint64_t lastp = 0;
    for (uint64_t p: positions){
        bv[p] = 1;
    }
    rrr_vector<256> rrr_bv(bv);
	return size_in_bytes(rrr_bv)/1024.0/1024.0;
}

bit_vector store_as_sdsl(vector<uint64_t>& positions, uint64_t bv_size, string filename){
    
    //bit_vector bv = bit_vector(bv_size, 0);
    bit_vector bv(bv_size, 0);
    uint64_t lastp = 0;
    for (uint64_t p: positions){
        bv[p] = 1;
    }
    if(filename=="rrr_main"){
        //if(DEBUG_MODE) debug2.fs<<bv;
    }
    rrr_vector<256> rrr_bv(bv);
    //cout << "rrr_MB_bv_mapping="<<size_in_bytes(rrr_bv_mapping)/1024.0/1024.0 << endl;
    store_to_file(rrr_bv, filename);	//"rrr_bv_mapping.sdsl"

    return bv;
}
void write_one(vector<uint64_t> & positions, uint64_t& b_it ){
    positions.push_back(b_it);
    b_it+=1;
}

void write_zero(vector<uint64_t> & positions, uint64_t& b_it ){
    b_it+=1;
}

uint64_t[] string_to_uint_arr(string uniq_ms_line, int C){
    uint64_t colors[int(ceil(SUPPORTED_COLOR/64))];
    uint64_t lo = std::stoull(uniq_ms_line.substr(0,std::min(64,int(C))), nullptr, 2) ; 
    // if(i==0){
    //     write_number_at_loc(positions_hd, lo, min(64, C), b_it_hd ); 
    // }

    uint64_t hi = 0;
    if(C > 64){
        string ss = uniq_ms_line.substr(64,(C-64));
        hi  = std::stoull(ss, nullptr, 2);
        // if(i==0){
        //     write_number_at_loc(positions_hd, hi, C-64, b_it_hd ); 
        // }
    }
    colors[1] = hi;
    colors[0] = lo;
    return colors;
}

void MST_global(int M, int C, ifstream& cmp_keys_fs, ifstream&  hds_fs){
    Graph g(M+1);

    string prev="";
    string curr="";
    int hdsum=0;

    uint64_t prev_num;
    uint64_t curr_num;

    string line;
    // InputFile cmp_keys("cmp_keys");
    // InputFile hds("hds");
    // InputFile uniq_ms("uniq_ms.txt");
   
    uint64_t zero_64bit = 0;
    g.addEdge(curr_num, prev_num, hd);
    // for(int i =0; i<M; i++){
    //     int hdzero = hammingDistance(hi, zero_64bit) + hammingDistance(lo, zero_64bit);
    //     g.addEdge(i, M, hd);
    // }
    
    //     g.addEdge(i, M, hd); //M th entry is all zero, 0 to M-1 entry is non zero
      
    uint64_t count = 0;
    while (getline(cmp_keys_fs,curr )){
        getline(hds_fs,hd_line );
        if(count != 0){
            if(curr!=prev){
                curr_num = std::stoi(curr);
                prev_num = std::stoi(prev);
                int hd = std::stoi(hd_line);
                hdsum+=hd;
                g.addEdge(curr_num, prev_num, hd);

                //uint64_t colors[2];
                //string_to_uint_arr(curr_num, C);
                //g.addEdge(curr_num, M, hammingDistance(0, colors[1]) + hammingDistance(0, colors[0]));

            }
            
        }
        count+=1;
        prev = curr;
	}

    // Function call
    g.kruskals_mst();
    uint64_t space_estimate = M*(log(1+M))+hdsum*(loc(C))+hdsum*1;
    
}
void test_decompression(int M, int C, string rrr_map_hd_filename = "rrr_map_hd", string rrr_map_hd_boundary_filename="rrr_map_hd_boundary"){
    // bit_vector<> b(10000000, 0);
    // b[b.size()/2] = 1;
    // sd_vector<> sdb(b);
    // store_to_file(sdb, "sdb.sdsl");
    // sdb = sd_vector<>();
    // cout << sdb.size() << endl; // 0
    // load_from_file(sdb, "sdb.sdsl");
    // cout << sdb.size() << endl; // 10000000 

    rrr_vector<256> rrr_map_hd;      
    load_from_file(rrr_map_hd, rrr_map_hd_filename);  //sdsl namespace
       
    rrr_vector<256> rrr_map_hd_boundary;      
    load_from_file(rrr_map_hd_boundary, rrr_map_hd_boundary_filename);  //sdsl namespace
     

    int lc = ceil(log2(C));
    stringstream buffer;
    buffer << rrr_map_hd;
    string hd_M = buffer.str();
    uint64_t b_it = 0;
    
     size_t ones = rrr_vector<256>::rank_1_type(&rrr_map_hd_boundary)(rrr_map_hd.size()); 
    rrr_vector<256>::select_1_type rrr_hd_sel(&rrr_map_hd_boundary);

    

// cout << rrr_map_hd << endl;
    if(ones!=M){
        cerr<<"Not matching"<<endl;
        exit(4);
    }

    size_t prev_begin = 0;
    for (size_t i=1; i <= M; ++i){
        size_t block_len = rrr_hd_sel(i) - prev_begin + 1;
        size_t numblocks = (block_len / lc);
        vector<int> flip_loc;
        while(numblocks){
            flip_loc.push_back(read_uint(hd_M, b_it, lc));
            numblocks--;
        }
        for (int f: flip_loc){
            //global_table[i-1] = string(col_vector);
        }
   
    }
   
    
    
        //for each select position
        //end = sel
        //read substr(0, end-start+1)
        //read bit<size
    }
    //decompressed global table
    //check diff of uniq_ms.txt, dec_global_table
}
void NONMST_global(string uniq_ms_filename, int M, int C){
    ifstream uniq_ms(uniq_ms_filename);
    string hd_line;
    
    vector<uint64_t> positions;
    uint64_t b_it = 0;

    vector<uint64_t> positions_hd;
    uint64_t b_it_hd = 0;

    uint64_t zero_64bit = 0;
    uint64_t prev_lo = 0;
    uint64_t prev_hi = 0;

    int hdsum = 0;
    for(int i = 0 ; i< M; i++){
        string uniq_ms_line;
        getline(uniq_ms, uniq_ms_line);
        uint64_t lo = std::stoull(uniq_ms_line.substr(0,std::min(64,int(C))), nullptr, 2) ; 
        // if(i==0){
		//     write_number_at_loc(positions_hd, lo, min(64, C), b_it_hd ); 
        // }

        uint64_t hi = 0;
        if(C > 64){
            string ss = uniq_ms_line.substr(64,(C-64));
            hi  = std::stoull(ss, nullptr, 2);
            // if(i==0){
		    //     write_number_at_loc(positions_hd, hi, C-64, b_it_hd ); 
            // }
        }

        int hd = hammingDistance(hi, zero_64bit) + hammingDistance(lo, zero_64bit);
        g.addEdge(i, M, hd); //M th entry is all zero, 0 to M-1 entry is non zero
        if(true){ //i!=0
            int hd_prev = hammingDistance(hi, prev_hi) + hammingDistance(lo, prev_lo);
            hdsum += hd_prev;
            int lc = ceil(log2(C));
            for (int i_bit = 0; i_bit < 64 && i_bit < C; i_bit += 1)
            {
                if (((prev_lo >> i_bit) & 1) != ((lo >> i_bit) & 1))
                {       
                    write_number_at_loc(positions_hd, i_bit, lc, b_it_hd); // i_bit is the different bit loc
                }
            }
            for (int i_bit = 64; i_bit < C; i_bit += 1)
            {
                int actual_i_bit = i_bit - 64;
                if (((prev_hi >> actual_i_bit) & 1) != ((hi >> actual_i_bit) & 1))
                {
                    write_number_at_loc(positions_hd, i_bit, lc, b_it_hd); // i_bit is the different bit loc
                }
            }
            // for(int ii = 0; ii< (hd_prev*lc)-1; ii++){
            //     write_zero(positions, b_it); // i_bit is the different bit loc
            // }
            // write_one(positions, b_it); // i_bit is the different bit loc
            if(hd_prev!=0){
                b_it+=hd_prev*lc-1;
                write_one(positions, b_it);
            }
            if(i!=0) g.addEdge(i, i-1, hd_prev);
        }
        prev_lo=lo;
        prev_hi=hi;
    }
    uniq_ms.close();

    //store_as_sdsl(positions_hd, b_it_hd, "rrr_map_hd" );	
    write_binary_bv_from_pos_vector( positions_hd, b_it_hd, "rrr_map_hd" );
    store_as_sdsl(positions, b_it, "rrr_map_hd_boundary" );	

}

int main (int argc, char* argv[]){
    cout<<"version: dup erase"<<endl;
	vector<string> args(argv + 1, argv + argc);
    int M, C;
    string uniq_ms_filename;
	uint64_t num_kmers=0;
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: tool -c <num-colors> -m <M> -i <input-dedup-sorted>" << endl;
            return 0;
        }else if (*i == "-c") {
            C = std::stoi(*++i);
        }else if (*i == "-m") {
            M = std::stoi(*++i);
        }else if (*i == "-i") {
            uniq_ms = (*++i);
        }
    }

    Graph g(M+1);

    string prev="";
    string curr;

    uint64_t prev_num;
    uint64_t curr_num;

    string line;
    NONMST_global(uniq_ms_filename, M, C);
   

    return EXIT_SUCCESS;
}