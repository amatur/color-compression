// C++ program for the above approach
  
#include <bits/stdc++.h>
#include<map>
#include<utility>
#include<vector>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;  

#define SUPPORTED_COLOR 128

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

    store_as_sdsl(positions_hd, b_it_hd, "rrr_map_hd" );	
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