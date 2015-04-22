#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <bitset>
#include <vector>
#include <ctime>
#include <ext/hash_map>
//#include <unordered_map>

using namespace std;
using namespace __gnu_cxx;

#define CODE_SIZE 16
typedef unsigned int CODE_TYPE; // 32 bit int for a code, == 16 bp nucleotide

char charset[4] = {'A','C','G','T'};
inline CODE_TYPE bitCode(char a){
    if(a=='C'||a=='c') return 0x1;
    else if(a=='G'||a=='g') return 0x2;
    else if(a=='T'||a=='t') return 0x3;
    return 0x0;
}
inline CODE_TYPE rc(CODE_TYPE code){
    CODE_TYPE rc_code=0;
    unsigned int total_bits = 8*sizeof(CODE_TYPE);
    for(int i=0; i<total_bits/2; i++)
        rc_code = rc_code | ((code>>(i*2) & 0x3) << (total_bits-2*i-2));
    return rc_code;
}
inline CODE_TYPE encodeSubread(const char* subread){
    unsigned int code=0;
    char a='A';
    unsigned int bitcode=0;
    for(int i=0; i<CODE_SIZE; i++){
        a=subread[i];
        bitcode=bitCode(a);
        code = code | (bitcode<<(2*i));
    }
    return code;
}
inline void decodeSubread(char* subread, unsigned int code){
    for(int i=0; i<CODE_SIZE; i++) subread[i] = charset[ ((code>>(2*i))&0x3) ];
}
inline string decodeSubread(unsigned int code){
    char nt[CODE_SIZE+1]; decodeSubread(nt, code);
    nt[CODE_SIZE]=0;
    return string(nt);
}

void encodeRead(string& read, CODE_TYPE (&codes)[400], int gap_size){
    int code_i, code_bit_shift;
    char a;
    for(int read_i=0; read_i<read.size(); read_i++){
        a = read[read_i];
        for(int gap_i=0; gap_i<=read_i && gap_i<gap_size; gap_i++){
            code_bit_shift = (read_i%CODE_SIZE) - gap_i; 
            if(read_i>=CODE_SIZE && (read_i%CODE_SIZE)<gap_i){
                code_i = (read_i/CODE_SIZE)*gap_size - gap_size + gap_i;
                code_bit_shift = (read_i%CODE_SIZE)+CODE_SIZE-gap_i;
            }else{
                code_i = (read_i/CODE_SIZE)*gap_size + gap_i;
            }
            CODE_TYPE& code = codes[code_i];
            code = code | (bitCode(a)<<(2*code_bit_shift));
        }
    }
    /*debug: dump the codes
    for(int i=0; i<5; i++){
        char de[16];
        decodeSubread(de, codes[i]);
        cout<<read.substr(i, 16)<<"\t"<<
            bitset<32>(codes[i])<<"\t"<<de<<endl;
    }*/
}

/**
 * quick hash table, adopted from gun c lib, best in performance
 */
#include "hashtab.h"
#define mix(a,b,c) \
{ \
    a -= b; a -= c; a ^= (c>>13); \
    b -= c; b -= a; b ^= (a<< 8); \
    c -= a; c -= b; c ^= ((b&0xffffffff)>>13); \
    a -= b; a -= c; a ^= ((c&0xffffffff)>>12); \
    b -= c; b -= a; b = (b ^ (a<<16)) & 0xffffffff; \
    c -= a; c -= b; c = (c ^ (b>> 5)) & 0xffffffff; \
    a -= b; a -= c; a = (a ^ (c>> 3)) & 0xffffffff; \
    b -= c; b -= a; b = (b ^ (a<<10)) & 0xffffffff; \
    c -= a; c -= b; c = (c ^ (b>>15)) & 0xffffffff; \
}
/* positions
 * 0xfca00000-0xffffffff, reserved for marking genomic repeating 16bp codes, up to 56,623,104 repeats(25% of codes are repeat)
 * ***********0xfca00001, repeated, lookup info in genomic_repeats[0]
 * ***********0xfca00002, repeated, lookup info in genomic_repeats[1]
 * 0xd0000000-0xfca00000, reserved for common variations in human genome, up to 748,683,008 snps
 * 0x00000001-0xd0000000, reserved for human genome positions, up to 3,489,660,927 genomic positions
 * 0x00000000, invalid position
 *
 * code = 0x00000000, invalid code, polyA was skept in current version
 */
struct CODE_POS_PAIR{
    CODE_TYPE code;
    unsigned int pos; /*0xffffff00-0xffffffff, repeat # 01-ff*/
    inline bool operator ==(const CODE_POS_PAIR& right){return this->code==right.code;}
    inline bool operator !=(const CODE_POS_PAIR& right){return this->code!=right.code;}
    inline bool operator ==(const unsigned int& right){return this->code==right;}
    inline bool operator !=(const unsigned int& right){return this->code!=right;}
    inline unsigned int operator %(const unsigned int& right){return this->code%right;}
};
struct DUP_POS{
    unsigned int size;
    unsigned int* entries;
    public:
    DUP_POS(){size=0; entries=0; }
    void dump(ostream& out){ cout<<"\nDUMP DUP_POS"<<endl; for(int i=0; i<size; i++){ out<<i<<"\t"<<entries[i]<<"\n"; } }
};
inline int element_eq(const void * p1, const void * p2){
    return ((CODE_POS_PAIR*)p1)->code==((CODE_POS_PAIR*)p2)->code;
}
void delete_element(void * ref){
    unsigned int* p=(unsigned int*)ref;
    (*p)=0;
}
unsigned int gcchash(const void* elem){
    const unsigned int v = ((const CODE_POS_PAIR*) elem)->code;
    unsigned a, b, c;
    a = b = 0x9e3779b9;
    a += v >> (sizeof (intptr_t) * CHAR_BIT / 2);
    b += v & (((intptr_t) 1 << (sizeof (intptr_t) * CHAR_BIT / 2)) - 1);
    c = 0x42135234;
    mix (a, b, c);
    return c;
}
/*wong's page on integers hash, second one in perfor
 * http://burtleburtle.net/bob/hash/integer.html */
inline uint32_t wonghash(const void* elem){
    uint32_t a = ((const CODE_POS_PAIR*) elem)->code;
    a = (a ^ 61) ^ (a >> 16);
    a = a + (a << 3);
    a = a ^ (a >> 4);
    a = a * 0x27d4eb2d;
    a = a ^ (a >> 15);
    return a;
}
inline unsigned int magichash(const void* elem) {
    uint32_t x = ((const CODE_POS_PAIR*) elem)->code;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x);
    return x;
}
/* fastest one
 * use 64 bit to generate 2 code once 
 * https://github.com/preshing/CompareIntegerMaps
 */
inline unsigned int preshhash(const void* elem){
    uint32_t h = ((const CODE_POS_PAIR*) elem)->code;
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}
inline unsigned int myhash(const void* elem){
    const CODE_POS_PAIR* pcode = ((const CODE_POS_PAIR*) elem);
    return (pcode->code>>4)+((pcode->code&0x4)<<24);
}
htab_t refGen_hash, refGen_hash_tier2;
DUP_POS genome_duplicates;
string genome;
string readFileAsString(string filename){
    FILE* pfile = fopen(filename.c_str(), "rb");
    fseek(pfile, 0, SEEK_END);
    unsigned long length = ftell(pfile);
    rewind(pfile);
    char* buffer = (char*)malloc(length);
    fread(buffer, 1, length, pfile);
    fclose(pfile);
    return string(buffer);
}
void hashRefGenome(string filename){
    stringstream ss( readFileAsString(filename) );
    string line;
    for(;getline(ss, line, '\n');){
        if(line[0]=='>') continue;
        genome+=line;
    }
    refGen_hash = htab_create_alloc(1024, preshhash, element_eq, delete_element, calloc, free);
    hash_map<CODE_TYPE, vector<unsigned int> > duplicates;
    unsigned int count_total_hit = 0;
    unsigned int count_multiple_hit = 0;
    unsigned int count_collisions = 0;
    unsigned int count_2nd_hit = 0;
    CODE_TYPE code = 0;
    for(unsigned int i=0; i<genome.size()-16; i+=CODE_SIZE/16){
        CODE_POS_PAIR* p_codepos = new CODE_POS_PAIR;
        p_codepos->code = encodeSubread(genome.substr(i, 16).c_str());
        p_codepos->pos = i;
        CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(refGen_hash, p_codepos);
        
        if(val==0){
            count_collisions = refGen_hash->collisions;
            CODE_POS_PAIR** addr = (CODE_POS_PAIR**) htab_find_slot(refGen_hash, p_codepos, INSERT);
            count_collisions = refGen_hash->collisions-count_collisions;
            (*addr) = p_codepos;
            if(count_collisions>=5) count_2nd_hit++;
            count_total_hit++;
        }else {
            count_multiple_hit++;
            if(duplicates.count( val->code )==0){ 
                duplicates.insert(make_pair<CODE_TYPE, vector<unsigned int> >(val->code, vector<unsigned int>())); 
                duplicates[val->code].push_back(val->pos);
            }
            duplicates[val->code].push_back(i);
        }
        // TODO, encoding reverse-complement
    }
    // to know where is the entry for a code
    cout<<"TOTAL:"<<count_total_hit<<"\tMULTI:"<<count_multiple_hit<<"\tTODO:"<<count_2nd_hit<<"\t";
    /*
     * code1_pos1|code1_pos2|0x00000000|code2_pos1|code2_pos2|...
     */
    hash_map<CODE_TYPE, unsigned int> dup_code_4_postion;
    genome_duplicates.size = duplicates.size()+count_multiple_hit*2;
    unsigned int* dup_pos = new unsigned int[genome_duplicates.size];
    memset(dup_pos, 0, sizeof(unsigned int)*genome_duplicates.size );
    unsigned int idx=0;
    for(hash_map<CODE_TYPE, vector<unsigned int> >::iterator ite = duplicates.begin(); ite!=duplicates.end(); ite++){
        dup_code_4_postion.insert(make_pair<CODE_TYPE, unsigned int>(ite->first, idx));
        for(int i=0; i<ite->second.size(); i++){
            dup_pos[idx]=ite->second[i];
            idx++;
        }
        dup_pos[idx++]=0;//0x00000000 seperator
    }
    genome_duplicates.entries = dup_pos;
    duplicates.clear();
    // point CODE_POS pair to DUP_POS
    for(int i=0; i<refGen_hash->size; i++){
        CODE_POS_PAIR* entry = (CODE_POS_PAIR*)refGen_hash->entries[i];
        if(entry==0) continue;
        CODE_TYPE code = entry->code;
        if( dup_code_4_postion.count(code) ) {
            entry->pos = dup_code_4_postion[code]+0xfca00001;
        }
    }
}
int make_codes_stat_4_genome(string filename){
    stringstream ss( readFileAsString(filename) );
    string line;
    for(;getline(ss, line, '\n');){
        if(line[0]=='>') continue;
        genome+=line;
    }
    htab_t stat_hash = htab_create_alloc(1024, preshhash, element_eq, delete_element, calloc, free);
    CODE_TYPE code = 0;
    for(unsigned int i=0; i<genome.size()-16; i+=CODE_SIZE){
        CODE_POS_PAIR* p_codepos = new CODE_POS_PAIR;
        p_codepos->code = encodeSubread(genome.substr(i, 16).c_str());
        p_codepos->pos = 1;
        CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(stat_hash, p_codepos);
        if(val==0){
            CODE_POS_PAIR** addr = (CODE_POS_PAIR**) htab_find_slot(stat_hash, p_codepos, INSERT);
            (*addr) = p_codepos;
        }else{
            val->pos++;
        }
    } 
    ofstream fout( (string(filename)+".kmerstat").c_str() );
    fout<<"code\tpos\n";
    for(int i=0; i<stat_hash->size; i++){
        CODE_POS_PAIR* val = (CODE_POS_PAIR*)stat_hash->entries[i];
        if(val==0) continue;
        fout<<val->code<<"\t"<<val->pos<<"\n";
    }
    fout.close();
    return 0;
}

int save_genome(const char* filename){
    FILE *pfile = fopen(filename, "wb");
    fwrite(refGen_hash, sizeof(htab), 1, pfile); 
    for(int i=0; i<refGen_hash->size; i++){
        CODE_POS_PAIR* pdata_entry = (CODE_POS_PAIR*) refGen_hash->entries[i];
        if( pdata_entry == HTAB_DELETED_ENTRY || pdata_entry==HTAB_EMPTY_ENTRY ){
            pdata_entry = new CODE_POS_PAIR;
            pdata_entry->code = 0;
            pdata_entry->pos = 0;
        }
        fwrite( pdata_entry, sizeof(CODE_POS_PAIR), 1, pfile );
    }
    fclose(pfile);

    pfile = fopen( (string(filename)+".dup").c_str(), "wb" );
    fwrite( &genome_duplicates, sizeof(DUP_POS), 1, pfile);
    fwrite( genome_duplicates.entries, sizeof(unsigned int)*genome_duplicates.size, 1, pfile);
    fclose(pfile);
    return 0;
}
unsigned int load_genome(const char* filename){

    cout<<"read in .genome file ["<<filename<<"]"<<endl;
    FILE *pfile = fopen( ((string(filename)+".genome").c_str()), "rb");
    if(pfile==0) { cout<<"no .gnome file. regenerating it"<<endl; return 0; }
    
    stringstream ss( readFileAsString(filename) );
    string line;
    for(;getline(ss, line, '\n');){
        if(line[0]=='>') continue;
        genome+=line;
    }
    refGen_hash = new htab;
    fread(refGen_hash, sizeof(htab), 1, pfile); 
    // restore call back functions
    // refGen_hash = htab_create_alloc(1024, preshhash, element_eq, delete_element, calloc, free);
    refGen_hash->hash_f = &preshhash;
    refGen_hash->eq_f = &element_eq;
    refGen_hash->del_f = &delete_element;
    refGen_hash->alloc_f = &calloc;
    refGen_hash->free_f = &free;
    // reset counters
    refGen_hash->searches = 0;
    refGen_hash->collisions= 0;
    // restore data
    refGen_hash->entries = new void *[refGen_hash->size];
    CODE_POS_PAIR* datablock = new CODE_POS_PAIR[refGen_hash->size];
    fread(datablock, sizeof(CODE_POS_PAIR), refGen_hash->size, pfile);
    cout<<"read in #"<<refGen_hash->size<<" entries"<<endl; 
    for(int i=0; i<refGen_hash->size; i++) {
        if(datablock[i].code==0) refGen_hash->entries[i] = 0;
        else refGen_hash->entries[i] = &datablock[i];
    }
    fclose(pfile);

    // loading dup file if existed
    pfile = fopen( ((string(filename)+".genome.dup").c_str()), "rb");
    if(pfile==0) { cout<<"no .dup file. continue without dup file "<<endl; return 1; }
    fread(&genome_duplicates, sizeof(DUP_POS), 1, pfile);
    cout<<"read in "<<genome_duplicates.size<<" dups"<<endl;
    genome_duplicates.entries = new unsigned int[genome_duplicates.size];
    fread(genome_duplicates.entries, genome_duplicates.size * sizeof(unsigned int), 1, pfile);
    
    cout<<"loading done"<<endl;
    return refGen_hash->size;
}

/* for debugging purpose
 * dump the whole hash table
 */
void htab_dump(htab_t htab){
    cout<<"\ndump the hash table"<<endl;
    printf("SIZE:%lu\n", htab->size);
    for(int i=0; i<htab->size; i++){
        cout<<"entry#"<<i<<"\t";
        CODE_POS_PAIR* pcode = (CODE_POS_PAIR*)htab->entries[i];
        if(pcode==0) {cout<<"empty\n"; continue; }
        cout<<bitset<32>(pcode->code)<<"\t"<<decodeSubread(pcode->code)<<"\t"<<
            pcode->pos<<"\t"<<(pcode->pos<=0xfca00000?' ':(pcode->pos-0xfca00001))<<"\n";
    }
}

/**
 * also encode vcf file into coding system for quick look ups
 */
struct SNP_ENTRY{
    unsigned int rsid; // rsid in dbSNP
    CODE_TYPE code; // rsid in dbSNP
    unsigned int pos; //genomic position
    unsigned char shift; // coding frame
    SNP_ENTRY(){rsid=0; code=0; pos=0; shift=0;}
};
vector<SNP_ENTRY> snp_entries;
void dump_snps(ostream& out){
    out<<"SNP DB SIZE = "<<snp_entries.size()<<endl;
    for(int i=0; i<snp_entries.size(); i++){
        out<<snp_entries[i].rsid<<"\t"<<decodeSubread(snp_entries[i].code)<<"\t"<<
            snp_entries[i].pos<<"\tshift="<<(unsigned int)snp_entries[i].shift<<"\t"
            <<genome.substr(snp_entries[i].pos-(unsigned int)snp_entries[i].shift, 16)<<"\n";
    }
}
int load_snps(string filename){
    FILE* pFile = fopen( filename.c_str(), "rb" );
    if(pFile==NULL) return 0;
    unsigned int len = 0;
    fread(&len, sizeof(unsigned int), 1, pFile);
    SNP_ENTRY* p_snp_entries = new SNP_ENTRY[len];
    fread(p_snp_entries, sizeof(SNP_ENTRY), len, pFile);
    for(int i=0; i<len; i++){
        SNP_ENTRY snp = p_snp_entries[i];
        snp_entries.push_back(snp);
        CODE_POS_PAIR* snp_code = new CODE_POS_PAIR;
        snp_code->code = snp.code;
        snp_code->pos = i+0xd0000000;
        CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(refGen_hash, snp_code);
        if(val!=NULL) { 
            cout<<"snp hit the genome, rebuild the whole SNP stack"<<endl;
            delete p_snp_entries; // release memory
            fclose(pFile);// close file
            return 0; 
        }
        CODE_POS_PAIR **entry = (CODE_POS_PAIR**)htab_find_slot(refGen_hash, snp_code, INSERT);
        *entry = snp_code;
    }
    fclose(pFile);
    delete p_snp_entries; // release memory
    return 1;
}
int encode_common_variation(string vcf_filename){
    char buffer[2000];
    if( load_snps(vcf_filename+".snp") ) return 1; // if snp already encoded
    FILE* pFile = fopen(vcf_filename.c_str(), "r");
    unsigned int c_large_ins = 0;
    unsigned int c_failure= 0;
    unsigned int lc = 0;
    for(;fgets(buffer,2000,pFile);){
        if(buffer[0]=='#') continue;
        string line(buffer);
        unsigned int sep1 = line.find('\t');
        unsigned int sep2 = line.find('\t', sep1+1);
        unsigned int sep3 = line.find('\t', sep2+1);
        unsigned int sep4 = line.find('\t', sep3+1);
        unsigned int sep5 = line.find('\t', sep4+1);
        const string& chr = line.substr(0, sep1);
        unsigned int pos = std::stoi(line.substr(sep1+1, sep2-sep1));
        const string& rsid = line.substr(sep2+1, sep3-sep2);
        const string& ref = line.substr(sep3+1, sep4-sep3-1);
        const string& alt = line.substr(sep4+1, sep5-sep4-1);
        if( (lc++)%100000==0 ) {cout<<lc<<endl;}
        //insertion large than CODE_SIZE, ignore
        if(alt.size()-1>CODE_SIZE){ c_large_ins++; continue; }
        bool inserted = false;
        // try different coding frame to identify a unique codes
        for(int i=1; i<CODE_SIZE; i++){
            CODE_POS_PAIR* snp = new CODE_POS_PAIR;
            snp->code = encodeSubread( (genome.substr(pos-i-1, i)
                        +alt
                        +genome.substr(pos+ref.size()-1,CODE_SIZE-i-alt.size()-1)).c_str() );
            snp->pos = 0xd0000000 + snp_entries.size(); // addr. space: 0xd0000000 - 0xfca0000000
            CODE_POS_PAIR* val = (CODE_POS_PAIR*) htab_find(refGen_hash, snp);
            if(val==0){
                CODE_POS_PAIR** entry = (CODE_POS_PAIR**) htab_find_slot(refGen_hash, snp, INSERT);
                *entry = snp;
                inserted=true;
                SNP_ENTRY* snp_entry = new SNP_ENTRY;
                snp_entry->rsid = std::stoi( rsid.substr(2,rsid.size()-2) );
                snp_entry->pos = pos;
                snp_entry->code = snp->code;
                snp_entry->shift = (unsigned char)i;
                snp_entries.push_back(*snp_entry);
                break;
            }
        }
        if(inserted == false) {c_failure++;}
        //too many SNPs? ignore the resting ones, or raise complains
        if( (lc+0xd0000000) > 0xfca0000000 ) {throw lc; break;}
    }
    fclose(pFile);

    // save the snp codes into flat binary file
    FILE* poutfile = fopen( (string(vcf_filename)+".snp").c_str(), "wb" );
    unsigned int len = snp_entries.size();
    fwrite( &len, sizeof(len), 1, poutfile);
    for(int i=0; i<snp_entries.size(); i++)
        fwrite( &snp_entries[i], sizeof(SNP_ENTRY), 1, poutfile);
    fclose(poutfile);

    return 0;
}
class BASKETS{
    unsigned int max_capacity; // 200 entry to max
    unsigned int size; //current size
    CODE_POS_PAIR codes[200];
    public:
    BASKETS(){max_capacity=200; size=0; }
    inline void insert(CODE_POS_PAIR& pair){ if(size<max_capacity) codes[size++] = pair; }
    inline void reset(){ std:memset(codes,0,sizeof(CODE_POS_PAIR)*200); size=0;  }
    inline unsigned int search_pos(unsigned int pos, unsigned short* snp_hits){
        unsigned int c_hit=0;
        for(int i=0; i<size; i++){
            if(codes[i].pos == pos) c_hit++;
            snp_hits[codes[i].pos-0xd0000000]++;
        }
        return c_hit;
    }
    inline unsigned int search_code(CODE_TYPE code){
        unsigned int c_hit=0;
        for(int i=0; i<size; i++){
            if(codes[i].code == code) c_hit++;
        }
        return c_hit;
    }
};
inline bool insert_basket(unsigned int* baskets, unsigned int pos){
    bool inserted = false;
    for(int b_i=0; b_i<10; b_i++){
        if( baskets[b_i*2]==0 
                || baskets[b_i*2]==pos ){
             baskets[b_i*2]=pos;
             baskets[b_i*2+1]++;
             inserted = true;
             //cout<<"insert into bsk#"<<b_i<<"\t"<<baskets[b_i*2]<<":"<<baskets[b_i*2+1]<<"\tread line#"<<lc<<"\tcode#"<<i<<endl;
             break;
        }
     }
    return inserted;
}
int main(int argc, char* argv[]){
    if( argc<=2 ){ cout<<"program <ref_genome> <fastaq>"<<endl; return 0; }

    // make_codes_stat_4_genome(argv[1]); return 0;

    cout<<"current time in ms\t"<<(std::time(0))<<endl;
    int ret = load_genome( argv[1] );
    if(ret==0) { 
        hashRefGenome(argv[1]);
        save_genome( (string(argv[1])+".genome").c_str() );
    }
    ret = encode_common_variation(argv[3]); // build snp information
    if(ret==0) save_genome( (string(argv[1])+".genome").c_str() );

    cout<<"current time in ms\t"<<(std::time(0))<<endl;
    //genome_duplicates.dump(cout);
    //htab_dump(refGen_hash);

    unsigned int baskets[20];
    CODE_TYPE codes[400];
    // read in the fastq file
    //stringstream is( readFileAsString(argv[2]) );
    ifstream is( argv[2] );
    FILE* p_input_file = fopen(argv[2], "r");
    char buffer[20000];
    ofstream fout((string(argv[2])+".subread").c_str());
    string line;
    unsigned int lc=0;
    unsigned int c_not_in_basket = 0;
    unsigned int c_hit_read = 0;
    unsigned int c_collisions = refGen_hash->collisions;
    unsigned int max_count, max_count_b_i, second_max_count, second_max_count_b_i;
    unsigned short* snp_hits = new unsigned short[snp_entries.size()];
    std::memset(snp_hits, 0, snp_entries.size()*sizeof(unsigned short));
    string read_name="";
    BASKETS snp_baskets;
    for(;fgets(buffer, 1999, p_input_file)!=NULL;){
        line = string(buffer);
        if( ((lc)%4)==0 ) read_name=line;
        if( ((++lc)%4)!=2 ) continue;
        if( line.size()>380 ) { cerr<<"WARN: read too long! take firt 384 bps"<<endl; line=line.substr(0,384); }
        std::memset(&baskets, 0, 20*sizeof(unsigned int));
        std::memset(&codes, 0, 400*sizeof(unsigned int));
        snp_baskets.reset();
        encodeRead(line, codes, CODE_SIZE);
        CODE_POS_PAIR* p_codepos = new CODE_POS_PAIR;
        for(int i=0; i<line.size()&i<400; ){
            // stop on empty code, or 16 poly'A's, or 'N's, current sequencers can not make accurate call on plyA
            if( codes[i]==0 ) { i++; continue; }
            p_codepos->code = codes[i];
            CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(refGen_hash, p_codepos);
            if( val==0 ) {i++; continue; }// no genome hit, continue
            // other wise, put it into basket
            unsigned int pos = val->pos<i?0:val->pos-i;// sometimes, query is longer than genome itself
            bool inserted = false;
            if(val->pos > 0xfca00001){ // run into non-unique code with multiple genomic hits
                for(unsigned int pos_i=val->pos-0xfca00001;;pos_i++){
                    pos = genome_duplicates.entries[pos_i];
                    if(pos==0) break;
                    //cout<<"debug: insert from multiple entry info "<<pos<<"\t"<<decodeSubread(val->code)<<"\t"<<genome.substr(pos,16)<<"\n";
                    inserted = insert_basket(baskets, pos<i?0:pos-i); // sometimes, query is longer than genome itself
                }
                //cout<<"hit multiple entries"<<endl;
            }else if(val->pos >= 0xd0000000){ // run into previous known variation
                unsigned int idx = val->pos-0xd0000000;
                pos = snp_entries[idx].pos-snp_entries[idx].shift;
                if(snp_entries[idx].pos>snp_entries[idx].shift) 
                    inserted = insert_basket(baskets, pos);
                snp_baskets.insert(*val);
            }else{
                //cout<<"debug: unique entry info "<<pos<<"\t"<<decodeSubread(val->code)<<"\t"<<genome.substr(pos,16)<<"\n";
                inserted = insert_basket(baskets, pos);
            }
            if( inserted == false ) c_not_in_basket++;
            i = (i/CODE_SIZE+1)*CODE_SIZE; // if inserted, jump to next window
        }
        max_count=0; second_max_count = 0;
        for(int b_i=0; b_i<10; b_i++){
            if( baskets[2*b_i+1] > max_count ){
                max_count = baskets[2*b_i+1];
                max_count_b_i = b_i;
            }
            if( baskets[2*b_i+1] < max_count && baskets[2*b_i+1] > second_max_count ){
                second_max_count = baskets[2*b_i+1];
                second_max_count_b_i = b_i;
            }
            //cout<<"basket#"<<b_i<<"\tpos="<<baskets[b_i*2]<<"\tcount="<<baskets[b_i*2+1]<<"\tmax_hit bkt#="<<max_count_b_i<<" count= "<<max_count<<endl;
        }
        c_hit_read++;
        //print alignment
        if(max_count>=4){ // 5*16=80bp, 1 mismatch in 100bp
            fout<<">READ ["<<read_name<<"] | hit offset ["<<baskets[2*max_count_b_i]<<"] | #hits="<<max_count<<"\n";
            fout<<line<<"\n"<<genome.substr( baskets[2*max_count_b_i],line.size() )<<"\n";
            // search for snp hits
            int ret = snp_baskets.search_pos(baskets[2*max_count_b_i], snp_hits);
        }
        if(second_max_count>=4){ // 5*16=80bp, 1 mismatch in 100bp
            fout<<">READ ["<<read_name<<"] | 2nd hit offset ["<<baskets[2*second_max_count_b_i]<<"] | #hits="<<second_max_count<<"\n";
            fout<<line<<"\n"<<genome.substr( baskets[2*second_max_count_b_i],line.size() )<<"\n";
        }
    }

    cout<<"current time in ms\t"<<std::time(0)<<endl;
    cout<<"done with hits:"<<c_hit_read<<endl;
    cout<<"hashtable collisions "<<refGen_hash->collisions-c_collisions<<endl;
    fout.close();
    fclose(p_input_file);

}

