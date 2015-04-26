#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <bitset>
#include <vector>
#include <ctime>
#include <algorithm>
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
            //cout<<read_i<<":"<<gap_i<<":"<<code_i<<"\t"<<bitset<32>(code)<<endl;
        }
    }
    /*debug: dump the codes
    cout<<read<<" in "<<gap_size<<endl;
    for(int i=0; i<read.size(); i++){
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
 * 0xf4200000-0xffffffff, reserved for marking genomic repeating 16bp codes, up to 199M repeats(25% of codes are repeat)
 * ***********0xf4200001, repeated, lookup info in genomic_repeats[0]
 * ***********0xf4200002, repeated, lookup info in genomic_repeats[1]
 * 0xd0000000-0xf4200000, reserved for common variations in human genome, up to 600M snps
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
    //debug, a highly repeat sequence on genome
    CODE_TYPE debug_code = encodeSubread("AACCCTAACCCTAACCCTAACCCTAACCCA"); 
    for(unsigned int i=0; i<genome.size()-16; i+=CODE_SIZE/16){
        CODE_POS_PAIR* p_codepos = new CODE_POS_PAIR;
        p_codepos->code = encodeSubread(genome.substr(i, 16).c_str());
        if(p_codepos->code==debug_code) cout<<"debug code:"<<i<<"\t"<<genome.substr(i, 16)<<endl;
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
    hash_map<CODE_TYPE, unsigned int> dup_code_2_postion;
    genome_duplicates.size = duplicates.size()*2+count_multiple_hit;
    unsigned int* dup_pos = new unsigned int[genome_duplicates.size];
    memset(dup_pos, 0, sizeof(unsigned int)*genome_duplicates.size );
    unsigned int idx=0;
    for(hash_map<CODE_TYPE, vector<unsigned int> >::iterator ite = duplicates.begin(); ite!=duplicates.end(); ite++){
        dup_code_2_postion.insert(make_pair<CODE_TYPE, unsigned int>(ite->first, idx));
        if(ite->first==debug_code) cout<<"debug ordering code idx="<<idx<<endl;
        for(int i=0; i<ite->second.size(); i++){
            dup_pos[idx]=ite->second[i];
            idx++;
            if(ite->first==debug_code) cout<<"debug ordering code "<<decodeSubread(ite->first)<<"\t"<<ite->second[i]<<endl;
        }
        dup_pos[idx++]=0;//0x00000000 seperator
    }
    genome_duplicates.entries = dup_pos;
    duplicates.clear();
    // point CODE_POS pair to DUP_POS
    for(int i=0; i<refGen_hash->size; i++){
        CODE_POS_PAIR* entry = (CODE_POS_PAIR*)refGen_hash->entries[i];
        if(entry==0) continue;
        if(entry->code==debug_code) cout<<"debug code "<<decodeSubread(entry->code)<<"\t"<<entry->pos<<endl;
        CODE_TYPE code = entry->code;
        if( dup_code_2_postion.count(code) ) {
            if( 0xffffffff-0xf4200001 < dup_code_2_postion[code] ) throw dup_code_2_postion[code]; // too many dups
            entry->pos = dup_code_2_postion[code]+0xf4200001;
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
            pcode->pos<<"\t"<<(pcode->pos<=0xf4200000?' ':(pcode->pos-0xf4200001))<<"\n";
    }
}

/**
 * also encode vcf file into coding system for quick look ups
 */
struct SNP_ENTRY{
    unsigned int rsid; // rsid in dbSNP
    CODE_TYPE code; // usable codon
    unsigned int pos; //genomic position
    unsigned char offset; // SNP offset on the code
    unsigned char type; // SNV 0b01, Ins 0b10(64 size info), Del 0b11(64 size info)
    SNP_ENTRY(){rsid=0; code=0; pos=0; offset=0;type=0b01;}
};
vector<SNP_ENTRY> snp_entries;
void dump_snps(ostream& out){
    out<<"SNP DB SIZE = "<<snp_entries.size()<<endl;
    for(int i=0; i<snp_entries.size(); i++){
        out<<snp_entries[i].rsid<<"\t"<<decodeSubread(snp_entries[i].code)<<"\t"<<
            snp_entries[i].pos<<"\tshift="<<(unsigned int)snp_entries[i].offset<<"\t"
            <<genome.substr(snp_entries[i].pos-(unsigned int)snp_entries[i].offset, 16)<<"\n";
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
    cout<<"load "<<len<<" snps"<<endl;
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
    unsigned int c_tries= 0;
    unsigned int lc = 0;
    cout<<"coding snps"<<endl;
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
        //if( (lc++)%20000==0 ) {cout<<c_failure<<"/"<<lc<<"/"<<c_tries<<"\t"<<time(0)<<endl;}
        //insertion large than CODE_SIZE, ignore
        if(alt.size()-1>CODE_SIZE){ c_large_ins++; continue; }
        bool inserted = false;
        // try different coding frame to identify a unique codes
        CODE_POS_PAIR* snp = new CODE_POS_PAIR;
        cout<<"debug snp info: "<<rsid<<"\t"<<pos<<":"<<ref<<"->"<<alt<<endl;
        cout<<"debug snp context: "<<genome.substr(pos-16,16)<<"|"<<genome.substr(pos,16)<<endl;
        const string& context = genome.substr(pos-CODE_SIZE,CODE_SIZE-1)+alt+genome.substr(pos,CODE_SIZE);
        for(int i=alt.size()-1; i<context.size()-CODE_SIZE-1; i++){
            c_tries++ ;
            snp->pos = 0xd0000000 + snp_entries.size(); // addr. space: 0xd0000000 - 0xf420000000
            snp->code = encodeSubread( context.substr(i,16).c_str() );
            cout<<"debug snp: "<<context.substr(i,16)<<"\tsht:"<<i<<endl;
            CODE_POS_PAIR* val = (CODE_POS_PAIR*) htab_find(refGen_hash, snp);
            if(val==0){
                CODE_POS_PAIR** entry = (CODE_POS_PAIR**) htab_find_slot(refGen_hash, snp, INSERT);
                *entry = snp;
                inserted=true;
                SNP_ENTRY* snp_entry = new SNP_ENTRY;
                snp_entry->rsid = std::stoi( rsid.substr(2,rsid.size()-2) );
                snp_entry->pos = pos;
                snp_entry->code = snp->code;
                snp_entry->offset = i>=CODE_SIZE?15:( (unsigned char) (15-i) );
                snp_entries.push_back(*snp_entry);
                cout<<"debug snp inserted into "<<snp->pos<<endl;
                break;
            }
        }
        delete snp; //release memory if no insertion
        if(inserted == false) {
            c_failure++;
            //cout<<"snp coding failure #"<<c_failure<<". "<<rsid<<endl;
            /*the reason of failure, conflict with genome entry? or a existing SNP?
             * it generates bunch of failures
             */
        }
        //too many SNPs? ignore the resting ones, or raise complains
        if( (lc+0xd0000000) > 0xf420000000 ) {throw lc; break;}
    }
    fclose(pFile);

    cout<<"save coded snps"<<endl;

    if(c_failure>0){ cout<<"Failure rate: "<<c_failure<<"/"<<lc<<endl; }

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
    inline bool insert(CODE_POS_PAIR& pair){ 
        if(size<max_capacity) {
            cout<<"debug: insert "<<pair.pos-0xd0000000<<":"<<bitset<32>(pair.code)<<endl;
            codes[size++] = pair; 
            return true;
        }else return false;
    }
    inline void reset(){ std:memset(codes,0,sizeof(CODE_POS_PAIR)*200); size=0;  }
    inline unsigned int search_pos(unsigned int pos, unsigned short* snp_hits){
        unsigned int c_hit=0;
        for(int i=0; i<size; i++){
            SNP_ENTRY& snp_entry = snp_entries[codes[i].pos-0xd0000000];
            cout<<"debug: search test code #"<<i<<"\t"<<codes[i].pos<<"\t"<<snp_entry.pos<<"\t"<<snp_entry.offset<<"\tvs. "<<pos<<"\n";
            // TODO: minus 1 due to db entry, different coordinate system
            if(snp_entry.pos-snp_entry.offset-1 == pos) c_hit++;
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
/* tag
 * 0b00------ : 00: known variation in dbSNP
 * 0b01001010 : 01: substititution, 001: C, 010, G
 * 0b10000001 : 10: deletion, 000001, 1 bp deleted
 * 0b11------ : 11: insertion, the following bits record the length of insertion
 */
struct Variation{
    unsigned char offset;
    unsigned char tag;
    public: Variation(){offset=0; tag=0;}
};
/* memory usage: 10B in total for recording a read mapping to ref
 */
struct ReadAlignment{
    string read; // memory heavy field, TODO
    unsigned int pos; // the ReadAlignment could be ordered by pos
    unsigned char size; // alignment started from pos
    Variation var[6]; // 6 variation for a read at most
    ReadAlignment(){pos=0;}
    public: int add(Variation _var){
                if(size>=6) return false;
                var[size++] = _var;
                return true;
            }
};
bool aligned_read_comparator(ReadAlignment& read1, ReadAlignment& read2){
    return read1.pos>read2.pos;
}
/* Variations 
 */
struct GenomeVariation{
    unsigned char tag; // the tag as used in Variation
    unsigned int c_evidence[3]; //3 slot for stroing different level of evidence
    public: GenomeVariation(){std::memset(c_evidence,0,3*sizeof(unsigned int));}
};
hash_map<unsigned int, GenomeVariation> variations;
inline int addVarToGenome(int read_coordinate, ReadAlignment& alignment, unsigned int level=0){
    unsigned int pos = 0;
    for(int i=0; i<6; i++){
        Variation var = alignment.var[i];
        if(var.tag==0) break;
        pos = read_coordinate + var.offset;
        if( variations.count(pos)== 0 ) {
            GenomeVariation genomevariation;
            genomevariation.tag = var.tag; 
            genomevariation.c_evidence[level]++;
            variations.insert( make_pair<unsigned int, GenomeVariation>(pos, genomevariation) );
            return 0;
        }
        if(variations[pos].tag!=0 && variations[pos].tag==var.tag){// same pos, same variation
            variations[pos].c_evidence[level]++;
        }else if( variations[pos].tag!=var.tag){//observed different variation on same pos
            // PUSH INTO TODO LIST
            cerr<<"multiple mutations in the same site"<<endl;
        }
        return 0;
    }
    return 0;
}
/*extend the hits to get a full alignment, discover indel and snv at the same time
 * SNP is also count into MISMATCHES
 * TODO issue: known IndelSNP hits would touble the extension step due to gaps
 * */
inline unsigned int extend_hits(unsigned int pos1, unsigned int pos2, 
        const string& read, unsigned int MAX_MISMATCH, ReadAlignment& alignment){
    // search start with smaller pos
    unsigned int pos = pos1<pos2?pos1:pos2;
    string context = genome.substr(pos1, read.size());
    string context2 = genome.substr(pos2, read.size());
    std::transform(context.begin(), context.end(),context.begin(), ::toupper);
    std::transform(context2.begin(), context2.end(),context2.begin(), ::toupper);
    cout<<"debug: pos1="<<pos1<<", pos2="<<pos2<<"\nREA:"<<read<<" "<<read.size()<<"\nREF:"<<context<<"\nREF:"<<context2<<" "<<context2.size()<<endl;
    unsigned int c_mismatch = 0;
    unsigned int c_perfect_hits_1 = 0;
    unsigned int c_perfect_hits_2 = 0;
    unsigned char* match_table = new unsigned char[read.size()];
    std::memset(match_table, 0, read.size());
    unsigned int pos_sum_1 = 0;
    unsigned int pos_sum_2 = 0;
    for(int i=0; i<read.size(); i++){
        if(read[i]==context[i] ){
            match_table[i] = (match_table[i] | 1);
            c_perfect_hits_1++;
            pos_sum_1 += i;
        } else if(read[i]==context2[i]){
            match_table[i] = (match_table[i] | 2);
            c_perfect_hits_2++;
            pos_sum_2 += i;
        }
    }
    if(c_perfect_hits_1>0) pos_sum_1=pos_sum_1/c_perfect_hits_1;
    if(c_perfect_hits_2>0) pos_sum_2=pos_sum_2/c_perfect_hits_2;
    // perfect alignment, no mismatch, no gap open
    if( c_perfect_hits_1==read.size() ) {
        addVarToGenome(pos1, alignment, 0); 
        return 0;
    }
    // with several mismatches from the major hit
    if( c_perfect_hits_1>=read.size()-MAX_MISMATCH) {
        for(int i=0; i<read.size(); i++){
            if( read[i]!=context[i] && read[i]!=context2[i] ){
                cout<<"debug "<<read[i]<<":"<<context[i]<<":"<<context2[i]<<endl;
                if( (++c_mismatch) > MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches during perfect matching"<<endl;
                    delete match_table; // release memory
                    return 1; // unmapped
                }
                Variation var;
                var.offset = i;
                var.tag = 0b01000000+(bitCode(context[i])<<3)+bitCode(read[i]);
                alignment.add( var );
            }
        }
        return 0; //success with < MAX_MISMATCH
    }
    if( pos1!=pos2 && pos_sum_2==pos_sum_1 ) { delete match_table; return 2; }//unmapped;
    if( pos1==pos2 ) { delete match_table; return 3; } //unmapped;

    // moving average to see the seperating point
    cout<<"debug pos_sum: #1="<<pos_sum_1<<"\t#2="<<pos_sum_2<<endl;
    unsigned int sum1=0;
    unsigned int sum2=0;
    int sep=0;
    for(sep=0; sep<read.size()-5; sep+=3){
        sum1 = sum2; sum2=0;
        for(int j=0; j<5; j++) sum2+=match_table[sep+j];
        if(( 2*sum1/5!=2*sum2/5) && sep>0) break;
    }
    for(int i=0; i<read.size(); i++) cout<<(int)match_table[i]<<"";
    cout<<"\ndebug roughly cut="<<sep<<endl;
    int error_count = 0;
    for(sep>3?(sep-=3):sep; sep<read.size(); sep++){
        if( context2[sep]==context[sep] && context[sep]==read[sep] ) continue;
        if( pos_sum_1<pos_sum_2 ){
            if( read[sep]!=context[sep] ) break;
            continue;
        }
        if( pos_sum_2<pos_sum_1 ){
            if( read[sep]!=context2[sep] ) break;
            continue;
        }
        cerr<<"can't find the cut position for double hits"<<endl; delete match_table; return 4;
    }
    cout<<"exact cut="<<sep<<endl;

    // reset mismatch and search for indels
    c_mismatch = 0;
    // find with insertions by pos_sum and pos info
    /*    ---------------III-------------
     *    pos1-----------
     *                     pos2----------
     *pos2---------- (-shift, woule be less than pos1)
     */
    if( pos1>pos2 && pos_sum_1 < pos_sum_2 ){
        // search from left
        for(int i=0; i<sep; i++){
            if(read[i]!=context[i]){
                if( (++c_mismatch) > MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding ins case #1.1"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = i;
                var.tag = 0b01000000+(bitCode(context[i])<<3)+bitCode(read[i]);
                alignment.add( var );
                cout<<"mismatch "<<__LINE__<<endl;
            }
        }
        for(int j=read.size()-1; j>=sep+pos1-pos2; j--){
            if(read[j]!=context2[j]){
                if( (++c_mismatch) > MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding ins case #1.2"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = j;
                var.tag = 0b01000000+(bitCode(context2[j])<<3)+bitCode(read[j]);
                alignment.add( var );
                cout<<"mismatch "<<__LINE__<<endl;
            }
        }
        Variation var;
        var.offset = sep;
        var.tag = 0b10000000+((pos1-pos2)>64?64:(pos1-pos2)); // if ins is large/equal than 64, kept 64
        alignment.add(var);
        cout<<"ins at "<<sep<<"   "<<__LINE__<<endl;
    }

    // the reverse of previous case
    if( pos2>pos1 && pos_sum_2 < pos_sum_1 ){
        // search from left
        for(int i=0; i<sep; i++){
            if(read[i]!=context2[i]){
                if( (++c_mismatch) >= MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding ins #2.1"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = i;
                var.tag = 0b01000000+(bitCode(context2[i])<<3)+bitCode(read[i]);
                alignment.add( var );
                cout<<"mismatch "<<__LINE__<<endl;
            }
        }
        for(int j=read.size()-1; j>=sep+pos2-pos1; j--){
            if(read[j]!=context[j]){
                if( (++c_mismatch) >= MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding ins #2.2"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = j;
                var.tag = 0b01000000+(bitCode(context[j])<<3)+bitCode(read[j]);
                alignment.add( var );
                cout<<"mismatch "<<__LINE__<<endl;
            }
        }
        Variation var;
        var.offset = sep;
        var.tag = 0b10000000+((pos2-pos1)>64?64:(pos2-pos1)); // if ins is large/equal than 64, kept 64
        alignment.add(var);
        cout<<"ins at "<<sep<<"   "<<__LINE__<<endl;
    }

    // find deletions 
    /*   ----------------
     *pos1-------pos2----
     *   pos2---- (-shift)
     */
    if( pos1<pos2 && pos_sum_1 < pos_sum_2 ){
        // search from left
        for(int i=0; i<sep; i++){
            if(read[i]!=context[i]){
                if( (++c_mismatch) >= MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding del #1.1"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = i;
                var.tag = 0b01000000+(bitCode(context[i])<<3)+bitCode(read[i]);
                alignment.add( var );
            }
        }
        for(int j=read.size()-1; j>=sep; j--){
            if(read[j]!=context2[j]){
                if( (++c_mismatch) >= MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding del #1.2 [offset="<<j<<","<<read[j]<<","<<context2[j]<<"]"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = j;
                var.tag = 0b01000000+(bitCode(context2[j])<<3)+bitCode(read[j]);
                alignment.add( var );
            }
        }
        Variation var;
        var.offset = sep;
        var.tag = 0b11000000+((pos2-pos1)>64?64:(pos1-pos2)); // if ins is large/equal than 64, kept 64
        alignment.add(var);
        cout<<"add a del in "<<(int)var.offset<<endl;
    }
    // the reverse of previous case
    if( pos1>pos2 && pos_sum_2 < pos_sum_1 ){
        // search from left
        for(int i=0; i<sep; i++){
            if(read[i]!=context2[i]){
                if( (++c_mismatch) >= MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding del #2.1"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = i;
                var.tag = 0b01000000+(bitCode(context2[i])<<3)+bitCode(read[i]);
                alignment.add( var );
            }
        }
        for(int j=read.size()-1; j>=sep; j--){
            if(read[j]!=context[j]){
                if( (++c_mismatch) >= MAX_MISMATCH) { // stop searching when there are too many mismatches
                    cerr<<"too many mismatches finding del #2.2"<<endl;
                    delete match_table; return 1;
                }
                Variation var;
                var.offset = j;
                var.tag = 0b01000000+(bitCode(context[j])<<3)+bitCode(read[j]);
                alignment.add( var );
            }
        }
        Variation var;
        var.offset = sep;
        var.tag = 0b11000000+((pos1-pos2)>64?64:(pos1-pos2)); // if ins is large/equal than 64, kept 64
        alignment.add(var);
    }
    
    delete match_table;
    return 0;
}
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
    //AACCCTAACCCTAACCCTAACCCTAACCCA
    CODE_POS_PAIR* tst = new CODE_POS_PAIR;
    tst->code = encodeSubread("AACCCTAACCCTAACCCTAACCCTAACCCA");
    CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(refGen_hash, tst);
    cout<<val->pos<<endl; 
    ret = encode_common_variation(argv[3]); // build snp information
    //if(ret==0) save_genome( (string(argv[1])+".genome").c_str() );

    cout<<"current time in ms\t"<<(std::time(0))<<endl;
    //genome_duplicates.dump(cout);
    //htab_dump(refGen_hash);

    vector<ReadAlignment> alignments;
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
        // to upper cases
        std::transform(line.begin(), line.end(),line.begin(), ::toupper);
        // remove blanks
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
        if( ((lc)%4)==0 ) read_name=line;
        if( ((++lc)%4)!=2 ) continue;
        if( line.size()>380 ) { cerr<<"WARN: read too long! take firt 384 bps"<<endl; line=line.substr(0,384); }
        std::memset(&baskets, 0, 20*sizeof(unsigned int));
        std::memset(&codes, 0, 400*sizeof(unsigned int));
        snp_baskets.reset();
        encodeRead(line, codes, CODE_SIZE);
        CODE_POS_PAIR* p_codepos = new CODE_POS_PAIR;
        for(int i=0; (i+15)<line.size()&i<400; ){
            // stop on empty code, or 16 poly'A's, or 'N's, current sequencers can not make accurate call on plyA
            if( codes[i]==0 ) { i++; continue; }
            cout<<"debug: start process #"<<i<<"\t"<<decodeSubread(codes[i])<<endl;
            p_codepos->code = codes[i];
            CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(refGen_hash, p_codepos);
            if( val==0 ) {i++; continue; }// no genome hit, continue
            // other wise, put it into basket
            cout<<"debug pos="<<val->pos<<endl;
            unsigned int pos = val->pos<i?0:val->pos-i;// sometimes, query is longer than genome itself
            bool inserted = false;
            if(val->pos > 0xf4200001){ // run into non-unique code with multiple genomic hits
                for(unsigned int pos_i=val->pos-0xf4200001;;pos_i++){
                    pos = genome_duplicates.entries[pos_i];
                    if(pos==0) break;
                    //cout<<"debug: insert from multiple entry info "<<pos<<"\t"<<decodeSubread(val->code)<<"\t"<<genome.substr(pos,16)<<"\n";
                    inserted = insert_basket(baskets, pos<i?0:pos-i); // sometimes, query is longer than genome itself
                }
                //cout<<"hit multiple entries"<<endl;
            }else if(val->pos >= 0xd0000000){ // run into previous known variation
                unsigned int idx = val->pos-0xd0000000;
                pos = snp_entries[idx].pos-snp_entries[idx].offset-1;
                if(snp_entries[idx].pos>snp_entries[idx].offset) {
                    inserted = insert_basket(baskets, pos);
                    snp_baskets.insert(*val);
                    // additional benifits from indel SNPs
                    unsigned char type = snp_entries[idx].type;
                    if( (type>>6) == 0b10 ){ //ins
                        insert_basket(baskets, pos-(type-0b10000000));
                    }else if( (type>>6) == 0b11 ){ //del
                        insert_basket(baskets, pos+(type-0b11000000));
                    }
                    cout<<"hit snp entries: "<<snp_entries[idx].rsid<<"\t"<<pos<<"\t"<<decodeSubread(snp_entries[idx].code)<<"\t"<<i<<":"<<line.substr(i,16)<<endl;
                }
            }else{
                //cout<<"debug: unique entry info "<<pos<<"\t"<<decodeSubread(val->code)<<"\t"<<genome.substr(pos,16)<<"\n";
                inserted = insert_basket(baskets, pos);
            }
            if( inserted == false ) c_not_in_basket++;
            else i = (i/CODE_SIZE+1)*CODE_SIZE; // if inserted, jump to next window
        }
        max_count=0; second_max_count = 0;
        for(int b_i=0; b_i<10; b_i++){
            if( baskets[2*b_i+1] > max_count ){
                second_max_count = max_count;
                second_max_count_b_i = max_count_b_i;
                max_count = baskets[2*b_i+1];
                max_count_b_i = b_i;
                continue;
            }
            if( baskets[2*b_i+1] > second_max_count ){
                second_max_count = baskets[2*b_i+1];
                second_max_count_b_i = b_i;
            }
        }
        if(max_count<1) continue; // 5*16=80bp, 1 mismatch in 100bp
        c_hit_read++;

        //print top-hit alignment
        fout<<">READ ["<<read_name<<"] | hit offset ["<<baskets[2*max_count_b_i]<<"] | #hits="<<max_count<<"\n";
        fout<<line<<"\n"<<genome.substr( baskets[2*max_count_b_i],line.size() )<<"\n";
        // search for snp hits
        int ret = snp_baskets.search_pos(baskets[2*max_count_b_i], snp_hits);
        cout<<"hit "<<ret<<" SNPs"<<endl;

        // print the second best match
        if(second_max_count>=4){ // 5*16=80bp, 1 mismatch in 100bp
            fout<<">READ ["<<read_name<<"] | 2nd hit offset ["<<baskets[2*second_max_count_b_i]<<"] | #hits="<<second_max_count<<"\n";
            fout<<line<<"\n"<<genome.substr( baskets[2*second_max_count_b_i],line.size() )<<"\n";
        }

        /** generate alignment based on the info in basket
         */
        ReadAlignment alignment;
        alignment.read=line; // memory heavy! TODO get a better way to save mem usage
        if( max_count+second_max_count>0 ){ //SETTO 0 temporarily
            int ret = extend_hits( baskets[2*max_count_b_i], baskets[2*second_max_count_b_i], line, 4, alignment);
            if( ret!=0 ){
                cout<<"un-mapped"<<endl;
            }
            char* pattern = new char[line.size()+1];
            char alphabet[4] = {'.', 'X','I','-'};
            for(int idx=0; idx<line.size(); idx++) pattern[idx]='.';
            for(int i=0; i<6; i++){
                Variation variant = alignment.var[i];
                if( variant.tag==0 ) break;
                pattern[ (unsigned short)variant.offset ] = alphabet[(variant.tag>>6)];
                if( (variant.tag>>6)==0b10 ){
                    for(int j=1; j< (int)(variant.tag-0b10000000); j++)
                        pattern[ (unsigned short)variant.offset + j ] = 'I';
                }
                cout<<(unsigned short)variant.offset<<":"<<bitset<8>(variant.tag)<<"\t";
            }
            pattern[line.size()]=0;
            cout<<"\nPTN:"<<pattern<<"\n";
        }
        alignments.push_back(alignment);
        // TODO, qual is missing, which should be used for variant quality calculation
    }
    unsigned int c_hit_snp=0;
    for(int i=0; i<snp_entries.size(); i++){
        if(snp_hits[i]>0) c_hit_snp++;
    }

    // coverage statistics for each variation
    for(int i=0; i<alignments.size(); i++){
        ReadAlignment& alignment = alignments[i];
        for(unsigned offset=0; offset<alignments[i].read.size(); offset++){
            unsigned int pos = offset+alignments[i].pos;
            if(variations.count(pos)>0) variations[pos].c_evidence[2]++;
        }
    }
    for(hash_map<unsigned int, GenomeVariation>::iterator ite=variations.begin(); 
            ite!=variations.end(); ite++){
        unsigned int pos = ite->first;
        GenomeVariation& v = ite->second;
        fout<<pos<<"\t"<<bitset<8>(v.tag)<<"\t"<<v.c_evidence[0]<<"/"<<v.c_evidence[2]<<"\n";
    }

    // sort & save all alignments
    std::sort(alignments.begin(), alignments.end(), aligned_read_comparator);

    cout<<"current time in ms: "<<std::time(0)<<endl;
    cout<<"done with hits:"<<c_hit_read<<". Potentially ("<<c_hit_snp<<" snps) discovered"<<endl;
    cout<<"hashtable collisions "<<refGen_hash->collisions-c_collisions<<endl;
    fout.close();
    fclose(p_input_file);

}

