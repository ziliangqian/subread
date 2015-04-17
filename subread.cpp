#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <bitset>
#include <ctime>

using namespace std;

#define CODE_SIZE 16
typedef unsigned int CODE_TYPE; // 32 bit int for a code, == 16 bp nucleotide

char charset[4] = {'A','C','G','T'};
inline CODE_TYPE bitCode(char a){
    if(a=='C'||a=='c') return 0x1;
    else if(a=='G'||a=='g') return 0x2;
    else if(a=='T'||a=='t') return 0x3;
    return 0x0;
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
                //cout<<"R";
            }else{
                code_i = (read_i/CODE_SIZE)*gap_size + gap_i;
                //cout<<"F";
            }
            CODE_TYPE& code = codes[code_i];
            code = code | (bitCode(a)<<(2*code_bit_shift));
            //cout<<"NT#"<<read_i<<"\t"<<a<<"\tCODE#"<<code_i<<"\tSHT:"<<code_bit_shift<<"\t"<<bitset<32>(code)<<endl;
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
void delete_element(void * ref){
    unsigned int* p=(unsigned int*)ref;
    (*p)=0;
}
inline int element_eq(const void * p1, const void * p2){
    return (*(unsigned int*)p1)==(*(unsigned int*)p2);
}
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
struct CODE_POS_PAIR{
    CODE_TYPE code;
    unsigned int pos;
};
unsigned int myhash2(const void* elem){
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
/*fastest one*/
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
htab_t refGen_hash;
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
    refGen_hash = htab_create_alloc(1024, 
            preshhash, 
            element_eq,
            delete_element, 
            calloc,
            free);

    unsigned int count_total_hit = 0;
    unsigned int count_multiple_hit = 0;
    CODE_TYPE code = 0;
    for(unsigned int i=0; i<genome.size()-16; i+=CODE_SIZE/16){
        CODE_POS_PAIR* p_codepos = new CODE_POS_PAIR;
        p_codepos->code = encodeSubread(genome.substr(i, 16).c_str());
        p_codepos->pos = i;
        CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(refGen_hash, p_codepos);
        if(val==0){
            CODE_POS_PAIR** addr = (CODE_POS_PAIR**) htab_find_slot(refGen_hash, p_codepos, INSERT);
            (*addr) = p_codepos;
            count_total_hit++;
        }else count_multiple_hit++;
    }
    cout<<"TOTAL:"<<count_total_hit<<"\tMULTI:"<<count_multiple_hit<<"\n";
}

int main(int argc, char* argv[]){
    if( argc<=2 ){ cout<<"program <ref_genome> <fastaq>"<<endl; return 0; }

    cout<<"current time in ms\t"<<(std::time(0))<<endl;
    hashRefGenome(argv[1]);
    cout<<"current time in ms\t"<<(std::time(0))<<endl;

    unsigned int baskets[20];
    CODE_TYPE codes[400];
    // read in the fastq file
    stringstream ss( readFileAsString(argv[2]) );
    string line;
    unsigned int lc=0;
    unsigned int c_not_in_basket = 0;
    unsigned int c_hit_read = 0;
    for(;getline(ss, line, '\n');){
        if( ((++lc)%4)!=2 ) continue;
        if( line.size()>380 ) { cerr<<"WARN: read too long! take firt 384 bps"<<endl; line=line.substr(0,384); }
        std::memset(&baskets, 0, 20*sizeof(unsigned int));
        std::memset(&codes, 0, 400*sizeof(unsigned int));
        encodeRead(line, codes, CODE_SIZE);
        CODE_POS_PAIR* p_codepos = new CODE_POS_PAIR;
        for(int i=0; i<line.size()&i<400; i++){
            p_codepos->code = codes[i];
            p_codepos->pos = i;
            CODE_POS_PAIR* val = (CODE_POS_PAIR*)htab_find(refGen_hash, p_codepos);
            if( val==0 ) continue; // no genome hit, continue
            // other wise, put it into basket
            unsigned int pos = val->pos-i;
            bool inserted = false;
            for(int b_i=0; b_i<10; b_i++){
                if( baskets[b_i*2]==0 
                        || baskets[b_i*2]==pos ){
                    baskets[b_i*2]=pos;
                    baskets[b_i*2+1]++;
                    inserted = true;
                    //cout<<"insert into bsk#"<<b_i<<"\t"<<baskets[b_i*2+1]<<"\tread line#"<<lc<<endl;
                    break;
                }
            }
            if( inserted == false ) c_not_in_basket++;
        }
        unsigned int max_count = baskets[0];
        unsigned int max_count_b_i = baskets[1];
        for(int b_i=1; b_i<10; b_i++){
            if( baskets[2*b_i+1] > max_count ){
                max_count = baskets[2*b_i+1];
                max_count_b_i = b_i;
            }
        }
        if(max_count>7){
            c_hit_read++;
        }
    }

    cout<<"current time in ms\t"<<std::time(0)<<endl;
    cout<<"done with hits:"<<c_hit_read<<endl;
}

