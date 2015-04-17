#include <iostream>
#include <fstream>
using namespace std;

#include "hashtab.h"
void delete_element(void * ref){
    unsigned int* p=(unsigned int*)ref;
    (*p)=0;
}
int element_eq(const void * p1, const void * p2){
    cout<<"compare "<<(*(unsigned int*)p1)<<"\t"<<(*(unsigned int*)p2)<<endl;
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
unsigned int myhash(const void* elem){
    const unsigned int v = *((const unsigned int*) elem);
    unsigned a, b, c;
     
    a = b = 0x9e3779b9;
    a += v >> (sizeof (intptr_t) * CHAR_BIT / 2);
    b += v & (((intptr_t) 1 << (sizeof (intptr_t) * CHAR_BIT / 2)) - 1);
    c = 0x42135234;
    mix (a, b, c);
    return c;
}

int main(){
    for(int i=0; i<1; i++) cout<<i<<endl;
    ifstream is("1.genome.fa0");
    is.seekg(0, std::ios::end);
    cout<<is.tellg()<<endl;
    is.close();

    //test the hash table function
    cout<<"prepare to create a hash"<<endl;
    htab_t ref_hash_gcc = htab_create_alloc(4l, 
            myhash, 
            element_eq,
            delete_element, 
            calloc,
            free);
    cout<<"created a hash"<<endl;
    unsigned int a[5]={1,1,1,1,10};
    for(int i=0; i<5; i++){
        unsigned int* val = (unsigned int *)htab_find(ref_hash_gcc, &a[i]);
        unsigned int pos;
        if( val!=0 ){
            cout<<"hit old entry "<<(*val)<<endl;
        }
        if( val==0 || ((*val))!=a[i]){
            unsigned int** addr = (unsigned int**)htab_find_slot(ref_hash_gcc, 
                &a[i], INSERT);
        cout<<"debug here:"<<addr<<endl;
        cout<<"new addr:"<<(&a[i])<<"\t"<<addr<<"\t"<<(*addr)<<endl;
        cout<<"new addr:"<<(&a[i])<<"\t"<<addr<<endl;
            (*addr) = &a[i];
            cout<<"***new entry for "<<a[i]<<"\tin pos\t"
                <<addr<<"\t"<<(**addr)<<endl;
        }
    }

    cout<<"testing query the hash map"<<endl;
    unsigned int b = 10;
    unsigned int* val = (unsigned int *)htab_find(ref_hash_gcc, &b);
    cout<<val<<endl;
}

