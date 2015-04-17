
unsigned int* encodeRefGenomeRaw(string filename){
    ifstream is("chr7.fa");
    //ifstream is("1.genome.fa");
    is.seekg(0, std::ios::end);
    int length = is.tellg();
    is.seekg(0, std::ios::beg);
    char* buffer = new char[length];
    is.read(buffer, length);
    is.close();

    string line;
    stringstream ss(buffer);
    while( getline(ss, line, '\n') ){
        if(line[0]=='>') continue;
        genome += line;
    }
    
    unsigned int* p = (unsigned int*)malloc(4294967296);
    std::memset(p, 0, line.size()*4);

    for(unsigned int i=0; i<genome.size()-32; i++){
        unsigned int code = encodeSubRead( genome.substr(i, 16) );
        p[code>>2]=i;
    }

    return (unsigned int*)p;
}

void encodeRefGenomeRawAndSave(string filename){
    ifstream is(filename.c_str());
    is.seekg(0, std::ios::end);
    int length = is.tellg();
    is.seekg(0, std::ios::beg);
    char* buffer = new char[length];
    is.read(buffer, length);
    is.close();

    string line;
    stringstream ss(buffer);
    while( getline(ss, line, '\n') ){
        if(line[0]=='>') continue;
        genome += line;
    }
    unsigned int* p = (unsigned int*)malloc(4294967296);
    for(int j=0; j<4; j++){
        std::memset(p, 0, line.size()*4);
        char a='0';
        if(j==1) a='1';
        if(j==2) a='2';
        if(j==3) a='3';
        ofstream fout( (filename+a).c_str(), ios::binary );
        for(unsigned int i=0; i<genome.size()-16; i++){
            unsigned int code = encodeSubRead( genome.substr(i, 16) );
            if( (code>>30)==j ) p[code>>2]=i;
        }
        fout.write( (char*)p, 4294967296);
        fout.close();
    }
}

unsigned int* readRefGenomeRawImage(string filename){
    /* ifstream is 5 times slower than FILE *
    ifstream is(filename.c_str());
    is.seekg(0, std::ios::end);
    long length = is.tellg();
    is.seekg(0, std::ios::beg);
    char* buffer = (char*)malloc(length);
    is.read(buffer, length);
    is.close();
    */

    // 2.8s for 4G file blob (SSD)
    FILE* pfile = fopen(filename.c_str(), "rb");
    fseek(pfile, 0, SEEK_END);
    long length = ftell(pfile);
    rewind(pfile);
    char* buffer = (char*)malloc(length);
    fread(buffer, 1, length, pfile);
    fclose(pfile);

    return (unsigned int*)buffer;
}

