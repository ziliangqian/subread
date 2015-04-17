#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <bitset>

using namespace std;

int main(){

    vector<string> inputs;
    string line;
    int lc = 0;
    while( getline(cin, line) ){
        lc++;
        if( (lc%4)==0 ) inputs.push_back(line);
    }

    for(vector<string>::iterator ite=inputs.begin(); ite!=inputs.end(); ite++){
        string line = *ite;
        string compressed = ""; //full 8bits qual codes
        char pre_char = line[0];
        unsigned short len = 1;
        for(int i=1; i<line.size(); i++){
            ostringstream converter;
            char cur_char = line[i];
            if( cur_char==pre_char ) len++;
            else {
                converter<<len;
                compressed = compressed + pre_char + converter.str();
                len=1;
                pre_char=cur_char;
            }
        }

        string compressed2 = ""; //1,2,3,4 2bits codes, 1/64-1
        for(int i=1; i<line.size(); i++){
            ostringstream converter;
            char cur_char = line[i];
            char pre_char_infoloss = (pre_char-'#'+4)/10;
            char cur_char_infoloss = (cur_char-'#'+4)/10;
            if( cur_char_infoloss==pre_char_infoloss ) len++;
            else {
                converter<<len;
                compressed2 = compressed2 + (char)(10*pre_char_infoloss+'#'-4) + converter.str()[0];
                pre_char=cur_char;
                len=1;
            }
        }
        double ratio = ((double)compressed2.size())/line.size();
        cout<<ratio<<endl;
//        if(ratio>0.8) cout<<line<<"\n"<<compressed2<<endl;
    }

    return 0;
}

