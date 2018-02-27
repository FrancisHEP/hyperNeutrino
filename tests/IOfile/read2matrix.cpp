 #include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
using namespace std;
int main(){
    ifstream ifile;               //定义输入文件
    ifile.open("example");     //作为输入文件打开
    int i=0;
    while(1){
        ifile>>i;                                 //由文件读入数据
        if(ifile.eof()!=0) break;            //当读到文件结束时，ifile.eof()为真
        cout<<setw(6)<<i<<setw(10)<<endl;     //屏幕显示
    }
    ifile.close();                 //关闭文件
    //remove("example");
    
    return 0;
}

