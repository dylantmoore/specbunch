// Minimal SPI host: no Stata executable, proprietary runtime, or estimation code.
#include "stplugin.h"
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <limits>
#ifdef _WIN32
#include <windows.h>
#define API extern "C" __declspec(dllexport)
#else
#include <dlfcn.h>
#define API extern "C" __attribute__((visibility("default")))
#endif
static std::vector<std::vector<double>> columns;
static std::map<std::string,std::string> macros;
static int rows=0, matrix_rows=0, matrix_cols=0;
static std::vector<double> matrix;
static int (*call_plugin)(int,char**)=nullptr;
static ST_plugin spi{};
static int nobs(){return rows;} static int nvars(){return (int)columns.size();}
static int one(){return 1;} static ST_boolean selected(int i){return i>=1&&i<=rows;}
static ST_boolean missing(double x){return !std::isfinite(x)||x>=spi.missval;}
static int get(int j,int i,double* v){if(j<1||j>nvars()||i<1||i>rows)return 198;*v=columns[j-1][i-1];return 0;}
static int put(int j,int i,double v){if(j<1||j>nvars()||i<1||i>rows)return 198;columns[j-1][i-1]=v;return 0;}
static int save(char* k,char* v){macros[k]=v;return 0;}
static int use(char* k,char* v,int len){auto it=macros.find(k);if(it==macros.end()||int(it->second.size())>len)return 198;std::strcpy(v,it->second.c_str());return 0;}
static int output(char* s){std::fputs(s,stderr);return 0;}
static int nr(char*){return matrix_rows;} static int nc(char*){return matrix_cols;}
static int matput(char*,int r,int c,double x){if(r<1||r>matrix_rows||c<1||c>matrix_cols)return 198;matrix[(r-1)*matrix_cols+c-1]=x;return 0;}
API int host_open(const char* path){
#ifdef _WIN32
 auto h=LoadLibraryA(path);if(!h)return 1;
 auto init=reinterpret_cast<int(*)(ST_plugin*)>(GetProcAddress(h,"pginit"));
 call_plugin=reinterpret_cast<int(*)(int,char**)>(GetProcAddress(h,"stata_call"));
#else
 auto h=dlopen(path,RTLD_NOW|RTLD_LOCAL);if(!h){std::fputs(dlerror(),stderr);return 1;}
 auto init=reinterpret_cast<int(*)(ST_plugin*)>(dlsym(h,"pginit"));
 call_plugin=reinterpret_cast<int(*)(int,char**)>(dlsym(h,"stata_call"));
#endif
 if(!init||!call_plugin)return 2;
 spi.nobs=nobs;spi.nvar=nvars;spi.nvars=nvars;spi.nobs1=one;spi.nobs2=nobs;spi.selobs=selected;
 spi.safevdata=get;spi.vdata=get;spi.safestore=put;spi.store=put;
 spi.macresave=save;spi.macuse=use;spi.spouterr=output;spi.spoutsml=output;
 spi.missval=8.98846567431158e307;spi.ismissing=missing;
 spi.safematstore=matput;spi.matstore=matput;spi.rowsof=nr;spi.colsof=nc;
 return init(&spi)==SD_PLUGINVER?0:3;
}
API void host_data(int n,int k){rows=n;columns.assign(k,std::vector<double>(n));macros.clear();}
API void host_column(int j,const double* p){std::copy(p,p+rows,columns.at(j).begin());}
API void host_matrix(int r,int c){matrix_rows=r;matrix_cols=c;matrix.assign(r*c,spi.missval);}
API int host_call(int n,char** args){return call_plugin?call_plugin(n,args):999;}
API const char* host_macro(const char* k){return macros[k].c_str();}
API double host_value(int i){return matrix.at(i);}
