#ifndef _WPSDLL_H_
#define _WPSDLL_H_
#ifdef WPSDLL_EXPORTS  
#define WPSDLL_API __declspec(dllexport)   
#else  
#define WPSDLL_API __declspec(dllimport)   
#endif  
WPSDLL_API int wpsmain(int argc, char* argv[]);
#endif

WPSDLL_API int wpsmain(int argc, char*argv[]);