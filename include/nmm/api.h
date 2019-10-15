#ifndef NMM_API_H_API
#define NMM_API_H_API

#ifdef _WIN32
  #ifdef NMM_API_EXPORTS
    #define NMM_API __declspec(dllexport)
  #else
    #define NMM_API __declspec(dllimport)
  #endif
#else
  #define NMM_API
#endif

#endif

