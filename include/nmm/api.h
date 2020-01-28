#ifndef NMM_API_H
#define NMM_API_H

#ifdef _WIN32
#ifdef NMM_EXPORTS
#define NMM_API __declspec(dllexport)
#else
#define NMM_API __declspec(dllimport)
#endif
#else
#define NMM_API __attribute__((visibility("default")))
#endif

#endif
