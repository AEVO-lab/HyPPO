#ifndef DEFINE_H
#define DEFINE_H

#ifdef WINDOWS
typedef __int64          int64;  //Portable signed long integer 8 bytes
typedef unsigned __int64 uint64; //Portable unsigned long integer 8 bytes
#else
typedef long long          int64; //Portable signed long integer 8 bytes
typedef unsigned long long uint64;//Portable unsigned long integer 8 byte
#endif

#endif // DEFINE_H
