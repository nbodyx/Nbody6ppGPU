//==================================================//
// Define some useful tools:                        
//-1. type transform:-------------------------------//
//  D: to<type T>(const type U& value);
//     tostr(const int/float/double& value);
//  I: transform value from type U to type T;
//     transform value to std::string;
//-2. exchange value: ------------------------------//
//  D: SWAP(T &a, T &b);
//  I: exchange a and b value;
//-3. quicksort :-----------------------------------//
//  D: 1) quicksort (T arr[], int m);
//     1.s) quicksort (vector<T> &arr, int m=-1); 
//     2) quicksort (T arr[], T brr[], int m);
//     2.s) quicksort (vector<T> &arr, vector<T> &brr,
//                    int m=-1); 
//  I: 1) quicksort arr with increasing order
//        m define the array region (0,m-1) to sort.
//        Should be carefule not out of array size!
//     1.s) same sort for vector, m is the region but
//          can be given any integer.
//     2) only sort arr as 1), and change the
//        order of brr synchronically.
//     (Press, et al. Numerical Recipe, 3rd. Edtion.
//        8.2 quicksort)
//-4. straight insertion sort=======================//
//  D: piksrt(T &arr[]);
//  I: straight insertion sort with increasing order.
//     Better for small number array.
//     (Press, et al. Numerical Recipe, 3rd. Edition.
//        8.1 straight insertion and shell's method)
//-5. find peak=====================================//
//  D: findpeak(T arr[], int m, compare cmp);
//  I: Find peak value index in arr based on compare
//-6. find max======================================//
//  D: findmax(T arr[], int m);
//  I: Find max value index in arr
//-7. find min======================================//
//  D: findmin(T arr[], int m);
//  I: Find min value index in arr
//-8. File IO=======================================//
//  D: mkdir(string path);
//  I: make path
////**  D: fileio(FILE* stream,string name,string opt);
////**  I: open file (name) to stream, opt: r/w/a+...
//==================================================//

#ifndef uftools_h
#define uftools_h

#include <iostream>
#include <cstdio>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <string.h>

/*
@(#)File:           $RCSfile: mkpath.c,v $
@(#)Version:        $Revision: 1.13 $
@(#)Last changed:   $Date: 2012/07/15 00:40:37 $
@(#)Purpose:        Create all directories in path
@(#)Author:         J Leffler
@(#)Copyright:      (C) JLSS 1990-91,1997-98,2001,2005,2008,2012
*/

/*TABSTOP=4*/

#include <errno.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */

//Declarion=========================================//
template <typename T>
T to(const int& value);

template <typename T>
T to(const double& value);

template <typename T>
T to(const std::string& value);

template <typename T>
T to(char*& value);

std::string tostr(const int& value);
std::string tostr(const double& value);
std::string tostr(const float& value);

template <typename T>
void SWAP(T &a, T &b);

template <typename T>
void quicksort(T arr[], int m);

template <typename T>
void quicksort(std::vector<T> &arr, int m=-1);

template <typename T>
void quicksort(T arr[], T brr[], int m);

template <typename T> 
void quicksort(std::vector<T> &arr, std::vector<T> &brr, int m=-1);

template <typename T>
void piksrt(T arr[],int m);

//typedef template <typename T>  bool (*compare) (const T &, const T &);

//template <typename T>
//int findpeak(T arr[], int m, compare cmp);

template <typename T>
int findmax(T arr[], int m);

template <typename T>
int findmin(T arr[], int m);

//1. type transform=================================//
template <typename T>
T to(const int& value)
{
  std::stringstream lw_sstr;
  T temp;
  lw_sstr.str("");
  lw_sstr.clear();
  lw_sstr<<value;
  lw_sstr>>temp;
  return temp;
}
template <typename T>
T to(const double& value)
{
  std::stringstream lw_sstr;
  T temp;
  lw_sstr.str("");
  lw_sstr.clear();
  lw_sstr<<value;
  lw_sstr>>temp;
  return temp;
}
template <typename T>
T to(const std::string& value)
{
  std::stringstream lw_sstr;
  T temp;
  lw_sstr.str("");
  lw_sstr.clear();
  lw_sstr<<value;
  lw_sstr>>temp;
  return temp;
}
template <typename T>
T to(char*& value)
{
  std::stringstream lw_sstr;
  T temp;
  lw_sstr.str("");
  lw_sstr.clear();
  lw_sstr<<value;
  lw_sstr>>temp;
  return temp;
}

std::string tostr(const int& value){
  std::stringstream lw_sstr;
  std::string temp;
  lw_sstr.str("");
  lw_sstr.clear();
  lw_sstr<<value;
  lw_sstr>>temp;
  return temp;
}

std::string tostr(const double& value){
  std::stringstream lw_sstr;
  std::string temp;
  lw_sstr.str("");
  lw_sstr.clear();
  lw_sstr<<value;
  lw_sstr>>temp;
  return temp;
}

std::string tostr(const float& value){
  std::stringstream lw_sstr;
  std::string temp;
  lw_sstr.str("");
  lw_sstr.clear();
  lw_sstr<<value;
  lw_sstr>>temp;
  return temp;
}

template <typename T>
void SWAP(T &a, T &b)
{
  T temp=a;
  a=b;
  b=temp;
}

template <typename T>
void quicksort(T arr[], int m)
{
  static const int M=7, NSTACK=64;
  int i,ir,j,k,jstack=-1,l=0;
  T a;
  int istack[NSTACK];
  ir=m-1;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        a=arr[j];
        for (i=j-1;i>=l;i--) {
          if (arr[i] <= a) break;
          arr[i+1]=arr[i];
        }
        arr[i+1]=a;
      }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a); 
        do j--; while (arr[j] > a); 
        if (j < i) break;
        SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      jstack += 2;
      if (jstack >= NSTACK) throw("NSTACK too small in sort."); if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

template <typename T>
void quicksort(std::vector<T> &arr, int m)
{
  static const int M=7, NSTACK=64;
  int i,ir,j,k,jstack=-1,l=0,n=arr.size();
  T a;
  std::vector<int> istack(NSTACK);
  if (m>0) n = std::min(m,n);
  ir=n-1;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        a=arr[j];
        for (i=j-1;i>=l;i--) {
          if (arr[i] <= a) break;
          arr[i+1]=arr[i];
        }
        arr[i+1]=a;
      }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a); 
        do j--; while (arr[j] > a); 
        if (j < i) break;
        SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      jstack += 2;
      if (jstack >= NSTACK) throw("NSTACK too small in sort."); if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

template <typename T>
void quicksort(T arr[], T brr[], int m)
{
  static const int M=7, NSTACK=64;
  int i,ir,j,k,jstack=-1,l=0;
  T a,b;
  int istack[NSTACK];
  ir=m-1;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        a=arr[j];
        b=brr[j];
        for (i=j-1;i>=l;i--) {
          if (arr[i] <= a) break;
          arr[i+1]=arr[i];
          brr[i+1]=brr[i];
        }
        arr[i+1]=a;
        brr[i+1]=b;
      }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1]);
      SWAP(brr[k],brr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
        SWAP(brr[l],brr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
        SWAP(brr[l+1],brr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
        SWAP(brr[l],brr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;) {
        do i++; while (arr[i] < a); 
        do j--; while (arr[j] > a); 
        if (j < i) break;
        SWAP(arr[i],arr[j]);
        SWAP(brr[i],brr[j]);
      }
      arr[l+1]=arr[j];
      brr[l+1]=brr[j];
      arr[j]=a;
      brr[j]=b;
      jstack += 2;
      if (jstack >= NSTACK) throw("NSTACK too small in sort."); if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

template <typename T> 
void quicksort(std::vector<T> &arr, std::vector<T> &brr, int m)
{
  static const int M=7, NSTACK=64;
  int i,ir,j,k,jstack=-1,l=0,n=arr.size();
  T a,b;
  std::vector<int> istack(NSTACK);
  if (m>0) n = std::min(m,n);
  ir=n-1;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        a=arr[j];
        b=brr[j];
        for (i=j-1;i>=l;i--) {
          if (arr[i] <= a) break;
          arr[i+1]=arr[i];
          brr[i+1]=brr[i];
        }
        arr[i+1]=a;
        brr[i+1]=b;
      }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1]);
      SWAP(brr[k],brr[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
        SWAP(brr[l],brr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
        SWAP(brr[l+1],brr[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
        SWAP(brr[l],brr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;) {
        do i++; while (arr[i] < a); 
        do j--; while (arr[j] > a); 
        if (j < i) break;
        SWAP(arr[i],arr[j]);
        SWAP(brr[i],brr[j]);
      }
      arr[l+1]=arr[j];
      brr[l+1]=brr[j];
      arr[j]=a;
      brr[j]=b;
      jstack += 2;
      if (jstack >= NSTACK) throw("NSTACK too small in sort."); if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

template <typename T>
void piksrt(T arr[], int m)
{
  int i,j;
  T a;
  for (j=1;j<m;j++)
  {
    a=arr[j];
    i=j;
    while (i>0 && arr[i-1]>a)
    {
      arr[i]=arr[i-1];
      i--;
    }
    arr[i]=a;
  }
}

template <typename T>
void piksrt(T arr[], T brr[], int m)
{
  int i,j;
  T a,b;
  for (j=1;j<m;j++)
  {
    a=arr[j];
    b=arr[j];
    i=j;
    while (i>0 && arr[i-1]>a)
    {
      arr[i]=arr[i-1];
      brr[i]=brr[i-1];
      i--;
    }
    arr[i]=a;
    brr[i]=b;
  }
}


//File IO===========================================//
typedef struct stat Stat;

static int do_mkdir(const char *path, mode_t mode)
{
    Stat            st;
    int             status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist. EEXIST for race condition */
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}

// /**
// ** mkpath - ensure all directories in path exist
// ** Algorithm takes the pessimistic view and works top-down to ensure
// ** each directory in path exists, rather than optimistically creating
// ** the last element and working backwards.
// */
 int mkpath(const char *path, mode_t mode)
 {
     char           *pp;
     char           *sp;
     int             status;
     char           *copypath = strdup(path);
     status = 0;
     pp = copypath;
     while (status == 0 && (sp = strchr(pp, '/')) != 0)
     {
         if (sp != pp)
         {
             /* Neither root nor double slash in path */
             *sp = '\0';
             status = do_mkdir(copypath, mode);
             *sp = '/';
         }
         pp = sp + 1;
     }
     if (status == 0)
         status = do_mkdir(path, mode);
     free(copypath);
     return (status);
 }

void mkdir(const std::string &path) {
  if(mkpath(path.c_str(),0777)==-1) {
    std::cerr<<"Error: "<<strerror(errno)<<std::endl;
    exit(1);
  }
}

// void fileio(FILE* stream, const std::string &name, const std::string &option) {
//   if ( (stream = fopen(name.c_str(),option.c_str())) ==NULL ) {
//     std::cerr<<"Error: Cannot open file "<<name<<"!\n";
//     exit(1);
//   }
// }
            
//Find max / min====================================//
// template <typename T>
// int findpeak(T arr[], int m,compare cmp){
//   int mindex=0;
//   for (int i=1;i<m;i++) {
//     if (cmp(arr[i],arr[mindex]))
//       mindex=i;
//   }
//   return mindex;
// }


template <typename T>
int findmax(T arr[], int m) {
  int mindex=0;
  for (int i=1;i<m;i++) {
    if (arr[mindex]<arr[i])
      mindex=i;
  }
  return mindex;
}

template <typename T>
int findmin(T arr[], int m) {
  int mindex=0;
  for (int i=1;i<m;i++) {
    if (arr[mindex]>arr[i])
      mindex=i;
  }
  return mindex;
}

#endif
