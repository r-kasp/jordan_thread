#ifndef THREAD_H
#define THREAD_H

#include <pthread.h>

class ARGS
{
public:
  int id = 0;
  int size = 0;
  double *a = nullptr;
  double *b = nullptr;
  double *rev = nullptr; //обратный блок
  int ind_max = -1;
  double min_norm = 1e308;
  double norm = 0;
  int n = 0;
  int m = 0;
  int res = 0;
  bool is_sing = false;
  double cpu_time = 0;
};

#endif
