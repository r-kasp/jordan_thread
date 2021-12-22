#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <sys/resource.h>
#include <sys/time.h>
#include "matrix.h"
#include "thread.h"
#include <unistd.h>

#define ok 0

double get_cpu_time()
{
  struct rusage buf;
  getrusage(RUSAGE_THREAD, &buf);
  return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e+6;
}


double get_full_time()
{
  struct timeval buf;
  gettimeofday(&buf, 0);
  return buf.tv_sec + buf.tv_usec / 1.e+6;
}

void reduce_sum (int p)
{
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int t_in = 0;
  static int t_out = 0;
  if (p <= 1)
    return;
  pthread_mutex_lock(&m);
  t_in++;
  if (t_in >= p)
    {
      t_out = 0;
      pthread_cond_broadcast (&c_in);
    }
  else
    {
      while (t_in < p)
        pthread_cond_wait (&c_in, &m);
    }
  t_out++;
  if (t_out >= p)
    {
      t_in = 0;
      pthread_cond_broadcast (&c_out);
    }
  else
    {
      while (t_out < p)
        pthread_cond_wait (&c_out, &m);
    }
  pthread_mutex_unlock(&m);
}

void solve(ARGS *arg)
{
  arg->cpu_time = get_cpu_time();
  int pcnt = arg->size;
  int p = arg->id;
  int m = arg->m;
  int n = arg->n;
  double *a = arg->a;
  double *b = arg->b;
  
  int s = n / m;
	int ost = n % m;
	double eps = 1e-16;
	eps *= arg->norm;
  double *rev = arg->rev;

  double *buf = new double[m*m];
  double *buf_mm = new double[m*m];
	double *buf_mm2 = new double[m*m];
	double *buf_ostm = new double[m*ost];
	double *buf_ostm2 = new double[m*ost];
	double *buf_ostost = new double[ost*ost];
	double *buf_ost = new double[ost];
	double *buf_m = new double[m];
	double *buf_m2 = new double[m];
	int *r_m = new int[m];
	int *r_ost = new int[ost];
  for (int k = 0; k < s; k++)
	{  
	  int ind_max;
	  int start;
	  if (is_thread(k, p, pcnt))
	    start = k;
	  else
	    start = get_index(k, p, pcnt);
	  reduce_sum(pcnt);
	  for (int ind = start; ind < s; ind += pcnt)
		{
      int ret = get_reversed_matrix_for_block(a, n, m, m, ind, k, buf, eps, buf_mm, r_m);
      if (ret == SING_MATRIX)
        continue;
      double check = norma_matrix(buf, m, m);
      
      if (check < arg->min_norm)
      {
        arg->min_norm = check;
        arg->ind_max = ind;
      }
		}
		
		reduce_sum(pcnt);
	  
	  for (int i = 0; i < arg->size; i++)
    {
      ARGS * bufarg = arg - arg->id + i;
      if (bufarg->min_norm <= arg->min_norm)
        {
          arg->min_norm = bufarg->min_norm;
          arg->ind_max = bufarg->ind_max;
        }
    }
    
    reduce_sum(pcnt);
    
    bool is_sing = true;
    for (int i = -p; p + i < pcnt; i++)
	    {
	      if ((arg + i)->ind_max >= 0)
	        {
	          is_sing = false;
	        }
	    }
	    
	  if (is_sing)
	  {
      for (int i = -p; p + i < pcnt; i++)
	    {
	      (arg + i)->is_sing = true;
	    }
	  }
	  else
	  {
	    for (int i = -p; p + i < pcnt; i++)
	    {
	      (arg + i)->is_sing = false;
	    }
	  }
	  
	  //была синх
		
		if (arg->is_sing)
		{
		  clean_memory(buf, buf_mm, buf_mm2, buf_ostm, buf_ostm2, buf_ostost, buf_ost, buf_m, buf_m2, r_m, r_ost);
		  arg->cpu_time = get_cpu_time() - arg->cpu_time;
      reduce_sum(pcnt);
      return;
		}
		
		//была синх
		reduce_sum(pcnt);
		if (arg->ind_max == -1)
		{
		  for (int i = -p; p + i < pcnt; i++)
	    {
	      if ((arg + i)->ind_max != -1)
	      {
	        arg->ind_max = (arg + i)->ind_max;
	        break;
	      }
	    }
	  }
		
		ind_max = arg->ind_max;
		//if (k > (n / pcnt) * pcnt && p < k % pcnt)
		//  continue;
		
		if (is_k_for_p(k, p, pcnt))
		{
	    get_reversed_matrix_for_block(a, n, m, m, ind_max, k, rev, eps, buf_mm, r_m);
	    //set_unitary_block(a, ind_max, k, n, m, m, m);
	  }
	  reduce_sum(pcnt); //ЖБ
	
		multiply_blocks_on_str(p, pcnt, a, rev, ind_max, k+1, n, m, buf_mm, buf_mm2, buf_ostm, buf_ostm2);
		
		if (p == pcnt - 1) // пусть вектора умножает первый поток
		  multiply_block_on_vector(rev, b, ind_max, m, m, buf_m);
		
		reduce_sum(pcnt); //ЖБ
		
		for (int ind = 0; ind < s; ind++)
		{
		  if (ind == ind_max)
				continue;
			get_block(a, ind, k, buf, n, m, m, m);
			if (norma_matrix(buf, m, m) < eps) //если на этой строчке итак элемент в k-ом столбце нулевой
			  continue;
			subtract_with_multiply(p, pcnt, a, ind, ind_max, k+1, buf, n, m, buf_mm, buf_mm2);
			if (p == 0)
			  subtract_with_multiply_vector(b, ind, ind_max, buf, n, m, buf_m, buf_m2, buf_ost);
		}
		if (ost != 0) 
		{
		  get_block(a, s, k, buf, n, m, ost, m);
			if (norma_matrix(buf, ost, m) < eps) //если на этой строчке итак элемент в k-ом столбце нулевой
			  continue;
			subtract_with_multiply(p, pcnt, a, s, ind_max, k+1, buf, n, m, buf_mm, buf_mm2);
			if (p == 0)
			  subtract_with_multiply_vector(b, s, ind_max, buf, n, m, buf_m, buf_m2, buf_ost);
		}
		//reduce_sum(pcnt);
		//была синх
		//change_str_blocks(p, pcnt, a, b, k, ind_max, s, m, n);
	  //была синх
	  
	  for (int i = 0; i < arg->size; i++)
    {
      ARGS * bufarg = arg - arg->id + i;
      bufarg->ind_max = -1;
      bufarg->min_norm = 1e308;
    }
    
    reduce_sum(pcnt);
    change_str_blocks(p, pcnt, a, b, k, ind_max, s, m, n, buf_mm, buf_mm2, buf_m, buf_m2);
	}
	if (ost != 0)
	{
	  //поделить эту строчку
	  if ((s + 1) % pcnt == p)
	  {
	    get_block(a, s, s, buf, n, m, ost, ost);
	    if (norma_matrix(buf, ost, ost) < eps)
	    {
	      for (int i = -p; p + i < pcnt; i++)
	      {
	        (arg + i)->is_sing = true;
	      }
	    }
	  }
	  reduce_sum(pcnt);
	  if (arg->is_sing)
		{
		  clean_memory(buf, buf_mm, buf_mm2, buf_ostm, buf_ostm2, buf_ostost, buf_ost, buf_m, buf_m2, r_m, r_ost);
      reduce_sum(pcnt);
      return;
		}
		if ((s + 1) % pcnt == p)
		{
      get_reversed_matrix_for_block(a, n, m, ost, s, s, rev, eps, buf_ostost, r_ost);
      set_unitary_block(a, s, s, n, m, ost, ost);
      multiply_block_on_vector(rev, b, s, m, ost, buf_ost);
      nill_last_column(a, b, n, m, buf_ostm, buf_m, buf_ost);
    }
	}
	reduce_sum(pcnt);
  clean_memory(buf, buf_mm, buf_mm2, buf_ostm, buf_ostm2, buf_ostost, buf_ost, buf_m, buf_m2, r_m, r_ost);
  arg->cpu_time = get_cpu_time() - arg->cpu_time;
  reduce_sum(pcnt);
  return;
}

static void* thread_func (void *args)
{
  ARGS *arg = (ARGS*) args;
  solve(arg);
  return ok;
}

int main(int argc, char *argv[])
{
    int n, m, p, r, s;
    char * filename = nullptr;
    if (!((argc == 6 || argc == 7) && sscanf(argv[1], "%d", &n) == 1 && 
        sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &p) == 1 && sscanf(argv[4], "%d", &r) == 1 && sscanf(argv[5], "%d", &s) == 1))
    {
        printf("Usage %s n m p r s (file) \n", argv[0]);
        return 1;
    }
    
    if (argc == 7)
        filename = argv[6];
    //Исходная матрица, которую не будем преобразовывать, чтобы сверить ответ
    double *a = new double[n*n];
    if (!a)
    {
        printf("Not enough memory\n");
        return 2;
    }
    //Вектор b
    double *b = new double[n];
    if (!b)
    {
        printf("Not enough memory\n");
        delete [] a;
        return 3;
    }
    //Копия вектора b
    double *c = new double[n];
    if (!c)
    {
        printf("Not enough memory\n");
        delete [] a;
        delete [] b;
        return 4;
    }
    
    if (filename)
    {
        int ret = read_matrix(a,n,filename);
        if (ret != SUCCESS)
        {
            switch (ret)
            {
                case ERROR_OPEN:
                    printf("Cannot open %s\n", filename);
                    break;
                case ERROR_READ:
                    printf("Cannot read %s\n", filename);
                    break;
                default:
                    printf("Unknown error %d in file %s\n", ret, filename);
            }
            delete [] a;
            delete [] b;
            delete [] c;
            return 5;
        }
    }
    else
        init_matrix(a,n,s);
    init_vector(b, a, n);
        
    print_matrix(a,n,n,r);
    
    copy_vector(c, b, n);
    
    int size = p;
    if (size > n / m + n % m)
      size = n / m + n % m;
    pthread_t *threads = new pthread_t[size - 1];
    ARGS *args = new ARGS[size];
    double *buf = new double[m*m];
    double norm = norma_matrix(a, n, n);
    for (int i = 0; i < size; i++)
      {
        args[i].id = i;
        args[i].a = a;
        args[i].b = b;
        args[i].size = size;
        args[i].ind_max = 0;
        args[i].n = n;
        args[i].m = m;
        args[i].rev = buf;
        args[i].norm = norm;
        args[i].cpu_time = 0;
      }
      
    double elapsed = get_full_time();
    for (int i = 0; i < size - 1; i++)
      {
        int ret = pthread_create (&threads[i], nullptr, thread_func, (void *)&args[i]);
        if (ret != ok)
          {
            printf("ERROR: CAN'T CREATE THREAD\n");
            delete [] a;
            delete [] b;
            delete [] c;
            delete [] buf;
            delete [] threads;
            delete [] args;
            return ERROR;
          }
      }
    ARGS *arg = &args[size - 1];
    solve(arg);    
    elapsed = get_full_time() - elapsed;
    
    int is_sing = false;
    for (int i = 0; i < size; i++)
    {
      if (args[i].is_sing)
        is_sing = true; 
    }
    if (filename)
    {
        int ret = read_matrix(a,n,filename);
        if (ret != SUCCESS)
        {
            switch (ret)
            {
                case ERROR_OPEN:
                    printf("Cannot open %s\n", filename);
                    break;
                case ERROR_READ:
                    printf("Cannot read %s\n", filename);
                    break;
                default:
                    printf("Unknown error %d in file %s\n", ret, filename);
            }
            delete [] a;
            delete [] b;
            delete [] c;
            delete [] buf;
            delete [] args;
            delete [] threads;
            return 7;
        }
    }
    else
        init_matrix(a,n,s);
    double *rsd = new double[n];
    if (!rsd)
    {
    	  printf("Not enough memory\n");
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] buf;
        delete [] args;
        delete [] threads;
        return 6;
    }
    multiply_matrix_and_vector(a, b, rsd, n, n);
    subtract_vectors(rsd, c, n);
    double residual = norma_vector(rsd, n) / norma_vector(c, n);
    printf("Result : \n");
    if (is_sing)
    {
        printf("Can't find solution. The matrix is singular\n");
        printf("Elapsed = %.2f\n", elapsed);
    }
    else
    {
        print_vector(b,n,r);
        printf("%s : residual = %e elapsed = %.2f s = %d n = %d m = %d p = %d\n", argv[0], residual, elapsed, s, n, m, p);
    }
    printf("CPU time for threads : \n");
    for (int i = 0; i < size; i++)
    {
      printf("thread : %d - %.2f\n", i, args[i].cpu_time);
    }
    delete [] a;
    delete [] b;
    delete [] rsd;
    delete [] c;
    delete [] buf;
    delete [] args;
    delete [] threads;
    return 0;
}
