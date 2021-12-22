#ifndef MATRIX_H
#define MATRIX_H

int read_matrix(double *a, int n, char *name);
void init_matrix(double *a, int n, int k);
void init_vector(double *b, double *a, int n);
void print_matrix(double *a, int m, int n, int r);
void copy_matrix(double *x, double *a, int n);
void copy_vector(double *c, double *b, int n);
double f(int s, int n, int i, int j);
int solve(int n, double *a, double *x, double *b, double *c);
double norma_vector(double *a, int n);
double norma_matrix(double *a, int m, int n);
void multiply_matrix_and_vector(double *a, double *b, double *c, int m, int n);
void multiply_blocks_m(double *a, double *b, double *c, int m);
void subtract_vectors(double *b, double *c, int n);
void print_vector(double *b, int n, int r);
int get_reversed_matrix_gauss(int n, double *a, double *x, int *r, double eps);


void fill_arranged_items(int *array, int size);
void set_block(double *a, int i, int j, double *b, int n, int m, int height, int width);
void get_block(double *a, int i, int j, double *b, int n, int m, int height, int width);
void get_vector(double *b, double *res, int ind, int m, int size);
void set_vector(double *b, double *buf, int ind, int m, int size);
void set_unitary_block(double *a, int i, int j, int n, int m, int width, int height);
void set_unitary_matrix(double *a, int m);
void set_null_block(double *a, int i, int j, int n, int m, int height, int width);
int get_reversed_matrix_for_block(double *a, int n, int m, int size, int i, int j, double *res, double eps, double *buf, int * r);
void change_str_blocks(int p, int pcnt, double *a, double *b, int k, int ind_max, int s, int m, int n, double *buf_mm, double *buf_mm2, double *buf_m, double *buf_m2);
void arrange_vector_items(double *b, int n, int m, int *numbers);
void multiply_block_on_vector(double *buf, double *b, int str_ind, int m, int block_sz, double *b_old);
void multiply_blocks_on_str(int p, int psize, double *a, double *mult, int i, int j, int n, int m, double *buf, double *res, double *buf2, double *res2);
void subtract_blocks(double *a, double *b, int n, int m);
void subtract_with_multiply(int p, int psize, double *a, int from_ind, int ind, int start_stl, double *mult, int n, int m, double *buf, double *res);
void subtract_with_multiply_vector(double *b, int from_ind, int ind, double *mult, int n, int m, double *buf, double *res, double *res2);
void fill_false(bool *a, int n);
int find_block_with_the_min_norm(double *a, int n, int m, int stl, int * numbers, double eps, double *buf, double *buf2, int *r);
void nill_last_column(double *a, double *b, int n, int m, double *buf, double *res, double *b_last);
bool is_k_for_p(int k, int p, int psize);
void clean_memory(double *buf, double *buf_mm, double *buf_mm2, double *buf_ostm, double *buf_ostm2, double *buf_ostost, double *buf_ost, double *buf_m, double *buf_m2, int *r_m, int *r_ost);
bool is_thread(int j, int p, int psize);
int get_index(int j, int p, int psize);

enum RETURN_CODES
{
    SUCCESS,
    ERROR_READ,
    ERROR_OPEN,
    SING_MATRIX,
    ERROR
};

#endif
