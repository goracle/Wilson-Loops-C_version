double** minors(double **M, int nn, int row, int column, double** minorr);

double* times(double* z, double* y, double* result);

double* getdet(double** M, int nn, double* temp_fin);

double getnorm(double** M, int col);

double getnorm2(double* a, double *result);

double* star(double* x, double* result);

double*** gen_rand_matrix();

int print_mat(double** M);

int testarg(int arg);

double**** update(double**** lattice, double*** container);

double*** findstap(int mu, int mua, int nu, int* x, int coor, int j, double** res1, double** res2, double** newM, double** newM2, double**** lattice, double ***staples, int updown);

double calculate_S(double ****lattice);

double plaq(double ****lattice, int coor, int nu, int mu, double result);

double **Times(double** M, double** N, double** res);

double trace_img(double **M);

double trace_real(double **M);

double** dagger(double** M, double** newM);

double**** initialize_lat();
  
int* findcoord(int coor, int* x);

