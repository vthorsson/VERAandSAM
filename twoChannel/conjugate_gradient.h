void frprmn(double*, int, double, int *, double *,
            double (*)(double []), void (*)(double [], double []));         
void linmin(double p[], double xi[], int n, double *,
                double (*func)(double []) );
double brent(double ax, double bx, double cx,
                double (*f)(double), double tol, double *xmin);
double f1dim(double x);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
                double *fc, double (*func)(double));            
