#ifndef Functions_h
#define Functions_h

double condition(double x);

void Jacobian (double **JJ, double dTC, double dmuBC, double angle1, double angle2);

void funcv(int n,double x[],double f[]);

double G (double R, double Theta);

double dGdmuBConT (double R, double Theta);

double dGdTConmuB (double R, double Theta);

double d2GdmuB2ConT (double R, double Theta);

double d2GdT2ConmuB (double R, double Theta);

double d2GdTdmuB (double R, double Theta);

double d4GdmuB4ConT (double R, double Theta);

double GausFunc (double x, double y, double sigmax, double sigmay);

void GausFilter2D (double **data, double sigmax, double sigmay, int lx, int ly);



#endif /* Functions_h */