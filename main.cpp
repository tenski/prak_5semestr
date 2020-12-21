#include <iostream>
#include <cmath>

double estimationError(double *y2, double *y, int N) {
    double res;

    N = N / 2;

    for (int i = 0; i < N + 1; i++) {
        if (res < abs(y[2 * i] - y2[i]) / 3) {
            res = abs(y[2 * i] - y2[i]) / 3;
        }
    }

    return res;
}


void realSolv(double *u, int N, double x0) {
    double * x;

    double h = 1.0 / N;
    
    double k = x0 + 1;
    double q = exp(x0);
    double f = exp(-x0*x0); 

    double C = f / (q * ((sqrt(k * q) - 1) * exp(-sqrt(q / k)) - (sqrt(k * q) + 1) * exp(sqrt(q / k))));
    x = new double [N + 1];

    for (int i = 0; i < N + 1; i++) {
        x[i] = i * h;
    }

    for (int i = 0; i < N + 1; i++) {
        u[i] = C * (exp(sqrt(q / k) * x[i]) + exp(-sqrt(q / k) * x[i])) + f / q;
    }

    delete [] x;
}


void sweepMethod(double *y, int N, double x0) {
    double *k;
    double *q;
    double *fi;
    double *x;
    
    double *alf;
    double *bet;
    double *a;
    double *b;
    double *c;

    double h = 1.0 / N;

    k = new double [N + 1];
    q = new double [N + 1];
    fi = new double [N + 1];
    x = new double [N + 1];

    alf = new double [N + 2];
    bet = new double [N + 2];
    a = new double [N + 2];
    b = new double [N + 2];
    c = new double [N + 2];

    for (int i = 0; i < N + 1; i++) {
        x[i] = i * h;

        if (x0 == 0.0) {
            k[i] = x[i] + 1;
            q[i] = exp(x[i]);
            fi[i] = exp(-x[i]*x[i]);
        } else {
            k[i] = x0 + 1;
            q[i] = exp(x0);
            fi[i] = exp(-x0*x0);
        }
    }

    for (int i = 1; i < N; i++) {
        a[i] = 1 / (2 * h * h) * (k[i - 1] + k[i]);
        b[i] = 1 / (2 * h * h) * (k[i] + k[i + 1]);
        c[i] = q[i] + 1 / (2 * h * h) * (k[i - 1] + 2 * k[i] + k[i + 1]);
    }

    alf[1] = (k[0] + k[1]) / (k[0] + k[1] + h * h * q[0]);
    bet[1] = h * h * fi[0] / (k[0] + k[1] + h * h * q[0]);

    for (int i = 1; i < N; i++) {
        alf[i + 1] = b[i] / (c[i] - a[i] * alf[i]);
        bet[i + 1] = (fi[i] + a[i] * bet[i]) / (c[i] - a[i] * alf[i]); 
    }

    double mu1 = 1 + 0.5 * h * q[N] + (k[N] + k[N - 1]) / (2 * h);
    double hi2 = (k[N] + k[N - 1]) / (2 * h * mu1);

    y[N] = (0.5 * h * fi[N] / mu1 + hi2 * bet[N]) / (1 - alf[N] * hi2);

    for (int i = N - 1; i >= 0; i--) {
        y[i] = alf[i + 1] * y[i + 1] + bet[i + 1];
    }

    delete [] k;
    delete [] q;
    delete [] fi;
    delete [] x;

    delete [] alf;
    delete [] bet;
    delete [] a;
    delete [] b;
    delete [] c;
}


void modelProblem(int N) {
    double *u;
    double *y;
    double *y2;

    int j = 0;

    double eps = 0.01, x0 = 0.5;

    while (N < 100000) {
        u = new double [N + 1];
        y = new double [N + 1];

        realSolv(u, N, x0);
        sweepMethod(y, N, x0);

        if (j == 1) {
            if (estimationError(y2, y, N) < eps) {
                for (int i = 0; i < N + 1; i += 1) {
                    std::cout << "x[" << i << "] = " << i * 1.0 / N << " , u[" << i << "] = " << u[i] << " , y[" << i << "] = " << y[i] << " , delta = " << abs(u[i] - y[i]) << std::endl;
                    /*std::cout << i << "&" << i * 1.0 / N << "&" << u[i] << "&" << y[i] << "&" << abs(u[i] - y[i]) << "\\\\" << std::endl; */
                }

                delete [] u;
                delete [] y;
                delete [] y2;

                break;
            } else {
                delete [] y2;

                y2 = new double [N + 1];

                for (int i = 0; i < N + 1; i++) {
                    y2[i] = y[i];
                }

                N = 2 * N;
            }
        } else {
            y2 = new double [N + 1];

            for (int i = 0; i < N + 1; i++) {
                y2[i] = y[i];
            }

            N = 2 * N;
        }

        delete [] u;
        delete [] y;

        j = 1;
    }
}


void generalProblem(int N) {
    double *u;
    double *y;
    double *y2;

    int j = 0;

    double eps = 0.01, x0 = 0.0;

    while (N < 100000) {
        y = new double [N + 1];

        sweepMethod(y, N, x0);

        if (j == 1) {
            if (estimationError(y2, y, N) < eps) {
                for (int i = 0; i < N + 1; i += 1) {
                    std::cout << "x[" << i << "] = " << i * 1.0 / N << " , y[" << i << "] = " << y[i] << std::endl;
                    /*std::cout << i << "&" << i * 1.0 / N << "&" << y[i] << "\\\\" << std::endl;*/
                }

                delete [] y;
                delete [] y2;

                break;
            } else {
                delete [] y2;

                y2 = new double [N + 1];

                for (int i = 0; i < N + 1; i++) {
                    y2[i] = y[i];
                }

                N = 2 * N;
            }
        } else {
            y2 = new double [N + 1];

            for (int i = 0; i < N + 1; i++) {
                y2[i] = y[i];
            }

            N = 2 * N;
        }

        delete [] y;

        j = 1;
    }
}

int main() {
    int N, model;

    std::cout << "Hello!" << std::endl;
    std::cout << "Please, choose which problem we will solve: model - 0, general - 1.\nEnter a number: ";
    std::cin  >> model;
    std::cout << "Please, enter N: ";
    std::cin  >> N;

    if (model == 0) {
        modelProblem(N);
    } else if (model == 1) {
        generalProblem(N);
    } else {
        std::cout << "When entering the first digit, choose between 0 and 1!\n Please, restart program!" << std::endl;
    }
    
    return 0;
}
