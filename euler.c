#include <stdio.h>
#include <math.h>

#define K 4.0       /* spring constant */
#define C 0.4       /* damping */
#define T_END 10.0  /* simulate from t = 0 to t = 10 */

/* Exact solution for x(0) = 1, v(0) = 0 (underdamped case). */
static void exact(double t, double *x, double *v) {
    double gamma = C / 2.0;
    double wd = sqrt(K - gamma * gamma);
    double e = exp(-gamma * t);
    *x = e * (cos(wd * t) + (gamma / wd) * sin(wd * t));
    *v = -e * (K / wd) * sin(wd * t);
}

/* The system as first-order equations: dx/dt = v, dv/dt = -K x - C v. */
static void deriv(double x, double v, double *dx, double *dv) {
    *dx = v;
    *dv = -K * x - C * v;
}

/* Forward Euler: n steps of size h, starting from x = 1, v = 0. */
static void euler(long n, double h, double *x, double *v) {
    double xi = 1.0, vi = 0.0;
    for (long i = 0; i < n; i++) {
        double a = -K * xi - C * vi;  /* acceleration from the OLD state */
        xi += h * vi;
        vi += h * a;
    }
    *x = xi;
    *v = vi;
}

/* Classic 4th-order Runge-Kutta: n steps of size h, same starting point. */
static void rk4(long n, double h, double *x, double *v) {
    double xi = 1.0, vi = 0.0;
    for (long i = 0; i < n; i++) {
        double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;
        deriv(xi, vi, &k1x, &k1v);
        deriv(xi + 0.5 * h * k1x, vi + 0.5 * h * k1v, &k2x, &k2v);
        deriv(xi + 0.5 * h * k2x, vi + 0.5 * h * k2v, &k3x, &k3v);
        deriv(xi + h * k3x, vi + h * k3v, &k4x, &k4v);
        xi += h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
        vi += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
    }
    *x = xi;
    *v = vi;
}

int main(void) {
    double x0, v0;
    exact(0.0, &x0, &v0);
    printf("sanity check: exact(0) = (%.3f, %.3f), should be (1.000, 0.000)\n\n", x0, v0);

    FILE *f = fopen("convergence.csv", "w");
    if (!f) { perror("fopen"); return 1; }
    fprintf(f, "h,euler_error,rk4_error\n");

    double xe, ve;
    exact(T_END, &xe, &ve);

    printf("%11s %13s %7s %13s %7s\n", "h", "euler err", "ratio", "rk4 err", "ratio");
    double pe = 0.0, pr = 0.0;
    for (int i = 0; i < 15; i++) {
        long n = 10L << i;
        double h = T_END / (double)n;
        double x, v;

        euler(n, h, &x, &v);
        double ee = sqrt((x - xe) * (x - xe) + (v - ve) * (v - ve));
        rk4(n, h, &x, &v);
        double er = sqrt((x - xe) * (x - xe) + (v - ve) * (v - ve));

        fprintf(f, "%.10e,%.10e,%.10e\n", h, ee, er);
        if (i == 0) {
            printf("%11.6f %13.4e %7s %13.4e %7s\n", h, ee, "-", er, "-");
        } else {
            printf("%11.6f %13.4e %7.2f %13.4e %7.2f\n", h, ee, pe / ee,
                   er, er > 0.0 ? pr / er : 0.0);
        }
        pe = ee;
        pr = er;
    }
    fclose(f);
    return 0;
}