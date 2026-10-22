#define _POSIX_C_SOURCE 200809L   /* needed for clock_gettime under -std=c99 */
#include <stdio.h>
#include <math.h>
#include <time.h>

#define K 4.0
#define C 0.4
#define T_END 10.0

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void exact(double t, double *x, double *v) {
    double gamma = C / 2.0;
    double wd = sqrt(K - gamma * gamma);
    double e = exp(-gamma * t);
    *x = e * (cos(wd * t) + (gamma / wd) * sin(wd * t));
    *v = -e * (K / wd) * sin(wd * t);
}

static void deriv(double x, double v, double *dx, double *dv) {
    *dx = v;
    *dv = -K * x - C * v;
}

static void euler(long n, double h, double *x, double *v) {
    double xi = 1.0, vi = 0.0;
    for (long i = 0; i < n; i++) {
        double a = -K * xi - C * vi;
        xi += h * vi;
        vi += h * a;
    }
    *x = xi;
    *v = vi;
}

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

static volatile double sink;   /* keeps the compiler from deleting the work */

int main(int argc, char **argv) {
    const char *out = argc > 1 ? argv[1] : "timing.csv";
    double xe, ve;
    exact(T_END, &xe, &ve);

    FILE *f = fopen(out, "w");
    if (!f) { perror("fopen"); return 1; }
    fprintf(f, "method,h,error,seconds\n");

    for (int m = 0; m < 2; m++) {
        const char *name = m == 0 ? "euler" : "rk4";
        for (int i = 0; i < 15; i++) {
            long n = 10L << i;
            volatile double hv = T_END / (double)n;   /* volatile: stops the compiler hoisting the call */
            long reps = 2000000L / n;
            if (reps < 1) reps = 1;

            double x = 0.0, v = 0.0, best = 1e30;
            for (int trial = 0; trial < 5; trial++) {   /* keep the fastest of 5 trials */
                double t0 = now();
                for (long r = 0; r < reps; r++) {
                    if (m == 0) euler(n, hv, &x, &v);
                    else        rk4(n, hv, &x, &v);
                    sink = x;
                }
                double dt = (now() - t0) / (double)reps;
                if (dt < best) best = dt;
            }
            double err = sqrt((x - xe) * (x - xe) + (v - ve) * (v - ve));
            fprintf(f, "%s,%.10e,%.10e,%.10e\n", name, (double)hv, err, best);
        }
    }
    fclose(f);
    printf("wrote %s\n", out);
    return 0;
}