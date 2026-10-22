#include <stdio.h>
#include <math.h>

#define K 4.0
#define H 0.05
#define T_END 500.0

static double energy(double x, double v) { return 0.5 * v * v + 0.5 * K * x * x; }

static void deriv(double x, double v, double *dx, double *dv) {
    *dx = v;
    *dv = -K * x;
}

static void euler_step(double *x, double *v, double h) {
    double a = -K * (*x);
    *x += h * (*v);
    *v += h * a;
}

static void rk4_step(double *x, double *v, double h) {
    double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;
    deriv(*x, *v, &k1x, &k1v);
    deriv(*x + 0.5 * h * k1x, *v + 0.5 * h * k1v, &k2x, &k2v);
    deriv(*x + 0.5 * h * k2x, *v + 0.5 * h * k2v, &k3x, &k3v);
    deriv(*x + h * k3x, *v + h * k3v, &k4x, &k4v);
    *x += h / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    *v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
}

/* Velocity Verlet: a symplectic method. */
static void verlet_step(double *x, double *v, double h) {
    double a = -K * (*x);
    *x += h * (*v) + 0.5 * h * h * a;
    double a_new = -K * (*x);
    *v += 0.5 * h * (a + a_new);
}

int main(void) {
    long n = (long)llround(T_END / H);
    double xe = 1.0, ve = 0.0, xr = 1.0, vr = 0.0, xv = 1.0, vv = 0.0;
    double E0 = energy(1.0, 0.0);
    double max_rk = 0.0, max_vl = 0.0;

    FILE *f = fopen("energy.csv", "w");
    if (!f) { perror("fopen"); return 1; }
    fprintf(f, "t,euler,rk4,verlet\n");

    for (long i = 0; i <= n; i++) {
        double Ee = energy(xe, ve) / E0;
        double Er = energy(xr, vr) / E0;
        double Ev = energy(xv, vv) / E0;
        if (fabs(Er - 1.0) > max_rk) max_rk = fabs(Er - 1.0);
        if (fabs(Ev - 1.0) > max_vl) max_vl = fabs(Ev - 1.0);
        if (i % 10 == 0) fprintf(f, "%.6f,%.12e,%.12e,%.12e\n", i * H, Ee, Er, Ev);
        if (i == n) {
            printf("final E/E0:   euler %.4e   rk4 %.10f   verlet %.10f\n", Ee, Er, Ev);
            break;
        }
        euler_step(&xe, &ve, H);
        rk4_step(&xr, &vr, H);
        verlet_step(&xv, &vv, H);
    }
    printf("max |E/E0 - 1|:   rk4 %.4e   verlet %.4e\n", max_rk, max_vl);
    fclose(f);
    return 0;
}