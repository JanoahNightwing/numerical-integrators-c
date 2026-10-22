#include <stdio.h>
#include <math.h>

#define K 4.0
#define H 0.05
#define T_END 100000.0
#define BOUND (K * H * H / 4.0)   /* bound velocity Verlet should stay within */

static double energy(double x, double v) { return 0.5 * v * v + 0.5 * K * x * x; }

static void deriv(double x, double v, double *dx, double *dv) {
    *dx = v;
    *dv = -K * x;
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

static void verlet_step(double *x, double *v, double h) {
    double a = -K * (*x);
    *x += h * (*v) + 0.5 * h * h * a;
    double a_new = -K * (*x);
    *v += 0.5 * h * (a + a_new);
}

int main(void) {
    long n = (long)llround(T_END / H);
    double xr = 1.0, vr = 0.0, xv = 1.0, vv = 0.0;
    double E0 = energy(1.0, 0.0);
    double max_rk = 0.0, max_vl = 0.0;
    double t_cross = -1.0;

    FILE *f = fopen("energy_long.csv", "w");
    if (!f) { perror("fopen"); return 1; }
    fprintf(f, "t,rk4,verlet\n");

    for (long i = 0; i <= n; i++) {
        double Er = energy(xr, vr) / E0;
        double Ev = energy(xv, vv) / E0;
        if (fabs(Er - 1.0) > max_rk) max_rk = fabs(Er - 1.0);
        if (fabs(Ev - 1.0) > max_vl) max_vl = fabs(Ev - 1.0);
        if (t_cross < 0.0 && fabs(Er - 1.0) > BOUND) t_cross = i * H;
        if (i % 200 == 0) fprintf(f, "%.6f,%.12e,%.12e\n", i * H, Er, Ev);
        if (i == n) {
            printf("final E/E0:   rk4 %.6f   verlet %.6f\n", Er, Ev);
            break;
        }
        rk4_step(&xr, &vr, H);
        verlet_step(&xv, &vv, H);
    }
    printf("max |E/E0 - 1|:   rk4 %.4e   verlet %.4e\n", max_rk, max_vl);
    if (t_cross >= 0.0)
        printf("RK4 energy error first exceeds %.4f at t = %.1f\n", BOUND, t_cross);
    else
        printf("RK4 energy error never exceeded %.4f\n", BOUND);
    fclose(f);
    return 0;
}