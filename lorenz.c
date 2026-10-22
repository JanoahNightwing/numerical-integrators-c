#include <stdio.h>
#include <math.h>

#define SIGMA 10.0
#define RHO   28.0
#define BETA  (8.0 / 3.0)
#define H     0.01
#define T_END 50.0

typedef struct { double x, y, z; } vec;

static vec f(vec s) {
    vec d = { SIGMA * (s.y - s.x), s.x * (RHO - s.z) - s.y, s.x * s.y - BETA * s.z };
    return d;
}

static vec axpy(vec s, double a, vec d) {
    vec r = { s.x + a * d.x, s.y + a * d.y, s.z + a * d.z };
    return r;
}

static vec rk4(vec s, double h) {
    vec k1 = f(s);
    vec k2 = f(axpy(s, 0.5 * h, k1));
    vec k3 = f(axpy(s, 0.5 * h, k2));
    vec k4 = f(axpy(s, h, k3));
    vec r = { s.x + h / 6.0 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x),
              s.y + h / 6.0 * (k1.y + 2 * k2.y + 2 * k3.y + k4.y),
              s.z + h / 6.0 * (k1.z + 2 * k2.z + 2 * k3.z + k4.z) };
    return r;
}

static double dist(vec a, vec b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

int main(void) {
    vec a = {1.0, 1.0, 1.0};
    for (int i = 0; i < 2000; i++) a = rk4(a, H);   /* burn-in: 20 time units to land on the attractor */
    vec b = a, c = a;
    b.x += 1e-10;
    long n = (long)llround(T_END / H);
    double t_b = -1.0, t_c = -1.0;

    FILE *fp = fopen("lorenz.csv", "w");
    if (!fp) { perror("fopen"); return 1; }
    fprintf(fp, "t,sep_perturb,sep_stepsize\n");

    for (long i = 0; i <= n; i++) {
        double sb = dist(a, b), sc = dist(a, c);
        if (t_b < 0.0 && sb > 1.0) t_b = i * H;
        if (t_c < 0.0 && sc > 1.0) t_c = i * H;
        if (i % 10 == 0) fprintf(fp, "%.4f,%.6e,%.6e\n", i * H, sb, sc);
        if (i == n) break;
        a = rk4(a, H);
        b = rk4(b, H);
        c = rk4(rk4(c, H / 2.0), H / 2.0);   /* two half-steps = one step of size H */
    }
    printf("1e-10 perturbation: separation first exceeds 1.0 at t = %.2f\n", t_b);
    printf("step size h vs h/2:  separation first exceeds 1.0 at t = %.2f\n", t_c);
    fclose(fp);
    return 0;
}