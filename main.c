#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#define m 1000
#define cx 0.25
#define A 2
#define p 1.1225
#define g 9.81

#define dt 0.00001
#define ROUTE_LENGTH 1000
#define V_MULTIPLIER 100

#define EXPECTED_TIME 100

const int SECIONS_SLOPES[] = {0, 2, -1, 1, 3, -2, 1, 2, 1, -1};

inline double compute_a(double t, double v, double alpha, int F)
    __attribute__((always_inline));

inline double compute_a(double t, double v, double alpha, int F) {
  return (F - (cx * p * A * v) - m * g * sin(alpha * M_PI / 180)) / m;
}

// Runke kutta 4
int main(int argc, char *argv[]) {
  double t = 0.0;
  double v = 0.0;
  double s = 0.0;
  int F = 1000;

  // Time stuff
  clock_t start, end;
  start = clock();

  float maybe_v = (ROUTE_LENGTH / 100.0f) - 0.909;

  double k1, k2, k3, k4;
  int i = 0;

  while (true) {
    if (s > (i > 0 ? i * 100 : i + 100)) {
      i++;
    } else if (i > sizeof(SECIONS_SLOPES) / sizeof(int)) {
      double absolute_error = fabs(t - EXPECTED_TIME);
      //       end = clock();
      //       printf("A: %f, V: %f, S: %f\n", k1, v, s);
      //       printf("T: %f, AbsoluteError: %f, MaxV: %f\n", t, absolute_error,
      //              maybe_v);
      //       printf("Time sint pent: %fs\n", ((double)(end - start)) /
      //       CLOCKS_PER_SEC);

      if (absolute_error < dt * V_MULTIPLIER * 10) {
        end = clock();
        printf("A: %f, V: %f, S: %f\n", k1, v, s);
        printf("T: %f, AbsoluteError: %f, MaxV: %f\n", t, absolute_error,
               maybe_v);
        printf("Total time: %fs\n", ((double)(end - start)) / CLOCKS_PER_SEC);

        break;
      } else {
        if (t < EXPECTED_TIME) {
          maybe_v -= dt * V_MULTIPLIER;
        } else if (t > EXPECTED_TIME) {
          maybe_v += dt * V_MULTIPLIER;
        }

        t = 0;
        s = 0;
        v = 0;
        i = 0;
      }
    }

    if (v > maybe_v) {
      F = 0;
    } else {
      F = 1000;
    }

    k1 = compute_a(t, v, SECIONS_SLOPES[i], F);
    k2 = compute_a(t + dt / 2, v, SECIONS_SLOPES[i], F);
    k3 = k2;
    k4 = compute_a(t + dt, v, SECIONS_SLOPES[i], F);

    v += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    s += v * dt;
    t += dt;
  }

  return EXIT_SUCCESS;
}
