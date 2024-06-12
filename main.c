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

#define dt 0.001
#define ROUTE_LENGTH 1000

#define EXPECTED_TIME 100

const int SECIONS_SLOPES[] = {0, 2, -1, 1, 3, -2, 1, 2, 1, -1};

inline double compute_a(double t, double v, double alpha, int F)
    __attribute__((always_inline));

inline double compute_a(double t, double v, double alpha, int F) {
  return (F - (cx * p * A * v) - m * g * sin(alpha * M_PI / 180)) / m;
}

// Runge kutta 4
int main(int argc, char *argv[]) {
  double t = 0.0;
  double v = 0.0;
  double s = 0.0;
  int F = 1000;

  // Time stuff
  clock_t start, end;
  start = clock();

  double acceleration_time =
      (double)(F - (cx * p * A * ROUTE_LENGTH / EXPECTED_TIME)) / m *
      ROUTE_LENGTH / EXPECTED_TIME;

  double k1, k2, k3, k4;
  int i = 0;

  while (true) {
    if (s > (i > 0 ? i * 100 : i + 100)) {
      i++;
    }

    if (t > acceleration_time) {
      F = 0;
    }

    if (t > 100) {
      if (fabs(s - ROUTE_LENGTH) < dt * 100) {
        if (s > ROUTE_LENGTH) {
          end = clock();

          printf("\n\n\n\n");
          printf("------------- Solved --------------\n\n");
          printf("A: %f, V: %f, S: %f\n", k1, v, s);
          printf("T: %f, AbsoluteError: %f, AccelerationTime: %f\n", t,
                 fabs(s - ROUTE_LENGTH), acceleration_time);
          printf("Total time: %fs\n", ((double)(end - start)) / CLOCKS_PER_SEC);

          break;
        }

      } else if (k1 < 0) {
        acceleration_time += dt;
        printf("A: %f, V: %f, S: %f\n", k1, v, s);
        printf("T: %f, AbsoluteError: %f, AccelerationTime: %fs\n", t,
               fabs(s - ROUTE_LENGTH), acceleration_time);
        printf("Total time: %fs\n", ((double)(end - start)) / CLOCKS_PER_SEC);
        t = 0;
        s = 0;
        v = 0;
        i = 0;
        F = 1000;
      } else if (s > ROUTE_LENGTH) {
        acceleration_time -= dt;
        printf("A: %f, V: %f, S: %f\n", k1, v, s);
        printf("T: %f, AbsoluteError: %f, AccelerationTime: %fs\n", t,
               fabs(s - ROUTE_LENGTH), acceleration_time);
        printf("Total time: %fs\n", ((double)(end - start)) / CLOCKS_PER_SEC);
        t = 0;
        s = 0;
        v = 0;
        i = 0;
        F = 1000;
      } else if (s < ROUTE_LENGTH) {
        acceleration_time += dt;
        printf("A: %f, V: %f, S: %f\n", k1, v, s);
        printf("T: %f, AbsoluteError: %f, AccelerationTime: %fs\n", t,
               fabs(s - ROUTE_LENGTH), acceleration_time);
        printf("Total time: %fs\n", ((double)(end - start)) / CLOCKS_PER_SEC);
        t = 0;
        s = 0;
        v = 0;
        i = 0;
        F = 1000;
      }
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
