#define _GNU_SOURCE
#include "mersenne_twister.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000            // 粒子数 
#define DT 0.005          // 時間刻み
#define SOFTENING 0.01    // ソフトニングパラメータ
#define TOTAL_TIME 10.0   // 時間

#define G_CONST 1.0
#define M_TOTAL 1.0
#define MASS (M_TOTAL / N)

typedef struct { double x, y, z; } Vec3;
typedef struct { Vec3 pos, vel; } Particle;

double ran_sym(void) {
    return genrand_real3() * 2.0 - 1.0;
}

// Box-Muller法で標準正規分布を与える
double ran_gauss(void) {
    double u1 = genrand_real3();
    double u2 = genrand_real3();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// qsort用の比較関数（半質量半径の計算に使用）
int compare_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

// 加速度の計算とポテンシャルエネルギー U の評価
double compute_acc_and_pot(const Particle p[N], Vec3 acc[N]) {
    for (int i = 0; i < N; i++) {
        acc[i] = (Vec3){0.0, 0.0, 0.0};
    }

    double U = 0.0;
    const double eps2 = SOFTENING * SOFTENING;
    const double gm = G_CONST * MASS;
    const double gm2 = G_CONST * MASS * MASS;

    for (int j = 0; j < N - 1; j++) {
        for (int k = j + 1; k < N; k++) {
            double dx = p[k].pos.x - p[j].pos.x;
            double dy = p[k].pos.y - p[j].pos.y;
            double dz = p[k].pos.z - p[j].pos.z;

            double r2 = dx * dx + dy * dy + dz * dz + eps2;
            double inv_r = 1.0 / sqrt(r2);
            double inv_r3 = inv_r * inv_r * inv_r;

            U -= gm2 * inv_r;

            double f = gm * inv_r3;
            acc[j].x += f * dx;
            acc[j].y += f * dy;
            acc[j].z += f * dz;

            acc[k].x -= f * dx;
            acc[k].y -= f * dy;
            acc[k].z -= f * dz;
        }
    }
    return U;
}

// 運動エネルギー K の計算
double compute_kinetic(const Particle p[N]) {
    double v2_sum = 0.0;
    for (int i = 0; i < N; i++) {
        v2_sum += p[i].vel.x * p[i].vel.x
                + p[i].vel.y * p[i].vel.y
                + p[i].vel.z * p[i].vel.z;
    }
    return 0.5 * MASS * v2_sum;
}

// 重心から測った半質量半径の計算
double compute_half_mass_radius(const Particle p[N]) {
    Vec3 cm = {0.0, 0.0, 0.0};
    for (int i = 0; i < N; i++) {
        cm.x += p[i].pos.x;
        cm.y += p[i].pos.y;
        cm.z += p[i].pos.z;
    }
    cm.x /= N;
    cm.y /= N;
    cm.z /= N;

    double dist[N];
    for (int i = 0; i < N; i++) {
        double dx = p[i].pos.x - cm.x;
        double dy = p[i].pos.y - cm.y;
        double dz = p[i].pos.z - cm.z;
        dist[i] = sqrt(dx * dx + dy * dy + dz * dz);
    }

    qsort(dist, N, sizeof(double), compare_double);

    // 全質量の半分（N/2 個分）を含む半径
    if (N % 2 == 0) {
        return 0.5 * (dist[N / 2 - 1] + dist[N / 2]);
    } else {
        return dist[N / 2];
    }
}

int main(void) {
    int steps = (int)(TOTAL_TIME / DT);
    double target_rv[] = {0.1, 0.2, 0.3, 0.4};
    int num_cases = sizeof(target_rv) / sizeof(target_rv[0]);

    Particle (*particles)[N] = malloc(sizeof(Particle) * (steps + 1) * N);
    if (particles == NULL) {
        fprintf(stderr, "Memory allocation error\n");
        return 1;
    }

    Vec3 acc[N];

    for (int c = 0; c < num_cases; c++) {
        double rv0 = target_rv[c];
        init_genrand(1245); // 各初期ビリアル比で同じ初期配置から比較できるようシードを固定

        // 1. 初期位置（半径 R=1 の球内に一様分布）と初期速度（正規分布）の生成
        Vec3 v_cm = {0.0, 0.0, 0.0};
        for (int i = 0; i < N; i++) {
            Vec3 p;
            do {
                p = (Vec3){ran_sym(), ran_sym(), ran_sym()};
            } while (p.x * p.x + p.y * p.y + p.z * p.z > 1.0);

            Vec3 v = (Vec3){ran_gauss(), ran_gauss(), ran_gauss()};
            particles[0][i].pos = p;
            particles[0][i].vel = v;

            v_cm.x += v.x;
            v_cm.y += v.y;
            v_cm.z += v.z;
        }
        v_cm.x /= N;
        v_cm.y /= N;
        v_cm.z /= N;

        // 2. 重心速度を差し引く
        for (int i = 0; i < N; i++) {
            particles[0][i].vel.x -= v_cm.x;
            particles[0][i].vel.y -= v_cm.y;
            particles[0][i].vel.z -= v_cm.z;
        }

        // 3. 初期ビリアル比 r_v = K / |U| が目標値になるよう全速度をスケーリング
        double U = compute_acc_and_pot(particles[0], acc);
        double K = compute_kinetic(particles[0]);
        double scale = sqrt((rv0 * fabs(U)) / K);

        for (int i = 0; i < N; i++) {
            particles[0][i].vel.x *= scale;
            particles[0][i].vel.y *= scale;
            particles[0][i].vel.z *= scale;
        }

        // 出力ファイルを開く
        char filename[64];
        snprintf(filename, sizeof(filename), "nbody_rv%.1f.dat", rv0);
        FILE *fp = fopen(filename, "w");
        if (fp == NULL) {
            fprintf(stderr, "File open error: %s\n", filename);
            free(particles);
            return 1;
        }
        fprintf(fp, "# t  U  K  E_total  virial_ratio(K/|U|)  half_mass_radius\n");

        // 4. シンプレクティック積分法（速度ベルレ法）による時間発展
        for (int step = 0; step <= steps; step++) {
            double t = step * DT;
            double K_now = compute_kinetic(particles[step]);
            double rv_now = K_now / fabs(U);
            double r_half = compute_half_mass_radius(particles[step]);

            fprintf(fp, "%.5f %.8e %.8e %.8e %.8e %.8e\n",
                    t, U, K_now, K_now + U, rv_now, r_half);

            if (step == steps) break;

            // 半ステップ速度更新 & 位置更新
            Vec3 v_half[N];
            for (int i = 0; i < N; i++) {
                v_half[i].x = particles[step][i].vel.x + 0.5 * acc[i].x * DT;
                v_half[i].y = particles[step][i].vel.y + 0.5 * acc[i].y * DT;
                v_half[i].z = particles[step][i].vel.z + 0.5 * acc[i].z * DT;

                particles[step + 1][i].pos.x = particles[step][i].pos.x + v_half[i].x * DT;
                particles[step + 1][i].pos.y = particles[step][i].pos.y + v_half[i].y * DT;
                particles[step + 1][i].pos.z = particles[step][i].pos.z + v_half[i].z * DT;
            }

            // 新しい位置での加速度とポテンシャルエネルギーを計算
            U = compute_acc_and_pot(particles[step + 1], acc);

            // 残り半ステップの速度更新
            for (int i = 0; i < N; i++) {
                particles[step + 1][i].vel.x = v_half[i].x + 0.5 * acc[i].x * DT;
                particles[step + 1][i].vel.y = v_half[i].y + 0.5 * acc[i].y * DT;
                particles[step + 1][i].vel.z = v_half[i].z + 0.5 * acc[i].z * DT;
            }
        }

        fclose(fp);
        printf("Finished rv = %.1f -> %s\n", rv0, filename);
    }

    free(particles);
    return 0;
}