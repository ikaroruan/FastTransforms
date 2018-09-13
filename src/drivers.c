// Driver routines for synthesis and analysis of harmonic polynomial transforms.

#include "fasttransforms.h"

void execute_sph_hi2lo(const RotationPlan * RP, double * A, const int M) {
    int N = RP->n;
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS()) {
        kernel_sph_hi2lo(RP, m, A + N*(2*m-1));
        kernel_sph_hi2lo(RP, m, A + N*(2*m));
    }
}

void execute_sph_lo2hi(const RotationPlan * RP, double * A, const int M) {
    int N = RP->n;
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS()) {
        kernel_sph_lo2hi(RP, m, A + N*(2*m-1));
        kernel_sph_lo2hi(RP, m, A + N*(2*m));
    }
}

void execute_sph_hi2lo_SSE(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_sph(A, B, N, M, 2);
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS())
        kernel_sph_hi2lo_SSE(RP, m, B + NB*(2*m-1));
    permute_t_sph(A, B, N, M, 2);
}

void execute_sph_lo2hi_SSE(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_sph(A, B, N, M, 2);
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS())
        kernel_sph_lo2hi_SSE(RP, m, B + NB*(2*m-1));
    permute_t_sph(A, B, N, M, 2);
}

void execute_sph_hi2lo_AVX(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    warp(A, N, M, 2);
    permute_sph(A, B, N, M, 4);
    for (int m = 2; m <= (M%8)/2; m++)
        kernel_sph_hi2lo_SSE(RP, m, B + NB*(2*m-1));
    #pragma omp parallel
    for (int m = (M%8+1)/2 + 4*FT_GET_THREAD_NUM(); m <= M/2; m += 4*FT_GET_NUM_THREADS()) {
        kernel_sph_hi2lo_AVX(RP, m, B + NB*(2*m-1));
        kernel_sph_hi2lo_AVX(RP, m+1, B + NB*(2*m+3));
    }
    permute_t_sph(A, B, N, M, 4);
    warp(A, N, M, 2);
}

void execute_sph_lo2hi_AVX(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    warp(A, N, M, 2);
    permute_sph(A, B, N, M, 4);
    for (int m = 2; m <= (M%8)/2; m++)
        kernel_sph_lo2hi_SSE(RP, m, B + NB*(2*m-1));
    #pragma omp parallel
    for (int m = (M%8+1)/2 + 4*FT_GET_THREAD_NUM(); m <= M/2; m += 4*FT_GET_NUM_THREADS()) {
        kernel_sph_lo2hi_AVX(RP, m, B + NB*(2*m-1));
        kernel_sph_lo2hi_AVX(RP, m+1, B + NB*(2*m+3));
    }
    permute_t_sph(A, B, N, M, 4);
    warp(A, N, M, 2);
}

void execute_sph_hi2lo_AVX512(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    int M_star = M%16;
    warp(A, N, M, 4);
    warp(A, N, M_star, 2);
    permute_sph(A, B, N, M, 8);
    for (int m = 2; m <= (M_star%8)/2; m++)
        kernel_sph_hi2lo_SSE(RP, m, B + NB*(2*m-1));
    for (int m = (M_star%8+1)/2; m <= M_star/2; m += 4) {
        kernel_sph_hi2lo_AVX(RP, m, B + NB*(2*m-1));
        kernel_sph_hi2lo_AVX(RP, m+1, B + NB*(2*m+3));
    }
    #pragma omp parallel
    for (int m = (M_star+1)/2 + 8*FT_GET_THREAD_NUM(); m <= M/2; m += 8*FT_GET_NUM_THREADS()) {
        kernel_sph_hi2lo_AVX512(RP, m, B + NB*(2*m-1));
        kernel_sph_hi2lo_AVX512(RP, m+1, B + NB*(2*m+7));
    }
    permute_t_sph(A, B, N, M, 8);
    warp(A, N, M_star, 2);
    warp_t(A, N, M, 4);
}

void execute_sph_lo2hi_AVX512(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    int M_star = M%16;
    warp(A, N, M, 4);
    warp(A, N, M_star, 2);
    permute_sph(A, B, N, M, 8);
    for (int m = 2; m <= (M_star%8)/2; m++)
        kernel_sph_lo2hi_SSE(RP, m, B + NB*(2*m-1));
    for (int m = (M_star%8+1)/2; m <= M_star/2; m += 4) {
        kernel_sph_lo2hi_AVX(RP, m, B + NB*(2*m-1));
        kernel_sph_lo2hi_AVX(RP, m+1, B + NB*(2*m+3));
    }
    #pragma omp parallel
    for (int m = (M_star+1)/2 + 8*FT_GET_THREAD_NUM(); m <= M/2; m += 8*FT_GET_NUM_THREADS()) {
        kernel_sph_lo2hi_AVX512(RP, m, B + NB*(2*m-1));
        kernel_sph_lo2hi_AVX512(RP, m+1, B + NB*(2*m+7));
    }
    permute_t_sph(A, B, N, M, 8);
    warp(A, N, M_star, 2);
    warp_t(A, N, M, 4);
}

void execute_tri_hi2lo(const RotationPlan * RP, double * A, const int M) {
    #pragma omp parallel
    for (int m = 1 + FT_GET_THREAD_NUM(); m < M; m += FT_GET_NUM_THREADS())
        kernel_tri_hi2lo(RP, m, A+(RP->n)*m);
}

void execute_tri_lo2hi(const RotationPlan * RP, double * A, const int M) {
    #pragma omp parallel
    for (int m = 1 + FT_GET_THREAD_NUM(); m < M; m += FT_GET_NUM_THREADS())
        kernel_tri_lo2hi(RP, m, A+(RP->n)*m);
}

void execute_tri_hi2lo_SSE(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_tri(A, B, N, M, 2);
    #pragma omp parallel
    for (int m = M%2+2*FT_GET_THREAD_NUM(); m < M; m += 2*FT_GET_NUM_THREADS())
        kernel_tri_hi2lo_SSE(RP, m, B+NB*m);
    permute_t_tri(A, B, N, M, 2);
}

void execute_tri_lo2hi_SSE(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_tri(A, B, N, M, 2);
    #pragma omp parallel
    for (int m = M%2+2*FT_GET_THREAD_NUM(); m < M; m += 2*FT_GET_NUM_THREADS())
        kernel_tri_lo2hi_SSE(RP, m, B+NB*m);
    permute_t_tri(A, B, N, M, 2);
}

void execute_tri_hi2lo_AVX(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_tri(A, B, N, M, 4);
    for (int m = M%2; m < M%8; m += 2)
        kernel_tri_hi2lo_SSE(RP, m, B+NB*m);
    #pragma omp parallel
    for (int m = M%8 + 4*FT_GET_THREAD_NUM(); m < M; m += 4*FT_GET_NUM_THREADS())
        kernel_tri_hi2lo_AVX(RP, m, B+NB*m);
    permute_t_tri(A, B, N, M, 4);
}

void execute_tri_lo2hi_AVX(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_tri(A, B, N, M, 4);
    for (int m = M%2; m < M%8; m += 2)
        kernel_tri_lo2hi_SSE(RP, m, B+NB*m);
    #pragma omp parallel
    for (int m = M%8 + 4*FT_GET_THREAD_NUM(); m < M; m += 4*FT_GET_NUM_THREADS())
        kernel_tri_lo2hi_AVX(RP, m, B+NB*m);
    permute_t_tri(A, B, N, M, 4);
}

void execute_tri_hi2lo_AVX512(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    int M_star = M%16;
    permute_tri(A, B, N, M, 8);
    for (int m = M%2; m < M_star%8; m += 2)
        kernel_tri_hi2lo_SSE(RP, m, B+NB*m);
    for (int m = M_star%8; m < M%16; m += 4)
        kernel_tri_hi2lo_AVX(RP, m, B+NB*m);
    #pragma omp parallel
    for (int m = M_star + 8*FT_GET_THREAD_NUM(); m < M; m += 8*FT_GET_NUM_THREADS())
        kernel_tri_hi2lo_AVX512(RP, m, B+NB*m);
    permute_t_tri(A, B, N, M, 8);
}

void execute_tri_lo2hi_AVX512(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    int M_star = M%16;
    permute_tri(A, B, N, M, 8);
    for (int m = M%2; m < M_star%8; m += 2)
        kernel_tri_lo2hi_SSE(RP, m, B+NB*m);
    for (int m = M_star%8; m < M%16; m += 4)
        kernel_tri_lo2hi_AVX(RP, m, B+NB*m);
    #pragma omp parallel
    for (int m = M_star + 8*FT_GET_THREAD_NUM(); m < M; m += 8*FT_GET_NUM_THREADS())
        kernel_tri_lo2hi_AVX512(RP, m, B+NB*m);
    permute_t_tri(A, B, N, M, 8);
}


void execute_disk_hi2lo(const RotationPlan * RP, double * A, const int M) {
    int N = RP->n;
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS()) {
        kernel_disk_hi2lo(RP, m, A + N*(2*m-1));
        kernel_disk_hi2lo(RP, m, A + N*(2*m));
    }
}

void execute_disk_lo2hi(const RotationPlan * RP, double * A, const int M) {
    int N = RP->n;
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS()) {
        kernel_disk_lo2hi(RP, m, A + N*(2*m-1));
        kernel_disk_lo2hi(RP, m, A + N*(2*m));
    }
}

void execute_disk_hi2lo_SSE(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_disk(A, B, N, M, 2);
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS())
        kernel_disk_hi2lo_SSE(RP, m, B + NB*(2*m-1));
    permute_t_disk(A, B, N, M, 2);
}

void execute_disk_lo2hi_SSE(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    permute_disk(A, B, N, M, 2);
    #pragma omp parallel
    for (int m = 2 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS())
        kernel_disk_lo2hi_SSE(RP, m, B + NB*(2*m-1));
    permute_t_disk(A, B, N, M, 2);
}

void execute_disk_hi2lo_AVX(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    warp(A, N, M, 2);
    permute_disk(A, B, N, M, 4);
    for (int m = 2; m <= (M%8)/2; m++)
        kernel_disk_hi2lo_SSE(RP, m, B + NB*(2*m-1));
    #pragma omp parallel
    for (int m = (M%8+1)/2 + 4*FT_GET_THREAD_NUM(); m <= M/2; m += 4*FT_GET_NUM_THREADS()) {
        kernel_disk_hi2lo_AVX(RP, m, B + NB*(2*m-1));
        kernel_disk_hi2lo_AVX(RP, m+1, B + NB*(2*m+3));
    }
    permute_t_disk(A, B, N, M, 4);
    warp(A, N, M, 2);
}

void execute_disk_lo2hi_AVX(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    warp(A, N, M, 2);
    permute_disk(A, B, N, M, 4);
    for (int m = 2; m <= (M%8)/2; m++)
        kernel_disk_lo2hi_SSE(RP, m, B + NB*(2*m-1));
    #pragma omp parallel
    for (int m = (M%8+1)/2 + 4*FT_GET_THREAD_NUM(); m <= M/2; m += 4*FT_GET_NUM_THREADS()) {
        kernel_disk_lo2hi_AVX(RP, m, B + NB*(2*m-1));
        kernel_disk_lo2hi_AVX(RP, m+1, B + NB*(2*m+3));
    }
    permute_t_disk(A, B, N, M, 4);
    warp(A, N, M, 2);
}

void execute_disk_hi2lo_AVX512(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    int M_star = M%16;
    warp(A, N, M, 4);
    warp(A, N, M_star, 2);
    permute_disk(A, B, N, M, 8);
    for (int m = 2; m <= (M_star%8)/2; m++)
        kernel_disk_hi2lo_SSE(RP, m, B + NB*(2*m-1));
    for (int m = (M_star%8+1)/2; m <= M_star/2; m += 4) {
        kernel_disk_hi2lo_AVX(RP, m, B + NB*(2*m-1));
        kernel_disk_hi2lo_AVX(RP, m+1, B + NB*(2*m+3));
    }
    #pragma omp parallel
    for (int m = (M_star+1)/2 + 8*FT_GET_THREAD_NUM(); m <= M/2; m += 8*FT_GET_NUM_THREADS()) {
        kernel_disk_hi2lo_AVX512(RP, m, B + NB*(2*m-1));
        kernel_disk_hi2lo_AVX512(RP, m+1, B + NB*(2*m+7));
    }
    permute_t_disk(A, B, N, M, 8);
    warp(A, N, M_star, 2);
    warp_t(A, N, M, 4);
}

void execute_disk_lo2hi_AVX512(const RotationPlan * RP, double * A, double * B, const int M) {
    int N = RP->n;
    int NB = ALIGNB(N);
    int M_star = M%16;
    warp(A, N, M, 4);
    warp(A, N, M_star, 2);
    permute_disk(A, B, N, M, 8);
    for (int m = 2; m <= (M_star%8)/2; m++)
        kernel_disk_lo2hi_SSE(RP, m, B + NB*(2*m-1));
    for (int m = (M_star%8+1)/2; m <= M_star/2; m += 4) {
        kernel_disk_lo2hi_AVX(RP, m, B + NB*(2*m-1));
        kernel_disk_lo2hi_AVX(RP, m+1, B + NB*(2*m+3));
    }
    #pragma omp parallel
    for (int m = (M_star+1)/2 + 8*FT_GET_THREAD_NUM(); m <= M/2; m += 8*FT_GET_NUM_THREADS()) {
        kernel_disk_lo2hi_AVX512(RP, m, B + NB*(2*m-1));
        kernel_disk_lo2hi_AVX512(RP, m+1, B + NB*(2*m+7));
    }
    permute_t_disk(A, B, N, M, 8);
    warp(A, N, M_star, 2);
    warp_t(A, N, M, 4);
}


void execute_spinsph_hi2lo(const SpinRotationPlan * SRP, double * A, const int M) {
    int N = SRP->n;
    kernel_spinsph_hi2lo(SRP, 0, A);
    #pragma omp parallel
    for (int m = 1 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS()) {
        kernel_spinsph_hi2lo(SRP, m, A + N*(2*m-1));
        kernel_spinsph_hi2lo(SRP, m, A + N*(2*m));
    }
}

void execute_spinsph_lo2hi(const SpinRotationPlan * SRP, double * A, const int M) {
    int N = SRP->n;
    kernel_spinsph_lo2hi(SRP, 0, A);
    #pragma omp parallel
    for (int m = 1 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS()) {
        kernel_spinsph_lo2hi(SRP, m, A + N*(2*m-1));
        kernel_spinsph_lo2hi(SRP, m, A + N*(2*m));
    }
}

void execute_spinsph_hi2lo_SSE(const SpinRotationPlan * SRP, double * A, double * B, const int M) {
    int N = SRP->n;
    int NB = ALIGNB(N);
    permute_spinsph(A, B, N, M, 2);
    kernel_spinsph_hi2lo(SRP, 0, B);
    #pragma omp parallel
    for (int m = 1 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS())
        kernel_spinsph_hi2lo_SSE(SRP, m, B + NB*(2*m-1));
   permute_t_spinsph(A, B, N, M, 2);
}

void execute_spinsph_lo2hi_SSE(const SpinRotationPlan * SRP, double * A, double * B, const int M) {
    int N = SRP->n;
    int NB = ALIGNB(N);
    permute_spinsph(A, B, N, M, 2);
    kernel_spinsph_lo2hi(SRP, 0, B);
    #pragma omp parallel
    for (int m = 1 + FT_GET_THREAD_NUM(); m <= M/2; m += FT_GET_NUM_THREADS())
        kernel_spinsph_lo2hi_SSE(SRP, m, B + NB*(2*m-1));
   permute_t_spinsph(A, B, N, M, 2);
}

void execute_spinsph_hi2lo_AVX(const SpinRotationPlan * SRP, double * A, double * B, const int M) {
    int N = SRP->n;
    int NB = ALIGNB(N);
    warp(A, N, M, 2);
    permute_spinsph(A, B, N, M, 4);
    kernel_spinsph_hi2lo(SRP, 0, B);
    for (int m = 1; m <= (M%8)/2; m++)
        kernel_spinsph_hi2lo_SSE(SRP, m, B + NB*(2*m-1));
    #pragma omp parallel
    for (int m = (M%8+1)/2 + 4*FT_GET_THREAD_NUM(); m <= M/2; m += 4*FT_GET_NUM_THREADS()) {
        kernel_spinsph_hi2lo_AVX(SRP, m, B + NB*(2*m-1));
        kernel_spinsph_hi2lo_AVX(SRP, m+1, B + NB*(2*m+3));
    }
   permute_t_spinsph(A, B, N, M, 4);
   warp(A, N, M, 2);
}

void execute_spinsph_lo2hi_AVX(const SpinRotationPlan * SRP, double * A, double * B, const int M) {
    int N = SRP->n;
    int NB = ALIGNB(N);
    warp(A, N, M, 2);
    permute_spinsph(A, B, N, M, 4);
    kernel_spinsph_lo2hi(SRP, 0, B);
    for (int m = 1; m <= (M%8)/2; m++)
        kernel_spinsph_lo2hi_SSE(SRP, m, B + NB*(2*m-1));
    #pragma omp parallel
    for (int m = (M%8+1)/2 + 4*FT_GET_THREAD_NUM(); m <= M/2; m += 4*FT_GET_NUM_THREADS()) {
        kernel_spinsph_lo2hi_AVX(SRP, m, B + NB*(2*m-1));
        kernel_spinsph_lo2hi_AVX(SRP, m+1, B + NB*(2*m+3));
    }
   permute_t_spinsph(A, B, N, M, 4);
   warp(A, N, M, 2);
}

void execute_spinsph_hi2lo_AVX512(const SpinRotationPlan * SRP, double * A, double * B, const int M) {
    int N = SRP->n;
    int NB = ALIGNB(N);
    int M_star = M%16;
    warp(A, N, M, 4);
    warp(A, N, M_star, 2);
    permute_spinsph(A, B, N, M, 8);
    kernel_spinsph_hi2lo(SRP, 0, B);
    for (int m = 1; m <= (M_star%8)/2; m++)
        kernel_spinsph_hi2lo_SSE(SRP, m, B + NB*(2*m-1));
    for (int m = (M_star%8+1)/2; m <= M_star/2; m += 4) {
        kernel_spinsph_hi2lo_AVX(SRP, m, B + NB*(2*m-1));
        kernel_spinsph_hi2lo_AVX(SRP, m+1, B + NB*(2*m+3));
    }
    #pragma omp parallel
    for (int m = (M_star+1)/2 + 8*FT_GET_THREAD_NUM(); m <= M/2; m += 8*FT_GET_NUM_THREADS()) {
        kernel_spinsph_hi2lo_AVX512(SRP, m, B + NB*(2*m-1));
        kernel_spinsph_hi2lo_AVX512(SRP, m+1, B + NB*(2*m+7));
    }
    permute_t_spinsph(A, B, N, M, 8);
    warp(A, N, M_star, 2);
    warp_t(A, N, M, 4);
}

void execute_spinsph_lo2hi_AVX512(const SpinRotationPlan * SRP, double * A, double * B, const int M) {
    int N = SRP->n;
    int M_star = M%16;
    int NB = ALIGNB(N);
    warp(A, N, M, 4);
    warp(A, N, M_star, 2);
    permute_spinsph(A, B, N, M, 8);
    kernel_spinsph_lo2hi(SRP, 0, B);
    for (int m = 1; m <= (M_star%8)/2; m++)
        kernel_spinsph_lo2hi_SSE(SRP, m, B + NB*(2*m-1));
    for (int m = (M_star%8+1)/2; m <= M_star/2; m += 4) {
        kernel_spinsph_lo2hi_AVX(SRP, m, B + NB*(2*m-1));
        kernel_spinsph_lo2hi_AVX(SRP, m+1, B + NB*(2*m+3));
    }
    #pragma omp parallel
    for (int m = (M_star+1)/2 + 8*FT_GET_THREAD_NUM(); m <= M/2; m += 8*FT_GET_NUM_THREADS()) {
        kernel_spinsph_lo2hi_AVX512(SRP, m, B + NB*(2*m-1));
        kernel_spinsph_lo2hi_AVX512(SRP, m+1, B + NB*(2*m+7));
    }
    permute_t_spinsph(A, B, N, M, 8);
    warp(A, N, M_star, 2);
    warp_t(A, N, M, 4);
}


void freeHarmonicPlan(HarmonicPlan * P) {
    freeRotationPlan(P->RP);
    VFREE(P->B);
    free(P->P1);
    free(P->P2);
    free(P->P1inv);
    free(P->P2inv);
    free(P);
}

HarmonicPlan * plan_sph2fourier(const int n) {
    HarmonicPlan * P = malloc(sizeof(HarmonicPlan));
    P->RP = plan_rotsphere(n);
    P->B = (double *) VMALLOC(ALIGNB(n) * (2*n-1) * sizeof(double));
    P->P1 = plan_leg2cheb(1, 0, n);
    P->P2 = plan_ultra2ultra(1, 0, n, 1.5, 1.0);
    P->P1inv = plan_cheb2leg(0, 1, n);
    P->P2inv = plan_ultra2ultra(0, 1, n, 1.0, 1.5);
    return P;
}

void execute_sph2fourier(const HarmonicPlan * P, double * A, const int N, const int M) {
    execute_sph_hi2lo_AVX512(P->RP, A, P->B, M);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+3)/4, 1.0, P->P1, N, A, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+2)/4, 1.0, P->P2, N, A+N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+1)/4, 1.0, P->P2, N, A+2*N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, M/4, 1.0, P->P1, N, A+3*N, 4*N);
}

void execute_fourier2sph(const HarmonicPlan * P, double * A, const int N, const int M) {
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+3)/4, 1.0, P->P1inv, N, A, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+2)/4, 1.0, P->P2inv, N, A+N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+1)/4, 1.0, P->P2inv, N, A+2*N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, M/4, 1.0, P->P1inv, N, A+3*N, 4*N);
    execute_sph_lo2hi_AVX512(P->RP, A, P->B, M);
}

HarmonicPlan * plan_tri2cheb(const int n, const double alpha, const double beta, const double gamma) {
    HarmonicPlan * P = malloc(sizeof(HarmonicPlan));
    P->RP = plan_rottriangle(n, alpha, beta, gamma);
    P->B = (double *) VMALLOC(ALIGNB(n) * n * sizeof(double));
    P->P1 = plan_jac2jac(1, 1, n, beta + gamma + 1.0, alpha, -0.5);
    double * P12 = plan_jac2jac(1, 1, n, alpha, -0.5, -0.5);
    alternate_sign(P12, n, n);
    alternate_sign_t(P12, n, n);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1.0, P12, n, P->P1, n);
    free(P12);
    P->P2 = plan_jac2jac(1, 1, n, gamma, beta, -0.5);
    double * P22 = plan_jac2jac(1, 1, n, beta, -0.5, -0.5);
    alternate_sign(P22, n, n);
    alternate_sign_t(P22, n, n);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1.0, P22, n, P->P2, n);
    free(P22);
    P->P1inv = plan_jac2jac(1, 1, n, -0.5, -0.5, alpha);
    double * P12inv = plan_jac2jac(1, 1, n, -0.5, alpha, beta + gamma + 1.0);
    alternate_sign(P->P1inv, n, n);
    alternate_sign_t(P->P1inv, n, n);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1.0, P12inv, n, P->P1inv, n);
    free(P12inv);
    P->P2inv = plan_jac2jac(1, 1, n, -0.5, -0.5, beta);
    double * P22inv = plan_jac2jac(1, 1, n, -0.5, beta, gamma);
    alternate_sign(P->P2inv, n, n);
    alternate_sign_t(P->P2inv, n, n);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1.0, P22inv, n, P->P2inv, n);
    free(P22inv);
    P->alpha = alpha;
    P->beta = beta;
    P->gamma = gamma;
    return P;
}

void execute_tri2cheb(const HarmonicPlan * P, double * A, const int N, const int M) {
    execute_tri_hi2lo_AVX512(P->RP, A, P->B, M);
    if ((P->beta + P->gamma != -1.5) || (P->alpha != -0.5))
        cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, M, 1.0, P->P1, N, A, N);
    if ((P->gamma != -0.5) || (P->beta != -0.5))
        cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit, N, M, 1.0, P->P2, N, A, N);
    chebyshev_normalization(A, N, M);
}

void execute_cheb2tri(const HarmonicPlan * P, double * A, const int N, const int M) {
    chebyshev_normalization_t(A, N, M);
    if ((P->beta != -0.5) || (P->gamma != -0.5))
        cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit, N, M, 1.0, P->P2inv, N, A, N);
    if ((P->alpha != -0.5) || (P->beta + P->gamma != -1.5))
        cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, M, 1.0, P->P1inv, N, A, N);
    execute_tri_lo2hi_AVX512(P->RP, A, P->B, M);
}

HarmonicPlan * plan_disk2cxf(const int n) {
    HarmonicPlan * P = malloc(sizeof(HarmonicPlan));
    P->RP = plan_rotdisk(n);
    P->B = (double *) VMALLOC(ALIGNB(n) * (4*n-3) * sizeof(double));
    P->P1 = plan_leg2cheb(1, 0, n);
    P->P2 = plan_jac2jac(1, 1, n, 1.0, 0.0, 0.5);
    alternate_sign(P->P2, n, n);
    alternate_sign_t(P->P2, n, n);
    double * P22 = plan_jac2jac(1, 1, n, 0.0, 0.5, -0.5);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1.0, P22, n, P->P2, n);
    free(P22);
    P->P1inv = plan_cheb2leg(0, 1, n);
    double * P22inv = plan_jac2jac(1, 1, n, 0.5, 0.0, 1.0);
    alternate_sign(P22inv, n, n);
    alternate_sign_t(P22inv, n, n);
    P->P2inv = plan_jac2jac(1, 1, n, -0.5, 0.5, 0.0);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, 1.0, P22inv, n, P->P2inv, n);
    free(P22inv);
    return P;
}

void execute_disk2cxf(const HarmonicPlan * P, double * A, const int N, const int M) {
    execute_disk_hi2lo_AVX512(P->RP, A, P->B, M);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+3)/4, 1.0, P->P1, N, A, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+2)/4, 1.0, P->P2, N, A+N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+1)/4, 1.0, P->P2, N, A+2*N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, M/4, 1.0, P->P1, N, A+3*N, 4*N);
    partial_chebyshev_normalization(A, N, M);
}

void execute_cxf2disk(const HarmonicPlan * P, double * A, const int N, const int M) {
    partial_chebyshev_normalization_t(A, N, M);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+3)/4, 1.0, P->P1inv, N, A, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+2)/4, 1.0, P->P2inv, N, A+N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, (M+1)/4, 1.0, P->P2inv, N, A+2*N, 4*N);
    cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, M/4, 1.0, P->P1inv, N, A+3*N, 4*N);
    execute_disk_lo2hi_AVX512(P->RP, A, P->B, M);
}

static void alternate_sign(double * A, const int N, const int M) {
    for (int j = 0; j < M; j++)
        for (int i = 0; i < N; i += 2)
            A[i+j*N] = -A[i+j*N];
}

static void alternate_sign_t(double * A, const int N, const int M) {
    for (int i = 0; i < N; i += 2)
        for (int j = 0; j < M; j++)
            A[j+i*M] = -A[j+i*M];
}

static void chebyshev_normalization(double * A, const int N, const int M) {
    A[0] *= M_1_PI;
    for (int i = 1; i < N; i++)
        A[i] *= M_SQRT2*M_1_PI;
    for (int j = 1; j < M; j++)
        A[j*N] *= M_SQRT2*M_1_PI;
    for (int j = 1; j < M; j++)
        for (int i = 1; i < N; i++)
            A[i+j*N] *= M_2_PI;
}

static void chebyshev_normalization_t(double * A, const int N, const int M) {
    A[0] *= M_PI;
    for (int i = 1; i < N; i++)
        A[i] *= M_SQRT1_2*M_PI;
    for (int j = 1; j < M; j++)
        A[j*N] *= M_SQRT1_2*M_PI;
    for (int j = 1; j < M; j++)
        for (int i = 1; i < N; i++)
            A[i+j*N] *= M_PI_2;
}

static void partial_chebyshev_normalization(double * A, const int N, const int M) {
    for (int j = 1; j < M; j += 4) {
        A[j*N] *= M_1_SQRT_PI;
        for (int i = 1; i < N; i++)
            A[i+j*N] *= M_SQRT2*M_1_SQRT_PI;
    }
    for (int j = 2; j < M; j += 4) {
        A[j*N] *= M_1_SQRT_PI;
        for (int i = 1; i < N; i++)
            A[i+j*N] *= M_SQRT2*M_1_SQRT_PI;
    }
}

static void partial_chebyshev_normalization_t(double * A, const int N, const int M) {
    for (int j = 1; j < M; j += 4) {
        A[j*N] *= M_SQRT_PI;
        for (int i = 1; i < N; i++)
            A[i+j*N] *= M_SQRT1_2*M_SQRT_PI;
    }
    for (int j = 2; j < M; j += 4) {
        A[j*N] *= M_SQRT_PI;
        for (int i = 1; i < N; i++)
            A[i+j*N] *= M_SQRT1_2*M_SQRT_PI;
    }
}
