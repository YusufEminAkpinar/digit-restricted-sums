#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpfr.h>
#include <gmp.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

// I know this is double eval, I just don't care rn. It is not a public library,
// only I use it and it's a small 'script'.
#define MIN(X, Y) ( ((X) < (Y)) ? (X) : (Y) ) 
#define C(X) ( ceil((double)X/2) )
#define GUARD 12
#define PREC 3000
#define BITS ceil( (PREC+1) / log10(2) )
#define L 3
#define LIMIT 1 + PREC/(L-1)
#define PRINTPERCALC 100

// #define LOG2 0.301029995663981195213738894724493026768189 
// #define STEP (int)(PREC/10.0) < 500 ? (int)(PREC/10.0) : 500
// #define STEP (int)(PREC/10.0)
#define STEP 500
#define RATIO (int)ceil((double) STEP / (L-1) / log10(2))
#define Mmax (int)floor( (BITS - (double)GUARD/2)/(L-1)*log10(2) )

typedef unsigned long int ulongi;

#define PRINT
// #define TIME

int binomMax = 10;
int binomPerLine = 10;

void um(mpfr_t u[LIMIT], mpz_t gammas[LIMIT], mpz_t *combi[binomPerLine], int decr);
void gamma_val(mpz_t gammas[LIMIT]);
void generate_binom_first(mpz_t *combi[binomPerLine]);
void get_binom(mpz_t binom, mpz_t *combi[binomPerLine], int x, int y);
void free_combi(mpz_t *combi[binomPerLine]);
bool isAdmissible(int x);
ulongi powi(int base, int exp);
void beta(mpfr_t betas[LIMIT], mpfr_t betaLimit, int decr);
void sommeFinie(mpfr_t value);
void sommeInfinite(mpfr_t value, mpfr_t u[LIMIT], mpfr_t beta[LIMIT]);


int main(int argc, char *argv[]){
	(void) argc;
	(void) argv;

	const int nbPrec = 1 + (int)floor( (BITS - (double)GUARD/2)/STEP );
	const int decr = 1 + (double)Mmax/nbPrec;
	printf("STEP: %d, NBPREC: %d, RATIO: %d, Mmax: %d, Decr: %d, PREC: %d\n", STEP, nbPrec, RATIO, Mmax, decr, PREC);
	
	mpfr_t ten; mpfr_init_set_ui(ten, 10, MPFR_RNDN);
	mpfr_t betaLimit; mpfr_init2(betaLimit, BITS+GUARD);
	mpfr_pow_si(betaLimit, ten, -1.1 * PREC, MPFR_RNDN);

	mpz_t *combi[binomPerLine];
	mpz_t gammas[LIMIT];
	mpfr_t u[LIMIT];
	mpfr_t betas[LIMIT];

	mpfr_t result;
	mpfr_init2(result, BITS+GUARD);
	mpfr_set_ui(result, 0, MPFR_RNDN);

	generate_binom_first(combi);
	mpz_t binom; mpz_init(binom);


	printf("Calculating gammas.\n");
	gamma_val(gammas);


	printf("Gammas finished.\nStarting um.\n");
	um(u, gammas, combi, decr);

	// for (int i = 0; i <= LIMIT; i++) {
	// 	printf("u[%d] = ", i);
	// 	mpfr_out_str(stdout, 10, 0, u[i], MPFR_RNDN);
	// 	printf("\n");
	// }

	printf("Calculating betas.\n");
	beta(betas, betaLimit, decr);
	printf("Betas finished.\nStarting finite sum.\n");

	sommeFinie(result);
	sommeInfinite(result, u, betas);


	printf("Calculation finished. For %f bit precision (%d digits), result is:\n", BITS, PREC);
	mpfr_out_str(stdout, 10, 0, result, MPFR_RNDN);
	printf("\n");


	printf("Calculated %d amount of u.\n", LIMIT);
	free_combi(combi);
	for (int i = 0; i < LIMIT-1; i++) {
		mpz_clear(gammas[i]);
		mpfr_clears(u[i], betas[i], NULL);
	}
	mpz_clear(binom);
	mpfr_clears(betaLimit, result, NULL);

	return 0;
}

void free_combi(mpz_t *combi[binomPerLine]) {

	for (int j = binomMax-binomPerLine; j < binomMax; j++) {
		for (int i = 0; i <= C(j)+1; i++) {
			mpz_clear(combi[j%binomPerLine][i]);
		}
		free(combi[j%binomPerLine]);
	}
}


void generate_binom_first(mpz_t *combi[binomPerLine]) {

	for (int i = 0; i < binomPerLine; i++) {
		combi[i] = malloc(sizeof(mpz_t) * (i + 2));
		memset(combi[i], 0, sizeof(mpz_t) * (i + 2));
	}

	mpz_init(combi[0][0]);
	mpz_init(combi[1][0]);
	mpz_init(combi[1][1]);

	mpz_set_ui(combi[0][0], 1);
	mpz_set_ui(combi[1][0], 1);
	mpz_set_ui(combi[1][1], 1);

	
	for (int j = 2; j < binomPerLine; j++) {
		mpz_init(combi[j][0]);
		mpz_set_ui(combi[j][0], 1);
		for (int i = 1; i <= C(j)+1; i++) {
			mpz_init(combi[j][i]);
			mpz_add(combi[j][i], combi[j-1][MIN(i, j-i)], combi[j-1][MIN(i-1, abs(j-i-1))]);
		}
	}
}

void get_binom(mpz_t binom, mpz_t *combi[binomPerLine], int x, int y) {
    if (x < binomMax) {
		mpz_set(binom, combi[x % binomPerLine][MIN(y, x-y)]);
		return;
	}


	for (int j = 0; j < binomPerLine; j++) {
		combi[j] = (mpz_t *)realloc(combi[j], sizeof(mpz_t) * (j + 2 + binomMax));
		if (combi[j] == NULL) {
			printf("\nCan not reallocate the binomial matrix. Exiting.\n");
			free_combi(combi);
			exit(0);
		}
		// j = 3500 -> calculating binomials for 3000-4000, So our binomMax=3000
		// For example at 100th line (which lies binomials for j = 2100),
		// we already initialized 
		// (100 + binomMax - binomPerLine)(=2100)*sizeof(mpz_t), byte of mem.
		// So for nth line we have n+binomMax elements already, we need to
		// initialize n+binomMax-binomPerLine -> n + binomMax
		// TODO: Actually, we only need C(j) + 2 but a bit of memory doesn't
		// hurt anybody at the beginning, right... right?? Fix it after it works.
		// maybe don't fix it? It just uses lots of memory,
		// doesn't read or write to it.
		for (int i = j + binomMax - binomPerLine; i < j + binomMax+2; i++) {
			mpz_init(combi[j%binomPerLine][i]);
		}
	}


	for (int j = binomMax; j < binomMax + binomPerLine; j++) {
		mpz_set_ui(combi[j%binomPerLine][0], 1);
		for (int i = 1; i <= C(j)+1; i++) {
			mpz_add(combi[j% binomPerLine][i], combi[(j-1)% binomPerLine][MIN(i, j-i)], combi[(j-1)% binomPerLine][MIN(i-1, j-i-1)]);
		}
	}
	binomMax += binomPerLine;
	mpz_set(binom, combi[x% binomPerLine][MIN(y, x-y)]);

}

void gamma_val(mpz_t gammas[LIMIT]) {
	mpz_t pow, res;
	mpz_inits(pow, res, NULL);
	mpz_init(gammas[0]);
	mpz_set_ui(gammas[0], 9);
	for (int k = 1; k < LIMIT; k++) {
		mpz_set_ui(res, 0);
		for (int i = 1; i < 9; i++) {
			mpz_ui_pow_ui(pow, i, k);
			mpz_add(res, res, pow);
		}
		mpz_init(gammas[k]);
		mpz_set(gammas[k], res);
	}
	mpz_clears(pow, res, NULL);
}


void um(mpfr_t u[LIMIT], mpz_t gammas[LIMIT], mpz_t *combi[binomPerLine], int decr) {

	int precisionBits = BITS+GUARD;

	mpz_t den; mpz_init(den);
	mpfr_t denf; mpfr_init2(denf, precisionBits);

	mpfr_init2(u[0], precisionBits); mpfr_set_ui(u[0], 10, MPFR_RNDN);

#ifdef TIME
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	double total_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
#endif

	for (int m = 1; m<LIMIT+1; m++) {

		if (m % decr == 0) {
			int newBits = BITS + GUARD - (int)((double)m/decr) * STEP;
			precisionBits = (newBits > 0) ? newBits : precisionBits;
		}
		mpfr_init2(u[m], precisionBits); mpfr_set_zero(u[m], 0);

		mpfr_t sumGlobal;
		mpfr_init2(sumGlobal, precisionBits);
		mpfr_set_ui(sumGlobal, 0, MPFR_RNDN);

#ifdef TIME
		struct timespec ts;
		clock_gettime(CLOCK_REALTIME, &ts);
		double start_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
		double end_time;
#endif

		#pragma omp parallel
		{
			mpz_t binomial; mpz_init(binomial);

			mpfr_t sumf; mpfr_init2(sumf, precisionBits);
			mpfr_t temp; mpfr_init2(temp, precisionBits);

			mpfr_set_ui(sumf, 0, MPFR_RNDN);

			
			#pragma omp for nowait
			for (int j = 1; j < m+1; j++) {
				#pragma omp critical
				{
					get_binom(binomial, combi, m, j);
				}
				mpz_mul(binomial, binomial, gammas[j]);
				mpfr_mul_z(temp, u[m-j], binomial, MPFR_RNDN);

				mpfr_add(sumf, sumf, temp, MPFR_RNDN);
			}

			#pragma omp critical
			{
				mpfr_add(sumGlobal, sumGlobal, sumf, MPFR_RNDN);
			}
			mpz_clear(binomial);
			mpfr_clears(sumf, temp, NULL);
		}

		mpz_ui_pow_ui(den, 10, m+1);
		mpz_sub_ui(den, den, 9);
		mpfr_set_z(denf, den, MPFR_RNDN);

		mpfr_div(u[m], sumGlobal, denf, MPFR_RNDN);

#ifdef PRINT	
			if (m%PRINTPERCALC == 0) {
#ifdef TIME
				clock_gettime(CLOCK_REALTIME, &ts);
				end_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
				double cpu_time_used = end_time - start_time;
#endif
				printf("u_%d is calculated with precision %d.\n", m, precisionBits);
#ifdef TIME
				printf("\rIt took %f seconds.\n", cpu_time_used);
				start_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
#endif
			}
#endif

		mpfr_clear(sumGlobal);
	}
#ifdef TIME
	clock_gettime(CLOCK_REALTIME, &ts);
	double end_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
	printf("\nTotally, betas took %f seconds.\n", end_time - total_time);
#endif
	mpz_clear(den);
	mpfr_clear(denf);
}


bool isAdmissible(int x) {
    while (x > 0) {
        if (x % 10 == 9) {
            return false;
        }
        x /= 10;
    }
    return true;
}


ulongi powi(int base, int exp) {
    ulongi result = 1;
    for (;;) {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }
    return result;
}

// beta(l, m) = sum of 1/n^m from n = 10^(l-1) to 10^l
void beta(mpfr_t betas[LIMIT], mpfr_t betaLimit, int decr) {

#ifdef TIME
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	double total_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
#endif
	mpfr_t one; mpfr_init_set_ui(one, 1, MPFR_RNDN);
	ulongi start = powi(10, L-1);
	ulongi end = start * 10;

	#pragma omp parallel 
	{
	int precisionBits = BITS+GUARD;
#ifdef TIME
		struct timespec ts;
		clock_gettime(CLOCK_REALTIME, &ts);
		double start_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
		double end_time;
#endif
		
		#pragma omp for
		for (int m = 0; m < LIMIT; m++) {
			if (m % decr == 0) {
				int newBits = BITS + GUARD - (int)((double)m/decr) * STEP;
				precisionBits = (newBits > 0) ? newBits : precisionBits;
			}
			mpfr_t term; mpfr_init2(term, precisionBits);
			mpz_t power; mpz_init(power);

			mpfr_init2(betas[m], precisionBits);
			mpfr_set_ui(betas[m], 0, MPFR_RNDN);

			for (ulongi n = start; n < end; n++) {

				if (isAdmissible(n)) {
					mpz_ui_pow_ui(power, n, m);
					mpfr_div_z(term, one, power, MPFR_RNDN);
					int comp = mpfr_cmp(term, betaLimit);
					if (comp < 0) {
						mpfr_add(betas[m], betas[m], term, MPFR_RNDN); // betas[m] += term
						// printf("Finished betas. Calculated beta[%d].\n", m);
						break;
					}
					mpfr_add(betas[m], betas[m], term, MPFR_RNDN); // betas[m] += term
				} 
			}
			mpfr_clear(term);
			mpz_clear(power);
#ifdef PRINT	
			if (m%PRINTPERCALC == 0) {
#ifdef TIME
				clock_gettime(CLOCK_REALTIME, &ts);
				end_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
				double cpu_time_used = end_time - start_time;
#endif
				printf("beta_%d is calculated with precision %d.\n", m, precisionBits);
#ifdef TIME
				printf("It took %f seconds.\n", cpu_time_used);
				start_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
#endif
			}
#endif
		}
	}
#ifdef TIME
	clock_gettime(CLOCK_REALTIME, &ts);
	double end_time = ts.tv_sec + ts.tv_nsec / 1.0e9;
	printf("\nTotally, betas took %f seconds.\n", end_time - total_time);
#endif
	mpfr_clears(one, NULL);
}


void sommeFinie(mpfr_t value) {
	mpfr_t s1; mpfr_init2(s1, BITS+GUARD);
	mpfr_t recip; mpfr_init2(recip, BITS+GUARD);
	mpfr_t one; mpfr_init_set_ui(one, 1, MPFR_RNDN);
	mpfr_t ten; mpfr_init_set_ui(ten, 10, MPFR_RNDN);

	mpfr_set_zero(s1, 0);
	ulongi start = powi(10, L-1);
	for (ulongi n = 1; n < start; n++) {
		if (isAdmissible(n)) {
			mpfr_div_ui(recip, one, n, MPFR_RNDN);
			mpfr_add(s1, s1, recip, MPFR_RNDN);
		}
	}

	mpfr_t s2; mpfr_init2(s2, BITS+GUARD);
	mpfr_set_zero(s2, 0);

	for (ulongi n = start; n < start*10; n++) {
		if (isAdmissible(n)) {
			mpfr_div_ui(recip, ten, n, MPFR_RNDN);
			mpfr_add(s2, s2, recip, MPFR_RNDN);
		}
	}

	mpfr_add(value, s1, s2, MPFR_RNDN);
	mpfr_clears(s1, recip, one, ten, s2, NULL);
}


void sommeInfinite(mpfr_t sum, mpfr_t u[LIMIT], mpfr_t betas[LIMIT]) {
	mpfr_t temp; mpfr_init2(temp, BITS+GUARD);

	for (int m = 1; m < LIMIT-1; m++) {
		mpfr_mul(temp, u[m], betas[m+1], MPFR_RNDN); //  temp = u * b
		
		mpfr_setsign(temp, temp, m%2, MPFR_RNDN); // t = (-1)**m * t

		mpfr_add(sum, sum, temp, MPFR_RNDN); // s =+ t
	}
	mpfr_clear(temp);
} 



