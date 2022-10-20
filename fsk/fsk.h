/*
 * Copyright (c) 2015-2021 David Rowe
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 *
 * 2015 - David Rowe designed modem in Octave
 * 2016 - Brady O'Brien converted code from Octave to C
 * 2021 - Steve Sampson redesigned demod for FFTW3 and real FFT
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#define FALSE       (0)
#define TRUE        (1)

#define MODE_2FSK   2
#define MODE_4FSK   4
#define MODE_M_MAX  4

#define MODE_BURST  0
#define MODE_STREAM 1

#define FDMDV_SCALE 750
    
/* default internal parameters */

#define FSK_DEFAULT_P    10     // Number of timing offsets, try to keep P >= 8
#define FSK_DEFAULT_NSYM 50     // See Nsym below
#define FSK_NONE         -1     // unused parameter

#ifndef M_PI
#define M_PI    3.14159265358979f
#endif

#ifndef TAU
#define TAU     (2.0f * M_PI)
#endif

#define cmplx(value) (cosf(value) + sinf(value) * I)
#define cmplxconj(value) (cosf(value) + sinf(value) * -I)
    
struct STATS {
    float f_est[MODE_M_MAX];
    float snr_est; // estimated SNR of rx signal in dB (3 kHz noise BW)
    float foff; // estimated freq offset in Hz
    float rx_timing; // estimated optimum timing offset in samples
    float clock_offset; // Estimated tx/rx sample clock offset in ppm
};

struct FSK {
    fftwf_plan plan;
    fftwf_complex *fftout;
    fftwf_complex *f_dc; // down converted samples

    complex float phi_c[MODE_M_MAX]; /* phase of each demod local oscillator */
    complex float tx_phase_c; /* TX phase, but complex */

    float *fftin;
    float *Sf; /* Average of magnitude spectrum */    
    float norm_rx_timing; /* Normalized RX timing */
    float EbNodB; /* Estimated EbNo in dB */
    float f_est[MODE_M_MAX]; /* Estimated frequencies (peak method) */
    float ppm; /* Estimated PPM clock offset */
    float SNRest; /* used for LLRs */
    float v_est; /* used for LLRs */
    float rx_sig_pow;
    float rx_nse_pow;

    int Ndft; /* freq offset est fft */
    int Fs; /* sample freq */
    int N; /* processing buffer size */
    int Rs; /* symbol rate */
    int Ts; /* samples per symbol */
    int Nmem; /* size of extra mem for timing adj */
    int P; /* oversample rate for timing est/adj */
    int Nsym; /* Number of symbols in averaging window */
    int Nbits; /* Number of bits in a processing frame */
    int f1_tx; /* f1 for modulator */
    int tone_spacing; /* Space between TX freqs for modulator (and option mask freq estimator) */
    int mode; /* 2FSK or 4FSK */
    int est_min; /* Minimum frequency for freq. estimator */
    int est_max; /* Maximum frequency for freq. estimator */
    int est_space; /* Minimum frequency spacing for freq. estimator */
    int nin; /* Number of samples to feed the next demod cycle */
    int stream_mode; /* stream: 1 burst: 0 */

    struct STATS *stats;
};

int fsk_get_error_num(void);

/*
 * Create a FSK modem
 *
 * int Fs - Sample frequency
 * int Rs - Symbol rate
 * int M  - 2 for 2FSK, 4 for 4FSK
 * int f1_tx - first tone frequency
 * int tone_spacing - frequency spacing (for modulator and optional "mask" freq estimator)
 */
struct FSK * fsk_create(int, int, int, int, int);

/*
 * Create a FSK modem - advanced version
 *
 * int Fs - Sample frequency
 * int Rs - Symbol rate
 * int M  - 2 for 2FSK, 4 for 4FSK
 * int P  - number of timing offsets to choose from (suggest >= 8)
 * int Nsym  - windows size for timing estimator
 * int f1_tx - first tone frequency
 * int tone_spacing - frequency spacing (for modulator and optional "mask" freq estimator)
 */
struct FSK * fsk_create_hbr(int, int, int, int, int, int, int);

/*
 * Clear the estimator states
 */
void fsk_clear_estimators(struct FSK *);

/*
 * Fills MODEM_STATS struct with demod statistics
 */
void fsk_get_demod_stats(struct FSK *, struct STATS *);

/*
 * Destroy an FSK state struct and free it's memory
 */
void fsk_destroy(struct FSK *);

/*
 * Modulates Nsym bits into N samples
 *
 * struct FSK *fsk   - FSK config/state struct, set up by fsk_create
 * float fsk_out[]   - Buffer for samples of modulated FSK, fsk->Ts*(nbits/(M>>1)) in length
 * uint8_t tx_bits[] - Buffer containing Nbits unpacked bits
 * int     nbits     - number of bits to transmit
 */
void fsk_mod(struct FSK *, float [], uint8_t [], int);

/*
 * Modulates Nsym bits into N complex samples
 *
 * struct FSK *fsk   - FSK config/state struct, set up by fsk_create
 * comp fsk_out[]    - Buffer for samples of modulated FSK, fsk->Ts*(nbits/(M>>1)) in length
 * uint8_t tx_bits[] - Buffer containing Nbits unpacked bits
 * int     nbits     - number of bits to transmit
 */
void fsk_mod_c(struct FSK *, complex float [], uint8_t [], int);

/*
 * Returns the number of samples needed for the next fsk_demod() cycle
 *
 * struct FSK *fsk - FSK config/state struct, set up by fsk_create
 * returns - number of samples to be fed into fsk_demod next cycle
 */
int fsk_nin(struct FSK *);

/*
 * Demodulate some number of FSK samples. The number of samples to be
 *  demodulated can be found by calling fsk_nin().
 *
 * struct FSK *fsk   - FSK config/state struct, set up by fsk_create
 * uint8_t rx_bits[] - Buffer for fsk->Nbits unpacked bits to be written
 * complex fsk_in[]    - nin samples of modulated FSK
 */
void fsk_demod(struct FSK *, uint8_t [], complex float []);

/*
 * Soft decision demodulation
 *
 * struct FSK *fsk - FSK config/state struct, set up by fsk_create
 * float rx_flit[] - M x Nsym array of filtermagnitude outputs
 * complex fsk_in[]  - nin samples of modualted FSK
 */
void fsk_demod_sd(struct FSK *, float [], complex float []);

/* Set the FSK modem in stream or burst demod mode */

void fsk_set_stream_mode(struct FSK *, int);

/*
 * Set the minimum and maximum frequencies at which the freq. estimator can find tones
 */
void fsk_set_freq_est_limits(struct FSK *, int, int);

#ifdef __cplusplus
}
#endif
