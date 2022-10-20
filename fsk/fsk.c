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
 * 2021 - Steve Sampson modified for use with complex syntax and FFTW3
 */

#include "fsk.h"

// imbed in binary
char const *copyright = "Copyright (c) 2015-2021 David Rowe, All rights reserved";

static float cnormf(fftwf_complex val) {
    float realf = crealf(val);
    float imagf = cimagf(val);

    return realf * realf + imagf * imagf;
}

static void fsk_demod_freq_est(struct FSK *, complex float []);
static void fsk_demod_core(struct FSK *, uint8_t [], float [], complex float []);
static struct FSK *fsk_create_core(int, int, int, int, int, int, int);

static char const *Wisdomf_file = "/etc/fftw/wisdomf";

static int error_num;

static struct FSK *fsk_create_core(int Fs, int Rs, int M, int P, int Nsym, int f1_tx, int tone_spacing) {
    struct FSK *fsk;

    /* Ts (Fs/Rs) must be an integer */
    /* Ts/P (Fs/Rs/P) must be an integer */
    /* P must be > 8 to have good choice of timing offsets to choose from */

    error_num = 0;

    fsk = (struct FSK *) calloc(1, sizeof (struct FSK));

    if (fsk == NULL) {
        error_num = 1;
        return NULL;
    }

    // Check for bogus mode

    switch (M) {
        case MODE_2FSK:
        case MODE_4FSK:
            fsk->mode = M;
            break;
        default:
            free(fsk);
            error_num = 2;
            return NULL;
    }

    // Need enough bins to within 10% of tone center

    int Ndft = powf(2.0f, ceilf(log2f((float) Fs / (0.1f * Rs)))); // 64 = 9600 fs 2400 rs 2FSK

    fsk->Fs = Fs;                           // Ex: 9600
    fsk->Rs = Rs;                           // Ex: 2400
    fsk->Ts = Fs / Rs;                      // Ex: 4 samples per symbol
    fsk->stream_mode = FALSE;
    fsk->P = P;                             // Ex: 10
    fsk->Nsym = Nsym;                       // Ex: 64 bits or symbols (8 bytes)
    fsk->N = fsk->Ts * fsk->Nsym;           // Ex: 256 samples
    fsk->Ndft = Ndft;                       // Ex: 64
    fsk->Nmem = fsk->N + (2 * fsk->Ts);     // Ex: 256 + 8 = 264
    fsk->f1_tx = f1_tx;                     // 900 Hz
    fsk->tone_spacing = tone_spacing;       // 2400 Hz
    fsk->nin = fsk->N;                      // Ex: 256
    fsk->Nbits = (fsk->mode == MODE_2FSK) ? fsk->Nsym : fsk->Nsym * 2;    // Ex: 2048 or 4096
    fsk->est_min = 0;
    fsk->est_max = Fs/2;                    // Nyquist 4800
    fsk->est_space = 0.75f * Rs;            // Ex: 1800 (Minimum frequency spacing)
    fsk->tx_phase_c = cmplx(0.0f);
    fsk->EbNodB = 0.0f;
    fsk->ppm = 0.0f;
    fsk->norm_rx_timing = 0;

    /* Set up rx state for each freq */

    for (int i = 0; i < fsk->mode; i++) {
        fsk->phi_c[i] = cmplx(0.0f);
        fsk->f_est[i] = 0.0f;
    }

    fsk->f_dc = (complex float *) calloc(fsk->mode * fsk->Nmem, sizeof (complex float));

    if (fsk->f_dc == NULL) {
        error_num = 3;
        return NULL;
    }

    fsk->Sf = (float *) calloc(Ndft, sizeof (float));

    if (fsk->Sf == NULL) {
        error_num = 4;
        return NULL;
    }

    fsk->stats = (struct STATS *) calloc(1, sizeof (struct STATS));

    if (fsk->stats == NULL) {
        error_num = 5;
        return NULL;
    }

    /*
     * Caution: if you compile FFTW for floats, then the file
     * we are looking for here is called /etc/fftw/wisdomf
     *
     * and the header must be like:
     * (fftw-3.3.10 fftwf_wisdom #xb8b78987 #x48270938 #xff2d15a5 #xefab3157)
     *
     * That is: PACKAGE "-" VERSION " " "fftwf_wisdom"
     *
     * Followed by four hashes and a /n
     *
     */
    int r = fftwf_import_system_wisdom();

    if (r != 1) {
        error_num = 6;
        return NULL;
    }

    r = fftwf_import_wisdom_from_filename(Wisdomf_file);

    if (r != 1) {
        error_num = 7;
        return NULL;
    }

    fsk->fftin = fftwf_malloc(sizeof(float) * fsk->Ndft);
    fsk->fftout = fftwf_malloc(sizeof(fftwf_complex) * fsk->Ndft);

    fsk->plan = fftwf_plan_dft_r2c_1d(fsk->Ndft, fsk->fftin, fsk->fftout, FFTW_ESTIMATE);

    return fsk;
}

int fsk_get_error_num() {
    return error_num;
}

/*---------------------------------------------------------------------------*\

  Call to free all memory.

\*---------------------------------------------------------------------------*/

void fsk_destroy(struct FSK *fsk) {
    fftwf_destroy_plan(fsk->plan);
    fftwf_free(fsk->fftout);
    fftwf_free(fsk->fftin);

    free(fsk->stats);
    free(fsk->Sf);
    free(fsk->f_dc);
    free(fsk);
}

/*---------------------------------------------------------------------------* \

  Create and initialize an instance of the FSK modem. Returns a pointer
  to the modem state/config struct. One modem config struct is used
  for both mod and demod.

  If you are intending to use only the demodulation functions, you can
  set f1_tx to FSK_NONE.

\*---------------------------------------------------------------------------*/

struct FSK *fsk_create(int Fs, int Rs, int M, int f1_tx, int tone_spacing) {
    return fsk_create_core(Fs, Rs, M, FSK_DEFAULT_P, FSK_DEFAULT_NSYM, f1_tx, tone_spacing);
}

/*---------------------------------------------------------------------------*\

  Alternate version of create allows user defined oversampling P and
  averaging window Nsym.  In the current version of the demod it's
  simply an alias for the default core function.

  P is the oversampling rate of the internal demod processing, which
  happens at Rs*P Hz.  We filter the tones at P different timing
  offsets, and choose the best one.

  P should be >=8, so we have a choice of at least 8 timing offsets.
  This may require some adjustment of Fs and Rs, as Fs/Rs/P must be
  an integer.

  Nsym is the number of symbols we average over demod parameters
  like the symbol timing.

\*---------------------------------------------------------------------------*/

struct FSK *fsk_create_hbr(int Fs, int Rs, int M, int P, int Nsym, int f1_tx, int tone_spacing) {
    return fsk_create_core(Fs, Rs, M, P, Nsym, f1_tx, tone_spacing);
}

/*---------------------------------------------------------------------------*\

  FSK modulator function, real output samples

\*---------------------------------------------------------------------------*/

void fsk_mod(struct FSK *fsk, float fsk_out[], uint8_t tx_bits[], int nbits) {
    complex float tx_phase_c = fsk->tx_phase_c; /* Current complex TX phase */
    complex float dosc_f[fsk->mode];        /* phase shift per sample */

    /* Init the per sample phase shift complex numbers */
    for (int m = 0; m < fsk->mode; m++) {
        dosc_f[m] = cmplx(TAU * ((float) (fsk->f1_tx + (fsk->tone_spacing * m)) / (float) fsk->Fs));
    }

    int bit_i = 0;
    int nsym = (fsk->mode == MODE_4FSK) ? nbits / 2 : nbits;

    for (int i = 0; i < nsym; i++) {
        int ibit = 0;

        if (fsk->mode == MODE_4FSK) {
            ibit = tx_bits[bit_i];
            bit_i++;
            ibit = (ibit << 1) | tx_bits[bit_i];
        } else {
            ibit = tx_bits[bit_i];
        }

        bit_i++;

        /* Look up symbol phase shift */
        complex float dph = dosc_f[ibit];

        /* Spin the oscillator for a symbol period */

        for (int j = 0; j < fsk->Ts; j++) {
            tx_phase_c *= dph;

            fsk_out[i * fsk->Ts + j] = crealf(tx_phase_c) * 2.0f;
        }
    }

    /* Normalize TX phase to prevent drift */
    tx_phase_c /= cabsf(tx_phase_c);

    /* save TX phase */
    fsk->tx_phase_c = tx_phase_c;
}

/*---------------------------------------------------------------------------*\

  FSK modulator function, complex output samples

\*---------------------------------------------------------------------------*/

void fsk_mod_c(struct FSK *fsk, complex float fsk_out[], uint8_t tx_bits[], int nbits) {
    complex float tx_phase_c = fsk->tx_phase_c; /* Current complex TX phase */
    complex float dosc_f[fsk->mode];            /* phase shift per sample */

    /* Init the per sample phase shift complex numbers */
    for (int m = 0; m < fsk->mode; m++) {
        dosc_f[m] = cmplx(TAU * ((float) (fsk->f1_tx + (fsk->tone_spacing * m)) / (float) fsk->Fs));
    }

    int bit_i = 0;
    int nsym = (fsk->mode == MODE_4FSK) ? nbits / 2 : nbits;

    for (int i = 0; i < nsym; i++) {
        int ibit = 0;

        if (fsk->mode == MODE_4FSK) {
            ibit = tx_bits[bit_i];
            bit_i++;

            ibit = (ibit << 1) | tx_bits[bit_i];
        } else {
            ibit = tx_bits[bit_i];
        }

        bit_i++;

        /* Look up symbol phase shift */
        complex float dph = dosc_f[ibit];

        /* Spin the oscillator for a symbol period */

        for (int j = 0; j < fsk->Ts; j++) {
            tx_phase_c *= dph;

            fsk_out[i * fsk->Ts + j] = tx_phase_c * 2.0f;
        }
    }

    /* Normalize TX phase to prevent drift */
    tx_phase_c /= cabsf(tx_phase_c);

    /* save TX phase */
    fsk->tx_phase_c = tx_phase_c;
}

/*---------------------------------------------------------------------------*\

  Call before each call to fsk_demod() to determine how many new
  samples you should pass in.  the number of samples will vary due to
  timing variations.

\*---------------------------------------------------------------------------*/

int fsk_nin(struct FSK *fsk) {
    return fsk->nin;
}

/*
 * Internal function to estimate the frequencies of the FSK tones.
 * This is split off because it is fairly complicated, needs
 * a bunch of memory, and probably takes more cycles than the rest
 * of the demod.
 *
 * Parameters:
 *
 * fsk - FSK struct from demod containing FSK config
 * fsk_in - block of samples in this demod cycles, must be nin long
 */
static void fsk_demod_freq_est(struct FSK *fsk, complex float fsk_in[]) {
    int Ndft2 = fsk->Ndft / 2;
    
    int fmask[Ndft2];

    int chunks = floorf((float) fsk->nin / Ndft2) - 1; // Ex: floor(8192/(64/2))-1 = 255
    int dftn = fsk->Ndft - 1;
    
    for (int j = 0; j < chunks; j++) {
        int offset = j * (fsk->Ndft / 2);     // Ex: 32, 64, 96...

        /*
         * Copy FSK buffer into FFT buffer and apply a hann window
         *
         * Each chunk starts reading the samples by an (Ndft / 2) Ex: 32 offset
         */

        for (int i = 0; i < fsk->Ndft; i++) {
            float hann = 0.5f - 0.5f * cosf(TAU * i / dftn);
            
            // Don't fly off the end
            if (i + offset < fsk->nin) {
                fsk->fftin[i] = fsk_in[i + offset] * hann;
            } else {
                fsk->fftin[i] = 0.0f;
            }
        }

        /* Do the real to complex FFT */

        fftwf_execute(fsk->plan);

        /* Find the mag^2 of each real FFT slot */

        for (int i = 0; i < Ndft2; i++) {
            float real_part = sqrtf(cnormf(fsk->fftout[i]));

            fsk->Sf[i] = (fsk->Sf[i] * .9f) + (real_part * .1f);

            // save a integer copy in fmask for peak detection below
            fmask[i] = (int) fsk->Sf[i];
        }
    }
    
    int freqi[fsk->mode];

    int f_bin = (fsk->est_max / Ndft2);  // Ex: 32 bins of 150 Hz = 4800 / 32 = 150
    
    int f_width = (int) ceilf((float) f_bin / ((float) fsk->est_space / Ndft2));  // Ex: 3

    /* Find the frequency peaks here */
        
    for (int i = 0; i < fsk->mode; i++) {
        int imax = 0;
        int max = 0;

        // Search the mask for a max
        // found peaks will be zeroed below
        // to remove from list while searching

        for (int j = 0; j < Ndft2; j++) {
            if (fmask[j] > max) {
                max = fmask[j];
                imax = j;
            }
        }

        // Blank out magnitude value Ex: +/- 450 Hz (+/- 3 bins

        int f_min = imax - f_width;              // Detect bin -3
        f_min = (f_min < 0) ? 0 : f_min;

        int f_max = imax + f_width;              // Detect bin + 3
        f_max = (f_max > Ndft2) ? Ndft2 : f_max;

        for (int j = f_min; j < f_max; j++) {
            fmask[j] = 0;
        }

        freqi[i] = imax;
    }

    // Insertion sort

    int ic = 1;
    while (ic < fsk->mode) {
        int x = freqi[ic];
        int jc = ic - 1;

        while (jc >= 0 && freqi[jc] > x) {
            freqi[jc + 1] = freqi[jc];
            jc--;
        }

        freqi[jc + 1] = x;
        ic++;
    }

    // Convert bin to frequency
    
    for (int i = 0; i < fsk->mode; i++) {
        fsk->f_est[i] = (float) freqi[i] * f_bin;
    }
}

/*
 * core demodulator function
 */
static void fsk_demod_core(struct FSK *fsk, uint8_t rx_bits[], float rx_filt[], complex float fsk_in[]) {
    /*
     * Estimate tone frequencies
     */
    fsk_demod_freq_est(fsk, fsk_in);

    int nold = (fsk->Nmem - fsk->nin);  // Ex: 8

    /* update filter (integrator) memory by shifting in nin samples */
    for (int m = 0; m < fsk->mode; m++) {
        for (int i = 0, j = (fsk->Nmem - nold); i < nold; i++, j++) {
            fsk->f_dc[m * fsk->Nmem + i] = fsk->f_dc[m * fsk->Nmem + j];
        }
    }

    /* freq shift down to around DC, ensuring continuous phase from last frame */

    for (int m = 0; m < fsk->mode; m++) {
        complex float dphi_m = cmplx(TAU * ((float) fsk->f_est[m] / (float) (fsk->Fs / 2)));

        for (int i = nold, j = 0; i < fsk->Nmem; i++, j++) {
            fsk->phi_c[m] *= dphi_m;

            fsk->f_dc[m * fsk->Nmem + i] = (fftwf_complex) (fsk_in[i] * conjf(fsk->phi_c[m]));
        }

        fsk->phi_c[m] /= cabsf(fsk->phi_c[m]);
    }

    // integrate over symbol period at a variety of offsets
    
    fftwf_complex f_int[fsk->mode][(fsk->Nsym + 1) * fsk->P];

    for (int i = 0; i < (fsk->Nsym + 1) * fsk->P; i++) {
        int st = (i * fsk->Ts) / fsk->P;
        int en = st + (fsk->Ts - 1);

        for (int m = 0; m < fsk->mode; m++) {
            f_int[m][i] = 0.0f;

            for (int j = st; j <= en; j++) {
                f_int[m][i] += fsk->f_dc[m * fsk->Nmem + j];
            }
        }
    }

    /*
     * Fine Timing Estimation
     */

    /*
     * Apply magic nonlinearity to f1_int and f2_int, shift down to 0,
     * extract angle, figure out how much to spin the oscillator to
     * extract magic spectral line
     */

    complex float dphift = cmplx(TAU * ((float) fsk->Rs / (float) (fsk->P * fsk->Rs)));

    complex float phi_ft = cmplx(0.0f);
    complex float t_c = 0.0f;

    for (int i = 0; i < (fsk->Nsym + 1) * fsk->P; i++) {
        /* Get abs^2 of fx_int[i], and add 'em */
        float absq = 0.0f;

        for (int m = 0; m < fsk->mode; m++) {
            absq += cnormf(f_int[m][i]);
        }

        /* Down shift and accumulate magic line */
        t_c += (phi_ft * absq);

        /* Spin the oscillator for the magic line shift */
        phi_ft *= dphift;
    }

    /* Check for NaNs in the fine timing estimate, return if found */
    /* otherwise segfaults happen */

    if (isnan(crealf(t_c)) || isnan(cimagf(t_c))) {
        return;
    }

    /* Get the magic angle */
    float norm_rx_timing = cargf(t_c) / TAU;

    float rx_timing = norm_rx_timing * (float) fsk->P;
    
    float old_norm_rx_timing = fsk->norm_rx_timing;
    fsk->norm_rx_timing = norm_rx_timing;

    /* Estimate sample clock offset */
    float d_norm_rx_timing = norm_rx_timing - old_norm_rx_timing;

    /* Filter out big jumps in due to nin change */
    if (fabsf(d_norm_rx_timing) < .2f) {
        float appm = 1e6f * d_norm_rx_timing / (float) fsk->Nsym;
        fsk->ppm = .9f * fsk->ppm + .1f * appm;
    }

    /*
     * Figure out how many samples are needed the next modem cycle
     *
     * In streaming mode, keep nin the same value. In burst packet
     * mode, adjust as we go.
     */

    if (fsk->stream_mode == FALSE) {
        // we're in burst mode, so adjust

        if (norm_rx_timing > 0.25f)
            fsk->nin = fsk->N + fsk->Ts / 4;
        else if (norm_rx_timing < -0.25f)
            fsk->nin = fsk->N - fsk->Ts / 4;
        else
            fsk->nin = fsk->N;
    }

    /* Re-sample the integrators with linear interpolation magic */

    int low_sample = (int) floorf(rx_timing);

    float fract = rx_timing - (float) low_sample;

    int high_sample = (int) ceilf(rx_timing);

    /* Vars for finding the max-of-4 for each bit */

    float tmax[fsk->mode];
    fftwf_complex t[fsk->mode]; /* complex number temps */

    float meanebno = 0.0f;
    float stdebno = 0.0f;

    float rx_nse_pow = 1E-12f;
    float rx_sig_pow = 0.0f;

    for (int i = 0; i < fsk->Nsym; i++) {
        /* resample at ideal sampling instant */
        int st = (i + 1) * fsk->P;

        for (int m = 0; m < fsk->mode; m++) {
            t[m] = (f_int[m][st + low_sample] * (1.0f - fract));
            t[m] += (f_int[m][st + high_sample] * fract);

            /* Figure mag^2 of each resampled fx_int */
            tmax[m] = cnormf(t[m]);
        }

        /* hard decision decoding of bits */
        float max = tmax[0];
        float min = tmax[0];

        int sym = 0; /* Index of maximum */

        for (int m = 0; m < fsk->mode; m++) {
            if (tmax[m] > max) {
                max = tmax[m];
                sym = m;
            }

            if (tmax[m] < min) {
                min = tmax[m];
            }
        }

        /*
         * Caller may not want any bits returned
         */
        if (rx_bits != NULL) {
            if (fsk->mode == MODE_2FSK) {
                rx_bits[i] = (sym == 1); /* 1 bit per symbol */
            } else if (fsk->mode == MODE_4FSK) {
                rx_bits[(i * 2)] = (sym >> 1) & 0x1; /* LSB */
                rx_bits[(i * 2) + 1] = sym & 0x1;    /* MSB 2 bit/symbol */
            }
        }

        /* Optionally output filter magnitudes for soft decision/LLR
           calculation.  Update SNRest always as this is a useful
           alternative to the earlier EbNo estimator below */
        float sum = 0.0f;

        for (int m = 0; m < fsk->mode; m++) {
            if (rx_filt != NULL)
                rx_filt[m * fsk->Nsym + i] = sqrtf(tmax[m]);

            sum += tmax[m];
        }

        rx_sig_pow += max;
        rx_nse_pow += (sum - max) / (fsk->mode - 1);

        /* Accumulate resampled int magnitude for EbNodB estimation */

        /* Accumulate the square of the sampled value */
        stdebno += max;

        /* Figure the abs value of the max tone */
        meanebno += sqrtf(max);
    }

    fsk->rx_sig_pow = rx_sig_pow = rx_sig_pow / fsk->Nsym;
    fsk->rx_nse_pow = rx_nse_pow = rx_nse_pow / fsk->Nsym;

    fsk->v_est = sqrtf(rx_sig_pow - rx_nse_pow);
    fsk->SNRest = rx_sig_pow / rx_nse_pow;

    /* Calculate mean for EbNodB estimation */
    meanebno = meanebno / (float) fsk->Nsym;

    /* Calculate the std. dev for EbNodB estimate */
    stdebno = (stdebno / (float) fsk->Nsym) - (meanebno * meanebno);

    /* trap any negative numbers to avoid NANs flowing through */
    if (stdebno > 0.0f) {
        stdebno = sqrtf(stdebno);
    } else {
        stdebno = 0.0f;
    }

    fsk->EbNodB = -6.0f + (20.0f * log10f((1e-6f + meanebno) / (1e-6f + stdebno)));

    /* Calculate and save SNR from EbNodB estimate */

    fsk->stats->snr_est = .5f * fsk->stats->snr_est + .5f * fsk->EbNodB;

    /* Save rx timing */
    fsk->stats->rx_timing = (float) rx_timing;

    /* Estimate and save frequency offset */
    float fc_avg = 0.0f;
    float fc_tx = 0.0f;

    for (int m = 0; m < fsk->mode; m++) {
        fc_avg += fsk->f_est[m];    // Ex: 900 + 3300
        fc_tx += (float) fsk->f1_tx + (float) (m * fsk->tone_spacing);  // Ex: 900 + 3300
    }

    fc_avg /= (float) fsk->mode;    // Ex: 2100
    fc_tx /= (float) fsk->mode;     // Ex: 2100

    fsk->stats->foff = fc_tx - fc_avg;

    for (int i = 0; i < fsk->mode; i++) {
        fsk->stats->f_est[i] = fsk->f_est[i];
    }
}

/*---------------------------------------------------------------------------*\

  FSK demodulator utility functions:

    fsk_demod...: complex samples in, bits out
    fsk_demos_sd: complex samples in, soft decision symbols out

\*---------------------------------------------------------------------------*/

void fsk_demod(struct FSK *fsk, uint8_t rx_bits[], complex float fsk_in[]) {
    fsk_demod_core(fsk, rx_bits, NULL, fsk_in);
}

void fsk_demod_sd(struct FSK *fsk, float rx_filt[], complex float fsk_in[]) {
    fsk_demod_core(fsk, NULL, rx_filt, fsk_in);
}

/* Set the FSK modem in stream or burst demod mode */

void fsk_set_stream_mode(struct FSK *fsk, int val) {
    fsk->nin = fsk->N;
    fsk->stream_mode = val;
}

/*
 * Clear freq estimator state
 */
void fsk_clear_estimators(struct FSK *fsk) {
    memset(fsk->Sf, 0, sizeof(float) * fsk->Ndft);
    fsk->nin = fsk->N;
}

void fsk_get_demod_stats(struct FSK *fsk, struct STATS *stats) {
    stats->clock_offset = fsk->ppm;
    stats->snr_est = fsk->stats->snr_est; // TODO: make this SNR not Eb/No
    stats->rx_timing = fsk->stats->rx_timing;
    stats->foff = fsk->stats->foff;
}

/*
 * Set the minimum and maximum frequencies at which the freq. estimator can find tones
 */
void fsk_set_freq_est_limits(struct FSK *fsk, int est_min, int est_max) {
    fsk->est_min = est_min;
    fsk->est_max = est_max;
}
