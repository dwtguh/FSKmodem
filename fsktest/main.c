/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 *
 * 2FSK engine test driver by S. Sampson 2021
 */

#include <stdio.h>
#include <complex.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <fsk.h>

#include "framing.h"

/*
 * return a uniform random value
 * in the range 0..n-1 inclusive
 */
static int randRange(int n) {
    int limit = RAND_MAX - (RAND_MAX % n);

    int r;

    while ((r = rand()) >= limit);

    return r % n;
}

static void outputSpectrum(float *data, int datalen, int loop) {
    printf("spectrum=[");

    for (int i = 0; i < datalen; i++) {
        printf("%g ", data[i]);
    }

    printf("];figure(\"Position\",[%d 200 560 420]);plot(spectrum);\n", 30 * loop);
}

/*
 * Simulate a 2FSK system with a 9600 sample rate and a 2400 symbol rate
 *
 * The frequency pairs are 900 Hz and 900 + symbol rate or 3300 Hz
 */
int main(int argc, char** argv) {
    struct FSK *fsk;
    FILE *fin;
    FILE *fout;
    bool plot = false;

    srand(time(0));

    if (argc == 2) {
        if (strcmp(argv[1], "-p") == 0) {   // plot
            plot = true;
        }
    }

    char *fname = "/tmp/mod1.raw";

    if ((fout = fopen(fname, "wb")) == NULL) {
        fprintf(stderr, "Error opening output modem sample file: %s\n", fname);
        exit(-1);
    }

    // Check to see that the framer works

    framer_create();

    int Fs = 9600;
    int Rs = 2400;
    int txfreq = 900;
    int txsep = 2400;
    int nsym = 64;

    // Instantiate our modem
    fsk = fsk_create_hbr(Fs, Rs, MODE_2FSK, FSK_DEFAULT_P, nsym, txfreq, txsep);

    if (fsk == NULL) {
        fsk_destroy(fsk);

        int error_num = fsk_get_error_num();
        fprintf(stderr, "Unable to instantiate FSK instance error_num = %d\n", error_num);

        exit(error_num);
    }

    /*
     * Transmit a packet stream rather than burst
     * So this sets the demodulator appropriately
     */
    fsk_set_stream_mode(fsk, MODE_STREAM);

    int16_t txbuf[fsk->N]; // Ex: (Fs / Rs) * nsym = 256
    float modulation[fsk->N];

    uint8_t codecData[fsk->Nbits / 8];
    uint8_t codecBits[fsk->Nbits];  // Codec frame is 64 bits with sync

    // Now send the real data

    for (int loop = 0; loop < 8; loop++) { // loop to give a good quantity of data
        for (int i = 0; i < (fsk->Nbits / 8); i++) {
            codecData[i] = (uint8_t) randRange(256); // create a simulated codec2 frame
        }

        getModemFrame(codecBits, codecData);

        // Modulate bits and write to a file

        fsk_mod(fsk, modulation, codecBits, fsk->Nbits);

        // Store the modem cycle audio in little-endian raw signed 16-bit PCM

        for (int i = 0; i < fsk->N; i++) {
            txbuf[i] = (int16_t) (modulation[i] * 8192.0f);
        }

        // write the PCM signed 16-bit samples to a file

        fwrite(txbuf, 2, fsk->N, fout);
    }

    fclose(fout);

    // open the test file for decoding

    if ((fin = fopen(fname, "rb")) == NULL) {
        fprintf(stderr, "Error opening input data file: %s\n", fname);
        exit(-1);
    }

    int16_t *rx_pcmbuf = (int16_t*) calloc(4 * (fsk->N + (fsk->Ts * 2)), sizeof (int16_t));
    uint8_t bits[fsk->Nsym];
    uint8_t frame[fsk->Nbits / 8];

    // Spectrum storage is sized to breathe a bit N + (samples per symbol * 2)
    // Generally it will use about nin samples on each call.

    complex float *spectrum = (complex float *) calloc((fsk->N + (fsk->Ts * 2)), sizeof (complex float));

    int loop = 0;

    while (fread(rx_pcmbuf, 2, fsk_nin(fsk), fin) == fsk_nin(fsk)) {
        for (int i = 0; i < fsk_nin(fsk); i++) {
            spectrum[i] = (((float) rx_pcmbuf[i]) / FDMDV_SCALE) + 0.0f * I;
        }

        fsk_demod(fsk, bits, spectrum);

        bool goodSync = getCodecFrame(fsk, frame, bits);

        // Plot Octave spectrum via pipe if plot boolean is true

        if (plot == true) {
            outputSpectrum(fsk->Sf, fsk->Ndft / 2, loop);
        } else {
            printf("\n--------------------------\n");
            printf("Frequency offset: %.4f\nPPM = %.4f\n\n", fsk->stats->foff, fsk->ppm);

            for (int i = 0; i < fsk->mode; i++) { // convert to Hz
                printf("Frequency %d: %.2f\n", i + 1, fsk->f_est[i]);
            }

            // We have bits now
            for (int i = 0; i < (fsk->Nbits); i++) {
                printf(frame[i] ? "1" : "0");
            }

            printf("\n\n%s Sync\n", goodSync == true ? "Good" : "Bad");

            printf("SNR = %.2f\n", fsk->stats->snr_est);
        }

        loop++;
    }

    fclose(fin);

    free(rx_pcmbuf);
    free(spectrum);

    fsk_destroy(fsk);

    fflush(stdout);

    // wait for key-press if not in plot mode

    if (plot == true) {
        getchar();
    }

    return (EXIT_SUCCESS);
}
