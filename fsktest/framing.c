/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 *
 * Data Format Definition
 *
 * Voice1[00:07]
 * Voice1[08:15]
 * Voice1[16:23]
 * Voice1[24:27] + Voice2[00:03]
 * Voice2[04:11]
 * Voice2[12:19]
 * Voice2[20:27]
 */

#include "fsk.h"
#include "framing.h"

static uint8_t SYNCWORD[] = {0, 1, 1, 0, 0, 1, 1, 1}; // 0x67
static uint8_t DATAWORD[] = {1, 1, 1, 1, 0, 0, 1, 0}; // 0xF2

static bool frameSync; // state is sync or not sync
//
static int frameIndex; // index into circular bit buffer
static int frameMissedCount; // How many Sync have been missed
static int framesSinceLastSync; // How many bits since the last Sync
//
static long frameTotalSyncBits; // Total RX-ed bits of SyncWord
static long frameTotalSyncErrors; // Total errors in UW bits
//
static uint8_t frameBits[64];   // FRAME_BITS

void framer_create() {
    frameSync = false;
}

// Get a single bit out of an MSB-first packed byte array

static uint8_t unpackMSB(uint8_t data[], int offset, int index) {
    return (data[offset + (index >> 3)] >> (7 - (index & 0x7)) & 0x1);
}

// See if the syncWord is where it should be, to within a tolerance

static bool matchSyncWord(struct FSK *fsk, uint8_t bits[], int tolerance, int *diffCount) {
    int diff = 0;

    *diffCount = 0;

    // Start bit pointer where Sync should be
    int ibit = frameIndex;

    if (ibit >= fsk->Nbits) {
        ibit -= fsk->Nbits;
    }

    // Walk through and match bits in frame with bits of syncWord

    for (int i = 0; i < SYNC_BITS; i++) {
        if (bits[ibit] != SYNCWORD[i]) {
            diff++;
        }

        ibit++;

        if (ibit >= fsk->Nbits) {
            ibit = 0;
        }
    }

    *diffCount = diff;

    return (diff <= tolerance);
}

static void extractCodecFrame(struct FSK *fsk, uint8_t codecBytes[], uint8_t bits[]) {
    int i;

    // Initialize the 8 bytes
    for (i = 0; i < CODEC_BYTES; i++) {
        codecBytes[i] = 0;
    }

    // Skip past the Sync Word
    int index = frameIndex + SYNC_BITS;

    if (index >= fsk->Nbits) {
        index -= fsk->Nbits;
    }

    // Extract and pack first voice codec 28-bit frame, MSB first

    for (i = 0; i < CODEC_BITS; i++) {

        codecBytes[i >> 3] |= (bits[index] & 0x1) << (7 - (i & 0x7));
        index++;

        if (index >= fsk->Nbits) {
            index = 0;
        }
    }

    index = frameIndex + SYNC_BITS + CODEC_BITS;

    if (index >= fsk->Nbits) {
        index -= fsk->Nbits;
    }

    // Extract and pack second voice codec 28-bit frame, MSB first

    for (i = 0; i < CODEC_BITS; i++) {

        codecBytes[4 + (i >> 3)] |= (bits[index] & 0x1) << (7 - (i & 0x7));
        index++;

        if (index >= fsk->Nbits) {
            index = 0;
        }
    }
}

/**
 * Method to return the modem data bit frame
 *
 * This modem data frame contains the Sync Word and two Codec2 700C Frames
 *
 * @param bits is the bit array containing the modem frame bits (64 bits)
 * @param codecBytes is the 8-bytes of codec2 frame bits
 */
void getModemFrame(uint8_t bits[], uint8_t codecBytes[]) {
    int i, ibit;

    // Fill out frame with prototype
    for (i = 0; i < SYNC_BITS; i++) {
        bits[i] = SYNCWORD[i];
    }

    // Fill out first 28-bit codec2 700C block
    ibit = 0;

    for (i = 0; i < CODEC_BITS; i++) {
        bits[i + SYNC_BITS] = unpackMSB(codecBytes, 0, ibit);
        ibit++;
    }

    // Fill out second 28-bit codec2 700C block
    ibit = 0;

    for (i = 0; i < CODEC_BITS; i++) {
        bits[i + (SYNC_BITS + CODEC_BITS)] = unpackMSB(codecBytes, 4, ibit);
        ibit++;
    }
}

/**
 * Method to return the codec2 frame bytes given a modem bit frame
 *
 * Locate the Sync Word if able, and extract the two Codec2 700C Frames
 *
 * @param codecBytes is the 8-bytes of codec2 frame bytes
 * @param modemFrameBits modem frame bits
 * @return boolean true if frame has been extracted successfully
 */
bool getCodecFrame(struct FSK *fsk, uint8_t codecBytes[], uint8_t modemFrameBits[]) {
    int diffCount;

    bool frameExtracted = false;

    for (int i = 0; i < fsk->Nbits; i++) {

        frameBits[frameIndex] = modemFrameBits[i];
        frameIndex++;

        if (frameIndex >= fsk->Nbits) {
            frameIndex -= fsk->Nbits;
        }

        if (frameSync == true) {
            // Already synchronized, just wait till syncWord is back where it should be

            framesSinceLastSync++;

            // syncWord should be here. We're sunk, so deframe anyway
            if (framesSinceLastSync == fsk->Nbits) {
                framesSinceLastSync = 0;

                if (!matchSyncWord(fsk, frameBits, 1, &diffCount)) {
                    frameMissedCount++;
                } else {
                    frameMissedCount = 0;
                }

                // If we go over the miss tolerance, go into no-sync
                if (frameMissedCount > 3) {
                    frameSync = false;
                }

                extractCodecFrame(fsk, codecBytes, frameBits);
                frameExtracted = true;

                frameTotalSyncBits += (long) SYNC_BITS;
                frameTotalSyncErrors += (long) diffCount;
            }
        } else if (matchSyncWord(fsk, frameBits, 0, &diffCount)) {
            // found sync

            frameSync = true;
            framesSinceLastSync = 0;
            frameMissedCount = 0;

            extractCodecFrame(fsk, codecBytes, frameBits);
            frameExtracted = true;

            frameTotalSyncBits += (long) SYNC_BITS;
            frameTotalSyncErrors += (long) diffCount;
        }
    }

    return frameExtracted;
}

int getFrameSizeInBits(struct FSK *fsk) {
    return fsk->Nbits;
}

int getVoiceCodecSizeInBytes() {
    return CODEC_BYTES; // 4 bytes per codec2 700 bit/s frame
}

bool IsFrameSync() {
    return frameSync;
}

long getFrameTotalSyncBits() {
    return frameTotalSyncBits;
}

long getFrameTotalSyncErrors() {
    return frameTotalSyncErrors;
}
