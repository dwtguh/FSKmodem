/*
 * Copyright (C) 1993-2017 David Rowe, Brady O'Brien
 *
 * All rights reserved
 *
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>

#define DATA_BITS 56
#define SYNC_BITS 8
#define CODEC_BITS 28
#define CODEC_BYTES 8

/**
 * Method to return the modem data bit frame
 *
 * This modem data frame contains the Sync Word and two Codec2 700C Frames
 *
 * @param bits is the bit array containing the modem frame bits
 * @param codecBytes is the 8-bytes of codec2 frame bytes
 */
void getModemFrame(uint8_t [], uint8_t []);

/**
 * Method to return the codec2 frame bytes given a modem bit frame
 *
 * Locate the Sync Word if able, and extract the two Codec2 700C Frames
 *
 * @param codecBytes is the 8-bytes of codec2 frame bytes
 * @param modemFrameBits modem frame bits
 * @return boolean true if frame has been extracted successfully
 */
bool getCodecFrame(struct FSK *, uint8_t [], uint8_t[]);

void framer_create(void);

int getFrameSizeInBits(struct FSK *);

int getVoiceCodecSizeInBytes(void);

bool IsFrameSync(void);

long getFrameTotalSyncBits(void);

long getFrameTotalSyncErrors(void);


#ifdef __cplusplus
}
#endif
