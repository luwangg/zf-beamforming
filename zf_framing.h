#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <gnuradio/gr_complex.h>
#include <liquid/liquid.h>
#include <assert.h>
#include <fftw3.h>
#include <string.h>
#include <vector>

#define USE_DEFAULT_MS  true

typedef int (*zf_framesync_callback)(gr_complex * _y,
                                     unsigned char * _p,
                                     unsigned int _M,
                                     void * _userdata);

void zf_frame_init_S0(unsigned char * _p,
                      unsigned int    _M,
                      msequence ms,
                      gr_complex * _S0,
                      gr_complex * _s0,
                      unsigned int *  _M_S0);

void zf_frame_init_S1(unsigned char * _p,
                      unsigned int    _M,
                      msequence ms,
                      gr_complex * _S1,
                      gr_complex * _s1,
                      unsigned int *  _M_S1);

struct zf_framegen_s {
    unsigned int M;         // number of subcarriers
    unsigned int cp_len;    // cyclic prefix length
    unsigned char * p;      // subcarrier allocation (null, pilot, data)

    // tapering/trasition
    unsigned int taper_len; // number of samples in tapering window/overlap
    float * taper;          // tapering window
    gr_complex *postfix; // overlapping symbol buffer

    // constants
    unsigned int M_null;    // number of null subcarriers
    unsigned int M_pilot;   // number of pilot subcarriers
    unsigned int M_data;    // number of data subcarriers
    unsigned int M_S0;      // number of enabled subcarriers in S0
    unsigned int M_S1;      // number of enabled subcarriers in S1

    // scaling factors
    float g_data;           //

    // transform object
    fftwf_plan ifft;          // ifft object
    gr_complex * X;      // frequency-domain buffer
    gr_complex * x;      // time-domain buffer

    // PLCP short
    gr_complex * S0;     // short sequence (frequency)
    gr_complex * s0;     // short sequence (time)

    // PLCP long
    gr_complex * S1;     // long sequence (frequency)
    gr_complex * s1;     // long sequence (time)

    // pilot sequence
    msequence ms_pilot;
};

typedef struct zf_framegen_s * zf_framegen;

// receiver state
typedef enum 
{
    OFDMFRAMESYNC_STATE_SEEKPLCP=0,   // seek initial PLCP
    OFDMFRAMESYNC_STATE_PLCPSHORT0,   // seek first PLCP short sequence
    OFDMFRAMESYNC_STATE_PLCPSHORT1,   // seek second PLCP short sequence
    OFDMFRAMESYNC_STATE_PLCPLONG,     // seek PLCP long sequence
    OFDMFRAMESYNC_STATE_RXSYMBOLS     // receive payload symbols
} zf_framesync_state_t;

struct zf_framesync_s {
    unsigned int M;         // number of subcarriers
    unsigned int M2;        // number of subcarriers (divided by 2)
    unsigned int cp_len;    // cyclic prefix length
    unsigned char * p;      // subcarrier allocation (null, pilot, data)

    // constants
    unsigned int M_null;    // number of null subcarriers
    unsigned int M_pilot;   // number of pilot subcarriers
    unsigned int M_data;    // number of data subcarriers
    unsigned int M_S0;      // number of enabled subcarriers in S0
    unsigned int M_S1;      // number of enabled subcarriers in S1

    // scaling factors
    float g_data;           // data symbols gain
    float g_S0;             // S0 training symbols gain
    float g_S1;             // S1 training symbols gain

    // transform object
    fftwf_plan fft;           // ifft object
    gr_complex * X;      // frequency-domain buffer
    gr_complex * x;      // time-domain buffer
    windowcf input_buffer;  // input sequence buffer

    // PLCP sequences
    gr_complex * S0;     // short sequence (freq)
    gr_complex * s0;     // short sequence (time)
    gr_complex * S1;     // long sequence (freq)
    gr_complex * s1;     // long sequence (time)

    // gain
    float g0;               // nominal gain (coarse initial estimate)
    gr_complex * G0;     // complex subcarrier gain estimate, S0[0]
    gr_complex * G1;     // complex subcarrier gain estimate, S0[1]
    gr_complex * G;      // complex subcarrier gain estimate
    gr_complex * B;      // subcarrier phase rotation due to backoff
    gr_complex * R;      //
    
    // receiver state
    zf_framesync_state_t state;

    // synchronizer objects
    nco_crcf nco_rx;        // numerically-controlled oscillator
    msequence ms_pilot;     // pilot sequence generator
    float phi_prime;        // ...
    float p1_prime;         // filtered pilot phase slope

    // timing
    unsigned int timer;         // input sample timer
    unsigned int num_symbols;   // symbol counter
    unsigned int backoff;       // sample timing backoff
    gr_complex s_hat_0;      // first S0 symbol metrics estimate
    gr_complex s_hat_1;      // second S0 symbol metrics estimate

    // detection thresholds
    float plcp_detect_thresh;   // plcp detection threshold, nominally 0.35
    float plcp_sync_thresh;     // long symbol threshold, nominally 0.30

    // callback
    ofdmframesync_callback callback;
    void * userdata;

    unsigned long int sync_index;
    unsigned long long int num_samples_processed;

#if DEBUG_OFDMFRAMESYNC
    int debug_enabled;
    int debug_objects_created;
    windowcf debug_x;
    windowf  debug_rssi;
    windowcf debug_framesyms;
    gr_complex * G_hat;  // complex subcarrier gain estimate, S1
    float * px;             // pilot x-value
    float * py;             // pilot y-value
    float p_phase[2];       // pilot polyfit
    windowf debug_pilot_0;  // pilot polyfit history, p[0]
    windowf debug_pilot_1;  // pilot polyfit history, p[1]
#endif
};

typedef struct zf_framesync_s * zf_framesync;

unsigned long int
zf_framesync_get_sync_index(zf_framesync _q);

unsigned long long int
zf_framesync_get_num_samples_processed(zf_framesync _q);

// generate short sequence symbols
//  _p      :   subcarrier allocation array
//  _M      :   total number of subcarriers
//  _S0     :   output symbol (freq)
//  _s0     :   output symbol (time)
//  _M_S0   :   total number of enabled subcarriers in S0
void zf_frame_init_S0(unsigned char * _p,
                       unsigned int    _M,
                       gr_complex * _S0,
                       gr_complex * _s0,
                       unsigned int *  _M_S0);

// generate long sequence symbols
//  _p      :   subcarrier allocation array
//  _M      :   total number of subcarriers
//  _S1     :   output symbol (freq)
//  _s1     :   output symbol (time)
//  _M_S1   :   total number of enabled subcarriers in S1
void zf_frame_init_S1(unsigned char * _p,
                       unsigned int    _M,
                       gr_complex * _S1,
                       gr_complex * _s1,
                       unsigned int *  _M_S1);

// initialize default subcarrier allocation
//  _M      :   number of subcarriers
//  _p      :   output subcarrier allocation array, [size: _M x 1]
//
// key: '.' (null), 'P' (pilot), '+' (data)
// .+++P+++++++P.........P+++++++P+++
//
void zf_frame_init_default_sctype(unsigned int _M,
                                   unsigned char * _p);

// validate subcarrier type (count number of null, pilot, and data
// subcarriers in the allocation)
//  _p          :   subcarrier allocation array, [size: _M x 1]
//  _M          :   number of subcarriers
//  _M_null     :   output number of null subcarriers
//  _M_pilot    :   output number of pilot subcarriers
//  _M_data     :   output number of data subcarriers
void zf_frame_validate_sctype(unsigned char * _p,
                               unsigned int _M,
                               unsigned int * _M_null,
                               unsigned int * _M_pilot,
                               unsigned int * _M_data);

// print subcarrier allocation to screen
//
// key: '.' (null), 'P' (pilot), '+' (data)
// .+++P+++++++P.........P+++++++P+++
//
void zf_frame_print_sctype(unsigned char * _p,
                            unsigned int    _M);
// create OFDM framing synchronizer object
//  _M          :   number of subcarriers, >10 typical
//  _cp_len     :   cyclic prefix length
//  _taper_len  :   taper length (OFDM symbol overlap)
//  _p          :   subcarrier allocation (null, pilot, data), [size: _M x 1]
//  _callback   :   user-defined callback function
//  _userdata   :   user-defined data pointer
zf_framesync zf_framesync_create(unsigned int           _M,
                                   unsigned int           _cp_len,
                                   unsigned int           _taper_len,
                                   unsigned char *        _p,
                                   zf_framesync_callback _callback,
                                   void *                 _userdata,
                                   msequence ms0,
                                   msequence ms1);

void zf_framesync_destroy(zf_framesync _q);

void zf_framesync_reset(zf_framesync _q);

zf_framesync_state_t
zf_framesync_execute(zf_framesync _q,
                          std::vector<gr_complex *> _x,
                          unsigned int _n);

void zf_framesync_get_G(zf_framesync _q,
                        gr_complex ** G);

// get receiver RSSI
float zf_framesync_get_rssi(zf_framesync _q);

// frame detection
void zf_framesync_execute_seekplcp(zf_framesync _q);

// frame detection
void zf_framesync_execute_S0a(zf_framesync _q);

// frame detection
void zf_framesync_execute_S0b(zf_framesync _q);

void zf_framesync_execute_S1(zf_framesync _q);
void zf_framesync_execute_rxsymbols(zf_framesync _q);

// compute S0 metrics
void zf_framesync_S0_metrics(zf_framesync _q,
                              gr_complex * _G,
                              gr_complex * _s_hat);

// estimate short sequence gain
//  _q      :   zf_framesync object
//  _x      :   input array (time), [size: M x 1]
//  _G      :   output gain (freq)
void zf_framesync_estimate_gain_S0(zf_framesync   _q,
                                    gr_complex * _x,
                                    gr_complex * _G);

// estimate long sequence gain
//  _q      :   zf_framesync object
//  _x      :   input array (time), [size: M x 1]
//  _G      :   output gain (freq)
void zf_framesync_estimate_gain_S1(zf_framesync _q,
                                    gr_complex * _x,
                                    gr_complex * _G);

// estimate complex equalizer gain from G0 and G1
//  _q      :   zf_framesync object
//  _ntaps  :   number of time-domain taps for smoothing
void zf_framesync_estimate_eqgain(zf_framesync _q,
                                   unsigned int _ntaps);

// estimate complex equalizer gain from G0 and G1 using polynomial fit
//  _q      :   zf_framesync object
//  _order  :   polynomial order
void zf_framesync_estimate_eqgain_poly(zf_framesync _q,
                                        unsigned int _order);

// recover symbol, correcting for gain, pilot phase, etc.
void zf_framesync_rxsymbol(zf_framesync _q);

// enable debugging
void zf_framesync_debug_enable(zf_framesync _q);
void zf_framesync_debug_disable(zf_framesync _q);
void zf_framesync_debug_print(zf_framesync _q,
                               const char * _filename);

// create OFDM framing generator object
//  _M          :   number of subcarriers, >10 typical
//  _cp_len     :   cyclic prefix length
//  _taper_len  :   taper length (OFDM symbol overlap)
//  _p          :   subcarrier allocation (null, pilot, data), [size: _M x 1]
zf_framegen zf_framegen_create(unsigned int    _M,
                                 unsigned int    _cp_len,
                                 unsigned int    _taper_len,
                                 unsigned char * _p,
                                 msequence ms0,
                                 msequence ms1);

void zf_framegen_destroy(zf_framegen _q);
void zf_framegen_print(zf_framegen _q);
void zf_framegen_reset(zf_framegen _q);

// write first PLCP short sequence 'symbol' to buffer
//
//  |<- 2*cp->|<-       M       ->|<-       M       ->|
//  |         |                   |                   |
//      +-----+-------------------+-------------------+
//     /      |                   |                   |
//    /  ..s0 |        s0         |        s0         |
//   /        |                   |                   |
//  +---------+-------------------+-------------------+-----> time
//  |                        |                        |
//  |<-        s0[a]       ->|<-        s0[b]       ->|
//  |        M + cp_len      |        M + cp_len      |
//
void zf_framegen_write_S0a(zf_framegen    _q,
                            gr_complex * _y);

void zf_framegen_write_S0b(zf_framegen _q,
                            gr_complex * _y);

void zf_framegen_write_S1(zf_framegen _q,
                           gr_complex * _y);

unsigned int
zf_framegen_write_sync_words(zf_framegen _q,
                             std::vector<gr_complex *> fg_buff);

// write OFDM symbol
//  _q      :   framing generator object
//  _x      :   input symbols, [size: _M x 1]
//  _y      :   output samples, [size: _M x 1]
void zf_framegen_writesymbol(zf_framegen    _q,
                              gr_complex * _x,
                              gr_complex * _y);

// write tail to output
void zf_framegen_writetail(zf_framegen    _q,
                            gr_complex * _buffer);

// generate symbol (add cyclic prefix/postfix, overlap)
//
//  ->|   |<- taper_len
//    +   +-----+-------------------+
//     \ /      |                   |
//      X       |      symbol       |
//     / \      |                   |
//    +---+-----+-------------------+----> time
//    |         |                   |
//    |<- cp  ->|<-       M       ->|
//
//  _q->x           :   input time-domain symbol [size: _q->M x 1]
//  _q->postfix     :   input:  post-fix from previous symbol [size: _q->taper_len x 1]
//                      output: post-fix from this new symbol
//  _q->taper       :   tapering window
//  _q->taper_len   :   tapering window length
//
//  _buffer         :   output sample buffer [size: (_q->M + _q->cp_len) x 1]
void zf_framegen_gensymbol(zf_framegen    _q,
                            gr_complex * _buffer);

gr_complex liquid_cexpjf(float theta);
