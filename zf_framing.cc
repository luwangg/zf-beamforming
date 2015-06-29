#include "zf_framing.h"

#define cargf   std::arg
#define conjf   std::conj
#define cabsf   std::abs
#define crealf  std::real
#define cimagf  std::imag

gr_complex liquid_cexpjf(float theta)
{
  return std::polar(1.0f, theta);
}

void zf_frame_init_S0(unsigned char * _p,
                      unsigned int    _M,
                      msequence ms,
                      gr_complex * _S0,
                      gr_complex * _s0,
                      unsigned int *  _M_S0)
{
    unsigned int i;

    unsigned int s;
    unsigned int M_S0 = 0;

    // short sequence
    for (i=0; i<_M; i++) {
        // generate symbol
        //s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            _S0[i] = 0.0f;
        } else {
            if ( (i%2) == 0 ) {
                // even subcarrer
                _S0[i] = s ? 1.0f : -1.0f;
                M_S0++;
            } else {
                // odd subcarrer (ignore)
                _S0[i] = 0.0f;
            }
        }
    }

    // ensure at least one subcarrier was enabled
    if (M_S0 == 0) {
        fprintf(stderr,"error: zf_frame_init_S0(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // set return value(s)
    *_M_S0 = M_S0;

    // run inverse fft to get time-domain sequence
    fft_run(_M, _S0, _s0, LIQUID_FFT_BACKWARD, 0);

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S0);
    for (i=0; i<_M; i++)
        _s0[i] *= g;
}

void zf_frame_init_S1(unsigned char * _p,
                      unsigned int    _M,
                      msequence ms,
                      gr_complex * _S1,
                      gr_complex * _s1,
                      unsigned int *  _M_S1)
{
    unsigned int i;

    unsigned int s;
    unsigned int M_S1 = 0;

    // long sequence
    for (i=0; i<_M; i++) {
        // generate symbol
        //s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            _S1[i] = 0.0f;
        } else {
            _S1[i] = s ? 1.0f : -1.0f;
            M_S1++;
        }
    }

    // ensure at least one subcarrier was enabled
    if (M_S1 == 0) {
        fprintf(stderr,"error: zf_frame_init_S1(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // set return value(s)
    *_M_S1 = M_S1;

    // run inverse fft to get time-domain sequence
    fft_run(_M, _S1, _s1, LIQUID_FFT_BACKWARD, 0);

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S1);
    for (i=0; i<_M; i++)
        _s1[i] *= g;
}

void zf_frame_init_S0(unsigned char * _p,
                       unsigned int    _M,
                       gr_complex * _S0,
                       gr_complex * _s0,
                       unsigned int *  _M_S0)
{
    unsigned int i;

    // compute m-sequence length
    unsigned int m = liquid_nextpow2(_M);
    if (m < 4)      m = 4;
    else if (m > 8) m = 8;

    // generate m-sequence generator object
    msequence ms = msequence_create_default(m);

    unsigned int s;
    unsigned int M_S0 = 0;

    // short sequence
    for (i=0; i<_M; i++) {
        // generate symbol
        //s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            _S0[i] = 0.0f;
        } else {
            if ( (i%2) == 0 ) {
                // even subcarrer
                _S0[i] = s ? 1.0f : -1.0f;
                M_S0++;
            } else {
                // odd subcarrer (ignore)
                _S0[i] = 0.0f;
            }
        }
    }

    // destroy objects
    msequence_destroy(ms);

    // ensure at least one subcarrier was enabled
    if (M_S0 == 0) {
        fprintf(stderr,"error: zf_frame_init_S0(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // set return value(s)
    *_M_S0 = M_S0;

    // run inverse fft to get time-domain sequence
    fft_run(_M, _S0, _s0, LIQUID_FFT_BACKWARD, 0);

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S0);
    for (i=0; i<_M; i++)
        _s0[i] *= g;
}

void zf_frame_init_S1(unsigned char * _p,
                       unsigned int    _M,
                       gr_complex * _S1,
                       gr_complex * _s1,
                       unsigned int *  _M_S1)
{
    unsigned int i;

    // compute m-sequence length
    unsigned int m = liquid_nextpow2(_M);
    if (m < 4)      m = 4;
    else if (m > 8) m = 8;

    // increase m such that the resulting S1 sequence will
    // differ significantly from S0 with the same subcarrier
    // allocation array
    m++;

    // generate m-sequence generator object
    msequence ms = msequence_create_default(m);

    unsigned int s;
    unsigned int M_S1 = 0;

    // long sequence
    for (i=0; i<_M; i++) {
        // generate symbol
        //s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            _S1[i] = 0.0f;
        } else {
            _S1[i] = s ? 1.0f : -1.0f;
            M_S1++;
        }
    }

    // destroy objects
    msequence_destroy(ms);

    // ensure at least one subcarrier was enabled
    if (M_S1 == 0) {
        fprintf(stderr,"error: zf_frame_init_S1(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // set return value(s)
    *_M_S1 = M_S1;

    // run inverse fft to get time-domain sequence
    fft_run(_M, _S1, _s1, LIQUID_FFT_BACKWARD, 0);

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S1);
    for (i=0; i<_M; i++)
        _s1[i] *= g;
}

void zf_frame_init_default_sctype(unsigned int _M,
                                   unsigned char * _p)
{
    // validate input
    if (_M < 6) {
        fprintf(stderr,"warning: zf_frame_init_default_sctype(), less than 4 subcarriers\n");
    }

    unsigned int i;
    unsigned int M2 = _M/2;

    // compute guard band
    unsigned int G = _M / 10;
    if (G < 2) G = 2;

    // designate pilot spacing
    unsigned int P = (_M > 34) ? 8 : 4;
    unsigned int P2 = P/2;

    // initialize as NULL
    for (i=0; i<_M; i++)
        _p[i] = OFDMFRAME_SCTYPE_NULL;

    // upper band
    for (i=1; i<M2-G; i++) {
        if ( ((i+P2)%P) == 0 )
            _p[i] = OFDMFRAME_SCTYPE_PILOT;
        else
            _p[i] = OFDMFRAME_SCTYPE_DATA;
    }

    // lower band
    for (i=1; i<M2-G; i++) {
        unsigned int k = _M - i;
        if ( ((i+P2)%P) == 0 )
            _p[k] = OFDMFRAME_SCTYPE_PILOT;
        else
            _p[k] = OFDMFRAME_SCTYPE_DATA;
    }
}

void zf_frame_validate_sctype(unsigned char * _p,
                               unsigned int _M,
                               unsigned int * _M_null,
                               unsigned int * _M_pilot,
                               unsigned int * _M_data)
{
    // clear counters
    unsigned int M_null  = 0;
    unsigned int M_pilot = 0;
    unsigned int M_data  = 0;

    unsigned int i;
    for (i=0; i<_M; i++) {
        // update appropriate counters
        if (_p[i] == OFDMFRAME_SCTYPE_NULL)
            M_null++;
        else if (_p[i] == OFDMFRAME_SCTYPE_PILOT)
            M_pilot++;
        else if (_p[i] == OFDMFRAME_SCTYPE_DATA)
            M_data++;
        else {
            fprintf(stderr,"error: zf_frame_validate_sctype(), invalid subcarrier type (%u)\n", _p[i]);
            exit(1);
        }
    }

    // set outputs
    *_M_null  = M_null;
    *_M_pilot = M_pilot;
    *_M_data  = M_data;
}

void zf_frame_print_sctype(unsigned char * _p,
                            unsigned int    _M)
{
    unsigned int i;

    printf("[");
    for (i=0; i<_M; i++) {
        unsigned int k = (i + _M/2) % _M;

        switch (_p[k]) {
        case OFDMFRAME_SCTYPE_NULL:     printf(".");    break;
        case OFDMFRAME_SCTYPE_PILOT:    printf("|");    break;
        case OFDMFRAME_SCTYPE_DATA:     printf("+");    break;
        default:
            fprintf(stderr,"error: zf_frame_print_default_sctype(), invalid subcarrier type\n");
            exit(1);
        }
    }

    printf("]\n");
}

zf_framesync zf_framesync_create(unsigned int           _M,
                                   unsigned int           _cp_len,
                                   unsigned int           _taper_len,
                                   unsigned char *        _p,
                                   zf_framesync_callback _callback,
                                   void *                 _userdata,
                                   msequence ms0,
                                   msequence ms1)
{
    zf_framesync q = (zf_framesync) malloc(sizeof(struct zf_framesync_s));

    // validate input
    if (_M < 8) {
        fprintf(stderr,"warning: zf_framesync_create(), less than 8 subcarriers\n");
    } else if (_M % 2) {
        fprintf(stderr,"error: zf_framesync_create(), number of subcarriers must be even\n");
        exit(1);
    } else if (_cp_len > _M) {
        fprintf(stderr,"error: zf_framesync_create(), cyclic prefix length cannot exceed number of subcarriers\n");
        exit(1);
    }
    q->M = _M;
    q->cp_len = _cp_len;

    // derived values
    q->M2 = _M/2;

    // subcarrier allocation
    q->p = (unsigned char*) malloc((q->M)*sizeof(unsigned char));
    if (_p == NULL) {
        zf_frame_init_default_sctype(q->M, q->p);
    } else {
        memmove(q->p, _p, q->M*sizeof(unsigned char));
    }

    // validate and count subcarrier allocation
    zf_frame_validate_sctype(q->p, q->M, &q->M_null, &q->M_pilot, &q->M_data);
    if ( (q->M_pilot + q->M_data) == 0) {
        fprintf(stderr,"error: zf_framesync_create(), must have at least one enabled subcarrier\n");
        exit(1);
    } else if (q->M_data == 0) {
        fprintf(stderr,"error: zf_framesync_create(), must have at least one data subcarriers\n");
        exit(1);
    } else if (q->M_pilot < 2) {
        fprintf(stderr,"error: zf_framesync_create(), must have at least two pilot subcarriers\n");
        exit(1);
    }

    // create transform object
    q->X = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->x = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->fft = fftwf_plan_dft_1d(q->M,
        reinterpret_cast<fftwf_complex *>(q->x),
        reinterpret_cast<fftwf_complex *>(q->X),
                               FFTW_FORWARD,
                               FFTW_ESTIMATE);
 
    // create input buffer the length of the transform
    q->input_buffer = windowcf_create(q->M + q->cp_len);

    // allocate memory for PLCP arrays
    q->S0 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->s0 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->S1 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->s1 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));

    // generate S0 and S1
    #if USE_DEFAULT_MS
      zf_frame_init_S0(q->p, q->M, q->S0, q->s0, &q->M_S0);
      zf_frame_init_S1(q->p, q->M, q->S1, q->s1, &q->M_S1);
    #else
      zf_frame_init_S0(q->p, q->M, ms0, q->S0, q->s0, &q->M_S0);
      zf_frame_init_S1(q->p, q->M, ms1, q->S1, q->s1, &q->M_S1);
    #endif

    // compute scaling factor
    q->g_data = sqrtf(q->M) / sqrtf(q->M_pilot + q->M_data);
    q->g_S0   = sqrtf(q->M) / sqrtf(q->M_S0);
    q->g_S1   = sqrtf(q->M) / sqrtf(q->M_S1);

    // gain
    q->g0 = 1.0f;
    q->G0 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->G1 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->G  = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->B  = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->R  = (gr_complex*) malloc((q->M)*sizeof(gr_complex));

#if 1
    memset(q->G0, 0x00, q->M*sizeof(gr_complex));
    memset(q->G1, 0x00, q->M*sizeof(gr_complex));
    memset(q->G , 0x00, q->M*sizeof(gr_complex));
    memset(q->B,  0x00, q->M*sizeof(gr_complex));
#endif

    // timing backoff
    q->backoff = q->cp_len < 2 ? q->cp_len : 2;
    float phi = (float)(q->backoff)*2.0f*M_PI/(float)(q->M);
    unsigned int i;
    for (i=0; i<q->M; i++)
        q->B[i] = liquid_cexpjf(i*phi);

    // set callback data
    q->callback = _callback;
    q->userdata = _userdata;

    // 
    // synchronizer objects
    //

    // numerically-controlled oscillator
    q->nco_rx = nco_crcf_create(LIQUID_NCO);

    // set pilot sequence
    q->ms_pilot = msequence_create_default(8);

#if OFDMFRAMESYNC_ENABLE_SQUELCH
    // coarse detection
    q->squelch_threshold = -25.0f;
    q->squelch_enabled = 0;
#endif

    // reset object
    zf_framesync_reset(q);

#if DEBUG_OFDMFRAMESYNC
    q->debug_enabled = 0;
    q->debug_objects_created = 0;

    q->debug_x =        NULL;
    q->debug_rssi =     NULL;
    q->debug_framesyms =NULL;
    
    q->G_hat = NULL;
    q->px    = NULL;
    q->py    = NULL;
    
    q->debug_pilot_0 = NULL;
    q->debug_pilot_1 = NULL;
#endif

    // return object
    
    q->sync_index = 0;
    q->num_samples_processed = 0;
    return q;
}

void zf_framesync_destroy(zf_framesync _q)
{
#if DEBUG_OFDMFRAMESYNC
    // destroy debugging objects
    if (_q->debug_x         != NULL) windowcf_destroy(_q->debug_x);
    if (_q->debug_rssi      != NULL) windowf_destroy(_q->debug_rssi);
    if (_q->debug_framesyms != NULL) windowcf_destroy(_q->debug_framesyms);
    if (_q->G_hat           != NULL) free(_q->G_hat);
    if (_q->px              != NULL) free(_q->px);
    if (_q->py              != NULL) free(_q->py);
    if (_q->debug_pilot_0   != NULL) windowf_destroy(_q->debug_pilot_0);
    if (_q->debug_pilot_1   != NULL) windowf_destroy(_q->debug_pilot_1);
#endif

    // free subcarrier type array memory
    free(_q->p);

    // free transform object
    windowcf_destroy(_q->input_buffer);
    free(_q->X);
    free(_q->x);
    fftwf_destroy_plan(_q->fft);

    // clean up PLCP arrays
    free(_q->S0);
    free(_q->s0);
    free(_q->S1);
    free(_q->s1);

    // free gain arrays
    free(_q->G0);
    free(_q->G1);
    free(_q->G);
    free(_q->B);
    free(_q->R);

    // destroy synchronizer objects
    nco_crcf_destroy(_q->nco_rx);           // numerically-controlled oscillator
    msequence_destroy(_q->ms_pilot);

    // free main object memory
    free(_q);
}

void zf_framesync_print(zf_framesync _q)
{
    printf("zf_framesync:\n");
    printf("    num subcarriers     :   %-u\n", _q->M);
    printf("    cyclic prefix len   :   %-u\n", _q->cp_len);
    //printf("    taper len           :   %-u\n", _q->taper_len);
}

void zf_framesync_reset(zf_framesync _q)
{
#if 0
    // reset gain parameters
    unsigned int i;
    for (i=0; i<_q->M; i++)
        _q->G[i] = 1.0f;
#endif

    // reset synchronizer objects
    nco_crcf_reset(_q->nco_rx);
    msequence_reset(_q->ms_pilot);

    // reset timers
    _q->timer = 0;
    _q->num_symbols = 0;
    _q->s_hat_0 = 0.0f;
    _q->s_hat_1 = 0.0f;
    _q->phi_prime = 0.0f;
    _q->p1_prime = 0.0f;

    // set thresholds (increase for small number of subcarriers)
    _q->plcp_detect_thresh = (_q->M > 44) ? 0.35f : 0.35f + 0.01f*(44 - _q->M);
    _q->plcp_sync_thresh   = (_q->M > 44) ? 0.30f : 0.30f + 0.01f*(44 - _q->M);

    // reset state
    _q->state = OFDMFRAMESYNC_STATE_SEEKPLCP;
}

zf_framesync_state_t 
zf_framesync_execute(zf_framesync _q,
                     std::vector<gr_complex *> _x,
                     unsigned int _n)
{
    unsigned int i;
    gr_complex x;
    for (i=0; i<_n; i++) {
        x = _x[0][i];

        // correct for carrier frequency offset
        if (_q->state != OFDMFRAMESYNC_STATE_SEEKPLCP) {
            nco_crcf_mix_down(_q->nco_rx, x, &x);
            nco_crcf_step(_q->nco_rx);
        }

        // save input sample to buffer
        windowcf_push(_q->input_buffer,x);

#if DEBUG_OFDMFRAMESYNC
        if (_q->debug_enabled) {
            windowcf_push(_q->debug_x, x);
            windowf_push(_q->debug_rssi, crealf(x)*crealf(x) + cimagf(x)*cimagf(x));
        }
#endif

        switch (_q->state) {
        case OFDMFRAMESYNC_STATE_SEEKPLCP:
            zf_framesync_execute_seekplcp(_q);
            break;
        case OFDMFRAMESYNC_STATE_PLCPSHORT0:
            zf_framesync_execute_S0a(_q);
            break;
        case OFDMFRAMESYNC_STATE_PLCPSHORT1:
            zf_framesync_execute_S0b(_q);
            break;
        case OFDMFRAMESYNC_STATE_PLCPLONG:
            zf_framesync_execute_S1(_q);
            break;
        case OFDMFRAMESYNC_STATE_RXSYMBOLS:
            zf_framesync_execute_rxsymbols(_q);
            break;
        default:;
        }

    } // for (i=0; i<_n; i++)
    return _q->state;
} // zf_framesync_execute()

void zf_framesync_get_G(zf_framesync _q,
                        gr_complex ** G)
{
  ;
}

unsigned long int
zf_framesync_get_sync_index(zf_framesync _q)
{
  return _q->sync_index;
}

unsigned long long int
zf_framesync_get_num_samples_processed(zf_framesync _q)
{
  return _q->num_samples_processed;
}

float zf_framesync_get_rssi(zf_framesync _q)
{
    return -10.0f*log10(_q->g0);
}

float zf_framesync_get_cfo(zf_framesync _q)
{
    return nco_crcf_get_frequency(_q->nco_rx);
}

void zf_framesync_execute_seekplcp(zf_framesync _q)
{
    _q->timer++;

    if (_q->timer < _q->M)
        return;

    // reset timer
    _q->timer = 0;

    //
    gr_complex * rc;
    windowcf_read(_q->input_buffer, &rc);

    // estimate gain
    unsigned int i;
    float g = 0.0f;
    for (i=_q->cp_len; i<_q->M + _q->cp_len; i++) {
        // compute |rc[i]|^2 efficiently
        g += crealf(rc[i])*crealf(rc[i]) + cimagf(rc[i])*cimagf(rc[i]);
    }
    g = (float)(_q->M) / g;

#if OFDMFRAMESYNC_ENABLE_SQUELCH
    // TODO : squelch here
    if ( -10*log10f( sqrtf(g) ) < _q->squelch_threshold &&
         _q->squelch_enabled)
    {
        printf("squelch\n");
        return;
    }
#endif

    // estimate S0 gain
    zf_framesync_estimate_gain_S0(_q, &rc[_q->cp_len], _q->G0);

    gr_complex s_hat;
    zf_framesync_S0_metrics(_q, _q->G0, &s_hat);
    s_hat *= g;

    float tau_hat  = cargf(s_hat) * (float)(_q->M2) / (2*M_PI);
#if DEBUG_OFDMFRAMESYNC_PRINT
    printf(" - gain=%12.3f, rssi=%12.8f, s_hat=%12.4f <%12.8f>, tau_hat=%8.3f\n",
            sqrt(g),
            -10*log10(g),
            cabsf(s_hat), cargf(s_hat),
            tau_hat);
#endif

    // save gain (permits dynamic invocation of get_rssi() method)
    _q->g0 = g;

    // 
    if (cabsf(s_hat) > _q->plcp_detect_thresh) {

        int dt = (int)roundf(tau_hat);
        // set timer appropriately...
        _q->timer = (_q->M + dt) % (_q->M2);
        _q->timer += _q->M; // add delay to help ensure good S0 estimate
        _q->state = OFDMFRAMESYNC_STATE_PLCPSHORT0;

//#if DEBUG_OFDMFRAMESYNC_PRINT
        printf("********** frame detected! ************\n");
        printf("    s_hat   :   %12.8f <%12.8f>\n", cabsf(s_hat), cargf(s_hat));
        printf("  tau_hat   :   %12.8f\n", tau_hat);
        printf("    dt      :   %12d\n", dt);
        printf("    timer   :   %12u\n", _q->timer);
//#endif
        //printf("exiting prematurely\n");
        //zf_framesync_destroy(_q);
        //exit(1);
    }

}

void zf_framesync_execute_S0a(zf_framesync _q)
{
    //printf("t : %u\n", _q->timer);
    _q->timer++;

    if (_q->timer < _q->M2)
        return;

    // reset timer
    _q->timer = 0;

    //
    gr_complex * rc;
    windowcf_read(_q->input_buffer, &rc);

    // TODO : re-estimate nominal gain

    // estimate S0 gain
    zf_framesync_estimate_gain_S0(_q, &rc[_q->cp_len], _q->G0);

    gr_complex s_hat;
    zf_framesync_S0_metrics(_q, _q->G0, &s_hat);
    s_hat *= _q->g0;

    _q->s_hat_0 = s_hat;

#if DEBUG_OFDMFRAMESYNC_PRINT
    float tau_hat  = cargf(s_hat) * (float)(_q->M2) / (2*M_PI);
    printf("********** S0[0] received ************\n");
    printf("    s_hat   :   %12.8f <%12.8f>\n", cabsf(s_hat), cargf(s_hat));
    printf("  tau_hat   :   %12.8f\n", tau_hat);
#endif

#if 0
    // TODO : also check for phase of s_hat (should be small)
    if (cabsf(s_hat) < 0.3f) {
        // false alarm
#if DEBUG_OFDMFRAMESYNC_PRINT
        printf("false alarm S0[0]\n");
#endif
        zf_framesync_reset(_q);
        return;
    }
#endif
    _q->state = OFDMFRAMESYNC_STATE_PLCPSHORT1;
}

void zf_framesync_execute_S0b(zf_framesync _q)
{
    //printf("t = %u\n", _q->timer);
    _q->timer++;

    if (_q->timer < _q->M2)
        return;

    // reset timer
    _q->timer = _q->M + _q->cp_len - _q->backoff;

    //
    gr_complex * rc;
    windowcf_read(_q->input_buffer, &rc);

    // estimate S0 gain
    zf_framesync_estimate_gain_S0(_q, &rc[_q->cp_len], _q->G1);

    gr_complex s_hat;
    zf_framesync_S0_metrics(_q, _q->G1, &s_hat);
    s_hat *= _q->g0;

    _q->s_hat_1 = s_hat;

#if DEBUG_OFDMFRAMESYNC_PRINT
    float tau_hat  = cargf(s_hat) * (float)(_q->M2) / (2*M_PI);
    printf("********** S0[1] received ************\n");
    printf("    s_hat   :   %12.8f <%12.8f>\n", cabsf(s_hat), cargf(s_hat));
    printf("  tau_hat   :   %12.8f\n", tau_hat);

    // new timing offset estimate
    tau_hat  = cargf(_q->s_hat_0 + _q->s_hat_1) * (float)(_q->M2) / (2*M_PI);
    printf("  tau_hat * :   %12.8f\n", tau_hat);

    printf("**********\n");
#endif

    // re-adjust timer accordingly
    float tau_prime = cargf(_q->s_hat_0 + _q->s_hat_1) * (float)(_q->M2) / (2*M_PI);
    _q->timer -= (int)roundf(tau_prime);

#if 0
    if (cabsf(s_hat) < 0.3f) {
#if DEBUG_OFDMFRAMESYNC_PRINT
        printf("false alarm S0[1]\n");
#endif
        // false alarm
        zf_framesync_reset(_q);
        return;
    }
#endif

    gr_complex g_hat = 0.0f;
    unsigned int i;
    for (i=0; i<_q->M; i++)
        g_hat += _q->G1[i] * conjf(_q->G0[i]);

#if 0
    // compute carrier frequency offset estimate using freq. domain method
    float nu_hat = 2.0f * cargf(g_hat) / (float)(_q->M);
#else
    // compute carrier frequency offset estimate using ML method
    gr_complex t0 = 0.0f;
    for (i=0; i<_q->M2; i++) {
        t0 += conjf(rc[i])       *       _q->s0[i] * 
                    rc[i+_q->M2] * conjf(_q->s0[i+_q->M2]);
    }
    float nu_hat = cargf(t0) / (float)(_q->M2);
#endif

#if DEBUG_OFDMFRAMESYNC_PRINT
    printf("   nu_hat   :   %12.8f\n", nu_hat);
#endif

    // set NCO frequency
    nco_crcf_set_frequency(_q->nco_rx, nu_hat);

    _q->state = OFDMFRAMESYNC_STATE_PLCPLONG;
}

void zf_framesync_execute_S1(zf_framesync _q)
{
    _q->timer--;

    if (_q->timer > 0)
        return;

    // increment number of symbols observed
    _q->num_symbols++;

    // run fft
    gr_complex * rc;
    windowcf_read(_q->input_buffer, &rc);

    // estimate S1 gain
    // TODO : add backoff in gain estimation
    zf_framesync_estimate_gain_S1(_q, &rc[_q->cp_len], _q->G);

    // compute detector output
    gr_complex g_hat = 0.0f;
    unsigned int i;
    for (i=0; i<_q->M; i++) {
        //g_hat += _q->G[(i+1+_q->M)%_q->M]*conjf(_q->G[(i+_q->M)%_q->M]);
        g_hat += _q->G[(i+1)%_q->M]*conjf(_q->G[i]);
    }
    g_hat /= _q->M_S1; // normalize output
    g_hat *= _q->g0;

    // rotate by complex phasor relative to timing backoff
    g_hat *= liquid_cexpjf((float)(_q->backoff)*2.0f*M_PI/(float)(_q->M));

#if DEBUG_OFDMFRAMESYNC_PRINT
    printf("    g_hat   :   %12.4f <%12.8f>\n", cabsf(g_hat), cargf(g_hat));
#endif

    // check conditions for g_hat:
    //  1. magnitude should be large (near unity) when aligned
    //  2. phase should be very near zero (time aligned)
    if (cabsf(g_hat) > _q->plcp_sync_thresh && fabsf(cargf(g_hat)) < 0.1f*M_PI ) {
        //printf("    acquisition\n");
        _q->state = OFDMFRAMESYNC_STATE_RXSYMBOLS;
        // reset timer
        _q->timer = _q->M + _q->cp_len + _q->backoff;
        _q->num_symbols = 0;

        // normalize gain by subcarriers, apply timing backoff correction
        float g = (float)(_q->M) / sqrtf(_q->M_pilot + _q->M_data);
        for (i=0; i<_q->M; i++) {
            _q->G[i] *= g;          // gain due to relative subcarrier allocation
            _q->G[i] *= _q->B[i];   // timing backoff correction
        }

#if 0
        // TODO : choose number of taps more appropriately
        //unsigned int ntaps = _q->M / 4;
        unsigned int ntaps = (_q->M < 8) ? 2 : 8;
        // FIXME : this is by far the most computationally complex part of synchronization
        zf_framesync_estimate_eqgain(_q, ntaps);
#else
        unsigned int poly_order = 4;
        if (poly_order >= _q->M_pilot + _q->M_data)
            poly_order = _q->M_pilot + _q->M_data - 1;
        zf_framesync_estimate_eqgain_poly(_q, poly_order);
#endif

#if 1
        // compute composite gain
        unsigned int i;
        for (i=0; i<_q->M; i++)
            _q->R[i] = _q->B[i] / _q->G[i];
#endif

        return;
#if 0
        printf("exiting prematurely\n");
        zf_framesync_destroy(_q);
        exit(1);
#endif
    }

    // check if we are stuck searching for the S1 symbol
    if (_q->num_symbols == 16) {
#if DEBUG_OFDMFRAMESYNC_PRINT
        printf("could not find S1 symbol. bailing...\n");
#endif
        zf_framesync_reset(_q);
    }

    // 'reset' timer (wait another half symbol)
    _q->timer = _q->M2;
}

void zf_framesync_execute_rxsymbols(zf_framesync _q)
{
    // wait for timeout
    _q->timer--;

    if (_q->timer == 0) {

        // run fft
        gr_complex * rc;
        windowcf_read(_q->input_buffer, &rc);
        memmove(_q->x, &rc[_q->cp_len-_q->backoff], (_q->M)*sizeof(gr_complex));
        fftwf_execute(_q->fft);

        // recover symbol in internal _q->X buffer
        zf_framesync_rxsymbol(_q);

#if DEBUG_OFDMFRAMESYNC
        if (_q->debug_enabled) {
            unsigned int i;
            for (i=0; i<_q->M; i++) {
                if (_q->p[i] == OFDMFRAME_SCTYPE_DATA)
                    windowcf_push(_q->debug_framesyms, _q->X[i]);
            }
        }
#endif
        // invoke callback
        if (_q->callback != NULL) {
            int retval = _q->callback(_q->X, _q->p, _q->M, _q->userdata);

            if (retval != 0)
                zf_framesync_reset(_q);
        }

        // reset timer
        _q->timer = _q->M + _q->cp_len;
    }

}

void zf_framesync_S0_metrics(zf_framesync _q,
                              gr_complex * _G,
                              gr_complex * _s_hat)
{
    // timing, carrier offset correction
    unsigned int i;
    gr_complex s_hat = 0.0f;

    // compute timing estimate, accumulate phase difference across
    // gains on subsequent pilot subcarriers (note that all the odd
    // subcarriers are NULL)
    for (i=0; i<_q->M; i+=2) {
        s_hat += _G[(i+2)%_q->M]*conjf(_G[i]);
    }
    s_hat /= _q->M_S0; // normalize output

    // set output values
    *_s_hat = s_hat;
}

void zf_framesync_estimate_gain_S0(zf_framesync   _q,
                                    gr_complex * _x,
                                    gr_complex * _G)
{
    // move input array into fft input buffer
    memmove(_q->x, _x, (_q->M)*sizeof(gr_complex));

    // compute fft, storing result into _q->X
    fftwf_execute(_q->fft);
    
    // compute gain, ignoring NULL subcarriers
    unsigned int i;
    float gain = sqrtf(_q->M_S0) / (float)(_q->M);

    for (i=0; i<_q->M; i++) {
        if (_q->p[i] != OFDMFRAME_SCTYPE_NULL && (i%2)==0) {
            // NOTE : if cabsf(_q->S0[i]) == 0 then we can multiply by conjugate
            //        rather than compute division
            //_G[i] = _q->X[i] / _q->S0[i];
            _G[i] = _q->X[i] * conjf(_q->S0[i]);
        } else {
            _G[i] = 0.0f;
        }

        // normalize gain
        _G[i] *= gain;
    }
}

void zf_framesync_estimate_gain_S1(zf_framesync _q,
                                    gr_complex * _x,
                                    gr_complex * _G)
{
    // move input array into fft input buffer
    memmove(_q->x, _x, (_q->M)*sizeof(gr_complex));

    // compute fft, storing result into _q->X
    fftwf_execute(_q->fft);
    
    // compute gain, ignoring NULL subcarriers
    unsigned int i;
    float gain = sqrtf(_q->M_S1) / (float)(_q->M);
    for (i=0; i<_q->M; i++) {
        if (_q->p[i] != OFDMFRAME_SCTYPE_NULL) {
            // NOTE : if cabsf(_q->S1[i]) == 0 then we can multiply by conjugate
            //        rather than compute division
            //_G[i] = _q->X[i] / _q->S1[i];
            _G[i] = _q->X[i] * conjf(_q->S1[i]);
        } else {
            _G[i] = 0.0f;
        }

        // normalize gain
        _G[i] *= gain;
    }   
}

void zf_framesync_estimate_eqgain(zf_framesync _q,
                                   unsigned int _ntaps)
{
#if DEBUG_OFDMFRAMESYNC
    if (_q->debug_enabled) {
        // copy pre-smoothed gain
        memmove(_q->G_hat, _q->G, _q->M*sizeof(gr_complex));
    }
#endif

    // validate input
    if (_ntaps == 0 || _ntaps > _q->M) {
        fprintf(stderr, "error: zf_framesync_estimate_eqgain(), ntaps must be in [1,M]\n");
        exit(1);
    }

    unsigned int i;

    // generate smoothing window (fft of temporal window)
    for (i=0; i<_q->M; i++)
        _q->x[i] = (i < _ntaps) ? 1.0f : 0.0f;
    fftwf_execute(_q->fft);

    memmove(_q->G0, _q->G, _q->M*sizeof(gr_complex));

    // smooth complex equalizer gains
    for (i=0; i<_q->M; i++) {
        // set gain to zero for null subcarriers
        if (_q->p[i] == OFDMFRAME_SCTYPE_NULL) {
            _q->G[i] = 0.0f;
            continue;
        }

        gr_complex w;
        gr_complex w0 = 0.0f;
        gr_complex G_hat = 0.0f;

        unsigned int j;
        for (j=0; j<_q->M; j++) {
            if (_q->p[j] == OFDMFRAME_SCTYPE_NULL) continue;

            // select window sample from array
            w = _q->X[(i + _q->M - j) % _q->M];

            // accumulate gain
            //G_hat += w * 0.5f * (_q->G0[j] + _q->G1[j]);
            G_hat += w * _q->G0[j];
            w0 += w;
        }

        // eliminate divide-by-zero issues
        if (cabsf(w0) < 1e-4f) {
            fprintf(stderr,"error: zf_framesync_estimate_eqgain(), weighting factor is zero\n");
            w0 = 1.0f;
        }
        _q->G[i] = G_hat / w0;
    }
}

void zf_framesync_estimate_eqgain_poly(zf_framesync _q,
                                        unsigned int _order)
{
#if DEBUG_OFDMFRAMESYNC
    if (_q->debug_enabled) {
        // copy pre-smoothed gain
        memmove(_q->G_hat, _q->G, _q->M*sizeof(gr_complex));
    }
#endif

    // polynomial interpolation
    unsigned int i;
    unsigned int N = _q->M_pilot + _q->M_data;
    if (_order > N-1) _order = N-1;
    if (_order > 10)  _order = 10;
    float x_freq[N];
    float y_abs[N];
    float y_arg[N];
    float p_abs[_order+1];
    float p_arg[_order+1];

    unsigned int n=0;
    unsigned int k;
    for (i=0; i<_q->M; i++) {

        // start at mid-point (effective fftshift)
        k = (i + _q->M2) % _q->M;

        if (_q->p[k] != OFDMFRAME_SCTYPE_NULL) {
            if (n == N) {
                fprintf(stderr, "error: zf_framesync_estimate_eqgain_poly(), pilot subcarrier mismatch\n");
                exit(1);
            }
            // store resulting...
            x_freq[n] = (k > _q->M2) ? (float)k - (float)(_q->M) : (float)k;
            x_freq[n] = x_freq[n] / (float)(_q->M);
            y_abs[n] = cabsf(_q->G[k]);
            y_arg[n] = cargf(_q->G[k]);

            // update counter
            n++;
        }
    }

    if (n != N) {
        fprintf(stderr, "error: zf_framesync_estimate_eqgain_poly(), pilot subcarrier mismatch\n");
        exit(1);
    }

    // try to unwrap phase
    for (i=1; i<N; i++) {
        while ((y_arg[i] - y_arg[i-1]) >  M_PI)
            y_arg[i] -= 2*M_PI;
        while ((y_arg[i] - y_arg[i-1]) < -M_PI)
            y_arg[i] += 2*M_PI;
    }

    // fit to polynomial
    polyf_fit(x_freq, y_abs, N, p_abs, _order+1);
    polyf_fit(x_freq, y_arg, N, p_arg, _order+1);

    // compute subcarrier gain
    for (i=0; i<_q->M; i++) {
        float freq = (i > _q->M2) ? (float)i - (float)(_q->M) : (float)i;
        freq = freq / (float)(_q->M);
        float A     = polyf_val(p_abs, _order+1, freq);
        float theta = polyf_val(p_arg, _order+1, freq);
        _q->G[i] = (_q->p[i] == OFDMFRAME_SCTYPE_NULL) ? 0.0f : A * liquid_cexpjf(theta);
    }

#if 0
    for (i=0; i<N; i++)
        printf("x(%3u) = %12.8f; y_abs(%3u) = %12.8f; y_arg(%3u) = %12.8f;\n",
                i+1, x_freq[i],
                i+1, y_abs[i],
                i+1, y_arg[i]);

    for (i=0; i<=_order; i++)
        printf("p_abs(%3u) = %12.8f;\n", i+1, p_abs[i]);
    for (i=0; i<=_order; i++)
        printf("p_arg(%3u) = %12.8f;\n", i+1, p_arg[i]);
#endif
}

void zf_framesync_rxsymbol(zf_framesync _q)
{
    // apply gain
    unsigned int i;
    for (i=0; i<_q->M; i++)
        _q->X[i] *= _q->R[i];

    // polynomial curve-fit
    float x_phase[_q->M_pilot];
    float y_phase[_q->M_pilot];
    float p_phase[2];

    unsigned int n=0;
    unsigned int k;
    gr_complex pilot = 1.0f;
    for (i=0; i<_q->M; i++) {

        // start at mid-point (effective fftshift)
        k = (i + _q->M2) % _q->M;

        if (_q->p[k]==OFDMFRAME_SCTYPE_PILOT) {
            if (n == _q->M_pilot) {
                fprintf(stderr,"warning: zf_framesync_rxsymbol(), pilot subcarrier mismatch\n");
                return;
            }
            pilot = (msequence_advance(_q->ms_pilot) ? 1.0f : -1.0f);
#if 0
            printf("pilot[%3u] = %12.4e + j*%12.4e (expected %12.4e + j*%12.4e)\n",
                    k,
                    crealf(_q->X[k]), cimagf(_q->X[k]),
                    crealf(pilot),    cimagf(pilot));
#endif
            // store resulting...
            x_phase[n] = (k > _q->M2) ? (float)k - (float)(_q->M) : (float)k;
            y_phase[n] = cargf(_q->X[k]*conjf(pilot));

            // update counter
            n++;

        }
    }

    if (n != _q->M_pilot) {
        fprintf(stderr,"warning: zf_framesync_rxsymbol(), pilot subcarrier mismatch\n");
        return;
    }

    // try to unwrap phase
    for (i=1; i<_q->M_pilot; i++) {
        while ((y_phase[i] - y_phase[i-1]) >  M_PI)
            y_phase[i] -= 2*M_PI;
        while ((y_phase[i] - y_phase[i-1]) < -M_PI)
            y_phase[i] += 2*M_PI;
    }

    // fit phase to 1st-order polynomial (2 coefficients)
    polyf_fit(x_phase, y_phase, _q->M_pilot, p_phase, 2);

    // filter slope estimate (timing offset)
    float alpha = 0.3f;
    p_phase[1] = alpha*p_phase[1] + (1-alpha)*_q->p1_prime;
    _q->p1_prime = p_phase[1];

#if DEBUG_OFDMFRAMESYNC
    if (_q->debug_enabled) {
        // save pilots
        memmove(_q->px, x_phase, _q->M_pilot*sizeof(float));
        memmove(_q->py, y_phase, _q->M_pilot*sizeof(float));

        // NOTE : swapping values for octave
        _q->p_phase[0] = p_phase[1];
        _q->p_phase[1] = p_phase[0];

        windowf_push(_q->debug_pilot_0, p_phase[0]);
        windowf_push(_q->debug_pilot_1, p_phase[1]);
    }
#endif

    // compensate for phase offset
    // TODO : find more computationally efficient way to do this
    for (i=0; i<_q->M; i++) {
        // only apply to data/pilot subcarriers
        if (_q->p[i] == OFDMFRAME_SCTYPE_NULL) {
            _q->X[i] = 0.0f;
        } else {
            float fx    = (i > _q->M2) ? (float)i - (float)(_q->M) : (float)i;
            float theta = polyf_val(p_phase, 2, fx);
            _q->X[i] *= liquid_cexpjf(-theta);
        }
    }

    // adjust NCO frequency based on differential phase
    if (_q->num_symbols > 0) {
        // compute phase error (unwrapped)
        float dphi_prime = p_phase[0] - _q->phi_prime;
        while (dphi_prime >  M_PI) dphi_prime -= M_2_PI;
        while (dphi_prime < -M_PI) dphi_prime += M_2_PI;

        // adjust NCO proportionally to phase error
        nco_crcf_adjust_frequency(_q->nco_rx, 1e-3f*dphi_prime);
    }
    // set internal phase state
    _q->phi_prime = p_phase[0];
    //printf("%3u : theta : %12.8f, nco freq: %12.8f\n", _q->num_symbols, p_phase[0], nco_crcf_get_frequency(_q->nco_rx));
    
    // increment symbol counter
    _q->num_symbols++;

#if 0
    for (i=0; i<_q->M_pilot; i++)
        printf("x_phase(%3u) = %12.8f; y_phase(%3u) = %12.8f;\n", i+1, x_phase[i], i+1, y_phase[i]);
    printf("poly : p0=%12.8f, p1=%12.8f\n", p_phase[0], p_phase[1]);
#endif
}

void zf_framesync_debug_enable(zf_framesync _q)
{
    // create debugging objects if necessary
#if DEBUG_OFDMFRAMESYNC
    if (_q->debug_objects_created)
        return;

    _q->debug_x         = windowcf_create(DEBUG_OFDMFRAMESYNC_BUFFER_LEN);
    _q->debug_rssi      = windowf_create(DEBUG_OFDMFRAMESYNC_BUFFER_LEN);
    _q->debug_framesyms = windowcf_create(DEBUG_OFDMFRAMESYNC_BUFFER_LEN);
    _q->G_hat           = (gr_complex*) malloc((_q->M)*sizeof(gr_complex));

    _q->px = (float*) malloc((_q->M_pilot)*sizeof(float));
    _q->py = (float*) malloc((_q->M_pilot)*sizeof(float));

    _q->debug_pilot_0 = windowf_create(DEBUG_OFDMFRAMESYNC_BUFFER_LEN);
    _q->debug_pilot_1 = windowf_create(DEBUG_OFDMFRAMESYNC_BUFFER_LEN);

    _q->debug_enabled   = 1;
    _q->debug_objects_created = 1;
#else
    fprintf(stderr,"zf_framesync_debug_enable(): compile-time debugging disabled\n");
#endif
}

void zf_framesync_debug_disable(zf_framesync _q)
{
    // disable debugging
#if DEBUG_OFDMFRAMESYNC
    _q->debug_enabled = 0;
#else
    fprintf(stderr,"zf_framesync_debug_enable(): compile-time debugging disabled\n");
#endif
}

void zf_framesync_debug_print(zf_framesync _q,
                               const char * _filename)
{
#if DEBUG_OFDMFRAMESYNC
    if (!_q->debug_objects_created) {
        fprintf(stderr,"error: zf_frame_debug_print(), debugging objects don't exist; enable debugging first\n");
        return;
    }

    FILE * fid = fopen(_filename,"w");
    if (!fid) {
        fprintf(stderr,"error: zf_frame_debug_print(), could not open '%s' for writing\n", _filename);
        return;
    }
    fprintf(fid,"%% %s : auto-generated file\n", DEBUG_OFDMFRAMESYNC_FILENAME);
    fprintf(fid,"close all;\n");
    fprintf(fid,"clear all;\n");
    fprintf(fid,"n = %u;\n", DEBUG_OFDMFRAMESYNC_BUFFER_LEN);
    fprintf(fid,"M = %u;\n", _q->M);
    fprintf(fid,"M_null  = %u;\n", _q->M_null);
    fprintf(fid,"M_pilot = %u;\n", _q->M_pilot);
    fprintf(fid,"M_data  = %u;\n", _q->M_data);
    unsigned int i;
    gr_complex * rc;
    float * r;

    // save subcarrier allocation
    fprintf(fid,"p = zeros(1,M);\n");
    for (i=0; i<_q->M; i++)
        fprintf(fid,"p(%4u) = %d;\n", i+1, _q->p[i]);
    fprintf(fid,"i_null  = find(p==%d);\n", OFDMFRAME_SCTYPE_NULL);
    fprintf(fid,"i_pilot = find(p==%d);\n", OFDMFRAME_SCTYPE_PILOT);
    fprintf(fid,"i_data  = find(p==%d);\n", OFDMFRAME_SCTYPE_DATA);

    // short, long, training sequences
    for (i=0; i<_q->M; i++) {
        fprintf(fid,"S0(%4u) = %12.4e + j*%12.4e;\n", i+1, crealf(_q->S0[i]), cimagf(_q->S0[i]));
        fprintf(fid,"S1(%4u) = %12.4e + j*%12.4e;\n", i+1, crealf(_q->S1[i]), cimagf(_q->S1[i]));
    }

    fprintf(fid,"x = zeros(1,n);\n");
    windowcf_read(_q->debug_x, &rc);
    for (i=0; i<DEBUG_OFDMFRAMESYNC_BUFFER_LEN; i++)
        fprintf(fid,"x(%4u) = %12.4e + j*%12.4e;\n", i+1, crealf(rc[i]), cimagf(rc[i]));
    fprintf(fid,"figure;\n");
    fprintf(fid,"plot(0:(n-1),real(x),0:(n-1),imag(x));\n");
    fprintf(fid,"xlabel('sample index');\n");
    fprintf(fid,"ylabel('received signal, x');\n");


    fprintf(fid,"s1 = [];\n");
    for (i=0; i<_q->M; i++)
        fprintf(fid,"s1(%3u) = %12.4e + j*%12.4e;\n", i+1, crealf(_q->s1[i]), cimagf(_q->s1[i]));


    // write agc_rssi
    fprintf(fid,"\n\n");
    fprintf(fid,"agc_rssi = zeros(1,%u);\n", DEBUG_OFDMFRAMESYNC_BUFFER_LEN);
    windowf_read(_q->debug_rssi, &r);
    for (i=0; i<DEBUG_OFDMFRAMESYNC_BUFFER_LEN; i++)
        fprintf(fid,"agc_rssi(%4u) = %12.4e;\n", i+1, r[i]);
    fprintf(fid,"agc_rssi = filter([0.00362168 0.00724336 0.00362168],[1 -1.82269490 0.83718163],agc_rssi);\n");
    fprintf(fid,"agc_rssi = 10*log10( agc_rssi );\n");
    fprintf(fid,"figure;\n");
    fprintf(fid,"plot(agc_rssi)\n");
    fprintf(fid,"ylabel('RSSI [dB]');\n");

    // write short, long symbols
    fprintf(fid,"\n\n");
    fprintf(fid,"S0 = zeros(1,M);\n");
    fprintf(fid,"S1 = zeros(1,M);\n");
    for (i=0; i<_q->M; i++) {
        fprintf(fid,"S0(%3u) = %12.8f + j*%12.8f;\n", i+1, crealf(_q->S0[i]), cimagf(_q->S0[i]));
        fprintf(fid,"S1(%3u) = %12.8f + j*%12.8f;\n", i+1, crealf(_q->S1[i]), cimagf(_q->S1[i]));
    }


    // write gain arrays
    fprintf(fid,"\n\n");
    fprintf(fid,"G0     = zeros(1,M);\n");
    fprintf(fid,"G1     = zeros(1,M);\n");
    fprintf(fid,"G_hat  = zeros(1,M);\n");
    fprintf(fid,"G      = zeros(1,M);\n");
    for (i=0; i<_q->M; i++) {
        fprintf(fid,"G0(%3u)    = %12.8f + j*%12.8f;\n", i+1, crealf(_q->G0[i]),   cimagf(_q->G0[i]));
        fprintf(fid,"G1(%3u)    = %12.8f + j*%12.8f;\n", i+1, crealf(_q->G1[i]),   cimagf(_q->G1[i]));
        fprintf(fid,"G_hat(%3u) = %12.8f + j*%12.8f;\n", i+1, crealf(_q->G_hat[i]),cimagf(_q->G_hat[i]));
        fprintf(fid,"G(%3u)     = %12.8f + j*%12.8f;\n", i+1, crealf(_q->G[i]),    cimagf(_q->G[i]));
    }
    fprintf(fid,"f = [0:(M-1)];\n");
    fprintf(fid,"figure;\n");
    fprintf(fid,"subplot(2,1,1);\n");
    fprintf(fid,"  plot(f, fftshift(abs(G_hat)),'sb',...\n");
    fprintf(fid,"       f, fftshift(abs(G)),'-k','LineWidth',2);\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  xlabel('subcarrier index');\n");
    fprintf(fid,"  ylabel('gain estimate (mag)');\n");
    fprintf(fid,"subplot(2,1,2);\n");
    fprintf(fid,"  plot(f, fftshift(arg(G_hat).*[abs(G0) > 1e-3]),'sb',...\n");
    fprintf(fid,"       f, fftshift(arg(G)),'-k','LineWidth',2);\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  xlabel('subcarrier index');\n");
    fprintf(fid,"  ylabel('gain estimate (phase)');\n");

    // write pilot response
    fprintf(fid,"\n\n");
    fprintf(fid,"px = zeros(1,M_pilot);\n");
    fprintf(fid,"py = zeros(1,M_pilot);\n");
    for (i=0; i<_q->M_pilot; i++) {
        fprintf(fid,"px(%3u) = %12.8f;\n", i+1, _q->px[i]);
        fprintf(fid,"py(%3u) = %12.8f;\n", i+1, _q->py[i]);
    }
    fprintf(fid,"p_phase(1) = %12.8f;\n", _q->p_phase[0]);
    fprintf(fid,"p_phase(2) = %12.8f;\n", _q->p_phase[1]);

    // save pilot history
    fprintf(fid,"p0 = zeros(1,M);\n");
    windowf_read(_q->debug_pilot_0, &r);
    for (i=0; i<DEBUG_OFDMFRAMESYNC_BUFFER_LEN; i++)
        fprintf(fid,"p0(%4u) = %12.4e;\n", i+1, r[i]);

    fprintf(fid,"p1 = zeros(1,M);\n");
    windowf_read(_q->debug_pilot_1, &r);
    for (i=0; i<DEBUG_OFDMFRAMESYNC_BUFFER_LEN; i++)
        fprintf(fid,"p1(%4u) = %12.4e;\n", i+1, r[i]);

    fprintf(fid,"figure;\n");
    fprintf(fid,"fp = (-M/2):(M/2);\n");
    fprintf(fid,"subplot(3,1,1);\n");
    fprintf(fid,"  plot(px, py, 'sb',...\n");
    fprintf(fid,"       fp, polyval(p_phase, fp), '-k');\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  legend('pilots','polyfit',0);\n");
    fprintf(fid,"  xlabel('subcarrier');\n");
    fprintf(fid,"  ylabel('phase');\n");
    fprintf(fid,"subplot(3,1,2);\n");
    fprintf(fid,"  plot(1:length(p0), p0);\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  ylabel('p0 (phase offset)');\n");
    fprintf(fid,"subplot(3,1,3);\n");
    fprintf(fid,"  plot(1:length(p1), p1);\n");
    fprintf(fid,"  grid on;\n");
    fprintf(fid,"  ylabel('p1 (phase slope)');\n");

    // write frame symbols
    fprintf(fid,"framesyms = zeros(1,n);\n");
    windowcf_read(_q->debug_framesyms, &rc);
    for (i=0; i<DEBUG_OFDMFRAMESYNC_BUFFER_LEN; i++)
        fprintf(fid,"framesyms(%4u) = %12.4e + j*%12.4e;\n", i+1, crealf(rc[i]), cimagf(rc[i]));
    fprintf(fid,"figure;\n");
    fprintf(fid,"plot(real(framesyms), imag(framesyms), 'x');\n");
    fprintf(fid,"xlabel('I');\n");
    fprintf(fid,"ylabel('Q');\n");
    fprintf(fid,"axis([-1 1 -1 1]*1.6);\n");
    fprintf(fid,"axis square;\n");
    fprintf(fid,"grid on;\n");

    fclose(fid);
    printf("zf_framesync/debug: results written to '%s'\n", _filename);
#else
    fprintf(stderr,"zf_framesync_debug_print(): compile-time debugging disabled\n");
#endif
}

zf_framegen zf_framegen_create(unsigned int    _M,
                               unsigned int    _cp_len,
                               unsigned int    _taper_len,
                               unsigned char * _p,
                               msequence ms0,
                               msequence ms1)
{
    // validate input
    if (_M < 2) {
        fprintf(stderr,"error: zf_framegen_create(), number of subcarriers must be at least 2\n");
        exit(1);
    } else if (_M % 2) {
        fprintf(stderr,"error: zf_framegen_create(), number of subcarriers must be even\n");
        exit(1);
    } else if (_cp_len > _M) {
        fprintf(stderr,"error: zf_framegen_create(), cyclic prefix cannot exceed symbol length\n");
        exit(1);
    } else if (_taper_len > _cp_len) {
        fprintf(stderr,"error: zf_framegen_create(), taper length cannot exceed cyclic prefix\n");
        exit(1);
    }

    zf_framegen q = (zf_framegen) malloc(sizeof(struct zf_framegen_s));
    q->M         = _M;
    q->cp_len    = _cp_len;
    q->taper_len = _taper_len;

    // allocate memory for subcarrier allocation IDs
    q->p = (unsigned char*) malloc((q->M)*sizeof(unsigned char));
    if (_p == NULL) {
        // initialize default subcarrier allocation
        zf_frame_init_default_sctype(q->M, q->p);
    } else {
        // copy user-defined subcarrier allocation
        memmove(q->p, _p, q->M*sizeof(unsigned char));
    }

    // validate and count subcarrier allocation
    zf_frame_validate_sctype(q->p, q->M, &q->M_null, &q->M_pilot, &q->M_data);
    if ( (q->M_pilot + q->M_data) == 0) {
        fprintf(stderr,"error: zf_framegen_create(), must have at least one enabled subcarrier\n");
        exit(1);
    } else if (q->M_data == 0) {
        fprintf(stderr,"error: zf_framegen_create(), must have at least one data subcarriers\n");
        exit(1);
    } else if (q->M_pilot < 2) {
        fprintf(stderr,"error: zf_framegen_create(), must have at least two pilot subcarriers\n");
        exit(1);
    }

    unsigned int i;

    // allocate memory for transform objects
    q->X = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->x = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->ifft = fftwf_plan_dft_1d(q->M,
        reinterpret_cast<fftwf_complex *>(q->X),
        reinterpret_cast<fftwf_complex *>(q->x),
                                FFTW_BACKWARD,
                                FFTW_ESTIMATE);

    // allocate memory for PLCP arrays
    q->S0 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->s0 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->S1 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));
    q->s1 = (gr_complex*) malloc((q->M)*sizeof(gr_complex));

    // generate S0 and S1
    #if USE_DEFAULT_MS
      zf_frame_init_S0(q->p, q->M, q->S0, q->s0, &q->M_S0);
      zf_frame_init_S1(q->p, q->M, q->S1, q->s1, &q->M_S1);
    #else
      zf_frame_init_S0(q->p, q->M, ms0, q->S0, q->s0, &q->M_S0);
      zf_frame_init_S1(q->p, q->M, ms1, q->S1, q->s1, &q->M_S1);
    #endif

    // create tapering window and transition buffer
    q->taper   = (float*)         malloc(q->taper_len * sizeof(float));
    q->postfix = (gr_complex*) malloc(q->taper_len * sizeof(gr_complex));
    for (i=0; i<q->taper_len; i++) {
        float t = ((float)i + 0.5f) / (float)(q->taper_len);
        float g = sinf(M_PI_2*t);
        q->taper[i] = g*g;
    }
#if 0
    // validate window symmetry
    for (i=0; i<q->taper_len; i++) {
        printf("    taper[%2u] = %12.8f (%12.8f)\n", i, q->taper[i],
            q->taper[i] + q->taper[q->taper_len - i - 1]);
    }
#endif

    // compute scaling factor
    q->g_data = 1.0f / sqrtf(q->M_pilot + q->M_data);

    // set pilot sequence
    q->ms_pilot = msequence_create_default(8);

    return q;
}

void zf_framegen_destroy(zf_framegen _q)
{
    // free subcarrier type array memory
    free(_q->p);

    // free transform array memory
    free(_q->X);
    free(_q->x);
    fftwf_destroy_plan(_q->ifft);

    // free tapering window and transition buffer
    free(_q->taper);
    free(_q->postfix);

    // free PLCP memory arrays
    free(_q->S0);
    free(_q->s0);
    free(_q->S1);
    free(_q->s1);

    // free pilot msequence object memory
    msequence_destroy(_q->ms_pilot);

    // free main object memory
    free(_q);
}

void zf_framegen_print(zf_framegen _q)
{
    printf("zf_framegen:\n");
    printf("    num subcarriers     :   %-u\n", _q->M);
    printf("      - NULL            :   %-u\n", _q->M_null);
    printf("      - pilot           :   %-u\n", _q->M_pilot);
    printf("      - data            :   %-u\n", _q->M_data);
    printf("    cyclic prefix len   :   %-u\n", _q->cp_len);
    printf("    taper len           :   %-u\n", _q->taper_len);
    printf("    ");
    zf_frame_print_sctype(_q->p, _q->M);
}

void zf_framegen_reset(zf_framegen _q)
{
    msequence_reset(_q->ms_pilot);

    // clear internal postfix buffer
    unsigned int i;
    for (i=0; i<_q->taper_len; i++)
        _q->postfix[i] = 0.0f;
}

void zf_framegen_write_S0a(zf_framegen    _q,
                            gr_complex * _y)
{
    unsigned int i;
    unsigned int k;
    for (i=0; i<_q->M + _q->cp_len; i++) {
        k = (i + _q->M - 2*_q->cp_len) % _q->M;
        _y[i] = _q->s0[k];
    }

    // apply tapering window
    for (i=0; i<_q->taper_len; i++)
        _y[i] *= _q->taper[i];
}

void zf_framegen_write_S0b(zf_framegen _q,
                            gr_complex * _y)
{
    unsigned int i;
    unsigned int k;
    for (i=0; i<_q->M + _q->cp_len; i++) {
        k = (i + _q->M - _q->cp_len) % _q->M;
        _y[i] = _q->s0[k];
    }

    // copy postfix (first 'taper_len' samples of s0 symbol)
    memmove(_q->postfix, _q->s0, _q->taper_len*sizeof(gr_complex));
}

void zf_framegen_write_S1(zf_framegen _q,
                           gr_complex * _y)
{
    // copy S1 symbol to output, adding cyclic prefix and tapering window
    memmove(_q->x, _q->s1, (_q->M)*sizeof(gr_complex));
    zf_framegen_gensymbol(_q, _y);
}

void zf_framegen_writesymbol(zf_framegen    _q,
                              gr_complex * _x,
                              gr_complex * _y)
{
    // move frequency data to internal buffer
    unsigned int i;
    unsigned int k;
    int sctype;
    for (i=0; i<_q->M; i++) {
        // start at mid-point (effective fftshift)
        k = (i + _q->M/2) % _q->M;

        sctype = _q->p[k];
        if (sctype==OFDMFRAME_SCTYPE_NULL) {
            // disabled subcarrier
            _q->X[k] = 0.0f;
        } else if (sctype==OFDMFRAME_SCTYPE_PILOT) {
            // pilot subcarrier
            _q->X[k] = (msequence_advance(_q->ms_pilot) ? 1.0f : -1.0f) * _q->g_data;
        } else {
            // data subcarrier
            _q->X[k] = _x[k] * _q->g_data;
        }

        //printf("X[%3u] = %12.8f + j*%12.8f;\n",i+1,crealf(_q->X[i]),cimagf(_q->X[i]));
    }

    // execute transform
    fftwf_execute(_q->ifft);

    // copy result to output, adding cyclic prefix and tapering window
    zf_framegen_gensymbol(_q, _y);
}

unsigned int
zf_framegen_write_sync_words(zf_framegen _q,
                             std::vector<gr_complex *> fg_buff)
{
  unsigned int num_streams = fg_buff.size();
  unsigned int sample_index = 0;
  unsigned int sps = (_q->M) + (_q->cp_len);
  unsigned int stream;

  gr_complex z(0.0, 0.0);
  gr_complex * zBuff = (gr_complex *) malloc
                       (sizeof(gr_complex)*sps);
  std::fill(zBuff, zBuff + sps, z);

  // write onto fg_buff[0] and put zeros on other streams
  zf_framegen_write_S0a(_q,
                        fg_buff[0] + sample_index);
  for(stream = 1; stream < num_streams; stream++){
    memmove(fg_buff[stream] + sample_index,
            zBuff,
            sizeof(gr_complex)*sps);
  }
  sample_index += (_q->M) + (_q->cp_len);
  zf_framegen_write_S0b(_q,
                        fg_buff[0] + sample_index);
  for(stream = 1; stream < num_streams; stream++){
    memmove(fg_buff[stream] + sample_index,
            zBuff,
            sizeof(gr_complex)*sps);
  }
  sample_index += (_q->M) + (_q->cp_len);
  zf_framegen_write_S1(_q,
                       fg_buff[0] + sample_index);
  for(stream = 1; stream < num_streams; stream++){
    memmove(fg_buff[stream] + sample_index,
            zBuff,
            sizeof(gr_complex)*sps);
  }
  sample_index += (_q->M) + (_q->cp_len);
  for(stream = 1; stream < num_streams; stream++){
    // S0a
    for(unsigned int chan = 0; chan < num_streams; chan++){
      memmove(fg_buff[chan] + sample_index,
              zBuff,
              sizeof(gr_complex)*sps);
    }
    sample_index += (_q->M) + (_q->cp_len);
    // S0b
    for(unsigned int chan = 0; chan < num_streams; chan++){
      memmove(fg_buff[chan] + sample_index,
              zBuff,
              sizeof(gr_complex)*sps);
    }
    sample_index += (_q->M) + (_q->cp_len);
    // S1
    for(unsigned int chan = 0; chan < num_streams; chan++){
      memmove(fg_buff[chan] + sample_index,
              zBuff,
              sizeof(gr_complex)*sps);
    }
    sample_index += (_q->M) + (_q->cp_len);
  }
  return sample_index;
}

// write tail to output
void zf_framegen_writetail(zf_framegen    _q,
                            gr_complex * _buffer)
{
    // write tail to output, applying tapering window
    unsigned int i;
    for (i=0; i<_q->taper_len; i++)
        _buffer[i] = _q->postfix[i] * _q->taper[_q->taper_len-i-1];
}

void zf_framegen_gensymbol(zf_framegen    _q,
                            gr_complex * _buffer)
{
    // copy input symbol with cyclic prefix to output symbol
    memmove( &_buffer[0],          &_q->x[_q->M-_q->cp_len], _q->cp_len*sizeof(gr_complex));
    memmove( &_buffer[_q->cp_len], &_q->x[               0], _q->M    * sizeof(gr_complex));
    
    // apply tapering window to over-lapping regions
    unsigned int i;
    for (i=0; i<_q->taper_len; i++) {
        _buffer[i] *= _q->taper[i];
        _buffer[i] += _q->postfix[i] * _q->taper[_q->taper_len-i-1];
    }

    // copy post-fix to output (first 'taper_len' samples of input symbol)
    memmove(_q->postfix, _q->x, _q->taper_len*sizeof(gr_complex));
}

