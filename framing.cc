#include "framing.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <volk/volk.h>
#include <volk/volk_malloc.h>
#include <string.h>

#define DEBUG_PRINT true

framegen::framegen(unsigned int _M,
                   unsigned int _cp_len,
                   unsigned int _num_streams,
                   std::vector<unsigned char *> const &_p,
                   std::vector<msequence> const &_ms_S0,
                   std::vector<msequence> const &_ms_S1,
                   std::vector<msequence> const &_ms_pilot) 
{
  // sanity checks for parameters
  assert(_p.size() == _num_streams);
  assert(_ms_S0.size() == _num_streams);
  assert(_ms_S1.size() == _num_streams);
  assert(_ms_pilot.size() == _num_streams);

  // assign values for variables
  M = _M;
  cp_len = _cp_len;
  symbol_len = M + cp_len;
  num_streams = _num_streams;

  // resize vectors to num streams
  p.resize(num_streams);
  M_null.resize(num_streams);
  M_pilot.resize(num_streams);
  M_data.resize(num_streams);
  M_S0.resize(num_streams);
  M_S1.resize(num_streams);
  g_data.resize(num_streams);
  W = (std::complex<float> **) malloc 
      (sizeof(std::complex<float>*)*num_streams);
  ifft.resize(num_streams);
  X.resize(num_streams);
  x.resize(num_streams);
  symbol.resize(num_streams);
  S0.resize(num_streams);
  s0.resize(num_streams);
  S1.resize(num_streams);
  s1.resize(num_streams);
  ms_pilot.resize(num_streams);

  // get volk alignment for volk_malloc
  size_t volk_alignment = volk_get_alignment();
  std::complex<float> I(1.0, 0.0);
  for(unsigned int i = 0; i < num_streams; i++)
  {
    // memory allocation
    // NOTE: all to be freed in destructor
    p[i]  = (unsigned char *) malloc
            (sizeof(unsigned char)*M);
    memmove(p[i], _p[i], sizeof(unsigned char)*M);
    S0[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    s0[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    S1[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    s1[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    // FIXME check if volk_malloc is to be used.
    symbol[i] = (std::complex<float> *) malloc 
                (sizeof(std::complex<float>)*symbol_len);
    W[i] = (std::complex<float> *) volk_malloc 
           (sizeof(std::complex<float>)*M, volk_alignment);
    // initialize the precoding to unity
    std::fill(W[i], W[i] + M, I);
    X[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    x[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    assert(volk_is_aligned(W[i]));
    assert(volk_is_aligned(X[i]));
    assert(volk_is_aligned(x[i]));
    ifft[i] = fftwf_plan_dft_1d(M,
    reinterpret_cast<fftwf_complex *>(X[i]),
    reinterpret_cast<fftwf_complex *>(x[i]),
                                FFTW_BACKWARD,
                                FFTW_ESTIMATE);
    // sanity check for subcarrier allocation
    ofdmframe_validate_sctype(p[i],
                              M,
                              &M_null[i],
                              &M_pilot[i],
                              &M_data[i]);
    // initialize S0 sequence
    ofdmframe_init_S0(p[i],
                      M,
                      S0[i],
                      s0[i],
                      &M_S0[i],
                      _ms_S0[i]);
    // initialize S1 sequence
    ofdmframe_init_S1(p[i],
                      M,
                      S1[i],
                      s1[i],
                      &M_S1[i],
                      _ms_S1[i]);
  }
}

void
framegen::set_W(std::complex<float> ** const _W)
{
  // TODO sanity checks on precoding vector
  for(unsigned int i = 0; i < num_streams; i++) {
    memmove(W[i], _W[i], sizeof(std::complex<float>)*M);
  }
}

unsigned int
framegen::write_sync_words(std::vector<std::complex<float> *> tx_buff)
{
  assert(num_streams == tx_buff.size());
  unsigned int sample_index = 0;
  // total length = 2 x S0 + 1 x S1
  unsigned int total_count = (2 + 1)*(M + cp_len)*num_streams;
  // write 2 copies of S0 followed by one copy S1
  // to each stream
  for(unsigned int stream = 0; stream < num_streams; stream++) {
    memmove(tx_buff[stream] + sample_index,
            s0[stream] + M - 2*cp_len,
            sizeof(std::complex<float>)*2*cp_len);
    sample_index += 2*cp_len;
    memmove(tx_buff[stream] + sample_index,
            s0[stream],
            sizeof(std::complex<float>)*M);
    sample_index += M;
    memmove(tx_buff[stream] + sample_index,
            s0[stream],
            sizeof(std::complex<float>)*M);
    sample_index += M;
    // add cyclic prefix
    memmove(tx_buff[stream] + sample_index,
            s1[stream] + M - cp_len,
            sizeof(std::complex<float>)*cp_len);
    sample_index += cp_len;
    // write ofdm symbol
    memmove(tx_buff[stream] + sample_index,
            s1[stream],
            sizeof(std::complex<float>)*M);
    sample_index += M;
  }
  // total count of sync samples
  assert(sample_index == total_count);
  return sample_index;
}

framegen::~framegen()
{
  for(unsigned int i = 0; i < num_streams; i++)
  {
    free(p[i]);
    free(S0[i]);
    free(s0[i]);
    free(S1[i]);
    free(s1[i]);
    free(symbol[i]);
    volk_free(W[i]);
    fftwf_free(X[i]);
    fftwf_free(x[i]);
    fftwf_destroy_plan(ifft[i]);
  }
  free(W);
}

void framegen::print()
{
    printf("ofdmframegen:\n");
    printf("    num subcarriers     :   %-u\n", M);
    printf("      - NULL            :   %-u\n", M_null[0]);
    printf("      - pilot           :   %-u\n", M_pilot[0]);
    printf("      - data            :   %-u\n", M_data[0]);
    printf("    cyclic prefix len   :   %-u\n", cp_len);
    printf("    ");
    ofdmframe_print_sctype(p[0], M);
}

framesync::framesync(unsigned int _M,
                     unsigned int _cp_len,
                     unsigned int _num_streams,
                     std::vector<unsigned char *> const &_p,
                     std::vector<msequence> const &_ms_S0,
                     std::vector<msequence> const &_ms_S1,
                     std::vector<msequence> const &_ms_pilot,
                     zf_callback _callback)
{
  // sanity checks for parameters
  assert(_p.size() == _num_streams);
  assert(_ms_S0.size() == _num_streams);
  assert(_ms_S1.size() == _num_streams);
  assert(_ms_pilot.size() == _num_streams);

  // assign values for variables
  M = _M;
  cp_len = _cp_len;
  symbol_len = M + cp_len;
  num_streams = _num_streams;
  callback = _callback;
  M2 = M/2;

  // resize vectors to num streams
  p.resize(num_streams);
  M_null.resize(num_streams);
  M_pilot.resize(num_streams);
  M_data.resize(num_streams);
  M_S0.resize(num_streams);
  M_S1.resize(num_streams);
  g_data.resize(num_streams);
  fft.resize(num_streams);
  X.resize(num_streams);
  x.resize(num_streams);
  symbol.resize(num_streams);
  S0.resize(num_streams);
  s0.resize(num_streams);
  S1.resize(num_streams);
  s1.resize(num_streams);
  ms_pilot.resize(num_streams);
  g0.resize(num_streams);
  G0.resize(num_streams);
  G1.resize(num_streams);
  G = (std::complex<float> **) malloc 
      (sizeof(std::complex<float>*)*num_streams);
  B.resize(num_streams);
  R.resize(num_streams);

  std::complex<float> Z(0, 0);
  for(unsigned int i = 0; i < num_streams; i++)
  {
    // memory allocation
    // NOTE: all to be freed in destructor
    p[i]  = (unsigned char *) malloc
            (sizeof(unsigned char)*M);
    memmove(p[i], _p[i], sizeof(unsigned char)*M);
    S0[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    s0[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    S1[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    s1[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    // FIXME volk_malloc is to be used?
    symbol[i] = (std::complex<float> *) malloc 
                (sizeof(std::complex<float>)*symbol_len);
    X[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    x[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    G0[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    G1[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*M);
    G[i] = (std::complex<float> *) malloc 
           (sizeof(std::complex<float>)*M);
    B[i] = (std::complex<float> *) malloc 
           (sizeof(std::complex<float>)*M);
    R[i] = (std::complex<float> *) malloc 
           (sizeof(std::complex<float>)*M);
    std::fill(G0[i], G0[i] + M, Z);
    std::fill(G1[i], G1[i] + M, Z);
    std::fill(G[i], G[i] + M, Z);
    std::fill(B[i], B[i] + M, Z);

    // timing backoff FIXME what is this??
    backoff = cp_len < 2 ? cp_len : 2;
    float phi = (float)(backoff)*2.0f*M_PI/(float)(M);
    unsigned int n;
    for (n = 0; n < M; n++)
        B[i][n] = std::polar(1.0f, i*phi);

    // initialize fft object
    fft[i] = fftwf_plan_dft_1d(M,
    reinterpret_cast<fftwf_complex *>(x[i]),
    reinterpret_cast<fftwf_complex *>(X[i]),
                               FFTW_FORWARD,
                               FFTW_ESTIMATE);
    // sanity check for subcarrier allocation
    ofdmframe_validate_sctype(p[i],
                              M,
                              &M_null[i],
                              &M_pilot[i],
                              &M_data[i]);
    // initialize S0 sequence
    ofdmframe_init_S0(p[i],
                      M,
                      S0[i],
                      s0[i],
                      &M_S0[i],
                      _ms_S0[i]);
    // initialize S1 sequence
    ofdmframe_init_S1(p[i],
                      M,
                      S1[i],
                      s1[i],
                      &M_S1[i],
                      _ms_S1[i]);
  }
  // numerically-controlled oscillator
  nco_rx = nco_crcf_create(LIQUID_NCO);
  input_buffer = windowcf_create((M + cp_len));
  sync_index = 0;
  num_samples_processed = 0;

  // volk malloc
  size_t volk_alignment = volk_get_alignment();
  volk_buff_fc1 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*
                   (M + cp_len)*3*num_streams,
                   volk_alignment);
  volk_buff_fc2 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*
                   (M + cp_len)*3*num_streams,
                   volk_alignment);
  volk_buff_fc3 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*
                   (M + cp_len)*3*num_streams,
                   volk_alignment);
  volk_buff_f1 = (float *) volk_malloc 
                 (sizeof(float)*
                  (M + cp_len)*3*num_streams,
                  volk_alignment);
  reset();
}

void framesync::print()
{
    printf("ofdmframegen:\n");
    printf("    num subcarriers     :   %-u\n", M);
    printf("      - NULL            :   %-u\n", M_null[0]);
    printf("      - pilot           :   %-u\n", M_pilot[0]);
    printf("      - data            :   %-u\n", M_data[0]);
    printf("    cyclic prefix len   :   %-u\n", cp_len);
    printf("    ");
    ofdmframe_print_sctype(p[0], M);
}

unsigned long int framesync::get_sync_index()
{
  return sync_index;
}

void framesync::get_G(std::complex<float> ** _G)
{
  for(unsigned int i = 0; i < num_streams; i++){
    memmove(_G[i], G[i],
            sizeof(std::complex<float>)*M);
  }
}

void framesync::reset()
{
  nco_crcf_reset(nco_rx);
  // reset timers
  timer = 0;
  num_symbols = 0;
  s_hat_0 = 0.0f;
  s_hat_1 = 0.0f;
  phi_prime = 0.0f;
  p1_prime = 0.0f;

  // set thresholds (increase for small number of subcarriers)
  // FIXME what is this??
//  plcp_detect_thresh = (M > 44) ? 0.35f : 0.35f + 0.01f*(44 - M);
//  plcp_sync_thresh   = (M > 44) ? 0.30f : 0.30f + 0.01f*(44 - M);
  plcp_detect_thresh = 0.35;
  plcp_sync_thresh   = 0.35;

  // reset state
  state = STATE_SEEK_SYNC;
}

unsigned long long int framesync::get_num_samples_processed()
{
  return num_samples_processed;
}

// NOTE: processes only stream in_buff[0]
framesync_states_t
framesync::execute(std::vector<std::complex<float>*>
                   const &in_buff,
                   unsigned int num_samples)
{
  unsigned int i = 0;
  std::complex<float> x;

  // TODO process a vector of inputs, rather than
  // one by one in a for loop.
  for(; i < num_samples; i++){
    x = in_buff[0][i];
    // FIXME correct for frequency offset
    if(state != STATE_SEEK_SYNC) {
      nco_crcf_mix_down(nco_rx, x, &x);
      nco_crcf_step(nco_rx);
    }

    // save input to buffer
    windowcf_push(input_buffer,x);
    num_samples_processed++;

    switch(state) {
      case STATE_SEEK_SYNC:
        execute_seek_sync();
        break;
      case STATE_SYNCING:
        execute_syncing();
        break;
      case STATE_ZF:
        execute_zf();
        break;
      default:
        std::cout << "Unknown state, exiting\n";
        exit(1);
    }
  }
  return state;
}

void framesync::execute_seek_sync()
{
  timer++;
  sync_index++;

  if (timer < M)
    return;

  // reset timer
  timer = 0;

  //
  std::complex<float> * rc;
  windowcf_read(input_buffer, &rc);

// copied from original source
    unsigned int i;
    float g = 0.0f;
    for (i=cp_len; i<M + cp_len; i++) {
        // compute |rc[i]|^2 efficiently
        g += std::real(rc[i])*std::real(rc[i]) + 
             std::imag(rc[i])*std::imag(rc[i]);
    }
    g = (float)(M) / g;

// try the code below if this works 

//  // estimate gain
//  memmove(volk_buff_fc1, rc,
//          sizeof(std::complex<float>)*(M + cp_len));
//  // TODO Write a volk kernel to do this more efficienlty
//  volk_32fc_x2_conjugate_dot_prod_32fc(volk_buff_fc3,
//                                       volk_buff_fc1,
//                                       volk_buff_fc1,
//                                       M + cp_len);
//  float g = std::real(volk_buff_fc3[0]);
//
//  // FIXME I believe it should be g / M
//  g = (float)(M) / g;
//  g = g / (float)(M);

  // FIXME squelch???

  // estimate S0 gain
  estimate_gain_S0(&rc[cp_len]);

  std::complex<float> s_hat;
  S0_metrics(&s_hat);
  s_hat *= g;

  float tau_hat  = std::arg(s_hat) * (float)(M2) / (2*M_PI);
  #if DEBUG_PRINT
    printf(" - gain=%12.3f, rssi=%12.8f, s_hat=%12.4f <%12.8f>, tau_hat=%8.3f\n",
            sqrt(g),
            -10*log10(g),
            std::abs(s_hat), std::arg(s_hat),
            tau_hat);
  #endif

    // save gain (permits dynamic invocation of get_rssi() method)
    g0[0] = g;

    // 
    if (std::abs(s_hat) > plcp_detect_thresh) {

        int dt = (int)roundf(tau_hat);
        // set timer appropriately...
        timer = (M + dt) % (M2);
        timer += M; // add delay to help ensure good S0 estimate
        state = STATE_SYNCING;

#if DEBUG_PRINT
        printf("********** frame detected! ************\n");
        printf("    s_hat   :   %12.8f <%12.8f>\n", std::abs(s_hat),
            std::arg(s_hat));
        printf("  tau_hat   :   %12.8f\n", tau_hat);
        printf("    dt      :   %12d\n", dt);
        printf("    timer   :   %12u\n", timer);
#endif
    }
}

void framesync::estimate_gain_S0(std::complex<float> * _x)
{
    // move input array into fft input buffer
    memmove(x[0], _x, M*sizeof(std::complex<float>));

    // compute fft, storing result into _q->X
    fftwf_execute(fft[0]);
    
    // compute gain, ignoring NULL subcarriers
    unsigned int i;
    float gain = sqrtf(M_S0[0]) / (float)(M);

    // FIXME do this with volk
    for (i=0; i<M; i++) {
        if (p[0][i] != OFDMFRAME_SCTYPE_NULL && (i%2)==0) {
        // NOTE : if cabsf(_q->S0[i]) == 0 then we can multiply 
        // by conjugate rather than compute division
        //_G[i] = _q->X[i] / _q->S0[i];
            G0[0][i] = X[0][i] * std::conj(S0[0][i]);
        }
        else {
            G0[0][i] = 0.0f;
        }

        // normalize gain
        G0[0][i] *= gain;
    }
}

void framesync::S0_metrics(std::complex<float> * _s_hat)
{
    // timing, carrier offset correction
    unsigned int i;
    std::complex<float> s_hat = 0.0f;

    // compute timing estimate, accumulate phase difference across
    // gains on subsequent pilot subcarriers (note that all the odd
    // subcarriers are NULL)
    for (i=0; i<M; i+=2) {
        s_hat += G0[0][(i+2)%M]*std::conj(G0[0][i]);
    }
    s_hat /= M_S0[0]; // normalize output

    // set output values
    *_s_hat = s_hat;
}

void framesync::execute_syncing()
{
  ;
}

void framesync::execute_zf()
{
  ;
}

framesync::~framesync()
{
  for(unsigned int i = 0; i < num_streams; i++)
  {
    free(S0[i]);
    free(s0[i]);
    free(S1[i]);
    free(s1[i]);
    free(symbol[i]);
    fftwf_free(X[i]);
    fftwf_free(x[i]);
    fftwf_destroy_plan(fft[i]);
    free(G0[i]);
    free(G1[i]);
    free(G[i]);
    free(B[i]);
    free(R[i]);
  }
  // numerically-controlled oscillator
  nco_crcf_destroy(nco_rx);
  windowcf_destroy(input_buffer);
  free(G);
  volk_free(volk_buff_fc1);
  volk_free(volk_buff_fc2);
  volk_free(volk_buff_fc3);
  volk_free(volk_buff_f1);
}

void ofdmframe_init_default_sctype(unsigned char * _p, unsigned int _M)
{
    // validate input
    if (_M < 6) {
        fprintf(stderr,"warning: ofdmframe_init_default_sctype(), less than 4 subcarriers\n");
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

void ofdmframe_validate_sctype(const unsigned char * _p,
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
            fprintf(stderr,"error: ofdmframe_validate_sctype(), invalid subcarrier type (%u)\n", _p[i]);
            exit(1);
        }
    }

    // set outputs
    *_M_null  = M_null;
    *_M_pilot = M_pilot;
    *_M_data  = M_data;
}

void ofdmframe_print_sctype(const unsigned char * _p, unsigned int _M)
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
            fprintf(stderr,"error: ofdmframe_print_default_sctype(), invalid subcarrier type\n");
            exit(1);
        }
    }

    printf("]\n");
}

void ofdmframe_init_S0(const unsigned char * _p,
                       unsigned int _M,
                       std::complex<float> * _S0,
                       std::complex<float> * _s0,
                       unsigned int *  _M_S0,
                       msequence ms)
{
    unsigned int i;

    unsigned int s;
    unsigned int M_S0 = 0;

    // short sequence
    for (i=0; i<_M; i++) {
        // generate symbol
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
        fprintf(stderr,"error: ofdmframe_init_S0(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // set return value(s)
    *_M_S0 = M_S0;

    // run inverse fft to get time-domain sequence
    // TODO make in independent of liquid
    fft_run(_M, _S0, _s0, LIQUID_FFT_BACKWARD, 0);

    // normalize time-domain sequence level
    // TODO Do this with volk
    float g = 1.0f / sqrtf(M_S0);
    for (i=0; i<_M; i++)
        _s0[i] *= g;
}

void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       unsigned int *  _M_S1,
                       msequence ms)
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
        fprintf(stderr,"error: ofdmframe_init_S1(), no subcarriers enabled; check allocation\n");
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
