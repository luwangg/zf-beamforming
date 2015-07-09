#include "framing.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <volk/volk.h>
#include <volk/volk_malloc.h>
#include <string.h>

#define DEBUG_PRINT true
#define USE_NCO     false
#define PLATEAU_THREASHOLD  0.95
#define ADD_NULL_CARRIERS   false
#define FILE_SC_OUT         "/tmp/f_sc_out"
#define FILE_s0_CORR        "/tmp/f_s0_corr"
#define DEBUG_LOG           true

#define FILE_AC11  "/tmp/f_ac11"
#define FILE_AC12  "/tmp/f_ac12"
#define FILE_AC13  "/tmp/f_ac13"
#define FILE_AC21  "/tmp/f_ac21"
#define FILE_AC22  "/tmp/f_ac22"
#define FILE_AC23  "/tmp/f_ac23"

gr_complex BPSK_CONSTELLATION[] = 
{
  gr_complex(-1.0, 0.0),
  gr_complex(1.0,  0.0),
};
#define BPSK_CONSTELLATION_SIZE 2
gr_complex QPSK_CONSTELLATION[] = 
{
  gr_complex(0.5, 0.5),
  gr_complex(-0.5, 0.5),
  gr_complex(-0.5, -0.5),
  gr_complex(0.5, -0.5)
};
#define QPSK_CONSTELLATION_SIZE 4
#define CONSTELLATION BPSK_CONSTELLATION
#define CONSTELLATION_SIZE BPSK_CONSTELLATION_SIZE

inline gr_complex
liquid_cexpjf(float theta)
{
  return std::polar(1.0f, theta);
}

inline float
cabsf(gr_complex z)
{
  return std::abs(z);
}

inline float
cargf(gr_complex z)
{
  return std::arg(z);
}

inline float
fabsf(float x)
{
  return std::abs(x);
}

inline gr_complex
conjf(gr_complex z)
{
  return std::conj(z);
}

framegen::framegen(unsigned int _M,
                   unsigned int _cp_len,
                   unsigned int _num_streams,
                   unsigned int _num_access_codes,
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
  num_access_codes = _num_access_codes;

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
            (sizeof(std::complex<float>)*
             M*num_access_codes);
    s1[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*
             M*num_access_codes);
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
                      num_access_codes,
                      S1[i],
                      s1[i],
                      &M_S1[i],
                      _ms_S1[i]);
    g_data[i] = 1.0f / sqrtf(M_pilot[i] + M_data[i]);
  }
  volk_buff_fc1 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_fc2 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_fc3 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_f1 = (float *) volk_malloc 
                 (sizeof(float)*
                  (M + cp_len)*3*num_streams,
                  volk_alignment);
}

unsigned int
framegen::get_num_streams()
{
  return num_streams;
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
  unsigned int total_count = (num_access_codes*num_streams + 1)*
                             (M + cp_len);
  // set tx_buff to 0;
  for(unsigned int stream = 0; stream < num_streams; stream++) {
    std::fill(tx_buff[stream],
              tx_buff[stream] + total_count,
              gr_complex(0.0, 0.0));
  }
  // write S0 onto ch0 and then S1 in TDMA
  memmove(tx_buff[0] + sample_index,
          s0[0] + M - cp_len,
          sizeof(std::complex<float>)*cp_len);
  sample_index += cp_len;
  memmove(tx_buff[0] + sample_index,
          s0[0],
          sizeof(std::complex<float>)*M);
  sample_index += M;
  for(unsigned int ac_id = 0; ac_id < num_access_codes; ac_id++) {
    for(unsigned int stream = 0; stream < num_streams; stream++) {
      // add cyclic prefix
      memmove(tx_buff[stream] + sample_index,
              s1[stream] + M*ac_id + M - cp_len,
              sizeof(std::complex<float>)*cp_len);
      sample_index += cp_len;
      // write ofdm symbol
      memmove(tx_buff[stream] + sample_index,
              s1[stream] + M*ac_id,
              sizeof(std::complex<float>)*M);
      sample_index += M;
    }
  }
  // total count of sync samples
  assert(sample_index == total_count);
  return sample_index;
}

void framegen::compute_W(gr_complex ** G)
{
  // compute the stronger stream
  std::vector<float> g;
  for(unsigned int stream = 0; stream < num_streams; stream++) {
    volk_32fc_magnitude_squared_32f(volk_buff_f1,
				    G[stream],
				    M);
    g.push_back(0.0f);
    for(unsigned int c = 0; c < M; c++) {
      g[stream] += volk_buff_f1[c];
    }
    g[stream] /= M_S1[stream];
  }
  // the following code assumes num_streams = 2
  // FIXME "fft shift"???
  if(g[0] > g[1]) { // channel 0 stronger
    // assign W[1][c] = 1;
    for(unsigned int c = 0; c < M; c++) {
      assert(p[0][c] == p[1][c]);
      if(p[0][c] == OFDMFRAME_SCTYPE_NULL) {
	W[0][c] = 0.0f;
	W[1][c] = 0.0f;
      }
      else {
	W[1][c] = 1.0f;
	W[0][c] = -(G[1][c] / G[0][c]);
	//W[0][c] = 0.0f;
	//W[1][c] = 0.0f;
      }
    }
  }
  else if(g[1] > g[0]) { // channel 1 stronger
    // assign W[0][c] = 1;
    for(unsigned int c = 0; c < M; c++) {
      assert(p[0][c] == p[1][c]);
      if(p[0][c] == OFDMFRAME_SCTYPE_NULL) {
	W[0][c] = 0.0f;
	W[1][c] = 0.0f;
      }
      else {
	W[0][c] = 1.0f;
	W[1][c] = -(G[0][c] / G[1][c]);
	//W[0][c] = 0.0f;
	//W[1][c] = 0.0f;
      }
    }
  }
  else { // this should not occure
    std::cout << "go to hell, I am done with you\n";
    exit(1);
  }
  std::cout << "******* printing W ***********\n";
  for(unsigned int sc = 0; sc < M; sc++) {
    printf("%12.8f + %12.8fi\t%12.8f + %12.8fi\n",
           std::real(W[0][sc]),
           std::imag(W[0][sc]),
           std::real(W[1][sc]),
           std::imag(W[1][sc])
           );
  }
}

unsigned int
framegen::write_zf_words_random(std::vector<gr_complex *> tx_buff)
{
  assert(num_streams == tx_buff.size());
  for(unsigned int i = 0; i < M; i++) { 
    volk_buff_fc2[i] = CONSTELLATION[rand()%CONSTELLATION_SIZE];
  }
  for(unsigned int chan = 0; chan < num_streams; chan++) {
    volk_32fc_s32fc_multiply_32fc(volk_buff_fc1,
                                  volk_buff_fc2,
                                  gr_complex(g_data[chan], 0.0f),
                                  M);
    volk_32fc_x2_multiply_32fc(X[chan],
                               volk_buff_fc1,
                               W[chan],
                               M);
    fftwf_execute(ifft[chan]);
    // cyclic prefix
    memmove(tx_buff[chan], x[chan] + M - cp_len,
            sizeof(gr_complex)*cp_len);
    memmove(tx_buff[chan] + cp_len,
            x[chan],
            sizeof(gr_complex)*M);
  }
  return M + cp_len;
}

unsigned int
framegen::write_words_random(std::vector<gr_complex *> tx_buff)
{
  assert(num_streams == tx_buff.size());
  for(unsigned int i = 0; i < M; i++) { 
    volk_buff_fc2[i] = CONSTELLATION[rand()%CONSTELLATION_SIZE];
  }
  for(unsigned int chan = 0; chan < num_streams; chan++) {
    volk_32fc_s32fc_multiply_32fc(X[chan],
                                  volk_buff_fc2,
                                  gr_complex(g_data[chan], 0.0f),
                                  M);
    fftwf_execute(ifft[chan]);
    // cyclic prefix
    memmove(tx_buff[chan], x[chan] + M - cp_len,
            sizeof(gr_complex)*cp_len);
    memmove(tx_buff[chan] + cp_len,
            x[chan],
            sizeof(gr_complex)*M);
  }
  return M + cp_len;
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
  volk_free(volk_buff_fc1);
  volk_free(volk_buff_fc2);
  volk_free(volk_buff_fc3);
  volk_free(volk_buff_f1);
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
                     unsigned int _num_access_codes,
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
  num_access_codes = _num_access_codes;
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
  s_hat_0.resize(num_streams);
  s_hat_1.resize(num_streams);

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
            (sizeof(std::complex<float>)*
             M*num_access_codes);
    s1[i] = (std::complex<float> *) malloc 
            (sizeof(std::complex<float>)*
             M*num_access_codes);
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
                      num_access_codes,
                      S1[i],
                      s1[i],
                      &M_S1[i],
                      _ms_S1[i]);
  }
  // numerically-controlled oscillator
  nco_rx = nco_crcf_create(LIQUID_NCO);
  input_buffer_len = (M + cp_len)*(num_access_codes*num_streams + 4);
  input_buffer = windowcf_create(input_buffer_len);
  sync_index = 0;
  num_samples_processed = 0;

  // Schmidl & Cox
  plateau_start = 0;
  plateau_end = 0;
  delay = wdelaycf_create(M2);
  float taps_delay_ma[M2];
  float taps_normalizer_ma[M];
  gr_complex taps_correlate_s0[M];
  for(unsigned int i = 0; i < M2; i++)
    taps_delay_ma[i] = -1.0f; // FIXME why not 1.0f
  for(unsigned int i = 0; i < M; i++)
    taps_normalizer_ma[i] = 0.5f;
  for(unsigned int i = 0; i < M; i++)
    taps_correlate_s0[i] = conj(s0[0][i]);
  delay_ma = firfilt_crcf_create(taps_delay_ma, M2);
  normalizer_ma = firfilt_rrrf_create(taps_normalizer_ma, M);

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
  corr_buff_s0 = (float *) malloc (sizeof(float)*(input_buffer_len - M));
  corr_buff_ac11 = (float *) malloc (sizeof(float)*(input_buffer_len - M));
  corr_buff_ac12 = (float *) malloc (sizeof(float)*(input_buffer_len - M));
  corr_buff_ac13 = (float *) malloc (sizeof(float)*(input_buffer_len - M));
  corr_buff_ac21 = (float *) malloc (sizeof(float)*(input_buffer_len - M));
  corr_buff_ac22 = (float *) malloc (sizeof(float)*(input_buffer_len - M));
  corr_buff_ac23 = (float *) malloc (sizeof(float)*(input_buffer_len - M));

  estimate_11 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  estimate_12 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  estimate_13 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  estimate_21 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  estimate_22 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  estimate_23 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  estimate_1  = (gr_complex *) malloc (sizeof(gr_complex)*M);
  estimate_2  = (gr_complex *) malloc (sizeof(gr_complex)*M);

  if(DEBUG_LOG) {
   f_sc_debug_out = fopen(FILE_SC_OUT, "wb");
   f_s0_corr = fopen(FILE_s0_CORR, "wb");
   corr_ac11 = fopen(FILE_AC11, "wb");
   corr_ac12 = fopen(FILE_AC12, "wb");
   corr_ac13 = fopen(FILE_AC13, "wb");
   corr_ac21 = fopen(FILE_AC21, "wb");
   corr_ac22 = fopen(FILE_AC22, "wb");
   corr_ac23 = fopen(FILE_AC23, "wb");
  }
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
  s_hat_0[0] = 0.0f;
  s_hat_1[0] = 0.0f;
  s_hat_0[1] = 0.0f;
  s_hat_1[1] = 0.0f;
  phi_prime = 0.0f;
  p1_prime = 0.0f;
  nco_freq = 0;

  // set thresholds (increase for small number of subcarriers)
  // FIXME what is this??
  plcp_detect_thresh = (M > 44) ? 0.35f : 0.35f + 0.01f*(44 - M);
  plcp_sync_thresh   = (M > 44) ? 0.30f : 0.30f + 0.01f*(44 - M);

  // reset state
  state = STATE_SEEK_PLATEAU;
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
  if(state == STATE_ZF)
    state = STATE_WAIT;

  // TODO process a vector of inputs, rather than
  // one by one in a for loop.
  for(; i < num_samples; i++){
    x = in_buff[0][i];

    // FIXME correct frequency offset
    switch(state) {
      case STATE_SEEK_PLATEAU:
	execute_sc_sync(x);
        break;
      case STATE_SAVE_ACCESS_CODES:
	execute_save_access_codes(x);
      case STATE_ZF:
	break;
      case STATE_WAIT:
        break;
      default:
        std::cout << "Unknown state, exiting\n";
        exit(1);
    }
    num_samples_processed++;
  }
  return state;
}

void framesync::execute_sc_sync(gr_complex _x)
{
  // save to input buffer
  windowcf_push(input_buffer, _x);
  gr_complex y, z;
  float w;
  wdelaycf_read(delay, &y);
  wdelaycf_push(delay, _x);
  firfilt_crcf_push(delay_ma, conj(y)*_x);
  firfilt_crcf_execute(delay_ma, &z);
  delay_magsquare = real(z*conj(z));
  normalizer_magsquare = real(_x*conj(_x));
  firfilt_rrrf_push(normalizer_ma, normalizer_magsquare);
  firfilt_rrrf_execute(normalizer_ma, &w);
  delay_normalize = delay_magsquare/(w*w);
  if(DEBUG_LOG) {
    fwrite(&delay_normalize, sizeof(float), 1, f_sc_debug_out);
  }
  if(delay_normalize > PLATEAU_THREASHOLD)
  {
    if(in_plateau)
      plateau_end = num_samples_processed;
    else {
      in_plateau = true;
      plateau_start = num_samples_processed;
    }
  }
  else {
    if((plateau_end - plateau_start > cp_len) && in_plateau)
      state = STATE_SAVE_ACCESS_CODES;
    in_plateau = false;
  }
}

void framesync::execute_save_access_codes(gr_complex _x)
{
  windowcf_push(input_buffer, _x);
  if(num_samples_processed - plateau_start == input_buffer_len - 
     M - cp_len) {
    state = STATE_ZF;
    printf("****** access codes saved ********\n");
    estimate_channel();
  }
}

void framesync::estimate_channel() {
  printf("***** locating s0 starting point ********\n");
  gr_complex * buff_ptr;
  gr_complex corr;
  float max_corr_s0   = 0.0;
  float max_corr_ac11 = 0.0;
  float max_corr_ac12 = 0.0;
  float max_corr_ac13 = 0.0;
  float max_corr_ac21 = 0.0;
  float max_corr_ac22 = 0.0;
  float max_corr_ac23 = 0.0;

  windowcf_read(input_buffer, &buff_ptr);
  // take M samples, compute FFT and correlate with
  // S0
  for(unsigned int i = 0; i < input_buffer_len - M; i++) {
    memmove(x[0], buff_ptr + i, sizeof(gr_complex)*M);
//   // find energy of samples in x[0]
//   volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
//					 x[0],
//					 x[0],
//					 M);
//   // scale to unit energy
//   float g = sqrtf(1/std::real(corr));
//   for(unsigned int k = 0; k < M; k++)
//     x[0][k] *= g;
    fftwf_execute(fft[0]);
    volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
					 X[0],
					 S0[0],
					 M);
    corr_buff_s0[i] = std::real(corr*std::conj(corr));
    if(corr_buff_s0[i] > max_corr_s0) {
      max_corr_s0 = corr_buff_s0[i];
      index_s0 = i;
    }
    volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
					 X[0],
					 S1[0] + 0*M,
					 M);
    corr_buff_ac11[i] = std::real(corr*std::conj(corr));
    if(corr_buff_ac11[i] > max_corr_ac11) {
      max_corr_ac11 = corr_buff_ac11[i];
      index_ac11 = i;
    }
    volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
					 X[0],
					 S1[0] + 1*M,
					 M);
    corr_buff_ac12[i] = std::real(corr*std::conj(corr));
    if(corr_buff_ac12[i] > max_corr_ac12) {
      max_corr_ac12 = corr_buff_ac12[i];
      index_ac12 = i;
    }
    volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
					 X[0],
					 S1[0] + 2*M,
					 M);
    corr_buff_ac13[i] = std::real(corr*std::conj(corr));
    if(corr_buff_ac13[i] > max_corr_ac13) {
      max_corr_ac13 = corr_buff_ac13[i];
      index_ac13 = i;
    }
    volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
					 X[0],
					 S1[1] + 0*M,
					 M);
    corr_buff_ac21[i] = std::real(corr*std::conj(corr));
    if(corr_buff_ac21[i] > max_corr_ac21) {
      max_corr_ac21 = corr_buff_ac21[i];
      index_ac21 = i;
    }
    volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
					 X[0],
					 S1[1] + 1*M,
					 M);
    corr_buff_ac22[i] = std::real(corr*std::conj(corr));
    if(corr_buff_ac22[i] > max_corr_ac22) {
      max_corr_ac22 = corr_buff_ac22[i];
      index_ac22 = i;
    }
    volk_32fc_x2_conjugate_dot_prod_32fc(&corr,
					 X[0],
					 S1[1] + 2*M,
					 M);
    corr_buff_ac23[i] = std::real(corr*std::conj(corr));
    if(corr_buff_ac23[i] > max_corr_ac23) {
      max_corr_ac23 = corr_buff_ac23[i];
      index_ac23 = i;
    }
  }
  fwrite(corr_buff_s0, sizeof(float), input_buffer_len - M, f_s0_corr);
  fwrite(corr_buff_ac11, sizeof(float), input_buffer_len - M, corr_ac11);
  fwrite(corr_buff_ac12, sizeof(float), input_buffer_len - M, corr_ac12);
  fwrite(corr_buff_ac13, sizeof(float), input_buffer_len - M, corr_ac13);
  fwrite(corr_buff_ac21, sizeof(float), input_buffer_len - M, corr_ac21);
  fwrite(corr_buff_ac22, sizeof(float), input_buffer_len - M, corr_ac22);
  fwrite(corr_buff_ac23, sizeof(float), input_buffer_len - M, corr_ac23);
  printf("%%%%   index_s0:  %d\n", index_s0);
  printf("%%%% index_ac11:  %d\n", index_ac11);
  printf("%%%% index_ac12:  %d\n", index_ac12);
  printf("%%%% index_ac13:  %d\n", index_ac13);
  printf("%%%% index_ac21:  %d\n", index_ac21);
  printf("%%%% index_ac22:  %d\n", index_ac22);
  printf("%%%% index_ac23:  %d\n", index_ac23);
  
  // estimate 11
  memmove(x[0], buff_ptr + index_ac11, sizeof(gr_complex)*M);
  fftwf_execute(fft[0]);
  for(unsigned int i = 0; i < M; i++) {
    if(p[0][i] == OFDMFRAME_SCTYPE_NULL)
      estimate_11[i] = 0.0f;
    else
      estimate_11[i] = X[0][i]/S1[0][0*M + i];
  }
  // estimate 12
  memmove(x[0], buff_ptr + index_ac12, sizeof(gr_complex)*M);
  fftwf_execute(fft[0]);
  for(unsigned int i = 0; i < M; i++) {
    if(p[0][i] == OFDMFRAME_SCTYPE_NULL)
      estimate_12[i] = 0.0f;
    else
      estimate_12[i] = X[0][i]/S1[0][1*M + i];
  }
  // estimate 13
  memmove(x[0], buff_ptr + index_ac13, sizeof(gr_complex)*M);
  fftwf_execute(fft[0]);
  for(unsigned int i = 0; i < M; i++) {
    if(p[0][i] == OFDMFRAME_SCTYPE_NULL)
      estimate_13[i] = 0.0f;
    else
      estimate_13[i] = X[0][i]/S1[0][2*M + i];
  }
  // estimate 21
  memmove(x[0], buff_ptr + index_ac21, sizeof(gr_complex)*M);
  fftwf_execute(fft[0]);
  for(unsigned int i = 0; i < M; i++) {
    if(p[1][i] == OFDMFRAME_SCTYPE_NULL)
      estimate_21[i] = 0.0f;
    else
      estimate_21[i] = X[0][i]/S1[1][0*M + i];
  }
  // estimate 22
  memmove(x[0], buff_ptr + index_ac22, sizeof(gr_complex)*M);
  fftwf_execute(fft[0]);
  for(unsigned int i = 0; i < M; i++) {
    if(p[1][i] == OFDMFRAME_SCTYPE_NULL)
      estimate_22[i] = 0.0f;
    else
      estimate_22[i] = X[0][i]/S1[1][1*M + i];
  }
  // estimate 23
  memmove(x[0], buff_ptr + index_ac23, sizeof(gr_complex)*M);
  fftwf_execute(fft[0]);
  for(unsigned int i = 0; i < M; i++) {
    if(p[1][i] == OFDMFRAME_SCTYPE_NULL)
      estimate_23[i] = 0.0f;
    else
      estimate_23[i] = X[0][i]/S1[1][2*M + i];
  }
  for(unsigned int i = 0; i < M; i++) {
    G[0][i] = estimate_11[i] +
              estimate_12[i] +
              estimate_13[i];
    G[1][i] = estimate_21[i] +
              estimate_22[i] +
              estimate_23[i];
  }
  callback(G);
  print_estimates();
}

void framesync::print_estimates() {
  printf("============ Channel Estimates ==============\n");
}

unsigned long int
framesync::get_plateau_start()
{
  return plateau_start;
}

unsigned long int
framesync::get_plateau_end()
{
  return plateau_end;
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
  
  // Schmidl & Cox
  wdelaycf_destroy(delay);
  firfilt_crcf_destroy(delay_ma);
  firfilt_rrrf_destroy(normalizer_ma);
  free(corr_buff_s0);
  free(corr_buff_ac11);
  free(corr_buff_ac12);
  free(corr_buff_ac13);
  free(corr_buff_ac21);
  free(corr_buff_ac22);
  free(corr_buff_ac23);
  free(estimate_11);
  free(estimate_12);
  free(estimate_13);
  free(estimate_21);
  free(estimate_22);
  free(estimate_23);
  free(estimate_1);
  free(estimate_2);

  if(DEBUG_LOG) {
    fclose(f_sc_debug_out);
    fclose(f_s0_corr);
    fclose(corr_ac11);
    fclose(corr_ac12);
    fclose(corr_ac13);
    fclose(corr_ac21);
    fclose(corr_ac22);
    fclose(corr_ac23);
  }
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
    unsigned int G = 0;
    if (ADD_NULL_CARRIERS) {
      G = _M / 10;
      if (G < 2) G = 2;
    }

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
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       unsigned int *  _M_S1,
                       msequence ms)
{
    unsigned int i, j;

    unsigned int s;
    unsigned int M_S1 = 0;

    // long sequence
    for(j = 0; j < _num_access_codes; j++) {
      for (i=0; i<_M; i++) {
        // generate symbol
        //s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            (_S1 + _M*j)[i] = 0.0f;
        } else {
            (_S1 + _M*j)[i] = s ? 1.0f : -1.0f;
            M_S1++;
        }
      }
      // run inverse fft to get time-domain sequence
      fft_run(_M, _S1 + _M*j, _s1 + _M*j, LIQUID_FFT_BACKWARD, 0);
    }

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S1/_num_access_codes);
    for (i=0; i<_M*_num_access_codes; i++)
        _s1[i] *= g;
    *_M_S1 = M_S1/_num_access_codes;
}
