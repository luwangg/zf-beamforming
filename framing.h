#ifndef FRAMING_H
#define FRAMING_H

#include <vector>
#include <fftw3.h>
#include <gnuradio/gr_complex.h>
#include <complex>
#include <liquid/liquid.h>

// callback
typedef void * (*zf_callback) (void *);

  // receiver state
typedef enum {
  STATE_SEEK_PLATEAU = 0,
  STATE_SAVE_ACCESS_CODES,
  STATE_SEEK_SYNC_CH1,
  STATE_PLCP_SHORT0_CH1,
  STATE_PLCP_SHORT1_CH1,
  STATE_PLCP_LONG_CH1,
  STATE_SEEK_SYNC_CH2,
  STATE_PLCP_SHORT0_CH2,
  STATE_PLCP_SHORT1_CH2,
  STATE_PLCP_LONG_CH2,
  STATE_ZF,
  STATE_WAIT
} framesync_states_t;

class framegen {
 private:
  // number of subcarriers
  unsigned int M;
  // cyclic prefix length
  unsigned int cp_len;
  // symbol length
  unsigned int symbol_len;
  // number of data streams
  unsigned int num_streams;
  // number of access codes
  unsigned int num_access_codes;
  // subcarrier allocation for each channel
  std::vector<unsigned char *> p;
  // number of null subcarriers for each channel
  std::vector<unsigned int> M_null;
  // number of pilot subcarriers for each channel
  std::vector<unsigned int> M_pilot;
  // number of data subcarriers for each channel
  std::vector<unsigned int> M_data;
  // number of enabled subcarriers in S0
  std::vector<unsigned int> M_S0;
  // number of enabled subcarriers in S1
  std::vector<unsigned int> M_S1;
  // frequency domain buffer
  std::vector<std::complex<float> *> X;
  // time domain buffer
  std::vector<std::complex<float> *> x;
  // scaling factors
  std::vector<float> g_data;
  // transmit beamformer
  std::complex<float> ** W;
  // transform object
  std::vector<fftwf_plan> ifft;
  // symbol
  std::vector<std::complex<float> *> symbol;
  // PLCP short
  std::vector<std::complex<float> *> S0;
  std::vector<std::complex<float> *> s0;
  // PLCP long
  std::vector<std::complex<float> *> S1;
  std::vector<std::complex<float> *> s1;
  // pilot sequences
  std::vector<msequence> ms_pilot;

  // volk buffer for volk operations
  std::complex<float> * volk_buff_fc1;
  std::complex<float> * volk_buff_fc2;
  std::complex<float> * volk_buff_fc3;
  float * volk_buff_f1;

 public:
  // constructor
  // _M: number of ofdm subcarriers
  // _cp_len: cyclic prefix length
  // _num_streams: number of data streams
  // p: subcarrier allocation vectors
  // _ms_S0: vector of pn sequence generator to generate
  // short PLCP sequence
  // _ms_S1: vector of pn sequence generator to generate
  // long PLCP sequence
  // _ms_pilot: vector of pn sequence generator to generate
  // pilot sequences.
  framegen(unsigned int _M,
           unsigned int _cp_len,
           unsigned int _num_streams,
           unsigned int _num_access_codes,
           std::vector<unsigned char *> const &_p,
           std::vector<msequence> const &_ms_S0,
           std::vector<msequence> const &_ms_S1,
           std::vector<msequence> const &_ms_pilot);
  // destructor
  ~framegen();
  void print();

  unsigned int
  write_sync_words(std::vector<std::complex<float> *> tx_buff);

  // set the precoding vector.
  void set_W(std::complex<float> ** const _W);
  gr_complex ** get_W();
  void compute_W(gr_complex ** G);
  unsigned int
  write_zf_words_random(std::vector<gr_complex *> tx_buff);
  unsigned int
  write_words_random(std::vector<gr_complex *> tx_buff);
  unsigned int get_num_streams();
};

class framesync {
 private:
  // number of subcarriers
  unsigned int M;
  // number of subcarriers divided by 2
  unsigned int M2;
  // cyclic prefix length
  unsigned int cp_len;
  // symbol length
  unsigned int symbol_len;
  // number of data streams
  unsigned int num_streams;
  // number of access codes per stream
  unsigned int num_access_codes;
  // subcarrier allocation for each channel
  std::vector<unsigned char *> p;
  // number of null subcarriers for each channel
  std::vector<unsigned int> M_null;
  // number of pilot subcarriers for each channel
  std::vector<unsigned int> M_pilot;
  // number of data subcarriers for each channel
  std::vector<unsigned int> M_data;
  // number of enabled subcarriers in S0
  std::vector<unsigned int> M_S0;
  // number of enabled subcarriers in S1
  std::vector<unsigned int> M_S1;
  // scaling factors
  std::vector<float> g_data;
  // transform object
  std::vector<fftwf_plan> fft;
  // frequency domain buffer
  std::vector<std::complex<float> *> X;
  // time domain buffer
  std::vector<std::complex<float> *> x;
  // symbol
  std::vector<std::complex<float> *> symbol;
  // PLCP short
  std::vector<std::complex<float> *> S0;
  std::vector<std::complex<float> *> s0;
  // PLCP long
  std::vector<std::complex<float> *> S1;
  std::vector<std::complex<float> *> s1;
  // pilot sequences
  std::vector<msequence> ms_pilot;
  // input buffer
  windowcf input_buffer;
  // index at which sync is found
  unsigned long int sync_index;
  // total number of samples processed
  unsigned long long int num_samples_processed;

  // gain
  std::vector<float> g0;                // nominal gain
  std::vector<std::complex<float>*> G0; // estimate from G0
  std::vector<std::complex<float>*> G1; // estimate from G1
  std::complex<float> ** G;             // CSI available in framesync 
  std::vector<std::complex<float>*> B;  // subcarrier phase rotation
                                        // due to backoff FIXME what?
  std::vector<std::complex<float>*> R;  // FIXME what?

  // synchronizer objects
  // FIXME only one of these??
  nco_crcf nco_rx;        // numerically-controlled oscillator
  float phi_prime;        // ...
  float p1_prime;         // filtered pilot phase slope
  float nco_freq;
  // timing
  unsigned int timer;         // input sample timer
  unsigned int num_symbols;   // symbol counter
  unsigned int backoff;       // sample timing backoff
  std::vector<gr_complex> s_hat_0;// first S0 symbol metrics estimate
  std::vector<gr_complex> s_hat_1;// second S0 symbol metrics estimate

  // detection thresholds
  // plcp detection threshold, nominally 0.35
  float plcp_detect_thresh;
  // long symbol threshold, nominally 0.30
  float plcp_sync_thresh;

  // callback
  zf_callback callback;
  void * userdata;

  // receiver state
  framesync_states_t state;

  // reset all
  void reset();
  // internals
  void execute_seek_sync(unsigned int chan);
  void execute_sc_sync(gr_complex _x);
  void execute_save_access_codes(gr_complex _x);
  void estimate_channel();
  void print_estimates();

  // volk buffer for volk operations
  std::complex<float> * volk_buff_fc1;
  std::complex<float> * volk_buff_fc2;
  std::complex<float> * volk_buff_fc3;
               float  * volk_buff_f1;

  // Schmidl & Cox
  wdelaycf delay;
  unsigned int input_buffer_len;
  firfilt_crcf delay_ma;
  firfilt_rrrf normalizer_ma;
  gr_complex delay_conjugate;
  gr_complex delay_corr;
  float delay_magsquare;
  float delay_normalize;
  float normalizer_magsquare;
  float normalizer_square;
  float peak_to_angle;
  unsigned long int plateau_start;
  unsigned long int plateau_end;
  bool in_plateau;
  FILE * f_sc_debug_out;
  FILE * f_s0_corr;
  FILE * corr_ac11;
  FILE * corr_ac12;
  FILE * corr_ac13;
  FILE * corr_ac21;
  FILE * corr_ac22;
  FILE * corr_ac23;
  float * corr_buff_s0;
  float * corr_buff_ac11;
  float * corr_buff_ac12;
  float * corr_buff_ac13;
  float * corr_buff_ac21;
  float * corr_buff_ac22;
  float * corr_buff_ac23;
  unsigned int index_s0;
  unsigned int index_ac11;
  unsigned int index_ac12;
  unsigned int index_ac13;
  unsigned int index_ac21;
  unsigned int index_ac22;
  unsigned int index_ac23;
  gr_complex * estimate_11;
  gr_complex * estimate_12;
  gr_complex * estimate_13;
  gr_complex * estimate_21;
  gr_complex * estimate_22;
  gr_complex * estimate_23;
  gr_complex * estimate_1;
  gr_complex * estimate_2;

 public:
  // constructor
  // _M: number of ofdm subcarriers
  // _cp_len: cyclic prefix length
  // _num_streams: number of data streams
  // p: subcarrier allocation vectors
  // _ms_S0: vector of pn sequence generator to generate
  // short PLCP sequence
  // _ms_S1: vector of pn sequence generator to generate
  // long PLCP sequence
  // _ms_pilot: vector of pn sequence generator to generate
  // pilot sequences.
  framesync(unsigned int _M,
            unsigned int _cp_len,
            unsigned int _num_streams,
            unsigned int _num_access_codes,
            std::vector<unsigned char *> const &_p,
            std::vector<msequence> const &_ms_S0,
            std::vector<msequence> const &_ms_S1,
            std::vector<msequence> const &_ms_pilot,
            zf_callback callback); 
  // destructor
  ~framesync();
  // method to expose CSI
  void get_G(std::complex<float> ** _G);
  void print();
  unsigned long int get_plateau_start();
  unsigned long int get_plateau_end();

  unsigned long int get_sync_index();
  unsigned long long int get_num_samples_processed();

  // seek PLCP, sync and compute CSI.
  // run callback to inform the caller thread.
  framesync_states_t 
  execute(std::vector<std::complex<float>*>
          const &in_buff,
          unsigned int num_samples);
};

// initialize default subcarrier allocation
//  _M      :   number of subcarriers
//  _p      :   output subcarrier allocation vector, [size: _M x 1]
//
// key: '.' (null), 'P' (pilot), '+' (data)
// .+++P+++++++P.........P+++++++P+++
//
void ofdmframe_init_default_sctype(unsigned char *_p, unsigned int _M);

// validate subcarrier type (count number of null, pilot, and data
// subcarriers in the allocation)
//  _p          :   subcarrier allocation array, [size: _M x 1]
//  _M_null     :   output number of null subcarriers
//  _M_pilot    :   output number of pilot subcarriers
//  _M_data     :   output number of data subcarriers
void ofdmframe_validate_sctype(const unsigned char * _p,
                               unsigned int _M,
                               unsigned int * _M_null,
                               unsigned int * _M_pilot,
                               unsigned int * _M_data);

// print subcarrier allocation to screen
//
// key: '.' (null), 'P' (pilot), '+' (data)
// .+++P+++++++P.........P+++++++P+++
//
void ofdmframe_print_sctype(const unsigned char * _p, unsigned int _M);

// generate short sequence symbols
//  _p      :   subcarrier allocation array
//  _S0     :   output symbol (freq)
//  _s0     :   output symbol (time)
//  _M_S0   :   total number of enabled subcarriers in S0
//  ms      :   pn_sequence generator to generate S0
void ofdmframe_init_S0(const unsigned char * _p,
                       unsigned int _M,
                       std::complex<float> * _S0,
                       std::complex<float> * _s0,
                       unsigned int *  _M_S0,
                       msequence ms);

// generate long sequence symbols
//  _p      :   subcarrier allocation array
//  _S1     :   output symbol (freq)
//  _s1     :   output symbol (time)
//  _M_S1   :   total number of enabled subcarriers in S1
//  ms      :   pn_sequence generator to generate S1
void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       unsigned int *  _M_S1,
                       msequence ms);

gr_complex liquid_cexpjf(float theta);
float cabsf(gr_complex z);
float cargf(gr_complex z);
float fabsf(float x);
gr_complex conjf(gr_complex z);

#endif // FRAMING_H
