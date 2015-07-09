/*
 * Copyright (c) 2014, 2015 Manu T S
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <uhd/usrp/multi_usrp.hpp>
#include <iostream>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <ctime>
#include <liquid/liquid.h>
#include <liquid/ofdm.h>
#include <volk/volk.h>
#include <volk/volk_malloc.h>
#include <unistd.h>
#include "framing.h"
#include "zf_framing.h"

// define default configurations
// USRP identities
#define N200_12             "addr=134.147.118.212"
#define N200_15             "addr=134.147.118.215"
#define X300A               "addr=134.147.118.216"
#define X300B               "addr=134.147.118.217"
#define B210A               "serial=308F955"
#define B210B               "serial=308F965"

// subdevice specifications
#define X300_SUBDEV_SPEC    "A:0 B:0"
#define N200_SUBDEV_SPEC    "A:0"
#define B210_SUBDEV_SPEC_TX "A:B A:A"
#define B210_SUBDEV_SPEC_RX "A:A A:B"

// tx/rx streamer configurations
#define CPU                 "fc32"
#define WIRE                "sc16"

// RF front end configurations
#define CENTER_FREQUENCY    5100e6
#define SAMPLING_RATE       1e6
#define TX_FRONTEND_GAIN    45.0
#define RX_FRONTEND_GAIN    45.0

// device synchronization configurations
typedef enum
{
  CLOCK_SOURCE_NONE = 0,
  CLOCK_SOURCE_EXTERNAL
} clock_source_type;

typedef enum
{
  TIME_SOURCE_NONE = 0,
  TIME_SOURCE_EXTERNAL
} time_source_type;
#define CLOCK_SOURCE        CLOCK_SOURCE_EXTERNAL
#define TIME_SOURCE         TIME_SOURCE_NONE

// OFDM configurations
#define NUM_SUBCARRIERS     64
#define CP_LENGTH           16
#define BASEBAND_GAIN       0.25
#define TAPER_LEN           0
#define NUM_ACCESS_CODES    3

// misc configurations
#define VERBOSITY           true
#define RUNTIME             10.0
#define SIMULATION          false
#define DECODE_ONLINE       true
#define NUM_STREAMS         2
#define LOG                 true
#define LOG_DIR             "/tmp/"
#define RX_RUNTIME          5.0
#define OFDM                true
// generator polynomials obtained from
// primitive_polys.pdf
#define LFSR_SMALL_LENGTH   12
#define LFSR_LARGE_LENGTH   13
#define LFSR_SMALL_0_GEN_POLY 010123
#define LFSR_SMALL_1_GEN_POLY 010151
#define LFSR_LARGE_0_GEN_POLY 020033
#define LFSR_LARGE_1_GEN_POLY 020047

namespace po = boost::program_options;

// global variablds
unsigned int pid;
#define PID_MAX               100
unsigned int num_frames_detected;
unsigned int num_valid_headers_received;
unsigned int num_valid_bytes_received;
unsigned int num_valid_packets_received;
time_t tx_begin, tx_end, rx_begin, rx_end;
float time_for_offline_decoding;


int ofdm_callback(unsigned char *  _header,
             int              _header_valid,
             unsigned char *  _payload,
             unsigned int     _payload_len,
             int              _payload_valid,
             framesyncstats_s _stats,
             void *           _userdata)
{
  // update global counters
  num_frames_detected++;

  if (_header_valid)
      num_valid_headers_received++;

  if (_payload_valid) {
      num_valid_packets_received++;
      num_valid_bytes_received += _payload_len;
  }
  return 0;
}
// NOTE: rxmd is declared global so that main can check for
// error flags in rxmd
uhd::rx_metadata_t rxmd;

// NOTE: Thread Synchronization
// --------------------------------------------------------
// Tx starts sending only when start_tx is true. Tx locks 
// mutex and checks for start_tx. If false it goes into
// condition wait which automatically unlocks the mutex.
// Rx thread is responsible for signalling Tx to start.
// After initialization, rx starts streaming, rx locks the 
// mutex. It signals the condition, and unlocks the mutex. 
// On receiving the condition, tx unlocks the mutex, and
// start transmitting the synchronization symbols.
bool start_tx = false;
bool stop_rx  = false;
bool start_zf = false;
pthread_mutex_t mutex_txrx;
pthread_cond_t condition_txrx;
pthread_cond_t condition_zf;

// NOTE: CSI Matrix and CSI Feedback
//                  CSI Matrix
//                  ----------
// SC-i refers to the ith subcarrier.
//
//            SC-1  SC-2  SC-3    ......    SC-M
//
// stream1    g11   g12   g13     ......    g1M
// 
// stream2    g21   g22   g23     ......    g2M
//
//   .         .     .     .      ......     .
//
//   .         .     .     .      ......     .
//
//   .         .     .     .      ......     .
//
// streamN    gN1   gN2   gN3     ......    gNM
// --------------------------------------------------------
// Entry gij in the above matrix refers to the channel gain
// of the jth subcarrier on the link from the ith transmit 
// antenna. The framesync has its own copy of CSI matrix.
// On receiving the sync word, the framesync computes the 
// CSI and updates its copy of G. The receive thread to be
// informed of the arrival of the sync word through the 
// value returned from the execute method in the framesync.
// Once rx thread is informed about the arrival of the sync
// word it copies the G matrix from framesync object, and
// runs a routine to compute the precoding matrix.
// NOTE: Space to hold CSI matrix to be allocated and freed
// in the main thread.
std::complex<float> ** G;

// NOTE: Precoding and Precoding Matrix
//                Precoding Matrix
//                ----------------
// SC-i refers to the ith subcarrier.
//
//            SC-1  SC-2  SC-3    ......    SC-M
//
// stream1    w11   w12   w13     ......    w1M
// 
// stream2    w21   w22   w23     ......    w2M
//
//   .         .     .     .      ......     .
//
//   .         .     .     .      ......     .
//
//   .         .     .     .      ......     .
//
// streamN    wN1   wN2   wN3     ......    wNM
// --------------------------------------------------------
// Entry wij in the above matrix refers to the coefficient 
// which should be multiplied to the subsymbol on the jth 
// subcarrier from ith antenna. A pointer to W is passed
// to the framegen during intialization. The framegen keeps
// its own copy of W. The initial value of the precoding 
// coefficients is set to one. Once the rx thread aquires 
// CSI information from the framesync, it runs a routine 
// which updates the W matrix, and signals the tx thread
// regarding update in W. Tx thread then asks the framegen
// to update its copy of W.
// NOTE: Space to hold W to be allocated and freed in main
std::complex<float> ** W;

// Compute precoding vector given CSI
void * design_zf_precoder (void *)
{
  std::cout << "Design zero-forcing Precoder\n";
  // TODO update W
  return NULL;
}

// call back
void * callback (void * __G)
{
  std::cout << "********  callback  **********\n";
  memmove(G[0], ((gr_complex **)__G)[0], sizeof(gr_complex)*NUM_SUBCARRIERS);
  memmove(G[1], ((gr_complex **)__G)[1], sizeof(gr_complex)*NUM_SUBCARRIERS);
  for(unsigned int sc = 0; sc < NUM_SUBCARRIERS; sc++) {
    printf("%12.8f + %12.8fi\t%12.8f + %12.8fi\n",
           std::real(G[0][sc]),
           std::imag(G[0][sc]),
           std::real(G[1][sc]),
           std::imag(G[1][sc])
           );
  }
  return NULL;
}

// structure to hold command line inputs
typedef struct
{
  double cent_freq;         // center frequency of transmission
  double samp_rate;         // usrp samping rate
  float dsp_gain;           // dsp gain
  double txgain;            // tx frontend gain
  double rxgain;            // rx frontend gain
  unsigned int M;           // number of subcarriers
  unsigned int cp_len;      // length of cyclic prefix
  bool verbose;             // verbosity of the app
} options;

// structure to hold usrp configurations
typedef struct
{
  double cent_freq;               // center frequency
  double samp_rate;               // sampling ratea
  float rf_gain;                  // sampling ratea
  clock_source_type clock_source; // clock source TODO List the sources
  time_source_type time_source;   // time source TODO List the sources 
} usrp_config;

// structure to load data to simulator
typedef struct
{
} simulator_data;

// structure to load data to tx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * tx;
  framegen * fg;
} tx_thread_data;

// structure to load data to rx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * rx;
  framesync * fs;
} rx_thread_data;

typedef struct
{
  uhd::usrp::multi_usrp::sptr * rx;
  liquid::ofdm::demodulator * dem;
} ofdm_dem_data;

void * ofdm_dem (void * _data)
{
  // vector of pointers to the modulator buffers.
  // TODO modifications for multichannel modulators
  std::complex<float> ** demodulator_buffer;
  // usrp buffer
  std::vector<std::complex<float> *> usrp_buffer;
  // usrp buffer len
  size_t usrp_buffer_len;
  // number of channels available on USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples received in each call of recv
  unsigned int num_samples_rcvd;
  // timeout
  float timeout = 0.2;
  // file to save data;
  std::vector<FILE *> fp;
  // total number of samples received
  unsigned int num_accumulated_samples;

  ofdm_dem_data * data = (ofdm_dem_data *)_data;

  num_channels = (*(data->rx))->get_rx_num_channels();

  for(size_t chan = 0; chan < num_channels; chan++)
    channels.push_back(chan);
  uhd::stream_args_t rx_stream_args(CPU, WIRE);
  rx_stream_args.channels = channels;
  uhd::rx_streamer::sptr rx_stream = 
    (*(data->rx))->get_rx_stream(rx_stream_args);

  usrp_buffer_len = rx_stream->get_max_num_samps();
  demodulator_buffer = (std::complex<float> **) 
    malloc(num_channels*sizeof(std::complex<float> *));
  for(size_t chan = 0; chan < num_channels; chan++) {
    usrp_buffer.push_back((std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*usrp_buffer_len));
    demodulator_buffer[chan] = usrp_buffer[chan];
    if(LOG) {
      fp.push_back(fopen(boost::str(
        boost::format("%srx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
    }
  }

  uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_cmd.stream_now = false;
  stream_cmd.time_spec = (*(data->rx))->get_time_now() 
                       + uhd::time_spec_t(0.1);
  rx_stream->issue_stream_cmd(stream_cmd);
  rx_begin = time(NULL);
  // lock mutex and signal tx to start
  std::cout << "Rx locking the mutex and signalling tx\n";
  pthread_mutex_lock(&mutex_txrx);
  start_tx = true;
  pthread_cond_signal(&condition_txrx);
  pthread_mutex_unlock(&mutex_txrx);

  num_accumulated_samples = 0;

  bool break_loop = false;
  while(1)
  {
    num_samples_rcvd = rx_stream->recv(usrp_buffer,
                                       usrp_buffer_len,
                                       rxmd,
                                       timeout);
    if(rxmd.error_code) {
      break;
    }
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite((void *)demodulator_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_rcvd,
                      fp[chan]) == num_samples_rcvd);
      }
    }
    (data->dem)->demodulate_samples(demodulator_buffer[0],
                                    num_samples_rcvd);
    num_accumulated_samples += num_samples_rcvd;

    // check of stop_rx is set
    pthread_mutex_lock(&mutex_txrx);
    if(stop_rx) {
      std::cout << "stop_rx is set\n";
      break_loop = true;
    }
    pthread_mutex_unlock(&mutex_txrx);
    if(break_loop)
      break;
  }

  rx_end = time(NULL);
  std::cout << "Exiting rx thread\n";
  // stop rx streaming
  stream_cmd.stream_mode = 
    uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
  rx_stream->issue_stream_cmd(stream_cmd);


  for(size_t chan = 0; chan < num_channels; chan++) {
    free(usrp_buffer[chan]);
  }
  free(demodulator_buffer);
  // close files
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++)
      fclose(fp[chan]);
  }
  pthread_exit(NULL);
}

typedef struct
{
  uhd::usrp::multi_usrp::sptr * tx;
  liquid::ofdm::modulator * mod;
} ofdm_mod_data;

void * ofdm_mod (void * _data)
{
  // number of subcarriers
  unsigned int M;
  // cyclic prefix length
  unsigned int cp_len;
  // number of OFDM symbols per packet
  unsigned int frame_len;
  // number of samples per frame at
  // the output of modulator
  unsigned int modulator_buffer_len;
  // vector of pointers to the modulator buffers.
  // TODO modifications for multichannel modulators
  std::complex<float> ** modulator_buffer;
  // usrp buffer
  std::vector<std::complex<float> *> usrp_buffer;
  // usrp buffer
  std::vector<std::complex<float> *> zero_buffer;
  // payload len
  unsigned int payload_len;
  // header len
  unsigned int header_len;
  // payload data
  unsigned char * payload;
  // header data
  unsigned char * header;
  // number of channels available with USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples sent in each call of send
  unsigned int num_samples_sent;
  std::vector<FILE *> fp;
  
  ofdm_mod_data * data = (ofdm_mod_data *)_data;

  M = (data->mod)->get_M();
  cp_len = (data->mod)->get_cp_len();
  frame_len = (data->mod)->get_frame_len();
  modulator_buffer_len = (M + cp_len)*frame_len;
  header_len = (data->mod)->get_h_usr_len();
  payload_len = (data->mod)->get_payload_dec_len();
  num_channels = (*(data->tx))->get_tx_num_channels();

  payload = (unsigned char *) malloc 
    (payload_len*sizeof(unsigned char));
  header = (unsigned char *) malloc 
    (header_len*sizeof(unsigned char));
  modulator_buffer = (std::complex<float> **) 
    malloc(num_channels*sizeof(std::complex<float> *));
  std::complex<float> Z(0.0, 0.0);
  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    modulator_buffer[chan] = (std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*modulator_buffer_len);
    usrp_buffer.push_back(modulator_buffer[chan]);
    zero_buffer.push_back((std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*modulator_buffer_len));
    std::fill(zero_buffer[chan], zero_buffer[chan] + 
              modulator_buffer_len,
              Z);
    if(LOG) {
      fp.push_back(fopen(boost::str(
        boost::format("%stx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
    }
  }

  uhd::stream_args_t tx_stream_args(CPU, WIRE);
  tx_stream_args.channels = channels;
  uhd::tx_streamer::sptr tx_stream = 
    (*(data->tx))->get_tx_stream(tx_stream_args);
  /*
   * Lock mutex and wait for signal.  
   * Note that the pthread_cond_wait routine will 
   * automatically and atomically unlock mutex while it waits.
   * Also, note that if start_tx is true before this 
   * routine is run by the tx thread, the loop will be 
   * skipped to prevent pthread_cond_wait from never returning. 
   */
  pthread_mutex_lock(&mutex_txrx);
    while (!(start_tx)) {
      // blocks tx thread
      std::cout << "Waiting for signal from rx\n";
      pthread_cond_wait(&condition_txrx, &mutex_txrx);
      std::cout << "Sigal received from rx thread\n";
    }
  pthread_mutex_unlock(&mutex_txrx);

  uhd::tx_metadata_t txmd;
  txmd.start_of_burst = true;
  txmd.end_of_burst = false;
  txmd.has_time_spec = true;
  txmd.time_spec = (*(data->tx))->get_time_now()
                 + uhd::time_spec_t(0.1);
  tx_begin = time(NULL);
  num_samples_sent = tx_stream->send(zero_buffer,
                                     modulator_buffer_len,
                                     txmd);
  assert(num_samples_sent == modulator_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(zero_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  txmd.start_of_burst = false;
  txmd.has_time_spec = false;

  unsigned int i;
  pid = 0;
  while(pid < PID_MAX)
  {
    header[0] = (pid >> 8) & 0xff;
    header[1] = (pid     ) & 0xff;
    for(i = 2; i < header_len; i++)
      header[i] = rand() & 0xff;

    for(i = 0; i < payload_len; i++)
      payload[i] = rand() & 0xff;

    (data->mod)->assemble_frame(header, payload);
    (data->mod)->assemble_output_samples(modulator_buffer[0]);
    num_samples_sent = tx_stream->send(usrp_buffer,
                                       modulator_buffer_len,
                                       txmd);
    assert(num_samples_sent == modulator_buffer_len);
    // logging
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite(usrp_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_sent,
                      fp[chan]) == num_samples_sent);
      }
    }
    pid++;
    // send some zeros
    num_samples_sent = tx_stream->send(zero_buffer,
                                       modulator_buffer_len,
                                       txmd);
    assert(num_samples_sent == modulator_buffer_len);
    // logging
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite(zero_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_sent,
                      fp[chan]) == num_samples_sent);
      }
    }
  }

  // transmission ends
  tx_end = time(NULL);
  usleep(1000000);
  pthread_mutex_lock(&mutex_txrx);
  stop_rx = true;
  pthread_mutex_unlock(&mutex_txrx);

  std::cout << "Exiting tx thread\n";
  // send the last packet
  txmd.end_of_burst = true;
  tx_stream->send(usrp_buffer,
                  0,
                  txmd);

  for(size_t chan = 0; chan < num_channels; chan++) {
    free(modulator_buffer[chan]);
    if(LOG) {
      fclose(fp[chan]);
    }
  }
  free(modulator_buffer);

  pthread_exit(NULL);
}

// function to read commandline options
// argc: argument count
// argv: argument list
// parameter options: pointer to a struct_options object
void read_options(int       argc,
                  char  **  argv,
                  options * d_options)
{
  //set the operating parameters
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "help message")
    ("freq,f",
     po::value<double>(&(d_options->cent_freq))->default_value(
                                            CENTER_FREQUENCY),
     "RF center frequency in Hz")
    ("rate,r",
     po::value<double>(&(d_options->samp_rate))->default_value(
                                            SAMPLING_RATE),
     "USRP Sampling rate")
    ("dsp_gain",
     po::value<float>(&(d_options->dsp_gain))->default_value(
                                          BASEBAND_GAIN),
     "TX DSP gain")
    ("tx_gain",
     po::value<double>(&(d_options->txgain))->default_value(
                                         TX_FRONTEND_GAIN),
     "TX Front end gain")
    ("rx_gain",
     po::value<double>(&(d_options->rxgain))->default_value(
                                         RX_FRONTEND_GAIN),
     "RX Front end gain")
    ("num_subcarriers",
     po::value<unsigned int>(&(d_options->M))->default_value(
                                          NUM_SUBCARRIERS),
     "Number of OFDM subcarriers")
    ("cp_len",
     po::value<unsigned int>(&(d_options->cp_len))->default_value(
                                               CP_LENGTH),
     "Cyclic Prefix Length")
    ("verbose,v",
     "Verbose")
    ("quite,q",
     "Quite")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  //print the help message
  if (vm.count("help")) {
    std::cout << boost::format("ofdmtxrx %s") % desc << std::endl;
    std::cout
      << std::endl
      << "Basic OFDM P2P Link\n"
      << std::endl;
    exit(0);
  }

  //check sanity of options
  if (vm.count("verbose") && vm.count("quite")) {
    std::cout << "Conflicting Options Verbose and Quite."
      << " Please use only one of those."
      << std::endl;
    exit(0);
  }

  if (vm.count("verbose"))
    d_options->verbose = true;
}

// initialize USRP
// parameter u: sptr to USRP
// parameter configs: values for freq, gain etc
// parameter is_tx: set to true if USRP is transmitting
// TODO set device time to 0 before tuning
void init_usrp(uhd::usrp::multi_usrp::sptr u,
               usrp_config * configs,
               bool is_tx)
{
  size_t num_chans;     // number of channels in USRP
  // set clock source
  switch(configs->clock_source) {
    case(CLOCK_SOURCE_NONE):
      break;
    case(CLOCK_SOURCE_EXTERNAL):
      u->set_clock_source("external");
      break;
    default:
      std::cout << "Clock source not known. Exiting\n";
      exit(1);
  }
  // set time source
  switch(configs->time_source) {
    case(TIME_SOURCE_NONE):
      break;
    case(TIME_SOURCE_EXTERNAL):
      u->set_time_source("external");
      break;
    default:
      std::cout << "Time source not known. Exiting\n";
      exit(1);
  }

  if(is_tx) {
    // set subdev specs
    if(u->get_mboard_name() == "X300") {
      u->set_tx_subdev_spec(
          uhd::usrp::subdev_spec_t(X300_SUBDEV_SPEC),
          uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if((u->get_mboard_name() == "N200") || 
            (u->get_mboard_name() == "N200r4")){
      u->set_tx_subdev_spec(
          uhd::usrp::subdev_spec_t(N200_SUBDEV_SPEC),
          uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if(u->get_mboard_name() == "B210") {
      u->set_tx_subdev_spec(
          uhd::usrp::subdev_spec_t(B210_SUBDEV_SPEC_TX),
          uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else {
      std::cout << "TX Motherboard not compatible\n"
                << "Subdevice specification for "
                << u->get_mboard_name()
                << " not known. Exiting\n";
      exit(1);
    }
    num_chans = u->get_tx_num_channels();
    // set freq, gain and antenna
    // TODO pass antenna as a parameter
    for (size_t chan = 0; chan < num_chans; chan++) {
      u->set_tx_rate(configs->samp_rate, chan);
      u->set_tx_gain(configs->rf_gain, chan);
      u->set_tx_antenna("TX/RX", chan);
    }
    // tune the channels simultaneously to obtain
    // same phase on all the channels
    uhd::time_spec_t cmd_time = u->get_time_now() + 
                                uhd::time_spec_t(0.1);
    u->set_command_time(cmd_time);
    uhd::tune_request_t tx_tune_request(configs->cent_freq);
    for (size_t chan = 0; chan < num_chans; chan++) 
      u->set_tx_freq(tx_tune_request, chan);
    u->clear_command_time();
  }
  else {
    // set subdev specs
    if(u->get_mboard_name() == "X300") {
      u->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
            X300_SUBDEV_SPEC),
            uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if((u->get_mboard_name() == "N200") ||
            (u->get_mboard_name() == "N200r4")) {
      u->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
            N200_SUBDEV_SPEC),
            uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else if(u->get_mboard_name() == "B210") {
      u->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
            B210_SUBDEV_SPEC_RX),
            uhd::usrp::multi_usrp::ALL_MBOARDS);
    }
    else {
      std::cout << "RX Motherboard not compatible\n"
                << "Subdevice specification for "
                << u->get_mboard_name()
                << " not known. Exiting\n";
      exit(1);
    }
    num_chans = u->get_rx_num_channels();
    // set freq, gain and antenna
    // TODO pass antenna as a parameter
    for (size_t chan = 0; chan < num_chans; chan++) {
      u->set_rx_rate(configs->samp_rate, chan);
      u->set_rx_gain(configs->rf_gain, chan);
      u->set_rx_antenna("TX/RX", chan);
    }
    // tune the channels simultaneously to obtain
    // same phase on all the channels
    uhd::time_spec_t cmd_time = u->get_time_now() + 
                                uhd::time_spec_t(0.1);
    u->set_command_time(cmd_time);
    uhd::tune_request_t rx_tune_request(configs->cent_freq);
    for (size_t chan = 0; chan < num_chans; chan++) 
      u->set_rx_freq(rx_tune_request, chan);
    u->clear_command_time();
  }
}

// simulator program
void * simulator (void * _data)
{
  return NULL;
}

// rx_worker thread
void * rx_worker (void * _data)
{
  rx_thread_data * data = (rx_thread_data *)_data;
  // rx buffer length
  unsigned int rx_buffer_len;
  // number of channels available on USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples read in each call of recv
  unsigned int num_samples_read;
  // usrp buffer
  std::vector<std::complex<float> *> rx_buffer;
  // log files
  std::vector<FILE *> fp;
  // state returned from framesync
  framesync_states_t sync_state = STATE_SEEK_PLATEAU;

  // FIXME is num_channels = NUM_STREAMS always??
  num_channels = (*(data->rx))->get_rx_num_channels();

  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    if(LOG) {
      fp.push_back(fopen(boost::str(
        boost::format("%srx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
    }
  }

  // initialize rx streamer
  uhd::stream_args_t rx_stream_args(CPU, WIRE);
  rx_stream_args.channels = channels;
  uhd::rx_streamer::sptr rx_stream = 
    (*(data->rx))->get_rx_stream(rx_stream_args);

  // initialilze rx_buffer_len, allocate memory
  rx_buffer_len = rx_stream->get_max_num_samps();
  for(size_t chan = 0; chan < num_channels; chan++) {
    // FIXME is volk malloc required here??
    rx_buffer.push_back((std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*rx_buffer_len));
  }

  // start receiving in 0.1 seconds
  uhd::stream_cmd_t stream_cmd(
      uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_cmd.stream_now = false;
  stream_cmd.time_spec = (*(data->rx))->get_time_now() 
                       + uhd::time_spec_t(0.1);
  rx_stream->issue_stream_cmd(stream_cmd);
  // lock mutex and signal tx to start
  std::cout << "Rx locking the mutex and signalling tx\n";
  pthread_mutex_lock(&mutex_txrx);
  start_tx = true;
  pthread_cond_signal(&condition_txrx);
  pthread_mutex_unlock(&mutex_txrx);

  rx_begin = time(NULL);
  // NOTE: timout should be larger if rx_buffer_length
  // is larger.
  float timeout = 0.2;
  // FIXME on what condition does rx stop streaming?
  // TODO receive sigal from tx to stop
  bool break_loop = false;
  while(1)
  {
    num_samples_read = rx_stream->recv(rx_buffer,
                                       rx_buffer_len,
                                       rxmd,
                                       timeout);
    if(rxmd.error_code)
      break; // TODO error code handling to be done in main

    // logging
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite(rx_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_read,
                      fp[chan]) == num_samples_read);
      }
    }

    // process the received samples through framesync
    if(sync_state < STATE_ZF) {
      sync_state = (data->fs)->execute(rx_buffer, num_samples_read);
      if(sync_state >= STATE_ZF) {
        // lock mutex and signal tx to start zf
        std::cout << "Rx signalling tx to start zf\n";
        pthread_mutex_lock(&mutex_txrx);
        start_zf = true;
        pthread_cond_signal(&condition_zf);
        pthread_mutex_unlock(&mutex_txrx);
      // TODO signal tx to start zero-forcing.
      }
    }

    // check of stop_rx is set
    pthread_mutex_lock(&mutex_txrx);
    if(stop_rx) {
      std::cout << "stop_rx is set\n";
      break_loop = true;
    }
    pthread_mutex_unlock(&mutex_txrx);
    if(break_loop)
      break;
  }

  rx_end = time(NULL);
  // stop rx streaming
  stream_cmd.stream_mode = 
    uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
  rx_stream->issue_stream_cmd(stream_cmd);

  // free memory
  for(size_t chan = 0; chan < num_channels; chan++) {
    free(rx_buffer[chan]);
  }
  // close files
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++)
      fclose(fp[chan]);
  }
  std::cout << "Exiting rx thread\n";
  pthread_exit(NULL);
}

// tx_worker thread
void * tx_worker (void * _data)
{
  tx_thread_data * data = (tx_thread_data *)_data;
  // tx_buffer_len = length of sync word
  unsigned int tx_buffer_len = 
    (CP_LENGTH + NUM_SUBCARRIERS)*
    (NUM_STREAMS*NUM_ACCESS_CODES + 1);
  // number of channels available with USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples sent in each call of send
  unsigned int num_samples_sent;
  // number of samples read from framegen
  unsigned int num_samples_read;
  // usrp buffer
  std::vector<std::complex<float> *> tx_buffer;
  // framegen buffer
  std::vector<std::complex<float> *> fg_buffer;
  // tx multiplier gain
  std::vector<std::complex<float> > bb_gain;
  // log files
  std::vector<FILE *> fp;

  // FIXME is num_channels = NUM_STREAMS always??
  num_channels = (*(data->tx))->get_tx_num_channels();

  size_t volk_alignment = volk_get_alignment();
  std::complex<float> Z(0.0, 0.0);
  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    tx_buffer.push_back((std::complex<float> *) volk_malloc 
      (sizeof(std::complex<float>)*tx_buffer_len,
       volk_alignment));
    fg_buffer.push_back((std::complex<float> *) volk_malloc 
      (sizeof(std::complex<float>)*tx_buffer_len,
       volk_alignment));
    std::fill(tx_buffer[chan], tx_buffer[chan] + tx_buffer_len,
              Z);
    if(LOG) {
      fp.push_back(fopen(boost::str(
        boost::format("%stx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
    }
    bb_gain.push_back(Z);
  }

  // reset baseband gains
  bb_gain[0] = BASEBAND_GAIN;
  bb_gain[1] = BASEBAND_GAIN;

  // initialize tx streamer
  uhd::stream_args_t tx_stream_args(CPU, WIRE);
  tx_stream_args.channels = channels;
  uhd::tx_streamer::sptr tx_stream = 
    (*(data->tx))->get_tx_stream(tx_stream_args);

  /*
   * Lock mutex and wait for signal.  
   * Note that the pthread_cond_wait routine will 
   * automatically and atomically unlock mutex while it waits.
   * Also, note that if start_tx is true before this 
   * routine is run by the tx thread, the loop will be 
   * skipped to prevent pthread_cond_wait from never returning. 
   */
  pthread_mutex_lock(&mutex_txrx);
  while (!(start_tx)) {
    // blocks tx thread
    std::cout << "Waiting for signal from rx\n";
    pthread_cond_wait(&condition_txrx, &mutex_txrx);
    std::cout << "Sigal received from rx thread\n";
  }
  pthread_mutex_unlock(&mutex_txrx);

  // begin tx
  uhd::tx_metadata_t txmd;
  txmd.start_of_burst = true;
  txmd.end_of_burst = false;
  txmd.has_time_spec = true;
  txmd.time_spec = (*(data->tx))->get_time_now()
                 + uhd::time_spec_t(0.1);

  tx_begin = time(NULL);

  // transmit zeros to flush out the buffers
  num_samples_sent = tx_stream->send(tx_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }

  num_samples_read = (data->fg)->write_sync_words(fg_buffer);
  assert(num_samples_read <= tx_buffer_len);
  // apply baseband gain
  // FIXME see if this multiply const can be done in place.
  // see if it is benificial to do this in separate
  // thread for each channel.
  for(size_t chan = 0; chan < num_channels; chan++) {
    volk_32fc_s32fc_multiply_32fc(tx_buffer[chan],
                                  fg_buffer[chan],
                                  bb_gain[chan],
                                  num_samples_read);
    // fill fg_buffer with zeros
    std::fill(fg_buffer[chan], fg_buffer[chan] + tx_buffer_len,
              Z);
  }
  // set tx metadata
  txmd.start_of_burst = false;
  txmd.has_time_spec = false;

  // send sync words
  for(unsigned int repeat = 0; repeat < 1; repeat++) {
    num_samples_sent = tx_stream->send(tx_buffer,
                                       num_samples_read,
                                       txmd,
                                       1.0);
    assert(num_samples_sent == num_samples_read);
    // logging
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite(tx_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_sent,
                      fp[chan]) == num_samples_sent);
      }
    }
    // send zeros in between
    for(unsigned int zBlock = 0; zBlock < 1; zBlock++) {
      num_samples_sent = tx_stream->send(fg_buffer,
                                         num_samples_read,
                                         txmd,
                                         1.0);
      assert(num_samples_sent == num_samples_read);
      // logging
      if(LOG) {
        for(size_t chan = 0; chan < num_channels; chan++) {
          assert(fwrite(fg_buffer[chan],
                        sizeof(std::complex<float>),
                        num_samples_sent,
                        fp[chan]) == num_samples_sent);
        }
      }
    }
  }
  // stop streaming
  txmd.end_of_burst = true;
  num_samples_sent = tx_stream->send(fg_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(fg_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  
  pthread_mutex_lock(&mutex_txrx);
  while (!(start_zf)) {
    // blocks tx thread
    std::cout << "Waiting for signal from rx to begin zf\n";
    pthread_cond_wait(&condition_zf, &mutex_txrx);
    std::cout << "Begin zf sigal received from rx thread\n";
  }
  pthread_mutex_unlock(&mutex_txrx);

  // update the zero forcing weights
  (data->fg)->compute_W(G);
  //************* ZERO FORCING ******************//
  std::cout << "********** begin zf ***********\n";
  txmd.start_of_burst = true;
  txmd.end_of_burst = false;
  txmd.has_time_spec = true;
  txmd.time_spec = (*(data->tx))->get_time_now()
                 + uhd::time_spec_t(0.001);
  num_samples_sent = tx_stream->send(fg_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(fg_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  txmd.start_of_burst = false;
  txmd.has_time_spec = false;

  for(unsigned int repeat = 0; repeat < 2; repeat++) {
    for(unsigned int words = 0; words < 40; words++) {
      num_samples_read = (data->fg)->write_zf_words_random(fg_buffer);
      assert(num_samples_read <= tx_buffer_len);
      // apply baseband gain
      // FIXME see if this multiply const can be done in place.
      // see if it is benificial to do this in separate
      // thread for each channel.
      for(size_t chan = 0; chan < num_channels; chan++) {
        volk_32fc_s32fc_multiply_32fc(tx_buffer[chan],
                                      fg_buffer[chan],
                                      bb_gain[chan],
                                      num_samples_read);
      }
      num_samples_sent = tx_stream->send(tx_buffer,
                                         num_samples_read,
                                         txmd,
                                         1.0);
      assert(num_samples_sent == num_samples_read);
      // logging
      if(LOG) {
        for(size_t chan = 0; chan < num_channels; chan++) {
          assert(fwrite(tx_buffer[chan],
                        sizeof(std::complex<float>),
                        num_samples_sent,
                        fp[chan]) == num_samples_sent);
        }
      }
    }
    for(unsigned int words = 0; words < 40; words++) {
      num_samples_read = (data->fg)->write_words_random(fg_buffer);
      assert(num_samples_read <= tx_buffer_len);
      // apply baseband gain
      // FIXME see if this multiply const can be done in place.
      // see if it is benificial to do this in separate
      // thread for each channel.
      for(size_t chan = 0; chan < num_channels; chan++) {
        volk_32fc_s32fc_multiply_32fc(tx_buffer[chan],
                                      fg_buffer[chan],
                                      bb_gain[chan],
                                      num_samples_read);
      }
      num_samples_sent = tx_stream->send(tx_buffer,
                                         num_samples_read,
                                         txmd,
                                         1.0);
      assert(num_samples_sent == num_samples_read);
      // logging
      if(LOG) {
        for(size_t chan = 0; chan < num_channels; chan++) {
          assert(fwrite(tx_buffer[chan],
                        sizeof(std::complex<float>),
                        num_samples_sent,
                        fp[chan]) == num_samples_sent);
        }
      }
    }
    for(size_t chan = 0; chan < num_channels; chan++)
      std::fill(tx_buffer[chan], tx_buffer[chan] + tx_buffer_len, Z);
    for(unsigned int words = 0; words < 10; words++) {
      num_samples_sent = tx_stream->send(tx_buffer,
                                         tx_buffer_len,
                                         txmd,
                                         1.0);
      assert(num_samples_sent == tx_buffer_len);
      // logging
      if(LOG) {
        for(size_t chan = 0; chan < num_channels; chan++) {
          assert(fwrite(tx_buffer[chan],
                        sizeof(std::complex<float>),
                        num_samples_sent,
                        fp[chan]) == num_samples_sent);
        }
      }
    }
  }
  // stop streaming
  txmd.end_of_burst = true;
  num_samples_sent = tx_stream->send(tx_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  
  // sleep for some time and signal rx to stop
  tx_end = time(NULL);
  usleep(1000000);
  pthread_mutex_lock(&mutex_txrx);
  stop_rx = true;
  pthread_mutex_unlock(&mutex_txrx);

  // free memory
  for(size_t chan = 0; chan < num_channels; chan++) {
    volk_free(tx_buffer[chan]);
    volk_free(fg_buffer[chan]);
  }
  // close files
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++)
      fclose(fp[chan]);
  }
  // return
  std::cout << "Exiting tx thread\n";
  pthread_exit(NULL);
}

// main
int UHD_SAFE_MAIN(int argc, char **argv)
{
  uhd::set_thread_priority_safe();

  // read command line options
  options d_options;
  read_options(argc, argv, &d_options);

  unsigned int M = NUM_SUBCARRIERS;
  unsigned int cp_len = CP_LENGTH;
  unsigned int num_streams = NUM_STREAMS;

  // allocate space for CSI and W
  unsigned int
    volk_alignment = volk_get_alignment();
  std::complex<float> I(1.0, 0.0);
  G = (std::complex<float> **) malloc
      (sizeof(std::complex<float> *)*num_streams);
  W = (std::complex<float> **) malloc
      (sizeof(std::complex<float> *)*num_streams);
  for(unsigned int i = 0; i < num_streams; i++)
  {
    G[i] = (std::complex<float> *) volk_malloc
      (sizeof(std::complex<float>)*M,
       volk_alignment);
    W[i] = (std::complex<float> *) volk_malloc
      (sizeof(std::complex<float>)*M,
       volk_alignment);
    std::fill(W[i], W[i] + M, I);
  }

  std::vector<unsigned char *> p;
  std::vector<msequence> ms_S0;
  std::vector<msequence> ms_S1;
  std::vector<msequence> ms_pilot;
  p.resize(num_streams);
  ms_S0.resize(num_streams);
  ms_S1.resize(num_streams);
  ms_pilot.resize(num_streams);
  for(unsigned int i = 0; i < num_streams; i++) {
    p[i] = (unsigned char *) malloc 
           (sizeof(unsigned char)*M);
    ofdmframe_init_default_sctype(p[i], M);
  }
  // NOTE: initialize the sequence generator for other streams as well
  // in case the number of streams is more that 2. Here we assume that
  // we only have 2 streams.
  ms_S0[0] = msequence_create(LFSR_SMALL_LENGTH, LFSR_SMALL_0_GEN_POLY, 1);
  ms_S0[1] = msequence_create(LFSR_SMALL_LENGTH, LFSR_SMALL_1_GEN_POLY, 1);
  ms_S1[0] = msequence_create(LFSR_LARGE_LENGTH, LFSR_LARGE_0_GEN_POLY, 1);
  ms_S1[1] = msequence_create(LFSR_LARGE_LENGTH, LFSR_LARGE_1_GEN_POLY, 1);

  framegen fg(M,
              cp_len,
              num_streams,
              NUM_ACCESS_CODES,
              p,
              ms_S0,
              ms_S1,
              ms_pilot);
  fg.set_W(W);

  msequence_reset(ms_S0[0]);
  msequence_reset(ms_S0[1]);
  msequence_reset(ms_S1[0]);
  msequence_reset(ms_S1[1]);
  
  framesync fs(M,
               cp_len,
               num_streams,
               NUM_ACCESS_CODES,
               p,
               ms_S0,
               ms_S1,
               ms_pilot,
               callback);
  // destroy msequence ojbects.
  msequence_destroy(ms_S0[0]);
  msequence_destroy(ms_S0[1]);
  msequence_destroy(ms_S1[0]);
  msequence_destroy(ms_S1[1]);

  usrp_config * tx_config;
  usrp_config * rx_config;
  tx_thread_data * tx_data;
  rx_thread_data * rx_data;
  ofdm_mod_data * mod_data;
  ofdm_dem_data * dem_data;
  uhd::usrp::multi_usrp::sptr tx;
  uhd::usrp::multi_usrp::sptr rx;

  // for ofdm tx rx loops
  void * user_data = (void *)&(d_options.samp_rate);
  std::complex<float> ofdm_mod_gain(BASEBAND_GAIN, 0.0);
  OFDMFRAME_STRUCT * frame_struct;
  frame_struct = (OFDMFRAME_STRUCT *) malloc 
    (sizeof(OFDMFRAME_STRUCT));
  frame_struct->VERSN =  104;
  frame_struct->H_USR =  2;
  frame_struct->H_LEN =  frame_struct->H_USR + 6;
  frame_struct->H_EC1 =  LIQUID_FEC_GOLAY2412;
  frame_struct->H_EC2 =  LIQUID_FEC_NONE;
  frame_struct->H_CRC =  LIQUID_CRC_32;
  frame_struct->H_MOD =  LIQUID_MODEM_BPSK;
  frame_struct->P_LEN =  512;
  frame_struct->P_EC1 =  LIQUID_FEC_NONE;
  frame_struct->P_EC2 =  LIQUID_FEC_NONE;
  frame_struct->P_CRC =  LIQUID_CRC_NONE;
  frame_struct->P_MOD =  LIQUID_MODEM_QPSK;

  liquid::ofdm::modulator mod(d_options.M,
                              d_options.cp_len,
                              TAPER_LEN,
                              NULL,
                              frame_struct);
  mod.set_tx_gain(ofdm_mod_gain);
  liquid::ofdm::demodulator dem(d_options.M,
                                d_options.cp_len,
                                TAPER_LEN,
                                NULL,
                                ofdm_callback,
                                user_data,
                                frame_struct);

  num_frames_detected         = 0;
  num_valid_headers_received  = 0;
  num_valid_bytes_received    = 0;
  num_valid_packets_received  = 0;

  pthread_t tx_thread, rx_thread;

  if(SIMULATION) {
  } 
  else {
    // initialize USRP
    rx = uhd::usrp::multi_usrp::make(
         uhd::device_addr_t(B210B));
    tx = uhd::usrp::multi_usrp::make(
         uhd::device_addr_t(B210A));
    tx_config = (usrp_config *) malloc 
                (sizeof(usrp_config));
    rx_config = (usrp_config *) malloc 
                (sizeof(usrp_config));
    // assign tx and rx sampling rate freq etc
    tx_config->cent_freq = d_options.cent_freq;
    rx_config->cent_freq = d_options.cent_freq;
    tx_config->samp_rate = d_options.samp_rate;
    rx_config->samp_rate = d_options.samp_rate;
    tx_config->rf_gain   = d_options.txgain;
    rx_config->rf_gain   = d_options.rxgain;
    // TODO Add command line options 
    // for clock and time source
    tx_config->clock_source = CLOCK_SOURCE;
    rx_config->clock_source = CLOCK_SOURCE;
    tx_config->time_source = TIME_SOURCE;
    rx_config->time_source = TIME_SOURCE;

    init_usrp(tx, tx_config, true);
    init_usrp(rx, rx_config, false);

    // initialize worker threads
    tx_data = (tx_thread_data *) malloc (sizeof(tx_thread_data));
    rx_data = (rx_thread_data *) malloc (sizeof(rx_thread_data));
    tx_data->tx = &tx;
    rx_data->rx = &rx;
    tx_data->fg = &fg;
    rx_data->fs = &fs;

    // ofdm modem
    mod_data = (ofdm_mod_data *) malloc (sizeof(tx_thread_data));
    dem_data = (ofdm_dem_data *) malloc (sizeof(rx_thread_data));
    mod_data->tx = &tx;
    dem_data->rx = &rx;
    mod_data->mod = &mod;
    dem_data->dem = &dem;

//    {
//      if(pthread_create(&tx_thread, NULL, ofdm_mod, (void *)mod_data)){
//        std::cout << "Error invoking tx thread\n";
//        return 1;
//      }
//      if(pthread_create(&rx_thread, NULL, ofdm_dem, (void *)dem_data)){
//        std::cout << "Error invoking rx thread\n";
//        return 1;
//      }
//      pthread_join(tx_thread, NULL);
//      pthread_join(rx_thread, NULL);
//    }
    {
      if(pthread_create(&tx_thread, NULL, tx_worker, (void *)tx_data)){
        std::cout << "Error invoking tx thread\n";
        return 1;
      }
      if(pthread_create(&rx_thread, NULL, rx_worker, (void *)rx_data)){
        std::cout << "Error invoking rx thread\n";
        return 1;
      }
      pthread_join(tx_thread, NULL);
      pthread_join(rx_thread, NULL);
    }
//    {
//      if(pthread_create(&tx_thread, NULL, tx_worker, (void *)tx_data)){
//        std::cout << "Error invoking tx thread\n";
//        return 1;
//      }
//      if(pthread_create(&rx_thread, NULL, ofdm_dem, (void *)dem_data)){
//        std::cout << "Error invoking rx thread\n";
//        return 1;
//      }
//      pthread_join(tx_thread, NULL);
//      pthread_join(rx_thread, NULL);
//    }
//    {
//      if(pthread_create(&tx_thread, NULL, ofdm_mod, (void *)mod_data)){
//        std::cout << "Error invoking tx thread\n";
//        return 1;
//      }
//      if(pthread_create(&rx_thread, NULL, rx_worker, (void *)rx_data)){
//        std::cout << "Error invoking rx thread\n";
//        return 1;
//      }
//      pthread_join(tx_thread, NULL);
//      pthread_join(rx_thread, NULL);
//    }

    // free malloc
    free(tx_data);
    free(rx_data);
    free(tx_config);
    free(rx_config);
    free(mod_data);
    free(dem_data);
  }
  for(unsigned int i = 0; i < num_streams; i++) {
    volk_free(G[i]);
    volk_free(W[i]);
    free(p[i]);
  }
  free(G);
  free(W);
  free(frame_struct);
  // print experiment output
  printf("    plateau width           : %6lu\n", fs.get_plateau_end() - 
                                                 fs.get_plateau_start() + 1);
  printf("    plateau start           : %6lu\n", fs.get_plateau_start());
  printf("    plateau end             : %6lu\n", fs.get_plateau_end());
  printf("    frames sync index       : %6lu\n", fs.get_sync_index());
  printf("    num samples processed   : %6llu\n", fs.get_num_samples_processed());
  printf("    frames transmitted      : %6u\n", pid);
  printf("    frames detected         : %6u\n", num_frames_detected);
  printf("    valid headers           : %6u\n", num_valid_headers_received);
  printf("    valid packets           : %6u\n", num_valid_packets_received);
  printf("    bytes received          : %6u\n", num_valid_bytes_received);
  printf("    tx run time             : %6lu\n", tx_end - tx_begin);
  printf("    rx run time             : %6lu\n", rx_end - rx_begin);
  printf("    bit rate                : %e bps\n",
         float(num_valid_bytes_received*8)/(rx_end - rx_begin));
  if(!(DECODE_ONLINE))
    printf("    decoding time           : %3.3f sec\n",
           time_for_offline_decoding);
  if(rxmd.error_code)
  {
    std::cout << "    rx thread exit with     : "
              << rxmd.strerror() << "\n";
  }
  return 0;
}
