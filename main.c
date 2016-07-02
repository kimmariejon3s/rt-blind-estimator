/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <sndfile.h>
#include <libswresample/swresample.h>
#include <libavutil/channel_layout.h>
#include <libavcodec/avcodec.h>
#include <liquid/liquid.h>
#include <nlopt.h>
#include "kiss_fft.h"
#include <svdlib.h>

//FIXME: make this 20 when I decide to do the split stuff
#define SPLIT		1
#define MAX_CHUNK_SIZE	1024 * 1024 * 1024
#define WINDOW_WIDTH	0.5
#define OVERLAP		0.98
#define OCT_FILT_ORDER	5
#define ENV_FILT_ORDER	4
#define N_OCT 		OCT_FILT_ORDER * 2
#define R_OCT 		N_OCT % 2
#define L_OCT 		(N_OCT - R_OCT) / 2

#define N_ENV 		ENV_FILT_ORDER
#define R_ENV 		N_ENV % 2
#define L_ENV 		(N_ENV - R_ENV) / 2
#define LOW_PASS_CUTOFF	80.0
#define DAMP_SEG_SIZE	0.5
#define DAMP_STEP_SIZE	0.002
#define MAX_RT		100.0	
#define MIN_RT		0.005	
#define	BEST_S_NUM	100000
#define PAR_UP_BOUND	0.99999
#define PAR_LOW_BOUND	0.99
#define SQP_STEP	5

#define SEG_LEN		3 * 60		/* 3 minutes */
#define RT_DECAY	25		/* 25 = T25 decay; 35 = T35 decay */

/* Functions */
int get_wav_data(void);
int compute_rt(int samp, int array_size, int *store_start, int *store_end, 
	double *dr, double *a, double *b, double *alpha, double *mean_rt,
	double *rt_sd);
int optimum_model(double a, double b, double alpha, double dr, int samp, int n,
	double *chan_2);
int optimum_model_2(double **chan, int n, int chan_sz, int samp, double *lin);
void hanning(int len, double *hann);
int process_wav_data(int band, float *wav_data, SF_INFO input_info, 
	SNDFILE *input, unsigned long long frames, double **p_alpha,
	double **p_a, double **p_b, double **p_dr, int **p_len, int **p_start,
	int **p_end);
float * resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, 
	int samp_freq, const uint8_t **wav_data, int band, int *resamp_frms);
float * oct_filt_data(float *output_data, float band, float samp_freq, 
	int resamp_frames);
float * get_envelope(float samp_freq, float *resampled_wav, 
	int resamp_frames);
int ** apply_polyfit(float samp_freq, int resamp_frames, float *env, 
	int *s_e_size);
int ** choose_segments(float *var, int step_size, int len, int not_dumped,
	int *seg_num_els);
int perform_ml(int **start_end, float *env, int s_e_size, 
	float *filtered_wav, int resamp_frames, int samp_freq,
	double **p_alpha, double **p_a_par, double **p_b_par, double **p_dr,
	int **p_len, int **p_start, int **p_end);
int ml_fit(float *data_seg, int len, float samp_freq, float max_abs_seg,
	double *a_par, double *b_par, double *alph_par);
double alpha_opt(int n, const double *a, double *grad, void *nldv);
double par_3_opt(int n, const double *a, double *grad, void *nldv);
int plot_wav(float *wav_data, int channels, sf_count_t frames, int samprate);

/* Structs */
struct nl_extra_data {
	float *data_seg;
	int len;
	double a_val;
	double b_val;
};


/* Globals */
int num_bands = 8;
int octave_bands[] = {63, 125, 250, 500, 1000, 2000, 4000, 8000};
int samp_freq_per_band[] = {3000, 3000, 3000, 3000, 3000, 6000, 12000, 24000};

#if OCT_FILT_ORDER >= ENV_FILT_ORDER
	float butt_b[3 * (L_OCT + R_OCT)] = {0};
	float butt_a[3 * (L_OCT + R_OCT)] = {0};
#else
	float butt_b[3 * (L_ENV + R_ENV)] = {0};
	float butt_a[3 * (L_ENV + R_ENV)] = {0};
#endif

/* Parameters as calculated by matlab */
#if 0
// 63 hz
float but_b[11] = {
		0.0000019065480398620, 0, -0.0000095327401993098,
		0, 0.0000190654803986197, 0, -0.0000190654803986197,
		0, 0.0000095327401993098, 0, -0.0000019065480398620
};
float but_a[11] = {1.0000000000000, -9.6137442795514, 4.16762001930964, -107.2819840395506, 181.6014605265677, -211.2216382366775, 170.09536018279864};

//250 Hz
float but_b[11] = {
   	0.000125964124667, 0, -0.000629820623333, 0, 0.001259641246665, 0, 		-0.001259641246665, 0, 0.000629820623333, 0, -0.000125964124667};

float but_a[11] = {
	1.000000000000000, -7.610710971362601, 27.056928005646160, 
	-58.995953452911763, 87.229825526402408, -91.305190165139038, 
	68.499525784864559, -36.381079989157115, 13.104095798414871, 
	-2.895830250311406, 0.299189780176664};

float butter_b[8][11] = {
	{1.0e-05 * 0.0191, 0, 1.0e-05 * -0.0953, 0, 1.0e-05 * 0.1907, 0, 
		1.0e-05 * -0.1907, 0, 1.0e-05 * 0.0953, 0, 1.0e-05 * -0.0191},
	{1.0e-04 * 0.0511, 0, 1.0e-04 * -0.2553, 0, 1.0e-04 * 0.5105, 0, 
		1.0e-04 * -0.5105, 0, 1.0e-04 * 0.2553, 0, 1.0e-04 * -0.0511},
	{0.0001, 0, -0.0006, 0, 0.0013, 0, -0.0013, 0, 0.0006, 0, -0.0001},
	{0.0026, 0, -0.0128, 0, 0.0257, 0, -0.0257, 0, 0.0128, 0, -0.0026},
	{0.0419, 0, -0.2094, 0, 0.4189, 0, -0.4189, 0, 0.2094, 0, -0.0419},
	{0.0419, 0, -0.2094, 0, 0.4189, 0, -0.4189, 0, 0.2094, 0, -0.0419},
	{0.0419, 0, -0.2094, 0, 0.4189, 0, -0.4189, 0, 0.2094, 0, -0.0419},
	{0.0419, 0, -0.2094, 0, 0.4189, 0, -0.4189, 0, 0.2094, 0, -0.0419}
};
float butter_a[8][11] = {
	{1.0000, -9.6137, 41.6762, -107.2820, 181.6015, -211.2216, 170.9536, 
		-95.0705, 34.7671, -7.5499, 0.7393},
	{1.0000, -9.0799, 37.4091, -92.0798, 149.9398, -168.7663, 132.9719,
		-72.4197, 26.0936, -5.6173, 0.5487},
	{1.0000, -7.6107, 27.0569, -58.9960, 87.2298, -91.3052, 68.4995,
		-36.3811, 13.1041, -2.8958, 0.2992},
	{1.0000, -3.6281, 7.9585, -11.8968, 13.6277, -11.9992, 8.3692, -4.4707,
		1.8221, -0.5004, 0.0848},
	{1.0000, 4.3339, 8.2183, 9.6168, 8.4562, 5.9140, 3.0377, 1.0937,
		0.3038, 0.0597, 0.0026},
	{1.0000, 4.3339, 8.2183, 9.6168, 8.4562, 5.9140, 3.0377, 1.0937,
                0.3038, 0.0597, 0.0026},
	{1.0000, 4.3339, 8.2183, 9.6168, 8.4562, 5.9140, 3.0377, 1.0937,
                0.3038, 0.0597, 0.0026},
	{1.0000, 4.3339, 8.2183, 9.6168, 8.4562, 5.9140, 3.0377, 1.0937,
                0.3038, 0.0597, 0.0026}
};
#endif


int main(void)
{
	double ret;

	ret = get_wav_data();

	if (ret < 0) {
		printf("Error returned during calculations...\n");
		return -1;
	} else
		printf("Success! Exiting...\n");

	return 0;
}


int get_wav_data(void) {
	SF_INFO input_info;
	float *wav_data;
	unsigned long long long_ret = 0;	
	int *len, *store_start, *store_end;
	int band, array_size, ret = 0;
	double *a, *b, *alpha, *dr, mean_rt, rt_sd;

	/* Create arrays large enough to store values */
	a = calloc(input_info.frames, sizeof(double));
	b = calloc(input_info.frames, sizeof(double));
	alpha = calloc(input_info.frames, sizeof(double));
	dr = calloc(input_info.frames, sizeof(double));
	len = calloc(input_info.frames, sizeof(int));
	store_start = calloc(input_info.frames, sizeof(int));
	store_end = calloc(input_info.frames, sizeof(int));

	/* Open file */
	input_info.format = 0;
	SNDFILE *input = 
		sf_open("/home/kim/wav_samples/mono_bach_partita_e_maj.wav",
		//sf_open("/home/kim/BlindRT/DS_Science_S6_L1.wav",
		SFM_READ, &input_info);

	/* Print header data */
	printf("Frames: %llu\nSample Rate: %d\nChannels: %d\nFormat: 0x%X\n"
		"Sections: %d\nSeekable: %d\n",(unsigned long long) 
		input_info.frames, input_info.samplerate, input_info.channels, 
		input_info.format, input_info.sections, input_info.seekable);

	/* FIXME: insert check to make sure that the wav file is in the 
	 *	expected format: mono and 16 bit signed??? */

	/* Store data in array */
	wav_data = (float *) calloc(input_info.frames * input_info.channels 
		/ SPLIT, sizeof(float));

	if (wav_data == NULL) {
		printf("Memory allocation error\nExiting...\n");
		sf_close(input);
		return -1;
	}

	/* Octave band filtering */
	for (band = 2; band < 3; band++) {	

		/* Single channel, so can call this way */
	//	for (i = 0; i <= SPLIT; i++) {
			long_ret = sf_read_float(input, wav_data, 
					floor(input_info.frames / SPLIT));

			printf("DEBUG Read: %llu\n", long_ret);


#if 0
	/* Plot the wav with Gnuplot, for the craic */	
	ret = plot_wav(wav_data, input_info.channels, input_info.frames, 
		input_info.samplerate);	

	for (i = 0, k = 0; i < WINDOW_WIDTH; i += WINDOW_WIDTH * (1 - OVERLAP),
	k++) {
		for (j = 0; j < input_info.frames * input_info.channels ; j++) {
			wav_data[j];
		}
	}
#endif

			/* Process the wav data */
			printf("DEBUG: Start process_wav_data\n");
			ret = process_wav_data(band, wav_data, input_info, 
				input, long_ret, &alpha, &a, &b, &dr, &len,
				&store_start, &store_end);
			printf("DEBUG: End process_wav_data\n");

			/* If ret <= 0, error */
			/* If > 0, ret is the size of the returned arrays */
			if (ret <= 0) {
				free(wav_data);
				sf_close(input);
				return -1;
			} else
				array_size = ret;
				
			// SUM THE ARRAYS IF SPLIT

			printf("DEBUG: Start compute_rt\n");
			/* Compute the reverberation time (RT) */
			ret = compute_rt(samp_freq_per_band[band],
				array_size, store_start, store_end, dr, a, b,
				alpha, &mean_rt, &rt_sd);
	//	}
		printf("Mean RT for %d band: %lf\n", octave_bands[band], mean_rt);
		printf("RT Standard Dev for %d band: %lf\n", octave_bands[band], rt_sd);
	}


	/* Clean up */
//	sf_close(input);
	free(wav_data);
	return ret;
}


/* 
 * Compute RT for room 
 */
int compute_rt(int samp, int array_size, int *store_start, int *store_end, 
		double *dr, double *a, double *b, double *alpha, 
		double *mean_rt, double *rt_sd) {
	int n_seg, file_len, pos_to, pos_from, seg_len_n, i, j, k, l_reg, ret_sz;
	int min_5dB_index, min_rtdB_index, nan_count = 0;
	int n_chan = 4 * samp;
	double **chan, chan_opt[n_chan], chan_opt_log[n_chan];
	double *sd_rev_time, max, min_5dB, min_rtdB, **mat_a;
	SVDRec svd_mat = malloc(sizeof(SVDRec));
	DMat dmat_a_mat = malloc(sizeof(DMat));

	mat_a = calloc(2, sizeof(double *));
	chan = calloc(array_size, sizeof(double *));

	for (i = 0; i < array_size; i++)
		chan[i] = calloc(n_chan, sizeof(double));

	/* Find the end position of the last decay phase by finding max
	 *	of store_end[] */
	file_len = store_end[0];
	for (i = 1; i < array_size; i++) {
		if (store_end[i] > file_len)
			file_len = store_end[i];
	}

	/* seg_len_n is the number of samples in the length of time
	 *	chosen, seg_len. n_seg is the number of chunks, For
	 *	example, if seg_len = 180 secs and file_len = 6200
	 *	samples, then n_seg = 3 chunks (2 chunks of 180 secs and
	 *	1 shorter chunk). */
	seg_len_n = SEG_LEN * samp;
	n_seg = ceil((double)file_len / seg_len_n);

	/* Create array to store standard deviation values */
	sd_rev_time = calloc(n_seg, sizeof(double));

	for (i = 0; i < n_seg; i++) {
		/* Get start and end positions of each chunk */
		pos_from = i * seg_len_n;
		pos_to = (i+1) * seg_len_n - 1;

		/* In case of partial chunk */
		if (i == n_seg - 1)
			pos_to = file_len;

		k = 0;
		for (j = 0; j < array_size; j++) {
			/* Find decay rates < -25 dB */
			if (dr[j] < -25.0) {
				/* Does decay rate occur in this chunk? */
				if (store_start[j] >= pos_from && 
						store_end[j] < pos_to) {
					printf("DR: %lf\n", dr[j]);
					optimum_model(a[j], b[j], alpha[j], dr[j],
							samp, n_chan, 
							&chan[k][0]);

					k++;		
				}
			}
		}

		if (k <= 0)
			return 0;

		printf("DEBUG: start optimum_model_2: %d\n", k);
		ret_sz = optimum_model_2(chan, n_chan, k, samp, chan_opt);

		/* Get cumulative sum */
		chan_opt_log[ret_sz - 1] = pow(chan_opt[ret_sz - 1], 2);
		for (j = ret_sz - 2; j >= 0; j--) {
			chan_opt_log[j] = pow(chan_opt[j], 2) +
				chan_opt_log[j + 1];
		}

		// DEBUG note: chan_opt_log[0] close to matlab

		/* Get max before getting log */
		max = fabs(chan_opt_log[0]);
		for (j = 1; j < ret_sz; j++) {
			if (fabs(chan_opt_log[j]) > max)
				max = fabs(chan_opt_log[j]);
		}

		for (j = 0; j < ret_sz; j++) {
			/* Convert to decibels */
			chan_opt_log[j] = 10 * log10(chan_opt_log[j] / max);

			if (j < 10)
				printf("CHANOPTLOG %d: %lf\n", j, chan_opt_log[j]);

			/* Find -5 dB and -RT_DECAY dB decay points */
			if (j == 0) {
				min_5dB = fabs(chan_opt_log[0] - (-5));
				min_rtdB = fabs(chan_opt_log[0] - (-RT_DECAY));
				min_5dB_index = 0;
				min_rtdB_index = 0;
			} else {
				if (fabs(chan_opt_log[j] - (-5)) < min_5dB) {
					min_5dB = fabs(chan_opt_log[j] - (-5));
					min_5dB_index = j;
				}

				if (fabs(chan_opt_log[j]-(-RT_DECAY)) < min_rtdB) {
					min_rtdB = fabs(chan_opt_log[j] -
						(-RT_DECAY));
					min_rtdB_index = j;
				}
			}
		}

		l_reg = min_rtdB_index - min_5dB_index + 1;

		/* Allocate space for mat_a */
		for (k = 0; k < 2; k++)
			mat_a[k] = calloc(l_reg, sizeof(double));

		/* Get the RT */
		for (j = 0; j < l_reg; j++) {
			mat_a[0][j] = 1.0;
			mat_a[1][j] = (double) (min_5dB_index + j);
		}

		/* Get the SVD of mat_a */
		dmat_a_mat->rows = 2;
		dmat_a_mat->cols = l_reg;
		dmat_a_mat->value = mat_a;				

		svd_mat = svdLAS2A(svdConvertDtoS(dmat_a_mat), 2);

		/* Get reciprocal of S */
		for (k = 0; k < 2; k++)
			svd_mat->S[k] = 1.0 / svd_mat->S[k];

		/* Multiply modified S val by V */
		for (k = 0; k < 2; k++)
			for (j = 0; j < 2; j++)
				svd_mat->Vt->value[k][j] *= svd_mat->S[k];

		/* Matrix multiplcation of modified V with U */
		for (j = 0; j < l_reg; j++) {
			mat_a[0][j] = svd_mat->Ut->value[0][j] *
				svd_mat->Vt->value[0][0];

			mat_a[1][j] = svd_mat->Ut->value[1][j] *
				svd_mat->Vt->value[1][1];

			mat_a[0][j] += svd_mat->Ut->value[1][j] *
					svd_mat->Vt->value[1][0];

			mat_a[1][j] += svd_mat->Ut->value[0][j] *
					svd_mat->Vt->value[0][1];
        	}
		printf("Foo\n");

		/* Multiply pseudoinverse by sig */
		for (k = 0; k < 2; k++)
			for (j = 0; j < l_reg; j++)
				mat_a[k][j] *= chan_opt_log[j + min_5dB_index];

		/* 
		 * If RT for segment != nan, add to rev_time
		 */
		if(!isnan(-60.0 / (samp * mat_a[1][0]))) {
			*mean_rt += -60.0 / (samp * mat_a[1][0]);
			sd_rev_time[i - nan_count] = -60.0 / (samp * mat_a[1][0]);
		} else
			nan_count++;

		/* Free array */
		for (k = 0; k < 2; k++)
			free(mat_a[k]);
	}

	/* Get mean RT and standard dev of RT */
	*mean_rt /= (n_seg - nan_count);

	for (i = 0; i < n_seg - nan_count; i++)
		*rt_sd = pow(fabs(sd_rev_time[i] - *mean_rt), 2);

	*rt_sd = sqrt(*rt_sd / (n_seg - nan_count - 1));

	/* Free array */
	free(mat_a);

	return 0;
}


int optimum_model(double a, double b, double alpha, double dr, int samp, int n,
		double *chan_2) {
	int i;
	double a_pow, b_pow, max;

	for (i = 0; i < n; i++) {
		a_pow = pow(a, i);
		b_pow = pow(b, i);
		chan_2[i] = alpha * a_pow + (1 - alpha) * b_pow;

		if (i == 0)
			max = fabs(chan_2[i]);
		else {
			if (fabs(chan_2[i]) > max)
				max = fabs(chan_2[i]);
		}
	}

	for (i = 0; i < n; i++)
		chan_2[i] /= max;

	return 0;
}


int optimum_model_2(double **chan, int n, int chan_sz, int samp, double *lin) {
	double winn = 0.1;
	double over = 0.5;
	int nn = winn * samp;
	int i, j, k, n1, min_index, l_reg; 
	int n0 = floor((1.0 - over) * nn);
	int n_sect = floor(((double)n - nn) / n0);
	double n2, min, chan_sum[chan_sz], last[n]; 
	double win[nn * n_sect];

	for (i = 0, n1 = 0; i < n_sect; i++) {
		n2 = n1 + nn;
		l_reg = n2 - n1;

		for (j = 0; j < chan_sz; j++) {
			chan_sum[j] = 0;
			for (k = n1; k < n2; k++) 
				chan_sum[j] += pow(chan[j][k], 2);

			/* Get the minimum */
			if (j == 0) {
				min = chan_sum[0];
				min_index = 0;
			}
			else {
				if (chan_sum[j] < min) {
					min = chan_sum[j];
					min_index = j;
				}
			}
		}

		// DEBUG note: min and min_index match matlab 250 Hz

		if (i == 0) {
			for (k = 0; k < n; k++) {
				if (k >= n1 && k < n2)
					lin[k] = chan[min_index][k];

				last[k] = chan[min_index][k];
			}
		} else {
			hanning(l_reg, win);

			for (k = l_reg / 2 - 1; k < l_reg; k++)
				win[k] = 1.0;

			for (k = 0; k < n; k++) {
				if (k >= n1 && k < n2) {
					lin[k] = last[k] * fabs(win[k] - 1.0) +
						chan[min_index][k] * win[k];
				}

				last[k] = chan[min_index][k];
			}
		}	
		n1 += n0;
	}
	//DEBUG note: lin[0-9] matches MATLAB

	return (int) n2;
}

void hanning(int len, double *hann) {
	int n;

	for (n = 0; n < len; n++)
		hann[n] = 0.5 * (1.0 - cos(2.0 * M_PI * n / (len - 1.0)));

	return;
}

int process_wav_data(int band, float *wav_data, SF_INFO input_info, 
		SNDFILE *input, unsigned long long frames, double **p_alpha,
		double **p_a, double **p_b, double **p_dr, int **p_len,
		int **p_start, int **p_end) {
	int ret, s_e_size, resamp_frames = 0, **start_end;
	float *resampled_wav, *filtered_wav, *env;
	double *alpha, *a, *b, *dr;
	int *len, *store_start, *store_end;

	printf("DEBUG: resample init\n");

	/* Specifications that apply to all subbands */
	SwrContext *resamp = swr_alloc();

	av_opt_set_int(resamp, "in_channel_layout", AV_CH_LAYOUT_MONO, 0);
	av_opt_set_int(resamp, "out_channel_layout", AV_CH_LAYOUT_MONO, 0);
	av_opt_set_int(resamp, "in_sample_rate", input_info.samplerate, 0);
	av_opt_set_sample_fmt(resamp, "in_sample_fmt", AV_SAMPLE_FMT_FLT, 0);
	av_opt_set_sample_fmt(resamp, "out_sample_fmt", AV_SAMPLE_FMT_FLT, 0);


	av_opt_set_int(resamp, "out_sample_rate", 
		samp_freq_per_band[band], 0); 

	ret = swr_init(resamp);

	if (ret < 0) {
		printf("Resample initialisation error\nExiting...\n");
		free(wav_data);
		sf_close(input);
		return -1;
	}

	/* Resample the audio */
	resampled_wav = resamp_wav_data(resamp, input_info.samplerate, 
		(uint64_t) frames, samp_freq_per_band[band],
		(const uint8_t **) &wav_data, octave_bands[band], 
		&resamp_frames);

	if (resamp_frames <= 0 || resampled_wav == NULL) {
		printf("Resample error\nExiting...\n");
		free(wav_data);
		sf_close(input);
		return -1;
	}

	/* Octave band filtering */
	filtered_wav = oct_filt_data(resampled_wav, 
		(float) octave_bands[band], (float) samp_freq_per_band[band],
		resamp_frames);
	
	/* Obtain signal envelope */
	env = get_envelope((float) samp_freq_per_band[band], filtered_wav, 
		resamp_frames);

	/* Polyfit algorithm & get decay segments */
	start_end = apply_polyfit((float) samp_freq_per_band[band], 
		resamp_frames, env, &s_e_size);

	/* Perform Maximum Liklihood Estimation */
	ret = perform_ml(start_end, env, s_e_size, filtered_wav, 
		resamp_frames, samp_freq_per_band[band], &alpha, &a, &b, &dr,
		&len, &store_start, &store_end);

	/* Return pointers to caller */
	*p_alpha = alpha;
	*p_a = a;
	*p_b = b;
	*p_dr = dr;
	*p_len = len;
	*p_start = store_start;
	*p_end = store_end;

	if (ret != 0)
		return -1;
	else
		return s_e_size;
}

float * resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, int samp_freq, const uint8_t **wav_data, int band, int *resamp_frms) {
	uint8_t *output_data;
	int ret, resamp_frames;

	resamp_frames = av_rescale_rnd(swr_get_delay(resamp, in_rate) + 
		num_frames, samp_freq, in_rate, AV_ROUND_UP);

	ret = av_samples_alloc(&output_data, NULL, 1, resamp_frames, 
		AV_SAMPLE_FMT_FLT, 0);

	if (ret < 0) {
		printf("Output buffer alloc failure\nExiting...\n");
		return NULL;
	}

	resamp_frames = swr_convert(resamp, &output_data, resamp_frames,
		wav_data, num_frames);

	if (resamp_frames < 0) {
		printf("Resampling failure\nExiting...\n");
		av_freep(&output_data);
		return NULL;
	}


#if 0
	/* Generate wav to check resampled data is not garbage */
	SF_INFO output_info;
	output_info.samplerate = samp_freq;
	output_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	output_info.channels = 1;

	SNDFILE* otpt = sf_open("/home/kim/wav_samples/Z.wav", 
		SFM_WRITE, &output_info);

	sf_write_float(otpt, (float *) output_data, resamp_frames);
#endif


	swr_free(&resamp);
//FIXME: free output_data ptr at end!
//	av_freep(&output_data);

	/* Store # frames in pointer */
	*resamp_frms = resamp_frames;

	/* Return pointer to beginning of data */
	return (float *) output_data;
}

float * oct_filt_data(float *resamp_data, float band, float samp_freq, int resamp_frames) {
	int i;
	float *filtered = calloc(resamp_frames, sizeof(float)); 

	memset(butt_b, 0, sizeof(butt_b));
	memset(butt_a, 0, sizeof(butt_a));

	/* Create bandpass filter - 7th and 8th parameter are ignored */
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_BANDPASS, 
		LIQUID_IIRDES_SOS, OCT_FILT_ORDER, 
		band / (sqrt(2.0) * samp_freq), 
		band / samp_freq, 1.0, 1.0, butt_b, butt_a);
	//printf("res: %d\n\n", iirdes_isstable(butt_b, butt_a, 3 * (L_OCT + R_OCT)));

//	for (i=0; i<(3 * (L_OCT + R_OCT)); i++)
//		printf("b: %f and a: %f\n", butt_b[i], butt_a[i]);		

	iirfilt_rrrf f_obj_oct = iirfilt_rrrf_create_sos(butt_b, butt_a, (L_OCT + R_OCT));


	/* Filter data */
	for (i = 0; i < resamp_frames; i++) 
		iirfilt_rrrf_execute(f_obj_oct, resamp_data[i], &filtered[i]);

	/* Free filter and reset coefficient arrays */
	iirfilt_rrrf_destroy(f_obj_oct);

	/* fiiltered data: return to previous function */
	return filtered;
}


float * get_envelope(float samp_freq, float *filtered_wav, 
		int resamp_frames) {
	float complex *filt_complex = calloc(resamp_frames, sizeof(float complex));
	float complex *tmp = calloc(resamp_frames, sizeof(float complex));
	float complex *tmp2 = calloc(resamp_frames, sizeof(float complex));
	float complex *hilb = calloc(resamp_frames, sizeof(float complex));
	float *env = calloc(resamp_frames, sizeof(float));
	int i, j;

	//FIXME add calloc NULL check	

	memset(butt_b, 0, sizeof(butt_b));
	memset(butt_a, 0, sizeof(butt_a));

	/* Create Low Pass Filter - 6th, 7th and 8th param are ignored */
	printf("DEBUG: Low pass filter and Hilbert transform\n");
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_LOWPASS, 
		LIQUID_IIRDES_SOS, ENV_FILT_ORDER, 
		LOW_PASS_CUTOFF / samp_freq, 0.1, 
		1.0, 1.0, butt_b, butt_a);

	iirfilt_crcf f_obj_env = iirfilt_crcf_create_sos(butt_b, butt_a, (L_ENV + R_ENV));

	fftplan fft_pln = fft_create_plan(resamp_frames, (float complex *) 
		filt_complex, (float complex *) tmp, LIQUID_FFT_FORWARD, 0);
	fftplan ifft_pln = fft_create_plan(resamp_frames, (float complex *) tmp,
 		(float complex *) hilb, LIQUID_FFT_BACKWARD, 0);

	/* Get the Hilbert transform - really it's the analytic signal as the
		Hilbert transform is stored in the imaginary part of the soln */
	for (i = 0; i < resamp_frames; i+= 1)
		filt_complex[i] = (float complex) filtered_wav[i];

	fft_execute(fft_pln);	

	for (i = 0; i < resamp_frames; i+= 1) {
		j = i + 1;

		if (resamp_frames % 2 == 0) {
			/* Even number of samples */
			/* Set if value should be non-zero (array contains zeros)*/
			if (j == 1 || j == ((resamp_frames / 2) + 1))
				tmp2[i] = 1.0;
			else if (j >= 2 && j <= (resamp_frames / 2))
				tmp2[i] = 2.0;
		} else {
			/* Odd number of samples */
			/* Set if value should be non-zero (array contains zeros)*/
			if (j == 1)
				tmp2[i] = 1.0;
			else if (j >= 2 && j <= (resamp_frames + 1) / 2)
				tmp2[i] = 2.0;
		}

		tmp[i] *= tmp2[i];
	}


	fft_execute(ifft_pln);

	for (i = 0; i < resamp_frames; i++) {
		/* Normalisation correction */
		hilb[i] *= (1.0 / resamp_frames);

		/* Get envelope */
		hilb[i] = csqrtf(filt_complex[i] * filt_complex[i] + 
			hilb[i] * hilb[i]);

		/* Apply low pass filter to envelope */
		iirfilt_crcf_execute(f_obj_env, hilb[i], &tmp[i]);
		env[i] = cabsf(tmp[i]);
//		if (i < 100)
//			printf("%d: %lf\n", i, env[i]);
	}


	//FIXME: why is this causing a fre pointer error???
	//fft_destroy_plan(fft_pln);
	//DONT FREE ENV
	free(tmp2);
	free(tmp);
	free(hilb);
	free(filt_complex);
	iirfilt_crcf_destroy(f_obj_env);

	/* envelope pointer: return to previous function */
	return env;
}


int ** apply_polyfit(float samp_freq, int resamp_frames, float *env, int *s_e_size) {
	int i, j, seg_size, step_size, num_segs, poly_coeff = 2, dump = 0;
	int **seg_index, invalid_coeff = 0, actual_seg_num_els = 0, max = 0;
	int **start_end; 
	float *poly_param, **poly_res, *poly_seg, max_poly_res, min_poly_res;
	float *var;
	float *log_env = calloc(resamp_frames, sizeof(float));
	//FIXME: add calloc NULL check

	for (i = 0; i < resamp_frames; i++)
		log_env[i] = log10f(env[i]);

	seg_size = (int) floorf(samp_freq * DAMP_SEG_SIZE);

	step_size = (int) floorf(samp_freq * DAMP_STEP_SIZE);

	num_segs = (int) floorf((resamp_frames - seg_size) / step_size + 1);

	//FIXME: add calloc NULL check
	poly_param = calloc(seg_size, sizeof(float));
	poly_res = calloc(num_segs, sizeof(float *));
	poly_seg = calloc(seg_size, sizeof(float));
	var = calloc(num_segs, sizeof(float));

	for (i = 0; i < num_segs; i++) 
		poly_res[i] = calloc(poly_coeff, sizeof(float));

	/* Fill in Polyfit function parameter */
	for (i = 0; i < seg_size; i++)
		poly_param[i] = i + 1;

	max_poly_res = -log10f(exp(6.91 / MAX_RT / samp_freq));
	min_poly_res = -log10f(exp(6.91 / MIN_RT / samp_freq));

	for (i = 0; i < num_segs; i++) {
		for (j = i * step_size; j < i * step_size + seg_size; j++)
			poly_seg[j - i * step_size] = log_env[j];
	
		polyf_fit(poly_param, poly_seg, seg_size, &poly_res[i][0], 
			poly_coeff);

		//TODO: CHECK MATLAB POLYFIT COEFFS
		//TODO: CHECK MATLAB VALUES UP TO HERE
		//DONE ABOVE - SEEMED FINE BUT POLY COEFFS REVERESED...
//		for (j = 0; j < poly_coeff; j++) {
//			printf("%lf", poly_res[i][j]);
//			if (j == 0)
//				printf(",");
//			else
//				printf("\n");
//		}

		for (j = 0; j < poly_coeff; j++) {
			//FIXME eeeh floating point comparisons are DODGEY in C
			if (poly_res[i][j] < min_poly_res || 
				poly_res[i][j] > max_poly_res) {
				invalid_coeff++;
			}
		}
		if (invalid_coeff  == poly_coeff) {
			var[i] = -1.0;
			dump++;
		}
		invalid_coeff = 0;		
	}

	/* Choose segments */
	seg_index = choose_segments(var, step_size, num_segs, num_segs - dump,
		&actual_seg_num_els);

	/* Number of decay segments should not be greater than BEST_S_NUM */
	if (actual_seg_num_els > BEST_S_NUM)
		actual_seg_num_els = BEST_S_NUM;

	start_end = calloc(2, sizeof(int *));
	for (i = 0; i < 2; i++)
		start_end[i] = calloc(actual_seg_num_els, sizeof(int));

	for (i = 0; i < actual_seg_num_els; i++) {
		max = 0;
		for (j = 0; j < actual_seg_num_els; j++) {
			if (seg_index[j][1] > seg_index[max][1])
				max = j;
		}

		start_end[0][i] = seg_index[max][0] * step_size;
		start_end[1][i] = seg_index[max][0] * step_size + 
			seg_size + seg_index[max][1] * step_size - 1;

		/* Remove longest segment - it has been chosen */
		seg_index[max][1] = 0;
		
	}

	/* Free 2d arrays */
	for (i = 0; i < num_segs; i++) 
		free(poly_res[i]);
	
	free(poly_seg);
	free(poly_res);
	free(var);
	free(poly_param);
	free(log_env);

	*s_e_size = actual_seg_num_els;

	return start_end;
}

int ** choose_segments(float *var, int step_size, int len, int not_dumped, 
		int *seg_num_els) {
	int i, k, j = 0;
	int **seg_index;

	/* Max possible size of seg_index */
	seg_index = calloc(not_dumped, sizeof(int *));

	for (i = 0; i < not_dumped; i++) 
		seg_index[i] = calloc(2, sizeof(int));

	/* Set decay segments and number of continuous segments in array */
	for (i = 0; i < len - 1; i += k) {
		k = 1;

		/* Floating point method of checking if var[i] == 0 */
		if (fabs(var[i]) < 0.0001) {
			seg_index[j][0] = i;
			seg_index[j][1] = 0;

			while (((i + k) < len) && fabs(var[i + k]) < 0.0001) {
				seg_index[j][1]++;
				k++;
			}
			j++;
			k += step_size;
		}
	}

	*seg_num_els = j;
	/* FIXME: don't free seg_index... what am I supposed to do about it? */
	return seg_index;
}

int perform_ml(int **start_end, float *env, int s_e_size, float *filtered_wav, 
		int resamp_frames, int samp_freq, double **p_alpha, 
		double **p_a_par, double **p_b_par, double **p_dr, int **p_len,
		int **p_start, int **p_end) {
	int sz, i, j, k, len, min_abs_filt_seg, max_abs_filt_seg, ret = 0;
	float *segment, *env_seg, *abs_filt_seg, *seg_seg2, max_seg2;
	float *abs_filtered = calloc(resamp_frames, sizeof(float));
	int *store_start = calloc(s_e_size, sizeof(int)); 
	int *store_end = calloc(s_e_size, sizeof(int));
	int *len_store = calloc(s_e_size, sizeof(int));
	double *a_par = calloc(s_e_size, sizeof(double));
	double *b_par = calloc(s_e_size, sizeof(double));
	double *alpha = calloc(s_e_size, sizeof(double));
	double *dr = calloc(s_e_size, sizeof(double));
	double *chan = calloc(6 * samp_freq, sizeof(double));
	double *y = calloc(6 * samp_freq, sizeof(double));
	double max_y;

	/* Get absolute vaule of filtered wav */
	for (i = 0; i < resamp_frames; i++)
		abs_filtered[i] = fabs(filtered_wav[i]);		

	/* Largest possible size of arrays - first pair of elements */
	sz = start_end[1][0] - start_end[0][0];
	segment = calloc(sz, sizeof(float));
	env_seg = calloc(sz, sizeof(float));
	abs_filt_seg = calloc(sz, sizeof(float));
	seg_seg2 = calloc(sz + sz - roundf((DAMP_SEG_SIZE / 2) * samp_freq), 
		sizeof(float));

	for (i = 0; i < s_e_size; i++) {
		len = abs(start_end[1][i] - start_end[0][i]);

		/* Store decay segment in new array */
		for (k = 0, j = start_end[0][i]; j < start_end[1][i]; j++, k++)
		{
			segment[k] = filtered_wav[j];
			env_seg[k] = env[j];
			abs_filt_seg[k] = abs_filtered[j];
		}
		// DEBUG: above values match matlab!

		/* Get location of minimum */
		min_abs_filt_seg = len - round((DAMP_SEG_SIZE / 2) * 
			samp_freq) - 1;
		
		for (j = min_abs_filt_seg; j < len; j++) {
			if (abs_filt_seg[j] < abs_filt_seg[min_abs_filt_seg])
				min_abs_filt_seg = j;
		}
		
		/* Get location of maximum */
		max_abs_filt_seg = 0;
		for (j = 0; j < round((DAMP_SEG_SIZE / 2) * samp_freq); j++) {
			if (abs_filt_seg[j] > abs_filt_seg[max_abs_filt_seg])
				max_abs_filt_seg = j;
		}

		/* Store new fine-tuned decay segment */
		for (k = 0, j = max_abs_filt_seg; j <= min_abs_filt_seg; k++, 
			j++)
			seg_seg2[k] = segment[j];

		/* NB: len_store[] contains INDEX to last element, not # of
			elements in array */
		len_store[i] = k - 1;
		
		k = 0;
		for (j = 0; j <= len_store[i]; j++) {
			if (fabs(seg_seg2[j]) > fabs(seg_seg2[k]))
				k = j; 
		}
		
		max_seg2 = fabs(seg_seg2[k]);	
		for (j = 0; j <= len_store[i]; j++)
			seg_seg2[j] /= max_seg2;

		/* Save fine-tuned start and end locations */
		// FIXME: are the -1's OK here? Compare to MATLAB
		store_start[i] = start_end[0][i] + max_abs_filt_seg;
		store_end[i] = start_end[0][i] + max_abs_filt_seg + 
				len_store[i] - 1;

#if 0
		/* ML fitting of decay model to the data */
		ret = ml_fit(seg_seg2, len_store[i], (float) samp_freq, 
			fabs(seg_seg2[k]), &a_par[i], &b_par[i], &alpha[i]);

                /* DEBUG: using matlab values to test */
                if (i == 0) {
                        a_par[i] = 0.999885193054;
                        b_par[i] = 0.995180243256;
                        alpha[i] = 1.000000000000;
                }
                if (i == 1) {
                        a_par[i] = 0.999878346327;
                        b_par[i] = 0.990000000000;
                        alpha[i] = 0.751467410761;
                }
                if (i == 2) {
                        a_par[i] = 0.9998735051930;
                        b_par[i] = 0.9900000000000;
                        alpha[i] = 0.7863959195535;
                }
                if (i == 3) {
                        a_par[i] = 0.999990000000;
                        b_par[i] = 0.998207538979;
                        alpha[i] = 0.667957424959;
                }
                /*******************************/
#endif
                if (i == 0) {
                        a_par[i] = 0.999094589804397;
                        b_par[i] = 0.994872826187869;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 1) {
                        a_par[i] = 0.999180097014578;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.212695672167877;
                }
                if (i == 2) {
                        a_par[i] = 0.999713323430732;
                        b_par[i] = 0.996342970347192;
                        alpha[i] = 0.374014803581545;
                }
                if (i == 3) {
                        a_par[i] = 0.999687612639516;
                        b_par[i] = 0.995605995258264;
                        alpha[i] = 0.091897918866414;
                }
                if (i == 4) {
                        a_par[i] = 0.994711785713134;
                        b_par[i] = 0.999024970844051;
                        alpha[i] = 0.0;
                }
                if (i == 5) {
                        a_par[i] = 0.995043336608471;
                        b_par[i] = 0.998847922767537;
                        alpha[i] = 0.0;
                }
                if (i == 6) {
                        a_par[i] = 0.994394826394472;
                        b_par[i] = 0.999019876183459;
                        alpha[i] = 0.0;
                }
                if (i == 7) {
                        a_par[i] = 0.999791266550412;
                        b_par[i] = 0.998029286558849;
                        alpha[i] = 0.696376397046217;
                }
                if (i == 8) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999509030932049;
                        alpha[i] = 0.0;
                }
                if (i == 9) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.999272637030537;
                        alpha[i] = 0.064216596905715;
                }
                if (i == 10) {
                        a_par[i] = 0.994031236911420;
                        b_par[i] = 0.999211181903717;
                        alpha[i] = 0.864570140222733;
                }
                if (i == 11) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.999339657895692;
                        alpha[i] = 0.225343158000227;
                }
                if (i == 12) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999494985340351;
                        alpha[i] = 0.319145306059731;
                }
                if (i == 13) {
                        a_par[i] = 0.999078051623354;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 14) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.998898257531995;
                        alpha[i] = 0.496185978717731;
                }
                if (i == 15) {
                        a_par[i] = 0.998182309046916;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.906716469849577;
                }
                if (i == 16) {
                        a_par[i] = 0.996974645165489;
                        b_par[i] = 0.999072395847762;
                        alpha[i] = 0.712963359594970;
                }
                if (i == 17) {
                        a_par[i] = 0.997944971151486;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.943110682342377;
                }
                if (i == 18) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.998357932928653;
                        alpha[i] = 0.140224235031150;
                }
                if (i == 19) {
                        a_par[i] = 0.999433444880068;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.376218010146496;
                }
                if (i == 20) {
                        a_par[i] = 0.999567015115998;
                        b_par[i] = 0.994720278601704;
                        alpha[i] = 0.144445669562915;
                }
                if (i == 21) {
                        a_par[i] = 0.999533092883626;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.104621032127595;
                }
                if (i == 22) {
                        a_par[i] = 0.999251569288705;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.172358075776195;
                }
                if (i == 23) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.998505804013255;
                        alpha[i] = 0.204383717664616;
                }
                if (i == 24) {
                        a_par[i] = 0.999433757942988;
                        b_par[i] = 0.994753977287462;
                        alpha[i] = 0.147747268686187;
                }
                if (i == 25) {
                        a_par[i] = 0.999158814183372;
                        b_par[i] = 0.995120948211300;
                        alpha[i] = 0.175261597302704;
                }
                if (i == 26) {
                        a_par[i] = 0.994708516594345;
                        b_par[i] = 0.998819933198468;
                        alpha[i] = 0.0;
                }
                if (i == 27) {
                        a_par[i] = 0.994142863363533;
                        b_par[i] = 0.999114614982346;
                        alpha[i] = -0.0;
                }
                if (i == 28) {
                        a_par[i] = 0.995798702641769;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.716848019811558;
                }
                if (i == 29) {
                        a_par[i] = 0.999495177542700;
                        b_par[i] = 0.997492620171249;
                        alpha[i] = 0.098172655772365;
                }
                if (i == 30) {
                        a_par[i] = 0.999255466614124;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.164653437368793;
                }
                if (i == 31) {
                        a_par[i] = 0.999627093659115;
                        b_par[i] = 0.997675558319822;
                        alpha[i] = 0.105743555931733;
                }
                if (i == 32) {
                        a_par[i] = 0.997746778185184;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.789756778534100;
                }
                if (i == 33) {
                        a_par[i] = 0.999417503812517;
                        b_par[i] = 0.999384258506497;
                        alpha[i] = 0.0;
                }
                if (i == 34) {
                        a_par[i] = 0.994566864788402;
                        b_par[i] = 0.999825342989669;
                        alpha[i] = 0.920800340297833;
                }
                if (i == 35) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.684330132699412;
                }
                if (i == 36) {
                        a_par[i] = 0.999161897452354;
                        b_par[i] = 0.999161900460670;
                        alpha[i] = 0.539573163424888;
                }
                if (i == 37) {
                        a_par[i] = 0.999780374941031;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.589498800657416;
                }
                if (i == 38) {
                        a_par[i] = 0.999890062058455;
                        b_par[i] = 0.999389970934506;
                        alpha[i] = 0.0;
                }
                if (i == 39) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.995706384699949;
                        alpha[i] = 0.533102880748713;
                }
                if (i == 40) {
                        a_par[i] = 0.997568226638431;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.874103377870471;
                }
                if (i == 41) {
                        a_par[i] = 0.990857493470198;
                        b_par[i] = 0.999799701598015;
                        alpha[i] = 0.799692841907514;
                }
                if (i == 42) {
                        a_par[i] = 0.999247828449652;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.538917883798136;
                }
                if (i == 43) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.998633369316183;
                        alpha[i] = 0.228572794714481;
                }
                if (i == 44) {
                        a_par[i] = 0.998995686598252;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.705943149406000;
                }
                if (i == 45) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.993488253215005;
                        alpha[i] = 0.332440637310982;
                }
                if (i == 46) {
                        a_par[i] = 0.999065049707584;
                        b_par[i] = 0.999064780161723;
                        alpha[i] = 0.306655707655123;
                }
                if (i == 47) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.992209081841129;
                        alpha[i] = 0.080058153911342;
                }
                if (i == 48) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999584857043679;
                        alpha[i] = 0.649346441588718;
                }
                if (i == 49) {
                        a_par[i] = 0.998270626706036;
                        b_par[i] = 0.998270692767028;
                        alpha[i] = 0.904538450398288;
                }
                if (i == 50) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.197876688359996;
                }
                if (i == 51) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999632236731992;
                        alpha[i] = 0.545485077299886;
                }
                if (i == 52) {
                        a_par[i] = 0.999979997411215;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.953197444323481;
                }
                if (i == 53) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.329032231385471;
                }
                if (i == 54) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 55) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.653258304572190;
                }
                if (i == 56) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 57) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.452280924341951;
                }
                if (i == 58) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.154125300792170;
                }
                if (i == 59) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 60) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.607410603268955;
                }
                if (i == 61) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.662868111795531;
                }
                if (i == 62) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999009381240013;
                        alpha[i] = 0.229121004431434;
                }
                if (i == 63) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 64) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.491347573920274;
                }
                if (i == 65) {
                        a_par[i] = 0.999822132974474;
                        b_par[i] = 0.999320103776189;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 66) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.633637992072218;
                }
                if (i == 67) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.616438010912253;
                }
                if (i == 68) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 0.155730142294404;
                }
                if (i == 69) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.521850366173101;
                }
                if (i == 70) {
                        a_par[i] = 0.999183098782412;
                        b_par[i] = 0.990045486335074;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 71) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999255589697719;
                        alpha[i] = 0.583290771847080;
                }
                if (i == 72) {
                        a_par[i] = 0.997895948312226;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.442090008013401;
                }
                if (i == 73) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.077965525936381;
                }
                if (i == 74) {
                        a_par[i] = 0.997895948312226;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.442090008013401;
                }
                if (i == 75) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 76) {
                        a_par[i] = 0.999990000000000;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 77) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999677231235847;
                        alpha[i] = 0.361290517355847;
                }
                if (i == 78) {
                        a_par[i] = 0.997614552490232;
                        b_par[i] = 0.990000000000000;
                        alpha[i] = 1.000000000000000;
                }
                if (i == 79) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999969121708210;
                        alpha[i] = 0.106665611482328;
                }
                if (i == 80) {
                        a_par[i] = 0.990000000000000;
                        b_par[i] = 0.999677231235847;
                        alpha[i] = 0.361290517355847;
                }
                if (i == 81) {
                        a_par[i] = 0.997462545276869;
                        b_par[i] = 0.999990000000000;
                        alpha[i] = 0.931579173775712;
		}

		/* FIXME: use ret to check for error */
		/* FIXME: fix all error cases - clean exit and all that */

		/* Cumpute decay curves using ML-calculated params */
		for (j = 0; j < 6 * samp_freq; j++) {
			chan[j] = alpha[i] * pow(a_par[i], j) + 
				(1.0 - alpha[i]) * pow(b_par[i], j);
		}

		/* Get the reverse cumsum of chan[] and assign to y in reverse */
		y[6 * samp_freq - 1] = pow(chan[6 * samp_freq - 1], 2);
		for (j = 6 * samp_freq - 2; j >= 0; j--) 
			y[j] = pow(chan[j], 2) + y[j + 1];

		/* Get max of fabs(y) */
		max_y = fabs(y[0]);
		for (j = 1; j < 6 * samp_freq; j++) {
			if (fabs(y[j]) > max_y)
				max_y = fabs(y[j]);		
		}

		/* Normalise y first value will be max; stores whole cumsum */
                for (j = 0; j < 6 * samp_freq; j++) 
			y[j] = 10.0 * log10(y[j] / max_y);

		dr[i] = y[len_store[i]];	
		printf("DEUBG DR %d: %lf\n",i, dr[i]);

		/* Reset to zero b/c arrays are longer than they need to be */
		memset(segment, 0, sz * sizeof(float));
		memset(seg_seg2, 0, (sz + sz - roundf((DAMP_SEG_SIZE / 2) * 
			samp_freq)) * sizeof(float));
		memset(env_seg, 0, sz * sizeof(float));
		memset(abs_filt_seg, 0, sz * sizeof(float));
	}

	/* Pointers used by caller cannot be freed. They are alpha[], a[],
		b[], len_store[], store_start[], store_end[] and dr[] */
	*p_alpha = alpha;
	*p_a_par = a_par;
	*p_b_par = b_par;
	*p_dr = dr;
	*p_len = len_store;
	*p_start = store_start;
        *p_end = store_end;

	/* Free pointers that are no longer needed */
	free(segment);
	free(seg_seg2);
	free(env_seg);
	free(abs_filt_seg);
	free(abs_filtered);
	free(chan);
	free(y);

	if (ret < 0)
		return -1;
	else
		return 0;
}

int ml_fit(float *data_seg, int len, float samp_freq, float max_abs_seg,
		double *a_par, double *b_par, double *alph_par) {
	double min, max, interval, j, gmax_val, x_fine[3], fine_val;
	double *coarse_grid, lb[3], ub[3];
	double like[SQP_STEP][SQP_STEP] = {0}, alpha[SQP_STEP][SQP_STEP] = {0};
	int gmax_pos[2], i, k, el_cg, ret = 0;
	struct nl_extra_data nld;

	lb[0] = 0.0;
	ub[0] = 1.0;

	nld.data_seg = data_seg;
	nld.len = len;

	nlopt_opt nl_obj1 = nlopt_create(NLOPT_LN_COBYLA, 1);
	nlopt_set_lower_bounds1(nl_obj1, lb[0]);
	nlopt_set_upper_bounds1(nl_obj1, ub[0]);
	nlopt_set_maxeval(nl_obj1, 300);
	nlopt_set_min_objective(nl_obj1, (nlopt_func) alpha_opt, (void *) &nld);


	min = (-6.91 / log(PAR_LOW_BOUND)) / 3000.0;
	min = exp(-6.91 / (samp_freq * min));

	max = (-6.91 / log(PAR_UP_BOUND)) / 3000.0;	
	max = exp(-6.91 / (samp_freq * max));
	
	interval = (max - min) / ((float) SQP_STEP - 1.0);
	el_cg = round((max - min) / interval) + 1;
	coarse_grid = calloc(el_cg, sizeof(double));

	/* Fill in coarse grid */
	for (i = 0, j = min; i < el_cg; i++, j += interval)
		coarse_grid[i] = j;

	// There is another max(abs(x)) division here in MATLAB
	// 	Makes no diff - dividing by 1.0
	for (i = 0; i < len; i++) 
		data_seg[i] /= max_abs_seg;

	/* Coarse grid minimisation search */
	for (i = 0; i < SQP_STEP; i++) {
		nld.b_val = coarse_grid[i];
 
		for (k = 0; k < SQP_STEP; k++) {
			alpha[k][i] = 0.5;
			nld.a_val = coarse_grid[k];

			ret |= nlopt_optimize(nl_obj1, &alpha[k][i], 
				&like[k][i]);

			like[k][i] *= -1;

			/* Store value and position ([][]) of global max */
			if (i == 0 && k == 0) {
				gmax_val = like[k][i];
				gmax_pos[0] = k;
				gmax_pos[1] = i;	
			} else {
				if (like[k][i] > gmax_val) {
					gmax_val = like[k][i];
					gmax_pos[0] = k;
					gmax_pos[1] = i;	
				}
			}

//			printf("k: %d, i: %d, ret: %d, a: %lf like: %le coarse: %lf\n", k+1, i+1, ret, alpha[k][i], like[k][i], coarse_grid[k]);
		}
	}

	/* No longer need nl_obj1; coarse search is done */
	nlopt_destroy(nl_obj1);

	/* Create new nlopt object for fine search */
	nlopt_opt nl_obj3 = nlopt_create(NLOPT_LN_COBYLA, 3);

	lb[0] = lb[1] = coarse_grid[0];
	ub[0] = ub[1] = coarse_grid[el_cg - 1];
	lb[2] = 0.0;
	ub[2] = 1.0;

	nlopt_set_lower_bounds(nl_obj3, lb);
	nlopt_set_upper_bounds(nl_obj3, ub);
	nlopt_set_maxeval(nl_obj3, 300);
	nlopt_set_min_objective(nl_obj3, (nlopt_func) par_3_opt, (void *) &nld);

#if 0
	// DEBUG TEST W MATLAB VALS
	x_fine[0] = 0.999878346327187;
	x_fine[1] = 0.990000000000000;
	x_fine[2] = 0.751467410761228;

	fine_val = par_3_opt(3, x_fine, NULL, (void *) &nld);	
	printf("DEBUG: %le\n", fine_val);
#endif

	/* Initial values for fine search */
	x_fine[0] = coarse_grid[gmax_pos[0]];
	x_fine[1] = coarse_grid[gmax_pos[1]];
	x_fine[2] = alpha[gmax_pos[0]][gmax_pos[1]];	
	printf("B4 XVAL: %lf %lf %lf %le\n", x_fine[0], x_fine[1], x_fine[2], -gmax_val); 

	ret |= nlopt_optimize(nl_obj3, x_fine, &fine_val);
	//printf("ret: %d XVAL: %lf %lf %lf FINE_VAL: %le\n", ret, x_fine[0], 
	//	x_fine[1], x_fine[2], fine_val);	
	
	/* Return a, b, alpha in pointers */
	*a_par = x_fine[0];
	*b_par = x_fine[1];
	*alph_par = x_fine[2];
	
	free(coarse_grid);

	if (ret < 0)
		return -1;
	else
		return 0;
}


/* Optimise with respect to alpha */
double alpha_opt(int n, const double *a, double *grad, void *nldv)
{
	/* FIXME ADD THE STUFF IN HERE TO DO WITH READING nld */
	int i;
	double alpha = *a;
	struct nl_extra_data *nld = (struct nl_extra_data *) nldv;	
	double sigma_tot = 0, *sigma = calloc(nld->len, sizeof(double));
	double like_a_tot = 0;
	double like_b_tot = 0, *like_b = calloc(nld->len, sizeof(double));

	for (i = 0; i < nld->len; i++) {
		/* Get sigma */
		sigma[i] = alpha * pow(nld->a_val, i) + 
			(1.0 - alpha) * pow(nld->b_val, i);

		sigma[i] = -1.0 / pow(sigma[i], 2);

		sigma[i] *= pow(nld->data_seg[i], 2);

		sigma_tot += sigma[i];
	}

	sigma_tot = sqrt(-sigma_tot / ( (double) nld->len ));

	for (i = 0; i < nld->len; i++) {
		/* Get likelihood*/
		like_a_tot += log( alpha * pow(nld->a_val, i) + 
			(1.0 - alpha) * pow(nld->b_val, i) );

		like_b[i] = alpha * pow(nld->a_val, i) + 
			(1.0 - alpha) * pow(nld->b_val, i);

		like_b[i] = pow(like_b[i], -2);

		like_b[i] *= pow(nld->data_seg[i], 2) / 
			(2.0 * pow(sigma_tot, 2));

		like_b_tot += like_b[i];	
	}

	like_b_tot = -like_a_tot - like_b_tot - 
		(double) nld->len * log(2 * M_PI * pow(sigma_tot, 2)) / 2.0;

	free(like_b);
	free(sigma);

	return -like_b_tot;
}

/* Optimise all three parameters */
double par_3_opt(int n, const double *a, double *grad, void *nldv)
{
	int i;
        struct nl_extra_data *nld = (struct nl_extra_data *) nldv;
	double a_val = a[0], b_val = a[1], alpha = a[2], sigma_tot = 0;
	double *env = calloc(nld->len, sizeof(double));
	double like = 0, x = 0;

	for (i = 0; i < nld->len; i++) {
		env[i] = alpha * pow(a_val, i) + (1.0 - alpha) * pow(b_val, i);

		sigma_tot -= (-1.0 / pow(env[i], 2)) * pow(nld->data_seg[i], 2);
	}

	sigma_tot = sqrt(sigma_tot / (double) nld->len);

	for (i = 0; i < nld->len; i++) {
		like -= log(env[i]);

		x += pow(env[i], -2) * pow(nld->data_seg[i], 2) / 
			(2.0 * pow(sigma_tot, 2));
	}

	like -= (x + (double) nld->len * log(2.0 * M_PI * 
			pow(sigma_tot, 2)) / 2.0);

	free(env);
	return -like;
}

/* Function uses Gnuplot to plot the wav envelope */
int plot_wav(float *wav_data, int channels, sf_count_t frames, int samprate)
{
	FILE *handle, *fp;
	int i;

	fp = fopen("foo.dat", "w");
	if (fp == NULL) {
		printf("Error opening file!\n");
		return -1;
	}

	
	fprintf(fp, "#Time\t\tAmplitude\n");
	for (i = 0; i < channels * frames; i++) {
	//for (i = 0; ((1/(double) samprate) * i) < 0.01; i++) {
		fprintf(fp, "%.10lf\t%f\n", (1/(double) samprate) * i, 
			wav_data[i]);
	}

	fclose(fp);

	handle = popen("gnuplot -persistent", "w");
	if (handle == NULL) {
		printf("Pipe open error\nExiting...");
		return -1;
	}

	fprintf(handle, "set term gif\n");
	fprintf(handle, "set output \"./wav.gif\"\n");
	fprintf(handle, "plot \"foo.dat\" with lines\n");	//test
	fprintf(handle, "set output\n");
	pclose(handle);

	return 0;
}
