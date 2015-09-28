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
#define LOW_PASS_CUTOFF	80
#define DAMP_SEG_SIZE	0.5
#define DAMP_STEP_SIZE	0.002
#define MAX_RT		100.0	
#define MIN_RT		0.005	
#define	BEST_S_NUM	100000
#define PAR_UP_BOUND	0.9999
#define PAR_LOW_BOUND	0.99
#define SQP_STEP	5

/* Functions */
int get_wav_data(void);
int process_wav_data(float *wav_data, SF_INFO input_info, SNDFILE *input);
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
int * perform_ml(int **start_end, float *env, int s_e_size, 
	float *filtered_wav, int resamp_frames, int samp_freq);
int ml_fit(float *data_seg, int len, float samp_freq, float max_abs_seg);
float alpha_opt(float alpha, float *data_seg, int len, float a_val, 
	float b_val);
int plot_wav(float *wav_data, int channels, sf_count_t frames, int samprate);

/* Globals */
int num_bands = 8;
int octave_bands[] = {63, 125, 250, 500, 1000, 2000, 4000, 8000};
int samp_freq_per_band[] = {3000, 3000, 3000, 3000, 3000, 6000, 12000, 24000};

#if OCT_FILT_ORDER >= ENV_FILT_ORDER
	float b[3 * (L_OCT + R_OCT)] = {0};
	float a[3 * (L_OCT + R_OCT)] = {0};
#else
	float b[3 * (L_ENV + R_ENV)] = {0};
	float a[3 * (L_ENV + R_ENV)] = {0};
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
	int ret;

	ret = get_wav_data();

	return ret;
}


int get_wav_data(void) {
	SF_INFO input_info;
	float *wav_data;
	unsigned long long long_ret = 0;	
	int ret;

	/* Open file */
	input_info.format = 0;
	SNDFILE *input = 
		sf_open("/home/kim/wav_samples/mono_bach_partita_e_maj.wav",
		SFM_READ, &input_info);

	/* Print header data */
	printf("Frames: %llu\nSample Rate: %d\nChannels: %d\nFormat: 0x%X\n"
		"Sections: %d\nSeekable: %d\n",(unsigned long long) 
		input_info.frames, input_info.samplerate, input_info.channels, 
		input_info.format, input_info.sections, input_info.seekable);

	/* FIXME: insert check to make sure that the wav file is in the 
	 *	expected format: mono and 16 bit signed??? */

	/* Store data in array */
	wav_data = (float *) calloc(input_info.frames * input_info.channels, 
		sizeof(float));

	if (wav_data == NULL) {
		printf("Memory allocation error\nExiting...\n");
		sf_close(input);
		return -1;
	}

	/* Single channel, so can call this way */
	long_ret = sf_read_float(input, wav_data, input_info.frames);		

	if (long_ret != input_info.frames) {
		printf("Error! Some frames not read\nExiting...");
		free(wav_data);
		sf_close(input);
		return -1;
	}


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
	ret = process_wav_data(wav_data, input_info, input);

	/* Clean up */	
	free(wav_data);
	sf_close(input);
	return ret;
}


int process_wav_data(float *wav_data, SF_INFO input_info, SNDFILE *input) {
	int i, ret, s_e_size, resamp_frames = 0, **start_end;
	float *resampled_wav, *filtered_wav, *env;

	printf("DEBUG: resample init\n");

	/* Specifications that apply to all subbands */
	SwrContext *resamp = swr_alloc();

	av_opt_set_int(resamp, "in_channel_layout", AV_CH_LAYOUT_MONO, 0);
	av_opt_set_int(resamp, "out_channel_layout", AV_CH_LAYOUT_MONO, 0);
	av_opt_set_int(resamp, "in_sample_rate", input_info.samplerate, 0);
	av_opt_set_sample_fmt(resamp, "in_sample_fmt", AV_SAMPLE_FMT_FLT, 0);
	av_opt_set_sample_fmt(resamp, "out_sample_fmt", AV_SAMPLE_FMT_FLT, 0);


	/* Loop through each octave band */
	for (i = 0; i < 1; i++) {
		av_opt_set_int(resamp, "out_sample_rate", 
			samp_freq_per_band[i], 0); 

		ret = swr_init(resamp);

		if (ret < 0) {
			printf("Resample initialisation error\nExiting...\n");
			free(wav_data);
			sf_close(input);
			return -1;
		}

		/* Resample the audio */
		resampled_wav = resamp_wav_data(resamp, input_info.samplerate, 
			(uint64_t) input_info.frames, samp_freq_per_band[i],
			(const uint8_t **) &wav_data, octave_bands[i], 
			&resamp_frames);

		if (resamp_frames <= 0 || resampled_wav == NULL) {
			free(wav_data);
			sf_close(input);
			return -1;
		}

		/* Octave band filtering */
		filtered_wav = oct_filt_data(resampled_wav, 
			(float) octave_bands[i], (float) samp_freq_per_band[i],
			resamp_frames);
		
		/* Obtain signal envelope */
		env = get_envelope((float) samp_freq_per_band[i], filtered_wav, 
			resamp_frames);

		/* Polyfit algorithm & get decay segments */
		start_end = apply_polyfit((float) samp_freq_per_band[i], 
			resamp_frames, env, &s_e_size);

		/* Preform Maximum Liklihood Estimation */
		perform_ml(start_end, env, s_e_size, filtered_wav, 
			resamp_frames, samp_freq_per_band[i]);
	}
	return 0;
}

float * resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, int samp_freq, const uint8_t **wav_data, int band, int *resamp_frms) {
	uint8_t *output_data;
	int ret, resamp_frames;
	float *x;

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

	/* Create bandpass filter - 7th and 8th parameter are ignored */
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_BANDPASS, 
		LIQUID_IIRDES_SOS, OCT_FILT_ORDER, 
		band / (sqrt(2.0) * samp_freq), 
		band / samp_freq, 1.0, 1.0, b, a);
	//printf("res: %d\n\n", iirdes_isstable(b, a, 3 * (L_OCT + R_OCT)));

//	for (i=0; i<(3 * (L_OCT + R_OCT)); i++)
//		printf("b: %f and a: %f\n", b[i], a[i]);		

	iirfilt_rrrf f_obj_oct = iirfilt_rrrf_create_sos(b, a, (L_OCT + R_OCT));


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

	memset(b, 0, sizeof(b));
	memset(a, 0, sizeof(a));

	/* Create Low Pass Filter - 6th, 7th and 8th param are ignored */
	printf("DEBUG: Low pass filter and Hilbert transform\n");
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_LOWPASS, 
		LIQUID_IIRDES_SOS, ENV_FILT_ORDER, 
		((float) LOW_PASS_CUTOFF) / samp_freq, 0.1, 
		1.0, 1.0, b, a);

	iirfilt_crcf f_obj_env = iirfilt_crcf_create_sos(b, a, (L_ENV + R_ENV));

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
	float *store, *var;
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
	store = calloc(num_segs, sizeof(float));
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

		store[i] = 6.91 / (resamp_frames * 
			logf(powf(10.0, -poly_res[i][0])));

		for (j = 0; j < poly_coeff; j++) {
			//FIXME eeeh floating point comparisons are DODGEY in C
			if (poly_res[i][j] < min_poly_res || 
				poly_res[i][j] > max_poly_res) {
				invalid_coeff++;
			}
		}
		if (invalid_coeff  == poly_coeff) {
			store[i] = 0.0;
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
	free(store);
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

int * perform_ml(int **start_end, float *env, int s_e_size, float *filtered_wav, int resamp_frames, int samp_freq) {
	int sz, i, j, k, len, min_abs_filt_seg, max_abs_filt_seg, max2;
	float *segment, *env_seg, *abs_filt_seg, *seg_seg2;
	float *abs_filtered = calloc(resamp_frames, sizeof(float));
	int *store_start = calloc(s_e_size, sizeof(int)); 
	int *store_end = calloc(s_e_size, sizeof(int));

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
		len = start_end[1][i] - start_end[0][i];

		/* Store decay segment in new array */
		for (k = 0, j = start_end[0][i]; j < start_end[1][i]; j++, k++)
		{
			segment[k] = filtered_wav[j];
			env_seg[k] = env[j];
			abs_filt_seg[k] = abs_filtered[j];
		}

		/* Get location of minimum */
		min_abs_filt_seg = len - roundf((DAMP_SEG_SIZE / 2) * 
			samp_freq) - 1;
		for (j = min_abs_filt_seg; j < len; j++) {
			if (abs_filt_seg[j] < abs_filt_seg[min_abs_filt_seg])
				min_abs_filt_seg = j;
		}

		/* Get location of maximum */
		max_abs_filt_seg = 0;
		for (j = 0; j < roundf((DAMP_SEG_SIZE / 2) * samp_freq); j++) {
			if (abs_filt_seg[j] > abs_filt_seg[max_abs_filt_seg])
				max_abs_filt_seg = j;
		}

		/* Store new fine-tuned decay segment */
		for (k = 0, j = max_abs_filt_seg; j <= min_abs_filt_seg; k++, 
			j++) {
			seg_seg2[k] = segment[j];
		}

		max2 = k;

		k = 0;
		for (j = 0; j <= max2; j++) {
			if (fabs(seg_seg2[j]) > fabs(seg_seg2[k]))
				k = j; 
		}
			
		for (j= 0; j <= max2; j++)
			seg_seg2[j] /= fabs(seg_seg2[k]);

		/* Save fine-tuned start and end locations */
		store_start[i] = start_end[0][i] + max_abs_filt_seg;
		store_end[i] = start_end[0][i] + max_abs_filt_seg + max2 - 1;

		/* ML fitting of decay model to the data */
		ml_fit(seg_seg2, max2, (float) samp_freq, fabs(seg_seg2[k]));

		/* Reset to zero b/c arrays are longer than they need to be */
		memset(segment, 0, sz * sizeof(float));
		memset(seg_seg2, 0, (sz + sz - roundf((DAMP_SEG_SIZE / 2) * 
			samp_freq)) * sizeof(float));
		memset(env_seg, 0, sz * sizeof(float));
		memset(abs_filt_seg, 0, sz * sizeof(float));
	}

	free(segment);
	free(seg_seg2);
	free(env_seg);
	free(abs_filt_seg);
	free(abs_filtered);
	free(store_start);
	free(store_end);
	return NULL;
}

int ml_fit(float *data_seg, int len, float samp_freq, float max_abs_seg) {
	float min, max, interval, j, alpha, *coarse_grid;
	nlopt_opt nl_obj1 = nlopt_create(NLOPT_LD_SLSQP, 3);
	int i, k;

	min = -6.91 / logf(PAR_LOW_BOUND) / 3000.0;
	min = exp(-6.91 / (samp_freq * min));

	max = -6.91 / logf(PAR_UP_BOUND) / 3000.0;	
	max = exp(-6.91 / (samp_freq * max));

	interval = (max - min) / ((float) SQP_STEP - 1.0);
	coarse_grid = calloc(max / interval, sizeof(float));

	/* Fill in coarse grid */
	for (i = 0, j = 0; i < max / interval; i++, j += interval)
		coarse_grid[i] = j;

	//FIXME: There is another max(abs(x)) division here in MTLB - typo?
	for (i = 0; i < len; i++)
		data_seg[i] /= max_abs_seg;

	for (i = 0; i < SQP_STEP; i++) {
		for (k = 0; k < SQP_STEP; k++) {
			alpha = 0.5;

			alpha_opt(alpha, data_seg, len, coarse_grid[k], 
				coarse_grid[i]);		
		}
	}


	free(coarse_grid);
	return 0;
}


/* Optimise with respect to alpha */
float alpha_opt(float alpha, float *data_seg, int len, float a_val, float b_val)
{
	int i;
	float sigma_tot = 0, *sigma = calloc(len, sizeof(float));
	float like_a_tot = 0;
	float like_b_tot = 0, *like_b = calloc(len, sizeof(float));

	for (i = 0; i < len; i++) {
		/* Get sigma */
		sigma[i] = alpha * powf(a_val, i) + 
			(1.0 - alpha) * powf(b_val, i);

		sigma[i] = -1.0 / powf(sigma[i], 2);

		sigma[i] *= powf(data_seg[i], 2);

		sigma_tot += sigma[i];

		/* Get likelihood*/
		like_a_tot += logf( alpha * powf(a_val, i) + 
			(1.0 - alpha) * powf(b_val, i) );

		like_b[i] = alpha * powf(a_val, i) + 
			(1.0 - alpha) * powf(b_val, i);

		like_b[i] = powf(like_b[i], -2);

		like_b[i] *= powf(data_seg[i], 2) / (2.0 * powf(sigma_tot, 2));

		like_b_tot += like_b[i];	
	}

	sigma_tot = sqrtf(-sigma_tot / ( (float) len ));

	like_b_tot = -like_a_tot - like_b_tot - 
		(float) len * logf(2 * M_PI * powf(sigma_tot, 2)) / 2.0;

	free(like_b);
	free(sigma);

	return -like_b_tot;
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
