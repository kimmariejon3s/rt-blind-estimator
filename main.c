/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sndfile.h>
#include <libswresample/swresample.h>
#include <libavutil/channel_layout.h>
#include <libavcodec/avcodec.h>
#include <liquid/liquid.h>
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

/* Functions */
int get_wav_data(void);
int process_wav_data(float *wav_data, SF_INFO input_info, SNDFILE *input);
int resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, 
	int samp_freq, const uint8_t **wav_data, int band);
int plot_wav(float *wav_data, int channels, sf_count_t frames, int samprate);

/* Globals */
int num_bands = 8;
int octave_bands[] = {63, 125, 250, 500, 1000, 2000, 4000, 8000};
int samp_freq_per_band[] = {3000, 3000, 3000, 3000, 3000, 6000, 12000, 24000};

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
	int i, ret = 0;
	float *resampled_wav;	

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
		ret = resamp_wav_data(resamp, input_info.samplerate, 
			(uint64_t) input_info.frames, samp_freq_per_band[i],
			(const uint8_t **) &wav_data, octave_bands[i]);

		if (ret < 0) {
			free(wav_data);
			sf_close(input);
			return -1;
		}

		/* Perform octave band filtering */

	}
}

int resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, int samp_freq, const uint8_t **wav_data, int band) {
	uint8_t *output_data;
	int ret, i, j, resamp_frames, half_len, flags;
	float *filtered, *x, *env, atten;

#if OCT_FILT_ORDER >= ENV_FILT_ORDER
	float b[3 * (L_OCT + R_OCT)] = {0};
	float a[3 * (L_OCT + R_OCT)] = {0};
#else
	float b[3 * (L_ENV + R_ENV)] = {0};
	float a[3 * (L_ENV + R_ENV)] = {0};
#endif

	resamp_frames = av_rescale_rnd(swr_get_delay(resamp, in_rate) + 
		num_frames, samp_freq, in_rate, AV_ROUND_UP);

	ret = av_samples_alloc(&output_data, NULL, 1, resamp_frames, 
		AV_SAMPLE_FMT_FLT, 0);

	if (ret < 0) {
		printf("Output buffer alloc failure\nExiting...\n");
		return -1;
	}

	resamp_frames = swr_convert(resamp, &output_data, resamp_frames,
		wav_data, num_frames);

	if (resamp_frames < 0) {
		printf("Resampling failure\nExiting...\n");
		av_freep(&output_data);
		return -1;
	}

	/* Bandpass filter - 7th and 8th parameter are ignored */
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_BANDPASS, 
		LIQUID_IIRDES_SOS, OCT_FILT_ORDER, 
		((float) band) / (sqrt(2.0) * ((float) samp_freq)), 
		((float) band) / ((float) samp_freq), 1.0, 1.0, b, a);
	//printf("res: %d\n\n", iirdes_isstable(b, a, 3 * (L_OCT + R_OCT)));

//	for (i=0; i<(3 * (L_OCT + R_OCT)); i++)
//		printf("b: %f and a: %f\n", b[i], a[i]);		

	iirfilt_rrrf f_obj_oct = iirfilt_rrrf_create_sos(b, a, (L_OCT + R_OCT));

	x = (float *) output_data;
	filtered = calloc(resamp_frames, sizeof(float)); 
//	x[0] = -1.331623101849397e-06;
//	x[1] = 6.919737431765877e-06;
//	x[2] = 5.191250643456192e-05;

	for (i = 0; i < resamp_frames; i++) {
		iirfilt_rrrf_execute(f_obj_oct, x[i], &filtered[i]);
	//	if (i < 100)
	//		printf("data%d: %lg %lg \n", i+1, (double) x[i], (double) filtered[i]);
	}

	iirfilt_rrrf_destroy(f_obj_oct);
	swr_free(&resamp);
	av_freep(&output_data);
	memset(b, 0, sizeof(b));
	memset(a, 0, sizeof(a));

#if 0
	/* Generate wav to check resampled data is not garbage */
	SF_INFO output_info;
	output_info.samplerate = samp_freq;
	output_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	output_info.channels = 1;

	SNDFILE* otpt = sf_open("/home/kim/wav_samples/Z.wav", 
		SFM_WRITE, &output_info);

	sf_write_float(otpt, (float *) output_data, resamp_frames);

//	FIXME: DO I want to return output_data or an int? Don't free it if I 
//		want to return it!
	return (float *)output_data;
#endif

	/* Low Pass Filter - 6th, 7th and 8th param are ignored */
	printf("DEBUG: Low pass filter and Hilbert transform\n");
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_LOWPASS, 
		LIQUID_IIRDES_SOS, ENV_FILT_ORDER, 
		((float) LOW_PASS_CUTOFF) / ((float)samp_freq), 0.1, 
		1.0, 1.0, b, a);

	iirfilt_rrrf f_obj_env = iirfilt_rrrf_create_sos(b, a, (L_ENV + R_ENV));

	//FIXME: do I want to memset x to zero? Tried, got seg fault
	env = calloc(resamp_frames, sizeof(float));

	flags = 0;
	
	float complex *filt_complex = calloc(resamp_frames, sizeof(float complex));
	float complex *tmp = calloc(resamp_frames, sizeof(float complex));
	float complex *tmp2 = calloc(resamp_frames, sizeof(float complex));
	float complex *hilb = calloc(resamp_frames, sizeof(float complex));

	fftplan fft_pln = fft_create_plan(resamp_frames, (float complex *) 
		filt_complex, (float complex *) tmp, LIQUID_FFT_FORWARD, flags);
	fftplan ifft_pln = fft_create_plan(resamp_frames, (float complex *) tmp,
 		(float complex *) hilb, LIQUID_FFT_BACKWARD, flags);

	fft_execute(ifft_pln);

	/* Get the Hilbert transform - really it's the analytic signal as the
		Hilbert transform is stored in the imaginary part of the soln */
	for (i = 0; i < resamp_frames; i+= 1)
		filt_complex[i] = (float complex)filtered[i];

	fft_execute(fft_pln);	

	for (i = 0; i < resamp_frames; i+= 1) {
		j = i + 1;
		tmp[i] *= (1.0);

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
	for (i = 0; i < 100; i++) {
		hilb[i] *= (1.0 / resamp_frames);
		printf("%d: %lg +i%lg\n\t%lg +i%lg\n", i, 
			creal(filt_complex[i]), cimag(filt_complex[i]), 
			creal(hilb[i]), cimag(hilb[i]));
	}

	//FIXME: why is this causing a fre pointer error???
	//fft_destroy_plan(fft_pln);
	free(tmp);
	free(tmp2);
	iirfilt_rrrf_destroy(f_obj_env);

	/* Get the full envelope */

	/* Apply low pass filter to envelope */
	for (i = 0; i < resamp_frames; i++) {
//		iirfilt_rrrf_execute(f_obj_env, x[i], &env[i]);
//		printf("data%d: %f %f \n", i, x[i], filtered[i]);
	}

	return 0;
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
