/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include "kiss_fft.h"
#include <libswresample/swresample.h>
#include <libavutil/channel_layout.h>
#include <libavcodec/avcodec.h>
#include <liquid/liquid.h>

#define MAX_CHUNK_SIZE	1024 * 1024 * 1024
#define WINDOW_WIDTH	0.5
#define OVERLAP		0.98
#define OCT_FILT_ORDER	5

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

#if 0
/* Parameters as calculated by matlab - currently not used */
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
	/* Hilbert Transform */
	printf("DEBUG: Hilbert transform\n");
}

int resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, int samp_freq, const uint8_t **wav_data, int band) {
	uint8_t *output_data;
	int ret, resamp_frames;
	float b[2 * OCT_FILT_ORDER + 1] = {0};
	float a[2 * OCT_FILT_ORDER + 1] = {0};

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

	int i;
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_BANDPASS, 
		LIQUID_IIRDES_TF, OCT_FILT_ORDER, 
		((float) band) / (sqrt(2.0) * ((float)samp_freq)),
		((float) band) / ((float)samp_freq), 1.0, 1.0, b, a);

	iirfilt_rrrf f_obj = iirfilt_rrrf_create(b, 2 *
		OCT_FILT_ORDER + 1, a, 2 * OCT_FILT_ORDER + 1);

	for (i = 0; i < 11; i++)
		printf("b: %f a: %f\n", b[i], a[i]);
	float y;
	float *x;
	x = (float *) output_data;

	for (i = 0; i < 100; i++) {
		iirfilt_rrrf_execute(f_obj, x[i], &y);
		printf("data%d: %f %f \n", i, x[i], y);
	}

	iirfilt_rrrf_destroy(f_obj);
	swr_free(&resamp);
	av_freep(&output_data);


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

	return 0;
}



#if 0
int octave_filter_wav() {
	int i;

	iirfilt_rrrf f_obj = iirfilt_rrrf_create(butter_b[0], 2 * 
		OCT_FILT_ORDER + 1, butter_a[0], 2 * OCT_FILT_ORDER + 1);
	float y[100];

	for (i = 0; i < 100; i++) {
		iirfilt_rrrf_execute(f_obj, (float) output_data[i], &y[i]);
		printf("data%d: %f %f\n", i, (float)output_data[i], y[i]);
	}
	av_freep(&output_data);
}
#endif
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
