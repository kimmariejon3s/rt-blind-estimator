/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include "kiss_fft.h"
#include "liquid.h"
#include <libswresample/swresample.h>
#include <libavutil/channel_layout.h>

#define MAX_CHUNK_SIZE	1024 * 1024 * 1024
#define WINDOW_WIDTH	0.5
#define OVERLAP		0.98
#define MAX_SAMP_FREQ	24000	/* Max freq in samp_freq_per_band[] */

//FIXME change this
#define CHAN	1

/* Globals */
int num_bands = 8;
int octave_bands[] = {63, 125, 250, 500, 1000, 2000, 4000, 8000};
unsigned int oct_filt_order = 5;
int samp_freq_per_band[] = {3000, 3000, 3000, 3000, 3000, 6000, 12000, 24000};
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

/* Functions */
int get_wav_data(void);
int plot_wav(short *wav_data, int channels, sf_count_t frames, int samprate);

//float test_b[11];
//float test_a[11];

int main(void)
{
	int ret;

	ret = get_wav_data();

	return ret;
}


int get_wav_data(void) {
	SF_INFO input_info;
	short *wav_data, *oct_filt_data;
	unsigned long long ret = 0;	
	uint8_t *output_data;
	int i, ret2, resamp_frames;

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
	wav_data = (short *) calloc(input_info.frames * input_info.channels, 
		sizeof(short));

	if (wav_data == NULL) {
		printf("Memory allocation error\nExiting...\n");
		sf_close(input);
		return -1;
	}

	ret = sf_readf_short(input, wav_data, input_info.frames);		

	if (ret != input_info.frames) {
		printf("Error! Some frames not read\nExiting...");
		free(wav_data);
		sf_close(input);
		return -1;
	}

	printf("DEBUG: resample\n");
	/* Specifications that apply to all subbands */
	SwrContext *resamp = swr_alloc();

	av_opt_set_int(resamp, "in_channel_layout", AV_CH_LAYOUT_MONO, 0);
	av_opt_set_int(resamp, "out_channel_layout", AV_CH_LAYOUT_MONO, 0);
	av_opt_set_int(resamp, "in_sample_rate", input_info.samplerate, 0);
	av_opt_set_sample_fmt(resamp, "in_sample_fmt", AV_SAMPLE_FMT_S16, 0);
	av_opt_set_sample_fmt(resamp, "out_sample_fmt", AV_SAMPLE_FMT_S16, 0);

	/* Octave band filtering */
	for (i = 0; i < num_bands; i++) {
		av_opt_set_int(resamp, "out_sample_rate", 
			samp_freq_per_band[i], 0); 

		ret2 = swr_init(resamp);

		if (ret2 < 0) {
			printf("Resample initialisation error\nExiting...\n");
			free(wav_data);
			sf_close(input);
                	return -1;
		}

		/* Resample the audio */
		resamp_frames = av_rescale_rnd(swr_get_delay(resamp, 
			input_info.samplerate) + input_info.frames, 
			samp_freq_per_band[i], input_info.samplerate, 
			AV_ROUND_UP);

		ret2 = av_samples_alloc(&output_data, NULL, 1, resamp_frames, 
			AV_SAMPLE_FMT_S16, 0);

		if (ret2 < 0) {
			printf("Output buffer alloc failure\nExiting...\n");
			free(wav_data);
			sf_close(input);
                	return -1;
		}

		resamp_frames = swr_convert(resamp, &output_data, resamp_frames,
			(const uint8_t **) &wav_data, input_info.frames);

		if (resamp_frames < 0) {
			printf("Resampling failure\nExiting...\n");
			av_freep(&output_data);
			free(wav_data);
			sf_close(input);
                	return -1;
		}

#if 0
		/* Generate wav to check resampled data is not garbage */
		SF_INFO output_info;
		output_info.samplerate = samp_freq_per_band[i];
		output_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
		output_info.channels = 1;

		SNDFILE* otpt = sf_open("/home/kim/wav_samples/Z.wav", 
			SFM_WRITE, &output_info);

		sf_writef_short(otpt, (short *) output_data, resamp_frames);

		/* Filter coefficient calculation using Liquid */
		liquid_iirdes(LIQUID_IIRDES_BUTTER,
			LIQUID_IIRDES_BANDPASS, LIQUID_IIRDES_TF, 
			oct_filt_order, 
			octave_bands[i] / (sqrt(2) * samp_freq_per_band[i]),
			octave_bands[i] / samp_freq_per_band[i], 1.0, 1.0,
			test_b, test_a);		
#endif		
		/* Octave band filtering */
		float *tmp = calloc(resamp_frames, 4);

		if (tmp ==NULL) {
			printf("Output buffer alloc failure\nExiting...\n");
			av_freep(&oct_filt_data);
			av_freep(&output_data);
			free(wav_data);
			sf_close(input);
                	return -1;
		}

		iirfilt_rrrf filt_obj = iirfilt_rrrf_create(&butter_b[i][0], 
			2 * oct_filt_order + 1, &butter_a[i][0], 
			2 * oct_filt_order + 1);	


		int j;
		for (j = 0; j < resamp_frames; j++)
			iirfilt_rrrf_execute(filt_obj, (float) output_data[j], 
				&tmp[j]);

		iirfilt_rrrf_destroy(filt_obj);

		av_freep(&output_data);
	}
	

	/* Hilbert Transform */
	printf("DEBUG: Hilbert transform\n");

	

#if 0	
	ret = plot_wav(wav_data, input_info.channels, input_info.frames, 
		input_info.samplerate);	

	for (i = 0, k = 0; i < WINDOW_WIDTH; i += WINDOW_WIDTH * (1 - OVERLAP),
	k++) {
		for (j = 0; j < input_info.frames * input_info.channels ; j++) {
			wav_data[j];
		}
	}

#endif
	free(wav_data);
	sf_close(input);
	return (int)ret;
}


/* Function uses Gnuplot to plot the wav envelope */
int plot_wav(short *wav_data, int channels, sf_count_t frames, int samprate)
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
		fprintf(fp, "%.10lf\t%hd\n", (1/(double) samprate) * i, 
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
