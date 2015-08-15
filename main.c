/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include "kiss_fft.h"
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
int samp_freq_per_band[] = {3000, 3000, 3000, 3000, 3000, 6000, 12000, 24000};

/* Functions */
int get_wav_data(void);
int plot_wav(short *wav_data, int channels, sf_count_t frames, int samprate);

int main(void)
{
	int ret;

	ret = get_wav_data();

	return ret;
}


int get_wav_data(void) {
	SF_INFO input_info;
	short *wav_data;
	unsigned long long ret = 0;	
	int i, ret2;

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
	for (i = 0; i < 1; i++) {
		av_opt_set_int(resamp, "out_sample_rate", 
			samp_freq_per_band[i], 0); 

		ret2 = swr_init(resamp);

		if (ret2 < 0) {
			printf("Resample initialisation error\nExiting...\n");
			free(wav_data);
			sf_close(input);
                	return -1;
		}

		uint8_t *output;
		int out_samples = av_rescale_rnd(swr_get_delay(resamp, 
			input_info.samplerate) + input_info.frames, 
			samp_freq_per_band[i], input_info.samplerate, 
			AV_ROUND_UP);

		av_samples_alloc(&output, NULL, 1, out_samples, 
			AV_SAMPLE_FMT_S16, 0);
		out_samples = swr_convert(resamp, &output, out_samples,
			(const uint8_t **) &wav_data, input_info.frames);

		SF_INFO output_info;
		output_info.samplerate = samp_freq_per_band[i];
		output_info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
		output_info.channels = 1;

		SNDFILE* otpt = sf_open("/home/kim/wav_samples/Z.wav", 
			SFM_WRITE, &output_info);


		sf_writef_short(otpt, (short *) output, out_samples);
		


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
