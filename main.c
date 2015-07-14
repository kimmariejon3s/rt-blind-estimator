/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>

/* Functions */
int get_wav_data(void);
int plot_wav();

int main(void)
{
	int ret;

	ret = get_wav_data();

	return ret;
}


int get_wav_data(void) {
	SF_INFO input_info;
	input_info.format = 0;
	short *wav_data;
	unsigned long long ret;	

	/* Open file */
	SNDFILE *input = 
		sf_open("/home/kim/wav_samples/mono_bach_partita_e_maj.wav",
		SFM_READ, &input_info);

	/* Print header data */
	printf("Frames: %llu\nSample Rate: %d\nChannels: %d\nFormat: 0x%X\n"
		"Sections: %d\nSeekable: %d\n",(unsigned long long) 
		input_info.frames, input_info.samplerate, input_info.channels, 
		input_info.format, input_info.sections, input_info.seekable);

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
	
	ret = plot_wav(wav_data, input_info.channels, input_info.frames);	

	free(wav_data);
	sf_close(input);
	return (int)ret;
}


int plot_wav(short **wav_data, int channels, sf_count_t frames)
{
	FILE *handle;

	handle = popen("gnuplot -persistent", "w");
	if (handle == NULL) {
		printf("Pipe open error\nExiting...");
		return -1;
	}

	fprintf(handle, "set term gif\n");
	fprintf(handle, "set output \"./wav.gif\"\n");
	fprintf(handle, "plot sin(x)/x\n");	//test
	fprintf(handle, "set output\n");
	pclose(handle);

	return 0;
}
