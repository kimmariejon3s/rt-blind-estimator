/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
*****************************************************************************/

#include <stdio.h>
#include <sndfile.h>

int main(void)
{
	SF_INFO input_info;
	input_info.format = 0;

	SNDFILE *input = sf_open("/home/kim/wav_samples/bach_partita_e_maj.wav",
				SFM_READ, &input_info);

	printf("Sample Rate: %d\nChannels: %d\nFormat: %d\nSections: %d\nSeekable: %d\n", input_info.samplerate, input_info.channels, input_info.format, input_info.sections, input_info.seekable);

	return 0;
}
