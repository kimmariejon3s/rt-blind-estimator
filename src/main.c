/*****************************************************************************
* This program allows the reverberation time of a venue to be determined
* blindly, using a sample of music or speech that is recorded in the venue
* C implementation of Paul Kendrick's Matlab algorithm:
*
* Kendrick, Paul, et al. "Blind estimation of reverberation parameters for
* non-diffuse rooms." Acta Acustica united with Acustica 93.5 (2007): 760-770.
*
* C implementation by Kim Jones <K.Jones20@edu.salford.ac.uk>
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <sndfile.h>
#include <libswresample/swresample.h>
#include <libavutil/channel_layout.h>
#include <libavcodec/avcodec.h>
#include <liquid.h>
#include <nlopt.h>
#include <svdlib.h>
#include "butter_params.h"
#include <libavutil/opt.h>
#include <windows.h>
#include <omp.h>

#define	MAX_PATH_SZ	512 
#define SPLIT		20		/* Process wav in SPLIT segment(s) */
#define MAX_CHUNK_SIZE	1024 * 1024 * 1024
#define WINDOW_WIDTH	0.5
#define OVERLAP		0.98
#define USE_HARD_BUTTER	1		/* Use butter_params.h if set */
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
#define	BEST_S_NUM	20000		/* Set to 100000 in Matlab */
#define PAR_UP_BOUND	0.99999
#define PAR_LOW_BOUND	0.99
#define SQP_STEP	5
//#define COARSE_DR_EST	1
#define SEG_LEN		3 * 60		/* 3 minutes */
#define RT_DECAY	25		/* 25 = T20 decay; 35 = T30 decay */
#define DISCARD_THRESH	0		/* Discard decay segs if <= this */

/* Functions */
int get_wav_data(char *filename, char *cwd);
int compute_rt(int samp, int array_size, int *store_start, int *store_end, 
	double *dr, double *a, double *b, double *alpha, double *mean_rt,
	double *rt_sd, int band, FILE *fh);
int optimum_model(double a, double b, double alpha, double dr, int samp, int n,
	double *chan_2);
int optimum_model_2(double **chan, int n, int chan_sz, int samp, double *lin);
void hanning(int len, double *hann);
int process_wav_data(int band, float *wav_data, SF_INFO input_info, 
	SNDFILE *input, unsigned long long frames, double *p_alpha,
	double *p_a, double *p_b, double *p_dr, int *p_start, int *p_end, int
	*filt_frames);
float * resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, 
	int samp_freq, const uint8_t **wav_data, int *resamp_frms);
float * oct_filt_data(float *output_data, int band_num, float samp_freq,
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
	int **p_start, int **p_end);
int ml_fit(float *data_seg, int len, float samp_freq, double *a_par, 
	double *b_par, double *alph_par);
double get_decay_region(int c_len, double a, double b, double alpha,
		int len_store);
double alpha_opt(int n, const double *a, double *grad, void *nldv);
double par_3_opt(int n, const double *a, double *grad, void *nldv);

/* Structs */
struct nl_extra_data {
	float *data_seg;
	int len;
	double a_val;
	double b_val;
};


/* Globals */
int seg_len_val = SEG_LEN;
int num_bands = 8;
int octave_bands[] = {63, 125, 250, 500, 1000, 2000, 4000, 8000};
int samp_freq_per_band[] = {3000, 3000, 3000, 3000, 3000, 6000, 12000, 24000};
const char fname[] = "blind_rts.log";

#if OCT_FILT_ORDER >= ENV_FILT_ORDER
	float butt_b[3 * (L_OCT + R_OCT)] = {0};
	float butt_a[3 * (L_OCT + R_OCT)] = {0};
#else
	float butt_b[3 * (L_ENV + R_ENV)] = {0};
	float butt_a[3 * (L_ENV + R_ENV)] = {0};
#endif

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine,
		int nCmdShow) 
{
	double ret;
	char filename[MAX_PATH_SZ] = "";
	char tmp_str[MAX_PATH_SZ]; 
	char cwd[MAX_PATH_SZ]; 
	int win_ret;
	char win_read[2];

	if(getcwd(cwd, sizeof(cwd)) == NULL)
		cwd[0] = '\0';

	OPENFILENAME ofn;
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = NULL;
	ofn.lpstrFilter = "WAV Files (*.wav)\0*.wav\0";
	ofn.lpstrFile = filename;
	ofn.nMaxFile = MAX_PATH_SZ;
	ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST;
	ofn.lpstrDefExt = "wav";

	snprintf(tmp_str, MAX_PATH_SZ, "Welcome to the Blind Reverberation Time Estimation "
		"Program!\n\nOutput for this program will be printed in a console window.\n"
		"\nData will also be logged to the file %s (in the same directory as the "
		"application).\n\nClick OK to continue...", fname);

	MessageBox(NULL, tmp_str, 
		"Blind Reverberation Time Estimation", MB_OK | MB_ICONASTERISK);
	MessageBox(NULL, "Please select a mono, 16-bit signed wav file when prompted.\n"
		"A file dialog box will open after you click OK...", 
		"Blind Reverberation Time Estimation", MB_OK | MB_ICONASTERISK);

	if (GetOpenFileName(&ofn)) {
		win_ret = MessageBox(NULL,
			"After the decay regions have been identified, the wav "
			"file will be analysed in segments.\n\nIdeally each segment should "
			"have 36 decays of at least -25 dB.\n\nThe default segment length "
			"is 3 minutes.\n\nDo you want to change the default segment "
			"length?\n\nIf you select \"No\", the default value will be used."
			"\n\nIf you are unsure what segment length to use, select \"No\" "
			"to keep the default value.\n\n"
			"If you select \"Yes\", you will be prompted to enter the desired "
			"segment length (in minutes) in the console window. Please press "
			"the enter key once you have typed the integer.",
			"Blind Reverberation Time Estimation", MB_YESNO | MB_ICONASTERISK);

		if (win_ret == IDYES) {
			printf("Please enter a value for segment length (in minutes).\n"
				"Input should be an integer between 1 and 9:\n");
			fgets(win_read, 2, stdin);
			if ((int) win_read[0] - (int) '0' > 0 &&
					(int) win_read[0] - (int) '0' < 10) {
				seg_len_val = (int) win_read[0] - (int) '0';

				snprintf(tmp_str, MAX_PATH_SZ, "Segment length is now "	
					"%d minute(s)!\nClick OK to continue...",
					seg_len_val);

				MessageBox(NULL, tmp_str,
					"Blind Reverberation Time Estimation", MB_OK | 
					MB_ICONASTERISK);	
			} else {
				 MessageBox(NULL, "Error! Invalid input for segment "
					"length!\n\nDefault of 3 minutes will be used...",
					"Blind Reverberation Time Estimation", MB_OK |
					MB_ICONEXCLAMATION);	
			}
		}
		chdir(cwd);
		ret = get_wav_data(filename, cwd);
	} else { 
		MessageBox(NULL, "Error Selecting File! Exiting...",
			"Blind Reverberation Time Estimation", MB_OK | MB_ICONEXCLAMATION);
		return 0;
	}

	if (ret < 0) {
		MessageBox(NULL, "Error returned during calculations...\nCheck console for"
			" more information.\nConsole will"
			" close after you click OK.",
			"Blind Reverberation Time Estimation", MB_OK | MB_ICONEXCLAMATION);
		return -1;
	} else {
		snprintf(tmp_str, MAX_PATH_SZ, "Results exported to log file: %s\n"
			"Console will close after you click OK.", fname);
		MessageBox(NULL, tmp_str, 
			"Blind Reverberation Time Estimation", MB_OK | MB_ICONASTERISK); 
        }
	return 0;
}


int get_wav_data(char *filename, char *cwd) {
	SF_INFO input_info;
	float *wav_data;
	unsigned long long long_ret = 0;	
	int *store_start, *store_end;
	int i, j, k, band, array_size, current_sz = 0, ret = 0;
	int filt_frm_tot, filt_frames;
	double *a, *b, *alpha, *dr, mean_rt, rt_sd;
	time_t ct = time(NULL);
	FILE *fh = fopen(fname, "a+");

	if (fh == NULL) {
		printf("Error opening log file!\n");
		return -1;
	}

	/* Open file */
	input_info.format = 0;
	SNDFILE *input = sf_open(filename, SFM_READ, &input_info);

	if (input == NULL) {
		MessageBox(NULL, "Error opening Wav file!",
			"Blind Reverberation Time Estimation", MB_ICONEXCLAMATION | MB_OK);
		return -1;
	}

	/* Print header data */
	printf("Frames: %llu\nSample Rate: %d\nChannels: %d\nFormat: 0x%X\n"
		"Sections: %d\nSeekable: %d\n",(unsigned long long) 
		input_info.frames, input_info.samplerate, input_info.channels, 
		input_info.format, input_info.sections, input_info.seekable);

	/* Check to make sure that the wav file is in the expected format:
		mono channel and 16 bit signed wav */
	if ((input_info.format & 0xFF000F) != 0x10002 &&
			(input_info.format & 0xFF000F) != 0x130002) {
		MessageBox(NULL, "Program requires 16-bit signed wav.\nMake sure format "
			"is 0x10002 or 0x130002.\nExiting...\n", 
			"Blind Reverberation Time Estimation", MB_ICONEXCLAMATION | MB_OK);
		sf_close(input);
		return -1;
	}
	if (input_info.channels != 1) {
		MessageBox(NULL, "Input must be a single channel wav file. Exiting...\n",
			"Blind Reverberation Time Estimation", MB_ICONEXCLAMATION | MB_OK);
		sf_close(input);
		return -1;
	}
	
	printf("Reading %s\n", filename);

	fprintf(fh, "=====================================================\n");
	fprintf(fh, "Reverberation Time calculations: %s\n", filename);
	fprintf(fh, "Time and Date: %s\n", ctime(&ct));

	/* Create arrays large enough to store values */
	a = calloc(BEST_S_NUM, sizeof(double));
	b = calloc(BEST_S_NUM, sizeof(double));
	alpha = calloc(BEST_S_NUM, sizeof(double));
	dr = calloc(BEST_S_NUM, sizeof(double));
	store_start = calloc(BEST_S_NUM, sizeof(int));
	store_end = calloc(BEST_S_NUM, sizeof(int));

	/* Store data in array */
	wav_data = (float *) calloc(input_info.frames / SPLIT, sizeof(float));

	if (wav_data == NULL) {
		printf("Memory allocation error\nExiting...\n");
		sf_close(input);
		return -1;
	}

	/* Octave band filtering */
	for (band = 0; band < 8; band++) {
		printf("\n\nCalculations of Reverberation Time for %d Hz "
			"octave band.\n", octave_bands[band]);
		fprintf(fh, "\nCalculations of Reverberation Time for %d Hz "
			"octave band.\n", octave_bands[band]);

		array_size = 0;
		current_sz = 0;
		filt_frm_tot = 0;
		filt_frames = 0;

		/* Single channel, so can call this way */
		for (i = 0; i < SPLIT; i++) {
			printf("\nPart %d of %d\n", i + 1, SPLIT);

			long_ret = sf_read_float(input, wav_data, 
					floor(input_info.frames / SPLIT));

			/* Process the wav data */
			ret = process_wav_data(band, wav_data, input_info,
				input, long_ret, alpha, a, b,
				dr, store_start, store_end, &filt_frames);

			/* If ret < 0, error */
			/* If >= 0, ret is the size of the returned arrays */
			if (ret < 0) {
				free(wav_data);
				sf_close(input);
				return -1;
			} else {
				array_size += ret;
				current_sz = ret;
				filt_frm_tot += filt_frames;
			}

			for (k = 0; k < current_sz; k++) {
				store_start[k] += (filt_frm_tot - filt_frames);
				store_end[k] += (filt_frm_tot - filt_frames);
			}

			a += current_sz;
			b += current_sz;
			alpha += current_sz;
			dr += current_sz;
			store_start += current_sz;
			store_end += current_sz;


			/* Reset bytes to zero */
			memset(wav_data, 0, sizeof(float) *
				(input_info.frames / SPLIT));
		}

		/* Return the pointers to zero position */
		if (sf_seek(input, 0, SEEK_SET) < 0) {
			printf("Error returning to start of wav!\n"
				"Exiting...\n");
			free(wav_data);
			sf_close(input);
			return ret;
		}

		a -= array_size;
		b -= array_size;
		alpha -= array_size;
		dr -= array_size;
		store_start -= array_size;
		store_end -= array_size;

		/* Compute the reverberation time (RT) */
		ret = compute_rt(samp_freq_per_band[band],
			array_size, store_start, store_end, dr, a, b,
			alpha, &mean_rt, &rt_sd, octave_bands[band], fh);

		fprintf(fh, "=============================================\n");

		/* Reset bytes to zero */
		memset(a, 0, sizeof(double) * BEST_S_NUM);
		memset(b, 0, sizeof(double) * BEST_S_NUM);
		memset(alpha, 0, sizeof(double) * BEST_S_NUM);
		memset(dr, 0, sizeof(double) * BEST_S_NUM);
		memset(store_start, 0, sizeof(int) * BEST_S_NUM);
		memset(store_end, 0, sizeof(int) * BEST_S_NUM);
	}


	/* Clean up */
	sf_close(input);
	free(wav_data);
	free(a);
	free(b);
	free(alpha);
	free(dr);
	free(store_start);
	free(store_end);
	fclose(fh);

	return ret;
}


/* 
 * Compute RT for room 
 */
int compute_rt(int samp, int array_size, int *store_start, int *store_end, 
		double *dr, double *a, double *b, double *alpha, 
		double *mean_rt, double *rt_sd, int band, FILE *fh) {
	int n_seg, file_len, pos_to, pos_from, seg_len_n, i, j, k, l_reg, ret_sz;
	int min_5dB_index, min_rtdB_index, nan_count = 0;
	int n_chan = 4 * samp, l = 0, discard = 0;
	double **chan, *chan_opt, *chan_opt_log, sd_tmp;
	double *sd_rev_time, max, min_5dB, min_rtdB, **mat_a, mat_a_mul[2];

	/* SVDLIBC stuff */
	SVDRec svd_mat = svdNewSVDRec();
	DMat dmat_a_mat;
	SVDVerbosity = 0;	/* 1: normal verbosity; 3: very verbose */

	/* Initialise */
	*mean_rt = 0.0;
	*rt_sd = 0.0;

	printf("Computing the Reverberation Time...\n");

	chan = calloc(array_size, sizeof(double *));
	chan_opt = calloc(n_chan, sizeof(double));
	chan_opt_log = calloc(n_chan, sizeof(double));

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
	 *	example, if seg_len_val = 180 secs and file_len = 6200
	 *	samples, then n_seg = 3 chunks (2 chunks of 180 secs and
	 *	1 shorter chunk). */
	seg_len_n = seg_len_val * samp;
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
					//printf("DR: %lf\n", dr[j]);
					optimum_model(a[j], b[j], alpha[j], dr[j],
							samp, n_chan, 
							&chan[k][0]);

					k++;
					l++;
				}
			}
		}

		printf("Segment: %d; Num Valid Decays: %d\n", i, k);
		fprintf(fh, "Segment: %d; Num Valid Decays: %d\n", i, k);

		if (k <= DISCARD_THRESH) {
			discard++;
			continue;
		}

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
		mat_a = calloc(l_reg, sizeof(double *));
		for (k = 0; k < l_reg; k++)
			mat_a[k] = calloc(2, sizeof(double));

		/* Get the RT */
		/* Note: due to (possible bug?) in SVDLIBC, I need to normalise
			the 2nd column of the matrix. I will multiply the RT
			at the end to compensate for the normalisation */
		for (j = 0; j < l_reg; j++) {
			mat_a[j][0] = 1.0;
			mat_a[j][1] = ((double) (min_5dB_index + j)) /
				min_rtdB_index;
		}

		/* Get the SVD of mat_a */
		dmat_a_mat = svdNewDMat(l_reg, 2);
		dmat_a_mat->value = mat_a;				

		svd_mat = svdLAS2A(svdConvertDtoS(dmat_a_mat), 0);

		/* Get reciprocal of S */
		for (k = 0; k < 2; k++)
			svd_mat->S[k] = 1.0 / svd_mat->S[k];

		/* Multiply modified S val by V */
		for (k = 0; k < 2; k++)
			for (j = 0; j < 2; j++)
				svd_mat->Vt->value[k][j] *= svd_mat->S[k];

		/* Matrix multiplcation of modified V with U */
		for (j = 0; j < l_reg; j++) {
			mat_a[j][0] = svd_mat->Ut->value[0][j] *
				svd_mat->Vt->value[0][0];

			mat_a[j][1] = svd_mat->Ut->value[1][j] *
				svd_mat->Vt->value[1][1];

			mat_a[j][0] += svd_mat->Ut->value[1][j] *
					svd_mat->Vt->value[1][0];

			mat_a[j][1] += svd_mat->Ut->value[0][j] *
					svd_mat->Vt->value[0][1];

			// DEBUG: MATLAB matches output when same values are
			//	passed to Matlab
        	}

		/* Multiply pseudoinverse by sig */
		mat_a_mul[0] = 0.0;
		mat_a_mul[1] = 0.0;
		for (j = 0; j < l_reg; j++) {
			mat_a_mul[0] += mat_a[j][0] *
				chan_opt_log[j + min_5dB_index];
			mat_a_mul[1] += mat_a[j][1] *
				chan_opt_log[j + min_5dB_index];
		}

		/* 
		 * If RT for segment != nan, add to rev_time
		 */
		if(!isnan(-60.0 / (samp * mat_a_mul[1]))) {
			sd_rev_time[i - nan_count - discard] =
				((double) (-60.0 * min_rtdB_index)) /
				(samp * mat_a_mul[1]);

			*mean_rt += sd_rev_time[i - nan_count - discard];
		} else
			nan_count++;

		/* Free array */
		for (k = 0; k < l_reg; k++)
			free(mat_a[k]);

		free(mat_a);
	}

	printf("n_seg: %d; nan_count: %d; discarded seg: %d\n",
		n_seg, nan_count, discard);
	printf("There are %d decays of at least -25 dB.\n", l);
	fprintf(fh, "n_seg: %d; nan_count: %d; discarded seg: %d\n",
		n_seg, nan_count, discard);
	fprintf(fh, "There are %d decays of at least -25 dB.\n", l);

	/* Get mean RT and standard dev of RT if num of valid RTs > 1 */
	if (n_seg - nan_count - discard > 1) {
		*mean_rt /= (n_seg - nan_count - discard);

		for (sd_tmp = 0.0, i = 0; i < n_seg - nan_count - discard; i++)
			sd_tmp += pow(fabs(sd_rev_time[i] - *mean_rt), 2);

		*rt_sd = sqrt(sd_tmp / (n_seg - nan_count - discard - 1));

		/* Get standard error from std dist (97.5 percentile) */
		*rt_sd = (*rt_sd / sqrt(n_seg - nan_count - discard)) * 1.96;
	}

	
	printf("Mean RT for %d band: %lf\n", band, *mean_rt);
	printf("RT Standard Dev for %d band: %lf\n", band, *rt_sd);
	fprintf(fh, "Mean RT for %d band: %lf\n", band, *mean_rt);
	fprintf(fh, "RT Standard Dev for %d band: %lf\n\n", band, *rt_sd);

	/* Free stuff */
	for (i = 0; i < array_size; i++)
		free(chan[i]);

	free(chan);
	free(chan_opt);
	free(chan_opt_log);
	free(sd_rev_time);
	svdFreeSVDRec(svd_mat);

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
	double n2, min, chan_sum[chan_sz], *last; 
	double win[nn * n_sect];

	last = calloc(n, sizeof(double));

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

	free(last);

	return (int) n2;
}

void hanning(int len, double *hann) {
	int n;

	for (n = 0; n < len; n++)
		hann[n] = 0.5 * (1.0 - cos(2.0 * M_PI * n / (len - 1.0)));

	return;
}

int process_wav_data(int band, float *wav_data, SF_INFO input_info, 
		SNDFILE *input, unsigned long long frames, double *p_alpha,
		double *p_a, double *p_b, double *p_dr,
		int *p_start, int *p_end, int *filt_frames) {
	int i, ret, s_e_size, resamp_frames = 0, **start_end;
	float *resampled_wav, *filtered_wav, *env;
	double *alpha, *a, *b, *dr;
	int *store_start, *store_end;

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
		return -1;
	}

	/* Resample the wav file */
	resampled_wav = resamp_wav_data(resamp, input_info.samplerate, 
		(uint64_t) frames, samp_freq_per_band[band],
		(const uint8_t **) &wav_data, &resamp_frames);

	swr_free(&resamp);
	if (resamp_frames <= 0 || resampled_wav == NULL) {
		printf("Resample error\nExiting...\n");
		return -1;
	}

	*filt_frames = resamp_frames;

	/* Octave band filtering */
	filtered_wav = oct_filt_data(resampled_wav, 
		band, (float) samp_freq_per_band[band],
		resamp_frames);

	/* Free memory no longer needed */
	av_free((uint8_t **) &resampled_wav[0]);
	resampled_wav = NULL;

	/* Obtain signal envelope */
	env = get_envelope((float) samp_freq_per_band[band], filtered_wav, 
		resamp_frames);

	/* Polyfit algorithm & get decay segments */
	start_end = apply_polyfit((float) samp_freq_per_band[band], 
		resamp_frames, env, &s_e_size);

	/* Perform Maximum Likelihood Estimation */
	ret = perform_ml(start_end, env, s_e_size, filtered_wav, 
		resamp_frames, samp_freq_per_band[band], &alpha, &a, &b, &dr,
		&store_start, &store_end);

	if (ret >= 0) {
		printf("Maximum Likelihood calculations have completed!\n%d "
			"potential segments analysed, %d may be valid\n",
			s_e_size, ret);

		/* Return pointers to caller */
		for (i = 0; i < ret; i++) {
			*(p_alpha + i) = alpha[i];
			*(p_a + i) = a[i];
			*(p_b + i) = b[i];
			*(p_dr + i) = dr[i];
			*(p_start + i) = store_start[i];
			*(p_end + i) = store_end[i];
		}
	}

	/* Free memory no longer needed */
	free(filtered_wav);
	filtered_wav = NULL;

	free(env);
	env = NULL;

	for (i = 0; i < 2; i++)
		free(start_end[i]);
	free(start_end);
	start_end = NULL;

	free(alpha);
	alpha = NULL;
	free(a);
	a = NULL;
	free(b);
	b = NULL;
	free(dr);
	dr = NULL;
	free(store_start);
	store_start = NULL;
	free(store_end);
	store_end = NULL;

	/* Return */
	if (ret < 0)
		return -1;
	else
		return ret;
}

float * resamp_wav_data(SwrContext *resamp, int in_rate, uint64_t num_frames, int samp_freq, const uint8_t **wav_data, int *resamp_frms) {
	uint8_t *output_data;
	int ret, resamp_frames;

	printf("Resampling wav...\n");
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
		av_free(&output_data[0]);
		output_data = NULL;
		return NULL;
	}

	/* Store # frames in pointer */
	*resamp_frms = resamp_frames;

	/* Return pointer to beginning of data */
	return (float *) output_data;
}

float * oct_filt_data(float *resamp_data, int band_num, float samp_freq, int resamp_frames) {
	int i, flag = 0;
	float *filtered = calloc(resamp_frames, sizeof(float)); 
	float band = octave_bands[band_num];

	memset(butt_b, 0, sizeof(butt_b));
	memset(butt_a, 0, sizeof(butt_a));

	printf("Applying octave band filter...\n");

#ifdef USE_HARD_BUTTER
	/* Only accurate when numbers are not tiny */
	if (fabs(butter_b[band_num][0]) > 1e-5)
		flag = 1;
#endif

	if (flag == 1) {
		iirfilt_rrrf f_obj_oct1 = iirfilt_rrrf_create(
			butter_b[band_num], 2 * OCT_FILT_ORDER + 1,
			butter_a[band_num], 2 * OCT_FILT_ORDER + 1);

		/* Filter data */
		for (i = 0; i < resamp_frames; i++)
			iirfilt_rrrf_execute(f_obj_oct1, resamp_data[i],
				&filtered[i]);

		/* Free filter and reset coefficient arrays */
		iirfilt_rrrf_destroy(f_obj_oct1);
	} else {
		/* Create bandpass filter - 7th and 8th parameter are ignored */
		liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_BANDPASS,
			LIQUID_IIRDES_SOS, OCT_FILT_ORDER,
			band / (sqrt(2.0) * samp_freq),
			band / samp_freq, 1.0, 1.0, butt_b, butt_a);
		//printf("res: %d\n\n", iirdes_isstable(butt_b, butt_a,
		//	3 * (L_OCT + R_OCT)));
		iirfilt_rrrf f_obj_oct2 = iirfilt_rrrf_create_sos(
			butt_b, butt_a, L_OCT + R_OCT);

		/* Filter data */
		for (i = 0; i < resamp_frames; i++)
			iirfilt_rrrf_execute(f_obj_oct2, resamp_data[i],
				&filtered[i]);

		/* Free filter and reset coefficient arrays */
		iirfilt_rrrf_destroy(f_obj_oct2);
	}

	/* filtered data: return to previous function */
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
	printf("Obtaining signal envelope...\n");
	liquid_iirdes(LIQUID_IIRDES_BUTTER, LIQUID_IIRDES_LOWPASS, 
		LIQUID_IIRDES_TF, ENV_FILT_ORDER,
		LOW_PASS_CUTOFF / samp_freq, 0.1, 
		1.0, 1.0, butt_b, butt_a);

	iirfilt_crcf f_obj_env = iirfilt_crcf_create(butt_b,
		ENV_FILT_ORDER+1, butt_a, ENV_FILT_ORDER+1);

	fftplan fft_pln = fft_create_plan(resamp_frames, (float complex *) 
		filt_complex, (float complex *) tmp, LIQUID_FFT_FORWARD, 0);
	fftplan ifft_pln = fft_create_plan(resamp_frames, (float complex *) tmp,
 		(float complex *) hilb, LIQUID_FFT_BACKWARD, 0);

	/* Get the Hilbert transform - really it's the analytic signal as the
		Hilbert transform is stored in the imaginary part of the soln */
	for (i = 0; i < resamp_frames; i++)
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
	}

	free(tmp2);
	free(tmp);
	free(hilb);
	free(filt_complex);
	iirfilt_crcf_destroy(f_obj_env);
	fft_destroy_plan(ifft_pln);
	fft_destroy_plan(fft_pln);

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

	printf("Applying polyfit algorithm...\n");

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

	/* Free arrays */
	for (i = 0; i < num_segs; i++) 
		free(poly_res[i]);
	free(poly_res);

	for (i = 0; i < num_segs - dump; i++)
		free(seg_index[i]);
	free(seg_index);
	seg_index = NULL;

	free(poly_seg);
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
		double **p_a_par, double **p_b_par, double **p_dr,
		int **p_start, int **p_end) {
	int sz, i, all_ret = 0;
	float *abs_filtered = calloc(resamp_frames, sizeof(float));
	int *store_start = calloc(s_e_size, sizeof(int)); 
	int *store_end = calloc(s_e_size, sizeof(int));
	int *len_store = calloc(s_e_size, sizeof(int));
	double *a_par = calloc(s_e_size, sizeof(double));
	double *b_par = calloc(s_e_size, sizeof(double));
	double *alpha = calloc(s_e_size, sizeof(double));
	double *dr = calloc(s_e_size, sizeof(double));

	printf("Performing Maximum Likelihood calculations\n"
		"This will take a while...\n");

	/* Get absolute vaule of filtered wav */
	for (i = 0; i < resamp_frames; i++)
		abs_filtered[i] = fabs(filtered_wav[i]);	

	/* Largest possible size of arrays - first pair of elements */
	sz = start_end[1][0] - start_end[0][0];

	#pragma omp parallel for schedule(static)
	for (i = 0; i < s_e_size; i++) {
		float *segment = calloc(sz, sizeof(float));
		float *env_seg = calloc(sz, sizeof(float));
		float *abs_filt_seg = calloc(sz, sizeof(float));
		float *seg_seg2 = calloc(sz + sz - roundf((DAMP_SEG_SIZE / 2) *
			samp_freq), sizeof(float));
		float max_seg2; 
		int k, j, len = abs(start_end[1][i] - start_end[0][i]);	
		int min_abs_filt_seg, max_abs_filt_seg, ret = 0;
	
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

		/* ML fitting of decay model to the data */
		ret = ml_fit(seg_seg2, len_store[i], (float) samp_freq, 
			&a_par[i], &b_par[i], &alpha[i]);

		/* Reset to zero b/c arrays are longer than they need to be */
		free(segment);
		free(seg_seg2);
		free(env_seg);
		free(abs_filt_seg);

		if (ret == -1) {
			/* DR being discarded */
			dr[i] = 0;
			continue;
		} else if (ret < -1) {
			/* Error */
			all_ret = -1;
		}

		/* If code gets to here, DR was possibly valid */
		/* Compute decay curves using ML-calculated params */
		dr[i] = get_decay_region(6 * samp_freq, a_par[i], b_par[i],
			alpha[i], len_store[i]);
	}


	/* Pointers used by caller cannot be freed. They are alpha[], a[],
		b[], store_start[], store_end[] and dr[] */
	*p_alpha = alpha;
	*p_a_par = a_par;
	*p_b_par = b_par;
	*p_dr = dr;
	*p_start = store_start;
        *p_end = store_end;

	/* Free pointers that are no longer needed */
	free(abs_filtered);
	free(len_store);

	if (all_ret != 0) {
		printf("Error during ML calculations!\n");
		return -1;
	} else
		return s_e_size;
}

double get_decay_region(int c_len, double a, double b, double alpha,
		int len_store) {
	double *chan = calloc(c_len, sizeof(double));
	double *y = calloc(c_len, sizeof(double));
	double dr;
	int j, max_y;

	/* Compute decay curves using ML-calculated params */
	for (j = 0; j < c_len; j++)
		chan[j] = alpha * pow(a, j) + (1.0 - alpha) * pow(b, j);

	/* Get the reverse cumsum of chan[] and assign to y in reverse */
	y[c_len - 1] = pow(chan[c_len - 1], 2);
	for (j = c_len - 2; j >= 0; j--)
		y[j] = pow(chan[j], 2) + y[j + 1];

	/* Get max of fabs(y) */
	max_y = fabs(y[0]);
	for (j = 1; j < c_len; j++) {
		if (fabs(y[j]) > max_y)
			max_y = fabs(y[j]);
	}

	/* Normalise y first value will be max; stores whole cumsum */
	for (j = 0; j < c_len; j++) 
		y[j] = 10.0 * log10(y[j] / max_y);

	dr = y[len_store];

	free(chan);
	free(y);

	return dr;
}

int ml_fit(float *data_seg, int len, float samp_freq, double *a_par, 
		double *b_par, double *alph_par) {
	double min, max, interval, j, gmax_val, x_fine[3], fine_val, dr;
	double *coarse_grid, lb[3], ub[3];
	double like[SQP_STEP][SQP_STEP] = {0}, alpha[SQP_STEP][SQP_STEP] = {0};
	int gmax_pos[2], i, k, el_cg, ret = 0;
	struct nl_extra_data nld;
	nld.data_seg = data_seg;
	nld.len = len;

	lb[0] = 0.0;
	ub[0] = 1.0;

	nlopt_opt nl_obj1 = nlopt_create(NLOPT_LN_COBYLA, 1);
	nlopt_set_lower_bounds1(nl_obj1, lb[0]);
	nlopt_set_upper_bounds1(nl_obj1, ub[0]);
	nlopt_set_maxeval(nl_obj1, 100);
	nlopt_set_min_objective(nl_obj1,
		 (nlopt_func) alpha_opt, (void *) &nld);

	min = (-6.91 / log(PAR_LOW_BOUND)) / 3000.0;
	min = exp(-6.91 / (samp_freq * min));

	max = (-6.91 / log(PAR_UP_BOUND)) / 3000.0;	
	max = exp(-6.91 / (samp_freq * max));
	
	interval = (max - min) / ((double) SQP_STEP - 1.0);
	el_cg = round((max - min) / interval) + 1;
	coarse_grid = calloc(el_cg, sizeof(double));

	/* Fill in coarse grid */
	for (i = 0, j = min; i < el_cg; i++, j += interval)
		coarse_grid[i] = j;

	// There is another max(abs(x)) division here in MATLAB
	// 	Makes no diff - dividing by 1.0 - leaving out
	
	/* Coarse grid minimisation search */
	for (i = 0; i < SQP_STEP; i++) {
		nld.b_val = coarse_grid[i];
		for (k = 0; k < SQP_STEP; k++) {
			nld.a_val = coarse_grid[k];
			alpha[k][i] = 0.5;

			ret = nlopt_optimize(nl_obj1, &alpha[k][i], 
				&like[k][i]);

			/* -4 is an acceptable nlopt return code */
			if (ret < 0 && ret != -4) {
				nlopt_destroy(nl_obj1);
				free(coarse_grid);
				return -2;
			}
			like[k][i] *= -1;
		}
	}

	for (i = 0; i < SQP_STEP; i++) {
		for (k = 0; k < SQP_STEP; k++) {	
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

			//printf("k: %d, i: %d, ret: %d, a: %lf like: %le coarse: %lf\n", k+1, i+1, ret, alpha[k][i], like[k][i], coarse_grid[k]);
		}
	}

	/* Not needed anymore */
	nlopt_destroy(nl_obj1);

#ifdef COARSE_DR_EST
	dr = get_decay_region(6 * samp_freq, coarse_grid[gmax_pos[0]],
		coarse_grid[gmax_pos[1]], alpha[gmax_pos[0]][gmax_pos[1]], len);

	if (dr > -1.0) {
		free(coarse_grid);
		return -1;
	}
#endif

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

	/* Initial values for fine search */
	x_fine[0] = coarse_grid[gmax_pos[0]];
	x_fine[1] = coarse_grid[gmax_pos[1]];
	x_fine[2] = alpha[gmax_pos[0]][gmax_pos[1]];	
	//printf("B4 XVAL: %lf %lf %lf %le\n", x_fine[0], x_fine[1], x_fine[2], gmax_val); 

	ret = nlopt_optimize(nl_obj3, x_fine, &fine_val);
	//printf("DEBUG FINE RESULTS - ret: %d XVAL: %lf %lf %lf FINE_VAL: %le\n", 
		//ret, x_fine[0], x_fine[1], x_fine[2], fine_val);
	
	/* Return a, b, alpha in pointers */
	*a_par = x_fine[0];
	*b_par = x_fine[1];
	*alph_par = x_fine[2];
	
	free(coarse_grid);
	nlopt_destroy(nl_obj3);

	if (ret < 0 && ret != -4)
		return -2;
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
