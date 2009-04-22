/*
 * Copyright (c) 2009 Matt Jibson <matt.jibson@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sndfile.h>

extern char *__progname;

/*
 * Returns the value of the zeroth order modified Bessel function of the first
 * kind at x.
 */
double
bessel_i0(double x)
{
	double i0 = 1, diff, ifac = 1;
	int i = 1;
	
	for (;;)
	{
		diff = pow(x / 2.0, i * 2) / pow(ifac, 2);

		if (diff < 1.0e-21)
			break;

		i0 += diff;
		i++;
		ifac *= i;
	}

	return i0;
}

/*
 * Returns the value of the Kaiser window at point n with shape alpha and
 * length M+1.
 */
double
kaiser(double alpha, double M, double n)
{
	if (n < 0 || n > M)
		return 0;

	return bessel_i0(alpha * sqrt(1.0 - pow(2.0 * n / M - 1.0, 2))) /
		bessel_i0(alpha);
}

/*
 * Returns a pointer to an array of length len representing one half of a
 * Kaiser sequence of size (len - 1) * 2 with shape alpha.
 */
double *
kaiser_table(double alpha, int len)
{
	double *ret;
	int i;
	int M = (len - 1) * 2;

	if ((ret = (double *)calloc(len, sizeof(double))) == NULL)
		return NULL;

	for (i = 0; i < len; i++)
		ret[i] = kaiser(alpha, M, len - i - 1);

	return ret;
}

/*
 * Returns the alpha (aka beta) to use as the Kaiser shape parameter when
 * given a sidelobe attenuation of dB, as specified by the MATLAB
 * documentation on its kaiser() function.
 */
double
kaiser_dB(double dB)
{
	if(dB > 50.0)
		return 0.1102 * (dB - 8.7);
	else if(dB >= 21.0)
		return 0.5842 * pow(dB - 21.0, 0.4) + 0.07886 * (dB - 21.0);
	else
		return 0;
}

/*
 * Bandlimited interpolation resampling algorithm based on work by
 * Julius O. Smith III: http://www-ccrma.stanford.edu/~jos/resample/.
*/
void
bi_resamp(char * in, char * out, int rate)
{
	int i, j;
	FILE * fd;
	SNDFILE * sf;
	SF_INFO sfinfo;

	const int n_c = sizeof(short) * 8; // word-length of stored impulse-response samples
	int n = 0; // samples into signal buffer (index-like)
	int l = 0; // index into filter coefficient table
	double eta = 0; // between 0 and 1 for doing linear interpolation between samples l and l+1 (initially) of filter table
	double t = 0; // desired time output
	const int n_n = sizeof(short) * 8; // word-length of n
	const int n_l = 1 + n_c / 2;
	const int n_eta = n_c / 2;
	const int N = 1 << n_n; // input signal buffer contains this many samples
	const int L = 1 << n_l; // filter table contains this many samples per zero-crossing
	const int N_z = 5; // number of zero-crossings
	const int N_h = L * N_z + 1; // filter table length

	if ((sf = sf_open(in, SFM_READ, &sfinfo)) == NULL)
		errx(1, "could not open input file: %s", in);

	const int Fs  = sfinfo.samplerate; // input rate
	const int Fsp = rate; // output rate
	double rho = (double)Fsp / Fs; // sampling-rate conversion factor
	short * x;
	short * dat;
	const int extra = rho < 1 ? ceil((double)(N_z * Fs) / Fsp) : N_z;
	const int dat_len = sfinfo.frames * sfinfo.channels;
	const int x_len = sfinfo.frames + 2 * extra;

	if ((x = (short *)calloc(x_len, sizeof(short))) == NULL)
		exit(1);

	if ((dat = (short *)calloc(dat_len, sizeof(short))) == NULL)
		exit(1);

	if(sf_readf_short(sf, dat, sfinfo.frames) != sfinfo.frames)
		exit(1);
	sf_close(sf);

	for(i = 0; i < x_len; i++)
		if(i < extra || i >= x_len - extra)
			x[i] = 0;
		else
			x[i] = dat[(i - extra) * sfinfo.channels];

	const double cutoff = kaiser_dB(80);
	double * h;
	double * hb;
	double sinc;

	if((h = kaiser_table(cutoff, N_h)) == NULL)
		exit(1);

	for(i = 1; i < N_h; i++)
	{
		sinc = (double)i / L * M_PI;
		h[i] *= sin(sinc) / sinc;
	}

	h[0] = 1.0;

	if ((hb = (double *)calloc(N_h, sizeof(double))) == NULL)
		exit(1);

	for(i = 1; i < N_h; i++)
		hb[i - 1] = h[i] - h[i - 1];
	hb[N_h - 1] = 0;

	short v;
	short * y;
	const int y_len = (double)(sfinfo.frames * Fsp / sfinfo.samplerate);

	if ((y = (short *)calloc(y_len, sizeof(short))) == NULL)
		exit(1);

	double xt, xtn;
	int step;

	for(j = 0; j < y_len; j++)
	{
		v = 0;

		n = t * Fs;
		xt = (double)n / Fs;
		xtn = (double)(n + 1) / Fs;
		eta = 1 - (xtn - t) / (xtn - xt);
		l = eta * L;
		step = L;

		n += extra;

		for(i = 0; l + i * step < N_h; i++)
			v += x[n - i] * (h[l + i * step] + eta * hb[l + i * step]);

		eta = 1 - eta;
		l = eta * L;
		step = L;

		for(i = 0; l + i * step < N_h; i++)
			v += x[n + 1 + i] * (h[l + i * step] + eta * hb[l + i * step]);

		y[j] = v;

		t += 1.0 / Fsp;
	}

	sfinfo.channels = 1;
	sfinfo.samplerate = Fsp;
	sfinfo.frames = y_len;
	if((sf = sf_open(out, SFM_WRITE, &sfinfo)) == NULL)
		errx(1, "could not open output file: %s", out);
	sf_close(sf);
}

int main(int argc, char *argv[])
{
	int rate;
	const char *errstr;

	if (argc != 4)
	{
	  fprintf(stderr,
	  	"usage: %s input rate output\n", __progname);
	  exit(1);
	}

	rate = strtonum(argv[2], 1, INT_MAX, &errstr);
	if (errstr)
		errx(1, "rate is %s: %s", errstr, argv[2]);

	bi_resamp(argv[1], argv[3], rate);

	return 0;
}
